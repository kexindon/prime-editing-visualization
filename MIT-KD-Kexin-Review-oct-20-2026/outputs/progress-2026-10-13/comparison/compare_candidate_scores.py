#!/usr/bin/env python3
"""Compare observed score rankings after a strict candidate identity join.

No experimental efficiencies are present. Agreement here does not measure
accuracy. The join rejects ambiguous/duplicate identities and unexpected
sequence geometry instead of using a nearest or approximate match.
"""
import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from candidate_identity import identity_hash, pegg_key, pridict_key
from optiprime_scaling import digest, save_json


def finite_score(value):
    result = float(value)
    if not np.isfinite(result):
        raise ValueError('Non-finite score in comparison input')
    return result


def rank_comparisons(df, columns, cohort):
    records = []
    for mutation_idx, part in df.groupby('mutation_idx', sort=True):
        for a, b in itertools.combinations(columns, 2):
            x, y = part[a].astype(float), part[b].astype(float)
            ra, rb = x.rank(method='average'), y.rank(method='average')
            defined = len(part) >= 2 and ra.nunique() > 1 and rb.nunique() > 1
            rho = float(ra.corr(rb)) if defined else None
            ids = part.candidate_identity_sha256
            top_a, top_b = set(ids[x.eq(x.max())]), set(ids[y.eq(y.max())])
            records.append({
                'cohort': cohort, 'mutation_idx': int(mutation_idx),
                'sequence_name': part.sequence_name.iloc[0], 'candidate_count': len(part),
                'score_a': a, 'score_b': b, 'spearman': rho,
                'undefined_reason': '' if defined else 'fewer than two candidates or constant scores',
                'a_tied_rows': int(x.duplicated(keep=False).sum()),
                'b_tied_rows': int(y.duplicated(keep=False).sum()),
                'a_top_ties': len(top_a), 'b_top_ties': len(top_b),
                'top_sets_intersect': bool(top_a & top_b),
                'top_sets_equal': top_a == top_b,
            })
    per_target = pd.DataFrame(records)
    summary = []
    for (a, b), part in per_target.groupby(['score_a', 'score_b'], sort=False):
        rhos = part.spearman.dropna()
        summary.append({
            'cohort': cohort, 'score_a': a, 'score_b': b,
            'target_count': len(part), 'defined_spearman_targets': len(rhos),
            'undefined_spearman_targets': len(part) - len(rhos),
            'median_within_target_spearman': float(rhos.median()) if len(rhos) else None,
            'mean_within_target_spearman': float(rhos.mean()) if len(rhos) else None,
            'min_within_target_spearman': float(rhos.min()) if len(rhos) else None,
            'max_within_target_spearman': float(rhos.max()) if len(rhos) else None,
            'targets_with_top_set_intersection': int(part.top_sets_intersect.sum()),
            'targets_with_equal_top_sets': int(part.top_sets_equal.sum()),
            'targets_with_tied_top_a': int(part.a_top_ties.gt(1).sum()),
            'targets_with_tied_top_b': int(part.b_top_ties.gt(1).sum()),
        })
    return per_target, summary


def runtime_agreement(runtime_dir):
    reference = pd.read_csv(runtime_dir / 'n1000_w1/scored_candidates.csv', keep_default_na=False)
    reference = reference.set_index('candidate_identity_sha256')
    records = []
    differences = []
    for run in json.loads((runtime_dir / 'runs.json').read_text()):
        if run['exit_code'] != 0:
            raise ValueError('Failed runtime experiment: ' + run['run'])
        if run['run'] == 'n1000_w1':
            continue
        other = pd.read_csv(runtime_dir / run['run'] / 'scored_candidates.csv', keep_default_na=False)
        other = other.set_index('candidate_identity_sha256')
        if other.index.has_duplicates or not other.index.isin(reference.index).all():
            raise ValueError('Invalid run candidate identity set')
        delta = other.OptiPrime_Score - reference.loc[other.index, 'OptiPrime_Score']
        if not np.isfinite(delta).all():
            raise ValueError('Nonfinite score differences')
        records.append({'run': run['run'], 'reference': 'n1000_w1', 'compared_rows': len(delta),
                        'max_absolute_difference': float(delta.abs().max()),
                        'median_absolute_difference': float(delta.abs().median()),
                        'rows_above_1e_6_absolute_difference': int(delta.abs().gt(1e-6).sum())})
        for cid, difference in delta.items():
            differences.append({'run': run['run'], 'candidate_identity_sha256': cid,
                                'difference_from_n1000_w1': difference})
    return records, pd.DataFrame(differences)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--runtime-dir', type=Path, required=True)
    ap.add_argument('--pridict-dir', type=Path, required=True)
    ap.add_argument('--mutation-lookup', type=Path, required=True)
    ap.add_argument('--historical-scores', type=Path)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    source = args.runtime_dir / 'n1000_w1/scored_candidates.csv'
    scored = pd.read_csv(source, keep_default_na=False)
    lookup = pd.read_csv(args.mutation_lookup, keep_default_na=False)
    if lookup.sequence_name.duplicated().any() or lookup.mutation_idx.duplicated().any():
        raise ValueError('Duplicate mutation lookup IDs')
    identities = dict(zip(lookup.sequence_name, zip(lookup.WT, lookup.ALT)))
    if len(set(identities.values())) != len(identities):
        raise ValueError('Duplicate target sequence identities')
    scored = scored.merge(lookup[['mutation_idx', 'sequence_name']], on='mutation_idx', validate='many_to_one')
    peg_keys = [pegg_key(row) for _, row in scored.iterrows()]
    if len(set(peg_keys)) != len(peg_keys):
        raise ValueError('Duplicate PEGG candidate keys')
    if [identity_hash(key) for key in peg_keys] != scored.candidate_identity_sha256.tolist():
        raise ValueError('Stored candidate hashes disagree with sequences')
    for name, key in zip(scored.sequence_name, peg_keys):
        if key[:2] != identities[name]:
            raise ValueError('PEGG target does not match lookup')

    manifest = args.pridict_dir / 'per_target_status.csv'
    status = pd.read_csv(manifest, keep_default_na=False)
    if set(status.sequence_name) != set(identities) or status.sequence_name.duplicated().any():
        raise ValueError('PRIDICT2 target manifest is not one-to-one with lookup')
    predictions, file_hashes = {}, {}
    for _, entry in status.iterrows():
        if entry.outcome != 'success':
            raise ValueError('Unsuccessful PRIDICT2 target')
        path = (args.pridict_dir / entry.candidate_output).resolve()
        if not path.is_relative_to(args.pridict_dir.resolve()):
            raise ValueError('Candidate output path escapes PRIDICT2 output directory')
        file_hashes[entry.candidate_output] = digest(path)
        count = 0
        with path.open(newline='') as handle:
            for line, row in enumerate(csv.DictReader(handle), start=2):
                count += 1
                key = pridict_key(row)
                if row['sequence_name'] != entry.sequence_name or key[:2] != identities[entry.sequence_name]:
                    raise ValueError('PRIDICT2 named target does not match fixture')
                if key in predictions:
                    raise ValueError('Ambiguous duplicate PRIDICT2 candidate key')
                predictions[key] = {
                    'PRIDICT2_HEK': finite_score(row['PRIDICT2_0_editing_Score_deep_HEK']),
                    'PRIDICT2_K562': finite_score(row['PRIDICT2_0_editing_Score_deep_K562']),
                    'pridict_spacer': row['Spacer-Sequence'],
                    'pridict_distance_to_nick': int(row['Editing_Position']),
                    'pridict_rha': int(row['RToverhanglength']),
                    'pridict_source_file': entry.candidate_output, 'pridict_source_csv_line': line,
                }
        if count != int(entry.candidate_rows):
            raise ValueError('PRIDICT2 count differs from manifest')

    joined = []
    for (_, row), key in zip(scored.iterrows(), peg_keys):
        pr = predictions.get(key)
        record = row.to_dict()
        if pr:
            if pr['pridict_spacer'] != row.Protospacer:
                raise ValueError('Matched native geometry has different synthesized spacer')
            if pr['pridict_distance_to_nick'] != int(row.Distance_to_nick) or pr['pridict_rha'] != int(row.RHA_size):
                raise ValueError('Matched sequence has inconsistent editing geometry')
            record.update(pr)
            record['match_status'] = 'exact'
        elif int(row.RHA_size) > 19:
            record['match_status'] = 'outside_pridict_rha_max19'
        elif int(row.RHA_size) < 3:
            record['match_status'] = 'outside_pridict_rha_min3'
        elif int(row.Distance_to_nick) > 25:
            record['match_status'] = 'outside_pridict_nick_distance_max25'
        else:
            raise ValueError('Unexplained absent candidate: ' + row.pegRNA_id)
        joined.append(record)
    joined = pd.DataFrame(joined)
    matched = joined[joined.match_status.eq('exact')].copy()
    joined.to_csv(args.output_dir / 'candidate_coverage.csv', index=False)
    matched.to_csv(args.output_dir / 'exact_matches.csv', index=False)
    per_target = joined.groupby(['mutation_idx', 'sequence_name', 'match_status']).size().unstack(fill_value=0)
    per_target['pegg_candidates'] = per_target.sum(axis=1)
    per_target['match_fraction'] = per_target['exact'] / per_target.pegg_candidates
    per_target.to_csv(args.output_dir / 'target_coverage.csv')
    summary = {
        'pegg_candidates': len(scored), 'pridict_candidates': len(predictions),
        'matched_candidates': len(matched), 'matched_fraction_of_pegg': len(matched) / len(scored),
        'targets': int(scored.mutation_idx.nunique()), 'matched_targets': int(matched.mutation_idx.nunique()),
        'minimum_matches_per_target': int(per_target['exact'].min()),
        'maximum_matches_per_target': int(per_target['exact'].max()),
        'match_status_counts': joined.match_status.value_counts().to_dict(),
        'matched_strands': matched.PAM_strand.value_counts().to_dict(),
        'matched_edit_types': matched.Variant_Type.value_counts().to_dict(),
        'duplicate_or_ambiguous_identities': 0,
        'interpretation': 'Exact target/guide geometry only. Models use different experimental metadata. No measured efficiencies; no biological accuracy conclusion.',
    }
    save_json(args.output_dir / 'coverage_summary.json', summary)
    baseline_columns = ['PEGG2_Score', 'RF_Score', 'OptiPrime_Score']
    frames, summaries = [], []
    for df, cols, cohort in ((joined, baseline_columns, 'all_1000_pegg_candidates'),
                             (matched, baseline_columns + ['PRIDICT2_HEK', 'PRIDICT2_K562'], 'exact_matched_subset')):
        if not np.isfinite(df[cols].astype(float).to_numpy()).all():
            raise ValueError('Missing/nonfinite score within declared comparison cohort')
        frame, stats = rank_comparisons(df, cols, cohort)
        frames.append(frame)
        summaries.extend(stats)
    pd.concat(frames, ignore_index=True).to_csv(args.output_dir / 'within_target_ranks.csv', index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / 'rank_summary.csv', index=False)
    save_json(args.output_dir / 'rank_summary.json', summaries)
    agreement, deltas = runtime_agreement(args.runtime_dir)
    save_json(args.output_dir / 'runtime_score_agreement.json', agreement)
    deltas.to_csv(args.output_dir / 'runtime_score_differences.csv', index=False)

    if args.historical_scores:
        old = pd.read_csv(args.historical_scores, keep_default_na=False)
        old['candidate_identity_sha256'] = [identity_hash(pegg_key(row)) for _, row in old.iterrows()]
        now = pd.read_csv(args.runtime_dir / 'n100_w1/scored_candidates.csv', keep_default_na=False)
        pair = old[['candidate_identity_sha256', 'OptiPrime_Score']].merge(
            now, on='candidate_identity_sha256', suffixes=('_historical', '_corrected'), validate='one_to_one')
        if len(pair) != len(old) or len(pair) != len(now):
            raise ValueError('Historical comparison candidate sets differ')
        pair['score_difference'] = pair.OptiPrime_Score_corrected - pair.OptiPrime_Score_historical
        pair['historical_blank_allele_affected'] = pair.REF.eq('') | pair.ALT.eq('')
        pair.to_csv(args.output_dir / 'historical_100_score_changes.csv', index=False)

    save_json(args.output_dir / 'provenance.json', {
        'input_scored_sha256': digest(source), 'mutation_lookup_sha256': digest(args.mutation_lookup),
        'pridict_manifest_sha256': digest(manifest), 'pridict_candidate_file_hashes': file_hashes,
        'script_sha256': digest(__file__),
        'identity_script_sha256': digest(Path(__file__).with_name('candidate_identity.py')),
        'historical_scores_sha256': digest(args.historical_scores) if args.historical_scores else None,
        'strict_key': ['full_forward_WT', 'full_forward_ALT', 'PAM_strand', 'PAM_oriented_start', 'native30', 'RTT', 'PBS'],
        'normalization': 'PRIDICT2 Edited_Sequence and RTrevcomp uppercased to remove introduced-allele case annotation; Rv full sequences reverse-complemented. No OptiPrime spacer normalization.',
        'within_target_spearman': 'Pearson correlation of average ranks, computed separately per mutation on the stated subset; constants undefined; summaries weight targets equally.',
        'top_agreement': 'Exact maximum score sets; report both nonempty intersection and equal sets, and retain tie counts.',
    })
    for path in (Path(__file__), Path(__file__).with_name('candidate_identity.py')):
        (args.output_dir / path.name).write_bytes(path.read_bytes())
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
