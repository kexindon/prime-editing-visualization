#!/usr/bin/env python3
"""Render the dated report from saved inference and comparison artifacts."""
import argparse
import json
from pathlib import Path

import pandas as pd


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(str(x) for x in row) + ' |' for row in rows])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workspace', type=Path, required=True)
    args = ap.parse_args()
    root = args.workspace.resolve()
    out = root / 'outputs/progress-2026-10-13'
    runs = json.loads((out / 'runtime/runs.json').read_text())
    coverage = json.loads((out / 'comparison/coverage_summary.json').read_text())
    agreement = json.loads((out / 'comparison/runtime_score_agreement.json').read_text())
    ranks = pd.read_csv(out / 'comparison/rank_summary.csv')
    historical = pd.read_csv(out / 'comparison/historical_100_score_changes.csv', keep_default_na=False)
    corrected = historical[historical.historical_blank_allele_affected]
    source = pd.read_csv(root / 'outputs/meeting-prep/pegg_100/pegg_candidates.csv', keep_default_na=False)
    expected_poly = source.RTT_PBS.str.contains('TTTT', regex=False) | source.Protospacer.str.contains('TTTT', regex=False)
    poly = source[['pegRNA_id', 'mutation_idx', 'Protospacer', 'RTT_PBS', 'contains_polyT_terminator']].copy()
    poly['expected_contains_polyT_terminator'] = expected_poly
    poly['annotation_changed'] = expected_poly.ne(source.contains_polyT_terminator)
    poly.to_csv(out / 'polyT_annotation_impact.csv', index=False)
    by_run = {r['run']: r for r in runs}
    single, parallel = by_run['n1000_w1'], by_run['n1000_w2']
    parallel_delta = next(r['max_absolute_difference'] for r in agreement if r['run'] == 'n1000_w2')
    runtime_table = table(
        ['Guides', 'Workers', 'Scorer call (s)', 'Whole process (s)', 'Finite scores', 'Largest child RSS (GB)'],
        [[r['candidate_count'], r['actual_workers'], f"{r['scorer_wall_seconds']:.2f}",
          f"{r['process_wall_seconds']:.2f}", f"{r['finite_scores']}/{r['candidate_count']}",
          f"{r['largest_reaped_child_peak_rss_bytes']/1e9:.3f}"] for r in runs])
    edit_table = table(['Guide', 'Edit', 'Historical score', 'Corrected score'],
                      [[r.pegRNA_id, r.Variant_Type, f'{r.OptiPrime_Score_historical:.9f}',
                        f'{r.OptiPrime_Score_corrected:.9f}'] for r in corrected.itertuples()])
    wanted = {('PEGG2_Score', 'OptiPrime_Score'), ('PEGG2_Score', 'PRIDICT2_HEK'),
              ('PEGG2_Score', 'PRIDICT2_K562'), ('OptiPrime_Score', 'PRIDICT2_HEK'),
              ('OptiPrime_Score', 'PRIDICT2_K562')}
    shown = ranks[ranks.cohort.eq('exact_matched_subset')]
    rank_table = table(['Score pair', 'Median within-target Spearman', 'Same top guide /100 targets'],
                       [[f'{r.score_a} / {r.score_b}', f'{r.median_within_target_spearman:.3f}',
                         r.targets_with_equal_top_sets] for r in shown.itertuples()
                        if (r.score_a, r.score_b) in wanted])
    all_pair = ranks[(ranks.cohort == 'all_1000_pegg_candidates') &
                     (ranks.score_a == 'PEGG2_Score') & (ranks.score_b == 'OptiPrime_Score')].iloc[0]
    report = f'''# MIT KD — PEGG progress, October 13, 2026

## Completed this round

- Added reproducible runtime and exact-candidate comparison scripts on local
  branch `codex/pegg-model-prep` (base PEGG commit `07eda419`).
- Scored all **1,000 existing PEGG guides across 100 targets**, with all five
  OptiPrime models. All five experiment processes exited 0; all returned scores
  were finite and candidate order was preserved.
- Joined **{coverage['matched_candidates']}/1,000 guides ({coverage['matched_fraction_of_pegg']:.1%})**
  exactly to saved PRIDICT2 predictions, covering **all 100 targets** with 2–9 matches
  each and no duplicate or ambiguous identities.
- Corrected three engineering issues: blank-allele handling, score assignment by
  row position, and polyT annotation. Model definitions and weights are unchanged.

**Status:** execution, scaling and candidate-comparison evidence are ready for
review. Kexin's mapping/metadata decision and measured sensor data are still
needed before interpreting these predictions as better experimental designs.

## 1. Runtime and numerical agreement

{runtime_table}

These are single observations on this Mac (arm64, macOS 27, 14 logical CPUs).
Fixed settings: CPU, numerical thread environment caps of 1, batch size 128, `Liu_HeLa`,
time 4, all five models. Candidate prefixes are nested and all cover 100 targets.
They are candidate counts, not counts of additional mutations.

The scorer call includes model worker startup, feature generation, inference and
output transport. Whole-process time also includes parent imports and archiving
the inputs. Time outside the scorer call was larger in the first process, so
100 guides took longer than 250 in whole-process time. Cache state and one repetition per configuration
prevent a firm general speed claim. These times exclude PEGG/PRIDICT2 design.

For 1,000 guides, two workers gave a **{single['scorer_wall_seconds']/parallel['scorer_wall_seconds']:.2f}× observed scorer speedup**.
The largest absolute prediction difference from the single-worker result was
**{parallel_delta:.3g}**; no difference exceeded 1e-6. The 100/250/500-guide results were
identical to their corresponding rows in the 1,000-guide single-worker run.
This measures numerical consistency for this fixture, not biological validity.

The RSS column is the largest individual child high-water mark. It is **not
aggregate memory for two simultaneous workers**. Parent peaks and byte values
are in [runtime.csv](runtime/runtime.csv).

## 2. Correction to the earlier 100-guide smoke

The old loader used default `pd.read_csv`, converting empty `REF`/`ALT` indel
alleles to NaN. The adapter then counted the three characters in `str(NaN)` as
bases when calculating context length. Both affected guides returned finite
scores, so the old 100/100 success claim demonstrated execution only.

{edit_table}

For `pegRNA_0` (deletion1), the unedited crop was 62 nt instead of 65 nt. For
`pegRNA_20` (insertion1), it was 55 nt instead of 52 nt. Edited crops were 55 nt in both.
The corrected experiment preserves empty strings with `keep_default_na=False`.
The package now rejects non-string alleles with a per-row explanation. The
maximum historical score change among the other 98 guides was
{historical.loc[~historical.historical_blank_allele_affected, 'score_difference'].abs().max():.3g}.

**The two old indel scores are superseded.** Original raw outputs and the sent
October 11 ZIP are preserved as historical artifacts. Use the dated corrected
scores and [full old/new comparison](comparison/historical_100_score_changes.csv).

## 3. Candidate coverage

| Match status | PEGG guides |
| --- | ---: |
| Exact match |553|
| Homology arm above PRIDICT2 maximum 19 nt |332|
| Homology arm below PRIDICT2 minimum 3 nt |109|
| Nick-to-edit distance 26 exceeds PRIDICT2 default 25 nt |6|

All 1,000 guides are accounted for. The unmatched 447 reflect the saved native
PRIDICT2 design search, not evidence that its underlying model cannot score them.
An arbitrary-candidate PRIDICT2 adapter remains later work after mapping review.

The strict join uses full forward WT/ALT, PAM strand, oriented PAM coordinate,
native30, RTT and PBS. Matched pairs additionally agree on synthesized spacer,
nick-to-edit distance and homology length. PRIDICT2 introduced-allele lowercase
annotation is normalized only in its edited context and RTT. Experimental
construct metadata are outside this sequence identity key.

Matches include 541 SNP, 5 deletion, 4 insertion and 3 equal-length replacement guides;
262 plus and 291 minus strand. [Per-target coverage](comparison/target_coverage.csv)
and [all 553 joined records](comparison/exact_matches.csv) retain source CSV lines.

## 4. Descriptive rank comparison

These comparisons use the same 553 matched guides, separately within each of the
100 targets. Each target has 2–9 guides. Targets are equally weighted when taking
the median; score maxima determine the top guide. There are no top-score ties
in the pairs shown below. All 100 within-target correlations were defined.

{rank_table}

Across all 1,000 PEGG guides (ten per target), PEGG2 versus OptiPrime had median
within-target Spearman **{all_pair.median_within_target_spearman:.3f}** and the same
top guide for **{all_pair.targets_with_equal_top_sets}/100 targets**. Baseline RF
comparisons, exact tie counts and each target's results are also saved in
[rank_summary.csv](comparison/rank_summary.csv) and
[within_target_ranks.csv](comparison/within_target_ranks.csv).

The predictions frequently rank candidates differently. This supplies a reason
to compare against measured sensor data; it does not identify which model is
more accurate. The subset already passed PEGG selection, most edits are
synthetic SNVs from a few loci, and models use different cellular/construct
contexts. No p-values, accuracy estimates or 20% efficiency thresholds are claimed.
The existing OptiPrime adapter also still supplies native 20 rather than PEGG's
synthesized G+19 spacer. Treat these rank comparisons as provisional until that
representation is reviewed; an exact target/guide join does not resolve it.

## 5. Package changes and their evidence

| Change | Evidence and limits |
| --- | --- |
| Preserve empty alleles and reject missing/non-string alleles | Actual corrected inference, archived exact transport and old/new score comparison. The invalid-input error branch has not been exercised by a unit test. |
| Assign scores with positional row IDs; reject invalid/duplicate returned IDs and flag missing/nonfinite values | Source contract review and all five real inference runs. Duplicate-index and failure branches have not been exercised by unit tests. |
| Check both extension and spacer for polyT | Source review found Python `pbs_rtt or proto` checked only the nonempty extension. Dataset analysis found {int(poly.annotation_changed.sum())} changed flags in this fixture ({int(expected_poly.sum())} positive), so this dataset does not exercise the missed-spacer case. |

No unit tests were added or run. Inference experiments and artifact analysis are
the execution evidence reported here. PEGG design was not rerun after the polyT
fix; benchmark candidates remain the frozen original fixture.
The four pre-existing sensor-length error flags remain in that fixture; these
scoring results do not establish that every candidate is ready for oligo synthesis.

## 6. Remaining lab decisions and team ownership

- **Kexin:** confirm synthesized spacer representation. The existing adapter
  supplies native 20; PEGG reports forced G+19, different for 871/1,000 guides.
  Confirm editor/scaffold/motif/group/time before scientific mapping changes.
- **Kexin:** provide sensor measurements with guide identity, target identity,
  replicate, cell line/editor and measurement time. Agree on ranking evaluation.
- **Sathvik:** owns integration and the shared benchmark workflow; use these
  results to choose the first lab-data run after those decisions.
- **Aditya:** prepare PRIDICT2 arbitrary-candidate feature mapping; the 447 missing
  candidates describe the required coverage beyond the native default designer.
- **Neeraj:** review the strict identity table and prepare sensor-data joins and
  within-target ranking analysis using stable target/guide identifiers.

## Files and provenance

- [Reproduction instructions](source/benchmarks/README.md).
- [Reviewable source patch](source/pegg-progress.patch), including new scripts.
- [Runtime source/weights/configuration hashes](runtime/provenance.json).
- [Comparison source/file hashes](comparison/provenance.json).
- [Prepared team update](TEAM_UPDATE.md).
- All raw selected inputs, transported inputs, predictions, exit codes and logs
  are in `runtime/`; all coverage, rank and historical-comparison tables are in
  `comparison/`. No new PRIDICT2 inference was run in this round.
- `SHA256SUMS.txt` records the dated evidence bundle. The earlier handoff was
  not rebuilt. The local branch is available for review; no remote merge occurred.
'''
    (out / 'README.md').write_text(report)
    update = f'''# PEGG team update — October 13

Following the cancelled meeting, I expanded the existing PEGG→OptiPrime run to
all 1,000 guides across our 100-target preparation fixture. Every guide scored.
The scorer took {single['scorer_wall_seconds']:.1f}s with one worker and
{parallel['scorer_wall_seconds']:.1f}s with two; this is one pilot run per setting.
Maximum score difference between those runs was {parallel_delta:.3g}.

There are 553 exact candidate matches to our existing PRIDICT2 outputs across
all 100 targets. The other 447 fall outside that native designer's default RHA or
nick-distance limits. Matched-guide rank comparisons are now in the report.

I also found a correction to the earlier 100-guide smoke: default CSV loading
turned empty indel alleles into NaN and changed two sequence crops. Those two
old predictions are superseded by the corrected scores. The package now reports
this input problem, preserves score assignment with duplicate DataFrame labels,
and checks spacer as well as extension for polyT.

Kexin: the main next decisions are still native 20 versus synthesized G+19
representation, experimental metadata, and a sensor-data sample. These results
measure execution and predicted-rank agreement; they do not establish accuracy.

Suggested next work remains: Sathvik integration/OptiPrime/shared benchmark,
Aditya arbitrary-candidate PRIDICT2 mapping, Neeraj identities/schema/sensor joins.

The full dated report, source patch, inputs, predictions and logs are included.
This is a prepared draft; it has not been sent.
'''
    (out / 'TEAM_UPDATE.md').write_text(update)


if __name__ == '__main__':
    main()
