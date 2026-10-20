#!/usr/bin/env python3
"""Measure the existing PEGG OptiPrime adapter on fixed nested guide sets.

Each configuration runs in a fresh PEGG process. This is an inference experiment,
not a biological validation or a unit-test suite. Existing run directories are
never overwritten. See README.md for scope and timing definitions.
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

from candidate_identity import identity_hash, pegg_key


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def git_snapshot(repo, destination, name):
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    patch = subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=repo)
    path = destination / (name + '-tracked.patch')
    path.write_bytes(patch)
    return {'commit': commit, 'tracked_diff_sha256': digest(path),
            'status': subprocess.check_output(['git', 'status', '--short'], cwd=repo, text=True),
            'tracked_patch': path.name}


def run_one(args):
    # Preserve the interpreter symlink: resolving .venv/bin/python would bypass
    # the environment when passed onward to the OptiPrime subprocess.
    sys.path.insert(0, str(args.pegg_repo.resolve()))
    import numpy as np
    import pandas as pd
    from pegg import optiprime

    rows = pd.read_csv(args.candidates, keep_default_na=False)
    if args.count > len(rows) or args.count <= 0:
        raise ValueError('Requested count outside available candidate range')
    rows = rows.sort_values(['pegRNA_rank', 'mutation_idx', 'pegRNA_id'], kind='stable')
    chosen = rows.head(args.count).reset_index(drop=True).copy()
    chosen['candidate_identity_sha256'] = [identity_hash(pegg_key(r)) for _, r in chosen.iterrows()]
    if chosen.candidate_identity_sha256.duplicated().any():
        raise ValueError('Duplicate candidate identities in selected set')
    chosen.to_csv(args.output_dir / 'selected_candidates.csv', index=False)
    # Archive the exact adapter transport before inference; useful for reviewing
    # crop lengths, spacer representation and indel handling without rerunning.
    transported, errors = optiprime._rows_from_pegdf(chosen)
    columns = ['spacer', 'rtt', 'pbs', 'full_unedited', 'full_edited', 'proto30']
    mapped = pd.DataFrame([r if r is not None else [None] * 6 for r in transported], columns=columns)
    mapped.insert(0, 'candidate_identity_sha256', chosen.candidate_identity_sha256)
    mapped['mapping_error'] = errors
    mapped.to_csv(args.output_dir / 'adapter_inputs.csv', index=False)
    start = time.perf_counter()
    scored = optiprime.score(
        chosen, optiprime_src=str(args.optiprime_repo.resolve()),
        optiprime_python=str(args.optiprime_python.absolute()), group='Liu_HeLa',
        time=4.0, batch_size=128, n_jobs=args.workers, quiet=False,
    )
    elapsed = time.perf_counter() - start
    scored.to_csv(args.output_dir / 'scored_candidates.csv', index=False)
    values = pd.to_numeric(scored.OptiPrime_Score, errors='coerce')
    good = np.isfinite(values)
    same_order = scored.candidate_identity_sha256.tolist() == chosen.candidate_identity_sha256.tolist()
    summary = {
        'candidate_count': len(chosen), 'target_count': int(chosen.mutation_idx.nunique()),
        'requested_workers': args.workers,
        'actual_workers': max(1, min(args.workers, (sum(r is not None for r in transported) + optiprime._MIN_SHARD - 1) // optiprime._MIN_SHARD)),
        'batch_size': 128, 'group': 'Liu_HeLa', 'model_time': 4.0,
        'scorer_wall_seconds': elapsed,
        'parent_peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024),
        'largest_reaped_child_peak_rss_bytes': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * (1 if sys.platform == 'darwin' else 1024),
        'memory_scope': 'Separate parent and largest child process high-water marks; not aggregate concurrent RSS.',
        'returned_rows': len(scored), 'finite_scores': int(good.sum()),
        'targets_with_finite_scores': int(scored.loc[good, 'mutation_idx'].nunique()),
        'preserved_candidate_order': same_order,
        'errors': scored.OptiPrime_error.fillna('').value_counts().to_dict(),
        'score_min': float(values[good].min()) if good.any() else None,
        'score_max': float(values[good].max()) if good.any() else None,
        'displayed_vs_adapter_spacer_mismatches': int(chosen.Protospacer.ne(chosen.Protospacer_30.str.slice(4, 24)).sum()),
        'selected_sha256': digest(args.output_dir / 'selected_candidates.csv'),
        'scored_sha256': digest(args.output_dir / 'scored_candidates.csv'),
        'adapter_inputs_sha256': digest(args.output_dir / 'adapter_inputs.csv'),
    }
    save_json(args.output_dir / 'summary.json', summary)
    print(json.dumps(summary, indent=2), flush=True)
    if len(scored) != len(chosen) or not good.all() or not same_order:
        raise RuntimeError('Incomplete inference; inspect scores, errors and logs')


def orchestrate(args):
    args.output_dir.mkdir(parents=True, exist_ok=False)
    scripts_dir = args.output_dir / 'scripts'
    scripts_dir.mkdir()
    for source in (Path(__file__), Path(__file__).with_name('candidate_identity.py')):
        (scripts_dir / source.name).write_bytes(source.read_bytes())
    env = dict(os.environ)
    caps = {k: '1' for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                            'NUMEXPR_NUM_THREADS', 'TF_NUM_INTRAOP_THREADS', 'TF_NUM_INTEROP_THREADS')}
    caps.update({'JAX_PLATFORMS': 'cpu', 'MPLBACKEND': 'Agg',
                 'XLA_FLAGS': '--xla_force_host_platform_device_count=1'})
    env.update(caps)
    env['MPLCONFIGDIR'] = str(args.output_dir.resolve() / 'mpl-cache')
    metadata = {
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'platform': platform.platform(), 'machine': platform.machine(), 'cpu_count': os.cpu_count(),
        'pegg_python': sys.executable, 'pegg_python_version': sys.version,
        'optiprime_python': str(args.optiprime_python.absolute()),
        'input_path': str(args.candidates.resolve()), 'input_sha256': digest(args.candidates),
        'configuration': {'group': 'Liu_HeLa', 'model_time': 4.0, 'batch_size': 128},
        'environment_overrides': caps,
        'pegg_source': git_snapshot(args.pegg_repo, args.output_dir, 'pegg'),
        'optiprime_source': git_snapshot(args.optiprime_repo, args.output_dir, 'optiprime'),
        'weight_files': {str(p.relative_to(args.optiprime_repo)): digest(p)
                         for p in sorted((args.optiprime_repo / 'weights').rglob('*')) if p.is_file()},
        'script_sha256': {p.name: digest(p) for p in scripts_dir.iterdir()},
        'interpretation': 'Single observations on a synthetic/public few-locus fixture; mapping unapproved; no accuracy or library-wide performance claim.',
    }
    save_json(args.output_dir / 'provenance.json', metadata)
    records = []
    for count, workers in ((100, 1), (250, 1), (500, 1), (1000, 1), (1000, 2)):
        name = 'n%d_w%d' % (count, workers)
        run_dir = args.output_dir / name
        run_dir.mkdir()
        command = [sys.executable, str(Path(__file__).absolute()), '--single',
                   '--pegg-repo', str(args.pegg_repo.resolve()),
                   '--optiprime-repo', str(args.optiprime_repo.resolve()),
                   '--optiprime-python', str(args.optiprime_python.absolute()),
                   '--candidates', str(args.candidates.resolve()), '--output-dir', str(run_dir.resolve()),
                   '--count', str(count), '--workers', str(workers)]
        save_json(run_dir / 'command.json', command)
        print('Starting ' + name, flush=True)
        start = time.perf_counter()
        with (run_dir / 'process.log').open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
        elapsed = time.perf_counter() - start
        record = {'run': name, 'exit_code': result.returncode, 'process_wall_seconds': elapsed}
        summary_path = run_dir / 'summary.json'
        if summary_path.exists():
            record.update(json.loads(summary_path.read_text()))
        save_json(run_dir / 'exit.json', {'exit_code': result.returncode, 'process_wall_seconds': elapsed})
        records.append(record)
        save_json(args.output_dir / 'runs.json', records)
        print('%s exit=%d process_seconds=%.3f' % (name, result.returncode, elapsed), flush=True)
        if result.returncode:
            raise RuntimeError('Inference failed in ' + name + '; raw log preserved')
    scalar_keys = [k for k, v in records[-1].items() if not isinstance(v, (dict, list))]
    with (args.output_dir / 'runtime.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=scalar_keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pegg-repo', type=Path, required=True)
    ap.add_argument('--optiprime-repo', type=Path, required=True)
    ap.add_argument('--optiprime-python', type=Path, required=True)
    ap.add_argument('--candidates', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--single', action='store_true', help=argparse.SUPPRESS)
    ap.add_argument('--count', type=int, default=100)
    ap.add_argument('--workers', type=int, default=1)
    args = ap.parse_args()
    try:
        run_one(args) if args.single else orchestrate(args)
    except Exception as exc:
        if args.single and args.output_dir.exists():
            save_json(args.output_dir / 'failure.json', {'error_type': type(exc).__name__, 'error': str(exc)})
        traceback.print_exc()
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
