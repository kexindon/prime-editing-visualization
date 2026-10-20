#!/usr/bin/env python3
"""Bundle the completed October13 evidence without rewriting historical results."""
import argparse
import hashlib
import json
import shutil
import subprocess
import zipfile
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workspace', type=Path, required=True)
    ap.add_argument('--uv', type=Path, required=True)
    args = ap.parse_args()
    root = args.workspace.resolve()
    repo = root / 'work/pegg/PEGG3.0'
    out = root / 'outputs/progress-2026-10-13'
    source = out / 'source'
    source.mkdir(parents=True, exist_ok=True)
    for part in ('benchmarks', 'preparation/2026-10-13'):
        for path in (repo / part).rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts:
                destination = source / path.relative_to(repo)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, destination)
    patch = subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=repo)
    for path in sorted(source.rglob('*')):
        if not path.is_file() or path.suffix not in ('.py', '.md'):
            continue
        relative = str(path.relative_to(source))
        diff = subprocess.run(['git', 'diff', '--no-index', '--binary', '--', '/dev/null', relative],
                              cwd=repo, capture_output=True)
        if diff.returncode != 1:
            raise RuntimeError('Could not construct addition patch for ' + relative)
        patch += diff.stdout
    (source / 'pegg-progress.patch').write_bytes(patch)
    env_dir = out / 'environments'
    env_dir.mkdir(exist_ok=True)
    for name, interpreter in (('pegg', root / 'work/pegg/.venv/bin/python'),
                              ('optiprime', root / 'work/optiprime/.venv/bin/python')):
        frozen = subprocess.check_output([str(args.uv), 'pip', 'freeze', '--python', str(interpreter)])
        (env_dir / (name + '-freeze.txt')).write_bytes(frozen)
    for name in ('compatibility-constraints.txt', 'compatibility.patch', 'lightgbm_build_log.patch', 'reproduce.md'):
        shutil.copy2(root / 'outputs/optiprime' / name, env_dir / ('optiprime-' + name))
    shutil.copy2(root / 'outputs/meeting-prep/CORRECTION-2026-10-13.md', out / 'HISTORICAL_CORRECTION.md')
    # In this bundle the correction note is next to the current report.
    correction = out / 'HISTORICAL_CORRECTION.md'
    correction.write_text(correction.read_text().replace('../progress-2026-10-13/README.md', 'README.md'))
    excluded = {'SHA256SUMS.txt', 'bundle.json', 'MIT-KD-PEGG-Progress-2026-10-13.zip'}
    files = sorted(p for p in out.rglob('*') if p.is_file() and p.name not in excluded
                   and '__pycache__' not in p.parts and 'mpl-cache' not in p.parts)
    manifest = out / 'SHA256SUMS.txt'
    manifest.write_text(''.join(sha(p) + '  ' + str(p.relative_to(out)) + '\n' for p in files))
    archive = out / 'MIT-KD-PEGG-Progress-2026-10-13.zip'
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zipped:
        for path in files + [manifest]:
            zipped.write(path, str(Path('MIT-KD-PEGG-Progress-2026-10-13') / path.relative_to(out)))
    summary = {'archive': archive.name, 'sha256': sha(archive), 'bytes': archive.stat().st_size,
               'files_in_archive': len(files) + 1, 'manifest_sha256': sha(manifest),
               'scope': 'New October13 outputs and source only; historical raw inputs/PRIDICT2 outputs remain in the earlier handoff.'}
    (out / 'bundle.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
