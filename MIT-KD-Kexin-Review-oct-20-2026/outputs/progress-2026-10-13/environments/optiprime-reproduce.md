# Reproduce the OptiPrime preparation runs

Working directory: `/Users/sathvikloke/Documents/Codex/2026-10-05/mit-x20`.

The isolated environment is `work/optiprime/.venv`; public source is `work/optiprime/repo`. The successful native PREDICT run scores all13 supplied mutations using all5 models. It is an execution check, not an external accuracy evaluation. **The filename parser overwrites the example's scaffold and motif to `BlpI_F+E` and `tevoPreQ1`, and sets `Liu_HEK293T` with effective time2.0 days.** Preserve that caveat when presenting these results.

## Re-run in the prepared environment

The logging harness only launches the native CLI and caps CPU library thread counts. It does not transform sequences, implement a new scoring adapter, or modify model architecture. Execute from the workspace:

```sh
work/optiprime/.venv/bin/python outputs/optiprime/run_logged.py predict13_repeat \
  "$PWD/work/optiprime/.venv/bin/python" PREDICT_PE.py \
  --data_dir ../inputs/native13 --rx_graph graphs/pe_model.rx \
  --weight_dirs weights/model_1 weights/model_2 weights/model_3 weights/model_4 weights/model_5 \
  --batch_size 13
```

The native PREDICT CLI writes a new timestamped directory beneath `work/optiprime/inputs/native13`. The harness writes stdout/stderr plus command, environment, exit code and elapsed time to `work/optiprime/logs`.

The native README design example uses this command (choose a fresh output directory; native DESIGN skips an already existing job directory):

```sh
work/optiprime/.venv/bin/python outputs/optiprime/run_logged.py design_col7a1_repeat \
  "$PWD/work/optiprime/.venv/bin/python" DESIGN_PE.py run \
  --run_data COL7A1_R185X.json --graph_rx graphs/pe_model.rx \
  --weight_dirs weights/model_1 weights/model_2 weights/model_3 weights/model_4 weights/model_5 \
  --out_path "$PWD/outputs/optiprime/native/design_repeat"
```

## Installation performed

The stock `pip install -r requirements.txt` equivalent failed when RS3's required LightGBM3.3.5 built from source and tried writing `~/LightGBM_compilation.log`. Exact error is retained in `logs/install.log`. We retained version3.3.5 and relocated that build log to its own source directory. Existing Homebrew CMake and libomp were used; nothing was installed globally. The compiler used upstream make's `-j4`; prediction processes use one thread for BLAS/OpenMP/TensorFlow and serial feature extraction.

`uv venv --python 3.11` selected existing `/opt/homebrew/opt/python@3.11/bin/python3.11`, version3.11.15. The requested `/opt/homebrew/bin/python3` is version3.14.7; it was used for file/provenance work rather than the model environment.

Core versions are JAX/jaxlib0.4.35, Flax0.10.2, Optax0.2.4, NumPy1.26.4, Pandas2.2.3, scikit-learn1.5.2, LightGBM3.3.5, RS3 0.0.16, ViennaRNA2.7.2 and TensorFlow2.18.1. Exact81-package freeze is `environment_freeze.txt`. `uv pip check` reports all installed packages compatible. The JAX/NumPy/Pandas pins also match the stack recommended in current PEGG3.0's OptiPrime notebook; the notebook's conflicting claim that this JAX stack needs NumPy≥2 was not followed.

Installation commands, shown from the workspace:

```sh
export UV_CACHE_DIR="$PWD/work/optiprime/cache/uv"
export UV_PYTHON_INSTALL_DIR="$PWD/work/optiprime/python"
/Users/sathvikloke/.local/bin/uv venv --python 3.11 work/optiprime/.venv
```

For a clean setup, fetch the LightGBM3.3.5 sdist described in `work/optiprime/lightgbm_source.json` (SHA256 in `source_checksums.json`), extract it to `work/optiprime/lightgbm-3.3.5`, and apply `lightgbm_build_log.patch` to its `setup.py`. Then:

```sh
CMAKE_POLICY_VERSION_MINIMUM=3.5 CMAKE_BUILD_PARALLEL_LEVEL=2 \
  /Users/sathvikloke/.local/bin/uv pip install --python work/optiprime/.venv/bin/python \
  work/optiprime/lightgbm-3.3.5
/Users/sathvikloke/.local/bin/uv pip install --python work/optiprime/.venv/bin/python \
  -r work/optiprime/repo/requirements.txt -c outputs/optiprime/compatibility-constraints.txt
/Users/sathvikloke/.local/bin/uv pip install --python work/optiprime/.venv/bin/python \
  tensorflow==2.18.1 -c outputs/optiprime/compatibility-constraints.txt
```

These commands reproduce the chosen constraints, while the full freeze records the exact resolved versions. To lock all transitive versions on a fresh setup after the local LightGBM build, use the full freeze, whose LightGBM line points to the local source path.

Apply `compatibility.patch` inside the source clone before running: it changes only the default `RxDataset` feature processing mode to serial. The source otherwise starts `cpu_count()-1` workers and changes process-global ViennaRNA RNA/DNA parameters in the parent; serial execution respects the intended pre/process/post sequence and avoids excess workers. No weights, architecture, forward pass, input sequence, or score formulas were changed. This is documented local execution, not an unmodified upstream CLI claim.

TensorFlow is an additional dependency absent from native `requirements.txt`: PREDICT's loader imports it at inference time. DESIGN does not need that loader. PE3's XGBoost secondary stage and twinPE were not exercised here.

## Evidence

- `provenance.json` and `source_checksums.json`: source commit, hardware, input/model hashes and sizes.
- `environment_freeze.txt`, `compatibility-constraints.txt`, `compatibility.patch`, `lightgbm_build_log.patch`: environment and exact changes.
- `logs/*.json` and `logs/*.log`: raw execution metadata and stdout/stderr.
- `native/predict13/`: all native prediction columns, joined inputs, transition matrices, annotated rows and verification counts.
- `input_output_mapping.md`: source-backed orientation, fields, defaults, groups, score units and review questions.

No input was taken from the lab, and no external biological validation is claimed. Distinct mutations and pegRNA length candidates are counted separately.
