# PEGG progress — October 13, 2026

The October 12 meeting was cancelled. The user requested continued work on the
existing PEGG project. This implements the already discussed runtime/coverage
work and focused engineering repairs. Work state is tracked in claude-mem's
`PEGG October13 progress` list.

## Design and scope

Three possible routes were considered: change the biological mapping now,
wait for the lab meeting, or make progress with the existing scorer and explicit
limitations. The third route supplies reviewable software and evidence while
preserving the team's agreement that Kexin reviews new mapping code.

Keep the original 100-target, 1,000-candidate fixture fixed. Call the existing
OptiPrime scorer under its existing sequence conversion and metadata. Repair
polyT annotation and assignment of returned scores to DataFrame row positions.
Source review also uncovered historical CSV conversion of empty alleles to NaN;
preserve strings in the new experiment and reject non-string alleles in the
adapter. No feature definition or weights are changed. Add a standalone benchmark script
under `benchmarks/` and analyze strict candidate overlap with the existing
PRIDICT2 exports. No PRIDICT2 wrapper or new biological mapping is implemented.

## Files and responsibilities

- `pegg/prime.py`: check the extension and spacer separately for polyT.
- `pegg/optiprime.py`: assign predictions by positional row IDs, validate returned
  row IDs, report missing/non-finite model outputs per row, and reject non-string
  alleles rather than interpreting NaN as three bases.
- `benchmarks/candidate_identity.py`: strict no-bystander sequence/geometry keys.
- `benchmarks/optiprime_scaling.py`: orchestrate isolated benchmark processes,
  preserve original candidates and IDs, and record raw outputs/metrics/logs.
- `benchmarks/compare_candidate_scores.py`: join only verified exact candidate
  identities, summarize coverage and descriptive ranking agreement.
- `benchmarks/README.md`: reproducible commands and interpretation limits.
- `preparation/2026-10-13/`: reviewable plan and final progress summary.
- Workspace `outputs/progress-2026-10-13/`: raw experimental evidence and reports.

## Execution sequence

1. Apply the two narrowly scoped fixes. Preserve other uncommitted preparation.
   Do not add or run unit tests; the user requested continued project work, and
   the session instructions require an explicit request for implementation tests.
2. Implement an isolated-process run mode and orchestration mode. Sort candidates
   by existing rank then target so nested prefixes cover all targets early.
   Score 100, 250, 500 and 1,000 candidates serially; repeat 1,000 with two worker
   processes. Use batch size 128, group `Liu_HeLa`, time 4, CPU and one thread per
   numerical library. No concurrent benchmark runs.
3. Save candidate hashes, source commit/diff hashes, returned IDs, scores/errors,
   elapsed whole-process and scorer-call time, parent/child peak memory and exit
   status. Preserve model errors as failures. No silent partial success or run
   directory overwrite.
4. Analyze complete supplied sequence/context/strand/RTT/PBS identities to join
   the pre-existing PRIDICT2 predictions; flag ambiguities rather than selecting
   an arbitrary candidate. Report excluded rows and matched target counts.
5. Compare baseline, OptiPrime and PRIDICT2 ranks only within the observed matched
   subset. Report ties and undefined correlations. No inferred biological accuracy.
6. Produce a dated report, save a source patch, and record remaining questions.

## Interpretation constraints

These are synthetic/public setup inputs from a few loci. Neither lab sensor data
nor a mapping approval has arrived. The existing adapter's native20 versus
synthesized-spacer discrepancy remains present and will be counted. Results
cannot establish better pegRNAs, calibrated thresholds or whole-proteome/library
performance. Model settings differ across tools. One repetition per configuration
is a pilot measurement, not a statistically established speed comparison.
Process RSS maxima are not aggregate memory across concurrent workers. Native
PRIDICT2 candidates and PEGG candidates are counted separately.
