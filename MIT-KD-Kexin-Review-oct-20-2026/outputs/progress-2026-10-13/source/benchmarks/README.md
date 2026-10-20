# Existing-adapter runtime and candidate comparison

These scripts measure the existing OptiPrime scorer and compare predictions for
exactly matched guide sequences. They are analysis tools; there are no measured
editing outcomes in the preparation fixture.

## Reproduce

Run from the workspace containing `work/` and `outputs/`, using the previously
prepared PEGG Python environment. Both output directories must be new. Paths can
be replaced with another checkout and the same source fixture.

```sh
work/pegg/.venv/bin/python work/pegg/PEGG3.0/benchmarks/optiprime_scaling.py \
  --pegg-repo work/pegg/PEGG3.0 \
  --optiprime-repo work/optiprime/repo \
  --optiprime-python work/optiprime/.venv/bin/python \
  --candidates outputs/meeting-prep/pegg_100/pegg_candidates.csv \
  --output-dir outputs/progress-rerun/runtime

work/pegg/.venv/bin/python work/pegg/PEGG3.0/benchmarks/compare_candidate_scores.py \
  --runtime-dir outputs/progress-rerun/runtime \
  --pridict-dir outputs/pridict2 \
  --mutation-lookup outputs/meeting-prep/pegg_100/mutation_lookup.csv \
  --historical-scores outputs/meeting-prep/optiprime_100/scored_candidates.csv \
  --output-dir outputs/progress-rerun/comparison
```

### Important input detail

Load PEGG CSVs with `pd.read_csv(path, keep_default_na=False)`. Empty `REF` or
`ALT` strings are valid indel alleles. Default pandas NA conversion caused two
incorrect context crops in the historical 100-guide smoke. The scorer now reports
non-string alleles as row errors. Do not replace arbitrary missing data with an
empty allele without knowing what it represents.

## Runtime experiment

- Sort existing candidates by PEGG rank, mutation index, and guide ID. Take nested
  prefixes of 100, 250, 500 and 1,000 rows; each prefix covers all 100 fixture targets.
- Run each prefix with one worker, and 1,000 rows again with two workers.
- Use all five released weight directories, group `Liu_HeLa`, time 4 and model
  batch size 128. These settings reproduce the existing adapter assumptions.
- Start each configuration in a fresh PEGG Python process. OptiPrime imports,
  featurization, model loading/JIT and prediction occur afresh in its subprocesses.
- Run configurations sequentially. Numerical-library thread environment variables
  are capped at one and JAX is set to CPU. This does not measure actual active
  thread counts or enforce operating-system CPU affinity.
- Save selected candidates, exact adapter inputs, scores/errors, command arguments,
  stdout/stderr logs, exit codes, source diffs and hashes, and weight hashes.
- Use one observation per configuration. Filesystem and OS caches may be warm;
  these runs do not establish a stable speed ratio or extrapolate to new loci.

`scorer_wall_seconds` covers the call to `optiprime.score`, including input
conversion, worker startup, featurization, inference, output parsing and row
assignment. It excludes parent imports and the separate input-archiving pass.
`process_wall_seconds` includes the whole fresh PEGG process. Neither number
includes the earlier PEGG design or PRIDICT2 design runs.

Memory is reported as the parent's RSS high-water mark and the largest reaped
child's RSS high-water mark. The latter is **not the sum of simultaneous workers**;
it cannot establish the memory required to run two workers together. macOS returns
RSS in bytes; Linux units are converted from KiB.

## Exact candidate join

`candidate_identity.py` builds a tuple of:

1. Full WT in input-forward orientation.
2. Full ALT in input-forward orientation.
3. PAM strand.
4. PAM start in the PAM-oriented full WT.
5. Native 30-base context.
6. Guide-oriented RTT.
7. Guide-oriented PBS.

For PEGG, PAM start is `PAM_start` on plus and `len(WT)-PAM_start` on minus.
For PRIDICT2, locate the unique `wide_initial_target` in `Original_Sequence`; PAM
start is that position plus 29. Native30 is `oriented_WT[p-24:p+6]`. Full reverse
strand PRIDICT2 contexts are reverse-complemented back to input orientation.
Only PRIDICT2 `Edited_Sequence` and `RTrevcomp` are uppercased to remove its
introduced-allele annotation. OptiPrime spacer case is never normalized.

Duplicate identities, ambiguous anchors, target identity mismatches, and
inconsistent PAM/spacer/geometry raise an error. All matched pairs must also
have identical synthesized spacers, nick-to-edit distances and homology lengths.
The hash identifies these guide and target sequences, **not** scaffold, editor,
motif or all other experimental metadata. The scripts reject silent-bystander
designs; this fixture has none. Unmatched-guide reason classification is specific
to the saved PRIDICT2 default search (PBS7–15, RHA3–19, nick distance at most25),
not a general claim about what the predictor can score internally.

Source basis: PEGG `prime.py:971–1046` at `07eda419b2c6d9b2f07b69ea048529e596c77edd`;
PRIDICT2 `pridict2_pegRNA_design.py:619–740` and defaults at lines20–28 at
`c133c35e205062e766bf515d1766169167e993b9`. Line numbers refer to those commits.

## Ranking comparisons

The analysis calculates average-rank Spearman correlations separately within each
mutation, then reports their median/mean with equal weight per target. Constant
scores and groups with fewer than two candidates are recorded as undefined.
Exact top-score ties are retained; both intersection and equality of top sets are
reported. There are no p-values or biological accuracy statistics.

Baseline/OptiPrime comparisons use all 1,000 PEGG candidates and also the matched
subset. PRIDICT2 comparisons use only the matched subset. This avoids comparing
the rank among ten PEGG candidates to the rank among hundreds of PRIDICT2 designs.

The same guide scores across nested batches and worker counts are compared by
their sequence identity. Differences above 1e-6 are counted descriptively; this
threshold is not a validation of biological equivalence.

## Limits still requiring lab review

- The existing adapter supplies native target20; PEGG synthesized spacers use G+19.
- Existing fixed editor/scaffold/motif/group/time settings need Kexin's review.
- PRIDICT2 HEK/K562 and OptiPrime HeLa predictions have different model contexts.
- The fixture has four upstream examples plus 96 synthetic SNVs at a few loci.
- The 1,000 candidates were already selected by PEGG. This does not evaluate the
  full design pool, `prime.run`'s OptiPrime prefilter, clinical variants or sensor
  accuracy. Lab sensor data are needed for ranking performance comparisons.

No unit-test suite was added or run for the package changes. Actual inference
experiments and source inspection are reported separately in the dated report.
