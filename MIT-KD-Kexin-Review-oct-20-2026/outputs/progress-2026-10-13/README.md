# MIT KD — PEGG progress, October 13, 2026

## Completed this round

- Added reproducible runtime and exact-candidate comparison scripts on local
  branch `codex/pegg-model-prep` (base PEGG commit `07eda419`).
- Scored all **1,000 existing PEGG guides across 100 targets**, with all five
  OptiPrime models. All five experiment processes exited 0; all returned scores
  were finite and candidate order was preserved.
- Joined **553/1,000 guides (55.3%)**
  exactly to saved PRIDICT2 predictions, covering **all 100 targets** with 2–9 matches
  each and no duplicate or ambiguous identities.
- Corrected three engineering issues: blank-allele handling, score assignment by
  row position, and polyT annotation. Model definitions and weights are unchanged.

**Status:** execution, scaling and candidate-comparison evidence are ready for
review. Kexin's mapping/metadata decision and measured sensor data are still
needed before interpreting these predictions as better experimental designs.

## 1. Runtime and numerical agreement

| Guides | Workers | Scorer call (s) | Whole process (s) | Finite scores | Largest child RSS (GB) |
| --- | --- | --- | --- | --- | --- |
| 100 | 1 | 16.60 | 29.15 | 100/100 | 0.991 |
| 250 | 1 | 23.96 | 25.09 | 250/250 | 0.978 |
| 500 | 1 | 41.73 | 42.96 | 500/500 | 1.000 |
| 1000 | 1 | 73.74 | 74.94 | 1000/1000 | 0.957 |
| 1000 | 2 | 40.50 | 41.70 | 1000/1000 | 0.944 |

These are single observations on this Mac (arm64, macOS 27, 14 logical CPUs).
Fixed settings: CPU, numerical thread environment caps of 1, batch size 128, `Liu_HeLa`,
time 4, all five models. Candidate prefixes are nested and all cover 100 targets.
They are candidate counts, not counts of additional mutations.

The scorer call includes model worker startup, feature generation, inference and
output transport. Whole-process time also includes parent imports and archiving
the inputs. Time outside the scorer call was larger in the first process, so
100 guides took longer than 250 in whole-process time. Cache state and one repetition per configuration
prevent a firm general speed claim. These times exclude PEGG/PRIDICT2 design.

For 1,000 guides, two workers gave a **1.82× observed scorer speedup**.
The largest absolute prediction difference from the single-worker result was
**4.54e-07**; no difference exceeded 1e-6. The 100/250/500-guide results were
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

| Guide | Edit | Historical score | Corrected score |
| --- | --- | --- | --- |
| pegRNA_0 | DEL | 0.016486328 | 0.100762998 |
| pegRNA_20 | INS | 0.002337609 | 0.139366888 |

For `pegRNA_0` (deletion1), the unedited crop was 62 nt instead of 65 nt. For
`pegRNA_20` (insertion1), it was 55 nt instead of 52 nt. Edited crops were 55 nt in both.
The corrected experiment preserves empty strings with `keep_default_na=False`.
The package now rejects non-string alleles with a per-row explanation. The
maximum historical score change among the other 98 guides was
2.62e-07.

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

| Score pair | Median within-target Spearman | Same top guide /100 targets |
| --- | --- | --- |
| PEGG2_Score / OptiPrime_Score | 0.000 | 6 |
| PEGG2_Score / PRIDICT2_HEK | 0.200 | 19 |
| PEGG2_Score / PRIDICT2_K562 | 0.300 | 25 |
| OptiPrime_Score / PRIDICT2_HEK | 0.800 | 48 |
| OptiPrime_Score / PRIDICT2_K562 | 0.779 | 35 |

Across all 1,000 PEGG guides (ten per target), PEGG2 versus OptiPrime had median
within-target Spearman **0.055** and the same
top guide for **1/100 targets**. Baseline RF
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
| Check both extension and spacer for polyT | Source review found Python `pbs_rtt or proto` checked only the nonempty extension. Dataset analysis found 0 changed flags in this fixture (85 positive), so this dataset does not exercise the missed-spacer case. |

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
