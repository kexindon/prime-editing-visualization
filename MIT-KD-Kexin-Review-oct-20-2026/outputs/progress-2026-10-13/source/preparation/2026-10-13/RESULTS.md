# October 13 progress summary

The five existing-adapter inference experiments all exited 0 with finite scores:
100, 250, 500 and 1,000 candidates with one worker, and 1,000 with two workers.
The 1,000-guide scorer call took 73.74 s serial and 40.50 s with two workers; maximum
score difference was 4.54e-7. These are single pilot observations on the existing
100-target fixture, not population-level throughput or accuracy measurements.

Strict sequence/coordinate matching found 553 exact overlaps with 55,845 saved
PRIDICT2 rows, spanning all 100 targets. The 447 absent designs are explained by
PRIDICT2 default RHA/nick-distance constraints. Median within-target Spearman
between OptiPrime and PRIDICT2 HEK is 0.800 on that subset; measured outcomes are
still needed to evaluate ranking accuracy.

Source changes fix blank-allele handling, positional score assignment and polyT
annotation. The historical 100-guide smoke had two incorrect indel crops; its
scores for pegRNA_0 and pegRNA_20 are superseded. New CSV reads preserve empty
strings. No unit tests were added/run; real inference, artifact analysis and
source inspection are the evidence collected this round.

The full report and portable dated ZIP are in workspace
`outputs/progress-2026-10-13/`. The earlier sent handoff remains unchanged.
[Benchmark instructions](../../benchmarks/README.md) describe reproduction.

Kexin's review remains necessary for native 20 versus synthesized G+19 and the
experimental metadata. Lab sensor measurements remain necessary for accuracy
claims and selecting a scoring default.
