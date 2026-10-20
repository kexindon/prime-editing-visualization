# PEGG team update — October 13

Following the cancelled meeting, I expanded the existing PEGG→OptiPrime run to
all 1,000 guides across our 100-target preparation fixture. Every guide scored.
The scorer took 73.7s with one worker and
40.5s with two; this is one pilot run per setting.
Maximum score difference between those runs was 4.54e-07.

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
