# Correction recorded October 13, 2026

The October9 100-guide OptiPrime smoke used default pandas CSV loading. Empty
REF/ALT alleles became NaN and the existing adapter counted `str(NaN)` as three
bases in its context crop calculation. Two predictions were affected:

| Guide | Target | Historical score | Corrected score |
| --- | --- | ---: | ---: |
| pegRNA_0 | deletion1 |0.016486328|0.100762998|
| pegRNA_20 | insertion1 |0.002337609|0.139366888|

The historical run still demonstrates execution; these two scores are superseded.
The corrected experiment uses `keep_default_na=False`; the package now rejects
non-string alleles rather than misinterpreting missing values. For all new runs,
use the scripts and exact model inputs in the [October13 report](../progress-2026-10-13/README.md).

Original raw files, the historical reproduction script and the October11 sent
ZIP are preserved. The old reproduction script should not be reused unchanged.
