# PRIDICT2 schema and model notes for Kexin review

Prepared from official repository commit `c133c35e205062e766bf515d1766169167e993b9`. This is a source map and a local setup exercise. No PEGG scoring wrapper has been implemented.

## Public input and candidate generation

The CLI accepts a CSV with `sequence_name` and `editseq`. Names must be unique. `editseq` is the full forward target context with a single parenthesized edit: `... (A/G) ...`, `... (+CTG) ...`, or `... (-TTC) ...` (without spaces). Slash indels with `-` are also parsed. Only A/C/G/T is accepted. Equal-length multi-base replacement is supported; mixed unequal-length replacement is rejected. Shared leading/trailing unchanged bases must sit outside the brackets. The README requests at least 100 bases on each side; the implementation checks `<99`, so it happens to admit 99. The fixture uses at least 100. [Parser](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict2_pegRNA_design.py#L66), [flank check](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict2_pegRNA_design.py#L526).

The CLI designs its own candidates on both strands with NGG PAMs, PBS lengths 7–15, RT overhang lengths 3–19, and a default nick-to-edit window of 25 bases. It has no public argument for a supplied list of predesigned spacers/PBS/RTTs. Its hard edit-length limit is 40 bases; source comments state that validation extends to 15 bases. These runs use unchanged default candidate settings, two model variants and fold 1 only. `--use_5folds` averages all five folds and was not exercised. [Defaults](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict2_pegRNA_design.py#L16).

The emitted `Spacer-Sequence` is a forced leading G followed by **19 genomic bases**, not a direct extraction of 20 genomic bases. `pegRNA` concatenates spacer, the fixed optimized scaffold, `RTrevcomp`, then `PBSrevcomp`. RTT/PBS guide strings are reverse complements of the candidate-oriented target sequence. A PEGG candidate match needs agreement about forced G, strand, RTT orientation, PBS orientation, nick coordinate and edit identity. Do not join scores by mutation name or spacer alone. [Candidate construction](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict2_pegRNA_design.py#L696).

## Internal scorer requirements

The internal `deeppridict(pegdataframe, models_lst_dict)` can score a supplied feature dataframe. This is a potential integration entry point, not a supported arbitrary-pegRNA CLI. Producing equivalent features for PEGG candidates remains work for a schema approved by Kexin. The public CLI passing an edit into PRIDICT2 would redesign candidates and therefore would not preserve an arbitrary PEGG candidate list. [Scorer](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict2_pegRNA_design.py#L315).

Its 23 selected raw columns are:

```text
wide_initial_target, wide_mutated_target, deepeditposition,
deepeditposition_lst, Correction_Type, Correction_Length,
protospacerlocation_only_initial, PBSlocation,
RT_initial_location, RT_mutated_location,
RToverhangmatches, RToverhanglength, RTlength, PBSlength,
RTmt, RToverhangmt, PBSmt, protospacermt, extensionmt,
original_base_mt, edited_base_mt, original_base_mt_nan, edited_base_mt_nan
```

`seq_id` is inserted internally. Location pairs and edit-position lists are converted to string list syntax (`[10, 29]`). The preprocessing path derives aligned initial/mutated sequences, correction one-hot categories, and `Correction_Length_effective` (number of actual mismatches for a replacement; indel length otherwise). It constructs 18 sequence-level features: six length/count values normalized by 50, seven Wallace melting temperatures normalized by 200, three correction-type indicators and two missing-allele-temperature indicators. MFE and optional additional features are disabled for this CLI. No global epigenetic or genomic-coordinate feature is supplied. [Feature normalization](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict/pridictv2/dataset.py#L107), [tensor features](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict/pridictv2/dataset.py#L383).

## Coordinates and sequence context

All positions below are local to each candidate-oriented sequence; none is a genomic coordinate. On `Rv`, full `Original_Sequence`, `Edited_Sequence`, `OriginalAllele` and `EditedAllele` are reverse-complemented relative to the supplied mutation string. `PRIDICT2_Format` stays in supplied forward orientation.

- The raw `wide_initial_target` and `wide_mutated_target` start ten bases before the 19-base genomic protospacer and are each sliced to 99 bases. The edited target uses lowercase introduced bases in CSV output.
- Location intervals are **zero-based, half-open**, before internal gap alignment. In default designs, genomic protospacer = `[10,29)`, nick/RT start = `26`, PBS = `[26-PBSlength,26)`, mutated RTT = `[26,26+RTlength)`.
- For replacement, initial RTT span equals mutated RTT span. For deletion, initial RTT end is longer by `Correction_Length`; for insertion, it is shorter by that amount.
- `deepeditposition` is the edit start in the 99-base raw target context; `Editing_Position = deepeditposition - 26`. The multi-replacement list contains each mismatching base position. Indels have one start position with length stored separately.
- The preprocessor then inserts gap symbols for indels and adjusts annotations; consumers should not pre-adjust raw intervals a second time. Tensor construction crops with a 20-base context setting and a 150-base internal span guard, so the raw 99-base columns are not an assertion that every tensor input is always exactly 99 positions.

[Raw coordinates](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict2_pegRNA_design.py#L737), [alignment/annotation](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict/pridictv2/data_preprocess.py#L43), [tensor window](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict/pridictv2/dataset.py#L220).

## Outputs and architecture

[columns.csv](columns.csv) lists **all 51 observed output headers**, their meanings, source lines at the pinned commit, and a literal value from the successful original example. Each row represents a candidate pegRNA, not a distinct mutation. The standard CSV exports intended-edit scores for HEK and K562, local rank per mutation, and percentile against the bundled Library-Diverse score table. The closest reference score supplies the percentile; percentile is not another efficiency percentage. The internal model also predicts unedited and indel outcomes, but the public candidate CSV does not expose those two probabilities. [Aggregation/export](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict2_pegRNA_design.py#L913).

The checkpoint configuration confirms two layers of bidirectional GRU encoders, nucleotide and region annotation embeddings, cell-specific local/global attention and sequence-feature embedders, and cell-specific three-outcome distribution decoders. Initial and mutated sequences use separate encoders. Prediction runs in evaluation mode, float32, with gradients disabled. Final intended-edit scores average the two model variants and multiply the edited-outcome fraction by 100. The saved fold-1 configuration and actual versions are in [environment_and_models.json](environment_and_models.json). The implementation chooses CUDA if present, otherwise CPU; this Mac run used CPU (no CUDA), one Torch thread, and no MPS override. [Model construction](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict/pridictv2/predict_outcomedistrib.py#L147).

The published model combines two training routes, both fine-tuned using HEK293T and K562 Library-Diverse data. HEK represents a mismatch-repair-deficient context and K562 a proficient one. The paper recommends choosing that context according to the experimental mismatch-repair setting. Chromatin effects require the separate ePRIDICT model. Local score agreement is not evidence of biological accuracy on these synthetic variants. [Mathis et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7617539/)

The protocol notes that PRIDICT2 does not model the extra nicking guide used by PE3/PE5 because these guides were absent from training. The current repository's optional nicking path imports `DeepCas9_TestCode`, which is absent in this checkout; it was not exercised. Primer generation was also not exercised. [Protocol](https://pmc.ncbi.nlm.nih.gov/articles/PMC7619000/), [lazy import](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict2_pegRNA_design.py#L259).

## Decisions for Kexin before a wrapper

1. Should PEGG's candidates be scored unchanged, or should PRIDICT2 redesign candidates for the same mutations?
2. Which shared schema fixes spacer leading G, strand, nick coordinate, PBS/RTT orientations, equal-length replacement intervals, indel gaps and feature calculations?
3. Which context score should the lab use (HEK/MMR-deficient or K562/MMR-proficient), and should production use five folds?
4. Which edit lengths, PBS/RTT ranges, editor/scaffold choices and variant types should the wrapper accept or reject?
5. Which known native PRIDICT2 examples must give equivalent features and scores before scoring arbitrary PEGG candidates?

These are review items, not claims that cross-tool scoring has already been validated.
