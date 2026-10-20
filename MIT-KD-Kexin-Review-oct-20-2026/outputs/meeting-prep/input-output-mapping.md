# PEGG, PRIDICT2 and OptiPrime: mapping for Kexin's review

Prepared October 9, 2026, for the October 12 check-in. **Draft for biological and
software-contract review; successful execution does not approve a mapping.**

## Versions and what is already implemented

| Tool | Inspected source revision | Relevant entry point |
|---|---|---|
| [PEGG3.0 repository](https://github.com/kexindon/PEGG3.0/tree/07eda419b2c6d9b2f07b69ea048529e596c77edd) | `07eda419`, package **3.1.0** | `prime.run`; existing `optiprime.score` |
| [PRIDICT2](https://github.com/uzh-dqbm-cmi/PRIDICT2/tree/c133c35e205062e766bf515d1766169167e993b9) | `c133c35` | Public design CLI; internal `deeppridict` scorer |
| [OptiPrime](https://github.com/alvin-hsu/optiprime-src/tree/475db8a43d2fb69b030a8a781e2d365306119cc5) | `475db8a` | `PREDICT_PE.py`; `DESIGN_PE.py run` |

The current PEGG repository already includes an OptiPrime adapter, example
notebook, and selectable scoring flags. This preparation inventories and executes
that existing code. No new PEGG/model integration wrapper was written. The local
model changes only permit serial execution in this environment; patches and logs
are saved with each model's deliverables.

## Inputs

| Tool / route | Input names | Meaning and requirements |
|---|---|---|
| PEGG `PrimeDesign` | `SEQ` | Full context with one `LEFT(REF/ALT)RIGHT` edit; insertion `(/ALT)`, deletion `(REF/)`. No genome download needed. |
| PEGG `WT_ALT` | `WT`, `ALT` | Full unedited and edited contexts; PEGG infers the edit through alignment. Suitable for the current simple setup fixture; complex indels need explicit checks. |
| PEGG `cBioPortal` | `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Type` | MAF-like records plus matching reference genome. Substitution/deletion coordinates are 1-based inclusive. Insertions use the source's flanking-boundary convention; validate on a supplied lab example before conversion. |
| PRIDICT2 public batch CLI | `sequence_name`, `editseq` | Unique target name and full sequence with `(REF/ALT)`, `(+INSERTION)`, or `(-DELETION)`; use at least 100 nt on each side. It designs and scores its own candidates. |
| OptiPrime direct candidate scorer | `spacer`, `rtt`, `pbs`, `full_unedited`, `full_edited` | Supplied pegRNA components plus allele contexts oriented to the PAM strand. Metadata, model group, and time also affect interpretation. |
| OptiPrime designer | JSON `name`, `unedited`, `edit_segments`, optional settings | Searches PBS/RTT lengths and optional edit combinations. Multiple candidates or silent-edit combinations are not multiple independent target mutations. |

PEGG's `WT_ALT`/`PrimeDesign` formatter drops custom input IDs and reconstructs
`mutation_idx` from row order. `pegRNA_id` is also regenerated at the end of each
run. Keep an external mutation lookup for this fixture; a future approved schema
needs durable target and candidate IDs. cBioPortal rows with a reference mismatch
can be dropped before `mutation_idx` is assigned.

Sources: [PEGG input parsing](https://github.com/kexindon/PEGG3.0/blob/07eda419b2c6d9b2f07b69ea048529e596c77edd/pegg/prime.py#L243),
[PEGG normalization](https://github.com/kexindon/PEGG3.0/blob/07eda419b2c6d9b2f07b69ea048529e596c77edd/pegg/prime.py#L1373),
[PRIDICT2 parsing](https://github.com/uzh-dqbm-cmi/PRIDICT2/blob/c133c35e205062e766bf515d1766169167e993b9/pridict2_pegRNA_design.py#L66),
[OptiPrime formatting](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/scripts/pe/pe_utils.py#L25).

## Candidate sequence fields

| PEGG output | Meaning | External-model implication |
|---|---|---|
| `wt_w_context`, `alt_w_context` | Full input-oriented reference and intended-edit contexts | `alt_w_context` excludes added bystanders; reconstruct the actual candidate allele when bystanders are present. |
| `PAM_strand` | `+` / `-` relative to input sequence | Orient both contexts to the candidate's PAM strand before using model-local positions. |
| `PAM_start`, `PAM_end` | Local boundaries in the forward input context | For `+`, interval is `[start,end)`; for `-`, endpoints are descending and interval is `[end,start)`. They are not genomic positions. |
| `Protospacer` | At the default size, forced uppercase G + 19 native bases; used in PEGG oligos | May differ from the native target20. This is a material issue for OptiPrime's spacer features. |
| `Protospacer_30` | Native PAM-oriented 4 nt flank + 20 nt target + 3 nt PAM + 3 nt flank | Native target20 is `[4:24]`; keep this separate from the synthesized spacer. |
| `RTT`, `PBS`, `RTT_PBS` | Guide-oriented DNA strings; extension is RTT followed by PBS | Complement checks should use the corresponding target strand. Do not reverse-complement them twice. |
| `Distance_to_nick`, `RHA_size` | Bases from nick to edit start, and homology after the edited allele | Used to reconstruct indels and fully edited candidate contexts. |
| `has_silent_bystander`, `bystander_positions` | Whether extra synonymous changes were generated; positions are offsets into reverse-complement(RTT) | Preserve candidate-specific changes. A declared reading frame alone does not establish a true biological coding context. |
| `sensor_wt`, `sensor_alt` | Cropped sensor constructs, with their own orientation | These are not universal full-context scorer inputs. Check `sensor_error` separately from candidate generation. |

Full observed PEGG field definitions: [pegg-columns.csv](pegg-columns.csv).
Source: [candidate construction](https://github.com/kexindon/PEGG3.0/blob/07eda419b2c6d9b2f07b69ea048529e596c77edd/pegg/prime.py#L971),
[edited-context reconstruction](https://github.com/kexindon/PEGG3.0/blob/07eda419b2c6d9b2f07b69ea048529e596c77edd/pegg/prime.py#L1085).

## Model-specific requirements

**OptiPrime:** reference context starts 4 nt upstream of the native20, PAM starts
at index24, and the nick boundary is21. Native design retains four bases after
the RTT endpoint: edited length `25 + RTT_length`, unedited length
`25 + RTT_length - allele_length_delta`. The existing PEGG adapter reconstructs
these contexts and supplies RTT/PBS unchanged. It uses `Liu_HeLa`, time4.0 by
default and fixes scaffold/editor/motif settings in its worker. Those settings
must be matched to the lab's constructs.

**Specific spacer discrepancy for review:** the existing adapter supplies native20
from `Protospacer_30[4:24]`, while PEGG synthesizes its `Protospacer` value. OptiPrime
uses the supplied spacer for pegRNA features and a lowercase leading `g` to mark
a forced-G mismatch. `pol3=True` does not repair that discrepancy. In this fixture,
871 of 1,000 returned PEGG candidates differ between the displayed spacer and
native20; the 100-candidate execution trial includes 86 such rows. Finite output
therefore establishes execution, not correct representation of those oligos.
Simply using an uppercase `Protospacer` would still omit the lowercase mismatch
feature. This needs a deliberate mapping decision before any correction.

Sources: [existing adapter](https://github.com/kexindon/PEGG3.0/blob/07eda419b2c6d9b2f07b69ea048529e596c77edd/pegg/optiprime.py#L253),
[OptiPrime spacer features](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/scripts/pe/pe_inputs.py#L219),
[native spacer construction](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/DESIGN_PE.py#L112).

**PRIDICT2:** the public CLI redesigns candidates; it cannot accept an arbitrary
PEGG candidate list. Option1 would need the internal scorer and equivalent raw
features. Its 23 selected fields include 99-nt initial/mutated target windows,
edit type/length/positions, PBS/RT intervals, sequence lengths/counts, and Wallace
melting temperatures. Default raw windows use protospacer interval `[10,29)` for
19 native bases and nick boundary26. Its model preprocessing handles indel gaps;
these coordinates must not be pre-adjusted twice. The public output spacer is
forced G + 19 native bases. Full detail and all 51 output columns are in
[the PRIDICT2 mapping](../pridict2/mapping.md) and [column dictionary](../pridict2/columns.csv).

The OptiPrime [detailed mapping](../optiprime/input_output_mapping.md) also records
native PREDICT's filename-driven metadata changes. Its supplied 13-target run uses
different experimental settings from the existing PEGG adapter. The resulting
numbers should not be treated as a same-condition comparison.

## Score columns and interpretation

| Column | Meaning |
|---|---|
| PEGG `PEGG2_Score` | Fixed weighted regression sum of 10 features; arbitrary model units, not a probability. |
| PEGG `RF_Score` | Raw random-forest prediction. Calling code does not establish its calibration or target units. |
| PEGG `OnTarget_Azimuth_Score` | Cas9 activity score scaled by100; it is not prime-editing efficiency. |
| PRIDICT2 `PRIDICT2_0_editing_Score_deep_HEK` / `..._K562` | Predicted intended-edit fraction multiplied by100 for the named model context. Current runs use fold1 and average two model variants. |
| PRIDICT2 `HEK_rank` / `K562_rank` | Descending rank within that mutation's own generated candidate set. |
| PRIDICT2 `HEK_percentile_to_librarydiverse` / `K562_percentile_to_librarydiverse` | Percentile relative to the supplied reference score distribution, not editing percentage. |
| OptiPrime native `model_0`…`model_4`, `mean_pred` | Individual model predictions and arithmetic ensemble mean; nominal fraction0–1 at the specified metadata/time. |
| Existing PEGG `OptiPrime_Score`, `OptiPrime_error` | Ensemble prediction or missing score with a per-row reason. A finite score does not certify mapping correctness. |

The lab's approximately20% **measured sensor editing** cutoff must not be assumed
to equal an unvalidated0.20 OptiPrime threshold. No such cutoff was applied here.

## Decisions requested before new integration code

1. Confirm the current PEGG3.1 snapshot is the intended baseline and whether the
   existing OptiPrime module should be reviewed/corrected as part of our work.
2. Agree the synthesized-spacer versus native-target convention, including
   lowercase forced-G mismatch representation in OptiPrime.
3. Confirm editor, scaffold, tevoPreQ1, model group and effective time for the lab.
   The weights contain more groups than the two contexts mentioned in the notes;
   available parameters alone do not establish validation for every context.
4. Confirm PRIDICT2 Option1 should preserve PEGG candidates through the internal
   feature/scorer path. Public-CLI redesign is a separate route.
5. Set supported edit types/lengths and PBS/RTT ranges, bystander policy, and
   stable target/candidate IDs. Include native-example feature parity before use.
6. Confirm the appropriate measured sensor dataset and evaluation split for the
   later ranking benchmark; the current fixture has no experimental outcomes.

One additional source discrepancy was found: `contains_polyT_terminator` currently
checks `RTT_PBS` for normal rows because the expression is `motif in (pbs_rtt or
proto)`; it does not check both strings as its comment suggests. This was documented
without changing PEGG behavior. Four of the 1,000 setup rows also have a reported
sensor-length error, which is separate from their successful guide generation.
