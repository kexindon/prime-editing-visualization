# OptiPrime input/output mapping for the October 12 check-in

This is source inspection and local execution preparation. Kexin has not validated the PEGG mapping. No new PEGG integration wrapper was written.

## Source and scorer choice

Source: [alvin-hsu/optiprime-src](https://github.com/alvin-hsu/optiprime-src/tree/475db8a43d2fb69b030a8a781e2d365306119cc5), commit `475db8a43d2fb69b030a8a781e2d365306119cc5`. The companion [optiprime-front](https://github.com/alvin-hsu/optiprime-front/tree/c643cdff2e56d789b57471c8d849f50ff4947078) is a web frontend and explicitly excludes the hosted scoring backend. It is not a replacement local scorer.

`PREDICT_PE.py` consumes existing pegRNAs; `DESIGN_PE.py run` searches PBS/RTT lengths and optional edit combinations. The README design command lacks the actual CLI's required `run` subcommand. `smoke` exits after imports and proves no prediction. The Dockerfiles refer to absent `graphs/pe_model_2.rx`, `rs3.patch`, and `sglearn.patch`; local execution uses the present `graphs/pe_model.rx`.

## Required pegRNA columns

All positions below are **zero-based, half-open Python slices** in a DNA context oriented to the PAM-containing strand. Convert minus-strand genomic context to this orientation before indexing.

| OptiPrime column | Required content | PEGG mapping to validate |
|---|---|---|
| `spacer` | Actual 5′→3′ spacer/transcribed guide, excluding scaffold. DNA T is converted to RNA U by `format_pe_df`. The ideal target protospacer is `full_unedited[4:24]`. | Actual `Protospacer` plus its explicit transcription convention; do not blindly equate a 20nt DNA target to the actual Pol III spacer. |
| `rtt` | Actual 5′→3′ reverse transcription template, as the pegRNA component, excluding PBS. Equivalent in the native designer to reverse complement of `full_edited[21:21+R]`, then T→U. | `RTT` already has the pegRNA orientation; verify it by reverse-complementing it against the edited DNA context. |
| `pbs` | Actual 5′→3′ PBS, excluding RTT. Native designer uses reverse complement of the last P bases of the first 17 protospacer bases: `full_unedited[21-P:21]` (native design P≤16). | `PBS` already has the pegRNA orientation; verify the reverse complement matches the reference strand immediately before the nick. |
| `full_unedited` | Uppercase reference DNA: 4 nt before 20nt protospacer; PAM begins at index24; nick boundary is21. Native design trims to length `25+R-delta`, where delta is edited minus unedited allele length. | Reconstruct PAM-oriented reference context with enough flanks; PEGG's `target_30mer` alone is insufficient for RTTs extending past index30. Preserve actual variant and candidate coordinates. |
| `full_edited` | Edited DNA in the same orientation and start coordinate. Native design trims to length `25+R`, retaining 4 nt beyond the RTT endpoint. Includes the intended edit plus any bystanders. | Derive from the allele actually encoded by that candidate's RTT; do not apply the requested mutation twice or omit silent bystanders. |

Supporting source: [DESIGN_PE.py lines113–142](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/DESIGN_PE.py#L113-L142), [rows generated at198–212](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/DESIGN_PE.py#L198-L212), [constants](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/scripts/pe/pe_constants.py#L1-L2), [format_pe_df](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/scripts/pe/pe_utils.py#L25-L76).

The final four bases matter: `POST_HOM_END=4` is used by the heteroduplex features. The source's assertion text says an extra post-homology base, but its constructed rows and feature code use four bases. RNA components are concatenated `spacer + scaffold + RTT + PBS + linker`. Motif type is a separate feature; the motif nucleotide sequence is not concatenated by `format_pe_df`.

`proto30` defaults to `full_unedited[:30]`; it must be 30 nt. An explicitly supplied `proto30` is supported for short contexts, but it must refer to the same target. Rule Set 3 overwrites positions25 and26 with G for its own feature, while the PAM features retain the actual PAM. The code models a four-base PAM context `full_unedited[24:28]`; the canonical NGG is the first three bases of that slice.

### Spacer transcription convention

Native `make_spacer(pol3=True)` uses the exact target20 if it begins G. If the preceding base is G it uses a 21nt spacer including that base. Otherwise it substitutes a lowercase `g` for the first target base, retaining 20nt total. Lowercase `g` marks a leading mismatch feature; uppercasing it loses that distinction. `pol3=False` uses the exact target20. The direct scorer preserves supplied spacer strings except T→U. Native `write_results` calls `make_spacer` without passing `pol3`, so its output spacer is not reliable for a `pol3=False` job without inspecting that inconsistency.

### Context checks before accepting a mapping

For each candidate, validate spacer/PAM orientation, reference and alternative identity, PBS complement, RTT complement, the 4nt upstream anchor and 4nt post-RTT context, and the edit delta. Reject ambiguous DNA bases: `seq_encoding` indexes the A/C/G/T vocabulary directly. Use separate fixtures for both strands, substitutions, insertions, deletions, PAM disruption, leading-G behavior, and bystanders. The graph has 50 synthesis steps and clips repetition lengths at50 (`rx_module.py`); candidates with RTT>50 require an explicit policy rather than treating the score as fully supported.

## Experimental metadata and defaults

| Column/setting | Base formatter default | Native DESIGN default | Native PREDICT behavior for `Liu_HEK293T_*_PE2.csv` |
|---|---|---|---|
| `scaffold_name` | `SpCas9_OG` | `OG_F+E` (`settings.scaffold`) | Overwritten to `BlpI_F+E` |
| `cas9_type` | `PE2-Cas9` | `PEmax-Cas9` | Overwritten to `PEmax-Cas9` |
| `cas9_pam` | `SpNGG` | `SpNGG` | Preserves supplied value |
| `rt_name` | `PE2-RT` | `PE2-RT` | Default when missing |
| `motif` | `tevoPreQ1` | `tevoPreQ1` | Overwritten to `tevoPreQ1` |
| `pe_type` | `PE2` | `PE2` | From filename field4: PE2 or PE4 |
| `pol3` | Feature routines use True if absent | True | Existing value preserved; absent means feature default True |
| `group` | `NO_GROUP` | `Liu_HeLa` via `settings.cell_type` | `Liu_HEK293T` from filename |
| `time` | 1.0 | 4.0 supplied directly to model | 3.0 days assigned by filename, then minus1→2.0 effective days |
| `linker` | empty | empty via formatter | Existing value preserved |

The same native PREDICT Liu parser uses5.0 days for HeLa, then subtracts1→4.0. Do not call scores in these different settings interchangeable. Source: [PREDICT preprocessing](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/PREDICT_PE.py#L39-L50), [filename metadata rules](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/scripts/pe/pe_datasets.py#L10-L75), [DESIGN settings](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/DESIGN_PE.py#L152-L181).

Native PREDICT only recognizes Liu, Schwank, Kim and YKim filename families. It mutates metadata according to each family, even when the input CSV already has some fields. In the bundled 13-row file, `scaffold=SpCas9_OG` is not the consumed `scaffold_name`, and `motif=none` is overwritten. Therefore that smoke run establishes execution; it does not establish scientifically appropriate metadata for those rows.

Present model weights carry group factors for `Kim_A549`, `Kim_DLD1`, `Kim_HCT116`, `Kim_HEK293T`, `Kim_HeLa`, `Kim_MDA-MB-231`, `Kim_NIH3T3`, `Liu_HEK293T`, `Liu_HeLa`, `Schwank_HEK293T`, `Schwank_K562`, and `Schwank_U2OS`. The list is verified from all45 public rate pickle files in `weight_metadata.json`. Unknown groups do not raise: their initialized group factors remain zero and only an INFO message records a missing factor. Cell-type one-hot encoding is separately inferred from group-name substrings; an unknown type receives an average encoding. HCT116 has group factors but is absent from the cell-type vocabulary, so that distinction also requires care. Source: [rate loader](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/reaction/rx_model.py#L342-L379), [cell type encoding](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/scripts/pe/pe_inputs.py#L335-L343).

`edited_frac`, `indel_frac`, `weight` default to0,0,1 to permit inference. They are observables/bookkeeping, not needed measured efficiencies. Rows with weight≤0 or missing required observables are removed. The native example's `PE2_editing`, `PE2_indel`, `PE3_editing`, and `PE3_indel` are not automatically mapped into those observables; no accuracy comparison should be inferred from a successful run.

## Outputs and score meaning

| Entry point/output | Interpretation |
|---|---|
| PREDICT `predictions.csv` | Saved row index, `model_0`…`model_4` and `mean_pred`. Output model_0 corresponds to first CLI weight directory (here model_1), not a source directory named model_0. Mean is arithmetic ensemble mean. |
| PREDICT `joined_df.csv` | Processed input rows and derived features, preserving the source identifiers but resetting index. Join outputs by preserved candidate identifier or explicit row key and check row count/order; do not rely on index across independently filtered tables. |
| PREDICT `gmats.npy` | Per-model per-row 55×55 continuous-time transition generator matrices, not additional probability scores. |
| DESIGN `full_results.txt.gz` | TSV columns: `pegRNA_name`, `edit_name`, `RTT_len`, `PBS_len`, `spacer`, `RTT`, `PBS`, `OptiPrime_score`. DNA T in output components. Ensemble mean only. |
| DESIGN `edit_mapping.txt.gz` | Maps edit-combination code to full unedited/edited sequences. |
| DESIGN `top1_by_editcombo.txt.gz` | Best scoring length combination per edit combination; ties may yield more than one row. |
| DESIGN `input.json` | Echoed per-job input and settings. |

Scores are dimensionless model predictions on a nominal0–1 scale. The reaction model exponentiates the rate matrix at the requested time, extracts the `unedited` and `edited` observables, and normalizes their sum to1. Thus the returned value is the normalized edited fraction among those modeled observables, not an indel probability, a confidence interval, or a calibrated guarantee for a new laboratory. Multiplying by100 converts the number to percentage units only. Source: [RxModule](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/reaction/rx_module.py#L65-L97), [PREDICT output](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/PREDICT_PE.py#L117-L147), [DESIGN output](https://github.com/alvin-hsu/optiprime-src/blob/475db8a43d2fb69b030a8a781e2d365306119cc5/DESIGN_PE.py#L338-L370).

## Feasibility of scoring PEGG candidates (Option1)

The public model supports existing-candidate scoring because the five required columns are explicit pegRNA and allele fields; no search is required. The practical gating work is validating exact orientation/geometry and experimental settings. Current PEGG3.0 already contains an OptiPrime adapter at commit `07eda419b2c6d9b2f07b69ea048529e596c77edd`. Executing or auditing that upstream adapter is distinct from authoring a new integration, but it still needs the requested biological mapping review.

Do not use the JSON designer to represent arbitrary already-selected PBS/RTT lengths: it performs its own search. For the common batch, count distinct target mutations separately from candidate rows, multiple PBS/RTT lengths, and silent edit alternatives. The native design example contains one mutation job; the bundled native PREDICT example contains13 distinct described mutations.
