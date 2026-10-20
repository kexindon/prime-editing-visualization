# Notes on the supplied PEGG slides

Read on October 9, 2026, for the October 12 check-in.

Source: `2026-09-14 PEGG3.0 Update (2).key`, supplied by Sathvik from Downloads.
The Keynote document contains four active slides with no speaker notes. Text and
slide order were extracted from the document archive; the full-resolution images
embedded in each slide were inspected. Keynote itself was not used.

| Slide | Content | Relevance to the preparation work |
|---|---|---|
| 1 | Prime editing components: protospacer, PBS, RTT, editor, and repair of the edited flap. The schematic labels PBS lengths 8–15 nt and RTT lengths 10–74 nt. | Distinguish the guide's individual sequence components from a complete expression construct. These labels are diagram values, not a verified contract for either prediction model. |
| 2 | Mutation inputs and design parameters enter PEGG; candidate generation is followed by ranking, filtering, oligo construction, and visualizations. | The agreed first route keeps PEGG's candidate generation and replaces or supplements the score at the ranking stage. |
| 3 | pegRNA expression and sensor architecture, with optional tevoPreQ1 motif, adapters, and restriction sites. | A sensor or synthesis oligo is not equivalent to the full genomic context required by a prediction model. Motif and editor configuration must match the interpretation of a model score. |
| 4 | PAM search; iteration over PBS/RTT lengths; ranking features including homology, GC, polyT, restriction sites, and on/off-target properties. Mutation inputs include chromosome coordinates or ClinVar IDs. | Record both the candidates generated and the subset actually scored. Candidate coverage can change through filtering even when model inference succeeds. |

The deck is an overview. It does not specify PRIDICT2/OptiPrime input columns,
index conventions, model weights, or a 100-mutation sample dataset. Those details
are being established from the versioned source code and execution results.

The repository named PEGG3.0 currently declares package version 3.1.0 and already
contains an OptiPrime scoring implementation. This is a source-code observation
at commit `07eda419b2c6d9b2f07b69ea048529e596c77edd`, not a claim made by the slides.
