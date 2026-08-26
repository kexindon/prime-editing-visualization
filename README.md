# prime-editing-visualization

Interactive web tool for visualizing prime editing guide RNA (pegRNA) design.

Paste a double-stranded DNA sequence and step through a design: set the reading
frame, pick a PAM and spacer, select the bases to modify, generate the PBS and
RTT, and compare the edited sequence and protein against the reference — with
optional silent bystander mutations.

All of the biology is done by [PEGG](https://github.com/samgould2/PEGG-2.0)
itself (`pegg.prime`, `pegg.bystander`). This repository is a visual front end
for it, not a reimplementation: PAM search, pegRNA geometry and silent bystander
enumeration are all delegated to PEGG so the two cannot drift apart.

## Features

**1. Sequence & reading frame** — paste DNA (FASTA headers, whitespace and
numbering are stripped). Choose the frame offset (0/1/2) and the strand to
translate; the protein is shown and each amino acid is drawn over its codon.

**2. PAM & spacer** — search any IUPAC PAM (`NGG`, `NG`, `NNGRRT`, …) on both
strands. Every hit lists its protospacer and nick site; click one to select it.

**3. Nucleotides to modify** — click and drag across the sequence map to select
bases, then type the replacement. Substitutions, insertions (zero-width
selection), deletions (empty replacement) and indels are all supported, and the
variant type is classified by PEGG.

**4. pegRNA** — generates the protospacer, PBS, RTT, the 3′ extension and the
full pegRNA, with the distance to the nick, the right homology arm, and whether
the edit disrupts the PAM or the protospacer. Everything is drawn onto the
sequence map in place.

**5. Silent bystanders** — optionally add synonymous mutations near the edit
(MMR evasion, and PAM knockout where they land in the PAM). Each option is
verified synonymous against the declared reading frame, listed with its base
change, mutation count, distance to the edit and whether it kills the PAM, and
highlighted on the map when selected.

**6. RTT comparison** — paste any RTT (or send one over from a design) and see
exactly what it installs: the edited sequence, the nucleotide differences, and
the resulting protein compared with the reference.

## Running it

```bash
pip install -r requirements.txt
python run.py
```

Then open <http://127.0.0.1:5000>.

PEGG must be importable. If it is not pip-installed, keep a `PEGG3.0` checkout
next to this repository (the default), or point at one explicitly:

```bash
PEGG_PATH=/path/to/PEGG3.0 python run.py
```

`/api/health` reports which copy of PEGG is actually in use.

## API

All endpoints take and return JSON. Offsets are 0-based and half-open, on the
forward strand of the sequence supplied.

| Endpoint | Purpose |
|---|---|
| `POST /api/translate` | Reading frame and protein for a sequence |
| `POST /api/pams` | All PAMs on both strands, with protospacers and nick sites |
| `POST /api/design` | pegRNA design(s) for an edit, with silent bystanders |
| `POST /api/apply_rtt` | Apply an RTT to the reference and report the result |
| `POST /api/compare` | Compare two sequences at nucleotide and protein level |
| `GET /api/health` | Liveness, and the PEGG install in use |

## Notes and limitations

**Silent bystanders require a forward-strand transcript.** PEGG's `orf` frame
mode declares the reading frame as an offset into the sequence you paste, which
cannot express a transcript running the other way. Paste the coding strand.
pegRNAs whose *PAM* is on the reverse strand are designed normally either way —
it is the transcript, not the PAM, that has to be forward. Reverse-strand
transcripts need genomic coordinates (PEGG's `cds` mode), which this tool does
not take.

**The sequence is assumed to be entirely coding** when bystanders are switched
on, as PEGG's `orf` mode requires. A sequence containing intronic or
untranslated bases would yield bystanders that are only apparently silent. CDS
membership and splice-site distance cannot be checked without coordinates.

**Designs are not scored.** PEGG's ranking models (`PEGG2_Score`, RF) operate on
a full library dataframe; this tool shows the design space for one edit rather
than ranking it. Use `pegg.prime.run()` for scored, library-scale design.

## Correctness

The bridge is checked against PEGG by round-trip: every generated RTT is applied
back onto the reference and must reproduce the intended edit exactly, and every
silent bystander is independently re-translated and must leave the protein
unchanged. Across randomized sequences, edit positions, variant types (SNP, ONP,
insertion, deletion, indel) and all three frames, on both PAM strands:

- 3147/3147 RTTs round-trip to the exact intended product
- 42430/42430 silent bystanders verified synonymous
