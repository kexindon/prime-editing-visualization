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

**[Full user guide →](DOCUMENTATION.md)** — install, what to paste, and a
step-by-step walkthrough.

## Features

**1. Sequence & reading frame** — paste DNA (FASTA headers, whitespace and
numbering are stripped). Choose the frame offset (0/1/2) and the strand to
translate; the protein is shown and each amino acid is drawn over its codon.
Residues are numbered from 1 by default, or from **First amino acid is #** when
the pasted sequence is a fragment — set it to 262 for an exon starting at G262
and an edit there is reported as `R273C` rather than `R12C`.

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

**Getting the sequence in** — `cds_for_visualizer.ipynb` turns a human or mouse
gene and a coding mutation coordinate into the exon to paste, in frame and on
the coding strand, and the position to select. See
[Getting a sequence to paste](#getting-a-sequence-to-paste).

## Running it

PEGG needs Python 3.9 or 3.10 -- its on-target scoring models are pickled with
scikit-learn 1.1.1, which does not build on newer Pythons.

```bash
conda create -n pegg python=3.9 -y
conda activate pegg
pip install -r requirements.txt
python run.py
```

Then open <http://127.0.0.1:5050>.

Already have the environment? Just `conda activate pegg && python run.py`, or
`pip install --upgrade pegg` first to pick up a new PEGG release.

Port 5050 rather than Flask's usual 5000, because macOS runs its AirPlay
Receiver on 5000 and Flask fails to bind there. Override with `PORT`, and use
`HOST=0.0.0.0` to let others on the network reach it:

```bash
PORT=8080 python run.py
HOST=0.0.0.0 python run.py      # reachable from other machines
```

PEGG must be importable. If it is not pip-installed, keep a `PEGG3.0` checkout
next to this repository (the default), or point at one explicitly:

```bash
PEGG_PATH=/path/to/PEGG3.0 python run.py
```

`/api/health` reports which copy of PEGG is actually in use, the Python version,
and whether cyvcf2 had to be stubbed out (see requirements.txt).

## Getting a sequence to paste

The tool searches for PAMs and builds the RTT **inside the string you paste**,
so that string has to be real genomic DNA. To go from a gene and a mutation
coordinate to that string, run
[`cds_for_visualizer.ipynb`](cds_for_visualizer.ipynb).

For a **coding variant** in a human or mouse gene it returns the **single exon
containing the mutation**, taken from the canonical transcript (the same curated
table H2M and PEGG use), in the transcript's direction and strand, trimmed at
both ends to whole codons so it is in frame from offset 0.

One exon is also exactly the region a pegRNA can reach: an RTT is one contiguous
stretch of genomic DNA, so it cannot cross an intron.

A variant is given either as a protein change or as a genomic coordinate:

```python
exon = get_exon('TP53', 'R175H')                      # protein change
exon = get_exon('TP53', 7577121, ref='G', alt='A')    # R273C hotspot, GRCh37
# TP53  ENST00000269305.4  - strand   exon 4 of 10, chr17:7577019-7577155
# 135 nt (2 trimmed to codon boundaries), 45 codons
# amino acids 262-306 of TP53 (G262 ... R306)
#
# paste the sequence; mutation at position 34 (1-based), reference base 'C'
#   that is R273 in the full protein (codon 12 of this exon)
#   expected 'C' on the transcript strand: matches
#   type 'T' into "Replace with"
#   room for bystanders in this exon: 33 nt 5', 101 nt 3'

print(exon['sequence'])
```

A protein change (`R175H`, `p.Arg175His`, `R306*`) is resolved against the
canonical transcript — the same transcript the numbering refers to — and the
reference residue is checked before anything else, so a build mismatch or a
wrong transcript is caught rather than silently designed against. Only
single-base routes are reported: where a substitution needs two or three bases
changed, that is said outright rather than guessed at, and where several single
-base routes exist they are all listed.

The returned dict carries the same values for programmatic use: `sequence`,
`protein`, `first_aa` / `last_aa` (the exon's span in the full protein),
`mutated_aa`, `position`, `genomic_position`, `ref_transcript` /
`alt_transcript`.

Note the allele flip: TP53 is on the `-` strand, so the genomic `G>A` is `C>T`
on the sequence you paste. `get_exon()` does that conversion and checks the
reference base agrees, which also catches a build mismatch.

Then, in the web tool:

1. Paste `exon['sequence']` into **Sequence & reading frame**.
2. Leave **Reading frame offset** on `Frame 0` and **Translate strand** on
   `Forward (+)` — the sequence is already in frame and already on the coding
   strand, which is the point of extracting it this way.
   Set **First amino acid is #** to the first number the notebook printed (262
   in the example above), so residues are numbered as they are in the gene.
3. In **Nucleotides to modify**, enter the reported position as both *Start* and
   *End* (both 1-based), and type `alt_transcript` into *Replace with*.
4. Pick a PAM in step 2 and design. Silent bystanders can be switched on in
   step 4; *Transcript strand* there is fixed at `Forward (+)`, which the
   extracted sequence already satisfies.

**Coding variants only.** An intronic or UTR position has no reading frame to be
silent against, and `get_exon()` raises rather than guessing. A mutation inside
a codon split across an intron is rejected for the same reason: that codon
cannot be edited from a single RTT.

The notebook needs a gffutils database and a reference FASTA of the **same
build**; paths are set in one cell near the top, with human GRCh37/GRCh38 and
mouse GRCm39 entries.

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

Both of these are what `cds_for_visualizer.ipynb` exists to handle: it returns a
single exon, reverse-complemented for a `-` strand gene and trimmed to whole
codons, so the sequence it hands you satisfies them **and** is genuine genomic
DNA — a spliced transcript would not be, and guides designed across a junction
would target sequence that does not exist.

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
