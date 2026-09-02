"""
Bridge between the web tool and PEGG.

Every design decision here is delegated to PEGG itself rather than being
reimplemented: PAM search, protospacer/PBS/RTT geometry and silent bystander
enumeration all call into pegg.prime / pegg.bystander. This module's job is only
to (a) build the inputs PEGG expects, (b) unpack its outputs into offsets the
browser can highlight, and (c) keep the two in the same coordinate system.

Coordinate system
------------------
The browser only ever deals with offsets into the forward (top) strand of the
sequence the user pasted. PEGG works in the PAM strand's orientation, which for
a '-' orientation pegRNA is the reverse complement. Every offset returned by
this module has been mapped back to forward-strand coordinates by to_forward();
nothing downstream needs to know which strand a design came from.
"""

import os
import sys

import Bio.Seq

#Import PEGG, preferring a sibling checkout over any pip-installed copy so the
#tool runs against the version of PEGG being developed alongside it. This is not
#merely a convenience: silent bystanders live in pegg.bystander, which older
#released versions of PEGG do not carry at all, and an installed copy would
#otherwise shadow the checkout and import cleanly right up until that point.
_HERE = os.path.dirname(os.path.abspath(__file__))
_PEGG_DEV = os.environ.get('PEGG_PATH') or os.path.abspath(
    os.path.join(_HERE, '..', '..', 'PEGG3.0'))

if os.path.isdir(os.path.join(_PEGG_DEV, 'pegg')):
    sys.path.insert(0, _PEGG_DEV)


def _stub_cyvcf2_if_broken():
    """
    pegg.prime imports cyvcf2 at module level, but uses it in exactly one
    function -- clinvar_VCF_translator(), for reading ClinVar VCFs -- which this
    tool never calls. cyvcf2 is a compiled extension linked against htslib, and
    a mismatched install raises ImportError on a missing symbol rather than
    simply being absent, which would take the whole app down over a feature it
    does not use.

    So: try the real thing first, and only if it is unusable put a stub in
    sys.modules so that PEGG's import succeeds. The stub's VCF raises if it is
    ever actually called, so this can silently disable VCF reading but can never
    silently return wrong data.
    """
    try:
        import cyvcf2  # noqa: F401
        return False
    except Exception:
        pass

    import types

    stub = types.ModuleType('cyvcf2')

    def _unavailable(*args, **kwargs):
        raise ImportError(
            'cyvcf2 is not usable in this environment, so PEGG cannot read VCF '
            'files here. This does not affect pegRNA design; only '
            'pegg.prime.clinvar_VCF_translator() needs it.')

    stub.VCF = _unavailable
    sys.modules['cyvcf2'] = stub
    return True


CYVCF2_STUBBED = _stub_cyvcf2_if_broken()

from pegg import prime, bystander


PEGG_SOURCE = os.path.dirname(prime.__file__)

#PEGG builds protospacers as 'G' + the proto_size nt 5' of the PAM.
DEFAULT_PROTO_SIZE = 19


class DesignError(ValueError):
    """Raised for input the user can fix, and reported as a 400 rather than a 500."""


#--- sequence helpers ---------

def revcomp(seq):
    return str(Bio.Seq.Seq(seq).reverse_complement())


def clean_sequence(seq):
    """
    Strips whitespace, digits and FASTA headers out of a pasted sequence and
    upper-cases it. Rejects anything that is not then a plain DNA string, since
    a silently-dropped character would shift every offset downstream of it.
    """
    if seq is None:
        raise DesignError('No sequence supplied.')

    lines = [i for i in str(seq).splitlines() if not i.startswith('>')]
    cleaned = ''.join(lines)
    cleaned = ''.join(i for i in cleaned if not i.isspace() and not i.isdigit())
    cleaned = cleaned.upper().replace('U', 'T')

    if len(cleaned) == 0:
        raise DesignError('No sequence supplied.')

    bad = sorted(set(i for i in cleaned if i not in 'ACGTN'))
    if bad:
        raise DesignError('Sequence contains non-DNA characters: %s' % ', '.join(bad))

    return cleaned


def translate(seq):
    """
    Translates a nucleotide string, ignoring any trailing 1-2 nt that do not
    complete a codon. Returns '' for anything shorter than one codon.
    """
    usable = len(seq) - (len(seq) % 3)
    if usable < 3:
        return ''
    return str(Bio.Seq.Seq(seq[:usable]).translate())


def reading_frame(seq, orf_start, strand='+'):
    """
    Translates `seq` in the frame beginning at `orf_start`, on the given strand,
    and returns both the protein and the codon layout the browser draws.

    Each codon carries its offsets in FORWARD-strand coordinates, so that a
    '-' strand reading frame highlights the same physical bases as a '+' one.

    Parameters
    -----------
    seq
        *type = str*

        Forward-strand sequence.

    orf_start
        *type = int*

        0, 1 or 2; offset at which the first codon begins, measured along the
        strand being translated.

    strand
        *type = str*

        '+' or '-'. On '-' the reverse complement is translated, and orf_start
        is measured from the 3' end of the forward sequence.
    """
    if orf_start not in (0, 1, 2):
        raise DesignError('Reading frame offset must be 0, 1 or 2.')

    work = seq if strand == '+' else revcomp(seq)
    L = len(seq)

    codons = []
    for i in range(orf_start, len(work) - 2, 3):
        codon = work[i:i + 3]
        aa = bystander.translate_codon(codon)

        #map the codon's three bases back onto the forward strand
        if strand == '+':
            positions = [i, i + 1, i + 2]
        else:
            positions = [L - 1 - i, L - 2 - i, L - 3 - i]

        codons.append({
            'codon': codon,
            'aa': aa if aa is not None else 'X',
            'index': len(codons),
            'start': min(positions),
            'end': max(positions) + 1,
            'positions': sorted(positions),
        })

    return {
        'strand': strand,
        'orf_start': orf_start,
        'protein': ''.join(i['aa'] for i in codons),
        'codons': codons,
    }


#--- PAM search ---------

#IUPAC codes PEGG's PAM_finder understands. Validated here because PAM_finder
#indexes its lookup table directly and raises a bare KeyError on anything else.
IUPAC = set('AGCTNRYSWKMBDHV')


def clean_pam(pam):
    """Validates and normalises a PAM string."""
    pam = (pam or '').strip().upper()
    if not pam:
        raise DesignError('No PAM sequence supplied.')
    bad = sorted(set(i for i in pam if i not in IUPAC))
    if bad:
        raise DesignError(
            'PAM contains characters that are not IUPAC codes: %s' % ', '.join(bad))
    return pam


def find_pams(seq, pam='NGG', proto_size=DEFAULT_PROTO_SIZE):
    """
    Finds every PAM on both strands using PEGG's own PAM_finder, and reports
    each one with its protospacer and nick site in forward-strand coordinates.

    A PAM is only returned when the full protospacer fits inside the supplied
    sequence; a truncated protospacer is not a real design choice.
    """
    seq = clean_sequence(seq)
    pam = clean_pam(pam)
    rc = revcomp(seq)
    L = len(seq)

    out = []
    for strand, work in (('+', seq), ('-', rc)):
        for start, end in prime.PAM_finder(work, pam):
            #protospacer occupies [start - proto_size, start) on `work`
            proto_start = start - proto_size
            if proto_start < 0:
                continue

            #PEGG nicks 3 nt 5' of the PAM, i.e. between work[start-3] and work[start-4]
            nick = start - 3

            if strand == '+':
                f_pam_start, f_pam_end = start, end
                f_proto_start, f_proto_end = proto_start, start
                f_nick = nick
            else:
                f_pam_start, f_pam_end = L - end, L - start
                f_proto_start, f_proto_end = L - start, L - proto_start
                f_nick = L - nick

            out.append({
                'id': '%s:%d' % (strand, f_pam_start),
                'strand': strand,
                'pam': work[start:end],
                'pam_start': f_pam_start,
                'pam_end': f_pam_end,
                'protospacer': 'G' + work[proto_start:start],
                'protospacer_start': f_proto_start,
                'protospacer_end': f_proto_end,
                'nick': f_nick,
            })

    out.sort(key=lambda x: (x['pam_start'], x['strand']))
    return out


#--- pegRNA design ---------

def _build_mutation(seq, edit_start, edit_end, alt_seq):
    """
    Builds a PEGG mutation object for replacing seq[edit_start:edit_end] with
    alt_seq, reusing PEGG's own PrimeDesign parser so that variant typing
    (SNP/ONP/INS/DEL/INDEL) is decided by PEGG and cannot drift from it.
    """
    if not (0 <= edit_start <= edit_end <= len(seq)):
        raise DesignError('Edit position lies outside the sequence.')

    ref_seq = seq[edit_start:edit_end]
    if ref_seq == alt_seq:
        raise DesignError('The edited sequence is identical to the reference.')

    pd_string = '%s(%s/%s)%s' % (seq[:edit_start], ref_seq, alt_seq, seq[edit_end:])
    var_type, wt, alt, left, right, ref, alt_s = prime.primedesign_formatter(pd_string)

    mut = prime.mutation(wt, alt, left, right, var_type, ref, alt_s)
    return mut


def design(seq, edit_start, edit_end, alt_seq, pam='NGG', pam_id=None,
           rtt_length=None, pbs_length=13, proto_size=DEFAULT_PROTO_SIZE,
           orf_start=None, transcript_strand='+', silent_bystander=False,
           bystander_window_nt=5, max_bystander_muts=2):
    """
    Designs the pegRNA(s) for one edit at one PAM, and returns everything the
    browser needs to draw them: the PBS/RTT/protospacer sequences, their offsets
    on the forward strand, the edited sequence, and the silent bystander options.

    The geometry mirrors prime.pegRNA_generator() exactly. It is re-derived here
    rather than called through run(), because run() takes a dataframe of variants
    and returns scored rows, while this tool needs the intermediate coordinates
    (nick site, left/right homology arm boundaries) that run() does not surface.
    Every rule it applies -- RTT start at PAM_start - 3, PBS taken 5' of the nick
    and reverse complemented, RTT reverse complemented, PAM disruption tested by
    re-running PAM_finder over RTT[3:3+len(PAM)] -- is taken from that function.
    """
    seq = clean_sequence(seq)
    pam = clean_pam(pam)
    alt_seq = '' if not alt_seq else clean_sequence(alt_seq)
    L = len(seq)

    mut = _build_mutation(seq, edit_start, edit_end, alt_seq)

    #Silent bystanders need a declared reading frame. This mirrors PEGG's 'orf'
    #frame mode: the frame is declared against the forward strand, and the
    #sequence is assumed to be entirely coding.
    frame_mode = 'orf'
    if silent_bystander:
        if orf_start is None:
            raise DesignError(
                'Silent bystanders need a reading frame. Set the ORF offset first.')

        #PEGG's 'orf' frame mode declares the reading frame as an offset into
        #the INPUT sequence, and has no way to express a transcript running the
        #other way: ORF_start alone cannot say where the reverse frame begins.
        #Accepting transcript_strand='-' here would silently translate in a
        #frame the user never specified and return "silent" changes that are not
        #silent. Reverse-strand transcripts need genomic coordinates (PEGG's
        #'cds' mode), which this tool does not take.
        if transcript_strand != '+':
            raise DesignError(
                'Silent bystanders are only supported for a transcript on the '
                'forward strand. The reading frame is declared as an offset into '
                'the sequence you pasted, so paste the coding strand and set the '
                'ORF offset. (pegRNAs on the reverse strand are still designed '
                'normally -- it is the transcript, not the PAM, that must be '
                'forward.)')

        mut.transcript_strand = transcript_strand

    candidates = find_pams(seq, pam, proto_size)
    if pam_id is not None:
        candidates = [i for i in candidates if i['id'] == pam_id]
        if not candidates:
            raise DesignError('No PAM matching %s.' % pam_id)

    designs = []
    for cand in candidates:
        d = _design_one(seq, mut, cand, pam, rtt_length, pbs_length, proto_size,
                        orf_start, transcript_strand, silent_bystander,
                        bystander_window_nt, max_bystander_muts, frame_mode)
        if d is not None:
            designs.append(d)

    return {
        'sequence': seq,
        'variant_type': mut.variant_type,
        'ref_seq': mut.ref_seq,
        'alt_seq': mut.alt_seq,
        'edit_start': edit_start,
        'edit_end': edit_end,
        'edited_sequence': seq[:edit_start] + alt_seq + seq[edit_end:],
        'designs': designs,
    }


def _design_one(seq, mut, cand, pam, rtt_length, pbs_length, proto_size,
                orf_start, transcript_strand, silent_bystander,
                window_nt, max_muts, frame_mode):
    """Designs the pegRNA for a single PAM. Returns None if it cannot be built."""
    L = len(seq)
    orientation = cand['strand']

    #switch into the PAM strand's frame of reference, exactly as
    #prime.pegRNA_generator() does
    if orientation == '+':
        work = mut.wt_forward
        left_len = len(mut.left_seq)
        ref_len = len(mut.ref_seq)
        alt_len = len(mut.alt_seq)
        alt = mut.alt_seq
    else:
        work = mut.wt_rc
        left_len = len(mut.left_seq_rc)
        ref_len = len(mut.ref_seq_rc)
        alt_len = len(mut.alt_seq_rc)
        alt = mut.alt_seq_rc

    def to_forward(offset, length=1):
        """Offset in `work` -> offset in the forward strand."""
        if orientation == '+':
            return offset
        return L - offset - length

    #PAM/protospacer offsets in `work` coordinates
    if orientation == '+':
        pam_start = cand['pam_start']
    else:
        pam_start = L - cand['pam_end']

    rtt_start = pam_start - 3
    if rtt_start < 0:
        return None

    left_rtt = work[rtt_start:left_len]
    #the nick must sit 5' of the edit on the PAM strand, or there is no left arm
    if len(left_rtt) < 0 or rtt_start > left_len:
        return None

    #choose an RTT length if the caller did not: PEGG's own default sweep is
    #10-20 nt of right homology arm, so take the shortest that clears the edit
    if rtt_length is None:
        rtt_length = len(left_rtt) + alt_len + 10

    remaining = rtt_length - (len(left_rtt) + alt_len)
    if remaining < 0 or (left_len + ref_len + remaining) > len(work):
        return None
    if rtt_start - pbs_length < 0:
        return None

    right_rtt = work[left_len + ref_len:left_len + ref_len + remaining]
    rtt_fwd = left_rtt + alt + right_rtt

    #PAM disruption, tested the way PEGG tests it
    pam_new = rtt_fwd[3:3 + len(pam)]
    pam_disrupted = len(prime.PAM_finder(pam_new, pam)) == 0
    proto_disrupted = len(left_rtt) < 3

    pbs_work = work[rtt_start - pbs_length:rtt_start]

    #--- silent bystanders ------------------------------------------------
    #Map an offset in the RTT onto the reference base it replaces. The RTT is
    #NOT a 1:1 overlay of the reference: it carries the alt allele, so for an
    #indel every base 3' of the edit is displaced by (alt_len - ref_len).
    #Treating the RTT as a straight overlay mis-places bystander highlights by
    #exactly that amount for insertions and deletions.
    def rtt_offset_to_ref(p):
        """RTT offset -> offset in `work`, or None if p falls inside the edit."""
        edit_lo = len(left_rtt)
        edit_hi = edit_lo + alt_len
        if p < edit_lo:
            return rtt_start + p                      #left homology arm: 1:1
        if p >= edit_hi:
            return rtt_start + p - alt_len + ref_len  #right arm: shifted
        return None                                   #inside the edit itself

    bystanders = []
    bystander_note = None
    if silent_bystander:
        #Phase of rtt_fwd[0] read in the transcript's direction. In 'orf' mode
        #the frame is declared against the forward strand, so a '-' orientation
        #offset is mapped back before asking PEGG for the phase.
        if orientation == '-':
            frame_ref = len(work) - 1 - rtt_start
        else:
            frame_ref = rtt_start

        frame0 = bystander.frame_of_offset(frame_ref, frame_mode, None, None, orf_start)

        #Correct the anchor when the transcript and the PAM are on opposite
        #strands AND the edit changes length.
        #
        #silent_bystanders() reads codons on reverse_complement(RTT_fwd) and
        #re-anchors the frame with reverse_frame_anchor(frame0, len(RTT_fwd)).
        #That assumes the RTT is the same length as the reference span it
        #replaces, which holds for a substitution but not for an indel: the RTT
        #carries the alt allele, so it is (alt_len - ref_len) longer than the
        #span. The first base of the reverse complement is therefore displaced
        #by that amount, and every codon boundary shifts with it -- yielding
        #changes that are synonymous in the shifted frame and non-synonymous in
        #the real one.
        #
        #Pre-compensating the anchor by the length change makes PEGG's own
        #re-anchoring land on the phase the reference span actually has. For a
        #substitution the term is zero and this is a no-op.
        if frame0 is not None and transcript_strand != orientation:
            frame0 = (frame0 + (alt_len - ref_len)) % 3

        if frame0 is None:
            bystander_note = 'No reading frame could be resolved at this position.'
        else:
            options = bystander.silent_bystanders(
                rtt_fwd, len(left_rtt), ref_len, alt_len,
                transcript_strand, orientation, frame0,
                RTT_genomic_positions=None, frame_map=None, boundaries=None,
                window_nt=window_nt, max_muts=max_muts,
                max_candidates=None)

            for opt in options:
                new_rtt_fwd = opt['RTT']
                new_pam = new_rtt_fwd[3:3 + len(pam)]

                #Several options routinely share a position and differ only in
                #which base they put there (e.g. the third position of a 4-fold
                #degenerate codon). Spell the change out so they can be told
                #apart -- a bare position list makes them look like duplicates.
                #Label each change by its position in the SEQUENCE (1-based),
                #not by its offset inside the RTT. Two numbering systems for the
                #same event -- "C39A" in the table but position 60 on the map --
                #read as two different sites.
                _fwd = {}
                for _p in opt['positions']:
                    _w = rtt_offset_to_ref(_p)
                    _fwd[_p] = None if _w is None else to_forward(_w)

                #Sorted by position in the sequence, so a multi-change label
                #reads in the same order as the coordinates beside it. On a '-'
                #strand design RTT order runs backwards along the sequence, which
                #otherwise listed "G66A, G60A" against positions "60, 66".
                changes = []
                for _p in sorted(opt['positions'],
                                 key=lambda q: (_fwd[q] is None, _fwd[q], q)):
                    _f = _fwd[_p]
                    _where = ('%d' % (_f + 1)) if _f is not None else ('RTT+%d' % _p)
                    changes.append('%s%s%s' % (rtt_fwd[_p], _where, new_rtt_fwd[_p]))

                bystanders.append({
                    'changes': changes,
                    'label': ', '.join(changes),
                    'rtt': revcomp(new_rtt_fwd),
                    'rtt_sense': new_rtt_fwd,
                    'n_muts': opt['n_muts'],
                    'dist_to_edit': opt['dist_to_edit'],
                    #offsets into the RTT, and onto the forward strand
                    'positions': opt['positions'],
                    'positions_forward': sorted(
                        to_forward(w) for w in
                        (rtt_offset_to_ref(p) for p in opt['positions'])
                        if w is not None),
                    'pam_disrupted': len(prime.PAM_finder(new_pam, pam)) == 0,
                })

            if not bystanders:
                bystander_note = (
                    'No synonymous change is available within %d nt of the edit.'
                    % window_nt)

    #Reference bases the RTT overwrites; needed both to complete the bystanders
    #below and to place the RTT on the forward strand further down.
    ref_span_ = len(left_rtt) + ref_len + len(right_rtt)

    #--- complete each bystander -------------------------------------------
    #The PBS is only known once the PBS loop above has run, so the parts of a
    #bystander that depend on it are filled in here rather than at construction.
    #Without this a bystander could be highlighted but not actually used: it
    #carried an RTT and nothing else, so there was no oligo to order and nothing
    #to draw annealed onto the target.
    pbs_pegRNA = revcomp(pbs_work)

    #the product this design installs, and its protein -- the same treatment the
    #bystanders get, so a design can be read without re-deriving it by hand
    _f0 = to_forward(rtt_start, ref_span_)
    design_edited = (seq[:_f0]
                     + (rtt_fwd if orientation == '+' else revcomp(rtt_fwd))
                     + seq[_f0 + ref_span_:])
    design_protein = (reading_frame(design_edited, orf_start, '+')['protein']
                      if orf_start is not None else None)

    for b in bystanders:
        b['pbs'] = pbs_pegRNA
        b['extension'] = b['rtt'] + pbs_pegRNA
        b['full_pegRNA'] = cand['protospacer'] + b['rtt'] + pbs_pegRNA

        #the product this bystander installs, and what it codes for, so that the
        #"silent" claim can be seen rather than taken on trust
        b_edited = (seq[:to_forward(rtt_start, ref_span_)]
                    + (b['rtt_sense'] if orientation == '+'
                       else revcomp(b['rtt_sense']))
                    + seq[to_forward(rtt_start, ref_span_) + ref_span_:])
        b['edited_sequence'] = b_edited
        if orf_start is not None:
            b['protein'] = reading_frame(b_edited, orf_start, '+')['protein']

    #--- map every feature back onto the forward strand --------------------
    #The reference span the RTT overwrites is NOT len(RTT): the RTT carries the
    #alt allele, so for an indel it is longer (insertion) or shorter (deletion)
    #than the reference it replaces. Report the span explicitly so that anything
    #re-installing the RTT does not have to infer it and get indels wrong.
    ref_span = ref_span_
    rtt_f_start = to_forward(rtt_start, ref_span)

    #The edit's own columns on the forward strand. In `work` coordinates it
    #occupies [left_len, left_len + ref_len); a pure insertion is the zero-width
    #point at left_len. to_forward() needs the span's length to reverse it, so a
    #zero-width edit is mapped via the following base and then collapsed.
    if ref_len:
        edit_f_start = to_forward(left_len, ref_len)
        edit_f_end = edit_f_start + ref_len
    else:
        edit_f_start = (left_len if orientation == '+'
                        else to_forward(left_len, 1) + 1)
        edit_f_end = edit_f_start
    pbs_f_start = to_forward(rtt_start - pbs_length, pbs_length)
    left_arm_f = to_forward(rtt_start, len(left_rtt)) if len(left_rtt) else rtt_f_start
    right_arm_start_work = left_len + ref_len
    right_arm_f = to_forward(right_arm_start_work, len(right_rtt))

    return {
        'id': cand['id'],
        'strand': orientation,
        'pam': cand['pam'],
        'pam_start': cand['pam_start'],
        'pam_end': cand['pam_end'],
        'protospacer': cand['protospacer'],
        'protospacer_start': cand['protospacer_start'],
        'protospacer_end': cand['protospacer_end'],
        'nick': cand['nick'],

        'pbs': revcomp(pbs_work),
        'pbs_length': pbs_length,
        'pbs_start': pbs_f_start,
        'pbs_end': pbs_f_start + pbs_length,

        'rtt': revcomp(rtt_fwd),
        'rtt_sense': rtt_fwd,
        'rtt_length': len(rtt_fwd),
        'rtt_start': rtt_f_start,
        'rtt_end': rtt_f_start + ref_span,
        #reference bases the RTT replaces; differs from rtt_length for indels
        'ref_span': ref_span,

        #Where the edit sits on the forward strand, and how many bases the RTT
        #carries there. The viewer needs both: across the edit the RTT is not
        #one base per reference column, so it has to be laid out in three runs
        #(left arm, edit, right arm) rather than as a single 1:1 overlay.
        'edit_start': edit_f_start,
        'edit_end': edit_f_end,
        'alt_len': alt_len,
        'ref_len': ref_len,

        'extension': revcomp(rtt_fwd) + revcomp(pbs_work),
        'full_pegRNA': cand['protospacer'] + revcomp(rtt_fwd) + revcomp(pbs_work),

        'left_arm': left_rtt,
        'left_arm_length': len(left_rtt),
        'left_arm_start': left_arm_f,
        'right_arm': right_rtt,
        'right_arm_length': len(right_rtt),
        'right_arm_start': right_arm_f,

        'distance_to_nick': len(left_rtt),
        'RHA': len(right_rtt),
        'PAM_disrupted': bool(pam_disrupted),
        'proto_disrupted': bool(proto_disrupted),

        'edited_sequence': design_edited,
        'protein': design_protein,

        'bystanders': bystanders,
        'bystander_note': bystander_note,
    }


#--- comparison ---------

def _align(a, b):
    """
    Global alignment of two strings, returned as the two gapped strings.

    Uses PEGG's own aligner so the gap penalties match the ones PEGG applies
    when it classifies a variant, rather than introducing a second opinion.
    Falls back to no alignment when either side is empty.
    """
    if not a or not b:
        return a or ('-' * len(b)), b or ('-' * len(a))

    al = prime.aligner.align(a, b)[0]
    return str(al[0]), str(al[1])


def compare(seq, edited, orf_start=None, strand='+'):
    """
    Compares reference and edited sequences at the nucleotide level, and -- when
    a reading frame is given -- at the protein level too.

    Differences are found by ALIGNING the two sequences, not by comparing them
    position by position. Comparing by offset only works when the two are the
    same length: after an insertion or deletion everything downstream is shifted
    by the length change, so a 3 nt indel was reported as ~45 differences -- one
    for nearly every base after the edit -- rather than as the one event it is.

    Returns nt_diffs as offsets into the REFERENCE (columns the viewer draws),
    plus a per-position alignment so an indel can be shown as inserted or
    deleted rather than as a wall of mismatches.
    """
    seq = clean_sequence(seq)
    edited = clean_sequence(edited)

    ref_aln, alt_aln = _align(seq, edited)

    #Walk the alignment, recording each difference against reference offsets.
    diffs = []          #substituted reference positions
    inserted = []       #{'after': ref offset, 'bases': str}
    deleted = []        #reference offsets absent from the edited sequence
    ri = -1             #last reference offset consumed
    for r, a in zip(ref_aln, alt_aln):
        if r == '-':
            if inserted and inserted[-1]['after'] == ri:
                inserted[-1]['bases'] += a
            else:
                inserted.append({'after': ri, 'bases': a})
            continue
        ri += 1
        if a == '-':
            deleted.append(ri)
        elif r != a:
            diffs.append(ri)

    out = {
        'ref': seq,
        'alt': edited,
        'nt_diffs': diffs,
        'nt_inserted': inserted,
        'nt_deleted': deleted,
        'ref_aligned': ref_aln,
        'alt_aligned': alt_aln,
        'length_change': len(edited) - len(seq),
    }

    if orf_start is not None:
        ref_frame = reading_frame(seq, orf_start, strand)
        alt_frame = reading_frame(edited, orf_start, strand)
        ref_p, alt_p = ref_frame['protein'], alt_frame['protein']

        #Aligned for the same reason as the nucleotides: an in-frame indel
        #shifts every residue after it, and comparing by index would report the
        #whole tail as changed when only one codon was gained or lost.
        rp_aln, ap_aln = _align(ref_p, alt_p)
        aa_diffs = []
        aa_inserted = []
        pi = -1
        for r, a in zip(rp_aln, ap_aln):
            if r == '-':
                #a residue the edit adds; it has no reference index of its own,
                #so record it as inserted after the last one consumed rather
                #than dropping it, which reported an in-frame insertion as no
                #protein change at all
                if aa_inserted and aa_inserted[-1]['after'] == pi:
                    aa_inserted[-1]['residues'] += a
                else:
                    aa_inserted.append({'after': pi, 'residues': a})
                continue
            pi += 1
            if a == '-':
                aa_diffs.append({'index': pi, 'ref': r, 'alt': '-'})
            elif r != a:
                aa_diffs.append({'index': pi, 'ref': r, 'alt': a})

        out['ref_protein'] = ref_p
        out['alt_protein'] = alt_p
        #the gapped forms, so the viewer can show an indel as a gap instead of
        #re-deriving a positional diff that is wrong past the edit
        out['ref_protein_aligned'] = rp_aln
        out['alt_protein_aligned'] = ap_aln
        out['aa_diffs'] = aa_diffs
        out['aa_inserted'] = aa_inserted
        out['silent'] = (len(aa_diffs) == 0 and len(ref_p) == len(alt_p))
        out['frameshift'] = (out['length_change'] % 3 != 0)

    return out
