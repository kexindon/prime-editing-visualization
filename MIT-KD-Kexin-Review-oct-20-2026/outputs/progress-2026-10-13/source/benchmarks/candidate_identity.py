"""Strict target/guide identity for the no-bystander preparation fixture.

This describes guide sequence geometry, not complete experimental constructs.
PRIDICT2's introduced-allele case annotation is normalized only where specified.
"""
import hashlib
import json


def reverse_complement(sequence):
    if not isinstance(sequence, str) or set(sequence) - set('ACGT'):
        raise ValueError('Expected uppercase unambiguous DNA')
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def identity_hash(key):
    return hashlib.sha256(json.dumps(
        ['target-guide-v1', *key], separators=(',', ':')
    ).encode()).hexdigest()


def pegg_key(row):
    wt, alt, strand = row['wt_w_context'], row['alt_w_context'], row['PAM_strand']
    if strand not in ('+', '-'):
        raise ValueError('Unexpected PEGG strand')
    if str(row.get('has_silent_bystander', False)).lower() not in ('false', '0', ''):
        raise ValueError('This comparison supports no-bystander candidates only')
    oriented = wt if strand == '+' else reverse_complement(wt)
    p = int(row['PAM_start']) if strand == '+' else len(wt) - int(row['PAM_start'])
    native30 = oriented[p - 24:p + 6]
    if len(native30) != 30 or native30 != row['Protospacer_30']:
        raise ValueError('PEGG native30 disagrees with context/PAM coordinates')
    if oriented[p:p + 3] != row['PAM']:
        raise ValueError('PEGG PAM disagrees with context')
    return (wt, alt, strand, p, native30, row['RTT'], row['PBS'])


def pridict_key(row):
    if row['Target-Strand'] not in ('Fw', 'Rv'):
        raise ValueError('Unexpected PRIDICT2 strand')
    strand = '+' if row['Target-Strand'] == 'Fw' else '-'
    oriented = row['Original_Sequence']
    edited = row['Edited_Sequence'].upper()
    window = row['wide_initial_target']
    h = oriented.find(window)
    if h < 0 or oriented.find(window, h + 1) != -1:
        raise ValueError('PRIDICT2 target window is absent or ambiguous')
    p = h + 29  # native19 starts at window offset10; PAM follows at29
    native30 = oriented[p - 24:p + 6]
    if len(native30) != 30 or native30 != window[5:35]:
        raise ValueError('PRIDICT2 native30 is inconsistent')
    if oriented[p + 1:p + 3] != 'GG':
        raise ValueError('Unexpected non-NGG PRIDICT2 PAM')
    if row['Spacer-Sequence'] != 'G' + oriented[p - 19:p]:
        raise ValueError('PRIDICT2 forced-G spacer is inconsistent')
    wt = oriented if strand == '+' else reverse_complement(oriented)
    alt = edited if strand == '+' else reverse_complement(edited)
    return (wt, alt, strand, p, native30, row['RTrevcomp'].upper(), row['PBSrevcomp'])
