"""
Flask app for the pegRNA design visualizer.

Thin HTTP layer over pegg_bridge: parse JSON, call the bridge, return JSON.
All biology lives in pegg_bridge / PEGG itself.
"""

import os

from flask import Flask, jsonify, render_template, request

from . import pegg_bridge
from .pegg_bridge import DesignError

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, '..'))

app = Flask(
    __name__,
    template_folder=os.path.join(_ROOT, 'templates'),
    static_folder=os.path.join(_ROOT, 'static'),
)


@app.errorhandler(DesignError)
def _design_error(err):
    """User-fixable input problems come back as 400 with a readable message."""
    return jsonify({'error': str(err)}), 400


def _body():
    return request.get_json(silent=True) or {}


def _int_or_none(val):
    if val is None or val == '':
        return None
    return int(val)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/translate', methods=['POST'])
def api_translate():
    """Reading frame for the pasted sequence, on either strand."""
    b = _body()
    seq = pegg_bridge.clean_sequence(b.get('sequence'))
    orf_start = int(b.get('orf_start', 0))
    strand = b.get('strand', '+')

    frame = pegg_bridge.reading_frame(seq, orf_start, strand)
    frame['sequence'] = seq
    frame['revcomp'] = pegg_bridge.revcomp(seq)
    return jsonify(frame)


@app.route('/api/pams', methods=['POST'])
def api_pams():
    """Every PAM on both strands, with its protospacer and nick site."""
    b = _body()
    pams = pegg_bridge.find_pams(
        b.get('sequence'),
        b.get('pam', 'NGG'),
        int(b.get('proto_size', pegg_bridge.DEFAULT_PROTO_SIZE)),
    )
    return jsonify({'pams': pams})


@app.route('/api/design', methods=['POST'])
def api_design():
    """pegRNA design(s) for one edit, optionally with silent bystanders."""
    b = _body()
    result = pegg_bridge.design(
        seq=b.get('sequence'),
        edit_start=int(b.get('edit_start', 0)),
        edit_end=int(b.get('edit_end', 0)),
        alt_seq=b.get('alt_seq', ''),
        pam=b.get('pam', 'NGG'),
        pam_id=b.get('pam_id'),
        rtt_length=_int_or_none(b.get('rtt_length')),
        pbs_length=int(b.get('pbs_length', 13)),
        proto_size=int(b.get('proto_size', pegg_bridge.DEFAULT_PROTO_SIZE)),
        orf_start=_int_or_none(b.get('orf_start')),
        transcript_strand=b.get('transcript_strand', '+'),
        silent_bystander=bool(b.get('silent_bystander', False)),
        bystander_window_nt=int(b.get('bystander_window_nt', 5)),
        max_bystander_muts=int(b.get('max_bystander_muts', 2)),
    )
    return jsonify(result)


@app.route('/api/compare', methods=['POST'])
def api_compare():
    """Reference vs edited, at nucleotide and (given a frame) protein level."""
    b = _body()
    result = pegg_bridge.compare(
        seq=b.get('sequence'),
        edited=b.get('edited'),
        orf_start=_int_or_none(b.get('orf_start')),
        strand=b.get('strand', '+'),
    )
    return jsonify(result)


@app.route('/api/apply_rtt', methods=['POST'])
def api_apply_rtt():
    """
    Applies a user-supplied RTT back onto the reference and reports what it
    would install, so a hand-written RTT can be checked against the design.

    The RTT is given in pegRNA orientation (as it is synthesised, 3'->5' along
    the target), so it is reverse complemented onto the PAM strand first -- the
    same orientation prime.pegRNA_generator() builds it in before its final
    reverse complement.
    """
    b = _body()
    seq = pegg_bridge.clean_sequence(b.get('sequence'))
    rtt = pegg_bridge.clean_sequence(b.get('rtt'))
    rtt_start = int(b.get('rtt_start'))
    strand = b.get('strand', '+')
    orf_start = _int_or_none(b.get('orf_start'))
    replaced_len = _int_or_none(b.get('replaced_length'))

    #RTT as it lies on the forward strand
    rtt_sense = pegg_bridge.revcomp(rtt) if b.get('rtt_is_pegRNA', True) else rtt
    if strand == '-':
        rtt_forward = pegg_bridge.revcomp(rtt_sense)
    else:
        rtt_forward = rtt_sense

    span = replaced_len if replaced_len is not None else len(rtt_forward)
    if rtt_start < 0 or rtt_start + span > len(seq):
        raise DesignError('RTT does not fit onto the sequence at that position.')

    edited = seq[:rtt_start] + rtt_forward + seq[rtt_start + span:]

    result = pegg_bridge.compare(seq, edited, orf_start, '+')
    result['rtt_forward'] = rtt_forward
    result['rtt_start'] = rtt_start
    result['edited_sequence'] = edited
    return jsonify(result)


@app.route('/api/health')
def api_health():
    return jsonify({'ok': True, 'pegg': pegg_bridge.PEGG_SOURCE})


def main():
    port = int(os.environ.get('PORT', 5000))
    app.run(host='127.0.0.1', port=port, debug=bool(os.environ.get('PEGVIZ_DEBUG')))


if __name__ == '__main__':
    main()
