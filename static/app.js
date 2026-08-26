/* pegRNA Design Visualizer — frontend.

   State lives in one object; every render is a pure function of it. All
   biology (PAM search, pegRNA geometry, silent bystanders, translation) is done
   server-side by PEGG — this file only collects inputs, draws results, and maps
   mouse selection onto sequence offsets.

   All offsets are 0-based, half-open, on the FORWARD strand, matching what the
   backend returns. The UI shows 1-based positions.                            */

const EXAMPLE =
  'ATGGCTAGCACCGGTAAGCTTCCCGGGTACCTGCAGGTCGACTCTAGAGGATCCCCGGG' +
  'CGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTA';

const LINE_WIDTH = 60;   // bases per row in the viewer

// Which feature colours a base when several overlap; most specific first.
// The PBS lies wholly inside the protospacer (it is complementary to the nicked
// strand just 5' of the nick), and the RTT covers the PAM, so the containing
// feature must lose to the contained one or it hides it completely.
const FEATURE_PRIORITY = ['f-edit', 'f-by', 'f-pam', 'f-pbs', 'f-rtt', 'f-spacer'];

const state = {
  sequence: '',
  revcomp: '',
  frame: null,          // {protein, codons, ...}
  pams: [],
  selectedPam: null,
  editStart: 0,
  editEnd: 0,
  altSeq: '',
  designs: [],
  activeDesign: null,
  activeBystander: null,
  dragAnchor: null,
};

const $ = (id) => document.getElementById(id);
const el = (tag, cls, txt) => {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (txt !== undefined) n.textContent = txt;
  return n;
};

function toast(msg, ms = 2600) {
  const t = $('toast');
  t.textContent = msg;
  t.classList.remove('hidden');
  clearTimeout(toast._t);
  toast._t = setTimeout(() => t.classList.add('hidden'), ms);
}

async function api(path, body) {
  const res = await fetch(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
  const data = await res.json().catch(() => ({ error: 'Malformed response from server.' }));
  if (!res.ok) throw new Error(data.error || ('Request failed (' + res.status + ')'));
  return data;
}

/* ---------- 1. sequence & reading frame ---------- */

function cleanLocal(s) {
  // Mirrors the server's cleaner so the viewer can render before any request
  // returns. The server remains the authority; this only avoids a flash of
  // unrendered text.
  return s.split('\n').filter(l => !l.startsWith('>')).join('')
          .replace(/[\s\d]/g, '').toUpperCase().replace(/U/g, 'T');
}

async function refreshFrame() {
  const raw = cleanLocal($('seq').value);
  if (!raw) {
    state.sequence = '';
    state.frame = null;
    render();
    return;
  }
  try {
    const data = await api('/api/translate', {
      sequence: raw,
      orf_start: +$('orfStart').value,
      strand: $('frameStrand').value,
    });
    state.sequence = data.sequence;
    state.revcomp = data.revcomp;
    state.frame = data;
    clampEdit();
    renderProtein();
    render();
  } catch (e) {
    state.sequence = '';
    state.frame = null;
    showError('proteinOut', e.message);
    render();
  }
}

function renderProtein() {
  const box = $('proteinOut');
  if (!state.frame || !$('showProtein').checked) {
    box.classList.add('hidden');
    return;
  }
  box.classList.remove('hidden');
  box.innerHTML = '';
  const strand = state.frame.strand === '+' ? 'forward (+)' : 'reverse (−)';
  box.appendChild(el('div', null,
    `Protein — ${strand} strand, frame ${state.frame.orf_start} · ${state.frame.protein.length} aa`));

  const line = el('div');
  for (const c of state.frame.protein) {
    const s = el('span', c === '*' ? 'stop' : null, c);
    line.appendChild(s);
  }
  box.appendChild(line);
}

function showError(containerId, msg) {
  const box = $(containerId);
  box.classList.remove('hidden');
  box.innerHTML = '';
  box.appendChild(el('div', 'error', msg));
}

/* ---------- 2. PAMs ---------- */

async function findPams() {
  if (!state.sequence) return toast('Paste a sequence first.');
  try {
    const data = await api('/api/pams', {
      sequence: state.sequence,
      pam: $('pam').value.trim().toUpperCase(),
      proto_size: +$('protoSize').value,
    });
    state.pams = data.pams;
    state.selectedPam = null;
    renderPams();
    render();
    toast(`${data.pams.length} PAM site${data.pams.length === 1 ? '' : 's'} found.`);
  } catch (e) {
    showError('pamList', e.message);
  }
}

function visiblePams() {
  const f = $('pamStrandFilter').value;
  return f === 'both' ? state.pams : state.pams.filter(p => p.strand === f);
}

function renderPams() {
  const box = $('pamList');
  box.innerHTML = '';
  const list = visiblePams();
  if (!list.length) {
    box.appendChild(el('p', 'note', 'No PAM sites — try a different PAM or sequence.'));
    return;
  }
  for (const p of list) {
    const chip = el('div', 'pam-chip' + (state.selectedPam === p.id ? ' active' : ''));
    chip.appendChild(el('div', 'seq', `${p.pam} @ ${p.pam_start + 1}`));
    chip.appendChild(el('div', 'meta',
      `${p.strand === '+' ? 'fwd' : 'rev'} · nick ${p.nick} · ${p.protospacer}`));
    chip.onclick = () => {
      state.selectedPam = (state.selectedPam === p.id) ? null : p.id;
      renderPams();
      render();
    };
    box.appendChild(chip);
  }
}

/* ---------- 3. edit selection ---------- */

function clampEdit() {
  const L = state.sequence.length;
  state.editStart = Math.max(0, Math.min(state.editStart, L));
  state.editEnd = Math.max(state.editStart, Math.min(state.editEnd, L));
  syncEditInputs();
}

function syncEditInputs() {
  // inputs are 1-based inclusive; state is 0-based half-open
  $('editStart').value = state.editStart + 1;
  $('editEnd').value = state.editEnd;
  $('refSeq').value = state.sequence.slice(state.editStart, state.editEnd) || '(none — insertion)';
  updateVarType();
}

function updateVarType() {
  const ref = state.sequence.slice(state.editStart, state.editEnd);
  const alt = cleanLocal($('altSeq').value);
  state.altSeq = alt;
  const b = $('varType');
  let t = '';
  if (!ref && !alt) t = '';
  else if (!ref) t = 'INS';
  else if (!alt) t = 'DEL';
  else if (ref.length === alt.length) t = ref.length === 1 ? 'SNP' : 'ONP';
  else t = 'INDEL';
  b.textContent = t;
  b.className = 'badge' + (t && t !== 'SNP' && t !== 'ONP' ? ' warn' : '');
}

function readEditInputs() {
  const s = Math.max(1, +$('editStart').value || 1);
  const e = Math.max(0, +$('editEnd').value || 0);
  state.editStart = s - 1;
  state.editEnd = Math.max(state.editStart, e);
  clampEdit();
  render();
}

/* ---------- sequence viewer ---------- */

/* Builds a per-base class map from the current selection and active design, so
   each base cell can be shaded with whatever features overlap it. */
function featureMap() {
  const L = state.sequence.length;
  const m = Array.from({ length: L }, () => ({ fwd: [], rev: [] }));
  const add = (from, to, cls, strand) => {
    for (let i = Math.max(0, from); i < Math.min(L, to); i++) {
      if (strand === '-') m[i].rev.push(cls); else m[i].fwd.push(cls);
    }
  };

  // PAM/spacer for every visible PAM, lightly; the selected one is emphasised
  const sel = state.pams.find(p => p.id === state.selectedPam);
  const d = state.activeDesign;

  if (d) {
    add(d.protospacer_start, d.protospacer_end, 'f-spacer', d.strand);
    add(d.pam_start, d.pam_end, 'f-pam', d.strand);
    add(d.pbs_start, d.pbs_end, 'f-pbs', d.strand);
    add(d.rtt_start, d.rtt_end, 'f-rtt', d.strand);
  } else if (sel) {
    add(sel.protospacer_start, sel.protospacer_end, 'f-spacer', sel.strand);
    add(sel.pam_start, sel.pam_end, 'f-pam', sel.strand);
  }

  // the intended edit, on both strands
  add(state.editStart, Math.max(state.editEnd, state.editStart + 1), 'f-edit', '+');
  add(state.editStart, Math.max(state.editEnd, state.editStart + 1), 'f-edit', '-');

  // silent bystander positions, if one is being previewed
  if (state.activeBystander) {
    for (const p of state.activeBystander.positions_forward) {
      if (p >= 0 && p < L) {
        m[p].fwd.push('f-by');
        m[p].rev.push('f-by');
      }
    }
  }

  return m;
}

function nickPositions() {
  const d = state.activeDesign || state.pams.find(p => p.id === state.selectedPam);
  return d ? [{ pos: d.nick, strand: d.strand }] : [];
}

function render() {
  const v = $('viewer');
  v.innerHTML = '';
  const seq = state.sequence;
  if (!seq) {
    v.appendChild(el('p', 'placeholder', 'Paste a sequence above to begin.'));
    renderLegend();
    return;
  }

  const fm = featureMap();
  const nicks = nickPositions();
  const showAA = $('showProtein').checked && state.frame;

  // codon lookup by forward-strand start offset, for the protein track
  const codonAt = new Map();
  if (showAA) for (const c of state.frame.codons) codonAt.set(c.start, c);

  for (let start = 0; start < seq.length; start += LINE_WIDTH) {
    const end = Math.min(start + LINE_WIDTH, seq.length);
    const block = el('div', 'block');

    // ruler: a tick every 10 bases
    const ruler = el('div', 'track ruler');
    for (let i = start; i < end; i++) {
      const c = el('span', 'cell');
      if ((i + 1) % 10 === 0) c.textContent = String(i + 1).slice(-3);
      ruler.appendChild(c);
    }
    block.appendChild(ruler);

    // protein track (above forward strand) — only for '+' strand frames
    if (showAA && state.frame.strand === '+') {
      const aaTrack = el('div', 'track');
      let i = start;
      while (i < end) {
        const c = codonAt.get(i);
        if (c && i + 3 <= end) {
          const cls = 'aa' + (c.aa === '*' ? ' stop' : '') + (isCodonEdited(c) ? ' changed' : '');
          aaTrack.appendChild(el('span', cls, c.aa));
          i += 3;
        } else {
          aaTrack.appendChild(el('span', 'cell'));
          i += 1;
        }
      }
      block.appendChild(aaTrack);
    }

    // forward strand
    block.appendChild(strandTrack(seq, start, end, fm, 'fwd', nicks, '+'));
    // reverse strand (complement, displayed 5'->3' left-to-right of the top strand)
    block.appendChild(strandTrack(complement(seq), start, end, fm, 'rev', nicks, '-'));

    // protein track below, for '-' strand frames
    if (showAA && state.frame.strand === '-') {
      const aaTrack = el('div', 'track');
      let i = start;
      while (i < end) {
        const c = codonAt.get(i);
        if (c && i + 3 <= end) {
          const cls = 'aa' + (c.aa === '*' ? ' stop' : '') + (isCodonEdited(c) ? ' changed' : '');
          aaTrack.appendChild(el('span', cls, c.aa));
          i += 3;
        } else {
          aaTrack.appendChild(el('span', 'cell'));
          i += 1;
        }
      }
      block.appendChild(aaTrack);
    }

    v.appendChild(block);
  }

  renderLegend();
}

function isCodonEdited(codon) {
  return codon.start < state.editEnd && codon.end > state.editStart &&
         state.editEnd > state.editStart;
}

function complement(s) {
  const M = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' };
  let o = '';
  for (const c of s) o += (M[c] || 'N');
  return o;
}

function strandTrack(seqStr, start, end, fm, kind, nicks, strandKey) {
  const track = el('div', 'track ' + kind);
  for (let i = start; i < end; i++) {
    const classes = ['cell', 'sel-target'];
    const feats = kind === 'fwd' ? fm[i].fwd : fm[i].rev;
    // A base can carry several features at once — the RTT spans the PAM, and
    // the edit sits inside the RTT. Only one background can show, so pick by
    // how much it tells the reader rather than by insertion order: the edit and
    // the bystanders are the point of the design, the PAM and the nick are what
    // place it, and the RTT/PBS are the broad extents those sit inside.
    // Ordering by insertion instead let the RTT paint over every PAM.
    if (feats.length) classes.push(FEATURE_PRIORITY.find(c => feats.includes(c)));
    if (i >= state.editStart && i < state.editEnd) classes.push('f-sel');

    const cell = el('span', classes.join(' '), seqStr[i]);
    cell.dataset.i = i;

    for (const n of nicks) {
      if (n.pos === i && ((n.strand === '+') === (kind === 'fwd'))) {
        cell.classList.add('nick');
        cell.title = 'nick site';
      }
    }
    track.appendChild(cell);
  }
  attachSelection(track);
  return track;
}

/* drag-to-select across base cells */
function attachSelection(track) {
  track.addEventListener('mousedown', (ev) => {
    const t = ev.target.closest('.cell');
    if (!t || t.dataset.i === undefined) return;
    ev.preventDefault();
    state.dragAnchor = +t.dataset.i;
    state.editStart = state.dragAnchor;
    state.editEnd = state.dragAnchor + 1;
    syncEditInputs();
    render();
  });
  track.addEventListener('mouseover', (ev) => {
    if (state.dragAnchor === null) return;
    const t = ev.target.closest('.cell');
    if (!t || t.dataset.i === undefined) return;
    const i = +t.dataset.i;
    state.editStart = Math.min(state.dragAnchor, i);
    state.editEnd = Math.max(state.dragAnchor, i) + 1;
    syncEditInputs();
    render();
  });
}
document.addEventListener('mouseup', () => { state.dragAnchor = null; });

function renderLegend() {
  const L = $('legend');
  L.innerHTML = '';
  const items = [
    ['--pam-bg', 'PAM'],
    ['--spacer-bg', 'Spacer / protospacer'],
    ['--pbs-bg', 'PBS'],
    ['--rtt-bg', 'RTT'],
    ['--edit-bg', 'Edit'],
    ['--bystand-bg', 'Silent bystander'],
  ];
  for (const [v, label] of items) {
    const s = el('span');
    const i = el('i');
    i.style.background = `var(${v})`;
    s.appendChild(i);
    s.appendChild(document.createTextNode(label));
    L.appendChild(s);
  }
  const n = el('span', null, '│ = nick site');
  L.appendChild(n);
}

/* ---------- 4. pegRNA design ---------- */

async function runDesign() {
  if (!state.sequence) return toast('Paste a sequence first.');
  const ref = state.sequence.slice(state.editStart, state.editEnd);
  const alt = cleanLocal($('altSeq').value);
  if (ref === alt) return toast('Choose bases to edit and a replacement.');

  const body = {
    sequence: state.sequence,
    edit_start: state.editStart,
    edit_end: state.editEnd,
    alt_seq: alt,
    pam: $('pam').value.trim().toUpperCase(),
    pam_id: state.selectedPam || null,
    pbs_length: +$('pbsLen').value,
    proto_size: +$('protoSize').value,
    rtt_length: $('rttLen').value === '' ? null : +$('rttLen').value,
    silent_bystander: $('silentBy').checked,
    bystander_window_nt: +$('byWindow').value,
    max_bystander_muts: +$('byMax').value,
    transcript_strand: $('txStrand').value,
  };
  if ($('silentBy').checked) body.orf_start = +$('orfStart').value;

  const btn = $('designBtn');
  btn.disabled = true;
  btn.textContent = 'Designing…';
  try {
    const data = await api('/api/design', body);
    state.designs = data.designs;
    state.activeDesign = data.designs[0] || null;
    state.activeBystander = null;
    renderDesigns(data);
    render();
    if (!data.designs.length) toast('No pegRNA can be built here — try another PAM or a longer RTT.');
  } catch (e) {
    showError('designOut', e.message);
  } finally {
    btn.disabled = false;
    btn.textContent = 'Design pegRNA';
  }
}

function renderDesigns(data) {
  const box = $('designOut');
  box.innerHTML = '';

  const head = el('p', 'note',
    `${data.variant_type} · ${data.ref_seq || '–'} → ${data.alt_seq || '–'} at ${data.edit_start + 1} · ` +
    `${data.designs.length} design${data.designs.length === 1 ? '' : 's'}`);
  box.appendChild(head);

  for (const d of data.designs) {
    const card = el('div', 'design-card');

    const h = el('h3');
    h.appendChild(document.createTextNode(
      `${d.strand === '+' ? 'Forward' : 'Reverse'} strand · PAM ${d.pam} @ ${d.pam_start + 1}`));
    if (d.PAM_disrupted) h.appendChild(el('span', 'badge ok', 'PAM disrupted'));
    if (d.proto_disrupted) h.appendChild(el('span', 'badge warn', 'spacer disrupted'));
    if (state.activeDesign && state.activeDesign.id === d.id) h.appendChild(el('span', 'badge', 'shown'));
    card.appendChild(h);

    const dl = el('dl', 'kv');
    const put = (k, v) => { dl.appendChild(el('dt', null, k)); dl.appendChild(el('dd', null, v)); };
    put('Spacer', d.protospacer);
    put('PBS', `${d.pbs}  (${d.pbs_length} nt)`);
    put('RTT', `${d.rtt}  (${d.rtt_length} nt)`);
    put('3′ extension', d.extension);
    put('Full pegRNA', d.full_pegRNA);
    put('Nick → edit', `${d.distance_to_nick} nt`);
    put('RHA', `${d.RHA} nt`);
    card.appendChild(dl);

    const show = el('button', 'ghost', 'Show on map');
    show.onclick = () => {
      state.activeDesign = d;
      state.activeBystander = null;
      renderDesigns(data);
      render();
    };
    card.appendChild(show);

    const useRtt = el('button', 'ghost', 'Send RTT to comparison');
    useRtt.onclick = () => {
      $('rttInput').value = d.rtt;
      $('rttOrient').value = 'peg';
      state.activeDesign = d;
      toast('RTT copied into the comparison box.');
      render();
    };
    card.appendChild(useRtt);

    if (d.bystander_note) card.appendChild(el('p', 'note', d.bystander_note));

    if (d.bystanders && d.bystanders.length) {
      const det = el('details');
      det.appendChild(el('summary', null,
        `${d.bystanders.length} silent bystander option${d.bystanders.length === 1 ? '' : 's'}`));

      const scroll = el('div', 'by-scroll');
      const tb = el('table', 'bystanders');
      const thead = el('thead');
      const hr = el('tr');
      for (const t of ['Change', '#', 'Dist. to edit', 'PAM KO', 'RTT']) hr.appendChild(el('th', null, t));
      thead.appendChild(hr);
      tb.appendChild(thead);

      const body = el('tbody');
      for (const b of d.bystanders) {
        const tr = el('tr');
        if (state.activeBystander === b) tr.classList.add('active');
        tr.appendChild(el('td', 'mono', b.label));
        tr.appendChild(el('td', null, String(b.n_muts)));
        tr.appendChild(el('td', null, `${b.dist_to_edit} nt`));
        tr.appendChild(el('td', null, b.pam_disrupted ? 'yes' : '—'));
        tr.appendChild(el('td', 'mono', b.rtt));
        tr.onclick = () => {
          state.activeDesign = d;
          state.activeBystander = (state.activeBystander === b) ? null : b;
          $('rttInput').value = b.rtt;
          renderDesigns(data);
          render();
        };
        body.appendChild(tr);
      }
      tb.appendChild(body);
      scroll.appendChild(tb);
      det.appendChild(scroll);
      card.appendChild(det);
    }

    box.appendChild(card);
  }
}

/* ---------- 5. RTT comparison ---------- */

async function runCompare() {
  const rtt = cleanLocal($('rttInput').value);
  if (!rtt) return toast('Paste an RTT to compare.');
  const d = state.activeDesign;
  if (!d) return toast('Design a pegRNA first, so the RTT has a nick to sit at.');

  try {
    const data = await api('/api/apply_rtt', {
      sequence: state.sequence,
      rtt: rtt,
      rtt_start: d.rtt_start,
      strand: d.strand,
      rtt_is_pegRNA: $('rttOrient').value === 'peg',
      // the reference span the RTT overwrites, which for an indel is not the
      // RTT's own length
      replaced_length: d.ref_span,
      orf_start: +$('orfStart').value,
    });
    renderCompare(data);
  } catch (e) {
    showError('compareOut', e.message);
  }
}

function renderCompare(c) {
  const box = $('compareOut');
  box.innerHTML = '';

  const summary = el('p');
  const badges = [];
  if (c.silent) badges.push(['ok', 'silent — protein unchanged']);
  else if (c.frameshift) badges.push(['warn', 'frameshift']);
  else if (c.aa_diffs && c.aa_diffs.length) badges.push(['warn', `${c.aa_diffs.length} aa changed`]);
  badges.push([null, `${c.nt_diffs.length} nt differ`]);
  if (c.length_change) badges.push(['warn', `${c.length_change > 0 ? '+' : ''}${c.length_change} nt`]);
  for (const [k, t] of badges) summary.appendChild(el('span', 'badge' + (k ? ' ' + k : ''), t));
  box.appendChild(summary);

  if (c.ref_protein !== undefined) {
    const wrap = el('div', 'design-card');
    wrap.appendChild(el('h3', null, 'Protein'));
    wrap.appendChild(proteinDiff('Reference', c.ref_protein, c.alt_protein));
    wrap.appendChild(proteinDiff('Edited', c.alt_protein, c.ref_protein));
    if (c.aa_diffs && c.aa_diffs.length) {
      const list = c.aa_diffs.map(d => `${d.ref}${d.index + 1}${d.alt}`).join(', ');
      wrap.appendChild(el('p', 'note', 'Changes: ' + list));
    }
    box.appendChild(wrap);
  }

  const seqBox = el('div', 'design-card');
  seqBox.appendChild(el('h3', null, 'Edited sequence'));
  const line = el('div', 'seqline');
  const diffs = new Set(c.nt_diffs);
  for (let i = 0; i < c.alt.length; i++) {
    const s = el('span', diffs.has(i) ? 'f-edit' : null, c.alt[i]);
    line.appendChild(s);
  }
  seqBox.appendChild(line);
  box.appendChild(seqBox);
}

function proteinDiff(label, seq, other) {
  const d = el('div');
  d.appendChild(el('div', 'note', label));
  const line = el('div', 'seqline');
  for (let i = 0; i < seq.length; i++) {
    const changed = other[i] !== undefined && other[i] !== seq[i];
    line.appendChild(el('span', changed ? 'f-edit' : (seq[i] === '*' ? 'stop' : null), seq[i]));
  }
  d.appendChild(line);
  return d;
}

/* ---------- wiring ---------- */

let seqTimer = null;
$('seq').addEventListener('input', () => {
  clearTimeout(seqTimer);
  seqTimer = setTimeout(refreshFrame, 250);
});
$('orfStart').addEventListener('change', refreshFrame);
$('frameStrand').addEventListener('change', refreshFrame);
$('showProtein').addEventListener('change', () => { renderProtein(); render(); });
$('loadExample').addEventListener('click', () => {
  $('seq').value = EXAMPLE;
  $('altSeq').value = 'A';
  state.editStart = 60; state.editEnd = 61;
  refreshFrame().then(findPams);
});

$('findPams').addEventListener('click', findPams);
$('pamStrandFilter').addEventListener('change', renderPams);

$('editStart').addEventListener('change', readEditInputs);
$('editEnd').addEventListener('change', readEditInputs);
$('altSeq').addEventListener('input', () => { updateVarType(); render(); });

$('designBtn').addEventListener('click', runDesign);
$('compareBtn').addEventListener('click', runCompare);

renderLegend();
