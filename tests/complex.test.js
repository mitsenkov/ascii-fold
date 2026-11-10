import test from 'node:test';
import assert from 'node:assert/strict';
import { parsePDB, parseCIF } from '../src/parse.js';
import { buildGeometry } from '../src/geometry.js';

function atom(serial, chain, num, x, y = 0, z = 0, name = 'CA', alt = ' ', insertion = ' ') {
  const line = Array(80).fill(' ');
  const put = (start, text) => [...text].forEach((c, i) => { line[start + i] = c; });
  put(0, 'ATOM  '); put(6, String(serial).padStart(5)); put(12, name.padStart(4));
  put(16, alt); put(17, 'ALA'); put(21, chain); put(22, String(num).padStart(4)); put(26, insertion);
  [x, y, z].forEach((v, i) => put(30 + 8 * i, v.toFixed(3).padStart(8)));
  put(60, ' 85.00');
  return line.join('');
}
function chain(id, offset, count = 8) {
  return Array.from({ length: count }, (_, i) => atom(i + 1, id, i + 1, i * 3.8, offset, Math.sin(i)));
}
function ssRecord(type, id, start, end) {
  const line = Array(80).fill(' ');
  const put = (i, s) => [...s].forEach((c, j) => { line[i + j] = c; });
  put(0, type === 'H' ? 'HELIX ' : 'SHEET ');
  put(type === 'H' ? 19 : 21, id);
  put(type === 'H' ? 21 : 22, String(start).padStart(4));
  put(33, String(end).padStart(4));
  return line.join('');
}
const pdb = [...chain('A', 0), 'TER', ...chain('B', 18), 'ENDMDL', ...chain('C', 90)].join('\n');

test('PDB retains both chains with overlapping residue numbers, first model only', () => {
  const p = parsePDB(pdb);
  assert.deepEqual(p.chains.map(c => [c.chain, c.rows.length]), [['A', 8], ['B', 8]]);
  assert.equal(p.rows.length, 16);
  assert.equal(p.chain, 'A');
});

test('PDB secondary structure is assigned independently by chain', () => {
  const p = parsePDB([ssRecord('H', 'A', 1, 8), ssRecord('E', 'B', 1, 8), pdb].join('\n'));
  assert.ok(p.chains[0].ss.every(s => s === 'H'));
  assert.ok(p.chains[1].ss.every(s => s === 'E'));
});

test('PDB keeps insertion codes and skips non-primary alternate locations', () => {
  const p = parsePDB([atom(1, ' ', 1, 0), atom(2, ' ', 1, 4, 0, 0, 'CA', ' ', 'A'), atom(3, ' ', 2, 8, 0, 0, 'CA', 'B')].join('\n'));
  assert.equal(p.chains[0].chain, ' ');
  assert.equal(p.rows.length, 2);
});

function cif() {
  const columns = ['group_PDB', 'label_atom_id', 'label_alt_id', 'auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'B_iso_or_equiv', 'pdbx_PDB_model_num'];
  const atoms = [];
  for (const [id, offset, model] of [['AA', 0, 1], ['BB', 18, 1], ['CC', 90, 2]])
    for (let i = 0; i < 8; i++) atoms.push(`ATOM CA . ${id} ${i + 1} ? ${i * 3.8} ${offset} 0 85 ${model}`);
  return ['data_synthetic', 'loop_', ...columns.map(c => `_atom_site.${c}`), ...atoms, '#'].join('\n');
}
test('mmCIF retains multi-character chain IDs and first model only', () => {
  const p = parseCIF(cif());
  assert.deepEqual(p.chains.map(c => [c.chain, c.rows.length]), [['AA', 8], ['BB', 8]]);
  assert.equal(p.rows.length, 16);
});

test('complex geometry preserves original relative positions without interchain ribbons', () => {
  const p = parsePDB(pdb), g = buildGeometry(p, { flat: 0 });
  const center = p.rows.reduce((a, r) => a.map((v, k) => v + r.CA[k] / p.rows.length), [0, 0, 0]);
  assert.equal(g.N, 16);
  for (const q of g.pts) {
    assert.ok([...q.p, ...q.side, ...q.nrm].every(Number.isFinite));
    if (q.t === 0) q.p.forEach((v, k) => assert.ok(Math.abs(v - (p.rows[q.res].CA[k] - center[k])) < 1e-9));
    assert.ok(Math.hypot(...q.p) + q.r <= g.R);
  }
  for (let i = 0; i + 1 < g.pts.length; i++)
    if (g.pts[i].chain !== g.pts[i + 1].chain) assert.equal(g.pts[i].gap, true);
  assert.deepEqual([...new Set(g.pts.map(q => q.chain))], ['A', 'B']);
});

test('single-chain geometry matches the legacy input shape', () => {
  const p = parsePDB(chain('A', 0).join('\n'));
  const legacy = { rows: p.rows, chain: p.chain, ss: p.ss, isAF: p.isAF };
  assert.deepEqual(buildGeometry(p), buildGeometry(legacy));
});

test('short peptide chains produce finite geometry', () => {
  const p = parsePDB([...chain('A', 0), ...chain('B', 18, 2)].join('\n'));
  const g = buildGeometry(p);
  assert.equal(g.N, 10);
  assert.ok(g.pts.every(q => [...q.p, ...q.nrm].every(Number.isFinite)));
});

test('complex confidence summary weights all residues and cylinder geometry stays finite', () => {
  const p = parsePDB([ssRecord('H', 'A', 1, 8), ssRecord('H', 'B', 1, 8), pdb].join('\n')); p.isAF = true;
  const g = buildGeometry(p, { hStyle: 'cylinder' });
  assert.equal(g.hasPlddt, true);
  assert.equal(g.meanP, 85);
  assert.ok(g.pts.every(q => [...q.p, ...q.nrm].every(Number.isFinite)));
});
