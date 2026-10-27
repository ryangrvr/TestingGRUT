"""Independent region ledger for one regenerated main-run snapshot (review of Codex 423d37b).

Own lattice energy, own FD4 force, own RK4 time step. Codex code is used only to read the binary
snapshot (pilot.snapshot) and to rebuild the published candidate mask (analyze.catalog), which is the
definition under audit. Region derivatives are checked by evolving the whole lattice +-h and
differencing, not from flux formulas.
usage: python ledger_check.py <codex pilot code dir> <snapshot.bin> <FINAL_AUDIT.json> <out.json>
"""
import sys, json
import numpy as np

code, snap, audit_path, out = sys.argv[1:5]
sys.path.insert(0, code)
import pilot, analyze

G, LX, LP, V = 0.1, 1e-3, 0.03, 10.0
assert (G, LX, LP) == (pilot.C.G, pilot.C.LX, pilot.C.LP) and abs(V / pilot.C.NU - 1) < 1e-14
V = pilot.C.NU  # 9.999999999999998 (v/m_X from GeV inputs); the snapshot was evolved with this value
BONDS = [(1, 4 / 3), (2, -1 / 12)]


def lap4(f, d):
    out = np.zeros_like(f)
    for ax in range(f.ndim):
        out += (16 * (np.roll(f, 1, ax) + np.roll(f, -1, ax)) - (np.roll(f, 2, ax) + np.roll(f, -2, ax)) - 30 * f) / (12 * d * d)
    return out


def force(q, d):
    x1, x2, ph = q
    r2 = x1 * x1 + x2 * x2
    mX = G * G * ph * ph + LX * r2
    return np.array([lap4(x1, d) - mX * x1, lap4(x2, d) - mX * x2,
                     lap4(ph, d) - LP * ph * (ph * ph - V * V) - G * G * ph * r2])


def densities(q, p, d):
    """Per-site X and phi energy; each bond's gradient energy split half/half between its ends."""
    x1, x2, ph = q
    r2 = x1 * x1 + x2 * x2
    eX = 0.5 * (p[0] ** 2 + p[1] ** 2) + 0.5 * G * G * ph * ph * r2 + LX / 4 * r2 * r2
    ep = 0.5 * p[2] ** 2 + LP / 4 * (ph * ph - V * V) ** 2
    for ax in range(3):
        for s, w in BONDS:
            bx = w / (2 * d * d) * ((np.roll(x1, -s, ax) - x1) ** 2 + (np.roll(x2, -s, ax) - x2) ** 2)
            bp = w / (2 * d * d) * (np.roll(ph, -s, ax) - ph) ** 2
            eX += 0.5 * (bx + np.roll(bx, s, ax))
            ep += 0.5 * (bp + np.roll(bp, s, ax))
    return eX, ep, x1 * p[1] - x2 * p[0]


def rk4(q, p, h, d):
    k1q, k1p = p, force(q, d)
    k2q, k2p = p + h / 2 * k1p, force(q + h / 2 * k1q, d)
    k3q, k3p = p + h / 2 * k2p, force(q + h / 2 * k2q, d)
    k4q, k4p = p + h * k3p, force(q + h * k3q, d)
    return q + h / 6 * (k1q + 2 * k2q + 2 * k3q + k4q), p + h / 6 * (k1p + 2 * k2p + 2 * k3p + k4p)


def region(q, p, d, m):
    dv = d ** 3
    eX, ep, j = densities(q, p, d)
    psi, pi = q[0] + 1j * q[1], p[0] + 1j * p[1]
    a2 = np.abs(psi[m]) ** 2
    Q, I = j[m].sum() * dv, a2.sum() * dv
    om = Q / I
    c = np.conj(psi[m]) * pi[m]
    K = 0.5 * (np.abs(pi[m]) ** 2).sum() * dv
    Kr = 0.5 * (c.real ** 2 / a2).sum() * dv
    Ks = 0.5 * ((c.imag - om * a2) ** 2 / a2).sum() * dv
    comp = dict(kphi=0.5 * (p[2][m] ** 2).sum() * dv, vphi=(LP / 4 * (q[2][m] ** 2 - V * V) ** 2).sum() * dv)
    comp['gradphi'] = ep[m].sum() * dv - comp['kphi'] - comp['vphi']
    return dict(cells=int(m.sum()), Q=Q, E_X=eX[m].sum() * dv, E_phi=ep[m].sum() * dv, omega_fit=om,
                kinetic_X=K, rot=Q * Q / (2 * I), kinetic_excess_fraction=(Kr + Ks) / K, phi_parts=comp)


def rates(q, p, d, m, h):
    """Central difference of region totals along the exact lattice flow."""
    qa, pa = rk4(q, p, h, d)
    qb, pb = rk4(q, p, -h, d)
    ra, rb = region(qa, pa, d, m), region(qb, pb, d, m)
    return {k: (ra[k] - rb[k]) / (2 * h) for k in ['E_X', 'E_phi', 'Q']}


N, L, t, arr = pilot.snapshot(snap)
q, p = np.asarray(arr[:3], float), np.asarray(arr[3:], float)
d = L / N
pub = json.load(open(audit_path))['cases']['mixed_N96']
assert pub['t'] == t and pub['N'] == N
_, labs = analyze.catalog(q, p, L, 4.)
ax = np.arange(N) * d
res = dict(t=t, N=N, L=L, regions=[])
for c in [x for x in pub['candidates'] if x['candidate']['compact']][:3]:
    cid, cen = c['candidate']['id'], c['candidate']['centre']
    dd = [(ax - x + L / 2) % L - L / 2 for x in cen]
    rr = np.sqrt(dd[0][:, None, None] ** 2 + dd[1][None, :, None] ** 2 + dd[2][None, None, :] ** 2)
    sph = [s for s in c['spheres'] if s['radius'] == 6.][0]
    for name, m, ref in [('mask', labs == cid, c['mask']), ('sphere6', rr < 6., sph)]:
        mine = region(q, p, d, m)
        mine['rates'] = {f'h={h}': rates(q, p, d, m, h) for h in [2e-3, 1e-3]}
        cmp = {k: (mine[k], ref[k], mine[k] / ref[k] - 1) for k in ['cells', 'Q', 'E_X', 'E_phi', 'omega_fit', 'kinetic_excess_fraction']}
        r = mine['rates']['h=0.001']
        cmp['dE_X/dt'] = (r['E_X'], ref['fixed_region_dE_X'], r['E_X'] / ref['fixed_region_dE_X'] - 1)
        cmp['dE_phi/dt'] = (r['E_phi'], ref['fixed_region_dE_phi'], r['E_phi'] / ref['fixed_region_dE_phi'] - 1)
        cmp['dQ/dt'] = (r['Q'], ref['fixed_region_dQ'], r['Q'] / ref['fixed_region_dQ'] - 1)
        res['regions'].append(dict(candidate=cid, region=name, mine=mine, compare_mine_published_relative=cmp))
        print(cid, name, {k: f'{v[2]:.2e}' for k, v in cmp.items()}, flush=True)
eX, ep, j = densities(q, p, d)
res['whole_box'] = dict(E_X=eX.sum() * d ** 3, E_phi=ep.sum() * d ** 3, Q=j.sum() * d ** 3,
                        released_false_vacuum=LP / 4 * V ** 4 * L ** 3)
json.dump(res, open(out, 'w'), indent=1, default=float)
