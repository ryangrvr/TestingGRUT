"""st_physics skeptic, check k4b: analysis of our physical-variable runs (k4_phys_tangent.py).

(1) Independent reproduction of Codex's published amplitude/charge responses (k = 2.0, 2.13, 2.4, 2.8) to tau 300.
(2) Growing-mode weight of a genuine 1-rad X rotation relative to one supplied preparation-phase radian.
(3) Vacuum-equivalent supplied radians per ln k, for four initial covariances at a = 1 (sensitivity test of the
    reviewer's single choice):
      A  adiabatic in the eigenbasis of M(0), w^2 = k^2 + m_i^2 (the reviewer's choice);
      B  massless Bunch-Davies for all six real fields, w = k;
      C  Codex's own modes.diagnostics regulator, w^2 = k^2 + max(m_i^2, 0) + 1;
      D  X only (massless), Y and Sigma unexcited.
(4) Band RMS (trapezoid in ln k over 1.8-2.8) first sampled crossings of 0.1, 0.3, 1 for amplitude and charge.

usage: python3 -B k4b_analysis.py <codex reality_framework dir> <k4 dir> <out.json>
"""
import sys, json
import numpy as np

RF, K4, OUT = sys.argv[1:4]
sys.path.insert(0, RF + '/grut_reheating_preparation_20261026/code')
import reheating as Rh  # noqa: E402

pkt = RF + '/grut_source_timing_20261026'
rr = {r['tag']: r for r in json.load(open(pkt + '/RESPONSE_RESULTS.json'))['rows']}
res = {}


def jac(zr, h=1e-2):
    E = np.eye(6) * h
    out = np.empty((6, 6))
    def f(v):
        z = v[0::2] + 1j * v[1::2]; _, F = Rh.potential_force(z)
        o = np.empty(6); o[0::2] = F.real; o[1::2] = F.imag; return o
    for j in range(6):
        out[:, j] = (8 * (f(zr + E[j]) - f(zr - E[j])) - (f(zr + 2 * E[j]) - f(zr - 2 * E[j]))) / (12 * h)
    return out


for tag, beta in (('quiet3', 3.), ('quiet5', 5.)):
    d = np.load(f'{K4}/k4_{tag}.npz')
    t, bg, T, colk, KV = d['times'], d['bg'], d['T'], d['colk'], d['KV']
    pub = np.load(f'{pkt}/{tag}.npz')
    a = np.exp(bg[:, 15]); H = bg[:, 16]
    zX = bg[:, 2] + 1j * bg[:, 3]; zXd = bg[:, 8] + 1j * bg[:, 9]
    qX = np.sqrt(2) * a * zX; pX = np.sqrt(2) * a * a * (zXd + H * zX)
    # background agreement with the published history at common times
    common = np.intersect1d(np.round(t, 6), np.round(pub['times'], 6))
    it = np.searchsorted(np.round(t, 6), common); ip = np.searchsorted(np.round(pub['times'], 6), common)
    bg_dx = float(np.max(abs(bg[it, 2:4] - pub['bg'][ip, 2:4])) / np.max(abs(pub['bg'][ip, 2:4])))
    bg_da = float(np.max(abs(bg[it, 15] - pub['bg'][ip, 15])))
    # conformal responses of every column
    dz = T[:, :, 0] + 1j * T[:, :, 1]; dzd = T[:, :, 6] + 1j * T[:, :, 7]
    dq = np.sqrt(2) * a[:, None] * dz; dp = np.sqrt(2) * (a * a)[:, None] * (dzd + H[:, None] * dz)
    A = np.real(dq * np.conj(qX)[:, None]) / abs(qX[:, None]) ** 2
    norm = abs(rr[tag]['late_NX_normalization_through80'])
    Q = (np.imag(np.conj(qX)[:, None] * dp) + np.imag(np.conj(dq) * pX[:, None])) / norm
    # supplied seed in the (q, p) basis at a = 1 (q = sqrt2 dz, p = q since H = 1, zdot = 0)
    v0 = d['v0']; dzs = d['dz_sup']
    qs = np.sqrt(2) * dzs; ps = qs.copy(); seed = np.r_[qs, ps]
    xr = np.zeros(12); xr[1] = np.sqrt(2) * v0[2]; xr[7] = xr[1]       # 1-rad X rotation, X real at a = 1
    out = dict(background_X_maxrel_vs_published=bg_dx, background_lna_maxabs_vs_published=bg_da,
               seed_components=dict(abs_chi_X=float(np.sqrt(2) * abs(v0[2] + 1j * v0[3])),
                                    dlnX=float(np.real(dzs[0] + 1j * dzs[1]) / v0[2]),
                                    dthetaX=float(np.imag((dzs[0] + 1j * dzs[1]) / v0[2])),
                                    dthetaY=float(np.imag((dzs[2] + 1j * dzs[3]) / (v0[4] + 1j * v0[5]))),
                                    dthetaS=float(np.imag((dzs[4] + 1j * dzs[5]) / (v0[6] + 1j * v0[7])))))
    late = (t >= 200) & (t <= 300)
    perk = []
    # initial mass matrix (own Jacobian; R from the model geometry)
    R0 = Rh.geometric(v0, beta, Rh.XI, Rh.ETA2)[6]
    M0 = jac(v0[2:8]) + Rh.XI * R0 * np.diag([1, 1, 0, 0, 0, 0.]) - R0 / 6 * np.eye(6)
    M0 = (M0 + M0.T) / 2
    lam, U = np.linalg.eigh(M0)
    out['M0_eigenvalues'] = lam.tolist()
    Pk = {c: np.zeros((len(t), len(KV))) for c in 'ABCD'}; PkQ = {c: np.zeros((len(t), len(KV))) for c in 'ABCD'}
    for ik, k in enumerate(KV):
        cols = np.flatnonzero(np.isclose(colk, k))[:12]
        Ak, Qk = A[:, cols], Q[:, cols]
        sup = Ak @ seed; supQ = Qk @ seed; xrot = Ak @ xr
        row = dict(k=float(k))
        jp = np.flatnonzero(np.isclose(pub['ks'], k))
        if len(jp):
            ref = pub['amplitude_per_radian'][:, jp[0]]; refQ = pub['charge_fraction_per_radian'][:, jp[0]]
            m = np.isin(np.round(pub['times'], 6), common) & (pub['times'] >= 20)
            mm = np.isin(np.round(t, 6), np.round(pub['times'][m], 6))
            row['amp_vs_published_maxrel'] = float(np.max(abs(sup[mm] - ref[m])) / np.max(abs(ref[m])))
            row['chg_vs_published_maxrel'] = float(np.max(abs(supQ[mm] - refQ[m])) / np.max(abs(refQ[m])))
        if np.isclose(k, 2.13):
            row['explicit_vs_combination_maxrel'] = float(np.max(abs(A[:, -2] - sup)) / np.max(abs(sup)))
            row['explicit_xrot_vs_combination_maxrel'] = float(np.max(abs(A[:, -1] - xrot)) / np.max(abs(xrot)))
        env_sup = np.max(abs(sup[late])); env_x = np.max(abs(xrot[late]))
        row['Xrot_1rad_over_supplied_1rad'] = float(env_x / env_sup)
        covs = {}
        w = np.sqrt(np.maximum(k * k + lam, 1e-12))
        covs['A'] = (np.block([[(U / (2 * w)) @ U.T, np.zeros((6, 6))], [np.zeros((6, 6)), (U * (w / 2)) @ U.T]]), bool(np.all(k * k + lam > 0)))
        covs['B'] = (np.diag(np.r_[np.full(6, 1 / (2 * k)), np.full(6, k / 2)]), True)
        wc = np.sqrt(k * k + np.maximum(lam, 0) + 1)
        covs['C'] = (np.block([[(U / (2 * wc)) @ U.T, np.zeros((6, 6))], [np.zeros((6, 6)), (U * (wc / 2)) @ U.T]]), True)
        CD = np.zeros((12, 12)); CD[0, 0] = CD[1, 1] = 1 / (2 * k); CD[6, 6] = CD[7, 7] = k / 2
        covs['D'] = (CD, True)
        pref = k ** 3 / (2 * np.pi ** 2)
        for c, (Cv, ok) in covs.items():
            P = pref * np.einsum('ti,ij,tj->t', Ak, Cv, Ak); PQ = pref * np.einsum('ti,ij,tj->t', Qk, Cv, Qk)
            Pk[c][:, ik] = P; PkQ[c][:, ik] = PQ
            row[f'vac_{c}_equiv_supplied_rad_per_lnk'] = float(np.max(np.sqrt(P[late])) / env_sup)
            row[f'vac_{c}_valid'] = ok
        row['naive_X_phase_per_lnk'] = float(k / (2 * np.pi) / out['seed_components']['abs_chi_X'])
        row['naive_times_weight'] = row['naive_X_phase_per_lnk'] * row['Xrot_1rad_over_supplied_1rad']
        perk.append(row)
    out['per_k'] = perk
    lk = np.log(KV)
    band = {}
    for c in 'ABCD':
        sa = np.sqrt(np.trapezoid(Pk[c], lk, axis=1)); sq = np.sqrt(np.trapezoid(PkQ[c], lk, axis=1))
        for name, s in (('amp', sa), ('chg', sq)):
            for lev in (0.1, 0.3, 1.0):
                wh = np.flatnonzero((s >= lev) & (t >= 1))
                band[f'{c}_{name}_{lev}'] = None if not len(wh) else float(t[wh[0]])
    out['band_rms_first_crossing_tau'] = band
    res[tag] = out
    print(tag, 'bg', bg_dx, bg_da, out['seed_components'], 'M0', np.round(lam, 3))
    for row in perk:
        print('  ', {k_: (round(v, 4) if isinstance(v, float) else v) for k_, v in row.items()})
    print('  band', band)
json.dump(res, open(OUT, 'w'), indent=1)
