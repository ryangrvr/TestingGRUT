"""Post-waterfall history with Pauli blocking of the axial channels: K, radiation fractions, SM plasma time [E on C].

The track's fluid (reheat.py) applies the half-period production rate Gamma_j = (m/pi) eps_j(A)/(S^2/2) at every
half-period, which assumes empty fermion modes each time. sk_axial.py shows (exact two-level evolution, both helicities,
constant amplitude) that the band saturates after ~3-6 half-periods at F_sat = 0.13 (chi->N2), 0.085 (N3), 0.055 (N4)
of rho_chi and 0.145 (sigma->N1) of rho_sigma at the slice amplitude; F_sat(A) is tabulated at A = 10, 20, 42.3 (3, 6,
12.3 for N1) in ../results/sk_axial_Fsat_vs_A.json and interpolated in ln A (slope 2 outside).

My fluid (independent code; the production tables eps_j(A) are the track's axial_tables.json, which I reproduced at
10^-7 at 12 points):
  condensates rho_sig, rho_chi (matter); per species j: N energy E_j, number n_j (mean energy E_j/n_j), in-band
  blocking energy B_j; SM radiation.
  production P_j = Gamma_j rho_c max(0, 1 - B_j/(F_sat,j rho_c)),
  dB_j/dt = P_j - (Gamma_lab,j + Gamma_th) B_j - (6H + 2 Gamma_c) B_j   (redshift + band-edge shrinkage),
  Gamma_lab,j = Gamma_rest,j M_j / max(M_j, <E>),  Gamma_rest = m_nu M^2/(8 pi v^2), m_nu = 0.05 eV (track's seesaw),
  N -> L H feeds the SM.  Gamma_th is an unknown N-thermalisation rate (0 baseline; 1 and 10 H_i brackets).
Validation: with blocking off ('unblocked') it reproduces the track's K = 2.05 starting from the bench state at t = 0.19.
Start: bench_HFa state at t = 0.19 H_i^-1 (eps_H ~ 0.45; the track's own friction has moved only 0.9% to radiation).
Usage: python3 -B sk_pauli.py -> ../results/sk_pauli.json
"""
import math, json, os
import numpy as np
from scipy.integrate import solve_ivp
import sk_common as C

mP = C.mP
HGEV = C.H_REF
MNU, VEW = 0.05e-9, 174.1
TAB = json.load(open(os.path.join(C.BENCH, 'axial_tables.json')))
FS = json.load(open(os.path.join(C.OUT, 'sk_axial_Fsat_vs_A.json')))
AXM = json.load(open(os.path.join(C.OUT, 'sk_axial.json')))['multi_period']


def fsat_table(tag, mh):
    pts = [(v['A'], v['Fsat']) for k, v in FS.items() if k.startswith(tag + '_')]
    key = {'N1': 'sigma->N1', 'N2': 'chi->N2', 'N3': 'chi->N3', 'N4': 'chi->N4'}[tag]
    rows = AXM[key]['rows']
    pts.append((AXM[key]['A'], float(np.mean([r['F_tot'] for r in rows[6:]]))))
    pts.sort()
    return np.log([p[0] for p in pts]), np.log([p[1] for p in pts])


class Chan:
    def __init__(self, tag, m, M, f):
        self.m, self.M, self.f = m, M, f
        key = min(TAB.keys(), key=lambda s: abs(float(s) - M / m))
        self.lA, self.le = np.log(TAB[key]['A']), np.log(np.maximum(TAB[key]['eps'], 1e-300))
        self.flA, self.flF = fsat_table(tag, M / m)
        beta2 = 1 - 4 * M * M / (m * m)
        self.Gp = (M / f) ** 2 * m * math.sqrt(beta2) / (16 * math.pi) if beta2 > 0 else 0.0
        self.Grest = MNU * (M * HGEV) ** 2 / (8 * math.pi * VEW ** 2) / HGEV

    def amp(self, rho):
        Sig = math.sqrt(2 * max(rho, 0.0)) / self.m
        return Sig / (2 * self.f), Sig / self.m

    def rate(self, rho):
        if rho <= 0:
            return 0.0
        A, S = self.amp(rho)
        lA = math.log(A)
        if lA < self.lA[0]:
            return self.Gp
        le = float(np.interp(lA, self.lA, self.le)) if lA <= self.lA[-1] else self.le[-1] + 4 * (lA - self.lA[-1])
        return (self.m / math.pi) * math.exp(le) / (0.5 * S * S)

    def fsat(self, rho):
        A, _ = self.amp(rho)
        lA = math.log(max(A, 1e-30))
        if lA < self.flA[0]:
            lF = self.flF[0] + 2 * (lA - self.flA[0])
        elif lA > self.flA[-1]:
            lF = self.flF[-1] + 2 * (lA - self.flA[-1])
        else:
            lF = float(np.interp(lA, self.flA, self.flF))
        return min(math.exp(lF), 1.0)

    def mu5(self, rho):
        return math.sqrt(2 * max(rho, 0.0)) / (2 * self.f)


def run(blocking=True, Gth=0.0, t_end=3e5, start=None, chans=None, pert_free=True):
    cs = [Chan('N1', 600., 250., 1200.)]
    cc = [Chan('N2', 250., 400., 2000.), Chan('N3', 250., 200., 2000.), Chan('N4', 250., 120., 2000.)]
    allc = [(c, 0) for c in cs] + [(c, 1) for c in cc]
    st = start or dict(t=0.19, a=1.209, H=0.9985, fs=0.054, fc=0.937, fr=0.009)
    rho0 = 3 * mP ** 2 * st['H'] ** 2
    a0 = st['a']
    nC = len(allc)

    def unpack(y):
        la = y[0]
        a = math.exp(la)
        rs, rc = max(y[1], 0) / a ** 3, max(y[2], 0) / a ** 3
        E = [max(y[3 + j], 0) / a ** 4 for j in range(nC)]
        n = [max(y[3 + nC + j], 1e-300) / a ** 3 for j in range(nC)]
        B = [max(y[3 + 2 * nC + j], 0) / a ** 4 for j in range(nC)]
        rsm = max(y[3 + 3 * nC], 0) / a ** 4
        return a, rs, rc, E, n, B, rsm

    def rhs(t, y):
        a, rs, rc, E, n, B, rsm = unpack(y)
        rho = rs + rc + sum(E) + rsm
        H = math.sqrt(rho / (3 * mP ** 2))
        P, Gl = [], []
        for j, (c, w) in enumerate(allc):
            rcond = rs if w == 0 else rc
            G = c.rate(rcond)
            if blocking:
                Es = c.fsat(rcond) * rcond
                bl = max(0.0, 1.0 - B[j] / Es) if Es > 0 else 0.0
                # perturbative products (k ~ m/2) are inside the band while mu5 > m/2: blocked alike
            else:
                bl = 1.0
            Pnp = G * rcond * bl
            # perturbative tail (2M < m): products at k = (m/2) beta; 'pert_free' treats them as unblocked
            # (N thermalised or removed fast at the shell), otherwise they are blocked like the band
            Pp = c.Gp * rcond * (1.0 if (pert_free or not blocking) else bl)
            P.append(max(Pnp, Pp))
            Emean = E[j] / n[j] if n[j] > 0 else c.M
            Gl.append(c.Grest * c.M / max(c.M, Emean))
        Gs = sum(P[j] for j, (c, w) in enumerate(allc) if w == 0) / max(rs, 1e-300)
        Gc = sum(P[j] for j, (c, w) in enumerate(allc) if w == 1) / max(rc, 1e-300)
        dy = [H, -Gs * rs * a ** 3, -Gc * rc * a ** 3]
        dE, dn, dB = [], [], []
        for j, (c, w) in enumerate(allc):
            rcond = rs if w == 0 else rc
            Einj = max(0.75 * c.mu5(rcond), c.M)
            Gband = Gl[j] if not blocking else c.Grest * c.M / max(c.M, 0.75 * c.mu5(rcond))
            dE.append((P[j] - Gl[j] * E[j]) * a ** 4)
            dn.append((P[j] / Einj - Gl[j] * n[j]) * a ** 3)
            Gcw = Gs if w == 0 else Gc
            dB.append((P[j] - (Gband + Gth + 2 * H + 2 * Gcw) * B[j]) * a ** 4)
        dsm = sum(Gl[j] * E[j] for j in range(nC)) * a ** 4
        return dy + dE + dn + dB + [dsm]

    y0 = [math.log(a0), st['fs'] * rho0 * a0 ** 3, st['fc'] * rho0 * a0 ** 3]
    Er0 = st['fr'] * rho0
    # initial radiation split as chi->N2 (dominant early channel)
    y0 += [(Er0 if j == 1 else 0.0) * a0 ** 4 for j in range(nC)]
    y0 += [(Er0 / max(allc[j][0].M, 0.75 * allc[j][0].mu5(st['fc'] * rho0)) if j == 1 else 1e-300) * a0 ** 3 for j in range(nC)]
    y0 += [(Er0 if j == 1 else 0.0) * a0 ** 4 for j in range(nC)]
    y0 += [0.0]
    ts = np.concatenate([np.linspace(st['t'], 3.0, 600), np.geomspace(3.0, t_end, 600)[1:]])
    sol = solve_ivp(rhs, (st['t'], t_end), y0, method='LSODA', rtol=1e-7, atol=1e-30, t_eval=ts, max_step=0.5)
    rec = dict(t=[], a=[], H=[], frad=[], fsm=[], fchi=[], fsig=[], K=[])
    for i, t in enumerate(sol.t):
        a, rs, rc, E, n, B, rsm = unpack(sol.y[:, i])
        rho = rs + rc + sum(E) + rsm
        H = math.sqrt(rho / (3 * mP ** 2))
        rec['t'].append(float(t)); rec['a'].append(a); rec['H'].append(H)
        rec['frad'].append((sum(E) + rsm) / rho); rec['fsm'].append(rsm / rho)
        rec['fchi'].append(rc / rho); rec['fsig'].append(rs / rho); rec['K'].append(a * a * H)
    T = np.array(rec['t'])

    def at(key, tt):
        return float(np.interp(tt, T, rec[key]))

    def first(key, thr):
        v = np.array(rec[key])
        i = int(np.argmax(v >= thr))
        return float(T[i]) if (v >= thr).any() else None
    summ = dict(blocking=blocking, Gamma_th=Gth, pert_free=pert_free, K=rec['K'][-1], frad_t03=at('frad', 0.3), frad_t06=at('frad', 0.6),
                frad_t1=at('frad', 1.0), frad_t3=at('frad', 3.0), frad_t10=at('frad', 10.0), frad_t100=at('frad', 100.0),
                t_SM_1em3=first('fsm', 1e-3), t_SM_half=first('fsm', 0.5), t_rad_99=first('frad', 0.99),
                fchi_t1=at('fchi', 1.0), fsig_t1=at('fsig', 1.0))
    return summ, rec


if __name__ == '__main__':
    out = dict(label='E on C (skeptic fluid with Pauli blocking)', note=__doc__, cases={})
    for nm, kw in [('unblocked_validation', dict(blocking=False)),
                   ('pauli_Gth0_pertfree', dict(blocking=True, Gth=0.0)),
                   ('pauli_Gth1_pertfree', dict(blocking=True, Gth=1.0)),
                   ('pauli_Gth10_pertfree', dict(blocking=True, Gth=10.0)),
                   ('pauli_Gth0_pertblocked', dict(blocking=True, Gth=0.0, pert_free=False))]:
        s, rec = run(**kw)
        out['cases'][nm] = s
        print(nm, json.dumps(s), flush=True)
        if nm == 'pauli_Gth0_pertfree':
            idx = np.unique(np.linspace(0, len(rec['t']) - 1, 120).astype(int))
            out['history_pauli_Gth0_pertfree'] = {k: [float('%.6g' % rec[k][i]) for i in idx] for k in rec}
    C.dump('sk_pauli.json', out)
