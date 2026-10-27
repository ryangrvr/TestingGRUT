#!/usr/bin/env python3
"""st_clocks_skeptic (Claude, own code): adversarial re-derivation of the decay widths used in Codex b73de64 and checked by
track st_clocks.  DIFFERENT METHOD from st_clocks/clocks.py:

* st_clocks wrote the amplitude by hand in the complex basis ((kappa/sqrt2) * 3! per channel, XXX and XbarXbarXbar) and
  checked only the Dalitz AREA by dblquad.
* Here: sympy expands the model potential (the rp/hw_gates spot-point potential, which is the model definition, not
  code) in REAL canonical fields X = (X1+iX2)/sqrt2, Y = (Y1+iY2)/sqrt2, Sigma = e^{i th}(h+s+iG)/sqrt2, diagonalises the Y
  mass matrix numerically, takes 4th/3rd derivatives along the eigen-directions to get every real-field vertex, and sums
  widths over every REAL final state with its own identical-particle factor.  Phase space by recursive factorisation
  Phi3 = int ds/(2pi) Phi2(M; 0, sqrt s) Phi2(sqrt s; 0, 0) (scipy quad), not the Dalitz area.
* Vector width by an explicit polarisation-vector sum and an angular quadrature in the A rest frame (no closed form).
* Two Sigma phases (th = -pi/2 as in the engine preparation, th = 0) to show phase independence.
No Codex or st_clocks code is imported.  Output: ../results/widths_skeptic.json
"""
import json, math
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'results' / 'widths_skeptic.json'
H_GEV = 133349619398.33566                         # [Q] source_interface meta.fixed_point.H_GeV
MX_GEV = 0.09240896373753621                       # [Q] fixed_point.mX_GeV
PAR = dict(lx=1e-3, ly=1.0, ls=1.0, kap=0.0034586847211042523, M=30.0, h=40.0, mu=7.159456159513794)   # [Q] spot point

X1, X2, Y1, Y2, s, G = sp.symbols('X1 X2 Y1 Y2 s G', real=True)
FIELDS = [X1, X2, Y1, Y2, s, G]


def potential(th):
    lx, ly, ls, kap, M, h, mu = [sp.Float(PAR[k], 30) for k in
                                 ('lx', 'ly', 'ls', 'kap', 'M', 'h', 'mu')]
    X = (X1 + sp.I * X2) / sp.sqrt(2)
    Y = (Y1 + sp.I * Y2) / sp.sqrt(2)
    S = sp.exp(sp.I * th) * (h + s + sp.I * G) / sp.sqrt(2)
    absq = lambda z: sp.expand(z * sp.conjugate(z))
    V = lx * absq(X) ** 2 + M ** 2 * absq(Y) + ly * absq(Y) ** 2 + ls * (absq(S) - h ** 2 / 2) ** 2 \
        + 2 * sp.re(sp.expand(kap * sp.conjugate(Y) * X ** 3 + mu * sp.conjugate(S) * Y * Y))
    return sp.expand(V)


def phi2_massless_pair(M, m1=0.0):
    """Lorentz-invariant two-body phase space for daughters (m1, 0): Phi2 = |p|/(16 pi^2 M) * 4 pi = |p|/(4 pi M);
    massless pair -> 1/(8 pi).  Identical-particle factors are applied by the caller."""
    p = (M * M - m1 * m1) / (2 * M)
    return p / (4 * math.pi * M)


def phi3_recursive(M):
    # Phi3 = int_0^{M^2} ds/(2 pi) Phi2(M; 0, sqrt s) Phi2(sqrt s; 0, 0)
    f = lambda ss: phi2_massless_pair(M, math.sqrt(ss)) * (1.0 / (8 * math.pi)) / (2 * math.pi) if ss > 0 else \
        phi2_massless_pair(M, 0.0) / (8 * math.pi) / (2 * math.pi)
    val, err = quad(f, 0.0, M * M, epsabs=0, epsrel=1e-13, limit=200)
    return val


def vertex_tensor(V, dirs_up, dirs_X, order):
    """derivative of V along direction vectors, evaluated at the vacuum (all fluctuations zero)."""
    expr = V
    for d in dirs_up + dirs_X:
        expr = sum(c * sp.diff(expr, f) for c, f in zip(d, FIELDS) if c != 0)
    return float(expr.subs({f: 0 for f in FIELDS}))


def unit(i):
    e = [0] * 6; e[i] = 1; return e


def analyse(th):
    V = potential(th)
    vac = {f: 0 for f in FIELDS}
    # Y mass matrix (real basis) at the vacuum
    MY = np.array([[float(sp.diff(V, a, b).subs(vac)) for b in (Y1, Y2)] for a in (Y1, Y2)])
    w, U = np.linalg.eigh(MY)
    mL, mH = math.sqrt(w[0]), math.sqrt(w[1])
    uL = [0, 0, U[0, 0], U[1, 0], 0, 0]; uH = [0, 0, U[0, 1], U[1, 1], 0, 0]
    ms2 = float(sp.diff(V, s, 2).subs(vac)); mG2 = float(sp.diff(V, G, 2).subs(vac))
    # Y_L -> three real X quanta: sum over unordered final states with identical-particle factors
    phi3 = phi3_recursive(mL)
    states = [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)]
    GL = 0.0; per_state = {}
    for st in states:
        T = vertex_tensor(V, [uL], [unit(i) for i in st], 4)
        sym = math.factorial(3) if len(set(st)) == 1 else (2 if len(set(st)) == 2 else 1)
        g = T * T * phi3 / (2 * mL) / sym
        per_state['X%d X%d X%d' % tuple(i + 1 for i in st)] = dict(vertex=T, symmetry=sym, gamma=g)
        GL += g
    # heavy mode -> 3X, same method
    phi3H = phi3_recursive(mH); GH3 = 0.0
    for st in states:
        T = vertex_tensor(V, [uH], [unit(i) for i in st], 4)
        sym = math.factorial(3) if len(set(st)) == 1 else (2 if len(set(st)) == 2 else 1)
        GH3 += T * T * phi3H / (2 * mH) / sym
    # condensate-assisted Y_L -> X X with one leg on a homogeneous condensate X_c (random phases): Gamma_2 / |X_c|^2
    rng = np.random.default_rng(1)
    ratios = []
    T4 = {}
    for i in range(2):
        for j in range(2):
            for k in range(2):
                T4[(i, j, k)] = vertex_tensor(V, [uL], [unit(i), unit(j), unit(k)], 4)
    for _ in range(6):
        ph = rng.uniform(0, 2 * math.pi); Xc = 1.0                         # |X_c| = 1 (complex canonical)
        xc = [math.sqrt(2) * Xc * math.cos(ph), math.sqrt(2) * Xc * math.sin(ph)]   # X = (X1+iX2)/sqrt2
        G2 = 0.0
        for (j, k) in [(0, 0), (0, 1), (1, 1)]:
            B = sum(T4[(i, j, k)] * xc[i] for i in range(2))
            sym = 2 if j == k else 1
            G2 += B * B * phi2_massless_pair(mL) / (2 * mL) / sym
        ratios.append(G2)
    c2 = float(np.mean(ratios)); c2_spread = float(np.ptp(ratios) / np.mean(ratios))
    # cascade Y_H -> Y_L + G (Goldstone of Sigma = longitudinal A at E >> m_A)
    VHLG = vertex_tensor(V, [uH], [uL, unit(5)], 3)
    pc = (mH ** 2 - mL ** 2) / (2 * mH)
    Gcasc = VHLG ** 2 * phi2_massless_pair(mH, mL) / (2 * mH)
    # Sigma radial decays
    msig = math.sqrt(ms2)
    VsGG = vertex_tensor(V, [unit(4)], [unit(5), unit(5)], 3)
    GsAA = VsGG ** 2 * phi2_massless_pair(msig) / (2 * msig) / 2
    VsLL = vertex_tensor(V, [unit(4)], [uL, uL], 3)
    pLL = math.sqrt(msig ** 2 / 4 - mL ** 2)
    GsLL = VsLL ** 2 * (pLL / (4 * math.pi * msig)) / (2 * msig) / 2
    return dict(theta_Sigma=float(th), mL=mL, mH=mH, m_sigma=msig, m_Goldstone2=mG2, eigvec_L=U[:, 0].tolist(),
                phi3_recursive_over_m2_256pi3=phi3 / (mL ** 2 / (256 * math.pi ** 3)),
                gamma_L=GL, gamma_L_over_formula=GL / (3 * PAR['kap'] ** 2 * mL / (256 * math.pi ** 3)),
                Gamma_L_GeV=GL * H_GEV, per_real_state=per_state, gamma_H_to_3X=GH3,
                cond_assisted_coeff=c2, cond_coeff_over_9k2_8pi_mL=c2 / (9 * PAR['kap'] ** 2 / (8 * math.pi * mL)),
                cond_phase_spread=c2_spread,
                vertex_YH_YL_G=VHLG, vertex_over_sqrt2_mu=VHLG / (math.sqrt(2) * PAR['mu']), gamma_cascade=Gcasc,
                Gamma_cascade_GeV=Gcasc * H_GEV, E_A_cascade=pc, gamma_sigma_AA=GsAA, gamma_sigma_YLYL=GsLL)


def vector_width(e, qS=6, h=PAR['h']):
    """A -> X Xbar: explicit polarisation vectors and angular quadrature in the A rest frame (H_i units)."""
    absS2 = h * h / 2                                   # |Sigma|^2 in vacuum (V has (|S|^2 - h^2/2)^2)
    mA = math.sqrt(2 * (qS * e) ** 2 * absS2)           # |D S|^2 = ... + (q e)^2 A^2 |S|^2 = (1/2) mA^2 A^2
    mx = MX_GEV / H_GEV
    p = math.sqrt(mA * mA / 4 - mx * mx); E = mA / 2
    eta = np.diag([1.0, -1.0, -1.0, -1.0])
    pols = [np.array([0, 1.0, 0, 0]), np.array([0, 0, 1.0, 0]), np.array([0, 0, 0, 1.0])]   # rest frame
    xg, wg = np.polynomial.legendre.leggauss(40)
    nph = 64; phis = 2 * math.pi * (np.arange(nph) + 0.5) / nph
    tot = 0.0
    for c, wc in zip(xg, wg):
        st = math.sqrt(1 - c * c)
        for ph in phis:
            n = np.array([st * math.cos(ph), st * math.sin(ph), c])
            p1 = np.r_[E, p * n]; p2 = np.r_[E, -p * n]
            amp2 = sum((eps @ eta @ (p1 - p2)) ** 2 for eps in pols) / 3.0 * e * e
            tot += wc * (2 * math.pi / nph) * amp2
    # Gamma = 1/(2 mA) * int dPhi2 |M|^2,  dPhi2 = |p| dOmega / (16 pi^2 mA)
    G = tot * p / (16 * math.pi ** 2 * mA) / (2 * mA)
    beta = 2 * p / mA
    return dict(e=e, mA=mA, mA_GeV=mA * H_GEV, gamma=G, Gamma_GeV=G * H_GEV,
                over_formula=G / (e * e * mA * beta ** 3 / (48 * math.pi)), e3h_over_8pi=e ** 3 * h / (8 * math.pi))


def main():
    res = dict(scalar={}, vector={})
    for name, th in [('theta_minus_pi_over_2', -sp.pi / 2), ('theta_0', sp.Integer(0))]:
        res['scalar'][name] = analyse(th)
        print(name, {k: v for k, v in res['scalar'][name].items() if k != 'per_real_state'}, flush=True)
    for e in (1e-7, 1e-6):
        res['vector']['e%g' % e] = vector_width(e)
        print(res['vector']['e%g' % e])
    g7, g6 = res['vector']['e1e-07']['gamma'], res['vector']['e1e-06']['gamma']
    res['vector_median_ratio'] = g6 / g7
    gL = res['scalar']['theta_minus_pi_over_2']['gamma_L']
    res['medians'] = dict(YL=math.log(2) / gL, A_e1e_7=math.log(2) / g7, A_e1e_6=math.log(2) / g6,
                          YL_1percent=-math.log1p(-0.01) / gL)
    res['codex'] = dict(Gamma_L_GeV=13413.681870567107, Gamma_A_e1e_7=2.1223251086667907e-10,
                        Gamma_A_e1e_6=2.1223251086668018e-07, tau_half_YL=6890793.564853669,
                        tau_half_A_e1e_7=4.3551721806072263e+20, tau_half_A_e1e_6=4.355172180607203e+17)
    res['rel_vs_codex'] = dict(
        Gamma_L=res['scalar']['theta_minus_pi_over_2']['Gamma_L_GeV'] / res['codex']['Gamma_L_GeV'] - 1,
        Gamma_A_e1e_7=res['vector']['e1e-07']['Gamma_GeV'] / res['codex']['Gamma_A_e1e_7'] - 1,
        Gamma_A_e1e_6=res['vector']['e1e-06']['Gamma_GeV'] / res['codex']['Gamma_A_e1e_6'] - 1,
        tau_half_YL=res['medians']['YL'] / res['codex']['tau_half_YL'] - 1,
        tau_half_A_e1e_7=res['medians']['A_e1e_7'] / res['codex']['tau_half_A_e1e_7'] - 1,
        tau_half_A_e1e_6=res['medians']['A_e1e_6'] / res['codex']['tau_half_A_e1e_6'] - 1)
    print(json.dumps(res['rel_vs_codex'], indent=1), res['vector_median_ratio'])
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2) + '\n')


if __name__ == '__main__':
    main()
