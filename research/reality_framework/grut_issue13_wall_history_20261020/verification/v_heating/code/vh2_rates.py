#!/usr/bin/env python3
"""v_heating check 2: gravitational (s-channel graviton) freeze-in SM bath -> real massless scalars.

Independent verifier code.  Metric (+,-,-,-).  Canonical graviton: g = eta + 2h/M_P (reduced M_P),
L_int = -(1/M_P) h_mn T^mn, de Donder propagator i P/k^2, P = (eta eta + eta eta - eta eta)/2.
=> M = [T_i.T_f - (1/2) trT_i trT_f] / (M_P^2 s)  (overall phase irrelevant).
Newton-limit check: two static masses give V = -m1 m2/(8 pi M_P^2 r).

Amplitudes computed numerically with explicit spinors (chiral rep, Weyl projector) and explicit transverse
polarization vectors, then integrated over the CM angle.  Rates use Bose/Fermi initial distributions.

Outputs: kappa_i  with  (1/2)_final * int dPhi2 Sum|M|^2 = kappa_i s^2/M_P^4 per initial species (all internal
states of one species summed: 1 real scalar, 1 Weyl field incl. its antiparticle, 1 gauge boson);
pair rate Gamma, energy rate, Monte Carlo <E>, <1/E> of produced quanta.
"""
import json, os
import numpy as np
from scipy import integrate, special

OUT = os.path.join(os.path.dirname(__file__), '..', 'results')
eta = np.diag([1.0, -1.0, -1.0, -1.0])
Z3, Z5 = special.zeta(3), special.zeta(5)

# ---------- Dirac algebra (chiral representation) ----------
s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
sig = [sx, sy, sz]
Zero = np.zeros((2, 2), dtype=complex)
gam = [np.block([[Zero, s0], [s0, Zero]])] + [np.block([[Zero, s], [-s, Zero]]) for s in sig]
g5 = np.block([[-s0, Zero], [Zero, s0]])
PL = (np.eye(4) - g5) / 2


def sqrtm_psd(M):
    w, U = np.linalg.eigh(M)
    w = np.clip(w, 0, None)
    return U @ np.diag(np.sqrt(w)) @ U.conj().T


def spinors(p):
    """u(p,s), v(p,s) for massless p (chiral rep, Peskin conventions)."""
    E = p[0]; pv = p[1:]
    ps = E * s0 - sum(pv[i] * sig[i] for i in range(3))      # p.sigma
    psb = E * s0 + sum(pv[i] * sig[i] for i in range(3))     # p.sigmabar
    A, B = sqrtm_psd(ps), sqrtm_psd(psb)
    us, vs = [], []
    for xi in (np.array([1, 0], complex), np.array([0, 1], complex)):
        us.append(np.concatenate([A @ xi, B @ xi]))
        vs.append(np.concatenate([A @ xi, -B @ xi]))
    return us, vs


def slash(p):
    return sum(eta[m, m] * p[m] * gam[m] for m in range(4))


def bar(w):
    return w.conj() @ gam[0]


def lower(T):
    return eta @ T @ eta


def contract(Ti, Tf):
    """T_i^{mn} P_{mnrs} T_f^{rs} with upper-index inputs."""
    TiL = lower(Ti)
    return np.sum(TiL * Tf) - 0.5 * np.trace(eta @ Ti) * np.trace(eta @ Tf)


def T_scalar(p1, p2, xi=0.0):
    """<0|T^{mn}|p1 p2> (or <k1 k2|T|0>) for a real scalar with nonminimal xi (improvement term)."""
    d = np.dot(p1, eta @ p2)
    P = p1 + p2; P2 = np.dot(P, eta @ P)
    etaU = np.linalg.inv(eta)
    return (np.outer(p1, p2) + np.outer(p2, p1) - etaU * d
            - 2 * xi * (np.outer(P, P) - etaU * P2))


def T_weyl(p1, p2, u1, v2):
    """(1/4) vbar(p2)[g^m (p1-p2)^n + g^n (p1-p2)^m] P_L u(p1)."""
    D = p1 - p2
    vb = bar(v2)
    cur = np.array([vb @ gam[m] @ PL @ u1 for m in range(4)])   # upper-index current
    return 0.25 * (np.outer(cur, D) + np.outer(D, cur))


def T_dirac(p1, p2, u1, v2):
    D = p1 - p2
    vb = bar(v2)
    cur = np.array([vb @ gam[m] @ u1 for m in range(4)])
    return 0.25 * (np.outer(cur, D) + np.outer(D, cur))


def Fmunu(p, e):
    return np.outer(p, e) - np.outer(e, p)     # upper indices


def T_vector(p1, e1, p2, e2):
    F1, F2 = Fmunu(p1, e1), Fmunu(p2, e2)
    # T^{mn} = F1^{m a} F2^{n}_{a} + F2^{m a} F1^{n}_{a} - (1/2) eta^{mn} F1^{ab} F2_{ab}
    F2lo = F2 @ eta          # F2^{n}_{a}
    F1lo = F1 @ eta
    t = F1 @ F2lo.T + F2 @ F1lo.T
    FF = np.sum(lower(F1) * F2)
    return t - 0.5 * np.linalg.inv(eta) * FF


def transverse(n):
    n = n / np.linalg.norm(n)
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(n, a); e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return [np.concatenate([[0.0], e1]), np.concatenate([[0.0], e2])]


def kin(E, th, ph=0.3):
    p1 = E * np.array([1, 0, 0, 1.0]); p2 = E * np.array([1, 0, 0, -1.0])
    n = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    k1 = E * np.concatenate([[1.0], n]); k2 = E * np.concatenate([[1.0], -n])
    return p1, p2, k1, k2


def M2sum(channel, th, xi_i=0.0, xi_f=0.0, E=1.0):
    """Sum over initial internal states of |T_i P T_f|^2 (M_P = 1, without the 1/s^2)."""
    p1, p2, k1, k2 = kin(E, th)
    Tf = T_scalar(k1, k2, xi_f)
    if channel == 'scalar':
        return abs(contract(T_scalar(p1, p2, xi_i), Tf))**2
    if channel in ('weyl', 'dirac'):
        us, _ = spinors(p1); _, vs = spinors(p2)
        tot = 0.0
        for u1 in us:
            for v2 in vs:
                Ti = T_weyl(p1, p2, u1, v2) if channel == 'weyl' else T_dirac(p1, p2, u1, v2)
                tot += abs(contract(Ti, Tf))**2
        return tot
    if channel == 'vector':
        tot = 0.0
        for e1 in transverse(p1[1:]):
            for e2 in transverse(p2[1:]):
                tot += abs(contract(T_vector(p1, e1, p2, e2), Tf))**2
        return tot
    raise ValueError


def kappa(channel, xi_i=0.0, xi_f=0.0):
    """(1/2) int dPhi2 Sum|M|^2 / s^2  with dPhi2 = dcos/(16 pi); s = 4E^2 with E = 1."""
    s = 4.0
    f = lambda c: M2sum(channel, np.arccos(c), xi_i, xi_f) / s**2
    val, _ = integrate.quad(f, -1, 1, epsabs=1e-14, epsrel=1e-12, limit=200)
    return 0.5 * val / (16 * np.pi) / s**2


def checks():
    out = {}
    # completeness sum_s u ubar = pslash (massless)
    p = np.array([2.0, 0.3, -1.1, np.sqrt(4 - 0.09 - 1.21)])
    us, vs = spinors(p)
    out['completeness_u'] = float(np.max(np.abs(sum(np.outer(u, bar(u)) for u in us) - slash(p))))
    out['completeness_v'] = float(np.max(np.abs(sum(np.outer(v, bar(v)) for v in vs) - slash(p))))
    # tracelessness of fermion and vector T, gauge invariance of vector T
    p1, p2, k1, k2 = kin(1.0, 0.7)
    us, _ = spinors(p1); _, vs = spinors(p2)
    out['trace_weyl'] = float(abs(np.trace(eta @ T_weyl(p1, p2, us[0], vs[1]))))
    e = transverse(p1[1:]); f = transverse(p2[1:])
    out['trace_vector'] = float(abs(np.trace(eta @ T_vector(p1, e[0], p2, f[1]))))
    out['gauge_inv_vector'] = float(np.max(np.abs(T_vector(p1, p1, p2, f[0]))))
    out['conservation_vector'] = float(np.max(np.abs((p1 + p2) @ eta @ T_vector(p1, e[0], p2, f[0]))))
    # Newton limit: static masses, t-channel with relativistic normalisation T^00 = 2m^2
    m1, m2 = 3.0, 5.0
    P1 = np.array([m1, 0, 0, 0.]); P2 = np.array([m2, 0, 0, 0.])
    T1 = np.outer(P1, P1) * 2 - np.linalg.inv(eta) * (m1**2 - m1**2)   # p p' + p' p - eta(p.p' - m^2), p=p'
    T2 = np.outer(P2, P2) * 2
    num = contract(T1, T2)           # = 2 m1^2 m2^2
    # M = num/(M_P^2 q^2) with q^2=-|q|^2 ; NR amplitude /(2m1 2m2): V(q) = -num/(4 m1 m2 |q|^2 M_P^2)
    coeff = num / (4 * m1 * m2)      # V(r) = -coeff/(4 pi r M_P^2) -> compare with G m1 m2 = m1 m2/(8 pi M_P^2)
    out['newton_ratio'] = float(coeff / (4 * np.pi) / (m1 * m2 / (8 * np.pi)))
    # analytic shapes: scalar t^2u^2/s^2, Weyl (t-u)^2 tu/(4 s^2), Dirac x2, vector 2 t^2u^2/s^2
    th = 1.1; s = 4.0
    t = -s * (1 - np.cos(th)) / 2; u = -s * (1 + np.cos(th)) / 2
    out['ratio_scalar_to_t2u2/s2'] = float(M2sum('scalar', th) / (t**2 * u**2 / s**2) * s**0)
    out['ratio_weyl_to_(t-u)^2tu/(4s^2)'] = float(M2sum('weyl', th) / ((t - u)**2 * t * u / (4 * s**2)))
    out['ratio_dirac_to_weyl'] = float(M2sum('dirac', th) / M2sum('weyl', th))
    out['ratio_vector_to_2t2u2/s2'] = float(M2sum('vector', th) / (2 * t**2 * u**2 / s**2))
    return out


def thermal_moments():
    J = {}
    J['B3'] = integrate.quad(lambda x: x**3 / np.expm1(x), 0, 200)[0]
    J['F3'] = integrate.quad(lambda x: x**3 / (np.exp(x) + 1), 0, 200)[0]
    J['B4'] = integrate.quad(lambda x: x**4 / np.expm1(x), 0, 200)[0]
    J['F4'] = integrate.quad(lambda x: x**4 / (np.exp(x) + 1), 0, 200)[0]
    return J


def sample_E(rng, n, stat):
    """sample x with density ~ x^3 f(x) (f Bose/Fermi), by inverse CDF on a fine grid."""
    x = np.linspace(1e-6, 60, 600001)
    f = 1 / np.expm1(x) if stat == 'B' else 1 / (np.exp(x) + 1)
    w = x**3 * f
    c = np.cumsum(w); c /= c[-1]
    return np.interp(rng.random(n), c, x)


def mc_final_energies(channel, stat, n, rng):
    """Produced-quanta energies (plasma frame, T = 1) for massless 2->2 graviton s-channel."""
    E1 = sample_E(rng, n, stat); E2 = sample_E(rng, n, stat)
    # relative angle weight (1-c)^2 on [-1,1]: sample y=1-c in [0,2] with density y^2 -> y = 2 U^{1/3}
    c12 = 1 - 2 * rng.random(n)**(1 / 3)
    # p1 along z, p2 at angle
    ph = 2 * np.pi * rng.random(n)
    s12 = np.sqrt(1 - c12**2)
    p1 = np.stack([E1, 0 * E1, 0 * E1, E1], 1)
    p2 = np.stack([E2, E2 * s12 * np.cos(ph), E2 * s12 * np.sin(ph), E2 * c12], 1)
    P = p1 + p2
    beta = P[:, 1:] / P[:, :1]
    # boost p1 to CM
    def boost(p, b):
        b2 = np.sum(b * b, 1); g = 1 / np.sqrt(1 - b2)
        bp = np.sum(b * p[:, 1:], 1)
        E = g * (p[:, 0] - bp)
        coef = (g - 1) * bp / np.where(b2 > 0, b2, 1) - g * p[:, 0]
        return np.concatenate([E[:, None], p[:, 1:] + coef[:, None] * b], 1)
    p1c = boost(p1, beta)
    nhat = p1c[:, 1:] / np.linalg.norm(p1c[:, 1:], axis=1)[:, None]
    Ecm = np.sqrt(np.maximum(P[:, 0]**2 - np.sum(P[:, 1:]**2, 1), 0)) / 2
    # final polar angle wrt beam from channel angular distribution (rejection sampling)
    if channel == 'scalar':
        g = lambda c: (1 - c**2)**2
    elif channel == 'weyl':
        g = lambda c: c**2 * (1 - c**2) / 0.25
    else:
        g = lambda c: (1 - c**2)**2
    c = np.empty(n); todo = np.arange(n)
    while todo.size:
        cand = 2 * rng.random(todo.size) - 1
        acc = rng.random(todo.size) < g(cand)
        c[todo[acc]] = cand[acc]; todo = todo[~acc]
    phi = 2 * np.pi * rng.random(n)
    a = np.where(np.abs(nhat[:, :1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    e1 = np.cross(nhat, a); e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(nhat, e1)
    sn = np.sqrt(1 - c**2)
    dirv = c[:, None] * nhat + sn[:, None] * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2)
    k3 = np.concatenate([Ecm[:, None], Ecm[:, None] * dirv], 1)
    k4 = np.concatenate([Ecm[:, None], -Ecm[:, None] * dirv], 1)
    E3 = boost(k3, -beta)[:, 0]; E4 = boost(k4, -beta)[:, 0]
    return np.concatenate([E3, E4]), np.concatenate([E1 + E2, E1 + E2])


def main():
    res = {'checks': checks()}
    kap = {ch: kappa(ch) for ch in ('scalar', 'weyl', 'dirac', 'vector')}
    kap['scalar_conformal_initial'] = kappa('scalar', xi_i=1 / 6)
    kap['scalar_xi_f=1/6'] = kappa('scalar', xi_f=1 / 6)
    kap['scalar_xi_f=-0.0034'] = kappa('scalar', xi_f=-0.0034)
    res['kappa'] = kap
    res['kappa_analytic'] = {'scalar': 1 / (480 * np.pi), 'weyl': 1 / (1920 * np.pi), 'vector': 1 / (240 * np.pi),
                             'scalar_conformal_initial': 1 / (2880 * np.pi)}
    # central-charge ratios with initial identical-particle factor S: conformal scalar (S=1/2), Weyl (S=1), vector (S=1/2)
    base = 0.5 * kap['scalar_conformal_initial']
    res['central_charge_ratios'] = [1.0, kap['weyl'] / base, 0.5 * kap['vector'] / base]
    J = thermal_moments()
    res['J'] = J
    res['J_closed'] = {'B3': np.pi**4 / 15, 'F3': 7 * np.pi**4 / 120, 'B4': 24 * Z5, 'F4': 22.5 * Z5}
    # rates per hidden real dof (M_P=1,T=1): Gamma = N S kappa J3^2/(3 pi^4) ; rhodot = N S kappa 2 J3 J4/(3 pi^4)
    chans = {'higgs': (4, 0.5, kap['scalar'], 'B'), 'weyl45': (45, 1.0, kap['weyl'], 'F'),
             'vector12': (12, 0.5, kap['vector'], 'B')}
    G = {}; R = {}
    for k, (N, S, kp, st) in chans.items():
        G[k] = N * S * kp * J[st + '3']**2 / (3 * np.pi**4)
        R[k] = N * S * kp * 2 * J[st + '3'] * J[st + '4'] / (3 * np.pi**4)
    G['total'] = sum(G[k] for k in chans); R['total'] = sum(R[k] for k in chans)
    res['pair_rate_per_dof'] = G; res['energy_rate_per_dof'] = R
    res['beta_E_4dof'] = 4 * R['total']
    res['pair_rate_4dof'] = 4 * G['total']
    res['particle_number_rate_4dof'] = 8 * G['total']
    res['mean_pair_energy_over_T'] = R['total'] / G['total']
    # literature: Clery et al. 2112.15214 eq.(86) beta0 and Bernal et al. 1803.01866 eq.(31) alpha
    beta0 = 3997 * np.pi**3 / 41472000
    alphaB = 3997 * np.pi**3 / 663552000
    res['clery_beta0'] = beta0
    res['ours_over_clery'] = G['total'] / beta0
    res['ours_with_doubled_weyl_over_clery'] = (G['total'] + G['weyl45']) / beta0
    res['bernal_alpha'] = alphaB
    res['ours_over_bernal'] = G['total'] / alphaB
    # Monte Carlo of produced quanta: <E/T>, <T/E> weighted over channels by pair rate
    rng = np.random.default_rng(20261021)
    n = 600000
    Es, Ws = [], []
    stats = {}
    for k, (N, S, kp, st) in chans.items():
        ch = {'higgs': 'scalar', 'weyl45': 'weyl', 'vector12': 'vector'}[k]
        E, Ein = mc_final_energies(ch, st, n, rng)
        stats[k] = {'mean_E': float(E.mean()), 'mean_invE': float(np.mean(1 / E)),
                    'energy_conservation_check': float(np.mean(E[:n] + E[n:] - Ein[:n]))}
        Es.append(E); Ws.append(np.full(E.size, G[k] / E.size))
    E = np.concatenate(Es); W = np.concatenate(Ws)
    res['mc_per_channel'] = stats
    res['mc_mean_E_over_T'] = float(np.sum(W * E) / np.sum(W))
    res['mc_mean_T_over_E'] = float(np.sum(W / E) / np.sum(W))
    res['analytic_mean_E_over_T'] = R['total'] / (2 * G['total'])
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'vh2_rates.json'), 'w') as f:
        json.dump(res, f, indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
