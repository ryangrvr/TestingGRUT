"""Track pt_relax (Claude, own code): conformal time to binding, equilibrium end point, parametric rates.

Items 1, 2 and 4 of the pt_relax brief. Nothing from Codex is imported. Inputs quoted [Q]:
  * clock a(s) = a_m + C s, a_m = 8.869218738806, C = 1.9857962750634526 (Codex 5f5513b MATHEMATICS.md);
  * H_i = 1.33349619398e11 GeV, v = 0.924089637 GeV = 10 m_X, m_X = g v = 92.41 MeV, g^2 = 0.01, lambda_phi = 0.03,
    lambda_X = 1e-3 (kin_source kincommon.py, Codex 225028b);  a_ref = 7.047e12, T_vis = 45.7 GeV (kin_source);
  * n_Q = 8112.665, rho = 105481.25 (IR12_s1_w1, conformal, H_i units; RUN_IR12_s1_w1.json).
Conventions: conformal fields chi = a phi_phys / H_i (canonical real components), conformal momenta and frequencies
in H_i units at a = 1; physical energy = conformal energy * H_i / a.  The reduced Hamiltonian has no explicit s
dependence (fieldstep.cpp; a'' = 0 in radiation), so a conformal rate Gamma_c is constant in s and the number of
relaxation e-folds before the onset is Gamma_c * s_onset; Gamma_phys / H = Gamma_c * a / C.
"""
import json, math, sys, os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, fsolve
from scipy.special import zeta

AM, C = 8.869218738806, 1.9857962750634526
H_GEV = 133349619398.33566
V_PHI = 0.924089637375362
GX2, LP, LX = 0.01, 0.03, 1e-3
MX = math.sqrt(GX2) * V_PHI
A_REF = 7.047e12
NQ, RHO = 8112.665383318339, 105481.24746711778
KAPPA = LP * (V_PHI / H_GEV) ** 2            # lambda_phi v^2 / H_i^2: tachyonic phi mass^2 coefficient of a^2
OUT = os.path.join(os.path.dirname(__file__), '..', 'results')


def s_of_a(a):
    return (a - AM) / C


def tau_of_s(s):
    return 20 + AM * s + 0.5 * C * s * s


def a_onset(r2, phi_self=0.0):
    """phi conformal mass^2 = g^2 <r^2> + phi_self - KAPPA a^2 crosses zero."""
    return math.sqrt((GX2 * r2 + phi_self) / KAPPA)


# ------------------------------------------------------------------ Bose integrals (conformal units)
def bose(x):
    return 1.0 / np.expm1(x) if x < 700 else 0.0


def n_th(T, M, mu):
    """number density of one bosonic degree with energy sqrt(k^2+M^2) and chemical potential mu (|mu| < M or M=0, mu<=0)"""
    f = lambda k: k * k * bose((math.sqrt(k * k + M * M) - mu) / T) / (2 * math.pi ** 2)
    return quad(f, 0, 60 * T + 10 * M, limit=400)[0]


def e_th(T, M, mu):
    f = lambda k: k * k * math.sqrt(k * k + M * M) * bose((math.sqrt(k * k + M * M) - mu) / T) / (2 * math.pi ** 2)
    return quad(f, 0, 60 * T + 10 * M, limit=400)[0]


def fluct(T, M, mu):
    """<delta^2> thermal part for one real-field degree: int f/omega (vacuum dropped)"""
    f = lambda k: k * k * bose((math.sqrt(k * k + M * M) - mu) / T) / math.sqrt(k * k + M * M) / (2 * math.pi ** 2)
    return quad(f, 0, 60 * T + 10 * M, limit=400)[0]


# ------------------------------------------------------------------ equilibrium models
def model0(nQ, rho, dof):
    """free massless BE; all charge in a circular quartic condensate E_c = (3/4) lambda^(1/3) nQ^(4/3)."""
    Ec = 0.75 * LX ** (1 / 3) * nQ ** (4 / 3)
    T = (30 * (rho - Ec) / (math.pi ** 2 * dof)) ** 0.25
    nbar = zeta(3) / math.pi ** 2 * T ** 3
    return dict(model=f'free massless BE, {dof} dof thermal', T=T, condensate_energy_fraction=Ec / rho,
                mu=(LX * nQ) ** (1 / 3), Xbar_per_net=nbar / nQ, r_c2=nQ / (LX * nQ) ** (1 / 3))


def model1(nQ, rho, phi_mode='equilibrium', phi_frac=None):
    """Self-consistent Hartree (one-loop thermal masses, vacuum parts dropped) with a circular condensate.
    Unknowns T, r_c^2. X quasiparticles: average mass M_X^2 = 2 l r_c^2 + 2 l <|dpsi|^2> + g^2 <phi^2>,
    chemical potential mu = condensate frequency, mu^2 = l (r_c^2 + 2 <|dpsi|^2>) + g^2 <phi^2>  (=> M_X^2 - mu^2 = l r_c^2).
    phi: M_phi^2 = g^2 (r_c^2 + <|dpsi|^2>) + 3 lp <phi^2>; equilibrium (mu_phi = 0, same T) or a fixed energy share.
    Energy (Gaussian/Hartree): 1/2 mu^2 r_c^2 + l/4 <r^4> + g^2/2 <phi^2><r^2> + 3 lp/4 <phi^2>^2
                               + sum (omega f) - 1/2 sum M^2 <delta^2>."""
    def solve_phi(T, rc2, dpsi, Tphi=None):
        P = 0.0
        for _ in range(200):
            Mp2 = GX2 * (rc2 + dpsi) + 3 * LP * P
            Pn = fluct(Tphi if Tphi else T, math.sqrt(Mp2), 0.0)
            if abs(Pn - P) < 1e-12 * (1 + P): break
            P = 0.5 * P + 0.5 * Pn
        return P, Mp2

    def state(T, rc2, Tphi=None):
        dpsi, P = 0.0, 0.0
        for _ in range(300):
            MX2 = 2 * LX * rc2 + 2 * LX * dpsi + GX2 * P
            mu = math.sqrt(LX * (rc2 + 2 * dpsi) + GX2 * P)
            M = math.sqrt(MX2)
            dn = fluct(T, M, mu) + fluct(T, M, -mu)          # <|dpsi|^2> = <d1^2>+<d2^2> = int (f+ + f-)/omega
            if phi_mode == 'none':
                Pn, Mp2 = 0.0, None
            else:
                Pn, Mp2 = solve_phi(T, rc2, dn, Tphi)
            if abs(dn - dpsi) < 1e-11 * (1 + dn) and abs(Pn - P) < 1e-11 * (1 + Pn):
                dpsi, P = dn, Pn; break
            dpsi, P = 0.5 * dpsi + 0.5 * dn, 0.5 * P + 0.5 * Pn
        MX2 = 2 * LX * rc2 + 2 * LX * dpsi + GX2 * P
        mu = math.sqrt(LX * (rc2 + 2 * dpsi) + GX2 * P); M = math.sqrt(MX2)
        npl, nmi = n_th(T, M, mu), n_th(T, M, -mu)
        charge = mu * rc2 + npl - nmi
        d = dpsi / 2
        r4 = rc2 ** 2 + 8 * rc2 * d + 8 * d * d
        E = 0.5 * mu * mu * rc2 + LX / 4 * r4 + GX2 / 2 * P * (rc2 + dpsi) + 3 * LP / 4 * P * P
        E += e_th(T, M, mu) + e_th(T, M, -mu) - 0.5 * MX2 * dpsi
        Ephi = 0.0
        if phi_mode != 'none':
            Mp2 = GX2 * (rc2 + dpsi) + 3 * LP * P
            Ephi = e_th(Tphi if Tphi else T, math.sqrt(Mp2), 0.0) - 0.5 * Mp2 * P
            E += Ephi
        return dict(T=T, rc2=rc2, mu=mu, M_X=M, M_phi=(math.sqrt(GX2 * (rc2 + dpsi) + 3 * LP * P) if phi_mode != 'none' else None),
                    dpsi=dpsi, phi2=P, n_plus=npl, n_minus=nmi, charge=charge, energy=E, Ephi_frac=Ephi / E)

    def resid(x):
        T, rc2 = math.exp(x[0]), math.exp(x[1])
        st = state(T, rc2)
        return [st['charge'] / nQ - 1, st['energy'] / rho - 1]
    Ec = 0.75 * LX ** (1 / 3) * nQ ** (4 / 3)
    T0 = (30 * max(rho - Ec, 0.05 * rho) / (math.pi ** 2 * 2.5)) ** 0.25
    x0 = [math.log(T0), math.log(nQ / (LX * nQ) ** (1 / 3))]
    sol, info, ier, msg = fsolve(resid, x0, xtol=1e-12, full_output=True)
    if ier != 1 or max(abs(v) for v in resid(sol)) > 1e-8:
        raise RuntimeError(('model1 did not converge', nQ, rho, phi_mode, msg))
    st = state(math.exp(sol[0]), math.exp(sol[1]))
    st.update(model=f'Hartree BE, circular condensate, phi {phi_mode}', Xbar_per_net=st['n_minus'] / nQ,
              residual=resid(sol), r2_total=st['rc2'] + st['dpsi'],
              condensate_charge_fraction=st['mu'] * st['rc2'] / nQ)
    return st


def model1_phi_share(nQ, rho, share):
    """phi out of equilibrium with a fixed energy share (its own temperature); X sector in Hartree equilibrium."""
    if share == 0:
        return model1(nQ, rho, 'none')
    # crude: phi treated as a separate free gas of mass ~ g r_c with energy share*rho; X sector gets (1-share)*rho
    st = model1(nQ, rho * (1 - share), 'none')
    st['model'] = f'Hartree BE, X only, phi holding {share:.3g} of the energy (decoupled)'
    return st


def classical_lattice_endpoint(nQ, rho, N=32, L=12., fields=3, mH2=None, mu=None):
    """Rayleigh-Jeans end point on the FD4 lattice: T_cl per oscillator; A_- = (1/V) sum_k T_cl/(omega_k + mu)."""
    m = np.fft.fftfreq(N) * N; u = 2 * np.pi * m / N
    f = 4 * np.sin(u / 2) ** 2 * (1 + np.sin(u / 2) ** 2 / 3) / (L / N) ** 2
    kl = (f[:, None, None] + f[None, :, None] + f[None, None, :]).ravel()[1:]
    V = L ** 3
    Ec = 0.75 * LX ** (1 / 3) * nQ ** (4 / 3)
    if mu is None: mu = (LX * nQ) ** (1 / 3)
    rc2 = nQ / mu
    if mH2 is None: mH2 = 2 * LX * rc2
    Tcl = (rho - Ec) * V / (fields * (N ** 3 - 1))
    om = np.sqrt(kl + mH2)
    Am = np.sum(Tcl / (om + mu)) / V / nQ
    Ap = np.sum(Tcl / (om - mu)) / V / nQ
    return dict(N=N, L=L, fields_sharing=fields, T_classical=Tcl, mu=mu, mH2=mH2, Xbar_action_per_net=Am,
                thermal_charge_fraction=(Ap - Am), note='harmonic RJ with Hartree mass; interaction energy neglected [E]')


def main():
    res = {}
    # ---------------- item 1: conformal time to binding
    rows = []
    cases = [('brief |chi_X| = 60', 60.0 ** 2), ('brief |chi_X| = 110', 110.0 ** 2)]
    if len(sys.argv) > 1:
        meas = json.loads(open(sys.argv[1]).read())
        for k, v in meas.items(): cases.append((k, v))
    eq = model1(NQ, RHO, 'equilibrium')
    cases.append(('equilibrium r_c^2 + <|dpsi|^2>', eq['rc2'] + eq['dpsi']))
    for name, r2 in cases:
        a = a_onset(r2)
        s = s_of_a(a)
        rows.append(dict(case=name, r2=r2, a_onset=a, s_onset=s, tau_onset=tau_of_s(s), Gamma_crit=1 / s,
                         Gamma_crit_over_Teq=1 / s / eq['T'], s_onset_over_34=s / 34.,
                         m_X_conformal_at_onset=MX * a / H_GEV))
    sref = s_of_a(A_REF)
    rows.append(dict(case='kin_source a_ref (heavy spot point, case-L coordinate)', a_onset=A_REF, s_onset=sref,
                     tau_onset=tau_of_s(sref), Gamma_crit=1 / sref, Gamma_crit_over_Teq=1 / sref / eq['T'],
                     s_onset_over_34=sref / 34., m_X_conformal_at_onset=MX * A_REF / H_GEV))
    res['item1_onset'] = dict(kappa=KAPPA, rows=rows,
        statement='Reduced conformal Hamiltonian is s-independent (no masses, R = 0): a constant conformal rate gives '
                  'Gamma_c * s_onset e-folds; Gamma_phys/H = Gamma_c a / C. Survival of the budget needs '
                  'Gamma_c < 1/s_onset.')
    # ---------------- item 4: equilibrium end point
    eqs = [model0(NQ, RHO, 2), model0(NQ, RHO, 3), eq, model1(NQ, RHO, 'none')]
    for sh in [0.003, 0.05, 0.10]:
        eqs.append(model1_phi_share(NQ, RHO, sh))
    res['item4_equilibrium'] = eqs
    # thermal X-bar after phi condenses (m_X switched on, conformal), at fixed conformal T and nQ
    post = []
    for name, a in [('a_onset (equilibrium r^2)', a_onset(eq['rc2'] + eq['dpsi'])), ('a_ref', A_REF)]:
        mconf = MX * a / H_GEV
        for T in sorted({round(eq['T'], 6), round(eqs[0]['T'], 6), round(eqs[1]['T'], 6)}):
            for mu_over_m in [0.0, 0.8611951308217225, 0.99]:
                nm = n_th(T, mconf, -mu_over_m * mconf)
                post.append(dict(at=name, a=a, T_conf=T, m_conf=mconf, T_over_mX=T / mconf, T_phys_GeV=T * H_GEV / a,
                                 mu_over_m=mu_over_m, Xbar_per_net=nm / NQ))
    res['item4_after_onset'] = post
    a_eq = a_onset(eq['rc2'] + eq['dpsi'])
    V0conf = LP * V_PHI ** 4 / 4 * (a_eq / H_GEV) ** 4
    res['item4_phi_vacuum_energy'] = dict(a=a_eq, V0_conformal=V0conf, V0_over_rho=V0conf / RHO,
                                          note='lambda_phi v^4/4 in conformal units at the onset; same order as the hidden energy [E]')
    # classical lattice end point(s)
    res['classical_lattice_endpoint'] = [classical_lattice_endpoint(NQ, RHO, fields=3),
                                         classical_lattice_endpoint(NQ, RHO, fields=2),
                                         classical_lattice_endpoint(NQ, RHO, N=48, fields=3)]
    # ---------------- item 2: parametric rates (conformal units) [E]
    T = eq['T']; rc = math.sqrt(eq['rc2'])
    G = 1 / s_of_a(a_eq)
    nX = zeta(3) / math.pi ** 2 * T ** 3
    rates = []
    def add(name, vertex4, cbr, formula):
        lo, hi = cbr[0] * vertex4 * T, cbr[1] * vertex4 * T
        rates.append(dict(process=name, formula=formula, Gamma_low=lo, Gamma_high=hi,
                          efolds_to_onset_low=lo / G, efolds_to_onset_high=hi / G))
    # 2<->4 in a thermal gas: c * V^4 * T, c from 4-body phase space 1/(4 pi)^5 times combinatorics (10-1e3) [E]
    cb = (1 / (4 * math.pi) ** 5 * 10, 1 / (4 * math.pi) ** 5 * 1e3)
    add('2<->4 X only (|Phi|^4 vertex 4 lambda_X), thermal f~1', (4 * LX) ** 4, cb, 'c (4 lambda_X)^4 T')
    add('2<->4 with two g^2 vertices (phi phi Phi Phi* vertex 2 g^2), thermal', (2 * GX2) ** 4, cb, 'c (2 g^2)^4 T')
    add('2<->4 mixed g^2 x lambda_phi (X Xbar -> 4 phi)', (2 * GX2) ** 2 * (6 * LP) ** 2, cb, 'c (2g^2)^2 (6 lambda_phi)^2 T')
    add('2<->4 phi self (vertex 6 lambda_phi), thermal phi', (6 * LP) ** 4, cb, 'c (6 lambda_phi)^4 T')
    # 2<->2 X Xbar <-> phi phi: Gamma ~ n <sigma v>, sigma ~ |M|^2/(32 pi s), |M| = 2 g^2, <s> ~ 18 T^2 (massless BE)
    g22 = nX * (2 * GX2) ** 2 / (32 * math.pi * 18 * T * T)
    rates.append(dict(process='2<->2 X Xbar <-> phi phi (number-conserving; converts Xbar to phi)',
                      formula='n_X (2g^2)^2/(32 pi <s>)', Gamma_low=g22 / 3, Gamma_high=g22 * 3,
                      efolds_to_onset_low=g22 / 3 / G, efolds_to_onset_high=g22 * 3 / G))
    # condensate-assisted 1<->2 with cubic vertex v3 = 2 lambda_X r_c (X) or 2 g^2 r_c (phi phi X) at omega ~ 3T
    for nm, v3 in [('condensate-assisted 1<->2, X cubic vertex ~ 2 lambda_X r_c', 2 * LX * rc),
                   ('condensate-assisted 1<->2, phi phi X vertex ~ 2 g^2 r_c', 2 * GX2 * rc)]:
        g = v3 ** 2 / (16 * math.pi * 3 * T)
        rates.append(dict(process=nm, formula='v3^2/(16 pi omega), omega = 3T, x (0.1-1) kinematic factor',
                          Gamma_low=0.1 * g, Gamma_high=g, efolds_to_onset_low=0.1 * g / G, efolds_to_onset_high=g / G,
                          caveat='with a circular (single-sign) condensate these conserve total number including '
                                 'condensate quanta; they equilibrate X with the condensate, not the X-bar number by themselves'))
    res['item2_rates'] = dict(T_eq=T, r_c=rc, Gamma_crit=G, rates=rates,
                              classical_over_occupied_note='In the present medium f ~ 1e2-1e3 at k <= 8: classical rates carry '
                              'Bose enhancement f^(n-1); measured directly in the long run instead of estimated.')
    # ---------------- item 4b: heavy spot point (kin_source release N_X, E_X) with injected energy R*E_X [E]
    NXs, EXs = 2788.1767940228406, 38386.25063336385
    inj = []
    for R in [0.0, 0.80557, 1.14445, 4.0538, 40.56262, 48.49912]:
        for mode in ['equilibrium', 'none']:
            st = model1(NXs, EXs * (1 + R), mode)
            inj.append(dict(R=R, phi=mode, T_conf=st['T'], T_over_mX_at_a_ref=st['T'] / (MX * A_REF / H_GEV),
                            T_phys_GeV_at_a_ref=st['T'] * H_GEV / A_REF, Xbar_per_net=st['Xbar_per_net'],
                            condensate_charge_fraction=st['condensate_charge_fraction']))
            print('INJ R=%.3f phi %s T=%.2f T/m(a_ref)=%.2f Xbar/net=%.3f' % (R, mode, st['T'], inj[-1]['T_over_mX_at_a_ref'], st['Xbar_per_net']))
    res['item4_heavy_spot_injection'] = dict(N_X=NXs, E_X=EXs, rows=inj,
        note='massless-phase (pre-onset) Hartree equilibrium at the heavy spot release with injected conformal energy R E_X; '
             'R values from kin_source (Y debris 0.806, mediators 1.144, vectors e=1e-6 4.05, e=1e-7 band 40.56 / broad 48.50)')
    json.dump(res, open(os.path.join(OUT, 'scales.json'), 'w'), indent=1, default=float)
    # short printout
    for r in rows: print('ONSET', r['case'], 'a=%.3e s=%.3e Gcrit=%.2e Gcrit/T=%.2e mconf=%.2f' % (r['a_onset'], r['s_onset'], r['Gamma_crit'], r['Gamma_crit_over_Teq'], r['m_X_conformal_at_onset']))
    for e in eqs: print('EQ', e['model'], 'T=%.3f Xbar/net=%.4f' % (e['T'], e['Xbar_per_net']), {k: round(e[k], 4) for k in ['mu', 'rc2', 'M_X', 'Ephi_frac'] if k in e and e[k] is not None})
    for p in post: print('POST', p['at'], 'T/m=%.3f Tphys=%.3f GeV mu/m=%.3f Xbar/net=%.4f' % (p['T_over_mX'], p['T_phys_GeV'], p['mu_over_m'], p['Xbar_per_net']))
    print('V0', res['item4_phi_vacuum_energy'])
    for c in res['classical_lattice_endpoint']: print('RJ', c)
    for r in rates: print('RATE', r['process'], '%.2e-%.2e efolds %.1e-%.1e' % (r['Gamma_low'], r['Gamma_high'], r['efolds_to_onset_low'], r['efolds_to_onset_high']))


if __name__ == '__main__':
    main()
