"""fx_observables item 1 (+ form factors for item 4): object family properties under GRUT-M1.  [C]

Own code (Claude, 1 Nov 2026).  Profiles from qball_m1_copy.py, an unchanged copy of
tracks/m1_vacuum/code/qball_m1.py (Numerov + Newton, fixed-Q mode).  For each Q:
  * Omega, E/(Q m_X), mass, r50/r90/r99 (charge radii), the phi-deficit half radius r_S (S = (1+S_c)/2);
  * scalar 'source' form factors in units where rho = m_X r:
        Fphi2(q) = 4 pi int rho^2 (S^2 - 1) j0(q rho) drho      (b |H|^2 phi^2 portal: source phi^2 - v^2)
        Fphi1(q) = 4 pi int rho^2 (S - 1)   j0(q rho) drho      (A phi |H|^2 trilinear portal: source phi - v)
        FX(q)    = 4 pi int rho^2 F^2       j0(q rho) drho      (lambda_HX |H|^2 |X|^2 portal: source |X|^2)
    and their first zeros q_0 (in m_X and GeV), recoil energies E_R = q_0^2 / 2 m_A for Xe131, Ge74, Ar40;
  * geometric object-object cross sections sigma = pi (2R)^2 for R = r50, r_S, r90, and sigma/M in cm^2/g;
  * merger energy release E(Q) + E(Q) - E(2Q).
Large Q (> 1e8) uses the bulk sphere R = (3Q / 4 pi n_B)^{1/3}, n_B = 449.31336788 m_X^3, and the fx_formation fit
E/Q = 0.775865 + 8.1429 Q^-1/3 + 125.79 Q^-2/3 - 2847/Q [Q own, fx_formation family.json fit3].
Writes results/objects.json and the profiles to $SP/wt2/fx_observables/profiles.npz (form factors are in the json).
"""
import json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qball_m1_copy as qb

RES = os.path.join(HERE, '..', 'results')
SP = os.environ.get('SP', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
SCR = os.path.join(SP, 'wt2', 'fx_observables')
os.makedirs(SCR, exist_ok=True)

MX = qb.MX_GEV                      # 0.0924089637 GeV
V = qb.V_GEV
HBARC_CM = 1.973269804e-14          # GeV cm
GEV_G = 1.78266192e-24              # g per GeV
CM_PER_INVMX = HBARC_CM / MX        # 2.1354e-13 cm
NB = 449.31336788
AMU = 0.9314941                     # GeV
TARGETS = {'Xe131': 131, 'Ge74': 74, 'Ar40': 40}
# nuclear masses with the mass excess neglected (A x amu), as in the usual recoil kinematics
MNUC = {k: A * AMU for k, A in TARGETS.items()}


def EQ_fit(Q):
    return 0.775865 + 8.1429 * Q**(-1 / 3) + 125.79 * Q**(-2 / 3) - 2847.0 / Q


def form_factors(r, F, S, h, qgrid):
    """4 pi int r^2 A(r) j0(q r) dr on the solver grid (trapezoid; integrands vanish at both ends)."""
    out = {}
    srcs = {'phi2': S * S - 1.0, 'phi1': S - 1.0, 'X2': F * F}
    for name, A in srcs.items():
        vals = []
        for q in qgrid:
            j0 = np.sinc(q * r / math.pi)  # sin(qr)/(qr)
            vals.append(4 * math.pi * h * np.sum(r * r * A * j0))
        out[name] = np.array(vals)
    return out


def first_zero(qgrid, Fq):
    s = np.sign(Fq)
    idx = np.where(s[1:] * s[:-1] < 0)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    # linear then secant refinement on the tabulated function
    q1, q2, f1, f2 = qgrid[i], qgrid[i + 1], Fq[i], Fq[i + 1]
    return float(q1 - f1 * (q2 - q1) / (f2 - f1))


def zero_refined(r, A, h, q_lo, q_hi):
    from scipy.optimize import brentq
    f = lambda q: 4 * math.pi * h * np.sum(r * r * A * np.sinc(q * r / math.pi))
    return brentq(f, q_lo, q_hi, xtol=1e-12)


def solve_Q(Q, L, seed, h=0.025):
    """Fixed-Q solve; same as qb.solve_fixed_Q but with the Newton tolerance raised above the round-off floor
    (~1e-13 max|u|, see the solver's own comment) for Q >= 5e7."""
    if seed is None or Q < 3e7:
        return qb.solve_fixed_Q(Q, L=L, h=h, seed=seed, Om0=0.88)
    sv = qb.Solver(L, h)
    s0, u0, w0, Om = seed
    u, w = sv.interp_from(s0, u0, w0)
    u, w, Om, res, it = sv.solve(u, w, Om, Qbar=qb.G2 * Q, tol=1e-12 * max(1.0, float(np.max(np.abs(u)))) * 10)
    ob = sv.observables(u, w, Om)
    ob.update(newton_residual=res, newton_iterations=it)
    return sv, u, w, Om, ob


def main():
    t0 = time.time()
    Qlist = [2e5, 1e5, 7.5e4, 4e5, 5e5, 1e6, 2e6, 1e7, 2e7, 3e7, 5e7, 7e7, 1e8]   # solve 2e5 first, continue down, then up
    Lfor = lambda Q: 60.0 if Q <= 1e6 else (90.0 if Q <= 2e7 else 140.0)
    qgrid = np.linspace(1e-4, 3.0, 3001)          # q in m_X
    rows = []
    prof = {}
    seed = None
    for Q in Qlist:
        L = Lfor(Q)
        if Q == 4e5:
            seed = seed_2e5
        sv, u, w, Om, ob = solve_Q(Q, L, seed)
        seed = (sv, u, w, Om)
        if Q == 2e5:
            seed_2e5 = seed
        r = sv.r
        F = u / r
        S = 1 + w / r
        Sc = ob['S0']
        # phi-deficit half radius
        target = 0.5 * (1 + Sc)
        rS = float(np.interp(target, S, r)) if np.all(np.diff(S[r < 3 * ob['r99']]) >= -1e-12) else float(r[np.argmin(np.abs(S - target))])
        ff = form_factors(r, F, S, sv.h, qgrid)
        row = dict(Q=Q, L=sv.L, h=sv.h, Omega=ob['Omega'], E_over_Qm=ob['E_over_Qm'], mass_GeV=ob['mass_GeV'],
                   S0=Sc, F0=ob['F0'], r50=ob['r50'], r90=ob['r90'], r99=ob['r99'], r_S_half=rS,
                   virial_residual=ob['virial_residual'], newton_residual=ob['newton_residual'])
        zeros = {}
        for name, A in (('phi2', S * S - 1), ('phi1', S - 1), ('X2', F * F)):
            qz = first_zero(qgrid, ff[name])
            if qz is not None:
                qz = zero_refined(r, A, sv.h, qz - 0.01, qz + 0.01)
            F0 = float(ff[name][0])
            zeros[name] = dict(F0_mX3units=F0, q0_mX=qz, q0_GeV=None if qz is None else qz * MX,
                               ER_keV={k: (None if qz is None else (qz * MX)**2 / (2 * m) * 1e6) for k, m in MNUC.items()})
        row['form_factor_zeros'] = zeros
        # charge check: Q = Omega int F^2 d^3rho / g^2
        row['charge_check'] = float(ob['Omega'] * ff['X2'][0] / qb.G2 / Q)
        # geometric cross sections
        M_g = ob['mass_GeV'] * GEV_G
        geo = {}
        for lab, R in (('r50', ob['r50']), ('r_S_half', rS), ('r90', ob['r90'])):
            sig_cm2 = math.pi * (2 * R * CM_PER_INVMX)**2
            geo[lab] = dict(R_mX=R, R_cm=R * CM_PER_INVMX, sigma_cm2=sig_cm2, sigma_over_M_cm2_g=sig_cm2 / M_g)
        row['geometric_obj_obj'] = geo
        row['bulk_sphere_R_mX'] = (3 * Q / (4 * math.pi * NB))**(1 / 3)
        rows.append(row)
        prof['Q%.3g' % Q] = np.vstack([r, F, S])
        prof['ff_Q%.3g' % Q] = np.vstack([qgrid, ff['phi2'], ff['phi1'], ff['X2']])
        print('Q=%.3g Om=%.6f E/Qm=%.6f M=%.1f GeV r50=%.2f r90=%.2f rS=%.2f  q0(phi2)=%.5f GeV ERXe=%.2f keV  sig/M(r90)=%.3e  chk=%.6f  t=%.0fs' % (
            Q, ob['Omega'], ob['E_over_Qm'], ob['mass_GeV'], ob['r50'], ob['r90'], rS,
            zeros['phi2']['q0_GeV'] or -1, zeros['phi2']['ER_keV']['Xe131'] or -1,
            geo['r90']['sigma_over_M_cm2_g'], row['charge_check'], time.time() - t0), flush=True)
    rows.sort(key=lambda z: z['Q'])
    np.savez_compressed(os.path.join(SCR, 'profiles.npz'), **prof)
    # mergers Q + Q -> 2Q from the solved family (where 2Q is solved) and from the fit
    byQ = {row['Q']: row for row in rows}
    mergers = []
    for Q in (1e5, 2e5, 5e5, 1e6, 1e7, 1e8):
        E1 = byQ[Q]['E_over_Qm'] * Q if Q in byQ else EQ_fit(Q) * Q
        E2 = byQ[2 * Q]['E_over_Qm'] * 2 * Q if (2 * Q) in byQ else EQ_fit(2 * Q) * 2 * Q
        dE = 2 * E1 - E2
        mergers.append(dict(Q=Q, dE_mX=dE, dE_GeV=dE * MX, frac_of_pair_mass=dE / (2 * E1),
                            source='solved' if (2 * Q) in byQ else 'fit', Omega_2Q=byQ[2 * Q]['Omega'] if (2 * Q) in byQ else None,
                            evap_cost_per_X_mX=(1 - byQ[2 * Q]['Omega']) if (2 * Q) in byQ else None))
    for Q in (1e12, 1e14, 1e17, 1e20):
        E1, E2 = EQ_fit(Q) * Q, EQ_fit(2 * Q) * 2 * Q
        mergers.append(dict(Q=Q, dE_mX=2 * E1 - E2, dE_GeV=(2 * E1 - E2) * MX, frac_of_pair_mass=(2 * E1 - E2) / (2 * E1), source='fit'))
    # large-Q objects (bulk spheres)
    large = []
    for Q in (1e8, 1e10, 1e12, 1e14, 1e17, 1e20):
        R = (3 * Q / (4 * math.pi * NB))**(1 / 3)
        Mg = EQ_fit(Q) * Q * MX * GEV_G
        sig = math.pi * (2 * R * CM_PER_INVMX)**2
        # phi^2 source of a uniform bulk sphere: (S_b^2 - 1) * volume; zero of 3 j1(x)/x at x = 4.4934
        q0 = 4.4934094579 / R
        large.append(dict(Q=Q, R_bulk_mX=R, R_cm=R * CM_PER_INVMX, mass_GeV=EQ_fit(Q) * Q * MX, mass_g=Mg,
                          sigma_cm2=sig, sigma_over_M_cm2_g=sig / Mg, q0_uniform_sphere_GeV=q0 * MX,
                          ER_Xe_keV=(q0 * MX)**2 / (2 * MNUC['Xe131']) * 1e6))
    out = dict(scope=__doc__.strip().splitlines()[0], units='m_X = %.10f GeV; 1/m_X = %.5e cm' % (MX, CM_PER_INVMX),
               rows=rows, mergers=mergers, large_Q_bulk_spheres=large,
               SIDM_reference='[Q/R] dwarfs/LSB 0.1-10 cm^2/g (Tulin & Yu 2018 arXiv:1705.02358; Kaplinghat, Tulin & Yu 2016 arXiv:1508.03339); clusters <~0.1-1 cm^2/g (Randall et al. 2008 arXiv:0704.0261 Bullet < 1.25; Harvey et al. 2015 arXiv:1503.07675 < 0.47; Sagunski et al. 2021 arXiv:2006.12515 < 0.35)',
               runtime_s=time.time() - t0)
    with open(os.path.join(RES, 'objects.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
    print('done %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
