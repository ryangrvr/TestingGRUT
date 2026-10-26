"""rp_portals / widths.py (Claude, own code).

Portal widths, branching, one-loop mixed charge trace and the isolated beta_b, for the heavy spot point
(e_X = 1e-7, h_Sigma = 40 H_i, H_i = 1.3334962e11 GeV, m_A = 6 e_X h_Sigma).

Independent derivation: the SM hypercharge content is enumerated field by field; the vector widths use the
textbook massless-final-state formulas Gamma(V -> psi psi) = g^2 m_V/(24 pi) per Weyl fermion and
Gamma(V -> phi phi*) = g^2 m_V/(48 pi) per complex scalar (times beta^3 for massive scalars); beta_b from the
Coleman-Weinberg log derivative with sympy.  Codex's numbers are quoted [Q] only for the comparison columns.

Usage: python3 -B widths.py <out.json>
"""
import json, math, sys
import sympy as sp

HI = 1.3334961939833566e11          # GeV, heavy spot point [Q own tracks / Codex INPUTS]
EX = 1e-7                            # hidden gauge coupling
HSIG = 40.0 * HI                     # h_Sigma (Sigma = h/sqrt2 e^{i theta}), hidden charge 6
MX = 0.1 * 0.924089637375362         # m_X = g_X v
GY_MZ = 0.3575                       # g'(M_Z), PDG-level value [R]
MZ = 91.1876

# --- SM field content above electroweak breaking: (name, hypercharge Y, multiplicity of Weyl fermions or complex scalars)
WEYL = []
for gen in range(3):
    WEYL += [("Q_L", 1/6, 6), ("u_R", 2/3, 3), ("d_R", -1/3, 3), ("L_L", -1/2, 2), ("e_R", -1.0, 1)]
SCALARS = [("H", 1/2, 2)]            # one complex doublet = two complex scalars
HIDDEN = [("X", 1, 0.0), ("Y", 3, 0.0), ("Sigma", 6, 0.0)]   # (name, hidden charge, hypercharge)

def sumY2(table):
    return sum(n * Y * Y for _, Y, n in table)

def gY_run(mu):
    # one-loop SM hypercharge running, b = 41/6 (no normalisation factor 3/5), from M_Z
    inv = 1/GY_MZ**2 - (41/6)/(8*math.pi**2) * math.log(mu/MZ)
    return 1/math.sqrt(inv)

def widths(eps, gY, mA):
    gsm = eps**2 * gY**2 * mA * (sumY2(WEYL)/(24*math.pi) + sumY2(SCALARS)/(48*math.pi))
    beta3 = (1 - 4*MX**2/mA**2)**1.5
    gx = EX**2 * mA/(48*math.pi) * beta3
    return gsm, gx, gx/(gx+gsm)

def main(out):
    mA = 6 * EX * HSIG
    sY2_f, sY2_s = sumY2(WEYL), sumY2(SCALARS)
    coeff = sY2_f/(24*math.pi) + sY2_s/(48*math.pi)         # Gamma_SM = eps^2 gY^2 mA * coeff
    # Codex's closed form 41/(96 pi)
    codex_coeff = 41/(96*math.pi)
    gY = 0.37
    gsm, gx, B = widths(1e-7, gY, mA)
    # sensitivity: running g_Y at mu = m_A
    gYrun = gY_run(mA)
    gsm_r, _, B_r = widths(1e-7, gYrun, mA)
    # mixed one-loop charge trace sum_i (hidden q_i)(hypercharge Y_i) over every field
    trace = sum(q*Y for _, q, Y in HIDDEN) + sum(0.0*Y for _, Y, n in WEYL + SCALARS)
    # isolated beta_b from the Coleman-Weinberg potential: V1 = (1/64pi^2) Str M^4 (ln M^2/mu^2 - 3/2)
    phi2, H2, m0, g, lam, mu = sp.symbols('phi2 H2 m0 g lam mu', positive=True)
    MX2 = m0**2 + g**2*phi2 + lam*H2            # field-dependent X mass^2, two real components
    V1 = sp.Rational(1, 64)/sp.pi**2 * 2 * MX2**2 * (sp.log(MX2/mu**2) - sp.Rational(3, 2))
    dV1_dlnmu = sp.simplify(mu*sp.diff(V1, mu))   # = -(1/32pi^2) Tr M^4
    # mu-independence of V_tree + V1: d b/d ln mu * phi2*H2 + ... + dV1/dlnmu = 0 -> beta_b = coefficient of phi2*H2 in -dV1/dlnmu
    poly = sp.Poly(sp.expand(-dV1_dlnmu), phi2, H2)
    beta_b_over = sp.simplify(poly.coeff_monomial(phi2*H2))
    beta_b_expected = g**2*lam/(8*sp.pi**2)
    res = dict(
        mA_GeV=mA, m_X_GeV=MX,
        sumY2_weyl=sY2_f, sumY2_complex_scalars=sY2_s,
        Gamma_SM_coefficient=coeff, codex_41_over_96pi=codex_coeff, coefficient_ratio=coeff/codex_coeff,
        gY=gY, Gamma_SM_GeV_eps1em7=gsm, Gamma_X_GeV=gx, B_X_eps1em7=B,
        codex_Gamma_SM=5.956199301217943e-10, codex_Gamma_X=2.122325108666791e-10, codex_B=0.26271197572541255,
        gY_run_at_mA=gYrun, Gamma_SM_with_running_gY=gsm_r, B_with_running_gY=B_r,
        Gamma_SM_over_Gamma_X_per_eps2_over_e2=(coeff*48*math.pi)*gY**2,   # = 20.5 gY^2
        one_loop_mixed_charge_trace=trace,
        beta_b_isolated=str(beta_b_over), beta_b_expected=str(beta_b_expected),
        beta_b_check=bool(sp.simplify(beta_b_over - beta_b_expected) == 0),
        beta_b_over_lambda_HX_at_gX0p1=float(0.01/(8*math.pi**2)),
        codex_beta_b_over_lambda_HX=0.00012665147955292222,
        B_grid={str(e): widths(e, gY, mA)[2] for e in [1e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 1e-5]},
        note=("Hidden fields carry no hypercharge and SM fields no hidden charge, so the one-loop B-A_X mixing trace vanishes. "
              "Stronger statement [A]: with X, Y, Sigma and the two couplings kappa Y*X^3, mu Sigma*Y^2, the two coupling phases are "
              "removable by rephasing (3 fields, 1 gauge direction), so hidden charge conjugation (X,Y,Sigma -> conjugates, A_X -> -A_X) "
              "is an exact symmetry of the stated tree Lagrangian including lambda_HX|X|^2|H|^2; B_mn F_X^mn is odd under it, so "
              "epsilon_BX is not generated at any loop order from these couplings alone. A CP-violating holder term (the supplied "
              "phase source) breaks that symmetry and could feed epsilon_BX radiatively only together with a hidden-SM bridge; "
              "not computable without the holder."),
    )
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps(res, indent=1))

if __name__ == '__main__':
    main(sys.argv[1])
