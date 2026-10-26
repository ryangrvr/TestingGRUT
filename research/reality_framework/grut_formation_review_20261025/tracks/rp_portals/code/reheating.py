"""rp_portals / reheating.py (Claude, own code).

Dependence of the thermal vector source on the reheating history.  A perturbatively decaying inflaton (matter,
width Gamma_I) feeds radiation; the hidden vector A is produced by SM inverse decays with the MB closure of
freezein.py on the actual T(a), decays with Gamma_tot, and the pair yield / conformal daughter energy are compared
with the hot radiation-era values (T_R >> m_A).  T_R := (90/(pi^2 g*))^{1/4} sqrt(Gamma_I M_r).  Constant g*,
no back-reaction of A on the bath, no entropy from A.  Also the two heavy-window histories: prompt radiation at
H_i (T_R = 3.08e14 GeV) and the finite releases Gamma_phi/H_i = 0.1..10 of hw_verify (T_R ~ 1e14 GeV) [Q own].

Usage: python3 -B reheating.py <out.json>
"""
import json, math, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import kve
from freezein import point, C_eq, FR, GS, MR, HI, MX, moment_solver

def TR_of_Gamma(G):
    return (90/(math.pi**2*GS))**0.25*math.sqrt(G*MR)

def Gamma_of_TR(TR):
    return TR**2/(MR*math.sqrt(90/(math.pi**2*GS)))

def run(TR, e=1e-7, eps=1e-7, Hstart_over_Gamma=1e4, rtol=1e-8):
    P = point(e, eps); mA, gx, gsm, gtot = P['mA'], P['gx'], P['gsm'], P['gtot']
    GI = Gamma_of_TR(TR)
    H0 = Hstart_over_Gamma*GI
    rhoI0 = 3*MR**2*H0**2
    C = C_eq()
    # state: ln rho_I, rho_r/rhoI0, N_A (comoving n_A a^3 / rhoI0^{3/4}), N_pair, E_pair (comoving energy / rhoI0)
    # units: energies in GeV; a = 1 at start.  Work with u = ln a.
    def T_of(rr):
        return (30*rr*rhoI0/(math.pi**2*GS))**0.25 if rr > 0 else 0.0
    def rhs(u, y):
        lrI, rr, NA, NP, EP = y
        rI = math.exp(lrI)*rhoI0
        H = math.sqrt((rI + rr*rhoI0)/(3*MR**2))
        a = math.exp(u); a3 = a**3
        T = T_of(rr)
        if T > 0:
            x = mA/T; k = kve(1, x)/kve(2, x)
            neq = 3*mA*mA*T/(2*math.pi**2)*kve(2, x)*math.exp(-x) if x < 700 else 0.0   # g_A m^2 T K2/(2 pi^2)
        else:
            x, k, neq = 1e9, 1.0, 0.0
        nA = NA/a3*rhoI0**0.75
        dNA = a3*(gsm*k*neq - gtot*k*nA)/H/rhoI0**0.75
        dNP = a3*gx*k*nA/H/rhoI0**0.75
        dEP = a**4*gx*mA*nA/H/rhoI0
        return [-3 - GI/H, -4*rr + GI*rI/H/rhoI0, dNA, dNP, dEP]
    # integrate until T has fallen to m_A/x_end with x_end as in the hot-RD solver and all A decayed
    xend = max(100.0, 40*math.sqrt(P['Hm']/gtot))
    # end when a is large enough: T_R a_R ~ T a; estimate a_end from entropy: T(a) ~ T_R (a_R/a); a_R ~ (H0/GI)^{2/3}
    aR = (H0/GI)**(2/3)
    a_end = aR*TR/(mA/xend)*3
    sol = solve_ivp(rhs, (0, math.log(a_end)), [0.0, 0.0, 0.0, 0.0, 0.0], method='LSODA', rtol=rtol, atol=1e-30,
                    dense_output=False)
    assert sol.success, sol.message
    lrI, rr, NA, NP, EP = sol.y[:, -1]
    a = math.exp(sol.t[-1]); T = T_of(rr)
    s = 2*math.pi**2/45*GS*T**3
    Ypair = NP*rhoI0**0.75/(a**3*s)
    J = EP*rhoI0/(a**4*s*T)
    # Tmax diagnostic
    Tmax = max(T_of(v) for v in sol.y[1])
    return dict(TR=TR, Gamma_I=GI, Tmax=Tmax, Y_pair=Ypair, J=J, R_thermal=J*FR, remaining_A_over_pairs=NA/NP if NP > 0 else None,
                xR=mA/TR, rhoI_end_over_rhor=math.exp(lrI)/rr)

def main(out):
    hot = moment_solver()
    rows = []
    for TR in [1e4, 3e4, 1e5, 2e5, 4e5, 7e5, 1e6, 2e6, 4e6, 1e7, 3e7, 1e8, 1e9]:
        r = run(TR)
        r['Y_pair_over_hotRD'] = r['Y_pair']/hot['Y_pair']; r['J_over_hotRD'] = r['J']/hot['J']
        rows.append(r); print('TR', TR, 'Y ratio', r['Y_pair_over_hotRD'], 'J ratio', r['J_over_hotRD'], flush=True)
    # power-law fit of the suppressed branch
    lo = [r for r in rows if r['TR'] <= 4e5]
    slope = np.polyfit(np.log([r['TR'] for r in lo]), np.log([r['Y_pair_over_hotRD'] for r in lo]), 1)[0]
    # heavy-window histories
    TR_prompt = (90/(math.pi**2*GS))**0.25*math.sqrt(3*MR**2*HI**2/(3*MR**2))*0+ (30*3*MR**2*HI**2/(math.pi**2*GS))**0.25
    heavy = dict(prompt_radiation_T_R=TR_prompt, prompt_xR=point()['mA']/TR_prompt,
                 finite_release={f'Gamma_phi/H_i={g}': TR_of_Gamma(g*HI) for g in [0.1, 1.0, 10.0]},
                 codex_reheating_runs={f'Gamma_I/H_i={g}': TR_of_Gamma(g*HI) for g in [0.03, 0.1, 0.3]})
    # the raw bound's minimum as a function of T_R: R_total >= R0 B^{3/2} + R_sat(T_R)(1-B)
    R0 = 40.56262288525944
    Rsat_hot = 5*hot['Y_pair_analytic']/(1-hot['B'])*FR      # 5 C (3pi/2) (Gamma_X/H_m) F_R
    bound = []
    for r in rows:
        Rs = Rsat_hot*r['J_over_hotRD']
        B = np.linspace(0, 1, 4001); f = R0*B**1.5 + Rs*(1-B)
        bound.append(dict(TR=r['TR'], R_sat=Rs, min_bound=float(f.min()), B_at_min=float(B[f.argmin()]),
                          outlet_can_lower_injection=bool(f.min() < R0*0.999)))
    res = dict(hot_RD=dict(Y_pair=hot['Y_pair'], J=hot['J'], R_thermal=hot['R_thermal']), rows=rows,
               low_TR_slope_dlnY_dlnTR=float(slope), heavy_window=heavy, bound_vs_TR=bound)
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != 'rows'}, indent=1))

if __name__ == '__main__':
    main(sys.argv[1])
