"""rp_portals / freezein.py (Claude, own code).

Thermal (inverse-decay) production of the hidden vector A through hypercharge kinetic mixing in a hot radiation
era, its decay into X Xbar (branching B) or SM (1-B), the injected pair yield Y_pair and the conformal daughter
energy J = rho_X/(s T), and the raw-source normalisation R = F_R J of Codex e160356 (F_R = 4 (M_r/H_i)^2/E_X with
E_X = 38386.25 H_i^4, the a = 1 X-condensate energy of the spot point; both [Q]).

Two independent solvers:
  (1) moment solver, Maxwell-Boltzmann closure (<m/E> = K1/K2 at the bath temperature), written in T with
      LSODA/BDF (Codex: Radau in ln x);
  (2) momentum-resolved solver: the occupation f(q, x) at fixed q = p/T is a linear ODE per bin with the exact
      inverse-decay gain Gamma_SM (m/E) e^{-E/T} and loss Gamma_tot (m/E) f; no shape assumption on the A population.
      This tests the "MB time-profile closure" of the packet.
Analytic controls: Y_pair(inf) = C (3 pi/2) (Gamma_X/H_m)(1-B), J >= 5 Y_pair (prompt bound), the Bessel integrals,
and the raw bound R_tot >= R_0 B^{3/2} + R_sat (1-B).

Usage: python3 -B freezein.py <out.json>
"""
import json, math, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import kve, kv

MR = 2.435e18
HI = 1.3334961939833566e11
GS = 106.75
GA = 3
MX = 0.1*0.924089637375362
GY = 0.37
EX_COND = 38386.25063336385      # [Q] a=1 conformal X energy, H_i^4 units
R0 = 40.56262288525944           # [Q Codex budget] primordial band injection at e = 1e-7, eps_BX = 0
FR = 4*(MR/HI)**2/EX_COND
HSIG = 40*HI

def C_eq():
    return 45*GA/(4*math.pi**4*GS)

def point(e=1e-7, eps=1e-7):
    mA = 6*e*HSIG
    gx = e**2*mA/(48*math.pi)*(1-4*MX**2/mA**2)**1.5
    gsm = 41*eps**2*GY**2*mA/(96*math.pi)
    Hm = math.sqrt(math.pi**2*GS/90)*mA**2/MR
    return dict(mA=mA, gx=gx, gsm=gsm, gtot=gx+gsm, B=gx/(gx+gsm), Hm=Hm)

def k12(x):
    return kve(1, x)/kve(2, x)

def yeq(x):
    # C x^2 K2(x), with the scaled Bessel to avoid underflow
    return C_eq()*x*x*kve(2, x)*np.exp(-x) if x < 700 else 0.0

def moment_solver(e=1e-7, eps=1e-7, x0=1e-4, rtol=1e-9, method='LSODA'):
    """Variables (Y_A, Y_pair, J) in w = ln x (x = m_A/T), own formulation."""
    P = point(e, eps)
    gx, gsm, gtot, Hm = P['gx'], P['gsm'], P['gtot'], P['Hm']
    if gsm == 0:
        return dict(**P, Y_pair=0.0, J=0.0, R_thermal=0.0, Y_A_end=0.0)
    xend = max(100.0, 40*math.sqrt(Hm/gtot))
    scale = C_eq()*(3*math.pi/2)*gsm/Hm            # total A ever produced (analytic); used only to scale variables
    def rhs(w, y):
        x = math.exp(w); k = k12(x); h = Hm/(x*x)
        prod = gsm*k*yeq(x)/h
        return [prod/scale - gtot*k*y[0]/h, gx*k*y[0]/h, gx*x*y[0]/h]   # all three in units of `scale`
    sol = solve_ivp(rhs, (math.log(x0), math.log(xend)), [0., 0., 0.], method=method, rtol=rtol, atol=rtol*1e-4)
    assert sol.success, sol.message
    ya, yp, J = sol.y[:, -1]*scale
    return dict(**P, Y_pair=float(yp), J=float(J), R_thermal=float(J*FR), Y_A_end=float(ya), x_end=xend,
                Y_pair_analytic=scale*P['B'], prompt_J=5*scale*P['B'], nfev=sol.nfev)

def momentum_solver(e=1e-7, eps=1e-7, nq=240, qmin=1e-2, qmax=60.0, x0=1e-4, rtol=1e-8):
    """Occupation f_j(x) per comoving-momentum bin q_j = p/T (constant under redshift in RD, constant g*).
    E/T = sqrt(q^2 + x^2), m/E = x/sqrt(q^2+x^2).  Y_A = C int q^2 f dq (reduces to C x^2 K2 for f = e^{-E/T})."""
    P = point(e, eps)
    gx, gsm, gtot, Hm = P['gx'], P['gsm'], P['gtot'], P['Hm']
    q = np.geomspace(qmin, qmax, nq); lnq = np.log(q); wq = np.gradient(lnq)*q**3   # int q^2 f dq = sum f q^3 dlnq
    C = C_eq()
    xend = max(100.0, 40*math.sqrt(Hm/gtot))
    scale = C*(3*math.pi/2)*gsm/Hm
    def rhs(w, y):
        x = math.exp(w); h = Hm/(x*x)
        E = np.sqrt(q*q + x*x); mE = x/E
        f = y[:nq]
        df = (gsm*mE*np.exp(-E) - gtot*mE*f)/h
        YA = C*np.sum(wq*f); YAk = C*np.sum(wq*mE*f)
        return np.concatenate([df, [gx*YAk/h, gx*x*YA/h]])
    sol = solve_ivp(rhs, (math.log(x0), math.log(xend)), np.zeros(nq+2), method='LSODA', rtol=rtol, atol=rtol*1e-6)
    assert sol.success, sol.message
    yp, J = sol.y[-2, -1], sol.y[-1, -1]
    YAend = C*np.sum(wq*sol.y[:nq, -1])
    return dict(Y_pair=float(yp), J=float(J), R_thermal=float(J*FR), Y_A_end=float(YAend),
                Y_pair_analytic=scale*P['B'], prompt_J=5*scale*P['B'], nq=nq, nfev=sol.nfev)

def thermal_check():
    """Does C int q^2 e^{-sqrt(q^2+x^2)} dq reproduce C x^2 K2(x)?  (quadrature identity used by the momentum solver)"""
    out = {}
    for x in [0.1, 1.0, 3.0, 10.0]:
        val = quad(lambda q: q*q*math.exp(-math.sqrt(q*q+x*x)), 0, 60)[0]
        out[str(x)] = val/(x*x*kv(2, x)) - 1
    return out

def instant_reheating_fraction(xR):
    """Fraction of the hot-RD pair yield produced if the radiation era starts (instantly) at T_R = m_A/xR."""
    num = quad(lambda x: x**3*kve(1, x)*math.exp(-x), xR, max(xR+50, 60))[0]
    return num/(3*math.pi/2)

def main(out):
    res = dict(FR=FR, C=C_eq(), inputs=point())
    # Bessel controls
    res['bessel'] = dict(q3=quad(lambda x: x**3*kve(1, x)*math.exp(-x), 0, 700)[0]/(3*math.pi/2) - 1,
                         q4=quad(lambda x: x**4*kve(2, x)*math.exp(-x), 0, 700)[0]/(15*math.pi/2) - 1,
                         thermal_identity=thermal_check())
    # analytic bound pieces at the spot point
    P = point(); Ysat = C_eq()*(3*math.pi/2)*P['gx']/P['Hm']; Rsat = 5*Ysat*FR
    res['analytic'] = dict(Y_pair_saturated=Ysat, R_sat=Rsat, codex_R_sat=39183.97395590286,
                           Rsat_over_1p5R0=Rsat/(1.5*R0),
                           bound_min_over_B=min(R0*B**1.5 + Rsat*(1-B) for B in np.linspace(0, 1, 2001)))
    # epsilon grid, both solvers
    rows = []
    for eps in [1e-11, 1e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 1e-5, 1e-4, 1e-3]:
        m = moment_solver(eps=eps); p = momentum_solver(eps=eps)
        B = m['B']
        rows.append(dict(eps_BX=eps, B=B, Gamma_SM=m['gsm'],
                         moment=dict(Y_pair=m['Y_pair'], J=m['J'], R_thermal=m['R_thermal'], delay=m['J']/m['prompt_J'],
                                     number_err=m['Y_pair']/m['Y_pair_analytic']-1, Y_A_end_over_Ypair=m['Y_A_end']/m['Y_pair']),
                         momentum=dict(Y_pair=p['Y_pair'], J=p['J'], R_thermal=p['R_thermal'], delay=p['J']/p['prompt_J'],
                                       number_err=p['Y_pair']/p['Y_pair_analytic']-1),
                         momentum_over_moment_J=p['J']/m['J'],
                         R_primordial=R0*B**1.5, R_total_moment=m['R_thermal']+R0*B**1.5,
                         R_total_momentum=p['R_thermal']+R0*B**1.5,
                         raw_lower_bound=R0*B**1.5+Rsat*(1-B),
                         f_rest_unprocessed=2*MX*p['Y_pair']/4.4e-10))
        print('eps', eps, 'B', B, 'R_th moment', m['R_thermal'], 'momentum', p['R_thermal'], 'Ypair', p['Y_pair'], flush=True)
    res['eps_rows'] = rows
    # tolerance / grid control at eps = 1e-7
    res['controls'] = dict(moment_rtol1em11_BDF=moment_solver(rtol=1e-11, method='BDF')['R_thermal'],
                           moment_x0_1em5=moment_solver(x0=1e-5)['R_thermal'],
                           momentum_nq480=momentum_solver(nq=480)['R_thermal'],
                           momentum_qmax120=momentum_solver(qmax=120.0, nq=320)['R_thermal'],
                           momentum_rtol1em10=momentum_solver(rtol=1e-10)['R_thermal'])
    # vector-mass (e) dependence at fixed eps_BX = 1e-7 and at fixed eps_BX/e = 1
    erows = []
    for e in [1e-7, 3e-7, 1e-6, 3e-6, 1e-5]:
        for mode, eps in [('fixed_eps_1e-7', 1e-7), ('fixed_ratio_eps_eq_e', e)]:
            m = moment_solver(e=e, eps=eps)
            erows.append(dict(e=e, mode=mode, eps_BX=eps, mA=m['mA'], B=m['B'], Y_pair=m['Y_pair'], J=m['J'],
                              R_thermal=m['R_thermal'], f_rest=2*MX*m['Y_pair']/4.4e-10,
                              R_band_no_portal_scaled_1_over_e=R0*1e-7/e,
                              R_primordial=R0*1e-7/e*m['B']**1.5))
    res['e_rows'] = erows
    # instant-reheating cutoff (radiation era starting at T_R)
    res['instant_reheating_TR'] = {f'{xR:g}': instant_reheating_fraction(xR) for xR in [0.1, 0.3, 1, 2, 3, 5, 8, 10, 15, 20]}
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ('eps_rows', 'e_rows')}, indent=1))

if __name__ == '__main__':
    main(sys.argv[1])
