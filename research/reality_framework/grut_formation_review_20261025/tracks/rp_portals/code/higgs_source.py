"""rp_portals / higgs_source.py (Claude, own code).

Scalar-portal source b_Hphi phi^2 H^dag H at the light-sector point (g_X = 0.1, lambda_phi = 0.03, v = 0.924 GeV):
Gamma(h -> ss), freeze-in of h -> ss -> 2 (X Xbar) from the thermal Higgs (vacuum masses), the injected pair yield,
its unprocessed rest-mass fraction against the 4.4e-10 GeV per-entropy benchmark, the prompt conformal energy
J_h = 5 Y_hss (h decays promptly: Gamma_h >> H) and the raw-source R_h = F_R J_h in the spot-point normalisation.
Independent of Codex's code: own freeze-in quadrature in T (not the closed form), a g*(T) sensitivity with
Boltzmann-suppressed heavy SM species, the ATLAS one-channel bounds, mixing angle, s -> X Xbar openness.

Usage: python3 -B higgs_source.py <out.json>
"""
import json, math, sys
import numpy as np
from scipy.integrate import quad
from scipy.special import kv, kve
from freezein import FR, MR, GS, MX

MH, VH, GAMH = 125.0, 246.0, 4.1e-3
GX, LAMPHI, V = 0.1, 0.03, 0.924089637375362
MS = math.sqrt(2*LAMPHI)*V
DM = 4.4e-10

def gstar_T(T):
    """SM g*_rho(T) with massless species 86.25 plus Boltzmann-suppressed W, Z, h, t (MB energy density ratio)."""
    g = 86.25
    for m, gi, ferm in [(80.4, 6, False), (91.2, 3, False), (125., 1, False), (173., 12, True)]:
        x = m/T
        # rho_massive/rho_massless (MB): (x^2/2)[x K1(x) + 3 K2(x)] / 3 ... use rho = g m^4/(2pi^2)[3K2/x^2 + K1/x]*... exact MB:
        # rho_MB(m,T) = g T^4/(2 pi^2) x^2 [3 K2(x) + x K1(x)]  ;  massless MB: g T^4 (3/pi^2)
        ratio = (x*x*(3*kve(2, x) + x*kve(1, x))*math.exp(-x)/2)/3
        g += gi*ratio*(7/8 if ferm else 1)
    return g

def width_hss(b):
    return b*b*VH*VH/(8*math.pi*MH)*math.sqrt(1-4*MS*MS/MH**2)

def freezein_Yh(b, gstar_const=None):
    """Y of h -> ss decays per entropy: int dT Gamma n_h^eq / (s H T) with n_h^eq = m^2 T K1... MB: Gamma <m/E> n_eq = Gamma m^2 T K1(x)/(2pi^2)."""
    G = width_hss(b)
    def integrand(lnT):
        T = math.exp(lnT); x = MH/T
        g = gstar_const if gstar_const else gstar_T(T)
        s = 2*math.pi**2/45*g*T**3
        H = math.sqrt(math.pi**2*g/90)*T*T/MR
        rate = G*MH*MH*T*kve(1, x)*math.exp(-x)/(2*math.pi**2)
        return rate/(s*H)            # dY/dlnT (T decreasing: dt = -dlnT/H)
    val = quad(integrand, math.log(MH/60), math.log(MH*50), limit=400)[0]
    return val

def energy_Jh(b, gstar_const=None):
    """Prompt conformal energy of the h -> ss -> X Xbar daughters: int dT Gamma m n_h^eq /(s T H T) (energy identity)."""
    G = width_hss(b)
    def integrand(lnT):
        T = math.exp(lnT); x = MH/T
        g = gstar_const if gstar_const else gstar_T(T)
        s = 2*math.pi**2/45*g*T**3
        H = math.sqrt(math.pi**2*g/90)*T*T/MR
        neq = MH*MH*T*kve(2, x)*math.exp(-x)/(2*math.pi**2)
        return G*MH*neq/(s*T*H)
    return quad(integrand, math.log(MH/60), math.log(MH*50), limit=400)[0]

def main(out):
    b = 5e-11
    G = width_hss(b)
    Hh = math.sqrt(math.pi**2*GS/90)*MH*MH/MR
    closed = 2*135/(8*math.pi**3*GS)*G/Hh                     # pairs: two s per h, one pair per s
    Yh_const = freezein_Yh(b, GS); Yh_run = freezein_Yh(b)
    Jh_const = energy_Jh(b, GS); Jh_run = energy_Jh(b)
    # g*s^{1/3} correction to R when the injection happens at g*(T) < g*(T_R)
    gT = gstar_T(MH/3)
    theta = 2*b*V*VH/(MH**2 - MS**2)
    Gs_XX = GX**4*V*V/(4*math.pi*MS)*math.sqrt(1-4*MX*MX/MS**2)
    Gs_SM_ee = theta**2*(0.000511**2*MS/(8*math.pi*VH**2))*(1-4*0.000511**2/MS**2)**1.5
    extra = GAMH*0.107/(1-0.107)
    bmax = math.sqrt(extra/(width_hss(1.0)))
    lmax = math.sqrt(extra/(VH*VH/(16*math.pi*MH)*math.sqrt(1-4*MX*MX/MH**2)))
    res = dict(
        m_s_GeV=MS, s_to_XX_open=MS > 2*MX, beta_s_XX=math.sqrt(1-4*MX*MX/MS**2),
        Gamma_hss_GeV=G, codex_Gamma_hss=4.8156786852335187e-20,
        H_at_mh_GeV=Hh,
        Y_pair_closed_form_gstar106p75=closed, codex_Y_pair=2.236557958585659e-08,
        Y_pair_quadrature_gstar106p75=2*Yh_const, Y_pair_quadrature_gstar_T=2*Yh_run,
        gstar_at_mh_over_3=gT,
        f_rest_unprocessed_gstar106p75=2*MX*2*Yh_const/DM, f_rest_unprocessed_gstar_T=2*MX*2*Yh_run/DM,
        codex_f_rest=9.394454695083644,
        J_h_prompt_gstar106p75=Jh_const, J_h_over_5Yhss=Jh_const/(5*Yh_const),
        R_h_gstar106p75=Jh_const*FR, R_h_gstar_T=Jh_run*FR*(GS/gT)**(1/3),
        R_h_scaling='R_h, Y_pair, f_rest all scale as b^2',
        mixing_angle=theta, codex_mixing_per_b=0.02909782991741602, mixing_angle_per_b=theta/b,
        Gamma_s_to_XX_GeV=Gs_XX, Gamma_s_to_ee_GeV=Gs_SM_ee, tau_s_if_only_ee_s=6.582e-25/Gs_SM_ee if Gs_SM_ee > 0 else None,
        Gamma_h_over_H_at_mh=GAMH/Hh,
        ATLAS_one_channel=dict(b_max=bmax, lambda_HX_max=lmax, codex_b_max=0.005050092595086073, codex_lambda_max=0.007141889923011246),
        b_for_f_rest_equal_0p01=b*math.sqrt(0.01/(2*MX*2*Yh_const/DM)),
        note=("Y_pair counts X Xbar pairs per visible entropy from h -> ss with each s -> X Xbar (open, beta = 0.58). "
              "Vacuum masses, MB statistics, thermal h, no erasure; the freeze-in integrand x^3 K1 peaks at T ~ m_h/2.5, in the broken phase."))
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps(res, indent=1))

if __name__ == '__main__':
    main(sys.argv[1])
