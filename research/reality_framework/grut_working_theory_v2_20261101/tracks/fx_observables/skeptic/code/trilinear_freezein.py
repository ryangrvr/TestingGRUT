"""fx_observables SKEPTIC (Claude, 1 Nov 2026): 2->2 estimate of phi freeze-in through h-phi mixing for the trilinear
portal A phi |H|^2.  [E] own code.

Reduced cross sections (massless quarks, scalar Yukawa y = theta m_t/v_H, colour/spin summed, thermal-mass cutoff
m_th^2 = g_s^2 T^2/6, L = ln(<s>/m_th^2), <s> ~ 18 T^2) derived here:
   q g -> q phi and qbar g -> qbar phi:  sigma_hat = 2 x 8 alpha_s y^2 (L - 3/2)      [Sum|M|^2 = 16 g^2 y^2 t^2/(s(-u))]
   q qbar -> g phi:                      sigma_hat = 16 alpha_s y^2 L                [Sum|M|^2 = 16 g^2 y^2 s^2/(t u)]
 gamma = sigma_hat T^4/(16 pi^4) (MB), x top Boltzmann factor (m/T)^2 K2(m/T)/2;  Y = int gamma/(s H T) dT,
 Delta S/Q = (3/r) Y/Y_Q with r = T_h/T_vis = 0.0059 (track).  Two thermal histories for the mixing:
   'fixed'     : theta = A v_H/m_h^2 (zero-T), m_t fixed, cut at T_EW = 160 GeV (the track's assumption)
   'meanfield' : v(T) = v sqrt(1 - T^2/T_c^2) floored at 0.3 v; m_t(T), m_h(T) ~ v(T)  -> y_eff = A m_t/(m_h^2 v(T)/v)
"""
import json, math, sys
import numpy as np
from scipy.special import kn

VH, MH, MT, MPL, YQ, R = 246.22, 125.20, 172.7, 1.22089e19, 5.07941423e-9, 0.0059
A_X = 2.737e-10       # GeV, XLZD 10 events (track)
TC = 160.0


def run(mode, A=A_X, gs=100.0, als=0.11):
    T = np.geomspace(5.0, TC, 4000)
    if mode == 'fixed':
        y = A * VH / MH**2 * MT / VH * np.ones_like(T)
        mt = MT * np.ones_like(T)
    else:
        vr = np.maximum(np.sqrt(np.clip(1 - (T / TC)**2, 0, 1)), 0.3)
        y = A * MT / (MH**2 * vr)
        mt = MT * vr
    g2 = 4 * math.pi * als
    L = np.log(18.0 / (g2 / 6))
    sh = 16 * als * y * y * (L - 1.5) + 16 * als * y * y * L
    x = mt / T
    supp = 0.5 * x * x * kn(2, x)
    gam = sh * T**4 / (16 * math.pi**4) * supp
    s = 2 * math.pi**2 / 45 * gs * T**3
    H = 1.66 * math.sqrt(gs) * T * T / MPL
    dY = gam / (s * H * T)
    Y = float(np.trapezoid(dY, T))
    Tpk = float(T[np.argmax(dY * T)])
    c_eff = None
    return dict(Y_phi=Y, dS_over_Q=3 / R * Y / YQ, T_peak_GeV=Tpk, L=float(L))


def main(fn):
    out = dict(A_GeV=A_X, fixed=run('fixed'), meanfield=run('meanfield'),
               track=dict(Y_top_c01=3.87e-12, Y_top_c1=3.87e-11, dS_over_Q_range=[0.4, 4.0]))
    for k in ('fixed', 'meanfield'):
        out[k]['equivalent_track_c'] = out[k]['Y_phi'] / 3.87e-11
        budget = (1.4, 2.2)
        out[k]['A_max_survival_GeV'] = [A_X * math.sqrt(b / out[k]['dS_over_Q']) for b in budget]
        out[k]['XLZD_events_at_A_max'] = [10 * (a / A_X)**2 for a in out[k]['A_max_survival_GeV']]
    print(json.dumps(out, indent=1))
    json.dump(out, open(fn, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
