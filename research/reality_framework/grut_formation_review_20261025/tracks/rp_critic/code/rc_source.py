#!/usr/bin/env python3
"""rp_critic: own momentum-resolved inverse-decay source for the kinetically mixed vector A (hot RD, MB, fixed mass).

Per comoving momentum q = p/T (constant in RD at constant g*), with x = m_A/T, E_q = sqrt(q^2+x^2), H = H_m/x^2:
   df_q/dln x = (x^2/H_m) (x/E_q) [ Gamma_SM e^{-E_q} - Gamma_tot f_q ]      (detailed balance against MB)
   Y_A = C int q^2 f_q dq,  C = 45 g_A/(4 pi^4 g*s), g_A = 3
   dY_pair/dln x = (Gamma_X/H_m) x^2 sum_q (x/E_q) f_q  (per-momentum decay rate)
   dJ/dln x      = (Gamma_X/H_m) x^3 Y_A                 (exact: rate x energy = Gamma m_A n_A)
Also the MB moment closure (Codex's) in the same code for the ratio. Written independently of rp_portals/freezein.py.
usage: rc_source.py OUT.json
"""
import sys, json, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import kn

M_R = 2.435e18
GSTAR = 106.75
G_Y = 0.37
E_X = 1e-7
H_SIG = 40 * 1.3334961939833566e11
M_A = 6 * E_X * H_SIG            # 3.2003908656e6 GeV
M_X = 0.1 * 0.924089637
GAMMA_X = E_X ** 2 * M_A / (48 * math.pi) * (1 - 4 * M_X ** 2 / M_A ** 2) ** 1.5
H_M = math.sqrt(math.pi ** 2 * GSTAR / 90) * M_A ** 2 / M_R
C_A = 45 * 3 / (4 * math.pi ** 4 * GSTAR)
EXN = 38386.250633                      # Codex's E_X normalisation (quoted)
H_I = 1.3334961939833566e11
F_R = 4 * (M_R / H_I) ** 2 / EXN

def gamma_sm(eps):
    return 41 * eps ** 2 * G_Y ** 2 * M_A / (96 * math.pi)

def momentum_resolved(eps, nq=400, qmin=1e-2, qmax=80.0, x0=1e-4, x1=None, rtol=1e-9, atol=1e-30):
    gsm = gamma_sm(eps); gtot = gsm + GAMMA_X
    if x1 is None:
        x1 = max(60.0, 40 * math.sqrt(H_M / gtot))   # well past the decay epoch x_d = sqrt(H_m/Gamma)
    lq = np.linspace(math.log(qmin), math.log(qmax), nq)
    q = np.exp(lq); w = np.gradient(lq) * q ** 3      # int q^2 f dq = sum q^3 f dln q
    def rhs(lnx, u):
        x = math.exp(lnx); f = u[:nq]
        E = np.sqrt(q * q + x * x)
        pref = x ** 2 / H_M * (x / E)
        df = pref * (gsm * np.exp(-E) - gtot * f)
        YA = C_A * np.dot(w, f)
        dYp = GAMMA_X / H_M * x ** 2 * C_A * np.dot(w, (x / E) * f)
        dJ = GAMMA_X / H_M * x ** 3 * YA
        return np.concatenate([df, [dYp, dJ]])
    u0 = np.zeros(nq + 2)
    sol = solve_ivp(rhs, (math.log(x0), math.log(x1)), u0, method='LSODA', rtol=rtol, atol=atol)
    assert sol.success, sol.message
    f_end = sol.y[:nq, -1]
    return dict(Y_pair=float(sol.y[nq, -1]), J=float(sol.y[nq + 1, -1]), remaining_A=float(C_A * np.dot(w, f_end)),
                x_end=x1, nfev=int(sol.nfev))

def closure(eps, x0=1e-4, x1=None, rtol=1e-10, atol=1e-30):
    gsm = gamma_sm(eps); gtot = gsm + GAMMA_X
    if x1 is None:
        x1 = max(60.0, 40 * math.sqrt(H_M / gtot))
    def rhs(lnx, u):
        x = math.exp(lnx); YA = u[0]
        k = kn(1, x) / kn(2, x) if x < 600 else 1 - 1.5 / x
        Yeq = C_A * x * x * kn(2, x) if x < 600 else C_A * math.sqrt(math.pi / 2) * x ** 1.5 * math.exp(-x) * (1 + 15 / (8 * x))
        return [x * x * k / H_M * (gsm * Yeq - gtot * YA), GAMMA_X / H_M * x * x * k * YA, GAMMA_X / H_M * x ** 3 * YA]
    sol = solve_ivp(rhs, (math.log(x0), math.log(x1)), [0.0, 0.0, 0.0], method='LSODA', rtol=rtol, atol=atol)
    return dict(Y_pair=float(sol.y[1, -1]), J=float(sol.y[2, -1]))

if __name__ == '__main__':
    out = dict(inputs=dict(m_A=M_A, Gamma_X=GAMMA_X, H_m=H_M, C=C_A, F_R=F_R), rows=[])
    for eps in [1e-9, 1e-8, 1e-7, 3e-7, 1e-6, 1e-5]:
        gsm = gamma_sm(eps); B = GAMMA_X / (GAMMA_X + gsm)
        ana = C_A * 1.5 * math.pi * GAMMA_X / H_M * (1 - B)
        mr = momentum_resolved(eps)
        mr2 = momentum_resolved(eps, nq=800, qmax=120.0, rtol=1e-10)
        cl = closure(eps)
        row = dict(eps=eps, B=B, Y_pair_analytic=ana, Y_pair_momentum=mr['Y_pair'], Y_pair_closure=cl['Y_pair'],
                   R_thermal_momentum=F_R * mr['J'], R_thermal_momentum_fine=F_R * mr2['J'], R_thermal_closure=F_R * cl['J'],
                   momentum_over_closure=mr['J'] / cl['J'], J_over_5Ypair_momentum=mr['J'] / (5 * mr['Y_pair']),
                   remaining_A=mr['remaining_A'], x_end=mr['x_end'], nfev=mr['nfev'])
        out['rows'].append(row)
        print(f"eps={eps:.0e} B={B:.6f} Ypair ana {ana:.6e} mom {mr['Y_pair']:.6e} clo {cl['Y_pair']:.6e} | "
              f"R_th mom {F_R*mr['J']:.1f} (fine {F_R*mr2['J']:.1f}) clo {F_R*cl['J']:.1f} ratio {mr['J']/cl['J']:.4f} | J/5Y {row['J_over_5Ypair_momentum']:.3f}")
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
