#!/usr/bin/env python3
"""rp_critic: independent homogeneous Jordan-frame reheating solver (own code).

Equations re-derived from the action S = int sqrt(-g)[F R/2 - (dI)^2/2 - sum|dz|^2 - m_I^2 I^2/2 - V_h + L_r],
F = M_r^2 - 2 xi |X|^2, written here from scratch in scipy (DOP853, H evolved, Friedmann checked).
Nothing is reused from Codex's reheat.c or from rp_reheating's rpreheat.c/rp_scipy.py except the quoted
Y free-mode diagnostic convention (hw_verify / Codex), which is implemented again here from its formula.

usage: rc_reheat.py OUT.json b gamma tau_end [rtol]
"""
import sys, json, math, cmath
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# spot point (H_i units for hidden fields; inflaton in M_r units)
LX, LY, LS = 1e-3, 1.0, 1.0
MM, HH = 30.0, 40.0
MU = 7.15945616
KAP = 0.003458684721
XI = -5.0 / 6.0
H_I = 1.3334961939833566e11      # GeV
M_R = 2.435e18                   # GeV
ETA2 = (H_I / M_R) ** 2
GSTAR = 106.75
Y_Q = 4.8294009077e-9

def potential_and_forces(x, y, s):
    x2, y2, s2 = abs(x) ** 2, abs(y) ** 2, abs(s) ** 2
    V = LX * x2 * x2 + MM * MM * y2 + LY * y2 * y2 + LS * (s2 - HH * HH / 2) ** 2 \
        + 2 * (KAP * (y.conjugate() * x ** 3) + MU * (s.conjugate() * y * y)).real
    Fx = 2 * LX * x2 * x + 3 * KAP * y * x.conjugate() ** 2
    Fy = MM * MM * y + 2 * LY * y2 * y + KAP * x ** 3 + 2 * MU * s * y.conjugate()
    Fs = 2 * LS * (s2 - HH * HH / 2) * s + MU * y * y
    return V, Fx, Fy, Fs

def ricci(b, i, j, x, y, s, vx, vy, vs):
    V, Fx, Fy, Fs = potential_and_forces(x, y, s)
    K = abs(vx) ** 2 + abs(vy) ** 2 + abs(vs) ** 2
    x2 = abs(x) ** 2
    f = 1 - 2 * XI * ETA2 * x2
    num = 2 * b * b * i * i - j * j + ETA2 * (4 * V - 2 * K + 12 * XI * abs(vx) ** 2
                                              - 12 * XI * (x.conjugate() * Fx).real)
    R = num / (f + 12 * XI * XI * ETA2 * x2)
    return R, V, Fx, Fy, Fs, K, f

def make_rhs(b, gam):
    def rhs(t, u):
        i, j = u[0], u[1]
        x = complex(u[2], u[3]); y = complex(u[4], u[5]); s = complex(u[6], u[7])
        vx = complex(u[8], u[9]); vy = complex(u[10], u[11]); vs = complex(u[12], u[13])
        r, lna, h = u[14], u[15], u[16]
        R, V, Fx, Fy, Fs, K, f = ricci(b, i, j, x, y, s, vx, vy, vs)
        ax = -3 * h * vx - Fx - XI * R * x
        ay = -3 * h * vy - Fy
        as_ = -3 * h * vs - Fs
        fd = -4 * XI * ETA2 * (x.conjugate() * vx).real
        fdd = -4 * XI * ETA2 * (abs(vx) ** 2 + (x.conjugate() * ax).real)
        hd = -(j * j + 2 * ETA2 * K + 4 * r / 3 + fdd - h * fd) / (2 * f)
        return [j, -(3 * h + gam) * j - b * b * i,
                vx.real, vx.imag, vy.real, vy.imag, vs.real, vs.imag,
                ax.real, ax.imag, ay.real, ay.imag, as_.real, as_.imag,
                gam * j * j - 4 * h * r, h, hd]
    return rhs

def prepare(b):
    """Conditional minimum: X real>0 radial-stationary with the actual R(0), Y free (complex),
    Sigma = -i*sigma radial-stationary; velocities zero; H=1; inflaton amplitude from Friedmann."""
    def inflaton_amp(V):
        return math.sqrt(max(2 * (3.0 - ETA2 * V) / (b * b), 0.0))
    def eqs(p):
        x = complex(p[0], 0.0); y = complex(p[1], p[2]); s = complex(0.0, -p[3])
        V, Fx, Fy, Fs = potential_and_forces(x, y, s)
        i0 = inflaton_amp(V)
        R, *_ = ricci(b, i0, 0.0, x, y, s, 0j, 0j, 0j)
        return [(Fx + XI * R * x).real, Fy.real, Fy.imag, (s.conjugate() * Fs).real / abs(s)]
    x0 = math.sqrt(-XI * 12.0 / (2 * LX))
    y0 = -KAP * x0 ** 3 / MM ** 2
    guess = [x0, y0, 0.0, HH / math.sqrt(2)]
    sol, info, ier, msg = fsolve(eqs, guess, xtol=1e-14, full_output=True)
    res = max(abs(v) for v in eqs(sol))
    x = complex(sol[0], 0); y = complex(sol[1], sol[2]); s = complex(0, -sol[3])
    V, *_ = potential_and_forces(x, y, s)
    i0 = inflaton_amp(V)
    R0, *_ = ricci(b, i0, 0.0, x, y, s, 0j, 0j, 0j)
    u0 = [i0, 0.0, x.real, x.imag, y.real, y.imag, s.real, s.imag, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 1.0]
    return u0, res, R0

def diagnostics(u, b):
    i, j = u[0], u[1]
    x = complex(u[2], u[3]); y = complex(u[4], u[5]); s = complex(u[6], u[7])
    vx = complex(u[8], u[9]); vy = complex(u[10], u[11]); vs = complex(u[12], u[13])
    r, lna, h = u[14], u[15], u[16]
    a = math.exp(lna); a3 = a ** 3; a4 = a3 * a
    R, V, Fx, Fy, Fs, K, f = ricci(b, i, j, x, y, s, vx, vy, vs)
    fd = -4 * XI * ETA2 * (x.conjugate() * vx).real
    rho = 0.5 * j * j + 0.5 * b * b * i * i + ETA2 * (K + V) + r
    fried = (3 * f * h * h + 3 * h * fd - rho) / rho
    NX = -2 * a3 * (x.conjugate() * vx).imag
    NY = -2 * a3 * (y.conjugate() * vy).imag
    NS = -2 * a3 * (s.conjugate() * vs).imag
    gauss = NX + 3 * NY + 6 * NS
    EX = a4 * (abs(vx) ** 2 + LX * abs(x) ** 4)
    eps = NX / (2 * math.sqrt(2 * LX) * (EX / (3 * LX)) ** 0.75)
    # Y free-mode diagnostic (quoted convention: subtract forced response with instantaneous Sigma)
    A = MM * MM; B = 2 * MU * s; J = KAP * x ** 3
    D = A * A - abs(B) ** 2
    yp = (-A * J + B * J.conjugate()) / D
    jd = 3 * KAP * x * x * vx + 2 * MU * vs * yp.conjugate()
    ypd = (-A * jd + B * jd.conjugate()) / D
    rot = cmath.exp(-0.5j * cmath.phase(s))
    theta = (s.conjugate() * vs).imag / abs(s) ** 2
    dy = (y - yp) * rot
    dv = (vy - ypd - 0.5j * theta * (y - yp)) * rot
    wh = math.sqrt(A + abs(B)); wl = math.sqrt(A - abs(B))
    NYH = a3 * (dv.real ** 2 + wh * wh * dy.real ** 2) / wh
    NYL = a3 * (dv.imag ** 2 + wl * wl * dy.imag ** 2) / wl
    # yield in GeV units
    rho_r = r * M_R ** 2 * H_I ** 2
    T = (30 * rho_r / (math.pi ** 2 * GSTAR)) ** 0.25
    svis = 2 * math.pi ** 2 * GSTAR * T ** 3 / 45
    Yx = NX * H_I ** 3 / (a3 * svis) if svis > 0 else float("nan")
    w_I = (0.5 * j * j - 0.5 * b * b * i * i) / max(0.5 * j * j + 0.5 * b * b * i * i, 1e-300)
    return dict(tau=None, a=a, H=h, R=R, w_I=w_I, rad_frac=r / rho, NX=NX, gauss=gauss, EX=EX, eps=eps,
                NYH=NYH, NYL=NYL, fried=fried, T_GeV=T, Y_over_YQ=Yx / Y_Q, dilution=(r * a4 / 3) ** 0.75,
                xiR=XI * R, absX=abs(x))

def run(b, gam, tau_end, rtol=1e-10, atol=1e-13):
    u0, res, R0 = prepare(b)
    rhs = make_rhs(b, gam)
    sol = solve_ivp(rhs, (0.0, tau_end), u0, method='DOP853', rtol=rtol, atol=atol, dense_output=True)
    assert sol.success, sol.message
    # late window: last 20% in tau, uniform grid (tau-weighted mean)
    tl = np.linspace(0.8 * tau_end, tau_end, 4001)
    late = [diagnostics(sol.sol(t), b) for t in tl]
    keys = ['eps', 'NX', 'NYH', 'NYL', 'Y_over_YQ', 'EX', 'a', 'H', 'rad_frac', 'dilution', 'T_GeV']
    mean = {k: float(np.mean([d[k] for d in late])) for k in keys}
    # checks over the whole run
    tg = np.linspace(0.0, tau_end, 20001)
    allq = [diagnostics(sol.sol(t), b) for t in tg]
    checks = dict(max_friedmann=float(max(abs(d['fried']) for d in allq)),
                  max_gauss=float(max(abs(d['gauss']) for d in allq)),
                  xiR_min=float(min(d['xiR'] for d in allq)), xiR_max=float(max(d['xiR'] for d in allq)),
                  nfev=int(sol.nfev), steps=int(len(sol.t)))
    # charge-generation timing: tau at which |NX| first reaches 50/90/99 % of the late mean
    NXt = np.array([d['NX'] for d in allq]); target = mean['NX']
    timing = {}
    for frac in (0.5, 0.9, 0.99):
        idx = np.where(np.sign(target) * NXt >= frac * abs(target))[0]
        timing[f'tau_{int(frac*100)}pct'] = float(tg[idx[0]]) if len(idx) else None
        if len(idx):
            timing[f'radfrac_at_{int(frac*100)}pct'] = float(allq[idx[0]]['rad_frac'])
    final = diagnostics(sol.y[:, -1], b)
    # window dependence of the late means (fractions of tau_end): [0.1,0.125], [0.16,0.2], [0.4,0.5], [0.8,1.0]
    windows = {}
    for lo, hi in [(0.1, 0.125), (0.16, 0.2), (0.4, 0.5), (0.8, 1.0)]:
        tw = np.linspace(lo * tau_end, hi * tau_end, 2001)
        dw = [diagnostics(sol.sol(t), b) for t in tw]
        windows[f'{lo*tau_end:g}-{hi*tau_end:g}'] = {k: float(np.mean([d[k] for d in dw])) for k in ('NX', 'eps', 'Y_over_YQ', 'NYL')}
    out = dict(b=b, gamma=gam, tau_end=tau_end, rtol=rtol, prep_residual=res, R0=R0,
               initial_state=u0, late_mean=mean,
               charge_ratio_to_instant=abs(mean['NX']) / 2788.18,
               yield_from_ratio_over_D=abs(mean['NX']) / 2788.18 / mean['dilution'],
               final=dict(a=final['a'], r=sol.y[14, -1], H=final['H'], rad_frac=final['rad_frac'], T_GeV=final['T_GeV']),
               timing=timing, checks=checks, windows=windows)
    return out

if __name__ == '__main__':
    outp = sys.argv[1]; b = float(sys.argv[2]); gam = float(sys.argv[3]); tau_end = float(sys.argv[4])
    rtol = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-10
    res = run(b, gam, tau_end, rtol)
    json.dump(res, open(outp, 'w'), indent=1)
    m = res['late_mean']
    print(f"b={b} gamma={gam} tau_end={tau_end}: eps={m['eps']:.6f} NX={m['NX']:.3f} NYL={m['NYL']:.6g} NYH={m['NYH']:.4g} "
          f"yield/YQ={m['Y_over_YQ']:.5f} (ratio/D {res['yield_from_ratio_over_D']:.5f}) a_f={res['final']['a']:.4f} "
          f"T_f={res['final']['T_GeV']:.4g} GeV | prep res {res['prep_residual']:.1e} R0={res['R0']:.4f} | "
          f"Friedmann {res['checks']['max_friedmann']:.1e} Gauss {res['checks']['max_gauss']:.1e} xiR [{res['checks']['xiR_min']:.2f},{res['checks']['xiR_max']:.2f}] "
          f"| timing {res['timing']} | nfev {res['checks']['nfev']}")
    for k, v in res['windows'].items():
        print(f"   window tau {k}: NX {v['NX']:.2f} eps {v['eps']:.5f} Y/YQ {v['Y_over_YQ']:.5f} NYL {v['NYL']:.5g}")
