"""Skeptic check of m1_transport T8/T13 fluid-toy numbers.

Part A (reproduction of the REVIEWER's toy, labelled): import the reviewer's mixture.py and test whether
run_two_stage(k, a_dec -> a_i) reproduces run(k, 'clustered'), which it must by construction; then
re-run the anomalous two-stage points (k = 1e6, sigma = 0.01 gave 0.0056) with a different integrator.

Part B (own code): an independent two-fluid integration (cold + gas fluid with c_s ~ 1/a) in the
reviewer's background, written from scratch, for the no-drag 'erased'/'clustered' cases and for the
proper two-stage problem.  Same physics as the toy, different code.
"""
import json, math, sys, os
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve()
TRACK = HERE.parents[2]
RES = HERE.parents[1] / 'results'
sys.path.insert(0, str(TRACK / 'code'))
import mixture as mx          # reviewer's toy (Claude review code, not Codex)
import transport_check as tc

bg = mx.bg; a_i = mx.a_i; f = mx.fgas; cs2 = mx.v_rms**2 / 3
out = {}

# ---------------- Part A
A = {}
for k in (1e5, 3e5, 1e6):
    ref = mx.run(k, cs2, 'clustered')['P_mix_over_P_cold_z0']
    two = mx.run_two_stage(k, cs2, a_i * 1.0001)['P_mix_over_P_cold_z0']
    A[f'k_{k:g}'] = {'run_clustered': ref, 'two_stage_adec_eq_ai': two, 'ratio': two / ref}
out['A_two_stage_consistency'] = A

# re-run the anomalous point with other integrators (copy of the reviewer's equations, different solver)
def two_stage_alt(k, a_dec, method='Radau', rtol=1e-10, atol=1e-14):
    N_end = math.log(1 / a_i); N_dec = math.log(a_dec / a_i)
    aHi = mx.aH_Mpc(a_i); a_k = a_i * aHi / k if k > aHi else a_i
    slope = 1.0 / (max(math.log(a_i / a_k), 0.0) + 1.0)
    def rhs1(N, y):
        a = a_i * math.exp(N); aH = mx.aH_Mpc(a); fric = 2 + mx.dlnH_dlna(a)
        Om_a = tc.Om / a**3 / float(bg.E2(a)); d, dp = y
        return [dp, -fric * dp + 1.5 * Om_a * d - (k * math.sqrt(f * cs2) * a_i / a / aH)**2 * d]
    def rhs2(N, y, gas):
        a = a_i * math.exp(N); aH = mx.aH_Mpc(a); fric = 2 + mx.dlnH_dlna(a)
        Om_a = tc.Om / a**3 / float(bg.E2(a)); dc, dcp, dg, dgp = y
        dm = (1 - f) * dc + f * dg if gas else dc
        src = 1.5 * Om_a * dm
        return [dcp, -fric * dcp + src, dgp, -fric * dgp + src - (k * math.sqrt(cs2) * a_i / a / aH)**2 * dg]
    s1 = solve_ivp(rhs1, (0, N_dec), [1.0, slope], method=method, rtol=rtol, atol=atol)
    d, dp = s1.y[0, -1], s1.y[1, -1]
    s2 = solve_ivp(rhs2, (N_dec, N_end), [d, dp, d, dp], args=(True,), method=method, rtol=rtol, atol=atol)
    sc = solve_ivp(rhs2, (0, N_end), [1.0, slope, 0, 0], args=(False,), method=method, rtol=rtol, atol=atol)
    y = s2.y[:, -1]
    return ((1 - f) * y[0] + f * y[2]) ** 2 / sc.y[0, -1] ** 2, {'s1': s1.success, 's2': s2.success, 'sc': sc.success, 'n2': int(s2.t.size)}

B = {}
for sig in (0.01, 1.0, 100.0):
    GH = 1.2141071422639381e-18 * sig * mx.v_rms * tc.mX_GeV / mx.H_i
    a_dec = a_i * math.sqrt(GH)
    for k in (3e5, 1e6):
        rev = mx.run_two_stage(k, cs2, a_dec)['P_mix_over_P_cold_z0']
        alt_r, info_r = two_stage_alt(k, a_dec, 'Radau')
        alt_d, info_d = two_stage_alt(k, a_dec, 'DOP853', 1e-10, 1e-14)
        B[f'sigma_{sig:g}_k_{k:g}'] = {'reviewer_two_stage': rev, 'Radau_tight': alt_r, 'DOP853_tight': alt_d, 'info': [info_r, info_d]}
        print('B', sig, k, rev, alt_r, alt_d, flush=True)
out['B_two_stage_resolvers'] = B

# ---------------- Part B: own two-fluid integration, own equations (ln a), reviewer's background only
def own_two_fluid(k, init='erased', stage_adec=None):
    """Return P_mix/P_cold at a = 1.  cold + gas (c_s = c_s,i a_i/a); optional tight-coupling stage with c_eff^2 = f c_s^2."""
    N_end = math.log(1 / a_i)
    aHi = mx.aH_Mpc(a_i); a_k = a_i * aHi / k if k > aHi else a_i
    slope = 1.0 / (max(math.log(a_i / a_k), 0.0) + 1.0)
    def coef(N):
        a = a_i * math.exp(N); aH = mx.aH_Mpc(a); fric = 2 + mx.dlnH_dlna(a)
        Om_a = tc.Om / a**3 / float(bg.E2(a)); kc = (k * math.sqrt(cs2) * a_i / a / aH) ** 2
        return fric, Om_a, kc
    def rhs(N, y, mode):
        fric, Om_a, kc = coef(N); dc, dcp, dg, dgp = y
        if mode == 'cold':
            return [dcp, -fric * dcp + 1.5 * Om_a * dc, 0, 0]
        if mode == 'tight':
            src = 1.5 * Om_a * dc
            return [dcp, -fric * dcp + src - f * kc * dc, dcp, -fric * dcp + src - f * kc * dc]
        src = 1.5 * Om_a * ((1 - f) * dc + f * dg)
        return [dcp, -fric * dcp + src, dgp, -fric * dgp + src - kc * dg]
    y0 = [1.0, slope, 1.0 if init == 'clustered' else 0.0, slope if init == 'clustered' else 0.0]
    if stage_adec is not None:
        N_dec = math.log(stage_adec / a_i)
        s1 = solve_ivp(rhs, (0, N_dec), [1.0, slope, 1.0, slope], args=('tight',), method='Radau', rtol=1e-10, atol=1e-14)
        y1 = s1.y[:, -1]; y1 = [y1[0], y1[1], y1[0], y1[1]]
        s2 = solve_ivp(rhs, (N_dec, N_end), y1, args=('two',), method='Radau', rtol=1e-10, atol=1e-14)
        y = s2.y[:, -1]
    else:
        s2 = solve_ivp(rhs, (0, N_end), y0, args=('two',), method='Radau', rtol=1e-10, atol=1e-14)
        y = s2.y[:, -1]
    sc = solve_ivp(rhs, (0, N_end), [1.0, slope, 0, 0], args=('cold',), method='Radau', rtol=1e-10, atol=1e-14)
    return ((1 - f) * y[0] + f * y[2]) ** 2 / sc.y[0, -1] ** 2

C = {}
for k in (3e4, 1e5, 3e5, 1e6):
    C[f'k_{k:g}'] = {'erased': own_two_fluid(k, 'erased'), 'clustered': own_two_fluid(k, 'clustered'),
                     'reviewer_isothermal': {3e4: (0.8857, 0.9905), 1e5: (0.8838, 0.9116), 3e5: (0.8691, 0.8881), 1e6: (0.7985, 0.7954)}[k],
                     'codex_BGK': {3e4: 0.979835, 1e5: 0.893210, 3e5: 0.849806, 1e6: 0.764169}[k]}
    print('C', k, C[f'k_{k:g}'], flush=True)
out['C_own_two_fluid_no_drag'] = C
D = {}
for sig in (1.0, 100.0, 1000.0):
    GH = 1.2141071422639381e-18 * sig * mx.v_rms * tc.mX_GeV / mx.H_i
    a_dec = a_i * math.sqrt(GH)
    D[f'sigma_{sig:g}'] = {f'k_{k:g}': own_two_fluid(k, 'clustered', stage_adec=a_dec) for k in (1e5, 3e5, 1e6)}
    D[f'sigma_{sig:g}']['x_dec'] = math.sqrt(GH)
    print('D', sig, D[f'sigma_{sig:g}'], flush=True)
out['D_own_two_stage'] = D
out['plateau_(1-f)^2'] = (1 - f) ** 2
RES.mkdir(parents=True, exist_ok=True)
(RES / 'sk_mixture_check.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
print(json.dumps(out, indent=1, default=float))
