#!/usr/bin/env python3
"""m1_holder track (Claude, own code) -- part 3 (charge change with the Option 1A holder) and the time-shift test of
part 4a, driving rpholder.c (our rp_reheating engine with mu_eff(R)).

Runs (all mode 0, m_I = 3 H_i, Gamma_I = 0.5 H_i, H_i t <= 800, late window = last 20 % as in Codex and rp_reheating):
  original      lock = 0, M1 preparation state (own Newton solve, rp_reheating.prepare)      -> Codex reference.json
  holder        lock = (10 i - mu)/12, own stationary point (results/stationary_point.json) -> Codex holder.json
  controls      rtol 1e-12; end 1000 (window 800-1000); half window
  response row  d ln|N_X| / d(delta r_a), d ln|N_X| / d(delta rdot_a) along the five physical Hessian eigenvectors
                (central differences, two step sizes)                                       -> Codex perturbations.json
  phase         d ln|N_X| / d(arg Sigma) at fixed amplitudes and velocities                 -> Codex 0.554384
  time shift    restart from the state at t0 and at t0 + eps: full shift, scalar-only shift, background-only shift
Reproduction label: none of Codex's code is run; Codex's numbers are quoted [Q] for comparison.
"""
import json, os, subprocess, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
SP = Path(os.environ.get('M1H_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/m1_holder'))
BIN = SP / 'build' / 'rpholder'
sys.path.insert(0, str(HERE.parent.parent / 'rp_reheating' / 'code'))
import rp                                                     # our own rp_reheating driver: prepare() only (read-only use)

H_GEV = 134284273491.30815; MPL = 2.435e18
ETA2 = (H_GEV / MPL) ** 2
P = dict(lx=1e-3, ly=1.0, ls=1.0, kap=0.0034586847211042523, M=30.0, h=40.0, mu=961402368973.181 / H_GEV, cH=10.0, xi=-5 / 6)
LOCK = (10j - P['mu']) / 12
B_I, GAM = 3.0, 0.5
COLS = rp.COLS
CODEX = dict(ref=-8120.739695645271, ref_half=-8118.6145409703995, holder=-8267.43703207837, holder_half=-8267.128654142336,
             ratio=1.0180645288398746, H_ratio=0.9881354099844972, phase_tangent=0.5543843049360061,
             phase_fd=[0.5543843710942298, 0.5543843215267685])


def build():
    BIN.parent.mkdir(parents=True, exist_ok=True)
    src = HERE / 'rpholder.c'
    if not BIN.exists() or BIN.stat().st_mtime < src.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-std=c11', '-o', str(BIN), str(src), '-lm'], check=True)


def run(tag, lock, z0=None, v0=None, state=None, end=800.0, rtol=1e-10, atol=1e-13, nsamp=1600, save=False):
    g = lambda v: '%.17g' % float(v)
    args = [str(BIN), '0', g(B_I), g(GAM), g(P['xi']), g(ETA2), g(P['cH']), g(P['lx']), g(P['ly']), g(P['ls']), g(P['kap']), g(P['M']), g(P['h']), g(P['mu']),
            g(lock.real), g(lock.imag), g(end), g(rtol), g(atol), str(int(nsamp))]
    if state is None:
        v0 = np.zeros(3, complex) if v0 is None else np.asarray(v0, complex)
        args += ['F'] + [g(v) for c in z0 for v in (c.real, c.imag)] + [g(v) for c in v0 for v in (c.real, c.imag)]
    else:
        args += ['S'] + [g(v) for v in state]
    t0 = time.time()
    out = subprocess.run(args, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(tag + ': ' + out.stderr)
    S, W, Y, END = [], {}, None, None
    for line in out.stdout.split('\n'):
        if line.startswith('S '):
            S.append([float(v) for v in line.split()[1:]])
        elif line.startswith('W'):
            v = line.split(); W[v[0]] = dict(span=float(v[1]), **dict(zip(COLS, map(float, v[2:]))))
        elif line.startswith('Y '):
            Y = [float(v) for v in line.split()[1:]]
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), rejects=int(v[2]), max_gauss=float(v[3]), max_H_consistency=float(v[4]), max_hidden_work=float(v[5]))
    S = np.array(S)
    if save:
        (SP / 'runs').mkdir(parents=True, exist_ok=True); np.save(SP / 'runs' / (tag + '.npy'), S)
    res = dict(tag=tag, end=end, rtol=rtol, lock=[lock.real, lock.imag], NX=W['W80']['NX'], NX_half=W['W90']['NX'], eps=W['W80']['eps'],
               NYL=W['W80']['NYL'], NYH=W['W80']['NYH'], R_range=[float(S[:, 3].min()), float(S[:, 3].max())], final_a=float(S[-1, 1]),
               rho_r_a4_end=float(S[-1, 16] * S[-1, 1] ** 4), checks=END, seconds=time.time() - t0)
    return res, S, Y


def main():
    build()
    out = {}
    stat = json.load(open(RES / 'stationary_point.json'))
    zs = np.array([complex(a, b) for a, b in stat['z_stationary']])
    U = np.array(stat['eigenvectors_columns']); ev = np.array(stat['eigenvalues'])
    phys = [k for k in range(6) if abs(ev[k]) > 1e-6]
    zprep, info = rp.prepare(B_I)                               # own solve of the M1 preparation point (lock = 0)
    out['preparation_state_own'] = [[c.real, c.imag] for c in zprep]
    out['preparation_residual'] = info['residual']

    # ---- original and holder
    ref, Sref, _ = run('original', 0j, zprep, save=True)
    hol, Shol, _ = run('holder', LOCK, zs, save=True)
    hol_t, _, _ = run('holder_rtol12', LOCK, zs, rtol=1e-12, atol=1e-15)
    ref_t, _, _ = run('original_rtol12', 0j, zprep, rtol=1e-12, atol=1e-15)
    hol_l, _, _ = run('holder_end1000', LOCK, zs, end=1000.0, nsamp=2000)
    ref_l, _, _ = run('original_end1000', 0j, zprep, end=1000.0, nsamp=2000)
    ratio = hol['NX'] / ref['NX']
    out['original'] = ref; out['holder'] = hol; out['controls'] = dict(holder_rtol12=hol_t, original_rtol12=ref_t, holder_end1000=hol_l, original_end1000=ref_l)
    out['comparison'] = dict(
        NX_original_own=ref['NX'], NX_original_codex=CODEX['ref'], original_rel_diff=ref['NX'] / CODEX['ref'] - 1,
        NX_half_original_own=ref['NX_half'], NX_half_original_codex=CODEX['ref_half'],
        NX_holder_own=hol['NX'], NX_holder_codex=CODEX['holder'], holder_rel_diff=hol['NX'] / CODEX['holder'] - 1,
        NX_half_holder_own=hol['NX_half'], NX_half_holder_codex=CODEX['holder_half'],
        charge_ratio_own=ratio, charge_ratio_codex=CODEX['ratio'], charge_change_percent_own=100 * (ratio - 1), charge_change_percent_codex=100 * (CODEX['ratio'] - 1),
        charge_ratio_rtol12=hol_t['NX'] / ref_t['NX'], charge_ratio_window_800_1000=hol_l['NX'] / ref_l['NX'], charge_ratio_half_window=hol['NX_half'] / ref['NX_half'],
        H_ratio_own=(1 / ratio) ** (2 / 3), H_ratio_codex=CODEX['H_ratio'], H_new_GeV_own=H_GEV * (1 / ratio) ** (2 / 3),
        H_ratio_fully_cooled_own=(1 / ratio / 1.008409936530420273) ** (2 / 3))
    # early profile of the two runs (for the notes): charge build-up, R, |Y|, phase
    def profile(S):
        keep = ['tau', 'a', 'R', 'NX', 'absX', 'absY', 'S_over_vev', 'radfrac', 'eps']
        idx = [COLS.index(c) for c in keep]
        early = S[S[:, 0] <= 20.0][::2][:, idx]
        return dict(columns=keep, rows=early.tolist())
    out['profiles'] = dict(original=profile(Sref), holder=profile(Shol))
    # fraction of final charge reached by tau = 1, 2, 3, 5, 10
    for nm, S in (('original', Sref), ('holder', Shol)):
        NXf = out[nm]['NX']
        out[nm]['charge_fraction_at_tau'] = {str(t): float(np.interp(t, S[:, 0], S[:, 4]) / NXf) for t in (0.5, 1, 2, 3, 5, 10, 20)}
        out[nm]['R_at_tau'] = {str(t): float(np.interp(t, S[:, 0], S[:, 3])) for t in (0.5, 1, 2, 3, 5)}

    # ---- response row along the physical eigenvectors (canonical real fields r = sqrt2 (Re, Im) pairs)
    def perturb(dr):
        dz = np.array([dr[0] + 1j * dr[1], dr[2] + 1j * dr[3], dr[4] + 1j * dr[5]]) / np.sqrt(2)
        return dz
    N0 = hol['NX']
    rows = {}
    for step in (1e-2, 1e-3):
        pos, vel = [], []
        for a in phys:
            e = U[:, a]
            Np, _, _ = run(f'pos{a}+', LOCK, zs + perturb(step * e)); Nm, _, _ = run(f'pos{a}-', LOCK, zs - perturb(step * e))
            pos.append((np.log(abs(Np['NX'])) - np.log(abs(Nm['NX']))) / (2 * step))
            Np, _, _ = run(f'vel{a}+', LOCK, zs, v0=perturb(step * e)); Nm, _, _ = run(f'vel{a}-', LOCK, zs, v0=-perturb(step * e))
            vel.append((np.log(abs(Np['NX'])) - np.log(abs(Nm['NX']))) / (2 * step))
        rows[str(step)] = dict(position=pos, velocity=vel)
    # Codex's canonical 12-row [Q perturbations.json] projected on OUR eigenvectors (same Hessian, sign-consistent)
    codex_row = np.array([0.008727474202496373, -0.03448797758673666, -0.0032932221196985086, 0.0010766259369597193, 0.014387590645153131, -0.00013868808357356228,
                          0.006290799639779052, 0.007656999917670617, -0.0009038792060610946, -0.0003518727880905591, -0.000578227596296392, -3.190277600716975e-05])
    Up = U[:, phys]
    proj = np.block([[Up, np.zeros_like(Up)], [Up, Up]])
    codex_long = codex_row @ proj
    own = np.r_[rows['0.001']['position'], rows['0.001']['velocity']]
    own2 = np.r_[rows['0.01']['position'], rows['0.01']['velocity']]
    out['response_row'] = dict(physical_eigenvalues=ev[phys].tolist(), own_fd_step_1e3=own.tolist(), own_fd_step_1e2=own2.tolist(),
                               own_richardson=((4 * own - own2) / 3).tolist(), codex_on_own_basis=codex_long.tolist(),
                               rel_diff_1e3=(own / codex_long - 1).tolist(), abs_diff_1e3=(own - codex_long).tolist(),
                               norm_rel_diff=float(np.linalg.norm(own - codex_long) / np.linalg.norm(codex_long)),
                               convention='position: d ln|N_X| / d(delta r) at zero velocity; velocity: d ln|N_X| / d(delta rdot) at zero displacement; '
                                          'r = sqrt2 (Re X, Im X, Re Y, Im Y, Re Sigma, Im Sigma) in H_i units, unit eigenvectors of the R = 12 Hessian')
    # gauge direction (should give zero response)
    gidx = [k for k in range(6) if abs(ev[k]) <= 1e-6][0]
    e = U[:, gidx]; step = 1e-2
    Np, _, _ = run('gauge+', LOCK, zs + perturb(step * e)); Nm, _, _ = run('gauge-', LOCK, zs - perturb(step * e))
    out['response_row']['gauge_direction_response'] = float((np.log(abs(Np['NX'])) - np.log(abs(Nm['NX']))) / (2 * step))

    # ---- Sigma phase sensitivity
    fd = {}
    for step in (1e-3, 5e-4):
        vals = []
        for sgn in (-1, 1):
            zz = zs.copy(); zz[2] *= np.exp(1j * sgn * step)
            r, _, _ = run(f'phase{sgn}{step}', LOCK, zz); vals.append(r['NX'])
        fd[str(step)] = float((np.log(abs(vals[1])) - np.log(abs(vals[0]))) / (2 * step))
    # the same for the original preparation (Codex d4747cb: -0.21996870)
    fdo = {}
    for step in (1e-3,):
        vals = []
        for sgn in (-1, 1):
            zz = zprep.copy(); zz[2] *= np.exp(1j * sgn * step)
            r, _, _ = run(f'ophase{sgn}{step}', 0j, zz); vals.append(r['NX'])
        fdo[str(step)] = float((np.log(abs(vals[1])) - np.log(abs(vals[0]))) / (2 * step))
    out['phase_sensitivity'] = dict(holder_own=fd, holder_codex_tangent=CODEX['phase_tangent'], holder_codex_fd=CODEX['phase_fd'],
                                    original_own=fdo, original_codex_d4747cb=-0.21996870110818634)

    # ---- time-shift (adiabatic mode) test: restart from the state at t0 and t0 + eps
    shift = {}
    for t0 in (3.0, 6.0, 12.0):
        for eps in (1e-2, 1e-3):
            _, _, Y0 = run('ts_a', LOCK, zs, end=t0, nsamp=10)
            _, _, Y1 = run('ts_b', LOCK, zs, end=t0 + eps, nsamp=10)
            Y0 = np.array(Y0); Y1 = np.array(Y1)
            dur = 800.0 - t0
            A, _, YA = run('ts_full0', LOCK, state=Y0, end=dur)
            Bf, _, YB = run('ts_full1', LOCK, state=Y1, end=dur)
            Yc = Y0.copy(); Yc[3:15] = Y1[3:15]                 # scalar fields and velocities shifted, background not
            C, _, YC = run('ts_scalar', LOCK, state=Yc, end=dur)
            Yd = Y1.copy(); Yd[3:15] = Y0[3:15]                 # background shifted, scalars not
            D, _, YD = run('ts_background', LOCK, state=Yd, end=dur)
            # comoving radiation entropy proxy at the end: (rho_r a^4)^(3/4) (g* constant)
            def ent(Y):
                return (Y[15] * np.exp(4 * Y[0])) ** 0.75
            shift[f't0={t0:g},eps={eps:g}'] = dict(
                NX_unshifted=A['NX'], NX_full_shift=Bf['NX'], NX_scalar_only_shift=C['NX'], NX_background_only_shift=D['NX'],
                dlnNX_deps_full=(Bf['NX'] / A['NX'] - 1) / eps, dlnNX_deps_scalar_only=(C['NX'] / A['NX'] - 1) / eps,
                dlnNX_deps_background_only=(D['NX'] / A['NX'] - 1) / eps,
                dlnS_deps_full=(ent(YB) / ent(YA) - 1) / eps, dlnS_deps_scalar_only=(ent(YC) / ent(YA) - 1) / eps, dlnS_deps_background_only=(ent(YD) / ent(YA) - 1) / eps,
                charge_entropy_perturbation_full=(Bf['NX'] / A['NX'] - 1) / eps - (ent(YB) / ent(YA) - 1) / eps,
                charge_entropy_perturbation_scalar_only=(C['NX'] / A['NX'] - 1) / eps - (ent(YC) / ent(YA) - 1) / eps,
                scalar_velocity_norm_at_t0=float(np.linalg.norm(Y0[9:15])), state_shift_norm=float(np.linalg.norm(Y1[3:15] - Y0[3:15])),
                sum_scalar_plus_background=(C['NX'] / A['NX'] - 1) / eps + (D['NX'] / A['NX'] - 1) / eps)
    out['time_shift_test'] = shift
    json.dump(out, open(RES / 'charge_transfer.json', 'w'), indent=1, default=float)
    print(json.dumps({k: out[k] for k in ('comparison', 'response_row', 'phase_sensitivity', 'time_shift_test')}, indent=1, default=float))
    print('checks original', ref['checks'], 'holder', hol['checks'], 'R range', ref['R_range'], hol['R_range'])
    print('charge fractions', out['original']['charge_fraction_at_tau'], out['holder']['charge_fraction_at_tau'])


if __name__ == '__main__':
    main()
