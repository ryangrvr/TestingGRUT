"""st_rerun_skeptic driver for skeptic_tangent.c (own code; see the C header for what is ours).

Generates hess_own.h from a sympy Hessian of the potential in six canonical real fields, builds the C engine in
scratch, prepares the seed from OUR preparation (rp.prepare, central difference +-1e-5, zero velocities), runs, and
detects first crossings by linear interpolation between uniform output samples.
Usage: python3 -B run_tangent.py <scratch_dir> <tag> <tau_end> <kset> <rtol> <dtout>
  kset: main (Codex's 18 k + 0.005, 0.01, 2.04, 2.16, 2.18, 2.2) | band (1.80..2.60 step 0.01) | k213 (2.13 only)
"""
import sys, json, time, subprocess, hashlib
from pathlib import Path
import numpy as np
import sympy as sp
scratch, tag, tau_end, kset, rtol, dtout = Path(sys.argv[1]), sys.argv[2], float(sys.argv[3]), sys.argv[4], float(sys.argv[5]), float(sys.argv[6])
HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3] / 'rp_reheating' / 'code'))
import rp  # noqa: E402  own preparation, parameters
P = rp.P
build = scratch / 'build'; build.mkdir(parents=True, exist_ok=True)

# ---- sympy Hessian -> C ----
f = sp.symbols('f0:6', real=True)
lx, ly, ls, kap, MM, hh, mu = sp.symbols('lx ly ls kap MM hh mu', real=True)
x = (f[0] + sp.I * f[1]) / sp.sqrt(2); y = (f[2] + sp.I * f[3]) / sp.sqrt(2); s = (f[4] + sp.I * f[5]) / sp.sqrt(2)
ab2 = lambda w: sp.expand(w * sp.conjugate(w))
V = sp.expand(lx * ab2(x) ** 2 + MM ** 2 * ab2(y) + ly * ab2(y) ** 2 + ls * (ab2(s) - hh ** 2 / 2) ** 2
              + 2 * sp.re(sp.expand(kap * sp.conjugate(y) * x ** 3 + mu * sp.conjugate(s) * y * y)))
Hm = sp.hessian(V, f)
lines = ['#ifndef M_SQRT2', '#define M_SQRT2 1.41421356237309504880', '#endif', 'static void hess_own(const double *f, double Hs[6][6]) {',
         '  const double f0=f[0],f1=f[1],f2=f[2],f3=f[3],f4=f[4],f5=f[5];']
for i in range(6):
    for j in range(i, 6):
        lines.append(f'  Hs[{i}][{j}] = Hs[{j}][{i}] = {sp.ccode(sp.simplify(Hm[i, j]))};')
lines.append('}')
(build / 'hess_own.h').write_text('\n'.join(lines) + '\n')
src = HERE.parent / 'skeptic_tangent.c'
key = hashlib.sha256(src.read_bytes() + (build / 'hess_own.h').read_bytes()).hexdigest()[:10]
exe = build / f'skeptic_tangent_{key}'
if not exe.exists():
    subprocess.run(['gcc', '-O3', '-std=c11', '-I', str(build), '-o', str(exe), str(src), '-lm'], check=True)

# ---- numerical check of the generated Hessian against a finite difference of rp.potential's force ----
def grad_phi(phi):
    z = np.array([phi[0] + 1j * phi[1], phi[2] + 1j * phi[3], phi[4] + 1j * phi[5]]) / np.sqrt(2)
    _, F = rp.potential(z); g = np.sqrt(2) * F
    return np.array([g[0].real, g[0].imag, g[1].real, g[1].imag, g[2].real, g[2].imag])
Hnum = sp.lambdify([f], Hm.subs({lx: P['lx'], ly: P['ly'], ls: P['ls'], kap: P['kap'], MM: P['M'], hh: P['h'], mu: P['mu']}), 'numpy')
phit = np.array([99.1, 5.3, -2.9, 1.2, .4, -39.7]); hstep = 1e-4
FD = np.array([(grad_phi(phit + hstep * e) - grad_phi(phit - hstep * e)) / (2 * hstep) for e in np.eye(6)]).T
hess_check = float(np.max(np.abs(np.array(Hnum(phit), float) - FD)) / np.max(np.abs(FD)))

# ---- seed ----
beta, decay = (3., .5) if tag == 'quiet3' else (5., 1.)
def zprep(dl):
    z, info = rp.prepare(beta, delta=dl); assert info['residual'] < 1e-8, info; return z
z0 = zprep(0.); zp = zprep(1e-5); zm = zprep(-1e-5); dz = (zp - zm) / 2e-5
seed = np.sqrt(2) * np.array([dz[0].real, dz[0].imag, dz[1].real, dz[1].imag, dz[2].real, dz[2].imag])
KS18 = [.02, .05, .1, .2, .4, .8, 1.2, 1.6, 2., 2.13, 2.4, 2.8, 3.2, 4., 5.6, 5.8, 6., 8.]
ks = {'main': np.array(KS18 + [.005, .01, 2.04, 2.16, 2.18, 2.2]), 'band': np.round(np.arange(1.80, 2.60001, .01), 4),
      'k213': np.array([2.13])}[kset]
nk = len(ks)
pars = np.array([beta, decay, rp.XI, rp.ETA2, P['lx'], P['ly'], P['ls'], P['kap'], P['M'], P['h'], P['mu'], tau_end, rtol, rtol * 1e-3,
                 dtout, nk, 0.])
label = f'{tag}_{kset}_{tau_end:g}_r{rtol:g}_d{dtout:g}'
inp = scratch / f'{label}.in'; outb = scratch / f'{label}.bin'
with inp.open('wb') as fh:
    pars.astype('<f8').tofile(fh)
    np.array([z0[0].real, z0[0].imag, z0[1].real, z0[1].imag, z0[2].real, z0[2].imag], '<f8').tofile(fh)
    ks.astype('<f8').tofile(fh); seed.astype('<f8').tofile(fh)
t0 = time.time()
p = subprocess.run([str(exe), str(inp), str(outb)], capture_output=True, text=True)
if p.returncode: raise SystemExit(f'engine failed {p.returncode} {p.stderr}')
secs = time.time() - t0
raw = np.fromfile(outb, '<f8').reshape(-1, 6 + 2 * nk)
T, A, H, NX, Rtr, Rge = raw[:, :6].T
amp = raw[:, 6::2]; dN = raw[:, 7::2]
w = (T >= 64) & (T <= 80)
norm_run = abs(np.trapezoid(NX[w], T[w]) / (T[w][-1] - T[w][0])) if w.sum() > 10 else None
norm = abs(json.loads((HERE.parents[1] / 'results' / 'bg_own.json').read_text())[tag]['NX_mean_64_80_dense'])
chg = dN / norm
np.savez_compressed(scratch / f'{label}.npz', t=T, a=A, H=H, NX=NX, amp=amp, charge=chg, ks=ks, seed=seed)
outb.unlink()

def first_cross(series, level, idx):
    m = np.max(np.abs(series[:, idx]), axis=1); wv = np.flatnonzero(m >= level)
    if not len(wv): return None
    i = int(wv[0]); j = idx[int(np.argmax(np.abs(series[i, idx])))]
    y1, y2 = abs(series[i - 1, j]), abs(series[i, j])
    tc = T[i - 1] + (level - y1) / (y2 - y1) * (T[i] - T[i - 1])
    return dict(tau=float(tc), k=float(ks[j]), a=float(np.interp(tc, T, A)), bracket=[float(T[i - 1]), float(T[i])])

grids = {'all': np.arange(nk)}
if kset == 'main': grids['18k'] = np.array([int(np.argmin(abs(ks - k))) for k in KS18])
res = dict(tag=tag, kset=kset, rtol=rtol, dtout=dtout, tau_end=tau_end, ks=ks.tolist(), engine=p.stdout.strip(), seconds=secs,
           hessian_vs_fd_of_rp_force=hess_check, seed=seed.tolist(), amplitude_per_radian_t0=float(amp[0, 0]),
           norm_NX_64_80=float(norm), norm_this_run=norm_run, max_R_trace_vs_geom_rel=float(np.max(np.abs(Rtr - Rge)) / np.max(np.abs(Rtr))),
           a_at={str(tt): float(np.interp(tt, T, A)) for tt in (56.46, 80., 134., 147.5, 772., 782., 800.) if tt <= tau_end},
           crossings={})
for gname, idx in grids.items():
    for sd in [.03, 1e-6, 1e-9]:
        res['crossings'][f'{gname}_{sd:g}'] = dict(
            amp10=first_cross(sd * amp, .1, idx), amp1=first_cross(sd * amp, 1., idx), charge1=first_cross(sd * chg, 1., idx),
            max_amp=float(sd * np.abs(amp[:, idx]).max()), max_charge=float(sd * np.abs(chg[:, idx]).max()))
    w80 = T <= 80
    sub = np.abs(amp[w80][:, idx]); ii = np.unravel_index(sub.argmax(), sub.shape)
    subc = np.abs(chg[w80][:, idx]); jj = np.unravel_index(subc.argmax(), subc.shape)
    res[f'through80_003_{gname}'] = dict(max_amp=float(.03 * sub.max()), tau=float(T[w80][ii[0]]), k=float(ks[idx][ii[1]]),
                                         max_charge=float(.03 * subc.max()), tau_c=float(T[w80][jj[0]]), k_c=float(ks[idx][jj[1]]))
out = HERE.parents[1] / 'results' / f'tangent_{label}.json'
out.write_text(json.dumps(res, indent=1, default=float) + '\n')
print(json.dumps({k: v for k, v in res.items() if k not in ('seed', 'ks')}, indent=1, default=float))
