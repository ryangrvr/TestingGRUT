"""en_rerun_skeptic -- own diagnostics and reference projections (no Codex import).
rfftn with explicit half-space weights (en_rerun used a full complex fftn; Codex uses rfftn).

Projections, all divided by the total net current density Q:
  n_-(omega) = 1/2 [ (|pi|^2 + omega^2 |chi|^2)/(2 omega) - Re(chi1* pi2 - chi2* pi1) ]  summed over k != 0
  references: massless omega^2 = k_eff^2; Codex isotropic Hartree 2 lX <r^2> + g2 <phi^2>; variants.
Floors (lowest n_- reachable inside a family of U(1)-covariant equal-time quadratic projections):
  'w'  : one omega per group                  -> 1/2 [ sqrt(sum A * sum B) - sum Qm ]
  'wg' : one omega and one chirp gamma per group, a = (omega chi + i(pi - gamma chi))/sqrt(2 omega)
         (pi -> pi - gamma chi is a canonical shear; the charge term is invariant under it)
                                              -> 1/2 [ sqrt(sum A sum B - (sum R)^2) - sum Qm ]
  with A = |pi1|^2+|pi2|^2, B = |chi1|^2+|chi2|^2, R = Re(chi1* pi1 + chi2* pi2), Qm as above.
  Groups: single Fourier mode; FD4 shell (equal FD4 symbol = cubic orbit, en_rerun's choice); isotropic |n|^2 shell.
"""
import numpy as np
from scipy import fft as sfft

LX, G2, LP = .001, .01, .03
_cache = {}


def fd4_1d(N, L, real_last=False):
    x = 2*np.pi*(np.fft.rfftfreq(N) if real_last else np.fft.fftfreq(N))
    return 4*np.sin(x/2)**2*(1 + np.sin(x/2)**2/3)/(L/N)**2


def fd4_symbol_r(N, L):
    a = fd4_1d(N, L); b = fd4_1d(N, L, True)
    return a[:, None, None] + a[None, :, None] + b[None, None, :]


def geometry(N, L):
    key = (N, L)
    if key in _cache:
        return _cache[key]
    sym = fd4_symbol_r(N, L)
    n = np.rint(np.fft.fftfreq(N)*N); nl = np.rint(np.fft.rfftfreq(N)*N)
    n2 = (n[:, None, None]**2 + n[None, :, None]**2 + nl[None, None, :]**2).astype(int)
    w = np.full(N//2 + 1, 2.); w[0] = 1.; w[-1] = 1.
    W = np.broadcast_to(w[None, None, :], sym.shape).copy(); W[0, 0, 0] = 0.
    kphys = 2*np.pi/L*np.sqrt(n2)
    nz = W > 0
    _, fd4_id = np.unique(np.round(sym[nz], 9), return_inverse=True)
    _, n2_id = np.unique(n2[nz], return_inverse=True)
    g = dict(sym=sym, W=W, kphys=kphys, nz=nz, fd4_id=fd4_id, n2_id=n2_id)
    _cache[key] = g
    return g


def diagnostics(q, p, L, sym, s, tau, proj=True):
    N = q.shape[-1]
    qh = sfft.rfftn(q, axes=(1, 2, 3), norm='forward', workers=2)
    ph = sfft.rfftn(p, axes=(1, 2, 3), norm='forward', workers=2)
    g = geometry(N, L); W = g['W']; nz = g['nz']
    Wf = W.copy(); Wf[0, 0, 0] = 1.  # include k=0 for totals
    r2 = q[0]**2 + q[1]**2; f2 = q[2]**2
    grad = .5*np.sum(Wf*sym*np.abs(qh)**2, axis=(1, 2, 3))
    j = q[0]*p[1] - q[1]*p[0]
    row = dict(eta=s, tau=tau,
               KX=float(.5*np.mean(p[0]**2 + p[1]**2)), Kphi=float(.5*np.mean(p[2]**2)),
               GX=float(grad[0] + grad[1]), Gphi=float(grad[2]),
               UX=float(np.mean(LX/4*r2*r2)), Uphi=float(np.mean(LP/4*f2*f2)), Umix=float(np.mean(G2/2*f2*r2)),
               nQ=float(np.mean(j)), charge_RMS=float(np.sqrt(np.mean(j*j))),
               q0=float(q[0].mean()), q1=float(q[1].mean()), qphi=float(q[2].mean()),
               p0=float(p[0].mean()), p1=float(p[1].mean()), pphi=float(p[2].mean()),
               X2=float(r2.mean()), phi2=float(f2.mean()))
    row['rho'] = row['KX'] + row['Kphi'] + row['GX'] + row['Gphi'] + row['UX'] + row['Uphi'] + row['Umix']
    if not proj:
        return row
    Qt = row['nQ']
    A = (np.abs(ph[0])**2 + np.abs(ph[1])**2)[nz]
    B = (np.abs(qh[0])**2 + np.abs(qh[1])**2)[nz]
    R = np.real(np.conj(qh[0])*ph[0] + np.conj(qh[1])*ph[1])[nz]
    Qm = np.real(np.conj(qh[0])*ph[1] - np.conj(qh[1])*ph[0])[nz]
    w = W[nz]; k2 = sym[nz]; kp = g['kphys'][nz]
    mbar2 = row['q0']**2 + row['q1']**2
    masses = dict(massless=0., hartree_codex=2*LX*row['X2'] + G2*row['phi2'],
                  hartree_no_zero_mode=2*LX*(row['X2'] - mbar2) + G2*row['phi2'],
                  tangential=LX*row['X2'] + G2*row['phi2'], radial=3*LX*row['X2'] + G2*row['phi2'])
    refs = {}
    for name, m2 in masses.items():
        om = np.sqrt(k2 + m2)
        nm = w*((A + om*om*B)/(2*om) - Qm)/2
        refs[name] = dict(m2=m2, total=float(nm.sum()/Qt), k_lt_1p8=float(nm[kp < 1.8].sum()/Qt),
                          k_lt_1=float(nm[kp < 1.].sum()/Qt), k_ge_1p8=float(nm[kp >= 1.8].sum()/Qt))
    floors = {}
    floors['mode_w'] = float(np.sum(w*(np.sqrt(A*B) - Qm)/2)/Qt)
    floors['mode_wg'] = float(np.sum(w*(np.sqrt(np.maximum(A*B - R*R, 0)) - Qm)/2)/Qt)
    for gname in ('fd4_id', 'n2_id'):
        gid = g[gname]
        As = np.bincount(gid, w*A); Bs = np.bincount(gid, w*B); Rs = np.bincount(gid, w*R); Qs = np.bincount(gid, w*Qm)
        floors[gname[:-3] + '_shell_w'] = float(np.sum((np.sqrt(As*Bs) - Qs)/2)/Qt)
        floors[gname[:-3] + '_shell_wg'] = float(np.sum((np.sqrt(np.maximum(As*Bs - Rs*Rs, 0)) - Qs)/2)/Qt)
        # IR part of the shell floors
        kg = np.bincount(gid, w*kp)/np.bincount(gid, w)
        floors[gname[:-3] + '_shell_w_k_lt_1p8'] = float(np.sum(((np.sqrt(As*Bs) - Qs)/2)[kg < 1.8])/Qt)
    row.update(mean_r2=mbar2, refs=refs, floors=floors,
               nonzero_charge_fraction=float(np.sum(w*Qm)/Qt))
    return row
