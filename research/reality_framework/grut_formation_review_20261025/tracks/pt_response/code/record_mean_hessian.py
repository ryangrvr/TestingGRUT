"""Background recorder.  ENGINE = Codex's fieldstep (copied 5f5513b tree, run from $PKC/code); DIAGNOSTICS = ours.

From the regenerated base48 state at s = 20 (bitwise equal to Codex's history), advance the nonlinear field to
s = 26.8 with the same step rule (hmax = 0.002) and record every Delta s = 0.01:
  * coherent zero mode (q1,q2,qphi,p1,p2,pphi),
  * spatial means needed for the spatially averaged Hessian Hbar_ab(s) of V (3x3: chi1, chi2, phi):
      Hbar_ab = <(lam r^2 + g2 phi^2) delta_ab + 2 lam chi_a chi_b>,  Hbar_a,phi = <2 g2 phi chi_a>,
      Hbar_phi,phi = <3 lamphi phi^2 + g2 r^2>
  * the same Hessian built from the coherent mode alone (bare condensate model),
  * Codex's Hartree_m2 = 2 lam <r^2> + g2 <phi^2>.
Output: $SP/review/pt_response/hessian_s20_26.8.npz
"""
import sys
sys.path.insert(0, '.')
from pt_common import *  # noqa  (Codex module from the copy)

LX_, G2_, LP_ = .001, .01, .03
N, L = 48, 12.
z = np.load(WORK / 'base48' / 'state_s20.npy'); z = np.ascontiguousarray(z)
set_box(L); lib = library()
s = 20.; ds = .01; nsamp = int(round(6.8 / ds)) + 1
out = dict(s=np.zeros(nsamp), mean=np.zeros((nsamp, 6)), Hbar=np.zeros((nsamp, 3, 3)), Hcoh=np.zeros((nsamp, 3, 3)),
           hf=np.zeros(nsamp), r2=np.zeros(nsamp), cov=np.zeros((nsamp, 3, 3)), rho=np.zeros(nsamp), Q=np.zeros(nsamp))

def record(i, s, z):
    c1, c2, ph = z[0], z[1], z[2]
    r2 = c1 * c1 + c2 * c2
    m = z.mean(axis=(1, 2, 3))
    H = np.zeros((3, 3))
    base = LX_ * r2.mean() + G2_ * (ph * ph).mean()
    H[0, 0] = base + 2 * LX_ * (c1 * c1).mean(); H[1, 1] = base + 2 * LX_ * (c2 * c2).mean()
    H[0, 1] = H[1, 0] = 2 * LX_ * (c1 * c2).mean()
    H[0, 2] = H[2, 0] = 2 * G2_ * (ph * c1).mean(); H[1, 2] = H[2, 1] = 2 * G2_ * (ph * c2).mean()
    H[2, 2] = 3 * LP_ * (ph * ph).mean() + G2_ * r2.mean()
    a, b, p = m[0], m[1], m[2]; rc2 = a * a + b * b
    Hc = np.array([[LX_ * rc2 + G2_ * p * p + 2 * LX_ * a * a, 2 * LX_ * a * b, 2 * G2_ * p * a],
                   [2 * LX_ * a * b, LX_ * rc2 + G2_ * p * p + 2 * LX_ * b * b, 2 * G2_ * p * b],
                   [2 * G2_ * p * a, 2 * G2_ * p * b, 3 * LP_ * p * p + G2_ * rc2]])
    d = z[:3] - m[:3, None, None, None]
    cov = np.einsum('aijk,bijk->ab', d, d) / d[0].size
    out['s'][i] = s; out['mean'][i] = m; out['Hbar'][i] = H; out['Hcoh'][i] = Hc
    out['hf'][i] = 2 * LX_ * r2.mean() + G2_ * (ph * ph).mean(); out['r2'][i] = r2.mean(); out['cov'][i] = cov
    out['Q'][i] = (c1 * z[4] - c2 * z[3]).mean()

record(0, s, z); prev = 0.
for i in range(1, nsamp):
    # same rule as Codex: call step for each sample interval with hmax 0.002 (here 5 steps of 0.002)
    step(lib, z, ds, .002, threads=2)
    record(i, 20. + i * ds, z)
    if i % 100 == 0: print(i, flush=True)
np.savez_compressed(Path(sys.argv[1]), **out)
print('saved', sys.argv[1])
