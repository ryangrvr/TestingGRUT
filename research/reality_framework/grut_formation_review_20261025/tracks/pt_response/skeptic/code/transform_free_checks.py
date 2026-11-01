"""SKEPTIC (pt_response_skeptic), INDEPENDENT code, no Codex code.
(1) Canonical-sum identity.  Codex's finite transform (delays 0..6 step .04, cos^2 taper, trapezoid, omega grid
    +-20 step .02, rho = -2 Im G(omega), Z = int w rho dw/2pi) applied to backgrounds with KNOWN response: free modes
    G = sin(w t)/w, two-line mixtures, damped lines, and a purely 'anomalous-converting' toy
    G_N = [sin(w1 t)/w1 + sin(w2 t)/w2]/2.  If Z_+ + Z_- - 1 has the same size as Codex's (<= 0.4%) for all of them,
    the sum tests only the transform.  Also the published normal arrays: Z_sum vs the discrete estimate of G_N'(0).
(2) Free control.  The exact discrete Yoshida-4 map (kick-drift-kick, w1 = 1/(2-2^{1/3})) for one FD4 mode
    (N=8, L=5, k = 2pi/5 along x), 10 substeps of h = .004 per sample, as a 2x2 matrix product in 50-digit mpmath;
    max |q(t) - sin(w t)/w| over 151 samples.
Usage: python3 -I transform_free_checks.py <packet dir> <out json>
"""
import sys, json
import numpy as np
import mpmath as mp
P = sys.argv[1].rstrip('/') + '/'; OUT = sys.argv[2]
t = np.arange(151) * .04; w = np.linspace(-20, 20, 2001); dw = w[1] - w[0]
taper = np.cos(np.pi * t / 12) ** 2; trap = np.full(151, .04); trap[[0, -1]] *= .5
TR = np.exp(-1j * w[:, None] * t[None, :])
def Zsum(g):
    rho = -2 * (TR @ (g * trap * taper)).imag
    return float(np.sum(w * rho) * dw / (2 * np.pi))
res = {}
free = {f'{om:g}': Zsum(np.sin(om * t) / om) - 1 for om in (1.25, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0)}
mix = {f'{a:g}_{b:g}': Zsum(.5 * (np.sin(a * t) / a + np.sin(b * t) / b)) - 1 for a, b in ((3.6, 4.4), (3.4, 8.6), (4.2, 5.0))}
damp = {f'{om:g}_G{G:g}': Zsum(np.exp(-G * t / 2) * np.sin(om * t) / om) - 1 for om in (4.0, 8.0) for G in (.3, 1.)}
cplx = {f'{om:g}': Zsum(np.exp(1j * 3.4 * t) * np.sin(om * t) / om) - 1 for om in (4.0, 8.0)}  # rotating-frame line
res['canonical_sum_minus_1'] = dict(free=free, two_lines=mix, damped=damp, complex_rotated=cplx)
pub = {}
for tag in ('s20_seed7', 's20.44_seed7', 's20_seed11'):
    d = np.load(P + f'RESPONSE_{tag}.npz'); nr = d['normal']
    zs = np.array([Zsum(nr[:, k]) for k in range(nr.shape[1])])
    slope = (nr[1] / .04).real
    pub[tag] = dict(Zsum_minus_1=(zs - 1).tolist(), initial_slope=slope.tolist(), max_abs=float(np.max(abs(zs - 1))))
res['published'] = pub
allf = [abs(v) for g in res['canonical_sum_minus_1'].values() for v in g.values()]
res['known_response_max_abs_Zsum_minus_1'] = max(allf)
print('known-response |Zsum-1| max %.4g' % max(allf)); print(json.dumps(res['canonical_sum_minus_1'], indent=1))
for k, v in pub.items(): print(k, 'max |Z-1| %.4g' % v['max_abs'])
# (2) exact discrete Yoshida map
mp.mp.dps = 50
N, L = 8, mp.mpf(5); dx = L / N; u = 2 * mp.pi / N
om2 = 4 * mp.sin(u / 2) ** 2 * (1 + mp.sin(u / 2) ** 2 / 3) / dx ** 2; om = mp.sqrt(om2)
w1 = 1 / (2 - mp.cbrt(2)); w0 = 1 - 2 * w1; h = mp.mpf('0.004')
def K(c): return mp.matrix([[1, 0], [-om2 * c * h, 1]])
def D(c): return mp.matrix([[1, c * h], [0, 1]])
step = K(w1 / 2) * D(w1) * K((w1 + w0) / 2) * D(w0) * K((w0 + w1) / 2) * D(w1) * K(w1 / 2)
step = step.T if False else step
# state column (q, p); apply right-to-left: first kick is the rightmost
S10 = step ** 10
y = mp.matrix([[0], [1]]); errs = []
for j in range(151):
    if j: y = S10 * y
    errs.append(abs(y[0] - mp.sin(om * mp.mpf(j) * mp.mpf('0.04')) / om))
res['free_control'] = dict(omega=float(om), max_abs_error_exact_discrete_map=float(max(errs)),
                           codex=2.161931467181688e-10, pt_response=2.161928e-10)
print('free control: omega %.16f  max error %.10e (Codex 2.161931e-10, pt_response 2.161928e-10)' % (float(om), float(max(errs))))
json.dump(res, open(OUT, 'w'), indent=1)
