"""REPRODUCTION driver (uses Codex's engine/initialiser from a COPY of the 5f5513b tree; run from $PKC/code).

Regenerates the base48 (N48, L12, seed 1, original band) nonlinear field at s = 20, 20.24, 20.44, 20.68 with exactly
pt_production's sampling/step rule (Delta s = 0.04 samples, hmax = 0.002), and checks the trace against the
published FIELD_TRACES.npz (base48) at every 25th sample.  Output: states in $GRUT_PARTICLE_WORKDIR/base48/.
"""
import sys, json, time
sys.path.insert(0, '.')
from pt_common import *            # noqa  (Codex module, copied tree)
from prepare import initialize     # noqa
from pt_production import observe, NAMES  # noqa

pub = sys.argv[1]  # published FIELD_TRACES.npz
out = WORK / 'base48'; out.mkdir(exist_ok=True)
N, L = 48, 12.
z, params = initialize(N, L, 1, None); z = np.ascontiguousarray(z)
set_box(L); lib = library()
ds = .04; samples = np.arange(851) * ds
kc2, kh2 = geometry(N, L, half=True); nbins = 96; edges = np.linspace(0, 24., nbins + 1)
bins = np.clip(np.digitize(np.sqrt(kc2).ravel(), edges) - 1, 0, nbins - 1)
weights = np.full(N // 2 + 1, 2.); weights[[0, -1]] = 1.; weights = np.broadcast_to(weights[None, None, :], kh2.shape).copy()
P = np.load(pub); pv = P['base48_values']
snap = {20., 20.24, 20.44, 20.68}; prev = 0.; worst = 0.; t0 = time.monotonic(); checks = []
for i, s in enumerate(samples):
    if s > 20.68 + 1e-9: break
    if i: step(lib, z, s - prev, .002, threads=2)
    prev = s
    if i % 25 == 0 or any(abs(s - x) < 1e-9 for x in snap):
        _, row, *_ = observe(z, L, kh2, weights, bins, nbins)
        rel = np.max(abs(row - pv[i]) / (1e-12 + abs(pv[i])) * (abs(pv[i]) > 1e-8))
        worst = max(worst, float(rel)); checks.append((float(s), float(rel)))
    if any(abs(s - x) < 1e-9 for x in snap):
        np.save(out / f'state_s{s:g}.npy', z); print('saved', s, 'worst rel trace diff so far', worst, round(time.monotonic() - t0, 1), flush=True)
json.dump(dict(worst_relative_trace_difference=worst, checks=checks, params={k: (v if not isinstance(v, np.generic) else v.item()) for k, v in params.items() if isinstance(v, (int, float, str, np.generic))}),
          open(out / 'regen_check.json', 'w'), indent=1)
print('done', worst)
