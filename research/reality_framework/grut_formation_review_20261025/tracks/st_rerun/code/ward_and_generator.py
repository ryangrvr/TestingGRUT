"""st_rerun audit helper (uses Codex's Python matrices; reproduction-type check, not independent physics).

(1) Ward identity behind the C projection: M g == ward, with ward = sqrt2 a^3 T F + a^2((xi on X) - 1/6) R g,
    evaluated on every published background sample of quiet3/quiet5 (tau 0..800).
(2) Direct-vector RHS == full generator column: compare Codex's Python generator(v,beta,k) @ x with the C formula
    written out in Python (Ward-projected), for random x, on the same samples.
(3) Gauss-law/canonical-momentum consistency of the published selected vectors: the canonical p and q' differ by
    (e^2/k^2) g (g.p - g'.q); report |g.p-g'.q| scale and the X-only part used as the local charge.
Usage: python3 -B ward_and_generator.py <codex_code_dir> <published_dir> <out_json>
"""
import sys, json
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from timing_common import M, np  # noqa
from reheating import geometric, XI, ETA2  # pinned 2cbc7f1/e160356 background helpers
pub = Path(sys.argv[2]); out = {}
TCH = M.TCH
rng = np.random.default_rng(1)

def c_style_rhs(v, beta, k, x, e=1e-7):
    Mm, g, gp, chi, chip = M.matrices(v, beta)
    z, vel, V, F, f, fd, R, hd = geometric(v, beta, XI, ETA2); a = np.exp(v[15])
    ward = np.empty(6)
    for i, q in enumerate([1, 3, 6]):
        extra = a * a * ((XI if i == 0 else 0) - 1 / 6) * R
        ward[2 * i] = -np.sqrt(2) * q * a ** 3 * F[i].imag + extra * g[2 * i]
        ward[2 * i + 1] = np.sqrt(2) * q * a ** 3 * F[i].real + extra * g[2 * i + 1]
    c = e * e / k ** 2; T = x
    current = g @ T[6:] - gp @ T[:6]; common = g @ T[:6] / (g @ g)
    mq = common * ward + Mm @ (T[:6] - common * g)
    d = np.empty(12); d[:6] = (T[6:] + c * g * current) / a; d[6:] = (-k * k * T[:6] - mq + c * gp * current) / a
    return d, Mm, g, ward

for tag, beta in [('quiet3', 3.), ('quiet5', 5.)]:
    z = np.load(pub / (tag + '.npz')); bg = z['bg']; times = z['times']
    wres = []; gres = []
    for n in range(0, len(bg), 7):
        v = bg[n]
        for k in (.02, 2.13, 5.8):
            x = rng.standard_normal(12)
            d, Mm, g, ward = c_style_rhs(v, beta, k, x)
            G = M.generator(v, beta, k)
            ref = G @ x
            gres.append(float(np.max(np.abs(d - ref)) / max(1e-300, np.max(np.abs(G)) * np.max(np.abs(x)))))
        wres.append(float(np.max(np.abs(Mm @ g - ward)) / (np.max(np.abs(Mm)) * np.max(np.abs(g)))))
    out[tag] = dict(samples=len(wres), ward_identity_max_relative=max(wres), ward_identity_at_tau800=wres[-1],
                    projected_vs_generator_max_relative=max(gres))
    print(tag, out[tag], flush=True)

# (3) canonical-momentum / velocity relation on the archived selected vectors (quiet3, k=2.13)
z = np.load(pub / 'quiet3.npz'); bg = z['bg']; sv = z['selected_vectors']; sk = z['selected_vector_ks']
rows = []
for n in [0, 200, 600, 1000, 1582]:
    v = bg[n]; Mm, g, gp, chi, chip = M.matrices(v, 3.)
    for j, k in enumerate(sk):
        x = sv[n, j]; cur = g @ x[6:] - gp @ x[:6]
        xcur = (TCH[:2, :2] @ chi[:2]) @ x[6:8] - (TCH[:2, :2] @ chip[:2]) @ x[:2]
        dN = x[0] * chip[1] - x[1] * chip[0] + chi[0] * x[7] - chi[1] * x[6]
        rows.append(dict(tau=float(z['times'][n]), k=float(k), total_current=float(cur), X_part_gp=float(xcur), responses_dN=float(dN),
                         gauge_velocity_shift_rel=float(1e-14 / k ** 2 * abs(cur) * np.linalg.norm(g) / max(1e-300, np.linalg.norm(x[6:])))))
out['local_charge_identity'] = dict(max_abs_diff_dN_vs_gX_dot_dp_minus_gpX_dot_dq=max(abs(r['responses_dN'] - r['X_part_gp']) / max(1., abs(r['X_part_gp'])) for r in rows),
                                    rows=rows)
print('local charge identity', out['local_charge_identity']['max_abs_diff_dN_vs_gX_dot_dp_minus_gpX_dot_dq'])
for r in rows[:6]: print(r)
Path(sys.argv[3]).write_text(json.dumps(out, indent=1) + '\n')
