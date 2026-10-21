"""Head start of the k=0 S mode over the lowest populated mode: audit's WKB-from-crossing number vs the actual
linear amplitude ratio (solver leapfrog, preparation initial velocity chi'=chi/40), and vs a sourced zero mode.
Usage: python3 -I head_start.py <codex_raw_dir> <out_json>"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import leapfrog_ql as lq

raw, outp = Path(sys.argv[1]), Path(sys.argv[2])
def wkb(c2, a):
    c = np.sqrt(c2); s = np.sqrt(max(a * a - c2, 0))
    return 40 * 0.5 * (a * s - c2 * np.log(max((a + s) / c, 1))) if a > c else 0.0
def lin_amp(k2, kap, a_end, dt=0.06):
    F2 = kap ** (2 / 3); x, p = 1.0, 1 / 40; steps = int(round((a_end - 1) * 40 / dt)); amp = []
    for it in range(1, steps + 1):
        a0, a1 = 1 + (it - 1) * dt / 40, 1 + it * dt / 40
        p += 0.5 * dt * (a0 * a0 - F2 - k2) * x; x += dt * p; p += 0.5 * dt * (a1 * a1 - F2 - k2) * x
    # envelope: |x| and sqrt(x^2 + (p/rate_local)^2) not needed after crossing; report |x|
    return abs(x)
res = {}
for run, (tab, N, dx, kap, a5) in {'main_s101': ('modes_L76.8_s101.txt', 192, .4, 2., 1.630), 'main_s202': ('modes_L76.8_s202.txt', 192, .4, 2., 1.578),
                                   'main_s303': ('modes_L76.8_s303.txt', 192, .4, 2., 1.617), 'kappa4_s101': ('modes_L76.8_s101.txt', 192, .4, 4., 1.929),
                                   'large_s101': ('modes_L115.2_s101.txt', 288, .4, 2., 1.635)}.items():
    t = np.loadtxt(raw / tab); n = t[:, :3]
    k2 = np.sum(4 * np.sin(np.pi * n / N) ** 2, axis=1) / dx ** 2
    order = np.argsort(k2); F2 = kap ** (2 / 3)
    r = {'a_rms05': a5, 'khat_min': float(np.sqrt(k2[order[0]]))}
    r['wkb_head_k0_vs_kmin'] = wkb(F2, a5) - wkb(F2 + k2[order[0]], a5)
    g0 = lin_amp(0.0, kap, a5); gk = [lin_amp(k2[i], kap, a5) for i in order[:6]]
    r['amp_head_k0_vs_kmin'] = float(np.log(g0 / gk[0]))
    r['amp_head_k0_vs_best_of_lowest6'] = float(np.log(g0 / max(gk)))
    r['lowest6_amp_efolds_rel_k0'] = [float(np.log(g / g0)) for g in gk]
    res[run] = r
    print(run, json.dumps({k: (round(v, 3) if isinstance(v, float) else [round(x, 2) for x in v] if isinstance(v, list) else v) for k, v in r.items()}))
outp.write_text(json.dumps(res, indent=1))
