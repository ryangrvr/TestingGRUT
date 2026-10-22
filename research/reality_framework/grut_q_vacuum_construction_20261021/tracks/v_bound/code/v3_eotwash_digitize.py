"""v_bound: independent re-digitization of the Eot-Wash 2020 alpha-lambda exclusion curve (Lee et al., PRL 124, 101101
(2020), arXiv:2002.11761, Fig. 5 bottom, 95% CL) at |alpha| = 1 and |alpha| = 1/3.

Method (different from the bound track's): the arXiv source figure fig5b1.pdf (downloaded to the scratchpad, NOT in the
repository) is converted to SVG with `pdftocairo -svg`; the plot frame is calibrated from the major tick marks (x: 1e-5,
1e-4, 1e-3 m; y: 1e-3 ... 1e6, ten ticks, least-squares fit); every stroked path segment is intersected with the
horizontal lines alpha = 1 and alpha = 1/3, and the crossings are grouped by stroke colour and width.
Usage:  python3 -I v3_eotwash_digitize.py /path/to/fig5b1.svg  > ../results/v3_eotwash_digitize.json
(The SVG path is an argument; no figure data are stored in the repository.)
"""
import json
import re
import sys
import numpy as np


def parse(svg_text):
    recs = []
    for body in re.findall(r'<path([^>]*)/>', svg_text):
        st = re.search(r'stroke="([^"]*)"', body)
        sw = re.search(r'stroke-width="([^"]*)"', body)
        d = re.search(r' d="([^"]*)"', body)
        if not (st and d):
            continue
        toks = re.findall(r'[MLCZ]|-?[\d.]+', d.group(1))
        pts, cmd, i = [], None, 0
        while i < len(toks):
            tk = toks[i]
            if tk in 'MLCZ':
                cmd = tk
                i += 1
                continue
            if cmd in ('M', 'L'):
                pts.append((float(toks[i]), float(toks[i + 1])))
                i += 2
            elif cmd == 'C':
                pts.append((float(toks[i + 4]), float(toks[i + 5])))
                i += 6
            else:
                i += 1
        recs.append((st.group(1), sw.group(1) if sw else None, pts))
    return recs


def main(path):
    recs = parse(open(path).read())
    black = [r for r in recs if r[0] == 'rgb(0%, 0%, 0%)' and len(r[2]) == 2]
    # major ticks: x ticks of length ~110 on the bottom axis; y ticks of length ~110 on the left axis
    xt = [r[2] for r in black if abs(r[2][0][0] - r[2][1][0]) < 1e-6 and 100 < abs(r[2][1][1] - r[2][0][1]) < 120]
    from collections import Counter
    ybot = min(y for y, n in Counter(round(min(p[0][1], p[1][1]), 1) for p in xt).items() if n >= 3)   # bottom axis
    xmaj = sorted({round(p[0][0], 3) for p in xt if abs(min(p[0][1], p[1][1]) - ybot) < 1})
    yt = [r[2] for r in black if abs(r[2][0][1] - r[2][1][1]) < 1e-6 and 100 < abs(r[2][1][0] - r[2][0][0]) < 120]
    xleft = min(x for x, n in Counter(round(min(p[0][0], p[1][0]), 1) for p in yt).items() if n >= 3)  # left axis
    ymaj = sorted({round(p[0][1], 3) for p in yt if abs(min(p[0][0], p[1][0]) - xleft) < 1})
    xv = np.array([-5.0, -4.0, -3.0])
    yv = np.arange(-3.0, 7.0)
    ax, bx = np.polyfit(xv, np.array(xmaj), 1)
    ay, by = np.polyfit(yv, np.array(ymaj), 1)
    lam_of = lambda X: 10**((X - bx) / ax)
    Y_of = lambda alpha: ay * np.log10(alpha) + by
    cross = {}
    for st, sw, pts in recs:
        for alpha in (1.0, 1.0 / 3.0):
            yc = Y_of(alpha)
            for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
                if (y1 - yc) * (y2 - yc) < 0:
                    xc = x1 + (yc - y1) / (y2 - y1) * (x2 - x1)
                    if xmaj[0] - 900 < xc < xmaj[-1] + 10:
                        cross.setdefault('%s | width %s' % (st, sw), []).append({'alpha': alpha, 'lambda_um': float(lam_of(xc) * 1e6)})
    purple = cross.get('rgb(33.300781%, 0%, 100%) | width 13.0392', [])
    a1 = sorted(c['lambda_um'] for c in purple if c['alpha'] == 1.0)
    a3 = sorted(c['lambda_um'] for c in purple if abs(c['alpha'] - 1 / 3) < 1e-9)
    out = {'source': 'Lee et al. 2020, arXiv:2002.11761, fig5b1.pdf (scratchpad copy), via pdftocairo -svg',
           'frame calibration': {'x major ticks (svg units)': xmaj, 'y major ticks (svg units)': ymaj,
                                 'svg units per decade x, y': [ax, ay]},
           'thick purple curves (Eot-Wash 2020 / 2007 / 2004), crossings in um': {'alpha=1': a1, 'alpha=1/3': a3},
           'all crossings by stroke': cross,
           'reading': ('Eot-Wash 2020 = leftmost thick purple curve (left edge of the green band); Eot-Wash 2007 = second; '
                       'the third is Eot-Wash 2004.')}
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
