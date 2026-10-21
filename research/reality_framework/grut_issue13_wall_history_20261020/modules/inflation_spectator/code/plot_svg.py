"""Hand-written SVG of p_-(H_inf) (no plotting library; numpy/scipy-only rule). Reads results/stochastic_1d.json
and results/final_sign.json; writes results/p_minus_vs_Hinf.svg."""
import json
import os
import numpy as np
from common import RESULTS, PC_SITE_SC

W, Hh = 780, 470
ML, MR, MT, MB = 74, 24, 70, 56
X0, X1 = 2.0, 13.0
Y0, Y1 = -40.0, 0.5


def sx(x):
    return ML + (x - X0) / (X1 - X0) * (W - ML - MR)


def sy(y):
    return MT + (Y1 - y) / (Y1 - Y0) * (Hh - MT - MB)


def path(xs, ys):
    pts = [(sx(x), sy(max(y, Y0))) for x, y in zip(xs, ys) if np.isfinite(y)]
    return "M" + " L".join("%.1f,%.1f" % p for p in pts)


def main():
    d = json.load(open(os.path.join(RESULTS, "stochastic_1d.json")))
    fs_path = os.path.join(RESULTS, "final_sign.json")
    fs = json.load(open(fs_path)) if os.path.exists(fs_path) else None
    rows = d["rows"]
    lx = np.log10([r["H_GeV"] for r in rows])
    series = [("N<tspan baseline-shift='sub' font-size='10'>*</tspan> = 50, fresh start at +v", "s1",
               [r["N50"]["delta_start"]["log10_p_minus"] for r in rows], "N*=50 fresh"),
              ("N<tspan baseline-shift='sub' font-size='10'>*</tspan> = 60, fresh start at +v", "s2",
               [r["N60"]["delta_start"]["log10_p_minus"] for r in rows], "N*=60 fresh"),
              ("N<tspan baseline-shift='sub' font-size='10'>*</tspan> = 50, relaxed + well", "s3",
               [r["N50"]["relaxed_start"]["log10_p_minus"] for r in rows], "N*=50 relaxed")]
    o = []
    o.append("<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 %d %d' width='%d' height='%d' role='img' "
             "aria-labelledby='t d'>" % (W, Hh, W, Hh))
    o.append("<title id='t'>Opposite-vacuum volume fraction versus inflationary Hubble rate</title>")
    o.append("<desc id='d'>Stochastic spectator phi of unchanged GRUT Candidate 1 (X=D=0 path, xi run to mu=H). "
             "Inflation is outside-action preparation data. Data: results/stochastic_1d.json.</desc>")
    o.append("""<style>
svg{--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--grid:#e4e3df;--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--ref:#8a8984}
@media (prefers-color-scheme: dark){svg{--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--grid:#34332f;--s1:#3987e5;--s2:#d95926;--s3:#199e70;--ref:#8f8e86}}
text{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;fill:var(--ink2);font-size:12px}
.ttl{fill:var(--ink);font-size:15px;font-weight:600}
.ln{fill:none;stroke-width:2;stroke-linejoin:round;stroke-linecap:round}
.ref{stroke:var(--ref);stroke-width:1}
.grid{stroke:var(--grid);stroke-width:1}
.hit{fill:transparent}
.hit:hover{fill:var(--ink2);fill-opacity:.35}
</style>""")
    o.append("<rect x='0' y='0' width='%d' height='%d' fill='var(--surface)'/>" % (W, Hh))
    o.append("<text class='ttl' x='%d' y='24'>Opposite-vacuum fraction p<tspan baseline-shift='sub' font-size='11'>&#8722;</tspan> "
             "after the last N<tspan baseline-shift='sub' font-size='11'>*</tspan> e-folds</text>" % ML)
    # gridlines and ticks
    for k in range(2, 14):
        x = sx(k)
        o.append("<line class='grid' x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f'/>" % (x, sy(Y1), x, sy(Y0)))
        o.append("<text x='%.1f' y='%.1f' text-anchor='middle'>10<tspan baseline-shift='super' font-size='9'>%d</tspan></text>"
                 % (x, sy(Y0) + 18, k))
    for yv in (0, -10, -20, -30, -40):
        y = sy(yv)
        o.append("<line class='grid' x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f'/>" % (sx(X0), y, sx(X1), y))
        lab = "1" if yv == 0 else "10<tspan baseline-shift='super' font-size='9'>%d</tspan>" % yv
        o.append("<text x='%.1f' y='%.1f' text-anchor='end'>%s</text>" % (sx(X0) - 8, y + 4, lab))
    o.append("<text x='%.1f' y='%d' text-anchor='middle'>H<tspan baseline-shift='sub' font-size='9'>inf</tspan> (GeV)</text>"
             % ((sx(X0) + sx(X1)) / 2, Hh - 12))
    o.append("<text transform='translate(16,%.1f) rotate(-90)' text-anchor='middle'>p<tspan baseline-shift='sub' font-size='9'>&#8722;</tspan> (log scale)</text>"
             % ((sy(Y0) + sy(Y1)) / 2))
    # reference lines
    # labels sit where no curve passes: p_c at the far left (curves are at the bottom there),
    # 1e-3 and 1e-10 at H ~ 1e5-1e7 GeV (curves are near the top there)
    for val, lab, xl in ((PC_SITE_SC, "p_c = 0.3116 (3D site percolation)", X0 + 0.06),
                         (1e-3, "10<tspan baseline-shift='super' font-size='9'>&#8722;3</tspan>", 5.0),
                         (1e-10, "10<tspan baseline-shift='super' font-size='9'>&#8722;10</tspan>", 5.0)):
        y = sy(np.log10(val))
        o.append("<line class='ref' x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f'/>" % (sx(X0), y, sx(X1), y))
        o.append("<text x='%.1f' y='%.1f'>%s</text>" % (sx(xl), y - 5, lab))
    o.append("<text x='%.1f' y='%.1f'>curves continue below 10<tspan baseline-shift='super' font-size='9'>&#8722;40</tspan> "
             "for H<tspan baseline-shift='sub' font-size='9'>inf</tspan> &lt; 165 GeV</text>" % (sx(2.6), sy(-37.5)))
    # series
    for name, cls, ys, short in series:
        o.append("<path class='ln' style='stroke:var(--%s)' d='%s'/>" % (cls, path(lx, ys)))
    if fs:
        xs = [np.log10(r["H_GeV"]) for r in fs["rows"]]
        ys = [np.log10(max(r["N50"]["p_minus_after_trapping"], 1e-300)) for r in fs["rows"]]
        for x, y in zip(xs, ys):
            o.append("<circle cx='%.1f' cy='%.1f' r='4' fill='var(--s1)' stroke='var(--surface)' stroke-width='2'>"
                     "<title>H_inf = %.3g GeV: p_- after post-inflation trapping (N*=50) = %.3g</title></circle>"
                     % (sx(x), sy(y), 10 ** x, 10 ** y))
    # hover targets
    for name, cls, ys, short in series:
        for x, y in zip(lx, ys):
            if np.isfinite(y) and y >= Y0:
                o.append("<circle class='hit' cx='%.1f' cy='%.1f' r='6'><title>%s: H_inf = %.4g GeV, p_- = %.3g</title></circle>"
                         % (sx(x), sy(y), short, 10 ** x, 10 ** y))
    # legend
    lxp = ML
    for name, cls, ys, short in series:
        o.append("<line x1='%d' y1='46' x2='%d' y2='46' class='ln' style='stroke:var(--%s)'/>" % (lxp, lxp + 22, cls))
        o.append("<text x='%d' y='50'>%s</text>" % (lxp + 28, name))
        lxp += 225
    if fs:
        o.append("<circle cx='%d' cy='%d' r='4' fill='var(--s1)'/>" % (sx(6.6), sy(-33)))
        o.append("<text x='%.1f' y='%.1f'>dots: N* = 50 after post-inflation trapping</text>" % (sx(6.6) + 10, sy(-33) + 4))
    o.append("</svg>")
    with open(os.path.join(RESULTS, "p_minus_vs_Hinf.svg"), "w") as fh:
        fh.write("\n".join(o))


if __name__ == "__main__":
    main()
