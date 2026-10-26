"""Create the result-centered report and precise scientific comparison figure."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import ROOT


def main():
    d=json.loads((ROOT/'COMBINED_RESULTS.json').read_text());b=json.loads((ROOT/'BENCHMARK_RESULTS.json').read_text())
    v=json.loads((ROOT/'VECTOR_RESULTS.json').read_text())
    fig,ax=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    colors=['#176d88','#62a5aa','#c27830']
    for j,r in enumerate(d['rows']):
        x=0.
        for k,color,name in zip(['R_Y','R_Sigma','R_heavy_cascade_central'],colors,['Y daughters','Sigma radial → vectors','heavy-Y cascade → vector']):
            ax[0].barh(j,r[k],left=x,height=.46,color=color,label=name if j==0 else None);x+=r[k]
        raw=next(q for q in r['thresholds'] if q['ymult']!=1.)
        lim=sum(raw['R_L_1040_bracket'])/2
        ax[0].plot([lim,lim],[j-.3,j+.3],c='#ad1736',lw=2,label='case-L threshold, raw yield' if j==0 else None)
        ax[0].text(lim,j+.34,f'{lim:.1f}',ha='center',color='#ad1736',fontsize=9)
        ax[0].text(x+4,j,f'{x:.1f}',va='center',fontsize=9)
    ax[0].set(yticks=[0,1],yticklabels=['candidate 10\nphase 30°','candidate 17\nphase 37.5°'],
        xlabel='Decay-averaged injection / original X energy',title='Both lower ledgers exceed the lump threshold',xlim=(0,330))
    ax[0].legend(fontsize=8,loc='center right')
    for j,r in enumerate(d['rows']):
        raw=r['ymult_raw'];low=b[f"{r['candidate']}_{raw:.6g}_R_lower"]['f_extra']
        full=next(q for q in d['full_central_benchmark_relics'] if q['candidate']==r['candidate'] and q['ymult']!=1.)['f_extra']
        ax[1].bar(j-.18,low,width=.34,color='#176d88',label='omit all cascade/primordial vectors' if j==0 else None)
        ax[1].bar(j+.18,full,width=.34,color='#c27830',label='central complete stated ledger' if j==0 else None)
        ax[1].text(j-.18,low*1.12,f'{low:.2f}',ha='center',fontsize=9)
        ax[1].text(j+.18,full*1.12,f'{full:.2f}',ha='center',fontsize=9)
    ax[1].axhline(.01,c='#ad1736',ls='--',lw=1.3,label='one-percent allowance')
    ax[1].set(yscale='log',ylim=(.006,45),xticks=[0,1],xticklabels=['candidate 10','candidate 17'],
        ylabel='Symmetric mass / calibrated dark-matter mass',title='Conditional equilibrium relic; raw charge density')
    ax[1].legend(fontsize=8,loc='lower left')
    for a in ax:a.grid(alpha=.18);a.set_axisbelow(True)
    fig.suptitle('Fixed low-scale candidates: Y and Sigma debris close the supplied case-L benchmark',fontsize=12)
    fig.savefig(ROOT/'low_scale_debris.png',dpi=170);plt.close(fig)
    table=[]
    for r in d['rows']:
        raw=next(q for q in r['thresholds'] if q['ymult']!=1.)
        lim=sum(raw['R_L_1040_bracket'])/2
        fr=next(q['f_extra'] for q in d['full_central_benchmark_relics'] if q['candidate']==r['candidate'] and q['ymult']!=1.)
        table.append(f"| {r['candidate']} | {np.degrees(r['delta']):g}° | {r['epsilon_net']:.6f} | {r['ymult_raw']:.3f} | {r['R_Y']:.2f} | {r['R_Sigma']:.2f} | {r['R_total_central_including_primordial']:.2f} | {lim:.2f} | {fr:.2f} |")
    text=f'''# The two smooth low-scale candidates fail the debris-inclusive lump benchmark

Codex / AInstein, 25 October 2026. Parent 373d37e; independent review
[f1aa7df](https://github.com/ryangrvr/TestingGRUT/commit/f1aa7df968cb3dc90ca32f4f32aeb0a819ddfce4).

**Both candidates fail.** Their prepared releases remain smooth, but the
decay-averaged Y and Sigma injection already exceeds the recomputed case-L
threshold, even if every heavy-cascade and primordial vector is omitted.
The complete central ledgers are R≈228.4 and 280.5, compared with thresholds
about 103.0 and 93.4 at the actual raw charge densities.

This closes the two remaining published low-scale releases in the stated
decay/equilibrium benchmark. It is not a theorem for every continuous
initial phase, altered history or fully interacting quantum formation model.

![Debris and relic comparison](low_scale_debris.png)

## Fixed cases and numerical result

Both use c_H=10 and H_i=4.348677982066179e11 GeV, unchanged couplings,
the full aligned radial minimum, a supplied common Sigma phase, zero
physical velocities and instantaneous radiation release. Each uses the
published 48-member sigma=0.001 phase ensemble; no favorable nominal
particle count is substituted for its average.

R denotes injected conformal radiation energy divided by the original
X-condensate energy. R_L is the 1% crossing of the supplied case-L proxy.
The final f is symmetric mass divided by the calibrated dark-matter mass,
assuming a later normalization of the net charge yield.

| Candidate | Phase | Net epsilon | Raw yield / Y_Q | R_Y | R_Sigma | Central total R | R_L at raw density | Final f |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(table)}

The conservative Y+Sigma ledgers, without any cascade vectors, are
R=199.40 and 249.69. Their final f values are 10.85 and 13.94: about
1085 and 1394 times a one-percent allowance. The central complete ledgers
give f=12.08 and 15.27. All of these are **benchmark calculations**, not
new observational confidence limits.

At a yield already normalized before hidden freeze-out, the thresholds
are lower still: R_L≈91.2 and 83.1. Using the raw factors 5.69/5.31 improves
the reaction time but does not rescue either case. No real dilution field
or reheating history is implemented by this comparison; without it the
raw net abundance remains too large.

## What supplies the injection

All 96 published canonical endpoints are projected with an independent
real-matrix method. Four new full scalar releases check long windows and
relative tolerances 2e-12 and 2e-14. The ensemble means are approximately:

| Candidate | N_X | N_YH | N_YL | N_Sigma radial |
|---:|---:|---:|---:|---:|
| 10 | 2693.96 | 16.78 | 22.18 | 24.98 |
| 17 | 2513.73 | 17.81 | 24.45 | 36.65 |

Every Y quantum is counted once as three X or three anti-X daughters.
Every radial Sigma produces two vectors. A heavy-to-light Y cascade adds
one vector without adding another Y quantum. Each vector ultimately gives
one X/anti-X pair. Sigma's stored angular charge is not counted as radial
particle debris.

Exponential decay is averaged over its actual radiation-era time
distribution. The late vectors carry matter-like energy for a long period
before decay, which explains why their injection matters despite their
small coupling. Redshifting momentum and time-dilated decay are checked
against an independent energy-work integral.

The heavy population during the early mixing/decay era is not solved
exactly. Zero to the full Y population is therefore retained as a cascade
bracket. Both cases fail at its **zero-cascade lower end**. The total
unprocessed pair counts are about 0.0465 and 0.0615 per net charge under
the central ledger, versus 0.00493 allowed at 1%, before the much larger
cold-X remainder is added.

## Vector withdrawal and source scope

The earlier approximately 1060 inflationary vector pairs per net charge
is **withdrawn**, as Claude states in f1aa7df. It used the constant vacuum
mass. A new Proca control with the declared changing-mass profile gives
only about 2.3–2.6e-6 primordial pairs per net charge for these cases,
or R≈0.011. That source cannot be used to argue against the point. Vectors
from Sigma/Y decays are a separate population and dominate this ledger.
Release-induced Goldstone production is not computed here.

## Thermal calculation and checks

The Bose gas and reaction integration are independently reimplemented.
Fixed quadratures match adaptive integrals; energy, entropy and charge
constraints close explicitly. Steps 260/520/1040, quadratures 96/128,
late endpoints and tenfold changes of the supplied 4->2 rates are retained.
The threshold is bracketed directly at each epsilon and each yield.

The main solver uses a monotone coarse entropy-degree table and an
implicit positive exponential midpoint update. A separate legacy result
keeps the reviewer table/update. The reviewer table's g_s switch can make
inferred time nonmonotone under refinement; this is not used in production.
The effect is small compared with the failure margin.

The calculation extends below the review's 1-keV endpoint. Even assigning
all net charge to ideal cold conversion targets for all later times cannot
erase the excess: the conservative post-tail lower bounds remain near
10–13 times the entire calibrated dark-matter mass, far above 0.01.

## What is and is not closed

The *published prepared releases* fail the frozen and supplied equilibrium
benchmarks, including the two remaining smooth candidates. A complete
"create and protect" cosmological point at these couplings is unavailable.
Charge protection and cold object binding retain their narrower results.

This postprocesses late particle projections on undamped homogeneous
releases. Early decay backreaction, a full coupled scalar/vector quantum
state, thermalization kinetics, fragmentation, hot lump stability and
solitosynthesis have not been certified. A quantitatively different
interacting formation history would be a new calculation. No claim of a
universal exclusion over all phases or an observational exclusion is made.

[MATHEMATICS.md](MATHEMATICS.md), [CORRECTIONS.md](CORRECTIONS.md),
[SOURCES.md](SOURCES.md), [REPRODUCE.md](REPRODUCE.md) and
[VERIFICATION_RESULTS.json](VERIFICATION_RESULTS.json) give definitions,
provenance, controls and the separate physical statuses. All source and
result files are included in MANIFEST.json. Claude review of this packet
is pending; same-session AI work is not outside replication.
'''
    (ROOT/'REPORT.md').write_text(text)


if __name__=='__main__':main()
