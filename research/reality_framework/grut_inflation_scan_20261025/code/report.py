"""Make the scientific figure and result-centered review handoff report."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import ROOT


def main():
    f=json.loads((ROOT/'FINDINGS.json').read_text());e=json.loads((ROOT/'ENSEMBLE_SUMMARY.json').read_text())
    y=json.loads((ROOT/'Y_SOURCE_RESULTS.json').read_text());g=json.loads((ROOT/'GRID_RESULTS.json').read_text())
    lc=list(csv.DictReader((ROOT/'LOCAL_WINDOW_RESULTS.csv').open()))
    fig,axs=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
    for c,color in [(1.,'#126e87'),(10.,'#b96a25')]:
        gg=[x for x in g['groups'] if x['cH']==c]
        for x in gg:
            hh=x['H_GeV']
            for p in x['phases']:
                axs[0,0].scatter(hh,abs(p['epsX']),s=13,c=color if p['provisional_regular'] else '#c8cdd1',alpha=.8)
        axs[0,0].plot([],[],marker='o',ls='',c=color,label=f'cH={c:g}; provisional regular')
    axs[0,0].scatter(f['best_supplied_phase_H_GeV'],.240661234,s=90,marker='*',c='#ad1d39',label='refined supplied-phase point',zorder=5)
    axs[0,0].axhline(.19,c='#777',ls='--',lw=1,label='optimistic X-only benchmark')
    axs[0,0].set(xscale='log',xlabel='Inflation scale H_i [GeV]',ylabel='Absolute X condensate fraction',ylim=(0,.75),title='Coarse grid and the missed refined point')
    axs[0,0].legend(fontsize=7,loc='upper left')
    hs=sorted(set(float(x['H_GeV']) for x in lc));ds=sorted(set(float(x['delta']) for x in lc))
    z=np.array([[abs(float(next(r['epsX'] for r in lc if float(r['H_GeV'])==h and float(r['delta'])==d))) for h in hs] for d in ds])
    mesh=axs[0,1].pcolormesh(np.array(hs)/1e13,np.array(ds)-np.pi/2,z,shading='nearest',cmap='viridis',vmin=0,vmax=.35)
    axs[0,1].contour(np.array(hs)/1e13,np.array(ds)-np.pi/2,z,levels=[.19],colors=['white'],linewidths=.9)
    axs[0,1].set(xlabel='H_i / 10^13 GeV',ylabel='Source-phase offset from pi/2 [rad]',title='Local supplied-phase map; white contour = 0.19')
    fig.colorbar(mesh,ax=axs[0,1],label='Absolute X-only charge fraction',shrink=.8)
    b=e['initial_ensembles'][0];a=e['candidate_ensembles'][0];aa=e['high_fraction_controls'][0]['larger_seed_ensemble']
    vals=[b['eps_net'],-a['eps_net'],-aa['eps_net']];errs=[b['eps_net_sem'],a['eps_net_sem'],aa['eps_net_sem']]
    axs[1,0].errorbar(range(3),vals,yerr=errs,fmt='o',capsize=5,c='#126e87')
    axs[1,0].axhline(.19,c='#777',ls='--',lw=1)
    axs[1,0].set(xticks=range(3),xticklabels=['old high scale\n3 phase seeds 0.001','new point\nsource seeds 0.001','new point\nsource seeds 0.03'],ylabel='Mean charge fraction\n(new-point sign reversed)',title='Homogeneous ensembles: mean ± sampling SEM',ylim=(-.045,.29))
    labels=['1% allowance','H=1.12e13, cH=1','H=1.58e12, cH=10','H=8.03e12, cH=1']
    vals=[y['pair_limit']]+[r['extra_Y_decay_pairs_per_abs_X_charge'] for r in y['rows']]
    axs[1,1].barh(labels,vals,color=['#777','#ad1d39','#b96a25','#126e87'])
    axs[1,1].set(xscale='log',xlabel='Y-decay extra pairs / absolute net charge',title='Conditional unthermalized Y-source budget')
    for i,v in enumerate(vals):axs[1,1].text(v*1.04,i,f'{v:.4g}',va='center',fontsize=8)
    axs[1,1].set_xlim(.003,.15)
    for ax in axs.ravel():ax.grid(alpha=.18);ax.set_axisbelow(True)
    fig.suptitle('Fixed couplings: a conditional charged release; preparation and relics remain open',fontsize=13)
    fig.savefig(ROOT/'inflation_scan.png',dpi=170);plt.close(fig)
    rows=[]
    for i in range(3):
        a=e['candidate_ensembles'][i];q=e['high_fraction_controls'][i]['larger_seed_ensemble'];p=y['rows'][i]
        rows.append(f"| {a['cH']:g} | {a['H_GeV']:.6g} | {a['delta']:.6f} | {a['nominal_epsX']:.6f} | {a['eps_net']:.6f} ± {a['eps_net_sem']:.6f} | {q['eps_net']:.6f} ± {q['eps_net_sem']:.6f} | {p['extra_Y_decay_pairs_per_abs_X_charge']:.5f} |")
    report=f'''# Inflation-scale scan: a conditional charged release survives the X-only test

Codex / AInstein, 25 October 2026. Parent b5174e5. Review applied:
[Claude 2793dab](https://github.com/ryangrvr/TestingGRUT/commit/2793dab800794e78ebd74c43a512eb15622a7532).

**The scan does not support a blanket “no inflation scale works” conclusion.**
At the unchanged couplings, H_i=1.1211167304e13 GeV, c_H=1 and an externally
supplied common source phase delta=pi/2 give a reproducible X charge
fraction 0.24066 and a yield about 66 times the required Y_Q. A small-noise
ensemble preserves it. This clears the proposed X-only necessary tests.
It does **not** supply the phase preparation, the calibrated final abundance,
or an acceptable full leftover population.

The same point has a nonzero Y excitation. If its late number survives to
perturbative Y decay and no number changing follows, it supplies about
**0.0834 extra pairs per net charge**, versus a one-percent allowance
**0.00493**. That conditional source alone is about **16.9 times** the
allowance. This is not a solved thermal relic or a proof that every possible
hidden-sector evolution fails.

![Numerical scan](inflation_scan.png)

## What was actually run

- {f['production_release_count']} production releases: 572 coarse runs
  (13 H values × 2 c_H values × 11 phases × 2 phase-shift choices),
  1,426 denser H/phase runs, 60 follow-up phase shifts, 205 local-map runs,
  2,048 phase-ensemble members and 15 tighter/longer/preparation controls.
  Eight additional mathematical-control releases are counted separately.
- c_H=1 and 10. Coarse H spans 3e11/sqrt(c_H) through 4.7e13 GeV;
  refinement covers the transition regions. Complete configurations are
  in inputs/; every grid row and every ensemble member is published.
- Full complex X/Y/Sigma equations, including both mediator interactions
  and quartics; no angular-only EFT. The inflationary Hubble mass is removed
  instantaneously into a prescribed radiation background. Phi=0,
  zero initial physical velocities, instant visible reheating.
- Couplings are the unchanged classical point, referenced explicitly to
  M_Y=1e12 GeV. This is not RG improvement or a quantum matching calculation;
  the review's high-field perturbativity concern remains.

The supplied mismatch rotates Sigma with the full *aligned* inflationary
radial amplitudes held fixed, leaving the kappa phase aligned. It is a
favorable preparation assumption, not an interaction found in the baseline
action. A preparation control instead minimizes radii with that phase held
externally; the leading point then gives only |epsilon_X|=0.0516. This
difference shows that the physical source/preparation matters.

## Three nominal high-fraction candidates

Signed fractions are retained; conjugating the preparation reverses the
sign. Sigma here is the phase-seed standard deviation, in radians. Errors
are Monte Carlo sampling SEMs, not physical or systematic uncertainties.

| c_H | H_i [GeV] | delta | Nominal epsilon_X | Net epsilon, sigma=0.001 | Net epsilon, sigma=0.03 | Conditional Y pairs / net charge |
|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

The leading point's 0.03-radian ensemble has 128 members. Its charge
fraction is -0.24909 ± 0.00333. The nominal value changes by less than
7e-9 between relative tolerances 2e-12 and 2e-14. The late charge fraction
remains about 0.241 when integration extends from a=1000 to a=3000.
The two other nominal candidates lose most of their net charge under the
larger seed experiment. The phase-noise distributions are supplied
sensitivity tests, not a computed inflationary stochastic state.

The local map resolves H within ±2% and five phases within ±0.06 radians.
It is structured and narrow: it must not be described as an entire smooth
H interval excluded or allowed by a handful of samples.

## What happens to the charge-fraction test when Y is included

For the leading point, the measured late quadratic excitation number is
J_Y≈13.269 in H_i units, with N_X≈-238.649. The number estimate subtracts
the late forced Y solution and is stable between late windows. The quoted
3-body width is Gamma_Y=3 kappa^2 M_Y/(256 pi^3)≈377.95 GeV.
Extrapolating a conserved Y number gives a_decay≈1.72e5 and Y excitation
energy about 275 times the X-condensate energy then.

Mapping that extra energy to the same quartic circular-charge scale gives
an **energy-equivalent** fraction about 0.00356. This is a warning that the
X-only 0.19 benchmark cannot simply be copied after the injection. It is
not a thermalization calculation, a new freeze-out law, or a legitimate
replacement for the review's coupled leftover calculation.

With unprocessed Y daughters, the source alone is about 17% of calibrated
dark-matter mass after normalizing the net yield to Y_Q. No dilution is
implemented here; the raw yield itself is about 66 Y_Q. In an equilibrating
history the source instead changes entropy and requires a new relic
calculation. Erasure through lumps is not assumed. Extra vector production,
Y/Sigma decay and the scalar transition are not silently omitted from a
claimed final population: that population has not been computed.

## Reproduction and numerical evidence

The independent cosmic-time integration agrees with the C conformal
integration at relative state error below 2.2e-10 in the short controls.
The old three-charge release reproduces to 1.6e-11. Force finite differences
agree below 6e-10; the gauge Ward identity below 1e-16. Exact real alignment
and CP reflection are preserved. An infinitesimal homogeneous phase seed
grows strongly, confirming the instability correction without inferring a
spatial domain spectrum from k=0 alone.

Across coarse/refined/local grids the maximum normalized Gauss residual
is {f['max_grid_Gauss_error']:.3g}, and the maximum canonical explicit-work
residual is {f['max_grid_work_error']:.3g}. Candidate controls are tighter.
These are numerical controls, separate from physical viability gates.

[REPRODUCE.md](REPRODUCE.md), [MATHEMATICS.md](MATHEMATICS.md),
[CORRECTIONS.md](CORRECTIONS.md), [SOURCES.md](SOURCES.md) and
[FINDINGS.json](FINDINGS.json) give methods, changes, assumptions and data.

## Decision and next review

The original high-scale aligned-release diagnosis is reproduced:
mean epsilon={f['baseline_high_scale_mean']:.5f} ± {f['baseline_high_scale_SEM']:.5f},
local rms {f['baseline_high_scale_local_rms']:.3f}. Lowering H can preserve
a supplied common asymmetry, so the homogeneous scan **does not rule out
the one-mediator link at these couplings by a universal charge-fraction
failure**. The existing action still does not implement that preparation.

No complete cosmological point passes all six requirements. The next
useful reviewer task is to reproduce the leading explicit point and add
its measured Y excitation to a decay/thermalization ledger, then test a
specified phase-preparation interaction if that ledger leaves room.
Longitudinal vectors require the evolving Higgs background; their initial
Goldstone-limit normalization is supplied in VECTOR_INITIALIZATION.json,
but no final occupation is claimed. Wall capture, string formation and
the full object spectrum remain open. Vacuum density and the C* observer
dictionary remain underived.

All computations here are same-author AI work. Claude's new review of this
packet is pending; a posted handoff is not evidence it has run.
'''
    (ROOT/'REPORT.md').write_text(report)


if __name__=='__main__':main()
