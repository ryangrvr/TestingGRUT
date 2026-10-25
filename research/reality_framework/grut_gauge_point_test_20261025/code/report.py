#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]


def main():
    pair=json.loads((ROOT/'PAIR_RESULTS.json').read_text());point=json.loads((ROOT/'POINT_RESULTS.json').read_text())
    screen=json.loads((ROOT/'SCREENING_RESULTS.json').read_text());release=json.loads((ROOT/'RELEASE_RESULTS.json').read_text())
    verify=json.loads((ROOT/'VERIFICATION_RESULTS.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.1))
    for d,color in [(5,'#225ea8'),(6,'#cb181d')]:
        r=next(r for r in pair['rows'] if r['power']==d and r['phase_points']==384)
        w=np.array(r['sampled_integrand'])
        axes[0].plot(r['sampled_q'],w/max(w),lw=2,color=color,label=fr'$X^{d}$')
    axes[0].set(xlim=(0,max(r['invariant_q_max'] for r in pair['rows'])),ylim=(0,1.06),xlabel=r'Total pair momentum $|\mathbf{p}_1+\mathbf{p}_2|/m_X$',ylabel='Vacuum Born pair weight / peak',title='Pairs can leave back to back')
    axes[0].legend(frameon=False)
    data=release['rows'][1]['samples']
    for key,label,color in [('comoving_X_charge','X','#225ea8'),('comoving_Y_charge','Y','#e6550d'),('comoving_Sigma_charge',r'$\Sigma$','#238b45')]:
        axes[1].plot(data['t'],data[key],label=label,color=color,lw=1.5)
    axes[1].axhline(0,color='#555555',ls='--',lw=1,label='Aligned source: all zero')
    axes[1].set(xlabel=r'Time $H_{rel}t$',ylabel=r'Comoving species charge / $H_{rel}^3$',title='Supplied misalignment transfers internal charge')
    axes[1].legend(frameon=False,fontsize=8,loc='best')
    for ax in axes:ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    fig.text(.5,.012,'Vacuum pair estimate and short homogeneous release diagnostic; no final abundance or complete lifetime is inferred.',ha='center',fontsize=8.5,color='#444444')
    fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(ROOT/'pair_and_charge_transfer.png',dpi=170);plt.close(fig)
    table='\n'.join(f"| {r['test']} | {r['baseline']} | {r['qualification']} |" for r in verify['six_physical_tests'])
    rates='\n'.join(f"| X^{r['power']} | {r['coherent_rate_per_second']:.4g} | {r['free_Born_pair_rate_per_second']:.4g} | {r['coefficient_max_for_one_percent']:.4g} |" for r in pair['summaries'])
    inf=point['inflation'];wall=point['wall_bias_variant'];g=pair['matched_gauge_point'];limits=point['symmetric_remainder']['extra_pair_budgets'][0]
    masses=np.sqrt(np.maximum(inf['angular_mass2_over_H2'][1:],0))
    s5=pair['summaries'][0]['S_min_if_abs_c_equals_exp_minusS']
    text=f'''# The gauge point survives late loss and binding; its genesis is not established

25 October 2026 · Codex / AInstein · parent 7737487 · response to 2e9f296

**The six-condition test does not yet yield a complete cosmological point.**
The missing pair-emission correction is accepted and independently
integrated. Cold classical binding has a large screened-gauge margin.
The early single-field Affleck-Dine window cannot be copied to this gauge
point, and its own phase alignment generates zero charge. Symmetric
leftovers and the full cosmological source remain open.

Eight author mathematical/numerical control groups pass. Those passes
are separate from the physical dispositions below.

## Stated point and six dispositions

Retain M_Y=1e12 GeV, h=1e11 GeV, mu=5e11 GeV, mediator y=0.001,
e=1e-7 and lambda_Y=lambda_Sigma=1. The original light target has
g_X=0.1, lambda_X=0.001, lambda_phi=0.03, v=0.924089637 GeV and Q=200000.
Here y is called kappa in the parent code and is not the Higgs portal.
No old yield, scale or empirical target is refit.

| Required test | Baseline disposition | Scope/result |
|---|---|---|
{table}

## The corrected survival comparison

The old forced-response numbers described coherent single-quantum
emission. We omitted spontaneous pair emission from the quadratic term
in X^d. With unit coefficient and reduced Planck suppression:

| Hypothetical operator | Coherent fractional rate (/s) | Vacuum Born pair fractional rate (/s) | Coefficient ceiling for the chosen 1%/13.8-Gyr target |
|---|---:|---:|---:|
{rates}

Our one-dimensional invariant phase-space integral differs from the
review's individual-momentum/angular integral by only a few ppm. This
checks the *same vacuum approximation* with different quadrature; it is
not an in-medium Bogoliubov calculation. No X5 term is added to the Z6 model.
The corrected c5 ceiling is about 1.80e-5. If c5=exp(-S) with unit
prefactor, its translation is S>={s5:.3f}; no instanton coefficient is predicted.

At the matched gauge point, c6 at reduced Planck normalization is
{g['effective_c6']:.6g}. Coherent plus vacuum-pair loss is
{g['included_channels_rate_per_second']:.6g}/s, giving an initial rate-times-age
of {g['included_channels_rate_times_age']:.6g}. This passes the engineering target
with a large margin. Other medium/multiparticle channels, evolution and
population averaging are uncomputed, so this is an included-channel estimate.

## Early matching and the phase obstruction

For the review's supplied H_inf=4.7e13 GeV and c_H=1, the vacuum
mediator mass/H is 0.0213, not a heavy-field limit. Its quadratic
elimination parameter |y||X|^2/M_Y^2 is about 1100 at the proposed
single-field starting amplitude. The full Y and Sigma quartics stabilize
the high-field potential and must remain in the calculation.

With the additional assumed -c_H H^2|X|^2 preparation, the full tree radial
minimum is approximately

| Field amplitude | GeV |
|---|---:|
| abs(X) | {inf['full_radial_amplitudes_GeV'][0]:.5g} |
| abs(Y) | {inf['full_radial_amplitudes_GeV'][1]:.5g} |
| abs(Sigma) | {inf['full_radial_amplitudes_GeV'][2]:.5g} |

The Higgs amplitude is about {inf['full_Higgs_amplitude_over_vacuum']:.1f} times its vacuum value.
The two physical angular masses are {masses[0]:.3f}H and {masses[1]:.3f}H;
one mode is light. Constant-coupling tree extrapolation and quantum running
are distinct questions; the latter is not computed. These masses do not
supply a CMB isocurvature prediction for the later X population.

More directly, real constant couplings preserve the real-field subspace.
Starting at the potential's own phase minimum with zero angular momenta
therefore gives exactly zero charge, throughout arbitrary real radial
evolution. Three full complex-field release histories and a tighter run
check the zero source, opposite conjugate sources and Gauss/stress ledgers.
The inflaton-dependent phase mismatch assumed in the review is not in the
stated action. Its construction remains necessary.

A **proposed additional** nonminimal coupling xi_X=-1/12 can realize the
c_H=1 curvature mass in de Sitter and suppress it in radiation domination.
It does not select the phase, and neither reheating nor the full
scalar-tensor background is solved. The short diagnostic ends at a~14,
well before a final asymmetry, heavy-particle decay or fragmentation.

## Cold gauge binding and compact discrete charge

The original object has a free-particle binding margin of
{screen['free_particle_binding_margin_GeV']:.6g} GeV. An admissible screened
classical trial adds at most {screen['added_gauge_energy_bound_GeV']:.6g} GeV,
only {screen['bound_over_existing_binding_margin']:.3g} of that margin. The
screening Fourier/Parseval control agrees to 6e-14. This supports binding;
it does not certify the entire coupled stability spectrum.

Q=200000 is **2 modulo 6**. Continuous Higgs screening is not the complete
discrete-gauge quantum-state construction. A compact whole can, for example,
group three such cores with total Q_X=600000, compensated by integer
Sigma number -100000. The quantum population state and its distribution
have not been constructed. The new vacuum heavy masses exceed this cold
object's total energy, closing their on-shell production at this point.

## Explicit wall-bias variant and leftover matter

The baseline lacks a sign bias. A separately proposed renormalizable tilt is

V_bias=Delta/4[2-3phi/v+(phi/v)^3], Delta=1e-22 GeV^4.

It preserves stationary +/-v, sets their energies to 0 and Delta, and
keeps the full potential positive. With stated scaling-network parameters,
the minimum estimated annihilation temperature is
{1000*wall['minimum_estimated_annihilation_temperature_GeV']:.3g} MeV, above the
adopted 1-MeV reference. This is a timing estimate, not a network simulation;
recent simulations exhibit different bias scaling. Its extra core trial
energy is {screen['proposed_bias_variant_added_trial_energy_GeV']:.3g} GeV.
The zero at +v is an explicit reference choice, not a Lambda solution.
Exact old source/virial identities change under this new operator.

At a chosen 1% allowance above the unchanged calibrated abundance, extra
cold free pairs are limited to Y_pair<={limits['extra_cold_pair_yield_max']:.4g},
or {limits['extra_pairs_per_net_charge_max']:.5g} pairs per net charge. Under
the review's approximate remainder model this needs epsilon>={limits['epsilon_min_in_quoted_remainder_model']:.5g}.
No annihilation/cannibalization or anti-object distribution is solved here.
A sudden-quench scalar-debris benchmark spans 0.34-4.2 times the calibrated
matter, but it is not a lower bound: adiabatic evolution may avoid it and
an existing total-energy ledger may already include it. Double counting
that energy would be incorrect.

Vacuum kinetic mixing is zero under the stated matching choice, no
bi-charged matter and real-coupling dark charge-conjugation symmetry.
This EFT assumption is not a quantum-gravity prediction; a new phase/CP
sector requires another mixing audit.

## Corrections and validation

Klinkhamer-Volovik already proposed the disequilibrium interpretation.
It is prior art; GRUT's closure link remains its wager. Competing observer
prescriptions give different assignments and GRUT has implemented none.
The precise corrections and primary read scopes are retained below.

Controls include an analytic Gaussian, six pair quadratures, nine
high-field starts, an independent phase Hessian, four full release histories,
three binding quadratures, a Fourier screening check and exact bias/mod-six
algebra. Gauss error is below 2.1e-14 and scalar stress error below 2.7e-13
in the short release. These controls do not turn the physical open gates
into passes. No complete external-review replay is claimed.

![Pair estimate and compensated charge diagnostic](pair_and_charge_transfer.png)

[Equations and limits](MATHEMATICS.md) · [corrections](CORRECTIONS.md) ·
[primary sources/read scopes](SOURCES.md) · [pair results](PAIR_RESULTS.json) ·
[point results](POINT_RESULTS.json) · [screening](SCREENING_RESULTS.json) ·
[full-field release](RELEASE_RESULTS.json) · [checks and six gates](VERIFICATION_RESULTS.json) ·
[reproduce](REPRODUCE.md) · [review handoff](CLAUDE_HANDOFF.md).
'''
    (ROOT/'REPORT.md').write_text(text)
    print('REPORT_AND_FIGURE_WRITTEN')


if __name__=='__main__':main()
