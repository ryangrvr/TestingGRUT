#!/usr/bin/env python3
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from boundary import ROOT
from manifest import main as manifest_main


def main():
    b=json.loads((ROOT/'BOUNDARY_RESULTS.json').read_text())
    a=json.loads((ROOT/'TRACE_AUDIT_RESULTS.json').read_text())
    v=json.loads((ROOT/'VERIFICATION_RESULTS.json').read_text())
    cases={r['name']:r for r in b['cases']};runs={r['name']:r for r in a['rows']}
    def reflection(name,k):
        return min(cases[name]['samples'],key=lambda r:abs(r['k_over_m']-k))['reflected_energy']
    def percent(value):
        return '<0.001%' if 100*value<.001 else f'{100*value:.4g}%'
    table='\n'.join(f"| {k:g} | {percent(reflection('original',k))} | {percent(reflection('old_larger_box',k))} | {percent(reflection('candidate_very_wide',k))} |" for k in [.01,.02,.03,.05,.1,.2,.5,1.])
    audit_table=[]
    for name,row in runs.items():
        x=row['variants'][1]
        audit_table.append(f"| {name} | {100*row['core60_late_mean_fraction']:.5f}% | {100*x['conditional_controlled_energy_fraction']:.3f}% | {100*x['conditional_all_band_antiparticle_number_fraction']:.4f}% |")
    base=runs['N0125'];weaker=runs['N0125_abs']
    core_shift=abs(weaker['core60_late_mean_fraction']/base['core60_late_mean_fraction']-1)
    pair_shift=weaker['variants'][1]['conditional_all_band_antiparticle_number_fraction']/base['variants'][1]['conditional_all_band_antiparticle_number_fraction']-1
    fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
    names=[('original','Original: D=90, eta=0.03','-'),
           ('old_larger_box','Old larger box: D=120, eta=0.03','-'),
           ('candidate_wide','Design only: D=360, eta=0.1','--'),
           ('candidate_very_wide','Design only: D=720, eta=0.1','--')]
    for name,label,style in names:
        data=cases[name]['samples'];ax.loglog([x['k_over_m'] for x in data],[100*x['reflected_energy'] for x in data],style,label=label)
    ax.axhline(1,color='black',ls=':',label='1% target')
    ax.set(xlabel=r'Test-wave momentum $p/m_X$',ylabel='Reflected energy (%)',
        title='The ordinary numerical layer fails a uniform soft-wave target',ylim=(1e-5,100))
    ax.grid(True,which='both',alpha=.15);ax.legend(fontsize=9)
    fig.savefig(ROOT/'boundary_reflection.png',dpi=180);plt.close(fig)
    text=f'''# The soft-field receiver has a measured boundary problem

23 October 2026 · Codex / AInstein · parent 59a1a71

**The old absorber is not uniformly transparent to the soft radiation we
need to measure.** An independent linear boundary calculation returns about
31.7% of test-wave energy at p/m_X=.02, and 14.8% at p/m_X=1, for its original
parameters. This is an additional reason the missing soft field cannot yet
be counted as asymptotic particles. It does not overturn the separately
measured finite-time core charge or supply a final cosmological population.

## What was built and checked

The actual nodal solver's cubic velocity damper and outer wall were reduced
to their continuum free radial Klein-Gordon scattering problem. The new
code computes reflection including both the layer entrance and wall return.
It checks that result against a separate finite-difference boundary-value
solver on six spacings, two closed-form problems, a positive flux-loss
integral and an ODE-tolerance control. No new full-field collapse simulation
was executed in this packet.

All six saved nodal traces were freshly loaded and matched to their original
float64 hashes. The measured signed energy/charge currents were reanalyzed
with actual nodal radii, window variations and the existing local free-field
cut. The two-probe pair has bitwise identical final fields and physical
histories, providing a radius control on the same physical realization.

The proposed receiver target is less than 1% reflected energy throughout
.01 <= p/m_X <= 1. Exploratory results preceded this engineering target;
it is not a blinded preregistration. All tested ordinary cubic layers fail
the target. Favorable narrow momentum ranges are not a complete receiver.

## The boundary benchmark

These are linear unit-test-wave coefficients, **not the fraction of the
formation debris that actually reflected**. The last column is a wider
ordinary-damper design, not a newly executed collapse. The free continuum
benchmark does not by itself certify the production discretization or a
nonlinear exterior.
Tiny computed coefficients are reported only as below 0.001%; their many
raw floating-point digits are not certified reflection floors.

| p/m_X | Original D=90, eta=.03 | Old larger box D=120, eta=.03 | Design D=720, eta=.1 |
|---:|---:|---:|---:|
{table}

The original layer returns 57.1% at p/m_X=.01. Even the tested D=720 design
returns 19.9% there. Stronger damping alone is not a reliable soft-wave cure:
it can increase the impedance mismatch at the entrance. The measured
positive damping-energy receiver still closes the numerical ledger, but
that ledger does not distinguish a true outgoing continuum from radiation
that returned and interacted again.

## Core convergence does not certify pair content

The following energy and antiparticle columns remain **conditional linear
mode diagnostics**. They are not physical source completeness or a measured
antiparticle abundance. The controlled-energy denominator is the inferred
outgoing linear-mode energy, which also has finite-time assumptions.

| Existing run | Late charge inside radius 60 | Conditional energy passing local cut | Conditional all-band antiparticle number fraction |
|---|---:|---:|---:|
{chr(10).join(audit_table)}

Halving absorber strength shifts the core fraction by only
{100*core_shift:.3f}% relatively, while the conditional all-band antiparticle
estimate shifts by {100*pair_shift:.2f}% relatively. That estimate is not
stable enough to supply the real anti-X population. At the original probe,
only about 25.1% of inferred outgoing linear-mode energy passes the existing
local cut; the remaining portion is not normalized into particles.

The 12000-unit record also has Fourier spacing about 5.24e-4. At p/m_X=.01,
the frequency lies only about 5e-5 above the vacuum mass. This FFT method
needs a characteristic duration around 1.26e5 to resolve that excess,
alongside appropriate arrival time, window, boundary and nonlinearity
controls. Smaller spatial cells alone cannot resolve this limitation.

## A concrete next receiver

A homogeneous-exterior complex-stretch reference is derived explicitly in
[MATHEMATICS.md](MATHEMATICS.md). In the continuum it has no entrance
reflection; a cubic stretch with D=90 and sigma_max=40 would suppress the
wall's returned energy below 2.4e-16 at p/m_X=.01. This is a frequency-domain
design reference from established numerical analysis. **No new time-domain
PML, stability certificate or complete formation source is implemented here.**

The deciding build must provide a causal massive-wave exterior, specify its
initial memory, pass discrete reflection/stability controls and preserve a
charge/energy plus time/momentum source interface. Its nonlinear matching
region and bound component need separate controls. A scalar multiplication
of the present spectrum by a reflection correction cannot supply that result.

## Verification and the remaining theory work

{v['passed']}/{v['total']} numerical control groups pass. The independent BVP's
finest complex reflection-amplitude error is about
{v['numerical_controls'][1]['evidence']['max_finest_complex_amplitude_error']:.2g};
its observed order is consistent with its declared second-order stencil.
Raw validation was actually performed for {v['raw_validation']['trace_count']}
traces and the matched two-probe fields. These are author checks, not an
external replication. **The original receiver target FAILS and the complete
asymptotic particle source remains uncertified.** That failure is separate
from the passing numerical controls.

The previous [conversion-history calculation](../grut_conversion_streaming_20261023/REPORT.md)
remains explicitly conditional on its supplied source and frozen objects.
This audit does not replace it with measured formation daughters. General
hot-core leakage, actual cosmological formation, evolving profiles,
mixed-matter perturbations and vacuum density remain unresolved.

![Linear reflection benchmark for the original layer and tested wider designs](boundary_reflection.png)

[Boundary results](BOUNDARY_RESULTS.json) · [trace audit](TRACE_AUDIT_RESULTS.json) ·
[controls](VERIFICATION_RESULTS.json) · [equations](MATHEMATICS.md) ·
[reproduce](REPRODUCE.md) · [Claude handoff](CLAUDE_HANDOFF.md) · [sources](SOURCES.md).
'''
    (ROOT/'REPORT.md').write_text(text)
    manifest_main()


if __name__=='__main__':main()
