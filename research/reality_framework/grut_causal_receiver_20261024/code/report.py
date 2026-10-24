#!/usr/bin/env python3
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from run import ROOT
from manifest import main as manifest_main


def read(name):return json.loads((ROOT/name).read_text())


def main():
    layer=read('DISCRETE_LAYER_RESULTS.json');controls=read('VERIFICATION_RESULTS.json')
    audit=read('SOURCE_AUDIT_RESULTS.json');interface=read('I025_RESULTS.json');bound=read('BOUND_MODE_DIAGNOSTICS.json')
    grid=[read(f'{name}_RESULTS.json') for name in ['P025_time','P0125','P00625']]
    refined=grid+[read('P025_remote_time_RESULTS.json')]
    near=grid[0];far=refined[-1];long=read('P025_remote_RESULTS.json')
    packets=[read(f'{name}_RESULTS.json') for name in ['free_X','free_X_time','free_s','free_soft']]
    max_rx=max(r['reflected_energy'] for r in layer['rows'] if r['mass']==1.)
    max_rs=max(r['reflected_energy'] for r in layer['rows'] if r['mass']!=1.)
    values=[r['core60_late_mean_fraction'] for r in grid];spread=max(values)/min(values)-1
    old_core=read('../grut_core_conversion_20261023/N00625_RESULTS.json')['core60_late_mean_fraction']
    grid_table='\n'.join(f"| {r['name']} | {r['parameters']['h']:g} | {r['parameters']['dt']:g} | {100*r['core60_late_mean_fraction']:.6f}% | {r['maximum_energy_plus_receiver_relative_error']:.3g} |" for r in refined)
    packet_table='\n'.join(f"| {r['name']} | {r['parameters']['tmax']:g} | {r['free_comparison']['exact_time_physical_energy_error_over_initial']:.3g} | {r['maximum_energy_plus_receiver_relative_error']:.3g} |" for r in packets)
    mode_table='\n'.join(f"| {r['name']} | {len(r['bound_modes'])} | {r['Qbar_in_instantaneous_subthreshold_modes']:.3f} | {r['Qbar_in_other_subthreshold_modes']:.3f} |" for r in bound['rows'])
    fig,axes=plt.subplots(1,3,figsize=(16,4.7),layout='constrained')
    old=read('../grut_soft_receiver_20261023/BOUNDARY_RESULTS.json')['cases'][0]['samples']
    axes[0].loglog([r['k_over_m'] for r in old],[100*r['reflected_energy'] for r in old],label='Old charged-wave layer')
    for mass,label in [(1.,'Causal charged layer'),(np.sqrt(6.),'Causal neutral layer')]:
        rows=[r for r in layer['rows'] if r['h']==.125 and r['mass']==mass]
        axes[0].loglog([r['k'] for r in rows],[100*r['reflected_energy'] for r in rows],'-o',ms=3,label=label)
    axes[0].axhline(1,color='black',ls=':',label='1% target')
    axes[0].set(xlabel=r'Test momentum $p/m_X$',ylabel='Reflected energy (%)',ylim=(1e-9,100),title='Actual grid and outer wall')
    axes[0].legend(fontsize=8);axes[0].grid(True,which='both',alpha=.15)
    axes[1].plot([r['parameters']['h'] for r in grid],[100*v for v in values],'-o')
    axes[1].axhline(100*old_core,color='gray',ls=':',label='Previous operator/layer')
    axes[1].set(xlabel='Spatial step h',ylabel='Charge inside radius 60 (%)',title='Same refined time step, finite-time means')
    axes[1].grid(True,alpha=.2);axes[1].legend(fontsize=8)
    history=np.array(interface['history']);axes[2].plot(history[:,0],history[:,10]/interface['Ebar_initial'],label='Net energy through cut')
    axes[2].plot(history[:,0],history[:,11]/interface['Qbar_initial'],label='Net charge through cut')
    axes[2].set(xlabel=r'Time $m_X t$',ylabel='Fraction of initial total',title='Exact classical extraction interface')
    axes[2].legend(fontsize=8);axes[2].grid(True,alpha=.2)
    fig.savefig(ROOT/'causal_receiver_results.png',dpi=170);plt.close(fig)
    body=f'''# The causal exterior receiver is implemented and tested

24 October 2026 · Codex / AInstein · parent 9636950

**The ordinary-damper reflection problem now has an executable causal
replacement.** Across {len(layer['rows'])} sampled combinations of grid,
momentum and mass, the maximum recorded reflected energy is about
{100*max_rx:.3g}% for the charged field and {100*max_rs:.4g}% for the neutral field. Four propagating
packet tests agree with independent free-field references. New prepared
formation runs and an exact conserving extraction interface were executed.
This is a numerical construction result, not completion of the physical
cosmological particle source or the working theory of reality.

## What changed

The new layer stores causal memory explicitly. It retains the physical
quartic force and uses a fourth-order staggered gradient/divergence in the
physical region. The old and new discrete spatial Hamiltonians are different;
the same action and smooth Gaussian preparation are directly refined again.
No physical damping interaction, portal width, initial charge-yield retuning
or new GRUT-RAI law is introduced.

The implemented direct second-order form, its zero-memory convention,
constant-coefficient stability argument, actual lattice leads, finite-wall
test and signed boundary work are derived in [MATHEMATICS.md](MATHEMATICS.md).
It follows established Klein-Gordon PML analysis, credited with exact read
scope in [SOURCES.md](SOURCES.md). The homogeneous matching result does not
become a nonlinear theorem by retaining the physical nonlinear force.

## Discrete and time-domain tests

The frequency calculation includes actual node/face memory factors, all
lattice spatial roots and the outer wall. A zero-layer closed-wall amplitude
checks its normalization independently. Complete representative C++ linear
eigenspectra have no growing modes at the recorded RK4 steps. Executable
Hamiltonian and Noether finite variations also check the receiver rates.

The following errors are physical-region error **energy divided by initial
packet energy**, not a pure reflection coefficient. Same-RK4 bulk comparisons
separately isolate boundary error and are recorded in each result. Reference
fast-mode tails are explicitly estimated.

| Packet test | Duration | Error versus exact-time free reference | Max field-plus-receiver energy drift |
|---|---:|---:|---:|
{packet_table}

The long soft-packet case has a deliberately coarser bulk time step. Its
ledger drift is disclosed and is not called physical emission or silently
required to meet the refined formation target. The complete physical-region
error remains below one percent by a wide margin.

## Prepared formation and boundary controls

These are new classical spherical runs of the same dense prepared Gaussian,
not a primordial-density 3D formation simulation. Late means use tau=9600
to 12000 at radius 60. Q=200000 remains the separate supplied stationary
scattering example and is not the object charge predicted by these runs.

| Run | h | dt | Late retained charge fraction | Max energy-plus-boundary error |
|---|---:|---:|---:|---:|
{grid_table}

The three spatial refinements have relative spread {100*spread:.4g}%.
Their mean values and observed differences are reported without a certified
Richardson extrapolation. At fixed h and dt, moving the layer from radius
360 to 720 shifts the finite-time mean core charge by
{100*abs(audit['boundary_radius_control']['core_charge_relative_shift']):.4g}% relatively.
The raw squared field difference inside radius 300 is
{audit['boundary_radius_control']['raw_complex_field_L2_squared_relative_difference_inside_300']:.3g}
of its norm; the separately disclosed phase-aligned value is not substituted
for that raw comparison.

The earlier five-point/ordinary-layer finest value was {100*old_core:.5f}%;
the new finest value is {100*values[-1]:.5f}%, a relative change of
{100*(values[-1]/old_core-1):.4g}%. Both values are preserved with their
numerical setups. The old small internal grid spread was not a complete
systematic uncertainty for a boundary that later failed its soft-wave test.

The two coarse time-step controls fail the stricter 1e-5 refined energy
target: their errors are {controls['coarse_budget_controls'][0]['energy_error']:.3g}
and {controls['coarse_budget_controls'][1]['energy_error']:.3g}.
They remain visible as controls; the maximum refined energy error is
{max(r['maximum_energy_plus_receiver_relative_error'] for r in refined):.3g}
and charge error is {max(r['maximum_charge_plus_receiver_relative_error'] for r in refined):.3g}.

The separate longer coarse run reaches tau=24000 and has late mean core
fraction {100*long['core60_late_mean_fraction']:.5f}%. Its different late-time
window is not substituted into the spatial comparison. This continuing
classical evolution is another reason not to label the radius-60 charge an
asymptotic cold-object or quantum lifetime measurement.

## An exact classical source interface

The instrumented I025 run adds the discrete hopping currents across radius
{interface['exact_cut_radius']:.3f}. Energy inside plus integrated net energy
through the cut closes to {interface['maximum_cut_budget_relative_error']['energy']:.3g};
charge closes to {interface['maximum_cut_budget_relative_error']['charge']:.3g}.
Physical fields, probe traces and layer memory are bitwise identical to the
otherwise identical uninstrumented P025_time run.

By tau=12000 this cut has passed net dimensionless energy
{interface['final_net_Eflux']:.6f} and charge {interface['final_net_Qflux']:.6f}.
No low-frequency cut is needed for that classical conservation statement.
The compact history records extraction time. It is not yet the production
time and momentum of a nonnegative particle ensemble, or a selected cosmic
formation epoch. Signed currents retain returning waves and coherent fields.

## Bound/coherent separation is still a physical source gate

The vacuum numerical boundary passes its declared tests. The nonlinear
field-to-particle interpretation remains separate. Layer mass shifts reach
about {near['maximum_layer_charged_mass2_shift']:.3g} at R=360 and
{far['maximum_layer_charged_mass2_shift']:.3g} at R=720, comparable to or
larger than the kinetic k^2 of the very slow modes. The existing local cut
therefore still leaves much of the inferred outgoing linear energy
uncontrolled. No missing energy is normalized into particles.

An instantaneous bound-mode projection is now implemented using the actual
snapshot charged operator, with the production stencil applied independently
to check its eigenvectors. This is a frozen, single-field diagnostic. It
does not replace the time-dependent coupled fluctuation problem or prove
future stable binding/free escape.

| Snapshot | Instantaneous subthreshold modes | Canonical charge in those modes | Charge in modes not 99% localized inside radius 60 |
|---|---:|---:|---:|
{mode_table}

The time/resolution dependence prevents banking that instantaneous split
as a completed bound-core subtraction. The exact classical source is now
measured; its complete asymptotic particle time/momentum spectrum remains
uncertified. The old conditional conversion-history horizons keep their
stated imposed-source and fixed-profile scope.

## Completed validation and review surface

{controls['passed']}/{controls['total']} declared numerical/finite-time control
groups pass. Raw validation checked {controls['raw_validation']['arrays']}
lossless numeric arrays. This is one author's construction and verification;
no second full replay, independent Claude rerun or empirical validation is
claimed. Six prepared runs, four free packets, the instrumented repeat and a
short smoke run are recorded; exploratory first-order pilot runs are excluded.

The compact packet includes code, run summaries, conserving histories,
diagnostics and hashes. Bulk lossless fields/traces remain reproducible and
are not uploaded here. The next physical calculation must separate the
remaining bound/coherent nonlinear field from a free outgoing source, then
use its actual production history in evolving-object reaction and
mixed-matter perturbation calculations. Actual cosmological/3D formation,
general hot-core leakage and vacuum density remain open.

![Causal reflection tests, spatial refinement and exact extraction history](causal_receiver_results.png)

[Review checks](CLAUDE_HANDOFF.md) · [reproduce](REPRODUCE.md) ·
[numerical validation](VERIFICATION_RESULTS.json) · [source audit](SOURCE_AUDIT_RESULTS.json) ·
[bound-mode diagnostics](BOUND_MODE_DIAGNOSTICS.json) · [exact interface](I025_RESULTS.json).
'''
    (ROOT/'REPORT.md').write_text(body);manifest_main()


if __name__=='__main__':main()
