"""Render the scientific report from the assessed measurements."""
from pathlib import Path
import json
P=Path(__file__).resolve().parent
def pct(v):return f'{100*v:.5f}%'
def yes(v):return 'pass' if v else 'FAIL / unresolved'
def main():
 s=json.loads((P/'SUMMARY.json').read_text());c=s['cases'];tr=s['direct_transverse'];e=json.loads((P/'EMPIRICAL_PREDICTION.json').read_text());geo=s['empirical_geometric']
 single='\n'.join(f"| {tag} | {c['single_generic_'+tag]['case']['h']} | {pct(c['single_generic_'+tag]['minimum_late_Q12_fraction'])} | {c['single_generic_'+tag]['late_core_EQ12_ranges'][0][0]:.6f}–{c['single_generic_'+tag]['late_core_EQ12_ranges'][0][1]:.6f} |" for tag in ['coarse','medium','fine'])
 pairs='\n'.join(f"| {name} | {c[name]['late_core_Q12_means'][0]:.6f} | {c[name]['late_core_Q12_means'][1]:.6f} | {pct(c[name]['minimum_late_Q12_fraction'])} | {yes(c[name]['retention_pass'])} |" for name in ['pair_reference_fine','pair_generic_fine','pair_negative_fine'])
 grids='\n'.join(f"| {name} | {pct(max(v['coarse_medium']['maximum_late_Q12_difference_over_initial_Q']))} | {pct(max(v['medium_fine']['maximum_late_Q12_difference_over_initial_Q']))} | {yes(v['charge_convergence_pass'])} |" for name,v in s['grid_comparisons'].items())
 control='\n'.join(f"| {name} | {pct(max(v['maximum_late_Q12_difference_over_initial_Q']))} | {yes(v['charge_control_pass'])} | "+(pct(v['maximum_late_relative_shift_amplification_relative_difference'])+' / '+yes(v['relative_shift_control_pass']) if 'relative_shift_control_pass'in v else 'not used')+' |' for name,v in s['time_boundary_controls'].items())
 modes='\n'.join(f"| {r['h']:.6f} | {r['endpoint_amplification'][2]:.6f} | {pct(r['endpoint_symmetry_residual'][0])} | {pct(r['endpoint_symmetry_residual'][1])} | {pct(r['late_core_shape'][0]['fraction_explained_by_translation'])} / {pct(r['late_core_shape'][1]['fraction_explained_by_translation'])} |" for r in sorted(s['large_domain_original_m1'],key=lambda r:-r['h']))
 impact='\n'.join(f"| {name} | {v['endpoint_x_over_impact'][0]:.6f} | {v['endpoint_x_over_impact'][1]:.6f} | {max(v['maximum_core_anisotropy_change_vs_reference']):.6g} |" for name,v in tr['cases'].items())
 impactchecks='\n'.join(f"| {name} | {pct(max(v))} | {yes(tr['checks'].get({'coarse_medium':'space','medium_fine':'space','half_seed':'linearization','box':'box','dt':'time'}[name],False)) if name!='coarse_medium' else 'reference difference'} |" for name,v in tr['relative_response_differences'].items())
 production='\n'.join(f"| {tag} | "+' | '.join(str(c['production_'+tag]['components'][str(t)]['count'])+' / '+pct(c['production_'+tag]['components'][str(t)]['selected_fraction_of_total_Q']) for t in [.02,.05,.1])+' |' for tag in ['coarse','medium','fine'])
 spectra='\n'.join(f"| {r['R']} | {r['dr']} | {r['phase_l0'][0]:.9g} | {r['translation_l1'][0]:.9g} | {r['fixed_charge_l0'][0]:.9f} | {r['shape_l2'][0]:.9f} |" for r in s['spectral']['radial_spectra'])
 maxQ=max(v['max_relative_Q_drift'] for v in c.values());maxE=max(v['max_relative_E_drift'] for v in c.values());nret=sum(v.get('retention_pass',False) for v in c.values())
 text=f'''# Candidate 1: finite objects in 3D and a conditional empirical exclusion

20 October 2026. Experimental research branch; independent review pending.

The tested isolated object and outgoing collision cores survive the declared
finite three-dimensional perturbations through tau160. The isolated branch also
passes a numerical fixed-charge curvature test with an analytic ordering of all
angular harmonics. The old transverse canonical norm cannot be interpreted as
an internal-instability rate: its spatial symmetry controls retain large errors.
Direct full-3D impact experiments separately test the outgoing trajectory response.

One common-scale cosmological preparation has been carried to a conservative
gravitational observable. It fails: the conserved-charge energy floor requires
H0 >= {geo['H0_predicted_lower_bound']:.2f} km/s/Mpc, incompatible with independent
geometric megamaser data. This is a conditional exclusion, not an observation
of the proposed new fields or a viable empirical theory of everything.

## 1. What was actually computed

The common law is the unchanged [Candidate-1 action](../grut_candidate_1_20261018/GRUT_CANDIDATE_1.md).
The new main campaign contains **{s['case_count']} nonlinear full Cartesian 3D runs**:
26 finite-object/control runs and three condensate-production runs. Spatial
symmetry is not imposed. Generic seeds contain oblique shape deformation,
spatial phase variation and velocity kicks; pair seeds include transverse
offsets and perturbed velocity directions. Both seed signs are tested.

The [main protocol](protocol.json) froze 23 cases; the [direct-impact supplement](TRANSVERSE_PROTOCOL.json)
froze six more after a specific concern emerged from the symmetry audit. There
are also three large-domain replays of the original linear m=1 code and a
separate radial spectral calculation. Neither is counted as a nonlinear 3D run.
All freezes were local, not public preregistrations. Additive/posthoc boundaries,
development failures and an interrupted/restarted run are disclosed in
[DEVELOPMENT_NOTES.md](DEVELOPMENT_NOTES.md).

The Cartesian solver evolves the leading three-field classical sector, D=0,
with conservative fixed-vacuum boundaries for objects and periodic boundaries
for production. It uses second-order differences, the matching lattice
Hamiltonian, and kick-drift-kick. The main object grids are h=.75,.5,.375;
dt=.05, with independent halving controls. Maximum charge drift across the
29 completed runs is {maxQ:.3g}; maximum Hamiltonian drift is {maxE:.3g}.
The phase tangent is exact to printed roundoff. A centered finite-difference
tangent check reduces its error by a factor 4 when its increment is halved.

## 2. What survives, and how accurately

Retention is defined before the long runs: all sampled tau120..160 must retain
at least 90% of initial Q in R12 core balls, with at least 15% in each of a pair
and disjoint pair balls. A single object additionally requires local E/abs(Q)<1
and one dominant full-grid connected component. R8 charge covariances and
2/5/10%-of-peak connectivity thresholds are retained as separate diagnostics.
**{nret}/26 finite-object runs meet their applicable retention criterion.**

Perturbed isolated object:

| Grid | h | Minimum late R12 charge fraction | Late local E/Q range |
|---|---:|---:|---:|
{single}

Finest outgoing objects (late mean charges are dimensionless Qbar):

| Case | Smaller core mean Q | Larger core mean Q | Minimum joint retained Q | Criterion |
|---|---:|---:|---:|---|
{pairs}

Local laboratory E/Q is reported in SUMMARY.json. It is a sufficient binding
check when below one, not a necessary one for a moving object. A core ball is
not proof that the outgoing state has settled to an exact stationary solution.

Maximum late charge differences, normalized to initial total Q:

| Family | Coarse–medium | Medium–fine | Charge refinement criterion |
|---|---:|---:|---|
{grids}

| Control | Maximum late charge change | Charge criterion | Relative-shift norm change / criterion |
|---|---:|---|---|
{control}

Whole-event rotations by .3 radians rotate both positions and velocities. The
largest late core-charge difference falls from
{pct(max(s['whole_event_rotation_control']['medium']['maximum_late_Q12_difference_over_initial_Q']))}
at h=.5 to {pct(max(s['whole_event_rotation_control']['fine']['maximum_late_Q12_difference_over_initial_Q']))}
at h=.375. This is a distinct check from shifting the sources while holding
their velocity directions fixed.

These are finite-time, finite-amplitude-family results. At the empirical 100 GeV
point, tau160 is about {e['physical_reference_object']['tau160_seconds']:.4g} seconds,
or about 31.35 reference R90 light-crossing times. A cosmological lifetime of
each individual composite has not been calculated.

## 3. The earlier transverse response

The original finite-volume cylindrical m=1 implementation was replayed at
R=180,Z=192 through tau160. Full translation and full rotation generators act
on all fields **and all canonical momenta**. The finest relative-shift
amplification reproduces the earlier value about 9.064, but its symmetry controls
are not sufficiently accurate to use it as an internal-growth exponent.

| h | Relative-shift canonical amplification | Translation residual | Rotation residual | Geometric response explained by each core's translation |
|---|---:|---:|---:|---:|
{modes}

The geometric fit uses real amplitude and S within each R8 core, removing a
local rigid transverse displacement. The large-domain controls show that
returned waves alone do not explain the old mode solver's discrepancy. Its
discrete m=0/m=1 operators and nonlinear chain rule do not preserve continuous
translation/rotation exactly. The residuals decrease with refinement but remain
material at the tested late times. A symmetry tangent's positive canonical norm
also need not remain constant even in the exact continuum encounter.

The smaller Cartesian boxes have their own late, nearly resolution-independent
global symmetry residual. Boundaries break continuous spatial symmetry;
increasing the box and inspecting local observables are therefore necessary.
The residuals are exposed in [SYMMETRY_FIGURE.png](SYMMETRY_FIGURE.png), not counted
as passing controls.

The deciding additional experiment displaces initial source centers by opposite
per-source delta=.02 (relative impact parameter b=.04), with all velocity directions held fixed, then evolves the actual
nonlinear fields in 3D. A delta=.01 control tests linear response. These are
the original relative-shift preparation up to an overall sign convention.

| Direct case | Final smaller-core x/delta | Final larger-core x/delta | Maximum change in core anisotropy vs matched reference |
|---|---:|---:|---:|
{impact}

| Response comparison | Maximum relative centroid-response change | Frozen criterion |
|---|---:|---|
{impactchecks}

These direct core observables, their late trajectory fits, charge retention and
shape checks test scattering deflection without relying on the old global
canonical norm. Their exact pass/fail values are retained in SUMMARY.json.
The response over this interval is chiefly motion of enduring cores. An
asymptotic exponential internal-instability rate has not been established.

## 4. Isolated-object energy curvature in every angular sector

At fixed Q the radial amplitude/S Hessian includes a **positive rank-one charge
constraint**. Omitting it would incorrectly retain the negative radial direction
of unconstrained E−Omega Q. The phase block has its exact global-phase zero
mode; the ell=1 amplitude/S block has its exact translation zero modes. The
angular barrier orders all ell>=2 amplitude sectors and ell>=1 phase sectors.
Derivations, assumptions and normalization are in [MATHEMATICS.md](MATHEMATICS.md).

| Outer r | dr | Phase zero (sampled) | Translation zero (sampled) | Lowest fixed-Q radial curvature | Lowest ell=2 curvature |
|---|---:|---:|---:|---:|---:|
{spectra}

The small negative sampled symmetry eigenvalues vanish as dr squared. The
positive radial levels just above .36 are finite-box scattering levels near
the vacuum threshold, not newly identified internal resonances. Re-solving
the profile gives dQ/dOmega = {s['spectral']['dQ_dOmega_centered']:.3f}.
The result supports a local fixed-charge energy minimum of the isolated
nodeless branch, modulo phase and translation, with all angular harmonics
covered by operator ordering. It is not a rigorous interval enclosure, a
global-minimizer certificate, or a theorem about the full time-dependent pair.
The additional D quadratic block is positive; finite D and metric perturbations
were not evolved in these 3D runs.

## 5. Production, abundance, survival and a named observable at one point

The common physical choice is gX=.1, lambda_phi=.01, lambda_X=.001,
v=1000 GeV, mX=100 GeV; mD0=100 GeV,hD=.01,lambdaD=.001. Leading direct
Standard-Model scalar portals are zero at matching in this EFT truncation.
This is a versioned empirical scale choice; it is not a measured GRUT constant.

**Production preparation:** a primordial nonthermal rotating condensate with
F=.5,S=sqrt(.75),Omega=sqrt(.775), perturbed by fixed fully 3D modes of RMS
amplitude .005, evolves and fragments on a periodic cube of side64. Its net
charge density is fixed to Qbar/volume=.2200852108. No creation of net charge
from an uncharged vacuum is derived. The charge asymmetry is an initial input.

Production endpoint counts / fraction of total signed Q selected at each threshold:

| Grid | 2% of peak F² | 5% of peak F² | 10% of peak F² |
|---|---:|---:|---:|
{production}

These are a connected density network and threshold-dependent concentrations,
not a converged soliton count or an asymptotic cosmological mass function. A
prediction based on identifying every counted component with a particle would
be unjustified. The abundance bound below does not make that identification.

**Abundance input and evolution:** impose YQ=YB, with YB calibrated from
independent deuterium+LUNA BBN, omega_b=.02233+-.00036 (no CMB baryon prior).
With the documented baryon-mass convention, YQ={e['charge_yield_YQ_equals_late_YB']:.7g},
and nQ(today)={e['present_charge_density_cm3']:.7g} per cubic centimetre. Equality
of the hidden yield and late baryon yield is an additional testable hypothesis;
the action does not force it or solve baryogenesis.

**Compatible microphysical epoch:** the initial visible temperature is about
{e['initial']['visible_temperature_GeV']:.5g} GeV with decoupled, much colder
hidden fields. H/mX={e['initial']['H_over_mX']:.5g}; H delta_t during tau240 is
{e['initial']['H_delta_t_for_tau240']:.5g}. The hidden energy is initially only
{e['initial']['rho_hidden_over_radiation']:.5g} of visible radiation. Local flat
production numerics are consistent with this point at that order.

**Survival:** within the adopted exact U(1) EFT, conserved visible entropy and
weak gravity, comoving net hidden charge persists. Fission or evaporation into
charged X carriers does not erase it. The classical bound E>=mu_star*mX*abs(Q),
mu_star=.5623413252, applies on an instantaneous flat FRW slice as well as locally.
It requires a present hidden energy density omega_hidden >= {e['omega_hidden_floor']:.7f},
regardless of the fraction in particular composites. Individual quantum
lifetimes, late charge violation, entropy dilution and strong collapse are not
calculated or included in this preparation.

**Primary named observable:** choose flat FRW, nonnegative other densities and
rho_Lambda=0 for the conservative point. The Friedmann constraint implies

H0 >= 100 sqrt(omega_b + omega_hidden_floor) =
**{geo['H0_predicted_lower_bound']:.2f} km/s/Mpc**.

This bound uses no late hidden equation of state. With three-sigma-low BBN
density and a 1% baryon-mass conversion allowance, it still gives
**{geo['H0_conservative_lower_bound']:.2f} km/s/Mpc**. Positive Lambda strengthens it.
Held-out geometric megamaser distances and velocities give **73.9+-3.0**,
with three-sigma upper edge **82.9**. The specified preparation is incompatible.
This is a comparison with a published inferred measurement; no new raw maser
likelihood or Gaussian exclusion significance is claimed.

If one additionally assumes a late cold collisionless relic, the same floor is
at least {s['empirical_cold_matter']['nominal_minimum_overproduction_factor']:.3f}
times Planck's omega_c=.1200+-.0012. That secondary comparison needs the late
cold-matter interpretation, which has not been derived from the production
network. The geometric H0 exclusion avoids that additional requirement.

The published targets were held out from parameter calibration, not blinded
to the analyst. The failing 100 GeV/YQ=YB point is retained. Changing mX, YQ
or cosmology is a new hypothesis. The action leaves those quantities free;
there is no unconditional successful abundance prediction or detection here.
For a viable empirical claim, a source/calibration must fix the yield and scale
independently, and a cosmological relic/observable calculation must then survive
new data without per-target refitting.

Primary sources: [Mossa et al., Nature 2020](https://www.pure.ed.ac.uk/ws/portalfiles/portal/166986970/mossa_nature_R2.pdf)
(BBN calibration); [Fixsen 2009](https://arxiv.org/abs/0911.1955) (CMB temperature);
[Pesce et al., ApJL 2020](https://arxiv.org/abs/2001.09213) (geometric H0);
[Planck 2018 VI](https://arxiv.org/html/1807.06209v4) (secondary cold-matter constraint);
[2022 CODATA](https://physics.nist.gov/cuu/Constants/Table/allascii.txt) (constants).
The exact values, roles and assumptions are in the calibration/held-out JSON files.

## 6. Reproduction and remaining boundary

After a clean clone of this branch:

```bash
python archive_io.py unpack
python analyze.py
python verify_packet.py
python empirical.py compare
python empirical_gravity.py
python make_figures.py
python report.py
```

Run those commands inside this directory. Re-evolving uses a separate fresh
copy before unpacking (existing raw runs are intentionally not overwritten):
`python run.py run`, `python run_transverse.py`, `node large_domain_modes.cjs`,
and `python spectral.py`. Set GRUT_BUILD_DIR to a writable location if needed.
The raw archives retain exact diagnostic histories, full-grid endpoint component
lists and subsampled float32 spatial volumes. Full evolving field histories
are not retained. `RAW_MANIFEST.json` verifies the archived bytes.

The finite-object conclusion is now substantially stronger than a bulk-frequency
match. Its scope is a declared perturbation family and finite interval, supported
by isolated-branch constrained-energy evidence. The original raw m=1 gain remains
numerically inadequate as an internal-instability rate. The empirical work has
produced a falsifiable conditional exclusion, while primordial yield generation,
an enduring composite population and positive empirical validation remain open.
Unpushed Claude outputs were unavailable and have not been incorporated or certified.
'''
 (P/'REPORT.md').write_text(text)
 print('Wrote REPORT.md')
if __name__=='__main__':main()
