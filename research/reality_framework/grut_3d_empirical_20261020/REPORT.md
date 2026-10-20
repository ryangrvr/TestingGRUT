# Candidate 1: finite objects in 3D and a conditional empirical exclusion

20 October 2026. Experimental research branch; independent review pending.

The tested isolated object and outgoing collision cores survive the declared
finite three-dimensional perturbations through tau160. The isolated branch also
passes a numerical fixed-charge curvature test with an analytic ordering of all
angular harmonics. The old transverse canonical norm cannot be interpreted as
an internal-instability rate: its spatial symmetry controls retain large errors.
Direct full-3D impact experiments separately test the outgoing trajectory response.

One common-scale cosmological preparation has been carried to a conservative
gravitational observable. It fails: the conserved-charge energy floor requires
H0 >= 116.66 km/s/Mpc, incompatible with independent
geometric megamaser data. This is a conditional exclusion, not an observation
of the proposed new fields or a viable empirical theory of everything.

## 1. What was actually computed

The common law is the unchanged [Candidate-1 action](../grut_candidate_1_20261018/GRUT_CANDIDATE_1.md).
The new main campaign contains **29 nonlinear full Cartesian 3D runs**:
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
29 completed runs is 4.91e-13; maximum Hamiltonian drift is 2.76e-05.
The phase tangent is exact to printed roundoff. A centered finite-difference
tangent check reduces its error by a factor 4 when its increment is halved.

## 2. What survives, and how accurately

Retention is defined before the long runs: all sampled tau120..160 must retain
at least 90% of initial Q in R12 core balls, with at least 15% in each of a pair
and disjoint pair balls. A single object additionally requires local E/abs(Q)<1
and one dominant full-grid connected component. R8 charge covariances and
2/5/10%-of-peak connectivity thresholds are retained as separate diagnostics.
**26/26 finite-object runs meet their applicable retention criterion.**

Perturbed isolated object:

| Grid | h | Minimum late R12 charge fraction | Late local E/Q range |
|---|---:|---:|---:|
| coarse | 0.75 | 99.85098% | 0.909894–0.909938 |
| medium | 0.5 | 99.84893% | 0.910952–0.910998 |
| fine | 0.375 | 99.84842% | 0.911321–0.911367 |

Finest outgoing objects (late mean charges are dimensionless Qbar):

| Case | Smaller core mean Q | Larger core mean Q | Minimum joint retained Q | Criterion |
|---|---:|---:|---:|---|
| pair_reference_fine | 259.867491 | 519.202919 | 97.48465% | pass |
| pair_generic_fine | 262.325148 | 513.556276 | 97.06750% | pass |
| pair_negative_fine | 263.677860 | 511.925937 | 96.98103% | pass |

Local laboratory E/Q is reported in SUMMARY.json. It is a sufficient binding
check when below one, not a necessary one for a moving object. A core ball is
not proof that the outgoing state has settled to an exact stationary solution.

Maximum late charge differences, normalized to initial total Q:

| Family | Coarse–medium | Medium–fine | Charge refinement criterion |
|---|---:|---:|---|
| single_reference | 0.00519% | 0.00098% | pass |
| single_generic | 0.00394% | 0.00181% | pass |
| pair_reference | 0.21427% | 0.06406% | pass |
| pair_generic | 0.33654% | 0.10203% | pass |

| Control | Maximum late charge change | Charge criterion | Relative-shift norm change / criterion |
|---|---:|---|---|
| single_box | 0.02229% | pass | not used |
| single_dt | 0.00015% | pass | not used |
| pair_box | 0.06678% | pass | 0.61816% / pass |
| pair_dt | 0.00271% | pass | 0.03364% / pass |

Whole-event rotations by .3 radians rotate both positions and velocities. The
largest late core-charge difference falls from
0.06077%
at h=.5 to 0.04318%
at h=.375. This is a distinct check from shifting the sources while holding
their velocity directions fixed.

These are finite-time, finite-amplitude-family results. At the empirical 100 GeV
point, tau160 is about 1.053e-24 seconds,
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
| 0.500000 | 5.394846 | 136.75577% | 89.31860% | 99.98563% / 99.65729% |
| 0.333333 | 7.970893 | 73.67583% | 47.80780% | 99.99723% / 99.77347% |
| 0.250000 | 9.063794 | 44.27692% | 28.66341% | 99.99872% / 99.80081% |

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
| impact_coarse | -13.726052 | 7.333882 | 0.000567155 |
| impact_medium | -14.114130 | 7.523045 | 0.000145112 |
| impact_fine | -14.241539 | 7.586161 | 0.000136989 |
| impact_half_fine | -14.239909 | 7.581834 | 8.78304e-05 |
| impact_box | -14.096267 | 7.529909 | 0.000157252 |
| impact_dt | -14.111692 | 7.522332 | 0.000160799 |

| Response comparison | Maximum relative centroid-response change | Frozen criterion |
|---|---:|---|
| coarse_medium | 2.74957% | reference difference |
| medium_fine | 0.89463% | pass |
| half_seed | 0.05704% | pass |
| box | 0.22715% | pass |
| dt | 0.03109% | pass |

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
| 50 | 0.1 | -6.49113225e-05 | -0.000170587901 | 0.364348002 | 0.285352436 |
| 50 | 0.05 | -1.62253604e-05 | -4.26265943e-05 | 0.364348691 | 0.285457169 |
| 50 | 0.025 | -4.05618498e-06 | -1.0655393e-05 | 0.364348863 | 0.285483333 |
| 60 | 0.05 | -1.62253605e-05 | -4.26265949e-05 | 0.362970677 | 0.285457169 |

The small negative sampled symmetry eigenvalues vanish as dr squared. The
positive radial levels just above .36 are finite-box scattering levels near
the vacuum threshold, not newly identified internal resonances. Re-solving
the profile gives dQ/dOmega = -4094.028.
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
| coarse | 315 / 85.22459% | 105 / 75.09988% | 29 / 68.84310% |
| medium | 281 / 85.15240% | 119 / 75.05493% | 36 / 68.58766% |
| fine | 248 / 86.11910% | 111 / 75.79647% | 39 / 68.99042% |

These are a connected density network and threshold-dependent concentrations,
not a converged soliton count or an asymptotic cosmological mass function. A
prediction based on identifying every counted component with a particle would
be unjustified. The abundance bound below does not make that identification.

**Abundance input and evolution:** impose YQ=YB, with YB calibrated from
independent deuterium+LUNA BBN, omega_b=.02233+-.00036 (no CMB baryon prior).
With the documented baryon-mass convention, YQ=8.675814e-11,
and nQ(today)=2.508369e-07 per cubic centimetre. Equality
of the hidden yield and late baryon yield is an additional testable hypothesis;
the action does not force it or solve baryogenesis.

**Compatible microphysical epoch:** the initial visible temperature is about
1.7563e+05 GeV with decoupled, much colder
hidden fields. H/mX=4.3336e-10; H delta_t during tau240 is
1.0401e-07. The hidden energy is initially only
6.2192e-14 of visible radiation. Local flat
production numerics are consistent with this point at that order.

**Survival:** within the adopted exact U(1) EFT, conserved visible entropy and
weak gravity, comoving net hidden charge persists. Fission or evaporation into
charged X carriers does not erase it. The classical bound E>=mu_star*mX*abs(Q),
mu_star=.5623413252, applies on an instantaneous flat FRW slice as well as locally.
It requires a present hidden energy density omega_hidden >= 1.3387081,
regardless of the fraction in particular composites. Individual quantum
lifetimes, late charge violation, entropy dilution and strong collapse are not
calculated or included in this preparation.

**Primary named observable:** choose flat FRW, nonnegative other densities and
rho_Lambda=0 for the conservative point. The Friedmann constraint implies

H0 >= 100 sqrt(omega_b + omega_hidden_floor) =
**116.66 km/s/Mpc**.

This bound uses no late hidden equation of state. With three-sigma-low BBN
density and a 1% baryon-mass conversion allowance, it still gives
**113.25 km/s/Mpc**. Positive Lambda strengthens it.
Held-out geometric megamaser distances and velocities give **73.9+-3.0**,
with three-sigma upper edge **82.9**. The specified preparation is incompatible.
This is a comparison with a published inferred measurement; no new raw maser
likelihood or Gaussian exclusion significance is claimed.

If one additionally assumes a late cold collisionless relic, the same floor is
at least 11.156
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
