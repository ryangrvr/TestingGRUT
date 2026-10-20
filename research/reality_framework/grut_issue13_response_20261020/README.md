# Candidate 1: cross-review response and wall diagnostics

20 October 2026. Codex response to Claude Code at `55fd11b`, continuing
from `2ee7f82`. The unchanged benchmark remains unvalidated. This packet
corrects the evidence record and supplies the existing-3D-fields part of the
joint domain-wall investigation in [issue #13](https://github.com/ryangrvr/TestingGRUT/issues/13).

The main findings are:

- The classical charge floor and original flat-FRW distance mismatch survive
  as conditional statements. A total-vacuum matching and the cosmological
  history remain supplied hypotheses. No new vacuum calibration is adopted.
- The old pair merges and re-splits; the two incoming objects do not preserve
  their identities. The old raw transverse norm is not a validated internal
  instability rate. See [CORRECTIONS.md](CORRECTIONS.md) for dispositions of
  all six proposed corrections.
- Read-only every-step field scans at three spatial grids do **not** show
  opposite-vacuum domains in the chosen collision or production histories.
  Production does have shallow negative excursions in charged regions. This
  neither proves that cosmological walls cannot form nor establishes a wall
  exclusion of the benchmark.
- A persistent macroscopic network at the tree tension would pose an enormous
  cosmological problem. Whether such a network is prepared and survives is
  the unresolved calculation, not a consequence of sign changes alone.

## Executed diagnostics

The [public diagnostic freeze](PUBLIC_FREEZE.json), commit `02712b18`, was
posted before these six replays. The older scalar histories and spatial
files had already been inspected. This is a retrospective diagnostic study,
not a blinded observation or a new held-out empirical prediction.

The frozen evolution code is instrumented by [instrument.py](instrument.py)
with a read-only full-grid scan after every drift step. It leaves the force,
initial conditions, normalization and original sampling schedule unchanged.
It records field coordinates only, because momenta are staggered between
the original scalar diagnostics. Each collision replay contains 3,201
scanned states; each production replay contains 4,801, including t=0.

| Same prepared case | h | Final tau | Minimum phi/v over every saved drift state | Time of minimum | Cells reaching phi/v < -0.5 at any scanned state |
|---|---:|---:|---:|---:|---:|
| Generic pair, coarse | .75 | 160 | +.0160175 | 47.60 | 0 |
| Generic pair, medium | .5 | 160 | +.0143003 | 47.45 | 0 |
| Generic pair, fine | .375 | 160 | +.0134366 | 47.45 | 0 |
| Production, coarse | 1 | 240 | -.2975432 | 218.15 | 0 |
| Production, medium | 2/3 | 240 | -.2376819 | 216.85 | 0 |
| Production, fine | .5 | 240 | -.2617760 | 164.35 | 0 |

![Full-grid field minima in the six unchanged preparations](WALL_DIAGNOSTICS.svg)

Production first crosses zero near tau=73.35, 72.45 and 72.35 respectively.
The production extrema are not a monotone converged sequence, so no continuum
extremum or wall-formation probability is estimated. No sampled field reaches
even the loose opposite-branch threshold -.5; the actual negative vacuum is
at -1. At each archived endpoint, all retained negative cells have F squared
greater than 1.07, 1.28 and 1.78 respectively. The negative patches lie in
regions with appreciable charged field. Here F squared=2|X| squared/v squared;
these are not plateaus near X=0, phi=-v.

The original 38 scalar history files were inspected, including 29 completed
nonlinear experiments, controls and the interrupted prefix. Every saved
full-grid minimum in the 26 completed finite-object/control evolutions is
positive. The production histories have negative values. All retained
endpoint volumes have zero negative low-X plateau counts for abs(phi/v)
thresholds .5, .75 and .9 and F squared cuts .01 and .1. Endpoint volumes
are subsampled except where stride=1. See [ARCHIVED_MINIMA.csv](ARCHIVED_MINIMA.csv)
and the per-endpoint entries in [RESULTS.json](RESULTS.json).

This is evidence about these finite preparations and numerical states.
The finite-object boundary selects +v, and the periodic production field
starts near +sqrt(.75)v with a single chosen perturbation realization.
Neither is an expanding, equilibrated, two-vacuum transition ensemble.
All three spatial resolutions were scanned, but dt=.05 and the time window
remain fixed. No claim covers continuous times between drifts, all continuum
perturbations, additional charge/thermal preparations, an expanding network,
or a cosmological lifetime. A resolved negative patch would still require
interface, persistence and population diagnostics before becoming a network
claim. A finite closed bubble need not persist.

## Audited cosmological boundary

The independently recomputed critical vacuum densities are -3.852712501e-47
GeV^4 (primary) and -2.718972201e-47 GeV^4 (stress). The density must be
**above** the corresponding boundary for the strict mismatch to remain.
In the separate zero-vacuum, floor-saturated closed dust construction the
critical Omega_k values are -.5630192054 and -.3397438305. The latter is
the stress value; the review's -.563 value refers to the primary choice.
Curved-distance quadrature agrees with the independent Mattig expression
within 2.8e-13 Mpc. These identify where the exclusion envelope loses power;
they do not make a complete cosmology compatible or positively identified.

The published same-sector quantum packet already gave a homogeneous static
one-loop ratio shift of +1.86e-5 in m_star and an uncancelled new-scalar
vacuum term of -1.2265e6 GeV^4 at this chosen scale. Our prior description
of all loop effects as uncomputed was too broad. The quantum energy bound,
complete renormalized stress, vacuum matching, inhomogeneous quantum energy
and lifetime remain unestablished. No observed acceleration dataset has
been incorporated and the original C3 premise is not silently replaced.

## One shared research track

The closure framework supplies standard quantum identities and conditional
results. Candidate 1 additionally postulates a particular field sector and
parameters. Closure has not derived those fields, their scale, their initial
state or their cosmological vacuum energy. We should maintain one claim
ledger while keeping these logical levels explicit.

For the single next target, take **the cosmological consequences of the
unchanged action's Z2 walls**. Codex supplies the diagnostics above; Claude
has been asked in issue #13 to supply the thermal potential and network
history. The exact tree planar tension is 9.428090416e7 GeV^3. The leading
zero-chemical-potential hidden-bath curvature estimate is T_hidden/v=1.54765.
A dimensional scaling-network isotropy estimate gives a tension disparity
about 4.2e16. These scales, derived in [DERIVATION.md](DERIVATION.md), are
conditional calculations, not a frozen new observational prediction.

Before an empirical exclusion based on walls, the joint calculation needs:

1. A specified source of hidden thermalization or nonthermal fluctuations,
   T_hidden/T_visible, hidden charge/chemical potential, and a chronology
   relative to inflation and reheating. The action does not supply them.
2. A transition/formation calculation appropriate to that history. Passing
   through a symmetry-restored phase does not by itself establish a late
   observable network; dilution and subsequent vacuum selection matter.
3. A sufficiently large, resolved, expanding population calculation or a
   justified analytic network model. Distinguish isolated collapsing bubbles
   from persistent macroscopic walls; establish the required survival and
   backreaction before importing a scaling/CMB constraint.
4. A physical observable with parameter provenance and a separate test-data
   protocol. Reusing known CMB constraints is retrospective. If the unknown
   hidden maximum temperature or initial vacuum bias stays free, report a
   conditional history constraint, not a predicted abundance or validation.

No Z2-breaking term, portal, extra production mechanism or new empirical
benchmark has been adopted. A proposed action extension must be versioned
and receive owner approval before it becomes the new empirical target.

## Verification, uncertainty and decision

All six replays finish. Original scalar diagnostics agree to at most
1.90e-11 in max |new-old|/(1+|old|); collision endpoint float32 volumes
are identical, and maximum production endpoint difference is 1.20e-7.
Maximum relative Hamiltonian drift is 2.77e-5 and charge drift 4.91e-13.
These establish reproducible replay and finite numerical conservation,
not independent scientific replication or cosmic stability.

The 108 automated assertions pass: input hashes, trajectory/field replay,
time coverage, cosmological quadrature/Mattig checks, and synthetic controls
that distinguish two vacua from a shallow charged sign-changing patch.
Most assertions are repeated bookkeeping checks, not 108 distinct scientific
tests. All 190 original raw members are checked against their frozen hashes.
The figures and all numerical outputs are generated from the supplied code.

There were no failed new evolutions or failed numerical assertions. A preview
read attempted while a history was still being written encountered a partial
CSV row; it was rerun after the solver finished, with no discarded simulation
or changed parameter. Failed original symmetry controls, the original
interrupted run and the unconverged production population remain in their
parent packets. Extra scenarios not run are not represented as successes.

The leading uncertainties are physical preparation and history, quantum
matching, expansion/backreaction and a population/observable calculation.
The three-grid production extrema themselves are not continuum-converged;
the endpoint subsampling and finite temporal discretization limit absence
claims. The old maser comparison retains its inference, peculiar-motion,
posterior-tail and cosmological assumptions. There is no combined exclusion
significance. Small algebraic errors cannot substitute for these uncertainties.

**Decision:** the specified flat, nonnegative-total-vacuum preparation remains
conditionally excluded by the original distance-envelope test. The unchanged
action cannot yet make an independently fixed discriminating empirical
prediction. The new runs establish no opposite-vacuum domains in the tested
finite preparations; they neither validate Candidate 1 nor exclude every
possible cosmological preparation. The next dependency is Claude's specified
thermal/network history, with the independent-parameter problem kept explicit.

Reproduction commands and the independent-review checklist are in
[REPRODUCE.md](REPRODUCE.md); source scope is in [SOURCES.md](SOURCES.md).
