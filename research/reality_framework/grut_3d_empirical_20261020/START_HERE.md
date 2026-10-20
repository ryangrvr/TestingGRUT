# Durable GRUT research handoff: 3D objects and an empirical exclusion

20 October 2026. Read REPORT.md, SUMMARY.json and VERIFICATION.json before
continuing. This checkpoint supersedes the earlier audit's *next-task list*,
while preserving every older frozen packet byte-for-byte.

## Owner's task and work boundary

The owner asked for two deciding advances: unrestricted 3D finite-object
perturbations with symmetry/convergence controls, and a common-parameter
production/abundance/survival/observable chain compared with held-out data.
Tighter fits to previously tested bulk frequencies do not fulfill either task.

Repository: `ryangrvr/TestingGRUT`.
Isolated research branch: `codex/grut-candidate1-two-object-20261019`.
Parent published commit: `9ae004a06868674487bca679697bb073a10b7c75`.
Use the branch's latest commit for this packet; do not infer current state from
the parent alone. Root AGENT_COORDINATION.md records the completed work.

Claude Code remains the main builder in the owner's other checkout. Its unpushed
tests were unavailable here. Do not claim to have read, incorporated or verified
them. No PR, merge, person-directed message or change to Claude's checkout was
part of this work. Preserve unrelated work and read/claim coordination before edits.

## Common physics

Candidate-1's action is in `../grut_candidate_1_20261018/GRUT_CANDIDATE_1.md`.
It adopts Einstein gravity, quantum rules, ordinary particle content, and adds
a quartic FLS/Q-ball scalar sector. a=1,b=.1,gX=.1 are unchanged. This is not
an established theory of everything or a consciousness result.

The new 3D computations evolve u,w,S with no imposed spatial symmetry, D=0 and
flat local gravity. The physical empirical point is v=1000 GeV,mX=100 GeV,
with a postulated primordial hidden yield YQ equal to the independently BBN-
calibrated late baryon yield. The original v/M_P=1e-5 illustrative scale is not
silently reused for this production epoch.

## What the evidence now supports

- The declared isolated and outgoing collision cores retain localized charge
  under the tested finite 3D shape/phase/velocity/impact perturbations through
  tau160. Core observables have grid, time, larger-box and whole-event rotation
  controls. Read their actual pass/fail values in SUMMARY.json.
- The isolated reference branch has positive fixed-charge radial curvature,
  correctly converging phase/translation zero modes, and analytic angular
  ordering for the remaining harmonics. This supports local energetic stability;
  it is not an interval-certified global-minimum theorem or a quantum lifetime.
- The original linear m=1 amplification about 9.064 at its finest old grid is
  reproduced. Its translation/rotation residuals remain material. The geometric
  response inside each outgoing core is overwhelmingly a local displacement.
  Direct nonlinear 3D opposite-source-offset runs independently measure the
  trajectory deflection, with a halved offset and time/box/grid controls.
- The production run makes a charged network/concentrations. Endpoint component
  counts depend on threshold and resolution; they are not a particle abundance
  or a converged cosmological soliton mass function.
- A complete **conditional exclusion bound** uses the conserved net yield,
  regardless of composite fission/evaporation: the energy-charge floor gives
  omega_hidden>=1.3387. Flat FRW with nonnegative other densities then requires
  H0>=116.66 km/s/Mpc, or >=113.25 under the documented conservative allowance.
  Independent geometric megamaser data give 73.9+-3.0. This preparation fails.
  The secondary Planck omega_c comparison additionally assumes late cold matter.

The published comparison targets were held out from fitting, not blinded to
the analyst. Primordial yield equality is an added hypothesis, not derived
baryogenesis or hidden-charge generation. No viable empirically validated GRUT
point or individual cosmological composite lifetime is established.

## Reproduce and inspect

Run in this directory:

```bash
python archive_io.py unpack
python analyze.py
python verify_packet.py
python empirical.py compare
python empirical_gravity.py
python make_figures.py
python report.py
```

For re-evolution, use a separate fresh copy before unpacking; runners refuse to
overwrite completed raw runs. `run.py run` executes the 23 main cases;
`run_transverse.py` executes six direct-impact cases; `large_domain_modes.cjs`
replays three original linear sectors; `spectral.py` re-solves the spectral
problem. Set GRUT_BUILD_DIR to a writable location if needed. Do not rerun the
development `preflight` command over the frozen PREFLIGHT.json.

The raw archives contain exact scalar histories, full-grid component lists,
and subsampled spatial volumes; complete evolving field histories were not
retained. All local freeze and execution boundaries are explicit. An interrupted
fine run's prefix is retained and compared with its unchanged-source restart.

## Next work must resolve a remaining scientific gap

For the original *global* linear m=1 norm, use a symmetry-consistent spatial
scheme or sufficient refinement/domain control to make generator residuals
small before quoting an instability exponent. Keep physical separate-core
motion distinct from intrinsic shape. Do not infer that every increasing norm
means a disintegrating particle.

For empirical success, independently fix or derive the hidden initial yield and
physical scale, follow an enduring relic distribution and survival history, and
predict an observable without fitting that same target. A different mass/yield
that merely matches a dark-matter density is a fit, not the missing prediction.
Preserve the rejected 100 GeV/YQ=YB point when proposing a new hypothesis.

When Claude results arrive, require source/diff, action/parameters, preparation,
raw run IDs, numerical controls and the exact claimed observable before combining
them. Keep finite-object durability, cosmological population, ordinary detection
and subjective experience as separate questions with their own evidence.
