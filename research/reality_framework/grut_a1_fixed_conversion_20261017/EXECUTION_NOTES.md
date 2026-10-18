# Execution and reproducibility disclosure

Baseline: 643ec7f9af41ec5d8d7816a5dcddbe24de2de7a7 (P1).
Charter/uncertainty freeze: ea45722e0a857cbbb74bac7cc2d8b15ad634d325.
Analytic and control-source freeze: 523f051cfcbc93bed95fe8095fd15593768a87d6.

The mathematical route was developed before the charter, explicitly disclosed
there. The protocol, parameters, thresholds, error method and budget were
frozen before controls. This is not blind preregistration of an unknown result.
The first primary execution passed; no physical variant, implementation
repair, budget relaxation or failed run occurred in this cycle.

Primary run: Python with OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1;
wall 0.987216965 s, peak resident memory 121.40625 MiB. The artifact byte
count at run start was 46,470; final packet storage is also below 100 MiB.
The optional verification replay is charged to the same 60-second total
small-control allowance. No paid compute, service or indefinite process is used.

The verification replay passed in 1.027149438 s, without replacing the
archived observations. Combined primary/replay numerical time was
2.014366403 s, below the original total ceiling. No failed execution or
physical revision was needed.

The path implementation retains every charge/field configuration in the exact
resource sector. The one-cell matrix and six-site kinetic norm are independently
written finite controls; all other code ancestry is disclosed. Full geometry
uses the archived P1 builder read-only with a two-mode diagnostic receiver,
then an exact invariant-resource block. All reverse transitions, storage-error
states and receiver occupations in that block remain. It is not an approximation
to the physical receiver or an N>=2^18 simulation.

Neutral relative blocks use cutoff 32, and finite physical paths have all
separations at their stated sizes. These check formulas; they do not prove
the complete virtual bounds or all-N statements. The latter rest on explicit
analysis, pending independent proof review. No residual-norm clipping,
outward rounding, automatic theorem prover, external reviewer or experiment
is claimed. Nine baseline source pins were checked before and after execution.

`python verify.py` checks hashes and recorded gates. Optional
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py --run-checks`
replays controls without overwriting CHECKS.json. It verifies reproducibility,
not independent mathematical validity. No background execution continues
after this cycle's decision.
