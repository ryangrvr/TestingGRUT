# S9 frozen decisions and stops

Freeze this file and PHYSICAL_TARGET.md before numerical controls or large
trajectories. Own this folder/branch only; preserve the S8 boundary and every
earlier failed route. No outside messages or subagents.

## A: one production sweep protocol; no fallback

Use S8-G1 physical integration, finite-receiver cut, operator audits and
global charge allocation, maximum total Schmidt rank 96. Keep complete
coherent history, exact charge, reverse operations and all retained receiver
occupations. Implement the direct in-range residual in SWEEP_CERTIFICATE.md.

One orthogonal-loss batch is one monotone field sweep. Reset on reversal,
boundary gate, normalization and every physical step. Add completed batches
linearly; never square-sum through returns. Boundary gates are singleton
batches. Proposed extra tiny drops use half the remaining batch allowance
divided by sqrt(maximum gates in this batch); hard rank choices still pay
their actual positive loss. Reserve remaining allowance over 8 batches per
remaining physical step. A pending batch counts against the stop gate.

Same .60 total state envelope as S8: Strang <=.20, orthogonal projection
<=.35, field factorization <=.02, all other costs <=.03, now including the
linear numerical in-range defect. Keep these quantities separate. Require
summed site Gram bounds <1e-6. Use the declared S8 IEEE/BLAS model; this
does not certify compiler/hardware correctness or external validation.

Before primaries: exact/direct small-matrix projector controls, a return
counterexample, full small physical wavefunctions and original N=3,L=576,
T=36 full-state reference. A failed mathematical or implementation check
stops A; do not fix a failed primary by changing budgets.

Then sequentially run both N=5 components (2,2) and (3,1), T=100,
physical L=1600, root weights 3/5 and 2/5. One fresh worker per component;
**1200 wall seconds**, including imports/setup/operator audits, soft peak
RSS 2048 MiB. This new predeclared limit pays for direct projection checks
and differs from S8's 900 and S7's 600. No equal-budget speed claim.
Measure whole worker wall/CPU/RSS; retain mid-sweep stops without event
intervals. Both components must complete within every gate, and their own
weighted macro lower bound must be >=1/2. No cross-method mixture passes.

If either fails, finish recording the other predeclared component, then
**close this numerical branch**. No rank increase, new allocation, solver
replacement, fallback or N=6 trajectory follows within S9.

## B: exact history-retaining formation analysis

Use only the class and time-sampling protocol in PHYSICAL_TARGET.md. First
derive finite physical-receiver correlation bounds and the complete coherent
formation-current identity. Count noncommuting append operators and both
creation and absorption commutators. A Markov-like drift is a comparison
only if its finite-time error is bounded back to the closed Hamiltonian.

Success: a uniform positive formation bound for explicit N>=N_0, or a
counterexample within the frozen physical class. A conditional theorem
must state its unproved condition and does not pass this gate. A bound
that is vacuous at accessible sizes must display that boundary. Small
numerics check identities/constants, not all-size proofs.

The initial analytic route suggests contact occupation can be controlled
without collapsing histories. This is a **preliminary insight**, not a
completed formation theorem at this freeze. If the route fails, report
the exact obstruction; do not silently assume independent emissions.

## Three ledgers and closure

Maintain PHYSICS_LEDGER.md, MATHEMATICS_LEDGER.md and METHODS_LEDGER.md.
Publish the finite successes/stops, proof/check status, supplied-input map
and integrated protocol if earned. A proof or finite recovery is not a
distinctive prediction about nature. No unlimited solver campaign is authorized.
