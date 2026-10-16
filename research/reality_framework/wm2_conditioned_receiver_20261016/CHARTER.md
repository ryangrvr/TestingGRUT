# WM2-S6: history-conditioned receiver transitions and residual control

Freeze before S6 numerical controls and trajectories. Parent S5:
98f1e861c3509c67a06032b7742c295069a9a273. Preserve every prior packet,
frozen CP-12, SCOUT-1, canonical documents, and the two S5 error ledgers.
This is a numerical construction for the unchanged S4 law, not S4 plus
a new spectral-gap claim. No subagents or external validation are used.

## First deliverable

Retain every ordered detuning-type word, original root weight, reverse
transition, and exact receiver charge. Permit each word u to use its
own orthonormal receiver frame W_u with every physical frequency retained
in its columns. Keep all k=|u| photon occupations within that frame.
There is no extra photon cutoff or radial/subset reduction.

Derive the complete append block from the rectangular Fock overlap
Gamma_k(W_v^dagger W_u), and use its adjoint for removal. Include
interference between all contributions arriving at the same history.
Derive a physical full-state error account for projected emission,
receiver-mode projections, physical-H Strang splitting, polynomial
exponentials, normalization and frame changes. An omitted component's
future returns must remain paid. Convergence by eye is insufficient.

Use a state-dependent Duhamel residual, not an unproved operator-norm
identification of different history frames. For B_P=PBP and R=(1-P)BP,
bound the first three residual Taylor coefficients and the full tail.
This remains valid for adaptively chosen numerical frames.

## Validation before size claims

Check random complex, unequal-rank Fock overlaps, adjoint append/pop,
and the complete residual against independently lifted physical matrices.
Include destructive/constructive interference controls. Check the common
frame limit against S5, and full propagated wavefunctions at N=3,4 with
small receivers. Compare a physical 576-mode N=3 receiver and S5's
archived controlled formation bounds where feasible. A code agreement
does not independently validate the mathematics or enclose roundoff.

Implement Fock overlaps both directly and through one-particle SVD plus
exact bosonic rotations; use direct blocks for <=65536 matrix entries,
otherwise matrix-free SVD. Record coefficient, frame and overlap costs.
Neither efficient rank growth nor polynomial history complexity is assumed.

## Frozen physical tests

Use the original neutral uniform-root preparation, receiver/probe vacuum,
L_N=64N^2, common g=.9, and one observation at T_N=4N^2. The macroscopic
event remains m>=max(2,ceil(N/2)); full assembly is separately counted.
Primary families: balanced +/-1/4 storage offsets, N=4,5,6, common
detuning 0 or 1/2. Keep root weights equal to original type counts/N.
Reuse only S5's proved sign identity and record every attempted component.
These finite families do not cover all independent disorder patterns.

## Frozen numerical construction and budgets

At each step evolve every frame with the full physical free receiver
unitary. Enlarge each nonvacuum word frame by contact; optionally borrow
one leading missing natural-orbital direction from its primitive parent
and children, weighted by their field density and coupling. This is a
numerical basis choice, not a physical edge or selected geometry.
Compress only by explicit word-conditioned orbital-vacuum projections.
All such projections act on the coherent global state and retain its
global normalization cost. Do not normalize each history independently.

Target state-error envelope .60: physical-H Strang <=.15, accumulated
orbital projections <=.20, projected-emission residual integral <=.20,
all other recorded terms <=.05. Use degree-8 projected exponentials,
three residual coefficients, rank ceiling 10 per history, total
coefficient ceiling 2,000,000 and cached overlap ceiling 8,000,000 complex
entries. Bound each instance by 180 wall seconds checked at step
boundaries. A one-step calculation may overrun that check; report it.
Do not silently relax any ceiling after observing formation.

Allocate unspent projection/residual budgets over remaining steps.
At a rank ceiling a projection may spend already-declared remaining
projection budget. Stop at rank, memory, time or exhausted error budget
and retain the intermediate state as UNRESOLVED at the target time.
No skipped residual may be declared zero. Tiny contact/frame thresholds
and observed normalization changes must be charged explicitly.

Use the maximum of zero, (sqrt(p_hat)-e)_+^2 and
1-[sqrt(1-p_hat)+e]^2. A complete finite-family pass requires probability
lower bound >=1/2 and all numerical gates at the original root weights.
Unknown components contribute zero lower bound without renormalization.
Ordinary double precision is not interval certification.

## Scientific scope

The first success is a correct, tested history-conditioned transition and
error calculus. Computational improvement and six-event completion are
separate empirical questions about the solver. A solver stop is neither
formation failure nor an all-size no-go. Derive actual cost dependence;
do not infer bounded active rank or a nonshrinking disorder interval.
Only earned formation bounds may be joined to S4, with every accepted
record and its actual receiver retained. Dimension, charged-matter
universality, a spatial phase and distinctive predictions remain open.
