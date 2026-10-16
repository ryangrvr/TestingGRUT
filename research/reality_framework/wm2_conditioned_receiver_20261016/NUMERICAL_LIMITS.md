# What the evaluated error account does and does not enclose

THEORY.md's identities and inequalities are exact-arithmetic statements
for the complete finite model. The implementation evaluates them with
ordinary NumPy/SciPy double precision. No packet status is a claim of
outward-rounded numerical certification or external mathematical review.

The complete residual is obtained as ||B Vx||^2-||V^dagger B Vx||^2.
It is nonnegative in exact arithmetic. Near cancellation, floating-point
subtraction can be slightly negative. The code records the largest
negative value, rejects a substantial violation, and clips tiny negative
values to zero for the evaluated square root. That is a disclosed
arithmetic convention, not a proof that every true residual is enclosed.
An interval evaluation or a stable rigorously enclosed residual Gram
factorization would be needed to close this arithmetic obligation.

Frame orthogonality is checked periodically. Polar frame changes and
observed global normalization changes are charged. Tiny directions not
added to a numerical frame remain included in the full physical-action
residual; they are not presumed absent from the physical Hamiltonian.
The analytic residual tail pays all higher uncomputed coefficients.

Natural-orbital projections use each history's unnormalized global
amplitude. One global renormalization is performed, and its full state
change is paid. The numerical procedure never gives a rare favorable
history weight one. Accepted formation must use original root weights
and every final receiver/ordering. Unknown target-time components keep
the trivial zero lower bound without renormalizing the remainder.

The original S6 trajectories were executed before caching their fixed
two-mode rotation matrices. The final code and G1 precompute those
identical matrices once per SVD/Givens plan. Additional controls compare
cached and uncached arithmetic bit for bit. This changes computation
reuse, not mathematical operations on a coefficient vector, frames,
physical parameters, error budgets or event definition. Wall-time
completion comparisons are therefore not clean measurements of the
residual-allocation policy alone. Reproductions can also stop at
different steps on different hardware.

Primary trajectories are run sequentially. Large verification work is
not deliberately run alongside a wall-limited primary instance. The
reference controls use independently constructed complete physical
matrices or the known exact quadratic uniform solution. These tests
establish implementation agreement in tested cases; they do not prove
all-size formation, uniformly sufficient ranks or physical acceptance.

S5's stronger finite formation results stay in their own unchanged
record. Replacing its solver, choosing a new numerical frame, or finding
a favorable intermediate population cannot upgrade an unresolved S6
target-time result. An exact transition theorem, a useful solver and a
successful physical formation theorem are separate accomplishments.

The prototype removes at most one orbital per history per compression
call. Its frame need not be minimal even when that history's current
one-body density has small rank. In particular a pure one-photon field
has rank-one density, but may retain unused orbitals as a workspace for
subsequent Hamiltonian action. A rank-ceiling stop therefore does not
prove that the physical receiver packet needs more than ten orbitals.
The action-density and exact Krylov module is a separately tested
construction; it was not used to rewrite either primary run.
