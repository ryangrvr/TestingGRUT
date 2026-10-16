# WM2-S3: a finite-receiver stability interval, with a spectral limitation

S3 removes the uniform-detuning formation obstruction in a declared
finite, closed model. It proves appreciable macroscopic assembly and
full-H operational transport throughout a nonshrinking parameter interval.
It does **not** retain S2's isolated real-energy consistency band:
the added receiver permits explicit coupled complement resonances.

This is a restricted constructive result. It does not establish a
generic spatial phase, dynamically selected dimension, common geometry
for original CP-12 matter, gravity or a distinctive empirical prediction.

## The mechanism and its inputs

S2 had nowhere to put a detuning mismatch. S3 adds a shared energy
receiver whose one-quantum spectrum is fixed in [2,6]. Each head-only
attachment releases binding energy into that receiver, with every
reverse transition retained. The released energy 4-delta lies inside
the receiver band for all tested detunings.

The Hamiltonian, contact operator and every initial register are given
in MODEL.md. Storage energy changes to 12+delta and assisted coefficient
changes to g/N; these are explicit new-law choices. The bare probe
backreaction and matter hopping remain. The receiver has
L_N=64 N^2 internal modes, initially vacuum, and exact capacity N-1
quanta from the conserved charge sum q-n_R=1.

All event pairs have identical primitive coefficients and no event
coordinates or pair affinities are supplied. The receiver's internal
spectral chain is a supplied auxiliary resource with no map to event
labels. The one-dimensional event class still comes from the single-end
grammar. Bosonic statistics, band shape, normalization, mode count,
probe species and seed are inputs.

## Fixed-interval formation theorem

The class was frozen before S3 controls:

- Uniform on-event storage detuning: -1/2 <= delta <= 1/2.
- Common emission coefficient: 9/10 <= g <= 11/10.
- Same initial vector, receiver, mode-count rule and observation window
  for the entire class; no parameter-dependent retuning.

For every N>=2 and every time in [4N^2,8N^2], the actual finite model obeys

\[
\boxed{\Pr[m\ge\max(2,\lceil N/2\rceil)]\ge\frac12.}
\]

The probability accepts every successful macroscopic word and root.
No favorable observation instant or reservoir packet is selected.
The interval is not inferred from a grid: MATHEMATICS.md proves it
continuously using a uniform spectral-density variation bound.

The vacuum dynamics admit an exact quadratic reduction. Uniform
rooted-word extensions act as annihilation of a bookkeeping pump
occupation. With N-1 initial bookkeeping quanta and one shared receiver,
the complete attachment distribution is binomial, governed by one
explicit Jacobi survival amplitude. The pump is an isometric
representation of free events, not an extra prepared resource.

A half-line comparison gives amplitude <=1000/(729 pi)<0.437 throughout
the window. Returning to the finite receiver changes that amplitude by
at most 2 exp(-120 N^2). The resulting conservative formation bound
is one-half. No master equation, dissipative rate, Markov limit or
irreversible bath is used.

The finite-resource numerical controls give:

| N | Actual finite receiver modes | Largest sampled pump survival | Smallest sampled macroscopic probability |
|---:|---:|---:|---:|
| 2 | 256 | 1.202511e-3 | 0.998797489 |
| 3 | 576 | 1.572319e-3 | 0.999997528 |
| 4 | 1024 | 1.703420e-3 | >0.999999 |

These floating-point samples support the reduction; the all-N theorem
does not rely on extrapolating them. Full-N connectivity at this window
and an all-time finite-system guarantee are not claimed.

## Exact escape from the previous detuning inequality

The old macroscopic compression had cost delta m with bounded kinetic
width. The new receiver has relative energy span [-2k,2k] at k
attachments, which accommodates delta k for every detuning in the
frozen interval. Its initial residual width is only
g sqrt(N-1)/N, so formation is not being purchased with an unexplained
extensive initial energy uncertainty.

The physical initial energy is 16N-4+delta. At mean attachment count k,

\[
16(N-1)=16(N-1-k)+(12+\delta)k
+\langle H_R\rangle+\langle H_{\rm int}\rangle.
\]

The receiver can hold at most 6(N-1) energy units, starts in vacuum and
stores the released binding energy. The maximum complete-ledger residual
in the archived finite controls is 6.409e-15.
The mode count grows as N^2 and the observation time as N^2, openly
replacing S2's tuned O(N log N) protocol with a larger resource budget.

Geometry recording has entropy
H(p_k)+sum p_k log2[(N-1)!/(N-1-k)!]+log2 N, growing as N log N in
the formed window. Its mean projective work is -<H_int>, bounded by
g. Probe injection is a separate preparation. No zero-cost physical
device, reservoir reset or irreversible lock is claimed.

## Complete operational protection survives

The prefix proof remains valid with every storage pattern and receiver
configuration included. Between originally attached events, the minimum
protected primitive route still equals their path separation. Every
P-to-Q-to-P excursion requires at least d+2 geometry/probe operations,
including arbitrary final word size and ordering.

The large receiver Hamiltonian is resummed within word/storage/probe
blocks. The remaining interaction has size-uniform norm at most 2.1,
including bosonic enhancement and all graph branches. Thus

\[
\|R_j e^{-itH}P_{i,d;j}\|
\le {\rm tail}_d(2.1|t|)
\le e^{2.1e^\mu|t|-\mu d}.
\]

R_j sums every final receiver, storage and geometry record. This is a
complete full-H propagation front, not a favorable reduced trajectory.
Its constants are uniform throughout the frozen detuning/coupling class.

For an injected adjacent two-site probe, the source variance is 1/4.
The number-phase receiver response has a product-of-tails distance
bound, and nearby response at t=0.005 is at least **0.0010279483**,
uniform in size, formed word and receiver state. The smallest nearby
response in the explicit small controls is approximately 0.00125.
The conditional receiver state in those controls is obtained from
neutral-seed evolution; it is not reset to a favorable packet.

Initial consistency has finite-time mismatch probability <=epsilon^2 t^2.
This is a Duhamel estimate, not S2's all-time leakage floor.
Conditional bulk-prefix failure is bounded by
tail_r(2.7|t|/sqrt(N))^2. A fixed-fraction interior of a formed
macroscopic component persists over probe times proportional to N.
An entire word is not irreversibly locked.

## The complete virtual bound has a resolution parameter

A new distance-weighted resolvent estimate sums every virtual receiver
and geometry path at z=E+i eta:

\[
\|R_j^P\Sigma(E+i\eta)P_{i,d;j}\|
\le\frac1{8\eta}
\left(1+\frac{5\eta}{21}\right)^{-d},
\quad \eta>0.
\]

It uses a weighted configuration-block resolvent, not a small global
norm for the extensive receiver. Sources may be coherent geometry and
receiver packets; receivers contain every consistent word of every size.
The operator bound is verified on 132 small all-record blocks using
complex Gaussian elimination and a Frobenius upper bound.

The eta parameter must remain visible. This result is **not**
S2's uniform real-energy Schur theorem.

## A coupled spectral obstruction is proved

Put all storage bits in the spent state, leave one event free, put the
probe on that free event and put all N-1 receiver quanta into a
far-contact modulated sine packet. The receiver can pay the 80-unit
free/spent defect cost by a small energy change per quantum.

For N>=64 and every allowed detuning there is a normalized Q state with

\[
\|QHQ\Psi_Q\|\le
r_N=\frac{2\pi(N-1)}{32N^2+1}\longrightarrow0.
\]

Its contact occupation is exactly zero. A spectral-projection argument
also proves a nonzero B coupling to at least one Q eigenvalue in
|E|<=2sqrt(N)r_N. These are coupled real poles of the full consistency-
sector Schur operator, inside the seed-energy window.

| N | Upper bound on distance from zero to Q spectrum | Coupled real pole exists within this absolute energy |
|---:|---:|---:|
| 64 | 3.020002e-3 | 4.832003e-2 |
| 128 | 1.521994e-3 | 3.443879e-2 |
| 512 | 3.827461e-4 | 1.732111e-2 |
| 4096 | 4.792520e-5 | 6.134425e-3 |

This rejects the S2-style complement-gap hypothesis and a bounded real
Schur operator throughout that window. It does **not** prove a shortcut in
actual matter propagation, nor divergence of every off-distance block.
The full-H front still holds. Absorbing those levels into a new
energy-defined sector would require a fresh metric and spectral proof.

A direct one-particle packet calculation through 1,048,576 receiver
modes agrees with the product-state residual formula. It is not a
million-mode many-body trajectory or a factorial Q diagonalization.

## What was executed

All archived numerical algorithms ran in the available JavaScript V8
runtime. The local shell disconnected during the campaign. The provided
dependency-free Node driver uses the same functions, but that driver
was not run under Node here. No Python simulation is reported.

- 1,487 complete small register/configuration states; 5,806 reverse checks
  and 4,634 virtual-exit distance checks.
- 36 unreduced vacuum wavefunction comparisons against the exact
  quadratic Fock prediction; maximum discrepancy
  2.177e-15.
- 90 continuum parameter sets and 270 observations through N=32;
  nine quadrature refinements.
- 27 actual L_N finite-receiver parameter sets and 81 observations
  through N=4; largest one-particle matrix dimension 1,025.
- 132 complex Schur blocks, 15 finite phase-pulse identities and eight
  full-H receiver propagation blocks.
- 18 analytic near-Q packet evaluations and two direct mode-packet checks.

The transport/register controls state their smaller L explicitly.
They are not presented as full-probe simulations with L_N. The exact
all-size bounds include that full receiver. Arithmetic is ordinary
double precision, not outward-rounded certification or external validation.

## Acceptance ledger and next boundary

| Gate | S3 result |
|---|---|
| Explicit finite energy receiver and reverse terms | Pass |
| Escape from S2 uniform-detuning obstruction | Pass |
| Nonshrinking formation interval | Pass in the frozen symmetric class |
| Full-H matter response and nearby transport | Pass in the declared probe sector |
| Conditional bulk-prefix persistence | Pass with supplied normalization |
| Isolated real-energy consistency sector | Fail: coupled Q poles |
| Arbitrary independent local perturbations | Not established |
| Irreversible locking or equilibrium phase | Not established |
| Dynamically selected dimension | Not earned; path grammar supplied |
| Common geometry for original CP-12 matter | Not established |
| Distinctive prediction about nature | None |

The next construction problem is an energy-defined operational sector
or a controlled nonequilibrium replacement whose receiver accessibility
is derived. A global spectral gap cannot simply be copied back in.
Keep the single-end grammar and the frozen branching counterexample
while that issue is resolved.
