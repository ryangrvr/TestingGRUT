# A constrained preparation kernel and the remaining phase freedom

This finite example makes a preparation slot concrete. It does not solve
the continuum gravitational constraints or derive the state's parameters.

Take three bosonic modes a,b,c with charges +1,+3,-6. The third represents
a conjugate charged mode, not an identification of the full Sigma field
with only one particle oscillator. With occupation cutoffs, project

\[
|\psi_\delta\rangle=\mathcal N P_0
|\sqrt6,\sqrt2,\sqrt2e^{i\delta}\rangle,
\quad P_0:\ n_a+3n_b-6n_c=0.
\]

The pre-projection coherent state is Gaussian in oscillator coordinates;
the projected state need not be Gaussian or a product of modes. It is
normalizable, positive as a density matrix and exactly gauge neutral.
The cutoffs `(18,8,6)` and `(24,12,9)` leave neutral dimensions31 and58.

Use the explicitly declared illustrative dimensionless Hamiltonian

\[
H=n_a+1.7n_b+2.3n_c
+.025(b^\dagger a^3+\mathrm{h.c.})
+.04(cb^2+\mathrm{h.c.}).
\]

Every transition conserves total charge. These frequencies and occupations
are toy choices; they are not the physical heavy-field mass spectrum or
the cosmological condensate abundance. The light-sector charge is
`L=n_a+3n_b`, with `d<L>/dt=i<[H,L]>`.

For phases `delta=+/-pi/2`, the conjugate wavefunctions have identical
energy, occupation probabilities and all single-mode occupation
entropies, while the torque changes sign. To exhibit zero initial signed
light charge, use the tensor product of the two conjugate copies, assign
opposite physical charges to the second copy, and observe `L_1-L_2`.
It starts at zero but its derivative is nonzero. Exchanging the supplied
phases reverses that derivative without changing the listed initial
scalar properties. The refined torque is `+/-2.71528939`; cutoff change
is .102%, while energy changes .00374%.

Append the same uncorrelated SM factor `rho_SM=diag(.2,.8)` to either
boundary state. Hidden/SM mutual information is zero for both; its largest
numerical residual is `4.45e-16`. This condition also fails to choose a
phase. The projection's internal mode entanglement should not be confused
with entanglement to SM radiation.

## Internal clock and physical observations

Let `U=exp(-.1i H)` and take five clock readings, with `W_t=U^t`. The
isometry and history are

\[
V=\frac1{\sqrt5}\sum_{t=0}^4|t\rangle\otimes W_t,
\quad|\Psi_\delta\rangle=V|\psi_\delta\rangle.
\]

The positive propagation operator is

\[
C=\frac12\sum_{t=0}^3
(|t+1\rangle\otimes1-|t\rangle\otimes U^\dagger)
(\langle t+1|\otimes1-\langle t|\otimes U).
\]

Equivalently each summand penalizes
`psi_(t+1)-U psi_t`; thus `CV=0`. This is an open finite clock with one
time-independent constraint, not a periodic clock or a claim of a cosmic
time endpoint. Let `P=VV^dagger`. For any matter event `Pi_a`,

\[
E_{a|t}=5P(|t\rangle\langle t|\otimes\Pi_a)P.
\]

These effects are positive, sum to P at fixed reading and commute with C.
Their expectation in the physical history is the ordinary conditional
Born probability in `W_t psi_delta`. The code directly checks the event
`n_a>=6` at reading2. It does not construct an autonomous physical reader
for every effect.

A state selector can be written as the positive physical kernel
`B_delta=V |psi_delta><psi_delta| V^dagger`, optionally with the same SM
factor. This provides a constrained preparation rule once delta is given.
It gives no reason to favor one delta, cosmic entropy or inflationary
holder over another. A real diagonal Lorentzian boundary phase preserves
configuration probabilities immediately after insertion, though it can
change their later evolution. An amplitude penalty, density kernel or
Euclidean prescription requires a distinct stated rule. Diffeomorphism
covariance of a boundary expression alone does not solve the continuum
Hamiltonian constraint or define its physical inner product.
