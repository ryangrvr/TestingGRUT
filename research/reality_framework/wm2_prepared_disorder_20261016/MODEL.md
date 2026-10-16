# S4 model and preparation account

S4 is S3 with independently assigned storage energies. The transition
grammar, receiver, normalizations and particle species are unchanged.
This is a separately declared disorder extension, not a change to CP-12.

For events u!=v let

    X_uv=b_uv^+ h_v^+ h_u^- f_v^-.

Bond bits b, head bits h and free-event bits f generate directed rooted
words: append a free event at the single head, or apply the adjoint pop.
Storage q_i is a bit, C_i=f_i+q_i-1, n_i=c_i^dagger c_i is the neutral
probe number and F_uv=n_u+n_v. The probe is either vacuum or in its
one-particle component. Receiver modes a_r are bosonic, r=0,...,L_N-1.

The complete law is

\[
\begin{aligned}
H={}&16\sum_i f_i+\sum_i(12+\delta_i)q_i+64\sum_i C_i^2+4n_R\\
&-\sum_{r=0}^{L_N-2}(a_{r+1}^\dagger a_r+\mathrm{h.c.})\\
&-\frac gN\sum_{u\ne v}(X_{uv}q_v^+a_0^\dagger+\mathrm{h.c.})\\
&-\frac1{4\sqrt N}\sum_{u\ne v}F_{uv}(X_{uv}+X_{uv}^\dagger)\\
&-\frac14\sum_{u\ne v}b_{uv}(c_v^\dagger c_u+\mathrm{h.c.}).
\end{aligned}
\]

L_N=64N^2, g in [0.9,1.1], delta_i=bar_delta+epsilon_i,
|bar_delta|<=1/2, |epsilon_i|<=1/4. No terms are deleted in the physical
calculation. The on-event detuning assignment is supplied; no event
coordinate, adjacency, pair affinity or ordering is supplied by it.
The internal receiver chain and its spectrum [2,6] are supplied resources.
Receiver mode distance is not the event metric.

## Initial state and energy

Tensor the frozen original CP-12 neutral ray with a uniform-root
bondless seed. In each root r component, h_r=q_r=1, f_r=0; all other
events have f=1, q=h=0. All bonds, receiver modes and the probe are
vacuum. CP-12 is a decoupled benchmark factor; it supplies no energy.

The conserved charge sum q_i-n_R=1 bounds n_R<=N-1. This is the only
excitation cap used, including in the sparse calculations. Within root
r the initial energy mean is 16N-4+delta_r and the energy variance is
g^2(N-1)/N^2. The original uniform-root superposition has the additional
root-energy variance Var_r(delta_r)<=1/16. It is bounded rather than
extensive. Computations subtract the conserved root energy separately
in each orthogonal root component; population observables are unaffected.

With the probe in vacuum, q equals word membership exactly, even under
independent disorder. The receiver has m-1 quanta for a word of length m.
Energy conservation accommodates each event's mismatch in the explicit
receiver; it does not imply that assembly succeeds at large N.

## Measurement and probe

At the declared formation observation time, record the complete root
and word, accept every m>=max(2,ceil(N/2)), and retain the actual
conditional receiver state. The record is a supplied projective
operation. It dephases the emission interaction, with mean work
-<H_emission>, magnitude at most g. Its entropy and device requirements
are not zero. The saved finite-size evidence includes record entropy.

Inject (|w_0>+|w_1>)/sqrt(2), a neutral probe on an occupied bond. This
preparation is an input, not a spontaneously derived matter species.
Its number variance at w_0 is 1/4 and its mean hopping energy is -J.
The preparation apparatus must account for that energy exchange. Apply
exp(-i theta n_(w_0)) and measure source-resolved number response under
the complete H, summing all final field/storage/geometry records.

The analytic interior theorem is uniform over every accepted word and
every admissible receiver state. Numerical conditional probes use the
actual field formed by their explicitly stated small finite receiver,
with no reset. They are distinct from the primary L_N formation runs.

The path grammar supplies the geometry class and its one-dimensional
ball growth. No branching or dimension-selection result is added.
Original charged CP-12 matter universality, gravity, equilibrium space,
observer mechanics and empirical predictions remain unestablished.
