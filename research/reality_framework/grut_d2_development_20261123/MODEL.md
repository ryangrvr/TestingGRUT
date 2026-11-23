# GRUT-D2: autonomous ordering and an internal reservoir record

23 November 2026. Exact finite model using adopted quantum mechanics.

The proposed universe contains a connected open chain of `n` physical spins,
`L` reservoir spins, and an inert diagnostic reference `R` storing the initial
defect count. The active Hilbert dimension is `2^(n+L)`; including the reference
it is `n*2^(n+L)`. Large calculations use an exact invariant-sector/free-fermion
reduction, not the full exponential matrix. This is a small local interaction
rule with a scalable family, not a proof of minimal Hilbert dimension.

Set hbar = g = 1 and Delta = 4. Physical spin indices are `0,...,n-1`.
Let `m=n-1`, `w_i=(1-Z_i Z_(i+1))/2`, and reservoir occupations
`b_j^dagger b_j=(1-Z_Bj)/2`. Here a lowering operator is
`b=(X+iY)/2=|0><1|`. In wall coordinates write `tau_i^- = |0><1|`.
The total time-independent Hamiltonian is

\[
 H=\Delta\left(\sum_{i=0}^{m-1}w_i+N_B\right)
 +g\sum_{i=0}^{m-2}(\tau_i^+\tau_{i+1}^-+\mathrm{h.c.})
 +g(\tau_{m-1}^-b_0^\dagger+\mathrm{h.c.})
 +g\sum_{j=0}^{L-2}(b_j^\dagger b_{j+1}+\mathrm{h.c.}).
\]

This is **not** a nonlocal interaction in the original spin coordinates.
The first term is an ordinary ferromagnetic Ising energy (plus a constant).
Adjacent wall hopping has physical form

\[
 \frac g2 X_{i+1}(1-Z_i Z_{i+2}),
\]

and the interface has physical form

\[
 \frac g2(X_{n-1}X_{B0}+Z_{n-2}Y_{n-1}Y_{B0}).
\]

The reservoir uses ordinary XX nearest-neighbor exchange. Every term has
support on at most three neighboring spins. The last displayed interface
implements resonant exchange of a domain wall with a reservoir excitation.
The energy preference, resonant interface and prepared reservoir are supplied;
calling them spontaneous design or a newly derived law would be incorrect.

## Preparation and operational quantities

Start the physical system in `I/2^n` and the reservoir in `|0...0>`.
Equivalently, the orientation `s_0` is a fair bit and the walls are independent
fair bits. The reference is initialized to their count `K`; it is a declared
diagnostic correlation, not dynamically generated memory. Its Hamiltonian is
zero and it does not control the evolution.

Organization is the **unconditional** probability of the family
`00...0` or `11...1`, equivalently `N_wall=0`. Its initial value is
`2^(-m)`. Both orientations remain equally likely, and `Z_0` is conserved.
The global spin flip on the physical chain commutes with H. Absolute
orientation is retained from a microscopic seed; alignment itself is favored
by the chosen interactions. There is no postselection.

The newly generated record is the reservoir's total excitation number `N_B`.
With `Q=N_wall+N_B=K` conserved, the decoder `K_guess=N_B` is correct exactly
when all walls have left the system. It stores a count, not the full original
wall pattern, a developmental program, a self-model, or consciousness. Access
to this collective reservoir observable is supplied. No observer or readout
apparatus emerges in this model.

`S2=sum(i*i for i in 1,...,m)`, `L=64*S2` and the scored window is
`8*S2 <= g*t <= 16*S2`. This intentionally buys a controlled lifetime by
increasing a cold reservoir and the elapsed time. Both grow as `O(n^3)`.
The memory resource has `L` initially pure spins. The total closed system
never exports entropy outside itself and eventually recurs.

## Time, energy and interventions

All interactions coexist in one fixed H. The propagation parameter t remains
the supplied Schrodinger time; this packet does not derive internal time or
supply a physical clock. Diagnostic flips and bath phase kicks use external
timing and coherent control. They are counterfactual controls with their own
work and initialization costs, not part of the autonomous baseline.

The vacuum has zero energy and H is nonnegative for Delta > 2g. Initially
`<H>=Delta*m/2`, `Var(H)=Delta^2*m/4 + g^2*m/2` for the connected fair-bit
preparation. The global active entropy is `n` bits, with `L` bits of bath
entropy capacity prepared. Removing the interface preserves the initial mean
energy but lowers the variance by `g^2/2`. This distinction is reported, not
called an energy-distribution-matched memory test.

## Interpretation gate

This construction replaces a sequence of corrective gates with a fixed local
law, retains both orientations, and supplies a controlled finite window.
It still inserts a preferred alignment relation and a reservoir preparation.
It is a concrete standard-physics mechanism for ordering and a passive record.
The stricter claim of emergence without a supplied target class, active
developmental feedback, any new GRUT prediction, consciousness and cosmic
predetermination remain unestablished. Candidate 1 is a separate program.
