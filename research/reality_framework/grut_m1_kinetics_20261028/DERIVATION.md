# Rates, energy and approximations

All energies and momenta below are in m_X units. The reference core has
Omega=0.88172049195 and Q0=200000. The canonical hXXbar vertex is g_h=2g_X=0.2,
m_h=sqrt(6), and Gamma_h=g_h^2 sqrt(1-4/m_h^2)/(16 pi m_h).
The three-channel radial coupling remains sqrt(2) F S; Argon's FS entry was
a normalization error in his written review, not in the committed operator.

Let N_+(E), N_-(E) be counts per supplied core. Time is
tau=n_core t/m_X^2, with constant n_core in this diagnostic. For a fixed
rotating-frequency channel nu, the energies are (nu+Omega, nu-Omega, nu).
If W_ij=sum_l (2l+1)|S_ij,l|^2, then a charged incident channel j has rate

    Gamma_ij / n_core = sigma_ij v_j = pi W_ij / (k_j E_j).

W_ij=W_ji follows from the real symmetric radial operator. The two charge
exchange directions shift gas energy and charge by (+/-2Omega,+/-2).
Neutral production shifts them by (+/-Omega,+/-1), after prompt vacuum decay
of the neutral quantum into one positive and one negative daughter.

For a neutral parent of energy E_h, each daughter's lab energy is uniform in

    E_h/2 +/- sqrt(E_h^2-m_h^2) sqrt(1-4/m_h^2)/2.

Both daughters are deposited; the neutral energy is counted once. If support
is truncated, retain the symmetric allowed angular interval and record its
excluded fraction. Never reject a whole decay because a tiny endpoint interval
is outside support, and never project an unsupported energetic particle to the
last energy node.

## Resonant opposite-charge redistribution

For isotropic incident particles, s=2+2(E1 E2-k1 k2 mu). The distinguishable
final-state cross section from the resonant diagram is

    sigma(s)=g_h^4 / [16 pi s ((s-m_h^2)^2 + m_h^2 Gamma_h^2)].

Multiply by the Moller relative velocity sqrt[s(s-4)]/(2 E1 E2), average over
mu with dmu/2, and use the narrow-width delta function. The result is

    <sigma v>_angles = pi g_h^2 / (8 E1 E2 k1 k2)

when |(E1 E2-2)/(k1 k2)|<1, and zero otherwise. The total parent energy is E1+E2;
its decay redistributes energy over the interval above. One X and one Xbar
enter and leave, so this process conserves both species counts and gas energy.
It is not an antiparticle sink or a heat sink. Its effects on core reactions
arise from changing the energy distributions.

The same normalization obeys the independent thermal identity

    gamma_inverse = Gamma_h m_h^2 T K1(m_h/T)/(2 pi^2).

The code checks this against an integrated two-particle phase-space expression
and compares the delta-function approximation with a resolved finite-width
Breit-Wigner integral. It does not add a separate inverse-decay term to the
eliminated-h description. This avoids counting the real intermediate state
twice. Off-shell opposite-charge diagrams and their interference are omitted,
not claimed to have been subtracted or calculated. See Ala-Mattinen et al.,
[Anatomy of real intermediate state-subtraction scheme](https://arxiv.org/abs/2309.16615),
for the general overlap issue; this packet uses only the stated dilute,
prompt-decay reduction.

## Ordinary identical-charge elastic collisions

The amplitude is

    M(XX->XX)=4 lambda_X + g_h^2/(t-m_h^2) + g_h^2/(u-m_h^2).

After azimuthal integration, d sigma/d cos(theta)=|M|^2/(64 pi s), including
the identical-final-particle factor 1/2. The operator integrates the initial
relative direction, CM scattering angle and scattering azimuth. It boosts
both daughters back to the gas frame and deposits both at their actual
energies. Incoming unordered pairs have an additional factor 1/2 for i=j.
The Xbar Xbar channel is charge conjugate. Number and energy conservation
hold per quadrature contribution, and gas-only Maxwell-Boltzmann stationarity
is tested separately. Angular and energy resolution errors are reported.

## What is conserved and what is uncomputed

With the frozen core slope, Q_core=Q0-N_++N_-, and

    I = E_gas - Omega(N_+-N_-)

is conserved. E_core,tangent=E(Q0)+Omega(Q_core-Q0). Kinetic energy is
E_gas-N_+-N_-, and remains explicitly in the ledger. The separate stationary
cold-mass diagnostic is E_stationary(Q_core)+N_++N_-; reaching it would require
actual cooling. Re-solving E(Q) checks the tangent error but does not turn the
frozen-profile trajectory into a self-consistent finite-Q evolution.

The original pulse also retains energy in its disturbed core. The difference
between the archived pulse energy and the stationary postpulse core plus
measured outgoing gas is retained as an unprocessed remainder, with matching
errors disclosed. Its later release is not evolved. A closed-world or fully
cooled final abundance is therefore not supplied by this experiment.

Repeated interpolation of E -> E+2Omega -> E on an arbitrary mesh causes
spurious diffusion even when number and mean energy are conserved. The
authoritative transport grid repeats a graded base interval at exact 2Omega
offsets. Supported forward/reverse exchange then maps between existing nodes.
Uniform aligned and unaligned quadratic grids are preserved as controls,
including their slow convergence and artificial late evolution.
