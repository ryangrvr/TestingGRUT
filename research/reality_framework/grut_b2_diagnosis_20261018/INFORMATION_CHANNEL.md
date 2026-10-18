# Interior-fuel information: routing, correlation and recovery

The channel classification is derived for the same color-blind B1 operators.
Numerical values below concern the declared three-event conditional collision
fixture, NOT B2's original eight-event two-fuel production protocol and NOT
the actual growing-geometry preparation. Independent proof review is pending.

## Collective variables replace individual carrier labels

At exactly two color carriers, every primitive color-blind motion induces
either I or SWAP on a declared ordering of color registers. Thus the full
Hamiltonian has H=A tensor I+B tensor S. In growing geometry A,B include
the geometry, receiver, matter, field and fuel configuration; they must not
be traced away before evolution. The representation is bookkeeping, not a
supplied physical metric or a choice of nonlocal energy eigenbasis.

Let P_T=(I+S)/2, P_S=(I-S)/2. A triplet input drives one environment evolution
phi_T; a singlet drives phi_S. Tracing that environment gives exactly

    E(rho)=P_T rho P_T+P_S rho P_S
             +gamma P_T rho P_S+gamma* P_S rho P_T,
    gamma=<phi_S,phi_T>.

For any color-blind accepted event R the UNNORMALIZED instrument instead
has diagonal coefficients ||R phi_T||^2, ||R phi_S||^2 and cross coefficient
<phi_S,R phi_T>. Its success probability can depend on input singlet/triplet
content. This avoids multiplying arrival and readability estimates.

The collective algebra B(C^3_triplet) direct-sum C_singlet is preserved.
A logical qubit (indeed a qutrit) encoded wholly in the triplet can retain
coherence even when the full arbitrary two-qubit space does not. That is
collective information preservation, not preservation of particle identities.
It can include a color stored in fuel and says nothing about whether the
carriers reach the required outgoing-matter receiver.

## What the frozen local fixture actually does

Same B1 registers/coefficient law, three events, middle fuel empty, two
initial co-pairs on outer events. The explicit color-collection convention
orders the registers as plus0, plus1, middle-fuel, plus2, sends the two
occupied color registers to output ports in that order, and then discards
the configuration ancilla. Implementing it requires coherent occupancy-
controlled routing; it is a stronger receiver than fixed outgoing matter
readout. A different configuration-dependent ordering changes the decoder
and can change gamma. The collective triplet statement is unaffected by
such permutations. At theta=6:

* gamma = 0.999995830217 + i 0.000004956901.
* |gamma| = 0.999995830229, a small attenuation rather than a large collapse
  of the internal information channel.
* After a known collective phase correction, reference-entangled-input
  fidelity is approximately 0.999998436336.

The inherited arithmetic accounts give gamma-modulus error <5.16e-8 plus
the stated scalar margin. `CHANNEL_ACCOUNT.json` records a conditional
interval wholly below one. These are evaluated accounts under B2's
binary64 contract, not outward-rounded or independently validated bounds.

The orthogonal color inputs |01> and |10> have output trace distance |gamma|.
If that environment is inaccessible, no deterministic color-only channel
can restore all arbitrary inputs exactly: trace distance cannot be increased
back to one by a subsequent channel. The deficit is small in this fixture.
An allowed exchange word alone did NOT establish a large observed information
failure, and this calculation should not be described that way.

The fuel-occupation distributions differ by about 8.82e-6 between singlet
and triplet input sectors. Fuel therefore carries part of their record;
position/field correlations and coherent fuel/matter amplitudes are included
in the complete Gram account. A coarse fuel-occupation record does not by
itself implement the full environment decoder.

## A specified receiver and explicit recovery operations

For return to the ORIGINAL outer co-pair configuration with fuel empty,
there is one Kraus operator M=a P_T+b P_S. Here

    p_return|T=0.999819322583,
    p_return|S=0.999810500709.

These are different input-dependent probabilities, not one universal
successful-arrival number. A known joint exchange-phase unitary corrects
arg a,arg b. The remaining positive operator is a singlet/triplet filter.
The evaluated maximally-entangled-input conditional fidelity after phase
correction is very close to one; its tiny last-digit deficit is not a
standalone numerical certification of an arbitrary-input decoder.

This high return probability is return to the same initial configuration.
It is not an outgoing transport experiment and cannot rescue B2's 1.93%
outgoing event. It does not prove that the fixture never moved in between.

Possible recovery resources are explicit:

1. **Color access only, environment discarded:** correct the known exchange
   phase; remaining |gamma|<1 attenuation cannot be perfectly undone for
   the full two-qubit input space. The triplet code remains exact.
2. **Full spatial configuration measured, feed-forward allowed:** each result
   r has M_r=a_r P_T+b_r P_S. Apply its joint phase correction and a heralded
   inverse filter. This yields an exact identity SUBCHANNEL with original-
   trial probability sum_r min(|a_r|^2,|b_r|^2). Its evaluated value here is
   0.999991094525. This protocol needs extra measurement, a collective gate
   and a filter-success flag; it is not autonomous and not the B2 protocol.
3. **Coherent access to all degrees of freedom:** the full evolution is an
   isometry on the prepared input. Its known inverse on that image can
   restore the information exactly, for example by reversing the entire
   dynamics. That requires control over fuel, positions, fields and any
   geometry/receiver involved. Global unitarity is not a claim that the
   declared receiver possesses that control.

These distinguish reversible routing, environment-correlated information,
and loss of recoverability in a restricted observable algebra. There is no
fundamental destruction of quantum information asserted. For the actual
growing B1 geometry, the same instrument FORM applies, but its coefficients
and physically affordable decoder remain uncomputed and unproved.
