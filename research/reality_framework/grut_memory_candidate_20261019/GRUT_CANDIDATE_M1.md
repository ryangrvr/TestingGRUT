# GRUT-M1-W: an explicit closed matter–memory model

This is a finite Ising repetition register with matter feedback and a
programmed quantum-carrier writer. Its proposed **state-dependent kinetic
saturation** rule is specified below, but the matched uniform weak-field
control retains the record better. Protection cannot be attributed to
facilitation in the demonstrated regime. The scope is a prepared finite
interacting register, not spacetime, Standard Model matter or quantum gravity.

## Theory datum

\[
\mathcal T_{M1W}=(\mathcal H,\mathfrak A,\mathcal C,H,\mathcal S_0,
\mathcal O,\Theta).
\]

For a ring of \(L\ge3\) cells,

\[
\mathcal H=L^2(\mathbb R,dQ)\otimes\mathbb C^2_R\otimes
\bigotimes_{i=0}^{L-1}(\mathbb C^2_{\sigma_i}\otimes\mathbb C^2_{\tau_i}).
\]

The \(\sigma\) spins are matter, the \(\tau\) spins are records, \(R\) is
the carrier's internal source bit, and \(Q,P\) are its position and
momentum, \([Q,P]=i\hbar\). Superscripts \(x,y,z\) denote Pauli matrices
with eigenvalues \(\pm1\), not spin operators divided by two. The
observable algebra consists of bounded operators on this space, together
with the spectral measurement procedures for \(Q,P,H\). There is no gauge
constraint or gravitational constraint. \(\mathcal C\) requires positive,
unit-trace states, finite initial kinetic energy and the stated operator
domains. The external Hamiltonian time \(t\) is fundamental in this model.

The core law is

\[
\begin{aligned}
H_{M1}={}&\sum_i\bigl[J(1+\alpha\tau_i^z)\sigma_i^z\sigma_{i+1}^z
+h_x\sigma_i^x+h_z\sigma_i^z-K\tau_i^z\tau_{i+1}^z
-\Delta F_i\tau_i^x\bigr],\\
F_i={}&\epsilon+(1-\epsilon)\frac{1-\tau_{i-1}^z\tau_{i+1}^z}{2}.
\end{aligned}
\]

All indices are modulo \(L\). Since \(F_i\) acts on neighbors, it
commutes with \(\tau_i^x\); the Hamiltonian is Hermitian. An ordered
neighborhood has flip amplitude \(\epsilon\Delta\); a mismatched
neighborhood has amplitude \(\Delta\). Matter sees the bond operator
\(J_i=J(1+\alpha\tau_i^z)\). The construction supplies this feedback and
the small ratio \(\epsilon\); neither is derived from the words “memory”
or “closure.”

The completed law includes the writer:

\[
H=\frac{P^2}{2M}+H_{M1}
+\Omega\frac{1-R^z}{2}\sum_i f_i(Q)\tau_i^x,
\]

\[
f_i(Q)=\begin{cases}
\cos^2[\pi(Q-i)/(2w)],&|Q-i|<w,\\
0,&|Q-i|\ge w.
\end{cases}
\]

These are fixed spatial couplings. A moving quantum particle meets them in
sequence; there is no prescribed time-dependent pulse. Its kinetic energy
changes when the register exchanges energy with it. Both transmitted and
reflected components belong to the state. The spin and interaction terms
are bounded, so their addition to \(P^2/(2M)\) gives a self-adjoint,
lower-bounded Hamiltonian on the free kinetic domain. Global evolution is
\(\rho(t)=e^{-iHt/\hbar}\rho(0)e^{iHt/\hbar}\).

## Postulates and their price

| Supplied postulate | Mathematical content and purpose | Obstruction addressed; dependent result |
|---|---|---|
| Quantum mechanics | The Hilbert space, unitary law and Born instruments above | Gives exact probability and energy bookkeeping. No derivation of quantum mechanics is claimed. |
| Fixed local substrate | Ring adjacency and two two-level systems per cell | Makes the common interaction neighborhood explicit. Dimension/locality are supplied, unlike the archived emergent-geometry ambition. |
| Kinetic saturation | The displayed \(F_i\), with \(0<\epsilon<1\) | Supplies slow flips on an ordered record and faster motion after defects. Its claimed retention benefit is rejected by the matched uniform-floor control. Its strength is a construction choice. |
| Physical feedback | The same \(\tau_i^z\) changes \(J_i\) | Prevents “memory” from being a disconnected label. The matter readout depends on this term. |
| Modeled preparation resources | A source qubit, massive particle, initial packet and blank register | Addresses the failed bare-writing protocol. An autonomous apparatus is included; spontaneous formation from generic resources is not inferred. |

Closed evolution, explicit energy accounting and local instruments have
independent physical motivation. The particular facilitation projector,
spin content, topology and parameter point are construction choices. Related
constrained quantum models already exist; this proposal has no established
priority claim. See [SOURCES.md](SOURCES.md).

## Preparations and calibration

The numerical benchmark uses \(\hbar=J=1\),

\[
(\alpha,h_x,h_z,K,\Delta,\epsilon)=(0.5,1.05,0.5,1,1,0.05).
\]

The autonomous writer uses \(L=4\), \(M=0.5\), \(p_0=100\),
\(q_0=-4\), position standard deviation \(0.35\), \(w=0.3\), and
\(\Omega=\pi(p_0/M)/(2w)\). The latter is the classical pulse-area
estimate for a flip; the quantum scattering calculation must determine
the actual error. The parameter is fixed before that calculation.

\(\mathcal S_0\) consists of finite-energy density operators, with these
specified experimental subclasses:

1. Prepared records: all memory spins \(|r_z\rangle\), matter spins
   \(|+_x\rangle\), \(r=\pm1\).
2. Core readout: the same records, matter \(|+_z\rangle\) except cell 1
   in \(|+_x\rangle\).
3. Autonomous writing: identical all-\(+z\) blank memories for both
   inputs; source \(R^z=r\); the core-readout matter preparation; Gaussian
   carrier packet with the stated \(q_0,p_0,\sigma_Q\).
4. Formation diagnostic: all memories \(|+_x\rangle\), matter in the
   single-spin ground state of \(h_xX+h_zZ\). This is a different claim
   from preparing a record and is separately tested.
5. Core equilibrium diagnostics: \(\rho_\beta=e^{-\beta H_{M1}}/Z\),
   \(\beta=1\). This is a state of the finite core, not an assertion that
   the writing experiment is in equilibrium. The free particle on
   \(\mathbb R\) has no normalized full-volume Gibbs density matrix.

Calibration would measure isolated spin fields, the two conditional bond
strengths, domain-wall energy and the active/inactive flip-amplitude ratio.
Packet and gate parameters are prepared apparatus inputs. No observed
cosmological or laboratory data has been fitted here.

## Operations and probabilities

For a memory or matter Pauli observable \(A\), use
\(P_\pm=(I\pm A)/2\), probability \(\operatorname{Tr}(P_\pm\rho)\),
and instrument \(\rho\mapsto P_\pm\rho P_\pm\). Equal prior source bits
are distinguished by measuring \(\tau_0^z\) or \(\sigma_1^y\). Their
guessing probability is

\[
P_A=\frac12+\frac{|\langle A\rangle_+-\langle A\rangle_-|}{4}.
\]

This is the success of the specified local measurement, not the optimal
global information-recovery probability. Core record encoding is one
classical bit repeated across \(L\) cells, not \(L\) protected independent
bits. A coherent source input evolves linearly into an entangled source–
register state; the model does not clone an arbitrary unknown quantum state.

The matched control replaces only \(F_i\) by the constant
\(\epsilon\). The low flip amplitude and the \(\alpha\) matter feedback
then remain, while facilitation is absent. It has better sampled retention
at all tested sizes and nearly equal matter readout. The original
\(F_i\equiv1\) control is a twentyfold-stronger transverse field, not
an isolated test of facilitation. The original M1 law is preserved rather
than silently replaced by the better control.

The writer is internal during the calculated evolution. Preparation and
terminal readout remain specified operations; a fully dynamical detector and
its reset apparatus have not been added. Nonselective readout energy changes
are explicitly calculated in [WORKED_REGIME.md](WORKED_REGIME.md).

## Scope of the proposed theory

This law supplies a concrete system in which records affect later
interactions and can persist while other observables spread. It has no
Einstein equation, graviton, horizon, Lorentz symmetry or Standard Model
limit. Changing a spin bond is not being called gravity. No original
Candidate 1 soliton, loop or pole-mass result is imported into this Hilbert
space. [COMPATIBILITY.md](COMPATIBILITY.md) identifies the remaining
obligations before any broader GRUT claim.
