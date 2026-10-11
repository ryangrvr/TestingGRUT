# CP-09 relational detectors: an explicit reference register and its physical back-action

Date: 2026-10-11. Independent continuation of [frozen CP-09](https://github.com/ryangrvr/TestingGRUT/blob/1dadbd60ae5719508ac40ce7633768549052691a/research/constructive_physics/attempt09/LAW.md), [low-energy response](https://github.com/ryangrvr/TestingGRUT/blob/b87c162d6a838f21e22d40a6f799fe6bfc65c0d5/research/reality_framework/cp09_low_energy_response/REPORT.md) and [disjoint relational roles](https://github.com/ryangrvr/TestingGRUT/blob/8bb18f4c14404d297b17c5748af4669ece7dbe29/research/reality_framework/cp09_relational_probes/REPORT.md). This branch forks from 8bb18f4. **It makes no changes to the frozen Hamiltonian, opens no CP-10 or Card 3, and makes no claim of a physical spacetime or clock.**

## Physical question

The earlier role-defined probes are mathematical block-diagonal observables on geometry × matter. They depend on *all* of C: the averaged source acts on every nonisolated minimal event and the averaged detector on every nonisolated maximal event (isolated events are removed from both role sets). They commute at the same external time and give a nonzero delayed ground-state Kubo response.

Can a separate physical reference register actually implement that operation? What happens if the register retains a record of which roles existed?

## Exact finite reference-register construction

Let f(C)=(S_C,D_C) denote the ordered pair of the nonisolated minimal set S_C and nonisolated maximal set D_C. Let F_N be the number of different values of f across all N-element labelled posets; define orthogonal ancilla states |f>. The ideal reference oracle is the unitary XOR permutation

\[
V:|C\rangle|\psi\rangle|r\rangle\mapsto
|C\rangle|\psi\rangle|r\oplus f(C)\rangle,
\]

with a reference of dimension 2^b, b=ceil(log2 F_N), after assigning distinct b-bit encodings to f. This is a **formally supplied global interaction**, not a physically derived local readout or free observer.

Define B_f as the average of X_i over the nonisolated minimal set S of tag f (or zero when empty). The source B on geometry × matter is

\[
B=\sum_C |C\rangle\langle C|\otimes B_{f(C)}.
\]

Then the following coherent reference-and-uncompute protocol implements the prior source exactly:

\[
\boxed{
V^\dagger\left[\sum_f |f\rangle\langle f|\otimes e^{i\varepsilon B_f}\right]
V\,|\psi\rangle|0\rangle
=
(e^{i\varepsilon B}|\psi\rangle)|0\rangle.
}
\]

The numerical control verifies an exactly cleaned reference and agreement with the ideal impulse (0 discrepancy at machine precision). This is an *implementability-in-principle* dilation, not an energy-conserving local gate network; the roles themselves may require global causal information.

The number of exactly orthogonal role tags has a strong general lower bound. Every proper nonempty subset S of N labelled events can be the lower layer of a **complete bipartite** order with upper layer its complement; such f(C)=(S,S^c) are distinct. Thus

\[
F_N\ge 2^N-2.
\]

Consequently any reference that **perfectly records every possible role identity in mutually orthogonal states** requires b >= ceil(log2(2^N-2)), which is N bits for N>=3. This is **not** a theorem that every implementation of B requires N stored ancilla bits: coherent controls may evaluate f reversibly with other space/time tradeoffs; it is a lower bound specifically for full orthogonal role-recording.

## Exact record-dephasing and excitation-gap bound

Suppose one applies the role-copy oracle to the CP-09 normalized pure ground state |g>, but **retains** a perfectly distinguishable role record instead of immediately uncomputing it. When the reference is unavailable, the reduced state of geometry and matter is

\[
\rho_{\rm rec}=\sum_fP_f|g\rangle\langle g|P_f,\quad
P_f=\sum_{C:f(C)=f}|C\rangle\langle C|\otimes I.
\]

Let p_f=<g|P_f|g>. Its ground-state fidelity after recording is

\[
\boxed{\langle g|\rho_{\rm rec}|g\rangle=\sum_f p_f^2.}
\]

For a nondegenerate ground state with excitation gap \(\Delta=E_1-E_0>0\), the spectral theorem gives the **universal within-this-model bound**

\[
\boxed{
\Delta E_{\rm rec}
:=\operatorname{tr}(H\rho_{\rm rec})-E_0
\ge\Delta\left(1-\sum_fp_f^2\right).
}
\]

The inequality follows directly from \(H\ge E_0|g\rangle\langle g|+E_1(I-|g\rangle\langle g|)\). This is an exact quantum measurement back-action result: a perfectly distinguishable persistent reference can change the CP-09 system energy and subsequent response. It does **not** mean total energy must fail conservation in a closed apparatus+clock+system model; the reference interaction has not been given its own Hamiltonian and energy ledger.

Uncomputing a coherent role tag *before CP-09 evolution* recovers the original pure system state and avoids this dephasing. Holding a role record throughout evolution is a physically different preparation.

## Independent N=2 and N=3 calculation

[check.py](check.py) independently enumerates all labelled posets and builds the unmodified CP-09 Hamiltonian, without importing frozen CP-08/09 code. It constructs the same disjoint role operators as 8bb18f4, checks equal-time commutativity, implements the coherent-copy/controlled-pulse/uncompute protocol, and computes the energy and subsequent response after a persistent record. NumPy/SciPy required.

| Observable | N=2 | N=3 |
|---|---:|---:|
| Number of posets | 3 | 19 |
| Distinct role tags F_N | 3 | 13 |
| Ancilla qubits for full role tag | 2 | 4 |
| Ground energy E0 | 0.0947552837 | 0.3109792266 |
| Gap Δ | 0.3182985418 | 0.1857122051 |
| Ground fidelity Σ p_f² after record | 0.3333427490 | 0.0863387479 |
| Rigorous lower bound on ΔE | 0.2121960309 | 0.1696780458 |
| Actual ΔE after retained record | 1.3112863611 | 3.0125206124 |
| Geometry TV at t=0.4, recorded vs unrecorded | 0.0013901530 | 0.1143498631 |
| Unrecorded source→detector response t=0.4 | 0.0525885135 | 0.1091901443 |
| Response following a retained role record | −0.0020149978 | 0.0467927797 |
| Coherent oracle source-pulse error | 0 | 0 |

The unrecorded responses match the previously published CP-09 disjoint role values without using their implementation. The *recorded response* means: apply the exact role-recording isometry at preparation time, keep the role reference inaccessible (trace it), then apply the same source impulse B and evolve under the same frozen H, measuring A at t=0.4. It is not a conditional response postselected on a particular record value. The geometry TV compares the recorded evolving mixture against the original stationary pure ground-state geometry distribution.

These tests are exact small-Hilbert-space numerical checks of a stipulated ordinary quantum model, not evidence of physical observers, quantum collapse or spacetime emergence.

## Scientific significance

A geometry-dependent intervention can be implemented with a coherently computed and erased reference, so the earlier Kubo result is not algebraically impossible to access. But it requires the very global geometry-role information whose physical accessibility CP-09 did not derive. If an observer keeps a distinguishable record of those roles, they change the subsequent dynamics: the reference is not a passive window onto geometry. This is expected measurement theory applied to the CP-09 model, not a novel law.

This distinguishes **three logically different constructions**: (1) an abstract relational operator; (2) an operationally stipulated coherent reversible oracle implementing that operator; and (3) a persistent physical reference record, which causes dephasing and measurable back-action. None provides a physical relational clock, energy-conserving detector apparatus, Lorentzian light cone or selected causal continuum.

## Next genuine construction obligation

Rather than add another abstract role-dependent projector, specify a **bounded-resource matter apparatus** with an actual interaction Hamiltonian and clock pointer, and test whether its record statistics and conserved total energy approximate the ideal probe without free global access to f(C). The same complete dynamics must preserve the CP-09 interference mechanism and overcome its all-pairs large-N problem, and must select the geometries supporting its inferred signal cones. The new exact gap–record bound can serve as an adversarial benchmark: any claimed nondisturbing *perfect* global role record of a nondegenerate ground state must either show explicitly how the system energy and total apparatus energy balance, or abandon one of its premises.

**Boundary:** no CP-10 is opened, no continuum dimension, spin-2 field, fundamental measurement outcome or novel empirical prediction earned. Earlier frozen packets and Card 3 unchanged.
