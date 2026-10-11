# CP09: positive configuration kinetics preserves interference and removes the negative ground-energy obstruction

2026-10-11. One construction, not an unrestricted search. Frozen law local commit c23185e preceded executable tests. Author-level calculations; independent scientific review outstanding. No Stage3 card or earlier frozen law changed.

## Actual mechanism

For each admissible single-relation transition C->C', let L_e psi = psi_C' - U_e psi_C. Postulate H = J sum L_e^dagger L_e + Omega sum(I-X_i), with J=1, Omega=1/4, U_ij=exp(-i pi Z_i X_j/4). This is the conventional positive connection-Laplacian kinetic principle applied to CP08's configuration graph. It preserves the same off-diagonal transport, and adds J times the allowed-transition degree on the diagonal. Configuration degree is NOT causal-link degree. No link cap is imposed.

The physical new assumption is the squared-difference kinetic form. It is not claimed derived from reality, causality or quantum consistency. This is a distinct model, not a correction inserted into frozen CP08. Ordinary quantum mechanics, external time, posets, qubits, all-pairs admissibility and couplings remain supplied.

## Exact energy result

For every state, <H>=J sum ||psi_C'-U_e psi_C||² + Omega sum <I-X_i> >=0. Thus the negative layered variational energy of CP08 cannot occur.

There is also an upper bound independent of numerical diagonalization. A total order has exactly N-1 admissible one-relation neighbors: removing adjacent cover relations. No additions are possible; removing a nonadjacent relation violates transitivity. Put all geometry support on one total order and prepare matter all X+. The matter energy is zero and hopping expectation zero; the diagonal kinetic energy is J(N-1). Hence

0 <= E0(N) <= J(N-1).

This supplies a uniformly bounded ground-energy-per-event range at fixed J. It does NOT prove convergence of that ratio, a thermodynamic state, finite excitation bandwidth, geometric selection or an operational continuum limit. Dense configurations still have large kinetic degrees. In particular the empty-order survival coefficient remains J²N(N-1): a diagonal term cannot remove that first departure probability. Ground-energy stability and dynamical scaling are separate problems.

## Interference survives

The new diagonal terms cannot connect the CP08 origin {0<2} to the full three-event chain in two steps. The same two transport routes therefore yield P_final(t)=3t^4/4+O(t^5). Independent full evolution gives P/t^4=0.75038019,0.75030342,0.75018003 for t=.02,.01,.005. We have removed a negative-energy obstruction without erasing the specific matter-induced interference.

## Executed controls

All labelled posets included at N=2,3,4. Ground energies 0.09475528,0.31097923,0.60375765. Independently summing squared differences agrees with matrix energy within 9.2e-15. Evolution norm checks pass. Turning transport off changes t=1 geometry probabilities by TV 0.00802,0.07970,0.07033 respectively. Thus the interaction remains operationally relevant.

Ground mean heights are 1.66416,2.25483,2.66993. These three sizes do not determine a dimension. At N=4, ground height-two probability is 0.41926 and height-three probability 0.47797. Neither a manifold phase nor suppression of the full layered ensemble follows.

At gamma=0, uniform geometry tensored with all-X+ is an exact zero mode. This is an important adversarial control: positive kinetic energy alone does not select geometry. Nontrivial transport can frustrate global covariant constancy, but that does not prove it favors continuumlike orders. The required class-level test remains open.

## Verdict and next physical work

CP09 is a successful finite positive-energy construction with retained route interference. It removes CP08's precise negative quadratic ground-energy problem. It does not remove the N-dependent short-time empty-state scale, generate a clock, recover Lorentzian response cones, select a manifold or produce Einstein gravity.

A defensible continuation must determine large-system excitation/response scaling and relational detector observables from this same frozen law. Adding preferred dimensions, graph caps or geometry penalties would be another hypothesis. No such changes are hidden here. This is a substantive construction result, not a completed theory of reality.

## Reproduction

python research/constructive_physics/attempt09/check.py (NumPy/SciPy). Imports frozen CP08 Hamiltonian-building definitions as a declared implementation dependency. New positive quadratic-form energy is evaluated independently; finite evolution and sparse eigenvectors give the reported controls. Terminal CP09_POSITIVE_CONFIGURATION_KINETIC_CHECKS_PASS. Results are author computations, not experiments. These are standard connection-Laplacian mathematics; originality of the assembled physical proposal is unestablished.
