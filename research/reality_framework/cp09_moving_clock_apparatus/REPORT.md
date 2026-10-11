# CP09 autonomous moving-clock pointer response: recovery and physical failure boundary

2026-10-11. One explicit apparatus benchmark, not a new fundamental generating law or CP10. Frozen CP09 unchanged. This bridges intervention and pointer recording in ONE time-independent Hamiltonian without matching the clock spectrum to CP09. It remains an idealized clock, not a physically completed observer.

## Joint Hamiltonian

Supply a clock coordinate q on the line with P=-i d/dq and speed v=1, source center q=0, detector center q=tau=.4. Let f_delta be the unit-area square function of width delta. Supply a pointer r initially |0> and use the previously declared disjoint relational-role source B and detector A. Postulate

H_total = v P_c + H09 - epsilon f_delta(q) B - zeta f_delta(q-tau) A tensor Y_r.

There is no externally scheduled pulse: both windows are fixed clock-coordinate potentials. Epsilon=zeta=1e-4 independently chosen weak-interaction areas. Only delta is varied as a declared approximation-convergence test, not a scan of candidate laws. Geometry/matter observables A,B still require supplied global role access; q is an apparatus coordinate, not derived physical space.

The initial system is CP09's ground state, pointer |0>, and any normalized clock packet with compact support wholly before the source. The linear clock equation translates that packet without dispersion. Once the packet has wholly passed both regions, all its trajectories have identical pointer statistics: initial ground-state evolution contributes a phase and final system-only evolution cannot change the pointer marginal. Thus no specially entangled history state or target-spectrum clock is required for this benchmark. A finite expectation-energy smooth initial packet can be chosen. The code calculates the trajectory-independent marginal using the exact method of characteristics, not a numerical clock-grid simulation.

## One complete intervention/readout experiment

Along a characteristic, traversing the source applies exp[-i(H09-epsilon B/delta)delta]; then free system evolution for tau-delta; then traversal of the pointer region applies exp[-i(H09-zeta A Y_r/delta)delta]. Reverse the source sign in an independently prepared repetition; read pointer X_r at the end. The estimator

-(<X_r>_+ - <X_r>_-)/(4 epsilon zeta)

approaches the ideal source-to-detector Kubo response. It is an actual pointer expectation in this declared apparatus, rather than differentiation between premanufactured stationary clock-history states. Both repetitions use the same clock Hamiltonian and windows; only the specified source sign differs.

| region width | pointer response | absolute error vs ideal .052588513539 |
|---|---:|---:|
| .2 | .052159977589 | .000428535950 |
| .1 | .052482739826 | .000105773713 |
| .05 | .052562155652 | .000026357886 |
| .025 | .052581929064 | .000006584475 |

The errors empirically decrease quadratically over these tested widths. This is not a uniform proof for arbitrary N or a full error bound; weak epsilon,zeta introduce additional finite-difference errors. In the delta->0 limit, finite-dimensional matrix exponentials analytically converge to the supplied ideal source and detector kicks, giving the Kubo derivative as epsilon,zeta->0. Peak interaction norms grow as 1/delta, so ideal recovery has an increasing apparatus-resource cost. No experimental feasibility is claimed.

## Conserved energy and decisive clock qualification

For finite delta the interaction potentials are bounded Hermitian operators. H_total is self-adjoint on the momentum domain and time independent; full joint evolution conserves its total energy expectation for states in its domain. This conservation is analytical; the code computes pointer trajectories, not numerical energy-component exchange of a sampled clock packet.

However vP has no lower spectral bound. Taking packets with arbitrarily negative mean momentum while keeping other terms bounded drives energy to minus infinity. Thus this cannot be credited as a stable physical clock model. Restricting to positive momenta alone is not a harmless repair: multiplication by the interaction windows generally mixes momentum sectors. A bounded-below dispersive clock, its wavepacket errors, reflection, back-action and energy budget must be independently constructed and tested. No such repair is silently included.

External Schrödinger time remains in the autonomous evolution, though the source and detector sequencing is set by clock position. Operational readout after packet passage does not derive time or establish equivalence between different clocks. No spacetime locality or physical implementation of the role oracle is earned.

## Result

A fixed non-spectrum-matched autonomous clock with product preparation and finite-duration source/pointer interactions RECOVERS the ideal response in a controlled idealization. This closes the purely representational gap at the benchmark level. The clock's unbounded energy and supplied global geometry access prevent calling it a healthy physical realization. Both the successful readout and exact failure must remain in the record.

Run python research/reality_framework/cp09_moving_clock_apparatus/check.py (NumPy/SciPy). Expected CP09_MOVING_CLOCK_POINTER_RESPONSE_PASS. Author model calculation, not experimental validation. This is conventional autonomous clock-controlled quantum dynamics; no fundamental novelty claimed.
