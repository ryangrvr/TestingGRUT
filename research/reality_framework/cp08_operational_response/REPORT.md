# CP08 operational matter response: causal order is not yet a signal cone

Author diagnostic, 2026-10-11. Frozen CP08 at e3cd89a16478ed3e36debee4373c7a5f33438108 is unchanged. No CP09 or Stage3 card is opened. This is a test of the existing law, not a replacement action.

## Exact result

Use its N=2 Hamiltonian, including all three labelled posets (empty, 0<1, 1<0), J=1, gamma=pi/4 and Omega=1/4. Choose independently declared matter observables B=X_0 (impulse generator) and A=Z_1 (readout). With external CP08 time, define

R_AB(t)=i theta(t) <[A(t),B]>.

Let P_empty project onto the initially empty geometry. Expand A(t)=A+it[H,A]-t²[H,[H,A]]/2+... . Direct Pauli/block multiplication gives

P_empty i[A(t),B] P_empty = -2 t² Y_0 Y_1 + O(t³).

The zeroth-order commutator vanishes; the first-order projected term vanishes because a single hopping step changes geometry. The second-order coefficient is a full matter operator identity, not a fitted scalar coefficient. A product preparation |Y+>_0 tensor |Y->_1 therefore gives R_AB(t)=2t²+O(t³). No initial entanglement or relation is needed. All intermediate geometry configurations and both relation orientations are included.

This response corresponds to the actual derivative of <A> after an impulse exp(+i epsilon B), followed by the same joint evolution. Independent dense matrix exponentiation and finite impulse differentiation agree within 1e-9. R/t² at t=.04,.02,.01 is 1.99726831,1.99931677,1.99982917. The gamma=0 comparator has identically commuting cross observables, because geometry and separate matter spins then factorize.

## Physical meaning and boundary

The initial order contains no relation between these events, yet the supplied global Hamiltonian produces an operational influence. Hence initial-poset incomparability does not by itself imply vanishing cross-response in CP08.

This is NOT a superluminal-signaling result. There is no reconstructed metric, physical spatial separation, relational clock or invariant light cone. The geometry changes during evolution, and paths can create a relation. A dynamically defined operational cone could include this influence. The exact result rules out identifying the initial poset with the complete operational cone without an additional correspondence argument; it does not rule out every evolving-geometry interpretation or continuum limit.

This short-time susceptibility is also not proof of controllable communication at an experimentally calibrated scale. The observables, state preparation, impulse and clock are mathematical inputs. Norms, units and physical implementation remain to be earned.

## Consequence for further construction

A new joint mechanism must define how a matter intervention, clock readings and detector records identify events and causal reachability while geometry fluctuates. Its propagation test must use that definition independently of the measured response; defining a cone afterward to contain every influence is not a prediction. A finite causal-link cap would not resolve this correspondence and would introduce its own Lorentzian difficulties.

The result does not repair CP08's separate layered-sector variational bound. A continuum proposal must address both operational propagation and the full geometry ensemble. A conditional local wave operator on a supplied sprinkling is a valuable comparator but cannot count as dynamically generated geometry.

## Reproduction

Run python research/reality_framework/cp08_operational_response/check.py from the repository. NumPy and SciPy required. It imports frozen attempt08/check.py without executing that file's test suite, regenerates results.json, checks the projected Pauli identity, compares Kubo and impulse responses, and checks gamma=0 factorization. Expected terminal CP08_OPERATIONAL_CROSS_RESPONSE_PASS.

No fundamental discovery, relativistic causality or new model is certified. This is a precise physical correspondence gap exposed by the existing model.
