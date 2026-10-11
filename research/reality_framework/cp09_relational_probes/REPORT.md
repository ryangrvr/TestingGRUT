# CP09 geometry-defined probe response

2026-10-11. Frozen Hamiltonian and couplings unchanged; no new candidate or Stage3 card. This extends b87c162's labelled-probe test with explicitly declared geometric roles. All labelled posets remain in the quantum state.

## Two operational definitions, retained together

For each configuration C let Min(C) and Max(C) be its minimal and maximal events. Define block-diagonal observables B_C=average X over Min(C), A_C=average Z over Max(C). These are invariant under simultaneous event/matter relabelling. No position labels or selected graph are needed.

However isolated events belong to BOTH sets. Thus A and B need not commute at equal time. Executed stationary-ground tests give immediate response i<[A,B]>=.13431237,.12060774,.09118620 for N=2,3,4. This is overlapping support, NOT propagation or acausal transmission. Merely naming probes source and sink does not make them separate physical detectors.

An explicit adversarial disjoint-support diagnostic removes the intersection from BOTH role sets, renormalizing each average by its remaining cardinality. If a set is empty its operator is zero, not a postselected or discarded configuration. This rule was added after the overlap result, is separately disclosed, changes only probe definitions and does not modify the law. With this rule [A,B]=0 exactly, configuration by configuration. Equal-time response vanishes within 1.2e-16.

## Actual delayed response

Use the same frozen CP09 ground state and full Hamiltonian evolution; calculate R(t)=i<[A(t),B]> with these geometry-conditioned observables. Independent impulse exp(i epsilon B) agrees within 1e-8 for both definitions. The impulse itself acts on the full joint state and is therefore a geometry-controlled global intervention, not an already realized local laboratory detector.

| N | disjoint R(.1) | disjoint R(.4) | disjoint R(1) |
|---|---:|---:|---:|
|2|.00412169|.05258851|.02909110|
|3|.00964701|.10919014|.08663352|
|4|.00919626|.09023540|.06801792|

This establishes nontrivial response between disjoint relational ROLE operators in the same fluctuating geometry state. It does not measure travel across a fixed distance: the roles and their membership change with configuration. No spatial scale, momentum, physical clock or signal velocity has been earned. A geometry-controlled operation may require nonlocal knowledge; mathematical relabelling invariance is not proof of physical implementability.

## Consequence

The symmetry obstacle for fixed labels does not prohibit relational responses. This explicit test realizes them without choosing a classical geometry, and shows why their support must be checked. It still falls short of detector localization and continuum propagation. Independently prepared physical reference records would have to realize these or better observables, with clock readings and a controlled approximation. No geometry preference, causal-link cap or target dimension was inserted.

Run python research/reality_framework/cp09_relational_probes/check.py (NumPy/SciPy). Expected CP09_RELATIONAL_ROLE_PROBE_PASS. Both original overlapping probes and disjoint diagnostic results are retained in results.json. The builder definitions are loaded from frozen CP09 without executing its original top-level suite. Author checks only; not experimental confirmation or an independently validated discovery.
