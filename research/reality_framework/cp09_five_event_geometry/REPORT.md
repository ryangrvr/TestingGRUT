# CP09 five-event full-ensemble low-energy test

2026-10-11. Frozen CP09 unchanged. Supporting clock work is not extended here. This is a foundational geometry/spectrum computation, not CP10 or a Stage3 card.

## Independent implementation and exhaustive ensemble

Generate all labelled strict posets by enumerating three choices for every unordered pair: absent, forward, backward, rejecting violations of transitivity. Independently assemble positive connection-Laplacian diagonal and transported hopping matrix from bit-level Pauli operations. No imported CP09 builder, Monte Carlo history sampling or selected geometries. N=3 and N=4 ground energies and geometry statistics reproduce earlier implementations before N=5 is evaluated.

N=5: all 4231 labelled orders, 135392 geometry-matter basis states, 2960192 sparse matrix entries. Lowest computed energies .953282910422, .995164440096, .995164440096, .995164440096, 1.087813046705. Ground energy per event .19065658. Gap .041881529674. Maximum eigenvector residual 4.4e-10. Sparse eigensolvers can miss degenerate multiplicities: no claim that the listed five modes determine all multiplicities or exhaust every low-energy symmetry sector. Residual controls do not constitute a rigorous global spectral certificate.

Together with previous sizes, gaps N=2,3,4,5 are .31829854,.18571221,.09650386,.04188153. This is a decreasing finite sequence, not a proof of gap closure, dispersion, a critical exponent or massless spin2 modes. No extrapolated fit is credited.

## Ground geometry versus counting

| height | N5 ground probability | uniform labelled-poset probability |
|---|---:|---:|
|1|.00024462|.00023635|
|2|.23715805|.19853463|
|3|.55679064|.54597022|
|4|.18528904|.22689671|
|5|.02051764|.02836209|

Mean ground height 2.98867702, mean relation count 5.44232609. Height three is dominant, BUT it is already dominant in the counting comparator. The ground dynamics enhances height-two weight and suppresses the tallest orders relative to the uniform ensemble at this size. Therefore the height-three peak cannot be presented as evidence that quantum dynamics has selected spacetime. Height is not a manifold classifier, and height-two/three categories are not equivalent to asymptotic entropy-dominating families or continuum classes at N=5.

This is a useful adverse finding: the new Hamiltonian's finite ground geometry remains fairly close to raw counting in this diagnostic. Positivity plus non-Abelian transport has not demonstrated class-level geometric selection. That does not prove asymptotic failure; it identifies what the present computation actually supports.

## Next obligation

Physical low-energy modes must be identified through independently defined relational observables, not their gaps alone. A controlled larger-N approximation must preserve the all-order ensemble and distinguish manifoldlike response from generic configuration fluctuations. Neither more accurate apparatus clocks nor one decreasing finite gap can supply that result.

The candidate remains a consistent supplied finite Hamiltonian with bounded ground-energy density, interference and observer controls, not established reality physics. This analysis exposes no new fundamental principle or empirical prediction. No law, couplings or prior packet changed.

Run python research/reality_framework/cp09_five_event_geometry/check.py, NumPy/SciPy. Expected CP09_FIVE_EVENT_GEOMETRY_PASS. Sparse arrays require sufficient memory. Results retain both ground and uniform distributions, energies, dimensions and residuals. Author numerical evidence; independent review outstanding.
