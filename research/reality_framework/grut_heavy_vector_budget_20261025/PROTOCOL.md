# Scope recorded before the new mode calculations

Parent 05f4c43; new input review f7d641a2570f148b67e4ad05bd8acbea5f97c960,
spot dictionary in issue comment 6406685915 and review comment 6407150523.

1. Reconstruct the specified conditional prepared state at M_Y/H=30,
   h/H=40, r=.45, c_H=10 and kappa from the given exact formula. Evolve
   the full unchanged scalar equations with the author's own DOP853 engine,
   using only a new runtime kappa argument and heavy-point output windows.
   This tests a supplied preparation, not the phase-holding action.
2. Compute the primordial longitudinal **Proca approximation** on the actual
   scalar Higgs norm. Standard de Sitter Bunch-Davies initial fluctuations,
   followed by the stated instant radiation history. Use the prepared norm
   during inflation; vary it to the aligned norm as a control.
3. The primary positive mode band is k/k_star in [.2,5], with k_star=sqrt(m_A/H),
   e=1e-7. Integrate it with 49 logarithmic nodes and refinement. It is a
   finite-band contribution, not a claim to have counted every vector.
   Keep signed numerical occupations and report any negative contribution.
4. Control the mode start, endpoint, tolerance, profile interpolation and
   constant/smooth-mass comparisons. e=1e-6 is one sensitivity check, not a
   newly adopted allowed cosmology or an optimization scan.
5. Compare the vector deposit plus the Y-only scalar ledger with the
   energy-conserving kinetic-equilibrium proxy in the pinned review. Credit
   the review's implementation explicitly if replayed. The original
   isentropic proxy and the corrected proxy are different premises.
6. Neither this Proca calculation nor a proxy relic replay is the full
   coupled scalar/longitudinal perturbation problem or a kinetic/lattice
   formation calculation. Additional release production must not be added
   twice to perturbative Sigma/Y decays. Initial-state and phase-holder
   assumptions remain visible. Report failure only within the tested domain.

The expected potential problem was identified before these computations:
the review's unfinished gate1 results give a primordial deposit near 48,
larger than its corrected threshold near 8. The result is being checked,
not treated as a completed independent critic verdict.

## Added common-mode transfer control

After the primary Proca band was computed, its finite-band R=40.56 exceeded
even the largest quoted corrected threshold. To test the unresolved
single-Higgs assumption, derive the full six-real-scalar quadratic
Coulomb-gauge Hamiltonian, including the nondynamical A0 constraint. Evolve
the initial common gauge-angle mode through the actual scalar release and
compare late matching at a=10,20,30,40 with the Proca result. This is an
additional physical-approximation check, not a new initial phase mechanism
or an evolution of every initial physical mode. Record the extra scope
before the coupled calculations. The primary verdict still belongs to its
stated quantum-initial-state and energy-conserving-proxy assumptions.
