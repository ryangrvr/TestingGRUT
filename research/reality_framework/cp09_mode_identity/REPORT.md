# CP09: identify the low-energy mode before calling it geometry

2026-10-11. Frozen CP09 unchanged. A direct foundational test of mode interpretation; no new law or Stage3 card.

At N=3 and N=4 compute 16 low-energy vectors using the independent full-ensemble builder. Probe the entire computed first-excited energy shell with (1) normalized total causal-relation count, (2) collective mean X, (3) labelled local X_0. These probes are declared independently of the eigenvectors. The first two are event-permutation invariant; the third is not. Ground expectation subtraction does not change off-diagonal excited-state weights.

| N | relation-count shell weight | mean-X shell weight | local-X0 shell weight |
|---|---:|---:|---:|
|3|1.33e-32|8.43e-32|.00261678656|
|4|1.76e-31|9.53e-31|.00386211812|

The gaps are .18571220505 and .09650386397. The lowest shell is visible to a labelled matter fluctuation but invisible to these two collective scalar probes. This is consistent with nontrivial permutation-sector excitations. It is NOT a complete representation-theoretic identification, nor a proof that no relational geometry operator can access the shell. Degenerate-shell weights are summed rather than relying on an arbitrary eigenvector basis; sparse eigenvalue multiplicities have not been rigorously certified.

The result eliminates a tempting interpretation: decreasing smallest energy gaps are not, by themselves, evidence for a physical geometry wave. An emergent spacetime claim requires a mode's operationally localized coupling, wavevector or relational separation dependence, dispersion and consistent causal response. None are established by CP09's gap sequence.

The broad theory-construction lesson is actionable: a mechanism must derive the coupling structure of its excitations, not simply create a small eigenvalue. Positive configuration kinetics and transport interference remain genuine finite achievements; they have not yet supplied such propagation. No extra geometric term or fitted projector is inserted to manufacture it.

Run python research/reality_framework/cp09_mode_identity/check.py. NumPy/SciPy. Expected CP09_MODE_IDENTITY_CONTROL_PASS. Results retain energies and absolute transition amplitudes for all computed modes (the field named all_mode_weights contains absolute amplitudes; gap_weights contains squared, shell-summed weights). Author-level control, independent review outstanding. No empirical discovery.
