# D1-01 — controlled effective histories of the existing interacting parent

This is a connected effective-description test of M1-01's E=12 finite parent. Its continuum matching failed and is not reopened here. Freeze this card and code before evaluating the new time protocols. The existing M1-01 outcomes are known; the new timings are a numerical check within the same adopted model, not untouched observational evidence.

## Parent and retained description

Use exactly the E=12, 128-node quadrature Hamiltonian, interacting ground state, compressed Weyl operator and absorbing overflow instrument of M1-01, m/g=0.5 and L=5. Reconstruct them with the frozen `bosonized.py`. First reproduce every archived parent history at (0), (0,1), (0,2), (0,1,2) to TV <1e-10.

Diagonalize the parent and retain its lowest K=4,8,16,32 energy eigenstates. This basis includes virtual free-boson components in the interacting eigenvectors; it does not change the parent. P_K commutes with the parent's Hamiltonian, so free evolution inside that subspace is exact. The ground-state preparation lies in every retained subspace. Compress the parent's ordinary Kraus operators again and send the additional discarded probability to the absorbing overflow register. Retain all coherences within P_K.

Evaluate the new protocols (0,0.7,1.9) and (0,0.6,1.5,2.4), with the same smearing, strength 2 and phase 0. The frozen primary reduction is K=32. It must satisfy TV <=1e-5 for both protocols, and the explicit bound below must be <=1e-3 and dominate the measured TV up to 1e-10 rounding tolerance. Report every K, not only the selected primary. A miss remains a miss; do not choose K afterward.

## Conditional finite-parent bound

Let M_s be an ordinary parent Kraus operation, including its preceding unitary interval. It is a contraction. For a branch with embedded reduced vector chi and known error bound d, the next ordinary branch satisfies

$$
\|M_s\psi-P_KM_s\chi\|
\le d+q_s,\qquad q_s=\|(1-P_K)M_s\chi\|.
$$

Thus recursively set d_child=d_parent+q_s, starting from d=0. Because P_K commutes with H and chi is retained, q_s=||(1-P_K)W chi||/2 for either sign. At the end,

$$
|p_s-\widehat p_s|\le 2\|\chi_s\|d_s+d_s^2.
$$

Both distributions have one aggregate absorbing-overflow probability equal to one minus the ordinary sum. Hence

$$
\mathrm{TV}(p,\widehat p)\le
\min\left(1,\sum_{s\in\{+,-\}^n}[2\|\chi_s\|d_s+d_s^2]\right).
$$

Compute q from the norm of the actual omitted vector, not a subtraction of nearly equal squared norms. This is a deterministic bound for the specified finite parent and reached branches, with numerical roundoff reported. It supplies no bound between E=12 and the continuum and is not uniform over other preparations or interventions. Rare postselection would need a separate conditional-error denominator.

## Comparator with the same calibration

Use all 1280 parent energy populations but dephase each retained ordinary branch in the interacting energy basis after every detector. Its transition matrices are |(I+s W)/2| squared entrywise. It has exactly the same initial preparation and first-outcome probabilities as the finite parent, and its dephasing leaves energy populations and mean energy unchanged. Subsequent time intervals cannot create coherence in that comparator. It is a deliberately changed instrument model, not an unphysical reset of energy or a claim that physical decoherence occurs spontaneously.

Require the single-time calibration to agree to TV <1e-10. The predeclared hypothesis is that at least one new history differs by TV >=0.001. If it does not, report no resolved need for coherence by this test. If it does, it rejects energy populations alone for these target histories, not every alternative reduced state or every physical decoherence model. Do not describe this standard quantum interference mechanism as a new GRUT law.

## Resources and outputs

One reconstruction, one BLAS thread, 300 seconds wall cap, no enlarged parent or refitted parameter. Preserve full parent, reduced and dephased probabilities, extra projection losses, the branchwise bounds, source hashes and numerical checks. Keep the baseline parent overflow separate from the additional reduction diagnostic. The publication pause and M1-01 verdict are unchanged.
