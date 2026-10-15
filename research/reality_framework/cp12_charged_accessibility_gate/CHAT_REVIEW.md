# Independent review handoff: charged accessibility

Read [REPORT.md](REPORT.md), [DERIVATION.md](DERIVATION.md),
[EXACT_ONSET.json](EXACT_ONSET.json), [NUMERICAL_AUDIT.json](NUMERICAL_AUDIT.json),
and [VERIFICATION.json](VERIFICATION.json).

## Claims to test

1. For Q=2N, the label-invariant target is a connected four-regular graph
   with every edge on a directed electric circuit and n_v>=ceil(N/4).
   The probability upper bound is proved for a larger event, so also bounds it.

2. The full unbounded Hamiltonian has bounded hopping:
   ||K|| <= |t| sqrt(2Q min(N,Q)).
   Check the rotor direct-integral/Fock argument and the self-adjoint-domain step.

3. For e_C>=L and n_v>=k, L1(E)>=2L+k. The neutral initial ray is an H0
   eigenvector and has first-hop norm sqrt(2Q)|t|. Check the finite unitary
   Duhamel remainder, especially the 1/min(N,Q) probability prefactor.
   The resulting fixed-time threshold is constant in N, not a spatial causal cone.

4. N=5,Q=10 has exactly 960 minimum-order charged full-circulation configurations.
   Check the regular-tournament and distinct-hole-endpoint argument. Rerun
   the forward/reverse exact signed amplitude sums, not a path-count surrogate:

   p_G5(T) = [1326415312 / 848779513441154588027109375] T^44 + O(T^46).

5. Check why P_G(H-5)^22 psi0 != 0 proves overlap with the original ray's
   actual cyclic subspace. The explicit overlap bound is an operator-norm
   bound and must not be interpreted as an appreciable transition probability.

6. The N=4 benchmark uses every retained mode and actual initial coefficients.
   Check the full-H residual correction to the absolute target-Gram envelope.
   This finite charged ring is not a second size of the degree-four family.

7. The N=5 budget32 run includes all minimum-order targets, and its outgoing
   residual controls omitted target support. The finite-time probability
   lower bounds all include zero; at T=1 only the upper bound is resolved.

## Reproduction status

All 960 labelled forward amplitudes were also checked directly. The exact core was executed independently in JavaScript; the individual target
also agreed with Python. All eight reverse checks and 680 independent unitary
sanity cases passed. The saved numerical session commands completed before
the container went offline. The packaged audit.py wrapper has not been rerun;
please execute it independently before calling this a complete reproduction.

No formal outward-rounded numerical certificate or external mathematical
certification is claimed.

## Scientific judgment requested

Does this establish strict finite charged-sector access from the original
neutral ray, with a controlled early-time collective upper bound? Are any
infinite-dimensional, permutation-sign, interference, or target-capacity steps
missing?

It does not establish appreciable large-N weight, persistent geometry, a chosen
dimension, a spatial causal cone, or a complete theory of reality. The next gate
requires charged weight and persistence at N>5 with graph choice, followed by
dimension and response in that same unfiltered state.

This file is a review handoff; its publication does not assert that another
chat has received it or independently reviewed it.
