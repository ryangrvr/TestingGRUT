# Independent two-star embedding audit and graph geometry of the low-energy trial

**2026-10-13. Additive follow-up to Q=N ground-sector theorem candidate.** Frozen CP-12 and canonical GRUT-RAI unchanged. Results are *independent physical-space numerical controls* and a *separate exact hypergeometric theorem* about the explicitly constructed coherent trial state, **not** a measurement of the unknown exact ground-state graph distribution.

## A. Physical Hilbert-space embedding of the 32-state two-star operator

Choose any adjacent degree-two vertices i,j with external edges L,R and mutual edge M. The pair-star operator h_i+h_j acts only on the five fermionic modes (a_i,a_j,b_L,b_M,b_R) and on electric incidence legs of L,M,R. It leaves fixed every other mode and the **outer half-leg electric field** z on L.

For each fixed configuration of spectators, the canonical anticommutation relations admit a unitary fermionic mode permutation (with its standard occupancy-dependent signs), moving the five active modes to the first five tensor factors. This turns the physical local hopping terms into the usual five-mode fermion hopping algebra. Different signs arising from external edge orientation can be removed by rephasing the external edge fermion modes; the two incidences of the shared edge necessarily have opposite signs, as in the certificate.

Enforce only the two central vertex Gauss constraints and the three edge-center Gauss constraints. They determine the remaining five electric half-legs from the active occupation bits and z, exactly via the map in \`cluster_bound.py\`; the *outside vertex Gauss constraints* and global fermion-number constraint select subspaces and may be dropped safely for a **lower** operator bound. Consequently, in every physical spectator sector the exact pair-star is the compression of the certificate's 32-state operator at integer z. When the two external edges meet the same third vertex (triangle), the outside vertex constraint is still only an additional restriction; there is no assumption that the two outer vertices are distinct.

This provides the missing algebraic embedding justification, subject to independent inspection of the mode-ordering unitary and frozen edge-orientation convention.

**Independent complete Gauss-space control:** An original physical-ring solver using the frozen mass/flux/cost/neutral-diagonal and global Fock signs constructs h_1+h_2 *directly on every physical basis state* for rings N=4,5 at electric cutoffs L=2,3,4. It checks that hopping never changes the exterior spectator bits or electric-flux labels, and computes the lowest eigenvalue:

| N | Cutoff L | dim | inf spec(h_1+h_2) at a=19/200 |
|---:|---:|---:|---:|
| 4 | 2 | 257 | -0.2155081441319937 |
| 4 | 3 | 397 | -0.2155081441319938 |
| 4 | 4 | 537 | -0.2155081441319936 |
| 5 | 2 | 874 | -0.2155081441319939 |
| 5 | 3 | 1378 | -0.2155081441319944 |

All exceed -27/125=-0.216. The same independently written full fixed-ring Hamiltonian reproduces Work's C6 energy at a_6=.06: L3 -0.6352035378471174, L4 -0.6352035378471245, and C9 at a_9=.075: L3 -0.8414821890327315. Thus the previously reported -0.5696 disagreement is reproducible as a *mismatch in the alternative solver*, not as a problem with the established target; its precise source in Claude's solver remains unverified.

These controls alone do **not** replace the all-N rational operator proof.

## B. NEW exact conclusion about the coherent trial: extensive isolated vertices

The Q=N coherent incidence-matching variational state used to bound E_0 is **not the exact ground state**. At excitation-number level k it is an equal-modulus superposition over all admissible matchings of k selected vertices to k *distinct* edges, and over every Q-edge graph containing those locked edges. Different matchings are orthogonal in their matter/electric occupations, so every *graph-diagonal* observable is equivalent to first sampling a matching and then sampling the remaining Q-k connections **uniformly without replacement** from the M-k allowed positions.

Fix any matching, and let U be any set of N-2k vertices not incident on its k locked edges (there must be at least that many). For v in U, the probability that it is isolated is **exactly**

\[
p_N=\frac{\binom{M-k-(N-1)}{Q-k}}{\binom{M-k}{Q-k}}.
\]

For distinct v,w in U, the joint isolation probability is

\[
p_{2,N}=\frac{\binom{M-k-(2N-3)}{Q-k}}{\binom{M-k}{Q-k}}.
\]

Let I_U count isolated vertices in U. Conditional on the matching,

\[
\mathbb E[I_U]=(N-2k)p_N,\qquad
\operatorname{Var}(I_U)=(N-2k)p_N(1-p_N)
+(N-2k)(N-2k-1)(p_{2,N}-p_N^2).
\]

At Q=nu N, k=xN+o(N) with 0<x<1/2 and nu>=x, we have p_N→exp[-2(nu-x)], and Var(I_U)=O(N). Hence by Chebyshev the total graph-isolate count I satisfies, *uniformly over matchings*,

\[
\boxed{\Pr\!\left[I/N\ge (1-2x)e^{-2(\nu-x)}-\varepsilon\right]\to1
\quad(\varepsilon>0).}
\]

For the Q=N trial tent wavepacket centered at k=N/7 and width floor(sqrt N), \`x=1/7\`, \`nu=1\`, so

\[
\boxed{(1-2x)e^{-2(1-x)}
=\frac57 e^{-12/7}=0.12863736581996596\ldots.}
\]

Thus the specific coherent variational energy witness has **at least approximately 12.86% isolated vertices asymptotically with probability tending to one** (minus any fixed small tolerance). In particular its connected-graph probability tends to zero. This fact is independent of the energy bound and holds despite the trial's vanishing positive rewiring-Laplacian energy.

**CRITICAL LIMIT:** Rayleigh–Ritz only proves that the *true* ground energy is <= the trial energy. It does NOT prove the true ground state resembles this trial, lies in the rewiring kernel, has extensive isolated vertices, or is necessarily fragmented. The Q=N degree-two exclusion gives no direct lower bound on the ground-state connectedness. Statements that the coherent matching state "wins" or "is the ground-state geometry" are unsupported.

## Reproduction and implications

- \`physical_two_star.py\` independently constructs the exact finite physical Gauss-law basis and two-star matrix, checking full global Fock signs and exterior invariant sectors. Requires NumPy/SciPy.
- \`trial_isolates.py\` uses stdlib gamma-function binomial-ratio formulas to check exact conditional first and second moments, the O(N) variance, and finite-size concentrations.
- Neither script changes or approximates the frozen Hamiltonian; the local physical-space check has a finite *numerical* electric cutoff but the cited exact rational 32x32 proof covers the unregulated z tail.

**Revised scientific ruling:** The Q=N near-degree-two *ground-state exclusion* has independent physical embedding controls; the low-energy coherent trial itself is rigorously nonspatial/disconnected. The full ground state's actual connectivity, graph distribution and operational spacetime properties **remain unknown**. A full-CAR algebraic audit of the two-star is advisable before describing its all-N proof as independently finalized. Do not infer a dynamical selection mechanism from the trial alone.
