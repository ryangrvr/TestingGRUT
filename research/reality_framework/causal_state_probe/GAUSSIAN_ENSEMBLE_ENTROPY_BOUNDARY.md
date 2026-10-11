# Gaussian causal-history weights: a bounded-selection obstruction

Date: 2026-10-11. Additive research control on \`codex/grut-causal-state-review\`. No CP-07, action principle for GRUT, frozen-packet edit, PR #12 change or claim of derived spacetime.

## Result and premises

Consider the **supplied** finite-dimensional Euclidean matter action with Hessian

\[
K_C=I+zS_C,\qquad S_C=\tfrac12(R_C+R_C^T),
\]

where \(R_C\) is the binary **cover/link matrix** of an \(N\)-element labelled causal order \(C\), with no diagonal entries. Assume the same labelling measure for all orders and positive-definite \(K_C\). Integrating a real scalar with a fixed Lebesgue measure gives the unnormalized Gaussian history factor

\[
w_C=(\det K_C)^{-1/2}.
\]

Neither this Hessian, this measure, the Wick rotation nor the state/boundary prescription is derived by causal order; the following is a **conditional counter-control**, not a proposed quantum-gravitational dynamics.

**Theorem (bounded finite-order reweighting).** If \(z=\alpha/N\) with fixed real \(0<\alpha<2\), then for *every* \(N\)-event finite causal order

\[
\boxed{\quad1\le w_C\le
\exp\left[\frac{\alpha^2}{32(1-\alpha/2)}\right].\quad}
\]

The bound is independent of \(N\) and \(C\).

**Proof.** The undirected Hasse graph has no triangles: a directed acyclic orientation of a triangle would force a transitive edge, which cannot be a cover. By Mantel's theorem its number of links \(L\) satisfies \(L\le\lfloor N^2/4\rfloor\). Since \(S=A_{\rm Hasse}/2\),

\[
\operatorname{tr}S=0,\qquad \operatorname{tr}S^2=L/2\le N^2/8.
\]

Every row has absolute off-diagonal sum at most \((N-1)/2\). Thus each eigenvalue \(y_i\) of \(zS\) obeys \(|y_i|<r=\alpha/2<1\), ensuring \(K_C>0\) throughout the labelled ensemble. Concavity of log gives \(\sum_i\log(1+y_i)\le\sum_i y_i=0\), whence \(\log w_C\ge0\). Moreover, for \(|y|\le r<1\),

\[
0\le-\log(1+y)+y\le\frac{y^2}{2(1-r)}.
\]

Therefore

\[
\log w_C=-\frac12\sum_i[\log(1+y_i)-y_i]
\le\frac{\sum_i y_i^2}{4(1-r)}
=\frac{z^2 \operatorname{tr}S^2}{4(1-r)}
\le\frac{\alpha^2}{32(1-\alpha/2)}.
\]

QED.

## Consequence for a whole history ensemble

Let \(U_N\) be the *uniform measure on labelled \(N\)-element posets*, and \(P_N\) its Gaussian-reweighted measure for this particular \(K_C\). For every set of histories \(A_N\),

\[
P_N(A_N)=\frac{\sum_{C\in A_N}w_C}{\sum_Cw_C}
\le \exp\left[\frac{\alpha^2}{32(1-\alpha/2)}\right]U_N(A_N).
\]

If \(U_N(A_N)\to0\), then \(P_N(A_N)\to0\). This Gaussian cannot turn an asymptotically rare class into a typical class by this reweighting. The uniform *labelled* ensemble is an independent supplied measure. For any different measure, this specific theorem must be interpreted relative to that measure; it is not a universal quantum gravity no-go result.

For a fixed partition of \(N=2m\) labelled events into two layers of size \(m\), each subset of the \(m^2\) potential cross-layer links is a distinct two-layer causal order. There are exactly \(2^{m^2}=2^{N^2/4}\) such labelled orders. The Gaussian partition contribution from this fixed family lies between \(2^{m^2}\) and \(e^{\alpha^2/[32(1-\alpha/2)]}\,2^{m^2}\); the positive bounded Gaussian weights leave this exponent unchanged. The better-known Kleitman–Rothschild generic orders typically have *three* layers, and their combinatorial dominance is established independently in the causal-set literature. Do not conflate that theorem with this fixed-partition two-layer counting identity.

## Exact maximally linked counter-control

For a balanced complete two-layer order with \(N=2m\), the Hasse graph is \(K_{m,m}\). The nonzero eigenvalues of \(S_C\) are \(\pm m/2=\pm N/4\). Hence

\[
\det K_C=1-\frac{z^2N^2}{16},\quad
w_C=\left(1-\frac{\alpha^2}{16}\right)^{-1/2}
\quad (z=\alpha/N).
\]

For any fixed nonzero \(z\), this Hessian eventually ceases to be positive definite: it has a zero eigenvalue at \(N=4/|z|\) and a negative eigenvalue afterward. Thus holding \(z\) fixed cannot give a well-defined **positive Euclidean Gaussian measure over all sizes and all causal orders**. This does not address oscillatory Lorentzian integrals or different kinetic operators.

## Independent algebraic/numerical controls

The accompanying \`gaussian_ensemble_entropy_control.py\` reproduces the same 180-event 4D sprinkling and random two-layer relation-matched control from \`supplied_probe.py\`, and compares a complete balanced two-layer order at \(z=0.01\) (\(\alpha=1.8\)).

| 180-event order | Relations | Links | \(\log w_C\) | \(w_C\) |
|---|---:|---:|---:|---:|
| Supplied 4D sprinkling | 1413 | 897 | 0.011222784465 | 1.011285996159 |
| Random height-two, matched relations | 1413 | 1413 | 0.017678763085 | 1.017835957384 |
| Complete balanced height-two | 8100 | 8100 | 0.113136722162 | 1.119785021912 |

The same fully linked weight 1.119785021912 follows exactly from the formula above for \(\alpha=1.8\). The uniform all-poset log-weight bound is 1.0125; it is deliberately loose. The log-count of the fixed-partition two-layer ensemble at \(N=180\) is \(8100\ln2\approx5614.49216\), dwarfing the \(O(1)\) weight preference. Exhausting all \(2^9=512\) bipartite orders on 3+3 events gives a Gaussian weight sum 539.831946987 and verifies the stated individual bounds. The **same mathematical computations** were independently checked in the analysis environment; the script itself is committed for partner reproducibility.

## What is not earned

- Neither the independent Euclidean scalar action nor the labelled-history measure follows uniquely from causal order.
- The link reward is not a microscopic derivation of Lorentzian propagation, manifoldlike geometry, stress-energy coupling, quantum outcomes or gravitational records.
- A different action, coupling scaling, complex interference or correlated physical measure can evade this toy obstruction. The published suppression of layered orders with a supplied discrete Einstein–Hilbert action is not contradicted.
- No comparison of two individual normalized weights proves the dominance of full history classes.

**Scientific ruling:** The retarded determinant offers no off-diagonal history discrimination. The associated supplied positive symmetric link-Gaussian does discriminate orders, but in the uniform labelled-ensemble/positive-for-all-orders regime its weight is uniformly bounded and cannot overturn a vanishing combinatorial fraction. These are separate, precise reasons *not* to mistake response-kernel mathematics for earned spacetime dynamics.

## Reproduce and read

\`\`\`bash
cd research/reality_framework/causal_state_probe
python gaussian_ensemble_entropy_control.py
\`\`\`

Earlier controls: [RETARDED_DETERMINANT_BOUNDARY.md](RETARDED_DETERMINANT_BOUNDARY.md), [REVIEW.md](REVIEW.md).

Prior art:
- S. Surya, *The causal set approach to quantum gravity*, Living Reviews in Relativity (2019), https://doi.org/10.1007/s41114-019-0023-1 (layered entropic abundance, dynamical challenge).
- S. Carlip, *Causal sets and an emerging continuum*, General Relativity and Gravitation (2024), https://doi.org/10.1007/s10714-024-03281-1 (distinct, supplied gravitational oscillatory action and conditional suppression).
- Mantel's extremal triangle-free graph bound (1907), standard combinatorics. The Hasse-graph application and Gaussian bounds above are direct deductions for the declared toy.
