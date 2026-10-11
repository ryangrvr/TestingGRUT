# How large a microscopic history weight must be?

2026-10-11. Additive analytical boundary, extending the verified Gaussian-ensemble closeout. This is *not* a new causal dynamics, quantum theory, CP candidate, Stage-3 action, or proof of spacetime emergence. Earlier frozen packets and PR #12 are unchanged.

## Question and definitions

Let \(\mathcal P_N\) be the labelled \(N\)-element partial orders, equipped **only for this comparison** with the uniform probability \(U_N\). Let \(\mathcal M_{N,d}^{\rm flat}\subseteq\mathcal P_N\) be the orders realizable by \(N\) labelled events in \(d\)-dimensional flat Minkowski spacetime, using the ordinary future-directed causal relation. The dimension \(d\) and Minkowski structure are **inputs in defining this comparison class**, not dynamically generated. This is not a theorem for all curved-manifoldlike causal sets.

## Bound 1: flat-spacetime-realizable orders are combinatorially rare

Take \(N\) events \(x_i=(t_i,\mathbf x_i)\). For each unordered pair \(i<j\), the causal relation is determined by the signs of the degree-1 time polynomial \(t_i-t_j\) and degree-2 interval polynomial
\[
(t_i-t_j)^2-\|\mathbf x_i-\mathbf x_j\|^2.
\]
There are \(m=N(N-1)\) polynomials of degree at most 2 in \(k=dN\) real variables. The standard sign-condition theorem (Warren/Alon, including zero signs) bounds the number of sign patterns by \((8e\,m\,2/k)^k\), for \(m\ge k\). Therefore, for \(N\ge d+1\),
\[
\boxed{|\mathcal M_{N,d}^{\rm flat}|\le
\left[\frac{16e(N-1)}{d}\right]^{dN}
=\exp(O_d(N\log N)).}
\]
This bound is deliberately loose and includes degenerate/lightlike sign configurations. Distinct causal orders can share a sign pattern, so the sign-pattern count is a safe upper bound.

For a fixed labelled bipartition of sizes \(\lfloor N/2\rfloor,\lceil N/2\rceil\), choosing each cross-layer relation independently gives \(2^{\lfloor N^2/4\rfloor}\) *distinct valid* labelled height-at-most-two posets. Thus
\[
|\mathcal P_N|\ge 2^{\lfloor N^2/4\rfloor}.
\]
It follows that
\[
\boxed{
U_N(\mathcal M_{N,d}^{\rm flat})
\le \frac{[16e(N-1)/d]^{dN}}{2^{\lfloor N^2/4\rfloor}}
=\exp[-(\ln2)N^2/4+O_d(N\log N)]\longrightarrow0.
}
\]
This uses only the simple two-layer lower bound, *not* the more detailed Kleitman–Rothschild asymptotic count. For \(d=2\) one can derive the same \(\exp(O(N\log N))\) geometric upper bound more simply from two orderings of null coordinates (\(\le(N!)^2\)); the sign-pattern method handles each fixed \(d\).

## Bound 2: any uniformly well-conditioned finite-species Gaussian is too weak (under a common positive measure)

Suppose every \(C\in\mathcal P_N\) has an *independently supplied*, real positive-definite Gaussian matter Hessian \(K_C\) of the **same** dimension \(m_N\le aN\), with all eigenvalues bounded between fixed \(0<\lambda_-<\lambda_+<\infty\), uniformly over \(C\) and \(N\). Assume the same Lebesgue normalization for every \(C\). The matter integral contributes
\[
w_C=(\det K_C)^{-1/2}.
\]
For any \(C,C'\),
\[
\frac{w_C}{w_{C'}}\le
\left(\frac{\lambda_+}{\lambda_-}\right)^{m_N/2}
\le \exp\left[\frac{aN}{2}\ln(\lambda_+/\lambda_-)\right].
\]
Writing \(P_N(A)=\sum_{C\in A}w_C/\sum_{C\in\mathcal P_N}w_C\), we obtain
\[
\boxed{
P_N(\mathcal M_{N,d}^{\rm flat})
\le \exp\left[-\frac{\ln2}{4}N^2+
O_d(N\log N)+
\frac{aN}{2}\ln(\lambda_+/\lambda_-)\right]\to0.
}
\]
So the earlier special link-Hessian bound of **\(O(1)\)** log-weights extends to a more general **\(O(N)\)** log-weight obstruction for uniformly well-conditioned, fixed-species Gaussian matter in the uniform labelled-poset ensemble. Neither result excludes a *genuinely derived* nonuniform history measure, oscillatory quantum interference, interactions beyond Gaussian theory, additional collective degrees, or a kinetic spectrum that becomes critical as \(N\) grows.

The Gaussian assumption is substantive: a quantum propagator/retarded commutator by itself is not a Euclidean Hessian and does not determine the path-integral measure. A \(\det K\) from a selected action is not a derivation of that action.

## A quantitative target (conditional, not a new law)

At fixed \(d\), with a *positive reweighting of a uniform labelled-poset ensemble*, reversal of the entropic disadvantage of the flat-Minkowski-realizable class requires a class-relative log-weight enhancement on the order of \(N^2\), rather than \(O(1)\) or \(O(N)\). This is a **necessary scale comparison**, not proof of sufficient selection. One possibility is large or critical collective matter effects; another is a genuine Lorentzian phase sum with structured destructive interference, in which positive-weight bounds cannot be applied. Either requires explicit microscopic dynamics, state and boundary conditions, and a common measure over complete history classes.

Do not generalize this counting bound from one fixed Minkowski geometry to **all** possible curved Lorentzian geometries with unrestricted metric degrees of freedom. Do not treat suppression of off-diagonal decoherence-functional terms as elimination of a class's diagonal probability: environmental records can suppress interference without selecting which histories occur.

## Research implication

The next scientific investigation should NOT repeat causal-pair estimates, rank counting, finite positive Gaussians, or retarded determinants. A candidate must give a microscopic evolution/action and a justified state/measure so that the **joint** matter–geometry amplitudes are defined. It must demonstrate class-level, state-independent (or declared-state-dependent) relative suppression while also yielding a physical retarded response, fluctuations and relational records, with the correct continuum limit. A numerical preference between two handpicked histories cannot decide the question.

## Sources and prior art

- Verified parent note: [GAUSSIAN_ENSEMBLE_ENTROPY_BOUNDARY.md](GAUSSIAN_ENSEMBLE_ENTROPY_BOUNDARY.md).
- S. Carlip, *Causal sets and an emerging continuum* (2024), https://doi.org/10.1007/s10714-024-03281-1 . Derives strong *oscillatory gravitational-action* suppression of broad layered families but not continuum emergence.
- Basu, Pollack & Roy, *An asymptotically tight bound on the number of semi-algebraically connected components of realizable sign conditions* (2006), https://arxiv.org/abs/math/0603256 . Sign-pattern counting methodology; the cruder all-sign Warren/Alon upper bound is used above.
- Warren/Alon sign-pattern bound cited in the literature, e.g. *Implicit representation conjecture for semi-algebraic graphs*, https://doi.org/10.1016/j.dam.2018.11.030 . The direct fixed-dimensional Minkowski count above is a deduction from this standard bound, **not a claim of new foundational physics**.
