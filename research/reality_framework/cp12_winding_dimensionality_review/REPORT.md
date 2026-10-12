# CP-12 independent winding response and dimensionality audit

2026-10-12. Independent additive mathematical/physical response to the frozen CP-12 induced Wilson-loop result [8e866d0](https://github.com/ryangrvr/TestingGRUT/commit/8e866d045f5c7419afc203211fc97282b0e8b90b). **No edit to frozen CP-12, no CP-13, no Stage-3 Card 3 and no claim to a theory of gravity.** [check.py](check.py) is independently implemented and reproducible with NumPy; equivalent script checks executed locally before publication.

## 1. The chosen torus families have different thermodynamic limits

The owner's supplied competitors were square-lattice torus graphs, all with Q=2N edges, degree four at each vertex, and **exactly N elementary square cycles**. The "thin" family is C_9 □ C_(N/9) with the narrow circumference fixed at **9**; the "balanced" family is C_(sqrt N) □ C_(sqrt N) for square N. For example:

| N | Thin | Balanced | Winding length thin / balanced |
| ---: | :---: | :---: | :---: |
| 324 | 9 x 36 | 18 x 18 | 9 / 18 |
| 729 | 9 x 81 | 27 x 27 | 9 / 27 |
| 2,916 | 9 x 324 | 54 x 54 | 9 / 54 |

This comparison is **not fixed aspect ratio** as N increases. Its two large-system limits are a *finite-width tube* and a genuinely two-dimensional square-torus sequence, respectively. The short nine-cycle winding channels occur at **density 1/9 per vertex** in the thin sequence; they are not a finite-volume correction that necessarily vanishes in its thermodynamic energy density.

Both candidate graphs are **supplied** competitors, not graph geometries dynamically formed by CP-12.

## 2. Independently exact propagation diagnostics

For the ordinary graph adjacency A on C_m □ C_n, calculate integer-valued closed-walk moments per start vertex. A walk of k<9 steps has no nonzero winding number on either 9x81 or 27x27; both moment sequences agree exactly. At k=9, the thin torus has two oriented nine-step windings from each vertex while the balanced torus has none. Thus

\[
\boxed{\operatorname{Tr}A_{9\times81}^{k}/729=
       \operatorname{Tr}A_{27\times27}^{k}/729\quad(0\le k\le8),}
\]
\[
\boxed{\frac{\operatorname{Tr}A_{9\times81}^{9}
             -\operatorname{Tr}A_{27\times27}^{9}}{729}=2.}
\]

An independent **integer dynamic-program** walk count checks every k=0..13. The difference at k=11 is 242 returns per start vertex; k=10 differs by zero. This fingerprint is **ordinary graph topology**, not a derivation of CP-12's gauge coefficient. It aligns with the owner's observation that physically novel flux states on nine-site windings first occur after 18 virtual half-link hops: topological route length 9, electric flux creation order 18.

The same edge geometry yields a directly propagation-defined test via the heat-return probability of the unweighted graph Laplacian:

\[
P_{m,n}(s)=\frac1{mn}\sum_{p,q}e^{-s[4-2\cos(2\pi p/m)-2\cos(2\pi q/n)]},
\quad
d_s(s)=-2\frac{d\log P}{d\log s}.
\]

| Supplied graph | d_s(20) at N=729 | d_s(20) at N=2,916 |
| :-- | --: | --: |
| Thin 9 x (N/9) | **1.009641** | **1.009641** |
| Balanced sqrtN x sqrtN | **2.003657** | **2.012825** |

At N=11,664 and s=200, thin 9x1296 has d_s=1.000627; balanced 108x108 has d_s=2.001196.

This is not a claim that heat diffusion is the full CP-12 gauge-coupled detector: it is a **static-graph correspondence test**. It shows that the chosen geometries are operationally distinct even though they have equal valence, edge count and square count.

The difference can be seen analytically from cycle heat kernels. Set

\[
p_m(s)=\frac1m\sum_{k=0}^{m-1}e^{-s[2-2\cos(2\pi k/m)]},
\qquad
p_\infty(s)=e^{-2s}I_0(2s).
\]

Then P_(m,n)(s)=p_m(s)p_n(s). In the respective N→∞ limits for every fixed s>0,

\[
P_{\rm thin}(s)\to p_9(s)p_\infty(s),\qquad
P_{\rm balanced}(s)\to p_\infty(s)^2.
\]

By the positive image representation p_9(s)=e^{-2s}\sum_{j\in\mathbb Z}I_{9j}(2s)>p_\infty(s) for s>0. Thus the two candidate families have provably different finite-scale transport spectra in their thermodynamic limits. Their large-s (but below the long-axis saturation scale) spectral dimensions approach one and two. **The dimensionality difference was already encoded by selecting the competitor families, not derived by CP-12 dynamics.**

## 3. Which way does the newly induced flux channel push?

The pinned owner's conditional Feshbach calculation gives a **positive** offdiagonal coupling for each simple length-ell cycle,

\[
\langle W_c0|H_{\rm eff}(E_0)|0\rangle
 =K_\ell(a)t_{\rm hop}^{2\ell}+O(t_{\rm hop}^{2\ell+2}),
\quad K_\ell(a)>0,
\]

where a=Gamma*(M-Q)/M includes CP-12's **actual graph-rewiring diagonal**. An uncharged unit of either circulating flux on this cycle has the diagonal energy offset

\[
\Delta_{\rm flux}=\ell(\epsilon_E-a)>0,
\]

since its 2ell half-link fields carry squared electric energy ell epsilon_E and its ell occupied edges cease to be removable, lowering the positive diagonal rewiring energy by ell a.

In the isolated, declared two-flux-state Schur-complement *channel* (one +flux and one -flux orientation), the formal leading vacuum-energy lowering is

\[
\boxed{\Delta E_{\rm channel}
  =-\frac{2K_\ell(a)^2}{\ell(\epsilon_E-a)}
    t_{\rm hop}^{4\ell}+\text{higher terms}.}
\]

This equation is a **channel-only accounting identity**, not the full vacuum-energy coefficient; other diagonal processes, energy dependence, interacting fluxes and graph fluctuations can contribute at the same or lower orders. The isolated channel sign is negative regardless of K's positive sign. Thus counting this new channel *alone* would **favor the extra short windings on the thin graph**, not the balanced torus. It cannot be invoked as evidence that nature selects balanced/extended 2D geometry.

Using the author's published K_9(a) values unchanged, and N/9 nine-cycle windings:

| N | K9(a) | Isolated thin total coefficient of t_hop^36 | Per-vertex coefficient |
| ---: | ---: | ---: | ---: |
| 324 | 0.0130893332 | -0.0090628124 | -0.0000279716 |
| 729 | 0.0132499634 | **-0.0209904954** | **-0.0000287935** |
| 2,916 | 0.0133471064 | -0.0854315666 | -0.0000292975 |

The per-site channel coefficient tends to a **nonzero finite value** if K_9(a) has its regular finite limit as a→Gamma=0.1 (the explicit positive path formula has finite denominators here). Since the balanced family's first winding length grows as sqrtN, it has no corresponding fixed-ell noncontractible channel. However, there are many other *contractible* loop channels on both graphs, and their complete linked contributions are not computed here.

### Perturbative convergence reality check

The owner's sufficient fixed-graph Neumann condition is 2Q*|t_hop|/delta<1, with delta=2m+epsilon_E/2-a. At N=729, Q=1458, this proves convergence only for

\[
|t_{\rm hop}|<
\frac{2+1/8-181/1820}{2916}
\approx\boxed{6.9463\times10^{-4}}.
\]

At the edge of that extremely conservative sufficient regime, the isolated channel's per-site t^36 contribution has magnitude less than about 10^-118. This is not a physical smallness proof at t_hop=1—the global norm criterion may be far from optimal—but it shows why **presently certified perturbation theory cannot rank the actual frozen-coupling geometries**. A size-uniform linked-cluster theorem or nonperturbative block free energies are still missing.

## 4. What the science gains

- The induced Wilson-loop result supplies an actual *matter-created* shape-sensitive gauge interaction. It is still an established style of strong-coupling physics, not an independently confirmed novel law.
- The tested comparison is **finite-width versus two-dimensional**, not two thermodynamically equivalent 2D aspect ratios. Its nine-cycle winding difference can contribute at nonzero **bulk density** in the finite-width family.
- Static spectral dimension gives an independent propagation observable capable of distinguishing the supplied competing geometry families; it is not a dynamic CP-12 spacetime-selection measurement.
- The positive Feshbach winding amplitude does **not** establish a favorable energy sign for balanced geometry. The isolated channel contribution is negative and could push toward the finite-width competitor.
- At frozen hopping t_hop=1, neither the local degree series nor the Wilson-flux series provides a controlled large-N free-energy ranking.

**Next falsifier with frozen CP-12:** evaluate at least one **common thermodynamic quantity**—a geometry-conditioned block ground-energy *per site*, or a matched canonical finite-temperature free-energy *density*—for families whose both dimensions scale together, PLUS the 9-wide competitor. Track full charged-sector diagonal corrections, flux mixing and rewiring diagonal, and compare the actual lowest-energy eigenstates with a spatially resolved source/detector response. Avoid treating a single vertex coefficient or a Laplacian heat kernel as the requested full result. Gauge sector/cutoff convergence, finite-size corrections and both topology families' scaling must be explicit.

## Reproducibility and provenance

Independent [check.py](check.py) enumerates lattice return walks by integer recursion, evaluates exact finite torus heat spectra, verifies the coefficient table algebraically from the owner's pinned K9(a) inputs, computes the isolated-channel formal energy shift and checks the frozen convergence threshold. A standalone identical implementation passed locally and writes results.json. No frozen CP-12 source file was imported or edited; full CP-12 itself was **not** diagonalized anew. The underlying maths of cycle heat kernels, spectral dimension and the Schur complement is standard. There is no claim of discovery of Einstein gravity or a new fundamental physical prediction.
