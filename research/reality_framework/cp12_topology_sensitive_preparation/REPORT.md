# CP-12: topology-sensitive response first appears at cubic order

14 October 2026. Additive research based on the independently reviewed checkpoint [`ce7894165f24b5eeec7ddc47e7bbdf074e9fc3d3`](https://github.com/ryangrvr/TestingGRUT/commit/ce7894165f24b5eeec7ddc47e7bbdf074e9fc3d3). Frozen CP-12, GRUT-RAI, PR #12 and Stage-3 Card 3 are unchanged.

## Outcome

The proposed adversarial test succeeds as a discrimination test. Two phase masks can have identical marked-edge count, identical degree at every labelled vertex, identical initial graph probabilities, purity, coherence, mean energy and energy variance, yet produce different **full-Hamiltonian graph responses at order T³**.

They remain indistinguishable through T², as the parent proved. Their degree-wedge responses also remain equal at T³. The distinction is not a degree-unevenness effect: it detects three-edge correlation topology, specifically the mask's triangles relative to open paths.

This is not spontaneous spatial formation. The mask is supplied, the leading distinction is weak at large N, and its sign reverses on conjugating the supplied phases. No physical principle choosing the mask or phase orientation has been derived.

## Parent audit

We read the full pinned parent derivation and ran its exact committed checker verbatim. It passed all five six-event invariant-observable controls. Its N=Q=6 masks have seven marked connections, identical degree sequence (3,3,2,2,2,2), but zero versus two marked triangles and two versus zero marked squares. The all-observable indistinguishability theorem through T² is sound under the declared neutral preparation and invariant graph-observable assumptions.

The parent's result is a statement about derivatives, not equality at later times. We now identify the first distinguishing derivative explicitly.

## New full-Hamiltonian reduction

Let c_j^H(O) be the coefficient of T^j in the expectation of a graph-diagonal observable under the full Hamiltonian, and c_j^R(O) the corresponding coefficient under neutral rewiring alone. For a neutral initial state,

\[
c_3^H(O)=c_3^R(O)-\frac23t_{\rm hop}^2c_1^R(O).
\]

For the flat-amplitude phase preparations and every event-permutation-invariant graph observable, c₁=0. Therefore all such coefficients through T³ agree with neutral rewiring. This strengthens the parent's order-two reduction, but does not permit substituting neutral evolution at finite times. The earlier four-event study already exhibits matter corrections at T⁴.

The proof retains the full mass, electric, hopping and positive rewiring terms. It uses the exact identities P T_hop² P=2Q t_hop²P and P T_hop R₁ T_hop P=2(Q−1)t_hop²R, where P is the neutral sea sector and R₁ acts on one-hop charged states. DERIVATION.md proves these by counting the removable edges in every allowed transition.

## Independent unregulated full-law certificate

As an auxiliary higher-density control, retain the parent's six-event masks but set **Q=12=2N**. The initial state coherently includes all 455 twelve-edge graphs, with phase z^k, z=−7/8+i√15/8. Both masks have exactly

\[
\langle H\rangle=\frac{156}{25},\qquad
\operatorname{Var}H=\frac{780063}{32500}.
\]

The independent reachable-state Hamiltonian implementation gives:

| Graph observable | Common initial expectation | Common T² coefficient | Square-rich mask T³ coefficient | Triangle-rich mask T³ coefficient |
|---|---:|---:|---:|---:|
| Triangle count | 880/91 | −29/39000 | 11√15/10920000 | 227√15/58500000 |
| Simple square count | 1485/91 | −87/130000 | 1243√15/273000000 | −37√15/9100000 |
| Degree-wedge count | 264/7 | −319/390000 | 41√15/19500000 | 41√15/19500000 |

The linear coefficients vanish. The first difference is cubic. The opposite cubic square coefficients do not mean opposite total initial slopes: both share a negative quadratic coefficient. These are expected **counts**, not probabilities of a square or connectedness.

`full_moments.py` uses exact integer arithmetic in Q(i√15), unbounded electric fields and the full frozen transition rule. One application reaches 11,375 states. The cubic contraction needs only the projection of H²ψ onto the one-action support; omitted states have zero overlap with the relevant bra vectors. This is an exact moment contraction, not a truncated time evolution or a complete many-body diagonalization. All retained states pass Gauss's law.

An additional non-invariant marked-edge-count control has a nonzero linear coefficient. Its full-minus-neutral cubic coefficient is exactly −4√15/375, confirming the correction term above rather than only testing cases where it vanishes.

## Main higher-density family: Q=3N throughout

For N=24,48,96,192, partition the event labels into six-event blocks. Place either parent mask in every block. Between selected pairs of blocks, mark **all** 36 possible cross-connections, using a predeclared round-robin ordering of block pairs. The same cross-connections are used for both masks. Their count makes the marked pool size a differ from M/2 by at most 18.

At every vertex the two masks have the same degree. Their edge-degree-product sums also agree. The triangle-rich mask has N/3 more marked triangles and N/3 fewer marked squares. These differences are extensive, not a single microscopic defect. The controller masks are not proposed spatial geometries; the actual initial graph probabilities are uniform over all Q-edge graphs.

Use the same energy-matching phase for both:

\[
\cos\theta=1-\frac{M(M-1)}{2a(M-a)},\qquad M=\binom N2.
\]

Their exact mean full energy is κQ+ΓQ(M−Q)/M. Their energy variances agree by the neutral Laplacian symmetry and the 2Q hopping norm identity. Both have min-entropy zero and neutral-basis relative-entropy coherence log₂ binom(M,Q). All results remain in this one Q=3N sector; the six-event Q=2N result is only an independent auxiliary certificate.

Let Δc₃ be the triangle-rich-minus-square-rich difference in the T³ coefficient. Lower-order differences vanish exactly.

| N | Δc₃ for triangle count | Δc₃ for square count | Δc₃ for degree wedges |
|---:|---:|---:|---:|
| 24 | +2.97827004×10⁻⁶ | −2.28334036×10⁻⁶ | 0 |
| 48 | +1.40558946×10⁻⁷ | −5.29914743×10⁻⁸ | 0 |
| 96 | +5.22941158×10⁻⁹ | −9.81806522×10⁻¹⁰ | 0 |
| 192 | +1.77213158×10⁻¹⁰ | −1.66191718×10⁻¹¹ | 0 |

These decimals evaluate exact rational coefficients times sinθ. `topology.py` computes colored motif embeddings and exact hypergeometric factorial moments. It independently agrees with closed edge-overlap formulas for the differences. `symbolic_check.py` proves, with formal polynomial arithmetic, that the cubic probability correction has degree three and verifies its leading coefficient for arbitrary M,Q,a,cosθ.

## Closed formula and controlled scaling

Write q₃(k)/sinθ=b₃k³+lower powers, and c=cosθ. With g=Γ/M,

\[
b_3=\frac43g^3(1-c)
\left[8(1-3c)-2M(1-13c)-3M^2(1+3c)+M^3(1+c)\right].
\]

For this degree- and edge-degree-matched family, with Δτ=N/3 extra marked triangles,

\[
\Delta c_3(C_3)=\sin\theta\,b_3\Delta\tau\,
\frac{6(Q)_3(M-Q)_3}{(M)_6},
\]

\[
\Delta c_3(C_4)=-\frac{3(Q-3)}{M-6}\Delta c_3(C_3).
\]

Falling factorials are denoted by (x)_r. These formulas follow from counting overlap of a marked three-edge pattern with a triangle or square. No energy fit or thermodynamic approximation is used.

For Q=3N, a=M/2+O(1), Γ fixed and this block family,

\[
\boxed{\Delta c_3(C_3)\sim36864\sqrt2\,\Gamma^3 N^{-5},\qquad
\Delta c_3(C_4)\sim-663552\sqrt2\,\Gamma^3 N^{-6}.}
\]

`asymptotic_check.py` validates the finite factorial identity independently over many pool sizes and evaluates the closed coefficients through N=6144. Those large-N entries are formula evaluations, not quantum evolution. The powers above are derived asymptotics, not fitted curves.

This suppression rules out interpreting the cubic signal itself as a finite macroscopic selection effect in this family. It does not rule out later collective evolution. No uniform Taylor remainder is proved, and substituting times growing with N into a fixed-N Taylor term would be unjustified.

## What the mechanism reads—and still does not select

Third order accesses three-marked-edge topology. The mask triangles are distinguishable from open paths after degree effects are removed. At this order, the square response is tied to that same three-edge information: it is not independent detection of a marked four-cycle. An independent four-cycle test would need masks with equal three-edge motif statistics and an order-four-or-higher calculation.

The frozen Hamiltonian is real in the declared basis. Conjugating the initial phases preserves graph probabilities, purity, coherence and **the entire distribution of energy measurements**, while reversing every odd-order graph coefficient. The conjugate pair has stronger energy matching than the distinct-topology pair, for which only the first two energy moments are matched. This follows from K H K=H and the corresponding equality of all spectral-projector probabilities; no many-body spectral scan is asserted. The exact checker independently verifies conjugation covariance of the full Hamiltonian action. Therefore the cubic steering direction comes from an externally supplied phase orientation. Recognizing supplied correlation topology is not generating the right topology or a cosmic time arrow.

This pass settles the proposed first-higher-order discrimination question positively, with an exact order and an exact size dependence. It does not establish loop packing, intrinsic three-dimensional volume growth, physical momentum, local propagation, universal gravity or a state-selection law. The calculations do not alter CP-12 to reward triangles or squares.

The next unresolved physical gate is a controller or boundary mechanism that produces relevant correlations without supplying their topology, followed by a same-state high-density finite-time formation and propagation test. The packet provides an adversarial benchmark for such a mechanism. It is not a working theory of reality.

## Reproduction

See README.md and `verify.py`. Verification runs the pinned parent checker and all four new controls afresh. Dependencies and source hashes are recorded in MANIFEST.json. Saved numerical tables alone are not the verification.
