# CP-12: equal-energy preparation phases change geometry response

**Research checkpoint, 14 October 2026.** Parent: [`a2f5b77f471e0c6c562030f559f9f7d69be79649`](https://github.com/ryangrvr/TestingGRUT/commit/a2f5b77f471e0c6c562030f559f9f7d69be79649). Additive study of the frozen CP-12 Hamiltonian. No new candidate action, state-selection law, or spacetime phase is asserted.

## Result

We constructed a preparation-controlled comparison in which the coherent and fully graph-dephased states have **identical initial graph probabilities and identical mean full-Hamiltonian energy**. Two coherent preparations also have identical purity, neutral-basis relative-entropy coherence, and full-Hamiltonian energy variance. Nevertheless, their subsequent graph responses differ, including opposite exact initial changes in square-loop probability.

This advances the proposed experiment beyond comparing a lower-energy coherent state with a higher-energy mixture. It also exposes a sharper preparation obligation: the amount of purity and the first two energy moments do not determine where interference directs the geometry. The phase pattern is additional supplied physical information.

We obtained three levels of evidence, kept distinct:

1. **Exact unregulated certificate:** the full frozen Hamiltonian gives opposite square responses for two four-event preparations with the same energy moments.
2. **Complete finite-time coupled evolution:** all physical states of the four-event, four-connection sector were included at five electric cutoffs, with fermionic signs, Gauss constraints, hopping and off-diagonal rewiring retained.
3. **One higher-density sector, Q=3N:** exact initial graph-response coefficients through T² were evaluated up to N=256. These are full-law derivatives, not finite-time many-body simulations or evidence of three-dimensional space.

## Parent theorem audit

The exact committed `cp12_preparation_minentropy_audit/check.py` was executed verbatim and passed. We read the arbitrary-spectrum extension added at `a2f5b77`. Its argument is valid with an O(N) initial excess-energy bound: conservation of the largest density-matrix eigenvalue alone would not control the unbounded electric tail.

The necessary deficit is measured relative to the neutral graph subspace, not a maximally mixed state of the infinite rotor Hilbert space. The theorem establishes a necessary resource for its microscopic loop-packing target; it does not make low min-entropy sufficient. The resemblance to a cosmological low-entropy boundary condition remains an analogy, not an identification with gravitational entropy. GRUT III's N1 and CP-04's supplied function remain distinct underdetermination problems.

The parent checker exhaustively counts six-event graphs, but its displayed small-N packing upper bounds are vacuous. Its numerical pass is a control on the formulas, not the proof of the asymptotic theorem. The generalized proof is analytical; the inherited checker was not expanded by the parent commit.

## A physically explicit preparation family

Let M=N(N−1)/2, D=binom(M,Q), and |G,0> denote a neutral sea on graph G. Choose a subset S of a=floor(M/2) possible connections and write k_S(G)=|G intersect S|. Prepare

\[
|\psi_{S,\theta}\rangle=D^{-1/2}\sum_{|G|=Q}e^{i\theta k_S(G)}|G,0\rangle.
\]

Every graph initially has probability 1/D. The gate

\[
U_{S,\theta}=\exp\left(i\theta\sum_{e\in S}n_e\right)
\]

preserves connection number and Gauss constraints. Its factors commute. Applied to the coherent uniform neutral sea, it changes only phases. Applied to the incoherent mixture, it changes nothing; it cannot purify that mixture.

With frozen Γ=0.1, choose

\[
\cos\theta=1-\frac{M(M-1)}{2a(M-a)}.
\]

Then exactly

\[
\langle H\rangle_{\psi}=\langle H\rangle_{\rm mix}
=\kappa Q+\frac{\Gamma Q(M-Q)}M.
\]

The phase angle uses connection count, not a fitted spectrum of CP-12. Different subsets S of the same size have the same mean energy and energy variance. The mixture shares the mean, but generally has a different variance. No equality of complete energy distributions is claimed.

This is an externally controlled preparation protocol. The coherent neutral input, the ability to address potential connections, and the chosen mask S are supplied. Starting from the coherent uniform input, the mean system work of this phase operation is ΓQ(M−Q)/M=O(N) at Q=νN. We have not constructed the controller Hamiltonian or proved that this state can be prepared autonomously. A generic mask has O(M) address data; a concise algorithmic mask is not automatically cheap to implement physically. No universal preparation-cost lower bound is inferred from the gate description.

## Exact full-law example

Use N=4, Q=4 and all frozen couplings. Both masks contain three potential connections:

- Concentrated: a star, S={(0,1),(0,2),(0,3)}.
- Spread: a path, S={(0,3),(1,2),(0,2)}.

Here cosθ=−2/3. All 15 graphs have initial probability 1/15. Three are four-cycles, so the initial square probability is 1/5. Both pure states have

\[
H_\infty=0,\qquad \langle H\rangle=\frac{32}{15},\qquad
\operatorname{Var}H=\frac{5402}{675}.
\]

The exact unregulated full-Hamiltonian results are

\[
P_\square^{\rm concentrated}(T)=\frac15-\frac{T^2}{810}
+\frac{2561}{1749600}T^4+O(T^5),
\]

\[
P_\square^{\rm spread}(T)=\frac15+\frac{T^2}{1215}
-\frac{2561}{2624400}T^4+O(T^5).
\]

The first and third derivatives vanish in these certificates. We do not need to assert vanishing of every odd derivative. `exact_check.py` evaluates the full Hamiltonian in exact arithmetic in Q(i sqrt(5)), with no electric cutoff. Its four Hamiltonian applications reach 1,173 physical basis states per preparation; every state passes Gauss's law.

The neutral-only fourth-order coefficients would instead be 1/349920 and −1/524880. Thus full matter/electric dynamics already changes this response at T⁴. A neutral Johnson-graph time curve must not be used as a proxy for later full CP-12 evolution.

All four-event Q=4 graphs are connected. This example tests square-versus-triangle-with-pendant organization, not formation of connectedness or macroscopic geometry.

## Full finite-time test and cutoff budget

We propagated the complete coupled Hamiltonian from both pure preparations and from the uniform neutral mixture. No graph transitions were dropped. The largest matrix has 10,275 states at electric cutoff L=5; the L=1,2,3,4 dimensions are 1,875, 3,975, 6,075 and 8,175.

| Time T | Concentrated phase: square probability | Spread phase: square probability | Dephased mixture |
|---:|---:|---:|---:|
| 0 | 0.200000000 | 0.200000000 | 0.200000000 |
| 0.4 | 0.199834595 | 0.200110362 | 0.200000066 |
| 2 | 0.198224505 | 0.201180840 | 0.200007771 |
| 8 | 0.180525467 | 0.214095971 | 0.200568322 |
| 20 | 0.116041270 | 0.256238869 | 0.201187493 |

These are L=5, dimensionless model results. The two pure preparations have the same variance; the mixture has variance 3601/450, so the pure/mixed comparison matches only the mean. The random-mask control happens to be vertex-isomorphic to the spread mask at N=4 and supplies no third independent example.

At T=20 the largest square-probability changes between successive cutoffs are:

| Cutoffs | Maximum absolute square-probability difference |
|---|---:|
| 1 → 2 | 2.56×10⁻³ |
| 2 → 3 | 6.68×10⁻⁴ |
| 3 → 4 | 4.99×10⁻⁶ |
| 4 → 5 | 1.25×10⁻⁹ |

The L=4→5 change in electric energy is larger, about 9.52×10⁻⁶. Observable convergence must be stated separately. Norm and energy conservation errors are below 10⁻⁹ in every run. This is strong numerical cutoff agreement for the reported square probabilities, not a rigorous infinite-cutoff finite-time certificate. The Taylor formulas above have that exact unregulated status independently.

## Higher-density test: Q=3N, one sector throughout

We used three explicit, unoptimized phase-mask algorithms: a fixed-seed random half of potential connections, a lexicographically concentrated half, and a nearly degree-balanced half obtained from round-robin matchings. None is a supplied spatial lattice. Each starts with the same uniform graph marginal and equal mean energy.

For any graph-only diagonal observable, the first two response derivatives of these neutral preparations receive no hopping contribution. This follows from occupations and Gauss-sector orthogonality, not small hopping. We therefore calculate the full-law initial coefficients using an exact two-color quotient, and independently verify them using rational hypergeometric factorial moments. We do not evolve that quotient and call it full CP-12.

Write <C₄(T)>=<C₄(0)>+c₂T²+O_N(T³), where C₄ counts all simple squares, not disjoint cycle packing.

| N | Initial <C₄> | Random c₂ | Concentrated c₂ | Spread c₂ |
|---:|---:|---:|---:|---:|
| 32 | 143.842040 | +3.83617×10⁻⁵ | +2.13610×10⁻² | −1.71565×10⁻³ |
| 128 | 157.075502 | +6.36917×10⁻⁷ | +2.23219×10⁻³ | −4.28404×10⁻⁵ |
| 256 | 159.503575 | +5.47852×10⁻⁷ | +6.06147×10⁻⁴ | −5.75654×10⁻⁶ |

The concentrated mask's favorable square coefficient accompanies a positive coefficient for the degree-wedge statistic sum_i binom(d_i,2). At N=256 these coefficients are respectively +0.000606147 and +0.004389293. The spread mask gives negative coefficients for both. Thus more loops in this control accompanies increased degree unevenness; it does not establish extended space.

The sign ordering also differs from the auxiliary N=4,Q=4 example. We do not merge their results into a single proposed phase. The Q=3N initial mean square count tends to 162, an O(1) number, not extensive loop packing. No uniform-in-N Taylor remainder or finite-time high-density formation is certified.

There is a useful exact simplification. For balanced color pools, the second-order reweighting is a quadratic function of k_S(G). For any event-permutation-invariant graph observable, its coefficient consequently depends on S only through the number of incident pairs of marked connections (the mask's degree-wedge count). At this order, the apparent loop control is therefore still tied to local degree correlations. It is not independent control of macroscopic topology. See DERIVATION.md.

## What this establishes for state selection

The min-entropy theorem supplies a necessary amount of initial spectral purity for the stated assembly event. This experiment supplies an explicit counterexample to treating that amount, even together with initial graph probabilities and two energy moments, as enough to determine geometric evolution.

The two pure states also have identical relative-entropy coherence log₂D in the neutral graph basis. The difference is the arrangement of phases relative to the frozen Hamiltonian. Changing the state while keeping H fixed is a physical change; changing both by a passive basis transformation would be a different comparison.

If event relabelling is implemented as a symmetry of H, permutation-twirling a preparation leaves every permutation-invariant observable's expectation unchanged. Such twirling is not projection onto the permutation-invariant Hilbert subspace and does not retain pure-state entropy. This packet uses the parent's full labelled ensemble and makes no gauge-quotient claim.

The surviving constructive task is a physical preparation mechanism that fixes correlations, followed by full high-density finite-time tests of loop packing, intrinsic volume growth, and gauge-invariant propagation in the same state. The current result gives an energy-controlled adversarial benchmark for that task. It neither supplies the selection mechanism nor shows it is impossible.

## Provenance and prior art

The quotient technique belongs to established symmetry reduction of quantum walks: Krovi and Brun, *Quantum walks on quotient graphs*, Phys. Rev. A 75, 062332 (2007), [doi:10.1103/PhysRevA.75.062332](https://doi.org/10.1103/PhysRevA.75.062332). Their discrete-time construction is a methodological antecedent, not a derivation of this CP-12 quotient.

Basis-dependent coherence measures are established quantum information theory: Baumgratz, Cramer and Plenio, *Quantifying Coherence*, Phys. Rev. Lett. 113, 140401 (2014), [doi:10.1103/PhysRevLett.113.140401](https://doi.org/10.1103/PhysRevLett.113.140401). No novelty is claimed for phase-sensitive interference or coherence as a resource. The contribution here is the explicit equal-energy construction, exact frozen-CP-12 certificates, and controlled comparisons.

See README.md for reproduction and MANIFEST.json for file hashes and dependencies. Verification executes the scripts afresh; it does not merely inspect saved result tables.
