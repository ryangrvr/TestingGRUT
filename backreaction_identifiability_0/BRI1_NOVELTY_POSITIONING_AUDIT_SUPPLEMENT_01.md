# BRI1 NOVELTY / POSITIONING AUDIT — SUPPLEMENT 01
## Hostile re-examination under `BRI1_PUBLICATION_VERIFICATION_CHARTER_01.md` §3 (V2 lane)
## Purpose: determine whether the omitted causal location-scale literature downgrades the V2 disposition. This is a hostile supplement, not a rescue exercise.

---

## Comparator A — Causal location-scale noise models (LOAD-BEARING)

**Source:** A. Immer, C. Utans, I. Khemakhem, B. Schölkopf, *On the Identifiability and Estimation of Causal Location-Scale Noise Models*, ICML 2023 (PMLR 202); arXiv:2210.09054. Primary source inspected (arXiv abstract + paper text via search extracts; core definitions and Theorem 1 structure verified).

**Core model:** Y = f(X) + g(X)·N, N ⊥ X, g > 0 (signed-strict-positive in the primary form; the NeurIPS 2023 companion, Strobl & Lasko line, uses the same construction).

**Question 1 — Is the LSNM class mathematically equivalent to the static/single-time form of BRI1's affine competitor?**

**Yes, at the static level.** E2± is F_q = M[q] + G[q]·ξ with one shared ξ-law. Mapping X := q (protocol/history variable), f := M, g := G, N := ξ gives exactly the LSNM functional form. The only structural differences: E2± allows g of either sign (BRI-O2), while primary LSNM assumes g > 0 (the sign-generalization is trivial — the standardized residual absorbs the sign); and E2±'s M, G are functionals of a *history* q_[0,t] (causal in the prefix), while LSNM's f, g are static functions of a covariate. At a single time, with q fixed, the two classes coincide.

**Question 2 — Is standardized-residual/noise-shape invariance explicit or immediate in their framework?**

**Immediate by construction.** In an LSNM, the standardized residual (Y − f(X))/g(X) = N, whose law does not depend on X — this is the *definition* of the model class, not a derived theorem. BRI1's ingredient H ("shape invariant surviving the affine quotient") is the same object: the standardized force. In the static LSNM setting, its invariance is assumed, not proven.

**Question 3 — What does their identifiability theorem identify?**

The causal *direction* (whether Y = f(X)+g(X)N or X = h(Y)+k(Y)N), up to pathological ODE cases. It is a theorem about direction identifiability *given the model holds in one direction* — not a theorem about rejecting the model class itself.

**Question 4 — Do they use multiple environments/interventions?**

Not in the Immer et al. theorem as inspected; the identifiability argument works from the joint law of (X, Y). The LSNM literature's multi-environment extension (e.g., the NeurIPS 2023 companion and related work on nonstationary causal models) treats "environments" as distribution shifts under which mechanisms are assumed invariant — again, invariance of the standardized noise across environments is part of the model's *definition* there, not a derived or tested property.

**Question 5 — Could their theorem subsume BRI1 after replacing X by a protocol/history variable?**

No, for a specific structural reason: the Immer theorem takes as *input* that an LSNM holds in one direction and concludes direction identifiability. BRI1-X1 is a *rejection* result: the physical force family (constructed from Hamiltonian dynamics) provably does **not** admit an LSNM representation with one shared N across the protocol family. Nothing in the LSNM identifiability machinery yields that rejection — it presupposes membership. The dynamical setting adds further structure absent from the static theory: the competitor involves a *shared causal process across intervention histories* (one ξ-indexed path law, one prefix-causal sign functional across all of 𝒳, per BRI-C2±), temporal correlations of the force *process*, and the degeneracy-set condition. The static theorem has no analogue of the sign-functional/prefix-consistency quotient.

**Question 6 — What remains different because BRI1 uses a shared causal process across an intervention family?**

Three things, each checked rather than asserted: (i) the competitor is a *process* law (fdd of a time-indexed force under each protocol), not a static conditional — so the invariance to be violated is across a family of protocols *on one environment*, with the same ξ-object serving every protocol, which has no static counterpart; (ii) the rejection is *constructive* — an explicit Hamiltonian system is proven to leave the class, not merely observed to fail a fit; (iii) the rejection is *quantified* (1/N_B rate) and tied to the reservoir restoration.

---

## Comparator B — Interventional SDE identifiability

**Source:** Zweig et al., *Towards Identifiability of Interventional Stochastic Differential Equations* (arXiv:2505.15987; UAI/PMLR). Abstract inspected.

- **Intervention concept:** stochastic interventions on the SDE system; identifiability of drift/diffusion components from interventional data.
- **Object identified:** model *parameters/components* within a specified SDE family — not non-representability against a competitor class.
- **Latent/shared noise:** the diffusion term is the noise object; no "one shared latent process across intervention environments, whose standardized shape must be invariant" structure appears.
- **Theorem type:** parameter recovery within a model family. **Does not subsume BRI1**, which is a non-representability theorem against a defined quotient class.
- Same verdict as original V2: positioning-relevant conceptual neighbour, not a threat.

## Comparator C — Makri 1999

**Source:** N. Makri, *The Linear Response Approximation and Its Lowest Order Corrections: An Influence Functional Approach*, J. Phys. Chem. B 103, 2823 (1999), DOI 10.1021/jp9847540. Search-level inspection (abstract/summary level; full text not accessible this session).

- Establishes the influence-functional route to linear response for a harmonic bath, with the lowest-order corrections from anharmonicity/nonlinearity captured by higher (multi-time) bath correlation functions.
- **Higher multitime force correlations:** yes, as the structure of corrections (this is the known origin of non-Gaussian corrections — consistent with ingredient B/C being KNOWN).
- **Finite nonlinear corrections / thermodynamic-limit harmonic mapping:** the harmonic mapping is exact in the thermodynamic/linear-response regime; corrections are organized as perturbative terms.
- **Coupling-size scaling / explicit inverse-N rate:** not found at the inspected level.
- **Protocol dependence / interventions:** not part of the framework; single-equilibrium-correlation setting.
- **Non-representability or identifiability statements:** none found.

Verdict: the **reservoir-limit physics** (Gaussian/effective-harmonic recovery; ingredients B, C, K) is prior art and remains so. No overlap with the identifiability theorem (I).

## Comparator D — Makri 2024

**Source:** N. Makri, *Parsing the Influence Functional: Harmonic Bath Mapping and Anharmonic Small Matrix Path Integral*, J. Phys. Chem. Lett. 15, 4616 (2024), DOI 10.1021/acs.jpclett.4c00908. Search-level inspection.

- New relative to 1999: a systematic decomposition ("parsing") of the influence functional into harmonic-mapping pieces plus explicitly non-harmonic remainders, with the small-matrix path-integral methodology for treating anharmonic contributions beyond the harmonic mapping.
- Does it approach BRI1's finite-bath distinction? Not at the inspected level: it is a computational/methodological advance for *including* anharmonic corrections, not an identifiability statement; no protocol family, no shared-competitor class, no escape rate. It reinforces that ingredient B/C/K are KNOWN and that the field *can* compute beyond-harmonic effects — which strengthens, not weakens, the claim that what is missing in the literature is the *competitor-class/identifiability* framing.

## Comparator E — Kiefer et al. (non-Gaussian GLE), arXiv:2505.15665

Original V2 analysis retained (full-text inspected in V2: projection-formalism comparison axis; no protocol family on a fixed environment; no shared-affine competitor class formulated; zero occurrences of "shared"/"exogenous"/"affine"; no 1/N axis; no theorem). **Additional cross-check against the location-scale literature demanded by this supplement:** Paper 1's implicit competitor is a *Markovian Gaussian embedding* — narrower than an LSNM-style affine class (it fixes the noise law to Gaussian, which BRI1's E2± deliberately does not). So Paper 1's "Gaussian approximation fails" finding is strictly weaker than a location-scale rejection even at the static level; it does not close the gap between its framing and the LSNM/causal literature. The "affine competitor" object is genuinely absent from the physical literature as a *defined* class.

## Comparator F — Carcarra & Akay, finite bath

**Source:** A. Carcaterra, A. Akay, *Fluctuation-dissipation and energy properties of a finite bath*, Phys. Rev. E 93, 032142 (2016), DOI 10.1103/PhysRevE.93.032142. Search-level inspection.

- Studies finite-bath FDT and energy properties: transient/non-asymptotic behaviour, energy exchange, deviations from ideal (infinite-bath) FDT.
- **Overlap with BRI1's scaling/result?** At the inspected level: finite-bath deviations from the reservoir ideal are studied, but as FDT/energy phenomenology — not as an escape rate from a defined shared-affine competitor class, no protocol family, no identifiability statement, no standardized-shape invariant. Ingredient J's *role* remains unclaimed there; finite-bath phenomenology (neighbour of A/J) is KNOWN-IN-NEARBY-FORM.

## Comparator G — older nonlinear/non-Gaussian Langevin literature

The strongest pre-2020 primary line inspected at this level: the Zwang/Mori orthogonal-dynamics literature and its explicit-solvent applications (the lineage that Paper 1 continues), plus heteroscedastic/regression statistics (GAMLSS/conditional transformation models — the statistical mainstream where "standardized residuals have one fixed shape, location and scale may vary with X" is the *standard modelling assumption*, decades old).

**Key finding:** in statistics, the statement "location and scale may depend on X while the (standardized) shape is one fixed distribution" is the *definition* of location-scale/conditional-transformation modelling — long predating LSNM causal work. This pushes ingredient G and H further toward KNOWN: the invariance logic is textbook in heterogeneous regression. What the statistics literature does not contain: the dynamical/process competitor, the multi-protocol rejection theorem, the constructive physical counterexample, the 1/N_B escape.

## Comparator H — shared-noise/intervention analogue search

Targeted searches for a theorem with: one common latent noise source; several interventions; invariance of standardized conditional distributions across interventions; rejection of a shared location-scale representation. Result: **no theorem of this shape found.** The nearest objects: (i) LSNM model criticism/misspecification testing (informal, not theorem-level, not multi-environment-with-shared-process); (ii) invariance-based causal discovery (invariant causal prediction and descendants — invariance of *mechanisms* across environments, used to find causal predictors, not to reject a shared-noise representation for a *process*); (iii) nonlinear ICA identifiability (latent variable recovery under distribution shifts). None is a first-class threat.

---

## Re-graded ingredient matrix (F–J)

| Ingredient | Original V2 | Re-graded | Basis |
|---|---|---|---|
| F — shared exogenous competitor | NOT-FOUND | **KNOWN-IN-NEARBY-FORM** | LSNM: one fixed N-law shared across all X-values is the model definition |
| G — affine/location-scale/sign quotient | NOT-FOUND | **KNOWN-IN-NEARBY-FORM** | LSNM is the static affine quotient; sign generalization trivial; textbook location-scale statistics |
| H — shape invariant after standardization | NOT-FOUND | **KNOWN-IN-NEARBY-FORM** | Standardized residual = N by construction in LSNM; conditional-transformation/GAMLSS mainstream |
| I — theorem-level failure of the representation | NOT-FOUND | **NOT-FOUND (narrowly)** | LSNM theorems presuppose membership and identify direction; no rejection theorem for a shared-process family found in either literature |
| J — O(1/N_B) escape rate | NOT-FOUND (as role) | **rate: KNOWN (structurally expected); role as escape rate from the defined quotient: NOT-FOUND (narrowly)** | CLT structurally predicts 1/N non-Gaussian corrections; no paper ties the rate to a defined competitor class |

The owner's stated expectation is confirmed on F, G, H; I and J remain unclaimed but *narrowly* — the gap is now the dynamical/process instantiation, not the invariance logic.

---

## The key equivalence test

**Is BRI1-X1 mathematically just an application of an existing location-scale-noise identifiability theorem to a Duffing bath?**

**No — but the distance is smaller than the original V2 claimed.** The exact theorem-level reason:

1. The LSNM identifiability theorem (Immer et al.) *presupposes* that the data admit an LSNM in one direction and derives direction identifiability up to pathological cases. It contains no mechanism for proving that a given family *fails* to admit an LSNM. BRI1-X1 is precisely such a failure proof. Applying Immer to (X=q, Y=F_q) would require assuming the LSNM holds — the negation of BRI1's conclusion.
2. The static LSNM object is a single conditional distribution; BRI1's competitor is a *shared causal process* across an intervention family (one ξ path law, one prefix-causal sign functional across 𝒳, degeneracy-set condition — the BRI-C2± quotient). The static theory has no prefix-consistency structure, no sign-functional coherence condition, and no process-level shape invariant.
3. The rejection is constructive and quantified: an explicit Hamiltonian system, a proven shape-invariant violation, an explicit 1/N_B rate, and reservoir restoration. None of these is derivable from the static literature.

**However**, the *invariance logic* at the heart of the witness (standardize, compare shape across conditions) is known in both the causal (LSNM) and statistical (location-scale regression) literatures — immediate by construction in the former, textbook in the latter. The novelty therefore does **not** reside in the invariance logic itself, and the original V2's "NOT-FOUND" for G and H was an overclaim.

---

## Final adjudication

The candidate contribution, re-assessed: the **competitor class formalization** (G) and **invariance witness** (H) are KNOWN-IN-NEARBY-FORM (static analogues well established). What remains genuinely unclaimed: the **dynamical/process-level theorem** (I, narrowly) — non-representability of a shared causal process across an interventional protocol family on a fixed physical environment — instantiated **constructively** from Hamiltonian dynamics, with the **rate** (J, narrowly) as the escape measure tied to the defined quotient and the reservoir restoration.

This is a genuine contribution, but its shape is: known invariance logic, applied in a new dynamical setting, with a constructive physical counterexample and quantified escape — an organizing/bridging result between the causal location-scale literature and finite-bath statistical physics, rather than a new mathematical theorem-shape.

**V2 disposition (supplemental): KNOWN-RESULT-NEW-FRAMING**

- It **SUPERSEDES** the disposition in 722ac2's V2 audit (`NOVEL-THEOREM-SHAPE`). The original audit's central factual error was rating G and H as NOT-FOUND; the causal location-scale literature and the location-scale regression mainstream make them KNOWN-IN-NEARBY-FORM.
- The original V2 document is preserved unchanged; this supplement governs.
- The remaining novelty is real but must be claimed as such: *first constructive, quantified, dynamical-instance of shared-location-scale rejection across an intervention family on a fixed environment* — with the static invariance logic credited to the causal/statistical literature.

**Implication for the charter's P3 (nearest prior literature materially weaker/different):** yes — materially different in kind (static parameter-identifiability vs. dynamical non-representability), but the paper must cite and position against LSNM identifiability explicitly.
