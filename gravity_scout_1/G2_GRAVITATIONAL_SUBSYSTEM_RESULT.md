# GRAVITY-SCOUT-1 · G2 — GRAVITATIONAL SUBSYSTEM STRUCTURE (result)

> **Repaired by GRAVITY REPAIR 02 (GR2-01 … GR2-06; `GRAVITY_CORRECTION_LEDGER.md`).** G2 is accepted provisionally
> after that repair. Where wording differs, the ledger takes precedence. Script and log outputs are kept as emitted.

**Question (owner, with GRAVITY REPAIR 01).** What replaces QFT split / type-I subsystem structure when gravitational
gauge constraints are imposed? Does the replacement remove any frozen residual information, or only relocate it into
charges, dressing, boundary observables, or resolution?

**Evidence:**
- the owner-verified sources (Donnelly–Giddings; Raju) and abstract-verified controls (Witten; CLPW);
- `G2_SUBSYSTEM_LEDGER.md`;
- the finite control `g2/g2_subsystem.py` + `.log`. This is type I throughout and an illustration only (GF-11): an
  independent code path, not an independent reviewer.

**G2-6 (no length-scale fishing).** No collar or Planck scale is sought here. G2 is not a second d gate.

## Verdict

| framework | outcome for the type-I intermediate 𝒩 | grade |
|---|---|---|
| flat AQFT (G = 0) | **𝒩 SURVIVES**. It is CONSTRAINED-NONUNIQUE, and canonical given a standard vector (Doplicher–Longo) | standard; DL bibliographically verified only |
| perturbative gravity (Donnelly–Giddings, O(κ)) | **The ordinary AQFT type-I interpolation is no longer the relevant localization structure. A charge-labelled gravitational splitting of Hilbert-space subspaces replaces its operational role** [GR2-01]. No continuum replacement algebra is claimed | **CONSTRAINED / CHARGE-SECTOR GRAVITATIONAL SPLITTING** + **RELOCATION into charges + dressing + perturbative order** |
| fine-grained gravity (Raju, scoped examples) | a Doplicher–Longo-style intermediate type-I factor implementing independent inside / outside specification is **BLOCKED / NOT APPLICABLE** [GR2-02]. No claim that no type-I factor of any kind can occur | **ORDINARY QFT SPLIT INDEPENDENCE FORBIDDEN IN CLASS**, at the verified scope only |
| crossed products (Witten; CLPW) | the type changes (II∞; II₁). No type-I interpolation, independent product preparation or universal gravitational subsystem algebra is inferred [GR2-05] | **NO TYPE-I INTERPOLATION ESTABLISHED** |

- **No unique gravitational analogue of 𝒩.** CONDITIONALLY SELECTED: none at G ≠ 0.
- **Source-backed:** subsystem independence is **OBSERVABLE-CLASS / COARSE-GRAINING / PERTURBATIVE-ORDER dependent**
  (G2-3) [GR2-03].
- **Quantitative physical A_resolution dependence: NOT YET DERIVED.** The ε toy is an illustration only
  (**A_resolution COUPLING CANDIDATE — ILLUSTRATION GRADE**) [GR2-03].
- **TRUE COMPRESSION: 0.** **Empirical payoff: none.**

The frameworks differ, and per G2-5 they are **not** forced into one answer.

## Revised G2 terminal [GR2-06]

- **Flat AQFT:** ordinary split independence / type-I interpolation is available under the split premises.
- **Perturbative gravity:** its operational role is replaced by a charge-labelled gravitational splitting, priced by
  total charges, dressing and perturbative order.
- **Fine-grained boundary-complete gravity** at Raju's verified scope: ordinary split independence fails.
- These statements concern **different observable algebras and levels of approximation**, and are not contradictory.
- **Source-backed dependence:** OBSERVABLE CLASS / COARSE GRAINING / PERTURBATIVE ORDER.
- **Quantitative physical A_resolution dependence: NOT YET DERIVED.** The ε toy is an illustration only.
- No unique gravitational analogue of 𝒩 has been selected.
- **TRUE COMPRESSION = 0.**

## G2-0 / G2-7 Flat and non-gravitational control

**Frozen AQFT** (QFT-SCOUT-1). A split inclusion A_in ⊂ 𝒩 ⊂ A_out with 𝒩 type I buys:
- **normal product extensions** of arbitrary normal marginals;
- **independent state specification** inside and outside the collar;
- **tensor-product implementation**: A_in ∨ A_out′ ≅ A_in ⊗̄ A_out′ (spatially).

**Same matter QFT with G = 0.** The ordinary split inclusion exists (under QP-4), and the type-I interpolation is recovered.

**Toy T0.** A_in = B(ℂ⁴) is a factor (center 1). A_out = 1 ⊗ B(ℂ⁵). Any marginals extend to a product.

**What changes at G ≠ 0:**

| item | G = 0 | G ≠ 0 |
|---|---|---|
| gauge dressing | none needed | gauge-invariant operators carry nonlocal gravitational dressings (QG-10) |
| asymptotic charges | not shared | total Poincaré / ADM charges are measurable outside (QG-8) |
| observable algebra | local net | dressed operators; the boundary Hamiltonian is a boundary observable |
| factorization / independence | type-I split, tensor product | charge-sector-wise (S1); fails in exact boundary algebras (S2, scoped) |

## G2-1 Perturbative gravitational splitting (Donnelly–Giddings; PRIMARY / SOURCE-TEXT VERIFIED by owner)

**Audit** (the owner's checklist):
- **Ordinary local gauge-invariant operators fail:** gravitational gauge invariance conflicts with tensor factorization and
  with a net of commuting local subalgebras already at leading order.
- **Dressings extend nonlocally.**
- **A gravitational splitting exists** for subspaces with fixed matrix elements of the total Poincaré charges.
- **Outside measurements** at the tested order do not resolve internal information beyond those charges.
- **Setup:** an arbitrary extended neighbourhood U_ε; no positive minimum ε at O(κ) (GR1-01).

**Toy T1** (the charge caricature). The outside also holds H_in:
- A_in becomes ⊕_E B(ℋ_E): dim 6, center 3 (the charge sectors). It is no longer a factor. **This direct sum is an
  algebraic caricature only** [GR2-01]. It illustrates shared charge labels, sector-wise independence and the loss of
  arbitrary product specification. It is **not** claimed to be the continuum algebra derived by Donnelly–Giddings.
- H_in is shared by both algebras. A naive ρ ⊗ σ fails to be a product for the dressed pair by exactly Var_ρ(H_in)
  (0.584 here).
- Charge-sharp marginals give products, and each sector is a type-I factor.

**Is it a genuine replacement or a relocation?** Both, and they are kept separate:
- **Replacement (CONSTRAINED).** The independence class genuinely shrinks: interior charge information is no longer
  independent of the exterior. In the toy, charge-indefinite product preparations across the dressed pair are forbidden.
  - This is a candidate H_cross class restriction **analogous** to G1-P1.
  - It is **toy-level**. No exact continuum no-product theorem is claimed. The continuum charges have continuous spectrum,
    and the source works with fixed charge matrix elements, not eigenstates.
- **Relocation.** What remains independent is labelled by supplied data:
  - the charges (QG-8);
  - the dressing prescription (QG-10);
  - the perturbative order (QG-9);
  - U_ε and the state class.

  No supplied information disappears.

**Grade: CONSTRAINED / CHARGE-SECTOR SPLITTING, with RELOCATION of the independence data into charges + dressing.**

## G2-2 Fine-grained split failure (Raju; PRIMARY / SOURCE-TEXT VERIFIED by owner, scoped)

**Established scope.** In the specific gravity settings treated (asymptotically flat / AdS), observables near the boundary
of a Cauchy slice fix the state on the entire slice (holography of information). The ordinary split property — the
ability to specify the state independently on a bounded subregion and its complement — fails there.

**Broader reading.** Not adopted:
- no extension to compact universes or to arbitrary quantum gravity;
- the paper scopes its examples.

**Question:** does fine-grained boundary completeness forbid independent inside / outside specification even with a
collar? **At the verified scope, yes**: the failure concerns the state on the whole slice, so a collar does not restore
independent specification.

**Terminal: ORDINARY QFT SPLIT INDEPENDENCE FORBIDDEN IN CLASS** (the scoped examples only) [GR2-02]. The specific QFT
split implementation (a Doplicher–Longo-style 𝒩 giving independent inside / outside specification) is
BLOCKED / NOT APPLICABLE there. No abstract claim excludes every type-I factor.

**Coarse-graining** [GR2-03]. The source also discusses whether coarse-graining the observable set restores an
approximate split description, and states that the answer / coarse-grained entropy is generally **state-dependent**.

**Toy T2** (algebraic skeleton only; the ingredient list is a summary, not re-read here):
- The outside holds the **exact** vacuum projector P₀ of H_tot (P₀ is a function of H_tot).
- With an outside-cyclic vacuum (λ ≠ 0), alg{1 ⊗ B(K), P₀} = B(H) (dim 400) and the interior commutant is ℂ.
- With a product vacuum (λ = 0), only vacuum-sector data is added (50 / 10).
- **Priced ingredients:** the boundary Hamiltonian (QG-11), vacuum uniqueness and cyclicity (QG-12), exact spectral
  data (QG-11).

## G2-3 Reconciling S1 and S2 (mandatory)

**Test.**
- S1 (Donnelly–Giddings) is **localization modulo total charges, at finite perturbative order**, for exterior measurements
  of the tested class.
- S2 (Raju) is **fine-grained distinguishability by the full boundary algebra**, including exact spectral projections.

**Result: consistent as stated.** The two concern different observable classes and resolutions, so they are not
contradictory.

**Toy illustration (T2 resolution ladder).** At fixed small λ, coarsening the numerical resolution ε steps the
computed outside algebra down through three stages:
- full B(H): no interior subsystem;
- an algebra whose dimensions equal T1's (75 / 6, at λ = 10⁻², ε = 10⁻⁴). This is a **dimension coincidence only, not an
  algebraic identification** [GR2-03];
- vacuum-sector data (50 / 10).

The step thresholds track the Schmidt scale of Ω (~λ). Exploiting fine-grained data requires amplification ~1/s_min.

**Caveat.** ε is a numerical proxy for observational resolution. This is an analogy, not a derivation of any
gravitational resolution limit. The literature debate between these positions is **not adjudicated** here.

**Recorded** [GR2-03]:
- **source-backed:** subsystem independence is **OBSERVABLE-CLASS / COARSE-GRAINING / PERTURBATIVE-ORDER dependent**;
- **illustration only:** the ε hierarchy (**A_resolution COUPLING CANDIDATE — ILLUSTRATION GRADE**). Gravity is **not**
  shown to physically determine A_resolution.

(Superseded wording, kept for the record: "subsystem structure is RESOLUTION / OBSERVABLE-CLASS / PERTURBATIVE-ORDER
dependent (a consistent reading; illustration grade for the resolution mechanism)".)

**Map to the frozen residual:**

| frozen component | gravitational reading |
|---|---|
| **Σ / A_partition** | **CONSTRAINED / OBSERVABLE-CLASS-PRICED** [GR2-04]. The physically useful subsystem decomposition depends on which gravitational observables are included |
| **A_resolution** | **POTENTIALLY LOAD-BEARING; NOT YET DERIVED** [GR2-04]. G2 motivates a coupling to Σ / A_partition but selects no physical resolution scale |
| **A_time** | not yet tested [GR2-04] |
| **A_interface** | **RELOCATION** [GR2-04]: dressing and boundary readout are supplied interface information |
| **gravitational charge sector** | **NEW SUPPLIED LABEL** (QG-8) [GR2-04] |
| **H_cross** | a candidate further class restriction at toy level: product preparations across a dressed pair require charge-sharp interiors |

## G2-4 Crossed-product controls (kept separate; PRIMARY ARXIV ABSTRACT VERIFIED)

| question | **Witten** (arXiv:2112.12828) | **CLPW** (arXiv:2206.10780) |
|---|---|---|
| setting | emergent large-N algebra of single-trace operators outside a black-hole horizon (from Leutheusser–Liu), with 1/N corrections | operators in a de Sitter static patch, gravitationally dressed to an observer's worldline |
| which algebra is type II? | the 1/N-corrected exterior algebra: **type II∞** | the dressed static-patch algebra: **type II₁** |
| what was crossed with what? | per the abstract, the type III₁ algebra **by its modular automorphism group** | the abstract states the dressing → II₁. The crossed-product mechanism is in the body, **not verified here** |
| trace? | yes (semifinite); entropy defined up to a state-independent additive constant | yes (finite); a maximum-entropy state (empty dS); entropy agrees with S_gen up to a constant |
| type-I interpolation? | **not reported** | **not reported** |
| independent inside / outside product preparation? | **not established here.** If the construction supplies an exact factor / commutant pair (e.g. the two exteriors), then G1-P1 forbids normal products for that pair. That pairing is in the body, not verified [GR1-06] | **not addressed** at abstract level |

**Rules observed:**
- No type-I factor is inferred from finite entropy (type II carries finite / renormalized entropy without type I).
- The crossed-product algebra is **not** identified with the Doplicher–Longo 𝒩.

## G2-5 Type-I-factor selection

The selection question was asked only after G2-1 – G2-4. The answer is per framework (see the Verdict).

- **There is no unique gravitational analogue of 𝒩.**
- The closest structure is the S1 charge-labelled gravitational splitting, which replaces 𝒩's *operational role*
  (caricatured as ⊕_charge 𝒩_E in the toy only [GR2-01]). It is **relocated** into
  the charge labels and the dressing choice, **not selected**.

## G2-8 Payoff firewall

- The structural changes to subsystem ontology are entered in the **residual ledger**, not the empirical ledger.
- No candidate passes the seven-criterion bar:
  - nothing is observable at accessible scales;
  - everything is standard to the quoted gravity literature;
  - nothing is residual-input invariant.

**ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS (preserved).**

## Scope and grades

| item | grade |
|---|---|
| Donnelly–Giddings; Raju | PRIMARY / SOURCE-TEXT VERIFIED (owner). The summaries here follow the owner's audit wording, and the ingredient list for S2 is **not re-read** in this environment |
| Witten; CLPW | PRIMARY ARXIV ABSTRACT VERIFIED; body-level claims flagged |
| Doplicher–Longo canonical 𝒩 | bibliographically verified only (carried from QFT-SCOUT-1) |
| Toy T0 – T2 | finite-dimensional illustration; ε is a numerical proxy, not a physical resolution |
| Charge-sharp product restriction | toy-level; no continuum theorem claimed |

All results are **auxiliary to canonical GRUT**: QG-1 – QG-13 and QP-5 are supplied.
