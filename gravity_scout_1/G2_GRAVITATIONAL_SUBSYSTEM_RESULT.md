# GRAVITY-SCOUT-1 · G2 — GRAVITATIONAL SUBSYSTEM STRUCTURE (result)

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
| perturbative gravity (Donnelly–Giddings, O(κ)) | **𝒩 REPLACED** by **charge-sector gravitational splitting**. Independence survives modulo total Poincaré charges | **CONSTRAINED / CHARGE-SECTOR SPLITTING**; information **RELOCATED** into charges + dressing + order |
| fine-grained gravity (Raju, scoped examples) | **𝒩 FORBIDDEN IN CLASS** | **QFT TYPE-I SPLIT STRUCTURE FORBIDDEN IN CLASS**, at the verified scope only |
| crossed products (Witten; CLPW) | **no type-I interpolation established**: the type changes, and no 𝒩 is reported | NOT APPLICABLE / NOT ESTABLISHED |

- **No unique gravitational analogue of 𝒩.** CONDITIONALLY SELECTED: none at G ≠ 0.
- **Subsystem structure is RESOLUTION / OBSERVABLE-CLASS / PERTURBATIVE-ORDER dependent.** This is the consistent
  reading of S1 versus S2 (G2-3).
- **TRUE COMPRESSION: 0.** **Empirical payoff: none.**

The frameworks differ, and per G2-5 they are **not** forced into one answer.

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
- A_in becomes ⊕_E B(ℋ_E): dim 6, center 3 (the charge sectors). It is no longer a factor.
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

**Terminal: QFT TYPE-I SPLIT STRUCTURE FORBIDDEN IN CLASS** (the scoped examples only).

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
- an algebra with the T1 charge-sector dimensions (75 / 6, at λ = 10⁻², ε = 10⁻⁴);
- vacuum-sector data (50 / 10).

The step thresholds track the Schmidt scale of Ω (~λ). Exploiting fine-grained data requires amplification ~1/s_min.

**Caveat.** ε is a numerical proxy for observational resolution. This is an analogy, not a derivation of any
gravitational resolution limit. The literature debate between these positions is **not adjudicated** here.

**Recorded:** **subsystem structure is RESOLUTION / OBSERVABLE-CLASS / PERTURBATIVE-ORDER dependent** (a consistent
reading; illustration grade for the resolution mechanism).

**Map to the frozen residual:**

| frozen component | gravitational reading |
|---|---|
| **Σ / A_partition** | The inside / outside partition is no longer fixed by net data alone. It is charge-sector-labelled (S1), or absent in the exact boundary algebra (S2). **A_partition becomes observable-class dependent** |
| **A_resolution** | **Load-bearing for whether a subsystem exists at all**: coarse exterior access gives charge-sector splitting, exact access gives none. A_resolution and Σ / A_partition are **coupled**, no longer independent residual entries |
| **A_interface** | The dressing prescription and boundary readout are interface data (A_readout at the boundary): **RELOCATION** into A_interface |
| **gravitational charge sector** | a **new supplied label** (QG-8) that the subsystem notion now depends on |
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
- The closest structure is the S1 charge-sector decomposition: 𝒩 → ⊕_charge 𝒩_E in the toy. It is **relocated** into
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
