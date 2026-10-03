# BRI0 — BACK-REACTION IDENTIFIABILITY 0: CHARTER 0 (charter only; NOT RUN)

**Status:** CHARTER 0. **No bath simulation, no numerical result, no candidate Hamiltonian, no witness protocols, no
bath-size grid.**

> **BRI0 SCOPE REPAIR 01 (§R at the end) supersedes the text below wherever they differ.**
> - **Owner ruling:** Charter 0 CONDITIONALLY ACCEPTED.
> - The primary competitor E₂ (G > 0) is replaced by **E₂±** (G ∈ ℝ \ {0}).
> - BRI-C2 is replaced by **BRI-C2±**.
> - The §8 witnesses are made reflection-safe.
> - BRI-UPPER is narrowed to its proved scope. Analytic statements are marked **PROVED HERE (INTERNALLY PROVED / NOT EXTERNALLY REVIEWED)**,
**KNOWN** (literature or standard technique, metadata only) or **CONJECTURE**. No PR, no merge.

**Branch:** `grut-backreaction-identifiability-0`, created from the frozen scientific parent
**`scout-0 @ ab2da47407bb670d94e2a52c87599fa13fd8ab99`**. That parent is not modified.

**Provenance:** BRI0 is scientifically downstream of SCOUT-0 (P-17, P-18 and the saturation map). It does **not**
descend from the DA0 / C3 architecture programme, and nothing here modifies:
- `grut-directed-autonomy-0`;
- `grut-reflexive-autonomy-0-frozen`;
- `scout-0`.

> **BRI0 IS A NEW-PREMISE CAMPAIGN.** It does not reopen, weaken or re-grade the SCOUT-0 saturation verdict.
> `SATURATION_CERTIFICATION_02.md` row 13 records nonlinear back-reaction parents as needing a new premise: "future
> campaigns, not unfinished SCOUT-0 business". BRI0 is such a campaign.

## 1. Inherited results (unaltered)

**P-17** (`playground/SCOUT_0/probes/P17_RESULT.md`, theorem). Within **class 𝓗** (system with any smooth potential;
harmonic bath; coupling linear in the bath coordinates **and** in the system coordinate, with counterterm; preparation
hypothesis (P): the free-force law does not depend on the system's initial state), the exact law of the reduced system
path equals that of the GLE with the same memory kernel, driven by an exogenous force with the same law.
**Reduced data identify the effective forcing law, not its microscopic ontology, within 𝓗.**

**P-18** (`P18_RESULT.md`). The same structural non-identifiability holds for a one-way deterministic chaotic driver (a
skew product with no back-reaction): exact at finite ε, and in the homogenisation limit. The record names the escape
condition: "back-reaction, where the bath's state depends on the system path beyond a deterministic memory functional".
It also names the residual risk: known back-reaction limits are again Markov diffusions that an exogenous
multiplicative-noise model reproduces.

**SCOUT-0 reopening key** (`SATURATION_CERTIFICATION_01.md`). The overturn condition is "a bath outside class 𝓗
(nonlinear back-reaction) whose reduced law no exogenous model with system-independent law plus deterministic memory
can match, exhibited explicitly". BRI0 makes "exogenous model" precise (§5) and **interventional** (§3 – §4).

## 2. Central question

> When **one fixed environment** is tested across a preregistered family of interventions on the retained system, can
> genuine back-reaction produce an environment-force law that cannot be represented by a **shared** exogenous
> location-scale process plus deterministic causal response?

**Shared** is the operative word. A separate exogenous process may not be fitted to each intervention.

## 3. Single-protocol vacuity (PROVED HERE; established before any physics)

**Setting.** A single protocol is one potential V, drive u(t), preparation and horizon [0, T], for a system
m q̈ = −V′(q) + u(t) + F_env(t). Define the realised environment force from the reduced path:

  **F_env(t) := m q̈(t) + V′(q(t)) − u(t).**

**PROP BRI-V (single-protocol vacuity).** Assume the deterministic initial data (q₀, q̇₀) are part of the protocol, and
that the ODE m q̈ = −V′(q) + u + f has a unique solution on [0, T] for every admissible forcing path f, with a measurable
solution map q = Φ(q₀, q̇₀, f). Then the additive exogenous model with **M ≡ 0** and **ξ with Law(ξ) := Law(F_env)**
reproduces the reduced path law of that protocol exactly.

**Proof.**
1. In the true system the realised path solves the ODE with f = F_env, so q = Φ(q₀, q̇₀, F_env) pathwise.
2. In the model, q = Φ(q₀, q̇₀, ξ).
3. Equal laws of the input give equal pushforwards: Law(q) = Φ(q₀, q̇₀, ·)_# Law(F_env) in both. ∎

**Remarks.**
- (i) Under a clamp protocol (§4) the statement is immediate: one protocol has one force law, and ξ := F_q.
- (ii) With **random** initial data, the same construction needs ξ jointly distributed with (q₀, q̇₀) as F_env is. An
  "exogenous" ξ independent of the initial data matches whenever F_env is independent of them. The point is unchanged:
  a single protocol never discriminates.

**Consequence: NO SINGLE-PROTOCOL IDENTIFIABILITY CLAIM IS ADMISSIBLE IN BRI0.** Every BRI0 discriminator compares **one
shared** environment model across **several** interventions.

## 4. Primary protocol family 𝒫 = 𝒳 (trajectory clamp; frozen)

**𝒳** is the class of externally prescribed system trajectories q: [0, T] → ℝ with:
- C² regularity, so the force is defined;
- **q(0) = 0** and **q̇(0) = 0**;
- q, q̇ and q̈ bounded on [0, T];
- **no protocol-dependent environment preparation:** the environment is prepared from one fixed law before every
  protocol.

A protocol is the prescribed history q(·). The environment receives it through one fixed microscopic coupling and
returns a random force process F_q(·). This is measured interventionally as the force needed to hold the prescribed
trajectory, minus the known inertial, potential and drive terms. No hidden environment variable is read.

**Why clamping:**
- it removes ambiguity from different system potentials;
- every protocol starts from the same environment preparation;
- differences in Law(F_q) are therefore caused by how the environment reacts to different system histories.

**The at-rest clamp q ≡ 0 belongs to 𝒳.** It is the natural reference protocol (§8).

**Witness subset.** A later, owner-reviewed candidate charter must fix a **finite witness subset 𝒫_test ⊂ 𝒳 before**
any candidate result is computed. **No witness trajectories are chosen in Charter 0.**

**Causality of the environment** (a standing hypothesis of the parent class, §6). The law of the force restricted to
[0, t] depends on q only through q_[0,t].

## 5. Competitor ladder E₁ ⊂ E₂ ⊂ E_univ (preregistered)

Every competitor object is **shared across all of 𝒳**:
- it may depend on the prescribed path q_[0,t];
- it may **not** depend on a protocol label;
- **no** finite-memory, Gaussian, Markov or stationarity restriction is imposed anywhere.

**E₁ — additive exogenous.** F_q(t) = M_t[q] + ξ(t).
- M is an arbitrary deterministic causal functional (nonlinear, unlimited memory).
- ξ is one exogenous scalar process with one fixed law, independent of q.
- The same (M, Law ξ) serves every q ∈ 𝒳.
- P-17's class 𝓗 lies in E₁ (PROP BRI-H below).

**E₂ — location-scale exogenous (PRIMARY).** F_q(t) = M_t[q] + G_t[q]·ξ(t).
- M is an arbitrary deterministic causal functional.
- **G_t[q] > 0** is an arbitrary deterministic causal scalar functional.
- ξ is one exogenous scalar process with one fixed, protocol-independent law.
- The same (M, G, Law ξ) serves every q.
- E₁ is the case G ≡ 1.
- **Restriction kept:** location-scale separability. System history may change the deterministic response and the
  *amplitude* of the random force, but not the **shape / copula class** of the random component. G may **not** be
  enlarged into a nonlinear transformation of ξ, which would collapse E₂ toward E_univ.

**E_univ — universal causal exogenous representation.** F_q(t) = 𝔉_t[q, U].
- U is one exogenous random object with a protocol-independent law.
- 𝔉 is an arbitrary deterministic functional that is causal in q (it depends on q_[0,t] and U).
- (Law U, 𝔉) is shared across 𝒳. There is no separability restriction.

### 5.1 Inclusions (PROVED HERE)

**PROP BRI-L.** E₁ ⊂ E₂ ⊂ E_univ, and both inclusions are **strict**.

**Proof of inclusion.**
- E₁ ⊂ E₂: take G ≡ 1.
- E₂ ⊂ E_univ: take U := ξ (a path-valued random object with a fixed law) and 𝔉_t[q, U] := M_t[q] + G_t[q]·U(t). This is
  causal in q and shared across protocols.

**Strictness**, with i.i.d. standard normal ξ(t) on a time grid (or any non-degenerate symmetric law):
- **(a) E₂ \ E₁.** Take F_q(t) = (1 + q(t)²)·ξ(t), which is in E₂. Var F_q(t) = (1 + q(t)²)² depends on q. In E₁, the
  centred law Law(F_q − E F_q) = Law(ξ − E ξ) is independent of q (see PROP BRI-C1). So F ∉ E₁.
- **(b) E_univ \ E₂.** Take F_q(t) = ξ(t) + q(t)·(ξ(t)² − 1), which is in E_univ. At any t with q(t) ≠ 0, the
  standardised skewness is non-zero, while at q ≡ 0 it is 0. In E₂, standardised skewness is protocol-invariant (PROP
  BRI-C2). So F ∉ E₂. ∎

### 5.2 Exact characterisations under the clamp family (PROVED HERE)

These go beyond the owner's "sufficient witnesses" (§8): under clamping and finite second moments, E₁ and E₂ membership
are **decidable by a protocol-invariance property of a single, canonically defined process**.

**Standing assumption.** E F_q(t)² < ∞ for all q and t.

**Definitions.**
- Centred force: F̊_q(t) := F_q(t) − E F_q(t).
- Standardised force: Z_q(t) := F̊_q(t) / sd F_q(t) wherever sd F_q(t) > 0.

**PROP BRI-C1.** F ∈ E₁ ⟺ the finite-dimensional distributions (fdd) of F̊_q are identical for all q ∈ 𝒳; equivalently,
equal to those of the at-rest clamp F̊₀.

**Proof.**
- (⇒) For clamped q, M_t[q] is a deterministic number, so F̊_q = ξ − E ξ for every q.
- (⇐) Set ξ := F̊₀ and M_t[q] := E F_q(t). M is causal because Law(F_q(t)) depends only on q_[0,t] (environment
  causality, §4). Then F_q = M[q] + ξ in law, for every q. ∎

**PROP BRI-C2.** F ∈ E₂ (with G > 0) ⟺ both of the following hold:
- (i) the degenerate-time set D_q := {t : sd F_q(t) = 0} is the same set D for all q;
- (ii) on every finite tuple of times outside D, the fdd of Z_q are identical for all q, i.e. equal to those of Z₀.

**Proof.**
- (⇒) For clamped q, F_q(t) = m_q(t) + g_q(t)·ξ(t) with deterministic m_q and g_q > 0. So sd F_q(t) = g_q(t)·sd ξ(t),
  which gives (i) with D = {t : sd ξ(t) = 0}. Because g_q cancels, Z_q(t) = (ξ(t) − E ξ(t))/sd ξ(t), giving (ii).
- (⇐) Set ξ := Z₀ off D and ξ := 0 on D, with M_t[q] := E F_q(t) and G_t[q] := sd F_q(t) off D (G := 1 on D). Both are
  causal by environment causality. Then F_q = M[q] + G[q]·ξ in fdd for every q. ∎

**COROLLARY (witness completeness).**
- At any non-degenerate time, standardised marginal shape (skewness, excess kurtosis, all standardised cumulant ratios)
  and the copula of any finite force vector are invariant under positive componentwise affine maps. So they are E₂
  invariants, and any protocol difference in them is a **sufficient** witness of E₂ escape (owner §8).
- Conversely, when all one-time marginals are continuous, Sklar's theorem gives fdd(Z_q) = (standardised marginals,
  copula). Then **the standardised marginals plus the copulas together are also necessary**: E₂ membership holds
  exactly when they match the at-rest clamp at every finite tuple.
- If marginals have atoms, the copula is not unique and the criterion is stated directly on fdd(Z_q).

**Scope of the characterisation.**
- It concerns **all** of 𝒳. A finite 𝒫_test can only **witness escape**: one tested q whose fdd(Z_q) differs from
  fdd(Z₀), beyond controlled uncertainty, proves no shared E₂ representation exists.
- Finite tests can **never** certify E₂ membership. That needs an analytic representation (BRI-E2-THEOREM, §12).
- Without second moments, replace mean / sd by any equivariant location / scale functionals (e.g. median and
  interquartile range); the proofs carry over.

**PROP BRI-H (P-17's class 𝓗 lies in E₁ under clamping).**
- From the P-17 proof with q(0) = 0 and q̇(0) = 0, the clamped force is F_q(t) = F_free(t) − ∫₀ᵗ γ(t − s) q̇(s) ds.
- The counterterm and initial-slip terms are either deterministic in q or vanish at q(0) = 0.
- F_free is a function of the environment's initial data alone, with a q-independent law.
- So M_t[q] = −∫γ q̇ and ξ = F_free. ∎
- Likewise, a one-way additive driver (P-18 skew product with additive coupling) lies in E₁.

## 6. E_univ bounding theorem (PROVED HERE; BRI-UPPER earned analytically in Charter 0)

**Parent class (classical, deterministic, causal).**
- The environment state Y ∈ ℝ^d obeys Ẏ = B(Y, q(t)).
- B is locally Lipschitz in Y, continuous in q, and no solution blows up on [0, T] for q ∈ 𝒳 (e.g. a global Lipschitz
  or energy bound).
- The force on the system is F_env = C(Y, q), with C measurable.
- The initial state is Y(0) = Ψ(q(0), U), with **one fixed map Ψ** and one exogenous random object U. Since
  q(0) = 0 on 𝒳, this is Ψ(0, U).
- **Reciprocal (back-reacting) coupling is included:** B may depend on q arbitrarily, and the system may feed energy into
  Y.

**THEOREM BRI-UPPER.** For every q ∈ 𝒳, F_q(t) = 𝔉_t[q, U] := C(𝒴_t[q, U], q(t)), where 𝒴 is the solution flow. The
pair (Law U, 𝔉) is shared across 𝒳, and 𝔉 is causal in q. **Hence every such environment's interventional force law
lies in E_univ.**

**Proof.**
1. **Existence and uniqueness:** Picard–Lindelöf plus the no-blow-up hypothesis give a unique solution Y on [0, T] for
   each (q, U).
2. **Causality:** if q and q′ agree on [0, t], then both solutions solve the same initial-value problem on [0, t] with
   the same Y(0), so by uniqueness they agree on [0, t]. Hence 𝒴_t depends only on q_[0,t].
3. **Measurability:** the solution depends continuously on the initial data, and Ψ is measurable, so U ↦ 𝒴_t[q, U] is
   measurable.
4. **Shared:** B, C, Ψ and Law U are fixed objects of the environment, not of the protocol.

∎

**Remark (KNOWN technique; scope statement).** E_univ is in fact the class of **all** non-anticipating protocol-indexed
force-law families, not just the deterministic-environment ones:
- In discrete time on standard Borel spaces, the sequential randomisation (transfer) lemma writes any non-anticipating
  kernel family as F_t = f_t(q_[0,t], F_{<t}, V_t), with V_t i.i.d. uniform and V independent of everything. This uses
  one shared U = (V_t)_t.
- Continuous time needs regularity conditions and is **not** claimed here.
- Stochastic (non-deterministic) classical environments are therefore also inside E_univ in the discrete-time setting.
- **Quantum environments are not covered** by this charter.

**Consequences (recorded):**
1. No classical, causal, back-reacting environment escapes E_univ.
2. This is a representation theorem, not a numerical result.
3. Back-reaction alone cannot establish primitive randomness or a distinct ontology: at the E_univ level, "the
   environment reacts to my path" and "a fixed random object processed by a causal response map" are the **same
   mathematical object**.
4. All scientific content lies in the boundary **E₂ vs E_univ**: whether the reaction stays location-scale separable.

## 7. Why E₂ is the primary rung

- **Escaping E₁ is weak.** History-dependent noise amplitude (multiplicative forcing) is enough, and that is
  known open-system structure. Control C1 (§10) is expected to escape E₁.
- **Escaping E₂ is non-trivial.** It means the system history changes the **shape or inter-time dependence** (the
  copula) of the environment's random response, beyond deterministic location and scale modulation of one common
  exogenous process. By BRI-C2, this is exactly a change in the law of the standardised clamped force.
- **Escaping E_univ is impossible** in the declared parent class (BRI-UPPER).

**E₂ ESCAPE IS THE PRIMARY LANE-3 TARGET.**

## 8. Location-scale invariants (preregistered sufficient witnesses)

For clamped protocols, M_q(t) and G_q(t) are deterministic, so E₂ imposes exact shared-law constraints.

**Invariants:**
- for any time t with non-zero variance: the standardised skewness, standardised excess kurtosis, and every standardised
  one-time cumulant ratio that exists are E₂-invariant;
- for any finite ordered tuple t₁ < … < t_k: the copula of (F_q(t₁), …, F_q(t_k)) is invariant under positive
  componentwise affine maps.

**Protocol differences** in standardised marginal shape or in the force-vector copula, with sampling / analytic
uncertainty controlled, are **sufficient** witnesses that no common positive location-scale (E₂) representation exists.
- They are **not necessary** in general.
- Under clamping with continuous marginals they **are** jointly necessary (COROLLARY, §5.2).
- Failing to observe them in a **finite** 𝒫_test never proves E₂ membership.

## 9. Bath-size axis N_B (mandatory for every future candidate)

**Hypothesis registered (CONJECTURE, not a result):**
- A finite anharmonic environment can be driven out of its initial statistical state by energy transferred from the
  system.
- For increasingly large, weakly perturbed environments, this effect may shrink toward an effectively fixed-reservoir
  (linear-response) law.
- **Any E₂ escape may be a finite / mesoscopic-bath effect that vanishes as N_B grows.**

**Requirements on any future candidate campaign:**
- at least **four** increasing bath sizes, fixed before computing the primary escape statistic;
- the microscopic coupling normalisation as a function of N_B stated explicitly;
- a report of whether each E₂-escape witness grows, saturates or vanishes with N_B;
- **no single-bath-size positive earns a general result.**

**The N_B grid is not chosen in Charter 0.**

## 10. Future model roles (registered structurally; NOT run; no Hamiltonian chosen)

**C1 — nonlinear-coupling harmonic control.** A harmonic bath with system coupling through A(q).
- Expected structure: eliminating the bath gives a nonlinear deterministic memory plus a multiplicative free-force term
  A′(q(t))·F_free(t).
- The role is to show that "multiplicative noise" alone is **not** an E₂ escape.

**AUDIT FINDING (owner decision required before any C1 computation).** With **G > 0** as chartered, C1 lies in E₂ only
when A′(q(t)) keeps a **strict constant sign** along every protocol considered. If A′(q(t)) changes sign or vanishes
along a clamped trajectory:
- the factor A′(q(t)) either flips the sign of a centred Gaussian F_free(t) at some times, or makes F_q(t) degenerate
  there;
- either changes the copula, or the degenerate-time set D_q, relative to the at-rest clamp;
- so C1 would **leave E₂ through a trivial sign / zero artefact.** That would be a false "escape" from known structure.

Two clean options, to be ruled on **before** computation. Neither is adopted silently here.
- **(a)** Keep E₂ (G > 0) and require C1 (and any candidate) to use monotone A, or a 𝒫_test along which A′(q(t)) never
  vanishes and never changes sign.
- **(b)** Define **E₂±**: G_t[q] ∈ ℝ, sign and zero allowed. BRI-C2 then becomes invariance of Z_q **up to a
  deterministic time-dependent sign**, with D_q allowed to vary.

**X1 — finite anharmonic back-reacting candidate.** A finite nonlinear oscillator environment with reciprocal coupling to
q. The system changes the environment's energy distribution and may therefore change the **shape** of the force
statistics it gets back. This is the primary future escape candidate.
- Its Hamiltonian, size grid and 𝒫_test require a **separate owner-reviewed candidate charter**, after the competitor
  classes are accepted.

## 11. Information firewall

**A valid E₂ escape may not come from:**
- different environment preparations per protocol;
- environment parameters changed between protocols;
- a ξ-law fitted per protocol;
- M / G functional forms changed between protocols;
- protocols selected after inspecting candidate results;
- conditioning on hidden environment variables unavailable to the interventional description;
- a finite-memory restriction imposed on the competitor;
- a Gaussianity restriction imposed on the competitor;
- a sign / zero artefact of G (§10 audit), unless option (b) is adopted.

The environment law and preparation are **one shared object** across the full protocol family.

## 12. Grades / terminals (chartered now; none assigned)

| ID | grade | meaning |
|---|---|---|
| **BRI-E1** | **E₁ ESCAPED; E₂ NOT ESCAPED** | history-dependent noise amplitude is required, but location-scale separability survives. Expected for C1 (under §10 option (a) or (b)). Known open-system structure; **not** GRUT distinctiveness |
| **BRI-E2+** | **E₂ ESCAPED AT FINITE BATH SIZE** | a preregistered interventional witness proves that no shared E₂ representation exists for the tested candidate and protocol family. **Must include the N_B scaling result.** Reading: **BACK-REACTION IS IDENTIFIABLE RELATIVE TO THE LOCATION-SCALE EXOGENOUS CLASS.** Not primitive randomness, not a unique ontology, not a GRUT-specific law |
| **BRI-E2−** | **NO E₂ ESCAPE FOUND IN THE TESTED BACK-REACTING CLASS** | all preregistered witnesses are consistent with E₂. Not a theorem unless an analytic E₂ representation is constructed |
| **BRI-E2-THEOREM** | **TESTED BACK-REACTING FAMILY PROVED TO LIE IN E₂** | requires an analytic shared representation (e.g. via BRI-C2), not numerical fitting. Extends P-17 / P-18 non-identifiability to the new family |
| **BRI-UPPER** | **E_UNIV NON-IDENTIFIABILITY BOUND** | THEOREM BRI-UPPER (§6). **Earned analytically in Charter 0** (INTERNALLY PROVED / NOT EXTERNALLY REVIEWED) |

## 13. Audit of the class definitions (hidden universality / arbitrary restrictions)

| check | finding |
|---|---|
| E₁ / E₂ secretly universal? | **No.** Both inclusions are strict (PROP BRI-L), and E₂ has exact checkable invariants (BRI-C2) |
| E_univ secretly restrictive? | **No.** It contains every classical deterministic causal environment (BRI-UPPER), and in discrete time every non-anticipating family (KNOWN randomisation). Its only exclusion is anticipation |
| Hidden finite-memory / Gaussian / Markov / stationarity restriction? | **None.** M and G are arbitrary causal functionals; ξ has an arbitrary path law |
| Causality of M, G a hidden restriction? | **No.** It is automatically satisfiable for causal environments (proofs of BRI-C1 / C2) |
| Single-protocol vacuity re-entering? | **Excluded.** Every class is shared across 𝒳 and contains the at-rest reference clamp |
| Scalar ξ, scalar G | Matches the scalar q and scalar force of this charter. A vector system would need a matrix G; **not** chartered |
| **G > 0 sign / zero artefact** | **Open: owner decision** (§10 options (a) / (b)) |
| Moment assumption | The characterisations use second moments, replaceable by equivariant location / scale functionals (§5.2) |
| Clamp realism | The force is measured interventionally under prescribed q, with no hidden variables. Free-evolution (unclamped) reduced data are **not** the primary object; any later use must re-derive the classes |
| Quantum environments | Out of scope |

## 14. Relation to GRUT (firewall)

**Even a BRI-E2+ would establish only this:** some finite back-reacting environments have interventional
reduced-force laws outside a shared location-scale exogenous representation.

**It would not establish:**
- a uniquely GRUT law or a GRUT prediction;
- primitive randomness;
- consciousness;
- TRUE COMPRESSION;
- escape from E_univ (impossible, BRI-UPPER).

**Its value:** it would identify a law-level invariant beyond the P-17 / P-18 separable quotient class, i.e. what
information genuine back-reaction can make observable.

## 15. Separations and hard stop

**Phase 0** (the SCOUT-0 owed independent reproductions) is **not** performed on this branch. It is verification /
packaging work under the frozen saturation programme, with its own future branch and ruling. **BRI0 is Lane 3 only.**

**Hard stop (Charter 0):**
- proved here: single-protocol vacuity, E₁ ⊂ E₂ ⊂ E_univ (strict), the E₁ / E₂ clamp characterisations, P-17's 𝓗 ⊂ E₁,
  and BRI-UPPER;
- class audit done;
- **no bath calculation, no anharmonic Hamiltonian, no 𝒫_test, no N_B grid.**

**Questions for owner review:**
1. G sign / zero (§10): option (a) or (b)?
2. Accept the clamp family 𝒳 as the primary 𝒫, with the at-rest clamp as the reference protocol?
3. Accept BRI-C2 as the exact operational meaning of E₂ membership under clamping, so that future witnesses target
   standardised-force fdd (marginal shape plus copula)?
4. Accept BRI-UPPER as earned (analytic, INTERNALLY PROVED / NOT EXTERNALLY REVIEWED)?


## §R — BRI0 SCOPE REPAIR 01 (owner review of `89b8236`; no candidate selected or run)

### BRI-O1 — Clamp family accepted

𝒫 = 𝒳 (§4) is BRI0's primary interventional protocol family. The scientific object is the **interventional
environment-force law**, not passive free-evolution reduced data. The following are kept:
- common q(0) = 0 and q̇(0) = 0;
- one shared environment preparation;
- a causal response;
- a finite 𝒫_test, to be preregistered later.

**Reference protocol.** The at-rest clamp q ≡ 0 is the default reference, but it may serve as the standardisation
reference **only at times where its force variance is non-zero** (BRI-O4). No witness trajectories are selected.

### BRI-O2 — Primary competitor: E₂± (signed location-scale / affine exogenous class)

**E₂±:** F_q(t) = M_t[q] + G_t[q]·ξ(t), where:
- M is an arbitrary deterministic causal functional;
- **G is an arbitrary deterministic causal functional with G_t[q] ∈ ℝ \ {0}**;
- ξ is one shared, protocol-independent path law.

Deterministic **sign** changes of the noise amplitude are allowed. **Arbitrary nonlinear transformations of ξ are not**:
the restriction remains affine (location-scale) separability.

The ladder becomes **E₁ ⊂ E₂± ⊂ E_univ**.
- **E₁:** deterministic nonlinear memory plus one fixed additive random law.
- **E₂±:** additionally, arbitrary deterministic causal magnitude and sign modulation.
- **E_univ:** an arbitrary causal transformation of one fixed random object.

**Why signs are allowed.** A deterministic factor such as A′(q(t)) changing sign does not change the underlying random
law. It is ordinary multiplicative response and must not count as the sought escape. The C1 control (§10) is therefore
expected to lie in E₂± wherever its multiplicative coefficient is non-zero. **The §10 sign audit is resolved in this
way. The positive-G class E₂ of §5 is retired as primary.**

### BRI-O3 — Zero-amplitude firewall

G = 0 is **not** allowed in E₂±. A zero multiplier makes the stochastic force degenerate and is not an invertible affine
map.

**Degeneracy set:** D_q := {t : Var F_q(t) = 0}.

**Exact consequence (PROVED HERE).** In E₂±, Var F_q(t) = G_t[q]²·Var ξ(t) with G ≠ 0. So D_q = D_ξ := {t : Var ξ(t) = 0}
for every q: **membership requires a common degeneracy set.**

**Rules for future finite witness campaigns:**
- primary witness times are chosen before results;
- the reference force must be non-degenerate at those times;
- a protocol-dependent zero of the force amplitude is reported as a **degeneracy change (BRI-DEG)**. Strictly, it is
  outside E₂±, but **it alone never earns BRI-E2+.**

**BRI-E2+ requires** a **non-degenerate** standardised-force law difference that cannot be removed by deterministic
affine location, magnitude and sign modulation. A zero of A′(q) therefore cannot become a trivial "success".

### BRI-O4 — BRI-C2± (replaces BRI-C2; PROVED HERE)

**Assumptions:** clamped protocols, E F_q(t)² < ∞, and environment causality (§4).

**Definition.** At non-degenerate times, Z_q(t) := (F_q(t) − E F_q(t)) / sd F_q(t).

**PROP BRI-C2±.** F ∈ E₂± ⟺ both of the following hold:
- (i) D_q = D, a common set, for all q ∈ 𝒳;
- (ii) for each q there is a **deterministic causal sign function** s_q(t) ∈ {+1, −1}, depending on q only through
  q_[0,t], such that the fdd of s_q·Z_q, on finite tuples of times outside D, are **common across all q ∈ 𝒳**.

**Proof (⇒).**
- For clamped q, F_q(t) = m_q(t) + g_q(t)·ξ(t), with deterministic, causal m_q and g_q, and g_q(t) ≠ 0.
- D_q = D_ξ (BRI-O3).
- Off D, Z_q(t) = sgn(g_q(t))·Z_ξ(t), where Z_ξ is the standardised ξ.
- Set s_q(t) := sgn(g_q(t)). This is causal because G is causal, and s_q·Z_q = Z_ξ has a q-independent law. ∎

**Proof (⇐), the converse construction.**
1. Let ℒ be the common fdd family of s_q·Z_q (off D).
2. Let ξ be a process with fdd ℒ off D and ξ(t) := 0 on D. Its law is protocol-independent by (ii).
3. Set M_t[q] := E F_q(t).
4. Set G_t[q] := s_q(t)·sd F_q(t) off D, and G_t[q] := 1 on D, so that G ≠ 0 everywhere.
5. M and G are causal: the law of F_q(t) depends only on q_[0,t] (environment causality), and s_q is causal by
   hypothesis.
6. Off D, M_t[q] + G_t[q]·ξ(t) has the same fdd as E F_q + sd F_q·s_q·(s_q Z_q) = F_q, because s_q² = 1. On D both
   sides equal the deterministic value E F_q(t).
7. Hence F_q has the same fdd as M[q] + G[q]·ξ for every q, with one shared Law(ξ). ∎

**Reference-independent form.** The statement is **common-orbit invariance across the family**; it does not depend on
any particular reference. The at-rest clamp is only a convenient representative of ℒ where it is non-degenerate. If a
later 𝒫_test uses times at which the at-rest force is degenerate, a different carrier reference must be fixed in that
candidate charter before results.

**Reflection-orbit form and its exact scope (correction to the owner text's "equivalently").** For a finite ordered
tuple τ = (t₁ < … < t_k) outside D, let L_q(τ) be the law of (Z_q(t₁), …, Z_q(t_k)), and let O(τ) be the orbit of the
reference law under the group R_k = {diag(s₁, …, s_k) : s_i ∈ {±1}}. Then:

- **(a) Necessity.** F ∈ E₂± ⟹ L_q(τ) ∈ O(τ) for every q and every τ. *Proof:* restrict the causal signs of (ii) to τ.
  ∎ **Consequently, a reflection-orbit violation at any single tuple is a sufficient witness of E₂± escape.** This is
  the primary witness (BRI-O5).
- **(b) Global, possibly non-causal signs.** If L_q(τ) ∈ O(τ) for every τ, then a **global** sign function
  s_q: [0, T] \ D → {±1} exists such that s_q·Z_q has the reference fdd. *Proof:* for each τ, the set of valid sign
  vectors is non-empty and finite. Marginalisation commutes with coordinate reflection, so valid vectors for τ′ ⊃ τ
  restrict to valid vectors for τ. These sets form an inverse system of non-empty compact sets, so by compactness of
  {±1}^{[0,T]} (Tychonoff) a consistent global selection exists. ∎
- **(c) The per-tuple orbit condition does NOT imply E₂± membership.** E₂± requires the selection in (b) to be
  **causal**, and orbit consistency does not force that.
  - **Counterexample (PROVED HERE).** Let ζ be standardised and symmetric with E ζ³ = 0 (e.g. standard normal), and
    η := ζ + ζ² − 1. Let the reference law at (t₁, t₂) be that of (ζ, η)/standardisation.
  - Take two protocols q and q′ that agree on [0, t₁] and differ on (t₁, t₂], with Z_{q′} = (ζ, η) and
    Z_q = (−ζ, η) (standardised), realised causally by an environment whose t₂-response sign depends on q after t₁.
  - Both lie in the reference orbit: q′ via (+, +) and q via (−, +). The latter is the **only** valid sign vector,
    because (−ζ, ±η) has law (ζ, ±(−ζ + ζ² − 1)), and neither sign reproduces (ζ, ζ + ζ² − 1).
  - So s_q(t₁) = −1 and s_{q′}(t₁) = +1 although q and q′ agree on [0, t₁]: no causal sign function exists, and F ∉ E₂±.
  - The family satisfies every per-tuple orbit condition, yet is outside E₂±.
- **(d) Registered consequence.** "Causal-sign inconsistency" is a **second, distinct** escape mode: orbit-consistent at
  every tuple, but with no causal sign selection. It is a genuine E₂± escape, but it involves future-dependent sign
  bookkeeping rather than a change in the standardised law's shape. **It may not earn BRI-E2+ without a separate owner
  ruling** in the candidate charter. The **primary** BRI-E2+ witness remains the reflection-orbit violation (a).

### BRI-O5 — Witness interpretation under E₂± (replaces §8)

**Invariant under E₂±**, at non-degenerate witness times, with no qualification:
- standardised **even** cumulants and their ratios (e.g. excess kurtosis);
- **absolute values** of standardised odd cumulants (e.g. |skewness|);
- any functional of the standardised law that is invariant under reflection;
- for a finite vector: the dependence law **modulo coordinatewise reflections**. Examples are |correlation| and
  the absolute value of every mixed standardised cumulant: a reflection multiplies κ(Z₁^{a₁}, …, Z_k^{a_k}) by
  ∏ s_i^{a_i}, so only its sign can change. The copula is invariant only up to the
  reflection action, under which a coordinate's u ↦ 1 − u.

**Not invariant:** the sign of standardised skewness, the sign of correlations, and the raw copula. **A raw copula-sign
or skew-sign change never counts as escape.**

**Primary mathematical witness:** the standardised finite-dimensional force law of some protocol lies **outside the
reference orbit** under every allowed coordinatewise reflection (BRI-O4 (a)). Marginal cumulants and copula diagnostics
are **secondary** ways of proving this.

### BRI-O6 — Exact-characterisation scope

BRI-C1 and BRI-C2± are exact under the clamp setup **as statements about finite-dimensional distributions.**
- fdd determine the law on the product (cylinder) σ-algebra.
- They determine the **path-space law** when the force processes have, for every q, versions with **continuous** (or
  càdlàg) paths. The cylinder σ-algebra then generates the Borel σ-algebra of C[0, T] (resp. D[0, T] with the Skorokhod
  topology), so fdd equality is path-law equality.
- For the converse construction to yield a path-regular ξ, it suffices that m_q(t) = E F_q(t) and sd F_q(t) are
  continuous and bounded away from 0 off D on the times used. Otherwise the representation is stated at the fdd level.
- **Fdd equality is not silently upgraded** to equality of arbitrary path-space measures.

**Finite tests:** a finite 𝒫_test can **prove escape**. It **cannot** prove membership in E₁ or E₂±. Membership requires
the full-family analytic characterisation or an explicit shared representation (BRI-E2-THEOREM).

### BRI-O7 — BRI-UPPER accepted, at its proved scope

**BRI-UPPER — E_UNIV NON-IDENTIFIABILITY BOUND: ACCEPTED, EARNED (INTERNALLY PROVED / NOT EXTERNALLY REVIEWED).**

**The proved continuous-time theorem:** every classical, deterministic, causal environment in the declared parent class
(Ẏ = B(Y, q), F = C(Y, q), Y(0) = Ψ(q(0), U), with the §6 regularity) lies in E_univ. Reciprocal and energy-absorbing
back-reaction is included.

**Kept separately:** the standard discrete-time randomisation result for general non-anticipating stochastic kernels on
standard Borel spaces (KNOWN technique).

**Wording repaired.** §6 consequence 1 ("No classical, causal, back-reacting environment escapes E_univ") and §13's
E_univ row are superseded by this narrower statement:
- **No environment in the declared continuous-time deterministic parent class escapes E_univ.**
- In discrete time, general non-anticipating stochastic environments are also representable.
- **A continuous-time statement for general stochastic environments is not claimed.** It would need separately stated
  regularity assumptions and a proof.
- Quantum environments remain out of scope.

### BRI-O8 — Strictness of the ladder under E₂± (PROVED HERE)

- **E₂± \ E₁.** F_q(t) = (1 + q(t)²)·ξ(t), with ξ i.i.d. standard normal. It is in E₂± (G > 0 is allowed). It is not in
  E₁, because Var F_q(t) = (1 + q(t)²)² varies with q, while E₁ requires a q-invariant centred law (BRI-C1). *Survives
  unchanged.*
- **E_univ \ E₂± (reflection-safe).** F_q(t) = ξ(t) + q(t)·(ξ(t)² − 1), with ξ i.i.d. standard normal. It is in E_univ.
  - Its third central moment is 6q + 8q³, where q = q(t), and its variance is 1 + 2q².
  - So the **absolute** standardised skewness is 0 at q(t) = 0 and **non-zero** whenever q(t) ≠ 0.
  - |skewness| is reflection-invariant (BRI-O5), so it is E₂±-invariant, and the at-rest clamp and any q with q(t) ≠ 0
    give different values.
  - Hence F ∉ E₂±. ∎

So **E₁ ⊊ E₂± ⊊ E_univ.**

### Accepted scientific reading

> A future **BRI-E2+** means: back-reaction changes the **standardised random response** in a way that **no shared
> deterministic affine modulation (location, magnitude and sign) of one exogenous process** can represent.

**Back-reaction itself is never the discriminator.** At the universal level it is always representable (BRI-UPPER). The
only interesting question is whether the response stays in the affine orbit M[q] + G[q]·ξ.

### Candidates (unchanged hold)

- **C1** remains a future structural control. Sign changes of A′(q) are now allowed; zeros are excluded from primary
  witness times and from primary escape grading (BRI-O3). No C1 computation, and A(q) is not chosen.
- **X1:** the exact anharmonic Hamiltonian, 𝒫_test, N_B grid, coupling normalisation and escape statistics are **not**
  chosen. They belong to the next candidate charter.

### Grades: additions under Repair 01

| ID | grade | meaning |
|---|---|---|
| **BRI-DEG** | **DEGENERACY CHANGE ONLY** | a protocol-dependent zero-variance set. Outside E₂± but **never** BRI-E2+ by itself |
| **BRI-CSI** | **CAUSAL-SIGN INCONSISTENCY ONLY** | orbit-consistent at every tuple, but with no causal sign selection (BRI-O4 (c)). An E₂± escape whose grade requires an owner ruling; **not** BRI-E2+ by default |

BRI-E1, BRI-E2+, BRI-E2−, BRI-E2-THEOREM and BRI-UPPER keep their meanings, with "E₂" read as **E₂±** throughout.
