# DA0 CORRECTION LEDGER (additive)

**Rule:** logs are kept as emitted. Where document wording differs, this ledger takes precedence.

## DA0 CHARTER REPAIR 01 (owner review of `bc9d3d2`; applied before any computation)

| ID | item | correction |
|---|---|---|
| **DC1-01** | cycle decomposition is not canonical | J, the cycle space and edge affinities are canonical. Any particular fundamental-cycle decomposition is **BASIS-PRICED / DIAGNOSTIC ONLY**, unless a basis-independent selector is proved. A "canonical cycle" must be invariantly identified from D; otherwise the terminal is **FLUX STRUCTURE ONLY — NO CANONICAL CYCLE**. EP = ½ Σ J log[(π_x G_xy)/(π_y G_yx)] ≥ 0 only for mutually supported edges; one-way edges need an explicit extended convention |
| **DC1-02** | real-part gap is not a sufficient non-normal cut certificate | Use r(λ) = −Re λ. Spectral sets must be conjugation-closed. Report E_N, ‖E_N‖, K_N(Γ) and conditioning. Perturbative stability needs the contour-resolvent product to vanish. Exploding conditioning ⇒ **SPECTRAL SEPARATION PRESENT — NONNORMAL RECOVERY UNSTABLE**, not a recovered macrostructure |
| **DC1-03** | orientation bit vs supplied non-reversible structure | P vs P* is one supplied bit. A and J are part of supplied D and carry more. The positive grade is **CONDITIONAL DERIVATION OF DIRECTED MACROSTRUCTURE FROM SUPPLIED NONREVERSIBLE D**. Test: covariance under P → P*; matched families with the same π and S but different A may differ (**A-PRICED**). R1 becomes a matched pair |
| **DC1-04** | partition recovery and directedness are separate | C1-P (partition without supplied inputs) and C1-D (directed macro-generator, conditional on C1-P) are separate questions. C1-D without C1-P earns no directed endogenous partition. For Q vs Q*: covariant partition, sign-reversed current and affinities, no supplied cycle basis. Final information accounting: derived / supplied / not earned |

## C1 note (process; no repair)

These changes were made **before the logged C1 run**:
- `macro_cycle` crashed on a round-off-level macro current in R1b. A declared computational zero (ZERO_J = 1e-13) and a
  guard requiring a cycle of length ≥ 3 were added.
- The Kato bound was first coded in the crude form ρR₀·δ/(1 − δ), which is O(√(η/g)). It was replaced by the
  residue-expanded form 2eaη + (e + ρa)δ²/(1 − δ), which is O(η/g) and is the one stated in THEOREM C1-T. Both are valid.
- The macro-current error constant was corrected from 2·q_max·m to **4**·q_max·m, as in PROP C1-D. All logged m are 0.
- An audit mapping blind macro-cycle labels to hidden labels, the R1 matched-pair check, and the DC1-02 defect-ring probe
  were added.
- The ad-hoc defect-ring probe that preceded this was run once outside the log. Its values are reproduced in the logged
  probe.

## DA0 C1 REPAIR 01 (owner review of `feb6bc3`)

The C1 terminal is preserved: **NONREVERSIBLE PARTITION DERIVATION**, graded **CONDITIONAL DERIVATION OF DIRECTED
MACROSTRUCTURE FROM SUPPLIED NONREVERSIBLE D — A-PRICED**. The original log is preserved.

| ID | item | correction |
|---|---|---|
| **C1R-01** | rank-K proof in THEOREM C1-T | The ‖E − E₀‖ < 1 argument is withdrawn; the hypotheses do not make the explicit bound < 1 at finite N. Replaced by a **homotopy** G_t = G₀ + tG₁ on the fixed circle |z| = ρ. For all t, ℓ_t ≤ ℓ < ρ and δ_t ≤ δ < 1, so the circle stays in the resolvent set and the rank is constant: rank E₁ = rank E₀ = K. The projector bound is kept for convergence only |
| **C1R-02** | macro-current recovery condition | N4 is replaced by **N4a** (j_min,N = min nonzero \|J_𝔅\|) and **N4b** (q_max,N·m_N / j_min,N → 0). Sufficient specialisation when j_min ≥ c_J η: q_max·η/g² → 0. N3 alone (η/g → 0) does not suffice when g → 0, e.g. η = g^{3/2}. **Audit:** R1a, R2 and R5 satisfy N4a / N4b asymptotically (g = O(1), η ∝ 1/M); **R4 fails both** (j_min/η → 0; q_max·η/g² grows ∝ M) |
| **C1R-03** | zero currents | "Support and sign pattern are recovered" is withdrawn. Under N4b, non-degenerate nonzero entries have their signs and relative magnitudes recovered. Hidden zeros are shown only to converge to zero, and an exact zero is never inferred from a vanishing sequence. Structural zeros by construction (R2's non-adjacent basins; R5's pair (1, 3)) are family-specific facts |
| **C1R-04** | macro EP | General macro-EP recovery is removed from PROP C1-D (log traffic ratios can be ill-conditioned). Kept: blind EP_Q, the numerical values, EP(G*) = EP(G), EP_Q ≤ EP |
| **C1R-05** | time-reversal certification | "G* satisfies N1 – N4 with the same constants" is withdrawn: the decoupled operator of G* need not be G₀*. Kept: J(G*) = −J(G) exactly, and J_Q(G*) = −J_Q(G) for any fixed partition. If G and G* are each certified for the same asymptotic partition, the recovered partitions agree and the currents reverse. The numerical covariance is evidence |
| **C1R-06** | K2 comparator scope | The theorem-like wording is replaced by a comparator observation: in the ideal first-order Markov description with predictively distinct macro-rows, causal states identify the macrostates. Verdict unchanged: **NOT EXCEEDED** |
| **C1R-07** | multi-cycle control | **R5 added**: K = 4, complete macro graph (cycle-space dim 3), supplied superposition of two independent circulations. No spanning tree, cycle basis or decomposition is used. Partition recovered blindly (misclassified 0 at M = 10 … 80); canonical J_Q reported as a matrix; J_Q is not a multiple of any single cycle (5 support edges on 4 vertices); time-reversal covariance holds to ≤ 1.4e-16. Terminal: **NONREVERSIBLE PARTITION DERIVATION + CANONICAL MACRO FLUX; NO CANONICAL CYCLE DECOMPOSITION** |
| **C1R-08** | non-normal instability | Remains an **untriggered possible terminal**; its absence is not claimed in general. N2 stays priced. No pathological family is required before C2 |

**Applied to:** `C1_NONREVERSIBLE_MACROSTRUCTURE.md`, `DA0_ZOOM_OUT_01.md`, `DA0_STATUS.md`. New numerics:
`c1/c1_repair.py`, `c1/c1_repair.log`.

## Owner ruling on `1c3e9ad` (no repair)

- **C1 ACCEPTED AFTER REPAIR 01. C2 APPROVED TO CHARTER.**
- **Wording discipline (standing).** At the logged finite sizes, the N4b sufficient quantity q_max·η/g² is still > 1 for
  R1a and R5 and ≈ 1 for R2 at the largest M. Their certification is **asymptotic, from the analytic 1/M scaling**, not
  from having entered a numerically small-error regime.

## DA0 C2 CHARTER REPAIR 01 (owner review of `92393d1`; before any C2 computation)

| ID | correction |
|---|---|
| **C2R-01** | Bounded description alone does not separate emergent from encoded architecture. **C2-F3′** requires both a DESCRIPTION TEST and an ARCHITECTURE-EXPLICITNESS TEST (no tree / level / basin / block / layer-count / level-indexed coupling / isomorphic latent variable). Otherwise the family is ARCHITECTURE EXPLICIT IN FAMILY, and recovery is graded ARCHITECTURE SUPPLIED BY FAMILY |
| **C2R-02** | GREM / CREM rejected as strong E1 (hierarchy explicit) and moved to the P1 comparator class. E1 is chosen after a recorded audit: **East model** (C2 §6a) |
| **C2R-03** | **C2-F6:** forbidden are supplied macro-units, subject units, subsystem blocks and recovery grains. Allowed but priced are microscopic sites / spins, locality and local rules. Microscopic sites may not be relabelled as emergent architecture |
| **C2R-04** | C2-A1 is split into **A1a** (abstract: K, p_min, C_V, s) and **A1b** (family → s). s = O(η/g) with a K-independent constant is not assumed. K-dependence of every constant is reported |
| **C2R-05** | Primary metric **m_rel = min_σ max_i π(B_i Δ B̂_σ(i))/π(B_i)**, with the permutation chosen on the audit side only. It is bounded explicitly from the idempotent errors. s·K^{3/2} → 0 is a candidate sufficient condition only |
| **C2R-06** | The "gap-free ⇒ maybe K-uniform" inference is withdrawn. Gap-free concerns coefficient separation, not dimension. Whether K^{3/2} is a proof artefact remains open |
| **C2R-07** | Nesting metric ν; no forced nesting, no hierarchical-clustering optimisation; exact numerical nesting is evidence only |
| **C2R-08** | **C2-F7** disorder firewall |
| **C2R-09** | H1 and H2 get dual grades (… and ARCHITECTURE SUPPLIED BY FAMILY). H1 must account for n·λ_s/g → 0 for the 2^n slow sector to separate |
| **C2R-10** | CONJECTURE C2-B is stated only under an explicit locality / rate / temperature class, and is not inferred from the product bound spread ≥ R^d |

## C2 note (process; no repair)

- **Exploratory runs.** The C2 script was run three times while adding (i) the A2b breakdown sweep, (ii) the P1 ζ-scaling
  check, and (iii) the East Δ_HS, L = 10 and modal-configuration diagnostics. The logged run is the final one. Earlier
  outputs were identical on the shared parts.
- **Numerical-diagnostic threshold flagged.** The East depth count classifies a cut as "diverging" when the fitted
  exponent in 1/q exceeds 0.5. This is a numerical diagnostic, not a certificate. Under C0-P4 it is never a modelling
  input.
- **Upper-bound m_rel.** m_rel is computed with a Hungarian (sum-optimal) matching, which upper-bounds the min-max m_rel
  of C2R-05.
- **Improved Lemma G5-3′.** Introduced in C2. It **supersedes** the anticipated s·K^{3/2} sufficient condition in the C2
  charter §3, with s·K^{3/4} under balance. The C2 charter's §3 remains as written (preregistration text).

## C2 SCOPE REPAIR 01 (owner review of `b403594`)

**Owner ruling:**
- **C2-A PASSES:** blind recovery with K_N → ∞ is proved under the stated regularity, with the balanced sufficient
  condition s_N·K_N^{3/4} → 0 (Lemma G5-3′ accepted).
- Growing architecture from a simple rule: **NOT FOUND**.
- Unbounded partition depth without supplied scaling: **NOT FOUND**.
- Intrinsic K-growth obstruction: **OPEN**.
- TRUE COMPRESSION: **0**.
- **C3 not opened.**

| ID | correction |
|---|---|
| **C2S-01** | "No K supplied" is too broad. Correct statement: **no K is supplied to the partition-recovery map once the canonical slow space is certified. The family-level certification of that slow space may use the hidden proof architecture.** The numerical controls used the construction rank to extract the first K modes. Firewall: blind recovery, conditional rank / cut certification (RA3-03) |
| **C2S-02** | The East terminal is evidence, not a theorem: **"No endogenous growing partition architecture found; the tested spectral hierarchy requires q → 0, consistent with known East theory."** It is not claimed that fixed-q East can never produce a partition-like hierarchy for all L |

**Applied to:** `C2_GROWING_ARCHITECTURE.md`, `DA0_ZOOM_OUT_02.md`, `DA0_STATUS.md`. **C2-E2** is preregistered in the C2
charter §6b; it is **not run**.

## C2S-03 / C2S-04 and C2 CLOSURE (owner review of `c394e03`)

| ID | correction |
|---|---|
| **C2S-03** | The registered E2 sequence 2 × 2 × L (L = N/4) keeps both transverse widths fixed. It is a finite-width spin-glass tube / ladder, effectively one-dimensional in the limit, not the cubic L × L × L sequence whose T_c ≈ 1.10 (finite-size scaling up to L = 40) was cited. T = 0.5 or 0.8 therefore does **not** place this sequence inside the 3D spin-glass phase. A positive result would show only finite-width / tube structure, a negative result would not rule out 3D architecture, and an indeterminate result adds nothing. Recorded: **E2 NOT RUN — ACCESSIBLE EXACT-SPECTRAL GEOMETRY DOES NOT PRESERVE THE TARGET 3D THERMODYNAMIC LIMIT.** It is **not** recorded as "E2 NO GROWING PARTITION ARCHITECTURE FOUND", and **not** counted as a failed physical test |
| **C2S-04** | "Global Z₂ produces exact paired sectors" is withdrawn. At finite N the heat-bath Glauber chain is irreducible; global spin flip commutes with the generator and gives parity structure and configuration pairing, not two disconnected Markov sectors. A K = 2 structure explainable solely by this supplied symmetry is non-novel for C2 |

**C2 FINAL OWNER RULING:** C2 COMPLETE — GROWING BLIND RECOVERY PROVED CONDITIONALLY; ENDOGENOUS GROWING ARCHITECTURE
FROM A SIMPLE FIXED RULE NOT FOUND; INTRINSIC K-GROWTH OBSTRUCTION OPEN. C3 not opened.

**Applied to:** `C2_CHARTER.md`, `C2_GROWING_ARCHITECTURE.md`, `DA0_ZOOM_OUT_02.md`, `DA0_STATUS.md`. Created
`C2_FINAL_HANDOFF.md`.

## C3 CHARTER REPAIR 01 (owner review of `f266be2`; charter only, no computation)

**Owner ruling:** C3 CHARTER CONDITIONALLY ACCEPTED — PRIMARY MODEL ACCEPTED, REPAIR BEFORE RUN. The pre-repair boundary
is `f266be2f6432e4c8de6c73a183ea5135fe816a8e`.

| ID | correction |
|---|---|
| **C3R-01** | "n = 1 ⇒ no architecture" and "architecture needs n ≠ 1" are **withdrawn**. At T = T′ the σ-marginal is uniform, but J \| σ ~ N(s_b/μ, T/μ) keeps bond / loop correlations (plaquette ⟨ΠJ_b⟩ = 1/μ⁴). CPS mean-field has an ordered q > 0 phase for n ≤ 2, including n = 1. A1′ is renamed **EQUILIBRIUM ADAPTIVE-COUPLING COMPARATOR**; A0 / A1 remain the nulls. For T ≠ T′, the two reservoirs are supplied, and entropy production / heat flow is to be calculated at finite ε, not asserted |
| **C3R-02** | Primary analytic object: **𝓕_T(J) = (μ/2)ΣJ_b² − T log Z_T(J)**, with dJ/dt = −ε∇𝓕_T + slow noise. Stationary points depend on (T, μ), **not on T′**. Hessian μδ − (1/T)Cov_J(s_b, s_b′). Exact linear instability of J = 0 at **T = 1/μ** |
| **C3R-03** | **C3-F2:** exact local gauge symmetry σ_i → η_iσ_i, J_ij → η_iη_jJ_ij. Architecture is counted modulo all exact symmetries, using gauge-invariant observables. No gauge fixing may generate K-growth. Generated gauge-invariant frustration is distinguished from gauge copies |
| **C3R-04** | The pairwise median-PR criterion is withdrawn: it certifies independent product bits (PR ~ N/2). It is replaced by **PR_edge on elementary transitions** (independent bits → O(1)), plus a non-vanishing contrast **liminf S_ab > 0** |
| **C3R-05** | **C3-M1 resolved model-specifically.** Detector = symmetry-inequivalent stable minima of 𝓕_T joined by canonical minimum-barrier saddles, with diverging barriers / exit times. No clustering, PCCA, chosen K, lag-time partitions or fitted macrostate count. C3-A1 requires growth **and** collective edges **and** growing barriers |
| **C3R-06** | **Analytic preflight first** (landscape → symmetries → stability → loop expansion → bond-local vs collective), on separate approval. No large-lattice simulation before review |
| **C3R-07** | **CPS motivates the mechanism; it does not validate the local lattice theory.** The analysed CPS model is mean-field / infinite-range; the nearest-neighbour version is a new construction |

**C3-B stays closed. No C3 code, no C3 calculation.**

## C3 CHARTER REPAIR 02 (owner review of `38febe1`; before the analytic preflight)

**Owner ruling:** Repair 01 accepted; primary model remains accepted. **Analytic preflight APPROVED after Repair 02. No
simulation.**

| ID | correction |
|---|---|
| **C3R-08** | A plaquette (or any bounded loop) is the first gauge-invariant **inter-bond interaction**, not C3 collectivity. Its support is bounded, so its participation is O(1), and it does not satisfy C3-A1. The preflight question is whether bounded loop interactions bootstrap into extended, symmetry-inequivalent metastable structures |
| **C3R-09** | The primary dimension is **fixed at d = 3 before results**. d = 1, 2 are analytic controls only |
| **C3R-10** | Mandatory ferromagnetic-envelope test: Z_T(J) ≤ Z_T(\|J\|), so 𝓕_T(J) ≥ 𝓕_T(\|J\|), with strictness. It constrains global minima only. The key question is frustrated **local** minima with diverging barriers |

## C3 analytic preflight note (process; no repair)

- **Secondary sanity check (`c3/c3_preflight_check.py`).** The first run reported a "frustrated" stable minimum on the
  cube at T = 1.1. Inspection showed it was the J = 0 minimum with round-off signs (|J_b| ~ 1e-10). A declared
  computational zero (|J_b| ≤ 1e-6 means absent) was added, and the logged run is the corrected one.
- The analytic results stand independently of this check.

## C3 PREFLIGHT SCOPE REPAIR 01 and PRIMARY MODEL 1 CLOSURE (owner review of `dde3b94`)

**Owner ruling:** **option (a)**. The window simulation is **not** authorised. The preflight is accepted after the scope
repairs. Primary Model 1 is closed as the primary C3 route, with its residuals preserved. This is **not PF-D**, and not a
universal no-go.

| ID | correction |
|---|---|
| **C3P-01** | "PF-B for T < T_lin/2" is narrowed to **PF-B — PROVED FOR SIGN / FRUSTRATION ARCHITECTURE ONLY** (full support, unfrustrated plaquettes, ≤ 8 winding sectors modulo gauge). For the full C3-A question, **NONUNIFORM MAGNITUDE ARCHITECTURE REMAINS OPEN**: Griffiths + Tarski give a uniform greatest fixed point but do not exclude intermediate nonuniform ones, and uniqueness is only claimed below an uncomputed T₀ |
| **C3P-02** | "The adaptive transition is first-order-like at T₁" is narrowed to: **the uniform stationary branch has a subcritical / first-order-like bifurcation structure. A full-space thermodynamic / adaptive transition and its location are not established.** Any uniform-ray free-energy crossing is a uniform-ray statement only |

**PRIMARY MODEL 1 TERMINAL:** PRIMARY MODEL 1 CLOSED — ADAPTIVE COUPLING GENERATES ORDER BUT NO GROWING COLLECTIVE
ARCHITECTURE MECHANISM WAS FOUND; SIGN/FRUSTRATION ARCHITECTURE IS O(1) IN THE PROVED LOW-T REGIME; FULL ARCHITECTURE
RETAINS INDETERMINATE RESIDUALS.

**Construction lesson.** Pure local Hebbian reinforcement is biased toward unfrustrated (Mattis-like) order. The next
law needs **endogenous local competition**. PM1 is retained as a C3 comparator.

**Primary Model 2:** a candidate audit / proposal only (`C3_PM2_CANDIDATE_AUDIT.md`). It is **not run** and **not
selected**. It adds a one-line structural result: **sign-blind constraints preserve the unfrustrated global envelope**.

**Boundaries preserved:** Repair 02 `31b5714e134b78dfe6b6ff75c814654c94f09cb2`; C2 `99428ff`; RA0-frozen `ab4fd86`.

## PM2 CHARTER REPAIR 01 (owner review of `80677c0`; before any PM2 computation)

| ID | correction |
|---|---|
| **PM2R-01** | Drive h_i = +1 (fixed source density), h_outlet = −(N − 1). The old normalisation made m_reroute ≤ 1, so ρ > 0 was impossible. PROP PM2-S (exact): h → ah maps stationary C → a^{2/(γ+1)}C and E → a^{2γ/(γ+1)}E, with topology unchanged. Trajectories are covariant only if the initial condition is rescaled. Reduced Ē = E/(N − 1)^{2/3} |
| **PM2R-02** | Pruning (C = 0) is absorbing, so fundamental swaps are **static** adjacencies and barriers are **static** diagnostics, not Stage-A dynamics. Lower-semicontinuous boundary convention. Stage-A ceiling: **PM2-A PARTIAL … DYNAMICAL BARRIERS / REVERSIBLE TRANSITIONS UNESTABLISHED** (not C3-A1) |
| **PM2R-03** | Seed-specific topologies and the number of distinct trees are not derived architecture. Only ensemble-stable structural laws can be earned. Q3 is replaced |
| **PM2R-04** | P1: per-tree CCDF fit on the window [√(N/10), √(10N)]; one τ per seed; t-intervals for PM2 and R-TREE; a single deterministic τ_SP; the pass rule is restated |
| **PM2R-05** | P2: L_branch(e) = longest upstream path within B_e to e's downstream endpoint; same rules |
| **PM2R-06** | P3: the full canonical swap set (all non-tree f × cycle edges e). n_reroute(f, e) = S_e exactly. Median per tree, fit ~N^ρ. Structural only. Static-barrier diagnostic on 10 preregistered random swaps per tree |
| **PM2R-07** | A1 is not an architecture-free null; A0 is the true null |

## PM2 STAGE-A IMPLEMENTATION DECLARATIONS (numerical only; written and committed BEFORE the Stage-A grid was run)

**Scope.** These are implementation choices. γ, ν, κ, δ, geometry, drive, sizes, seeds and statistics are unchanged.
Before this entry, only timing and validation runs had been made, at L ∈ {8, 12, 16, 24, 32, 48}, seeds 0 – 7. The only
quantities inspected were final topologies, energies, step counts and wall time. **No Stage-A statistic (τ, η_H, P3,
Strahler, barriers, controls) was computed or inspected.**

| ID | declaration |
|---|---|
| **PM2-I0** | **Abandoned exploration (no result used).** A semi-implicit (per-edge backward-Euler) integrator in C was tried at L = 8. It hit the 400 000-step cap with residual 7×10⁻², and its topology disagreed with explicit Euler. It was removed from the code |
| **PM2-I1** | **Exact spanning-tree early stop.** Stop as soon as the active set is a connected set of N − 1 edges, and set C to its exact limit C*_e = (S_e²/νγ)^{1/(γ+1)}. **Proof:** on a spanning-tree support, Kirchhoff fixes Q_e = ±S_e ≠ 0 independently of C. Each active edge then obeys an autonomous 1-D ODE with a unique attracting root C* > 0 (dC/dt > 0 as C → 0⁺), so no active edge can be pruned; pruned edges never revive (PM2R-02). The final topology is therefore already fixed. **Validation** (`c3/pm2_earlystop_check.py`): against the registered full-Euler run, L = 8 seeds 0 – 7 and L = 12 seeds 0 – 3 give **12/12 the same topology**, max\|ΔC\|/max C ≤ 3×10⁻⁹, \|ΔE\|/E ≤ 8×10⁻¹⁵ |
| **PM2-I2** | **Production integrator: explicit Euler in y = √C** (exact rewrite for γ = 1/2, using Q_e = C_e Δp_e): dy_e/dt = (κ/2)(y_e Δp_e² − νγ). This is the **same ODE**. In C, a dying edge has relative rate ∝ C^{−1/2}, which throttles the global step; in y, extinction is linear in time, so the variable is non-stiff. Step: max_e \|Δy_e\|/max(y_e, y_ref·max y) ≤ max_rel_change, with y_ref = 10⁻³. An edge whose y crosses 0 is pruned (finite-time extinction of the exact ODE), and the 10⁻¹² floor is kept. Stop rules are as registered, plus PM2-I1. **Production max_rel_change = 0.0125**, replacing the registered C-Euler 0.05 |
| **PM2-I3** | **Why PM2-I2 (step-refinement check, `c3/pm2_dt_check.py`).** Final topologies were compared against y-Euler at 0.00625 (the finest reference), L = 8 seeds 0 – 5 and L = 12 seeds 0 – 3. **Perturbed seeds:** registered C-Euler at 0.05 differs by 2 – 10 edges in 9/9 runs, and still differs at 0.0125 in 3/9. y-Euler at 0.05 differs in 3/9 (2 – 4 edges), and **at 0.0125 in 0/9**. The registered C-Euler at 5% is therefore **not step-converged at the level of individual trees**; y-Euler at 0.0125 is the most converged affordable choice. Cross-integrator checks at L ≤ 24 (`c3/pm2_yvar_check.py`) showed per-run topology differences of 2 – 46 edges, with energies agreeing to ≲ 1 % |
| **PM2-I4** | **The per-run tree is discretisation-sensitive; Stage A is graded on ensembles only.** Consistent with PM2R-03, seed-specific trees are not used. **Uniform start (seed 0):** its topology does **not** converge under step refinement (differences of 4 – 30 edges at every step size). The symmetric start breaks symmetry through round-off, so **the uniform-start tree is reported as round-off-selected and is not used for any criterion.** **Preregistered robustness requirement (Stage-A positive item 5):** the full grid is also run with y-Euler at max_rel_change = 0.05. Every pass/fail verdict (P1, P2, ρ > 0, Strahler growth) must be **identical at both step sizes**, or that criterion is graded *not robust* |
| **PM2-I5** | **Parallelism.** The 21 PM2 runs per L execute in a 4-process pool. Results are deterministic and independent of scheduling |

## PM2 STAGE A — RUN RECORD (after the grid; additive)

| ID | record |
|---|---|
| **PM2-S1** | **Grid run exactly as preregistered.** L ∈ {16, 24, 32, 48, 64}; uniform start plus seeds 1 – 20 (δ = 0.01); R-TREE seeds 1 – 20; SP-TREE; A0 / A1 / A2. No change to γ, ν, κ, δ, geometry or drive. No noise. Production integrator per PM2-I2. The robustness pass per PM2-I4. Step cap never reached; 105 / 105 production runs ended on spanning trees |
| **PM2-S2** | **Operational restart (no effect on results).** The first production launch was killed after L = 24 (about 14 min) because it would have exceeded the 2 h background-job limit. Old pool workers were terminated. The grid was relaunched from scratch as a detached process, and the logged run is the complete relaunch. The runs are deterministic, so the discarded partial output was identical in content (the L = 16 summary of the first launch, inspected before the kill, matches the relaunch exactly) |
| **PM2-S3** | **Declared reading of Stage-A positive item 4** (no numeric rule in the charter), fixed in `c3/pm2_stageA_summary.py`. It was written after the L = 16 numbers were visible and before L ≥ 24 were. The PM2 Strahler slope interval must lie entirely above the R-TREE slope interval and above the SP slope, and at L = 48 and 64 the PM2 interval must lie above the R-TREE interval and the SP value. The outcome (FAIL) holds under any reading, since PM2 never exceeds R-TREE. P3 reading: the 95% lower bound of pooled-OLS ρ must be > 0 |
| **PM2-S4** | **Exact readings added for interpretation** (INTERNALLY PROVED / NOT EXTERNALLY REVIEWED). **PROP PM2-U:** for even L, the exact uniform-start trajectory cannot reach a spanning tree. The diagonal reflection fixes no edge, so σ-symmetric active sets have even size, while N − 1 is odd. Uniform-start trees are therefore round-off-selected and used for nothing. **Corollary of n_reroute = S_e:** ρ > 0 is generic to every growing tree family; the controls give ρ_R = 0.78, ρ_A2 = 0.75 and ρ_SP = 0.51 |
| **PM2-S5** | **Stage-A outcome** (`C3_PM2_STAGE_A.md`). Production verdicts: P1 PASS (finite-size; τ_PM2 drifts 0.99 → 0.60 toward the R-TREE band ≈ 0.38 – 0.44), P2 FAIL, P3 PASS (ρ = 0.674 [0.651, 0.698], generic and below the controls), Strahler-beyond-controls FAIL. **The Stage-A positive is not met; PM2-A PARTIAL is not reached.** Proposed terminal: **PM2-B — GENERATED TREE / NETWORK — GENERIC TOPOLOGICAL HIERARCHY ONLY** (owner to rule). Robustness pass (step 0.05): all verdicts identical (item 5 met) |

## PM2 STAGE A — OWNER RULING (review of `2139cc4`)

**Ruling:** **PM2 DETERMINISTIC STAGE A ACCEPTED — PM2-B TERMINAL ACCEPTED WITH SCOPE DISCIPLINE.**

**Terminal (precise scope; supersedes the shorter PM2-B wording in `C3_PM2_STAGE_A.md` §4 and in PM2-S5):**

> **PM2-B — GENERATED TREE / NETWORK; REGISTERED HIERARCHY IS GENERIC / RANDOM-TREE-LIKE; NO DISTINCT GROWING COLLECTIVE
> ARCHITECTURE ESTABLISHED.**

**Deterministic PM2 Stage A is CLOSED.** The reviewed Stage-A scientific boundary is
**`grut-directed-autonomy-0 @ 2139cc436e79782c231286f4437bdf62402140e0`**. This is recorded as a boundary on the branch;
no new branch or tag is created.

| ID | correction / record |
|---|---|
| **PM2-O1** | **Owner reading accepted:** trees are generated from the perturbed ensemble; flow feedback matters (A1 selects no tree); P1 passes its finite-size rule but drifts toward R-TREE; P2 fails; P3 ρ > 0 is generic, with PM2 below R-TREE / A2; Strahler does not exceed R-TREE; reduced static barriers do not grow; the robustness pass preserves every verdict; PM2-A PARTIAL is not reached; **TRUE COMPRESSION = 0** |
| **PM2-O2** | **Scope firewall.** PM2-B must **not** be read as "adaptive flow networks have no hierarchy" or as any universal negative. The supplied hierarchy A2 is itself not cleanly separated from R-TREE by the registered P1 / P2 / Strahler suite. **Earned negative:** no distinctive growing collective architecture was established **under the registered invariant suite** |
| **PM2-O3** | **R1 (asymptotic τ) is OPEN and NOT ACTIONABLE at present.** No larger-L run will be made to chase τ. P1 is the lone survivor, and its separation shrinks with size. A larger-L distinction would not repair the failures of depth, rerouting distinctiveness or barriers |
| **PM2-O4** | **PM2-BIS is NOT opened.** Fluctuating sources and edge revival remain physically plausible future mechanisms. Adding them before the architecture detector is qualified would make any later positive uninterpretable |
| **PM2-O5** | **Next gate: C3-D0 — ARCHITECTURE DETECTOR QUALIFICATION** (`C3_D0_ARCHITECTURE_DETECTOR.md`). Charter and candidate audit only, calibrated on controls only. PM2 Stage-A outputs are frozen and may not be used to tune a detector. No D0 computation is made until the detector choice and the second positive hierarchy are preregistered and reviewed |

**Preserved:** PM1 terminal; all previous boundaries; RA0-frozen `ab4fd86`. **C3-B closed. No PR. No merge.**

## D0 CHARTER REPAIR 01 (owner review of `e52a276`; before any D0 computation)

**Ruling:** C3-D0 CHARTER CONDITIONALLY ACCEPTED. After this repair is committed, the D0 control computation is approved
to run without another owner stop. The repaired text is `C3_D0_ARCHITECTURE_DETECTOR.md` §R.

| ID | correction |
|---|---|
| **D0-O1** | P1 (nested 3 × 3 block tree) accepted as the second supplied hierarchy; ARCHITECTURE SUPPLIED; never a C3 result |
| **D0-O2** | Φ accepted as the primary. The √N / √S exponents are frozen. **Scope:** an operational detector of pathwise recurrence of collective coarse branching motifs; a failure rejects this operational definition only |
| **D0-O3** | **Wording error corrected.** Randomising below N^{1/3} is **not** a vanishing fraction of the log-mass range; it is about **one-third**. The ε-controls keep about **two-thirds** of the supplied recursive range (a growing majority). Q6 and the ε-controls are kept as a substantial robustness test |
| **D0-O4** | Q3 replaced by an adjacent-pair non-shrinking rule: Δ_{j+1} ≥ Δ_j − [h_R(L_j) + h_R(L_{j+1})]. Q6 is made exact: Δlow_j = lo_{Fε} − up_R > 0 for all j; lo_{Fε} > Φ_{N1}; adjacent non-shrinking with the combined Fε and R half-widths. Interpretation: **scale-structural, non-shrinking over the preregistered finite grid**, not asymptotically proved |
| **D0-O5** | Grids fixed exactly: G2 = {16, 32, 64, 128}, G3 = {27, 81, 243} |
| **D0-O6** | Implementation firewall: detector plus unit tests (relabelling, children-order and isomorphism invariance; analytic N1 Φ = 0; path Φ = 0; complete binary depth 6 Φ = 6/7; construction validity) **before** any control value. No inspection of control Φ while modifying code, except for exact-test violations, which must be ledgered |
| **D0-O7** | Grades: **D0-A-N** (noise-tolerant; the only grade that may support applying Φ to a stochastic generator) / **D0-A-T** (template-only; an instrument result, not a failure) / D0-B / D0-INDETERMINATE |
| **D0-O8** | Run only N0, N1, P0, P1, P0-ε, P1-ε on the fixed grids and seeds. PM2 trees are not scored |

**Preserved:** PM2 Stage-A boundary `2139cc4`; PM2-B and its scope firewalls; all C2 / C3 boundaries; RA0-frozen
`ab4fd86`. No PM2-BIS, no C3-B, no PR, no merge.
