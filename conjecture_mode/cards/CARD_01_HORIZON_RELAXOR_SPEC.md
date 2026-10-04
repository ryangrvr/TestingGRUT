# CARD #1 — HORIZON RELAXOR — C0 SPECIFICATION (v1)

> **POSTULATED HEURISTIC TOY LAW — NOT DERIVED FROM THE GRUT EFFECTIVE STRESS TENSOR.**
> The frozen record itself labels this form "heuristic / illustrative of STRUCTURE; the exact form needs the full Calzetta–Hu in-in stress tensor" (`calc/wz_dark_energy.py`, docstring). Nothing in this card upgrades that status.

| | |
|---|---|
| **Branch** | `grut-conjecture-mode-1`, from `grut-program-governance-1-frozen @ 5baa8afd83b8a2dbce31b8f3c4bd484417afbc11` |
| **Governing law** | Conjecture Mode Charter 01, including Annex B (Card #1 requirements) and PROGRAM_CAMPAIGN_GATE_01 (Door C) |
| **Commit class** | **CARD-01 PRE-DATA SPEC BOUNDARY**. Pre-data only |
| **State** | **NOT EXECUTABLE.** Owner threshold ruling pending (`CARD_01_OWNER_THRESHOLD_RULING.md` does not exist). No network preflight run. No DESI product touched |
| **Executor C1 verdict** | **CARD-01-C1-PASS (with notes N1–N6)**, subject to owner review (§C1) |

---

## C0.1 — Residual target

**Instantiated entry.** The **inserted second, slow cosmological relaxation scale and its response structure**: the "two-scale IR mode for w(z)". The program records it as an **assumed / supplied input**, not a derived one.

| Frozen source | Location | Content |
|---|---|---|
| scout-1 constitution | `scout-1 @ a2987fed5880da62c874857faa0a93099d23f64d`, `CHARTER.md` line 53 ("assumed (name each at the rung it enters)") | "… the **two-scale IR mode for w(z)** …" listed among the **assumed** inputs |
| Rung-7 record | same commit, `calc/RESULTS_wz.md` §(A)–(B), §Status | Observable w(z) evolution "requires a second, slow scale τ₂ ∼ 1/H₀". The two-scale vacuum is "a new named structural input". "Ledger +2 (ε + the two-scale commitment)" |
| Emergence chain | same commit, `EMERGENCE_CHAIN.md` §10 | "+3: the amplitude, the two-scale commitment, the single-departure shape — all priced" |
| Booked ansatz | same commit, `GRUT_PREDICTION_GATE_GAMMA_T.md` horn (b) | χ(ω) = A/(1 − iωτ_c) + B/(1 − iωτ₂), τ₂ ∼ 1/H₀ — "the BOOKED two-scale ansatz" |
| Frozen residual ledger (governance base) | `residual_synthesis/MASTER_RESIDUAL_LEDGER.md` row **physical scale** (line 52): "scale-priced (τ, ε …)", final grade **SUPPLIED**; `residual_synthesis/SUPPLIED_INPUT_MASTER_TABLE.md` row **SP-4** ("objectives, thresholds, scales (… τ, ε)"), **SUPPLIED** | The nearest frozen-ledger rows: a supplied scale τ and a supplied amplitude ε |

**Granularity note (N6).** The governance-base master residual ledger has **no dedicated row** for the cosmological IR relaxation mode. It is carried only inside the generic **physical scale / SP-4** rows. The precise named input is the scout-1 "two-scale IR mode for w(z)". This card instantiates that supplied input. It **does not derive it** and does not claim to eliminate it.

## C0.2 — Origin / GRUT rationale

**Admissible origin categories used:** (b) a declared fork, and (d) the core responsiveness / memory principle. No other category is claimed.

1. **Finite-memory / responsive vacuum (d).** GRUT's vacuum has a causal frequency-dependent response χ(ω). Out of equilibrium, its effective equation of state can leave w = −1 (`calc/RESULTS_wz.md`, Goal; W1 in `calc/RESULTS_wz_sign.md`).
2. **Rung-7 requirement of a slow IR scale.** The confirmed UV-cutoff memory (H₀τ_c ∼ 10⁻⁴⁰) gives w = −1 flat to about 80 decimals. Any observable w(z) evolution requires response power at ω ∼ H(z), i.e. a second slow scale (`RESULTS_wz.md` §(A)).
3. **Declared two-scale fork (b).** Single-pole tabletop behaviour (rung 3) and an evolving w(z) coexist only by scale separation: a UV cutoff plus an IR horizon-scale mode. The record names this as an explicit commitment, a "named, falsifiable structural input" (`RESULTS_wz.md` §(B); `README.md` lines 139–144).
4. **Horizon motivation for τ₂ ∼ H₀⁻¹.** The IR scale is "horizon-motivated (de Sitter / Gibbons–Hawking) … possibly natural rather than tuned", with the cosmic-coincidence question explicitly left open (`RESULTS_wz.md` §(B)).
5. **Branch (dissipative).** The second-law (bulk-viscosity, ζ ≥ 0) reading fixes the dissipative branch on the **side** w ≤ −1. It does **not** fix the slope w_a (`calc/RESULTS_wz_sign.md`, "The amended boundary").

**What remains arbitrary / POSTULATED (not dictated by the origin):**
- the identification τ₂ = H₀⁻¹ **exactly** (the origin gives only τ₂ ∼ H₀⁻¹);
- the single-relaxor storage response shape R(x) = x²/(1+x²);
- the use of a **fixed reference background** E_ref (0.31, 0.69) inside R;
- the amplitude ε;
- the choice of the dissipative reading over the record's reactive / phase-lag reading. The latter sits on w ≥ −1 and is the ε > 0 control. The branch is "fixed by origin" only through the record's second-law **interpretation**.

**C-F9 origin audit — PASS.**
- **The DESI anomaly is not the origin.** The toy was written to test structural differentiation (w(z) evolution being impossible for a static-Λ family). Its own §(D) reports a **mismatch** with the DESI pattern rather than fitting it.
- **The branch is not chosen to approach the data.** The dissipative branch gives `w0 < −1` under the registered projection (§C0.5). That is the side opposite the published headline preference for `w0 > −1` (disclosed in §PED). The GRUT-justified branch is tested regardless (C-F9).
- **Disclosure.** The historical record repeatedly discussed DESI's direction, including a retracted "wrong sign" reading and an overseer self-correction of a "toward-DESI" over-claim (`RESULTS_wz_sign.md`). These are recorded under §PED. None of them altered R(x), τ₂ or E_ref.

## C0.3 — Exact POSTULATED law (transcribed, not reinvented)

**Historical sources** (all at `scout-1 @ a2987fed5880da62c874857faa0a93099d23f64d`, read via `git show`; the scout-1 checkout was not touched):

| File | Blob |
|---|---|
| `calc/wz_dark_energy.py` | `cd900011717730d68c0fbe39bd1386e05e67fa04` |
| `calc/RESULTS_wz.md` | `d9d2b317dee742351005455537c115e09158a90c` |
| `calc/RESULTS_wz_sign.md` | `70713240c30977069608a125de2194006aed9a2b` |
| `PHYSICS_LEDGER/RUNG7_TWO_POLE_COMPARISON.md` | `f7ccd77931712856f042ee7646874b80f9814b79` |

**Historical implementation** (`calc/wz_dark_energy.py`):
- `OM = 0.31`, `OL = 1 − OM`
- `E(z) = sqrt(OM (1+z)^3 + OL)`
- `x = E(z) * Htau`
- `w = −1 + eps * x²/(1 + x²)`
- with `Htau = 1` for the IR horizon scale

**Card #1 v1 law:**

```
E_ref(z)² = 0.31 (1+z)³ + 0.69          (fixed reference function; NOT the varied cosmology)
x(z)      = H_ref(z) τ₂ = E_ref(z) · (H₀ τ₂),   τ₂ = H₀⁻¹  ⇒  x = E_ref(z)
R(z)      = E_ref(z)² / [1 + E_ref(z)²]
w(z; ε)   = −1 + ε R(z)
```

**Model semantics (owner decision, frozen).** The response is a **prescribed one-parameter shape `w(z; ε)`**, not an implicit equation for the self-consistent H.
- The reference function E_ref, which enters only the postulated shape, is fixed at (0.31, 0.69).
- The physical parameters (ω_b, ω_c, H₀ or θ, A_s, n_s, τ_reio, nuisances) are varied by the likelihood code and **do not** feed back into R.
- ρ_DE follows from separate covariant conservation: ρ_DE(z)/ρ_DE,0 = exp[3ε ∫₀ᶻ R(z′)/(1+z′) dz′]. Its present value is fixed by flatness closure.
- Replacing E_ref by the self-consistent H would be **a new card or version**, priced separately.

**Branches.**
- **Primary (preregistered GRUT claim):** ε < 0, with ε = 0 as the nested ΛCDM null.
- **ε > 0: CONTROL — NOT A GRUT CLAIM.** It cannot affect the Card #1 verdict.

**Domain.** z ∈ (−1, ∞) for the background. Data contact uses 0 ≤ z ≲ 10⁹ (early universe through today). The registered CPL window is 0 ≤ z ≤ 2.

### Analytic sign proof (Annex B §2.2)

For z > −1, E_ref² = 0.31(1+z)³ + 0.69 > 0. Hence R = E_ref²/(1+E_ref²) ∈ (0, 1).

R is strictly increasing in z: dR/dz = (dE_ref²/dz)/(1+E_ref²)² = 0.93(1+z)²/(1+E_ref²)² > 0.

For z ≥ 0, E_ref ≥ 1, so R ∈ [½, 1). In particular:
- R(0) = ½;
- R(2) = 0.900596…;
- R → 1 as z → ∞;
- R → 0.69/1.69 = 0.40828… as z → −1.

Therefore, for **ε < 0**: w(z) − (−1) = εR(z) < 0 for **every** z > −1. On z ≥ 0 this gives the tighter bound

  −1 + ε < w(z) ≤ −1 + ε/2 < −1.

**ε < 0 ⇒ w(z) ≤ −1 (strictly < −1) on the whole declared domain.** The sign proof **passes**, so CARD-01-REPAIR is not triggered.

The dissipative branch therefore never touches w = −1 (no phantom-divide crossing). This is consistent with the record's no-crossing structure (`RESULTS_wz_sign.md` (i); `RUNG7_TWO_POLE_COMPARISON.md`).

## C0.4 — Information price

| Item | Status | Count |
|---|---|---|
| τ₂ = H₀⁻¹ (exact identification; the origin gives only τ₂ ∼ H₀⁻¹) | **POSTULATED relation** | 1 dimensional-scale relation |
| ε | **one fitted scalar amplitude** | 1 parameter |
| R(x) = x²/(1+x²) (single-relaxor storage response) | **POSTULATED functional form** | 1 function |
| dissipative sign (ε < 0) | **fixed by origin** (record's second-law interpretation); not derived by this card. If that interpretation is withdrawn, it becomes one priced branch bit | 0 (conditional 1) |
| reference background (Ω_m^ref, Ω_Λ^ref) = (0.31, 0.69) | **projection / input choice inherited** from the historical toy | 2 constants (one if flatness is imposed) |
| perturbation closure for the effective fluid: CAMB **PPF** (Fang–Hu–Lewis), default settings | **projection convention, physically load-bearing** for CMB / ISW at the perturbation level | 1 convention |
| w(a) table: 3000 log-spaced points, a ∈ [10⁻⁹, 1], CAMB internal interpolation | numerical tolerance (verified: background H agrees with independent quadrature to ≤ 1.3 × 10⁻⁹) | 0 |
| D3 standard parameters (ω_b, ω_c, θ or H₀, A_s, n_s, τ_reio, Σm_ν fixed at 0.06 eV) and likelihood nuisances | standard ΛCDM-sector parameters, shared with the null and the CPL comparator | as in the null |

**Card #1 is not parameter-free.** It is a **one-parameter, one-sided shape family** (ε ≤ 0). Relative to ΛCDM it adds:
- 1 fitted parameter (ε);
- 2 postulated structures (τ₂ = H₀⁻¹ and R);
- 2 inherited reference constants;
- 1 perturbation closure.

Any low-dimensional shape relation (§C0.5) is **conditional on the frozen projection**.

## C0.5 — Projections

### (a) Registered CPL diagnostic (diagnostic only; never the kill test)

**Frozen:**
- the window `0 ≤ z ≤ 2`;
- `w_CPL(z) = w0 + wa z/(1+z)`;
- the single **unweighted global least-squares** projection in redshift, `J(w0, wa) = ∫₀² [w_model(z) − w_CPL(z)]² dz`.

No alternative weighting is permitted in v1. The local z = 0 derivative map is **not** the registered projection.

**Model-only result** (`code/card01_cpl_projection.py`, log `code/card01_cpl_projection.log`; mpmath at 30 digits, scipy cross-check agrees to 10⁻¹²). Because w_model + 1 = εR and J is quadratic:

```
(1 + w0, wa) = ε · (c0, ca),   c0 = 0.433262795583067,   ca = 0.674115607186169
wa / (1 + w0) = ca / c0 = 1.55590467046443   — EXACT ε-independent constant; the locus is EXACTLY linear
```

Fractional residual of the fit: ∫(R − fit)² / ∫(R − R̄)² = 0.0225.

| Question (Annex B §2.5) | Model-only answer |
|---|---|
| Is the locus exactly linear? | **Yes.** It is a straight line through ΛCDM (−1, 0) with slope 1.5559 |
| Frozen global-fit ratio | **wa/(1+w0) = 1.55590…** (not 1.37) |
| Difference from the old local-derivative map | The historical value-plus-slope map gives c0 = R(0) = ½ and ca = R′(0) = 3Ω_m^ref/4 = 0.2325, so the **ratio is 0.465** (= 3Ω_m^ref/2). This is labeled **UNREGISTERED HISTORICAL DIAGNOSTIC**. The registered global ratio is 3.35× larger |
| Provenance of "1.37" | **Not produced by either map, and not found in the frozen scout-1 record** (searched at a2987fe). Its origin is the earlier reconnaissance outside this record. It stays unbanked |
| Can the global CPL fit give w0 > −1 although pointwise w ≤ −1? | **No, for this law and projection.** 1 + w0 = ε c0 with c0 = 0.4333 > 0, so **w0 < −1 for every ε < 0** |
| Sign of wa on the primary branch | wa = ε ca with ca > 0, so **wa < 0 for every ε < 0** |
| Example | At ε = −0.1: (w0, wa) = (−1.04333, −0.06741), with pointwise w ∈ [−1.0901, −1.05] on [0, 2] |

**Structural fact (model-derived, not a forecast).** The primary-branch CPL image is the open ray {w0 < −1, wa < 0} along the line wa = 1.5559(1 + w0). The ε > 0 control image is the opposite ray {w0 > −1, wa > 0}. Neither ray enters the quadrant {w0 > −1, wa < 0}. As in Annex B §4, this is a diagnostic geometry statement. **The verdict is decided only by D3.**

### (b) D3 primary model interface (decisive route) — frozen implementation

**Stack:**
- **CAMB 2.0.4** (pip), Boltzmann and background solver;
- **Cobaya 3.6.2** (pip), the likelihood framework. This is the same framework family DESI DR2 cosmology analyses used, per general knowledge and to be confirmed (§C0.6).

**Interface:** `code/card01_cobaya_theory.py`, class `Card01CAMB`.
- It subclasses Cobaya's CAMB theory.
- It consumes one sampled parameter, `card01_eps`.
- It passes the **exact** w(a) = −1 + εR(z) to `DarkEnergyPPF.set_w_a_table` on a table of 3000 log-spaced points over a ∈ [10⁻⁹, 1].
- At ε = 0 it uses CAMB's cosmological constant.
- **The w(z) shape is not modified for software convenience and is never replaced by its CPL projection.**

**Model-only verification** (no likelihood, no data):

| Check | Result | Log |
|---|---|---|
| **F1.** CAMB background H(z) vs an independent quadrature of ρ_DE at the same physical densities, z ∈ {0 … 1100}, ε ∈ {−0.05, −0.2, −0.5, +0.2} | max relative difference ≤ **1.3 × 10⁻⁹** | `code/card01_d3_feasibility.log` |
| **F2.** Continuity, ε = −10⁻⁶ vs ΛCDM | max \|ΔH/H\| = 1.4 × 10⁻⁷ | `code/card01_d3_feasibility.log` |
| **F3.** Full lensed TT/EE spectra compute on both branches | computed without error, ℓ_max = 2500 | `code/card01_d3_feasibility.log` |
| **F4.** Cobaya plumbing: sampled `card01_eps` reaches CAMB | Cobaya-returned H(z) and D_A(z) at z = 0.51, 1.32, 2.33 equal direct CAMB to **0.0** relative difference; r_drag and C_ℓ returned | `code/card01_cobaya_plumbing_test.log` |
| Cost | one full CAMB C_ℓ evaluation (ℓ_max 2500, lensing accuracy 1) ≈ **1.7 s** on this 4-core container | — |

**D3 procedure (frozen; numerical thresholds are owner-lock slots, §S):**

1. **Likelihood combination:** the owner-approved S6 primary. Each robustness combination is run separately.
2. **Primary-branch profile likelihood.** For each ε on the frozen grid
   ε ∈ {0, −0.005, −0.01, −0.015, −0.02, −0.03, −0.04, −0.05, −0.06, −0.08, −0.10, −0.12, −0.15, −0.20, −0.25, −0.30, −0.40, −0.50, −0.70, −1.00} ∩ [ε_min, 0],
   compute χ²_min(ε) by Cobaya's minimizer (default backend). All non-ε parameters are free under the standard priors of the S6 likelihood set.
   - **Minimizer protocol:** 3 independent starts per grid point. Accept the lowest χ². Flag the point if the starts disagree by more than 0.2 in χ².
   - After the grid, run a free-ε minimization on [ε_min, 0], seeded from the best grid point.
3. **Comparators under the identical likelihood (C-F7):**
   - **ΛCDM** = ε = 0, the nested null, from step 2;
   - **w0waCDM (CPL)** χ²_min, using Cobaya CAMB's built-in `w`, `wa` with PPF. This is needed only to define the "targeted signal" magnitude (S3).
4. **Control (labeled):** the same grid mirrored to ε > 0, reported as `CONTROL — NOT A GRUT CLAIM`. It cannot affect the result state.
5. **Posterior (secondary report only):** an MCMC on ε ∈ [ε_min, 0] with Cobaya's MCMC. Stop at Gelman–Rubin R−1 < 0.02. This is reported, not decisive.
6. **D1/D2 diagnostics:** official DESI DR2 w0wa chains, **only if reachable at preflight**. If not, D1 is reported as unavailable and **no substitute** is used (§4B). D2 overlays the frozen CPL line from (a). All are labeled `DIAGNOSTIC — NOT THE KILL TEST`.

If the full CMB likelihood cannot be evaluated (blocked data), D3 is **not** performed with a compressed-CMB substitute without a prior owner downgrade ruling (Charter §4B). The card is then graded CARD-01-ACCESS-BLOCKED, or CARD-01-C2R by ruling.

## C0.6 — Dataset proposal (pre-data; nothing opened or downloaded)

**Documentation access.** `www.desi.lbl.gov` (the official DESI site, documentation page) is **blocked by this container's egress proxy**; a WebFetch of the DR2 chains release page returned EGRESS_BLOCKED. That was a documentation request, not product access. The product reachability preflight was **not** run.

The proposal below therefore rests on:
- (i) the owner's statement that the official site confirms DR2 cosmology chains and posterior-maximization products are released;
- (ii) the **likelihood definitions shipped inside the Cobaya 3.6.2 pip package**. Only module names and the data-file **names** in YAML definitions were read; no data values, no data installation;
- (iii) the executor's general knowledge of the published analysis set-up, flagged **TO CONFIRM against official documentation at preflight**.

**Proposed S6 hierarchy (CLAUDE-PROPOSAL, non-binding; see §S6):**

| Role | Combination | Cobaya modules (definitions present locally; data not installed) |
|---|---|---|
| **Primary** | **DESI DR2 BAO + CMB + DES-Y5 SN** | `bao.desi_dr2.desi_bao_all` (data files `bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt`, from the Cobaya `bao_data` GitHub release v2.6); CMB = Planck PR4 CamSpec high-ℓ TTTEEE (`planck_NPIPE_highl_CamSpec`) + Planck 2018 low-ℓ TT (`planck_2018_lowl.TT`) + low-ℓ EE (`planck_2018_lowl.EE`) + CMB lensing (Planck PR4 + ACT DR6; external `act_dr6_lenslike` package) **[composition TO CONFIRM]**; SN = `sn.desy5` |
| Robustness 1 | DESI DR2 BAO + CMB + Pantheon+ | `sn.pantheonplus` |
| Robustness 2 | DESI DR2 BAO + CMB + Union3 | `sn.union3` |
| Not primary (owner decision) | DES-Y5 recalibrated ("Dovekie") | `sn.desdovekie` exists in Cobaya 3.6.2. It is offered as a robustness row only if the owner wants it. The principal published DR2 comparison used the original DES-Y5 (TO CONFIRM) |

**No DR1 substitution.** `bao.desi_2024_*` (DR1) modules are present in the package and are **excluded**.

**Provenance question for the owner (O3).** The CMB and SN likelihoods are official products of Planck, ACT and DES, not DESI. The DESI BAO Gaussian likelihood files are distributed through the Cobaya `bao_data` repository. The owner must rule whether these distributions count as "official", conditional on SHA-256 matching against each collaboration's release at preflight.

## C0.7 — Kill form (governance nested-null states; no forecast)

The result states are exactly those of Charter §4C and Annex B §5:
- EXPLANATORY-SURVIVES
- EXPLANATORY-KILLED
- MODEL-KILLED
- INCONCLUSIVE
- ACCESS-BLOCKED
- RECONNAISSANCE
- PROCEDURE-VOID

Their numerical boundaries are the **owner-lock slots S1–S8**. They are currently **EMPTY**; executor proposals are in §S.

**Operational decision tree (proposed with S1–S5; binding only once owner-approved).** Let:
- V = {ε ∈ [ε_min, 0] : χ²(ε) − χ²_min,branch ≤ Δ_V} be the viable set;
- Δχ²(ε) = χ²(ε) − χ²(0);
- Δ_sig be the improvement threshold;
- F(ε) = Δχ²(ε) / Δχ²_w0wa, where Δχ²_w0wa = χ²_min,w0wa − χ²_min,ΛCDM is the targeted-signal magnitude.

| Order | Condition | Result state |
|---|---|---|
| 1 | Δχ²_w0wa > −Δ_sig: no targeted dynamical signal in this likelihood | **INCONCLUSIVE** (nothing to explain) |
| 2 | V ⊂ (−ε_null, 0] | **EXPLANATORY-KILLED** (only the null neighborhood remains viable) |
| 3 | ∃ ε < −ε_null with Δχ²(ε) ≤ −Δ_sig **and** F(ε) ≥ F_min | **EXPLANATORY-SURVIVES** |
| 4 | ∃ ε < −ε_null with Δχ²(ε) ≤ −Δ_sig, but F < F_min everywhere on the branch | **EXPLANATORY-KILLED** (the deformation cannot reach the targeted signal) |
| 5 | otherwise: V extends past the null but no significant improvement | **INCONCLUSIVE** |
| overriding | χ²_min,branch − χ²_min,w0wa ≥ Δ_model (the whole ε ≤ 0 manifold, including ΛCDM, conflicts with the data relative to the flexible comparator) | **MODEL-KILLED** |

## C1 — Minimal firewall (run before any data)

| # | Check | Verdict | Basis |
|---|---|---|---|
| 1 | Dimensions | **PASS** | ε, H₀τ₂ and E_ref are dimensionless; w is dimensionless |
| 2 | Conservation / background consistency | **PASS** | A non-interacting DE component, separately covariantly conserved: ρ_DE ∝ exp[3ε∫R/(1+z)dz], closed form, finite. Matter and radiation are untouched |
| 3 | Well-posed w(z) cosmology | **PASS** | w is explicit, smooth and bounded (−1+ε, −1+ε/2] on z ≥ 0. CAMB/PPF evaluates the background and perturbations (F1–F3). There is no w = −1 crossing on either branch |
| 4 | Initial / boundary data | **PASS** | No new datum. ρ_DE,0 is fixed by flatness closure. The reference background is a declared, priced constant (C0.4) |
| 5 | Equivalence registry | **PASS** | Not equivalent to EQ-01, CPL or wCDM (`EQUIVALENCE_REGISTRY.md`) |
| 6 | Early-universe behaviour (**hostile check**) | **PASS (primary)**; control acceptable | At high z, w → −1+ε. On ε < 0, ρ_DE falls toward the past at least as fast as (1+z)^{3ε/2} on z ≥ 0 (since R ≥ ½), asymptotically as (1+z)^{3ε}, so the **DE fraction is suppressed below ΛCDM's**. Computed Ω_DE(z = 1100) = 1.24 × 10⁻⁹ (ΛCDM), 4.6 × 10⁻¹⁰ (ε = −0.05), 2.3 × 10⁻¹¹ (ε = −0.2), 5.9 × 10⁻¹⁴ (ε = −0.5). At z = 10⁵: ≤ 7 × 10⁻¹⁷. **No early-dark-energy behaviour on the primary branch.** Control ε = +0.2 gives Ω_DE(1100) = 6.7 × 10⁻⁸, also negligible. r_drag is unchanged to 10⁻⁶ (147.069 Mpc) |
| 7 | Parameter identifiability | **PASS (pre-data)** | ε changes H(z) and D_A(z) at BAO redshifts by percent-level amounts for \|ε\| ∼ 0.2 (plumbing log) and shifts θ* at fixed H₀ (feasibility log). It is partially degenerate with H₀ / Ω_m in CMB alone. BAO + SN shape breaks the degeneracy in principle. The actual constraint is data-dependent and not asserted |
| 8 | No post-selection | **PASS** (with disclosure) | The law is transcribed verbatim from a frozen record (2026-06-25 toy). The branch comes from the record's second-law reading. The CPL weighting was chosen by the owner before any model-only computation was shown. Thresholds are owner-locked. Historical DESI-direction discussion is disclosed (§PED) |
| 9 | Nearest standard comparator | **NOTED** | Nearest classes: phantom-side dynamical DE with a fixed shape; CPL w0waCDM (comparator in D3); constant-w phantom wCDM. Distinctive content: a **one-parameter, fixed-shape** family whose CPL image is a line of slope 1.5559 through ΛCDM. The GRUT-load-bearing pieces (τ₂ = H₀⁻¹ and R) are both POSTULATED. Distinctiveness from generic phantom DE is **modest** and is not claimed as a GRUT result |
| 10 | Information price | **PASS** | Fully priced in C0.4; not parameter-free |
| 11 | Origin audit / C-F9 | **PASS** | C0.2 |

**Notes (non-fatal, carried to C3 if ever reached):**
- **N1 — phantom side / NEC.** w < −1 violates the null energy condition for a perfect-fluid reading. The record's reading is an **effective viscous stress** (Π = −3ζH ≤ 0), not a fundamental ghost field. Stability of a microscopic realization is not addressed by this card.
- **N2 — future behaviour.** Because R uses E_ref, |1+w| → 0.408|ε| ≠ 0 as z → −1. This is a constant phantom w in the far future: a finite-time "big rip" for any ε < 0, at times far beyond the data. It is also in tension with the record's statement that (w+1) → 0 at the de Sitter attractor. It lies outside the data domain and is not a C1 kill. A self-consistent-H version would differ, and that is a new card.
- **N3 — perturbation closure.** PPF is a priced projection convention (C0.4). A different closure (e.g. a fluid with c_s² = 1) is a different projection, not a repair.
- **N4 — reference-background semantics.** The prescribed-shape reading is internally consistent in CAMB (F1–F4). There is **no D3 technical inconsistency**.
- **N5 — slope not second-law-fixed.** The record states the w_a slope is not fixed by the second law. In this card the slope is fixed by the **postulated** R and E_ref, not by the origin.
- **N6 — residual-ledger granularity** (C0.1).

**Executor C1 verdict: CARD-01-C1-PASS (with notes N1–N6), subject to owner review.**

## PED — Prior-Exposure Disclosure (Charter §4D; all participants)

| Participant | Exposure |
|---|---|
| **Scientific owner (human)** | Knows the published qualitative DESI DR2 dynamical-DE preference (w0 > −1, wa < 0 in CPL). Directed the CPL weighting choice before any model-only Card #1 computation was shown |
| **Reviewer** | Knows the published qualitative DR2 preference (stated in the owner ruling). Cited the official DESI DR2 chains release page |
| **Executor / orchestrator (Claude Code)** | Has general prior knowledge of the published DESI DR2 BAO + CMB + SN headline (the w0 > −1, wa < 0 preference and its reported significance range) and of the published analysis set-up. Has **not** opened any DESI likelihood, chain or data file. Read only Cobaya package YAML definitions (module and file **names**) |
| **Sub-agents** | **None used** in this turn |
| **Historical reconnaissance (disclosed)** | Card #1's historical grade is **C2-R reconnaissance** against approximate DESI + CMB + SN information (governance Annex A). The scout-1 record compared the toy against "DESI 2024–25 hints w0 ≈ −0.8, wa ≈ −0.6" (`RESULTS_wz.md` §(D)), with sign retractions recorded in `RESULTS_wz_sign.md` |

**No participant has opened the Card #1 target chains or products under this protocol.** This disclosure is a record of exposure, **not a forecast**. Under ER-3, this spec contains no forecast of its own verdict.

## S — Owner-lock slot proposals

> Every entry below is a **CLAUDE-PROPOSAL (non-binding)**. The executor cannot authorize any slot. Final values require the scientific owner's explicit approval in a committed `CARD_01_OWNER_THRESHOLD_RULING.md` that predates any data access. The reviewer may audit and recommend.
>
> The proposals use generic statistical conventions and model-only quantities. **None is derived from the published DESI preferred point.**

### S1 — `epsilon_null`: CLAUDE-PROPOSAL (non-binding)
- **Value:** **ε_null = 0.03**. Null neighborhood: −0.03 < ε ≤ 0.
- **Rationale:** a model-intrinsic "percent-level indistinguishability from Λ". At |ε| = 0.03:
  - max pointwise |1+w| on the registered window [0, 2] is 0.9 × 0.03 = 0.027 (below 3 %);
  - the CPL image is |1+w0| = 0.013, |wa| = 0.020.
  A deviation that small is at or below conventional percent-level dark-energy systematics floors, whatever any particular dataset prefers.
- **Sensitivity:** the explanatory verdict is directly sensitive to it. Tightening (e.g. 0.01) shrinks the null neighborhood and makes it easier to escape the kill with a tiny ε. Loosening (e.g. 0.05) demands a larger deviation before SURVIVES is possible.
- **If tightened:** more borderline outcomes move from EXPLANATORY-KILLED to INCONCLUSIVE or SURVIVES.
- **If loosened:** outcomes move toward EXPLANATORY-KILLED.

### S2 — Primary statistic: CLAUDE-PROPOSAL (non-binding)
- **Definition:** the **profile likelihood** on the primary branch, Δχ²(ε) = χ²_min(ε) − χ²_min(ε = 0). All other parameters are minimized, under the identical S6 likelihood (§C0.5(b), D3 steps 2–3).
- Because ε = 0 sits on the boundary of [ε_min, 0], significance uses the **Chernoff ½χ²₀ + ½χ²₁** reference.
- The Bayesian posterior on ε is a secondary report only.
- **Rationale:** prior-independent; directly comparable with the nested null and the CPL comparator (C-F7); well defined at the boundary.
- **Sensitivity:** depends on minimizer convergence. Hence the 3-start protocol and the 0.2 χ² flag.
- **If replaced by a posterior criterion:** results would depend on the S7 prior volume (the Bayesian Occam penalty grows with |ε_min|).

### S3 — EXPLANATORY-SURVIVES threshold: CLAUDE-PROPOSAL (non-binding)
- **Value.** SURVIVES requires a point ε < −ε_null with both:
  - (a) **Δχ²(ε) ≤ −4** (Δ_sig = 4: one-sided p ≈ 0.023 under the Chernoff mixture, about 2σ), **and**
  - (b) **F(ε) = Δχ²(ε)/Δχ²_w0wa ≥ 0.5**, i.e. it recovers at least half of the likelihood improvement the two-parameter CPL model achieves over ΛCDM in the same likelihood.
- **Viable set:** Δ_V = 3.84 (95 %, 1 dof profile).
- **Rationale:**
  - (a) is the conventional 2σ evidence level for one added parameter;
  - (b) operationalizes "meaningfully accounts for the targeted signal" without reference to any particular (w0, wa) point. It compares like with like (χ² improvements), and the card has one parameter against CPL's two.
- **Sensitivity:** (b) depends on Δχ²_w0wa, which is itself data-determined. If the targeted signal is weak, F is noisy, which is why decision-tree row 1 routes "no signal" to INCONCLUSIVE.
- **If tightened** (Δ_sig = 9 or F_min = 0.8): SURVIVES becomes harder, and more outcomes fall to INCONCLUSIVE / KILLED.
- **If loosened** (Δ_sig = 1 or F_min = 0.25): a marginal improvement could be labeled SURVIVES. That would be weak evidence and open to misreading as support.

### S4 — MODEL-KILLED threshold: CLAUDE-PROPOSAL (non-binding)
- **Value:** **χ²_min,branch − χ²_min,w0wa ≥ 25**. The best point of the whole ε ≤ 0 manifold, ΛCDM included, is worse than the flexible CPL comparator by Δχ² ≥ 25 (≈ 4.6σ for the 1-dof difference).
- **Rationale:** MODEL-KILLED is reserved for the much stronger case (Charter §4C). Since ΛCDM is nested, this outcome would simultaneously mean ΛCDM itself is excluded at that level in this likelihood. A 5σ-class convention is appropriate for a claim of that weight.
- **Sensitivity:** absolute goodness of fit with CMB likelihoods is ill-defined, so a relative criterion is used. That makes the outcome depend on the CPL comparator's adequacy.
- **If tightened** (e.g. 16): MODEL-KILLED could be triggered by a ~4σ ΛCDM tension.
- **If loosened** (e.g. 36): even stronger tension is needed.

### S5 — INCONCLUSIVE band: CLAUDE-PROPOSAL (non-binding)
- **Definition:** INCONCLUSIVE iff either:
  - (i) **Δχ²_w0wa > −4** (no targeted dynamical signal at 2σ in this likelihood); or
  - (ii) V extends below −ε_null but no ε < −ε_null reaches Δχ² ≤ −4 (decision-tree rows 1 and 5).
- **Rationale:** separates "the data cannot discriminate" from "the data reject the explanation". Without a targeted signal, the explanatory question has no object.
- **Sensitivity:** shares Δ_sig with S3. Raising Δ_sig widens the band.
- **If tightened:** fewer INCONCLUSIVE outcomes, and more forced verdicts on weak evidence.
- **If loosened:** more INCONCLUSIVE outcomes.

### S6 — Dataset combination(s): CLAUDE-PROPOSAL (non-binding)
- **Value:** primary = **DESI DR2 BAO (all tracers) + CMB (Planck PR4 CamSpec high-ℓ TTTEEE + Planck 2018 low-ℓ TT + low-ℓ EE + Planck PR4 + ACT DR6 lensing) + DES-Y5 SN**.
  - Robustness: the same with **Pantheon+**, and with **Union3**.
  - The DES-Y5 "Dovekie" recalibration is optional robustness, at owner discretion.
  - The verdict is decided on the **primary** only. Robustness rows are reported and flagged if their result state differs.
- **Rationale:** matches the principal published dynamical-DE comparison set-up (per general knowledge, **TO CONFIRM** against official documentation at preflight, because the documentation host is blocked here). It keeps DESI's own SN-robustness structure without choosing the SN sample after the fact.
- **Sensitivity:** the dynamical-DE signal strength is known (published headline, disclosed) to vary across SN samples. That is why the primary is fixed now and robustness rows cannot change the verdict.
- **If the CMB composition cannot be confirmed:** the owner must rule on the composition before data access, or the run is ACCESS-BLOCKED.
- **Provenance:** the owner must rule on question O3.

### S7 — ε prior / ε_min: CLAUDE-PROPOSAL (non-binding)
- **Value:** **ε ∈ [−1, 0]**, uniform (used only for the secondary posterior); profile grid down to ε_min = −1. Control: [0, +1], separate.
- **Rationale:** |ε| ≤ 1 keeps the toy's "deviation amplitude ≤ the relaxor's full storage response", bounding w ∈ (−2, −1) at all z. The historical toy explored |ε| ≤ 0.6. Data constraints are expected well inside this range, so the boundary should not bind.
- **Sensitivity:** the profile-likelihood verdict (S2–S5) is prior-independent unless the best fit hits ε_min. The posterior report depends on the range.
- **If tightened** (e.g. [−0.3, 0]): the risk is truncating a best fit, which would need to be flagged.
- **If loosened** (e.g. [−2, 0]): it admits w < −2 at high z. Physically extreme, but harmless to the profile.

### S8 — CPL diagnostic credible level: CLAUDE-PROPOSAL (non-binding)
- **Value:** report the **68 % and 95 %** highest-posterior-density contours of the official DESI DR2 w0wa chains for the S6 primary combination, if reachable. The overlay statement uses **95 %**. Also report the posterior mass in the quadrant {w0 > −1, wa < 0} and the minimum posterior-density level crossed by the frozen Card #1 line.
- **Rationale:** standard reporting levels. The CPL overlay is diagnostic only and cannot change the result state (Annex B §4).
- **Sensitivity:** none for the verdict, by construction.
- **If changed:** only the wording of the diagnostic changes.

## Technical obstacles to D3 (pre-data)

| ID | Obstacle | Status |
|---|---|---|
| **O1** | **Network.** The official DESI documentation host is egress-blocked (WebFetch to `www.desi.lbl.gov`). The likelihood data (Cobaya `bao_data` GitHub release, Planck PR4 / PR3 likelihood files, ACT DR6 lensing package and data, DES-Y5 / Pantheon+ / Union3 files) all require downloads that have **not** been attempted. **ACCESS-BLOCKED is a live risk** for at least the D1 chains | open; resolved only by the preflight (after the owner ruling) |
| **O2** | **Compute.** The full CMB likelihood costs ≈ 1.7 s per CAMB evaluation. The profile requires ~20 grid points × 3 starts × O(10³) evaluations, plus a CPL comparator and robustness rows: roughly a day of wall time on 4 cores. It must run detached (setsid / nohup), as with PM2. This is feasible, not blocking | noted |
| **O3** | **"Official product" definition** for non-DESI likelihoods and for the Cobaya-distributed DESI BAO files | **owner ruling needed** |
| **O4** | **D1 needs DESI's own w0wa chains**, a separate product from the likelihoods. If unreachable, D1 is reported as unavailable with no substitute. D3 is unaffected | noted |
| **O5** | **Accuracy settings.** CAMB accuracy (lens_potential_accuracy ≥ 1, ℓ_max per the CamSpec requirements) must be set from each likelihood's documented requirements. They are to be frozen in the run configuration **before** any likelihood value is inspected | to freeze at run-configuration commit (post-ruling, pre-data) |

**No technical inconsistency** was found between the frozen prescribed-shape model and a standard CAMB + Cobaya stack (N4).

## Commit-boundary statement

This document, the registries, the ledger and `code/` (model-only scripts and logs) constitute the **CARD-01 PRE-DATA SPEC BOUNDARY**.

**After this commit (all forbidden):**
- no network preflight;
- no DESI product access;
- no chain download;
- no D3 run;
- no `CARD_01_OWNER_THRESHOLD_RULING.md`;
- no self-approval of any proposal.

The next valid step is the owner threshold ruling.
