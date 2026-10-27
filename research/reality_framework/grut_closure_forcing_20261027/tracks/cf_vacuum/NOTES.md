# Track cf_vacuum: closure as vacuum equilibrium (P3), frozen and confronted with the GRUT branches

27 October 2026. Claude, closure-forcing parallel pass (issue #13 context). Brief: `../../BRIEF.md`.
Records built on: `grut_q_vacuum_construction_20261021/tracks/vacua/NOTES.md` ("vacua §n") and its verifier
`tracks/v_vacua/VERIFY.md` ("v_vacua Cn"), both draft/unreviewed; Codex's foundation gate (2908e91) and exact
architecture (04c53f5); the heavy-window and portal records in `grut_formation_review_20261025/`.

**Labels.** [D] derived here; [C] computed by `code/` (output in `results/`); [Q] stated by the cited source (read level in
`results/sources.json`); [E] estimate; [A] argued; [R] recalled.

**Scope.** No validation, no exclusion of the real world, no new theory. Every comparison with an existing measurement
is **retrospective**. "Held out" means data that do not exist yet. Every postulate below is an **added assumption**; none
is a consequence of the word "closed".

**Reproduction.** `cd code && python3 -B run_all.py` (about 8 min). Needs the SMDR v1.3 build that v_vacua left in the
scratchpad (`.../scratchpad/qtheory/v_vacua_dl/smdr/SMDR-1.3`; `smdr_common.py` and `smdr_driver.c` are v_vacua's,
copied unchanged; `sm_rge.py` is track vacua's 3-loop code, copied unchanged). Downloads (EEGLS and Lebedev TeX, PDG
2026 QCD review) are in `.../scratchpad/closure/cf_vacuum/dl/`, not in the repository.

---------------------------------------------------------------------------------------------------------------------

## 0. Answer in one screen

| # | Result | Label |
|---|---|---|
| 1 | **Frozen relation (PDG 2026 inputs).** Under P_q read as "the final equilibrium is the global minimum" (or the worldline end-point, which reduces to it when Γ > 0; vacua §3), with the SM valid to M_Pl: **M_t ≤ 171.00 + 0.494 (M_h − 125.13) + 0.355 (α_s − 0.1180)/0.0009 GeV (±0.11 matching, ±0.03 criterion)**, μ* = 5.0e17 GeV. Pipeline reproduces v_vacua's 171.026 at the PDG 2025 point to 1e-5 GeV and Hiller et al.'s 171.04. Written to `results/frozen_relation.json`. | [C] |
| 2 | **Retrospective tension.** PDG 2026 direct M_t = 172.60 ± 0.27: Δ = 1.60 GeV, **z = 3.45** (p = 3e-4); with a 0.3 GeV MC-to-pole error 2.90; with 0.6 GeV 2.11. ATLAS+CMS 172.52 ± 0.33: 3.04 (2.61 with 0.3 MC). Cross-section pole mass 172.1 ± 0.6: **z = 1.56**. FLAG α_s = 0.1183 ± 0.0007: 3.65. z = 0 needs m_pole − m_MC = −1.60 GeV, or α_s = 0.1221, or M_h = 128.43 GeV. | [C] |
| 3 | **The brief's formula is not Buttazzo et al.'s.** "171.53 + 0.15(M_h−125.15) − 0.34 Δα_s/0.0007" is not in 1307.3536 and has the wrong sign of the α_s term; it gives 171.72 at PDG 2026. The paper's own linear formula, inverted, is M_t < 173.34 + (M_h − 129.6)/2.0 + 0.25 Δα_s/0.0007 (= 170.96 here); its quoted "171.53 ± 0.42" does not follow from it (vacua §5.3, v_vacua row 26). | [Q][C] |
| 4 | **Instability scales (M_t = 172.60).** MSbar λ = 0 at **Λ_λ = 4.9e10 GeV** (2.0e10 to 1.4e11 for ±1σ in M_t; 3.9e11 with the cross-section mass); Landau-gauge λ_eff = 0 at Λ_I = 3.5e11 (2-loop); barrier top h_max = 2.7e11; Λ_V = 2.6e10 GeV. V_eff(2.435e18) = −7.5e70 GeV⁴ = 3e117 ρ_Λ. | [C] |
| 5 | **P_q forces Higgs-coupled new physics below ~Λ_λ if metastability is confirmed.** Tree-threshold singlet (EEGLS 1203.0237): stability needs δ = λ_HS²/λ_S ∈ [δ_min, δ_max] ≈ [0.007–0.16, 0.12–0.36] and M_S ≲ 4.5e10 GeV (none at 1e11); then sin θ · M_S ≈ √(2δ) v ∈ **[≈30, ≈210] GeV**, independent of the singlet vev. Measurable only if M_S ≲ few TeV. | [C] |
| 6 | **GRUT's own fields do not provide the escape at their reference values.** φ and X (loop-only, light): portal needed 0.117 and 0.175, against invisible-width bounds 5.1e-3 and 7.1e-3 (factor 23–24; Γ(h→ss) ≈ 0.27–0.29 GeV ≈ 70 Γ_SM). Y: no vev, M_Y = 4.0e12 ≫ Λ_λ. Σ at the heavy reference (λ_Σ = 1, M_Σ = 7.5e12 = 153 Λ_λ): no window for any perturbative portal. | [C] |
| 7 | **A Σ window exists only far from the reference.** With the stated vev h = 5.33e12 GeV: M_Σ ≲ 4e10 GeV ⇒ **λ_Σ ≲ 3e-5** (reference 1), λ_HΣ = 7e-12 … 6e-3, θ = 4e-1 … 3e-9 (θ < 2e-4 for M_Σ ≥ 1e6). | [C] |
| 8 | **P_q then bites the hidden sector.** Absolute stability of the (Σ, Y) potential at tree level needs λ_Σ ≥ 9.05e-4 (λ_Y = 1, λ_ΣY = 0) because the trilinear μΣ*Y² makes Y tachyonic at \|Σ\| > 2.22 \|Σ\|_vac. The reference passes (margin 1100); every Σ window point fails, with a Y-condensed vacuum 2e98–2e100 ρ_Λ deeper, **unless λ_ΣY ≥ 0.047–0.057** (≈ μ²/M_Y²) is added. | [D][C] |
| 9 | **Costs of the consistent corner** [E]: λ_Σ tuned against the one-loop Y threshold (−1.37e-5; with λ_ΣY ≈ 0.057 total \|δλ_Σ\| ≈ 1.3–3.4e-4) at 0.1 to 1e-11; m_Σ² tuned at 5e-3 to 1e-12; Higgs mass parameter shifted by λ_HΣ h²/2 (tuning 1e-11 to 1e-19, absent when λ_HΣ = 0). Σ's radial mode becomes lighter than H_i (M_Σ ≤ 0.34 H_i), so the heavy-window results computed at m_Σ² = 3210 H_i² are not established there. | [E][A] |
| 10 | **Standard-physics obligation found on the way (independent of P_q).** The heavy window has H_i = 1.33e11 GeV ≈ 0.5 h_max at the central M_t (random-walk h_rms ≈ 1.6e11 GeV over 60 e-folds), so it needs a Higgs–curvature coupling ξ_H ≳ 0.05–0.07 (2011.03763 [Q abstract]) or a stabiliser. ACTION.md allows ξ_H but does not state it. | [C][E] |
| 11 | **Verdict (gate).** P_q(a) + SM: a strict, testable restriction (MIXED: testable, but the relation is the standard stability bound and the reading is an added assumption; retrospective tension 3.0–3.5σ, 1.6σ with the cross-section mass; held out FCC-ee ≈ 12σ, HL-LHC + lattice ≈ 5.5σ if central values stay). P_q(a) + GRUT at its reference: inherits the same tension, no internal escape. P_q(a) + GRUT Σ window: consistent only in a tuned corner, no measurable signature. It is the best *testable* candidate found, but it is not GRUT-native and not forced by closure. | [A] |

---------------------------------------------------------------------------------------------------------------------

## 1. The postulates (gate item 1)

The q-theory action and P_q are Klinkhamer–Volovik's (0711.3170, 0806.2805; KV 1102.3152 §6 for the "global
equilibrium" choice). QV1: S = ∫√−g[½K(q)R − ε(q)] + S_m; one global integration constant μ; P_q sets μ = μ₀[ε_∞]
for "the final T → 0 vacuum" (HANDOFF.md). With several matter vacua, Λ_i = F(V_i − ε*(μ)), F(0) = 0, F′ > 0
(vacua §1, CONFIRMED by v_vacua rows 1–4).

- **P3a (global/worldline reading).** μ is fixed by equilibrium of the **global minimum of the full renormalized
  effective potential of all sectors** (SM, hidden, and any UV completion); the observed ρ_Λ is a separate small offset
  that may not absorb vacuum-energy differences. Equivalent, when Γ > 0, to "the vacuum every worldline ends in"
  (vacua §3.3, v_vacua §5). Constrains: the ordering of vacuum energies, i.e. (M_t, M_h, α_s) in the SM and the scalar
  couplings of any extension.
- **P3d (coexistence/MPP).** All final vacua coexist at P = 0: equality M_t = M_t,crit (Froggatt–Nielsen hep-ph/9511371).
  Not implied by the action (vacua §3.4).
- **P3v (volume reading, de Sitter).** The final vacuum is ours by volume (Γ/H⁴ < 9/4π): **no SM relation** (vacua §3.3).

None of these readings is chosen by "no exterior": GRUT's closure does not decide global versus final, nor worldline
versus volume measure (vacua §6.5; FOUNDATION_GATE: "vacuum equilibration [is an] extra assumption").

## 2. What P3a forces, and the baseline it forbids (gate items 2–3)

**Compatibility [Q/C via vacua, v_vacua].** The q-equations pass a Bianchi check; μ is conserved and continuous across
matter walls; P_q is a final-state (boundary) condition, not enforced by dynamics (HANDOFF; ledger §1.3). Fixed point
under nucleation: with Γ > 0 the worldline reading has the unique self-consistent fixed point "true vacuum" (vacua §3.3).

**Derivation of the restriction [D].** Under P3a, Λ_ours = F(V_ours − V_true) (+ the separate offset). If a deeper
vacuum exists anywhere in the validity range, Λ_ours ≥ |V_true − V_ours|. In the SM at PDG 2026 central values the
potential at h = 2.435e18 GeV is −7.5e70 GeV⁴ [C], so Λ_ours ≥ 3e117 ρ_Λ, which is retrospectively excluded. The
tolerance is δλ ~ 4ρ_Λ/h*⁴ ~ 1e-117, so P3a ⇔ "no deeper vacuum": λ(μ) ≥ 0 up to M_Pl, the boundary being
λ(μ*) = β_λ(μ*) = 0 (Bednyakov et al. 1507.08833).

**Forbidden baseline [D].** The SM valid to M_Pl with (M_t, M_h, α_s) on the metastable side — e.g. the PDG 2026
central point — is a legitimate standard-physics model: its lifetime is ≫ the age of the universe (vacua §2.3) and
nothing in standard physics excludes it. P3a excludes it. The restriction S_C ⊊ S_base is strict, and it is not a
relabelling of the same distributions: it removes a region of measurable parameter space.

## 3. The frozen SM relation (gate item 4) — `code/c1_fit_grid.py`, `code/c2_sm_relation.py`

### 3.1 Inputs [Q]

| Input | Value | Source |
|---|---|---|
| M_h | 125.13 ± 0.11 GeV (S = 1.5) | PDG 2026 S126M (pdglive, 2026 edition) |
| M_t direct (MC-template) | 172.60 ± 0.27 GeV (S = 1.5) | PDG 2026 Q007TP |
| M_t pole from σ(tt̄) | 172.1 ± 0.6 GeV | PDG 2026 Q007TP4 |
| M_t ATLAS+CMS | 172.52 ± 0.33 GeV | arXiv:2402.08713 (abstract) |
| α_s(M_Z) | 0.1180 ± 0.0009 (FLAG 2024: 0.1183 ± 0.0007) | PDG QCD review rev. Aug 2025, eqs. (9.25), (9.23) |
| Other SM inputs | SMDR's PDG-2025 reference values | Martin–Robertson 1907.02500 |

PDG 2026 moved M_h from 125.20 to 125.13 and M_t from 172.56 ± 0.31 to 172.60 ± 0.27; α_s unchanged.

### 3.2 Method and validation [C]

SMDR full matching at Q = M_t on 11 fits (quadratic in M_t, linear in M_h, α_s; interpolation residuals ≤ 1e-5),
SMDR running at loop order 3 and 5, criterion min_{μ≤M_Pl} λ(μ) = 0 — the verifier's pipeline (S). Checks: on
v_vacua's own fit grid it returns 171.0264 at (125.20, 0.1180) (v_vacua C1: 171.026); at (125.7, 0.1185) it gives
171.47 (Bednyakov et al. 171.44). Loop 2 − loop 3: −0.110 GeV; loop 5 − loop 3: +0.005 GeV.

### 3.3 The line

    M_t ≤ 170.997 + 0.494 (M_h − 125.13) + 0.355 (α_s − 0.1180)/0.0009 GeV,   ±0.113 (matching) ±0.032 (criterion)
    μ* = 4.96e17 GeV.   Equivalently M_h ≥ 128.43 GeV at M_t = 172.60, or α_s ≥ 0.1221.

The equality is the MPP/criticality form (P3d).

### 3.4 Tension, and the top-mass interpretation (retrospective) [C]

σ_tot = σ(M_t) ⊕ 0.494σ(M_h) ⊕ 0.355σ(α_s)/0.0009 ⊕ 0.113 ⊕ extras.

| M_t input | extras | Δ (GeV) | σ_tot | z | p (one-sided) |
|---|---|---|---|---|---|
| direct 172.60 ± 0.27 | none | 1.603 | 0.463 | **3.46** | 2.7e-4 |
| direct | criterion + 0.3 MC→pole | 1.603 | 0.553 | **2.90** | 1.9e-3 |
| direct | criterion + 0.6 MC→pole | 1.603 | 0.759 | 2.11 | 0.017 |
| direct | criterion + 0.088 matching shift (v_vacua C3) | 1.603 | 0.472 | 3.39 | 3.5e-4 |
| ATLAS+CMS 172.52 ± 0.33 | criterion | 1.523 | 0.501 | 3.04 | 1.2e-3 |
| ATLAS+CMS | criterion + 0.3 MC | 1.523 | 0.584 | 2.61 | 4.6e-3 |
| σ(tt̄) pole 172.1 ± 0.6 | criterion | 1.103 | 0.709 | **1.56** | 0.060 |
| direct, FLAG α_s 0.1183 ± 0.0007 | criterion | 1.485 | 0.407 | 3.65 | — |

**MC versus pole.** Treat m_pole = m_MC + d. z(d) = 3.45 + 2.15 d/GeV: z = 2.16 at d = −0.6, 4.75 at d = +0.6, and z = 0
only at d = −1.60 GeV. Butenschoen et al. 1608.01318 [Q abstract] find m_MC differs from the pole mass by 900/600 MeV
(NLL/NNLL) and agrees with the MSR mass at R = 1 GeV; the sign relevant here is not settled at the read level used,
so the 0.3 and 0.6 GeV rows are bracketing errors, not shifts. The direct and cross-section masses differ by 0.5 GeV,
and that difference is the main reason the tension spans 1.6σ to 3.5σ.

**Held out (illustrative only; data do not exist).** FCC-ee (δm_t ≈ 40 MeV, δα_s = 0.00012, δM_h = 4 MeV; precisions
from vacua §5.6 [Q/R]): σ_tot = 0.133 GeV, z = 12.1 if the direct central value stays (8.3 at the cross-section value),
17.9 with halved matching errors. HL-LHC + lattice (δm_t ≈ 0.2 [R], δα_s = 0.00045): z = 5.5.

**Newer literature [Q abstracts].** arXiv search 2024-01 to 2026-10 (31 hits; `results/sources.json`) found no update of
the SM line beyond Hiller et al. 2401.08811 (PDG 2024: 171.04 λ-criterion, 171.10 effective; 5σ needs errors 2–3×
smaller). Contrast, also credited: the "metastability bound" programme (2408.10297, 2503.03825, 2503.22787) proposes
the **opposite** selection — that the Higgs sits near the metastability edge — so a confirmed metastable vacuum would
favour that rival principle over P3a.

## 4. The escape P3a forces if metastability is confirmed — `code/c3_threshold.py`

### 4.1 Instability scales (M_t = 172.60) [C]

| Scale | Value (GeV) | ±1σ M_t (172.33 / 172.87) | σ(tt̄) mass 172.1 |
|---|---|---|---|
| Λ_λ (MSbar λ = 0) | 4.9e10 | 1.4e11 / 2.0e10 | 3.9e11 |
| Λ_I (λ_eff = 0, Landau, 2-loop) | 3.5e11 | 1.0e12 / 1.4e11 | 3.1e12 |
| h_max (barrier top, 2-loop) | 2.7e11 | 8.0e11 / 1.0e11 | 2.4e12 |
| Λ_V = V_max^{1/4} | 2.6e10 | 7.4e10 / 1.0e10 | 2.1e11 |

At Buttazzo's inputs the same code gives Λ_λ/Λ_V = 2.0 and Λ_I/Λ_V = 12.8 (theirs: ≈ 2 and ≈ 13), log10 Λ_V = 9.61
(theirs 9.5, the 0.1 being the 2013-versus-SMDR matching shift).

### 4.2 Generic tree-level threshold (EEGLS 1203.0237; Lebedev 1203.0156) [Q][D][C]

EEGLS normalisation V = λ_H(H†H − v²/2)² + λ_S(S†S − w²/2)² + 2λ_HS(H†H − v²/2)(S†S − w²/2); integrating S out gives
λ_SM = λ_H − δ, δ = λ_HS²/λ_S, M_S² = 2λ_S w². Their one-loop RGEs (scalar parts) were re-derived here from the
general formula 16π² dV/d ln μ = ½ Tr M⁴ with sympy: 24λ_H² + 4λ_HS², 8λ_HS² + 20λ_S², 4λ_HS(3λ_H + 2λ_HS + 2λ_S)
— identical [D]. Real-singlet and complex-X versions are in `threshold.json` (A).

**Stability test [C].** RG-improved tree potential at μ = h: SM λ ≥ 0 below M_S; above M_S, V minimised over |S|
(valley: (λ_H − δ)x²; S = 0 branch: λ_H x² + λ_S w⁴/4 − λ_HS x w²) must be ≥ 0 up to M_Pl; all quartics < 4π to M_Pl.
SM part 3-loop, singlet part 1-loop, tree matching at M_S.

**One useful identity [D].** Above M_S the valley quartic λ_v = λ_H − δ runs as β_SM(λ_v) + δ(24λ_v + 8δ + 8λ_S −
16λ_HS)/16π². So with a small portal (h_c = w√(λ_S/λ_HS) ≫ M_S) the valley stays SM-like apart from the 8δ² term; the
tree shift alone helps only once the singlet stops tracking (|S| → 0 at h_c). The numerics include both effects.

Results: windows exist for M_S from 300 GeV to 4e10 GeV, none at 1e11 (λ_SM < 0 below M_S); δ_min = 0.007 (M_S =
TeV, small portal) rising to 0.15 at 3e10; δ_max = 0.12–0.36 from perturbativity; λ_HS = 1 has no perturbative window.
EEGLS: "the mechanism can be operative even for a very heavy singlet, as long as its mass is smaller than the instability
scale" [Q TeX]. Using the λ_eff criterion instead of MSbar λ would raise the mass edge by about ×7 (Λ_I/Λ_λ).

**What would be measurable [C].** In the heavy limit θ ≈ λ_HS v/(λ_S w) = √(2δ) v / M_S, independent of w (exact 2×2
diagonalisation agrees to 4e-9 relative at M_S = 1e6). So the window gives **sin θ · M_S ∈ [≈30, ≈210] GeV** (lower edge rising with M_S).
Higgs couplings shift by θ²/2: measurable at HL-LHC (|sin θ| ≳ 0.15–0.2 [E/R]) only for M_S ≲ 1 TeV, at FCC-ee
(|sin θ| ≳ 0.05 [E/R]) for M_S ≲ 0.6–4 TeV. Above that the escape is invisible (Lebedev: "mixing … can be as small as
1e-7" [Q TeX]). P3a therefore does **not** predict a measurable Higgs signal; it predicts a band that is empty of
signal unless the stabiliser is light.

### 4.3 GRUT's fields one by one [C][Q]

| Field (reference) | Mechanism | Requirement | GRUT value or bound | Verdict |
|---|---|---|---|---|
| φ (v_φ = 0.924 GeV, λ_φ = 0.03) | v_φ ≪ v: no tree threshold (mixing raises m_h, so it destabilises slightly) [D]; loops only | b ≥ 0.117 (b φ²H†H) | b ≤ 5.05e-3 (ATLAS B_inv < 0.107; rp_portals [C]); detector-level target 1–1.5e-10 | excluded by 23× (Γ(h→φφ) = 0.27 GeV) |
| X (m_X = 92 MeV, no vev) | loops only | λ_HX ≥ 0.175 | λ_HX ≤ 7.14e-3 | excluded by 24× |
| Y (M_Y = 4.0e12, no vev) | loops above M_Y only | M_Y < Λ_λ | M_Y = 81 Λ_λ | cannot act on (Λ_λ, M_Y) |
| Σ at reference (λ_Σ = 1, h = 5.33e12) | tree threshold | M_Σ ≲ 4e10 | M_Σ = 7.54e12 = 153 Λ_λ | no window for any portal |
| Σ, λ_Σ lowered | tree threshold | see §5 | — | window, with costs |

The inflaton, right-handed neutrinos or other singlets are not in GRUT's stated field content (THEORY.md §1); using
them would be supplied physics.

## 5. The only GRUT-internal window: Σ with a small quartic — `code/c3_threshold.py`, `code/c4_hidden_and_costs.py`

### 5.1 The window at the stated vev h = 5.334e12 GeV [C]

λ_Σ = M_Σ²/(2h²), λ_HΣ (GRUT normalisation λ_HΣ|H|²|Σ|²) = 2λ_HS, δ = λ_HΣ²/(4λ_Σ).

| M_Σ (GeV) | λ_Σ | δ range (+ portal) | δ_min (− portal) | λ_HΣ (+) | θ | h_c (GeV) |
|---|---|---|---|---|---|---|
| 1e3 | 1.8e-20 | 0.0078–0.136 | 0.038 | 2.3e-11 – 9.8e-11 | 0.031 – 0.13 | 2.1e8 |
| 1e6 | 1.8e-14 | 0.0102–0.211 | 0.069 | 2.7e-8 – 1.2e-7 | 3.5e-5 – 1.6e-4 | 6.1e9 |
| 1e8 | 1.8e-10 | 0.0143–0.265 | 0.092 | 3.2e-6 – 1.4e-5 | 4.2e-7 – 1.8e-6 | 5.6e10 |
| 1e9 | 1.8e-8 | 0.065–0.296 | 0.107 | 6.8e-5 – 1.4e-4 | 8.9e-8 – 1.9e-7 | 1.2e11 |
| 1e10 | 1.8e-6 | 0.120–0.332 | 0.129 | 9.2e-4 – 1.5e-3 | 1.2e-8 – 2.0e-8 | 3.3e11 |
| 4e10 | 2.8e-5 | 0.162–0.357 | 0.158 | 4.3e-3 – 6.3e-3 | 3.5e-9 – 5.2e-9 | 6.1e11 |
| 1e11, 7.54e12 (ref.) | 1.8e-4, 1 | none | none | — | — | — |

So P3a + GRUT requires **λ_Σ ≲ 3e-5** — 4.5 to 20 orders below the reference — and a portal λ_HΣ that the reference
leaves at zero (and that REVIEW_REHEATING's portal ledger never considered).

### 5.2 P3a applied to the hidden sector itself [D][C]

P3a is a statement about the global minimum of the *whole* potential, so it also constrains the hidden scalars. With
aligned phases, μΣ*Y² + h.c. → −2μ|Σ||Y|², so the Y mass² along Σ is m²(s) = M_Y² + λ_ΣY s² − 2μs. At our vacuum
m²/M_Y² = 1 − r = 0.55 (r = 0.45, as in hw_gates). For s > M_Y²/(2μ) = 2.22 s_vac, Y condenses and
V_min(s) = λ_Σ(s² − h²/2)² − m⁴(s)/(4λ_Y). Requiring V_min ≥ 0 for all s:

| λ_Y | λ_ΣY = 0 | 0.01 | 0.03 | 0.05 | μ²/M_Y² = 0.0570 |
|---|---|---|---|---|---|
| 1 | λ_Σ ≥ 9.05e-4 | 6.1e-4 | 2.0e-4 | 1.3e-5 | 0 |
| 2 | 4.5e-4 | 3.1e-4 | 1.0e-4 | 6.7e-6 | 0 |
| 4π | 7.2e-5 | 4.9e-5 | 1.6e-5 | 1.1e-6 | 0 |

- The reference (λ_Σ = λ_Y = 1) satisfies it with margin ≈ 1100. P3a in the hidden sector, like gate 4's κ bound, is
  then a consistency condition the reference already meets: it adds no test there.
- Every point of the §5.1 window violates it unless λ_ΣY ≥ 0.047 (M_Σ = 4e10) to 0.057 (M_Σ ≤ 1e8), with λ_Y = 1. The
  Y-condensed vacuum would be 2e98–2e100 ρ_Λ below ours.
- So the P3a + GRUT corner is (λ_Σ ≲ 3e-5, λ_ΣY ≳ 0.05, δ in the §5.1 band), plus κ within gate 4's bound.

### 5.3 What the corner costs [E][A]

- **Quartic.** One-loop Y threshold: δλ_Σ = −μ⁴/(24π²M_Y⁴) = −1.37e-5 [D: Coleman–Weinberg expansion with Y masses
  M_Y² ± 2μ|Σ|]; with λ_ΣY = 0.057 the λ_ΣYμ²/(8π²M_Y²) and λ_ΣY² log pieces bring |δλ_Σ| to 1.3–3.4e-4. Against
  λ_Σ = 1.6e-5 … 1.8e-14 this is a tuning of 0.1 to 1e-11.
- **Mass.** δm_Σ² ≈ μ² ln(M_Y²/M_Σ²)/(8π²) ≈ 1–5e23 GeV² against M_Σ² = 9e20 … 1e12: tuning 5e-3 to 1e-12.
- **Higgs.** λ_HΣ h²/2 = 3e14 … 7e22 GeV² against m_h²/2: tuning 1e-11 … 1e-19. With λ_HΣ = 0 (the reference;
  radiatively consistent because no coupling links the sectors) this tuning is absent, so P3a reintroduces a
  hierarchy cost that GRUT's decoupled hidden sector did not have.
- **Heavy-window dynamics.** M_Σ ≤ 4e10 = 0.34 H_i: the Σ radial mode is light during inflation and reheating, whereas
  hw_gates and the yield 1.00 Y_Q used m_Σ² = 3210 H_i² (frozen radial mode). Those results are not established in
  the corner, and the light radial mode is a new degree of freedom for the release, the strings (gate 3 used β =
  2λ_Σ/36e² = 5.6e12; at λ_Σ = 1e-8 it is 5.6e4, still type II) and debris. This is an obligation, not a refutation.
- **Measurable signature.** None for M_Σ ≥ 1e5 GeV (θ < 2e-3, θ²/2 < 2e-6). The light corner (M_Σ ~ TeV, θ =
  0.03–0.13) would be visible at HL-LHC/FCC-ee, but needs λ_Σ ≈ 2e-20, a tuning of about 1e-15 [E].

### 5.4 Inflation scale against the barrier (standard physics, independent of P_q) [C][E][Q]

| M_t | h_max (GeV) | H_i/h_max | H_i/Λ_λ |
|---|---|---|---|
| 172.60 | 2.7e11 | 0.50 | 2.7 |
| 172.87 (+1σ) | 1.0e11 | 1.3 | 6.6 |
| 172.1 (σ(tt̄)) | 2.4e12 | 0.056 | 0.34 |

With H_i = 1.33e11 GeV, a massless Higgs random-walks to √N H/2π ≈ 1.6e11 GeV over 60 e-folds ≈ 0.6 h_max [E].
Unless the Higgs gets a Hubble-scale mass, the heavy window falls into the AdS-patch problem (Espinosa et al.
1505.04825 [Q abstract]). The usual remedy is ξ_H R|H|² with ξ_H ≳ 0.051–0.066 (Mantziris–Markkanen–Rajantie
2011.03763 [Q abstract], for their inflation models). In ACTION.md's sign convention (F ⊃ −2ξ_H H†H) that means
ξ_H > 0 — the opposite sign to the ξ_X = −5/6 chosen for X. The low-scale candidates (H_i = 4.35e11) need it more.
This obligation exists whether or not P_q is adopted; under P3a with the Σ corner the Higgs direction is stable and
it is met automatically.

## 6. Gate summary (Codex's five items) and verdicts

| Item | P3a + SM to M_Pl | P3a + GRUT reference | P3a + GRUT Σ corner | P3d (MPP) | P3v (volume) |
|---|---|---|---|---|---|
| 1 Postulate | §1 | same + GRUT fields with reference couplings | same + λ_Σ ≲ 3e-5, λ_ΣY ≳ 0.05, λ_HΣ in band | equality | final = ours by volume |
| 2 Compatibility | yes (vacua §1–3) | yes | yes at tree level; one-loop tuned | yes, assumption not derived | yes |
| 3 Forbidden baseline | metastable SM (e.g. PDG central) | same; GRUT has no internal escape | the reference itself (λ_Σ = 1) **and** the λ_ΣY = 0 hidden sector | everything off the line | none |
| 4 Discriminator | z = (M_t − M_t,crit)/σ, frozen line §3.3 | same | none measurable for M_Σ ≥ 1e5; sin θ·M_S band only if light | same z, two-sided | — |
| 5 Supplied / failure | §7 / z > 5 with clean short-distance m_t | same | heavy-window re-solve; failure: no consistent Y/Σ/Higgs vacuum ordering with Y_Q | z ≠ 0 at 5σ; FN μ* = M_Pl form: M_h,crit = 128.6–132 GeV vs 125.13 (vacua §5.4 [Q]) | — |
| Verdict | **MIXED** (testable restriction; standard relation; 3.0–3.5σ / 1.6σ retrospective) | **INCONSISTENT WITH GRUT BRANCH at reference** (inherits the tension) | **UNTESTABLE** (tuned, invisible) | **MIXED** (same tension; FN double form already off) | **STANDARD** (no content) |

## 7. Credit, novelty and what remains supplied

**Credit.** Relation: Degrassi et al. 1205.6497, Buttazzo et al. 1307.3536, Bednyakov et al. 1507.08833, Hiller et
al. 2401.08811 (stability bound); Froggatt–Nielsen hep-ph/9511371, Shaposhnikov–Wetterich 0912.0208, Klinkhamer
1302.1496 (criticality/MPP). Reading: KV 1102.3152 §6 and Arkani-Hamed–Hall–Kolda–Murayama astro-ph/0005111 (global
minimum); Kawai–Okada 1110.2303 and Hamada–Kawai–Kawana 1509.05955 (final-state condition → MPP). Escape: EEGLS
1203.0237, Lebedev 1203.0156. Code: SMDR (1907.02500); the vacua/v_vacua pipeline.

**Novelty.** Recomputation and freezing at PDG 2026; the observation that P3a constrains the GRUT hidden sector
(λ_Σλ_Y ≳ 9e-4 or λ_ΣY ≳ μ²/M_Y²) and that GRUT's only stabiliser candidate conflicts with that constraint unless
λ_ΣY is added; the fixed-vev window for Σ; the ξ_H obligation of the heavy window. These are bookkeeping on known
physics. **Nothing is GRUT-native**: the q-sector, P_q, the reading and the relation are all published.

**Still supplied.** (1) the q-theory action (ε(q), K(q), a 4-form) and P_q itself, a final-state condition the
dynamics does not enforce; (2) the reading of "final vacuum" (global vs worldline vs volume) — closure does not pick it;
(3) SM validity to M_Pl, or a specified UV completion with **no deeper vacua** (generic landscapes violate this [A]);
(4) the residual ρ_Λ (QV1-cal offset, treated as separate); (5) the top-mass interpretation; (6) in GRUT, λ_HΣ, λ_ΣY and
ξ_H, none fixed by the reference.

**Is this the best candidate for distinctive closure content?** It is the only route in the program that turns an
added closure-motivated postulate into a strict restriction on measured quantities with a pre-registered, decisive
held-out test (an e+e− threshold top mass plus α_s). In that sense it is the best *testable* candidate. It is not
distinctive of GRUT: P3a + SM is the standard absolute-stability bound with a mandate attached, it is currently on
the disfavoured side (3.0–3.5σ with the MC mass, 1.6σ with the cross-section mass), and GRUT's own fields offer no
escape at their reference values. If FCC-ee confirms metastability, P3a survives in GRUT only through the tuned Σ
corner, whose consequences are unmeasurable.

## 8. Precision and limits

- SM line: SMDR matching (2-loop + leading 3-loop EW, 4-loop QCD M_t↔y_t), loop order 5 running; ±0.11 GeV matching
  kept from vacua; MSbar λ criterion ('eff' differs by +0.05–0.06 GeV, vacua/v_vacua).
- Other SMDR inputs (M_Z, G_F, α, light masses) are SMDR's PDG-2025 values; M_W effect < 1e-3 GeV (vacua).
- Threshold scans: SM 3-loop + singlet 1-loop, tree matching, RG-improved tree potential; δ_min to ≈0.002, mass
  edges to a factor ≈2 (×7 between the MSbar and λ_eff criteria). Gauge coupling of Σ (6e = 6e-7) neglected.
- Hidden-sector bound: tree level, aligned phases, X = 0. Gate 4 covers the κ (X–Y) direction separately; mixed
  X–Y–Σ directions with small λ_Σ were not scanned.
- Radiative costs and the inflation statement are estimates. ξ_H bounds and MC–pole numbers are quoted at abstract level.
- Same-session AI work; not outside replication.

## 9. Files

- `code/c1_fit_grid.py` → `results/c1_fit_grid.json` (11 SMDR fits at PDG 2026).
- `code/c2_sm_relation.py` → `results/sm_relation.json` (line, slopes, tension, MC scan, scales, held out),
  `results/frozen_relation.json` (the frozen relation, dated 2026-10-27).
- `code/c3_threshold.py` → `results/threshold.json` (sympy betas; generic EEGLS; GRUT Σ window; φ, X; Y).
- `code/c4_hidden_and_costs.py` → `results/hidden_and_costs.json` (hidden-sector P3a bound, λ_ΣY needs, costs, inflation).
- `code/smdr_common.py`, `code/smdr_driver.c` (v_vacua), `code/sm_rge.py` (vacua): reused unchanged.
- `results/sources.json`: inputs, sources, read levels, search record.
