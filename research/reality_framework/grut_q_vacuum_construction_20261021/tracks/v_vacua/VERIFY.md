# Verification of track `vacua` (independent, adversarial)

21-22 October 2026. Verifier track `v_vacua` of the stage-2 construction pass (`../../BRIEF.md`,
`../../BRIEF_STAGE2.md`). Object of verification: `../vacua/NOTES.md`, `../vacua/code/`, `../vacua/results/` as they
stood on this date (NOTES.md 573 lines, results written 23:15-23:23).

**Labels.** [D] derived here; [C] computed by `code/` (outputs in `results/`); [Q] stated by a source at the read
level given; [E] estimate; [A] argued; [R] recalled. Verdicts: CONFIRMED / CORRECTED (with the corrected value) /
REFUTED / NOT CHECKED.

**Independence.** None of the track's code is used to produce any number here. The checks use different formulations
and tools:
- **q-sector.** The track's code is not read for method. The static de Sitter equations are re-derived from an
  explicitly computed Einstein tensor. The series is solved with generic Taylor data. The junction conditions come
  from the A″-coefficients of G_ab and R for a wall metric. The thick wall uses a sparse finite-difference solve.
- **SM running.** RG running and beta functions come from **SMDR v1.3** (Martin and Robertson, arXiv:1907.02500;
  3-loop SM betas of Bednyakov et al. and Chetyrkin-Zoller, with 4- and 5-loop pieces at higher loop order). It is
  driven by our own C driver (`code/smdr_driver.c`).
- **Matching.** Two matchings are used:
  - (B) the Buttazzo et al. 2013 fit formulas, re-transcribed from the TeX;
  - (S) SMDR's own full fit from pole masses: 2-loop plus leading 3-loop electroweak, and 4-loop QCD in M_t ↔ y_t.
- **Audit-only import.** The track's `sm_rge.py` is imported only in `v3b` part A, to compare its beta functions with
  SMDR's term by term. Bytecode writing is disabled, so nothing is written into `../vacua/`.

SMDR, arXiv TeX, PDG files and abstracts are in the scratchpad
(`.../qtheory/v_vacua_dl/`, `.../qtheory/vacua_dl/` read-only). None are in the repository.

**Scope.** No new theory, no empirical validation and no confirmatory test is claimed. Every comparison with
measured (M_t, M_h, α_s) is **retrospective**. "Held out" is used only for FCC-ee and HL-LHC data that do not exist
yet.

---------------------------------------------------------------------------------------------------------------------

## 0. Summary

The track's load-bearing structure survives:
- the action-level formula Λ_i = F(V_i − ε*(μ));
- the junction behaviour;
- the nucleation and crunch numbers;
- the fixed-point reading of P_q;
- the stability-inequality form of the SM relation, including its numbers as a reproduction of the Buttazzo-matching
  pipeline (170.938 against 170.937 GeV).

**No item is REFUTED.**

**Main correction.** The relation's central value is biased by the 2013 matching it uses. With current matching (SMDR,
which reproduces Hiller et al. 2024 to 0.01 GeV) the stability/criticality line is

    M_t ≤ 171.03 + 0.494 (M_h − 125.20) + 0.355 (α_s − 0.1180)/0.0009 GeV   (±0.11 th)        [C]

against the track's 170.94.

**Retrospective tension.**

| M_t input | z (verifier) | z (track) |
|---|---|---|
| direct M_t | 3.15 | 3.34 |
| with the 0.3 GeV MC-to-pole error | 2.68 | 2.84 |
| cross-section pole mass | 1.73 | 1.84 |

**Held-out reach.** FCC-ee 11.6σ (track 12.2σ); HL-LHC 5.25σ (track 5.6σ).

The qualitative verdict is unchanged: (a)/(c)-worldline with the SM to M_Pl is in retrospective tension, and that is
not a decisive failure.

**Other corrections are minor:**
- the m_q w = 30 tracking deviation;
- the "< 0.5 %" bound;
- the λ = 0 crossing scale;
- M_h,crit;
- missing credit.

---------------------------------------------------------------------------------------------------------------------

## 1. Verdict table

Machine-readable: `results/verdicts.json` (built by `code/v4_verdicts.py` from the result files).

| # | Claim (track §) | Track | Verifier | Verdict |
|---|---|---|---|---|
| 1 | Static dS equations E1: ε′ − 2K′Λ/K = μ, E2: Λ = V + ε − μq (1.1) | as stated | Einstein tensor recomputed, G = −3H²g, R = 12H² | CONFIRMED |
| 2 | a2 = 2K′²/(K²ε″) = 4/(3K m_q²); a3 = 8K′²(3K″ε″ − K′ε‴)/(3K³ε″³); dq/dx = 2K′/(Kε″); dK/K = a2 x (1.3) | as stated | identical (sympy, generic Taylor data) | CONFIRMED |
| 3 | dΛ/dV = K²D/(K²D − 4ΛK′²) (1.3) | as stated | equal **on shell** (E1 used); pole at (3/8)K m_q² | CONFIRMED |
| 4 | KV-type ε, static roots (1.3) | 0.010190 (K′=1, x=0.01); 0.0100177 (K′=0.3) | 0.010190; 0.0100177 | CONFIRMED |
| 5 | "correction < 0.5 % down to \|V\| ~ 1e71" (0, 1.3) | < 0.5 % | a2x = 4.27e-3 at 7.5e70, **5.7e-3 at 1e71**; 0.5 % at 8.8e70 | CORRECTED |
| 6 | q-coupled matter: Λ_i = V_i − ε*(μ − θ_i); ΔΛ = ΔV + q_μΔθ (1.4) | as stated | sympy residual 0 | CONFIRMED |
| 7 | Light-q thin wall: [Kx] = 0, [n·∂K] = S/3 = −σ, −K[Kx_ab] = S_ab − h_ab S/3 (1.5) | as stated | δ-parts of G_ab (2[A′/A]h_ab) and R (−6[A′/A]) recomputed | CONFIRMED |
| 8 | Linear wall: kink −0.994σ (m w = 0.01), −0.997σ (0.1); tracking deviation 1.8 % (m w = 10), 0.12 % (m w = 30) (1.5) | as stated | −0.9949, −0.9974; 1.86 %; **0.220 %** (analytic 2/(m w)² = 0.222 %) | CORRECTED |
| 9 | Past light-cone volumes: π t⁴/3; MD 81π/1485 t⁴; ΛCDM 0.1230/H0⁴, η0 = 3.18/H0, t0 = 13.79 Gyr (2.1) | as stated | 3π/55 = 81π/1485; 0.1230; 3.181; 13.79 | CONFIRMED |
| 10 | Volume criterion Γ/H⁴ < 9/(4π); late dS slope (4π/3)/H³; Buttazzo's H t − 11/6 (2.2) | as stated | 9/(4π) = 0.7162; −11/6 recovered | CONFIRMED |
| 11 | Rates (1e-577 Gyr⁻¹Gpc⁻³): log Γ = −741.3, log P = −574.8, log Γ/H_L⁴ = −573.6; 1e574 Gyr (dS), 1e145 Gyr (Minkowski) (2.3) | as stated | −741.26, −574.80, −573.57; 1e574.2, 1e144.6 | CONFIRMED |
| 12 | AdS crunch t_c = π√(3K/\|Λ\|): 3e-41 s, 2e-38 s, 55 Gyr (2.4) | as stated | 3.18e-41 s, 1.95e-38 s, 55.1 Gyr; a = l sin(t/l) solves open FRW | CONFIRMED |
| 13 | CDL quenching moves the boundary by < 1e-3 GeV (2.5) | ε/thr = 5934 at +0.001 GeV | 5947; ε = thr at δM_t ≈ 3e-6 GeV | CONFIRMED (ratios at ≥ 0.01 GeV not meaningful, §4.6) |
| 14 | Track's 1-3 loop β functions and γ_h (5.2) | own transcription | = SMDR to 1e-13 (1-2 loop), ≤ 2e-5 relative (3 loop, rounded coefficients), 4-loop QCD 1.5e-6 | CONFIRMED |
| 15 | At 1307.3536 inputs, λ(M_Pl) = −0.01437; g1, g2, g3, y_t = 0.61537, 0.50546, 0.48731, 0.38248 (5.3) | as stated | −0.01435; 0.61537, 0.50546, 0.48732, 0.38246 | CONFIRMED |
| 16 | **M_t,crit = 170.94 GeV at (125.20, 0.1180), μ* = 4.7e17 GeV (5.4)** | 170.94 | (B) 170.938, μ* 4.73e17; **(S) 171.026 (3-loop), 171.031 (loop order 5)**, μ* 5.1e17 | **CORRECTED → 171.03** |
| 17 | Slopes 0.495 /GeV (M_h), 0.353 per 0.0009 (α_s) (5.4) | as stated | (B) 0.495, 0.354; (S) 0.494, 0.355 | CONFIRMED |
| 18 | Shifts: δλ +0.072, δy_t −0.087, 2-loop − 3-loop −0.109; 4-loop "several times smaller" [E] (5.4) | as stated | +0.072, −0.087, −0.111; loop order 4 (5) − 3: +0.005 (+0.005) | CONFIRMED |
| 19 | Matching error 0.113 GeV covers the matching (5.4, 7) | ±0.0003 (λ), ±0.0005 (y_t) | S − B at M_t: y_t −0.00094, g3 −0.0017, λ +0.00012 → +0.088 GeV | CORRECTED (§4.3) |
| 20 | M_h,crit = 128.54 at M_t = 172.56 (5.4, 5.5) | 128.54 | (B) 128.53; **(S) 128.36** | CORRECTED |
| 21 | FN double condition (M_t, M_h) at μ* = 1.22e19 / 2.435e18 / 1e18 (5.4) | 174.20/132.03; 172.61/128.64; 171.71/126.77 | 174.21/132.03; 172.61/128.65; 171.71/126.78 (B) | CONFIRMED |
| 22 | V_eff(2.435e18) = −7.5e70 GeV⁴ = 3e117 ρ_Λ; λ = 0 at 4.2e10 GeV (5.4, 3.1) | as stated | −7.55e70 (1-loop λ_eff), 3.0e117; **λ = 0 at 4.7e10 GeV** (6.2e9, not 5.7e9, at the 1307.3536 point) | CORRECTED (crossing only) |
| 23 | z = 3.3 (direct), 2.8 (+0.3 MC), 1.8 (xsec) (5.5) | 3.34, 2.84, 1.84 | (B) 3.33, 2.83, 1.84; **(S) 3.15, 2.68, 1.73** | CORRECTED |
| 24 | Held out: FCC-ee σ = 0.13 GeV, z ≈ 12; HL-LHC z ≈ 5.6 (5.6) | 0.1325, 12.25, 5.57 | (B) 0.1325, 12.24, 5.56; **(S) 11.6, 5.25** | CORRECTED |
| 25 | 'eff' criterion 171.00; eff − msbar = +0.064 (5.3, 5.4) | as stated | 1-loop λ_eff: 170.989; +0.051 | CONFIRMED (2-loop piece not recoded) |
| 26 | Literature cross-checks: Bednyakov 171.44 at (125.7, 0.1185); 1307.3536 formula gives 171.12, not its quoted 171.53 (5.3) | own 171.38; inconsistency recorded | (B) 171.38, (S) 171.47-171.48; formula 171.115 | CONFIRMED |
| 27 | PDG 2025 inputs (5.1) | 125.20±0.11; 172.56±0.31; 172.4±0.7; 0.1180±0.0009; 80.3692±0.0133 | identical in pdg-2025-v0.2.3.sqlite and the QCD review eq. (9.25) | CONFIRMED |
| 28 | Readings (a)-(d), fixed point, Volovik counting (3) | see track §3.5 | re-derived; two interpretive notes (§5) | CONFIRMED |
| 29 | Credit and novelty (4, 6.4) | as listed | five omissions (§6) | CORRECTED (credit only) |
| — | Λ_V = 10^9.497 GeV, Λ_I/Λ_V = 12.7; 2-loop λ_eff values; FCC/HL-LHC precisions; AFS and CMS rates beyond abstract | as stated | — | NOT CHECKED |

---------------------------------------------------------------------------------------------------------------------

## 2. Corrections (every one, with the corrected value)

C1 (row 16). **Stability/criticality line: 170.94 → 171.03 GeV** at (M_h, α_s) = (125.20, 0.1180).
- Running: SMDR loop order 5 (3-loop: 171.026).
- Matching: SMDR full fit (S).
- The track's pipeline is reproduced exactly with its own (2013) matching: 170.938.
- The difference is the matching (C3).
- 171.03 agrees with Hiller et al. 2401.08811 (171.04, same criterion, PDG 2024 = same M_h and α_s). The track lists
  that number as consistent without resolving the 0.1 GeV gap.
- μ* = 5.1e17 GeV (S) against 4.7e17 (B).
- [C]

C2 (row 20). **M_h,crit at M_t = 172.56: 128.54 → 128.36 GeV** (S). At M_t = 172.4 the value is 128.02 (Hiller:
127.97). The α_s needed at M_t = 172.4 is 0.12148 (Hiller: 0.1215). [C]

C3 (row 19). **Matching error budget.** At the PDG point the 2013 fits differ from SMDR's full fit by:
- y_t: −0.00094, i.e. 1.9 times the quoted ±0.0005. The source is the 4-loop QCD term in M_t ↔ m_t; 1307.3536 §3 itself
  notes that its y_t omits it.
- g3(M_t): −0.0017. This term is not in the track's budget at all. Pure-QCD 4-loop running gives 1.16305 and SMDR
  1.16343 (+0.0004 electroweak), against the fit's 1.16516. Bednyakov et al. 1507.08833 likewise have 1.16540, which
  is the pure-QCD value 1.16499 + 0.0004. The 1307.3536 TeX carries a commented-out "1.1645" beside 1.1666
  (`results/v3c_g3_qcd.json`).
- λ: +0.00012.

Effect on M_t,crit, one coupling at a time:

| Coupling | Shift (GeV) |
|---|---|
| y_t | +0.163 |
| g3 | −0.111 |
| λ | +0.029 |
| g, g′ | +0.008 |
| **net** | **+0.088** |

The net shift lies inside the quoted 0.113 GeV only because the two large pieces cancel. The budget should name the
4-loop y_t term and the g3 input. [C]

C4 (rows 23, 24). **Tension and reach** with C1 (same error budget as the track):
- direct M_t: z = 3.15 (p = 8e-4);
- with the criterion and the 0.3 GeV MC-to-pole error: 2.68 (p = 3.7e-3);
- cross-section pole mass: 1.73 (p = 0.042);
- adding |S − B| = 0.088 GeV as a further matching error: 3.10, 2.64, 1.72.

Held out, illustrative only:
- FCC-ee: σ_tot = 0.1326 GeV, z = 11.6 (direct central value), 10.4 (cross-section);
- HL-LHC + lattice: z = 5.25.

The verdict "retrospective tension, not decisive" stands. [C]

C5 (row 8). **Tracking deviation at m_q w = 30: 0.12 % → 0.22 %.** The leading-order deviation for a sech² wall is
ρ″/(m²ρ) = 2/(m w)² = 0.222 %. The track's 0.12 % contradicts its own stated scaling, so it is a resolution artifact.
The values at m_q w = 10 (1.8 %) and the kinks are confirmed. [C]

C6 (row 5). **"< 0.5 % for |V| down to 1e71" → < 0.6 % (0.57 % at 1e71; 0.5 % is reached at 8.8e70 GeV⁴).** This is
trivial; the conclusion is unchanged. [C]

C7 (row 22). **λ = 0 crossing at PDG central values: 4.2e10 → 4.7e10 GeV.** At the 1307.3536 point it is 6.2e9 GeV,
not 5.7e9. The track takes the grid point below the crossing on a 60-point grid. [C]

C8 (row 29). **Credit**: see §6.

---------------------------------------------------------------------------------------------------------------------

## 3. q-sector checks (`code/v1_qsector.py`, `results/v1_qsector.json`)

### 3.1 Static vacua at a shared μ [D][C]

- **de Sitter check.** For the de Sitter metric (flat slicing), sympy gives G_μν = −3H² g_μν and R = 12H². With
  T = −V g, the (g) equation of BRIEF.md gives 3K H² = V + ε − μq ≡ Λ. The (q) equation gives ε′ − 2K′Λ/K = μ. This
  is the track's E1 and E2, signs included.
- **Series.** Expand about q_μ (ε′(q_μ) = μ) with generic Taylor data (K0..K3, e2..e4). Solving order by order gives
  the track's a2 and a3 exactly, together with dq/dx = 2K′/(Kε″) and dK/K = a2 x. These match equations §5.2 for a
  constant shift.
- **Expansion rate.** Expressed with K0, the expansion rate is H² = Λ/(3K(q)) = x/(3K0) + O(x³), with no x² term.
  This reproduces equations §5.2 and is not a contradiction. The track's O(x/(K m_q²)) correction applies to Λ in
  local-G units, and in H² it cancels. "Differences gravitate fully" therefore holds even more closely for the
  expansion rate.
- **Slope.** Implicit differentiation gives dΛ/dV = (K²ε″ − 2ΛKK″ + 2ΛK′²)/(same + 2(μ − ε′)KK′). Using E1 this equals
  the track's K²D/(K²D − 4ΛK′²). The formula holds on the static family only, which is how the track uses it.
- **Numerical roots** (KV ε = (−q² + q⁴/3)/2, K = 1 + K′(q − 1)):

  | K′ | x | Λ (exact) | Three-term series |
  |---|---|---|---|
  | 1 | 0.01 | 0.01019005 | 0.01018933 |
  | 0.3 | 0.01 | 0.01001772 | 0.01001771 |
  | 0 | any | x exactly | x |

### 3.2 Matter coupled to the 3-form, coexistence [D][C]

With the matter term −θ_i q and K′ = 0:
- ε′ + θ_i = μ;
- Λ_i = V_i + ε(q_i) − q_i(μ − θ_i) = V_i − ε*(μ − θ_i) (sympy residual 0);
- q_i − q_j = (θ_j − θ_i)/ε″;
- the first-order difference term is q_μ(θ_i − θ_j).

At θ = 0, P_i − P_j = V_j − V_i for every μ, so there is no Maxwell point. With θ_i ≠ θ_j a coexistence value μ_c
exists, and P = 0 there is one further condition involving θ. The k-form Jacobian argument (identical rows q_a) is
correct. Volovik's counting ν ≤ k was checked against hep-ph/0309144 (TeX, lines 300-314) and the quote is accurate.

### 3.3 Junction conditions [D][C]

Take ds² = dn² + A(n)²η_ab dx^a dx^b, a symmetric pure-tension wall:
- the A″ coefficient in G_ab is 2h_ab/A, i.e. G_ab ⊃ −([Kx_ab] − h_ab[Kx])δ with Kx_ab = (A′/A)h_ab;
- in R it is −6/A, i.e. R ⊃ −2[Kx]δ;
- G_nn has none.

With ∇∇K ⊃ n n[∂_nK]δ, the δ-part of (g) is −K([Kx_ab] − h_ab[Kx]) + h_ab[∂_nK] = S_ab. Since R = 2(ε′ − μ)/K′ has no δ
when q is continuous, [Kx] = 0. Then the trace gives [∂_nK] = S/3 and the traceless part gives −K[Kx_ab] = S_ab − h_ab
S/3, which are the track's statements.

Linear wall, sech² profile, finite-difference solve:
- kink −0.9949σ (m w = 0.01) and −0.9974σ (m w = 0.1);
- local-tracking deviation 42.9 % (m w = 1), 1.86 % (10), **0.220 %** (30) → C5.

### 3.4 Sizes [C]

- Natural m_q = √(2/3) M = 1.99e18 GeV.
- ξ_max = m_h²/M² = 2.6e-33.
- m_h² h*² = 3.5e39 to 9.3e40 GeV⁴ (track ~1e40).
- m_q w = 1.6e16 for w = 1/m_h; m_q R = 19.9 for R = 1/(1e17 GeV).
- δλ tolerance 4ρ_Λ/h*⁴ = 2e-117 (h* = 4.7e17), 3e-120 (2.435e18). The track's "~1e-118" is within that range.

---------------------------------------------------------------------------------------------------------------------

## 4. Nucleation and SM checks

### 4.1 Light-cone volumes, rates, crunches (`code/v2_nucleation.py`) [D][C]

- **Analytic volumes.** Minkowski V4 = πt⁴/3. Matter domination gives V4/t⁴ = 3π/55 = 81π/1485 = 0.17136 and
  V4 H⁴ = 16π/1485. De Sitter from a = 1 gives (3H⁴/4π)V4 − Ht → −11/6 (Buttazzo's form).
- **ΛCDM** (Ω_m = 0.315, Ω_r = 9e-5, H0 = 67.4), by conformal-time quadrature:
  - V4 = 0.12299/H0⁴, η0 = 3.181/H0, t0 = 13.79 Gyr;
  - Buttazzo's "≈ 0.15" corresponds to Ω_m ≈ 0.27 (0.143);
  - the late slope approaches 4π/(3H_L³) (7.377 at a = 1e3 against 7.390).
- **Rates.** 1 Gyr⁻¹Gpc⁻³ = 5.455e-165 GeV⁴. Every row of the track's table is reproduced to 0.01 in the log. Rates
  were checked against the CMS 1707.09301 abstract (1e-577; 1e-295 to 1e-1465) and the AFS 1707.08124 abstract (1e139
  yr; > 1e58 yr at 95 %).
- **Crunch.** The open-slicing a = l sin(t/l) satisfies H² = 1/a² − 1/l² exactly (sympy). The crunch times are
  reproduced.

### 4.2 Audit of the track's RGEs against SMDR (`v3b` part A) [C]

At six coupling points the loop-by-loop pieces were compared, for d(g_i²)/d ln μ², d(y_t²)/d ln μ², dλ/d ln μ² and γ_h
(track) against d ln v/d ln Q (SMDR, Landau gauge):
- 1-2 loop: agreement to ≤ 1.4e-13 relative;
- 3 loop: ≤ 1.2e-11 for the gauge couplings, ≤ 3e-6 for y_t, ≤ 1.8e-5 for λ, ≤ 8e-5 for γ_h. These are rounding of
  the 1307.3536 App. A numerical coefficients;
- the 4-loop pure-QCD g3 term (−2472.28 G3⁵/(16π²)⁴) agrees with SMDR to 1.5e-6.

The transcription is correct.

### 4.3 The SM line, two matchings (`v3b` parts B, C, F; `v3a`; `v3c`) [C]

| Quantity | Track | (B) 2013 fits + SMDR running | (S) SMDR fit + SMDR running | Literature |
|---|---|---|---|---|
| λ(M_Pl) at 1307.3536 inputs | −0.01437 | −0.01435 (loop order 5: −0.01433) | — | −0.0143 |
| M_t,crit (125.15, 0.1184) | 171.07 | 171.07 | 171.16 | 171.12 (formula), 171.53 (quoted) |
| M_t,crit (125.7, 0.1185) | 171.38 | 171.38 | 171.47 (171.48, loop order 5) | 171.44 (Bednyakov) |
| **M_t,crit (125.20, 0.1180)** | **170.94** | **170.938** | **171.026 / 171.031** | **171.04** (Hiller, λ > 0), 171.10 (eff) |
| μ* (GeV) | 4.7e17 | 4.73e17 | 5.06e17 | — |
| 2-loop − 3-loop RGE | −0.109 | −0.111 | −0.109 | — |
| M_h,crit at 172.56 | 128.54 | 128.53 | 128.36 | — |
| α_s for stability at M_t = 172.4 | — | — | 0.12148 | 0.1215 (Hiller) |

**Interpolation.** The SMDR couplings are interpolated in M_t (quadratic, 6 fits) and linearly in M_h and α_s. The
residuals at the 16 fitted points are ≤ 1e-5 within ±0.6 GeV of M_h, and up to 4e-4 in λ at M_h = 132 (not used for
the line).

The criterion is min_{μ ≤ M_Pl} λ(μ) = 0 (as the track's 'msbar'). The minimum is found by a parabola on a 241-point
grid in ln Q from 1e8 GeV to M_Pl. The cut at the reduced Planck mass gives the same result, since μ* < 2.4e18.

### 4.4 Tension and held-out (`v3b` part D) [C]

The error budget is the track's: σ(M_t); the propagated σ(M_h) and σ(α_s); 0.113 matching; 0.032 criterion; and
optionally 0.3 MC.
- (B) reproduces the track: 3.33 / 2.83 / 1.84, FCC-ee 12.24, HL-LHC 5.56.
- (S) gives C4.

The held-out scenarios are correctly labelled as data that do not exist yet. The FCC-ee and HL-LHC precision inputs
themselves are NOT CHECKED (abstract-level and [R] in the track).

### 4.5 Potential depth (`v3b` part B) [C]

At PDG central values (M_t = 172.56, (B) matching), with my own 1-loop Landau-gauge λ_eff and e^{4Γ} from SMDR's v(Q):

| Scale | λ | λ_eff (1-loop) | V_eff (GeV⁴) |
|---|---|---|---|
| 1e17 GeV | −0.01046 | −0.00983 | −2.05e65 |
| 2.435e18 GeV | −0.01041 | −0.01007 | −7.55e70 = 3.0e117 ρ_Λ |

These agree with the track's 2-loop values (−2.03e65, −7.48e70) to 1 %. The track's 1-loop λ_eff formula and sign
convention for Γ were checked against the Coleman-Weinberg form (Sher's ξ = e^{−∫γ}).

### 4.6 CDL quench (`v3b` part E) [E]

With the track's definition (σ = ∫√(2V) over V > 0 below h_min), ε/threshold is:

| δM_t above M_t,crit (GeV) | ε/threshold |
|---|---|
| 1e-6 | 0.38 |
| 1e-5 | 4.3 |
| 1e-4 | 64 |
| 1e-3 | 5947 (track: 5934) |

ε reaches the threshold at δM_t ≈ 3e-6 GeV, so "gravity moves the boundary by < 1e-3 GeV" is confirmed with a wide
margin.

Caveat: for δM_t ≥ 0.01 GeV the true minimum moves to ≥ 2e18 GeV and the thin-wall picture fails. A full-path tension
∫√(2|V|) gives ε/threshold = 4 at 0.01 GeV and 0.17 at 0.05 GeV, against the track's 7.6e8 and 2.9e17. The track's
large-δ ratios are therefore not meaningful, but no conclusion uses them.

---------------------------------------------------------------------------------------------------------------------

## 5. Readings of P_q, logic [D][A]

Re-derived:
- (a) gives Λ_ours = V_ours − V_true ≥ 0.
- (c) with a worldline measure has the unique fixed point f = true when Γ > 0. With f = ours, true-vacuum interiors
  are AdS and crunch, every worldline (V4 → ∞) ends there, and no T → 0 final state exists. So (c) reduces to (a).
- (c) with a volume measure in de Sitter keeps f = ours when Γ/H⁴ < 9/(4π).
- (d) needs V_i = V_j, and nothing in the q-sector supplies it, so it is FN's principle added.
- Shared μ is automatic ([μ] = 0 across matter walls).

None of these is contradicted by the stage-1 records (equations §5, ledger §1.3-1.4).

Two notes, neither a correction of a conclusion:
1. Reading (b) "self-consistent only if Γ = 0" imposes the end-point condition of (c) on (b). It is an interpretive
   step, [A], and the track's [D] label on §3.2 overstates it.
2. Under (a) or (c)-worldline with **uncalibrated** QV1, reproducing ρ_Λ,obs > 0 is possible in only two ways. Either
   our vacuum is the global minimum, so Λ_ours = 0 and QV1 has its known residual problem. Or the deeper vacuum lies
   exactly ρ_Λ below ours, which is the equality M_t = M_t,crit to ~1e-118 in δλ, i.e. AHHKM's tuned splitting. Both
   satisfy the track's inequality, so its relation stands.

The track states this as "uses only the smallness of ρ_Λ". Making the equality branch explicit would help the gate
assembly.

---------------------------------------------------------------------------------------------------------------------

## 6. Credit and novelty

The credits the track gives check out at the stated read levels:
- FN hep-ph/9511371;
- Klinkhamer 1302.1496 (FN-type MPP with gravity, not q-theory, as the track says);
- Volovik hep-ph/0309144 (quotes verified in the TeX);
- the Buttazzo stability formula and its 2.8σ;
- Bednyakov's table (g_s = 1.16540, y_t = 0.93517, λ = 0.12714);
- the Hiller et al. table and their "< 2σ / 3σ" statements.

Missing credit, verified at abstract level (scratchpad `v_vacua_dl/arxiv/`):

| Work | Why it belongs |
|---|---|
| Shaposhnikov-Wetterich, arXiv:0912.0208 (PLB 683 (2010) 196) | λ(M_Pl) ≈ 0 with a fixed point at zero gives m_H ≈ 126 GeV. This is the best-known criticality prediction and belongs beside FN for the equality form (d). |
| Holthausen-Lim-Lindner, arXiv:1112.2415 | Planck-scale boundary conditions (λ = 0, β_λ = 0, ...) give Higgs-mass predictions. |
| Bezrukov-Kalmykov-Kniehl-Shaposhnikov, arXiv:1205.2893 | NNLO-level stability bound, contemporaneous with Degrassi et al. 1205.6497. |
| Kawai-Okada, arXiv:1110.2303; Hamada-Kawai-Kawana, arXiv:1509.05955 | Coupling constants, including Λ, fixed by a condition on the far-future state of the universe, which yields MPP-type conditions. This is the closest prior art to the track's chain "a final-state condition fixes the vacuum-energy reference, so a stability/degeneracy condition follows". The mechanism differs (wormholes or a multi-local action, not a 4-form), so the track's narrower "not found in searched papers" claim can stand with this citation. |
| Martin-Robertson, arXiv:1907.02500 (SMDR); Alam et al. (2022) | The matching behind Hiller's 171.04, and the source of C1-C3. |

Novelty verdict: the track's own description, "bookkeeping on known physics", is accurate. Nothing is GRUT-native
(agreed).

---------------------------------------------------------------------------------------------------------------------

## 7. Inputs and provenance

- M_h = 125.20 ± 0.11 (S126M), M_t = 172.56 ± 0.31 (Q007TP, S = 1.6), 172.4 +0.68 −0.64 (Q007TP4) and M_W = 80.3692 ±
  0.0133 (S043M). All were read from `pdg-2025-v0.2.3.sqlite` (track's download, scratchpad).
- α_s = 0.1180 ± 0.0009 (PDG 2025 QCD review eq. (9.25), text grep).
- SMDR's other inputs are the PDG 2025 reference values (`smdr_pdg_2025.h`): G_F, M_Z = 91.1876, α, Δα_had = 0.02783,
  light masses. With these, SMDR's M_W output is 80.350 (PDG convention), which affects the line by < 1e-3 GeV
  (track).
- Buttazzo fit formulas: 1307.3536 TeX `artv2.tex` lines 1183-1225 (re-transcribed).
- Cosmology: H0 = 67.4, Ω_m = 0.315, Ω_r = 9e-5 [R], as in the track. ρ_Λ and M are from BRIEF.md.
- Calibration: none. The SM line uses no observed vacuum-energy value. QV1-cal appears only in readings where the
  track labels it.

## 8. Limits

- Matching (S) is SMDR's: 2-loop plus leading 3-loop electroweak/QCD and 4-loop QCD for the top mass. The quoted
  ±0.11 th is kept from the track; I did not re-estimate it. The 0.088 GeV matching shift measures what the 2013
  error bars missed.
- The 'eff' criterion was checked at 1 loop only. Λ_V and Λ_I were not checked.
- The CDL estimate is thin-wall (§4.6).
- The interpretation of the MC mass (±0.3 GeV) is taken from the track and not recomputed.

## 9. Files and reproduction

- `code/smdr_driver.c`: own C driver for SMDR (modes: betas, run, fit, fitrun).
- `code/v3_build_smdr.sh`: downloads and builds SMDR in the scratchpad.
- `code/smdr_common.py`: wrapper.
- `code/v1_qsector.py`: §3.
- `code/v2_nucleation.py`: §4.1.
- `code/v3a_smdr_fit_grid.py`: 16 SMDR fits, cached in `results/v3a_smdr_fit_grid.json`.
- `code/v3b_sm_check.py`: parts A-F.
- `code/v3c_g3_qcd.py`: pure-QCD g3.
- `code/v4_verdicts.py`: `results/verdicts.json`.
- `code/run_all.py`.

Reproduce (about 4 minutes; numpy 2.x, scipy, sympy, gcc):

    cd code && sh v3_build_smdr.sh && python3 -B run_all.py
