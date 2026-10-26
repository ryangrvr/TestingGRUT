# Track rp_portals: portal sources, matching checks and cold leakage in Codex e160356 (handoff items 3–4)

26 October 2026 · Claude, issue #13 · same-session AI work, not outside replication. Packet read only through
`git archive e160356 research/reality_framework/grut_reheating_preparation_20261026` into the scratchpad
(CLAUDE_HANDOFF.md, REPORT.md, MATHEMATICS.md, BOUNDARY.md, PORTALS.md, PORTAL_RESULTS.json, SOURCES.json, code/portals.py
read; **no Codex code was run**). Parent inputs read at 04c53f5 (heavy_vector_budget INPUTS/BUDGET/BACKGROUND,
gauge_point_test PAIR_RESULTS, heavy_scaling and kinetic_contract MATHEMATICS, the 8b0d96b / fa3bc97 portal reports).
Nothing was committed; nothing was written outside this directory and the scratchpad. The kinetic tracks (kin_*) were
not touched; nothing here is a kinetic or f_sym statement.

**Scope.** One supplied tree-level working model at the heavy spot point (H_i = 1.3334962e11 GeV, M_Y = 30 H_i, h = 40 H_i,
μ = 7.15945616 H_i, κ = 3.458684721e-3, r = 0.45, e = 1e-7, c_H = 10; light sector g_X = 0.1, λ_X = 1e-3, λ_φ = 0.03,
v = 0.924089637 GeV, Q = 200000, Y_Q = 4.8294e-9). Raw sources only: no f_sym, no plasma kinetics, no observational
statement, no validation or exclusion. Nothing here is evidence that X, Y, Σ, φ or the vector exist.
Labels: [D] derived, [C] computed here, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

**Independence.** Every [C] number is from this track's own code in `code/` (own enumeration of the SM hypercharge
content, own moment solver written in a different variable/integrator, an own *momentum-resolved* solver that Codex does
not have, own reheating-era solver, own Higgs freeze-in quadrature with a g*(T) sensitivity, sympy derivation of β_b).
Reuse of our earlier tracks, stated: the unit-coefficient X⁶ rates of `tracks/xcharge` (forced.py, pairs.py) and the
`tracks/hw_gates` gate-5 survival numbers are read as data for the leakage check; the case-L threshold variants and the
existing debris numbers (R_med = 1.1437, R_A = 48.47 own modes / 64.0 GMR) are quoted from `tracks/hw_critic` and
`tracks/hw_gates` [Q own]. Codex's E_X = 38386.25 H_i⁴ and R_0 = 40.5626 (the finite-band primordial injection) are [Q]
normalisations, not recomputed here (our own band gives 48.5; both are carried).

## 0. Verdict table

| # | Codex claim (e160356) | Codex [Q] | This track | Verdict |
|---|---|---|---|---|
| 1 | Γ_SM = 41 ε² g_Y² m_A/(96π) (ΣY² Weyl = 10, Higgs 1/2); Γ_X = e² m_A β³/(48π); B_X at ε_BX = 1e-7 | 5.9562e-10, 2.12233e-10 GeV; B = 0.262712 | explicit field-by-field sum gives the same coefficient (ratio 1.000000); same widths; B = 0.262712 [D, C]. Sensitivity: one-loop-run g_Y(m_A) = 0.380 gives B = 0.252 | CONFIRMED |
| 2 | Thermal pair yield Y_pair = 1.663e-7 (MB closure); exact identity Y_pair = C(3π/2)(Γ_X/H_m)(1−B) | 1.66294e-7 | 1.66294e-7 from both solvers; identity holds to 6e-10 (moment) and ≤ 4e-4 (momentum-resolved, ε ≥ 1e-10); equals the HJMRW decay freeze-in formula with g_A = 3, Γ = Γ_SM, times B [C, D, R] | CONFIRMED |
| 3 | Hidden energy injection R_thermal ≈ 9.84e5, R_total = 983809.9 at ε_BX = 1e-7 | 983804.4 + 5.46 | moment solver (same closure, own code): 983804.44; **momentum-resolved: 968381** (R_total 968387), −1.57 %; the MB time profile overstates the delay factor (34.05 → 33.52) because the free-streaming A keep p/T ≈ 3 while the closure assigns them the thermal p/T ≈ √(2x) [C] | CORRECTED (small): R_total = 9.68e5 |
| 4 | Raw bound R_total ≥ 40.5626 B^{3/2} + 39184 (1−B) ≥ 40.5626; derivation from branching and free on-shell redshift (J ≥ 5 Y_pair, R_prim = R_0 B^{3/2}) | as stated | re-derived [D]; R_sat = 39183.974 own evaluation; R_sat/(1.5 R_0) = 644 so f(B) decreases, min = R_0 at B = 1; J/J_prompt ≥ 1.00002 at every ε in the momentum-resolved runs; Bessel integrals 3π/2, 15π/2 to 2e-12 [C] | CONFIRMED |
| 5 | Source is a hot-RD statement; low T_R or other vector mass needs separate treatment | caveat | T_R scan (decaying inflaton, own solver): ratio to hot RD = 1.00 (T_R ≥ 1e7), 0.93 (4e6), 0.26 (1e6), 0.010 (4e5), 1.4e-6 (1e5) GeV. The visible outlet lowers the raw injection below R_0 only if T_R ≲ 2–3e5 GeV (Γ_I ≲ 1e-18 H_i). Every history on the table is hot: prompt T_R = 3.08e14, finite release 1e14–1e15, Codex's six reheating runs 5.3e13–1.7e14 GeV (x_R ~ 1e-8): ratio = 1 [C] | CONFIRMED and quantified: the caveat is inactive for the heavy window |
| 6 | (vector mass) | — | at fixed ε_BX = 1e-7: R_thermal = 9.8e5, 2.2e5, 4.6e4, 7.8e3 for e = 1e-7, 1e-6, 3e-6, 1e-5 (band R_0 ∝ 1/e: 40.6 → 0.41); at fixed ε_BX/e = 1: R_thermal ∝ e, 5.8e6 at e = 3e-6 [C] | NEW |
| 7 | Higgs-portal source b_Hφ = 5e-11: Γ_hss = 4.81568e-20 GeV, Y_pair = 2.237e-8, f_rest = 9.39 | as stated | 4.81568e-20; Y_pair = 2.23656e-8 (closed form and quadrature agree to 6e-7); f_rest = 9.394. With g*(T) (95.5 at m_h/3) instead of 106.75: Y_pair = 2.65e-8, f_rest = 11.1 (+18 %). Added: prompt conformal energy J_h = 5 Y_hss = 5.59e-8 ⇒ R_h = 1943 (2450 with g*(T)) in the spot normalisation; Γ_h/H = 1.9e11 (prompt); s → X X̄ open, β = 0.577, Γ = 1.7e-5 GeV [C] | CONFIRMED; +18 % g* sensitivity |
| 8 | θ/b = 0.0291; ATLAS one-channel b_max = 5.0501e-3, λ_HX,max = 7.1419e-3 | as stated | 0.0290978; 5.05009e-3; 7.14189e-3 [C] | CONFIRMED |
| 9 | One-loop mixed gauge-charge trace = 0; ε_BX allowed, not forced | 0 | Σ_i q_i Y_i = 0 field by field [C]. Sharper [A]: the two coupling phases (κ, μ) are removable (3 fields, 1 gauge direction), so hidden charge conjugation is an exact symmetry of the stated tree Lagrangian including λ_HX|X|²|H|²; B_μν F_X^μν is odd under it ⇒ ε_BX is not generated at **any** loop order from these couplings; only a CP-violating holder plus a hidden–SM bridge could feed it | CONFIRMED and sharpened |
| 10 | Isolated β_b = g_X² λ_HX/(8π²) | 1.2665e-4 λ_HX | sympy Coleman–Weinberg log-derivative with Tr M_X⁴ = 2(m₀² + g_X²φ² + λ_HX|H|²)²: identical [D]; 1.26651e-4 λ_HX | CONFIRMED |
| 11 | Cold leakage with the full heavy denominator: α₆ = 2.10886e-31 GeV⁻², 1/(1−r²) = 1.25392, Γt₀ = 7.08368e-17, pair/coherent 2653 | as stated | α₆ 2.10886e-31, c6 = 1.25039e6, ratio 1.25392; Γt₀ = 7.0837e-17 from our xcharge coherent (3.9196e-50 /s) + free-Born pair (1.03997e-46 /s) unit rates; LDA variants 6.90–7.15e-17; hw_gates' own earlier value 7.0835e-17; leading denominator would give 4.51e-17 [C, Q own] | CONFIRMED |
| 12 | Heavy-photon contact: σ_Xe(0) = 1.855e-51 cm² at ε_γ = 1, 5.05e-13 events/t·yr | as stated | 1.8550e-51 cm², 5.046e-13 [C, arithmetic] | CONFIRMED (not asked; side check) |
| 13 | Six reheating runs, preparation state, Weyl/clock statements (items 1–2) | — | outside this track | NOT CHECKED |

No numerical error was found in Codex's portal ledger beyond the 1–2 % closure effect in item 3.

## 1. Widths, branching, trace, β_b (`code/widths.py`)

Above electroweak breaking the massive A acquires the coupling ε g_Y Y to the hypercharge current (field shift B → B − εA
at O(ε)) [D]. Enumerating Q_L (6 Weyl, Y = 1/6), u_R (3, 2/3), d_R (3, −1/3), L (2, −1/2), e_R (1, −1) per generation gives
ΣY² = 10/3 per generation, 10 for three; the Higgs doublet gives 1/2. With Γ(V → ψψ̄) = g²m_V/(24π) per Weyl fermion and
Γ(V → φφ*) = g²m_V/(48π) per complex scalar [R], Γ_SM = ε²g_Y²m_A[10/(24π) + 1/(96π)] = 41ε²g_Y²m_A/(96π) exactly. With
g_Y = 0.37, m_A = 6eh = 3.2003909e6 GeV: Γ_SM(1e-7) = 5.9562e-10 GeV, Γ_X = 2.12233e-10 GeV (β³ = 1 − 5e-9), B = 0.262712.
Γ_SM/Γ_X = 20.5 (ε g_Y/e)² = 2.806 (ε/e)², so B depends only on ε_BX/e. Sensitivity: g_Y run at one loop to m_A is 0.380
(b = 41/6), giving Γ_SM 5.6 % higher and B = 0.252.

Mixed trace: no field carries both charges, Σ q_i Y_i = 0 [C]. β_b: V₁ = (1/64π²)·2M_X⁴(ln M_X²/μ² − 3/2) with
M_X² = m₀² + g_X²φ² + λ_HX H†H; μ∂V₁/∂μ = −(1/32π²)·2M_X⁴, whose φ²H†H coefficient is 4g_X²λ_HX/(32π²) = g_X²λ_HX/(8π²)
(sympy, `beta_b_check = true`) [D]. This is one additive contribution; the b-proportional pieces (b·λ_H, b·λ_φ, b·g_X²…) are
multiplicative, so b = λ_HX = 0 at a matching scale is radiatively consistent with the SM-blind hidden sector [A].

## 2. Thermal vector source (`code/freezein.py`)

**Model (Codex's, re-derived).** x = m_A/T, H = H_m/x², H_m = √(π²g*/90) m_A²/M_r = 1.43919e-5 GeV, Y_eq = C x²K₂(x),
C = 45·3/(4π⁴g*s) = 3.24569e-3. Inverse-decay gain Γ_SM⟨m/E⟩n_eq = Γ_SM m_A²T K₁/(2π²) [D, exact for MB]; loss Γ_tot⟨m/E⟩;
energy to daughters Γ_X m_A n_A (exact: the rate Γm/E times the energy E) [D]; J = ρ_X/(sT) = 4ρ_X/(3ρ_r) at constant g*;
R = F_R J with F_R = 4(M_r/H_i)²/E_X = 3.47455e10 [Q normalisation], i.e. R = ρ_X a⁴/(E_X H_i⁴) with a = 1 at the release.

**Solver (1), moment, MB closure, own formulation** (LSODA/BDF in ln x; controls rtol 1e-9/1e-11, x₀ 1e-4/1e-5):
R_thermal = 983804.44 at ε = 1e-7, the 19-point Codex grid reproduced to ≤ 1e-7 relative at every ε (this is an
independent code solving the same closure, so it is a reproduction of the closure, not of Codex's program).

**Solver (2), momentum-resolved, no closure.** Per comoving momentum q = p/T (constant in RD at constant g*),
∂f/∂ln x = (x²/H_m)(x/√(q²+x²))[Γ_SM e^{−√(q²+x²)} − Γ_tot f]; Y_A = C∫q²f dq (reduces to Cx²K₂ for a thermal f, checked to
1e-14); 240–480 log bins, q ∈ [0.01, 60–120], LSODA rtol 1e-8/1e-10: R_thermal = 968381.24 ± 0.01 at ε = 1e-7.
Ratio momentum/moment of J: 0.976 (1e-11, number tolerance-limited), 0.991 (1e-10 … 1e-8), 0.984 (1e-7), 0.965 (3e-7),
0.925 (1e-6), 0.951 (1e-5), 0.998 (1e-4), 1.000 (1e-3). Interpretation: the closure assigns the A population the thermal
⟨m/E⟩ = K₁(x)/K₂(x) at the current x, i.e. p ≈ √(2m_AT); the actual free-streaming A produced at x ≈ 2–3 keep p/T ≈ 3, so
they are colder, decay slightly earlier and deposit slightly less comoving energy. The effect is largest where the decay
epoch x_d = √(H_m/Γ_tot) is smallest (ε ≈ 1e-6: x_d ≈ 20). The pair number is closure-independent (identity).

| ε_BX | B | Y_pair | R_thermal (closure / momentum) | R_total (momentum) | raw bound | f_rest |
|---|---|---|---|---|---|---|
| 1e-10 | 0.999997 | 6.33e-13 | 7.24 / 7.18 | 47.7 | 40.67 | 2.7e-4 |
| 1e-9 | 0.99972 | 6.33e-11 | 724 / 718 | 758 | 51.5 | 0.027 |
| 1e-8 | 0.9727 | 6.16e-9 | 69501 / 68902 | 68941 | 1109 | 2.59 |
| 1e-7 | 0.2627 | 1.663e-7 | 983804 / 968381 | 968387 | 28895 | 69.9 |
| 1e-6 | 0.00355 | 2.247e-7 | 176058 / 162892 | 162892 | 39045 | 94.4 |
| 1e-4 | 3.6e-7 | 2.255e-7 | 39325 / 39261 | 39261 | 39184 | 94.7 |

Small-ε scaling: R_thermal = 7.18e20 ε² (momentum), f_rest = 2.66e16 ε².

**The bound.** R_prim = R_0 B^{3/2} [D]: a cold band decays at ⟨a⟩ ∝ Γ_tot^{−1/2} with a fraction B into pairs. J ≥ 5Y_pair
[D]: aE = √(q² + a²m_A²) cannot decrease, so immediate decay minimises the daughters' conformal energy, and the prompt
integral ∫x⁴K₂ = 15π/2 is five times ∫x³K₁ = 3π/2. R_sat = 5C(3π/2)(Γ_X/H_m)F_R = 39183.974 (own). Since
R_sat > 1.5 R_0 (ratio 644), f(B) = R_0B^{3/2} + R_sat(1−B) is decreasing and its minimum over B ∈ [0, 1] is R_0 at B = 1
(grid check: 40.5626). The "visible outlet" cannot lower the raw hidden injection in this source model. CONFIRMED.

**Literature (credit, read levels).** Hall–Jedamzik–March-Russell–West, arXiv:0911.1120 (abstract read; the decay freeze-in
formula Y = 135 g Γ M_Pl/(8π³·1.66 g*^{3/2} m²) recalled [R]): Codex's C(3π/2)Γ/H_m is identical with g = 3, Γ = Γ_SM.
Chu–Hambye–Tytgat, arXiv:1112.0493 (abstract read: kinetic-mixing freeze-in regime). Redondo–Postma arXiv:0811.0326,
An–Pospelov–Pradler arXiv:1302.3884, Dvorkin–Lin–Schutz arXiv:1902.08623 (abstracts read): in-medium effects, including the
resonance where the thermal (Debye) mass of the mixed gauge boson equals m_A (here at T ≈ 2 m_A for hypercharge, m_D² ≈
(11/6)g'²T²) and plasmon decay, add to the vacuum inverse-decay rate. Codex's source omits these and the O(α) 2→2
channels; they increase the source, so the conclusions above only strengthen [A]. Fabbrichesi–Gabrielli–Lanfranchi is
context [Q Codex].

## 3. Reheating-temperature and vector-mass dependence (`code/reheating.py`, `freezein.py` e-rows)

Decaying inflaton (Γ_I) → radiation, H₀ = 10⁴Γ_I (T_max ≈ 6 T_R), constant g*, inverse decays on the actual T(a), no
back-reaction. Y_pair and J relative to hot RD: 1.000 (T_R = 1e9–1e7 GeV), 0.9997 (3e7), 0.994 (1e7), 0.928 (4e6), 0.687 (2e6),
0.259 (1e6), 0.099 (7e5), 0.0100 (4e5), 1.65e-4 (2e5), 1.4e-6 (1e5), 9e-12 (3e4). Below T_max ≈ m_A the suppression is
exponential (fitted d ln Y/d ln T_R ≈ 14 between 1e4 and 4e5 GeV for this T_max/T_R; for T_max ≫ m_A ≫ T_R it would be a
power law, Giudice–Kolb–Riotto 2001 [R]). Instant reheating at T_R: fraction ∫_{x_R}^∞x³K₁/(3π/2) = 0.95, 0.52, 0.18,
5.1e-3 at x_R = 1, 3, 5, 10. The outlet can lower the raw injection below R_0 only when R_sat(T_R) < 1.5R_0, i.e. between
T_R = 2e5 (yes, R_sat = 6.4) and 4e5 GeV (no, 390). The heavy window's radiation history has T_R = 3.08e14 GeV (prompt),
1e14–1e15 (finite release Γ_φ/H_i = 0.1–10) and Codex's six runs 5.3e13–1.7e14 GeV: x_R ≈ 1e-8, ratio = 1 to < 1e-6. A
T_R ≲ 3e5 GeV would be Γ_I ≲ 1e-18 H_i, a different cosmology whose entropy and Y_Q normalisation would have to be redone.

Vector mass: m_A = 6eh ∝ e, Γ_X/H_m ∝ e/(6h), B fixed by ε/e. At fixed ε_BX = 1e-7 the visible branch closes as e rises
(B = 0.76, 0.97, 0.997, 0.9997 at e = 3e-7, 1e-6, 3e-6, 1e-5) and R_thermal falls 9.8e5 → 9.4e5 → 2.2e5 → 4.6e4 → 7.8e3,
while the no-portal band R_0 ∝ 1/e falls 40.6 → 0.41. At fixed ε_BX/e = 1, R_thermal ∝ e: 1.7e6, 3.2e6, 5.8e6, 1.1e7
(f_rest 210–7000). The e ≳ 3e-6 remedy for the inflationary-vector gate therefore does not relax the ε_BX obligation.

## 4. Higgs-portal source (`code/higgs_source.py`)

L ⊃ −b v_H h s² from bφ²H†H ⇒ Γ(h → ss) = b²v_H²β/(8πm_h) = 4.81568e-20 GeV at b = 5e-11 (β = 0.99999). Thermal h, vacuum
masses, MB: Y_hss = 135Γ/(8π³g*s H(m_h)) [D, = HJMRW], two pairs per h (m_s = 0.2264 GeV > 2m_X = 0.1848, s → XX̄ open,
Γ = 1.7e-5 GeV, prompt): Y_pair = 2.23656e-8 (closed form 2.2365580e-8; own quadrature in T, 2.2365567e-8). The integrand
peaks at T ≈ m_h/2.5 (broken phase); replacing g* = 106.75 by g*(T) with MB-suppressed W, Z, h, t (g* = 95.5 at m_h/3):
Y_pair = 2.648e-8, f_rest = 11.1 (+18 %). f_rest = 2m_X Y_pair/4.4e-10 = 9.394 (Codex 9.394). Γ_h/H(m_h) = 1.9e11, so the
decay is prompt and the conformal daughter energy is J_h = 5Y_hss = 5.5914e-8 (quadrature/5Y = 1.0000003); in the spot
normalisation R_h = F_R J_h = 1943 (2450 with g*(T), including the (g*s,R/g*s)^{1/3} conversion). θ/b = 0.0290978 [D];
ATLAS B_inv < 0.107 one-channel: b_max = 5.0501e-3, λ_HX,max = 7.1419e-3 [C]. The vertex b v_H h φ² exists in either φ
phase, so the h-decay source is phase-independent; only the daughters' identity (s → XX̄ vs φ quanta) is phase-dependent
[A] (Codex's "thermal phases not certified" caveat applies to the daughters, not to the rate).

## 5. Cold leakage (`code/leakage.py`)

α₆ = κ²μh/(√2(M_Y⁴ − 2μ²h²)) = 2.10886e-31 GeV⁻² (leading 1.68182e-31), r = 0.45, exact/leading = 1/(1−r²) = 1.253918,
c6 = α₆M_r² = 1.250390e6 (hw_gates 1.250374e6, identical convention α₆ = μκ²⟨Σ⟩/(M_Y⁴ − 4μ²⟨Σ⟩²)). Our xcharge
Q = 200000 unit rates (c6 = 1 at M_r): coherent 3.91960e-50 /s, two-quantum anti-X pair 1.03997e-46 (free Born),
1.04901e-46 / 1.01252e-46 (LDA variants); Codex's "included" 1.040362e-46 = coherent + free-Born pair. Γt₀ (t₀ = 4.3550e17 s)
= 7.0837e-17 (free Born), 7.145e-17 / 6.897e-17 (LDA), 4.505e-17 with the leading denominator; pair/coherent 2653.27;
1 % ceiling c6 = 1.4856e13 (margin 1.19e7). As xcharge stated, the pair channel is a free-Born/LDA estimate with O(1)
medium corrections uncomputed; no in-medium lifetime is inferred.

## 6. Item 5: the window between detectability and injection (`code/headroom.py`, `code/headroom_e.py`)

Raw-source statement: the portal injection must fit under headroom = R_L − R_med − R_A, with R_L the case-L threshold
variants on record [Q own hw_critic thresholds.json: isentropic proxy 87.85 (table) / 81.5–82.0 (direct) / 80 (Codex
calibration); energy-conserving ke 8.26 (μ_c = 0.861 m) / 6.97 (0.923 m); ke + s route 19.6], R_med = 1.1437 [Q own] and
R_A the inflationary vector at e = 1e-7 (Codex band 40.56 [Q], own modes 48.47, GMR formula 64.0 [Q own]).

| e | threshold | headroom | ε_BX,max (R_thermal ≤ headroom) | b_Hφ,max (R_h ≤ headroom) |
|---|---|---|---|---|
| 1e-7 | isentropic proxy (80–87.85) | 14.9–46.1 | 1.4e-10 – 2.5e-10 | 4.4e-12 – 7.7e-12 |
| 1e-7 | ke 8.26 / 6.97 | −33 to −58 | none (the point already fails) | none |
| 1e-7 | ke + s 19.6 | −22 to −46 | none | none |
| 3e-6 | isentropic | 78.5–79.0 | 4.2e-9 | 1.0e-11 |
| 3e-6 | ke 8.26 | 5.3–5.8 | 1.1e-9 (ε/e ≤ 3.7e-4) | 2.6e-12 – 2.7e-12 |
| 3e-6 | ke + s 19.6 | 16.6–17.1 | 1.9e-9 – 2.0e-9 | 4.6e-12 – 4.7e-12 |
| 1e-5 | isentropic / ke / ke + s | 80 / 6.8–6.9 / 18.1–18.3 | 1.0e-8 / 3.0e-9 / 5.0e-9 | 1.0e-11 / 3.0e-12 / 4.8e-12 |

Rest-mass-only statements (f_rest ≤ 0.01 allowance, no energy processing at all): ε_BX ≤ 6.1e-10 (e = 1e-7), b ≤ 1.63e-12
(any e). For the Higgs source the rest-mass criterion is the tighter one because its pairs are soft (R_h/f_rest = 207
against 1.4e4 for the vector source).

**Consequence for the empirical branch.** The fa3bc97/8b0d96b proposal needed κ = 2b at the detector scale: κ_90 =
5.0e-10 (XENONnT 20 t·y projection median), 6.0–7.7e-10 (two LZ estimates), 2.5e-9 (XENON1T) for the Q = 39944, v = 1 GeV
object [Q 8b0d96b], i.e. b = 2.5e-10 – 1.25e-9; rescaled to the heavy-point object (E_Q = 17064.5 GeV, v = 0.924; rate per
κ² ∝ E_Q/(λv²)², ×6.4; form factor not redone) b ≈ 1.0–1.5e-10 (XENONnT/LZ level) to 4.9e-10 [E]. The 8b0d96b history was
cold by construction (T_max = 20 MeV), so this source was absent there; the heavy window's bath (T_R ≈ 3e14 GeV) passes
through T ≈ m_h, where it is unavoidable. At b = 1.0–1.5e-10 the raw source is Y_pair = 0.9–2.0e-7, f_rest = 37–87 and
R_h = 7.6e3–1.8e4, against the best raw-source allowance b ≤ 7.7e-12 (energy headroom, isentropic proxy, Codex band) or
b ≤ 1.6e-12 (rest mass alone). **The detectable Higgs portal overshoots the heavy-point raw-source window by a factor
13–60 in b (170–4000 in injection; ≥ 2000 in unprocessed rest mass versus the 1 % allowance).** Only a kinetic erasure of
≳ 99.9 % of the injected pairs (kin_* territory, not assumed here) could reopen it; under the energy-conserving thresholds
at e = 1e-7 there is no window regardless of the portal. The kinetic-mixing portal has no detector relevance at this
m_A (5e-13 events/t·yr even at ε_γ = 1) and its "visible outlet" cannot reduce the hidden injection (item 4); it is
constrained to ε_BX ≲ 2e-10 (isentropic) or ≲ 1–5e-9 at e = 3e-6–1e-5 (ke variants).

## 7. Corrections and notes for Codex

1. R_thermal(1e-7) = 9.68e5, R_total = 9.68e5 (momentum-resolved), not 9.84e5: the MB time-profile closure is 1–2 % high
   at ε ≤ 1e-7 and up to 7.5 % high at ε = 1e-6; the sign is as expected (free-streaming A are colder than thermal). No
   conclusion changes.
2. "Thermal phases are not certified" for the Higgs source should be narrowed: the h → φφ rate b²v_H²/(8πm_h) holds in
   either φ phase; only the daughters' subsequent identity is phase-dependent.
3. The ε_BX statement can be sharpened from "not forced at one loop" to "forbidden at all orders by hidden C of the stated
   tree Lagrangian (phases of κ, μ removable)"; the loophole is a CP-violating holder with a hidden–SM bridge.
4. The hot-RD caveat is inactive for every history in the packet and in our heavy-window tracks (T_R ≥ 5e13 GeV); the
   outlet would need T_R ≲ 3e5 GeV.
5. The g* = 106.75 constant in the Higgs source is 18 % low in Y_pair relative to g*(T ≈ m_h/3).
6. g_Y = 0.37 vs the one-loop-run 0.38 at m_A changes B from 0.263 to 0.252 (immaterial).

## 8. Limits of this track

Same-session AI work. Fixed vector mass, constant g* in the vector source, MB statistics, vacuum widths, no 2→2 or
in-medium production (both would raise the source), no hidden inverse reactions, no plasma kinetics, no f_sym, no
formation. The headroom numbers inherit the equilibrium-proxy thresholds and the e = 1e-7 vector estimates of our earlier
tracks, and E_X, R_0 are Codex's normalisations. The detector-level b for the heavy object is a rescaling [E], not a recast.
Nothing here is an observational statement.

## Reproduction (~25 s, 4 cores)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/rp_portals/code
cd $T
python3 -B widths.py ../results/widths.json
python3 -B freezein.py ../results/freezein.json
python3 -B reheating.py ../results/reheating.json
python3 -B higgs_source.py ../results/higgs_source.json
python3 -B leakage.py ../results/leakage.json
python3 -B headroom.py ../results/freezein.json ../results/higgs_source.json ../results/headroom.json
python3 -B headroom_e.py ../results/higgs_source.json ../results/headroom_e.json
```
Inputs read as data: `tracks/xcharge/results/{forced,pairs}_Q200000.json`, `tracks/hw_gates/results/gate5_light.json`,
`tracks/hw_critic/results/thresholds.json`. Python 3, numpy 2.5.3, scipy 1.18.1, sympy 1.14.
