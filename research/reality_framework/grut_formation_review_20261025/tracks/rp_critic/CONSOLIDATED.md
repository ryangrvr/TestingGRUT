# Consolidated review of Codex e160356 (physical reheating, constrained preparation, portal sources) — track rp_critic

26 October 2026 · Claude (issue #13) · consolidates tracks `rp_reheating`, `rp_preparation` and `rp_portals` (their NOTES.md read
in full), cross-checks them against each other and against the packet, and recomputes the load-bearing numbers with own code
(`code/`, outputs in `results/`). Packet: Codex e160356 (`origin/codex/grut-curvature-boundary-20261026`,
`research/reality_framework/grut_reheating_preparation_20261026/`), read only through `git archive e160356` into the scratchpad;
no Codex code was run by this track; nothing was committed; the kinetic tracks (kin_*) were not touched.

**Scope.** One supplied tree-level working model at the heavy-mediator spot point (H_i = 1.3334962e11 GeV, M_Y = 30 H_i, h = 40 H_i,
μ = 7.15945616 H_i, κ = 3.458684721e-3, r = 0.45, e = 1e-7, c_H = 10 ↔ ξ_X = −5/6; light sector g_X = 0.1, λ_X = 1e-3, λ_φ = 0.03,
v = 0.924089637 GeV, Q = 200000; Y_Q = 4.8294e-9), a quadratic-inflaton surrogate for reheating, a three-mode toy for the
preparation question and hot-RD source ledgers for the portals. Nothing here is evidence that X, Y, Σ, φ, the vector or the
inflaton exist. No observational statement, no validation, no exclusion, no GRUT prediction; Candidate 1 is unchanged and
paused; the vacuum density remains underived. All four agents (Codex and the three tracks) and this one are same-session AI
work: this is not outside replication. Labels: [D] derived, [C] computed (own code), [Q] quoted, [E] estimate, [A] argued.

## 0. Bottom line

Every number in e160356 that was checked reproduces, by three independent engines for the reheating runs (Codex's C DOP853
with H evolved; rp_reheating's C DP5(4) with H reconstructed from the constraint, plus a SciPy formulation; this track's
own SciPy DOP853 with H evolved), by three independent builds of the boundary model, and by two independent momentum-resolved
source solvers. **No numerical error was found in the packet.** Four things are corrected in substance or scope: (i) the
yield deficit is a property of Codex's slow-reheating grid (Γ_I/m_I ≤ 0.15), not of the model — a matter era of about one
Hubble time (Γ_I ≈ 0.2–0.3 m_I) gives 1.01–1.05 Y_Q with no hidden coupling changed; (ii) the thermal vector source is
R_total(ε_BX = 1e-7) = 9.68e5, not 9.84e5 (the MB closure is 1.6% high); (iii) the "phase not selected" result is a theorem
for conjugation-even conditions, but the symmetry is antiunitary conjugation in the number basis (time-reversal type), not
C/CP, and the mutual-information row is an assumption, not a test; (iv) the matter-era curvature mass is −2.5H² on average
(oscillating −10 … +5 H²), which corrects the physical content of our own hw_verify finite-release row. The window's status:
it survives a physical reheating on a line in the surrogate's two parameters, keeps its supplied phase and holder, and
acquires a new obligation for any portal that would make it detectable (§5).

## 1. Verdicts over Codex's four handoff items

| # | Codex claim (CLAUDE_HANDOFF.md) | Tracks | This track (own code) | Verdict |
|---|---|---|---|---|
| 1a | Jordan-frame equations (Friedmann with 3HḞ, Raychaudhuri, trace-eliminated R, transfer Γİ², work identities, Bianchi) | rp_reheating re-derived line by line; no missing term, no sign error [D] | own formulation from MATHEMATICS.md; own preparation solve lands on Codex's initial state to 2e-11 (x = 76.39978266, y = −2.11150985 − 0.93932904 i, σ = 28.29314213, R(0) = 12.000) [C] | CONFIRMED |
| 1b | Six runs: ε −0.191 … −0.219 / −0.180 … −0.181; N_YL 0.00352 / 1.61–1.74; yield 0.333–0.784 / 0.185–0.518 Y_Q | identical to printed digits (worst 1.4e-4 on N_YH) | b2_g0.1: ε −0.198769, yield 0.55956, N_YL 0.0035065; b2_g0.3: −0.218682, 0.78427; b10_g0.1: −0.181528, 0.32505, N_YL 1.70713 — identical to Codex's digits [C] | CONFIRMED |
| 1c | Radau/DOP853 3.2e-12; Friedmann ≤ 1.27e-9; long run 0.014%/0.22%; phase ±1e-8 → 4.4e-10; ξ = 0 and RD controls | all reproduced (rp_reheating) | not repeated (controls, not load-bearing) | CONFIRMED [Q rp_reheating] |
| 1d | "None of these six unchanged-scale runs produces the required yield; no retuning performed" | true for the six; same model gives 0.94–1.15 Y_Q along Γ_I ≈ 0.2–0.3 m_I, m_I = 3–50 H_i | (m_I, Γ_I)/H_i = (3, 0.5): 1.041; (5, 1): 1.041; (10, 3): 1.014 Y_Q at τ = 320–400, window spread ±0.5% for τ ≥ 40 (rp_reheating's 1.045/1.047/1.019 are the τ = 64–80 window) [C] | CONFIRMED as stated; **CORRECTED in scope** (grid, not model) |
| 1e | m_I/H_i = 10 "produces more Y excitation" | 2m_I ≈ ω_YL = 22.25 H_i resonance, peak N_YL 2.22 at m_I = 12 | N_YL(m_I = 12, Γ = 0.1) = 2.2233 [C] | CONFIRMED and explained |
| 1f | Preparation: conditional minimum at the actual R(0), phases supplied, no holder | identical state; holder off implicitly at τ = 0; 0.048 light Y quanta | identical state (above) | CONFIRMED as conditional |
| 1g | Finite-k, decays, fragmentation, thermal bath not certified | not computed | not computed | NOT CHECKED (open, §6) |
| 2a | Neutral projected coherent state: dims 31/58; equal E, occupations, mode entropies at δ = ±π/2; torques ∓1.3563/∓1.3576; doubled ∓2.7153 (cutoff change 0.102%) | rp_preparation: all digits; converged torque 1.35764502; Codex's "refined 2.71528939" is the (24,12,9) value, converged 2.7152900 | own third build: dims 31/58/88; E(0, ±π/2, π) = 14.311195, 14.084921, 13.858647 (cutoff 2); T(+π/2) = −1.35626541 / −1.35764469 / −1.35764502; eigenstate torques 0 [C] | CONFIRMED |
| 2b | Clock POVM: CV = 0, effects positive, commute with C, Born match 0.7184909 | reproduced; spectrum 1 − cos(πk/5), gap 0.190983 | not repeated | CONFIRMED [Q rp_preparation] |
| 2c | Hidden/SM mutual information 0 "also fails to choose a phase" | imposed by the product form ρ_h ⊗ diag(.2,.8); correlated SM bit gives T_net = (2p − 1)T₊ | agreed [A] | CORRECTED (assumption, not test) |
| 2d | "Neutrality, these entropy conditions and an uncorrelated SM factor do not select the charge-producing phase" | theorem: K = complex conjugation in the number basis, KHK = H, Kψ_δ = ψ_{−δ}; all K-even conditions are blind to the sign; T = 6 dE/dδ so energy extremisation selects δ = π with T = 0; stationary rules give T = 0 | T = 6 dE/dδ checked to 2e-9; T(δ) = −12·0.04·2√2 sin δ = −1.357645 sin δ and E(δ) = E₀ + 0.226274 cos δ exact (E(0) − E(π))/2 = 0.2262742 at cutoff 2 [C, D] | CONFIRMED as theorem; label **CORRECTED** (conjugation of time-reversal type, not C/CP: no antiparticle modes) |
| 2e | Weyl ≡ 0 in FLRW; real boundary phase keeps probabilities; scope sentences | trivial; agreed | — | CONFIRMED (trivial); boundary term otherwise NOT CHECKED |
| 3a | Γ_SM = 41ε²g_Y²m_A/(96π), Γ_X, B = 0.262712 at ε_BX = 1e-7 | field-by-field ΣY² = 10 + ½ | B = 0.26271198 [C] | CONFIRMED |
| 3b | Y_pair = 1.663e-7 = C(3π/2)(Γ_X/H_m)(1 − B) | both solvers | 1.662942e-7 analytic, momentum-resolved, closure [C] | CONFIRMED |
| 3c | R_thermal ≈ 9.84e5, R_total = 983809.9 (MB time-profile closure) | momentum-resolved 968381 (−1.6%) | own momentum solver 968381.2; own closure 983804.4; ratio 0.9843; delay J/5Y_pair 33.52 vs 34.05 [C] | **CORRECTED (small): R_total = 9.68e5** |
| 3d | R_total ≥ 40.5626 B^{3/2} + 39184(1 − B) ≥ 40.5626 | re-derived; R_sat/1.5R_0 = 644 | J ≥ 5Y_pair holds in every own run (J/5Y = 1.13 at ε = 1e-5) [C] | CONFIRMED |
| 3e | Hot-RD, fixed-mass source; low T_R or other m_A needs separate treatment | T_R scan: ratio 1.00 for T_R ≥ 1e7 GeV; every history on the table has T_R ≥ 5e13 GeV | T_RH prefactor 0.5406; six runs 5.3e13–1.7e14 GeV; the line 2.2e14–5.3e14 GeV [C] | CONFIRMED; caveat inactive |
| 4a | Higgs source b = 5e-11: Γ_hss 4.81568e-20 GeV, Y_pair 2.237e-8, f_rest 9.39 | same; g*(T) raises Y_pair 18%; R_h = 1943 | Γ_hss 4.81568e-20, Y_pair 2.23656e-8, f_rest 9.394, R_h 1942.8 (own arithmetic) [C] | CONFIRMED |
| 4b | θ/b, ATLAS b_max 5.05e-3, λ_max 7.14e-3; mixed trace 0; β_b = g_X²λ_HX/(8π²) | identical; sympy CW derivation | not repeated | CONFIRMED [Q rp_portals] |
| 4c | Cold leakage α₆ 2.10886e-31 GeV⁻², 1/(1 − r²) = 1.25392, Γt₀ = 7.08368e-17 | 7.0837e-17 from our xcharge rates; LDA 6.90–7.15e-17 | not repeated | CONFIRMED [Q rp_portals, xcharge] |

## 2. Cross-checks between the tracks and against e160356 (own evaluation)

| Quantity | Codex [Q] | rp_* track [Q] | This track [C] | Status |
|---|---|---|---|---|
| Initial state x, Re y, Im y, Im s | 76.39978266076902, −2.111509854283515, −0.9393290404896636, −28.293142125868524 | identical (own Newton solve) | 76.39978265923669, −2.1115098541668, −0.9393290405025, −28.29314212586876 (residual 4e-12) | agree ≤ 2e-11 |
| Yield bookkeeping yield/Y_Q = (\|N_X\|/2788.18)/D, D = (ρ_r a⁴/3)^{3/4} | yield from s_vis(ρ_r,final) | decomposition | D from Codex's recorded final states reproduces every Codex yield to 5 digits; direct Y in GeV units identical; RD reference 0.99980 Y_Q, so Y_Q is the instant-release yield | agree |
| Dilution of the six runs | — | 13.4 / 7.86 / 5.51; 6.43 / 3.66 / 2.31 | 13.395 / 7.861 / 5.507; 6.425 / 3.657 / 2.311 | agree |
| Charge timing | "charge fixed early" (implicit) | 90% by τ = 2.6 (m_I = 2) / 0.5 (10) | 90% at τ = 2.51 (radiation 13% of ρ) / 0.40 (3%); 99% at τ = 3.1 / 2.4 | agree |
| ξR/H² range realised | qualitative | −10 … +4.3 (m_I = 2), +4.7 (10) | −10 … +0.65 (m_I = 2, Γ = 0.1), +3.16 (10, 0.1), +3.40 (12), +1.85 (10, 3) over the sampled grid (20001 points; peaks between samples not resolved) | consistent (range [−10, +5] [D]) |
| Yield line (m_I, Γ_I)/H_i = (3, 0.5) / (5, 1) / (10, 3) | not run | 1.045 / 1.047 / 1.019 (τ = 64–80) | τ = 40–50: 1.037 / 1.036 / 1.005; 64–80: 1.045 / 1.047 / 1.019; 160–200: 1.042 / 1.042 / 1.013; 320–400: 1.041 / 1.041 / 1.014. With Codex's Γt = 25 convention (τ_end = 50 / 25 / 8.3): 1.037 / 1.034 / **1.070** | agree; late means are window-dependent at ±0.5% (τ ≥ 40) and by +5% at Γt = 25 for (10, 3) |
| Y resonance peak | "much larger" at m_I = 10 | N_YL = 2.22 at m_I = 12, Γ = 0.1 | 2.2233 | agree |
| T_RH on the line | T_vis,final quoted (end-of-run temperature) | 2.2e14 (Γ = 0.5 H_i) – 5.3e14 GeV (3 H_i) | 0.5406 √(Γ_I M_r): 2.18e14, 3.08e14, 5.34e14 GeV at Γ_I = 0.5, 1, 3 H_i | agree |
| Vector source at ε_BX = 1e-9 / 1e-8 / 1e-7 / 3e-7 / 1e-6 / 1e-5 | closure 724.1 / 69501 / 983804 / 502714 / 176058 / 46597 | momentum 718 / 68902 / 968381 / — / 162892 / — | momentum 717.9 / 68902 / 968381 / 485009 / 162892 / 44326; closure 724.1 / 69501 / 983804 / 502714 / 176058 / 46597 | agree; closure 0.9–7.5% high |
| Boundary model: projection probability, occupations | 0.0410786 / 0.0408848; 5.8335, 1.7850, 1.8647 | same; untruncated 0.0408829 | 0.04107863 / 0.04088476 / 0.04088286 (30,15,12); occupations identical | agree |
| ⟨cb²⟩_proj = α_cα_b² = 2√2 e^{iδ} | — | 7e-14 at (36,18,15) | 2.8255529 i (18,8,6), 2.8284264 i (24,12,9), 2.8284271 i (30,15,12) vs 2√2 = 2.8284271 | agree (truncation at cutoff 1) |
| Higgs source | 4.81568e-20; 2.23656e-8; 9.394 | same; R_h 1943 | same to all digits; R_h 1942.8 | agree |
| Debris R_central of the m_I = 10 runs and on the line (5.2–8.0; 0.008–2.2) | not in packet | ledger reuse [E] | not recomputed | NOT CHECKED (ledger proxy) |

Own integration checks (`results/rc_*.json`): Gauss ≤ 7e-8 absolute on N ~ 1e4; relative Friedmann residual ≤ 1e-8 for
τ_end ≤ 83 and 7e-6 (τ = 250) to 2e-5 (τ = 400) from the late-time drift of the independently evolved H at atol 1e-13
(ρ has fallen by 1e-6 to 1e-8 by then; a_f, ε, N_X and the yield agree with the constraint-solved rp_reheating engine and
with Codex to ≤ 1e-9, so the drift is a diagnostic of my code, not of the physics); preparation residual 4e-12.

## 3. Corrections for Codex (consolidated)

1. **R_total(ε_BX = 1e-7) = 9.68e5, not 9.84e5.** The MB time-profile closure assigns the vector population the thermal
   ⟨m/E⟩ at the current x; the free-streaming A keep p/T ≈ 3, decay slightly earlier and deposit 1.6% less conformal energy
   (7.5% at ε = 1e-6). Two independent momentum-resolved solvers (rp_portals, this track) agree to 1e-6. Nothing downstream changes.
2. **The yield deficit is the grid's, not the model's.** Charge is fixed by τ ≲ 3 (90% while radiation is ≤ 13% of ρ);
   the entropy produced afterwards dilutes it by D = 2.3–13.4. With Γ_I ≈ 0.2–0.3 m_I (matter era ≲ 1 Hubble time) the same
   model with the same hidden couplings, scale and phase gives 1.01–1.05 Y_Q for m_I = 3–50 H_i, T_RH ≈ 2–5e14 GeV. "No
   retuning performed" is accurate; "none of the six reaches Y_Q" should carry the cause and the grid label.
3. **The m_I = 10 Y excitation is a resonance,** 2m_I ≈ ω_YL = (M_Y² − 2μ|Σ|)^{1/2} = 22.25 H_i (heavy mode 36.1 H_i at
   m_I ≈ 19), peaking at N_YL = 2.22 for m_I = 12; it is set by τ ≲ 2 and then conserved (no Y decay in the equations).
   Under our ledger it makes the m_I = 6–16 runs marginal debris cases (R ≈ 2–8 vs 1.14 instant [E]); m_I ≲ 5 H_i is quiet.
4. **Late-mean convention.** Γt = 25 is adequate for the six runs (τ_end ≥ 83) but too short for Γ_I ≳ H_i: at (10, 3) the
   τ = 6.7–8.3 mean is 5% above the τ ≥ 40 value because N_X still relaxes after reheating (through the excited Y). The
   τ-weighted last-20% mean oscillates by ±0.5% between windows for τ ≥ 40; quote line yields to two digits.
5. **Preparation labels.** The symmetry used is antiunitary conjugation in the occupation basis (time-reversal type; no mode
   has an antiparticle partner, so C is not an operator on the space), under which E, occupations, entropies, MI and Q are
   even and the torque is odd; the statement is a theorem for any such condition set under a real H, not an example. The
   MI = 0 row is imposed by the product form, not tested. Energy minimisation does select a phase — the torque-free δ = π
   (T = 6 dE/dδ) — so "energy does not select" should read "energy selects a torque-free phase". Stationary rules (eigenstates,
   Gibbs, equivariant F(H)) give T = 0 unconditionally. The "refined torque 2.71528939" is the (24,12,9) value; converged
   2.7152900. SOURCES.json pins none of Page–Wootters, conditional-probability clocks or projected-coherent-state literature.
6. **Portal notes.** The h → φφ rate b²v_H²β/(8πm_h) holds in either φ phase (the daughters' identity is phase-dependent);
   ε_BX is forbidden at all orders by the hidden charge conjugation of the stated tree Lagrangian (κ, μ phases removable),
   so "not forced at one loop" understates it — the loophole is a CP-violating holder plus a hidden–SM bridge [A]; g*(T) at
   T ≈ m_h/3 raises the Higgs source 18%; the hot-RD caveat would need T_R ≲ 3e5 GeV; g_Y run to m_A gives B = 0.252 (immaterial).

## 4. Overreach flags on our own tracks

- **rp_reheating, item 10.** The three-digit line yields are window-dependent (±0.5%, §2); quote 1.01–1.05. "Without
  retuning" is true for every hidden coupling, the scale and the phase, but the line is the locus where the yield crosses
  Y_Q in the surrogate's own two parameters; since yield(Γ_I = 0.3 m_I) ≥ 1 for m_I ≥ 3 and the yield is monotonic in Γ_I,
  a crossing exists for every m_I ≥ 3 — the content is "the model's perturbative range spans 0.18–1.15 Y_Q", not a prediction
  of Y_Q. Γ_I/m_I = 0.2–0.3 is the edge of perturbative decay (no preheating modelled); T_RH = 2–5e14 GeV coincides with the
  heavy window's prompt-release value, which is consistent, not independent. Debris on the line is ledger reuse [E].
- **rp_reheating, item 12 (our hw_verify row).** Accepted: the finite-release row "Γ_φ/H_i = 1 → ε −0.279, 1.16 Y_Q, R 0.03"
  stands in its model but not in the covariant one (−2.5H² average, not −10H²); hw_critic's "a finite release is favourable
  (R falls to 0.03 and |ε| rises to 0.38)" should be withdrawn in favour of "|ε| stays 0.15–0.22; R stays small only off the
  Y resonance".
- **rp_portals, §6 (headroom and detectability).** The headroom R_L − R_med − R_A is a debris proxy built for hard pairs
  injected *before onset* (mediator and vector decays complete 220× before the case-L coordinate a_ref, T_vis ≈ 46 GeV, where
  the φ transition and formation begin [Q kin_source]). The thermal vector source decays on the same clock as the primordial
  vector and is inside that domain. The Higgs source is not: h → ss peaks at T ≈ m_h/2.5 ≈ 50 GeV, i.e. *at onset*, during
  the φ transition, so the b_Hφ headroom (≤ 4–8e-12) and "overshoots by 13–60 in b" compare a source to a threshold outside
  its domain, and "only ≳ 99.9% kinetic erasure (kin_*) could reopen it" misattributes: kin_* as specified has all sources
  complete before onset. The robust part is the unprocessed rest mass: f_rest = 9.4 at b = 5e-11 and 37–87 at the
  detector-level b ≈ 1.0–1.5e-10 [E rescaling], so a detectable Higgs portal requires ≥ 97–99% erasure of pairs injected at
  onset — an uncomputed obligation on the empirical branch, not a closure. The vector-source statement ε_BX ≤ 1.4–2.5e-10
  (isentropic proxy) is in domain but inherits the undecided threshold (no window at e = 1e-7 under the energy-conserving proxy).
- **rp_preparation.** No overreach found; "T = 6 dE/dδ ⇒ torque-free extremum" holds along the L-rotation direction, as stated.
- **All three tracks and this one:** same-session AI work with shared conventions (yield, ledger, thresholds); agreement to
  1e-9 between engines written in the same session is a check on implementation, not on the model.

## 5. What e160356 changes about the heavy-mediator window

| Aspect | Before (hw_critic/CONSOLIDATED.md) | After this packet and review |
|---|---|---|
| Release history | prompt radiation (RD) or a fluid inflaton with −c_H H² through the matter era | a physical oscillating inflaton with the covariant ξ_X R mass (−10 … +5 H², matter average −2.5H²); the quiet transfer survives; the charge is fixed by τ ≲ 3 and then diluted |
| Yield | 1.00 Y_Q by construction at the prompt release; 1.16 at Γ_φ = H_i (now withdrawn) | 0.18–0.78 Y_Q on Codex's grid; 1.01–1.05 Y_Q along Γ_I ≈ 0.2–0.3 m_I. **New stated requirement:** reheating within about one Hubble time (T_RH ≈ 2–5e14 GeV) and m_I ≲ 5 H_i (or ≳ 30) to avoid the Y resonance; m_I = 6–16 H_i gives R ≈ 2–8 of mediator debris [E] |
| Preparation | phase supplied through a holder; held-phase attractor; gate 2: equilibrium gives P(π/2)/P(0) = e^{−43,277} | unchanged in substance, now with a theorem: no conjugation-even condition (energy, neutrality, entropies, MI, stationarity) can select the sign; energy extremisation picks the torque-free phase. The holder with a supplied CP phase remains the only route; its switch-off is still assumed to coincide with the end of de Sitter |
| Debris and threshold | mediator 1.14; vector 40.6–64 at e = 1e-7; decisive calculation kin_* | unchanged; plus the resonance debris if m_I is in 6–16 |
| Portals | b_Hφ = 5e-11 detector proposal (8b0d96b/fa3bc97); ε_BX unconstrained | **new obligation:** any ε_BX adds a thermal vector source R ≈ 7.2e20 ε² (momentum-resolved), so ε_BX ≤ 1.4–2.5e-10 under the isentropic proxy and none at e = 1e-7 under the energy-conserving one; the "outlet" cannot lower the raw band (bound). The Higgs portal at a detectable b injects f_rest = 37–87 at onset: unresolved. ε_BX = 0 and b = 0 are radiatively consistent choices for the no-portal window [A] |
| Status line | exists under homogeneous classical release; requires supplied phase, supplied holder, one undecided threshold | exists also under a physical reheating on the line above; requires in addition Γ_I ≳ 0.2 m_I with m_I ≲ 5 H_i; a detectable portal adds an uncomputed onset-epoch erasure requirement. Not a cosmological point; no observational evidence; Candidate 1 unchanged; vacuum density underived; closure supplied no new law (Codex's sentence is accurate) |

## 6. Open

- The phase value and the holder are supplied; the de Sitter stage is not evolved; the holder is assumed to switch off at τ = 0.
- Homogeneous, classical, tree level: no finite-k modes, particle production, mediator decays in the dynamics, backreaction,
  thermal corrections or fragmentation on the reheating backgrounds (Codex: `physical_modes_stability_certified: false`).
- Debris on the line and for m_I = 6–16 is ledger-and-proxy [E]; the threshold is kin_*'s (running, not touched).
- Portal sources: constant g*, MB statistics, vacuum widths, no 2→2 or in-medium production (both raise the source), no
  hidden inverse reactions; the φ phase at T ≈ m_h is not certified; the onset-epoch fate of the Higgs-source pairs is uncomputed.
- The preparation theorem is about the listed conditions in a toy with no gravity, gauge field or antiparticles; whether any
  cosmological preparation rule GRUT adopts is conjugation-even is a question about that rule.
- Everything is same-session AI work; nothing is observational.

## 7. The single most decisive next calculation (other than kin_*)

**Linear finite-k stability of X, Y, Σ and the longitudinal vector on the physical reheating backgrounds** — Floquet/mode
equations on the recorded homogeneous solutions in the quiet-and-sufficient corner (m_I = 3–5 H_i, Γ_I = 0.2–0.3 m_I) and on
Codex's m_I = 10 runs, through the first few inflaton oscillations where the X curvature mass swings between −10H² and +5H² at
frequency 2m_I with |X| ≈ 30–76 H_i, and through the 2m_I ≈ ω_YL resonance. It decides whether the homogeneous quiet transfer
(fixed by τ ≲ 3) and the 1.0 Y_Q yield survive inhomogeneous growth — the gate the packet itself leaves uncertified and the one
that could remove the window independently of the threshold; it is cheap (linear modes on stored backgrounds) and does not
overlap kin_*. Runner-up, for the empirical branch and as an extension of the kin_* solver rather than a duplicate: the
onset-epoch fate of the Higgs-source pairs (injected at T_vis ≈ 50 GeV during the φ transition), which decides whether any
detector-level b is compatible with the dark-matter budget.

## Reproduction (own code; ≈ 12 min on 4 cores)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/rp_critic/code
cd $T
python3 -B rc_reheat.py ../results/rc_b2_g0.1.json 2 0.1 250          # Codex b2_g0.1 (≈ 2 min); likewise 2 0.3 83.3333, 10 0.1 250, 12 0.1 250
python3 -B rc_reheat.py ../results/rc_b3_g0.5_t400.json 3 0.5 400     # yield line with window diagnostics; likewise 5 1 400, 10 3 400
python3 -B rc_reheat.py ../results/rc_b3_g0.5.json 3 0.5 50           # Gamma*t = 25 convention; likewise 5 1 25, 10 3 8.3333, 3 0.9 27.8
python3 -B rc_source.py ../results/rc_source.json                      # momentum-resolved vs closure vector source (30 s)
python3 -B rc_boundary.py ../results/rc_boundary.json                  # three-mode boundary model (20 s)
python3 -B rc_arith.py <scratch>/.../REHEATING_RESULTS.json ../results/rc_arith.json
python3 -B rc_summary.py                                               # -> results/summary.json
```
Codex's packet is taken from `git archive e160356 research/reality_framework/grut_reheating_preparation_20261026` into the
scratchpad. Reuse from our earlier tracks: none in code; the Y free-mode diagnostic convention (hw_verify/Codex) is
re-implemented from its formula; thresholds, R_med, R_A and the onset coordinate are quoted from hw_critic, hw_gates and
kin_source. Credits: Dine–Randall–Thomas and Affleck–Dine [R]; scalar–tensor FLRW equations (standard) [R]; Hall–Jedamzik–
March-Russell–West decay freeze-in [R]; Page–Wootters [R]; Sakharov [R]; Penrose [R]; Codex e160356 [Q].
