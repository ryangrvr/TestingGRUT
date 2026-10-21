# CRITIC: completeness and consistency of the issue #13 wall-history packet (verifier key v_critic)

21 October 2026, state of the packet at 04:20 UTC. I read BRIEF.md, all seven module NOTES.md, their results
headline JSON, the six finished verifier reports (v_thermal, v_charge, v_heating, v_inflation, v_network, v_cosmo),
and Codex's FEASIBILITY.md and README.md at 135d535 (`git show`, read-only). The module `charged_exit_sim` is still
running (2 workers since 03:07; 14 of 34 runs written; no NOTES.md yet). I treat its results as **preliminary** and
read only its `results/summary.json`, `validation.json` and the docstrings of `cesim.py`/`runs.py` to learn what it
simulates. `v_relic` was still running and has no report.

I did not recompute the modules. My own code (`critic_checks.py` -> `critic_results.json`, < 5 s, numpy/scipy) does
only the bookkeeping needed to reconcile modules: time-temperature conversions, the charged-exit temperature
using the charge_phase n_local(T_h) column (read as quoted data), a thin-wall thermal-nucleation estimate, the X
bound state on the wall, a spectator-condensate relic estimate, and the reference-history numbers. Labels: [D]
derived here, [C] computed here, [Q] quoted, [E] order-of-magnitude estimate, [P] outside-action preparation datum.

## 1. Cross-module numerical inconsistencies and the value to adopt

| # | Quantity | Values in the packet | Adopt | Why |
|---|---|---|---|---|
| I1 | Hidden T_c (zero charge) | thermal_potential 1.5571-1.5641 v (resummed AE/Parwani/DR, 4 scales; v_thermal CONFIRMED); charge_phase T_c0 = 1.5575 v (leading order with screened X loop); hidden_heating 1.5476 v "±30%"; network_sim 1.5476 v (illustrative); relic_yield both | **T_c = 1.557-1.564 v** (±0.3% soft-loop systematic), conditional on μ_N = 0 | Only resummed value. The ±30% band of hidden_heating is 100x wider than the actual uncertainty; its bounds move by only +0.35-0.6% (T_RH ∝ T_c^{4/7}). charge_phase's 1.5575 lies inside the range. The perturbative order is scheme dependent (v_thermal: gap scheme gives spurious φ_c/T_c = 0.026); the network-formation conclusion does not depend on the order (p = 1/2 either way). |
| I2 | Wall-domination epoch | network_sim t_w = 3.9e4 s; cosmo_observable T_dom = 6.367 keV, t_dom = 2.71e4 s; inflation_spectator T_dom = 8.0 keV | **T_dom = 6.37 keV (5.45-7.34), t_dom = 2.7e4 s** | Not a contradiction but three definitions. t_w = 3 M_P²/(4Aσ) = 3.88e4 s is the no-backreaction time and corresponds to 5.82 keV [C]; cosmo_observable includes backreaction (walls speed up H). 8.0 keV is σ/(M_P²H) = 1, a one-wall-per-Hubble-patch criterion [C: 7.97 keV], an O(1) convention appropriate only to isolated late-re-entering domains. The final report should quote 6.37 keV and label t_w as the no-backreaction estimate. |
| I3 | De Sitter thresholds | hidden_heating: equilibrium restoration H > 436 GeV (origin locally stable 379 GeV), finite duration 0.57 v (60 e-folds), "no macroscopic wrong-vacuum regions for H ≲ 184-275 GeV"; inflation_spectator: p_- < 1e-3 at 368/355 GeV (N*=50/60), N_reg < 1e-3 at 144/138 GeV, no re-entry after 1 MeV at 183-197 GeV, percolation (regime D) 1784/1607 GeV at end of inflation, 1.43/1.35 TeV after trapping | **inflation_spectator thresholds, with verifier corrections: H_perc 1792/1617 GeV (v_inflation); with X fluctuating from zero 1.70/1.47 TeV (v_heating); "no domain re-entering after BBN" ≈ 184 GeV; finite-duration X restoration 0.60 v at N = 60 (v_heating)** | The two modules answer different questions. 436 GeV needs X equilibrated (~200-300 e-folds) and is not a percolation threshold for an N_obs-scale window from a +v start. hidden_heating's 275 GeV end is a Δ = 5 e-fold count, not conservative (v_heating). |
| I4 | n_c usage (charge thresholds) | relic_yield uses the cold lower spinodal n_c = 0.0316 v³ as a "PLACEHOLDER" for both (a) the exit of a charge-held φ = 0 phase and (b) "no charge-driven restoration iff n_Q(T_Q) < n_c" for supply into a cold +v sector; charge_phase: n_c (φ = 0 locally stable), n_t = 0.0401 (homogeneous global switch), n_s1 = 0.0428 (φ ≠ 0 branch ceases), n_* = 0.178 (homogeneous φ = 0 is ground state), and n_c(T_h) rising to ≈ 0.042 v³ near T_h ≈ 1.5 v | (a) **exit at n = n_local(T_h) at the exit hidden temperature**: T_SM,exit = 203-241 TeV for r = 1 and 575-578 TeV for r = 0.05 (ξ = 5e-4 to 7.58e-4) [C], versus the cold formula 198 r^(-1/3) TeV (2-21% low). (b) **T_Q < 198 r^(-1/3) TeV is sufficient for no restoration; T_Q > 219 r^(-1/3) TeV (n > n_s1) forces it; the window in between is unresolved** (branch 1 is FLS-unstable there, so fragmentation is expected, but it was not simulated: charged_exit_sim runs κ = 0.9, 2, 4 only, while n_s1/n_c = 1.35). | n_t is a homogeneous bookkeeping marker with no dynamical role. Between n_* and n_c the homogeneous φ = 0 phase is metastable to phase separation, but nucleation is negligible: T = 0 S_4 ≥ 2.5e5 (charge_phase, v_charge) and my thin-wall thermal estimate S_3/T ≥ 2.6e3 for T_h ≤ 0.4 v [E] (critic_results.json). The spinodal exit at n_local(T_h) stands. |
| I5 | When does the charge, not the temperature, end the restored phase? | relic_yield: "for r ≳ 1e-3 the charge is Bose-condensed at T_c and keeps φ = 0 ... until T_SM = 198 r^(-1/3) TeV"; charge-holds criterion ξ < 7.8e-3 r^(1/3). charge_phase: on its ξ grid the charge route holds at ξ ≤ 0.005 and the exit is thermal at ξ = 0.0075 for r = 1 | **There is no sharp threshold: n_local(T_h)/T_h³ rises monotonically from 0 at T_c0 as T_h falls, so the exit happens where it meets the constant n_Q/T_h³ = 4.06e-9 r/ξ³.** At ξ = 7.58e-4: r = 1e-3 is not condensed (relic's own table: condensed fraction 0; r_BEC = 1.22e-3) and exits thermally; T_h,exit ≤ 1.0 v needs r ≥ 4.3e-3 and ≤ 0.75 v needs r ≥ 0.012 [C]. relic_yield's ξ < 7.8e-3 r^(1/3) uses n_c,cold/T_c³ = 0.0085 instead of n_local/T³ = 0.0125-0.019 near T_h = 1.3-1.5 v, so it overstates the charge-route ξ range by ~1.5x (consistent with charge_phase's 0.005-0.0075 bracket). | Both are far above the relic-compatible ξ ≤ 7.6e-4, so the practical consequence is small, but relic_yield's "r ≳ 1e-3" sentence contradicts its own capacity table. |
| I6 | ξ thresholds | kinetic: 6.6e-5 (thermal_potential, φ elastic), 4.3e-5 (hidden_heating), (3.6-5.4)e-5 (v_thermal exact kinematics); X kinetic 1.45e-4; D-φ exchange 0.0148; 2→4 chemical 0.068 (band 0.0084-0.27, [E]); broken-phase 3→2 ξ_32 = 3.48e-3 [C]; relic ξ_max = 7.58e-4 (8.36e-4 without D); ξ_BEC = 7.09e-3 r^(1/3) | **Kinetic window (4-7)e-5 (v/TeV)^(1/2) ≲ ξ ≤ 7.6e-4 (v/TeV)^(-1/3)** | No contradiction, but the final report must say that in the whole relic-compatible window the bath is never in number (chemical) equilibrium (both 2→4 and 3→2 thresholds are above ξ_max) and D is not in contact with φ (0.0148 > ξ_max), so μ_N and the D population are production inputs and T_c = 1.557-1.564 v is conditional on μ_N = 0. |
| I7 | Formation hierarchy | network_sim: "H^-1/δ = 2.1e13" at T = 1.5476 v and "formation hierarchy H^-1/ξ_f ~ 1e13" (C2) | **H^-1/δ ≈ 5.2e6 at the reference thermal exit (ξ = 5e-4, T_SM = 3.1e6 GeV) and 1.5e8 for the charged exit at r = 0.05 [C]; H^-1/ξ_f = N_H^(1/3) ≈ 10^4.6-10^5.2 (thermal, ξ = 5e-4, z = 1-2) and ≈ 2e6 (charged exit, charge_phase ξ_KZ ≈ 850 ℓ)** | 2.1e13 is the ξ = 1 (T_SM = T_h) value, which the relic bound excludes. network_sim C2 conflates the wall width δ with the domain size ξ_f. The extrapolation from the simulated H^-1/δ ≤ 72 (3D physical) is therefore ~5 decades for the reference, not 11. |
| I8 | Scaling amplitude A | network_sim 2D PRS 0.652 ± 0.002 (v_network: 0.69 ± 0.06 over 5 realizations); 3D PRS 0.777 ± 0.013; 3D physical 0.884 ± 0.062 (box-affected); cosmo_observable A = 0.8 ± 0.1 (HKS); preliminary charged 2D runs A_wall 1.37-2.35 vs neutral 0.68-1.06 | **A = 0.8 ± 0.1 for the neutral network** | Consistent. If the charged network kept A ≈ 1.4-2.4 (2D, not 3D, not scaling), T_dom ∝ A^(1/2) rises by ≤ 1.7x: the observable conclusion only strengthens. |
| I9 | Domain counts | thermal_potential N_H per H^-3; charge_phase per Hubble volume, τ_Q = 1/(3H), mean-field ν = 1/2; charged_exit_sim ξ_KZ 962-968 ℓ vs charge_phase 850 ℓ | state conventions | 0.62 dex volume convention and O(1) τ_Q choices; N ≫ 1 in every convention. |
| I10 | R_H range | cosmo_observable 1.14e9-2.06e9 (n = 1), 2.32e5-2.76e5 (n = 2), "at least 5.4 orders" | **1.04e9-2.30e9 and 2.17e5-3.44e5 (joint variation; v_cosmo); at least 10^5.3** | One-at-a-time ranges; conclusion unchanged. |
| I11 | Y_B, ω_b | Y_B = 8.674e-11 (charge_phase, relic) vs Codex 8.6758e-11; ω_b 0.02233 (BRIEF, BBN) vs 0.02237 (cosmo, Planck) | either | ≤ 0.7% convention spread (v_charge); irrelevant. |
| I12 | Symmetric-relic initial condition in the charged branch | relic_yield starts its Boltzmann system at T_h = 4 m_X = 0.4 v in the +v broken phase | flag | In the charged branch φ = 0 persists to T_h,exit ≈ 0.1-0.45 v [C], overlapping the start. Hidden number is frozen for ξ < ξ_32, so Ω_sym changes little, but the spectrum/conversions used are not those of that history. |
| I13 | Hidden energy budget at large r | ξ (thermal) treated as the whole hidden energy | flag | The branch-2 condensate energy at T_c is 5.9x (r = 1, ξ = 5e-4) or 1.1x (r = 1, ξ = 7.58e-4) the thermal bath energy, but only 0.02-0.11 for r = 0.05 [D]. For r ≲ 0.09 (retrospective DM: ω_asym ≥ 1.338 r) the thermal ξ is adequate. |

## 2. Codex's requirements: are they met?

**R1 Explicit history with all free inputs** (thermalisation mechanism, T_h/T_vis, hidden charge and chemical potential,
chronology). **PARTIALLY MET.** Every module names its inputs and labels outside-action data, and hidden_heating shows
the action itself couples the sectors only through gravity, which cannot thermalise the hidden sector
(Γ₂₂/H ~ 1e-26). Missing: (a) no single assembled history with numerical values exists in the packet (I supply one
in §5); (b) the producer of a μ_N ≈ 0 hidden bath at ξ ≤ 7.6e-4 is not exhibited. relic_yield's example "B_h =
(g_h/g*)ξ⁴ ≲ 1.2e-14" is the thermal-energy form; a direct branching of a heavy inflaton produces quanta with
E ≈ m_I/2 and Y_N ≈ 1.5 B_h T_RH/m_I [D], i.e. a strongly underpopulated bath (μ_N ≪ 0). By relic_yield's own
monotonicity result (restoration needs n ≥ n_th(T_c) = 1.8 v³ for any μ_N ≤ 0) such a bath restores later or not at
all, so the example is not self-consistent; μ_N (or Y_N) must be declared as a separate input. (c) The memory of an
inflationary φ condensate through thermal restoration is argued (P-invariant thermal state), not computed; it needs
Γ_φ > H, which holds for ξ above the kinetic threshold.

**R2 Transition/formation appropriate to the history** (resummation, nonzero-charge phase structure). **PARTIALLY
MET.** Neutral: met (resummed T_c, X/D never condense, p(+) = p(-) proof, KZ count, percolation). Charged: the T = 0
structure is exact and verified; the finite-T (T, n) diagram is leading order only (no daisy/CW, a ~10%
gas-condensate partition ambiguity per v_charge); the charged exit was argued (spinodal, p = 1/2, KZ ~1e2-1e3 ℓ)
and only now is being simulated (charged_exit_sim, 2D, preliminary: walls form in all 7 charged-exit runs with
p(+) at formation 0.493-0.510, mean 0.501).

**R3 Formation AND persistence of a macroscopic network, expansion and backreaction before any scaling/CMB import.**
**PARTIALLY MET.** Neutral unbiased (X = D = 0) network: formation and persistence inside the box (3D physical to
H^-1/δ = 72; PRS to τ = L/2; spanning wall ≥ 97.4% of area in 3D), scaling approached but not established
(μ = -0.90 to -0.975). Extrapolation to the physical hierarchy (5e6 for the reference, 1.5e8 for the charged exit)
rests on the scaling attractor (literature-supported, not proven here). Backreaction is included at the background
level in cosmo_observable with radiation-era A; there is no dynamics through domination (n = 1/2 bracket). The
subdominant CMB template is correctly not used. **Charged exit: persistence is NOT established** (see §2.1).

**R4 Observable with parameter provenance, magnitude fixed before comparison, retrospective data labelled, conditional
constraints.** **MET.** R_H = H0_pred/H0_meas is fixed by σ (action), A (literature), the n bracket, T_CMB and
ρ_Λ = 0; H0 data are labelled retrospective; everything is phrased as conditional history constraints. Caveats: no
freeze artefact (hash/timestamp) records the observable definition before comparison, unlike Codex's
PUBLIC_FREEZE; v = 1 TeV is the benchmark, not independently fixed (FEASIBILITY requires that for a "predictive
magnitude"), so the result is a conditional rejection of a history, which the module says.

**R5 No action extension; inflation/reheating labelled outside the action.** **MET.** The ξR φ² terms are inside the
EFT matching (zero at μ0, generated by running) and are used as such.

**R6 Explicit σ in GeV³; literature scale labels reconciled.** **MET.** σ = 9.428090416e7 GeV³ throughout; the LMS
"0.93/0.96 MeV" labels reconciled with σ^(1/3) via a kg/cm² reading (an interpretation, confirmed numerically by
v_cosmo, not an erratum).

### 2.1 Is the network from the CHARGED exit shown to persist?

No, not by the completed modules. charge_phase (H1) and relic_yield (§8, R3) state formation only and say
"persistence is not addressed". network_sim simulated only the neutral X = D = 0 quench. The charge_phase sign-memory
toy shows closed −v regions inside a connected +v background collapse, which is irrelevant to an unbiased p = 1/2
exit. What exists:

- **Symmetry argument (sound):** P (φ → −φ) commutes with H(t) and with the U(1)_X charge, so a P-invariant charged
  state gives p(+) = p(−) and no energetic bias (charge_phase §4, v_charge CONFIRMED). The wall-removal mechanisms in
  the literature (vacuum-energy bias, initial-population bias) are absent.
- **Expected direction of charge effects (argument, [E]):** the planar wall is a Pöschl-Teller well for X with
  l(l+1) = 2g_X²/λ_φ = 2, so l = 1 exactly at the benchmark: X is reflectionless and has one wall bound state at
  ω = (m_X² − λ_φ v²/2)^(1/2) = 70.7 GeV [C], above m* = 56.2 GeV. Charge therefore prefers Q-balls but can sit on
  walls. Lumps attached to walls add friction and pin walls, which slows coarsening (larger A, earlier domination).
  Charge trapped on a closed wall may stabilise small closed walls (extra relics), but cannot remove the infinite
  wall cluster.
- **Preliminary simulation (charged_exit_sim, in progress):** 2D z-invariant lattice, φ and complex X, radiation era,
  exit from a homogeneous branch-2 condensate at n_i = 2 n_c or 4 n_c; charge conserved to 2e-16, energy residual
  converges as dt². Walls form in 7/7 charged-exit runs, p(+) at formation 0.501 ± 0.005 (seed std), walls present at every
  output to a = 2.5-3, 5-30% of the charge in lumps with roughly 2/3 of large lumps attached to walls, and the wall
  length exceeds the neutral partner's by 1.3-3.4x at matched times (slower coarsening). Subcritical supply into a +v
  background (κ = 0.9): no walls. Forced supply above n_s1 (κ = 4): walls form with p(+) = 0.47.
- **Limits of that evidence:** 2D (Q-balls are Q-strings; 2D p = 1/2 sits at the percolation threshold), H/m_X at
  the n_c crossing 4e-3 to 1.6e-2 against the physical 4.7e-9 (r = 0.05) or 5.5e-10 (r = 1), only a factor 2-3 in a
  after formation, no NOTES.md, the 3D runs, resolution pair, de-charged controls and H → 0 continuations not yet
  finished.

**Verdict:** formation from the charged exit is supported (symmetry proof plus preliminary 2D lattice). Persistence
needs (i) completion of charged_exit_sim including 3D and the de-charged controls, and (ii) a written argument that
connects the late-time charged network to the neutral scaling network (lumps detaching or co-moving, A_charged ≥
A_neutral), plus a statement on charge-stabilised closed walls as an additional relic. Until then the charged branch
should be reported as "formation established at p = 1/2; persistence expected by symmetry and preliminary lattice
evidence, not established".

## 3. Overclaims (location, issue, correction)

1. thermal_potential §0/§5: "In every resummed scheme the perturbative transition is second order". A gap-equation
   soft-φ resummation gives a spurious weak first-order transition, φ_c/T_c = 0.026 (v_thermal). Say "in the AE,
   Parwani and DR schemes". §11 C1 "re-breaks it continuously" must carry "(argued from 3D Ising universality, not
   computed)". The packet nowhere computed a nonperturbative order; any "nonperturbatively second order" in the final
   report would be an overclaim. Formation (p = 1/2, N_H ≫ 1) is independent of the order.
2. charge_phase H1 (bold): "A random-sign domain pattern spanning the horizon forms." It is a conditional formation
   statement that rests on homogeneous super-horizon supply, undamped-memory assumptions and a mean-field KZ estimate
   (v_charge). Also "1e-9 of a Hubble time" should be (3-6)e-8, and the 248-315 "components" are F² components, not
   negative-φ patches (v_charge).
3. relic_yield §0/§8: "For r ≳ 1e-3 the hidden charge is Bose-condensed at T_c and keeps φ = 0 ... until
   T_SM = 198 r^(-1/3) TeV." Its own table gives condensed fraction 0 at r = 1e-3 (r_BEC = 1.22e-3 at ξ_max); a
   distinctly charge-driven exit needs r ≳ 4e-3-1e-2, and the exit T_SM is 2-21% above the cold formula (I5).
4. relic_yield §9: "No charge-driven restoration occurs if and only if T_Q < T_n(r)". Only "if": restoration is forced
   only above n_s1 (T_Q > 219 r^(-1/3) TeV); the window n_c-n_s1 is unresolved (I4).
5. relic_yield §7: the example hidden reheating "B_h ≲ 1.2e-14" presumes a thermal (μ_N = 0) number, inconsistent with
   a heavy-inflaton branching (R1(b)).
6. relic_yield §8 summary: "No work analyses a persistent wall network of an exact Z2 in an FLS-type formation
   history." The search covered eight papers, several read only at abstract/introduction or "searched" level (OCR
   scans). Restrict to "none of the eight papers checked".
7. hidden_heating §7: "H_inf ≳ 0.44 v restores Z2_φ during inflation through hidden-sector gravity alone." True only in
   the stochastic equilibrium (~200-300 e-folds); from a +v start over N_obs = 50-60 e-folds percolation needs
   H ≳ 1.5-1.8 TeV (v_heating). "Wrong-vacuum regions of macroscopic size are avoided only for H ≲ 184-275 GeV": the
   275 GeV end is not conservative. "Moves T bounds by 8%": T_RH,max moves by 4.7% (v_heating).
8. network_sim §6 C2: "the network disappeared in every tested case" is contradicted by its own n_seeds_gone_by_end
   (1/3 at p = 0.45 in 3D PRS and 2D physical; 2/3 at p = 0.40 in 2D physical) (v_network); "for a formation
   hierarchy H^-1/ξ_f ~ 1e13" conflates δ with ξ_f (I7) and is an extrapolation of K10 by many decades. §1 "isolated
   closed walls always collapse" holds for the neutral X = D = 0 sector only; charged closed walls were not studied.
   2D PRS "A = 0.652 ± 0.002" understates the realization scatter (0.69 ± 0.06).
9. cosmo_observable §0/§4: "at least 5.4 orders" is 5.3 with joint variation; "BBN only mildly affected" omits the
   D/H epoch (ΔN_eff-equivalent 0.05-0.28 at 70-30 keV, v_cosmo) and no BBN network was run; "residual is a common
   conversion-constant offset" is speculative.
10. inflation_spectator: H_perc 1784/1607 GeV should be 1792/1617; the small-B Λ₁ column is the restricted-equilibrium
    MFPT, mislabelled QSD, overstated up to 24%; the X entropic column is an over-correction to be used as a bound;
    "non-log de Sitter terms are O(1/ln(H/μ0))" is wrong at ln(H/μ0) ≈ 1 (v_inflation). Its classification row D
    "a macroscopic wall network forms at re-entry" is formation only; persistence is outside the module (stated
    there).

No module claims "excluded" for Candidate 1 or calls any number a prediction; that framing is consistent with the
BRIEF. The packet-level statement must likewise say "this rejects a history, not the action" and that a wall
signal would not identify Candidate 1 among Z2 sectors (FEASIBILITY).

## 4. Minimal additional calculations, ranked

1. **Finish charged_exit_sim and write its NOTES.md** (in progress; hours): 3D runs, resolution pair, de-charged
   controls, H → 0 continuation; report p(+) statistics, wall persistence and coarsening exponent against the neutral
   partners, charge fraction on walls, and whether charge stabilises closed walls. Add the persistence argument of
   §2.1. Closes R3 for the charged branch.
2. **Assemble the single reference history** (§5) with all inputs and recompute the dependent numbers consistently
   (resummed T_c, ξ, Y_N/μ_N, r, H_inf, exit temperature from n_local(T_h), H^-1/δ at formation, T_dom with
   backreaction, R_H). Bookkeeping only; my critic_results.json gives the numbers. Closes R1(a).
3. **Production self-consistency:** either exhibit an outside-action producer that gives μ_N ≈ 0 at ξ ≈ 5e-4 (e.g.
   decay of a light hidden-sector-only mediator with m ≲ few T_h) or rederive the restoration condition for an
   underpopulated bath (relic_yield's n ≥ 1.8 v³ criterion) as a function of the producer's (B_h, m_I, T_RH).
   Analytic; closes R1(b).
4. **Physical-hierarchy support for persistence:** a longer 3D physical-equation run or a VOS calibration showing that
   A is independent of δ/H^-1, so that the ~5-decade (thermal) or ~8-decade (charged) extrapolation is justified.
   Medium cost.
5. **Collapsing minority domains in inflation regimes B/C** (H_inf ≈ 140 GeV-1.6 TeV): combine inflation_spectator's
   dp_-/dn with cosmo_observable's relic bound (T_end ≳ 60 MeV × f c_b p_-) and the post-T_dom self-domination/PBH
   flag. Analytic. My rough reading [E]: with a relaxed start and f c_b ~ 1, injection continuing to late times needs
   Λ₁/2 ≲ 1e-11, i.e. H_inf ≲ 240 GeV; a fresh start suppresses late re-entry and weakens this. This is Codex's
   "quantified consequence of collapsing bubbles", currently left as a free-input formula.
6. **Spectator-condensate relic for H_inf ≫ v:** after inflation the light X, φ, D condensates (⟨U⟩ = c H⁴,
   c = 0.04-0.09) fragment into stable quanta with no number-changing. My estimate [E] with instant reheating:
   Y_N ≈ 1e-9 (H_inf/1e13 GeV)^(3/2), above the relic bound 4.3e-12 for H_inf ≳ 3e11 GeV. Needed for any branch
   that invokes high-scale inflation (critic_results.json, spectator_condensate_relic).
7. Lower priority: next-to-leading (T, n) phase diagram near the exit (shifts T_n by a few %); finite-T bounce in the
   n_* > n > n_c window (my thin-wall S_3/T ≥ 2.6e3 says irrelevant); dynamics through wall domination (even n = 2
   gives R_H ≥ 2e5); 3D lattice order of the transition (formation does not depend on it); BBN D/H shift; a freeze
   record of the observable definition.

## 5. Reference history and complementary branch

### 5.1 Reference H★: hot, relic-compatible hidden sector with negligible charge (network forms)

| Input | Value | Status |
|---|---|---|
| Action | unchanged benchmark: v = 1 TeV, g_X = 0.1, λ_φ = 0.01, λ_X = 0.001, h_D = 0.01, λ_D = 0.001, m_D0 = 100 GeV; MSbar matching at μ0 = 100 GeV, ξ_a(μ0) = 0, portals 0; ρ_Λ = 0 (base) | action / BRIEF; ρ_Λ free |
| Inflation | H_inf = 1e8 GeV, instant reheating of the SM to T_RH ≈ 8e12 GeV; φ sign memory irrelevant (erased by thermal restoration) | [P]; H_inf ≤ 1e10 GeV keeps the spectator relic ≲ 1e-2 of the bound [E] |
| Hidden bath | outside-action producer leaves a P-invariant hidden population of all four dof (φ, q1, q2, D), kinetically thermal below T_SM ≈ 1.8e8 GeV, with ξ = T_h/T_SM = 5.0e-4, μ_N = 0 (Y_N ≈ 1.3e-12), hence T_h,max ≫ T_c | [P]; thermal-form B_h equivalent 2.3e-15 (only if the products are thermal-like, R1(b)) |
| Hidden charge | r = Y_Q/Y_B = 1e-4, homogeneous, present with the bath (below r_BEC = 3.5e-4: not condensed) | [P] |
| Entropy | separate conservation in both sectors; g* = 106.75, g_h = 4; no late entropy injection | [P] |

Consequences in the modules: Z2 restored thermally; transition at T_h = T_c = 1.557-1.564 v (resummed;
perturbatively second order in AE/Parwani/DR, order not computed nonperturbatively) at T_SM = T_c/ξ = 3.1e6 GeV,
where H/m_X = 1.4e-7 and H^-1/δ = 5.2e6 [C]. Charge negligible (μ_Q/T ≈ 0.01 < M/T = 0.034). p(+) = p(−) = 1/2
exactly; KZ N_H ≈ 10^13.7 (z = 1) to 10^15.7 (z = 2) per H^-3 [interpolated from thermal_potential]; both vacua
percolate. Relic: ω_sym ≈ 0.034, ω_asym = 1.3e-4 (≤ 0.12 retrospective); maser product m_X r + 79.28 ω_sym = 2.7 GeV
< 63.93 GeV. If the network persists: T_dom = 6.37 keV, R_H = 1.0e9-2.3e9 (n = 1) or 2.2e5-3.4e5 (n = 2), so this
history is inconsistent with the known H0 by ≥ 10^5.3 (retrospective).

**Persistent network: YES according to the modules** (formation computed; persistence supported by the neutral,
unbiased lattice networks of network_sim plus the HKS/Leite-Martins scaling literature, extrapolated ~5 decades in
H^-1/δ). This history is rejected, conditionally on its inputs. It is not an exclusion of the action.

*Charged variant H★Q* (same, r = 0.05): BEC at T_c (condensed fraction ≈ 0.99); φ = 0 held by the condensate until
T_h ≈ 0.29 v, T_SM ≈ 578 TeV (cold formula 538 TeV), H/m_X = 4.7e-9, H^-1/δ = 1.5e8 [C]; spinodal exit with
p = 1/2; ω_asym = 0.067, total ≈ 0.10. **Network forms; persistence not established** (§2.1).

### 5.2 Complementary branch C★: cold hidden sector with low-scale inflation (no network)

| Input | Value | Status |
|---|---|---|
| Action | as above | |
| Inflation | H_inf = 120 GeV; N_* ≈ 48 (instant reheating); the observable patch is at φ = +v at the start of the last N_* e-folds, with N_tot − N_* ≳ N_roll + N_lag ≈ 17 + 22 ≈ 40 e-folds | [P]; pre-inflationary vacuum selection is preparation data |
| Reheating | instant, visible only, T_RH ≈ 9.2e9 GeV; direct hidden branching B_h = 0 | [P] |
| Hidden heating | gravity only (in-action): freeze-in below the formal 1.35e10 GeV and nonthermal 9.1e12 GeV thresholds; condensate scattering far below its thresholds | computed in hidden_heating |
| Hidden charge | r = 0.05, supplied homogeneously into the cold +v sector at T_Q = 100 TeV (n/n_c = 6.4e-3, well below the 538 TeV threshold) | [P] |

Consequences: during inflation H_inf = 120 GeV lies in regime A (N_reg < 1e-3 below 138-144 GeV; p_- ~ 1e-72 or
less), so no opposite-vacuum region is expected in the observable patch; no hidden restoration afterwards; the charge
sits on branch 1 and fragments into Q-balls with φ ≈ 0 cores inside a +v background (charge_phase H2; Codex's 3D
production replays at n = 0.70 n_c are this class; charged_exit_sim κ = 0.9 run, preliminary: no walls). Relic:
ω_asym ≥ 0.067, ω_sym negligible; maser product m_X r = 5 GeV < 63.93 GeV.

**Persistent network: NO according to the modules.** Caveats: the stochastic treatment is marginal at m_φ ~ H
(m_φ/H = 1.2 at 120 GeV), N_* and the +v start are preparation data, and the charge-supply mechanism is unspecified.

Between the two branches lie (i) inflation regimes B/C (rare or closed minority domains; consequences only as a
formula, gap 5 above) and (ii) H_inf ≳ 1.5-1.8 TeV or any hidden restoration, which lead back to H★-type networks.

## 6. Own numbers (critic_results.json)

- t_w = 3.882e4 s ↔ 5.822 keV; σ/(M_P²H) = 1 ↔ 7.97 keV [C].
- Charged exit (n_local(T_h) from charge_phase, linear interpolation): r = 1: T_SM,exit = 203 TeV (ξ = 5e-4), 241 TeV
  (ξ = 7.58e-4); r = 0.05: 578 / 575 TeV; r = 0.0897: 491 / 478 TeV; r = 0.64: 285 / 256 TeV [C].
- Thin-wall vacuum-bubble nucleation in the charged φ = 0 phase: S_4 = 2.5e5 at n_c (reproduces charge_phase);
  S_3/T = 6.9e3 (T_h = 0.15 v) and 2.6e3 (0.4 v) at n_c, larger at higher n [E].
- X wall bound state ω = 70.71 GeV (Pöschl-Teller l = 1, reflectionless); D: l = 0.0196, binding ≲ 0.01 GeV [C].
- Spectator condensate relic: Y_N/4.3e-12 ≈ 7e-6, 7e-3, 0.23, 230 at H_inf = 1e8, 1e10, 1e11, 1e13 GeV [E].
- Condensate/thermal energy at T_c: 0.11 (r = 0.05, ξ = 5e-4), 5.9 (r = 1, ξ = 5e-4) [D].

## 7. Not checked

The modules' internal computations (verified by the per-module verifiers); the literature statements read at
abstract level by modules (FGGK and the other seven FLS papers in particular); the charged_exit_sim runs themselves
(in progress); v_relic (not finished).
