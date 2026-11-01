# m1_ledger: the M1 mass-normalisation chain and the frozen mass test against computed populations

1 November 2026. Claude review track for issue #13 (brief: [BRIEF_M1.md](../../BRIEF_M1.md), key `m1_ledger`). Codex's
packets were read through `git show` only (d4747cb grut_matched_candidate_20261027, 1fcc6b7 grut_five_outputs_20261027,
ed73791 grut_m1_dynamics_20261027, bf8ab3b grut_m1_scattering_20261028, c481530 grut_m1_kinetics_20261028, b5174e5
grut_gauge_point_test_20261025, d71cc01 grut_portal_debris_window_20261023). Same-session AI work on one supplied
model; nothing here validates, excludes or predicts. No commits; nothing written outside this track and the scratchpad.

**Independent vs reproduction.** `code/frozen_test.py` (algebra, own), `code/scale_chain.py` Part A (reconstruction
from Codex's quoted JSON numbers, [D]) and Part B (runs of **our own** reheating engine `tracks/rp_reheating/code/rp.py` +
`rpreheat.c`, [C]), `code/born_ratio.py` (own derivation and own quadrature of the vacuum-Born pair channel on both
profiles, [C]; the coherent channel, 6e-5 of the total, is taken from Codex's numbers [Q]). No Codex code was run.
One thread, foreground; total compute about 7 s.

## 1. The frozen test [D] — CONFIRMED

With the three-reservoir ledger ρ_dark/s = Y_Q[η ε_Q + (1−η) m_X + 2 m_X p] and the calibration Y_Q ε_Q = D:

(ρ/s)/D = η + (1−η) m_X/ε_Q + 2p m_X/ε_Q = **1 + (m_X/ε_Q − 1)(1−η) + (2 m_X/ε_Q) p**.

| Quantity | Ours | Codex |
|---|---:|---:|
| m_X = g v | 0.0924089637 GeV | same |
| ε_Q = 17324.83235/200000 | 0.0866241617 GeV = 0.93739999 m_X | 0.93739999 |
| Y_Q = D/ε_Q | 5.079414233e-9 (Y_Q ε_Q = 4.40000e-10 GeV) | 5.07941423e-9 |
| (1 − ε/m)/(ε/m) | 0.06678047 | 0.0667805 |
| 2/(ε/m) | 2.13356093 | 2.13356 |
| η_min at p = 0 | 0.8502556 | 0.85026 |
| p_max at η = 1 | 0.0046870 | 0.004687 |

The charge constraint N₊/Q_net = 1 − η + p ≥ 0 allows η > 1 when p > 0 (ed73791's remark holds).

**The old target** (what "old Y_Q = 4.8294e-9" meant). D/Y_Q,old = 0.0911086 GeV = **0.98592 m_X** is not an object E/Q: it
is the v-calibrated mass per net charge of the fixed-yield witness (d71cc01 FIXED_YIELD_RESULTS: late mass per object
18221.72 GeV against the object's 17064.54 GeV). The old object (E/Q = 0.923316 m_X) therefore delivered 0.936494 D at
Y_Q,old (Codex 1fcc6b7: 0.936494) and needed ×1.0678122 (Codex 1.0678122). The same number fixes the kin_solve
convention: f = 2 m_X P/(D/Y_Q,old) = 2.02855 P/n_Q (REVIEW_KINETIC's 2.0286).

## 2. The scale chain [D, C] — CONFIRMED, with the reading corrected

The brief's objection (5.07941423/4.8294 = 1.0518 would need +3.42%, not +0.70%) assumes the old heavy point delivered
Y_Q,old. It did not. Codex's own MATCHED_RESULTS rows and our engine agree on the chain:

1. Supplied Y_Q,old = 4.8294009e-9 with the v-calibrated 0.98592 m_X per charge (old target).
2. Old object 0.923316 m_X → needs 5.156893e-9 (×1.0678).
3. The (3, 0.5) history at H_old = 1.3335e11 GeV actually yields **1.040807 Y_Q,old = 5.026476e-9** (τ = 800, late window);
   (5, 1): 1.042251 Y_Q,old = 5.033450e-9. Our engine (own code, τ = 800): 5.0264755160e-9 (rel. −8.6e-12) and
   5.0334502218e-9 (−4.8e-11).
4. Old-object repair: λ = (5.156893/5.026476)^{2/3} = **1.0172235 (+1.722%)**; (5, 1): 1.016277 (+1.628%). These are the
   "nine homogeneous controls … 1.722%/1.628%" of d4747cb's REPORT: the old-profile targets, superseded by M1.
5. M1 object 0.9374 m_X → target 5.079414e-9; λ = (5.079414/5.026476)^{2/3} = **1.00700905 (+0.7009%)**; (5, 1):
   1.0060786 (+0.6079%). From our τ = 800 run: 1.0070090496 and 1.0060785791 (agree to 6e-12 and 3e-11); our engine at
   Codex's λ gives Y/target = 1 − 7e-12.
6. So 1.0518 = 1.040807 (what the old point already gave) × 1.010532 (= 1.00700905^{3/2}). No exponent other than 3/2
   and no hidden 1.7% repair is involved.

**Exponent.** Y = N_X H_i³/(a³ s) with s ∝ T³, ρ_r = r M_red² H_i² ⇒ T ∝ (H_i M_red)^{1/2} ⇒ **Y ∝ H_i^{3/2}** at fixed
dimensionless history; the history depends on the scale only through η² = (H_i/M_red)², i.e. the hidden sector's
back-reaction on the constraint. Own-engine measurement (all heavy scales co-scaled, η² scaled, τ = 800, rtol 1e-10):

| history | λ | Y(λ)/Y(1) | d ln Y / d ln H_i | N_X, a, ρ_r a⁴ ratios |
|---|---:|---:|---:|---|
| (3, 0.5) | 1.00700905 / 1.0172235 / 1.05 | 1.0105320 / 1.0259462 / 1.0759298 | 1.50000000017 / 1.50000000020 / 1.50000000020 | 1 ± 1e-11 |
| (5, 1) | same | 1.0105320 / 1.0259462 / 1.0759298 | 1.50000000020 / 1.50000000020 / 1.50000000019 | 1 ± 1e-11 |

Codex's rows imply the same: exponent 1.5000000002 (Y) and 0.5000000004 (T_vis) for all three matched rows. The
"actual ODE rerun" is therefore a unit conversion of an unchanged dimensionless state (hidden back-reaction ≤ 1e-11 on
N_X); it cannot fail, and it did not.

**Our earlier 1.0454.** `rp_reheating` reported 1.0454206 Y_Q for yr_b3_g0.5 at τ = 80 (last-20% window τ 64–80);
at τ = 800 the same engine gives 1.0408073 (Codex 1.0408073). The ratio 1.004432 equals N_X(64–80)/N_X(640–800): it is
the homogeneous charge estimator's clock dependence (Codex's charge_clock: −8147.1 at τ = 80, −8116.1 at τ = 800), not
an engine difference. CORRECTED: the comparable number is 1.0408 at τ = 800; the 0.44% clock drift is the largest
uncertainty in the whole scale arithmetic (the half-window change is 2.6e-4).

## 3. The frozen test against the populations we have [E, A]

p in the frozen test is the **leftover** pair ratio Y_X̄/Y_Q after annihilation/conversion; kin_solve's f is the pair
mass over the old budget per net charge, so p = f/2.02855 and the M1 term is 2.13356 p = **1.05177 f**. Rows
(`results/frozen_test.json`; term = 0.0667805(1−η) + 2.13356 p, allowance 0.01):

| Scenario (source) | η | p | term | vs 0.01 |
|---|---:|---:|---:|---|
| pt_relax equilibrium X̄/charge at the pilot energy, no annihilation [E] | 1 | 0.069–0.110 | 0.147–0.235 | fails ×15–23 (needs 4–7% survival) |
| same with e = 1e-7 vector heat thermalised [E] | 1 | 1.9–2.9 | 4.05–6.19 | ×405–619 (needs 0.16–0.25% survival) |
| classical plateau 1.75 (relaxes at ~1e-3 per conformal unit over ~1e11–1e12 units) [C] | 1 | 1.75 | 3.73 | ×373 |
| kin_solve gas-only, e = 1e-6 band R = 4.97 (f 0.0045855) / R = 5 (0.0047584) [C] | 1 | 0.00226 / 0.00235 | 0.00482 / 0.00500 | **passes**; room (1−η) ≤ 0.077 (η ≥ 0.923) |
| kin_solve gas-only R = 7 (f 0.0288) [C] | 1 | 0.01422 | 0.0303 | ×3.0 |
| kin_solve e = 1e-7 band R = 41.7 (f 1.492) [C] | 1 | 0.736 | 1.569 | ×157 |
| cold dilute masks 30–32% read as η (Codex: not η) [C] | 0.30–0.32 | 0 | 0.0454–0.0467 | ×4.5–4.7 |
| dense primary masks 49.4% [C] | 0.494 | 0 | 0.0338 | ×3.4 |
| hot elliptic control 0.6% (64³) / 3.45% (96³), not converged [C] | 0.006 / 0.0345 | 0 | 0.0664 / 0.0645 | ×6.6 / 6.4 |
| pulses, measured outgoing X̄/Q₀ as raw p (5/10/20/40%) [C] | 1 | 0.00209 / 0.00845 / 0.0345 / 0.146 | 0.0045 / 0.0180 / 0.0735 / 0.311 | 5% passes; ×1.8 / 7.4 / 31 |
| Codex's conditional cold completions (all X̄ converted) frozen formula [C] | 0.99782 / 0.99105 / 0.96087 / 0.79519 | 0 | 0.00015 / 0.00060 / 0.00261 / 0.01368 | first three pass; 40% ×1.4 |
| same with the re-solved E(Q_f)/Q_f (ours = Codex) | | | 0.000275 / 0.001127 / 0.004903 / 0.02482 | frozen formula understates ×1.9 |
| 10% pulse kinetics τ = 1000 (Q_core/Q₀ 0.94186, 0.559 X̄) [C] | 0.94186 | 2.8e-6 | 0.00389 (frozen) / 0.00727 (Codex, re-solved E/Q) | passes only if the gas's 1.3475% kinetic energy leaves (2.07% total otherwise) |

The kin_solve crossing: f = 0.01 at R = 5.82 (log-linear between R = 5 and 7); the M1 allowance f ≤ 0.009508 moves it to
R = 5.77. The frozen formula's fixed ε_Q is itself a systematic: at the reduced Q of the completions E/Q rises
(0.93752–0.95053) and the true excess is 1.9× the frozen one; our re-evaluation with Codex's re-solved E/Q reproduces
its 0.0275/0.1127/0.4903/2.482% to all digits.

**Verdict on item 3.** No computed formation-type population satisfies the test. The passes are all conditional or
partial: the gas-only forced-thermal moments at R ≲ 5.8 (no objects, η assumed ≥ 0.92), the 5% pulse read raw, and
Codex's idealised completions (every X̄ converted, free X never recaptured, heat gone). Everything with a measured
X̄ inventory before annihilation (pt_relax, plateau, e = 1e-7) fails by 15–620×, and every clumping fraction read as η
fails by 3.4–6.6×. Two calculations would decide it: **η** needs a formation ensemble — expanding 3D evolution from the
matched release spectrum with the deformed cubic, returning the charge-weighted object spectrum y_Q and a core/tail
membership rule (masks of 30–52% are not it); **p** needs chemical relaxation of the X̄ inventory (2↔4 and
condensate-assisted channels, pt_relax's classical drift) followed by annihilation/conversion/capture with the actual
heat R (∝ 1/e through the vector daughters) — kin_solve covers the gas-only limb, c481530 one fixed-profile pulse.

## 4. Charge-loss normalisation Γt₀ = 6.27e-17 [D, C] — CONFIRMED

The light operator is A X⁶ + h.c. (dimension 6) with A = κ²μ(h/√2)/(M_Y⁴(1−r²)), the remnant of the dimension-7
κ²μ X⁶Σ*/M_Y⁴ after Σ → h/√2 (2|A|/B = r = 0.45 is scale-invariant). Under the coordinated change κ fixed, μ, h, M_Y ∝ λ:
**A ∝ λ⁻²** → A_new = 2.10886e-31/1.00700905² = 2.0796057e-31 GeV⁻² (Codex 2.0796055e-31; relative change 0.9861279 both).

Own Born derivation (X⁶ → 15 X₀⁴χ², two X̄ from four condensate quanta, W = 4Ω, two-body phase space (π/2)β,
identical-particle ½): Γ = (900/32π³)|A|²(v/√2)⁸ m_X⁻³ ∫₀^{q_max} q²|F̃⁴(q)|²β dq — the same prefactor 225/(8π³) as Codex,
evaluated with a different substitution (q = q_max(1−u²), Gauss–Legendre) on a different radial grid (dr = 0.01, Simpson):

| | old profile (Ω 0.861195, W 3.44478) | M1 profile (Ω 0.881720, W 3.52688) |
|---|---:|---:|
| dimensionless integral I (ours / Codex) | 3149626.59 / — | 2868466.51 / 2868466.51 (4e-12) |
| Born fractional core rate at c₆ = 1 (s⁻¹) | 1.039970e-46 / 1.039970e-46 | 9.471343e-47 / 9.471343e-47 |
| coherent channel [Q Codex] | 3.92e-50 (0.038%) | 6.10e-51 (0.0064%) |
| c₆ = A M_red² | 1.250391e6 | 1.233045e6 |
| **Γt₀** (included × c₆² × 4.3549e17 s) | **7.0837e-17** (our review's 7.08e-17) | **6.27163e-17** (Codex 6.27163e-17) |

Ratio new/old = 0.88536 = **0.91045 (profile: the F⁴ source shrinks and W rises) × 0.97245 (coupling λ⁻⁴)**; from the
profile change alone 7.0837e-17 × 0.91045 = 6.4493e-17, and the coordinated scale brings it to 6.2716e-17. The "modified
neutral response" (the −1.2(s−1) term in the s-channel potential) enters only the coherent channel, 6e-5 of the total,
so it cannot move Γt₀ at the quoted precision. Limits: vacuum Born, cold ground profile, included channels only — as Codex
states; a hot or excited core (ed73791's 0.2–15% X̄ emission per pulse) is a different, far larger loss channel.

## Verdict table

| # | Claim | Label | Codex | Ours | Verdict |
|---|---|---|---|---|---|
| 1 | Frozen test coefficients and endpoints | [D] | 0.0667805, 2.13356; η ≥ 0.85026, p ≤ 0.004687 | 0.06678047, 2.13356093; 0.8502556, 0.0046870 | CONFIRMED |
| 2 | Y_Q = D/ε_Q with ε_Q = 0.9374 m_X | [D] | 5.07941423e-9 | 5.079414233e-9; Y_Q ε_Q = 4.4000e-10 GeV | CONFIRMED |
| 3 | Old Y_Q 4.8294e-9 ↔ old objects give 0.936494 D | [D] | 0.936494; ×1.0678 | D/Y_Q,old = 0.98592 m_X (v-calibrated, not an E/Q); 0.9364942; ×1.0678122 | CONFIRMED (origin clarified) |
| 4 | +0.7009% coordinated scale from the old point gives Y_Q,M1 | [C] | 1.00700905 | 1.0070090496 from our τ = 800 run; the old point gives 1.040807 Y_Q,old, so 1.0518 = 1.040807 × 1.00700905^{3/2} | CONFIRMED (brief's reading corrected) |
| 5 | 1.722%/1.628% repairs; 0.6079% for (5, 1) | [C] | 1.0172235, 1.016277, 1.0060786 | 1.0172235042, 1.016277, 1.0060785791 (old-object targets, superseded) | CONFIRMED |
| 6 | Y ∝ H_i^{3/2} on a fixed dimensionless history | [D, C] | used as root estimate, verified by rerun | exponent 1.5000000002 (own engine, λ = 1.007–1.05, both histories); dimensionless state unchanged to 1e-11 | CONFIRMED |
| 7 | Baseline yields at the old point | [C] | 5.026475516e-9 (3, 0.5); 5.033450222e-9 (5, 1) | 5.026475516e-9 (−9e-12); 5.033450222e-9 (−5e-11), own engine τ = 800 | CONFIRMED |
| 8 | Our rp_reheating 1.0454 Y_Q for (3, 0.5) | [C] | 1.0408 | 1.0454 at τ = 80 → 1.0408 at τ = 800; ratio 1.00443 = charge-clock drift | CORRECTED (ours) |
| 9 | Frozen test vs computed populations | [E, A] | η, p "remain to compute" | no formation-type population passes; fails ×1.8–619; conditional passes only (kin R ≲ 5.8, 5% pulse, idealised completions) | NEW |
| 10 | Frozen formula with fixed ε_Q vs re-solved E(Q) | [D] | 0.0275/0.1127/0.4903/2.482% | frozen 0.0146/0.060/0.261/1.368%; with Codex's E(Q_f) 0.0275/0.1127/0.4903/2.482% | CONFIRMED; frozen form understates ×1.9 |
| 11 | A → 0.98613 A under the scale change | [D] | 2.0796055e-31 | 2.0796057e-31 (λ⁻²) | CONFIRMED |
| 12 | Born pair channel on the M1 profile | [C] | I = 2868466.51; 9.4713e-47 s⁻¹ | 2868466.51 (4e-12); 9.4713e-47; profile factor 0.91073 | CONFIRMED |
| 13 | Γt₀ = 6.27e-17 (from 7.08e-17) | [C] | 6.2716e-17 | 6.27163e-17 = 7.0837e-17 × 0.91045 × 0.97245 | CONFIRMED |

## What remains uncertified

- The scale chain is arithmetic on a homogeneous charge estimator whose late value drifts 0.44% between τ = 80 and 800
  and 0.026% across the last half-window; the ±0.7% calibration is smaller than that systematic's natural scale.
- The 3/2 exponent assumes the hidden sector stays a spectator (true to 1e-11 here) and g* = 106.75 at T_vis ~ 7.7e12 GeV.
- The frozen test's 1% is a program benchmark; η and p are not computed by any packet; the Q distribution changes ε_Q.
- Γt₀ is the vacuum-Born cold-ground-state channel only; excited cores emit X̄ at 10⁻³–10⁻¹ of Q per pulse.

Reproduce: `cd code; OMP_NUM_THREADS=1 python3 -B frozen_test.py; python3 -B born_ratio.py; python3 -B scale_chain.py`
(the last builds `rpreheat.c` into the scratchpad; profiles are read from the scratchpad copies saved by `git show`
d4747cb:…/PROFILE_M1.csv (sha c08176e9…) and d71cc01:…/profile_Q200000_lambda0.03_refined.csv (sha 7a737d67…)).
