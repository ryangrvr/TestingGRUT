# SKEPTIC of track dq_entropy (end-to-end S/Q of the M1r release, X ↔ Σ exchange to freeze-out)

2 November 2026 · Claude, adversarial skeptic in the decisive pass (BRIEF_DECISIVE.md; issue #13 context). This is
same-session AI work on one supplied model, not outside replication. Nothing here is observational evidence, a
validation, an exclusion of the real world or a new theory. Every input of M1r is a **model choice**, not a
consequence of closure. Git was read-only and nothing was committed. Labels: [D] derived, [C] computed with our own
code, [E] estimate, [Q] quoted, [A] argued, [R] recalled.

Code is in `skeptic/code/` (all own, written fresh unless stated). Results are in `skeptic/results/*.json` (total
under 100 kB). Scratch runs are in `$SP/wt2/dq_entropy_skeptic/`.

## Verdict in one paragraph

**I could not refute any dq_entropy number.** Every key result was recomputed with a different method or code:
- the frozen release (N_X, N_Σ, E_X) with an independent conformal-time symplectic integrator, run 4× further in a;
- S/Q(R), R_crit and the onset with an independent Bessel-series Hartree solver;
- the latent heat with an independent minimisation of the homogeneous functional;
- the debris ledger, the readings A/B, ν and H_i re-typed from their definitions.

All agree to the quoted precision or better.

The closed form S/Q(W) is, to 0.1–0.25%, the analytic massless three-dof ideal-gas law, S/Q = 0.19046 (W − 1)^{3/4}
[D]. The "exact universality" is a scaling identity, not an independent test.

The new issues are small:
- the frozen N_X is −8047.826 ± 0.003 (dq_entropy: −8047.84 ± 0.22);
- the Planck-condition runs use a rounded ln(M_P/M_Y) = 13.31, and the self-consistent value is 13.537;
- a computed dephasing toy reshapes the fragmentation systematic;
- dq_entropy's window means drop the last 0.25-wide bin.

With dq_heat's directly computed R(e) fed through these mappings:
- S/Q(10⁻⁶) = 4.954–4.957, free X with a margin of ×1.43;
- e = 2×10⁻⁶ is **free X** (S/Q = 3.515–3.518, a 1.3% margin), not Q-matter as in dq_entropy's reading-B band table;
- e_crit = 2.058–2.062×10⁻⁶ (bulk).

## 1. Frozen release: independent late-time integrator [C]

**Method** (`sk_late.c`, `sk_late_driver.py`, `sk_state1000.py`).
- I ran dq_entropy's own rpdq.c from its source to H_i t = 1000. It reproduces N_X(400–1000) = −8045.485140302, as stated.
- I took the full state there (a = 63.0146, η = 31.4663, K = a²H = 1.985805) and continued with a different method:
  - conformal fields χ = aZ in conformal time;
  - an exact radiation background a(η) = a₀ + K(η − η₀), with no Friedmann back-reaction (it is 5×10⁻¹¹) and no ξR
    (R ~ η² ~ 10⁻¹⁵);
  - a fourth-order Yoshida velocity-Verlet integrator, symplectic and charge-conserving, with step h = c_h/(40a);
  - the exact conformal energy |χ′|² + λ|χ|⁴ recorded alongside the engine's definition.
- I ran to η = 1300, a = 2582 (3.4×10⁸ steps), against dq_entropy's a = 630.

**Checks.**
- Step-size convergence: c_h = 0.1 against 0.2 changes the window means of N_X by ≤ 1.3×10⁻⁶.
- Gauss drift ≤ 2×10⁻⁹.
- Against rpdq on its own windows the agreement is **4×10⁻⁷ in N_X** (5×10⁻¹¹ relative) and 3×10⁻⁴ in E_X. This needed
  one correction; see N4.

| quantity | dq_entropy (rpdq, a ≤ 630, a⁻² fit) | skeptic (sk_late, a ≤ 2582, a⁻² fit over a ≥ 854) |
|---|---:|---:|
| N_X frozen | −8047.84 ± 0.22 | **−8047.8259 ± 0.0030** |
| N_Σ frozen | +1341.307 ± 0.037 | +1341.30432 ± 0.0005 |
| E_X late (conformal) | 188737.7 ± 1.1 | **188737.6639 ± 0.0006** (exact conformal energy, flat to 10⁻⁶ from a ≈ 300) |
| N_X + 6N_Σ | – | 1.4×10⁻⁵ |
| \|N_Y\| late | < 10⁻⁶ | < 1.9×10⁻⁶ (a > 300) |
| d ln std(N_X)/d ln a | −2.03 | −2.013 (a ≤ 2582), −2.028 (a ≤ 630) |
| ε, E/E_c, E_X/N_X | 0.127374, 15.6036, 23.452 | 0.1273739, 15.60364, 23.4520 |
| Y/Y_Q (M1 H_i) | 0.991023 ± 0.000027 | 0.9910212 |

Over a ≤ 630, the same a⁻² fit on my windows gives −8047.839, which is theirs. The extension to a = 2582 moves the
value by 1.8×10⁻⁶ and shrinks the error 75-fold. **Item 1 is confirmed** [C].

## 2. S/Q(R): an independent Hartree solver [C, D]

**Method** (`sk_hartree_bessel.py`). The physics model is the same, but every numerical step differs:
- Bose integrals as exponentially scaled Bessel series (`kve`), validated against quadrature to 3×10⁻¹⁶;
- entropy from s = (ε + P − μn)/T per quasiparticle species, not from occupation numbers;
- (D, P) gaps by fixed point, nested inside a two-unknown root in (T, r_c²), instead of their four-unknown hybrid root.

| R | 0 | 1 | 2 | **2.1442** | 3 | 3.6484 | 4 | 5 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| S/Q skeptic | 1.419326 | 2.449756 | 3.348887 | **3.471619** | 4.173092 | 4.679156 | 4.946010 | 5.680494 | 7.721549 |
| S/Q dq_entropy | 1.4193 | 2.4498 | 3.3489 | (R_crit) | 4.1731 | – | 4.9460 | 5.6805 | 7.7215 |
| ideal massless gas 0.19046 (W − 1)^{3/4} | 1.4228 | 2.4541 | 3.3538 | 3.4766 | 4.1784 | 4.6847 | 4.9517 | 5.6865 | 7.7283 |

Thresholds, with W_crit independent of the release:

| threshold | W_crit | R_crit | dq_entropy R_crit |
|---|---:|---:|---:|
| bulk, 3.4716 | 49.0609 | 2.144193 | 2.1442 |
| Q = 2×10⁵ | 39.3127 | 1.519459 | 1.5195 |
| binding | 32.9156 | 1.109482 | 1.1095 |

The onset at the bulk R_crit is T_h/m_X = 7.4155, n_Q = 154.099 m_X³ and a_on = 5.2301×10¹² (at H_i = 1.293×10¹¹), as
stated.

**Universality [D/C].** Solving at the same W on a different release normalisation (the M_P-portal release) gives the
same S/Q to 7×10⁻¹⁶.
- This is the dimensional identity of a conformal plasma: S/N = f(E/N^{4/3}).
- It is not an independent physical test. The "direct against spot-normalised to 10⁻¹⁵" check confirms only that the
  code is scale-covariant.
- The fitted closed form 0.190255 (W − 0.97494)^{0.75013} − 0.00428 is the analytic law
  S/Q = (4/3)(3π²/30)^{1/4}(3/4)^{3/4}λ_X^{1/4}(W − 1)^{3/4} = 0.19046 (W − 1)^{3/4}.
  - This is a three-dof massless gas plus a circular condensate holding all the charge.
  - The Hartree masses lower it by 0.1–0.25%.
  - The ideal-gas law alone gives W_crit = 48.97 (bulk).

**Item 2 is confirmed** [C].

## 3. Post-relaxation entropy [C/E]

**Latent heat** (`sk_homfunc.py`). I minimised the homogeneous T = 0 M1 functional directly at fixed n:
e(n) = min over (|X|², φ) of [n²/4|X|² + g²φ²|X|² + λ|X|⁴ + V_φ].

| quantity | skeptic | dq_entropy / fx_formation |
|---|---:|---:|
| n_b (m_X³) | 449.313 | 449.31 |
| Ω_b (m_X) | 0.7758647 | 0.775865 |
| φ/v in the bulk | 0.151 | – |
| n_sp (m_X³) | 233.208 | 233.2 |
| E/Q_hom(n_sp) (m_X) | 0.827277 | 0.82727 |
| φ_eq(n_sp)/v | 0.336 | 0.34 |
| ΔE/Q | **0.051412** | 0.0514 |
| Δ(S/Q) at R_crit, f_c ΔE/T_sp (a_sp/a_on = 0.8519, T_sp = 8.705) | **0.005526** | 0.0055 |

- The thermal-mass variant (0.036, an [E] value) was not recomputed.
- The surface bound 1.3×10⁻³ and the R_crit shift (2.1362–2.1395) re-type correctly.
- The direction is one-sided: all of these raise S/Q.

**Debris epochs** (`sk_debris_epoch.py`). The ledger uses the late (a ≥ 315) values of N_YH and N_σ, but Y_H decays at
⟨a⟩ ≈ 9.4 and σ at a ≈ 2.2.
- At a = 9.4: N_YH = 0.010376 (late value 0.010304) and N_YL = 0.56670 (late value 0.56642).
- At a = 2.2–5: N_σ = 2–7×10⁻⁵, against 2.0×10⁻⁵ late; the diagnostic is ambiguous for a ≲ 2.
- The effect on R is ≤ 3×10⁻⁴ at e = 10⁻⁶. It is negligible [C].

## 4. Readings, ledger, H_i, and the branch with the direct heat (`sk_ledger.py`) [C]

**Ledger and readings, re-typed from their definitions.**
- Ledger: m_L = 22.2486, Γ_L = 3κ²m_L/256π³ = 1.0059×10⁻⁷ (I re-derived the width: |M|² = 18κ² per channel, two
  channels), ⟨a_L⟩ = 5568.65, **R_Y = 0.378584**, R_casc + R_Σ = 0.018414 (10⁻⁶/e)^{1/2}.
- Transfer: ν = 0.95954.
- e = 10⁻⁶: S/Q_B,band = 4.6788 + 0.0055 and S/Q_A = 5.6565 + 0.0055.
- e = 2×10⁻⁶: S/Q_B,band = 3.3635 and S/Q_B,broad = 3.6131 (before increments).
- H_i ledger: free X 1.29297×10¹¹; Q-matter 1.4990–1.5245×10¹¹; Q = 2×10⁵ convention 1.35094×10¹¹.

All of these reproduce.

**Direct heat.** dq_heat's R(e) is in units of its 400–1000 E_X = 188706.29 at N_X = 8045.485, not of the frozen 188737.7.
So I used W = (E_X,frozen + R·E_X,win)/E_c(N_frozen), which lowers S/Q by 0.0011 at 10⁻⁶. The S/Q column below is the
independent Hartree value; the branch includes the item-3 increments.

| e | R_total (dq_heat) | W (frozen frame) | S/Q | branch (incl. +0.004–0.007) |
|---|---:|---:|---:|---|
| 10⁻⁶ | 4.0058 | 78.099 | 4.9499 | free X (×1.43) |
| 1.5×10⁻⁶ | 2.7957 | 59.220 | 4.0091 | free X |
| **2×10⁻⁶** | 2.1912 | 49.788 | **3.5110** | **free X, margin 1.3%** (reading B band said Q-matter, 3.370) |
| 3×10⁻⁶ | 1.5862 | 40.349 | 2.9877 | Q-matter (bulk); Q = 2×10⁵ objects do not survive |
| 5×10⁻⁶ | 1.1028 | 32.808 | 2.5466 | Q-matter |

e_crit with the direct heat:

| threshold | Hartree | with increments |
|---|---:|---:|
| bulk | 2.053×10⁻⁶ | 2.058–2.062×10⁻⁶ |
| Q = 2×10⁵ | 3.174×10⁻⁶ | 3.187–3.196×10⁻⁶ |
| binding | 4.953×10⁻⁶ | 4.98–5.00×10⁻⁶ |

These match dq_heat's 2.06, 3.20 and 4.96×10⁻⁶.

The direct W_vec is 56.31 E_c per charge at 10⁻⁶. That is 1.110× the band transfer and 0.939× the broad one, so the §8
failure condition (below 0.5) is not met.

## New issues

**N1. Frozen-value precision [C].** N_X = −8047.8259 ± 0.0030, N_Σ = 1341.3043 and E_X = 188737.6639. dq_entropy's
±0.22 is 75× conservative, and its central value is 1.8×10⁻⁶ off. This is irrelevant to the physics.

**N2. Planck condition: rounding and self-consistency [C]** (`sk_portal_summary.py`).
- dq_entropy uses c(M_Y) = −13.31/4π² exactly. dq_heat and the fx_holder skeptic use ln(2.435×10¹⁸/30H_M1) = 13.3121.
  The difference accounts for the whole 5.05×10⁻⁵ gap between dq_entropy's Planck N_X(400–1000) = −10605.624 and their
  −10606.160. I reproduced both to 10⁻⁹.
- The self-consistent value at that release's own free-X H_i is L = ln(M_P/30H_i) = 13.537.
  - It is not "13.49": that number used the superseded 1.124×10¹¹.
  - It gives N_X +0.55% (d ln N_X/dL = 2.46% per unit here, against 1.70% at c = 0), E/E_c −0.043 and R_Y ×1.008.
  - **H_i(free X, M_P) = 1.0724×10¹¹ (−0.37%)**, not −0.2%.
  - S/Q changes by about −0.01.

**N3. Fragmentation-dephasing systematic, a computed toy [C on a toy]** (`sk_dephase.py`).
- Setup: start from the rpdq state at H_i t_f = 50–80 and switch κ (the only X-number-violating coupling) off with a cos²
  ramp over 0–16 X orbits.

  | κ switch-off | N_X shift relative to frozen |
  |---|---|
  | sudden | −2.17% … −0.64% (the instantaneous dressing) |
  | 1 orbit | −1.82% … +0.28% |
  | 2 orbits | −1.12% … +0.64% |
  | 4 orbits | −0.04% … +0.58% (+0.36 … +0.55% at t_f = 60–70) |
  | 8 orbits | +0.08 … +0.15% |
  | 16 orbits | +0.02 … +0.03% |

- Readings:
  - dq_entropy's "plausible −0.5% … +0.2%" (4-orbit running means) is **not** what a 4-orbit dephasing gives.
  - Running means are not finite-time switch-off responses.
  - The overall band is −2.2% … +0.64% in N_X, i.e. H_i −0.4% … +1.5%. For t_f = 60–70 with 1–4-orbit ramps it is
    −1.1% … +0.55%.
- This is a toy. Real fragmentation keeps the local X³ source; a lattice with Y and Σ is still needed.

**N4. Window construction [C].**
- `freeze_analysis.windows()` keeps only 0.25-bins with η₁ ≤ hi + 10⁻⁹, and rpdq's bins close on the first step past
  each edge. So every "Δη = 20" window is really [lo, hi − 0.25]. Rebuilt that way, my windows match rpdq's to 4×10⁻⁷.
- The same-edge comparison differs by ±0.01 in N_X and up to ±9 in the engine E_X, because the engine E_X carries a ±400
  orbit-phase cross term.
- Harmless for the fits.

**N5. The closed form is the ideal-gas law [D].** S/Q = 0.19046 (W − 1)^{3/4} to within 0.25%. Universality is a scaling
identity, so the 10⁻¹⁵ agreement does not test the physics. What remains to test is Hartree against the true
(Bogoliubov) spectrum: the rotating-condensate quasiparticles at k ≲ M ≈ 0.8 m_X ≪ T ≈ 7.4 m_X carry about (M/T)³
~ 10⁻³ of the entropy [A].

**N6. The branch at e = 2×10⁻⁶ flips with the direct heat [C on dq_heat C].** S/Q = 3.515–3.518 is free X, against
reading B band's Q-matter 3.370.
- The margin, 0.04, is comparable to dq_heat's single-Higgs-direction uncertainty: ±0.035 in R at 2×10⁻⁶, i.e. ±0.03 in
  S/Q.
- **e = 2×10⁻⁶ is therefore marginal.** The statement "Q-matter needs e ≳ 2×10⁻⁶" should read e > 2.06×10⁻⁶ (bulk),
  > 3.19×10⁻⁶ (Q = 2×10⁵) and > 4.98×10⁻⁶ (binding).

**N7. Critic hand-off frame [C].** dq_entropy says to give R "in units of E_X = 188737.7". dq_heat's R is per
E_X,win = 188706.29 at N_X,win.
- Feeding it into W = 15.6036(1 + R) overstates S/Q by 0.0006–0.0011.
- dq_heat's own W_threshold of 49.0568 (against 49.0609) comes from its E/E_c taken from the bin-averaged ε (15.6060)
  rather than from window means (15.6071).
- Both are negligible.

**N8. H_i follows the S/Q used [C].** The free-X H_i = 1.2930×10¹¹ uses p at the reading-B S/Q of 4.685. At the direct
4.955, p ∝ (S/Q)^{≈9.5} gives p ≈ 9.7×10⁻⁴, M = 1.0019 and H_i = 1.2923×10¹¹ (−0.05%) [E on the p mapping].

## Omitted processes looked for (none found that lowers S/Q at the percent level) [A/E]

- **CP-odd Y_L → XXX against X̄X̄X̄.** All vacuum phases are removable at late times (u = μ real), so there is no
  asymmetry.
- **Condensate-stimulated Y_L → XX.** Rate ~ 4×10⁻¹¹ against Γ_L = 10⁻⁷. Even a fully asymmetric decay moves N_X by
  ≤ 2×10⁻⁴.
- **Bose enhancement of vector or Y_L products.** None: E ≫ T.
- **X-condensate Higgsing of A.** Order e²|X|² ≪ e²h².
- **Gravitational-wave loss from fragmentation.** ≲ 10⁻¹² of the hidden energy.
- **Thermalisation of the hard vector products.** First collision at Γ/H ≈ 0.3 at the vector decay (a ≈ 1.4×10⁹),
  growing ∝ a, so complete by a ~ 10¹⁰–10¹¹, long before a_on = 5.2×10¹².
- **Incomplete number-changing relaxation.** This would *lower* the pre-onset entropy below the Hartree maximum. It
  remains the assumed [Q REVIEW_TRANSPORT] input, and it is the one direction not covered by "everything omitted raises
  S/Q".

## Reproduction (one or two processes; about 4 min)

```
export OMP_NUM_THREADS=1 SK=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/dq_entropy_skeptic
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/dq_entropy/skeptic/code
python3 -B sk_hartree_bessel.py ../results/sk_hartree_bessel.json                      # 10 s
python3 -B sk_homfunc.py ../results/sk_homfunc.json                                    # 5 s
DQ_SCRATCH=$SK python3 -B sk_state1000.py base_1000 0 1000 $SK/runs/state_base_1000.json
python3 -B sk_late_driver.py $SK/runs/state_base_1000.json 1300 0.2 $SK/runs/late_c02_1300.json   # 85 s
python3 -B sk_late_driver.py $SK/runs/state_base_1000.json 320 0.1 $SK/runs/late_c01_320.json     # 11 s
python3 -B sk_late_summary.py $SK/runs/late_c02_1300.json $SK/runs/late_c01_320.json ../results/sk_late.json
for c in "planck1331 -0.33714623856987896" "planck_dqheat -0.33719832310397263" "planck_sc -0.342806"; do set -- $c
  DQ_SCRATCH=$SK python3 -B sk_state1000.py $1 $2 1000 $SK/runs/state_$1.json; done
python3 -B sk_portal_summary.py $SK/runs ../results/sk_portal.json
python3 -B sk_dephase.py ../results/sk_dephase.json 50 55 60 62 64 66 68 70 75 80     # 40 s
python3 -B sk_debris_epoch.py ../results/sk_debris_epoch.json
python3 -B sk_ledger.py ../results                                                     # 15 s
```

The rpdq.c build reads the dq_entropy source read-only, and dq_common.py is imported read-only. sk_late.c is new.
