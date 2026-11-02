# SKEPTIC of track inf_potential

3 November 2026 · Claude, adversarial pass on `tracks/inf_potential` (BRIEF_INFLATON.md; issue #13 context). This is
same-session AI work on one supplied model; nothing here is a validation, an exclusion or evidence for any sector.
**Every potential, slice, calibration point, width and ledger is a model choice.** Data comparisons are retrospective.
Labels: [D] derived, [C] computed with my own code (or the track's engine re-run by me), [E] estimate, [Q] quoted from a
reviewed track, [A] argued, [R] recalled or literature.

Code and results are in `skeptic/code/` and `skeptic/results/`. Git was read-only and nothing was committed. Every run was
in the foreground on one core, and the scratch area was `$SP/wt2/inf_potential_skeptic/`. The track's files were only read.

## 0. Verdict in brief

1. **The inflationary side holds [C, independent].** My own cosmic-time Mukhanov–Sasaki code, with its own N_* matching,
   agrees with the track on all three benchmarks: n_s to 1×10⁻⁴, r to 0.4%, α_s to 1%, A_s to 0.4%, N_* to 0.03–0.11,
   m_I exactly, and H_*/H(S1) = 1.0003 / 1.00007. The Floquet chart for T1 is reproduced exactly. The holder-induced
   Y-pair Born count for HT6t is confirmed within a factor ≈ 0.8–0.4 by a different, time-domain method.
2. **The branch conclusions do not hold [C engine, E ledger].** The track's own engine copy (rpinf.c MODE 4, run from S1)
   leaves a homogeneous Y_H remnant of N_YH = 73.8 (HT6t), 82.1 (T1) and 9.2 (E1) per H_i³. M1r's is 0.0104, so these are
   7116, 7910 and 883 times larger. It also leaves N_σ at 535–12300 times M1r's value.
   - The NOTES' smoke table and the rough re-ledger in `summary.py` omit both. The value is stored in
     `results/engine_test.json → smoke.*.NYH`.
   - Each Y_H yields one long-lived Y_L and one cascade vector A_L (the dq_heat ledger). On the benchmarks' own release,
     S/Q(e = 10⁻⁶) is 149 / 561 / 285 (HT6t / T1 / E1), and S/Q(e → ∞) is 58 / 292 / 148.
   - With the R_YY channel removed and Y_L deposits set to zero (the most favourable bookkeeping), S/Q is still 113 / 90 /
     23 at 10⁻⁶.
   - **No benchmark has a Q-matter branch at any e, under either calibration.** "What decides the branch is where the
     holder is calibrated" (digest, R6, R7) is therefore not supported.
3. **The cause is the inflaton frequency, not the potential shape [C].** I scanned M1r's own engine (MODE 0, quadratic,
   Γ_I = 0.5) over m_I. S/Q(10⁻⁶) is 9.4, 5.7, **4.97**, 5.8, 7.5, 13 and 125 at m_I = 2, 2.5, **3**, 4, 5, 6 and 8 H_i. It
   is ≥ 140 for m_I ≥ 11 and 800–1500 for m_I = 30–42.
   - A Q-matter branch (S/Q < 3.47 at e ≤ 10⁻⁵) exists only for m_I ≈ 2–4 H_i.
   - The brief's criterion ("m_I ≲ 5 or ≳ 30 H_i") and the track's design target (m_I ≲ 5, or a nearly harmonic
     30–44.5 H_i) are both insufficient.
   - No benchmark family here reaches m_I ≲ 10 H_i at H_i ≈ 1.3×10¹¹ GeV (the track's own scan_H).
4. **Secondary findings.**
   - `exit_factor.py` applies the K shift to the spectrum's k but not to the R_vec ∝ K² amplitude. Applied consistently,
     S1 gives S/Q = 3.58 for HT6t (margin ×1.03, e_crit = 1.08×10⁻⁶) instead of 4.71.
   - The Send "Q-matter" label is contradicted by the track's own R_YY(Send): on the M1r release, S/Q = 38 / 3461 / 2246.
   - ACT DR6 (P-ACT-LB, n_s = 0.9743 ± 0.0034) puts T1/E1 at −3.7σ and HT6t at −2.8σ [R].
   - The engine width Γ_I = 0.5 H_i used for every N_*, K, R_YY and smoke number is blocked for the brief's singlet-fermion
     channel by a factor 10²–10⁴ in y|φ_I|/(m_I/2).
   - Two of WORKING_THEORY's own H_i systematics put HT6t outside its "clean window".

## 1. What was recomputed and holds

**CMB [C] (`sk_cmb.py`).**
- Method: background in cosmic time; scalar modes as the flat-gauge Q with the a⁻³d(a³φ̇²/H)/dt mass term; R = HQ/φ̇;
  tensors as massless modes. Bunch–Davies start at k/aH = 250, read-out at 2×10⁻³.
- My N_* model is a cycle-averaged matter phase with constant Γ, then entropy conservation (g_* = 106.75, g_s0 = 3.909).
- Only the shape parameters were taken from bench_*.json. V0 was set by H(S1) = 1.29×10¹¹ GeV, as in the track.

| | N_* (Γ = 0.5) | A_s at track V0 | n_s | α_s | r | m_I/H_i | H_*/H(S1) |
|---|---|---|---|---|---|---|---|
| HT6t skeptic / track | 52.82 / 52.79 | 2.102 / 2.100 ×10⁻⁹ | 0.96493 / 0.96490 | −9.53 / −9.54 ×10⁻⁴ | 2.707 / 2.710 ×10⁻⁷ | 42.372 / 42.372 | 1.00027 / 1.00027 |
| T1 | 52.80 / 52.69 | 2.109 / 2.100 | 0.96185 / 0.96177 | −7.23 / −7.26 | 2.697 / 2.709 | 128.08 / 128.08 | 1.00007 / 1.00007 |
| E1 | 52.80 / 52.69 | 2.109 / 2.100 | 0.96185 / 0.96177 | −7.23 / −7.26 | 2.697 / 2.709 | 256.16 / 256.16 | 1.00007 |

- The 0.1 difference in N_* for T1/E1 comes from my matter-only reheating (it ignores ⟨w⟩ < 0 in the first
  oscillations). It also accounts for the 0.4% in A_s.
- Instant reheating gives 53.02 / 53.00 in both codes.
- The pure-hilltop tilts agree with n_s = 1 − 2(p−1)/[(p−2)N] [D].
- r = 2H²/(π²A_sM_P²) is attained, not just a floor, as the track says.

**Floquet [C] (`sk_floquet.py`, monodromy over one exact period, m_I = 1).**
- T1 at A₁ = 1.746 M: μ_max/m = 0.2487 at k = 0.34 m; then 0.208 / 0.055 / 0.015 at 0.5 / 0.2 / 0.1 A₁. This is
  identical to the track.
- HT6t: 0.160 (track 0.142). E1: 0.115 (track 0.131).
- The differences come from how the amplitude of the asymmetric orbit is defined (I use the steep-side turning point and the
  equal-energy partner). They agree to the same order. Fragmentation of T1/E1 and homogeneity of HT6t are not challenged.

**Y-pair Born rates [C] (`sk_yy_born.py`).**
- Method: first-order Bogoliubov integrals, β_k = (i/2)∫dm²/ω e^{−2i∫ω} and the LH analogue, on my own post-inflation
  background with Γ = 0.5. Every harmonic, the amplitude decay and the expansion are automatic. This is not the per-period
  FFT.
- HT6t to 12 H_i⁻¹: N_YL = 2.40 / 1.99 / 1.22 and N_YH = 1.63 / 1.37 / 0.92 for switch-on windows of 0.15 / 0.3 / 0.6
  H_i⁻¹. The track has 2.95 and 2.00.
- The rate formulas, thresholds and couplings (|c| = g_I|cos α||Σ|, |c′| = 2g_I|sin α||Σ|) are verified [D].
- Caveat: the count sits in the first ≈ 0.3 H_i⁻¹ after ε_H = 1, when the engine shows the Y condensate still releasing
  (\|Y\| = 2.5–3.7 H_i; §2). The late-vacuum masses behind the thresholds are least valid there, so R_YY(HT6t) is
  uncertain by ≳ ×2.

**Engine and inputs [C].**
- I rebuilt rpinf.c from the track source and reproduced the smoke numbers: HT6t N_X = −3768 (track −3781), E_X = 3.67×10⁴,
  ε = −0.204. MODE 0 reproduces M1r: N_YL = 0.5666, N_YH = 0.01038.
- The MODE 4 changes (potential, force, trace, Ṙ) are correct by inspection: dR/dt = 6v′I′ + 2(3H + Γ)I′² follows from
  R = 4v − I′².
- Tolerance rtol 10⁻¹² against 10⁻¹⁰ changes N_YH by < 10⁻⁴.
- G, α_I = 2.192 − π for HT6t, and μ_b = −2.584 + 13.609i check arithmetically.
- The charge sign is the same as M1r's in all three.

**CP bookkeeping [A, checked].**
- arg(G/μ_late) = α_I is the physical phase in every benchmark. In HT6t the complex μ_b is a frame artefact, not a tuning.
- In HT6t the second vacuum (x ≈ −1) has μ_late = −12.3 + 27.2i, so m_L² = 900 − √2·29.9·40 < 0 and Y condenses there.
  This is a different hidden vacuum. Within the observable patch it is irrelevant (φ_* − φ_top ≈ 0.055 M_P ≫ H/2π),
  as the track argues.

**Induced couplings [E].**
- The Coleman–Weinberg shift of the inflaton slope from φ_I-dependent Y masses is ΔV′ ≈ 2.5×10⁻⁴ H_i³, against
  V′_* ≈ 1×10⁴ H_i³.
- The holder tadpole is 2Re(GΣ*Y²) ≈ 7×10⁻³ H_i³.
- Both are negligible for the CMB numbers.

## 2. The benchmarks' own release: an omitted homogeneous Y_H / σ remnant (main refutation)

`sk_engine.py` re-runs rpinf.c MODE 4 exactly as `engine_test.py` does (S1 start, Γ_I = 0.5). Numbers are given at
a = 1 at ε_H = 1; M1r's are at its own slice.

| late (τ ≥ 2/3 τ_end) | N_X | E_X | N_YL | **N_YH** | N_σ | N_YH / M1r |
|---|---|---|---|---|---|---|
| M1r (MODE 0, m_I = 3) | −7880 (τ ≤ 60) | 1.891×10⁵ | 0.5666 | **0.01038** | 2.0×10⁻⁵ | 1 |
| HT6t | −3768 | 3.67×10⁴ | 5.23 | **73.84** | 0.0108 | **7116** |
| T1 | −5579 | 4.10×10⁴ | 11.63 | **82.09** | 0.248 | **7910** |
| E1 | −4696 | 4.02×10⁴ | 0.474 | **9.16** | 0.0191 | **883** |

**When the remnant forms.**
- In HT6t it forms at τ = 1.05–1.6 from S1, i.e. 0.05–0.6 H_i⁻¹ after ε_H = 1. At that time \|Y\| falls from 3.2 to 0.1,
  while u swings around μ at ω ≈ 30 H_i, close to m_H = 36.1.
- In T1 and E1 it forms 0.05–0.15 H_i⁻¹ after the end. That is **before** their fragmentation (0.35–0.43 H_i⁻¹), so the
  fragmentation caveat does not remove it.
- The remnant is a stable adiabatic invariant thereafter (spread 3×10⁻⁴).

**Ledger.** dq_heat's ledger [Q] counts every Y_H as one Y_L parent plus one A_L. The A_L becomes non-relativistic and
decays with the primordial vectors. N_σ adds 2 A_L. On the benchmark's own release (`sk_reledger.py`, scenario E) this
gives:

| S/Q, holder at S1 | e = 10⁻⁶ | 10⁻⁵ | e → ∞ | no R_YY, 10⁻⁶ / ∞ | no R_YY, no Y_L deposit, 10⁻⁶ / 10⁻⁵ |
|---|---|---|---|---|---|
| HT6t | **148.7** | 89.8 | 58.4 | 144.8 / 55.8 | 113.2 / 47.4 |
| T1 | 560.9 | 383.3 | 291.8 | 116.4 / 46.2 | 89.9 / 37.6 |
| E1 | 284.8 | 194.5 | 148.2 | 28.4 / 10.1 | 23.0 / 9.1 |

For comparison:
- The track's own rough re-ledger (scenario F) reproduces 19.26 / 511 / 276.
- In HT6t the hidden heat is 56 E_X without the Y_H remnant and 873 E_X with it.
- The last column is ledger-independent at e ≲ 10⁻⁵. It sets inf_reheat's prompt-radiation ledger for Y_L to zero, and
  even then the cascade A_L alone keeps S/Q ≫ 3.47.

**Consequences [C/E].**
- Under the track's own engine and the program's ledger, S/Q ≫ 3.47 in all three benchmarks at every e. There is no
  e_crit (R6 refuted).
- The release is not M1r-like, so "HT6t is the M1r-faithful benchmark" fails.
- At S/Q ~ 10², the X̄ pair fraction p is far outside its computed range: d ln p/d ln(S/Q) ≈ 9.5 near 5 [Q WORKING_THEORY
  §3]. The asymmetric-relic H_i calibration then does not apply [A].
- The calibration loop is not closed either. The benchmark's release has N_X × 0.47–0.69 and K × 0.76–0.93 relative to M1r.
  At fixed dimensionless history (Y ∝ N_X H_i^{3/2}/K^{3/2}) the yield moves by 0.70 / 0.77 / 0.65, so H_i ≈ 1.63 / 1.54 /
  1.72 ×10¹¹ GeV before any pair correction [E] (`sk_kinematics.py`).

## 3. Why: the release is non-adiabatic once m_I is comparable to the Y masses (control scan)

`sk_mi_scan.py` runs M1r's own engine (MODE 0 is byte-identical to rpdq.c). It uses the quadratic oscillator started at rest
with μ_eff = 10i, Γ_I = 0.5 and τ ≤ 40, and varies only m_I. The ledger is as in §2, with the static exit factor
(f_eff = 1). At m_I = 3 this reproduces dq_heat's M1r S/Q (4.97 against 4.951).

| m_I/H_i | 1.5 | 2 | 2.5 | **3** | 3.5 | 4 | 5 | 6 | 8 | 11 | 20 | 30 | 42.4 | 128 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N_YL + N_YH | 0.026 | 0.093 | 0.28 | **0.58** | 1.02 | 2.16 | 4.13 | 14.1 | 184 | 547 | 1417 | 2291 | 5286 | 596 |
| K | 5.59 | 2.97 | 2.28 | 1.99 | 1.84 | 1.75 | 1.66 | 1.61 | 1.57 | 1.54 | 1.52 | 1.52 | 1.52 | 1.52 |
| S/Q(10⁻⁶) | 28.0 | 9.38 | 5.67 | **4.97** | 5.12 | 5.83 | 7.51 | 13.2 | 125 | 141 | 193 | 795 | 1505 | 416 |
| S/Q(10⁻⁵) | 5.62 | 2.81 | 2.21 | 2.21 | 2.51 | 3.34 | 4.79 | 10.4 | 84 | 124 | 171 | 490 | 997 | 249 |

Reading [C engine, E ledger, A interpretation]:
- **M1r's m_I = 3 sits at the minimum of S/Q(m_I).** Below it, K grows and with it the vector heat (∝ K²). Above it, the
  Y condensate (\|Y\| = 6.9 H_i at the slice) can no longer follow the release adiabatically.
- The Q-matter branch (S/Q < 3.47 at e ≤ 10⁻⁵) exists only for m_I ≈ 2–4 H_i.
- The brief's "Y resonance" window (m_I ≲ 5 or ≳ 30) and the track's design target (§5) are not sufficient conditions.
  m_I = 5 is already free X at 10⁻⁵, and every m_I ≥ 8 gives S/Q ≥ 80.
- No single-field plateau or hilltop at H_i ≈ 1.3×10¹¹ GeV with r ≈ 3×10⁻⁷ reaches m_I ≲ 10 H_i (track scan_H). The
  M1r surrogate needs a different construction, e.g. a multi-field ending or the holder on a light spectator [A].
- This also slightly extends WORKING_THEORY §4's [E] range: m_I = 2 gives 9.4 at Γ_I = 0.5.

## 4. The exit factor: the K shift is applied to k but not to the amplitude

`exit_factor.py` shifts the M1r per-mode spectrum by k → k(K_bench/K_M1r)^{1/2}, with K_bench = K_end e^{2N_geo}. It then
multiplies M1r's R_vec(k ≤ 1) by f_eff only. R_vec ∝ K² at a fixed release, as WORKING_THEORY §4 and dq_heat state: the
number scales as K^{3/2} and the deposit as K^{1/2}.

Treating K consistently means using the K of the same physical slice as M1r's (the oscillation start, a = 1 at ε_H = 1):
K = 1.516 for HT6t and 1.854 for T1, against 1.986. This rescales R_vec by (K/K_M1r)² and the Y deposits by
(K/K_M1r)^{1/2} (scenario B):

| S1, M1r release | HT6t | T1 | E1 |
|---|---|---|---|
| track (A): S/Q(10⁻⁶), e_crit | 4.711, 1.90×10⁻⁶ | 4.596, 1.83×10⁻⁶ | 4.559, 1.81×10⁻⁶ |
| K-consistent (B) | **3.575, 1.08×10⁻⁶** | 4.270, 1.57×10⁻⁶ | 4.248, 1.56×10⁻⁶ |

On the M1r release itself, HT6t at S1 is therefore marginal (×1.03), not ×1.36. Section 2 supersedes both rows.

## 5. The Send calibration: the "Q-matter" label contradicts the track's own R_YY

R7 (S/Q = 2.16 / 2.04 / 2.01, exit factor only) is reproduced (scenario A).

Adding the track's own Send-calibrated Born R_YY (g_I × 3.3–5.7; R_YY = 74 / 3.1×10⁴ / 1.7×10⁴ in M1r units, from
`yy_harmonics.json`) on the same M1r release gives S/Q = **38.2 / 3461 / 2246** (scenario C). At S1 the same addition
gives 9.3 / 363 / 166 (D).

Both calibrations are therefore free X once the track's own new heat source is included. The claim that the calibration
point decides the branch rests only on isolating the exit factor.

## 6. CMB target has moved since Planck 2018 [R]

ACT DR6 combined with Planck (P-ACT) gives n_s = 0.9709 ± 0.0038. Adding CMB lensing and DESI DR1 BAO (P-ACT-LB) gives
0.9743 ± 0.0034 (Louis et al. arXiv:2503.14452, as quoted e.g. in arXiv:2510.18656 and 2509.23329). A DESI-DR2 combination
quoted at 0.9752 ± 0.0030 is unverified here. SPT-3G D1 (arXiv:2506.20707) combined with ACT/Planck and DESI likewise
disfavours Planck-preferred models at > 2σ (arXiv:2512.05108).

Pulls with my n_s:

| | HT6t | T1 / E1 | PL2 (track, 0.9717) |
|---|---|---|---|
| P-ACT-LB | −2.8σ | −3.7σ | −0.8σ |
| P-ACT | −1.6σ | −2.4σ | +0.2σ |

The brief's Planck 2018 target (0.9649 ± 0.0042) is met. The track's "too blue" PL2 family is the one the newer data prefer.
HT6t's β would have to be re-tuned to n_s ≈ 0.974. A longer matter era (§7) pushes every benchmark redder. These are
retrospective tensions, not exclusions.

## 7. Γ_I = 0.5 H_i is not available to the brief's reheating channel [D/E]

Γ = N_f y²m_I/8π = 0.5 H_i needs the following. The blocking factor is y\|φ_I\|/(m_I/2) at the first turning point; a_kin
is when redshift alone (A ∝ a^{−3/2}) lifts it.

| | y (N_f = 1, 3) | blocking factor | a_kin | matter e-folds | t_kin (H_i⁻¹) | K if decay at a_kin |
|---|---|---|---|---|---|---|
| HT6t | 0.54, 0.31 | 1.8×10⁴, 1.1×10⁴ | 690, 480 | 6.5, 6.2 | 1.3×10⁴, 7.3×10³ | 25, 21 (×17, ×14) |
| T1 | 0.31, 0.18 | 3.1×10³, 1.8×10³ | 210, 150 | 5.4, 5.0 | 2.1×10³, 1.2×10³ | 14, 12 |
| E1 | 0.22, 0.13 | 211, 122 | 35, 25 | 3.6, 3.2 | 142, 82 | 5.9, 4.9 |

- Thermal masses yT ≫ m_I at T_RH ~ 10¹⁴ GeV block it further (cf. Kolb–Notari–Riotto hep-ph/0307241) [R].
- In fragmented T1/E1 the oscillon amplitudes stay ~M, so the blocking persists inside them [A].
- Only non-perturbative fermion production at the zero crossings remains, which needs its own computation (Greene–Kofman
  hep-ph/9807339; Giudice–Peloso–Riotto–Tkachev hep-ph/9905242; Peloso–Sorbo hep-ph/0003045) [R].
- Consequences of a perturbative history:
  - the Higgs-plasma deadline (τ ≈ 8–10 H_i⁻¹) is missed by 10–10³;
  - N_* drops by 0.8–1.6, so n_s falls by about 0.6–1.1×10⁻³;
  - K rises ×3–17. R_vec ∝ K² at a fixed release, which raises S/Q further.
- Every Γ_I = 0.5 number in the track (N_*, K, R_YY, the smoke release) is conditional on a width the stated channel cannot
  deliver perturbatively.

## 8. Smaller points

- **HT6t's clean window is fragile [C, interpolation of scan_H].**
  - Its lower edge (m_I = 44.5 H_i) is at H_i = 1.20×10¹¹, only 7% below the free-X calibration.
  - WORKING_THEORY §5's portals-at-M_P calibration (1.08×10¹¹) and the −16% history shift both give m_I/H_i = 47.9. That
    opens the fundamental φ_I → Y_LY_L, as in T1 (R_YY ~ 10³). (§2–3 make the window moot anyway.)
- **R_YY units.** "R_YY ≈ 7–9 against M1r's budget 4.0" is in M1r E_X units. Against the benchmark's own E_X
  (3.6×10⁴ instead of 1.9×10⁵, with K = 1.52) it is ≈ 31 for HT6t (track summary D_YY/E_X).
- **The "larger homogeneous Y_L remnant" (5.2 against 0.57).** The comparison is correctly normalised, but it is the small
  part. The remnant is dominated by N_YH, which the table omits.
- **The exit-factor solver** (dq_model_copy) was not re-derived here. The S1/Send f_eff values were taken as given; §2
  makes them subdominant.

## 9. Checks of the track's results

| id | verdict | reason |
|---|---|---|
| R1 HT6t CMB | holds | independent MS code: n_s 0.96493, r 2.707×10⁻⁷, α_s −9.53×10⁻⁴, A_s 2.102×10⁻⁹, N_* 52.82, m_I 42.37 |
| R2 T1/E1 CMB | holds | n_s 0.96185, r 2.697×10⁻⁷, m_I 128.08/256.16; ACT DR6 tension −3.7σ (retrospective) |
| R3 pure hilltops / PL2 | holds (analytic) | n_s = 1 − 2(p−1)/[(p−2)N] gives 0.943/0.953/0.956 |
| R4 N_* vs history | holds as conditional | 53.02 instant, 52.80–52.82 for Γ = 0.5 (my model). The Γ = 0.5 history is kinematically blocked for the singlet-fermion channel (§7) |
| R5 s(N) | holds | my background gives the same φ_I(S1) and φ_I(end) to 6 digits |
| R6 S1 exit factor and S/Q | **refuted** | e_crit and near-threshold S/Q are artefacts of M1r's release: own release gives S/Q 149/561/285 with no e_crit (§2). On M1r's release, consistent K² gives 3.58 for HT6t (§4) |
| R7 Send → Q-matter | **refuted** | the track's own R_YY(Send) gives S/Q 38/3461/2246 even on M1r's release (§5) |
| R8 R_YY | holds (Born) | time-domain Born gives 0.4–0.8 of the track's HT6t count (window-dependent); couplings and thresholds verified. Understated ×4.5 against own E_X and subdominant to the N_YH remnant |
| R9 Floquet and fragmentation | holds | T1 μ_max/m = 0.2487 at k = 0.34 m (identical); HT6t/E1 same order |
| R10 Y resonance and clean window | **refuted** ("clean") | inside the window HT6t has N_YH ×7100 M1r and S/Q ≈ 150; MODE-0 m_I = 30–42 gives S/Q 800–1500; only m_I ≈ 2–4 keeps M1r's structure (§3) |
| R11 engine | holds, incomplete report | MODE 0 and MODE 4 reproduced; the smoke table omits N_YH (and N_σ) |
| R12 inputs | holds | G, α_I, μ_b, g_I and the S1 timing reproduced |

## 10. Not checked

- The dq_heat ledger itself (long-lived Y_L at vacuum rate). The §2 last column shows the conclusion does not depend on it
  for e ≲ 10⁻⁵.
- Finite-k Y with the new inflaton (an rpmodes.c copy).
- A lattice for T1/E1.
- The exit-factor stationary solver.
- The second-order (φφ → YY) process beyond the track's ×1.3 control.

## Code (`skeptic/code/`, own; reproduce in ≈ 3 min on one core)

| file | purpose |
|---|---|
| sk_cmb.py | independent background, MS spectra and N_* → `sk_cmb.json` |
| sk_floquet.py | monodromy Floquet chart → `sk_floquet.json` |
| sk_yy_born.py | time-domain Born Y-pair count → `sk_yy_born.json` |
| sk_engine.py | re-run of the track's rpinf.c (built in scratch) with the N_YH/N_YL/N_σ time series → `sk_engine.json` |
| sk_mi_scan.py | M1r engine m_I scan with the ledger → `sk_mi_scan.json` |
| sk_reledger.py | scenarios A–F and E0/E1 → `sk_reledger.json` |
| sk_kinematics.py | kinematic blocking, window, yield shift → `sk_kinematics.json` |
| sk_summary.py | → `sk_summary.json` |

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/inf_potential/skeptic/code; export OMP_NUM_THREADS=1
python3 -B sk_cmb.py HT6t T1 E1; python3 -B sk_engine.py M1r HT6t 60; python3 -B sk_engine.py T1 E1 40
python3 -B sk_reledger.py; python3 -B sk_mi_scan.py; python3 -B sk_yy_born.py HT6t; python3 -B sk_floquet.py T1 HT6t E1
python3 -B sk_kinematics.py; python3 -B sk_summary.py
```

Read-only imports: `tracks/dq_heat/skeptic/code/dq_model_copy.py` (the stationary preparation for the engine start),
`tracks/dq_heat/results/heat_table.json`, and the track's results files.

## Credits

Planck 2018 X arXiv:1807.06211; BICEP/Keck arXiv:2110.00483; ACT DR6 arXiv:2503.14452 (n_s values as quoted in
arXiv:2510.18656, 2509.23329); SPT-3G D1 arXiv:2506.20707; arXiv:2512.05108; Kallosh–Linde arXiv:1306.5220;
Boubekeur–Lyth hep-ph/0502047; Amin et al. arXiv:1106.3335; Lozanov–Amin arXiv:1608.01213; Kolb–Notari–Riotto
hep-ph/0307241; Greene–Kofman hep-ph/9807339; Giudice–Peloso–Riotto–Tkachev hep-ph/9905242; Peloso–Sorbo
hep-ph/0003045 [R]. Our own tracks, read-only: inf_potential (rpinf.c, bench files), dq_heat (ledger, dq_model_copy),
dq_entropy (S/Q closed form).
