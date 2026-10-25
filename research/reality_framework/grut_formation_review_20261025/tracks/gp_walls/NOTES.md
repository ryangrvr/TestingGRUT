# Track gp_walls: the proposed wall tilt at the gauge point (handoff item 4 of b5174e5)

25 October 2026. Claude, issue #13. Packets read through `git archive` into the scratchpad:
- b5174e5, `grut_gauge_point_test_20261025/`: REPORT, CLAUDE_HANDOFF, MATHEMATICS, CORRECTIONS, SOURCES, `code/point.py`, `code/screening.py`, POINT and SCREENING results.
- 7737487, `grut_z6_gauge_bridge_20261025/`: CONSTRUCTION, MATHEMATICS.
- d71cc01, `grut_portal_debris_window_20261023/`: MATHEMATICS (widths and portal), the FIXED_YIELD row, and the Q = 200000 profile file. The profile file is used only as a comparison.

All the code here is mine. I read Codex's `point.py`; I did not run it. Codex's timing formula is re-evaluated in
my own code and labelled **reproduction (formula)**. Every number is conditional on supplied fields and inputs. None
is a GRUT prediction, a validation or an exclusion. Comparisons with PTA and LVK data are retrospective.

## 1. Model and inputs [Q]

The light sector is U = λ/4 (φ²−v²)² + g²φ²|X|² + λ_X|X|⁴ + ½(m_D0² + h_D²φ²)D² + λ_D D⁴/4. The Higgs portal is
κ/2 (φ²−v²)(|H|² − v_H²/2) with κ ≲ 9.6e-10. Parameters: g = 0.1, λ_X = 0.001, λ = λ_φ = 0.03, v = 0.924089637 GeV,
Q = 200000. Derived quantities: m_X = 92.41 MeV, m_s = 226.35 MeV, m_D = 92.87 MeV. The tilt is
V_bias = Δ/4 (2 − 3x + x³), with x = φ/v and Δ = 1e-22 GeV⁴. The heavy Z6 sector (M_Y, h, μ, y, e) contains no φ, so
it does not enter here.

## 2. Algebra of the tilt [D] (`code/bias_algebra.py`, sympy)

- V′(±v) = 0 exactly, with V(+v) = 0 and V(−v) = Δ.
- Factorisation: V0 + V_bias = (x−1)²/4 [λv⁴(x+1)² + Δ(x+2)], with zero residual.
- The bracket's minimum is Δ − Δ²/(4λv⁴) > 0 if and only if Δ < 4λv⁴. The full potential is therefore ≥ 0.
- Higgs portal: minimising over |H|² replaces λ by λ_eff = λ − κ²/(4λ_H). Positivity then needs
  Δ < 4λ_eff v⁴ = 0.0875 GeV⁴, which leaves a margin of 8.8e20.
- Curvatures: m²(±v) = 2λv² ± 3Δ/(2v²). The shift is 1.76e-22 GeV².
- **Prior form.** V_bias equals Saikawa's eq. (2.23), εvφ(φ²/3 − v²) with εv⁴ = 3Δ/4, plus the constant Δ/2. The
  residual is zero.
- **Tension and profile.** V_bias(x) + V_bias(−x) = Δ, so V_bias − Δ/2 is odd. The first-order tension shift of a
  kink centred at φ = 0 therefore vanishes exactly; numerically it is 9e-26 GeV³, a quadrature edge effect. The
  second-order relative shift is about (Δ/λv⁴)² ≈ 2e-41.
  - σ = (2√2/3)√λ v³ = 0.128863 GeV³, so σ^(1/3) = 505.1 MeV. The wall width is √(2/λ)/v = 8.84 GeV⁻¹.
  - No static wall exists. The pressure accelerates the wall at Δ/σ = 7.8e-22 GeV.
- **False vacuum.** Thin-wall nucleation of +v bubbles inside −v regions has bounce exponent
  B = 27π²σ⁴/(2Δ³) = 3.7e64. False domains therefore disappear only through wall motion.

## 3. Core energy [C] (`code/core_profile.py`)

My solver treats Ω as a boundary-value parameter and imposes the charge directly. It gives Ω = 0.8611951,
E = 17064.5375 GeV, a binding margin of 1417.255 GeV and S(0) = 0.0060. It agrees with Codex's profile file to 4e-11.

The tilt adds ∫V_bias d³x = **1.7382e-17 GeV** on the core, matching Codex's 1.738e-17. Because V_bias − Δ/2 is odd,
an object in a −v domain is lighter by this amount. The core's half-depletion radius is 46.4 GeV⁻¹.

## 4. Annihilation timing [C] (`code/timing_gw.py`)

**Method.** t_ann = C_ann·A·σ/Δ with H = 1/(2t) (Saikawa eqs. 2.27–2.28). g*(T) and g*s(T) are computed
self-consistently from free e, μ and π gases.

- **Reproduction (formula).** Codex's fixed-g* grid gives a minimum of 4.2607 MeV, as Codex states. That minimum
  uses g* = 61.75 at T ≈ 4 MeV, which is not consistent: g* ≈ 10.75 there. The error is in the conservative
  direction.
- **Self-consistent results.**

  | Inputs | T_ann |
  |---|---|
  | HKS/Saikawa A = 0.8 ± 0.1, C_ann = 2–5 | 13.9–24.3 MeV |
  | Full grid, A = 0.5–2, C_ann = 1–10 | 6.60–39.3 MeV |

- Δ_min for T_ann = 1 MeV is 4.6e-25 GeV⁴ (A = 0.8, C = 5) or 2.3e-24 GeV⁴ (A = 2, C = 10). So Δ/Δ_min = 218 or 44.
- Saikawa's no-domination bound (eq. 2.29) is 1.2e-38 GeV⁴.
- [A] The 1 MeV reference is an engineering choice. The products are dark (X X̄, §6), so BBN does not set it
  directly. The wall-only deadline is much later: domination near 0.1 eV, and the CMB bound on surviving walls
  (Zel'dovich–Kobzarev–Okun).

**Bias-scaling variants [E].** Babichev et al. (2504.07902) fit τ_ann ≈ 3τ_i C(k)(λv/ε)^(1/3), with bias εχ³ and
V_bias ≈ 2εv³. They fitted ε/(λv) = 0.0035–0.05 and started at H_i = √λ v. Here ε/(λv) = 2.3e-21, which is 18
decades below the fitted range. Three ways to transfer the fit:
- B1, as normalised (their eq. 27): T_ann = 11.5–17.9 GeV. Not applicable, because formation here is much later.
- B2, their tentative dust-like reading (ρ_wall ∝ a⁻³ from formation at T_f = 40–86 GeV): 0.50–0.64 GeV.
- B3, eq. 26 taken literally with τ_i at formation: **4.9–10.6 keV**, which fails the 1 MeV target.

The authors caution against extrapolation and say the standard law may be restored at very small ε. The pass is
therefore robust to A, C_ann and g*, but not to the scaling law. The formation temperature of 40–86 GeV is taken from
the one-field history in tracks/affleck_dine; it is not established at the gauge point.

## 5. Gravitational waves [C]

Emission is ρ_gw = ε̃ G A²σ² at t_ann, with ε̃ = 0.7 ± 0.4. I redshift it with my own code and cross-check against
Saikawa eqs. 3.9 and 3.10, which agree within 1–2%.
- Ω_gw h² today: **6e-41 to 2.5e-36**, peaking at **0.75–4.7 nHz**. Central case (A = 0.8, C = 5): 1.0e-38 at 1.7 nHz.
- NANOGrav 15-yr (A = 2.4e-15 at 1/yr, γ = 13/3 [R]) corresponds to Ω h² ≈ 5e-10 at 1.7 nHz, which is 26–31 decades
  above this signal.
- The LVK O3 bound (Ω < 5.8e-9 [R]) is in a band ~10¹⁰ times higher in frequency.
- Babichev et al.'s ×10 production enhancement does not change this. There is no signal and no constraint, and the
  walls cannot account for the PTA signal (retrospective).

## 6. Energy at annihilation [C]

- **Size.** Terminal ρ = Aσ/t + Δ/2 = (6–15)e-23 GeV⁴. That is 1.4e-17 to 9.0e-15 of radiation, and ≤ 1.0e-7 of the
  calibrated ρ_DM/s. Codex gives 1.85e-7 using the full Δ, which agrees in order.
- **Where it goes.** Into s quanta.

  | s channel | Value |
  |---|---|
  | Width to X X̄ | 1.733e-5 GeV |
  | BR to X X̄ | 0.99995 |
  | BR to DD | 4.9e-5 |
  | BR to visible (κ = 9.6e-10, θ = 1.4e-11) | 8.6e-28 |
  | Lifetime | 3.80e-20 s (matches d71cc01) |
  | X momentum from s at rest | 65 MeV |

- **Leftover pairs.** Wall-made pairs, including the network's earlier radiation, are Y ≤ 2.3e-16. That is ≤ 9.5e-6
  of Codex's 1% budget (Y_pair ≤ 2.38e-11). The tilt therefore adds a stable but negligible symmetric X X̄ component.

## 7. Objects and walls: new open item [C/D/E] (`code/wall_core.py`, `timing_gw.py`)

- **Attraction [C].** The core holds φ ≈ 0, so a wall through its centre saves tension.
  - A product trial S = S_Q(r) tanh(z/w), with X frozen and so fixed charge, gives U(0) = **−1398.7 GeV**. This is a
    variational lower bound on the attraction.
  - Relaxing φ on an axisymmetric grid gives −1405 (h = 0.25) and −1409 GeV (h = 0.15).
  - The range is about 75 GeV⁻¹: U(d) = −974, −351 and −50 GeV at d = 3, 5 and 7 m_X⁻¹.
  - The well depth is about 99% of the object's binding margin of 1417 GeV. It is comparable to the wall-frame
    kinetic energy (γ−1)M = 1554 GeV at v_w = 0.4 [R for v_w].
- **Inertia [C].**
  - ρ_DM/ρ_wall = 5.5e12 at 86 GeV and 3.3e8 at T_ann.
  - The objects sweep and pin walls once the capture probability per crossing exceeds σ/(ρ_DM v_w t) ≈ 6e-13
    (formation) to 9e-9 (T_ann).
  - Each object meets about 6 walls before annihilation.
  - Capture needs an inelastic loss above ½Mv². With losses of 1%, 10% or 100% of the well depth, crossings slower
    than v = 0.04, 0.13 or 0.40 are captured.
  - Two inelastic channels have no gap [A]: the wall's ripple modes (translation Goldstone) and the bound X modes
    below.
  - Fully loaded walls are a snowplow. False domains then clear only when ρ_DM ≲ Δ, at **T ≈ 50 keV** [E, crude].
- **Wall-bound X [D].** With a Pöschl–Teller well, the X gap on the wall is μ₀ = 0.82834 m_X (finite differences
  give 0.82834; one bound mode). This is below Ω = 0.86120. Charge can move from a core into the wall at a gain of
  3.0 MeV per charge, so a core touching a wall is not stationary.
- **Packing.** At 86 GeV the cores take up about 24% of the volume (2.4% at 40 GeV). The walls form inside the
  fragmenting condensate.
- **Object preference for the false vacuum.** Each object gains 2δE = 3.5e-17 GeV in a false domain. Summed over
  objects this offsets 0.20Δ at 86 GeV, but only 1e-13 Δ at T_ann.

**Consequence.** The scaling timing assumes a free relativistic network. With this population that assumption is not
established. Object capture could delay annihilation past the 1 MeV target, concentrate objects, move charge onto
walls, and reshape dN/dQ. A coupled wall–object simulation is needed.

## 8. Explicit Z2 breaking: new channels? [D/E]

About +v, the tilt adds no tadpole. Its cubic shift is Δ/(4v³) = 3.2e-23 GeV, which is 1.1e-21 of the existing λv.
- **Self-decay.** s → ss is kinematically closed. All Z2-odd s couplings (s|X|², sD², s|H|²) already exist through
  ⟨φ⟩ = v, and the tilt changes them by about 1e-21.
- **Portal.** A one-loop induced linear portal A ≈ 1.2e-33 GeV gives mixing 1.8e-35, which is 1e-24 of the existing
  θ = 1.4e-11. The loop tadpole is 4e-4 of the tilt's linear term, which a renormalisation condition absorbs.
- **Other symmetries.** D parity, the X gauge charge and dark C are all untouched: V_bias is a real function of φ
  alone. This confirms Codex's kinetic-mixing statement for the bias.
- **Planck suppression.** v⁵/M_red = 2.8e-19 GeV⁴, so Δ corresponds to a coefficient of 3.6e-4 [E].

There is no significant new decay or mixing channel.

## 9. Verdicts

| Claim (b5174e5) | Codex | Mine | Verdict |
|---|---|---|---|
| ±v stationary, energies 0 and Δ | exact | sympy, exact | CONFIRMED [D] |
| Full potential positive | bracket min > 0 | Δ − Δ²/4λv⁴; with portal, margin 8.8e20 | CONFIRMED [D] |
| False-vacuum curvature | 2λv² − 3Δ/2v² | same | CONFIRMED [D] |
| Tension and profile unchanged (implicit) | uses σ unperturbed | first order exactly 0; ~1e-41 | CONFIRMED [D] |
| Extra core energy | 1.74e-17 GeV | 1.7382e-17 GeV | CONFIRMED [C] |
| Minimum T_ann | 4.26 MeV | 4.2607 (formula); self-consistent g*: 6.60; central 13.9–24.3 MeV | CORRECTED (conservative) |
| Exceeds 1 MeV across controls | yes | yes for t ∝ 1/Δ; Babichev-type transfers span 5 keV–18 GeV | CONFIRMED for the standard law only |
| Terminal energy ≪ DM | 1.85e-7 | ≤ 1.0e-7; to X X̄ (BR 0.99995); pairs ≤ 1e-5 of budget | CONFIRMED [C] |
| GW | not computed | Ω h² ≤ 2.5e-36 at ~nHz | NEW [C] |
| Annihilation unaffected by objects | not addressed | attraction 1.40 TeV; pinning threshold 1e-12–1e-8; X wall gap < Ω | NOT CHECKED dynamically; open |
| Bias keeps C and kinetic mixing at zero | yes | yes; no new channel | CONFIRMED [D/E] |

## 10. Limits

- Tree level and classical throughout.
- The wall–object numbers freeze X and use a planar static wall. There is no collision simulation, capture
  probability or wall-velocity distribution.
- The g* model uses free gases and neglects the QCD transition (T < 150 MeV is used).
- The Babichev transfers are readings of fits made far outside their range.
- The formation temperature comes from the one-field AD history.
- The GW estimate is a standard scaling estimate.

## 11. Credit and read levels

- **Read.** Saikawa, arXiv:1703.02576: eqs. 2.18, 2.19, 2.23, 2.24, 2.27–2.29, 3.5, 3.9 and 3.10 read in HTML [Q].
- **Read in summary.** Babichev, Dankovsky, Gorbunov, Ramazanov and Vikman, arXiv:2504.07902: abstract and eqs.
  25–32 read through an HTML summary, including their extrapolation caveat [Q].
- **Abstract only.**
  - Hiramatsu, Kawasaki and Saikawa, JCAP 02 (2014) 031, arXiv:1309.5001.
  - Kitajima et al., arXiv:2306.17146.
- **Recalled [R].**
  - HKS, JCAP 05 (2010) 032; Kawasaki–Saikawa, JCAP 09 (2011) 008.
  - Gelmini–Gleiser–Kolb, PRD 39 (1989) 1558; Larsson–Sarkar–White, PRD 55 (1997) 5129.
  - Zel'dovich–Kobzarev–Okun (1974); Vilenkin, PRD 23 (1981) 852.
  - Coleman, PRD 15 (1977) 2929 (thin wall); Pöschl–Teller well.
  - Dvali–Liu–Vachaspati, PRL 80 (1998) 2281 (walls sweeping solitons).
  - Q-ball formation: Kusenko–Shaposhnikov; Enqvist–McDonald; Kasuya–Kawasaki.
  - NANOGrav 15-yr, ApJL 951 (2023) L8; LVK O3 isotropic, PRD 104 (2021) 022004.

## 12. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/gp_walls
S=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/gp_walls
python3 code/bias_algebra.py results/bias_algebra.json
python3 code/core_profile.py results/core_profile.json $S/debris/.../data/profile_Q200000_lambda0.03_refined.csv $S/core_profile.npz
python3 code/wall_core.py $S/core_profile.npz results/wall_core.json      # ~90 s
python3 code/timing_gw.py results/timing_gw.json
```

The CSV argument is optional and is used only for the comparison. Extract it with `git archive d71cc01
research/reality_framework/grut_portal_debris_window_20261023`.
