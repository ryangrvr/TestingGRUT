# Track dq_critic: the decisive tracks combined into verdicts

2 November 2026 · Claude, decisive pass (BRIEF_DECISIVE.md of this directory; issue #13 context). This is same-session AI
work on one supplied model (M1r), not outside replication. Nothing here is observational evidence, a validation, an
exclusion of the real world or a new theory. **Every input is a model choice, not a consequence of closure.** That
covers M1r, e, the inflaton history, the portal renormalisation condition, ξ_H and the kinetic models. Data comparisons
are retrospective. Git was read-only and nothing was committed. At most one worker process ran at a time, in the foreground. Labels: [D] derived, [C] computed with our own
code, [E] estimate, [Q] quoted from a reviewed track, [A] argued, [R] recalled or literature.

Inputs, each read together with its SKEPTIC.md:
- dq_heat: R(e);
- dq_entropy: S/Q(R) and the frozen release;
- qm_micro: relic microphysics;
- qm_transfer: the P(k) map;
- higgs_reheat: Higgs survival.

Where a skeptic corrected or refuted a number, the corrected number is used. Every contradiction between tracks was
settled by a quick check of my own (§6).

## 0. Verdict in one paragraph

**On the static slice that M1r specifies, the declared e = 10⁻⁶ is free X**: S/Q = 4.955, outer band 4.915–4.994, a
margin of ×1.43 [C].
- The switch to Q-matter is at e_crit = 2.06×10⁻⁶ (bulk), 3.19×10⁻⁶ (Q = 2×10⁵) and 4.99×10⁻⁶ (binding objects).
- If the induced portals vanish at M_P instead, the bulk switch moves to 1.46×10⁻⁶.

**The branch at the declared coupling is not decided by M1r as specified.**
- The vector heat is set by the holder value when k_* left the horizon, 3–4 e-folds before the slice.
- Free X at 10⁻⁶ needs f_eff ≥ 0.469 for that exit factor.
- Every rolling-inflaton embedding tried gives f_eff = 0.35–0.44, hence S/Q = 3.11–3.37: Q-matter.
- The supplied decay history also enters through K² (S/Q(10⁻⁶) = 4.34–8.03 over K = 1.75–2.97) [E].
- WORKING_THEORY's "robust in direction" claim is therefore false.

**The Q-matter branch has no testable P(k) signature at any e** (§2).
- At the benchmark e = 3×10⁻⁶ (S/Q = 2.993), P/P_ΛCDM − 1 is between −4×10⁻⁵ and −10⁻⁶ at k = 100/Mpc, and smaller at
  lower k. That is 110–210× in ζ below a 20 keV-class reach, and far from current bounds.
- What is distinctive is the relic itself: compact objects of 0.2–7 Planck masses, ~0.3 per m² per year at Earth. Only
  proposed gravitational direct detection could see them.

**H_i and r** [C on Q/E inputs]:

| Branch | H_i (GeV) | r |
|---|---|---|
| free X | 1.290–1.292×10¹¹ | 2.7×10⁻⁷ |
| Q-matter | 1.48–1.50×10¹¹ | 3.6×10⁻⁷ |
| M_P portals, free X | 1.072×10¹¹ | 1.9×10⁻⁷ |

The r values are floors if H at CMB exit exceeds H_i.

**Higgs: ξ_EW = 0.1 is inside the allowed window only below a critical top mass** [C runs; E criterion].
- Inflation is safe at every corner.
- Gravitational-only reheating survives for M_t ≤ 172.84 GeV (free-X H_i), 172.79 (Q-matter H_i) or 172.89 (M_P
  portals), each ±0.09 GeV from the criterion's conventions.
- Against PDG 2024 (M_t = 172.57 ± 0.29 GeV) these are +0.9σ, +0.8σ and +1.1σ.
- Above that top mass a Higgs-coupled plasma is needed by τ ≈ 8–10 H_i⁻¹, and M1r supplies none.

## 1. Item 1: S/Q(e) of the M1r release with the computed heat [C]

**Method** (`code/sq_branch.py`).
- dq_heat's R(e) components are interpolated log-log in e:
  - piece A (frozen-angle vectors, ∝ e^−1.007);
  - the release-era band plus the separate-universe angle, and the cascade (∝ e^−1/2);
  - R_Y = 0.3792.
- Frame. dq_heat's R is per its window E_X (188706.29). The plasma is set by the frozen release (dq_entropy SKEPTIC:
  N_X = 8047.8259, E_X = 188737.6639). So W = (E_X,f + R·E_X,w)/E_c(N_X,f), with E_c = (N/c)^{4/3} and c = 6.977568.
- S/Q(W) comes from dq_entropy's closed form, checked by a direct solve with an unchanged copy of its Hartree solver
  (`code/hartree_copy.py`): 4.94991 / 3.51097 / 2.98766 at e = 1 / 2 / 3×10⁻⁶, equal to 5 digits.
- Item-3 entropy is added: +0.004–0.007.
- Bands on R:
  - computed: the multi-field check of the release-era band, |ΔR| ≤ 0.003 (10⁻⁶/e)^{1/2};
  - outer: dq_heat's full-vs-noX difference, ±0.05 (10⁻⁶/e)^{1/2};
  - plus the omitted channels from dq_heat SKEPTIC N4 (−1.7×10⁻³ … +2.6×10⁻⁴).

| e | R_total | W | S/Q (outer band) | branch | S/Q, M_P portals | S/Q on rolling paths (fixed-H drift / quadratic / plateau) |
|---|---:|---:|---|---|---:|---|
| 10⁻⁷ | 36.72 | 588.5 | 22.72 (22.65–22.80) | free X | 17.09 | 9.66 / 12.49 / 10.04 |
| **10⁻⁶** | 4.006 | 78.10 | **4.955 (4.915–4.994)** | **free X** | 4.089 | **3.24 / 3.37 / 3.11 (Q-matter)** |
| 1.5×10⁻⁶ | 2.796 | 59.22 | 4.015 (3.979–4.049) | free X | 3.434 (Q-matter) | 2.85 / 2.91 / 2.74 |
| 2×10⁻⁶ | 2.191 | 49.79 | 3.516 (3.484–3.548) | free X, margin 0.4–2.2% | 3.091 | 2.64 / 2.66 / 2.55 |
| 2.5×10⁻⁶ | 1.828 | 44.12 | 3.205 | Q-matter (bulk) | 2.878 | 2.51 / 2.51 / 2.42 |
| **3×10⁻⁶** | 1.586 | 40.35 | **2.993 (2.965–3.020)** | Q-matter; Q = 2×10⁵ objects do not survive | 2.734 | 2.42 / 2.41 / 2.34 |
| 5×10⁻⁶ | 1.103 | 32.81 | 2.552 | Q-matter; binding objects survive above 4.99×10⁻⁶ | 2.436 | 2.22 / 2.19 / 2.16 |
| 10⁻⁵ | 0.741 | 27.16 | 2.204 | Q-matter | 2.205 | 2.05 / 2.03 / 2.01 |
| 10⁻⁴ | 0.420 | 22.15 | 1.880, near the floor 1.838 (e → ∞) | Q-matter | 1.990 | 1.87 |

**Critical couplings** (S/Q = threshold; `results/sq_branch.json`):

| threshold | central | computed band | outer band | M_P portals | fixed-H drift / quadratic / plateau |
|---|---|---|---|---|---|
| bulk (3.4716) | **2.060×10⁻⁶** | 2.054–2.065 | 2.016–2.103 | **1.46×10⁻⁶** | 0.82 / 0.93 / 0.75 ×10⁻⁶ |
| Q = 2×10⁵ (2.9284) | **3.19×10⁻⁶** | 3.177–3.201 | 3.105–3.276 | 2.36×10⁻⁶ | 1.37 / 1.47 / 1.20 ×10⁻⁶ |
| binding (2.5531) | 4.99×10⁻⁶ | 4.96–5.02 | 4.82–5.16 | 3.97×10⁻⁶ | 2.32 / 2.34 / 1.97 ×10⁻⁶ |

**Notes on the critical couplings.**
- **Exit factor.** At e = 10⁻⁶ free X needs f_eff ≥ **0.469** (dq_heat SKEPTIC: 0.470).
- **Rolling-path values.** My path e_crit are 5–8% below the skeptic's (0.86 / 0.95 / 0.81). The skeptic interpolated
  ln W linearly in ln e across the coarse 10⁻⁷–10⁻⁶ interval; I interpolated the components. Mine are the more accurate.
- **K sensitivity [E]** (`code/verdict.py`). R_vec carries K² at fixed E_X.

  | K | S/Q(10⁻⁶) | e_crit (bulk) |
  |---|---:|---:|
  | 1.75 | 4.34 | 1.60×10⁻⁶ |
  | 2.97 | 8.03 | 4.57×10⁻⁶ |

  The E_X of those inflaton histories is not computed.

**Settled with quick checks.**
1. **e = 2×10⁻⁶.** dq_entropy reading B said Q-matter (3.370); the direct heat says free X (3.516). The direct heat
   supersedes the transfer. The branch stays free X even at the outer band's low edge, 3.484, but only by 0.4%.
2. **The H-growth cancellation.** dq_heat estimated a factor 1.06–1.15; its skeptic found f = 0.62–0.46 [C]. I re-solved
   dq_heat's own held minimum (`code/exit_check.py`, read-only import of dq_model.py).
   - f = 0.623 at s = h_r = 1.5 and 0.507 at 2, against the estimate's 1.055 and 0.911.
   - H alone gives f = 1.048 (h_r = 1.5) and 1.043 (h_r = 2): F/H_exit stays at 270 because the ξ_X R hold scales |X|
     with H.
   - The skeptic is right.
3. **Readings A and B of dq_entropy and qm_micro's e(S/Q)** are superseded by the direct heat.
   - e(S/Q) = 5.41×10⁻⁶ (2.5), 2.98×10⁻⁶ (3.0), 2.51×10⁻⁶ (3.2), 2.06×10⁻⁶ (3.4716).
   - qm_micro's reading-B values were 4.9 / 2.7 / 2.27×10⁻⁶.
4. **M_P condition.**
   - The direct M_P budget gives S/Q(10⁻⁶) = 4.09 after the self-consistent-L shift of −0.01 [E]. This supersedes
     dq_entropy's transferred 4.01 / 4.34.
   - Free-X H_i = 1.0724×10¹¹, at L = 13.537 (dq_entropy SKEPTIC N2).

## 2. Item 2: the Q-matter branch at the benchmark e = 3×10⁻⁶ [C/E]

**Relic** (`code/qm_benchmark.py`; qm_micro with its skeptic's corrections, interpolated to S/Q = 2.993) [C on E
models].

| Quantity | Perfect interior trapping | Beliaev-only trapping |
|---|---|---|
| ⟨Q⟩ (charge-weighted; Bose + latent-heat-limited exchange) | 1.16–1.26×10²⁰ | ~0.2–1.1×10²¹ (not recomputed with corrections) |
| Vapour fraction (energy-conserving), charge / mass | 0.119–0.121 / 0.149–0.151 | 0.175 / 0.214 |
| M/m_X | 0.802 | 0.815 |
| T_f (m_X) | 0.0089 | 0.0125 |
| Kinetic freeze-out T_vis | 0.54–0.99 MeV | 10.9–12.9 MeV |
| ζ_eff | 1.0–1.1×10⁻⁷ | 3.6–5.6×10⁻⁹, for every Q ≪ 10²⁴ (qm_transfer SKEPTIC) |
| σ_p today (qm_transfer convention) | 1.7–2.0×10⁻¹¹ | 1.7–1.8×10⁻¹² |
| λ_fs | 0.18–0.33 kpc | 0.025 kpc |

In both cases s̄(T_f) = 0.47 and p ≈ 0.

**P(k)/P_ΛCDM at z = 0.** qm_transfer's direct runs at qm_micro's points sit on the code's floor (|1 − P| ≤ 1.2×10⁻³).
I therefore extracted the physical deficit from qm_transfer's own grid. The family is X_fc0.1_Tf0.01 (f_v = 0.150 by
mass, T_f = 0.010, matching the benchmark). The steps:
1. **Remove the floor.** Divide by the ζ = 10⁻⁷ member. The floor pattern repeats to ≤ 2×10⁻⁵ across all ζ = 10⁻⁷
   models. Then correct for that member's own deficit.
2. **Extrapolate in σ_p.**
   - Upper estimate: the local exponent measured between ζ = 10⁻⁶ and 10⁻⁵. It is 1.95 at k = 12 and 1.43 at k = 96, and
     it rises toward 2 as σ falls.
   - Lower estimate: n = 2 from ζ = 10⁻⁵.
3. **Correct and scale.**
   - Multiply the upper estimate by the converged K_rsa = 800/200 deficit ratio, 1.02–1.07 (qm_transfer SKEPTIC).
   - Scale linearly with f_v.

Labels: [C] on the grid, [E] for the extrapolation.

| k (1/Mpc) | 10 | 20 | 50 | 100 |
|---|---|---|---|---|
| **benchmark, perfect trapping: 1 − P/P_ΛCDM** | 3–8×10⁻⁸ | 1.1–3.0×10⁻⁷ | 0.5–5.8×10⁻⁶ | 1.2–42×10⁻⁶ |
| benchmark, Beliaev trapping | 0.4–1.3×10⁻⁹ | 1.5–4.5×10⁻⁹ | 0.6–17×10⁻⁸ | 0.16–2.0×10⁻⁶ |
| thermal relic 5.7 keV (Viel fit; current Lyα bound) | 0.051 | 0.218 | 0.827 | 0.998 |
| thermal relic 20 keV (a "20 keV-class" future reach) | 0.0023 | 0.011 | 0.082 | 0.327 |

**Excluded already?** No. The benchmark deficit is 4–6 decades below what the current Lyman-α (5.7 keV, arXiv:2309.04533),
Milky Way satellite (6.5 keV, arXiv:2008.00022) and lensing (6.1 keV, arXiv:2405.01620) bounds allow. qm_transfer's
scaling law, extrapolated beyond its 1.5–80 keV fit range [E], gives an effective thermal mass of:

| Probe | Perfect trapping | Beliaev trapping |
|---|---|---|
| Lyα | ≈ 340 keV | ≈ 940 keV |
| satellites | ≈ 250 keV | ≈ 630 keV |
| lensing | ≈ 190 keV | ≈ 390 keV |

**In reach of upcoming probes?** No.
- A probe reaching 10–20 keV sees this vapour only for ζ ≳ 1.2–2.3×10⁻⁵ (f_v = 10–20%). The benchmark lies 110–210× below
  in ζ.
- Reaching the band would need ⟨Q⟩ ≲ 10¹³–6×10¹³ with perfect trapping, 10⁷ times smaller than ripening gives.
- With Beliaev trapping, no object size reaches it.

**Is there an e range with a distinctive P(k) prediction?** No.
- Across the whole Q-matter range (e from 2.06×10⁻⁶ up to e → ∞, S/Q = 3.47 → 1.84), ⟨Q⟩ stays at 0.4–1.5×10²⁰. It is
  set by ripening and nearly independent of S/Q.
- f_v rises only to ~0.13–0.19 (charge) near the threshold, and ζ stays ~10⁻⁷ or ~5×10⁻⁹.
- WORKING_THEORY §6's "O(10%) deficit near k ≈ 10–50/Mpc for ⟨Q⟩ ~ 10¹²" does not occur in M1r.

**What is distinctive is the composition** [D on C]. At the benchmark:
- objects of M = (8.4–9.1)×10¹⁸ GeV = 1.5–1.6×10⁻⁵ g = **0.68–0.74 M_Pl**;
- radius 0.85 nm;
- σ/m = 1.5×10⁻⁹ cm²/g;
- a flux of 0.24–0.35 m⁻² yr⁻¹ at Earth (ρ = 0.3–0.4 GeV cm⁻³, v = 230 km/s).

Over the branch the masses run from ~0.2 M_Pl (S/Q → 1.84) to ~7 M_Pl (the Beliaev upper end).
- **What could see them.** With zero portals only gravity couples. This is the target of proposed gravitational direct
  detection with large arrays of quantum-limited impulse sensors (Carney, Ghosh, Krnjaic, Taylor, arXiv:1903.00492;
  Windchime roadmap, arXiv:2203.07242, with a flux scale of ~1 m⁻² yr⁻¹ × M_Pl/m). No scheduled experiment reaches it [R].
- **Why it matters.** A gravitational detection at that mass and flux would separate the Q-matter branch from free X,
  which is a 92 MeV particle that cannot be seen gravitationally. That is a branch discriminator in principle, not an
  upcoming test.

## 3. Item 3: H_i and r per branch [C on Q/E inputs] (`code/hi_calibration.py`)

**Formula.** H_i = H_M1 (0.9374/(M·Y/Y_Q))^{2/3} with the frozen Y/Y_Q = 0.9910212.
- Free X: M = 1 + 2p. p(S/Q) comes from wt_critic's map [E], whose ×2.5 disagreement with REVIEW_KINETIC case L gives the
  range.
- Q-matter: M = 0.802 (perfect trapping) to 0.815 (Beliaev).
- r = 2H²/(π²A_s M_P²), with A_s = 2.099×10⁻⁹ [R Planck 2018].

| branch (e) | M/m_X | H_i (GeV) | with fragmentation-dephasing systematic | M_P portals | r (M_Y / M_P) |
|---|---|---|---|---|---|
| free X (10⁻⁶, static) | 1.0019–1.0047 | **1.2899–1.2923×10¹¹** | 1.284–1.312×10¹¹ | 1.0724×10¹¹ | 2.71×10⁻⁷ / 1.87×10⁻⁷ |
| free X (1.5–2×10⁻⁶) | 1.0001–1.0008 | 1.2933–1.2939×10¹¹ | – | Q-matter: 1.229–1.244×10¹¹ | 2.73×10⁻⁷ / 2.47×10⁻⁷ |
| Q-matter (10⁻⁶ on rolling paths) | 0.804–0.818 | 1.479–1.497×10¹¹ | – | – | 3.6×10⁻⁷ |
| **Q-matter (3×10⁻⁶ benchmark)** | 0.802–0.815 | **1.483–1.499×10¹¹** | 1.477–1.521×10¹¹ | 1.236–1.248×10¹¹ | 3.58–3.66×10⁻⁷ / 2.49×10⁻⁷ |
| Q-matter (10⁻⁵) | 0.795–0.804 | 1.497–1.508×10¹¹ | – | 1.244–1.253×10¹¹ | 3.65×10⁻⁷ |
| free X (10⁻⁷) [E] | 2.68–2.80 (p ≈ 0.84–0.9) | 6.5–6.7×10¹⁰ | – | 7.6×10¹⁰ | 7×10⁻⁸ |

How to read the table:
- The dephasing systematic is N_X −2.2% … +0.64% (dq_entropy SKEPTIC N3).
- The r values assume H_CMB = H_i. A rolling history has H_CMB ≥ H_i, so they are floors [A].
- The supplied-history sensitivities of WORKING_THEORY §5 still apply on top [Q]: H_i ∝ φ_I(slice)^{−2.23}, and the
  inflaton history moves H_i by −16% or +21%.

## 4. Item 4: Higgs survival with ξ_H = 0.1 [C runs, E criterion] (`code/higgs_check.py`)

**Method.** The higgs_reheat SKEPTIC's normalisation-fixed copy of the mode code (`modes_canon.py`, read-only, unchanged)
is run at this pass's H_i. Reheating is gravitational only, N = 60, with the scale-resolved criterion. ξ_EW = 0.1 runs to
ξ(H_i) = 0.0878.

| M_t (GeV) | H_i | D₆₀ (c_H = ½ / 2) | verdict |
|---|---|---|---|
| 172.60 | 1.2923×10¹¹ (free X) | −181 (−209 / −222) | safe |
| 172.87 | 1.2923×10¹¹ | **+12.7** (+6.6 / +6.0) | fails |
| 172.60 | 1.4908×10¹¹ (Q-matter) | −136 | safe |
| 172.87 | 1.4908×10¹¹ | **+25.8** (+20.5 / +19.7) | fails |
| 172.87 | 1.0724×10¹¹ (M_P free X) | −7.9 | safe |

**Critical top mass** (bisection to 4 MeV):

| Branch | M_t,crit (GeV) | Distance above PDG 2024 (172.57 ± 0.29) |
|---|---|---|
| free X | **172.84** | +0.91σ |
| Q-matter | **172.79** | +0.76σ |
| M_P free X | **172.89** | +1.10σ |

- Each value carries ±0.09 GeV from the criterion's conventions (±30 in D, higgs_reheat SKEPTIC §3, divided by the
  computed slope dD/dM_t ≈ 335 per GeV near the crossing).
- If M_t were Gaussian, M_t < M_t,crit would hold with probability 0.78–0.86 [E].

**Answer.** Is ξ_H = 0.1 inside the window? Yes during inflation: ξ_EW,min = 0.044–0.055, a margin of ×1.8–2.3. Through
reheating, gravitational only, yes for M_t below these critical values; the lower edge at central M_t is 0.072–0.085. No
above them: the edge is 0.096–0.112 at +1σ.

Three conditions on that answer:
1. **Criterion.** It rests on our argued scale-resolved criterion. The literature's end-of-inflation patch count
   (Kearney–Yoo–Zurek arXiv:1503.05193; Espinosa et al. arXiv:1505.04825) would call even central M_t "fails" without a
   plasma.
2. **Plasma rescue.** It needs Higgs-coupled states with T ≳ 10¹¹ GeV by τ ≈ 8–10 H_i⁻¹, which is 3.5–5×10⁻³⁵ s. M1r's
   reflection-even inflaton forbids φ_I|H|², and its singlet-fermion route decays too late (~5×10⁴ H_i⁻¹) [A, higgs_reheat
   SKEPTIC].
3. **Thermal decay** is negligible: S3/T ≥ 478.

## 5. Final verdict table

| e | S/Q (static slice) | branch | relic | H_i (GeV) | r | observable | status |
|---|---|---|---|---|---|---|---|
| 10⁻⁷ | 22.7 | free X | cold asymmetric X gas with a large symmetric part (p ≈ 0.8) | 6.5–6.7×10¹⁰ | 7×10⁻⁸ | none | heat [C]; p, H_i [E] |
| **10⁻⁶ (declared)** | **4.955** (4.915–4.994) | **free X** on the static slice; **Q-matter** (3.11–3.37) on every rolling embedding | free X: p ≈ 1–2.4×10⁻³, M = 1.002–1.005 m_X | 1.290–1.292×10¹¹ (1.072 with M_P portals; 1.48–1.50 if Q-matter) | 2.7×10⁻⁷ | none (CDM-like; σ/m 6×10⁻⁸ cm²/g) | **conditional on the unsupplied pre-slice inflaton history** [C on E paths] |
| 1.5×10⁻⁶ | 4.015 | free X (Q-matter with M_P portals) | free X | 1.293×10¹¹ | 2.7×10⁻⁷ | none | [C] |
| 2×10⁻⁶ | 3.516 (3.484–3.548) | free X, margin 0.4–2.2% | free X | 1.294×10¹¹ | 2.7×10⁻⁷ | none | [C], marginal |
| **2.06×10⁻⁶** | 3.4716 | bulk e_crit (1.46×10⁻⁶ with M_P portals; 0.75–0.93×10⁻⁶ rolling; 1.6–4.6×10⁻⁶ over K) | – | – | – | – | [C] |
| **3×10⁻⁶ (benchmark)** | **2.993** | **bulk Q-matter** | ⟨Q⟩ = 1.2×10²⁰ (0.7 M_Pl objects), vapour 12–17% of the charge, M = 0.802–0.815 m_X, p ≈ 0 | **1.483–1.499×10¹¹** | 3.6×10⁻⁷ | P/P_ΛCDM ≥ 1 − 4×10⁻⁵ at k ≤ 100/Mpc; Planck-mass objects (gravitational detection only) | [C] heat; [C/E] relic; [E] P(k) |
| 3.19×10⁻⁶ | 2.928 | Q = 2×10⁵ e_crit | – | – | – | – | [C] |
| 5×10⁻⁶ | 2.552 | Q-matter (binding objects too, above 4.99×10⁻⁶) | ⟨Q⟩ ≈ 8×10¹⁹ | 1.491–1.504×10¹¹ | 3.6×10⁻⁷ | as at 3×10⁻⁶ | [C/E] |
| 10⁻⁵ … ∞ | 2.20 → 1.84 | Q-matter | ⟨Q⟩ 4–7×10¹⁹ (0.2–0.4 M_Pl) | 1.50–1.51×10¹¹ | 3.7×10⁻⁷ | as above | [C/E] |

**Higgs, every branch.** ξ_EW = 0.1 is safe during inflation. It is safe through gravitational reheating for
M_t ≤ 172.79–172.89 GeV, and needs an unsupplied prompt Higgs-coupled plasma above that.

## 6. Contradictions settled

| Contradiction | Settled as | How |
|---|---|---|
| Branch at e = 2×10⁻⁶: Q-matter (dq_entropy reading B) vs free X (direct heat) | **free X, margin 0.4–2.2%** | direct heat supersedes the transfer; outer band low edge 3.484 > 3.4716 [C] |
| H growing with φ_I "nearly cancels" (dq_heat) vs f = 0.46–0.77 (skeptic) | **skeptic** | own re-solve: f = 0.623 / 0.507 at s = h_r = 1.5 / 2; f = 1.04–1.05 for H alone [C] |
| e_crit on rolling paths 0.81–0.95×10⁻⁶ (skeptic) | **0.75–0.93×10⁻⁶** | component interpolation instead of ln W interpolation [C] |
| e(S/Q) used by qm_micro | 5.41 / 2.98 / 2.51×10⁻⁶ for S/Q = 2.5 / 3.0 / 3.2 | direct heat [C] |
| Free-streaming length: qm_micro 0.34–0.61 kpc vs qm_transfer "~10 pc" | **0.18–0.33 kpc** (perfect trapping); 25 pc (Beliaev) | own integral ∫√3σ_p/(a³H) da [C]; qm_transfer's "~10 pc" is ~20× low; neither matters |
| Vapour fraction: qm_micro 0.08–0.11, qm_transfer 0.14–0.15, qm_micro SKEPTIC 0.12 (S/Q = 3) | **0.12 (perfect) – 0.175 (Beliaev)**, charge | energy-conserving, area-corrected [Q] |
| Higgs plasma deadline: τ ≈ 15 (track) vs 8–10 (skeptic) | **8–10** | consistent normalisation [Q] |
| Higgs at +1σ: "marginal, D = +1" (track) vs +13 (skeptic) | **fails under the criterion**, +12.7 at H_i = 1.2923×10¹¹ | own run [C] |
| M_P free-X H_i 1.076 vs 1.0724×10¹¹ | **1.0724×10¹¹** | self-consistent L = 13.537 [C] |
| Free-X H_i 1.2930 vs 1.2923×10¹¹ | **1.2899–1.2923×10¹¹** | p at the direct S/Q, with the p-map spread [C on E] |
| Onset a ≈ 3.8×10¹² (brief) vs 5.25×10¹² | **5.25×10¹² in the M1r frame** | dq_heat [Q]; 3.8×10¹² was the spot frame |

## 7. Remaining open inputs, ranked by what they move

1. **The inflaton history before the slice** (holder value and H when k_* left the horizon, 3–4 e-folds before the slice).
   - It decides the branch at the declared e = 10⁻⁶: free X needs f_eff ≥ 0.469; the rolling embeddings give 0.35–0.44.
   - It also sets H at CMB exit, hence r, and the slice value φ_I, hence H_i ∝ φ_I^{−2.23}.
   - M1r's (m_I = 3H_i, φ_I = 0.8165 M_P) is a quadratic potential released at rest, which cannot itself be the
     slow-roll history [A].
2. **The inflaton decay history and channel** (m_I, Γ_I).
   - K = 1.75–2.97 changes R_vec by 0.78–2.24 and e_crit to 1.6–4.6×10⁻⁶ [E]; E_X for those histories is uncomputed.
   - The decay channel decides whether a Higgs-coupled plasma forms by τ ≈ 8–10 H_i⁻¹.
3. **The induced-portal renormalisation condition.** e_crit 2.06 → 1.46×10⁻⁶, H_i −17%.
4. **The top mass, together with the Higgs decay criterion.** M_t,crit lies 0.8–1.1σ above the PDG central value, with
   ±0.09 GeV from the criterion. The standard patch count disagrees.
5. **Fragmentation-time dephasing of X ↔ Σ.** H_i +1.5% … −0.4%. Needs a lattice with Y and Σ.
6. **Completeness of chemical relaxation before onset** [Q]. It is the one omitted process that lowers S/Q, by up to ~0.1
   if only the massless channels act.
7. **Q-matter interior trapping (perfect vs Beliaev).** f_v 0.12 vs 0.175, M 0.802 vs 0.815, H_i ±0.5%, ζ 10⁻⁷ vs
   5×10⁻⁹. No observable consequence either way.
8. **p(S/Q).** The ×2.5 spread between methods moves free-X H_i by 0.2%, and by more at e ≲ 3×10⁻⁷.
9. Small: the single-Higgs-direction treatment of sub-horizon vectors (≤ 0.003 computed, ±0.05 outer in R); Hartree vs
   Bogoliubov spectrum (~10⁻³ of S); g* at onset.

## 8. Corrections to WORKING_THEORY.md (for the synthesis)

- **§0 / §4.** "S/Q ≈ 4.7–5.1 at e = 10⁻⁶" becomes **4.955 on the static slice**.
- **"Robust in direction" is false.** The pre-slice inflaton history gives 3.11–3.37 (Q-matter), and the decay history
  gives 4.34–8.03.
- **"Q-matter needs e ≳ 2×10⁻⁶"** becomes e > 2.06×10⁻⁶ (bulk), 3.19×10⁻⁶ (Q = 2×10⁵) or 4.99×10⁻⁶ (binding); 1.46×10⁻⁶
  with M_P portals.
- **§8's failure condition is met** on the rolling embeddings (W_vec 0.35–0.42 of the broad transfer) and not on the
  static slice (0.94).
- **Q-matter relic.**
  - ⟨Q⟩ = 0.4–1.5×10²⁰: Planck-mass objects, not "10¹²–10²⁰".
  - The vapour is 10–19% of the charge at S/Q = 2.5–3.2 (10–13% perfect trapping, 14–19% Beliaev), not 3–11%.
  - M = 0.80–0.82 m_X, not 0.78–0.80.
  - **The §6 small-scale signature does not occur** (≤ 4×10⁻⁵ at k = 100/Mpc).
- **H_i and r.** H_i = 1.290–1.292×10¹¹ (free X) and 1.48–1.50×10¹¹ (Q-matter); r = 2.7 and 3.6×10⁻⁷, as floors.
- **Higgs.** "Not computed through reheating" becomes: computed, safe gravitationally for M_t ≲ 172.8–172.9 GeV, and
  dependent on an unsupplied plasma above that.

## Code and reproduction (one process, ~7 min; `code/`)

| File | Item | Reuse |
|---|---|---|
| `sq_branch.py` | 1 | reads dq_heat, dq_entropy and skeptic results; `hartree_copy.py` is an **unchanged copy** of dq_entropy/code/hartree_dq.py |
| `qm_benchmark.py` | 2 | reads qm_transfer grid, qm_micro_eval and the skeptic K_rsa runs; Viel fit [R] |
| `hi_calibration.py` | 3 | p map from wt_critic; M from qm_micro SKEPTIC |
| `higgs_check.py` | 4 | imports higgs_reheat/skeptic/code/modes_canon.py read-only |
| `exit_check.py` | §6 check | imports dq_heat/code/dq_model.py read-only |
| `verdict.py` | 5 | assembles `results/verdict.json`, `results/summary.json`; K sensitivity |

```
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/dq_critic/code; export OMP_NUM_THREADS=1
python3 -B sq_branch.py && python3 -B qm_benchmark.py && python3 -B hi_calibration.py && python3 -B exit_check.py
python3 -B higgs_check.py points            # 1 min
for H in 1.2923e11 1.4908e11; do python3 -B higgs_check.py mtcrit $H; done; python3 -B higgs_check.py mtcrit 1.0724e11 172.87 173.10
python3 -B verdict.py
```

The Higgs runs were made from the scratch directory `$SP/wt2/dq_critic`. Results are in `results/*.json` (under 110 kB in
total).
