# SKEPTIC of track qm_transfer (Q-matter + kinetic vapour transfer function)

2 November 2026 · Claude, adversarial pass (key qm_transfer_skeptic; BRIEF_DECISIVE.md of this directory; issue #13
context). This is same-session AI work on one supplied model (M1r, Q-matter branch), not outside replication.
- Every input is a model choice, not a consequence of closure. Data comparisons are retrospective, and the probe
  criteria are proxies.
- Git was read-only and nothing was committed. At most 2 worker processes, all in the foreground.
- The track's code was imported read-only and never modified.
- Labels: [D] derived, [C] computed with our own code, [E] estimate, [A] argued, [Q] quoted from a reviewed track,
  [R] literature (arXiv abstracts read on 2 Nov 2026).
- Code: `skeptic/code/`. Results: `skeptic/results/` (192 kB in total).

## Verdict in one paragraph

Five things hold:
- the code itself;
- the qualitative vapour physics (re-emission after capture freeze-out, the lag, kinetic freeze-out 1.7× early);
- the kinetic-vs-fluid gap (R3);
- the physical-line limit, ζ_max ≈ 3×10⁻³ at f_v ≈ 3%, within about ±20%;
- the qm_micro verdict (CDM-like).

Four things do not hold as stated:
1. **Production numerics.** The production switch K_rsa = 200 is not converged at k ≳ 30/Mpc. At K_rsa = 800 (which
   equals 1600 to ≤ 1% at k ≤ 96):
   - P(96) is 0.282 instead of 0.335 at ζ = 3×10⁻³, and 0.070 instead of 0.117 at 10⁻²;
   - the thermal relic's "4–6% colder than CLASS" disappears: k_hm agrees with Vogel & Abazajian's CLASS fit to ≤ 0.7%.
2. **The heating coefficient.** "σ_p,thermo = 0.038 ζT_f, 3× sudden" is an artefact of the 24-bin q grid, which ends
   at 25 q_f,eq. The rms sits in a < 1% tail. With the grid extended, the coefficient is 0.079–0.083 (5.6–5.9×
   sudden). The robust mean ⟨q⟩ is 2.3× sudden.
3. **Where the thermostat works.** An object can act as the T_o thermostat only if it can hold one latent heat:
   S(Δ) > Δ/T_o. That fails for Q ≲ 3×10⁸–10⁹ (ζ ≳ 10⁻³ at s = 1), where re-emission conserves energy (N_e = N_a).
4. **Object size of the band (new, decisive).** qm_micro's own interior-trapping estimate (Beliaev decay) gives
   ζ_eff = 0.216 s̄/ℓ ≈ 4–6×10⁻⁹ for every object with R ≪ ℓ ≈ 2×10⁷ m_X⁻¹, independent of Q. The "observable band
   at ⟨Q⟩ ≈ 4×10⁷–3×10¹², s ≈ 1" then has no realisation in M1r as modelled. Reaching it needs a thermalisation length
   5×10⁴–2×10⁶ times shorter than Beliaev's.

## 1. Checks of the track's results

| ID | Claim | Recomputation | Verdict |
|---|---|---|---|
| R1 | Validation (EH98 ±3%; k_hm within 4–6% of the CLASS-based fit; ≤ 1% vs two-fluid; controls; Murgia δA) | VA22 fit (2210.10753 Table 2: a = 0.0437, b = −1.188, ν = 1.049, θ = 2.012, η = 0.2463), own evaluation: k_hm = 28.2 / 46.2 / 50.4 / 58.9 at 3.5 / 5.3 / 5.7 / 6.5 keV; Viel 25.9 / 41.1 / 44.6 / 51.6. The track's 28.6 / 46.3 / 50.9 / 59.5 and 26.3 / 41.7 / 44.8 / 52.3 agree to ≤ 1.5%. The own code at **K_rsa = 800** (`sk_wdm_K.py`): k_hm = **50.56 (5.7 keV), 59.28 (6.5 keV)** against VA22 50.35 / 58.85, and P(34, 48, 68) within 1% of VA22. At the production 200 it is 53.9 / 62.7 | **holds**; the 4–6% offset is the K_rsa = 200 bias, not "colder like CLASS" [C] |
| R2 | Energy balance forces re-emission; σ_p,thermo = 0.038(ζT_f)^0.99 vs 0.0124 sudden; lag 0.57; a_f 1.7× earlier | Independent background (`sk_vapour_bg.py`: own cosmology, 160–320 bins up to 1800–8000 q_f, exponential integrator). Reproduced: sudden σ_p 6.85×10⁻⁹ (track 6.93); a_f,kin 1.497×10⁻⁷ (1.502); a_f,eq/a_f,kin = 1.72; T_v/T_o = 0.556 at Γ/H = 2. **Thermo:** with the track's grid edge (25 q_f,eq) σ_p = 2.45×10⁻⁸ (track 2.53). It grows with the edge: 3.05 / 3.51 / 3.85 / 4.06×10⁻⁸ at 100 / 400 / 1800 / 8000 q_f, so C = 0.079–0.083 (×5.6–5.9 sudden). ⟨q⟩, q50 and q90 change by ≤ 2%: ⟨q⟩ ×2.3, q50 ×1.29, q90 ×2.5 over sudden. Mechanism (object MC, below): N_e/N_a = 0.963–0.969, equal to the canonical value, but only for Q ≳ 6×10⁸ | **refuted as a number** (grid-truncated rms; the mechanism holds only for Q ≳ 10⁹) [C] |
| R3 | Kinetic suppression ≪ two-fluid where tight coupling fails | Rerun at K_rsa = 800 (`sk_twofluid_K.py`): Q = 10⁸: P(50) = 0.715 (vs 0.521), P(100) = 0.559 (vs 0.210). Q = 10¹⁰: P(100) = 0.787 (vs 0.720) | **holds**; the gap shrinks by 2–6% [C] |
| R4 | P(k) on S/Q = 1.235 at k = 5…96 | K_rsa = 800 (K1600 equal to ≤ 0.3% at k ≤ 68, 1% at 96). ζ = 10⁻³: P(48, 96) = **0.746, 0.567** (track 0.762, 0.602). 3×10⁻³: **0.528, 0.282** (0.559, 0.335). 10⁻²: **0.294, 0.070** (0.336, 0.117). k ≤ 24 within 1.6% | **refuted at k ≥ 48 for ζ ≥ 10⁻³** (5–40% too little suppression); k ≤ 24 holds [C] |
| R5 | m_eff,Lyα = 19.3 keV (σ_p/10⁻⁸)^−0.49 (f_v/0.1)^−0.45 | The fit is to the track's truncated σ_p. Plugging in the converged σ_p (3.85×10⁻⁸) at ζ = 10⁻⁴ gives 15.8 keV against the computed 18.2. It mixes thermo (tail-inflated) and sudden models; S at ζ = 10⁻³ deviates by 0.15 dex | holds only with the track's σ_p convention; not portable [C] |
| R6 | ζ_max: physical S/Q = 1.235: 3.0×10⁻³ (Lyα), 3.8×10⁻³ (sat) | `sk_limits.py`: production 3.03 / 3.61×10⁻³. K800 thermo: **2.81 / 3.23×10⁻³**. K800 sudden (the right physics at Q = 3.7×10⁷, §2b): **3.74 / 3.48×10⁻³**. Lyα with k_up = ∞ for model and reference: 4.5–5.5×10⁻³ | **holds to ≈ ±20%**; the k_up choice is worth ±60% [C on E] |
| R7 | Observable band; ⟨Q⟩ ≈ 4×10⁷–3×10¹² at s ≈ 1 | The ζ band shifts by ≤ 5% (K_rsa). The ⟨Q⟩ translation needs perfect thermalising capture. Under qm_micro's Beliaev trapping, ζ_eff ≤ 3.6×10⁻⁹ (T = 0.004) … 5.6×10⁻⁹ (0.009) for every Q (§2c). With perfect trapping and s̄(0.004) = 0.35, Q must be 23× smaller than at s = 1 | ζ band holds; **⟨Q⟩ band not realisable** [C on E] |
| R8 | CDM-like at qm_micro's point; f_v 1.3–1.7× qm_micro's (b) | ζ_eff ≲ 10⁻⁷, so σ_p ≲ 10⁻¹⁰ by both of my backgrounds' scaling; no K_rsa issue at the 10⁻³ floor. Object MC at Q ≫ Q_c reproduces the canonical thermostat (N_e/N_a = 0.96–0.97). qm_micro's SKEPTIC §6 independently finds energy-conserving f_v 42–49% above variant (b) [Q]. X–X scattering is 4–44× the object rate there (§2d): no observable change | **holds** [C] |
| R9 | Consistent with Boyarsky (any mass for f < 0.35) and Baur (f ≲ 0.14 m^1.37); quote ζ_max ≈ 2–3×10⁻³ | 0812.0010 abstract: "within the mass range 5 keV < m_NRP < ∞" (thermal ≳ 1.1 keV), F < 0.6 (99.7%), F < 0.35 (Bayesian 95%). The vapours at the edge have m_eq ≈ 0.2–0.6 keV, outside that range. 1706.03118 abstract: F ≈ 0.35 (m/keV)^1.37 (BOSS) and 0.20 (with high-z data), 3σ, not 0.14. With 0.20 and the sudden σ_p appropriate at Q = 3.7×10⁷ (m_eq = 0.60 keV): F ≤ 0.099 ≫ f_v = 0.029, so no tightening | **refuted** (misquotes; the tightening to 2–3×10⁻³ is unsupported) [R/E] |

## 2. New issues

### a. Radiation-streaming switch [C] (`sk_grid_K.py`, `sk_krsa_metrics.py` → `krsa_metrics.json`, `krsa_raw_runs.json`)

Nine production models were rerun on the full 12-point k grid at K_rsa = 800 (one also at 1600), each with its own
CDM reference.
- The ratio converges at 800. 1600 changes it by ≤ 0.3% at k ≤ 68.
- The production ratio is biased high, i.e. too little suppression. The track's caveat quoted ≤ 9% at k = 68. In fact:
  - k = 68: 10%;
  - k = 96: 16% (ζ = 3×10⁻³) and 40% (ζ = 10⁻²).
- The metrics shift little: m_eff,Lyα 5.72 → 5.60 keV and m_eff,sat 6.83 → 6.64 keV at ζ = 3×10⁻³. So ζ_max moves by
  −7 to −11%.
- The same bias makes the own-code thermal relics colder than CLASS (R1). At 800 they reproduce VA22.
- Consequence: `judge_own` ("+30% tighter") is effectively the CLASS-calibrated reference, not an own-code oddity.

### b. Thermostat validity [C on E] (`sk_object_mc.py` → `object_mc.json`; `sk_trap_xx.py`)

**Model.** A microcanonical Gillespie simulation of one object in the track's own post-freeze-out vapour
(`sk_probe_bg.py`):
- Weisskopf evaporation from ρ(U) = e^{S(U)};
- internal modes: bulk phonons (c_s = 0.5467 [Q]) plus capillary ripplons (σ = 98.8, e + P = 348.6 [Q]);
- perfect internal thermalisation, which is the track's most favourable case.

**Canonical limit.** In this limit the kernel is exactly the track's thermostat. For Q = 10⁹–10¹² it reproduces T_o
(to ≤ 2%), N_e/N_a (to ≤ 0.1%) and ⟨E_out⟩ = 2T_o (to ≤ 3%). This is an independent confirmation of the mechanism.

**Failure for small objects.** An object can sit at the threshold holding one latent heat only if S(Δ) > Δ/T_o.
- Ripplons dominate S below Q ≈ 5×10¹⁵.
- Q_c = 2.5×10⁸ (T_o = 0.005), 5.6×10⁸ (0.004) and 1.5×10⁹ (0.003).
- In the scan at fixed vapour (Δ/T_o = 57.4), the heating efficiency is:

| Q | 10⁷–10⁸ | 3×10⁸ | ≥ 6×10⁸ |
|---|---|---|---|
| heating efficiency | 0.000 (N_e/N_a = 1.000) | 0.30 | 0.98–1.01 |

- On the physical line, ζ = 3×10⁻³ (Q = 3.7×10⁷) and ζ = 10⁻² (Q = 10⁶) are fully non-thermostat. ζ = 10⁻³
  (Q = 10⁹) is canonical.
- There the 'sudden' family is the right model: no heating, no net condensation, quasi-elastic re-emission.
- Effect: ζ_max = 3.5–3.7×10⁻³ (§1 R6). The track's f_v = 1% entry (ζ ≈ 6×10⁻³, Q ≈ 5×10⁶) is in this regime too.
  The S/P difference at 3×10⁻³–10⁻² (m_eff 6.2 vs 5.6 keV) suggests +20–40% in ζ_max there [E].
- Lowest ripplon quanta (6.5×10⁻³ m_X at Q = 10⁸) make small objects even more elastic than the continuum ρ(U).

**Numerical note.** A first, energy-grid version of the simulation stalled at U = Δ and produced spurious heating.
The continuous-U version is the one reported.

### c. Interior trapping makes the band's object sizes unreachable [C on qm_micro's E] (`trap_xx.json`)

qm_micro's §1f [Q] has two ingredients:
- an absorbed X becomes one phonon (k ≈ 0.41);
- unless it decays within ~R, it is re-emitted at the same energy, i.e. an elastic event.

So thermalising capture is s_eff = s̄ (1 − e^{−8R/3ℓ}), with ℓ = 2.05×10⁷ (T = 0.004–0.005) and 1.82×10⁷ (0.009)
m_X⁻¹ [Q]. Then:
- ζ_eff = s̄ P_trap Q^{−1/3} → 8 s̄ (3/4π n_B)^{1/3}/(3ℓ) = **3.6×10⁻⁹ (T = 0.004), 4.0×10⁻⁹ (0.005), 5.6×10⁻⁹ (0.009)**
  for all Q ≲ 10²⁴;
- for example, Q = 10⁹ gives P_trap = 1.0×10⁻⁵ and ζ_eff = 3.6×10⁻⁹ (3.5×10⁻⁴ with perfect trapping).

Reaching ζ = 7×10⁻⁵ / 10⁻⁴ / 3×10⁻⁴ / 10⁻³ / 3×10⁻³ needs both of the following:
- Q ≤ 1.2×10¹¹ / 4.1×10¹⁰ / 1.5×10⁹ / 4.1×10⁷ / 1.5×10⁶, even with perfect trapping (s̄ = 0.35);
- a thermalisation length ≲ R, which is 5×10⁴ / 7×10⁴ / 2×10⁵ / 7×10⁵ / 2×10⁶ times shorter than Beliaev's.

Landau damping is about 40× weaker [Q]. Planar phonon-to-ripplon conversion is kinematically closed (k_∥ = 0.41 <
k_ripplon(w) ≈ 0.56) [D]. The qm_micro skeptic puts the Beliaev inputs at ×2–3 [Q]; that leaves ζ_eff ≲ 10⁻⁸.

In M1r as modelled, therefore, no ⟨Q⟩ reaches the observable band, and the Q-matter branch is CDM-like for every
object size. Ripening sizes (10²⁰) are not needed for this conclusion.

### d. X–X scattering [C] (`trap_xx.json`)

Inputs: σ_XX = 6.0×10⁻⁸ cm²/g × m_X = 2.17×10⁻⁷ m_X⁻² (wt_critic [Q]). The ratio of X–X to object collisions is
N_v σ_XX √2 / (K Q^{−1/3} C_o^{2/3}):

| ζ (P line) | 10⁻² | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ | 10⁻⁷ | qm_micro points |
|---|---|---|---|---|---|---|---|
| X–X / object | 3.3×10⁻⁵ | 3.7×10⁻⁴ | 4.2×10⁻³ | 0.05 | 0.6 | 7.7 | 3.8–44 |

- The track's "≈ 10⁻⁴" holds only near ζ = 10⁻². The constraint-relevant region is safe.
- At qm_micro's point the vapour is self-collisional until a ≈ (2–7) a_f (a fluid, not free-streaming). With
  σ_p ~ 10⁻¹¹ there is no observable change.

### e. Smaller items
- **Lyα criterion.** The k_up = 75 h/Mpc calibration (to Murgia's δA_REF) moves ζ_max by +60% if k_up = ∞ is used
  for both model and reference, because the mixed-model plateau above 50/Mpc enters only in the latter [C]. It is a
  defensible calibration, but it is the largest single systematic in R6.
- **Mapping to ⟨Q⟩ with qm_micro's sticking.** It uses s ≈ 0.47. At the band's T_f ≈ 0.004–0.005, s̄ = 0.35–0.38
  [Q], so with perfect trapping Q must be 19–23× (not 10×) below the s = 1 values.
- **Objects' own thermal motion** (omitted, checked): at Q = 10⁶, ζ = 10⁻², the all-DM-equivalent mass is m_eq ≈ 37 keV
  (≈ 270 keV at Q = 3.7×10⁷). This is negligible next to the vapour [E].
- **f_v.** My f_v is 2–5% below the track's at the same model (0.0267–0.0277 vs 0.0280, ζ = 10⁻⁴). qm_micro's skeptic
  separately finds the single-⟨Q⟩ area makes it ~5% high [Q]. Neither matters for any verdict.

## 3. What survives, for the synthesis

- Kinetic code: the code is correct when run converged (K_rsa ≥ 800). Production P(k) at k ≥ 48 is biased toward CDM.
- Re-emission after freeze-out: correct for objects with S(Δ) > Δ/T_o (Q ≳ 10⁹). Its σ_p coefficient is
  grid-dependent; ⟨q⟩ is 2.3× sudden.
- Limits: ζ_max ≈ 3–4×10⁻³ at f_v ≈ 3% (Lyα ≈ satellites), with ±60% from the Lyα k_up choice.
- **Physical realisation:** with qm_micro's trapping, ζ_eff ≈ 4–6×10⁻⁹ for any Q. In M1r as modelled there is no
  P(k) signature on any scale computed (k ≤ 136/Mpc), for any object size. The parametric band needs an interior
  thermalisation mechanism 10⁴–10⁶ times faster than Beliaev decay, which is outside the computed model.

## Reproduction (OMP_NUM_THREADS=1; scratch $SP/wt2/qm_transfer_skeptic; run from it; S = tracks/qm_transfer/skeptic)

```
python3 -I $S/code/sk_grid_K.py gridK_a.json 800 540 P_SQ1.235_z3e-03,P_SQ3.4_z3e-04   # ~4 min per model; also b–f:
#   b: P_SQ1.235_z1e-03,P_SQ3.4_z1e-03 (800)  c: P_SQ1.235_z3e-03 (1600)  d: S_SQ1.235_z3e-03,S_SQ1.235_z1e-03 (800)
#   e: P_SQ1.235_z1e-02 (800)  f: S_SQ1.235_z1e-02 (800)
python3 -I $S/code/sk_krsa_metrics.py $S/results/krsa_metrics.json gridK_[a-f].json
python3 -I $S/code/sk_limits.py $S/results/krsa_metrics.json $S/results/limits_recomputed.json
python3 -I $S/code/sk_wdm_K.py wdmK_57.json 5.7 ; python3 -I $S/code/sk_wdm_K.py wdmK_65.json 6.5      # ~2.5 min each
python3 -I $S/code/sk_twofluid_K.py tf_K_1e8.json 1e8 50,100 ; python3 -I $S/code/sk_twofluid_K.py tf_K_1e10.json 1e10 100
python3 -I $S/code/sk_vapour_bg.py $S/results/vapour_bg.json                                          # 20 s
python3 -I $S/code/sk_probe_bg.py $S/results/probe_bg.json                                            # 1 min
python3 -I $S/code/sk_object_mc.py $S/results/probe_bg.json $S/results/object_mc.json 20000           # 75 s
python3 -I $S/code/sk_trap_xx.py $S/results/trap_xx.json
```
`sk_krsa.py` is the first timing probe (superseded by `sk_grid_K.py`).

Literature read for this pass (abstracts):
- arXiv:0812.0010 (Boyarsky et al.);
- arXiv:1706.03118 (Baur et al.);
- arXiv:2210.10753 (Vogel & Abazajian, Table 2 via ar5iv);
- arXiv:2309.04533 (Iršič et al., 5.7 keV; 4.1 keV with k_max < 0.1 s/km).
