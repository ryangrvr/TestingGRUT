# e2e_higgs: Higgs survival on the final histories

3 November 2026 · Claude, end-to-end pass (BRIEF_E2E.md, track e2e_higgs; GitHub issue #13 context).

This is same-session AI work on one supplied model. It is not outside replication. Nothing here is a validation, an
exclusion of the real world, a new theory or a GRUT prediction.

**Every choice is a model choice, not a consequence of closure.** That includes:
- HFa;
- ξ_H = 0.1, read as ξ(M_t);
- the argued criteria (Hawking–Moss, the scale-resolved danger exponent D, the medium masses);
- the closures of the histories (local vs box-mean Pauli blocking, the splice at t_sw = 3);
- the H_i range.

PDG inputs (M_t = 172.57 ± 0.29, M_h = 125.20 ± 0.11, α_s = 0.1180 ± 0.0009) are recalled [R], and comparisons with them
are retrospective.

Labels: [D] derived, [C] computed here, [Q] quoted from a track, [E] estimate or argued criterion, [A] argued,
[R] recalled or literature.

## 0. What was done (lean, as the brief asks)

**Code reuse.** The mf_higgs code is imported read-only and unchanged:
- `mh_modes.run`: Higgs modes, Bunch–Davies start, RK4, Hartree;
- `mh_modes.danger`: the criterion D;
- `mh_inflation`: HM counts;
- `mh_common`: the N_obs bookkeeping;

on `inf_higgs/ih_common.SMTab` and `higgs_reheat/smpot`.

**What is new** (code/eh_common.py):
1. A background class `FinalBG` with the same interface as `mh_common.Background`, built from the final histories.
2. The medium-mass column, with the mf_higgs skeptic's N-gas term.
3. The exact ε_H = 1 endpoint in the instantaneous count (code/eh_infl.py).

**Regression [C].** On the track's own export background the same code gives D(172.57 / 172.86) = −298.9 / −7.78. The
track has −299 / −7.8.

**Histories** (dimensionless, H_i = 1 and a = 1 at σ = σ_c). The raw homogeneous bench R is used up to t = 0.19 in every
case, then:

| tag | t > 0.19 | K (own) |
|---|---|---|
| KIN | e2e_kinetic fluid_HFa (cycle-mean R), homogeneous parent, self-blocked shell | 16.35 |
| LAT_P0 | e2e_lattice P0 box means (a, H, R_cyc, R_rms) to t = 3, then e2e_kinetic's ε_H(t) integrated from the lattice end state [A] | 14.6 (splice) |
| LAT_GM | same with the lattice's box-mean saturation closure (P0GM) | 17.6 (splice) |
| LATF_P0 | lattice P0 then the lattice's own HB.stage2 fluid | 3.6 |

Construction choices:
- **ln a.** It is rebuilt from ∫H dt, so that the friction and a³ agree (Wronskian error ≤ 1.2×10⁻⁴). The tabulated a of
  the kinetic export differs by ≤ 0.019 in ln a, a row-sampling effect that moves D by 0.06.
- **H_i.** The range is 5×10¹⁰–2.4×10¹¹ GeV. e2e_chain had no results when this ran, so no H_i was settled. At each H_i the
  dimensionless history is the design one [A]; N_obs uses e2e_kinetic's K(H_i) = 9.19 / 16.35 / 17.49 [Q].

**Media [E]** (the mf_higgs skeptic's model):
- **N gas.** m² = Σ_j y_j² ρ_Nj / E_N², using e2e_sm's S1 Yukawas y² = (0.054, 0.156, 0.043, 0.026) (∝ H_i). The
  per-species shares come from e2e_kinetic, or from the lattice's f_rad − f_SM.
- **E_N.** It is held at the injection bound ¾ μ5,max = 8.2×10³ H_i (conservative).
- **SM.** On KIN: e2e_kinetic's own column k_H T_cal². On the lattice histories: the hard gas 0.0228 ρ_SM/(E_N/2)².

## 1. During inflation (code/eh_infl.py → results/inflation.json) [C on E]

N0 is computed with the kinetic K. The waterfall to ε_H = 1 is common to all histories, since the lattice is homogeneous
there to about 10⁻⁵.

| H_i (GeV) | ξ_min const-H, tree (172.57 / 172.86) | loop | cap at ξ = 0.1, const-H | instantaneous with exact ε_H = 1 end: cap (α_s 0.1171 / 0.1189); ξ_min |
|---|---|---|---|---|
| 5×10¹⁰ | 0.0456 / 0.0505 | 0.036 / 0.043 | > 176 | 173.93 (173.60 / 174.27); 0.064 / 0.074 |
| 1.3×10¹¹ | 0.0487 / 0.0531 | 0.042 / 0.047 | > 176 | 173.73 (173.39 / 174.07); 0.071 / 0.079 |
| 2.4×10¹¹ | 0.0503 / 0.0546 | 0.044 / 0.050 | > 176 | 173.61 (173.27 / 173.95); 0.074 / 0.082 |

- **Regression.** At 1.3×10¹¹ these reproduce the mf_higgs skeptic's numbers: B_end = 425.13, log₁₀⟨N⟩ = −113.5 / −62.2
  / −12.9 / +10.6 against its −112.1 / … The offset of −1.4 is exactly the K = 16.35 shift in N0.
- **Reading.** As the skeptic showed, the instantaneous count is not a bracket, because it diverges 0.003 H⁻¹ past ε_H = 1.
  It is listed only as the corrected form of the old "pessimistic bracket".
- **Verdict.** The meaningful inflationary verdict is the constant-H one: ξ_H = 0.1 is safe, with a cap above 176 GeV, at
  every H_i.

## 2. After inflation, gravitational only (code/eh_post.py main → results/post.json) [C on E]

M_t,crit is where D crosses 0, on a 0.05 GeV grid. The band covers the criterion variants c_H = 0.5 / 2 and k_c × 0.5 / 2.

| history @ H_i | N_obs | D(172.57 / 172.86) | **M_t,crit** (band) | α_s 0.1171 / 0.1189 | M_h 125.09 / 125.31 |
|---|---|---|---|---|---|
| KIN @ 1.3e11 | 57.01 | −342 / −21.5 | **172.914** (172.884–173.005) | 172.586 / 173.240 | 172.851 / 172.975 |
| LAT_P0 @ 1.3e11 | 57.07 | −336 / −19.4 | **172.909** (172.877–173.000) | 172.583 / 173.237 | 172.847 / 172.971 |
| LAT_GM @ 1.3e11 | 56.97 | −349 / −22.5 | **172.916** (172.886–173.006) | 172.589 / 173.242 | 172.854 / 172.978 |
| LATF_P0 @ 1.3e11 | 57.77 | −319 / −12.5 | **172.893** (172.847–172.990) | 172.566 / 173.219 | 172.832 / 172.951 |
| KIN @ 5e10 | 56.82 | −1325 / −283 | **173.221** (173.186–173.329) | 172.900 / 173.547 | 173.161 / 173.284 |
| KIN @ 2.4e11 | 57.28 | −105 / **+33.8** | **172.736** (172.708–172.820) | **172.408** / 173.065 | 172.677 / 172.798 |

- **The lattice and the plasma closure barely matter: ±0.02 GeV at fixed H_i.**
  - The lattice histories at 5×10¹⁰ give 173.198–173.225.
  - At 2.4×10¹¹ they give 172.717–172.738.
  - The N_obs choice (own K, kinetic K or lattice-fluid K, 3.6–17.6) moves M_t,crit by ≤ 0.003.
  - Running to t_end = 400 instead of 80 changes it by ≤ 0.0003 (results/control_tend.json).
  - The roll at the critical mass is at τ_c ≈ 18–21 H_i⁻¹.
- **H_i is the lever.**
  - At 1.3×10¹¹ the critical mass is +1.2σ above PDG.
  - At 5×10¹⁰ it is +2.2σ.
  - At 2.4×10¹¹ it is +0.6σ, and D(172.86) = +34 fails at +1σ (the mf_higgs skeptic had +34.5 on mf_chain's Pauli
    history).
- **α_s ∓ 1σ moves M_t,crit by −0.33 / +0.33 GeV.**
  - At α_s = 0.1171, central M_t is marginal at 1.3×10¹¹: D = −6.0 (KIN), −4.1 (LAT_P0), −6.7 (GM), +2.5 (LATF_P0).
  - At 2.4×10¹¹ it fails: D = +40.3.
- **M_h ∓ 1σ moves it by ∓0.06 GeV.**
- **Deadlines** (gravitational only, D(τ_th) = 0) at 1.3×10¹¹ give the time by which a Higgs-coupled medium must exist:
  - 173.0 GeV: 8.3–11.3 H_i⁻¹;
  - 173.15: 5.6–7.3;
  - 173.45: 3.6–4.4;
  - 174.6: 1.47–1.79.
  - At 2.4×10¹¹: 4.3–5.2 at 173.15 and 1.33–1.61 at 174.6.

## 3. The lattice's ⟨R⟩ and its fluctuations [C, with A]

- **⟨R⟩/H² after ε_H = 1.**
  - The local-blocking lattice (P0) is more radiation-like: R_cyc/H² ≈ 1.1 at t = 1–3, against 2.0–2.3 for KIN and GM.
  - It still moves M_t,crit by only −0.005 (splice) to −0.02 GeV (own fluid). D is set mostly at τ ≈ 5–20, after the
    lattice block.
- **R fluctuations, pessimistic local variant.** Here R̄ is replaced by R̄ − n·R_rms(t), with R_rms the instantaneous
  spatial rms of R(x) over the box: 2.6 / 3.5 (P0 / GM) at t = 0.5, and 0.10 / 0.35 at t = 3.
  - n = 1: no change (P0), or −0.05 GeV (GM) at 1.3×10¹¹.
  - n = 2: M_t,crit = 172.565–172.574 (−0.34 GeV), and 172.41–172.42 at 2.4×10¹¹.
- **Why the variant is a bound only [A].**
  - The fluctuations sit at k ~ a m_χ (m_χ ≈ 480 H_i) and oscillate at about ω ~ m_χ.
  - A rolling region has size 1/k_c, with k_c ≈ 0.3–0.5 at a ≈ 10–20, so it averages over about (a m_χ/k_c)³ ~ 10¹³–10¹⁴
    fragments and over many periods.
  - Its R is then the mean to about R_rms × 10⁻⁷.
  - So only the box mean enters, and the lattice's R fluctuations do not change the result.
- **Upper edge [E].** A local tachyonic patch would need ξ|δR| ≳ m_χ², i.e. ξ ≳ 480²/4.6 ≈ 5×10⁴. That is far above the
  upper edge ξ_EW ≈ 111 (21) from the coherent oscillation at 172.57 (172.86) [Q mf_higgs skeptic; not recomputed here,
  because the final histories carry cycle-mean R after t = 0.19].

## 4. With the media [C on E]

On every history and every H_i, D < 0 for all M_t ≤ 174.6 GeV with N gas + SM. Smallest factor on the medium m² that
still gives D(174.6) < 0:

| | KIN (1.3e11 / 5e10 / 2.4e11) | LAT_P0 | LAT_GM | LATF_P0 |
|---|---|---|---|---|
| N + SM | 1e-5 / 3e-6 / 3e-5 | 1e-3 | 3e-3 | 3e-4 |
| N gas alone (e2e_sm y) | 0.03 / 0.1 / 0.03 | 0.01–0.03 | 0.03–0.1 | 0.01–0.03 |
| N gas, smallest y only | 0.1 / 0.1 / 0.03 | 0.03 | 0.03–0.1 | 0.03–0.1 |

**When the media are in place.**
- The N gas reaches 2.4×10³ H_i² (KIN) or 1.1×10⁵ H_i² (lattice) by t = 0.2.
- Its peak is 2–4×10⁵ at t = 0.25, against −V″ ≤ 0.5 H_i² inside the abundance limit.
- The SM column exceeds 10³ H_i² from t ≈ 0.25 (KIN). So the media precede the earliest deadline (1.3–1.8 H_i⁻¹) by more
  than 1 H_i⁻¹.

**Margins.**
- The N gas alone has only ×10–100 margin. Its m² falls as the N decay (f_N ≈ 0.005 at t = 100), while −V″ stays fixed.
- Later safety rests on the SM medium. On KIN that is e2e_kinetic's k_H T_cal². On the lattice histories it is the hard-gas
  formula on the lattice's f_SM, which is 10–20× smaller, hence the 10⁻³ margin.
- The SM hard-gas sign carries the skeptic's gauge caveat (N3). After t_bath_thermal = 6 H_i⁻¹, the thermal mass
  (κ_T > 0) applies [Q e2e_kinetic].

**Lower ξ edge** (results/xiscan.json):
- Gravitational only, after inflation, at 1.3×10¹¹ on all four histories: ξ_EW = 0.066–0.067 (172.57), 0.093–0.096
  (172.86) and 0.111–0.115 (173.15).
- At 5×10¹⁰: < 0.05 / 0.068 / 0.093. At 2.4×10¹¹: 0.083 / 0.106 / 0.121.
- So at 2.4×10¹¹, ξ_H = 0.1 is inside the gravitational-only window only below about 172.8 GeV.
- With the media, the ξ_EW = 0.05 runs give D(174.6) = −∞ (KIN) and −6554 (LAT_P0). The post-inflation edge then drops
  below 0.05, so the window's lower edge is the inflationary one.

## 5. Verdict for ξ_H = 0.1 on the final histories (model choices as above)

**ξ window.** The lower edge is 0.046–0.050 (central M_t) and 0.051–0.055 (+1σ), set during inflation over H_i = 5×10¹⁰–
2.4×10¹¹, tree. The loop schemes are lower. The upper edge is ≈ 111 (21) [Q]. ξ_H = 0.1 sits about ×2 above the lower edge.

- **During inflation.** Safe at every H_i in range: ξ_min = 0.046–0.055, cap > 176 GeV.
- **After inflation, gravitational only.** M_t,crit = 172.91 ± 0.02 (closures) at 1.3×10¹¹; 173.20–173.23 at 5×10¹⁰;
  172.72–172.74 at 2.4×10¹¹.
  - The criterion band is −0.03 / +0.09 GeV.
  - α_s ∓ 1σ moves it by ∓0.33 GeV; M_h ∓ 1σ by ∓0.06 GeV.
  - The lattice, its two blocking closures and its R(x) fluctuations change this by ≤ 0.02 GeV (the n = 2 local bound
    aside).
- **With the media.** Safe to ≥ 174.6 GeV on every history and H_i, with ×10³–10⁵ margin (×10–100 from the N gas alone).
  This does not depend on the gravitational-only marginality.

## 6. Not established / caveats

- e2e_chain's self-consistent H_i was not available, so the 5×10¹⁰–2.4×10¹¹ range is reported instead. At each H_i the
  dimensionless design history is reused [A]; only N_obs, ξ(H_i) and the SM inputs change.
- e2e_chain appeared while this ran, with scan files only (results/scan_sHFa_*.json, "H_i_req" 2.1–6.5×10¹¹ GeV,
  not settled). For use there: M_t,crit (gravitational only, KIN) falls by about 0.30 GeV per e-fold of H_i
  (173.221 / 172.914 / 172.736 at 5e10 / 1.3e11 / 2.4e11). Extrapolated [E], that gives ≈ 172.61 at 3.7×10¹¹ and
  ≈ 172.45 at 6.5×10¹¹, i.e. about central M_t. The media verdict does not depend on this.
- The history skeptic was not available (null in the inputs). LAT_P0/GM continue with e2e_kinetic's ε_H(t) after t = 3
  [A]: there is a jump in R/H² (1.16 → 1.84) for P0, and the own-fluid LATF_P0 brackets it.
- The final histories carry cycle-mean R after t = 0.19. The R lines and the R < 0 windows of the coherent stage are
  dropped (the track showed they move D by ≤ 0.04). The upper ξ edge was not recomputed.
- The media are hard-gas forward-scattering estimates with fixed injection energies, with no back-reaction. The kinetic
  history used m̃ = 0.05 eV for N2, while e2e_sm uses 0.0912 eV.
- D is the argued criterion of higgs_reheat / inf_higgs. Thermal tunnelling is not computed (it is minor per
  1608.02555 [Q mf_higgs skeptic]).

## 7. Files and reproduction

- `code/eh_common.py`: FinalBG, media, scratch-safe imports.
- `code/eh_post.py`: `main` → `results/post.json`; `xiscan` → `results/xiscan.json`.
- `code/eh_infl.py` → `results/inflation.json`.
- `code/eh_summary.py` → `results/summary.json`.
- `results/control_tend.json`: t_end = 400 control (inline driver using `run_cached(..., t_end=400)` and
  `eh_post.crit_scan`).

Mode runs and SM tables are cached in
`/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_higgs/` (about 1.5 s per
mode run; the whole track takes about 15 min on one core).

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/e2e_higgs/code; export OMP_NUM_THREADS=1
python3 -B eh_infl.py; python3 -B eh_post.py main; python3 -B eh_post.py xiscan; python3 -B eh_summary.py
```

## Literature (verified by web search this session unless marked)

- Mantziris, Markkanen, Rajantie, arXiv:2011.03763: vacuum-decay bounds on ξ during inflation (ξ_EW ≳ 0.06 in field
  models).
- Herranen, Markkanen, Nurmi, Rajantie, arXiv:1506.04065: curvature-driven Higgs production at the end of inflation for
  ξ ≳ 1.
- Figueroa, Rajantie, Torrentí, arXiv:1709.00398: lattice upper bound ξ ≳ 4–5 at M_t ≈ 173.3.
- Not re-verified here [Q mf_higgs skeptic]: Markkanen–Nurmi–Rajantie–Stopyra, arXiv:1804.02020 (curved-space one-loop
  potential); Salvio–Strumia–Tetradis–Urbano, arXiv:1608.02555 (thermal decay).
