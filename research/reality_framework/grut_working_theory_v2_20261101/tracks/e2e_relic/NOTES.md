# Track e2e_relic: the dark-matter relic machinery for M2

Claude, 3 Nov 2026. Brief: `BRIEF_E2E.md`, section e2e_relic. Labels: [D] derived, [C] computed, [Q] quoted from a track,
[E] estimate, [A] assumption, [R] literature. Every recipe, start convention, fit and threshold below is a **model choice**,
not a consequence of closure. Nothing here validates or excludes anything, and nothing here is a GRUT prediction. Data
comparisons are retrospective.

## Result in one paragraph

1. **The pair-fraction map gap (6.74 < S/Q < 17.1) is filled.** I ran the reviewed case-L network (fx_chain's `m1_kin`,
   read-only) at 12 points from 5 to 21.8 [C].
   - The ×1.85 "junction" the mf_chain skeptic found at S/Q = 6.74 is not a physical feature. The inherited maps
     **change start convention** there: the M1r gas map below uses the relaxed-onset start, and inf_chain and the
     skeptic above use forced equilibrium at a_ref.
   - Run with the relaxed start, the network reproduces the gas map at 6.742 (p = 0.016619 against 0.016617). It
     joins the a_ref rows from above, the gap between them closing to 3% at 21.8.
   - The corrected map uses one convention, relaxed, throughout. It is exported as `code/pmap.py`, with the table in
     `results/pmap.json`.
   - In the gap it sits above the inherited bridge by ×1.06–1.84 at S/Q = 7–16, and 6% below it at 17.1. Relative
     to the a_ref network it is 7–29% lower at S/Q = 8–14.
2. **The relic is now a function** of (S/Q, H_i, r_hv = s_h/s_vis, e): `code/relic.py`.
   - It reproduces qm_micro's M1r Q-matter benchmark, as corrected by its skeptic: at S/Q = 3.0, ⟨Q⟩ = 1.267×10²⁰,
     f_v = 0.1197 / 0.1751 and M = 0.8027 / 0.8152 m_X.
   - It reproduces the free-X point M = 1.002 (S/Q = 4.955) and the c0 = 1.40 value M = 0.807 [C].
   - New here: ⟨Q⟩ ∝ Y_Q^1.00 [C]. f_v depends only weakly on Y_Q: −7 to −8% per decade with perfect trapping, −18 to
     −20% per decade with Beliaev trapping.
   - At fixed S/Q, e does not enter, because the hidden U(1) is Higgsed with m_A = 240 e H_i ≫ m_X [D].
3. **Observables** [C on the model; local halo values R]:
   - At the Q-matter points the objects are 8–19 μg (0.37–0.86 M_Pl), 0.69–0.91 nm in radius.
   - Their local density is 2.4–5.7×10⁻¹⁴ m⁻³, and they cross a 1 m² area 0.27–0.67 times a year.
   - Their self-interaction is ≤ 7×10⁻⁹ cm²/g. The vapour gives m_eff,Lyα ≈ 160–250 keV, so P(k) is CDM-like.
   - Free X: σ/m = 6.0×10⁻⁸ cm²/g (XX) and 7.1×10⁻⁷ cm²/g (XX̄), with no annihilation channel to light states.

## 1. The corrected pair-fraction map (`code/pmap_runs.py`, `code/pmap.py`, `results/pmap_runs.json`, `results/pmap.json`)

**Method** [C on a conditional recipe].
- Engine: fx_chain/code/kin/m1_kin.py, the M1-patched copy of Codex 34f23f8's moment_controls. It is imported
  read-only, with its work directory redirected to my scratch.
- Gas-only rule, M1 inputs. The case-L network is the same recipe as inf_chain's `pmap_kin` and the mf_chain
  skeptic's `sk_pmap`, with R = W(S/Q)/13.0412 − 1.
- Two start conventions:
  - *relaxed*: x0 = ln(a_on/a_ref), linearly interpolated over fx_chain's heat_onset rows (R = 0…49.6, φ in
    equilibrium). The gap R = 8–30 lies between the computed rows R = 5.70 and 41.7, so this is interpolation.
  - *a_ref*: x0 = 0.
- Charge errors are ≤ 1.1×10⁻⁶ and energy–work errors ≤ 2×10⁻⁸ on every run.

| S/Q | 6.742 | 7.0 | 7.5 | 8 | 9 | 10 | 12 | 14 | 16 | 17.131 | 21.797 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| p, relaxed (this track) | **0.01662** | 0.02201 | 0.03486 | 0.05044 | 0.08759 | 0.1299 | 0.2229 | 0.3215 | 0.4228 | 0.4808 | 0.7231 |
| p, a_ref (sk_pmap / inf_chain) | 0.0308 | – | – | 0.0711 | – | 0.1539 | 0.2472 | 0.3453 | 0.4458 | 0.5034 | 0.7439 |
| relaxed / a_ref | 0.54 | – | – | 0.71 | – | 0.84 | 0.90 | 0.93 | 0.95 | 0.955 | 0.972 |
| inherited bridge (chain_mfi) | 0.0166 | 0.0191 | 0.0246 | 0.0312 | 0.0481 | 0.0708 | 0.138 | 0.244 | 0.399 | 0.513 | 0.762 |

- **The relaxed run at 6.742 equals the M1r gas map** (0.016619 against 0.016617) [C]. So the gas map *is* the
  relaxed network, and the ×1.85 junction is the start-convention switch, not a discontinuity in p.
- The same switch is ×3.5 at S/Q = 4.94 (fx_chain's own a_ref row, 0.00322 against 0.00091) [Q fx_chain].
- **Above 21.8** no relaxed rows exist. I use inf_chain's a_ref rows times 1 − 0.0280 (21.797/S)², a fit to the
  computed relaxed/a_ref ratio at 14–21.8 [E]. The correction is ≤ 0.8% for S/Q ≥ 41.6.
- **Model choice.** The central map is the relaxed start, the convention M1r's own map uses below the gap. The
  WORKING_THEORY §7 basis is "the plateau relaxes chemically before onset". `p_aref()` gives the a_ref edge as a
  start-placement band, which is 7–29% at S/Q = 8–14 and ≤ 5% above 16.
- **Yield dependence (new)** [C].
  - Writing f = Y_Q/5.0794×10⁻⁹, the network gives β = d ln p/d ln Y_Q = −0.383 (S/Q = 5), −0.192 (7), −0.083 (10)
    and −0.035 (21.8), from runs at f = 0.1 and 10.
  - inf_chain's self-consistent runs give −0.021 to −0.019 at 155–958.
  - `p_selfconsistent()` solves the abundance-calibrated f = 0.9374/(1 + 2p). That replaces inf_chain's sc_ratio,
    and agrees with it to ≤ 0.7% where both exist.
  - β is held constant outside 5–958 [E].
- **Slopes.** d ln p/d ln(S/Q) is 7.1 / 5.2 / 3.35 / 2.37 / 2.05 at S/Q = 7 / 8 / 10 / 12.2 / 14.7. The a_ref values
  are 4.9 / 4.2 / 3.0 / 2.2 / 1.96, and the bridge's is 3.68. For the isocurvature derivative in the gap, use the map
  (d ln M/d ln(S/Q) = 0.77 at 12.2 and 0.87 at 14.7) rather than any fixed slope.

**Consequences at fixed history** (M = 1 + 2p self-consistent; H_req ∝ M^−2/3) [C]:

| S/Q | 7 | 8 | 10 | 12.2 | 14.7 | 15.9 |
|---|---|---|---|---|---|---|
| M, inherited bridge | 1.038 | 1.062 | 1.142 | 1.294 | 1.584 | 1.780 |
| **M, corrected** | 1.045 | 1.103 | 1.267 | 1.479 | 1.737 | 1.867 |
| H_req, corrected / inherited | 0.996 | 0.975 | 0.933 | 0.915 | 0.940 | 0.969 |

- The exported/unblocked hybrid relics (S/Q 12.2–15.9) get M = 1.48–1.87 m_X, not 1.29–1.78. Their H_req falls by
  3–9%, inside the 2–13% that mf_critic flagged.
- Against the skeptic's a_ref fill (M = 1.77 at S/Q = 14.7), the corrected value is 2% lower.
- Pauli-point relics (S/Q ≈ 61.5) are unchanged: M = 6.85 against 6.8.
- The M2-X point (S/Q ≈ 4.16) has M = 1.0004, unchanged.

## 2. The relic as a function (`code/relic.py`, `code/relic_runs.py`, `results/relic_table.json`, `results/relic_runs.json`)

`relic(SQ, H_i=None, r_hv=None, e=1e-6, Y_Q=None, trapping='both')`:
- Y_Q = n_Q/s = r_hv/(S/Q).
- If neither r_hv nor Y_Q is given, the abundance-calibrated Y_Q = D/(M m_X) is solved with M (D = 4.4×10⁻¹⁰ GeV,
  supplied).
- The function returns the branch, p or ⟨Q⟩, f_v (perfect and Beliaev trapping), M (central and range), objects per
  charge, Ω/Ω_target, m_A, and flags for every extrapolation.

**Branch.** Q-matter iff S/Q < 3.4716, the bulk threshold [Q dq_entropy / fx_chain skeptic]. Ripened objects
(Q ~ 10²⁰) are bulk-like, so the Q = 2×10⁵ threshold of 2.928 no longer applies [E].

**Free X.** p from §1; M = 1 + 2p; X̄/X = p/(1 + p).

**Q-matter, inherited tables** [C on E models, Q qm_micro + skeptic]:
- ⟨Q⟩ from qm_micro's ripening with the skeptic's Bose and latent-heat factors (coarsen_heat bose_heat_sφ1 / sφ0).
- f_v from the skeptic's energy-conserving freeze-out (sk_freeze).
- M = 0.775866 + 0.224135 f_v; p ≈ 0.

| S/Q (at Y_Q,ref = 6.065×10⁻⁹) | 1.5 | 2.5 | 3.0 | 3.2 |
|---|---|---|---|---|
| ⟨Q⟩ (sφ = 1 / 0) | 3.82 / 3.61×10¹⁹ | 7.96 / 7.67×10¹⁹ | 1.267 / 1.17×10²⁰ | 1.45 / 1.50×10²⁰ |
| f_v perfect / Beliaev | 0.0555 / 0.0827 | 0.0971 / 0.1414 | 0.1197 / 0.1754 | 0.1292 / 0.189 |

**New: dependence on the hidden-to-visible ratio** (`relic_runs.py`; qm_micro's qm_vapour, the skeptic's coarsen_heat
and sk_freeze run unchanged with QM_MASS and SK_SCR redirected; S/Q = 3.0) [C on E models]:

| Y_Q / Y_Q,ref | ⟨Q⟩ | f_v perfect | f_v Beliaev | T_vis at Γ_a = H (perfect / Beliaev) |
|---|---|---|---|---|
| 0.1 | 1.240×10¹⁹ | 0.1304 | 0.2179 | 3.6 / 79 MeV |
| 1 (regression) | **1.267×10²⁰** (skeptic 1.27×10²⁰) | **0.1197** (0.1197) | **0.1751** (0.175) | 0.85 / 12.8 MeV |
| 10 | 1.272×10²¹ | 0.1111 | 0.1437 | 0.23 / 1.6 MeV |

- **⟨Q⟩ ∝ Y_Q^1.00** [C]. This is the ripening law [D]: at a fixed hidden state, T_vis ∝ Y_Q^−1/3 and H ∝ Y_Q^−2/3,
  so R² ∝ K_W t gives Q ∝ t^3/2 ∝ Y_Q.
- f_v ∝ Y_Q^−0.035 (perfect) and Y_Q^−0.09 (Beliaev).
- These exponents are applied at every S/Q [E: measured at one S/Q].
- At the calibrated abundance Y_Q is 5.9×10⁻⁹, so the exponents matter only if the chain's history misses the
  abundance.
- H_i enters only through Y_Q and the m_A check [D on the model].
- e enters only through S/Q(e), which the chain computes. The hidden U(1) is Higgsed: m_A = 6 e h = 240 e H_i =
  2×10⁷ GeV at e = 10⁻⁶. That gives no Coulomb bound on Q and no light vector. The code flags m_A < 10 m_X, which
  happens for e < 4×10⁻¹⁴ at H_i = 10¹¹ GeV.

**Regression** (`python3 -B relic.py`) [C]:
- S/Q 2.99 (M1r Q-matter benchmark, e = 3×10⁻⁶): ⟨Q⟩ = 1.26×10²⁰, f_v = 0.119 / 0.175, M = 0.8026–0.8150. The
  reference is ⟨Q⟩ 1.17–1.27×10²⁰, f_v 0.120 / 0.175, M 0.8023–0.8152.
- S/Q 2.88, calibrated: M_central = 0.8075, against WORKING_THEORY's 0.807.
- S/Q 4.955: M = 1.0019, against 1.002.

**Flags.**
- S/Q between 3.2 and 3.47 is a linear extrapolation of the tables [E]. The evaporation dip deepens toward the
  threshold (η_eq,min = 0.078 at 3.2), so ⟨Q⟩ there is least certain.
- S/Q < 1.5 is clamped. On M1r's heat, S/Q < 1.84 is unreachable.

## 3. Observables (`observables()` in relic.py, `results/relic_observables.json`)

Local halo: ρ = 0.3 GeV/cm³, v0 = 238 km/s, v_E = 250.5 km/s, v_esc = 544 km/s, following the conventions of Baxter et
al., arXiv:2105.00599 [R; the paper is verified, the values are as recommended there to my knowledge]. The Earth-frame
mean speed is 354 km/s [C]. The mass fraction in objects is 1 − f_v/M ≈ 0.85.

| benchmark (S/Q, calibrated Y_Q) | ⟨Q⟩ | f_v (perf/Bel) | M/m_X | object mass | R | n_local | crossings of 1 m² per yr | σ/m (contact) | m_eff,Lyα |
|---|---|---|---|---|---|---|---|---|---|
| M1r, e = 3×10⁻⁶ (2.99) | 1.22×10²⁰ | 0.119 / 0.175 | 0.803–0.815 | 15.6 μg = 0.72 M_Pl (3.6 reduced) | 0.86 nm | 2.9×10⁻¹⁴ m⁻³ | 0.33 | 5.9×10⁻⁹ cm²/g | 160–230 keV |
| HFa c0 = 1.40, 10⁻⁶ (2.88) | 1.12×10²⁰ | 0.114 / 0.168 | 0.802–0.814 | 14.3 μg = 0.66 M_Pl | 0.83 nm | 3.2×10⁻¹⁴ | 0.36 | 6.1×10⁻⁹ | 160–235 |
| HFa c0 = 1.40, 3×10⁻⁶ (2.12) | 6.3×10¹⁹ | 0.081 / 0.119 | 0.794–0.803 | 8.0 μg = 0.37 M_Pl | 0.69 nm | 5.7×10⁻¹⁴ | 0.67 | 7.4×10⁻⁹ | 170–250 |
| HFa c0 = 1.36, 3×10⁻⁶ (3.28) [E extrap.] | 1.47×10²⁰ | 0.133 / 0.195 | 0.806–0.820 | 18.8 μg = 0.86 M_Pl | 0.91 nm | 2.4×10⁻¹⁴ | 0.27 | 5.6×10⁻⁹ | 160–230 |

**Object size and flux.**
- Mass m_obj = 0.775866 ⟨Q⟩ m_X. Radius R = (3Q/4π n_B)^1/3 with n_B = 449.3 m_X³ [Q qm_micro].
- The ⟨Q⟩ uncertainty from φ absorption is ≤ 8%. The model uncertainty in ⟨Q⟩ is a factor ≲ 2 [Q qm_micro skeptic].
- The mean separation is about 32 km. The flux through area A is n⟨v⟩A: 0.3 per m² per year, 3 per (3 m)² per year
  [C]. This agrees with WORKING_THEORY's "0.3 per m² per year".

**Detectors** [R, verified by search; abstracts only].
- Carney, Ghosh, Krnjaic & Taylor, "Proposal for gravitational direct detection of dark matter", arXiv:1903.00492
  (PRD 102, 072003). They propose an array of quantum-limited impulse sensors and state that a meter-scale apparatus
  could detect any candidate "around the Planck mass or heavier".
- The Windchime Collaboration, Snowmass 2021 white paper, arXiv:2203.07242. It describes the programme of large
  mechanical-sensor arrays aimed at near-Planck-scale masses.
- The objects here are 0.4–0.9 M_Pl (2–4 reduced Planck masses), so they lie at the lower edge of that stated target.
- A 1 m² array would see 0.27–0.67 tracks per year, each crossing in about 3 μs/m. A 1 mg sensor at b = 1 mm receives
  an impulse I = 2 G m_obj m_s/(b v) ≈ 5×10⁻²⁷ kg m/s [C].
- I have **not** assessed whether any array reaches that impulse. No scheduled experiment is claimed.

**Self-interaction.**
- Q-matter: the contact (merger) cross-section π(2R)²/m_obj is 5.6–7.4×10⁻⁹ cm²/g. That is an upper estimate,
  7 decades below SIDM-relevant values (~0.1–1 cm²/g [E]).
- Free X: from tree-level λ_X|X|⁴ plus φ exchange (trilinear 2g²vφ|X|², m_φ² = 6 m_X²),
  σ/m(XX) = 6.0×10⁻⁸ cm²/g (reproducing WORKING_THEORY's 6×10⁻⁸) and σ/m(XX̄) = 7.1×10⁻⁷ cm²/g [D/C].
  The mixture is 6.0×10⁻⁸ (S/Q 4.16), 2.8×10⁻⁷ (14.7) and 3.8×10⁻⁷ (61.5).

**P(k) via qm_transfer.**
- ζ = s_eff ⟨Q_cap⟩^−1/3 with s_eff(T_f) = 0.47 gives ζ ≈ 1.0–1.05×10⁻⁷. Then σ_p = C(ζT_f)^0.99, with C = 0.038
  (track) or 0.081 (skeptic, converged), and f_v,mass ≈ 0.14–0.15.
- qm_transfer's fit then gives m_eff,Lyα ≈ 160–250 keV [Q fit, used at the lower edge of its ζ range; the skeptic says
  the fit is not portable to ≲ 15%].
- qm_transfer's direct runs at these ζ give |ΔP/P| ≤ 10⁻³ (the numerical floor) at k ≤ 136/Mpc [Q].
- So P(k) is CDM-like in the Q-matter branch, 3 decades in ζ below a 20 keV-class reach.

**Free-X relic for M2-X** [C via pmap; others D/Q]:
- Particle mass 92.4 MeV; local density 3.25 cm⁻³; flux 1.1×10⁸ cm⁻² s⁻¹.
- X̄/X = 2.1×10⁻⁴ at S/Q 4.16 (M = 1.0004) and 9.7×10⁻⁴ at 4.955 (M = 1.0019). In the gap it is 0.27 at 14.7
  (M = 1.74); at the Pauli point it is 0.75 (M = 6.85).
- X X̄ → φφ is closed (m_φ = 2.45 m_X) and there is no light hidden state, so there is no annihilation signal and
  ΔN_eff ≈ 0 [D, Q WORKING_THEORY §6].
- P(k) is CDM-like, with free streaming 0.3–0.5 pc [Q]. There is no direct detection at zero portals [Q].

## Code, reuse and reproduction

| file | what | reuse |
|---|---|---|
| `code/pmap_runs.py` | case-L network rows, relaxed or a_ref start, optional visible-clock Y_Q factor (inf_chain's patch) | imports fx_chain/code/kin/m1_kin.py read-only (work dir → scratch) |
| `code/pmap.py` | corrected map p(S/Q, Y_Q), p_aref, β, p_selfconsistent, M_freeX; writes results/pmap.json | data rows of wt_critic/fx_chain (gas map), inf_chain pmap_kin and mf_chain skeptic sk_pmap, embedded; chain_mfi imported read-only for comparison |
| `code/relic_runs.py` | Y_Q dependence of ⟨Q⟩, f_v | imports qm_micro/skeptic/code/sk_freeze.py (and through it coarsen_heat) read-only; histories from qm_micro/code/qm_vapour.py with QM_MASS = 0.785/f |
| `code/relic.py` | relic(), observables(), regression() | own; tables from qm_micro + skeptic results |
| `code/relic_tables.py` | writes relic_table.json, relic_observables.json, relic_runs.json | own |

Reproduction (OMP_NUM_THREADS=1, SCR = $SP/wt2/e2e_relic, about 8 min):
```
cd tracks/e2e_relic/code
python3 -B pmap_runs.py relaxed 5.0 6.742 7.0 7.5 8 9 10 12 14 16 17.131 21.797      # ~16 s per row
E2E_YQ_FACTOR=0.1 python3 -B pmap_runs.py relaxed 5.0 7.0 10 21.797; E2E_YQ_FACTOR=10 python3 -B pmap_runs.py relaxed 5.0 10
for f in 1 0.1 10; do mkdir -p $SCR/yq$f
  QM_MASS=$(python3 -c "print(0.785/$f)") python3 -B -W ignore ../../qm_micro/code/qm_vapour.py $SCR/yq$f/hist_3.0.json 3.0
  SK_SCR=$SCR/yq$f SK_CASES=bose_heat_sphi1 python3 -B -W ignore ../../qm_micro/skeptic/code/coarsen_heat.py $SCR/yq$f/cheat_3.0.json 3.0
  SK_SCR=$SCR/yq$f python3 -B -W ignore relic_runs.py $f $SCR/yq$f/relic_$f.json; done
python3 -B pmap.py; python3 -B relic_tables.py; python3 -B relic.py
```

## Limits

- **The network is a conditional recipe** [C, inherited]: a forced chemical-equilibrium start, no objects and MB
  fugacities. The start convention is a model choice. The a_ref band is 7–29% in p at S/Q = 8–14. Above 21.8 the
  relaxed/a_ref fit is [E].
- **Q-matter.** ⟨Q⟩ and f_v come from qm_micro's E-level models with the skeptic's corrections. They are not
  recomputed on a hybrid thermal history. The only history dependence the relic has at fixed S/Q is the visible clock
  (Y_Q), which I computed; the hybrid's post-release hidden history enters only through S/Q and Y_Q [D on the model].
  The Y_Q exponents were measured at S/Q = 3.0 only. Beliaev-trapping ⟨Q⟩ was not recomputed. S/Q 3.2–3.47 is
  extrapolated.
- **Detectors.** Only abstracts were verified. Rates are geometric (n⟨v⟩A). I have not checked any sensor-noise or
  sensitivity claim.
- **Not done (budget):** relaxed rows above S/Q = 21.8, Y_Q exponents at other S/Q, the hybrid-history ⟨Q⟩, and a 3D
  coarsening.
