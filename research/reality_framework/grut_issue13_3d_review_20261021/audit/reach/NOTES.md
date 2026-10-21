# Track `reach`: empirical reach of unchanged Candidate 1 in the classes that pass the wall tests

21 October 2026. Audit of Codex `a69bd86` (`grut_expanding_charge_20261021`) against Claude `c7b8480`
(`grut_issue13_wall_history_20261020`), issue #13. This note does not exclude, validate or identify Candidate 1. It asks
which observations could test the unchanged action, coupled to visible matter only through gravity
(hidden_heating), in the two classes that the wall tests do not reject:

- **cold histories:** no hidden restoration after inflation; charge supplied below n_c onto one vacuum sign;
- **biased early-decaying exits.**

Labels:

| Label | Meaning |
|---|---|
| [C] | computed by the code here |
| [D] | closed form |
| [E] | estimate |
| [Q] | quoted from the packets |
| [P] | preparation input |
| [RETRO] | retrospective datum, known before this analysis |
| [R] | recalled literature value, not fetched in this session; used only for orders of magnitude |

## Reproduce

```
cd research/reality_framework/grut_issue13_3d_review_20261021/audit/reach/code
OMP_NUM_THREADS=2 python3 qball_family.py     # ~1 min; reads the archived profile read-only -> results/qball_family.json
OMP_NUM_THREADS=2 python3 reach.py            # ~2 min -> results/reach.json
```

- `reach.py` reads `results/qball_family.json`.
- It quotes the Codex `479f111` production endpoint counts. To recompute the fine-grid 10% statistic, set
  `CODEX479_RAW` to a directory holding the extracted `raw/production_*.csv.components.csv`. These were extracted
  read-only with `git show 479f111:.../raw_diagnostics.tar.gz`.
- Verdicts are in `results/verdicts.json`.

## 1. What gravity alone lets the hidden sector change

With no portal, every visible effect goes through the stress tensor. That covers:

- the background expansion (total pressureless density, distances, H0);
- the growth of structure, including any departure from cold collisionless behaviour (self-interaction, free
  streaming, discreteness);
- extra radiation (N_eff, BBN);
- gravitational waves;
- compact objects and lensing;
- the initial perturbation type (isocurvature).

Direct detection, colliders, indirect detection and BBN injection are absent at this order. Graviton-induced
portals are ≲ 3e-34 [Q].

## 2. The "hidden charge is all of the dark matter" calibration [C, RETRO]

- Ω h² = (E/Q)·r·Y_B·s0/ρ_c,100 = 0.023799 GeV⁻¹ × (E/Q)·r.
- With ω_c = 0.120 [RETRO, Planck 2018] and Y_B = 8.674e-11 [Q, BBN], the calibration fixes one product:
  **e·g_X·v·r = 5.042 GeV**, where e = E/(m_X Q). It does not fix v, r, T_Q or the preparation separately.
- r at the benchmark v = 1 TeV:

  | Matter form | e | r |
  |---|---|---|
  | Q-matter floor | m*/m_X = 0.5623 | 0.0897 |
  | Codex 479f111 endpoint lumps (Q ≈ 5e5) | 0.71 | 0.071 |
  | archived soliton (Q = 39,944) | 0.911 | 0.0553 |
  | free X quanta | 1 | 0.0504 |

- **relic_yield's retrospective r_max = 0.0897 is CONFIRMED** (0.08966; ω_asym = 1.3383 r). It is the upper edge of
  the all-DM band, valid only in the Q-matter limit. Finite cold-branch Q-balls give r ≈ 0.055–0.071.
- Charge thresholds: T_n = 443–536 TeV. Fragmentation needs T_Q < T_n.
- Codex's maser product holds by construction: m_X r = 5.0–9.0 GeV < 63.93 GeV, and the floor
  H0 ≥ 100√(ω_b+ω_c) = 37.7 km/s/Mpc. The floor bounds ω_h from above, so with ω_c used to calibrate, distances and H0
  are the ΛCDM consistency test. They do not test Candidate 1.

## 3. Cold-branch Q-ball population, in units of v [C, E]

**Q-ball family.** `qball_family.py` reproduces the archived soliton exactly: Q = 39,944, E/(m_X Q) = 0.91108652,
R90 = 5.1034 ℓ. It then continues the branch:

| Q | E/(m_X Q) | R90 [ℓ] |
|---|---|---|
| 1.65e4 | 1.006 | 5.42 |
| **1.77e4** | **1.000** | 5.21 |
| 4.0e4 | 0.911 | 5.10 |
| 1.15e5 | 0.806 | 6.24 |
| 5.8e5 | 0.699 | 9.58 |

- The classical binding threshold is **Q_b = 1.77e4**. charge_phase's thin-wall estimate was ≈ 1.3e4 [E].
  **CORRECTED** (+36%).
- dQ/dΩ < 0 along the whole computed range.
- Large-Q extrapolation: e → m* + c₁Q^(−1/3) + c₂Q^(−2/3), with radii R ≈ (3Q/(4π n_*))^(1/3).
- Units: ℓ = 1/(100 GeV) (TeV/v) = 1.97e-16 cm (TeV/v). Mass M = e·g_X v·Q.

**Codex 479f111 production endpoint** [Q] (fine grid, unexpanding replay at 0.70 n_c, τ = 240; Codex's caveat:
threshold-dependent, not a converged mass function):

- At the 10% threshold, 39 components hold 68% of the charge. The 15 above Q_b carry most of it.
- Charge-weighted mean Q ≈ 6.2e5; the largest lump at 5% has Q ≈ 1.5e6. At 2% a single network holds 85%.

**Dependence on the supply density.**

- Linear fragmentation follows the charge_phase cubic. It is reproduced, then extended to small x:
  - γ_max = 0.450 x · 100 GeV (v/TeV);
  - k_max = 0.949 √x / ℓ.
- Per fastest-growing cell, Q_frag ≈ x(2π/k_max)³/λ_φ ≈ 2.9e4 x^(−1/2) [E]:
  - ≈ 1e5 near n_c;
  - growing for more dilute supply.
- The supplied condensate fragments only if γ_max ≳ 30 H [E]. At v = 1 TeV under the all-DM calibration, this sets a
  floor T_Q ≳ 0.35–0.6 GeV, where Q_frag ≈ 6e13–1e14.
- Below the floor the condensate stays homogeneous and is CDM-like. γ/H falls as 1/a, so it never fragments later by
  this instability.
- **Mass function [E, P].** The cold-class mass function therefore spans Q ≈ 1.8e4 → 1e14:
  - M ≈ 1.8e6 GeV → 5e15 GeV (v/TeV), i.e. 1.6e-51 → 5e-42 M_⊙ (v/TeV);
  - radii R ≈ 1e-15 → 1e-13 cm (TeV/v).
- Where in that range the population lies is set by T_Q and the supply form [P], not by the action.

## 4. Channel by channel, benchmark v = 1 TeV, all-DM calibration

| Channel | Candidate 1 value [C/E] | Reach / bound [R or RETRO] | Gap | Free inputs it depends on |
|---|---|---|---|---|
| ω_c (total pressureless density) | calibrates e g_X v r = 5.04 GeV | Planck 0.120 ± 0.001 | consumed by the calibration | r, v, e (preparation) |
| H0 / distances / maser floor | ΛCDM-identical once ω_c is fit; floor 37.7 km/s/Mpc | masers 73.9 ± 3.0 | not a C1 test | ρ_Λ (free) |
| σ/m, Q-balls | 3e-13 – 4e-12 cm²/g (contact π(R₁+R₂)², Q = 1.8e4 – 1.5e6) | SIDM bounds ~0.1–1, forecasts ~0.01–0.1 cm²/g | ≥ 10.4 decades | ∝ v⁻³ Q^(−1/3) |
| σ/m, elementary X / h | 7e-16 / 1.1e-14 cm²/g (threshold Born; XX amplitude 4\|g_X⁴/λ_φ − λ_X\| = 0.036) | same | ≥ 13 decades | ∝ v⁻³ |
| Free streaming / Lyman-α | λ_fs ≈ 2e-12 – 8e-12 Mpc for T_Q ≈ T_n; ≤ (1–6)e-6 Mpc at the fragmentation floor (β_i = 0.5–0.9) | λ_fs ≲ 0.03–0.1 Mpc | ≥ 3.7 decades | T_Q [P], β_i (not computed), v |
| Smallest halos / kinetic decoupling | never in kinetic equilibrium with the SM; free-streaming cutoff M_cut ≤ 4.5e-6 M_⊙ (≲ 1e-23 M_⊙ near T_n) | substructure probes ~1e6–1e8 M_⊙ | ≥ 11 decades in mass | as above |
| Lensing / compact objects | M ≤ 5e-42 M_⊙ (v/TeV) | femtolensing ≥ 1e-16, microlensing ≥ 1e-11 M_⊙ | ≥ 25 decades | v, Q |
| N_eff | ρ_h/ρ_rad = 1.3e-15 (T_Q = 440 TeV), 5.8e-13 (1 TeV). Even full conversion gives ΔN_eff ≤ 1.7e-12. The hidden spectrum has no massless state. | Planck ±0.17; CMB-S4 σ ≈ 0.03 | ≥ 10 decades | T_Q |
| BBN | hidden matter is pressureless, ρ/ρ_rad = 5.8e-7 at 1 MeV (as for any CDM); no injection | D/H, Y_p | none distinguishable | — |
| GW, cold-class fragmentation | emission fraction ≲ (ρ_h/ρ_r)² ≲ 3e-25 for T_Q ≥ 1 TeV (2e-30 at 440 TeV) | any detector | unobservable | T_Q |
| GW, early-decaying biased networks | see §5 | PTA (NANOGrav-15 ~1e-9–1e-8), SKA/LISA/ET ~1e-13–1e-12, BBO/DECIGO ~1e-17–1e-16 | ≥ 4 decades | m₀ [P], f, v |
| PBH / late minority domains | only cold histories with H_inf ≈ 210–560 GeV; M ≥ 4πM_P⁴/σ = 4.2e9 M_⊙ (v/TeV)⁻³; abundance changes ~1 decade per 1–3% in H_inf [Q minority_relic]; BH fraction not computed | SMBH/PBH limits (not imported) | possible but not calibrated | H_inf, N_*, start [P]; r_t ≈ 4e-25–5e-24, unmeasurable |
| Isocurvature | none from the action at a measurable level (spectator relic Y ~ 1e-25 [Q]); otherwise set by the [P] charge source | β_iso ≲ 0.035 | tests the source, not the action | source [P] |

**Biased class, extra components.**

- The symmetric relic from a bath at ξ_kin ≈ (4–7)e-5 gives Ω h² ≈ 2e-5 – 9e-5. It is a negligible fraction.
- Network debris decaying near T_safe with energetic quanta (p/m = 1/f) has
  λ_fs = 5e-5 / 3.6e-3 / 0.32 Mpc for f = 1 / 0.1 / 0.01.
- The f = 0.01 corner is Lyman-α relevant if that debris is all of the DM. That corner is tuned by m₀ [P], and f is
  not computed.

## 5. Gravitational waves from early-decaying biased networks [D, E]

**HKS radiation-era relation.** Ω_gw(t_ann) = (3ε̃/32π)(ρ_w/ρ_r)², with ε̃ = 0.7.

**A cap from the relic condition.** The relic condition that defines the surviving class,
f·ρ_w(T_ann)/s ≤ ρ_DM/s = 0.437 eV [RETRO], caps ρ_w/ρ_r **independently of v**:

Ω_gw h² (today) ≲ 2.5e-19 f⁻² (MeV/T_ann)²

| T_ann | f_peak today | Cap on Ω_gw h² |
|---|---|---|
| 0.3 MeV | 3e-11 Hz | 4e-18 f⁻² |
| 10 MeV | 1.1e-9 Hz (PTA) | 2.5e-21 f⁻² |
| 100 MeV | 1.2e-8 Hz | 2.3e-23 f⁻² |
| 10 TeV | 1.7e-3 Hz (LISA) | 1e-33 f⁻² |

**At the benchmark tension** (HKS, at the quoted T_safe band):

| f | T_safe | Ω_gw h² | Peak |
|---|---|---|---|
| 1 | 30.6 MeV | 3e-22 | 3.5e-9 Hz |
| 0.1 | 3.3 MeV | 2e-18 | 3.7e-10 Hz |
| 0.01 | 0.31 MeV | 6.5e-14 | 3e-11 Hz |

- Scaling: Ω ∝ f⁻⁴ (v/TeV)⁻⁶ and f_peak ∝ f (v/TeV)³.
- Networks that decay at formation are horizon-formula values of 2e-52 (500 TeV) and 3e-59 (2.6e7 GeV). The true
  sub-horizon values are lower still.
- Even the most favourable corner (f = 0.01, latest decay) lies below the PTA band. Inside the PTA band the cap is
  ≤ 2.5e-17. The LISA, ET and DECIGO bands are hopeless.
- **A detectable wall GW signal requires the network to decay close to its own domination.** The relic condition
  forbids that in Candidate 1, because every hidden species is massive and stable and couples only by gravity. Any
  such signal would also be generic to exact-Z₂ wall sectors.

## 6. Answer

**(i) Separate observational tests available without an action extension.**

| Observable | Discriminating power |
|---|---|
| ω_c | none: it is consumed by the calibration |
| H0 / distances / maser floor | none beyond ΛCDM (passes) |
| SIDM σ/m | none: ≥ 10 decades below reach |
| Lyman-α / small-scale / minimum halo mass | none in the bulk of both classes: ≥ 3.7 decades. Only a tuned debris corner (f ~ 0.01) is relevant, and there f is uncomputed and m₀ is [P] |
| Lensing / compact objects | none: ≥ 25 decades |
| N_eff / BBN | none: ≤ 1e-11 |
| GW | none: capped as in §5 |
| Isocurvature | tests the [P] source, not the action |
| PBHs from late minority domains | the one channel that could be large, but only for H_inf in a ~210–560 GeV window that cannot be calibrated independently (r_t ~ 1e-24). The BH formation fraction is not computed, and the signal is generic to Z₂ walls. At most it would calibrate (H_inf, v); it would not test them. |

**(ii) Stated plainly: no current or planned observation distinguishes the benchmark, calibrated so the hidden charge
is all of the dark matter, from collisionless cold dark matter.**

- The all-DM calibration consumes ω_c and fixes only e·g_X·v·r.
- No gravity-only observable within reach depends on the remaining parameters.
- This strengthens Codex's REVIEW verdict: the model cannot yet make a discriminating prediction. The block is not
  only missing calibrations. Under gravity-only coupling at v ~ TeV, the predictions themselves sit 4–25 decades below
  reach.

**What would be required** (not adopted; each needs an owner-approved, versioned change):

1. **An independent fix of the charge yield.** A versioned charge-source extension that ties Y_Q to something
   measured, e.g. Y_Q = Y_B, would let ω_c predict v. For r = 1 this gives e·g_X·v = 5.04 GeV, i.e. v ≈ 50–90 GeV.
   That v is still unobservable under gravity-only coupling, so (2) or (3) is needed as well.
2. **A scale change of the benchmark** (same action form):
   - SIDM: v ≈ 0.07–0.35 GeV gives Q-ball σ/m ≈ 0.1–1 cm²/g for Q = 1.8e4–6e5. All-DM then needs r ≈ 370 and
     T_Q < T_n ≈ 4 GeV.
   - This only becomes a prediction if an independently fixed preparation pins the mass function. σ/m ∝ Q^(−1/3) spans
     ~5 decades over the allowed Q.
   - The cold class at that scale needs H_inf ≲ 0.14–0.56 v [E, Claude's thresholds rescaled, not recomputed].
   - Lyman-α instead needs v ≲ 5–40 GeV together with the latest allowed supply T_Q (fragmentation floor), and the
     fragmentation velocities, which can be computed within the action.
3. **A nonzero visible-sector coupling.** Changing a matching choice (e.g. a Higgs portal κ|H|²φ²) makes v and the
   couplings measurable. But κ ≳ 1e-6 [E] thermalizes the hidden sector above T_c unless T_RH < T_c ≈ 1.56 v. That
   re-enters the rejected hot (H★-type) history, so any such extension must come with a low-reheating history and be
   rechecked against the wall tests.

The owner's statement, "connect a converged preparation to independently calibrated physical inputs and a separate
observational test", is necessary but not sufficient for these classes. At the benchmark scale a converged preparation
plus the all-DM calibration still yields no discriminating observable.

## 7. Verdicts (claims checked)

| Claim | Verdict | Value here |
|---|---|---|
| relic_yield r_max = 0.0897 (ω_h ≤ 0.12) | CONFIRMED | 0.08966; it is the Q-matter edge. Finite Q-balls give 0.055–0.071, free X 0.050 |
| ω_asym = 1.338 r (relic_yield) | CONFIRMED | 1.3383 r |
| "r ≈ 0.09 at the benchmark v" (task framing) | CORRECTED | 0.050–0.090. Codex-like cold Q-balls (Q ~ 4e4–6e5) need 0.055–0.071 |
| Thin-wall binding threshold Q ≈ 1.3e4 [E] (charge_phase) | CORRECTED | 1.77e4, from the computed family |
| Archived soliton Q = 39,944, E/Q = 0.91108652, R90 = 5.1034 | CONFIRMED | reproduced to 1e-9 |
| charge_phase growth rates 0.0045 / 0.045 / 0.112 / 0.226 and k = 0.404 | CONFIRMED | 0.0045 / 0.0448 / 0.1121 / 0.2264; k = 0.4035 |
| Codex maser product compatible with all-DM | CONFIRMED | m_X r = 5.0–9.0 GeV < 63.93 GeV |
| Codex: no discriminating empirical prediction yet | CONFIRMED | strengthened: no gravity-only observable within reach |
| Claude report: a wall signal would not distinguish C1 from other Z₂ sectors | CONFIRMED | GW §5; PBH channel generic |
| Claude report: GWs of early-decaying networks not computed | addressed [E] | capped at ≲ 2.5e-19 f⁻² (MeV/T_ann)², v-independent |

## 8. Limits

- **Q-ball family.** Flat-space, classical, D = 0, one branch. Beyond Q = 5.8e5 it is extrapolated (thin-wall form).
  Geometric contact cross-sections are [E]; inelastic merging and velocity dependence are not modelled.
- **Fragmentation.** Linear theory plus a γ ≳ 30H criterion. Fragment charge per fastest cell is [E]. Formation
  velocities are not computed; brackets are used, and the energy bound gives β ≤ 0.55–0.84.
- **Background and data.** g*(T) is a step table [R]. Detector and astrophysical sensitivities are recalled [R], not
  fetched, and are used only as orders of magnitude.
- **GW cap.** Uses the HKS radiation-era relation with A = 0.8 and ε̃ = 0.7 [Q]. The decay spectrum of a
  population-biased (pressure-free) network is not computed. f is [P] in [0.01, 1].
- **PBH channel.** Quoted from minority_relic. No black-hole formation fraction and no PBH limits were computed.
- **Out of scope.** Inflation above percolation and the n_c–n_s1 window are not covered here (unresolved in the
  packet).
- **No likelihoods.** No likelihood was computed. Every "gap" is an order-of-magnitude ratio, not a significance.
- **Data provenance.** All data are retrospective, and no held-out comparison was made.
