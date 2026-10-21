# Verifier `v_reach`: adversarial check of audit track `reach`

21 October 2026, issue #13. Verifies `audit/reach/` (NOTES.md, code/reach.py, code/qball_family.py, results/*).
Written only in `audit/v_reach/`. Nothing here excludes, validates or identifies Candidate 1. Labels: [C] computed
here, [D] closed form, [E] estimate, [P] preparation input, [Q] quoted, [R] recalled literature value (not fetched,
order of magnitude only), [RETRO] retrospective datum.

## Reproduce

```
cd research/reality_framework/grut_issue13_3d_review_20261021/audit/v_reach/code
OMP_NUM_THREADS=2 python3 v_qball_min.py ../results/v_qball_min_raw.json      # ~8 min, two resolutions
OMP_NUM_THREADS=1 python3 v_qball_large.py ../results/v_qball_large.json     # ~4 min, Q = 6.19e5, 1e6, 1.52e6
python3 v_summary.py ../results                                                # Richardson, threshold -> v_qball_summary.json
OMP_NUM_THREADS=1 python3 v_dispersion.py ../results/v_dispersion.json         # ~1 min
python3 v_cosmo.py ../results/v_qball_summary.json ../results/v_cosmo.json     # ~1 min
# Codex 479f111 component lists, extracted read-only with git show ... raw_diagnostics.tar.gz into the scratchpad:
python3 -I v_codex479.py <extracted raw dir> ../results/v_codex479.json
```

## Method (independent of the track)

| Track result | Track method | Verifier method |
|---|---|---|
| Q-ball family, binding threshold, archived soliton | `solve_bvp` continuation in Omega / F0 from the archived profile | **Direct constrained energy minimisation at fixed charge** (E[F,S] + Qbar²/2I, L-BFGS, analytic gradient), own tanh initial guess, cell-centred radial grid, h = 0.08 and 0.04, Richardson extrapolation. Does not read the archived profile. |
| Fragmentation growth rates | charge_phase cubic, `np.roots` | (A) Linearisation re-derived by hand into a **6×6 first-order system**, eigenvalues; (B) **direct 1D nonlinear lattice evolution** of the full field equations seeded at k_max; (C) closed-form NR-condensate limit. |
| Calibration, maser, T_n, PBH, σ/m, N_eff | `reach.py` | Own constants (s0 and Y_B recomputed), own script. |
| Free streaming, fragmentation floor, GW | `reach.py`, step-table g* | Own integral; **g* below 20 MeV computed** from Fermi–Dirac e± integrals with entropy-conserving T_ν/T; GW spectrum shape f³ / f⁻¹ around the peak [R: HKS 2014 / Saikawa 2017] added. |
| Codex 479f111 endpoint statistics | quoted counts | Recomputed from the raw component CSVs (read-only extract). |

## Results

### Q-ball family by constrained minimisation [C]

Richardson-extrapolated (h = 0.08, 0.04; the two differ by ≤ 1e-5 in e):

| Q | e = E/(m_X Q) | Ω/m_X | R90 [ℓ] (cell cumulant) |
|---|---|---|---|
| 1.50e4 | 1.01338 | 0.9503 | 6.22 |
| 1.65e4 | 1.00610 | 0.9207 | 5.42 |
| 1.74e4 | 1.00138 | 0.9094 | 5.25 |
| **1.77e4** | **0.99979** | 0.9061 | 5.21 |
| 2.4e4 | 0.96823 | 0.8581 | 4.93 |
| **39,944** | **0.9110861** | **0.80000** | 5.098 |
| 1.0e5 | 0.81802 | 0.7282 | 6.04 |
| 2.5e5 | 0.74687 | 0.6807 | 7.58 |
| 6.19e5 | 0.69575 | 0.6481 | 9.75 |
| 1.0e6 | 0.67477 | 0.6348 | 11.23 |
| 1.52e6 | 0.65927 | 0.6250 | 12.75 |

- **Binding threshold Q_b = 1.766e4** (e = 1 crossing). Track: 1.77e4. charge_phase thin-wall
  (m* + 3σ_Q/(n_* R) = m_X with σ_Q = 0.68, n_* = 1.778): 1.34e4. The track's correction is **CONFIRMED**.
- Archived soliton: e = 0.9110861 vs 0.91108652 (4e-7), Ω = 0.80000. R90 = 5.098 with this code's cell cumulant
  versus 5.1034 (convention difference ≈ 0.1%). **CONFIRMED.**
- Ω decreases monotonically with Q over 1.5e4–1.52e6 (dQ/dΩ < 0), a wider range than the track's 1.65e4–5.8e5.
  The branch continues below Q_b (Q = 1.5e4 is still localised, e = 1.013); its lower end was not located.
- Track e values agree: e.g. reach Ω = 0.68 → (2.54e5, 0.74573); here (2.5e5, 0.74687) at Ω = 0.6807. The track's
  thin-wall extrapolation e(6.185e5) = 0.69632 vs computed 0.69575 (8e-4, immaterial).

### Fragmentation (charge_phase cubic) [C]

| x | γ_max 6×6 | γ 1D lattice | k_max | track |
|---|---|---|---|---|
| 0.01 | 0.004498 | 0.004498 | 0.0943 | 0.0045 |
| 0.1 | 0.044849 | 0.044853 | 0.2820 | 0.0448 |
| 0.25 | 0.112140 | 0.112150 | 0.4035 | 0.1121, k = 0.4035 |
| 0.5 | 0.226400 | 0.226403 | 0.4605 | 0.2264 |

Small-x closed form [D]: a dilute φ ≠ 0 condensate is a nonrelativistic attractive Bose gas with
n dμ/dn = −(a² − b)x/2, so γ_max = (a² − b)x/2 = **0.450 x** and k_max = √((a² − b)x) = **0.9487 √x**. The 6×6
system gives 0.449998 x and 0.94869 √x at x = 1e-8. The track's 0.450 x / 0.949 √x is **CONFIRMED**.

### Calibration "hidden charge = all DM" [C, RETRO ω_c]

- s0 = 2891.2 cm⁻³, Y_B = 8.6734e-11 (proton-mass convention), Ω h² per GeV·r = 0.023799; **(E/Q)·r = 5.0422 GeV**.
- r at v = 1 TeV: free X (e = 1) 0.05042; homogeneous condensate at the Codex density (e = 0.9442 [Q charge_phase])
  0.05340; archived soliton 0.05534; Q = 6.19e5 (Codex fine-10% charge-weighted) 0.0725; Q = 1.52e6 (largest Codex
  lump) 0.0765; Q-matter floor (e = m*) **0.089664**. ω_asym = 1.33832 r. T_n = 442.8–536.5 TeV.
- Maser: m_X r = 5.04–8.97 GeV < 63.93 GeV; floor 100√(ω_b + ω_c) = 37.73 km/s/Mpc. **CONFIRMED.**

**Missing case in the band (finding V1).** In Candidate 1 with gravity-only coupling, binding energy released when a
condensate fragments stays in the hidden sector: h is stable (m_h = 141 GeV < 2 m_X, no portal), XX̄ → hh is closed
at rest (2 m_h > 2 m_X), and there is no massless hidden state. Only the kinetic part of the radiated quanta
redshifts. The e that enters the calibration is therefore the late-time hidden rest energy per net charge, which lies
between the final Q-ball value and the supplied state's value. For the Codex preparation (e_supply = 0.9442) this
gives r ∈ [0.0534, 0.0725] for the charge-weighted lumps, not 0.055–0.071. The overall band's upper edge 0.0897 is
rigorous (e_eff ≥ m* always); its lower edge 0.0504 holds only if the hidden density is the net charge alone — a
[P] supply that carries a symmetric XX̄ or h component pushes r below it.

### Codex 479f111 endpoint [C, from raw CSV]

Fine grid, 10% threshold: 39 components, top Qbar 12460.1/5775.0/4589.2/3454.5/3043.8, charge-weighted mean Qbar
6185.3 (Q = 6.19e5), 15 components above Q_b — all **reproduced**. Charge fraction is 0.690 (track 0.680); 5%: 0.758
(track 0.744). The charge-weighted mean is grid- and threshold-dependent: 3.6e5 (coarse 10%) to 1.07e6 (fine 5%), and
4.8e6 at 2% where one network holds 86%. Codex's caveat (not a mass function) stands.

### Cosmology channels [C, E]

| Quantity | Track | Verifier | Status |
|---|---|---|---|
| fragmentation floor T_Q (γ = 30H), v = 1 TeV | 0.35–0.60 GeV | 0.338–0.601 GeV | confirmed |
| λ_fs at floor, β = 0.5/0.9 | ≤ (1–6)e-6 Mpc | 9.9e-7–5.9e-6 | confirmed |
| λ_fs at T_Q = 4.4e5 GeV | 2e-12–8e-12 | 2.1e-12–7.6e-12 | confirmed |
| biased debris λ_fs (f = 1/0.1/0.01) | 4.9e-5 / 3.6e-3 / 0.32 Mpc | 4.9e-5 / 3.6e-3 / 0.28 | confirmed (g*) |
| σ/m contact, Q_b / archived | 4e-12 / 2e-12 | 4.21e-12 / 1.96e-12 | confirmed |
| σ/m XX / hh (threshold Born, M = 0.036 / 24 λ_φ) | 7e-16 / 1.1e-14 | 7.04e-16 / 1.11e-14 | confirmed |
| M_cut at λ = 6e-6 Mpc | 4.5e-6 M⊙ | 4.49e-6 | confirmed |
| ΔN_eff full conversion at T_Q = 1 TeV | 1.7e-12 | 1.67e-12 | confirmed (T_Q ≥ 1 TeV only) |
| ρ_DM/ρ_r at 1 MeV | 5.8e-7 | 5.83e-7 | confirmed |
| PBH M_min = 4π M_P⁴/σ | 4.20e9 M⊙ | 4.203e9 | confirmed |
| r_t at H_inf = 210–560 GeV | 4e-25–5e-24 | 7.2e-25–5.1e-24 | confirmed (lower edge was for 150 GeV) |
| portal κ thermalising at T_c | ≳ 1e-6 | 1.05e-6 | confirmed [E] |
| GW cap coefficient | 2.5e-19 f⁻²(MeV/T)² | 2.47e-19 | confirmed (reach.py's text string says 1.6e-19: cosmetic) |
| benchmark Ω_gw h², f = 1 / 0.1 | 3e-22 / 2e-18 | 2.97e-22 / 2.20e-18 | confirmed |
| benchmark Ω_gw h², f = 0.01 (0.31 MeV) | 6.5e-14 at 3.0e-11 Hz | **3.5e-14 at 3.4e-11 Hz** | corrected: track's step table has g* = 5.3 at 0.31 MeV; computed g* = 8.98. Its 6.5e-14 exceeds its own relic cap (3.7e-14) at that T. |
| Ω_gw h² inside the PTA band | "cap ≤ 2.5e-17" | **≤ 4.8e-16 at 2e-9 Hz (f = 0.01); 6e-16 for the benchmark f = 0.01 network** | corrected: the track used peak-in-band networks only; the f⁻¹ UV tail of networks peaking at 3e-11 Hz dominates the band. General: ≤ 4.8e-20 f⁻² at 2e-9 Hz, 9.6e-21 f⁻² at 1e-8 Hz. |
| lensing: max Q-ball mass | 5e-42 M⊙ (per-cell estimate) | 2.9e-42–5.1e-42 (same estimate) | estimate confirmed; **causal ceiling at the floor = 1.7e-10–9.4e-10 M⊙**, inside microlensing reach |
| Q-ball radii, top of range (Q = 1e14) | "1e-13 cm" | **1.0e-12 cm** (R = 5.1e3 ℓ) | corrected (text slip; no consequence) |

## Findings

V1. Calibration band: see above. Codex-like lumps need r ≈ 0.053–0.077, not 0.055–0.071 (energy retention plus
    the actual charge-weighted Q = 6.2e5 and the largest lump Q = 1.5e6; the track used e(5e5) = 0.7065).
V2. GW in the PTA band is ≈ 20–40× larger than the track's "≤ 2.5e-17" once the f⁻¹ spectral tail of later-decaying
    networks is included (≤ 5e-16 at the benchmark). Still ≈ 6 decades below the NANOGrav-15 level (~1e-9 [R]).
    Relative to SKA-era forecasts [R, uncertain, quoted in the literature roughly 1e-16–1e-13] the margin at the
    tuned f = 0.01 corner is 0 to ~3 decades, not "≥ 4". For v below the benchmark (lower allowed T_ann) the
    v-independent cap with tail rises further (∝ 1/T_ann), so the v-independence of the cap does not mean a
    v-independent in-band ceiling.
V3. Free streaming "≥ 3.7 decades in the bulk of both classes" holds only for the cold class with a cold-condensate
    supply and β_i ≤ 0.9. Biased-class debris near T_safe is 2.8 (f = 1), 0.9 (f = 0.1) and −1 (f = 0.01) decades
    from 0.03 Mpc; the debris fraction of DM there is ≈ T_safe/T_ann, so relevance concentrates at the latest
    allowed decay. A charge supplied as energetic X quanta [P] (not a condensate) is unbounded by the track's
    estimate; it tests the supply, not the action.
V4. The r = 1 (Y_Q = Y_B) extension, v ≈ 50–90 GeV, is called "still unobservable". At its latest-supply corner
    (r = 1: v = 50.4 GeV with e = 1, v = 89.7 GeV with e = m*) the computed λ_fs is 1.2e-3–3.3e-3 Mpc for
    β_i = 0.5 (near the track's energy bound 0.55–0.84), up to 1.1e-2 Mpc for β_i = 0.9, and 2.0e-5–5.7e-5 Mpc for
    β_i = 0.01. That is 1–1.4 decades (β = 0.5) or 0.4 decades (β = 0.9) below 0.03 Mpc. The margin is set by the
    uncomputed formation velocity, so "unobservable" is likely but not established at this corner.
V5. Lensing "≥ 25 decades" is an [E] resting on the per-fastest-cell fragment charge with no hierarchical merging.
    The only hard ceiling computed (all DM in one horizon at the fragmentation floor) is 1.7e-10–9.4e-10 M⊙, which
    lies inside microlensing reach. The estimate is physically plausible (merging by ~25 decades is implausible),
    but it should be labelled an estimate, not a bound.
V6. SIDM gap: ≥ 10.4 decades is relative to the current 0.1 cm²/g bound; relative to forecasts of ~0.01 cm²/g [R] it
    is ≥ 9.4. M_cut ≥ 11 decades is relative to current substructure probes (~1e6 M⊙ [R]).
V7. N_eff ≤ 1.7e-12 and the fragmentation-GW bound are quoted for T_Q ≥ 1 TeV; supply can be as late as the
    fragmentation floor (or later, homogeneous). At T_Q = 0.5 GeV ρ_h/ρ_r ≈ 1e-9: still negligible, and the
    hidden spectrum has no massless state, so no conclusion changes.
V8. Internal: reach.json's "Codex endpoint" row uses e(5e5) while the text quotes Q ≈ 6.2e5; the GW "scaling" string
    says 1.6e-19 while NOTES says 2.5e-19 (the computed rows give 2.47e-19).

## Verdicts

| # | Claim (track `reach`) | Verdict | Verifier value / reason |
|---|---|---|---|
| 1 | relic_yield retrospective r_max = 0.0897 (omega_h <= 0.12), CONFIRMED by reach | **CONFIRMED** | 0.089664. Own s0 = 2891.2 cm^-3, Y_B = 8.6734e-11: (E/Q) r = 5.0422 GeV; with e = m*/m_X = 0.56234, r = 0.089664; omega_asym = 1.33832 r. It is a rigorous upper edge because any hidden state has e_eff >= m*. |
| 2 | reach: all-DM calibration fixes only e*gX*v*r = 5.042 GeV; r in [0.0504, 0.0897]; Codex-like cold Q-balls need 0.055-0.071 (correcting the 'r ~ 0.09' framing) | **CORRECTED** | r in [0.0504, 0.0897]; Codex-like 0.053-0.077. Product and band edges reproduced (5.0422 GeV; 0.05042-0.08966), so 'r ~ 0.09' is indeed only the Q-matter edge. But (i) released binding energy stays as stable hidden h / X quanta (m_h < 2 m_X, XX->hh closed, no massless state), so e_eff lies between the Q-ball value and the supplied state's e (0.9442 at the Codex density); (ii) the charge-weighted Codex lump is Q = 6.19e5 (e = 0.6958, r = 0.0725) and the largest Q = 1.52e6 (e = 0.6593, r = 0.0765), while the track used e(5e5). Codex-like range is 0.053-0.077. The lower edge 0.0504 holds only if the hidden density is the net charge alone. |
| 3 | reach: thin-wall binding threshold Q ~ 1.3e4 [E] (charge_phase) corrected to 1.77e4; dQ/dOmega < 0 over 1.65e4-5.8e5 | **CONFIRMED** | 1.766e4. Independent fixed-charge energy minimisation (no archived profile, two resolutions, Richardson): e = 1 at Q = 1.766e4; thin-wall formula gives 1.34e4. Omega decreases monotonically over Q = 1.5e4-1.52e6. |
| 4 | Archived soliton Q = 39,944, E/(m_X Q) = 0.91108652, R90 = 5.1034 | **CONFIRMED** | e = 0.9110861, Omega = 0.80000. Minimisation at fixed Qbar = 399.44 gives e = 0.9110861 (4e-7 off) and Omega = 0.80000; R90 5.098 with a cell-based energy cumulant (convention difference ~0.1%). |
| 5 | Branch-1 fragmentation gamma_max = 0.0045/0.0448/0.1121/0.2264 at x = 0.01/0.1/0.25/0.5; k_max = 0.4035; small-x law 0.450 x, 0.949 sqrt(x) | **CONFIRMED** | 0.004498/0.044849/0.11214/0.22640; k = 0.40351; 0.45 x, 0.9487 sqrt(x). Linearisation re-derived by hand as a 6x6 first-order system (eigenvalues) and checked by direct 1D nonlinear lattice evolution of the full field equations (rates agree to 1e-5). Small-x closed form: attractive NR condensate, gamma = (a^2-b) x/2, k = sqrt((a^2-b) x) = 0.9487 sqrt(x). |
| 6 | Codex maser product m_X|r| < 63.93 GeV and H0 floor compatible with hidden charge = all DM; not a separate test once omega_c calibrates | **CONFIRMED** | 5.04-8.97 GeV; 37.73 km/s/Mpc. m_X r = 5.04-8.97 GeV; floor 100 sqrt(0.02233+0.120) = 37.73 km/s/Mpc. The floor and the maser envelope bound omega_h only from above. |
| 7 | Codex a69bd86 REVIEW: the model cannot yet make a discriminating empirical prediction because identified inputs or calculations are missing | **CONFIRMED** | REVIEW.md states this verbatim; nothing found here yields a Candidate-1-specific observable. Every channel that comes near reach (biased-class debris free streaming, PTA-band GW tail, latest-supply free streaming) is controlled by [P] inputs (f, m0, T_ann, supply form, beta_i) and is generic to Z2 walls. |
| 8 | reach strengthening: every gravity-only channel at v ~ 1 TeV sits 3.7-25 decades below current or planned reach (SIDM >= 10.4, lambda_fs >= 3.7, M_cut >= 11, lensing >= 25, N_eff <= 1.7e-12, GW <= 2.5e-17 in PTA band) | **CORRECTED** | cold class: 3.7-25 dec (estimates); biased class near T_safe: 2.8 to -1 dec in free streaming; GW tail <= 5e-16. Overstated as a uniform margin. Biased-class debris near T_safe: lambda_fs = 4.9e-5 / 3.6e-3 / 0.28 Mpc for f = 1/0.1/0.01, i.e. 2.8 / 0.9 / -1 decades from 0.03 Mpc. PTA-band GW including the f^-1 tail: <= 4.8e-16 (f = 0.01) at 2e-9 Hz, not 2.5e-17. Lensing >= 25 decades is an estimate from per-cell fragment charge; the hard causal ceiling at the fragmentation floor (1.7e-10-9.4e-10 Msun) lies inside microlensing reach. SIDM margin vs forecasts (0.01 cm^2/g) is 9.4 decades. Cold-class numbers (lambda_fs, sigma/m, M_cut, N_eff) reproduce for a cold-condensate supply. |
| 9 | reach estimate (Claude c7b8480 left GW of early-decaying networks uncomputed): Omega_gw h^2 <~ 2.5e-19 f^-2 (MeV/T_ann)^2; benchmark 3e-22 / 2e-18 / 6.5e-14 for f = 1/0.1/0.01; <= 2.5e-17 inside the PTA band | **CORRECTED** | cap 2.47e-19 f^-2 (MeV/T)^2; f=0.01 benchmark 3.5e-14; in-band <= 4.8e-16 (f=0.01). REPORT.md line 301 does say GWs were not computed. Cap coefficient reproduced (2.47e-19) and f = 1, 0.1 benchmarks reproduced (2.97e-22, 2.20e-18). f = 0.01: the track's step table has g* = 5.3 at 0.31 MeV (computed 8.98), giving 6.5e-14, which exceeds its own relic cap there; computed value 3.5e-14 at 3.4e-11 Hz. The PTA-band ceiling ignores the f^-1 UV tail of networks peaking below the band: <= 4.8e-20 f^-2 at 2e-9 Hz (benchmark f = 0.01 network: 6e-16). reach.py's text string says 1.6e-19. |
| 10 | Claude c7b8480: a wall signal would not distinguish Candidate 1 from other Z2 sectors | **CONFIRMED** | 4.203e9 Msun. GW amplitude and frequency depend only on sigma, A, eps and T_ann; the minority-domain mass scale 4 pi M_P^4/sigma (reduced Planck mass) = 4.203e9 Msun at v = 1 TeV. |
| 11 | reach: Q-ball mass function at v = 1 TeV spans Q ~ 1.8e4 -> 1e14, M 1.6e-51 -> 5e-42 Msun, radii 1e-15 -> 1e-13 cm | **CORRECTED** | R(1e14) = 1.0e-12 cm. Masses and Q range reproduce as an [E] (floor Q_frag = 5.7e13-1.0e14). The top radius is (3Q/(4 pi n*))^(1/3) = 5.1e3 ell = 1.0e-12 cm, not 1e-13 cm. No downstream consequence (contact sigma/m ~1e-15 cm^2/g there). |
| 12 | reach extension: Y_Q = Y_B fixes v ~ 50-90 GeV, which is still unobservable under gravity-only coupling | **CORRECTED** | lambda_fs 1.2e-3-1.1e-2 Mpc (beta 0.5-0.9). v = 50.4-89.7 GeV reproduced. At the latest-supply (fragmentation-floor) corner lambda_fs = 1.2e-3-3.3e-3 Mpc for beta_i = 0.5 and up to 1.1e-2 Mpc for beta_i = 0.9, only 0.4-1.4 decades below 0.03 Mpc. 'Unobservable' is likely but depends on the uncomputed formation velocity; not established. |
| 13 | reach: Codex 479f111 fine-grid 10% endpoint has 39 components holding 68% of the charge, 15 above Q_b, charge-weighted mean Q ~ 6.2e5 | **CONFIRMED** | 39; 15; 6.19e5; 0.690. Recomputed from the raw component CSVs: 39 components, 15 above Q_b, top Qbar and charge-weighted mean 6185.3 identical. Charge fraction 0.690 (quoted 0.680), 5%: 0.758 (0.744). Charge-weighted mean ranges 3.6e5-1.07e6 across grids/thresholds. |

## Limits

- Q-ball family: classical, flat space, D = 0, spherical, one branch, same benchmark couplings. Minimisation finds
  local minima from a tanh guess and continuation; it does not prove global minimality. R90 here uses a cell-based
  energy cumulant (≈ 0.1% different from the track's convention).
- 1D lattice growth check is one Fourier mode in a one-wavelength periodic box; it checks the linear rate only.
- g* above 20 MeV is a step table [R]; below, e± are computed with instantaneous ν decoupling at 2 MeV.
- GW spectral shape (f³ below, f⁻¹ above the peak) is recalled [R]; detector sensitivities are recalled [R], not
  fetched. No likelihoods. The relic condition uses κ = 1 (κ = 2 lowers every cap by 4).
- PBH, isocurvature and minority-domain numbers were only re-derived where closed form (M_min, r_t); the
  minority_relic abundance statements were not rechecked.
- All data are retrospective; no held-out comparison.
