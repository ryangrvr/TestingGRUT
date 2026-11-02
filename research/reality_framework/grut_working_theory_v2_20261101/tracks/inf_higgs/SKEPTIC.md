# SKEPTIC of track inf_higgs (key inf_higgs_skeptic)

3 November 2026 · Claude, adversarial pass on `tracks/inf_higgs` (BRIEF_INFLATON.md; GitHub issue #13 context). This is
same-session AI work on one supplied model. It is not outside replication. Nothing here is a validation, an exclusion of
the real world or a new theory. **Every choice in the track and here is a model choice:** ξ_H, the benchmarks, the decay
history, the S1 channel (y, M_N), the A-term and portal renormalisation, and the survival criterion. None follows from
closure. Git was read-only and nothing was committed. Every run was in the foreground, with at most two processes.

Labels: [D] derived, [C] computed with my own code (`skeptic/code/`, output in `skeptic/results/`), [E] estimate,
[A] argued, [Q] quoted, [R] recalled.

**Method policy.**
- The SM input is shared on purpose: `higgs_reheat/code/smpot.py` (3-loop running, flat 2-loop λ_eff, Γ, SMDR matching),
  imported read-only. The brief says not to redo it.
- Everything else is re-written:
  - potentials from the formulas;
  - background pre-run in cosmic time with LSODA;
  - Higgs modes in the friction-free variable u_k = a^{3/2}h_k, with u'' + [k²/a² + M² − 9H²/4 − 3Ḣ/2]u = 0,
    co-integrated with the background by a vectorised fixed-step RK4;
  - a Bunch–Davies start from scipy Hankel functions;
  - my own implementation of the track's danger criterion D, and my own period mean.
- The track uses h_k with 3Hḣ, DOP853 and mpmath.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict on one screen

| Track result | Verdict | Recomputed | Label |
|---|---|---|---|
| R1 M1r regression | **stands** | own code: +12.46 (Hartree) / +13.3 (no Hartree; the higgs_reheat skeptic's code had +13.2) at (172.87, 1.29e11, M_h = 125.13); −186.8 / −185.7 at 172.60; τ_c 14.88 | C |
| R2 inflation with H(N) | **numbers stand**; the criterion is then violated inside the track's own claimed range (§3.2) | own background and barrier: B_end 684 / 430 / 310, B_const 2425 / 1580 / 1160; ξ_EW,min(tree) 0.0595 / 0.0701 / 0.0788 (constant H 0.043 / 0.048 / 0.053); one loop 0.067 at 172.86 | C |
| R3 Gaussian tail | conclusion stands; value is sampling-dependent | The tail rises steeply through the last 0.01 H⁻¹: log₁₀ = −761 at τ = 0.991 (R/H² = 7.85), ≈ −560 at ε_H = 1, −417 at τ = 1.011. The track's −698 is its last save before the end. 10⁻¹⁰² even at 174.6 GeV | C |
| R4 R < 0 windows | **stands** | R_min = −6.0 H²; windows last 0.026 / 0.010 / 0.005 H⁻¹ with exponents 0.015 / 0.0056 / 0.0028 (HT6t / T1 / E1); the surrogate gives −3.84 H², 0.224 and 0.102; fraction of time with R < 0 is 0.374 / 0.350 / 0.351 | C |
| R5 D at ξ_EW = 0.1 | **stands** | D₆₀ at 172.57 / 172.86: HT6t MD −374.1 / −24.2; T1 MD −371.6 / −24.2; E1 MD −371.9 / −24.1; HT6t G05 −314.4 / −5.7; surrogate −306.3 / −5.4 | C |
| R6 M_t,crit | **stands as computed**; it is conditional on α_s and on the input edition (§4, N5 and N7) | HT6t MD 172.917 (N = 60) / 172.937 (N_obs); T1 and E1 MD 172.919 / 172.938; HT6t G05 172.874 / 172.884; surrogate 172.874. The HT6t background rescaled to 1.5×10¹¹ gives 172.873 (the track's recalibrated potential gives 172.871) | C |
| R7 deadlines | **stands** | MD: 16.5 (172.95), 10.85 (173.05), 8.55 (173.15), 5.62 (173.45). G05: 12.4 / 9.8 / 7.54 / 6.42 / 4.78 at 172.90–173.45 | C |
| R8 lower ξ edge | **stands** | D₆₀(ξ_EW = 0.095, 172.86) = +3.8; D₆₀(0.067, 172.57) = +0.3 | C |
| R9 upper edge | not rechecked | — | — |
| R10 A bound during inflation | **stands as an inflation-only number, but it is not the binding bound** (§3.3) | own count: 17.1 / 12.1 TeV (HT6t, 172.57 / 172.86), 42 / 39 (T1), 45 / 41 (E1). The track's T1 and E1 values are ~15% conservative (its plateau tail uses φ_I at N = 50) | C |
| R11 S1 plasma timing | stands at the model's [E] level | Analytic rate dm²/dt ≈ c y_ν² n_GK/2π gives τ(10⁻³) ≈ 1.5 for HT6t at y = 3×10⁻⁵, M_N = m_I/4 (track 1.10) [E]. Composition: a pure L + H gas gives c = 0.031–0.034 through the longitudinal W/Z, 1.4–1.5× the track's 0.0228 (favourable) [A/C] | E |
| R12 D with the S1 plasma | **refuted** (§3.1) | The track's −inf entries are a truncation of its h̄ grid at 20 h_max. With the grid extended to 1000 h_max:<br>• HT6t, y = 3×10⁻⁵: crosses 0 at 173.78 GeV (M_N = m_I/16) and 174.52 GeV (m_I/4)<br>• T1, m_I/4, y = 10⁻⁵: crosses at 173.96<br>• E1, y = 10⁻⁵: crosses at 174.48 (m_I/4) and 173.78 (m_I/16)<br>The track lists all of these as "> 174.6". In addition, the inflationary H(N) count fails above 173.75 GeV (tree) / 174.15 GeV (one loop) | C |
| R13 natural A_ind vs the bound | **refuted** (§3.3, §4 N4) | Gravitational only, survival at +1σ needs \|Aφ_I(S1)\| ≲ 0.07 H_i²: A ≲ 1.9 TeV (HT6t), ≈ 5.3 TeV (T1, E1), 6–7× below R10. At R10's own bound, M_t,crit = 172.60–172.61 on all three benchmarks. With the S1 plasma, the "compatible" case (HT6t, y = 10⁻⁴, M_N = m_I/16, A = 7.7 TeV) gives M_t,crit = 173.50. The same loop's R-even piece is larger by yφ_I/2M_N and is not controlled | C/D |
| R14 R-even portal | **conclusion stands**; values are grid-truncated | With the extended grid at λ_min: D = −2992 (172.86) and −1441 (174.6), against "≤ −4200 / −3300". At 0.1 λ_min: −480 at 174.6, against "≤ −820" | C |
| R15 numerical controls | **stands** | the independent code agrees to ΔD ≤ 9 at 172.57 and ≤ 0.5 at 172.86 | C |

**Headline.**
- The gravitational-only numbers of the track are reproduced by an independent code to ±1 in D and ±0.002 GeV in M_t,crit.
- What does not hold is the plasma-era statement "safe to M_t ≥ 174.6 GeV". Its −inf entries are a truncated field grid.
  The track's own inflationary criterion fails from 173.75 GeV (tree) or 174.15 GeV (one loop), and a plasma cannot act
  during inflation.
- The S1 channel's induced Higgs couplings are more dangerous than the track's A-term "catch" says:
  - the A-term matters mainly through the end-of-inflation variance, which sets the post-inflation danger, so its bound
    is 6–7 times tighter than the Hawking–Moss one;
  - the R-even piece of the same loop is 10–10³ H_i² for y ≥ 10⁻⁴, with a sign set by an unstated renormalisation
    condition.

---------------------------------------------------------------------------------------------------------------------

## 1. What was recomputed, and how

| Check | Track | Skeptic | File |
|---|---|---|---|
| background | e-fold pre-run + DOP853, inf_common potentials | own potential formulas; LSODA pre-run from slow roll ~12 e-folds before the end, then RK4 | `sk_common.py` |
| modes | h_k with 3Hḣ, DOP853 (rtol 1e-8), mpmath BD | u_k = a^{3/2}h_k, fixed-step RK4 (160 steps per inflaton period), scipy BD; Wronskian error 6×10⁻⁸ at the end | `sk_common.evolve` |
| Hartree | 6λ e^{−4Γ}σ², adiabatic subtraction | same convention (it is a convention); also run without | `sk_common.evolve` |
| ⟨R⟩ | running mean over the local period from turning points | running mean over the period between successive maxima of φ_I | `sk_common.period_mean` |
| D | `ih_modes.danger` | own implementation from the higgs_reheat skeptic's definition; optional extended h̄ grid | `sk_common.danger` |
| inflation | `ih_inflation.py` | own background history, own vectorised barrier on the SM splines; one loop via smpot.V_dS | `sk_infl.py` |

**Cross-check of the criterion code** (`sk_xcheck.py`). My danger function on the **track's** cached HT6t mode run, with
the track's plasma curve, gives −289.0 / −87.7 (y = 3×10⁻⁵, m_I/4, at 172.86 / 173.45). The track has −287.6 / −87.7. My
own mode run gives −286.9 / −88.0.

**A bug of my own, found and fixed.** My npz loader turned `frag = None` into 0.0, which switched off the period mean for
reloaded runs. All numbers here were produced after the fix. The fresh-run regressions (R5, R6, the A-term mode runs) were
never affected.

## 2. Reproductions (all [C])

| (M_t = 172.57 / 172.86) | track D₆₀ | skeptic D₆₀ | track M_t,crit (60 / N_obs) | skeptic |
|---|---|---|---|---|
| HT6t MD | −365 / −23.8 | −374.1 / −24.2 | 172.917 / 172.936 | 172.917 / 172.937 |
| T1 MD | −371 / −23.9 | −371.6 / −24.2 | 172.917 / 172.936 | 172.919 / 172.938 |
| E1 MD | −371 / −24.0 | −371.9 / −24.1 | 172.917 / 172.936 | 172.919 / 172.938 |
| HT6t G05 | −312 / −6.2 | −314.4 / −5.7 | 172.875 / 172.884 | 172.874 / 172.884 |
| M1r surrogate (Γ = 0.5) | −310 / −5.1 | −306.3 / −5.4 | 172.874 | 172.874 |

**Further reproductions.**
- D(τ_th) saturates by τ ≈ 30, so the τ ≤ 80 window is sufficient. At 172.95: −28.7 at τ_th = 10, +5.9 at 20, +11.2
  from 30 to 80.
- The R < 0 windows, the deadlines and the lower ξ edges are in §0.

## 3. Refutations

### 3.1 R12: the "−inf" entries are a truncated field grid [C] (`sk_plasmaD.py`, `results/sk_plasmaD.json`)

**Why the grid matters.**
- `ih_modes.danger` scans h̄ ∈ [h_max, 20 h_max] only.
- At higher M_t the barrier falls below the end-of-inflation spread. In canonical units, h_max is 0.905 H_i (172.86),
  0.159 (173.45), 0.039 (174.05) and 0.013 (174.6), while σ(k < aH) ≈ 0.16 H_i.
- Regions with h̄ ≫ 20 h_max are therefore abundant (x ≪ 3N).
- When none of the scanned h̄ rolls, the track returns D = −inf and reads it as safe. Larger h̄ do roll.

| D₆₀, S1 plasma (track's curves) | 172.86 | 173.45 | 174.05 | 174.6 | M_t,crit, track grid | **extended grid** |
|---|---|---|---|---|---|---|
| HT6t m_I/4, y = 3×10⁻⁵ | −288 | −88 | −24.5 (track −inf) | +3.3 (−inf) | > 174.6 | **174.52** |
| HT6t m_I/16, y = 3×10⁻⁵ | −192 | −29 | +15.4 (−inf) | +33.6 (−inf) | > 174.6 | **173.78** |
| T1 m_I/4, y = 10⁻⁵ | −219 | −43 | +5.8 (−inf) | +26.4 (−inf) | > 174.6 | **173.96** |
| E1 m_I/4, y = 10⁻⁵ | −292 | −87 | −22.6 (−inf) | +4.8 (−inf) | > 174.6 | **174.48** |
| E1 m_I/16, y = 10⁻⁵ | −191 | −29 | +15.7 (−inf) | +33.7 (−inf) | > 174.6 | **173.78** |
| HT6t m_I/4, y = 10⁻⁴ | −417 | −165 | −77 | −37.5 | > 174.6 | > 174.6 |
| y ≥ 10⁻⁴, all benchmarks | | | | < 0 | > 174.6 | > 174.6 |

The three crossings the track does report (HT6t y = 10⁻⁵: 173.39 and 173.13; T1 m_I/16: 173.44) are reproduced
(173.39 / 173.13 / 173.44).

### 3.2 The track's own inflationary criterion caps the plasma claim [C] (`sk_infl.py`, `results/sk_infl.json`)

The track evaluated the H(N) Hawking–Moss count only at 172.28 / 172.57 / 172.86. Extended (HT6t, N₀ = N_obs + 1 = 55.8;
T1 agrees to 0.0002):

| M_t | 172.86 | 173.15 | 173.45 | 173.75 | 174.05 | 174.6 |
|---|---|---|---|---|---|---|
| ξ_EW,min, tree / one loop | 0.079 / 0.067 | 0.086 / — | 0.094 / 0.084 | **0.100** / 0.091 | 0.106 / 0.098 | 0.116 / **0.109** |
| log₁₀⟨N⟩ at ξ_EW = 0.1, tree / one loop | −63 / −139 | −34 / — | −14 / −41 | **+0.1** / −19 | +10 / −3.9 | +23 / **+14** |

**What this means.**
- Under the criterion the track adopts for inflation, ξ_EW = 0.1 fails during inflation from **173.75 GeV (tree)** and
  **≈ 174.15 GeV (one loop)**, i.e. +4.1σ and +5.4σ (PDG 2024 σ = 0.29).
- A plasma cannot act during inflation. So "safe to M_t ≥ 174.6 GeV with the S1 plasma" (digest, R12, verdict) is
  inconsistent with R2.
- The non-equilibrium Gaussian count would pass to 174.6 (10⁻¹⁰²). Using it, however, abandons R2's ξ_min as well.

### 3.3 The A-term: the inflationary bound is not the binding one [C] (`sk_aterm.py`, `results/sk_aterm*.json`)

Mode runs on the MD backgrounds with m_H² += Aφ_I(t), present in the modes and, as a period mean, in the criterion.
Destabilising sign: Aφ_I < 0 during inflation.

| gravitational only | A (TeV) | Δm²(S1)/H² | D₆₀(172.57) | D₆₀(172.86) | M_t,crit |
|---|---|---|---|---|---|
| HT6t | 0 / 0.8 / 1.5 / 3 / 6 / 7.7 / **12** | 0 / −0.030 / −0.056 / −0.112 / −0.224 / −0.287 / −0.448 | −374 / −330 / −296 / −231 / −131 / −88 / **−12** | −24 / −13 / −3.9 / +12.6 / +36.6 / +46.3 / **+62.5** | 172.917 / 172.891 / 172.870 / 172.827 / 172.748 / 172.698 / **172.601** |
| T1 | −5 / **−34** | −0.070 / −0.477 | −290 / −18 | −2.8 / **+60.3** | 172.868 / **172.615** |
| E1 | **−35** | −0.458 | −15 | **+61.6** | **172.607** |
| HT6t, T1 (stabilising sign) | −12 / +34 | +0.45 / +0.48 | −1299 / −1307 | −531 / −380 | — |

**Reading.**
- At R10's own bound (12 / 34 / 35 TeV), every benchmark fails at +1σ with D ≈ +61, and M_t,crit drops to
  172.60–172.61.
- The D = 0 bound at 172.86 is |Aφ_I(S1)| ≈ 0.07 H_i², i.e. **A ≲ 1.9 TeV (HT6t), ≈ 5.3 TeV (T1, E1)**. At central M_t it
  is about 12.7 TeV (HT6t).
- **The mechanism is mainly inflationary.** T1 has no DC piece (⟨φ_I⟩ = 0) and shows the same failure. A lower Higgs mass
  in the last e-folds raises the super-horizon variance that the post-inflation criterion inherits; Hawking–Moss does
  not see that.
- **There is also a DC piece** on the asymmetric potentials. In the matter era ⟨φ_I⟩ ∝ Φ² ∝ H², so
  Aφ_I = κ_A H² with κ_A = −5.4×10⁻³ per TeV (HT6t) and +1.0×10⁻³ per TeV (E1), of the same sign as during inflation;
  T1 has none. At the R10 bound that is 25% (HT6t) or 14% (E1) of 3ξ.

**With the S1 plasma** (the track's curves, extended grid):

| case | A_ind (track's estimate, without / with log) | A used | M_t,crit |
|---|---|---|---|
| HT6t, y = 10⁻⁴, m_I/16 (R13 "compatible") | 0.24 / 7.7 TeV | 7.7 | **173.50** (A = 0: > 174.6) |
| HT6t, y = 3×10⁻⁵, m_I/16 | 0.07 / 2.3 | 3 | 173.40 (A = 0: 173.78) |
| HT6t, y = 3×10⁻⁵ / 10⁻⁴ / y_nat, m_I/4 | 1.2–34 / 3.9–112 / 11–319 | 12 | 172.95 / 173.35 / 173.78 |
| T1, y = 3×10⁻⁵, m_I/4 | 11–284 | −34 | 173.57 |
| T1, y = 10⁻⁵, m_I/16 | 0.2–6.5 | −34 / −5 | 172.72 / 173.30 (A = 0: 173.44) |

### 3.4 R14 values [C] (`sk_xcheck.py`)

The portal's conclusion (D ≪ 0) holds. Its quoted values come from the truncated grid; the maximum sits at h̄ > 20 h_max
even at 172.86.

| | track grid (20 h_max) | extended grid (1000 h_max) |
|---|---|---|
| λ_min, at 172.86 / 173.75 / 174.6 | −7282 / −21289 / −18385 | −2992 / −1464 / −1441 |
| 0.1 λ_min, at 172.86 / 173.75 / 174.6 | −2542 / −884 / −1326 | −1576 / −775 / −480 |

## 4. New issues

- **N1 [C] Field-grid truncation in the danger function.**
  - **Where.** `ih_modes.danger` (h̄ ≤ 20 h_max), and my own default.
  - **Symptom.** Whenever h_max ≲ σ (M_t ≳ 173.4 here), D = −inf is returned as safe.
  - **Fix.** Scan h̄ up to the abundance limit x = 3N (φ̄ ≈ (6Nσ²)^{1/2} ≈ 3–4 H_i), or flag −inf as "no roll inside
    the grid".
  - **Effect.** It changes R12 and R14 (§3.1, §3.4).

- **N2 [C] The plasma rescue is capped by inflation.**
  - The H(N) criterion fails at ξ_EW = 0.1 for M_t ≥ 173.75 GeV (tree) or ≳ 174.15 GeV (one loop) (§3.2).
  - The overall reach with S1 is therefore min(173.75–174.15, the post-inflation M_t,crit of §3.1/§3.3), not
    "≥ 174.6".

- **N3 [C] The A-term bound must come from the post-inflation D, not from Hawking–Moss.**
  - Gravitational only at +1σ: |Aφ_I(S1)| ≲ 0.07 H_i², i.e. 1.9 / ≈ 5.3 / ≈ 5.3 TeV (HT6t / T1 / E1) (§3.3).
  - With the S1 plasma it loosens, but A at the track's natural sizes still costs 0.3–1.5 GeV in M_t,crit.

- **N4 [D/C] S1's own loop induces an R-even Higgs mass that the track sets to zero.** (`sk_induced.py`)
  - **The term.** The N–ν mass matrix [[0, y_νh/√2], [y_νh/√2, M]] with M = |M_N + yφ_I| gives at one loop
    δm_h²(φ_I) = −(y_ν²/8π²) M² [ln(M²/μ²) − 1].
  - **Why it dominates.** Expanding M² gives the A-term the track bounds, plus an R-even λ_ind φ_I²|H|² that is larger by
    yφ_I/2M_N (≈ 65 at y_nat on HT6t). It does not depend on the CP branch.

  At S1 (H_i² units):

  | | c M² | best μ: max \|δm²\| over the last 5 e-folds | μ = M_N | μ = M_P |
  |---|---|---|---|---|
  | HT6t m_I/4, y = 10⁻⁵ / 3×10⁻⁵ / 10⁻⁴ / y_nat | 0.04 / 0.51 / 6.3 / 53 | 0.016 / 0.19 / 2.3 / 19 | −0.06 / −2.1 / −42 / −457 | +1.1 / +13 / +140 / +1050 |
  | HT6t m_I/16 | 0.015 / 0.14 / 1.6 / 13 | 0.005 / 0.053 / 0.60 / 4.9 | −0.07 / −1.0 / −15 / −152 | +0.4 / +3.5 / +36 / +266 |
  | T1 m_I/4 | 0.22 / 0.64 / 3.9 / 57 | 0.063 / 0.23 / 1.5 / 22 | +0.02 / −0.63 / −11 / −309 | +5.8 / +17 / +93 / +1200 |
  | E1 m_I/4 | 1.1 / 2.3 / 9.4 / 112 | 0.22 / 0.76 / 3.5 / 44 | +0.6 / −0.37 / −15 / −456 | +29 / +57 / +221 / +2370 |

  - **Under "vanish at μ = M_N"** (the analogue of the program's "vanish at M_Y"), the Higgs is tachyonic during
    inflation for y ≳ 3×10⁻⁵. That is fatal, whatever the plasma does.
  - **Under "vanish at M_P"**, the loop is a built-in R-even portal, which is stabilising.
  - **No choice of μ** keeps |δm²| below the ~0.07 H² that §3.3 tolerates, for y ≥ 3×10⁻⁵ at M_N = m_I/4, or for
    y ≥ 10⁻⁴ at m_I/16.
  - **Consequence.** The y range where S1 supplies the plasma in time (y ≥ 3×10⁻⁵) is the range where S1's own loop
    decides the Higgs during inflation. Every S1 result needs a stated renormalisation condition for λ_IH(μ), like the
    hidden-sector portals.

- **N5 [C] α_s is not a small systematic.**
  - α_s = 0.1180 ∓ 0.0009 moves M_t,crit to 172.594 / > 173.15 (MD; D(173.15) = −45 at 0.1189), and to
    172.560 / > 173.15 (G05).
  - At α_s = 0.1171, central M_t gives D = −7.1 (MD) and **+8.5 (G05)**.
  - That is ±0.33 GeV, 5× the M_h effect, and as large as the whole margin between central and critical M_t.
  - It is consistent with Buttazzo et al. 1307.3536 [R]: Δα_s = 0.0009 ↔ ΔM_t ≈ 0.43 GeV in Λ_I.
  - The track's "α_s not varied" hides it.

- **N6 [C] H_i must be calibrated before the gravitational-only verdict means anything.**
  - On the HT6t background with the SM at a different H_i, M_t,crit is:

    | H_i (GeV) | 1.29×10¹¹ | 1.5×10¹¹ | 2×10¹¹ | 3×10¹¹ | 5×10¹¹ | 10¹² |
    |---|---|---|---|---|---|---|
    | M_t,crit | 172.917 | 172.873 | 172.790 | 172.684 | 172.561 | 172.447 |

  - Central M_t fails for H_i ≳ 5×10¹¹.
  - inf_reheat (and its skeptic) find that slow decay pushes the dark-matter calibration of an M1r-like release to
    ~5×10¹² GeV [Q].
  - The MD history at H_i = 1.29×10¹¹ is therefore a fixed-H_i stand-in, not a consistent benchmark (the track says so in
    a caveat). If the calibration moves as inf_reheat found, the gravitational-only verdict flips at central M_t.

- **N7 [Q/C] The σ language depends on the input edition.**
  - The track switched to PDG 2024 (172.57 ± 0.29, M_h = 125.20). higgs_reheat used 172.60 ± 0.27 and M_h = 125.13.
  - PDG's 2025 summary has 172.56 ± 0.31 (direct), and a 2026 listing gives 172.60 ± 0.27 [Q, pdg.lbl.gov].
  - At (172.60 ± 0.27, 125.13): MD gives M_t,crit = 172.878 (+1.03σ) with D(172.87) = −2.9; G05 gives 172.838
    (+0.88σ) with **D(172.87) = +12.2**.
  - The "+1.2σ" headline is thus +1.0σ at the inputs the program used one pass earlier.
  - The MC-versus-pole top-mass interpretation (≈ 0.5 GeV; Hoang arXiv:2004.12915 [R]) is not counted either.

- **N8 [C/E] Fragmentation equation of state (T1, E1).**
  - The track's MDf assumes w = 0 exactly after fragmentation.
  - A radiation-like post-fragmentation fluid (RDf, w = 1/3: a "void"/scalar-radiation bracket) gives, for T1,
    D(172.86) = −3.0 and M_t,crit = 172.869 / 172.889 (N_obs).
  - So T1/E1 are bracketed at 172.87–172.92, not 172.913–172.917. Only a lattice can say where the oscillon gas sits.

- **N9 [C] The MD histories are not CMB-calibrated.**
  - At T_RH = 10¹¹ GeV: N_* = 50.35; n_s = 0.9625 (HT6t) and 0.9600 (T1, E1); A_s at fixed shape −8.5% / −8.8%
    (from inf_potential's tables).
  - Planck 2018 pulls are −0.58σ / −1.17σ; against P-ACT-LB (0.9743 ± 0.0034, as quoted by the inf_potential skeptic)
    they are −3.5σ / −4.2σ (retrospective, not exclusion).
  - r = 2.7×10⁻⁷ ≪ 0.036 is unaffected.
  - Restoring A_s raises H_i by ≈ 4% at fixed shape, i.e. ΔM_t,crit ≈ −0.01 GeV. This is minor.

- **N10 [A] CP-sign bookkeeping.**
  - If |A_ind| exceeds the D bound of N3, only the branch with Aφ_I > 0 survives. The observed asymmetry sign is then
    (α_I, sign φ_I) correlated with sign(A_ind) ∝ sign(yM_N): a third CP-relevant input the track does not state.
  - The R-even δm² of N4 is branch-blind, so it cannot be dodged by the branch choice.
  - HT6t has no R-symmetry at all (φ_I < 0 throughout inflation, and the potential is not even). A tree-level Aφ_I|H|²
    is then not forbidden, and "natural size" arguments for A do not apply there.

- **N11 [A/C] Hard-gas composition (favourable to the track).**
  - The track's caveat says an L + H gas couples only through g², g′² and λ.
  - In the region field h̄, the four Higgs-doublet quanta are one radial mode (3λh̄²) and three longitudinal W±, Z
    (equivalence theorem: m_W² = g²h̄²/4, m_Z² = (g² + g′²)h̄²/4); leptons couple only through y_ℓ.
  - So m_ind² = (n_SM/2E)·[3λ + (3g² + g′²)/4]/4, i.e. c = 0.031–0.034 at μ = 10¹¹–10¹³, against the track's
    thermal-normalised 0.0228.
  - Gauge bosons and tops are not needed for the sign.

- **N12 [D] Small conventions.**
  - The region count omits 4π/3 and the comoving radius of the last-scattering sphere (≈ 3.2/a₀H₀): together +2.1 in D.
  - The Mathieu q of the 2ω part of ξR is 9ξH²/2ω², not /4.
  - Both are immaterial.

## 5. Omitted-effects checklist (brief)

| effect | status |
|---|---|
| preheating of the Higgs by ξR | none: per-window exponent ≤ 0.015, q ≲ 10⁻⁴ (reproduced) |
| inflaton fragmentation | T1/E1 bracket 172.87–172.92 (N8); HT6t homogeneous |
| thermal / hard-gas corrections | composition favourable (N11). The thermal decay rate at T_RH ≈ 10⁹–10¹¹ GeV is not recomputed; higgs_reheat's skeptic found S₃/T ≥ 478 at higher T |
| induced couplings | A-term (N3); R-even S1 loop mass (N4); λ_IH convention as stated |
| kinematics, perturbativity | S1 blocking and open fractions taken from the track [E]. The portal has q ≤ 5×10⁻³ (reproduced from the formula) |
| CP-sign bookkeeping | N10 |
| M1r assumptions | H_i not recalibrated (N6); the decay history Γ_I ≈ 0 is not M1r's (3, 0.5). G05 equals the M1r surrogate (reproduced) |
| CMB | N9; the benchmarks' n_s and r themselves were checked by the inf_potential skeptic [Q] |

## 6. Corrected verdict for ξ_EW = 0.1 (model choices as in the track)

- **During inflation.** Safe at the PDG points under the track's H(N) Hawking–Moss criterion: ξ_min = 0.070 / 0.079
  (tree), reproduced. It fails from M_t ≈ 173.75 GeV (tree) / 174.15 GeV (one loop). The non-equilibrium Gaussian count
  is far safer (10⁻¹⁰² at 174.6).
- **After inflation, gravitational only.**
  - Safe at central M_t; at +1σ it passes by tens in D. M_t,crit = 172.917–172.938 is reproduced.
  - This holds only at H_i ≈ 1.3×10¹¹ (it falls to 172.56 at 5×10¹¹), at α_s = 0.1180 (±0.33 GeV) and with A = 0.
  - At the newer inputs it is +1.0σ (MD) and +0.9σ (G05).
- **With the S1 plasma.**
  - Post-inflation M_t,crit is 173.1–174.5 for y ≤ 3×10⁻⁵ and > 174.6 for y ≥ 10⁻⁴. The overall reach is capped at
    173.75–174.15 GeV by inflation.
  - This holds only if the S1 loop's A-term is ≲ 2–5 TeV, or the plasma compensates (M_t,crit 172.95–173.8 at
    A = 3–12 TeV).
  - It also needs the R-even S1 loop mass to be renormalised to a non-negative value over the observable trajectory.
    That fails for "vanish at M_N" when y ≳ 3×10⁻⁵.
- **With λ_IH ≥ λ_min.** Safe throughout the computed range (D ≤ −1441 at 174.6, extended grid). If the S1 loop is
  renormalised at M_P, it acts like such a portal.

## 7. Not established here

- The upper ξ edge (R9) and the Hartree runaway were not recomputed.
- The plasma model (GK refill, open fraction, hard-gas mass) was checked only at the order-of-magnitude level and for
  composition.
- No lattice: T1/E1 after fragmentation is bracketed, not computed.
- The post-inflation effect of the R-even S1 loop mass (an oscillating M(t) through m_N zeros) was not run through the
  mode code; only its inflationary size was computed.
- The SM input (smpot) is shared with the track. A different SM code was not used.

## 8. Files and reproduction

All under `research/reality_framework/grut_working_theory_v2_20261101/tracks/inf_higgs/skeptic/`. Runs are cached in
`/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/inf_higgs_skeptic/`. The plasma
curves and portal runs are read-only from the track's scratchpad (`.../wt2/inf_higgs/`).

| code | content → results |
|---|---|
| `sk_common.py` | SM splines (smpot read-only), own potentials, background, RK4 modes, period mean, danger |
| `sk_bench.py run <name> <hist>` | D(M_t), M_t,crit for HT6t/T1/E1 MD, HT6t G05, T1 RDf, M1r → `sk_bench.json` |
| `sk_aterm.py run / plasma` | A-term mode runs; A-term + S1 plasma → `sk_aterm.json`, `sk_aterm_plasma.json` |
| `sk_plasmaD.py` | S1 plasma D with the track's and the extended h̄ grid → `sk_plasmaD.json` |
| `sk_infl.py main / loop / Abound` | H(N) Hawking–Moss to 174.6 GeV, tree and one loop; inflation-only A bound → `sk_infl.json` |
| `sk_induced.py` | S1 loop-induced R-even Higgs mass along the trajectory → `sk_induced.json` |
| `sk_sens.py Hi / alphas / tauth` | H_i and α_s sensitivity, D(τ_th) saturation → `sk_sens.json` |
| `sk_inputs.py` | newer PDG inputs; N_*, n_s, A_s of the MD histories → `sk_inputs.json` |
| `sk_misc.py` | hard-gas composition; Gaussian tail vs M_t → `sk_misc.json` |
| `sk_repro.py` | R4 windows, R7 deadlines, R8 edges → `sk_repro.json` |
| `sk_xcheck.py` | criterion code on the track's run; portal with the extended grid → `sk_xcheck.json` |
| `sk_summary.py` | → `sk_summary.json` |

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/inf_higgs/skeptic/code; export OMP_NUM_THREADS=1
for j in "HT6t MD" "HT6t G05" "T1 MD" "T1 RDf" "E1 MD" "M1r G05"; do python3 -B sk_bench.py run $j 0; done   # ~10 min
for A in 0.8 1.5 3 6 7.7 12 -12; do python3 -B sk_aterm.py run HT6t $A; done
python3 -B sk_aterm.py run T1 -5; python3 -B sk_aterm.py run T1 -34; python3 -B sk_aterm.py run T1 34
python3 -B sk_aterm.py run E1 -35
python3 -B sk_aterm.py plasma HT6t 7.7 "MN=0.0625 y=0.0001"   # etc., see SKEPTIC §3.3
python3 -B sk_plasmaD.py; python3 -B sk_infl.py main; python3 -B sk_infl.py loop; python3 -B sk_infl.py Abound
python3 -B sk_induced.py
for H in 1.5e11 2e11 3e11 5e11 1e12; do python3 -B sk_sens.py Hi $H; done
python3 -B sk_sens.py alphas; python3 -B sk_sens.py tauth
python3 -B sk_inputs.py; python3 -B sk_misc.py; python3 -B sk_repro.py; python3 -B sk_xcheck.py; python3 -B sk_summary.py
```

**Literature** (abstract level unless stated):
- Buttazzo et al. 1307.3536 [R]
- Mantziris–Markkanen–Rajantie 2011.03763 and 2207.00696 [Q, via the track]
- Markkanen–Nurmi–Rajantie–Stopyra 1804.02020 [R]
- Kearney–Yoo–Zurek 1503.05193; Espinosa et al. 1505.04825 [R]
- Greene–Kofman hep-ph/9807339 [R]
- Lozanov–Amin 1608.01213 [R]
- Vissani hep-ph/9709409 (seesaw threshold correction to the Higgs mass) [R]
- Hoang 2004.12915 [R]
- PDG 2025 / 2026 top-quark listings, pdg.lbl.gov [Q]
- Planck 2018 X 1807.06211; BICEP/Keck 2110.00483 [R]
- ACT DR6 2503.14452 (as quoted by the inf_potential skeptic) [Q]

Our tracks, read-only: inf_higgs (code, results, scratchpad runs), higgs_reheat (smpot) and its skeptic,
inf_potential, inf_reheat and their skeptics.
