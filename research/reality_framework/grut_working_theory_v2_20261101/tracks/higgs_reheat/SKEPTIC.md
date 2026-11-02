# SKEPTIC of track higgs_reheat (Higgs survival in M1r with ξ_H = 0.1)

2 November 2026 · Claude, adversarial skeptic, decisive pass (GitHub issue #13 context). This is same-session AI work
on one supplied model. It is not outside replication, not a validation and not an exclusion. Code is in `skeptic/code/`
and output in `skeptic/results/`. The SM input (2-loop λ_eff and the anomalous dimension Γ from cf_vacuum's pipeline)
was reused, since the brief says not to redo it. Everything else was recomputed with a different formulation or a
different code.

Labels: [C] computed here, [E] estimate, [A] argued, [Q] quoted.

## 0. One-screen result

| Track result | Skeptic verdict | Recomputed |
|---|---|---|
| R1 SM barrier, Λ_I | not re-derived (brief: do not redo) | — |
| R2 ξ(M_t) = 0.1 → ξ(H_i) = 0.0878 | **confirmed** [C] | own 1-loop RGEs with the Buttazzo matching: factor 1.189 against 1.1837, so ξ(H_i) = 0.0874 (0.5%) |
| R3 inflation: m² = 0.937 H² "one loop"; exact variance 0.0322 | **mislabelled / partly wrong** [C] | 0.937 = 12ξ e^{2Γ}, the tree value in the M_t-normalised field; the genuine one-loop shift is only +0.0006. The canonical mass is 1.055 H². The exact super-horizon variance is 0.0280 H² (canonical) or 0.0315 H² (h₀ units), not 0.0322. B and h_b stand |
| R4 ξ_min (inflation) | **confirmed** [C] | own tree barrier search: B = 2096 / 2016 / 1461 / 1418 / 1113 / 1087 and ξ_EW,min(tree) = 0.0455–0.0548, identical. B does not depend on the field normalisation |
| R5 end-of-inflation patches above h_max | not refuted | Fokker–Planck tail is normalisation-invariant; light-formula x = 53.7 reproduced |
| R6 background | **confirmed** [C] | independent e-fold-time integration: R_min = −1.0636 at τ = 0.704; equality at τ = 1.786; T(0.05) = 4.03 / 4.36×10¹³; T_max = 1.59 / 1.72×10¹⁴ GeV at τ = 0.81 |
| R7 nothing grows at ξ = 0.088 | not refuted | own conformal-time modes: σ²(k<1, τ=0) = 0.0280 against 0.0287 |
| R8 gravitational-only lower edges | **shifted** [C] | 0.077 / 0.082 / 0.103 / 0.107 / 0.047 (canonical), against 0.074 / 0.079 / 0.100 / 0.104 / 0.043. Criterion band ±7% (region size, window), not < 3% |
| R9 benchmark D₆₀ | **refuted at +1σ** [C] | −185 / −132 / **+13.2** / **+26.6** / −928 (own code). The track's own code with the normalisation fixed gives −182 / −131 / **+12.5** / +25.6 / −912; at (172.87, 1.29e11), cH = ½ and 2 give +6.4 / +5.8 |
| R10 plasma must appear by τ ≈ 15 | **refuted** [C] | ≈ 8 (worst corner) and ≈ 10 (172.87, 1.29e11), for a sharp onset with consistent normalisation. The track's own code crosses D = 0 at τ_th ≈ 13, not 15; its tanh ramp acts about 1.8 H_i⁻¹ before τ_th |
| R11 thermal S3/T = 581–896 | **values refuted, conclusion stands** [C] | the same e^{4Γ} mismatch: S3/T = 478–738 (−18%). Bubbles 10⁻¹²⁰ (172.87) and 10⁻¹⁶⁰ (172.60), not 10⁻¹⁶⁶ and 10⁻²¹³. Still negligible against 215 |
| R12 verdict | **partly refuted** | see §6 |

## 1. The main error: two field normalisations mixed [C]

**The two normalisations.**
- `smpot.lam_eff_flat(h)` multiplies by e^{4Γ(h)}: V = λ_eff h₀⁴/4, with h₀ the field renormalised at M_t (Buttazzo
  convention).
- So h_max = 2.69×10¹¹ GeV and V″(h₀) are in h₀ units.
- The canonically normalised fluctuation field at μ ≈ H_i is φ = e^{Γ(H_i)} h₀, with e^{Γ} = 0.941–0.942. Its kinetic
  term in h₀ is Z = e^{2Γ} = 0.887.

**What modes.py does.** It evolves canonical modes with mass ξR, which is the correct mass for φ. It then compares
their variance with h_max in h₀ units, and it uses d²V/dh₀² as the tachyonic rate.

**The consistent version.**
- x = φ_max²/(2σ_φ²) is smaller by e^{2Γ} = 0.887.
- m_c² = −d²V/dφ² is larger by 1.128.
- Both make the track's D too negative, that is, too safe.

**inflation.py.** It finds the same factor and calls it a loop effect. m² = 2ΔV/h₀² = 0.937 is 12ξ e^{2Γ} = 0.9362 plus a
genuine one-loop +0.0006 (`sk_infl.py`). The light formula 0.0406 H² happens to be right in h₀ units, because
σ² ∝ 1/m². The "exact" 0.0322 uses ν(0.937) and is neither normalisation.

**thermal.py.** It has the same mismatch: a canonical bounce field carries the quartic λ_eff e^{4Γ} (e^{4Γ} = 0.79–0.81),
together with masses g²u²/4. Removing e^{4Γ} lowers S3/T by 17–18% (`sk_thermal.py`). The own Bessel-K thermal
functions and own shooting reproduce the track's S3/T to 5 digits in the track's convention. Derrick holds to ≤ 6×10⁻⁹.

**Benchmark ξ_EW = 0.1, N = 60, gravitational only** (`sk_modes.py`: conformal time, u = a h, scipy Hankel BD start,
own background):

| (M_t, H_i) | track's D | own code, track normalisation | own code, canonical | track code with the fix (`modes_canon.py`) |
|---|---|---|---|---|
| 172.60, 1.29e11 | −229 | −229.4 | −185.1 | −182.2 |
| 172.60, 1.51e11 | −168 | −167.0 | −132.0 | −131.2 |
| 172.87, 1.29e11 | +1.1 | +2.5 | **+13.2** | **+12.5** (cH = ½: +6.4; cH = 2: +5.8) |
| 172.87, 1.51e11 | +18 | +18.1 | **+26.6** | +25.6 |
| 172.33, 1.51e11 | −1087 | −1090.5 | −928.0 | −911.8 |

The independent code reproduces the track in its own convention to within ΔD ≤ 3.8 and the lower edges to ≤ 0.002. The
correction therefore comes from the normalisation, not from code differences.

## 2. Plasma deadline [C]

D(τ_th) is the danger exponent counting only regions whose roll begins before a sharp plasma onset (`sk_modes.py criteria`):

| | track normalisation | canonical |
|---|---|---|
| 172.87, 1.51e11 | zero at τ ≈ 11 | **zero at τ ≈ 8** (+0.4 at 8, −34 at 5) |
| 172.87, 1.29e11 | zero at τ ≈ 15.5 | **zero at τ ≈ 10** (−0.9 at 10, +6.2 at 12) |

**The track's own code** (`results/track_code_plasma_deadline.json`) at the worst corner:
- D = −4.9 at τ_th = 12 and +4.6 at 14, so the crossing is at ≈ 13. The quoted "≈ 15" was a linear interpolation
  between τ_th = 10 and 20.
- The tanh ramp, (1 + tanh((τ − τ_th)/0.2))/2 multiplied by m_T² ~ 10⁵ H_i², already exceeds |V″| ~ 10⁻³ H_i² at
  τ_th − 1.7.

**Corrected statement.** A Higgs-coupled population must give m_p² ≳ 10⁻³ H_i² by τ ≈ 8–10 (3.5–5×10⁻³⁵ s), not by 15.

**The magnitude needed is tiny [E].** For the thermal mass to act at h ≈ h_max one needs T ≳ g h_max/2 ≈ 3×10¹⁰–10¹¹ GeV
in Higgs-coupled states. That is about 10⁻¹² of ρ_tot at τ ≈ 8. Any prompt branching into SM states suffices; the
issue is timing, not amount.

## 3. The criterion's own spread at +1σ [C]

D is steep in ξ (dD/dξ_EW ≈ 3×10³). The O(1) conventions inside the argued criterion therefore move D by tens:
- region wavenumber k_c = f a m_c, with f = ½, 1, 2;
- a Gaussian window instead of a sharp cut.

| | D range, canonical |
|---|---|
| 172.60, 1.29e11 | −342 to −123 |
| 172.60, 1.51e11 | −257 to −83 |
| 172.87, 1.29e11 | −26 to +29 |
| 172.87, 1.51e11 | −5 to +39 |

The canonical lower edges (`sk_edges_band.py`) range as follows:
- 0.072–0.080 (172.60, 1.29e11);
- 0.076–0.085 (172.60, 1.51e11);
- 0.096–0.108 (172.87, 1.29e11);
- 0.099–0.112 (172.87, 1.51e11).

So central M_t is safe under every variant. The +1σ corners are undetermined by this criterion: they fail for the
track's own choice (f = 1, sharp) and pass only for f = ½. The track's "ξ_EW ≥ 0.11 removes the marginal corner in
every case" fails for f = 2 at (172.87, 1.51e11), where the edge is 0.112.

**Inflaton history [C].** Γ_I = 0.1–2 H_i at m_I = 3 moves D by −4.5 to +1.3 at +1σ. A slower decay, with a longer
matter era and ⟨R⟩ > 0, is slightly safer. No value of Γ_I tried rescues the +1σ corners (D = +8.7 at Γ_I = 0.1).

## 4. Confirmations [C]

- **ξ running.** Own 1-loop SM RGEs from the Buttazzo et al. 1307.3536 fits (not the SMDR grid) give the factor
  1.1886–1.1892 against 1.1834–1.1838.
- **Background.** Different variables (e-fold time, P = dI/dN, H from the constraint) give the same R(τ), equality
  and T values.
- **Hawking–Moss tree exponents and ξ_min(tree).** Reproduced exactly. B is independent of the normalisation, so
  R4's conclusion that "inflation is safe at every corner, ξ_EW,min = 0.044–0.049" stands.
- **κ_T formula.** 3g²/16 + g′²/16 + y_t²/4 + λ/2 is reproduced from the small-y limit of the J functions.

## 5. New issues

1. **[C] Normalisation mismatch** (§1). It affects modes.py, edges.py, controls.py and thermal.py. Fix: put h_max and
   V″ into φ = e^{Γ}h₀ units, or evolve h₀ with Z = e^{2Γ}. Upper edges (19–26, mean field 2.2–11) shift too; not
   recomputed here.
2. **[C] Plasma deadline** τ ≈ 8–10, not 15 (§2). The tanh ramp tail should not be read as the onset time.
3. **[C] The +1σ status is criterion-dominated** (§3). "Marginal / fails" should read "fails under the track's
   criterion; undetermined within ±40 in D".
4. **[A] M1r as written has no renormalizable channel that heats Higgs-coupled SM states.**
   - The reflection-even inflaton forbids φ_I|H|² (WORKING_THEORY §2), and the sector is declared to have "zero SM
     portals".
   - A gauge-singlet φ_I has no renormalizable Yukawa to the chiral SM fermions. The "fermionic channel" must go to
     singlet fermions (for example N), or through dimension 5.
   - For N with seesaw-compatible Yukawas (y² = 2m_ν M_N/v² ≈ 3×10⁻⁴ at M_N ≈ 1.9×10¹¹ GeV, m_ν = 0.05 eV), the
     N lifetime is about 5×10⁴ H_i⁻¹ ≫ 8.
   - So the prompt-plasma rescue needs an unstated extension. The gravitational-only case is the default for M1r as
     stated.
5. **[A] The supplied perturbative (3, 0.5) decay is not what any Yukawa channel does at these amplitudes.**
   - Γ_I/m_I = 1/6 needs N_f y² = 4.2.
   - With φ_I(0) = 0.82 M_P, the fermion mass y φ_I ≈ 4×10¹⁸/√N_f GeV ≫ m_I = 3.9×10¹¹ GeV. The decay is kinematically
     closed except near zero crossings, so it is non-perturbative.
   - ρ_I/m_I⁴ ≈ 1.3×10¹³: soft fermionic products would be Pauli-blocked by about 10¹⁵ unless they thermalise
     instantly.
   - The engine's T = 4×10¹³ GeV at τ = 0.05 is therefore an idealisation of the supplied history. D depends on Γ_I
     only weakly (§3); the plasma timing depends on the channel strongly.
6. **[Q/A] Literature framing.** The standard end-of-inflation patch count (Kearney–Yoo–Zurek 1503.05193; Espinosa et
   al. 1505.04825, who also note thermal rescue of overshooting patches) gives "fails" at M_t ≥ 172.6 without plasma.
   The track's safe-at-central-M_t verdict without plasma rests on its own argued, scale-resolved criterion, which I
   find physically better motivated. The track reports both readings.
7. **[C] Stated numerical controls understate the uncertainty.** "m_c-onset variants move the lower edge < 3%"
   holds, but region-size and window conventions move it ±7%.

## 6. Corrected verdict for ξ_EW = 0.1

**During inflation: safe** at every corner. This is unchanged and confirmed.

**Oscillating-R era: nothing grows.** Unchanged.

**Reheating, gravitational only** (no Higgs-coupled plasma before τ ≈ 10):
- **Safe at M_t ≤ 172.6.** D ≤ −83 under every convention tried.
- **Fails at M_t = 172.87 for both H_i calibrations** under the track's own criterion with consistent normalisation:
  D = +13 (free-X calibration) and +27 (Q-matter calibration).
- **Undetermined within the criterion's conventions** at that M_t: D from −26 to +39.

**With a plasma:** safe if Higgs-coupled states with T ≳ 10¹¹ GeV exist by τ ≈ 8–10 H_i⁻¹. M1r does not supply such a
channel at the renormalizable level [A].

**Thermal decay:** negligible (S3/T ≥ 478).

ξ_EW ≈ 0.11–0.12 would clear every variant computed here at +1σ. That is an observation, not a proposed change.

## 7. Files

All under `research/reality_framework/grut_working_theory_v2_20261101/tracks/higgs_reheat/skeptic/`.

**Code (`code/`)**
- `sk_modes.py`: independent modes, background and danger function. Modes: `bench`, `edges <i>`, `criteria`,
  `history`.
- `modes_canon.py`: a copy of the track's modes.py with the normalisation fix; edits are marked `# SK`.
- `sk_edges_band.py`: lower edges under the criterion variants.
- `sk_infl.py`: ξ running, the curvature-mass decomposition, variances and tree Hawking–Moss.
- `sk_thermal.py`: Bessel-K thermal functions and shooting, in both normalisations.
- `sk_bubbles.py`: thermal bubble counts.
- `sk_summary.py`: collects the headline numbers.

**Results (`results/`)**
- `sk_modes.json`, `sk_edges_*.json`, `sk_edges_band_*.json`, `modes_canon_benchmark.json`
- `sk_infl.json`, `sk_thermal.json`, `sk_bubbles.json`
- `track_code_plasma_deadline.json`
- `sk_summary.json`

**Reproduction** (in `skeptic/code`, one thread each; about 6 minutes in total):

```
python3 -B sk_modes.py bench
python3 -B sk_modes.py criteria
python3 -B sk_modes.py history
python3 -B sk_modes.py edges 0 1
python3 -B sk_modes.py edges 2
python3 -B sk_modes.py edges 3
python3 -B sk_modes.py edges 4
python3 -B sk_edges_band.py 0 1
python3 -B sk_edges_band.py 2 3
python3 -B sk_infl.py
python3 -B sk_thermal.py
python3 -B sk_bubbles.py
python3 -B sk_summary.py
```

`run_modes_canon.py` runs `modes_canon.py` at the five corners and writes `modes_canon_benchmark.json` (about 1 minute).
