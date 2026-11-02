# Track inf_potential: explicit inflaton-potential benchmarks for M1r

3 November 2026 · Claude, inflaton-sector pass (BRIEF_INFLATON.md of this directory; issue #13 context). This is
same-session AI work on one supplied model. It is not outside replication. Nothing here is evidence for an inflaton, the
hidden sector or the holder, and nothing is a validation or an exclusion. **Every potential, slice convention, holder
calibration and decay width below is a model choice. None of them follows from closure.** Data comparisons (Planck 2018
n_s, A_s; BICEP/Keck r) are retrospective.

Labels: [D] derived, [C] computed with our own code, [E] estimate, [Q] quoted from a reviewed track, [A] argued,
[R] recalled or literature.

Git was read-only and nothing was committed. Every run was in the foreground, with at most one worker process at a time.
One early scan exceeded the 2-minute default call limit. The harness moved it to the background, and I stopped it before
any further work and re-ran it in the foreground. The reheating track's T_RH did not exist when this was done. N_* is
therefore given for the engine's own width Γ_I = 0.5 H_i, and for stated (w_reh, T_RH) values.

## 0. Result in brief

**Benchmarks.** H_i = 1.29×10¹¹ GeV, the free-X calibration. The 1.50×10¹¹ (Q-matter calibration) rows are in
`results/summary.json`.

| | **HT6t** (tilted hilltop, p = 6) | **T1** (α-attractor T-model) | **E1** (α-attractor E-model) |
|---|---|---|---|
| V/V0 | (1 − βx − x⁶)²/f_top², x = φ/v | tanh²(φ/√(6α)) | (1 − e^{−φ/√(3α/2)})² |
| parameters | v = 0.34694 M_P, β = 1.5529×10⁻⁵ | α = 6.0968×10⁻⁵ | α = 6.0969×10⁻⁵ |
| V0^{1/4} (GeV) | 7.377×10¹⁴ | 7.376×10¹⁴ | 7.376×10¹⁴ |
| N_* (Γ_I = 0.5 H_i) | 52.79 | 52.69 | 52.69 |
| A_s | 2.100×10⁻⁹ | 2.100×10⁻⁹ | 2.100×10⁻⁹ |
| n_s | **0.9649** (tuned by β) | 0.9618 (−0.75σ) | 0.9618 (−0.75σ) |
| α_s | −9.5×10⁻⁴ | −7.3×10⁻⁴ | −7.3×10⁻⁴ |
| r | 2.71×10⁻⁷ | 2.71×10⁻⁷ | 2.71×10⁻⁷ |
| n_t | −3.5×10⁻⁸ | −3.5×10⁻⁸ | −3.5×10⁻⁸ |
| H_*/H_i, H_end/H_i | 1.0003, 0.957 | 1.0001, 0.987 | 1.0001, 0.987 |
| m_I (GeV); m_I/H_i; m_I/H_end | 5.47×10¹²; **42.4**; 44.3 | 1.65×10¹³; 128; 130 | 3.30×10¹³; 256; 259 |
| holder field φ_I = φ − φ_min at S1 / at ε_H = 1 (M_P) | −0.2549 / −0.0775 | 0.0960 / 0.0213 | 0.0894 / 0.0157 |
| s(3.8 e-folds before the end), S1 slice | 1.126 | 1.171 | 1.184 |
| g_I (μ_eff = 10i at S1 / at ε_H = 1) | 2.56×10⁻⁶ / 8.40×10⁻⁶ | 6.79×10⁻⁶ / 3.06×10⁻⁵ | 7.29×10⁻⁶ / 4.15×10⁻⁵ |
| first turning points of φ_I (M_P) | +0.0375, −0.0751 | −0.0334, +0.0278 | −0.0065, +0.0270 |
| oscillation ω (H_i), first 12 / 6 H_i⁻¹ | 29.5 → 42.3 | 49 → 127 | 86 → 255 |
| Y resonance (m_I ≲ 5 or ≳ 30 H_i) | avoided | avoided | avoided |
| φ_I → Y_LY_L open (m_I > 44.50 H_i)? | closed; 2nd harmonic open | **open** | **open** |
| Floquet μ_max/H_i at the first amplitude | 6.0 | 32 | 33 |
| linear backreaction (fragmentation) | none by 2.5 e-folds | **0.35 e-fold after the end** | **0.30 e-fold** |
| exit factor f_eff, S/Q(10⁻⁶), M1r release: holder fixed at S1 | 0.910, 4.71 (free X) | 0.867, 4.60 (free X) | 0.854, 4.56 (free X) |
| same, holder fixed at ε_H = 1 | 0.061, 2.16 (Q-matter) | 0.027, 2.04 (Q-matter) | 0.016, 2.01 (Q-matter) |
| holder-induced Y heat R_YY (M1r E_X units, Γ_I = 0.5, S1) [E] | 6.8 (Born; linear modes ×1.3) | ≈ 1.5×10³ | ≈ 5.3×10² |
| CP bookkeeping | reflection-odd; α_I = −0.949 (= 2.192 − π); μ_b = −2.58 + 13.61i in the hilltop frame | reflection-even potential; α_I = 2.192, φ_I > 0 | potential not reflection-even; α_I = 2.192, inflation only on φ > 0 |

All rows are [C] given the model choices, except the R_YY row, which is [E] (Born rates [C], ledger [E]).

**Five findings**
1. **The inflationary side is easy [C].** Every family reaches H_i = 1×10¹¹ – 1×10¹² GeV with A_s matched and
   r = 1.6×10⁻⁷ – 1.7×10⁻⁵.
   - The α-attractors give n_s = 1 − 2/N_* ≈ 0.962.
   - Pure small-field hilltops are too red: n_s = 0.943 (p = 4), 0.952 (p = 6), 0.955 (p = 8). A linear tilt β ≈ 10⁻⁵
     brings p = 6 to 0.9649.
   - The polynomial plateau V0φ²/(φ² + μ²) gives 0.972 (+1.6σ).
   - The CMB pivot sees H_* = H_i to 3×10⁻⁴. The WORKING_THEORY r values (2.7 and 3.6×10⁻⁷) are attained, not just
     floors.
2. **The branch is decided by where the holder is calibrated, not by the potential shape [C on the M1r release].**
   - On every benchmark the inflaton is a plateau field. In the last 5 e-folds of slow roll it changes by only 13–20%, so
     H/F at the horizon exit of the relic vectors stays close to its static value **if** μ_eff = 10i is imposed at the end
     of slow roll (S1, one e-fold before ε_H = 1). Then S/Q(10⁻⁶) = 4.56–4.71: free X, with e_crit = 1.80–1.90×10⁻⁶.
   - If μ_eff = 10i is imposed at ε_H = 1, the hold value 3–5 e-folds earlier is 3.3–6.7 times larger. Then f_eff is
     0.016–0.06 and S/Q(10⁻⁶) is about 2.0–2.2: Q-matter.
   - The dq_heat SKEPTIC's rolling paths (S/Q = 3.1–3.4) lie between these two calibrations.
3. **A new heat source the static engine could not see: the holder makes the inflaton decay into Y pairs
   [C rates, E heat].**
   - The trilinear −g_Iφ_IΣ*Y² couples φ_I to Y_LY_L, Y_LY_H and Y_HY_H, with thresholds 44.50, 58.35 and 72.24 H_i.
   - α-attractors have m_I ≈ π√(6A_s)M_P/N_* × (1 for T, 2 for E) = 1.65 / 3.3×10¹³ GeV [D, C]. That is 128 / 256 H_i, so
     the channel is open. Because the Y_L lives long (Γ_L = 10⁻⁷ H_i), even a branching ratio of ~10⁻¹¹ deposits
     R_YY ~ 10³ E_X. M1r's whole budget is 4.0.
   - For HT6t (m_I = 42.4 H_i) the fundamental is closed. The anharmonic second harmonic (19% of the fundamental at first)
     is open, and gives R_YY ≈ 7–9.
   - M1r's quadratic surrogate (m_I = 3 H_i, no harmonics) gives zero [D].
4. **T1 and E1 fragment within ~0.3–0.35 e-fold of the end of inflation [C linear, E criterion].**
   - μ_max/H_i ≈ 32 in the first oscillations, which agrees with Lozanov–Amin's M ≪ M_P regime.
   - The homogeneous engine therefore describes them only for ≈ 0.4 H_i⁻¹ after the end.
   - HT6t stays homogeneous.
5. **Engine.** `code/rpinf.c` is the M1r engine with an explicit-potential MODE 4.
   - Its MODE 0 output is byte-identical to rpdq.c.
   - MODE 4 reproduces our Python background to ≤ 6×10⁻⁷.
   - A smoke run from S1 (Γ_I = 0.5) gives the same-sign charge for all three benchmarks: N_X = −3781 / −5589 / −4715,
     normalised at ε_H = 1, against M1r's −8045. It also gives ε = 0.21–0.28 (M1r: 0.127) and a larger homogeneous Y_L
     remnant (5.2 for HT6t, against 0.57). This is not a release analysis.

**Recommendation for the chain track: run HT6t first, then T1 (§8).**
- HT6t is the M1r-faithful benchmark. Its n_s is on target, its oscillation is homogeneous, its m_I lies between the
  forced Y resonance and the pair threshold, and its holder is nearly static. Its clean window (Y resonance avoided and
  φ_I → Y_LY_L closed) is H_i ≈ 1.2–2.2×10¹¹ GeV, which is exactly the calibrated range.
- T1 is the reflection-even textbook attractor and the stress test.
- E1 is optional. Its inflation is identical to T1's.

## 1. Families, credits and calibration (`code/inf_common.py`, `calibrate.py`)

| key | potential | credit |
|---|---|---|
| T | V0 tanh²(φ/√(6α)) | Kallosh–Linde, arXiv:1306.5220 (T-models; r = 12α/N²) [R] |
| E | V0(1 − e^{−√(2/3α)φ})²; Starobinsky at α = 1 | Kallosh–Linde–Roest arXiv:1311.0472; Starobinsky 1980 [R] |
| HT | V0(1 − βx − x^p)²/f_top², x = φ/v, quadratic minimum at f = 0 | Boubekeur–Lyth hep-ph/0502047 (squared hilltop); Kohri–Lin–Lyth arXiv:0707.3826 (hilltop with a quadratic term). The linear tilt β is our deformation |
| PL | V0φ^q/(φ^q + μ^q), q = 2 | brane/KKLTI form (Kachru et al. hep-th/0308055) and polynomial attractors [R] |

**Method [C].**
- The background is the full equation in e-fold time, φ_NN = −(3 − ε_H)(φ_N + V′/V), integrated from ≥ 75 e-folds to
  ε_H = 1.
- The spectra are exact Mukhanov–Sasaki and tensor modes, from k/aH = 200 to 10⁻³. n_s and α_s come from a 5-point
  cubic fit in ln k.
- Check: Starobinsky at N_* = 55 gives n_s = 0.9649 and r = 3.55×10⁻³, the standard numerical values [R]. Slow roll is
  0.001 lower in n_s.
- The slow-roll estimate of A_s is 4% low at small α, so the exact spectrum is used throughout.

**Calibration.** V0 is set so that H(S1) = H_i. The shape parameter is then fixed by A_s = 2.100×10⁻⁹ at k = 0.05/Mpc
(for HT6t, (v, β) by A_s together with n_s = 0.9649). N_* is iterated to consistency:
N_* = ln(H_*/k_*) + ln(a_end/a_0).

**History for N_*.** For Γ_I = 0.5 H_i (the engine's convention), the homogeneous oscillation with a constant width is
integrated until radiation holds all but 10⁻⁶ of the energy, followed by entropy conservation with g_* = 106.75 and
T_0 = 2.7255 K. This gives T_inst = 3.0×10¹⁴ GeV.

| history (T1; HT6t) | N_* | n_s at fixed shape |
|---|---|---|
| instant | 53.00; 53.02 | 0.9620; 0.9651 |
| Γ_I = 0.1 / 0.5 / 2 H_i | 52.44 / 52.69 / 52.86; 52.56 / 52.79 / 52.93 | 0.9616–0.9619; 0.9647–0.9650 |
| w = 0, T_RH = 10¹⁴ / 10¹² / 10¹⁰ / 10⁸ GeV | 52.64 / 51.10 / 49.57 / 48.03; 52.66 / 51.12 / 49.59 / 48.05 | 0.9617 / 0.9606 / 0.9594 / 0.9581; 0.9648 / 0.9633 / 0.9617 / 0.9599 |
| w = 1/3 (any T_RH) | 53.00; 53.02 | as instant |

How to use the table:
- ΔN_* = −(1/3) ln(T_inst/T_RH) for a matter-like phase [D].
- n_s moves by ≈ 7×10⁻⁴ per e-fold. A_s at fixed shape moves too, by −6% at T_RH = 10¹² (re-calibrate when the
  reheating track supplies T_RH).
- Fragmented n = 1 attractors behave as matter (w → 0; Lozanov–Amin arXiv:1608.01213) [R], so for them use the w = 0
  rows with the actual T_RH.

**Reach in H_i [C]** (`scan_H.py`, Γ_I = 0.5):

| H_i (GeV) | 1×10¹¹ | 1.29×10¹¹ | 1.5×10¹¹ | 2×10¹¹ | 3×10¹¹ | 5×10¹¹ | 1×10¹² |
|---|---|---|---|---|---|---|---|
| T1: α | 3.64×10⁻⁵ | 6.10×10⁻⁵ | 8.27×10⁻⁵ | 1.48×10⁻⁴ | 3.36×10⁻⁴ | 9.46×10⁻⁴ | 3.87×10⁻³ |
| T1: m_I/H_i | 166 | 128 | 110 | 82 | 55 | 33 | 16 |
| HT6t: v (M_P) | 0.292 | 0.347 | 0.385 | 0.468 | 0.618 | 0.878 | 1.42 |
| HT6t: m_I/H_i | 50.4 | 42.4 | 38.2 | 31.4 | 23.8 | 16.8 | 10.5 |
| r | 1.6×10⁻⁷ | 2.7×10⁻⁷ | 3.7×10⁻⁷ | 6.5×10⁻⁷ | 1.5×10⁻⁶ | 4.1×10⁻⁶ | 1.7×10⁻⁵ |

Reading the table:
- For α-attractors m_I ≈ 1.6×10¹³ GeV (T) at every α [D: m_I = H_*/√α with H_* = πM_P(6A_sα)^{1/2}/N_*].
- The clean window (30 < m_I/H_i < 44.5) is H_i ≈ 3.7–5.5×10¹¹ GeV for T1, and ≈ 1.2–2.2×10¹¹ GeV for HT6t.
- References: HT4 has v = 0.741 and m_I = 13.3 H_i, and oscillates at 11.5–13.3 H_i, on the forced Y resonance. It is also
  5.3σ red. HT8 has m_I = 57.8 H_i (fundamental open). PL2 has μ = 7.3×10⁻⁵ M_P and m_I = 3.3×10⁴ H_i, and is not
  pursued.

## 2. The last 0–15 e-folds (needed for the vectors' horizon-exit factor) [C]

Full tables are in `results/bench_<name>.json → by_Hi[i].history`. Each row gives:
- N before ε_H = 1;
- φ and the holder field φ_I (M_P);
- H in GeV and in units of H_i;
- ε_H and η_V;
- s_S1 = φ_I/φ_I(S1) and s_end = φ_I/φ_I(ε_H = 1);
- h_r relative to S1 and to the end.

HT6t (H_i = 1.29×10¹¹):

| N before end | 0 | 0.5 | 1 (S1) | 2 | 3 | 3.8 | 5 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| φ_I (M_P) | −0.0775 | −0.2304 | −0.2549 | −0.2738 | −0.2826 | −0.2871 | −0.2918 | −0.2990 | −0.3020 |
| H (10¹¹ GeV) | 1.2340 | 1.2892 | 1.2900 | 1.2903 | 1.2903 | 1.2903 | 1.2903 | 1.2903 | 1.2903 |
| ε_H | 1 | 3.1×10⁻³ | 4.8×10⁻⁴ | 6.9×10⁻⁵ | 2.2×10⁻⁵ | 1.2×10⁻⁵ | 5.4×10⁻⁶ | 1.5×10⁻⁶ | 8.5×10⁻⁷ |
| s_S1 | 0.304 | 0.904 | 1 | 1.074 | 1.108 | 1.126 | 1.145 | 1.173 | 1.185 |

T1 (H_i = 1.29×10¹¹):

| N before end | 0 | 0.5 | 1 (S1) | 2 | 3 | 3.8 | 5 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| φ_I (M_P) | 0.0213 | 0.0863 | 0.0960 | 0.1048 | 0.1097 | 0.1124 | 0.1155 | 0.1206 | 0.1230 |
| H (10¹¹ GeV) | 1.2731 | 1.2899 | 1.2900 | 1.2901 | 1.2901 | 1.2901 | 1.2901 | 1.2901 | 1.2901 |
| ε_H | 1 | 4.2×10⁻⁴ | 8.9×10⁻⁵ | 1.9×10⁻⁵ | 7.6×10⁻⁶ | 4.5×10⁻⁶ | 2.5×10⁻⁶ | 8.9×10⁻⁷ | 5.5×10⁻⁷ |
| s_S1 | 0.222 | 0.900 | 1 | 1.092 | 1.143 | 1.171 | 1.203 | 1.257 | 1.281 |

E1 matches T1 to 1% in s_S1. The difference is a field shift.

Two points on these histories:
- **The time from S1 to ε_H = 1 is 1.0005 H_i⁻¹ (T1, E1) and 1.0033 H_i⁻¹ (HT6t).** Anything quoted from the M1r slice,
  for example the Higgs plasma deadline τ ≈ 8–10, maps to about 1 H_i⁻¹ later if S1 is the slice.
- **On every plateau the field falls by a factor 3–6 in the last e-fold.** The holder is then released during the final
  e-fold of inflation and the first oscillation, not by an oscillator starting from rest.

## 3. End of inflation and the oscillation [C]

**Turning points and frequency.** These come from `osc_gamma0` / `osc_gamma05` in the bench files. The turning-point
tables start at a = 1 at ε_H = 1 and give t in H_i⁻¹. The windowed w columns are averages over 40 windows.
- **T1.** First turning point φ_I = −0.0334 M_P, where V/V0 = 0.885: still on the plateau, about 1.75 M from the minimum.
  ω starts at 0.38 m_I = 49 H_i and rises toward m_I as the amplitude falls. ⟨w⟩ is −0.34 at first, then ≈ 0.
- **E1.** It is asymmetric:
  - steep side φ_I = −0.0065 M_P, plateau side +0.0270 M_P;
  - ω = 86 → 255 H_i;
  - its harmonics are strong (2nd 0.34, 3rd 0.16 of the fundamental) because the potential is not even.
- **HT6t.**
  - Turning points φ_I = +0.0375 M_P (steep side, V/V0 = 0.73) and −0.0751 M_P (hill side).
  - ω = 29.5 → 42.3 H_i.
  - Harmonics at the first period: 0.19 (2nd), 0.05 (3rd), 0.02 (4th).
  - ⟨w⟩ oscillates between −0.2 and +0.3 in the first windows.

**Y resonance (2m_I ≈ 22.25 H_i, i.e. m_I ≈ 11 H_i).** All three avoid it. Their lowest instantaneous frequencies are
29.5, 49 and 86 H_i. HT4 does not avoid it (ω = 11.5–13.3 H_i).

**Floquet chart, Minkowski, at the turning-point amplitudes of the real background** (`modes.py floquet`; cf. Amin et al.
arXiv:1106.3335) [C].

| | A = A₁ | 0.5 A₁ | 0.2 A₁ | 0.1 A₁ |
|---|---|---|---|---|
| T1 μ_max/m_I (μ/H_i) | 0.249 (31.8), k = 0.34 m_I | 0.208 (26.6) | 0.055 (7.0) | 0.015 (1.9) |
| E1 μ_max/m_I (μ/H_i) | 0.131 (33.5) | 0.139 (35.7) | 0.069 (17.7) | 0.027 (6.9) |
| HT6t μ_max/m_I (μ/H_i) | 0.142 (6.0) | 0.074 (3.1) | 0.019 (0.8) | 0.006 (0.24) |

**Fragmentation on the expanding background** (`modes.py resonance`) [C linear, E criterion]. Linear inflaton modes are
started from vacuum at ε_H = 1, with the width acting on the modes as on the background. The criterion is the time when
the excess variance reaches 0.1–1 of A(t)².
- **T1:** at 0.35 e-fold (τ = 0.43 H_i⁻¹), at k/a ≈ 0.35 m_I.
- **E1:** at 0.30 e-fold (τ = 0.35 H_i⁻¹).
- **HT6t:** ratio ≤ 3×10⁻⁶ after 2.5 e-folds with Γ = 0 (n_k ≤ 10⁶, narrow and fading). It decays away with Γ = 0.5.

Consequences:
- T1 and E1 fragment immediately. With a quadratic minimum, the result is oscillons or transients, then w → 0
  (Lozanov–Amin arXiv:1608.01213, 1710.06851; Krajewski–Turzyński–Wieczorek arXiv:1801.01786) [R]. Confirming that needs
  a lattice, which is not run here.
- HT6t keeps a homogeneous condensate. Antusch–Nolde–Orani (arXiv:1503.06075) found hill crossing only for much smaller
  hilltop scales [R].

## 4. The holder on a real inflaton: calibration conventions and the exit factor (`code/exit_factor.py`) [C on the M1r release]

M1r fixes g_I by μ_eff = μ + g_Ie^{iα_I}φ̂_I = 10i H_i "at the slice". With a real potential this is a choice of where:

| convention | μ_eff = 10i at | geometric slice (a = 1, H = H_i) | meaning |
|---|---|---|---|
| **S1** (primary) | 1 e-fold before ε_H = 1, where \|η_V\| ≈ 1: the end of slow roll | S1 | the static preparation is accurate there (ε_H = 5×10⁻⁴ / 9×10⁻⁵); the hold over the earlier e-folds is within 13–20% |
| Send | ε_H = 1 | ε_H = 1 | literal "end of inflation". The field is falling at φ̇ = √2 H M_P, so the static preparation is not adiabatic there |
| Cexit | 3.6 e-folds before the end (exit of the dominant relic vector modes) | S1 | the static value is reproduced at the relic's exit |

**Method.** The dq_heat SKEPTIC's method (sk_exit_path.py) is applied to our paths:
f(N) = [h_r F(1, 1)/F(s, h_r)]², with F re-solved at u = μ + s(10i − μ) and R = 6(2 − ε)h_r². The solver is
dq_model_copy.py, imported read-only. Two further inputs:
- The M1r per-mode spectrum is shifted by k → k(K/K_M1r)^{1/2} for the benchmark's K = a²H at the geometric slice. K is
  1.85 (T1, E1) and 1.52 (HT6t) at ε_H = 1, with Γ = 0.5.
- The release (N_X, E_X, mediators) is M1r's.

| e = 10⁻⁶ | HT6t | T1 | E1 | HT6 (pure, n_s −3σ) |
|---|---|---|---|---|
| **S1**: f_eff, S/Q, e_crit | 0.910, **4.71**, 1.90×10⁻⁶ | 0.867, **4.60**, 1.82×10⁻⁶ | 0.854, **4.56**, 1.80×10⁻⁶ | 0.895, 4.67, 1.87×10⁻⁶ |
| Send | 0.061, **2.16**, 1.75×10⁻⁷ | 0.027, **2.04**, < 10⁻⁷ | 0.016, **2.01**, < 10⁻⁷ | 0.062, 2.17, 1.80×10⁻⁷ |
| Cexit | 1.003, 4.96, 2.07×10⁻⁶ | 1.004, 4.96, 2.08×10⁻⁶ | 1.004, 4.96, 2.08×10⁻⁶ | 1.003, 4.96, 2.07×10⁻⁶ |

**Reading [C on the release, A on relevance].**
- In the S1 and Cexit conventions the exit states remain heavy: the lightest m²/H² at N = 3.8 is 8.8–9.9.
- Under Send they become Y-condensed, with lightest m²/H² of 200–800 and F much larger.
- The convention changes g_I by φ_I(S1)/φ_I(end) = 3.3 (HT6t) to 5.7 (E1). The two are therefore different Lagrangians,
  and with them the release changes too. The table isolates the exit factor only.
- A rolling inflaton does **not** force Q-matter at 10⁻⁶ by itself. The static result needs the hold value at the relic's
  exit to be within ≈ 20% of the calibrated one, and every plateau here supplies that if g_I is calibrated at the end of
  slow roll.

**Charge sign.**
- **T1.** The potential is even. The holder breaks φ → −φ (as in M1r), so α_I and sign(φ_I) together carry the sign, and
  α_I ≡ α_I + π holds if the decay sector is reflection-even.
- **E1.** The potential is not even, and inflation happens only on φ > 0. α_I = 2.192 is then the full CP phase; there is
  no π ambiguity.
- **HT6t.** The holder field is φ − v, with φ_I < 0 throughout inflation, so α_I = 2.192 − π = −0.949. In the field
  measured from the hilltop, the Lagrangian reads μ_bΣ*Y² + g_Ie^{iα_I}φΣ*Y² with complex μ_b = −2.58 + 13.61i H_i: the
  CP-odd input sits in a bare complex μ.
  - If the inflaton were to roll into the other vacuum (−v; the potential is not even when β ≠ 0), the late μ would differ.
    This is not a domain-wall issue within the observable patch (φ_* ≫ H/2π) [A].

## 5. The holder makes the inflaton decay into Y pairs (`code/yy_harmonics.py`, `modes.py ycontrol`)

**Couplings [D].** At the late vacuum (u = μ real, |Σ| = h/√2), −(g_Ie^{iα_I}φ_IΣ*Y² + h.c.) gives:

| coupling | value | channel threshold |
|---|---|---|
| φ_I y_H² | \|c\| = g_I\|cos α_I\|\|Σ\| | 2m_H = 72.24 H_i |
| φ_I y_L² | the same | 2m_L = 44.50 H_i |
| φ_I y_L y_H | 2g_I\|sin α_I\|\|Σ\| | m_L + m_H = 58.35 H_i |

The masses are m_L = 22.25 and m_H = 36.12 H_i.

**Rates [C].** A homogeneous oscillation φ_I = Σφ_n cos(nωt) produces pairs at the Born rates c²φ_n²β/16π (identical
pairs) and c′²φ_n²β/32π. Each full period of the computed background is Fourier analysed for n = 1–10.
- **T1 and E1 (fundamental open).** The perturbative widths are Γ_YY = 2.0×10⁻¹¹ (T1) and 1.3×10⁻¹¹ H_i (E1) with g_I
  fixed at S1.
- **HT6t (fundamental closed).** The 2nd and 3rd harmonics are open.
- **Narrow-resonance check.** q = δm²/4ω² is 0.002–0.02 (S1) and 0.01–0.07 (Send), and q²ω/H ≤ 0.2, so the
  perturbative Born treatment holds.
- **Control.** Linear light-Y modes on the HT6t background, LL channel only, 500 k up to 250 H_i, Γ = 0.5 to 12 H_i⁻¹:
  N = 1.63 against Born 1.26. **Born is low by ×1.3.**

**Heat ledger [E].** Each Y_L deposits 0.6565 E_X (dq_heat ledger). Each Y_H deposits 0.6565 + 1.78 (10⁻⁶/e)^{1/2} E_X,
through Y_H → Y_L + A_L [Q]. With Γ_I = 0.5 H_i:

| R_YY (units of M1r's E_X) | g_I fixed at S1 | at ε_H = 1 |
|---|---|---|
| HT6t | 6.8 (×1.3 linear → ≈ 9) | 74 |
| T1 | 1.5×10³ | 3.1×10⁴ |
| E1 | 5.3×10² | 1.7×10⁴ |

M1r's whole budget at e = 10⁻⁶ is R = 4.0, and its quadratic surrogate gives R_YY = 0 [D].

Without decay (Γ = 0), R_YY keeps growing: HT6t reaches 13 by 6 H_i⁻¹. **For T1 and E1 these homogeneous numbers are
indicative only**, because the inflaton fragments at 0.3–0.35 e-fold.

**Consequence [E].** On the M1r release, R_YY ≈ 7–9 alone gives W = W0(1 + 0.38 + 7–9) and S/Q ≈ 7–8 at any e. **No
Q-matter branch survives in HT6t unless the harmonic Y production is far smaller than Born and the linear modes give.**
For α-attractors S/Q reaches ~10², which means free X with a large symmetric (X̄) component. That changes M = 1 + 2p and
hence H_i; not computed here.

The design target this implies:
- m_I ≲ 5 H_i, like M1r's surrogate. No benchmark here reaches that at H_i ≈ 10¹¹ GeV, because the minimum's curvature
  is set by the inflationary field range.
- Or a nearly harmonic minimum with 30 ≲ m_I < 44.5 H_i and (φ₂/φ₁)² ≲ 10⁻³ in the first oscillations.

## 6. The modified engine (`code/rpinf.c`, `engine_test.py`, `engine_inputs.py`)

**What rpinf.c is.** A copy of `tracks/dq_entropy/code/rpdq.c` with a MODE 4.
- I is the holder field φ_I = φ − φ_min (M_r units), and v(I) = V0·shape(I + φ_min)/η² in H_i²M_r².
- The equations: I″ + (3H + Γ)I′ + v′(I) + √η² ∂V_h/∂φ̂ = 0, ρ_I = I′²/2 + v, trace −I′² + 4v, and
  Ṙ = 2(3H + Γ)I′² + 6v′I′ (the same approximation as mode 0).
- The potential comes from the environment: RPINF_KIND = 1 (T), 2 (E) or 3 (HT), and RPINF_AMP, P1, P2, P3, FTOP,
  PHIMIN.
- The F spec takes RPINF_I0 and RPINF_SIGN. |I′₀| is fixed by the H = 1 constraint.
- The argument b is used only as the step cap (pass m_I/H_i).
- The command-line interface is unchanged, so dq_model.run-style drivers work as they are.

**Tests [C].**
1. MODE 0 output is byte-identical to rpdq.c on the M1r prepared state (λ = 0.5, u = 10i, (3, 0.5), τ ≤ 20).
2. MODE 4 against the Python background, from S1 to 3 H_i⁻¹: max |ΔH/H| = 1.1×10⁻⁸ (HT6t), 7.0×10⁻⁸ (T1) and
   6.7×10⁻⁸ (E1); max |Δφ̂|/max φ̂ ≤ 6.2×10⁻⁷.

**Inputs at S1** (`results/engine_inputs.json`). η² = (H_i/M_r)² = 2.8066×10⁻¹⁵ at 1.29×10¹¹.

| | KIND | AMP = V0/η² | P1 | other | I₀ = φ_I(S1) | sign I′₀ | b (cap) | G = g_Ie^{iα_I} |
|---|---|---|---|---|---|---|---|---|
| HT6t | 3 | 3.0016355 | v = 0.34694215 | P2 = 6, P3 = 1.5529239×10⁻⁵, FTOP = 1.000000988, PHIMIN = 0.34694126 | −0.25492939 | + | 42.37 | 1.48782×10⁻⁶ − 2.07812×10⁻⁶ i |
| T1 | 1 | 3.0004369 | M = 0.019126074 | PHIMIN = 0 | 0.095970662 | − | 128.08 | −3.95214×10⁻⁶ + 5.52017×10⁻⁶ i |
| E1 | 2 | 3.0004361 | M = 0.0095631301 | PHIMIN = 0 | 0.089363035 | − | 256.16 | −4.24437×10⁻⁶ + 5.92834×10⁻⁶ i |

The hidden fields start from the static preparation at u = 10i, R = 12, as in M1r. At S1, ε_H ≤ 5×10⁻⁴, so R = 12 is
accurate to 3×10⁻⁴.

**Smoke test** (Γ_I = 0.5, τ = 200–300 from S1; rescaled to a = 1 at ε_H = 1). This is not a converged release: the
window and the ledger have not been checked.

| | N_X | E_X | ε | homogeneous N_YL |
|---|---|---|---|---|
| HT6t | −3781 | 3.64×10⁴ | −0.206 | 5.2 |
| T1 | −5589 | 4.09×10⁴ | −0.279 | 11.6 |
| E1 | −4715 | 3.99×10⁴ | −0.239 | 0.47 |
| M1r | −8045 | 1.89×10⁵ | −0.127 | 0.57 |

A rough re-ledger on these numbers (`summary.py → rough_reledger_E`) [E], for HT6t:
- S/Q(10⁻⁶) = 5.5 with the vectors only;
- 10.8 with the homogeneous Y_L added;
- 19 with R_YY added;
- 17 at e → ∞.

The direction is free X with margin. **This is for the chain track to compute properly.**

## 7. Limits

- Single-field, homogeneous background with linear fluctuations. Metric perturbations are dropped in the post-inflation
  modes.
- Fragmentation is a linear-theory criterion. No lattice was run, so T1/E1 oscillon formation is [E] from the Floquet
  rates and the literature.
- The Y-pair estimate treats Σ and u as at the late vacuum, and uses vacuum Y. The homogeneous Y/Σ condensates in the
  first H_i⁻¹ are ignored. The ledger is dq_heat's (decay at rest, full deposit).
- The exit factor and the S/Q values use M1r's release and per-mode spectrum, apart from the K shift. They isolate one
  effect.
- The holder's induced portals, inflaton sequestering and boundedness are as in fx_holder.
  - With a plateau the inflaton has no large-field quartic. The joint potential is bounded only for |φ_I| ≲ M_P, where
    the holder-induced −g_I⁴φ_I⁴/λ³ is small against V0 [A]. That holds over the whole inflationary range here
    (|φ_I| ≤ 0.32 M_P).
- No reheating channel was built (that is the reheating track's job). N_* is given for stated histories and must be
  re-iterated with the real T_RH.
- Γ_I = 0.5 H_i is the engine's supplied value. For these masses it is Γ/m_I = 0.012 (HT6t) and 0.004 (T1).

## 8. Recommendation for the chain track

1. **HT6t, the primary benchmark: H_i = 1.29×10¹¹ (free-X calibration) and 1.50×10¹¹ (Q-matter calibration), holder
   fixed at S1.** Run it end to end:
   - rpinf MODE 4 from S1;
   - Proca modes on the recorded background;
   - finite-k Y with the new inflaton (an rpmodes.c copy with v(I)), to settle the harmonic R_YY ≈ 7–9;
   - Hartree S/Q, and Higgs.

   Bracket the result with the Send calibration (g_I × 3.3).
2. **T1, the stress test.** Use H_i = 1.29×10¹¹ with the holder at S1.
   - It fragments at 0.35 e-fold, so the homogeneous release is valid only for τ ≲ 0.4 after the end.
   - The open φ_I → Y_LY_L must be added (R_YY ~ 10³).
   - Expect free X with a large symmetric component. If the chain wants a reflection-even benchmark with a clean Y
     window, T1 at H_i ≈ 4–5×10¹¹ GeV (m_I = 30–44 H_i) is the alternative, but it is outside the dark-matter
     calibration.
3. **E1, optional.** It has the same inflation as T1 and a reflection-odd potential (no α_I + π ambiguity), with
   m_I = 256 H_i.

## 9. Cross-check with the sibling reheating track (`tracks/inf_reheat`, read at the end; read-only)

inf_reheat appeared while this track was being finished. It used stand-in potentials, not these benchmarks, and its
slice is H_i := H_end, which is our Send convention.

**Where the two tracks agree.**
- Its holder-leakage rate above 2m_YL, 2–3×10⁻¹⁰ H_i (its §0 item 9), matches ours with g_I fixed at ε_H = 1:
  Γ_YY = 4.1×10⁻¹⁰ H_i for T1 [Q / C].
- Its R-spurion statement (G ↔ −G) agrees with our CP rows in §4.

**Where they differ: the hidden-heat ledger.**
- inf_reheat writes R_inf = 3M_P²K²·BR/E_X. That treats the deposit as prompt radiation.
- We use dq_heat's ledger: each Y_L is long-lived (Γ_L = 10⁻⁷ H_i), and its conformal energy grows as a until it decays
  at ⟨a_L⟩ ≈ 5.6×10³. The deposit per quantum is then m_L⟨a_L⟩ against roughly (m_I/2)·a_prod, about 10³ times larger.
- For T1 (Send, Γ_I = 0.5) this gives R_YY ≈ 9×10³ by our ledger, against O(10–60) by theirs.

**Where we add a channel.** inf_reheat did not consider anharmonic harmonics: its stand-ins are quadratic or a T-model.
The HT6t channel in §5 (R_YY ≈ 7–9 with the fundamental closed) is new here.

**T_RH.** Its T_RH ≈ 3×10¹² GeV, with a matter era of ≈ 9 e-folds, belongs to its large-α stand-ins. If that T_RH were
adopted here (w = 0), N_* would be ≈ 51.5 and n_s(T1) ≈ 0.961. HT6t would need β re-tuned (n_s at fixed shape ≈ 0.964).

A combiner should settle which ledger applies [A: ours is the one the S/Q chain has used throughout].

## Code (`code/`, all own)

| file | purpose |
|---|---|
| inf_common.py | potentials, background, Mukhanov–Sasaki spectra, post-inflation |
| calibrate.py | A_s / H_i / n_s calibration, N_* |
| benchmarks.py | benchmark records and tables → `results/bench_*.json` |
| scan_H.py | reach in H_i → `scan_H.json` |
| exit_factor.py | exit factor; read-only import of dq_heat/skeptic/code/dq_model_copy.py and dq_heat results → `exit_factor.json` |
| yy_harmonics.py | Born Y-pair production → `yy_harmonics.json` |
| modes.py | Floquet, self-resonance and Y-mode control → `modes_*.json` |
| rpinf.c | copy of dq_entropy/code/rpdq.c plus MODE 4 |
| engine_test.py | regression, MODE 4 check and smoke runs → `engine_test.json` |
| engine_inputs.py | → `engine_inputs.json` |
| summary.py | → `summary.json` |

## Reproduce (one core; about 25 min)

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/inf_potential/code; export OMP_NUM_THREADS=1
python3 -B benchmarks.py HT6t T1 && python3 -B benchmarks.py E1 HT6 && python3 -B benchmarks.py HT4 HT8 PL2   # ~10 min
python3 -B scan_H.py                                                   # 5.5 min
python3 -B exit_factor.py HT6t T1 && python3 -B exit_factor.py E1 HT6  # ~6 min
python3 -B yy_harmonics.py HT6t T1 E1                                  # ~3 min
for b in HT6t T1 E1; do python3 -B modes.py $b floquet resonance; done; python3 -B modes.py HT6t ycontrol
python3 -B engine_test.py && python3 -B engine_inputs.py && python3 -B summary.py   # binaries built in $SP/wt2/inf_potential/build
```

## Credits

Kallosh–Linde arXiv:1306.5220; Kallosh–Linde–Roest arXiv:1311.0472; Boubekeur–Lyth hep-ph/0502047; Kohri–Lin–Lyth
arXiv:0707.3826; Lozanov–Amin arXiv:1608.01213 and 1710.06851; Amin–Easther–Finkel–Flauger–Hertzberg arXiv:1106.3335;
Krajewski–Turzyński–Wieczorek arXiv:1801.01786; Antusch–Nolde–Orani arXiv:1503.06075; Planck 2018 X arXiv:1807.06211;
BICEP/Keck arXiv:2110.00483; N_* matching after Liddle–Leach astro-ph/0305263 [R]. Our own tracks reused, read-only:
- dq_entropy (rpdq.c, the S/Q closed form);
- dq_heat and its skeptic (the per-mode vector spectrum, the exit method, the ledger, dq_model_copy);
- fx_holder (the holder conventions).
