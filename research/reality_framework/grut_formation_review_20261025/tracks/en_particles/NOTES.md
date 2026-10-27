# en_particles — export, particle-count references and a first measured-dispersion count

27 October 2026. Claude review track for Codex 7d9c93b (`grut_early_nonlinear_20261027`), issue #13.
Scope: the reduced scalar-limit X/φ pilot only. Nothing here is a physical cosmological abundance, a validation or an
exclusion. Everything is classical field data from one reduced model.

**Independence labels.**
- *Own code*: everything in `code/` (`export_counts.py`, `snapshot_counts.py`, `spectral_window.py`, `response.py`,
  `subwindow_totals.py`, `zero_mode_orbit.py`, `summarize_particles.py`) and all analysis numbers.
- *Reproduction-type*: the field histories. They come from Codex's unmodified `early3d.cpp` and `initial.py`, run in a
  copy of the packet (`$SP/review/en_particles/pkt`, 2 threads). `window_run.py` and `make_initial.py` are thin
  drivers around them.
- The regenerated N48 history reproduces Codex's published τ = 800 state. ρ agrees to 1×10⁻¹⁵ relative, and the
  export rebuilt with own binning agrees with SPECTRAL_EXPORT.npz to 1.8×10⁻¹¹.

## 1. Export and reference counts [C]/[D]

**Charge and energy.** These come from SPECTRAL_EXPORT.npz alone (`results/export_counts.json`):
- Q₀ = 8556.674 and Q_nonzero = −443.197, so Q_total/nQ − 1 = 3.6×10⁻¹⁵.
- Q₀/nQ = 1.054625.
- K = 70564.77 (of which 32937.44 is zero-mode kinetic), G = 8075.02 and U = 26827.65, so E/ρ − 1 = 2.4×10⁻¹⁵.
- Hartree m² = 2λ⟨|χ|²⟩ + g²⟨φ²⟩ = 19.2636, which reproduces Codex's value.

**Counts.**
- Bin-level from the export (50 bins): massless 2.243 with ω_b² = Σk²|χ|²/Σ|χ|², or 2.171 with bin-mean k². Hartree
  1.1797.
- Per mode, exact, with own ψ = χ₁ + iχ₂ formulation on the regenerated τ = 800 snapshot (`results/snapshot_counts.json`):
  **massless 2.237515, Hartree 1.179962**. Codex has 2.2375154 and 1.1799624.

**Derivation [D].** Take a_j = (ωχ_j + iπ_j)/√(2ω) and a_± = (a₁ ∓ ia₂)/√2. Summed over ±k, the terms that are odd in k
cancel. That leaves:
- n₊ + n₋ = Σ(|π₁|² + |π₂|² + ω²(|χ₁|² + |χ₂|²))/(2ω);
- n₊ − n₋ = Re(χ₁*π₂ − χ₂*π₁), which is the code's charge density χ₁π₂ − χ₂π₁.

So Codex's ½[… ± …] is right in sign and in its factors of ½.
- With the forward FFT (field = Σ_k c_k e^{ikx}), each sum is a number density per code volume. That is equivalent to
  the 1/√V convention.
- The rfft half-space weights (2 in the interior, 1 on the kz = 0 and Nyquist planes) are correct.
- Check: at τ = 20 the seeded band (k ≈ 3.5–7.5) shows 0.51–0.76 per mode, which is the half quantum after the τ 0–20
  transfer.

**Zero point.**
- Seeded half-quanta: ½(8³ − 1.8³)/(6π²) = 4.27 per code volume per species, which is **5.3×10⁻⁴ per net charge**.
- Half a quantum in every lattice mode would be 3.9×10⁻³.
- At τ = 20 the nonzero-mode n₋/Q is 0.0019 (massless) or 0.0024 (Hartree).
- At τ = 800, modes with occupation below 1 contribute 9.3×10⁻⁴. The UV modes (k ≳ 9.5, occupation below ½) are outside
  classical validity (Epelbaum–Gelis–Wu, arXiv:1402.0115), but they are negligible in the count.
- So the zero point is at most 0.4% of either reference, and it does not cause the 2.24 versus 1.18 difference. The
  difference comes from the choice of ω at low k, where the occupation per mode is 10⁴–10⁵.

**Lower bounds.**
- The best single ω per bin gives 1.119 per net charge.
- A per-mode Cauchy–Schwarz floor gives 0.891.
- The symplectic (basis-free) charged block is strongly U(1)-anisotropic (0.63), so its eigenvalues do not split into
  n±.

## 2. Measured dispersion [E]

**Conformal time matters.** Since dτ/ds = a ≈ 52–56, the range τ 700–800 spans only **Δs = 1.83**, and the whole run
spans s = 23.9. Frequencies must therefore come from windows in s.

**How the windows were made.** The engine's clock only labels outputs, because the Hamiltonian is autonomous in s.
Restarting a saved 64-bit state with (τ₀ = s, a₀ = 1, C = 0) gives snapshots uniform in s (Δs = 0.05, h = 0.004). This
restart reproduces Codex's τ = 800 state to 9×10⁻⁷ from τ 700, or 4×10⁻⁶ from τ 340.

**Windows.**
- First: s 0–6 (τ 20–115).
- Early: s 5.56–11.56 (τ 100–255).
- Late: s 14.03–23.88 (τ 340–798).
- Extension: s 23.92–33.92 (τ 800–1460). This runs the same reduced Hamiltonian past Codex's range.
- Controls: seed 2 at N48, and L = 16 at N64, over τ 340–800.

**Estimator.** Use X(ω) = Σ_j h_j ψ_k(s_j) e^{−iωs_j} with a Hann taper h, and ω > 0 for positive charge. Then
n_± = Σ_{ω≷0}|ω|P_k(ω). Its charge, Σ ωP, matches the equal-time charge to better than 10⁻⁸ in every window.

**Retarded response.** Two twin runs kicked by ±ε in π_ψ start at τ 500, with T = 5. The results agree for ε = 0.05 and
0.2. They give G = −½{ψ, ψ*} and ρ = −2 Im G̃:
- Z₊ = 0.500–0.512 and Z₋ = 0.489–0.504 in every shell.
- The sum rule Z₊ + Z₋ = 1.001–1.008.
- So the normalisation is free-like, and the response-normalised count is n̄₋ = n_count/(2Z₋).

**Dispersion (late window).**
- The ρ peaks follow ω² ≈ k² + m² with **m² = 18.9–19.4 for k = 4–7**. The window-averaged Hartree value is 20.1.
- F peaks give m² = 17.4 ± 3. Widths (FWHM) are 1.4–1.9 at k ≥ 4, against a Hann resolution of 0.91.
- At k ≲ 2, F peaks at ±ω₀ = ±3.49, the zero-mode frequency, while ρ peaks at ±(4.0–4.5). So the coherent, non-FDT
  condensate modulation sits there.
- Branches are **symmetric about ω = 0** to within 0.17 (ρ) and 0.13 (F). They are not symmetric about ω₀, as a Bogoliubov
  spectrum of a circular condensate would be (Kapusta-type, [R]). So the co-rotating frame is not a natural basis here.

**Counts.** n₋ per net charge, nonzero modes:

| Window | Spectral | Response-normalised | Hartree equal-time (mean, range) | Massless equal-time (mean, range) |
|---|---:|---:|---|---|
| τ 100–255 | 0.700 | – | 0.684 (0.55–0.93) | 0.861 (0.43–1.52) |
| τ 340–798 | **1.145** | – | 1.152 (0.91–1.25) | 1.778 (0.91–2.77) |
| τ 547–798 | **1.158** | – | 1.166 | 1.860 |
| τ 500–750 | 1.166 | **1.177** | 1.144 | 1.815 |
| seed 2, τ 544–798 | 1.167 | – | 1.125 | 1.763 |
| L16 N64, τ 549–801 | 1.231 | – | 1.199 | 1.979 |
| τ 800–1460 (extension) | 1.273 | – | 1.298 | 2.414 |

**The measured dispersion selects the Hartree-like reference (to within 3–6%), not the massless one.** The massless
equal-time count is a phase sample:
- across Codex's own τ 700–800 snapshots it is 2.17, 1.31, 2.64, 1.40, 2.16 and 2.24;
- Hartree over the same snapshots is 1.15, 1.22, 1.08, 1.20, 1.17 and 1.18.

**The zero mode.** The spatial mean is an eccentric orbit, with lines at ±ω₀ holding 99.7% of its power.
- Its counter-rotating part is 1.05 (τ 100–255), 0.61 (τ 340–798) and 0.48 (τ 800–1460) per net charge.
- **Zero-mode plus nonzero-mode n₋ is 1.741–1.753 in all 22 sub-windows from τ 20 to 1460, and in both controls.**
- The homogeneous τ = 20 orbit alone gives 1.739.

So in this lab-frame quasiparticle basis, the X̄ number is set by the eccentricity of the matched orbit. It is
approximately conserved, as expected when number-conserving 2↔2 scattering dominates (Micha–Tkachev; Berges–Sexty
[R]). Fragmentation moves it from the coherent orbit into waves, with an inverse cascade toward k ≈ 1.

**φ quanta per net charge.** 0.0030 (spectral), against 0.0030–0.0032 (Hartree) and 0.0036 (massless).

## 3. Comparison with REVIEW_TIMING [E]

- **Gas at ω, 1.05–2.48 for (3,.5).** The nonzero-mode values 1.15–1.27 sit at its low end.
- **Excess energy into pairs, 1.36–2.62.** The nonzero-mode count by τ 800 is below this range. The total including the
  zero mode, 1.75, is inside it. Energy outside the mean is 63%, against the ≤ 88% budget.
- **Fully equilibrated massless end point, 0.10–0.12.** It is not approached by τ 1460, because the count is conserved.
  Reaching it requires number-changing processes that are absent on these time scales.
- **kin_source bracket, 0.13–2.6.** It contains every measured value.

## Verdicts

| # | Claim | Codex | Ours | Verdict |
|---|---|---|---|---|
| 1 | Massless reference n₋/net at τ 800, N48 [C] | 2.238 | 2.237515 per mode; 2.243 / 2.171 from bins | CONFIRMED |
| 2 | Isotropic Hartree reference [C] | 1.180 | 1.179962; 1.1797 from bins | CONFIRMED |
| 3 | Export reconstructs charge and energy [C] | 4.5e-16 / 2.3e-16 | 3.6e-15 / 2.4e-15 (export only); rebuild 1.8e-11 | CONFIRMED |
| 4 | n± formula: sign, ½ factors, normalisation [D] | as written | derived; vacuum check 0.5 per mode | CONFIRMED |
| 5 | The ~1.9× split exceeds grid/time-step changes [C] | yes | yes at fixed τ, but the massless value swings 1.31–2.64 between τ 700–800 snapshots | CORRECTED (phase-dependent) |
| 6 | Mean charge / net charge at τ 800 [C] | 1.0546 | 1.054625; oscillates 1.047–1.169 over τ 340–798 | CONFIRMED (phase sample) |
| 7 | First measured-dispersion count [E] | open | 1.15–1.18 (N48), 1.17 (seed 2), 1.23 (L16); Hartree-like | NEW |
| 8 | Zero-mode X̄ content and conserved total [E] | not counted | 0.59–0.61 now; total 1.74–1.75 | NEW |
| 9 | Seeded zero point per net charge [C] | not given | 5.3e-4 (band), ≤ 3.9e-3 (lattice) | NEW |
| 10 | Our own "pairs near k_peak" 1.36–2.62 [E] | – | nonzero 1.15–1.23 by τ 800; total 1.75 | CORRECTED (nonzero part below) |

## Limits (not certified)

- Classical-statistical runs with single realisations: two seeds, two boxes, N48/N64.
- At k ≲ 2.6, which holds about 90% of n₋, F and ρ peak at different ω. The classical FDT therefore fails there, and the
  particle reading is weakest where it matters most.
- The ± split of the coherent zero mode is a choice of basis.
- The response was measured at one s₀, with a π kick only.
- The extension window runs the reduced Hamiltonian past τ 800.
- Not included: gauge, heavy-field and quantum effects, collision terms, and the binding transition.
- The physical n̄, f_sym and threshold remain open.

## Reproduce

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad; R=$SP/review/en_particles; C=tracks/en_particles/code
export OMP_NUM_THREADS=2 GRUT_EARLY_WORKDIR=$R/ework   # plus the env of BRIEF_EARLY.md
python3 -I $C/export_counts.py $P/SPECTRAL_EXPORT.npz $P/FIELD_RESULTS.json out.json
(cd $R/pkt/grut_early_nonlinear_20261027/code && python3 -B -c "from run import run; run('N48',N=48,threads=2)")
python3 -B $C/window_run.py $R/pkt/grut_early_nonlinear_20261027/code $R/ework/N48/state_s0000016.bin $R/W_late 9.85 0.05 2.0
python3 -I $C/spectral_window.py $R/W_late/window.npz out.json
python3 -B $C/window_run.py ... state_s0000024.bin $R/K_0.05 5.0 0.05 2.5 2 0.004 0.05 7   # and -0.05
python3 -I $C/response.py $R/K_0.05/window.npz $R/K_-0.05/window.npz $R/K_0.05/kick.npy out.json stat.json
SUB=100 STEP=20 python3 -I $C/subwindow_totals.py out.json $R/W_*/window.npz
```

Credits: Khlebnikov–Tkachev (classical lattice method); Micha–Tkachev (turbulent kinetics, particle numbers); Aarts–Berges
(classical spectral functions) [R]; Berges–Sexty (number-conserving condensation) [R]; Epelbaum–Gelis–Wu arXiv:1402.0115
(classical UV limits).
