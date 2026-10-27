# Skeptic of pt_rerun: Codex 5f5513b, issue #13

27 October 2026. Claude, key `pt_rerun_skeptic`. This is an adversarial check of [NOTES.md](NOTES.md). Brief:
[BRIEF_TRANSPORT.md](../../BRIEF_TRANSPORT.md).

This is same-session AI work on one supplied, reduced classical scalar model. It is not outside replication. Nothing
here is a validation, an exclusion, a GRUT prediction or new theory.

## Which numbers are which

Every number below is **independent**: our own code, written for this review, in `skeptic/code/`. Results are in
`skeptic/results/`. No Codex script was run.

The inputs are of three kinds:
- Codex's published, git-identical files: `PREPARATION_L12/L24.npz`, `FIELD_TRACES.npz`, `AVERAGE_*.npz`,
  `WINDOW_*.json`, `COHERENT_MATCHING.json` and `MATCH_L12.npz`;
- pt_rerun's rerun of IR12_s1_w1 (`psi.npy`, `trace.npz`);
- `RESULTS.json`, read through `git show 5f5513b:`.

Methods were chosen to differ from both Codex and pt_rerun:
- **Coherent orbit.** We use the u = r² cubic: T_r in closed form via ellipk, and J_r by 30-digit tanh-sinh quadrature.
  The lab statistic comes from our own 8th-order Kahan–Li symplectic integrator.
- **Infrared.** No transfer integrator is used. The transfer rows are recovered from Codex's factors and splined in k²,
  and the singular state factor is applied exactly.
- **Spectra.** The DFT is an explicit matrix with absolute-time phases, plus a 4× zero-padded grid.
- **Coherent part.** All eight histories are taken from the published float64 traces.

## Verdicts on pt_rerun's claims

| pt_rerun claim | pt_rerun | Skeptic (own code) | Refuted? |
|---|---|---|---|
| C1a. IR12_s1_w1 budget 1.7505527 / 1.7414468 | reproduction, 1 of 8 histories | Explicit DFT of rerun ψ: 1.7505527 / 1.7414468 (Δ 2.6e-14 / 8e-15). Coherent part of **all 8** histories × 2 windows from published traces: ≤ 4.6e-9 | No |
| C1b. Σ_{ω<0}\|ω\|P = m\|B\|² | synthetic tones; orbit ±5e-6 | 4× zero padding changes the total by 7.8e-7 / 4.0e-6 | No |
| C1c (NEW). Totals including the k > 8 tail | 1.7530–1.7546; **1.7501–1.7513**; refinement 0.048%/0.056% | 1.7530–1.7546; **1.7501–1.7518** (IR12_s2_w1: 1.741795 + 0.009964 = 1.751759 was missed). Refinement (fine − coarse) is **−0.048% / −0.056%**: the sign flips from +0.052% / +0.160% | **Yes**, in part: one bound and the sign. The substance holds (see N4) |
| C2a. Eigenvalues; "15 / 60 unstable oscillators" | −2.832376, −1.999649; 1.682967, 1.414090 | Same. The 15/60 are (shell, direction) pairs, Codex's own label. Lattice vectors: L12 146 + 80 = 226; L24 1140 + 690 = 1830 | No (label) |
| C2b. Pure canonical state and normalisation | Williamson ½ ± 3e-15 | Canonical pair algebra: c_k = (X+iY)/√(2V) gives ⟨\|c_k\|²⟩ = 1/(2ΩV). The width-3 factor rebuilt from the width-1 transfer agrees to 1.1e-15 | No |
| C2c. Continuum handoff RMS | 4.437% / 2.235% (w1); 4.390% / 2.211% (w3) | No integrator: **4.43726% / 2.23501%** and **4.39021% / 2.21098%**, agreeing to ~1e-9 relative. Spline leave-half-out error ≤ 1.3e-6 | No |
| C2c. "Insensitive to widths" (σ 0.1–30 → 4.25–4.40%) | — | True only inside Codex's σ family. σ = 300 gives 4.665%. Other pure rules on k < 1.8 give 4.325% (clamp Ω ≥ 1) and 4.664% (Ω = √\|k²+λ\|) | No (qualified) |
| C2d (NEW). "Jumps ×2.4 at 1.414 and ×5.8 at 1.683"; "(1.683, 1.8) ≈ 11%" | — | No finite jump. On the stable side, the per-mode τ = 20 X variance diverges as c(k−k₀)^−½ (c = 2.47 and 13.76). ×2.39 and ×5.82 are the ratios between the nearest **L24** shells; L12 gives ×2.35 and ×3.36. With a panel break at 1.8, the share is **11.66%** (σ = 1) and 11.92% (σ = 3); pt_rerun's base quadrature gave 11.14%, its doubled-node run 11.67% | No (wording and number refined) |
| C3a. J_r 1.763725279; lab 1.737788353 | 7e-16; 1.4e-10 | J_r: 1.7637252793969964 (Δ −6.7e-16). T_r and Δθ agree to ≤ 2e-15. Lab statistic, KL8 at M = 2048–8192: Δ ≤ 1.5e-13 | No |
| C3b. Slow matching → J_r | smoothstep ≤ 1.7e-5 (T = 80); ≤ 2.2e-6 (T = 200) | C^∞ ramp at m_out = 30, 8 phases. The deviation falls faster than any power of T: 0.14 (T = 5), 2.9e-3 (10), 9.6e-7 (20), **2.9e-10 (40)**, 2.9e-11 (80, tolerance floor). The effective radial well is single for Q > 0, so there is no separatrix (Q₀/Q ≥ 0.95 over s 14–34 in the medium as well) | No (stronger) |
| C3c. Sudden range (CORRECTED by pt_rerun to 0.013–18.30) | 64 phases at m_out = 30 | n₋ depends only on r² (energy conservation), so the exact extremes sit at the turning points: **0.00535–18.312** at m_out = 30. At m_out = 25.92 (= Q/r_min²) the minimum is **0** (3e-11). m_out = 10: 0.2445–5.7795; 100: 0.5293–62.197; 3.678 (= Δθ/T_r): 1.2973–2.0284 | No: the CORRECTED verdict stands, with exact values |
| Qualification: "the 1.5% lab-vs-J_r gap indicates the size of the late ambiguity" | [A] | **Underestimate.** See N1 below: the late gap is +2…+6% (bare quartic J_r) or −8…−10% (Hartree J_r) | **Yes** |
| C5. Export: coherent part 113–129% of Q, "28–33% of the energy" | IR12_s1_w1 only | All 8 histories × 2 windows × 2 rules: ⟨v₀v₀ᵀ⟩ from traces ≤ 1.5e-15. Coherent charge share **101–129%** (base48 s 24–34: 101%). 26–33% is the coherent **kinetic** share; its own quartic adds 7–11%. The squared mean keeps 0.35–0.49% of Q | No (range widened) |
| Design note: IR12_s1 ≡ base48 seed band | code reading | Confirmed. `initial.py` (7d9c93b) and `prepare.initialize` use the same SeedSequence([seed, 100L, i+128, j+128, k+128]) and 14-normal draws, and the old-band factors are copied | No |
| C4 (matrix response), C6 (chemical relaxation) | NOT CHECKED | not examined | — |

## New issues (N)

**N1. The late coherent component has a matching ambiguity of several per cent, not 1.5%** [C/E].

We made two radial-action variants for the zero mode of every history and window:
- bare-quartic J_r of the zero mode alone;
- J_r with the fluctuation Hartree mass M² = m²_H − 2λ|ψ₀|² ≈ 9–11.

We Hann-averaged them and compared them with the lab statistic (0.40–0.68 per net charge):
- J_q − lab = **+0.013 to +0.023**, a ratio of 1.02–1.06;
- J_H − lab = **−0.039 to −0.055**, a ratio of 0.90–0.92.

On the total of ≈ 1.745 this is **−3.2% to +1.3%**. That is 10–40× the eight-history control spread (0.07% / 0.24%).

Within each 10-unit window the instantaneous J_q swings by ×3–5 (for example 0.17–0.87). The zero mode exchanges
action with the waves and is not on an adiabatic central orbit, so a "slow-matched" coherent value is not defined
without that coupling.
- The anomalous ⟨δψ²⟩ψ₀* term is not central and is omitted.
- Neither variant is "the" answer; they bracket the definitional freedom.
- Equal-time harmonic counts are at least the frequency statistic for quasi-free modes. For coherent plus waves at
  m² = 7, 12, 20 and 26 they give 1.82, 1.76, 1.83 and 1.92 (IR12_s1_w1, s 24–34), against 1.741.

**N2. The tail-included total still declines, and the residual grows** [C]. These are 5-unit Hann windows of the rerun ψ:
- k ≤ 8: 1.7517 → 1.7388 (−0.73%);
- k ≤ 8 + tail: 1.75403, 1.75400, 1.75371, 1.75312, 1.75212, 1.75112, 1.75026 (−0.22%).

The successive drops grow from 3e-5 to about 1e-3 per 2.5 units of s. So "most of the decline is cascade into k > 8"
holds, at 71% of it. But the remainder is not zero and is not constant.

**N3. The tail is grid content.** The k > 8 tail differs by ×1.6–2.1 between N32 and N48 of the same realization:
0.00336 vs 0.00160 (s 14–24) and 0.00966 vs 0.00590 (s 24–34). On N32 the "k > 8" modes are cube corners with FD4
k_eff ≈ 6–10.7. Read the tail-included sum as "all lattice modes of that grid", not as continuum k > 8.

**N4. pt_rerun's tail arithmetic.** The correct ranges are 1.75299–1.75459 (s 14–24) and 1.75013–1.75176 (s 24–34). The
window-to-window change is −0.145% to −0.170%. base48 → IR12_fine including the tail is +0.004% / −0.019%.

**N5. The handoff depends on the state rule at the −2.5% / +5% level.** Keeping the old band fixed, pure Gaussian rules for k < 1.8 give a
continuum X RMS of 4.33–4.67% and a momentum RMS of 2.21–2.33%. The k < 1.8 share runs 15–27%. Base48 → IR12_fine
moves A₋ by only −0.01% / −0.07%, so this matters for the stated handoff, not for the budget.

## What stands

pt_rerun's reproduction and its CONFIRMED verdicts survive every independent recomputation tried here. The changes are
the ones listed above:
- one wrong bound (1.7513 → 1.7518);
- a sign;
- the sudden-switch extremes (exact, with a zero minimum attainable);
- the jump wording (a divergence, not a jump);
- the export range across all 8 histories;
- the size of the late coherent matching ambiguity (N1), which is the substantive change.

Not checked: the wave part of seven histories (no ψ files), the response probes, chemical relaxation, and any
physical n̄, f_sym or threshold.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
P=$SP/review/pt_rerun/published; R=$SP/review/pt_rerun/pwork/IR12_s1_w1
T=research/reality_framework/grut_formation_review_20261025/tracks/pt_rerun/skeptic
E=$SP/f55/research/reality_framework/grut_early_nonlinear_20261027
python3 -I $T/code/sk_coherent.py $E/MATCH_L12.npz $P/COHERENT_MATCHING.json $T/results/sk_coherent.json
python3 -I $T/code/sk_ir.py $P/PREPARATION_L24.npz $P/PREPARATION_L12.npz $T/results/sk_ir.json
python3 -I $T/code/sk_ir_alt.py $P/PREPARATION_L24.npz $T/results/sk_ir_alt.json
python3 -I $T/code/sk_traces.py $P/FIELD_TRACES.npz $P $T/results/sk_traces.json
python3 -I $T/code/sk_spectrum.py $R $P/WINDOW_IR12_s1_w1_14_24.json $P/WINDOW_IR12_s1_w1_24_34.json $P/FIELD_TRACES.npz $T/results/sk_spectrum.json
python3 -I $T/code/sk_subwindows.py $R $T/results/sk_subwindows.json
```

The tail arithmetic (N4) comes from `git show 5f5513b:…/RESULTS.json` and is recorded in `skeptic/results/sk_verdicts.json`.
