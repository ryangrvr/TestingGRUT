# Brief: review of Codex 5f5513b (averaged classical source, infrared controls, coherent matching), issue #13

27 October 2026. The Codex packet is `research/reality_framework/grut_particle_transport_20261027/` at commit
5f5513b; read it through git only. It accepts our review fd71d3a (REVIEW_EARLY.md) and continues the same reduced
conformal X/φ model:
- H = ½π² + ½(∇χ)² + (λ_X/4)r⁴ + (g²/2)χ_φ²r² + (λ_φ/4)χ_φ⁴, with λ_X = 10⁻³, g² = 0.01, λ_φ = 0.03;
- FD4 lattice with Yoshida stepping;
- clock τ = 20 + a_m s + Cs²/2, with a_m = 8.869218738806 and C = 1.9857962750634526.

## Codex's claims

1. **Eight histories to s = 34 (τ ≈ 1469) give a classical negative-frequency action budget A₋.**
   - The budget is 1.7505–1.7517 per net charge over s 14–24 and 1.7413–1.7454 over s 24–34.
   - The histories cover: base48 and its half step; two infrared seeds each at L = 12 (N32) and L = 24 (N64); a
     second unstable-mode width; and an L12 N48 refinement.
   - The coherent/wave split moves while the total stays near 1.74–1.75.
   - A₋ = Σ_{k, ω<0} |ω| P_k(ω), with Hann-windowed P_k(ω) = |Σ_j h_j ψ_k(s_j) e^{−iω s_j}|² / (M Σ h²) and
     ψ = χ₁ + iχ₂.
2. **Infrared preparation.** The initial six-field Hessian has eigenvalues −2.83238 and −1.99965, so there is no
   instantaneous ground state below k ≈ 1.683 and 1.414. In the Hessian eigenbasis the added modes use
   C_qq = ½ diag(1/Ω), C_pp = ½ diag(Ω), with Ω = √(k² + σ²) and σ = 1 or 3 in the unstable directions. At the τ = 20
   handoff, the X perturbation RMS is about 4.3% and the momentum RMS about 2.2%.
3. **Coherent counting.**
   - In a free massive basis, ψ = Ae^{+ims} + Be^{−ims} gives n₊ = m\|A\|² and n₋ = m\|B\|², so n₋ equals the radial
     action J_r.
   - The starting homogeneous orbit has a lab-frequency statistic of 1.737788 and J_r = 1.763725 per coherent charge.
   - Slow free matching, m(s) = 30f(s/T) with λ ramped to 0, returns J_r. A fixed-quartic mass-only ramp agrees to
     5×10⁻⁷.
   - Sudden matching gives 0.077–18.18, depending on phase.
   - So no extra counting postulate is needed.
4. **The medium response is a matrix.** Three tangent probes give a conjugate/normal RMS ratio of 0.44–1.70 in the
   count-dominant infrared. The canonical sum holds within 0.4%. A free control matches to 2.2×10⁻¹⁰ with zero mixing.
5. **Averaged export.** The coherent outer product is kept separately. Charge and energy reconstruct to 6×10⁻¹⁵ and
   1.4×10⁻¹⁵. Both uniform-conformal and uniform-cosmic averages are supplied.
6. **Next.** Chemical relaxation before binding is now the deciding calculation.

## Codex's review request (CLAUDE_HANDOFF.md)

Check:
- the initial-infrared width prescription and its canonical covariance;
- the averaged export's coherent bookkeeping;
- above all, the distinction between the radial action and free matching.

A short engine test plus one IR12 history is enough to start. There is no need to rerun the old seven.

The rules of [BRIEF.md](BRIEF.md) apply:
- write only in `tracks/<key>/` (nothing over 2 MB);
- git read-only, and no commits;
- label independent code versus reproduction;
- give a verdict per claim, with claim labels;
- never claim validation, exclusion or new theory;
- return a digest of at most about 700 words.

## Prepared environment

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
# Codex tree (git archive of 5f5513b, whole research/reality_framework): $SP/f55/research/reality_framework/
export GRUT_KINETIC_INPUT_ROOT=$SP/pinned/research/reality_framework/grut_formation_review_20261025
export GRUT_KINETIC_WORKDIR=$SP/kwork GRUT_ONSET_WORKDIR=$SP/owork
export GRUT_TIMING_WORKDIR=$SP/review/<key>/twork GRUT_PARTICLE_WORKDIR=$SP/review/<key>/pwork
export GRUT_PARTICLE_WORKERS=1 GRUT_PARTICLE_THREADS=2
cd $SP/f55/research/reality_framework/grut_particle_transport_20261027/code && python3 -B -c "import pt_common"
```

- Codex scripts write outputs into the packet directory. Only `pt_rerun` may run them in place in `$SP/f55`, after
  copying the published outputs to `$SP/review/pt_rerun/published/`.
- Every other track that needs a Codex script copies the whole `$SP/f55/research/reality_framework` tree to
  `$SP/review/<key>/rf/` and runs from the copy.
- The machine has 4 cores and another workflow is running, so use at most 2 threads.
- Our own earlier code is in `tracks/en_indep/code` (an independent FD4/SRKN integrator), `tracks/en_particles` and
  `tracks/en_critic`.
