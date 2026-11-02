# Brief: a multi-field inflaton sector for M1r

Read [WORKING_THEORY.md](WORKING_THEORY.md) §0, §4a, §7–§9 and [BRIEF_INFLATON.md](BRIEF_INFLATON.md) first; the
earlier rules apply (end of this file). All new tracks are `tracks/mf_<key>/`.

## Where we are

The inflaton pass (tracks `inf_*`, combiner `inf_critic`) found that no single-field plateau or hilltop realises M1r's
assumed history. That history has three parts:
- **quiet:** the holder field oscillates with m ≈ 2–4 H_i, so it does not drive Y_L;
- **fast:** radiation domination within about a Hubble time, K = a²H ≈ 2;
- **low-scale:** H_i ≈ 1.2–1.5×10¹¹ GeV, r ≈ 3×10⁻⁷.

On every explicit benchmark the relic is free X X̄ (S/Q = 155–987), with no e_crit. The open item is a sector beyond a
single field. This pass builds and tests one.

## The external proposal (Gemini, relayed by the owner; treat it as input to test, not as fact)

**Two-field hybrid:**
- V = λ_χ/4 (χ² − v_χ²)² + ½ g_hyb² σ² χ² + V_SR(σ).
- σ is the slow-roll clock. χ is held at 0 during inflation and waterfalls at σ_c.
- λ_χ v_χ⁴/4 = 3 M_P² H_i² with H_i ∈ [1.25, 1.55]×10¹¹ GeV, and m_χ = √(2λ_χ) v_χ = 3 H_i.
- χ reheats fast through y_N χ N̄N with y_N ≈ 1.5–2, plus a Higgs portal λ_χH χ²|H|², with Γ_eff/H_i ≈ 0.5–1.
- V_SR(σ) is chosen to give n_s ≈ 0.967–0.973, to remove the ACT tension.
- The claims: S/Q < 3.47, Q-matter at Q ~ 10²⁰, Higgs rescued by the plasma, and leptogenesis kept.

**Premises in the proposal that our record contradicts (do not build on them):**
- **"Any realistic rolling history gives Q-matter."** Refuted in the inflaton pass. That result changed only the
  horizon-exit factor on M1r's own release. On the benchmarks' own releases S/Q = 155–987.
- **What M1r's supplied history gives at e = 10⁻⁶.** With (m_I, Γ_I) = (3, 0.5) H_i on the static slice it gives
  free X (S/Q = 4.955). Q-matter needs e > 2.06×10⁻⁶ or less heat. Recovering M1r's history recovers the branch
  point in e; it does not by itself give Q-matter at 10⁻⁶.
- **The width formula.** For −(y/2)φN̄N (Majorana) it is Γ = y²mβ³/16π; y²m/8π is the Dirac normalisation
  (inf_reheat).
- **The scale.** V0^{1/4} = (3 M_P² H_i²)^{1/4} = 7.4×10¹⁴ GeV at H_i = 1.3×10¹¹ (reduced M_P), not 2.4×10¹⁵.
- **The tilt.** A V0-dominated hybrid with a convex V_SR(σ) has a blue tilt. A red tilt needs negative curvature
  (hilltop) or a log or plateau term; it is not automatic.

## Claude's quick estimates (labelled [E]; verify or refute each, do not assume)

- **G1, energy and amplitude.** A field of mass m oscillating with the full energy V0 = 3M_P²H² has amplitude
  √6 H M_P/m: 0.82 M_P at m = 3H. A quiet field that carries the energy therefore has a Planckian amplitude.
  - Fast perturbative decay is then blocked: y|φ| ≫ m/2.
  - The one-loop potential from an O(1) Yukawa swamps the tree potential.
  - This is the single-field obstruction again. It suggests the holder must sit on a **subdominant** field.
- **G2, the proposal as written.**
  - m_χ = 3H and λ_χ v_χ⁴/4 = V0 give v_χ = √(8/3) M_P = 1.63 M_P and λ_χ = 4.8×10⁻¹⁵.
  - With y = 2, y⁴/16π² = 0.10 ≫ λ_χ: the Coleman–Weinberg term destroys χ's potential, by 2×10¹³.
  - y v_χ/(m_χ/2) = 4×10⁷: kinematically blocked.
  - The tachyonic mass² at χ = 0 is −m_χ²/2 = −4.5 H², so the waterfall is "mild". It lasts about 17 e-folds from
    H/2π to v_χ. That raises the CMB and small-scale curvature (PBH) questions.
  - A real χ gives domain walls with v ~ M_P.
  - Holder placement: χ = 0 during inflation, so a holder on χ is off during inflation and cannot prepare the state.
- **G3, a heavy waterfall can reheat perturbatively.** If N gets its mass from the χ vev (M_N = y v_χ), the decay
  needs y < √(λ_χ/2), so Γ_χ ≤ λ_χ m_χ/32π.
  - With V0 fixed, Γ_χ ≥ H needs m_χ ≳ 3–5×10³ H_i, with v_χ ≈ 1–1.6×10⁻³ M_P and λ_χ ≈ 5×10⁻³–4×10⁻².
  - Then y⁴/16π² ≪ λ_χ, so the loop is safe.
  - N with M_N ~ 10³ H_i decays to LH within the Higgs deadline (seesaw y_ν ~ 0.5) [E].
- **G4, the clock σ after the waterfall.** On the valley χ² = v²(1 − σ²/σ_c²), V = ½m_σ²σ² − m_σ²σ⁴/(4σ_c²), with
  m_σ = g_hyb v_χ and V(σ_c) = V0. So σ, not χ, ends up holding about V0 and must decay fast too.
  - σ → −σ makes σ stable unless the symmetry is broken.
  - Hilltop example (η = −0.0175, r = 2.75×10⁻⁷): σ_* = 0.0106 M_P, σ_c ≈ 0.028 M_P, m_σ ≈ 125 H_i.
  - A Yukawa on σ is capped near 2×10⁻³ by flatness, giving Γ_σ ~ 10⁻⁵ H.
  - The σ decay route (odd terms in V_SR, σ → χχ, σσ → χχ, or tachyonic preheating plus χ decays) is a central
    question.
- **G5, a spectator holder.** Put the holder on a field s:
  - the holder term is −(g_s e^{iα} s Σ*Y² + h.c.);
  - s is light during inflation, so it is frozen at s_* and the preparation is stationary;
  - s gains m_s ≈ 3 H_i only after inflation.
  - Two ways to get that mass:
    - (a) from the waterfall: ½ g_s² χ² s², so m_s = g_s v_χ;
    - (b) from a non-minimal coupling −½ξ_s R s² that cancels a bare m_s² = 9 H_i² during inflation, with ξ_s ≈ −3/4
      for R = 12H². This needs H nearly constant, which a V0-dominated hybrid gives.
  - The CP sign is α plus sign(s_*). The reflection s → −s.
  - **Isocurvature.** In the asymmetric branch, Ω_DM depends on s_*. Planck's β_iso bound (arXiv:1807.06211) roughly
    requires s_* ≳ 2×10¹⁵ GeV × |d ln Ω_DM/d ln s_*|.
  - **Energy.** s's energy fraction at reheating is ≈ 1.5 (s_*/M_P)². It is matter-like until s decays, so check
    early s domination, entropy dilution, and any s decay into the hidden sector (which would add hidden heat).
- **G6, Higgs.** A portal linear in a fluctuation about the vacuum gives a Higgs mass² that oscillates in sign. Examples
  are λ_χH(χ² − v²)|H|² ≈ 2λ_χH v δχ |H|², or Aσ|H|². These are inf_reheat's tachyonic windows. A positive χ²|H|² or
  (v² − χ²)|H|² term can instead keep the Higgs heavy during inflation (cf. Lebedev–Westphal arXiv:1210.6987). Check
  the sign structure through the waterfall.
- **G7, heat.** R_vec ∝ K² and K ≥ 1 for radiation from the slice; M1r's (3, 0.5) gives 1.99. Fast reheating with
  K → 1–1.3 could cut the vector heat by up to ×4. But the release itself changes when the holder switches off by
  mass turn-on and redshift instead of decay. **Do not prejudge the branch; compute it.**

## Tracks

### mf_gemini: audit the proposal as written

- Verify or refute G1 and G2 with your own code: the background and the waterfall, both linear and δN.
- Find whether **any** choice of (m_χ, y_N, λ_χH, g_hyb, V_SR, holder placement) meets the proposal's own stated
  requirements: m_χ = 3H_i, Γ_eff/H_i ≥ 0.5, H_i ∈ [1.25, 1.55]×10¹¹, n_s 0.967–0.973, no defects, holder working.
  Include tachyonic preheating as a non-perturbative route.
- Check the claimed n_s range for the hybrid family.
- If the proposal fails, state the minimal change that would make it work, if any.

### mf_hybrid: build the repaired architecture and export benchmark histories

**What to build.** A heavy, small-field hybrid with a spectator holder (G3–G5), in at least two benchmarks, and a
variant with the ξ_s spectator if it works.

**Requirements:**
- **CMB.** Planck 2018 n_s = 0.9649 ± 0.0042, A_s = 2.10×10⁻⁹; r < 0.036 (BICEP/Keck 2021, arXiv:2110.00483).
  Report the ACT DR6 pull (n_s = 0.9743 ± 0.0034, arXiv:2503.14452) retrospectively. Give α_s.
- **Waterfall.** Fast, with sub-Hubble tachyonic growth, or show that a mild waterfall is safe. The curvature power
  on small scales: PBH.
- **Defects.** Avoid them: a tilted or "smooth" hybrid (Lazarides–Panagiotakopoulos), a biased waterfall, or show that
  they are inflated away.
- **Radiative stability.** The full one-loop potential along the trajectory.
- **Decay channels.** For σ, for χ and for s, kinematically open, with the CP-sign bookkeeping kept.
- **The post-inflation history.** Equation of state, K, T_RH, and when a Higgs-coupled plasma appears.
- **The spectator s.** Its value s_*, its isocurvature, its energy fraction, and its decay.
- **Gravitational effects on the hidden sector** at the abrupt end:
  - R oscillating at 2m_χ and 2m_σ;
  - ξ_X = −5/6 production;
  - the Y resonance 2m ≈ ω_YL = 22.25 H_i.

**Lattice.** A 3D lattice of the waterfall and preheating is allowed here, at 64³ or below.

**Exports for the chain.** For each benchmark write `tracks/mf_hybrid/results/bench_<name>.json` (< 2 MB) holding the
history from at least 5 e-folds before the slice to K convergence:
- t (in H_i⁻¹), N, a, H, ε_H, R;
- the holder field s(t), or m_s(t) and s_*;
- the energy fractions;
- every model parameter.

Note H_i is the Hubble rate at the M1r slice. Define the slice as M1r does: the end-of-inflation slice, with a = 1 at
H = H_i. State the convention.

**Interim export.** Write an early `bench_*.json` as soon as one benchmark history exists, and mark it interim in
`results/STATUS.md`.

### mf_ideal: the chain on idealised fast histories (starts at once, independent of mf_hybrid)

**The engine.** Copy the engine (`tracks/inf_chain/code/rpchain.c`, or `tracks/inf_potential/code/rpinf.c` MODE 4)
and add a mode with two features:
- a tabulated or parametrised background: radiation after the slice with a short matter-like interval chosen to give
  K = 1.0, 1.3, 1.6, 2.0 and 3.0;
- the holder on a spectator s. The holder is μ_eff = μ + G s(t)/H_i, with s frozen during inflation and the mass
  m_s ∈ {2, 3, 4} H_i turned on at the slice over a time t_on ∈ {0.1, 1} H_i⁻¹. s then evolves by its own equation of
  motion (no decay, or a small decay).

**Regression first.** Reproduce M1r's (3, 0.5) release in the copy (dq_entropy / inf_chain benchmark).

**The full chain** on each idealised history:
- release, with N_X, ε, E/E_c and the Y_L / Y_H remnants;
- vectors: the `tracks/dq_heat` Proca code, or the inf_chain `vmodes.c`, with the frozen-angle exit factor on a
  constant-H inflationary history;
- S/Q from the closed form;
- the scan in e, e_crit, the branch, and the relic;
- the H_i calibration;
- the isocurvature derivative d ln Ω_DM/d ln s_*.

**Output.** The map S/Q(K, m_s, t_on, e). This tells mf_chain and the critic what a fast spectator history needs.

### mf_chain: the full chain on mf_hybrid's benchmarks (starts after mf_hybrid and mf_ideal)

- Use mf_ideal's engine mode, reading `tracks/mf_hybrid/results/bench_*.json`.
- Run everything as in mf_ideal, plus:
  - the heavy remnants;
  - any hidden heat from s and σ decays;
  - gravitational production at the abrupt end.
- **Calibrate H_i self-consistently.** In the hybrid, V0 is free, so rescale and iterate until the H_i the dark matter
  needs equals the benchmark's own. Report the iteration.
- Give the branch at e = 10⁻⁶ and at 3×10⁻⁶, e_crit, and the relic. For Q-matter, give ⟨Q⟩ and M from the qm_micro
  machinery if needed.
- Give r, n_s and the isocurvature.

### mf_higgs: Higgs survival on mf_hybrid's benchmarks (starts after mf_hybrid)

- Use `tracks/inf_higgs/code` together with its skeptic's corrections:
  - N1: scan the field grid up to the abundance limit, not 20 h_max;
  - N2: cap with the inflationary H(N) criterion;
  - N5: vary α_s ± 0.0009;
  - N6: the H_i dependence.
- Cover:
  - inflation, with H nearly constant and any portal masses;
  - the waterfall and preheating, with the portal sign windows (G6) and R oscillations at 2m_χ;
  - plasma timing from the actual decay chain.
- Give the critical M_t and the ξ_H window.

## Rules

- Write only in `research/reality_framework/grut_working_theory_v2_20261101/tracks/<key>/` (NOTES.md, code/,
  results/; nothing over 2 MB; no compiled binaries or `__pycache__` left in the tree). Scratch goes in
  `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/<key>/`. Never write
  anywhere else, including `/`.
- Git is read-only. No commits.
- Use at most 2 threads or worker processes; other agents share the 4 cores. Run in the foreground (Bash timeout up to
  600000 ms per call; chunk long runs with restart files). Never end your turn with a background job running.
- **Labels and claims:**
  - Claim labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.
  - Every choice is a model choice, not a consequence of closure.
  - Never claim validation, exclusion of the real world, a new theory or a GRUT prediction.
  - Data comparisons are retrospective.
- **Citations.** Cite primary literature with arXiv numbers, and verify each one; for current limits use WebSearch via
  ToolSearch "select:WebSearch". Candidates to check:
  - hybrid inflation: Linde astro-ph/9307002; Copeland et al. astro-ph/9401011;
  - smooth and shifted hybrid: hep-ph/9506325, hep-ph/0002151;
  - tachyonic preheating: hep-ph/0012142, hep-th/0106179, hep-ph/0208228;
  - mild waterfall and PBH: 1006.4522, 1102.5612, 1501.07565;
  - the equation of state after preheating: hep-ph/0507096, 1608.01213;
  - the Higgs in inflation and preheating: 1210.6987, 1407.3141, 1602.07719.
- Return a digest of at most about 700 words.
