# Brief: complete the end-to-end model ("M2")

The owner asked us to "complete the end-to-end model". Read [WORKING_THEORY.md](WORKING_THEORY.md) §0, §2, §4b, §7–§9,
[BRIEF_MULTIFIELD.md](BRIEF_MULTIFIELD.md) and `tracks/mf_critic/NOTES.md` first; the earlier rules apply (end of
this file). All new tracks are `tracks/e2e_<key>/`.

## What "complete" means here

One fully specified model, called **M2**, is M1r's hidden sector with the multi-field inflaton sector of the
multi-field pass (HFa/HFb). It needs three things.

**1. A parameter card.** Every field, every coupling and every initial condition, each marked:
- **supplied**, a model choice;
- **fitted**, to a named datum;
- **derived** or **computed**.

Nothing may be left unspecified. The holder strength c0 and phase, e and s_* are supplied numbers. Say so, and give
the window each must sit in.

**2. A computed history from before CMB horizon exit to today,** with no uncomputed link:
- inflation;
- the waterfall and preheating;
- reheating and thermalisation;
- the hidden release;
- relaxation and the branch;
- objects or free X;
- the Standard Model (SM) side: neutrino masses, the baryon asymmetry, BBN.

**3. An observables ledger.** For each quantity give the model's value and the relevant bound or measurement, quoted
retrospectively with arXiv numbers:
- CMB: n_s, r, α_s, isocurvature, f_NL;
- dark matter: Ω_DM, the relic type, the object mass, size and flux, the self-interaction, P(k);
- Σm_ν, the baryon asymmetry η_B, ΔN_eff, BBN;
- gravitational waves from the waterfall and preheating;
- the absence of defects and PBHs;
- Higgs stability.

**Two benchmarks:**
- **M2-Q:** Q-matter at the declared e = 10⁻⁶.
- **M2-X:** M1r's two-branch picture, free X at 10⁻⁶ with e_crit ≈ 2×10⁻⁶.

Both are tuned choices of c0. Report the tuning honestly. Completing the model does not make it a prediction of
GRUT: closure fixes none of the supplied numbers.

## Where things stand (multi-field pass)

- **The sectors and their history:**
  - HFa / HFb: `tracks/mf_hybrid/results/bench_HFa.json`, `bench_HFb.json`; code in `tracks/mf_hybrid/code`.
  - The Pauli-blocked post-waterfall fluid: `tracks/mf_hybrid/skeptic/code/sk_pauli.py`.
  - The history builder and H_i iteration: `tracks/mf_chain/code/hist_build.py`, `iterate.py`, `cmb_rescale.py`.
- **The engine** (MODE 5, spectator holder, tabulated background): `tracks/mf_ideal/code/rpmfi.c`, plus `bench_run.py`,
  `vectors_mfi.py`, `vmodes.c` and README.
- **The holder-strength scans:** `tracks/mf_critic/code/hyb_c0.py`, `m1r_c0.py`, and
  `tracks/mf_chain/skeptic/code/sk_c0chain.py`.
- **The strong-holder points** (Pauli Γ_th = 0 history, HFa):
  - c0 = 1.40: H_i = 8.40×10¹⁰ GeV, Q-matter at 10⁻⁶, S/Q = 2.88;
  - c0 = 1.36: H_i ≈ 8.85×10¹⁰ GeV, M1r's picture.
  - Window: ±3% in |μ_eff| and ±5° in phase.
  - At Γ_th = 10 H_i the points move to c0 ≈ 1.25–1.30 at H_i ≈ 5×10¹⁰ GeV, near the axial effective theory's cutoff.
- **Open links (mf_critic §5; both skeptics):**
  1. The post-waterfall history needs a lattice: walls, fragmentation, K. It also needs the thermalisation rate Γ_th
     and the Pauli blocking computed, not bracketed.
  2. The pair-fraction map has an unsampled gap at S/Q 6.7–17. ⟨Q⟩ and the vapour must be recomputed on the hybrid's
     thermal history.
  3. The SM-side completion:
     - Σm_ν = 0.20 eV for four singlets at m̃ = 0.05 eV;
     - the N_s soft mass;
     - the UV completion of the axial couplings: periodic couplings on non-periodic fields;
     - the anomalous SU(2) couplings of the rotating axial phases;
     - the baryon asymmetry, which nothing has computed yet.
  4. The induced holder-portal condition (c_L ≲ 1/16π², or vanishing at M_Y).
  5. Isocurvature: s_* ≥ 1.4×10¹⁶ GeV at the strong holder.

## Tracks

### e2e_lattice: the post-waterfall history on a 3D lattice

Lattices are allowed in this track, up to 64³; a 96³ convergence check is allowed if it fits in about 2 h per run on
2 threads.

**Fields.** Evolve σ and the real biased χ of HFa (and HFb if the parameters differ) in an expanding background with
self-consistent a(t), from just before σ_c until the coherent stage ends. Include:
- the spectator s (HFa: mass g_s χ; HFb: ξ_s R), so that its mass and value can become inhomogeneous;
- the energy drain into the axial singlet fermions as a local, Pauli-limited transfer term. Bracket it with the
  mf_hybrid skeptic's saturation model, and state the closure.

**Outputs:**
- Do walls form and percolate? Measure the fraction of the −v vacuum and the wall area per volume. Extrapolate with
  the bias pressure to their decay or dominance.
- Fragmentation and the gradient-energy fraction.
- The equation of state w(t) and the energy fractions.
- ⟨R⟩(t) and its fluctuation.
- The spread of s and m_s over the box.
- A gravitational-wave estimate: peak frequency and Ω_GW today. Cite verified literature (e.g. Dufaux et al. on
  preheating GWs).

Export a history `results/hist_<name>.json` that the chain can splice: t, a, H, R, the fractions, and s(t) or m_s(t)
with its variance.

**Convergence.** Check box size, lattice spacing and time step.

### e2e_kinetic: thermalisation and K

**What to compute.** Take the singlet-ν reheating of HFa/HFb and compute the real transfer of energy into an SM
thermal bath. Do not bracket it. Cover:
- Pauli-blocked production by the oscillating axial phases;
- N decays N → LH with the seesaw Yukawas of `e2e_sm`'s neutrino completion. Use the mf_hybrid values if e2e_sm is
  not ready; the dependence must be stated;
- the SM cascade and thermalisation, with the LPM-suppressed splitting rates (Harigaya–Mukaida; Kurkela–Moore; verify
  the arXiv numbers).

**Outputs:** Γ_th(t), f_rad(t), K, T_RH, and the time the SM bath reaches 50% and thermal equilibrium. Write
`results/fluid_<name>.json` in a form the chain can combine with the lattice's coherent-stage history.

### e2e_sm: the Standard Model side

**1. Neutrinos.**
- Choose a seesaw completion that gives the measured Δm²₂₁ and Δm²₃₁ with Σm_ν below current bounds. Quote NuFIT and
  DESI/Planck retrospectively, verified.
- Use a complex Casas–Ibarra R, or fewer or differently weighted singlets.
- Keep the N decay fast enough for reheating, and keep the Higgs medium mass that mf_higgs's skeptic relied on.
- Fix N_s: give its soft mass and decay, or state that it is dark radiation with its ΔN_eff.

**2. UV completion of the axial couplings.** For example, each N gets its mass from a complex scalar
Φ = (ρ/√2) e^{iφ/f}. Say what φ is, and how σ and χ, the fields actually coupled, enter. Then:
- check that the completion keeps the one-loop flatness and the CP bookkeeping, with the holder as the only s-odd term;
- check that it introduces no new light states or defects.

**3. Baryogenesis.** Compute the baryon asymmetry with Boltzmann equations and verified rates. Two sources:
- leptogenesis from the N decays, non-thermal or thermal, given T_RH and M_N, with an explicit CP phase;
- the anomalous SU(2) couplings of the rotating phases, which can give spontaneous lepto- or baryogenesis.

Give η_B, compare it with the measured value, and say whether it is reached, overproduced or underproduced and what
fixes it.

**4. Other checks.** BBN; ΔN_eff; any new long-lived state; washout; sphaleron conversion.

### e2e_relic: the dark-matter relic machinery

**1. The pair-fraction gap.** Fill the gap in p(S/Q) between 6.7 and 17 with the reviewed case-L kinetic network, and
give the corrected map as a function. See `tracks/mf_chain/skeptic/code/sk_pmap.py` and the inf_chain p-map.

**2. The Q-matter relic as a function of the inputs.** Inputs: S/Q, H_i, the hidden-to-visible entropy ratio and
e. Compute ⟨Q⟩, the vapour fraction and the mass per charge M, using `tracks/qm_micro` and its skeptic. Give a
callable function or table that the chain can evaluate on the final history.

**3. Observables of the relic:**
- object mass, radius and number density;
- the flux on Earth and the event rates for proposed gravitational or impulse detectors, cited and verified;
- the self-interaction;
- P(k) via `tracks/qm_transfer`;
- the free-X relic properties for M2-X.

### e2e_chain: the final self-consistent chain (after e2e_lattice, e2e_kinetic, e2e_relic and the history skeptic)

**1. The final history.** Splice the lattice's coherent stage with the kinetic track's thermalisation into final
histories for HFa and HFb. Apply the history skeptic's corrections where they are concrete.

**2. The holder scan.** On each final history, scan c0 and the phase with the full chain:
- the mf_ideal engine;
- the vectors with the benchmark's own H(N) and F(N);
- the ledger, including the leak with the c_L condition stated;
- the closed form;
- e2e_relic's corrected p-map.

Find the self-consistent points M2-Q and M2-X. "Self-consistent" means the H_i the dark matter needs equals the
benchmark's own, through the V0 rescaling and with N_* re-solved. Give the window in c0 and phase for each, and the
H_i, r and CMB capacity at each.

**3. At M2-Q and M2-X:**
- S/Q(e) and e_crit;
- the relic, using e2e_relic's function;
- the heavy remnants;
- the gravitational production;
- the isocurvature, with the s_* it requires;
- the s energy fraction and decay;
- the n_s / r / α_s from Mukhanov–Sasaki modes at the final scale;
- a hydrostatic check that the spectator spread from the lattice does not break the release (inhomogeneous m_s).

**4. Is c0 natural?** If s_* comes from de Sitter equilibrium (the Starobinsky–Yokoyama stochastic distribution) for
the residual inflationary mass of s, what is the distribution of c0, and how likely is the window? If equilibrium
cannot reach s_* ~ 10¹⁶ GeV, say what then sets s_* (an initial condition from before the observable e-folds).

### e2e_higgs: Higgs survival on the final histories (after e2e_lattice, e2e_kinetic and the history skeptic)

Use `tracks/mf_higgs/code` with both skeptic corrections: the exact ε_H = 1 endpoint, and the N-gas medium mass. Cover
inflation at the final H_i, the lattice's ⟨R⟩ and its fluctuations, and the kinetic track's plasma. Give the critical
M_t, the ξ_H window, and the α_s and M_h dependence.

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
  - Data comparisons are retrospective. Cite primary literature with arXiv numbers you have verified; use WebSearch via
    ToolSearch "select:WebSearch".
- Return a digest of at most about 700 words.
