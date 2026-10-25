# Track `conversion`: independent check of 59a1a71 (conversion streaming) and the capture channel of bb211d1

Claude review, 25 October 2026. Read-only on Codex's branches. No commits.

Packets were read with `git archive` into the scratchpad:
- 59a1a71: `research/reality_framework/grut_conversion_streaming_20261023/`.
- bb211d1: `grut_core_conversion_20261023/`, plus its parents `grut_portal_debris_window_20261023/` and
  `grut_daughter_core_20261023/` (non-data files), at the same commit.

**Independence.** Every number marked [C] comes from my own code in `code/`, using methods different from
Codex's. The model inputs I take from the packets are quoted [Q] and listed in section 1. I did not run Codex's
code. The only Codex-derived numbers that pass through my code are its 32-row rate table: my transport was run
once on that table, labelled "Codex table", to separate table effects from transport effects. My own
independent rate table is labelled "my table".

Labels: [D] derived, [C] computed here, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

Scope. This is one model with supplied fields: the imposed s-decay source, frozen cold Q = 200000 ground
objects, the original Y_Q and a supplied background. Nothing here is evidence that the fields exist. The
structure comparison in section 8 is retrospective and approximate. It is not a likelihood, an exclusion or a
validation.

## 0. Verdict summary

All numbers in "Mine" are independent [C]. I ran no Codex code, so there are no reproduction runs. The "Codex
table" runs feed Codex's quoted 32-row rate table into my own transport code.

| # | Claim (packet) | Codex | Mine | Verdict |
|---|---|---|---|---|
| a | ω/m_X ≠ E/(Q m_X) at Q = 200000 (59a1a71) | 0.861195 / 0.923316 | 0.8611951308 / 0.9233161354, from my own BVP solve of the profile equations. dE/dQ = 0.8611951 by finite difference | CONFIRMED |
| cap | Conversion σv, cold and E = √1.5 (bb211d1) | 1.13039e-16 / 8.12250e-15 cm³/s | 1.13038e-16 / 8.12250e-15 from my own three-channel solver; unitarity 3e-10 | CONFIRMED |
| cap | 32-node rate table between nodes | PCHIP | 741-point scan: the narrow ℓ = 5 resonance (p = 0.508, σv = 7.69e-14) is under-sampled, and the table overestimates the integral over 0.45–0.56 by 1.43× | CORRECTED |
| b | Initial Γ/H, z = 10⁶ | ≈ 27 | 27.0066 | CONFIRMED |
| b | Half first-converted, z = 10⁶ | z ≈ 974487 | 974383 (my table) / 974487 (Codex table) | CONFIRMED |
| b | First-event survival separate from net anti-X | separated | separated: survival e^{−49.7} ≈ 1e-22, net anti-X ≈ 2.6e-7 | CONFIRMED |
| b | Net anti-X today, z = 10⁶ | 7.9e-7 (59a1a71), 1.08e-6 (bb211d1) | ≈ 2.6e-7; Codex's own grids converge towards ~1e-7 | CORRECTED |
| b | Median exit z, z = 10⁶ | 890452 | 889,170–889,455 (five MC runs, two tables) | CORRECTED (source-bin bias, minor) |
| b | Forward events, 20 MeV (bb211d1/59a1a71) | 39600–39831 | 34336 | CORRECTED (diagnostic only) |
| c | E_esc = 1 + 2Ω and u_esc | 2.722390 / 2.532076 | same; the MC mean q at escape matches to 0.02% | CONFIRMED |
| d | Horizons to today: 20 MeV / 1 keV / z = 10⁶ | 194.95 pc / 2.36919 / 8.37832 Mpc | 194.91 pc / 2.36921 / 8.3780 ± 0.001 Mpc | CONFIRMED. With a g*(T) background the 20 MeV value is ≈ 199 pc |
| d | Interpretation as streaming length | collisionless; flagged as uncertain | elastic σ_T keeps X coupled (Γ_T/H up to 76) until z ≈ 1.7e4 (z = 10⁶ source); the effective length is ≈ 4.2 Mpc | Codex's caveat CONFIRMED and quantified: the horizons are upper bounds |
| e | Mass fractions | 3.175 / 3.359 / 6.718 / 5.469% | 3.17513 / 3.35882 / 6.71764 / 5.46881%. Evolving ground state: 5.4781% (release +0.17%) | CONFIRMED |
| screen | Structure (new, approximate) | not done | z_source ≳ 3×10⁶ allowed; ≈ 10⁶ marginal-to-excluded; ≲ 3×10⁵ excluded (approx.) | see §8 |

## 1. Quoted inputs [Q]

| Input | Value | Source |
|---|---|---|
| Field equations (dimensionless), a = λ_X/g² = 0.1, b = λ_φ/g² = 3 | ψ̈ − ∇²ψ + (S² + a\|ψ\|²)ψ = 0; S̈ − ∇²S + [b(S²−1) + \|ψ\|²]S = 0 | bb211d1 MATHEMATICS.md §2 |
| g_X, v_φ, m_X = g v | 0.1; 0.924089637375362 GeV; 0.0924089637375362 GeV | PARAMETER_POINT.json |
| Charge and energy normalisation | g²Q = 4πΩ∫ρ²F²; E = (v/g)·ē | portal_debris profiles.py (definitions only) |
| Y_Q, s0, Q | 4.8294009076988175e-9; 2891.2 cm⁻³; 200000 | 59a1a71 streaming.py |
| Cosmology | H0 = 67.4, Ω_r = 9.20917977964354e-5, Ω_m = 0.315, flat | 59a1a71 streaming.py |
| Source epochs | a(20 MeV) = (T0/20 MeV)(3.909/10.75)^{1/3}; a(1 keV) = T0/1 keV; z = 10⁶ | 59a1a71 streaming.py |
| Late X per object; object mass; late total mass per object | 12521.802854958796; 17064.537454786467 GeV; 18221.721841256607 GeV | FIXED_YIELD_RESULTS.json (row Q = 200000) |
| Rate table (only for the "Codex table" runs) | 32 rows of σv(p/m_X) | bb211d1 RATE_RESULTS.json |

## 2. Profile, ω and E/(Q m_X): item (a)

**Method [C]** (`code/profile.py`). I solved the profile equations myself; I did not use the packet profile. The
solve is a fixed-frequency boundary-value problem in the regular variables u = ρF and z = ρ(S−1), so there is no
1/ρ term. The boundary conditions are u(0) = z(0) = 0, u′ + κu = 0 and z = 0 at ρ = 60, with tolerance 1e-10.
Charge and energy come from Simpson quadrature, and Ω for Q = 200000 from root finding. dE/dQ is a symmetric
finite difference of (E(Ω), Q(Ω)). The packet's profile CSV was integrated separately, as a data check only.

| Quantity | Codex [Q] | Mine [C] | Verdict |
|---|---|---|---|
| Ω = ω/m_X at Q = 200000 | 0.8611951308217225 | 0.8611951308218845 | CONFIRMED |
| E/(Q m_X) | 0.9233161354 | 0.9233161353725 | CONFIRMED |
| dE/dQ / m_X (finite difference, δΩ = 5e-5; error falls as δΩ²) | (claims = Ω) | 0.86119510 (differs from Ω by 3e-8) | CONFIRMED: ω = dE/dQ ≠ E/Q |
| Object energy | 17064.537454786 GeV | 17064.537454385 GeV | CONFIRMED |
| Virial (Derrick) residual | 3e-10 | 2.4e-15 | — |
| Packet CSV by my quadrature (data check) | — | Q = 200000.0000, E/(Qm) = 0.92331613537 | consistent |

So the chemical energy released per removed charge is ω, not E/Q [D, C]. Codex's correction from 0.923316 to
0.861195 is right.

**Frozen-profile check [C].** I also solved the ground state at Q′ = Q − 2N = 187478.197, where N = 6260.9014 is
the pair count per object.
- Ω rises to 0.8641633.
- E(Q) − E(Q′) = 998.202 GeV, against the frozen value 2ωN = 996.512 GeV (+0.170%).
- The outgoing energy at the end of complete conversion is 1 + 2Ω′ = 2.72833, against 2.72239 (+0.22%).
  u_esc rises to 2.53846, against 2.53208 (+0.25%).

These bound the frozen-profile error in the ground-state energy ledger at the 0.2% level, if the cores relax to
ground states. Evolving excited cores are not covered.

## 3. Capture/conversion channel of bb211d1: independent rates

**Method [C]** (`code/scatter.py`).
- Profile: my own Q = 200000 profile.
- Linearisation: I re-derived it from the field equations above, with δψ = e^{−iΩτ}[u e^{−iντ} + v* e^{iντ}],
  δS = w e^{−iντ} + c.c. and U = ρ(u, v, √2 w).
  - The u-equation gets aF²v + 2SF w; the w-equation gets SF(u+v).
  - With W = √2 w the coupling matrix is real symmetric, with V13 = V23 = √2 SF. This reproduces the matrix in
    bb211d1 MATHEMATICS.md §2 [D].
  - The √2 follows from equal number flux. For the model Lagrangian ½|∂ψ|² + ½(∂S)², a mode of S carries twice
    the energy density of a ψ mode of equal amplitude [D].
- Charge convention: an incoming anti-X of energy e is the negative-frequency part of the v channel
  (ν = e + Ω). The X goes out in the u channel with energy e + 2Ω. s is closed for e < √6 − Ω = 1.588 [D].
- Integration: adaptive DOP853 of the three regular solutions together. A regular series start uses
  M0ρ²/(2(2ℓ+3)). QR re-orthonormalisation every 1.5 units in ρ stabilises against the growing closed channel.
- Matching at ρ_m = 40: Riccati–Bessel decomposition of the open channels, and removal of the growing
  modified-Bessel component of s. The S matrix comes from flux-normalised in/out amplitudes, with
  σ = π(2ℓ+1)|S_uv|²/p² and v = p/e.
- This method differs from Codex's relaxation BVP.

| p/m_X | Codex σv [Q] (cm³/s) | Mine [C] | Rel. diff. |
|---|---|---|---|
| 1e-5 (cold) | 1.130389e-16 | 1.130384e-16 | 4e-6 |
| 0.1 | 1.363309e-14 | 1.363309e-14 | 3e-7 |
| 0.50474 | 6.386916e-14 | 6.386897e-14 | 3e-6 |
| 0.70711 = √0.5 (E = √1.5) | 8.122497e-15 | 8.122497e-15 | 2e-8 |

Per-ℓ probabilities agree too: e.g. at E = √1.5, P0 = 0.0749289 and P1 = 0.1332008 for both.

Numerical controls [C]:
- unitarity |S†S − 1| < 3e-10 over all scans;
- reciprocity |S_uv − S_vu| < 1e-11;
- moving ρ_m from 40 to 30, or rtol from 1e-11 to 1e-9, changes σv by < 5e-9.

Verdict: the **rate table entries are CONFIRMED** at the sampled momenta.

**Dense scan [C]** (709 momenta in steps of 0.001, plus 0.0005 steps across the ℓ = 5 resonance; fine ℓ = 6 and
7 scans in steps of 0.0002 — see section 3b):
- Between the table nodes, Codex's PCHIP interpolation agrees with my dense scan to within 0.05–0.15% in the
  integral of σv over 0 < p < 0.45 and over 0.56 < p < 0.708.
- Pointwise across p = 0.45–0.56 the interpolated table is off by factors 0.36 to 1.23. The narrow ℓ = 5 resonance
  peaks at p = 0.508, with σv = 7.69e-14, P5 = 0.80 and FWHM ≈ 0.012.
- Codex's single node at 0.5047 sits on the resonance flank, and PCHIP spreads it across the whole interval. The
  table therefore overestimates the integral of σv over 0.45–0.56 by a factor of **1.43**.
- Codex disclosed this risk: "narrow resonance interpolation not uniformly certified". It is now quantified.
  Effects on outputs are in sections 4 and 5.

## 3b. High-ℓ resonance search

I scanned single partial waves in steps of 0.0002 [C]: ℓ = 6 over 0.30–0.708 and ℓ = 7 over 0.40–0.708
(`code/ell_scan.py`, results/rate_scan.json). Neither has an interior maximum. Each rises only towards the band
edge: P6 = 0.0068 and P7 = 2.8e-5 at p = 0.708. For ℓ ≥ 8, P < 1e-9.

Within the source band the only narrow structure is the ℓ = 5 resonance at p = 0.508, and the 0.0005 grid
resolves it. A resonance narrower than 0.0002 in ℓ ≥ 6 is not excluded by these scans. Any such resonance would
carry P ≪ 1e-2 of rate weight.

## 4. First conversion, inverse conversion and the anti-X tail: item (b)

**Methods [C]** (`code/transport_mc.py`, `code/anti_tail.py`).
- **First-event survival.** Adaptive quadrature of τ = ∫ n_obj σv/H d ln a along the never-converted
  trajectory u = u0 a_f/a, split at the table knots. Inverse reactions cannot enter a first-event probability.
- **Population.** A grid-free Monte Carlo of individual particles in (species, comoving momentum q = a p):
  - exponential waiting times inside ln a steps, with ds = 0.004/max(Γ/H) and several events per step allowed;
  - each event conserves ν: E_+ = E_− + 2Ω;
  - the inverse rate follows from detailed balance;
  - a positive particle escapes at the analytic crossing a = q/u_esc.

  There is no energy grid and no remapping, so the method is unlike Codex's δ grid. Convergence: ds_fac 0.004
  against 0.001, and seeds 1, 2, 3.
- **20 MeV case** (Γ/H ≈ 3.2e6). A local-equilibrium drift ODE in ν: N_+/N_− = p_+E_+/(p_−E_−) and
  dν/d ln a = −⟨e − 1/e⟩.
- **Late net anti-X.** A first-order estimate f_anti ≈ ⟨∫ λ_+(δ) S_−(a, p_−(δ)) d ln a⟩ over the last approach
  of each escaping particle to the threshold. Here S_− is the probability that the produced anti-X survives to
  today. I take the whole band as an upper envelope.

| Quantity | Codex [Q] | Mine [C] | Verdict |
|---|---|---|---|
| Γ/H at source, z = 10⁶ | 27.0 | 27.0066 | CONFIRMED |
| Γ/H at source, 1 keV | 115 | 115.138 | CONFIRMED |
| Γ/H at source, 20 MeV | 3.23e6 | 3.2275e6 | CONFIRMED |
| Half first-converted, z = 10⁶ source | z ≈ 974487 | 974383 (my table); 974487 (Codex table, PCHIP); 974562 (Codex table, log-linear) | CONFIRMED (0.01% in z) |
| Total first-event optical depth to today, z = 10⁶ | 54.61 | 49.74 (my table); 54.61 (Codex table) | CORRECTED (−9%, from the ℓ = 5 resonance); survival is ~1e-22 either way |
| MC never-converted fraction against the analytic S0 | — | 0.99561 against 0.99601 at Δln a = 1.5e-4 (n = 1e5, first recorded step) | consistent |
| Median exit from the reactive band, z = 10⁶ | 890452 | 889167–889455 (five MC runs; 889167 with my table) | CORRECTED, minor: Codex's value carries a source-bin bias. Its own sequence 893742 → 891527 → 890452 extrapolates to ≈ 889400 |
| Median exit, 1 keV | 3854182 | 3849489–3849504 | as above (Codex sequence 3866442 → 3858788 → 3854182) |
| Escape, 20 MeV | median z = 1.0847e11 | a_esc = 9.2276e-12 (z = 1.0837e11), universal Δln a = 0.09612 in the equilibrium limit | CONFIRMED |
| Forward events per initial pair, z = 10⁶ | 1.30644 | 1.30885 (MC, Codex table PCHIP, n = 1e5); 1.2975 (Codex table, log-linear); **1.2657 (my table)**; 1 keV: 2.3846 (Codex) against 2.3508 (Codex table, log-linear) and **2.2029 (my table)** | Transport CONFIRMED on the same table (+0.18%, within MC noise). The table's over-broad resonance raises the counts by about 3% (z = 10⁶) and 8% (1 keV): a diagnostic CORRECTION |
| Forward events per pair, 20 MeV | 39600–39784 (grid-dependent); 39830.6 in bb211d1 | 34336 (my table, equilibrium); 39537 (Codex table PCHIP, equilibrium) | CORRECTED: about −13%. Diagnostic only; the excess comes from the over-broad ℓ = 5 resonance in the table |
| Net anti-X / initial anti-X today, z = 10⁶ | 7.9e-7 (59a1a71 finest); 1.08e-6 (bb211d1) | ≈ 2.6e-7 (first-order upper envelope; 2.61e-7 Codex table, 2.62e-7 my table). MC: 0/2e5 with the Codex table; 1/1e5 with my table, whose tail weight is 1.4e-7 | CORRECTED, see text |
| Net anti-X today, 1 keV | 2.7e-8 (finest; 8.5e-8 and 2.3e-7 on coarser grids) | ≈ 1.0e-8 | CORRECTED (same grid artefact) |

**Net anti-X tail.** Codex's value falls with grid resolution: 2.86e-6 → 1.47e-6 → 7.9e-7 over its three grids.
That is roughly first-order convergence, and a linear extrapolation gives ≈ 1e-7.

My estimate has two parts [C]:
- MC: no anti-X remain out of 2e5 once all positives have escaped, so the fraction is below 1.5e-5.
- First-order integral: ≈ 2.6e-7 per escaping particle.

The integrand peaks where an X close to threshold (δ ≈ 1e-4–3e-4) back-converts. That makes an anti-X with
p_− ≈ 0.014–0.025, whose survival is S_− ≈ 0.03–0.24, because σv falls steeply below p ≈ 0.02.

So the physical late anti-X fraction at the z = 10⁶ source is of order 1e-7. It is not ~1e-6; the ~1e-6 level is
numerical diffusion of the δ grid. Codex itself warned that this tail must not be read as a relic fraction. The
first-event survival to today is e^{−τ} ≈ 1e-22. It is a different, far smaller number, so Codex's separation of
the two quantities is right.

## 5. Kinematics and the escape receiver: item (c)

Inputs: my Ω (section 2), plus energy conservation with recoil neglected. The recoil is E_X/M_obj ≈ 1.5e-5.

| Quantity | Codex [Q] | Mine [C/D] | Verdict |
|---|---|---|---|
| Outgoing energy, cold anti-X: E_X/m_X = 1 + 2Ω | 2.7223903 (251.573 MeV) | 2.7223903 (251.573 MeV) | CONFIRMED |
| Escape momentum u_esc = √((1+2Ω)² − 1) | 2.532076 | 2.5320760 | CONFIRMED |
| Warm input E = √1.5: K_out, p_out | 1.947135, 2.772292 | 1.947135, 2.772292 | CONFIRMED |
| Comoving momentum at escape q/m = a_esc u_esc: mean, z = 10⁶ | 2.87809e-6 | 2.87752e-6 / 2.87795e-6 / 2.87817e-6 (three MC runs) | CONFIRMED (≤ 0.02%) |
| Mean q/m, 1 keV | 6.59499e-7 | 6.59505e-7 / 6.59487e-7 | CONFIRMED |
| q/m, 20 MeV | 2.33704e-11 (finest; coarser 2.33866, 2.33758) | 2.33650e-11 (equilibrium) | CONFIRMED (Codex converges towards mine) |
| Relative spread of q, z = 10⁶ / 1 keV | — | 3.6% / 1.0% | — |

[D] Inverse conversion closes exactly when the positive particle's energy falls below 1 + 2Ω, because the
anti-X it would make needs E ≥ m_X. Below that energy, this model's other inelastic channels are also closed:
- one-quantum absorption with s emission needs E_X − Ω ≥ √6, i.e. E_X ≥ 3.31;
- pair channels need more.

So u_esc is the same for every escaping particle, and q = a_esc u_esc is exact under free redshift. Codex's code
computes a_esc = a·u/u_esc inside the half-step (`streaming.py`, half_drift), which is the same exact crossing.

Number is conserved exactly in my MC, since every particle is followed. The core charge change is −2 × (net
forward events) and E_out = E_in + 2Ω per event, both by construction.

The frozen-profile caveat is quantified in section 2. If cores relax to ground states, u_esc rises by up to 0.25%
by the end of complete conversion, and horizons rise by at most that much.

## 6. Post-escape horizons: item (d)

**Method [C].** λ = c ∫ (q/a)/√(1+(q/a)²) da/(a²H), by my own quadrature for each escaped particle. The number-
weighted mean uses the per-particle escape times from the MC, or the equilibrium escape at 20 MeV. Controls:
- in pure radiation domination my quadrature equals the asinh closed form to 1e-15;
- the light-cone prefix is χ = c ∫_{a_f}^{a_esc} da/(a²H).

| Source | Codex mean to today [Q] | Mine to today [C] | Codex to z = 5 | Mine to z = 5 | Codex upper (with prefix) | Mine: mean + prefix |
|---|---|---|---|---|---|---|
| 20 MeV | 194.95 pc | 194.91 pc | 194.46 pc | 194.42 pc | 195.34 pc | 195.30 pc |
| 1 keV | 2.36919 Mpc | 2.36921 / 2.36915 Mpc | 2.35546 | 2.35548 / 2.35542 | 2.38105 | 2.38107 |
| z = 10⁶ | 8.37832 Mpc | 8.37695 / 8.37797 / 8.37844 / 8.37852 Mpc | 8.31841 | 8.31705 / 8.31806 | 8.44154 | 8.4401 / 8.4412 |

Verdict: **CONFIRMED** (agreement ≤ 0.02%, within MC scatter). With my own rate table the horizons are 8.3794 Mpc (z = 10⁶), 2.36929 Mpc (1 keV) and 194.91 pc (20 MeV). The resonance correction does not move the horizons (≤ 0.02%), because escape times barely change.

Sensitivity [C].
- **Background.** Codex uses the extrapolated Ω_r a⁻⁴. With entropy-conserving g*(T) through e⁺e⁻
  annihilation, ρ_r at 20 MeV is 0.830 times the extrapolated value, and the 20 MeV horizon becomes 199.1 pc
  (+2.2%). The 1 keV and z = 10⁶ cases are unchanged (< 1e-5).

  Codex disclosed this background choice, so this is a small CORRECTION to the 20 MeV number only: ≈ 199 pc.
- **Original positive daughters.** With q = a_f √0.5, their horizons are 51.1 pc (20 MeV), 0.644 Mpc (1 keV) and
  2.268 Mpc (z = 10⁶).
- **Epoch dependence.** Results/epoch_scan_codex.json gives the mean converted-X horizon against source
  redshift:

  | z_source | 1e5 | 3e5 | 1e6 | 3e6 | 1e7 | 3e7 | 1e8 | 1e9 |
  |---|---|---|---|---|---|---|---|---|
  | Horizon (Mpc) | 63.6 | 23.9 | 8.38 | 3.22 | 1.12 | 0.419 | 0.141 | 0.0171 |

  At z = 1e5 conversion is incomplete: 0.7% of the anti-X stay unconverted (Γ/H = 2.7 at the source).
- **Post-escape elastic scattering.** See section 6b. This is the largest open caveat on these horizons.

## 6b. Post-escape elastic scattering (new sensitivity; not in the packets)

**Method [C]** (`code/elastic.py`, `code/elastic_epochs.py`).
- After escape the anti-X channel is closed, so X–object scattering is single-channel and elastic. The v and s
  channels are closed and their growing parts are removed. The solver is the same as in section 3.
- I computed the momentum-transfer cross section σ_T = (4π/k²) Σ (ℓ+1) sin²(δ_ℓ − δ_ℓ+1) on the frozen
  Q = 200000 object.

| p/m_X | 0.01 | 0.1 | 0.2 | 0.4 | 0.7 | 1.0 | 1.4 | 1.8 | 2.2 | 2.5 |
|---|---|---|---|---|---|---|---|---|---|---|
| σ_T (m_X⁻²) | 341 | 243 | 167 | 289 | 120 | 35.1 | 9.0 | 3.1 | 1.5 | 2.0 |
| σ_el (m_X⁻²) | 343 | 332 | 390 | 534 | 484 | 291 | 149 | 84 | 50 | 36 |

At low momentum σ_el → 343 m_X⁻² ≈ 4π(5.2/m_X)², a hard-sphere value matching the object radius. That serves as
a plausibility check.

Γ_T/H = n_obj σ_T v/H along each escaped cohort [C]:

| Source | Γ_T/H at escape | Max Γ_T/H | Elastic decoupling z | Collisionless horizon | From decoupling | Diffusion length while coupled | Effective (combined in quadrature) | q_eff/q |
|---|---|---|---|---|---|---|---|---|
| z = 10⁵ | 0.62 | 5.8 | 6574 | 61.7 Mpc | 35.5 | 13.9 | 38.3 | 0.58 |
| z = 3×10⁵ | 2.2 | 21 | 10301 | 23.4 Mpc | 12.0 | 4.75 | 12.9 | 0.52 |
| z = 10⁶ | 7.7 | 76 | 17016 | 8.30 Mpc | 4.02 | 1.34 | 4.23 | 0.48 |
| 1 keV | 33 | 332 | 32394 | 2.36 Mpc | 1.11 | 0.31 | 1.15 | 0.47 |
| z = 10⁷ | 78 | 782 | 49031 | 1.11 Mpc | 0.52 | 0.13 | 0.54 | 0.46 |
| 20 MeV | 9.4e5 | — | 5.25e6 | 195 pc | 94 pc | — | — | — |

[E] In the frozen-object model the escaped X are **not** collisionless at first. While Γ_T/H > 1 they random-walk
among the objects, and they decouple elastically only at z ≈ 10⁴–10⁵. The comoving streaming length that
controls small-scale structure is then about half of Codex's collisionless horizon.

This is a sensitivity estimate, not a transport solution:
- angular kernels are summarised by σ_T;
- the random-walk/streaming split is crude;
- object number density and profile are frozen.

Codex explicitly said the horizons are "not a proof of all-interaction decoupling", and this confirms that the
concern is quantitatively important. The structure screen (section 8) therefore runs both a collisionless case and
an elastic case with q_eff.


My own arithmetic from the quoted per-object budget (`code/mass_fractions.py`, results/mass_fractions.json):

| Quantity | Codex [Q] | Mine, frozen profile [C] | Mine, evolving ground state [C] | Verdict |
|---|---|---|---|---|
| Initial anti-X rest mass / old calibrated matter | 3.175% | 3.17513% | — | CONFIRMED |
| Converted X / remaining rest matter | 3.359% | 3.35882% | 3.35915% | CONFIRMED |
| All free X / remaining rest matter | 6.718% | 6.71764% | 6.71829% | CONFIRMED |
| Core chemical energy released / old budget | 5.469% (996.512 GeV) | 5.46881% (996.512 GeV) | 5.47809% (998.202 GeV) | CONFIRMED (the frozen-profile value is low by 0.17%) |
| Core charge transferred | 6.26% | 6.2609% | — | CONFIRMED |

The budget closes: object 17064.537 + free X 1157.127 + D 0.058 = 18221.722 GeV. "Remaining rest matter" means
the whole dark-sector rest mass after the kinetic energy has redshifted away; baryons are not included. The
abundance calibration therefore loses 5.47–5.48% unless the yield is retuned. Codex did not retune it.

## 8. First structure-data screen (retrospective, approximate)

This is a screen, not a likelihood. It is retrospective and does not claim any exclusion. It assumes:
- the converted X are a dark-matter fraction f = 3.359% (or 6.718% with the colder original daughters), with
  the computed escape-momentum distribution;
- linear theory and adiabatic initial conditions;
- no post-escape interactions in the "collisionless" case; the elastic case mimics coupling with q_eff;
- the imposed source and frozen objects.

### 8a. Published constraints used (read levels)

| Constraint | Result used [Q] | Read level |
|---|---|---|
| Garcia-Gallego, Iršič, Haehnelt, Viel, Bolton 2025, PRD 112, 043502 (arXiv:2504.06367): Lyα, UVES/HIRES, z = 4.2–5.0, mixed CDM+thermal WDM | f_WDM < 0.16 / 0.35 / 0.50 / 0.67 (2σ) at m = 1 / 2 / 3 / 4 keV; **no constraints below 1 keV** (emulator range) | abstract + results text (HTML) |
| Boyarsky, Lesgourgues, Ruchayskiy, Viel 2009, JCAP 05, 012 (arXiv:0812.0010): WMAP5 + SDSS/VHS Lyα, mixed | any mass allowed if F_WDM < 0.35 (Bayesian 95%) or < 0.6 (frequentist 99.7%), for 5 keV < m_NRP < ∞ only | abstract |
| Diamanti, Ando, Gariazzo, Mena, Weniger 2017, JCAP 06, 008 (arXiv:1701.03128): Planck15 + BAO + MW satellite count; thermal relic with T = T_ν, so only m/T matters | 2σ (fermion): f ≤ 0.29 (1–10 keV), 0.43 (10–100 keV) [text]; ≈ 0.16 (100 eV–1 keV) and ≈ 0.04 (10–100 eV) [read from their Fig. 5 at 300 dpi, ±20%] | full text §4–5 + figure |
| Tan, Dekker, Drlica-Wagner 2025, PRD 111, 063079 (arXiv:2409.18917): DES + PS1 Milky Way satellites, mixed | f_WDM bound plateaus at ≲ 0.45 for m ≲ 1.5 keV; pure WDM m > 6.6 keV | abstract |
| Iršič et al. 2024, PRD 109, 043511 (arXiv:2309.04533); Nadler et al. 2021, PRL 126, 091101 | pure thermal WDM m > 5.7 keV (Lyα), > 6.5 keV (satellites) | abstract (search summary) |
| Murgia, Merle, Viel, Totzauer, Schneider 2017, JCAP 11, 046 (arXiv:1704.07838) | area criterion δA (definition used here); δA_ref = 0.38 (3.5 keV), 0.21 (5.3 keV) in their pipeline | §5.2 (HTML) |
| Rogers & Poulin 2025 (Phys. Rev. Research; arXiv:2311.16377v3): eBOSS DR14 P1D compressed | Δ²_* = 0.310 ± 0.020, n_* = −2.340 ± 0.006, ρ = 0.512 at k_* = 0.009 s/km, z = 3; 4.9σ from Planck ΛCDM; a mixed WDM with f ≈ 0.022, m ≈ 10 eV (T/T_cmb ≈ 0.19) eases it | full text (HTML, summarised) |
| DESI DR1 P1D cosmology (arXiv:2601.21432, 2026) | Δ²_* = 0.379 ± 0.032, n_* = −2.309 ± 0.019 (consistent with Planck); correlation not used | abstract (search summary) |
| Decaying DM to a massive daughter with a kick: Wang, Croft, Peter, Zentner, Purcell 2013, PRD 88, 123515 (arXiv:1309.7354); Mau et al. 2022, ApJ 932, 128; Fuß & Garny 2023, JCAP 10, 020 (arXiv:2210.06117) | v_kick ≲ 30–70 km/s for τ ≲ 10 Gyr (SDSS Lyα); τ > 18 (29) Gyr for v_kick = 20 (40) km/s (satellites, 95%); BOSS DR14 bound τ ≳ 18 Gyr for small splittings | abstracts |
| Decant, Heisig, Hooper, Lopez-Honorez 2022, JCAP 03, 041 (arXiv:2111.09321): superWIMP / freeze-in Lyα reinterpretation | method context for decay-produced (kicked) DM, including mixed cases | abstract |

[A] How the decaying-DM bounds map here. Those bounds are for late decays, τ ≳ Gyr, of essentially all the DM,
with kicks of tens of km/s. In this model a 3.4% fraction receives a relativistic kick (p/m = 2.5) at t ≲ 1 yr
for the z = 10⁶ source. The kick then redshifts, to a velocity today of c·q = 0.86 km/s (z = 10⁶) or 0.20 km/s
(1 keV).

So the decaying-DM limits do not apply directly; the right analogue is a mixed warm+cold component with a
non-thermal, nearly monochromatic momentum distribution. Linear CLASS P(k) with that distribution is used below,
rather than a thermal-mass relabelling.

Equivalent thermal mass (matching ⟨p⟩/m, T = T_ν, as in Diamanti et al.) [C]: m_eq = 3.151 T_ν / ⟨q/m⟩. The
values below are for the collisionless case; the elastic case roughly doubles m_eq.

| Source | m_eq |
|---|---|
| z = 10⁵ | 14.5 eV |
| z = 3×10⁵ | 51 eV |
| z = 10⁶ | 184 eV |
| 1 keV | 804 eV |
| z = 10⁷ | 1.9 keV |
| z = 10⁸ | 19 keV |
| 20 MeV | 23 MeV |

### 8b. CLASS screen (own code: `code/class_screen.py`; CLASS v3.3.2 via classy, Lesgourgues 2011; Blas, Lesgourgues & Tram 2011; ncdm: Lesgourgues & Tram 2011)

**Pipeline and calibration [C].**
- Base cosmology: h = 0.674, ω_b = 0.0224, ω_dm = 0.1207 (total), N_ur = 3.046, A_s = 2.1e-9, n_s = 0.965.
- The converted X are one ncdm species with a tabulated, strictly positive psd: Gaussian at 10% width for narrow
  samples, a kernel density for broad ones. P(k) changes by < 3e-4 between 3.6% and 10.7% width, or between 60
  and 120 momentum bins.
- The original daughters (q = a_f√0.5) are a second ncdm species.
- P_k_max = 120 h/Mpc, with a power-law tail in the P1D integral.

Murgia's δA, evaluated with my pipeline on thermal relics, gives:
- pure WDM: 0.465 (3.5 keV), 0.332 (5.3 keV), 0.312 (5.7 keV);
- the 2σ mixed bound of Garcia-Gallego et al., 1 keV and f = 0.16: 0.367 when the warm temperature is set as if
  it were all of the DM, 0.281 when it is set by its own abundance;
- 2 keV, f = 0.35: 0.419 / 0.337.

My absolute δA runs about 0.08–0.12 above Murgia's quoted 0.38 and 0.21 (truncation and k-range choices). I
therefore compare against my own calibrators.

Compressed Lyα parameters at k_* = 0.009 s/km ≈ 1.03 h/Mpc, z = 3: my ΛCDM gives Δ²_* = 0.3613 and
n_* = −2.3037. Each model's shift is scored against DESI DR1 with χ² (correlation not used) and reported as
Δχ² relative to ΛCDM. The eBOSS Δχ² is listed for context only, because Planck ΛCDM is already 4.9σ off eBOSS.

**Classification rule (approximate):**
- allowed: δA < δA(5.7 keV) and Δχ²_DESI < 4 and f below the mixed-WDM bound for m_eq;
- marginal: δA up to δA(3.5 keV), or 4 ≤ Δχ² < 9, or f within 1.5× of the bound;
- excluded (approx.): anything beyond those.

Full rows are in results/structure_screen_summary.json.

**Collisionless after escape** (the converted X only, f = 3.359%):

| z_source (T) | q/m today | m_eq | δA | P(z=3) ratio at k = 1 / 5 / 20 h/Mpc | Δ²_* ratio | Δn_* | Δχ²_DESI | Δχ²_eBOSS | Screen |
|---|---|---|---|---|---|---|---|---|---|
| 1.2e11 (20 MeV) | 2.3e-11 | 23 MeV | ≈ 0 (195 pc) | 1, 1, 1 | 1 | 0 | 0 | 0 | allowed (cold) |
| 1e9 | 2.8e-9 | 190 keV | 0.003 | 1.000, 1.000, 1.000 | 1.000 | 0.000 | 0.0 | 0.0 | allowed |
| 1e8 | 2.8e-8 | 19 keV | 0.031 | 1.000, 1.000, 0.991 | 1.000 | 0.000 | 0.0 | −0.1 | allowed |
| 3e7 | 9.3e-8 | 5.7 keV | 0.082 | 1.000, 0.995, 0.935 | 1.000 | −0.0002 | 0.0 | −0.5 | allowed |
| 1e7 | 2.8e-7 | 1.9 keV | 0.128 | 0.999, 0.966, 0.899 | 0.999 | −0.002 | 0.0 | −5.0 | allowed |
| 4.26e6 (1 keV) | 6.6e-7 | 804 eV | 0.167 | 0.995, 0.912, 0.857 | 0.994 | −0.013 | +0.2 | −21.5 | allowed |
| 3e6 | 9.4e-7 | 564 eV | 0.184 | 0.990, 0.906, 0.839 | 0.989 | −0.025 | +1.1 | −30.6 | allowed |
| 1e6 (235 eV) | 2.9e-6 | 184 eV | 0.232 | 0.938, 0.856, 0.783 | 0.935 | −0.090 | +21.3 | +92.8 | excluded (approx.), by the compressed P1D |
| 3e5 | 1.0e-5 | 51 eV | 0.273 | 0.880, 0.791, 0.739 | 0.879 | −0.071 | +15.2 | +11.2 | excluded (approx.) |
| 1e5 | 3.7e-5 | 14.5 eV | 0.299 | 0.810, 0.741, 0.721 | 0.808 | −0.069 | +18.2 | −2.8 | excluded (approx.) |

**With elastic coupling (q_eff, section 6b):**
- z = 10⁶: δA = 0.201, Δn_* = −0.046, Δχ²_DESI = +4.9 → **marginal**.
- z = 3×10⁵: δA = 0.253, Δχ² = +4.0 → marginal (borderline).
- z = 10⁵: Δχ² = +16 → excluded (approx.).
- All z ≥ 3×10⁶: allowed.

**All free X** (converted plus original daughters, f = 6.718%), collisionless:
- z ≥ 3×10⁶: allowed. δA rises to 0.25 (1 keV) and 0.28 (3×10⁶), near the 5.7 keV calibrator 0.31.
- z ≤ 10⁶: excluded (approx.), with Δχ²_DESI = +29 at 10⁶.

With elastic coupling the all-free-X case is allowed for z ≥ 3×10⁶, marginal at 10⁶ (Δχ² = +6.0, δA = 0.31),
and excluded at ≤ 3×10⁵.

Other observables [C]:
- σ8 changes by ≤ 0.14% for z ≥ 10⁶, and by −1.1% (3×10⁵) and −4% (10⁵).
- CMB lensing C_ℓ^φφ is reduced by ≤ 0.1% at ℓ ≤ 400 and by 1.8% at ℓ = 2000 for z = 10⁶. The 1 keV case gives
  0.04% and 0.15%; the z = 3×10⁵ case 0.8% at ℓ = 400 and 7.5% at ℓ = 2000.
- Lensed TT changes by < 4e-4.

Planck and ACT lensing are measured to about 2–3%, so CMB/lensing does not constrain any case except possibly
z ≤ 3×10⁵ at high ℓ.

Comparison with the literature bounds, by f and m_eq:
- f = 3.4–6.7% is far below the satellite plateau (f ≲ 0.45; Tan et al. 2025) at every epoch.
- It is below the Lyα mixed bounds where they apply, m ≥ 1 keV (Garcia-Gallego et al. 2025: 0.16 at 1 keV).
- It is below Diamanti et al.'s ≈ 0.16 for the 100 eV–1 keV decade.
- At 10–100 eV, f = 3.4% sits at Diamanti's ≈ 0.04 limit, so the z ≈ 3×10⁵ source is marginal and the
  z ≈ 10⁵ source (m_eq ≈ 15–26 eV) is excluded.
- No published mixed bound exists for the non-thermal, nearly monochromatic shape. The z = 10⁶ step in P(k)
  falls right at k ≈ 1 h/Mpc, which is why the compressed P1D slope n_* is the most sensitive quantity here
  (Δn_* = −0.09 against σ(n_*) = 0.019).

### 8c. Screen outcome by production epoch (retrospective, approximate)

| Imposed source epoch | Converted X only (3.36%) | All free X (6.72%) |
|---|---|---|
| z_source ≳ 3×10⁶ (T ≳ 0.7 keV), including the 1 keV and 20 MeV sources | **allowed** | **allowed** (near the δA boundary at 1 keV–3×10⁶) |
| z_source ≈ 10⁶ (T ≈ 0.24 keV) | **excluded (approx.)** if collisionless; **marginal** with elastic coupling | excluded (approx.) collisionless; marginal with elastic coupling |
| z_source ≈ 3×10⁵ | excluded (approx.) collisionless; marginal with elastic coupling | excluded (approx.) |
| z_source ≲ 10⁵ (T ≲ 25 eV) | **excluded (approx.)** | excluded (approx.) |

These labels depend on:
- a two-number compressed P1D comparison that ignores the DESI Δ²–n correlation, applied to a feature that is
  not ΛCDM-shaped;
- the linear δA screen;
- the frozen-object elastic estimate.

They are not exclusions.

eBOSS, which itself sits 4.9σ from Planck ΛCDM, would move the z ≳ 3×10⁶ cases towards agreement
(Δχ² < 0) and the z = 10⁶ collisionless case away from it. This is context, not evidence for the model.


## 9. Reproduction

Set S=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/conversion,
then extract the packets read-only:

    git archive 59a1a71 research/reality_framework/grut_conversion_streaming_20261023 | tar -x -C $S/p59
    git archive bb211d1 research/reality_framework/grut_core_conversion_20261023 | tar -x -C $S/pbb
    git archive bb211d1 research/reality_framework/grut_portal_debris_window_20261023 | tar -x -C $S/parents

Then run these, with C = this track's code/:

    python3 -I $C/profile.py $S/parents/.../profile_Q200000_lambda0.03_refined.csv   # section 2; writes my profile to $S
    python3 -I $C/scatter.py spot                                                    # section 3 spot rates
    python3 -I $C/scatter.py scan 0.001 0.236 236 a   (also b: 0.237 0.472 236; c: 0.473 0.708 236; refine1: 0.500 0.516 33)
    python3 -I $C/ell_scan.py 6 0.30 0.708 2041 ; python3 -I $C/ell_scan.py 7 0.40 0.708 1541
    python3 -I $C/transport_mc.py all codex ; python3 -I $C/transport_mc.py all mine quick
    python3 -I $C/anti_tail.py codex ; python3 -I $C/epoch_scan.py codex
    python3 -I $C/mass_fractions.py
    python3 -I $C/elastic.py ; python3 -I $C/elastic_epochs.py
    $S/venv/bin/python -I $C/class_screen.py cal|model|elastic      # pip install classy==3.3.2.0 in a scratch venv
    python3 -I $C/summarize_screen.py

Python 3.13, NumPy 2.5.3, SciPy 1.18.1, classy 3.3.2.0. Monte Carlo seeds are fixed. Large intermediates
(escape-time samples, scans, CLASS psd files) are in $S. Results JSON are in results/.

## 10. What remains uncertified

- **The source.** The actual soft formation source is still missing: the imposed instantaneous s-decay source,
  Y_Q and Q = 200000 are supplied. General hot-core leakage and vacuum density are also open.
- **Frozen versus evolving objects.** Ground-state relaxation is bounded at 0.2–0.25% (section 2). Excited
  newborn cores and the mixture of object charges are not treated.
- **Post-escape transport.** Elastic scattering roughly halves the streaming length (section 6b). A real transport
  calculation with angular kernels and a coupled-fluid perturbation treatment is needed.
- **Structure.** Section 8 is a linear, approximate screen. There is no Lyα hydrodynamic likelihood, no
  marginalisation over the thermal history, and no named full data likelihood. The compressed-parameter
  comparison assumes a ΛCDM-like shape near k_*.
- **Rate table.** Resonances in ℓ ≥ 6 narrower than 0.0002 in p/m_X are not excluded. Weak-wave rates on a cold
  ground object are not rates for hot newborn cores.
