# m1_critic: completeness critic for Claude's review of Codex's GRUT-M1 chain (issue #13)

1 November 2026 · Claude · same-session AI work, not outside replication. Inputs: BRIEF_M1.md, the four review tracks
(`m1_vacuum`, `m1_ledger`, `m1_holder`, `m1_transport`) with their skeptic passes, and Codex's packets read through
`git show` only (d4747cb CLAUDE_HANDOFF/CANDIDATE_M1/CONSEQUENCES; ed73791 REPORT; c481530 REPORT; 94f660d REPORT;
f8186fa REPORT/FOLLOWUP_NOTES/STABILIZED_CANDIDATE). Nothing checked out, nothing committed. Own quick checks in
`code/critic_checks.py` (algebra and 1D quadrature, 0.2 s, one thread); results in `results/critic_checks.json`,
`results/checkpoint_kurtosis_repro.json` (labelled reproduction from Codex's npz) and `results/verdicts.json`.
Remote re-checked today (`git ls-remote origin`): f8186fa is the newest `codex/grut-m1-*` head; the brief's "newer
isocurvature packet" (4e-12 / 3.02 / 3.435e-70 / 38 checks) does not exist on the remote.

Everything below concerns one supplied working model. Nothing here validates or excludes GRUT-M1 or any theory, and
all data comparisons are retrospective. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued,
[R] recalled.

## 1. Coverage of Codex's bounded request and of the later load-bearing claims

Codex's request (d4747cb): "refute the one-vacuum argument and scalar BVP/Hessian first; then check the new
mass-normalization chain and included charge-loss normalization."

| Item | Covered by | How well | Gap / settled here |
|---|---|---|---|
| One-vacuum argument (V', discriminant, threshold, V ≥ 0, hidden tree bound, no walls) | m1_vacuum + skeptic | exact rationals, brute-force bound test, wall BVP | complete. Skeptic N5: Codex's own ε grid straddles the one-vacuum threshold (ε = 0.001 has a metastable second minimum); immaterial for ε = 0.002 |
| Scalar BVP and Hessian (Q = 2e5 object, first law, fixed-Q Hessian, zero modes) | m1_vacuum + skeptic | two independent solvers agree with Codex to 1e-13; Hessian rebuilt two ways | complete, with one correction (the "positive eigenvalue" is the free-X box mode; the right statement is VK, which holds). NEW: Q_min = 5.67e4, binding threshold E/Q = m_X at Q = 7.13e4 |
| Mass-normalisation chain (Y_Q, λ = 1.00700905, Y ∝ H_i^{3/2}, old/new targets) | m1_ledger + skeptic | own engine at τ = 800 and 50-digit arithmetic on Codex's rows | complete. The brief's "+3.42 % needed" was a misreading (the old point already gave 1.0408 Y_Q,old). Precision overstated by Codex and track: λ is defined only to ±1e-2 (clock, history, D rounding; §3) |
| Charge-loss normalisation Γt₀ = 6.27e-17 | m1_ledger + skeptic | Born integral by three quadratures (2868466.51); A ∝ λ⁻² | complete for the vacuum-Born channel; coherent channel (6e-5) quoted [Q] |
| Holder Hessian and charge change (f8186fa) | m1_holder + skeptic | three independent codes on the ODE; mpmath point and Hessian | complete; cause of the late drift settled here (§2) |
| Isocurvature cancellation | m1_holder + skeptic | theorem written out; own nonlinear restarts; exact Hankel power table | the argument and the f8186fa table are checked; the specific 4e-12 / 3.02 / 3.435e-70 numbers are NOT CHECKED because the packet is absent |
| Transport forecast (σ, rates, free streaming, BGK curve, observability) | m1_transport + skeptic | own cross sections, rates, backgrounds; fluid toy for the BGK magnitudes; coupled X–core solve | σ, rates, free streaming complete; BGK hierarchy itself NOT CHECKED (toy agrees 1–3 %); clock-scaling sign settled here (§2); observability assessed [R] |
| **Gas source** (ed73791 3D pulses 0.002/0.008/0.030/0.095 Q; c481530 chain 1695 X̄ → 13189 free X, ×7.37; bf8ab3b three-channel X̄ scattering) | nobody (3D excluded by the budget; the kinetics chain only re-read as ledger arithmetic) | — | **uncertified root of the whole gas/transport chain**: the 6.6 % and v_rms = 0.068 c rest on a supplied pulse and an unreviewed chain |
| 94f660d loose ends: H/m_X = 1.9301e-21 at 20 MeV, refit −0.567645 %, sextic +9.927 % charge, kurtosis 1.917 → 1.673 | here (arithmetic [D]); m1_holder (sextic point only) | H/m_X 1.930099844e-21, refit −0.5676449 %, FOLLOWUP −1.18646 %/−1.7366 % all reproduce (`results/critic_checks.json` B); kurtosis initial 1.91708 reproduced from the npz | sextic +9.927 % charge NOT CHECKED (different operator, no engine); final kurtosis 1.673 NOT CHECKED |
| Finite-k controls, CLASS, 64/96-bin dependence, Codex's C-code conservation checks, portal-extension node (Xe 27.56 keV etc.), e = 1e-6 vector clocks, 6c504de/8f67a83/12b5696/4fe963f packets | nobody | — | NOT CHECKED; none is load-bearing for the owner's question except as quoted inputs |

## 2. Contradictions between tracks and skeptics, settled

1. **Holder late-charge drift (m1_holder 11/20 vs its skeptic).** The track attributed the ±0.15 % window dependence of
   the +1.806 % charge change to 3.08 light-Y quanta exchanging charge with X; the skeptic found the partner is Σ.
   Settled for the skeptic by Gauss-constraint arithmetic on the skeptic's own saved window means [C, data Q]:
   N_X + 3N_Y + 6N_Σ = 0 to 3e-8 in every window; between the 640–800 and 800–1000 windows ΔN_X = −10.286 (holder)
   while 6ΔN_Σ = +10.290 and 3ΔN_Y = −0.0039, i.e. Σ carries 100.04 % and Y 0.04 % of the change (original run:
   ΔN_X = +2.2394, 6ΔN_Σ = −2.2394, 3ΔN_Y = 5e-6). The Y excitation (×540) is real but is not the channel. The
   running-window ratio spans 1.01749–1.02082 (±0.17 %), so Codex's "+1.80645 %" and "H_i = 1.32691046e11 GeV" are
   meaningful to four digits; the exchange has not frozen at τ = 1000.
2. **EFT control of the Option 1A holder (m1_holder 3 "CONFIRMED for Y, margin 5" vs skeptic "not established").**
   Λ_eff = M_Y/ξ = 0.9757 H_i = 1.31e11 GeV [D]; the amplitudes are 5.9 (Y), 29 (Σ), 89 (X) Λ_eff and R/Λ_eff² = 12.6
   at R = 12H_i²; the coefficient relative to the declared 1e15 GeV cutoff is 7632; the holder term is 1.72 μ. The
   skeptic is right: no EFT margin exists relative to the operator's own scale; "margin 5" measured against M_Y is
   not the relevant statement. Also the boundedness proof needs 0.8 λ_X > 6.23e-4 (margin 1.6) and fails under any
   tree-level singlet completion with g₁ ≤ 4π (skeptic [E]); so f8186fa's "bounded below at fixed R" is CONFIRMED as
   algebra of the written polynomial and conditional as physics.
3. **Formation-clock scaling of the free-streaming length (m1_transport T11 vs its skeptic).** Own quadrature [E]
   (p ∝ 1/a from release at T_f with p_f = 0.684 m_X, crude g_*s table above 2 MeV, standard background below):
   λ_fs to equality = 388 pc (T_f = 2 MeV), 32.6 pc (20 MeV), 1.3 pc (200 MeV), 0.075 pc (2 GeV). Direction ∝ 1/T_f
   up to logarithms: the skeptic is right, the track's "λ_fs ∝ T_f (200 MeV → 366 pc)" is inverted. The track's
   conclusion that the scale is an input survives; an earlier clock makes the dip smaller and less observable.
4. **Core interior for a slow X (m1_transport T13 mechanism).** Own derivation of the rotating-frame Bogoliubov
   relation for the charged condensate, (k² − ν²)(k² − ν² + 2aF²) = 4Ω²ν² [D]: gapless branch, c_s² = aF²/(2Ω² + aF²)
   = 0.308 at the centre, and a free X with k = 0.0684 (ν = 0.1206 in the rotating frame) propagates inside with
   k_in = 0.216 m_X. The reviewer's "classically forbidden interior" (diagonal S² + 2aF² − 1 = +0.40) ignores the
   anomalous coupling; the skeptic's mechanism is right. The consequence is unchanged: both methods give a geometric
   cross section (hard sphere 160–265/m_X², coupled solve 332/m_X²), Γ/H ~ 1e6 at 2 MeV, momentum decoupling at a
   few keV, so the k ≳ 3e5/Mpc part of Codex's curve and its 0.62 MeV decoupling are not robust inside the model.
5. **M1V-10 "CORRECTED" (m1_vacuum vs skeptic).** The box-mode identification is right; Codex's DERIVATION.md already
   carries the finite-grid caveat but REPORT/HANDOFF do not. Settled as CORRECTED (wording of REPORT/HANDOFF).
6. **T9 "CORRECTED" on "testable" (m1_transport vs skeptic).** Codex wrote "a concrete conditional forecast that can
   now be challenged ... not yet an earned, distinctive GRUT prediction"; "testable" was the reviewer's word. Settled
   as CONFIRMED plus NEW (the observability analysis: k = 1e4–1e6/Mpc is M_hm = 5 to 5e-6 M_sun, 3.5 decades in k
   below any current probe; planned minihalo probes are O(1) in abundance).
7. **λ precision (m1_ledger vs skeptic).** Both agree on the arithmetic; the skeptic's band stands: the supplied
   D = 4.4e-10 GeV is 0.61 % above the Planck-2018 value 4.373e-10 GeV [R] (own recomputation from Ω_c h² = 0.1200,
   ρ_c/h² = 1.05368e-5 GeV cm⁻³, s₀ = 2891.2 cm⁻³), which alone moves λ by −0.40 %; clock choice ±0.2 %, history
   −0.09 %. Every legal reading lies in λ = 0.994–1.011 against "+0.7009 %" to eight digits.

## 3. Verdicts stronger than their evidence

- Codex: "fixed-Q radial Hessian positive (0.2233–0.2244)" → the eigenvalue is the free-X box level converging to
  1 − Ω² = 0.22257; the object's radial stability rests on VK (dQ/dω < 0), which does hold.
- Codex: "+0.7009 % by actual ODE rerun", λ = 1.0070090495735164 → arithmetic on a homogeneous charge estimator with
  0.44 % clock drift; the "rerun" cannot fail by construction (Y ∝ H_i^{3/2} exactly at fixed dimensionless history);
  meaningful to ±1 %.
- Codex: "+1.80645 %", "H_i = 1.32691046e11 GeV", "0.554384305" → four, four and five significant digits respectively
  (±0.17 % window drift; Codex's slope is 4.3e-6 high).
- Codex: "V bounded below at fixed R" → true for the written polynomial; fragile (λ_X margin 1.6) and not stable under
  the only kind of completion anyone has written down.
- Codex: "71/71" and "38/38 checks" → same-author controls, as Codex itself says; not independent.
- m1_holder 3 "CONFIRMED for Y" → not established (item 2 above). m1_holder 11/20 cause → refuted (item 1).
- m1_transport T11 clock direction → inverted (item 3); T13 mechanism → wrong premise, right size (item 4);
  T9 "CORRECTED" → CONFIRMED + NEW (item 6); T1's "not a contact cross-section" corrects the task text, not Codex.
- m1_vacuum M1V-10 "CORRECTED" → CORRECTED (wording) (item 5); roll time 1.9059 → 1.9079/m_X (0.1 %).
- m1_ledger: "(5,1): 1.016277" → 1.0162836 (slip); pulse rows 2–3 % high (charge-conserving η > 1; own
  recomputation 0.00433/0.01755/0.07209/0.3114, `results/critic_checks.json` G).
- Brief: "μ = 7.2 H_i" → 7.159456 H_i; "+3.42 % needed" → misreading; "newer isocurvature packet" → not on the remote.

## 4. Plain answer to "how is the theory coming together"

**What M1 establishes (as checked by this review, within its stated limits: classical, tree level, homogeneous or
spherical, radial stability only).**
- The tilted potential has one vacuum, V ≥ 0, no static walls; the hidden tree potential is bounded with slack. Exact.
- The Q = 200000 object exists with E/(Q m_X) = 0.93740 (binding 6.26 %), Ω = 0.881720, mass 17324.83 GeV; it is
  radially stable in the VK sense; objects are bound for Q > 7.13e4 and exist for Q > 5.67e4; E/Q → 0.77586 slowly
  (0.78 only near Q ~ 1e10). Confirmed to 1e-13 by two more solvers.
- A dilute X condensate below n = 233 m_X³ has an imaginary sound speed: it fragments (standard Q-ball
  fragmentation; Kusenko–Shaposhnikov, Enqvist–McDonald, Kasuya–Kawasaki). This is the clumping mechanism; it does
  not yet give a Q distribution.
- The calibration arithmetic is self-consistent: Y_Q = D/ε_Q = 5.079e-9, H_i = 1.3428e11 GeV × (1 ± 0.01).
- The Option 1A holder has a stable prepared point at R = 12 H_i² with the wanted phase; it raises the late charge by
  1.8 % (±0.2 %) and costs 1.2 % in H_i. The heavy spectator sector is adiabatic at k → 0 (Weinberg's theorem) and
  its intrinsic charge power decays as e^{−3N}; the "no resolved isocurvature signal" is the generic heavy-spectator
  expectation, not a model feature, and not P_CDI.
- Given a 6.6 % free-X gas at v = 0.068 c at 2 MeV, its XX self-scattering (σ = 2.17e-7/m_X², contact + φ exchange),
  7.9 collisions per particle, 0.62 MeV viscosity crossing and 35 pc collisionless free streaming are right.
- The cold vacuum-Born charge loss Γt₀ = 6.27e-17 is right for that channel.

**What M1 assumes (every supplied input found in the packets).**
Light sector: g_X = 0.1, λ_X = 1e-3, λ_φ = 0.03, v = 0.92409 GeV (m_X = 92.41 MeV), the new tilt **ε = 0.002**, the
reference charge **Q = 200000**. Heavy sector: M_Y/H_i = 30, h/H_i = 40, μ/H_i = 7.159, κ = 3.4587e-3, λ_Y = λ_Σ = 1,
e = 1e-6, ξ_X = −5/6, cutoff 1e15 GeV, zero SM portals, **ξ_H = 0.1** (stabilising choice, survival not computed).
Phase holder: **ξ = 30.7467 and α = 2.1921** (one complex coefficient fixed by two conditions at R = 12 H_i²; α is
the supplied CP-odd sign; the operator's own scale is ~H_i); the sextic ξ₆ = 4e8, Λ = 1e15 GeV is superseded.
Reheating: quadratic inflaton with m_I/H = 3, Γ_I/H = 0.5 (alternative 5, 1), a supplied prepared initial state
(amplitudes, zero velocities, arg Σ − 6 arg X = −π/2), g* = 106.75, homogeneous; H_i itself is the calibrated output.
Targets: **D = ρ_DM/s = 4.4e-10 GeV** (rounded; Planck 4.373e-10 [R]), the vacuum energy ρ_vac = 2.52e-47 GeV⁴
(listed as an input, underived), the matching hypothesis η = 1, p = 0, and the **1 % tolerance** (a program
benchmark, as Codex says, not a confidence interval). Formation and gas: **one 10 % radial Gaussian pulse on one
prepared core** (5/20/40 % as controls), the **20 MeV formation clock** with n_core = 1.21e-15 m_X³ and a static
background for τ ≤ 1000, the inherited 96-bin checkpoint, the 2 MeV BGK matching (equal density and velocity
perturbations, an assumed temperature perturbation, everything else cold, radiation one fluid), gas = 5.877 % of
matter. Also the fixed dimensionless history (Y ∝ H_i^{3/2}) and the spectator assumption (true to 1e-11 here).

**Is the parsec-scale curve a distinctive, testable prediction?** No, on both counts, and Codex says the first
itself. Amplitude: the plateau (1 − f)² = 0.886 is set by f, and f = 0.845 % × 7.37 comes from the supplied pulse
amplitude (f ∝ A²) times an unreviewed catalytic chain at the supplied clock. Scale: λ_fs ∝ 1/T_f, set by the
supplied 20 MeV clock. Shape: the standard mixed cold + warm curve (Boyarsky et al. 2009); a two-parameter mixed-DM
fit absorbs both inputs. Observability: the dip lives at k = 1e4–1e6/Mpc, M_hm = 5 to 5e-6 M_sun, 3.5 decades in k
below the smallest scale reached by Lyman-α, satellites, lensing or streams; planned minihalo probes (PTA, caustic
microlensing, 21-cm) are order-unity in abundance, not 10 %. Inside the model the k ≳ 3e5/Mpc part is not robust:
gas–core scattering (Γ/H ~ 1e6, omitted by Codex) locks the gas to the cores until T ~ 2–3 keV and turns the
free-streaming cutoff into a phase-sensitive acoustic feature. Nothing at Lyman-α scales is affected (|ΔP/P| ≲ 1e-8).

**What the frozen mass test needs, against what has been computed.** The test 0.0668(1 − η) + 2.134 p ≤ 0.01
(η = charge fraction in Q = 2e5 objects, p = surviving pairs per charge) is confirmed as algebra. No formation-type
population computed by anyone satisfies it: pt_relax equilibrium X̄ inventories 0.07–0.11 per charge give 0.15–0.23
(×15–23); the coherent-orbit inventories under the M1 tilt (pt_m1tilt, slow limit: 1.49–1.52 X̄ + 0.16–0.18 φ per
charge for the τ-20 orbit, 0.11–0.15 + 0.22–0.25 for the late orbit; the M1 tilt removes the chaotic split but the
slow limit is still lossy) give p ≈ 1.67 and 0.35 before annihilation, i.e. ×360 and ×75; the 3D mask fractions read
as η (0.31, 0.49) give ×4.6 and ×3.4 (Codex rightly says masks are not η); the pulses give 0.0043/0.0175/0.072/0.311
with charge-conserving η (5 % passes, 10 % fails ×1.75). The passes are conditional: kin_solve gas-only annihilation
at total heat R ≲ 5.8 (the e = 1e-6 band R ≈ 4.97 passes narrowly, 0.0050), the 10 % pulse after Codex's idealised
τ = 1000 kinetics (0.0039 frozen, 0.0073 with the re-solved E(Q), and 1.35 % of the mass still in gas heat), and the
cold completions with every X̄ converted. What is needed: **η** from a formation ensemble (an expanding evolution from
the matched coherent release under the tilted potential, returning the charge-weighted object spectrum y_Q with
membership by E/Q < m_X, not by density masks), and **p** from the X̄ inventory after the φ transition (pt_m1tilt)
followed by annihilation/conversion/capture at the actual heat R(e) (kin_solve covers the gas-only limb; c481530 one
fixed-profile pulse). One more point the tracks make only implicitly: a failure of the 1 % test is not a
falsification, because the dimensionless history is unchanged and H_i refits by [1 + 0.0668(1 − η) + 2.134p]^{−2/3}
(0.91 at p = 0.07, 0.69 at p = 0.35, 0.36 at p = 1.67). What a large p would change is the composition: symmetric
pairs would be 13 %, 43 % and 78 % of the dark matter in those three cases, and the "asymmetric charged objects"
picture would be gone unless annihilation at the actual heat removes them.

## 5. The single most valuable next calculation and its failure condition

**A formation ensemble under the M1 potential**: evolve the matched coherent release (the pt_m1tilt/pt_matching
state with the deformed cubic, expanding background, the actual density clock) in three dimensions far enough past
fragmentation to return (i) the charge-weighted object spectrum y_Q and η = Σ_{Q ≥ 7.13e4} Q y_Q / Q_net with
membership by E/Q < m_X, (ii) the X̄ inventory and its fate through annihilation/conversion/capture at the actual
heat R(e), hence p, and (iii) the free-X fraction f and velocity distribution that replace the supplied pulse and
clock. Everything downstream (ε_Q, the calibration, the gas, the transfer curve) is a function of these three.
**Failure condition** (Codex's own, made explicit): with the re-solved E(Q)/Q, 0.0668(1 − η) + 2.134 p > 0.01 means
the M1 calibration at (Q = 2e5, η = 1, p = 0) is wrong and the composition is not what M1 claims; the hard version is
either no objects above Q = 7.13e4 at all (nothing bound forms) or p ≳ 0.1 surviving (dark matter mostly symmetric
pairs). This is a 3D run and outside this review's budget; two cheap ODE precursors are worth doing first: follow the
X ↔ Σ exchange past τ = 1000 to freeze-out (it fixes the "late charge" input to better than 0.2 %), and solve the
X–core three-channel elastic problem at k ~ 0.07 m_X (it fixes the gas decoupling temperature).

## 6. Consolidated verdict table (item · verdict · ours · Codex · method)

| # | Item | Verdict | Ours | Codex | Method |
|---|---|---|---|---|---|
| 1 | One-vacuum potential: discriminant, threshold, V ≥ 0, no walls | CONFIRMED | −9/62500 exact; ε > 0.001716; V(−v) = 160 m_X⁴ | −0.000144; 0.001715729 | exact algebra, two passes (m1_vacuum) |
| 2 | Hidden tree bound λ_Y,eff = 0.758516, margin 2445.04 | CONFIRMED | same; brute force slack 1.03–3.5 | same | algebra + minimisation |
| 3 | Bulk Ω = 0.77586471, n = 449.31 m_X³ | CONFIRMED | 30-digit identical; spinodal 233.2 | same | algebra |
| 4 | Q = 2e5 object: Ω, E/Qm, mass | CONFIRMED | 0.8817204919524, 0.9373999905, 17324.832346 GeV (two solvers, 1e-13) | same | own Numerov and cell-centred BVP |
| 5 | Fixed-Q radial Hessian "positive 0.2233–0.2244" | CORRECTED (wording) | box mode → 1 − Ω² = 0.22257; VK dQ/dω < 0 holds; Q_min 5.67e4, bound for Q > 7.13e4 | 0.2243/0.2233 | own Hessian, secular test |
| 6 | Cold-fluid growth F = 1: Γ 0.11428 at k 0.42 | CONFIRMED | 0.114280 at 0.41923 = −c_s² below the spinodal | same | own 6×6 / cubic |
| 7 | Frozen-test coefficients 0.0667805, 2.13356 | CONFIRMED | same to 2e-15 | same | algebra |
| 8 | Calibration λ = 1.00700905, Y_Q = 5.07941423e-9 | CONFIRMED (arithmetic); precision CORRECTED | λ band 0.994–1.011 (clock, history, D rounding) | 8 digits | own engine τ = 800; 50-digit arithmetic |
| 9 | Y_Q ∝ H_i^{3/2} "verified by ODE rerun" | CONFIRMED (by construction) | 1.5000000002 | as stated | own engine; Codex control rows |
| 10 | Γt₀ = 6.27e-17 (from 7.08e-17) | CONFIRMED | 6.27163e-17; I = 2868466.51 (3 quadratures) | 6.2716e-17 | own Born derivation |
| 11 | Holder (ξ, α), stationary point, Hessian, bound −285370 | CONFIRMED (bound fragile) | identical to 1e-14; λ_X margin 1.6; completion with g₁ ≤ 4π breaks B > 0 | as stated | mpmath / sympy |
| 12 | EFT control of the holder | NEW (not established) | Λ_eff = M_Y/ξ = 0.976 H_i; amplitudes 6–89 Λ_eff | "dimension 5, not UV-complete" | algebra (here, skeptic) |
| 13 | Holder charge +1.80645 % | CONFIRMED to 4 digits | 1.0180646 (640–800); 1.0175–1.0208 running windows | 1.01806453 | three ODE codes |
| 14 | Cause of the late drift | SETTLED: X ↔ Σ exchange | 6ΔN_Σ = −ΔN_X to 4e-4; Y share 4e-4 | not reported | Gauss arithmetic (here) |
| 15 | H_i refit 1.32691046e11 GeV | CONFIRMED to 4 digits | 1.3269 ± 0.0015e11 | 9 digits | arithmetic |
| 16 | Phase slope 0.554384305 | CORRECTED | 0.554382 (two codes) | 0.5543843 | own finite differences |
| 17 | Time shift = adiabatic mode; residual 4e-12 / 3.02 | theorem CONFIRMED; numbers NOT CHECKED | −9e-7 = drift rate; scalar-only −0.48 | packet absent | own restarts; Weinberg 2003 |
| 18 | Spectator power 1.236e-10 … 1.443e-70 | CONFIRMED | 1.4e-9 with Codex's row; e^{−3N} to 1e-4 | as stated | exact Hankel, x-ODE |
| 19 | Hidden-scalar power 3.435e-70 | NOT CHECKED | 1.443e-70 under the f8186fa definition | 3.435e-70 (absent) | — |
| 20 | σ(XX) = 2.16627561e-7/m_X²; Γ(φ → XX̄) = 1.875659e-4 m_X | CONFIRMED | same (contact + φ exchange; exchange 3.3× contact) | same | own Feynman rules ×2 |
| 21 | Γ_coll/H 15.72; 7.91 collisions; 0.62156 MeV; 35.51 pc | CONFIRMED | 16.0 / 7.90 / 0.6227 MeV / 36.6 pc (Maxwell); repro 15.81 / 7.81 / — / 35.50 | as stated | own rates, background; repro |
| 22 | Kurtosis 1.91713 → 1.67301 | initial CONFIRMED (repro); final NOT CHECKED | 1.91708 from the npz | as stated | reproduction (here) |
| 23 | BGK curve 0.9979 … 0.7642 | NOT CHECKED (magnitudes agree) | fluid toy 0.88–0.91 (1e5), 0.78–0.80 (1e6) | as stated | two fluid toys |
| 24 | "Conditional forecast, not yet distinctive" | CONFIRMED + NEW | unobservable: M_hm 5–5e-6 M_sun, 3.5 decades below probes; f ∝ A², λ_fs ∝ 1/T_f | as stated | halo-mass conversion; literature [R] |
| 25 | Gas source 6.595 %, v_rms 0.0684 c | CORRECTED (interpretation); root NOT CHECKED | 0.845 % (supplied pulse) × 7.37 (unreviewed chain) at a supplied clock | 6.5946 % | ledger arithmetic |
| 26 | Gas–core elastic scattering | NEW | σ 160–332/m_X², Γ/H ~ 1e6, decoupling 2–3 keV; interior propagating (k_in 0.22) | omitted | hard sphere, coupled solve, dispersion (here) |
| 27 | 94f660d: H/m_X 1.9301e-21; refit −0.567645 %; FOLLOWUP −1.186 %/−1.74 % | CONFIRMED | 1.930100e-21; −0.567645 %; −1.18646 %/−1.7366 % | same | arithmetic (here) |
| 28 | Sextic holder δ = −2.09993, lightest 10.669; +9.927 % charge | point CONFIRMED; charge NOT CHECKED | −2.099927, 10.66925 | as stated | own point solve |
| 29 | ed73791 3D: object survives (E/Q 0.943), masks 30–52 %, pulse emission 0.002–0.095 | NOT CHECKED (3D) | — | as stated | — |
| 30 | c481530 kinetics chain 1695 X̄ → 0.559; free X ×7.37; Q_core/Q₀ 0.9419 | NOT CHECKED | ledger arithmetic only | as stated | — |
| 31 | Frozen test vs computed populations | NEW | no formation-type population passes (×1.75–360); passes conditional | "η, p remain to compute" | ledger + pt_m1tilt + kin_solve |

## 7. Files and reproduction

```
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd research/reality_framework/grut_formation_review_20261025/tracks/m1_critic/code
python3 -I critic_checks.py        # -> ../results/critic_checks.json (0.2 s; item C reads the holder skeptic's JSON,
                                   #    item H reads the m1_transport scratch copy of Codex's npz if present)
```
`results/checkpoint_kurtosis_repro.json` is a labelled reproduction from Codex's 6dc9556 npz; `results/verdicts.json`
is the table above in machine-readable form. Scratch: /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/m1_critic/.
Credit: Weinberg 2003 and Wands–Malik–Lyth–Liddle 2000 (adiabatic mode); Vakhitov–Kolokolov / Lee–Pang / Friedberg–Lee–Sirlin
(Q-ball stability); Kusenko–Shaposhnikov, Enqvist–McDonald, Kasuya–Kawasaki (fragmentation); Boyarsky et al. 2009 and the
mixed/warm dark matter literature cited in m1_transport (observability, [R]).
