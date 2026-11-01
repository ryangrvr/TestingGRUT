# Review of Codex's GRUT-M1 chain (d4747cb through f8186fa)

1 November 2026. Claude, issue #13. Brief: [BRIEF_M1.md](BRIEF_M1.md). Four tracks, each with an adversarial skeptic,
then a completeness critic: `tracks/m1_vacuum`, `m1_ledger`, `m1_holder`, `m1_transport`, `m1_critic`. The transition
check under the M1 potential (`tracks/pt_m1tilt`) is included where relevant and is **preliminary** until its skeptic
reports. Each track's NOTES.md and SKEPTIC.md give methods, commands and full numbers.

Claude was unavailable 27–31 October; the later M1 packets were handed to a third AI ("Argon"/Gemini) and had no
independent review until now. This is same-session AI work on one supplied model. It is not outside replication.
Nothing here is observational evidence, a validation, an exclusion or a GRUT prediction. All data comparisons are
retrospective.

## Result

**Every Codex number we checked reproduces**, most of them with our own code to machine precision: the one-vacuum
algebra, the Q = 200000 object, the bulk state, the cold-fluid dispersion, the frozen mass test, the calibration chain
and its exponent, the charge-loss rate, the holder's stationary point, Hessian, bound and charge change, the adiabatic
cancellation and the de Sitter power table, the gas cross sections, rates and free-streaming distances. Codex's
bounded request (one-vacuum argument, soliton and Hessian, mass normalisation, charge-loss normalisation) is fully
covered and confirmed.

Six interpretations change:
1. The soliton's "positive fixed-charge radial Hessian" is the lowest free-X box mode; what holds is the
   Vakhitov–Kolokolov criterion, and it holds for every Q > 5.67×10⁴ (new).
2. The calibration λ = 1.00700905 is right as arithmetic, but every legal reading of the same calculation spans
   0.994–1.011, and the target D = 4.4×10⁻¹⁰ GeV is a rounded supplied value 0.6% above the Planck value.
3. The holder's late charge has not frozen: it exchanges slowly with Σ, so +1.806% is good to about ±0.17% and
   H_i = 1.3269 ± 0.0015×10¹¹ GeV. Its EFT control is not established: the operator's own scale is M_Y/ξ ≈ H_i.
4. The parsec-scale power dip is unobservable with any current or planned probe, has the generic mixed cold+warm
   shape, and its amplitude and scale are inputs (a supplied pulse and a supplied clock). Its sub-30 pc part is not
   robust even inside the model, because gas–core scattering, which Codex omitted, is 10⁵ times faster than gas–gas.
5. The adiabatic cancellation is the generic heavy-spectator result, not an isocurvature bound.
6. No population computed so far passes the frozen mass test; the passes are conditional recipes.

## Verdicts

| Item | Codex | Ours | Verdict |
|---|---|---|---|
| One vacuum: V′ ∝ (u−1)[λu² + (λ−3ε)u + 3ε], discriminant −0.000144; V(v) = V′(v) = 0, V″ unchanged | as stated | exact: −9/62500; V/v⁴ = (u−1)²[λ(u+1)²/4 − ε(u−1)] ≥ 0, zero only at +v; one-vacuum window ε ∈ (0.0017157, 0.0583); V‴(v) down 6.67% | CONFIRMED |
| Hidden tree bound: λ_Y,eff = 0.758516, large-Σ margin 2445.04 | as stated | identical; brute-force minimisation of the actual potential shows slack (×1.03 at s_c, ×3.5 at 100 H_i). The real caveat is loop running (the Landau pole near 6×10¹⁵ GeV at λ_Y = λ_Σ = 1) | CONFIRMED |
| Bulk state S = 0.15117, F = 2.40648, Ω = 0.775865, n = 449.31 m_X³; no second minimum at fixed charge density | as stated | identical to 12 digits; S(F) unique in (0, 1) for every F; fixed-n Hessian positive everywhere; the only instability is thermodynamic (dΩ/dn < 0 below the spinodal n = 233.2 m_X³). The tilt moves the interior from φ ≈ 0 to S_c = 0.1258 and Ω_bulk from 0.7401 to 0.7759 | CONFIRMED |
| Q = 200000 object: Ω = 0.88172049, E/(Qm_X) = 0.93739999, S_c = 0.12577, mass 17324.83235 GeV; dE/dQ = ω to 2.9×10⁻⁷; dΩ/dQ < 0 | as stated | two own BVP solvers: Ω = 0.8817204919524, E/Qm = 0.9373999905262, mass 17324.832346 GeV (10⁻¹³); the 2.9×10⁻⁷ is Codex's finite-difference step; dΩ/dQ = −2.0500×10⁻⁷ | CONFIRMED |
| Fixed-Q radial Hessian positive (0.2233–0.2244) | as stated | reproduced exactly, but the lowest eigenvalue runs 0.2308 → 0.2227 with box size toward the continuum edge 1 − Ω² = 0.22257 as (π/L)²: it is the lowest free-X box mode. L₊ has exactly one negative eigenvalue (−0.095), lifted by the charge term. Correct statement: no discrete radial fixed-charge mode below threshold; Vakhitov–Kolokolov dQ/dΩ = −4.9×10⁶ < 0 | CORRECTED (wording) |
| Object family | Q = 200000 is "an assumed reference" | E(Q)/Q = 0.9790 (10⁵), 0.9374 (2×10⁵), 0.8928 (5×10⁵), 0.8670 (10⁶), 0.8161 (10⁷); binding E/Q = m_X at Q = 7.13×10⁴; smallest object Q_min = 5.67×10⁴ (Ω = 0.9727); every Q > Q_min is VK-stable; the bulk 0.7759 is reached only near Q ~ 10¹⁰ (surface tension 99 m_X³) | NEW |
| Cold-fluid growth for dilute densities: Γ = 0.11428 m at k = 0.42 (F = 1) | as stated | 0.114280 at 0.41923; band k < 0.615; F = 0.5: 0.0292; F = 1.5: 0.1838; bulk and core densities stable. At k → 0 it is the imaginary sound speed of a condensate below the spinodal: standard fragmentation into objects (Kusenko–Shaposhnikov; Enqvist–McDonald; Kasuya–Kawasaki), not an instability of the objects | CONFIRMED |
| Walls removed | "removes the two-vacuum φ sign-wall mechanism" | V(−v) = 8εv⁴ = 160 m_X⁴, V(0) = 95; V monotone below v, so no static wall exists; negative-φ patches cost ≥ 160 m_X⁴ per volume and roll back in 1.9/m_X. Note: Codex's deformation grid straddles the threshold (ε = 0.001 still has a metastable minimum at u = −0.77) | CONFIRMED |
| Frozen test 0.0667805(1−η) + 2.13356p ≤ 0.01; η ≥ 0.85026 or p ≤ 0.004687 | as stated | coefficients (m_X/ε_Q − 1) and 2m_X/ε_Q with ε_Q = 0.0866242 GeV: 0.0667805, 2.133561; Y_Q ε_Q = 4.4000×10⁻¹⁰ GeV. The "old" Y_Q = 4.8294×10⁻⁹ was never an object calibration: D/Y_Q,old = 0.98593 m_X = the old object 0.923316 m_X + 2 × 0.031306 m_X of scalar-release pairs (the origin of kin_solve's 2.0286) | CONFIRMED |
| Calibration: +0.700905% coordinated scale, Y_Q ∝ H_i^{3/2}, Y_Q = 5.07941423×10⁻⁹ | as stated, eight digits | the heavy point on (3, 0.5) delivered 1.040807 × 4.8294×10⁻⁹ = 5.0265×10⁻⁹ (our engine at τ = 800: 5.02648×10⁻⁹); λ = (5.07941/5.02648)^{2/3} = 1.00700905; exponent measured 1.5000000002 (Codex's own factor-4 controls already fix it; the "ODE rerun" cannot fail by construction). Our earlier 1.0454 was the τ = 80 window. **Precision:** Codex's own charge clocks give λ = 1.0048–1.0083; the (5, 1) history gives 1.00608; and D = 4.4×10⁻¹⁰ GeV is a rounded supplied value, 0.61% above Planck 2018's 4.373×10⁻¹⁰ [R] with a ±1% band. Every legal reading lies in λ = 0.994–1.011 | CONFIRMED; precision CORRECTED |
| Charge loss Γt₀ = 6.27×10⁻¹⁷ (was 7.08×10⁻¹⁷) | as stated | 6.27163×10⁻¹⁷ = old × 0.9105 (profile) × 0.9725 (A ∝ λ⁻²); three Born quadratures agree to 10⁻⁹ | CONFIRMED |
| Option 1A holder: dimension 5; ξ = 30.74672874, α = 2.19214407; μ_eff(12H²) = 10i H_i | as stated | ξ and α reconstructed exactly (μ = 7.159456 H_i); the pair is a tuned cancellation at R = 12 H_i² exactly; dδ_eq/dlnR = 0.7159456. α is a rephasing-invariant CP phase: the supplied sign our preparation theorem requires. Codex states no selection law is claimed | CONFIRMED |
| Stationary point (87.19556, 5.75425, 28.38718) H_i, phase −π/2, Hessian (13.0914, 14.5892, 538.147, 1546.66, 3239.89), bound −285370 H_i⁴ | as stated | three own codes agree to 10⁻¹⁴; 600-start global search finds no other minimum; the bound chain is exact but its λ_X margin is only 1.6 | CONFIRMED |
| EFT control of the holder | "dimension-5 … not UV-complete" | **not established:** the operator's own scale is Λ_eff = M_Y/ξ = 0.976 H_i, below every field amplitude (\|Y\|, \|Σ\|, \|X\| = 5.9, 29, 89 Λ_eff) and below √R (R/Λ_eff² = 12.6); relative to the declared 10¹⁵ GeV cutoff the coefficient is 7632; a loop would give ξ ~ 10⁻⁴–10⁻³. Any tree-level singlet completion (g₁ M_S S R + g₂ S Σ*Y²) induces −(1693/g₁²) λ_Y\|Y\|⁴ at the prepared point, so the fixed-R boundedness proof needs g₁ > 57 and fails for any g₁ ≤ 4π [E] | NEW |
| Holder charge +1.80645% → H_i = 1.32691046×10¹¹ GeV | as stated, eight digits | two independent codes: ratio 1.0180646 (Codex's digits beyond the 7th are not reproducible). The late N_X has **not frozen**: 80-wide windows over H_i t = 400–1000 give 1.0175–1.0208 (±0.17%). The partner is Σ, through the Y-mediated sextic κ²uΣ*X⁶/M_Y⁴ (6ΔN_Σ = −ΔN_X to 4×10⁻⁴), kicked once per late X orbit (period ~300–360 H_i t, decaying as a⁻²). H_i refit: 1.3269 ± 0.0015×10¹¹ GeV | CONFIRMED (4 digits); drift NEW |
| Phase slope 0.554384305 | as stated | 0.5543820 (two codes; Codex's is 4×10⁻⁶ high); d4747cb's −0.21997 is the re-prepared conditional minimum (also reproduced); the two definitions differ | CORRECTED (digits) |
| Adiabatic single-clock cancellation; spectator charge power 1.236×10⁻¹⁰ … 1.443×10⁻⁷⁰ | as stated | theorem written out (Weinberg 2003; Wands–Malik–Lyth–Liddle 2000): δy = εẏ solves the tangent system including the −V_rR δR forcing; our nonlinear restarts give −9×10⁻⁷ (= the late drift rate); a scalar-only shift gives −0.48 at t₀ = 3 and is cancelled by the background shift. Power table reproduced to 1.4×10⁻⁹ by exact Hankel modes; ratios e^{−3ΔN} to 10⁻⁴. "No resolved correlated signal" is the generic expectation for a spectator with m² ≥ 13 H² > 9H²/4, not an isocurvature bound | CONFIRMED |
| Newer isocurvature packet (residual 4×10⁻¹², 3.02, charge power 3.435×10⁻⁷⁰, 38 checks) | reported in chat | **not on the remote** (f8186fa is the newest M1 head on 1 November). 3.435×10⁻⁷⁰ does not follow from f8186fa's response row under any convention tested | NOT CHECKED |
| σ(XX → XX) = 2.16627561×10⁻⁷/m_X²; Γ(φ → XX̄) = 1.875659×10⁻⁴ m_X | as stated | both exact; σ is contact (−4λ_X) plus t/u-channel φ exchange (+0.01333), as Codex's derivation states; the exchange is 3.3× the contact amplitude | CONFIRMED |
| Γ_coll/H = 15.7 at 2 MeV; 7.91 collisions over a factor 10; viscosity proxy crosses H at 0.62 MeV; Maxwellisation | as stated | 16.0 (Maxwellian) / 15.8 (Codex's shape); 7.90; 0.623 MeV (the total rate crosses at 0.51 MeV); initial kurtosis 1.91708 reproduced | CONFIRMED |
| Free streaming to equality 35.5 / 39.9 / 34.1 pc; z_eq = 3420.5 | as stated | 35.50 / 39.86 / 34.06 with Codex's distribution in our background; 36.6 / 39.7 for a Maxwellian; z_eq = 3417 | CONFIRMED |
| BGK transfer: P_mix/P_cold = 0.893 at 63 pc, 0.764 at 6.3 pc | as stated | hierarchy NOT CHECKED; our two-fluid toy gives 0.88–0.91 at k = 10⁵/Mpc: the generic cold+warm plateau (1−f)² = 0.886 plus a growth deficit | NOT CHECKED / consistent |
| "A concrete conditional forecast that can now be challenged" | as stated | k = 10⁵/Mpc is M_hm ≈ 5×10⁻³ M_⊙ (k = 10⁴: 5 M_⊙). Current probes stop at k ≈ 25–50/Mpc (M_hm 4×10⁷–3×10⁸ M_⊙: Lyman-α arXiv:1702.01764; satellites 2008.00022; lensing 1908.06983; streams 1911.02663); planned minihalo probes (PTA 2005.03030, 2012.09857; microlensing 1908.01773; 21-cm 1911.11144) are O(1) in abundance, not 10%. The shape is the standard mixed cold+warm one (Boyarsky et al. 0812.0010). Lyman-α puts no constraint on a 6–7% fraction with v ≈ 0.07c at 2 MeV (\|ΔP/P\| ≤ 10⁻⁸ at 10–30 h/Mpc) | NEW: unobservable, not distinctive |
| Gas source: 6.595% free X, v_rms = 0.068c | as stated | = 0.845% direct X̄ emission from one supplied 10% pulse on one core, times a ×7.37 catalytic chain (c481530), at a supplied 20 MeV clock; emission ∝ A^{2.0–2.1} over the 5–40% pulses; p_rms(20 MeV) = 0.68 m_X from the pulse's 2.3% energy excess. Both the amplitude (f ∝ A²) and the scale (λ_fs ∝ 1/T_f: 2 MeV → ≈400 pc, 200 MeV → ≈5 pc) are inputs. The 3D pulse runs and the kinetics chain are the uncertified root | CONFIRMED as Codex states; interpretation NEW |
| Gas–core elastic scattering | omitted | coupled three-channel solve: σ ≈ 332/m_X² at k = 0.068 m_X (geometric; a resonance 1333/m_X² near k = 0.15), Γ/H ≈ 1.4×10⁶ at 2 MeV, against 16 for gas–gas; momentum decoupling at T ≈ 2–3 keV, not 0.62 MeV. The core interior is a propagating Bogoliubov medium, not a hard wall, but the cross section is geometric anyway. Codex's own X̄ results implied σ_tot ≈ 210/m_X². In a two-stage fluid toy the k ≥ 3×10⁵/Mpc part of the curve is not robust (0.3–0.7 at 3×10⁵; 0.006–0.6 at 10⁶, phase sensitive) | NEW |
| 94f660d: H/m_X = 1.930×10⁻²¹, refit −0.5676%; sextic holder δ = −2.0999, +9.93% charge | as stated | 1.930100×10⁻²¹, −0.567645%, FOLLOWUP's −1.18646% / −1.7366% reproduce; sextic point and δ reproduce; the +9.93% charge NOT CHECKED | CONFIRMED |

Not checked anywhere: the ed73791 3D runs (object survival, masks, pulse emission 0.002–0.095 Q), the c481530
kinetics chain, the BGK code, the final kurtosis, the CLASS attempts, and the 38-check packet.

## The frozen mass test against what has been computed

With η the charge fraction in reference objects and p the leftover pairs per charge **after** annihilation, the test
is 0.0667805(1−η) + 2.13356p ≤ 0.01. No formation-type population computed so far passes it; the passes are
conditional recipes:

| Scenario | Source | Test term | Margin |
|---|---|---:|---:|
| Equilibrium X̄ at the pilot's own energy, 0.069–0.110 per charge, before annihilation | REVIEW_TRANSPORT, pt_relax | 0.147–0.235 | ×15–23 |
| Same with the e = 10⁻⁷ vector heat (2.3–2.9) | pt_relax | 4.9–6.2 | ×490–620 |
| The 1.75 plateau | 5f5513b | 3.73 | ×370 |
| Coherent-orbit inventory after the M1 transition, slow limit: 1.67 (τ-20) / 0.35 (late) pairs per charge (preliminary) | pt_m1tilt | 3.56 / 0.75 | ×356 / ×75 |
| Cold dilute masks 30–32% of charge, read as η | ed73791 | 0.045–0.047 | ×4.5 |
| Hot elliptic control 0.6% | ed73791 | 0.066 | ×6.6 |
| Pulses 5/10/20/40%, charge-conserving η, emitted X̄ and φ as pairs | bf8ab3b | 0.0043 / 0.0176 / 0.072 / 0.311 | 5% passes; others ×1.8–31 |
| 10% pulse kinetics (η = 0.942, p = 3×10⁻⁶) | c481530 | 0.0039 frozen / 0.0073 re-solved E(Q), plus 1.35% of the reference mass still in gas kinetic energy | conditional |
| Gas-only annihilation at total heat R = 5 / 7 / 41 (e = 10⁻⁷) | kin_solve | 0.0050 / 0.030 / 1.57 | R = 5 passes (leaves η ≥ 0.923); crossing R ≈ 5.4–5.8 |

Failure of the test is not a falsification: H_i refits by [1 + 0.0668(1−η) + 2.134p]^{−2/3}, but pairs would then be
13–78% of the dark matter. What the test needs is a formation ensemble for η (masks are not η; membership by
E/Q < m_X) and, for p, the X̄ inventory after the φ transition followed by annihilation at the actual heat R(e).

## What M1 establishes and what it assumes

Established, at tree level and classically: one vacuum and no walls; a bound, VK-stable Q = 200000 object with
E/Q = 0.9374 m_X; objects bound above Q = 7.13×10⁴ and existing above 5.67×10⁴; a dilute condensate that fragments
(negative compressibility below 233 m_X³); a self-consistent calibration H_i = 1.343×10¹¹ GeV (±1%); a locally stable
phase holder that costs 1.2% in H_i; heavy-spectator adiabaticity; correct transport numbers for a given gas.

Supplied, every one: g_X = 0.1, λ_X = 10⁻³, λ_φ = 0.03, v = 0.924 GeV, **ε = 0.002**, **Q = 200000**; M_Y/H = 30,
h/H = 40, μ/H = 7.16, κ = 3.46×10⁻³, λ_Y = λ_Σ = 1, e = 10⁻⁶, ξ_X = −5/6, a 10¹⁵ GeV cutoff, zero portals,
**ξ_H = 0.1**; **ξ = 30.75 and α = 2.192** (the CP sign, in an operator whose own scale is ~H_i); the inflaton
(3, 0.5) with a supplied prepared state and g* = 106.75; **D = 4.4×10⁻¹⁰ GeV** (rounded), ρ_vac = 2.52×10⁻⁴⁷ GeV⁴
(underived), η = 1 and p = 0, the **1% tolerance** (a program benchmark); **one 10% pulse on one core at a 20 MeV
clock**, a static kinetics window, and the 2 MeV matching of the transport hierarchy.

The parsec-scale curve is neither distinctive nor testable: amplitude from the pulse, scale from the clock, generic
mixed shape, 3.5 decades in k below any probe, and not robust below 30 pc inside the model.

## The transition under the M1 potential (preliminary; `tracks/pt_m1tilt`)

Under V_φ with the ε = 0.002 tilt, the φ-condensation transition of the coherent X orbit is deterministic (the untilted
chaotic split disappears) and adiabatic for ramp rates ε_rate ≥ 10⁻⁴ (τ-20 orbit) / 10⁻³ (late), returning J_r to
2×10⁻⁴. Below a sharp threshold it deterministically loses radial action, growing as the rate slows: at
ε_rate = 10⁻⁵–10⁻⁶ the τ-20 orbit gives 1.49–1.52 X̄ (0.85 J_r) plus 0.16–0.18 φ quanta per charge; the late orbit
0.11–0.15 (0.23–0.30 J_r) plus 0.22–0.25. A tilt-driven φ oscillation slaved to the radial frequency releases a free φ
component where that regime ends. The physical rate is ~10⁻¹², ten decades beyond the scan [E]. The adversarial check
is still running; this paragraph will be updated.

## Most valuable next calculation

The formation ensemble under M1: expanding 3D evolution from the matched release with the deformed cubic, returning
y_Q, η (membership by E/Q < m_X), the velocity distribution, and p after annihilation at R(e). Failure condition:
0.0668(1−η) + 2.134p > 0.01 with the re-solved E(Q), or no objects above Q = 7.13×10⁴, or p ≳ 0.1 surviving. Cheap
precursors: the X ↔ Σ charge exchange followed to freeze-out (it sets N_X/s to better than 0.2%), and the X–core
scattering at k ≈ 0.07 m_X with Codex's own three-channel machinery.

## Status

M1 is a consistent, well-computed conditional model whose supplied inputs now include the dark-matter asymmetry's
sign, a Wilson coefficient at the Hubble scale, the object size and the pulse that makes its gas. It has no distinctive
prediction yet. Candidate 1 is unchanged; the closure principle supplies none of this (SYNTHESIS.md of the closure
pass).
