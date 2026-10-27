# Skeptic review of track cf_state (closure of the state, of time, and of the dark-matter history)

27 October 2026. Claude, key `cf_state_skeptic`, parallel pass for issue #13 (brief: `../../BRIEF.md`). This is
same-session AI work, not an outside replication. Git was used read-only and nothing was committed.

- Code: `skeptic/code/skeptic_checks.py`, independent of `code/`.
- Results: `skeptic/results/skeptic_checks.json`, with stdout in `skeptic_checks.stdout.txt`.
- Downloaded papers are kept in the scratchpad, not the repository.

**Labels.** [D] derived, [C] computed, [Q] quoted from a source read in this pass, [E] estimate, [A] argued,
[R] recalled.

**Scope.** Every data comparison is **retrospective**. Nothing here validates or excludes a model of the real world.
Every postulate is an added choice, not a consequence of "closed".

## 0. Bottom line

**The headline survives and gets stronger.** No postulate in this family gives GRUT content that is both
distinctive and testable. Every verdict label stands.

**The track reproduces exactly.** Its four scripts, rerun from a scratch copy, give byte-identical JSON.

**Refuted or wrong as stated:**
1. **P5, the one positive claim.** The track's "[D] closure-flavoured consequence" is not forced: "energy closure forbids
   unbounded white-noise CSL heating, so the X-ray bounds do not apply". See §2.1.
2. **P5 data status is out of date.** XENONnT (PRL 136, 120201 (2026), arXiv:2506.05507) gives
   - λ/r_C² < 3.0×10⁻³ s⁻¹m⁻² at 90% CL (3.7×10⁻³ at 95%);
   - **the GRW point is excluded at 9.1σ** for white noise [Q].

   The track says "the GRW point is still allowed" and quotes Majorana 2022 as the current bound. Both are superseded.
   The explicit baseline the track says P5 forbids is therefore already excluded by data.
3. **P6c, the forbidden baseline.** The track says P6c forbids "the externally switched −c_H H(t)² model". That model
   has an internal representative, so P6c forbids nothing. See §2.4.

**Wrong comparator or misapplied number (the verdict is unchanged):**
- **P6b.** The Sr-clock discriminator (2.5×10⁻²⁶ against 0.70) and the "Milburn τ < 9.7×10¹⁰ t_P" bound use a
  measurement that is blind, by construction, to one-clock (common-mode) dephasing (§2.3).
- **P7a.** The spontaneous-CP window "a/b ∈ (0.929, 1.083) at c/b = 13" applies to no actual program point (§2.5).
- **P7a.** The "correction 2" attribution to REVIEW_REHEATING item 3 is wrong. Only rp_preparation §4 says "holds for
  either" (§2.5).
- **P5 Part A.** The 8.9×10⁻¹⁶ "closed whole reproduces Lindblad" figure is a self-comparison, on a degenerate initial
  state (§2.1).

**Missed postulates (§3):**
- the no-boundary state;
- final or two-time boundary states;
- a CPT-symmetric whole (Boyle–Finn–Turok);
- semiclassical gravity sourced by ⟨T⟩ (Page–Geilker).

## 1. Reproduction and independent recomputation

| Quantity | Track | Skeptic | Method (independent) |
|---|---|---|---|
| All four result JSON files | – | byte-identical on rerun | scratch copy, same code |
| CSL rate, 10⁻¹⁴ kg sphere, GRW | 2.97245×10⁷ s⁻¹ (34 ns) | 2.97247×10⁷ s⁻¹ (ratio 1.000000) | position space: ball distance distribution p(s), not the k-space form factor |
| Same, at XENONnT 95% bound | – | 1.10×10⁷ s⁻¹ (91 ns) | same |
| Fein molecule CSL exponent (GRW / Adler) | 3.44×10⁻¹⁰ / 0.138 | 3.443×10⁻¹⁰ / 0.1377. XENONnT 95%: 1.27×10⁻¹⁰ | point formula |
| DP: E_G, τ (10⁻¹⁴ kg, 250 μm) | 7.8×10⁻³³ J, 13.6 ms | 7.767×10⁻³³ J, 13.58 ms. Granular part at R₀ = 4.5×10⁻¹⁰ m: 5.9×10⁻⁴¹ J | analytic |
| White-CSL heating, GRW | 5.0×10⁻⁴⁴ W per nucleon | 5.02×10⁻⁴⁴ W, which is **1.4×10⁻⁷ eV per nucleon in 13.8 Gyr** | – |
| GPP exponents (Sr 37 s / molecule / sphere) | 2.5×10⁻²⁶ / 1.1×10⁻⁴⁷ / 3.0×10⁻²³ | 2.467×10⁻²⁶ / 1.138×10⁻⁴⁷ / 3.04×10⁻²³. TeV needed for D = 10⁻³ in 1 s: 1.19 | GPP eq. (6) re-read in hep-th/0406260 [Q] |
| X̄ per Q at τ = 20 | 1.73778835311 | 1.73778835311 (relative 1.2×10⁻¹³). Parseval 1 − 2×10⁻¹⁶. E/E_min = 8.60007 | polar ODE (r, p_r, θ) with conserved Q, Radau, Gauss–Legendre Fourier quadrature (not FFT) |
| Release patches in today's Hubble volume | 3×10⁷⁵–5×10⁷⁶ | 1.30×10⁷⁶ (g* = 106.75), 1.27×10⁷⁶ (116.75). Instant reheating T_RH = 3.0–3.1×10¹⁴ GeV; residual 8.8×10⁻³⁹ | own constants (H₀ = 67.4, g_s0 = 3.937) |
| Spontaneous-CP criterion | 597/597 | **500/500 random, 200/200 within 0.5–3% of the boundary** | brute-force 721² or 1441² grid plus Nelder–Mead, no BFGS |
| Physical phases (dimension ≤ 4) | 0 (Z₂_φ exact) / 1 (broken) | same; **plus 1 from a renormalisable inflaton coupling χΣ*Y², even with Z₂_φ exact** | sympy rank |
| P6c representatives | ξR −2.54; ρ −10; V −5.0; K −1490 | ρ-representative ≡ −10H² in every era. **Switched representative −10H²·Θ(H − H_end): −10 during inflation, exactly 0 after** | own inflaton run (m²χ², χ₀ = 15) |

Primary values re-checked at source in this pass [Q]:
- **Majorana v5 text (arXiv:2202.01343).** λ/r_C² < 0.494 s⁻¹m⁻² (coherent nuclei, 95%), 17.4 (quasi-free electrons),
  R_DP > 2.54×10⁻¹⁰ m. The track's quotes are correct, but superseded.
- **XENONnT (arXiv:2506.05507 v2; PRL 136, 120201 (2026)).**
  - λ/r_C² < 3.0×10⁻³ (90%), 3.7×10⁻³ (95%), so λ < 3.0–3.7×10⁻¹⁷ s⁻¹ at r_C = 10⁻⁷ m;
  - the GRW parameter set deviates from the best fit by 9.1σ;
  - R₀ > 4.9×10⁻¹⁰ m (90%), 4.5×10⁻¹⁰ m (95%);
  - white (Markovian) noise is assumed.
- **Bothwell et al. (arXiv:2109.12238) full text.** "By comparing two uncorrelated regions within our atomic sample, we
  determine the contrast and relative phase difference between the two subensembles". The second π/2 pulse has "a
  random phase relative to the first", and an ellipse is fitted (Fig. 2; Methods).
- **GPP (hep-th/0406260).**
  - eq. (6): log(ρ₁₂(T)/ρ₁₂(0)) = −(3/2) t_P^{4/3} T^{2/3} ω₁₂²;
  - "The effect is too small to be observed in the lab";
  - the derivation evolves clock and system "with respect to an ideal time t".
- **Zhu et al. (arXiv:2605.26860, abstract).** λ < 1.27×10⁻⁵ s⁻¹ at r_C = 10⁻⁵ m. Confirmed.

## 2. Postulate-by-postulate audit

### 2.1 P5 (global purity plus a unitary whole). Verdict STANDARD stands. The positive claim is refuted.

1. **Energy closure does not force a finite internal bath [D, C].** Four independent reasons:
   - **Scope.** The derivation uses THM 2.1 with H bounded below and finite ⟨H⟩ for a whole that "includes gravity as a
     quantum subsystem". FORMAL_CORE itself rules this out:
     - §4 row "global energy conservation in cosmology: SCOPE";
     - DEF 1.6: "in general, no conserved total energy";
     - COND 5.4.

     On the compact T³ branch, the gravitational Hamiltonian is a constraint (H ≈ 0) and is unbounded below in the
     conformal mode. In open FLRW the total energy is not finite. The premise fails in the very setting P5 specifies.
   - **No finite-time content.** At GRW parameters, white CSL deposits **1.4×10⁻⁷ eV per nucleon over 13.8 Gyr** [C].
     The CMB alone holds 1.0×10⁶ eV per baryon, 7.6×10¹² times more [C]. Any internal bath with more capacity than
     this reproduces white-noise statistics at every finite time. "Unbounded heating" is a t → ∞ statement with no
     observable consequence. This fails Codex's gate: "Merely rewriting the same distributions using ...
     environment labels does not establish strictness" [Q 2908e91].
   - **The X-ray non-sequitur.** The keV bounds test the noise spectrum at ω ≳ 10¹⁸–10²⁰ s⁻¹ over years of
     exposure, not the long-time energy budget. A finite bath hotter than about 10 keV, with a cutoff above about
     10¹⁹ s⁻¹, emits the same X-rays.
     - Bilardello et al. (arXiv:1601.03672) state that the fate of the X-ray bound under dissipative CSL is unclear
       [Q via search summary].
     - The coloured-noise escape needs a frequency cutoff below about 10¹⁸ Hz [Q via search summary].
     - Closure fixes neither the temperature nor the cutoff. "The keV X-ray bounds then do not apply" should read
       "may not apply, depending on supplied bath parameters".
   - **Wrong relabelling.** Under P5-law there is no collapse at all. An internal bath that reproduces a CSL-like
     master equation is ordinary hidden-sector decoherence. Smirne–Bassi dissipative CSL is a collapse model with a
     classical noise. The closer credit for "the noise must be a physical field with its own energy" is the
     energy-conserving collapse literature (Pearle; Bassi–Ippoliti–Vacchini, J. Phys. A 38, 8017 (2005) [R]).
2. **The explicit forbidden baseline is already excluded.** XENONnT excludes white CSL at (10⁻¹⁶ s⁻¹, 10⁻⁷ m) at 9.1σ
   [Q]. So the track's discriminator separates P5 from a model the data already reject. P5 adds nothing. The data
   status line "The GRW point is still allowed" must be withdrawn.
3. **Part A is a weak test [C].**
   - The initial state |+⟩ is an eigenstate of H = σ_x/2 (‖[H, ρ₀]‖ = 0). The averaged dynamics is therefore pure
     dephasing, x(t) = e^{−2γt}, and the tunnelling never acts.
   - The 8.9×10⁻¹⁶ compares the collision model with a map that applies the same cos 2θ factor per step. Symbolically
     the single collision gives exactly r₀₁ cos 2θ, so the figure checks only the tensor contraction.
   - The track's "6×10⁻³ … is step error" is a θ-calibration error (θ² = γτ). With matched cos 2θ = e^{−2γτ} it falls to
     4×10⁻¹⁶ on the degenerate state.
   - **On a site-localised state |0⟩** the closed dilation converges to the exact Lindblad solution:
     - at order τ with Lie splitting, and τ² with matched θ plus Strang;
     - 8.8×10⁻⁴ for the full K = 12 dilation;
     - 5×10⁻⁸ at K = 1536.
   - The QSD and random-unitary averages agree within sampling error: 6.6×10⁻³ and 1.1×10⁻² against 1/√N = 1.1×10⁻²
     at N = 8000.
   - The claim survives. The number should be replaced.
4. **Reversibility is not operational [A].** The fidelity 1 − 10⁻¹⁵ applies the global inverse as one exterior
   operation. DEF 1.5(a) says an externally timed intervention opens the whole. Inside Ω, a reversal needs a
   controller and clock, subject to COND 2.13 and 2.20 and, under P6b, to GPP decoherence. "Every decoherence
   reversible in principle" is a property of the map, not of anything an internal observer can do.
5. **Also confirmed [C].**
   - DP τ = 13.6 ms;
   - CSL 34 ns at GRW;
   - the Fein exponents;
   - the heating figure.

### 2.2 P6a (ideal internal clocks). STANDARD, confirmed.

Nothing to refute. An ideal Page–Wootters clock reproduces standard quantum mechanics by construction. The
Salecker–Wigner algebra was checked in the track's sympy file.

### 2.3 P6b (GPP / Ng–van Dam clock bound). UNTESTABLE stands. The comparator is invalid.

1. **The Bothwell 37 s observable is blind to the effect [Q, D, C].**
   - **The GPP map is one-clock and common-mode.** It is ρ(T) = ∫dt P_T(t) U(t)ρU(t)†: a mixture of *common* time
     translations of the whole sample. (Milburn's [H,[H,ρ]] with total H is the same.)
   - **Bothwell measures a relative quantity.** They extract the relative phase and contrast of two sub-ensembles,
     with the laser phase randomised. That removes common-mode phase noise by design.
   - **Algebra [D, C].** For O = σ₊⁽¹⁾σ₋⁽²⁾, Tr(O·(−[H,[H,ρ]])) = −(ω₁ − ω₂)² Tr(Oρ). This was checked numerically: the
     ratio is −0.0049 against −(0.07)². For a single atom against the clock it is −ω².
   - **Monte Carlo [C]**, 2000 atoms per region, 4000 shots, single-atom D = 2 (coherence 0.135):

     | Jitter | Ellipse contrasts | Relative phase |
     |---|---|---|
     | common-mode (the GPP/Milburn form) | 0.992 / 0.995 | recovered (cos Δφ 0.821 against 0.825) |
     | independent per atom | 0.136 / 0.136 | – |

   - **Consequence.** GPP's prediction for the measured observable is D ≤ (Δω/ω)² × 2.5×10⁻²⁶ ≈ 2.5×10⁻⁵⁶ for
     Δω/ω ≤ 10⁻¹⁵ [E], not 2.5×10⁻²⁶. The "gap 2.8×10²⁵" and the **retrospective Milburn bound τ < 5.2×10⁻³³ s are
     not valid**: that measurement cannot see this class of dephasing.
   - The Ng–van Dam floor (3.8×10⁻³²) against 7.6×10⁻²¹ is likewise a synchronous differential comparison. It is
     also common-mode, so it is not a test.
2. **The valid comparator** references the atomic phase to an independent oscillator, the clock that defines T: a
   laser-referenced Ramsey experiment, limited by the laser coherence of about 1–10 s [E]. There, D_GPP ≈ 10⁻²⁷–10⁻²⁶.
   The verdict UNTESTABLE is unchanged and, if anything, reinforced. GPP themselves write "too small to be observed in
   the lab" [Q].
3. **The worst-case vertical ΔE is not what actual experiments do.** The track already flags this, but the table
   presents these rows as comparators. In the horizontal interferometers that were actually run, the branches have
   equal total energy, so the GPP prediction is 0.
4. **Tension with P5 [A].** GPP present realistic-clock non-unitarity as information loss for black holes. P5 lists
   "black-hole unitarity" as forced, and P6 forbids the external t on which P5-law's U(t) is built. Taken together,
   P5 + P6b predict effective irreversibility relative to every physical clock. P5's reversibility and black-hole
   statements hold only with the external time that P6 removes. The track treats the two separately and calls them
   compatible.

### 2.4 P6c (no externally timed switches). The forbidden baseline is refuted. Verdict: STANDARD ("forces nothing").

[D] In an expanding universe ρ̇ = −3Hφ̇² ≤ 0, so ρ, and with it H² = ρ/3M², is non-increasing. Any switching time t_s
can therefore be written as Θ(H − H(t_s)), a function of internal fields.

[C] In our own inflaton run (m²χ², χ₀ = 15):
- **The switched model.** The representative −10H²·Θ(H − H_end) equals −10 throughout inflation. It is exactly 0 after
  t_end = 17.42 m⁻¹, at χ_end = 0.995 M. ρ is monotone at every step.
- **The never-switched model.** This is c ρ/M² ≡ −10H² in every era (Friedmann).

Both the "externally switched" and the "never switched" histories therefore satisfy P6c, so S_C = S_base. The
track's own table already shows that the representative is a free function. The forbidden-baseline column should
read "none". This matches Codex's gate sentence on internal-clock rewriting.

### 2.5 P7a (preparation theorem). STANDARD stands. Two scope corrections.

1. **The spontaneous-CP window is applied at the wrong point [C].**
   - The track's ratio c/b = 13 is gp_phase's, computed at the **7737487 gauge point** (H_inf = 4.7×10¹³ GeV,
     M_Y = 10¹² GeV, h = 10¹¹ GeV, κ = 10⁻³), not at the heavy window.
   - At that point (gp_phase's quoted minimum):
     - a = 2κ|Y||X|³ and b = 2μ|Σ||Y|² give **a/b = 2517**;
     - c/b = 12.970 reproduces gp_phase's 12.972;
     - the triangle criterion fails, so there is no canted minimum.
   - At the **heavy spot point** (own minimum |X|, |Y|, |Σ| = 80.202, 3.4425, 28.311 H, reproducing hw_gates):
     - **a/b = 2.56**;
     - natural c/b = 1.1×10⁻³ (Σ*YX³/Λ, |coefficient| = 1, Λ = M_red);
     - the criterion fails.
   - A canted minimum there needs c/b ∈ (0.72, 1.64), that is Λ ∈ (1.7, 3.9)×10¹⁵ GeV. That is a tuned, supplied
     dimension-5 coefficient.
   - gp_phase §4 adds a dynamical objection [Q]. A supplied DRT-type mismatch of 0.1–0.3 rad "leaves no resolved net
     charge … the multi-field dynamics forget which side of the real subspace they started on". Ensemble means stayed
     |⟨ε⟩| ≤ 0.034 even with explicit CP phases.
   - So the proposed widening of REVIEW_REHEATING's "only route" is not supported at either program point. At most
     it is a logical possibility [A] that needs a tuned Λ and an untested holder.
2. **Misattribution.**
   - REVIEW_REHEATING item 3 already says the symmetry is "antiunitary conjugation … of time-reversal type, not C or
     CP" [Q]. It does not say "holds for either".
   - That sentence is in rp_preparation §4 [Q], which itself cites gp_phase's C version.
   - The track's substantive finding stands and is new [C]: the K-theorem holds only at the selection instant, and
     ⟨L⟩ moves from 11.19 to 11.65 for a K-invariant state.
   - The correction should be addressed to rp_preparation, not to REVIEW_REHEATING.
3. **"Exactly C-invariant" is a supplied truncation, not a property of the stated architecture.**
   - A renormalisable inflaton coupling χΣ*Y² adds a physical phase even with Z₂_φ exact [C]. gp_phase §3 already
     noted it.
   - Codex's THEORY.md keeps L_higher and states: "Additional allowed coefficients cannot be declared zero at every
     scale" [Q 04c53f5]. gp_phase counts 8 phases at dimension ≤ 5.
   - The forced negative is therefore conditional on supplied matching data. That is still Sakharov, so the verdict
     is STANDARD.

### 2.6 P7b, P7c, P7d. STANDARD and "not forced" stand.

- **P7b.** Nothing to refute.
- **P7c.**
  - Reproduced to 1.2×10⁻¹³ by a different formulation (§1).
  - New [C]: the 1.74 is **definition dependent**. An adiabatic split with the orbit-averaged mean-field frequency
    (ω_eff = 3.08) gives X̄/Q = **1.80** instead of 1.74, 3.8% higher.
  - This strengthens "not forced": closure cannot fix a number whose value depends on how particles are defined
    for a nonlinear orbit.
- **P7d.**
  - Order of magnitude confirmed (D/S ≈ 8×10⁻¹⁶).
  - The susceptibility uses g = 2 for one complex scalar, whose number susceptibility is T²/3. Correcting this
    doubles D/S to about 1.6×10⁻¹⁵. That is still 15 orders from binding.

## 3. Missed postulates in this family (closure of the state)

FORMAL_CORE §7 lists these as "known kinds outside" the table. No cf_* track covers them, apart from cf_vacuum's P_q.

1. **No-boundary state of the whole (Hartle–Hawking 1983 [R]).**
   - The unweighted weight is ∝ exp(8π²M²/H²). The exponent is 2.6×10¹⁶ at the heavy window's H_i and 2.1×10¹¹ at
     4.7×10¹³ GeV [E, C].
   - It prefers the lowest scale the potential allows, so a supplied heavy-window inflation is suppressed by
     exp(−O(10¹⁶)) against any lower available scale.
   - Volume weighting (Hartle–Hawking–Hertog, arXiv:0711.4630 [R]) reverses this.
   - Measure dependent. It is likely UNTESTABLE, or INCONSISTENT WITH the heavy-window branch under the bare
     measure. It is a relabelling.
2. **Final or two-time boundary states.**
   - Sources: Aharonov–Bergmann–Lebowitz 1964; Gell-Mann–Hartle (gr-qc/9304023 [R]); Horowitz–Maldacena
     (hep-th/0310281 [R]).
   - These are the most literal "closure of the state in time" strengthenings, and they can change Born statistics.
   - Not treated here.
3. **A CPT-symmetric whole** (Boyle–Finn–Turok, PRL 121, 251301 (2018), arXiv:1803.08928 [Q via search]).
   - This is a published closure-of-state symmetry postulate that **is testable**:
     - lightest neutrino massless;
     - light neutrinos Majorana;
     - no primordial long-wavelength gravitational waves;
     - dark matter a 4.8×10⁸ GeV right-handed neutrino.
   - It is the published realisation of the track's P7a loophole: a globally symmetric whole with an asymmetric
     branch.
   - Its dark-matter identification is not GRUT's X, so as stated it is INCONSISTENT WITH THE GRUT BRANCH. It is
     not GRUT-native.
4. **Unitary whole with classical gravity sourced by ⟨T⟩** (Page–Geilker, PRL 47, 979 (1981) [R]).
   - This is the one "pure unitary whole without collapse" variant already disfavoured retrospectively.
   - P5 avoids it only because it declares gravity quantum.
   - Struyve (arXiv:2510.20991 [Q abstract via search]) notes that some semiclassical theories give no GIE. That bears
     on P5's failure condition "absence of GIE".

## 4. Corrections the track should make

1. **P5 data status.** Replace "GRW still allowed; λ < 4.9×10⁻¹⁵ (Majorana)" with the XENONnT values (§1).
   - At the new bound, the 10⁻¹⁴ kg sphere's CSL time is 91 ns.
   - The Fein exponent at the new bound is 1.3×10⁻¹⁰.
   - Update the DP bound to R₀ > 4.5×10⁻¹⁰ m (95%).
2. **P5 "forces".** Drop the [D] energy-closure consequence, or relabel it [A], scope-limited, with no finite-time
   content. Drop "X-ray bounds then do not apply".
3. **P5 Part A.** Replace 8.9×10⁻¹⁶ with the convergence statement on a non-eigenstate initial state.
4. **P6b.**
   - Withdraw the Sr comparator, the Milburn τ bound and the Ng–van Dam versus 7.6×10⁻²¹ comparison as tests.
   - Use a laser-referenced Ramsey experiment [E].
5. **P6c.** Set the forbidden baseline to "none".
6. **P7a.**
   - Re-scope the spontaneous route as above.
   - Re-address correction 2 to rp_preparation §4.
   - Add the inflaton-coupling and L_higher conditions to "exactly C-invariant".
7. **P7c.** Note that the 1.74 depends on the particle definition: 1.80 with an adiabatic split.

## 5. Reproduction

```
cd research/reality_framework/grut_closure_forcing_20261027/tracks/cf_state/skeptic/code
python3 -B skeptic_checks.py ../results/skeptic_checks.json   # ~140 s; part G (grid brute force) dominates
```
