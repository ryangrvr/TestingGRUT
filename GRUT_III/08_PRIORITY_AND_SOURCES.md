# 08 — Priority and sources: what the literature already contains (equation-level, where the text could be read)

Method: each paper was fetched in full text (arXiv PDF → text) by an independent agent asked a specific question
and required to quote equations. "FULLTEXT" means the body was read; verdicts are the agent's, checked by the author.
Entries marked *pending* were still in progress when this file was written and are appended as they arrive.

## Pérez-Nadal, Roura, Verdaguer (2009), "Stress tensor fluctuations in de Sitter spacetime", JCAP 05 (2010) 036 — FULLTEXT

- **Treats:** the free **minimally coupled** (ξ = 0) scalar of arbitrary mass in Bunch–Davies on dS_D, in the
  dS-invariant bitensor basis F = P nnnn + Q(…) + R(…) + S(…) + T gg (eq. 28), with closed-form coefficient functions
  (eqs. 40–43), elementary-function limits (small mass eqs. 48–50; strictly massless eqs. 55–56 in the Allen–Folacci
  vacua; large distance eqs. 58–63) and the timelike continuation Z → Z ∓ iε (Sec. 7).
- **Does not treat** the conformally coupled field: "In this paper we have considered only minimally coupled fields,
  i.e., ξ = 0" (Sec. 6.2); general ξ is deferred (Sec. 8). No worldline restriction, no proper-time Fourier transform,
  no KMS discussion.
- **Bearing:** revision 1's and R9's statement that the conformal-field dS correlator "is in PNRV 2009" was
  **inaccurate** and is corrected here: PNRV supply the general dS-invariant framework and the ξ = 0 case (relevant
  to D10's massless sector, where the worldline propagator is Allen–Folacci's). **Priority verdict: NOT THIS RESULT.**

## Cho & Hu (2014), "Noise kernels of stochastic gravity in conformally-flat spacetimes", CQG 32 (2015) 055006 — FULLTEXT

- **Contains:** closed analytic noise kernels of a conformally coupled scalar in Minkowski, Einstein and open
  Einstein spaces (Phillips–Hu formula eqs. 18–21; coefficient functions C_ij eqs. 23–28; Minkowski set eq. 29), with
  flat/closed/open FRW for **arbitrary a(η)** by the conformal rule (C_ij)_fFRW = a⁻²(η)(C_ij)_M a⁻²(η′) (eq. 30);
  de Sitter in flat slicing with Bunch–Davies "conformal to the Minkowski one" (Sec. III.A), global and static dS
  (eqs. 37–40). The conformal rule Ñ = Ω(x)⁻² N Ω(x′)⁻² (eq. 17) is credited by them to **Eftekharzadeh, Bates, Roura,
  Anderson, Hu, PRD 85, 044037 (2012)**. Thermality at T = H/2π is stated.
- **Equation-level check:** (C₁₁)_M at Δs = 0 equals 3Δη⁴/(12π⁴Δη¹²) = **1/(4π⁴Δη⁸)**, the same constant as `v1` §1.4;
  with a = −1/(Hη) and the frame conversion, GRUT III's worldline kernel (1/4π⁴)[H/(2 sinh(HΔt/2))]⁸ follows in two
  lines (verified symbolically to zero difference by the checking agent).
- **Does not contain:** the worldline sinh⁻⁸ form written out, any proper-time Fourier transform, the odd-polynomial
  spectrum, or the 4H decay rate (zero hits for "Fourier/spectrum/worldline/proper time/10080/560").
- **Priority verdict: CLOSELY RELATED; the worldline kernel is a corollary of their eqs. (29)+(30) and cannot be
  claimed as new; the proper-time spectrum G⁺_T(ω) is not printed there.** D2 and R9 now say exactly this.

## Satin, Cho, Hu (2016), "Conformally-related Einstein–Langevin equations …", PRD 94, 064019 — FULLTEXT (arXiv:1603.06206)

- **Contains:** the conformal transformation of the influence action with the anomaly terms (eqs. 6–8, a = 1/1920π²,
  b = −1/5760π²); of ⟨T^ab⟩ (eqs. 15–16, agreeing with Horowitz–Wald); the **noise-kernel rule** Ñ^abcd(x,y) =
  e^{6ω(x)} N^abcd e^{6ω(y)} (eq. 29, credited to Eftekharzadeh et al. 2012); the **dissipation-kernel rule**
  H̃_A = e^{6ω(x)} H_A e^{6ω(y)} (eq. 40, new there); and the statement that the fluctuation–dissipation relation is
  preserved under conformal transformation (eqs. 41–43).
- **Bearing:** this is exactly the mapping step GRUT III uses ("T̂ = a⁻⁴T_flat on each leg", D2) and the reason the
  noise/dissipation pairing of `07` T1 is inherited from Minkowski in dS. Revision 1 listed SCH 2016 among sources for
  the dS position-space kernel; that was inaccurate — SCH give the *rules*, not the dS kernel. **Priority verdict:
  CLOSELY RELATED for the mapping step; NOT THIS RESULT for the worldline spectrum.**

## Masanes, Galley, Müller (2019), "The measurement postulates of quantum mechanics are operationally redundant", Nat. Commun. 10, 1361 — FULLTEXT

- **Assume:** complex separable Hilbert space (finite or countably infinite d), unitaries as reversible
  transformations, tensor-product composition; operational closure of outcome-probability functions under mixtures,
  unitaries and composition (eqs. 7–9); a ⋆-product for bi-local outcomes that is **associative** (eq. 14, motivated by
  observer-independence of subsystem partition); and "possibility of state estimation" (finite-dimensional systems
  are characterised by finitely many outcomes). They explicitly do **not** assume non-contextuality, operators for
  outcomes, or density-matrix mixed states.
- **Derive:** the POVM structure, the Born rule and the state-update rule.
- **Bearing:** GRUT III's P5 (non-contextuality → Gleason/Busch) is a *weaker* route than MGM's; 01 P5 and 05 O3 now
  say so. "None derives probabilities from causality alone" is correct for MGM: their inputs are Hilbert space,
  unitarity, composition, associativity and state estimation. **No priority conflict** (GRUT III claims no derivation).

## The Stage-3 charter rule the reviewer cited (located by the branch-digest agent)

Not in this repository. Found in the owner's `vsctestinggrut` repository at `PROGRAM/STAGE3/STAGE3_CHARTER.md`
(blob `8dd7b05…`, commit `c1275a8`). Relevant text as quoted by the agent:
- §13.1: "The standard constraint set 𝒮 is imposed **before** any freedom is counted. Relations implied by 𝒮 earn
  nothing."
- PS-5: "Gibbs/KMS/FDT structure … any input from which a standard theorem … derives a KMS state, detailed balance …"
- NR-12(b): "If any premise is a PS-5 structure …, whether supplied, hidden or generated by 𝒦, then every claim whose
  derivation passes through it …" [earns nothing].
- §13.6 "Precedents and RECOVERY": Kubo/FDT relations fixed by P0's connected correlators lie in 𝒮.
- **DEF-16 (self-referential flows):** "KMS, detailed-balance or FDT forms relative to a flow constructed from the state
  itself (modular or thermal-time constructions …). KMS-type credit or RECOVERY needs a time evolution certified
  independently."
- §14 SCOPE-Q (G2-12 item 4): "scope-Q cards target a platform where RH-KMS is certified to fail."
- CT-75: "Generated-premise laundering: 𝒦 generates detailed balance or KMS … and ℛ★ is its FDT/Onsager consequence
  → STANDARD-IMPLIED."

**Bearing:** these rules say, in the owner's own earlier words, what revision 2 says: modular-flow thermality (P1) is
self-referential and earns nothing (DEF-16); FDT consequences are standard-implied (NR-12(b), CT-75); and a distinctive
claim needs a platform where KMS is certified to fail (SCOPE-Q) — the quasi-de Sitter program of `05` O5.
The revision-2 texts that said "§13.6 not found in this repository" are correct as stated and are now supplemented
by this location.

## Constructive-physics branch (CP-01 … CP-04; digest), and CP-05 (`codex/grut-causal-response-05`)

Common pattern (digest): geometry and quantum structure supplied; one mechanism derived; a complete conventional
comparator with all EFT operators reproduces every soft observable; verdicts "CONDITIONAL ONLY". The branch's standing
rule: "Deriving properties of an input is not deriving the physical input." CP-01: curvature-susceptibility scalar,
cosh potential, λ₄ = (2/3)m²/M_P², family not closed under its own scalar loop. CP-02: shared induced-geometry
mechanism, C₇ = 3κ²/8π³, no Einstein limit. CP-03: clockless scalar determinant mechanism, real-signature obstruction.
CP-04: global vacuum-constraint construction. CP-05: F = M² + 2g_cφ scalar–tensor, (□ − μ²)φ ∝ T, covariant memory =
the scalar's retarded propagator with the J₁ tail — "WORKING LEADING SCALAR–TENSOR MODEL; NO NEW GRUT LAW". Its review's
three qualifications (Cassini as a measurement not a tolerance; Eöt-Wash 38.6 µm only at gravitational strength;
"every causal stable kernel is an integrated-out field" unproved in the strong sense) are incorporated in `03` and
`07` T6.

## DESI Collaboration (2025), "DESI DR2 Results II: BAO and cosmological constraints", arXiv:2503.14738 — FULLTEXT

- **Significance of w₀w_aCDM over ΛCDM** (Table VI; Δχ²_MAP, two extra d.o.f., Wilks): DESI alone 1.7σ; DESI+CMB 3.1σ
  (Δχ² = −12.5, driven mainly by CMB lensing); DESI+CMB+Pantheon+ 2.8σ; +Union3 3.8σ; +DESY5 4.2σ. "No combination
  exceeds 4.2σ." Best fits: DESI+CMB w₀ = −0.42 ± 0.21, w_a = −1.75 ± 0.58; +Pantheon+ w₀ = −0.838 ± 0.055,
  w_a = −0.62; +Union3 w₀ = −0.667 ± 0.088, w_a = −1.09; +DESY5 w₀ = −0.752 ± 0.057, w_a = −0.86. Favoured w(z) is
  > −1 at low z with a phantom crossing above z ≈ 0.4. Within ΛCDM, DESI BAO and CMB parameters are in 2.3σ tension.
- **DR3:** the string does not occur in the paper; only "future DESI analyses". Revision 1's "DR3 decides" was wrong.
- **Bearing:** data source for F1 only; GRUT III's revision-2 wording is accurate. Note that the favoured history
  *crosses* w = −1, which the record's own single-relaxor analysis already excluded for GRUT I.

## Milgrom (1999), "The modified dynamics as a vacuum effect", Phys. Lett. A 253, 273 — FULLTEXT

- Takes T(a) = (1/2π)(a² + Λ/3)^{1/2} from Narnhofer–Peter–Thirring [27] and Deser–Levin [28] (eq. 6) and the
  5-acceleration a₅ = (a² + Λ/3)^{1/2} (eq. 7); reads ΔT = T(a) − T_Λ as MOND inertia; derives 2πΔT = a·μ̂(a/â₀) with
  μ̂(x) = [1 + (2x)⁻²]^{1/2} − (2x)⁻¹, **â₀ = 2(Λ/3)^{1/2}** (eqs. 8–9), and a variant with μ(x) = x/(1+x²)^{1/2},
  **a₀ = (Λ/3)^{1/2}** (eqs. 10–11; abstract). "I can offer no specific mechanism that begets MOND inertia from vacuum
  effects." Numerically c√(Λ/3) = 5.4–5.9 × 10⁻¹⁰ m s⁻² (4.5–4.9 × observed a₀); 2c√(Λ/3) is 9–10 ×.
- **Bearing:** `06` corrected (coefficients are definite and overshoot; the cH₀/2π form is Milgrom's later canonical
  statement). The a₀ ∝ H(z) test is implicit in any identification with the Hubble scale rather than Λ and has been
  discussed in the MOND literature; `06` records it as context only. **Priority: Milgrom's, as revision 2 states.**

## Han, Lykken, Zhang (1999), "On Kaluza–Klein states from large extra dimensions", PRD 59, 105006 — FULLTEXT

- Conventions: ĝ = η + κ̂ĥ with κ̂² = 16πG (after eq. 2), kinetic coefficient 1/4 (eq. 2), polarisation sum without the
  1/2 (eq. A.5; footnote [19] corrects an earlier version) — equivalent to the canonical convention of `v4`.
- Widths (Sec. 3.1.1): massless vectors Γ = Nκ²m³/160π (eq. 46); massive vectors δκ²m³(1−4r)^{1/2}(13/12 + (14/39)r +
  (4/13)r²)/80π (eq. 49); Dirac fermions N_cκ²m³(1−4r)^{3/2}(1+8r/3)/320π (eq. 52); complex scalar
  κ²m³(1−4r)^{5/2}/960π (eq. 53) ⇒ real scalar 1/1920π. Ratios real scalar : Dirac : vector = 1 : 6 : 12.
- **Bearing:** D5's threshold laws and weights are **contained** in HLZ; `v4` reproduces them from explicit amplitudes
  (a verification, not a discovery, as `v4` and D5 say). The massive-vector polynomial is quoted above and not computed.

## Jacobson (1995), PRL 75, 1260; Eling, Guedens, Jacobson (2006), PRL 96, 121301 — FULLTEXT

- Jacobson: an arbitrary point p; a 2-surface element P whose null congruence has vanishing expansion and shear at p;
  "Through any spacetime point there are local Rindler horizons in all null directions"; χ^a the approximate boost
  Killing field, T = ħκ/2π "where κ is the acceleration of the Killing orbit on which the norm of χ^a is unity";
  δQ = −κ∫λT_kk dλdA (eq. 2), δA = −∫λR_kk dλdA (eq. 5); "δQ = TdS … can only be valid if T_kk = (ħη/2π)R_kk for all
  null k^a"; the κ scale ambiguity is explicitly noted as harmless. The **only** validity condition is curvature radius ≫
  cutoff length; **no condition on any acceleration appears** — confirming R1.
- **Bearing:** D3 is **entirely contained** in Jacobson 1995 and is labelled IMPORTED. EGJ 2006 (shear/viscosity) read;
  the identification of the horizon shear viscosity underlies D4's membrane statement.

## Narnhofer, Peter, Thirring (1996); Deser & Levin (1997) — FULLTEXT (Deser–Levin), abstract (NPT)

- 2πT = (Λ/3 + a²)^{1/2} = (R⁻² + a²)^{1/2} (NPT abstract; Deser–Levin abstract and eq. 8, who write "For dS, this
  recovers a known result" and credit [2] = NPT). Deser–Levin add the 5-acceleration reading, the Tagirov argument that
  the conformally coupled massless Wightman function in dS equals the flat 5-space one, and the AdS extension
  2πT = (a² − R⁻²)^{1/2}. **Bearing:** `06`'s formula is correct and correctly attributed; priority NPT 1996.

## Independent review of revision 2 (owner's reviewers, 2026-10-11)

- Published check: `research/reality_framework/grut_iii_independent_review/REVIEW_V11_R2.md` and
  `check_v11_powerlaw_expansion.py` on branch `codex/grut-causal-state-review` (commit `0aca9229`). Result: at
  x = 1 − i the exact power-law correlator differs from W_dS(H(1−ε)) by 2.77e-4, 6.93e-5, 1.74e-5, 4.34e-6 at
  ε = 0.04, 0.02, 0.01, 0.005 (quartering with halving ε: second-order scaling); the leading coefficient
  ε²[x²/4 − (x³/12)coth(x/2)] reproduced numerically. Second reviewer: an independent implementation reproduces
  v11's p = 100, 200 values to all digits; p = 400, 800 + Richardson give D/ε → 2πω/H̄ to 3e-5.
- Independent reruns: all eleven scripts pass; v7 matches the smeared-detector result; v9 reproduces Solodukhin's
  k₁ (scalar 1/6 − ξ, Dirac +1/3); v10 reproduces D = H³/4π². The literature entries in this file were not
  independently checked by them.
- Accepted in full; `v11` rewritten (rev. 3); see RESPONSE R15.

## Verification status of the new scripts (honest record)

The three-lens adversarial verification of `v7`–`v10` scheduled in the round-two workflow **did not run**: the four
verifier agents failed at launch with a usage-limit error. In their place the author (i) re-ran every script
(`RESULTS.md`), (ii) re-derived by hand the ρ_T polynomial expansion used in `v8` ((ω²+H²)(ω²+4H²)(ω²+9H²) =
ω⁶ + 14H²ω⁴ + 49H⁴ω² + 36H⁶), the Maxwell heat-kernel coefficient k₁ = −2/3 (4·(1/6) − 1 − 2·(1/6)) and the species
sums (4 + 45 − 48)/6 = 1/6 and (0 + 45 − 48)/6 = −1/2 used in `v9`, and the Gaussian-convolution and lens-volume
steps of `v7` §A, and (iii) checked `v7` §B's two independent closed-form derivations against each other as the
script reports (agreement to 1e-26). `v7`'s own caveats (the 1/t edge; L⁻⁸ scaling) were adopted into D9. `v10` and
`v11` were author-written; `v11`'s revision-2 headline was then **refuted** by the independent reviewers and corrected
(rev. 3). The reviewers' reruns cover v7, v9, v10 and the power-law part of v11; `v8` and the general-a(t) and ΛCDM
parts of `v11` rev. 3 remain author-level.
