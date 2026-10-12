# 02 — Derivations (revision 2)

Units ħ = c = k_B = 1 unless shown. Each item names its tier: IMPORTED (standard result; zero credit),
DERIVED (computed in `calc/` from the postulates plus imports; still standard physics unless marked),
RECORD CORRECTION (closes an item in the GRUT record). Scripts: `calc/v1`–`v11`; outputs in `RESULTS.md`.

## D1 — The medium's state is thermal at the horizon temperature (IMPORTED; checked in `v1`)

For a conformally coupled scalar in the Bunch–Davies state on dS₄ (a = −1/(Hη)), along one comoving worldline

    g(τ) = −(H²/16π²) / sinh²(H(τ − iε)/2),      because  a(η)a(η′)(η−η′)² = 4 sinh²(HΔt/2)/H²  exactly.

The image-sum identity Σₙ(z + inβ)⁻² = (π/β)²/sinh²(πz/β) shows this *is* the flat thermal correlator at
β = 2π/H (Gibbons–Hawking 1977). Shifting the contour to Im τ = −β/2 (where the integrand is smooth) gives

    G⁺(ω) = (ω/2π)/(1 − e^{−2πω/H}),    S_sym(ω) = G⁺(ω) + G⁺(−ω) = (ω/2π) coth(πω/H).

**RECORD CORRECTION (proved).** `calc/worldline_reduction.py` defines its "symmetrised" spectrum as
S(w) = 2 Re ∫₀^∞ e^{iwτ} g(τ) dτ. For any Wightman function g(−τ) = g(τ)*, hence 2 Re ∫₀^∞ = ∫_{−∞}^{∞} = G⁺(w):
the record computed the **one-sided** Planck spectrum and compared it with two **symmetrised** closed forms. Its
value 0.03393 at w = 0.1H is G⁺(0.1H) = 0.034116 up to a 0.55% ε-extrapolation error (ε ∈ {0.04 … 0.005}). The
"UNRESOLVED distinct floored spectrum" is the Planck spectrum; its floor is the classical Nyquist floor
2T·(Im χ/ω)(0), not "zero memory".

## D2 — The stress-tensor noise kernel, in closed form (DERIVED; corollary of standard results)

Explicit Wick contraction of the improved (ξ = 1/6) stress tensor (`v1` 1.4; general point in `v7` B.0):

    ⟨T₀₀(x) T₀₀(0)⟩_c = C_T (I₀₀² − ¼)/(x²)⁴,  C_T = 1/(3π⁴),  I₀₀ = 1 − 2t²/x²;   on the time axis 1/(4π⁴ t⁸).

Mapping with D1 (orthonormal-frame components T̂ = a⁻⁴ T_flat on each leg) the geodesic observer in dS sees

    ⟨T̂₀₀(t) T̂₀₀(0)⟩_c = (1/4π⁴) [H/(2 sinh(H(t − iε)/2))]⁸,

and, by the contour shift plus the recursion F_{2n+2} = F_{2n}(4n²a² + ω²)/(2n(2n+1)a²) for ∫e^{iωu}sech^{2n}(au)du,

    G⁺_T(ω) = ω (ω² + H²)(ω² + 4H²)(ω² + 9H²) / ( 10080 π³ (1 − e^{−2πω/H}) ).

Checks (`v1` 1.5: quadrature agreement ≤ 1e-16 and the 4H decay; `v8` B: FDT/KMS algebra on the closed form): KMS G⁺_T(ω)/G⁺_T(−ω) = e^{2πω/H}; flat limit
ω⁷/(10080π³); G⁺_T(0) = H⁷/(560π⁴); asymptotic decay rate of the time-domain kernel 4H.

**Reading, with the labels the reviews required.** (i) Occupation is Planckian at T = H/2π. (ii) The
commutator spectrum ρ_T = G⁺_T(ω) − G⁺_T(−ω) = ω(ω²+H²)(ω²+4H²)(ω²+9H²)/(10080π³) is an **odd polynomial**: the
thermal denominator cancels exactly, the pointwise retarded response is a finite sum of δ-derivatives (Huygens),
and all thermal structure sits in N_T = ½coth(πω/H)ρ_T. (iii) **4H is the asymptotic decay *rate*** of the pointwise
correlator (first surviving pole after the polynomial cancels the Matsubara rungs n = 1, 2, 3); it is *not* an
autocorrelation time — the pointwise correlator diverges as t⁻⁸ and has no normalisable one (R2). The thermal scale
β/2π = 1/H is a different quantity. "Memory" is defined only for smeared observables (D9, `07`). (iv) For
spatially extended or homogeneous modes the response integrates the light cone to the horizon and acquires the
record's power-law tails (48/πΔ⁵ …); those are not a medium "with memory" either.

**Priority (R9; checked at equation level in `08`).** The position-space conformal-field noise kernel in
conformally flat backgrounds is **Cho–Hu (2014)**, eqs. (29)+(30) — their (C₁₁)_M at Δs = 0 is exactly 1/(4π⁴Δη⁸), the
constant found independently in `v1` — with the conformal rule Ñ = Ω⁻²NΩ⁻² due to Eftekharzadeh–Bates–Roura–Anderson–Hu
(2012) and the dissipation-kernel rule and FDR preservation in Satin–Cho–Hu (2016). (Pérez-Nadal–Roura–Verdaguer 2009
treat the minimally coupled field only.) The worldline sinh⁻⁸ form is a two-line corollary and is not new; the
proper-time spectrum G⁺_T(ω), the odd-polynomial dissipative part and the 4H rate are not printed in those papers and
are recorded as corollaries/placements. "Memory lives in the noise kernel" is Hu–Verdaguer stochastic gravity.

## D3 — Einstein's equation as the medium's equation of state (IMPORTED: Jacobson 1995; reproduced in `v2`)

At a point p take any null k^a and the local Rindler horizon through p with boost field χ^a = −κλk^a; Raychaudhuri
with θ(p) = σ(p) = 0 gives δA = −∫λR_kk dλ dA; the heat flux is δQ = −κ∫λT_kk dλ dA; Clausius with T = κ/2π and
S = ηA gives T_kk = (η/2π)R_kk — **κ cancels**, and nothing refers to any observer's physical acceleration (R1).
The null-cone lemma (SVD nullspace = metric, `v2` 2.2) and the Bianchi identity give

    G_ab + Λ g_ab = 8πG T_ab,   with  η = 1/4G inserted,  Λ an integration constant.

GR is **recovered** at the price of P2 plus the datum G. Λ is not a response of the medium.

## D4 — Horizon transport and the thermal scale (IMPORTED; numbers in `v3`)

The membrane paradigm's **surface** shear viscosity of a horizon is η_s = c³/(16πG) = 8.03 × 10³³ **kg s⁻¹**
(a 2-d fluid quantity; revision 1 mislabelled it Pa·s, R3), the entropy per area is s = k_B c³/(4Għ), and
η_s/s = ħ/(4πk_B). This describes the horizon as an effective 2-d fluid; it does not make empty spacetime a bulk
fluid. The thermal scale τ_β = ħ/(2πk_B T) = 1/κ equals 1/H₀ = 4.58 × 10¹⁷ s for the cosmological horizon and
4GM/c³ = 19.7 µs for a solar-mass black hole. Static-patch identities (exact): E = ρ_Λ V = 1/(2GH) = TS, so
F = 0; |E_Komar| = 1/(GH) = ½Nk_BT with N = A/ℓ_P². These are standard.

## D5 — The matter-sector absorptive kernel from the spectrum (DERIVED; `v4`)

Unitarity with the vertex (κ/2)h_μν T^μν gives |Im Π(s)| = ½Σ_f S_f ∫dΦ₂ |M|², averaged over the five spin-2
polarisations. Explicit amplitudes: real scalar Σ|ε·T|² = 8p⁴/3 (law β⁵); Dirac fermion, massless ratio to the real
scalar **6.000000000**, law **β³(1 + 8r/3)** to 2e-16 (Han–Lykken–Zhang's form, reproduced); photon pair, ratio
**12.000000000**. Hence

    Im Π(s) = c_eff(s) · G s²,     c_eff(s) = (1/120) Σᵢ wᵢ Fᵢ(4mᵢ²/s),    w = 1 : 3 : 12 (real scalar : Weyl : vector),

the frequency-resolved Weyl-anomaly coefficient (c_SM = 283/120 above 2m_t). **Against the primary text (`08`):**
Han–Lykken–Zhang (1999) eqs. (46), (52), (53) give Γ(h̃→VV) = Nκ²m³/160π, Γ(h̃→ff̄) = N_cκ²m³(1−4r)^{3/2}(1+8r/3)/320π,
Γ(h̃→HH̄) = κ²m³(1−4r)^{5/2}/960π for a *complex* scalar (so 1/1920π per real scalar), with their κ² = 16πG and a
non-canonical kinetic term (coefficient 1/4, |ε|² = 2) — physically identical to the canonical convention used here;
the ratios 1 : 6 : 12 and the two threshold laws derived in `v4` agree. Their massive-vector factor is
(1−4r)^{1/2}(13/12 + (14/39)r + (4/13)r²) (eq. 49), **not computed here**. **Scope (R4):** the variable is the
timelike invariant mass √s, not a GW frequency — a free on-shell GW has s = 0 and *no* species is open at any
frequency; the Standard-Model staircase in `v4` is **illustrative** (massive-vector threshold polynomial not computed,
hadronic region 0.3–2 GeV schematic). The derived content is the weights, the two threshold laws, and the edge data:
locations 2mᵢ, exponents 5/2, 3/2, 1/2. The FRW remainder Γ ∼ c_eff t_P² H₀⁴/2ω ≈ 5 × 10⁻¹⁶² s⁻¹ reproduces and
explains the record's "invisible" verdicts.

## D6 — The anomaly sector (DERIVED numbers; `v5`)

a = (N_s + 11N_D + 62N_V)/360, c = (N_s + 6N_D + 12N_V)/120; SM: a = 2.765, c = 2.358. Anomaly-driven
self-consistent dS (Starobinsky 1980): H² = πM_P²/a → H = 1.07 M_P at SM content, excluded (r < 0.036 needs
a ≳ 2 × 10¹¹). GRUT III does not claim anomaly-driven inflation.

## D7 — The outcome rule (IMPORTED given P5)

Non-contextuality ⇒ Tr(ρE) (Gleason; Busch). The revision-1 route "P5 from P0 + P1 via HJW" is withdrawn (R5).

## D8 — Where the arrow enters (IMPORTED; located)

sign(β) = the choice of the future horizon as the traced-out one (P4). Not derived.

## D9 — Memory from access geometry (DERIVED; `v7`) — Theorem T3 of `07`

Free massless scalar, G_R = θ(t)δ(t − r)/4πr, Gaussian profile of width L: the convolution f_L ∗ f_L is the
Gaussian of width √2L and the radial integral against the light-cone delta gives, exactly,

    χ_L(t) = θ(t) · t · e^{−t²/4L²}/(8π^{3/2}L³),   ∫₀^∞ χ_L dt = 1/(4π^{3/2}L) (the smeared Coulomb potential),
    ⟨t⟩ = √π L,   ρ_L(ω) = (ω/2π) e^{−ω²L²}  (pointwise ω/2π times the form factor);

for a top-hat profile of radius R the smeared response vanishes identically for t > 2R (lens-volume formula; checked
at t/R = 0.5 … 3). For the conformal stress tensor the same smearing turns the pointwise polynomial response into

    χ_L(t) = −(t⁶ − 30L²t⁴ + 180L⁴t² − 120L⁶) e^{−t²/4L²} / (184320 π^{7/2} L¹³ t),

obtained two independent ways (position-space distribution and momentum-space CFT form) and confirmed by
quadrature to 1e-26 at nine values of t/L: a kernel of support ∼ L with zero crossings and a Gaussian tail, from an
observable whose pointwise response has no memory at all. **Nuance found by `v7`:** spatial smearing alone leaves a
1/t edge at small t (χ_L → 1/(1536π^{7/2}L⁷t), because the dissipative spectrum of a spatially smeared T₀₀ tends to the
constant √π/(1536π³L⁷) at high frequency — pairs of arbitrarily high energy with bounded total momentum), so the static
susceptibility ∫χ_L dt diverges logarithmically and an integrable kernel needs a *temporal* smearing as well; the
scaling is χ_L = L⁻⁸g(t/L) (F_L has dimension 4). KMS for smeared observables, N_L = ½coth(βω/2)ρ_L, verified to
1e-13 at two temperatures. *Corollary:* apparent memory of order L/c can always be manufactured by the access
profile and is never evidence for a vacuum relaxation law.

## D10 — The gravitational sector: stationary thermal process plus a Brownian zero mode (IMPORTED form; `v10`) — Theorem T5

Massless minimally coupled sector (each graviton polarisation) in the dS-breaking Bunch–Davies-type state
(Allen–Folacci 1987):

    G(t,t′) = (H²/8π²)[ 1/(2 sinh²(HΔt/2)) − ln(2 sinh²(HΔt/2)) + H(t + t′) + c ],   1 − Z = 2 sinh²(HΔt/2) on the worldline.

∂_t∂_{t′}G depends on Δt only (every velocity correlator and commutator is stationary); the stationary part is
periodic under Δt → Δt + 2πi/H (checked; full KMS analyticity in the strip is not verified here); the non-stationary
part **at large separation** (HΔt ≫ 1) is the Brownian covariance D·min(t, t′) + const with **D = H³/4π²** — an
asymptotic statement, not an exact covariance at all separations — i.e. d⟨φ²⟩/dt = H³/4π² (Vilenkin–Ford; Linde; Starobinsky), and for tensors
d⟨h_{ij}h^{ij}⟩/dt = 2H³/(π²M̄_P²) = H·P_T. **RECORD CORRECTION:** the record's "graviton channel has no stationary
reduced object / no derived memory time" (tt_worldline) becomes: all response data are stationary; the zero mode
random-walks with a derived diffusion constant. "Infinite memory" is the random walk; "finite memory" is thermal.

## D11 — The conventional KMS structure under slow roll: exact first-order thermality, O(ε²) defect, present-day baseline (DERIVED; `v11` rev. 3)

**Correction of revision 2.** The earlier D11 reported a first-order temperature-shift coefficient 1.01 and a
"1%-of-ε frequency dependence" from finite-ε tables without extrapolating. Two independent reviewers (one with a
published check, `research/reality_framework/grut_iii_independent_review/` on `codex/grut-causal-state-review`) showed
that this was O(ε²) contamination. They were right; here is what is true.

**Theorem (first order; any a(t) whose slow-roll parameters ε = −Ḣ/H² and η = ε̇/Hε are both first-order small).** For a conformal field in the conformal vacuum along a comoving worldline,
W = −1/(4π²a(t)a(t′)(Δη − i0)²) exactly. Expanding about the midpoint t̄ to first adiabatic order (ε = −Ḣ/H², with
η = ε̇/Hε also first-order small), the correction to the de Sitter form is a **pure shift of the Hubble rate**:

    W(t,t′) = W_dS(Δ; H_eff) + O(ε²),     H_eff = H̄ + Ḣ̄/H̄ = H̄(1 − ε̄) = (d/dt) ln ȧ |_{t̄},

so the state is **exactly KMS to first order** at T_eff = (ä/ȧ)/2π (a statement at fixed x = HΔ, i.e. for ω ≫ εH; the
expansion is not uniform in Δ; and a statement about the fixed-mean-time Wigner spectrum — an instantaneous one-sided
Unruh–DeWitt rate differs at O(ε) by a frequency-dependent Ḣ term, see D13 "Convention" and `11`) (`v11` A: for the power law a ∝ t^p,
ln[a₁a₂Δη²/(4 sinh²(x/2)/H_eff²)] has no 1/q term, q = p − 1; `v11` B: the general-a(t) construction to second order
has vanishing first-order remainder at random parameter points, to 1e-15).

**Second order.** For the power law, W/W_dS(H_eff) = 1 + ε²[x²/4 − (x³/12)coth(x/2)] + O(ε³) with x = H_effΔ
(symbolic, and exact-expression check at q = 10⁵ to 1e-11); f₂(x − 2πi) ≠ f₂(−x), so **KMS is first broken at
O(ε²)**. For general a(t) the second-order remainder has an additional εη term (printed by `v11` B).

**Spectrum (Richardson-extrapolated, p = 200…1600).** The first-order defect relative to H̄ is exactly
D/ε → 2πω/H̄ (c = 1 to 3e-5), i.e. a temperature shift and nothing else. The genuine second-order defect relative to
the exact H_eff, D₂(ω) = lim[ln(G⁺(ω)/G⁺(−ω)) − 2πω/H_eff]/ε², is 2.50, 11.00, 31.47, 52.14 at ω/H = 0.5, 1, 2, 3;
divided by 2πω/H it rises from 0.80 to 2.77: a **frequency-dependent** KMS violation at O(ε²), the first conventional
one. Consequence for `03`/`07`: the conformal sector is KMS through O(ε), so the FDT boundary applies through first
order; a distinctive residual in slow roll would have to be O(ε²) or live in a sector where this expansion does not
hold.

**The present universe is not slow-roll.** Today ε₀ = (3/2)Ω_m = 0.47, so H_eff(t₀) = 0.53 H₀ and no expansion applies.
Along the comoving worldline a correlator exists only for |Δ| < 2t₀ = 1.9/H₀ (the big bang), so a Fourier spectrum at
t₀ requires a window and is not computed. In the time domain the exact ΛCDM correlator deviates from the first-order
thermal form W_dS(Δ; H_eff(t₀)) by 0.12%, 0.56%, 3.6% and 19% at ΔH₀ = 0.25, 0.5, 1, 1.5 — factors 4–18 above the
η = 0 second-order form ε₀²f₂ (the audit `11` notes that including the εη term of `v11` B, with η₀ = −2.06, makes it 4–22;
the conclusion is unchanged): the present-day comoving conformal vacuum is **not thermal at any single temperature** and
must be treated with the exact kernel. (A formal baseline for the vacuum medium; the physical photon field at these
scales is dominated by the 2.7 K CMB bath.)

## D12 — The induced-gravity / entanglement matching behind P2 (IMPORTED theorems, numbers in `v9`) — O1 made precise

In the heat-kernel scheme (cutoff ε), per species: minimal scalar δ(1/4G) = S_ent/A = 1/(48πε²); Dirac
1/(24πε²); Weyl 1/(48πε²); conformal scalar δ(1/4G) = 0 while S_ent/A = 1/(48πε²), the mismatch −ξ/(8πε²) being
exactly the Wald contact term −2πξ⟨φ²⟩_ε (Solodukhin eqs. 141–144); Maxwell: S_ent/A = 1/(24πε²) plus Kabat's
contact term −1/(8πε²) gives δ(1/4G) = −1/(12πε²). Species sums for the Standard Model: Σk₁ = +1/6 (ξ = 0) or
−1/2 (ξ = 1/6), i.e. **1/G_induced is tiny or negative** (the 12 vectors, −2/3 each, nearly cancel 45 Weyl at
+1/6 each); with three ν_R the ξ = 1/6 sum vanishes exactly. Sakharov's one-loop induced gravity with this content
does *not* give G (Visser 2002 reaches the same conclusion). Coefficients are scheme-dependent (heat-kernel 1/48π vs
Susskind–Uglum 1/360π); the scheme-free statement is Cooperman–Luty's S_ent = ¼δ(M_P²)A. **P2 is consistent, not
derived; the finite G needs the UV completion.**

## D13 — The massless minimally coupled sector under slow roll: first-order thermality of the fixed-mean-time spectrum at (ä/ȧ)/2π, σ-free (DERIVED; `v14`; audited, `11`)

**Convention (from the audit).** Every statement here is about the fixed-mean-time (Wigner) two-point function W(t̄, Δ)
Fourier-transformed in Δ at ω ≠ 0 — the same convention D11 uses for the conformal sector. It is not a statement about an
instantaneous Unruh–DeWitt rate: that rate is one-sided in time, W(τ, τ−s) = W(t̄ = τ − s/2, s), and differs from the Wigner
spectrum at O(ε) by a frequency-dependent term proportional to Ḣ (in both sectors; structurally it is s·∂_H W, not a
temperature shift; the auditor's one-sided numbers for the conformal sector are c(E) = 0.24, 0.93, 3.97 at E/H = 0.25, 0.5, 1,
tending to the bare H̄ reading as E → 0 — not re-derived by the author). For this sector a linearly coupled detector is
non-thermal already at ε = 0 because of the secular zero-mode growth (Garbrecht–Prokopec 2004), which the Wigner convention
places at ω = 0 (a derivative-coupled detector removes it). "Which temperature at first order" is therefore a question with
a convention in it; the Wigner convention is the one the record's two-time kernels use (`07` §6), and it is the one answered.

**Representation.** By the `v12` 12.4 lemma the power-law massless modes are the de Sitter modes of index
ν = 3/2 + ε/(1−ε) ≡ 3/2 + δ up to a ratio of scale factors, so the comoving two-point function is the Bunch–Davies de Sitter
function of index ν (continued to ν > 3/2) times that ratio: W(t,t′) = (1−ε)²H(t)H(t′)Γ(3+δ)Γ(−δ)₂F₁(3+δ, −δ; 2; w)/16π² + C_IR,
w = 1 + (τ−τ′)²/4ττ′, τ = −1/((1−ε)aH); the IR-cutoff state differs from the continued function by the constant C_IR
(frozen super-horizon modes). Checked against a direct mode integral (`v14` B: 1e-11 where the integral converges; the
cutoff differences at ν = 1.75 to 1e-10, the remainder scaling as k_IR^{5−2ν}). Around the pole, ₂F₁ = 1 − δS₁(w) − δ²S₂(w) + O(δ³) with S₁ = w/2(1−w) − ln(1−w)
(Allen–Folacci's massless function) and S₂ = −½ln(1−w) + Li₂(w) + ln(1−w)/2w + ½ − ¾w/(1−w); Γ(3+δ)Γ(−δ) = −(2/δ)(1 + 3δ/2 + …).

**First-order theorem** (Richardson in ε against the exact closed form, `v14` C): at fixed mean time t̄,
W(t̄,Δ) = (1 + 3ε/2)·W_mmc^{dS}(Δ; ä/ȧ) + ε(H̄²/8π²)S₂(cosh²(H̄Δ/2)) + c₀ + c₂Δ² + O(ε²), c₂ = −(1−ε)³εH̄⁴/32π² (c₀ depends on t̄:
it is the O(1) secular growth of the equal-time variance, (H̄³/4π²)t̄ + C_IR, Δ-independent at fixed t̄). The conformal-factor
combination is D11's pure shift to H_eff = ä/ȧ; the only new first-order term is the dilogarithm on the zeroth-order
argument. Using H̄ or the Kodama–Hayward rate H̄(1 − ε/2) instead, or dropping S₂, leaves an O(ε) residual.

**KMS.** Any function of cosh²(H̄Δ/2) is symmetric about the strip midline Im Δ = −π/H̄ (its value there is a function of
−sinh²(H̄s/2), even in s), hence KMS at 2π/H̄ by itself; c₀ and c₂Δ² live at ω = 0. So the first-order KMS defect of the
ω ≠ 0 part at β = 2π/(ä/ȧ) is O(ε²), non-uniformly in Δ (for Δ ≪ 1/(εH̄)); `v14` D″ checks this on the exact closed form,
whose defect after removal of the ω = 0 polynomial scales as ε² at β_eff and as ε at the Kodama–Hayward and bare rates:
**the fixed-mean-time spectrum of the massless sector is thermal at T = (ä/ȧ)/2π + O(ε²), like the conformal sector's,
and the Kodama–Hayward rate is excluded at O(ε) — with no coarse-graining scale anywhere in the computation** (the
fixed-mean-time convention is a choice, but not a free parameter).
The O(ε) physics is a *shape* change of the Wigner spectral density, not a temperature change:
ρ(ω) = ρ_mmc(ω; H_eff)[1 + ε r(ω)], ρ_mmc = (ω² + H²)/2πω, with r(0) = 2 − 2 ln 2 = 0.614, r(H̄) = 0.461, r(4H̄) = 0.085, r → 0
in the UV (the S₂ pole 3/8π²Δ² cancels the 3ε/2 amplitude exactly: the Hadamard coefficient is unchanged, as it must be);
large-ω asymptote r → 3(2ω⁴ − ω² − 1)/(4ω⁴(ω² + 1)) ≈ 3/2ω² (auditor's derivation, checked in `v14` at ω = 8 to 1e-5).

**Byproduct (σ-free diffusion coefficient).** The large-Δ expansion of the fixed-mean-time function in 1/H̄ ≪ Δ ≪ 1/(εH̄) is
W = const − (D/2)|Δ| − (3εH̄⁴/32π²)Δ² + …, with D = (H̄³/4π²)[1 − (1 + 2 ln 2)ε] + O(ε²) = (H_eff³/4π²)[1 + (2 − 2 ln 2)ε]; D is
equivalently the coefficient of the 1/ω² Brownian tail of the Wigner spectrum. The Δ² term is the expansion of the secular
growth of the equal-time variance (ω = 0) and exceeds the O(ε) correction to the |Δ| term for H̄Δ ≳ 3, so D is extracted by
analytic separation, not from a plateau; reading D as "the zero-mode diffusion coefficient" is the usual stochastic
identification and is an interpretive step (a Markovian walk with slowly varying D has no such Δ² term). With that step, the
stochastic A₀(σ)H³ of `v13` reproduces D at σ* = e^{1−γ_E} = 1.526 (one value, by monotonicity): a calibration of the
coarse-graining scale for this observable at first order, not a derivation of σ (no γ_E survives; one ln 2 does). This is
the σ-free answer to the question `v13` posed badly.

**Scope.** Comoving worldline, fixed-mean-time spectrum at ω ≠ 0; massless minimally coupled scalar; power-law FRW with
constant ε. The extension to general slow-roll a(t) is **not** established here: the D11 argument covers only the
conformal-factor combination, and the new S₂ term comes from the constant index shift δ = ε/(1−ε), which has no direct
analogue when ε varies (the modes are then not Hankel functions of a fixed index). Gravitons share this field's mode equation
at linear order in dS and power-law FRW, so the same statement follows for the TT-gauge tensor two-point function in the same
convention; no tensor computation was done, and the question of a physical (Riemann-coupled) detector is not addressed.
Interacting (loop) corrections are O5 item 2. All of this is conventional QFT; it establishes nothing about P1 beyond
consistency.
