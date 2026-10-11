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

Checks (`v1` 1.5; `v8` B): quadrature agreement ≤ 1e-16; KMS G⁺_T(ω)/G⁺_T(−ω) = e^{2πω/H}; flat limit
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

the frequency-resolved Weyl-anomaly coefficient (c_SM = 283/120 above 2m_t). **Scope (R4):** the variable is the
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
quadrature to 1e-26 at eight values of t/L: a smooth kernel of support ∼ L with zero crossings, from an observable
whose pointwise response has no memory at all. KMS for smeared observables, N_L = ½coth(βω/2)ρ_L, verified to
1e-13 at two temperatures. *Corollary:* apparent memory of order L/c can always be manufactured by the access
profile and is never evidence for a vacuum relaxation law.

## D10 — The gravitational sector: stationary thermal process plus a Brownian zero mode (IMPORTED form; `v10`) — Theorem T5

Massless minimally coupled sector (each graviton polarisation) in the dS-breaking Bunch–Davies-type state
(Allen–Folacci 1987):

    G(t,t′) = (H²/8π²)[ 1/(2 sinh²(HΔt/2)) − ln(2 sinh²(HΔt/2)) + H(t + t′) + c ],   1 − Z = 2 sinh²(HΔt/2) on the worldline.

∂_t∂_{t′}G depends on Δt only (every velocity correlator and commutator is stationary); the stationary part is
KMS-periodic (period 2πi/H); the non-stationary part at large separation is exactly the Brownian covariance
D·min(t, t′) with **D = H³/4π²**, i.e. d⟨φ²⟩/dt = H³/4π² (Vilenkin–Ford; Linde; Starobinsky), and for tensors
d⟨h_{ij}h^{ij}⟩/dt = 2H³/(π²M̄_P²) = H·P_T. **RECORD CORRECTION:** the record's "graviton channel has no stationary
reduced object / no derived memory time" (tt_worldline) becomes: all response data are stationary; the zero mode
random-walks with a derived diffusion constant. "Infinite memory" is the random walk; "finite memory" is thermal.

## D11 — The conventional KMS defect under slow roll (DERIVED baseline; `v11`)

For a conformal field in the conformal vacuum of a power-law FRW universe a ∝ t^p, ε = 1/p, the exact worldline
Wightman function W = −1/(4π²a(t)a(t′)(Δη − i0)²) has mean-time Wigner spectrum G⁺(ω; t̄) = H̄·g(ω/H̄; ε)
(self-similar), recovering the Planck spectrum as ε → 0. The KMS defect D(ω) = ln[G⁺(ω)/G⁺(−ω)] − 2πω/H̄ is odd,
vanishes identically at ε = 0, and is O(ε) with D/ε → f(ω/H) = 3.17, 6.37, 12.79, 19.21 at ω/H = 0.5, 1, 2, 3:
predominantly an **effective-temperature shift** T_eff ≈ (H̄/2π)(1 − 1.01ε) (the Hubble rate about one Hubble time
after the mean time), plus a **convention-free** frequency dependence [β_eff(3H) − β_eff(0.5H)]/(εβ) ≈ 0.010 — no
single temperature fits. This is conventional QFT in FRW and fixes the O(ε) violation completely for free conformal
fields; it is the baseline against which any non-KMS GRUT residual must be reported (`07` §6, `05` O5).

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
