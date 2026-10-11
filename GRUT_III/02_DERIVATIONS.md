# 02 — Derivations

Units ħ = c = k_B = 1 unless shown. Scripts: `calc/v1`–`v6`; outputs in `RESULTS.md`.

## D1 — The medium's state is thermal at the horizon temperature (THEOREM; `v1` §1.1–1.3)

For a conformally coupled scalar in the Bunch–Davies state on dS₄ (flat slicing, a = −1/(Hη)), the Wightman
function at two events on one comoving worldline is

    g(τ) = −(H²/16π²) / sinh²(H(τ − iε)/2),        τ = proper-time separation,

because a(η)a(η')(η−η')² = 4 sinh²(HΔt/2)/H² exactly (`v1` 1.1). The image-sum identity
Σₙ (z + inβ)⁻² = (π/β)²/sinh²(πz/β) (`v1` 1.3) shows this *is* the flat-space thermal correlator at β = 2π/H.
Its Fourier transform, obtained by shifting the contour to Im τ = −β/2 where the integrand is smooth,

    G⁺(ω) = (ω/2π) / (1 − e^{−2πω/H}),    S(ω) ≡ G⁺(ω) + G⁺(−ω) = (ω/2π) coth(πω/H),    S(0) = H/2π²,

verified by quadrature to 1e-17 (`v1` 1.2). *Record correction:* `calc/RESULTS_worldline_reduction.md` reports
0.03393 at ω = 0.1H and calls the spectrum "UNRESOLVED" after two closed forms fail; the exact one-sided value is
G⁺(0.1H) = 0.034116 (the record's quadrature undershoots by 0.5%), and the two rejected closed forms were
symmetrised objects. The floor S(0) = H/2π² is the classical Nyquist floor 2T·(Im χ/ω)(0); it is not
"zero memory" — the kernel decays as e^{−H|τ|}.

## D2 — The noise kernel and the memory, in closed form (DERIVED; `v1` §1.4–1.5)

Explicit Wick contraction of the improved (ξ = 1/6) stress tensor gives, for one real conformal scalar in flat
Euclidean space along the time axis (and tracelessness ⟨T_μμ T₀₀⟩ = 0 at generic points as a check),

    ⟨T₀₀(t) T₀₀(0)⟩_c = (3/4) C_T / t⁸ = 1/(4π⁴ t⁸),     C_T = 1/(3π⁴).

Mapping with the identity of D1 (orthonormal-frame components T̂ = a⁻⁴ T_flat on each leg), the de Sitter
geodesic observer sees

    ⟨T̂₀₀(t) T̂₀₀(0)⟩_c = (1/4π⁴) [ H / (2 sinh(H(t − iε)/2)) ]⁸      — exactly stationary.

Contour shift plus the recursion F_{2n+2} = F_{2n}(4n²a² + ω²)/(2n(2n+1)a²) for ∫e^{iωu} sech^{2n}(au) du gives
the Wightman spectrum in closed form:

    G⁺_T(ω) = ω (ω² + H²)(ω² + 4H²)(ω² + 9H²) / ( 10080 π³ (1 − e^{−2πω/H}) ).

Checks (`v1` 1.5, all to ≤1e-16 unless stated): quadrature agreement at seven frequencies; KMS detailed balance
G⁺_T(ω)/G⁺_T(−ω) = e^{2πω/H}; flat limit ω⁷/(10080π³) (ratio 1.0004 at ω = 200H); classical floor
G⁺_T(0) = H⁷/(560π⁴); time-domain decay rate 4H (4.0002 measured at τ = 10–12/H).

**Reading.** (i) Occupation is Planckian at T = H/2π: *derived*. (ii) The dissipative part
G⁺_T(ω) − G⁺_T(−ω) ∝ ω(ω²+H²)(ω²+4H²)(ω²+9H²) is an **odd polynomial**: the retarded response at a point is a
finite sum of δ-function derivatives — *no memory in the response*. This is Huygens' principle for a conformal
medium and it holds in flat space too (ω⁷). (iii) All memory is in the **noise**, with correlation time
1/(4H): the Planck denominator has zeros at ω = inH, the polynomial cancels n = 1, 2, 3, and the first surviving
pole at 4iH sets the decay. This is the "Matsubara-ladder inheritance" the record fenced, exhibited. (iv) For a
spatially extended (k ≠ 0, or homogeneous k = 0) perturbation the response integrates the light cone and acquires
the power-law tails the record's T3 chain found (48/(πΔ⁵), −4/(πΔ³), 2/(πΔ)); those tails are the light cone
reaching to the horizon, not a property of a medium "with" memory. The record's internal tension — "finite
memory" (GRUT_I_What_Survived) versus "dS forces infinite scale-free memory" (PROGRAM_FREEZE) — is resolved:
both were true of different objects (noise at a point vs response of a homogeneous mode).

## D3 — Einstein's equation as the medium's equation of state (DERIVED from P1+P2; `v2`)

At a point p take any null k^a and the local Rindler horizon through p with boost Killing field χ^a = −κλk^a.
Raychaudhuri with θ(p) = σ(p) = 0 gives θ = −λ R_ab k^a k^b + O(λ²), hence a horizon-area change
δA = −∫λ R_kk dλ dA per unit area; the heat flux is δQ = ∫ T_ab χ^a dΣ^b = −κ ∫ λ T_kk dλ dA; the KMS temperature
of P1 is T = κ/2π. Clausius δQ = T·η δA gives T_kk = (η/2π) R_kk for every null k (`v2` 2.1). A symmetric
tensor whose null-cone values vanish is proportional to the metric (`v2` 2.2: SVD nullspace one-dimensional,
spanned by η_ab), so R_ab − (2π/η) T_ab = f g_ab; the contracted Bianchi identity with ∇T = 0 forces
f = R/2 + Λ (`v2` 2.3, FRW check). With η = 1/4G:

    G_ab + Λ g_ab = 8πG T_ab,     Λ = integration constant.

GR is hosted. Λ is *not a response of the medium*: it is the one constant the equation of state leaves free,
which is why it is listed as a datum and why GRUT III predicts w ≡ −1 (see 03).

## D4 — Transport coefficients and the Planckian memory time (THEOREM/DERIVED; `v3`)

The medium's Markov-limit response to shear is the membrane paradigm's: η = c³/(16πG) = 8.03×10³³ Pa·s,
entropy density s = k_B c³/(4Għ) = 1.32×10⁴⁶ J K⁻¹ m⁻² per unit area, so **η/s = ħ/(4πk_B)** (KSS value): the
vacuum medium is a perfect Planckian fluid, and in GR this dissipation is exactly the shear term of Raychaudhuri
(tidal heating). The thermal memory time τ = ħ/(2πk_B T) equals 1/κ in every KMS frame: 1/H₀ = 4.58×10¹⁷ s for
the cosmological horizon, 4GM/c³ = 19.7 µs for a solar-mass black hole, c/a for Rindler. Static-patch identities
(exact rationals): E = ρ_Λ V = 1/(2GH) = TS, so the medium's **free energy vanishes**, F = 0; and
|E_Komar| = 1/(GH) = ½ N k_B T with N = A/ℓ_P² (holographic equipartition, Padmanabhan). These are the
"quotient values" of the gravitational medium, and they are functions of (ħ, c, G, k_B) alone.

## D5 — The matter-sector absorptive kernel from the spectrum (DERIVED; `v4`)

With g = η + κh and coupling −(κ/2)h_μν T^μν, unitarity gives |Im Π(s)| = ½ Σ_f S_f ∫dΦ₂ |M|², averaged over the
five spin-2 polarizations with Σ_λ ε_μν ε*_ρσ = ½(P_μρ P_νσ + P_μσ P_νρ) − ⅓P_μν P_ρσ. Explicit amplitudes
(`v4`): real scalar Σ_λ|ε·T|² = 8p⁴/3 exactly (threshold law β⁵); Dirac fermion, from the trace
Tr[(p̸₁+m)V^{μν}(p̸₂−m)V̄^{ρσ}] with V^{μν} = ¼[γ^μ(p₁−p₂)^ν + (μ↔ν)]: massless ratio to the real scalar
**6.000000000**, threshold law **β³(1 + 8r/3)** reproduced to 2e-16 (the Han–Lykken–Zhang form, from scratch);
massless vector (photon pair, physical polarizations, Maxwell T^μν): ratio **12.000000000**. Hence

    Im Π(s) = c_eff(s) · G s²,     c_eff(s) = (1/120) Σᵢ wᵢ Fᵢ(4mᵢ²/s),

with w = 1 (real scalar), 3 (Weyl), 6 (Dirac), 12 (vector) — the Weyl-anomaly weights — and F the threshold
functions (β⁵, β³(1+8r/3), β·[polynomial]). The massless limit is the anomaly coefficient itself:
Im Π = c G s², c_SM = 283/120 = 2.358 (45 Weyl). The Standard-Model staircase (`v4` 4.4) rises from
120c_eff = 12 (photon only, below 2m_ν) through 21 (below 2m_e), 27 (below 2m_μ), to 283 above 2m_t.

**Consequences.** The edge locations are 2mᵢ, the edge exponents are 5/2, 3/2, 1/2 by spin, the weights are
fixed by statistics: *this is the quotient object of SCOUT_0's ZOOM_OUT_03, computed*. The Källén–Lehmann
time-domain kernel K_R(t) = θ(t)∫ds ρ(s) sin(√s t)/√s then has power-law memory with oscillations at the
thresholds: t^{−7/2}cos(2mt) (scalars), t^{−5/2} (fermions), t^{−3/2} (vectors) — finite memory with derived
structure, not a single pole. On shell, Im Π(0) = 0: a free gravitational wave is not absorbed by the vacuum in
flat space (massless quanta cannot pair-produce). The FRW remainder is Γ ∼ c_eff t_P² H₀⁴/2ω ≈ 5×10⁻¹⁶² s⁻¹ at
100 Hz: the record's "invisible" verdicts are reproduced and *explained* (the medium's leading response is D3,
not the loop).

## D6 — The anomaly sector (DERIVED; `v5`)

a = (N_s + 11N_D + 62N_V)/360, c = (N_s + 6N_D + 12N_V)/120; SM: a = 2.765, c = 2.358 (C_T = 40c/π⁴ = 0.968).
Checks: conformal scalar ⟨T^μ_μ⟩_dS = −3aH⁴/2π² = −H⁴/240π²; C_T(scalar) = 1/3π⁴ as in D2. The self-consistent
de Sitter solution driven by the anomaly alone (Starobinsky 1980) has H² = πM_P²/a → H = 1.07 M_P at SM content:
Planckian, excluded by r < 0.036 (which would require a ≳ 2×10¹¹, i.e. ∼10¹² species). GRUT III therefore does
*not* claim anomaly-driven inflation; the R² scale is supplied, as Starobinsky concluded. (This retires the
record's conformalon/second-mode line definitively: w = +1/3 was the symptom, Planckian H is the cause.)

## D7 — The outcome rule (THEOREM given P5)

Decomposition-independence of outcome statistics (SCOUT_0 P-15's characterisation of the martingale branch
coordinate) is non-contextuality; Gleason's theorem then gives Tr(ρP) in dimension ≥ 3, and the POVM extension
(Busch) in all dimensions. Nothing else about measurement is claimed: GRUT III does not derive P5.

## D8 — Where the arrow enters (THEOREM; unchanged from the record, located)

The direction of dissipation is sign(β). In P1, β > 0 is the statement that the *future* horizon is traced
out, equivalently that the state is Hadamard with positive frequencies toward the future. The record's result —
existence of dissipation intrinsic, direction imported — stands; GRUT III adds only that the import is a single
geometric choice (P4), not a property of a bath.
