# 07 — The access-resolved response: the object a GRUT claim must be about

Revision 2's central structural change. GRUT I asked "what is the vacuum's memory kernel?" and supplied one.
GRUT III (rev. 1) said "memory lives in the noise, with time 1/(4H)" and was rightly told that it had named two
different quantities and that a pointwise correlator has no autocorrelation time. Revision 2 replaces the word
"memory" by a computable object attached to an *observable and a state*, and states what a distinctive GRUT claim
must look like in that language.

## 1. The object

Let O be a causally complete region, A(O) its local algebra (P0), ω the state restricted to A(O) (P1), and
f^{μν} a smooth test tensor supported in O defining the smeared observable F_f = ∫ f^{μν} T_{μν}. Define

    W_f(t)   = ω( F_f(t) F_f(0) )                          (Wightman function along the observer's time),
    G⁺_f(ω)  = ∫ dt e^{iωt} W_f(t),
    ρ_f(ω)   = G⁺_f(ω) − G⁺_f(−ω)                           (commutator spectrum = dissipative response),
    N_f(ω)   = ½ [G⁺_f(ω) + G⁺_f(−ω)]                      (symmetrised noise),
    χ_f(t)   = i θ(t) ω([F_f(t), F_f(0)])                   (retarded kernel; χ̃_f = Hilbert transform of ρ_f).

**T1 (FDT, algebraic).** If ω is KMS at β with respect to the chosen time, then G⁺_f(−ω) = e^{−βω} G⁺_f(ω) and
therefore N_f(ω) = ½ coth(βω/2) ρ_f(ω) for every f. The noise is never an independent datum once the state is
thermal: it is the dissipative response times the Planck factor. (Verified numerically for smeared observables in
v7 §C and for the pointwise conformal correlator in v8 §B.)

**Definition (memory of an observable).** The *response memory* of (f, ω) is the support of χ_f(t); the *noise
memory* is the decay of N_f in the time domain. Both are computed from (f, ω). Neither is a property of "the
vacuum" without f.

**Definition (residual).** For a proposed theory X and the conventional theory (GR + QFT in the same state), with
the same observable f, state ω and preparation,

    Δρ_f(ω) = ρ_f^X(ω) − ρ_f^conv(ω),      ΔN_f = N_f^X − N_f^conv,      Δκ_n,f = κ_n^X − κ_n^conv  (n ≥ 3).

**T2 (trichotomy).** If ω is KMS for both theories at the same β, then ΔN_f = ½coth(βω/2) Δρ_f. Hence a theory
with the same state and the same observables can differ from the conventional one only by (a) a different
dissipative response Δρ_f ≠ 0, (b) a different state (non-KMS, or KMS at a different β, or non-stationary), or
(c) higher cumulants Δκ_n ≠ 0. Anything else is a reinterpretation of the same correlators. (This is the reviewer's
observation, adopted verbatim as the organising principle.)

## 2. Theorems about where memory lives (all imported/standard; computed in `calc/`)

**T3 (access geometry makes memory; v7).** Free massless scalar, G_R = θ(t)δ(t − r)/4πr. For a detector with
normalised Gaussian profile of width L the doubly-smeared retarded kernel is exactly

    χ_L(t) = θ(t) · t · e^{−t²/4L²} / (8 π^{3/2} L³),

a smooth function of support ∼ L/c, although the unsmeared response is on the light cone. For a profile supported in
a ball of radius R the smeared response vanishes for t > 2R/c. For the conformal stress tensor the same smearing
turns the pointwise polynomial ρ_T (no memory) into a kernel of support ∼ L with a Gaussian tail (v7 §B), keeping a
1/t edge at small t (spatial smearing bounds momentum, not energy); an integrable kernel needs temporal smearing too.
*Corollary:* an apparent
memory time can always be manufactured by the observer's access profile and can never be credited to the vacuum.

**T4 (pointwise conformal response is local; v1, v8).** For a conformal field in dS along a geodesic,
ρ_T(ω) = ω(ω² + H²)(ω² + 4H²)(ω² + 9H²)/(10080π³): an odd polynomial, so χ_T(t) is a finite sum of δ-derivatives.
The thermal denominator cancels in ρ_T and reappears only in N_T. The pointwise correlator's asymptotic decay rate
is 4H; it has no normalisable autocorrelation time (t⁻⁸ at short times).

**T5 (the gravitational sector is a stationary thermal process plus a Brownian zero mode; v10).** For the
massless minimally coupled sector (each graviton polarisation), the worldline two-point function is
(H²/8π²)[1/(2 sinh²(HΔt/2)) − ln(2 sinh²(HΔt/2)) + H(t + t') + c]. Every derivative correlator and every commutator
is stationary and 2πi/H-periodic; the non-stationary part is, at large separation, the Brownian covariance
D·min(t, t') with D = H³/4π² per canonically normalised field (asymptotic, not exact at all separations), i.e. d⟨h_{ij}h^{ij}⟩/dt = 2H³/(π²M̄_P²) = H·P_T.
The record's "graviton channel has no stationary reduced object" is therefore: *all* its response data are
stationary; only the zero-mode position random-walks, with a derived diffusion constant (Starobinsky's stochastic
inflation). "Infinite memory" is the random walk; "finite memory" is thermal.


**T6 (realizability: every causal passive kernel is an integrated-out field, at the two-point level; IMPORTED).**
Let χ be a causal linear response kernel whose spectral function ρ(ω) satisfies ω ρ(ω) ≥ 0 (passivity/stability).
Then (i) classically there is a bath of harmonic oscillators, linearly coupled, with spectral density proportional to
ρ, whose elimination reproduces χ exactly (Ford–Kac–Mazur 1965; Caldeira–Leggett 1983; Mori–Zwanzig; the record's
P-17 is the converse: reduced data identify the forcing law, not the bath); (ii) in relativistic QFT a covariant
two-point function with positive spectral measure dσ(s) is that of a *generalized free field*, a superposition of free
fields of mass √s (Källén–Lehmann), so the kernel is literally "an integrated-out field". Scope: two-point/Gaussian
level only; it says nothing about higher cumulants or back-reaction; and the realizing field content is **not unique**
(many local Lagrangians give the same σ — P-17's non-identifiability); and both statements are *mathematical
representations* — a harmonic bath for a passive temporal kernel, a generalized free field for a covariant positive
spectral measure (Källén–Lehmann needs Lorentz covariance and spectral positivity, which an arbitrary causal kernel
does not supply) — not the existence of an actual local field theory with those excitations. Consequence for constructions: eliminating one
specified local field (the constructive branch's CP-05, with its massive retarded kernel and J₁ tail inside the light
cone) exhibits *a* realization of a passive kernel; the kernel cannot *select* that field. "Every causal stable kernel
is an integrated-out field" is therefore proved in this sense and false in the stronger sense (it does not single out a
sector) — which is the distinction the CP-05 review drew.

## 3. The memory dictionary (replaces every earlier use of the word)

| object | what is computed | result | source |
|---|---|---|---|
| pointwise conformal T, dS geodesic | ρ_T, N_T | ρ_T polynomial (local); N_T Planckian; decay rate 4H | v1, v8 |
| smeared observable, width L, any state | χ_f, N_f | support ∼ L/c (response); N_f = ½coth·ρ_f if KMS | v7 |
| homogeneous mode, k → 0 | χ(t) | power-law tails (48/πΔ⁵, …): light cone to the horizon | record T3 |
| matter thresholds | ρ(s) | edges 2mᵢ; exponents 5/2, 3/2, 1/2; Im Π = c_eff G s² | v4 |
| massless/graviton sector, dS | W(t,t') | stationary thermal + Brownian zero mode, D = H³/4π² | v10 |

## 4. What a distinctive GRUT claim must now be

A **residual report** (template implemented in v8 §D): the observable f (with its support), the state ω and
preparation, the conventional ρ_f^conv and N_f^conv, the proposed ρ_f^X (or state, or κ_n), the residual, and the
*mechanism* that produces it. A residual with no mechanism is what SCOUT_0 called "supplied". A mechanism that is
an additional field or operator in P0's content is conventional physics with a price, not GRUT. A mechanism that is
*forced* by P0–P5 and excluded in the conventional comparison would be a GRUT result. **None exists at present.**

Two further requirements, learned from the constructive branch's CP-05 review: (i) **no refitting** — a residual test
must connect *separately measured* quantities (e.g. a static range, a force strength, and a propagating response)
without re-fitting a parameter between them; a relation such as μ² = m²(1 − 3η) is not an independent test while m is
unmeasured; (ii) **the price of a mechanism is displayed** — if the mechanism is an added light field coupled to T_μν,
the existing fifth-force bounds apply (`03`, F1 note) and must be quoted with their confidence levels.

## 5. Applied to the record (v8 §C)

The register's rung-3 noise model (J ∝ ω³ with a cutoff, times coth) and the S5-1 Lorentzian-exponent family are
internally KMS-consistent, but their ρ differs in *shape* from the conformal medium's (low-frequency exponents 3 and
0–1 against 7; a log-log slope residual of −4 for ω³): they assert Δρ ≠ 0 relative to the conventional vacuum
with no mechanism. That is the precise content of the record's own verdict "supplied", now as a number.

## 6. Scope and the one live opening

T1–T2 are two-point, stationary statements. The non-stationary case is handled by two-time kernels (the Wigner
mean-time form the record's dS paper already used), and T5 shows the only non-stationary object in the free
gravitational sector is the zero-mode random walk. Under slow roll the conformal sector stays KMS **through first
order** at T_eff = (ä/ȧ)/2π (D11, corrected after review); the first conventional violation is O(ε²) and frequency
dependent, and today (ε₀ = 0.47) the exact kernel is needed. The apparent-horizon temperature (H/2π)(1 − ε/2) is
not a bath temperature in P1's sense (`09`): it labels a non-Killing horizon and closes Hayward's first law as an identity. The genuinely open regime is the *interacting* gravitational
infrared in de Sitter, where secular growth of loop corrections is unresolved in the literature and where a
non-stationary residual (ΔW(t, t') beyond the Brownian term) would be the first place a "responsive" correction
could appear without being inserted. That is a research problem in quantum gravity, stated here as the one opening;
it is not a GRUT result.
