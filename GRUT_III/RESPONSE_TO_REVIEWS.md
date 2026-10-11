# Response to the three external reviews of GRUT III (revision 1 → revision 2)

Three reviews of commit `862e666` arrived on 2026-10-11 (two long-form, one with an independent re-derivation of
the species weights and the fermion threshold law). Every claim was audited against the code and the physics;
where a reviewer's statement was itself checkable it was checked. Verdicts below use **CORRECT / PARTLY /
INCORRECT** for the reviewer's claim, followed by what changed. Nothing was argued away.

| # | Reviewer claim | Verdict | Change |
|---|---|---|---|
| R1 | a₀ "Jacobson fails below cH" argument is unfounded | **CORRECT** | Retracted; `06` demoted to a coincidence appendix |
| R2 | 1/(4H) and 1/H are different quantities; no finite autocorrelation time without smearing | **CORRECT** | Relabelled everywhere; memory now defined only for smeared observables (`07`, v7) |
| R3 | η = c³/16πG has units kg/s (surface viscosity), not Pa·s | **CORRECT** | v3 fixed; D4 reworded (horizon is a 2-d fluid; spacetime is not a bulk fluid) |
| R4 | SM staircase illustrative; √s is not GW frequency | **CORRECT** | v4 relabelled; D5 and 03 reworded |
| R5 | HJW cannot derive Born from causality; check MGM 2019 etc. | **CORRECT** | O3 rewritten; claim withdrawn |
| R6 | w ≡ −1 is ΛCDM's prediction too; extra light fields allowed; no DR3 | **CORRECT** | 03 F1 reframed as a falsifier of the *minimal* framework, shared with ΛCDM; DR2 numbers cited |
| R7 | +16→+5 is reorganisation; tag imports as recovery | **CORRECT** (rule §13.6 not found here; principle adopted) | 01, 04 relabelled; "compression" withdrawn |
| R8 | Type III algebras: no tensor factorisation; P1 needs modular form | **CORRECT** | P1 restated in Tomita–Takesaki terms; "influence functional" explicitly not supplied by P1 |
| R9 | Priority: PNRV 2009, Cho–Hu 2014, SCH 2016, Hu–Verdaguer | **CORRECT** (checked at equation level, `08`; PNRV is ξ = 0 only — the conformal kernel is Cho–Hu's) | "Not on the record" downgraded to "corollary of Cho–Hu (2014); placement" |
| R10 | 0.0339 vs 0.0341 is 0.55%, and only an explanation if the record's precision is coarser | **CORRECT**, and the explanation is now proved | D1 rewritten with the proof |
| R11 | N = ½coth(βω/2)ρ; Δρ is the decisive object | **CORRECT**; adopted with scope notes | New `07`; v8 implements it |
| R12 | Gaussian-smeared free field: χ_L(t) = θ(t) t e^{−t²/4L²}/(8π^{3/2}L³); access geometry makes apparent memory | **CORRECT** (verified analytically; v7 numerically) | Theorem T3 in `07` |
| R13 | G is inserted via P2/P3; nothing yet distinguishes GRUT III from GR+QFT; some sentences overclaim | **CORRECT** | Overclaiming sentences listed below and removed |

## R1 — the a₀ mechanism

The reviewer is right and I was wrong. In Jacobson's construction the local Rindler horizon is an auxiliary
boost horizon through an *arbitrary* point p; its normalisation κ appears in both the heat flux (χ^a = −κλk^a) and
the temperature (κ/2π) and cancels (this is visible in `calc/v2_jacobson_eos.py` §2.1, where κ drops out of the
solved T_kk). Nothing in the construction refers to any observer's physical acceleration; the patch is taken small
compared with the curvature radius, which is always possible. The statement "for a_galaxy ≲ cH the derivation is
not licensed" confused the auxiliary boost with a star's orbit. **Retracted.** `06` now carries only the dimensional
coincidence (which is Milgrom's, 1999) and the a₀ ∝ H(z) discriminator (also proposed before), marked "no
mechanism within GRUT III".

## R2 — memory times

Correct. Revision 2 uses three distinct names and never conflates them: τ_β = β/2π = 1/κ (the thermal scale);
γ_T = 4H (the asymptotic decay *rate* of the pointwise ⟨T̂T̂⟩ correlator, from the first surviving pole after the
polynomial cancels the n = 1, 2, 3 Matsubara rungs); and the *memory of an observable*, defined only for a smeared
observable (f, ω) through the support of its retarded kernel and the decay of its noise kernel (`07`). The
pointwise correlator diverges as t⁻⁸ and has no normalisable autocorrelation time; v7 shows what smearing does.

## R3 — units

Correct: [c³/G] = kg s⁻¹. The quantity is the membrane paradigm's *surface* shear viscosity of a horizon
(a 2-d fluid), and η/s = ħ/4πk_B is the consistent ratio. The sentence "the vacuum medium is a perfect Planckian
fluid" overreached: the *horizon* is described as a 2-d fluid in an established effective description; empty
spacetime is not asserted to be a bulk fluid. v3 and D4 corrected.

## R4 — the Standard-Model staircase

Correct on both counts. (i) The massive-vector threshold polynomial was not computed (leading β used) and the
hadronic region is a proxy: the staircase is **illustrative**, with the derived ingredients being the massless weights
1 : 3 : 12, the scalar β⁵ law and the fermion β³(1 + 8r/3) law (the last two from explicit amplitudes, agreeing
with Han–Lykken–Zhang). (ii) The threshold variable is the timelike invariant mass √s; a free on-shell GW has
s = 0 at every frequency, so *no* species is open for it; the old "at LIGO frequencies" line conflated energy with
invariant mass and is gone.

## R5 — Born rule

Correct. Hughston–Jozsa–Wootters realises ensemble decompositions by measurements on a purification *inside*
standard quantum mechanics; invoking it to derive Born from causality assumes what is to be derived. The
suggestion that P5 might follow from P0 + P1 is withdrawn. Masanes–Galley–Müller (2019) derive the measurement
postulates from the remaining quantum postulates plus operational assumptions (including finite-dimensional state
characterisation); Gleason/Busch, Zurek's envariance and Deutsch–Wallace each trade the postulate for other
assumptions. GRUT III adds nothing to this literature and now says so (O3).

## R6 — w ≡ −1

Correct. Λ as an integration constant is what GR gives; ΛCDM predicts w = −1 too. The honest statement: the
*minimal* GRUT III (P0–P5 with the known field content and no additional light fields) predicts w ≡ −1; a
model-independent ≥5σ detection of evolving dark energy would refute the minimal framework exactly as it would
refute ΛCDM, and would force an additional light field into P0's content (allowed, and priced). It is a falsifier,
not a discriminator. DESI DR2 (2025) prefers w₀w_aCDM over ΛCDM at ~2.8–4.2σ depending on the supernova sample
(numbers checked against the paper in `08`); there is no DR3 cosmology release. "DR3 decides" is removed.

## R7 — the ledger

Correct in substance. The "+16 → +5" sentence read as if eleven inputs had been derived. What happened is: seven
register nodes were replaced by *imported theorems* (Bisognano–Wichmann/Kay–Wald/Gibbons–Hawking thermality;
Jacobson's equation of state; Huygens for conformal fields; the membrane transport coefficients; the dS
identities) — these earn GRUT **zero credit** and are now tagged IMPORTED; two became data (G, ħ) by the charter's
own datum ruling; two became identifications (the split and the flow). GRUT III has *five framework postulates*,
and the word "compression" is withdrawn. The reviewer's Stage-3 charter "§13.6" was not in this repository; the branch-digest agent located it in the
owner's `vsctestinggrut` repository (`PROGRAM/STAGE3/STAGE3_CHARTER.md`, see `08`): §13.1 "relations implied by 𝒮 earn
nothing", NR-12(b)/CT-75 on KMS/FDT premises, §13.6 RECOVERY, and **DEF-16** — KMS relative to a modular flow
constructed from the state itself earns no credit — which is precisely the caveat revision 2 attaches to P1.

## R8 — type III and P1

Correct and useful. Local algebras of relativistic QFT are type III₁; there is no tensor factorisation
H = H_O ⊗ H_{O'} and no reduced density matrix. P1 is restated (01): the observer's medium is the triple
(A(O), ω|_{A(O)}, σ_t^ω) — the local von Neumann algebra, the restricted state, and its Tomita–Takesaki modular
flow; the bath is the commutant A(O)' (= A(O') where Haag duality holds). Thermality in the modular parameter is
then automatic; the *physics* is whether the modular flow is geometric (Bisognano–Wichmann wedges, Kay–Wald
Killing horizons, the dS static patch; Hislop–Longo for double cones in CFTs), which is exactly when a temperature
is a number. P1 supplies no influence functional and no dynamics; those come from P0 (the theory) and P3 (the
probe coupling), after which the Feynman–Vernon functional is *computed*. This is stated explicitly now.

## R9 — priority

Checked at equation level (`08`). Correction to my own first answer: Pérez-Nadal–Roura–Verdaguer (2009) treat only
the **minimally coupled** scalar and do not contain the conformal case; the position-space conformal-field kernel in
dS is **Cho–Hu (2014)**, eqs. (29)+(30), whose Minkowski coefficient at Δs = 0 is exactly 1/(4π⁴Δη⁸) — the same
constant as `v1` — with the conformal rule due to Eftekharzadeh–Bates–Roura–Anderson–Hu (2012); Satin–Cho–Hu (2016)
give the dissipation-kernel rule and the preservation of the FDR. GRUT III's worldline kernel is a two-line corollary
of Cho–Hu and is not new; the proper-time spectrum G⁺_T(ω) and the odd-polynomial reading are not printed in any of
these and are recorded as corollaries. "Memory lives in the noise kernel" is Hu–Verdaguer stochastic gravity.

## R10 — the record correction

The gap is 0.55% (0.034116 vs 0.03393); "within half a percent" is replaced by the number. More importantly the
explanation is now a proof rather than a conjecture: the record's `calc/worldline_reduction.py` *defines* its
"symmetrised" spectrum as S(w) = 2 Re ∫₀^∞ e^{iwτ} g(τ) dτ, and for any Wightman function g(−τ) = g(τ)*, so this
equals ∫_{−∞}^{∞} e^{iwτ} g dτ = G⁺(w) exactly — the one-sided Planck spectrum, not the symmetrised one. The two
"candidate closed forms" it tested were symmetrised objects. The 0.55% residual is its ε-extrapolation error
(ε ∈ {0.04, …, 0.005}).

## R11 — the FDT residual

Correct, and adopted as the central object of revision 2 (`07`, v8). Scope notes added: the identity
N = ½coth(βω/2)ρ holds for stationary KMS states at the two-point (Gaussian) level; the framework is extended with
(i) cumulant residuals Δκ_n for non-Gaussian claims (the record's "C-B cubic drift" was one), (ii) a two-time form
for non-stationary states (the gravitational zero mode, FRW), and (iii) explicit observable dependence (ρ_f for a
test tensor f). The reviewer's trichotomy — change ρ, change the state, or reinterpret — is the organising
statement of what a GRUT claim must be.

## R12 — smearing

Correct; verified analytically (the convolution of two Gaussians of width L is a Gaussian of width √2 L; the radial
integral against δ(t − r)/4πr gives θ(t) t e^{−t²/4L²}/(8π^{3/2}L³)) and numerically in v7, together with the
compact-support statement (no response after 2R/c) and the stress-tensor analogue. Consequence adopted: apparent
memory of order L/c arises from access geometry alone, so access geometry can never be credited as a vacuum
relaxation law. This is now Theorem T3 of `07`.

## R13 — circularity and overclaiming

Correct. G enters through P2 (η = 1/4G) and the coupling through P3 (κ); the matter-loop results are then
gravitational QFT, and nothing in GRUT III distinguishes it from GR + QFT at present. Sentences removed or
rewritten in revision 2: README "the kernel is nowhere supplied, GR is hosted" (→ "no kernel is supplied; GR is
imported via Jacobson with G inserted"); 00 "GR goes from borrowed to hosted" (→ "GR is recovered, at the price
of P2 and the datum G"); 00/README "the frontier question is dissolved" (→ "the frontier question is relocated:
P1 names the bath algebra; it supplies no dynamics"); README "ledger falls from +16 to +5" (→ "five framework
postulates; imports tagged; no compression claimed"); 02 D2 "none of this is supplied" kept only with "and all of it
is standard"; README "Discipline" paragraph "not on the GRUT record" (→ "corollary; placement").

## What revision 2 adds (not just subtracts)

- `07_ACCESS_RESOLVED_RESPONSE.md` + `calc/v8`: the Δρ residual framework, with the record's old kernels run
  through it (they assert Δρ ≠ 0 by fiat — which is what SCOUT_0 meant by "supplied").
- `calc/v7`: memory from access geometry — the reviewer's result, its stress-tensor analogue, and the KMS identity for
  smeared observables.
- `calc/v9`: the induced-gravity/entanglement-entropy matching that makes "where G comes from" (O1) a precise,
  species-by-species consistency condition with its known obstructions.
- `calc/v10`: the gravitational sector's non-stationarity in closed form — a Brownian zero mode with derived diffusion
  constant H³/4π² (tensor: 2H³/π²M̄_P²), all derivative correlators stationary and KMS — closing the record's
  "no stationary reduced object" finding with a formula.
- `08_PRIORITY_AND_SOURCES.md`: what the literature already contains, checked at equation level where the texts
  could be read.

## R14 — Addendum (owner's note, 2026-10-11): KMS by construction means standard by construction

**Claim.** GRUT III imposes KMS via Bisognano–Wichmann; by the FDT boundary any extra noise at the same β is forced
extra dissipation; hence in its own regime GRUT III can only be standard. Distinctive locks need platforms where KMS
fails (the record's "Scope-Q" ruling reached the same point). Compute Δρ in a non-KMS state; first candidate
quasi-de Sitter with ε = −Ḣ/H², where the violation is O(ε) and conventional stochastic gravity predicts it too; a
GRUT claim must be a restriction conventional theory leaves free.

**Verdict: CORRECT, and adopted as the organising statement of revision 2** (`03` first section, `05` O5, `07`).
What was done about it: `v11` computes the conventional O(ε) KMS defect exactly for conformal matter in power-law
inflation (contour-shifted Wigner spectrum of the exact worldline Wightman function). Result: the defect is odd,
vanishes at ε = 0, scales as ε, and is predominantly an effective-temperature shift T_eff ≈ (H/2π)(1 − 1.01ε) with a
convention-free frequency dependence of about 1% of ε. Conventional theory fixes it completely for free conformal
fields; the only sector it leaves open is the interacting gravitational infrared, which `05` O5 now sets out as the
single remaining program. The "Scope-Q / G2-12" ruling was not found in this repository; the convergence is recorded
on the owner's word.

**a₀.** Recorded as Milgrom's (1999) vacuum-effect conjecture with the Narnhofer–Peter–Thirring / Deser–Levin
temperature √(a² + H²)/2π, not as GRUT's (`06`).


## R15 — Addendum (second round of reviews, 2026-10-11): the v11 headline was wrong

**Claim (two independent reviewers; one published with code on `codex/grut-causal-state-review`).** At first order in
ε the exact power-law correlator equals the de Sitter correlator with H → H̄(1 − ε): W_ε = W_dS[1 + ε(x coth(x/2) − 2)]
+ O(ε²) = W_dS(Δ; H(1−ε)) + O(ε²); the first-order detailed-balance relation is thermal at T_eff = (H/2π)(1 − ε); v11's
c ≈ 1.014 and its "1%-of-ε frequency dependence" were finite-ε (O(ε²)) contamination; the second-order coefficient is
W_ε/W_dS(H_eff) = 1 + ε²[x²/4 − (x³/12)coth(x/2)]; Richardson extrapolation gives D/ε → 2πω/H̄ to 3e-5.

**Verdict: CORRECT.** My convergence test accepted a 5% drift between p = 100 and 200 as "converged"; the drift was the
O(ε²) term. `v11` rev. 3 now proves the first-order theorem symbolically for the power law (no 1/q term) and for
general a(t) to second adiabatic order (η = ε̇/Hε counted as first-order small), confirms the reviewers' f₂(x) exactly
and at q = 10⁵ numerically, extrapolates D/ε → 2πω/H̄ (c = 1 to 3e-5), and computes the true O(ε²) frequency-dependent
defect D₂(ω). The physics consequence strengthens the framework's own conclusion: the conformal sector is KMS through
first order, so the FDT boundary holds through O(ε). Also adopted: the late universe is not slow-roll (ε₀ = 0.47); the
exact present-day kernel deviates from any thermal form by up to 19% at Δ = 1.5/H₀ (`v11` D).

**Also accepted from this round:** T6's scope (a harmonic-bath or generalized-free-field representation is a
mathematical representation, not an actual local field theory; Källén–Lehmann needs covariance and positivity) and
T5's scope (the Brownian covariance is asymptotic; imaginary-time periodicity does not settle full KMS analyticity).
Both reworded in `07` and D10. The reviewers' verification of v7, v9, v10 against the smeared-detector result,
Solodukhin's coefficients and the standard D = H³/4π² is recorded in `08`.


## R16 — Addendum (third round, 2026-10-11): which temperature is the bath's?

**Claim.** T_det = (ä/ȧ)/2π = (H/2π)(1 − ε) and the apparent-horizon T_KH = (H/2π)(1 − ε/2) differ at O(ε); GRUT III
must choose (equal / different by a stated rule / silent); this is a first-order window in the non-thermal regime.

**Verdict: the question is correct and answered — option (b′), rule derived (`09`, `v12`).** P1 assigns a temperature
only where the restricted state's modular flow is geometric; the apparent horizon is not a Killing horizon for ε ≠ 0
(Ṙ_A = ε), so P1 assigns it none; T_det is the state's KMS parameter and the detector reading; T_KH is Hayward's
dynamical surface gravity, which makes the horizon first law an identity of the Einstein equations D3 already supplies
(verified symbolically: Cai–Kim with T₀ = 1/2πR_A; Akbar–Cai with the signed κ). No heat-flow law is implied; the O(ε)
mismatch is the apparent horizon's expansion rate. The definite statement — a comoving detector coupled to a conformal
field thermalises at (ä/ȧ)/2π, not at T_KH — is conventional QFT and not a GRUT discriminator. Caveat accepted: this is
for conformal fields in the conformal vacuum; the massless/graviton sector is the next computation, with the
mode-equation identity of `v12` 12.4 as the tool. Attribution note: arXiv:0812.3006 is Jiang–Ke–Peng–Feng; the κ formula
is Hayward's, used by Cai–Kim (2005), Akbar–Cai (2007) and Cai–Cao–Hu (2009); Tian–Booth (2014) is the "temperature
confusion" reference, read at abstract level here.
