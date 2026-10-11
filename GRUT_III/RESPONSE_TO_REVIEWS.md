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
convention-free frequency dependence of about 1% of ε. [**Superseded — see R15:** the coefficient is exactly 1 and the
frequency dependence is O(ε²); the sentence is kept as the record of what revision 2 claimed.] Conventional theory fixes it completely for free conformal
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


## R17 — Addendum (fourth round, 2026-10-11): revision 2.2 accepted; two cautions

Accepted as stated: (1) the `v12` lemma identifies mode equations in conformal time only; comoving worldline spectra in
power-law FRW require the power-law t(η) and cannot borrow de Sitter Wightman functions (`09` §4 amended); (2) item 1 of
`05` O5 is likely textbook (Abbott–Wise 1984; Lucchin–Matarrese 1985; Garbrecht–Prokopec 2004), so the literature is
checked first and effort moves to item 3 — whether the secular-infrared resummations disagree — which is the only item
that can open a window rather than close one. `10_SECULAR_IR_SURVEY.md` records that check.

## R18 — Addendum (fifth round, 2026-10-11): the item-3 survey and the decisive `v13` computation

The fourth-round reprioritisation (R17) was followed. **Survey (`10`):** the scalar-sector resummations agree at leading log
(one Fokker–Planck equation and equilibrium from stochastic, diagrammatic, wavefunction, static-patch, soft-EFT and Euclidean
methods) and at NLO; no opening. The graviton-sector disagreement is real but is a dispute over which quantities are
gauge-invariant observables in perturbative quantum gravity on dS, not between competing resummations; none of P0–P5 speaks to
it. Item 1 is published except for the O(ε) worldline KMS defect of the massless sector (conventional, small). **`v13`:** the
published O(ε) Fokker–Planck coefficients of the massless zero mode in power-law FRW (Cho–Kim–Kitamoto 2015, eqs. 3.14–3.16,
confirmed at equation level) define T_zm = (3−ε)A₀H⁴/2 = 1/Vol(S⁴) + O(ε). Its O(ε) coefficient 2/3 − 2γ_E − 2 ln 2 = −1.874 is
the Euclidean weight of no geometric horizon temperature: the required rate H_X = H(1 − 0.4685 ε) is neither the detector value
(c = 1), nor Kodama–Hayward (c = 1/2), nor H. Kodama–Hayward is numerically close (diffusion coefficient −1.50 vs −1.54), but
the gap 5/2 − 2γ_E − 2 ln 2 is exact and transcendental — a near-coincidence that the script now checks explicitly so that no
later reader mistakes it for a match. **Consequence:** P1 assigns no temperature for ε ≠ 0 (non-geometric modular flow) and adds
nothing; conventional theory fixes the coefficients. O5 items 1 and 3 close. With `v11` (conformal sector, O(ε²)), `v12`
(horizon first laws are identities) and `v13` (massless sector, O(ε)), every regime computed on this branch is fixed by GR + QFT
and P0–P5 add no selection. `05` O5 now states the fork this leaves: accept GRUT III as the thermodynamic reading of GR + QFT, or
propose a sixth, parameter-free postulate that breaks KMS-by-construction and test it with the `07` instrument against the
baselines fixed here. Owed: independent review of `v12` and `v13` (`08`, O7).

## R19 — Addendum (sixth round, 2026-10-11): `v13`'s inference refuted; restated

The reviewer's point is accepted in full. The stochastic noise is evaluated at k = σaH; the published coefficients are at σ = 1;
`v13` rev. 2 derives from the Bunch–Davies mode function that the white-noise amplitude at general σ is A₀(1)σ^{−2ε/(1−ε)}H³,
so c₁(σ) = 1 − 2γ_E − 2 ln(2σ) and c_X(σ) = 0.4685 + ½ ln σ — σ-independent only in exact dS. The reviewer's table is reproduced
(c_X = 0 at σ = 0.39, 1/2 at σ = 1.07, −0.68 at σ = 0.1, −1.83 at σ = 0.01) and the "transcendental gap" of rev. 1 is γ_E + ln(2σ)
at σ = 1. The O(ε) zero-mode coefficient is not a physical invariant and cannot test horizon temperatures; rev. 1's inference
is withdrawn, its arithmetic stands. **What survives:** items 1 and 3 close on P1 (no geometric modular flow → no P1
temperature) and on conventional QFT fixing the σ-independent correlators — the grounds the reviewer named. **On the suggested
σ-independent quantity:** the equal-time correlator at fixed physical separation is σ-free but, for a tilted scale-free
spectrum, comparing it with an instantaneous H still requires choosing the separation, which shifts the O(ε) coefficient by
the logarithm of that choice; the quantity without that ambiguity is frequency-resolved — the comoving detector's Wigner
spectrum at ω ≠ 0 and its KMS ratio at O(ε), exactly what `v11` computed for the conformal sector. That is `v14` (via the `v12`
12.4 lemma and the dS massive two-point function expanded around ν = 3/2). **Governance:** agreed that a sixth postulate
breaking KMS is a Stage-3 𝒦 card, to be written as Card 1 under the frozen charter rather than as new governance; `05` O5 says
so. **Owed verifications:** a workflow of independent re-derivation and adversarial-audit agents over `v8`, the general-a(t) and
ΛCDM parts of `v11`, `v12` and `v13` rev. 2 was run this round (`11_VERIFICATION_ROUND.md`).

## R20 — Addendum (sixth round, 2026-10-11): the σ-free computation (`v14`)

Done the same round. The massless minimally coupled sector's comoving two-point function in power-law FRW is, by the `v12`
lemma, the Bunch–Davies de Sitter function of index 3/2 + ε/(1−ε) times a ratio of scale factors (checked against a direct
mode integral); expanding around the massless pole gives, at first order, (1 + 3ε/2)·W_mmc^{dS}(Δ; ä/ȧ) + ε(H̄²/8π²)S₂(cosh²(H̄Δ/2))
+ ω = 0 terms, S₂ an explicit dilogarithm, verified by Richardson in ε. Any function of cosh²(H̄Δ/2) is KMS at 2π/H̄ by
midline symmetry, so the first-order KMS defect at T = (ä/ȧ)/2π vanishes: **the massless sector is thermal at the detector
temperature at first order, like the conformal sector; Kodama–Hayward is excluded at O(ε); nothing depends on a
coarse-graining scale.** The first-order effect is a shape change of the spectral density (r(0) = 2 − 2 ln 2, r(∞) = 0: the
Hadamard coefficient is unchanged), and the σ-free zero-mode diffusion coefficient is D = (H̄³/4π²)[1 − (1 + 2 ln 2)ε], which
the stochastic A₀(σ)H³ reproduces only at σ* = e^{1−γ_E} = 1.53 — the exact function calibrates σ for this observable. On the
reviewer's proposed σ-free quantity (the equal-time correlator at fixed physical separation): the frequency-resolved
Wightman function was used instead because the tilted spectrum leaves any single-separation comparison with a logarithm of
the chosen separation; the two are the same information organised differently, and the frequency-domain version has the
KMS structure in it directly. P1 adds nothing to any of this; D13 records it as conventional QFT. Owed: independent review
of `v14` (next round), and a literature check for r(ω) and the −(1 + 2 ln 2) coefficient.
