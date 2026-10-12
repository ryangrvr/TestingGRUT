# 09 — Two first-order temperatures in FRW, and which one GRUT III's principle calls a temperature

**Question (reviewer, 2026-10-11).** In a slowly evolving FRW universe the comoving detector reads
T_det = (ä/ȧ)/2π = (H/2π)(1 − ε) (D11), while the apparent horizon carries Hayward's dynamical surface gravity
κ = −(1/R_A)(1 − Ṙ_A/(2HR_A)), i.e. T_KH = (H/2π)(1 − ε/2) in flat FRW (Hayward 1998; Cai–Kim 2005; Cai–Cao–Hu 2009).
They agree in de Sitter and differ at O(ε). A principle that makes the bath the horizon's causal complement and gets
gravity from horizon thermodynamics must say which is the bath temperature: (a) equal, (b) different by a stated rule,
or (c) silent.

**Answer: (b′), with the temperature derived (`v11`, `v14`) and the exclusion of the horizon labels following from P1's definition.** The two numbers are different objects, and P1 says which is a
temperature.

## 1. What each quantity is (`calc/v12`)

- **T_det = (ä/ȧ)/2π** is the KMS parameter of the *state*: for a conformal field in the conformal vacuum, the comoving
  worldline's two-point functions are, to first adiabatic order, exactly those of a thermal state at this temperature
  (D11, theorem; any a(t) with ε and η first-order small — not the present universe, where η₀ ≈ −2). It is the KMS parameter
  of the fixed-mean-time (Wigner) spectrum at first order. An Unruh–DeWitt detector's *instantaneous* rate is one-sided in time
  and differs at O(ε) by a frequency-dependent term from Ḣ (audit, `11`; D13 "Convention"), so "what the detector reads" is
  convention-laden at this order; the Wigner convention is the record's.
- **T_KH = |κ|/2π** is the surface gravity of the apparent horizon with Kodama normalisation. For ε ≠ 0 the apparent
  horizon is **not a Killing horizon** (the Kodama vector's Lie derivative of the metric is −a²rḢ; `v12` 12.2, now
  computed). P1 assigns a temperature only where the modular flow of the restricted state is geometric (the list in `01`:
  wedges, Killing horizons, conformal double cones); the apparent-horizon region for ε ≠ 0 is in none of these classes, so
  for ε ≠ 0 **P1 assigns no temperature to the region inside R_A** — an application of P1's definition, not a computation.
  T_KH is a geometric label *for the comoving bath*. [Correction, `13`: it is more than a label at the horizon itself.
  Kurpicz, Pinamonti & Verch (Lett. Math. Phys. 111 (2021) 110) prove that Hadamard states become KMS with respect to the
  Kodama flow at the Hayward surface gravity in a near-horizon scaling limit, for spherically symmetric outer trapping
  horizons. Whether this extends to the cosmological apparent horizon is open. So T_KH is the near-horizon scaling-limit
  temperature and T_det the comoving detector's bath temperature: different observables, both conventional.] The literature's own candidate, T₀ = 1/(2πR_A) = H/2π (Cai–Kim 2005; Cai–Cao–Hu 2009;
  adopted by Tian–Booth 2014), is likewise a geometric label; revision 2.3 misattributed T_KH to the first two (audit, `11`).
- The first-order mismatch is half the apparent horizon's expansion rate relative to T_KH, T_det/T_KH − 1 = −Ṙ_A/2 + O(ε²),
  and the full rate relative to T₀, T_det/T₀ − 1 = −Ṙ_A = −ε exactly. Equal iff ε = 0.

## 2. Why no heat-flow law is implied

GRUT III derives gravity in D3 from **local Rindler horizons** (Jacobson), whose boost normalisation cancels; no
cosmological temperature enters, and the Friedmann equations follow from the resulting Einstein equations. The
apparent-horizon "first laws" are then **identities** of those equations, not inputs (`v12` 12.3, symbolic, flat FRW
with a perfect fluid):
- Cai–Kim: the energy flux through the apparent horizon in dt equals T₀ dS with T₀ = 1/(2πR_A) — identity;
- Akbar–Cai / Hayward: dE_MS = (κ/2π) dS + W dV with the **signed** κ, E_MS = ρV, W = (ρ − p)/2 — identity.
The *signed* κ/2π = −T_KH is the bookkeeping temperature that closes Hayward's unified first law (with +T_KH it does not
close); Tian–Booth (2014) compare T₀ with Hayward's κ/2π and adopt T₀, their positive-heat-out convention fixing the sign of
the Gibbs relation, not the temperature (corrected after the audit, `11`). Option (b) in the reviewer's list — a
non-equilibrium horizon–bath heat-flow law at O(ε) — would be needed only if T_KH were a bath temperature. It is not.

## 3. What is testable, and what it is not

The definite first-order statement is: **for a conformal field in the conformal vacuum, the comoving worldline's
fixed-mean-time two-point spectrum is KMS, to first adiabatic order, at (ä/ȧ)/2π, not at either apparent-horizon label**
(and, by `v14`, the same holds for the massless sector's ω ≠ 0 spectrum). This is conventional quantum field theory in FRW,
checkable by anyone, and it is not a GRUT discriminator — it is GR + QFT read correctly. It does not apply to the present
universe (ε₀ = 0.47, η₀ = −2; D11); the order of magnitude of any such temperature today, H₀/2π ∼ 10⁻³⁰ K, is unobservable
regardless. Its value for this program is that it closes a door: the O(ε) window
the reviewer identified is not a non-thermal window for the conformal sector; the first conventional non-thermality
remains O(ε²) (D11), and the bath temperature's meaning in FRW is now fixed by P1 rather than left to a choice between
literature conventions.

## 4. The massless sector, done (`v14`, D13)

The caution below was followed: the worldline function was built with the power-law t(η), from the de Sitter function of
index ν = 3/2 + ε/(1−ε) times the ratio of scale factors, and checked against a direct mode integral. Result: at first order
the fixed-mean-time (Wigner) spectrum of the massless minimally coupled sector along a comoving worldline is, at ω ≠ 0, KMS
at **T_det = (ä/ȧ)/2π + O(ε²)**, the same temperature as the conformal sector's in the same convention; the Kodama–Hayward
rate (H/2π)(1 − ε/2) is excluded at O(ε). The computation contains no coarse-graining scale (contrast `v13`). So the answer
(b′) of §3 is now established for both sectors at first order in the Wigner convention (constant ε). Caveat from the audit:
a linearly coupled Unruh–DeWitt detector's instantaneous rate for this sector is non-thermal already in exact de Sitter
(the secular zero mode gives it an H³/4π²E² term; Garbrecht–Prokopec), so "the detector temperature" of the massless sector
is a statement about the Wigner spectrum or a derivative-coupled detector, not about the linear instantaneous rate. The first-order effect in the massless sector is a shape change of the spectral density
(r(0) = 2 − 2 ln 2, r(∞) = 0), not a temperature change.

## 5. Caveats and priority (historical; the computation is now `v14`)

- T_det was derived for conformal fields in the conformal vacuum. For the massless minimally coupled sector (and
  gravitons) the state is not conformally related to Minkowski and the result is not yet known. **Tool (`v12` 12.4):**
  in power-law FRW, a ∝ (−η)^{1/2−ν}, the massless mode equation v″ + (k² − (ν² − ¼)/η²)v = 0 is identical to that of a
  field of mass m in de Sitter with ν² = 9/4 − m²/H², with ν = 3/2 + ε/(1−ε) (the dS analogue is tachyonic, so the IR
  needs the usual care; the velocity correlator is IR-finite). The massless worldline correlator in power-law FRW is
  therefore the dS massive one times a ratio of scale factors, and the first-order temperature of its stationary part is
  computable from known dS functions — **with one caution (reviewer):** the lemma identifies the *mode equations in
  conformal time*; a comoving detector's proper time maps to η differently in power-law FRW than in de Sitter, so the
  worldline spectra must be built with the power-law t(η), not by transplanting de Sitter Wightman functions. And the
  power-law massless modes are textbook (Abbott–Wise 1984; Lucchin–Matarrese 1985), with comoving-detector responses for
  massless fields in dS in Garbrecht–Prokopec (2004): if the item is largely published, effort goes to `05` O5 item 3
  instead.
- Priority: the detector-vs-horizon distinction and the sign conventions are discussed in Tian & Booth,
  arXiv:1411.6547; the Kodama–Hayward temperature in FRW is Cai–Kim (2005), Akbar–Cai (2007), Cai–Cao–Hu (2009);
  the reviewer's arXiv:0812.3006 is Jiang–Ke–Peng–Feng (tunnelling and the unified first law), which uses the same
  κ. The first-order detector temperature (ä/ȧ)/2π for conformal fields is derived here; a literature check for it is
  owed (`05` O7).
