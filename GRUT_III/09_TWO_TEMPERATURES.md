# 09 — Two first-order temperatures in FRW, and which one GRUT III's principle calls a temperature

**Question (reviewer, 2026-10-11).** In a slowly evolving FRW universe the comoving detector reads
T_det = (ä/ȧ)/2π = (H/2π)(1 − ε) (D11), while the apparent horizon carries Hayward's dynamical surface gravity
κ = −(1/R_A)(1 − Ṙ_A/(2HR_A)), i.e. T_KH = (H/2π)(1 − ε/2) in flat FRW (Hayward 1998; Cai–Kim 2005; Cai–Cao–Hu 2009).
They agree in de Sitter and differ at O(ε). A principle that makes the bath the horizon's causal complement and gets
gravity from horizon thermodynamics must say which is the bath temperature: (a) equal, (b) different by a stated rule,
or (c) silent.

**Answer: (b′), with the rule derived, not stated.** The two numbers are different objects, and P1 says which is a
temperature.

## 1. What each quantity is (`calc/v12`)

- **T_det = (ä/ȧ)/2π** is the KMS parameter of the *state*: for a conformal field in the conformal vacuum, the comoving
  worldline's two-point functions are, to first adiabatic order, exactly those of a thermal state at this temperature
  (D11, theorem; any a(t)). It is what an Unruh–DeWitt detector on the comoving worldline reads at first order.
- **T_KH = |κ|/2π** is the surface gravity of the apparent horizon with Kodama normalisation. For ε ≠ 0 the apparent
  horizon is **not a Killing horizon** (Ṙ_A = ε ≠ 0; `v12` 12.2). P1 assigns a temperature only where the modular flow
  of the restricted state is geometric — Killing horizons and wedges — so for ε ≠ 0 **P1 assigns no temperature to the
  region inside R_A**. T_KH is a geometric label.
- The first-order mismatch is exactly the apparent horizon's expansion rate: T_det/T_KH = 1 − ε/2 + O(ε²) with
  Ṙ_A = ε. Equal iff ε = 0.

## 2. Why no heat-flow law is implied

GRUT III derives gravity in D3 from **local Rindler horizons** (Jacobson), whose boost normalisation cancels; no
cosmological temperature enters, and the Friedmann equations follow from the resulting Einstein equations. The
apparent-horizon "first laws" are then **identities** of those equations, not inputs (`v12` 12.3, symbolic, flat FRW
with a perfect fluid):
- Cai–Kim: the energy flux through the apparent horizon in dt equals T₀ dS with T₀ = 1/(2πR_A) — identity;
- Akbar–Cai / Hayward: dE_MS = (κ/2π) dS + W dV with the **signed** κ, E_MS = ρV, W = (ρ − p)/2 — identity.
T_KH is the bookkeeping temperature that closes Hayward's unified first law; its sign convention is the "temperature
confusion" that Tian–Booth (2014) resolve by a positive-heat-out convention. Option (b) in the reviewer's list — a
non-equilibrium horizon–bath heat-flow law at O(ε) — would be needed only if T_KH were a bath temperature. It is not.

## 3. What is testable, and what it is not

The definite first-order statement is: **a comoving detector coupled to a conformal field in the conformal vacuum
thermalises, to first order, at (ä/ȧ)/2π, not at the apparent-horizon temperature.** This is conventional quantum
field theory in FRW, checkable by anyone, and it is not a GRUT discriminator — it is GR + QFT read correctly. Its
magnitude today (T ∼ 10⁻³⁰ K) is unobservable. Its value for this program is that it closes a door: the O(ε) window
the reviewer identified is not a non-thermal window for the conformal sector; the first conventional non-thermality
remains O(ε²) (D11), and the bath temperature's meaning in FRW is now fixed by P1 rather than left to a choice between
literature conventions.

## 4. The massless sector, done (`v14`, D13)

The caution below was followed: the worldline function was built with the power-law t(η), from the de Sitter function of
index ν = 3/2 + ε/(1−ε) times the ratio of scale factors, and checked against a direct mode integral. Result: at first order
the massless minimally coupled sector seen by a comoving detector is KMS at **T_det = (ä/ȧ)/2π + O(ε²)**, the same
temperature as the conformal sector; the Kodama–Hayward rate (H/2π)(1 − ε/2) is excluded at O(ε). The computation contains
no coarse-graining scale (contrast `v13`). So the answer (b′) of §3 is now established for both sectors at first order,
not only for the conformal one. The first-order effect in the massless sector is a shape change of the spectral density
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
