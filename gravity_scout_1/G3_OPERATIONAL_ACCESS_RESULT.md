# GRAVITY-SCOUT-1 · G3 — OPERATIONAL ACCESS, FINITE TIME AND GRAVITATIONAL RESOLUTION (result)

**Central question (owner).** Does gravity itself constrain or select the observable class / time access / resolution
at which a subsystem description exists? Or does G2 merely relocate subsystem structure into still-supplied A_interface,
A_resolution and A_time?

**Evidence:**
- `G3_ACCESS_LEDGER.md` (access vectors, the per-setting ledger, source grades);
- `g3/g3_timeband.py` + `.log` (G3-2 conditioning illustration; an independent code path, not an independent reviewer).

**Not a d gate.** No length scale is sought.

## Verdict (per setting — not collapsed, per the owner rule)

| setting | terminal |
|---|---|
| **AdS, boundary time band** (CPR protocol; WdW perturbative holography) | **ACCESS CONSTRAINED + CONDITIONAL ON ACCESS DATA.** Gravity excludes an interior subsystem hidden from boundary time-band observers, *given* declared ε_t, Λ, operator class and a unique vacuum. **Gravity does not determine ε_t, Λ or the operator class** |
| **AdS, robustness** (G3-2 illustration) | **EXACT ALGEBRAIC ACCESS ≠ ROBUST OPERATIONAL ACCESS.** Every ε_t > 0 works algebraically, but the conditioning grows like ε_t^{−(4N−2)}. The robust time band grows with the energy cutoff. This is **kinematic**, not gravitational |
| **Asymptotically flat, finite retarded time** (BCHW) | **A_time CONSTRAINED-NONUNIQUE.** The Bondi mass is not in the algebra of any finite portion of null infinity (under the entropy-bound premises). A class of finite-time measurements is forbidden |
| **Coarse-graining** (Raju) | **COARSE-GRAINING SUPPLIED.** No gravity-derived principle uniquely selects the coarse-grained algebra; the result is generally state-dependent |
| **Perturbative order** | **PERTURBATIVE ORDER SUPPLIED.** No canonical nested hierarchy is supplied by the sources used |

- **A_time ↔ A_interface coupling: ASYMPTOTICS-DEPENDENT.**
  - In AdS, boundary energy data is usable in any time band.
  - In flat space, even the gravitational charge labels (Bondi mass) are not finite-time accessible.
- **CONDITIONALLY SELECTED: none. TRUE COMPRESSION: 0. Empirical payoff: none.**

## G3-1 AdS finite-time boundary protocol (Chowdhury–Papadoulaki–Raju, arXiv:2008.01740)

**Owner-stated source-text facts** (recorded; **not re-read here**, source located only):
- observers remain near the boundary;
- their detectors operate only in a small time band [0, ε];
- they use simple low-energy unitaries and low-energy measurements;
- gravitational backreaction and vacuum entanglement are essential;
- within the stated low-energy Hilbert space, the protocol identifies the bulk state;
- the same protocol fails in nongravitational QFT.

**Priced:**
- global AdS geometry;
- a unique vacuum;
- a low-energy cutoff Λ with Λ ≪ M_Pl;
- the time band ε;
- the allowed simple-operator degree;
- projective / energy measurements;
- semiclassical gravity.

**Does gravity determine ε, Λ or the operator class? No.**
- ε is arbitrary: the protocol works for a small band, not a selected one.
- Λ is declared (Λ ≪ M_Pl is a premise, not an output).
- "Simple low-energy operators" is a declared class.

The exact reconstruction is therefore **CONDITIONAL ON ACCESS DATA**: a relocation into A_time / A_resolution /
A_interface, not a selection.

**What gravity does remove.** Given those access data, an interior subsystem whose low-energy state is hidden from every
boundary time band is excluded. Without gravity the same protocol fails, so the combination

> (independent hidden interior) + (boundary time-band access, low-energy sector)

is admissible at G = 0 and **excluded** at G ≠ 0 at the stated scope.

**Terminal: ACCESS CONSTRAINED** (a combination removed from the admissible residual class), **conditional on the access
data.**

## G3-2 Time-band shrinking / conditioning control (ILLUSTRATION ONLY)

**Model.** The l = 0 sector of a free scalar in global AdS₄, Δ = 3, truncated to N modes with ω_n = 3 + 2n.
- Extracting a_m from the boundary operator on [0, ε_t] is a 2N-exponential moment problem.
- Its minimum-norm solution has ‖f‖² = (G⁻¹)_mm.
- At ε_t = π, G = π·I (perfect conditioning).
- Mode normalizations c_n are set to 1. They rescale the result polynomially and do not change the ε_t-dependence.
- This is the standard AdS spectrum, **not** a transcription of the source's own truncated construction, which was not
  re-read here.

**Results** (160-digit arithmetic):
- **The condition number scales as ε_t^{−(4N−2)}.** It reaches 1.3 × 10⁸ (N = 2) and 7.2 × 10⁴⁵ (N = 8) at ε_t = π/32.
- **The lowest-mode coefficient norm** reaches 2.0 × 10⁴ and 2.8 × 10²² there (relative to the full-period value).
- **The fixed-robustness band grows with the cutoff.** The smallest ε_t/π with cond ≤ 10⁶ is
  0.070 / 0.196 / 0.314 / 0.476 / 0.580 for N = 2 / 3 / 4 / 6 / 8.

**Establishes only:** **EXACT ALGEBRAIC ACCESS ≠ ROBUST OPERATIONAL ACCESS.**
- The states stay mathematically reachable for every ε_t > 0, but reconstruction becomes ill-conditioned.
- A_time and Λ (A_resolution) are **coupled kinematically**. The construction contains no gravity, and with gravity off
  it is identical (G3-9).

**Not** identified with a physical measurement precision (firewall).

## G3-3 Asymptotically flat finite-time control (Bousso–Chandrasekaran–Halpern–Wall, arXiv:1709.08632)

**Abstract-verified statement.**
- The Bondi mass cannot be observed in finite retarded time, so it is not contained in the algebra on any finite portion
  of null infinity.
- This follows from recently discovered asymptotic entropy bounds.
- Attempts to measure a conserved charge at arbitrarily large radius in fixed retarded time are thwarted by quantum
  fluctuations.

**Priced:**
- asymptotically flat geometry;
- null infinity;
- the asymptotic entropy-bound assumptions (body not re-read);
- the definition of the observable algebra on a finite null-infinity portion.

**Reading.** This is a genuine **gravity-related access constraint**: it forbids a class of finite-time measurements
(finite-window charge readout). It selects no window.

**Terminal: A_time CONSTRAINED-NONUNIQUE (asymptotically flat).**

**Comparison with AdS.** In AdS, boundary energy data enters a protocol that works in an arbitrarily small time band. At
null infinity, the analogous total charge is not available in any finite window. Notably, the **charge labels of G2's
perturbative gravitational splitting are themselves not finite-time accessible at null infinity.**

**Recorded: A_time → A_interface coupling is ASYMPTOTICS-DEPENDENT.** The two geometries are not forced into one answer.

## G3-4 Wheeler–DeWitt perturbative control (Chowdhury–Godet–Papadoulaki–Raju, arXiv:2107.14802)

**Abstract-located.**
- At leading nontrivial order around AdS, the WdW and diffeomorphism constraints force correlations between a component
  of the asymptotic metric and energetic excitations of matter / gravitons.
- This gives a perturbative version of holography.
- **Owner-stated (not re-read):** two states / density matrices agreeing on the boundary for an infinitesimal time
  interval agree in the bulk, at the stated perturbative scope.

**Classification: both.**
- **Independent support** that gravity constrains the observable algebra. The correlation is forced by the gravitational
  constraints themselves, not by a chosen protocol: **CONSTRAINED** (Σ / A_partition: no independent bulk sector hidden
  from boundary data at that order).
- **And RELOCATION.** It consumes the supplied AdS boundary structure, the perturbative order and the state class. It
  selects no time band (any infinitesimal one) and no resolution.

## G3-5 Coarse-graining (Raju, arXiv:2110.05470)

The source distinguishes three things:
- fine-grained boundary completeness;
- possible coarse-grained observable sets;
- a state-dependent coarse-grained entropy / approximate split.

**Is there a gravity-derived principle that uniquely chooses the coarse-grained algebra? None found in the sources used.**

**Terminal: A_resolution / observable class REMAINS SUPPLIED (COARSE-GRAINING SUPPLIED).** This is the key selector test,
and it fails to select.

## G3-6 Perturbative-order test

- Donnelly–Giddings (O(κ)) and WdW (leading nontrivial order) each give **one** observable class at **their** order.
- No source used here supplies a canonically nested sequence 𝒪₁ ⊂ 𝒪₂ ⊂ … converging to the fine-grained boundary algebra.
- None is invented.

**Terminal: PERTURBATIVE ORDER SUPPLIED.**

## G3-7 Source-backed vs toy resolution

| | source-backed access hierarchy | G2 ε toy |
|---|---|---|
| variable | observable algebra, asymptotics (∂), time band ε_t, energy cutoff Λ, perturbative order k | a numerical threshold for discarding small generators |
| physical? | each entry is a physically defined access datum, but supplied | no |
| primary source mapping physical resolution ↔ toy ε? | — | **none found → NO** |

The G2 ε ladder is preserved as a **conceptual illustration only**. The G3-2 ε_t is a physical time band, but its
conditioning number is not a precision.

## G3-8 Selector bar

**Required for CONDITIONALLY SELECTED:** gravitational premises + state / asymptotics → a unique operational observable
algebra / time resolution, **without freely declaring the target cutoff**.

**Not met.**
- Every reconstruction theorem used needs a declared ε_t, Λ, operator class, asymptotic region or order. That is
  **RELOCATION / CONDITIONAL ACCESS**.
- **Earned:**
  - **A_time CONSTRAINED-NONUNIQUE** (asymptotically flat): a class of finite-time measurements is forbidden;
  - **ACCESS CONSTRAINED** (AdS): a hidden-interior combination is excluded, conditional on access data.

## G3-9 Flat / nongravitational controls

| control | result |
|---|---|
| AdS protocol, gravity off | fails (owner-stated source result). The hidden-interior combination is admissible at G = 0 |
| null infinity, gravity off | no Bondi mass; the finite-window statement has no gravity-off analogue here |
| WdW, gravity off | no constraint, so no forced boundary–bulk correlation |
| conditioning toy | identical with gravity off (kinematic) |
| ordinary QFT | AQFT locality / split independence is the comparison; independent interior specification is available under the split premises |

## G3-10 Orientation firewall

The following are **supplied** in every setting used:
- retarded time;
- future null infinity;
- the positive AdS Hamiltonian / unique vacuum;
- the time band's [0, ε] direction.

**No orientation is inferred.** Orientation remains supplied and is deferred.

## Payoff

The seven-criterion bar is unchanged.
- Tomography / state reconstruction, finite-time charge non-measurability and conditioning numbers are structural or
  standard to the quoted literature.
- None is a GRUT-distinctive observable.

**ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS (preserved).** G3 affects the residual ledger only.

## Residual after G3

    D_dyn ⊕ [Σ ⊗ (H_marginals, H_cross, H_epoch)] ⊕ [A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time]

with the gravity refinements:

| component | refinement |
|---|---|
| Σ / A_partition | CONSTRAINED / OBSERVABLE-CLASS-PRICED (G2). AdS: a hidden-interior combination excluded given the access data (G3) |
| A_interface | RELOCATION (dressing, boundary detectors, null-infinity window) |
| A_resolution | SUPPLIED (Λ, operator class, coarse-graining). Kinematically coupled to A_time |
| A_time | asymptotics-dependent. AdS: conditional access. Flat: CONSTRAINED-NONUNIQUE |
| charge sector | NEW SUPPLIED LABEL; not finite-time accessible at null infinity |
| orientation | supplied |

**TRUE COMPRESSION = 0.** All results are auxiliary to canonical GRUT.
