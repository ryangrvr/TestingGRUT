# GRAVITY-SCOUT-1 · G3 — OPERATIONAL ACCESS, FINITE TIME AND GRAVITATIONAL RESOLUTION (result)

> **Repaired by GRAVITY REPAIR 03 (GR3-01 … GR3-06; `GRAVITY_CORRECTION_LEDGER.md`).** G3 is accepted provisionally
> after that repair. Where wording differs, the ledger takes precedence. Script and log outputs are kept as emitted.

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
| **AdS, boundary time band** (CPR protocol; WdW perturbative holography) | **ACCESS CONSTRAINED — CONDITIONAL ON ACCESS DATA.** At the verified low-energy / perturbative AdS scope, gravity excludes a bulk-state degree of freedom that remains independent of the specified near-boundary access data [GR3-04]. **Gravity does not determine ε_t, Λ or the operator class** |
| **AdS, robustness** (G3-2 illustration) | **KINEMATIC ILLUSTRATION OF ALGEBRAIC ACCESS VS ROBUSTNESS** [GR3-03]. *In the truncated free-field toy only:* every ε_t > 0 is algebraically invertible, the conditioning grows like ε_t^{−(4N−2)} for the tested construction, and fixed conditioning needs a wider band as N grows. Not a gravity theorem |
| **Asymptotically flat, finite retarded time** (BCHW) | **A_time ⊗ A_interface(Bondi-charge readout): CONSTRAINED-NONUNIQUE — FINITE-TIME CHARGE ACCESS FORBIDDEN IN CLASS** [GR3-01]. Finite windows remain admissible, and Bondi-mass readout remains a meaningful target. Only the combination is forbidden. **A_time itself: SUPPLIED / NOT SELECTED** |
| **Coarse-graining** (Raju) | **COARSE-GRAINING SUPPLIED.** No gravity-derived principle uniquely selects the coarse-grained algebra; the result is generally state-dependent |
| **Perturbative order** | **PERTURBATIVE ORDER SUPPLIED.** No canonical nested hierarchy is supplied by the sources used |

- **A_time ↔ A_interface coupling: ASYMPTOTICS-DEPENDENT.**
  - In AdS, boundary energy data is usable in any time band.
  - In flat space, even the gravitational charge labels (Bondi mass) are not finite-time accessible.
- **CONDITIONALLY SELECTED: none. TRUE COMPRESSION: 0. Empirical payoff: none.**

## G3-1 AdS finite-time boundary protocol (Chowdhury–Papadoulaki–Raju, arXiv:2008.01740)

**PRIMARY ARXIV ABSTRACT VERIFIED** [GR3-02] for exactly the following (details beyond the abstract keep their prior
owner-audited / source-text grade):
- observers live near the boundary of global AdS and do not leave the near-boundary region;
- they use simple low-energy unitaries and make measurements in a small interval of time;
- gravitational backreaction and vacuum entanglement are essential;
- the low-energy bulk state can be completely identified;
- the protocol fails in theories without gravity, including nongravitational gauge theories.

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

**What gravity does remove** [GR3-04]. At the verified low-energy / perturbative AdS scope, gravity excludes a bulk-state
degree of freedom that remains independent of the specified near-boundary access data. This is **not** generalized to
arbitrary-energy states, nonperturbative AdS quantum gravity, arbitrary observer algebras, or all asymptotics. Without gravity the same protocol fails, so the combination

> (bulk-state degree of freedom independent of the access data) + (specified near-boundary time-band access, low-energy sector)

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

**Establishes only** [GR3-03]: **KINEMATIC ILLUSTRATION OF ALGEBRAIC ACCESS VS ROBUSTNESS** (exact algebraic access ≠ robust
operational access *in this truncated free-field model*). The sources give "small time band" (CPR) and "infinitesimal
interval" (WdW) at their scopes. The statements "every ε_t > 0 works" and the exponent 4N − 2 belong to this toy, not to
any gravity theorem.
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

**Terminal** [GR3-01]: **A_time ⊗ A_interface(Bondi-charge readout): CONSTRAINED-NONUNIQUE — FINITE-TIME CHARGE ACCESS
FORBIDDEN IN CLASS** (asymptotically flat).
- Finite time windows remain admissible.
- A_time itself stays **SUPPLIED / NOT SELECTED**.
- This is a residual **coupling** constraint, not an elimination of A_time.

(Superseded: "A_time CONSTRAINED-NONUNIQUE" as a standalone statement.)

**Comparison with AdS.** In AdS, boundary energy data enters a protocol that works in an arbitrarily small time band. At
null infinity, the analogous total charge is not available in any finite window. Notably, the **charge labels of G2's
perturbative gravitational splitting are themselves not finite-time accessible at null infinity.**

**Recorded: A_time → A_interface coupling is ASYMPTOTICS-DEPENDENT.** The two geometries are not forced into one answer.

## G3-4 Wheeler–DeWitt perturbative control (Chowdhury–Godet–Papadoulaki–Raju, arXiv:2107.14802)

**PRIMARY ARXIV ABSTRACT VERIFIED** [GR3-02].
- At leading nontrivial order in Newton's constant about AdS, the WdW and diffeomorphism constraints force correlations between a component
  of the asymptotic metric and energetic excitations of matter / gravitons.
- Strictly localized excitations are disallowed.
- Two states or density matrices coinciding at the boundary for an infinitesimal interval of time coincide everywhere in
  the bulk.
- This gives perturbative holography at the stated scope. The infinitesimal-time statement is now abstract-verified, no
  longer only owner-stated.

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

**Grade** [GR3-05]:
- Raju's primary abstract verifies only that coarse-graining possibilities are discussed.
- The statements about state-dependent coarse-grained entropy keep their prior **owner source-text verified** grade. They
  are not independently re-verified here.
- **A_resolution: SUPPLIED / NOT SELECTED.**

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
  - **A_time ⊗ A_interface(Bondi-charge readout) CONSTRAINED-NONUNIQUE** (asymptotically flat): finite-time charge access
    is forbidden in class [GR3-01];
  - **ACCESS CONSTRAINED** (AdS): a bulk degree of freedom independent of the specified access data is excluded at the
    verified scope, conditional on access data [GR3-04].

## G3-9 Flat / nongravitational controls

| control | result |
|---|---|
| AdS protocol, gravity off | fails in theories without gravity, including nongravitational gauge theories (primary abstract verified [GR3-02]). The independent-bulk combination is admissible at G = 0 |
| null infinity, gravity off | no Bondi mass; the finite-window statement has no gravity-off analogue here |
| WdW, gravity off | no constraint, so no forced boundary–bulk correlation |
| conditioning toy | identical with gravity off (kinematic) |
| ordinary QFT | AQFT locality / split independence is the comparison; independent interior specification is available under the split premises |

## G3-10 Orientation firewall

**Owner scope ruling (with GRAVITY REPAIR 03): no orientation gate is run.** Terminal: **ORIENTATION REMAINS SUPPLIED /
NOT SELECTED** (a campaign scope ruling, not a theorem).

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

| component | refinement (final G3 mapping, GR3-06) |
|---|---|
| Σ / A_partition | **CONSTRAINED / OBSERVABLE-CLASS-PRICED** |
| A_interface | **RELOCATION / LOAD-BEARING** (dressing, charge / boundary readout, asymptotic-observable choice) |
| A_resolution | **SUPPLIED / NOT SELECTED**. Coupled kinematically to time-band robustness in the toy; gravity does not select it |
| A_time | **SUPPLIED / NOT SELECTED**. Specific combinations with observable targets are constrained: asymptotically flat, finite time × Bondi-charge readout is forbidden; AdS, small / infinitesimal boundary-time data can determine low-energy / perturbative bulk information under the stated assumptions |
| perturbative order | **SUPPLIED** |
| asymptotic structure | **SUPPLIED and load-bearing** |
| charge sector | NEW SUPPLIED LABEL; not finite-time accessible at null infinity |
| orientation | **SUPPLIED / NOT SELECTED** (campaign scope ruling) |

**General statement** [GR3-06]: **gravity constrains admissible combinations of subsystem, interface, time and
observable-class data. It does not select the individual residual entries.**

**TRUE COMPRESSION = 0.** All results are auxiliary to canonical GRUT.
