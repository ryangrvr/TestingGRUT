# G3 ACCESS LEDGER — access data per setting, what is supplied, and what gravity contributes

> **Repaired by GRAVITY REPAIR 03 (GR3-01 … GR3-06; `GRAVITY_CORRECTION_LEDGER.md`).** G3 is accepted provisionally
> after that repair. Where wording differs, the ledger takes precedence. Script and log outputs are kept as emitted.

## 1. Operational resolution is not one ε (G3-0)

Resolution is recorded as a **vector of access data**. Each entry is supplied unless a source derives it:

| access datum | symbol |
|---|---|
| accessible observable algebra / operator class (incl. "simple" operator degree) | 𝒪 |
| energy cutoff (low-energy Hilbert space) | Λ |
| time band / retarded-time window | ε_t |
| perturbative order | k |
| asymptotic region (AdS boundary; null infinity) | ∂ |
| measurement precision (only where physically defined) | δ |

The G2 toy's ε is **none of these**: it is a numerical threshold for discarding small generators (G3-7).

## 2. Setting-by-setting ledger

| setting / source | access data **supplied** | what **gravity** contributes | gravity off (G3-9) | terminal |
|---|---|---|---|---|
| **AdS boundary-observer protocol** (Chowdhury–Papadoulaki–Raju, arXiv:2008.01740) | global AdS geometry (∂); a unique vacuum; low-energy cutoff Λ ≪ M_Pl; time band [0, ε_t]; simple low-energy unitaries and low-energy / projective measurements (𝒪); semiclassical gravity | gravitational backreaction / boundary energy data plus vacuum entanglement let near-boundary observers identify the bulk state within the stated low-energy Hilbert space | the protocol **fails** in theories without gravity, including nongravitational gauge theories (primary abstract) | **ACCESS CONSTRAINED — CONDITIONAL ON ACCESS DATA** (a bulk degree of freedom independent of the specified access data is excluded at the verified low-energy scope [GR3-04]; ε_t, Λ, 𝒪 declared, not derived) |
| **Time-band conditioning toy** (G3-2; `g3/g3_timeband.py`) | the AdS₄ l = 0 spectrum ω_n = 3 + 2n; mode cutoff N (a Λ proxy); time band ε_t | **none** (kinematic Fourier analysis) | **identical** | **KINEMATIC ILLUSTRATION OF ALGEBRAIC ACCESS VS ROBUSTNESS** (toy-only; not a gravity theorem) [GR3-03] |
| **Asymptotically flat finite time** (Bousso–Chandrasekaran–Halpern–Wall, arXiv:1709.08632) | asymptotically flat geometry; null infinity; a finite retarded-time window; the asymptotic entropy bounds the result rests on; the definition of the observable algebra on a finite portion of null infinity | the Bondi mass (a gravitational charge) is **not** in the algebra on any finite portion of null infinity; attempts at large radius in fixed retarded time are thwarted by quantum fluctuations | the Bondi mass itself is absent; no finite-time statement transfers | **A_time ⊗ A_interface(Bondi-charge readout): CONSTRAINED-NONUNIQUE — FINITE-TIME CHARGE ACCESS FORBIDDEN IN CLASS** [GR3-01]; A_time itself supplied / not selected; conditional on the entropy-bound premises |
| **Wheeler–DeWitt perturbative holography** (Chowdhury–Godet–Papadoulaki–Raju, arXiv:2107.14802) | an AdS background (∂); leading nontrivial perturbative order (k); the WdW / diffeomorphism constraints; the state class | the constraints force correlations between a component of the asymptotic metric and bulk energetic excitations. Per the owner's audit, states agreeing on the boundary for an infinitesimal time interval agree in the bulk at the stated scope | no constraint ⇒ no forced correlation | **CONSTRAINED** (independent perturbative support for a gravity-constrained observable algebra) **+ RELOCATION** into the supplied AdS boundary and order |
| **Coarse-graining** (Raju, arXiv:2110.05470) | which coarse-grained observable set is used | the fine-grained boundary completeness; the coarse-grained entropy / approximate split is **generally state-dependent** | ordinary split | **COARSE-GRAINING SUPPLIED** |
| **Perturbative order** (Donnelly–Giddings O(κ); WdW leading order) | k, separately in each source | different observable classes at "leading order" | — | **PERTURBATIVE ORDER SUPPLIED** (no canonical nested hierarchy O₁ ⊂ O₂ ⊂ … → fine-grained algebra in the sources used) |

## 3. G3-2 numbers (illustration; `g3/g3_timeband.log`, 160-digit arithmetic)

**Condition number** of the time-band moment system (rows N, columns ε_t):

| N | ε_t = π | π/2 | π/4 | π/8 | π/16 | π/32 |
|---|---|---|---|---|---|---|
| 2 | 1 | 5.7 | 2.4e2 | 2.6e4 | 2.0e6 | 1.3e8 |
| 4 | 1 | 7.0e2 | 3.2e7 | 8.5e11 | 2.0e16 | 3.6e20 |
| 8 | 1 | 1.6e8 | 1.1e18 | 2.6e27 | 5.2e36 | 7.2e45 |

- The minimum-norm smearing for the lowest mode grows from 1 (at ε_t = π) to 2.0e4 (N = 2) and 2.8e22 (N = 8) at π/32.
- The condition number scales roughly as ε_t^{−(4N−2)} at small ε_t, **for the tested construction only** [GR3-03].

**Fixed-robustness trade-off** (smallest ε_t/π with cond ≤ 10⁶):

| N | 2 | 3 | 4 | 6 | 8 |
|---|---|---|---|---|---|
| ε_t/π | 0.070 | 0.196 | 0.314 | 0.476 | 0.580 |

The admissible robust time band **grows with the mode cutoff**. This is a **kinematic** A_time ↔ Λ coupling, present with
gravity off.

**Firewall.** The condition number and coefficient norm are **not** identified with a physical measurement precision.

## 4. Residual entries after G3

| component | status |
|---|---|
| Σ / A_partition | **CONSTRAINED / OBSERVABLE-CLASS-PRICED** [GR3-06]. AdS: a bulk degree of freedom independent of the specified access data is excluded at the verified scope [GR3-04] |
| A_interface | **RELOCATION / LOAD-BEARING** |
| A_resolution | **SUPPLIED / NOT SELECTED** (kinematic coupling to time-band robustness in the toy only) |
| A_time | **SUPPLIED / NOT SELECTED**. Combination constraints: flat, finite time × Bondi-charge readout forbidden; AdS, small / infinitesimal boundary-time data determines low-energy / perturbative bulk information under the stated assumptions |
| perturbative order | SUPPLIED |
| asymptotic structure | SUPPLIED and load-bearing |
| orientation | **SUPPLIED / NOT SELECTED** (campaign scope ruling; no orientation gate) |

## 5. Source grades

| source | grade | note |
|---|---|---|
| Bousso–Chandrasekaran–Halpern–Wall, PRD 97, 046014 (2018), arXiv:1709.08632 | **PRIMARY ARXIV ABSTRACT VERIFIED** [GR3-02]: the Bondi mass is not measurable in finite retarded time; it is not in the algebra of any finite portion of future null infinity; large-radius / fixed-retarded-time attempts are obstructed by quantum fluctuations; the derivation is tied to asymptotic entropy bounds | detailed body assumptions unverified |
| Chowdhury–Godet–Papadoulaki–Raju, JHEP 03 (2022) 019, arXiv:2107.14802 | **PRIMARY ARXIV ABSTRACT VERIFIED** [GR3-02]: leading nontrivial order in G about AdS; forced asymptotic-metric / energetic-excitation correlations; strictly localized excitations disallowed; boundary coincidence for an infinitesimal interval ⇒ bulk coincidence; perturbative holography | — |
| Chowdhury–Papadoulaki–Raju, SciPost Phys. 10, 106 (2021), arXiv:2008.01740 | **PRIMARY ARXIV ABSTRACT VERIFIED** [GR3-02]: near-boundary observers in global AdS; simple low-energy unitaries; measurements in a small time interval; observers stay near the boundary; backreaction and vacuum entanglement essential; low-energy bulk state completely identified; fails without gravity, incl. nongravitational gauge theories | details beyond the abstract keep their prior owner-audited grade. Caution: one search summary attributed to this paper a coarse-grained time-band algebra with a non-trivial commutant in states containing a macroscopic bulk observer. This may conflate it with later work (e.g. "Holographic observers for time-band algebras", JHEP 06 (2025) 242). **Not used as evidence** |
| Raju, arXiv:2110.05470 | PRIMARY / SOURCE-TEXT VERIFIED (owner) | [GR3-05] the primary abstract verifies only that coarse-graining possibilities are discussed. The state-dependence statements keep the prior owner source-text grade and are not re-verified here |
