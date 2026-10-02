# G3 ACCESS LEDGER — access data per setting, what is supplied, and what gravity contributes

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
| **AdS boundary-observer protocol** (Chowdhury–Papadoulaki–Raju, arXiv:2008.01740) | global AdS geometry (∂); a unique vacuum; low-energy cutoff Λ ≪ M_Pl; time band [0, ε_t]; simple low-energy unitaries and low-energy / projective measurements (𝒪); semiclassical gravity | gravitational backreaction / boundary energy data plus vacuum entanglement let near-boundary observers identify the bulk state within the stated low-energy Hilbert space | the protocol's bulk-identification mechanism **fails** (owner-stated source result) | **ACCESS CONSTRAINED** (a hidden-interior combination is excluded) **+ CONDITIONAL ON ACCESS DATA** (ε_t, Λ, 𝒪 declared, not derived) |
| **Time-band conditioning toy** (G3-2; `g3/g3_timeband.py`) | the AdS₄ l = 0 spectrum ω_n = 3 + 2n; mode cutoff N (a Λ proxy); time band ε_t | **none** (kinematic Fourier analysis) | **identical** | **EXACT ALGEBRAIC ACCESS ≠ ROBUST OPERATIONAL ACCESS** (illustration) |
| **Asymptotically flat finite time** (Bousso–Chandrasekaran–Halpern–Wall, arXiv:1709.08632) | asymptotically flat geometry; null infinity; a finite retarded-time window; the asymptotic entropy bounds the result rests on; the definition of the observable algebra on a finite portion of null infinity | the Bondi mass (a gravitational charge) is **not** in the algebra on any finite portion of null infinity; attempts at large radius in fixed retarded time are thwarted by quantum fluctuations | the Bondi mass itself is absent; no finite-time statement transfers | **A_time CONSTRAINED-NONUNIQUE** (asymptotically flat); conditional on the entropy-bound premises |
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
- The condition number scales roughly as ε_t^{−(4N−2)} at small ε_t.

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
| Σ / A_partition | CONSTRAINED / OBSERVABLE-CLASS-PRICED (G2). **AdS:** an interior subsystem hidden from every boundary time band is excluded *given* the access data (ACCESS CONSTRAINED, conditional) |
| A_interface | **RELOCATION:** the boundary detector / operator class, dressing, and the null-infinity window are supplied |
| A_resolution | **SUPPLIED** (Λ, 𝒪, coarse-graining). Kinematically coupled to A_time (conditioning), not selected by gravity |
| A_time | **AdS:** CONDITIONAL ACCESS (any ε_t > 0 suffices algebraically; robustness degrades). **Asymptotically flat:** CONSTRAINED-NONUNIQUE (finite windows exclude the Bondi mass). **ASYMPTOTICS-DEPENDENT** |
| orientation | supplied (retarded time, future null infinity, the positive AdS Hamiltonian); not attacked (G3-10) |

## 5. Source grades

| source | grade | note |
|---|---|---|
| Bousso–Chandrasekaran–Halpern–Wall, PRD 97, 046014 (2018), arXiv:1709.08632 | **PRIMARY ABSTRACT VERIFIED** (abstract text retrieved via search) | body / exact assumptions not re-read |
| Chowdhury–Godet–Papadoulaki–Raju, JHEP 03 (2022) 019, arXiv:2107.14802 | **ABSTRACT LOCATED** (summary retrieved) | the "infinitesimal time interval" statement is **owner-stated**, not re-read here |
| Chowdhury–Papadoulaki–Raju, SciPost Phys. 10, 106 (2021), arXiv:2008.01740 | **SOURCE LOCATED** | all protocol facts in §2 are **owner-stated source-text facts**, not re-read here. Caution: one search summary attributed to this paper a coarse-grained time-band algebra with a non-trivial commutant in states containing a macroscopic bulk observer. This may conflate it with later work (e.g. "Holographic observers for time-band algebras", JHEP 06 (2025) 242). **Not used as evidence** |
| Raju, arXiv:2110.05470 | PRIMARY / SOURCE-TEXT VERIFIED (owner) | coarse-graining / state-dependence statement per the owner's audit |
