# S2-Σ (+ S2-G1) RESULT — factorization and local-dimension selection without hidden weights

**Charter:** `probes/PROBE_CHARTERS.md`, Wave 2 §S2-Σ. Pre-registered at `a4d433e`, before any run.

**Files:**
- `probes/S2-Sigma/s2_sigma.py` (+ `.log`): Σ-1 … Σ-7;
- `probes/S2-Sigma/s2_sigma_scope.py` (+ `.log`): Σ-8 scope check.

**Wall attacked:** Σ, and G1 (local dimension).

**Setup.**
- n = 6 qubit slots, N = 64.
- **1414 candidates**, all set partitions of the 6 slots into ≥ 2 blocks. That gives 202 groupings with local dimensions 2 … 32, including mixed ones, and all 10 local-dimension types. Each grouping is crossed with 7 frames:
  - the identity;
  - 4 random Clifford circuits;
  - a Haar unitary;
  - the commutant frame W = e^{−iH·0.7}.
- **Criteria vector:** V = (L, Q, R, P, M). Each criterion is LU-invariant and oriented so that larger is better:
  - L = locality;
  - Q = Carroll–Singh-type entanglement-growth penalty;
  - R = record redundancy;
  - P = autonomy;
  - M = description length.
- **Models:**
  - M1: clustered pairs plus weak bridges;
  - M2: TFIM ring (translation-symmetric);
  - M3: record-forming star;
  - M4: generic random 2-local chain (the positive control).
- **States:**
  - product;
  - random local product;
  - ground;
  - mid-spectrum eigenstate;
  - Haar;
  - record-forming;
  - Gibbs β = 1.

**Labels:**
- **PARETO NONUNIQUENESS** (outcome B);
- **DIMENSION-PRICED** (C);
- **SCALE-PRICED** (D);
- **DEGENERACY-PROTECTED NONUNIQUENESS** (Σ-7);
- **STATE-PRICED (→ H)**;
- **LOCAL DIMENSION REMAINS A Σ PRIMITIVE**;
- outcome A (**ROBUST TPS SELECTOR**): **NOT OBTAINED** in the open candidate set. It is obtained only in restricted scopes (Σ-8).

## 0. Headline

> **No candidate weakly dominates in any of 58 runs**: 28 (model × state) runs plus 30 scale runs.
> - Every Pareto front holds between 3 and 148 inequivalent TPSs.
> - In the unrestricted runs, the front holds **5 to all 10 local-dimension types**.
> - Choosing among them needs weights, normalizations, a scale and a state.
>
> The only dominance found anywhere appears after **supplying the local dimension and number of factors** (CPR scope).
> Even then it is unique only **up to the commutant of H**. Adding ψ breaks that degeneracy only by choosing a
> **time origin**.

## Σ-1 Pareto fronts (τ = 1, δ = 0.1, ε = 0)

| model | front sizes over 7 states | local-dim types on front | weakly dominant |
|---|---|---|---|
| M1 clustered | 23 – 49 | 7 – 10 | none (0/7) |
| M2 TFIM ring | 18 – 148 | 7 – 9 | none (0/7) |
| M3 record star | 19 – 45 | 7 – 10 | none (0/7) |
| M4 generic chain | 9 – 40 | 5 – 10 | none (0/7) |

**Systematic trade-off:**
- L and P are maximized by **coarse** groupings, e.g. the 16×4 grouping (01)(2345), where they reach k_eff = 1.0003 and
  autonomy 0.9997 for M1.
- M is maximized only by the **finest** grouping (qubits).
- Q and R are maximized by state-specific frames.

The criteria disagree, as functions of local dimension, before any weighting enters.

## Σ-2 Monotone-transformation hostile

- **Exact statement:** Pareto dominance is invariant under every strictly increasing componentwise transform. The
  fronts above are therefore objective-scale-invariant by construction. This is the only invariant object.
- **Scalar winners are not invariant.** Across 4 transforms (identity, log1p, square, sqrt) × 3 normalizations
  (min-max, z-score, rank), equal weights give **1 – 5 distinct winners per run**.
- Only 4 of the 28 main runs have a single transform-stable winner. In 3 of those 4 the winner sits in the commutant
  frame, i.e. it is epoch-dependent (Σ-7b).

## Σ-3 Weight simplex (Dirichlet(1) on the 4-simplex, 20 000 samples; the sampling measure is supplied)

- **Min-max normalization:** 4 – 26 distinct winners per run.
- **Rank normalization:** 12 – 58 distinct winners per run.
- A top share **> 0.9 occurs once in 28 main runs**: M4 with the random local product state, 0.91, a commutant-frame
  winner.
- Large weight regions are occupied by **different local dimensions**. Example, M1 Haar state:
  - 0.83 for 4×4×4 (01)(23)(45);
  - 0.07 for 4×4×2×2;
  - 0.07 for 8×8.
- Calling any weight vector "natural" is forbidden. None is singled out.

## Σ-4 Local dimension as a variable (G1)

- The fronts contain qubits, ququarts, octits and mixed (8×4×2, 16×2×2, …) TPSs at the same time.
- The local-dimension choice is a **different** question from the TPS choice: for fixed d, frames still compete (Σ-8).
- The same objective family cannot jointly select d: L/P against M is a direct dimension conflict.
- → **LOCAL DIMENSION REMAINS A Σ PRIMITIVE.**

## Σ-5 Scale hostile (M1, M4; record-forming state)

Over the grid τ ∈ {0.2, 1, 5}, δ ∈ {0.1, 0.3}, fragment ∈ {single, pairs} and GUE perturbation ε ∈ {0, 0.05, 0.3}
(the perturbation measure is supplied):
- dominance **never** appears;
- the front size moves between 3 and 58;
- the top simplex winner changes across τ, ε, δ and the fragment definition. Example for M4 at ε = 0:
  - 8×4×2 at τ = 0.2 and 1;
  - 16×2×2 at τ = 5;
  - 4×2×2×2×2 with pair fragments.

No stable scale window exists. → **SCALE-PRICED.**

## Σ-6 State dependence

For fixed H, the top scalar / simplex winner changes with ψ in every model. Example, M4:

| state | top winner |
|---|---|
| product | commutant 8×4×2 |
| ground | id 16×2×2 / 8×2×2×2 (tie 0.24) |
| mid-spectrum | id 4×4×4 |
| Gibbs | id 8×4×2 |

In eigenstates and Gibbs states Q is identically 0 for every candidate (no dynamics). State-dependent criteria then
lose all discriminating power. → **STATE-PRICED**: the cost moves into **H** (the universe's state picks its
subsystems).

## Σ-7 Symmetry / degeneracy

**(a) Translation ring (M2).** Translation-related pairings (01)(23)(45) and (12)(34)(05) have *identical* vectors.
110 of the 148 front members share an exactly identical vector with another candidate.

**(b) Commutant frames, in every model.** For W = e^{−iHs}:
- the H-only criteria (L, P, M) are **identical** to the identity frame, with max difference ≤ 4.4·10⁻¹⁶;
- the state criteria equal the identity-frame values **at the time-shifted state ψ(−s)**, with max difference 0.0;
- they differ from the identity frame at ψ by up to 2.87.

So every H-only selector is blind to the whole commutant {e^{−if(H)}}. For a non-degenerate spectrum that is an
N-torus of inequivalent TPSs. Adding ψ breaks the degeneracy only by picking an **epoch**. That is why commutant frames
often "win": they are the identity frame at a better-chosen time origin.

→ **DEGENERACY-PROTECTED NONUNIQUENESS**, both discrete (translation) and continuous (commutant).

## Σ-8 Literature conflict, settled by scope in these toy models

Scope check (`s2_sigma_scope.py`): candidate sets are restricted step by step.

| scope | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| d = 2, n = 6 fixed, H-only (L, P, M) | none | **id ≡ all commutant frames** | none | **id ≡ all commutant frames** |
| d = 2, n = 6, H + state | none | none | none | single: **commutant s = 0.7** (an epoch) |
| d = 4 pairings fixed, H-only | **id (01)(23)(45) ≡ commutant frames** | none | none | none |

| Work | Scope classification |
|---|---|
| **Cotler–Penington–Ranard** (arXiv:1702.06142, owner-cited) | **CONFIRMED IN SCOPE.** With the local dimension, the number of factors and the locality framework supplied, and when H is genuinely local in that d (M2, M4 for d = 2; M1 for d = 4), locality selects the TPS **up to the commutant of H**. Outside its native d, or with the wrong d supplied, no dominance. The theorem's inputs (k, d) are exactly the Σ primitives found here |
| **Carroll–Singh** (PRA 103, 022213, owner-cited) | **SCOPE-PRICED.** Quasiclassicality (Q) needs a state and a time horizon. It breaks the commutant degeneracy only by choosing an epoch, and its winner changes with ψ, τ and d. A preferred factorization "from H" via quasiclassicality is in fact from (H, ψ, τ, d, criterion), and it is not robust to adding the other natural criteria (L, R, P, M) without weights |
| **Stoica** (arXiv:2103.15104, owner-cited) | **CONFIRMED IN SCOPE.** H-only structure is invariant under the commutant (Σ-7b, exact), so it cannot be unique. Adding ψ makes it unique only relative to a time origin plus supplied criteria |

All three references are **owner-cited**. Fetching arXiv and APS was blocked by the egress proxy, so the abstracts were
**not verified** in this session.

The three claims are **compatible**. They differ in which Σ items they assume.

## Σ-9 Information accounting

| entry | status |
|---|---|
| H | given |
| ψ | given; **changes the winner** (Σ-6) |
| total dimension N | given |
| **local dimension d** | **supplied**: not jointly selectable; L/P vs M conflict |
| **number of factors** | **supplied**: ≥ 2 constraint; trivial factorization otherwise |
| **objective family** | **supplied**: 5 criteria disagree |
| **weights** | **supplied**: 4 – 58 simplex winners |
| **normalization / units** | **supplied**: 1 – 5 scalar winners across transforms |
| **time scale τ** | **supplied**: winner changes |
| **coarse scale δ** | **supplied**: winner changes |
| **fragment definition** | **supplied**: winner changes |
| **perturbation measure / ε** | **supplied**: winner changes |
| **epoch / time origin** | **supplied**: commutant degeneracy is broken only by it |
| **candidate family** (measure over TPSs) | **supplied**: a finite family; dominance could only get *harder* with more candidates |
| commutant of H, symmetries | **gauge for H-only selectors** (exact degeneracy) |

**Σ DERIVED requires every entry except H and ψ to be eliminated or gauge. Here none is eliminated.** → Σ is **not
derived**.

**Status: S2-Σ + G1 COMPLETE.**
- PARETO NONUNIQUENESS in 58/58 runs.
- DIMENSION-PRICED; LOCAL DIMENSION REMAINS A Σ PRIMITIVE.
- SCALE-PRICED; STATE-PRICED.
- DEGENERACY-PROTECTED (translation; the commutant of H, exact).
- A robust selector appears only in CPR scope (d, number of factors and locality supplied; H native to that d). It is
  unique up to the commutant, and ψ breaks that only by choosing an epoch.
