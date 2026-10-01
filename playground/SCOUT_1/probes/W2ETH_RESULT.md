# SCOUT-1 W2-ETH RESULT — KMS without a supplied passive preparation? (closed-system thermalization)

**Charter:** `PROBE_CHARTERS.md` §W2-ETH. Preregistered outcomes: ATTRACTOR SELECTION (prices) / NO.

**Files:** `w2eth_thermalization.py`, with log `w2eth_thermalization.log`.

**Setup:**
- Model: mixed-field Ising ring `H = Σ ZZ + gΣX + hΣZ`, L = 8, 10, 12, with exact diagonalization.
- Initial states: product states rotated by θ.
- Late-time state: the infinite-time average (diagonal ensemble).
- Comparison: against the Gibbs state at the β that reproduces `⟨H⟩`, using the trace distance of the
  2-site reduced density matrices.
- The states sit above the infinite-temperature energy of this sign convention, so β < 0 throughout. For
  `−H` the identical numbers hold with β > 0.

**Labels:**
- ATTRACTOR SELECTION (premise-priced: genericity);
- KNOWN RESULT IMPORT: ETH (Deutsch 1991; Srednicki 1994; Rigol–Dunjko–Olshanii 2008); GGE for integrable chains —
  STANDARD-TEXTBOOK, ED-indicated ✓ (finite-size);
- NON-DISTINCTIVE.

## 0. Verdict

> **ATTRACTOR SELECTION (shape) — W1-P's preparation price is traded for a genericity price.**
>
> **What emerges.** In the non-integrable chain, unitary dynamics alone drives the subsystem toward a
> single-β Gibbs state from non-thermal product states. No passive preparation is supplied. The distance to
> Gibbs falls with L for every θ:
>
> | θ | L = 8 | L = 10 | L = 12 | initial distance |
> |---|---|---|---|---|
> | π/2 | 0.114 | 0.086 | 0.079 | 0.65 |
> | 3π/4 | 0.157 | 0.108 | 0.094 | 0.75 |
> | 0 | 0.0209 | 0.0206 | 0.0198 | |
> | π/4 | 0.0436 | 0.0412 | 0.0398 | |
>
> **The integrable control is not cleanly thermal.** Two of four θ drift *away* from Gibbs as L grows
> (0.034 → 0.039 → 0.042), and two drift slowly down. This is consistent with relaxation to a GGE rather than
> Gibbs. It is finite-size-limited, and the import is standard.
>
> **Prices:**
> - **Non-integrability of the supplied Hamiltonian** (generic, measure-one, but supplied).
> - **Subsystem restriction** — the global state stays pure.
> - **ETH itself**, which is a hypothesis supported numerically, not a theorem.
> - **The temperature is NOT selected.** β is fixed by the supplied initial energy density (−0.72, −0.63, −0.30,
>   −0.20 as θ varies). This is W1-C again: only the shape (single β) is attracted to.

## 1. What this means

- **First dynamical replacement of a Wave-1 price.** The *preparation* premise of W1-P (complete passivity)
  is replaced by **dynamics + genericity**. That is a different and weaker kind of premise: it holds for
  almost every Hamiltonian, rather than being a special supplied value.
- **A sixth premise class appears: genericity.** "Holds for almost all members" (non-integrability) is
  neither a symmetry nor a preparation nor a tuning. It is the closest thing to an *earned* condition found
  so far, because it **removes** information instead of adding it.
- **The GRUT record's own baths are mostly quadratic (integrable).** C1-a/L0-1 chains, P-17 oscillator baths
  and the S6 parent are free systems, and these are exactly the class that relaxes to GGE, not KMS. So within
  the record's declared parents, **genericity is absent by construction**, and the KMS shape is not attracted
  to. **Cross-layer NO-GO (scoped):** the quadratic/free scope that makes the record's reductions exact (E-3
  linearity ⇒ exact reduction) is the scope in which dynamical thermalization fails.
- **Hostile notes:**
  - L ≤ 12; the decay rates are not extrapolated.
  - The diagonal ensemble ignores degenerate-block coherences. Product states are translation- and
    parity-symmetric, so this is a small effect, but it is not quantified.

**Status: W2-ETH COMPLETE — ATTRACTOR SELECTION of the KMS *shape* under a genericity price (non-integrability).
T stays supplied (initial energy). Record's quadratic parents excluded from the mechanism (E-3 tension).**

## Addendum — does genericity fix a universal *number*? (`w2eth_level_stats.py`, `.log`)

**Setup:** level-spacing ratio `⟨r⟩`, L = 11 open chain, random site fields (no spatial symmetry), three
realizations per case.

| Case | `⟨r⟩` per realization | Mean | Reference |
|---|---|---|---|
| generic real | 0.509, 0.487, 0.521 | 0.506 | GOE 0.5307 |
| integrable random TFIM | 0.403, 0.373, 0.354 | 0.377 | Poisson 0.3863 |
| "complex" fields (added Y) | 0.527, 0.523, 0.524 | 0.525 | GOE-like |

- **Genericity does fix a universal weight-0 number, but only within a supplied symmetry class.** The value
  of `⟨r⟩` is set by the antiunitary class (Dyson's threefold way).
- **Hostile catch (recorded post-hoc in the log).** The "complex" control was meant to be GUE but is **not**.
  Single-site X/Y fields can be rotated into the x–z plane by local z-rotations, which commute with ZZ. The
  model therefore has a hidden antiunitary and is GOE. This is a live example of hostile test #2: *the
  symmetry class is a supplied datum, and it can be hidden.*
- **Pattern:** genericity × (supplied symmetry class) × (imported theorem: RMT universality). TP-1 again,
  with the genericity class as the principle.
