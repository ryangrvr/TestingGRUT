# S2-D-arrow RESULT — can a closed unitary universe select an arrow?

**Charter:** `probes/PROBE_CHARTERS.md` §S2-D-ARROW. Pre-registered at `04b0152`, theorem-first.

**Files:**
- `results/S2-D-arrow_D0_THEOREM.md` (D0: theorem + Proposition 2);
- `probes/S2-D-arrow/s2_darrow.py` (+ `.log`): D1 – D10;
- `probes/S2-D-arrow/d2_check.py` (+ `.log`): large-sample D2.

**Main model:**
- the real (time-reversal-symmetric) chaotic mixed-field Ising chain H = Σ ZZ + 0.9045 Σ X + 0.809 Σ Z;
- n = 10 for pure states and n = 8 for mixed states;
- S = the first qubit(s) of the computational TPS unless stated.

**Labels:**
- outcome **A (CLOSED-UNITARY ARROW SELECTED): EXCLUDED** (theorem, in scope);
- **B (ARROW NEEDS SPECIAL GLOBAL STATE)**: holds, sharpened to **H_CORR + FRESHNESS**, not low entropy;
- **C (ARROW NEEDS Σ + H)**: holds;
- **D (TYPICALITY ONLY)**: typicality gives **equilibrium, not an arrow**;
- **E (TWO-SIDED / JANUS)**: holds. The special middle state can be fixed by **(H, Σ) + a selector functional**, never
  by H alone. **TIME ORIENTATION NOT SELECTED.**

## STRONG verdict (every admissible state)

> **EVERY-STATE ARROW = IMPOSSIBLE** (D0, Corollary 1). This covers every continuous functional: subsystem entropy,
> coarse-grained entropy, mutual information, record diagnostics. In a finite closed unitary system no such functional
> is monotone along any orbit forever (Theorem 1, recurrence). For every non-conserved functional there is an
> admissible state whose forward evolution decreases it.

Numerical confirmations:
- **D1 (Corollary 2, exact).** The time-reversed state Kψ_prod(t* = 8) sends S_S from 1.238 to 0.000 going
  forward. The mirror identity holds to 8.9·10⁻¹⁶.
- **Recurrence (Lemma 1).** With 3 qubits, fidelity reaches 0.99999 at t = 46 251.7, and S(qubit 0) returns to 10⁻⁴.
- **Eigenstate.** S_S = 1.883, constant.

## WEAK verdict (typical states)

| measure | ⟨S_S(0)⟩ (max 2) | ⟨ΔS⟩ over Δt = 1 | P(increase) |
|---|---|---|---|
| Haar on ℂ¹⁰²⁴ (150 samples) | 1.990 | +0.0002 | 0.52 |
| Haar on the energy shell \|E\| < 1 (1500 samples) | 1.978 | +0.00004 ± 0.00019 | **0.499 ± 0.013** |

(A first 150-sample shell run gave 0.61. That is a fluctuation: invariance of the measure under U_t and K forces 0.5.)

**EQUILIBRIUM TYPICALITY ≠ ARROW.**
- Typical states are already at equilibrium at t = 0.
- There is no entropy budget and no direction.
- What typicality selects is the *equilibrium appearance*. It is MEASURE-, TPS- and constraint-priced (the shell).

## What the effective arrows that do occur actually need

### D3 / D6 — correlation, not bath entropy (n = 8, S = qubit 0)

States C and D have **identical marginals**. The bath starts nearly maximally mixed (6.902 of 7 bits). Entries are
S_S at t = 0, 1, 2, 3, 5:

| state | S_S(t) | I(S:E)(0) |
|---|---|---|
| A: product \|0⟩⟨0\| ⊗ 𝟙/128 | 0 → 0.825 → 0.999 → 0.955 → 0.979 (rises) | 0 |
| C: correlated, Kρ_A(τ = 3)K | 0.955 → 0.999 → 0.825 → **0.000** → 0.999 (**recoheres to purity**) | 0.857 |
| D: product of C's marginals | 0.955 → 0.997 → 0.986 → 0.968 → 0.993 (relaxes, stays near max) | 0 |
| λC + (1−λ)D, λ = 0.3 | dips only to 0.857 | — |
| λC + (1−λ)D, λ = 0.7 | dips to 0.532 | — |

- With the **same marginals and an almost maximally mixed bath**, one state relaxes and the other runs an
  anti-arrow.
- The anti-arrow strength grows with the correlation fraction λ.
- **The relevant boundary datum is the S–E correlation pattern (H_CORR), not low environment entropy.**

### D4 — collision models: bath entropy is not the price; bath correlation breaks Markovianity

Partial-swap collisions with θ = 0.6. The system TD is between preparations |0⟩ and |+⟩.

| fresh ancillas | system TD sequence | monotone? | S_sys(a) | record in ancilla 1 |
|---|---|---|---|---|
| pure \|0⟩ | 0.535 → 0.135 | yes | 0 → 0 (fixed point) | 0.324 |
| thermal p₀ = 0.8 | 0.502 → 0.068 | yes | 0.342 → 0.694 | 0.265 |
| **maximally mixed** | 0.482 → 0.048 | **yes** | **0.633 → 0.997** | **0.225** |
| classically correlated (marginals maximally mixed) | 0.482, 0.263, 0.170, **0.202**, 0.203, 0.150, 0.072 | **no** | 0.633 → 0.997 | 0.225 |
| GHZ ancillas (marginals maximally mixed) | same as classically correlated | **no** | 0.633 → 0.997 | 0.225 |

The five questions, answered separately:

| question | answer |
|---|---|
| (1) attractor | exists in all cases |
| (2) its value | fixed by the **ancilla marginal** (final ρ₀₀ = 1.000 / 0.814 / 0.534) |
| (3) Markovianity | holds for **uncorrelated** fresh ancillas of *any* purity; **fails with correlated ancillas**, even with maximally mixed marginals |
| (4) entropy arrow | present even with a **maximally mixed** bath |
| (5) records | formed even with a maximally mixed bath (0.225) |

→ A **high-entropy bath works**. What it needs is **independence** (no bath correlations) and **freshness**.

### D5 — freshness firewall (maximally mixed bath, 21 collisions)

| bath | system TD | monotone? | S_sys(a) at end | record in ancilla 1 at end |
|---|---|---|---|---|
| fresh | monotone (above) | yes | — | — |
| M = 1 reused | 0.482, 0.037, 0.693, 0.170 … 0.706 | no | **0.633 → 0.007** | 0.001 |
| M = 2 reused | non-monotone | no | 0.436 | 0.032 |
| M = 4 reused | non-monotone | no | 0.938 | 0.081 |
| M = 7 reused | non-monotone | no | 0.985 | 0.165 |

With M = 1 the system **re-purifies**. Removing freshness removes semigroup behaviour.
→ **ARROW FRESHNESS-PRICED.** Freshness is physically a **global state in which not-yet-interacted degrees of freedom
are uncorrelated with the system**. That is an H_CORR statement about the future-interacting part of the universe.

### D7 — record arrow (S + 6 env; the antiunitary symmetry Θ = X_S K is Σ-local)

| | t: 0 → 0.866 | redundancy |
|---|---|---|
| forward from \|+⟩\|0…⟩ | ⟨I(S:E_k)⟩ 0 → 0.221 → 0.485 → 0.733 → 0.921 → 0.985 | 0 → 6 |
| from Θψ(t_f) | ⟨I(S:E_k)⟩ 0.985 → … → **0.000** | **6 → 0** (records un-form) |

A Haar global state gives ⟨I(S:E_k)⟩ = 0.064 with redundancy 0: no records at all.
→ **RECORD ARROW H-PRICED.** It needs a product-like boundary relative to Σ.

### D8 — Σ-dependence (same global state ψ_prod(4), same H)

| frame | S_S(t) for t = 0 … 8 | reading |
|---|---|---|
| computational | 0.930 … 1.238 | mild changes |
| Clifford | 1.94 – 1.98 | flat, at equilibrium |
| Householder (maps the state to a product) | **0.000** → 1.536 → … | an "arrow from a product state" |

Whether a given global state is "special" (low correlation) **depends on the factorization**. In the Householder frame
the state is a product state, but H is non-local in that frame.

→ **"Low-entropy / low-correlation environment" is not TPS-invariant. The arrow price is H_CORR relative to Σ: H + Σ
jointly.** The D0 mirror (Corollary 2) also needs a **Σ-local** time reversal.

### D9 — access

- Global TD between two nearby product starts is **0.0158 at every t**: exact information is preserved.
- S_S rises (0 → 1.013) from the product start and falls (→ 0.000) from its time reverse.
- The coarse magnetization entropy is **non-monotone** in both cases (2.683 → 2.707 → 2.585), and is exactly mirrored.

Coarse / local access exposes arrows but does not create them. Every arrow seen is **A_res-PRICED** (it needs reduced
or coarse access) **and** H_CORR-priced.

### D10 — Janus, and a structurally fixed special state

- **Real product state:** |S_S(t) − S_S(−t)| = 0.0 exactly. → **TIME ORIENTATION NOT SELECTED** (GTZ-type two-sided
  arrow).
- **Proposition 2:** any state fixed covariantly by H alone is stationary. The ground state has S_S = 0.406 constant.
  → **no arrow from an H-only boundary condition.**
- **Fixed by (H, Σ): the minimal-energy product state (mean field).**
  - E = −11.918 (ground: −12.442).
  - Found as **2 reflection-related minimizers**, so it is unique modulo Sym(H).
  - Both give Janus arrows: S_S = 0.566, 0.584, 0.325, **0**, 0.325, 0.584, 0.566 at t = −4 … 4, with half-chain
    entropy 1.033 at t = 4.
  - → **BOUNDARY CONDITION FIXED STRUCTURALLY BY (H, Σ) mod Sym(H); TIME ORIENTATION NOT SELECTED.**
- **Selector hostile.** Other (H, Σ)-covariant extremal principles pick *different* special states, all with overlap
  0.000 with the mean-field state:
  - the max-energy product state is almost stationary (S_S ≤ 0.009);
  - the min-variance product state gives a Janus arrow up to 0.903.

  So the middle condition is fixed by (H, Σ) **only after a selector functional is chosen**: preference-priced, as in
  S2-Σ.

## D12 — accounting summary

See `ledgers/ARROW_ORIGIN_LEDGER.md`.

**The minimal supplied item for every effective arrow found is:**

> **H_arrow = a low-correlation / independence / freshness boundary condition, relative to a factorization Σ, at a
> chosen time.**
>
> It is **not** low environment entropy:
> - a maximally mixed bath gives a Markovian arrow and records;
> - a 6.9-bit bath with correlations gives an anti-arrow.

**Status: S2-D-arrow COMPLETE.**
- Outcome A is excluded by theorem (finite closed unitary).
- Effective arrows need H_CORR + freshness, relative to Σ (B + C), plus reduced / coarse access (A_res).
- Typicality explains equilibrium, not the arrow (D).
- A special middle state can be fixed structurally by (H, Σ) + a selector functional. Never by H alone (Proposition 2).
  The time orientation is not selected (E).
