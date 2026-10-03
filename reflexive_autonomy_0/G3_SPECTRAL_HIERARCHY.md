# RA0 · G3 — CANONICAL ASYMPTOTIC SPECTRAL HIERARCHY (result)

> **Repaired by RA0 REPAIR 02** (RA2-01 … 05; `RA0_CORRECTION_LEDGER.md`). G3 is accepted provisionally after this repair.
> Where wording differs, the ledger takes precedence. Logs are kept as emitted; all numerical values are preserved.

**Central question.** Can a microscopic family define a complete multiscale hierarchy of dynamically separated sectors,
without choosing an ε, an observable channel, an eigenvalue cutoff or a preferred macro scale?

**Evidence:**
- `g3/g3_spectral.py` + `g3/g3_spectral.log` (independent code path, not independent reviewer; NUMERICAL ILLUSTRATION);
- propositions below (labelled).

## G3.1 Setup (allowed data only)

**Operator.**
- Positive relaxation operator: L_N = I − P_N for the discrete-time chains (A, B), or −G_N for continuous-time generators
  (C, D).
- All families are reversible, so the spectrum is real.

**Spectral data and degeneracy.**
- Ordered **distinct levels** 0 = λ₀ < λ₁ < λ₂ < … (multiplicities merged).
- Cuts sit only *between* distinct levels, so every projector E = 1_{[0, λ_k]}(L_N) is **basis-free**. Degeneracy needs no
  eigenvector choice.

**Cross-N firewall (G3.4).** Only the following are used:
- ordered levels;
- ranks (dim E);
- adjacent ratios r_k = λ_{k+1}/λ_k on the nonzero spectrum (the zero eigenspace is handled separately);
- spectral counting functions.

No eigenvector is matched across sizes.

**Divergence (G3.2).** A cut is **diverging** iff its ratio grows without bound along a **rank-indexed** sequence:
- ranks are allowed data, so this is permitted;
- no finite threshold or constant is used;
- a *large but bounded* gap is **not** a cut.

Identifying a slow projector with named structure (basins, density modes) is a **diagnostic using supplied labels**. It is
not part of the definition.

**Scope of canonicality [RA2-02].**
- *Earned, every finite N:* ordered spectral levels are intrinsic; the spectral projector at a specified gap is basis-free;
  degeneracies need no eigenbasis choice.
- *Conditional, across N:* a cross-N filtration is canonical only if the sequence of gap ranks is itself fixed by spectral
  data without externally matching modes. This is established in the constructed fixed-rank controls D1, D2 and the
  ladder with β > 2. **No general theorem for arbitrary families is claimed.**

## G3.5 Controls

### A — independent particles (N two-state particles)

| N | 2 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|
| largest adjacent ratio (rank) | 1.65 (3) | 1.65 (5) | 1.65 (7) | 1.65 (9) | 1.65 (11) |

**PROP G3-A.** The levels are 1 − θ^k (θ = 1 − a − b, k = 0 … N). The ratios (1 − θ^{k+1})/(1 − θ^k) ≤ 1 + θ are bounded
uniformly in N. ∎

**Terminal: NO NONTRIVIAL ASYMPTOTIC SPECTRAL SEPARATION.**

### B — SSEP ring (N = L/2; full many-particle kernel)

| L (states) | top adjacent ratios (ratio, rank) | spread λ_max/λ₁ | counting N(xλ₁), x = 1, 4, 9, 16 |
|---|---|---|---|
| 8 (70) | (1.619, 3), (1.248, 17), (1.209, 8) | 19.3 | 2, 16, 45, 66 |
| 10 (252) | (1.716, 3), (1.261, 5), (1.261, 21) | 36.7 | 2, 20, 60, 151 |
| 12 (924) | (1.775, 3), (1.315, 5), (1.182, 25) | 62.6 | 2, 24, 79, 229 |
| 14 (3432) | (1.814, 3), (1.275, 5), (1.187, 8) | 98.6 | 2, 21, 103, 330 |

**One-particle sector (exact formula), L = 100 / 1000 / 10 000:**
- largest ratio 3.996 / 4.000 / 4.000 at rank 2;
- spread 1.0e3 / 1.0e5 / 1.0e7;
- counting N(xλ₁) for x = 1, 4, 9, 16 is exactly **2, 4, 6, 8**.
- [RA2-03] For fixed x and large L, λ_q/λ₁ → q², so the ±q one-particle count approaches **2⌊√x⌋** (away from
  finite-size / Nyquist qualifications). At the tested perfect squares this is exactly 2, 4, 6, 8; it is **not** 2√x
  for general x.

**PROP G3-B (no canonical cut in the slow window).** For SSEP on the ring:
1. By self-duality, the span of {1, η_x} is invariant and carries the one-particle levels λ_q^{(1)} ∝ sin²(πq/L). This is
   KNOWN / standard (cited from the standard literature).
2. The many-body gap equals λ₁^{(1)} (Caputo–Liggett–Richthammer; abstract located).
3. Consecutive one-particle ratios are ≤ 4: sin t / t is decreasing on (0, π], so for 1 ≤ q < q + 1 ≤ L/2,
   sin(π(q+1)/L) ≤ ((q+1)/q)·sin(πq/L) ≤ 2·sin(πq/L), and squaring gives ratio ≤ 4.
4. So in the window [λ₁, λ_max^{(1)}] there is always a many-body level within a factor 4 above any level. **Every
   adjacent many-body ratio there is ≤ 4, uniformly in L.** ∎ (assembled from 1 – 3; not externally reviewed)

**Above λ_max^{(1)}** (the window [4/L, 2]): the bound is not proved. Numerically every ratio is ≤ 1.81 for L ≤ 14.

**Reading.**
- SSEP has an intrinsic **scale hierarchy**: the spread grows ~L². It is a **continuum / tower**, not a sharp split.
- The one-particle counting function has an L-independent **scaling profile** (2⌊√x⌋ [RA2-03]). That is a canonical cross-N object
  using counting functions only, and it encodes spectral dimension 1, consistent with frozen SCOUT S2-G2 (spectral
  dimension = CONDITIONAL DERIVATION from graph + probe).
- The **many-body** counting function grows with L at fixed x. The "slow sector below xλ₁" is not size-stable.

**Terminal: MULTISCALE CONTINUUM — NO CANONICAL CUT.**

### C — two-lane ladder (one particle), inter-lane rate γ_L = L^(−β)

Largest adjacent ratio @ rank:

| β | L = 50 | 100 | 200 | 400 | verdict |
|---|---|---|---|---|---|
| 0 | 3.98 @ 3 | 4 @ 3 | 4 @ 3 | 4 @ 3 | continuum (bounded) |
| 1 | 2.54 @ 3 | 4 @ 3 | 4 @ 3 | 4 @ 3 | continuum (bounded) |
| 2 | 19.7 @ 2 | 19.7 @ 2 | 19.7 @ 2 | 19.7 @ 2 | **large but bounded: no cut** (no-ε rule) |
| 3 | 986 @ 2 | 1.97e3 @ 2 | 3.95e3 @ 2 | 7.9e3 @ 2 | **diverging (∝ L): canonical projector of rank 2** (constants + lane imbalance) |

**Terminal: FAMILY / SCALING PRICED.**
- The same geometry has a canonical cut **only** for β > 2.
- Note β = 2: a visibly large, constant gap of 19.7 is **not** admitted. Admitting it would need a chosen C, which is
  exactly the A_resolution relocation the no-ε rule forbids.

### D — metastable positive controls

**[RA2-04] Interpretation.** D1 – D3 were **built** with metastable scale separation. They establish that the no-ε
spectral criterion (i) detects true asymptotic scale separation and (ii) returns the correct basis-free projector /
filtration in positive controls. They do **not** show that generic dynamics generate metastability: the rule detects
objectively supplied separation; it does not create it.

**D1 — 3 symmetric basins (complete graphs of size M), inter-basin rate ∝ 1/M.**

| M | 5 | 10 | 20 | 40 | 80 |
|---|---|---|---|---|---|
| top ratio @ rank | 4.00 @ 3 | 7.33 @ 3 | 14.0 @ 3 | 27.3 @ 3 | 54.0 @ 3 |

- The ratio ≈ 0.67 M, diverging. The rank-3 projector = constants + 2 inter-basin modes.
- **Terminal: CANONICAL SPECTRAL PROJECTOR.**
- Diagnostic: the projector equals the span of the basin indicators to 1e-15.

**D2 — nested: 4 basins in 2 super-basins (rates 1 : M⁻¹ : M⁻²).**

| M | 5 | 10 | 20 | 40 |
|---|---|---|---|---|
| top two cuts (ratio, rank) | (3.0, 2), (2.7, 4) | (5.5, 2), (5.1, 4) | (10.5, 2), (10.0, 4) | (20.5, 2), (20.0, 4) |

- **Two** diverging cuts, nested: rank 2 ⊂ rank 4.
- **Terminal: CANONICAL NESTED SPECTRAL FILTRATION.** The output is the whole filtration, not a preferred level.

**D3 — non-symmetric metastable control (random intra-basin rates, inter ∝ 1/M).**

| M | 5 | 10 | 20 | 40 | 80 |
|---|---|---|---|---|---|
| top ratio @ rank | 1.65 @ 3 | 2.80 @ 3 | 5.06 @ 3 | 11.06 @ 3 | 23.49 @ 3 |

The diverging cut is at rank 3. **G3 succeeds in the tested positive controls, which were built with a genuine scale
separation [RA2-04].**

[RA2-01] Each M in D3 uses a **newly generated** random matrix; the tested objects are not a nested realization of one
defined family, and no limit has been proved. G4 replaces D3 with a deterministic family rule (`g4/`).

## G3.6 Algebra test (spectral subspace ≠ observable algebra)

| subspace | contains constants | P-invariant | product-closed | generated algebra |
|---|---|---|---|---|
| D1 slow projector (symmetric basins) | yes | yes (spectral) | **yes** (residual ~1e-15) | = functions of basin label (proper) |
| D3 slow projector (non-symmetric) | yes | yes | residual 0.324 → 0.068 → 0.021 → 0.0067 → 0.0025 (M = 5 … 80). **[RA2-01] Numerical evidence for asymptotic algebraization, not a proof of convergence to zero** (fresh random matrix per M; residual measured on eigenbasis pairs, Euclidean norm) | at finite M the exact subspace is not an algebra |
| SSEP lowest eigenspace + constants (L = 10) | yes | yes | **no** (residual 1.000) | functions of the E₁ coordinates: 121 level sets of 252 states (**proper**), but **not semigroup-closed** (one-step lumpability defect 0.40) |

**Reading.**
- A canonical spectral projector was observed to approach an observable algebra **only** in the metastable controls
  (exactly for symmetric basins; numerical evidence only for D3 [RA2-01]).
- In the hydrodynamic case the lowest eigenspace is not an algebra, and the algebra it generates is proper but not
  autonomous. **ALGEBRA CLOSURE FAILED** (SSEP).

## G3.7 Terminal(s), per family (not collapsed)

| family | terminal |
|---|---|
| A independent particles | **NO NONTRIVIAL ASYMPTOTIC SPECTRAL SEPARATION** |
| B SSEP (hydrodynamic) | **MULTISCALE CONTINUUM — NO CANONICAL CUT**; ALGEBRA CLOSURE FAILED (lowest eigenspace); canonical counting-function scaling profile |
| C ladder | **FAMILY / SCALING PRICED** (cut only for β > 2) |
| D1 / D3 metastable | **CANONICAL SPECTRAL PROJECTOR** (algebra exact for symmetric basins; D3: numerical evidence for asymptotic algebraization, not proved [RA2-01]) |
| D2 nested metastable | **CANONICAL NESTED SPECTRAL FILTRATION** |
| overall [RA2-05, preserved] | **CONDITIONAL DYNAMICAL HIERARCHY DERIVATION.** When the supplied family carries a genuine diverging scale separation, the dynamics alone define a canonical (nested) projector hierarchy, with no ε, channel, cutoff or chosen level. When it does not (hydrodynamic continuum, independent particles), no canonical cut exists |

## G3.8 Selector firewall

Even the perfect filtration (D2) selects **none** of the following:
- Σ;
- A_partition;
- A_resolution;
- consciousness;
- orientation (all families are reversible, so the spectra are time-symmetric).

It derives dynamically distinguished **timescale structure**, conditional on the supplied family.

## G3.9 Consciousness entry condition

| situation | perspective-independent hierarchy supplied by G3? | consciousness question formulable without a supplied cut? |
|---|---|---|
| metastable-type families (D1 – D3, ladder β > 2) | **yes** (canonical projector / nested filtration) | **formulable** ("can an internally defined structure inhabit one level of this filtration without an external subsystem boundary?") |
| hydrodynamic continuum (SSEP), independent particles, ladder β ≤ 2 | **no** | **blocked**: "occupying a level" would require choosing a cut, i.e. A_resolution |

**G3 therefore does not justify a general consciousness layer.** It identifies the precise class in which the owner's
next question is well-posed. That class is fixed by the **supplied family**, which still has to carry the separation.
