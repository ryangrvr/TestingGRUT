# P-02 — SF-1 CONTROL MAP (fermionic ingredient → bosonic replacement)

Extracted from the record, not reconstructed from memory: `SF1_FORMATION_CHARTER_01.md` §§1–7 (frozen
`8a4f71c`), `SF1_FORMATION_VERDICT_01.md`, `SF1_OWNER_RULING_02.md`, `SFG0_OWNER_RULING_01.md`.
**The only intended structural change is fermionic → bosonic statistics.**

| SF-1 ingredient (record) | P-02 main model | Changed? |
|---|---|---|
| Parent `H = −Σ_j (c_j†c_{j+1} + h.c.)`, ring `ℤ_L`, periodic, hopping 1, spinless, L even, full Fock space | `H_B = −Σ_j (b_j†b_{j+1} + h.c.)`, same ring, BC, hopping, full bosonic Fock space | **statistics only** |
| One-particle levels `ε(k) = −2cos k`, `k = 2πj/L` | identical | no |
| No chemical potential, no sector-specific term, one H builder | identical | no |
| Selector: exactly conserved `N = Σ n_j`, `[H, N] = 0` | identical (`[H_B, N] = 0`) | no |
| Reference state: lowest-energy eigenstate within sector `N` (a state declaration; no relaxation claim) | identical | no |
| Odd-N convention; A-ODD (symmetric fill, positive Fermi-edge gap, odd holes) | odd N kept. **A-ODD is fermion-specific** (it describes a Fermi sea) and cannot transport; replaced by **A-BOSE**: the sector ground state is unique and equals `(b_{k=0}†)^N|0⟩/√N!`, checked on exact Fock blocks | **replacement (forced by statistics; reason stated)** |
| Families: **D** `N = L/2` (L ≡ 2 mod 4), **E** `N = 1`; controls D-¼, D-¾ `N = L/4, 3L/4` (L ≡ 4 mod 8), E-3 `N = 3`, Ē `N = L−1`; C-6 `N = nearest odd √L` | identical sector rules and L sequences. **Ē and D-¾ lose their particle–hole meaning** (no PH symmetry for bosons); they are kept as literal N-rules | sectors transport literally; PH reading lost |
| Readout `ρ_q = Σ_j e^{−iqj} n_j`, `q = 2πm/L`; object `S_{L,N}(q, ω)` support, weight > 10⁻¹², grouped over degenerate eigenspaces | identical | no |
| Prescription P (L → ∞ at fixed q, then q → 0, `z_P`), P′ (`q = 2π/L`, `z_{P′}`); invariants I-z, I-q (soft set mod `q ~ −q ~ q+2π`); class = (I-z, I-q); SF-0 quotient | identical | no |
| V-FOCK (full Fock, L ∈ {6, 10, 12}) | **V-BFOCK**: exact bosonic sector blocks (occupation basis, `n_j ≤ N`), L ∈ {6, 10}, N ∈ {1, 3, 5}; bosonic Fock space is infinite, so the sector block is the exact finite object | replacement (bosonic Fock space is infinite-dimensional; the sector block is exact) |
| Terminal ladder: FORMATION-OF-LAW-CLASS / STATE-NOT-LAW / EDGE-ONLY-NONROBUST / SECTOR-SMUGGLED / UNFORMULABLE | mapped onto the auditor's P-02 labels: SAME-SPLIT / BOSON-COLLAPSE / NEW-BOSONIC-SPLIT / UNCONTROLLED | — |
| **Supplied inputs** (record): the sector (conserved boundary/state data); fermionic statistics (CA-1 F, Tier C, "conditional on supplied conservative fermionic statistics"); the reference-state declaration; the coarse-graining | same, with bosonic statistics supplied instead | — |

**Declared diagnostics (controls only; they never overwrite the main verdict):**

- **HCB (hard-core bosons)**: `n_j ≤ 1` with bosonic (commuting) operators. **NEW ASSUMPTION** (an
  infinite on-site repulsion). Purpose: separate **exclusion** from **antisymmetry**. Not a rescue of
  the free-boson result.
- **FPB (Fermi pattern in bosons)**: the free-boson eigenstate with one boson in each of the N
  lowest levels. Purpose: separate "which occupation is the ground state" from "which transitions are
  blocked". It is not a sector ground state, so it is not a reference state, by construction.
