# mf_hybrid results: status of the exports

2 November 2026 · Claude (mf_hybrid). All values are model choices; see ../NOTES.md.

| file | status | content |
|---|---|---|
| bench_HFa.json | **final** for its model choices | hybrid (m_χ, m_σ) = (250, 600) H_i, n_s = 0.9649 (Planck-centred), spectator mass from the waterfall (m_s = g_s χ), axial singlet-neutrino reheating, K = 2.048 |
| bench_HFb.json | **final** for its model choices | same sector, n_s = 0.9743 (ACT-centred), spectator mass from ξ_s R (ξ_s = −3/4), K = 2.047 |
| axial_check.json | final [C] | exact half-period fermion production against Landau–Zener |
| axial_tables.json | final [C] | ε(A) production tables used by the reheating rates |
| diagnostics.json | final [C/E] | defects, waterfall ζ, R spectrum, decay channels, gravitational production, spectator |
| structure.json | final [D/C] | parameter relations, one-loop e-fold capacity scan, perturbative-reheating mass bound |

**What "final" means here.**
- The CMB numbers (exact MS), the waterfall, the defect check and the homogeneous history to t = 0.6 H_i⁻¹ are [C].
- The reheating history after t = 0.6 is a fluid model whose rates are [E], built on the [C] half-period production. It has
  no Pauli blocking, fragmentation or N thermalisation dynamics.
- K, T and the plasma time inherit that [E].

**Interim exports.** None remain. An earlier HFa export, with K = 9.8 (two χ species at f = 4000 H_i), was superseded
and overwritten.

**Slice convention.** a = 1 and H = H_i at σ = σ_c (end of slow roll, R = 12H_i²). ε_H = 1 follows at t = 0.19 H_i⁻¹, with
H = 0.997 H_i.
