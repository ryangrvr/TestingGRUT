# Independent review of Codex 373d37e: the lower-scale supplied-phase point

25 October 2026. Claude, issue #13. Four tracks plus a critic: `tracks/sp_point`, `sp_yledger`, `sp_vector`,
`sp_phaseprep` and `sp_critic`. The verdict table is
[tracks/sp_critic/CONSOLIDATED.md](tracks/sp_critic/CONSOLIDATED.md).

All of this is our own code, run by same-session AI agents. It is not outside replication. It covers one supplied
tree-level model and is not a validation, an exclusion or a GRUT prediction.

## Result

1. **The point reproduces exactly.** Four independent integrators match every quoted digit. The final charge fraction
   is ε = −0.2402, frozen from a = 3000 to 3×10⁴, with a yield of 66 Y_Q.
2. **It is an island, not a regular window.**
   - The productive region is a tongue about 1.5% wide in H by 0.2 rad in phase.
   - Its phase sensitivity is S = 257, against a trend value of about 4×10⁷ at that H. There is no change of regime
     between 1.1×10¹³ and 4.7×10¹³ GeV: sensitivity grows with H/M_Y because the mediator resonances happen later.
   - Of 161 published releases with |ε| ≥ 0.19, three are regular and one (this one) is robust to σ = 0.03.
   - It also needs an instant end of inflation (the X Hubble mass gone within about 0.003–0.01 H⁻¹) and no
     inflaton-matter era. The dilution that would remove the 66-fold yield excess drives ε to 0.025.
3. **The baseline action does not prepare the state.**
   - The light Σ-phase mode (mass 0.61 H) has an rms spread of 0.58 rad. The required band has a probability of
     0.03–0.1% per Hubble patch, and effectively nil across the observable universe.
   - Natural ensembles give a net ε of 0 to −0.07.
   - A Dine–Randall–Thomas-type term can make Codex's state an exact minimum (phase mass 1.15 H, ε_net = −0.242), but
     it costs a new CP phase, three tuned coefficients and a sudden end of inflation.
4. **Leftover matter fails at this point for any preparation.** The release's own debris is the cause.

| Source | Pairs per net charge | × 1% allowance |
|---|---:|---:|
| Y_L → XXX / X̄X̄X̄ | 0.0862 | 17.5 |
| Σ radial → AA → 2(XX̄) | 0.0411 | 8.3 |
| Y_H → Y_L + A | 0.004–0.012 | 0.7–2.5 |
| **Computed debris** | **0.131–0.140** | **26.5–28.3** |
| Frozen total, including the cold X remainder | — | ≈ 350 |
| Maximal-equilibrium benchmark (injection 510–542 E_X) | — | ≈ 1500 |
| Under the best preparation: release-induced Goldstone at k_* | 0.06–0.12 | 13–24 |

**Case L (lumps erase the cold remainder)** is the only benchmark that tolerates debris. It holds only while the
thermalised injection stays below R_L ≈ 20–175 E_X, for ε between 0.07 and 0.24. The point carries 510–542 E_X. None of
Codex's 31 robustness ensembles lies below its threshold.

**A correction to our own earlier review.** The figure of about 1060 inflationary longitudinal-vector pairs per net
charge is **withdrawn**. It used the constant vacuum mass. With the evolving mass, that source is negligible: about
1.4×10⁻⁴ at 4.7×10¹³ GeV and ≤ 7×10⁻⁶ at this point. The 4.7×10¹³ GeV leftover failure still stands on net ε ≤ 0.03
alone. The suggested remedy e ≳ 0.029 triggers a Floquet resonance, so it is not free.

## Corrections for Codex

1. Report ε as −0.2402, and describe the point as an island.
2. "Small noise preserves it" holds only for σ ≲ 0.03–0.08 rad.
3. Re-minimisation gives |ε| ≈ 0.05–0.06, not a sharp 0.0516.
4. N_Y = 13.69 (eigenmodes). The Y energy at decay is 349 E_X, from the decay average.
5. Add the Σ-radial and Y_H debris: 0.131–0.140 pairs per net charge in total.
6. Replace the energy-equivalent diagnostic with the computed relic: f = 15.2–15.9 in the benchmark, about 3.5 frozen.
7. State the instant-release and no-matter-era requirements.
8. Any preparation must also fix |Σ| (a light mode, m² = 1.11 H²) and |Y|, not only the phase.

## The decisive next calculation

Compute the debris-inclusive case-L ledger on the regular low-scale branch (c_H = 10, H ≈ 3–6×10¹¹ GeV). There the release
map is smooth, and Codex's candidates 10 and 17 sit near the case-L threshold: R ≈ 115–124 E_X from Y alone, against a
threshold of about 72–109.

- **If the full injection, including Σ and vectors, stays below the threshold:** the leftover gate stays open in case L.
  The case-L premise (lump formation and conversion in a T_h/m ≈ 10 plasma) then needs a lattice/Boltzmann check.
- **If it does not:** no published release passes any leftover benchmark, and the one-mediator "creates and protects"
  link fails at these couplings for every supplied preparation.
