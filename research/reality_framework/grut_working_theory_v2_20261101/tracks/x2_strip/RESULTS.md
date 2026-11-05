# M2-X: no strict strip within the simulated range

4 November 2026. Both registered grids are complete: 135 keys at each of
`hLF84_mid` and `hLF_mid`. The 13 already finished rows at each node are
preserved exactly; 244 missing rows were computed with the original engine
and scan code. The grids were sealed before evaluating any roots. The
acceptance rule is the unchanged [pre-registration](PREREGISTRATION.md)
from `6d836ca`. The result tests the declared lattice-fermion closure and
relic prescription on this grid, not a model-independent dark matter bound.

**Strict result: 0/135 M2-X points, hence no strip.** There are 37
self-consistent free-X roots inside the two-node bracket before applying
the threshold band. Of those, 32 have no crossing, three have a crossing
at just one node, and two have a numeric interpolated bulk `e_crit`.
Those two thresholds are `2.18356e-5` and `2.52898e-5`. Even the lower is
5.30 times the registered upper limit `4.12e-6`. Thus free-X roots exist
in the bracket, but none has the required transition strength there.

The Q-matter calculation uses the separate registered condition:
self-consistent, bracketed Q roots. It gives two isolated points:

| c0 | delta (degrees) | H* (GeV) | S/Q at e=1e-6 | Bulk e_crit |
|---:|---:|---:|---:|---:|
| 1.380 | 0 | 8.49995e10 | 3.33421 | 8.35670e-7 |
| 1.395 | +2 | 8.49152e10 | 3.37198 | 8.87874e-7 |

Neither Q root is in the M2-X threshold band. The table does not establish
that the real universe forms a Q-matter population or that either root
survives a more complete fermion closure.

**Outside the bracket.** Extrapolating only the computed pair gives 15
self-consistent free-X roots in the band, all below the lower computed node,
at H* from `7.34632e10` to `8.08294e10 GeV`. This is a diagnostic for where
additional lattice coverage matters. It is excluded from the strict result.
The fresh four-node extended scan, with the two rescaled lattice proxies,
has not been run. Its histories are restored and archived separately.

**Replay and recovery.** Four fermion table summaries and all six checks
on each of the four rebuilt histories match their archived values exactly.
Five baseline scan quantities on each computed node match with maximum
relative difference `2.42e-16`. An interrupted scratch trace stopped part
of the first pass after 226/270 rows. A fresh scratch directory reused the
verified input tables and every completed row. All 200 rows computed after
Claude's pause and before that restart were compared and retained exactly;
only the 44 unfinished rows were recomputed. No unreadable trace was treated
as a failed physical point or dropped from the grid.

The machine-readable results are in [blind_verdict.json](results/blind_verdict.json)
and the two `completed_scan_*.json.gz` archives. [REPLAY.md](REPLAY.md) gives
the commands, and [the history provenance](results/restored_histories/provenance.json)
records the input hashes, checks and runtime. [LOWER_H_TEST.md](LOWER_H_TEST.md)
specifies the separate lower-scale scalar lattice check.

## Separate lower-scale lattice control

The preregistered control ran the `64^3` scalar lattice at `5.5e10 GeV`
through `t=3.0000625` (9204 steps), retaining the H84 drain table, grid,
time step, cutoff and seed. Both histories are compared using their last
diagnostic near `t=3` (the H55 row is at `2.9980937`). The main differences:

| Quantity near the hand-over | Fresh H55 lattice | Reused H84 lattice |
|---|---:|---:|
| Radiation energy fraction | 0.669585 | 0.832143 |
| Chi condensate energy fraction | 0.317327 | 0.157878 |
| SM energy fraction | 0.031440 | 0.068986 |
| Cycle-mean R / H_i^2 | 0.026988 | 0.013215 |
| Mean spectator m_s^2 / H_i^2 | 9.65821 | 9.15866 |
| Maximum absolute budget diagnostic | 0.013304 | 0.009704 |

The radiation difference is −16.26 percentage points, and the surviving
chi fraction is 2.01 times the proxy value. The early radiation ordering
reverses: at `t=0.5` it is 0.587 (H55) versus 0.443 (H84). This is a real
change in the computed hand-over on the held drain prescription. It does
not tell us the late K or the relic branch on the new history. The H55 run
has one scalar resolution and no recalibration of its fermion drain.

The [lower-scale result](results/lower_H_lattice.json) and
`results/lower_H_diag.txt.gz` preserve the run. The engine completed normally;
an initial export guard had incorrectly assumed diagnostics were written
every step. The corrected guard checks the actual completion flag and the
six-step output cadence. No lattice dynamics or parameters were changed
for that export correction.

The next numerical dependencies are a kinetic splice on an actually
computed lower-scale history, and the distinct full 3D fermion spacing
test. The latter has not been run. The existing `lf3d` scalar lattice plus
line fermions cannot stand in for a full Dirac field calculation. For scale,
the reported H84 hot-spot comoving momenta reach `1.2e5 H_i`; representing
those in a uniform periodic `L=0.42/H_i` box would need about `1.6e4` cells
per axis, even before a convergence margin. A 64-to-96 scalar-grid change
alone does not solve that resolution problem.
