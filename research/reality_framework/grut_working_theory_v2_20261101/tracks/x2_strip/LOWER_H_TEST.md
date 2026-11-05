# Follow-up: compute the scalar lattice at the lower scale

This check is separate from the fixed M2-X verdict. The original `hLF55_mid`
history uses a lattice block computed at `8.4e10 GeV`; its kinetic block was
computed at `5.5e10 GeV`. Replace the reused lattice block with a fresh
HFa scalar lattice at `5.5e10 GeV`, after the two-node blind scan is complete.

Keep the exact grid, time step, vacuum seed and momentum cutoff from the
archived `KFmid_H84_nl` run. Keep its H84-calibrated `kappa_j(t)` table and all
other `lf3d` closure settings. This isolates the change in the physical scale
under the *same transferred closure*. It does not calibrate Pauli blocking
at H55. Use the unchanged scalar lattice engine through `t=3` and archive the
diagnostics. Compare the radiation fractions at `t=0.3, 0.5, 1, 1.5, 2, 3`,
maximum gradient energy fraction, first `f_coh < 0.1` time if present, energy
budget and the hand-over scale factor, H, curvature and spectator mass. No
cutoff for declaring equivalence is chosen. Report the differences.

This is a three-dimensional **scalar** lattice with the existing fermion
drain approximation. The independent full 3D Dirac calculation on matched
backgrounds at varied spacing is a different check. None of this follow-up's
rows can be added to the registered two-node verdict. A later full H55 relic
test must also rebuild the kinetic splice on this new lattice and state its
own analysis rule before running the scan.
