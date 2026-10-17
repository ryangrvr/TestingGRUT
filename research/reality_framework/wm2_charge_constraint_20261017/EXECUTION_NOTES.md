# Execution disclosure

The initial diagnostic code used a linear spread of storage detunings
inside the larger absolute box |delta_i|<=3/4. That spread need not fit
the frozen common-plus-independent-disorder family. Before publication,
the default was corrected to bar_delta=.23 with independently specified
offsets spread across [-.25,.25], and every archived control was rerun.
The final CHECKS.json records those detunings. No candidate coupling,
constraint, theoretical bound, acceptance gate or physical receiver was
changed. Earlier scratch diagnostics are not a second formation protocol.

The implementation also tests the unconditioned local receiver alongside
the prefix-weighted readout, and adds separate nearby charged protocols.
These address the existing frozen internal-channel and viable-matter
requirements. There are no formation runs or fallback trajectories.

Final results use ordinary double precision, sparse matrix exponentials
and source-pulse finite differences. Integer constraints and rational
constant checks are exact. No interval enclosure of all floating-point
operations or independent human proof review is claimed. The source
hashes in EXECUTED_SOURCES.json identify the final executed implementation.
