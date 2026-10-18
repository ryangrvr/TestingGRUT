# Challenge, integration and review handoff

These are functions performed by one assistant with shared context, not four
independent agents or an external review. The uncertainty map was recorded
before the detailed proof artifacts, but the author had already developed
the analytic route. No blind Challenger or independent replication is claimed.

| Function | Challenge or action | Disposition |
|---|---|---|
| Builder | Specify fixed H, state, resource conservation and return mechanism | MODEL.md, RETENTION.md |
| Challenger | Can extensive outside dynamics wash out the isolated-cell return? | H_cut retains complete outside evolution; only five local hops differ |
| Challenger | Did zero-flux truncation or dropping interactions create the return? | All physical shifts/fields and density terms remain; empty guards make cross-density vanish in reference only |
| Challenger | Is a global O(N) perturbation treated as small? | No: cut norm <=5/48, independent of N |
| Integrator | Do preparation, measurement and resource count compose? | Initial record/fuel/program supplied; P=2(K-F) exact; no later postselection |
| Integrator | Is the old number-conserving P1 theorem silently reused? | No: q_x=n_x+2n_f,x supplies the new weighted moment and J_0=6 |
| Recorder | Preserve P1, freeze before controls, log source identity and output | SOURCE_PINS/EXECUTED_SOURCES/CHECKS and separate branch |

## Independent review priorities

1. Recheck the five-state cell and exact fuel Green function. Negative hopping
   and reverse conversion signs must match Gauss at both link boundaries.
2. Prove the guard-factorization claim with the full outside evolution left
   on. Outside charge configurations are spectators in the cut norm estimate;
   each corresponding kinetic matrix is a compression of ordinary six-site
   hard-core hopping. This is not a nonlocal primitive implementation.
3. Check that squaring the projected Duhamel amplitude is justified because
   the reference fuel is exactly full/empty at the designated times. Summing
   local probabilities uses expectation linearity, not physical independence.
4. Recheck all boundary drives in the moving-prefix comparison with fuel as
   part of the local resource. Conversion sources cannot occur in the empty
   tail under the static program, and full H retains them after excursions.
5. Check complete Neumann/parity remainders and energy-dependent corrections
   in the neutral bands. The first excited two branches require separate
   complement gaps above and below their energy, not the ground-band bound.
6. Verify the Sylvester estimate on adjacent-state low-band weight and its
   one-pair scope. It cannot be exported to a finite-density gas by relabeling.
7. Check the two-pair interaction bound. A diagonal Schur interaction is not
   an on-shell scattering amplitude; the packet leaves that open.

## What the negative conclusion does and does not cover

It covers the frozen EVERY-TIME threshold/window and the class containing
lambda=1, Omega=13/4. It does not exclude delayed windows, time averages,
other couplings, explicitly prepared clocks, enlarged autonomous reservoirs,
or a different microscopic law. Those are different hypotheses and require
new frozen campaigns. No assumption may be altered to repaint this gate.

Locality stays supplied at the generator level. A1's failure is conditional
on an already recorded path and prepared fuel; it neither refutes S9 nor
establishes autonomous formation from the original empty seed. The most
promising architecture is not inferred merely from one candidate failing.
