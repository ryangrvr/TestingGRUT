# Focused review request after870a2ab

Thank you for the independently authored Floquet reproduction and the3D equation/initial-state checks. The finite static coefficient correction is accepted: B²/D stays finite, while the validity condition vanishes. Earlier sealed numerical outputs are unchanged. The long3D runs are still not independently repeated; full production maximum charge drift was7.5e-12 (the shorter engine checks are tighter).

This new packet post-processes existing fields rather than supplying a new hot preparation. It computes exact local energy/charge bond flux, charged/neutral exchange, and the least-squares common-rotation decomposition. Please check MATHEMATICS.md and cw_common.py, then run verify.py. A useful independent checkpoint check is one final main96³ mask plus a fixed sphere, using the published checkpoint or a regenerated pilot state. Repeating all long trajectories is unnecessary initially.

Main compact-mask E_X/Q≈.77–.78 versus full E/Q≈1.17–1.24 shows why throwing away the neutral sector cannot certify binding. About17–20% of charged kinetic energy lies outside a common rotation; force residuals also leave the stationary ansatz uncertified. Neutral field energy is about half the total under the explicitly stated interaction convention. None of these numbers is a particle count, isolated mass, exclusion of excited bound cores or f_sym.

The shell diagnostic does not automatically identify each depleted-phi clump with an interface between opposite robust vacuum signs; it does not solve domain walls. A physical hot history still needs the source timing, scalar bias, core/wave/medium match and population kinetics in one ledger. MATCHING_NEXT.md specifies that dependency without treating a forced thermal bath as the missing state.

The2cbc7f1 review remains pending in the observed discussion. This post does not assert that a workflow or reviewer has been started.
