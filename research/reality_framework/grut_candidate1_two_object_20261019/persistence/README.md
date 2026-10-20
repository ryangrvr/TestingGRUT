# Post-collision persistence and azimuthal response

This is an additive follow-up on the isolated two-object research branch. The physical scalar couplings and original kernel are unchanged. A generated extension synchronizes linear m=1 and m=2 azimuthal perturbations with the same background Verlet integrator. Perturbations do not backreact on that background.

The protocol is frozen before the time-160/240 evolutions. A short timing pilot and successful implementation controls preceded the freeze and are preserved. The selected canonical norm can grow through translations, deflections and phase shifts; growth alone is not a nonlinear stability or chaos theorem.

Run `python3 build_kernel.py`, `node controls.cjs`, then `node run.cjs RUN_ID` for an id in protocol.json. `kernel_origin.json` pins the source identities. The final report will identify completed outcomes and their numerical limits.
