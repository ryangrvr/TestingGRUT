# Protocol amendment 01 — before any physical trajectory

**18 October 2026**

The original `PROTOCOL.md` is retained. Before executing a physical evolution, the Constructor supplied an analytically controlled perturbation of the same action: deform both archived profiles by

\[
M=\operatorname{diag}(e^\epsilon,e^{-\epsilon},1),\quad
A_0(y)=F(|My|),\quad S_0(y)=S_*(|My|).
\]

Use this family instead of the originally proposed quadrupole-amplitude/phase formula. The original family was never evolved; there is no discarded unfavorable trajectory. This is a pre-execution improvement to preparation control, not a change of physical law or parameters. All seven run names, grids, epsilons, time steps, horizon and diagnostics remain as specified.

Because `det M=1`, the continuum charge norm, potential energy and rotation energy are exactly unchanged for `A_dot=−i Omega A`. For the spherical reference, angular integration gives

\[
\Delta\bar E=\left[\frac{\operatorname{tr}(MM^T)}3-1\right]T_{grad}
=\frac43T_{grad}\sinh^2\epsilon,
\quad T_{grad}=66.558341736314.
\]

This is a smooth, genuinely three-dimensional nonradial distortion with known injected energy. No phase kick is included in this amendment. Normalize `Omega_grid` on the actual periodic grid to the same target Q, compute the matching spectral energy, and compare it with the analytic continuum preparation energy. The difference is a separate spatial/preparation discretization diagnostic, not physical dissipation.
