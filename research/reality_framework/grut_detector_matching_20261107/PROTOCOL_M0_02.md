# M0-02 — a separately frozen endpoint-matching hypothesis

M0-01 missed its registered windows. Its data are known and preserved. This follow-up changes the finite-a boundary dictionary on an analytic basis; it is not the original test and is not an untouched test of a rule chosen independently of M0-01. No continuous parameter is fitted to the histories. Freeze this card and implementation before calculating the new histories.

## Derivation and hypothesis

The free open hopping chain has single-particle levels

$$
\varepsilon_\ell=a^{-1}\cos\frac{\pi\ell}{N+1},\qquad \ell=1,\ldots,N.
$$

At half filling its lowest neutral excitation is

$$
\Delta_{\rm hop}=\frac2a\sin\frac{\pi}{2(N+1)}
=\frac\pi{(N+1)a}+O\!\left(\frac1{a(N+1)^3}\right).
$$

Thus the long-wavelength free chain measures the distance between its two ghost endpoints, `(N+1)a`, whereas M0-01 fixed the occupied-site span `(N-1)a`. The difference is a leading order-a boundary correction. **Hypothesis:** fixing the ghost-endpoint distance also improves the matched Schwinger histories after the declared bulk mass improvement. The free-chain identity does not prove this for the interacting model.

Fix the same target length `L=5`, but set `a=5/(N+1)`. Interpret sites as `z_j=(j+1)a` and link midpoints as `z_(j+1/2)=(j+1.5)a`. Keep the standard Gauss law, theta zero, even neutral sector, `m_lat=-a/8`, target boundary values `phi/sqrt(pi)=1/4`, Gaussian smearing `sigma=0.8`, strength 2, phase 0, and all M0-01 times. The continuum mean and covariance are unchanged. Link quadrature omits the unoccupied edge regions; this remaining detector-discretization error is not assumed absent.

The primary observable remains the weighted raw link field. A predeclared control at N=10 and 20 removes only its known alternating c-number:

$$
L_j=\tfrac12\sum_{i\le j}Z_i+\tfrac14[1+(-1)^j],
\qquad
F_{\rm smooth}=F_{\rm raw}-\delta_N,
\quad \delta_N=\tfrac14\sum_j f_j(-1)^j.
$$

Implement this control by an instrument phase `-c*delta_N`, not by fitting a measured mean. It does not replace a failed primary verdict.

## Frozen execution

Use primary N=10,14,18,20 and a tighter N=14 replay. The controls are the two specified smoothed observables; bare-zero results from M0-01 are retained but not rerun here. All original M0-01 numerical tolerances, TV windows, independent reference, maximum mode count and resource caps are unchanged and encoded in `protocol_m0_02.json`. The same analysis is reused without new fits or extrapolation. The wrapper reuses the frozen lattice algorithm at an occupied-site span `5*(N-1)/(N+1)` and translates the origin; the relative Gaussian positions coincide exactly with the ghost-endpoint coordinates above.

Before production, check the analytic free hopping gap against independent one-body diagonalization for the four N values to absolute error below `1e-12`. This checks the stated identity, not the interacting hypothesis.

Passing would support the improved finite-a dictionary for the declared massless ideal-instrument protocols. It would not prove a unique physical endpoint convention, a massive error bound, a continuum supremum theorem, causal implementation of the extended detector, or new GRUT dynamics. A missed gate remains a miss and is not repaired with the smoothed control. Zenodo remains paused.
