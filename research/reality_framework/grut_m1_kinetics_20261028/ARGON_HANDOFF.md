# Independent review target: M1 fixed-profile transport

Input: bf8ab3b96c07691410145c3b8c21ca18b009d66d. New packet:
research/reality_framework/grut_m1_kinetics_20261028/.

Ryan's latest shared page has been read. This is a reconstruction from committed
inputs; the stalled implementation itself was not accessible. No review request
has been sent and no independent-author numerical replication is claimed.

Keep M1 fixed. The corrected canonical radial off-diagonal entry is sqrt(2)FS,
as already implemented in bf8ab3b. The width is Gamma_h/m_X=0.0001875658992.

At optical time tau=1000, the primary graded translated grid gives:

- Initial Nbar=1695.0095; final Nbar=0.558941.
- Final Nplus=11629.3892 and Qcore=188371.1698.
- Gas kinetic energy=2526.2895m_X; 8.17589 positive particles remain above
  the reverse threshold. No heat has been removed.
- Re-solved stationary cold-mass excess=0.726635%, conditional on later cooling
  with unchanged composition. It is not a final mass-budget pass.
- An additional 419.35m_X of original pulse/core/matching energy is explicitly
  retained but its release has not been evolved.

Please check the following independently, preferably with another implementation:

1. Derive the isotropic resonant collision coefficient
   pi g_h^2/(8 E1 E2 k1 k2) and the kinematic condition
   |(E1 E2-2)/(k1 k2)|<1. Check its equality to thermal inverse-decay reaction
   density and that the real h intermediate state has been counted once.
2. Check the identical-charge amplitude 4lambda_X+g_h^2/(t-6)+g_h^2/(u-6),
   the identical-final-state 1/2, and the separate unordered-incoming-pair factor.
3. Check that both boosted daughters are deposited with the correct lab energy,
   and that omitted support is recorded rather than turned into cooling.
4. Reproduce the measured pulse source, the primary translated-grid refinement,
   the no-diffusion two-state exchange control, and the gas-only thermal residual.
   Passing conservation alone is insufficient. Earlier unaligned-grid histories
   and the failed endpoint handling remain available for comparison.
5. Audit the core-energy remainder against the original pulse energy and verify
   the finite-Q stationary masses. The fixed profile is an assumption: sampled
   cross-section changes reach 7.5% at the final charge. Do not label this
   calculation a self-consistent expanding-universe solution.

Report which scripts and data you actually executed, numerical tolerances,
differences, and any failure. Analytical verification should remain explicitly
distinct from numerical replication. See REPRODUCE.md for commands and
INPUT_PINS.json for input hashes.
