# Bounded radiation seed audit

Freeze the 1 meV transfer, \(r_{0.05}=0.002765\), the fitted inflation
parameters and both \(w_s=3/5\) smooth conversion histories in the
[replication packet](../replication/README.md). The audit diagnoses the
local fluid closure, with \(f=\rho_r/\rho\), \(Q=(\sigma-0.8)H\rho_s\rho_r/\rho\),
and \(df/dN=\sigma f(1-f)\). Use the end of each inflation blend
\(N_i=\delta\) for the seed inventory. There is no new parameter search
or observational holdout in this audit.

1. Recover \(f_i\) and \(H_i\) from the frozen input and give the energy
   density, energy and causal-quantum count in a Hubble sphere. For
   conventional on-shell radiation compare \(E_{\mathrm{sphere}}/H_i\)
   with one. Also compare the initial inventory in the causal volume
   of the nominal transition and on the pivot scale. Treat the
   equilibrium-temperature conversion as formal.
2. Keep the calculated hard tensor and vector populations separate
   from the putative visible bath. Take their resolved \(8\le k/k_{\rm end}
   \le16\) band at \(N=4\). In a completion with one minimally gravitating
   massless Maxwell field, evaluate the phase-averaged kinetic estimate
   for two hard longitudinal vector quanta producing a Maxwell pair
   over \(4\le N\le4.1\). Use only this band and relative incident-angle
   cosines in \([-0.5,0.5]\). Report the direct injection before any
   nonlinear amplification. Keep the bath identification and
   thermalization assumptions explicit.
3. Evaluate the exact long-wavelength separate-universe response to
   independent initial logit perturbations on a uniform-total-density
   slice. Compare the derivative with \((0.2/\sigma)(1-f_i)\).
   Illustrate root-power amplitudes \(10^{-4},10^{-3},10^{-2}\) at
   the frozen pivot. Use 10% of the frozen scalar power as a comparison
   budget, with zero cross correlation for these examples. This is
   a consistency diagnostic for the frozen calibration.
4. Distinguish an effective homogeneous history, a kinetic seed, and
   a microscopic reheating model. A kinetic seed that cannot populate
   a causal volume leaves the fluid realization open. A produced
   species enters \(Q\) only when the completion identifies it as
   the relevant component.
5. Validate the exact reproduction block from the distributed ZIP in a
   fresh Python environment. Use the extracted packet and its declared
   dependencies. Record whether a file correction is needed and verify
   all eight mode outputs and both comparison ledgers.

The production check is restricted to the stated field and state
assumptions. In particular the ideal thermal degrees of freedom in
the transfer calculation do not specify curvature couplings, a visible
vacuum, or the source of the bath.
