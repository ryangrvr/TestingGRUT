# Independent review: M1 expansion, local source, and energy budget

Input parent `c4815306753ece00dad2b03412f0ed9fe89b8f92`. No change to M1
couplings, supplied yield, reference charge, or the 1% comparison. This is a
same-author calculation; no independent replication is claimed.

1. Verify the event-counter correction in `channels64.json`. Explain the
   difference among gross neutral events, net charge loss from the neutral
   channel, and forward-minus-reverse exchange. Check the stored stoichiometry
   against the parent matrix; do not infer gross cycles from endpoints alone.
2. Reproduce the `h=0.025,dt=0.003125` radial pulse through `m_X t=10000`,
   including the same-lattice stationary core at the measured charge, signed
   spectra and finite-window accounting. Is the residual near `197.54 m_X`
   supported, and what physical wave, matching, or off-shell components remain?
   Quantify separately finite-radius, time-step and spatial-grid errors.
3. Audit the seven computed `W(k,Q)` tables, especially the moving ell=4
   resonance, thresholds and reciprocity. Recalculate selected populated
   k/Q points with a second method. Verify that the sampled holdouts do not
   hide a large error in the **net** reaction rate when forward and reverse
   contributions nearly cancel.
4. Check the energy-space FRW equation for counts per core, the `tau(A)`
   mapping, the fixed-frequency exchange control and the explicitly measured
   interpolation energy defect. Reproduce 64/96 energy-grid and frequency
   spacing controls. Quantify how off-grid exchange diffusion changes tails;
   exact number and mean-energy conservation do not suffice.
5. Confirm the source-specific component ledger at `a/a0=10`: the fine-source
   96-node run gives cold excess 0.8052%, kinetic energy 30.30 `m_X`, held
   unassigned energy 199.37 `m_X`, and total modeled excess 0.9277%. State
   conditions under which total energy bounds all later cold compositions,
   including the pressure of an excited core and possible energy transfers.
6. Analyze the pre-mixing cloud. The simultaneous-source causal-ball fraction
   by `m_X t=10000` is at most 0.00509 at the supplied mean density, before
   allowing for the small initial core radius. Estimate or compute collisions
   with the actual radial/angular spectra, finite h decay length and return to
   the parent core. Do not apply the global mean-density collision clock to a
   single un-mixed plume. If the scalar h is an explicit kinetic species,
   remove the overlapping real-intermediate-state contribution from the
   reduced `X Xbar -> h -> X Xbar` kernel.

Report scripts actually run, physical assumptions and numerical limitations.
The experiment does not derive the initial pulse, core abundance, formation
epoch, a quantum lifetime, or an observational prediction.
