# What was checked and corrected during execution

- No M1 parameter, charge yield, source-pulse amplitude or mass tolerance was
  retuned. The separate protocols distinguish linear scattering, nonlinear packet
  controls and the later explicitly conditional cold-completion diagnostic.
- The initial dispatch's canonical decay width and low-energy neutral threshold
  were checked before a capture interpretation. The old kinetic reference already
  has the correct width; no old solver was changed.
- Two scratch traces were incomplete when read, including a previously generated
  pulse20 trace and one regenerated trace. No malformed row was skipped. Regenerated
  all required traces under a separate temporary work directory, checked every
  history column against published NPZ data (zero difference), then extracted
  spectra in the same sequential command chain. Original inputs remain unchanged.
- The first complex solve_bvp formulation exhausted its mesh budget. Rewrote
  the same linear boundary problem as real and imaginary equations; it converges
  with residuals below2e-8 and flux closure about1e-12. No physical threshold or
  boundary condition was relaxed.
- The first packet comparison applied the core response to the no-core control
  and clipped very small spectral tails to a local interpolation range. Corrected
  the control to vacuum and used a common explicitly covered momentum interval.
  The raw nonlinear histories did not change. The actual packet spectra, not
  only the central momentum, enter every comparison.
- A plotted ell4 resonance near k=.340 prompted a denser local momentum grid,
  a radial refinement and additional FFT padding. It changes the broadest-source
  covered average by about.6%, which is retained. The packet grids still differ
  by about4% relatively in their small fast-X probability; this is disclosed.
- The96-node angular verification missed the1e-10 check by1.7e-9 at high energy.
  Increased only the independent quadrature to192 nodes; the unchanged kernels
  then meet the numerical check. No science acceptance tolerance was loosened.
- One combined patch failed its exact-context check and wrote nothing. Reapplied
  the new cold-completion files and clock edit against the actual context before
  executing them. This was an editing failure, not a numerical result.

All calculations and separate-method checks in this packet are same-author.
Argon's supplied text is a proposal and algebraic commentary, not independent
execution of these new numerical runs. No external-model workflow was started.

Before committing the packet, the kinetic interpretation prompted a direct
inverse-channel check. Generalized only the incident column of the same matrix;
all archived forward probabilities recheck unchanged. Full open-channel
unitarity/reciprocity and momentum-squared detailed balance are checked. The
fast-X reverse source is now explicit in the report and handoff. The ideal cold
completion still assumes successful net processing and remains conditional.
