# Sources and read scope

- The actual implementation is the parent
  [radial_nodal.cpp](../grut_core_conversion_20261023/code/radial_nodal.cpp).
  Its velocity multiplier, cubic layer, wall, probe radius and output
  normalization were read directly. Its six trace hashes and run parameters
  are pinned in TRACE_AUDIT_RESULTS.json. Reusing these author runs is not
  independent replication.
- The previous
  [exterior projection](../grut_core_conversion_20261023/code/exterior.py)
  and [field spectra](../grut_daughter_core_20261023/code/spectra.py) supply
  the conventions and disclosed local-free-field cut. Their conditional
  particle interpretation is not promoted here.
- Steven G. Johnson, [Notes on Perfectly Matched Layers](https://math.mit.edu/~stevenj/18.369/spring14/pml.pdf).
  Read the retrieved PDF text for the absorber/interface distinction,
  homogeneous linear exterior assumptions, complex coordinate stretch and
  discretization limitations (sections 1, 3 and 4). The massive-wave formulas
  in this packet are explicitly derived here; no Maxwell or massless-wave
  time-domain implementation is silently transferred into GRUT.

This is an author numerical-boundary investigation. No experimental dataset,
cosmological likelihood, new installation or external reviewer execution
was used.
