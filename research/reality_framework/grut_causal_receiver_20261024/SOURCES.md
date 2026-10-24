# Sources and scope

- Xavier Antoine and Xiaofei Zhao, [Pseudospectral methods with PML for
  nonlinear Klein-Gordon equations in classical and non-relativistic regimes](https://arxiv.org/html/2101.02528v1),
  arXiv:2101.02528, later JCP 448 (2022), 110728. Read the retrieved
  introduction and section 2: causal frequency shifting, auxiliary-memory
  realization, absorption profiles, numerical controls and the limitation
  of exact matching to the linear problem. This is established numerical
  analysis. This packet derives and tests its own staggered-grid realization;
  it does not copy that paper's pseudo-spectral integrator or import its
  nonlinear error curves as a GRUT certificate.
- The [parent boundary audit](../grut_soft_receiver_20261023/REPORT.md)
  fixes the failed ordinary-damper comparison and its action/trace conventions.
- The [parent radial action and preparation](../grut_core_conversion_20261023/code/radial_nodal.cpp)
  and [initial data construction](../grut_core_conversion_20261023/code/refinement.py)
  were read directly. The new spatial gradient operator is identified
  explicitly rather than mislabeled as the same discrete Hamiltonian.

An institutional copy of the separate Baara-Diaz-Tlemcani paper was sought
but its document endpoint returned a bot-check page; no full-text read or
implementation based on that inaccessible paper is claimed. No new external
package, experimental dataset or independent reviewer execution was used.
