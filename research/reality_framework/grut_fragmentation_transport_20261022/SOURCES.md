# Sources, credit and retrieval scope

22 October 2026 (America/Chicago). Primary sources; novelty is not asserted.

## Formation

* R. Friedberg, T. D. Lee and A. Sirlin, **Class of scalar-field soliton
  solutions in three space dimensions**, Phys. Rev. D 13, 2739 (1976),
  [doi:10.1103/PhysRevD.13.2739](https://doi.org/10.1103/PhysRevD.13.2739).
  Model-family credit already present in the unified-working-theory packet.
  Not newly rederived as a GRUT mechanism; no new paper read is claimed here.
* A. Kusenko and M. Shaposhnikov, **Supersymmetric Q-balls as dark matter**,
  hep-ph/9709492v3, Phys. Lett. B 418, 46 (1998),
  [primary PDF](https://arxiv.org/pdf/hep-ph/9709492v3).
  Read abstract, formation introduction and instability/fastest-mode discussion
  including the numerical linear-amplification estimates. Credit for condensate
  fragmentation as a cosmological route. This paper does not supply a final
  full nonlinear distribution for this action.
* K. Enqvist and J. McDonald, **Q-Balls and Baryogenesis in the MSSM**,
  hep-ph/9711514v5, Nucl. Phys. B 538, 321 (1999),
  [primary PDF](https://arxiv.org/pdf/hep-ph/9711514v5).
  Read abstract/introduction and section 4.3 on negative pressure and growth.
  Credit for the instability-to-fragments reasoning. Its potential, primordial
  preparation and charge estimates are different. It is not our lattice engine.
* S. Kasuya and M. Kawasaki, **Q-ball Formation through Affleck-Dine
  Mechanism**, hep-ph/9909509v3, Phys. Rev. D 61, 041301 (2000),
  [primary PDF](https://arxiv.org/pdf/hep-ph/9909509v3).
  Read initial setup, dimensional qualification and nonlinear spectra/size
  discussion. Credit for nonlinear 1D/2D/3D lattice formation; planar/string
  runs cannot be reinterpreted as 3D object distributions.
* S. Kasuya and M. Kawasaki, **Q-ball formation in the gravity-mediated
  SUSY breaking scenario**, hep-ph/0002285,
  [primary record](https://arxiv.org/abs/hep-ph/0002285).
  Bibliographic follow-on identified; full calculation not read or reproduced
  here. We do not consume a numerical coefficient from this paper.
* T. Hiramatsu, M. Kawasaki and F. Takahashi, **Numerical study of Q-ball
  formation in gravity mediation**, 1003.1779v2,
  [primary full text](https://arxiv.org/html/1003.1779v2).
  Read introduction, setup/resolution criteria, discussion and numerical/
  identification appendices. Credit for resolution-dependent charge-distribution
  cautions, higher-order integration and object identification. Our code instead
  uses spectral Strang splitting and simple periodic connected components;
  those components are not identified as relaxed Q-balls.

The MSSM/AD mechanisms in these papers are not additional operators silently
inserted into our potential. None fixes the supplied primordial charge here.
The preceding packet explicitly credited Hardy et al. 1504.05419 and
Gresham--Lou--Zurek 1707.02316 / 1805.04512 for composite detection and
formation context. Those citations and the earlier detector qualification stand.

## Spectral distortion and cosmological response

* J. Chluba, A. Ravenni and S. K. Acharya, **Thermalization of large energy
  release in the early Universe**, 2005.11325v2,
  [primary PDF](https://arxiv.org/pdf/2005.11325v2), MNRAS 498, 959 (2020).
  Read abstract and visibility discussion/equation (49), plus its large-release
  qualification. The historical FIRAS energy screen ~6e-5 and small-distortion
  Jbb=exp[-(z/1.98e6)^(5/2)] come from this account of established methods.
  This is not a latest-limit fit, and the large-energy cases are not precise
  spectra computed with its nonlinear thermalization solver.
* Background point H0=67.4, Omega_m=.315, T0=2.7255 K and N_eff=3.046:
  adopted illustrative standard-cosmology inputs, not new measurements or a
  refit. Omega_L is supplied by flatness, not derived by GRUT.
* Scalar lifetimes and the two photon-response envelopes are consumed from
  the pinned `grut_portal_empirical_test_20261022/HISTORY_RESULTS.json`
  at parent bd2cac3. We do not claim a new SM scalar-decay precision calculation.
* Binding budgets, supplied charge, linear-onset scales, wall qualification
  and static-profile inputs are consumed from the pinned
  `grut_formation_wall_audit_20261022/` packet at parent bd2cac3.

The explicit post-annihilation entropy/photon conversion, collisionless
decay/redshift kernels, threshold amplitudes, finite-momentum angular
integrals and irreversible-channel example are calculated in this packet.
They are standard classical/kinetic methods applied to the declared model;
no literature priority claim accompanies them. The earlier BBN nuclear
calculation is not rerun and is not treated as covering the formation source.
