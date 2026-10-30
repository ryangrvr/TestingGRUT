# Sources and provenance

The parent repository commit and exact input files are pinned in [INPUT_PINS.json](INPUT_PINS.json). Numerical outputs, scalar lower bounds, the Option 1B counterexample, implementation choices and reported control results in this packet are our calculations.

* [Ma & Bertschinger (1995), Cosmological Perturbation Theory](https://arxiv.org/abs/astro-ph/9506072): standard conformal-Newtonian Einstein, fluid and collisionless hierarchy conventions. The paper's full species treatment is broader than the simplified test cosmology here; it does not validate our formation or BGK assumptions.
* [Yunis, Argüelles & López Nacir (2020), Boltzmann hierarchies for self-interacting warm dark matter scenarios](https://arxiv.org/abs/2002.05778): primary treatment of momentum-dependent elastic collision hierarchies and relaxation-time approximations. Abstract and scope consulted; their exact kernels are not imported as a validation of our chosen BGK closure.
* [White, Minamitsuji & Sasaki (2012), Curvature perturbation in multi-field inflation with non-minimal coupling](https://arxiv.org/abs/1205.0656): primary context for the distinction between scalar-field fluctuations and cosmological curvature/isocurvature observables. Abstract and scope consulted; no claim to reproduce that paper's calculations.
* [CLASS source commit](https://github.com/lesgourg/class_public/tree/814ce59c380a5cccd0bed1747a384ff7bfc1f9cf), v3.4.1: external cosmological solver used only for the stated collisionless controls. Configuration, successful spectra and failure tails are retained.

No author above endorses GRUT. No new observational likelihood or empirical validation is performed in this packet.
