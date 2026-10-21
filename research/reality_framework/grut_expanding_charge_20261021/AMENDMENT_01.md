# Pre-production numerical clarifications

2026-10-21, after protocol freeze `188b92b`, before any production run.

1. The Gaussian Fourier polynomial contains 64 distinct nonzero modes, with wavevectors
   drawn with physical comoving rms components 1/ell and cutoff |k|<=2.5/ell. Cosine and
   sine coefficients are saved and each field is normalized analytically to unit rms.
   This is a supplied perturbation family, not a thermal or quantum spectrum.
2. The formation criterion in the protocol is still the raw S>0.5 and S<-0.5 fractions.
   To avoid calling an excursion within a restored X-rich core a genuine vacuum, the
   topology masks and basin seeds additionally require physical F^2<0.5. Both raw and
   genuine-vacuum fractions are retained. Basin labels extend from these seeds by
   six-neighbor graph distance; the seed order is geometric and Z2-equivariant.
3. Topology cadence is 4.8 rather than 4 conformal units so every dt in the suite divides
   it exactly. Scalar cadence is 1.2. The endpoint is always included.
4. Core coarea and core-energy estimators are retained, but a sharp |S|<0.5 mask has
   sizeable grid-position errors in the coarse boosted planar control. Add a smooth
   coarea estimator with delta(S)=3(1-S^2)/4 on |S|<1. A smooth phi-energy proxy uses
   weight 1-S^2 and divides by 0.8, the exact planar-kink energy fraction. Both calibrate
   on an unoccupied planar kink; neither is asserted to be the full energy of a charged
   network. The all-field density, pressure and conformal work remain the conserved
   energy accounting. No unique wall/receiver partition has been derived.
5. First validation had one failed exact-mirror check: a 3.66e-13 field difference.
   Parallel reduction order in the initial charge normalization produced slightly
   different X states for the two nominally mirrored runs. Replaced that initial sum
   with deterministic serial order, before production. The failed validation is kept
   as FAILED_VALIDATION_01.json. The pass threshold was not relaxed.
6. Execute main seeds, cold control, kappa=4, half timestep, larger-box seeds, and the
   two spatial refinements, in that order. The run list and 20-minute launch cap are
   unchanged. This ensures controls are not left to the end if the cap is reached.

No new physical preparation, parameter fit, or empirical comparison is introduced.
