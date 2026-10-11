# CP09 clock error budget: discrete pulse-area correction

2026-10-11. Exact committed independent review 78013a6 was executed and passed. Its width .15 comparator at the SAME weak areas epsilon=zeta=.001 is .052349192601; the .052349242 figure uses smaller weak areas. Difference approximately 5e-8 is immaterial to the reported percent discrepancy but should be counted consistently.

## Numerical calibration issue and controlled rerun

Published bounded-clock simulation sampled square windows as indicator/width. At finite grid this does not have exactly unit discrete integral. Source-area / detector-area / product:

256: .9375 / 1.04166667 / .9765625
512: .98958333 / .98958333 / .97927517
1024: 1.015625 / .98958333 / 1.00504557
2048: 1.00260417 / 1.00260417 / 1.00521512

The leading weak source-pointer signal scales with the area product. Thus grid-dependent area is an additional calibration error, separate from intrinsic clock effects. This does not retroactively change the frozen simulator. A new explicit control normalizes each sampled window by its own number of points times dx and uses that same coefficient in evolution AND energy accounting. Finite-grid support widths/centers still differ from the continuum; normalization alone does not remove every discretization error.

Same mass320, width .15, preparation, source/detector and all other controls:

| grid | time step | unit-area response | residual vs instantaneous ideal |
|---|---:|---:|---:|
|1024|.0025|.053683215575|+.001094702036|
|2048|.00125|.053683095065|+.001094581526|

Responses differ by 1.21e-7 on this combined refinement, compared with 8.98e-6 before normalization. Energy drift <=1.2e-10. Residual ideal-response discrepancy is now 2.08%, rather than 2.61%. No claim of complete convergence from two grids.

At the finest grid, subtracting the matched finite-width translation benchmark gives .001333902464. Therefore the residual budget is:

finite-width translation shift = -.000239320937;
additional bounded-clock/numerical/interaction difference = +.001333902464;
net instantaneous-ideal difference = +.001094581526.

The second term is NOT uniquely dispersion or back-action; numerical differences and the finite pulse shapes remain mixed. Larger mass/width studies need the same pulse calibration and separate spatial/time refinements. The result confirms that narrowing windows alone cannot be credited with fixing the clock.

## Scope

This improves an explicitly supplied conventional apparatus test. It does not derive global geometry-role access, locality, a relational time origin, continuum geometry or gravity. CP09 and its prior results remain frozen. No candidate scan or CP10. Run python research/reality_framework/cp09_clock_area_control/check.py. It loads the original simulator definitions and applies an openly specified calibration control without editing that file. Results retained; author numerical evidence only.
