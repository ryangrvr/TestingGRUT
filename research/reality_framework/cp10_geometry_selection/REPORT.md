# CP10 reference feedback on ground geometry

2026-10-11. Frozen CP10 unchanged. Full-ensemble ground-state diagnostic at N2,3,4; no new candidate, parameters, selected graphs or Stage3 card.

Independently assemble CP10 from the published CP09 sparse builder, comparability reference Laplacian and reference-position matter interaction. All labelled orders included. Compare ground geometry probabilities with CP09 at the same N. Ground energies CP10: .260465360127,.482868701245,.780965486648. N4 joint dimension14016. Ground vector residual <=7.6e-15. Smallest computed gaps .21949772,.17686689,.09929163; no continuum extrapolation or complete degenerate-shell certification.

Geometry total variation CP10 vs CP09 is .00183661,.00415258,.00682705 respectively. N4 height probabilities:

|height|CP09|CP10|
|---|---:|---:|
|1|.00452363|.00474920|
|2|.41926436|.42190483|
|3|.47796601|.47580858|
|4|.09824601|.09753738|

This test distinguishes acquiring geometry-reference correlations from substantially selecting a new geometry distribution. The frozen reference coupling causes nonzero feedback, but only a small change at tested sizes. It neither shows a continuumlike class wins nor proves asymptotic failure. Height classes are not manifold classifiers. Low excited spectra are not operational propagation laws.

CP10 earns a jointly invariant positive dynamical reference; it does not resolve CP09's geometry selection. Additional reference accuracy or recording sophistication should not be credited with that task absent evidence. A next generating mechanism must physically explain selection and propagation together, with all new premises accounted for. No amended coupling or extra penalty has been inserted here.

Run python research/reality_framework/cp10_geometry_selection/check.py, NumPy/SciPy. Expected CP10_GROUND_GEOMETRY_CONTROL_PASS. Results preserve ground heights, full-configuration total variation and eigen residuals. Author computations; independent review outstanding; no new physics claimed.
