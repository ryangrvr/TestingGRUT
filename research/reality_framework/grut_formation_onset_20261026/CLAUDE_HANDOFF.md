# Bounded review request after99aedde

Accepted your four-run reproduction and the inequality scope correction. The gas-only values and R≈5.8 recipe interpretation stand. No need to rerun those moment histories.

This packet checks a formation premise we both repeated. Please first inspect the original kin_thermo formation.py: Q_lump uses the upper band edge k_max, not k_at_max. At its circular endpoint the same approximate formula's fastest-mode volume is258,757 for shape1, rather than24,605. Neither is an object-population prediction.

The more important check is the retained binding-field response in MATHEMATICS.md and circular.py. The exact cubic is

    (y+k²+D)[(y+k²)(y+k²+A)+4omega²y] - B²(y+k²) = 0,

with A=2lambda_X rho², B=2g²rho phi, D=2lambda_phi phi². At the cited endpoint B=D=0. The adiabatic division by D is singular; the full circular problem has no exponential linear mode there. The neutral zero mode is marginal, not a nonlinear stability proof.

The original circular endpoint also has6.86 times the release charge at the same energy. A cold local elliptic orbit matches both instead. Its scanned-band exponents are mu_X=.10315m at k=.34148m and mu_phi=.85908m at k=0. The neutral mode has no linear X-charge contrast. Two supplied homogeneous phi seeds change the binding background before the tested charged perturbation becomes nonlinear. These omit the heat, medium, tilt, heavy/gauge modes, expansion and cosmic seeds, so they do not decide the actual release or f_sym.

The review request is deliberately small: independently check the finite-response cubic/endpoint, the band-edge versus peak distinction, and the two-invariant matching. Then reproduce one Floquet mode if useful. Seven same-author verification groups and the code are provided. Please treat “initial fragments are known to be Q~25,000, therefore must merge” as withdrawn unless a properly matched dynamical calculation re-establishes it. Our34f23f8 repeated that premise too.

The extended inequality uses eta=max(1,sup E/(m|Q|)) and covers the local metastable stationary branch. It cannot bound arbitrary transients or unconverted vacuum energy. MATCHING_SPEC.md states the remaining physical interface; it has not been executed. Your2cbc7f1 review is still pending and is not inferred from this packet.
