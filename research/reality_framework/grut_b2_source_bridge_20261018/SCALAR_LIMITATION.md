# A recurrent measure satisfying the currently banked scalar summaries

Author-derived mathematical countermodel to an implication. This is NOT
the spectral measure of the unchanged source h, not an added candidate,
and not a violating B2 length/time. It shows exactly what cannot be inferred
using only the named scalar summaries. The actual resolvent and primitive
Hamiltonian contain additional information not imposed on this countermodel.

## 1. One declared measure, no fitted time or grid

On the dimensionless energy x=(E-Omega)/J define

    dmu_*= (23/791) delta_0
             +(768/791)(1/225) sum_(n=-112)^112 delta_(n/64).  (1)

All weights are positive, with the zero atom including both contributions.
Normalization and symmetry give mass one and mean zero. Since

    (1/225) sum_(n=-112)^112 (n/64)^2
        =112*113/(3*64^2)=791/768,

its variance is EXACTLY one. Therefore its physical energy measure has
the same mean Omega and variance J^2 as the original fuel preparation.
All its support is in [-7/4,7/4]. B2-SS's proven bulk margins d0/J and
dp/J are both >7/4, so every atom is strictly inside the exact reference
ground-band energy interval. In particular there is no atom at/below e0,
the ground-band interval mass is one (>79.9%), and the mass above 2.5 is
zero. These statements concern compatibility with energy summaries, not
the existence of this measure for the original source.

The central interval [-3/2,3/2] includes exactly n=-96,...,96. Thus

    mu_*(central)=23/791+(768/791)(193/225)
                 =51133/59325 >5/9.                (2)

These energies are also in the velocity-compatible interval of the
reference bulk dispersion. No claim about the countermodel's own bulk
dispersion or spatial primitive moves is made.

## 2. Continuous early-time depletion

After removing the irrelevant phase exp(-iOmega t), its exact return is

    F_*(theta)=23/791+(768/791)
       sin(225 theta/128)/[225 sin(theta/128)],       (3)

with the continuous value at removable singularities. For theta in [2,16],
theta/128 lies in [1/64,1/8]. The sine denominator is positive, and

    sin(theta/128)>=sin(1/64)
          >=1/64-1/(6*64^3)>0.

Bounding the numerator by one yields the CONTINUOUS inequality

    |F_*| <=23/791+(768/791)/
                      [225(1/64-1/(6*64^3))]
           <31/100<13/20.                           (4)

No time grid is used. It satisfies the banked |G|<.65 early-window
consequence more strongly, along with the moment and energy facts above.
This does not assert that it satisfies every inequality in B1/B2-SO.

## 3. An exact return in a declared late window

At the predeclared theta_*=128*pi every phase exp(-i theta_* n/64)
equals one. Consequently

    |F_*(128*pi)|=1.                                (5)

Take the predeclared L=3072 in 24N, L>=2016. Its frozen interval in
dimensionless time is [384,512]. From 3<pi<4, theta_* lies strictly
inside that interval. Thus the listed scalar summaries, even together,
do NOT imply |G|<=2/3 on the late interval. They permit complete return.

A finite Hilbert realization is multiplication by Omega+Jx on
L2(mu_*), with normalized cyclic seed f_*=1. This has 225 distinct
eigenvalues. The identity

    h_* f_*=Omega f_*-J c_*,
    c_*=-(h_*-Omega)f_*/J

has c_* normalized and orthogonal to f_* by the exact moments, matching
the first fuel-column identity alone. At theta_* the complete state is a
phase times f_*. ANY chosen detector orthogonal to that seed has zero
arrival probability then. This demonstrates that scalar fuel facts do
not in themselves enforce detector arrival either.

## 4. Exact scope of the negative result

The invalid inference is:

    specified energy/moment weights + no bottom atom
       + continuous early fuel bound
       ==> uniform late fuel bound or receiver capture.

The witness refutes THIS implication. It does not refute source operation,
the comparison theorem, source-locality, a scattering proof, or all possible
analytical methods. It is wholly pure point inside the permitted energy
interval, precisely a spectral type not yet excluded by the banked facts.
It is not claimed to obey the original infinite lattice resolvent, exact
charge-coordinate support, full essential spectrum, or the small weighted
generator correction. Those extra structures may still exclude recurrence
for the actual h; using them to prove (15) or a late Fourier estimate is
the unresolved physical analysis.

This is the permitted 'demonstrated limitation of the analytical method'
outcome: no extrapolation from scalar summaries can close the growing
window. No actual violating source point is reported.
