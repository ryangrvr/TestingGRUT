# Actual-source closed-channel kernel: classification and a real-energy bound

Analytic addendum to e9773891e0dc5986331e99b9ed7aa1a0d8ab95cb.
Author-derived; independent mathematical review pending. No Hamiltonian,
preparation, detector, parameter or time window is changed. This is not
a new numerical campaign or an operating pass. No source dynamics ran.

## Exact operator class

Use B2-AS's h on C f plus ell2(N_0 x N_0), with jp=1/64,jm=1/96,
A=13/46080,J=1/3840,Omega=2-2A. The exact min-coordinate representation
is a constant INFINITE-transverse-width block Jacobi operator plus the
fuel/root defect. On a finite physical charge quadrant 0<=i,j<L, layer
n=min(i,j) has 2(L-n)-1 relative states. A fixed-width compression is not
the finite physical host. The full source is unbounded, so it is not
unitarily a bounded free scalar Jacobi operator plus a bounded finite-rank
perturbation. A scalar low-energy cyclic representation would additionally
require its coefficients and transformed physical detectors to be derived.

## Closed channels can be eliminated on the real axis

Let P span f and every co-pair |n,n>, and Q=1-P. Set B=QhP. On Q write

    QhQ=D_Q+V_Q,
    (D_Q)_(i,j)=2+|i-j|+(1/4)1_(|i-j|=1), i!=j.

Then D_Q>=13/4 and, from the full primitive hopping norm,

    ||B||,||V_Q|| <=2(jp+jm)=5/96.

For Re z<=5/2 let C(z)=D_Q-z. It satisfies ||C^-1||<=4/3 and

    ||C^-1 V_Q|| <=q=5/72<1.

Thus the EXACT open-sector resolvent is

    P(z-h)^-1P=[z-PhP-Sigma(z)]^-1,
    Sigma(z)=B*(z-QhQ)^-1 B
      =-B* sum_(k>=0) (-C^-1 V_Q)^k C^-1 B.        (1)

This converges in operator norm uniformly throughout that closed
half-plane, including real energies. It retains every virtual path.
It does not apply a 1/Im z bound to the open propagation sector.

Each primitive charge hop flips i+j parity, while D_Q preserves parity.
B maps the even co-pair sector into the odd Q sector; Bf=0. Therefore
all odd-k terms in (1) vanish as operators. Let Sigma^[m] retain k=0,2,...,2m.
Summing the omitted even geometric series gives

    ||Sigma-Sigma^[m]||
      <= (5/96)^2 (4/3) q^(2m+2)/(1-q^2)
       = (75/20636)(5/72)^(2m+2).                 (2)

This bound is uniform in spatial size and Re z<=2.5; it does not grow
with the propagation time. It controls the closed-channel self-energy,
NOT the boundary value of the full open resolvent after inversion.

## Spatial structure is retained

A k-term path in (1) contains k+2 constituent hops. Transferring a
co-pair by d sites needs at least 2d such hops. Consequently Sigma^[m]
has co-pair bandwidth at most m+1. Each retained entry is an explicit
finite sum of path amplitudes, with its energy-dependent denominators.
Away from a finite root boundary layer these coefficients are translation
invariant; boundary paths are retained, not replaced by bulk paths.
The full Sigma is generally infinite range. For real E<=2.5 all relevant
path weights in -Sigma(E) are nonnegative. Two explicit paths from
|n,n> to |n+d,n+d> alternately use relative separation 1 and 2, remaining
in Q until the last hop. For every n>=0 and integer d>=1 they imply

    -Sigma_(n+d,n)(E)
      >= 2(jp*jm)^d /
         [(13/4-E)^d (4-E)^(d-1)] >0.             (3)

The first hop is plus (or minus), then d-1 repetitions of the leading
charge hop followed by the trailing charge hop, then the trailing charge
closes the pair. The two orientations are distinct and have identical
denominator products. There can be no cancellation of this lower bound
at these real energies, since every other Neumann path has the same sign.

Hence the eliminated operator is not exactly tridiagonal or a finite
boundary alteration of a nearest-neighbor scalar chain in the physical
co-pair basis. The reduction is energy dependent with controlled spatial
tails. This closes the theorem-applicability/transverse-truncation part
of the actual-kernel analysis, not the positive receiver-capture estimate.

The k=0 term is explicit:

    Sigma^[0]_(n,n)(z)
       =-(2-delta_(n,0))(jp^2+jm^2)/(13/4-z),
    Sigma^[0]_(n+1,n)(z)=-2jp*jm/(13/4-z),
    Sigma^[0]_(f,x)=Sigma^[0]_(x,f)=0.              (4)

At the algebraic reference point z=2, PhP+Sigma^[0](2) is exactly the
uniform fuel/co-pair half-line chain with diagonal Omega and hopping -J.
This recovers the old leading reference, but is not an exact replacement:
the full energy dependence and the remainder (2) remain in the source.
In particular z=2 is only a reference point, not a shifted fuel energy
or a retuned operating band.

## What this does and does not close

The source's unbounded closed-channel part now has a real-energy,
distance-sensitive, quantitatively controlled path expansion. It can be
reduced to a finite-band energy-dependent kernel with the explicit
remainder (2), while keeping the retained spatial phases.

The remaining task is to invert z-PhP-Sigma(z) near its open band and
control physical receiver kernels uniformly for n~L,t~L/J. Small Sigma
error alone does not bound this inverse near spectrum. One cannot
exponentiate the energy-dependent Schur operator as a time-independent
Hamiltonian. Threshold/pole analysis, source preparation, coherent
inversion and the existing finite-host comparison are still required.
No new late-return or arrival lower bound, actual violating point, or
completed B2 operating theorem is claimed.

The new calculation was analytic; elementary Fraction checks verify its
constants. Those checks are not a preregistered physical experiment or
a spectral classification, and do not replace independent proof review.
