# What actually moves in P1: neutral bound-pair bands

This is analysis of **unchanged P1 during its controlled operating stage**,
with conversion off as P1 specifies. It is not an effective theory of A1
with conversion on, not a claim that P1 prepares a ground-band gas, and not
relativistic matter. Results are author-derived, independent review pending.

## 1. Exact one-pair Hamiltonian

On a fixed recorded path with one plus and one minus, basis |i,j> labels
their positions. Gauss determines one unit of field on the interval between
them and zero elsewhere. Capacity one discards no one-pair position states.
The exact Hamiltonian, keeping every separation and both orientations, is

    H|i,j>=[2+|i-j|+kappa 1_(|i-j|=1)]|i,j>
             -j_+ sum_(i'~i)|i',j>-j_- sum_(j'~j)|i,j'>.     (1)

Co-location is allowed. Thus constituents are confined; their low kinetic
energy does not describe independent free charges. On an infinite uniform
recorded path, Fourier-transform the center coordinate and remove a relative
phase. At center momentum k the relative-coordinate Jacobi operator is

    H_k=2+|r|+kappa 1_(|r|=1)-alpha(k)(T+T^*),
    alpha(k)^2=j_+^2+j_-^2+2j_+j_- cos k.                    (2)

k is a conditional effective-path variable, not primitive coordinates.
The diverging relative potential has discrete relative bound states for
each k. Their center dispersion determines composite motion.

Let tau=2(j_++j_-)=5/96, Delta_0=1+kappa=5/4. On a finite path the kinetic
norm is <=tau, uniformly in length; (2) has the same bound. With P_0 the
co-located subspace, Q_0=1-P_0, Q_0HQ_0>=2+Delta_0-tau. Exact Schur at
E=2 and the complete Neumann series, INCLUDING every separated state, give

    H_eff(2)=2P_0 - A sum_x deg(x)|x><x|
                    -t_0 sum_(x~x+1)(|x+1><x|+h.c.)+R_0,
    A=(j_+^2+j_-^2)/Delta_0=13/46080,
    t_0=2j_+j_-/Delta_0=1/3840.                              (3)

There are two virtual orderings of the charge hops; their amplitudes add.
No elementary neutral-pair hopping was inserted into P1.

## 2. Complete virtual error, not only second order

Each constituent hop flips parity of i+j. Co-located endpoints have even
parity, so all odd total hopping orders vanish. With r_0=tau/Delta_0=1/24,
the complete remainder after the second-order term obeys

    ||R_0|| <= (tau^2/Delta_0) r_0^2/(1-r_0^2)
             =tau^2/(575 Delta_0).                          (4)

The low band has E<=2 by min-max and E>=2-b_0,
b_0=tau^2/(Delta_0-tau). Resolvent identity bounds the additional
energy-dependence correction by tau^4/(Delta_0-tau)^3. Therefore the exact
lowest bound-pair dispersion satisfies

    |E_0(k)-[2-2A-2t_0 cos k]| <=d_0,
    d_0=tau^2/(575Delta_0)+tau^4/(Delta_0-tau)^3 <8.1e-6.     (5)

The leading bandwidth 4t_0=1/960 exceeds 2d_0. In particular the exact
band is not flat: neutral bound objects are mobile. This is an eigenvalue
and virtual-matrix result, not license to exponentiate H_eff(E) as an exact
time Hamiltonian. A long-time packet approximation needs a low-band isometry
and its own time-error estimate; that is not asserted here.

For the isolated low spectral projector Pi_0, a Sylvester/Laplace-integral
bound on Q_0 Pi_0 gives

    ||Q_0 Pi_0|| <=tau/(Delta_0-tau)=1/23.                   (6)

Indeed QHQ is above 2+Delta_0-tau and the low spectrum is <=2;
integrating the off-diagonal Sylvester equation costs this spectral separation.
Thus low-band states have small separated-charge weight, but an initially
adjacent pair (a vector in Q_0) has low-band weight at most 1/529. P1's
adjacent-pair production is NOT a ground-band preparation mechanism. This
statement is for one pair; no independent-pair approximation at density is made.

## 3. The initially relevant adjacent-dipole bands

P1's one-pair preparation instead starts near separation |r|=1, with energy
E_1=2+Delta_0=13/4 after its density interaction is on. The two orientations
are degenerate before hopping. Retain r=+1,-1; the r=0 intermediate lies
Delta_0 below, and r=+/-2 lies Delta_2=3/4 above. In the relative gauge of
(2), the complete second-order block is

    H_1^(2)(k)=E_1 I+alpha^2[(-8/15) I+(4/5) sigma_x].      (7)

The two parity branches therefore have

    E_even^(2)=E_1+(4/15)alpha^2,
    E_odd^(2) =E_1-(4/3)alpha^2.                             (8)

Their leading center hopping coefficients in the convention -2t cos k are
t_even=-1/23040 and t_odd=1/4608, different from t_0. Orientation is a
physical internal degree, measurable by the sign of the intervening field.
In an ungauged orientation basis the off-diagonal acquires the corresponding
k-dependent phase; no common independent spin qubit is implied.

For an explicit full-sector error let g_1=3/4, r_1=tau/g_1=5/72,
c_1=g_1-tau, b_1=tau^2/(g_1-2tau). The complement has spectrum separated
from E_1 by at least c_1. Initial Weyl bounds |E-E_1|<=tau followed by the
Schur equation improve this to |E-E_1|<=b_1. Parity again eliminates odd
hopping orders. The uniform error for either branch in (8) is at most

    d_1=(tau^2/g_1)r_1^2/(1-r_1^2)
               +tau^2 b_1/[c_1(c_1-b_1)] <4.2e-5.           (9)

The even leading bandwidth 1/5760 and odd bandwidth 1/1152 each exceed
2d_1. Both exact adjacent-dipole branches are nonflat. Thus the produced
neutral excitations can move in more than one internal band; their speeds
and mixing depend on that band and the supplied couplings. This is not a
universal relativistic cone or observed particle identification.

The finite relative-cutoff matrices in the code are diagnostics of (5),(9).
Proofs (4),(9) sum complete virtual sectors using norm convergence and do
not replace the infinite potential by a small receiver. No cutoff trajectory
is promoted to an all-size formation or scattering result.

## 4. Two pairs: an interaction, not a completed scattering theory

For two plus/minus pairs, let P_c contain configurations with two co-located
neutral pairs, with hard-core pair occupancy d_x. On that subspace,

    P_c H P_c=4P_c+2kappa sum_x d_x d_(x+1).                 (10)

Adjacent co-pairs cost an additional 1/2 compared with separated ones. All
Q_c configurations have at least one field unit, so their diagonal is >=5.
The full two-pair kinetic norm is <=tau_2=2tau=5/48. At E_*=9/2,

    ||Sigma_2(E_*)|| <=tau_2^2/(1/2-tau_2)=:B_2<0.028.       (11)

Consequently the exact Schur diagonal difference between an adjacent and a
separated co-pair configuration is at least 1/2-2B_2>0.444. This is a
nonzero short-range collision interaction in the neutral sector, not a free
gas of independent composite particles. It retains all charged virtual paths.

An on-shell scattering matrix, many-pair kinetic theory, and robust transport
of an encoded internal qubit are **unresolved**. Equations (7),(10),(11)
identify physical orientation and collision channels, not their complete
information lifetime or asymptotic scattering amplitudes. Reporting them as
an established interacting neutral phase would exceed this bounded result.

The primary A1 retention gate already fails; this cycle stops rather than
opening an unbounded scattering/continuum or replacement-architecture campaign.
