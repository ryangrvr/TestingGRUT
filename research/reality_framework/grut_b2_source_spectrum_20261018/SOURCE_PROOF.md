# Source spectrum, preparation weights, and a finite retention interval

Author-derived; independent proof review pending. Fixed law, original
one-source fuel preparation, all charge separations and reverse processes.
The frozen joint comparison is used unchanged. Neither full source gate
is settled by these partial results.

## 1. Exact source spectral objects

On fuel f and all |i,j>, i,j>=0, the spatial source matrix is

    h_f,f=Omega=2-2A, h_f,(0,0)=-J,
    h_(i,j),(i,j)=2+|i-j|+(1/4)1_(|i-j|=1)-A 1_(i=j=0),
    plus hopping=-j+, minus hopping=-j-.

Here j+=1/64, j-=1/96, A=13/46080, J=1/3840. It is the self-adjoint
diagonal multiplication operator plus bounded hopping/conversion; finite
support is a form core. The scalar fuel spectral measure

    nu(B)=||1_B(h)f||^2

has total mass 1, mean Omega and variance J^2. These follow immediately
from h f=Omega f-J|0,0>.

Write its exact Lebesgue decomposition as

    nu=sum_b w_b delta_(E_b)+nu_sc+rho(E)dE,
    G(t)=sum_b w_b exp(-itE_b)
                   +integral exp(-itE)dnu_sc
                   +integral exp(-itE)rho(E)dE.       (5)

The decomposition is exact, but its point/singular-continuous weights
inside and above the first band are not solved here. A finite matrix's
discrete spectrum cannot identify them. The receiver probability A_L also
contains cross terms between spectral pieces; it is not the sum of their
separate probabilities at finite times.

Block inversion gives, without a memoryless approximation,

    <f,(z-h)^-1 f>=1/[z-Omega-J^2 m(z)],
    m(z)=<0,0|(z-H_ch)^-1|0,0>.                       (6)

Every pole, threshold and continuous part remains in this expression.

## 2. The exact bulk ground band brackets the fuel energy

Set Delta=5/4, b0=j++j-=5/192, bp=j+-j-=1/192. For simultaneous-translation
momentum k, the exact relative-coordinate fiber is

    h(k)=diag(2+|r|+(1/4)1_(|r|=1))-b(k)(S+S*),
    b(k)=sqrt(j+^2+j-^2+2j+j- cos k), r in Z.

It has a simple strictly positive even ground state; write e0=e(0),
epi=e(pi). The confining diagonal gives compact fiber resolvent. Replacing
a ground vector by its absolute value cannot increase its energy; the
irreducible nearest-neighbor eigen-equation makes it strictly positive.
The prior all-order derivative proof gives e'(k)>0 for 0<k<pi, so the
lowest bulk band is [e0,epi]. No relative cutoff enters this definition.

A normalized Rayleigh trial supported on r=0,+1,-1 with coefficients
1,b0/Delta,b0/Delta gives

    e0 <= 2-Q0,
    Q0=(2 b0^2/Delta)/(1+2 b0^2/Delta^2),
    Omega-e0>=d0:=Q0-2A>3J/2.                       (7)

At k=pi compress away r=0. That complement is >=2+Delta-2bp, while
epi<=2. Schur elimination of the complement and ||B||^2=2bp^2 imply

    epi>=2-2bp^2/(Delta-2bp),
    epi-Omega>=dp:=2A-2bp^2/(Delta-2bp)>3J/2.         (8)

Thus the ACTUAL fuel energy lies inside the exact mobile band with explicit
margins. This does not prove its spectral weight is all radiative.

## 3. No source bound state at or below e0

Let u_r>0 be the even, normalized bulk k=0 ground vector and d=Omega-e0.
Define a positive generalized vector on the ENTIRE source graph by

    w_(i,j)=u_(j-i),    w_f=J u_0/d.

It is not square-summable, because it is translation-constant along the
diagonal. Away from boundaries, (h-e0)w=0. Missing hops at the half-line
axes give

    ((h-e0)w)_(0,j)=j+ u_(j+1)>0, j>0,
    ((h-e0)w)_(i,0)=j- u_(i+1)>0, i>0.

The fuel residual is zero. At the root co-pair, the bulk equation
(2-e0)u0=2b0 u1 gives

    ((h-e0)w)_(0,0)
      =[b0 u1/u0-A-J^2/d]u0
      =[d/2-J^2/d]u0 >=(J/12)u0>0,                  (9)

using d>3J/2. Both source coupling directions and the degree counterterm
are included; removing conversion would be a different calculation.

For finitely supported phi, the exact ground-state representation is

    <phi,(h-e0)phi>
      =(1/2)sum_(x,y) t_xy w_x w_y
                           |phi_x/w_x-phi_y/w_y|^2
       +sum_x [((h-e0)w)_x/w_x]|phi_x|^2 >=0,         (10)

where h_xy=-t_xy for x!=y and the first sum counts both orientations.
Form-core closure gives h>=e0. If a square-summable vector were an
eigenvector at e0, all nonnegative terms vanish. The strictly positive
root residual forces its root component to vanish, and connectedness
forces every component to vanish. Hence

    spectrum(h) contains nothing below e0,
    point_spectrum(h) intersects (-infinity,e0] trivially. (11)

Bulk wave packets translated far from the source, with increasing center
support and a converging relative tail cutoff, are Weyl sequences at e0.
Thus e0 is the source spectrum's bottom, not an isolated source eigenvalue.
The same argument embeds every bulk ground-band energy in the essential
spectrum. This does NOT rule out point or singular-continuous measure
inside [e0,epi], or localized modes ABOVE that band.

## 4. Actual fuel preparation: quantitative spectral weights

Because the mean is Omega and the variance is J^2, one-sided Cantelli
with the exact dp in (8) yields

    nu((epi,infinity))<=J^2/(J^2+dp^2)<201/1000.

Combining this with h>=e0 gives

    nu([e0,epi])>799/1000.                            (12)

This is interval spectral weight of the ORIGINAL fuel state, not a
claim of 79.9% outgoing probability or an absolutely continuous weight.

Define the central energy interval C=[Omega-3J/2,Omega+3J/2]. Equations
(7)-(8) put it strictly inside the ground band, and Chebyshev gives

    nu(C)>=1-J^2/(3J/2)^2=5/9.                       (13)

These energies correspond to useful BULK velocities. To check that statement
without a time-growing approximation, apply B2-SO's exact weighted graph
argument on the full-line one-resource bulk. Its A_P is zero, hence obeys
the same conservative constants. Translation covariance of the unique
graph fixed point makes its co-pair generator a scalar Fourier symbol;
the exact correction bound ||E||_mu<J/12 implies

    |e(k)-(Omega-2J cos k)|<J/12.

For e(k) in C, |cos k|<19/24. The prior exact derivative lower bound
e'(k)>(16928/8825)J sin k then gives, for the outward k in (0,pi),

    e'(k)>(16928/8825)J sqrt(215)/24 > (11/10)J.       (14)

The exact bulk upper speed is <(11/5)J. Thus C lies in the velocity range
which fits the frozen receiver rectangle. Still, (13) does not identify
nu(C) as outgoing: embedded/singular source components and finite-time
coherent cross terms must be controlled before converting it into A_L.

The inherited high spectral bound remains

    nu([2+Delta-5/96,infinity))<=36/27613^2.

It does not pay for unclassified localized components inside the low band.

## 5. A continuous original-source return bound through theta=16

Use the existing B1 return-inclusive intertwiner bound, with theta=Jt.
It is for the SAME source h, with all virtual states and reverse conversion.
Its initial error e_init=1/115200 and residual rho<1/105000 give

    |G(t)-exp(-iOmega t)J_1(2theta)/theta|
       <=e_init+(rho/J)theta,
    |G(t)|<=1/theta+1/115200+(32/875)theta.           (15)

The full source is the uniform finite-path limit; B1's constants are length
independent and the finite-support source operators converge strongly on
the common form core. No full-geometry theorem is used. The scalar |J_1|
bound follows from its unit-modulus Fourier integral. The right side of
(15) is convex on theta>0, so its maximum on [2,16] is at an endpoint.
Exact rational comparison of both endpoints gives

    |G(t)|<13/20<2/3 for EVERY 2<=Jt<=16.             (16)

This is a continuous interval, not sampled favorable times. OBSERVABLES.md
explains why |G|^2 is the ENTIRE one-source fuel occupation uniformly over
input qubits. No conditioning on successful emission or arrival is used.

For the full physical two-source model, retain B2-SO's unchanged delta<0.001.
At L>=2016 and theta<=16, its one-source cut tail is <10^-6. With w<10^-6,
|g_low-G|<=w and the exact retained-matter identity,

    m_M>=1-10^-6-(13/20+10^-6)^2-10^-12.

The existing joint comparison therefore implies

    p_(q=2)^(L)(t) >=[m_M-1/1000]^2 >33/100,
                   2<=Jt<=16, L>=2016 in 24N.        (17)

q=2 counts two plus/minus pairs' worth of matter; constituents may be
separated. This is not two co-pairs in the distant receivers. The physical
I_L begins at theta=L/8>=252, so (16)-(17) cover only an initial segment
of the required retention interval and do NOT pass the operating gate.

## 6. Source verdict and exact open input

The source has no trapped point below or at its ground-band bottom; its
fuel energy lies inside that band; substantial initial spectral interval
weight is bounded; and the original fuel-return inequality holds through
theta=16. These are scoped partial source results, with proof review pending.

Unresolved: the actual point/singular/outgoing decomposition INSIDE and
above the ground band, its preparation weights, the outgoing generalized
eigenfunctions with physical co-pair detector, and their oscillatory
remainder uniformly for all t in I_L. No bound

    inf_(L,t in I_L) A_L(t)>=17/50

has been derived or refuted. Nor has |G|<=2/3 for theta from 16 through
L/6. Existing information gives no positive joint lower constant sufficient
for the target in that growing window. No explicit violating L or t is
claimed. The source sufficient route and the operating gate remain
mathematically unresolved; this cycle stops without changing the law.
