# Phase-sensitive open inverse and the exact remaining theorem

Author-derived; independent external proof review pending. Use (z-h)^-1.
All formulas are for the actual nominal source and original fuel seed.

## 1. Bulk Toeplitz inverse with its phases retained

Let kappa(z,k)=z-2-sigma(z,k), with sigma derived in BULK_BOUNDARY.md.
For Im z>0 and Re z<=5/2, the closed-channel resolvent has nonpositive
imaginary part, so Im kappa(z,k)>=Im z>0 on real k. Thus kappa is
nonzero on the circle with zero winding. Its principal logarithm has
an absolutely convergent Fourier series for each FIXED z:

    ell_j(z)=(1/2pi) integral_(-pi)^pi
                     log(kappa(z,k)) exp(-ijk) dk.

Indeed kappa is analytic in a k strip and is bounded away from zero on
the compact real circle at fixed z; it has a nonzero analytic neighborhood.
No width uniform as Im z approaches zero is claimed.

For w=exp(ik) define the scalar factors

    kappa_+(z,w)=exp[ell_0/2+sum_(j>=1)ell_j w^j],
    kappa_-(z,w)=exp[ell_0/2+sum_(j<=-1)ell_j w^j].  (1)

They and their inverses have absolutely summable one-sided coefficients,
and kappa=kappa_- kappa_+. With Toeplitz entries a_(n-m), coefficient
convolution gives EXACTLY

    A(z):=T_+(kappa)=T_+(kappa_-) T_+(kappa_+),
    g_T(z):=A(z)^-1=T_+(1/kappa_+) T_+(1/kappa_-). (2)

For the first identity, missing convolution indices l<0 vanish because
kappa_- has no positive coefficients. Each triangular Toeplitz factor
has its one-sided reciprocal inverse. This proves (2) without importing
a scalar finite-defect Jacobi theorem.

Let alpha_j=[w^j]1/kappa_+, j>=0, and beta_j=[w^j]1/kappa_-, j<=0.
Then the SPATIALLY RESOLVED bulk-compressed open inverse is

    (g_T)_(n,m)=sum_(l=0)^min(n,m) alpha_(n-l) beta_(l-m),
    (g_T)_(n,0)=alpha_n exp(-ell_0/2).              (3)

There is no complex conjugation in (3); this is an analytic resolvent
factorization, not a positive spectral decomposition. The actual symbol,
its zeros and its Fourier coefficients carry the propagation phase.
A norm bound is not substituted for these matrix elements.

## 2. The actual boundary and fuel coupling

Let U=W-A|0><0| be the derived trace-class boundary kernel. The charged
co-pair inverse and full fuel inverse are exactly

    S(z)=[A(z)-U(z)]^-1
        =[I-g_T(z)U(z)]^-1 g_T(z),
    F(z)=[z-Omega-J^2<0,S(z)0>]^-1,
    F_n(z)=-J F(z)<n,S(z)0>.                      (4)

This is equivalent to B2-AS's block-Weyl formula. It retains the
original fuel; no favorable outgoing state or fitted emission weight
has been supplied. For Im z>0 all inverses are the actual Schur
resolvents. The boundary product is trace class there, so its Fredholm
determinant can diagnose failure of inversion where the stated operator
realizations exist. An unweighted L2 real-axis determinant on the band
is NOT silently defined by (4).

Equations (1)-(4), together with the explicit bulk continued fraction,
are an actual spatially resolved open-channel representation. They do
not evaluate its boundary values or settle the operating gate.

## 3. What must be proved to pass to the physical amplitude

Let nu_f and kappa_n be the actual fuel and fuel-to-co-pair spectral
measures. At regular boundary values their densities come from
-Im F/pi and -Im F_n/pi. The entire amplitudes also contain all atom
residues and singular continuous measures:

    G(t)=integral exp(-itE)dnu_f(E),
    v_n(t)=integral exp(-itE)d kappa_n(E),
    A_L(t)=sum_(n=L/8)^(3L/8-1)|v_n(t)|^2.         (5)

The exact unsolved theorem is quantitative real-axis control of (3)-(4)
and their coherent Fourier inversion, including:

* the roots of kappa(E,k) and the behavior of its factors at e0,epi;
* the localized inverse I-g_T(E+i0)U(E), with source-visible exceptional
  poles or singular components either included or bounded;
* the fuel denominator in (4), its preparation residues and any resonance;
* the n~L, t~L/J oscillatory remainder, uniformly over the frozen receiver
  interval and continuous observation window.

The target is A_L>=17/50 throughout I_L and |G|<=2/3 on the complete
[2/J,L/(6J)] interval, or a sharper direct route to the physical event.
A possible violation of either SUFFICIENT source inequality is not itself
a certified physical failure. The B2-AS review packet supplies physical
upper event bounds needed for that conclusion.

## 4. Why the new small boundary norm does not finish the theorem

The derived ||U(E)||_1<J/500 on I_b does not imply
||g_T(E+i0)U(E)||<1 as an unweighted L2 estimate. g_T has propagating
spectrum there. Weighted/localized kernel estimates, or equivalent
explicit boundary scattering estimates, are still required. In particular
trace class alone does not exclude embedded or singular source components.
No general continuous-spectrum theorem has been imported without checking
its hypotheses for this energy-dependent operator pencil.

For a truncated self-energy, the user's inverse-error condition remains
essential: if R_m=K_m^-1, ||R_m||=M and the complete kernel error is
epsilon with M epsilon<1, then

    ||K^-1-R_m||<=M^2 epsilon/(1-M epsilon).

This follows by the resolvent identity, but is not uniform on the open
band merely because epsilon is small. No discarded Hamiltonian is
exponentiated as exact dynamics and no Im z resolution is suppressed.

## 5. Finite-host composition is retained, not re-inferred

The whole-line operator defines only the bulk symbol. W restores the
physical half-line source boundary exactly. The full finite TWO-source
H_L is then related to the original prepared half-line packets by the
existing B2-SO same-space defect estimate, not by equality of spectra:

    delta_L=12/27613+1500 exp(-L/105)<1/1000,
           0<=t<=L/(6J), L in 24N,L>=2016.

That estimate includes finite-host boundary crossings, collisions,
reabsorption and coherent returns on its stated horizon, pending proof
review. No additional finite-boundary approximation is introduced here.
If the missing source inequalities were earned, the original-trial
joint lower bound remains >458329/4000000=0.11458225, not 34%.
No such probability is established by the current kernel decomposition.

## Stop

The actual bulk-symbol and localized-boundary representation are derived,
with explicit complete remainders. The next step requires the named
real-axis factor/scattering theorem. This packet stops there under the
user's instruction, with the operating gate unresolved. It does not
replace that theorem by scalar moments, tiny errors, or another universe.
