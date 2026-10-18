# Derivative and preparation-tail verification

Author-derived; independent external mathematical review pending.
Use the actual h and the (z-h)^-1 convention pinned in SOURCE_PINS.json.
P spans fuel and all co-pairs, Q=1-P, B=QhP, C=QhQ.

The diagonal separated-charge energy is >=13/4 and its hopping norm
is <=5/96, hence C>=307/96. For real E<=5/2,

    C-E>=67/96,
    Sigma(E)=B*(E-C)^-1B <=0,
    ||Sigma(E)|| <=(5/96)^2/(67/96)=25/6432.

Norm-differentiating the bounded closed-sector resolvent gives

    Sigma'(E)=-B*(E-C)^-2B <=0,
    ||Sigma'(E)|| <=(5/96)^2/(67/96)^2=25/4489.

The fundamental theorem of calculus in operator norm proves the user's
Lipschitz estimate. For K(E)=E-PhP-Sigma(E),

    I<=K'(E)<= (1+25/4489)I.

These are valid closed-sector estimates. They do not bound K(E)^-1 on
its propagating spectrum. No continuous-spectrum theorem is inherited
from an isolated-eigenvalue perturbation result.

## The proposed high tail is valid but weaker than the existing one

Since (h-Omega)f=-J|0,0>, the fuel spectral measure has mean Omega and
variance J^2. Therefore, also including the cutoff point,

    nu_f([5/2,infinity)) <=J^2/(5/2-Omega)^2
                          =36/133010089,
    sqrt(tail)<=6/11533.

This checks the user's arithmetic independently. It is not used to replace
the stronger B2-SO value. Here that stronger separation can be rederived:
PhP <=(2+J)P. For 2+J<E<307/96, QhQ-E is positive and invertible, and
the Schur complement

    PhP-E-B*(QhQ-E)^-1B

is strictly negative and boundedly invertible. Thus the full h has no
spectrum in that interval. The cutoff 5/2 is inside it, so every high
spectral value is >=307/96. Applying the same variance estimate gives

    omega=nu_f([5/2,infinity))
          <=J^2/(307/96-Omega)^2=36/27613^2,
    sqrt(omega)<=6/27613.

The existing stronger bound and the existing joint error ledger are
preserved. Neither tail controls possible singular/point contributions
inside the low-energy sector. No physical filtering is performed.
