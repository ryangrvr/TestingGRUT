# Exact remaining actual-source scattering theorem

Author-derived kernel results are available, but the following theorem
has not been proved. This is a prepared mathematical problem, not a claim
that an external reviewer was contacted or that mathematical work must end.

## Inputs already derived for the actual operator

Nominal jp=1/64,jm=1/96,kappa=1/4,A=13/46080,J=1/3840,Omega=2-2A.
The bulk closed-channel symbol is the actual continued fraction

    v(k)=jp^2+jm^2+2jp jm cos k,
    sigma(z,k)=2v(k)/[z-13/4-v(k)/[z-4-v(k)/[z-5-...]]].

Set kappa(z,k)=z-2-sigma(z,k), A(z)=T_+(kappa), g_T=A^-1. The
fixed-Im z>0 spatial kernel is explicitly given by the Fourier-log
Wiener-Hopf factors in SOURCE_KERNEL.md. The actual co-pair boundary
correction U=W-A|0><0| (here scalar A=13/46080 in the second expression)
is trace class, exponentially localized and satisfies

    ||U(E)||_1<J/500, E in I_b=[2-5J,2].

I_b contains the known first bulk band, not necessarily every source
atom below 2.5. The fuel is NOT absorbed into this small correction:

    S(z)=[I-g_T(z)U(z)]^-1 g_T(z),
    F(z)=[z-Omega-J^2<0,S(z)0>]^-1,
    F_n(z)=-J F(z)<n,S(z)0>.

All conventions are (z-h)^-1. The upper/lower jump, residues and any
singular measures reconstruct the actual source amplitudes.

## The theorem needed, rather than another norm estimate

1. Establish quantitative real-axis limits for the actual scalar factors
   and the localized inverse I-g_T(E+i0)U(E), including first-band edges
   and source-visible exceptional energies. An unweighted L2 inverse
   bound on the band is not assumed.
2. Bound the fuel denominator and every source-visible pole residue or
   singular contribution. A.e. regular density alone cannot exclude
   exceptional singular measure. Tiny weight above 2.5 cannot pay for
   unknown atoms inside the low-energy sector.
3. Derive source-to-site amplitudes with their outgoing phase, with
   controlled remainder uniform for n~L,t~L/J. A possible sufficient
   interior form to PROVE is <n,S(E+i0)0>=a_out(E)exp(ik(E)n)+r_n(E),
   with explicit edge terms and derivative/remainder constants. This
   form is not claimed as an established result.
4. Integrate coherently before summing the detector probabilities. Prove

       inf_(L in 24N,L>=2016;t in I_L)
          sum_(n=L/8)^(3L/8-1)|integral exp(-itE)d kappa_n(E)|^2 >=17/50,
       sup_(L in family;2/J<=t<=L/(6J))
          |integral exp(-itE)dnu_f(E)| <=2/3,

   or establish the same physical operating gates by a sharper direct
   route. The 96% central-capture sufficient condition is not mandatory.

The original finite-host comparison, including returns, remains
delta_L=12/27613+1500 exp(-L/105)<.001 on [0,L/(6J)]. If the above
sufficient source conditions hold, physical joint output exceeds
0.11458225 and q=2 retention exceeds .25, with exact conditional two-qubit
quality under the endpoint law. No factorization of physical successes
or incoherent mixture of energies is assumed.

To certify a physical failure instead, use controlled UPPER event bounds
for the original experiment. A failure of one sufficient source condition,
of the small-kernel Neumann route, or of a proposed asymptotic form alone
would not refute the operating gate.

## Method scope

The explicit actual kernel closes the bulk/source-boundary decomposition.
It does not make the remaining inverse well conditioned. The Feshbach
algebraic invertibility relation is distinct from a quantitative
continuous-spectrum/scattering theorem. For methodological context see
Griesemer and Hasler, [On the Smooth Feshbach-Schur Map](https://web.ma.utexas.edu/mp_arc/c/07/07-102.pdf).
No propagation theorem from that paper is claimed for this source.
