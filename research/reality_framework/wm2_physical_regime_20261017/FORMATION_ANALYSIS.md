# S9 formation theorem — complete histories and finite receiver

**Builder derivation, not independently reviewed.** The frozen preliminary
route is preserved at `85fc1bd`. No radial, subset-only, independent-emission,
Markov-bath or classical-history approximation enters the theorem below.
Small matrices and rational constants are checked in FORMATION_CONTROLS.json;
they do not constitute an all-N numerical test or external validation.

## Result and its unusually large finite-size boundary

For the unchanged WM2/S4 probe-vacuum law, every root, every independent
pattern delta_v=bar_delta+e_v with |bar_delta|<=1/2, |e_v|<=1/4,
g in [9/10,11/10], physical L=64N^2, and declared neutral/vacuum seed,

\[
 N\ge 2^{28},\quad 4N^2\le t\le8N^2
 \quad\Longrightarrow\quad
 \Pr[m\ge\lceil N/2\rceil]>\tfrac12.                 \tag{1}
\]

Consequently uniform time sampling on the **frozen** interval has probability
at least 1/2, above the frozen 1/4 target. This gives no useful formation
lower bound at five or six events, or uniformly over 6<N<2^28. It is a
conservative asymptotic theorem, not a tractable experiment or an efficient
simulation claim. It concerns a connected macroscopic partial path, not
full-N assembly, permanent locking, selected dimension or spacetime.

## 1. Exact word algebra and physical energy convention

Fix one root, M=N-1. S_v appends unused non-root event v to an ordered word;
S=sum_v S_v, K=word length excluding root, F=M-K. Unique parents imply

\[
 S_v^\dagger S_w=\delta_{vw}P_{v\,\mathrm{free}},\qquad
 S^\dagger S=F,\qquad \|S\|\le\sqrt M.              \tag{2}
\]

More generally A_z=sum_v exp(-i delta_v z)S_v satisfies A_z^dagger A_z=F.
Its norm is at most sqrt(M). In particular

\[
 \|[S,A_z]\|,\ \|[S^\dagger,A_z]\|\le2M.            \tag{3}
\]

**[S,A_z] is generally nonzero.** Removing it would silently replace the
last-in/first-out word grammar by a different emitter algebra.

In the consistent probe-vacuum formation component the exact relative law is

\[
 H=D+d\Gamma(h_L)-\lambda(Sa_0^\dagger+S^\dagger a_0),\quad
 D=\sum_v\delta_v P_{v\,\mathrm{attached}},\quad\lambda=g/N. \tag{4}
\]

h_L is the actual L-site path matrix with off-diagonal -1. Physical receiver
frequencies are 4+sigma(h_L) in [2,6]; (4) subtracts the conserved
4(M-K+n_R) and a constant. It changes no physical evolution in the prepared
charge sector n_R=K. To use creation/annihilation identities, extend (4) to
the ordinary full Fock space; all vectors used lie in finitely many adjacent
charge sectors. This introduces no additional bath, energy, or cutoff.
The decoupled CP-12 benchmark factor remains inert. All reverse terms stay.

## 2. A uniform correlation bound for the physical receiver

Let C_L(u)=<0|exp(-i h_L u)|0>. The half-line sine transform gives

\[
 C_\infty(u)=\int_{-2}^2 e^{-ixu}\rho(x)\,dx,\qquad
 \rho(x)=\frac{\sqrt{4-x^2}}{2\pi}.                 \tag{5}
\]

The density follows by x=-2cos(k) from the normalized boundary weight
(2/pi)sin(k)^2 dk, 0<k<pi; it integrates to one. For u>=1 split off the
two endpoint intervals of width 1/u. Their total absolute mass is at most
4/(3pi)u^(-3/2). On [-a,a], a=2-1/u, integrate twice by parts. Since
rho''<0 and |rho'(a)|<=sqrt(u)/(pi sqrt(3)), the interior absolute bound is

    [2/pi + 4/(pi sqrt(3))] u^(-3/2).

The total coefficient is <2 (already pi>3, sqrt(3)>5/3 suffice).
For u<1 use the spectral mass one. Thus

\[
 |C_\infty(u)|\le\min(1,2u^{-3/2}),\quad
 \int_0^\infty|C_\infty|\le5,\quad
 \int_0^\tau u|C_\infty(u)|du\le4\sqrt\tau\ (\tau\ge1),
\]
\[
 \int_\tau^\infty|C_\infty|du\le4/\sqrt\tau\ (\tau\ge1),
 \qquad \int_t^\infty|C_\infty|du\le5/\sqrt t\ (t>0). \tag{6}
\]

The finite and half-line closed walks agree through order 2L-1. Both
one-particle matrices have norm <=2. The Taylor tail therefore gives,
for 0<=u<=T_max=8N^2 and L=64N^2,

\[
 |C_L(u)-C_\infty(u)|\le
 2\sum_{n\ge2L}\frac{(2u)^n}{n!}
 \le\frac{16}{7}(e/8)^{128N^2}
 <\varepsilon_N:=2^{-120N^2}.                       \tag{7}
\]

Use (n!)>=(n/e)^n, term ratio <1/8, and e<3; the last deliberately loose
inequality holds for N>=1. In particular

\[
 \int_0^t|C_L(u)|du\le5+T_{\max}\varepsilon_N\le6. \tag{8}
\]

No finite receiver state is dropped: (7) compares its free correlation
inside an exact retarded identity, not a substituted many-body evolution.

## 3. Exact contact occupation and the history-retaining retardation bound

Let O(t)=exp(iHt)O exp(-iHt). Solving the linear receiver equation and
using the **initial receiver vacuum** gives the exact vector identity

\[
 a_0(t)\psi_0=i\lambda\int_0^t C_L(u)S(t-u)\psi_0\,du. \tag{9}
\]

By (2),(8), throughout the full observation horizon,

\[
 \|a_0(t)\psi_0\|\le b_N:=6\lambda\sqrt M,
 \qquad \langle n_0(t)\rangle\le36g^2M/N^2.          \tag{10}
\]

The field has not been reset; remote quanta and coherent returns are present.
Moreover ||a_0^dagger(t)psi_0||=sqrt(1+<n_0(t)>), which is **not** small.
For fixed u<=t differentiate [A_(u-v)](t-v), 0<=v<=u. The D commutator
cancels its explicit phase derivative, and dGamma(h_L) commutes with A.
The remaining derivative is -i[V,A_(u-v)](t-v). At equal times geometry
and receiver operators commute. Equations (3),(10) imply

\[
 \|(S(t-u)-A_u(t))\psi_0\|\le B_Nu,\qquad
 B_N=2\lambda M(\sqrt{1+b_N^2}+b_N).                \tag{11}
\]

For N>=64, b_N<=33/40 and B_N<187/40<5. Also the left side is at most
2sqrt(M). Split the memory integral at tau=sqrt(N). With (6), for every
t>0 (the same bound is loose but valid when t<tau),

\[
 \left\|\int_0^t C_\infty(u)(S(t-u)-A_u(t))\psi_0du\right\|
 \le4B_N\sqrt\tau+8\sqrt M/\sqrt\tau\le28N^{1/4}.   \tag{12}
\]

The creation commutator in (3) pays for the full history noncommutativity.
Neither independent emissions nor an irreversible channel was assumed.

## 4. Controlled current comparison, not an assumed master equation

The absolutely convergent half-transform obeys

\[
 2\operatorname{Re}\int_0^\infty C_\infty(u)e^{-i\delta u}du
 =2\pi\rho(-\delta)=\sqrt{4-\delta^2}.              \tag{13}
\]

For example obtain this by Abel damping eta>0, integrate the Lorentz
kernel against (5), then take eta down to zero; (6) permits dominated
convergence. Here every delta is strictly inside the band.

From (4), d<K>/dt=2lambda Im <S(t)psi_0,a_0(t)psi_0>. In (9) replace
only its retarded vector by the instantaneous phased vector, retaining an
explicit error. The three vector errors are at most

    lambda*28N^(1/4)              [history retardation],
    lambda*sqrt(M)*5/sqrt(t)      [initial missing tail],
    lambda*sqrt(M)*T_max*epsilon_N [finite receiver correlation].

Using S^dagger A_u=sum_v exp(-i delta_v u)P_(v free), one obtains the
**exact current with a bounded remainder**

\[
 \frac{d\langle K\rangle}{dt}
 =\lambda^2\sum_v\sqrt{4-\delta_v^2}
      \langle P_{v\,\mathrm{free}}(t)\rangle+r_N(t),
\]
\[
 |r_N(t)|\le56g^2N^{-5/4}+\frac{10g^2}{N\sqrt t}
                    +16g^2N\varepsilon_N,
 \quad N\ge64,\ 0<t\le8N^2.                       \tag{14}
\]

This is a finite-time first-moment comparison back to the full closed
Hamiltonian. It is **not** a trace-norm approximation by a Lindblad law,
a factorized state, or an independent-emitter probability distribution.
Its O(N^(-1/4)) macroscopic error survives the complete coherent history.

## 5. From the current to a macroscopic probability

Let f(t)=<F(t)>. Since |delta_v|<=3/4,
sqrt(4-delta_v^2)>=sqrt(55)/4>gamma_0=9/5. Equations (14) give

    f'(t)<=-gamma_0*(g/N)^2*f(t)+|r_N(t)|, f(0)=M.

Gronwall, retaining damping on the constant first remainder and bounding
only the other two integrals by their absolute values, yields

\[
 f(t)\le M e^{-\gamma_0(g/N)^2t}
     +\frac{280}{9}N^{3/4}+73+155N^3\varepsilon_N.  \tag{15}
\]

For 4N^2<=t<=8N^2, g>=9/10 gives gamma_0(g/N)^2t>=729/125.
An exact rational partial Taylor sum verifies exp(729/125)>256.
Therefore the first term is <M/256.

If macro formation fails, K<=ceil(N/2)-2, so F>=floor(N/2)+1>=N/2.
Markov's inequality applies to this **quantum observable distribution**;
it does not replace dynamics by a classical process. Thus

\[
 1-p_N(t)\le\frac{2f(t)}N
 <\frac1{128}+\frac{560}{9}N^{-1/4}
                   +\frac{146}N+310N^2\varepsilon_N. \tag{16}
\]

For N>=2^28 the right side is at most

    1/128 + 35/72 + 146/2^28 + 2^(-100) < 1/2.

The final receiver-tail inequality follows from monotonicity of
N^2*2^(-120N^2) for N>=1. This proves (1). No root-dependent constant
appears, so averaging the **original root weights** preserves it.
Time averaging on the originally frozen interval preserves it as well.

## Attempted counterexamples and limits

Arbitrary signed/ramped/repeated/distinct detunings inside the frozen class
are included by |delta|<=3/4; no permutation symmetry is used. The explicit
creation-commutator control refutes an absorption-only proof shortcut.
Returning projection errors refute an all-time square-sum numerical shortcut.
These are counterexamples to proposed reasoning, not counterexamples to
formation. No in-class physical counterexample was found or is claimed.

The estimate becomes vacuous at accessible sizes and does not settle the
unresolved six-event target. The resource sequence is huge: at N=2^28,
L=2^62 receiver modes, capacity N-1 quanta, and observation t>=2^58
in receiver-hopping units. Coupling to each event is g/N; asymptotic times
are paid for that normalization. A typical formed path can be used only
conditionally for an O(N) interior-probe window, versus O(N^2) formation.
This is not a permanently stable phase. See INTEGRATED_PROTOCOL.md.

The proof uses a specific supplied receiver bandwidth and vacuum plus a
single-end grammar. It establishes a model property, not a unique ontology,
a selected physical dimension, universal original charged matter, gravity,
or a distinctive empirical prediction. Independent mathematical review
of (7), (11)-(14) and the numerical sweep certificate remains desirable.
