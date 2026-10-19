# Calculations deciding the proposed closure fork

These calculations test particular inferences in the review of `8fdd674`.
They do not define a new cosmological quantum state.

## Classical Candidate 1 channels

The original dimensionless scalar potential is

\[
\bar U(F,S)=\tfrac14(S^2-1)^2+\tfrac12S^2F^2+.025F^4.
\]

For each prescribed integer charge N solve the same regular radial equations,
using an auxiliary integral J and Omega as a parameter:

\[
F''+2F'/r=(S^2+.1F^2-\Omega^2)F,\quad
S''+2S'/r=(S^2-1+F^2)S,\quad J'=r^2F^2.
\]

Boundary conditions are `F'(0)=S'(0)=J(0)=0`, `F(L)=0`, `S(L)=1`,
`4 pi Omega J(L)=g^2 N`. The singular origin terms are handled by the BVP
solver's regular singular matrix, not by discarding the origin. The energy is

\[
\frac{E_0(N)}{gv}=\frac{4\pi}{g^2}\int dr\,r^2
\left[\tfrac12(F'^2+S'^2)+\tfrac12\Omega^2F^2+\bar U(F,S)\right],\quad g=.1.
\]

At tree order `mX/(gv)=1`. The results at L50/tolerance 1e-10 are

| N | Omega | E0/(gv) |
|---:|---:|---:|
| 39944 | .799999626333 | 36392.42286136366 |
| 39943 | .800002069464 | 36391.62286051577 |
| 19972 | .885074195209 | 19729.91009896411 |

Thus `E0(N)-E0(N-1)-mX=-.199999152115 gv`, and
`E0(N)-2E0(N/2)=-3067.397336565 gv`. All six L40/L50 solves converge;
energy changes under refinement are below `1.4e-9 gv`. These comparisons use
the tree particle mass throughout; no tree/pole-order mixture is called a
quantum channel calculation.

The exact classical branch identity is `dE0/dQ=omega`. A finite integer step
is its integral `E0(N)-E0(N-1)=integral_(N-1)^N omega(q)dq`, not exactly the
endpoint value. Here the finite-step result `.800000847885 gv` differs from
`omega(N)` by about `1.22e-6 gv`, negligible for the `.2 gv` emission gap.

## The first slow-roll statement, with its actual scope

Use a free, massless, conformally coupled scalar in the conformal vacuum of
spatially flat FLRW. For a comoving worldline and noncoincident real proper
times, the vacuum two-point function is

\[
W(t_1,t_2)=-\frac{1}{4\pi^2a(t_1)a(t_2)[\eta(t_1)-\eta(t_2)-i0]^2}.
\]

Set the mean time to zero and normalize H(0)=a(0)=1. Let
`epsilon=-Hdot/H^2`, with the remaining slow-roll hierarchy appropriately
small. Write `s=t1-t2`, `z=s/2`. To first order at fixed separation,

\[
a(t)=\exp(t-\tfrac12\epsilon t^2)+O(\epsilon^2),
\quad a(s/2)a(-s/2)=1-\epsilon z^2+O(\epsilon^2),
\]

\[
\Delta\eta=2\sinh z+\epsilon[(z^2+2)\sinh z-2z\cosh z]+O(\epsilon^2).
\]

Consequently

\[
a_+a_-(\Delta\eta)^2=4\sinh^2z[1+2\epsilon(1-z\coth z)]+O(\epsilon^2),
\]

\[
W(0,s)=W_{\rm dS}(s;H=1)[1+2\epsilon(z\coth z-1)]
=W_{\rm dS}(s;H_{\rm eff}=1-\epsilon)+O(\epsilon^2),
\]

where `W_dS(s;h)=-h^2/[16 pi^2 sinh^2(h(s-i0)/2)]`. Restoring H gives

\[
H_{\rm eff}=H+\dot H/H=H(1-\epsilon)=\ddot a/\dot a.
\]

For constant epsilon the exact test uses
`a(t)=(1+epsilon t)^(1/epsilon)` and
`eta(t)=-(1+epsilon t)^((epsilon-1)/epsilon)/(1-epsilon)`.
The exact finite-separation residual after the temperature shift scales as
epsilon squared and approaches the independently derived coefficient

\[
\lim_{\epsilon\to0}\frac{W_\epsilon/W_{\rm dS}(H_{\rm eff})-1}{\epsilon^2}
=\frac{s^2}{4}-\frac{s^3}{12}\coth(s/2).
\]

This rechecks the conformal part of GRUT III v11. It is a local-in-mean-time
finite-separation statement; a Fourier/Wigner interpretation additionally
requires the frequency and adiabatic qualifications in the archived result.
An integral over arbitrarily long detector times is not justified by a
finite-separation expansion alone.

## Why the detector qualification matters at the same order

An endpoint-anchored detector uses `W(t,t-s)` rather than `W(t+s/2,t-s/2)`.
At t=0 its midpoint lies at `-s/2`; for power-law slow roll,
`H(-s/2)=1+epsilon s/2+O(epsilon^2)`. Therefore

\[
\frac{W_\epsilon(0,-s)}{W_{\rm dS}(s;1-\epsilon)}-1
=\epsilon\frac{s}{2}\left[2-s\coth(s/2)\right]+O(\epsilon^2).
\]

The coefficient is nonzero and depends on separation. The code verifies it
against the exact kernel. It also verifies first-order change under translating
both times. Thus the full state is nonstationary in comoving proper time;
a thermal form at a fixed midpoint does not establish a global KMS state or
an exactly thermal, finite-duration detector experiment. This calculation
does not recompute an entire detector response or v14's minimally coupled IR
sector. The latter's stored claim explicitly excludes its zero-frequency
polynomial contributions and specifies the Wigner convention.

Obadia's [2008 paper](https://arxiv.org/pdf/0804.2890), eqs. 20–22, gives the
conformal comoving local coefficient `Tloc=sqrt(H^2+2 Hdot)/(2 pi)`. Its
slow-roll expansion is the same `H(1-epsilon)/(2 pi)` at first order. The
paper distinguishes that coefficient from actual thermality. Our stronger
finite-separation recheck is not declared novel on the strength of a search
that found no identical sentence.

## Closedness and a finite curvature measurement

For k=+1 FLRW, `Omega_K=-1/(a^2 H^2)` in units c=1. Positive spatial curvature
fixes the sign; it supplies no nonzero lower bound on its present magnitude.
The ordinary no-boundary de Sitter saddle is an explicit limiting example:

\[
a(t)=H_\Lambda^{-1}\cosh(H_\Lambda t),\quad
H(t)=H_\Lambda\tanh(H_\Lambda t),\quad
N=\log[a(t)/a(0)].
\]

\[
\Omega_K=-\frac1{\sinh^2(H_\Lambda t)}=-\frac1{e^{2N}-1}.
\]

At N5 it is `-4.5402e-5`, and at N10 it is `-2.06115e-9`, while spatial
sections remain three-spheres. The cap is a regular Euclidean hemisphere,
continued from `a_E=sin(H_Lambda tau)/H_Lambda`. This exact background is
not a fit to our hot-big-bang history; it disproves the proposed implication
that a finite near-flatness bound alone excludes all closed no-boundary
histories. More generally, during inflation

\[
\frac{|\Omega_{K,f}|}{|\Omega_{K,i}|}
=e^{-2N}\left(\frac{H_i}{H_f}\right)^2.
\]

A reheating model and a distribution/bound on the duration of inflation are
needed to turn a curvature sign into a present numerical prediction. A model
can be statistically disfavored by its curvature distribution even when
arbitrarily small negative values remain in its support. That is different
from ruling out the entire sign class with a finite error bar.

Compactness alone is weaker still: periodically identify x,y,z in the flat
metric `dell^2=a(t)^2(dx^2+dy^2+dz^2)`. The spatial volume is finite and its
local scalar curvature is zero. This topological example is not asserted to
satisfy every no-boundary or global sequestering equation.

## A closed-geometry check on the proposed temperature transfer

There is an additional matching obligation when joining a closed no-boundary
saddle to the spatially flat v11 calculation. In the same global de Sitter
example, a comoving worldline is geodesic. For a conformal scalar in the
Euclidean state its exact geodesic two-point function has the de Sitter form
at `T=H_Lambda/(2 pi)`, independent of its midpoint. This uses the standard
Gibbons-Hawking result as an input; it is not a new mode-sum calculation here.

But the proposed scale-factor expression on this closed slicing is

\[
\frac{\ddot a}{\dot a}=H_\Lambda\coth(H_\Lambda t),\qquad
\epsilon_H=-\frac{\dot H}{H^2}=\Omega_K.
\]

The ratio to the actual geodesic temperature is
`sqrt(1+|Omega_K|)`, with discrepancy `|Omega_K|/2+O(Omega_K^2)`.
Thus extending the flat formula to a closed saddle without accounting for
curvature already fails in this exactly soluble limit. A valid slow-roll
bridge must compute the curved-space correlator or impose a matching limit
that makes curvature errors smaller than the retained slow-roll order. This
does not contradict v11's explicitly spatially flat result.

## What additional sequestering assumptions can do

An original global sequestering equation
`MP^2 Gmunu=Tmunu-(1/4)gmunu <T>` implies `<R>=0`. In homogeneous FLRW,

\[
\int dt\,a^3 R=6[a^3H]_{i}^{f}
+6\int dt\,(ka-a^3H^2).
\]

With finite integrals and a vanishing endpoint term, this would require
`k integral a dt=integral a^3 H^2 dt>0` for a nonstatic history. That is a
conditional dynamical result, with explicit boundary assumptions. It does
not follow from the phrase "nothing outside". Stiff singular endpoints and
their cutoff need separate treatment. The published local four-form version
of sequestering changes the residual constraint and does not require finite
spacetime volume; the original global restriction cannot be assigned to all
sequestering theories. Neither derivation supplies a minimum observable
|Omega_K| for a specified present epoch without a complete cosmological model.
