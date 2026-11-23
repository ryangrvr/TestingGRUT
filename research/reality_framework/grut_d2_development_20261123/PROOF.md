# GRUT-D2: exact reduction and finite-window bounds

23 November 2026. These are mathematical statements about the supplied finite
Hamiltonian and preparation in `MODEL.md`. They establish autonomous ordering,
a passive count record, and recovery after one specified intervention. They do
not establish an unsupplied target class, developmental feedback, a new GRUT
law, consciousness, or cosmological selection. No thermodynamic or Markovian
limit is used.

## 1. A local spin Hamiltonian is exactly a wall-and-reservoir path

Let the physical spins be indexed by `0,...,n-1`, put `m=n-1`, and define the
computational-basis change

\[
 r=s_0,\qquad u_i=s_i\mathbin\oplus s_{i+1}\quad(0\le i<m).
\]

It is a bijection: `s_j = r xor u_0 xor ... xor u_(j-1)`. Consequently it is a
unitary change of basis on the physical system, with no projection or extra
degrees of freedom. In these coordinates

\[
 \tau_i^z=Z_iZ_{i+1},\qquad
 \tau_i^x=\prod_{j=i+1}^{n-1}X_j,\qquad
 \tau_i^y=i\tau_i^x\tau_i^z,\qquad
 w_i=\frac{1-\tau_i^z}{2}.
\]

The suffix product flips exactly wall `i` and leaves the orientation `r`
unchanged. Although an individual wall Pauli operator has a string, the
Hamiltonian exchanges adjacent walls and its strings cancel:

\[
 \begin{aligned}
 \tau_i^x\tau_{i+1}^x&=X_{i+1},\\
 \tau_i^y\tau_{i+1}^y&=-X_{i+1}Z_iZ_{i+2},\\
 g(\tau_i^+\tau_{i+1}^-+\mathrm{h.c.})
 &=\frac g2X_{i+1}(1-Z_iZ_{i+2}).
 \end{aligned}
\]

At the exit, `tau_(m-1)^x=X_(n-1)` and
`tau_(m-1)^y=Z_(n-2)Y_(n-1)`. Hence

\[
 g(\tau_{m-1}^-b_0^\dagger+\mathrm{h.c.})
 =\frac g2\left(X_{n-1}X_{B0}+Z_{n-2}Y_{n-1}Y_{B0}\right).
\]

Thus all stated physical terms act on at most three neighboring spins. The
orientation `Z_0` is conserved. Global physical-spin reversal commutes with the
Hamiltonian, so this law favors alignment without preferring either absolute
orientation. A conserved boundary seed is nevertheless supplied.

Order the logical path as the `m` walls followed by the `L` reservoir sites.
The Jordan-Wigner transformation gives exactly

\[
 H=\Delta Q+\sum_{j,k=0}^{m+L-1}T_{jk}c_j^\dagger c_k,
 \qquad
 T_{jk}=g(\delta_{j,k+1}+\delta_{j+1,k}),\qquad
 Q=\sum_jc_j^\dagger c_j.
\]

The orientation is an independent inert two-state factor. `T` is the finite
open-path hopping matrix, with norm at most `2g`; therefore
`H >= (Delta-2g)Q >= 0` when `Delta>2g`. For the registered `Delta=4g`, its
vacuum family is exactly the two aligned physical configurations with an
empty reservoir. This energetic preference is an input to the construction.

## 2. Exact propagation and unconditional organization

Write `N=m+L`, `P` for the first `m` sites, and

\[
 U(t)=e^{-iTt},\qquad A(t)=PU(t)P\big|_{P\mathcal H}.
\]

The on-site term contributes a common one-particle phase `exp(-i Delta t)`,
which cancels from the number probabilities used below. The exact finite
path eigenmodes give

\[
 A_{ij}(t)=\frac2{N+1}\sum_{k=1}^{N}
 \sin\frac{ik\pi}{N+1}\sin\frac{jk\pi}{N+1}
 e^{-2igt\cos[k\pi/(N+1)]},\qquad 1\le i,j\le m.
\]

Initially the physical spins are maximally mixed, the orientation is fair,
the wall occupations are independent with `p=1/2`, and the reservoir is
vacuum. The wall-and-reservoir state is a number-conserving fermionic
Gaussian state with one-particle covariance `C(0)=pP`. Therefore

\[
 C_S(t)=pA(t)A(t)^\dagger,\qquad
 \mathbb E N_{\rm wall}(t)=p\|A(t)\|_F^2.
\]

Diagonalizing this restricted covariance turns its number distribution into
a sum of independent Bernoulli occupations. Its vacuum probability is

\[
 \boxed{P_{\rm order}(t)=\det[I-pA(t)A(t)^\dagger].}
\]

This is the probability that all physical spins agree, accepting both
orientations. It is unconditional. Initially it is `2^(-m)`. Since the
indicator of any remaining wall is bounded by the wall number,

\[
 \boxed{P_{\rm order}(t)\ge1-p\|A(t)\|_F^2.}
\]

The analytic bound, rather than a sampled trajectory, will certify an entire
finite interval. The numerical mode sum evaluates the same finite model; it
does not approximate a bath continuum or truncate a many-particle sector.

## 3. The reservoir writes a count record

The reference `R` is a supplied inert diagnostic copy of the initial wall
count `K`, with alphabet `0,...,m`. It is not newly generated memory. Resolving
this classical preparation into each value of `K`, the active state remains
in the sector

\[
 Q=N_{\rm wall}+N_B=K,
\]

because `[H,Q]=0`. Thus the reservoir decoder `K_guess=N_B` succeeds if and
only if no wall remains:

\[
 \boxed{P(K_{\rm guess}=K)=P(N_B=K)
       =P(N_{\rm wall}=0)=P_{\rm order}.}
\]

This identity involves no postselection and establishes a newly written
reservoir count record. It does not establish retention of the original
spatial wall pattern. For the specified decoder with error `delta`, Fano's
inequality also gives

\[
 I(K:N_B)\ge H(K)-h_2(\delta)-\delta\log_2m,
\]

where `K` has the initial binomial distribution and the last term is zero
for `m=1`. Initially the empty reservoir has zero count information.

The record is passive. The Hamiltonian contains no controller that reads
this count and chooses a subsequent action. The collective readout is an
access rule supplied by the modeler, and `N_B` is not exactly conserved
while exchange with the walls continues.

## 4. A uniform half-line estimate

Use a half-line only as an analytic comparison, not as the physical system.
Its hopping propagator is

\[
 U^{(\infty)}_{ij}(t)=\frac2\pi\int_0^\pi
 \sin(ik)\sin(jk)e^{-2igt\cos k}\,dk,
 \qquad i,j\ge1.
\]

For positive integers `i,j`, set

\[
 F_{ij}(k)=\frac{\sin(ik)\sin(jk)}{\sin k}
 =\sum_{r=0}^{\min(i,j)-1}
 \sin[(i+j-1-2r)k].
\]

The last equality is a finite trigonometric identity. `F` vanishes at both
integration endpoints. Since
`d exp(-2igt cos k)/dk = 2igt sin(k) exp(-2igt cos k)`, integration by parts
gives

\[
 |U^{(\infty)}_{ij}(t)|
 \le\frac1{\pi gt}\int_0^\pi|F'_{ij}(k)|\,dk
 \le\frac{2ij}{\pi gt}.
\]

For the last inequality, each positive-integer cosine has absolute integral
`2` over `[0,pi]`, and the positive frequencies in the sum add to `ij`.
Writing

\[
 S_2=\sum_{i=1}^{m}i^2=\frac{m(m+1)(2m+1)}6,
\]

the system-block Frobenius norm obeys

\[
 \|U^{(\infty)}_{SS}(t)\|_F\le\frac{2S_2}{\pi gt}.
\]

## 5. A finite-reservoir certificate with exact rational gates

For any `i,j` in the first `m` sites, a walk contributing to a matrix power
cannot distinguish the finite path from the half-line until it crosses the
missing edge and returns. The shortest such walk has length at least

\[
 K_*=2(L+1).
\]

Consequently the `ij` moments of their hopping matrices agree through order
`K_*-1`. Using the power series of the exponential and the norm bound `2g`
for both matrices gives, uniformly on `0<=t<=T`,

\[
 |U^{(L)}_{ij}(t)-U^{(\infty)}_{ij}(t)|
 \le2\sum_{q=K_*}^{\infty}\frac{x^q}{q!},\qquad x=2gT.
\]

For `x<K_*`, the exponential/Chernoff estimate of this tail is

\[
 2\sum_{q=K_*}^{\infty}\frac{x^q}{q!}
 \le2\left(\frac{ex}{K_*}\right)^{K_*}.
\]

The registered family chooses

\[
 L=64S_2,\qquad gt\in[8S_2,16S_2],\qquad gT=16S_2.
\]

Here `x=32S2` and `K_*=128S2+2`, so `4x<K_*`. The elementary bound `e<3`
therefore supplies the exact rational upper bound

\[
 \epsilon:=2\left(\frac34\right)^{K_*}.
\]

There are `m^2` system-block entries, so their error matrix has Frobenius
norm at most `m epsilon`. Combining the estimates,

\[
 \boxed{\|A(t)\|_F\le a(t)
 :=\frac{2S_2}{\pi gt}+m\epsilon,}
\]

\[
 \boxed{P_{\rm order}(t)=P_{\rm count\ correct}(t)
 \ge1-\frac12a(t)^2.}
\]

Replacing `pi` by the strict rational lower bound `3.14` makes the right-hand
amplitude bound weaker and entirely rational at the registered integer
times. `model.bounds` compares these rational numbers directly with the
frozen gates. For every scored size below, throughout its complete window,

\[
 \boxed{P_{\rm order}>0.9968,\qquad
        P_{\rm count\ correct}>0.9968.}
\]

The rounded lower display is `0.9968305002231328`; rounding is not used for
acceptance. These are single-read bounds at every chosen time in a finite
window, not a joint probability of never leaving the ordered family under
continuous monitoring. Repeated reads require their measurement backaction
and apparatus to be included. The bounds do not extend to arbitrarily
late time. As an additional first-member check, the sharper geometric tail
bound is an exact rational comparison:

\[
 \frac{2\,32^{130}}{130!}\frac1{1-32/131}
 <2\times10^{-24}.
\]

It is stronger than the coarse `2(3/4)^130` used in the general gate.

| Physical spins n | Walls m | S2 | Reservoir spins L | Certified gt window | One-particle sites |
|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 64 | 8 to 16 | 65 |
| 3 | 2 | 5 | 320 | 40 to 80 | 322 |
| 4 | 3 | 14 | 896 | 112 to 224 | 899 |
| 5 | 4 | 30 | 1920 | 240 to 480 | 1924 |
| 8 | 7 | 140 | 8960 | 1120 to 2240 | 8967 |

The proof applies to this connected family; these are not independent
copies of a two-spin system. Reservoir size and certified waiting time grow
as `O(n^3)`. The full active Hilbert dimension remains `2^(n+L)`; the small
propagator calculation is an exact free-fermion reduction, not a claim of
small physical resources.

## 6. A specified endpoint disturbance and autonomous recovery

At `tau=8S2/g`, apply the physical endpoint flip `X_(n-1)`. It is logical
`X` at wall site `k=m-1`; the fixed Hamiltonian then acts for another
`s=8S2/g`. This is an externally supplied disturbance, not part of the
autonomous baseline.

Under this spin flip, Jordan-Wigner annihilation operators transform as

\[
 c_j\mapsto c_j\ (j<k),\qquad
 c_k\mapsto c_k^\dagger,\qquad
 c_j\mapsto-c_j\ (j>k).
\]

Before the flip the state conserves particle number. Let `C=pU(tau)PU(tau)^dagger`
be its normal covariance and set

\[
 D_0=\operatorname{diag}(\underbrace{1,\ldots,1}_{m-1},0,-1,\ldots,-1).
\]

Direct operator conjugation gives the exact normal covariance immediately
after the flip:

\[
 C'=D_0CD_0+(1-C_{kk})|k\rangle\langle k|.
\]

The flip can introduce anomalous pair correlations. Accordingly, the
post-flip proof uses expected wall number and the union bound; it does not
apply the earlier number-conserving vacuum determinant without justification.
Number-conserving propagation evolves the normal covariance by `U(s)C'U(s)^dagger`
even when anomalous correlations are present.

Write `D_0=-I+PBP`, where

\[
 B=\operatorname{diag}(2,\ldots,2,1),\qquad\|B\|\le2.
\]

The group law for `U` now gives the exact expression

\[
 \boxed{\begin{aligned}
 \mathbb E N_{\rm wall}^{\rm flip}(\tau+s)
 ={}&p\|-A(\tau+s)+A(s)BA(\tau)\|_F^2\\
 &+(1-C_{kk})\|A(s)|k\rangle\|^2.
 \end{aligned}}
\]

In particular, the first squared norm contains an interference cross term.
Its sign is not fixed, and it must not be dropped or treated as a positive
causal contribution. A conservative triangle inequality retains it:

\[
 \mathbb E N_{\rm wall}^{\rm flip}(\tau+s)
 \le p[a(\tau+s)+2a(s)a(\tau)]^2+a(s)^2.
\]

For `tau=s=8S2/g`, the entire calculation lies within the same finite-bath
comparison horizon `T=16S2/g`. Thus

\[
 \boxed{P_{\rm order}^{\rm flip}(2\tau)
 \ge1-\tfrac12[a(2\tau)+2a(\tau)^2]^2-a(\tau)^2
 >0.9922.}
\]

The conservative rational evaluation displays `0.9922835622546656`. This
proves recovery after one specified physical endpoint flip; it is not a
claim about arbitrary disturbances, repeated attacks, or unlimited repair.

The original decoder also needs the intervention accounted for. The flip
changes `Q` by `+1` if its wall was empty and by `-1` otherwise. Hence the
known-intervention decoder `K_guess=N_B-1` is correct whenever the wall was
empty at the flip and the system is ordered at the final readout. Its success
is at least the recovery lower bound minus `p a(tau)^2`, approximately
`0.9891140624777984` using the same conservative constants. The subtraction
and its external knowledge are part of this diagnostic; the unperturbed
identity `K_guess=N_B` is not silently transferred to a changed charge sector.

## 7. A count-preserving bath phase echo

This counterexample tests whether the written count is a sufficient physical
description of future behavior. Let the full hopping-path staggering be

\[
 D_{jj}=(-1)^j,\qquad DTD=-T.
\]

At time `tau`, apply only the reservoir part of this diagonal phase operation:
`D_B` equals the identity on system sites and equals `D` on reservoir sites.
It leaves the complete physical-system reduced state unchanged, as well as
the full distribution of reservoir occupation readouts and their joint
distribution with the inert count reference. Every individual occupation
marginal and the `N_B` count record are therefore unchanged at the kick.
It can change off-diagonal coherences of the complete reservoir reduced
state; this is not a comparison with both subsystem density matrices fixed.

Let `D_S` be the restriction of `D` to the system sites. Since
`U(tau)DU(tau)=D` and `D_B=D+P(I-D_S)P`, the exact system block after a further
time `tau` is

\[
 \boxed{A_{\rm echo}=D_S+A(\tau)(I-D_S)A(\tau).}
\]

The undisturbed block at this time would be `A(2tau)`. They need not agree.
For `m=1`, `D_S=1` and `A_echo=1` exactly: the system's wall population returns
to its initial `p=1/2`, even though its high-accuracy count record at `tau`
was unchanged by the phase kick. This supplies an explicit failure of
count-only future prediction. It does not establish a record-reading agent.

The phase kick is external coherent control. It preserves `Delta Q`. For
this initial diagonal preparation, nearest-neighbor covariances are purely
imaginary under the bipartite hopping evolution, so the mean hopping energy
is zero before and after the real diagonal sign kick. Its mean work is zero,
but its energy variance and microscopic coherences can change; zero mean
work is not a fully resource-matched intervention.

## 8. Energy, entropy capacity, and controls

The active initial state has `n=m+1` bits of entropy, all in the physical
spins; the reservoir is pure. These are entropy capacities and preparations
inside the closed system, not an external entropy sink. Its initial energy
moments are

\[
 \langle H\rangle=\Delta m/2,\qquad
 \operatorname{Var}(H)=\Delta^2m/4+g^2m/2.
\]

To see the second expression, the diagonal term contributes
`Delta^2 m p(1-p)`. A system bond contributes `2g^2 p(1-p)` and the interface
contributes `g^2 p`; other reservoir bonds annihilate the vacuum. Cross terms
have zero expectation in the initial occupation-diagonal state. Setting
`p=1/2` gives the displayed result.

Removing only the interface retains all degrees of freedom, preparations and
the mean energy, but reduces this variance by `g^2/2`. Wall number is then
conserved within the physical chain, so its probability of zero walls stays
`2^(-m)`. This identifies the energy/information transport channel; changing
that channel is not an isolated causal test of a memory readout. Making the
whole active system maximally mixed instead leaves it invariant under every
unitary and cannot produce the stated ordering. That control also changes
the bath entropy and energy resources.

For a direct entropy accounting, let `lambda_i` be the eigenvalues of
`A A^dagger`. The nonzero reservoir covariance eigenvalues are
`p(1-lambda_i)`. At `p=1/2`,

\[
 S_S=1+\sum_i h_2(p\lambda_i),\qquad
 S_B=\sum_i h_2[p(1-\lambda_i)],\qquad
 S_{SB}=m+1.
\]

These formulas show the ordering process transfers entropy into the prepared
reservoir while global entropy remains fixed. More generally, if ordering
has error at most `delta<1/2`, its two-dimensional ordered family implies

\[
 S_S\le1+h_2(\delta)+\delta\log_2(2^m-1).
\]

Subadditivity and `S_SB=m+1` then require

\[
 \boxed{S_B\ge m-h_2(\delta)-\delta\log_2(2^m-1).}
\]

The initially pure reservoir provides this capacity explicitly. The inert
reference has its own initially prepared classical correlation and zero
Hamiltonian; it neither supplies a dynamically generated controller nor
changes the active-state energy/entropy accounting above. External kicks
also require preparation, timing and control resources. The endpoint flip
can inject energy; the model contains no accounting rule that makes it free.

## 9. Finite closure excludes an exact attractor or perfect plateau

For any finite Hamiltonian and fixed observable, the expectation is a finite
sum of exponentials of energy differences. If it has a late-time limit, all
nonzero-frequency coefficients must vanish: multiply by a frequency phase
and take a long-time average, or average the squared deviation from the
putative limit. The expectation must then have been constant at every time.
Finite unitary evolution also has arbitrarily close recurrences, returning
all these expectations close to their initial values. The initial ordering
probability here is `2^(-m)<1`, so permanent attraction to probability one is
impossible.

The expectation is real analytic in time. A probability exactly equal to
one on any open time interval would therefore equal one identically, also
contradicting the initial state. The proven high-probability finite windows
avoid both obstructions; neither a perfect plateau nor permanent memory is
claimed. Recurrence bounds the interpretation of the reservoir record.

The external Schrödinger parameter, time origin, boundary geometry, local
interaction rule, cold reservoir, conserved orientation seed and allowed
readouts remain specified inputs. Mathematical completion of the scoped
finite-window gates does not discharge the stricter target-free emergence,
active developmental feedback, novel-physics or cosmological gates.
