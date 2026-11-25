# D7 exact recovery bound and finite implementation

The recovery task has a fixed input qubit Q, an inaccessible reference
R and an accessible programme register. The hidden condition b changes
the sign of a phase rotation on Q. The register available to the
recovery has state $\rho_b$. The payload and programme are prepared
independently conditional on b, and the programme gets no access to
R. These conditions matter for the bound.

## Recovery theorem

Write $s_b=(-1)^b$, let the source have probabilities $\pi_b$, and
put $\bar\rho=\sum_b\pi_b\rho_b$. The phase disturbance is
$e^{-is_b\phi Z_Q}$. Define the programme operator

$$
C=\sum_b\pi_b e^{-2is_b\phi}\rho_b .
$$

For any recovery channel acting on Q and the accessible programme,
the maximum overlap of the output and R with the original Bell state
is

$$
\boxed{F_{\rm opt}=\frac{1+\|C\|_1}{2}.}
$$

The same formula holds for any finite mixture of diagonal qubit
phase rotations by using their respective off-diagonal factors in C.
It allows coherent processing of the programme and arbitrary
discarded recovery ancillas.

**Upper bound.** In the two reference-labelled components 00 and 11,
the incoming joint density operator has diagonal programme blocks
$\bar\rho/2$ and off-diagonal blocks $C/2,C^\dagger/2$. Let V be a
Stinespring isometry for the recovery from Q and the programme to
output O and an environment. Let $V_i=\langle i|_O V |i\rangle_Q$.
The overlap after recovery is

$$
F=\frac14\left[
\operatorname{Tr} V_0\bar\rho V_0^\dagger+
\operatorname{Tr} V_1\bar\rho V_1^\dagger+
2\operatorname{Re}\operatorname{Tr}(V_1^\dagger V_0 C)
\right].
$$

The other output blocks do not contribute to this Bell projection.
Both $V_i$ are contractions, since V is an isometry. The two
diagonal terms are each at most 1, and the trace inequality gives
$\left|\operatorname{Tr}(V_1^\dagger V_0 C)\right|\le\|C\|_1$.

**Attainment.** Write $C=W|C|$ with W extended to a unitary on the
programme space. Preserve the Q basis and apply $W^\dagger$ to the
programme for Q=0 and identity for Q=1. After tracing the programme,
both diagonal terms equal 1 and the last trace equals $\|C\|_1$.
Equivalently, measure the spectral projectors of $W^\dagger$ on the
programme and apply the corresponding phase to Q=0. This yields the
same recovered channel. The optimizing programme measurement can
depend on the phase ensemble.

## Application to the D6-type record

Let the source be balanced and put $D=\sin(2\theta)$ with
$0\le\theta\le\pi/4$. The one-use assay states
$|\psi_b\rangle=e^{-is_b\theta Z}|+\rangle$ have overlap
$c=\cos(2\theta)$. Their Y expectations are $s_bD$. The prepared
$R_x(\pi/2)$ rotates the Y measurement onto the Z pointer. Copying
to M and dumping the old probe leaves a binary memory with
$\Pr(g=b\mid b)=(1+D)/2$ in the intact main arm. Its mutual
information with b is $1-h_2((1+D)/2)$ under the fair prior.

For a classical programme consisting of this pointer, C is diagonal
with entries $(\cos(2\phi)\mp iD\sin(2\phi))/2$. If the entire
unmeasured probe instead serves as the programme, its C has rank at
most two. In the span of the probes,
$\operatorname{Tr} C^\dagger C=(1+c^2\cos(4\phi))/2$
and $|\det C|=(1-c^2)/4$. The square of the sum of the two
singular values is $\operatorname{Tr} C^\dagger C+2|\det C|$, including
the rank-one limit, both programme uses give

$$
F_{\rm opt}(\theta,\phi)
=\frac12\left[1+\sqrt{\cos^2(2\phi)+D^2\sin^2(2\phi)}\right].
$$

Thus the specific D6 pointer already attains the coherent-probe
optimum in this family.

The resulting compensation is local. Conditional on b and the
latched setting g, the payload rotates by
$e^{-i(s_b\phi-s_g\beta) Z}$, leaving Bell fidelity
$\cos^2(s_b\phi-s_g\beta)$. Averaging the intact record gives

$$
F(\beta)=\frac12\left[
  1+\cos(2\phi)\cos(2\beta)+D\sin(2\phi)\sin(2\beta)
\right].
$$

Maximizing aligns the two-dimensional unit vector
$(\cos(2\beta),\sin(2\beta))$ with
$(\cos(2\phi),D\sin(2\phi))$, giving the frozen calibrated angle
and the upper bound above. When the latter vector vanishes all
angles tie. The blind optimum uses identical programme states for
both source values, so $F_{\rm blind}=(1+|\cos(2\phi)|)/2$.
The blind comparator already allows *any* downstream channel on
the phase-disturbed qubit and source-independent ancillas.
It has the same one-qubit input interface as the assisted arm.

For the weak and shallow setting $D=1/2,\phi=\pi/8$ the calibrated
angle is $\tfrac12\arctan(1/2)$. The fidelity is
$(1+\sqrt{5/8})/2$. The full choice $\beta=\phi$ achieves 7/8.
Inversion replaces D by −D while holding β fixed. The fair
scramble and passive selector replace D by 0. A disabled selector
is fixed at zero, which yields the same balanced score at that
angle. These three comparisons that hold β fixed do not replace
the optimized blind bound. The direct port copies b and achieves
unit fidelity with $\beta=\phi$.

For the deep setting $\phi=\pi/4$, blind averaging removes
the entire off-diagonal part of Q. Its channel is complete
dephasing, and the Bell score cannot exceed 1/2. The weak assayed
record yields 3/4. This score certifies entanglement, since a
separable state has overlap with a Bell state at most 1/2.
The fixed-replacement null prepares some receiver state σ wholly
independent of R. It gives exactly
$\langle\Phi^+|\sigma\otimes I_R/2|\Phi^+\rangle=1/4$.

## Exact exchange-clock model

In the one-head sector set $U^{(0)}=I$ and
$U^{(j+1)}=U_j U^{(j)}$. The block-diagonal unitary
$G_{\rm clock}=\sum_j |j\rangle\langle j|\otimes U^{(j)}$
conjugates the physical Hamiltonian to $\omega J_x\otimes I$
with spin $L/2$. In particular each physical link has the
same coefficient as the ordered spin ladder, bounded by J.
From a localized head at 0 the probability of site j at
physical time t is

$$
\binom{L}{j} p^j(1-p)^{L-j},\qquad
p=\sin^2(\omega t/2).
$$

To see this, represent the spin as the totally symmetric
subspace of L spin-half systems, with initial state all zero.
The X rotation is a product, and counting flipped systems
gives the binomial distribution. At $T=\pi/\omega$ the head
arrives at L with probability one and the data unitary is
exactly $U^{(L)}$ up to phase. The data evolution at the
endpoint has the stated commitment, payload phase, compensation
and release order.

For every j≥8 the receiver has the endpoint payload and
reference joint marginal. For j<8 its own marginal is the
initial |0⟩, independent of the reference, whose marginal is
$I/2$. Tracing the head removes inter-position coherences.
The measured receiver probability at any time is therefore

$$
F_{\rm receiver}(t)=\frac14+
\Pr(j\ge8)\left(F_{\rm endpoint}-\frac14\right).
$$

This includes all clock branches. In $|\omega t-\pi|\le0.2$,
$1-p\le\sin^2(0.1)<0.01$, so the probability of an unfinished
release is at most $2^L(0.01)^{L-7}$. This yields
$6.5536\times 10^{-14}$ at L=16 and decreases with L. The
advantage over the matched blind clock scales by the same
ready probability. The unpadded L=8 clock has probability
$\cos^{16}(0.1)$ of being ready at the flanks and falls below
the frozen 0.95 gate. At $2T$ the spin and all data return
to the initial density operator, and the receiver Bell
score is 1/4.

## Resources and allowable information

The one-head sector has energy values $E_k=\omega(k-L/2)$.
The initial weights of the localized spin in the X basis
are $2^{-L}\binom Lk$ in every arm. The same product-spin
construction used for position probabilities proves this.
Mean energy is zero with the declared reference, the
variance is $\omega^2 L/4$, and the sector norm is below
$2J$. The time-independent evolution takes no external
work after preparation. An additional source-dependent control,
new query, alternate routing or access to the isolated R
would change the recovery problem.

The occupied head has 48 labels. The stationary Hilbert
dimensions at sites 1, 2 and 5 are 48, 16 and 16, respectively.
The receiver and isolated reference contribute two each.
The full dimension is $49^{L+1}\cdot49{,}152$, and the
one-head dimension is $(L+1)48\cdot49{,}152$. The common
inputs include clean registers, two mixed coins, and a
prepared Bell pair. For balanced b the initial ensemble
entropy is three bits, including the two coin bits, and is
preserved. Any discarded record remains in the finite dump,
disable sink or coin correlations. There is no dissipative bath.
Engineered exchanges, the prepared clock and Bell pair, control
calibration, readout and reset have separate laboratory costs.
