# Exact derivations and controlled limits

## 1. Ancillary channels from the finite unitary

For a degenerate cell in $I/2$, the cross term of
$U=\cos\alpha\,I-i\sin\alpha\,X\otimes X$ vanishes on
tracing the cell. The channel is
$(1-p)\rho+pX\rho X$, $p=\sin^2\alpha$.
For the resonant exchange with a cell of excited probability $q$,
the zero- and two-excitation basis states are unchanged and the
one-excitation subspace rotates by α. Direct trace of that unitary
gives
$$
 x\mapsto (1-p)x+pq,\qquad
 \rho_{01}\mapsto\sqrt{1-p}\,\rho_{01}.
$$
The mean bare heat to the cell is $\Delta p(x-q)$. The sign
uses the pre-contact payload population. The entire unitary
preserves excitation energy on the resonant link.

Let $\delta=T/n$, $P_\delta=I+\kappa\delta(X-I)$, and
$D_\delta=\operatorname{diag}(e^{-2ig\delta},e^{2ig\delta})$.
For initial sign $+$, define
$$
 z_n=\boldsymbol 1^\top(P_\delta D_\delta)^n e_+,\quad
 d_n=D_0(1-2\mu\delta)^n,\quad
 \lambda_n=(1-\gamma_h\delta)^n(1-\gamma_c\delta)^{n_c}.
$$
The other initial sign yields $\overline z_n$. The complex
factor tracks the entire source history. Replacing the history
by its final state would generally change the phase. The record
contrast about the *current* sign is $d_n(1-2\kappa\delta)^n$.
Population transitions are independent of the signs and record.
The conditional decoherence factor includes $\sqrt{\lambda_n}$.

## 2. D7 and the joint recovery optimum

Take a qubit channel whose two population transition probabilities
$u=\Pr(1|0)$, $v=\Pr(0|1)$ are independent of the hidden
alternative, with $\lambda=1-u-v\ge0$. Its off-diagonal
multiplier is $\alpha_b$, and the accessible programme
has density $\rho_b$ independent of the input conditional on b.
Define $\bar\rho=\sum_b\pi_b\rho_b$ and
$C=\sum_b\pi_b\alpha_b\rho_b$.
The Bell input's reference-labelled diagonal blocks before
recovery have Q probabilities $(1-u,u)$ and $(v,1-v)$
and programme $\bar\rho$. Its off-diagonal block is $C/2$
between $|00\rangle$ and $|11\rangle$.

For any Stinespring isometry $V:Q\otimes M\to O\otimes E$,
put $V_i=\langle i|_O V |i\rangle_Q$.
Isometry on each input Q sector gives the fidelity
$$
 F=\tfrac14\left[
 u+v+\lambda\,\operatorname{Tr}\bar\rho(V_0^\dagger V_0+V_1^\dagger V_1)
 +2\operatorname{Re}\operatorname{Tr}(V_1^\dagger V_0 C)\right].
$$
Both $V_i$ are contractions. When $\lambda\ge0$ the diagonal
terms are bounded by $1+\lambda$, and the coherent term by
$2\|C\|_1$. A controlled polar unitary on the programme
and unchanged computational populations attain both bounds.
The partial isometry in the polar decomposition of C may be
extended to a unitary. Hence
$$
 \boxed{F_{\mathrm{opt}}=(1+\lambda+2\|C\|_1)/4.}
$$
With $\lambda=1$ this is the exact D7 theorem
$\boxed{F_{\mathrm{opt}}=(1+\|C\|_1)/2}$.
Averaging the programme when it is inaccessible replaces
$\|C\|_1$ by $|\operatorname{Tr}C|$. Neither optimum
uses the raw probe, current source, fresh access to past cells
or the Bell reference. This comparison allows arbitrary output
channels at the specified interface.

In I1 let $a=\operatorname{Re}z_n,\ b=\operatorname{Im}z_n$.
The pointer M has two equally probable values. Its weighted
coherence operator C is diagonal with entries
$\sqrt{\lambda_n}(a\pm i d_n b)/2$.
The known phase correction conditional on M attains the trace
norm. In particular
$$
 F_M=\tfrac14(1+\lambda_n)+
     \tfrac12\sqrt{\lambda_n}\sqrt{a^2+d_n^2 b^2},\qquad
 F_0=\tfrac14(1+\lambda_n)+\tfrac12\sqrt{\lambda_n}|a|.
$$
The blind optimum is idle for the registered phase range because
each trajectory has absolute total rotation at most $\pi/4$.
Every trajectory's off-diagonal phase has nonnegative real part.
With fixed actuator error ±ε, the off-diagonal correction term
is multiplied by $\cos(2\epsilon)$. On this grid it is
nonnegative. Each pointer value has the same benefit, and
holding an independent idle action available yields
$$
 R_{\rm usable}=\tfrac12\sqrt{\lambda_n}\,
     [\cos(2\epsilon)\sqrt{a^2+d_n^2b^2}-|a|]_+ .
$$
This describes the registered restricted control. A different
noise model or actuator set changes that optimum.

The ideal gain is quadratic in small d when $a\ne0$, with
leading term $\sqrt\lambda\,d^2b^2/(4|a|)$. At exact $a=0$
it is linear, $\sqrt\lambda\,|d\,b|/2$. The small-contrast
crossover depends on independently determining $a$, since
a small detuning of the deep benchmark restores the quadratic
regime. The zero-gain boundary for the restricted actuator is an
ordinary optimal-policy switch.

## 3. Continuous-time evolution and error

Let $\mathcal X\rho=X\rho X$ and
$\Gamma_q=(1-q){\cal D}_{\sigma_-}+q{\cal D}_{\sigma_+}$, with
${\cal D}_L\rho=L\rho L^\dagger-\{L^\dagger L,\rho\}/2$.
The exact flip and exchange channels above are respectively
$\exp[-\tfrac12\log(1-2p)(\mathcal X-I)]$
and $\exp[-\log(1-p)\Gamma_q]$.
Indeed these give the corresponding affine population and
coherence factors. In the small-contact limit the sum generator is
$$
 {\cal L}(t)\rho=-ig[Z_BZ_Q,\rho]
 +\kappa(\mathcal X_B-I)\rho+\mu(\mathcal X_M-I)\rho
 +\gamma_h\Gamma_{q_h}\rho+\gamma_c(t)\Gamma_{q_c}\rho .
$$
Here $\gamma_c(t)=\gamma_c$ before the supply stop and zero
afterwards. All rates are nonnegative, so each time segment is
GKSL and the product is causal and CPTP on B,M,Q. R is isolated.
The Q-only process need not be CP-divisible after hiding B.
The continuous phase factor is
$$
 z(T)=\boldsymbol1^\top e^{T[\kappa(X-I)-2igZ]}e_+ .
$$
Writing $\nu^2=\kappa^2-4g^2$ gives the equivalent
$e^{-\kappa T}[\cosh(\nu T)+(\kappa-2ig)\sinh(\nu T)/\nu]$,
with the continuous limit at $\nu=0$. The classical thermal
population obeys $\dot x=\gamma_h(q_h-x)+\gamma_c(t)(q_c-x)$.
Within a segment $x(t)=x_*+(x(0)-x_*)e^{-(\gamma_h+\gamma_c)t}$,
$x_*=(\gamma_hq_h+\gamma_cq_c)/(\gamma_h+\gamma_c)$ when
the denominator is positive. Integrating this law gives the
signed heat, which is separate from the Bell score.

For the operator error, take diamond norms of superoperators.
$\|\mathcal X-I\|_\diamond,\|\Gamma_q\|_\diamond\le2$
and $\|[-igZ_BZ_Q,\cdot]\|_\diamond\le2|g|$.
For $p<1/2$, the rate replacements use
$$
 -\tfrac12\log(1-2p)-p\le p^2/(1-2p),\quad
 -\log(1-p)-p\le p^2/[2(1-p)].
$$
Multiply each by the corresponding generator norm. For a
fixed segment, the difference between a product of ideal local
semigroups and the exponential of the sum is bounded per step by
$\delta^2\sum_{i<j}\|{\cal L}_i\|_\diamond\|{\cal L}_j\|_\diamond$.
Duhamel's identity twice gives the commutator bound
$\delta^2\|[{\cal L}_i,{\cal L}_j]\|_\diamond/2$
for each pair. Products and convex combinations of the CPTP
semigroups are contractive, then the commutator norm is at
most twice the product of norms. Telescoping the errors across
contacts and segments gives the implemented $O(T\delta)$
bound when rates are fixed and the stop is aligned. Scores
with effects in [0,1], including their optima, differ by at
most half the channel bound. The optimizer may change among
the fixed feasible set on either input. No bound for reusing
correlated ancillas follows from this product approximation.

As $\delta\to0$ the angle of each collision is $O(\sqrt\delta)$.
A shrinking physical pulse interval can require a growing
peak strength, and the number of prepared cells per unit time
grows as $1/\delta$. The fine grid is a convergence check
of the declared family. It is not a fixed finite-resource
limit operating forever.

## 4. Entropy and matching

For a finite block S and initially product Gibbs bath cells
$\otimes_j\tau_j$, unitary evolution gives the exact identity
$$
 \Delta S(S)+\sum_j\beta_j Q_j
 =D\left(\rho'_{S\mathcal E}\,\middle\|\,\rho'_S\otimes_j \tau_j\right)
 \ge0 .
$$
Here $Q_j$ is the increase in a bath's bare energy.
It follows by expanding the relative entropy and using
total entropy invariance. The right side contains both
system–bath correlations and the deviation of the baths
from their initial product. For the Q marginal alone the
fresh-cell identity applies contact by contact, because a
new cell is a product with Q before it interacts. Phase and
bit contacts do not change the marginal x or its entropy.
The stored marginal ledger is that sum for Q. It differs
from the entropy of B,M,Q,R or the total closed system.
The separate density reviewer also verifies the local
entropy identity including the retained block.

With interactions off at the endpoints, all registered
unitaries preserve the sum of bare energies, so bare
net work is zero in that ideal episode. Supplying the
Gibbs cells, fresh gaps, mixed controls, Bell pair, clock
and low-entropy memory need not be work free. At hot
temperature the non-equilibrium free energy of one clean
cold Gibbs cell relative to a hot one is
$\beta_h^{-1}D(\tau_c\|\tau_h)$. A stock of $N_c$
cold qubits accepts at most $N_c(1-q_c)\Delta$ net energy.
Reusing a small stock changes the channel, and no bath
can retain an indefinitely positive net energy flux at
fixed finite capacity. Scoped statements about finite
entropy production do not price a universal $R$ unit.

## 5. Organization, collective decay and network access

For N independent copies with final damage probability
x, exact restoration has probability $(1-x)^N$.
A coarse damage-density criterion $\alpha$ has
$\Pr[\operatorname{Bin}(N,x)\le\lfloor N\alpha\rfloor]$.
For $x\ne\alpha$, concentration bounds its wrong-side
probability by $e^{-2N(x-\alpha)^2}$, with the usual
choice of open or closed inequality. The response as x
crosses α is a readout cutoff. The independent channel and
mean x are smooth. No collective interaction produced it.

For two physical qubits and vacuum cells, a common
lowering operator $\sigma_-^1+\sigma_-^2$ has a bright
eigenmode with exchange amplitude $\sqrt2$ and a
decoupled dark one in the single-excitation sector.
An equal one-site defect prior is a half-bright/half-dark
mixture. In the calibrated mean-rate match its survival is
$[1+(1-2p)^n]/2$. Jensen's inequality bounds this below
by the independent survival $(1-p)^n$, so the rate-matched
aggregate cannot beat independent repair. Pure bright
preparation can beat it. The peak-rate match weakens the
common arm without changing the dark obstruction. Two
vacuum cells are budgeted per step in each arm. The
interaction norms are distinct and the calibration angles
in the output expose that distinction.

For the XY study, the source bit is vacuum vs a single
local excitation. Each bond is bounded below by
$-|J_{ij}|(n_i+n_j)$, so the chosen gap $2J_G$ makes
the vacuum a ground state in all compared networks. Let $u(t)$ be the one-body hopping
propagator. On a region S the alternative to vacuum
has an orthogonal excitation component of probability
$p_S(t)=\sum_{i\in S}|u_{i0}(t)|^2$. Its trace distance
from vacuum is exactly $p_S(t)$. The guessing advantage
above $1/2$ is $p_S(t)/2$. If a spectral eigenvalue
has projector $P_e$, the long-time Cesàro probability
at i is $\sum_e|(P_e)_{i0}|^2$. For the uniform ring
each plane-wave energy is at most twofold degenerate,
so the average of $p_S$ is at most $2|S|/N$.
The infinite uniform lattice has the Fourier dispersion
$E(k)=\Delta_G+J_G\cos(ka)$. Fixed-site amplitudes
decay by oscillatory integration. The quadratic expansion
near an extremum has remainder at most
$J_G|ka|^4/24$ relative to that expansion point.
This yields a conditional low-wave-number envelope
under an already specified lattice.

For the complete graph the symmetric subspace has
eigenvalue $J_G$ and its orthogonal complement has
$-J_G/(N-1)$. A specified remote-site transfer is at
most $4/N^2$. The initial site retains high probability
from the degenerate sector. These exact comparisons
show that preparation-specific retention and transmission
can give different rankings at the same incident coupling
sum. They do not select geometric locality. [Results](RESULTS.json)
