# D6 exact finite-system result

Let $s_b=(-1)^b$, $p_\theta=(1+\sin(2\theta))/2$, $m=3n+4$ and
$L=2m$. For the primary intact arm of [MODEL.md](MODEL.md) at
$T=\pi/\omega$, the probability of arriving with both target
layouts and a correct physical forecast is $p_\theta$ in each
open condition $b,z$. This achieves the optimal response given
the declared single phase query. In the sealed control the
structural score cannot exceed $1/2$ for the fair source prior.
Both statements include all trials.

## Exchange clock and terminal map

In the conserved one-head sector the $j$-th forward exchange is
$|j+1\rangle\langle j|\otimes U_j$. Set $V_0=I$,
$V_{j+1}=U_jV_j$ and $G=\sum_j|j\rangle\langle j|\otimes V_j$.
All $U_j$ commute with $K$, since the source and strength labels
are stationary and none of the exchanges acts on W. Therefore

$$
G^\dagger H G=\omega J_x\otimes I+I\otimes K,
$$

where $J_x$ acts on a spin of length $L/2$ in its ordered $J_z$
basis. The matrix element $\langle j+1|J_x|j\rangle$ is
$\sqrt{(j+1)(L-j)}/2$. For each pure preparation in the coin
mixture, the state is

$$
\sum_{j=0}^L a_j(\phi)|j\rangle V_j e^{-itK}|\chi_0\rangle,\qquad
a_j(\phi)=\sqrt{\binom Lj}[-i\sin(\phi/2)]^j[\cos(\phi/2)]^{L-j},
\qquad \phi=\omega t .
$$

Represent the spin as the symmetric sector of $L$ qubits. The
initial spin is the all-zero state and its $J_x$ rotation is a
product of $L$ identical $X$ rotations. The coefficients follow
from the binomial multiplicities. At $\phi=\pi$ the entire
amplitude is at $j=L$. Output probabilities are those of the
ordered reversible map. At intermediate times several positions
can be occupied coherently. The clock law is independent of
condition, control arm and tape data.

## Phase acquisition, witness and later use

The state of P just past the probe link is
$|\psi_b\rangle=e^{-is_b\theta Z}|+\rangle$. Its Y expectation
is $s_b\sin(2\theta)$ and its overlap with the other source's
state is $\cos(2\theta)$. The rotation $R_x(\pi/2)$ maps the Y
decision basis to the Z pointer basis. The copy and subsequent
swap prepare

$$
\sum_{r=0,1}\alpha_{b,r}|r\rangle_M|r\rangle_D,\qquad
|\alpha_{b,b}|^2=p_\theta,
$$

with P restored to zero. The dump keeps the original decision
distribution through the downstream interventions. Tracing
it out gives a diagonal pointer on M. Averaging the fair
source yields $I(b:M)=1-h_2(p_\theta)$ before the intervention.
This is 0.1887218755 bits at weak strength.

Conditioned on an open source $b$, W independently evolves
as $e^{-it\kappa s_b Z_W}|+\rangle_W$. Its Y outcome decodes
to $b$ with probability $(1+\sin(\omega t/2))/2$. At $T$ this
probability is one. In the sealed arm W decodes to 0 at $T$
for both $b$. No head exchange acts on W.

Both banks use the same decision $d$. Bank B also uses the
independent challenge $z$, loaded after A, and a parity
computed by toggling the reset P. After both banks are full,
both words are correct exactly when $d=b$. Thus their errors
are shared. In open arms the event also gives the correct
W forecast at $T$. For the main port, averaging the fair
source, the intact, scrambled, inverted and disabled
structural scores are $p_\theta$, $1/2$, $1-p_\theta$ and
$1/2$. All forward response links preserve the post-intervention
pointer. The passive port keeps M but selects a fair independent
bit and scores $1/2$. The open direct port uses a copy of
$b$ loaded at the source and scores 1 even when M is disabled.

## Accessible input and optimality

At ordered position 3, tracing out the source station and
intervention station leaves the head and developmental
substrate. In the main arm only M depends on $b$. A and Q
are blank, P was reset by its swap with D, and the future
challenge is independent of $b$. The mask, original dump
value, and any disabled memory remain at the stations.
The complete response input thus has trace distance
$\sin(2\theta)$ between the source hypotheses when intact
or inverted, and zero when scrambled or disabled.

The same factorization includes coherent response positions.
For $j\ge3$ write $V_j=W_jV_3$, where $W_j$ never acts on
the two stations. The map forming
$\sum_{j\ge3}a_j|j\rangle W_j$ has squared norm operator
$p_{\ge3}I$, where $p_{\ge3}$ is independent of $b$.
Dividing this map by $\sqrt{p_{\ge3}}$ gives a fixed isometry,
including if the data is entangled with an excluded station.
The proof applies to any subset of response positions with
nonzero clock probability. It does not require a measurement
of the head to run the protocol.

For equal priors and accessible input states $\rho_0,\rho_1$,
the mutually exclusive successful output effects of any
fixed response channel induce a decision measurement on
the input. Completing failure outcomes into guesses gives

$$
P_{\rm success}\le P_{\rm guess}
 =\frac12+\frac14\|\rho_0-\rho_1\|_1.
$$

In particular the success of a binary guess effect $E$ is
$\tfrac12[1+\mathrm{Tr}\,E(\rho_0-\rho_1)]$, maximized by the
positive eigenspace projector. The declared single phase
query has trace distance $\sin(2\theta)$. Even a pure probe
entangled with an initially independent reference cannot
improve on this for $0\le\theta\le\pi/4$: the magnitude
of the conditional state overlap is
$|\cos(2\theta)+i\sin(2\theta)\langle Z_P\rangle|$,
which is at least $\cos(2\theta)$. Convexity extends this
bound to mixed inputs. The registered decoder attains it.

In the sealed sector the whole Hamiltonian factors as
$I_B\otimes H_{\rm sealed}$, and all degrees of freedom other
than B start independently of $b$. Every operation on
that accessible complement, including repeated internal
probes, produces identical accessible histories. A
physical coupling to B would grant new access. This includes the perfectly discriminating
alternative furnished by the open direct port.

## Finite time, stability and recurrence

In $|\omega t-\pi|\le0.2$ let
$q=\cos^2(\omega t/2)\le \sin^2(0.1)<0.01$.
An unfinished bank requires at least $m+1$ head-count
failures. Using the binomial law,

$$
\Pr(j<m)\le 2^L q^{m+1}
 \le 0.01(0.04)^m = \varepsilon_n
 \le 1.048576\times10^{-16}.
$$

The weak intact joint structural and forecast probability,
including every unfinished trial, is at least
$\frac34(1-\varepsilon_n)\frac{1+\cos(0.1)}{2}>0.747$.

To bound a measurement that also accesses the earlier clock
sectors, project each pure conditional state onto $j\ge m$
and normalize. Its trace distance from the original is at most
$\sqrt{\varepsilon_n}$. Mixing coins and reducing to accessible
degrees cannot increase this difference. The binary guessing
advantage of the entire clock and the response degrees over
the completed sector is at most $\sqrt{\varepsilon_n}$ in this
window, with the two stations excluded. This statement
concerns the unperturbed history. Earlier clock sectors can
carry probe information, and new interventions can re-enter
the probe station.

For the intermediate diagnostic take $p=(3n+7)/(2L)$. The
A-only probability is exactly
$\sum_{j=n+3}^{2n+4}\binom Lj p^j(1-p)^{L-j}$.
The n=2 value is 0.660416206084, below the 0.90 gate.
The n=64 value is 0.999867119390. A Chebyshev lower bound is
$1-4Lp(1-p)/(n+3)^2$, already above 0.90 at n=64.

At $2T$ the head returns to zero and the data to $V_0=I$.
Both tapes are blank. W is $|-\rangle$ up to phase, so its
Y distribution returns but its state differs. At $4T$ the
witness propagator is $-I$ and the whole density operator
recurs exactly. The apparatus supplies no lasting output
or accumulated memory across these cycles.

## Energy and entropy

The localized clock has binomial weights over $J_x$.
The witness starts along +X and therefore assigns equal
independent weight to its energy levels $\pm\kappa$, for
each source. The joint law is
$E_{k,w}=\omega(k-L/2)+w\kappa$ with probability
$2^{-(L+1)}\binom Lk$. All arms have mean zero and variance
$\omega^2L/4+\kappa^2$. Each exchange has norm at most $J$,
and the one-head Hamiltonian has norm below $2J$.
There is no externally performed work after preparation.
Correlation transfers use the explicit finite D, E and coin
registers, the global entropy is conserved, and the clock
can undo the episode. Preparation and reset costs are priced
separately in [RESOURCES.md](RESOURCES.md).
