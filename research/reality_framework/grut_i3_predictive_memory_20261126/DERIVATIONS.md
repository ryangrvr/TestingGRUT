# Finite-time guarantees

The delayed two-state benchmark, the first-hit comparison, the
rate-distortion condition and the analytical primary floors are derived
in the frozen feasibility gate. Here the exact finite dynamics yield the
remaining values and specify how model error enters the test.
[Feasibility gate](FEASIBILITY.md)

## Catalytic yield and optimal use of the record

Let relative state $0$ mean a match. After the orientation transfer,
$B\oplus A$ is a telegraph process of rate $r=q+\nu$, independent
of the bond transition noise. Along a relative path $z(t)$, the
conditional bond probability $x(t)$ solves
$$
 \dot x=c_{z(t)}(\pi-x),\qquad x(0)=0.
$$
Set $Q=\left(\begin{smallmatrix}-r&r\\r&-r\end{smallmatrix}\right)$
and $D=\operatorname{diag}(c_h,c_l)$. Killing at the time-dependent
rate $c_{z(t)}$ gives the survival probabilities
$$
 k_j=\mathbf{1}^{\mathsf T}e^{(Q-D)T_r}e_j,\qquad
 s_h=\pi(1-k_0),\qquad s_l=\pi(1-k_1).
$$
To order the paths for the two initial conditions, couple a mismatched
pair $(0,1)$ to jump together into $(1,1)$ at rate $r$ and into
$(0,0)$ at rate $r$. After coalescence, use common flips. Each
marginal has the correct rate, and the first path spends at least as
much time matched. There is strict inequality in the survival integrals
on positive-probability paths when $c_h>c_l$ and $T_r>0$.

At the transfer, a symmetric record with posterior match probability
$p$ gives pre-noise advantage $(2p-1)(s_h-s_l)$ for the majority
orientation. The independent transfer flip reduces this advantage by
$1-2\eta$. Thus the record-driven permutation implements the optimal
decision over all locked orientation choices available through $M$
under the frozen access model. The score is
$$
 S_M=s_l+(s_h-s_l)p_A.
$$
A policy that reads more registers or acts through a new response
interface has a larger decision space. The no-record comparison already
grants a larger timing and latching space than the sensing arm uses.

At the primary row the required locked match is .807226 and the required
information for the action is .202927 nats. The delivered actuator
matches with probability .890263 and has .347182 nats of mutual
information with the current conformation at entry. The inverted arm
has the same marginal information and an adverse decision rule.
Information quantity and causal use are separately testable.
[Results](RESULTS.json)

## Independent scaling

For $N$ isolated compartments using the same local policy and
independent sources, baths and internal noise, write $K$ for the
number of completed material bonds. For uniform success probability
$p$, the exact reliability is
$\Pr(K\ge\lceil r_\star N\rceil)
=\sum_{j=\lceil r_\star N\rceil}^{N}{N\choose j}p^j(1-p)^{N-j}$.
For independent null policies with different scores each at most $p_0$,
couple each outcome with an independent uniform draw to dominate $K$
by a binomial with parameter $p_0$. This gives the reported null
upper tail. The weak law gives convergence to one for the predictive
probability and to zero for the null when $p>r_\star>p_0$. Exponential
tails follow by optimizing the Bernoulli moment generating function.
The number of ports, independent environmental conditions, and fresh
preparations grows linearly with $N$. Correlated compartments and a
fixed budget shared by increasing $N$ are outside this guarantee.
[Scaling results](RESULTS.json)

## Perturbations beyond the registered noise

For two jump descriptions on the same states with rates $w$ and
$\widetilde w$, couple equal transitions at the minimum of the two
rates as long as their paths agree. A phase of duration $t$ adds at
most $t\sup_i\sum_{j\ne i}|w_{ij}-\widetilde w_{ij}|$ to the probability
that the paths differ. Initial total variation, the phase terms, and
the maximum conditional total variations of altered gate kernels add
over the episode. This is a conservative finite-time bound on the
error of the bond score. A short duration error can be included by
bounding any extra jump by the maximum outgoing rate times the duration
error. The bound assumes a jump description with the declared
observable states. Resolving bath memory, transit-time information or
additional observables requires a measured model or an additional
coupling argument.
