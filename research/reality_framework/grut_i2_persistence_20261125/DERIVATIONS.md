# Derivations for GRUT-I2

Work throughout in units of $k_B T_b$ and use the fixed
time-homogeneous jump law in [MODEL.md](MODEL.md). The
jump model assumes an overdamped thermal solvent and rare
activated hops. Local detailed balance fixes rate ratios.
The barrier prefactors, gap and concentration regime are
physical inputs. The model does not derive their values.

## Microscopic reversibility and finite fuel

For a particle–vacancy exchange at positions i,i+1, only the
two exterior contacts can change. Their old and new counts
give $d=B(y)-B(x)\in\{-1,0,1\}$. At the thermal exchange the
energy change is $-Jd$. A common barrier prefactor and an
activation energy symmetric about the initial and final
energies give the stated rates with ratio $e^{Jd}$. The
patch conformation changes only the common prefactor for the
reverse pair. A local fuel-assisted jump also takes one
specific fuel molecule from F to W for $d=+1$, with the reverse
using W. Total energy change is $-(J+A)d$. The individual
molecule's rates therefore have ratio $e^{(J+A)d}$.
At a sensor jump the energy change is $-A d_s$. The reverse
retains the same e and uses the opposite $d_s$. Source
and thermal receptor flips have zero energy change.

The full labelled chain obeys
$$
 \frac{w(z\to z')}{w(z'\to z)}=e^{E(z)-E(z')}.
$$
The activation factors are an ordinary reversible reaction
model. Only one local reactive patch samples the fuel and
receptor. There is no feedback from the score, no spatial
broadcast and no external switch clock.

All states with a given particle arrangement, e, m and
fuel count have the same aggregate rates to every other
count state. Summing the rates over the $k$ F or $K-k$
W choices proves strong Markov lumpability. Exchangeability
of the initial labelled stock implies uniform states within
each count class at later times. The count equilibrium is
$$
 \pi(z)=Z^{-1}g(k)e^{-E(z)},\qquad g(k)=\binom Kk.
$$
For an assisted $d=+1$ hop with the same e and m in both
directions, the count ratio is
$$
 \log \frac{w_{\rm a}(x,k\to y,k-1)}
               {w_{\rm a}(y,k-1\to x,k)}
 =J+A+\log\frac{k}{K-k+1}.
$$
A nonzero actuator coefficient changes the activity but
cancels in the ratio. The same counting law holds for
sensor conversion with $A$ replacing $J+A$. Every pair
satisfies detailed balance even when the controller is
reading the receptor. With finite stock, any steady
distribution on communicating components has zero stationary
jump currents. Transient work comes from the charged initial
stock. Changing barriers at fixed energy does not add a
new affinity.

## Observables and propagation

Let $Q$ act on probability columns. For $p_0$ over the
retained states, $p(t)=e^{tQ}p_0$. Let $O$ be the shape
class and $Q_{OO}$ its principal block. Its diagonal
keeps the original escape rate. Summing the mass that
has never left O gives
$$
  S_O(t)=\boldsymbol 1^\top e^{tQ_{OO}}p_{0,O}.
$$
The endpoint mass $\boldsymbol 1^\top p_O(t)$ also counts
reentry. The work budget is computed on the full $p(t)$,
including failures, and never on the survivors alone.
Bonds and longest component are physical configuration
functions, invariant under rotations of the ring.
The threshold treats many arrangements as the same shape.

For a signed increment h on a directed jump, put
$r_h(z)=\sum_{z'}w(z\to z')h(z,z')$. The expected accumulated
increment is $\int_0^t r_h^\top p(v)dv$. A block with rows
$r_h^\top$ below Q yields the integrated fuel use of each channel.
The joint history statistic is obtained by resolving the
unmeasured ensemble at the tag time $t_*$ by e, then
propagating both pieces. It gives the actual
$I(m_t:e_{t_*}\mid e_t)$ under the trajectory law. This
tag is an analysis variable. The model reads the current
physical m alone. It never consults a future event.

The labelled check uses positive uniformization on the
species-resolved generator. With $\nu$ at least the
maximum escape rate, $P=I+Q/\nu$ is stochastic, and
$P_{OO}=I+Q_{OO}/\nu$ is substochastic. The propagator is
the Poisson-weighted sum of their powers. Truncating at
N omits at most $\Pr[\operatorname{Pois}(\nu t)>N]$ in
mass. Each signed reward grows in absolute value at
most linearly in the number of uniformized steps, with
coefficient $\|r_h\|_\infty/\nu\le1$. The corresponding
tail is bounded by the Poisson first tail moment.
The tolerance is separate from floating-point error.

## Cluster scaling without a phase transition

Write $n=L/2$. A configuration with c occupied runs has
c positive particle-run lengths summing to n and c positive
gap lengths summing to n. Choosing a labelled start and
dividing by the c choices of start for the same pattern gives
$$
 N_c=\frac Lc\binom{n-1}{c-1}^2 .
$$
It has $B=n-c$. Thus the canonical adjacent probability ratio
is
$$
 \frac{\pi(c+1)}{\pi(c)}
   =e^{-J}\frac{(n-c)^2}{c(c+1)} .
$$
For fixed finite J the mode fraction approaches
$f_*=(1+e^{J/2})^{-1}$. Away from $f_*$ by any fixed
fraction the consecutive ratio is bounded away from one
on the appropriate side. Comparing with a point halfway
back to the mode gives exponential suppression at each
fixed displacement. The equilibrium density has a smooth
limit at every finite J. A sharp histogram cutoff does
not mark a singularity in the interactions.

Conditional on c, a uniformly selected run has the
composition length law. The probability that at least
one occupied run has length at least m is bounded by
$$
 c\,\frac{\binom{n-m}{c-1}}{\binom{n-1}{c-1}} ,
$$
clipped at one and with zero numerator outside support.
For $c/n$ bounded from zero this decays exponentially
in m up to the factor c. The largest run is $O(\log L)$
with high probability at any fixed finite J. This is a
correlated nearest-neighbor system. No independent-site
assumption enters the run calculation.

For the driven sizes, one patch and one stock remain
fixed. This changes the fuel capacity per site. A general
fixed-time comparison uses shared Poisson opportunities
on every background bond, each with rate ceiling
$\bar r=r_{\rm bg}e^{J/2}$. Couple the acceptance draws
in two variants. Initially mark the two patch sites.
A background bond can disagree only when its four-site
rate support intersects a marked site. Each affected
opportunity marks at most two further sites, and each
marked site touches at most four background supports.
Therefore the expected marked count is at most
$2e^{8\bar r t}$. A marked site can change at most two
bonds. At identical initial populations the difference
in expected bonds **per ring site** is at most
$\min(1,4e^{8\bar r t}/L)$. This bounds what one local
patch can change asymptotically in fixed time. It does
not give a driven bulk phase diagram.

## Cost frontier for the frozen controls

For each setting $j$ the exact Markov solution supplies
first-exit reliability $S_j$ and reset work $W_j$ under a
declared reset convention. A trial-level mixture with
retained independent setting label has expected
$(S,W)=\sum_j p_j(S_j,W_j)$. Minimizing work at reliability
at least r is a two-constraint linear program, with
normalization and one payoff constraint. At an optimum
one support point suffices when the payoff constraint is
slack. When it binds, an extreme solution has at most two
nonzero weights. A direct proof uses a null vector of
the two rows $(1,S_j)$ on any support of greater size, then
moves to a lower-cost boundary until one weight is zero.
The registered pairwise interpolation enumerates all such
extremes and the pure settings. The same result holds
for maximizing survival at a work cap. The formulas give
the exact constrained optimum when the inputs are exact.
The numerical evaluation has the precision and scope in
[REPRODUCTION.md](REPRODUCTION.md).

This optimization is restricted to one frozen finite
rate menu, the specified preparations, access and reset
conventions. There is no strictly positive model-independent
energy requirement set by duration and reliability alone.
For example, increasing a common activation barrier
can lengthen survival without supplying an ongoing
chemical force. The passive hop rates and temperature
are held fixed inside the registered cost comparison.
The thermodynamic infima are derived in
[THERMODYNAMICS.md](THERMODYNAMICS.md).
