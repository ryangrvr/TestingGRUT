# GRUT-I3 independent internal review

This review and replay were authored by a separate internal AI. They follow
the freeze at b78f6199ec8a240e75e808cfd7204c915f5a6e54. No external RAI
session or hardware validation took place. The replay imports neither the
production solver nor its outputs.

## Replay and quantitative verdict

[critical_replay.py](critical_replay.py) propagates all 128 labelled states
of $B,M,F_s,A,C,X,F_r$ by forward uniformization. Each finite generator is
built by enumerating its allowed transitions. The two orientation pulses act
as a bijection. It propagates the two conditional catalytic yields from
separate labelled starting states. No renormalization or generator reduction
is used.

~~~bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 critical_replay.py --output /tmp/grut-i3-critic.json
~~~

The 32 rows and 160 physical arms took 0.317 seconds in one run and used 653
nontrivial propagations. The maximum omitted Poisson tail per propagation
was $9.029\times10^{-15}$. The maximum missing mass was
$1.80\times10^{-14}$. No negative component or violation of $X+F_r=1$
was observed. Log edge-ratio error was below $8.89\times10^{-16}$.
The fixed-environment acquisition and yield limits agree within
$4.45\times10^{-15}$. The yield mixture and retention predictions agree
within $5.89\times10^{-15}$ and $8.55\times10^{-15}$, respectively.
The maximum ledger identity error was $1.34\times10^{-15}$.
The registered comparison against the independent production implementation
passed all 2,624 scalar checks. The largest difference was
$6.28\times10^{-13}$ for separate reset work in the disabled arm at
grid point 04,04. [Comparison](CROSSCHECK.json)

| Row | Acquisition match | Stored match | Actuator match | Predictive bond | Feedback ceiling plus cue | Result at .70 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Primary | .949099 | .907175 | .890263 | .759906 | .669591 | Pass |
| Moderate | .936040 | .869790 | .846769 | .720035 | .662519 | Pass |
| Severe | .867029 | .670007 | .628278 | .557014 | .623012 | Miss |

Fourteen of thirty sweep points reach .70 with the configured actuator.
The primary disabled, passive and shuffled arms each yield .478354, and
the inverted arm yields .196801. These are configured controls. The stronger
relaxed feedback ceiling is the useful comparison for the necessity claim.

## What Jensen bounds

A controller can learn from the absence of a previous bond. Fix its independent
random seed and follow the unique hypothetical trajectory on which no bond
has happened yet. Along that trajectory time, its prior actions and its
observed material history specify an open-loop orientation schedule. Given the
exogenous conformer path, survival to the deadline without a hit is the
exponential of minus the integrated forward hazard on this schedule. The
expectation is over the unconditioned path law, even though a posterior inferred
from survival need not stay balanced.

At each time its exogenous marginal is balanced. The single-port rate budget
makes the expectation of the integrated forward hazard at most
$\pi\overline c T_r$. Jensen therefore bounds the first-hit probability,
and final bonding implies an earlier hit. The half split and perfect latch
attain the ceiling in the granted relaxation. They provide a stronger null
than the fixed arms or the actual reversible bond. The same argument for an
open-loop schedule without latching uses the two directions' common rate
multiplier and gives $\pi(1-e^{-\overline c T_r})$.

The cue correction applies to total variation of the joint source path and
extra signal against the corresponding product law. Its cost is added to an
expected payoff bounded between zero and one. A measured accuracy for a
single bit does not bound a general stream of control-dependent cues.
Reagent dwell, passive orientation-dependent scattering, early loading
times, and invisible sensing by the actuator would each need a measured
access or leakage model. Source backaction would also break the bound
if it changed the exogenous path law.

The stated detachment matters. For example, the fair C control destroys the
marginal predictor in A while A and C together reconstruct the transferred
record. Conditional sensing fuel state may also encode useful reliability
information in combination with the orientation. The admissible actuator
does not reopen these registers after transfer. A final correlation also
need not be a usable predictive record. The disabled primary arm has
final source–device mutual information .290519 nats after reaction, despite
an initial actuator match of .5. The material history itself can create
correlations. The ceiling handles precisely the declared observations.

## Thermodynamics and work

The detailed-balance ratio on each physical edge is the exponential of the
negative bare-energy change. Forward sensing spends one sensor fuel if it
occurs, and a reverse jump can restore it. Each net bond spends one distinct
assembly fuel. Since the pulses only permute the degenerate registers, heat
into the bath is $Q=a_s\langle1-F_s\rangle+(a_r+J)P(X=1)$. Prescribed changes
of barrier and coupling can still have implementation costs even when the
net reversible pulse work on the bits is zero.

Reset is priced on the device marginal and lacks access to B. The full
entropy production is $\Sigma=Q+S_{\rm full}^{\rm end}
-S_{\rm full}^{\rm init}$ in the adopted units. The source entropy stays
$\ln 2$ and its initial correlation with the device vanishes. As a result,
joint device reset work satisfies $W_{\rm reset}^{\rm joint}=\Sigma+I(B;d)$.
The separate reset charge adds total device correlation. Joint use of
device correlations assumes full abstract control on those registers
at reset. Neither ideal work includes a free recovery of the source–device
correlation.

For the primary predictive row, the fresh-preparation price is
$8.848712\,k_BT_b$, ideal joint reset is $5.877154\,k_BT_b$, and
separate reset is $6.443609\,k_BT_b$. The conservative fresh-preparation
bound including the assumed overhead is $9.946752\,k_BT_b$.
This verifies the ensemble budget within the frozen model. The benchmark
does not establish a lower actual work than every lower-yield null. The
same allowance covers work ignored by the microscopic switching description
and must be independently justified for any apparatus.

## Physical boundary

A real device would have to realize the common forward/reverse barrier
speed-up at approximately fixed binding affinity, a conformer-selective
receptor transition without appreciable source backaction, and shielding
through the transfer interval. The model uses a bath and activated-jump
rates from the outset. An attempt rate large enough for the stated
barriers and rapid intrawell equilibration, a calibrated conformation
leakage bound, and a work budget for contacts and reset remain conditions
to be established for the actual preparation. A controller that operates
a second full-rate port or changes the affinity defines a different benchmark.

The successful primary and moderate cases and registered replay agreement
earn **B, a predictive advantage explained by established physics**.
Necessity concerns the .70 reliability target in the screened, single-port
interface. It follows from ordinary jump kinetics, estimation, and
thermodynamic accounting. Independent compartments can show sharper
reliability scaling without a collective transition. The conventional
analytic explanation leaves no basis here for a distinctive physical law.
