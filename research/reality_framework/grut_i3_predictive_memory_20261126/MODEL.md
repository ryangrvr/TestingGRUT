# I3 — Screened conformer-dependent assembly

## Physical system and scope

One material unit can form a reversible bond in a catalytic pocket. A passing reagent has
two thermodynamically degenerate conformers, \(B=0,1\). The pocket accepts one of two
catalyst orientations, \(A=0,1\). A match lowers the common activation barrier of bond
formation and breakup. Both conformers have the same binding free energy. The two favorable
responses therefore compete for the same catalytic port over a finite interval.
A separate upstream contact can couple the reagent to a bistable receptor \(M\). A
screened interval separates receptor contact, orientation transfer and catalytic contact.

The model resolves classical chemical microstates and uses a time-dependent contact
schedule in a single thermal bath at \(T_b\). Its kinetic description applies when
intrawell and bath relaxation are fast relative to the residence times and reaction
rates. The activation barriers are specified below. The model takes the activated
jump law and the phase preparation as inputs. It describes a class of screened
chemical or bistable-well devices, with no calibrated reagent or apparatus assumed.

## States, energies and supply

The state is \((B,M,F_s,A,C,X,F_r)\) in \(\{0,1\}^7\).
\(F_s\) is a sensor fuel molecule, \(F_r\) a chemically distinct assembly fuel molecule,
\(X\) the presence of a physical material bond, and \(C\) an independent orientation-error
bit. A charge of one in a fuel bit denotes one available molecule.
In units of \(k_BT_b\) the bare energy is
\[
  E = a_s F_s + a_r F_r - JX .
\]
\(B,M,A,C\) have degenerate levels. Initially \(B,M\) are independent fair bits,
\(C\) is an independent Bernoulli bit with parameter \(\eta\), and
\(F_s=F_r=1,\ X=0,\ A=0\). All arms start with the same physical carriers and
capacity. The shuffled-record control uses a fair \(C\) and is charged for
its actual preparation. Chemical selectivity keeps \(F_s\) out of the material
port. The production trajectory stays in \(X+F_r=1\). Separate loading contacts
for reset can reach both values of \(X\) and \(F_r\) independently.

Time is in units of one nominal catalyst residence time. All rates are inverse
times in those units. \(B\) flips symmetrically at rate \(q\) throughout.
Two degenerate wells and a transition-state barrier
\(k_BT_b\ln(\kappa_B/q)\) above their minima with attempt rate \(\kappa_B>q\)
give this rate in the activated jump law. \(q=0\) denotes a closed barrier.
The thermal flips of \(M\) and \(A\) specified below have the same construction
with their respective barriers. The exit rates of the reagent do not depend on
memory, actuator or binding. This excludes backaction on the source in the rate
model and is a physical condition on a realization.

## Phase 1 — Receptor contact, duration \(t_s\)

The source and receptor undergo their independent flips at \(q,\mu\).
The only other edge, for each \(b\), is
\[
 (B=b,M=1-b,F_s=1)
 \mathrel{\mathop{\rightleftarrows}^{s e^{a_s/2}}_{s e^{-a_s/2}}}
 (B=b,M=b,F_s=0).
\]
A conformer-sensitive transition state opens this edge and selects which receptor
flip it catalyzes. With attempt rate \(\kappa_s\) its energy above the mean of
the endpoint energies is \(k_BT_b\ln(\kappa_s/s)\). The forward-to-reverse
ratio is \(e^{a_s}\). The reagent is returned in the same conformer during
a sensing jump. Uncatalyzed receptor flips and reverse sensing can undo a record.
The sensor has one molecule of capacity and no refill. The same molecule may
be recharged by the reverse step. The sensor is then detached from the reagent.

## Phase 2 — Screening and orientation transfer, total delay \(\tau\)

\(B\) continues to flip at \(q\). \(M\) flips at \(\mu\) for \(\tau-\delta\).
The source remains screened from both receptor and actuator. During the last
\(\delta\), \(M,A,C\) are isolated from their thermal flips. Two equal-duration
local pulses implement \(\mathrm{SWAP}_{MA}\) followed by a flip of \(A\)
controlled by \(C\). The resulting classical permutation has
\(M'=A\), \(A'=M\oplus C\). \(B\) keeps flipping during both pulses. The
inverted control complements the actuator map by using \(1-C\) as the
active control state. The passive arm skips the map and retains its record.
The disabled arm removes only the catalyzed sensing edge and performs the map.
The shuffled arm performs the map with independent fair \(C\).

A microscopic realization of the classical permutations on degenerate registers
uses isolated pulses \(H_{\rm int}(t)=g(t)G\) with \(G=\mathrm{SWAP}_{MA}\)
or \(G=\lvert 1\rangle\!\langle1\rvert_C\otimes X_A\), and
\(\int g(t)dt/\hbar=\pi/2\) for each pulse. Pulse phases are immaterial to
the diagonal input and to the subsequent population dynamics. Bare energies
and total entropy are unchanged. For reversible drive of this form the net
work on the isolated bits vanishes because \(\langle G\rangle\) is constant
during each pulse and \(g\) starts and ends at zero. Transient pulse energy
returns to the drive. The common control allowance in the protocol covers
unrecovered control, timing, routing and reset overhead within the declared
model budget. A realization must separately verify that overhead and the
assumed shielding and pulse rates.

## Phase 3 — Catalysis, duration \(T_r\)

\(B\) and \(A\) flip independently at \(q,\nu\). The receptor and its fuel are
detached. Only
\[
 (X=0,F_r=1)\mathrel{\mathop{\rightleftarrows}^{c\pi}_{c(1-\pi)}}
 (X=1,F_r=0),\qquad
 \pi=(1+e^{-(J+a_r)})^{-1}
\]
can change the material. \(c=c_h\) when \(A=B\) and \(c=c_l\) otherwise.
These rates follow from a common activated prefactor
\(c/[2\cosh((J+a_r)/2)]\) and the two opposite half-affinity factors.
In particular both directions speed up by the same factor. There is
at most one net assembly fuel molecule used. Both source labels
are interchangeable. The physical output is the event \(X=1\) at the fixed
deadline, with no conditioning on the reagent, fuel or data readout.
This is a finite-time assembly reliability observable.

## Access and matched controls

The initial fixtures, environmental path distribution, initial reagent, reaction
times, fuel capacity and maximum binding rate are common. The acquired information
available in the response phase is physically the orientation transferred from
\(M\). Other device registers and the supplied clock have no reagent-conformation
information on the branch prior to the first bond. Failed attempts leave no
observable trace other than the time elapsed without a bond. The catalyst
changes the transition rate and cannot separately signal the conformation.
A null may sense at the upstream contact. Any causal remnant that it can use
later, including orienting an actuator early or caching a filtered arrival
time, counts as retained information. The no-record null has the same present
material state and same physical access to the reaction outcome. It has no
source-correlated remnant from upstream and no extra current-state sensor.
Relative actuator noise can be ignored or compensated at zero cost in the
relaxed null, which strengthens the comparison.

The bound grants either orientation at any time and even a fractional split
\(u\in[0,1]\) of the single port. The total kinetic budget gives
\[
  c_{B=0}(u)=c_l+(c_h-c_l)u,\quad
  c_{B=1}(u)=c_h-(c_h-c_l)u,\quad
  \overline c=(c_h+c_l)/2.
\]
Any open-loop contact modulation that turns the port off is also admissible.
The memoryless material-feedback bound additionally grants cost-free, immediate
latching upon first binding. Independent randomization and any fixed control
schedule are allowed. Both kinds of null may use the entire common work
budget. A second full-rate port, a different affinity, and any extra
conformation-reporting channel change the comparison. A free perfect cue
at the start of the response provides a separate diagnostic at the same
reaction and timing parameters. The additional information channel is explicit.
