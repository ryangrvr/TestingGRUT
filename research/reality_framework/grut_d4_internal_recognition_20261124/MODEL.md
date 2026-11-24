# GRUT-D4 internal discrepancy detection and historical repair

24 November 2026. One fixed Hamiltonian writes two local replicas from a source
bit, then responds to discrepancies among the three. The historical bit is
prepared as a fair classical mixture. It supplies a score for the experiment
and does not enter the dynamics as a separate reference.

## Physical model

The sites are a bit S, two three-level memories M1 and M2 with labels 0, 1
and b (blank), a writer spin W, and two finite spin chains B1 and B2. All
interactions are time independent. Initially S is r, both memories are blank,
W is 0, and both chains are vacuum. At most four sites participate in each
interaction. All sites and their geometry are supplied.

On the eight memory-ready basis states define

$$
 r=\mathrm{majority}(s,m_1,m_2),\qquad
 x_i=s\oplus m_i\quad(i=1,2).
$$

This change of coordinates is bijective on those eight states. Its inverse
is `s=r xor (x1 and x2)` and `m_i=s xor x_i`. The historical score uses the
initial S value. The coordinate r in the Hamiltonian is always the majority
of the current state. No register in the interaction is set by the initial
value. The parity bits correspond to two physical comparisons.

Let `p_i` project onto ready states with x_i=1. Let `A_i` lower x_i from 1
to 0 at fixed r and fixed other parity, summed over those two labels and
both settings of the other parity. These operators vanish on states that
contain any blank. For each chain define `b_(i,j)=|0><1|` and
`N_i=sum_j b_(i,j)^† b_(i,j)`. Set hbar = g = 1, Delta = 4g, and
`kappa=g/65536`. The Hamiltonian is

$$
\begin{aligned}
H={}&\Delta\big(P_{\mathrm{both\ blank}}+n_W\big)
+\kappa\sum_{r=0}^1\big(
 |r,r,r,1\rangle\langle r,b,b,0|+\mathrm{h.c.}\big)\\
&+\sum_{i=1}^2\bigg[
 \Delta(p_i+N_i)+g\big(A_i b_{i,1}^\dagger+\mathrm{h.c.}\big)
+g\sum_{j=1}^{L_i-1}\big(b_{i,j}^\dagger b_{i,j+1}+\mathrm{h.c.}\big)
 \bigg] .
\end{aligned}
$$

The writer kets are in S,M1,M2,W order. The first exchange writes both
memories while moving one excitation from the blank to W. It can run in
reverse. All other terms dissipate a parity discrepancy into the finite
spin chains by reversible transport. Label inversion of S, M1 and M2
commutes with H, with b left fixed. S is distinguished by the comparison
geometry. Neither historical label has a privileged energy or rule.

The charges `Q_W=P_both_blank+n_W` and `Q_i=p_i+N_i` are conserved between
faults. The writer exchanges a blank state and a ready **unanimous** state
only. The corrected sectors preserve the majority. On blank states this
same coordinate is assigned the present S value. With vacuum correction
chains and Q_i=0, writing is the exact two-state Rabi evolution and is
perfect at `kappa*w=pi/2`. The proof also scores a window
`|kappa*w-pi/2| <= 1/512`. The nominal writing time is
`32768*pi/g`. The slow coupling and long preparation are explicit resources.

## Disturbances, scores and clock

Faults are basis-label exchanges X on the selected site. The memory
exchange fixes b. Faults are chosen independently of r and do not set
an erasure symbol, feed a syndrome to H, or touch W or either chain.
An ideal laboratory pulse specifies when and where a test fault occurs.
After writing the same H continues to act. Its correction exchange
responds to x_i immediately upon a discrepancy, with no timed trigger.
The model supplies no independent physical clock or pulse source.

For the separated-fault experiment, begin with a pulse after writing and
use up to four single-site pulses, each followed by a wait of
`128/g` to `144/g`. The next pulse can act on any site including the
previous one. Both B chains remain in the system throughout. We score
one terminal read in a run, with no intermediate measurement, reset, or
conditioning on success. Interactions run continuously through each wait.

- Historical restoration is S equal to the initially prepared bit.
- Record retention is both memories equal to that bit.
- Joint success is the whole ready triple equal to the historical codeword.
- Detection uses the parity pair, with 00 for no discrepancy. In the
  single-fault class S, M1 and M2 faults give 11, 10 and 01, respectively.
- A joint resource score requires joint success and the expected cumulative
  count in each bath given the applied pulse masks. These counts can store
  the cumulative parity excitations. They do not encode fault order in
  general.
- The undamaged score tests spontaneous repair-chain activity and changes
  of the historical bit. We track writer reversion to blank separately.

## Proof and numerical scope

The registered large chains each have 2304 spins. The full physical
dimension is `36 * 2^(L_1+L_2)`. For analysis, the ready subspace without
the writer exchange is equivalent to two independent parity-source XX
paths, with an inert logical majority between faults. During faults that
coordinate has an exact transformation. The small production matrices
retain the writer and every other physical term. The large numerical
comparison is exact for the ready model with the writer exchange
omitted after acquisition. Its norm error relative to the physical
writer-on trajectory is bounded analytically. This is a comparison for
the proof, rather than a switch used to implement the experiment.

Standard unitary dynamics, the parity comparisons, majority correction,
blank memories, vacuum baths, pulse type and timing, and the assumption
of sufficiently separated single faults are declared. Two simultaneous
flips can leave a discrepant majority of the opposite color. Three flips
produce the opposite unanimous word. The historical value is then
ambiguous without additional evidence.
