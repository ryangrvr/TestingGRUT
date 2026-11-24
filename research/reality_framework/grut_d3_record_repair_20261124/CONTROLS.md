# Matched controls and what they establish

The historical bit r is fair. The common timed SWAP moves r from the system
S into the inert fault register D and prepares S in its erased state e.
No dynamical term acts on D or the coin C. All three conditions start
from the same physical preparation and undergo this same fault.

| Condition | Local operation on M after the fault | Information in the memory label | Historical joint success |
|---|---|---|---|
| Intact | Identity | Correct label whenever written | J |
| Scrambled | Flip the two labels when the same fair coin is 1 | Independent of past r when C is unread | J/2 |
| Wrong | Flip the two labels in both coin states | Inverted label whenever written | 0 |

The flip leaves a blank memory blank. It changes the correlation with D
without discarding a degree of freedom. The wrong-record condition retains
the association in an inverted code. The scrambled condition hides it behind
the unobserved coin. The isolated coin remains in every condition. Together
with a nonblank M, it recovers the label. The Hamiltonian has no such access.

## Present state and energy matching

Let A denote all evolving degrees of freedom S, M, W and B. Given r, the
state just after the fault and before the memory control is

$$
 |\phi_r\rangle =
 |e\rangle_S\!
 \left(a_0|\mathrm{vac}\rangle_M|0\rangle_W
      +\sum_{x=1}^{L_W}a_x|r\rangle_M|x\rangle_W\right)|0\rangle_B .
$$

The amplitudes depend on w and are the same for either r. The memory flip
maps `phi_r` to `phi_(1-r)`. It therefore preserves
`(rho_phi0+rho_phi1)/2`, including conditional on either coin value.
Tracing out D gives exactly the same full density on A and C in every
condition. The joint state of D and S just after the operation also agrees
across conditions. Probes confined to A have identical unconditioned
statistics during the subsequent evolution. Historical recovery is a joint
readout that includes D and is recorded without postselection.

On this preparation F_M acts exactly as G=F_S F_M because S=e. The full
color symmetry commutes with H. For each full-Hamiltonian spectral projector
`Pi_E`, each r, and each coin value, we have

$$
 \langle\phi_r|F_M^\dagger\Pi_E F_M|\phi_r\rangle
 =\langle\phi_r|G^\dagger\Pi_E G|\phi_r\rangle
 =\langle\phi_r|\Pi_E|\phi_r\rangle .
$$

This includes energy degeneracies and all moments. The memory-only flip
need only coincide with the symmetry on this input. One could realize the
equivalent symmetry operation on S and M on the same input, or make it
conditional on C. Ordinary energy measurement probabilities are matched.
Implementation of the diagnostic unitary still uses external control.

Immediately before the fault, mean energy is Delta and variance is g^2.
The common SWAP increases Q_B from zero to one. In the writer path
the even and odd positions have a relative imaginary phase from the
bipartite hopping, so the mean of every hopping edge is zero before and
after the fault. Mean energy rises by Delta=4g. The numerical ledger
records the post-fault variance, which is shared by all arms at each w.
Energy distribution matching compares the arms after the common fault.

## Accessible record and resource accounting

With D read as a fair bit and M read in its pointer basis, the intact
pre-repair mutual information is `1-|a_0|^2` bits. A blank outcome conveys
no bit and has equal probability for both earlier labels. A nonblank value
determines r. The wrong-label condition has the same information value and
opposite code. The scrambled pointer has zero mutual information with D.
The entire active state conditional on r is also independent of r in that
condition once C is unread. Hence any guess of a fair r made solely from
A has success at most one half.

In the intact and wrong conditions the full A states conditional on the
two r values overlap by `|a_0|^2`, since they share the blank component.
For their equal mixture the two possible nonzero eigenvalues are
`(1+|a_0|^2)/2` and `(1-|a_0|^2)/2`. Sector calculations for each
conditional label respect this shared physical corner. The total
experiment starts with two bits of von Neumann entropy from the fair
historical label and fair independent coin. Unitary evolution, the fault
SWAP and the memory control preserve that entropy. D, the memory blank
and the reservoir vacua require preparation. The demonstrated episode
uses 2688 W spins and 64 B spins with vacuum initialization. Their sizes
are conservative choices for the analytic bounds.

## Negative results

- **Same ordinary configuration probability.** The three conditions have
  identical probability of reaching the allowed family {0,1}. A measurement
  of generic order therefore has zero ability to distinguish which history
  was recovered. The change is in the joint D, S, M outcomes.
- **Short chains.** With one spin per chain, the specified timing restores
  S with certainty and returns M to the blank state. Joint success is
  exactly zero. It is necessary to test retention independently of repair.
- **Unflagged disturbance.** A flip 0↔1 on S with B still empty leaves Q_B=0.
  It creates no erasure flag and S is frozen in the other configuration.
  The repair mechanism only addresses the specified erasure.
- **Finite windows.** The bounds apply to a single terminal read at any
  allowed pair of times. Repeated measurements require an instrument model.
  Further faults, shrinking the reservoirs, removing the external controller
  and permanent repair require separate analysis.

The comparisons test causal use of a local classical record. The local
copy and repair rules, initial configuration family, fault symbol and
reservoir preparation enter as inputs. Autonomy here is the fixed
Hamiltonian between preparation, the common disturbance and readout.
The model's rules already choose the kind of structure that it can restore.
