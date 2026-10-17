# Why global allocation differs from a charge-by-charge cap

This analysis follows the original S8 stop and precedes G1 trajectories.
It does not reclassify an original failure or assert a new formation result.

The total Schmidt rank is the sum of ranks in the orthogonal charge blocks.
With a fixed total budget, the best Frobenius approximation keeps the globally
largest singular values. Reserving at most 16 per block is an additional
restriction: it can discard information while other blocks leave total-rank
capacity unused. S8's observed total rank 55 at a per-charge-16 stop makes
this a concrete numerical question. G1 changes that allocation rule only,
with a separate preregistered gate and unchanged total rank 96.

## An exact allowed-state capacity witness

For one root component, enumerate the nonempty type histories u. For each,
put |u|-1 photons in a receiver mode left of a chosen cut and one photon in
a distinct orthogonal right-hand mode b_u. There are enough physical receiver
modes in the original L_N law. Form the normalized coherent state

\[
|\Psi\rangle={1\over\sqrt n}\sum_{u\ne\varnothing}
 |u, (|u|-1)\hbox{ left photons}\rangle\otimes|1_{b_u}\rangle.
\]

Every term has the exact conserved total M: history charge M-|u| plus
|u| photons. Left charge is the same M-1 in every term. Left histories
and right photon modes are orthogonal, so the single charge block has n
equal Schmidt values 1/sqrt(n). Each *conditional history* has Schmidt
rank one across the same cut.

For N=5 remaining counts (2,2), n=18. A 16-vector charge cap necessarily
discards weight 2/18=1/9 (unnormalized norm 1/3). For N=6 counts (2,3),
n=33, so it discards weight 17/33 (norm about .7177). A total-rank-96
allocation represents either witness without truncation.

These are **allowed Hilbert-space packets, not states shown to form from
the neutral preparation**. They prove a capacity distinction, not why a
particular numerical trajectory stopped, a general no-go, or a physically
selected geometry. G1 must still demonstrate target-time performance.

`global_controls.py` checks these exact counts and compares its allocation
proposal against globally sorted singular spectra. The actual algorithm
accepts a truncation only after the same positive tail/arithmetic account.
Singular-value sorting is a proposal, not a separate accuracy certificate.
