# Neutral-origin recorded-history prototype

The original coherent neutral CP-12 ray evolves unforced until time 5. A rare
charged, fully circulating ring sector then has weight

\[
p_G=0.05502210766596276,\qquad
p_G\in[0.0550221026799,0.0550221126521].
\]

The state was not initialized with a loop filter, supplied host or flux pattern.
The following operation is an explicit diagnostic intervention. Record both the
labelled graph and the event flag with the complete channel

\[
\Phi(\rho)=\sum_g\sum_{b=0}^1P_{g,b}\rho P_{g,b},\quad
P_{g,1}=P_gG,\quad P_{g,0}=P_g(1-G).
\]

Here G means every active connection belongs to electric circulation and at
least one vertex fermion is occupied. At N=4,Q=4 this forces a spanning ring.
Every record outcome is retained; the G=1 outcomes are analyzed conditionally.
This recording changes the physics by removing coherences. It is not an
autonomously generated classical selection principle.

After the record, a gauge-invariant source exp(−iθn₀) is applied and the unchanged
Hamiltonian evolves for lag 0.5. Each graph record identifies the vertex opposite
source 0 on its initial ring. Its local density is measured independently of the
final graph. The result is not a final-distance reclassification.

| Quantity, conditional on initial G record | Result | Controlled interval |
|---|---:|---:|
| Same labelled host at final time | 0.999999997486 | [0.999999816238, 1] |
| Full circulation at final time | 0.893146518185 | [0.893146351599, 0.893146684773] |
| Opposite-receiver density susceptibility | 0.000131168601283 | [0.000131012542, 0.000131324661] |

The unnormalized opposite-receiver susceptibility is 7.21717290219×10⁻⁶.
Restricting its final detector jointly to the same initial host gives
7.21716051069×10⁻⁶. This confirms the response is carried by the retained-host
sector. The conditional response denominator is the *initial* G probability.
It has exactly zero source derivative because n₀ commutes with all record
projectors, so no final-event quotient derivative is being omitted.

The complete diagnostic record has a measured-system entropy of
4.08240090830 bits, with continuity error estimate 6.75397×10⁻⁷ bits, and adds
energy

\[
\Delta E_{\rm record}=0.0654739727471\ \pm\ 2.46046\times10^{-7}.
\]

This system-energy change is not a complete apparatus or reset-work budget.
A coherent apparatus can preserve purity of the enlarged system while the
measured system becomes mixed. Its physical construction is not supplied by
this calculation.

## Exact symmetry and numerical reduction

The waiting state evolves in the exact signed S₃ sector fixing source 0,
dimension 2,242. A labelled P_g is not invariant under that group and cannot be
implemented by selecting canonical representatives. The calculation therefore
constructs the complete physical E²≤256 basis, dimension 13,413, and embeds the
waiting state by a signed isometry J. The checked residuals are

\[
\|J^\dagger J-I\|=1.0514\times10^{-14},\qquad
\|H_{\rm phys}J-JH_{S_3}\|=1.1357\times10^{-14}.
\]

After recording, all branch vectors evolve in that complete physical basis,
including source-symmetry sectors which the record can populate. Every outgoing
transition omitted by E²≤256 is retained in the residual operator. This is a
numerical projection and never a modified physical law.

## Error control and resource accounting

The waiting wavefunction has full unbounded-law error e_w≤1.062315×10⁻⁸.
The map ψ↦(P_gGψ)_g is a contraction into branch-column Frobenius norm, as is
its composition with n₀. Thus the inherited waiting error is e_w, not the
square root of the number of branches times e_w. The lag propagation supplies
separate grouped Frobenius bounds for X and Y. For

\[
R=2\operatorname{Im}\sum_g\langle x_g,D_gy_g\rangle,
\]

the error obeys

\[
|\delta R|\le2(e_x\|\widetilde Y\|_F+
e_y\|\widetilde X\|_F+e_xe_y).
\]

Positive-probability amplitude intervals and division by the positive interval
for p_G control the normalized values. Survival intervals are intersected with
[0,1]. A separate 5×10⁻¹² numerical reporting tolerance is added. Calculations
use floating point and are not formal interval-arithmetic certificates.

The energy error does **not** bound an unbounded energy expectation by a mere
wavefunction norm. Electric, mass and diagonal rewiring terms commute with the
record. Therefore the actual cost observable is bounded:

\[
A=\Phi(H)-H=\Phi(K+R_{\rm off})-(K+R_{\rm off}),
\]

\[
\|A\|\le2\left(\sqrt{2NQ}|t|+
\frac{\Gamma Q(M-Q)}M\right).
\]

Its expectation error is at most ||A||(2e_w+e_w²), plus the stated reporting
tolerance. Because a pure state's complete projective branch vectors are
orthogonal, the measured-system entropy equals the finite-outcome Shannon
entropy. A conservative total variation bound δ≤2e_w+e_w² gives
h₂(δ)+δlog₂29 as the entropy error estimate.

This prototype combines naturally reached charged circulation, a retained host
and gauge-invariant density response in one finite recorded history. The event
is rare, the explicit record is an intervention, the graph diameter is two, and
the result establishes neither a dominant autonomous geometric phase nor a
macroscopic causal cone.

Reproduce with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python neutral_record.py`.
The trusted quotient cache is optional; a complete quotient is rebuilt if the
local cache is unavailable. The physical postrecord Hamiltonian is rebuilt and
the isometry/intertwiner checks are always rerun.
