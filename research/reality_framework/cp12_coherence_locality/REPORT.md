# CP-12: a controlled prepared local-response regime

**New positive result:** frozen CP-12 supports a controlled regime in which
coherent geometry preparations, physical Gauss constraints, and local
matter–gauge detector response coexist. On a supplied three-dimensional
bounded-degree host, the full-response error tends to zero over observation
times growing with the host's linear size.

This is an additive continuation of
[`d2fbd05`](https://github.com/ryangrvr/TestingGRUT/commit/d2fbd05c420a209faae618c88fbe4d1902bc999c).
The ground-state exclusion did not end the investigation. It showed that
ground-phase selection and the survival of a prepared physical geometry
are different questions. Here we derive and test the latter, and quantify
why recording or calibrating geometry changes the former.

The complete proof is in [DERIVATION.md](DERIVATION.md). CP-12's Hamiltonian,
all couplings, its physical constraints, prior candidates, canonical
GRUT-RAI, PR #12, and Stage-3 Card 3 remain unchanged. No CP-13 is opened.

## 1. Full dynamics approaches the local detector theory

Let \(F\) be an independently specified neighborhood pool with
\(S=Q+r_N\) possible short bonds. Prepare any physical state with all its
\(Q\) occupied bonds inside \(F\). The state can be coherent across many
different graphs. The pool's \(r_N\) vacancies need not be fixed or
classically recorded.

For graph-diagonal physical probes with detector norm at most one and
source norm at most one on the prepared sector,

\[
\boxed{|R_{\rm full}(T)-R_{\rm graph\ blocks}(T)|
\le4\frac{\Gamma^2Q(M-S)(r_N+1)}{M^2}T^2
+4\frac{\Gamma Qr_N}{M}|T|.}
\]

The graph blocks retain the **actual rewiring diagonal**. This is not a
comparison against a more convenient modified Hamiltonian. The estimate
covers the entire physical matter/electric space and arbitrary initial
coherence within the declared pool.

With fixed \(Q/N\), the error is
\(O((r_N+1)T^2/N+r_N|T|/N)\). If the pool has bounded incidence degree,
each block has bounded local hopping interactions and on-site electric
terms. Even, gauge-invariant probes satisfy a Lieb-Robinson-type
commutator cone in that independently declared incidence metric. The
full response inherits this cone with the stated error.

Unbounded on-site electric terms and fermion anticommutation are explicitly
handled in the derivation. The mathematical locality argument is established
physics, with primary antecedents in
[Nachtergaele–Sims, 1410.8174](https://arxiv.org/abs/1410.8174) and
[Nachtergaele–Sims–Young, 1705.08553](https://arxiv.org/abs/1705.08553).
The new application is the controlled comparison to CP-12's full evolving
geometry. A Lieb-Robinson bound is not a Lorentzian light cone.

### An actual source and detector, not a geometry-role oracle

The source modulates **all** frozen hopping terms at one vertex; the
detector is a vertex fermion occupation. Both preserve Gauss's law.
The source does not query pool membership or select its occupied
neighbors: the Hamiltonian's presence factors activate existing bonds.
Only its scalar strength is calibrated to bound the prepared-sector
norm. The pool and detector separation are declared before evolution.
No global minimum/maximum-role query is supplied. The source's larger
norm on other graph sectors is recorded and does not enter this
initially pool-supported linear-response estimate.

Complete 1,032-state full-law controls compare full evolution, pool
compression, and complete independent graph blocks. For the coherent
four-cycle pool, at \(T=0.4\):

| Detector's declared incidence distance from source support | Full response | Graph-block response |
| ---: | ---: | ---: |
| 1 | 0.0279111093 | 0.0279388028 |
| 3 | 0.0000148752 | 0.0000146160 |

Independent finite source-pulse derivatives agree within
\(1.6\times10^{-10}\). The analytical error bounds hold in both the cycle
pool and a branched-pool control. These finite readings demonstrate the
specified model operation; they are not empirical measurements.

An additional operator-level control diagonalizes complete conditional
path/cycle blocks. Integer nested commutators vanish exactly through every
order shorter than the interaction chain connecting the probes. The
commutator norms satisfy the independently derived chain bound at
cutoffs 3 and 4. Uniformity of the bound does not require the
cutoff-dependent operator norms themselves to converge uniformly.

## 2. The controlled window can be macroscopic

A supplied periodic diamond-type host has \(N=2L^3\) vertices, \(Q=2N\)
bonds, and bounded degree. Additional disjoint short bonds form a pool
with \(r_N\) coherent vacancies and maximum degree five. Its spatial
growth is supplied, not selected by CP-12. Sparse Laplacian spectra
independently agree with the host's exact Bloch formula.

For observation time \(T=L\), the response bound tends to zero if
\(r_N=o(L)\). This includes an increasing vacancy inventory
\(r_N=\lfloor\sqrt L\rfloor\). Thus the result extends beyond one fixed
graph and beyond fixed observation times.

For the illustrative single-vacancy preparation:

| \(L\) | Events \(N\) | Observation time | Absolute full-response error upper bound |
| ---: | ---: | ---: | ---: |
| 16 | 8,192 | 16 | 0.013122 |
| 64 | 524,288 | 64 | 0.002696 |
| 256 | 33,554,432 | 256 | 0.000638 |

These are **analytical bounds evaluated at large sizes**, not full
many-body simulations and not relative percentage errors. No physical
speed is fitted. The microscopic external time remains an input.

The neutral uniform-pool preparation is physically allowed and has

\[
\langle H\rangle=\kappa Q+\Gamma Q(M-S)/M,
\qquad\operatorname{Var}(H)=2t^2Q+b_0^2.
\]

Its energy density approaches 1.2 in the frozen units, and its relative
energy width tends to zero. Complete small-system matrices independently
confirm both moments. It is a finite-energy-density excited preparation,
not a tailored clock history or a claimed ground state. Its neighborhood
and its preparation mechanism remain supplied physical inputs.

## 3. Prepared geometry and ground geometry behave differently

For any initially pool-supported physical state, including coherent
states, full CP-12 obeys

\[
\boxed{\langle k(T)\rangle\le
\Gamma^2 Q(M-S)S\,T^2/M^2,}
\]

where \(k\) counts outside bonds. At fixed extensive \(Q,S\), the
coefficient is bounded independently of \(N\). A macroscopic density of
outside bonds therefore requires at least order \(\sqrt N\) time.
This density theorem alone does not prevent individual shortcuts; the
source–detector estimate above is the separate operational result.

That same estimate places a limit on back-reaction in the recovered
window. Every graph-only readout has exactly zero response under the
independent graph blocks. Its full response to a bounded graph-diagonal
matter intervention is therefore at most the same error bound, and
vanishes in the controlled limit. Direct connection-presence readouts
pass this check. Stable matter propagation and a finite gravitational
metric response have not been derived together.

The ground-state theorem gives the contrasting constraint. At \(Q=2N\),
the graph-changing energy contribution must satisfy

\[
\liminf\frac{-\langle V\rangle}{\|V\|}\ge0.2815756658\ldots .
\]

For any independently specified neighborhood pool of \(O(N)\) bonds,
the ground state's mean fraction of outside bonds is at least this
same asymptotic lower bound. At finite \(N=32,768\), an outward
certificate for a pool of \(4N\) bonds gives
\(\langle k\rangle\ge17,520.77\); the bound is about 26.73% of its
\(Q=65,536\) bonds. The asymptotic lower limit is about 28.16%.

This fixed-calibration theorem is not a universal exclusion of spatial
ground phases. A metric/reference system that moves with and interacts
with geometry is outside its assumptions. Defining the neighborhood
after observing a graph also changes the question. The supplied-host
prepared regime does not establish that missing dynamical reference.

## 4. Recording geometry is a physical operation with back-action

An ideal retained full-graph record dephases only the graph configuration.
Its system-energy increase is exactly

\[
\Delta E_{\rm record}=-\langle V\rangle.
\]

The preceding full-H bounds require an asymptotic ground-state cost of
at least 0.056315 per vertex at \(Q=2N\). For independently encoded edge
records with distinguishability \(D\), the tested architecture gives
\(\Delta E=D^2[-\langle V\rangle]\). An explicit joint system–reference
state calculation confirms it.

| Complete finite ground control | Full graph-record system-energy increase |
| --- | ---: |
| \(N=4,Q=3,L=3\) | 0.1202468379 |
| \(N=5,Q=5,L=3\) | 0.2028193065 |

Separately recording each potential edge obeys the exact sum rule
\(\sum_e\Delta E_e=2\Delta E_{\rm full\ record}\). A single microscopic
record can be cheap even when a complete macroscopic geometry record
has extensive cost. This is not a heat-dissipation bound, a consciousness
effect, or a complete autonomous apparatus. The reference Hamiltonian
and its energy exchange must still be modeled.

## 5. Required irregularity does not automatically destroy organization

The complete five-vertex ground control is connected with probability
0.949036509, compared with the exact uniform graph baseline
\(37/42=0.880952381\). Its degree variance is 0.536714, below the
uniform baseline \(2/3\). Matter and gauge dynamics therefore produce
measurable finite-system organization even though exact regularity is
not predominant.

This is the separate \(Q=N\) finite control, not the \(Q=2N\) asymptotic
theorem. Cutoffs 2 and 3 change connectedness by
\(3.65\times10^{-10}\), coherence energy by \(1.02\times10^{-10}\), and
the reported geometry statistics by at most \(1.63\times10^{-9}\).
These checks do not justify a thermodynamic dimensional extrapolation.

## Reproduction and research boundary

From the repository root, with Python, NumPy, SciPy, and NetworkX:

```bash
python research/reality_framework/cp12_coherence_locality/verify.py
python research/reality_framework/cp12_coherence_locality/locality_pool.py
python research/reality_framework/cp12_coherence_locality/prepared_response.py
python research/reality_framework/cp12_coherence_locality/light_cone_control.py
python research/reality_framework/cp12_coherence_locality/three_dimensional_host.py
python research/reality_framework/cp12_coherence_locality/full_ground_records.py --n 4 --q 3 --L 3
python research/reality_framework/cp12_coherence_locality/full_ground_records.py --n 5 --q 5 --L 2
python research/reality_framework/cp12_coherence_locality/full_ground_records.py --n 5 --q 5 --L 3
```

Set BLAS threads to one for predictable resource use. The last command
retains all 427,986 physical states. All scripts were executed; the
manifest records source hashes and the complete parameter sets.

We have earned a working **conditional local matter–gauge response sector**
inside the unchanged evolving-geometry model. The main scientific
opening is a long-lived excited regime, rather than assuming physical
space must be the absolute ground phase.

We have not earned host formation, physical reference dynamics, Lorentz
invariance, a gravitational constraint, universal spin-2 coupling, or a
complete theory of reality. A Lieb-Robinson cone exists on expanders too;
it does not select spatial growth. The next mechanism must produce the
prepared spatial/reference regime through physical dynamics and then
test its metric response, without supplying the desired geometry as an
answer.

The addressed vertices and source strength are calibrated microscopic
probes. A gauge interpretation of event relabelling would additionally
require relational matter tags/dressing; their physical implementation
is not provided by the labelled response calculation.
