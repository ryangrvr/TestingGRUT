# WM2-S1: a single-end assembly Hamiltonian

This is a new, bounded construction within WM2, not an alteration of
CP-12 or WM1. S1 uses no coordinates, Euclidean distance, ordering of
event labels or preferred pair adjacency. It does supply a single-end
assembly grammar, which restricts reachable geometries to paths. That
one-dimensional bias is an input, not a derived dimension of nature.

## Primitive registers and support

For each of N event labels i introduce a head bit h_i, a free bit f_i
and a neutral probe mode c_i. For each ordered pair u!=v introduce a
directed bond bit b_uv. Raising/lowering symbols below are hard-core
qubit operators; c is fermionic. Define

\[
X_{uv}=b_{uv}^{+}h_v^{+}h_u^{-}f_v^{-},\qquad
K=\sum_i f_i,
\]

\[
H_{S1}=UK-\frac{\epsilon}{\sqrt N}
\sum_{u\ne v}(X_{uv}+X_{uv}^{\dagger})
-J\sum_{u\ne v} b_{uv}
(c_u^\dagger c_v+c_v^\dagger c_u).
\]

In the last sum b_uv denotes its occupation projector, not its raising
operator. Each assembly term acts on four qubits (one bond, two head
bits, one free bit); a probe term acts on one bond and two probe modes.
The diagonal uses one free bit at a time. Every ordered pair occurs
with the same coefficient. Labels carry no spatial information. The
epsilon/sqrt(N) normalization is an explicit supplied rule: it pays for
the growing number of possible attachments and is not an earned constant.

Append X_uv moves the head from u to a free v and adds u->v; its adjoint
removes the incoming edge of the current head. No degree test involving
all incident edges and no stored spatial order is used in a primitive
term. Number of heads and number of probes are conserved. K plus total
bond number is conserved. All terms are Hermitian in the stated sums.

## Initial state and complete invariant configuration component

Let |seed_r> have no bonds, h_r=1, f_r=0 and f_i=1 for i!=r, with other
head bits zero. Use the permutation-invariant seed

\[
|s_0\rangle=N^{-1/2}\sum_r|seed_r\rangle.
\]

The root sectors are invariant. The calculations use representative
r=0; permutation covariance supplies the others. An unseeded state
with all heads zero and no bonds is dark: the seed is a real additional
preparation resource, explicitly supplied rather than derived.

Induction under the actual register terms proves that the entire
reachable component consists of directed words w=(r,v_2,...,v_m) with
distinct labels, with head at the last entry and every absent label
free. Every partial ordering is included. Append any unused label;
pop the last label; hop a probe across any currently occupied bond.
No other register transition leaves this component. Its orthogonal
complement is invariant because the Hamiltonian is Hermitian. The
tuple representation in model.py is therefore an exact representation
of this component, not a cutoff eliminating virtual geometries.

The graph of allowed configurations has states (w,p), including free
probe positions p not in w. Such a free probe is stationary until its
event attaches. The graph degree can grow with N during assembly; only
the protected **event** graph below has degree at most two. Finite
primitive support is not confused with bounded configuration degree.

## Energy-defined sector and generated operational metric

P0=1_(K=0) within the reachable component includes **every** fully
assembled ordering. No preferred word is selected. The exact low
spectral projector Pi is separately derived in MATHEMATICS.md.

One primitive transition remaining inside P0 is a probe hop along a
bond of its particular assembled word. There is no primitive P0-to-P0
geometry reconnection. Thus, after a complete geometry record w,
d_op^w(i,j) is the shortest number of these protected primitive hops.
It equals the separation of i and j in w, a derived description of the
occupied bonds rather than coordinates supplied to H. This metric is
fixed before the full-Q resolvent is evaluated. The theorem compares it
to **all** virtual words and **all** final assembled orderings.

These conditional event metrics have degree <=2, diameter N-1 and
linear ball growth. Their adjacency depends on which events attach;
their one-dimensional class follows from the grammar. That distinction
is a standing limitation of this construction.

## Original-neutral embedding, recording and probe preparation

For the requested original neutral preparation, an explicit benchmark
embedding is

\[
|\Psi_{\rm init}\rangle=|\psi_{\rm CP12,neutral}\rangle
\otimes|s_0\rangle\otimes|0_c\rangle,
\qquad H_{\rm embed}=H_{\rm CP12}\otimes I+I\otimes H_{S1}.
\]

The old CP-12 factor remains complete and untouched; no old charged
state is removed. This deliberately minimal product embedding does
**not** establish a common geometry for its charged a,b excitations or
a full-system isolated low-energy band. It permits an exact test of S1
assembly from the original ray with all new registers specified.
There is no reservoir, annealing schedule or hidden energy sink. The
seed has auxiliary energy U(N-1) and energy width epsilon sqrt((N-1)/N).
The resulting formation obstruction is proved for all times.

The later probe controls inject one c fermion after recording a fully
assembled branch. They are conditional transport controls, not a
successful formation experiment. Their source is a local c-number
phase with variance 1/4. Probe number conservation means neither the
probe nor this packet is produced automatically from probe vacuum.
No universal low-energy matter claim follows.

A complete geometry record includes every partial word and every root,
not just successful assemblies. Its energy change is minus the mean
assembly kinetic energy. The full radial entropy and record work are
reported; postselection probabilities are never hidden. This ideal
instrument calculation does not model a physical measurement device.
