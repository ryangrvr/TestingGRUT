# WM2-S3 microscopic law and declared resources

S3 is a separate finite closed Hamiltonian. It preserves the single-end
grammar and replaces the S2 storage-assisted term. It changes both
storage energy and append normalization openly. No frozen model is edited.
Use hbar=1 and the supplied model energy units.

## Registers and Hamiltonian

For N events retain bond bits b_uv, one head bit h_i per event, free bits
f_i, storage bits q_i and the neutral probe c_i. Define

\[
X_{uv}=b_{uv}^{+}h_v^{+}h_u^{-}f_v^{-},\quad
C_i=f_i+q_i-1,\quad M_{\rm mis}=\sum_i C_i^2,\quad
F_{uv}=n_{c,u}+n_{c,v}.
\]

Add a shared receiver of L bosonic modes a_r, initially vacuum:

\[
H_R=\Delta n_R-\sum_{r=1}^{L-1}
(a_{r+1}^{\dagger}a_r+a_r^{\dagger}a_{r+1}),\qquad
n_R=\sum_r a_r^\dagger a_r,\quad \Delta=4.
\]

Its one-quantum spectrum lies in [2,6]. The complete law is

\[
\begin{aligned}
H_{\delta,g}={}&U\sum_i f_i+(U-\Delta+\delta)\sum_i q_i
 +\Lambda M_{\rm mis}+H_R\\
&-\frac gN\sum_{u\ne v}
(X_{uv}q_v^+a_1^\dagger+a_1q_v^-X_{uv}^\dagger)\\
&-\frac{\epsilon}{\sqrt N}\sum_{u\ne v}
F_{uv}(X_{uv}+X_{uv}^\dagger)\\
&-J\sum_{u\ne v}b_{uv}
(c_u^\dagger c_v+c_v^\dagger c_u).
\end{aligned}
\]

Values: U=16, Lambda=64, epsilon=J=1/4, delta in [-1/2,1/2],
g in [9/10,11/10]. All ordered event pairs have identical coefficients.
A primitive append acts on the original four geometry bits, q_v and the
contact oscillator. Reverse terms are retained. Receiver mode motion
acts on two oscillator registers. Probe-assisted bare append and matter
hopping retain their S2 supports.

The bath spectrum and contact vector are inputs. The Jacobi mode index
is internal receiver structure, with no assignment to event labels or
event coordinates. The same receiver couples only through the movable
end. A unitary mode-basis change can remove the receiver's chain form
without changing event-pair symmetry or the results. No event distance
is read from that internal mode chain.

The receiver is nevertheless a **shared additional resource**, not a
derived local field. Causality is proved for the declared neutral probe
and grammar; arbitrary new couplings or species could invalidate it.

## Initial preparation and exact excitation capacity

For each root r use no bonds, h_r=1, f_r=0, q_r=1, all other h=0,
f=1,q=0, and receiver/probe vacuum. Superpose roots uniformly.
No formed geometry or incoming reservoir packet is supplied.

The law conserves sum q - n_R=1 on this component. Consequently

\[
0\le n_R=\sum_iq_i-1\le N-1.
\]

This finite excitation capacity follows from the exact charge, including
after probe injection. It is not a boson-number truncation. All allowed
storage patterns and all occupations of all L modes are retained.

The physical implementation in the theorem uses **L_N=64 N^2 modes**.
Its initial receiver energy is zero and its maximum possible receiver
energy on the component is at most 6(N-1). Mode count grows quadratically;
an infinite receiver is only an analytic comparison with an explicit
finite-model error bound. There is no resetting or irreversible sink.

The original CP-12 neutral ray remains a separate product benchmark
under H_CP12 tensor I + I tensor H_S3. No original term is removed and
no energy flows from that factor into S3. Common geometry for its
original charged species, a combined low-energy band and a continuum
gravitational law are not established.

## Energy ledger

The initial auxiliary energy is UN-Delta+delta. It includes free-event
energy U(N-1) and the supplied root-storage energy U-Delta+delta.

| Process | Free-event change | Storage change | Receiver change |
|---|---:|---:|---|
| Emission-assisted append | -U | U-Delta+delta | Create a contact quantum |
| Reverse assisted pop | U | -(U-Delta+delta) | Annihilate a contact quantum |
| Bare append at consistency | -U | 0 | Unchanged; mismatch penalty +Lambda |
| Bare pop at consistency | U | 0 | Unchanged; mismatch penalty +Lambda |

The created contact quantum is not an energy eigenstate. The receiver's
actual positive spectrum provides frequencies throughout [2,6], so it
can accept the released Delta-delta in [3.5,4.5]. No coefficient is
adjusted to delta. The internal hopping supplies the full spectral
dynamics and all reabsorption paths.

In the consistent vacuum component, let k=m-1 be net attachments.
After subtracting the seed energy, the exact Hamiltonian has diagonal
delta k plus the receiver hopping and assisted interaction. Thus

\[
\delta\langle k\rangle+\langle H_R-\Delta n_R\rangle
+\langle H_{\rm int}\rangle=0.
\]

Equivalently, initial binding energy equals remaining binding energy,
new storage energy, physical receiver energy and interaction energy.
MATHEMATICS.md and the finite-model controls account for all four.

## Reference sector and metric

P0 denotes **storage consistency**, all C_i=0, at every assembly size,
with every allowed receiver state. In probe vacuum it is exactly
invariant. After injection the bare interaction can leave it.
P0 must not be described as an isolated exact low-energy band in S3.

A recorded geometry is a rooted distinct-label word. Operational event
distance is the minimum protected primitive operations, allowing real
assembly, reverse assembly and receiver motion, between its attached
events. The prefix proof gives their path separation; all receiver
states are included. The one-dimensional class is still grammar-supplied.

## Frozen stability class and remaining boundaries

The positive interval covers uniform permutation-symmetric storage
detuning and a common emission coefficient in the intervals above.
It does not cover independent site disorder, a changed receiver
spectral measure, arbitrary added interactions or finite probe density.

The field vacuum, bosonic statistics, receiver spectrum, growing mode
count, normalization, root seed and path grammar are inputs. The finite
model is reversible and eventually permits recurrences. The theorem is
a scalable observation-window result, followed by conditional bulk
persistence, not an irreversible equilibrium phase. No empirical
prediction or fundamental constant is derived.
