# GRUT-WM1: proposed law and input ledger

GRUT-WM1 means **walk mediator 1**. It is a separately named transport
proposal built over the frozen CP-12 geometry/matter Hamiltonian, not a
revision of CP-12 and not a candidate already meeting the spatial-phase
gate.

Introduce two neutral fermion modes per event, L_i and H_i. In a sector
with one new probe, they are a low register and a detuned mediator
register. With the original CP-12 active-edge bit g_e and electric fields,
define the gauge-invariant through-field projector

\[
 F_e=\mathbf1_{E_{ei}E_{ej}<0},\qquad B_e=g_eF_e.
\]

The **complete proposed microscopic law** is

\[
\begin{split}
H_{\rm WM1}={}&H_{\rm CP12}\otimes I+
\Delta\sum_i H_i^\dagger H_i
-g\sum_i(L_i^\dagger H_i+H_i^\dagger L_i)\\
&-v\sum_{e=(i,j)}B_e(H_i^\dagger H_j+H_j^\dagger H_i).
\end{split}
\]

The notation g for the conversion amplitude is distinct from the edge
bit g_e. Here Delta>0 and g,v>0 are new supplied energies. All original
CP-12 hopping and rewiring terms remain. There is no deletion of an
unwanted path. New probes can change register at an event; only their
mediator register traverses an edge.

Each new term is Hermitian, commutes with the original Gauss operators
and preserves the new total fermion number. In the one-probe sector the
added off-diagonal operator has norm at most |g|+|v|(N-1); the detuning
has norm Delta. Thus at each finite N it is a bounded perturbation of
the original Hamiltonian on its original domain. This is not a uniform
gap proof for the full coupled geometry.

## Formation identity and its limits

Every added term annihilates the new-probe vacuum. For every original
state, including the original equal coherent neutral preparation,

\[
e^{-itH_{\rm WM1}}(|\psi_0\rangle\otimes|0\rangle_{LH})
=e^{-itH_{\rm CP12}}|\psi_0\rangle\otimes|0\rangle_{LH}.
\]

This preserves **exactly the old formation probabilities, successes and
failures in probe vacuum**. It does not improve them. It does not create
a low-energy probe from that vacuum: number conservation forbids that.
Injecting a probe, recording a geometry and preparing its dressed state
are additional operations. Their cost and effect must be included in a
future complete original-neutral experiment.

The original charged a,b modes still propagate under CP-12. Calling them
"only geometry" would not establish their decoupling or universal matter
locality. That gate remains open.

## The exact solvable module

Conditioning on a **static supplied** through-field graph G with adjacency
A gives, in the one-probe sector,

\[
H_G=\begin{pmatrix}0&-gI\\-gI&\Delta I-vA\end{pmatrix}.
\]

The cycle and square tori in the controls are supplied graphs. Degree
bound z=4, their topology and their distances are inputs. No result on
these fixtures demonstrates graph formation, a selected dimension, a
continuum or gravity. A static-module theorem is kept separate from a
theorem about WM1 with all CP-12 geometry states retained.

## Independent geometry-change control

GC1 is an independent toy Hamiltonian, **not CP-12 and not a hidden
truncation advertised as WM1**. A geometry bit q=0,1 selects a cycle A
or that cycle plus one antipodal chord S. It has four bands, ordered
(0,L),(0,H),(1,L),(1,H):

\[
H_{\rm GC1}=\begin{pmatrix}
0&-gI&-\lambda I&0\\
-gI&\Delta I-vA&0&-\lambda I\\
-\lambda I&0&\delta I&-gI\\
0&-\lambda I&-gI&(\Delta+\delta)I-v(A+S)
\end{pmatrix}.
\]

Distance is fixed to the q=0 cycle **before** evaluating any coupling.
This deliberately stringent control has a uniform gap, bounded degree
(at most three), small geometry excursions and all virtual paths
retained. Its role is to test whether those properties alone protect
the original metric. They do not. It is not a claim that every CP-12
phase has the GC1 dynamics or the same exact all-orders lower bound.

## Input ledger

| Supplied input | Derived consequence in this packet | Still unearned |
|---|---|---|
| Quantum/Fock structure, external time and frozen CP-12 | Probe-vacuum identity, Gauss and number preservation | Quantum mechanics and time from GRUT |
| New neutral L,H species, Delta,g,v and flux gate | WM1 is explicit; an exact activation word survives | Universal protected matter |
| Static graph G, degree bound z | Complete Schur sum, gap, localized dressing and transport bounds | Dynamically selected spatial phase |
| Cycle/torus fixtures and source packet | Full-H density-phase responses with positive variance | Original-neutral source preparation |
| GC1 geometry bit, delta,lambda and chord | Uniform gap plus a rigorously nondecaying virtual shortcut | A replacement law solving geometry protection |

No supplied parameter is an earned quotient value. SCOUT-1 remains frozen.
