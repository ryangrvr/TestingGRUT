# The explicit microscopic law being tested — unchanged S4/WM2

This restates the existing law, not a successor Hamiltonian. No new
coefficients are fitted in S9.
Units set hbar and the receiver hopping energy to one. These are unit
conventions and supplied scales, not derived fundamental constants.

## Registers and preparation

There are N labeled events, with no supplied adjacency between event pairs.
An orthogonal geometry state w=(w_0,...,w_(m-1)) is a simple ordered word,
the reachable encoding of the single-end bond/head/free-event grammar.
All choices of a free next label occur with the same coupling; label order
is not a spatial input. The grammar does, however, supply the path class.

q is an N-bit storage mask, q_(w_0)=1. The receiver is a **supplied**
L=64N^2 mode chain, with exact charge n_R=|q|-1 and no extra photon cutoff.
Its vacuum is part of the declared seed. A neutral probe has vacuum or one
particle at an event p. Probe creation is an explicit later intervention.

Initially w consists only of a root, q is its single bit, and receiver/probe
are vacuum. The original root is uniformly superposed (or equivalently
weighted for the block-diagonal observables). It is never optimized.

## Relative Hamiltonian

Let mask(w) denote the occupied word labels. Include the full diagonal

    D=16(|q|-m)+64|q xor mask(w)|
         +sum_v delta_v(q_v-1_(v=root)) + dGamma(h_L),

where h_L has off-diagonal -1 and zero diagonal. delta_v=bar_delta+e_v,
|bar_delta|<=1/2, |e_v|<=1/4; g in [.9,1.1]. The physical receiver
frequencies are 4+sigma(h_L) in [2,6]; the displayed convention subtracts
the original conserved-charge energy and a constant. The neutral seed has
relative mean energy zero. All storage/receiver energy remains accounted.

Put epsilon=J=1/4. Add all the following terms and their exact adjoints:

* **Assisted append:** for v not in w with q_v=0, append v, set q_v=1,
  create a contact-mode photon, with coefficient -g/N. In occupation
  matrix elements the coefficient is -(g/N)sqrt(n_0+1). It does not depend
  on the probe position. Its reverse pops the end, removes q_v and absorbs
  the contact photon. Every receiver occupation allowed by charge remains.
* **Bare append:** for v not in w, append v without changing q or the
  receiver, with coefficient -epsilon/sqrt(N), but only if the probe is
  at the current head or at v. Its reverse has the same head/parent
  condition. It is zero in probe vacuum. These potentially dangerous
  geometry transitions are retained physically throughout S9.
* **Matter hopping:** with coefficient -J, move the probe across each
  consecutive occupied bond of w. No other fixed event-pair adjacency
  is supplied. All geometry changes can alter these bonds dynamically.

Thus H=D+V_assisted+V_bare+V_matter. This is precisely the state-space law
of the frozen S4 sparse construction, including inconsistent storage states.
The decoupled CP-12 benchmark factor remains a spectator; its original
charged matter is not identified with this neutral probe.

## Size-uniform changing-interaction bound

Each assisted child has a unique parent word. A directed assisted append
has squared source action at most

    (g/N)^2 (N-|q|)(n_0+1)
      <=(g/N)^2 (N-|q|)|q| <=g^2/4.

The Hermitian assisted interaction therefore has norm <=g. A directed bare
append has at most N children of amplitude epsilon/sqrt(N), with unique
parent/probe input for each target, so its norm <=epsilon; including its
adjoint gives <=2epsilon. A path's degree is <=2, hence matter norm <=2J.
Consequently **on the complete physical Hilbert space**

    ||V_assisted+V_bare+V_matter||<=g+2epsilon+2J<=2.1.

The configuration-preserving D may have an extensive norm; all full-H
word-count and comparison bounds use its exact interaction picture.
No gap, deleted transition, limited receiver energy, or nonlocal eigenbasis
is used. Before probe injection, storage stays consistent and n_R=m-1;
the exact formation reduction is equation (4) of FORMATION_ANALYSIS.md.
