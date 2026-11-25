# D5 exact development and pulse theorem

25 November 2026. This proof refers to the frozen model and its scored
initial states. Each conditional initial state is a computational-basis
seed. The fair mixtures in the protocol follow by linearity.

## Physical sector and transport

On route sites j and j+1 an exchange takes the block with a head at
j and vacuum at j+1 to the block with the head at j+1 and vacuum at
j. Its matrix between these blocks is the local unitary U_j.
The adjoint supplies reverse exchange. Set the term to zero on both
vacant and both occupied blocks. Number of heads is conserved. A
work-dependent term acts on the two head sites and the affected tape
site, which includes its three qutrits. Its exchange matrix has norm
J_j because the two active off-diagonal blocks are a unitary and
its adjoint times J_j. Inserting the frozen J_j and using
`(j+1)(L-j) <= (L+1)^2/4` gives J_j<=J. All onsite dimensions,
support sizes and graph degrees have bounds independent of N.

Let V_0=I and V_(j+1)=U_j V_j. In the one-head sector let

$$
D=\sum_{j=0}^{L}|j\rangle\langle j|\otimes V_j.
$$

Conjugating a bare nearest-neighbor hop by D puts
`V_(j+1) V_j^dagger=U_j` on its forward edge. Hence

$$
H=D(\omega J_x\otimes I)D^\dagger,
$$

where J_x is the spin-L/2 matrix in the J_z basis. This is an
identity on the full one-head work sector, not just the intended
path. Representing that clock matrix as the symmetric subspace of
L independent spins gives, for seed s,

$$
|\Psi(t)\rangle=\sum_{j=0}^{L} a_j(t)|j\rangle V_j|s\rangle,
\qquad
a_j(t)=(-i)^j\sqrt{\binom Lj}
\cos^{L-j}(\omega t/2)\sin^j(\omega t/2).
$$

The clock marginal is binomial with p=sin^2(omega t/2). At
T=pi/omega the state is exactly `(-i)^L |L> V_L|s>`. Every
one-head state recurs at 2T up to the common phase (-1)^L.
The eigenvalues are omega(k-L/2). The preparation at j=0 assigns
them weights `binomial(L,k)/2^L` for any work preparation. Thus all
internal arms have the same energy law, mean zero and variance
`omega^2 L/4`, despite different diagnostic selectors and data.
The same spectrum and input law apply to the matched null.
This is the known engineered spin-chain transfer construction.
[Christandl et al.](https://arxiv.org/abs/quant-ph/0309131)

## Local rule, record and terminal word

The assay is a swap of two labels controlled by h. It produces
M=h from the blank without changing the seed. Intact growth uses
this newly acquired bit. After i forward visits, P=i mod 2.
The first k=N/4 growth visits fill only that prefix, and the
first m=N/2 visits still leave the suffix blank. A color flip
leaves the number of occupied cells unchanged. There are no
further depositions during the return path.

Since m is even, the first rebuild begins with P=0. Each
forward rebuild prepares the requested symbol in the previously
blank B_i and swaps it into A_i, storing its old symbol. The
suffix has never been filled in the earlier ordered states, and
the last rebuild visit is the first one at which the whole tape
is filled for intact arms. Before the first reverse visit the
phase is N mod 2. The pre-toggle makes the phase at site N-1
equal to (N-1) mod 2, and each subsequent pre-toggle gives i
mod 2 at site i. Replacing A from the fresh C bank therefore
maintains the same pattern and records the previous A, even for
the alternating word at the turn.

These statements give terminal A_i=M xor (G*(i mod 2)) for
valid G,M at the central edge. A disabled register requests
blanks. On an intact arm the output is the original intended word.
A fault at F=1 or 2 leaves one site wrong while at least
half the tape is still unbuilt. The corresponding B stores the
flipped symbol, and the subsequent complete word is new in the
ordered path. A fault at F=3 corrupts a completed word, whose
wrong symbol is acquired in C during maintenance. Every archived
flip and original word is present at the exact endpoint.

These are ordered-coordinate statements. Unitary evolution is a
coherent walk with backward edges and recurrence, not a classical
history observed to move monotonically. No clock projection is
applied during the internal runs.

## Causal controls at one damaged substrate

The control edge has identical tape action across K. Only the
chosen head field and, on a swap, its backup can change. The
damaged A and both other tapes are identical at that point for
fixed h,g and F=2. The subsequent deposition rule has no access
to the diagnostic metadata K, F, c or the backups. Inverting G
changes half the final colors. Inverting M changes every color.
Disabling either input makes the final tape blank. The fair
scramble is exactly the average of unchanged and inverted
whole-word success. This gives the exact scores 1, 0, 0 and 1/2
for intact, inverted, disabled and scrambled, respectively, for
both inputs. These arms retain the same dimension, prepared
ancillas, entropy and energy law.

The null reads h instead of M. In every intact history M=h at
all write coordinates, so the two models have the same entire
intact orbit. The null's tape and head statistics and acquired
record match there. The null ignores the M changes under
intervention. The observed contrast identifies which stored
variable the interacting dynamics uses. The pilot and backup
records remain physical information. See [LIMITS.md](LIMITS.md).

## A laboratory pulse propagated with the full Hamiltonian

Take F=0 and K=intact and apply one tape permutation X at time
tau. This state need not remain in the original legal history
span. In the D frame the pulse at position j is
`O_j=V_j^dagger X V_j`. On the appropriate idle band, V_j is
constant, giving one common O. Let q be the binomial probability
of this band. Replacing the exact pulse by I tensor O changes
the full state vector by at most `2 sqrt(1-q)`, since the
two coincide on the band and the difference of unitaries has
norm at most two. The remaining evolution preserves this norm
bound. The ideal comparison has the target, flipped archive
and head at the endpoint exactly. The complement of this
joint event then has probability at most `4(1-q)` in the
actual experiment. All pulse outcomes count.

There is also a coherent exact endpoint formula. For
`U_c(t)=exp(-i omega J_x t)`,

$$
\langle L|U_c(T-\tau)|j\rangle
=(-i)^L\,\overline{a_j(\tau)}.
$$

Indeed the full clock transfer maps |0> to (-i)^L|L>
and all states perpendicular to |0> to the perpendicular
complement of |L>. Compose the resulting bra with
U_c(tau)^dagger. Thus the endpoint work vector, dropping the
common phase, is

$$
\sum_j |a_j(\tau)|^2
    V_L V_j^\dagger X V_j|s\rangle.
$$

Collect the weights for identical **complete** work basis
states before taking squared moduli. Uncombined distinct
work states are orthogonal. For a given final work state with
aggregate weight w its unconditional endpoint probability
is w^2. The missed endpoint weight is the sum of pair products
2 w_i w_j over distinct unordered groups. This identity keeps
extremely small failed probabilities visible without a
subtraction from one. All j are included, even those outside
the intended stage. A common successful group of mass at
least q also gives the stronger direct lower bound q^2.

## Strict finite certificate and its scope

The three fixed bands, in inclusive head coordinates, are

$$
[1+k,L/4],\qquad [L/4+1+k,L/2],\qquad
[L/2+1+m+N,3L/4].
$$

For p within 1/64 of 1/8, 3/8 or 5/8, the respective mean
Lp is at least L/10 from either boundary for the frozen
L=256N and N>=8. The binomial Hoeffding bound gives
`1-q <= 2 exp(-L/50)`. From e>8/3 an exact rational
upper bound is `epsilon = 2(3/8)^floor(L/50)`.
At N=8, `4 epsilon < 7.32e-17`, and it decreases
thereafter. These inequalities cover intervals of single
lab pulses as well as the frozen sample times.

With no pulse the late idle interval preserves the complete
output and archives on all positions after the last
maintenance visit. On the first transit and early return,
whenever p>=7/8, the binomial mean is at least L/10
above the start of that interval. It has probability at
least 1-epsilon. The same ideal pulse comparison bounds
output-plus-archive failure by 9 epsilon in that interval
after any single certified pulse, when the read scores data
without requiring the precise head position L. Finite
recurrence eventually removes this persistence.

The early bands give strong unfinished-stage probabilities
at the probe time, not a claim that wave amplitude at later
coordinates is exactly zero at earlier laboratory times.
The internal path proves the order of its writes, and the
pulse theorem proves a joint unconditioned outcome and
archive score in real time.
See [RESOURCES.md](RESOURCES.md) for the storage and energy
costs and [CONTROLS.md](CONTROLS.md) for the tested contrasts.
