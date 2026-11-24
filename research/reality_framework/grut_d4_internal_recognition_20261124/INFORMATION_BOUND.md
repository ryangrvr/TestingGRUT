# How much of the past can be restored?

24 November 2026. These bounds concern a classical historical configuration,
with all systems accessible to the repair dynamics included in the input.

## Historical restoration is bounded by state discrimination

Let R be the historical label with prior probabilities p_r. Let rho_r be
the complete state accessible to the repair mechanism at the moment it
responds, including every memory, bath and clock degree it can use. Suppose
a channel E maps that state to the output, where orthogonal projectors P_r
score the restored configurations. The adjoint of E maps the P_r to a
sub-POVM on the accessible input. Assign any remaining outcome to a guess
to complete a POVM. Hence

$$
 P_{\mathrm{restore}}
 =\sum_r p_r\operatorname{Tr}[P_r\mathcal E(\rho_r)]
 \le \max_{\{F_r\}}\sum_r p_r\operatorname{Tr}[F_r\rho_r]
 = P_{\mathrm{guess}}(R\mid\mathrm{accessible}).
$$

In conditional min-entropy notation, the last quantity is
`2^(-H_min(R|accessible))`. The channel may include arbitrarily long
finite evolution and a measurement or a fixed coherent apparatus. It
obeys the same bound. [König, Renner and Schaffner](https://arxiv.org/abs/0807.1338)

For d equiprobable historical labels whose accessible states have common
support of dimension D, use `rho_r <= I_D` inside each POVM expectation.
Summing traces gives `P_restore <= min(1,D/d)`. An independent ancillary
state can be absorbed into the POVM and does not increase D for this
bound. A b-bit record with no other information about the history can
therefore support exact recovery of at most 2^b equiprobable alternatives.
For error probability epsilon, the mutual information of the scored guess
and R is at least `log2(d) - h2(epsilon) - epsilon log2(d-1)` by conditioning
on whether the guess failed. This separates storage of a past value from
storage of the syndrome needed to repair it.

## Two copies detect a mismatch but cannot identify the changed copy

Take two distinct labels a and b. From initial (a,a), change the second
copy to b. From (b,b), change the first copy to a. Both histories give
(a,b). If the accessible ancillas agree too, every repair effect has the
same expectation under the two inputs. The probabilities of returning
(a,a) and (b,b) sum to at most one. The equal-prior success ceiling and the
worst-case ceiling are both one half. A coin guess of the two histories
attains this ceiling on that pair.

This is the obstruction to repairing both unflagged system and memory
errors with one system and one memory. If a bath or clock is correlated
with which copy changed, include that additional information in rho_r
and apply the same discrimination bound. A located erasure carries
such information in its erasure symbol.

For n replicas with substitutions at at most t sites, let f_r(y) count
copies of label r in a current word y. The set of historical candidates is

$$
 C(y)=\{r:f_r(y)\ge n-t\}.
$$

Distinct initial codewords are distance n apart. For d>=2, all admissible
histories are unambiguous exactly when n>2t. At a word consistent with
k histories and identical ancillary records, the minimax ceiling is
1/k. For t<n the largest possible k is
`min(d, floor(n/(n-t)))`. Each of the k labels needs at least n-t copies.
Conversely, one can allocate that many positions to each compatible label
and use any remaining positions to complete the word.

Three replicas suffice for one substitution. They do not distinguish the
current word 100 that came from 000 by changing S from the same word that
came from 111 by changing M1 and M2. D4 has an exact collision even when
its writer spin and both reservoirs are included. For d>=3 and two possible
substitutions, three different current values permit three histories, giving a
one-third worst-word ceiling. The code assumes a bound on simultaneous damage.

## Capacity for the information released by correction

Fix one historical label. There are K=1+n(d-1) orthogonal possibilities
consisting of no error or one changed site and value. Exact unitary repair
maps each logical input to the same pure codeword. The environment outputs
must therefore distinguish those K inputs by orthogonal states. For a fixed
initial environment of rank r_E, the K input subspaces have dimension
K r_E in total. Their image lies in the space of the one codeword times the
environment, so `dim(environment) >= K r_E`.

With uniform error input and an independent initial environment state sigma_E,
exact repair leaves the codeword pure and preserves global entropy. The
environment's increase in entropy is then exactly log2(K). Its unused
capacity must be at least this large. With a pure environment and uniform
input, approximate average codeword probability p also obeys
`p <= min(1,dim(environment)/K)` by the rank of the output codeword projector.
The capacity statement concerns the information released by a correction,
separate from the label needed for historical identity.

For D4's binary three-copy case, the two chains represent the four
possible count signatures 00, 11, 10 and 01. For an isolated fault with
initially empty chains, these match no fault, S, M1 and M2 upon complete
repair. The same chains can add further quanta across later rounds.
Terminal counts only preserve the accumulated parity totals. They need
not identify the exact sequence. There is no assumption that all ordered
sequences are stored in orthogonal reservoir states.

| Historical labels d | Bits needed for perfect history discrimination | One bit alone: maximum success | Orthogonal outputs for no or one fault in three copies |
|---:|---:|---:|---:|
| 2 | 1 | 1 | 4 |
| 3 | 1.585 | 2/3 | 7 |
| 4 | 2 | 1/2 | 10 |
| 8 | 3 | 1/4 | 22 |
| 16 | 4 | 1/8 | 46 |
| 64 | 6 | 1/32 | 190 |

For orientation, a d-label single-fault realization can use the same
equal-label writer and one finite drain for each of the 3(d-1) ordered
site and nonzero displacement pairs (j,k). Let `E_(r,j,k)` be the string
obtained from (r,r,r) by changing site j to r+k modulo d. Exchange
each `E_(r,j,k)` with (r,r,r) while creating an excitation in port (j,k),
summed over r with identical couplings. Global cyclic label shift is a
symmetry. A single substitution selects one charge-one finite XX path.
The finite path proof therefore gives a first-repair guarantee by the
same source survival estimate, with the writer perturbation accounted
for. This extension shows how the number of syndrome ports can grow,
while the general lower bound counts their total distinguishable
environment states. The d>2 case has no registered repeated physical
trial. The scaling table and tuple enumeration test the information
bounds. [RESULTS.json](RESULTS.json)
