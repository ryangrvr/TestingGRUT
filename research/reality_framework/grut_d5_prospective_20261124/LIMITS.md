# Prospective indistinguishability and a matched null

## Information required before completion

Let Y label a proposed future target and let rho_y give the complete
conditional physical state on all degrees of freedom available to the
development mechanism before the action. This includes hidden records,
seed, substrate, archives, coupling to controllers and any clock or
reservoir carrying information. Let E be the common channel from that
state to the output, allowing any subsequent dynamics that receives
the same further inputs. The projectors Q_y for mutually exclusive
whole-word successes have sum at most identity. Thus

$$
P_{\mathrm{success}}
=\sum_y \pi_y\operatorname{Tr}(Q_y\mathcal E(\rho_y))
\le P_{\mathrm{guess}}(Y\mid\mathrm{accessible}).
$$

The pulled-back effects E*(Q_y) can be completed to a POVM by assigning
the failure effect to a guess. This proves the inequality by optimal
state discrimination. For identical conditional states, its bound is
the largest prior. For d equal-prior states supported within dimension
D, positivity and rho_y <= I on that support give a bound of min(1,D/d).
The connection of conditional min-entropy and guessing probability is
standard. [König, Renner and Schaffner](https://arxiv.org/abs/0807.1338)

If the future choice is independent and unavailable, exact success is
at most 1/4 for four equally likely distinct words. With known phase
and two equally likely rule choices it is at most 1/2. Using a per-cell
reward instead allows 3/4 when phase is fixed, since uniform and
alternating agree at half the cells. For overlapping reward effects
the correct independent-choice baseline is the largest eigenvalue of
the prior-weighted average effect. This is why D5 scores whole words.

A later trace-preserving map cannot signal an independent choice back
to the earlier read. If I_a is the earlier outcome instrument and E_y
is the later map, Tr[E_y(I_a(rho))] equals Tr[I_a(rho)] for every y.
A later observed outcome can be statistically correlated with the
earlier one. Conditioning on that outcome is a different experiment
and its selection probability matters.

D5's target label is computed from the prepared seed and acquired
record. It is available through existing physical correlations. The
disturbed visible tape by itself is an incomplete description of the
state. The complete head still has h, and a disabled register's backup
holds its old value. One cannot infer an information-theoretic erasure
or identical complete states from the record controls. They identify a
causal dependence of the supplied write rule on a particular memory.

## Null that reproduces shape and record acquisition

In the null, change the writer's argument from M to h and leave
everything else in the model specified in MODEL.md fixed. This null
has its own explicit microscopic Hamiltonian. It is the same exchange
formula, with the same clock envelope and with the amended local
color selection on each depositing edge. It still acquires M=h at
the assay. Intact data states follow the same ordered path at each
head position, so the entire intact distribution, including clock
coherence and correlations with M, agrees. A physical A-only
disturbance leaves the two write inputs equal in the intact arms.

Inverting or disabling the acquired M changes the main model's final
word and leaves the null's word intact. Scrambling M gives the main
model half its historical whole-word score and the null full score.
Applying the corresponding control to G affects both models. The
null demonstrates that even an organized pattern with an acquired,
accurate record is insufficient evidence that the record has a
functional role. The perturbation differentiates the two mechanisms.

## Where the program resides

The starting tape contains no colored cells or sitewise copy, and
every deposition uses the same phase and rule selection. The two-bit
choice of member of the family is in h,G. The local color formula,
the prearranged size and geometry, the forward and reverse visits,
their overwrite slots, and the engineered transport weights are all
supplied. There is no rule for selecting this family from macroscopic
dynamics independent of those inputs. The full overwrite acts on
undamaged cells too. This implementation repairs by regeneration
during a scheduled sweep, without detecting a discrepancy as D4 did.

At an internal clock coordinate before the last rebuild the target
has yet to occupy all cells on the ordered legal path. The actual
unobserved quantum head is a superposition of such coordinates.
Coordinate snapshots are compatible with reversals and cannot be
interpreted as observed monotone trajectories. At the tested early
laboratory times the unfinished branches have probability near one,
with bounded nonzero tails. The model does not establish strict
absence of any completed amplitude at all positive earlier times.

The quantum walk realizes a finite local program with no runtime
gate sequencer, while the stage route and diagnostic flags are
physical inputs selected by the experimenter. External unflagged
flips act directly on a tape cell at specified laboratory times.
No head measurement selects fault trials. The complete word and the
authenticated archive are scored jointly with head arrival.
