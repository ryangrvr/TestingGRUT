# Preparation and records are physical resources

The original neutral graph ray has a simple constructive external preparation;
it does not require a lookup table of desired spatial graphs or edge phases.
Encode M possible graph flags and start a reversible register at r=Q. At stage
j, with M-j remaining flags, perform the count-controlled rotation

\[
|0\rangle\mapsto
\sqrt{1-r/(M-j)}|0\rangle+\sqrt{r/(M-j)}|1\rangle,
\]

and decrement r coherently when the flag is1. At the last stage r=0 for every
branch. Every Q-edge graph has amplitude1/sqrt(binomial(M,Q)). Vertex vacuum,
the filled potential-edge fermion sea and E=0 complete the supplied input and
satisfy Gauss law on every branch. `preparation.py` checks the exact rational
branch probabilities on the finite control cases, including455 graph parents.

This uses M controlled stages and an O(logQ)-bit count register. A conservative
table implementation has O(MQ) count-controlled rotations plus reversible
arithmetic; finite-precision gate synthesis adds its own overhead. The angles
depend only on counts, not topology. This is polynomial circuit-description
complexity, not a demonstration that the frozen Hamiltonian autonomously makes
its initial ray. The controller, initialization, quantum gate model, energy,
precision and elapsed preparation time remain supplied physical resources.

Likewise, the neutral-origin recorded prototype retains every outcome of the
complete {P_gG,P_g(1-G)} channel. There are30 formal outcomes and18 nonzero
outcomes here, requiring a memory with at least18 distinguishable records.
The induced measured-system entropy is4.0824bits; system energy rises by
0.065474 in the frozen units. These are computed costs of the specified channel,
not a complete hardware work/reset budget or a free selection principle.
Sampling the G records at probability0.0550221 takes about18.17 trials per
accepted event on average. Unconditional probabilities must remain reported.

The record is used to define the receiver's **pre-source intrinsic distance**.
It avoids redefining distance after a geometric change, while paying for the
intervention that makes the history operational. No consciousness variable
selects an outcome; a unitary implementation correlates the system with physical
memory, and a subsequent conditional analysis uses the ordinary supplied Born
probabilities. Macroscopic stable internal records generated without this
external diagnostic still need to be demonstrated.
