# The supplied preparation already contains a giant component

The diagonal graph distribution of the coherent initial state is exactly a
uniform Q-edge graph on N-1 vertices, together with the isolated reference.
Coherence changes subsequent dynamics, but not this initial graph measurement.

The classical giant-component result of Erdős and Rényi applies when
Q/(N-1) -> nu > 1/2. Its limiting largest-component fraction is the nonzero
solution S of S=1-exp(-2nu S). At nu=1 this is approximately .796812.
The isolated-vertex fraction tends to exp(-2nu), about .135335 at nu=1.
Thus the same preparation can have a macroscopic connected component while
the probability that **every** vertex is connected tends to zero.

Primary source: Erdős and Rényi, *On the Evolution of Random Graphs* (1960),
Theorem 9b and the following discussion, printed pp.56–57; degree probabilities
in §10. [Original paper](https://static.renyi.hu/~p_erdos/1960-10.pdf).

An independent sampling control of this exact initial distribution gives:

| N | Mean largest-component fraction | Mean sampled graph distance within largest component | Mean maximum degree |
|---:|---:|---:|---:|
| 64 | .79656 | 4.90 | 5.79 |
| 256 | .79559 | 7.08 | 7.00 |
| 1024 | .79609 | 8.96 | 7.97 |
| 4096 | .79625 | 11.07 | 8.72 |

One hundred independent graphs per size, Q=N-1, seed12012026. Mean distance uses
eight sampled source vertices per graph, not all pairs. Standard errors and
sampled lower bounds on diameter are in `preparation_geometry_results.json`.
This is sampling of the initial quantum state's diagonal probabilities,
**not** a simulation or approximation to full quantum time evolution.

The partner's global-connectivity theorem remains valid. Its physical target
must be distinguished from a spatial macroscopic component: demanding one
component containing absolutely everything is stronger than needed for a
large connected region in which observers might reside. Conversely, observing
a giant component here would not establish formation, since it was already
typical of the supplied preparation.

The needed target is geometric and operational: a macroscopic component with
controlled volume-versus-radius growth, excitation/transport scaling, and
independently tested detector locality. The sampling table does not establish
any particular continuum geometry or dimension. The full-dynamics packet's
path/ring controls are finite-size diagnostics only; neither they nor global
connectedness can substitute for that target.
