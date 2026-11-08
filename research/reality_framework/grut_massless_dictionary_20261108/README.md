# B1-01 — measured massless histories improve; dictionary control still misses

The unfitted finite-volume detector reduces the N20 maximum discrepancy from
the ideal smooth massless target from **0.00167772 to 0.000420236**. First-readout,
sequential-history and disturbance comparisons meet their frozen 0.001 windows.
However the independently calculated cell-detector-to-smooth-target difference
is **0.00101075**, above its 0.00025 allowance. Other discrepancies partly cancel
it. The full registered test therefore remains **massless dictionary matching
unresolved**, and the corrected massive experiment was **not evaluated**.

This improves a declared experiment without fitting a physical parameter. It
does not establish a controlled detector dictionary for transfer, a continuum
certificate, or a new GRUT law. P1's interacting remainder stays separate.

## Derivation and chronology

[DERIVATION.md](DERIVATION.md) and [PROTOCOL_B1_01.md](PROTOCOL_B1_01.md), including
the code, analysis and [parameter card](protocol_b1_01.json), were frozen at
[`871a14c41a8afcdd678ffffc125d51aa712d2efa`](https://github.com/ryangrvr/TestingGRUT/commit/871a14c41a8afcdd678ffffc125d51aa712d2efa)
before corrected physical histories. Production used a clean checkout of that
public commit. Earlier M0/P1 data were already known: this is a new registered
dictionary hypothesis, not an untouched holdout or a retroactive repair of M0.

The exact Gauss-law identity is L_j=S_j+1/4+(-1)^j/4. With link cells I_j and
the original normalized Gaussian f, set v_j=integral_Ij f and b=1-sum_j v_j.
The rule is

\[
F_a^{lat}=\sum_jv_j[L_j-(-1)^j/4]+b/4=1/4+\sum_jv_jS_j.
\]

It integrates the prescribed smearing, removes the explicit alternating
c-number, and retains the boundary contribution. Interior weights are not
renormalized. Length, mass prescription, sigma, detector strength, times,
Hamiltonian and ground-state preparation rule remain fixed. The matrix W=exp(2iF)
is unitary on the finite lattice, and K±=(I±W)/2 implements the full disturbance.

The cell-average interpretation is an explicit finite-resolution assumption,
not a unique consequence of the spin Hamiltonian. The free hopping modes derive
the ghost endpoint convention; an exact interacting finite-spacing boundary
matching remains unproved. A smooth-classical-profile reconstruction bound is
O(a^2), with its norms and domain stated. It is **not** applied to unsmeared quantum
vacuum fluctuations. The latter are compared directly using two analytic Gaussian
references: the actual cell smearing with boundary completion, and the original
smooth target. No coefficient was selected to minimize their discrepancy.

## Frozen comparison results

All numbers below are complete outcome-distribution total variation, except the
disturbance comparison, which is half the L1 difference of signed marginal
responses to the earlier readout.

| N20 quantity | Result | Frozen allowance | Verdict |
|---|---:|---:|---|
| First readout versus smooth target | 0.00000801062 | 0.001 | pass |
| Maximum sequential history versus smooth target | 0.000420236 | 0.001 | pass |
| Lattice versus matched cell reference | 0.000832336 | 0.001 | pass |
| Cell reference versus smooth target | 0.001010753 | 0.00025 | **miss** |
| Signed disturbance difference | 0.000428247 | 0.001 | pass |
| N18→N20 history change | 0.000104770 | 0.00025 | pass |
| Refined N14 solver change | 4.68e-15 | 1e-8 | pass |

| Times | Lattice→smooth target | Lattice→cell reference | Cell→smooth target |
|---|---:|---:|---:|
| (0) | 0.000008011 | 0.000249779 | 0.000241768 |
| (0,1) | 0.000420236 | 0.000295986 | 0.000605547 |
| (0,2) | 0.000295638 | 0.000832336 | 0.000540604 |
| (0,1,2) | 0.000420236 | 0.000811200 | 0.001010753 |

The gate miss is not a failure of floating-point precision or a resource stop.
The analytic reference's maximum registered tail bound is 5.50e-7, much smaller
than the geometry miss. Cell-reference mode doubling changes histories by at
most 8.81e-11 at N20; target-reference mode doubling by 7.76e-14. These reference
controls are ordinary floating-point evaluations, not interval enclosures.

The N10/N14/N18/N20 target discrepancies are 0.00203565, 0.000902741,
0.000522558 and 0.000420236. Their corresponding cell-to-target differences are
0.00369392, 0.00198518, 0.00123584 and 0.00101075. Both sequences decrease in the
tested range. No convergence exponent or infinite-resolution remainder is fitted.

## Why the apparent target success is insufficient for transfer

For the longest history, the two comparison differences have half-L1 sizes
0.000811200 and 0.001010753, but their vector sum has size only 0.000420236.
Their directions partly cancel. This is a derived statement about the saved
probability vectors, not a claim that a fixed percentage of error has one cause.
The separate geometry gate was frozen specifically to require control of this
component rather than acceptance through an aggregate cancellation.

The static field mean at N20 is 0.12731175 and variance 0.06215489, compared with
the ideal continuum values 0.12722990 and 0.06216349. These close static values
do not by themselves establish sequential response. The first-readout TV is
about 52 times smaller than the largest sequential TV.

The preregistered four-part signed diagnosis uses two *auxiliary* Gaussian
summaries, one retaining the lattice's measured two-point kernel and another
retaining only its static moments with continuum time dependence. It finds,
for (0,1,2):

| Signed contrast | Half-L1 magnitude |
|---|---:|
| Beyond the selected two-point description | 0.000487702 |
| Temporal two-point response | 0.000135710 |
| Static moments relative to the cell reference | 0.001115475 |
| Cell-to-smooth spatial reconstruction | 0.001010753 |
| Sum as a signed probability vector | **0.000420236** |

The contrasts close exactly at stored arithmetic precision. Their norms cannot
be added as causal percentages. The auxiliaries use measured moments solely
for diagnosis; neither is the physical candidate, a fitted correction, or a
massive prediction. The finite-spacing lattice's departure from this Gaussian
summary is compatible with a Gaussian massless continuum limit; it does not
demonstrate interacting continuum non-Gaussianity.

The unresolved term that is directly identified is the chosen detector's finite
spatial reconstruction, including its edge treatment. The remaining
lattice-to-cell difference still combines state, propagation and operator
matching effects. This test does not uniquely identify a missing microscopic
boundary operator or mass correction.

## Costs, verification and scope

N remains at or below 20; no interacting parent cutoff was increased. The five
massless jobs completed in about 61 seconds including subprocess overhead,
within the 600-second total and 180-second per-point caps, using one BLAS thread.
The conditional-transfer barrier recorded
`not_evaluated_massless_gates_not_passed`; there are no corrected massive output
arrays. The earlier raw instrument was replayed at every N, agreeing with its
archive within TV 5.84e-15. That checks the unchanged baseline; it does not make
this a separately implemented physical replication.

Actually checked: analytic cell coefficients against independent quadrature,
free hopping gap against a one-body diagonalization, synthetic ordered Gaussian
words, exact spin/Gauss detector rewrite, clean frozen-checkout production,
solver refinement, complete probabilities, checkpoint dimensions and norms,
16 frozen-source/input hashes, all frozen verdicts, and signed decomposition.
The 120 saved distributions include physical, reference, replay and diagnostic
groups; this is not 120 independent physical tests. Their largest normalization
discrepancy is 1.51e-13. A forged acceptance flag, forged massive-stage verdict and
bypassed transfer barrier are rejected. See [VALIDATION.json](VALIDATION.json),
[raw analysis](raw/analysis.json), [events](raw/events.json),
[provenance](raw/provenance.json) and [transfer decision](raw/massive_stage.json).

## Reproduction

From the repository root of this result checkout:

```bash
python3 -m venv /tmp/grut-b1-venv
/tmp/grut-b1-venv/bin/python -m pip install -r research/reality_framework/grut_detector_matching_20261107/requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /tmp/grut-b1-venv/bin/python research/reality_framework/grut_massless_dictionary_20261108/verify_results.py --negative-controls
```

To repeat production, use fresh paths. The output directory must be outside the
clean frozen checkout:

```bash
git worktree add --detach /tmp/grut-b1-frozen 871a14c41a8afcdd678ffffc125d51aa712d2efa
cd /tmp/grut-b1-frozen
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /tmp/grut-b1-venv/bin/python research/reality_framework/grut_massless_dictionary_20261108/run.py --output-dir /tmp/grut-b1-replay
```

The postproduction verifier is in the result commit, not the earlier freeze.
Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0 were used. Raw probability files and
ground-state checkpoints are included; no older packet is overwritten.

The bounded test stops with its geometry miss intact. No alternate smearing,
coefficient adjustment, endpoint fit or transfer to the massive model was made.
P1's interacting parent remainder is still uncontrolled at the requested
accuracy; E2's distinct finite-parent success remains intact. No new GRUT law,
interacting continuum certification, merge, release or Zenodo upload follows.
