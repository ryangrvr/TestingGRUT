# Constructive GRUT research sprint — 15 October 2026

The frozen CP-12 model now supports a tested conditional construction connecting
charged network formation, virtual gauge interactions, persistent relational
records and matter response. It is not a completed theory of spacetime.

Start with [the executive report](EXECUTIVE_REPORT.md) or the full
[working-theory manuscript](MANUSCRIPT.md). A compact independent-review assignment
is in [CHAT_REVIEW.md](CHAT_REVIEW.md).

| Component | Main evidence | Scope |
|---|---|---|
| Charged accessibility | [Derivation](charged_access/DERIVATION.md), [exact outputs](charged_access/RESULTS.json) | Exact leading amplitudes, N5–N7, graph choice and signed interference |
| Induced gauge/coordination interactions | [Derivation](induced_gauge/DERIVATION.md), [outputs](induced_gauge/RESULTS.json) | Controlled heavy-matter limit; energy-dependent frozen-parameter check |
| Locality bridge | [Synthesis](theory_context/synthesis.md), [checks](theory_context/LOCALITY_BRIDGE_CHECK.json) | Exact conditional graph-block bound; supplied bounded-degree hosts |
| Joint protection and response | [Report](protected_geometry/REPORT.md), [outputs](protected_geometry/RESULTS.json) | Full finite dynamics from supplied charged ring; two cutoffs |
| Neutral-origin recorded prototype | [Derivation](protected_geometry/NEUTRAL_RECORD.md), [outputs](protected_geometry/NEUTRAL_RECORD.json) | Explicit recording costs; final event weight retained |
| Preparation resources | [Accounting](PREPARATION.md), [exact checks](PREPARATION.json) | Topology-neutral externally controlled preparation |
| Earlier pending reproduction | [Report](reproduction/REPRODUCTION_REPORT.md), [comparison](reproduction/COMPARISON.json) | Published1801afe, fresh Node/Python, no cache, all65 checks |
| Internal independent review | [Review](independent_review/review.md), [algebra](independent_review/ALGEBRA_CHECK.json) | Separate worker review; not external certification |

## Reproduce

Use the repository at this packet's commit. Python3.12, NumPy, SciPy and Matplotlib
are sufficient for the new work; Node24 is used for the prior exact suite.
The tested package versions are in `requirements.txt`.

```bash
cd research/reality_framework/constructive_sprint_20261015
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py --rerun
```

The command re-executes all eight new analytical/numerical programs in a temporary
mirror, compares their scientific outputs, validates source hashes and checks
saved prior reproduction logs. It does not replace saved scientific outputs.
The root verification record is regenerated. Per-program commands also work,
for example `python charged_access/run.py` and
`python protected_geometry/verify.py --rerun`; these replace their own outputs.

Add `--baseline` to the root command to repeat the older Node/Python suite as
well. That suite completed during this sprint and took about12 minutes including
both runtimes. The new exact N7 calculation is roughly100 seconds on the tested
single-thread runtime; the complete new suite takes several minutes. Calculations
run sequentially to avoid overlapping their workspaces. Platform timings vary.
No formation cache is needed.

All additions are confined to this directory. The parent is1801afee27096c03c2a3b64cc811ac0be0d2a4a7,
continuing the requested a2f447dc37ea0a71f62c6e990a4630d727ebd85e baseline.
The [manifest](MANIFEST.json) records packet hashes and source provenance.

## Main physical limit

No result establishes a dominant autonomous three-dimensional phase, emergent
external time, a Lorentz-invariant continuum, Einstein dynamics or Standard Model
matter. The manuscript identifies the explicit matching conditions and the next
joint growing-size test. The original Hamiltonian and earlier findings are preserved.
