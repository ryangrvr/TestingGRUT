# GRUT builder: source-resolved propagation

Continuation of `4e7d6c5a91fe96737ab9a8a140ec2b50eb140ebc`.

**Physical result:** the tested six-site recorded-ring metric fails its
finite source-resolved benchmark. At lag 1/4 a receiver at ring distance 3,
with an active direct chord, responds more than 550 times as strongly as
a receiver at ring distance 2. The original ring has over 90.4% retention.
The original-neutral ring formation probability also decreases across the
auxiliary sizes. This does not establish a spatial phase or an empirical
prediction, and is not a no-go theorem for all CP-12 phases.

The Hamiltonian, its gauge-invariant hopping, all external connections
and all graph rewiring are unchanged. The matched `Q=N+1` sequence is
explicitly auxiliary; a separate five-site `Q=2N` result is included.
No large-size scaling of the primary `Q=2N` sector is claimed.

## Read and review

- [REPORT.md](REPORT.md): results, interval tables, formation and record effects.
- [DERIVATION.md](DERIVATION.md): source response, complete instruments and bounds.
- [TEST_PROTOCOL.md](TEST_PROTOCOL.md): prospectively fixed six-site comparison;
  four- and five-site pilots were already inspected.
- [PROGRAM_TESTS.md](PROGRAM_TESTS.md): five GRUT gates and milestone decisions.
- [ARCHITECTURE_CHARTER.md](ARCHITECTURE_CHARTER.md): obligations of a replacement
  metric/transport hypothesis, without altering the frozen experiment.
- [HANDOFF.md](HANDOFF.md): concise continuation record for the reviewer.
- [SOURCE_RESPONSE.svg](SOURCE_RESPONSE.svg): static scientific figure from data.
- [ASSESSMENT.json](ASSESSMENT.json), [VERIFICATION.json](VERIFICATION.json) and
  [manifest.json](manifest.json): operational decisions, checks and provenance.

## Reproduce

Requirements used here: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0,
Matplotlib 3.10.8, a C++ compiler for the inherited signed-symmetry engine,
and approximately 10 GB RAM. NumPy/SciPy versions are recorded as the
tested environment, not an assertion that other versions cannot work.
Large numerical jobs run serially; simultaneous large jobs can exceed RAM.
Cold-cache reconstruction is substantially slower than output verification.

From the repository root:

```sh
python research/reality_framework/cp12_source_resolved_propagation_20261015/reproduce.py
python research/reality_framework/cp12_source_resolved_propagation_20261015/reproduce.py --full
```

The first command verifies the saved evidence and hashes. The second
regenerates the calculations and independent controls before verification.
Disposable caches use `GRUT_SHORTCUT_CACHE` and `GRUT_SOURCE_CACHE`, or
temporary directories by default. Caches are not required evidence files
and are not committed. Full reproduction depends only on repository
sources and the listed numerical/compiler dependencies.

`neutral_prepare.py` performs original full-basis neutral evolution.
`neutral_krylov.py` accelerates the large neutral calculation while keeping
its complete electric boundary and finite-H relation residual. Its
controls compare against full-basis evolution and a larger six-site
Krylov space. `source_basis.py` implements signed rooted lifts and full-H
Galerkin extension. `response_run.py` evolves both source vectors and
reports individual receivers, shell means, core retention and explicitly
final-conditioned derivatives. `refined_record.py` includes every record
outcome and raw rewiring coherence. The labelled N=4 controls check the
source calculation and actual nonselective recording disturbance.

All interval claims use analytic truncation/integration bounds evaluated
in ordinary floating point with stated allowances. They are not formal
outward-rounded certificates. Numerical consistency checks passing is
distinct from the physical acceptance test passing; the latter fails
for the specified candidate branch.
