# Reproduction and continuation

Python, NumPy, SciPy, Matplotlib, gcc and g++ with OpenMP are required. No package installation, paid compute or AI subagents were used. Ordinary numerical processes run the histories. Runtime and thread counts are recorded per case; they are not universal timing estimates.

Start from this branch. A shallow/partial clone may need the pinned review/build objects:

```bash
git fetch origin fd71d3a416db929252561d2cf296e915dd710f37
git fetch origin e71730c24fe6a0b28ddd72cd30982c5b39c9aedf
```

The inherited kinetic bootstrap extracts its pinned completed interfaces read-only. The active fields still use the same couplings and fixed radiation clock. The new working directory defaults to `../grut_particle_transport_work` relative to the repository; `GRUT_PARTICLE_WORKDIR` overrides it. Inherited mode/reheating scratch directories are set under it by default. No old packet output is rewritten.

Run these scripts from the repository root, with `OPENBLAS_NUM_THREADS=1`:

```bash
python -B research/reality_framework/grut_particle_transport_20261027/code/package.py --sources
python -B research/reality_framework/grut_particle_transport_20261027/code/verify.py
python -B research/reality_framework/grut_particle_transport_20261027/code/prepare.py
python -B research/reality_framework/grut_particle_transport_20261027/code/pt_production.py
python -B research/reality_framework/grut_particle_transport_20261027/code/analyze.py
python -B research/reality_framework/grut_particle_transport_20261027/code/coherent.py
python -B research/reality_framework/grut_particle_transport_20261027/code/response.py
python -B research/reality_framework/grut_particle_transport_20261027/code/check_response_free.py
python -B research/reality_framework/grut_particle_transport_20261027/code/summarize.py
python -B research/reality_framework/grut_particle_transport_20261027/code/check_results.py
python -B research/reality_framework/grut_particle_transport_20261027/code/plot.py
```

`pt_production.py` accepts one or more case tags for a bounded rerun, such as `base48` or `IR12_s1_w1`. Default production uses two processes and three OpenMP threads each. `GRUT_PARTICLE_WORKERS` and `GRUT_PARTICLE_THREADS` override those counts. The archived infrared runs used two threads per process; the original-band runs used three. The Hamiltonian, time steps, grids and sampling are independent of that scheduling choice. The two volume boxes have equal spacing; the fine L12 run checks the coarse infrared result.

The retarded stage requires the main `base48` snapshots at s=20 and 20.44. Default response uses one numerical process with two OpenMP threads. All compiled objects use `-ffp-contract=off`; a content hash selects the stepping library. Temporary compiler outputs are replaced atomically. The saved classical arrays use double precision; selected temporal Fourier fields use complex64, followed by double-precision analysis. Publication retains full scalar traces, compact frequency spectra and averaged covariance, not every working temporal coefficient.

For the review-comparison check, extract the focused immutable inputs to `../review_fd71d3a`, retaining paths beneath `research/reality_framework/grut_formation_review_20261025`. `code/package.py --sources` can do this from the fetched commit and writes exact source pins. A partial clone may download missing blobs; a full clone already has them. This extraction is read-only with respect to the Git repository and does not merge the review branch.

The final check `code/package.py` verifies source hashes, finite arrays, links, manifests and the eight numerical groups. `--seal` records an intentional regenerated package after figure inspection. Cross-platform raw hashes can differ; compare the stated mathematical/numerical quantities before treating a byte difference as a physical failure.

The separate checkpoint named in [DATASET.json](DATASET.json) contains four full double-precision final states at s=34, plus JSON metadata. For example:

```python
import json, numpy as np
d = np.load('GRUT_particle_transport_checkpoint_20261027.npz')
metadata = json.loads(str(d['metadata_json']))
fields = d['IR24_s1_w1_fields'].copy()
# fields order: chi1, chi2, phi, pi1, pi2, pi_phi
# Read L, couplings and units from metadata before advancing this state.
```

These data permit continuation without repeating the early evolution. They contain a regulated classical field, not the universe's quantum state or a certified particle distribution. No physical abundance or threshold is reconstructed by choosing a different Fourier reference.
