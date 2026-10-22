# Sources and reproduction record

22 October 2026 UTC. Sources are credited to their authors. The review
does not claim the known mechanisms as GRUT discoveries. Source read
levels below describe this reconciliation; earlier packets retain their
own, sometimes fuller, source records.

| Source | Read/check level in this reconciliation | Use |
|---|---|---|
| Claude checkpoint [78352e2](https://github.com/ryangrvr/TestingGRUT/commit/78352e2c24d2a9edee113b0d789b984eff77ac55) | Handoff, stage-2 brief, all four build dispositions and completed verifier dispositions; relevant code/results inspected | Exact object reviewed; original files preserved |
| Codex checkpoint [ca2c867](https://github.com/ryangrvr/TestingGRUT/commit/ca2c86776a715b646c6fd3b60af10518dd9c84ac) | Prior same-author full build, derivations and recorded checks, reread here | Independently formulated comparison, including Einstein-frame check and false-phase example |
| [KV 0806.2805](https://arxiv.org/abs/0806.2805) and [1102.3152](https://arxiv.org/abs/1102.3152) | Relevant equations/state reading via the two archived derivations; 0806 also directly read in the earlier Codex packet | Four-form action, scalar-tensor map, final-state matching cost |
| [KV16a 1601.00601v6](https://arxiv.org/pdf/1601.00601v6) | Primary PDF: relevant source, backreaction, asymptotic and appendix equations; archived notes/code compared | Phenomenological curvature-production law; new constant-K finite-history analysis |
| [KSV16 1601.04676v4](https://arxiv.org/pdf/1601.04676v4) | Primary abstract and pertinent equations through archived transcription; exact static family rederived | Dissipation vanishing for static de Sitter; no new full basin-grid reproduction |
| [K22 2207.03453v4](https://arxiv.org/pdf/2207.03453v4) | Primary action, cosmology, analytic solution, table and endpoint discussion; new table integration | Separate brane realization; source and original-action distinction |
| [Lee et al. 2002.11761v1](https://arxiv.org/pdf/2002.11761v1) | Primary paper and fresh source figure; figure rendered and visually checked; archived SVG digitizer rerun | Unscreened alpha=1/3 range, approximately 54 micrometres |
| [Martin–Robertson 1907.02500](https://arxiv.org/abs/1907.02500) and [SMDR v1.3](https://www.niu.edu/spmartin/SMDR/) | Official program description; source archive downloaded, build files/driver inspected, library compiled locally and executed | Pole/MSbar matching and beta functions for a new minimum-finding calculation |
| [Hiller et al. 2401.08811v3](https://arxiv.org/html/2401.08811v3) | Primary section II and Table 1 | Cross-check near 171.04 GeV; mass-definition and correlation caveats |
| [Buttazzo et al. 1307.3536](https://arxiv.org/abs/1307.3536), [Froggatt–Nielsen hep-ph/9511371](https://arxiv.org/abs/hep-ph/9511371) | Relevant formulas and credits in prior derivations and archived track; not fully reread in this review | Existing SM stability and multiple-point conditions, not a new prediction |
| Hu–Sawicki, Starobinsky, Appleby–Battye and Wang–Hui–Khoury | Archived bound/verifier analyses and their source read levels | Scope of the specific screened examples; no new class-wide theorem |

The PDG values used for the retrospective arithmetic are the snapshot
already exposed in Claude's record. The review does not present them
as a fresh survey of current world averages. Likewise, the chosen
rho_Lambda, age and thermal inputs are reference values, not held-out
measurements. No newly computed observational likelihood is reported.

## Computational reproduction

From this construction directory:

```sh
python tracks/codex_review/code/review.py
python tracks/relate/code/freeze.py
python tracks/codex_review/code/verify_package.py
```

The first needs numpy, scipy and sympy. The second performs the offline
relation replay from archived SMDR outputs; it cannot recreate the
library computation without a build. Regenerating files changes their
manifest hashes, so use a separate checkout for reproduction or retain
the published originals for comparison. The package verifier checks
published hashes and dispositions without rewriting them.

For the full SMDR check, obtain the official v1.3 archive and verify
SHA256 `398b3947a6e5e02244d5adcd45911a3825e056a1b6f65bcf203975da0fc91b41`.
Build its libraries in a separate directory, compile the archived
`tracks/v_vacua/code/smdr_driver.c` against them, and copy its
`benchmark_data/ReferenceModel.dat` to that build's working directory.
Then run:

```sh
python tracks/relate/code/freeze.py --smdr-home /absolute/path/to/SMDR-1.3
```

The new minimizer roots interpolated beta_lambda on an 801-point flow
and brackets the top mass. The recorded three-/four-/five-loop boundary
values reproduce Claude's 241-point parabolic-minimum results to better
than 0.1 micro-GeV. A direct full SMDR fit at the critical point checks
the matching interpolation locally; this does not validate all possible
future input points or the physical theory error.

The official library is reused, rather than independently recoding its
beta functions. Its GPL code is not vendored here. The archive,
reference-model, driver and compiled-library hashes are recorded. No
system-wide installation or `make install` was run. The first library
build hit archive ownership metadata errors; safe extraction of the
nested dependencies corrected that environment issue. The first fit
needed the explicit reference-model file. Neither required a physics
or input change.

The review's first Python run found SymPy missing. SymPy 1.14.0 and its
dependency were installed into an isolated scratch directory, then the
complete new calculation passed. Numerical library versions appear in
`results/review_results.json`.

The Eot-Wash replay uses the source `fig5b1.pdf`, `pdftocairo -svg` and
the preserved `tracks/v_bound/code/v3_eotwash_digitize.py`. It is an
independent download with the same digitization algorithm, not a new
independent algorithm. Raw crossings, the alpha=1 anchor correction,
source hashes and converter version are in `results/eotwash_replay.json`.
The anchor-corrected alpha=1/3 crossing is 54.0236 micrometres, or
3.65261 meV. Digitization precision is not a new confidence interval.

## Verification boundary

Newly run: 17 symbolic/numerical review dispositions, five relation
replay/reproduction dispositions, fresh laboratory curve replay and
package integrity checks. Prior Codex build checks remain in its
unchanged packet. Original Claude checks are credited as archived
results except where the new calculations explicitly reproduce them.
No full CMB/SN/BAO refit, entire old basin grid, quantum tunneling,
microscopic dissipation, all-order effective potential or outside
replication is claimed.
