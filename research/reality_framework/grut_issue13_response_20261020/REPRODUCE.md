# Reproduction and independent review

This branch adds a response packet to `2ee7f82`; it does not merge Claude's
branch. The exact Claude review is linked by commit in CORRECTIONS.md.
No original frozen science file is modified. The coordination log is updated.

From repository root, with Python 3.12, numpy, scipy and matplotlib:

```bash
python research/reality_framework/grut_3d_empirical_20261020/archive_io.py unpack
python research/reality_framework/grut_issue13_response_20261020/archive.py unpack
python research/reality_framework/grut_issue13_response_20261020/analyze.py
python research/reality_framework/grut_issue13_response_20261020/plot.py
```

The first command authenticates and extracts all 190 original raw members.
The second restores all 36 exact new raw outputs: 30 files are archived,
and six endpoint volumes are reconstructed from their immutable parent
volumes plus 83 explicitly stored byte replacements. All source and target
hashes are checked. No new measurements are silently replaced by old data.
See RAW_MANIFEST.json. The archive is deterministic, with zero timestamps.

To repeat the evolutions, using g++ and OpenMP:

```bash
python research/reality_framework/grut_issue13_response_20261020/instrument.py
python research/reality_framework/grut_issue13_response_20261020/analyze.py
```

This reruns the original six commands using current paths and overwrites
the extracted `runs/` copies. Preserve a separate copy before rerunning
if exact output-byte comparison is wanted. OpenMP reduction order can change
scalar histories at roundoff and a few endpoint float32 least significant
bits, as already quantified. A rerun should be assessed by the trajectory
tolerances and physical diagnostics, not advertised as a byte-identical
independent experiment. `archive.py unpack` restores the published originals.
`archive.py verify` checks those exact restored bytes.

RUN_ENVIRONMENT.json records compiler options, original and instrumented
source hashes, executable identity and runtime. Each execution record has
the full command, case, exit status, stdout, stderr, elapsed time and the
parent execution identity. Compiler products are generated outside the
packet and are not required as opaque binaries. The recorded environment
uses numpy 2.3.5, scipy 1.17.0, g++ 13.3.0 and 8 OpenMP threads per run.
The new build omits the parent's `-march=native` optimization for portability;
both keep `-ffp-contract=off` and omit fast-math. Exact build commands are
recorded. The replay errors quantify numerical agreement across these builds.

RESULTS.json contains all cosmological roots, all endpoint/old-history
diagnostics, replay errors, conservation measures, analytic scales, the
108 assertions and the analysis environment. PACKET_MANIFEST.json inventories
the response packet apart from itself. Its hashes establish file integrity,
not an independent scientific reproduction.

Independent-review checklist:

- Read the unchanged action and the previous quantum matching report. Check
  that the total gravitational vacuum constant is distinct from the subtracted
  excitation energy, without double-counting the polynomial constant.
- Derive the charge floor with unit charge, gradients and D positivity; check
  the normalization and the physical baryon-density conversion.
- Rebuild the six-row comparison independently before importing our functions.
  Disclose familiarity with these public measurements. No blind test is claimed.
- Verify the strict inequality/equality distinction at the vacuum boundary.
  Check the curved sine factor, positivity of H squared, pre-antipodal
  monotonicity and the distinction between omega_k and Omega_k.
- Inspect the tiny source instrumentation. Confirm that it writes no state
  and uses no staggered momentum for per-step energy or charge inference.
- Check the original +v boundary conditions and single-vacuum production
  preparation. Neither supplies a thermal network ensemble.
- Recover raw data from the archive and endpoint patches, compare hashes,
  and rerun evolution on a separate implementation if possible.
- Check that sign changes and negative-vacuum plateaus are distinguished.
  Inspect threshold controls and their limitations; do not infer topology
  or a mass function from threshold counts.
- Check finite time, fixed dt, all three spatial grids, nonmonotone production
  extrema, unresolved population and physical-history uncertainty.
- Re-derive tree wall tension and the high-T trace coefficient; keep the
  zero-chemical-potential equilibrium assumptions explicit.
- Before applying a CMB constraint, supply a network history, tension/source
  normalization, expansion/backreaction treatment, calibrated physical
  parameters and a separate observational protocol.
- Confirm that the data justify only the claim in each row of
  CLAIM_EVIDENCE.md. No field identification, external scientific replication,
  cosmological lifetime or finished theory has been established.
