# WM2-S4: prepared bulk dynamics and independent disorder

Read [REPORT.md](REPORT.md) for the result, [MATHEMATICS.md](MATHEMATICS.md)
for proofs, [MODEL.md](MODEL.md) for the law/preparation/energy ledger,
and [HANDOFF.md](HANDOFF.md) for the next bounded question.
[S3_PROOF_CHECK.md](S3_PROOF_CHECK.md) records the same-Builder review
and its independence limit. [CHARTER.md](CHARTER.md) was committed before
S4 numerical trajectories.

The pre-execution charter commit was local. Shell GitHub push lacked
authentication, so API publication uses its own commit metadata.
`CHARTER_FREEZE.json` preserves the original raw commit object and
charter blob hash; the verifier checks both. This preserves integrity
and local execution provenance, not independent timestamp attestation.

S4 proves a conditional bulk response-feedback window under independent
storage detuning. Nonshrinking-disorder, all-size formation remains open.
It preserves S3's coupled real-energy poles and does not claim an
isolated low-energy spatial phase.

From the repository root, with Python, NumPy, SciPy and Node on PATH:

```sh
python research/reality_framework/wm2_prepared_disorder_20261016/reproduce.py
```

This verifies saved result gates and every packet content hash. To rerun:

```sh
OPENBLAS_NUM_THREADS=1 python research/reality_framework/wm2_prepared_disorder_20261016/reproduce.py --recompute
```

Recomputation includes the entire archived S3 Node suite in a temporary
copy, independent Python sparse matrices, exact type-word isometries,
first-error configuration searches, physical-receiver corner tests,
a direct 576-mode many-boson comparison, small-mode controls, full-H
phase probes and large-mode packet residual checks. It rewrites only
this S4 packet's evidence and hashes. S3 and all earlier files stay
unchanged. Runtime is ordinarily minutes rather than a large simulation.
The largest full vacuum matrix here is dimension 333,505; the largest
full one-probe matrix is 15,600. All charge-allowed occupations are used.

`RESULTS.json` contains raw finite results and execution disclosure;
`ASSESSMENT.json` separates proved, observed and open claims;
`VERIFICATION.json` summarizes checks; `manifest.json` hashes the packet.
Ordinary double precision is used. Analytic receiver-cut errors exclude
roundoff. Neither saved hashes nor a successful run certifies all-N
proofs or independent external validation.
