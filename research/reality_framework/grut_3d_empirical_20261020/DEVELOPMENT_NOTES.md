# Source and execution boundaries

The main Cartesian 23-case protocol and solver were locally frozen before any
tau>6 Cartesian result was inspected. They were not publicly preregistered.
The later spectral, large-domain original-mode, geometric-H0 and direct-impact
protocols explicitly record their additive, posthoc scope. Their hypotheses and
criteria were saved before their respective new results were generated. No
earlier frozen packet, kernel, profile or result has been rewritten.

Development pilots consisted of a T=3/6 U(1) symmetry check, a centered
finite-difference tangent check, a rectangular-indexing check and timing tests.
A first T=4,dt=.075 pilot was rejected by the integer-step check before stepping.
The development runner first failed JSON serialization of a Path after a short
C++ run; recording was corrected. A short L=12 cube versus Lz=18 test gave a
relative R8 charge difference 1.0717926821e-7, just beyond its 1e-7 criterion.
This test also changed the physical boundary. It was repeated with transverse
L=18 and longitudinal Lz=24, keeping the core causally better separated from
the boundary, and returned zero at printed precision. The final short controls
and their diagnostics are retained. The earlier overwritten development-pilot
histories were not retained; no failed long-run science criterion was erased.

The first multi-case tool session lost its execution transport while the fine
reference pair had reached tau65. It had no completed execution record. The
partial CSV and interruption note are retained. The case was rerun with the
same frozen source, parameters and criteria. `verify_packet.py` compares the
overlapping diagnostic prefix. The host cgroup reported zero OOM events at the
interruption; its underlying transport cause was not established.

Main execution uses sequential cases within each runner. Independent later
control batches may run concurrently, with separate build directories to avoid
replacing an executable in use. Sources, field inputs, grid spacing, timestep
and criteria remain fixed; every exact command and execution duration is saved.
Wall time is not a physical observable.

The complete raw scalar diagnostics, endpoint full-grid component lists and
float32 subsampled spatial volumes are archived. Full field/momentum histories
are not archived. The programs and frozen preparations reproduce the runs.
`run.py preflight` is a development/freezing command and rewrites PREFLIGHT.json;
run it only in a separate working copy. For analysis of the frozen evidence,
unpack and use `verify_packet.py` instead.

The common physical point, source boundary, unpublished-Claude limit and the
remaining empirical gaps are summarized in START_HERE.md and REPORT.md. This
work does not claim to include Claude's unpushed tests or to supersede the
repository's historical canonical charter with an experimental branch result.
