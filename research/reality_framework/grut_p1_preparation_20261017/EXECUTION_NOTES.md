# Execution disclosure

Charter freeze: 44f8e82d046c141713e6cc6d61698d06c2abc88e.
Analytical specification and initial control-source freeze:
44b2748390c2e0113be29d090669bd1bd6f91f36.
Baseline D1: 830746b180ac8af6d5dc40a9e7777740981ceee7.

The design and elementary analytic estimates preceded the charter; this is
not represented as blind preregistration of an unknown mathematical result.
Parameters and the primary architecture were frozen before executable controls.
The finite diagnostic plan was written before those controls ran.

Primary execution 1 failed in final formula evaluation at N=1024 because a
large integer factorial could not be converted to float. Earlier controls
had returned. No CHECKS.json was produced. FIRST_RUN.txt archives that failure.
The correction uses log-gamma above order 170 with a positive exp(-700) floor
for smaller-than-floating-range tail terms, instead of an apparent exact
zero. This is an ordinary-double formula evaluation, not a rigorous enclosure.

The entire original bounded suite then passed, unchanged in physical law,
times, sizes, parameters and error/resource gates. CHECKS.json records that
execution: 2.459034583 seconds, 117.38671875 MiB peak resident memory. The
60-second and 1-GiB limits were not relaxed. OPENBLAS_NUM_THREADS=1 and
OMP_NUM_THREADS=1 prevent unnecessary parallel numerical workers.

`EXECUTED_SOURCES.json` retains both the original freeze hashes and the
corrected control-source hashes. Source pins for seven archived baseline
files were checked before and after execution. All size-dependent theorem
claims are analytical; evaluations at 256,512,1024 are not wavefunction
simulations of those sizes. The inherited large threshold is only recorded,
not simulated. No experimental data, external reviewer, automatic proof
checker or independent all-size validation exists in this packet.

The independent two-site tensor builder checks the path implementation's
entries; the fixed-path and complete-geometry builders otherwise share their
explicit transition code. They must not be described as wholly independent
implementations. The small geometry receiver has two modes and all allowed
occupations. It is not substituted for the physical receiver in a formation
claim. All physical charge sectors at that finite size are retained, including
more pairs than the selected pump can populate from vacuum.

The record-energy paragraph and original-trial signed-response corollary were
added to the exposition after the controls; neither changes the executed
Hamiltonians, gates or primary success condition. They are analytical claims
for review, not additional measured quantities.

Reproduce without overwriting the archived observations:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py --run-checks

Default `python verify.py` only verifies owned-file hashes, execution-source
hashes, unchanged baseline pins, scope fields and archived control gates.
It does not independently validate the mathematics or roundoff contract.
