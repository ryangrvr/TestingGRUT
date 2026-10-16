# WM2-S3: finite-receiver detuning interval

Separate candidate based on WM2-S2
9cbfef8b70f90c875df7b87b59b945c7bc99d02e.

S3 adds an explicit finite shared bosonic receiver and retains every
reverse term. For uniform storage detuning in [-1/2,1/2] and common
emission coefficient in [9/10,11/10], the actual finite closed law has
macroscopic formation probability >=1/2 throughout [4N^2,8N^2].
The receiver uses 64N^2 internal modes and at most N-1 quanta.

Full-H operational response fronts, nearby response and conditional
bulk-prefix persistence survive throughout this class. The real-energy
spectral gate fails for the full consistency sector: explicit coupled
Q poles approach the seed-energy window. A new complex-energy virtual
bound keeps eta>0 and must not be read as S2's real-energy theorem.

Read REPORT.md, MODEL.md and MATHEMATICS.md, then HANDOFF.md.
RESULTS.json records what ran. ASSESSMENT.json keeps physical gates
separate from numerical consistency in VERIFICATION.json.

From a repository checkout:

    node research/reality_framework/wm2_detuning_receiver_20261016/reproduce.cjs
    node research/reality_framework/wm2_detuning_receiver_20261016/reproduce.cjs --recompute

Default verifies saved evidence and content hashes. The second command
reruns every bounded control before refreshing the manifest. No package
installation or remote service is required. The same numerical functions
were executed in JavaScript V8 during construction; the Node filesystem
driver itself was not executed here because local shell access failed.
Last-digit differences between engines are not physical differences.

This is a restricted finite-window operational module, not a generic
spatial phase or complete WM2. The receiver spectrum, growing mode count,
normalization, bosonic statistics and one-dimensional grammar are inputs.
Independent site disorder, selected dimension, original charged-matter
universality, gravity and empirical novelty are unearned.

S2, S1, WM1, frozen CP-12, SCOUT-1 and canonical GRUT documents are
preserved through the parent tree. The publication was constructed
directly through the repository API after the local environment
disconnected; a new checkout/fetch is needed for the new local files.
