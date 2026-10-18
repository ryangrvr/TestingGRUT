# Bounded execution and reproduction

Scope was frozen at c616af32c94507d2151f263563fa7421d2ddf644.
Proof/code were frozen at 4bea0bdb0c7d77a1d965908bce0ac3f65f13e3ba
before primary execution. EXECUTED_SOURCES.json hashes all nine frozen
files. No source or proof was changed after observing results.

`python run_controls.py primary` ran once, followed by the unchanged
`python run_controls.py replay` once. Both exited zero and emitted
byte-identical CONTROLS.json with SHA-256
8b2905d41c5884856b60118cc1c3282e8f381c002b8d29d07fdc764759b67e55.

Primary child wall time: 0.169137815 seconds; replay: 0.175596973 seconds;
combined: 0.344734788 seconds. Peak child RSS: 13,568 KiB (13.25 MiB).
Limits: 25-second CPU/wall per child, 60-second combined envelope,
1-GiB address space and 25-MiB files. No fallback or relaxed budget.

All 31 named checks passed. Seven check historical pins; 24 check finite
primitive algebra/margins. The N=4 whole-word P control has dimension 84
and 160 boundary image states. Every B*B matrix entry is compared with
the formula, and reverse geometry entries are checked. Minimum-word
neighborhood counts are 1,3,9,27 for H^3 and 1,2,6,12,25 for H^4.

Deciding arithmetic uses integer/Fraction and finite exact square-root
sums, including bosonic square-root coefficients. The controls implement
off-diagonal geometry and matter primitives. Diagonals and receiver-only
hopping cannot enter these minimum-length matter/geometry targets; their
absence is not a truncated dynamics approximation. The source-defect and
storage-energy checks separately use the actual diagonal formulas.

The receiver test vectors have contact/spectator Fock occupation, with
all other actual physical modes in vacuum. No time evolution restricted
to those occupied modes is performed or claimed. The mathematical
theorem retains arbitrary actual conditional sigma_w. There are zero
dynamics runs, zero operating-length simulations, no interval propagation
and no achieved late-time comparison or information-error estimate.

Run `python verify.py` to check the final packet, frozen source, seven
input pins, 18 parent-manifest entries, replay identity and resource
ceilings. The saved runner refuses duplicate metadata. For a separate
reproduction, use a fresh checkout of the source-freeze commit, add the
unchanged EXECUTED_SOURCES.json and run primary/replay there.

MANIFEST.json covers all packet files except itself; the final commit
pins the manifest. Publication stays on the isolated right-boundary branch.
No canonical merge, historical overwrite, installation, paid computation,
persistent service, external reviewer contact or physical experiment took
place. Internal shared-context challenge and exact finite verification
do not discharge independent proof review.
