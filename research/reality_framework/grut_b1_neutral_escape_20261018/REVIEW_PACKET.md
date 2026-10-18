# Short adversarial review route

Status: author-derived proof and internal numerical checks. No external
reviewer, independent implementation, independent validation or empirical
test. Builder/challenger/integrator/recorder are roles performed by one
assistant sharing context and code. No AI agreement is counted as evidence.

## First inspect the claim, then try to break it

Claim: `GRUT_CANDIDATE.md`, exactly one initial resource, permitted fuel
interval, no later switch, >=64-event fixed path. Matter probability >1/2
for theta in [2,4]; at theta=4, co-pair arrival x>=4 with color success >3/5.
The conditional moving-geometry extension additionally requires the large
accepted record specified there. A1 remains closed. There is no density claim.

Read `CENTRAL_RESULT.md` in this order:

1. Rebuild the complete (m^2+1)-state one-resource matrix from the primitive
   charge, electric, density, fuel and conversion operators. Check all signs
   and positive E^2; do not begin from the effective walk.
2. Verify B^*D^-1 B, including both virtual hop orders, the end diagonal,
   the supplied counterterm, and the exact nominal K of equation (1).
3. Multiply out HZ-ZK from (2). Check parity cancellation B^*x2=0 and every
   term in (3). This is the highest-value short algebraic challenge.
4. Check residual norms and initial Z|f>-|f>; then derive (4) by Duhamel.
   Do not assume Z is isometric or claim closeness of the full state to
   the *undressed* walk. Only P-supported amplitudes use RZ=R.
5. Verify the half-line normalization, finite-end drive estimate and exact
   rational enclosure of the first five Bessel terms in (8). Check that
   distance x and walk index n differ by one, because n=0 is fuel.
6. Check parameter-perturbation norm, continuous-time survival bound, and
   conversion of arrival amplitude to probability. No time-grid inference.
7. Independently rederive (11) for the actual moving geometry. This is the
   longest dependency: preserve the actual conditional receiver, all q
   errors, all charged separations, boundary degree correction and fuel-
   supported pop gate. Prior P1/A1 proofs motivate its form, not its truth.
8. Check exact one-resource color factorization, erasure accounting and
   energy bounds. Only then compose recording, fuel loading and propagation.

The finite checks do not validate step 7 or the all-N mathematics. If step 7
fails but 1-6 hold, retain the fixed-path theorem and mark the full-geometry
extension unproved. If S9 fails, discard only the original-seed acceptance
corollary; the conditional record theorem does not depend on S9.

## Exact geometry-generator reference

For clarity, the inherited geometric adjacency operator sums over every
rooted simple word w and every free label v, mapping |w> to |wv> with head
and directed bond updated. The child must be empty of matter/fuel and the
new/removed link must have E=0. A pop also requires the removed head empty.
Assisted append multiplies q_v^+ b_contact^dagger and -g/N; its adjoint pops.
Bare append multiplies -(n_head+n_child)/(4sqrt(N)), leaves q and receiver
unchanged and has its exact adjoint. These terms and the diagonal in
`GRUT_CANDIDATE.md` define H_geo over all storage-error configurations.
There is no event-label ordering, preferred pair matrix, external geometry
schedule, or projected q=a evolution. Root and endpoint are relational
flags of the recorded word. On vacuum fuel/matter all new B1 terms vanish,
so the old formation block is an exact embedding, not a module analogy.

## Source compatibility ledger

| Source | What may carry forward | What may not carry forward |
|---|---|---|
| S9 `82bcab3` | Same fuel/matter-vacuum formation block, pending review | Fuel loading, B1 lifetime or validation |
| IG1 `650104c` | Necessity of declaring a generator support rule | Derived universal locality |
| D1 `830746b` | Warning about preparation and charged energy | Its externally populated density theorem |
| P1 `643ec7f` | Positive E^2, bounded links, confinement, physical gates | Switched-off conversion operation as B1 dynamics |
| A1 `22916ef` | Closed negative gate; virtual pair algebra as starting identities | Its transient preparation, gate, source law or times |
| B1 | New neutral conversion, fuel level, counterterm, color and proof | A theory of gravity or finite-density world |

All cited local files are byte-pinned in `SOURCE_PINS.json`. The five core
specification/proof/code sources and pins were hashed in
`EXECUTED_SOURCES.json` before primary execution.

## Reproduction

From the repository root, Python with numpy/scipy installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python research/reality_framework/grut_b1_neutral_escape_20261018/verify.py
```

That checks integrity without rerunning trajectories. To rerun the frozen
small controls, replace `verify.py` by `checks.py`; JSON is printed to stdout.
No packages are installed by either script. The executed environment used
Python 3.12.14, numpy 2.3.5, scipy 1.17.0. Fraction arithmetic checks scalar
inequalities exactly; matrix diagonalization uses double precision and has
no global roundoff enclosure. Execution consistency is not mathematical review.
