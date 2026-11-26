# Reproduction

Run from this packet directory using Python 3 and the standard
library. The reference run used Python 3.12. No random seed,
external data or floating-point tolerance is needed. Each
calculation completes within the frozen 30-second budget.

~~~bash
run_dir=$(mktemp -d)
timeout 30s python replay_obstruction.py --output "$run_dir/moments.json"
timeout 30s python critical_replay_f1.py --output "$run_dir/resolvents.json"
cmp RESULTS.json "$run_dir/moments.json"
cmp CRITICAL_RESULTS.json "$run_dir/resolvents.json"
~~~

The first program checks exact conjugacy, port moments, graph
support, rank and the separation certificate. The independent
program checks projected inverses at $2n$ distinct shifts of
$V+zI$ per fixture and verifies rank from sampled inverse
columns. The degree of each entry's difference numerator is at
most $2n-1$ by the cofactor formula. Positive stiffness ensures
no singularity at the test shifts. Identical samples therefore
certify an identity. The observed equality of a finite number of
moments from the builder is likewise sufficient with the exact
conjugacy check. The two representations agree on all ranks,
graph invariants and the added-probe threshold. The dark control
has the expected deficient rank.
[Exact moments](RESULTS.json) · [Exact resolvents](CRITICAL_RESULTS.json)

The archived tree digest uses Git's NUL-delimited path list.
From the repository root it can be reproduced by the following
read-only commands. The index includes the exact text search and
selection expression for the documentary audit, with source Git
blob identifiers.

~~~bash
git ls-tree -r -z --name-only 14cf11c7fed1413d63f9139f25be9e1569a97f8c | sha256sum
git diff --name-status 14cf11c7fed1413d63f9139f25be9e1569a97f8c HEAD
~~~

The continuum result and the positive classification are analytic
proofs. A replay tests the finite witnesses and their controls.
Historical evidence is cited at its frozen status.
[Evidence audit](FOUNDATIONAL_EVIDENCE_AUDIT.md)
