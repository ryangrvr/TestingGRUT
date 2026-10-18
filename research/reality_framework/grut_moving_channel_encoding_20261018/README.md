# MC-E1: what the moving channel can preserve

Start with [REPORT.md](REPORT.md), then [RECOVERY_THEOREM.md](RECOVERY_THEOREM.md)
and [FINITE_DIAGNOSTIC.md](FINITE_DIAGNOSTIC.md). [PREDICTIVE_INTERFACE.md](PREDICTIVE_INTERFACE.md)
specifies measured quantities and held-out tests. [CHARTER.md](CHARTER.md)
freezes scope; [INTERNAL_CHALLENGE.md](INTERNAL_CHALLENGE.md) and
[EXECUTION_NOTES.md](EXECUTION_NOTES.md) distinguish evidence categories.

New: an optimal color-only recovery bound and an exact full-law N4 finite-order
diagnostic with nonparallel successful environment vectors on a short interval.
Inherited: identity/SWAP closure and the triplet noiseless code. MC-E1-T is
an explicit qutrit-input protocol using that code, not a new law or an
arbitrary-two-qubit operating pass. Its actual output u and retention u_T
remain unbounded on the macroscopic moving windows.

Run `python verify.py` for saved hashes, exact predictions, primary/replay
identity and historical preservation. This is archive verification, not
another physical run. The bounded runner refuses duplicate primary/replay
execution; reproducing it requires a fresh source-freeze checkout and the
same EXECUTED_SOURCES.json. Ordinary replay is not independent validation.
The mathematical arguments and exact finite results are author-derived;
independent review and physical validation have not occurred.
