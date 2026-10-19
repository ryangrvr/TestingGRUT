# GRUT closed-whole working theory (packet `grut_closed_whole_theory_20261019`)

19 October 2026. Owner: D. Ryan Grover. Branch `codex/grut-closure-fork-20261019`.

This packet is the working theory requested by the handoff in
`../grut_closed_whole_handoff_20261019/`. Its parent object is a question, not an action:

> What does a closed system, whose every subsystem's bath is the rest of the universe,
> force beyond known physics?

Read in this order:

| File | What it is |
|---|---|
| [WORKING_THEORY.md](WORKING_THEORY.md) | The readable report: what GRUT now says, what calculates, what closure constrains, which proposed additions did not survive, what remains supplied |
| [FORMAL_CORE.md](FORMAL_CORE.md) | Definitions, exact global-to-local rules, what closure forces and does not force (with countermodels), the theory datum, what an additional law would have to be |
| [CALCULATIONS.md](CALCULATIONS.md) | One closed whole carried from preparation to measurement (calc 1) and the nested/overlap compatibility check (calc 2), with every number |
| [CAMPAIGN_LEDGER.md](CAMPAIGN_LEDGER.md) | Every sector as a campaign: question, assumptions, executed evidence, answer, classification, remaining freedom, failures retained, source |
| [ledger/](ledger/) | The detailed, line-cited ledger entries behind the table, one per sector, each adversarially re-read against its sources |
| [OPEN_CLAIMS.md](OPEN_CLAIMS.md) | Exact claim boundaries, retained failures, unexecuted tasks |
| [SOURCES.md](SOURCES.md) | Repository sources (paths, commits, branches) and primary literature actually read |
| `code/` | `closed_whole_demo.py`, `marginal_compatibility.py`, `verify.py`, `finalize.py` |
| `results/` | JSON outputs, the M1 matched-control re-run record |
| `MANIFEST.json` | sha256 of every packet file |

Reproduce: `python3 code/closed_whole_demo.py` (~1 min), `python3 code/marginal_compatibility.py`
(<1 s), then `python3 code/verify.py --rerun`. Requirements: Python ≥ 3.12, numpy, scipy, cvxpy
(CLARABEL solver, bundled with cvxpy).

Status in one line: **a working closed-system framework of known physics with conditional
model instances and computed boundaries; no additional law or state-selection rule has
survived; nothing here is validated universal physics.**
