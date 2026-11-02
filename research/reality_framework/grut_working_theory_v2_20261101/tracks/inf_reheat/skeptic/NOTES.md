# inf_reheat skeptic (key inf_reheat_skeptic): pointer

The write-up is `../SKEPTIC.md`. This folder holds its code (`code/sk_tmodel.py`, `code/sk_channels.py`,
`code/sk_summary.py`) and output (`results/sk_tmodel.json`, `results/sk_channels.json`, `results/sk_summary.json`).

Run order, in `code/`, one thread: `sk_tmodel.py` → `sk_channels.py` → `sk_summary.py`. The track's
`../results/*.json` are read-only inputs.

Same-session AI work on one supplied model. No validation or exclusion is claimed; every choice is a model choice.
