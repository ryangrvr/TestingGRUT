# inf_chain SKEPTIC: index

The write-up is `../SKEPTIC.md`. Code is in `code/`. It is my own and independent of the track's code paths, with three
exceptions:
- the patched, read-only-sourced engine copy in `sk_gamma_switch.py`;
- the track's driver, `ch_common.py`, imported read-only by that script for the slice preparation;
- the shared SM running, `higgs_reheat/code/smpot.py`, imported read-only by `sk_higgs_TB.py`.

Results are in `results/*.json`, with the headline numbers in `results/sk_summary.json`. Labels and caveats are as in
SKEPTIC.md.

Every choice is a model choice. Nothing here is a validation or an exclusion.

Process note: a stray 2.7 KB copy of `results/sk_gamma_switch.json` was left at `/sk_gamma_switch.json` by a mistyped
`cp`. The harness's safety check blocked its removal, and the owner can delete it.
