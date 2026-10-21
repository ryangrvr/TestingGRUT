# Retained implementation failure

The first expanded archive check selected all scalar CSVs, including validation
controls, then required exactly twelve production runs. It stopped before writing
new results: `runs: 31`, maximum pressure-identity relative error
`2.36871690980499e-16`. This was a selector/count failure, not a failed field identity.
The selector now excludes the archived `v_` and `g_v_` validation controls and
requires the twelve production records. No threshold was loosened. No new lattice
run, empirical comparison, or physical parameter adjustment was attempted here.
