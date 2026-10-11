# GRUT III — changelog

## Revision 2 (2026-10-11, after three external reviews and an owner's note)

Corrections (all accepted; see `RESPONSE_TO_REVIEWS.md`):
- η = c³/16πG relabelled as the horizon's surface viscosity, units kg s⁻¹ (was Pa·s).
- "Memory time 1/(4H)" relabelled as the asymptotic decay rate of the pointwise ⟨T̂T̂⟩; memory now defined only for
  smeared observables (`07`, `v7`).
- a₀ Jacobson-route argument retracted; `06` is a coincidence appendix crediting Milgrom (1999).
- w ≡ −1 reframed as a falsifier of the minimal framework shared with ΛCDM; "DESI DR3" removed; DR2 figures cited.
- P1 restated in modular (Tomita–Takesaki) form; no factorisation or influence functional claimed from it.
- Ledger relabelled: 5 postulates, imports tagged IMPORTED (zero credit), data; "compression" withdrawn.
- P5 derivation via HJW withdrawn; reconstruction literature cited.
- v4 wording: √s vs frequency; staircase marked illustrative; on-shell s = 0.
- Priority: dS/conformally-flat noise kernels credited to Pérez-Nadal–Roura–Verdaguer and Cho–Hu; "not on the record"
  downgraded to "corollary; placement".
- Record correction D1 now proved (the record's S(w) is G⁺(w) by its own definition); gap quoted as 0.55%.

Additions:
- `07_ACCESS_RESOLVED_RESPONSE.md`: the Δρ residual framework (FDT identity, trichotomy, memory dictionary).
- `calc/v7`: memory from access geometry (reviewer's smeared-field formula; smeared stress-tensor kernel closed form;
  KMS for smeared observables).
- `calc/v8`: residual framework; the record's kernels as Δρ ≠ 0 by fiat.
- `calc/v9`: induced-gravity / entanglement-entropy matching per species; SM species sums (induced 1/G tiny or negative).
- `calc/v10`: gravitational sector = stationary thermal process + Brownian zero mode, D = H³/4π².
- `calc/v11`: conventional O(ε) KMS defect under slow roll (the first non-KMS baseline).
- `calc/requirements.txt`; `run_all.py` extended to v1–v11.

## Revision 1 (2026-10-11, commit 862e666)

Initial principle-first reconstruction: P0–P5, D1–D8, v1–v6.
