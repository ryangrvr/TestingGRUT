# Execution and evidence contract

Date: 2026-10-18 UTC. This was a bounded consolidation/local-algebra cycle,
not a new physical dynamics experiment.

## Frozen execution

1. Scope and historical input hashes were published at
   43bf3afe6b91f7aff0a0d74746c96447de6095fd.
2. The theorem, moving-law diagnosis, review dossier, challenge disclosure
   and Python source were published at
   e6d2b6a35236ae3cc6596ecda4c06cb5309495b0 before execution.
3. EXECUTED_SOURCES.json records SHA-256 for all ten frozen files.
4. `python run_controls.py primary` ran once; then the unchanged
   `python run_controls.py replay` ran once. No fallback ran.

Both workers exited zero. Their CONTROLS.json output was byte-identical,
SHA-256 c702ffba8953ad2cb6a25b7f98ddcd49db2238fc416d2ee38c618be7d87cc87d.
Primary wall time: 0.047508571 seconds; replay: 0.062356902 seconds;
combined: 0.109865473 seconds. Peak child RSS: 10,496 KiB.
Each child had a 25-second CPU/wall ceiling, 1-GiB address-space ceiling
and 25-MiB file-size ceiling, within the combined 60-second envelope.
No code or proof was changed in response to execution.

## What was checked

All 39 checks passed: 24 local algebra/margin checks plus 15 pinned
historical inputs. The exact finite local matrix has 67 scalar states.
It checks closure, Hermiticity, resource/field constraints, unique outward
two-operation coefficients, the degree-energy difference and initial
energy moments. Separate elementary gate checks exhibit forward append
past a fueled head with an admissible full-receiver occupation. Rational
arithmetic checks operating times and the CONDITIONAL future error budget.
The append check is a direct reading of declared gates, not an independently
generated complete moving-Hamiltonian matrix.

All deciding arithmetic uses integers/Fraction; decimal outputs are
display only. Exact local checks do not enclose the archived spectral
integrals or constitute an all-size analytical proof. There were zero
trajectory and zero operating-length simulations. The hypothetical
epsilon_geo=1/1000 was not measured or proved.

## Reproduction and integrity

The saved runner deliberately refuses duplicate primary/replay metadata.
To check the delivered evidence, run `python verify.py`. It checks the
current manifest, frozen source, 15 historical pins, all 19 files in the
archived B2-RA manifest, identical worker output and resource ceilings.
For a separate replay, use a fresh checkout of the source-freeze commit,
add the unchanged EXECUTED_SOURCES.json record and run primary/replay
there, leaving this historical record untouched.

MANIFEST.json covers every final packet file other than itself. Git commit
identity additionally pins that manifest. Publication is on the isolated
consolidation branch; no canonical merge or historical overwrite occurs.

Internal shared-context challenge, execution consistency, independent
mathematical review and empirical validation remain distinct. Only the
first two were performed here; independent review is pending and no
physical experiment was performed.
