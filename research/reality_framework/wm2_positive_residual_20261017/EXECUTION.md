# Execution disclosure

The pre-trajectory charter, physical controls, five-event recovery, fresh
S5 baselines and gated six-event worker ran in that order. Both five-event
root components completed. Every primary worker had the frozen 600-second
budget and ran sequentially. The six-event run stopped due to rank/error
constraints before that limit, not because it failed to form.

`EXECUTED_SOURCES.json` pins the numerical source used for the controls and
all primary workers. Those files were not changed after the first primary
component. The whole-worker timers begin before imports and report through
result construction, before the final JSON serialization/print. Peak RSS is
the Linux process high-water mark at that point. Native rebuild/validation
costs are separate. No installed package, implicit bath or outside process
was used to alter a physical trajectory.

The six-event worker was invoked through `reproduce.py --six`, which checked
the saved recovery gate first and wrote `replay/s7_N6_root0.json`. After it
terminated, that exact JSON was copied to `results/s7_N6_root0.json` for the
final archive. It was not rerun, edited to change status, or promoted to a
target-time estimate. The other root would use S5's proved mirror identity;
an unresolved component does not earn a family formation claim.

After the campaign ended, NEXT_RESIDUAL.md derives a state-weighted positive
residual bound. Its separate small controls do not affect the frozen runs.
The first control harness passed a strided matrix column to the raw native
ABI (which requires contiguous buffers), causing its assertion to fail.
That harness was corrected to copy columns contiguously. The final test
also constructs its physical B using the independent immutable S5 Python
operator. No primary-worker buffer follows that strided construction; their
vectors are contiguous by construction. No production source or primary
trajectory was changed for the post-campaign check. Raw native interfaces
must be called with their declared dtype/contiguous-buffer contract.

S7's primary controls and their archived outputs are unchanged. Saved-record
verifiers passed for S5/S6. Full all-size proofs, generic rank growth and the
compiler/BLAS assumptions have not received independent external validation.
The weighted control is an ordinary-double diagnostic; primary arithmetic
enclosures have the narrower declared model in ARITHMETIC.md.
