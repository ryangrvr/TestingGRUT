# GRUT-I2 independent review packet

Prepared for an external reviewer. This packet includes a separate
computational implementation of the registered small cases.
The implementations were prepared in the same project and do
not constitute independent external replication. There has
been no laboratory run or dispatched external review.

## Immutable inputs and replay

Parent is draft PR #21 at ca33f29f75e2f2d83e0a984b67af23a8adf35598.
The I2 freeze is 4dcc33c4720969b2898d7a0320a4a9ec3de9ea39.
Take the equations, thresholds, and axes from its files and
derive a comparator independently of the stored outputs.
[Microscopic model](MODEL.md) · [Frozen protocol](PROTOCOL.md) ·
[Parameters](PARAMETERS.json) · [Replay](REPRODUCTION.md)

1. Label each fuel molecule and derive the reverse rates when one
   label flips. Check that the count generator is a lumping and that
   the hidden count entropy is included. Challenge the solvent
   and mixing assumptions separately from algebraic balance.
2. Verify neighbor support at the patch, invariant gas energy,
   source switching and the receptor's missing read port in the
   sensor-off arms. Check whether the choice of shape favors a
   specific site or layout. Check the flat formation control.
3. Reproduce both first-exit and endpoint membership from the
   same start. Confirm that failures and paths that reenter
   contribute to the thermodynamic ensemble, even though they
   fail the primary survival score.
4. Track the net chemical conversions for the two channels,
   restoration entropy and residual correlations. Derive the
   ideal reset costs and state which later controls are needed
   to attain those infima. Price any additional apparatus
   separately. Test that the no-fuel curve agrees with passive
   survival and that reverse fuel moves remain at positive capacity.
5. Recompute the bounded control frontier in both reset cases
   using linear programming or an independent convex hull.
   Challenge the zero incremental benefit of an acquired record.
   The null should include uninformed stochastic activation and
   equal-rate sensing. Treat the direct cue as a stronger access
   class. Distinguish conditional history information from
   a benefit at matched work.
6. Derive the equilibrium run count and its scaling. Look for
   prior dependence, finite-size artifacts and false coarse
   thresholds. Do not equate the small one-patch experiment
   with uniform bulk drive.
7. Verify the conditional rate ratio without using future outcomes
   or an unexplained affinity. Compare the proposed laboratory
   contrast against ordinary chemical kinetics and state what
   additional evidence would be needed to claim any discrepancy.

Return the input hashes, reproduction conditions, survival and
cost residuals, your access and thermodynamic assumptions,
and the strongest verdict warranted. Mark whether the result
is an independent computation, derivation or hardware test.
[Derivations](DERIVATIONS.md) · [Full thermodynamic ledger](THERMODYNAMICS.md)
