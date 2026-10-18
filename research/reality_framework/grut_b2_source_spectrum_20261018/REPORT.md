# B2-SS: partial source results; the scaling gate stays open

**Verdict:** the focused source analysis does not yet prove or refute the
two growing-window source estimates. It adds exact observable definitions,
a source trapping exclusion, preparation-weight bounds and a continuous
finite return interval. B2-SO's comparison is unchanged.

## The probability accounting is clarified

The 34% condition is ONE-source co-pair probability. The comparator is an
unnormalized product embedded in the full physical space, and its joint
event weight is the square of the one-source weight. With all reconstruction
and comparison errors paid, that condition would imply joint output
probability above 11.458225%, sufficient for the frozen 10% gate. It would
not imply a joint probability above 33.8%. The condition itself remains open.

Fuel return is also precisely defined: in the exact one-resource sector,
occupied fuel is the sole spatial/electric configuration f. The return map
is G(t)I_2 on the input qubit, including an entangled reference, and total
one-source fuel occupation is exactly |G(t)|^2. This is unconditioned and
specific to the declared fuel seed. It is not the full two-source marginal
fuel occupation, whose subspace includes many matter configurations.

## New source results

All results below are author-derived, independent proof review pending.

* **No source eigenstate below or at the mobile band's bottom.** A positive
  supersolution on the entire configuration graph, including fuel and every
  charge separation, gives an exact ground-state representation. Its strict
  root residual excludes a threshold eigenstate as well as lower poles.
  Embedded and above-band localized components are not excluded.
* **Actual preparation weights are bounded.** The nominal fuel energy lies
  inside the exact ground band with margins greater than 1.5J. More than
  79.9% of its spectral weight is in that band's energy interval; at least
  5/9 is in a central interval compatible with the required bulk velocities.
  These are energy weights, not demonstrated outgoing probabilities.
* **Continuous return control extends through theta=16.** The same original
  fuel preparation obeys |G(t)|<0.65 for every 2<=Jt<=16. Through the frozen
  comparison, original two-source q=2 probability exceeds 33% throughout
  that interval for L>=2016 in the declared family. q=2 permits separated
  constituents; this is not co-pair arrival in distant receivers.

## What prevents an operating verdict

At the smallest allowed length the receiver window begins at Jt=252.
The finite return lemma ends at 16. Late return, embedded/source-localized
weight, outgoing receiver amplitudes, and their coherent spectral cross
terms remain unresolved. Neither a sampled favorable time nor the interval
spectral weight can fill that gap.

The exact remaining conditions are A_L(t)>=17/50 throughout every I_L and
|G(t)|<=2/3 throughout [2/J,L/(6J)]. No violating point was found, but no
uniform proof was completed. The joint event remains linked to these
quantities by the existing full-H theorem, without assuming physical
independent sources. No additional success probability is claimed.

## Verification and scope

The frozen rational checks and finite source controls at K=8,12,16 passed;
all five data files were byte-identical on replay. Combined workers took
1.18 seconds, peak RSS 51.79 MiB. No time evolution or operating population
was numerically run. Finite source spectra and a relative-coordinate bulk
control are diagnostics with ordinary double precision, not infinite-source
classification or outward-rounded certification. Code ancestry is shared.

No source law, physical window, threshold, or comparison constant was
changed. B2's old failure, A1's closure and the closed solver branch remain.
The scope is prepared fixed geometry with two supplied fuels and nominal
couplings, not a robust autonomous matter-bearing world. This cycle stops
with the source conditions mathematically unresolved.

Read OBSERVABLES.md for the exact conditioning/input quantifiers,
SOURCE_PROOF.md for derivations, and SOURCE_VERDICT.md for the remaining
spectral kernel and time-window estimate.
