# Fresh check of the S3 arguments

This is a same-Builder mathematical check and a second implementation,
not an independent external referee report. The user-supplied assessment
is treated as a guide; its referenced reviews/grut_review_report.md
was not available here and is not claimed as an inspected source.
The actual S3 source packet at 80f5d6a was read and preserved.

## Uniform formation proof

1. The symmetric word vectors have unique parents, giving exactly
   sqrt(N-1-k) for appending from k attachments. This validates the
   bosonic bookkeeping pump only for common detuning and coupling.
2. For lambda=g/N and N>=2, lambda^2<=121/400. The contact Weyl
   denominator is monotone outside the receiver band. Both edge signs
   exclude an exterior pump pole; the interior density has no real pole.
3. Completing the density denominator gives the dimensionless minimum
   factor at least 254/279>9/10 for |delta|<=1/2. Its inverse has total
   variation at most twice the inverse minimum. Together with edge
   numerator variation four this gives 40/(9 pi lambda^2).
4. Fourier integration by parts then gives A=1000/(729 pi) at
   t>=4N^2, so 1-2A^2=0.6186926... . This has a strict margin over 0.6.
5. The finite-receiver return path must use 2L_N+2 hops. The stated
   exponential comparison follows since
   16+128[1-log(8)] < -120. Resumming diagonal phases is essential.
6. The binomial pump expectation and Markov inequality yield the
   conservative macroscopic probability 1/2 for the full continuous
   symmetric interval and every time in the declared window. A grid
   of successful trajectories is not the proof of that interval.

No missing implication was identified in these scoped steps. This
does not establish independent-disorder formation: the invariant
subspace in step 1 fails, as S4 explicitly checks.

## Coupled complement poles

The all-spent, one-free-event configuration with probe on that free
event has diagonal 80+delta(N-1)+dGamma(h_rel). Its matter motion is
zero. Its sole bare append leaves Q; assisted pop annihilates a contact
quantum. The far-supported packet has zero contact occupation, so its
Q residual reduces to the displayed many-boson packet residual. Its
bare image in P has norm epsilon/sqrt(N). Spectral projection onto the
stated shrinking interval loses at most 1/(2sqrt(N)) of the witness;
the triangle inequality leaves nonzero B-dagger coupling. Thus a
coupled pole exists in that window for the declared finite receiver.

This establishes neither occupation by the physical preparation nor
divergence of every source-receiver block. S4 retains this distinction.

## Execution cross-checks

The exact archived Node driver was run, then fully recomputed in a
disposable copy under Node v24.19.0. All original files stayed unchanged.
S4 separately constructs sparse matrices in Python/SciPy, compares
complete matrix entries against the archived JavaScript construction,
and compares its exact type-word reduction against unreduced vacuum
wavefunctions. These are useful reproducibility checks. They do not
replace an external mathematical review or outward-rounded arithmetic.
