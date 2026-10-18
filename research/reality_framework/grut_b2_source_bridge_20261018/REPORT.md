# B2-SB: the missing bridge is explicit; scalar inference cannot supply it

**Decision: demonstrated limitation of the specified scalar-inference
route. Actual unchanged-source operation remains mathematically unresolved.**
No replacement time, Hamiltonian, source trajectory or B3 was introduced.
All new mathematical arguments are author-derived, with independent
review pending. The parent is B2-SS ac4bf51.

The exact decomposition is now in [DECOMPOSITION.md](DECOMPOSITION.md),
including spectral-type and energy cross terms, the physical detector
kernel, high-state reconstruction, the spatial cut, and the existing
unnormalized joint comparison. It gives, for the actual source,

    p_out >= [[sqrt(A_L)-sigma]_+^2-delta_L]_+^2,
    p_q2 >= [1-omega-(|G|+omega)^2-tau_L^2-delta_L]_+^2.

Here sigma<=6/27613, omega<=sigma^2,
delta_L=12/27613+1500 exp(-L/105)<.001, and
tau_L<1.01 exp(.3-L/105). These inherited bounds are uniform through
t=L/(6J). The missing physical functions are A_L and G on the frozen
intervals, rather than another refinement of those comparison constants.

## Two sufficient estimates that would connect the existing results

1. **Arrival.** Let c(t)=exp(-ith)1_[Omega-1.5J,Omega+1.5J](h)f and
   w_C=||c||^2>=5/9. A sharp orthogonality-aware inequality proves

       A_L >= [sqrt(w_C eta_C)
                    -sqrt((1-w_C)(1-eta_C))]_+^2,
       eta_C=||R_left,L c||^2/w_C.

   Consequently, capture eta_C>=24/25 throughout the frozen I_L would
   imply A_L>32041/90000>.34, even allowing coherent interference from
   every other component. This 96% capture is NOT established. Bulk
   velocity compatibility alone does not give it. The detector's finite
   ray margin is at least L/120-1, but its scattering/capture remainder
   remains unbounded at the needed 4% level.

2. **Retention.** Write the dimensionless low-energy fuel measure as
   a singular part of mass beta plus a density r(x), x=(E-Omega)/J.
   An oscillatory bound D(theta) for the density gives

       |G|<=omega+beta+D(theta).

   One sufficient route is a zero-extended BV density of variation V
   with omega+beta+V/16<=2/3. It would cover every theta>=16 and join
   the existing [2,16] theorem. Neither the singular-mass bound nor the
   density regularity/variation is currently proved. This route is not
   necessary; failure of its hypotheses would not refute operation.

## A precise negative result about the current scalar evidence

[SCALAR_LIMITATION.md](SCALAR_LIMITATION.md) constructs ONE fixed atomic
measure with the same fuel mean and variance, no weight at/below the
reference band bottom, all weight inside the first-band energy interval,
central weight 51133/59325 (about 86.19%), and a continuous early return
bound below .306 on 2<=theta<=16. Yet its fuel return is EXACTLY one at
theta=128*pi, inside the frozen [384,512] window for L=3072.

This refutes the implication from those scalar summaries to late
retention. At that recurrence every detector orthogonal to the fuel seed
also has zero arrival. **It is not the spectral measure of the actual h,
does not satisfy its full microscopic specification, and is not a B2
violating point.** The unchanged resolvent/primitive structure may still
exclude this behavior; an actual proof must use that extra information.

## Accounting and stop

The 34% condition remains ONE-source arrival. If proved, its unchanged
joint consequence is >11.458225%, not 33.8%. No source conditions have
been established on the late window by this packet. B2's failed eight-event
endpoint remains failed, and neither failure nor success of the scaling
family is inferred from the logical witness.

All 40 frozen controls passed, and the single replay produced byte-identical
output. Rational expressions are exact; ordinary-double two-dimensional
and trigonometric checks are diagnostics. Total child-worker wall time was
0.082285 seconds; peak child RSS was 10.125 MiB. No source dynamics, new
finite-size formation calculation, paid resources or external review ran.

This bounded cycle is closed. The next actual-source proof obligation is
the central-packet receiver kernel/capture and late fuel Fourier control,
or another controlled route to the same operating gate. The current
evidence does not license extrapolation across theta=16 to theta>=252.
