# One bounded interval factorization correction

The first run enclosed the predictor residuals and a negative radial trial for all three charges. Direct interval block LDL widened until the pivot determinant interval included zero at block 123 for every registered shift and charge. Thus the original run has a method failure and gives no continuum index. The raw result and code are retained in the preceding commit. [Initial result](INITIAL_ATTEMPT.json)

Freeze this amendment before the sole corrected production invocation. Keep exactly the same profiles, mesh, exterior, shifts, acceptance criteria and arithmetic precision. Keep the original 180-second wall cap and five spectral counts per charge for this invocation. The direct propagated interval pivots are superseded by a test of a point factorization with outward residual and inverse bounds.

Construct a point unit block bidiagonal triangular matrix $T$ and a symmetric block diagonal $D$. Their stored multiprecision values define exact dyadic point factors. Let $\widehat A=TDT^{\mathsf T}$. Enclose the pivot inertia of $D$ and the block tridiagonal residual of the original nodal matrix $A-\widehat A$ independently with outward intervals. Bound $\|D^{-1}\|_\infty$ from the small pivots. For $T^{-1}$ bound the vector of absolute row sums by the recurrence $r_i=\mathbf 1+|T_{i,i-1}|r_{i-1}$. Applying the reverse transpose recurrence bounds its absolute column sums. Hence
\[
\|\widehat A^{-1}\|_\infty\le
\|T^{-1}\|_1\|D^{-1}\|_\infty\|T^{-1}\|_\infty.
\]
Require the product of this upper bound and the interval residual norm upper bound to be strictly below one. The symmetric path $\widehat A+t(A-\widehat A)$ is then invertible for every $0\le t\le 1$ by the Neumann criterion, so the exact nodal matrix has the inertia of $D$. The continuum nodal-comparison proof and every physical gate are unchanged. Failure of the corrected test ends this campaign at its achieved scope.
