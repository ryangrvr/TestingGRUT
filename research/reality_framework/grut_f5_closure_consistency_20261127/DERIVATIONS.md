# Closure obstruction and exact benchmarks

## 1. Arbitrary sampled laws survive the consistency requirement

**Proposition.** Let $\mathcal K$ be a matter Hilbert space and let $I,W_1,\ldots,W_L$ be any finite list of unitary prefixes on it. They may represent successive products of arbitrary unitary step gates. In the class frozen for F5 there is a positive physical state space with two good ideal clock charts, a state-dependent clock rate and a nonseparable clock interaction, containing this sequence in the sector $r=0$. The initial matter density in that sector is arbitrary. Every joint physical event, including events from imperfect clock POVMs and retained records, has the same joint and conditional probabilities in the two complete descriptions.

**Construction.** Choose an odd prime $N>L$ with $N\ge5$. For both $s$, set $V_{s,0,j}=W_j$ at the required ticks and fill any unused prefixes arbitrarily, keeping the specified wrap to $V_{s,0,0}=I$. In the $r\ne0$ sectors keep $\mathcal K$ unchanged and set only the prefix at $r=1,a=1$ to $-I$. Then $A=j$ in $r=0$ yields $W_j\rho_{\mathcal K}W_j^\dagger$ after discarding the remaining degrees of freedom if the original seed factor was $\rho_{\mathcal K}$. Choose it freely, pure or mixed. The clock $\nu_s$ and the interaction in the other relative sector are independent of the target prefixes.

At fixed $s=0$, strip the common one-step translations of both clocks. Use rows $a=0,2$ and columns $b=1,3$. The four relative labels on this subarray are $1,3,4,1$ for $N=5$, and more generally $1,3,N-1,1$ for $N\ge5$. All are outside $r=0$. The residual phase is $-1$ only at the first entry. The sign matrix
\[
\begin{pmatrix}-1&1\\ 1&1\end{pmatrix}
\]
has determinant $-2$, ruling out a product of separate clock phases in the declared factorization even with matter held fixed. In the full shift, the rate of B depends on matter $S$ and the relative clock sector controls the matter gate table. The coupling witness is a property of the specified clock factors and effects. A simultaneous unitary transformation of all factors and their observables can change its representation.

The [physical construction](MATHEMATICAL_FOUNDATIONS.md) gives two unitary reductions and carries each physical effect $F_{\mathcal L}$ to the corresponding chart:
\[
\mathrm{tr}\,[U_A\rho U_A^\dagger\, U_A F_{\mathcal L} U_A^\dagger]
=\mathrm{tr}\,[\rho F_{\mathcal L}]
=\mathrm{tr}\,[U_B\rho U_B^\dagger\, U_B F_{\mathcal L} U_B^\dagger].
\]
Use this twice for a conditional event, once with its joint effect and once with its filter. This holds for every supplied unitary history and preparation. With the event algebra unrestricted, changing the prefixes can be represented as a unitary change of the kinematic description of a shift. The counterexamples hold the physical clock and matter readout interface fixed, so the compressed effects of the same bare record procedure are generally different across different laws. The no-write comparison below witnesses an operational difference between laws under that fixed interface. It also shows that agreement between descriptions of each law survives that difference.

**Scope.** This is a constructive obstruction to using matched internal-clock probabilities as a law- or state-selection rule. It applies to arbitrary finite samples and arbitrary seed density in the admitted class. A finite history can sample a Hamiltonian at any finite list of ticks, with its cyclic end supplied separately. The theorem does not require an unmodified Hamiltonian to continue around the cycle, a common microscopic locality prescription, or a fixed complexity or energy budget as $N$ grows. The proof also allows separable infinite-dimensional $\mathcal K$ with unitary prefixes. Only the reference degrees and the number of sampled ticks have to be finite. The absence of external time therefore does not itself require finite Hilbert-space dimension within this framework.

## 2. A quantitative same-event criterion

Transport both descriptions to a common physical Hilbert space. Suppose two procedures are specified by effect pairs $0\le E_i\le D_i\le I$, where $D_i$ is the condition and $E_i$ is the joint event and condition. For $p_i(\rho)=\mathrm{tr}(\rho E_i)/\mathrm{tr}(\rho D_i)$ with both denominator probabilities at least $d_*>0$,
\[
|p_A(\rho)-p_B(\rho)|\le
\min\!\left(1,\frac{\|E_A-E_B\|+\|D_A-D_B\|}{d_*}\right).
\]
To obtain the second term, use a common denominator $d_A=\mathrm{tr}(\rho D_A)$ in the difference and subtract the term $p_B\,\mathrm{tr}[\rho(D_A-D_B)]/d_A$. Born positivity gives $0\le p_B\le1$. Each expectation difference is bounded by the operator norm. The inequality quantifies the effect of measurement calibration and rare conditioning. It is a consequence of the ordinary probability rule.

There is also an exact finite-dimensional converse at the level of a binary conditional probability. For nonzero $D_i$, suppose the conditional functions agree on all full-rank density matrices and are nonconstant. Choose $\rho,\sigma$ with different values of the common conditional function. Evaluate the two functions on a convex mixture and cross-multiply. The product of the nonzero difference of conditional values and
\[
\mathrm{tr}(\rho D_A)\,\mathrm{tr}(\sigma D_B)
-\mathrm{tr}(\sigma D_A)\,\mathrm{tr}(\rho D_B)
\]
vanishes. It follows that the ratio of the two filter expectations is equal on $\rho,\sigma$. Every other full-rank state has a conditional value different from at least one of these, so that ratio is one positive constant $\lambda$ on the full-rank density matrices. Their span then gives $D_B=\lambda D_A$, and agreement gives $E_B=\lambda E_A$. This characterization is of the conditional probability functions. A common factor changes the overall rate of the condition, and identical probabilities alone do not identify the post-measurement instruments. Constant conditional functions can agree even for inequivalent filters.

## 3. The registered history and physical backreaction

With seed $|+\rangle_S|r=0\rangle|00\rangle_{DR}$, write $d,q$ for the data and record bits. The supplied source is the same for both $s$ components. Explicit prefixes give

| $a$ | $b$ at $s=0$ | $b$ at $s=1$ | $(d,q)$ |
| ---: | ---: | ---: | :---: |
| 0 | 0 | 0 | 00 |
| 1 | 1 | 2 | 10 |
| 2 | 2 | 4 | 11 |
| 3 | 3 | 1 | 01 |
| 4 | 4 | 3 | 01 |

The write at $a=2$ copies the current bit into $R$. The bit stays written when the data changes at $a=3$, through $a=4$. The inverse wrap returns the two bits to their seed values. A stationary global state therefore has relational records with a finite operational history. Replacing the seed by $I_{\mathcal L}/40$ makes every complete ideal slice maximally mixed, giving record probability $1/2$ at all $a$. A blank record and a specific developmental sequence are preparations. They do not follow from imposing stationarity.

At $a=0$ both B labels coincide, leaving $\langle X_S\rangle=1$. At $a=2$ the bare B reading distinguishes $s=0$ from $s=1$. Tracing out B then sets $\langle X_S\rangle=0$. The two orthogonal $S$ preparations $|+\rangle,|-\rangle$ have trace distance one at full access. In the A slice at $a=2$ both give the same mixture on $S,D,R$ if B is inaccessible. They also give the same mixture on $S,D,R$ for a noisy A report of 2, since its second raw tick $a=1$ has distinct B labels on the two branches as well. These limits are exact losses of operational distinguishability under access restriction. With the other clock retained the ideal chart maps preserve it.

## 4. Physical probabilities and clock changes

The record probabilities at ideal A and B readings of 2 are respectively $1$ and $1/2$. B at 2 selects $(s,a)=(0,2),(1,1)$ with equal weights. With the frozen one-tick reporting errors $e_A=1/4$, $e_B=1/3$, the separate reports of 2 give $3/4$ (A) and $1/2$ (B). These are different filters on the cyclic orbit and need not select the same record. They compare procedures using only one clock report at a time.

For the joint noisy report $\alpha=\beta=2$ the three contributing rows are $(s,a)=(0,2),(0,1),(1,1)$. Each had prior weight $1/10$. Including the error weights gives $1/20,1/120,1/60$. The record is one only on the first row. The filter probability is $3/40$, the joint filter-and-record probability is $1/20$, and the conditional probability is $2/3$. Each is the same whether the exact A or B reduction represents the event. The joint filter is different from the product $1/25$ of marginal filter probabilities. Omitting the writer with identical clock settings and updating the inverse wrap by the same rule leaves the filter probability unchanged and makes the record probability zero. That model also satisfies matched-event consistency.

The phase-probe control uses the $r=1$ plus seed and retains the full Hilbert space of the other registers. Its physical projector is the projector onto that seed. At A anchor 2 the supported basis states of the remaining registers have B readings 3 and 0 for $s=0,1$ and both amplitudes are positive. At B anchor 2 the A readings for those branches are 1 and 3, giving signs $-$ and $+$. The transported effect gives probability one in both descriptions. Erasing the relative sign in the B effect while retaining the true state makes the overlap exactly zero. Thus matching only the diagonal clock populations would miss a physical comparison.

## 5. Calibration and boundary controls

If the B stride is changed to zero for $s=1$, its raw reading equals $r$ along that branch. A fixed B section sees only one of the five relative sectors there and all five on the $s=0$ branch. Its rank is $4+20=24$, compared with 40 for a complete ideal chart. A remains a good clock and the physical probabilities of actual joint POVMs remain consistent. A B-based coordinate reduction here cannot describe every physical preparation.

The reciprocal wrap cancels the net supplied product of matter gates. In a more general cyclic constraint with an arbitrary net unitary $W$ on the matter sector, an invariant history requires its seed to satisfy $W|\psi\rangle=|\psi\rangle$. In the frozen sign-flipped wrap control, $T'^5=-I$, and hence $\ker(T'-I)=0$. Its gauge group has order ten and the original five-step twirl is no longer the projector. This control identifies a dependence on the chosen boundary and gauge sector. A different gauge character or wrap phase would be other supplied data. There is no universal boundary condition selecting physical preparations in the proposition.
