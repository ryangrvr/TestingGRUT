# Exact standard controls and hostile witnesses

**Issue #3; builder-authored; external review pending.** These are STANDARD
RESULT / REPRODUCTION / COUNTEREXAMPLE statements. No originality claim, new
physical law, GRUT candidate, selection score, information price, or card.
Numerical fixtures in `RESULTS.json` illustrate the proofs; they do not prove
the general statements. Source IDs refer to `SOURCES.json`.

## C01–C02. Probability, angle, and units

For a supplied complex qubit and a supplied Z measurement,
`psi = cos(theta/2)|0> + exp(i phi)sin(theta/2)|1>` gives
`P0 = cos²(theta/2)` and `P1 = sin²(theta/2)` by the Born rule. The formula
does not derive the state space, Born rule, dynamics, or physical meaning of
the basis. At theta=pi/2, phi=0 and phi=pi give identical Z probabilities but
orthogonal X eigenstates. Their X-plus probabilities are 1 and 0. Thus one
angular marginal does not determine the accessible quantum state or response.
Global phase and simultaneous rotation of state and measurement leave
probabilities unchanged. A Bloch direction becomes a spatial direction only
after an additional physical representation/calibration is given. [Q03,F01]

For a pure qubit the Fubini–Study line element is
`ds² = (dtheta² + sin²(theta)dphi²)/4`; the pure-state quantum Fisher metric
is four times this convention. This is geometric quantum mechanics, not
spacetime curvature or a gravitational field equation. [F01,F03]

In Q01's early reciprocal projection formula, the real unit vectors
`psi=(3/5,4/5)` and `q=(1,0)` yield `||psi||/(psi dot q)=5/3`.
This cannot be a cosine. Restricting to this real example also avoids any
ambiguity about conjugation of complex amplitudes.
With dimensionless theta, `omega` in T^-1, and physical x in L, Q02's
`(1/omega)partial_x theta` is T/L, not L/T. Its
`(1/omega)laplacian theta` is T/L², not the Newtonian potential's L²/T².
Q01's `theta''` is T^-2, not physical acceleration L/T².
Natural units or rescaled coordinates require an explicit scale/conversion;
these are objections to the stated physical identifications, not to all
possible alternative conventions. Q03's later `a=c theta'(-sin theta,cos
theta)` repairs the acceleration units and is a valid derivative of its
postulated velocity vector. The older error is not imputed to that revision.

## C03. The plus-sign determinant cannot have real-frequency zeros

Let H be a finite-dimensional invertible self-adjoint matrix, with real
nonzero eigenvalues lambda_j. Define the **actual Q03 section 171** finite
version by

```
Z_H(z;omega) = Tr[(omega²+H²)^(-z/2)]
log Det_*(omega²+H²) = -partial_z Z_H(0;omega).
```

The exponent -z/2 matters. Direct differentiation gives
`log Det_* = (1/2)sum_j log(omega²+lambda_j²)`, so

```
Det_*(omega²+H²)/Det_*(H²)
    = product_j (1+omega²/lambda_j²)^(1/2) >= 1
```

for every real omega. It is strictly positive. The ordinary determinant
uses exponent -z and has the same conclusion without the half power.
Q03 asks for this ratio to equal
`|xi(1/2+i omega)|²/|xi(1/2)|²`. That target has zeros at real ordinates of
critical-line zeta zeros; their existence does not require the Riemann
hypothesis. At any such ordinate the two sides cannot agree. [Q03,M01]

The exact fixture is D=diag(1,2), H=[[0,D],[D,0]], giving eigenvalues
{-2,-1,1,2}. Its Q03 ratio is
`(1+omega²)(1+omega²/4)`; the ordinary ratio is its square. No numerical
approximation to zeta is used.

**Scope of obstruction.** This disproves the stated finite spectral equality
and its finite discretization as an exact identification. An infinite
regularized determinant defined as the exponential of a finite regular
zeta derivative is also nonzero wherever that definition exists. Hence the
infinite claim requires either a failure of that regularity, a different
determinant/argument, or a different target at the zeros. This audit does not
prove impossibility of every infinite spectral construction, disprove the
Riemann hypothesis, or reject a corrected Hilbert–Pólya approach. If H has a
zero mode, the stated denominator is undefined; it is not a successful case.

## C04. Self-adjointness does not supply a thermal trace or vacuum residual

For a densely defined closed D, the standard block operator
`H=[[0,D*],[D,0]]` on `dom(D) direct-sum dom(D*)` is self-adjoint. Its grading
`Gamma=diag(I,-I)` satisfies `Gamma H Gamma=-H` on that domain. Consequently
the nonzero spectrum is paired under lambda -> -lambda. In every finite
block truncation preserving that structure, `Tr H=0` exactly. A residual
signed sum therefore requires a specified symmetry breaking, nonpaired
cutoff, or different observable. A signed operator trace also is not by
itself a derivation of physical vacuum stress energy. [Q03]

On l² of paired basis states, take the self-adjoint diagonal operator with
eigenvalues +/-n, n>=1, and its natural square-summability domain. For beta>0,
the partial thermal traces are `Z_N=2 sum_(n=1)^N cosh(beta n)>=exp(beta N)`.
Thus `Tr exp(-beta H)` is infinite. Self-adjointness and unitarity alone do
not establish Q03's Gibbs state; lower boundedness and trace-class/control
conditions are additional assumptions. A finite cutoff makes Z finite but
does not establish the unregulated physical theory.

## C05. Isometric encoding has a capacity condition

For U: H_bulk -> H_screen with `U*U=I`, U is injective and
`dim(H_bulk)<=dim(H_screen)`. If the source has dimension 3 and the screen
has dimension 2, rank(U*U)<=2, so U*U cannot be I3. The fixture that keeps
the first two coordinates gives diag(1,1,0), with squared Frobenius error 1.
An infinite bulk cannot be isometrically encoded into a finite-dimensional
screen. A restricted finite code subspace can be encoded; its selection is
then a premise requiring specification. Entangling bulk and screen in a
Schmidt state and tracing out bulk does not provide a pure screen isometry
for arbitrary bulk inputs. [Q03,F05]

## C06. Monotonic quantum information geometry is not unique

At the full-rank state rho=diag(p,q), p=1/4,q=3/4, take traceless Hermitian
tangent A with A12=A21=1 and diagonal entries zero. For a metric with
Morozova–Chentsov coefficient c, `g(A,A)=2c(p,q)`. The SLD coefficient is
`2/(p+q)`, giving 4. The Bogoliubov–Kubo–Mori coefficient is
`(log p-log q)/(p-q)`, giving `4 log 3`, which differs from 4. Both belong
to the standard contractive quantum metric family. Contractivity alone
does not select one mixed-state metric; a selected metric still does not
select a Hamiltonian, spacetime, or gauge content. [F01,F02]

## C07. A Hilbert space alone does not identify subsystems

In C² tensor C², let `|B>=(|00>+|11>)/sqrt(2)`. Its first marginal is I2/2,
with purity 1/2. A controlled-NOT unitary maps it to `|+> tensor |0>`, whose
first marginal is pure. Under a simultaneous conjugation of the state and
all observables, expectation values are unchanged. Reading the conjugated
observable algebras as an alternative factorization gives different
subsystem entanglement for the same abstract description. [F16]

This is **not** a local symmetry of an already physically fixed partition:
CNOT is an entangling global operation in the original factorization. The
result says the abstract Hilbert space does not fix that partition. Local
interactions, a local observable algebra, spatial structure, and accessible
operations can fix physically relevant subsystem relations. A GRUT claim of
unique identity must specify its invariants under such representation changes,
not equate arbitrary factorizations of a fixed laboratory to physically
identical intervention repertoires.

## C08–C09. Explicit alternatives delimit quantum reconstruction premises

The PR box is `p(a,b|x,y)=1/2` if `a xor b = xy`, and 0 otherwise, for four
binary variables. It is positive, normalized, and each local marginal is
1/2 independently of the other input. Its correlations are (1,1,1,-1), so
CHSH=4. Enumerating all deterministic local assignments gives CHSH<=2.
Quantum mechanics gives the Tsirelson bound 2 sqrt(2). No-signalling alone
therefore does not derive quantum probabilities. [F08,F09]

In real quantum theory, local symmetric 2x2 observables are spanned by I,X,Z.
The matrix Y tensor Y is real symmetric and squares to I4. Set
`rho_+/-=(I4 +/- Y tensor Y)/4`. Each has trace 1 and squares to rho/2, hence
has nonnegative eigenvalues 0 and 1/2. All nine expectations of products of
I,X,Z agree, but the global Y tensor Y expectation is +1 or -1. Real quantum
theory is not locally tomographic for these systems. The unnormalized real
state-space dimensions are 3 locally and 10 jointly, exceeding 3*3; complex
dimensions are 4 and 16=4*4. These are standard alternatives. Local tomography
excludes this real theory but does not by itself exclude classical theory.
Purification excludes a strictly classical simplex theory with ordinary
classical composites, because a pure joint classical state has pure marginals;
a Bell state supplies a mixed quantum marginal. A full reconstruction needs
its other explicitly stated operational axioms as well. [F07]

## C10. Reproduce the existing kinetic control; do not expand it

For labelled states s=-1,+1 define rates
`q_(-,+)=exp(g(u)+u)/2`, `q_(+,-)=exp(g(u)-u)/2` and
`pi(s)=exp(us)/(2 cosh u)`. Detailed balance follows because each stationary
flux is `exp(g(u))/(4 cosh u)`. The gap is `exp(g(u))cosh u`. With the matrix
J having every row pi,
`P(t,u)=J+exp[-t exp(g(u))cosh u](I-J)`.

The two supplied conventional models g=0 and g=u² have the same stationary
family and the same unforced generator at u=0, yet differ at u=t=1. This is a
minimal reproduction of the supplied theorem package's two-state control,
not a new family or a fundamental activity selector. The ambiguity disappears
when the microscopic transition law is specified; it is about incomplete
input data, not evidence that nature lacks kinetics. Analytic finite-jet and
smooth infinite-jet theorems are reported from the supplied record, not newly
reproduced in this fixture. The supplied three-state isospectral extension
was not independently rerun here. [G06]

## C11. A genuine selector, with its premises visible

Assume SU(3)xSU(2)xU(1), one generation of the SM representations, one Higgs
doublet, no right-handed neutrino, and all charged-fermion Yukawa couplings.
Write LH-doublet hypercharges q,l; RH-singlet charges u,d,e; Higgs charge h.
Yukawa invariance requires `u=q+h`, `d=q-h`, `e=l-h`.
The SU(2)²U(1) anomaly requires `3q+l=0`; the mixed gravitational anomaly
requires `6q-3u-3d+2l-e=l+h=0`. Thus

```
(q,u,d,l,e,h) = q*(1,4,-2,-3,-6,3).
```

The SU(3)²U(1) and U(1)³ anomaly equations then vanish identically. For q=1/6
the usual hypercharges result. This is a conditional loss of charge-ratio
freedom: after Yukawa relations only q,l,h remain; these two independent
linear constraints leave one normalization. It does **not** determine the
representations, number of generations, Higgs selection, gauge coupling,
masses, or charge normalization apart from its conventional absorption into
the coupling. Adding a right-handed neutrino changes the admissible anomaly
problem and permits additional directions such as B-L. This audit does not
claim a universal charge uniqueness theorem. [F17]

## C12–C13. Prime labels and projection do not establish new physics

The first four primes multiply to 210, not 30. The units modulo 30 are
{1,7,11,13,17,19,23,29}; their complement occupies 22/30=11/15 of residue
classes. The composite 49 is itself a unit modulo 30. A residue-class count
does not specify a stress-energy density, dark-sector particle population,
interaction, or cosmological mass fraction. [Q05]

For the completely positive trace-preserving X-dephasing channel
`E_X(rho)=(rho+X rho X)/2`, the symmetric |+> state is fixed, whereas
`|0><0|` becomes I2/2. This reverses a claimed basis-independent protection
of asymmetry over symmetry. Protection depends on the supplied noise and
coupling basis, not an integer's prime label. It does not preclude a specified
microscopic noise model that protects a chosen state.

Q04 simultaneously defines `vx=cos(omega t)`, `vy=sin(omega t)` and
`v²=vx²+vy²`, then identifies `v=omega` at c=1. The first equations require
v²=1 for every omega; omega=1/2 gives the exact contradiction 1=1/4. If the
vector denotes a unit clock motion and omega denotes a different massive
velocity ratio, these must be distinct variables with a further derived
map; the displayed equations do not provide that map.
