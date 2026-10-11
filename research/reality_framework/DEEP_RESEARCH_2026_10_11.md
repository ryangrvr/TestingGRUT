# GRUT: research toward a predictive synthesis

Date: 2026-10-11. Scope: owner's requested research pass while new notes are
pending. No CP-07 action is proposed or evaluated; CP-06 and earlier packets
remain preserved. These are source-backed findings and analytical synthesis,
not independent validation or an established theory of reality.

## Result of the pass

The strongest near-term construction route is a **single physical account of
gravitational response, quantum fluctuations and relational records**. Use
quantum gravitational EFT as its correspondence baseline, an influence
functional to calculate reduced dynamics, and physical gravitationally dressed
probes to define what is observed. This combination can make several previously
separate calculations mutually accountable without inventing another scalar
potential. It does not itself explain the microscopic origin of gravity.

The strongest specific microscopic alternative found here is the 2025 quantum
entropic-force construction of Carney and collaborators. It actually calculates
a force and its associated noise. Its unresolved relativistic completion is
where a substantive contribution could occur; duplicating its Newtonian limit
would not be a discovery.

These are my assessments, not conclusions asserted by the cited papers. The
search was focused rather than exhaustive. No evidence found here selects
a unique ontology, establishes consciousness-dependent dynamics, or completes
the quantum-gravity problem.

## 1. Start with physical subsystems, not coordinate-labelled boxes

Donnelly and Giddings explicitly construct dressed observables and prove an
obstruction to strictly local subsystem algebras in perturbative gravity.
A compactly supported operator carrying nontrivial Poincare charges needs
asymptotic metric dependence once made diffeomorphism invariant. This already
occurs at leading gravitational order. [1,2]

This is a direct reason to revisit GRUT's record/access vocabulary: observer,
source and environment cannot simply be declared exact tensor factors of the
physical gravitational Hilbert space. Approximate operational subsystems may
still be useful. Long-range dressing is not a license for superluminal messages.
The result has an asymptotically flat, perturbative domain and does not supply
an arbitrary cosmological observer algebra.

**Construction implication (our inference):** define a record as an outcome
of a specified relational or dressed detector observable. A coordinate value
of h_ab is not, by itself, an experimental record. A source-probe calculation
should state reference systems, physical dressing and which gravitational
degrees of freedom remain accessible. This can correct the framework without
assuming consciousness has a special physical status.

## 2. Response, memory and fluctuations belong in one calculation

Feynman-Vernon reduction expresses the effect of ignored quantum degrees of
freedom through an influence functional. For a Gaussian environment and a
linear source coupling, its two-point data determine both retarded response
and noise. The relation is exact within that Gaussian specification; an
interacting environment requires higher connected correlations. [3]

With H_int=-j B, define

    D_R(t,t')=(i/hbar) theta(t-t') <[B(t),B(t')]>,
    N(t,t')=(1/2)<{delta B(t),delta B(t')}>.

For two prescribed source histories j_plus and j_minus, set
j_delta=j_plus-j_minus. The Gaussian influence magnitude is

    |F[j_plus,j_minus]| = exp[-(1/(2 hbar^2))
                              integral dt dt' j_delta(t) N(t,t') j_delta(t')].

This is a visibility factor for the specified history experiment, not a
general solution of the measurement-outcome problem. Retarded response uses
D_R and the prepared state; the real influence phase also depends on both
histories. Positivity applies to N on physical, suitably smeared probes.

Hu and Verdaguer's stochastic-gravity framework adds stress fluctuations to
semiclassical back-reaction through an Einstein-Langevin equation. Its metric
correlations contain both initial gravitational fluctuations and matter-induced
fluctuations. It does not replace the complete quantum theory with classical
noise, and its leading treatment omits higher gravitational correlations. [4]

**Construction implication (our inference):** do not fit a gravitational
memory kernel, detector noise and decoherence separately. Derive them from
one action, state and physical observable specification. Disagreement between
those outputs would then be a meaningful failure of the assembled model.
Mean response alone cannot fix its quantum state or all fluctuation data.

## 3. There is a relevant advance beyond the old Markov approximation

The 2013 Anastopoulos-Hu derivation starts from weak-field GR and QFT and
produces a gravitational decoherence master equation, while identifying gauge
pitfalls in interpretations of coordinate fluctuations. Its one-particle,
nonrelativistic, Markov treatment gives energy-basis decoherence and depends
on environmental data, not only a Planck scale. [5]

A 2025 paper, Non-Markovian quantum master and Fokker-Planck equation for
gravitational systems and gravitational decoherence, extends related analysis
to a two-mass graviton environment at arbitrary temperatures. Its abstract
reports qualitatively different low-temperature temporal behavior from the
Markov case. This is a relevant comparator for any vacuum-memory construction;
a single relaxation rate cannot automatically substitute for it. [6]

The publisher abstract of [6], not its full derivation, was inspected in this
pass. Its predictions cannot be imported into a different apparatus or state
without matching its assumptions. No new convergence test of the old Gaussian
radial bubble is warranted.

**Construction implication (our inference):** retain the prepared gravitational
state, initial correlations and finite-time dependence before using a Markov
or thermal approximation. A low-temperature response may have no single
decay time, and reduced decoherence need not mean fundamental irreversibility.

## 4. A real microscopic entropic-force model exists

Carney et al. (PRX 2025) specify masses, quantum mediators and a thermalizing
environment. Their spin/oscillator constructions produce a Newtonian mean
force and calculate noise/decoherence. They provide distinct model versions,
not a full relativistic gravity theory. In the nonlocal version a parameter
limit can suppress decoherence and reproduce ordinary perturbative-gravity
behavior within the studied nonrelativistic regime. They do not establish
a bound excluding that limit. [7]

This is more useful than an unspecified vacuum bath: dynamics and quantum
matter are present. Its free parameters, thermal setting and spatial coupling
remain physical assumptions. A label such as "local model" does not establish
a relativistic causal completion. The paper leaves that problem open.

**Priority assessment:** investigate what its specified mediator coupling and
thermalization would require in a constraint-compatible relativistic theory.
The target is not simply recovering an inverse-square force; it is recovering
gravitational response to all stress-energy, the physical tensor sector and
consistent quantum fluctuations from the same mechanism.

### Analytical diagnostic: do not mistake thermal stiffness for fixed noise

We can make an independent equilibrium calculation before choosing a new law.
Consider commuting independent occupations n_a=0,1 with
H(r)=sum_a epsilon_a(r)n_a, held at a specified temperature T. Then

    F(r)=-T sum_a log[1+exp(-epsilon_a/T)],
    p_a=1/[1+exp(epsilon_a/T)].

Direct differentiation gives

    F'=sum_a p_a epsilon_a',
    F''=sum_a p_a epsilon_a''
          -(1/T)sum_a p_a(1-p_a)(epsilon_a')^2.

The force operator is -H', so

    Var(force)=sum_a p_a(1-p_a)(epsilon_a')^2,
    F''=<H''>-Var(force)/T.

This identity is ordinary equilibrium physics, not a GRUT law or a universal
gravitational fluctuation theorem. It shows that mean force and stiffness
need not fix force variance unless the microscopic curvature term <H''> is
also specified. Noncommuting Hamiltonians replace the simple variance with
an imaginary-time correlation integral; dynamic noise additionally requires
the relaxation spectrum.

Even if a model matches F=-G m1 m2/r, it must account for <H''>, its state
and its force spectrum before claiming a parameter-free noise prediction.
This gives an immediate test of whether a proposed synthesis actually links
observables or merely fits them separately.

## 5. Keep the two competing coherence pictures explicit

Galley, Giacomini and Selby's reversibility theorem constrains a classical
gravitational sector coupled to fully nonclassical matter with back-reaction.
At least one declared property must fail; it does not determine which nature
chooses. [8]

Oppenheim and collaborators derive a decoherence-diffusion trade-off using
complete positivity in their classical-quantum dynamical framework. It is a
comparison restriction, not an arbitrary non-Markovian quantum-gravity law.
The distinction between classical stochastic back-reaction and quantum
reduction matters even when mean forces agree. [9]

**Construction implication (our inference):** choose and disclose the status
of gravitational variables. A quantum state traced over inaccessible modes,
a classical phase-space distribution coupled to quantum matter, and a
semiclassical expectation-value equation are different theories. Combining
their attractive equations without a common state/evolution law is not a
consistent synthesis. Compare joint response, coherence and noise, rather
than treating a matched mean force as a full equivalence.

## 6. Entanglement offers a bridge, not a completed origin story

Jacobson's entanglement-equilibrium argument supplies a conditional link to
Einstein dynamics; conformal and nonconformal matter have different obligations.
Casini, Galante and Myers compute relevant-operator effects that challenge
a key nonconformal assumption. Thus the simple small-ball argument should
not be asserted for the complete massive Standard Model without further
analysis. [10,11]

Takayanagi's 2025 essay describes holographic connections between quantum
information and spacetime and identifies emergence of time and extension
to realistic cosmologies as open targets. It is an expert research outlook,
not experimental evidence that our universe is a tensor network. [12]

**Priority assessment:** holography is a serious longer-term microscopic lead;
an arbitrary finite code or lattice that reproduces an area law does not
earn Lorentzian spacetime, universal coupling or observed cosmology. Its
explanatory power must be tested through a specified duality and observables.

## Construction priorities from this pass

| Priority | Physical work | What would advance GRUT |
|---|---|---|
| 1 | Constraint-compatible influence calculation with dressed/relational probes | One state and dynamics predict response, fluctuations and physical records |
| 2 | Relativistic completion of a specified microscopic entropic mechanism | A derived tensor response and universal stress coupling, with linked noise |
| 3 | Entanglement/holographic microscopic dynamics | Earned causal geometry and measurable correspondence beyond an area-law toy |

This is a ranked assessment, not three authorized candidate attempts. New
notes may change the ranking. Nothing here reinstates the closed Leiden
apparatus investigation, consumes Card 3 or silently opens another CP law.

The synthesis should be judged by common physical dynamics and predictions,
not by whether it removes every input. GRUT can include measured constants
and established quantum structures. Its additional content must be specified
and tested. No consciousness premise is required to define physical records;
any future consciousness hypothesis needs its own operational evidence.

## Sources and inspection scope

1. Donnelly & Giddings, Diffeomorphism-invariant observables and their nonlocal
   algebra, PRD 93 024030 (2016), https://arxiv.org/abs/1507.07921.
   Publisher abstract and arXiv record inspected; erratum noted by publisher.
2. Donnelly & Giddings, Observables, gravitational dressing, and obstructions
   to locality and subsystems, PRD 94 104038 (2016),
   https://arxiv.org/pdf/1607.01025. PDF theorem and proof assumptions inspected.
3. Feynman & Vernon, Annals of Physics 24 118-173 (1963),
   https://authors.library.caltech.edu/records/x858p-sg795.
   Author-repository abstract of the reprint inspected; Gaussian formulas here
   use explicitly stated conventions, not a claim to audit the original paper.
4. Hu & Verdaguer, Stochastic Gravity: Theory and Applications,
   https://arxiv.org/abs/0802.0658. Full text functional approach, noise and
   initial/induced-correlation sections inspected.
5. Anastopoulos & Hu, https://arxiv.org/abs/1305.5231.
   Full text constraints, reduced dynamics and approximation discussion inspected.
6. Non-Markovian quantum master and Fokker-Planck equation for gravitational
   systems and gravitational decoherence, PRD (2025),
   https://doi.org/10.1103/b4n2-426h. Publisher abstract inspected only.
7. Carney, Karydas, Scharnhorst, Singh & Taylor, On the Quantum Mechanics of
   Entropic Forces, PRX 15 031038 (2025), https://arxiv.org/pdf/2502.17575.
   PDF microscopic models, free-energy expansion, quantum reduction and
   limitations inspected. Publisher DOI https://doi.org/10.1103/y7sy-3by1.
8. Galley, Giacomini & Selby, Quantum 7 1142 (2023),
   https://quantum-journal.org/papers/q-2023-10-16-1142/.
   Publisher theorem summary and assumptions inspected.
9. Oppenheim et al., Nature Communications 14 7910 (2023),
   https://www.nature.com/articles/s41467-023-43348-2.
   Full text class, positivity and trade-off discussion inspected.
10. Jacobson, https://arxiv.org/abs/1505.04753; PDF examined in preceding pass.
11. Casini, Galante & Myers, JHEP 03 (2016) 194,
    https://doi.org/10.1007/JHEP03(2016)194; publisher and author-repository
    abstract inspected. Detailed nonconformal calculations not independently rerun.
12. Takayanagi, PRL 134 240001 (2025), https://doi.org/10.1103/pg4r-fy8n;
    publisher abstract and research outlook inspected.

No new experimental claims, complete literature coverage, independently checked
derivations or successful fundamental selection are asserted by this pass.
