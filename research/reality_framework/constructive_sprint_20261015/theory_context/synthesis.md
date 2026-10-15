# A constructive architecture for the strongest presently supported GRUT

Author synthesis, 15 October 2026 UTC. The frozen CP-12 law and earlier results
are preserved. This note supplies a conditional locality bridge and a coherent
physical interpretation; it does not declare a completed theory of nature.

## 1. The working mathematical construction

The strongest supported microscopic construction is a quantum theory of
connection resources, compact gauge transport and fermionic matter. Its physical
space is the exact Gauss-law subspace

\[
\mathcal H_{N,Q}^{\rm phys}
=\bigoplus_{\substack{g\subset K_N\\|g|=Q}}
\{\psi_g:G_i\psi_g=G_e\psi_g=0\},
\]
\[
G_i=\sum_{e\ni i}E_{ie}-n_{a_i},\qquad
G_e=-E_{ie}-E_{je}-n_{b_e}+1.
\]

Inactive connections have filled edge fermions and zero fields. Active
connections carry two integer electric fields with \([E,U]=U\). The frozen law is

\[
H=\kappa Q+2mN_v+\frac{\epsilon_E}{2}\sum_{e,i\in e}E_{ie}^2
+t\sum_{e=(i,j)}n_e[(a_i^\dagger U_{ie}-a_j^\dagger U_{je})b_e+\mathrm{h.c.}]
+\frac{\Gamma}{M}\sum_r L_r^\dagger L_r,
\quad M=\binom N2.
\]

Here \(\kappa=1/2,m=t=1,\epsilon_E=1/4,\Gamma=1/10\). A neutral removable
edge can exchange with any absent edge. The supplied initial ray is the coherent
equal superposition of all Q-edge graphs with empty vertices, filled edge
fermions and zero fields. No spatial host is encoded in this ray.

For finite N,Q the diagonal electric operator is self-adjoint and bounded below
on its usual square-summable integer-rotor domain. Hopping and rewiring are
bounded perturbations: the coarse bounds \(\|H_{\rm hop}\|\le2|t|Q\) and
\(\|R\|\le2\Gamma Q\) suffice. Thus H has a well-defined self-adjoint,
bounded-below realization, with exact unitary evolution. Q, total fermion number,
Gauss constraints, energy and vertex-permutation symmetry are conserved. This
is a consistent finite-resource quantum model, rather than a potential added to
an unspecified substrate.

Complex Hilbert space, unitary evolution, tensor composition, Fermi statistics,
compact U(1), the Born rule and external Hamiltonian time are inputs. None follows
merely from writing Gauss constraints. They cannot be advertised as explained
quantum rules. N,Q and the coherent initial preparation are additional supplied
physical resources; their cosmological origin is not derived.

## 2. One internal mechanism, with separate physical stages

The constructive mechanism is now specific:

1. Matter hopping creates vertex occupation, edge holes and electric flux while
   preserving Gauss law.
2. Nonzero flux or an edge hole makes that edge unavailable to the specified
   neutral rewiring operator.
3. Quantum amplitudes around a closed path generate gauge-invariant electric
   circulation and virtual Wilson-loop vertices. These amplitudes retain
   fermionic interference and are not classical counts of possible paths.
4. Persistent connection records would support a local Hamiltonian on the
   recorded network. The next section turns this last statement into a
   quantitative theorem, including the error from actual graph motion.

The earlier packets establish stages 1–3 in defined finite systems, together
with exact onset and resource constraints. Large-network formation remains a
different dynamical question. The a2f447d low-matter obstruction and charged
energy-matched comparators show why charge must be included in that question.
The present sprint advances access beyond K5: its independent charged-access
calculation sums signed amplitudes from the original N6,Q12 neutral ray over
all 455 graph parents and all 15 labelled four-regular octahedral target hosts,
and obtains

\[
p_{G_6}(T)=
\frac{28937809141133489}{173110482737718477580671905802392578125}
T^{52}+O(T^{54}).
\]

This is positive actual-ray access in a sector with graph choice, rather than
only energetic admissibility. Its target and exact integer controls are in the
companion charged-access packet. At N7,Q14, the same shortest-path formation
kernel distinguishes two exhaustively enumerated four-regular graph families:
complements of C7 and of disjoint C3∪C4. Their per-host coefficient ratio is
approximately 0.5979550202 in that order. Consequently the C7-complement fraction
within this very rare leading-order charged target is approximately 0.6721452827,
versus its uniform-graph-family prior 24/31. This is an internal constructive
topology discrimination during charged-network formation. It does not identify
a preferred dimension, an appreciable unconditional phase or a long-time state.

The earlier induced-loop packet is a useful genuine interaction bridge. In a
declared fixed-graph Feshbach compression its vacuum-to-one-flux matrix element
is \(K_\ell(a)t^{2\ell}+O_g(t^{2\ell+2})\), with positive minimum-path
coefficient and \(a=\Gamma(M-Q)/M\). This is one derived piece of an effective
gauge theory. Hopping-generated loop interactions are established lattice-gauge
physics; the CP-12 coefficients and their validity domain are model-specific
calculations [S6]. It is not yet an energy-independent, uniformly convergent
large-N magnetic Hamiltonian at frozen t=1, nor a photon or gravitational mode.

## 3. A new constructive bridge: local dynamics with a full-law error

Let \(P_g\) project onto one specified Q-edge graph, retaining **all** physical
matter and integer electric configurations on it. Its true compression is

\[
H_g=P_gHP_g=D_g+K_g+a\sum_{e\in g}r_e,
\quad a=\Gamma(M-Q)/M,
\]

where \(r_e\) projects onto the neutral removable state of edge e. The last
term is essential: removing off-diagonal rewirings while retaining their escape
diagonal leaves a physical local edge potential. Dropping it compares against a
different Hamiltonian and changes O(Γ) local physics.

Only R leaves the graph block. Different pairs \((e\in g,f\notin g)\)
lead to different graphs \(g-e+f\). Their outgoing ranges are orthogonal.
Therefore for every \(v\in P_g\mathcal H\), without an electric cutoff,

\[
\boxed{\|(1-P_g)HP_gv\|^2
=\left(\frac\Gamma M\right)^2(M-Q)
\sum_{e\in g}\|r_ev\|^2.}
\tag{1}
\]

Let \(v_\psi(s)=e^{-iH_gs}\psi\), with normalized
\(\psi\in P_g\mathcal H\). Duhamel's identity gives

\[
\|e^{-iHT}\psi-v_\psi(T)\|
\le\varepsilon_\psi(T):=
\int_0^{|T|}\|(1-P_g)HP_gv_\psi(s)\|\,ds
\le C_N|T|,
\]
\[
\boxed{C_N=\frac\Gamma M\sqrt{Q(M-Q)}.}\tag{2}
\]

For Q=νN at fixed ν, \(C_N=\Gamma\sqrt{2\nu/N}[1+O(N^{-1})]\).
Thus every supplied fixed host, with arbitrary allowed coherent matter/flux,
remains close to its **true local compression** for \(T=o(\sqrt N)\).
Flux protection can improve this state-dependent estimate; it is not needed for
the uniform bound. This attributes part of conditional host persistence directly
to the frozen 1/M graph kinetic normalization, rather than to flux alone.

For a normalized state with \(|E_{ie}|,|E_{je}|\ge L\) initially on every
present edge, the edge must undergo at least 2L of its own hopping operations
before becoming removable. The full two-leg hopping operator on that edge has
norm \(\sqrt2|t|\). All remaining terms in H_g preserve its fields. A unitary
Dyson remainder therefore gives the refinement

\[
\varepsilon_\psi(T)\le C_N\min\left\{|T|,
\frac{(\sqrt2|t|)^{2L}|T|^{2L+1}}{(2L+1)!}\right\}.
\tag{3}
\]

This is a supplied-state bound. Larger flux costs electric energy at least
\(\epsilon_EQ L^2\); it is not a free means of selecting a geometry from the
neutral ray. In particular, using L growing without bound at fixed extensive
available energy is not justified.

For a bounded Hermitian source B that preserves the graph and detector O, define
\(\chi_H=-i\langle[O_H(T),B]\rangle_\psi\), corresponding to the source
pulse \(e^{-i\theta B}\), and use H_g for \(\chi_g\).
The evolved vectors for both \(\psi\) and \(B\psi\) enter the response.
Writing \(\varepsilon_{B\psi}\) for the same residual integral with the
unnormalized source vector gives the exact error estimate

\[
\boxed{|\chi_H(T)-\chi_g(T)|
\le2\|O\|\big(\|B\psi\|\varepsilon_\psi(T)
+\varepsilon_{B\psi}(T)\big)
\le4\|O\|\|B\|C_N|T|.}\tag{4}
\]

The first inequality follows by adding and subtracting one cross matrix element
between \(U\psi\) and \(UB\psi\); it controls actual linear response, not
only graph survival. A source must satisfy the same protected-field assumption
to use (3) on Bψ. Vertex density sources do so. A detector-specific dual estimate
can improve (4); it cannot invalidate the uniform statement.

### Why the compressed model is local

Associate one site to each vertex fermion and one site to each edge center,
including its edge fermion and two rotors. H_g acts on the barycentric graph of
g: hopping only joins a vertex with one adjacent edge center; electric terms,
mass terms and the retained escape diagonal are onsite. For maximum graph degree
bounded independently of N, hopping interactions have uniformly bounded norm.
The unbounded electric terms are onsite; passing to their interaction picture
leaves the interaction norms bounded. Standard finite-range Lieb–Robinson
arguments then bound commutators of bounded, gauge-invariant even local observables
outside a finite-velocity graph cone [S1,S2,S10]. The bound can be established on
the enlarged local rotor/CAR algebra and restricted to the invariant physical
Gauss subspace; a factorization of the constrained physical space is not needed.
S10 supplies the unbounded-onsite version, including the required strong-topology
interaction-picture argument. No relativistic Lorentz symmetry is implied by a
Lieb–Robinson upper bound.

Combining that local cone with (4) produces a full CP-12 response bound consisting
of the usual distance-decaying local term **plus a controlled graph-motion
floor** of order \(|T|/\sqrt N\). This is an explicit route from a persistent
network to operational approximate locality under the unchanged Hamiltonian.
It neither assumes nor proves that the original unfiltered ray generates that
network.

If a supplied family has diameter O(N^(1/d)) and transports ballistically on its
local compression, (2) is sufficient for small graph-motion error over a crossing
time when d>2. This proves a useful conditional regime for d=3 and also every
larger d. The bound being inconclusive for d≤2 is **not** instability, exclusion
or dimension selection. Random bounded-degree hosts with logarithmic diameter
are also covered and must remain adversarial comparators.

### Scope of graph conditioning and records

Equations (1)–(4) are single-host results. For a coherent superposition of many
graphs, amplitudes leaving different starting graphs can reach the same final
graph and interfere. Summing the single-host probabilities or their rms bounds
does not retain the 1/√N gain. A full coherent state may instead be compared with
a declared direct-sum trial Hamiltonian using its **actual full residual**;
that residual can be extensive.

An actual graph/record measurement also has a resource cost. Copying orthogonal
graph records into a memory can entangle system and apparatus and remove
interference from the reduced system. It is not an energy-neutral physical
selection principle. Final readout conditioning alone is not proof that a
structure existed throughout a propagation history. Initial-distance and
final-distance detectors must be distinguished in the joint evolving state.

The companion `locality_bridge_check.py` checks (1), the escape diagonal and
the response comparison on every fixed-graph block of the raw regulated
N4,Q4 physical Hamiltonian. It is a finite numerical control of the algebra;
the unbounded-H statements above are analytical, not inferred by cutoff fitting.

## 4. How dimension would be earned

The candidate spatial state must have appreciable probability under the
unfiltered dynamics, persistent connected records, bounded local coordination,
and two compatible intrinsic dimensions across a size-growing interval:

\[
V_g(x,r)=|B_g(x,r)|\sim r^{d_H},\qquad
P_g(s)=N^{-1}\operatorname{Tr}e^{-s\Delta_g}\sim s^{-d_s/2}.
\]

The graph Laplacian here is only a structural diagnostic. The matter response
must independently exhibit propagation with distance on that same state; its
propagator cannot be replaced by an assumed Laplacian. A 3D interpretation needs
d_H≈d_s≈3 on ranges away from lattice and finite-size cutoffs, approximately
isotropic response, and state weight/persistence scaling. Four regularity and
N+1 independent cycles are compatible with supplied one-, two-, three-dimensional
and nonspatial networks. They do not fix dimension.

## 5. A conditional continuum completion with explicit matching obligations

The economical continuation is a **conditional infrared description**, not a
second microscopic law. If state-selected relational geometry, local transport,
a gap hierarchy and a common approximately Lorentzian cone are established,
integrating out modes can organize an IR action as

\[
S_{\rm IR}=\int d^4x\sqrt{-g}\left[
\frac{M_{\rm eff}^2}{2}R-\Lambda_{\rm eff}
-\frac{1}{4g_{\rm eff}^2}F_{\mu\nu}F^{\mu\nu}
+\mathcal L_{\rm matter}+\sum_i\frac{c_i\mathcal O_i}{\Lambda_*^{d_i-4}}
\right].
\]

Every symbol is a **matching target**, not a previously derived CP-12 output.
The dimensionality and Lorentzian form may be used here only after the proposed
limits have been established. If derived, the matching calculation would predict
coefficients and departures rather than freely appending Einstein gravity.
Known gravitational EFT is a controlled comparison at low energy; it is not a
GRUT discovery [S3].

There are two concrete gravity bridges, both conditional:

* **Induced action:** once a common covariant matter operator has been derived,
  its functional determinant contains curvature terms. For a real scalar,
  \(\Gamma_E=\tfrac12\operatorname{Tr}\log(-\nabla^2+m^2+\xi R)\).
  In four dimensions its proper-time expansion contains
  \((4\pi s)^{-2}\int\sqrt g[1+s(1/6-\xi)R+\cdots]\).
  Matching the R term generates an Einstein-action coefficient, with
  \(M_{\rm ind}^2=(1/6-\xi)(4\pi)^{-2}
  \int_{\Lambda^{-2}}^\infty ds\,s^{-2}e^{-m^2s}\) for this chosen regulator.
  Its value is regulator/matching sensitive and allowed bare counterterms
  remain; neither the physical Planck scale nor small cosmological constant is
  fixed by writing a cutoff. This is the Sakharov/heat-kernel mechanism [S4,S5].
* **Entanglement equilibrium:** if a Lorentzian continuum, local vacuum modular
  Hamiltonian and entropy-area coefficient \(\alpha\) have been derived,
  stationary generalized entropy of small balls gives
  \(\delta G_{00}=2\pi\delta\langle T_{00}\rangle/(\alpha\hbar)\).
  This is Jacobson's conditional result, with additional nonconformal-matter
  obligations; it does not generate its own geometric premises [S7].

Neither bridge is currently a derivation of CP-12 gravity. Compact U(1) Gauss
constraints are not the gravitational Hamiltonian/momentum constraints. Event
permutation multiplets are not Lorentz spin. A valid graviton claim additionally
needs a physical two-helicity massless tensor response and the appropriate Ward
identities. If a Lorentzian massless spin-two S-matrix exists, the soft-gauge
identity requires a common coupling to all matter [S8]. That conditional theorem
cannot manufacture the required tensor mode.

The Weinberg–Witten theorem must also be addressed precisely: an exactly
Lorentz-covariant conserved physical stress tensor with the relevant matrix
elements obstructs a massless spin-two carrier [S9]. CP-12 lacks microscopic
Lorentz invariance, so its premises are not automatically present. This is not
an automatic evasion in a proposed Lorentzian IR limit; the physical stress and
gravitational constraint algebra would have to be constructed.

The prior CP03 determinant model establishes a conditional large-N pole in a
supplied four-dimensional Lorentzian expansion with tuning. Its real classical
positive-Gram branch cannot be Lorentzian. No map from CP-12 physical states to
that model's fields, measure or constraints has been derived. CP01's scalar-
tensor action and CP02's embedding EFT are likewise independently supplied
continuum constructions, not CP-12 consequences. Combining their successful
formulas by naming the result GRUT would not join the mechanisms.

## 6. Matter, clocks and observers within this architecture

The current matter degrees of freedom are vertex and edge spinless fermions
with one compact gauge symmetry and a supplied background charge. This is
internally consistent but not the Standard Model. Relativistic spin, chiral
representations, anomaly cancellation, three families, nonabelian interactions,
Higgs physics, masses and observed couplings have not been derived. Adding
SU(3)×SU(2)×U(1) matter by hand would create a declared extension, not an
explanation of the observed inventory.

Time T in CP-12 is already present. Operational clocks can be modeled as finite
subsystems correlated with processes; that constructs readings rather than the
origin of external time. Earlier CP09 bounded-below apparatus calculations show
how an autonomous clock/record can approximate a source–detector response, but
that different Hamiltonian's apparatus cannot be inherited without a coupling
and energy/preparation calculation. A finite closed unitary model also does not
derive a fundamental entropy arrow: total pure-state entropy remains zero,
while reduced entropy and physical records can grow from low-correlation inputs.

Observers can consistently be finite matter/flux systems with accessible
records. A record's distinguishability and lifetime are ordinary dynamical
quantities, as are its disturbance and erasure costs. The conditional state
after a physical record is read is computed by a specified instrument; it is
not a consciousness-dependent modification of H. No present calculation derives
phenomenal consciousness, a single-outcome selection law or Born probabilities.
Their absence does not invalidate the model's ordinary conditional quantum
predictions, and no unsupported consciousness variable is added to rescue it.

The numerical couplings are specified model units, not measured constants of
nature. No physical edge length, second, Planck scale, or emergent value of ℏ
has been fixed. A proposed phenomenology must independently calibrate those
quantities and then hold out an observable relation for testing; changing units
cannot count as a new empirical prediction.

## 7. Comparison and constructive verdict

| Question | Actual construction | Status | What would advance it |
|---|---|---|---|
| Quantum mathematical rules | Complex physical Hilbert space and unitary H | Input | An operational reconstruction would still need its independent axioms |
| Matter–geometry interaction | Gauge-covariant hopping changes removability and graph response | Derived finite-model mechanism | Joint larger charged-network formation test |
| Closed gauge interaction | Sign-resolved virtual loop vertices | Derived in stated effective calculation | Uniform frozen-coupling matching and gapless gauge modes |
| Charged extensive access | Exact N6 actual-ray onset with graph choice; N7 topology-dependent formation kernel | Finite positive access and short-time bias; macroscopic weight open | Finite-time network weight and lifetime at growing size |
| Persistent local dynamics | Equations (1)–(4), including local escape potential | Exact conditional host theorem | Apply residual to an internally formed unfiltered network |
| Three spatial dimensions | Intrinsic graph/response tests defined | Hypothesis, not selected | Competing-dimensional state weights, growth and transport |
| Time/causal cone | External time; conditional graph cone plus leakage floor | Input time and conditional locality | Common operational clocks and approximate Lorentzian scaling |
| Einstein gravity | Explicit matching route through IR action or entropy equilibrium | Conditional continuation, not derived | Two tensor modes, constraints, common stress coupling and matching |
| Ordinary particles/forces | Spinless fermions and compact U(1) | Partial model, not SM | Derive a consistent low-energy spectrum and representations |
| Physical observers | Dynamical subsystem records and instruments | Conventional consistent account | Autonomous accessible stable records in formed geometry |
| Empirical novelty | No calibrated beyond-established-physics prediction yet | Open | Fix scales and compare a held-out measurable relation |

There is no mathematical need at this stage to alter the frozen Hamiltonian.
It has not been proved unable to make a charged relational spatial phase. Its
all-pairs rewiring, canonical graph entropy and preparation dependence remain
concrete challenges. A successor should be introduced only to address an
identified failure with an independently motivated new ingredient; inserting
a preferred lattice, dimension or plaquette reward merely to obtain the target
would not explain it. A local adjacency-mediated graph-move law would be a
possible separately labelled locality hypothesis, but would impose a new rule
and can preserve unwanted disconnected sectors. It is not adopted here.

The strongest working theory actually constructed is therefore a predictive
finite quantum graph–gauge model with a derived internal protection mechanism,
accessibility tests, and a rigorous conditional route to local dynamics.
The continuum/gravity/Standard-Model portion is an explicit matching program,
not a completed unified theory. The positive locality theorem is constructive
progress: it identifies exactly how an actual persistent network could become
physically consequential, and supplies a measurable error budget for that claim.

## 8. Discriminating tests and the next integrated experiment

The model makes falsifiable **model-specific** statements even before it makes
a new empirical prediction about nature. For t=0 the neutral ray only acquires
a phase: no matter, flux, protection or density response can form. At nonzero
t, hole-supported protected triangles and occupied-edge flux triangles have
different exact short-time onsets. Geometry-only resource or energy-moment
comparators cannot reproduce all actual-ray phase and circulation responses.
Equations (1)–(4) add a full-law relation between measured removable-edge weight,
fixed-host residual, and matter-response error. These distinguish mechanisms
inside prospective quantum simulations; they do not yet exclude conventional
quantum physics, which also implements this Hamiltonian.

The next decisive experiment is a finite-time **unfiltered** Q=2N evolution with
graph choice beyond the exact onset kernel. Measure charged network weight, persistent protected circulation,
intrinsic ball growth, spectral geometry and spatial density response in the
same state. Retain the pre-intervention distance/record and compare it with
final-distance classification. Use the actual coherent full-state residual or
explicitly costed physical records. Compare random regular, thin, fragmented
and supplied-dimensional controls at matched resources. A nonzero small-system
signal alone is insufficient; a controlled size-dependent joint window would
justify the continuum matching calculation.

## Sources

Source records, inspected primary URLs, scope and retrieval date are in
`SOURCES.json`. S1–S10 denote those records. Earlier repository packets are the
primary source for the frozen law and model-specific theorems; their preserved
files and immutable commits should accompany the integrated manuscript.
