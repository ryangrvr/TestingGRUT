# CP12: fixed connection number cures density, but not every entropy obstruction

This is a theorem for the frozen labelled-graph physical Hilbert space and canonical Gibbs state, not for all measures, preparations or emergent-gravity programs. No new candidate or thermal scaling is introduced.

## Exact sparsity

[H,Q]=0, so in the Q=νN sector the edge density is identically Q/M~2ν/N at EVERY temperature. CP11's fixed-temperature dense-edge failure is genuinely removed within this supplied conserved-resource premise. No calculation determines ν or proves nature has this inventory.

## Physical partition functions are finite

Let D=κQ+2m N_vertex+(ε_E/2)Σ E_l² on the exact Gauss-law subspace with inactive legs frozen. At finite N,Q and β>0, the trace of exp(-βD) is finite. For any graph G let z_G be its diagonal partition function. Its filled-edge, E=0 state gives z_G≥exp(-βκQ). Discarding Gauss constraints gives the valid overcount

 z_G≤exp(-βκQ) 2^(N+Q) θ(βε_E/2)^(2Q),

where θ(a)=Σ_{z∈Z} exp(-a z²)<infinity. The 2^(N+Q) factor overcounts active fermion occupations; it includes no independent states on absent connections. These bounds also hold uniformly for electric cutoffs, using the full θ as an upper bound.

There are K_N=binomial(M,Q) labelled graph configurations. For a graph class A, in the reference Gibbs state σ=exp(-βD)/Z_D,

 σ(A) ≤ (|A|/K_N) exp[C_β(N,Q)],
 C_β=(N+Q)ln2+2Q ln θ(βε_E/2).

## Quantum mixing cannot be treated as a classical weight — use relative entropy

Write H=D+V. There are 2Q active fermion-hopping terms, each norm≤t, so ||H_hop||≤2tQ. The reversible-rewiring Laplacian has maximum degree≤Q(M-Q), giving ||H_rewire||≤2ΓQ. Thus

 ||V||≤2(t+Γ)Q.

For ρ=exp(-βH)/Z_H, the exact quantum relative entropy is

 S(ρ||σ)=-β Tr(ρV)+ln Z_D-ln Z_H ≤2β||V||≤4β(t+Γ)Q.

The partition inequality used here follows from eigenvalue bounds for a bounded perturbation; it remains valid although H and D do not commute. Measurement of the projector onto A and relative-entropy monotonicity yield

 S(ρ||σ) ≥ p ln(1/q)-ln2,
 p=Tr(P_Aρ), q=Tr(P_Aσ).

Whenever the denominator is positive,

 p ≤ [4β(t+Γ)Q+ln2] /
     [ln(K_N/|A|)-(N+Q)ln2-2Q ln θ(βε_E/2)].

This proof handles the full dynamical gauge–fermion–connection Gibbs state. It does not replace the quantum hopping by a static optimized flux.

## Restricted spatial-family implication

Suppose Q=νN+O(1), ν>1 fixed, and a target family has at most N! exp(sN+o(N)) labelled graphs, for fixed s. This covers relabellings of a fixed lattice sequence and any explicitly bounded exponential number of variants. It does NOT automatically cover every curved spacetime discretization or all possible spatial graphs.

 ln K_N = νN ln N + O(N),
 ln |A| ≤ N ln N + O(N).

For fixed β, the denominator is (ν-1)N ln N-O_β(N), while the numerator is O_β(N). Therefore

 P_β(A) = O_β(1/ln N) -> 0.

Resource conservation removes the O(N²) dense-edge entropy, but an O(N ln N) counting disadvantage still defeats these target families under finite-strength extensive interactions at fixed temperature. The bound is conservative; it is not a practical finite-size crossover estimate.

For ν=1 this argument gives no asymptotic exclusion: both leading counting powers coincide. For β depending on N, the fixed-temperature conclusion cannot simply be retained. We adopt neither temperature scaling nor a changed measure as a repair. At zero temperature the relative-entropy bound above is not decisive.

## Priority and interpretation

Superextensive graph entropy and vanishing transition temperatures are already discussed by Si Chen and Steven S. Plotkin, *Statistical mechanics of graph models and their implications for emergent spacetime manifolds*, arXiv:1210.3372 / Phys. Rev. D 87, 084011 (2013). The model-specific quantum relative-entropy bound here is an author derivation, pending review, not a claim that the underlying entropy problem is new.

Quantum coherence, gauge constraints and conserved resources all have physical content in CP12. They do not automatically make a canonical finite-temperature state concentrate on spatial geometry. A nonequilibrium state could behave differently, but would require its own actual preparation and evolution; it cannot be assumed to evade this result.
