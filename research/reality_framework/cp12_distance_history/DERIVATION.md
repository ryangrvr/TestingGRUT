# CP-12: distance histories, measurement backaction, and shortcut capacity

This packet uses the frozen CP-12 law and full labelled physical Hilbert space. Numerical projections bound total integer electric square; they are not changes to the law. Sources, detectors and records commute with Gauss generators. The input pure neutral ray and the chosen Q sector remain supplied.

## 1. Three distinct experiments

Choose source vertex 0 and receiver vertex 1 without supplying any graph or embedding. Let S=n₀, T=n₁, U(t)=exp(−iHt), and v=U(w)ψ₀, where ψ₀ is the neutral equal superposition over all Q-edge graphs. Let F_r project onto source–receiver graph distance r at the intervention time. The disconnected class is denoted infinity in prose and N in the data. These graph projectors satisfy Σ_r F_r=I and commute with S, T and the electric-budget projector Π.

An operational initial-distance record is the isometry

    M v = Σ_r F_r v ⊗ |r⟩.

Ignoring its orthogonal record states gives the nonselective Lüders state Σ_r F_r|v⟩⟨v|F_r. It is generally different from |v⟩⟨v|. The record is not a free, noninvasive classical label. It is an ideal diagnostic involving a global graph property; an autonomous local apparatus performing it has not been derived.

Write v_r=F_r v and p_r=||v_r||². Since S commutes with F_r, the source kick exp(−iηS) preserves p_r. The conditional linear susceptibility is therefore the branch-weighted susceptibility divided by p_r; no derivative of p_r is needed.

The joint signed response, classified both before and after intervention, is

    C_rs = 2 Im ⟨ U(τ)v_r | F_s T | U(τ)S v_r ⟩.

It is the derivative with respect to η, not a transition probability. Its row and column sums give

    χ_past(r) = Σ_s C_rs,
    χ_post,recorded(s) = Σ_r C_rs.

These are comparable classifications of THE SAME recorded experiment. Separately, the unrecorded state has

    χ_post,coherent(s) = 2 Im ⟨ U(τ)v | F_s T | U(τ)S v ⟩.

Comparing this last quantity directly to χ_past conflates distance reclassification and recording backaction. The exact difference between the unrecorded and recorded total responses is the sum of cross-record coherences:

    χ_total,coherent − χ_total,recorded
      = 2 Im Σ_(r≠r') ⟨ U(τ)v_r | T | U(τ)S v_r' ⟩.

The recorded shortcut therefore belongs to the INSTRUMENTED state. Recording can inject rewiring energy and disrupt cancellations in a coherent geometry. It must not be attributed unchanged to the unrecorded input. A candidate phase would need either a derived physical metric record or demonstrated small readout backaction for the observables used to define its causal region.

The ordered-pair quotient fixes vertices 0 AND 1 individually, retaining all permutations of other vertices with the parent fermionic signs. A source-only quotient cannot implement an individual receiver's distance record. Before recording, symmetry makes this receiver representative of the other vertices; after recording, receiver exchange symmetry is broken. Recorded receiver density must not be silently equated to the parent's average over all receivers.

## 2. Positive geometry histories and classification ambiguity

Define positive branch weights

    q_rs = ||F_s U(τ)v_r||²,
    b_rs = ||F_s U(τ)S v_r||².

They obey Σ_s q_rs=p_r and Σ_s b_rs=||S v_r||². The raw off-diagonal sum Σ_(r≠s)q_rs measures distance change in this recorded experiment, not identity or persistence of individual edges.

Cauchy–Schwarz gives |C_rs|≤2√(q_rs b_rs). Consequently

    χ_past(s)−χ_post,recorded(s)
      = Σ_(j≠s)(C_sj−C_js).

Its absolute value is bounded by the corresponding sum of 2√(q b) terms. For projected forward vectors with norm errors e_x,r and e_y,r, replace √q by √q_projected+e_x,r and √b by √b_projected+e_y,r. This yields the saved full-law upper bound from geometry-transition weights.

Joint response entries are signed and are not themselves conditional probabilities. A response at a final distance does not identify the distance or path used by matter during evolution.

## 3. Recording resources

For pure v, the orthogonal vectors v_r imply that the nonselective system state has nonzero eigenvalues p_r. Its entropy and min-entropy, in bits, are

    S_recorded = −Σ_r p_r log₂p_r,
    H∞,recorded = −log₂(max_r p_r).

A unitary recording on a pure apparatus need not increase TOTAL entropy: these are system/record entanglement entropies when the record is ignored. The apparatus, memory, clock and global distance computation remain external resources, not a derived selection mechanism.

The source-free neutral input has exactly ⟨H⟩=Q/2 and Var(H)=2Q. At the record time, define

    ΔE_record = Σ_r ⟨v_r,Hv_r⟩−⟨v,Hv⟩.

Only rewiring contributes, because graph projectors commute with mass, electric energy and hopping. Let R be the positive frozen rewiring operator. Then ΔE_record is the expectation of the bounded operator D_R=ΣF_r R F_r−R. A safe full-law bound is

    ||D_R|| ≤ ||R|| ≤ 2Γ Q(M−Q)/M.

The first inequality follows from 0≤R≤rI and 0≤ΣF_rRF_r≤rI; their difference lies between −rI and rI. Therefore a waiting-state norm error ε_w bounds the record-work error by

    |ΔE_record,true−ΔE_record,projected|
       ≤ [2Γ Q(M−Q)/M] (2ε_w+ε_w²).

Thus record work can be controlled despite unbounded electric energy; one must not justify it merely by a wavefunction norm error applied to an unbounded H. Conditional energies shown in the separate configuration-dephasing audit are finite-projection resource diagnostics and are not independently certified unbounded-energy expectations.

An abstract record requires a reversible graph-distance computation and at least enough distinguishable memory states for the available classes. No gate-count estimate for an actual CP-12 apparatus is established here. Such a global diagnostic operation does not constitute a physical preparation mechanism generated by H.

## 4. Full-law error control for recorded paths

Retain a Π-invariant signed source–receiver sector, with H_Π=ΠHΠ and B=(1−Π)HΠ. All full-law rewiring diagonals and outgoing destinations are preserved. Recording maps Π to the direct sum of Π copies. Source and receiver densities and all distance projectors commute with Π.

For recorded post-distance responses, use purification matrices X=[U_Π(τ)v_r] and Y=[U_Π(τ)S v_r]. Their inner product is Frobenius. The adjoints of the waiting-time insertions are

    M†[z_r] = Σ_r F_r z_r,
    M† S[z_r] = S Σ_r F_r z_r.

M is an isometry and S a contraction. Propagate the detector-weighted outputs backward through the post-source segment, apply these adjoints, and propagate backward through the pre-source segment. The preceding packet's residual orthogonality argument applies unchanged to this enlarged Hilbert space:

    |χ_true−χ_projected| ≤ 2(e_y d_y+e_x d_x+e_x e_y).

For each past-distance branch, the waiting-time insertions are F_r and S F_r. For unrecorded paths they are I and S. Recorded Frobenius residuals are integrated as matrix norms directly. Taking the square root of summed squared per-column integrals would not generally be an upper bound on the integral of the Frobenius norm; the code does not make that substitution.

## 5. Accurate leakage integration, not a physical Taylor extrapolation

With A=H_Π−Q/2 and mesh width h, expand the leakage VECTOR B exp(−iAs)x about each midpoint. Integrating its norm by the triangle inequality gives the valid upper sum

    Σ_midpoints Σ_(k=0)^5 [h^(k+1)/(2^k (k+1)!)] ||B A^k x_mid||
      + Σ_midpoints [h^7/(2^6 7!)] ||B|| ||A^6 x||.

Use ||B||≤sqrt(||B||₁||B||∞). Unitarity preserves ||A^6x||. The same formula holds for Frobenius norms and reversed evolution. The code also evaluates the parent's first-degree upper sum and takes the smaller of the two valid bounds.

The time-dependent states come from finite-time matrix exponentials. The polynomial expansion bounds local quadrature error in the residual integral; it is not an expansion of the measured response extrapolated to long times. Bounds are evaluated in ordinary floating point. Roundoff is tested independently, including an exact finite-pulse check, but no interval-arithmetic certificate is claimed. In particular, tiny analytic truncation bounds below machine-scale uncertainty must not be marketed as rounded certified response intervals.

Branch probability error is safely at most 2ε_w+ε_w². Conditional response intervals divide the unconditional interval by the resulting positive probability interval, including both signs correctly. The tables distinguish conditional suppression from mere rarity of far-distance pairs.

## 6. Exact mechanism controls

Hold the full-law waiting state and distance record fixed, and change terms only in diagnostic post-source comparator evolutions.

* Without rewiring, every graph sector is invariant. If source and receiver lie in different components in that graph, their even density algebras evolve independently and commute. Therefore χ_past(infinity)=0 exactly for EVERY input state supported on those graphs, including entangled inputs. Initial correlations alone cannot evade this operator statement.
* Without hopping, the remaining Hamiltonian commutes with S. The source-kick response of these commuting density/graph detectors vanishes identically, although rewiring can still change distance.

The global electric-square projection can artificially couple otherwise independent components. Accordingly, the projected no-rewiring response need not be treated as an exact physical null without its full-law error bound. The saved controls check compatibility with zero after bounding that error.

In the basis builder, rewiring preserves Π, so its outgoing electric-budget block is zero. The no-rewiring comparator has the same B as the full law. The no-hopping comparator has B=0. The code reconstructs the comparator rewiring operator from exact transitions rather than deleting entries according to representative graph labels: different physical rewirings can fold to the same canonical representative.

## 7. An exact edge-pool coupling norm

Let a specified pool contain L potential graph edges, and let Z project onto graphs containing NONE of them. Assume 1≤Q≤M−L, so that this subspace is nonempty. A single swap takes Z only into the sector with exactly one pool edge.

Mass, electric energy and hopping preserve graph bits. Thus the entire off-diagonal block (I−Z)HZ is a rewiring block. Each physical source configuration has at most QL legal cross-pool destinations. Each target with exactly one pool edge has at most M−L−Q+1 inverse destinations. All raw cross amplitudes have magnitude Γ/M. The row/column Schur bound gives

    ||(I−Z)HZ||² ≤ (Γ²/M²) Q L (M−L−Q+1).

This is equality, not merely an upper estimate. Restrict to neutral seas with all electric fields zero. Every edge is removable. The incidence matrix is biregular, with source degree QL and target degree M−L−Q+1. Uniform normalized vectors on the two graph strata attain its Schur bound. They are physical Gauss states, so the same lower bound applies to the full unbounded electric-field Hilbert space. Therefore

    g_L² = ||(I−Z)HZ||²
         = Γ² Q L (M−L−Q+1)/M².

No electric cutoff, weak-hopping limit, or time-series approximation enters this norm identity. Incidence enumeration checks its exact degrees and rational squared norm in nine finite cases, including a higher-density Q=2N example. Large-N rows in EDGE_SPEED.json are formula evaluations only.

For L=1,

    g_1² = Γ² Q(M−Q)/M².

If the selected edge is initially absent, differentiating its population p gives |p'|≤2g₁√(p(1−p)). Integrating from p(0)=0 yields the full-law bound

    P(edge present at τ) ≤ sin²(min(g₁|τ|,π/2)).

The same bound applies to losing an initially present edge. Mixtures follow by purification or convexity. The checker tests the corresponding weighted transition bounds in the recorded data. This is a bound on one specified edge, not on the entire change of graph distance.

At Q~νN, M~N²/2:

    g_1² ~ 2νΓ²/N,
    L~λN ⇒ g_L → Γ√(2νλ).

Thus normalization makes each specified channel weaker but does NOT make the collective operator capacity for O(N) candidate edges vanish. The cross SUBMATRIX that also keeps d fixed near edges has source degree (Q−d)L and the same target degree, leaving this limit unchanged for fixed d. This counting corollary still concerns a restricted capacity subspace, not the actual spatial phase.

The norm limit does not prove that a selected state realizes it. For a single complete physical configuration with no pool edges and r removable occupied edges, all Lr outgoing destinations are orthogonal, so its squared cross norm is only Γ² Lr/M². Coherent superpositions can enhance the cross norm. This difference is why operator capacity, actual-state exposure, and a finite-time causal cone must remain distinct questions.

## 8. Actual-state exposure and the next gate

For each recorded nonadjacent class, the audit computes ||F₁ H F_r v||²/p_r. Since hopping and diagonals preserve distance, this is a direct-edge rewiring matrix-element diagnostic. It is not a finite-time probability or a new time-series prediction.

Complete dephasing in the FULL physical configuration basis, including graph, fermions and rotors, gives the corresponding diagnostic Γ² Σ_s |(F_r v)_s|² r_removable(s)/(M² p_r). It accounts for orbit multiplicities in dephased entropy, and reports changed energies rather than presenting the comparison as free or energy matched. This is stronger dephasing than graph-only dephasing.

Its bounded-operator intervals use ||F₁HF_r||≤g₁ and the waiting-state error. For the dephased diagnostic, the diagonal norm is at most Γ²Q/M². These are evaluated numerically and retain the same floating-point qualification as the other bounds.

There is also a constructive state-dependent bound. Let r_rem count removable occupied edges, and let a normalized state have no pool edges. Each target has at most d_in=M−L−Q+1 parents. Apply Cauchy–Schwarz to the amplitude sum at each target and then sum targets. Each source has L r_rem outgoing destinations. This proves

    ||(I−Z)H Z ψ||²
       ≤ (Γ²/M²) L(M−L−Q+1) ⟨r_rem⟩_ψ.

For L=1 this implies that the actual-state/dephased-configuration coupling ratio is at most M−Q. The saved state audit checks that bound. For Q~νN and L~λN, the upper estimate becomes asymptotically 2Γ²λ ⟨r_rem⟩/N. Thus ⟨r_rem⟩=o(N) is a SUFFICIENT state-dependent condition to suppress collective pool coupling. It is not necessary: phases or orthogonal records can suppress interference even with many removable edges. Nor does pinning arbitrary edges select a spatial graph.

There is a finite-time consequence with an explicit dynamical condition. If the pool is initially absent and p(t) is the probability that it becomes nonempty, the flux identity gives

    |p'(t)| ≤ 2√p(t) [|Γ|√(L(M−L−Q+1))/M] √⟨r_rem⟩_t.

Integrating the bound on d√p/dt proves

    p(T) ≤ min(1, [Γ² L(M−L−Q+1)/M²]
                    [∫₀ᵀ √⟨r_rem⟩_t dt]²).

This is not a short-time Taylor approximation. If the removable-edge expectation is o(N) UNIFORMLY over a fixed time interval at Q~νN and L~λN, pool creation vanishes in that limit. The early numerical states here do not satisfy the required suppression, and suppression of geometry updates alone would still not prove a three-dimensional matter cone.

Gauss's edge equation gives a physical interpretation. An occupied edge with both electric half-fields zero must have its center fermion filled, so it is removable. Conversely, any nonremovable occupied edge has nonzero integer electric square. Consequently the diagonal operator inequality

    H_E ≥ (ε_E/2)(Q−r_rem)

holds exactly. With ε_E=1/4, suppressing removability requires at least 1/8 electric-energy unit per pinned edge. For a completely charge-neutral circulation state the stronger factor ε_E applies, because the half-fields are (f,−f). This energy must be accounted for together with graph cost, mass and hopping; it is not free preparation purity. A physical mechanism might create persistent gauge/matter records that inhibit rewiring, but that mechanism and its simultaneous spatial organization have NOT been demonstrated by these early-time runs.

The reported evolution series uses Q=N−1 as an AUXILIARY propagation test: connected graphs in that sector are trees and do not provide an extensive conventional cycle-space gauge capacity. It is not evidence for a three-dimensional photon phase. The pool-norm theorem also applies to the higher-density sectors, but their large-size source–detector dynamics is not calculated here.

A candidate spatial phase must demonstrate STATE-DEPENDENT suppression of collective long-distance updates and a controlled growing range of propagation while preserving its gauge capacity and organization. These results identify a concrete obstacle and a quantity to control; they do not establish a no-go against every state or derive the state-selection principle.

## 9. Exact finite-pulse extraction

For the projector source S²=S, exp(−iηS)=I+(exp(−iη)−1)S. Let f(η) be a target population after that pulse and the subsequent evolution, with the same initial-distance branch. Expanding only this TWO-EIGENVALUE source operation gives

    f(η)−f(−η)=2 sin(η) χ.

Thus χ=[f(π/2)−f(−π/2)]/2 EXACTLY, without a small-η finite-difference approximation. This identity holds for the recorded mixed experiment by linearity as well. The pulse checker computes the two populations independently with matrix exponentials and compares them with the bilinear susceptibility. It also reports the pulses' conditional work. The finite interventions have different energy costs; they are not being claimed as resource-matched preparations.
