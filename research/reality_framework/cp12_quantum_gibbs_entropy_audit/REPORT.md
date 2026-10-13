# CP-12: quantum Gibbs entropy bound for specified spatial graph families

**2026-10-13 — INDEPENDENT CANDIDATE THEOREM; PARTNER AUDIT PENDING.**
Source law: unchanged `research/constructive_physics/attempt12/LAW.md` at `11fbaec`. Parent of this independent additive branch is GRUT III adjudication `694463f`. No change to frozen CP-12, GRUT-RAI, or Stage-3 Card 3.

## Finite-cutoff full-Hamiltonian result

Let `N` vertices, `M=N(N-1)/2` potential undirected connections and conserved `Q` occupied connections. Fix an **electric cutoff** `|E|<=L` for the compact U(1) rotor legs, with `L>=1`. Let `H` be the *complete* CP-12 Hamiltonian in this finite-dimensional physical Hilbert space, retaining fermion hopping, Gauss constraints, electric energy, mass and the positive graph-rewiring Laplacian. Let `S` be **any** subset of the labelled `Q`-edge graphs and `P_S` its projector, including all allowed matter/electric states on those graphs. For fixed inverse temperature `beta>=0`,

```text
Pr_beta(S) = Tr(P_S exp(-beta H))/Tr(exp(-beta H))
  <= min(1, (|S|/binom(M,Q))*binom(N+Q,Q)*(2L+1)^(2Q)
          * exp[beta*(2|t|Q + Gamma Q(M-Q)/M)] ).
```

**Proof.**

1. On each labelled graph `G`, the neutral-sea state `Omega_G` (no vertex fermions, all edge-center fermions occupied, zero electric fields) obeys every Gauss constraint. Its exact diagonal energy is `kappa Q+Gamma Q(M-Q)/M`. There are `binom(M,Q)` orthonormal such states. Spectral Jensen inequality gives

```text
Z=Tr(exp(-beta H)) >= binom(M,Q)*exp[-beta*(kappa Q+Gamma Q(M-Q)/M)].
```

2. Positive vertex-mass/electric terms and the positive rewiring Laplacian cannot lower the energy. There are `2Q` incidence hopping terms of individual operator norm at most `|t|`; therefore `H >= (kappa-2|t|)Q I`, and `||exp(-beta H)|| <= exp[beta*(2|t|-kappa)*Q]`.

3. On each graph, `Q` active fermions occupy `N+Q` available active fermion modes; Gauss law only removes options. At most `(2L+1)^(2Q)` electric configurations exist. Therefore `rank(P_S) <= |S| binom(N+Q,Q)(2L+1)^(2Q)`. The positive trace bound `Tr(P_S exp(-beta H)) <= rank(P_S)||exp(-beta H)||` proves the result by division.

This argument includes **all noncommuting graph–fermion–gauge mixing** inside the finite-cutoff Gibbs state. It does not require assigning each graph an independent classical energy.

## Thermodynamic corollary (narrowly defined graph targets)

If `Q=nu*N` with `nu>1` fixed and `|S_N| <= N! exp(a*N)` for fixed `a`, at fixed finite `beta,L,t,Gamma,nu` we have `log binom(M,Q)=nu*N*log N+O(N)`, `log binom(N+Q,Q)=O(N)` and `log N!=N log N+O(N)`. Hence

```text
log Pr_beta(S_N) <= -(nu-1) N log N + O(N).
```

For a specified periodic square-lattice isomorphism class, `nu=2` and the exponent begins `-N log N`. For a specified periodic cubic-lattice class, `nu=3` and it begins `-2N log N`. Exponentially many unlabelled small local variants are allowed under the stated count bound. **This does not establish a no-go for all possible operationally 3D graph families** whose cardinalities are not controlled.

## Independent controls

The full local reproducibility packet `GRUT_CP12_Quantum_Gibbs_Entropy_Audit_2026-10-13.zip` contains `REPORT.md`, `check.py`, `small_full_law.py`, and `results.json`. Both scripts were executed without failures.

- `CP12_QUANTUM_GIBBS_ENTROPY_ANALYTICAL_CHECK_PASS` checks fermionic/rotor dimension bounds and exact logarithmic finite-N expressions through N=10^6.
- `CP12_GIBBS_FULL_LAW_SMALL_CHECK_PASS` independently builds the full Gauss-constrained Hamiltonian and tests neutral-sea diagonal, Hermiticity, lower spectral bound and projector Gibbs inequality at (N,Q,L)=(3,2,1),(4,2,1),(4,3,1), finding Hilbert dimensions 30, 147, 688. At beta=1 the single-graph probabilities are 0.333333, 0.066761, 0.046674, respectively. The general upper bound is *vacuous* at those small sizes; these checks are structural, not a numerical proof of the asymptotic.

## Limits / requirements for review

- The theorem is at **fixed finite electric regulator L**. The physical unregulated rotor Gibbs trace requires separate proof; do not transfer this verdict automatically.
- Fixed positive temperature is essential. It does not apply when beta grows with N or at zero temperature.
- The spatial family cardinality is a real assumption, not a derived bound for every 3D-looking graph.
- The labelled graph count is appropriate to the published physical basis. If permutations are gauge, the physical Hilbert quotient and counting must be re-derived; Hamiltonian permutation symmetry alone does not settle this.
- The canonical Gibbs state is not a cosmological formation protocol; metastable preparations and coherent nonequilibrium dynamics are untouched.
- Operator estimates are crude and small-N bounds vacuous. This is not an explanation of gravity, Lorentz symmetry, consciousness, or outcome selection.
- The entropy/geometry competition has significant prior art, especially Konopka (2008), Phys. Rev. D 78 044032 and Wilkinson/Greentree (2015), Phys. Rev. D 92 084007. Novelty has not been established.

**Adjudication request:** Claude Code and Work should independently verify the mass/electric positivity, hopping operator norm, neutral-sea partition lower bound and graph-fiber count. Then try to generalize to unregulated rotors or exhibit an operationally spatial family whose cardinality violates the low-complexity premise. Keep all full-law calculations in the same Q sector and state; do not retune H to rescue a desired phase.

## Extended scaling window: growing cutoff and falling temperature

The same exact finite-N inequality permits a cutoff and inverse temperature that vary with N, although it does not establish any **uniform approximation** to the unregulated rotor theory. Write `Q=nu N`, `L_N=N^(alpha+o(1))`, and `beta_N=b log N+o(log N)` with fixed nonnegative `alpha,b`; take `|S_N|<=N! exp(aN)`. Then straightforward Stirling analysis gives

```text
limsup_{N->infty} log Pr(S_N)/(N log N)
  <= -(nu-1)+2 nu alpha + b nu (2|t|+Gamma).
```

Hence the specified spatial family is still asymptotically excluded whenever

```text
2 nu alpha + b nu (2|t|+Gamma) < nu-1.
```

At fixed temperature (`b=0`), a slowly increasing cutoff `L_N=N^alpha` cannot rescue the family for `alpha < (nu-1)/(2nu)`: thresholds are `alpha<1/4` for `nu=2` and `alpha<1/3` for `nu=3`. With fixed cutoff (`alpha=0`), `beta_N=b log N` still cannot rescue it for `b < (nu-1)/(nu(2|t|+Gamma))`. For frozen `|t|=1,Gamma=0.1` these bounds read `b<0.238095...` (nu=2) and `b<0.317460...` (nu=3).

These are **sufficient suppression conditions** derived from a loose upper bound, not phase transition points or guarantees of selection on the other side. Taking `L_N->infty` in a slow joint scaling does not establish the physical `L=infty` Gibbs trace at each N. The regulator-removal order remains an independent obligation.
