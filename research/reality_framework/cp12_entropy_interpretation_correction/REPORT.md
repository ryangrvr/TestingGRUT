# CP-12 — Entropy Interpretation Corrections (2026-10-13)

**Independent additive audit; requires Claude/Work review.** Based on the unregulated CP-12 Gibbs entropy theorem at commit `4bb0fa5`; no change to frozen CP-12, GRUT-RAI, or Stage-3 Card 3. Detailed offline calculation and exact finite enumerations are in the attached research packet.

## 1. Uniform sparse graphs are typically symmetric, not asymmetric

For uniform labelled G(N,Q) with Q=nu N, M=binom(N,2), the isolation probability for one vertex is EXACTLY

```
p0 = binom(binom(N-1,2),Q) / binom(M,Q).
p00 = binom(binom(N-2,2),Q) / binom(M,Q).
```

For I=# isolated vertices, E I=N p0 and Var I=N p0(1-p0)+N(N-1)(p00-p0^2). At fixed nu, p0=e^(-2nu)+O(1/N), p00=e^(-4nu)+O(1/N), Var I=O(N). Therefore I/N converges in probability to e^(-2nu)>0. Every permutation among isolated vertices is a graph automorphism, so |Aut(G)|>=I!, with high probability. **Dividing both totals by N! on the ground that most sparse graphs are asymmetric is incorrect.** Random fixed-k regular graphs (k>=3) are different: they are typically asymmetric, connected, and nonspatial.

Permutations as gauge still require a full physical Hilbert-space quotient with fermionic signs/gauge states. This counting correction does NOT prove gauging selects spacetime.

## 2. Generalized full-quantum Gibbs entropy-exponent corollary

The accepted regulator-free CP-12 theorem at fixed beta>0 supplies, for ANY labelled Q-edge graph class S_N,

```
P_beta(S_N) <= O(1/log N)
 + exp(log |S_N| - log binom(M,Q) + O(N log log N)).
```

Define sigma=limsup log|S_N|/(N log N). At Q=nu N, nu>1, log binom(M,Q)=nu N log N+O(N). Therefore sigma<nu => P_beta(S_N)->0.

**Critical correction:** This does not require an admissible family to include almost every graph. A fraction exp(-cN) is still exponentially SMALL and is not excluded by the present estimate. All k-regular graphs, k=2nu, have log count=(k/2)N log N+O(N)=nu N log N+O(N) (McKay/Wormald), so their graph-count probability is exp[-Theta(N)] and their sigma=nu. The theorem alone CANNOT exclude this family. Being regular, however, does not imply physically three-dimensional. Any broader operational family needs its own rigorous count and dynamical observables.

## 3. Temperature, superselection, and locality

For sigma<nu, beta_N=o(log N) remains insufficient to overcome the entropy estimate; beta_N=Omega(log N) is a NECESSARY SCALE TO EVADE THE BOUND, not a critical temperature or evidence of order. The Q=N ring is one-dimensional and its ground-state fragmentation test is highly informative **for that candidate**, not decisive for the Q=2N or Q=3N sectors or physically explained nonequilibrium states. Full CP-12 conserves Q.

Uniform G(N,nu N) is disconnected with high probability (a positive density of isolates); it is **not an expander**. Random k>=3 regular graphs are typically expander-like, but those form a distinct constrained ensemble. Fixed-background Lieb-Robinson graph-distance bounds cannot be inserted into full CP-12's all-pairs graph-changing Hamiltonian without additional state/detector assumptions. Effective physical edge suppression remains an untested hypothesis.

## Next scientific obligation

1. Audit the **complete full-H** frozen-coupling unregulated Q=N ring-versus-fragments ground-state projector and spectral energy-density gap as N grows.
2. Independently evaluate Q=2N/3N candidate spatial-versus-random-regular and fragmented competitors, with operational dimension, physical responses, and preparation budgets in the SAME Q and state.
3. Analyze full permutation-gauge quotient separately; do not infer it from graph automorphism averages.
4. Do not call frozen CP-12 globally falsified by a failed spanning-ring candidate.

## Exact checker and prior art

The accompanying `check.py` independently evaluates exact p0,p00, variance and the fixed-k regular-graph asymptotic. It is NOT a rerun of the full unregulated quantum Hamiltonian.
McKay/Wormald, Combinatorica 11 (1991), https://doi.org/10.1007/BF01275671 . Sparse random graph isolated-vertex asymptotics are classical.

**Verdict:** Claimed Gibbs theorem preserved; two stronger inferences in the partner discussion corrected; no emergent GR, no full-CP12 no-go.
