# External RAI review packet

**Status.** Ready for independent review. The review and second implementation included here are internal to the GRUT project. The result has no external replication at publication.

## Question and scope

Test whether the absence of exact connected ordinary internal global phases on light local fields suffices to select the electric charge assignment. The criterion is motivated by a restricted consequence of the AdS/CFT no-global-symmetry result. The test class is a weakly coupled four-dimensional spin Wilsonian gauge EFT. The gauge group is $(SU(3)\times SU(2)\times U(1)_y)/\mathbb Z_6$ with its standard diagonal generator. Keep one positive Higgs weight $y_H=3$, three family-universal standard nonabelian matter representations, one $n^c$ per family, generic four Dirac Yukawa matrices and arbitrary finite-dimension local polynomials. Respect the cutoff and declare all special preparations or matching assumptions. [Frozen protocol](FROZEN_PROTOCOL.md) · [Parameters](PARAMETERS.json)

## Result to try to refute

The family-universal Yukawa and anomaly conditions leave the weights
\[
y=(1,-4,2,-3,6,0,3)+6k(1,-1,-1,-3,3,3,0),\qquad k\in\mathbb Z.
\]
The global form fixes the allowed integer lattice. Modulo gauge, the renormalizable phase freedom has directions $B$ and $L$, of which $B-L$ is anomaly-free. The $QQQL$ invariant has $(\Delta B,\Delta L)=(1,1)$. The operator $n^c n^c$ multiplied by $3|k|$ powers of $C_\mu C^\mu$ or its conjugate has $(0,-2)$ with $C_\mu=H^T i\sigma_2D_\mu H$. The charge of $C_\mu$ is $6$. Choose the conjugate for $k>0$ and the original for $k<0$. These terms remove both connected phase directions modulo gauge. For $k=0$ omit the dressing.

Color triality and Lorentz parity give $|\Delta[3(B-L)]|\ge6$ when nonzero. The two single-field inequalities $|y_0-s\,3(B-L)|\le3 d_{\rm eng}$ for $s=1/2$ and $3/2$, together with gauge neutrality, are claimed to give the lower bounds $12|k|-1$ for $k<0$ and $12k+3$ otherwise. The upper witness has dimension $3+18|k|$. The physical observable $q_n/|q_e|=3k/|1+3k|$ is independent of charge-normalization convention, assuming an identifiable long-lived $udd$ state and ordinary monopole measurement. The minimum noncanonical magnitude is $3/4$. [Derivation](DERIVATIONS.md) · [Empirical comparator](REALITY_TEST.md)

## Required challenges

1. Recompute the gauge kernel on the covering group and trace perturbative anomalies for arbitrary integer $k$. Check whether the use of the spin-bordism theorem fully addresses the stated global gauge anomaly scope. Identify extra assumptions for nonspin and quantum-gravity contexts.
2. Classify the connected flavor centralizers with the frozen genericity conditions. Test whether a residual ordinary global phase survives modulo the gauge subgroup, and whether either explicit breaking vertex vanishes identically or is redundant. Reconstruct the charged vector contact or another on-shell test.
3. Challenge the lower bound using derivatives, conjugates, flavor contractions, integration by parts, equation-of-motion relations and gauge-compatible field redefinitions. Distinguish the engineering dimension above electroweak breaking from an observable rate after matching.
4. Check local causal perturbation theory, radiative stability and the absence of a uniform ultraviolet claim. Ask which independent premise would remove the noncanonical solutions. If a stronger conjecture is invoked, derive its applicability and physical input requirements separately.
5. Verify that the physical charge difference cannot be removed by normalization, renormalization scheme or source calibration. Audit the actual neutron species, the assumptions of the monopole comparison, both reported uncertainty evaluations and the measurement's response range.
6. Compare priority against established charge dequantization and high-dimension accidental symmetry examples. In particular compare any new alleged implication to Fichet and Saraswat's quantitative rate conjecture and the conventional Majorana condition. Do not count a different notation or a matched formula as new physics. [Novelty audit](NOVELTY_AUDIT.md)
7. Run the registered commands. Report algebraic objections, differing results, and the resulting A/B/C/D classification with scope and supporting sources.

## Inputs and expected verification

- Parent is [PR #26](https://github.com/ryangrvr/TestingGRUT/pull/26) at 95ad5c0f158b128d337da8fe862894015b5d8f8c.
- [Freeze commit](https://github.com/ryangrvr/TestingGRUT/commit/2f6a734a456c23272b8a413afa880eb4234f83f5) precedes production. Run only the registered five fixtures when assessing the coded result.
- [Reproduction](REPRODUCTION.md) includes two exact-arithmetic implementations, the retained first-run check failure and archived outputs.
- Primary references include [Babu and Mohapatra](https://doi.org/10.1103/PhysRevD.41.271), [Davighi et al.](https://arxiv.org/abs/1910.11277), [Harlow and Ooguri](https://doi.org/10.1103/PhysRevLett.122.191601), [Fichet and Saraswat](https://doi.org/10.1007/JHEP01(2020)088), and [Baumann et al.](https://doi.org/10.1103/PhysRevD.37.3107).
