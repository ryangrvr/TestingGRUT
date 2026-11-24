# GRUT-D3 record-dependent repair: exact finite proof

24 November 2026. Protocol `GRUT-D3-01` concerns **classical historical-bit
restoration with a retained pointer after one declared fault**. The local
writing and repair rules, two degenerate configurations, ready memory,
vacuum chains, fault alphabet, timing and readout are supplied assumptions.
The proof establishes the registered finite-window mathematical gates under
standard unitary quantum mechanics.

## 1. Charges and acquisition of the record

Every term of the Hamiltonian in `MODEL.md` commutes with

$$
Q_W=P_{\mathrm{vac}}^M+N_W,
\qquad Q_B=P_e^S+N_B.
$$

An interface lowers its site's blank or erased occupation while creating one
excitation in the corresponding chain, or performs the reverse transition.
Chain hopping preserves excitation number. The other interface preserves
the same charge through its configuration projector.

Initially `Q_W=1` and `Q_B=0`. Positivity of both summands of `Q_B` keeps S
in an ordinary configuration and B in its vacuum. The writer preserves S's
initial label `r`. For this fixed label its exact basis is

$$
\begin{aligned}
|0\rangle_W&=|r\rangle_S|\mathrm{vac}\rangle_M|\mathrm{vac}\rangle_W,\\
|x\rangle_W&=|r\rangle_S|r\rangle_M|1_x\rangle_W,
\qquad 1\le x\le L_W.
\end{aligned}
$$

Apart from the constant energy `Delta`, its Hamiltonian is the adjacency
matrix of an open path with edge weight g. Let its hopping propagator be
`U_W(t)`, and write

$$
w_x(t)=\langle x|U_W(t)|0\rangle,
\qquad a_W(t)=w_0(t),\qquad c_W(t)=w_1(t).
$$

The memory equals the historical label with probability

$$
P(M=r\text{ at }w)=1-|a_W(w)|^2.
$$

This record is acquired under the same fixed Hamiltonian later used for
repair. Its two possible values are determined by the initial label.

The common fault swaps S with the prepared register D. Afterwards S is e,
D contains r, and `Q_W=Q_B=1`. D remains inert during the ensuing evolution.
The fault changes the charge sector through a declared external operation.

## 2. Exact repair sector and its dark corner

Fix the selected memory label `q`. In the intact arm `q=r`. In the wrong arm
`q=1-r`. Label basis states by `(x,y)`, where

$$
\begin{array}{c|cc}
 &0&\text{positive coordinate}\\ \hline
x&M=\mathrm{vac},\ W\text{ vacuum}&M=q,\ W\text{ excitation at }x\\
y&S=e,\ B\text{ vacuum}&S=q,\ B\text{ excitation at }y.
\end{array}
$$

For `x,y>0`, both sites contain q. At `(0,0)` both sites lack a label, so
neither interface can act. At `y=0`, the writer interface is blocked by S=e.
At `x=0`, the repair interface is blocked by M blank. All other edges follow
the corresponding uniform path exchange.

Thus the exact hopping matrix is

$$
T=T_0-V,
\qquad T_0=T_W\otimes I+I\otimes T_B,
$$

$$
V=g\bigl(|00\rangle\langle10|+|00\rangle\langle01|
          +|10\rangle\langle00|+|01\rangle\langle00|\bigr).
$$

The two edges at the origin are deleted and the corner is isolated. Each
label gives an isomorphic invariant subspace. They share this one dark
corner and are orthogonal on the remaining basis states. They are therefore
described as overlapping invariant subspaces, rather than an orthogonal
direct sum including two copies of the corner.

The on-site term is `2 Delta` throughout this sector and contributes only a
common phase. The post-fault initial vector is exactly

$$
|\psi(0)\rangle
=\sum_{x=0}^{L_W}w_x(w)|x,0\rangle.
$$

Joint restoration and retention is the event `x>0,y>0` in the intact arm.
It means S=D and M=D at the terminal read, with all other outcomes counted
as failures.

## 3. Uniform finite-path bounds

For a half-line with source at its endpoint, the hopping propagator has the
spectral integral

$$
U^{\infty}_{ij}(t)=\frac2\pi\int_0^\pi
\sin(ik)\sin(jk)e^{-2igt\cos k}\,dk,
\qquad i,j\ge1.
$$

The finite trigonometric identity

$$
\frac{\sin(ik)\sin(jk)}{\sin k}
=\sum_{r=0}^{\min(i,j)-1}\sin[(i+j-1-2r)k]
$$

has a derivative with absolute integral at most `2ij`. Its endpoint values
vanish. Integration by parts therefore gives

$$
|U^{\infty}_{ij}(t)|\le\frac{2ij}{\pi gt}.
$$

In particular, the first two amplitudes from the source obey

$$
|a^{\infty}(t)|\le\frac2{\pi gt},
\qquad |c^{\infty}(t)|\le\frac4{\pi gt}.
$$

For a finite chain of L reservoir sites plus its source, a walk from the
source to either of the first two sites must make at least `2L` hops to
cross the missing far-boundary edge and return. The finite and half-line
matrix powers agree below that conservative threshold. Both hopping
matrices have norm at most `2g`, which gives the uniform amplitude error

$$
\epsilon_L(T)\le2\sum_{k=2L}^{\infty}\frac{(2gT)^k}{k!}
\le2\left(\frac{e(2gT)}{2L}\right)^{2L}.
$$

The last estimate applies when `2gT<2L`. The registered resources use

$$
\begin{aligned}
L_W&=2688, &gw&\in[640,656], &g(w+u)&\le672,\\
L_B&=64,   &gu&\in[8,16].
\end{aligned}
$$

Both comparisons have `(2gT)/(2L)=1/4`. Using `e<3`, set the exact rational
bounds

$$
\epsilon_W=2(3/4)^{5376},
\qquad \epsilon_B=2(3/4)^{128}.
$$

The equality of this ratio to `1/4` is allowed because the strict inequality
needed here is `e<3`. Put

$$
\beta=\frac2{3.14}=\frac{100}{157},
\qquad A=\frac\beta{640}+\epsilon_W.
$$

The mathematical bound `pi>3.14` gives, throughout all required writer times,

$$
|a_W|\le A,\qquad |c_W|\le2A.
$$

For dimensionless elapsed repair time `v=gs>=1`,

$$
|a_B(s)|\le\frac\beta v+\epsilon_B,
\qquad |c_B(s)|\le\frac{2\beta}v+\epsilon_B.
$$

At shorter times the exact finite-path unitarity bound
`|a_B|^2+|c_B|^2<=1` will suffice. In particular, record acquisition has the
uniform certificate

$$
P(M=r\text{ at }w)\ge1-A^2>0.99999.
$$

## 4. A controlled comparison with product evolution

The product-path comparison vector is

$$
|\phi(s)\rangle
=U_W(w+s)|0\rangle\otimes U_B(s)|0\rangle.
$$

It has the same initial vector as the exact cut-corner evolution. Duhamel's
formula bounds their state-vector distance by

$$
\|\psi(u)-\phi(u)\|
\le\int_0^u\|V\phi(s)\|\,ds.
$$

The perturbation acts on three coordinates. Its exact norm is

$$
\frac{\|V\phi\|^2}{g^2}
=|c_Wa_B+a_Wc_B|^2+2|a_Wa_B|^2.
$$

This expression retains the interference cross term. The following bounds
use absolute values and apply for either sign of that term.

For `0<=v<=1`, the writer estimates and finite-path unitarity give

$$
\frac{\|V\phi\|^2}{g^2A^2}
\le6|a_B|^2+4|a_Bc_B|+|c_B|^2
\le7.
$$

The last inequality bounds the largest eigenvalue of
`[[6,2],[2,1]]` by 7. Since `sqrt(7)<8/3`, the contribution of this interval
to Duhamel's integral is at most `8A/3`.

For `1<=v<=16`, the three amplitude bounds above imply

$$
\frac{\|V\phi\|}{g}
\le A\sqrt{(4\beta/v+3\epsilon_B)^2
                    +2(\beta/v+\epsilon_B)^2}
\le A\left(\frac{9\beta}{2v}+\frac72\epsilon_B\right).
$$

Here the triangle inequality for the two real vectors uses
`sqrt(18)<9/2` and `sqrt(11)<7/2`. Integrating with `ln(16)<3` gives the
single rational error bound

$$
E:=A\left(\frac83+\frac{27}2\beta+\frac{105}2\epsilon_B\right),
\qquad \|\psi(u)-\phi(u)\|\le E.
$$

For the product vector, failure of the joint event occurs at `x=0` or `y=0`.
Its failure amplitude is at most `|a_W|+|a_B|`. For every registered `w,u`,

$$
\|P_{\mathrm{failure}}\psi(u)\|
\le\eta:=A+\frac\beta8+\epsilon_B+E.
$$

All constants in this expression are rational. Their exact comparison gives

$$
\boxed{
P_{\mathrm{intact}}(S=D,M=D)
\ge1-\eta^2
>0.9915>0.99.
}
$$

The displayed approximation is `1-eta^2 = 0.9915682364924134`. The gate
uses the rational inequality rather than rounded propagation samples.
It covers the complete rectangle `gw in [640,656]`, `gu in [8,16]`.
Either single-site historical success probability is at least this joint
probability. The reference retains its initial label throughout repair.

## 5. Causal record controls with matched energy distributions

Let `F_S` and `F_M` exchange configurations 0 and 1 while preserving the
third state. Simultaneous exchange

$$
G=F_SF_M
$$

commutes with H, as is seen by relabeling the two interface sums. At the
post-fault checkpoint, S=e, so applying `F_M` to the state acts exactly as
G. This includes the coherent blank-corner amplitude. If `rho_r` denotes
the checkpoint state conditioned on the historical label r, then

$$
F_M\rho_rF_M=G\rho_rG.
$$

Every spectral projector of H commutes with G. Thus every energy measurement
probability is identical across the memory arms, conditional on the actual
historical label and on the coin value. The comparison matches the whole
energy distribution on these states.

Fair averaging over r makes the unconditional active checkpoint state
invariant under `F_M`. The fair coin is present in every arm. Consequently
the complete unconditional active-plus-coin state is identical in the
intact, wrong and scrambled arms. The arms differ in their correlations
with the historical label held in D. The current system S is e in each arm.
`CONTROLS.md` specifies the operational resource comparison.

In the wrong arm, the invariant selected-label sector contains only S=e or
S=1-r. Its restoration and joint historical scores are exactly zero at all
post-fault times. The coin-scrambled arm is the equal mixture of this arm
and the intact arm. Therefore, exactly for every writing and repair time,

$$
\begin{aligned}
P_{\mathrm{wrong}}(S=D)&=P_{\mathrm{wrong}}(S=D,M=D)=0,\\
P_{\mathrm{scrambled}}(S=D)&=\tfrac12P_{\mathrm{intact}}(S=D),\\
P_{\mathrm{scrambled}}(S=D,M=D)
&=\tfrac12P_{\mathrm{intact}}(S=D,M=D).
\end{aligned}
$$

It follows throughout the registered window that

$$
\boxed{
P_{\mathrm{intact}}(S=D,M=D)
-P_{\mathrm{scrambled}}(S=D,M=D)>0.495.
}
$$

The probability of reaching either ordinary S configuration is invariant
under G and is exactly the same in all three arms. The intervention changes
which historical configuration is restored, rather than the overall rate
of leaving the erasure state. The inverted pointer preserves information
while reversing its association with history. The coin mask removes the
accessible association from the active state and retains correlations in
the full state including D and C.

During writing the mean energy is Delta. Hopping has zero expectation for a
path endpoint evolved by its own hopping matrix. Its endpoint-to-neighbor
coherence is purely imaginary on this bipartite path, so blocking that
interface at the common fault leaves zero mean hopping energy. The fault
changes `Q_B` from zero to one and costs mean work Delta. Afterwards the
mean energy is `2 Delta` in every arm. The memory interventions preserve
the energy distribution on the declared checkpoint states. Their coherent
control and timing remain supplied resources.

## 6. Exact limitations witnessed by smaller or different faults

For the registered short-chain control `L_W=L_B=1`, a writing time
`gw=pi/2` puts the record at `x=1` with certainty. Following the fault, the
nontrivial component of the graph is the three-site path

$$
|10\rangle\longleftrightarrow|11\rangle
\longleftrightarrow|01\rangle.
$$

Starting at `|10>`, its amplitudes after time u are

$$
\begin{aligned}
|\psi(u)\rangle={}&
\frac{1+\cos(\sqrt2gu)}2|10\rangle
-\frac{i\sin(\sqrt2gu)}{\sqrt2}|11\rangle\\
&+\frac{\cos(\sqrt2gu)-1}2|01\rangle.
\end{aligned}
$$

At `gu=pi/sqrt(2)` the state is `-|01>`. S is restored with probability one
and M is blank with probability one. The joint repair-and-retention score
is zero. Finite reversible transport can restore the system while consuming
the pointer, which explains the separate retention gate and the supplied
long writer chain.

An unflagged configuration flip gives a different exact failure. Replace
the erasure fault by `F_S`, leaving S in an ordinary but wrong configuration
and B in its vacuum. Then `Q_B=0` remains conserved. Both `P_e^S` and `N_B`
stay zero, so the repair interface stays inactive. The writer terms preserve
S's ordinary label. S therefore retains the flipped label and never returns
to its historical value under this H. The accepted fault class is the
specified flagged erasure.

## 7. Mathematical verdict and supplied physics

The registered gates establish **classical historical-bit restoration with
a retained pointer after one declared fault**. A local record acquired
during writing selects the restored configuration under the same fixed H.
The intact, wrong and coin-scrambled comparisons match the current active
state, dimensions and conditional energy law while changing historical
correlations. The retention and advantage certificates hold uniformly over
the finite writing and reading windows.

The Hamiltonian explicitly supplies the conditional copying and repair
relation. The ready memory and finite reservoirs supply its initialization
resources. The common fault and its pulse timing are external, and the
Schrodinger parameter remains supplied without an internal clock. Repeated
repair lifetime, an emergent controller, a novel physical law and a
cosmological selection rule remain open questions. The result belongs to
standard quantum mechanics and leaves Candidate 1 unchanged.
