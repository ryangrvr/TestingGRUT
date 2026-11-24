# GRUT-D3 record-dependent repair

24 November 2026. The physical question is whether memory acquired in one
part of an autonomous evolution can subsequently select which configuration
is restored after damage. D2 supplied a passive reservoir count. Here a
separate local three-state memory records one of two equally treated
configurations and controls which transition repairs a designated fault.

The active system consists of two three-level sites S and M and two finite
spin chains W and B. S has states 0, 1 and e (erased). M has states 0, 1 and
vac (blank). W receives the excitation emitted during writing. B receives
the excitation emitted during repair. The full experiment also includes
an inert three-level fault register D with S's basis and an inert fair
coin C. Its dimension is 54 times 2^(L_W+L_B). A classical record suffices for
the experiment, so the initial S state is a uniform incoherent mixture of
0 and 1. M starts blank, both chains start in their all-zero basis states, and
D starts in e. C is maximally mixed and appears in every condition.

For a chain J with occupations n_(J,j) and length L_J, use
`b_(J,j)=|0><1|`, `N_J=sum(n_(J,j))`, and adjacent XX exchange. Define the
conserved charges `Q_W=P_vac^M+N_W` and `Q_B=P_e^S+N_B`, present during
every uninterrupted interval. Put hbar = g = 1 and Delta = 4g.
The time-independent Hamiltonian is

$$
\begin{aligned}
H={}&\Delta(Q_W+Q_B)
+g\sum_{J=W,B}\sum_{j=1}^{L_J-1}
(b_{J,j}^\dagger b_{J,j+1}+\text{h.c.})\\
&+g\sum_{r=0}^1
\big(P_r^S\otimes |r\rangle\langle\mathrm{vac}|_M\,b_{W,1}^\dagger
  +P_r^M\otimes |r\rangle\langle e|_S\,b_{B,1}^\dagger
  +\text{h.c.}\big).
\end{aligned}
$$

These terms act on at most three sites. The two chain ends attach to S and
M through the corresponding three-site exchanges. The two designated
configurations are degenerate and exchanged by the simultaneous flip
`G=F_S F_M`. Each F swaps 0 and 1 and preserves the third state. G commutes
with H. The local writing and repair rules, fault alphabet, ready memory
and vacuum chains are supplied inputs. The selected label comes from the
earlier configuration.

## Experimental sequence

1. Evolve for a writing time w under H. For initial configuration r, S
   stays r because Q_B=0. Q_W=1, so the writer is exactly an open path of
   L_W+1 sites. Its origin is M blank and W vacuum. At x>=1, M=r and
   the single W excitation is at x.
2. At w, apply one externally timed SWAP of S and D in every condition.
   This erases S to e and stores r in the inert D. This operation is a
   reversible disturbance with a prepared fault register. The dynamics
   acts as the identity on D and C outside the specified interventions.
3. Leave M as written, flip its labels with F_M, or apply F_M when C=1.
   The coin mask removes accessible information about r from the active
   system while correlation with C retains it in the total state. The
   unconditional active density matrix and the entire energy distribution
   match across these comparisons. After the fault F_M acts on the occupied
   support exactly as G, even for the coherent blank component of the record.
4. Continue under the same H for a duration u. Score a single terminal
   read of S and M relative to D with all outcomes included. H has no
   coupling to the historical label in D or the coin in C. The supplied
   Schrodinger time and pulse timing remain external to the model.

After the disturbance both charges equal one. For a fixed selected memory
label m, identify x with the writer excitation location and y with the
repair excitation location. At x=0 M is blank, and at x>0 M=m. At y=0 S=e,
and at y>0 S=m. The invariant sector is the rectangular product of paths
of sizes L_W+1 and L_B+1 with the two edges touching (0,0) deleted. When
both sites are blank or erased there is no label to transfer. This isolates
the common corner. Other edges follow the two always-present hopping rules.
The identical sector matrices for m=0 and m=1 give the matched controls.

The primary outcome is the unconditional probability of S=D and M=D at
the final read. The erasure symbol and blank never count as historical
success. A second score counts S=D regardless of M. The probability that
S is in either of its two ordinary configurations is separately measured.
It can stay unchanged across controls even when the historical match fails.
