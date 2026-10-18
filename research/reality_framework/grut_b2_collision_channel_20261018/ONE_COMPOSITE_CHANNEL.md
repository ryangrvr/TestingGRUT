# B1 input-to-output channel and spectral account

This clarifies the unchanged B1 result; it is not independent validation.

## Exact preparation-to-output instrument

Let Q be the input qubit and Z an arbitrary external reference, initially
rho_QZ. Let xi be the actual accepted word/storage/receiver state, independent
of QZ. Prepare the one root fuel excited with color Q, matter/fields vacuum.
The one-resource physical space is canonically H_spatial tensor C^2: the
unique color resides in fuel if F=1 and in the plus if F=0. All B1 gates,
energies, hopping and conversion are color-blind, including geometry moves.
Therefore U_B1(t)=U_spatial(t) tensor I_Q on this reachable space.

Let R be the ONE common success projection: retained recorded prefix intact,
co-located plus/minus at recorded distance >=4, and fuel empty. Sum over
all accepted final positions, all tail geometries and all actual receiver
fields. A decoder reads plus color while retaining position/environment
only as discarded ancillary information. Then

    E_success(rho_QZ)=p_R rho_QZ,
    E_failure(rho_QZ)=(1-p_R)|erase><erase| tensor rho_Z.

p_R=Tr[R U_spatial xi_spatial U_spatial^*] is independent of the qubit.
The successful conditional channel has ZERO information error in the exact
color-symmetric model, also for coherence and entanglement with Z. It is
stronger than distinguishing two fixed inputs. It is NOT a product of an
arrival estimate with a separately favorable readability probability.

B1's >3/5 lower bound is a bound on this joint event, conditional on the
specified prepared geometry and fuel. If S9 is correct, its >1/2 formation
acceptance may be composed by conditional probability after recording and
fuel loading, giving >3/10 per seed trial. Those interventions remain inputs.
No claim here prepares the fuel or input qubit autonomously.

## Spectral components and what is actually bounded

In a fixed m-event B1 record, P contains fuel and co-pairs (dimension m+1),
Q contains ALL separated-charge states. QHQ>=2+Delta-tau, Delta=5/4,
tau=5/96. By min-max, the m+1 low eigenvalues are <=sup PHP<=2+lambda,
and every remaining eigenvalue is >=2+Delta-tau. Thus a genuine cluster
gap exists in this ONE-RESOURCE FIXED-PATH matter matrix. It is not a
restored storage-consistency gap of the full receiver/geometry theory.

For initial fuel |f>, <H>=Omega and Var(H)=lambda^2 exactly. The spectral
measure of every state above the cluster satisfies

    w_high <= lambda^2/(2+Delta-tau-Omega)^2.

At nominal parameters this is 36/27613^2 <4.73e-8. This weight is conserved
by exact evolution. It bounds contributions of the entire high cluster,
not selected individual eigenvectors, and all such components remain in
B1's time-domain bound. This spectral statement uses the same fuel-on law;
it does not import the conversion-off pair-band overlap calculation.

The NOMINAL leading half-line walk has a sine-transform spectral density,
absolutely continuous on [Omega0-2J,Omega0+2J], and no trapped eigencomponent.
Its endpoint return amplitude is J_1(2theta)/theta. This statement is about
the reference walk. The complete charged half-line operator's possible
localized low-band states are not classified in B1 or this packet.

That absence theorem is not needed for the declared finite-time guarantees:
the exact intertwiner residual controls the WHOLE physical state against a
dressed reference vector, including any trapped, returning, or reflected
component. For every successful/fuel projection R in P, RZ=R. The error is
e0+rho*t plus explicit finite-end, parameter and moving-geometry corrections.
No spectral component or reverse amplitude is thrown away. A long-time
irreversible-emission or infinite-volume no-bound-state claim remains open.

The prescribed m=16 control archives all eigenvalues and initial-fuel
residues, their exact finite-volume coherent reconstruction and cluster
weights. Finite-volume states all have point spectra; inverse participation
or endpoint weight alone would not prove permanent trapping. No such
classification is inferred from that control.

## Geometry and return obligations

B1's initial consistency q=a and n_receiver=m-1 are consequences of the
declared vacuum formation/recording block, not a choice of a favorable
receiver eigenvector. Its conditional receiver is retained. The common-
space reflected-prefix reference requires N>=2^28, m>=N/2 and t<=15360;
its boundary drive, fuel-support gates and degree-term discrepancy are
spelled out in the pinned B1 proof/challenge. Those arguments are author-
derived and still need independent review. The new finite spectral control
does not check the enormous receiver or validate that extension.

Every reverse conversion and boundary reflection is already in the fixed-
path matrix and in the Duhamel comparison. An exact finite path can recur
later; neither the channel identity nor its single specified arrival time
implies permanently retained matter or increasing-distance propagation.
