# Final operating claim and complete dependency

Author-derived algebra; inherited comparison proofs remain independently
unreviewed. The two missing source bounds below are NOT assumed proved.

## The claim to decide

For every L in 24N, L>=2016, set T_L=L/(6J), I_L=[L/(8J),T_L].
For every density operator rho on the two ordered fuel-color qubits and
any external reference, prepare both endpoint fuels, matter/field vacuum,
and rho under the frozen physical H_L. Require

    inf_(t in I_L) p_out,L(t)>=1/10,
    inf_(2/J<=t<=T_L) p_q2,L(t)>=1/4,
    successful output trace-distance error<=1/100.     (1)

Output means both fuels empty and co-pairs in the original two receiver
intervals; q=2 permits separated constituents. The information statement
uses the same output event and local ordered-color decoding. The source
comparison covers this L range only, with nominal fixed parameters.

## Dependency chain

The source is a mathematical reference to the original preparation, not
a physical replacement or extra measurement. Arrows are proved comparison
maps, not a sequential intervention protocol.

| Arrow | Exact relation / estimate | Quantifiers and range | Status and error |
|---|---|---|---|
| Original seed to source reference | F_L embeds f tensor f under reflected endpoint ordering | Fixed vacuum/resource preparation; arbitrary color/reference input; all L in family | Exact seed identification; source factorization is NOT asserted for physical evolution |
| Source preparation to actual half-line evolution | v(t)=exp(-ith)f, all |i,j>, i,j>=0 | All t; same nominal jp,jm,kappa,A,J,Omega | Exact definition; source spectral solution missing |
| v to low source u | Pi=1_(-infinity,2.5)(h), u=Pi v, omega=||(1-Pi)f||^2<=36/27613^2 | All t; color tensor identity | Derived B2-SO; sigma=sqrt(omega)<=6/27613; no physical filtering/postselection |
| u to cut source C_M u | M=L/2-2; C_M keeps fuel and max(i,j)<=M | 0<=t<=T_L | Derived uniform tail tau<=1.01 exp(.3-L/105) |
| Cut reference to physical evolution | Phi=E_L(C_Mu tensor C_Mu), ||Psi-Phi||<=delta_L | Same Hilbert space; every color/reference input; entire 0<=t<=T_L | Derived B2-SO; delta_L=12/27613+1500 exp(-L/105)<.001; includes both seed discrepancy and all later coherent returns/interactions |
| Actual source to filtered receiver | A_L=||R_left v||^2, a_L=||R_left u||^2; a_L>=[sqrt(A_L)-sigma]_+^2>=A_L-2sigma | R_left=[L/8,3L/8-1]; t in I_L | Exact projector inequality; full A_L still unbounded below at target |
| Filtered receiver to joint reference | ||R_joint Phi||=a_L; joint reference probability=a_L^2 | Same original outgoing projector; no renormalization | Exact embedding identity; ONE-source 34% is NOT joint 34% |
| Joint reference to physical probability | p_out>=[a_L-delta_L]_+^2 | Entire I_L; all inputs | Exact unitary norm transfer; no independent-success assumption |
| Actual return to retained matter | G=<f,v>, g=<f,u>, |g-G|<=omega; m=1-omega-|g|^2-tau^2 | All source t, cut bound on 0<=t<=T_L | Exact spectral measure and sector identity; late G still uncontrolled |
| Retained matter to physical q=2 | p_q2>=[m-delta_L]_+^2 | 2/J<=t<=T_L | Exact event transfer; no scalar double-fuel-overlap substitution |
| Physical output to information | H_colored=H_scalar tensor I_4; channel on success is p_out rho | All two-qubit inputs/reference states, fixed endpoint law | Existing exact ordering result; conditional error zero if p_out>0; does not supply probability |

## Numerical margin, without changing the physical gates

One sufficient ORIGINAL-source interface remains

    A_L(t)>=17/50 on every I_L,
    |G(t)|<=2/3 on every [2/J,T_L].                    (2)

The source upper/lower bounds themselves are missing. If established,
2sigma<1/2000, delta_L<1/1000 give exactly

    p_out>(17/50-1/2000-1/1000)^2
           =458329/4000000=0.11458225>1/10.           (3)

Using omega<10^-6 and tau^2<10^-6, the retention margin is

    p_q2>[1-2/10^6-(2/3+1/10^6)^2-1/1000]^2
          =24909781540810162540000081 /
             81000000000000000000000000 >1/4.       (4)

Controls verify these rational inequalities and the all-L error constants
at the worst allowed L=2016. They use exp(x)>=sum_(k=0)^64 x^k/k!, x>0,
to upper-bound exp(-x) rationally. Since both error envelopes decrease
with L, no grid of lengths is used. This checks composition only; it
does not enclose A_L or G for the actual source.

## Where the complete proof stops

B2-SS proved |G|<.65 only for 2<=Jt<=16. I_L starts at Jt>=252.
B2-SB shows the scalar summaries do not imply the missing late estimates.
The actual M(z) and spatial resolvent elements in ACTUAL_RESOLVENT.md
therefore have to supply (2), or a stronger direct joint route. Their
real-axis boundary values and uniform oscillatory inversion are the
unresolved arrows. A failure of the 96% central-capture sufficient route
alone would not refute (1).

No eventual-decay theorem, small finite-host eigenspectrum, color quality,
or independently favorable marginal probability closes these arrows.
