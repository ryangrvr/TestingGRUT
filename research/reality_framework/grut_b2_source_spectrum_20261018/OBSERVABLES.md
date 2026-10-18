# What exactly is being bounded?

The physical gate and B2-SO comparison remain unchanged. All all-size
mathematics is author-derived pending independent review.

## One source

The exact conserved resource-one space is

    H_source=(span{f} direct-sum ell^2{(i,j):i,j>=0}) tensor C^2_color.

Fuel f is the ONLY spatial/electric configuration with occupied fuel:
resource conservation forces zero matter, hence Gauss forces zero field.
All other spatial configurations contain one plus and one minus, possibly
separated, with empty fuel. Source h=h_spatial tensor I_2. Define

    V(t)=exp(-ith_spatial),   G(t)=<f,V(t)f>,
    P_F=|f><f| tensor I_2.

Then the complete fuel-return amplitude map for the DECLARED fuel input is

    P_F exp(-ith)(|f> tensor I_2)=|f> tensor [G(t)I_2].  (1)

Its operator norm is |G(t)|. For every input color density matrix, including
a purification in an inaccessible reference, total fuel occupation is
EXACTLY |G(t)|^2. This is not conditional on receiver arrival, low energy,
or a measurement. It is not a bound uniform over arbitrary spatial seeds.

The other exact unfiltered source observable is

    A_L(t)=||P_co,Rleft V(t)f||^2,
    P_co,Rleft=sum_(x=L/8)^(3L/8-1)|x,x><x,x|.         (2)

It is a ONE-source probability, with the actual fuel preparation and
physical co-location detector, not the full joint outgoing probability.

## Two physical sources

The full resource-two H_L has a much larger fuel-occupied subspace. For
example one chosen fuel is occupied in one double-fuel configuration and
in L^2 one-pair configurations with the other fuel empty. Its scalar rank
is 1+L^2, versus rank one for the single source. Its marginal occupation
cannot be obtained by squaring <FF,exp(-itH_L)FF>. B2's 70.46% marginals
and 48.89% joint double-fuel occupation remain distinct.

B2-SO's comparison uses u=V(t)Pi f, UNNORMALIZED, cuts it before the center,
and embeds the product as Phi_L. If a_L is its ONE-source receiver
probability, then

    ||R_L Phi_L||^2=a_L^2, NOT a_L.

The normalization obeys ||Phi_L||<=1. For normalized physical Psi_L and
||Psi_L-Phi_L||<=delta, reverse triangle gives

    p_out>= [a_L-delta]_+^2,
    |p_out-a_L^2|<=delta(1+||Phi_L||)<=2delta.          (3)

Thus the general 2delta estimate is applicable even to this unnormalized
comparator, but its reference event is a_L^2. It cannot turn a one-source
34% condition into a joint physical probability above 33.8%.

With s=6/27613, a_L>=A_L-2s, 2s<1/2000, delta<1/1000,
the sufficient A_L>=17/50 gives the CORRECT joint bound

    p_out>(17/50-1/2000-1/1000)^2
          =458329/4000000=0.11458225>1/10.             (4)

This multiplication is only inside B2-SO's explicitly controlled reference;
it is not an independent-source assumption for the physical evolution.
All two-source returns and interactions are retained in its state error.
No such lower bound is earned until the source estimate is proved.

## Frozen quantifiers

The sufficient source questions are:

* Every L in 24N, L>=2016, and EVERY t in I_L=[L/(8J),L/(6J)]: A_L(t)>=17/50.
* Every same L and EVERY t in [2/J,L/(6J)]: |G(t)|<=2/3.

Color input is arbitrary and may be entangled with a reference. Parameters
are the single nominal set, with no disorder robustness claim. H_L, its
seed, R_L, and all physical acceptance thresholds are unchanged.

Proving a shorter fuel-return interval is a partial source result. Finding
a violation of one sufficient condition would not alone refute p_out>=1/10
or q=2 retention by every possible proof. No violated source inequality
or operating counterexample has been established in this cycle.
