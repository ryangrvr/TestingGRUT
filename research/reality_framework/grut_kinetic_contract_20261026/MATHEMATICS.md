# Source, inventories and energy moments

Units are natural. Set a=1 at the end of inflation. Quantities N_A and N_X
in the source file are comoving numbers divided by H_i^3; E_X is conformal
energy divided by H_i^4. There is one longitudinal polarization. Each vector
decay supplies one X and one anti-X, with one bosonic degree per species.

## Exponential decays in the supplied radiation history

Let gamma=Gamma_A/H_i, H(a)=H_i/a^2. For a nonrelativistic parent,

\[
{\cal N}_A(a)={\cal N}_{A,0}
 \exp[-\gamma(a^2-1)/2],\qquad
\frac{d{\cal N}_\pm}{d\ln a}=\gamma a^2{\cal N}_A(a).
\]

Here calligraphic numbers are dimensional comoving densities, so
cal N_A,0=N_A H_i^3. In the parent rest frame,
p_star=sqrt(m_A^2/4−m_X^2). Define q=a p, q0=p_star/sqrt(gamma), x=q/q0.
Up to the negligible lower cutoff a=1 for this late decay,

\[
\frac{d{\cal N}_\pm}{dq}
 ={\cal N}_{A,0}\frac{q}{q_0^2}e^{-q^2/(2q_0^2)},\qquad
\langle q\rangle=q_0\sqrt{\pi/2},\qquad
\langle q^2\rangle=2q_0^2.
\]

The exact finite-start expression has q≥p_star and an additional
exp(gamma/2) factor. Here gamma≈1.6e-21 at e=1e-7, so the distinction is
irrelevant to the reported precision. The earlier parent-mode momenta are
negligible compared with m_A a over the appreciably populated decay interval.
Even at the upper band edge, the estimated fraction decayed before parent
p/m_A falls below .1 is less than 1e-10 for both supplied values of e.
Do not use this rest-parent source if that hierarchy ceases to hold.

The phase-space convention is n_i=integral dp p^2 f_i/(2 pi^2).
The physical source for each daughter is

\[
(\partial_t-Hp\partial_p)f_\pm\big|_A
 =\frac{2\pi^2}{p^2}\Gamma_A n_A\delta(p-p_\star).
\]

Its integrated number moment is Gamma_A n_A for each sign, and its total
daughter energy moment is Gamma_A m_A n_A. At fixed comoving momentum,

\[
\partial_t f_\pm(q,t)\big|_A
 =\frac{2\pi^2\Gamma_A {\cal N}_A}{q^2}
   \delta(q-a p_\star).
\]

After the source has finished and before scattering,
f_pm(q)=2 pi^2 q^(-2) d cal N_pm/dq. This distribution is not an equilibrium
gas. For these hard daughters the conformal radiation moment gives

\[
R_A=\frac{2N_A\langle q\rangle}{H_i E_X}
 \simeq\frac{N_A(m_A/H_i)}{E_X}
       \sqrt{\frac{\pi H_i}{2\Gamma_A}}.
\]

The first expression uses the exact daughter momentum; the second replaces
2p_star by m_A. They agree to 1.6e-15 for the stated inputs. A kinetic solver
uses the exact energy moment when daughters become nonrelativistic.

The source's unprocessed extra rest mass relative to the inherited dark
matter normalization is f_rest=2 m_X Y_Q N_A/(DM_per_s |N_X|). It is .0015025,
whereas the inherited thermal replay gives f=.67368. Neither value by itself
is the full kinetic answer. Elastic scattering can redistribute energy;
inverse inelastic reactions are needed to turn it into extra particle number.

## Fast inverse decay carries a pair inventory

Use species order (X, anti-X, s) and stoichiometric columns

\[
\nu_{s\to X\bar X}=(1,1,-1),\quad
\nu_{ss\to X\bar X}=(1,1,-2),\quad
\nu_{2X+2\bar X\to X+\bar X}=(-1,-1,0).
\]

The charge vector (1,−1,0) annihilates every column. The inventory
P=n_-+n_s is conserved by inverse decay alone, and decreases by one per
forward ss reaction or forward 4→2 event. When the scalar eventually
decays, its population contributes to the surviving pair count.

For a fixed-temperature Maxwell–Boltzmann fast-inverse-decay approximation,

\[
\mu_s=\mu_++\mu_-,\qquad n_s=K(T)n_+n_-,\qquad
\left.\frac{\partial P}{\partial n_-}\right|_{T,n_Q}
 =1+K(n_++n_-),\quad n_+=n_Q+n_-.
\]

This is the denominator already present in the reviewed s-channel rate.
For a changing temperature and net gas charge the complete differential is

\[
dP=[1+K(n_++n_-)]dn_-+K n_-dn_Q+K'(T)n_+n_-dT.
\]

It must be solved together with the energy equation, rather than retaining
only its fixed-temperature derivative. Condensed objects add their charge
and number variables. For Bose gases the simple MB product is approximate;
an explicit s distribution or the correct Bose susceptibilities can replace it.

For any on-shell channel the gain-minus-loss integrand includes
product_in f product_out(1+f) minus product_out f product_in(1+f).
It vanishes for a common-temperature Bose state when the channel conserves
energy and sum mu_in=sum mu_out. Direct checks use s decay, ss conversion
and a charge-conserving 4→2 process. These do not determine thermal rates.

## Energy and the finite object population

For a homogeneous isolated sector after the sources finish,

\[
\dot\rho_h+3H(\rho_h+P_h)=0.
\]

With parent decays the gas receives Gamma_A m_A n_A, and the parent loses
the identical energy. The combined equation still conserves energy, with
the appropriate parent pressure. Portal transfer, if retained, appears
with opposite signs in the two sectors. At fixed medium masses, the sum of
all internal collision energy moments vanishes. If masses depend on a
mean field or temperature, its energy/work must be included consistently.

For a monodisperse zero-temperature population rho_L=n_L E(Q),

\[
d\rho_L=\omega\,d(n_LQ)+(E-Q\omega)\,dn_L,
\qquad\omega=\frac{dE}{dQ}.
\]

At Q=200000 the archived E/Q=.9233161354 m_X and omega=.8611951308 m_X
give E−Q omega=1148.1075 GeV per object. Their difference is physical.
The second term matters when objects form, merge or fragment. In a
distributed population evolve n_L(Q,t), include transport in Q and the
energy of every formation/capture event; do not infer the distribution from
the single cold reference object. The derivative test in verify.py uses
an arbitrary smooth test energy, not a new fit to the objects.

In the proposed erasure channel anti-X+L_Q→X+L_(Q−2), the emitted particle
has K_out=K_in+E(Q)−E(Q−2), approximately K_in+2 omega. The extra rest-mass
inventory falls but this energy remains in the hidden gas unless an explicit
escape channel carries it away. Its nonthermal source and inverse process
belong in a kinetic treatment of the erasure premise.
