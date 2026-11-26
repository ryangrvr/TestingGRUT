# Analytical feasibility gate

## Prediction and decision

For a symmetric two-state environment, \(\Pr(B_{t+\tau}=B_t)\) solves
\(\dot p=q(1-p)-qp\) with \(p(0)=1\).
A symmetric observation at \(t\) with error \(\epsilon\) therefore gives
\[
 P_{\mathrm{correct}}
 =(1-\epsilon)\frac{1+e^{-2q\tau}}{2}
   +\epsilon\frac{1-e^{-2q\tau}}2
 =\frac12+\left(\frac12-\epsilon\right)e^{-2q\tau}.
\]
The current belief is a sufficient statistic for a Markov source if the
future conditional on that belief has the declared transition law.
An independent thermal record flip rate \(\mu\) during the unshielded
delay, and transfer error \(\eta\), change the decision accuracy to
\[
 p_A=\frac12+\frac{1-2\epsilon}{2}
              e^{-2q\tau-2\mu(\tau-\delta)}(1-2\eta).
\]
Learning about the past can have zero decision value when this contrast
vanishes, when accessible present data already supply the same belief,
or when the best action is common across conditions.

Use \(s_h,s_l\) for bond success conditional on an initial match or mismatch
at catalyst contact. Given a relative source/actuator path, the scalar
birth-death equation gives \(S_X=\pi(1-\exp[-\int_0^{T_r}c\,dt])\).
The relative state flips at \(q+\nu\). The expected exponentials are
computed by a two-state killed generator with killing rates \(c_h,c_l\).
A monotone coupling of relative states starting matched and mismatched
gives \(s_h>s_l\) when \(c_h>c_l\) and \(T_r>0\).
Thus a single orientation chosen using the more probable class is optimal
among locked actions. To meet a reliability \(r_\star\) in that class,
\[
 p_A\ge p_\star=\frac{r_\star-s_l}{s_h-s_l}.
\]
This is informative when \(1/2<p_\star\le1\). A negative acquired
contrast can be flipped when the channel orientation is known.
The required positive contrast \(2p_\star-1\) specifies the retention
condition on the product of sensing, delay and transfer factors above.
For an action selected from observations \(Y\), binary Fano and data
processing require \(I(B;Y)\ge\ln 2-h(1-p_\star)\), where \(h\) is the
binary entropy in nats. The symmetric optimal channel attains the
Hamming rate-distortion curve. This is a task information condition.
Physical work depends on the memory implementation and reset access.

## Optimized matched null and resource condition

For a time-dependent open-loop policy, condition on its independent
randomness and average over the undisturbed balanced environment path.
The rate sum fixes the average integrated rate at \(\overline cT_r\),
or less if the port is closed. Convexity of \(e^{-x}\) yields
\[
 S_{\rm OL}\le \pi(1-e^{-\overline cT_r}).
\]
A half split attains this bound in the relaxed port family. For arbitrary
material-feedback policies, the only possible observation trajectory
prior to first binding is \(X=0\) throughout. Along this no-hit branch the
policy, conditioned on its independent coins, is an open-loop schedule.
The first-hit survival conditional on the environment uses integrated forward
hazard \(\pi\int c\,dt\). The same Jensen argument gives
\[
 S_{\rm FB,end}\le \Pr(\text{hit by }T_r)\le
       1-e^{-\pi\overline cT_r}.
\]
A half split with a perfect immediate latch attains the bound in a relaxation.
This argument grants the null all permissible control of the port without
charging control work, to strengthen the comparison. A new material
event that reveals \(B\) before binding would invalidate the no-hit branch
argument for that device.

For an additional passive cue, use total variation between its joint law
with the exogenous response source path and the product of its marginal
laws. Under the product comparison the cue may randomize the null
schedule. Any fixed policy's score in \([0,1]\) can increase by at most
this total variation. This bound includes a cue at delivery and applies
to an entire extra signal only if its joint path law satisfies the same
independence error bound. A single noisy direct observation of current
\(B\) with balanced symmetric accuracy \(.5+\ell\) has joint TV \(\ell\)
with the response path. Data processing gives the bound and inclusion of
\(B\) at the start of the path gives equality.

For a matched work cap \(W_\star\), a sufficient condition for
thermodynamically budgeted capability is
\[
  1-e^{-\pi\overline cT_r}+\ell < r_\star \le s_l+(s_h-s_l)p_A,\qquad
  W_{\rm prep}^{\rm upper}+W_{\rm common}^{\rm max}\le W_\star .
\]
This establishes a reliability advantage at equal resource caps.
It does not require the sensing arm to spend less than an arm completing
a smaller fraction of assemblies. In a priced objective \(S-\lambda\beta W\),
the relevant difference instead compares the sensing arm with the maximum
over the null work/reliability frontier. For any particular null with the
same fixed architecture the sensing score must exceed its score by more
than \(\lambda\) times the excess work.

## Feasibility certificate before production

For the primary row, conditioning on no thermal flip of \(B\) or \(M\)
in the acquisition interval gives the conservative floor
\[
 L_s=e^{-(q+\mu)t_s}
 \left[\frac12+\frac12\frac{e^{a_s/2}}{e^{a_s/2}+e^{-a_s/2}}
           (1-e^{-s(e^{a_s/2}+e^{-a_s/2})t_s})\right].
\]
Insert \(1-L_s\) for \(\epsilon\) in the decision formula to get a
lower bound \(L_A\). For a matched start, keep only the contribution
of paths with no relative flip, giving
\(L_h=e^{-(q+\nu)T_r}\pi(1-e^{-c_hT_r})\).
For a mismatched start,
\(L_l=\pi(1-e^{-c_lT_r})\) follows from the pointwise rate floor.
They give \(L_l+(L_h-L_l)L_A>.71\) at the primary input. The no-record
upper bound including the declared cue is below \(.67\).
The discard-and-reprepare device budget, even charging \(C\) as pure, is
\(2\ln(1+e^3)+\ln(1+e^{1.2})+2\ln2<9\). The total is below \(10\)
with the specified common allowance. These are bounds on the declared
finite-state model. They are sufficient to justify the bounded replay
and the registered noise and access controls.

## Conventional classification

The feasibility gate closes the additional-principle claim for this
construction. Each constraint follows by ordinary stochastic dynamics,
decision theory and nonequilibrium free energy. Published conventional
comparators cover the mechanisms at issue.

| Question | Conventional result and scope |
| --- | --- |
| Why use a delayed record? | [Still, Sivak, Bell and Crooks (2012)](https://doi.org/10.1103/PhysRevLett.109.120604) distinguish stored predictive information from obsolete information in driven systems. [Becker, Mugler and ten Wolde (2015)](https://doi.org/10.1103/PhysRevLett.115.258103) derive physically implementable temporal predictors in biochemical signaling. |
| How general is the advantage? | [Blackwell (1953)](https://doi.org/10.1214/aoms/1177729032) compares decision experiments by access to more informative data. [Smallwood and Sondik (1973)](https://doi.org/10.1287/opre.21.5.1071) optimize finite-horizon decisions under partial observation. Here the locked binary action is a one-stage Bayesian decision and the feedback null is a bounded partially observed control problem. |
| When can an organism benefit from a changing environment? | [Rivoire and Leibler (2011)](https://doi.org/10.1007/s10955-011-0166-2) quantify decision value under population growth with environmental variation. I3 uses a single finite-time material bond and a different utility. |
| How should acquisition time and fuel be counted? | [Govern and ten Wolde (2014)](https://doi.org/10.1073/pnas.1411524111) identify receptor sampling, readout and energy constraints in cellular sensing. [Malaguti and ten Wolde (2021)](https://doi.org/10.7554/eLife.62574) treat integration time under a changing input. |
| What pays for feedback and correlation storage? | [Sagawa and Ueda (2012)](https://doi.org/10.1103/PhysRevE.85.021104) account for information in nonequilibrium feedback. [Rao and Esposito (2016)](https://doi.org/10.1103/PhysRevX.6.041064) derive chemical-network thermodynamics with free-energy and detailed-balance accounting. |
| What information can a binary reliability target require? | [Shannon (1959)](https://gwern.net/doc/cs/algorithm/information/1959-shannon.pdf) gives lossy coding with a fidelity criterion. The binary error calculation here is the Hamming special case. |

Positive feasibility is conditional on interface restrictions, a finite
response time, a single catalytic port and the work allowance. The small
registered replay tests that conditional physical example and its failures.
It supplies no evidence of an additional constraint on ordinary dynamics.
