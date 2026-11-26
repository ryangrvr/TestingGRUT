# Thermodynamic accounting

Work and heat below are in units of $k_BT_b$, and entropies are in nats.
All enabled thermal edges obey
$$
 \frac{w(i\to j)}{w(j\to i)}=e^{E_i-E_j}.
$$
The kinetics in each phase preserve the product Gibbs law, possibly with
disconnected sectors. The response has one occupied sector
$X+F_r=1$, and each forward bond transition spends a fuel molecule with
a reverse refund. The isolated transfer pulses are permutations on
degenerate levels. The contact timing and actuator pulses are supplied
physical controls. An ideal pulse returns its transient drive energy
at the end. Operational friction, clock and drive losses, transport,
fixture return and loading/reset overhead fall under the separately
declared common work allowance.

## Preparation and finite stock

The reference for device reset is the product canonical law
$\gamma_d$ of the bare energies. The reset contacts allow fuel and
material bits to be loaded individually. The reagent remains an
exogenous thermal fair bit, independent of the supplied initial device.
The device begins with one molecule of each fuel, the material unbound,
a fair receptor, a ready orientation and the error bit. The preparation
free energy relative to thermal equilibrium is
$$
 W_{\rm prep}=D_{\rm KL}(p_{d,0}\|\gamma_d)
 =\ln(1+e^{a_s})+\ln(1+e^{a_r})+\ln(1+e^J)
    +2\ln2-h(\eta_{\rm arm}).
$$
This prices the orientation reset and error state, even in controls
that leave those fixtures idle. The swapped-out receptor is retained
in the ledger. One may discard the endpoint state and prepare this
ensemble again, so the full preparation charge is a conservative
sufficient budget within ideal isothermal work accounting. Use either
fresh preparation or credited reset per cycle. The worst preparation
charge sets $h(\eta_{\rm arm})=0$. It is $8.946752\,k_BT_b$ plus
up to $1\,k_BT_b$ for common overhead, for a total allowance of
$9.946752\,k_BT_b$ against the $10\,k_BT_b$ cap. Reset duration is
free in the protocol. The work budget is for the ensemble per trial.

The two fuel stocks are bounded in every trial. The sensor can exhaust
its single fuel and only a reverse transition can restore it within the
acquisition phase. The assembly stock is exhausted whenever $X=1$.
No storing or resetting of measurement results is required by the
autonomous transitions within each supplied contact window. The
bond measurement is an assay of the unconditional output and supplies
no fuel, timing cue or control information to the device.

## Heat, storage and reset

Conservation of energy and the zero modeled bare-energy work during the
phases give
$$
 Q=a_s\,\mathbb{E}[1-F_s]+(a_r+J)\Pr(X=1).
$$
This is heat into the bath, with a net upper bound $7.2\,k_BT_b$
per episode at the registered parameters. Degenerate flips and record
decay can increase entropy while transferring zero bare-energy heat.
For the joint distribution including $B$,
$\Sigma=Q+S_{B,d}^T-S_{B,d}^0\ge0$ follows by Gibbs relative entropy
contraction on each thermal interval. A Gibbs-preserving pulse leaves
that relative entropy invariant. The reference bath supplies rapid
thermalization at the stated temperature. Its heat and entropy intake
are counted. Temperature drift in a finite apparatus has to fit its
rate-error budget.

With final device distribution $p_{d,T}$, ideal joint device reset
without access to the reagent costs
$$
 W_{\rm joint}
 =D_{\rm KL}(p_{d,0}\|\gamma_d)-D_{\rm KL}(p_{d,T}\|\gamma_d)
 =Q+S_d^T-S_d^0
 =\Sigma+I(B;d)_T .
$$
For separate marginal resets the joint correlation free energy is
discarded. The cost is
$W_{\rm sep}=W_{\rm joint}+\sum_i S_{d_i}^T-S_d^T$.
Reset to the same independent input restores both prepared molecules
and the material bit. These work identities include the memory and
its correlations without adding a separate erase-the-bit charge.
The ideal joint result assumes reversible control on the device
registers and no work recovery conditional on source or assay.
The preparation budget avoids relying on extraction from any endpoint
correlation.

| Primary physical arm | Net sensor fuel | Net assembly fuel | Bath heat | Ideal joint reset | Ideal separate reset |
| --- | ---: | ---: | ---: | ---: | ---: |
| Predictive | .485631 | .759906 | 4.648499 | 5.877154 | 6.443609 |
| Disabled | 0 | .478354 | 2.009085 | 2.701295 | 3.393505 |
| Passive record | .485631 | .478354 | 3.465978 | 4.719817 | 5.598676 |
| Inverted | .485631 | .196801 | 2.283458 | 3.456117 | 3.968063 |
| Shuffled | .485631 | .478354 | 3.465978 | 4.666636 | 5.543132 |

Primary full entropy production is $5.477455$ for the predictive arm.
It leaves .399699 nats correlated with the reagent and .566454 nats of
total device correlation. The sensor and actuator errors are physically
populated in every arm. The shuffled control uses a fair $C$, with a
lower actual preparation charge and the same conservative resource cap.
The null ceilings deliberately grant free allocation control and free
latching, so charging those operations could only reduce the null set.
[All ledger entries](RESULTS.json)

The demonstrated advantage is attaining a reliability at a common
work cap. At the primary settings the disabled arm uses less work
and forms fewer bonds than the predictive arm. A realized device must
measure the work of barrier modulation, shielding, transit,
initialization, readout and reset in the full cycle. The $1\,k_BT_b$
allowance is a stated assumption covering overhead in the model.

For a priced objective $S-\lambda W$, suppose the complete null cycle
returns its thermodynamic resources and therefore has nonnegative net work.
Using the primary relaxation plus cue as a bound on the best null score,
the predictive arm strictly wins whenever
$$
0\le\lambda <
 \frac{.7599061758-.6695914360}{9.946751532}
 \simeq .009080
$$
per thermal unit. This is a sufficient guarantee and uses the conservative
work budget for the predictive arm. The optimum priced trade-off outside
that interval is a different task.
