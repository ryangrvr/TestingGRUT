# Nonminimal homogeneous reheating equations

Natural units, reduced Planck mass `M_r=2.435e18 GeV`, spatially flat FLRW,
metric convention as in the parent action. For a homogeneous neutral
background the gauge field can be zero when the total gauge current is
zero. This does not omit the relative charged-field phases.

Use canonical real inflaton I, complex fields X,Y,Sigma and

\[
S=\int\!\sqrt{-g}\,[\tfrac12 F R-\tfrac12(\partial I)^2
-\sum_z|D z|^2-\tfrac12m_I^2 I^2-V_h+\mathcal L_r],
\quad F=M_r^2-2\xi |X|^2,\quad\xi=-5/6.
\]

The early potential, with phi held at zero and its tiny low-energy scale
neglected relative to H_i, is

\[
V_h=\lambda_X|X|^4+M_Y^2|Y|^2+\lambda_Y|Y|^4
+\lambda_\Sigma(|\Sigma|^2-h_\Sigma^2/2)^2
+2\operatorname{Re}(\kappa Y^*X^3+\mu\Sigma^*Y^2).
\]

The old hidden parameters are used; no new cross quartic is silently
claimed to vanish at every scale. This is a tree-level matched model.

## Dimensionless evolution

Let `tau=H_i t`, `i=I/M_r`, `j=i'`, `z=(X,Y,Sigma)/H_i`, `v=z'`,
`r=rho_r/(M_r^2 H_i^2)`, `eta^2=H_i^2/M_r^2`, `b=m_I/H_i`,
`gamma=Gamma_I/H_i`, `h_c=H/H_i`, and `U=b^2 i^2/2`.
Primes are tau derivatives. Write the rescaled potential as V and
`F_z=partial V/partial z*`, `K=sum |v|^2`. The code's `H` column is h_c,
not the Higgs field; the constant `HH=40` is h_Sigma/H_i.

\[
\begin{aligned}
f&=F/M_r^2=1-2\xi\eta^2|x|^2,&
f'&=-4\xi\eta^2\operatorname{Re}(x^*v_x),\\
\mathcal R\equiv R/H_i^2
&=\frac{4U-j^2+\eta^2[4V-2K+12\xi|v_x|^2
-12\xi\operatorname{Re}(x^*F_x)]}
{f+12\xi^2\eta^2|x|^2},\\
z''&=-3h_cz'-F_z-\xi\mathcal R x\,\mathbf e_x,&
i''&=-(3h_c+\gamma)j-b^2i,\\
r'&=\gamma j^2-4h_cr,& (\ln a)'&=h_c,\\
f''&=-4\xi\eta^2[|v_x|^2+\operatorname{Re}(x^*x'')],\\
h_c'&=-\frac{j^2+2\eta^2K+4r/3+f''-h_cf'}{2f}.
\end{aligned}
\]

The inflaton friction transfers exactly `gamma j^2` into radiation. The
hidden canonical energy exchanges work with the nonminimal term. Checks:

\[
\begin{aligned}
3fh_c^2+3h_cf'&=\tfrac12j^2+U+\eta^2(K+V)+r\equiv\rho,\\
\mathcal R&=6(2h_c^2+h_c'),\\
\rho'&=-3h_c[j^2+2\eta^2K+4r/3]+\tfrac12 f'\mathcal R,\\
(K+V)'&=-6h_cK-2\xi\mathcal R\operatorname{Re}(x^*v_x),\\
Q_{\rm tot}&=a^3[-2\operatorname{Im}(x^*v_x)
-6\operatorname{Im}(y^*v_y)-12\operatorname{Im}(s^*v_s)]=0.
\end{aligned}
\]

The geometry formula follows by combining the trace, Raychaudhuri equation
and X equation. A separately assembled three-by-three linear system tests
the elimination at 50 random states. H is evolved, not reconstructed from
the Friedmann constraint. Work integrals are extra evolution variables.

Initially `H/H_i=1`, all physical field velocities vanish, X is positive
real and Sigma has supplied phase `-pi/2`. Conditional radial/X-Y stationarity
is solved using the actual initial Ricci scalar; the inflaton amplitude is
then fixed by Friedmann. Relative phases are not minimized. This is not an
inflationary evolution or physical holder construction.

## Charge, yield and diagnostic occupation

The reported quantities are `N_X=-2a^3 Im(x* v_x)` and
`E_X=a^4(|v_x|^2+lambda_X |x|^4)`. The inherited quartic normalization is

\[
\epsilon_Q=\frac{N_X}{2\sqrt{2\lambda_X}
(E_X/(3\lambda_X))^{3/4}}.
\]

Late means are tau-weighted over the last 20% of each run. The radiation
reference is separately compared in the old a-weighted interval 30 to40.
Visible entropy uses the final radiation density and `g*=g*s=106.75`:

\[
Y_X=\frac{N_XH_i^3}{a^3s_{\rm vis}},\qquad
T=(30\rho_r/(\pi^2g_*))^{1/4},\qquad
s_{\rm vis}=2\pi^2g_{*s}T^3/45.
\]

Y diagnostics subtract the algebraic forced solution with instantaneous
Sigma and split the two real Y eigenfrequencies. They do not simulate
decays. The small nonlinear Y correction and an actual time-dependent
decay ledger are not certified by the late occupations. Sigma's radial
diagnostic is likewise not an integrated decaying particle population.

The cold origin criterion during ideal matter domination is
`m_X,vac^2 > 2.5 H^2`. During the early quartic condensate one needs the
coupled perturbation equations instead. In the massless linear
matter-origin control, `H=2/(3t)`, so

\[
\ddot X+2\dot X/t-10X/(9t^2)=0,
\qquad X\propto t^{2/3}\propto a
\]

is the growing solution. This diagnostic was reproduced numerically; it
does not establish instability of the nonlinear populated background.

## Verification limits

The independent Radau comparison ends at tau10 and tests all 19 evolved
variables. Full-length tolerance and endpoint controls use the primary
engine. These are different numerical implementations by the same author,
not an independent research-group replication. There are no finite-k
inhomogeneities, actual decay daughters, loop-corrected profiles or
thermal stress tensors in these runs.
