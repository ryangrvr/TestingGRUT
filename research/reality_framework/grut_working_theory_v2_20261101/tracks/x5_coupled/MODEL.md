# Model and work identity

Use cosmic time, $H_i=1$ units, a flat homogeneous metric and a fixed
comoving radial cutoff $Q$. The fields are $\sigma,\chi,s$. The four
Majorana species use the archived H55 masses and phase scales. The
phase $\theta_j=\phi_{i(j)}/f_j$ depends on the parent $\sigma$ for N1
and $\chi$ for N2–N4. The spectator potential is
$V_s=g_s^2\chi^2s^2/2$. Its contribution to both scalar equations is
obtained from this same potential.

## States and normalization

For each helicity $h=\pm1$ and comoving momentum $q$, let $p=q/a$.
Use the unrotated two-level Hamiltonian

$$
\mathsf H_h^{\rm u} =
M\cos\theta\,\sigma_x-M\sin\theta\,\sigma_y+hp\,\sigma_z .
$$

The state has Bloch vector $r$. Rotating the state by the inverse
chiral phase gives the rotation axis $b=(M,0,hp+\dot\theta/2)$ and
$\dot r=2 b\times r$. The rotation is a change of numerical basis
within the same finite regulator. The physical mass axis in that basis
is $b_0=(M,0,hp)$ and $E=|b_0|$.

Write $\sum w$ for $\sum_h\int_0^Q q^2dq/(2\pi^2)$. The subtracted
static vacuum contributes zero to the following stress and force.
The factor one-half multiplying the two-level trace uses Majorana
counting. With $f=(1+b_0\cdot r/E)/2$ there is one particle occupation
per helicity.

$$
\rho_N=\frac{\sum w\,E f}{a^3},\quad
P_N=\frac{\sum w\,hp(r_z+hp/E)}{6a^3},\quad
J_i=-\sum_{j\to i}\frac{\sum w\,M_j r_y}{2 f_j a^3}.
$$

Here $f_j$ in the force is the model's phase scale. The source has
$\dot\rho_N+3H(\rho_N+P_N)=\sum_i J_i\dot\phi_i$.
The scalar equation carries $-J_i$. Pressure uses the full two-level
state. Away from an adiabatic slice it generally differs from the
kinetic pressure obtained by dropping the coherences in the physical
mass basis.

## Second-order subtraction

With $n=\mathsf H^{\rm u}/E$ interpreted as a unit Bloch axis, the
first adiabatic correction to the negative vacuum is
$\delta r=(n\times\dot n)/(2E)$. Normalization then gives the second
order physical energy per helicity

$$
\epsilon^{(2)}=\frac{|\dot n|^2}{16E}
 =\frac{M^2\dot\theta^2}{16E^3}+\frac{M^2 H^2p^2}{16E^5}.
$$

Define $Z_i(a)$ by summing $M^2/(8 f_j^2 E^3)$ over species with
parent $i$ and the physical measure $\sum w/a^3$. Define $Z_g(a)$ by
summing $M^2p^2/(8 E^5)$ over all species with the same measure. In the
continuum radial integrals through $p_{\max}=Q/a$, the primitives are
$\operatorname{asinh}(p_{\max}/M)-\beta$ and
$\operatorname{asinh}(p_{\max}/M)-\beta-\beta^3/3$, respectively,
where $\beta=p_{\max}/\sqrt{p_{\max}^2+M^2}$.

The finite prescription subtracts these whole second-order terms,
including their finite parts, from the homogeneous action. The
classical unit kinetic coefficients and $M_P$ are the prescribed
matched parameters. The implemented matter kinetic coefficient is
$A_i(a)=1-Z_i(a)$. Its contribution to density is $A_i\dot\phi_i^2/2$
and to pressure is
$A_i\dot\phi_i^2/2+(\partial_{\ln a}A_i)\dot\phi_i^2/6$. The same $A_i$
enters the canonical momentum $\Pi_i=a^3A_i\dot\phi_i$. The gravity
term moves the constraint denominator from $3M_P^2$ to
$3M_P^2+Z_g/2$. The logged matter stress includes the scalar subtraction.
The gravitational subtraction is recorded separately. Its constraint
contribution follows by varying the lapse in the homogeneous action.

At leading derivative order, expanding the two helicities of the
rotated force gives
$J_i^{(2)}=Z_i\ddot\phi_i+H(3Z_i+\partial_{\ln a}Z_i)\dot\phi_i$.
The scalar inertia and friction corrections therefore accompany the
energy subtraction. This scheme maintains the finite-model work
identity. The comoving cutoff and the choice of finite terms are part
of the calculation. A covariant stress calculation requires matching
the additional metric counterterms and their finite parts.

## Discrete exchange

At fixed scale factor, a midpoint velocity $\bar v_i$ gives a constant
rotated axis for each mode. Its unitary step conserves $b\cdot r$.
Consequently the change of normal-ordered comoving fermion energy per
mode is $-\bar v_i w\,\Delta r_z/(4 f_j)$, exactly the opposite of the
work of the scalar impulse $w\,\Delta r_z/(4 f_j)$. The update is

$$
\Delta\Pi_i =
 -a^3\Delta t\,\overline{\partial_i V}
 +\sum_{j\to i}\frac{\sum w\,\Delta r_z}{4 f_j}.
$$

The discrete gradient contracts with the field displacement to give the
change of the potential. This also gives an exact finite-regulator shift
charge $\Pi_i-\sum_{j\to i}\sum w\,r_z/(4 f_j)$ when the potential is
independent of that parent and collisions are absent.

Symmetric expansion substeps preserve the canonical momenta and spin
state while integrating the Friedmann constraint. The pressure-work
integral is evaluated separately. Spectator damping transfers its
kinetic energy loss to dark radiation. The optional fermion decay
control uses amplitude damping in the instantaneous free-mass basis at
rate $\Gamma_{\rm rest}M/E$ and transfers that energy to the visible
radiation variable. The latter model can damp reversible polarization
during strong phase motion. It is a sensitivity test for the
collisionless result and leaves the physical collision kernel open.
