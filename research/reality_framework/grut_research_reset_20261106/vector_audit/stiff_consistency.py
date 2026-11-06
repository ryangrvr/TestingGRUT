"""Matched smooth histories and gravitational excitation of hard modes."""
import argparse
import json
from math import comb
import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import betainc, hankel2, h2vp

from inflation import Background, FINE
from transfer import A_S, DM_PER_ENTROPY, HISTORIES, M, MPL, RESULTS, write

LATE = HISTORIES["late_stiff"]
P_STIFF = LATE.p
NAMES = {"narrow": (.06, 24.), "wide": (.15, 12.)}
FIT = json.loads((RESULTS / "late_stiff_extension.json").read_text())["result"]["fine"]
H_END = FIT["H_end_gev"]
H_R = LATE.h_reheat*M
N_R = np.log(H_END/H_R)/P_STIFF
X_END = (M/H_END)**(1/P_STIFF)
Q_END = X_END*H_END/M
K_MIN = .001
K_MAX = 64.
TENSOR_PIVOT = FIT["tensor_ratio"]*A_S


class Smooth:
    def __init__(self, name):
        self.width, self.rate = NAMES[name]
        self.inflation = Background(FIT["alpha"], FINE.background_tolerance)
        b = np.sqrt(2/(3*FIT["alpha"]))
        log_h, eps, _ = self.inflation.track(0.)
        self.log_h0 = log_h
        phi = -np.log1p(-np.sqrt(3-eps)*np.exp(log_h))/b

        def derivative(_n, y):
            speed = y[1]
            return [speed, -(3-speed**2/2)*(speed+2*b/np.expm1(b*y[0]))]

        post = solve_ivp(derivative, (0., self.width), [phi, -np.sqrt(2*eps)],
                         method="DOP853", rtol=2e-13, atol=2e-14, dense_output=True)
        self.post = post.sol
        self.b = b
        self.end = brentq(lambda n: self.background(n)[1]-1,
                           -self.width, self.width)

    def inflation_values(self, n):
        if n <= 0.:
            log_h, eps, _ = self.inflation.track(n)
            return log_h-self.log_h0, eps
        phi, speed = self.post(n)
        eps = .5*speed*speed
        log_h = np.log(-np.expm1(-self.b*phi))-.5*np.log(3-eps)
        return log_h-self.log_h0, eps

    def background(self, n):
        """Return log(H/H_end) and epsilon for the conserved background."""
        if n <= -self.width:
            return self.inflation_values(n)
        u = n-N_R
        stiff = -P_STIFF*n
        if n < self.width:
            infl, eps_inf = self.inflation_values(n)
            z = (n+self.width)/(2*self.width)
            blend = betainc(5., 5., z)
            slope = 630*(z*(1-z))**4/(2*self.width)
            return ((1-blend)*infl+blend*stiff,
                    (1-blend)*eps_inf+blend*P_STIFF-slope*(stiff-infl))
        softplus = np.logaddexp(0., self.rate*u)
        fraction = np.exp(-np.logaddexp(0., -self.rate*u))
        return stiff+(P_STIFF-2)*softplus/self.rate, P_STIFF-(P_STIFF-2)*fraction

    def horizon(self, n):
        return np.exp(n+self.background(n)[0])

    def potential(self, n):
        log_h, eps = self.background(n)
        return np.exp(2*(n+log_h))*(2-eps)

    def basis(self, k, n):
        u = k/((P_STIFF-1)*self.horizon(n))
        order = 1/(P_STIFF-1)-.5
        amp = np.sqrt(np.pi*u/2)
        h = hankel2(order, u)
        return amp*h, (P_STIFF-1)*amp*(.5*h+u*h2vp(order, u))

    def vacuum(self, k, n):
        def frequency(time):
            d = .001
            u = self.potential(time)
            left, right = self.potential(time-d), self.potential(time+d)
            b = self.horizon(time)
            u_second = b*b*((right-2*u+left)/d**2+
                            (1-self.background(time)[1])*(right-left)/(2*d))
            return k-u/(2*k)+(u_second-u*u)/(8*k**3)

        w = frequency(n)
        d = .001
        slope = (frequency(n+d)-frequency(n-d))/(2*d)
        amp = np.sqrt(k/w)
        return [complex(amp), amp*(-1j*w/self.horizon(n)-slope/(2*w))]


def hard_mode(model, k, start_ratio, tolerance, times):
    end = model.width+.025
    start_ceiling = -model.width-.05
    target = min(k/start_ratio, model.horizon(start_ceiling))
    start = brentq(lambda n: model.horizon(n)-target, -26., start_ceiling,
                   xtol=1e-12)

    def derivative(n, state):
        v, speed = state
        h, eps = model.background(n)
        return speed, -(1-eps)*speed-(k*k*np.exp(-2*(n+h))-2+eps)*v

    sol = solve_ivp(derivative, (start, end), model.vacuum(k, start),
                    method="DOP853", rtol=tolerance, atol=tolerance/20,
                    dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    v, speed = sol.y[:, -1]
    f, df = model.basis(k, end)
    norm = 2j*k/model.horizon(end)
    beta = (f*speed-df*v)/norm
    alpha = (df.conjugate()*v-f.conjugate()*speed)/norm
    occ = float(abs(beta)**2)
    # Below the WKB start the excited state difference is omitted.
    stress = np.zeros((2, len(times)))
    for index, n in enumerate(times):
        if n < start:
            continue
        b = model.horizon(n)
        if n >= end and k/b > 40.:
            stress[:, index] = (4*occ, 4*occ/3)
            continue
        if n <= end:
            value, slope = sol.sol(n)
        else:
            f, df = model.basis(k, n)
            value = alpha*f+beta*f.conjugate()
            slope = alpha*df+beta*df.conjugate()
        field = abs(value)**2
        kinetic = abs((slope-value)*b/k)**2
        # Local order-0/2/4 adiabatic subtraction with an IR window.
        # The window dependence is a local stress diagnostic.
        z = (b/k)**2
        d = .001
        eps = model.background(n)[1]
        u = model.potential(n)
        un = (model.potential(n+d)-model.potential(n-d))/(2*d)
        unn = (model.potential(n+d)-2*u+model.potential(n-d))/d**2
        a = -u/2
        a_n = -un/2
        c = (b*b*(unn+(1-eps)*un)-u*u)/8
        z_k = 1/k**2
        window = np.clip((k/b-3.)/3., 0., 1.)**2
        e_ad = 2+window*(z+(a*a+b*b*(a_n-a))*z_k**2)
        p_ad = (2/3+window*((4*a/3+b*b)*z_k+
                             (4*c/3-a*a/3+b*b*(a_n-a))*z_k**2))
        stress[:, index] = (kinetic+field-e_ad, kinetic-field/3-p_ad)
    wronskian = abs(abs(alpha)**2-occ-1)/max(1., occ)
    return occ, stress*k**4, float(wronskian)


def potential_derivatives(model, n):
    """Differentiate the background ODE and the matching polynomial analytically."""
    from sympy import lambdify, symbols

    velocity, slope, b = symbols("velocity slope b")
    velocity_n = -(3-velocity**2/2)*(velocity+slope)
    slope_n = -slope*(slope/2+b)*velocity
    eps = velocity**2/2
    expressions = [eps]
    for _ in range(4):
        eps = eps.diff(velocity)*velocity_n+eps.diff(slope)*slope_n
        expressions.append(eps)
    evaluate = lambdify((velocity, slope, b), expressions, "numpy", cse=True)

    in_inflation = n < model.width
    h_inf, eps_inf = np.array([model.inflation_values(x) for x in n[in_inflation]]).T
    field_factor = np.sqrt(3-eps_inf)*np.exp(h_inf+model.log_h0)
    eps_jets = evaluate(-np.sqrt(2*eps_inf), 2*model.b*(1-field_factor)/field_factor,
                        model.b)
    log_h = np.zeros((6, len(n)))
    log_h[0] = -P_STIFF*n
    log_h[1] = -P_STIFF
    infl = np.array([h_inf, *[-value for value in eps_jets]])
    correction = infl-log_h[:, in_inflation]
    z = np.clip((n[in_inflation]+model.width)/(2*model.width), 0., 1.)
    rise = np.polynomial.Polynomial([0., 0., 0., 0., 630., -2520.,
                                     3780., -2520., 630.])
    weight = [1-betainc(5., 5., z)]
    for order in range(1, 6):
        value = -rise.deriv(order-1)(z)/(2*model.width)**order
        weight.append(np.where(z > 0, value, 0.))
    for order in range(6):
        log_h[order, in_inflation] += sum(comb(order, j)*weight[j]*correction[order-j]
                                         for j in range(order+1))

    ell = log_h.copy()
    ell[0] += n
    ell[1] += 1.
    horizon = np.exp(ell[0])
    horizon_sq = [horizon**2]
    for order in range(1, 5):
        horizon_sq.append(sum(2*comb(order-1, j-1)*ell[j]*horizon_sq[order-j]
                              for j in range(1, order+1)))
    mass = log_h[1:].copy()
    mass[0] += 2.
    u, d1, d2, d3, d4 = [
        sum(comb(order, j)*mass[j]*horizon_sq[order-j] for j in range(order+1))
        for order in range(5)]
    l1, l2, l3 = ell[1:4]
    u1 = horizon*d1
    u2 = horizon**2*(d2+l1*d1)
    u3 = horizon**3*(d3+3*l1*d2+(l2+2*l1*l1)*d1)
    u4 = horizon**4*(d4+6*l1*d3+(4*l2+11*l1*l1)*d2
                     +(l3+7*l1*l2+6*l1**3)*d1)
    return horizon, u, u1, u2, u3, u4


def uv_envelope(model, cutoff, mesh_scale=1):
    """Bound the excitation using both absolute and oscillatory WKB estimates."""
    n = np.unique(np.r_[np.linspace(-15, -model.width, 1300*mesh_scale+1),
                        np.linspace(-model.width, model.width, 2000*mesh_scale+1),
                        np.linspace(model.width, 6., 1600*mesh_scale+1)])
    # Outside this interval the inflation and stiff residuals decay
    # with at least the third power of the comoving horizon.
    b, u, u1, u2, u3, u4 = potential_derivatives(model, n)
    a, a1, a2 = -.5*u, -.5*u1, -.5*u2
    c, c1, c2 = (u2-u*u)/8, (u3-2*u*u1)/8, (u4-2*(u1*u1+u*u2))/8
    z = 1/cutoff**2
    w_min = 1-abs(a)*z-abs(c)*z*z
    if np.min(w_min) <= 0.:
        raise ValueError("The WKB ultraviolet lower frequency bound failed.")
    residual = (abs(a2)*(abs(a)+abs(c)*z)/(2*w_min)+2*abs(a*c)
                +c*c*z+abs(c2)/(2*w_min)+3*(abs(a1)+abs(c1)*z)**2/(4*w_min**2))
    # Integration in conformal time, with a factor two safety margin.
    c5 = 2*simpson(residual/(2*w_min*b), x=n)
    at_cut = c5/cutoff**5
    factor = np.sinh(at_cut)/at_cut
    absolute = (factor*c5)**2/(6*cutoff**6)

    # For W2=sqrt(k²-U), integrate the leading residual excitation by
    # parts once. The remaining exact Bogoliubov terms cost at most 3 J²
    # when J=integral |residual/(2 W2)| is small.
    w2_min = np.sqrt(1-max(abs(u))*z)
    j3 = 2*simpson((abs(u2)/(8*w2_min**3)
                    +5*u1*u1*z/(32*w2_min**5))/b, x=n)
    if j3/cutoff**3 >= .1:
        raise ValueError("The oscillatory ultraviolet remainder is outside its bound.")
    oscillatory4 = 2*simpson((abs(u3)/(16*w2_min**4)
                              +9*abs(u1*u2)*z/(32*w2_min**6)
                              +15*abs(u1)**3*z*z/(64*w2_min**8))/b, x=n)
    oscillatory4 += 3*j3*j3*z
    oscillatory = oscillatory4**2/(4*cutoff**4)
    return dict(absolute_energy=float(absolute), oscillatory_energy=float(oscillatory),
                coherent_pressure_integral=float(factor*c5/cutoff),
                fourth_order_residual=float(c5))


def soft_vector(model, count, tolerance):
    q = np.geomspace(.0005, 128., count)
    q2 = q*q
    n = np.r_[np.linspace(-16., 2.5, 151), np.log(16.), np.log(24.)]
    x0 = np.exp(n[0])
    small = q2*x0**(2*(P_STIFF-1))
    xh = x0**(1-P_STIFF)
    kinetic = x0*x0/(x0*x0+q2)
    y0 = np.r_[1-small/(2*(P_STIFF-1)*(P_STIFF+1)),
               -xh*kinetic*small/(P_STIFF+1)]

    def derivative(nx, data):
        x = np.exp(nx)
        h = H_END/M*np.exp(model.background(nx-np.log(X_END))[0])
        y, momentum = data.reshape(2, -1)
        kinetic = x*x/(x*x+q2)
        return np.r_[momentum/(x*kinetic*h), -x*y/h]

    sol = solve_ivp(derivative, (n[0], n[-1]), y0, method="DOP853", t_eval=n,
                    rtol=tolerance, atol=tolerance/50)
    if not sol.success:
        raise RuntimeError(sol.message)
    x = np.exp(n)[:, None]
    h = H_END/M*np.exp([model.background(z-np.log(X_END))[0] for z in n])
    amplitude, momentum = np.split(sol.y.T, 2, axis=1)
    kinetic = x*x/(x*x+q2)
    potential_energy = x*x*amplitude**2
    kinetic_energy = momentum**2/kinetic
    # The pivot spectrum gives an upper envelope for this later-exiting band.
    prefactor = TENSOR_PIVOT/(48*x[:, 0]**4*h*h)
    energy = prefactor*simpson(q2*(kinetic_energy+potential_energy), x=np.log(q), axis=1)
    pressure = prefactor*simpson(q2*((3-2*kinetic)*kinetic_energy-potential_energy)/3,
                                x=np.log(q), axis=1)
    c = 1-kinetic
    b = np.sqrt(kinetic)
    canonical = b*amplitude
    derivative = momentum/b+x*h[:, None]*c*canonical
    frequency = np.sqrt(x*x+q2+3*(x*h[:, None])**2*c*(1-c))
    weight = q2*(derivative**2+frequency**2*canonical**2)/frequency
    return q, weight[-1], n, energy, pressure


def fit_frozen_power(q, model):
    from inflation import LOG_PIVOT, mode
    log_k = np.log(LATE.k_star*q)
    offset = LOG_PIVOT+FIT["pivot_efolds"]-model.inflation.track(-FIT["pivot_efolds"])[0]
    nodes = np.linspace(log_k[0], log_k[-1], 41)
    values = []
    for node in nodes:
        n = model.inflation.exit(node-offset)
        log_h = model.inflation.track(n)[0]
        values.append(2*log_h+np.log(mode(model.inflation.track, n, False, FINE)))
    return 2*FIT["V0_over_Mpl4"]/np.pi**2*np.exp(CubicSpline(nodes, values)(log_k))


def run(name, fine):
    model = Smooth(name)
    times = np.unique(np.r_[np.linspace(-14., -.35, 37),
                            np.linspace(-.35, .45, 101), np.linspace(.45, 3., 31)])
    k = np.geomspace(K_MIN, K_MAX, 257 if fine else 129)
    start_ratio, tolerance = (92., 8e-12) if fine else (46., 4e-10)
    rows = [hard_mode(model, w, start_ratio, tolerance, times) for w in k]
    occ = np.array([row[0] for row in rows])
    stress = np.array([row[1] for row in rows])
    integral = simpson(k**4*occ, x=np.log(k))
    above_end = k >= 1.
    fraction_above = simpson((k**4*occ)[above_end], x=np.log(k[above_end]))/integral
    end_rho = H_END**2/(24*np.pi**2*MPL**2)
    bg = np.array([model.background(n) for n in times])
    stress_ratio = end_rho*np.exp(-4*times-2*bg[:, 0])*simpson(
        stress, x=np.log(k), axis=0)
    q, weight, nx, soft_energy, soft_pressure = soft_vector(
        model, 641 if fine else 321, 3e-10 if fine else 4e-9)
    q_tail = np.geomspace(q[-1], 1e8, 97)[1:]
    power, coefficient = LATE.uv_power
    full_q = np.r_[q, q_tail]
    weighted = np.r_[weight, coefficient*q_tail**(-power)]*fit_frozen_power(full_q, model)
    prefactor = M*M*MPL*MPL/(16*LATE.entropy_reference*DM_PER_ENTROPY)
    density = prefactor*(simpson(weighted, x=np.log(full_q))+.5*weighted[0])
    sample = np.linspace(-model.width, model.width, 3001)
    eps = [model.background(n)[1] for n in sample]
    file = RESULTS/f"consistency_{name}_{'fine' if fine else 'coarse'}.npz"
    np.savez_compressed(file, k_end_units=k, occupation=occ, early_efolds=times,
                        early_hard_energy_over_background=stress_ratio[0],
                        early_hard_pressure_over_background=stress_ratio[1],
                        soft_vector_log_x=nx, soft_energy_over_background=soft_energy,
                        soft_pressure_over_background=soft_pressure, q=full_q,
                        vector_tensor_weight=weighted)
    return dict(blend_half_width=model.width, reheating_rate=model.rate,
                epsilon_end_blend_min=float(min(eps)), epsilon_end_blend_max=float(max(eps)),
                physical_end_shift_efolds=model.end,
                physical_pivot_to_end_efolds=FIT["pivot_efolds"]+model.end,
                hard_integral=float(integral), hard_fraction_above_original_end=float(fraction_above),
                max_relative_wronskian_error=max(row[2] for row in rows),
                max_early_hard_energy_abs_fraction=float(max(abs(stress_ratio[0]))),
                max_early_hard_pressure_abs_fraction=float(max(abs(stress_ratio[1]))),
                max_soft_energy_fraction=float(max(soft_energy)),
                max_soft_pressure_abs_fraction=float(max(abs(soft_pressure))),
                vector_rest_density_over_frozen_target=float(density),
                radiation_asymptote_H_ratio=float(np.exp(model.background(N_R+3)[0]+
                                                           6+P_STIFF*N_R)),
                soft_data_file=file.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history", choices=NAMES)
    parser.add_argument("--fine", action="store_true")
    args = parser.parse_args()
    result = run(args.history, args.fine)
    label = "fine" if args.fine else "coarse"
    write(RESULTS/f"consistency_{args.history}_{label}.json", result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
