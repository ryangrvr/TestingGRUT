"""Fixed-input background and mode replication. Imports no archive routines."""
from argparse import ArgumentParser
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import brentq
from scipy.special import expit, gamma, jv, jvp, yv, yvp

ROOT = Path(__file__).resolve().parent
FROZEN = json.loads((ROOT / "frozen_inputs.json").read_text())
UNITS = FROZEN["units"]
MASS = UNITS["vector_mass_GeV"]
PLANCK = UNITS["Mbar_GeV"]
MPC = UNITS["GeV_to_Mpc_inverse"]
PIVOT = UNITS["k_pivot_Mpc_inverse"]
TARGET = UNITS["target_DM_energy_over_entropy_GeV"]
G_RHO = UNITS["g_energy"]
G_S = UNITS["g_entropy"]
S_TODAY = (2*np.pi**2/45)*UNITS["g_entropy_today"]*UNITS["T_today_GeV"]**3
T_STAR_RD = (90*(PLANCK*MASS)**2/(np.pi**2*G_RHO))**.25
K_STAR_RD = MASS*MPC*(S_TODAY/(2*np.pi**2/45*G_S*T_STAR_RD**3))**(1/3)
K_HARD_MIN, K_HARD_MAX = .001, 64.


class Expansion:
    """Canonical Friedmann/Klein–Gordon integration in cosmic time."""

    def __init__(self, name, fine):
        self.name = name
        self.case = case = FROZEN["histories"][name]
        self.p = case["pre_radiation_H_power"]
        self.width = case["inflation_end_blend_halfwidth"]
        self.rate = case["reheat_rate_per_efold"]
        self.rtol = 8e-12 if fine else 3e-10
        b = np.sqrt(2/(3*case["alpha"]))
        scale = PLANCK*np.sqrt(case["potential_over_Mbar4"])

        def friedmann(phi, speed):
            energy = .5*speed**2+(-np.expm1(-b*phi))**2
            return np.sqrt(energy/3)

        def evolution(_time, state):
            phi, speed, _ = state
            slope = -2*b*np.exp(-b*phi)*np.expm1(-b*phi)
            return speed, -3*friedmann(phi, speed)*speed-slope, friedmann(phi, speed)

        def end_event(_time, state):
            phi, speed, _ = state
            return .5*(speed/friedmann(phi, speed))**2-1

        end_event.terminal = True
        end_event.direction = 1
        phi0 = np.log1p(240*b*b)/b
        speed0 = -2*b*np.exp(-b*phi0)/np.sqrt(3)
        inflation = solve_ivp(evolution, (0., 500.), (phi0, speed0, 0.),
                              events=end_event, dense_output=True, rtol=self.rtol,
                              atol=self.rtol/25, method="RK45", max_step=.25)
        if len(inflation.t_events[0]) != 1:
            raise RuntimeError(f"{name} did not reach the end of inflation")
        stop = inflation.t_events[0][0]
        last = inflation.y_events[0][0]
        self.h_end = float(scale*friedmann(last[0], last[1]))
        times = np.linspace(0., stop, 30001 if fine else 16001)
        solution = inflation.sol(times)
        if self.width:
            def blend_end(_time, state):
                return state[2]-last[2]-self.width

            blend_end.terminal = True
            post = solve_ivp(evolution, (stop, stop+4.), last, events=blend_end,
                             dense_output=True, rtol=self.rtol, atol=self.rtol/25,
                             max_step=.005)
            if len(post.t_events[0]) != 1:
                raise RuntimeError("The continuation did not cover the blend")
            post_times = np.linspace(stop, post.t_events[0][0], 501)[1:]
            solution = np.column_stack((solution, post.sol(post_times)))
        phi, speed, efolds = solution
        hubble = friedmann(phi, speed)
        epsilon = .5*(speed/hubble)**2
        self.inflation = CubicHermiteSpline(efolds-last[2],
                                            np.log(scale*hubble/self.h_end), -epsilon)
        self.inflation_slope = self.inflation.derivative()

        t_reh = case["nominal_reheat_temperature_GeV"]
        if t_reh is None:
            self.h_reh = self.h_end
            t_reh = (90*(PLANCK*self.h_end)**2/(np.pi**2*G_RHO))**.25
        else:
            self.h_reh = np.sqrt(np.pi**2*G_RHO/90)*t_reh**2/PLANCK
        entropy_density = 2*np.pi**2/45*G_S*t_reh**3
        a_reh = (S_TODAY/entropy_density)**(1/3)
        self.n_reh = np.log(self.h_end/self.h_reh)/self.p
        self.a_end = a_reh*np.exp(-self.n_reh)
        self.k_end = self.a_end*self.h_end*MPC

        if self.h_reh >= MASS:
            self.a_star = a_reh*np.sqrt(self.h_reh/MASS)
        else:
            # Use the archived convention for a_star, from the pure stiff asymptote.
            self.a_star = self.a_end*(self.h_end/MASS)**(1/self.p)
        self.n_star = np.log(self.a_star/self.a_end)
        self.k_star = self.a_star*MASS*MPC
        self.s_star = S_TODAY/self.a_star**3

        if self.width:
            def conservation(n, _y):
                f = expit(self.rate*(n-self.n_reh))
                return [-self.p*(1-f)-2*f]

            end = self.n_reh+24.
            ns = np.linspace(self.width, end, 6501 if fine else 3501)
            solved = solve_ivp(conservation, (self.width, end), [-self.p*self.width],
                               t_eval=ns, rtol=self.rtol, atol=self.rtol/50,
                               max_step=.04)
            if not solved.success:
                raise RuntimeError(solved.message)
            self.conversion = CubicHermiteSpline(ns, solved.y[0], conservation(ns, None)[0])
            self.physical_end = brentq(lambda n: self.state(n)[1]-1,
                                       -self.width, self.width)
        else:
            self.physical_end = 0.

    def state(self, n):
        """Return ln(H/H_end) and -d ln(H)/dN."""
        if n < 0 and (not self.width or n <= -self.width):
            return float(self.inflation(n)), -float(self.inflation_slope(n))
        if self.width and n < self.width:
            z = (n+self.width)/(2*self.width)
            blend = z**5*(126+z*(-420+z*(540+z*(-315+70*z))))
            derivative = 630*z**4*(1-z)**4/(2*self.width)
            first = float(self.inflation(n))
            second = -self.p*n
            return ((1-blend)*first+blend*second,
                    -(1-blend)*float(self.inflation_slope(n))+blend*self.p
                    -derivative*(second-first))
        if self.width:
            fraction = expit(self.rate*(n-self.n_reh))
            return float(self.conversion(n)), self.p*(1-fraction)+2*fraction
        if n < self.n_reh:
            return -self.p*n, self.p
        return -self.p*self.n_reh-2*(n-self.n_reh), 2.

    def horizon_exit(self, k):
        return brentq(lambda n: np.log(self.k_end/k)+n+self.state(n)[0],
                      -100., min(-.05, self.physical_end))


def canonical_wave(bg, K, start, stop, origin):
    """Solve the canonically normalized massless wave in linear scale factor."""
    a0, a1 = np.exp(np.array([start, stop])-origin)

    def wkb_ratio(n):
        h, eps = bg.state(n)
        return np.sqrt(1-(2-eps)*(np.exp(n+h)/K)**2)

    ratio = wkb_ratio(start)
    phase_rate = K*np.exp(-start-bg.state(start)[0])/a0
    dn = .0003
    dlogw = np.log(wkb_ratio(start+dn)/wkb_ratio(start-dn))/(2*dn*a0)
    norm = 1/np.sqrt(ratio)

    def derivative(a, state):
        n = np.log(a)+origin
        h, eps = bg.state(n)
        horizon_ratio = K*np.exp(-n-h)
        u, du = state
        return du, -(2-eps)*du/a-(horizon_ratio**2-2+eps)*u/a**2

    run = solve_ivp(derivative, (a0, a1),
                    (complex(norm), norm*(-.5*dlogw-1j*ratio*phase_rate)),
                    method="DOP853", rtol=bg.rtol, atol=bg.rtol/40)
    if not run.success:
        raise RuntimeError(run.message)
    return run.y[:, -1]


def tensor_power(bg, k, fine):
    exit_at = bg.horizon_exit(k)
    K = k/bg.k_end
    inside = 112. if fine else 56.
    start = brentq(lambda n: n+bg.state(n)[0]-np.log(K/inside),
                   exit_at-8., exit_at)
    stop = min(exit_at+(9. if fine else 7.), bg.physical_end-.12)
    if stop-exit_at < 6:
        raise ValueError("Requested primordial mode did not freeze")
    amplitude, _ = canonical_wave(bg, K, start, stop, exit_at)
    h_exit = bg.h_end*np.exp(bg.state(exit_at)[0])
    return float(2/np.pi**2*(h_exit/PLANCK)**2*abs(amplitude*np.exp(exit_at-stop))**2)


def primordial(bg, fine, qmin, qmax, additional_range=None):
    lo, hi = bg.k_star*qmin, bg.k_star*qmax
    if additional_range:
        lo, hi = min(lo, additional_range[0]), max(hi, additional_range[1])
    nodes = np.linspace(np.log(lo), np.log(hi), 65 if fine else 33)
    power = np.array([tensor_power(bg, np.exp(k), fine) for k in nodes])
    interpolated = PchipInterpolator(nodes, np.log(power), extrapolate=False)

    def spectrum(k):
        return np.exp(interpolated(np.clip(np.log(k), nodes[0], nodes[-1])))

    return spectrum


def longitudinal(bg, fine):
    """Second-order longitudinal equation and time-averaged stress readout."""
    q_min, q_max = (.0005, 128.) if fine else (.001, 64.)
    q = np.geomspace(q_min, q_max, 643 if fine else 323)
    q2 = q*q
    x0 = 3e-8 if fine else 1e-7
    sample = np.unique(np.r_[np.geomspace(x0, 19.9, 151),
                            np.linspace(20., 24., 257), np.linspace(28., 32., 257)])

    def elements(x):
        n = bg.n_star+np.log(x)
        log_h, eps = bg.state(n)
        ahat = x*x*bg.h_end/MASS*np.exp(log_h)
        return x*x/(q2+x*x), ahat, eps

    def equation(x, state):
        field, slope = state.reshape(2, -1)
        kinetic, ahat, eps = elements(x)
        curvature = -(4-eps-2*kinetic)*slope/x-(x*x+q2)*field/ahat**2
        return np.r_[slope, curvature]

    # The second-order equation needs A and dA/dx, both continuous at a
    # change in epsilon with continuous H. The early-matter step is explicit.
    segments = [x0, 32.]
    boundary = np.exp(bg.n_reh-bg.n_star)
    if not bg.width and x0 < boundary < 32.:
        segments.insert(1, boundary)
    state = np.r_[np.ones(len(q)), np.zeros(len(q))]
    states = []
    for start, stop in zip(segments[:-1], segments[1:]):
        run = solve_ivp(equation, (start, stop), state, method="DOP853",
                        dense_output=True, rtol=bg.rtol, atol=bg.rtol/40)
        if not run.success:
            raise RuntimeError(run.message)
        mask = (sample >= start) & ((sample < stop) | (stop == 32.))
        states.extend(run.sol(sample[mask]).T)
        state = run.y[:, -1]

    rows = []
    for x, state in zip(sample, states):
        field, slope = state.reshape(2, -1)
        kinetic, ahat, eps = elements(x)
        dfield = ahat*slope
        electric = kinetic*dfield*dfield
        mass_energy = x*x*field*field
        energy = electric+mass_energy
        pressure = ((3-2*kinetic)*electric-mass_energy)/3
        curve = equation(x, state)[len(q):]
        ddfield = ahat*(curve+(2-eps)*slope/x)
        denergy = (2*kinetic*(1-kinetic)*dfield**2/x+2*kinetic*dfield*ddfield
                   +2*x*field*field+2*x*x*field*slope)
        conservation = max(abs(x*denergy-energy+3*pressure)/energy)
        weight = q2*energy/np.sqrt(x*x+q2)
        rows.append((weight, q2*energy/ahat**2, q2*pressure/ahat**2, conservation))

    weight, stress_energy, stress_pressure = (np.array([row[i] for row in rows])
                                               for i in range(3))

    def readout(low, high):
        band = (sample >= low) & (sample <= high)
        return simpson(weight[band], x=sample[band], axis=0)/(high-low)

    return (q, readout(20., 24.), readout(28., 32.), sample, stress_energy,
            stress_pressure, max(row[3] for row in rows))


def high_q_tail(bg):
    """Large-q frozen mode matching to the massless regular Bessel solution."""
    p = bg.p if bg.h_reh < MASS else 2.
    index = 1/(p-1)-.5
    scale = (2*gamma(1+index))**2*2**(2*index-1)*(p-1)**(2*index+1)/np.pi
    return 2*index, scale


def fold(bg, q, response, spectrum):
    power, amplitude = high_q_tail(bg)
    tail = np.geomspace(q[-1], bg.case["q_tail_stop"], 121)
    soft = response*spectrum(q*bg.k_star)
    tail_weight = amplitude*tail**(-power)*spectrum(tail*bg.k_star)
    band = simpson(soft, x=np.log(q))+.5*soft[0]
    high = simpson(tail_weight, x=np.log(tail))
    normalization = MASS**2*PLANCK**2/(16*bg.s_star*TARGET)
    leftover = amplitude*tail[-1]**(-power)*spectrum(tail[-1]*bg.k_star)/power
    return {"ratio": float(normalization*(band+high)),
            "analytic_tail_share": float(high/(band+high)),
            "post_cut_upper_if_nonincreasing": float(normalization*leftover),
            "UV_coefficient_ratio_at_direct_edge": float(response[-1]/
                                                         (amplitude*q[-1]**(-power)))}


def hard_occupation(bg, fine):
    """Compute hard transfer by matching real Bessel coefficients."""
    k = np.geomspace(K_HARD_MIN, K_HARD_MAX, 233 if fine else 117)
    occ = []
    norm_errors = []
    stop = bg.width+.055
    yf = np.exp(stop)
    for K in k:
        start = min(-bg.width-.07, brentq(lambda n: n+bg.state(n)[0]-
                                         np.log(K/(104 if fine else 52)), -24., -.16)
                    if K < .7*(104 if fine else 52) else -bg.width-.07)
        value, derivative = canonical_wave(bg, K, start, stop, 0.)
        p = bg.p
        s = K*np.exp((p-1)*stop)/(p-1)
        order = 1/(p-1)-.5
        prefactor = np.sqrt(np.pi*s/2)
        bessel = np.array([jv(order, s), yv(order, s)])
        bessel_slope = np.array([jvp(order, s), yvp(order, s)])
        basis = prefactor*np.vstack((bessel, (p-1)/yf*(.5*bessel+s*bessel_slope)))
        j_coeff, y_coeff = np.linalg.solve(basis, [value, derivative])
        alpha, beta = .5*(j_coeff+1j*y_coeff), .5*(j_coeff-1j*y_coeff)
        occ.append(float(abs(beta)**2))
        norm_errors.append(float(abs(abs(alpha)**2-abs(beta)**2-1)/max(1., abs(beta)**2)))
    return k, np.array(occ), max(norm_errors)


def measure(name, fine):
    bg = Expansion(name, fine)
    q, early, late, xs, energy, pressure, ward = longitudinal(bg, fine)
    if bg.width:
        shared_range = (q[0]*K_STAR_RD, bg.case["q_tail_stop"]*bg.k_star)
        fitted_power = primordial(bg, fine, q[0], bg.case["q_tail_stop"], shared_range)
        matched_power = fitted_power
    else:
        fitted_power = primordial(bg, fine, q[0], bg.case["q_tail_stop"])
        matched_bg = Expansion("narrow", fine)
        q_stop = matched_bg.case["q_tail_stop"]
        matched_power = primordial(matched_bg, fine, q[0], q_stop,
                                  (q[0]*K_STAR_RD, q_stop*matched_bg.k_star))
    fitted = fold(bg, q, late, fitted_power)
    previous = fold(bg, q, early, fitted_power)
    matched = fold(bg, q, late, matched_power)
    flat = fold(bg, q, late, lambda k: 0*k+UNITS["A_s_pivot"])
    pivot_tensors = tensor_power(bg, PIVOT, fine)
    weighted_energy = simpson(energy*fitted_power(q*bg.k_star), x=np.log(q), axis=1)/48
    weighted_pressure = simpson(pressure*fitted_power(q*bg.k_star), x=np.log(q), axis=1)/48
    result = dict(case=name, fine=fine, H_end_GeV=bg.h_end, a_end=bg.a_end,
                  entropy_per_a_star_cubed_GeV3=bg.s_star, k_star_Mpc_inverse=bg.k_star,
                  pivot_efolds=bg.physical_end-bg.horizon_exit(PIVOT),
                  r_at_frozen_As=pivot_tensors/UNITS["A_s_pivot"],
                  fitted=fitted, matched_stiff_power=matched, common_flat_power=flat,
                  late_window_shift=abs(previous["ratio"]/fitted["ratio"]-1),
                  longitudinal_conservation_relative_residual=ward,
                  max_resolved_soft_energy_fraction=float(max(weighted_energy)),
                  max_resolved_soft_pressure_abs_fraction=float(max(abs(weighted_pressure))))
    output = ROOT / "measured"
    output.mkdir(exist_ok=True)
    label = f"{name}_{'fine' if fine else 'coarse'}"
    arrays = dict(q=q, stress_x=xs, stress_energy=weighted_energy, stress_pressure=weighted_pressure,
                  soft_response=late, frozen_spectrum=fitted_power(q*bg.k_star))
    if bg.width:
        k, occ, norm_error = hard_occupation(bg, fine)
        hard = simpson(k**4*occ, x=np.log(k))
        result.update(hard_K4_occupation_integral=float(hard),
                      hard_max_relative_wronskian_error=norm_error)
        arrays.update(hard_K=k, hard_occupation=occ)
    np.savez_compressed(output / (label+".npz"), **arrays)
    (output / (label+".json")).write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=FROZEN["histories"])
    parser.add_argument("--fine", action="store_true")
    options = parser.parse_args()
    measure(options.case, options.fine)
