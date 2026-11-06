"""Background, vacuum modes, and vector abundance for the frozen histories."""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import json

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import gamma

from transfer import (A_S, DM_PER_ENTROPY, G, G_S0, GEV_TO_INVERSE_MPC, HISTORIES,
                      K_PIVOT, M, MPL, RESULTS, T_CMB, write)

LOG_PIVOT = np.log(K_PIVOT)


@dataclass(frozen=True)
class Precision:
    inside: float
    after: float
    mode_tolerance: float
    background_tolerance: float
    tensor_samples: int
    tilt_step: float


COARSE = Precision(72., 8., 2e-9, 2e-10, 17, .08)
FINE = Precision(144., 10., 2e-11, 2e-12, 33, .04)


class Background:
    def __init__(self, alpha, tolerance):
        b = np.sqrt(2/(3*alpha))

        def derivative(_n, y):
            phi, speed = y
            eps = .5*speed**2
            return np.array([speed, -(3-eps)*(speed+2*b/np.expm1(b*phi))])

        def end(_n, y):
            return .5*y[1]**2-1

        end.terminal = True
        end.direction = 1
        slow_end = np.log1p(np.sqrt(2.)*b)/b

        def slow_duration(phi):
            return ((np.exp(b*phi)-np.exp(b*slow_end))/(2*b*b)
                    -(phi-slow_end)/(2*b))

        phi_start = brentq(lambda phi: slow_duration(phi)-100, slow_end, slow_end+50)
        speed_start = -2*b/np.expm1(b*phi_start)
        sol = solve_ivp(derivative, (0., 125.), [phi_start, speed_start], events=end,
                        dense_output=True, method="DOP853", rtol=tolerance,
                        atol=tolerance/30, max_step=.2)
        if not sol.success or not len(sol.t_events[0]):
            raise RuntimeError("Inflation did not reach epsilon_H=1.")
        self.duration = float(sol.t_events[0][0])
        times = np.linspace(0, self.duration, 1+int(self.duration/.01))
        phi, speed = sol.sol(times)
        eps = .5*speed**2
        eps2 = 2*derivative(times, (phi, speed))[1]/speed
        logh = np.log(-np.expm1(-b*phi))-.5*np.log(3-eps)
        self.track = CubicSpline(times-self.duration, np.column_stack((logh, eps, eps2)))

    def exit(self, log_k_internal):
        return brentq(lambda n: n+self.track(n)[0]-log_k_internal,
                      -self.duration+1, -.02, xtol=4e-13)


def mode(track, exit_n, scalar, precision):
    log_h, eps_exit, _ = track(exit_n)
    log_k = exit_n+log_h

    def log_ratio(n):
        return log_k-n-track(n)[0]

    start = brentq(lambda n: log_ratio(n)-np.log(precision.inside),
                   exit_n-10, exit_n)
    stop = min(exit_n+precision.after, -.05)
    if stop-exit_n < 6:
        raise ValueError("A mode reached the end of inflation before the frozen readout.")

    def mass(n):
        _, eps, eps2 = track(n)
        z_slope = 1+eps2/2 if scalar else 1.
        z_curve = track(n, 1)[2]/2 if scalar and isinstance(track, CubicSpline) else 0.
        return z_curve+z_slope**2+(1-eps)*z_slope

    _, eps_i, eps2_i = track(start)
    q_i = precision.inside
    mu = mass(start)
    mu_prime = (mass(start+.001)-mass(start-.001))/.002
    w_ratio_sq = 1-mu/q_i**2
    w_ratio_derivative = -(mu_prime+2*(1-eps_i)*mu)/q_i**2
    amplitude = np.exp(exit_n-start)/w_ratio_sq**.25
    z_slope = 1.
    if scalar:
        amplitude *= np.sqrt(eps_exit/eps_i)
        z_slope += eps2_i/2
    derivative = amplitude*(-1j*q_i*np.sqrt(w_ratio_sq)
                            -w_ratio_derivative/(4*w_ratio_sq)-z_slope)

    def rhs(n, state):
        lh, eps, eps2 = track(n)
        value, speed = state
        friction = 3-eps+(eps2 if scalar else 0.)
        return (speed, -friction*speed-np.exp(2*(log_k-n-lh))*value)

    sol = solve_ivp(rhs, (start, stop), [complex(amplitude), derivative], t_eval=[stop],
                    method="DOP853", rtol=precision.mode_tolerance,
                    atol=precision.mode_tolerance/30)
    if not sol.success:
        raise RuntimeError(sol.message)
    return float(abs(sol.y[0, -1])**2)


def mode_controls():
    rows = []
    for eps in (.0, .02):
        def track(n):
            return -eps*n, eps, 0.
        nu = (3-eps)/(2*(1-eps))
        exact = (2**(nu-1.5)*gamma(nu)/gamma(1.5))**2*(1-eps)**(2*nu-1)
        computed = mode(track, -30., False, FINE)
        scalar = mode(track, -30., True, FINE) if eps else None
        rows.append(dict(epsilon=eps, expected=float(exact), measured=computed,
                         relative_error=float(abs(computed/exact-1)),
                         scalar_relative_error=abs(scalar/exact-1) if eps else None))
    return rows


class Calibrated:
    def __init__(self, alpha, history, precision):
        self.alpha, self.history, self.precision = alpha, history, precision
        self.background = Background(alpha, precision.background_tolerance)
        self.log_v = brentq(self.amplitude_residual, -35., -16., xtol=1e-10)
        self.log_offset, self.end_a, self.reheat_temperature, self.H_end = self.mapping(self.log_v)
        self.pivot_exit = self.background.exit(LOG_PIVOT-self.log_offset)
        self.scalar_pivot = np.exp(self.scalar_log_power(LOG_PIVOT, self.log_v,
                                                         self.log_offset))

    def mapping(self, log_v):
        h_end = self.background.track(0.)[0]
        log_h_scale = np.log(MPL)+log_v/2
        h_end = np.exp(log_h_scale+h_end)
        rho_end = 3*MPL*MPL*h_end*h_end
        instantaneous = (30*rho_end/(np.pi**2*G))**.25
        reheat_temperature = self.history.reheat_temperature or instantaneous
        if reheat_temperature > instantaneous*(1+1e-13):
            raise ValueError("The specified radiation temperature exceeds the end energy.")
        a_reh = T_CMB/reheat_temperature*(G_S0/G)**(1/3)
        a_end = a_reh*(reheat_temperature/instantaneous)**(4/(3*(1+self.history.w)))
        return (np.log(a_end*GEV_TO_INVERSE_MPC)+log_h_scale,
                float(a_end), float(reheat_temperature), float(h_end))

    def scalar_log_power(self, log_k, log_v, offset):
        n = self.background.exit(log_k-offset)
        lh, eps, _ = self.background.track(n)
        response = mode(self.background.track, n, True, self.precision)
        return log_v+2*lh+np.log(response/(8*np.pi**2*eps))

    def amplitude_residual(self, log_v):
        offset = self.mapping(log_v)[0]
        return self.scalar_log_power(LOG_PIVOT, log_v, offset)-np.log(A_S)

    def tensor(self, log_k):
        n = self.background.exit(log_k-self.log_offset)
        log_h = self.background.track(n)[0]
        response = mode(self.background.track, n, False, self.precision)
        return np.log(2/np.pi**2)+self.log_v+2*log_h+np.log(response)

    def abundance(self):
        history, precision = self.history, self.precision
        with np.load(RESULTS / (history.name+"_transfer.npz")) as packet:
            q, weight = packet["q"], packet["weight"]
        power, amplitude = history.uv_power
        log_k_star = np.log(history.k_star)
        # Read every inflation mode at least eight e-folds after exit.
        # Bound the positive, decreasing-power contribution beyond the last mode.
        end_limit = -8.2+self.background.track(-8.2)[0]+self.log_offset-log_k_star
        q_stop = min(np.exp(end_limit), 1e8 if power < 1 else 1e6)
        if q_stop <= q[-1]:
            raise ValueError("The directly integrated transfer range includes late exit.")
        tail_q = np.geomspace(q[-1], q_stop, 97)[1:]
        log_q = np.log(np.r_[q, tail_q])
        w = np.r_[weight, amplitude*tail_q**(-power)]
        nodes = np.linspace(log_q[0], log_q[-1], precision.tensor_samples)
        log_power = np.array([self.tensor(log_k_star+node) for node in nodes])
        tensor = np.exp(CubicSpline(nodes, log_power)(log_q))
        integral = (simpson(w*tensor, x=log_q)+.5*w[0]*tensor[0])
        window = (np.exp(log_q) >= .25) & (np.exp(log_q) <= 4.)
        finite = simpson(w[window]*tensor[window], x=log_q[window])
        prefactor = M*M*MPL*MPL/(16*history.entropy_reference*DM_PER_ENTROPY)
        return dict(density_ratio=float(prefactor*integral),
                    finite_band_density_ratio=float(prefactor*finite),
                    tail_density_ratio_upper=float(prefactor*tensor[-1]*
                                                   amplitude*q_stop**(-power)/power),
                    tail_q_stop=float(q_stop),
                    nonincreasing_sampled_tensor=bool(np.max(np.diff(log_power)) <= 0),
                    log_k_star=log_k_star,
                    vector_peak_exit=-float(self.background.exit(
                        np.log(history.k_star*1.375)-self.log_offset)))

    def summary(self, abundance, full):
        tp = self.tensor(LOG_PIVOT)
        result = dict(alpha=self.alpha, V0_over_Mpl4=float(np.exp(self.log_v)),
                      scalar_power_pivot=float(self.scalar_pivot),
                      H_pivot_gev=float(MPL*np.exp(self.log_v/2+
                                                 self.background.track(self.pivot_exit)[0])),
                      H_end_gev=self.H_end, pivot_efolds=-float(self.pivot_exit),
                      tensor_ratio=float(np.exp(tp)/self.scalar_pivot),
                      reheat_temperature_gev=self.reheat_temperature,
                      a_end=self.end_a, **abundance)
        if full:
            d = self.precision.tilt_step
            y = np.array([self.scalar_log_power(LOG_PIVOT+i*d, self.log_v, self.log_offset)
                          for i in (-2, -1, 0, 1, 2)])
            slope = (y[0]-8*y[1]+8*y[3]-y[4])/(12*d)
            curvature = (-y[0]+16*y[1]-30*y[2]+16*y[3]-y[4])/(12*d*d)
            result.update(n_s=float(1+slope), running=float(curvature))
        return result


def calibrate(history, precision, band_only=False):
    cache = {}
    key = "finite_band_density_ratio" if band_only else "density_ratio"

    def get(alpha):
        if alpha not in cache:
            model = Calibrated(alpha, history, precision)
            cache[alpha] = model, model.abundance()
        return cache[alpha]

    def residual(alpha):
        return get(alpha)[1][key]-1

    edge = [residual(1.)+1, residual(40.)+1]
    if (edge[0]-1)*(edge[1]-1) > 0:
        return dict(root_bracketed=False, calibration_key=key, bracket_densities=edge)
    alpha = brentq(residual, 1., 40., xtol=2e-5)
    model, abundance = get(alpha)
    return dict(root_bracketed=True, calibration_key=key, bracket_densities=edge,
                **model.summary(abundance, full=True))


def fitted(history):
    coarse, fine = calibrate(history, COARSE), calibrate(history, FINE)
    row = dict(coarse=coarse, fine=fine)
    if fine["root_bracketed"] and coarse["root_bracketed"]:
        row["r_refinement"] = abs(fine["tensor_ratio"]/coarse["tensor_ratio"]-1)
        row["tilt_refinement"] = abs(fine["n_s"]-coarse["n_s"])
        row["running_refinement"] = abs(fine["running"]-coarse["running"])
    print(history.name, fine, flush=True)
    return row


def fixed_shape(alpha):
    model = Calibrated(alpha, HISTORIES["prompt"], FINE)
    return model.summary(model.abundance(), full=True)


def main():
    qualification = json.loads((RESULTS / "transfer.json").read_text())
    if not qualification["passed"]:
        raise SystemExit("Qualify the independent transfer first.")
    controls = mode_controls()
    if max(max(row["relative_error"], row["scalar_relative_error"] or 0)
           for row in controls) > 2e-5:
        raise SystemExit("The constant-epsilon mode control missed.")
    print("Vacuum mode controls", controls, flush=True)
    with ProcessPoolExecutor(max_workers=3) as executor:
        fitted_jobs = {name: executor.submit(fitted, history)
                       for name, history in HISTORIES.items()}
        fixed_jobs = {str(alpha): executor.submit(fixed_shape, alpha) for alpha in (1., 4.)}
        finite_job = executor.submit(calibrate, HISTORIES["late_stiff"], FINE, True)
        results = {name: job.result() for name, job in fitted_jobs.items()}
        fixed = {name: job.result() for name, job in fixed_jobs.items()}
        finite = finite_job.result()
    convergence = all(row.get("r_refinement", 1.) < .002
                      and row.get("tilt_refinement", 1.) < 2e-4
                      and abs(row["fine"]["density_ratio"]-1) < .001
                      for row in results.values())
    result = dict(fixed_mass_eV=M*1e9, scalar_power_input=A_S,
                  density_per_entropy_input_eV=DM_PER_ENTROPY*1e9,
                  controls=controls, fitted_histories=results,
                  fixed_shape_controls=fixed, late_stiff_finite_band_fit=finite,
                  fitted_refinements_pass=bool(convergence))
    write(RESULTS / "inflation.json", result)


if __name__ == "__main__":
    main()
