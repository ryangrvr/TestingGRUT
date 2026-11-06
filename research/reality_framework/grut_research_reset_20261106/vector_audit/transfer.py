"""Independent longitudinal transfer in cosmic log scale factor."""
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.special import gamma

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
M = 1e-12  # 1 meV in GeV
MPL = 2.435e18
G = 106.75
G_S0 = 3.91
T_CMB = 2.34865e-13
GEV_TO_INVERSE_MPC = 1.563738305e38
A_S = 2.1e-9
DM_PER_ENTROPY = .44e-9
K_PIVOT = .05
T_STAR = np.sqrt(M*MPL/(np.pi*np.sqrt(G/90)))
IR_RADIATION = 2*gamma(.75)**2/np.pi


@dataclass(frozen=True)
class History:
    name: str
    w: float = 1/3
    reheat_temperature: float | None = None

    @property
    def p(self):
        return 1.5*(1+self.w)

    @property
    def h_reheat(self):
        if self.reheat_temperature is None:
            return None
        return (self.reheat_temperature/T_STAR)**2

    @property
    def x_reheat(self):
        h = self.h_reheat
        if h is None:
            return None
        return h**(-.5 if h >= 1 else -1/self.p)

    @property
    def entropy_reference(self):
        if self.reheat_temperature is None:
            return entropy(T_STAR)
        return self.x_reheat**3*entropy(self.reheat_temperature)

    @property
    def k_star(self):
        return M*GEV_TO_INVERSE_MPC*(entropy(T_CMB, G_S0)/self.entropy_reference)**(1/3)

    @property
    def uv_power(self):
        p = self.p if self.h_reheat is not None and self.h_reheat < 1 else 2.
        nu = (3-p)/(2*(p-1))
        amplitude = 2**(2*nu+1)*gamma(nu+1)**2/np.pi*(p-1)**(2*nu+1)
        return 2*nu, amplitude


HISTORIES = {
    "prompt": History("prompt"),
    "early_matter": History("early_matter", 0., 1e9),
    "late_stiff": History("late_stiff", .6, 200.),
}


def entropy(temperature, degrees=G):
    return 2*np.pi**2/45*degrees*temperature**3


def write(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")
    temporary.replace(path)


def solve(q, history, start, stops, rtol):
    """A and P=x K h dA/dln x, with h=H/m and K=x²/(x²+q²)."""
    q2 = q*q
    transition = history.x_reheat
    p = history.p if transition is not None else 2.
    if transition is not None:
        start = min(start, transition/100)
        pre_scale = history.h_reheat*transition**p
        rad_scale = history.h_reheat*transition**2
    else:
        pre_scale = rad_scale = 1.
    h_start = pre_scale*start**(-p)
    small = q2/(start*h_start)**2
    coefficient = 1/(2*(p-1)*(p+1))
    kinetic = start**2/(start**2+q2)
    initial = np.r_[1-coefficient*small,
                    -start*kinetic*h_start*small/(p+1)]

    def segment(left, right, output, state, index, scale):
        def rhs(n, data):
            x = np.exp(n)
            h = scale*x**(-index)
            a, momentum = data.reshape(2, -1)
            kinetic = x*x/(x*x+q2)
            return np.r_[momentum/(x*kinetic*h), -x*a/h]

        sol = solve_ivp(rhs, (np.log(left), np.log(right)), state,
                        t_eval=np.log(output), method="DOP853",
                        rtol=rtol, atol=rtol/50)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y.T

    if transition is not None:
        initial = segment(start, transition, [transition], initial, p, pre_scale)[-1]
        start = transition
    states = segment(start, max(stops), stops, initial, 2., rad_scale)
    twice_action = []
    for x, state in zip(stops, states):
        a, momentum = state.reshape(2, -1)
        b = x/np.sqrt(x*x+q2)
        c = q2/(x*x+q2)
        ah = rad_scale/x
        canonical = b*a
        derivative = momentum/b+ah*c*canonical
        frequency = np.sqrt(x*x+q2+3*ah**2*c*(1-c))
        twice_action.append((derivative**2+frequency**2*canonical**2)/frequency)
    return np.array(twice_action)


def evaluate(history, q_min, q_max, count, start, stops, rtol):
    q = np.geomspace(q_min, q_max, count)
    action = solve(q, history, start, stops, rtol)
    spectrum = action*q*q
    power, amplitude = history.uv_power
    uv_tail = amplitude*q_max**(-power)/power
    integrals = (simpson(spectrum, x=np.log(q), axis=1)
                 +.5*action[:, 0]*q_min**2+uv_tail)
    band = (q >= .25) & (q <= 4)
    row = dict(q_min=q_min, q_max=q_max, points=count, late_readouts=stops,
               integral=integrals.tolist(), analytic_UV_tail=uv_tail,
               relative_UV_coefficient_error=float(
                   abs(spectrum[-1, -1]/(amplitude*q_max**(-power))-1)),
               infrared_coefficient=float(action[-1, 0]),
               finite_band_integral=float(simpson(spectrum[-1, band], x=np.log(q[band]))))
    return row, q, spectrum[-1]


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    previous = json.loads((HERE.parent / "bounded_route/results/vector_transfer.json").read_text())
    rows = {}
    for name, history in HISTORIES.items():
        coarse, _, _ = evaluate(history, .001, 64., 321, 1e-5, (12., 18.), 4e-9)
        fine, q, w = evaluate(history, .0005, 128., 641, 5e-6, (16., 24.), 3e-10)
        value = fine["integral"][-1]
        r_flat = 16*DM_PER_ENTROPY*history.entropy_reference/(M*M*MPL*MPL*A_S*value)
        flat_band_upper = r_flat*value/fine["finite_band_integral"]
        row = dict(coarse=coarse, fine=fine, entropy_gev3=history.entropy_reference,
                   k_star_Mpc_inverse=history.k_star, H_reheat_over_m=history.h_reheat,
                   x_reheat=history.x_reheat, flat_r=r_flat,
                   flat_r_using_only_finite_band=flat_band_upper,
                   relative_grid_shift=abs(coarse["integral"][-1]/value-1),
                   relative_late_shift=abs(fine["integral"][0]/value-1))
        rows[name] = row
        np.savez_compressed(RESULTS / (name+"_transfer.npz"), q=q, weight=w)
    reference = previous["fine"]["transfer_integrals"][-1]
    p = rows["prompt"]
    error = abs(p["fine"]["integral"][-1]/reference-1)
    early_error = abs(rows["early_matter"]["fine"]["integral"][-1]/
                      p["fine"]["integral"][-1]-1)
    prompt_ir_error = abs(p["fine"]["infrared_coefficient"]/IR_RADIATION-1)
    checks = dict(independent_prompt=error < .005,
                  infrared=prompt_ir_error < .015,
                  ultraviolet=max(row["fine"]["relative_UV_coefficient_error"]
                                  for row in rows.values()) < .015,
                  early_transfer=early_error < .005,
                  quadrature=max(row["relative_grid_shift"] for row in rows.values()) < .005,
                  late=max(row["relative_late_shift"] for row in rows.values()) < .005,
                  late_counterexample_band=rows["late_stiff"]["flat_r_using_only_finite_band"] <
                      rows["prompt"]["flat_r"])
    result = dict(passed=all(checks.values()), checks={k: bool(v) for k, v in checks.items()},
                  previous_prompt_integral=reference, relative_prompt_error=error,
                  early_relative_transfer_error=early_error, prompt_ir_error=prompt_ir_error,
                  radiation_T_star_gev=T_STAR, histories=rows)
    write(RESULTS / "transfer.json", result)
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit("The independent transfer qualification missed.")


if __name__ == "__main__":
    main()
