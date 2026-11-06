"""Conservative relic bound for the frozen 40 GeV Higgs-portal comparator."""
import argparse
import hashlib
import json
from math import exp, log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import kv, kve

from match import DIRECTORY
from model import write_json

M = 40.
HIGGS = 125.09
VEV = 246.22
SM_WIDTH = .00410
MIN_POLE_WIDTH = .0040
WIDTH_ENVELOPE = .020
ATLAS_LIMIT = .107
T0 = 2.
X0 = M/T0
GPLANCK = 1.22089e19
G_MAX, G_MIN = 106.75, 3.90
S_OVER_H_GMAX = 12.
DM_YIELD = .44e-9/M
WIDTH_SOURCE = ("https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/"
                "9da32dcda3b65d4763b4af191d553018247adcfa/data/lhc-hxswg/sm/br/BR.txt")


def sigma_per_lambda_squared(energy, width, pole_width):
    denominator = (energy**2-HIGGS**2)**2+HIGGS**2*pole_width**2
    return 2*VEV**2*width/(energy*denominator)


def thermal_upper_weight(energy):
    """Uniform over T <= T0 for energies >= 95 GeV."""
    y = energy/T0
    bessel = sqrt(pi/(2*y))*(1+3/(8*y))/(pi/(2*X0))
    return energy**3*sqrt(energy**2-4*M*M)/(8*M**4*T0)*bessel*exp(-(energy-2*M)/T0)


def exact_thermal_weight(energy, x):
    temperature = M/x
    return (energy**3*sqrt(energy**2-4*M*M)/(8*M**4*temperature)*
            kve(1, energy/temperature)/kve(2, x)**2*exp(-(energy-2*M)/temperature))


def integrals():
    low = 2*VEV**2*WIDTH_ENVELOPE/(2*M*(HIGGS**2-95**2)**2)
    breakpoints = (95., 110., 120., HIGGS-.05, HIGGS, HIGGS+.05, 130., 140.)
    middle = sum(quad(lambda e: thermal_upper_weight(e)*sigma_per_lambda_squared(
        e, WIDTH_ENVELOPE, MIN_POLE_WIDTH), a, b, epsabs=1e-15, epsrel=1e-9)[0]
        for a, b in zip(breakpoints, breakpoints[1:]))
    # A deliberately loose cross-section-times-velocity ceiling, 1 GeV^-2,
    # above 140 GeV. The Boltzmann weight renders its contribution small.
    tail = quad(thermal_upper_weight, 140, 300, epsabs=1e-19, epsrel=1e-8)[0]
    tail += quad(thermal_upper_weight, 300, np.inf, epsabs=1e-40)[0]
    return low, middle, tail


def relic_bound(table):
    mass, width = table[:, 0], table[:, -1]
    selected = (mass >= 2*M) & (mass <= 140)
    assert mass[0] <= 2*M and max(width[selected]) < WIDTH_ENVELOPE
    decay_factor = VEV**2/(32*pi*HIGGS)*sqrt(1-4*M*M/HIGGS**2)
    allowed = SM_WIDTH*ATLAS_LIMIT/(1-ATLAS_LIMIT)
    lambda_max = sqrt(allowed/decay_factor)
    low, mid, tail = integrals()
    rate = lambda_max**2*(low+mid)+tail
    entropy_factor = 1+log(G_MAX/G_MIN)/3
    # H = sqrt(4*pi^3/45) sqrt(g_rho) T^2 / M_Planck.
    integral_coefficient = ((2*pi*pi/45)/sqrt(4*pi**3/45)*
                            GPLANCK*S_OVER_H_GMAX*T0*entropy_factor)
    y_equilibrium_floor = 45*X0**2*kv(2, X0)/(4*pi**4*G_MAX)
    y_floor = 1/(1/y_equilibrium_floor+integral_coefficient*rate)
    needed = sqrt((1/DM_YIELD-1/y_equilibrium_floor-integral_coefficient*tail)/
                  (integral_coefficient*(low+mid)))
    required_width = decay_factor*needed**2
    controls = {}
    for x in (20., 30., 50.):
        value = quad(lambda e: exact_thermal_weight(e, x)*lambda_max**2*
                     sigma_per_lambda_squared(e, np.interp(e, mass, width), SM_WIDTH+allowed),
                     2*M, 140., points=[HIGGS-.05, HIGGS, HIGGS+.05],
                     epsabs=1e-17, epsrel=1e-6)[0]
        controls[str(x)] = value
    assert max(controls.values()) < rate
    return dict(mass_gev=M, higgs_mass_gev=HIGGS, sm_width_gev=SM_WIDTH,
                max_reference_width_80_to_140_gev=float(max(width[selected])),
                allowed_invisible_br_95=ATLAS_LIMIT, max_abs_lambda=lambda_max,
                low_and_mid_rate_per_lambda_squared=[low, mid],
                high_tail_rate_bound=tail, max_thermal_rate_gev_minus2=rate,
                thermal_rate_controls=controls, equilibrium_yield_floor_at_x20=y_equilibrium_floor,
                late_annihilation_integral_upper=integral_coefficient*rate,
                late_yield_floor=y_floor, all_dark_matter_yield=DM_YIELD,
                abundance_over_target_lower=y_floor/DM_YIELD,
                necessary_abs_lambda_lower=needed,
                necessary_invisible_br_lower=required_width/(SM_WIDTH+required_width),
                tested_point_excluded=bool(y_floor > DM_YIELD))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("width_table", type=Path)
    args = parser.parse_args()
    decision = json.loads((DIRECTORY / "route_decision.json").read_text())
    if decision["decision"] != "close_ALP_lane":
        raise SystemExit("The conditional replacement is still closed.")
    table = np.loadtxt(args.width_table, skiprows=1)
    result = relic_bound(table)
    result.update(width_table=WIDTH_SOURCE, width_sha256=hashlib.sha256(
        args.width_table.read_bytes()).hexdigest())
    write_json(DIRECTORY / "portal.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
