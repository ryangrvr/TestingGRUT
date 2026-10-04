#!/usr/bin/env python3
"""MODEL-ONLY plumbing test of Card01CAMB inside Cobaya.  NO data likelihood: the only "likelihood" is a
synthetic external function that requests H(z), angular-diameter distance, rdrag and CMB Cls and returns 0.
It checks that the sampled `card01_eps` reaches CAMB and changes the predictions exactly as direct CAMB does.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import camb
from cobaya.model import get_model
import card01_cobaya_theory as T

ZS = [0.51, 1.32, 2.33]
captured = {}


def dummy(_self=None):
    return 0.0


def make_model():
    def lik(_self):
        pr = _self.provider
        captured["H"] = pr.get_Hubble(ZS)
        captured["DA"] = pr.get_angular_diameter_distance(ZS)
        captured["rdrag"] = pr.get_param("rdrag")
        captured["tt220"] = pr.get_Cl(ell_factor=True, units="muK2")["tt"][220]
        return 0.0
    info = {
        "theory": {"card01": {"external": T.Card01CAMB, "extra_args": {"lens_potential_accuracy": 0}}},
        "likelihood": {"synthetic_zero": {"external": lik,
                                           "requires": {"Hubble": {"z": ZS}, "angular_diameter_distance": {"z": ZS},
                                                        "rdrag": None, "Cl": {"tt": 300}}}},
        "params": {"H0": 67.5, "ombh2": 0.0224, "omch2": 0.12, "mnu": 0.06, "tau": 0.055,
                   "As": 2.1e-9, "ns": 0.965, "card01_eps": {"prior": {"min": -1, "max": 1}},
                   "rdrag": None},
    }
    return get_model(info)


def direct(eps):
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.5, ombh2=0.0224, omch2=0.12, mnu=0.06, omk=0, tau=0.055)
    p.InitPower.set_params(As=2.1e-9, ns=0.965)
    if eps != 0:
        de = camb.dark_energy.DarkEnergyPPF(); de.set_w_a_table(T.A_TAB, T.w_table(eps)); p.DarkEnergy = de
    bg = camb.get_background(p)
    return np.array([bg.hubble_parameter(z) for z in ZS]), bg.angular_diameter_distance(ZS)


def main():
    m = make_model()
    for eps in [0.0, -0.2, 0.2]:
        m.logposterior({"card01_eps": eps})
        Hd, Dd = direct(eps)
        print(f"eps={eps:+.2f}: cobaya H={np.round(captured['H'],4)} DA={np.round(captured['DA'],3)} "
              f"rdrag={captured['rdrag']:.4f} TT220={captured['tt220']:.2f} | "
              f"max rel diff vs direct CAMB: H {np.max(abs(captured['H']/Hd-1)):.1e}, DA {np.max(abs(captured['DA']/Dd-1)):.1e}")


if __name__ == "__main__":
    main()
