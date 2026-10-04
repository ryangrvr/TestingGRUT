"""CARD-01 v1 Cobaya theory wrapper (FROZEN IMPLEMENTATION CANDIDATE; model-only test below).

Cobaya's stock CAMB theory accepts CPL (w, wa) but not a sampled tabulated w(a).  This module subclasses
cobaya.theories.camb.CAMB so that one extra sampled parameter `card01_eps` is consumed and the exact
Card #1 v1 law  w(a) = -1 + eps * R(z),  R = E_ref^2/(1+E_ref^2),  E_ref^2 = 0.31(1+z)^3 + 0.69
is passed to CAMB's DarkEnergyPPF via set_w_a_table.  The w(z) shape is NOT altered for software
convenience and is never replaced by its CPL projection.  eps == 0 uses CAMB's cosmological constant.

Usage in a Cobaya input:  theory: {card01_cobaya_theory.Card01CAMB: {python_path: <this dir>, ...}}
"""
import numpy as np
from cobaya.theories.camb import CAMB
from cobaya.theories.camb import camb as _cambmod

OMR, OLR = 0.31, 0.69
A_TAB = np.logspace(-9, 0, 3000)


def w_table(eps):
    z = 1.0 / A_TAB - 1.0
    e2 = OMR * (1 + z) ** 3 + OLR
    return -1.0 + eps * e2 / (1 + e2)


_orig = _cambmod.CambTransfers.get_can_support_params


def _transfers_support(self):
    s = _orig(self)
    if isinstance(self.cobaya_camb, Card01CAMB):
        s = set(s) | {"card01_eps"}
    return s


_cambmod.CambTransfers.get_can_support_params = _transfers_support


class Card01CAMB(CAMB):
    def set(self, params_values_dict, state):
        pv = dict(params_values_dict)
        eps = float(pv.pop("card01_eps", 0.0))
        p = super().set(pv, state)
        if p and eps != 0.0:
            de = self.camb.dark_energy.DarkEnergyPPF()
            de.set_w_a_table(A_TAB, w_table(eps))
            p.DarkEnergy = de
        return p
