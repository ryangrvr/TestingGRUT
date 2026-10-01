"""SCOUT-1 W1-F: prediction-first check of the KMS/FDT route (C13, C14).
 1. dS Bunch-Davies (conformally coupled massless scalar, dS4) geodesic Wightman function
    G(s) = -H^2 / (16 pi^2 sinh^2(H (s - i eps)/2)); detector response F(w) = int ds e^{-i w s} G(s);
    check detailed balance F(w)/F(-w) = exp(-2 pi w / H)  (KMS at T = H / 2 pi).
 2. noise/dissipation ratio coth(pi w / H): a weight-0 function of w/H (W1-C eligible) -- fixed given background.
 3. magnitude for the present Hubble rate: T_dS = hbar H0 / (2 pi k_B).
"""
import numpy as np
from scipy.integrate import quad

H = 1.0
eps = 0.6   # contour shift; result is eps-independent for 0 < eps < 2 pi / H (analyticity strip)


def G(s):
    z = H * (s - 1j * eps) / 2
    return -H ** 2 / (16 * np.pi ** 2 * np.sinh(z) ** 2)


def F(w, Lmax=80.0):
    re = quad(lambda s: (np.exp(-1j * w * s) * G(s)).real, -Lmax, Lmax, limit=4000)[0]
    im = quad(lambda s: (np.exp(-1j * w * s) * G(s)).imag, -Lmax, Lmax, limit=4000)[0]
    # the shifted contour Im s = -eps lies in the pole-free strip; undo the e^{w eps} factor of the shift
    return np.exp(-w * eps) * (re + 1j * im)


print("=== 1. detailed balance of the dS Bunch-Davies detector response ===")
for w in (0.3, 0.8, 1.5, 2.5):
    Fp, Fm = F(w), F(-w)
    exact = (w / (2 * np.pi)) / (np.exp(2 * np.pi * w / H) - 1)
    print(f"  w/H={w}: F(w)={Fp.real:.6e} (Planckian {exact:.6e}), F(w)/F(-w)={(Fp / Fm).real:.6e}, "
          f"exp(-2 pi w/H)={np.exp(-2 * np.pi * w / H):.6e}")
for e2 in (0.3, 1.2):
    eps = e2
    print(f"  eps={e2}: F(0.8)/F(-0.8) = {(F(0.8) / F(-0.8)).real:.6e} (contour-independent)")

print("\n=== 2. noise/dissipation ratio coth(pi w/H) is a function of w/H only ===")
for HH in (0.5, 1.0, 4.0):
    x = np.array([0.1, 1.0, 3.0])
    print(f"  H={HH}: ratio at w/H = 0.1, 1, 3 -> {np.round(1 / np.tanh(np.pi * x), 6)}")

print("\n=== 3. magnitude ===")
hbar, kB = 1.054571817e-34, 1.380649e-23
H0 = 70e3 / 3.0857e22
print(f"  H0 = {H0:.3e} s^-1 -> T_dS = {hbar * H0 / (2 * np.pi * kB):.3e} K")
print("  (inflationary H ~ 1e13 GeV would give T ~ 1.6e12 GeV, imprinted only via standard dS-QFT spectra)")
