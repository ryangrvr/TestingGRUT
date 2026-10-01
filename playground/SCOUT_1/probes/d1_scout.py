"""SCOUT-1 D1-SCOUT (pin = 0 + extended nonlinear dynamics), SCOUT-ONLY premise experiment.

Part A (D1-4 / D1-6): the INHERITED S2 C-B drift dx = [-K x - 4 beta x^3] dt + B dW extended on a pin-free chain
  (K = Laplacian), uniform T (equilibrium). Exact 1D transfer operator of exp(-V/T),
  V = sum (x_{i+1}-x_i)^2 / 2 + beta sum x_i^4:  correlation length xi = 1/ln(l0/l1).
  beta = 0: zero mode, xi = infinity (EW interface, IR visible). beta > 0: does the on-site quartic regenerate a mass?

Part B (D1-1 / D1-3): conserved scalar u on a ring, bond (conserved) noise,
  du/dt = -d_x F,  F_{i+1/2} = -(u_{i+1} - u_i) + J(ubar) + sqrt(2T/dt) eta + [Model-B term] ,
  operator content scanned:  J = 0 | eps*u (affinity, constant mobility) | eps*u^2 (lambda_2, density-dependent
  mobility) | eps*u^3 (Z2-odd cubic current) | Model B: F += -d_x(4 beta u^3) (inherited quartic made conservative).
  Measured: growth exponent of the integrated height h_i = sum_{j<=i} u_j  (EW 1/4, KPZ 1/3).

Part C (D1-4): pin > 0 (non-conserving relaxation -r u) under the same lambda_2 dynamics: width saturation.
"""
import numpy as np

rng = np.random.default_rng(41)


# ---------------- Part A ----------------
def transfer_xi(beta, T, xmax=8.0, n=1201):
    x = np.linspace(-xmax, xmax, n)
    dx = x[1] - x[0]
    X, Y = np.meshgrid(x, x, indexing="ij")
    Kmat = np.exp(-((X - Y) ** 2 / 2 + beta * (X ** 4 + Y ** 4) / 2) / T) * dx
    ev = np.sort(np.abs(np.linalg.eigvalsh(Kmat)))[::-1]
    return 1 / np.log(ev[0] / ev[1])


def partA():
    print("=== A. inherited C-B quartic on a pin-free chain (exact transfer operator, T = 1) ===")
    for beta in (1.0, 0.1, 0.01, 0.001, 1e-4):
        xm = max(8.0, 4.0 * (1 / beta) ** 0.25 if beta > 0 else 8.0)
        xm = min(xm, 60.0)
        xi = transfer_xi(beta, 1.0, xmax=xm, n=1601)
        print(f"  beta = {beta:7.4f}: correlation length xi = {xi:9.3f}   (xi * beta^(1/3) = {xi * beta ** (1 / 3):.3f})")
    print("  beta = 0: Gaussian random-walk interface, xi = infinity (pin = 0 exposes the zero mode)")
    print("  => the inherited on-site quartic regenerates a mass xi ~ beta^(-1/3) (1D: no finite-T transition):")
    print("     with the S2 nonlinearity present, pin = 0 does NOT produce a critical IR at all.")


# ---------------- Part B / C ----------------
def conserved_run(N, tmax, dt, T, J, eps, beta=0.0, r=0.0, reps=4):
    times = np.unique(np.round(np.logspace(0, np.log10(tmax), 24) / dt).astype(int))
    W2 = np.zeros(len(times))
    for _ in range(reps):
        u = np.zeros(N)
        k = 0
        for s in range(1, times[-1] + 1):
            up = np.roll(u, -1)
            ub = 0.5 * (u + up)
            F = -(up - u) + eps * J(ub) + np.sqrt(2 * T / dt) * rng.standard_normal(N)
            if beta:
                mu = 4 * beta * u ** 3
                F += -(np.roll(mu, -1) - mu)
            u = u - dt * (F - np.roll(F, 1))
            if r:
                u = u - dt * r * u + np.sqrt(2 * T * r * dt) * 0.0
            if s == times[k]:
                h = np.cumsum(u)
                W2[k] += np.var(h)
                k += 1
    return times * dt, np.sqrt(W2 / reps)


def exponents(t, W, wins=((5, 50), (50, 500), (500, 2000))):
    out = []
    for a, b in wins:
        m = (t >= a) & (t <= b)
        if m.sum() >= 3:
            out.append(np.polyfit(np.log(t[m]), np.log(W[m]), 1)[0])
        else:
            out.append(np.nan)
    return out


def partB():
    print("\n=== B. conserved scalar, operator content scanned (N = 4096, T = 1, dt = 0.02, t <= 2000) ===")
    cases = [
        ("J = 0 (reciprocal, E-1)", lambda v: 0 * v, 0.0, 0.0),
        ("J = eps u (affinity E-4, constant mobility)", lambda v: v, 1.0, 0.0),
        ("Model B: inherited quartic made conservative", lambda v: 0 * v, 0.0, 0.25),
        ("affinity + Model-B quartic", lambda v: v, 1.0, 0.25),
        ("J = eps u^3 (Z2-odd cubic current, SUPPLIED)", lambda v: v ** 3, 1.0, 0.0),
        ("J = eps u^2 (lambda_2, SUPPLIED mobility)", lambda v: v ** 2, 1.0, 0.0),
    ]
    for name, J, eps, beta in cases:
        t, W = conserved_run(4096, 2000, 0.02, 1.0, J, eps, beta=beta)
        e = exponents(t, W)
        print(f"  {name:48s}: growth exponent [5-50, 50-500, 500-2000] = " + ", ".join(f"{x:.3f}" for x in e))
    print("  references: EW 1/4, KPZ 1/3")


def partC():
    print("\n=== C. pin role under the same lambda_2 dynamics ===")
    for r in (0.0, 0.01, 0.1):
        t, W = conserved_run(4096, 2000, 0.02, 1.0, lambda v: v ** 2, 1.0, r=r, reps=3)
        print(f"  r = {r:5.2f}: W(t) at t = 10, 100, 1000, 2000 = {np.round(np.interp([10, 100, 1000, 2000], t, W), 3)}")
    print("  pin > 0 cuts the scaling off (W saturates); pin = 0 only lets the supplied operator's class become visible")


if __name__ == "__main__":
    import sys
    part = sys.argv[1] if len(sys.argv) > 1 else "ABC"
    if "A" in part: partA()
    if "B" in part: partB()
    if "C" in part: partC()


def cubic_rerun():
    """Z2-odd cubic current with a Lax-Friedrichs (monotone) flux; the centred flux of part B was unstable."""
    print("\n=== B'. J = eps u^3 rerun with Lax-Friedrichs flux (eps = 0.2), and lambda_2 reference with the same flux ===")
    for name, J, Jp in (("J = 0.2 u^3 (cubic, Z2-odd)", lambda v: 0.2 * v ** 3, lambda v: 0.6 * v ** 2),
                        ("J = 0.2 u^2 (lambda_2)", lambda v: 0.2 * v ** 2, lambda v: 0.4 * np.abs(v))):
        N, dt, T, tmax, reps = 4096, 0.02, 1.0, 2000, 4
        times = np.unique(np.round(np.logspace(0, np.log10(tmax), 24) / dt).astype(int))
        W2 = np.zeros(len(times))
        for _ in range(reps):
            u = np.zeros(N); k = 0
            for s in range(1, times[-1] + 1):
                up = np.roll(u, -1)
                a = np.maximum(Jp(u), Jp(up))
                FJ = 0.5 * (J(u) + J(up)) - 0.5 * a * (up - u)
                F = -(up - u) + FJ + np.sqrt(2 * T / dt) * rng.standard_normal(N)
                u = u - dt * (F - np.roll(F, 1))
                if s == times[k]:
                    W2[k] += np.var(np.cumsum(u)); k += 1
        t = times * dt; W = np.sqrt(W2 / reps)
        print(f"  {name:30s}: growth exponent [5-50, 50-500, 500-2000] = " + ", ".join(f"{x:.3f}" for x in exponents(t, W)))
