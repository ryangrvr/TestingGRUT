import numpy as np
n = 240; L = 6.0
xg, wg = np.polynomial.legendre.leggauss(n)
xs = L * xg; wx = L * wg
X0 = np.broadcast_to(xs[:, None], (n, n)).copy()
P0 = np.broadcast_to(xs[None, :], (n, n)).copy()
rho = np.exp(-0.5*P0**2 - 0.5*X0**2 - 0.25*X0**4)
W = wx[:, None] * wx[None, :]
Wn = W * rho; Wn /= Wn.sum()
def s(u): return 10*u**3 - 15*u**4 + 6*u**5
def q_p1(t): return s(t/np.pi) if t <= np.pi else 1.0
def integrate(X0, P0, T, eps, dt):
    steps = max(1, int(round(T/dt))); h = T/steps
    x = X0.copy(); v = P0.copy()
    for k in range(steps):
        q = q_p1(k*h)
        k1v = -x - x**3 + eps*q; k1x = v
        x2 = x + 0.5*h*k1x; v2 = v + 0.5*h*k1v
        qm = q_p1(k*h + 0.5*h)
        k2v = -x2 - x2**3 + eps*qm; k2x = v2
        x3 = x + 0.5*h*k2x; v3 = v + 0.5*h*k2v
        k3v = -x3 - x3**3 + eps*qm; k3x = v3
        x4 = x + h*k3x; v4 = v + h*k3v
        q1 = q_p1(k*h + h)
        k4v = -x4 - x4**3 + eps*q1; k4x = v4
        x = x + (h/6)*(k1x+2*k2x+2*k3x+k4x)
        v = v + (h/6)*(k1v+2*k2v+2*k3v+k4v)
    return x
dt = 5e-4
for T in [0.25, 0.5, 0.75, 1.0]:
    row = []
    for nb in [4, 8, 16, 32, 64, 128]:
        eps = 1.0/np.sqrt(nb)
        xT = integrate(X0, P0, T, eps, dt)
        m1 = (xT*Wn).sum(); m2 = (xT**2*Wn).sum(); m3 = (xT**3*Wn).sum()
        var = m2 - m1*m1; k3 = m3 - 3*m1*m2 + 2*m1*m1*m1
        k3F = k3/np.sqrt(nb); g1F = k3F/var**1.5
        row.append(g1F)
    print(f't={T}: ' + ' '.join(f'{g:.6e}' for g in row))
gam = []
for nb in [4, 8, 16, 32, 64, 128]:
    eps = 1.0/np.sqrt(nb)
    xT = integrate(X0, P0, 0.5, eps, dt)
    m1 = (xT*Wn).sum(); m2 = (xT**2*Wn).sum(); m3 = (xT**3*Wn).sum()
    var = m2 - m1*m1; k3 = m3 - 3*m1*m2 + 2*m1*m1*m1
    gam.append(k3/np.sqrt(nb)/var**1.5)
gam = np.array(gam)
nb_arr = np.array([4.,8.,16.,32.,64.,128.])
p, a = np.polyfit(np.log(nb_arr), np.log(np.abs(gam)), 1)
print(f'Fitted p at t*=0.5: {-p:.6f}')
print(f'N_B*gamma1: ' + ' '.join(f'{nb*g:.4e}' for nb,g in zip(nb_arr,gam)))
# odd-epsilon check
for nb in [4, 16]:
    eps = 1.0/np.sqrt(nb)
    xT_p = integrate(X0, P0, 0.5, eps, dt)
    xT_m = integrate(X0, P0, 0.5, -eps, dt)
    m1p=(xT_p*Wn).sum(); m2p=(xT_p**2*Wn).sum(); m3p=(xT_p**3*Wn).sum()
    k3p = m3p-3*m1p*m2p+2*m1p*m1p*m1p
    m1m=(xT_m*Wn).sum(); m2m=(xT_m**2*Wn).sum(); m3m=(xT_m**3*Wn).sum()
    k3m = m3m-3*m1m*m2m+2*m1m*m1m*m1m
    print(f'odd-eps nb={nb}: k3(+eps)={k3p:.6e}, k3(-eps)={k3m:.6e}, sum={k3p+k3m:.2e}, ratio={k3m/k3p:.6f}')
