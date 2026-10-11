"""Independent CP-09 stationary-clock spectral budget control (NumPy/SciPy).

Construct N=2 CP-09 without importing its code. This tests the stationary
noninteracting Page-Wootters history construction, NOT preparation by the
previous autonomous clock-pointer apparatus.
"""
import math
import numpy as np
from scipy.linalg import eigh, expm

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def cp09_n2():
    gam = math.pi / 4
    om = .25
    P01, P10 = np.kron(Z, X), np.kron(X, Z)
    U01 = math.cos(gam)*np.eye(4) - 1j*math.sin(gam)*P01
    U10 = math.cos(gam)*np.eye(4) - 1j*math.sin(gam)*P10
    onsite = om*(2*np.eye(4)-np.kron(X, I)-np.kron(I, X))
    H = np.zeros((12, 12), complex)
    for k, U in [(1, U01), (2, U10)]:
        H[:4, :4] += np.eye(4)
        H[4*k:4*(k+1), 4*k:4*(k+1)] += np.eye(4)
        H[4*k:4*(k+1), :4] -= U
        H[:4, 4*k:4*(k+1)] -= U.conj().T
    for k in range(3):
        H[4*k:4*(k+1), 4*k:4*(k+1)] += onsite
    A = np.zeros_like(H)
    B = np.zeros_like(H)
    # Geometry: empty, 0<1, 1<0. Disjoint relational source/detector.
    A[4:8, 4:8] = np.kron(I, Z)
    A[8:12, 8:12] = np.kron(Z, I)
    B[4:8, 4:8] = np.kron(X, I)
    B[8:12, 8:12] = np.kron(I, X)
    return H, A, B


def conditioned_clock(psi, t, E, V, chosen):
    """Stationary joint history with q EXACTLY tuned clock energy levels."""
    q = len(chosen)
    amps = V[:, chosen].conj().T @ psi
    joint = amps[:, None]*V[:, chosen].T
    Ec = -E[chosen]
    residual = np.linalg.norm(
        joint @ (V@np.diag(E)@V.conj().T).T + Ec[:, None]*joint)
    assert residual < 1e-11
    # clock ket |t> = q^(-1/2) sum exp(-i Ec t)|a>.
    # bra has phase exp(+i Ec t).
    raw = (np.exp(1j*Ec*t) @ joint)/np.sqrt(q)
    p = np.vdot(raw, raw).real
    return raw/np.sqrt(p), p, residual


def check():
    H, A, B = cp09_n2()
    E, V = eigh(H)
    assert abs(E[0]-.09475528369770218)<1e-11
    assert abs(E[1]-E[0]-.31829854178778927)<1e-11
    g = V[:, 0]
    b = V.conj().T @ (B@g)
    modes = [k for k in range(1, 12) if abs(b[k])>1e-9]
    assert len(modes)==6, modes
    assert all(abs(b[k])<1e-11 for k in range(1, 12) if k not in modes)
    pop = np.array([abs(b[k])**2 for k in modes])
    assert abs(pop.sum()-.4406182061422456)<1e-11
    assert abs(pop.max()-.1783721202157285)<1e-11
    min_tangent_error = math.sqrt(float(pop.sum()-pop.max()))
    assert abs(min_tangent_error-.5120996835836915)<1e-11
    print('occupied first-order system energies:', [0]+modes)
    print('excited weights:', [(k, round(float(abs(b[k])**2), 10)) for k in modes])
    print('two-level clock, keeping ground: minimum excluded tangent norm=',min_tangent_error)
    # For the selected impulse, B^2 and all powers retain the same sector.
    Ps = V[:, [0]+modes]@V[:, [0]+modes].conj().T
    assert np.linalg.norm((np.eye(12)-Ps)@(B@B@g)) < 1e-11
    eps=1e-5
    plus=expm(1j*eps*B)@g
    minus=expm(-1j*eps*B)@g
    assert np.linalg.norm((np.eye(12)-Ps)@plus)<1e-11
    assert np.linalg.norm((np.eye(12)-Ps)@minus)<1e-11
    for t in (.1,.4,1.):
        plus_cond,pplus,resplus = conditioned_clock(plus,t,E,V,[0]+modes)
        minus_cond,pminus,resminus = conditioned_clock(minus,t,E,V,[0]+modes)
        directplus=expm(-1j*H*t)@plus
        directminus=expm(-1j*H*t)@minus
        assert abs(pplus-1/7)<1e-11 and abs(pminus-1/7)<1e-11
        assert np.linalg.norm(plus_cond-directplus)<1e-11
        assert np.linalg.norm(minus_cond-directminus)<1e-11
        response=(np.vdot(plus_cond,A@plus_cond).real -
                  np.vdot(minus_cond,A@minus_cond).real)/(2*eps)
        U=expm(-1j*H*t)
        At=U.conj().T@A@U
        ideal=float((1j*np.vdot(g,(At@B-B@At)@g)).real)
        assert abs(response-ideal)<1e-8
        print(f'q=7 clock: t={t}, response={response:.12f}, ideal={ideal:.12f}, p={pplus:.12f}')
    # Explicit q=2 spectral truncation, NOT a bound on response error:
    # retain ground plus largest B|g> excited sector.
    best=modes[int(np.argmax(pop))]
    P2=V[:, [0,best]]@V[:, [0,best]].conj().T
    plus2=P2@plus;plus2/=np.linalg.norm(plus2)
    minus2=P2@minus;minus2/=np.linalg.norm(minus2)
    t=.4; U=expm(-1j*H*t)
    R2=(np.vdot(U@plus2,A@(U@plus2)).real-
        np.vdot(U@minus2,A@(U@minus2)).real)/(2*eps)
    print(f'q=2 spectral truncation (ground + strongest): response at .4={R2:.12f}, '
          f'vs ideal 0.052588513539')
    assert abs(R2-.12983643596)<1e-8
    print('PASS: seven matched clock energies recover exact response; '
          'two cannot reproduce the full first-order clock history')


if __name__ == '__main__':
    check()
