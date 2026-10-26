#!/usr/bin/env python3
"""kin_rates (Claude, own code): tree amplitudes of the light sector for any number of legs, vectorised over
kinematic configurations, with complex (regulated) propagator masses and per-line propagator switches.

Method: Berends-Giele off-shell recursion in the real field basis (x1, x2, s) with the exact cubic and quartic
vertex tensors of U (model.py).  The method is the one used by our gp_leftovers track (amp42.py); the code here is
rewritten: batch dimension, complex masses, propagator switches for the real-intermediate-state subtraction, and
analytic 2->2 references for every channel used below.  Convention: all momenta incoming, metric (+,-,-,-),
M = -prod(V) prod 1/(p^2 - m^2 + i m Gamma) over tree diagrams; an outgoing particle is an incoming antiparticle
with momentum -k.  Legs 0..n-2 are recursed, leg n-1 is the off-shell leg; an internal line is a subset S of
{0..n-2} with 2 <= |S| <= n-2.
"""
import itertools
import numpy as np
import model as M

METRIC = np.array([1.0, -1.0, -1.0, -1.0])


def minkowski_sq(p):
    """p: (..., 4) -> p^2 (...)."""
    return p[..., 0]**2 - p[..., 1]**2 - p[..., 2]**2 - p[..., 3]**2


def dot(p, q):
    return p[..., 0] * q[..., 0] - p[..., 1] * q[..., 1] - p[..., 2] * q[..., 2] - p[..., 3] * q[..., 3]


def _partitions2(S):
    """Unordered partitions of bitmask S into two non-empty parts (A contains the lowest bit)."""
    low = S & -S
    sub = (S - 1) & S
    out = []
    while sub:
        if sub & low:
            out.append((sub, S ^ sub))
        sub = (sub - 1) & S
    return out


def _partitions3(S):
    out = []
    for A, R in _partitions2(S):
        for B, C in _partitions2(R):
            out.append((A, B, C))
    return out


class Amplitude:
    """Amplitude evaluator for a fixed external leg list.

    legs: list of species strings 'X', 'Xb', 's' for incoming legs and ('X','out') etc. for outgoing; momenta are
    given physical (positive energy) and the class flips outgoing momenta internally.
    m2c: complex mass^2 array for (x1, x2, s): m^2 - i m Gamma_reg.
    """

    def __init__(self, legs, m2c=None, v=M.V, widths=(0.0, 0.0)):
        self.legs = legs
        self.n = len(legs)
        m2, V3, V4 = M.get_vertices(v)
        self.m2 = m2
        self.set_vertices(V3, V4)
        if m2c is None:
            gx, gs = widths
            m2c = np.array([m2[0] - 1j * np.sqrt(m2[0]) * gx, m2[1] - 1j * np.sqrt(m2[1]) * gx,
                            m2[2] - 1j * np.sqrt(m2[2]) * gs])
        self.m2c = np.asarray(m2c, dtype=complex)
        self.pol = []
        self.sign = []
        self.charge = []
        for leg in legs:
            if isinstance(leg, tuple):
                sp, io = leg
            else:
                sp, io = leg, 'in'
            if io == 'in':
                self.pol.append(M.POL_IN[sp]); self.sign.append(1.0); self.charge.append(M.CHARGE[sp])
            else:
                self.pol.append(M.POL_OUT[sp]); self.sign.append(-1.0); self.charge.append(-M.CHARGE[sp])
        self.sign = np.array(self.sign)
        nrec = self.n - 1
        self.full = (1 << nrec) - 1
        # subsets by size, with partitions precomputed
        self.subsets = {}
        for S in range(1, self.full + 1):
            k = bin(S).count('1')
            if k >= 2:
                self.subsets[S] = (_partitions2(S), _partitions3(S))
        self.order = sorted(self.subsets, key=lambda S: bin(S).count('1'))
        self.internal = [S for S in self.order if bin(S).count('1') <= self.n - 2]
        # charge of each subset and the propagating field (0/1 -> x, 2 -> s, None)
        self.subset_charge = {S: sum(self.charge[i] for i in range(nrec) if S >> i & 1) for S in self.order}

    def set_vertices(self, V3, V4):
        """Store the vertex tensors as lists of their non-zero entries (7 cubic, 21 quartic) so that each
        contraction is a handful of vectorised products: K[n,a] = sum V3[a,b,c] J1[n,b] J2[n,c]."""
        self.V3, self.V4 = V3, V4
        self.nz3 = {a: [(b, c, V3[a, b, c]) for b in range(3) for c in range(3) if V3[a, b, c] != 0] for a in range(3)}
        self.nz4 = {a: [(b, c, d, V4[a, b, c, d]) for b in range(3) for c in range(3) for d in range(3)
                        if V4[a, b, c, d] != 0] for a in range(3)}

    def _c3(self, J1, J2):
        K = np.zeros_like(J1)
        for a, lst in self.nz3.items():
            for b, c, v in lst:
                K[:, a] += v * J1[:, b] * J2[:, c]
        return K

    def _c4(self, J1, J2, J3):
        K = np.zeros_like(J1)
        for a, lst in self.nz4.items():
            for b, c, d, v in lst:
                K[:, a] += v * J1[:, b] * J2[:, c] * J3[:, d]
        return K

    def propagating(self, S):
        q = self.subset_charge[S]
        if q == 0:
            return 's'
        if abs(q) == 1:
            return 'X'
        return None

    def __call__(self, momenta, switch=None, return_invariants=False):
        """momenta: (N, n, 4) physical momenta of the legs; switch: dict {S: alpha} scaling propagator of line S.
        Returns complex M (N,)."""
        p = np.asarray(momenta, dtype=float) * self.sign[None, :, None]
        N = p.shape[0]
        nrec = self.n - 1
        J = {}
        for i in range(nrec):
            J[1 << i] = np.broadcast_to(self.pol[i], (N, 3)).astype(complex)
        PS = {}
        for i in range(nrec):
            PS[1 << i] = p[:, i, :]
        inv = {}
        for S in self.order:
            # momentum of subset
            if S not in PS:
                low = S & -S
                PS[S] = PS[low] + PS[S ^ low]
            P2 = minkowski_sq(PS[S])
            inv[S] = P2
            K = np.zeros((N, 3), dtype=complex)
            parts2, parts3 = self.subsets[S]
            for A, B in parts2:
                K += self._c3(J[A], J[B])
            for A, B, C in parts3:
                K += self._c4(J[A], J[B], J[C])
            if S == self.full:
                amp = -np.einsum('a,na->n', self.pol[nrec], K)
                if return_invariants:
                    return amp, inv
                return amp
            D = 1.0 / (P2[:, None] - self.m2c[None, :])
            if switch and S in switch:
                D = D * switch[S]
            J[S] = K * D


# ------------------------------------------------------------------ analytic 2->2 references [D]
def mandelstam(p1, p2, k1, k2):
    s = minkowski_sq(p1 + p2); t = minkowski_sq(p1 - k1); u = minkowski_sq(p1 - k2)
    return s, t, u


def M_XXb_XXb(s, t, u, mX2=M.MX**2, ms2c=M.MS**2, lam3=M.LAM3, lX=M.LX):
    return -4 * lX - lam3**2 * (1 / (s - ms2c) + 1 / (t - ms2c))


def M_XX_XX(s, t, u, mX2=M.MX**2, ms2c=M.MS**2, lam3=M.LAM3, lX=M.LX):
    return -4 * lX - lam3**2 * (1 / (t - ms2c) + 1 / (u - ms2c))


def M_Xs_Xs(s, t, u, mX2c=M.MX**2, ms2c=M.MS**2, lam3=M.LAM3, lam_s3=M.LAM_S3, g2=M.GX**2):
    return -2 * g2 - lam3**2 * (1 / (s - mX2c) + 1 / (u - mX2c)) - lam3 * lam_s3 / (t - ms2c)


def M_ss_ss(s, t, u, ms2c=M.MS**2, lam_s3=M.LAM_S3, lam_s4=M.LAM_S4):
    return -lam_s4 - lam_s3**2 * (1 / (s - ms2c) + 1 / (t - ms2c) + 1 / (u - ms2c))


def M_ss_XXb(s, t, u, mX2c=M.MX**2, ms2c=M.MS**2, lam3=M.LAM3, lam_s3=M.LAM_S3, g2=M.GX**2):
    return -2 * g2 - lam_s3 * lam3 / (s - ms2c) - lam3**2 * (1 / (t - mX2c) + 1 / (u - mX2c))


# ------------------------------------------------------------------ kinematics helpers
def cm_momenta(sqrt_s, m1, m2, m3, m4, cos_th, phi=0.0):
    """2->2 CM momenta for arbitrary arrays of cos_th (N,). Returns p1,p2,k1,k2 each (N,4)."""
    s = sqrt_s**2
    p = np.sqrt(kallen(s, m1 * m1, m2 * m2)) / (2 * sqrt_s)
    k = np.sqrt(kallen(s, m3 * m3, m4 * m4)) / (2 * sqrt_s)
    E1 = (s + m1 * m1 - m2 * m2) / (2 * sqrt_s); E2 = sqrt_s - E1
    E3 = (s + m3 * m3 - m4 * m4) / (2 * sqrt_s); E4 = sqrt_s - E3
    c = np.asarray(cos_th, dtype=float); N = c.shape[0]
    sn = np.sqrt(np.clip(1 - c * c, 0, None))
    z = np.zeros(N); o = np.ones(N)
    p1 = np.stack([E1 * o, z, z, p * o], 1); p2 = np.stack([E2 * o, z, z, -p * o], 1)
    k1 = np.stack([E3 * o, k * sn * np.cos(phi), k * sn * np.sin(phi), k * c], 1)
    k2 = np.stack([E4 * o, -k * sn * np.cos(phi), -k * sn * np.sin(phi), -k * c], 1)
    return p1, p2, k1, k2


def kallen(a, b, c):
    return a * a + b * b + c * c - 2 * a * b - 2 * a * c - 2 * b * c


# ------------------------------------------------------------------ self-test
def _test():
    rng = np.random.default_rng(3)
    out = {}
    cases = [
        ('XXb_XXb', ['X', 'Xb', ('X', 'out'), ('Xb', 'out')], (M.MX, M.MX, M.MX, M.MX), M_XXb_XXb),
        ('XX_XX', ['X', 'X', ('X', 'out'), ('X', 'out')], (M.MX, M.MX, M.MX, M.MX), M_XX_XX),
        ('Xs_Xs', ['X', 's', ('X', 'out'), ('s', 'out')], (M.MX, M.MS, M.MX, M.MS), M_Xs_Xs),
        ('ss_ss', ['s', 's', ('s', 'out'), ('s', 'out')], (M.MS, M.MS, M.MS, M.MS), M_ss_ss),
        ('ss_XXb', ['s', 's', ('X', 'out'), ('Xb', 'out')], (M.MS, M.MS, M.MX, M.MX), M_ss_XXb),
    ]
    for name, legs, ms, fan in cases:
        A = Amplitude(legs)
        sq = max(ms[0] + ms[1], ms[2] + ms[3]) * (1.01 + 2 * rng.random(7))
        c = rng.uniform(-1, 1, 7)
        res = []
        for sqs, cc in zip(sq, c):
            p1, p2, k1, k2 = cm_momenta(sqs, *ms, np.array([cc]), phi=0.4)
            mom = np.stack([p1, p2, k1, k2], 1)
            m = A(mom)[0]
            s, t, u = mandelstam(p1[0], p2[0], k1[0], k2[0])
            res.append(abs(m - fan(s, t, u)) / abs(fan(s, t, u)))
        out[name] = max(res)
    # charge violation X X -> X Xbar
    A = Amplitude(['X', 'X', ('X', 'out'), ('Xb', 'out')])
    p1, p2, k1, k2 = cm_momenta(3 * M.MX, M.MX, M.MX, M.MX, M.MX, np.array([0.3]))
    out['charge_violating'] = abs(A(np.stack([p1, p2, k1, k2], 1))[0])
    # 1->2 decay amplitude
    A = Amplitude(['s', ('X', 'out'), ('Xb', 'out')])
    pst = np.sqrt(M.MS**2 / 4 - M.MX**2)
    mom = np.array([[[M.MS, 0, 0, 0], [M.MS / 2, 0, 0, pst], [M.MS / 2, 0, 0, -pst]]])
    out['s_decay_M_over_lam3'] = complex(A(mom)[0]) / M.LAM3
    # threshold 4->2 (reproduction of our gp_leftovers numbers; own code)
    m = M.MX; kf = np.sqrt(3) * m
    rest = np.array([m, 0, 0, 0.])
    k1 = np.array([2 * m, 0, 0, kf]); k2 = np.array([2 * m, 0, 0, -kf])
    Aa = Amplitude(['X', 'X', 'Xb', 'Xb', ('X', 'out'), ('Xb', 'out')])
    Ab = Amplitude(['X', 'X', 'X', 'Xb', ('X', 'out'), ('X', 'out')])
    mom = np.array([[rest, rest, rest, rest, k1, k2]])
    out['M_a_threshold'] = complex(Aa(mom)[0]).real
    out['M_b_threshold'] = complex(Ab(mom)[0]).real
    # lambda_X-only check: 12 lX^2 / m^2 (hand diagram sum of our earlier track) via a modified vertex set
    import copy
    Al = copy.copy(Aa)
    V3 = np.zeros((3, 3, 3)); V4 = np.zeros((3, 3, 3, 3))
    V4[0, 0, 0, 0] = V4[1, 1, 1, 1] = 6 * M.LX
    for perm in set(itertools.permutations([0, 0, 1, 1])):
        V4[perm] = 2 * M.LX
    Al.set_vertices(V3, V4)
    out['M_a_lambdaX_only'] = complex(Al(mom)[0]).real
    out['M_a_lambdaX_only_expected_12lX2_over_m2'] = 12 * M.LX**2 / m**2
    return out


if __name__ == '__main__':
    import json
    print(json.dumps(_test(), indent=1, default=str))
