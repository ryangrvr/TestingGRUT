"""Reproduce Card 1 and hostile checks; all results originate here.

Run from the repository: python research/card01/audit.py
Exact arithmetic uses stdlib Fraction, not a floating symbolic surrogate.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import expm

from law import channel, commutator, generator, spectral_blocks

ROOT = Path(__file__).resolve().parent
TAU = F(1, 5)
LEVELS = [F(0), F(1), F(3)]
NOISE_LEVELS = [F(0), F(2), F(1)]
PAIRS = [(0, 1), (0, 2), (1, 2)]
TIMES = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0]


@dataclass(frozen=True)
class Q:
    r: F = F(0)
    i: F = F(0)

    def __add__(self, other):
        if not isinstance(other, Q):
            other = Q(F(other))
        return Q(self.r + other.r, self.i + other.i)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.r, -self.i)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not isinstance(other, Q):
            other = Q(F(other))
        return Q(self.r * other.r - self.i * other.i,
                 self.r * other.i + self.i * other.r)

    __rmul__ = __mul__


def qmat(n):
    return [[Q() for _ in range(n)] for _ in range(n)]


def qmul(a, b):
    n = len(a)
    return [[sum((a[i][k] * b[k][j] for k in range(n)), Q())
             for j in range(n)] for i in range(n)]


def qcomm(a, b):
    ab, ba = qmul(a, b), qmul(b, a)
    return [[ab[i][j] - ba[i][j] for j in range(len(a))]
            for i in range(len(a))]


def qscale(s, a):
    return [[s * z for z in row] for row in a]


def qadd(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a))]
            for i in range(len(a))]


def exact_basis_controls():
    a = qmat(3)
    for i, level in enumerate(LEVELS):
        a[i][i] = Q(level)
    cells = []
    for i in range(3):
        for j in range(3):
            e = qmat(3)
            e[i][j] = Q(F(1))
            c = qcomm(a, e)
            actual = qadd(qscale(Q(F(0), F(-1)), c),
                          qscale(-TAU / 2, qcomm(a, c)))
            gap = LEVELS[i] - LEVELS[j]
            expected = qscale(Q(-TAU * gap**2 / 2, -gap), e)
            if actual != expected:
                raise AssertionError((i, j, actual, expected))
            cells.append({"i": i, "j": j,
                          "coefficient_real": str(actual[i][j].r),
                          "coefficient_imag": str(actual[i][j].i),
                          "off_cell_entries_zero": all(actual[k][l] == Q()
                             for k in range(3) for l in range(3) if (k, l) != (i, j))})
    return cells


def unit(n, i, j):
    e = np.zeros((n, n), complex)
    e[i, j] = 1
    return e


def supermatrix(fn, n):
    return np.column_stack([fn(unit(n, i, j)).reshape(-1)
                            for i in range(n) for j in range(n)])


def choi(fn, n):
    return sum((np.kron(unit(n, i, j), fn(unit(n, i, j)))
                for i in range(n) for j in range(n)), np.zeros((n*n, n*n), complex))


def baseline_channel(a, noise, x, t):
    gaps = a[:, None] - a[None, :]
    distances2 = sum((v[:, None] - v[None, :])**2 for v in noise)
    return np.exp((-1j*gaps - distances2/2)*t) * x


def collision_dilation(a, noise, t):
    """Exact Schur channel with one d-level zero-energy pure ancilla per interval.

    Return the controlled unitary, its system-ancilla isometry and Gram factor.
    Gram representation is exact analytically; eigendecomposition is numerical.
    Each controlled block's first column is the corresponding environment state.
    """
    n = len(a)
    c = baseline_channel(a, noise, np.ones((n, n)), t)
    eig, v = np.linalg.eigh(c)
    if eig.min() < -1e-12:
        raise AssertionError("not a CP Schur kernel")
    s = (v * np.sqrt(np.maximum(eig, 0))) @ v.conj().T
    # For row s_i, <e_j|e_i> = (s s^dagger)_ij = c_ij.
    states = s.T
    blocks = []
    for j in range(n):
        e = states[:, j]
        e = e / np.linalg.norm(e)
        columns = [e]
        for basis in np.eye(n, dtype=complex).T:
            x = basis.copy()
            for col in columns:
                x -= col * np.vdot(col, x)
            if np.linalg.norm(x) > 1e-10:
                columns.append(x / np.linalg.norm(x))
            if len(columns) == n:
                break
        blocks.append(np.column_stack(columns))
    u = sum((np.kron(unit(n, j, j), blocks[j]) for j in range(n)),
            np.zeros((n*n, n*n), complex))
    embedding = np.kron(np.eye(n), np.eye(n)[:, [0]])
    return u, u @ embedding, s


def partial_env(x, n):
    return np.trace(x.reshape(n, n, n, n), axis1=1, axis2=3)


def numerical_controls():
    a = np.diag([0., 1., 3.])
    tau = float(TAU)
    av = np.diag(a)
    models = {
        "candidate": [math.sqrt(tau) * av],
        "hostile": [math.sqrt(tau) * np.array([0., 2., 1.])],
    }
    matrix = supermatrix(lambda x: generator(a, x, tau), 3)
    worst = {k: 0.0 for k in ["expm", "trace", "energy", "semigroup", "unitary",
         "energy_commutator", "dilation_channel", "detailed_balance", "basis_covariance",
         "energy_shift", "ramsey", "fixed_algebra_projector"]}
    smallest_choi = 1.0
    rows = []
    fixture_rho = np.ones((3, 3), complex) / 3
    sigma = np.diag([0.2, 0.3, 0.5])  # a test weight, never an initial Gibbs premise
    for name, noise in models.items():
        for t in TIMES:
            phi = lambda x: baseline_channel(av, noise, x, t)
            c = choi(phi, 3)
            smallest_choi = min(smallest_choi, float(np.linalg.eigvalsh(c).min()))
            u, w, gram = collision_dilation(av, noise, t)
            worst["unitary"] = max(worst["unitary"], np.max(abs(u.conj().T@u-np.eye(9))))
            h = np.kron(a, np.eye(3))
            worst["energy_commutator"] = max(worst["energy_commutator"], np.max(abs(commutator(h, u))))
            for i in range(3):
                for j in range(3):
                    e = unit(3, i, j)
                    exact_channel = phi(e)
                    actual = partial_env(w @ e @ w.conj().T, 3)
                    worst["dilation_channel"] = max(worst["dilation_channel"], np.max(abs(actual-exact_channel)))
                    worst["trace"] = max(worst["trace"], abs(np.trace(exact_channel)-np.trace(e)))
                    worst["energy"] = max(worst["energy"], abs(np.trace(a@exact_channel)-np.trace(a@e)))
                    worst["semigroup"] = max(worst["semigroup"], np.max(abs(phi(phi(e))-baseline_channel(av, noise, e, 2*t))))
                    if name == "candidate":
                        via_expm = (expm(t*matrix) @ e.reshape(-1)).reshape(3, 3)
                        worst["expm"] = max(worst["expm"], np.max(abs(via_expm-channel(a, e, tau, t))))
            for i, j in PAIRS:
                rho = (unit(3,i,i)+unit(3,j,j)+unit(3,i,j)+unit(3,j,i))/2
                x = unit(3,i,j)+unit(3,j,i)
                y = -1j*unit(3,i,j)+1j*unit(3,j,i)
                r = phi(rho)
                z = np.trace(r@x)+1j*np.trace(r@y)
                rate = sum((v[i]-v[j])**2 for v in noise)/2
                expected = np.exp((-rate-1j*(av[j]-av[i]))*t)
                worst["ramsey"] = max(worst["ramsey"], abs(z-expected))
                rows.append({"model":name,"time":t,"pair":[i,j],"frequency":float(av[j]-av[i]),
                             "decay_rate":float(rate),"z_real":float(z.real),"z_imag":float(z.imag)})
        # Check dissipative GNS self-adjointness on all labelled matrix-unit pairs.
        noise_d = [np.diag(v) for v in noise]
        dissip = lambda x: -sum((commutator(d,commutator(d,x))/2 for d in noise_d), np.zeros_like(x))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        x,y = unit(3,i,j),unit(3,k,l)
                        err = abs(np.trace(sigma@x.conj().T@dissip(y))-np.trace(sigma@dissip(x).conj().T@y))
                        worst["detailed_balance"] = max(worst["detailed_balance"],err)

    # Same numeric assignment under exact representation moves, not a new law.
    v = np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
    ar = v@a@v.conj().T
    for t in TIMES:
        want = v@channel(a,fixture_rho,tau,t)@v.conj().T
        worst["basis_covariance"] = max(worst["basis_covariance"], np.max(abs(channel(ar,v@fixture_rho@v.conj().T,tau,t)-want)))
        worst["energy_shift"] = max(worst["energy_shift"], np.max(abs(channel(a+7*np.eye(3),fixture_rho,tau,t)-channel(a,fixture_rho,tau,t))))
    # Extractors use the initial A: no law evaluation or fitted cutoff occurs here.
    blocks = spectral_blocks(a)
    for _, p in blocks:
        worst["fixed_algebra_projector"] = max(worst["fixed_algebra_projector"], np.max(abs(generator(a,p,tau))))
    no_decay = supermatrix(lambda x: -1j*commutator(a,x),3)
    nullity = lambda x: int(np.sum(np.linalg.svd(x,compute_uv=False)<1e-11))
    if nullity(matrix)!=3 or nullity(no_decay)!=3:
        raise AssertionError("fixed-algebra ablation mismatch")
    for p in blocks:
        if np.max(abs(commutator(a,p[1])))>1e-12:
            raise AssertionError("identity decoder did not survive ablation")
    scalar = 2*np.eye(3)
    if len(spectral_blocks(scalar))!=1 or np.max(abs(generator(scalar,fixture_rho,tau)))>1e-12:
        raise AssertionError("scalar rejection fixture incorrect")

    # Same tau, second finite input: composite energies 0,1,1,2.
    # This extra evaluated instance is disclosed in IP12, never a candidate repair.
    aq = np.diag([0.,1.]); at = np.kron(aq,np.eye(2))+np.kron(np.eye(2),aq)
    rhoq = np.ones((2,2))/2; rhot = np.kron(rhoq,rhoq)
    t = math.pi/2
    rcol = channel(at,rhot,tau,t)
    rind = np.kron(channel(aq,rhoq,tau,t),channel(aq,rhoq,tau,t))
    xq = np.array([[0.,1.],[1.,0.]])
    xx = np.kron(xq,xq)
    collective = float(np.trace(rcol@xx).real)
    independent = float(np.trace(rind@xx).real)
    expected = (1-math.exp(-tau*math.pi))/2
    if abs(collective-expected)>1e-12 or abs(independent)>1e-12:
        raise AssertionError("composition identity failure")
    coherence = unit(4,1,2)
    if np.max(abs(generator(at,coherence,tau)))>1e-12:
        raise AssertionError("equal-energy coherence not protected")
    if max(worst.values())>1e-10 or smallest_choi<-1e-10:
        raise AssertionError((worst,smallest_choi))
    return {"worst_errors":{k:float(v) for k,v in worst.items()},
            "smallest_choi_eigenvalue":smallest_choi,"matrix_unit_interval_checks":108,
            "ramsey_cells":rows,"fixed_algebra_nullity_candidate":3,
            "fixed_algebra_nullity_without_dissipation":3,
            "scalar_input_blocks":1,"scalar_generator_zero":True,
            "composite_block_ranks":[int(round(np.trace(p).real)) for _,p in spectral_blocks(at)],
            "composition":{"time":t,"collective_XX":collective,"independent_XX":independent,
                 "analytic_collective_XX":expected,"protected_01_10_coherence":True,
                 "is_superluminal_signaling_claim":False},
            "finite_dilation":{"ancillas_per_interval":1,"ancilla_dimension":3,
                "ancilla_hamiltonian_zero":True,"ancilla_preparation":"pure blank, supplied",
                "same_resource_budget_both_models":True,"arbitrary_interventions_between_intervals":True}}


def exact_relations():
    rows=[]
    for i,j in PAIRS:
        omega=LEVELS[j]-LEVELS[i]
        gc=TAU*omega**2/2
        gh=TAU*(NOISE_LEVELS[j]-NOISE_LEVELS[i])**2/2
        rows.append({"pair":[i,j],"omega":str(omega),"candidate_rate":str(gc),
                     "countermodel_rate":str(gh),"candidate_tau_estimate":str(2*gc/omega**2),
                     "countermodel_tau_estimate":str(2*gh/omega**2)})
    w=[F(r['omega']) for r in rows]
    gc=[F(r['candidate_rate']) for r in rows]
    gh=[F(r['countermodel_rate']) for r in rows]
    rc=[gc[0]*w[k]**2-gc[k]*w[0]**2 for k in (1,2)]
    rh=[gh[0]*w[k]**2-gh[k]*w[0]**2 for k in (1,2)]
    if rc!=[F(0),F(0)] or rh!=[F(7,2),F(3,2)]:
        raise AssertionError((rc,rh))
    return {"rows":rows,"candidate_residuals":list(map(str,rc)),
            "countermodel_residuals":list(map(str,rh)),
            "baseline_rate_cone_dimension_at_fixed_A":3,"candidate_rate_ray_dimension":1,
            "conditional_codimension":2,"charter_earned_coupling_bits":0}


def main():
    freeze=json.loads((ROOT/'FREEZE.json').read_text())
    for name,digest in freeze['files'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError('frozen file changed: '+name)
    out={"status":"CARD1_ARITHMETIC_CONTROLS_PASS_CANDIDATE_KILL",
         "scientific_status":"AUTHOR_KILL_PENDING_INDEPENDENT_REVIEW",
         "budget":{"slots_consumed":1,"slots_remaining":2},
         "exact_basis_controls":exact_basis_controls(),"exact_restriction":exact_relations(),
         "numerical_controls":numerical_controls(),
         "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
         "frozen_files_verified":True}
    (ROOT/'results.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'])
    print(json.dumps({"exact_restriction":out['exact_restriction'],
                      "worst_errors":out['numerical_controls']['worst_errors']},indent=2))


if __name__=='__main__':
    main()
