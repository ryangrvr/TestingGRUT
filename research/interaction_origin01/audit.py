"""Finite conventional completion null. No GRUT candidate is searched or scored."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent
TOL = 4e-11
I = np.eye(9, dtype=complex)
b = np.array([[0, 1, 0], [0, 0, np.sqrt(2)], [0, 0, 0]], complex)
b1, b2 = np.kron(b, np.eye(3)), np.kron(np.eye(3), b)
n1, n2 = b1.conj().T @ b1, b2.conj().T @ b2
N = n1 + n2
Nvals = np.rint(np.diag(N).real).astype(int)
P = [0, 3, 1]  # |00>, |10>, |01>; physical labels are fixed
Ps = [np.diag((Nvals == n).astype(float)) for n in range(5)]
K = np.array([[1.3, -1], [-1, 1.3]], complex)
OMEGA = 5.0
US = [0.0, 0.25, 0.5, 1.0]
TIMES = [0.0, 0.125, 0.5, 1.0, 2.0]
DELTAS = [0.125, 0.03125, 0.0078125]
dark1 = np.zeros(9, complex); dark1[3] = dark1[1] = 1/np.sqrt(2)
dark2 = np.zeros(9, complex)
dark2[6] = dark2[2] = 0.5; dark2[4] = 1/np.sqrt(2)
X21 = np.outer(dark2, dark1.conj())
base_jumps = [np.sqrt(0.6)*b1, np.sqrt(0.6)*b2,
              np.sqrt(2.0)*(b1-b2)]
jumps = [(mu, n, L @ Ps[n]) for mu, L in enumerate(base_jumps)
         for n in range(1, 5)]
all_checks = []
cells = []

def check(name, error, group, **inputs):
    error = float(abs(error))
    item = dict(name=name, group=group, inputs=inputs, error=error,
                pass_=bool(error <= TOL))
    all_checks.append(item)
    if not item['pass_']:
        raise AssertionError(item)

def norm(X):
    return float(np.linalg.norm(X, 2))

def liouvillian(H, js):
    G = -1j*(np.kron(I, H)-np.kron(H.T, I))
    for L in js:
        M = L.conj().T @ L
        G += np.kron(L.conj(), L)-0.5*(np.kron(I, M)+np.kron(M.T, I))
    return G

def apply(S, X):
    return (S @ X.reshape(-1, order='F')).reshape(9, 9, order='F')

def matrix_unit(i, j):
    X = np.zeros((9, 9), complex)
    X[i, j] = 1
    return X

def low_formula(X, t):
    A = X[np.ix_(P, P)]
    T = expm(-(K+1j*OMEGA*np.eye(2))*t)
    R = T @ A[1:, 1:] @ T.conj().T
    B = np.zeros((3, 3), complex)
    B[0, 0] = A[0, 0] + np.trace(A[1:, 1:]) - np.trace(R)
    B[1:, 1:] = R
    B[1:, 0] = T @ A[1:, 0]
    B[0, 1:] = A[0, 1:] @ T.conj().T
    Y = np.zeros((9, 9), complex)
    Y[np.ix_(P, P)] = B
    return Y

def choi(S):
    C = np.zeros((81, 81), complex)
    for i in range(9):
        for j in range(9):
            C[i*9:(i+1)*9, j*9:(j+1)*9] = apply(S, matrix_unit(i, j))
    return C

def collision_channel(L, H, gap, delta):
    plus = np.array([[0, 0], [1, 0]], complex)
    HA = np.diag([0, gap])
    C = np.kron(L, plus) + np.kron(L.conj().T, plus.conj().T)
    Hfree = np.kron(H, np.eye(2)) + np.kron(I, HA)
    # Each of twelve collisions lasts delta/12; all free phases total delta.
    W = expm(-1j*(delta/12*Hfree + np.sqrt(delta)*C))
    Kraus = [W[a::2, 0::2] for a in range(2)]
    S = sum(np.kron(A.conj(), A) for A in Kraus)
    return S, C, Hfree, W

# Exact integer/rational controls are independent of floating-point fixtures.
pin_exact = Fraction(3,10)
diag_exact = 1+pin_exact
assert 2*diag_exact == pin_exact+(pin_exact+2)
assert diag_exact**2-1 == pin_exact*(pin_exact+2)
sector_sizes_exact = [sum(a+c == n for a in range(3) for c in range(3))
                      for n in range(5)]
exact = {
    'K_characteristic_roots': [str(Fraction(3, 10)), str(Fraction(23, 10))],
    'number_sector_dimensions': sector_sizes_exact,
    'number_conserving_completion_dimension': sum(d*d for d in [3, 2, 1]),
    'quartic_pair_completion_dimension': 3*3,
    'real_exchange_symmetric_quartic_dimension': 4,
    'addition_energy': [dict(u=str(Fraction(str(u))),
                             E0='0', E1='5',
                             E2=str(Fraction(10)+Fraction(str(u))),
                             chi=str(Fraction(str(u)))) for u in US],
}
assert exact['number_conserving_completion_dimension'] == 14
assert exact['quartic_pair_completion_dimension'] == 9
assert [int(np.trace(Q).real) for Q in Ps] == exact['number_sector_dimensions']
assert all(Fraction(row['E2'])-2*Fraction(row['E1']) == Fraction(row['chi'])
           for row in exact['addition_energy'])

Msum = sum(L.conj().T @ L for _, _, L in jumps)
check('loss_metric', norm(Msum-2*sum(K[i,j]*[b1,b2][i].conj().T @ [b1,b2][j]
                                   for i in range(2) for j in range(2))), 'structure')

# Quartic annihilators restrict to the identity coordinate map on N=2.
pair_ops = [b1@b1/np.sqrt(2), b1@b2, b2@b2/np.sqrt(2)]
two = [6, 4, 2]  # |20>, |11>, |02>
quartic_basis = []
for a in range(3):
    for c in range(3):
        Q = pair_ops[a].conj().T @ pair_ops[c]
        check('quartic_vanishes_low', norm(Q[:, P])+norm(Q[P, :]),
              'quartic', pair_a=a, pair_c=c)
        expected = np.zeros((3,3), complex); expected[a,c] = 1
        check('quartic_pair_coordinates', norm(Q[np.ix_(two,two)]-expected),
              'quartic', pair_a=a, pair_c=c)
        quartic_basis.append(Q.reshape(-1))
assert np.linalg.matrix_rank(np.array(quartic_basis), tol=1e-10) == 9
check('central_quartic_identity', norm(sum(A.conj().T@A for A in pair_ops)
                                      -0.5*(N@N-N)), 'quartic')

# Full low-sector superoperator, not selected state trajectories only.
G0 = None
collision_ref = {}
for u in US:
    H = OMEGA*N + u/2*(N@N-N)
    G = liouvillian(H, [L for _,_,L in jumps])
    if G0 is None: G0 = G
    # Unitary rotations of jump representations leave the physical generator
    # unchanged; unlike u, their angle is not an observable interaction dial.
    for theta in [0.2, 0.7]:
        R = np.array([[np.cos(theta),-np.sin(theta),0],
                      [np.sin(theta),np.cos(theta),0],[0,0,1]])
        rotated = []
        for sector in range(1,5):
            old = [L for _,n,L in jumps if n == sector]
            rotated.extend(sum(R[a,c]*old[c] for c in range(3)) for a in range(3))
        check('jump_representation_rotation', norm(liouvillian(H,rotated)-G),
              'representation',u=u,theta=theta)
    check('number_symmetry', norm(H@N-N@H), 'admissibility', u=u)
    check('positive_H', max(0, -np.linalg.eigvalsh(H).min()), 'admissibility', u=u)
    check('finite_energy_ceiling', max(0, norm(H)-26), 'admissibility', u=u)
    check('trace_preserving_generator', norm(np.eye(9).reshape(-1, order='F')@G),
          'admissibility', u=u)
    energy_derivative = 1j*(H@H-H@H)
    passive_rhs = np.zeros((9,9), complex)
    for mu,n,L in jumps:
        gap = OMEGA+u*(n-1)
        M = L.conj().T@L
        energy_derivative += L.conj().T@H@L-0.5*(M@H+H@M)
        passive_rhs -= gap*M
        check('lowering_energy', norm(H@L-L@H+gap*L), 'admissibility', u=u, mu=mu, n=n)
    check('energy_drain_identity', norm(energy_derivative-passive_rhs), 'admissibility', u=u)
    check('no_system_energy_gain', max(0,np.linalg.eigvalsh(energy_derivative).max()),
          'admissibility', u=u)
    for i in P:
        for j in P:
            X = matrix_unit(i,j)
            check('generator_identical_low', norm(apply(G-G0,X)), 'low_generator', u=u, i=i,j=j)
    for t in TIMES:
        S = expm(t*G)
        for i in P:
            for j in P:
                X = matrix_unit(i,j)
                check('whole_low_channel', norm(apply(S,X)-low_formula(X,t)),
                      'low_channel', u=u,t=t,i=i,j=j)
        if t == 0.5:
            C = choi(S)
            check('whole_channel_CP', max(0,-np.linalg.eigvalsh((C+C.conj().T)/2).min()),
                  'admissibility', u=u,t=t)
        # Two-to-one dark-mode coherence has an exactly shifted decay pole.
        X = X21
        pole = -3*0.3-1j*(OMEGA+u)
        expected = np.exp(pole*t)*X
        observed = apply(S,X)
        check('two_particle_pole', norm(observed-expected), 'spectroscopy', u=u,t=t)
        amplitude = np.trace(X.conj().T@observed)
        cells.append(dict(u=u,t=t,coherence_real=float(amplitude.real),
                          coherence_imag=float(amplitude.imag),
                          omega01=OMEGA,omega12=OMEGA+u,chi=u,
                          invariant_ratio_chi_over_kappa=u/0.3))
    # Full semigroup detects that the physical model has changed outside P.
    if u:
        difference = norm(apply(G-G0,X21))
        check('physical_difference_equals_u', abs(difference-u), 'spectroscopy', u=u)

    for delta in DELTAS:
        Sround = np.eye(81, dtype=complex)
        energy_commutators = []
        for mu,n,L in jumps:
            gap = OMEGA+u*(n-1)
            S,C,Hfree,W = collision_channel(L,H,gap,delta)
            check('collision_energy_commutator', norm(C@Hfree-Hfree@C),
                  'collision_energy', u=u,delta=delta,mu=mu,n=n)
            total_number = np.kron(N,np.eye(2))+np.kron(I,np.diag([0,1]))
            check('collision_number_commutator', norm(C@total_number-total_number@C),
                  'collision_energy', u=u,delta=delta,mu=mu,n=n)
            check('collision_unitary', norm(W.conj().T@W-np.eye(18)),
                  'collision_energy', u=u,delta=delta,mu=mu,n=n)
            Sround = S @ Sround
        if delta not in collision_ref: collision_ref[delta] = Sround
        for i in P:
            for j in P:
                check('finite_collision_low_equality',
                      norm(apply(Sround-collision_ref[delta],matrix_unit(i,j))),
                      'collision_low',u=u,delta=delta,i=i,j=j)
        cells.append(dict(u=u,delta=delta,
                          full_one_round_semigroup_error=norm(Sround-expm(delta*G)),
                          generator_approximation_error=norm((Sround-np.eye(81))/delta-G),
                          ancillas_per_round=12,round_duration=delta,
                          common_interaction_norm_ceiling=48/np.sqrt(delta)))

# Hypothesis relaxation controls: two-excitation access distinguishes models;
# a genuine full-field linear closure is spoiled by the supplied interaction.
field_controls = []
for u in US:
    H = OMEGA*N+u/2*(N@N-N)
    comm = 1j*(H@b1-b1@H)
    nonlinear = comm+1j*OMEGA*b1
    check('nonlinear_commutator', norm(nonlinear+1j*u*b1@(N-I)),
          'boundary', u=u)
    field_controls.append(dict(u=u,full_field_nonlinearity=norm(nonlinear),
                               low_sector_nonlinearity=norm(nonlinear[:,P])))

groups = {}
for x in all_checks:
    group = groups.setdefault(x['group'],dict(count=0,max_error=0.0))
    group['count'] += 1
    group['max_error'] = max(group['max_error'],x['error'])
out = dict(status='STANDARD_INTERACTION_COMPLETION_CONTROLS_PASS',
           scope='finite conventional null; no selector or Card 2',
           runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
           inputs=dict(omega=OMEGA,pin=0.3,edge_weight=1,local_levels=3,u=US,
                       times=TIMES,collision_deltas=DELTAS,low_basis=P),
           exact_controls=exact,groups=groups,spectroscopy_and_collision_cells=cells,
           field_boundary_controls=field_controls,checks=all_checks)
(ROOT/'results.json').write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
table = ['# Generated interaction-completion controls','',
         '| u | one-particle frequencies | two-particle addition shift chi | chi / pin |',
         '|---:|---:|---:|---:|']
for u in US: table.append(f'| {u:g} | 5 | {u:g} | {u/0.3:.8g} |')
table += ['', '| Check group | Cells | Maximum error |','|---|---:|---:|']
for group,x in groups.items(): table.append(f"| {group} | {x['count']} | {x['max_error']:.3e} |")
table += ['', 'Finite collision errors are convergence diagnostics, not exact semigroup identities.',
          'The report supplies a whole-domain analytic O(delta) limit bound and prices the diverging pulse ceiling.']
(ROOT/'TABLES.md').write_text('\n'.join(table)+'\n')
print(json.dumps(dict(status=out['status'],checks=len(all_checks),groups=groups),sort_keys=True))
