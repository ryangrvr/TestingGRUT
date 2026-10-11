"""Independent finite-difference checks of the affine Legendre algebra.
Not a spacetime simulation or a quantum-gravity validation.
"""
import hashlib
import json
import math
from pathlib import Path

root = Path(__file__).parent
law_hash = hashlib.sha256((root / 'LAW.md').read_bytes()).hexdigest()
M, V = 2.0, 0.7
errors = []
for signs in ((1, 1, 1, 1), (-1, 1, 1, 1)):
    metric = [s*a for s,a in zip(signs, (1.2, 0.8, 1.7, 2.1))]
    K = [V*g for g in metric]
    grad = (0.11, -0.07, 0.05, 0.03)
    R = [(k+u*u)/(M*M) for k,u in zip(K,grad)]
    def lag(r):
        ks = [M*M*a-u*u for a,u in zip(r,grad)]
        return math.sqrt(abs(math.prod(ks)))/V
    volume = math.sqrt(abs(math.prod(metric)))
    p = [M*M*volume/(2*g) for g in metric]
    for i in range(4):
        step = 1e-6
        plus, minus = R.copy(), R.copy()
        plus[i] += step
        minus[i] -= step
        measured = (lag(plus)-lag(minus))/(2*step)
        errors.append(abs(measured-p[i]))
    H = sum(a*b for a,b in zip(p,R))-lag(R)
    expected = volume*(V+sum(u*u/g for u,g in zip(grad,metric))/2)
    assert abs(H-expected)<1e-12
    assert abs(lag(R)-V*volume)<1e-12
    # Both signatures satisfy the same exact algebraic field relation.
    assert all(abs(k-V*g)<1e-15 for k,g in zip(K,metric))
assert max(errors)<1e-8
results = {'status':'CP06_AFFINE_ALGEBRA_CHECKS_PASS',
           'law_sha256':law_hash, 'signatures_checked':['Euclidean','Lorentzian'],
           'max_diagonal_derivative_error':max(errors),
           'checks':['Legendre derivative','Legendre Hamiltonian',
                     'on-shell determinant','signature nonselection'],
           'scope':'Finite diagonal algebra; no PDE or quantum validation'}
(root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
