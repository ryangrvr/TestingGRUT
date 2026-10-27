"""Twin (Lyapunov) test: start from own tau=100 snapshot with a random relative perturbation eps of each field's
fluctuation RMS (all lattice modes, q only), write a resumable checkpoint for indep3d.py."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import indep3d as I
src, out, eps = sys.argv[1], sys.argv[2], float(sys.argv[3])
am, C = 8.869218738806, 1.9857962750634526
z = np.load(src); q, p, s, tau = z['q'].copy(), z['p'].copy(), float(z['s']), float(z['tau'])
N = q.shape[-1]; L = 12.0
rng = np.random.default_rng(20261027)
for c in range(3):
    rms = np.std(q[c]); q[c] += eps * rms * rng.normal(size=q[c].shape)
os.makedirs(out, exist_ok=True)
lat = I.Lattice(N, L)
rec = [I.diagnostics(lat, q, p, s, tau, am + C * s, 0.0)]
np.savez(os.path.join(out, 'checkpoint.npz'), q=q, p=p, s=s, tau=tau, N=N, L=L)
np.save(os.path.join(out, 'records.npy'), np.array(rec))
print('perturbed start at tau', tau, 's', s, 'eps', eps)
