"""M2b: the noise notation of the 'snapshot vs push' section (reviewers: if A = -(D + sQ) Pi, the stationary precision is 2 Pi only
when the diffusion matrix equals the same D). Prospectively specified. 200 random cases, 6 nodes, seed 2222:
Pi SPD with Pi_IE = 0 (I = {0,1}, B = {2,3}, E = {4,5}); Q antisymmetric; D diagonal U(0.5, 2); s = 0.3.
Predictions: (i) diffusion = D: stationary precision equals 2 Pi (max relative error < 1e-9), so its I-E block is 0;
(ii) diffusion = I (while D != I in the drift): relative error > 1e-2 in >= 95% of cases, and the precision's I-E block is non-zero
(> 1e-3 relative) in >= 95%."""
import os, json
import numpy as np
from scipy.linalg import solve_continuous_lyapunov
HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(2222); I, B, E = [0, 1], [2, 3], [4, 5]
e1, e2, z2 = [], [], []
for _ in range(200):
    M = rng.uniform(-1, 1, (6, 6)) * (rng.random((6, 6)) < 0.6); M = (M + M.T) / 2
    M[np.ix_(I, E)] = 0; M[np.ix_(E, I)] = 0; np.fill_diagonal(M, 0)
    Pi = M + (0.5 - np.linalg.eigvalsh(M).min()) * np.eye(6)
    Q = rng.standard_normal((6, 6)); Q = Q - Q.T; D = np.diag(rng.uniform(0.5, 2, 6)); A = -(D + 0.3 * Q) @ Pi
    P1 = np.linalg.inv(solve_continuous_lyapunov(A, -D)); P2 = np.linalg.inv(solve_continuous_lyapunov(A, -np.eye(6)))
    e1.append(np.abs(P1 - 2 * Pi).max() / np.abs(2 * Pi).max()); e2.append(np.abs(P2 - 2 * Pi).max() / np.abs(2 * Pi).max())
    z2.append(np.abs(P2[np.ix_(I, E)]).max() / np.abs(P2).max())
out = dict(diffusion_D_max_relerr=float(max(e1)), diffusion_I_frac_relerr_gt=float(np.mean(np.array(e2) > 1e-2)),
           diffusion_I_frac_IE_nonzero=float(np.mean(np.array(z2) > 1e-3)), diffusion_I_median_IE_rel=float(np.median(z2)))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M2b.json'), 'w'), indent=1)
