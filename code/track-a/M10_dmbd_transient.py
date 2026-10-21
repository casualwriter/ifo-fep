"""M10: is DMBD-class residual dependence transient (start effect) or steady (boundary structure)? Prospectively specified
(written before running). Same model class as M9 (s=2, b=4, z=2, diagonal unit noise, no s-z blocks), spectral radius <= 0.5,
T in {8, 16, 32}, 100 models per case, seed 9696. Measure: I(s path; z path | b path) in nats.
Predictions: (G) general boundary, x0 known: conditional MI grows with T, MI(T=32)/MI(T=8) >= 3 (median of per-model ratio)
-> a steady per-step leak (structural). (S) split boundary, stationary start: MI(T=32)/MI(T=8) <= 1.5 -> a start effect (transient).
(K) split, x0 known: < 1e-10 for all T."""
import os, json
import numpy as np
from scipy.linalg import solve_discrete_lyapunov
HERE = os.path.dirname(os.path.abspath(__file__))
s, b1, b2, z = [0, 1], [2, 3], [4, 5], [6, 7]; b = b1 + b2; n = 8
rng = np.random.default_rng(9696)
def gen(split):
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if {i, j} & set(s) and {i, j} & set(z): continue
            if split and ((i in b1 and j in z) or (i in b2 and j in s)): continue
            if rng.random() < 0.6: A[i, j] = rng.uniform(0.1, 0.4) * rng.choice([-1, 1])
    r = np.abs(np.linalg.eigvals(A)).max(); return A * (0.5 / r) if r > 0.5 else A
ld = lambda X: np.linalg.slogdet(X)[1]
def cmi(A, stationary, T):
    V = solve_discrete_lyapunov(A, np.eye(n)) if stationary else np.zeros((n, n))
    Vs = [V]
    for t in range(T): Vs.append(A @ Vs[-1] @ A.T + np.eye(n))
    P = [np.eye(n)]
    for t in range(T): P.append(A @ P[-1])
    C = np.zeros(((T + 1) * n, (T + 1) * n))
    for t in range(T + 1):
        for u in range(t + 1):
            blk = P[t - u] @ Vs[u]; C[t*n:(t+1)*n, u*n:(u+1)*n] = blk; C[u*n:(u+1)*n, t*n:(t+1)*n] = blk.T
    t0 = 0 if stationary else 1
    idx = lambda g: [t * n + k for t in range(t0, T + 1) for k in g]
    Sp, Zp, Bp = idx(s), idx(z), idx(b); L = lambda X: ld(C[np.ix_(X, X)])
    return 0.5 * (L(Sp + Bp) + L(Zp + Bp) - L(Bp) - L(Sp + Zp + Bp))
out = {}
for key, split, stat in [('G_general_known', False, False), ('S_split_stationary', True, True), ('K_split_known', True, False)]:
    vals = []
    for _ in range(100):
        A = gen(split); vals.append([cmi(A, stat, T) for T in (8, 16, 32)])
    v = np.array(vals)
    out[key] = dict(median_T8_16_32=np.median(v, 0).tolist(), ratio32_8_median=float(np.median(v[:, 2] / np.maximum(v[:, 0], 1e-300))),
                    max_abs=float(np.abs(v).max()))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M10.json'), 'w'), indent=1)
