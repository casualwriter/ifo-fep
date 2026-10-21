"""M13: does a long blanket history before time 0 remove the start condition? (Theorem A, long-lived particle). Prospectively specified.
Discrete-time stationary process, same model class and generator as M10 (s=2, b=4 split as b1 hears s / b2 hears z, z=2, unit diagonal
noise, spectral radius <= 0.5, seed 9696), 100 models per case. The observer conditions on the blanket on [-L, T] and asks about the
s and z paths on [0, T], T = 8. Measure: I(s_{0:T}; z_{0:T} | b_{-L:T}) in nats, L in {0, 2, 4, 8, 16}.
Predictions: split boundary (M = 0): median MI falls with L and is < 1e-8 at L = 16 (start effect forgotten).
General boundary (M != 0): median MI at L = 16 stays >= 0.5 x its value at L = 0 (the M leak does not fade)."""
import os, json
import numpy as np
from scipy.linalg import solve_discrete_lyapunov
HERE = os.path.dirname(os.path.abspath(__file__))
s, b1, b2, z = [0, 1], [2, 3], [4, 5], [6, 7]; b = b1 + b2; n, T = 8, 8
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
def mi(A, L):
    Sg = solve_discrete_lyapunov(A, np.eye(n)); N = L + T + 1        # times -L..T -> index 0..N-1; stationary
    P = [np.eye(n)]
    for _ in range(N): P.append(A @ P[-1])
    C = np.zeros((N * n, N * n))
    for t in range(N):
        for u in range(t + 1):
            blk = P[t - u] @ Sg; C[t*n:(t+1)*n, u*n:(u+1)*n] = blk; C[u*n:(u+1)*n, t*n:(t+1)*n] = blk.T
    idx = lambda g, t0: [t * n + k for t in range(t0, N) for k in g]
    Sp, Zp, Bp = idx(s, L), idx(z, L), idx(b, 0); Lg = lambda X: ld(C[np.ix_(X, X)])
    return 0.5 * (Lg(Sp + Bp) + Lg(Zp + Bp) - Lg(Bp) - Lg(Sp + Zp + Bp))
out = {}
for key, split in [('split_M0', True), ('general_Mne0', False)]:
    v = np.array([[mi(A, L) for L in (0, 2, 4, 8, 16)] for A in [gen(split) for _ in range(100)]])
    out[key] = dict(L=[0, 2, 4, 8, 16], median=np.median(v, 0).tolist(), q90=np.quantile(v, 0.9, 0).tolist())
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M13.json'), 'w'), indent=1)
