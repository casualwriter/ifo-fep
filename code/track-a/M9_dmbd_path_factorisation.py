"""M9: does the path-level blanket factorisation p(s, z | b) = p(s | b) p(z | b) (Beck & Ramstead 2025, Eq. 7) hold in the
DMBD model class? Discrete-time linear Gaussian x_{t+1} = A x_t + w_t, w ~ N(0, I), with the s-z blocks of A zero (as in DMBD),
boundary dynamics depending on s, b and z. Exact Gaussian: condition on the whole b path (t = 0..T) and measure the conditional
dependence between the s path and the z path: I(s path; z path | b path) (nats). Predictions fixed before running:
D1 general DMBD-type A (b hears s and z), x_0 fixed (known): I > 1e-3 in >= 95% of 200 models.
D2 split (each b node hears s or z, not both), x_0 fixed: I < 1e-10 in all models.
D3 split, x_0 drawn from the stationary law (not known): I > 1e-6 in >= 90% of models (second collider route via the start).
Sizes: s = 2, b = 4 (b1 = {2,3} hears s, b2 = {4,5} hears z in the split), z = 2; T = 8; weights U(0.1, 0.4) +-, prob 0.6;
spectral radius scaled to <= 0.9. Seed 9595."""
import os, json
import numpy as np
from scipy.linalg import solve_discrete_lyapunov
HERE = os.path.dirname(os.path.abspath(__file__))
s, b1, b2, z = [0, 1], [2, 3], [4, 5], [6, 7]; b = b1 + b2; n, T = 8, 8
rng = np.random.default_rng(9595)
def gen(split):
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if {i, j} & set(s) and {i, j} & set(z): continue                 # no direct s-z coupling (both directions)
            if split and ((i in b1 and j in z) or (i in b2 and j in s)): continue   # b1 hears s only, b2 hears z only
            if rng.random() < 0.6: A[i, j] = rng.uniform(0.1, 0.4) * rng.choice([-1, 1])
    r = np.abs(np.linalg.eigvals(A)).max(); return A * (0.9 / r) if r > 0.9 else A
ld = lambda X: np.linalg.slogdet(X)[1]
def cmi(A, stationary):
    S0 = solve_discrete_lyapunov(A, np.eye(n)) if stationary else np.zeros((n, n))
    Phi = [np.linalg.matrix_power(A, t) for t in range(T + 1)]
    C = np.zeros(((T + 1) * n, (T + 1) * n))
    for t in range(T + 1):
        for u in range(t + 1):
            V = Phi[u] @ S0 @ Phi[u].T + sum(Phi[k] @ Phi[k].T for k in range(u))     # Cov(x_u)
            blk = Phi[t - u] @ V; C[t * n:(t + 1) * n, u * n:(u + 1) * n] = blk; C[u * n:(u + 1) * n, t * n:(t + 1) * n] = blk.T
    idx = lambda g, t0: [t * n + k for t in range(t0, T + 1) for k in g]
    t0 = 0 if stationary else 1                                              # with x_0 known, x_0 carries no randomness
    Sp, Zp, Bp = idx(s, t0), idx(z, t0), idx(b, t0)
    L = lambda X: ld(C[np.ix_(X, X)])
    return 0.5 * (L(Sp + Bp) + L(Zp + Bp) - L(Bp) - L(Sp + Zp + Bp))
out = {}
g1 = [cmi(gen(False), False) for _ in range(200)]; g2 = [cmi(gen(True), False) for _ in range(200)]; g3 = [cmi(gen(True), True) for _ in range(200)]
out['D1_general_x0known'] = dict(frac_gt_1e3=float(np.mean(np.array(g1) > 1e-3)), median=float(np.median(g1)))
out['D2_split_x0known'] = dict(max_abs=float(np.max(np.abs(g2))), median=float(np.median(g2)))
out['D3_split_stationary'] = dict(frac_gt_1e6=float(np.mean(np.array(g3) > 1e-6)), median=float(np.median(g3)))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M9.json'), 'w'), indent=1)
