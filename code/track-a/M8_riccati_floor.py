"""M8: exact floor for an observer who knows x(0) and the whole B path, via the Kalman-Bucy Riccati equation
(plan/20261016-k1-next-steps-v1.md section 2). Predictions fixed before running:
R1 M = 0 (Loop 3 split nets, eps = 0): Riccati floor < 1e-10 (median and max).
R2 eps = 0.1: Riccati floor matches the discrete computation at m = 32 within 5% (median |ratio - 1|).
R3 Riccati floor / eps^2 at eps = 0.01 within 5% of that at eps = 0.001 (median of per-net ratio).
State xi = (u, x_I, x_E); F = diag(0, A_II, A_EE); Q = diag(0, D_I, D_E); H = [0, A_BI, A_BE]; R = D_B.
Sigma' = F S + S F' + Q - S H' R^-1 H S. Nets: Loop 3 generator (seed 9191, 200 nets), A[4,6] = eps, push 0.1, h = 2."""
import os, sys, json
import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M3_loop2.py"), encoding="utf-8").read().split("\nout = {}")[0])
I, S, A_, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + A_
hears = {**{k: I + S + A_ for k in I + A_}, **{k: E + S + A_ for k in S + E}}
rng = np.random.default_rng(9191); su = 0.1; H_ = 2.0
def gen():
    M = -np.eye(8)
    for i in range(8):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: M[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
    M[4, 0] = 0.5; return M
def stab(M):
    m = np.linalg.eigvals(M).real.max(); return M - (m + 0.1) * np.eye(8) if m > -0.1 else M
ld = lambda X: np.linalg.slogdet(X)[1]
def riccati_ratio(A, src, tgt, D=None):
    D = np.eye(8) if D is None else D; p, q = len(src), len(tgt); k = 2 * p + q
    F = np.zeros((k, k)); F[p:2 * p, p:2 * p] = A[np.ix_(src, src)]; F[2 * p:, 2 * p:] = A[np.ix_(tgt, tgt)]
    Q = np.zeros((k, k)); Q[p:2 * p, p:2 * p] = D[np.ix_(src, src)]; Q[2 * p:, 2 * p:] = D[np.ix_(tgt, tgt)]
    Hm = np.hstack([np.zeros((len(B), p)), A[np.ix_(B, src)], A[np.ix_(B, tgt)]]); G = Hm.T @ np.linalg.solve(D[np.ix_(B, B)], Hm)
    S0 = np.zeros((k, k)); S0[:p, :p] = np.eye(p); S0[:p, p:2 * p] = su * np.eye(p); S0[p:2 * p, :p] = su * np.eye(p); S0[p:2 * p, p:2 * p] = su ** 2 * np.eye(p)
    f = lambda t, y: (lambda S_: (F @ S_ + S_ @ F.T + Q - S_ @ G @ S_).ravel())(y.reshape(k, k))
    Sh = solve_ivp(f, (0, H_), S0.ravel(), rtol=1e-11, atol=1e-14, method='DOP853').y[:, -1].reshape(k, k); Sh = (Sh + Sh.T) / 2
    Suu, SuE, SEE = Sh[:p, :p], Sh[:p, 2 * p:], Sh[2 * p:, 2 * p:]
    num = 0.5 * (ld(Suu) - ld(Suu - SuE @ np.linalg.solve(SEE, SuE.T)))
    Sig = solve_continuous_lyapunov(A, -D); Ph = expm(A * H_); W = Sig - Ph @ Sig @ Ph.T
    Sg = Ph[np.ix_(tgt, src)] * su; Wt = W[np.ix_(tgt, tgt)]
    den = 0.5 * (ld(Sg @ Sg.T + Wt) - ld(Wt))
    return num / den
def riccati_floor(A, D=None): return max(riccati_ratio(A, I, E, D), riccati_ratio(A, E, I, D))
def discrete_floor(A, m=32):
    o = Obs(A, H_, m); P = su * np.eye(8); return np.nanmax([o.lc(P[:, I], I, B, E), o.lc(P[:, E], E, B, I)])
nets = [gen() for _ in range(200)]; r1, r2, r3 = [], [], []
for M0 in nets:
    A0 = stab(M0.copy()); r1.append(riccati_floor(A0))
    M = M0.copy(); M[4, 6] = 0.1; A1 = stab(M); r2.append(riccati_floor(A1) / discrete_floor(A1))
    a, b = M0.copy(), M0.copy(); a[4, 6] = 0.01; b[4, 6] = 0.001
    r3.append((riccati_floor(stab(a)) / 0.01 ** 2) / (riccati_floor(stab(b)) / 0.001 ** 2))
out = dict(R1_median=float(np.median(np.abs(r1))), R1_max=float(np.max(np.abs(r1))),
           R2_median_abs_dev=float(np.median(np.abs(np.array(r2) - 1))), R2_q90_abs_dev=float(np.quantile(np.abs(np.array(r2) - 1), 0.9)),
           R3_median_abs_dev=float(np.median(np.abs(np.array(r3) - 1))))
print(json.dumps(out, indent=1)); json.dump(dict(summary=out, r2=r2, r3=r3), open(os.path.join(HERE, 'M8.json'), 'w'), indent=1)
