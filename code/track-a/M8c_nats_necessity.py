"""M8c: (1) the floor in nats (numerator, denominator) for M8's networks at eps = 0.1 and eps = 0, and the number of networks
whose denominator is ~0 (< 1e-12); (2) algebra check of the no-push necessity term: x(0) known, no push, Riccati on (x_I, x_E)
given the blanket path: Sigma_IE(t) = -(t^3/3) D_I M D_E + O(t^4). Prospectively specified prediction for (2): with D = I,
median relative error of Sigma_IE(t) against -(t^3/3) M at t = 0.01 below 1%."""
import os, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, solve_continuous_lyapunov
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M8_riccati_floor.py"), encoding="utf-8").read().split("\nnets = [gen()")[0])
def nats(A, src, tgt):
    p, q = len(src), len(tgt); k = 2 * p + q; D = np.eye(8)
    F = np.zeros((k, k)); F[p:2*p, p:2*p] = A[np.ix_(src, src)]; F[2*p:, 2*p:] = A[np.ix_(tgt, tgt)]
    Q = np.zeros((k, k)); Q[p:2*p, p:2*p] = np.eye(p); Q[2*p:, 2*p:] = np.eye(q)
    Hm = np.hstack([np.zeros((len(B), p)), A[np.ix_(B, src)], A[np.ix_(B, tgt)]]); G = Hm.T @ Hm
    S0 = np.zeros((k, k)); S0[:p, :p] = np.eye(p); S0[:p, p:2*p] = S0[p:2*p, :p] = su * np.eye(p); S0[p:2*p, p:2*p] = su ** 2 * np.eye(p)
    f = lambda t, y: (lambda S_: (F @ S_ + S_ @ F.T + Q - S_ @ G @ S_).ravel())(y.reshape(k, k))
    Sh = solve_ivp(f, (0, H_), S0.ravel(), rtol=1e-11, atol=1e-14, method='DOP853').y[:, -1].reshape(k, k); Sh = (Sh + Sh.T) / 2
    Suu, SuE, SEE = Sh[:p, :p], Sh[:p, 2*p:], Sh[2*p:, 2*p:]
    num = 0.5 * (ld(Suu) - ld(Suu - SuE @ np.linalg.solve(SEE, SuE.T)))
    Sig = solve_continuous_lyapunov(A, -D); Ph = expm(A * H_); W = Sig - Ph @ Sig @ Ph.T
    Sg = Ph[np.ix_(tgt, src)] * su; Wt = W[np.ix_(tgt, tgt)]
    return num, 0.5 * (ld(Sg @ Sg.T + Wt) - ld(Wt))
def nopush_cross(A, t):
    p, q = len(I), len(E); k = p + q
    F = np.zeros((k, k)); F[:p, :p] = A[np.ix_(I, I)]; F[p:, p:] = A[np.ix_(E, E)]
    Hm = np.hstack([A[np.ix_(B, I)], A[np.ix_(B, E)]]); G = Hm.T @ Hm
    f = lambda tt, y: (lambda S_: (F @ S_ + S_ @ F.T + np.eye(k) - S_ @ G @ S_).ravel())(y.reshape(k, k))
    St = solve_ivp(f, (0, t), np.zeros(k * k), rtol=1e-12, atol=1e-20, method='DOP853').y[:, -1].reshape(k, k)
    M = A[np.ix_(B, I)].T @ A[np.ix_(B, E)]
    return St[:p, p:], -(t ** 3 / 3) * M
nets = [gen() for _ in range(200)]; r = {0.1: [], 0.0: []}; err = []
for M0 in nets:
    for eps in r:
        M = M0.copy(); M[4, 6] = eps; A = stab(M)
        r[eps].append([nats(A, I, E), nats(A, E, I)])
    M = M0.copy(); M[4, 6] = 0.1; A = stab(M)
    got, pred = nopush_cross(A, 0.01)
    if np.abs(pred).max() > 0: err.append(float(np.linalg.norm(got - pred) / np.linalg.norm(pred)))
out = {}
for eps, v in r.items():
    v = np.array(v)                                   # nets x dir x (num, den)
    out[f'eps_{eps}'] = dict(num_nats_median_IE_EI=np.median(v[:, :, 0], 0).tolist(), den_nats_median_IE_EI=np.median(v[:, :, 1], 0).tolist(),
                             n_den_lt_1e12=int((v[:, :, 1] < 1e-12).sum()))
out['necessity_t3_rel_err_median'] = float(np.median(err)); out['necessity_t3_rel_err_q90'] = float(np.quantile(err, 0.9))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M8c.json'), 'w'), indent=1)
