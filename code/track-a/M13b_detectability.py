"""M13b: does the long-history theorem (A2) hold through detectability, also when the interior alone is unstable? Prospectively specified.
Continuous time, 8 nodes: I = {0,1}, S = {2,3} (hears E, S, A), A_ = {4,5} (hears I, S, A), E = {6,7}; wiring cut and split blanket (M = 0);
unit diagonal noise. Observer: has read the blanket on [0, L], does not know the start (stationary start). Filtered covariance of
(x_I, x_E) given the blanket path: Riccati with F = diag(A_II, A_EE), Q = I, G = H' H, H = [A_BI, A_BE], started at the stationary
conditional covariance given x_B(0). Score: I(x_I(L); x_E(L) | blanket path on [0, L]) from the filtered covariance.
Lemma (proof in the paper): with a wiring cut and A stable, (F, H) is detectable (PBH), so the Riccati solution converges to the
stabilising solution, which is block-diagonal when M = 0; the cross block decays exponentially.
Predictions:
 (a) 'unstable interior' networks (A stable, A_II has an eigenvalue with real part > 0; 30 networks): MI(L = 40) < 1e-6 * MI(0) in all.
 (b) 'random' split networks (30): same.
 (c) 'hidden slow mode': interior node 0 drives nothing (unobservable through B), self-rate -0.02: forgetting is slow,
     MI(40)/MI(0) > 1e-3, and the log-slope of MI over L in [10, 40] lies in [-0.12, -0.02] (rate set by the hidden mode, ~2 x 0.02).
Seed 1313."""
import os, json
import numpy as np
from scipy.linalg import solve_continuous_lyapunov
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
I, S, Ac, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + Ac; n = 8
hears = {**{k: I + S + Ac for k in I + Ac}, **{k: E + S + Ac for k in S + E}}
rng = np.random.default_rng(1313)
ld = lambda X: np.linalg.slogdet(X)[1]
def gen():
    A = -np.eye(n)
    for i in range(n):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: A[i, j] = rng.uniform(0.2, 0.9) * rng.choice([-1, 1])
    return A
def stable(A, margin=0.05): return np.linalg.eigvals(A).real.max() < -margin
def mi_curve(A, Ls=(0, 5, 10, 20, 40)):
    IE = I + E; k = len(IE)
    Sg = solve_continuous_lyapunov(A, -np.eye(n))
    S0 = Sg[np.ix_(IE, IE)] - Sg[np.ix_(IE, B)] @ np.linalg.solve(Sg[np.ix_(B, B)], Sg[np.ix_(B, IE)])
    F = np.zeros((k, k)); F[:2, :2] = A[np.ix_(I, I)]; F[2:, 2:] = A[np.ix_(E, E)]
    H = np.hstack([A[np.ix_(B, I)], A[np.ix_(B, E)]]); G = H.T @ H
    f = lambda t, y: (lambda P: (F @ P + P @ F.T + np.eye(k) - P @ G @ P).ravel())(y.reshape(k, k))
    sol = solve_ivp(f, (0, max(Ls)), S0.ravel(), t_eval=list(Ls), rtol=1e-11, atol=1e-16, method='DOP853')
    out = []
    for j in range(len(Ls)):
        P = sol.y[:, j].reshape(k, k); P = (P + P.T) / 2
        out.append(max(0.5 * (ld(P[:2, :2]) + ld(P[2:, 2:]) - ld(P)), 0.0))
    return out
res = {'unstable_interior': [], 'random': [], 'hidden_slow': []}
while len(res['unstable_interior']) < 30 or len(res['random']) < 30:
    A = gen()
    if not stable(A): continue
    if np.linalg.eigvals(A[np.ix_(I, I)]).real.max() <= 0:
        A2 = A.copy(); A2[0, 0] = rng.uniform(0.1, 0.6)      # make the interior alone unstable
        if stable(A2) and np.linalg.eigvals(A2[np.ix_(I, I)]).real.max() > 0 and len(res['unstable_interior']) < 30:
            res['unstable_interior'].append(mi_curve(A2))
        if len(res['random']) < 30: res['random'].append(mi_curve(A))
while len(res['hidden_slow']) < 30:
    A = gen(); A[:, 0] = 0.0; A[0, 0] = -0.02                 # node 0 drives nothing: hidden from B, slow
    A[0, :] = 0.0; A[0, 0] = -0.02
    for j in [1] + S + Ac:                                    # node 0 is driven by the blanket and node 1, so it shares the start
        if rng.random() < 0.6: A[0, j] = rng.uniform(0.2, 0.9) * rng.choice([-1, 1])
    if stable(A, margin=0.01): res['hidden_slow'].append(mi_curve(A))
out = {}
Ls = [0, 5, 10, 20, 40]
for k, v in res.items():
    v = np.array(v); ratio = v[:, -1] / np.maximum(v[:, 0], 1e-300)
    sl = [np.polyfit(Ls[2:], np.log(np.maximum(r[2:], 1e-300)), 1)[0] for r in v]
    out[k] = dict(n=len(v), median_MI=np.median(v, 0).tolist(), max_ratio_40_0=float(ratio.max()), min_ratio_40_0=float(ratio.min()),
                  median_logslope_10_40=float(np.median(sl)), slope_range=[float(np.min(sl)), float(np.max(sl))])
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M13b.json'), 'w'), indent=1)
