"""M16: can direct inside-outside couplings cancel in the path likelihood? (Step 4, the endpoint-term argument). Prospectively specified.
Euler discretisation of dx = A x dt + dW (dt = 0.02, 50 steps, h = 1), start known (x0 = 0). 8 nodes: I = {0,1}, S = {2,3}, A_ = {4,5},
E = {6,7}, split blanket (M = 0). Exact Gaussian path MI: I(I path; E path | B path). 40 networks per case, seed 1616.
Cases: (a) wiring cut (no I-E edges): MI < 1e-10 in all.
       (b) 'symmetric' direct couplings A_EI = c R, A_IE = c R' (with D = I this gives G = C', the cancellation the reviewers raised;
           the integral term vanishes, the endpoint term does not), c = 0.2: MI > 1e-4 in all.
       (c) one-way direct coupling A_EI = c R only, c = 0.2: MI > 1e-4 in all (sanity)."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
I, S, Ac, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + Ac; n = 8
hears = {**{k: I + S + Ac for k in I + Ac}, **{k: E + S + Ac for k in S + E}}
rng = np.random.default_rng(1616); dt, K = 0.02, 50
ld = lambda X: np.linalg.slogdet(X)[1]
def gen():
    A = -np.eye(n)
    for i in range(n):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: A[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
    return A
def path_mi(A):
    Phi = np.eye(n) + A * dt; Qd = dt * np.eye(n)
    P = [np.eye(n)]
    for _ in range(K): P.append(Phi @ P[-1])
    V = [np.zeros((n, n))]
    for _ in range(K): V.append(Phi @ V[-1] @ Phi.T + Qd)
    C = np.zeros((K * n, K * n))
    for a in range(1, K + 1):
        for b in range(1, a + 1):
            blk = P[a - b] @ V[b]; C[(a-1)*n:a*n, (b-1)*n:b*n] = blk; C[(b-1)*n:b*n, (a-1)*n:a*n] = blk.T
    idx = lambda g: [(t - 1) * n + k for t in range(1, K + 1) for k in g]
    Ip, Ep, Bp = idx(I), idx(E), idx(B); L = lambda X: ld(C[np.ix_(X, X)])
    return 0.5 * (L(Ip + Bp) + L(Ep + Bp) - L(Bp) - L(Ip + Ep + Bp))
out = {}
for case in ['a_wiring_cut', 'b_symmetric', 'c_one_way']:
    v = []
    while len(v) < 40:
        A = gen()
        if case != 'a_wiring_cut':
            R = rng.uniform(0.2, 0.8, (2, 2)) * rng.choice([-1, 1], (2, 2))
            A[np.ix_(E, I)] = 0.2 * R
            if case == 'b_symmetric': A[np.ix_(I, E)] = 0.2 * R.T
        if np.linalg.eigvals(A).real.max() < -0.05: v.append(path_mi(A))
    v = np.array(v); out[case] = dict(min=float(v.min()), median=float(np.median(v)), max=float(v.max()))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M16.json'), 'w'), indent=1)
