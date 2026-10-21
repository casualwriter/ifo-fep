"""M17: does the internal observer need the future of the blanket? Causal (filtering) reader. Prospectively specified.
Answers the Vector review of fep-blanket v0.6 (inbox/20261019-fep-blanket-v0.6-vector-review.md): "conditioning on the whole
blanket path, including the future, is not something the particle can do".
Setup as M13b: continuous time, 8 nodes, I = {0,1}, S = {2,3}, A_ = {4,5}, E = {6,7}, wiring cut (A_IE = A_EI = 0), unit diagonal noise.
Causal reader: has read the blanket on (-inf, t] only (no future), stationary. Its uncertainty about (x_I(t), x_E(t)) is the
stabilising solution P of the filter Riccati  F P + P F' + Q - P G P = 0,  F = diag(A_II, A_EE), Q = I, G = H' D_B^-1 H, H = [A_BI, A_BE].
Lemma (proof, 3 lines): the I-E block of G is M = A_BI' D_B^-1 A_BE. If P_IE = 0, the I-E block of the Riccati reads -P_I M P_E = 0,
so M = 0 (P_I, P_E > 0). If M = 0, the Riccati splits into two blocks, and the unique stabilising P is block-diagonal.
So, with a wiring cut: current-state screening for a causal reader  <=>  M = 0. The same condition as for paths; no future needed.
Score (a): I(x_I(t); x_E(t) | blanket on (-inf, t]) from P.
Predictions (a):
 (a1) split networks (M = 0; 30): MI < 1e-10 in all.
 (a2) collider networks (one blanket node hears I and E, so M != 0; 30): MI > 1e-4 in all.
 (a3) max over all 60 of |(Riccati residual) I-E block + P_I M P_E| at P_IE set to 0 is < 1e-9 (the identity used in the proof).
Exploratory (b), no threshold, about a different question (does x_I(t) itself encode the reader's estimate of x_E(t)? The FEP's
"synchronisation" claim; Vector's suggested test 1): R2 = share of the variance of E[x_E(t) | blanket past] that is linear in x_I(t).
Direction expected: R2 < 1 in all networks (screening does not imply encoding).
Seed 1717."""
import os, json
import numpy as np
from scipy.linalg import solve_continuous_are, solve_continuous_lyapunov
HERE = os.path.dirname(os.path.abspath(__file__))
I, S, Ac, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + Ac; IE = I + E; n = 8
split = {**{k: I + S + Ac for k in I + Ac}, **{k: E + S + Ac for k in S + E}}   # S hears E, A_ hears I: M = 0
rng = np.random.default_rng(1717)
ld = lambda X: np.linalg.slogdet(X)[1]
ix = np.ix_

def gen(collider):
    hears = {k: list(v) for k, v in split.items()}
    if collider: hears[2] = E + S + Ac + [0]          # sensory node 2 also hears interior node 0
    A = -np.eye(n)
    for i in range(n):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: A[i, j] = rng.uniform(0.2, 0.9) * rng.choice([-1, 1])
    if collider: A[2, 0] = rng.uniform(0.3, 0.9) * rng.choice([-1, 1]); A[2, 6] = rng.uniform(0.3, 0.9) * rng.choice([-1, 1])
    return A

def analyse(A):
    k = 4
    F = np.zeros((k, k)); F[:2, :2] = A[ix(I, I)]; F[2:, 2:] = A[ix(E, E)]
    H = np.hstack([A[ix(B, I)], A[ix(B, E)]]); G = H.T @ H            # D_B = I
    M = A[ix(B, I)].T @ A[ix(B, E)]
    P = solve_continuous_are(F.T, H.T, np.eye(k), np.eye(len(B)))     # filter Riccati (dual form)
    P = (P + P.T) / 2
    mi = max(0.5 * (ld(P[:2, :2]) + ld(P[2:, 2:]) - ld(P)), 0.0)
    P0 = P.copy(); P0[:2, 2:] = 0; P0[2:, :2] = 0                      # identity check at P_IE = 0
    R = F @ P0 + P0 @ F.T + np.eye(k) - P0 @ G @ P0
    ident = np.abs(R[:2, 2:] + P0[:2, :2] @ M @ P0[2:, 2:]).max()
    # (b) augmented system z = (x, xhat), xhat = estimate of (x_I, x_E) given the blanket past
    K = P @ H.T                                                        # gain, D_B = I
    AB_IE = np.vstack([A[ix(I, B)], A[ix(E, B)]]); ABB = A[ix(B, B)]
    Z = np.zeros((n + k, n + k)); Z[:n, :n] = A
    # d xhat = (F xhat + A_(IE)B x_B) dt + K (dx_B - (A_BB x_B + H xhat) dt), dx_B = A_B. x dt + dW_B
    Z[n:, n:] = F - K @ H
    Z[ix(range(n, n + k), B)] += AB_IE
    Z[n:, :n] += K @ A[B, :]
    Z[ix(range(n, n + k), B)] -= K @ ABB
    N = np.zeros((n + k, n)); N[:n] = np.eye(n); N[n:, B] = K            # noise loading
    C = solve_continuous_lyapunov(Z, -N @ N.T)
    xI = I; mE = [n + 2, n + 3]
    Cm = C[ix(mE, mE)]; Cmi = C[ix(mE, xI)]; Cii = C[ix(xI, xI)]
    r2 = float(np.trace(Cmi @ np.linalg.solve(Cii, Cmi.T)) / np.trace(Cm))
    # sanity: estimation error covariance from the augmented system equals P
    D = np.zeros((k, n + k)); D[:, IE] = np.eye(k); D[:, n:] = -np.eye(k)
    perr = np.abs(D @ C @ D.T - P).max()
    return mi, float(ident), r2, float(perr), float(np.abs(M).max())

res = {'split': [], 'collider': []}
while len(res['split']) < 30 or len(res['collider']) < 30:
    for key in res:
        if len(res[key]) >= 30: continue
        A = gen(key == 'collider')
        if np.linalg.eigvals(A).real.max() < -0.05: res[key].append(analyse(A))
out = {}
for key, v in res.items():
    v = np.array(v)
    out[key] = dict(n=len(v), MI_min=float(v[:, 0].min()), MI_median=float(np.median(v[:, 0])), MI_max=float(v[:, 0].max()),
                    max_M=float(v[:, 4].max()), R2_min=float(v[:, 2].min()), R2_median=float(np.median(v[:, 2])),
                    R2_max=float(v[:, 2].max()), max_identity_err=float(v[:, 1].max()), max_P_check_err=float(v[:, 3].max()))
out['pred_a1_split_MI_lt_1e-10'] = out['split']['MI_max'] < 1e-10
out['pred_a2_collider_MI_gt_1e-4'] = out['collider']['MI_min'] > 1e-4
out['pred_a3_identity_lt_1e-9'] = max(out['split']['max_identity_err'], out['collider']['max_identity_err']) < 1e-9
out['expl_b_R2_lt_1_all'] = max(out['split']['R2_max'], out['collider']['R2_max']) < 1 - 1e-9
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M17.json'), 'w'), indent=1)
