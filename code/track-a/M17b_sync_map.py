"""M17b: a fair test of the FEP's encoding claim (synchronisation map). Prospectively specified.
M17 (b) asked whether the internal STATE x_I(t) carries the reader's estimate of x_E(t). The FEP claims less: the internal
EXPECTATION given the blanket, eta = E[x_I(t) | blanket], determines the external expectation mu = E[x_E(t) | blanket]
through a map mu = sigma(eta) (the synchronisation map). In linear-Gaussian systems both are linear in the blanket data, so a map
exists iff mu is linear in eta, i.e. R2(mu on eta) = 1 (share of the variance of mu explained linearly by eta).
Two readings of "given the blanket":
 (S) snapshot: the current blanket state x_B(t) only (stationary density);
 (H) history: the past blanket path (-inf, t] (causal reader of M17; Kalman filter).
Networks: as M17 (8 nodes, I = {0,1}, S = {2,3}, A_ = {4,5}, E = {6,7}, wiring cut, unit noise), split (M = 0) and collider (M != 0), 30 each.
Positive control (code check): 6 nodes, I = {0,1}, B = {2,3} (one sensory, one active), E = {4,5}: under (S), eta and mu are both
linear maps of the same 2-dim x_B, so if eta's map is invertible the synchronisation map exists and R2 = 1.
Predictions:
 (p1) control, (S): R2 > 1 - 1e-9 in all 30.
 (p2) 8-node networks, (S): R2 < 0.999 in >= 90% (split and collider separately).
 (p3) 8-node networks, (H): R2 < 0.999 in >= 90% (split and collider separately).
Death clause: if p2 or p3 fails, the synchronisation map generically exists in these networks and our 'screening is not encoding'
reading of M17 is wrong for the FEP's actual claim.
Seed 1718."""
import os, json
import numpy as np
from scipy.linalg import solve_continuous_are, solve_continuous_lyapunov
HERE = os.path.dirname(os.path.abspath(__file__))
ix = np.ix_
rng = np.random.default_rng(1718)

def r2(Cmm, Cme, Cee):
    # share of the variance of m that is linear in e
    return float(np.trace(Cme @ np.linalg.solve(Cee, Cme.T)) / np.trace(Cmm))

def gen(n, I, S, Ac, E, collider):
    hears = {**{k: I + S + Ac for k in I + Ac}, **{k: E + S + Ac for k in S + E}}
    if collider: hears[S[0]] = hears[S[0]] + [I[0]]
    A = -np.eye(n)
    for i in range(n):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: A[i, j] = rng.uniform(0.2, 0.9) * rng.choice([-1, 1])
    if collider: A[S[0], I[0]] = rng.uniform(0.3, 0.9) * rng.choice([-1, 1]); A[S[0], E[0]] = rng.uniform(0.3, 0.9) * rng.choice([-1, 1])
    return A

def snapshot(A, I, B, E):
    Sg = solve_continuous_lyapunov(A, -np.eye(len(A)))
    Kb = np.linalg.inv(Sg[ix(B, B)])
    WI = Sg[ix(I, B)] @ Kb; WE = Sg[ix(E, B)] @ Kb          # eta = WI x_B, mu = WE x_B
    Cbb = Sg[ix(B, B)]
    return r2(WE @ Cbb @ WE.T, WE @ Cbb @ WI.T, WI @ Cbb @ WI.T)

def history(A, I, B, E):
    n = len(A); IE = I + E; k = len(IE)
    F = np.zeros((k, k)); F[:len(I), :len(I)] = A[ix(I, I)]; F[len(I):, len(I):] = A[ix(E, E)]
    H = np.hstack([A[ix(B, I)], A[ix(B, E)]])
    P = solve_continuous_are(F.T, H.T, np.eye(k), np.eye(len(B))); P = (P + P.T) / 2
    K = P @ H.T
    Z = np.zeros((n + k, n + k)); Z[:n, :n] = A
    hat = list(range(n, n + k))
    Z[n:, n:] = F - K @ H
    Z[ix(hat, B)] += np.vstack([A[ix(I, B)], A[ix(E, B)]]) - K @ A[ix(B, B)]
    Z[n:, :n] += K @ A[B, :]
    N = np.zeros((n + k, n)); N[:n] = np.eye(n); N[ix(hat, B)] = K
    C = solve_continuous_lyapunov(Z, -N @ N.T)
    hI = hat[:len(I)]; hE = hat[len(I):]
    return r2(C[ix(hE, hE)], C[ix(hE, hI)], C[ix(hI, hI)])

def stable(A): return np.linalg.eigvals(A).real.max() < -0.05

out = {}
# positive control
ctl = []
while len(ctl) < 30:
    A = gen(6, [0, 1], [2], [3], [4, 5], False)
    if stable(A): ctl.append(snapshot(A, [0, 1], [2, 3], [4, 5]))
out['control_S_R2_min'] = float(min(ctl))
I, S, Ac, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + Ac
for key in ['split', 'collider']:
    rs, rh = [], []
    while len(rs) < 30:
        A = gen(8, I, S, Ac, E, key == 'collider')
        if stable(A): rs.append(snapshot(A, I, B, E)); rh.append(history(A, I, B, E))
    rs, rh = np.array(rs), np.array(rh)
    out[key] = dict(S_median=float(np.median(rs)), S_range=[float(rs.min()), float(rs.max())], S_frac_lt_0999=float(np.mean(rs < 0.999)),
                    H_median=float(np.median(rh)), H_range=[float(rh.min()), float(rh.max())], H_frac_lt_0999=float(np.mean(rh < 0.999)))
out['p1_control'] = out['control_S_R2_min'] > 1 - 1e-9
out['p2_snapshot'] = all(out[k]['S_frac_lt_0999'] >= 0.9 for k in ['split', 'collider'])
out['p3_history'] = all(out[k]['H_frac_lt_0999'] >= 0.9 for k in ['split', 'collider'])
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M17b.json'), 'w'), indent=1)
