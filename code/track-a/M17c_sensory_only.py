"""M17c: two follow-ups from the round-3 review of fep-blanket v0.7. Prospectively specified.
(a) Vector: the interior drives the active states; does a reader need them? Causal reader of the SENSORY states only, (-inf, t],
    in split networks (M = 0, wiring cut). Hidden: (x_I, x_A, x_E). Score: I(x_I(t); x_E(t) | past of S) from the stationary filter.
    Control: reader of the whole blanket (S and A) gives 0 (Proposition 1c).
    Prediction: sensory-only MI > 1e-4 in >= 90% of 30 split networks; whole-blanket MI < 1e-10 in all.
    Death: if sensory-only MI is ~0, a reader does not need its active states, and the efference-copy reading is wrong.
(b) Nova: show that encoding does not need screening. 6-node networks, I = {0,1}, B = {2,3}, E = {4,5}, blanket state 2 hears both
    sides (M != 0). Snapshot synchronisation map: R2(mu on eta).
    Prediction: R2 > 1 - 1e-9 in all 30 (map exists), while causal-reader MI > 1e-4 in all 30 (no screening).
    With M17b (split networks: screening without encoding), this gives both non-implications by example.
Seed 1719.
Code fix after run 1 (disclosed): in (b) the generator did not force node 2 to hear E, so some networks had M = 0 (b_M_min = 0)
and the (b) prediction failed for that reason. Run 1 output kept in M17c_run1_buggy.json. Fix: also force A[2, 4] != 0.
Second fix (run 2 crashed): W_I can be rank-deficient, so the tool now uses a pseudo-inverse; all networks are kept, and the
rank of W_I is reported (a map needs ker W_I in ker W_E, which a rank-deficient W_I can break)."""
import os, sys, json
import numpy as np
from scipy.linalg import solve_continuous_are
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'tools'))
from blanket_report import causal_leak, sync_r2, M_matrix
ix = np.ix_
rng = np.random.default_rng(1719)
ld = lambda X: np.linalg.slogdet(X)[1]

def gen(n, hears, extra=()):
    A = -np.eye(n)
    for i in range(n):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: A[i, j] = rng.uniform(0.2, 0.9) * rng.choice([-1, 1])
    for (i, j) in extra: A[i, j] = rng.uniform(0.3, 0.9) * rng.choice([-1, 1])
    return A

def stable(A): return np.linalg.eigvals(A).real.max() < -0.05

def sensory_only(A, I, S, Ac, E):
    Z = I + Ac + E; F = A[ix(Z, Z)]; H = A[ix(S, Z)]
    P = solve_continuous_are(F.T, H.T, np.eye(len(Z)), np.eye(len(S))); P = (P + P.T) / 2
    a, e = list(range(len(I))), list(range(len(I) + len(Ac), len(Z)))
    return max(0.5 * (ld(P[ix(a, a)]) + ld(P[ix(e, e)]) - ld(P[ix(a + e, a + e)])), 0.0)

out = {}
# (a)
I, S, Ac, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + Ac
split = {**{k: I + S + Ac for k in I + Ac}, **{k: E + S + Ac for k in S + E}}
so, wb = [], []
while len(so) < 30:
    A = gen(8, split)
    if stable(A): so.append(sensory_only(A, I, S, Ac, E)); wb.append(causal_leak(A, np.eye(8), I, B, E))
so, wb = np.array(so), np.array(wb)
out['a_sensory_only'] = dict(min=float(so.min()), median=float(so.median() if hasattr(so, 'median') else np.median(so)), max=float(so.max()), frac_gt_1e4=float(np.mean(so > 1e-4)))
out['a_whole_blanket_max'] = float(wb.max())
# (b)
I6, B6, E6 = [0, 1], [2, 3], [4, 5]
h6 = {0: [0, 1, 2, 3], 1: [0, 1, 2, 3], 3: [0, 1, 2, 3], 2: [2, 3, 4, 5], 4: [2, 3, 4, 5], 5: [2, 3, 4, 5]}   # 2 hears E (sensory), 3 hears I (active)
r2s, cl, Ms, ranks = [], [], [], []
while len(r2s) < 30:
    A = gen(6, h6, extra=[(2, 0), (2, 4)])               # sensory 2 hears interior 0 and exterior 4: M != 0
    if stable(A):
        r2s.append(sync_r2(A, np.eye(6), I6, B6, E6, 'snapshot')); cl.append(causal_leak(A, np.eye(6), I6, B6, E6))
        Ms.append(np.abs(M_matrix(A, np.eye(6), I6, B6, E6)).max())
        from scipy.linalg import solve_continuous_lyapunov
        Sg = solve_continuous_lyapunov(A, -np.eye(6)); ranks.append(int(np.linalg.matrix_rank(Sg[ix(I6, B6)], tol=1e-10)))
r2s, cl = np.array(r2s), np.array(cl)
out['b_snapshot_R2_min'] = float(r2s.min()); out['b_causal_MI_min'] = float(cl.min()); out['b_causal_MI_median'] = float(np.median(cl))
out['b_M_min'] = float(min(Ms)); out['b_rank_WI'] = ranks; out['b_R2'] = [float(x) for x in r2s]
out['pred_a'] = out['a_sensory_only']['frac_gt_1e4'] >= 0.9 and out['a_whole_blanket_max'] < 1e-10
out['pred_b'] = out['b_snapshot_R2_min'] > 1 - 1e-9 and out['b_causal_MI_min'] > 1e-4
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M17c.json'), 'w'), indent=1)
