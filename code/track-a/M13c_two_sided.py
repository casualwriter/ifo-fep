"""M13c (POST HOC, after M13b prediction (c) failed: one hidden slow interior mode did NOT slow forgetting).
Explanation from the Riccati equation: with M = 0 the cross block obeys dS_IE/dt = C_I S_IE + S_IE C_E', with closed-loop matrices
C_I = A_II - S_II G_II and C_E = A_EE - S_EE G_EE, so S_IE decays at (slowest rate of C_I) + (slowest rate of C_E).
Forgetting is slow only if BOTH sides hide a slow mode from the blanket. Test: hidden slow node 0 in I AND hidden slow node 6 in E
(each drives nothing, self-rate -0.02), 30 networks: expect MI(40)/MI(0) > 1e-3 and log-slope of MI over L in [10, 40] near
-2 x (0.02 + 0.02) = -0.08 (range [-0.16, -0.04]). Control: hidden slow node only in I (as M13b (c)): fast."""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M13b_detectability.py"), encoding="utf-8").read().split("\nres = {")[0])
rng = np.random.default_rng(1314)
def hide(A, node, others):
    A[:, node] = 0.0; A[node, :] = 0.0; A[node, node] = -0.02
    for j in others:
        if rng.random() < 0.6: A[node, j] = rng.uniform(0.2, 0.9) * rng.choice([-1, 1])
    return A
Ls = (0, 5, 10, 20, 40); res = {'both_sides': [], 'one_side': []}
while len(res['both_sides']) < 30 or len(res['one_side']) < 30:
    A = gen()
    A1 = hide(A.copy(), 0, [1] + S + Ac)
    A2 = hide(A1.copy(), 6, [7] + S + Ac)
    if len(res['one_side']) < 30 and stable(A1, 0.01): res['one_side'].append(mi_curve(A1, Ls))
    if len(res['both_sides']) < 30 and stable(A2, 0.01): res['both_sides'].append(mi_curve(A2, Ls))
out = {}
for k, v in res.items():
    v = np.array(v); r = v[:, -1] / np.maximum(v[:, 0], 1e-300)
    sl = [np.polyfit([L for L, m in zip(Ls[2:], row[2:]) if m > 1e-14], np.log([m for m in row[2:] if m > 1e-14]), 1)[0]
          for row in v if sum(m > 1e-14 for m in row[2:]) >= 2]
    out[k] = dict(n=len(v), median_MI=np.median(v, 0).tolist(), median_ratio_40_0=float(np.median(r)),
                  median_logslope=float(np.median(sl)) if sl else None, n_slope=len(sl))
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M13c.json'), 'w'), indent=1)
