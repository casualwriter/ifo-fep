"""M4b (loop 3 follow-up, referee's discriminating test; predictions fixed before running).
Thm 9 says the floor depends on the system only through M = A_BI' D_B^-1 A_BE. Loop 3 nets (seed 9191, 200), with A[2,7] = 0 and
A[2,6] = 0.5 (so A_2E is along e_6), eta = 0.1 (D[2,4] = D[4,2]). Then M = A_4I' [Dinv_44 eps e_6' + Dinv_42 A_2E] = 0 at
eps* = -Dinv_42 * A[2,6] / Dinv_44 (exact). Conditions: (a) eta only; (b) eps = eps* with eta = 0; (c) both (cancel, M = 0);
(d) eps = -eps* with eta (M doubled).
Predictions: plateau m = 16 vs 32 equal within 5% for (a), (b), (d); (c) median floor < 5% of (a) and within 10x of the eps = eta = 0
baseline; (d) median floor / (a) in [3, 5] (|M| doubled -> 4x). Death: (c) > 20% of (a) -> floor is not a function of M alone."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "M3_loop2.py"), encoding="utf-8").read().split("\nout = {}")[0])
I, S, A_, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + A_
hears = {**{k: I + S + A_ for k in I + A_}, **{k: E + S + A_ for k in S + E}}
rng = np.random.default_rng(9191); su = 0.1
def gen():
    M = -np.eye(8)
    for i in range(8):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: M[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
    M[4, 0] = 0.5; return M
def stab(M):
    m = np.linalg.eigvals(M).real.max(); return M - (m + 0.1) * np.eye(8) if m > -0.1 else M
def floor(M, D, m):
    o = Obs(M, 2.0, m, D=D); P = su * np.eye(8)
    return np.nanmax([o.lc(P[:, I], I, B, E), o.lc(P[:, E], E, B, I)])
eta = 0.1; D = np.eye(8); D[2, 4] = D[4, 2] = eta; Di = np.linalg.inv(D[np.ix_(B, B)])   # B order 2,3,4,5 -> idx 2 = node 4, 0 = node 2
eps_star = -Di[2, 0] * 0.5 / Di[2, 2]
res = {k: {16: [], 32: []} for k in ['base', 'a_eta', 'b_eps', 'c_cancel', 'd_double']}
for g in range(200):
    M0 = gen(); M0[2, 7] = 0.0; M0[2, 6] = 0.5
    for k, (eps, DD) in dict(base=(0, None), a_eta=(0, D), b_eps=(eps_star, None), c_cancel=(eps_star, D), d_double=(-eps_star, D)).items():
        M = M0.copy(); M[4, 6] = eps; M = stab(M)
        for m in (16, 32): res[k][m].append(floor(M, DD, m))
med = {k: {m: float(np.nanmedian(v)) for m, v in d.items()} for k, d in res.items()}
c_over_a = float(np.nanmedian(np.array(res['c_cancel'][32]) / np.array(res['a_eta'][32])))
d_over_a = float(np.nanmedian(np.array(res['d_double'][32]) / np.array(res['a_eta'][32])))
out = dict(eps_star=float(eps_star), median=med, c_over_a_median=c_over_a, d_over_a_median=d_over_a,
           plateau_16_vs_32={k: med[k][16] / med[k][32] for k in med})
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'M4b.json'), 'w'), indent=1)
