"""M2 (loop 1 of plan/20261016-k1-math-deep-v1.md). Predictions fixed before running.
L1a Thm 6 thickness. 8 nodes: I = {0,1}, B1 = {2,3}, B2 = {4,5}, E = {6,7}, B = B1 u B2. 'thick': edges only inside a group and between
  adjacent layers (I-B1, B1-B2, B2-E), so d_E = 3. 'thin': any edge except I-E (d_E = 2). Edge prob 0.6, weights U(0.2,0.8) +-, diag -1,
  stability shift. 200 nets each, seed 9090. Predictions: thick median log-log slope of leak (h = 0.005 -> 0.02) in [3.8, 4.2];
  median |pred/exact - 1| at h = 0.02 < 0.05; median tau_B(thick) > median tau_B(thin) (tau_B = first h with leak > 0.5, grid to 20).
L1b Thm 7 no-go. X2 networks (seed 21). Prediction: among nets with an I-E path, 0 nets with leak < 1e-12 at any h in [0.01, 10].
L1c Thm 3' complex slowest pair. T4 networks (seed 5050, first 50) whose slowest mode is complex. Ratio of period-averaged sums at
  h = H + [0, 2 pi / omega) with H = max(20, 15 / gap), vs the |v| formula. Prediction: median over nets of max |diff| over candidates < 0.02.
L1d Thm 5' noise. (i) diagonal D ~ U(0.5, 2): median pred/exact at s = 5e-4 within 1%. (ii) full SPD D: leak at s = 0 above 1e-3 in >= 95%."""
import os, sys, json
import numpy as np
from scipy.linalg import expm
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'causal-partition-v0_6'))
import expC
from T1_graded_blankets import n, net
from T4_observer_dependence import rho_from
out = {}
def ldir(R2, src, mid, tgt):
    s = lambda t: R2[np.ix_(t, src)].sum(); d = s(mid) + s(tgt)
    return np.nan if d < 1e-300 else s(tgt) / d
def leak(R2, I, B, E): return np.nanmax([ldir(R2, I, B, E), ldir(R2, E, B, I)])
def stab(A):
    m = np.linalg.eigvals(A).real.max()
    return A - (m + 0.1) * np.eye(len(A)) if m > -0.1 else A

# L1a
rng = np.random.default_rng(9090); I, B1, B2, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = B1 + B2
grp = {k: g for g, ks in enumerate([I, B1, B2, E]) for k in ks}
def gen(thick):
    A = -np.eye(8)
    for i in range(8):
        for j in range(8):
            if i == j: continue
            gi, gj = grp[i], grp[j]
            ok = abs(gi - gj) <= 1 if thick else not ({gi, gj} == {0, 3})
            if ok and rng.random() < 0.6: A[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
    return stab(A)
hs = np.arange(0.05, 20.0001, 0.05); res = {}
for thick in [True, False]:
    sl, er, ta, t1 = [], [], [], []
    for g in range(200):
        A = gen(thick); L = lambda h: leak(expm(A * h) ** 2, I, B, E)
        if thick:
            A3 = np.linalg.matrix_power(A, 3)
            c = max((A3[np.ix_(E, I)] ** 2).sum() / max((A[np.ix_(B, I)] ** 2).sum(), 1e-300),
                    (A3[np.ix_(I, E)] ** 2).sum() / max((A[np.ix_(B, E)] ** 2).sum(), 1e-300))
            l1, l2 = L(0.005), L(0.02)
            if c > 0 and l1 > 0: sl.append(np.log(l2 / l1) / np.log(4)); er.append(abs((0.02 ** 2 / 6) ** 2 * c / l2 - 1))   # fix: skip nets without a length-3 path
        lf = np.array([L(h) for h in hs]); ta.append(float(hs[np.argmax(lf > 0.5)]) if (lf > 0.5).any() else np.inf)
        t1.append(float(hs[np.argmax(lf > 0.1)]) if (lf > 0.1).any() else np.inf)   # POST HOC (after run 1): tau at leak 0.1
    res['thick' if thick else 'thin'] = dict(slope_median=float(np.median(sl)) if sl else None, err_median=float(np.median(er)) if er else None,
                                             tau_median=float(np.median(ta)), tau_inf=int(np.isinf(ta).sum()), n_slope=len(sl),
                                             posthoc_tau01_median=float(np.median(t1)), posthoc_tau01_inf=int(np.isinf(t1).sum()))
out['L1a'] = res; print('L1a', res, flush=True)

# L1b
expC.rng = np.random.default_rng(21); bad, withpath, nopath = 0, 0, 0
for g in range(200):
    A, I_, B_, E_ = expC.random_planted(0.0); W = A[:6, :6]
    reach = np.linalg.matrix_power(np.abs(W) + np.eye(6), 6) > 0
    if not (reach[np.ix_(E_, I_)].any() or reach[np.ix_(I_, E_)].any()): nopath += 1; continue
    withpath += 1
    lf = np.array([leak(expm(W * h) ** 2, I_, B_, E_) for h in np.arange(0.01, 10.0001, 0.01)])
    bad += int(np.nanmin(lf) < 1e-12)
out['L1b'] = dict(with_path=withpath, no_path=nopath, nets_with_zero_leak=bad); print('L1b', out['L1b'], flush=True)

# L1c
rng = np.random.default_rng(5050); diffs = []
for g in range(50):
    A = net(rng, 0.0); rng.choice(n, 2, replace=False); Ao = A[:n, :n]
    lam, V = np.linalg.eig(Ao); k1 = np.argmax(lam.real)
    if abs(lam[k1].imag) < 1e-9: continue
    om = abs(lam[k1].imag); re = np.sort(np.unique(np.round(lam.real, 12)))[::-1]; gap = re[0] - re[1]
    H = max(20.0, 15 / gap); S = Ao - lam.real[k1] * np.eye(n)
    avg = np.mean([expm(S * h) ** 2 for h in H + np.arange(24) * (2 * np.pi / om) / 24], axis=0)
    r_avg = rho_from(avg, range(n)); r_f = rho_from(np.outer(np.abs(V[:, k1]) ** 2, np.ones(n)), range(n))
    ok = ~np.isnan(r_avg) & ~np.isnan(r_f); diffs.append(dict(maxdiff=float(np.abs(r_avg - r_f)[ok].max()), gap=float(gap), H=H,
                                                             spearman=float(spearmanr(r_avg[ok], r_f[ok]).statistic)))
out['L1c'] = dict(nets=len(diffs), maxdiff_median=float(np.median([d['maxdiff'] for d in diffs])),
                  spearman_median=float(np.median([d['spearman'] for d in diffs])), rows=diffs); print('L1c', {k: v for k, v in out['L1c'].items() if k != 'rows'}, flush=True)

# L1d
rng = np.random.default_rng(8181); I, B, E = [0, 1], [2, 3], [4, 5]; rat, lk0 = [], []
for g in range(200):
    M = rng.uniform(-1, 1, (6, 6)) * (rng.random((6, 6)) < 0.6); M = (M + M.T) / 2
    M[np.ix_(I, E)] = 0; M[np.ix_(E, I)] = 0; np.fill_diagonal(M, 0); Pi = M + (0.5 - np.linalg.eigvalsh(M).min()) * np.eye(6)
    Q = rng.standard_normal((6, 6)); Q = Q - Q.T; D = np.diag(rng.uniform(0.5, 2, 6)); s = 5e-4
    DP, QP = D @ Pi, Q @ Pi
    if (DP[np.ix_(B, I)] ** 2).sum() == 0 or (DP[np.ix_(B, E)] ** 2).sum() == 0: continue
    c = max((QP[np.ix_(E, I)] ** 2).sum() / (DP[np.ix_(B, I)] ** 2).sum(), (QP[np.ix_(I, E)] ** 2).sum() / (DP[np.ix_(B, E)] ** 2).sum())
    rat.append(s ** 2 * c / leak((-(D + s * Q) @ Pi * (1 - np.eye(6))) ** 2, I, B, E))
    G = rng.standard_normal((6, 6)); Df = G @ G.T / 6 + 0.3 * np.eye(6)
    lk0.append(leak((-(Df @ Pi) * (1 - np.eye(6))) ** 2, I, B, E))
out['L1d'] = dict(diag_pred_over_exact_median=float(np.median(rat)), full_D_leak_s0_above_1e3=float(np.mean(np.array(lk0) > 1e-3)),
                  full_D_leak_s0_median=float(np.median(lk0)))
print('L1d', out['L1d'])
json.dump(out, open(os.path.join(HERE, 'M2.json'), 'w'), indent=1, default=float)
