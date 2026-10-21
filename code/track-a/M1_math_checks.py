"""M1: numerical checks of the derivations in plan/20261016-k1-math-derivations-v1.md (written before this run).
Linear system dx = A x dt + noise, response to a push after time h: R(h) = expm(A h). Leak 1 - rho as in T1/X2.

C1 Theorem 1 (leak grows as h^2). Valid blanket (A_IE = A_EI = 0): 1 - rho(h) = (h^2/4) kappa + O(h^3),
   kappa = max over directions of ||A_EB A_BI||_F^2 / ||A_BI||_F^2. Ratio form rho ~ 1 / (1 + h^2 kappa / 4) gives tau_B ~ 2 / sqrt(kappa).
   Networks: X2 networks (6 nodes, seed 21). Predictions: log-log slope of leak on h in [0.01, 0.05] within [1.9, 2.1] (median);
   median |pred/exact - 1| at h = 0.1 below 0.05. tau_B pred/exact: reported only (leading order, no prediction).
C2 Theorem 3 (slow-observer limit). If the slowest mode of A is real and simple (right eigenvector v),
   rho(h -> inf) = ||v_B||^2 / (||v_B||^2 + max(||v_I||^2, ||v_E||^2)): the slow observer ranks candidates by where the slow mode lives,
   not by the wiring. Networks: T4 networks (12 nodes, seed 5050, first 50). Predictions: median Spearman(rho_h, rho_inf) over
   candidates rises with h and is >= 0.9 at h = 30 (networks with a real slowest mode); at h = 3 (T4's slow observer)
   Spearman(rho_3, rho_inf) > Spearman(rho_3, rho_0) in the median.
   Post hoc (T4 already seen): candidates with internal positive feedback put a larger share of the slow mode on I.
   Constructive: the candidate that maximises min(rho_0, rho_inf) ("robust blanket") -- report its min(rho_0, rho_3).
C3 Theorem 4 (no sham -> false alarms grow with K). If the reference model has a fixed error b on the start states,
   the no-sham z statistic is N(b sqrt(K) / sigma, 1), so P(reject) -> 1 as K grows; the sham-difference statistic is N(0, 1).
   Model: 2 observed nodes, no coupling 0 -> 1, both driven by a hidden slow drive; start states as in F5 (hidden part redrawn).
   Prediction: empirical false rejection (200 repeats per K, same data) within 0.08 of the formula for K in {250, 1000, 4000};
   sham false rejection <= 0.08 for every K.
C4 Power rule. Per-trial slope SD = sqrt(2.5 / h1) / u for unit noise (h1 = 0.05: 7.07 / u). Prediction: empirical SD within 10%.
   Consequence: sham test needs u^2 K ~ 100 (z_c + z_beta)^2 / a^2 (h1 = 0.05).
C5 Theorem 5 (snapshot blanket, push observer). A = -(I + s Q) Pi, Pi_IE = 0, Q antisymmetric: leak = O(s^2),
   leak_IE ~ s^2 ||(Q Pi)_EI||^2 / ||Pi_BI||^2. Prediction: median log-log slope in s in [1.9, 2.1]; median pred/exact within 10% at s = 0.05.
Death: any prediction outside its range -> the derivation is wrong or its assumptions do not hold; report and fix the note."""
import os, sys, json
import numpy as np
from scipy.linalg import expm
from scipy.stats import spearmanr, norm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'causal-partition-v0_6'))
import expC
from T1_graded_blankets import n, net, PARTS
from T4_observer_dependence import rho_from, pos_cycle
out = {}

def leak_dir(R2, src, mid, tgt):
    s = lambda t: R2[np.ix_(t, src)].sum()
    d = s(mid) + s(tgt)
    return np.nan if d < 1e-300 else s(tgt) / d
def leak(R2, I, B, E): return np.nanmax([leak_dir(R2, I, B, E), leak_dir(R2, E, B, I)])

# ---------- C1 ----------
expC.rng = np.random.default_rng(21); rows = []
for g in range(200):
    A, I, B, E = expC.random_planted(0.0); W = A[:6, :6]
    ks = []
    for src, tgt in [(I, E), (E, I)]:
        d = (W[np.ix_(B, src)] ** 2).sum()
        if d > 0: ks.append(((W[np.ix_(tgt, B)] @ W[np.ix_(B, src)]) ** 2).sum() / d)
    if not ks or max(ks) == 0: continue
    k = max(ks); L = lambda h: leak(expm(W * h) ** 2, I, B, E)
    l1, l2 = L(0.01), L(0.05); slope = np.log(l2 / l1) / np.log(5)
    hs = np.linspace(0.01, 10, 1000); lf = np.array([L(h) for h in hs[::10]])
    tau = float(hs[::10][np.argmax(lf > 0.5)]) if (lf > 0.5).any() else None
    rows.append(dict(slope=slope, err01=abs((0.01 * k / 4) / L(0.1) - 1), tau=tau, tau_pred=2 / np.sqrt(k)))
tr = [r['tau_pred'] / r['tau'] for r in rows if r['tau']]
out['C1'] = dict(networks=len(rows), slope_median=float(np.median([r['slope'] for r in rows])),
                 err_h01_median=float(np.median([r['err01'] for r in rows])),
                 tau_pred_over_exact_median=float(np.median(tr)), tau_pred_over_exact_q10_q90=np.quantile(tr, [0.1, 0.9]).tolist(), n_tau=len(tr))
print('C1', out['C1'], flush=True)

# ---------- C2 ----------
rng = np.random.default_rng(5050); rows = []; fbI, noI = [], []
for g in range(50):
    A = net(rng, 0.0); rng.choice(n, 2, replace=False)               # keep the T4 random stream
    Ao = A[:n, :n]; W = Ao * (1 - np.eye(n))
    lam, V = np.linalg.eig(Ao); k1 = np.argmax(lam.real); real = abs(lam[k1].imag) < 1e-9
    gap = lam.real[k1] - np.sort(lam.real)[-3 if not real else -2]
    v2 = np.abs(V[:, k1]) ** 2
    S = Ao - lam.real[k1] * np.eye(n)                                    # shift: ratios do not change, no underflow
    r0 = rho_from(W ** 2, range(n)); rinf = rho_from(np.outer(v2, np.ones(n)), range(n))
    rh = {h: rho_from(expm(S * h) ** 2, range(n)) for h in [1, 3, 10, 30]}
    ok = ~np.isnan(r0) & ~np.isnan(rinf)
    for h in rh: ok &= ~np.isnan(rh[h])
    sp = lambda a, b: float(spearmanr(a[ok], b[ok]).statistic)
    rob = np.argmax(np.where(ok, np.minimum(r0, rinf), -1)); b1 = np.argmax(np.where(ok, r0, -1))
    row = dict(real=bool(real), gap=float(gap), **{f'sp_inf_h{h}': sp(rh[h], rinf) for h in rh},
               sp_0_h3=sp(rh[3], r0), best3_eq_bestinf=bool(np.argmax(np.where(ok, rh[3], -1)) == np.argmax(np.where(ok, rinf, -1))),
               best30_eq_bestinf=bool(np.argmax(np.where(ok, rh[30], -1)) == np.argmax(np.where(ok, rinf, -1))),
               robust_min_r0_r3=float(min(r0[rob], rh[3][rob])), o1best_min_r0_r3=float(min(r0[b1], rh[3][b1])))
    rows.append(row)
    if real and g < 20:                                                  # post hoc: slow-mode share on I, with / without internal loop
        for p in rng.choice(np.nonzero(ok & (r0 >= 0.5))[0], 300, replace=False):
            I = list(PARTS[p][0]); (fbI if pos_cycle(W, I) else noI).append(float(v2[I].sum() / v2.sum()))
R = [r for r in rows if r['real']]
med = lambda key, rr=R: float(np.median([r[key] for r in rr]))
out['C2'] = dict(networks=len(rows), real_slowest=len(R), **{k: med(k) for k in ['sp_inf_h1', 'sp_inf_h3', 'sp_inf_h10', 'sp_inf_h30', 'sp_0_h3', 'gap']},
                 best3_eq_bestinf=float(np.mean([r['best3_eq_bestinf'] for r in R])), best30_eq_bestinf=float(np.mean([r['best30_eq_bestinf'] for r in R])),
                 robust_min_r0_r3=med('robust_min_r0_r3', rows), o1best_min_r0_r3=med('o1best_min_r0_r3', rows),
                 posthoc_slowmode_share_on_I=dict(with_loop=float(np.median(fbI)), without=float(np.median(noI)), n=(len(fbI), len(noI))))
print('C2', out['C2'], flush=True)

# ---------- C3 + C4 ----------
rng = np.random.default_rng(7070); dt, h1 = 0.01, 0.05
A = np.array([[-1.0, 0, 0.8], [0, -1.0, 0.8], [0, 0, -0.1]])           # nodes 0, 1 observed; 2 = hidden slow drive; no 0 -> 1 coupling
x = np.zeros(3); X = np.empty((10000, 2))
for _ in range(3000): x = x + A @ x * dt + rng.standard_normal(3) * np.sqrt(dt)
for k in range(10000):
    for _ in range(5): x = x + A @ x * dt + rng.standard_normal(3) * np.sqrt(dt)
    X[k] = x[:2]
F1 = np.linalg.lstsq(X[:-1], X[1:], rcond=None)[0]; F2 = np.linalg.lstsq(X[:-2], X[2:], rcond=None)[0]
def trials(K, u):
    idx = rng.integers(0, len(X) - 2, K); s = np.zeros((K, 3)); s[:, :2] = X[idx]; s[:, 2] = rng.standard_normal(K) * np.sqrt(0.5 / 0.1)
    x0 = s[:, :2].copy(); s[:, 0] += u
    for st in range(10):
        s = s + s @ A.T * dt + rng.standard_normal((K, 3)) * np.sqrt(dt)
        if st == 4: r1 = s[:, :2] - x0 @ F1
    return ((4 * r1 - (s[:, :2] - x0 @ F2)) / (2 * h1))[:, 1]            # target node 1, not divided by u
u = 3.0; big = trials(200000, 0.0); b = big.mean(); sig = big.std()     # reference-model error b on the start states, per-trial SD
c3 = {}; zc = norm.ppf(1 - 0.05 / 2)
for K in [250, 1000, 4000]:
    ns, sh = [], []
    for rep in range(200):
        p = trials(K, u) / u; q = trials(K, 0.0) / u
        ns.append(abs(p.mean() / (p.std(ddof=1) / np.sqrt(K))) > zc)
        sh.append(abs((p.mean() - q.mean()) / np.sqrt(p.var(ddof=1) / K + q.var(ddof=1) / K)) > zc)
    m = (b / u) * np.sqrt(K) / (sig / u); pred = norm.sf(zc - m) + norm.cdf(-zc - m)
    c3[K] = dict(no_sham=float(np.mean(ns)), formula=float(pred), sham=float(np.mean(sh)))
out['C3'] = dict(bias_b=float(b), sd=float(sig), by_K=c3)
out['C4'] = dict(sd_per_trial_u1=float(sig), formula=float(np.sqrt(2.5 / h1)))
print('C3', out['C3'], '\nC4', out['C4'], flush=True)

# ---------- C5 ----------
rng = np.random.default_rng(8080); rows = []
I, B, E = [0, 1], [2, 3], [4, 5]
for g in range(200):
    M = rng.uniform(-1, 1, (6, 6)) * (rng.random((6, 6)) < 0.6); M = (M + M.T) / 2
    M[np.ix_(I, E)] = 0; M[np.ix_(E, I)] = 0; np.fill_diagonal(M, 0)
    Pi = M + (0.5 - np.linalg.eigvalsh(M).min()) * np.eye(6)
    Q = rng.standard_normal((6, 6)); Q = Q - Q.T
    L = lambda s: leak((-(np.eye(6) + s * Q) @ Pi * (1 - np.eye(6))) ** 2, I, B, E)
    QP = Q @ Pi
    pred = 0.05 ** 2 * max((QP[np.ix_(E, I)] ** 2).sum() / (Pi[np.ix_(B, I)] ** 2).sum(), (QP[np.ix_(I, E)] ** 2).sum() / (Pi[np.ix_(B, E)] ** 2).sum())
    if (Pi[np.ix_(B, I)] ** 2).sum() == 0 or (Pi[np.ix_(B, E)] ** 2).sum() == 0: continue
    rows.append(dict(slope=np.log(L(0.01) / L(0.001)) / np.log(10), ratio=pred / L(0.05)))
out['C5'] = dict(networks=len(rows), slope_median=float(np.median([r['slope'] for r in rows])),
                 pred_over_exact_median=float(np.median([r['ratio'] for r in rows])))
print('C5', out['C5'])
json.dump(out, open(os.path.join(HERE, 'M1.json'), 'w'), indent=1, default=str)
