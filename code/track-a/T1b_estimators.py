"""T1b: fix the crossing-share estimator near rho = 1. T1 found estimates biased low near 1 (noise floor; death clause 2).
Predictions and decision rule fixed before running: plan/20261016-manifest-probability-crossing-share-v1.md, section 5.2.

Same design as T1, new network seed 4040, 4000 push and 4000 sham trials per node. Four estimators on the same networks and data:
 E0  max(a^2 - v, 0), first 1000 trials (T1 method)
 Ea  max(a^2 - v, 0), all 4000 trials
 Eb  a^2 - v untruncated, first 1000; numerator and denominator summed first, rho clipped to [0, 1]; denominator <= 0 -> failure (NaN)
 Ec  a^2 if |a|/sqrt(v) > Bonferroni threshold (132 couplings per network, alpha 0.05) else 0, first 1000
Main criterion (each estimator, both c): control rho_hat >= 0.95 in >= 80%; rank-correct not more than 5 points below E0; MAE <= 0.10.
Decision: among those meeting it, the lowest MAE becomes the standard estimator.
Death clause: at c = 0 no estimator reaches control >= 0.95 in >= 80% -> absolute rho near 1 not measurable with K <= 4000."""
import os, sys, json, time
import numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import T1_graded_blankets as T
from T1_graded_blankets import n, dt, U, H1, BINS, MASK, net, simulate, J2_obs, rho_all
from stats import wilson
KMAX, KLOW = 4000, 1000
ZB = norm.ppf(1 - 0.05 / (2 * n * (n - 1)))
def slopes(A, X, j, u, rng, F1, F2, K):
    idx = rng.integers(0, len(X) - 2, K); x = np.zeros((K, n + 1)); x[:, :n] = X[idx]
    x[:, n] = rng.standard_normal(K) * np.sqrt(0.5 / abs(A[n, n]))
    x0 = x[:, :n].copy(); x[:, j] += u; sq = np.sqrt(dt)
    for s in range(10):
        x = x + x @ A.T * dt + rng.standard_normal((K, n + 1)) * sq
        if s == 4: r1 = x[:, :n] - x0 @ F1
    r2 = x[:, :n] - x0 @ F2
    return (4 * r1 - r2) / (2 * H1)
def ests(A, X, rng):
    F1 = np.linalg.lstsq(X[:-1], X[1:], rcond=None)[0]; F2 = np.linalg.lstsq(X[:-2], X[2:], rcond=None)[0]
    a = {k: np.zeros((n, n)) for k in ['lo', 'hi']}; v = {k: np.zeros((n, n)) for k in ['lo', 'hi']}
    for j in range(n):
        p = slopes(A, X, j, U, rng, F1, F2, KMAX); s = slopes(A, X, j, 0.0, rng, F1, F2, KMAX)
        for k, K in [('lo', KLOW), ('hi', KMAX)]:
            a[k][:, j] = (p[:K].mean(0) - s[:K].mean(0)) / U
            v[k][:, j] = (p[:K].var(0, ddof=1) / K + s[:K].var(0, ddof=1) / K) / U ** 2
    z = np.abs(a['lo']) / np.sqrt(v['lo'])
    return dict(E0=np.maximum(a['lo'] ** 2 - v['lo'], 0), Ea=np.maximum(a['hi'] ** 2 - v['hi'], 0),
                Eb=a['lo'] ** 2 - v['lo'], Ec=np.where(z > ZB, a['lo'] ** 2, 0.0))
def rho_untrunc(J2):
    J2 = J2 * (1 - np.eye(n))
    s = lambda x, y: np.einsum('pt,ts,ps->p', MASK[x], J2, MASK[y])
    out = []
    for num, leak in [(s('B', 'I'), s('E', 'I')), (s('B', 'E'), s('I', 'E'))]:
        den = num + leak
        with np.errstate(invalid='ignore', divide='ignore'):
            out.append(np.where(den > 0, np.clip(num / den, 0, 1), np.nan))
    return np.minimum(out[0], out[1])
def run(c, G, seed=4040, max_tries=400):
    rng = np.random.default_rng(seed + int(10 * c)); rows = []; skipped = 0; tries = 0
    while len(rows) < G and tries < max_tries:
        tries += 1
        A = net(rng, c); J2 = A[:n, :n] ** 2; r = rho_all(J2); pick = {}
        for b, (lo, hi) in BINS.items():
            ok = np.nonzero((r >= lo) & (r < hi))[0]
            if len(ok): pick[b] = int(rng.choice(ok))
        ctrl = np.nonzero(r == 1)[0]
        if len(pick) < 3: skipped += 1; continue
        if len(ctrl): pick['control'] = int(rng.choice(ctrl))
        X = simulate(A, rng); E = ests(A, X, rng)
        R = {k: (rho_untrunc(E[k]) if k == 'Eb' else rho_all(E[k])) for k in E}; R['obs'] = rho_all(J2_obs(X))
        rows.append({b: dict(true=float(r[p]), **{k: float(R[k][p]) for k in R}) for b, p in pick.items()})
    return rows, skipped
def summary(rows):
    out = {}
    for m in ['E0', 'Ea', 'Eb', 'Ec', 'obs']:
        rank = [bool(r['low'][m] < r['mid'][m] < r['high'][m]) for r in rows]
        err = np.array([abs(r[b][m] - r[b]['true']) for r in rows for b in BINS])
        ctl = [bool(r['control'][m] >= 0.95) for r in rows if 'control' in r]
        bias = {b: round(float(np.nanmean([r[b][m] - r[b]['true'] for r in rows if b in r])), 3) for b in list(BINS) + ['control']}
        out[m] = dict(rank_correct=wilson(int(sum(rank)), len(rank)), mae=round(float(np.nanmean(err)), 3), n_failed=int(np.isnan(err).sum()),
                      control_ge_095=wilson(int(sum(ctl)), len(ctl)), n_control=len(ctl), bias=bias)
    return out
if __name__ == "__main__":
    G = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    out = {}; t0 = time.time()
    for c in [0.0, 0.4]:
        rows, skipped = run(c, G)
        out[f'c={c}'] = dict(summary=summary(rows), kept=len(rows), skipped=skipped, rows=rows)
        print(c, 'kept', len(rows), 'skipped', skipped, round(time.time() - t0), 's', flush=True)
        for m, s in out[f'c={c}']['summary'].items():
            print('  ', m, 'rank', [round(x, 3) for x in s['rank_correct']], 'mae', s['mae'], 'fail', s['n_failed'],
                  'ctl>=.95', [round(x, 3) for x in s['control_ge_095']], 'n_ctl', s['n_control'], 'bias', s['bias'], flush=True)
    json.dump(out, open(os.path.join(HERE, f'T1b_G{G}.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
