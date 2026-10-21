"""F5b: small pushes with the SHAM-controlled test (F3b), which is the version the paper reports (3% / 9% false rejection).
F5 used the F3 test without sham; under the hidden drive its false rejection was 12-24% and grew with K, so F5 does not
measure the paper's test. F5b repeats the F5 grid with F3b's test (code copied from F3b_sham.py, which has no main guard).

Fixed before running: same networks as F5 (seed 2020, 50 per c), u in {0.5, 1, 2, 3}, K in {1000, 4000}; sham = u = 0
pushes from the same kind of start states; two-sample z per I-E pair, Bonferroni 0.05.
Prediction: false rejection <= 10% in every cell (the sham cancels the hidden-drive bias); detection depends on u^2 K only:
 (u=1, K=4000) ~ (u=2, K=1000) and (u=0.5, K=4000) ~ (u=1, K=1000).
Death clauses: (1) false rejection > 15% in any cell -> sham test not calibrated at that push size.
 (2) a matched pair with non-overlapping 95% Wilson intervals for detection -> u^2 K scaling fails."""
import sys, os, json, time
import numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'causal-partition-v0_6'))
import F3_verify20 as F
from stats import wilson
n, dt = F.n, F.dt
def trials(A, X, j, u, K, rng, F1, F2):
    idx = rng.integers(0, len(X) - 2, K); x = np.zeros((K, n + 1)); x[:, :n] = X[idx]; x[:, n] = rng.standard_normal(K) * np.sqrt(0.5 / abs(A[n, n]))
    x0 = x[:, :n].copy(); x[:, j] += u; sq = np.sqrt(dt)
    for s in range(10):
        x = x + x @ A.T * dt + rng.standard_normal((K, n + 1)) * sq
        if s == 4: r1 = x[:, :n] - x0 @ F1
    r2 = x[:, :n] - x0 @ F2
    return (4 * r1 - r2) / (2 * 0.05)
def sham_test(A, X, cand, K, u, rng, alpha=0.05):
    I, B, E = cand; F1 = np.linalg.lstsq(X[:-1], X[1:], rcond=None)[0]; F2 = np.linalg.lstsq(X[:-2], X[2:], rcond=None)[0]; z = []
    for src, tgt in [(I, E), (E, I)]:
        for j in src:
            p = trials(A, X, j, u, K, rng, F1, F2); s = trials(A, X, j, 0.0, K, rng, F1, F2)
            for t in tgt:
                d = p[:, t].mean() - s[:, t].mean(); se = np.sqrt(p[:, t].var(ddof=1) / K + s[:, t].var(ddof=1) / K); z.append(d / se)
    return bool(np.max(np.abs(z)) > norm.ppf(1 - alpha / (2 * len(z))))
U, KS, G = [0.5, 1.0, 2.0, 3.0], [1000, 4000], 50
if __name__ == "__main__":
    out = {}; t0 = time.time()
    for c in [0.0, 0.4]:
        F.rng = np.random.default_rng(2020)
        nets = [F.net(c) for _ in range(G)]
        data = [F.simulate(A, 10000, 5) for A, *_ in nets]
        for u in U:
            for K in KS:
                rng = np.random.default_rng(int(1000 * u) + K)
                fr, det = [], []
                for (A, true, bad, isbad), X in zip(nets, data):
                    fr.append(sham_test(A, X, true, K, u, rng))
                    if isbad: det.append(sham_test(A, X, bad, K, u, rng))
                key = f'c={c}_u={u}_K={K}'
                out[key] = dict(false_reject=wilson(int(sum(fr)), len(fr)), detect=wilson(int(sum(det)), len(det)), n_bad=len(det))
                print(key, out[key], round(time.time() - t0), 's', flush=True)
    json.dump(out, open(os.path.join(HERE, 'F5b.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
