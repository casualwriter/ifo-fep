"""R4 (PRE-REGISTERED): graded boundary "tunnel vs room" on new real data, seen by observers with different time scales.
Plan: plan/20261016-r4-graded-leak-prereg.md (fixed before running).

Data: Causal Chambers wt_walks_v1, loads_hatch_mix_{fast,slow}_run_1..5 (10 files). Observer: push load_in; read pressure_ambient
(outside) and pressure_upwind (inside, reference); time scale h in {1,3,10,30,100} rows; raw sensor grain.
Per h: non-overlapping windows; drop windows with any sampling gap > 1.5 s. OLS dv = b dload_in + c1 dload_out + c2 dhatch + file
intercepts, HC1 SE. Leak L(h) = |b_ambient| / |b_upwind|; rho_O(h) = 1 - L(h). Bootstrap over windows (1000), 95% percentile CI.
Predictions: (1) L(1) upper < 0.005; (2) L(100) - L(1) CI > 0; (3) L(100) in [0.005, 0.10]; (4) signal_1 not significant at any h.
Death clauses: (1) L(1) upper >= 0.005; (2) L(100) - L(1) CI includes <= 0; (3) b_upwind never significant -> test void;
(4) negative control significant at any h -> design biased, result not used."""
import os, glob, json
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'data', 'causalchamber', 'wt_walks_v1')
FILES = [f'loads_hatch_mix_{s}_run_{i}.csv' for s in ['fast', 'slow'] for i in range(1, 6)]
HS, B, MARGIN = [1, 3, 10, 30, 100], 1000, 0.005
def windows(df, h, fid):
    t = df['timestamp'].to_numpy(); dt = np.diff(t); n = len(df)
    starts = np.arange(0, n - h, h)
    ok = np.array([dt[s:s + h].max() <= 1.5 for s in starts])
    s = starts[ok]; d = lambda c: df[c].to_numpy(float)[s + h] - df[c].to_numpy(float)[s]
    return dict(fid=np.full(len(s), fid), din=d('load_in'), dout=d('load_out'), dhat=d('hatch'),
                amb=d('pressure_ambient'), up=d('pressure_upwind'), sig=d('signal_1'), sec=float(np.median(dt)) * h)
def design(W):
    F = np.eye(W['fid'].max() + 1)[W['fid']]
    return np.column_stack([W['din'], W['dout'], W['dhat'], F])
def ols_b(X, y, se=False):
    XtXi = np.linalg.pinv(X.T @ X); b = XtXi @ X.T @ y
    if not se: return b[0]
    e = y - X @ b; n, k = X.shape; V = XtXi @ (X.T * e ** 2) @ X @ XtXi * n / (n - k)
    return b[0], np.sqrt(V[0, 0])
if __name__ == "__main__":
    dfs = [pd.read_csv(os.path.join(DATA, f)) for f in FILES]
    rng = np.random.default_rng(4040); res = {}; boots = {}
    thr_neg = norm.ppf(1 - 0.05 / (2 * len(HS)))
    for h in HS:
        parts = [windows(df, h, i) for i, df in enumerate(dfs)]
        W = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if k != 'sec'}
        sec = float(np.median([p['sec'] for p in parts])); X = design(W)
        b_amb, se_amb = ols_b(X, W['amb'], True); b_up, se_up = ols_b(X, W['up'], True); b_sig, se_sig = ols_b(X, W['sig'], True)
        L = abs(b_amb) / abs(b_up); n = len(W['din']); Lb = np.empty(B)
        for k in range(B):
            i = rng.integers(0, n, n); Xi = X[i]
            Lb[k] = abs(ols_b(Xi, W['amb'][i])) / abs(ols_b(Xi, W['up'][i]))
        boots[h] = Lb
        res[h] = dict(seconds=round(sec, 2), windows=n, b_ambient=float(b_amb), z_ambient=round(float(b_amb / se_amb), 2),
                      b_upwind=float(b_up), z_upwind=round(float(b_up / se_up), 2), L=float(L),
                      L_ci=[float(np.quantile(Lb, 0.025)), float(np.quantile(Lb, 0.975))], rho=float(1 - L),
                      z_signal1=round(float(b_sig / se_sig), 2))
        print(h, res[h], flush=True)
    diff = boots[100] - boots[1]; dci = [float(np.quantile(diff, 0.025)), float(np.quantile(diff, 0.975))]
    up_sig = any(abs(r['z_upwind']) > 1.96 for r in res.values()); neg = [h for h, r in res.items() if abs(r['z_signal1']) > thr_neg]
    verdict = dict(
        void=not up_sig, negative_control_failed=neg,
        prediction_1=res[1]['L_ci'][1] < MARGIN, prediction_2=dci[0] > 0,
        prediction_3=MARGIN <= res[100]['L'] <= 0.10, prediction_4=len(neg) == 0,
        death_clause_1='FIRED' if res[1]['L_ci'][1] >= MARGIN else 'not fired',
        death_clause_2='FIRED' if dci[0] <= 0 else 'not fired',
        death_clause_3='VOID' if not up_sig else 'ok', death_clause_4='FIRED' if neg else 'not fired',
        L100_minus_L1_ci=dci)
    print('VERDICT', verdict)
    json.dump(dict(results={str(k): v for k, v in res.items()}, verdict=verdict), open(os.path.join(HERE, 'R4.json'), 'w'), indent=1)
