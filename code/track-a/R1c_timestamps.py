"""R1c (POST-HOC, found after R1): sampling is not uniform. The step into row 0 of each block (the sham) takes ~0.93 s,
into row 5 (the push) ~0.17 s, other steps ~0.15 s. R1 extrapolated per sample, so the sham baseline used a slope over an
uneven time span. Here the same test uses real time: s = (v[t] - v[t-2]) / (tau[t] - tau[t-2]),
prediction v[t+k] = v[t] + s (tau[t+k] - tau[t]). Windows 4 (R1) and 8 samples. Question: do R1's conclusions hold?"""
import os, glob, json
import numpy as np, pandas as pd
from scipy.stats import norm
import R1_wt_impulse as R
def resp(v, tau, t, w):
    s = (v[:, t] - v[:, t - 2]) / (tau[:, t] - tau[:, t - 2])
    return (v[:, t + 1:t + 1 + w] - (v[:, [t]] + s[:, None] * (tau[:, t + 1:t + 1 + w] - tau[:, [t]]))).mean(1)
out = {}
for f in sorted(glob.glob(os.path.join(R.DATA, '*.csv'))):
    df = pd.read_csv(f); B = df['flag'].nunique()
    X = {k: df[k].to_numpy(float).reshape(B, 50) for k in R.VARS + ['counter', 'load_in', 'intervention', 'timestamp']}
    ok = np.r_[False, (X['counter'][1:, 0] - X['counter'][:-1, -1]) % 100000 == 1] & (X['intervention'][:, 0] == 1) & (X['intervention'][:, 5] == 1)
    ok &= np.r_[False, X['load_in'][:-1, -1] < 0.5] & (X['load_in'][:, 0] < 0.5) & (X['load_in'][:, 5] > 0.5)
    idx = np.nonzero(ok)[0]; n = len(idx); thr = norm.ppf(1 - 0.05 / (2 * len(R.VARS)))
    tau = np.hstack([X['timestamp'][idx - 1, -2:], X['timestamp'][idx]])
    name = os.path.basename(f); out[name] = {}
    for w in [4, 8]:
        z = {}
        for k in R.VARS:
            ext = np.hstack([X[k][idx - 1, -2:], X[k][idx]])
            d = resp(ext, tau, 7, w) - resp(ext, tau, 2, w)
            z[k] = round(float(d.mean() / (d.std(ddof=1) / np.sqrt(n))), 2)
        out[name][f'w={w}'] = dict(n=n, z=z, missed=[k for k in R.DESC if abs(z[k]) <= thr], false=[k for k in R.NONDESC if abs(z[k]) > thr])
        print(name, w, 'n', n, 'missed', out[name][f'w={w}']['missed'], 'false', out[name][f'w={w}']['false'],
              {k: z[k] for k in ['current_out', 'pressure_downwind'] + R.NONDESC})
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R1c.json'), 'w'), indent=1)
