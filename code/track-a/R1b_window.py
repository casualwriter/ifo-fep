"""R1b (POST-HOC, chosen after seeing R1): response window length. R1 used 4 samples after the event (0.6 s) and missed
current_out and pressure_downwind. Here: windows of 4 (R1), 8 and 16 samples. Question: are the misses slow responses,
and does calibration (no non-descendant significant) hold with longer windows? Same blocks and test as R1.
Note: the push at row 5 reverses at row 25, so windows up to 16 stay inside the push phase."""
import os, glob, json
import numpy as np
import R1_wt_impulse as R
from scipy.stats import norm
import pandas as pd
def responses(v, t, w):
    s = (v[:, t] - v[:, t - 2]) / 2
    return v[:, t + 1:t + 1 + w].mean(1) - (v[:, t] + (w + 1) / 2 * s)
out = {}
for f in sorted(glob.glob(os.path.join(R.DATA, '*.csv'))):
    df = pd.read_csv(f); B = df['flag'].nunique()
    X = {k: df[k].to_numpy().reshape(B, 50) for k in R.VARS + ['counter', 'load_in', 'intervention']}
    ok = np.r_[False, X['counter'][1:, 0] == X['counter'][:-1, -1] + 1] & (X['intervention'][:, 0] == 1) & (X['intervention'][:, 5] == 1)
    ok &= np.r_[False, X['load_in'][:-1, -1] < 0.5] & (X['load_in'][:, 0] < 0.5) & (X['load_in'][:, 5] > 0.5)
    idx = np.nonzero(ok)[0]; n = len(idx); thr = norm.ppf(1 - 0.05 / (2 * len(R.VARS)))
    name = os.path.basename(f); out[name] = {}
    for w in [4, 8, 16]:
        z = {}
        for k in R.VARS:
            ext = np.hstack([X[k][idx - 1, -2:], X[k][idx]])
            d = responses(ext, 7, w) - responses(ext, 2, w)
            z[k] = round(float(d.mean() / (d.std(ddof=1) / np.sqrt(n))), 2)
        out[name][f'w={w}'] = dict(z=z, missed=[k for k in R.DESC if abs(z[k]) <= thr], false=[k for k in R.NONDESC if abs(z[k]) > thr])
        print(name, w, 'missed', out[name][f'w={w}']['missed'], 'false', out[name][f'w={w}']['false'],
              {k: z[k] for k in ['current_out', 'pressure_downwind', 'pressure_ambient', 'signal_1', 'signal_2']})
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R1b.json'), 'w'), indent=1)
