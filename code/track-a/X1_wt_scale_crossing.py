"""X1 (EXPLORATORY, not pre-registered): does the outside air pressure cross into the wind tunnel, and at which time scale?
Same data as R1 (wt_intake_impulse_v1). Changes of pressure_ambient vs changes of each tunnel barometer, over time scales from
one sample (0.15 s) to many blocks (hours). Within a block: quiet rows 27-49 only (no push phase). Across blocks: mean of rows 30-49
of each block, differenced over m blocks (one block ~ 8.6 s incl. the ~0.9 s gap). Reported: correlation of changes and the
regression slope (tunnel change per unit outside change; 1 = full crossing)."""
import os, glob, json
import numpy as np, pandas as pd
import R1_wt_impulse as R
BAROS = ['pressure_intake', 'pressure_upwind', 'pressure_downwind']
def stats(a, b):
    a, b = a.ravel(), b.ravel(); c = np.corrcoef(a, b)[0, 1]; s = np.cov(a, b)[0, 1] / a.var()
    return round(float(c), 3), round(float(s), 3)
out = {}
for f in sorted(glob.glob(os.path.join(R.DATA, '*.csv'))):
    df = pd.read_csv(f); B = df['flag'].nunique(); name = os.path.basename(f)
    X = {k: df[k].to_numpy(float).reshape(B, 50) for k in BAROS + ['pressure_ambient', 'timestamp']}
    blk_dt = float(np.median(np.diff(X['timestamp'][:, 30])))
    res = {}
    for L in [1, 5, 20]:                                         # within quiet rows
        t = np.arange(27, 50 - L)
        da = X['pressure_ambient'][:, t + L] - X['pressure_ambient'][:, t]
        res[f'{L * 0.15:.2f}s'] = {k: stats(da, X[k][:, t + L] - X[k][:, t]) for k in BAROS}
    M = {k: X[k][:, 30:50].mean(1) for k in BAROS + ['pressure_ambient']}
    for m in [1, 10, 100, 1000]:                                 # across blocks
        da = M['pressure_ambient'][m:] - M['pressure_ambient'][:-m]
        res[f'{m * blk_dt / 60:.1f}min'] = {k: stats(da, M[k][m:] - M[k][:-m]) for k in BAROS}
    out[name] = res
    print(name, 'block', round(blk_dt, 2), 's')
    for sc, v in res.items(): print('  ', sc.rjust(8), {k.replace('pressure_', ''): v[k] for k in BAROS})
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'X1.json'), 'w'), indent=1)
