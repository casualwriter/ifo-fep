"""R3b (POST-HOC, after R3's death clause 1 fired): is the outside barometer's response to fan pushes a real slow physical effect?
Same data and model as R3. Checks: (1) effect size of the outside barometer vs the tunnel barometers (Pa per unit load);
(2) response windows further after the change (rows t+5..t+9, t+20..t+24, t+50..t+54, t+90..t+94; ~1-15 s after the first
post-change row); a real physical effect should persist, a fluke should not; (3) placebo: the same regression with d shifted by one
event (d of the next event), which must give no effect."""
import os, json
import numpy as np, pandas as pd
import R3_wt_bernoulli_prereg as R
out = {}
for f in ['random_loads_intake.csv', 'random_loads_both.csv']:
    df = pd.read_csv(os.path.join(R.DATA, f)); c = df['counter'].to_numpy()
    iv = np.nonzero(df['intervention'].to_numpy() == 1)[0]
    iv = np.array([t for t in iv[1:] if t + 95 <= len(df) and np.all(np.diff(c[t - 5:t + 5]) % 100000 == 1)])
    d = (df['load_in'].to_numpy()[iv] - df['load_in'].to_numpy()[iv - 1])[:, None]
    res = {}
    for s in ['pressure_ambient', 'pressure_upwind', 'pressure_downwind', 'pressure_intake']:
        v = df[s].to_numpy(float); base = np.array([v[t - 5:t].mean() for t in iv]); row = {}
        for k in [0, 5, 20, 50, 90]:
            r = np.array([v[t + k:t + k + 5].mean() for t in iv]) - base
            b, se = R.ols_hc1(d, r); row[f't+{k}'] = (round(float(b[0]), 2), round(float(b[0] / se[0]), 2))
        r0 = np.array([v[t:t + 5].mean() for t in iv]) - base
        b, se = R.ols_hc1(np.roll(d, -1)[:-1], r0[:-1]); row['placebo_next_d'] = (round(float(b[0]), 2), round(float(b[0] / se[0]), 2))
        res[s] = row
    out[f] = res
    print(f)
    for s, row in res.items(): print('  ', s.ljust(18), row)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R3b.json'), 'w'), indent=1)
