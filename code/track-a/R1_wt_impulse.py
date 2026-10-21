"""R1: push + sham test on real hardware (Causal Chambers wind tunnel, dataset wt_intake_impulse_v1, Gamella et al. 2025).

Design found in the data: each 50-row block (column 'flag') sends a command at row 0 that sets load_in to 0.01 while
load_in is already 0.01 (a real sham: command sent, no change), and at row 5 a push load_in 0.01 -> 1.0.
Ground truth (causalchamber.ground_truth, wind tunnel 'standard'): load_in has 8 sensor descendants
(rpm_in, rpm_out, current_in, current_out, pressure_intake, pressure_upwind, pressure_downwind, mic)
and 3 sensor non-descendants (pressure_ambient, signal_1, signal_2). All barometers share a hidden common cause
(atmospheric pressure), which drives pressure_ambient and the tunnel barometers together.

Fixed before running:
 - response of variable v to an event at row t: r = mean(v[t+1..t+4]) - (v[t] + 2.5 s), s = (v[t] - v[t-2]) / 2
   (deviation from linear extrapolation of the last 3 samples). Sham: t = row 0 of the block (rows t-2, t-1 come
   from the previous block; blocks used only if contiguous in 'counter'). Push: t = row 5.
 - per block d = r_push - r_sham; z = mean(d) / (sd(d) / sqrt(n)); two-sided, Bonferroni over the 11 sensors, alpha 0.05.
 - primary file: load_out_0.5_osr_downwind_4.csv; the other 4 files are replications.
 - observational contrast: correlation of 1-step increments and of levels between pressure_ambient and
   pressure_downwind, outside push windows (rows 30-49 of each block).

Prediction: all 8 descendants respond (significant); none of the 3 non-descendants does. The observational level
correlation between pressure_ambient and pressure_downwind is high (shared atmosphere), so an observational
criterion would link two sensors that the push separates.
Death clauses: (1) any non-descendant significant in the primary file -> the test is not calibrated on this hardware.
 (2) more than one descendant missed in the primary file -> the test lacks power at this setting."""
import os, sys, json, glob
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'data', 'causalchamber', 'wt_intake_impulse_v1')
DESC = ['rpm_in', 'rpm_out', 'current_in', 'current_out', 'pressure_intake', 'pressure_upwind', 'pressure_downwind', 'mic']
NONDESC = ['pressure_ambient', 'signal_1', 'signal_2']
VARS = DESC + NONDESC
def responses(v, t):
    s = (v[:, t] - v[:, t - 2]) / 2
    return v[:, t + 1:t + 5].mean(1) - (v[:, t] + 2.5 * s)
def analyse(path):
    df = pd.read_csv(path)
    B = df['flag'].nunique(); X = {k: df[k].to_numpy().reshape(B, 50) for k in VARS + ['counter', 'load_in', 'intervention']}
    # prepend the last 2 rows of the previous block, so the sham at row 0 has t-2, t-1
    ok = np.r_[False, X['counter'][1:, 0] == X['counter'][:-1, -1] + 1]
    ok &= (X['intervention'][:, 0] == 1) & (X['intervention'][:, 5] == 1)
    ok &= np.r_[False, X['load_in'][:-1, -1] < 0.5] & (X['load_in'][:, 0] < 0.5) & (X['load_in'][:, 5] > 0.5)
    idx = np.nonzero(ok)[0]; n = len(idx); thr = norm.ppf(1 - 0.05 / (2 * len(VARS)))
    res = {}
    for k in VARS:
        ext = np.hstack([X[k][idx - 1, -2:], X[k][idx]])                    # columns: prev48, prev49, row0..row49
        d = responses(ext, 7) - responses(ext, 2)                           # push at row 5 (=col 7), sham at row 0 (=col 2)
        z = float(d.mean() / (d.std(ddof=1) / np.sqrt(n)))
        res[k] = dict(z=round(z, 2), significant=bool(abs(z) > thr), descendant=k in DESC)
    quiet = slice(30, 50)
    a, b = X['pressure_ambient'][idx, quiet], X['pressure_downwind'][idx, quiet]
    obs = dict(level_corr=round(float(np.corrcoef(a.ravel(), b.ravel())[0, 1]), 3),
               increment_corr=round(float(np.corrcoef(np.diff(a, axis=1).ravel(), np.diff(b, axis=1).ravel())[0, 1]), 3))
    miss = [k for k in DESC if not res[k]['significant']]; false = [k for k in NONDESC if res[k]['significant']]
    return dict(blocks_used=n, z_threshold=round(float(thr), 2), sensors=res, missed_descendants=miss,
                false_rejections=false, observational_ambient_vs_downwind=obs)
if __name__ == "__main__":
    out = {}
    for f in sorted(glob.glob(os.path.join(DATA, '*.csv'))):
        name = os.path.basename(f); out[name] = analyse(f)
        r = out[name]
        print(name, 'n', r['blocks_used'], 'missed', r['missed_descendants'], 'false', r['false_rejections'], r['observational_ambient_vs_downwind'])
        print('   ', {k: v['z'] for k, v in r['sensors'].items()})
    json.dump(out, open(os.path.join(HERE, 'R1.json'), 'w'), indent=1)
