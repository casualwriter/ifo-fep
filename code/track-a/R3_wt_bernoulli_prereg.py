"""R3 (PRE-REGISTERED, gate G1): hardware redo of the push response test on NEW data. Plan: plan/20261016-g1-hardware-redo-prereg.md.

Data: Causal Chambers wt_bernoulli_v1: random_loads_intake.csv (intake-fan load set to a random level at 2000 interventions, one every
100 rows; exhaust load fixed) and random_loads_both.csv (intake and exhaust loads random). About 8.3 s gap before each intervention row;
the system reaches its new steady state during the gap.
Differences from R1 (stated before running): level shift ~8-9 s after the change, not an initial slope; dose-response instead of a sham
group (intercept = response at zero push = sham baseline); identical timing for all events; random push sizes, so slow drifts are
unrelated to the push size.
Analysis: events = intervention rows (skip the first, and any with non-contiguous counters within 5 rows). d_in = load_in[t] - load_in[t-1]
(d_out likewise). Response r = mean(v[t..t+4]) - mean(v[t-5..t-1]). OLS r = a + b_in d_in (+ b_out d_out), HC1 robust SE, z = b/SE.
Bonferroni within file (11 tests intake file, 22 both file), alpha 0.05. Ground truth: each load has 8 sensor descendants and 3
non-descendants (pressure_ambient, signal_1, signal_2).
Predictions: (1) 0 of 9 actuator x non-descendant tests significant; (2) >= 7/8 descendants detected per actuator-file combination.
Death clauses: (1) any non-descendant significant -> calibration claim fails on this hardware; (2) < 6/8 descendants in any combination
-> power insufficient at this setting.
AMENDMENT (before any results; first run stopped on a singular design): in random_loads_both.csv both loads are identical in every row,
so the both file uses one joint push d = change of the common load; 11 tests per file; 2 combinations, 6 non-descendant tests."""
import os, json
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'data', 'causalchamber', 'wt_bernoulli_v1')
DESC = ['rpm_in', 'rpm_out', 'current_in', 'current_out', 'pressure_intake', 'pressure_upwind', 'pressure_downwind', 'mic']
NONDESC = ['pressure_ambient', 'signal_1', 'signal_2']
SENS = DESC + NONDESC
def ols_hc1(X, y):
    X = np.c_[np.ones(len(y)), X]; XtXi = np.linalg.inv(X.T @ X); b = XtXi @ X.T @ y; e = y - X @ b
    n, k = X.shape; V = XtXi @ (X.T * e ** 2) @ X @ XtXi * n / (n - k)
    return b[1:], np.sqrt(np.diag(V))[1:]
def analyse(fname, acts):
    df = pd.read_csv(os.path.join(DATA, fname)); c = df['counter'].to_numpy()
    iv = np.nonzero(df['intervention'].to_numpy() == 1)[0]
    iv = [t for t in iv[1:] if t - 5 >= 0 and t + 5 <= len(df) and np.all(np.diff(c[t - 5:t + 5]) % 100000 == 1)]
    iv = np.array(iv)
    D = np.column_stack([df[a].to_numpy()[iv] - df[a].to_numpy()[iv - 1] for a in acts])
    m = len(acts) * len(SENS); thr = norm.ppf(1 - 0.05 / (2 * m)); res = {}
    for s in SENS:
        v = df[s].to_numpy(float)
        r = np.array([v[t:t + 5].mean() - v[t - 5:t].mean() for t in iv])
        b, se = ols_hc1(D, r)
        for j, a in enumerate(acts):
            res[f'{a}->{s}'] = dict(b=float(b[j]), z=round(float(b[j] / se[j]), 2), significant=bool(abs(b[j] / se[j]) > thr), descendant=s in DESC)
    return dict(events=len(iv), z_threshold=round(float(thr), 2), tests=res)
if __name__ == "__main__":
    out = {'random_loads_intake.csv': analyse('random_loads_intake.csv', ['load_in']),
           'random_loads_both.csv': analyse('random_loads_both.csv', ['load_in'])}   # joint push: load_in == load_out
    summary = {}
    for f, r in out.items():
        for a in ['load_in']:
            t = {k.split('->')[1]: v for k, v in r['tests'].items() if k.startswith(a + '->')}
            det = [s for s in DESC if t[s]['significant']]; fal = [s for s in NONDESC if t[s]['significant']]
            summary[('intake push' if 'intake' in f else 'joint push') + f' ({f})'] = dict(detected=f'{len(det)}/8', missed=[s for s in DESC if s not in det], false=fal,
                                         z={s: t[s]['z'] for s in SENS})
        print(f, 'events', r['events'], 'threshold', r['z_threshold'])
    n_false = sum(len(v['false']) for v in summary.values()); min_det = min(int(v['detected'][0]) for v in summary.values())
    verdict = dict(false_calls=f'{n_false}/6', min_detected=min_det,
                   death_clause_1='FIRED' if n_false > 0 else 'not fired', death_clause_2='FIRED' if min_det < 6 else 'not fired',
                   prediction_1='met' if n_false == 0 else 'failed', prediction_2='met' if min_det >= 7 else 'failed')
    for k, v in summary.items(): print(' ', k, v)
    print('VERDICT', verdict)
    json.dump(dict(results=out, summary=summary, verdict=verdict), open(os.path.join(HERE, 'R3.json'), 'w'), indent=1)
