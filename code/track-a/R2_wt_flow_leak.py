"""R2: does the observational flow test link sensors that only share a hidden slow driver? Real data, same files as R1.
pressure_ambient and the tunnel barometers (intake, upwind, downwind) share atmospheric pressure (hidden common cause);
the ground truth has no causal path between them, and R1 shows that pushes on the intake fan do not move pressure_ambient.

Fixed before running:
 - series: only quiet, evenly spaced rows: rows 27-49 of each 50-row block (after the push ends at row 25; the step
   into row 0 of a block takes ~0.93 s, so blocks are not joined). Spacing 0.15 s. (Changed before the first complete
   run: the first attempt stopped on a check and showed the uneven spacing.)
 - flow test (as F3.flow_test): Gaussian conditional MI, CMI(x_{t+L} ; y_t | x_t, x_{t-1}) and the reverse direction,
   likelihood-ratio 2N*CMI ~ chi2(1); "linked" if min p < 0.05/2. Lags L in {1, 5}.
 - pairs: pressure_ambient vs each tunnel barometer (shared hidden cause, no causal path).
   Negative control: pressure_ambient vs signal_1 and vs signal_2 (no shared cause, no path).
   Positive control: rpm_in vs pressure_upwind (both driven by load_in; airflow links them).
Prediction: the flow test declares ambient-barometer pairs linked (leak from the shared slow driver, as in paper result 4),
 does not link the negative controls, and links the positive control.
Death clause: if ambient-barometer pairs are NOT linked at both lags in the primary file, the real-data leak claim fails
 (the shared atmosphere is too slow or too weak to fool the flow test here)."""
import os, glob, json
import numpy as np, pandas as pd
from scipy.stats import chi2
import R1_wt_impulse as R
def cmi(C, i, e, b):
    ld = lambda idx: np.linalg.slogdet(C[np.ix_(idx, idx)])[1]
    return 0.5 * (ld(i + b) + ld(e + b) - ld(b) - ld(i + e + b))
def flow(x, y, L):
    """x, y: (blocks, 50) arrays. Uses t with t-1 >= 27 and t+L <= 49 inside each block."""
    x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std()
    t = np.arange(28, 50 - L)
    def one(a, b):                                  # does b_t predict a_{t+L} given a_t, a_{t-1}?
        Z = np.c_[a[:, t + L].ravel(), b[:, t].ravel(), a[:, t].ravel(), a[:, t - 1].ravel()]
        C = np.cov(Z.T); N = len(Z); m = cmi(C, [0], [1], [2, 3])
        return float(chi2.sf(2 * N * m, 1)), float(m)
    p1, m1 = one(x, y); p2, m2 = one(y, x)
    return dict(p_min=min(p1, p2), cmi=[round(m1, 6), round(m2, 6)], linked=min(p1, p2) < 0.025)
PAIRS = [('pressure_ambient', 'pressure_intake', 'shared'), ('pressure_ambient', 'pressure_upwind', 'shared'),
         ('pressure_ambient', 'pressure_downwind', 'shared'), ('pressure_ambient', 'signal_1', 'negative'),
         ('pressure_ambient', 'signal_2', 'negative'), ('rpm_in', 'pressure_upwind', 'positive')]
if __name__ == "__main__":
    out = {}
    for f in sorted(glob.glob(os.path.join(R.DATA, '*.csv'))):
        df = pd.read_csv(f); B = df['flag'].nunique()
        name = os.path.basename(f); out[name] = {}
        for a, b, kind in PAIRS:
            for L in [1, 5]:
                r = flow(df[a].to_numpy(float).reshape(B, 50), df[b].to_numpy(float).reshape(B, 50), L)
                out[name][f'{a}|{b}|L={L}'] = dict(r, kind=kind)
        print(name, {k.replace('pressure_', 'p_'): (v['kind'][:3], v['linked'], f"{v['p_min']:.1e}") for k, v in out[name].items()})
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R2.json'), 'w'), indent=1)
