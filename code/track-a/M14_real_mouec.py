"""M14: real case. Is the sparse-coupling (wiring) blanket of a brain system a blanket for an observer who reads it continuously?
Prospectively specified (written before running).

Data: Gilson et al. (2018, NeuroImage 180:534) via mb-BCA/notebooks_review2019 (Apache-2.0): BOLD, 22 subjects x 5 runs
(2 rest, 3 movie) x 66 Desikan ROIs x 300 TR (TR = 2 s); SC_anat (DTI). Fit: pyMOU (mb-BCA/pyMOU f0626be, unchanged) with the
authors' notebook settings: band-pass 0.01-0.2 Hz (Butterworth order 3, filtfilt), centre, rescale FC, mask_EC = SC > 0 plus
homotopic links, diagonal Sigma, algo 'NeNe2020', eta_C = 0.005, eta_Sigma = 0.01. pyMOU's J has J_ij = i -> j, so A = J.T.

Partitions (fixed in advance): I in {R-visual: rCUN rLING rPCAL rLOCC; L-visual: lCUN lLING lPCAL lLOCC;
R-sensorimotor: rPREC rPSTC rPARC; L-sensorimotor: lPREC lPSTC lPARC}. B = all ROIs linked to I in mask_EC (either direction);
E = the rest. Then the fitted A has no direct I-E term (wiring cut), so the wiring observer reports an exact blanket.
Observer: start known, blanket read continuously, push sigma_u = 0.1, h = 2 TR. Score: residual leak (Riccati, blanket_report.cont_leak),
worst direction.
Null (truth M = 0 by construction): for each fit and I set, split the blanket: each b in B keeps inputs only from the side (I or E)
with the larger summed |A|; simulate the split model with the fitted Sigma (Euler dt = 0.02 TR, burn-in 200 TR, 300 samples at 1 TR),
apply the same filter/centre/rescale/fit pipeline, and compute the same leak for the refitted model. Seed 1414.
Predictions: (P1) median real leak > 0.01; (P2) real leak > the null leak of the same (subject, run, I set) in >= 80% of cases.
Death: P2 < 50% -> on this data, the M leak cannot be told apart from estimation noise (inconclusive, not evidence for M = 0).
Exploratory: |M_hat| (Frobenius, both directions), leak rest vs movie."""
import os, sys, json, time, warnings
import numpy as np
import scipy.signal as spsg
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..', '..')
sys.path.insert(0, os.path.join(ROOT, 'code', 'pymou_f0626be')); sys.path.insert(0, os.path.join(ROOT, 'code', 'tools'))
from mou_model import MOU
from blanket_report import cont_leak, M_matrix
warnings.filterwarnings('ignore')
D_ = os.path.join(ROOT, 'data', 'mouec')
ts = np.load(os.path.join(D_, 'ts_emp.npy')); SC = np.load(os.path.join(D_, 'SC_anat.npy'))
labels = [s.strip() for s in np.load(os.path.join(D_, 'ROI_labels.npy'), allow_pickle=True)]
n_sub, n_run, N, T = ts.shape
b_f, a_f = spsg.iirfilter(3, [0.01 / 0.25, 0.2 / 0.25], btype='bandpass', ftype='butter')
mask = SC > 0
for i in range(N): mask[i, i] = False; mask[i, N - 1 - i] = True
mask_S = np.eye(N, dtype=bool)
def prep(X):                                    # X: N x T -> filtered, centred
    Y = spsg.filtfilt(b_f, a_f, X, axis=-1); return Y - Y.mean(1, keepdims=True)
# rescale factor as in the notebook (from rest runs, zero-lag variance)
Yall = np.array([[prep(ts[s, r]) for r in range(n_run)] for s in range(n_sub)])
var0 = np.mean([np.mean(np.var(Yall[s, r], axis=1)) for s in range(n_sub) for r in range(2)])
scale = np.sqrt(0.5 / var0)
def fit(Y):
    m = MOU(); m.fit((Y * scale).T, mask_C=mask, mask_Sigma=mask_S, algo_version="NeNe2020", eta_C=0.005, eta_Sigma=0.01)
    return m.J.T.copy(), m.Sigma.copy(), m.d_fit['correlation'] ** 2
SETS = {'R-visual': ['rCUN', 'rLING', 'rPCAL', 'rLOCC'], 'L-visual': ['lCUN', 'lLING', 'lPCAL', 'lLOCC'],
        'R-sensorimotor': ['rPREC', 'rPSTC', 'rPARC'], 'L-sensorimotor': ['lPREC', 'lPSTC', 'lPARC']}
PART = {}
for k, names in SETS.items():
    I = [labels.index(x) for x in names]
    nb = (mask[:, I].any(1) | mask[I, :].any(0)); nb[I] = False
    B = list(np.nonzero(nb)[0]); E = [j for j in range(N) if j not in I and j not in B]
    PART[k] = (I, B, E)
def split(A, I, B, E):
    A2 = A.copy()
    for b in B:
        if np.abs(A[b, I]).sum() >= np.abs(A[b, E]).sum(): A2[b, E] = 0.0
        else: A2[b, I] = 0.0
    m = np.linalg.eigvals(A2).real.max()
    if m > -0.01: A2 -= (m + 0.01) * np.eye(N)
    return A2
rng = np.random.default_rng(1414)
def simulate(A, S, n=300, dt=0.02, burn=200):
    L = np.linalg.cholesky(S + 1e-12 * np.eye(N)); x = np.zeros(N); out = np.empty((N, n)); k = int(round(1 / dt))
    for _ in range(int(burn / dt)): x = x + A @ x * dt + np.sqrt(dt) * L @ rng.standard_normal(N)
    for t in range(n):
        for _ in range(k): x = x + A @ x * dt + np.sqrt(dt) * L @ rng.standard_normal(N)
        out[:, t] = x
    return out
rows = []; t0 = time.time()
for s in range(n_sub):
    for r in range(n_run):
        A, S, R2 = fit(Yall[s, r])
        for k, (I, B, E) in PART.items():
            real = cont_leak(A, S, I, B, E, 2.0)
            Mn = max(np.linalg.norm(M_matrix(A, S, I, B, E)), np.linalg.norm(M_matrix(A, S, E, B, I)))
            As = split(A, I, B, E); Ysim = prep(simulate(As, S)) / scale
            An, Sn, _ = fit(Ysim)
            null = cont_leak(An, Sn, I, B, E, 2.0)
            rows.append(dict(sub=s, run=r, rest=r < 2, set=k, R2=float(R2), real=float(real), null=float(null), M=float(Mn),
                             nB=len(B), nE=len(E)))
        print(s, r, round(R2, 3), [(x['set'][:5], round(x['real'], 4), round(x['null'], 4)) for x in rows[-4:]], round(time.time() - t0), 's', flush=True)
real = np.array([x['real'] for x in rows]); null = np.array([x['null'] for x in rows]); ok = ~np.isnan(real) & ~np.isnan(null)
out = dict(n=int(ok.sum()), P1_median_real=float(np.median(real[ok])), median_null=float(np.median(null[ok])),
           P2_frac_real_gt_null=float(np.mean(real[ok] > null[ok])),
           by_set={k: dict(median_real=float(np.nanmedian([x['real'] for x in rows if x['set'] == k])),
                           median_null=float(np.nanmedian([x['null'] for x in rows if x['set'] == k])),
                           frac_gt=float(np.mean([x['real'] > x['null'] for x in rows if x['set'] == k])),
                           nB=PART[k][1].__len__(), nE=len(PART[k][2])) for k in PART},
           rest_vs_movie=dict(rest=float(np.nanmedian([x['real'] for x in rows if x['rest']])), movie=float(np.nanmedian([x['real'] for x in rows if not x['rest']]))),
           median_R2=float(np.median([x['R2'] for x in rows])))
print(json.dumps(out, indent=1))
json.dump(dict(summary=out, rows=rows, partitions={k: [list(map(int, v)) for v in p] for k, p in PART.items()}),
          open(os.path.join(HERE, 'M14.json'), 'w'), indent=1)
