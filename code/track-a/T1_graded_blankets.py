"""T1: graded blankets. Can push + sham measure the crossing share rho of several coexisting candidate blankets?
Design and death clauses: plan/20261016-manifest-probability-crossing-share-v1.md, section 5 (fixed before running).

crossing share for a candidate (I, B, E), J = direct couplings (off-diagonal of A):
  rho_IE = sum_{i in I, b in B} J_bi^2 / (sum J_bi^2 + sum_{i in I, e in E} J_ei^2), rho_EI the same with I and E swapped,
  rho = min(rho_IE, rho_EI). rho = 1: no direct I-E coupling.
Setup: 12 observed nodes + 1 hidden slow drive (c in {0, 0.4}); edge probability 0.25, weights U(0.2, 0.8) with random sign;
network seed 3030; 50 kept networks per c. In each network, enumerate all 4/4/4 partitions and pick one candidate with true rho in
low [0.20, 0.40], mid [0.65, 0.85], high [0.95, 1); a control with rho = 1 if one exists. Networks without all three are skipped (counted).
Estimators: (a) push + sham: push every node (u = 3, K = 1000, sham u = 0), initial slope a = (4 r(h1) - r(2 h1)) / (2 h1 u),
  h1 = 0.05; J^2 estimate = max(a^2 - var(a), 0), a = mean push slope - mean sham slope. (b) observational: Jacobian from data,
  logm(F)/h with F from regressing x_{t+h} on x_t (N = 10000, h = 0.05), then J^2.
Prediction: push + sham ranks the three candidates correctly in >= 90% of networks at c = 0 and 0.4; MAE <= 0.10; control rho_hat >= 0.95
  in >= 90%. Observational: similar at c = 0; rank-correct drops by >= 20 points at c = 0.4.
Death clauses: (1) push + sham rank-correct < 70% at c = 0; (2) control rho_hat >= 0.95 in < 80%; (3) observational within 5 points of
  push + sham at c = 0.4."""
import os, sys, json, time, itertools
import numpy as np
from scipy.linalg import logm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'causal-partition-v0_6'))
from stats import wilson
n, dt, U, K, H1 = 12, 0.01, 3.0, 1000, 0.05
BINS = {'low': (0.20, 0.40), 'mid': (0.65, 0.85), 'high': (0.95, 1.0)}
I_IDX = list(itertools.combinations(range(n), 4))
PARTS = [(I, E, tuple(sorted(set(range(n)) - set(I) - set(E))))
         for I in I_IDX for E in itertools.combinations(sorted(set(range(n)) - set(I)), 4) if I[0] < E[0]]   # I/E unordered
MASK = {k: np.zeros((len(PARTS), n)) for k in 'IEB'}
for p, (I, E, B) in enumerate(PARTS):
    MASK['I'][p, list(I)] = 1; MASK['E'][p, list(E)] = 1; MASK['B'][p, list(B)] = 1
def rho_all(J2):
    """rho for every partition; J2[j, i] = squared coupling i -> j (diagonal ignored). NaN if a side has no outgoing coupling."""
    J2 = J2 * (1 - np.eye(n))
    s = lambda a, b: np.einsum('pt,ts,ps->p', MASK[a], J2, MASK[b])     # sum of J2[target in a, source in b]
    with np.errstate(invalid='ignore', divide='ignore'):
        r_ie = s('B', 'I') / (s('B', 'I') + s('E', 'I')); r_ei = s('B', 'E') / (s('B', 'E') + s('I', 'E'))
    return np.minimum(r_ie, r_ei)
def net(rng, c):
    A = -np.eye(n + 1); A[n, n] = -0.2
    M = (rng.random((n, n)) < 0.25) & ~np.eye(n, dtype=bool)
    A[:n, :n] += M * rng.uniform(0.2, 0.8, (n, n)) * rng.choice([-1, 1], (n, n))
    A[:n, n] = c * rng.uniform(0.5, 1.5, n)
    m = np.linalg.eigvals(A).real.max()
    if m > -0.1: A[:n, :n] -= (m + 0.1) * np.eye(n)
    return A
def simulate(A, rng, nsamp=10000, every=5):
    x = np.zeros(n + 1); out = np.empty((nsamp, n)); sq = np.sqrt(dt)
    for _ in range(3000): x = x + A @ x * dt + rng.standard_normal(n + 1) * sq
    for k in range(nsamp):
        for _ in range(every): x = x + A @ x * dt + rng.standard_normal(n + 1) * sq
        out[k] = x[:n]
    return out
def slopes(A, X, j, u, rng, F1, F2):
    idx = rng.integers(0, len(X) - 2, K); x = np.zeros((K, n + 1)); x[:, :n] = X[idx]
    x[:, n] = rng.standard_normal(K) * np.sqrt(0.5 / abs(A[n, n]))
    x0 = x[:, :n].copy(); x[:, j] += u; sq = np.sqrt(dt)
    for s in range(10):
        x = x + x @ A.T * dt + rng.standard_normal((K, n + 1)) * sq
        if s == 4: r1 = x[:, :n] - x0 @ F1
    r2 = x[:, :n] - x0 @ F2
    return (4 * r1 - r2) / (2 * H1)                               # per trial, not yet divided by u
def J2_push(A, X, rng):
    F1 = np.linalg.lstsq(X[:-1], X[1:], rcond=None)[0]; F2 = np.linalg.lstsq(X[:-2], X[2:], rcond=None)[0]
    J2 = np.zeros((n, n))
    for j in range(n):
        p = slopes(A, X, j, U, rng, F1, F2); s = slopes(A, X, j, 0.0, rng, F1, F2)
        a = (p.mean(0) - s.mean(0)) / U; v = (p.var(0, ddof=1) / K + s.var(0, ddof=1) / K) / U ** 2
        J2[:, j] = np.maximum(a ** 2 - v, 0)
    return J2
def J2_obs(X, h=0.05):
    F = np.linalg.lstsq(X[:-1], X[1:], rcond=None)[0].T            # x_{t+h} = F x_t
    return np.real(logm(F) / h) ** 2
def run(c, G, seed=3030, max_tries=400):
    rng = np.random.default_rng(seed + int(10 * c)); rows = []; skipped = 0; tries = 0
    while len(rows) < G and tries < max_tries:
        tries += 1
        A = net(rng, c); J2 = A[:n, :n] ** 2; r = rho_all(J2)
        pick = {}
        for b, (lo, hi) in BINS.items():
            ok = np.nonzero((r >= lo) & (r < hi))[0]
            if len(ok): pick[b] = int(rng.choice(ok))
        ctrl = np.nonzero(r == 1)[0]
        if len(pick) < 3: skipped += 1; continue
        if len(ctrl): pick['control'] = int(rng.choice(ctrl))
        X = simulate(A, rng)
        rp, ro = rho_all(J2_push(A, X, rng)), rho_all(J2_obs(X))
        rows.append({b: dict(true=float(r[p]), push=float(rp[p]), obs=float(ro[p])) for b, p in pick.items()})
    return rows, skipped
def summary(rows):
    out = {}
    for m in ['push', 'obs']:
        rank = [r['low'][m] < r['mid'][m] < r['high'][m] for r in rows]
        mae = [abs(r[b][m] - r[b]['true']) for r in rows for b in BINS]
        ctl = [r['control'][m] >= 0.95 for r in rows if 'control' in r]
        out[m] = dict(rank_correct=wilson(int(sum(rank)), len(rank)), mae=round(float(np.nanmean(mae)), 3),
                      control_ge_095=wilson(int(sum(ctl)), len(ctl)) if ctl else None, n_control=len(ctl))
    return out
if __name__ == "__main__":
    G = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    out = {}; t0 = time.time()
    for c in [0.0, 0.4]:
        rows, skipped = run(c, G)
        out[f'c={c}'] = dict(summary=summary(rows), kept=len(rows), skipped=skipped, rows=rows)
        print(c, 'kept', len(rows), 'skipped', skipped, json.dumps(out[f'c={c}']['summary'], default=float), round(time.time() - t0), 's', flush=True)
    json.dump(out, open(os.path.join(HERE, f'T1_G{G}.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
