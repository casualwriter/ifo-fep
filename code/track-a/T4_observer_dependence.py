"""T4 (PRE-REGISTERED): the blanket comes from the observer. Does the best blanket change with the observer, and are blankets with
internal positive feedback more observer-invariant? Plan: plan/20261016-t4-observer-dependence-prereg.md (fixed before running).

12-node linear networks (T1 generator, no hidden drive), 200 networks, seed 5050, exact. Candidates: all 4/4/4 partitions.
Observers: O1 direct couplings (h -> 0), all nodes; O2 response expm(A*1); O3 response expm(A*3); O4 like O1 but 2 random nodes can be
neither pushed nor read. rho per observer as in T1 (squared response shares, min over directions). Agreement range = max - min of
rho over the 4 observers. Internal positive feedback: a directed cycle (length 2-4) inside I with positive weight product.
Predictions: (1) best blanket under O1 and O3 differs in >= 30% of networks; (2) among candidates with rho_O1 >= 0.5, median range is
smaller with internal positive feedback (one-sided Mann-Whitney p < 0.05).
Death clauses: (1) best blanket same in >= 90% of networks; (2) prediction 2 fails."""
import os, sys, json, itertools
import numpy as np
from scipy.linalg import expm
from scipy.stats import mannwhitneyu
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from T1_graded_blankets import n, net, MASK, PARTS
def rho_from(R2, keep):
    """R2[j, i] = squared response of j to a push on i; nodes outside `keep` are neither pushed nor read."""
    k = np.zeros(n); k[list(keep)] = 1
    R2 = R2 * (1 - np.eye(n)) * np.outer(k, k)
    s = lambda a, b: np.einsum('pt,ts,ps->p', MASK[a], R2, MASK[b])
    with np.errstate(invalid='ignore', divide='ignore'):
        r1 = s('B', 'I') / (s('B', 'I') + s('E', 'I')); r2 = s('B', 'E') / (s('B', 'E') + s('I', 'E'))
    return np.fmin(r1, r2)                      # a side with no visible outgoing influence is skipped (fmin ignores NaN)
def pos_cycle(W, nodes):
    for L in range(2, len(nodes) + 1):
        for cyc in itertools.permutations(nodes, L):
            if cyc[0] != min(cyc): continue
            w = np.prod([W[cyc[(i + 1) % L], cyc[i]] for i in range(L)])
            if w > 0: return True
    return False
if __name__ == "__main__":
    rng = np.random.default_rng(5050)
    differs, ranges_fb, ranges_no = [], [], []
    for g in range(200):
        A = net(rng, 0.0); Ao = A[:n, :n]; W = Ao * (1 - np.eye(n))
        hide = rng.choice(n, 2, replace=False); allk = range(n); part = [i for i in range(n) if i not in hide]
        rhos = {'O1': rho_from(W ** 2, allk), 'O2': rho_from(expm(A)[:n, :n] ** 2, allk),
                'O3': rho_from(expm(3 * A)[:n, :n] ** 2, allk), 'O4': rho_from(W ** 2, part)}
        R = np.vstack([rhos[k] for k in ['O1', 'O2', 'O3', 'O4']])
        ok = ~np.isnan(R).any(0)
        best1 = np.nanargmax(np.where(ok, R[0], -1)); best3 = np.nanargmax(np.where(ok, R[2], -1))
        differs.append(bool(best1 != best3))
        rng_ = R.max(0) - R.min(0)
        cand = np.nonzero(ok & (R[0] >= 0.5))[0]
        if len(cand) > 300: cand = rng.choice(cand, 300, replace=False)     # cap per network for run time
        for p in cand:
            I = list(PARTS[p][0])
            (ranges_fb if pos_cycle(W, I) else ranges_no).append(float(rng_[p]))
    frac = float(np.mean(differs))
    u = mannwhitneyu(ranges_fb, ranges_no, alternative='less')
    out = dict(best_differs_O1_O3=frac, n_fb=len(ranges_fb), n_no=len(ranges_no),
               median_range_fb=float(np.median(ranges_fb)), median_range_no=float(np.median(ranges_no)), p_one_sided=float(u.pvalue),
               prediction_1=frac >= 0.30, prediction_2=bool(u.pvalue < 0.05),
               death_clause_1='FIRED' if (1 - frac) >= 0.90 else 'not fired', death_clause_2='FIRED' if u.pvalue >= 0.05 else 'not fired')
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, 'T4.json'), 'w'), indent=1)
