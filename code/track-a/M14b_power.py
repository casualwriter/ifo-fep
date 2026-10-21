"""M14b (POST HOC, after M14's death clause fired: real leak > null in only 18%). How much data is needed to tell M != 0 from M = 0?
For subjects 0-9, run 0 (rest), using the fitted model A (M != 0, 'true') and its split version A_split (M = 0) as ground truths:
simulate T in {300, 1200, 4800} TR, run the same filter / fit pipeline, and compute the continuous-reader leak (h = 2 TR) for the four
partitions. Separation = fraction of (subject, partition) pairs with leak(refit of A) > leak(refit of A_split); also median leaks.
Expectation (stated before this run): separation rises with T, and the refit leak of A_split falls with T."""
import os, sys, json, time, warnings
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "M14_real_mouec.py"), encoding="utf-8").read().split("\nrows = []; t0")[0]
exec(src)
rng = np.random.default_rng(1415)
res = {T_: [] for T_ in (300, 1200, 4800)}; t0 = time.time()
for s in range(10):
    A, S, _ = fit(Yall[s, 0])
    for k, (I, B, E) in PART.items():
        As = split(A, I, B, E); truth = cont_leak(A, S, I, B, E, 2.0)
        for T_ in res:
            Yt = prep(simulate(A, S, n=T_)) / scale; Ys = prep(simulate(As, S, n=T_)) / scale
            At, St, _ = fit(Yt); An, Sn, _ = fit(Ys)
            res[T_].append(dict(sub=s, set=k, truth=float(truth), leak_true=float(cont_leak(At, St, I, B, E, 2.0)),
                                leak_split=float(cont_leak(An, Sn, I, B, E, 2.0))))
    print(s, round(time.time() - t0), 's', flush=True)
out = {}
for T_, v in res.items():
    lt = np.array([x['leak_true'] for x in v]); ls = np.array([x['leak_split'] for x in v])
    out[str(T_)] = dict(separation=float(np.mean(lt > ls)), median_leak_true_refit=float(np.median(lt)),
                        median_leak_split_refit=float(np.median(ls)), median_truth=float(np.median([x['truth'] for x in v])))
print(json.dumps(out, indent=1)); json.dump(dict(summary=out, rows=res), open(os.path.join(HERE, 'M14b.json'), 'w'), indent=1)
