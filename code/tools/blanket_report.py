"""blanket_report.py - how much is a proposed Markov blanket a blanket, for which observer?

A Markov blanket is treated here as a relation between a system and an observer, not as a property of the system alone.
For a linear stochastic system dx = A x dt + D^(1/2) dW and a proposed partition (I, B, E), this tool reports the blanket
as seen by several observer classes, and which observers agree.

Observer classes (each gives a leak score; 0 = exact blanket for that observer):
  wiring      : direct couplings only (h -> 0). Score = share of squared direct I-E coupling (noise units).
  total(h)    : pushes I (or E), does NOT read B, reads the other side after time h. Score = total-effect share.
  read(h, d)  : knows the start state, reads B every d, reads the other side at h. Score = conditional push-information ratio.
  cont(h)     : knows the start state, reads B continuously (d -> 0). Score = exact residual leak (Kalman-Bucy Riccati).
Structural check: M = A_BI' D_B^-1 A_BE (both directions). M = 0 is the condition under which every observer of class
cont(h) sees an exact blanket (path-level, start known; a wiring cut is also needed).

Defaults of report(): start state known; isotropic push sigma_u = 0.1 on all of I (and on all of E for the other direction), in the units
of A and D as given; time scales hs = (0.5, 2.0); read interval d = 0.1; tol = 1e-6 (a leak below tol is shown as 'exact', meaning
numerically zero; a theorem-certified zero needs a wiring cut and M = 0); a direction with zero denominator gives nan and is skipped.
No stability or stationarity check: give a stable A.

Usage:  python blanket_report.py            (runs a demo)
        from blanket_report import report;  report(A, D, I, B, E)
Needs numpy and scipy only. Units: give A and D in your own units; scores marked [units] change with units, others do not."""
import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov, solve_continuous_are
from scipy.integrate import solve_ivp

def _ld(X):
    return np.linalg.slogdet(X)[1]

def M_matrix(A, D, I, B, E):
    """Structural cross term for I -> E: A_BI' D_B^-1 A_BE (|I| x |E|)."""
    return A[np.ix_(B, I)].T @ np.linalg.solve(D[np.ix_(B, B)], A[np.ix_(B, E)])

def wiring_leak(A, I, B, E):
    """[units] share of squared direct coupling that bypasses B, worst direction."""
    def one(src, tgt):
        d = (A[np.ix_(B, src)] ** 2).sum() + (A[np.ix_(tgt, src)] ** 2).sum()
        return np.nan if d == 0 else (A[np.ix_(tgt, src)] ** 2).sum() / d
    return np.nanmax([one(I, E), one(E, I)])

def total_leak(A, I, B, E, h):
    """[units] total-effect share after time h (observer does not read B), worst direction."""
    R2 = expm(A * h) ** 2
    def one(src, tgt):
        a, b = R2[np.ix_(tgt, src)].sum(), R2[np.ix_(B, src)].sum()
        return np.nan if a + b == 0 else a / (a + b)
    return np.nanmax([one(I, E), one(E, I)])

def read_leak(A, D, I, B, E, h, d, sigma_u=0.1):
    """Conditional push-information ratio: start known, B read every d, other side read at h. Unit-free."""
    n = len(A); m = max(1, int(round(h / d))); ts = h * np.arange(1, m + 1) / m
    Sig = solve_continuous_lyapunov(A, -D); Phi = [expm(A * t) for t in ts]
    W = [Sig - P @ Sig @ P.T for P in Phi]; CN = np.zeros((n * m, n * m))
    for a in range(m):
        for b in range(a + 1):
            blk = expm(A * (ts[a] - ts[b])) @ W[b]
            CN[a*n:(a+1)*n, b*n:(b+1)*n] = blk; CN[b*n:(b+1)*n, a*n:(a+1)*n] = blk.T
    def one(src, tgt):
        P = sigma_u * np.eye(n)[:, src]; Mx = np.vstack([Ph @ P for Ph in Phi])
        z = [a * n + k for a in range(m) for k in B]; y = [(m - 1) * n + k for k in tgt]
        info = lambda S: 0.5 * (_ld(Mx[S] @ Mx[S].T + CN[np.ix_(S, S)]) - _ld(CN[np.ix_(S, S)]))
        tot = info(y)
        return np.nan if tot < 1e-14 else (info(y + z) - info(z)) / tot
    return np.nanmax([one(I, E), one(E, I)])

def cont_leak(A, D, I, B, E, h, sigma_u=0.1):
    """Exact residual leak for continuous reading of B (start known), via the Kalman-Bucy Riccati equation. Unit-free.
    Assumes D block-diagonal over (I, B, E) and a wiring cut; returns nan if the push has no route to the other side."""
    def one(src, tgt):
        p, q = len(src), len(tgt); k = 2 * p + q
        F = np.zeros((k, k)); F[p:2*p, p:2*p] = A[np.ix_(src, src)]; F[2*p:, 2*p:] = A[np.ix_(tgt, tgt)]
        Q = np.zeros((k, k)); Q[p:2*p, p:2*p] = D[np.ix_(src, src)]; Q[2*p:, 2*p:] = D[np.ix_(tgt, tgt)]
        H = np.hstack([np.zeros((len(B), p)), A[np.ix_(B, src)], A[np.ix_(B, tgt)]])
        G = H.T @ np.linalg.solve(D[np.ix_(B, B)], H)
        S0 = np.zeros((k, k)); S0[:p, :p] = np.eye(p); S0[:p, p:2*p] = S0[p:2*p, :p] = sigma_u * np.eye(p)
        S0[p:2*p, p:2*p] = sigma_u ** 2 * np.eye(p)
        f = lambda t, y: (lambda S: (F @ S + S @ F.T + Q - S @ G @ S).ravel())(y.reshape(k, k))
        S = solve_ivp(f, (0, h), S0.ravel(), rtol=1e-10, atol=1e-14, method='DOP853').y[:, -1].reshape(k, k)
        S = (S + S.T) / 2; Suu, SuE, SEE = S[:p, :p], S[:p, 2*p:], S[2*p:, 2*p:]
        num = 0.5 * (_ld(Suu) - _ld(Suu - SuE @ np.linalg.solve(SEE, SuE.T)))
        Sig = solve_continuous_lyapunov(A, -D); Ph = expm(A * h); Wt = (Sig - Ph @ Sig @ Ph.T)[np.ix_(tgt, tgt)]
        Sg = Ph[np.ix_(tgt, src)] * sigma_u; den = 0.5 * (_ld(Sg @ Sg.T + Wt) - _ld(Wt))
        return np.nan if den < 1e-14 else max(num, 0.0) / den
    return np.nanmax([one(I, E), one(E, I)])

def _filter(A, D, I, B, E):
    """Stationary Kalman filter for (x_I, x_E) given the past of x_B. Returns F, H, P (error covariance), K (gain)."""
    IE = list(I) + list(E); ix = np.ix_
    F, H = A[ix(IE, IE)], A[ix(B, IE)]
    P = solve_continuous_are(F.T, H.T, D[ix(IE, IE)], D[ix(B, B)]); P = (P + P.T) / 2
    return F, H, P, P @ H.T @ np.linalg.inv(D[ix(B, B)])

def causal_leak(A, D, I, B, E):
    """Causal reader (stationary, has read only the past of B): I(x_I(t); x_E(t) | past of B) in nats. Unit-free.
    With a wiring cut and D block-diagonal, it is 0 iff M = 0 (companion Proposition 1c)."""
    p = len(I); P = _filter(A, D, I, B, E)[2]
    return max(0.5 * (_ld(P[:p, :p]) + _ld(P[p:, p:]) - _ld(P)), 0.0)

def sync_r2(A, D, I, B, E, reading="history"):
    """Encoding check (FEP synchronisation map). R2 = share of the variance of mu = E[x_E | blanket] that is linear in
    eta = E[x_I | blanket]. A map mu = sigma(eta) exists iff R2 = 1. reading = 'snapshot' (current x_B) or 'history' (past of B)."""
    ix = np.ix_; n = len(A); r2 = lambda Cmm, Cme, Cee: float(np.trace(Cme @ np.linalg.pinv(Cee) @ Cme.T) / np.trace(Cmm))   # pinv: eta may be degenerate
    if reading == "snapshot":
        Sg = solve_continuous_lyapunov(A, -D); Kb = np.linalg.inv(Sg[ix(B, B)]); Cbb = Sg[ix(B, B)]
        WI, WE = Sg[ix(I, B)] @ Kb, Sg[ix(E, B)] @ Kb
        return r2(WE @ Cbb @ WE.T, WE @ Cbb @ WI.T, WI @ Cbb @ WI.T)
    IE = list(I) + list(E); k = len(IE); F, H, P, K = _filter(A, D, I, B, E)
    hat = list(range(n, n + k)); Z = np.zeros((n + k, n + k)); Z[:n, :n] = A
    Z[n:, n:] = F - K @ H                                   # d xhat = (F xhat + A_(IE)B x_B) dt + K (dx_B - (A_BB x_B + H xhat) dt)
    Z[ix(hat, B)] += A[ix(IE, B)] - K @ A[ix(B, B)]
    Z[n:, :n] += K @ A[B, :]
    N = np.zeros((n + k, n + k)); N[:n, :n] = D; N[:n, n:] = D[:, B] @ K.T; N[n:, :n] = N[:n, n:].T; N[n:, n:] = K @ D[ix(B, B)] @ K.T
    C = solve_continuous_lyapunov(Z, -N); hI, hE = hat[:len(I)], hat[len(I):]
    return r2(C[ix(hE, hE)], C[ix(hE, hI)], C[ix(hI, hI)])

def report(A, D, I, B, E, hs=(0.5, 2.0), d=0.1, tol=1e-6, sigma_u=0.1, show=True):
    """Leak of partition (I, B, E) for each observer class, and which observers see an exact blanket (leak < tol).
    Also reports encoding (synchronisation map R2), which is a separate question from screening."""
    A, D = np.asarray(A, float), np.asarray(D, float)
    rows = [("structure", "|M| (I->E, E->I)", max(np.abs(M_matrix(A, D, I, B, E)).max(), np.abs(M_matrix(A, D, E, B, I)).max())),
            ("wiring", "direct couplings [units]", wiring_leak(A, I, B, E))]
    for h in hs:
        rows.append((f"total(h={h})", "does not read B [units]", total_leak(A, I, B, E, h)))
        rows.append((f"read(h={h},d={d})", "reads B every d", read_leak(A, D, I, B, E, h, d, sigma_u)))
        rows.append((f"cont(h={h})", "reads B continuously", cont_leak(A, D, I, B, E, h, sigma_u)))
    rows.append(("causal", "reads past of B only [nats]", causal_leak(A, D, I, B, E)))
    enc = [("sync map, snapshot", "R2 (1 = map exists)", sync_r2(A, D, I, B, E, "snapshot")),
           ("sync map, history", "R2 (1 = map exists)", sync_r2(A, D, I, B, E, "history"))]
    if show:
        print(f"partition I={list(I)} B={list(B)} E={list(E)}")
        for name, what, v in rows:
            tag = "" if name == "structure" else ("  exact" if v < tol else "")
            print(f"  {name:<18} {what:<28} {v:10.3g}{tag}")
        agree = [r[0] for r in rows[1:] if r[2] < tol]
        print("  observers that see an exact blanket:", agree if agree else "none")
        print("  encoding (separate from screening):")
        for name, what, v in enc:
            print(f"  {name:<18} {what:<28} {v:10.3g}{'  map exists' if v > 1 - tol else ''}")
    return rows + enc

if __name__ == "__main__":
    rng = np.random.default_rng(1)
    I, S, Ac, E = [0, 1], [2, 3], [4, 5], [6, 7]; B = S + Ac
    hears = {**{k: I + S + Ac for k in I + Ac}, **{k: E + S + Ac for k in S + E}}   # split blanket: S hears E, Ac hears I
    A = -1.5 * np.eye(8)
    for i in range(8):
        for j in hears[i]:
            if i != j and rng.random() < 0.6: A[i, j] = rng.uniform(0.2, 0.8) * rng.choice([-1, 1])
    A[4, 0] = 0.5
    m = np.linalg.eigvals(A).real.max()
    if m > -0.1: A -= (m + 0.1) * np.eye(8)
    D = np.eye(8)
    print("Case 1: split blanket (sensory hears outside, active hears inside), independent noise")
    report(A, D, I, B, E)
    A2 = A.copy(); A2[4, 6] = 0.3
    print("\nCase 2: same system, but active state 4 also hears external state 6 (blanket state hears both sides)")
    report(A2, D, I, B, E)
    D3 = D.copy(); D3[2, 4] = D3[4, 2] = 0.3
    print("\nCase 3: split wiring, but sensory 2 and active 4 share noise")
    report(A, D3, I, B, E)
