"""F5b.0 (plan §8.12, step 2): leak forecast WITHOUT Park's law, evaluated open loop.

The replica's leak is Park's law (julia/SMDFRHRM.jl, `leak_flow`):
    L(t) = min(-d + A d (t - 500)^(1/B), 0)
The twin may not use it (it generates the data). It sees noisy measurements y = L + v,
v ~ N(0, sigma^2), every h = 0.5 s, and forecasts L(tau + j h), j = 1..K = 26 (the N2
horizon), at every query tau.

Forecasters
  persist   L_j = y_tau
  trend2    L_j = y_tau + j (y_tau - y_{tau-h})
  huber     robust line (Huber, IRLS) on the last W = 10 s, extrapolated
  kalman    local-linear-trend Kalman filter (level, slope); q_l, q_s calibrated offline
            by the mean multi-step log score on training realisations, r = sigma^2
  kalman_r  the same, with an intervention: if the normalised innovation exceeds 5, the
            jump is absorbed in the level (x_level = y, P_ll += nu^2), the slope is kept
The Kalman filters give their own Sigma_L (multi-step covariance); the other forecasters get
an empirical sigma_j per horizon from the training set.

Tests (out of sample): the nominal and extreme families (draw ranges of
julia/f5a2_monte_carlo.jl), measurement noise sigma in {0, 0.1, 0.3} kg/s, and break jumps
(1-2 per run, Delta ~ -U(1, 5) kg/s at random times). Also the 14 closed-loop realisations of
F5a.3 (results/f5a2_monte_carlo/runs_extreme.csv).

Outputs: results/f5b0_leak_forecast/{summary.csv, kalman_params.json}.
Usage: .venv/Scripts/python.exe python/f5b0_leak_forecast.py
"""
import csv
import json
import os
import itertools
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = os.path.join(ROOT, "results", "f5b0_leak_forecast")
os.makedirs(OUT, exist_ok=True)
H, K, TF = 0.5, 26, 500.0
BURN = 60.0          # the filter starts 60 s before the first query
SPAN = 240.0         # queries over 240 s (481 queries, as in F5a.3)
W = 10.0             # Huber window [s]
SIGMAS = (0.0, 0.1, 0.3)


def leak(t, A, B, d, jumps=()):
    L = -d + A * d * np.clip(t - TF, 0, None) ** (1 / B)
    for tj, dj in jumps:
        L = L + dj * (t >= tj)
    return np.minimum(L, 0.0)


def draw(rng, family, n, jumps=False):
    """Realisations: (A, B, d, t_first_query, jumps)."""
    out = []
    for _ in range(n):
        if family == "extreme":
            A, B, d = 0.1 + 0.1 * rng.random(), 10 + 5 * rng.random(), 15 + 10 * rng.random()
            t0 = 660 + 90 * rng.random()
        else:
            A, B, d = 0.2 + 0.2 * rng.random(), 9 + 6 * rng.random(), 8 + 2 * rng.random()
            t0 = 1140 + 100 * rng.random()
        t0 = np.round(t0 / H) * H
        J = []
        if jumps:
            for _ in range(rng.integers(1, 3)):
                J.append((t0 + SPAN * rng.random(), -(1 + 4 * rng.random())))
        out.append((A, B, d, t0, tuple(J)))
    return out


def series(reals, sigma, rng):
    """Measurement grid (from t0 - BURN to t0 + SPAN), measurements and truth incl. horizon."""
    n = int(round((BURN + SPAN) / H)) + 1
    m = n + K
    T = np.zeros((len(reals), m)); Ltrue = np.zeros_like(T)
    for i, (A, B, d, t0, J) in enumerate(reals):
        T[i] = t0 - BURN + H * np.arange(m)
        Ltrue[i] = leak(T[i], A, B, d, J)
    Y = Ltrue[:, :n] + sigma * rng.standard_normal((len(reals), n))
    q0 = int(round(BURN / H))           # index of the first query
    return Y, Ltrue, q0, n


def truth_windows(Ltrue, q0, n):
    """E[r, k, j] = L(tau_k + j h), k over queries, j = 1..K."""
    idx = np.arange(q0, n)[:, None] + np.arange(1, K + 1)[None, :]
    return Ltrue[:, idx]


# ------------------------------------------------------------------ forecasters
def fc_persist(Y, q0, n):
    return np.repeat(Y[:, q0:n, None], K, axis=2)


def fc_trend2(Y, q0, n):
    s = Y[:, q0:n] - Y[:, q0 - 1:n - 1]
    return Y[:, q0:n, None] + s[:, :, None] * np.arange(1, K + 1)[None, None, :]


def fc_huber(Y, q0, n, sigma):
    w = int(round(W / H)); x = np.arange(-w + 1, 1, dtype=float)
    F = np.zeros((Y.shape[0], n - q0, K))
    for k, i in enumerate(range(q0, n)):
        Yw = Y[:, i - w + 1:i + 1]
        X = np.column_stack([np.ones(w), x])
        beta = np.linalg.lstsq(X, Yw.T, rcond=None)[0]               # (2, R)
        for _ in range(5):                                           # IRLS
            res = Yw - (X @ beta).T
            s = 1.4826 * np.median(np.abs(res - np.median(res, axis=1, keepdims=True)), axis=1)
            s = np.maximum(s, max(sigma, 1e-6))
            c = 1.345 * s[:, None]
            wt = np.where(np.abs(res) <= c, 1.0, c / np.maximum(np.abs(res), 1e-12))
            for r in range(Y.shape[0]):
                Xw = X * wt[r][:, None]
                beta[:, r] = np.linalg.solve(X.T @ Xw + 1e-12 * np.eye(2), Xw.T @ Yw[r])
        F[:, k, :] = beta[0][:, None] + beta[1][:, None] * np.arange(1, K + 1)[None, :]
    return F


def kalman(Y, q0, n, ql, qs, r, reset=False):
    """Local linear trend: x = (level, slope per step). Returns forecasts (R, nq, K), their
    variances (R, nq, K) and, for the last query, nothing else. P is common to all
    realisations unless `reset` (then per realisation)."""
    R_ = Y.shape[0]
    Fm = np.array([[1.0, 1.0], [0.0, 1.0]]); Q = np.diag([ql, qs])
    x = np.zeros((R_, 2)); x[:, 0] = Y[:, 0]
    P = np.tile(np.diag([r + 1.0, 1.0]), (R_, 1, 1))
    jj = np.arange(1, K + 1)
    proc = np.cumsum(ql + (jj - 1) ** 2 * qs)                         # sum_{m<j} (q_l + m^2 q_s)
    Fc = np.zeros((R_, n - q0, K)); Vc = np.zeros_like(Fc)
    loglik = 0.0
    for i in range(1, n):
        x = x @ Fm.T
        P = Fm @ P @ Fm.T + Q
        S = P[:, 0, 0] + r
        nu = Y[:, i] - x[:, 0]
        if i > 20:
            loglik += np.sum(-0.5 * (np.log(2 * np.pi * S) + nu ** 2 / S))
        Kg = P[:, :, 0] / S[:, None]
        x = x + Kg * nu[:, None]
        P = P - Kg[:, :, None] * P[:, 0, :][:, None, :]
        if reset:                       # intervention: absorb a jump of the break in the level
            big = np.abs(nu) / np.sqrt(S) > 5
            if big.any():
                x[big, 0] = Y[big, i]
                P[big, 0, 0] += nu[big] ** 2
        if i >= q0:
            k = i - q0
            Fc[:, k, :] = x[:, [0]] + x[:, [1]] * jj[None, :]
            Vc[:, k, :] = (P[:, 0, 0][:, None] + 2 * jj[None, :] * P[:, 0, 1][:, None]
                           + jj[None, :] ** 2 * P[:, 1, 1][:, None] + proc[None, :])
    return Fc, Vc, loglik


def logscore(E, V):
    return np.mean(-0.5 * (np.log(2 * np.pi * V) + E ** 2 / V))


def calibrate(Y, Lt, q0, n, sigma):
    """(q_l, q_s) maximising the mean multi-step log score on the training set."""
    r = max(sigma ** 2, 1e-4)          # floor: sensor resolution 0.01 kg/s
    T = truth_windows(Lt, q0, n)
    best = None
    for lql, lqs in itertools.product(np.arange(-12, -1.5, 1.0), np.arange(-16, -5.5, 1.0)):
        Fc, Vc, _ = kalman(Y, q0, n, 10 ** lql, 10 ** lqs, r)
        s = logscore(T - Fc, Vc)
        if best is None or s > best[0]:
            best = (s, lql, lqs)
    # refine on a half-decade grid around the best
    s0, a, b = best
    for lql, lqs in itertools.product(a + np.array([-0.5, 0, 0.5]), b + np.array([-0.5, 0, 0.5])):
        Fc, Vc, _ = kalman(Y, q0, n, 10 ** lql, 10 ** lqs, r)
        s = logscore(T - Fc, Vc)
        if s > best[0]:
            best = (s, lql, lqs)
    return dict(q_l=10 ** best[1], q_s=10 ** best[2], r=r, logscore=best[0])


def jump_masks(reals, q0, n):
    """(R, nq) masks: `cross` = a jump falls inside the forecast window (tau, tau + K h]
    (unforeseeable); `post` = a jump fell in (tau - 10 s, tau] (recovery)."""
    nq = n - q0
    cross = np.zeros((len(reals), nq), bool); post = np.zeros_like(cross)
    for i, (A, B, d, t0, J) in enumerate(reals):
        tau = t0 + H * np.arange(nq)
        for tj, _ in J:
            cross[i] |= (tj > tau) & (tj <= tau + K * H)
            post[i] |= (tj <= tau) & (tj > tau - 10.0)
    return cross, post


def metrics(name, E, V, mask=None):
    """E, V: (R, nq, K) errors and predictive variances; `mask` (R, nq) selects queries."""
    if mask is not None:
        E = E[mask][None]; V = np.broadcast_to(V, mask.shape + (K,))[mask][None]
    a = np.abs(E)
    z = E / np.sqrt(V)
    from scipy.stats import norm, kstest
    pit = norm.cdf(z[:, :, [0, 12, 25]]).ravel()
    return dict(forecaster=name,
                maxerr_p99=np.quantile(a.max(axis=2), 0.99),
                maxerr_max=a.max(),
                rmse_j1=np.sqrt(np.mean(E[:, :, 0] ** 2)), rmse_j13=np.sqrt(np.mean(E[:, :, 12] ** 2)),
                rmse_j26=np.sqrt(np.mean(E[:, :, 25] ** 2)),
                cov90_j26=np.mean(np.abs(z[:, :, 25]) <= 1.645), cov95_all=np.mean(np.abs(z) <= 1.96),
                pit_ks=kstest(pit, "uniform").statistic,
                logscore=logscore(E, V))


def main():
    rng = np.random.default_rng(20261116)
    rows, params = [], {}
    for fam in ("nominal", "extreme"):
        train = draw(rng, fam, 60)
        for sigma in SIGMAS:
            Ytr, Ltr, q0, n = series(train, sigma, rng)
            par = calibrate(Ytr, Ltr, q0, n, sigma)
            params[f"{fam}_sigma{sigma}"] = par
            Ttr = truth_windows(Ltr, q0, n)
            # empirical sigma_j of the simple forecasters (training residuals)
            emp = {}
            for nm, F in (("persist", fc_persist(Ytr, q0, n)), ("trend2", fc_trend2(Ytr, q0, n)),
                          ("huber", fc_huber(Ytr, q0, n, sigma))):
                emp[nm] = np.maximum(np.mean((Ttr - F) ** 2, axis=(0, 1)), 1e-12)
            for jumps in (False, True):
                test = draw(rng, fam, 60, jumps)
                Y, Lt, q0, n = series(test, sigma, rng)
                T = truth_windows(Lt, q0, n)
                cand = {"persist": fc_persist(Y, q0, n), "trend2": fc_trend2(Y, q0, n),
                        "huber": fc_huber(Y, q0, n, sigma)}
                cross, post = jump_masks(test, q0, n)
                sets = [("all", None)] if not jumps else [("no_jump_in_window", ~cross), ("recovery_10s", post & ~cross)]
                res = {nm: (T - F, np.broadcast_to(emp[nm], F.shape)) for nm, F in cand.items()}
                for nm, rs in (("kalman", False), ("kalman_r", True)):
                    Fc, Vc, _ = kalman(Y, q0, n, par["q_l"], par["q_s"], par["r"], reset=rs)
                    res[nm] = (T - Fc, Vc)
                for sname, msk in sets:
                    for nm, (E, V) in res.items():
                        rows.append(dict(family=fam, sigma=sigma, set=sname, **metrics(nm, E, V, msk)))
                print(f"{fam:8s} sigma={sigma:.1f} jumps={jumps!s:5s} done")
    with open(os.path.join(OUT, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    with open(os.path.join(OUT, "kalman_params.json"), "w") as f:
        json.dump(params, f, indent=2)
    hdr = f"{'family':8s} {'sig':>4s} {'set':18s} {'forecaster':9s} {'p99max':>8s} {'max':>8s} {'rmse1':>8s} {'rmse26':>8s} {'cov90@26':>8s} {'cov95':>6s} {'KS':>5s} {'logS':>7s}"
    print(hdr)
    for r in rows:
        print(f"{r['family']:8s} {r['sigma']:4.1f} {r['set']:18s} {r['forecaster']:9s} {r['maxerr_p99']:8.4f} {r['maxerr_max']:8.3f} "
              f"{r['rmse_j1']:8.4f} {r['rmse_j26']:8.4f} {r['cov90_j26']:8.3f} {r['cov95_all']:6.3f} {r['pit_ks']:5.2f} {r['logscore']:7.2f}")
    print(json.dumps(params, indent=1))


if __name__ == "__main__":
    main()
