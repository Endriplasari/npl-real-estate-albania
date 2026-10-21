"""
Size and power of the ARDL bounds F-test under the M1 design.

Two things the placebo regressions and the standard-error-based minimum
detectable effect do NOT establish, and that are therefore computed here:

  (1) the empirical rejection rate of the bounds F-test under the null of
      no level relationship, evaluated against the finite-sample critical
      values actually used in the paper;
  (2) the rejection rate of the same test under alternatives in which a
      level relationship exists with a known speed of adjustment.

The conditional ECM is coded directly so the simulation is independent of
the statsmodels routine that produced the reported statistics; agreement
between the two is itself the check.

Design mirrors M1: 41 quarterly observations, four regressors, ARDL(2,1)
in unrestricted error-correction form, unrestricted intercept (Case III),
one fixed dummy. 39 effective observations, 12 parameters, 27 residual df.
"""
import pathlib, json
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "output"

N_OBS = 41          # sample length, as in M1
K = 4               # regressors entering in levels
P_LAG, Q_LAG = 2, 1 # ARDL orders selected by AIC for M1
NREP_SIZE = 5000
NREP_POWER = 2000

# finite-sample critical values reported for M1 (statsmodels, 20,000 sims, Case III)
CV = {"lower": {90.0: 2.667634, 95.0: 3.216148, 99.0: 4.457236},
      "upper": {90.0: 3.849789, 95.0: 4.588766, 99.0: 6.214031}}

rng = np.random.default_rng(20261020)


def bounds_F(y, X, dummy, p=P_LAG, q=Q_LAG):
    """Bounds F for H0: phi = 0 and beta = 0 in the conditional ECM."""
    T, k = X.shape
    dy = np.diff(y)
    dX = np.diff(X, axis=0)
    start = max(p, q)                      # first usable index into the differenced series
    n = T - 1 - start + 1
    cols = [np.ones(n)]                    # unrestricted intercept
    lev = []
    # levels at t-1
    y_lag = y[start - 1:T - 1]
    cols.append(y_lag); lev.append(len(cols) - 1)
    for j in range(k):
        cols.append(X[start - 1:T - 1, j]); lev.append(len(cols) - 1)
    # short-run: dy lags 1..p-1
    for i in range(1, p):
        cols.append(dy[start - 1 - i:T - 1 - i])
    # short-run: dX lags 0..q-1
    for j in range(k):
        for i in range(0, q):
            cols.append(dX[start - 1 - i:T - 1 - i, j])
    cols.append(dummy[start:T])            # fixed dummy
    Z = np.column_stack(cols)
    yv = dy[start - 1:]
    if Z.shape[0] != yv.shape[0]:
        m = min(Z.shape[0], yv.shape[0]); Z = Z[:m]; yv = yv[:m]
    b, *_ = np.linalg.lstsq(Z, yv, rcond=None)
    e = yv - Z @ b
    rss_u = e @ e
    dfr = len(yv) - Z.shape[1]
    keep = [c for c in range(Z.shape[1]) if c not in lev]
    Zr = Z[:, keep]
    br, *_ = np.linalg.lstsq(Zr, yv, rcond=None)
    er = yv - Zr @ br
    rss_r = er @ er
    return ((rss_r - rss_u) / len(lev)) / (rss_u / dfr), dfr


def make_dummy(T):
    """COVID moratorium dummy, as in the data: four consecutive quarters."""
    d = np.zeros(T); d[17:21] = 1.0
    return d


def draw_regressors(T, integ):
    """integ: 'I1' all unit root, 'I0' all stationary, 'mixed' 2+2 as in M1."""
    e = rng.standard_normal((T, K))
    if integ == "I1":
        return np.cumsum(e, axis=0)
    if integ == "I0":
        return e
    X = np.empty((T, K))
    X[:, :2] = np.cumsum(e[:, :2], axis=0)   # exposure shares: I(1)
    X[:, 2:] = e[:, 2:]                      # GDP growth, real policy rate: I(0)
    return X


# ----------------------------------------------------------------- SIZE
def size(integ, nrep=NREP_SIZE):
    D = make_dummy(N_OBS)
    F = np.empty(nrep)
    for r in range(nrep):
        y = np.cumsum(rng.standard_normal(N_OBS))     # no level relationship
        X = draw_regressors(N_OBS, integ)
        F[r], dfr = bounds_F(y, X, D)
    out = {"n_rep": nrep, "mean_F": float(F.mean())}
    for lvl in (90.0, 95.0, 99.0):
        a = (100 - lvl) / 100
        out[f"reject_upper_{int(lvl)}"] = float((F > CV["upper"][lvl]).mean())
        out[f"reject_lower_{int(lvl)}"] = float((F > CV["lower"][lvl]).mean())
        out[f"nominal_{int(lvl)}"] = a
    out["q95_simulated"] = float(np.quantile(F, 0.95))
    return out


# ----------------------------------------------------------------- POWER
def power(phi, sigma, nrep=NREP_POWER, integ="mixed"):
    """DGP: y error-corrects towards y = beta'x with speed |phi|.

    beta loads on the first regressor only, so the long-run relationship is
    y = x1 + equilibrium error with standard deviation sigma.
    """
    D = make_dummy(N_OBS)
    beta = np.zeros(K); beta[0] = 1.0
    crit = CV["upper"][95.0]
    rej = 0
    for r in range(nrep):
        X = draw_regressors(N_OBS, integ)
        eq = X @ beta
        y = np.empty(N_OBS); y[0] = eq[0]
        e = rng.standard_normal(N_OBS) * sigma
        dX = np.diff(X, axis=0, prepend=X[:1])
        for t in range(1, N_OBS):
            y[t] = y[t - 1] + phi * (y[t - 1] - eq[t - 1]) + 0.5 * (dX[t] @ beta) + e[t]
        F, _ = bounds_F(y, X, D)
        rej += F > crit
    return rej / nrep


if __name__ == "__main__":
    res = {"design": {"n": N_OBS, "k": K, "p": P_LAG, "q": Q_LAG,
                      "case": "III (unrestricted intercept, no trend)",
                      "critical_values": CV,
                      "cv_source": "statsmodels UECM.bounds_test(asymptotic=False, nsim=20000), as reported for M1"}}

    print("=" * 78)
    print("EMPIRICAL SIZE OF THE BOUNDS F-TEST")
    print(f"{NREP_SIZE} replications, no level relationship imposed")
    print("=" * 78)
    print(f"{'regressors':<26}{'mean F':>9}{'vs upper 5%':>13}{'vs lower 5%':>13}{'vs upper 10%':>14}")
    res["size"] = {}
    for integ, lab in [("I1", "all I(1)"), ("I0", "all I(0)"), ("mixed", "two I(1), two I(0)")]:
        s = size(integ)
        res["size"][integ] = s
        print(f"{lab:<26}{s['mean_F']:9.3f}{s['reject_upper_95']:13.4f}"
              f"{s['reject_lower_95']:13.4f}{s['reject_upper_90']:14.4f}")
    print("\nUnder the null the rejection rate at the I(1) upper bound should be at or")
    print("below the nominal level; at the I(0) lower bound it should be close to it.")

    print()
    print("=" * 78)
    print("POWER OF THE BOUNDS F-TEST AGAINST ERROR-CORRECTING ALTERNATIVES")
    print(f"{NREP_POWER} replications, 5% level, I(1) upper bound")
    print("=" * 78)
    sigmas = [0.25, 0.50, 1.00]
    print(f"{'|phi|':>7}{'half-life (q)':>15}" + "".join(f"{'sigma='+str(s):>12}" for s in sigmas))
    res["power"] = {}
    for phi in (-0.10, -0.20, -0.30, -0.40, -0.50, -0.60, -0.75):
        hl = np.log(0.5) / np.log(1 + phi)
        row = {}
        line = f"{abs(phi):7.2f}{hl:15.1f}"
        for s in sigmas:
            pw = power(phi, s)
            row[str(s)] = pw
            line += f"{pw:12.3f}"
        res["power"][f"{abs(phi):.2f}"] = {"half_life_q": float(hl), **row}
        print(line)
    print("\nsigma is the standard deviation of the equilibrium error relative to unit")
    print("innovations in the regressors: 0.25 is a tight long-run relationship.")

    json.dump(res, open(str(OUT / "size_power.json"), "w"), indent=1)
    print("\nsaved -> output/size_power.json")
