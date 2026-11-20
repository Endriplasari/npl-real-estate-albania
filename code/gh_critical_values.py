"""
Gregory-Hansen (1996a) residual-based cointegration test with one endogenous
break, and Monte Carlo critical values for the sample sizes used in the paper.

Revision 2: the test statistic and its critical values are now produced by
the same function (`gh_stat`, ADF* with one augmentation lag, 15% trimming),
the regressor sets contain I(1) variables only (m = 1, 2), and the simulation
is validated against the asymptotic values of Gregory and Hansen (1996a,
Table 1) for the level-shift (C) and regime-shift (C/S) models.
"""
import pathlib, json, warnings
import numpy as np

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "output"

SEED = 19960303
NREP_FINITE = 5000
NREP_VALID = 2000

# Gregory and Hansen (1996a), Journal of Econometrics 70, Table 1, ADF* statistic
GH96A = {
    ("C", 1): {1: -5.13, 5: -4.61, 10: -4.34},
    ("CS", 1): {1: -5.47, 5: -4.95, 10: -4.68},
    ("C", 2): {1: -5.44, 5: -4.92, 10: -4.69},
    ("CS", 2): {1: -5.97, 5: -5.50, 10: -5.23},
}


def adf_nc(e, p=1):
    """ADF t-statistic on residuals, no deterministic terms, p lags."""
    de = np.diff(e)
    Z = [e[p:-1]]
    for i in range(1, p + 1):
        Z.append(de[p - i:-i])
    Z = np.column_stack(Z)
    y = de[p:]
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    r = y - Z @ b
    s2 = r @ r / (len(y) - Z.shape[1])
    XtXi = np.linalg.pinv(Z.T @ Z)
    return b[0] / np.sqrt(s2 * XtXi[0, 0])


def gh_stat(y, X, model, trim=0.15):
    """Minimum ADF* over candidate break dates. Returns (statistic, break index)."""
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    n = len(y)
    lo, hi = int(trim * n), int((1 - trim) * n)
    best, tb_best = np.inf, None
    ones = np.ones(n)
    for tb in range(lo, hi):
        D = (np.arange(n) >= tb).astype(float)
        if model == "C":
            Z = np.column_stack([ones, D, X])
        else:                                   # C/S: regime shift
            Z = np.column_stack([ones, D, X, X * D[:, None]])
        b, *_ = np.linalg.lstsq(Z, y, rcond=None)
        s = adf_nc(y - Z @ b)
        if s < best:
            best, tb_best = s, tb
    return float(best), tb_best


def mc_gh(n, m, model, nrep, rng):
    out = np.empty(nrep)
    for r in range(nrep):
        e = rng.standard_normal((n, m + 1))
        out[r] = gh_stat(np.cumsum(e[:, 0]), np.cumsum(e[:, 1:], axis=0), model)[0]
    return {q: float(np.quantile(out, q / 100)) for q in (1, 5, 10)}


if __name__ == "__main__":
    rng = np.random.default_rng(SEED)
    res = {"seed": SEED, "trim": 0.15, "lags": 1,
           "published_asymptotic_GH1996a_Table1": {f"{k[0]}_{k[1]}": v for k, v in GH96A.items()}}

    print("=" * 78)
    print("VALIDATION at n = 250 against Gregory and Hansen (1996a, Table 1)")
    print("=" * 78)
    val, maxdev = {}, 0.0
    print(f"{'model':>6}{'m':>3}{'level':>7}{'simulated':>11}{'published':>11}{'diff':>7}")
    for model in ("C", "CS"):
        for m in (1, 2):
            sim = mc_gh(250, m, model, NREP_VALID, rng)
            val[f"{model}_{m}"] = sim
            for q in (1, 5, 10):
                dev = sim[q] - GH96A[(model, m)][q]
                maxdev = max(maxdev, abs(dev))
                print(f"{model:>6}{m:3}{q:6}%{sim[q]:11.2f}{GH96A[(model, m)][q]:11.2f}{dev:7.2f}")
    res["validation_n250"] = val
    res["validation_max_abs_deviation"] = maxdev
    print(f"maximum absolute deviation: {maxdev:.2f}")

    print()
    print("=" * 78)
    print(f"FINITE-SAMPLE CRITICAL VALUES ({NREP_FINITE} replications, 15% trimming, one lag)")
    print("=" * 78)
    fin = {}
    for n, m in ((41, 2), (39, 1)):
        for model in ("C", "CS"):
            r = mc_gh(n, m, model, NREP_FINITE, rng)
            fin[f"{model}_m{m}_n{n}"] = r
            print(f"  n={n} m={m} {model:>2}: " + "  ".join(f"{q}%: {r[q]:6.2f}" for q in (1, 5, 10)))
    res["finite_sample"] = fin
    json.dump({str(k): v for k, v in res.items()}, open(str(OUT / "gh_cv.json"), "w"), indent=1, default=str)
    print("\nsaved -> output/gh_cv.json")
