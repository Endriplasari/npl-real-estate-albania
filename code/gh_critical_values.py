import pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW  = ROOT / "data" / "raw"
DER  = ROOT / "data" / "derived"
OUT  = ROOT / "output"
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json

rng = np.random.default_rng(19960303)


def adf_nc(e, p=1):
    """ADF t-statistic on residuals, no deterministics, p lags."""
    de = np.diff(e)
    n = len(de) - p
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
    n = len(y)
    lo, hi = int(trim * n), int((1 - trim) * n)
    best = np.inf
    ones = np.ones(n)
    for tb in range(lo, hi):
        D = (np.arange(n) >= tb).astype(float)
        if model == "C":
            Z = np.column_stack([ones, D, X])
        elif model == "CS":
            Z = np.column_stack([ones, D, X, X * D[:, None]])
        else:                                    # C/S/T : regime and trend shift
            t = np.arange(n, dtype=float)
            Z = np.column_stack([ones, D, t, X, X * D[:, None]])
        b, *_ = np.linalg.lstsq(Z, y, rcond=None)
        e = y - Z @ b
        s = adf_nc(e)
        if s < best:
            best = s
    return best


def mc_gh(n, m, model, nrep=2000):
    out = np.empty(nrep)
    for r in range(nrep):
        e = rng.standard_normal((n, m + 1))
        y = np.cumsum(e[:, 0])
        X = np.cumsum(e[:, 1:], axis=0)
        out[r] = gh_stat(y, X, model)
    return {q: float(np.quantile(out, q / 100)) for q in (1, 2.5, 5, 10)}


print("=" * 86)
print("VALIDATION: simulated large-sample values for the regime-and-trend-shift model")
print("against the published asymptotic values of Gregory and Hansen (1996b, Table 1)")
print("=" * 86)
PUB = {1: {1: -6.02, 2.5: -5.72, 5: -5.50, 10: -5.24},
       2: {1: -6.45, 2.5: -6.17, 5: -5.96, 10: -5.72},
       3: {1: -6.89, 2.5: -6.65, 5: -6.32, 10: -6.16},
       4: {1: -7.31, 2.5: -7.06, 5: -6.84, 10: -6.58}}
print(f"{'m':>3}{'level':>8}{'simulated (n=250)':>20}{'published':>12}{'diff':>8}")
val = {}
for m in (3, 4):
    sim = mc_gh(250, m, "CST", nrep=600)
    val[m] = sim
    for q in (1, 5, 10):
        print(f"{m:3}{q:7}%{sim[q]:20.2f}{PUB[m][q]:12.2f}{sim[q]-PUB[m][q]:8.2f}")

print()
print("=" * 86)
print("TABLE. Finite-sample critical values for the Gregory-Hansen ADF* statistic")
print("Monte Carlo, 2,000 replications, n = 42, 15% trimming, one lag")
print("=" * 86)
cv = {}
for model, name in [("C", "level shift (C)"), ("CS", "regime shift (C/S)")]:
    for m in (3, 4):
        r = mc_gh(42, m, model, nrep=2000)
        cv[f"{model}_{m}"] = r
        print(f"  {name:22} m={m}   "
              + "   ".join(f"{q}%: {r[q]:7.3f}" for q in (1, 2.5, 5, 10)))

json.dump({"validation_CST_n250": {str(k): v for k, v in val.items()},
           "published_CST_asymptotic": {str(k): v for k, v in PUB.items()},
           "finite_sample_n42": cv},
          open(str(OUT / "gh_cv.json"), "w"), indent=1)
print("\nsaved -> gh_cv.json")
