"""
Calibrated housing-downturn scenario, revision 2.

Changes relative to the version deposited with the original submission:

  * Incremental losses only. Revision 1 applied the stressed impairment rate
    to total exposure, which counted the existing, already provisioned stock
    of impaired loans as a new loss. The stock and the flow are now separated,
    and existing provisions are netted against the stock.
  * Recovery and workout costs: legal and administrative costs as a share of
    proceeds, and discounting over the time to recovery.
  * The impairment rate depends on the equity position: exposures in negative
    equity after the price fall carry the multiple kappa (Haughwout and Okah,
    2009, report roughly twice the delinquency rate under negative equity).
  * Loan-to-value ratios are distributed within buckets instead of being
    represented by one value per bucket; the latest observed bucket share is
    reported alongside the survey average.
  * One-at-a-time sensitivity (replacing the two-way grid), a Monte Carlo loss
    distribution with tail measures, and a reverse stress test.
"""
import math
import numpy as np, pandas as pd

from common import load, jdump

SEED = 20261114
NDRAW = 200_000
d = pd.read_csv(str(__import__("common").DER / "dataset_quarterly_2016Q1_2026Q2.csv"), index_col=0)
S = {}

# ---------------- 1. Decomposition of the NPL ratio: d ln(ratio) = d ln(stock) - d ln(credit)
def decomp(a, b):
    r0, r1 = d.loc[a, "npl"], d.loc[b, "npl"]
    s0, s1 = d.loc[a, "stok_npl"], d.loc[b, "stok_npl"]
    c0, c1 = d.loc[a, "kredi_tot"], d.loc[b, "kredi_tot"]
    dr = math.log(r1 / r0); ds = math.log(s1 / s0); dc = math.log(c1 / c0)
    return dict(period=f"{a}-{b}", npl_from=r0, npl_to=r1, dln_ratio=dr, dln_stock=ds, dln_credit=-dc,
                share_numerator=ds / dr, share_denominator=-dc / dr,
                stock_from_bn=s0 / 1000, stock_to_bn=s1 / 1000, credit_from_bn=c0 / 1000, credit_to_bn=c1 / 1000)
S["decomp"] = [decomp("2016Q1", "2026Q2"), decomp("2016Q1", "2021Q4"), decomp("2021Q4", "2026Q2")]

# ---------------- 2. Inputs (2026Q2, ALL million, Bank of Albania)
last = d.loc["2026Q2"]
TOT = float(last["kredi_tot"]); MORT = float(last["hipoteka"]); BRE = float(last["pp_biznes"])
NPL_STOCK = float(last["stok_npl"]); R0 = float(last["npl"]) / 100.0
RWA = 1194029.43           # FSI 1.1.2, June 2026
CAP = 235630.63            # FSI 1.1.1, June 2026 (regulatory capital)
CAR0 = 100 * CAP / RWA
NET_NPL_TO_CAP = 6.14      # FSI 1.4.2, June 2026: NPLs net of provisions / regulatory capital, %
COV = 1 - (NET_NPL_TO_CAP / 100 * CAP) / NPL_STOCK          # provision coverage of the NPL stock
ROA = 1.48                 # FSI 2.5, June 2026: pre-tax return on assets, annualised, %
CAP_TO_ASSETS = 9.93       # FSI 1.3.2, June 2026: regulatory capital / total assets, %
PROFIT = ROA / 100 * CAP / (CAP_TO_ASSETS / 100)             # implied annual pre-tax profit
MINCAR = 12.0

LTV_LOW_SHARES = {"2021H2": 0.50, "2022H1": 0.86, "2022H2": 0.80, "2023H1": 0.85, "2023H2": 0.60,
                  "2024H1": 0.65, "2025H1": 0.62, "2025H2": 0.38}
SL_MEAN = float(np.mean(list(LTV_LOW_SHARES.values())))
SL_LATEST = LTV_LOW_SHARES["2025H2"]

BASE = dict(sL=SL_MEAN,                 # share of the residential book with LTV <= 60%
            low=(0.40, 0.60),           # LTV range, low bucket (uniform)
            high=(0.60, 0.90),          # LTV range, high bucket (uniform)
            corp=(0.50, 0.90),          # LTV range, corporate real estate (uniform)
            h_res=0.20, h_com=0.30,     # forced-sale haircuts
            w=0.10,                     # workout cost, share of proceeds (Doing Business 2020: 10% of estate)
            T=2.0,                      # years to recovery (Doing Business 2020: 2.0 years)
            i=0.05,                     # discount rate
            kappa=2.0,                  # impairment multiple under negative equity relative to positive equity
            cov=COV)                    # provisions already held against the impaired stock
S["inputs"] = dict(total_loans=TOT, mortgage=MORT, mortgage_share=100 * MORT / TOT, business_re=BRE,
                   business_re_share=100 * BRE / TOT, npl_ratio=100 * R0, npl_stock=NPL_STOCK, rwa=RWA,
                   capital=CAP, car0=CAR0, provision_coverage=COV, implied_pretax_profit=PROFIT,
                   loss_to_min_car_bn=(CAP - MINCAR / 100 * RWA) / 1000, sL_mean=SL_MEAN,
                   sL_latest=SL_LATEST, ltv_low_shares=LTV_LOW_SHARES, base=BASE)

NG = 240     # grid points per bucket


def bucket(E, lo, hi, s, m, h, P):
    """Losses on a book of size E with LTV uniform on [lo, hi]."""
    ell = lo + (np.arange(NG) + 0.5) * (hi - lo) / NG
    delta = (1 + P["i"]) ** (-P["T"])
    lgd = np.maximum(0.0, 1 - (1 - s) * (1 - h) * (1 - P["w"]) * delta / ell)
    neg = ell > (1 - s)                                  # negative equity after the price fall
    r = np.minimum(1.0, np.where(neg, m * P["kappa"], m) * R0)
    new_flow = np.maximum(0.0, r - R0)                   # newly impaired share of the book
    loss_new = new_flow * lgd
    loss_old = R0 * np.maximum(0.0, lgd - P["cov"])      # shortfall on the existing impaired stock
    return dict(E=E, neg=float(neg.mean()), lgd=float(lgd.mean()),
                lgd_neg=float(lgd[neg].mean()) if neg.any() else 0.0,
                flow=float(E * new_flow.mean()), loss_new=float(E * loss_new.mean()),
                loss_old=float(E * loss_old.mean()))


def run(s, m, **over):
    P = dict(BASE); P.update(over)
    b = [bucket(MORT * P["sL"], *P["low"], s, m, P["h_res"], P),
         bucket(MORT * (1 - P["sL"]), *P["high"], s, m, P["h_res"], P),
         bucket(BRE, *P["corp"], s, m, P["h_com"], P)]
    loss_new = sum(x["loss_new"] for x in b); loss_old = sum(x["loss_old"] for x in b)
    loss = loss_new + loss_old
    flow = sum(x["flow"] for x in b)
    res_E = b[0]["E"] + b[1]["E"]
    return dict(shock=s, mult=m, mult_neg=m * P["kappa"],
                rate_pos=100 * min(1, m * R0), rate_neg=100 * min(1, m * P["kappa"] * R0),
                neg_res=100 * (b[0]["E"] * b[0]["neg"] + b[1]["E"] * b[1]["neg"]) / res_E,
                neg_corp=100 * b[2]["neg"],
                lgd_res=100 * (b[0]["E"] * b[0]["lgd"] + b[1]["E"] * b[1]["lgd"]) / res_E,
                lgd_corp=100 * b[2]["lgd"],
                new_impaired_bn=flow / 1000, loss_new_bn=loss_new / 1000, loss_old_bn=loss_old / 1000,
                loss_bn=loss / 1000, loss_pct_capital=100 * loss / CAP, car_after=100 * (CAP - loss) / RWA,
                npl_after=100 * (NPL_STOCK + flow) / TOT,
                loss_corp_share=100 * (b[2]["loss_new"] + b[2]["loss_old"]) / loss if loss > 0 else 0.0,
                loss_pct_profit=100 * loss / PROFIT)


def run_v1(s, m, sL=SL_MEAN):
    """Revision 1 calculation, kept for the reconciliation reported to the reviewers."""
    lgd = lambda ltv, h: max(0.0, 1 - (1 - s) * (1 - h) / ltv)
    rs = min(1.0, R0 * m)
    loss = MORT * rs * (sL * lgd(0.50, 0.20) + (1 - sL) * lgd(0.75, 0.20)) + BRE * rs * lgd(0.70, 0.30)
    return dict(loss_bn=loss / 1000, car_after=100 * (CAP - loss) / RWA)


# ---------------- 3. Scenario grid
grid = [run(s, m) for s in (0.0, 0.20, 0.30) for m in (2, 3, 4)]
S["grid"] = grid
S["grid_v1"] = [dict(shock=s, mult=m, **run_v1(s, m)) for s in (0.0, 0.20, 0.30) for m in (2, 3, 4)]
S["grid_latest_ltv"] = [run(s, m, sL=SL_LATEST) for s in (0.20, 0.30) for m in (2, 3, 4)]
print(f"provision coverage {COV:.3f}; implied pre-tax profit {PROFIT/1000:.1f} bn; CAR0 {CAR0:.2f}")
print(f"{'shock':>6}{'m':>3}{'neg res':>8}{'neg corp':>9}{'LGD res':>8}{'LGD corp':>9}{'new NPL':>8}{'L new':>7}"
      f"{'L old':>7}{'loss':>7}{'%cap':>6}{'CAR':>7}{'NPL':>6}{'corp%':>7}{'v1 loss':>9}")
for g, v in zip(grid, S["grid_v1"]):
    print(f"{g['shock']:6.0%}{g['mult']:3}{g['neg_res']:8.1f}{g['neg_corp']:9.1f}{g['lgd_res']:8.1f}{g['lgd_corp']:9.1f}"
          f"{g['new_impaired_bn']:8.1f}{g['loss_new_bn']:7.1f}{g['loss_old_bn']:7.2f}{g['loss_bn']:7.1f}"
          f"{g['loss_pct_capital']:6.1f}{g['car_after']:7.2f}{g['npl_after']:6.1f}{g['loss_corp_share']:7.1f}{v['loss_bn']:9.1f}")

# ---------------- 4. One-at-a-time sensitivity around (30% fall, multiple 3)
ref = run(0.30, 3)
tests = [
    ("Price fall", "20%", "40%", lambda: run(0.20, 3), lambda: run(0.40, 3)),
    ("Impairment multiple, positive equity", "2", "4", lambda: run(0.30, 2), lambda: run(0.30, 4)),
    ("Negative-equity multiple (kappa)", "1", "3", lambda: run(0.30, 3, kappa=1.0), lambda: run(0.30, 3, kappa=3.0)),
    ("Forced-sale haircut, residential / commercial", "10% / 20%", "30% / 40%",
     lambda: run(0.30, 3, h_res=0.10, h_com=0.20), lambda: run(0.30, 3, h_res=0.30, h_com=0.40)),
    ("Corporate LTV range", "40-80%", "60-100%",
     lambda: run(0.30, 3, corp=(0.40, 0.80)), lambda: run(0.30, 3, corp=(0.60, 1.00))),
    ("Corporate haircut alone", "20%", "50%", lambda: run(0.30, 3, h_com=0.20), lambda: run(0.30, 3, h_com=0.50)),
    ("Share of residential book with LTV <= 60%", "86% (highest survey)", "38% (latest survey)",
     lambda: run(0.30, 3, sL=0.86), lambda: run(0.30, 3, sL=SL_LATEST)),
    ("Upper limit of the high-LTV bucket", "80%", "100%",
     lambda: run(0.30, 3, high=(0.60, 0.80)), lambda: run(0.30, 3, high=(0.60, 1.00))),
    ("Workout cost, share of proceeds", "5%", "15%", lambda: run(0.30, 3, w=0.05), lambda: run(0.30, 3, w=0.15)),
    ("Time to recovery", "1 year", "4 years", lambda: run(0.30, 3, T=1.0), lambda: run(0.30, 3, T=4.0)),
    ("Provision coverage of existing NPLs", "70%", "50%", lambda: run(0.30, 3, cov=0.70), lambda: run(0.30, 3, cov=0.50)),
]
sens = []
for name, lo_lab, hi_lab, flo, fhi in tests:
    a, b = flo(), fhi()
    sens.append(dict(parameter=name, low_label=lo_lab, high_label=hi_lab, loss_low=a["loss_bn"], loss_high=b["loss_bn"],
                     car_low=a["car_after"], car_high=b["car_after"], range_bn=abs(b["loss_bn"] - a["loss_bn"])))
sens.sort(key=lambda r: -r["range_bn"])
S["sensitivity"] = dict(reference=ref, rows=sens, min_car=min(min(r["car_low"], r["car_high"]) for r in sens))
print(f"\nreference (30%, x3): loss {ref['loss_bn']:.1f} bn, CAR {ref['car_after']:.2f}")
for r in sens:
    print(f"  {r['parameter']:46} {r['low_label']:>22} -> {r['loss_low']:6.1f} bn | {r['high_label']:>20} -> "
          f"{r['loss_high']:6.1f} bn | range {r['range_bn']:5.1f} | CAR {r['car_low']:.2f} / {r['car_high']:.2f}")

# ---------------- 5. Monte Carlo loss distribution
rng = np.random.default_rng(SEED)
tri = lambda a, c, b, n: rng.triangular(a, c, b, n)

def simulate(shock_draw, n=NDRAW):
    P = dict(m=tri(1.5, 3.0, 5.0, n), kappa=tri(1.0, 2.0, 3.0, n), h_res=tri(0.10, 0.20, 0.35, n),
             h_com=tri(0.20, 0.30, 0.50, n), sL=tri(0.38, SL_MEAN, 0.86, n), hi_top=tri(0.80, 0.90, 1.00, n),
             corp_mid=tri(0.60, 0.70, 0.80, n), w=tri(0.05, 0.10, 0.15, n), T=tri(1.0, 2.0, 4.0, n),
             i=tri(0.03, 0.05, 0.08, n), cov=tri(0.50, COV, 0.70, n), s=shock_draw(n))
    loss = np.empty(n)
    u = (np.arange(40) + 0.5) / 40                       # 40-point grid within each bucket
    delta = (1 + P["i"]) ** (-P["T"])
    def book(E, lo, hi, h):
        ell = lo[:, None] + u[None, :] * (hi - lo)[:, None]
        lgd = np.maximum(0.0, 1 - ((1 - P["s"]) * (1 - h) * (1 - P["w"]) * delta)[:, None] / ell)
        neg = ell > (1 - P["s"])[:, None]
        r = np.minimum(1.0, np.where(neg, (P["m"] * P["kappa"])[:, None], P["m"][:, None]) * R0)
        return E * (np.maximum(0.0, r - R0) * lgd + R0 * np.maximum(0.0, lgd - P["cov"][:, None])).mean(axis=1)
    ones = np.ones(n)
    loss = (book(MORT * P["sL"], 0.40 * ones, 0.60 * ones, P["h_res"])
            + book(MORT * (1 - P["sL"]), 0.60 * ones, P["hi_top"], P["h_res"])
            + book(BRE, P["corp_mid"] - 0.20, P["corp_mid"] + 0.20, P["h_com"]))
    return loss

def summarise(loss):
    q = lambda p: float(np.quantile(loss, p))
    es = lambda p: float(loss[loss >= np.quantile(loss, p)].mean())
    car = lambda L: 100 * (CAP - L) / RWA
    return dict(mean_bn=float(loss.mean()) / 1000, median_bn=q(0.5) / 1000, p90_bn=q(0.90) / 1000,
                p95_bn=q(0.95) / 1000, p99_bn=q(0.99) / 1000, es95_bn=es(0.95) / 1000, es99_bn=es(0.99) / 1000,
                max_bn=float(loss.max()) / 1000,
                mean_pct_cap=100 * float(loss.mean()) / CAP, p95_pct_cap=100 * q(0.95) / CAP,
                p99_pct_cap=100 * q(0.99) / CAP, es99_pct_cap=100 * es(0.99) / CAP,
                car_mean=car(float(loss.mean())), car_p95=car(q(0.95)), car_p99=car(q(0.99)),
                car_es99=car(es(0.99)), car_min=car(float(loss.max())),
                prob_car_below_15=float((car(loss) < 15).mean()), prob_car_below_12=float((car(loss) < MINCAR).mean()),
                prob_loss_above_profit=float((loss > PROFIT).mean()))

MC = {"n": NDRAW, "seed": SEED,
      "distributions": "triangular (min, mode, max): m (1.5, 3, 5); kappa (1, 2, 3); residential haircut (0.10, 0.20, 0.35); "
                       "commercial haircut (0.20, 0.30, 0.50); low-LTV share (0.38, survey mean, 0.86); upper limit of high "
                       "bucket (0.80, 0.90, 1.00); corporate LTV midpoint (0.60, 0.70, 0.80) with half-width 0.20; workout "
                       "cost (0.05, 0.10, 0.15); years to recovery (1, 2, 4); discount rate (0.03, 0.05, 0.08); provision "
                       "coverage (0.50, FSI-implied, 0.70)"}
MC["fall_20"] = summarise(simulate(lambda n: np.full(n, 0.20)))
MC["fall_30"] = summarise(simulate(lambda n: np.full(n, 0.30)))
MC["fall_uncertain"] = summarise(simulate(lambda n: tri(0.10, 0.30, 0.50, n)))
S["monte_carlo"] = MC
print()
for k in ("fall_20", "fall_30", "fall_uncertain"):
    v = MC[k]
    print(f"MC {k:15} mean {v['mean_bn']:5.1f}  median {v['median_bn']:5.1f}  P95 {v['p95_bn']:5.1f}  P99 {v['p99_bn']:5.1f}  "
          f"ES99 {v['es99_bn']:5.1f} bn | CAR mean {v['car_mean']:.2f} P95 {v['car_p95']:.2f} P99 {v['car_p99']:.2f} "
          f"min {v['car_min']:.2f} | P(CAR<15) {v['prob_car_below_15']:.4f} P(L>profit) {v['prob_loss_above_profit']:.3f}")

# ---------------- 6. Reverse stress test: what breaches the 12% minimum?
def solve_m(s, target=MINCAR, **over):
    lo, hi = 1.0, 1 / R0
    if run(s, hi, **over)["car_after"] > target:
        return None
    for _ in range(60):
        mid = (lo + hi) / 2
        if run(s, mid, **over)["car_after"] > target:
            lo = mid
        else:
            hi = mid
    return hi

rev = []
for s in (0.20, 0.30, 0.40, 0.50):
    for lab, over in (("baseline assumptions", {}), ("latest LTV distribution, haircuts 30% / 40%",
                                                      dict(sL=SL_LATEST, h_res=0.30, h_com=0.40))):
        m = solve_m(s, **over)
        if m is None:
            rev.append(dict(shock=s, variant=lab, mult=None))
            continue
        g = run(s, m, **over)
        rev.append(dict(shock=s, variant=lab, mult=m, rate_pos=g["rate_pos"], rate_neg=g["rate_neg"],
                        npl_after=g["npl_after"], loss_bn=g["loss_bn"], new_impaired_bn=g["new_impaired_bn"]))
S["reverse"] = dict(target_car=MINCAR, loss_needed_bn=(CAP - MINCAR / 100 * RWA) / 1000, rows=rev,
                    peak_npl_2014=25.0, sample_peak_npl=float(d["npl"].max()))
print(f"\nreverse stress test: loss needed to reach {MINCAR:.0f}% = {S['reverse']['loss_needed_bn']:.1f} bn")
for r in rev:
    if r["mult"] is None:
        print(f"  fall {r['shock']:.0%} {r['variant']}: not reached"); continue
    print(f"  fall {r['shock']:.0%} {r['variant']:44}: multiple {r['mult']:5.1f} -> impairment rate {r['rate_pos']:.1f}% "
          f"(positive equity) / {r['rate_neg']:.1f}% (negative equity), system NPL ratio {r['npl_after']:.1f}%")

# ---------------- 7. Benchmark against the IMF staff stress test (Koosakul and Topi, 2025)
# Their credit-risk test raises NPLs on the whole loan book by 100%, 200% and 440% with a 65% provisioning rate.
bench = []
for inc in (1.0, 2.0, 4.4):
    L = inc * NPL_STOCK * 0.65
    bench.append(dict(npl_increase_pct=100 * inc, loss_bn=L / 1000, car_after=100 * (CAP - L) / RWA,
                      car_change=100 * (CAP - L) / RWA - CAR0))
S["benchmark_imf_assumptions"] = bench
print("\nIMF-style benchmark (all loans, 65% provisioning):",
      [(b["npl_increase_pct"], round(b["loss_bn"], 1), round(b["car_after"], 2)) for b in bench])

jdump(S, "scenario.json")
print("\nsaved -> output/scenario.json")
