import pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW  = ROOT / "data" / "raw"
DER  = ROOT / "data" / "derived"
OUT  = ROOT / "output"
import json, math, numpy as np, pandas as pd
d = pd.read_csv(str(DER / "dataset_quarterly_2016Q1_2026Q2.csv"), index_col=0)
S = {}

# ---------------- 1. Decomposition of NPL-ratio decline: d ln(ratio) = d ln(stock) - d ln(credit)
def decomp(a, b):
    r0, r1 = d.loc[a, "npl"], d.loc[b, "npl"]
    s0, s1 = d.loc[a, "stok_npl"], d.loc[b, "stok_npl"]
    c0, c1 = d.loc[a, "kredi_tot"], d.loc[b, "kredi_tot"]
    dr = math.log(r1 / r0); ds = math.log(s1 / s0); dc = math.log(c1 / c0)
    return dict(period=f"{a}-{b}", npl_from=r0, npl_to=r1, dln_ratio=dr, dln_stock=ds, dln_credit=-dc,
                share_numerator=ds / dr, share_denominator=-dc / dr,
                stock_from_bn=s0/1000, stock_to_bn=s1/1000, credit_from_bn=c0/1000, credit_to_bn=c1/1000)
S["decomp"] = [decomp("2016Q1", "2026Q2"), decomp("2016Q1", "2021Q4"), decomp("2021Q4", "2026Q2")]
for x in S["decomp"]:
    print(f"{x['period']}: NPL {x['npl_from']:.2f}->{x['npl_to']:.2f} | numerator {x['share_numerator']*100:.0f}% | denominator {x['share_denominator']*100:.0f}%")

# ---------------- 2. Inputs (2026Q2, million ALL, BSH)
last = d.loc["2026Q2"]
tot = last["kredi_tot"]; mort = last["hipoteka"]; bre = last["pp_biznes"]
npl_stock = last["stok_npl"]; npl = last["npl"]
RWA = 1194029.0        # BSH FSI 1.1.2, June 2026
CAP = 235631.0         # BSH FSI 1.1.1, June 2026
CAR0 = 100 * CAP / RWA
S["inputs"] = dict(total_loans=tot, mortgage=mort, mortgage_share=100*mort/tot, business_re=bre,
                   business_re_share=100*bre/tot, npl_ratio=npl, npl_stock=npl_stock, rwa=RWA, capital=CAP, car0=CAR0)

# LTV structure from BoA real-estate surveys: share of credit-financed purchases with loan <= 60% of value
ltv_low_shares = {"2021H2":0.50, "2022H1":0.86, "2022H2":0.80, "2023H1":0.85, "2023H2":0.60,
                  "2024H1":0.65, "2025H1":0.62, "2025H2":0.38}
sL = float(np.mean(list(ltv_low_shares.values())))     # baseline share in low-LTV bucket
S["ltv_low_shares"] = ltv_low_shares; S["sL_baseline"] = sL

def lgd(ltv, shock, haircut):
    cov = (1 - shock) * (1 - haircut) / ltv
    return max(0.0, 1 - cov)

def run(shock, pd_mult, sL=sL, h_res=0.20, h_com=0.30, ltv_low=0.50, ltv_high=0.75, ltv_bus=0.70):
    pd_base = npl / 100.0
    pd_s = min(1.0, pd_base * pd_mult)
    # residential book
    lgd_res = sL * lgd(ltv_low, shock, h_res) + (1 - sL) * lgd(ltv_high, shock, h_res)
    loss_res = mort * pd_s * lgd_res
    # business real-estate-purpose book
    lgd_bus = lgd(ltv_bus, shock, h_com)
    loss_bus = bre * pd_s * lgd_bus
    loss = loss_res + loss_bus
    new_npl_stock = npl_stock + (mort + bre) * (pd_s - pd_base)
    return dict(shock=shock, pd_mult=pd_mult, pd_stress=100*pd_s, lgd_res=100*lgd_res, lgd_bus=100*lgd_bus,
                loss_bn=loss/1000, loss_pct_capital=100*loss/CAP, car_after=100*(CAP - loss)/RWA,
                npl_after=100*new_npl_stock/tot,
                share_neg_equity_res=100*(1-sL) if lgd(ltv_high, shock, 0) > 0 else 0.0)

grid = []
for shock in (0.0, 0.20, 0.30):
    for pm in (2, 3, 4):
        grid.append(run(shock, pm))
S["grid"] = grid
for g in grid:
    print(f"shock {g['shock']:.0%} PDx{g['pd_mult']}: PD={g['pd_stress']:.1f}% LGDres={g['lgd_res']:.1f}% LGDbus={g['lgd_bus']:.1f}% loss={g['loss_bn']:.1f}bn ({g['loss_pct_capital']:.1f}% cap) CAR {CAR0:.2f}->{g['car_after']:.2f} NPL->{g['npl_after']:.1f}%")

# sensitivity: haircut and low-LTV share, at shock 30%, PDx3
sens = []
for h in (0.10, 0.20, 0.30):
    for s in (0.50, sL, 0.80):
        r = run(0.30, 3, sL=s, h_res=h, h_com=h+0.10)
        sens.append(dict(haircut=h, sL=s, loss_bn=r["loss_bn"], car_after=r["car_after"]))
S["sensitivity"] = sens
# capital buffer: how large a loss (bn) brings CAR to 12% (regulatory min) 
S["loss_to_12pct_bn"] = (CAP - 0.12 * RWA) / 1000
print("Loss absorbable before CAR hits 12%:", round(S["loss_to_12pct_bn"],1), "bn ALL")
json.dump(S, open(str(OUT / "scenario.json"), "w"), indent=1, default=float)
