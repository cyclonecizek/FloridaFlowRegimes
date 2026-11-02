#!/usr/bin/env python3
"""
Cape Canaveral lightning lookup from the AMU First Strike Timing Tool.

Reads the hidden RawData sheet of "First Strike Timing Tool Excel2010.xlsm"
(one row per warm-season day: XMR morning 1000-700 mb mean wind direction and
speed, Thompson index, lightning occurrence on KSC/CCAFS and the local hour of
the first strike) and writes data/cape_ltg.json:

  * probability of lightning by month group x 8-sector XMR flow x Thompson bin,
    shrunk toward the parent cell so sparse cells do not swing wildly
    (beta-binomial / empirical Bayes with K pseudo-days)
  * first-strike hour counts (local time) by month group and sector

Usage:
  python scripts/build_cape_lightning.py "First Strike Timing Tool Excel2010.xlsm"
Requires openpyxl (pip install openpyxl).
"""
import json, math, os, sys
from collections import defaultdict
from datetime import datetime, timezone

K = 10  # pseudo-days of prior weight toward the parent cell
SECTORS = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
TI_BINS = ["<25", "25-34", ">=34"]
GROUPS = {"May": [5], "JJA": [6, 7, 8], "Sep": [9], "Oct": [10], "All": [5, 6, 7, 8, 9, 10]}
HOURS = list(range(7, 24))  # first strike local hour 07-23


def sector8(d):
    """45 WS FRN sectors used by the tool (N >335 or <=20, then 45-46 deg wedges)."""
    d = d % 360
    if d > 335 or d <= 20: return "N"
    for name, lo, hi in (("NE", 20, 66), ("E", 66, 110), ("SE", 110, 155), ("S", 155, 200),
                         ("SW", 200, 245), ("W", 245, 290), ("NW", 290, 335)):
        if lo < d <= hi: return name


def ti_bin(t):
    if t is None: return None
    return "<25" if t < 25 else ("25-34" if t < 34 else ">=34")


def shrink(yes, n, parent):
    return (yes + K * parent) / (n + K) if n + K > 0 else parent


def main(path):
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    rows = list(wb["RawData"].iter_rows(values_only=True))
    hdr = rows[0]
    col = {h: i for i, h in reversed(list(enumerate(hdr))) if h}  # first occurrence wins
    days = []
    for r in rows[1:]:
        d = r[col["Date1"]]
        if not isinstance(d, datetime): continue
        wdir, spd, ti = r[col["WDir"]], r[col["Speed"]], r[col["TI"]]
        if not isinstance(wdir, (int, float)): continue
        ltg = r[col["LtgOccur"]] == "Yes"
        fh = r[col["FirstStrikeLocalHour"]]
        days.append({"m": d.month, "sec": sector8(float(wdir)), "spd": float(spd) if isinstance(spd, (int, float)) else None,
                     "ti": float(ti) if isinstance(ti, (int, float)) else None, "ltg": ltg,
                     "fh": int(fh) if ltg and isinstance(fh, (int, float)) else None, "y": d.year})
    years = sorted({x["y"] for x in days})
    out = {"meta": {"source": "AMU First Strike Timing Tool (Excel 2010), RawData sheet", "file": os.path.basename(path),
                    "days": len(days), "years": [years[0], years[-1]], "months": "May-October",
                    "predictand": "Lightning on KSC/CCAFS during the day; first strike hour is local time (07-23)",
                    "predictors": "XMR morning 1000-700 mb mean wind (8 sectors) and Thompson index",
                    "shrinkage_pseudo_days": K, "built": datetime.now(timezone.utc).isoformat(timespec="minutes")},
           "sectors": SECTORS, "tiBins": TI_BINS, "hours": HOURS, "groups": {}}
    for g, months in GROUPS.items():
        sub = [x for x in days if x["m"] in months]
        n, yes = len(sub), sum(x["ltg"] for x in sub)
        base = yes / n
        G = {"months": months, "n": n, "yes": yes, "p": round(base, 4), "sector": {}, "timing": {}}
        tim_all = [0] * len(HOURS)
        for x in sub:
            if x["fh"] in HOURS: tim_all[HOURS.index(x["fh"])] += 1
        G["timing"]["all"] = tim_all
        for s in SECTORS:
            ss = [x for x in sub if x["sec"] == s]
            sn, sy = len(ss), sum(x["ltg"] for x in ss)
            ps = shrink(sy, sn, base)
            cell = {"n": sn, "yes": sy, "p": round(ps, 4), "ti": {}}
            for b in TI_BINS:
                tb = [x for x in ss if ti_bin(x["ti"]) == b]
                tn, ty = len(tb), sum(x["ltg"] for x in tb)
                cell["ti"][b] = {"n": tn, "yes": ty, "p": round(shrink(ty, tn, ps), 4)}
            tim = [0] * len(HOURS)
            for x in ss:
                if x["fh"] in HOURS: tim[HOURS.index(x["fh"])] += 1
            G["timing"][s] = tim
            G["sector"][s] = cell
        out["groups"][g] = G
    os.makedirs("data", exist_ok=True)
    with open("data/cape_ltg.json", "w") as f: json.dump(out, f, separators=(",", ":"))
    print(f"wrote data/cape_ltg.json from {len(days)} days, {years[0]}-{years[-1]}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) < 2: sys.exit(__doc__)
    main(sys.argv[1])
