#!/usr/bin/env python3
"""
Florida Flow Regime - climatology and recent-sounding builder.

Pulls Florida radiosonde archives from the Iowa Environmental Mesonet (IEM),
computes the same parameters the dashboard computes live (1000-700 hPa
layer-mean wind + NWS TBW regime, PWAT, K index, ML lifted index, Thompson
index, MLCAPE, 500 hPa temp, 700-500 lapse rate, 850-500 mean RH), and writes:

  data/climo.json   day-of-year percentiles (+/-15 day window), regime
                    frequencies and within-regime percentiles per station
  data/recent.json  the last N days of observed soundings for the season log

Usage
  python scripts/build_climatology.py climo  [--start 1996 --end 2025]
  python scripts/build_climatology.py recent [--days 60]
  python scripts/build_climatology.py climo --csv-dir path/  # use local IEM-format CSVs
  python scripts/build_climatology.py import-table --table xmr.csv --station KXMR --start 1996 --end 2025
      # wide per-launch parameter table -> data/climo_KXMR.json, which the
      # climo run then uses in place of downloading that station

Pure standard library: runs on a stock GitHub Actions runner.
"""
import argparse, csv, io, json, math, os, sys, time, urllib.request, urllib.error
from collections import defaultdict
from datetime import datetime, timedelta, timezone

STATIONS = {  # IEM id -> label
    "KTBW": "Tampa Bay (Ruskin)",
    "KMFL": "Miami",
    "KJAX": "Jacksonville",
    "KXMR": "Cape Canaveral",
    "KKEY": "Key West",
    "KTLH": "Tallahassee",
}
IEM = "https://mesonet.agron.iastate.edu/cgi-bin/request/raob.py?station={st}&sts={s}&ets={e}"
VARS = ["pwat", "ki", "li", "ti", "mlcape", "t500", "lr75", "rhmid", "mwspd"]
PCTS = list(range(0, 101, 5))
DOY_STEP = 3
WINDOW = 15
LIGHT, STRONG = 5.0, 10.0

RD, CP, LV, EPS, G = 287.04, 1005.7, 2.501e6, 0.622, 9.80665
KAPPA = RD / CP


def isnum(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


# ---------------- thermodynamics (mirror of src/core.js) ----------------
def esat(tc): return 6.112 * math.exp(17.67 * tc / (tc + 243.5))
def td_from_e(e):
    l = math.log(e / 6.112); return 243.5 * l / (17.67 - l)
def mixr(p, tdc):
    e = esat(tdc); return EPS * e / (p - e)
def td_from_w(p, w): return td_from_e(w * p / (EPS + w))
def rh_of(tc, tdc): return 100 * esat(tdc) / esat(tc)
def theta(p, tk): return tk * (1000 / p) ** KAPPA
def tvirt(tk, w): return tk * (1 + w / EPS) / (1 + w)
def to_uv(d, s):
    r = math.radians(d); return -s * math.sin(r), -s * math.cos(r)
def from_uv(u, v): return (math.degrees(math.atan2(-u, -v)) + 360) % 360, math.hypot(u, v)


def normalize(levels):
    seen = {}
    for q in levels:
        p = q.get("p")
        if not isnum(p) or p < 50 or p > 1100: continue
        k = round(p * 10)
        cur = seen.setdefault(k, {"p": p})
        for f in ("t", "td", "dir", "spd"):
            if isnum(q.get(f)) and not isnum(cur.get(f)): cur[f] = q[f]
    out = sorted(seen.values(), key=lambda q: -q["p"])
    for q in out:
        if isnum(q.get("t")) and isnum(q.get("td")) and q["td"] > q["t"]: q["td"] = q["t"]
        if isnum(q.get("dir")) and isnum(q.get("spd")): q["u"], q["v"] = to_uv(q["dir"], q["spd"])
    return out


def interp(prof, p, key):
    lo = hi = None
    for q in prof:
        if not isnum(q.get(key)): continue
        if q["p"] >= p: lo = q
        else: hi = q; break
    if lo is not None and abs(lo["p"] - p) < 1e-6: return lo[key]
    if lo is None or hi is None: return None
    f = math.log(lo["p"] / p) / math.log(lo["p"] / hi["p"])
    return lo[key] + f * (hi[key] - lo[key])


def layer_mean_wind(prof, pbot=1000, ptop=700):
    if not prof: return None
    bot = min(pbot, prof[0]["p"])
    ub, vb, ut, vt = interp(prof, bot, "u"), interp(prof, bot, "v"), interp(prof, ptop, "u"), interp(prof, ptop, "v")
    if not all(isnum(x) for x in (ub, vb, ut, vt)): return None
    pts = [(bot, ub, vb)] + [(q["p"], q["u"], q["v"]) for q in prof if isnum(q.get("u")) and ptop < q["p"] < bot] + [(ptop, ut, vt)]
    su = sv = 0.0; gap = 0.0
    for a, b in zip(pts, pts[1:]):
        dp = a[0] - b[0]; gap = max(gap, dp)
        su += (a[1] + b[1]) / 2 * dp; sv += (a[2] + b[2]) / 2 * dp
    if gap > 160: return None
    u, v = su / (bot - ptop), sv / (bot - ptop)
    d, s = from_uv(u, v)
    return {"u": u, "v": v, "dir": d, "spd": s}


def ml_parcel(prof, depth=100):
    if not prof: return None
    psfc = prof[0]["p"]; ptop = psfc - depth
    sth = sw = 0.0; n = 0; p = psfc
    while p >= ptop - 1e-6:
        t, td = interp(prof, p, "t"), interp(prof, p, "td")
        if isnum(t) and isnum(td):
            sth += theta(p, t + 273.15); sw += mixr(p, td); n += 1
        p -= 5
    if n < 3: return None
    th, w = sth / n, sw / n
    tk = th * (psfc / 1000) ** KAPPA
    tdk = min(td_from_w(psfc, w) + 273.15, tk)
    return {"p": psfc, "tk": tk, "tdk": tdk, "w": w, "theta": th}


def lcl(tk, tdk, p):
    tl = 1 / (1 / (tdk - 56) + math.log(tk / tdk) / 800) + 56
    return tl, p * (tl / tk) ** (1 / KAPPA)


def moist_dtdp(p, tk):
    es = esat(tk - 273.15); rs = EPS * es / (p - es)
    return ((RD * tk + LV * rs) / p) / (CP + (LV * LV * rs * EPS) / (RD * tk * tk))


def moist_step(tk, p0, p1):
    n = max(1, math.ceil(abs(p1 - p0) / 5)); h = (p1 - p0) / n; p, t = p0, tk
    for _ in range(n):
        k1 = moist_dtdp(p, t); k2 = moist_dtdp(p + h / 2, t + h / 2 * k1)
        k3 = moist_dtdp(p + h / 2, t + h / 2 * k2); k4 = moist_dtdp(p + h, t + h * k3)
        t += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4); p += h
    return t


def parcel_path(pc, pressures):
    tl, pl = lcl(pc["tk"], pc["tdk"], pc["p"])
    tm, pm = tl, pl; out = []
    for p in pressures:
        if p >= pl:
            out.append((p, pc["theta"] * (p / 1000) ** KAPPA, pc["w"]))
        else:
            tm = moist_step(tm, pm, p); pm = p
            es = esat(tm - 273.15); out.append((p, tm, EPS * es / (p - es)))
    return out


def cape_of(prof, pc):
    temps = [q for q in prof if isnum(q.get("t"))]
    if not temps: return None
    ptop = max(temps[-1]["p"], 100)
    ps = []; p = pc["p"]
    while p >= ptop: ps.append(p); p -= 10
    if len(ps) < 5: return None
    path = parcel_path(pc, ps)
    B = []
    for (p, tk, w) in path:
        te, tde = interp(prof, p, "t"), interp(prof, p, "td")
        if not isnum(te): B.append(None); continue
        we = mixr(p, tde) if isnum(tde) else 0
        B.append(tvirt(tk, w) - tvirt(te + 273.15, we))
    _, pl = lcl(pc["tk"], pc["tdk"], pc["p"])
    lfc = el = -1
    for i, p in enumerate(ps):
        if isnum(B[i]) and B[i] > 0 and p <= pl + 1e-6:
            if lfc < 0: lfc = i
            el = i
    if lfc < 0: return 0.0
    cape = 0.0
    for i in range(lfc, el):
        if not (isnum(B[i]) and isnum(B[i + 1])): continue
        b = (B[i] + B[i + 1]) / 2
        if b > 0: cape += RD * b * math.log(ps[i] / ps[i + 1])
    return cape


def classify(d, s, light=LIGHT, strong=STRONG):
    if not (isnum(d) and isnum(s)): return None
    if s < light: return 1
    d = d % 360
    if d < 120: m, st = 2, 3
    elif d < 190: m, st = 6, 7
    elif d < 290: m, st = 4, 5
    else: m, st = 8, 9
    return st if s > strong else m


def metrics(levels):
    prof = normalize(levels); out = {}
    if len(prof) < 4: return out
    mw = layer_mean_wind(prof)
    if mw:
        out["mwdir"], out["mwspd"] = mw["dir"], mw["spd"]
        out["u"], out["v"] = mw["u"], mw["v"]
        out["regime"] = classify(mw["dir"], mw["spd"])
    T = lambda p: interp(prof, p, "t"); TD = lambda p: interp(prof, p, "td")
    t850, td850, t700, td700, t500, td500 = T(850), TD(850), T(700), TD(700), T(500), TD(500)
    if isnum(t500): out["t500"] = t500
    if all(isnum(x) for x in (t850, td850, t700, td700, t500)): out["ki"] = (t850 - t500) + td850 - (t700 - td700)
    if isnum(t700) and isnum(t500):
        dz = (RD / G) * ((t700 + t500) / 2 + 273.15) * math.log(700 / 500)
        out["lr75"] = (t700 - t500) / dz * 1000
    rhs = [rh_of(a, b) for a, b in ((t850, td850), (t700, td700), (t500, td500)) if isnum(a) and isnum(b)]
    if len(rhs) >= 2: out["rhmid"] = sum(rhs) / len(rhs)
    moist = [q for q in prof if isnum(q.get("td")) and q["p"] >= 300]
    if len(moist) >= 4 and moist[-1]["p"] <= 400:
        s = 0.0
        for a, b in zip(moist, moist[1:]):
            s += (mixr(a["p"], a["td"]) + mixr(b["p"], b["td"])) / 2 * (a["p"] - b["p"]) * 100
        out["pwat"] = s / G
    pc = ml_parcel(prof)
    if pc:
        if isnum(t500):
            out["li"] = t500 + 273.15 - parcel_path(pc, [500])[0][1]
            if isnum(out.get("ki")): out["ti"] = out["ki"] - out["li"]
        out["mlcape"] = cape_of(prof, pc)
    return out


# ---------------- IEM data ----------------
def fnum(x):
    try:
        v = float(x); return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def parse_iem_csv(text):
    """Return {valid_iso: [levels]} from IEM RAOB CSV text."""
    rdr = csv.DictReader(io.StringIO(text))
    cols = {c.lower(): c for c in (rdr.fieldnames or [])}
    def col(*names):
        for n in names:
            if n in cols: return cols[n]
        return None
    cv, cp = col("validutc", "valid"), col("pressure_mb", "pres", "pressure")
    ct, ctd = col("tmpc", "temperature"), col("dwpc", "dewpoint")
    cd, cs = col("drct", "direction"), col("speed_kts", "sknt", "speed")
    out = defaultdict(list)
    if not (cv and cp): return out
    for r in rdr:
        out[r[cv]].append({"p": fnum(r[cp]), "t": fnum(r.get(ct)) if ct else None, "td": fnum(r.get(ctd)) if ctd else None,
                           "dir": fnum(r.get(cd)) if cd else None, "spd": fnum(r.get(cs)) if cs else None})
    return out


def fetch(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "florida-flow-regime-dashboard"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read().decode("utf-8", "replace")
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"  retry {i + 1}: {e}", file=sys.stderr); time.sleep(5 * (i + 1))
    return ""


def parse_valid(s):
    s = s.strip().replace("Z", "").replace("T", " ")
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y%m%d%H%M"):
        try: return datetime.strptime(s[:19] if "-" in s else s, fmt).replace(tzinfo=timezone.utc)
        except ValueError: pass
    return None


def soundings_for(st, start, end, csv_dir=None):
    """Yield (datetime, cycle, metrics) for station between datetimes."""
    texts = []
    if csv_dir:
        for fn in sorted(os.listdir(csv_dir)):
            if fn.lower().endswith(".csv") and st.lower() in fn.lower():
                texts.append(open(os.path.join(csv_dir, fn), encoding="utf-8", errors="replace").read())
    else:
        cur = start
        while cur < end:  # one request per calendar year chunk
            nxt = min(datetime(cur.year + 1, 1, 1, tzinfo=timezone.utc), end)
            url = IEM.format(st=st, s=cur.strftime("%Y-%m-%dT%H:%MZ"), e=nxt.strftime("%Y-%m-%dT%H:%MZ"))
            print(f"  {st} {cur:%Y-%m-%d} -> {nxt:%Y-%m-%d}", file=sys.stderr)
            texts.append(fetch(url)); time.sleep(1.5); cur = nxt
    for text in texts:
        for valid, levels in parse_iem_csv(text).items():
            dt = parse_valid(valid)
            if not dt or not (start <= dt < end): continue
            hr = dt.hour + dt.minute / 60
            cycle = "12" if 10 <= hr <= 14 else ("00" if hr >= 22 or hr <= 2 else None)
            if not cycle: continue
            if cycle == "00" and hr >= 22: dt = dt + timedelta(hours=3)  # assign late launches to the next 00Z date
            m = metrics(levels)
            if m: yield dt, cycle, m


def percentiles(vals):
    xs = sorted(v for v in vals if isnum(v))
    if len(xs) < 20: return None
    out = []
    for p in PCTS:
        k = (len(xs) - 1) * p / 100; f = math.floor(k); c = min(f + 1, len(xs) - 1)
        out.append(round(xs[f] + (xs[c] - xs[f]) * (k - f), 1))
    return out


def doy(dt): return dt.timetuple().tm_yday


def station_climo(by, label, doys):
    """by[cycle][doy] -> list of metric dicts. Returns the station entry for climo.json."""
    n12 = sum(len(v) for v in by["12"].values())
    if n12 < 200:
        print(f"  skipped: only {n12} 12Z soundings", file=sys.stderr); return None
    sres = {"label": label, "n12": n12}
    for cyc in ("12", "00"):
        cres = {v: [] for v in VARS}
        freq = []
        for d in doys:
            pool = []
            for off in range(-WINDOW, WINDOW + 1):
                pool.extend(by[cyc].get((d - 1 + off) % 366 + 1, []))
            for v in VARS: cres[v].append(percentiles([m.get(v) for m in pool]))
            if cyc == "12":
                cnt = [0] * 9; tot = 0
                for m in pool:
                    r = m.get("regime")
                    if r: cnt[r - 1] += 1; tot += 1
                freq.append([round(c / tot, 3) for c in cnt] if tot >= 20 else None)
        sres[cyc] = cres
        if cyc == "12": sres["regimeFreq"] = freq
    br = {v: {} for v in VARS}
    wet = [m for d, ms in by["12"].items() if 152 <= d <= 273 for m in ms]
    for r in range(1, 10):
        sub = [m for m in wet if m.get("regime") == r]
        for v in VARS:
            q = percentiles([m.get(v) for m in sub])
            if q: br[v][str(r)] = q
        br.setdefault("_n", {})[str(r)] = len(sub)
    sres["byRegime"] = br
    return sres


def import_table(args):
    """Wide per-sounding parameter table (one row per launch) -> data/climo_<ST>.json.
    Uses the table's full-resolution PWAT, ML lifted index, Thompson, MLCAPE and
    1000-700 hPa mean vector; mid-level RH is rebuilt from the 850/700/500 columns
    so it matches the dashboard's definition."""
    def f(row, k):
        v = fnum(row.get(k)); return v
    best = {}  # (date, cycle) -> (distance from synoptic hour, metrics)
    nrows = 0
    with open(args.table, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            nrows += 1
            hr = fnum(row.get("nominal_hour_utc"))
            if hr is None or not row.get("date"): continue
            if 11 <= hr <= 13: cyc, dist = "12", abs(hr - 12)
            elif hr >= 23 or hr <= 1: cyc, dist = "00", min(abs(hr - 0), abs(hr - 24))
            else: continue
            d = datetime.strptime(row["date"][:10], "%Y-%m-%d").replace(tzinfo=timezone.utc)
            if hr >= 23: d += timedelta(days=1)  # late launch belongs to the next 00Z date
            if not (args.start <= d.year <= args.end): continue
            m = {"pwat": f(row, "pwat_mm"), "ki": f(row, "k_index"), "li": f(row, "lifted_index"),
                 "ti": f(row, "thompson_index"), "mlcape": f(row, "mlcape"), "t500": f(row, "t500"),
                 "lr75": f(row, "lr_700_500")}
            rh = [f(row, k) for k in ("rh850", "rh700", "rh500")]
            rh = [x for x in rh if isnum(x)]
            if len(rh) >= 2: m["rhmid"] = sum(rh) / len(rh)
            u, v = f(row, "u_1000_700"), f(row, "v_1000_700")
            if isnum(u) and isnum(v):
                wd, ws = from_uv(u, v)
                m.update({"u": u, "v": v, "mwdir": wd, "mwspd": ws, "regime": classify(wd, ws)})
            key = (d.strftime("%Y-%m-%d"), cyc)
            if key not in best or dist < best[key][0]: best[key] = (dist, d, m)
    by = {"12": defaultdict(list), "00": defaultdict(list)}
    for (ds, cyc), (_, d, m) in best.items(): by[cyc][doy(d)].append(m)
    doys = list(range(1, 367, DOY_STEP))
    sres = station_climo(by, args.label or STATIONS.get(args.station, args.station), doys)
    if not sres: sys.exit("Not enough 12Z soundings in the table.")
    yrs = sorted({k[0][:4] for k in best})
    sres["source"] = {"file": os.path.basename(args.table), "rows": nrows, "years": [int(yrs[0]), int(yrs[-1])],
                      "n00": sum(len(x) for x in by["00"].values())}
    os.makedirs("data", exist_ok=True)
    out = os.path.join("data", f"climo_{args.station}.json")
    with open(out, "w") as fh: json.dump(sres, fh, separators=(",", ":"))
    print(f"wrote {out}: {sres['n12']} 12Z and {sres['source']['n00']} 00Z soundings, {yrs[0]}-{yrs[-1]}", file=sys.stderr)
    if args.merge:
        path = os.path.join("data", "climo.json")
        if os.path.exists(path): base = json.load(open(path))
        else:
            base = {"meta": {"source": "Station override files", "years": [int(yrs[0]), int(yrs[-1])], "window_days": WINDOW,
                             "built": datetime.now(timezone.utc).isoformat(timespec="minutes"), "regime_thresholds_kt": [LIGHT, STRONG]},
                    "pcts": PCTS, "doys": doys, "doyStep": DOY_STEP, "stations": {}}
        base["stations"][args.station] = sres
        with open(path, "w") as fh: json.dump(base, fh, separators=(",", ":"))
        print(f"merged into {path}", file=sys.stderr)


def build_climo(args):
    start = datetime(args.start, 1, 1, tzinfo=timezone.utc)
    end = datetime(args.end + 1, 1, 1, tzinfo=timezone.utc)
    doys = list(range(1, 367, DOY_STEP))
    result = {"meta": {"source": "Iowa Environmental Mesonet RAOB archive", "years": [args.start, args.end],
                       "window_days": WINDOW, "built": datetime.now(timezone.utc).isoformat(timespec="minutes"),
                       "units": {"pwat": "mm", "ki": "C", "li": "C", "ti": "C", "mlcape": "J/kg", "t500": "C",
                                 "lr75": "C/km", "rhmid": "%", "mwspd": "kt"},
                       "regime_thresholds_kt": [LIGHT, STRONG]},
              "pcts": PCTS, "doys": doys, "doyStep": DOY_STEP, "stations": {}}
    for st, label in STATIONS.items():
        override = os.path.join("data", f"climo_{st}.json")
        if os.path.exists(override) and not args.csv_dir:
            result["stations"][st] = json.load(open(override))
            print(f"{st} {label}: using {override}", file=sys.stderr)
            continue
        print(f"{st} {label}", file=sys.stderr)
        by = {"12": defaultdict(list), "00": defaultdict(list)}
        for dt, cyc, m in soundings_for(st, start, end, args.csv_dir):
            by[cyc][doy(dt)].append(m)
        sres = station_climo(by, label, doys)
        if sres: result["stations"][st] = sres
    if not result["stations"]:
        sys.exit("No station had enough soundings. Check the IEM request URL or --csv-dir input.")
    os.makedirs("data", exist_ok=True)
    with open("data/climo.json", "w") as f: json.dump(result, f, separators=(",", ":"))
    print("wrote data/climo.json", file=sys.stderr)


def build_recent(args):
    end = datetime.now(timezone.utc) + timedelta(hours=1)
    start = (end - timedelta(days=args.days)).replace(hour=0, minute=0, second=0, microsecond=0)
    out = {"built": datetime.now(timezone.utc).isoformat(timespec="minutes"), "days": args.days, "stations": {}}
    for st in STATIONS:
        rows = []
        url = IEM.format(st=st, s=start.strftime("%Y-%m-%dT%H:%MZ"), e=end.strftime("%Y-%m-%dT%H:%MZ"))
        text = fetch(url); time.sleep(1.0)
        for valid, levels in parse_iem_csv(text).items():
            dt = parse_valid(valid)
            if not dt: continue
            hr = dt.hour
            cycle = "12" if 10 <= hr <= 14 else ("00" if hr >= 22 or hr <= 2 else None)
            if not cycle: continue
            if cycle == "00" and hr >= 22: dt += timedelta(hours=3)
            m = metrics(levels)
            if not m: continue
            row = {"date": dt.strftime("%Y-%m-%d"), "cycle": cycle}
            row.update({k: (round(v, 2) if isinstance(v, float) else v) for k, v in m.items()})
            rows.append(row)
        rows.sort(key=lambda r: (r["date"], r["cycle"]))
        out["stations"][st] = rows
        print(f"{st}: {len(rows)} soundings", file=sys.stderr)
    if not any(out["stations"].values()):
        sys.exit("No recent soundings returned from IEM.")
    os.makedirs("data", exist_ok=True)
    with open("data/recent.json", "w") as f: json.dump(out, f, separators=(",", ":"))
    print("wrote data/recent.json", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["climo", "recent", "import-table", "selftest"])
    ap.add_argument("--start", type=int, default=1996)
    ap.add_argument("--end", type=int, default=2025)
    ap.add_argument("--days", type=int, default=60)
    ap.add_argument("--csv-dir", default=None)
    ap.add_argument("--table", help="wide per-sounding parameter CSV for import-table")
    ap.add_argument("--station", default="KXMR")
    ap.add_argument("--label", default=None)
    ap.add_argument("--merge", action="store_true", help="also merge into data/climo.json")
    a = ap.parse_args()
    if a.mode == "climo": build_climo(a)
    elif a.mode == "recent": build_recent(a)
    elif a.mode == "import-table": import_table(a)
    else:
        levels = json.load(open(sys.argv[3] if len(sys.argv) > 3 else "/tmp/synth.json"))
        print(json.dumps({k: round(v, 3) if isinstance(v, float) else v for k, v in metrics(levels).items()}))
