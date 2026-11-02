# Florida Flow Regime

An automated summertime sea breeze flow-regime dashboard for Florida, built on the
NWS Tampa Bay regime scheme (12Z 1000–700 mb layer-mean vector wind, regimes 1–9),
with moisture and instability compared against sounding climatology.

## What's on the page

- **Regime wheel and call.** The nine TBW regimes drawn as the polar partition they
  are (direction sectors by speed bands). Today's call comes from the Tampa Bay 12Z
  sounding once it is in, otherwise from the model consensus.
- **Seven-day outlook.** GFS, ECMWF IFS, ECMWF AIFS, ICON, GEM and HRRR regimes in a
  Tampa Bay analysis box, plus consensus and agreement.
- **Ensemble regime probabilities.** Every GEFS (31) and ECMWF ENS (51) member is
  classified individually, giving regime odds per day for each system and combined
  (equal weight per system). The detail view adds flow-direction odds, ensemble-mean
  flow, ridge-axis spread and ridge-position regime odds. Members also plot on the
  regime wheel, and the ingredients panel can show ensemble means with member ranges.
- **Layer flow and ridge axis.** 1000–700 mb mean wind across the state (12Z, 18Z,
  00Z) with the subtropical ridge axis traced where easterlies turn westerly.
- **Coverage ingredients vs normal.** PWAT, Thompson, K, lifted index, MLCAPE,
  500 mb temperature, 700–500 mb lapse rate and mid-level RH for TBW, XMR, MFL, JAX,
  KEY and TLH, as percentiles for the date and within the current regime.
- **Regime rain climatology.** The TBW PoP graphics for the regime and time of day.
- **Florida soundings.** Latest launches with mini skew-Ts and computed parameters.
- **Season so far.** Each morning's Tampa Bay regime and PWAT against normal.

## Cape Canaveral page (`cape.html`)

Same data and code as the statewide page, focused on KSC and CCSFS:

- **AMU flow regimes** (Lambert 2007, NASA/CR-2007-214732): 1000–700 mb sectors at
  MFL, TBW and JAX set SW-1, SW-2, SE-1, SE-2, NE, NW or Other, and the XMR morning
  flow (above 4 kt) resolves Other/Missing and the SE-1/SW-2 boundary. The same rule
  runs on the soundings, every model and every GEFS and EC ENS member.
- **XMR flow wheel** with the AMU local sectors, and a station strip showing the
  four soundings' sectors behind the call.
- **Lightning odds from XMR flow and Thompson index**, from the AMU First Strike
  Timing Tool record (3,892 May–October days, 1989–2011): chance of lightning on
  KSC/CCSFS by 8-sector XMR 1000–700 mb flow and TI bin, looked up separately for the
  XMR sounding, the selected model and every GEFS and EC ENS member, plus the chance
  the first strike has happened by each local hour. Rebuild with
  `python scripts/build_cape_lightning.py "First Strike Timing Tool Excel2010.xlsm"`
  (writes `data/cape_ltg.json`; the workbook itself doesn't need to be in the repo).
- **Lightning odds** from the AMU flow-regime table (≥1 CG in the 10 KSC/CCSFS 5 n mi
  warning circles, 0700–0000 EDT, May–September 1989–2005), for the regime call and
  weighted by the ensemble regime odds.
- **Local map** with 5/10/20/30 n mi rings around the SLF and the Planning Forecast
  points.

The two pages link to each other from the header.

## Setup (GitHub Pages)

1. Create a repo (for example `FloridaFlowRegime`) and upload everything in this
   folder, keeping the `.github/workflows` and `scripts` paths.
2. Settings → Pages → deploy from the `main` branch, root folder.
3. Actions → **Soundings - climatology and season log** → Run workflow.
   The first run downloads 1996–2025 soundings for six Florida sites from the Iowa
   Environmental Mesonet and writes `data/climo.json` (expect roughly 30–60 minutes).
   After that it refreshes `data/recent.json` after each 12Z and 00Z cycle.

The page works before the Action finishes; percentiles and the season log appear
once the data files exist.

## Running locally

Open `index.html` directly. Live model and sounding data load from Open-Meteo and
IEM. Browsers block reading `data/*.json` from `file://`, so in Settings point the
climatology and season log URLs at your GitHub Pages copy, for example
`https://cyclonecizek.github.io/FloridaFlowRegime/data/climo.json`.

## Using your own sounding archive

`build_climatology.py` reads IEM-format CSVs (columns `validUTC, pressure_mb, tmpc,
dwpc, drct, speed_kts`) from a folder whose file names include the station ID:

    python scripts/build_climatology.py climo --start 1994 --end 2025 --csv-dir my_soundings/

## Cape Canaveral (XMR) climatology

XMR uses your own 1996–2025 archive rather than IEM. `data/climo_KXMR.json` was
built from the per-launch parameter table with:

    python scripts/build_climatology.py import-table --table xmr_allyear_soundings.csv --station KXMR --merge

Whenever the Action rebuilds the climatology, any `data/climo_<STATION>.json` file is
used in place of downloading that station. You don't need to commit the raw CSV.

## Method notes

- Regime thresholds default to the current TBW Quick Reference (regime 1 below 5 kt,
  strong above 10 kt). Older TBW material used 4 kt; change it in Settings.
- The model analysis box defaults to 27.0–28.5°N, 81.5–83.0°W as a stand-in for the
  TBW RAP analysis box; it is editable.
- The browser (`index.html`) and `scripts/build_climatology.py` share identical
  thermodynamics, so live values and climatology are computed the same way.
- Model profiles are compared against sounding climatology, so model biases carry
  through to the percentiles.
- The coverage signal is an experimental heuristic, not a calibrated probability.
- Ensemble-mean vectors shrink toward regime 1 when members disagree on direction,
  so each ensemble votes in the consensus with its most likely member regime.
- Ensemble requests are the heaviest on Open-Meteo's free limits, so processed results
  are cached in the browser for 3 hours. Refresh clears the cache.
- If a model is missing, hover its red chip in the status line to see Open-Meteo's
  error message.

## Sources

NWS Tampa Bay Thunderstorm Climatology and Quick Reference; Chavez et al. (2022),
*Weather and Forecasting* 37, 1867–1883; Lericos et al. (2002), *Weather and
Forecasting* 17, 83–98. Data: Open-Meteo, Iowa Environmental Mesonet.
