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
- **Layer flow and ridge axis.** 1000–700 mb mean wind across the state (12Z, 18Z,
  00Z) with the subtropical ridge axis traced where easterlies turn westerly.
- **Coverage ingredients vs normal.** PWAT, Thompson, K, lifted index, MLCAPE,
  500 mb temperature, 700–500 mb lapse rate and mid-level RH for TBW, XMR, MFL, JAX,
  KEY and TLH, as percentiles for the date and within the current regime.
- **Regime rain climatology.** The TBW PoP graphics for the regime and time of day.
- **Florida soundings.** Latest launches with mini skew-Ts and computed parameters.
- **Season so far.** Each morning's Tampa Bay regime and PWAT against normal.

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

## Sources

NWS Tampa Bay Thunderstorm Climatology and Quick Reference; Chavez et al. (2022),
*Weather and Forecasting* 37, 1867–1883; Lericos et al. (2002), *Weather and
Forecasting* 17, 83–98. Data: Open-Meteo, Iowa Environmental Mesonet.
