# Blind Transit-Detection Pipeline with Error-Correcting-Code Contamination Screen

**Status:** exploratory research tool, validated end-to-end on 2 real TESS
targets. See `METHODS.md` for the full validation history, including bugs
found and fixed during testing. Not yet benchmarked against established
tools (BLS/TLS/`wotan`) -- see Limitations in `METHODS.md`.

## What this does

Given a TESS light curve (SPOC 2-minute cadence CSV + the original FITS
files), this pipeline:

1. Detrends per sector (rolling median).
2. Screens for contamination using **two merged candidate sources** --
   TESS's own SPOC quality flags (precise but category-limited) and a blind
   self-calibrated scan (broader coverage, lower precision) -- then decides
   whether to spare or excise each candidate by checking it against the
   pipeline's own bootstrap-validated signal, so real transits aren't
   mistaken for contamination.
3. Runs a blind box-least-squares-style period search.
4. Calibrates significance with a bootstrap null, not a fixed SNR cutoff.

## Requirements

- Python 3.9+
- `numpy`, `pandas`

No `astropy` dependency -- `fits_reader.py` is a minimal, self-contained
FITS binary-table parser (useful in sandboxed/offline environments; if you
have `astropy` available, feel free to swap it in for `fits_reader.read_lightcurve_fits`).

```bash
pip install numpy pandas
```

## Getting data

You need, per target:
- A cleaned CSV with columns `time_btjd, pdcsap_flux_normalized, quality, sector`
  (quality==0 rows only is fine; the pipeline reads the *original* FITS
  separately for the full quality bitmask).
- The original SPOC light-curve FITS file(s) for the same sectors, unmodified.

Both are downloadable from the [MAST Portal](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html):
search your target, filter Mission to TESS, Product Type to Timeseries.

## Usage

### Command line

```bash
cd src
python pipeline.py \
    --csv path/to/target_clean.csv \
    --fits path/to/sector1_lc.fits path/to/sector2_lc.fits \
    --period-min 0.5 --period-max 15.0 \
    --label "My Target"
```

Exits 0 if a significant signal was detected, 1 otherwise. Prints the full
screening and search log to stdout.

### As a module

```python
from pipeline import run_pipeline

result = run_pipeline(
    csv_path="target_clean.csv",
    fits_paths=["sector1_lc.fits", "sector2_lc.fits"],
    period_range=(0.5, 15.0),
)

print(result["period"], result["depth"], result["snr"], result["verdict"])
```

`result` is a dict: `period` (days), `depth` (fraction, multiply by 1e6 for
ppm), `duration_hours`, `snr`, `threshold` (the bootstrap significance
threshold that was cleared), `verdict` (`"DETECTED"` or `"NOT SIGNIFICANT"`),
`candidates` (every flagged region and its spare/excise decision), and
`n_excised` (cadences removed by the screen).

## Key parameters (all validated empirically -- see METHODS.md before changing)

| Parameter | Default | Meaning |
|---|---|---|
| `DETREND_WINDOW_DAYS` | 0.5 | Per-sector rolling-median window |
| `BLIND_BLOCK` | 91 | ~3h blocks for the blind contamination scan |
| `EXCISE_MARGIN_DAYS` | 3.0/24 | Window removed around each excised point -- validated to be necessary to cover a real event's recovery tail, not just its peak |
| `--n-periods` | 5000 | Period grid density |
| `--n-boot` | 20 | Bootstrap trials for the final significance threshold |

## Examples

`examples/` contains the exact scripts used to validate this pipeline
against WASP-18 b and Pi Mensae c, reproducing the numbers in `METHODS.md`.

## Validation summary

| Target | Recovered period | Published period | Depth | Verdict |
|---|---|---|---|---|
| WASP-18 b | 0.9415-0.9416 d | 0.94145 d | ~9,900-10,000 ppm | DETECTED |
| Pi Mensae c | 6.266-6.270 d | 6.2679 d | 147-315 ppm | DETECTED |

Both recovered blind (no ephemeris used before the search), on real public
TESS data, with the significance threshold self-calibrated from the data
rather than assumed.

## Known limitations

See `METHODS.md` Section 9 for the full, honest list. In short: n=2 real
targets, n=1 real contamination case used to validate the screen, no
comparison against established tools, coarse (box-fit) depth precision, and
the structural periodicity fix requires an already-validated signal to
reference -- unsolved for single-transit targets.

## License / citation

Research code, not a released package. If you use this, cite it as an
unpublished exploratory tool pending the validation work described in
`METHODS.md`.
