# README

## Usage

One-shot fit + plot:

For double crystal ball:
`python src/fit/dcb_driver.py data.root [--outdir out]`

For gaussian:
`python src/fit/gauss_driver.py data.root [--outdir out]`

Each driver runs the fit, saves `results_dcb.json` / `results_gauss.json`,
and writes `dcb_mass_fit.png` / `gauss_mass_fit.png` into `--outdir` (default: `out/`).

## Step-by-step usage

For double crystal ball:

Fit: `python src/fit/fit_dcb.py data.root --outfile results.json`
Plot: `python src/fit/plot_dcb.py results.json`

For gaussian:

Fit: `python src/fit/fit_gauss.py data.root --output results.json`
Plot: `python src/fit/plot_gauss.py results.json`

You can adjust the configurations in the corresponding `plot_*.py` file's `CONFIG` dict, which includes:

- `description`
- `output`
- `title`
- `xlim` (optional)
- `pull_ylim` (optional)
- `build_text`

## Dalitz plot fit (bin-by-bin)

Fit the eta mass peak separately in each (m12, m23) Dalitz bin to get a
background-subtracted signal yield per bin, instead of raw candidate counts.
Requires mean/sigma/c0/c1 from a prior inclusive fit (same `--xmin`/`--xmax`)
so each bin only has to float its yields:

```
python src/fit/fit_gauss.py data.root --output results_gauss.json
python src/fit/fit_dalitz.py data.root \
    --mean 547.9 --sigma 8.2 --c0 -0.05 --c1 0.02 \
    --output results_dalitz.json
python src/fit/plot_dalitz_fit.py results_dalitz.json
```

(Substitute the actual fitted `mean`/`sigma`/`c0`/`c1` from `results_gauss.json`.)
`fit_dalitz.py` outputs a JSON with fitted `n_sig`/`n_bkg` (+ errors) per
Dalitz bin; `plot_dalitz_fit.py` renders it as a yield-map PNG (bins that
weren't fit or didn't converge are masked white, and the physical signal
boundary is drawn as reference lines).

Default grid (`--m12-min/max`, `--m23-min/max`, 25x25 bins) is padded around
the physical eta -> mu mu gamma signal window rather than spanning the full
kinematic triangle -- see `fit_dalitz.py`'s docstring and the `M12_PHYS`/
`M23_PHYS` derivation in `plot_dalitz_fit.py` for where those numbers come
from. Override them if you're fitting a different decay or want a coarser/
finer grid.
