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
