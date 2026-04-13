# README

## Usage

For double crystal ball:

Fit: `python src/fit/fit_dcb.py data.root --output results.json`
Plot: `python src/fit/plot_dcb_driver.py results.json`

For gaussian:

Fit: `python src/fit/fit_gauss.py data.root --output results.json`
Plot: `python src/fit/plot_gaussian_driver.py results.json`

You can adjust the configurations in the corresponding `plot_*_driver.py` file, which includes:

- `description`
- `output`
- `title`
- `xlim` (optional)
- `pull_ylim` (optional)
- `build_test`
