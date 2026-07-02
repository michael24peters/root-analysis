"""
gauss_driver.py

Run the Gaussian + Chebyshev mass fit on a ROOT file and produce the fit
plot in one step:

    fit_gauss → results_gauss.json → gauss_mass_fit.png

Usage:
    python gauss_driver.py input.root [--outdir out]
                          [--xmin 480] [--xmax 620] [--nbins 80]
                          [--mean-init 548] [--sigma-init 10]
"""

import sys
import os
import json
import argparse

from fit_gauss import fit_gauss
from plot_gauss import CONFIG as PLOT_CONFIG
from plot_mass_fit import plot_fit_data

parser = argparse.ArgumentParser(
    description="Fit and plot a Gaussian + Chebyshev mass fit in one step",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
parser.add_argument("infile", help="Input ROOT file")
parser.add_argument("--outdir", default="out",
                    help="Output directory for results_gauss.json and gauss_mass_fit.png")
parser.add_argument("--xmin", type=float, default=480.0, help="Lower mass bound [MeV]")
parser.add_argument("--xmax", type=float, default=620.0, help="Upper mass bound [MeV]")
parser.add_argument("--nbins", type=int, default=80, help="Number of histogram bins")
parser.add_argument("--mean-init", type=float, default=547.9, help="Initial mean [MeV]")
parser.add_argument("--sigma-init", type=float, default=10.0, help="Initial sigma [MeV]")
args = parser.parse_args()

result = fit_gauss(
    args.infile,
    xmin=args.xmin, xmax=args.xmax, nbins=args.nbins,
    mean_init=args.mean_init, sigma_init=args.sigma_init,
)

os.makedirs(args.outdir, exist_ok=True)
json_path = os.path.join(args.outdir, "results_gauss.json")
with open(json_path, "w") as f:
    json.dump(result, f, indent=2)
print(f"[done] JSON saved to: {json_path}", file=sys.stderr)

plot_config = dict(PLOT_CONFIG)
plot_config["output"] = os.path.join(args.outdir, "gauss_mass_fit.png")
plot_fit_data(result, plot_config)
