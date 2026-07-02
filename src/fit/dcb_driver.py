"""
dcb_driver.py

Run the asymmetric DCB + Chebyshev mass fit on a ROOT file and produce the
fit plot in one step:

    fit_dcb → results_dcb.json → dcb_mass_fit.png

Usage:
    python dcb_driver.py input.root [--outdir out]
                       [--xmin 480] [--xmax 620] [--nbins 80]
                       [--mean-init 548] [--sigma-init 10]
                       [--alphaL-init 1.5] [--nL-init 2]
                       [--alphaR-init 1.5] [--nR-init 2]
"""

import sys
import os
import json
import argparse

from fit_dcb import fit_dcb
from plot_dcb import CONFIG as PLOT_CONFIG
from plot_mass_fit import plot_fit_data

parser = argparse.ArgumentParser(
    description="Fit and plot an asymmetric DCB + Chebyshev mass fit in one step",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
parser.add_argument("infile", help="Input ROOT file")
parser.add_argument("--outdir", default="out",
                    help="Output directory for results_dcb.json and dcb_mass_fit.png")
parser.add_argument("--xmin", type=float, default=480.0, help="Lower mass bound [MeV]")
parser.add_argument("--xmax", type=float, default=620.0, help="Upper mass bound [MeV]")
parser.add_argument("--nbins", type=int, default=80, help="Number of histogram bins")
parser.add_argument("--mean-init", type=float, default=547.9, help="Initial mean [MeV]")
parser.add_argument("--sigma-init", type=float, default=10.0, help="Initial sigma [MeV]")
parser.add_argument("--alphaL-init", type=float, default=1.5, help="Initial left tail α_L")
parser.add_argument("--nL-init", type=float, default=2.0, help="Initial left tail n_L")
parser.add_argument("--alphaR-init", type=float, default=1.5, help="Initial right tail α_R")
parser.add_argument("--nR-init", type=float, default=2.0, help="Initial right tail n_R")
args = parser.parse_args()

result = fit_dcb(
    args.infile,
    xmin=args.xmin, xmax=args.xmax, nbins=args.nbins,
    mean_init=args.mean_init, sigma_init=args.sigma_init,
    alphaL_init=args.alphaL_init, nL_init=args.nL_init,
    alphaR_init=args.alphaR_init, nR_init=args.nR_init,
)

os.makedirs(args.outdir, exist_ok=True)
json_path = os.path.join(args.outdir, "results_dcb.json")
with open(json_path, "w") as f:
    json.dump(result, f, indent=2)
print(f"[done] JSON saved to: {json_path}", file=sys.stderr)

plot_config = dict(PLOT_CONFIG)
plot_config["output"] = os.path.join(args.outdir, "dcb_mass_fit.png")
plot_fit_data(result, plot_config)
