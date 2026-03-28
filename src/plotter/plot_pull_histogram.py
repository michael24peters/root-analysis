#!/usr/bin/env python3
"""
plot_pull_histogram.py  –  Histogram of per-bin pull values from a DCB mass fit
─────────────────────────────────────────────────────────────────────────────────
Reads a .npy or whitespace-delimited .txt file of pull values (one per bin)
produced by plot_dcb_mass_fit.py --dump-pulls, then plots a pull-value histogram
overlaid with a unit Gaussian N(0, 1) for reference.

Usage
─────
  python plot_pull_histogram.py pulls.npy              [options]
  python plot_pull_histogram.py pulls.txt              [options]

Options
───────
  --output PATH     Output PNG path           (default: out/pull_histogram.png)
  --nbins  INT      Histogram bins            (default: 20)
  --title  STR      Plot title                (default: "Pull distribution")
  --xmin   FLOAT    x-axis lower bound        (default: auto from data)
  --xmax   FLOAT    x-axis upper bound        (default: auto from data)

─────────────────────────────────────────────────────────────────────────────────
How to produce the pull file from plot_dcb_mass_fit.py
─────────────────────────────────────────────────────────────────────────────────
Add  --dump-pulls  to your fit command:

  python plot_dcb_mass_fit.py data.root fit.png --dump-pulls

This writes  out/pulls.npy  alongside the existing JSON and PNG outputs.
To add this flag, insert the following two lines into plot_dcb_mass_fit.py:

  # In the argparse block (near the other --output-json argument):
  g_out.add_argument("--dump-pulls", action="store_true",
                     help="Save per-bin pull values to out/pulls.npy")

  # After the `pulls` list is fully populated (after the per-bin loop):
  if args.dump_pulls:
      import numpy as np
      pulls_path = os.path.join("out", "pulls.npy")
      np.save(pulls_path, np.array(pulls))
      print(f"[done] Pull values saved to: {pulls_path}")
"""

import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
from scipy.optimize import curve_fit

# ─── CLI ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Histogram of pull values from a DCB mass fit.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument("pulls_file",
                    help="Pull values file: .npy (numpy) or .txt (whitespace-delimited)")
parser.add_argument("--output", default="out/pull_histogram.png",
                    help="Output PNG path (default: out/pull_histogram.png)")
parser.add_argument("--nbins",  type=int,   default=20,
                    help="Number of histogram bins (default: 20)")
parser.add_argument("--title",  default="Pull distribution",
                    help='Plot title (default: "Pull distribution")')
parser.add_argument("--xmin",   type=float, default=None,
                    help="x-axis lower bound (default: auto)")
parser.add_argument("--xmax",   type=float, default=None,
                    help="x-axis upper bound (default: auto)")
args = parser.parse_args()


# ─── Load pull values ─────────────────────────────────────────────────────────
if not os.path.exists(args.pulls_file):
    sys.exit(f"File not found: {args.pulls_file}")

pulls = np.load(args.pulls_file)
pulls = pulls.astype(float)

n_pulls = len(pulls)
print(f"[info] Loaded {n_pulls} pull values from: {args.pulls_file}")


# ─── Define Gaussian function ─────────────────────────────────────────────────
def gaussian(x, mu, sigma, amp):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ─── Get histogram values ─────────────────────────────────────────────────────
# Each pull is defined as (y_data - y_fit) / sigma_data for a given bin. If the
# fit is good then each pull is drawn from a normal distribution. There are
# three resulting cases:
# 1) sigma ~= 1: The fit is good.
# 2) sigma < 1: Overestimated or too free; too good of a fit.
# 3) sigma > 1: Underestimated or too constrained; poor fit.
mean = np.mean(pulls)

# Standard deviation defined by
# sigma = sqrt( (1/(N-1)) * sum((x_i - x_bar)^2) ), 
# where x_bar is the sample mean. 
#
# We're computing a sample std and are trying to estimate the true std of the 
# underlying distribution. We also only know the sample mean, not the true mean,
# so Bessel correction is used (where we multiply by N/(N-1) yielding 1/(N-1) 
# instead of 1/N) to get an unbiased estimator.
# 
std = np.std(pulls, ddof=1)

# x range for the plot
xmin = args.xmin if args.xmin is not None else min(-4.0, mean - 3.5 * std)
xmax = args.xmax if args.xmax is not None else max( 4.0, mean + 3.5 * std)

# Histogram the pulls to get bin counts for the Gaussian amplitude seed
# Returns the counts in each bin and the edges of the bins
counts, edges = np.histogram(pulls, bins=args.nbins, range=(xmin, xmax))
bin_width     = edges[1] - edges[0]

# ─── Plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))

# Histogram of pull values
ax.hist(
    pulls,
    bins=args.nbins,
    range=(xmin, xmax),
    color="#e07b39",   # orange — matches kOrange+7 from the fit script
    alpha=0.75,
    edgecolor="white",
    linewidth=0.7,
    label=f"Pull values  ($N_{{\\mathrm{{bins}}}}={n_pulls}$)",
    zorder=2,
)

# Overlay a unit Gaussian for reference
x_ref  = np.linspace(xmin, xmax, 400)
# Find the area under the histogram to scale the Gaussian PDF
# Each bar's area is count * bin_width, so total histogram area is 
# sum(counts) * bin_width. A normal gaussian PDF has area 1, so multiplying by
# the histogram's total area exactly scales it to match.
#
# If we assume the histogram of pull values should follow a normal distribution,
# then the area and standard deviation we get 
area   = n_pulls * bin_width
# calculates pdf for norm and scales it to match the histogram area
y_ref  = area * stats.norm.pdf(x_ref)
ax.plot(x_ref, y_ref,
        color="steelblue", linewidth=1.6, linestyle="--",
        label=r"Normal Gaussian (scaled to hist area)", zorder=3)

# Axis labels, title, legend
ax.set_xlabel(r"Pull  $(y_{\mathrm{data}} - y_{\mathrm{fit}})\,/\,\sigma_{\mathrm{data}}$",
              fontsize=12)
ax.set_ylabel("Bins", fontsize=12)
ax.set_title(args.title, fontsize=13)
ax.set_xlim(xmin, xmax)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.tick_params(which="both", direction="in", top=True, right=True)
ax.legend(fontsize=9, framealpha=0.9, frameon=False)
# Show mean and std in top-left corner of the plot
ax.text(
    0.03, 0.97,
    f"$\\mu = {mean:+.3f}$\n$\\sigma = {std:.3f}$",
    transform=ax.transAxes,      # coordinates are 0–1 in axes fraction
    ha="left", va="top",
    fontsize=9,
    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8),
)

fig.tight_layout()

# ─── Save ─────────────────────────────────────────────────────────────────────
os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
fig.savefig(args.output, dpi=150)
print(f"[done] Pull histogram saved to: {args.output}")

# ─── Print summary ────────────────────────────────────────────────────────────
W = 48
print()
print("─" * W)
print(f"  {'Pull distribution summary':^{W-4}}")
print("─" * W)
print(f"  N bins            : {n_pulls}")
print(f"  Bin width         : {bin_width:.4f}")
print(f"  Sample mean       : {mean:+.4f}")
print(f"  Sample std        : {std:.4f}")
print("─" * W)
