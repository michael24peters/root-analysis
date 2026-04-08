"""
fit_gauss.py
────────────
Fit a Gaussian signal + 2nd-order Chebyshev background to the eta candidate
mass distribution from a ROOT file, using iminuit ExtendedBinnedNLL.

Outputs fit results as JSON (stdout by default, or --output file.json).
The JSON is self-contained: it includes histogram data, smooth curve points,
and pull values so that plot_fit_result.py never needs the ROOT file.

Usage
─────
    python fit_gauss.py input.root [--output result.json]
                        [--xmin 480] [--xmax 620] [--nbins 80]
                        [--mean-init 548] [--sigma-init 10]

Pipeline mode (no intermediate file):
    python fit_gauss.py input.root | python plot_fit_result.py - output.png
"""

import sys
import os
import json
import argparse
import numpy as np

# Ensure src/utils is importable regardless of working directory
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', 'utils'))
import fitting

# ─── CLI ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Gaussian + Chebyshev mass fit → JSON",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
parser.add_argument("infile", help="Input ROOT file")
parser.add_argument("--output", default=None,
                    help="Output JSON file (default: stdout)")
parser.add_argument("--xmin",       type=float, default=480.0,
                    help="Lower mass bound [MeV] (default: 480)")
parser.add_argument("--xmax",       type=float, default=620.0,
                    help="Upper mass bound [MeV] (default: 620)")
parser.add_argument("--nbins",      type=int,   default=80,
                    help="Number of histogram bins (default: 80)")
parser.add_argument("--mean-init",  type=float, default=547.9,
                    help="Initial mean [MeV] (default: 547.9)")
parser.add_argument("--sigma-init", type=float, default=10.0,
                    help="Initial sigma [MeV] (default: 10.0)")
args = parser.parse_args()
xmin, xmax, nbins = args.xmin, args.xmax, args.nbins

# ─── Load data ────────────────────────────────────────────────────────────────
if not os.path.exists(args.infile):
    sys.exit(f"[error] File not found: {args.infile}")

centers, counts, errors = fitting.load_histogram(
    args.infile, xmin=xmin, xmax=xmax, nbins=nbins)
bin_edges = np.linspace(xmin, xmax, nbins + 1)
N = counts.sum()
print(f"[info] {int(N)} entries in [{xmin}, {xmax}] MeV across {nbins} bins",
      file=sys.stderr)

# ─── Fit ──────────────────────────────────────────────────────────────────────
cost = fitting.make_gauss_cost(counts, bin_edges, xmin, xmax)

init = dict(
    mean  = args.mean_init,
    sigma = args.sigma_init,
    n_sig = N * 0.5,
    n_bkg = N * 0.5,
    c0    = 0.0,
    c1    = 0.0,
)
limits = dict(
    mean  = (xmin, xmax),
    sigma = (1.0, 50.0),
    n_sig = (0.0, N * 2),
    n_bkg = (0.0, N * 2),
    c0    = (-1.0, 1.0),
    c1    = (-1.4, 1.4),
)

m = fitting.run_fit(cost, init, limits)
p, e = m.values, m.errors

print(f"[fit] valid={m.valid}  "
      f"mean={p['mean']:.3f}±{e['mean']:.3f} MeV  "
      f"σ={p['sigma']:.3f}±{e['sigma']:.3f} MeV",
      file=sys.stderr)

# ─── Bin predictions (CDF-based, matches fit model exactly) ───────────────────
bw = (xmax - xmin) / nbins

sig_cdf_edges = fitting.gauss_cdf(bin_edges, p['mean'], p['sigma'], xmin, xmax)
bkg_cdf_edges = fitting.cheb_cdf(bin_edges, xmin, xmax, p['c0'], p['c1'])
bin_sig = np.diff(sig_cdf_edges) * p['n_sig']
bin_bkg = np.diff(bkg_cdf_edges) * p['n_bkg']
bin_fit = bin_sig + bin_bkg

# Check cdf for sig
print("shape:", sig_cdf_edges.shape)
print("min/max:", sig_cdf_edges.min(), sig_cdf_edges.max())
print("monotonic:", np.all(np.diff(sig_cdf_edges) >= -1e-12))
# Repeat for bkg
print("shape:", bkg_cdf_edges.shape)
print("min/max:", bkg_cdf_edges.min(), bkg_cdf_edges.max())
print("monotonic:", np.all(np.diff(bkg_cdf_edges) >= -1e-12))

# ─── Goodness of fit (Pearson χ², matching RooFit's chiSquare convention) ─────
mask          = bin_fit > 0
chi2          = float(np.sum((counts[mask] - bin_fit[mask])**2 / bin_fit[mask]))
ndof          = int(mask.sum()) - m.nfit
chi2_per_ndof = chi2 / ndof if ndof > 0 else float('nan')

# Data-based pulls: (d − f) / √d  (matching RooFit pull plots; 0 for empty bins)
pulls = np.where(counts > 0,
                 (counts - bin_fit) / np.sqrt(counts),
                 0.0)

print(f"[fit] χ²/ndf = {chi2_per_ndof:.4f}  ({chi2:.2f} / {ndof})",
      file=sys.stderr)

# ─── Smooth curve (200 pts) for plotting ──────────────────────────────────────
x_curve   = np.linspace(xmin, xmax, 200)
curve_sig = fitting.gauss_pdf(x_curve, p['mean'], p['sigma'], xmin, xmax) * p['n_sig'] * bw
curve_bkg = fitting.cheb_bkg_pdf(x_curve, xmin, xmax, p['c0'], p['c1']) * p['n_bkg'] * bw
curve_fit = curve_sig + curve_bkg

# ─── Build JSON ───────────────────────────────────────────────────────────────
result = {
    "meta": {
        "input_file": args.infile,
        "model":      "gauss",
        "xmin":  xmin,
        "xmax":  xmax,
        "nbins": nbins,
        "symmetric": False,
    },
    "fit": {
        "converged":      bool(m.valid),
        "valid":          bool(m.valid),
        "fval":           round(float(m.fval), 4),
        "chi2":           round(chi2, 4),
        "ndof":           ndof,
        "chi2_per_ndof":  round(chi2_per_ndof, 6),
        "n_free_params":  m.nfit,
        "parameters": {
            "mean":  {"value": round(float(p['mean']),  6),
                      "error": round(float(e['mean']),  6), "unit": "MeV"},
            "sigma": {"value": round(float(p['sigma']), 6),
                      "error": round(float(e['sigma']), 6), "unit": "MeV"},
            "n_sig": {"value": round(float(p['n_sig']), 2),
                      "error": round(float(e['n_sig']), 2)},
            "n_bkg": {"value": round(float(p['n_bkg']), 2),
                      "error": round(float(e['n_bkg']), 2)},
            "c0":    {"value": round(float(p['c0']),    6),
                      "error": round(float(e['c0']),    6),
                      "note": "Chebyshev linear coefficient (slope)"},
            "c1":    {"value": round(float(p['c1']),    6),
                      "error": round(float(e['c1']),    6),
                      "note": "Chebyshev quadratic coefficient (curvature)"},
        },
    },
    "histogram": {
        "bin_centers": centers.tolist(),
        "bin_counts":  counts.tolist(),
        "bin_errors":  errors.tolist(),
        "bin_fit":     bin_fit.tolist(),
        "bin_sig":     bin_sig.tolist(),
        "bin_bkg":     bin_bkg.tolist(),
        "bin_pulls":   pulls.tolist(),
    },
    "curve": {
        "x":   x_curve.tolist(),
        "fit": curve_fit.tolist(),
        "sig": curve_sig.tolist(),
        "bkg": curve_bkg.tolist(),
    },
}

# ─── Output ───────────────────────────────────────────────────────────────────
json_str = json.dumps(result, indent=2)

if args.output:
    outdir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(outdir, exist_ok=True)
    with open(args.output, "w") as f:
        f.write(json_str)
    print(f"[done] JSON saved to: {args.output}", file=sys.stderr)
else:
    print(json_str)
