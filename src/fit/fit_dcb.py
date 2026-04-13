"""
fit_dcb_asym.py
────────────────
Fit an asymmetric Double Crystal Ball (DCB) signal + 2nd-order Chebyshev background
to the eta candidate mass distribution from a ROOT file using iminuit ExtendedBinnedNLL.

Tail parameters for left and right sides are independent: α_L, n_L, α_R, n_R.

Outputs fit results as JSON (stdout by default, or --output file.json).

Usage
─────
    python fit_dcb_asym.py input.root [--output result.json]
                               [--xmin 480] [--xmax 620] [--nbins 80]
                               [--mean-init 548] [--sigma-init 10]
                               [--alphaL-init 1.5] [--nL-init 2]
                               [--alphaR-init 1.5] [--nR-init 2]

Pipeline mode:
    python fit_dcb_asym.py input.root | python plot_dcb.py - output.png
"""

import sys
import os
import json
import argparse
import numpy as np

# Ensure utils/fitting.py is importable
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', 'utils'))
import fitting

# ─── CLI ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Asymmetric Double Crystal Ball + Chebyshev mass fit → JSON",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
parser.add_argument("infile", help="Input ROOT file")
parser.add_argument("--output", default=None, help="Output JSON file (default: stdout)")
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

# ─── Build cost function ──────────────────────────────────────────────────────
cost = fitting.make_dcb_cost(counts, bin_edges, xmin, xmax)  # asymmetric version

# ─── Initial parameters and limits ───────────────────────────────────────────
init = dict(
    mean   = args.mean_init,
    sigma  = args.sigma_init,
    alphaL = args.alphaL_init,
    nL     = args.nL_init,
    alphaR = args.alphaR_init,
    nR     = args.nR_init,
    n_sig  = N * 0.5,
    n_bkg  = N * 0.5,
    c0     = 0.0,
    c1     = 0.0,
)
limits = dict(
    mean   = (xmin, xmax),
    sigma  = (1.0, 50.0),
    alphaL = (0.1, 10.0),
    nL     = (0.1, 50.0),
    alphaR = (0.1, 10.0),
    nR     = (0.1, 50.0),
    n_sig  = (0.0, N*2),
    n_bkg  = (0.0, N*2),
    c0     = (-1.0, 1.0),
    c1     = (-1.4, 1.4),
)

# ─── Run fit ──────────────────────────────────────────────────────────────────
m = fitting.run_fit(cost, init, limits)
p, e = m.values, m.errors

print(f"[fit] valid={m.valid}  "
      f"mean={p['mean']:.3f}±{e['mean']:.3f} MeV  "
      f"σ={p['sigma']:.3f}±{e['sigma']:.3f} MeV",
      file=sys.stderr)
print(f"[fit] α_L={p['alphaL']:.3f}±{e['alphaL']:.3f}  n_L={p['nL']:.3f}±{e['nL']:.3f}",
      file=sys.stderr)
print(f"[fit] α_R={p['alphaR']:.3f}±{e['alphaR']:.3f}  n_R={p['nR']:.3f}±{e['nR']:.3f}",
      file=sys.stderr)

# ─── Compute bin predictions ──────────────────────────────────────────────────
bw = (xmax - xmin) / nbins
sig_cdf_edges = fitting.dcb_cdf(bin_edges, p['mean'], p['sigma'],
                                p['alphaL'], p['nL'], p['alphaR'], p['nR'],
                                xmin, xmax)
bkg_cdf_edges = fitting.cheb_cdf(bin_edges, xmin, xmax, p['c0'], p['c1'])
bin_sig = np.diff(sig_cdf_edges) * p['n_sig']
bin_bkg = np.diff(bkg_cdf_edges) * p['n_bkg']
bin_fit = bin_sig + bin_bkg

# ─── Goodness of fit ─────────────────────────────────────────────────────────
mask = bin_fit > 0
chi2 = float(np.sum((counts[mask] - bin_fit[mask])**2 / bin_fit[mask]))
ndof = int(mask.sum()) - m.nfit
chi2_per_ndof = chi2 / ndof if ndof > 0 else float('nan')

# Returns pull or 0 if count is 0 for each bin.
# Pull: (y_data - y_fit) / sigma_data, where sigma_data = sqrt(y_data)
pulls = np.where(counts > 0,
                 (counts - bin_fit) / np.sqrt(counts),
                 0.0)

print(f"[fit] χ²/ndf = {chi2_per_ndof:.4f}  ({chi2:.2f} / {ndof})",
      file=sys.stderr)

# ─── Smooth curves ───────────────────────────────────────────────────────────
x_curve = np.linspace(xmin, xmax, 200)
curve_sig = fitting.dcb_pdf(x_curve, p['mean'], p['sigma'],
                             p['alphaL'], p['nL'], p['alphaR'], p['nR'],
                             xmin, xmax) * p['n_sig'] * bw
curve_bkg = fitting.cheb_bkg_pdf(x_curve, xmin, xmax, p['c0'], p['c1']) * p['n_bkg'] * bw
curve_fit = curve_sig + curve_bkg

# ─── Build JSON ─────────────────────────────────────────────────────────────
result = {
    "meta": {
        "input_file": args.infile,
        "model": "dcb_asym",
        "xmin": xmin,
        "xmax": xmax,
        "nbins": nbins,
        "symmetric": False,
    },
    "fit": {
        "converged": bool(m.valid),
        "valid": bool(m.valid),
        "fval": round(float(m.fval), 4),
        "chi2": round(chi2, 4),
        "ndof": ndof,
        "chi2_per_ndof": round(chi2_per_ndof, 6),
        "n_free_params": m.nfit,
        "parameters": {
            "mean":  {"value": round(float(p['mean']), 6),
                      "error": round(float(e['mean']), 6), "unit": "MeV"},
            "sigma": {"value": round(float(p['sigma']), 6),
                      "error": round(float(e['sigma']), 6), "unit": "MeV"},
            "alphaL":{"value": round(float(p['alphaL']), 6),
                      "error": round(float(e['alphaL']), 6)},
            "nL":    {"value": round(float(p['nL']), 6),
                      "error": round(float(e['nL']), 6)},
            "alphaR":{"value": round(float(p['alphaR']), 6),
                      "error": round(float(e['alphaR']), 6)},
            "nR":    {"value": round(float(p['nR']), 6),
                      "error": round(float(e['nR']), 6)},
            "n_sig": {"value": round(float(p['n_sig']), 2),
                      "error": round(float(e['n_sig']), 2)},
            "n_bkg": {"value": round(float(p['n_bkg']), 2),
                      "error": round(float(e['n_bkg']), 2)},
            "c0":    {"value": round(float(p['c0']), 6),
                      "error": round(float(e['c0']), 6),
                      "note": "Chebyshev linear coefficient (slope)"},
            "c1":    {"value": round(float(p['c1']), 6),
                      "error": round(float(e['c1']), 6),
                      "note": "Chebyshev quadratic coefficient (curvature)"},
        },
    },
    "histogram": {
        "bin_centers": centers.tolist(),
        "bin_counts": counts.tolist(),
        "bin_errors": errors.tolist(),
        "bin_fit": bin_fit.tolist(),
        "bin_sig": bin_sig.tolist(),
        "bin_bkg": bin_bkg.tolist(),
        "bin_pulls": pulls.tolist(),
    },
    "curve": {
        "x": x_curve.tolist(),
        "fit": curve_fit.tolist(),
        "sig": curve_sig.tolist(),
        "bkg": curve_bkg.tolist(),
    },
}

# ─── Output JSON ─────────────────────────────────────────────────────────────
json_str = json.dumps(result, indent=2)
if args.output:
    outdir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(outdir, exist_ok=True)
    with open(args.output, "w") as f:
        f.write(json_str)
    print(f"[done] JSON saved to: {args.output}", file=sys.stderr)
else:
    print(json_str)