"""
fit_gauss.py

Fit a Gaussian signal + 2nd-order Chebyshev background to the eta candidate
mass distribution from a ROOT file, using iminuit ExtendedBinnedNLL.

Outputs fit results as JSON (stdout by default, or --output file.json).
The JSON is self-contained: it includes histogram data, smooth curve points,
and pull values so that plot_fit_result.py never needs the ROOT file.

Usage:
    python fit_gauss.py input.root [--output result.json]
                        [--xmin 480] [--xmax 620] [--nbins 80]
                        [--mean-init 548] [--sigma-init 10]
"""

import sys
import os
import json
import argparse
import numpy as np

# Ensure src/utils is importable regardless of working directory
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', 'utils'))
import fit_utils


def fit_gauss(infile, xmin=480.0, xmax=620.0, nbins=80,
              mean_init=547.9, sigma_init=10.0):
    """Run the Gaussian + Chebyshev fit and return the result dict."""
    if not os.path.exists(infile):
        sys.exit(f"[error] File not found: {infile}")

    centers, counts, errors = fit_utils.load_histogram(
        infile, xmin=xmin, xmax=xmax, nbins=nbins)
    bin_edges = np.linspace(xmin, xmax, nbins + 1)
    N_CAND = counts.sum()
    print(f"[info] {int(N_CAND)} entries in [{xmin}, {xmax}] MeV across {nbins} bins",
          file=sys.stderr)

    # Fit
    cost = fit_utils.make_gauss_cost(counts, bin_edges, xmin, xmax)

    init = dict(
        mean  = mean_init,
        sigma = sigma_init,
        n_sig = N_CAND * 0.5,
        c0    = 0.0,
        c1    = 0.0,
    )
    limits = dict(
        mean  = (xmin, xmax),
        sigma = (1.0, 50.0),
        n_sig = (0.0, N_CAND),
        c0    = (-1.0, 1.0),
        c1    = (-1.4, 1.4),
    )

    m = fit_utils.run_fit(cost, init, limits)
    p, e = m.values, m.errors

    # n_bkg is not a free parameter: n_cand = N_CAND is fixed, and n_sig + n_bkg = N_CAND
    # by unitarity, so n_bkg is derived from n_sig. Since N_CAND carries no
    # uncertainty, n_bkg's error equals n_sig's error (n_bkg = N_CAND - n_sig).
    n_bkg = N_CAND - p['n_sig']
    n_bkg_err = e['n_sig']

    print(f"[fit] valid={m.valid}  "
          f"mean={p['mean']:.3f}±{e['mean']:.3f} MeV  "
          f"σ={p['sigma']:.3f}±{e['sigma']:.3f} MeV",
          file=sys.stderr)

    # Bin predictions
    bw = (xmax - xmin) / nbins

    sig_cdf_edges = fit_utils.gauss_cdf(bin_edges, p['mean'], p['sigma'], xmin, xmax)
    bkg_cdf_edges = fit_utils.cheb_cdf(bin_edges, xmin, xmax, p['c0'], p['c1'])
    bin_sig = np.diff(sig_cdf_edges) * p['n_sig']
    bin_bkg = np.diff(bkg_cdf_edges) * n_bkg
    bin_fit = bin_sig + bin_bkg

    # Check cdf for sig
    print("sig cdf shape:", sig_cdf_edges.shape)
    print("sig cdf min/max:", sig_cdf_edges.min(), sig_cdf_edges.max())
    print("sig cdf monotonic:", np.all(np.diff(sig_cdf_edges) >= -1e-12))
    # Repeat for bkg
    print("bkg cdf shape:", bkg_cdf_edges.shape)
    print("bkg cdf min/max:", bkg_cdf_edges.min(), bkg_cdf_edges.max())
    print("bkg cdf monotonic:", np.all(np.diff(bkg_cdf_edges) >= -1e-12))

    # Goodness of fit
    mask = bin_fit > 0
    chi2 = float(np.sum((counts[mask] - bin_fit[mask])**2 / bin_fit[mask]))
    ndof = int(mask.sum()) - m.nfit
    chi2_per_ndof = chi2 / ndof if ndof > 0 else float('nan')

    # Data-based pulls: (d − f) / √d  (matching RooFit pull plots; 0 for empty bins)
    pulls = np.where(counts > 0,
                     (counts - bin_fit) / np.sqrt(counts),
                     0.0)

    print(f"[fit] χ²/ndf = {chi2_per_ndof:.4f}  ({chi2:.2f} / {ndof})",
          file=sys.stderr)

    # Smooth curve
    x_curve   = np.linspace(xmin, xmax, 200)
    curve_sig = fit_utils.gauss_pdf(x_curve, p['mean'], p['sigma'], xmin, xmax) * p['n_sig'] * bw
    curve_bkg = fit_utils.cheb_bkg_pdf(x_curve, xmin, xmax, p['c0'], p['c1']) * n_bkg * bw
    curve_fit = curve_sig + curve_bkg

    # Build JSON
    result = {
        "meta": {
            "input_file": infile,
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
                "n_bkg": {"value": round(float(n_bkg), 2),
                          "error": round(float(n_bkg_err), 2),
                          "note": "derived: n_cand - n_sig (not a free fit parameter)"},
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
    return result


def _parse_args():
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
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result = fit_gauss(
        args.infile,
        xmin=args.xmin, xmax=args.xmax, nbins=args.nbins,
        mean_init=args.mean_init, sigma_init=args.sigma_init,
    )

    json_str = json.dumps(result, indent=2)
    if args.output:
        outdir = os.path.dirname(os.path.abspath(args.output))
        os.makedirs(outdir, exist_ok=True)
        with open(args.output, "w") as f:
            f.write(json_str)
        print(f"[done] JSON saved to: {args.output}", file=sys.stderr)
    else:
        print(json_str)
