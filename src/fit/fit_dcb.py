"""
fit_dcb.py

Fit an asymmetric Double Crystal Ball (DCB) signal + 2nd-order Chebyshev background
to the eta candidate mass distribution from a ROOT file using iminuit ExtendedBinnedNLL.

Tail parameters for left and right sides are independent: α_L, n_L, α_R, n_R.

Outputs fit results as JSON (stdout by default, or --outfile file.json).

Usage:
    python fit_dcb.py input.root [--outfile file.json]
                               [--xmin 480] [--xmax 620] [--nbins 80]
                               [--mean-init 548] [--sigma-init 10]
                               [--alphaL-init 1.5] [--nL-init 2]
                               [--alphaR-init 1.5] [--nR-init 2]
"""

import sys
import os
import json
import argparse
import numpy as np

# Ensure utils/fit_utils.py is importable
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', 'utils'))
import fit_utils


def fit_dcb(infile, xmin=480.0, xmax=620.0, nbins=80,
            mean_init=547.9, sigma_init=10.0,
            alphaL_init=1.5, nL_init=2.0,
            alphaR_init=1.5, nR_init=2.0):
    """Run the asymmetric DCB + Chebyshev fit and return the result dict."""
    if not os.path.exists(infile):
        sys.exit(f"[error] File not found: {infile}")

    centers, counts, errors = fit_utils.load_histogram(
        infile, xmin=xmin, xmax=xmax, nbins=nbins)
    bin_edges = np.linspace(xmin, xmax, nbins + 1)
    N_CAND = counts.sum()
    print(f"[info] {int(N_CAND)} entries in [{xmin}, {xmax}] MeV across {nbins} bins",
          file=sys.stderr)

    # Cost function
    cost = fit_utils.make_dcb_cost(counts, bin_edges, xmin, xmax)  # asymmetric version

    # Initial parameters and limits
    init = dict(
        mean   = mean_init,
        sigma  = sigma_init,
        alphaL = alphaL_init,
        nL     = nL_init,
        alphaR = alphaR_init,
        nR     = nR_init,
        n_sig  = N_CAND * 0.5,
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
        n_sig  = (0.0, N_CAND),
        c0     = (-1.0, 1.0),
        c1     = (-1.4, 1.4),
    )

    # Run fit
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
    print(f"[fit] α_L={p['alphaL']:.3f}±{e['alphaL']:.3f}  n_L={p['nL']:.3f}±{e['nL']:.3f}",
          file=sys.stderr)
    print(f"[fit] α_R={p['alphaR']:.3f}±{e['alphaR']:.3f}  n_R={p['nR']:.3f}±{e['nR']:.3f}",
          file=sys.stderr)

    # Compute bin predictions
    bw = (xmax - xmin) / nbins
    sig_cdf_edges = fit_utils.dcb_cdf(bin_edges, p['mean'], p['sigma'],
                                    p['alphaL'], p['nL'], p['alphaR'], p['nR'],
                                    xmin, xmax)
    bkg_cdf_edges = fit_utils.cheb_cdf(bin_edges, xmin, xmax, p['c0'], p['c1'])
    bin_sig = np.diff(sig_cdf_edges) * p['n_sig']
    bin_bkg = np.diff(bkg_cdf_edges) * n_bkg
    bin_fit = bin_sig + bin_bkg

    # Goodness of fit
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

    # Smooth curves
    x_curve = np.linspace(xmin, xmax, 200)
    curve_sig = fit_utils.dcb_pdf(x_curve, p['mean'], p['sigma'],
                                 p['alphaL'], p['nL'], p['alphaR'], p['nR'],
                                 xmin, xmax) * p['n_sig'] * bw
    curve_bkg = fit_utils.cheb_bkg_pdf(x_curve, xmin, xmax, p['c0'], p['c1']) * n_bkg * bw
    curve_fit = curve_sig + curve_bkg

    # Build JSON
    result = {
        "meta": {
            "input_file": infile,
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
                "n_bkg": {"value": round(float(n_bkg), 2),
                          "error": round(float(n_bkg_err), 2),
                          "note": "derived: n_cand - n_sig (not a free fit parameter)"},
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
    return result


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Asymmetric Double Crystal Ball + Chebyshev mass fit → JSON",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("infile", help="Input ROOT file")
    parser.add_argument("--outfile", default=None, help="Output JSON file (default: stdout)")
    parser.add_argument("--xmin", type=float, default=480.0, help="Lower mass bound [MeV]")
    parser.add_argument("--xmax", type=float, default=620.0, help="Upper mass bound [MeV]")
    parser.add_argument("--nbins", type=int, default=80, help="Number of histogram bins")
    parser.add_argument("--mean-init", type=float, default=547.9, help="Initial mean [MeV]")
    parser.add_argument("--sigma-init", type=float, default=10.0, help="Initial sigma [MeV]")
    parser.add_argument("--alphaL-init", type=float, default=1.5, help="Initial left tail α_L")
    parser.add_argument("--nL-init", type=float, default=2.0, help="Initial left tail n_L")
    parser.add_argument("--alphaR-init", type=float, default=1.5, help="Initial right tail α_R")
    parser.add_argument("--nR-init", type=float, default=2.0, help="Initial right tail n_R")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result = fit_dcb(
        args.infile,
        xmin=args.xmin, xmax=args.xmax, nbins=args.nbins,
        mean_init=args.mean_init, sigma_init=args.sigma_init,
        alphaL_init=args.alphaL_init, nL_init=args.nL_init,
        alphaR_init=args.alphaR_init, nR_init=args.nR_init,
    )

    json_str = json.dumps(result, indent=2)
    if args.outfile:
        args.outfile = 'out/' + args.outfile if not args.outfile.startswith('out/') else args.outfile
        outdir = os.path.dirname(os.path.abspath(args.outfile))
        os.makedirs(outdir, exist_ok=True)
        with open(args.outfile, "w") as f:
            f.write(json_str)
        print(f"[done] JSON saved to: {args.outfile}", file=sys.stderr)
    else:
        print(json_str)
