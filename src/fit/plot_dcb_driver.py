"""
Plot DCB fit results from JSON file
Author: Michael Peters
Usage: python plot_dcb_driver.py fit_result.json
"""

from plot_mass_fit import plot_fit, val_err

def build_text(params, chi2_per_ndof):
    mean, mean_err   = val_err(params["mean"])
    sigma, sigma_err = val_err(params["sigma"])
    n_sig, n_sig_err = val_err(params["n_sig"])
    n_bkg, n_bkg_err = val_err(params["n_bkg"])
    c0, c0_err       = val_err(params["c0"])
    c1, c1_err       = val_err(params["c1"])

    # DCB tail parameters (asymmetric)
    alphaL, alphaL_err = val_err(params.get("alphaL", {}))
    nL, nL_err         = val_err(params.get("nL", {}))
    alphaR, alphaR_err = val_err(params.get("alphaR", {}))
    nR, nR_err         = val_err(params.get("nR", {}))

    plot_text = (
        f"$\\mu$ = {mean:.2f} ± {mean_err:.2f} MeV\n"
        f"$\\sigma$ = {sigma:.2f} ± {sigma_err:.2f} MeV\n"
        f"$\\alpha_L$ = {alphaL:.2f} ± {alphaL_err:.2f}, "
        f"$n_L$ = {nL:.1f} ± {nL_err:.1f}\n"
        f"$\\alpha_R$ = {alphaR:.2f} ± {alphaR_err:.2f}, "
        f"$n_R$ = {nR:.2f} ± {nR_err:.2f}\n"
        f"$n_{{sig}}$ = {n_sig:.0f} ± {n_sig_err:.0f}\n"
        f"$n_{{bkg}}$ = {n_bkg:.0f} ± {n_bkg_err:.0f}\n"
        f"$c_0$ = {c0:.2f} ± {c0_err:.2f}, $c_1$ = {c1:.2f} ± {c1_err:.2f}\n"
        f"$\\chi^2/ndof$ = {chi2_per_ndof:.2f}"
    )
    term_text = (
        f"μ = {mean:.2f} ± {mean_err:.2f} MeV\n"
        f"σ = {sigma:.2f} ± {sigma_err:.2f} MeV\n"
        f"α_L = {alphaL:.2f} ± {alphaL_err:.2f}, "
        f"n_L = {nL:.1f} ± {nL_err:.1f}\n"
        f"α_R = {alphaR:.2f} ± {alphaR_err:.2f}, "
        f"n_R = {nR:.2f} ± {nR_err:.2f}\n"
        f"n_sig = {n_sig:.0f} ± {n_sig_err:.0f}\n"
        f"n_bkg = {n_bkg:.0f} ± {n_bkg_err:.0f}\n"
        f"c_0 = {c0:.2f} ± {c0_err:.2f}, c_1 = {c1:.2f} ± {c1_err:.2f}\n"
        f"χ²/ndof = {chi2_per_ndof:.2f}"
    )
    return plot_text, term_text

plot_fit({
    "description"    : "Plot DCB fit results from JSON file",
    "output"         : "out/dcb_mass_fit.png",
    "title"          : "Eta Mass",
    # "xlim"           : (420, 680),
    # "pull_ylim"      : (-3, 3),
    "build_text"     : build_text,
})