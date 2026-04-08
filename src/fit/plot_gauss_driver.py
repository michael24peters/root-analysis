"""
Plot Gaussian fit results from JSON file
Author: Michael Peters
Usage: python plot_gauss_driver.py fit_result.json
"""

from plot_mass_fit import plot_fit, val_err

def build_text(params, chi2_per_ndof):
    mean,  mean_err  = val_err(params["mean"])
    sigma, sigma_err = val_err(params["sigma"])
    plot_text = (
        f"$\\mu$ = {mean:.2f} ± {mean_err:.2f} MeV\n"
        f"$\\sigma$ = {sigma:.2f} ± {sigma_err:.2f} MeV\n"
        f"$\\chi^2/ndof$ = {chi2_per_ndof:.2f}"
    )
    term_text = (
        f"μ = {mean:.2f} ± {mean_err:.2f} MeV\n"
        f"σ = {sigma:.2f} ± {sigma_err:.2f} MeV\n"
        f"χ²/ndof = {chi2_per_ndof:.2f}"
    )
    return plot_text, term_text

plot_fit({
    "description"    : "Plot Gaussian fit results from JSON file",
    "output"         : "out/gauss_mass_fit.png",
    "title"          : "Eta Mass",
    # "xlim"           : (420, 680),
    # "pull_ylim"      : (-3, 3),
    "build_text"     : build_text,
})