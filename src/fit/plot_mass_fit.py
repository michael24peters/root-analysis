"""
plot_fit.py  -  Library for plotting fit results from JSON files.

Usage: python plot_mass_fit.py input.json [output.png]
"""

import argparse
import json
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description='Plot fit results from JSON file')
    parser.add_argument("input", help="Input JSON file")
    return parser.parse_args()

def load_data(infile):
    with open(infile) as f: return json.load(f)

def val_err(param):
    return param.get("value", np.nan), param.get("error", np.nan)

def print_summary(text):
    W = 48
    print()
    print("─" * W)
    print(f"{'Fit summary':^{W-4}}")
    print("─" * W)
    print(text)
    print("─" * W)

def plot_fit(config):
    # Parse arguments
    args = parse_args()
    output = config["output"] if "output" in config else "out/fit_result.png"
    print(f"[INFO] Reading from {args.input} and writing to {output}.")

    # Load JSON data
    data = load_data(args.input)

    # Extract histogram data
    hist = data["histogram"]
    x = np.array(hist["bin_centers"])
    y = np.array(hist["bin_counts"])
    yerr = np.array(hist["bin_errors"])
    fit = np.array(hist["bin_fit"])
    sig = np.array(hist["bin_sig"])
    bkg = np.array(hist["bin_bkg"])
    pulls = np.array(hist["bin_pulls"])
    # LaTeX font
    plt.rcParams["font.family"] = "STIXGeneral"
    plt.rcParams["font.size"] = 16
    # Create figure with two panels
    fig, (ax, ax_pull) = plt.subplots(
        2, 1, figsize=(12, 12),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True
    )
    # Title and label
    ax.set_title(config["title"])

    # Top panel: data + fit
    # Set x-axis limits based on config or data range
    if config.get("x_lim") is not None: xmin, xmax = config["x_lim"]
    else: xmin, xmax = x.min() - 0.5, x.max() + 0.5
    ax.set_xlim(xmin, xmax)
    ax.errorbar(x, y, yerr=yerr,
                fmt='o',
                label="Data",
                color="black",
                markersize=5,
                markeredgewidth=1,
                elinewidth=0.8,
                capsize=2
    )
    # Plot fit, signal, background
    ax.plot(x, fit, label="Total fit", color="red")
    ax.plot(x, sig, '--', label="Signal", color="blue")
    ax.plot(x, bkg, '--', label="Background", color="green")
    ax.axvline(x=547.86, color='black', linestyle='--', label='PDG Mass')

    # Grid, legend, and label
    ax.grid(alpha=0.3)
    ax.legend(frameon=False)
    ax.set_ylabel("Candidates")

    # Add fit parameters to plot
    params = data["fit"]["parameters"]
    chi2_per_ndof = data["fit"]["chi2_per_ndof"]
    # Get text from config function, with error handling for compatibility
    try: 
        plot_text, term_text = config["build_text"](params, chi2_per_ndof)
    except: 
        plot_text = config["build_text"](params, chi2_per_ndof)
        term_text = None
    ax.text(0.025, 0.95, plot_text,
            transform=ax.transAxes,
            verticalalignment='top',
    )

    # Bottom panel: pulls
    bin_width = x[1] - x[0]
    ax_pull.bar(x, pulls, width=bin_width, color="red", edgecolor='none')
    ax_pull.axhline(0, color='black', lw=0.8)
    # y axis config
    # Set y axis limits based on config or pull range
    if config.get("pull_ylim") is not None: ymin, ymax = config["pull_ylim"]
    else: ymin, ymax = int(np.floor(pulls.min())) - 1, int(np.ceil(pulls.max())) + 1
    ax_pull.set_ylim(ymin, ymax)
    ax_pull.set_yticks(np.arange(ymin, ymax + 1, 1))  # TODO: test ticks
    ax_pull.tick_params(axis='y', direction='in')
    # Grid and labels
    ax_pull.grid(alpha=0.3)
    ax_pull.set_ylabel("Pulls")
    ax_pull.set_xlabel("Mass (MeV/c²)")
    # Show x labels on both panels, comment out for only bottom panel
    # ax.tick_params(axis='x', which='both', labelbottom=True)

    # Save plot
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    print(f"[DONE] Plot saved to {output}.")
    # Print summary
    if term_text: print_summary(term_text)
    else: print_summary(plot_text)
