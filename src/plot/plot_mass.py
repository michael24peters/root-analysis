################################################################################
# Script to make mass plots from root file.                                    #
# Usage: python plot_mass.py input.root [output.png]                           #
# Author: Michael Peters                                                       #
################################################################################

import argparse
from anaroot.src.utils import fit_utils
import uproot
import numpy as np
import matplotlib.pyplot as plt

# IO
parser = argparse.ArgumentParser(description='Plot mass distributions from ROOT file')
parser.add_argument("input", 
                    help="Input ROOT file")
parser.add_argument("output", nargs="?", default="mass_plot.png",
                    help="Output PNG file (default: mass_plot.png)")
# Flag to only plot "best" candidate (one candidate per event), 
# using fit_utils.find_best_candidate() to select the best candidate.
parser.add_argument("--best", action="store_true", default=False,
                    help="Only plot best candidate per event (default: False)")
# Flag to plot full mass range (400-1200 MeV) instead of zoomed-in range 
# (450-600 MeV).
parser.add_argument("--full-range", action="store_true", default=False,
                    help="Plot full mass range (400-1200 MeV) instead of zoomed-in range (450-600 MeV) (default: False)")
args = parser.parse_args()
print(f"Reading from {args.input} and writing to {args.output}.")

# Get mass info as numpy array
with uproot.open(args.input) as f:
    tree = f["tree"]
    # Raw, per-event jagged array of candidate masses (always needed for the
    # candidate-multiplicity diagnostics below, regardless of --best).
    arr = tree["tag_dtf_m"].array(library="np")

    # Get lengths of arrays in tag_dtf_m to see how many candidates per event
    lengths = np.array([len(x) for x in arr])
    print(f"{'─' * 48}")
    print(f"Events with 0 candidates: {np.sum(lengths == 0)}")
    print(f"Events with 1 candidate:  {np.sum(lengths == 1)}")
    # TODO: Any events with 2+ candidates need to be handled in the future, as
    # flattening will mix candidates from different events. For now, we just print
    # them out and don't worry about it, but in the future we will want to do "best
    # candidate" selection (e.g., the one closest to the eta mass) or somehow index
    # candidates by event to distinguish this.
    print(f"Events with 2+ candidates: {np.sum(lengths >= 2)}")
    print(f"Events with 3+ candidates: {np.sum(lengths >= 3)}")
    print(f"Max candidates in an event: {lengths.max()}")
    # Candidates plotted with --best=False is every candidate (lengths.sum());
    # with --best=True it's one per event that has at least one candidate.
    # Report both so the effect of --best is visible regardless of which
    # mode this run actually uses.
    n_all_candidates = int(lengths.sum())
    n_best_candidates = int(np.sum(lengths >= 1))
    print(f"Candidates plotted (--best=False): {n_all_candidates}")
    print(f"Candidates plotted (--best=True):  {n_best_candidates}")
    print(f"Difference (candidates dropped by --best selection): {n_all_candidates - n_best_candidates}")
    print(f"{'─' * 48}")

    if args.best:
        # Select the best candidate for each event
        # Use metric "tag_dtf_chi2" to select the best candidate (lowest chi2)
        metric_arr = tree["tag_dtf_chi2"].array(library="np")
        # find_best_candidate(candidates, metrics=None, min=True) -> best_candidate, best_idx
        # Returns (None, None) for events with 0 candidates, so those are
        # filtered out below rather than fed into the histogram.
        best = [fit_utils.find_best_candidate(candidates, metrics=metrics, min=True)[0]
                for candidates, metrics in zip(arr, metric_arr)]
        tag_m = np.array([x for x in best if x is not None])
    else:
        # Flatten array of arrays
        tag_m = np.concatenate(arr)

print(f"Read {len(tag_m)} entries from {args.input}")
print(f"shape: {tag_m.shape}")
print(f"dtype: {tag_m.dtype}")
print(f"min: {tag_m.min():.2f}")
print(f"max: {tag_m.max():.2f}")
print(f"{'─' * 48}")

# Plot
print("Plotting histogram...")
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 16
plt.figure(figsize=(12, 9))
if args.full_range: hist_range = (400, 1200)
else: hist_range = (450, 600)
plt.hist(tag_m, bins=200, range=hist_range, histtype='step', color='black', label='Candidates')
plt.axvline(x=547.86, color='tab:red', linestyle='--', label='Eta')
if args.full_range:
    plt.axvline(x=770., color='tab:orange', linestyle='--', label='Rho')
    plt.axvline(x=782.65, color='tab:pink', linestyle='--', label='Omega')
    plt.axvline(x=957.78, color='tab:purple', linestyle='--', label='Eta Prime')
    plt.axvline(x=1019.46, color='tab:brown', linestyle='--', label='Phi')
# plt.xlim(300, 700)
plt.xlabel('Mass [MeV]')
plt.ylabel('Candidates')
plt.title('Mass Distribution of Candidates')
plt.grid(alpha=0.3)
plt.legend(frameon=False)
plt.savefig(f"{args.output}", dpi=300)
print(f"Saved plot to {args.output}")
