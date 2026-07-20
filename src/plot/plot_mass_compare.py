################################################################################
# Script to compare "all candidates" vs "best candidate per event" mass plots  #
# on one overlaid plot.                                                        #
# Usage: python plot_mass_compare.py input.root [output.png]                   #
# Author: Michael Peters                                                       #
################################################################################

import argparse
from anaroot.src.utils import fit_utils
import uproot
import numpy as np
import matplotlib.pyplot as plt

# IO
parser = argparse.ArgumentParser(description='Compare all-candidate vs best-candidate mass distributions')
parser.add_argument("input",
                    help="Input ROOT file")
parser.add_argument("output", nargs="?", default="mass_plot_compare.png",
                    help="Output PNG file (default: mass_plot_compare.png)")
parser.add_argument("--full-range", action="store_true", default=False,
                    help="Plot full mass range (400-1200 MeV) instead of zoomed-in range (450-600 MeV) (default: False)")
args = parser.parse_args()
print(f"Reading from {args.input} and writing to {args.output}.")

# Get mass info as numpy arrays
with uproot.open(args.input) as f:
    tree = f["tree"]
    arr = tree["tag_dtf_m"].array(library="np")
    metric_arr = tree["tag_dtf_chi2"].array(library="np")

# All candidates, flattened
tag_m_all = np.concatenate(arr)

# Best candidate per event (lowest tag_dtf_chi2). find_best_candidate returns
# (None, None) for events with 0 candidates, so those are dropped.
best = [fit_utils.find_best_candidate(candidates, metrics=metrics, min=True)[0]
        for candidates, metrics in zip(arr, metric_arr)]
tag_m_best = np.array([x for x in best if x is not None])

print(f"{'─' * 48}")
print(f"Candidates plotted (all):  {len(tag_m_all)}")
print(f"Candidates plotted (best): {len(tag_m_best)}")
print(f"Difference (candidates dropped by best selection): {len(tag_m_all) - len(tag_m_best)}")
print(f"{'─' * 48}")

# Plot
print("Plotting histogram...")
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 16
plt.figure(figsize=(12, 9))

bins = 200
if args.full_range: hist_range = (400, 1200)
else: hist_range = (450, 600)

# All candidates: red low-opacity fill + red full-opacity outline.
# Best candidate: blue low-opacity fill + blue full-opacity outline.
# Where the two fills overlap, the alpha blending reads as purple.
plt.hist(tag_m_all, bins=bins, range=hist_range, histtype='stepfilled',
         color='red', alpha=0.25, label='All candidates')
plt.hist(tag_m_all, bins=bins, range=hist_range, histtype='step',
         color='red', alpha=1.0, linewidth=1.5)

plt.hist(tag_m_best, bins=bins, range=hist_range, histtype='stepfilled',
         color='blue', alpha=0.25, label='Best candidate')
plt.hist(tag_m_best, bins=bins, range=hist_range, histtype='step',
         color='blue', alpha=1.0, linewidth=1.5)

plt.axvline(x=547.86, color='tab:green', linestyle='--', label='Eta')

if args.full_range:
    plt.axvline(x=770., color='tab:orange', linestyle='--', label='Rho')
    plt.axvline(x=782.65, color='tab:pink', linestyle='--', label='Omega')
    plt.axvline(x=957.78, color='tab:purple', linestyle='--', label='Eta Prime')
    plt.axvline(x=1019.46, color='tab:brown', linestyle='--', label='Phi')

plt.xlabel('Mass [MeV]')
plt.ylabel('Candidates')
plt.title('Mass Distribution: All vs Best Candidate')
plt.grid(alpha=0.3)
plt.legend(frameon=False)
plt.savefig(f"{args.output}", dpi=300)
print(f"Saved plot to {args.output}")
