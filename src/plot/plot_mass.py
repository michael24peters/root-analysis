################################################################################
# Script to make mass plots from root file.                                    #
# Usage: python plot_mass.py input.root [output.png]                           #
# Author: Michael Peters                                                       #
################################################################################

import argparse
import uproot
import numpy as np
import matplotlib.pyplot as plt

# IO
parser = argparse.ArgumentParser(description='Plot mass distributions from ROOT file')
parser.add_argument("input", 
                    help="Input ROOT file")
parser.add_argument("output", nargs="?", default="mass_plot.png",
                    help="Output PNG file (default: mass_plot.png)")
args = parser.parse_args()
print(f"Reading from {args.input} and writing to {args.output}.")

# Get mass info as numpy array
with uproot.open(args.input) as f:
    tree = f["tree"]
    # Convert to numpy array
    tag_m = tree["tag_dtf_m"].array(library="np") 

# Get lengths of arrays in tag_dtf_m to see how many candidates per event
lengths = np.array([len(x) for x in tag_m])
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
print(f"{'─' * 48}")

# Flatten array of arrays
tag_m = np.concatenate(tag_m)  # flatten array of arrays
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
plt.hist(tag_m, bins=200, range=(450, 600), histtype='step', color='black', label='Candidates')
plt.axvline(x=547.86, color='tab:red', linestyle='--', label='Eta')
# plt.axvline(x=770., color='tab:green', linestyle='--', label='Rho')
# plt.axvline(x=782.65, color='tab:pink', linestyle='--', label='Omega')
# plt.axvline(x=957.78, color='tab:blue', linestyle='--', label='Eta Prime')
# plt.axvline(x=1019.46, color='tab:purple', linestyle='--', label='Phi')
# plt.xlim(300, 700)
plt.xlabel('Mass [MeV]')
plt.ylabel('Candidates')
plt.title('Mass Distribution of Candidates')
plt.grid(alpha=0.3)
plt.legend(frameon=False)
plt.savefig(f"{args.output}", dpi=300)
print(f"Saved plot to {args.output}")
