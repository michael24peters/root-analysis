# Counts TIS and TOS candidates and calculates the efficiency of the trigger
# selection. Writes individual counts and efficiency to the terminal (or a 
# json file, optionally).)
# Usage: python src/analysis/trig_eff.py input.root

import argparse
import uproot
import awkward as ak
import numpy as np
from pathlib import Path

# IO
parser = argparse.ArgumentParser(description='Calculate trigger efficiency from ROOT file')
parser.add_argument("input", help="Input ROOT file")
parser.add_argument("--output", "-o", default=None,
                    help="Output JSON file (optional)")
args = parser.parse_args()

# Establish output path for json file if specified. Should always write to out/
# directory, and the output argument is the filename to save in out/
if args.output:
    output_path = Path("out") / args.output
    print(f"Reading from {args.input}, writing to {output_path}.")
else: print(f"Reading from {args.input}.")

# Get trigger info as numpy array. Triggers:
# - L0: tag_l0_tos0, tag_l0_tos1
# - HLT1: tag_hlt1_tos0, tag_hlt1_tos1
# - HLT2: tag_hlt2_tos0, tag_hlt2_tos1, tag_hlt_tos2, tag_hlt_tis
with uproot.open(args.input) as f:
    tree = f["tree"]
    # Load as awkward arrays (branches are jagged: variable candidates per event),
    # then reduce to one boolean per event with ak.any(..., axis=1)
    hlt2_tos0 = ak.any(tree["tag_hlt2_tos0"].array(library="ak"), axis=1)
    hlt2_tos1 = ak.any(tree["tag_hlt2_tos1"].array(library="ak"), axis=1)
    hlt2_tos2 = ak.any(tree["tag_hlt2_tos2"].array(library="ak"), axis=1)
    hlt2_tis  = ak.any(tree["tag_hlt2_tis"].array(library="ak"),  axis=1)

# Trigger efficiency is defined as number of TOS given TIS divided by total
# number of TIS events, i.e., eff = N(TOS & TIS) / N(TIS). Calculate for each
# trigger.
effs = {}
print(f"{'─' * 48}")
print(f"Trigger Efficiencies:")
print(f"{'─' * 48}")

n_tis = np.sum(hlt2_tis)
print(f"TIS events: {n_tis}")
n_tos_given_tis_0 = np.sum(hlt2_tos0 & hlt2_tis)
eff0 = n_tos_given_tis_0 / n_tis if n_tis > 0 else 0
print(f"TOS0 & TIS: {n_tos_given_tis_0} (efficiency: {eff0:.4f})")
n_tos_given_tis_1 = np.sum(hlt2_tos1 & hlt2_tis)
eff1 = n_tos_given_tis_1 / n_tis if n_tis > 0 else 0

print(f"TOS1 & TIS: {n_tos_given_tis_1} (efficiency: {eff1:.4f})") 
n_tos_given_tis_2 = np.sum(hlt2_tos2 & hlt2_tis)
eff2 = n_tos_given_tis_2 / n_tis if n_tis > 0 else 0
print(f"TOS2 & TIS: {n_tos_given_tis_2} (efficiency: {eff2:.4f})")

# Then calculate the overall efficiency of the trigger, which is
# eff = N((L0 & HLT1 & HLT2) TOS & TIS) / N((L0 & HLT1 & HLT2) TIS).
overall_tos = (hlt2_tos0 | hlt2_tos1 | hlt2_tos2)
n_overall_tos_given_tis = np.sum(overall_tos & hlt2_tis)
overall_eff = n_overall_tos_given_tis / n_tis if n_tis > 0 else 0
print(f"{'─' * 48}")
print(f"Overall Trigger Efficiency (L0 & HLT1 & HLT2): {overall_eff:.4f} ({n_overall_tos_given_tis} / {n_tis})")
print(f"{'─' * 48}")