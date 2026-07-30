"""
trig_eff.py

Counts TIS and TOS candidates and calculates the efficiency of the trigger
selection. Writes individual counts and efficiency to the terminal (or a
json file, optionally).

Usage:
    python trig_eff.py input.root [--outfile out/output.json]
"""

import argparse
import uproot
import awkward as ak
import numpy as np
from pathlib import Path

# IO
parser = argparse.ArgumentParser(description='Calculate trigger efficiency from ROOT file')
parser.add_argument('input', help='Input ROOT file')
parser.add_argument('--outfile', default=None,
                    help='Output JSON file (optional)')
args = parser.parse_args()

# Establish output path for json file if specified. Should always write to out/
# directory, and the output argument is the filename to save in out/
if args.outfile:
    output_path = Path(args.outfile)
    print(f'Reading from {args.input}, writing to {output_path}.')
else: print(f'Reading from {args.input}.')

# Get trigger info as numpy array. Triggers:
# - L0: tag_l0_tos0, tag_l0_tos1
# - HLT1: tag_hlt1_tos0, tag_hlt1_tos1
# - HLT2: tag_hlt2_tos0, tag_hlt2_tos1, tag_hlt_tos2, tag_hlt_tis
with uproot.open(args.input) as f:
    tree = f['tree']
    # Load as awkward arrays (branches are jagged: variable candidates per event),
    # then reduce to one boolean per event with ak.any(..., axis=1)
    hlt2_tos0 = ak.any(tree['tag_hlt2_tos0'].array(library='ak'), axis=1)
    hlt2_tos1 = ak.any(tree['tag_hlt2_tos1'].array(library='ak'), axis=1)
    hlt2_tos2 = ak.any(tree['tag_hlt2_tos2'].array(library='ak'), axis=1)
    hlt2_tis0 = ak.any(tree['tag_hlt2_tis0'].array(library='ak'), axis=1)
    hlt2_tis1 = ak.any(tree['tag_hlt2_tis1'].array(library='ak'), axis=1)
    hlt2_tis2 = ak.any(tree['tag_hlt2_tis2'].array(library='ak'), axis=1)
    hlt2_tos2 = ak.any(tree['tag_hlt2_tos2'].array(library='ak'), axis=1)
    hlt2_tis_topo = ak.any(tree['tag_hlt2_tis_topo'].array(library='ak'), axis=1)
    hlt2_tos_topo = ak.any(tree['tag_hlt2_tos_topo'].array(library='ak'), axis=1)

print(f"{'─' * 48}")
print(f'Trigger Counts:')
print(f"{'─' * 48}")
# Print number of TIS and TOS events for each trigger. Use ak.sum to count
# True values, then put out of total number of events (len(hlt2_tos0)).
print(f'TOS0: {np.sum(hlt2_tos0)} / {len(hlt2_tos0)} events')
print(f'TIS0: {np.sum(hlt2_tis0)} / {len(hlt2_tis0)} events')
print(f'TOS1: {np.sum(hlt2_tos1)} / {len(hlt2_tos1)} events')
print(f'TIS1: {np.sum(hlt2_tis1)} / {len(hlt2_tis1)} events')
print(f'TOS2: {np.sum(hlt2_tos2)} / {len(hlt2_tos2)} events')
print(f'TIS2: {np.sum(hlt2_tis2)} / {len(hlt2_tis2)} events')
print(f'TOS Topo: {np.sum(hlt2_tos_topo)} / {len(hlt2_tos_topo)} events')
print(f'TIS Topo: {np.sum(hlt2_tis_topo)} / {len(hlt2_tis_topo)} events')

# Trigger efficiency is defined as number of TOS given TIS divided by total
# number of TIS events, i.e., eff = N(TOS & TIS) / N(TIS). Calculate for each
# trigger.
effs = {}
print(f"{'─' * 48}")
print(f'Trigger Efficiencies:')
print(f"{'─' * 48}")

n_tis0 = np.sum(hlt2_tis0)
n_tis1 = np.sum(hlt2_tis1)
n_tis2 = np.sum(hlt2_tis2)
n_tis_topo = np.sum(hlt2_tis_topo)
print(f'TIS events: {n_tis0} (TIS0), {n_tis1} (TIS1), {n_tis2} (TIS2), {n_tis_topo} (TIS Topo)')
n_tos_given_tis_0 = np.sum(hlt2_tos0 & hlt2_tis0)
eff0 = n_tos_given_tis_0 / n_tis0 if n_tis0 > 0 else 0
print(f'TOS0 & TIS: {n_tos_given_tis_0} (efficiency: {eff0:.4f})')
n_tos_given_tis_1 = np.sum(hlt2_tos1 & hlt2_tis1)
eff1 = n_tos_given_tis_1 / n_tis1 if n_tis1 > 0 else 0
print(f'TOS1 & TIS: {n_tos_given_tis_1} (efficiency: {eff1:.4f})') 
n_tos_given_tis_2 = np.sum(hlt2_tos2 & hlt2_tis2)
eff2 = n_tos_given_tis_2 / n_tis2 if n_tis2 > 0 else 0
print(f'TOS2 & TIS: {n_tos_given_tis_2} (efficiency: {eff2:.4f})')
n_tos_given_tis_topo = np.sum(hlt2_tos_topo & hlt2_tis_topo)
eff3 = n_tos_given_tis_topo / n_tis_topo if n_tis_topo > 0 else 0
print(f'TOS Topo & TIS: {n_tos_given_tis_topo} (efficiency: {eff3:.4f})')

# Then calculate the overall efficiency of the trigger, which is
# eff = N((L0 & HLT1 & HLT2) TOS & TIS) / N((L0 & HLT1 & HLT2) TIS).
overall_tos = (hlt2_tos0 | hlt2_tos1 | hlt2_tos2 | hlt2_tos_topo)
overall_tis = (hlt2_tis0 | hlt2_tis1 | hlt2_tis2 | hlt2_tis_topo)
n_tis = np.sum(overall_tis)
n_tos_given_tis = np.sum(overall_tos & overall_tis)
overall_eff = n_tos_given_tis / n_tis if n_tis > 0 else 0
print(f"{'─' * 48}")
print(f'Overall HLT2 Trigger Efficiency: {overall_eff:.4f} ({n_tos_given_tis} / {n_tis})')

# Repeat except without including TIS2/TOS2 (i.e., displaced trigger), since
# this one seems to have problems
overall_tos_no_2 = (hlt2_tos0 | hlt2_tos1 | hlt2_tos_topo)
overall_tis_no_2 = (hlt2_tis0 | hlt2_tis1 | hlt2_tis_topo)
n_tis_no_2 = np.sum(overall_tis_no_2)
overall_eff_no_2 = np.sum(overall_tos_no_2 & overall_tis_no_2)
print(f'Overall HLT2 Trigger Efficiency: {overall_eff_no_2 / n_tis_no_2 if n_tis_no_2 > 0 else 0:.4f} ({overall_eff_no_2} / {n_tis_no_2})')
print(f'  (no TIS2/TOS2)')
print(f"{'─' * 48}")

# Write results to JSON if requested
if args.outfile:
    import json
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results = {
        'TOS0': int(np.sum(hlt2_tos0)),
        'TIS0': int(np.sum(hlt2_tis0)),
        'TOS1': int(np.sum(hlt2_tos1)),
        'TIS1': int(np.sum(hlt2_tis1)),
        'TOS2': int(np.sum(hlt2_tos2)),
        'TIS2': int(np.sum(hlt2_tis2)),
        'TOS Topo': int(np.sum(hlt2_tos_topo)),
        'TIS Topo': int(np.sum(hlt2_tis_topo)),
        'efficiency_TOS0_given_TIS0': float(eff0),
        'efficiency_TOS1_given_TIS1': float(eff1),
        'efficiency_TOS2_given_TIS2': float(eff2),
        'efficiency_TOS_Topo_given_TIS_Topo': float(eff3),
        'overall_efficiency': float(overall_eff),
        'overall_efficiency_no_2': float(overall_eff_no_2 / n_tis_no_2 if n_tis_no_2 > 0 else 0),
    }
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=4)
    print(f'Wrote results to {output_path}.')
