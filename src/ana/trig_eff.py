"""
trig_eff.py

Counts TIS and TOS candidates and calculates the efficiency of the trigger
selection. Writes individual counts and efficiency to the terminal (or a
json file, optionally).

Usage:
    python trig_eff.py input.root [out/output.json]
"""

import argparse
import uproot
import awkward as ak
import numpy as np
from pathlib import Path

# IO
parser = argparse.ArgumentParser(description='Calculate trigger efficiency from ROOT file')
parser.add_argument('input', help='Input ROOT file')
parser.add_argument('output', default=None,
                    help='Output JSON file (optional)')
args = parser.parse_args()

# Establish output path for json file if specified. Should always write to out/
# directory, and the output argument is the filename to save in out/
if args.output:
    output_path = Path(args.output)
    print(f'[INFO] Reading from {args.input}, writing to {output_path}.')
else: print(f'[INFO] Reading from {args.input}.')

# Get trigger info as numpy array
with uproot.open(args.input) as f:
    tree = f['tree']
    # Load as awkward arrays (branches are jagged: variable candidates per event),
    # then reduce to one boolean per event with ak.any(..., axis=1)
    # L0
    l0_tos0 = ak.any(tree['tag_l0_tos0'].array(library='ak'), axis=1)
    l0_tos1 = ak.any(tree['tag_l0_tos1'].array(library='ak'), axis=1)
    l0_tis0 = ak.any(tree['tag_l0_tis0'].array(library='ak'), axis=1)
    l0_tis1 = ak.any(tree['tag_l0_tis1'].array(library='ak'), axis=1)
    print(f'[INFO] Loaded {len(l0_tos0)} events for L0.')
    # Hlt1
    hlt1_tos0 = ak.any(tree['tag_hlt1_tos0'].array(library='ak'), axis=1)
    hlt1_tos1 = ak.any(tree['tag_hlt1_tos1'].array(library='ak'), axis=1)
    hlt1_tis0 = ak.any(tree['tag_hlt1_tis0'].array(library='ak'), axis=1)
    hlt1_tis1 = ak.any(tree['tag_hlt1_tis1'].array(library='ak'), axis=1)
    print(f'[INFO] Loaded {len(hlt1_tos0)} events for Hlt1.')
    # Hlt2
    hlt2_tos0 = ak.any(tree['tag_hlt2_tos0'].array(library='ak'), axis=1)
    hlt2_tos1 = ak.any(tree['tag_hlt2_tos1'].array(library='ak'), axis=1)
    hlt2_tos2 = ak.any(tree['tag_hlt2_tos2'].array(library='ak'), axis=1)
    hlt2_tis0 = ak.any(tree['tag_hlt2_tis0'].array(library='ak'), axis=1)
    hlt2_tis1 = ak.any(tree['tag_hlt2_tis1'].array(library='ak'), axis=1)
    print(f'[INFO] Loaded {len(hlt2_tos0)} events for Hlt2.')

def calc_efficiency(tos, tis):
    """Calculate trigger efficiency using the formula:
    eff = N(eta2mumu | reconstructed and TIS and TOS) / N(eta2mumu | reconstructed and TIS)
    where N is the number of events satisfying the conditions.
    """
    n_tis = np.sum(tis)
    n_tos_given_tis = np.sum(tos & tis)
    eff = n_tos_given_tis / n_tis if n_tis > 0 else 0
    return eff, n_tos_given_tis, n_tis

def print_counts(tos, tis, trigger_name):
    """Print the number of TIS and TOS events for a given trigger."""
    print(f"{'─' * 60}")
    print(f'Trigger Counts for {trigger_name}:')
    print(f'TOS: {np.sum(tos)} / {len(tos)} events')
    print(f'TIS: {np.sum(tis)} / {len(tis)} events')

def print_efficiency(eff, n_tos_given_tis, n_tis, trigger_name):
    """Print the trigger efficiency for a given trigger."""
    print(f"{'─' * 60}")
    print(f'Trigger Efficiency for {trigger_name}:')
    print(f'TIS events: {n_tis}')
    print(f'TOS & TIS events: {n_tos_given_tis}')
    print(f'Efficiency: {eff:.4f}')

# L0 calculations
eff0_l0, n_tos_given_tis_0_l0, n_tis0_l0 = calc_efficiency(l0_tos0, l0_tis0)
eff1_l0, n_tos_given_tis_1_l0, n_tis1_l0 = calc_efficiency(l0_tos1, l0_tis1)
total_eff_l0, n_tos_given_tis_total_l0, n_tis_total_l0 = calc_efficiency(l0_tos0 | l0_tos1, l0_tis0 | l0_tis1)
print_counts(l0_tos0, l0_tis0, 'L0DiMuonDecision')  # tos0
print_counts(l0_tos1, l0_tis1, 'L0MuonDecision')  # tos1
print_efficiency(eff0_l0, n_tos_given_tis_0_l0, n_tis0_l0, 'L0DiMuonDecision')  # tos0
print_efficiency(eff1_l0, n_tos_given_tis_1_l0, n_tis1_l0, 'L0MuonDecision')  # tos1
print_efficiency(total_eff_l0, n_tos_given_tis_total_l0, n_tis_total_l0, 'L0')  # total

# Hlt1 calculations
eff0_hlt1, n_tos_given_tis_0_hlt1, n_tis0_hlt1 = calc_efficiency(hlt1_tos0, hlt1_tis0)
eff1_hlt1, n_tos_given_tis_1_hlt1, n_tis1_hlt1 = calc_efficiency(hlt1_tos1, hlt1_tis1)
total_eff_hlt1, n_tos_given_tis_total_hlt1, n_tis_total_hlt1 = calc_efficiency(hlt1_tos0 | hlt1_tos1,
                                                                               hlt1_tis0 | hlt1_tis1)
print_counts(hlt1_tos0, hlt1_tis0, 'Hlt1DiMuonNoIPDecision')  # tos0
print_counts(hlt1_tos1, hlt1_tis1, 'Hlt1DiMuonLowMassDecision')  # tos1
print_efficiency(eff0_hlt1, n_tos_given_tis_0_hlt1, n_tis0_hlt1, 'Hlt1DiMuonNoIPDecision')  # tos0
print_efficiency(eff1_hlt1, n_tos_given_tis_1_hlt1, n_tis1_hlt1, 'Hlt1DiMuonLowMassDecision')  # tos1
print_efficiency(total_eff_hlt1, n_tos_given_tis_total_hlt1, n_tis_total_hlt1, 'Hlt1')  # total

# Hlt2 calculations
eff0, n_tos_given_tis_0, n_tis0 = calc_efficiency(hlt2_tos0, hlt2_tis0)
eff1, n_tos_given_tis_1, n_tis1 = calc_efficiency(hlt2_tos1, hlt2_tis1)
total_eff_hlt2, n_tos_given_tis_total_hlt2, n_tis_total_hlt2 = calc_efficiency(hlt2_tos0 | hlt2_tos1, 
                                                                hlt2_tis0 | hlt2_tis1)
print_counts(hlt2_tos0, hlt2_tis0, 'Hlt2ExoticaPrmptDiMuonTurbo')  # tos0
print_counts(hlt2_tos1, hlt2_tis1, 'Hlt2ExoticaDiMuonNoIPTurbo')  # tos1
print_counts(hlt2_tos0 | hlt2_tos1, hlt2_tis0 | hlt2_tis1, 'Hlt2')  # total
print_efficiency(eff0, n_tos_given_tis_0, n_tis0, 'Hlt2ExoticaPrmptDiMuonTurbo')  # tos0
print_efficiency(eff1, n_tos_given_tis_1, n_tis1, 'Hlt2ExoticaDiMuonNoIPTurbo')  # tos1
print_efficiency(total_eff_hlt2, n_tos_given_tis_total_hlt2, n_tis_total_hlt2, 'Hlt2')  # total

# ------------------------------------------------------------------------------

# Write results to JSON if requested
if args.output:
    import json
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results = {
        'L0': {
            'TOS0': int(np.sum(l0_tos0)),
            'TIS0': int(np.sum(l0_tis0)),
            'TOS1': int(np.sum(l0_tos1)),
            'TIS1': int(np.sum(l0_tis1)),
            'efficiency_TOS0_given_TIS0': float(eff0_l0),
            'efficiency_TOS1_given_TIS1': float(eff1_l0),
            'overall_efficiency': float(total_eff_l0),
        },
        'HLT1': {
            'TOS0': int(np.sum(hlt1_tos0)),
            'TIS0': int(np.sum(hlt1_tis0)),
            'TOS1': int(np.sum(hlt1_tos1)),
            'TIS1': int(np.sum(hlt1_tis1)),
            'efficiency_TOS0_given_TIS0': float(eff0_hlt1),
            'efficiency_TOS1_given_TIS1': float(eff1_hlt1),
            'overall_efficiency': float(total_eff_hlt1),
        },
        'HLT2': {
            'TOS0': int(np.sum(hlt2_tos0)),
            'TIS0': int(np.sum(hlt2_tis0)),
            'TOS1': int(np.sum(hlt2_tos1)),
            'TIS1': int(np.sum(hlt2_tis1)),
            'efficiency_TOS0_given_TIS0': float(eff0),
            'efficiency_TOS1_given_TIS1': float(eff1),
            'overall_efficiency': float(total_eff_hlt2),
        }
    }
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=4)
    print(f'Wrote results to {output_path}.')
