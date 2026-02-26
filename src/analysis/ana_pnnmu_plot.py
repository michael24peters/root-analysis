# Create plot of relative efficiency drop from baseline (probnn_mu > 0.4)
# So baseline should be 1.0 and the rest showing percent drop from that.
# Do this for signal and background, put on same plot vs probnn_mu cut value.
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import argparse
from pathlib import Path
import re

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile',
                    help = 'Input ROOT file')
parser.add_argument(
    '-d', '--dir',
    help='Directory with text outputs to parse'
)
args = parser.parse_args()
try: infile = args.infile
except: raise ValueError('Input file must be specified with -i flag.')
try: indir = args.dir
except: raise ValueError('Input directory must be specified with -d flag.')

# --- Get candidate and signal counts from text files ---
# Get path
results = {}
indir_path = Path(indir)
if not indir_path.is_dir():
    raise ValueError(f'Input directory not found: {indir}')

# Get text files
txt_files = sorted(indir_path.glob('*.txt'))
if not txt_files:
    raise ValueError(f'No .txt files found in directory: {indir}')

# Regular expressions to match candidate and signal lines
candidates_re = re.compile(r'^Total candidates processed:\s*(\d+)\s*$')
signal_re = re.compile(r'^Total signal candidates:\s*(\d+)\s*$')

# Parse each text file to extract candidate and signal counts
for txt_file in txt_files:
    candidates = None
    signal = None
    with txt_file.open('r', encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if candidates is None:
                match = candidates_re.match(line)
                if match:
                    candidates = int(match.group(1))
                    continue
            if signal is None:
                match = signal_re.match(line)
                if match:
                    signal = int(match.group(1))
                    continue

    # Store results
    results[txt_file.name] = {
        'candidates': candidates,
        'signal': signal,
    }

# Define arrays
pnnmu_cuts = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975])
can_counts, sig_counts = [], []

# Extract values
print('Parsed candidate/signal counts:')
for name, counts in results.items():
    print(f'{name}: candidates={counts["candidates"]}, signal={counts["signal"]}')
    can_counts.append(counts['candidates'])
    sig_counts.append(counts['signal'])

# Convert to numpy arrays for calculations
can_counts = np.array(can_counts)
sig_counts = np.array(sig_counts)
bkg_counts = can_counts - sig_counts

# Import ROOT file to get total number of gen-level candidates
tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

# Get total number of gen-level candidates
ngen = 0
for entryIdx in range(0, tree.GetEntries()):
    # NOTE: Count will be slightly off since in previous code we only allowed
    # one gen-level candidate per event, here we are counting all candidates.
    tree.GetEntry(entryIdx)
    mc_pid = getattr(tree, 'mc_pid')
    mc_pid = [int(pid) for pid in mc_pid]
    for i, pid in enumerate(mc_pid):
        if pid == 221: ngen += 1

# Absolute rates
sig_abs_rate = sig_counts / can_counts
bkg_abs_rate = bkg_counts / can_counts

# Relative rates
sig_rel_rate = sig_counts / sig_counts[0]
bkg_rel_rate = bkg_counts / bkg_counts[0]

# Efficiencies
eff_abs = can_counts / ngen
sig_eff_abs = sig_counts / ngen

# Relative efficiencies
eff_rel = eff_abs / eff_abs[0]
sig_eff_rel = sig_eff_abs / sig_eff_abs[0]

# Plot sig_rel_rate and bkg_rel_rate vs pnnmu_cuts
plt.figure(figsize=(10, 6))
plt.plot(pnnmu_cuts, sig_rel_rate, marker='o', color='blue', label='Signal Relative Rate')
plt.plot(pnnmu_cuts, bkg_rel_rate, marker='o', color='black', label='Background Relative Rate')
plt.title('Relative Signal and Background Rates vs PROBNNmu Cut')
plt.xlabel('PROBNNmu')
plt.ylabel('Relative Rate')
plt.ylim(0, 1.1)
plt.grid(True)
plt.legend()
plt.savefig('out/pnnmu_cuts.png')
plt.close()
print('pnnmu_cuts.png written.')

plt.figure(figsize=(10, 6))
plt.plot(pnnmu_cuts, eff_rel, marker='o', color='black', label='Efficiency Relative')
plt.plot(pnnmu_cuts, sig_eff_rel, marker='o', color='blue', label='Signal Efficiency Relative')
plt.title('Relative Efficiencies vs PROBNNmu Cut')
plt.xlabel('PROBNNmu')
plt.ylabel('Relative Efficiency')
plt.ylim(0, 1.1)
plt.grid(True)
plt.legend()
plt.savefig('out/pnnmu_cuts_eff.png')
plt.close()
print('pnnmu_cuts_eff.png written.')