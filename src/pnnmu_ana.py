# Create plot of relative efficiency drop from baseline (probnn_mu > 0.4)
# So baseline should be 1.0 and the rest showing percent drop from that.
# Do this for signal and background, put on same plot vs probnn_mu cut value.
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile',
                    help = 'Input ROOT file')
args = parser.parse_args()
try: infile = args.infile
except: raise ValueError('Input file must be specified with -i flag.')

tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

# Get total number of gen-level candidates
ngen = 0
for entryIdx in range(0, tree.GetEntries()):
    tree.GetEntry(entryIdx)
    mc_pid = getattr(tree, 'mc_pid')
    mc_pid = [int(pid) for pid in mc_pid]
    for i, pid in enumerate(mc_pid):
        if pid == 221: ngen += 1

pnnmu_cuts = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975])
can_counts = np.array([250, 181, 133, 113, 82, 42, 29, 11])
sig_counts = np.array([3, 2, 2, 2, 2, 2, 2, 0])
bkg_counts = can_counts - sig_counts

sig_abs_rate = sig_counts / can_counts
bkg_abs_rate = bkg_counts / can_counts
# eff = ncan / ngen
eff_abs = can_counts / ngen
# sig eff = nsig / ngen
sig_eff_abs = sig_counts / ngen

sig_rel_rate = sig_counts / sig_counts[0]
bkg_rel_rate = bkg_counts / bkg_counts[0]
eff_rel = eff_abs / eff_abs[0]
sig_eff_rel = sig_eff_abs / sig_eff_abs[0]

# print(f'pnnmu_cuts = {pnnmu_cuts}')  # debug
# print(f'can_counts = {can_counts}')  # debug
# print(f'sig_counts = {sig_counts}')  # debug
# print(f'bkg_counts = {bkg_counts}')  # debug
# print(f'sig_abs_rate = {sig_abs_rate}')  # debug
# print(f'bkg_abs_rate = {bkg_abs_rate}')  # debug
# print(f'eff_abs = {eff_abs}')  # debug
# print(f'sig_eff_abs = {sig_eff_abs}')  # debug
# print(f'eff_rel = {eff_rel}')  # debug
# print(f'sig_eff_rel = {sig_eff_rel}')  # debug

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