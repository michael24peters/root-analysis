###############################################################################
# Step 3 of 4                                                                 #
# Script to make mass histograms from reduced ntuple and save to root file.   #
# Author: Michael Peters                                                      #
###############################################################################
'''
- Inside one event, we want to make sure prt_idx_mom == 0 for all daughters.
- Then we want to make sure prt_idx_gen points to the correct MC pids:
    - The 0 row in mc_pid should be 221.
    - For mu-: prt_idx_gen should point to mc_pid == 13, and its
    corresponding mc_idx_mom should == 0. 
    - For mu+: prt_idx_gen should point to mc_pid == -13, and its
    corresponding mc_idx_mom should == 0.
    - For gamma: prt_idx_gen should point to mc_pid == 22, and its
    corresponding mc_idx_mom should == 0.
'''

import ROOT
from utils.create_histograms import create_histograms
import sys

sig_file = False
if 'sig' in sys.argv[1:]: sig_file = True

if sig_file:
    infile = 'MC_2018_Signal/eta2MuMuGamma_mc_20251208.root'
    outfile = 'hist/sig_hist_m.root'
else:
    # infile = 'red/reduced.root'
    infile = 'red/reduced_fiducial_cuts.root'
    outfile = 'hist/hist_m.root'

print(f'Reading from {infile}, writing to {outfile}.')

arr_sig, arr_bkg, arr_tot = [], [], []  # arrays for signal, background, total
nsig, nbkg, ntot = 0, 0, 0  # counters for signal, background, total
ntags, ncan = 0, 0  # debug counters

# Combine files to create single histogram
tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

# Event loop
for entryIdx in range(0, tree.GetEntries()):
    tree.GetEntry(entryIdx)

    # Reconstructed tag info
    tag_pid = getattr(tree, 'tag_pid')  # eta pid, type vector<double>
    tag_m = getattr(tree, 'tag_m')  # eta mass

    # Reconstructed daughter info
    prt_pid = getattr(tree, 'prt_pid')  # daughter pids
    # MC-match info for daughters to MC-matched daughters (-1 if no match)
    prt_idx_gen = getattr(tree, 'prt_idx_gen')
    # MC-match info for daughters to candidate tag (-1 if mother not candidate
    # gen-level eta)
    prt_idx_mom = getattr(tree, 'prt_idx_mom')
    
    # Gen-level MCParticle info
    mc_pid = getattr(tree, 'mc_pid')  # all mc particle pids
    # MC-match info for gen-level index to its mother (-1 if mother is not eta
    # candidate)
    mc_idx_mom = getattr(tree, 'mc_idx_mom')

    # Reformat above lists for easier handling
    tag_m_new = []
    if len(tag_m) > 0: 
        for i in range(len(tag_m)): tag_m_new.append(float(tag_m[i]))
        tag_m = tag_m_new
    prt_idx_gen = [int(idx) for idx in prt_idx_gen]
    prt_idx_mom = [int(idx) for idx in prt_idx_mom]
    mc_idx_mom = [int(idx) for idx in mc_idx_mom]

    # Skip empty events
    ntags = len(tag_pid)
    if ntags == 0: continue
    ncan += ntags

    # print(f'tag_pid: {tag_pid}, prt_pid: {prt_pid}, mc_pid: {mc_pid}')  # debug

    # Loop over all candidates in event. Use indexing information to find
    # correctly MC-matched events vs incorrectly MC-matched events.
    # TODO: create sub-categories of background for mis-matched photons vs 
    # mis-matched muons vs mis-matched dimuon pairs.
    for i in range(ntags):
        is_signal = True
        # Check that there is a MC eta candidate. If not, we immediately know this
        # is background.
        # TODO: This doesn't handle multiple eta candidates (very rare).
        if mc_pid[0] != 221: is_signal = False
        if tag_pid[i] != 221: continue  # skip failed reco/non-eta candidates
        # Loop for each daughters of candidate
        for j in range(i*3, i*3+3):
            if prt_idx_mom[j] != i:  break  # Skip failed reco, shouldn't happen
            
            if prt_idx_gen[j] == -1: is_signal = False  # no MC match
            elif mc_pid[prt_idx_gen[j]] != prt_pid[j]: is_signal = False  # wrong MC match
            elif mc_idx_mom[prt_idx_gen[j]] != i: is_signal = False  # wrong mother match
            if not is_signal: break  # TODO: continue loop, classify background type
    
        # Fill arrays and increment counters
        # Append flattened elements of tag_m to arrays
        arr_tot.append(tag_m[i]); ntot += 1
        if is_signal: arr_sig.append(tag_m[i]); nsig += 1
        else: arr_bkg.append(tag_m[i]); nbkg += 1

# Print summary statistics
print(f'Number of events processed: {tree.GetEntries()}')
print(f'Number of signal candidates: {nsig}')
print(f'Number of background candidates: {nbkg}')
print(f'Number of total candidates: {ntot}')

# Close TFile
tfile.Close()

# Create histogram variables
binwidth = 10  # MeV
arrays = [arr_sig, arr_bkg, arr_tot]
names = ['sig', 'bkg', 'tot']
binwidths = [binwidth] * len(arrays)

# Create histograms and save to ROOT file
create_histograms(outfile, binwidths, arrays, names)
print(f'Done: wrote histograms to {outfile}')
