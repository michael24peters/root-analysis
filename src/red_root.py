###############################################################################
# Step 1 of 4                                                                 #
# Script to make a reduced root file from many large, sparse root files.      #
# Author: Michael Peters                                                      #
###############################################################################
'''This is step one of four in generating plots from ntuple data. This script
reads in multiple ntuple root files, extracts all non-empty events and fills
all information into a new ntuple root file.
'''

import ROOT
import os

pre = '/data/home/michael24peters/anaroot/ntuple/MC_2018_MinBias_100M/'
infiles = [
    pre + 'magdown/00334331_00000001_1.etamumugamma.root',
    pre + 'magdown/00334331_00000002_1.etamumugamma.root',
    pre + 'magup/00334330_00000001_1.etamumugamma.root',
    pre + 'magup/00334330_00000002_1.etamumugamma.root'
]
outfile = 'red/reduced.root'

# Ensure output directory exists
os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)

print(f'Reading from {len(infiles)} files:')
for f in infiles:
    print(f'  - {f}')
print(f'Writing to: {outfile}')

# Create TChain from all input files
chain = ROOT.TChain('tree')
for file in infiles:
    chain.Add(file)

# Create reduced TFile and TTree
tfile = ROOT.TFile.Open(outfile, "RECREATE")
tree = chain.CloneTree(0)  # structure of original tree only

# Loop variables
check_interval = 1000000  # print status every n events
branch_names = ['tag_pid', 'prt_pid', 'mc_pid']  # empty event indicators

# Loop over all entries in chain and fill only non-empty events
for entryIdx in range(0, chain.GetEntries()):
    # Print status
    if entryIdx % check_interval == 0 and entryIdx > 0:
        print(f'  - Processed {entryIdx:,d} events, kept {tree.GetEntries()}...')

    chain.GetEntry(entryIdx)

    # Check if event is empty by looking at all branches
    is_empty = True
    for branch_name in branch_names:
        branch = getattr(chain, branch_name)
        if len(branch) > 0:
            is_empty = False
            break
    if is_empty: continue

    # print("Event passed selection. Filling event...")  # debug

    # Fill all branches for this entry
    tree.Fill()

# Write to TFile
tree.Write()
tfile.Close()

print(f'Processed {chain.GetEntries()} events, kept {tree.GetEntries()}...')
print(f'Done: wrote tree to {outfile}.')
