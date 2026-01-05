###############################################################################
# Step 3 of 4                                                                 #
# Script to make reconstructed histograms from ntuple and save to root file.  #
# Author: Michael Peters                                                      #
###############################################################################

import ROOT
from utils.create_histograms import create_histograms
import sys

sig_file = False
if 'sig' in sys.argv[1:]:
    sig_file = True

if sig_file:
    infile = 'MC_2018_Signal/eta2MuMuGamma_mc_20251208.root'
    outfile = 'hist/sig_hist_rec.root'
else:
    infile = 'red/reduced_fiducial_cuts.root'
    outfile = 'hist/hist_rec.root'

print(f'Reading from {infile}, writing to {outfile}:')

# Arrays to hold histogram data
arr_tag_pid, arr_prt_pid = [], []
arr_tag_pt, arr_prt_pt = [], []
arr_tag_p, arr_prt_p = [], []
arr_tag_pz, arr_prt_pz = [], []
arr_tag_m = []
# Counters
ntag, nprt = 0, 0

# Combine files to create single histogram
tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

# Event loop
for entryIdx in range(0, tree.GetEntries()):
    tree.GetEntry(entryIdx)

    tag_pid = getattr(tree, 'tag_pid')  # type vector<double>
    tag_p = [getattr(tree, 'tag_px'),
                getattr(tree, 'tag_py'),
                getattr(tree, 'tag_pz')]
    tag_e = getattr(tree, 'tag_e')
    prt_pid = getattr(tree, 'prt_pid')
    prt_p = [getattr(tree, 'prt_px'),
                getattr(tree, 'prt_py'),
                getattr(tree, 'prt_pz')]
    
    # Skip empty events
    if len(tag_pid) == 0: continue

    # Compute momentum
    tag_mom = [ (px**2 + py**2 + pz**2)**0.5 for px, py, pz in zip(tag_p[0], tag_p[1], tag_p[2]) ]
    prt_mom = [ (px**2 + py**2 + pz**2)**0.5 for px, py, pz in zip(prt_p[0], prt_p[1], prt_p[2]) ]

    # Compute transverse momentum
    tag_pt = [ (px**2 + py**2)**0.5 for px, py in zip(tag_p[0], tag_p[1]) ]
    prt_pt = [ (px**2 + py**2)**0.5 for px, py in zip(prt_p[0], prt_p[1]) ]

    # Extract z component of momentum
    tag_pz = tag_p[2]
    prt_pz = prt_p[2]

    # Compute mass
    tag_m = [ (e**2 - (px**2 + py**2 + pz**2))**0.5 for e, px, py, pz in zip(tag_e, tag_p[0], tag_p[1], tag_p[2]) ]

    # Extract and fill tag information
    for pid, pt, pz, mom, m in zip(tag_pid, tag_pt, tag_pz, tag_mom, tag_m):
        arr_tag_pid.append(float(pid))
        arr_tag_pt.append(float(pt))
        arr_tag_pz.append(float(pz))
        arr_tag_p.append(float(mom))
        arr_tag_m.append(float(m))
        ntag += 1

    # Extract and fill particle information
    for pid, pt, pz, mom in zip(prt_pid, prt_pt, prt_pz, prt_mom):
        arr_prt_pid.append(float(pid))
        arr_prt_pt.append(float(pt))
        arr_prt_pz.append(float(pz))
        arr_prt_p.append(float(mom))
        nprt += 1

# Close TFile
tfile.Close()

# Create histogram variables
binwidths = [1, 1,  # pid bins
             10,  # mass bins (MeV)
             1000, 100, 1000,  # momentum, pt, pz bins (MeV)
             1000, 100, 1000]  # momentum, pt, pz bins (MeV)
arrays = [arr_tag_pid, arr_prt_pid, 
          arr_tag_m, 
          arr_tag_p, arr_tag_pt, arr_tag_pz,
          arr_prt_p, arr_prt_pt, arr_prt_pz]
names = ['tag_pid', 'prt_pid',
         'tag_m',
         'tag_p', 'tag_pt', 'tag_pz',
         'prt_p', 'prt_pt', 'prt_pz']
         
# Create histograms and save to output file
create_histograms(outfile, binwidths, arrays, names)

print('Number of reconstructed tags = ', ntag)
print('Number of reconstructed daughters = ', nprt)
print(f'Done: wrote histograms to {outfile}')
