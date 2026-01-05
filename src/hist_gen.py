###############################################################################
# Step 3 of 4                                                                 #
# Script to make gen-level histograms from ntuple and save to root file.      #
# Author: Michael Peters                                                      #
###############################################################################

import ROOT
from utils.create_histograms import create_histograms
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '-o', '--outfile',
    help='Output ROOT file'
)
parser.add_argument(
    '-s', '--sig',
    action='store_true',
    help='Use signal file'
)

args = parser.parse_args()

sig_file = args.sig
if 'sig' in sys.argv[1:]:
    sig_file = True

if sig_file:
    infile = 'MC_2018_Signal/eta2MuMuGamma_mc_20251121.root'
    def_outfile = 'hist/sig_hist_gen.root'
else:
    infile = 'red/reduced_fiducial_cuts.root'
    def_outfile = 'hist/hist_gen.root'

outfile = ('hist' + args.outfile) if args.outfile else def_outfile
print(f'Reading from {infile}, writing to {outfile}.')

# Arrays to hold histogram data
arr_mc_pid = []
arr_mc_pt = []
arr_mc_p = []
arr_mc_pz = []
arr_mc_m = []
# Counters
ntag = 0

# Combine files to create single histogram
tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

# Event loop
for entryIdx in range(0, tree.GetEntries()):
    tree.GetEntry(entryIdx)

    # Extract gen-level tag and particle information
    mc_pid = getattr(tree, 'mc_pid')  # type vector<double>
    mc_px = getattr(tree, 'mc_px')
    mc_py = getattr(tree, 'mc_py')
    mc_pz = getattr(tree, 'mc_pz')
    mc_e = getattr(tree, 'mc_e')
    
    # Skip empty events
    if len(mc_pid) == 0: continue

    for i, pid in enumerate(mc_pid):
        pid = int(pid)
        
        # Compute momentum
        p = (mc_px[i]**2 + mc_py[i]**2 + mc_pz[i]**2)**0.5
        
        # Compute transverse momentum
        pt = (mc_px[i]**2 + mc_py[i]**2)**0.5
        
        # Compute mass for eta only
        if pid == 221: m = (mc_e[i]**2 - p**2)**0.5
        
        # Fill arrays
        arr_mc_pid.append(float(pid))
        arr_mc_pt.append(float(pt))
        arr_mc_p.append(float(p))
        arr_mc_pz.append(float(mc_pz[i]))
        if pid == 221: arr_mc_m.append(float(m))
        ntag += 1

# Close TFile
tfile.Close()

# Create histogram variables
binwidths = [1,  # pid bins
             1,  # mass bins (MeV)
             2000, 100, 2000]  # momentum, pt, pz bins (MeV)
arrays = [arr_mc_pid, 
          arr_mc_m, 
          arr_mc_p, arr_mc_pt, arr_mc_pz]
names = ['mc_pid',
         'mc_m',
         'mc_p', 'mc_pt', 'mc_pz']

# Create histograms and save to file
create_histograms(outfile, binwidths, arrays, names)

print('Number of generated tags = ', ntag)
print(f'Done: wrote histograms to {outfile}')
