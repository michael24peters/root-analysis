###############################################################################
# Step 2 of 4                                                                 #
# Script to apply LHCb requirements on generator-level particles.             #
# Author: Michael Peters                                                      #
###############################################################################

import ROOT
import array
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
    def_outfile = 'root/sig_reduced_fiducial_cuts.root'
else:
    infile = 'red/reduced.root'
    def_outfile = 'red/reduced_fiducial_cuts.root'

outfile = ('red' + args.outfile) if args.outfile else def_outfile
print(f'Reading from {infile}, writing to {outfile}.')

tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

new_tfile = ROOT.TFile.Open(outfile, "RECREATE")
new_tfile.cd()
new_tree = tree.CloneTree(0)

# Create new branch to indicate which gen particles failed fiducial cuts
mc_req_pass = ROOT.std.vector('int')()
new_tree.Branch('mc_req_pass', mc_req_pass)

def pseudorapidity(px, py, pz):
    """
    Calculate pseudorapidity using same method as ROOT.TVector3::PseudoRapidity(),
    except using try / except to catch division by zero.
    Source: https://root.cern.ch/doc/master/classTVector3.html#aedc6fc6f5f6f3f3d4e4f4e4f4f4f4f4
    Code: https://root.cern.ch/doc/master/TVector3_8cxx_source.html#l00345
    """
    p = (px**2 + py**2 + pz**2)**0.5
    cosTheta = pz / p if p != 0 else 1.0  # Avoid division by zero
    if cosTheta * cosTheta < 1:
        return -0.5 * ROOT.TMath.Log((1.0 - cosTheta) / (1.0 + cosTheta))
    if pz == 0: return 0.0
    elif pz > 0: return 1e10
    else: return -1e10

print(f'entries: {tree.GetEntries()}')
for entryIdx in range(0, tree.GetEntries()):
    # Print status
    check_interval = 100000
    if entryIdx % check_interval == 0 and entryIdx > 0:
        print(f'  - Processed {entryIdx:,d} events, kept {new_tree.GetEntries()}...')
    
    tree.GetEntry(entryIdx)
    
    prt_pid = getattr(tree, 'prt_pid')  # Reconstructed particle pids
    mc_pid = getattr(tree, 'mc_pid')  # MC-matched daughter pids
    px = getattr(tree, 'mc_px')  # MC-matched daughter px
    py = getattr(tree, 'mc_py')  # MC-matched daughter py
    pz = getattr(tree, 'mc_pz')  # MC-matched daughter pz 
    
    # Loop over mc_pid and apply cuts based on pid
    has_rec = len(prt_pid) > 0
    has_mc = len(mc_pid) > 0
    passed = True
    for i, pid in enumerate(mc_pid):
        pid = int(pid)
        p = (px[i]**2 + py[i]**2 + pz[i]**2)**0.5  # momentum
        pt = (px[i]**2 + py[i]**2)**0.5  # transverse momentum
        
        # Prevent division by zero (wouldn't pass cuts anyway)
        if p == 0: passed = False; break
        eta = pseudorapidity(px[i], py[i], pz[i])
        
        # Apply cuts
        if pid == 221:
            mc_req_pass.push_back(-1)
            continue
        elif abs(pid) == 13:
            # LHCb requirements
            passed = passed and (2.0 < eta < 4.5)
            # Reconstruction requirements
            # passed = passed and (pt > 500) and (p > 3000) and (2.0 < eta < 4.5)
        elif pid == 22:
            # LHCb requirements
            passed = passed and (2.0 < eta < 4.5)
            # Reconstruction requirements
            # passed = passed and (pt > 500) and (2.0 < eta < 4.5)

        # Save info on failed particles in event
        if passed: mc_req_pass.push_back(1)
        else: 
            mc_req_pass.push_back(0)

    if has_rec: new_tree.Fill()
    elif has_mc and passed: new_tree.Fill()
    
    # Clear vector for next event
    mc_req_pass.clear()

print(f'Total kept entries: {new_tree.GetEntries()}')

# Close input file
tfile.Close()

# Write new tree to output file 
new_tree.Write()
new_tfile.Close()

print(f'Done: wrote reduced tree with fiducial cuts to {outfile}.')