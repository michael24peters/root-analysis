################################################################################
# Script to apply fiducial requirements (also calculates efficiencies).        #
# Author: Michael Peters                                                       #
################################################################################
# TODO: make this work for eta -> 4mu and eta -> 2mu 2e cases as well. Requires
# more generalized code.

import ROOT
import argparse
import logging

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

def passes_reqs(pid, px, py, pz):
    """Check if a particle with given pid and momentum passes fiducial
    requirements.

    Fiducial requirements:
    - Pseudorapidity (eta) in [2, 4.5]
    - Muon pT > 500 MeV
    - Muon P > 3 GeV
    - Photon pT > 500 MeV
    """
    p = (px**2 + py**2 + pz**2)**0.5  # momentum
    pt = (px**2 + py**2)**0.5  # transverse momentum
    
    # Prevent division by zero (wouldn't pass cuts anyway)
    if p == 0: return False
    eta = pseudorapidity(px, py, pz)
    
    # Apply requirements
    if abs(pid) == 13: return (2.0 < eta < 4.5) and (pt > 500) and (p > 3000)
    elif pid == 22: return (pt > 500) and (2.0 < eta < 4.5)
    elif pid == 221: return True
    else: return False

#-------------------------------------------------------------------------------

def apply_fiducial_reqs(tree, gamma_flag):
    """Apply fiducial cuts to generator-level particles.
    
    Returns a new tree with only events that pass the fiducial cuts.
    """

    new_tree = tree.CloneTree(0)

    logging.info(f'Entries: {tree.GetEntries()}')
    for entryIdx in range(0, tree.GetEntries()):
        # Print status
        check_interval = 100000
        if entryIdx % check_interval == 0 and entryIdx > 0:
            logging.info(f'Processed {entryIdx:,d} events, kept {new_tree.GetEntries()}...')
        
        tree.GetEntry(entryIdx)
        
        prt_pid = getattr(tree, 'prt_pid')  # Reconstructed particle pids
        mc_pid = getattr(tree, 'mc_pid')  # MC-matched daughter pids
        px = getattr(tree, 'mc_px')  # MC-matched daughter px
        py = getattr(tree, 'mc_py')  # MC-matched daughter py
        pz = getattr(tree, 'mc_pz')  # MC-matched daughter pz 

        prt_pid = [int(pid) for pid in prt_pid]
        mc_pid = [int(pid) for pid in mc_pid]
        
        passed = True
        # This is not the most elegant way to do this but is straightforward.
        # eta -> mu+ mu- gamma case
        if gamma_flag:
            for i in range(0, len(mc_pid) - 3, 4):
                pids = mc_pid[i:i + 4]
                px4 = px[i:i + 4]
                py4 = py[i:i + 4]
                pz4 = pz[i:i + 4]
                if len(pids) < 4: continue

                if pids[0] != 221 or pids[1] != -13 or pids[2] != 13 or pids[3] != 22:
                    continue

                passed = all(passes_reqs(pids[j], px4[j], py4[j], pz4[j]) 
                            for j in range(4))
        # eta -> mu+ mu- case
        else:
            for i in range(0, len(mc_pid) - 2, 3):
                pids = mc_pid[i:i + 3]
                px3 = px[i:i + 3]
                py3 = py[i:i + 3]
                pz3 = pz[i:i + 3]
                if len(pids) < 3: continue

                if pids[0] != 221 or pids[1] != -13 or pids[2] != 13:
                    continue

                passed = all(passes_reqs(pids[j], px3[j], py3[j], pz3[j]) 
                            for j in range(3))
            
        if passed: new_tree.Fill()

    return new_tree

#-------------------------------------------------------------------------------

# Set up logging to stderr only
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
# Also print to stderr
logging.getLogger().addHandler(logging.StreamHandler())

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Input ROOT file")
# eta -> mu+ mu- (gamma) flag
parser.add_argument(
    '-g', '--gamma',
    action='store_true',
    help='eta -> mu+ mu- gamma flag'
)
args = parser.parse_args()

infile = args.infile
# Put outfile in same directory as infile with fixed name
if infile is not None:
    outfile = '/'.join(infile.split('/')[:-1]) + '/fiducial_requirements.root'
else: raise ValueError('Input file must be specified with -i or --infile flag.')
# Output file name
logging.info(f'Reading from {infile}, writing to {outfile}.')

tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

new_tfile = ROOT.TFile.Open(outfile, "RECREATE")
new_tfile.cd()

# Apply fiducial requirements
new_tree = apply_fiducial_reqs(tree, args.gamma)

logging.info(f'Total kept entries: {new_tree.GetEntries()}')

# Close input file
tfile.Close()

# Write new tree to output file 
new_tree.Write()
new_tfile.Close()

logging.info(f'Done: wrote reduced tree with fiducial requirements to {outfile}.')
