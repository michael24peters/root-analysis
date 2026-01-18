################################################################################
# Methods to apply fiducial requirements and calculate efficiencies.           #
# Author: Michael Peters                                                       #
################################################################################
# TODO: Split into two files
# TODO: Might replace offline_gen_cuts.py entirely with this

import ROOT
import array
import sys
import argparse

#===============================================================================


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


#===============================================================================


def passes_reqs(pid, px, py, pz):
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


#===============================================================================


def apply_fiducial_reqs(tree):
    """Apply fiducial cuts to generator-level particles.
    
    Returns a new tree with only events that pass the fiducial cuts.
    Fiducial requirements:
    - Pseudorapidity (eta) in [2, 4.5]
    - Muon pT > 500 MeV
    - Muon P > 3 GeV
    - Photon pT > 500 MeV
    """

    new_tree = tree.CloneTree(0)

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

        prt_pid = [int(pid) for pid in prt_pid]
        mc_pid = [int(pid) for pid in mc_pid]
        
        passed = True
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
            
        if passed: new_tree.Fill()

    return new_tree


#===============================================================================


def calc_efficiency(tree):
    """Calculate efficiency with fiducial requirements in place."""
    nreco, ngen = 0, 0
    # Loop through tree, count number of reconstructed candidates and number of
    # generator level candidates.
    for entryIdx in range(0, tree.GetEntries()):
        tree.GetEntry(entryIdx)
        
        tag_pid = getattr(tree, 'tag_pid')
        mc_pid = getattr(tree, 'mc_pid')

        tag_pid = [int(pid) for pid in tag_pid]
        mc_pid = [int(pid) for pid in mc_pid]
        
        # Loop over mc_pid and count candidates
        nreco += 1 if len(tag_pid) > 0 else 0
        # Each generator level candidate has 4 particles, and all the junk has
        # already been thrown out, so we can just divide by 4.
        ngen += len(mc_pid) // 4

    return nreco / ngen if ngen > 0 else 0.0


#===============================================================================


def calc_sig_efficiency(tree):
    """Calculate signal efficiency with fiducial requirements in place."""
    nreco_matches, ngen = 0, 0
    # Loop through tree, count number of reconstructed decays which match to
    # generator level decays and count number of generator level decays.
    for entryIdx in range(0, tree.GetEntries()):
        tree.GetEntry(entryIdx)
        
        tag_pid = getattr(tree, 'tag_pid')
        prt_pid = getattr(tree, 'prt_pid')
        prt_idx_gen = getattr(tree, 'prt_idx_gen')
        mc_pid = getattr(tree, 'mc_pid')
        mc_idx_mom = getattr(tree, 'mc_idx_mom')

        tag_pid = [int(pid) for pid in tag_pid]
        prt_pid = [int(pid) for pid in prt_pid]
        prt_idx_gen = [int(idx) for idx in prt_idx_gen]
        mc_pid = [int(pid) for pid in mc_pid]
        mc_idx_mom = [int(idx) for idx in mc_idx_mom]
        
        # Loop through prt_pid and check if each prt_idx_gen points to a matching
        # mc_pid which is part of a signal decay. If all 3 reconstructed daughters
        # match to generator level daughters from the same signal decay, count
        # this as a reconstructed signal decay matching to generator level.
        ngen += len(mc_pid) // 4
        for i in range(0, len(prt_pid), 3):
            pids = prt_pid[i:i + 3]
            if len(pids) < 3: continue
            if pids[0] != -13 or pids[1] != 13 or pids[2] != 22: continue
            # This is a reco signal candidate, check if all daughters match to
            # generator level signal daughters.
            passed = True
            for j in range(1, 3):
                gen_idx = prt_idx_gen[i + j]
                if gen_idx < 0 or gen_idx >= len(mc_pid):
                    passed = False
                    break
                mcp = mc_pid[gen_idx]
                mcp_mom_idx = mc_idx_mom[gen_idx]
                if mcp_mom_idx == -1:
                    passed = False
                    break
                mc_mom_pid = mc_pid[mcp_mom_idx]
                if not (mc_mom_pid == 221 and 
                        ((pids[j] == -13 and mcp == -13) or
                         (pids[j] == 13 and mcp == 13) or
                         (pids[j] == 22 and mcp == 22))):
                    passed = False
                    break
            if passed: nreco_matches += 1

    return nreco_matches / ngen if ngen > 0 else 0.0


#===============================================================================

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
    def_outfile = 'root/sig_reduced_fiducial_reqs.root'
else:
    infile = 'red/reduced.root'
    def_outfile = 'red/reduced_fiducial_reqs.root'

outfile = ('red' + args.outfile) if args.outfile else def_outfile
print(f'Reading from {infile}, writing to {outfile}.')

tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

new_tfile = ROOT.TFile.Open(outfile, "RECREATE")
new_tfile.cd()

new_tree = apply_fiducial_reqs(tree)

print(f'Total kept entries: {new_tree.GetEntries()}')
eff = calc_efficiency(new_tree)
sig_eff = calc_sig_efficiency(new_tree)

# Close input file
tfile.Close()

# Write new tree to output file 
new_tree.Write()
new_tfile.Close()

print(f'Done: wrote reduced tree with fiducial requirements to {outfile}.')

print(f'Efficiency with fiducial requirements: {eff:.6f}')
print(f'Signal efficiency with fiducial requirements: {sig_eff:.6f}')