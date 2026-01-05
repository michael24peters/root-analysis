###############################################################################
# Script to analyze background.                                               #
# Author: Michael Peters                                                      #
###############################################################################

import ROOT
import sys
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true', 
                    help='Enable verbose output')
parser.add_argument('-s', '--sig', action='store_true',
                    help='Analyze signal file instead of minbias file')
args = parser.parse_args()

verbose = args.verbose
sig_file = args.sig

if sig_file:
    infile = 'MC_2018_Signal/eta2MuMuGamma_mc_20251208.root'
    outfile = 'hist/sig_hist_m.root'
else:
    # infile = 'red/reduced.root'
    infile = 'red/reduced_fiducial_cuts.root'
    # outfile = 'hist/TODO.root'

print(f'Reading from {infile}.')

# Counters
ncan = 0  # Total number of candidates
nsig = 0  # Total number of signal particles
nbkg = 0  # Total number of misreconstructed particles
npho = 0  # Number of misreconstructed photons
nmu = 0  # Number of misreconstructed muons
ndimu = 0  # Number of misreconstructed dimuons

bkgd = []  # List of misreconstructed particles with different PIDs
decay_bkgd = []  # List of reconstructed particles and their MC-matched particles per decay mode
fid_fail = []  # Particles failing LHCb fiducial cuts

# Combine files to create single histogram
tfile = ROOT.TFile.Open(infile, 'READ')
tree = tfile.Get('tree')

# Event loop
for entryIdx in range(0, tree.GetEntries()):
    
    tree.GetEntry(entryIdx)
    
    # Reconstructed particle information
    tag_pid = getattr(tree, 'tag_pid')
    prt_pid = getattr(tree, 'prt_pid')
    # MC-matching index information
    prt_idx_gen = getattr(tree, 'prt_idx_gen')
    prt_idx_mom = getattr(tree, 'prt_idx_mom')
    # Generator particle information
    mc_pid = getattr(tree, 'mc_pid')
    mc_idx_mom = getattr(tree, 'mc_idx_mom')
    # Passed LHCb fiducial cuts
    mc_req_pass = getattr(tree, 'mc_req_pass')

    # Reformat above lists for easier handling
    prt_pid = [int(pid) for pid in prt_pid]
    prt_idx_gen = [int(idx) for idx in prt_idx_gen]
    prt_idx_mom = [int(idx) for idx in prt_idx_mom]
    mc_pid = [int(pid) for pid in mc_pid]
    mc_idx_mom = [int(idx) for idx in mc_idx_mom]
    mc_req_pass = [int(val) for val in mc_req_pass]

    # Skip empty events
    ntags = len(tag_pid)
    if ntags == 0: continue
    ncan += ntags

    # Loop over all candidates in event. Use indexing information to find
    # correctly MC-matched events vs incorrectly MC-matched events.
    for i in range(ntags):
        dimu_err = [False, False]
        is_signal = True
        mcds = []  # MC-matched daughters collected for this candidate
        # Display event number at beginning of mcds list
        mcds.append(entryIdx)
        # Check that there is a MC eta candidate. If not, we immediately know this
        # is background.
         # TODO: This doesn't handle multiple eta candidates (very rare).
        if mc_pid[0] != 221: is_signal = False
        if tag_pid[i] != 221: continue  # skip failed reco/non-eta candidates
        # Loop for each daughters of candidate
        for j in range(i*3, i*3+3):
            if prt_idx_mom[j] != i:  break  # Skip failed reco, shouldn't happen
            
            # Get mc_pid of matched generator particle
            mcpid = mc_pid[prt_idx_gen[j]]
            
            # Check background conditions
            if prt_idx_gen[j] == -1: is_signal = False  # no MC match
            elif mc_pid[prt_idx_gen[j]] != prt_pid[j]: is_signal = False  # wrong MC match
            elif mc_idx_mom[prt_idx_gen[j]] != i: is_signal = False  # wrong mother match

            # Classify background type
            if not is_signal: 
                mcds.append([mcpid, prt_pid[j], True])
                if mcpid == -13: nmu += 1; dimu_err[0] = True
                elif mcpid == 13: nmu += 1; dimu_err[1] = True
                elif mcpid == 22: npho += 1
                else: bkgd.append([entryIdx, prt_pid[j], prt_idx_gen[j], mcpid, mc_idx_mom[prt_idx_gen[j]]])
                if all(dimu_err): ndimu += 1; dimu_err = [False, False]
                if mc_req_pass[prt_idx_gen[j]] == 0: fid_fail.append([entryIdx, prt_pid[j], prt_idx_gen[j], mcpid, mc_idx_mom[prt_idx_gen[j]]])
            else: mcds.append([mcpid, prt_pid[j], False])
    
        # Fill arrays and increment counters
        # Append flattened elements of tag_m to arrays
        if is_signal: nsig += 1
        else: decay_bkgd.append(mcds); nbkg += 1
        # arr_tot.append(tag_m[i])
        # if is_signal: arr_sig.append(tag_m[i]); nsig += 1
        # else: arr_bkg.append(tag_m[i]); nbkg += 1

# Extract set of unique misreconstructed MC particle pids
mcps = list(set([mcp[3] for mcp in bkgd]))

# Close TFile
tfile.Close()

# Print results
with open ('out/decay_bkgd.txt', 'w') as f:
     f.write('')  # Clear file contents
     output = '\n========================= Background Analysis Results =========================='
     output += f'\nTotal candidates processed: {ncan}'
     output += f'\nTotal signal candidates: {nsig}'
     output += f'\nTotal misreconstructed particles: {nbkg}'
     output += f'\n  - Misreconstructed photons: {npho}'
     output += f'\n  - Misreconstructed muons: {nmu}'
     output += f'\n  - Misreconstructed dimuons: {ndimu}'
     output += f'\nSet of unique misreconstructed MC-matched pids:\n  {mcps}'
     if verbose:
         if decay_bkgd:
             output += f'\nDecay modes of misreconstructed particles:'
             output += f'\n  Format:\n  Event #: [mcpid0, pid0, flag], ' + \
                       '[mcpid1, pid1, flag], [mcpid2, pid2, flag],\n  ' + \
                       'where !! is displayed if the two particles come from ' + \
                       'different gen-level\n  decays (non-signal). So, for example, ' + \
                       '[-13, -13, !!] might occur even\n  though they have the same ' + \
                       'pid but come from different decays.'
             for dec in decay_bkgd:
                 evtIdx = dec[0]
                 parts = []
                 for part in dec[1:]:
                     mcpid, pid, flag = part
                     if flag: parts.append(f'[{mcpid:5d}, {pid:4d}, !!]')
                     else: parts.append(f'[{mcpid:5d}, {pid:4d},   ]')
                 joined = ', '.join(parts)
                 output += f'\n  - Event {evtIdx:5d}: {joined}'
         # Probably redundant with above printout
         # if bkgd:
         #     output += f'\nOther misreconstructed particle pids:'
         #     for evt, pid, idx, mc_pid, mom in bkgd:
         #         output += f'\n  - Event {evt:5d}, pid {pid:4d}, idx {idx:2d}, mc_pid: {mc_pid:5d} mc_idx_mom {mom:2d}'
         if fid_fail: 
             output += f'\nParticles failing LHCb fiducial cuts but still MC-matched:'
             for evt, pid, idx, mc_pid, mom in fid_fail:
                 output += f'\n  - Event {evt:5d}, pid {pid:4d}, idx {idx:2d}, mc_pid: {mc_pid:5d} mc_idx_mom {mom:2d}'
     output += '\n===============================================================================\n'
     print(output)
     f.write(output)
