################################################################################
# Methods to calculate efficiencies.                                           #
# Author: Michael Peters                                                       #
################################################################################

import ROOT

def calc_ratio(tree):
    """Calculate efficiency as ratio with fiducial requirements in place."""
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

    return (nreco, ngen)


#===============================================================================

def calc_efficiency(tree):
    """Calculate efficiency with fiducial requirements in place."""
    nreco, ngen = calc_ratio(tree)
    if ngen == 0: return 0.0
    return nreco / ngen

#===============================================================================


def calc_sig_ratio(tree):
    """Calculate signal efficiency as ratio with fiducial requirements in place."""
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

    return (nreco_matches, ngen)

#===============================================================================

def calc_sig_efficiency(tree):
    """Calculate signal efficiency with fiducial requirements in place."""
    nreco_matches, ngen = calc_sig_ratio(tree)
    if ngen == 0: return 0.0
    return nreco_matches / ngen
