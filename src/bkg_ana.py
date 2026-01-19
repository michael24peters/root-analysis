################################################################################
# Script to analyze background.                                                #
# Author: Michael Peters                                                       #
################################################################################

from __future__ import annotations

import ROOT
import argparse
from dataclasses import dataclass
from enum import Enum
from collections import Counter
from utils.calculate_efficiency import calc_ratio, calc_sig_ratio
from utils.calculate_efficiency import calc_efficiency, calc_sig_efficiency

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true', 
                    help='Enable verbose output')
parser.add_argument('-s', '--sig', action='store_true',
                    help='Analyze signal file instead of minbias file')
parser.add_argument('-o', '--outfile', action='store_true',
                    help='Write to output text file (default: none)')
args = parser.parse_args()

verbose = args.verbose
is_sig_file = args.sig
write_to_outfile = args.outfile

if is_sig_file:
    infile = 'MC_2018_Signal/eta2MuMuGamma_mc_20251208.root'
else:
    # infile = 'red/reduced.root'
    # infile = 'red/reduced_fiducial_cuts.root'
    infile = 'red/reduced_fiducial_reqs.root'

if write_to_outfile: print(f'Reading from {infile}, writing to out/bkg_ana.txt.')
else: print(f'Reading from {infile}.')

# Possible error categories for a decay candidate
class ErrorType(str, Enum):
    MUP_PID_MISMATCH = 'MUP_PID_MISMATCH'
    MUM_PID_MISMATCH = 'MUM_PID_MISMATCH'
    DIMUON_PID_MISMATCH = 'DIMUON_PID_MISMATCH'
    PHOTON_PID_MISMATCH = 'PHOTON_PID_MISMATCH'
    MUP_ERROR = 'MUP_ERROR'
    MUM_ERROR = 'MUM_ERROR'
    DIMUON_ERROR = 'DIMUON_ERROR'
    PHOTON_ERROR = 'PHOTON_ERROR'
    OTHER_ERROR = 'OTHER_ERROR'

@dataclass
class DaughterMatch:
    prt_pid: int
    prt_idx_gen: int
    mc_pid: int | None
    mc_idx_mom: int | None
    err_type: ErrorType

@dataclass
class Candidate:
    evt: int
    can_idx: int
    dtrs: list[DaughterMatch]
    has_dimu_mismatch: bool
    has_dimu_err: bool

fid_fail = []  # Particles failing LHCb fiducial cuts
ncan, nsig, nbkg = 0, 0, 0  # Total candidates, signal, and background counters
candidates: list[Candidate] = []  # List of all candidates
mup_mismatches = []  # List of MC pids causing mu+ PID mismatches
mum_mismatches = []  # List of MC pids causing mu- PID mismatches
pho_mismatches = []  # List of MC pids causing photon PID mismatches
other_mismatches = []  # List of MC pids causing other mismatches

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

    # Reformat above lists for easier handling
    prt_pid = [int(pid) for pid in prt_pid]
    prt_idx_gen = [int(idx) for idx in prt_idx_gen]
    prt_idx_mom = [int(idx) for idx in prt_idx_mom]
    mc_pid = [int(pid) for pid in mc_pid]
    mc_idx_mom = [int(idx) for idx in mc_idx_mom]

    # Skip empty events
    ntags = len(tag_pid)
    if ntags == 0: continue
    ncan += ntags
    
    for i in range(ntags):
        if tag_pid[i] != 221: continue  # skip failed reco/non-eta candidates

        is_signal = True
        dtrs: list[DaughterMatch] = []
        dimu_mismatch = [False, False]
        dimu_err = [False, False]

        # Requires at least one MC eta candidate
        # Note: cannot handle multiple eta candidates (very rare)
        if mc_pid[0] != 221: is_signal = False

        for j in range(i*3, i*3+3):
            if prt_idx_mom[j] != i: break  # Skip failed reco, shouldn't happen
            if tag_pid[i] != 221: continue  # Skip failed reco/non-eta candidates

            is_pid_mismatch, is_from_eta = False, True

            # Particle has no MC match
            if prt_idx_gen[j] == -1:
                is_signal = False
                is_from_eta = False
                err_type = ErrorType.OTHER_ERROR
            # Particle matches to a MC particle with a different pid
            elif mc_pid[prt_idx_gen[j]] != prt_pid[j]:
                is_signal = False
                is_from_eta = False
                is_pid_mismatch = True
            # Particle is correct pid but didn't come from eta candidate
            elif mc_idx_mom[prt_idx_gen[j]] != i:
                is_signal = False
                is_from_eta = False
            
            # Background
            if not is_from_eta:
                # Classify PID mismatch type
                if is_pid_mismatch:
                    # MUON_PID_MISMATCH condition
                    if prt_pid[j] == -13:
                        err_type = ErrorType.MUP_PID_MISMATCH
                        mup_mismatches.append(mc_pid[prt_idx_gen[j]])
                        dimu_mismatch[0] = True
                    # MUON_PID_MISMATCH condition
                    elif prt_pid[j] == 13:
                        err_type = ErrorType.MUM_PID_MISMATCH
                        mum_mismatches.append(mc_pid[prt_idx_gen[j]])
                        dimu_mismatch[1] = True
                    # PHOTON_PID_MISMATCH condition
                    elif prt_pid[j] == 22:
                        err_type = ErrorType.PHOTON_PID_MISMATCH
                        pho_mismatches.append(mc_pid[prt_idx_gen[j]])
                    else:
                        err_type = ErrorType.OTHER_ERROR
                        other_mismatches.append(mc_pid[prt_idx_gen[j]])
                # Classify daughter error type
                elif not is_from_eta:
                    if prt_pid[j] == -13:
                        err_type = ErrorType.MUP_ERROR
                        dimu_err[0] = True
                    elif prt_pid[j] == 13:
                        err_type = ErrorType.MUM_ERROR
                        dimu_err[1] = True
                    elif prt_pid[j] == 22:
                        err_type = ErrorType.PHOTON_ERROR
                # Some other error which should be studied
                else: err_type = ErrorType.OTHER_ERROR
            else: err_type = None  # Correctly matched dtr of signal candidate

            try: 
                dtrs.append(DaughterMatch(prt_pid=prt_pid[j],
                    prt_idx_gen=prt_idx_gen[j],
                    mc_pid=mc_pid[prt_idx_gen[j]],
                    mc_idx_mom=mc_idx_mom[prt_idx_gen[j]],
                    err_type=err_type))
            except: 
                dtrs.append(DaughterMatch(prt_pid=prt_pid[j],
                    prt_idx_gen=prt_idx_gen[j],
                    mc_pid=None,
                    mc_idx_mom=None,
                    err_type=err_type))
                print('Warning: Could not assign mc_pid or mc_idx_mom for daughter.')

        if is_signal: nsig += 1
        else: nbkg += 1
        candidate = Candidate(evt=entryIdx,
                              can_idx=i,
                              dtrs=dtrs,
                              has_dimu_mismatch=all(dimu_mismatch),
                              has_dimu_err=all(dimu_err))
        candidates.append(candidate)
        
# Collect analytics
ERROR_TYPES = list(ErrorType)
err_counters = {err: 0 for err in ERROR_TYPES}
for can in candidates:
    # Dimuon errors can only be observed at candidate-level
    if can.has_dimu_mismatch: err_counters[ErrorType.DIMUON_PID_MISMATCH] += 1
    if can.has_dimu_err: err_counters[ErrorType.DIMUON_ERROR] += 1
    for dtr in can.dtrs:
        # Increment daughter counters
        for err in ERROR_TYPES:
            if dtr.err_type == err:
                err_counters[err] += 1

#-------------------------------------------------------------------------------


def format_pid_freq_table(rows: list[tuple[int, int]]) -> str:
    """Return a plain text table for (pid, count) rows."""
    if not rows:
        return '(none)\n'

    pid_w = max(len("PID"), max(len(str(pid)) for pid, _ in rows))
    cnt_w = max(len("#"), max(len(str(cnt)) for _, cnt in rows))

    out = ''
    out += '+' + '-' * (pid_w + 2) + '+' + '-' * (cnt_w + 2) + '+\n'
    out += f'| {"PID":>{pid_w}} | {"#":<{cnt_w}} |\n'
    out += f'| {"-" * pid_w} | {"-" * cnt_w} |\n'
    for pid, count in rows:
        out += f'| {pid:>{pid_w}d} | {count:<{cnt_w}d} |\n'
    out += '+' + '-' * (pid_w + 2) + '+' + '-' * (cnt_w + 2) + '+\n'
    return out


#-------------------------------------------------------------------------------


def get_analytics():
    output = '='*25 + ' Background Analysis Results ' + '='*26 + '\n'
    # Key to explain counters
    output += '*_MISMATCH: Daughter has MC match but reco pid does not match gen pid.\n'
    output += '*_ERROR: Daughter has MC match but did reco did not match to candidate gen dtr.\n'
    output += 'Note: DIMUON_* errors do not overwrite single MU*_* errors.\n'
    output += '-'*80 + '\n'
    # Summary statistics
    output += f'Total candidates processed: {ncan}\n'
    output += f'Total signal candidates: {nsig}\n'
    output += f'Total background candidates: {nbkg}\n'
    output += '-'*80 + '\n'
    # Error counts
    output += 'Background error counts:\n'
    for err in err_counters:
        # Remove ErrorType. prefix for display
        output += f'- {err.strip("ErrorType.")}: {err_counters[err]}\n'
    output += '-'*80 + '\n'
    # List of MC pids causing mismatches
    output += 'List of PID mismatches (ranked by frequency):\n'
    # Sort list by frequency
    c_mup = Counter(mup_mismatches)
    c_mum = Counter(mum_mismatches)
    c_pho = Counter(pho_mismatches)
    c_other = Counter(other_mismatches)
    # Each Counter.most_common() returns [(pid, count), ...]
    # Formatted tables
    output += '\n--- MU+ ---\n'
    output += format_pid_freq_table(c_mup.most_common())

    output += '\n--- MU- ---\n'
    output += format_pid_freq_table(c_mum.most_common())

    output += '\n--- PHOTON ---\n'
    output += format_pid_freq_table(c_pho.most_common())

    output += '\n--- OTHER ---\n'
    output += format_pid_freq_table(c_other.most_common())
    output += '-'*80 + '\n'

    # Calculate efficiencies with fiducial requirements in place
    eff_ratio = calc_ratio(tree)
    sig_eff_ratio = calc_sig_ratio(tree)
    eff = calc_efficiency(tree)
    sig_eff = calc_sig_efficiency(tree)

    output += f'Efficiency with fiducial requirements: {eff_ratio[0]}/{eff_ratio[1]} = {eff:.4f}\n'
    output += f'Signal efficiency with fiducial requirements: {sig_eff_ratio[0]}/{sig_eff_ratio[1]} = {sig_eff:.4f}\n'
    output += '-'*80 + '\n'
    
    # Verbose output
    verbose_output = ''
    if verbose:
        verbose_output += 'Verbose candidate information:\n'
        for can in candidates:
            verbose_output += f'\nEvent {can.evt}, Candidate {can.can_idx}:\n'
            for dtr in can.dtrs:
                verbose_output += f'  Daughter PID {dtr.prt_pid:3d}, '
                verbose_output += f'Gen idx {dtr.prt_idx_gen:2d}, '
                verbose_output += f'MC PID {dtr.mc_pid:5d}, '
                verbose_output += f'MC mom idx {dtr.mc_idx_mom:2d}, '
                verbose_output += f'Error type: {dtr.err_type}\n'

    return output, verbose_output


#-------------------------------------------------------------------------------

# Print or write analytics to file
output, verbose_output = get_analytics()
print(output)
if write_to_outfile:
    with open('out/bkg_ana.txt', 'w') as f:
        f.write('') # Clear file contents
        f.write(output)
        if verbose:
            f.write('\n' + verbose_output)
    print('Background analysis results written to out/bkg_ana.txt file.')
elif verbose: 
    print('Verbose output not written. Use -o flag to write to file.')
    print('-' * 80)
