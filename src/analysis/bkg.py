# ╔════════════════════════════════════════════════════════════════════════════╗
# ║ Script to analyze background.                                              ║
# ║ Author: Michael Peters                                                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝

from __future__ import annotations

import ROOT
import argparse
from dataclasses import dataclass
from enum import Enum
from collections import Counter
import os
import sys
import json

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', required=True,
                    help = 'Input ROOT file')
parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='Enable verbose output')
parser.add_argument('-o', '--outfile', default='out/bkg_ana',
                    help='Output text file')
parser.add_argument('--pnnmu_cut', type=float, default=None,
                    help='PROBNNmu cut value applied to input file')
# Store number of daughters per candidate, 2 for eta -> mu+ mu-, 3 for 
# eta -> mu+ mu- gamma, etc.
parser.add_argument('--dtrs', action='store', default=2, type=int,
                    help='Number of daughters per candidate (default: 2)')
args = parser.parse_args()

# Handle arguments, set up IO.
infile = args.infile
os.makedirs('out', exist_ok=True)
outfile = args.outfile
verbose = args.verbose
pnnmu_cut = args.pnnmu_cut
if pnnmu_cut is not None:
    print(f'[info] Applying PROBNNmu cut: {pnnmu_cut}')

print(f'[info] Reading from {infile}, writing to {outfile}.')

# Possible error categories for a decay candidate
class ErrorType(str, Enum):
    MUP_PID_MISMATCH = 'MUP_PID_MISMATCH'
    MUP_ONLY_PID_MISMATCH = 'MUP_ONLY_PID_MISMATCH'
    MUM_PID_MISMATCH = 'MUM_PID_MISMATCH'
    MUM_ONLY_PID_MISMATCH = 'MUM_ONLY_PID_MISMATCH'
    DIMUON_PID_MISMATCH = 'DIMUON_PID_MISMATCH'
    PHOTON_PID_MISMATCH = 'PHOTON_PID_MISMATCH'
    MUP_ERROR = 'MUP_ERROR'
    MUP_ONLY_ERROR = 'MUP_ONLY_ERROR'
    MUM_ERROR = 'MUM_ERROR'
    MUM_ONLY_ERROR = 'MUM_ONLY_ERROR'
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
    
    # Print status every 500,000 events
    check_interval = 500000
    if entryIdx % check_interval == 0 and entryIdx > 0:
        print(f'[info] Processed {entryIdx:,d} events...')
    
    # Explicitly define types for easier handling
    # Reconstructed particle information
    tag_pid = [int(pid) for pid in tree.tag_pid]
    prt_pid = [int(pid) for pid in tree.prt_pid]
    # MC-matching index information
    prt_idx_gen = [int(idx) for idx in tree.prt_idx_gen]
    prt_idx_mom = [int(idx) for idx in tree.prt_idx_mom]
    # Generator particle information
    mc_pid = [int(pid) for pid in tree.mc_pid]
    mc_idx_mom = [int(idx) for idx in tree.mc_idx_mom]
    # PROBNNmu branch
    probnn_mu = [float(prob) for prob in tree.prt_pnn_mu]

    # Skip empty events
    ntags = len(tag_pid)
    if ntags == 0: continue
    
    error_log = ''
    for i in range(ntags):
        if tag_pid[i] != 221: continue  # skip failed reco/non-eta candidates

        # NOTE: can only handle one gen-level eta candidate per event, improve
        # later if needed, or just select the first one and ignore the others.
        # If so, then this check is needed to prevent out-of-bounds errors.
        #
        # Define j0 and j1 depending on gamma_flag
        ndtrs = args.dtrs  # Number of daughters per candidate
        j0, j1 = i * ndtrs, i * ndtrs + ndtrs
        if j1 > len(prt_pid) or j1 > len(prt_idx_gen) or j1 > len(prt_idx_mom):
            error_log += f'[warning] Event {entryIdx}, candidate {i} skipped, has incomplete daughter info.\n'
            continue

        # Apply PROBNNmu cut if specified
        if pnnmu_cut is not None:
            passed = True
            for j in range(j0, j1):
                if abs(prt_pid[j]) == 13 and probnn_mu[j] < pnnmu_cut:
                    passed = False
                    break
            # If candidate fails probnn_mu cut, skip it.
            if not passed: continue
        ncan += 1  # Count candidate only after passing pnnmu cut

        is_signal = True
        dtrs: list[DaughterMatch] = []
        dimu_mismatch = [False, False]
        dimu_err = [False, False]

        # Requires at least one MC eta candidate
        # NOTE: cannot handle multiple eta candidates (very rare)
        if mc_pid[0] != 221: is_signal = False

        for j in range(j0, j1):
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
                print(f'[warning] Could not assign mc_pid or mc_idx_mom for daughter.')

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
mup_err_only_count, mum_err_only_count = 0, 0
mu_AND_dimu_err_count = 0
for can in candidates:
    # Dimuon errors can only be observed at candidate-level
    if can.has_dimu_mismatch: 
        err_counters[ErrorType.DIMUON_PID_MISMATCH] += 1
    if can.has_dimu_err: 
        err_counters[ErrorType.DIMUON_ERROR] += 1
    for dtr in can.dtrs:
        # Increment daughter counters
        for err in ERROR_TYPES:
            if dtr.err_type == err:
                err_counters[err] += 1
                if err == ErrorType.MUP_PID_MISMATCH and not can.has_dimu_mismatch:
                    mup_err_only_count += 1
                elif err == ErrorType.MUM_PID_MISMATCH and not can.has_dimu_mismatch:
                    mum_err_only_count += 1
                elif err == ErrorType.MUP_ERROR and not can.has_dimu_err:
                    mu_AND_dimu_err_count += 1
                elif err == ErrorType.MUM_ERROR and not can.has_dimu_err:
                    mu_AND_dimu_err_count += 1

# Add MUON_ONLY_* counters
err_counters[ErrorType.MUP_ONLY_PID_MISMATCH] = mup_err_only_count
err_counters[ErrorType.MUM_ONLY_PID_MISMATCH] = mum_err_only_count
err_counters[ErrorType.MUP_ONLY_ERROR] = mup_err_only_count
err_counters[ErrorType.MUM_ONLY_ERROR] = mum_err_only_count


#───────────────────────────────────────────────────────────────────────────────
def format_pid_freq_table(label, rows: list[tuple[int, int]]) -> str:
    """
    Return a plain text table for (pid, count) rows.
    """
    if not rows:
        return '(none)\n'

    out = ''
    # Determine column widths based on label and data
    label_w = len(label)
    pid_w = max(len(" PID "), max(len(str(pid)) for pid, _ in rows))
    cnt_w = max(len(" # "), max(len(str(cnt)) for _, cnt in rows))
    
    # Label should look like '\n───── MU- ─────\n' but with dynamic width based
    # on label
    # Label
    out += f'\n{"═" * ((pid_w + cnt_w + 5 - label_w) // 2)} {label} ' + \
           f'{"═" * ((pid_w + cnt_w + 5 - label_w) // 2)}\n'
    # Header
    out += f"┌─{' PID ':>{pid_w}}─┬─{' # ':<{cnt_w}}─┐\n"
    out += f"├─{'─' * pid_w}─┼─{'─' * cnt_w}─┤\n"
    # Data
    for pid, count in rows:
        out += f'│ {pid:>{pid_w}d} │ {count:<{cnt_w}d} │\n'
    out += '└' + '─' * (pid_w + 2) + '┴' + '─' * (cnt_w + 2) + '┘\n'
    return out


#───────────────────────────────────────────────────────────────────────────────
def get_analytics(verbose=False):
    """
    Return a formatted string of background analysis results, including counts,
    error rates, and list of PID mismatches ranked by frequency for each
    particle type. If verbose=True, also include candidate-level information and
    error log. Designed for printing to console or writing to text file.
    """
    output = ''
    W = 48  # Line width
    output += '═' * W + '\n'
    output += f"{'Background analysis results':^{W}}\n"
    output += '═' * W + '\n'
    # Key to explain counters
    output += '  *_MISMATCH: Daughter has MC match but reco' + '\n' + \
              '    pid does not match gen pid.\n'
    output += '  *_ERROR: Daughter has MC match but reco did' + '\n' + \
              '    not match to candidate gen dtr.\n'
    output += '─' * W + '\n'
    # Summary statistics
    output += f'  {"Total candidates processed:":<31} {ncan:5d}\n'
    output += f'  {"Total signal candidates:":<31} {nsig:5d}\n'
    output += f'  {"Total background candidates:":<31} {nbkg:5d}\n'
    # Error counts
    output += '─' * W + '\n'
    output += f"{'Background error counts':^{W}}\n"
    output += '─' * W + '\n'
    for err in err_counters:
        # Remove ErrorType. prefix for display
        label = err.name + ':'
        output += f'  {label:<31} {err_counters[err]:5d}\n'
    # Error rates
    output += '─' * W + '\n'
    output += f"{'Background error rates':^{W}}\n"
    output += '─' * W + '\n'
    # P(dimuon error | muon error) = N(muon & dimuon) / N(muon)
    prob_dimu_given_mu = err_counters[ErrorType.DIMUON_ERROR] / \
                         max(err_counters[ErrorType.MUP_ERROR], 
                         err_counters[ErrorType.MUM_ERROR], 1)
    # P(muon error)
    prob_mu = (err_counters[ErrorType.MUP_ERROR] / ncan,
               err_counters[ErrorType.MUM_ERROR] / ncan)
    # P(photon error)
    prob_pho = err_counters[ErrorType.PHOTON_ERROR] / ncan
    output += f'  P(dimuon error | muon error) = {prob_dimu_given_mu:.4f}\n'
    output += f'  P(mu+ error)                 = {prob_mu[0]:.4f}\n'
    output += f'  P(mu- error)                 = {prob_mu[1]:.4f}\n'
    output += f'  P(photon error)              = {prob_pho:.4f}\n'
    # List of MC pids causing mismatches
    output += '─' * W + '\n'
    output += f"{'List of PID mismatches (ranked by frequency)':^{W}}\n"
    output += '─' * W + '\n'
    # Sort list by frequency
    c_mup = Counter(mup_mismatches)
    c_mum = Counter(mum_mismatches)
    c_pho = Counter(pho_mismatches)
    c_other = Counter(other_mismatches)
    # Each Counter.most_common() returns [(pid, count), ...]
    # Formatted tables for each particle type
    output += format_pid_freq_table('MU+', c_mup.most_common())

    output += format_pid_freq_table('MU-', c_mum.most_common())

    output += format_pid_freq_table('PHOTON', c_pho.most_common())

    output += format_pid_freq_table('OTHER', c_other.most_common())
    output += '─' * W + '\n'

    # Verbose output
    W = 120  # Wider line width for verbose output
    verbose_output = ''
    if verbose:
        verbose_output += '═' * W + '\n'
        verbose_output += f"{'Verbose candidate information':^{W}}\n"
        verbose_output += '═' * W + '\n'
        for can in candidates:
            verbose_output += f'\nEvent {can.evt}, Candidate {can.can_idx}:\n'
            for dtr in can.dtrs:
                verbose_output += f'  Daughter PID {dtr.prt_pid:3d}, '
                verbose_output += f'Gen idx {dtr.prt_idx_gen:2d}, '
                verbose_output += f'MC PID {dtr.mc_pid:5d}, '
                verbose_output += f'MC mom idx {dtr.mc_idx_mom:2d}, '
                verbose_output += f'Error type: {dtr.err_type}\n'
        verbose_output += '─' * W + '\n'

    return output, verbose_output


#───────────────────────────────────────────────────────────────────────────────
def get_pid_mismatch_analytics(rows: list[tuple[int, int]]):
    """
    Get list of PID mismatches ranked by frequency for a given particle type and
    store to a dict for JSON output.
    """
    sub_analytics = {}
    for pid, count in rows: sub_analytics[pid] = count
    return sub_analytics

#───────────────────────────────────────────────────────────────────────────────
def get_json_analytics(verbose=False):
    """
    Return analytics in a JSON-serializable dict format.
    """
    # Basic counts
    analytics = {
        'total_candidates': ncan,
        'total_signal': nsig,
        'total_background': nbkg
    }
    # Background error counts
    for err in err_counters:
        # Remove ErrorType. prefix for display
        label = str(err).strip('ErrorType.') + ':'
        counter = err_counters[err]
        analytics[label] = counter
    # Background error rates
    # P(dimuon error | muon error) = N(muon & dimuon) / N(muon)
    prob_dimu_given_mu = err_counters[ErrorType.DIMUON_ERROR] / \
                         max(err_counters[ErrorType.MUP_ERROR], 
                         err_counters[ErrorType.MUM_ERROR], 1)
    # P(muon error)
    prob_mu = (err_counters[ErrorType.MUP_ERROR] / ncan,
               err_counters[ErrorType.MUM_ERROR] / ncan)
    # P(photon error)
    prob_pho = err_counters[ErrorType.PHOTON_ERROR] / ncan
    analytics['error_rates'] = {
        'P(dimuon error | muon error)': prob_dimu_given_mu,
        'P(mu+ error)': prob_mu[0],
        'P(mu- error)': prob_mu[1],
        'P(photon error)': prob_pho
    }
    # List of MC pids causing mismatches
    c_mup = Counter(mup_mismatches)
    c_mum = Counter(mum_mismatches)
    c_pho = Counter(pho_mismatches)
    c_other = Counter(other_mismatches)
    analytics['pid_mismatches'] = {
        'MU+': get_pid_mismatch_analytics(c_mup.most_common()) if c_mup else None,
        'MU-': get_pid_mismatch_analytics(c_mum.most_common()) if c_mum else None,
        'PHOTON': get_pid_mismatch_analytics(c_pho.most_common()) if c_pho else None,
        'OTHER': get_pid_mismatch_analytics(c_other.most_common()) if c_other else None
    }

    # Verbose candidate-level information
    verbose_analytics = None
    if verbose:
        verbose_analytics = {}
        verbose_analytics['candidates'] = []
        for can in candidates:
            can_info = {
                'event': can.evt,
                'candidate_index': can.can_idx,
                'has_dimu_mismatch': can.has_dimu_mismatch,
                'has_dimu_err': can.has_dimu_err,
                'daughters': []
            }
            for dtr in can.dtrs:
                dtr_info = {
                    'prt_pid': dtr.prt_pid,
                    'prt_idx_gen': dtr.prt_idx_gen,
                    'mc_pid': dtr.mc_pid,
                    'mc_idx_mom': dtr.mc_idx_mom,
                    'err_type': str(dtr.err_type)
                }
                can_info['daughters'].append(dtr_info)
            verbose_analytics['candidates'].append(can_info)

    return analytics, verbose_analytics


#───────────────────────────────────────────────────────────────────────────────
# Print or write analytics to file
output, verbose_output = get_analytics(verbose=verbose)
analytics, verbose_analytics = get_json_analytics(verbose=verbose)
print(output)
# Write text output to file
with open(outfile + '.txt', 'w') as f:
    f.write('') # Clear file contents
    f.write(output)
    if verbose: f.write('\n' + verbose_output)
with open(outfile + '.json', 'w') as f:
    json.dump(analytics, f, indent=2)
    if verbose: json.dump(verbose_analytics, f, indent=2)
print(f'[done] Background analysis results saved to: {outfile}')
