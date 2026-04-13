# ╔════════════════════════════════════════════════════════════════════════════╗
# ║ Script to plot eta mass distributions broken down by background category.  ║
# ║ Each candidate is assigned to exactly one category based on the            ║
# ║ combination of daughter error types from bkg.py classification.            ║
# ║ Author: Michael Peters                                                     ║
# ╚════════════════════════════════════════════════════════════════════════════╝

from __future__ import annotations

import ROOT
import argparse
import logging
import os
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

from anaroot.src.analysis.bkg import classify_candidate, ErrorType, Candidate

# ─── Category definitions ─────────────────────────────────────────────────────
# Ordered for legend display: signal first, then mismatch group, error group,
# then catch-all.
CATEGORIES = [
    'Signal',
    'Muon PID mismatch only',
    'Dimuon PID mismatch only',
    'Photon PID mismatch only',
    'Combination PID mismatch',
    'Muon error only',
    'Dimuon error only',
    'Photon error only',
    'Combination error',
    'Other error',
]


def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot eta mass distributions by background error category.')
    parser.add_argument('-i', '--infile', required=True,
                        help='Input ROOT file')
    parser.add_argument('-o', '--outfile', default='out/bkg_mass',
                        help='Output file prefix (default: out/bkg_mass)')
    parser.add_argument('--dtrs', type=int, default=2,
                        help='Number of daughters per candidate (default: 2)')
    parser.add_argument('--pnnmu_cut', type=float, default=None,
                        help='PROBNNmu cut applied to muon daughters')
    parser.add_argument('--log', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='Logging level (default: INFO)')
    return parser.parse_args()


# ─── Per-candidate label ──────────────────────────────────────────────────────
def get_candidate_label(can: Candidate) -> str:
    """
    Map a Candidate to exactly one category label based on the full set of
    daughter error types. Candidates with a mix of PID-mismatch and error
    daughters (or any OTHER_ERROR) fall into 'Other error'.
    """
    # Signal candidates
    if can.is_sig: return 'Signal'

    has_mup_mm = any(d.err_type == ErrorType.MUP_PID_MISMATCH for d in can.dtrs)
    has_mum_mm = any(d.err_type == ErrorType.MUM_PID_MISMATCH for d in can.dtrs)
    has_pho_mm = any(d.err_type == ErrorType.PHOTON_PID_MISMATCH for d in can.dtrs)
    has_mup_err = any(d.err_type == ErrorType.MUP_ERROR for d in can.dtrs)
    has_mum_err = any(d.err_type == ErrorType.MUM_ERROR for d in can.dtrs)
    has_pho_err = any(d.err_type == ErrorType.PHOTON_ERROR for d in can.dtrs)
    has_other = any(d.err_type == ErrorType.OTHER_ERROR for d in can.dtrs)

    has_any_mm = has_mup_mm or has_mum_mm or has_pho_mm
    has_any_err = has_mup_err or has_mum_err or has_pho_err

    # Catch-all for other errors and mixed categories not captured
    if has_other or (has_any_mm and has_any_err):
        return 'Other error'

    # PID mismatch errors
    if has_any_mm:
        if has_pho_mm and (has_mup_mm or has_mum_mm): return 'Combination PID mismatch'
        if has_mup_mm and has_mum_mm: return 'Dimuon PID mismatch only'
        if has_pho_mm: return 'Photon PID mismatch only'
        return 'Muon PID mismatch only'

    # Non-PID errors
    if has_any_err:
        if has_pho_err and (has_mup_err or has_mum_err): return 'Combination error'
        if has_mup_err and has_mum_err: return 'Dimuon error only'
        if has_pho_err: return 'Photon error only'
        return 'Muon error only'

    return 'Other error'


# ─── Event loop ───────────────────────────────────────────────────────────────
def run_event_loop(args) -> dict[str, list[float]]:
    """
    Loop over the input ROOT tree, classify each candidate, and collect
    tag_m values grouped by category label.
    """
    masses: dict[str, list[float]] = defaultdict(list)

    tfile = ROOT.TFile.Open(args.infile, 'READ')
    tree = tfile.Get('tree')

    for entryIdx in range(tree.GetEntries()):
        tree.GetEntry(entryIdx)

        if entryIdx % 500_000 == 0 and entryIdx > 0:
            logging.info(f'Processed {entryIdx:,d} events...')

        tag_pid = [int(pid) for pid in tree.tag_pid]
        tag_m = [float(m) for m in tree.tag_m]
        prt_pid = [int(pid) for pid in tree.prt_pid]
        prt_idx_gen = [int(idx) for idx in tree.prt_idx_gen]
        prt_idx_mom = [int(idx) for idx in tree.prt_idx_mom]
        mc_pid = [int(pid) for pid in tree.mc_pid]
        mc_idx_mom = [int(idx) for idx in tree.mc_idx_mom]
        probnn_mu = [float(p) for p in tree.prt_pnn_mu]

        # Skip events with no candidates
        if not tag_pid: continue

        for i in range(len(tag_pid)):
            if tag_pid[i] != 221:
                continue
            candidate = classify_candidate(entryIdx, i, prt_pid, prt_idx_gen, 
                prt_idx_mom, mc_pid, mc_idx_mom, probnn_mu, args)
            if candidate is None: continue
            label = get_candidate_label(candidate)
            masses[label].append(tag_m[i])

    tfile.Close()
    return masses


# ─── Plot ─────────────────────────────────────────────────────────────────────
def plot_masses(masses: dict[str, list[float]], outfile: str) -> None:
    """
    Overlay eta mass histograms for each category on a single matplotlib
    figure and save to {outfile}.png.
    """
    bins = np.arange(400, 710, 10)  # 400–700 MeV, 10 MeV bins
    colors = plt.get_cmap('tab10').colors  # 10 distinct colours

    # Font settings
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['font.size'] = 14
    # Creature plot
    plt.figure(figsize=(12, 9))

    for idx, label in enumerate(CATEGORIES):
        data = masses.get(label, [])
        if not data:
            continue
        arr = np.array(data)
        plt.hist(arr, bins=bins, histtype='step',
                color=colors[idx % len(colors)],
                label=f'{label} (n={len(arr):,})')

    plt.axvline(x=547.86, color='black', linestyle='--', linewidth=1,
               label='PDG $m_\\eta$ = 547.86 MeV')
    plt.xlabel('Mass [MeV]')
    plt.ylabel('Candidates')
    plt.title('Eta Mass Signal vs Background Categories')
    plt.grid(alpha=0.3)
    plt.legend(frameon=False, fontsize=11)

    os.makedirs(os.path.dirname(outfile) or '.', exist_ok=True)
    plt.savefig(outfile + '.png', dpi=300, bbox_inches='tight')
    logging.info(f'Saved plot to {outfile}.png')


# ─── Main ─────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log),
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    logging.info(f'Reading from {args.infile}, writing to {args.outfile}.png')
    if args.pnnmu_cut is not None:
        logging.info(f'Applying PROBNNmu cut: {args.pnnmu_cut}')

    masses = run_event_loop(args)

    total = sum(len(v) for v in masses.values())
    for label in CATEGORIES:
        n = len(masses.get(label, []))
        if n:
            logging.info(f'  {label + ":":<35} {n:6,d}  ({100*n/total:.1f}%)')

    plot_masses(masses, args.outfile)
    logging.info('DONE')


if __name__ == '__main__':
    main()
