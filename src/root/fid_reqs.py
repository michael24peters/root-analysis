"""
fid_reqs.py

Script to apply fiducial requirements to generator-level particles in a ROOT
file. Uses uproot + awkward for I/O (no PyROOT/ROOT dependency).

TODO: make this work for eta -> 4mu and eta -> 2mu 2e cases as well. Requires
more generalized code.

Usage:
    python fid_reqs.py input.root [--gamma]
"""

import argparse
import logging

import awkward as ak
import numpy as np
import uproot

# -------------------------------------------------------------------------------


def pseudorapidity(px, py, pz):
    """
    Vectorized pseudorapidity, matching ROOT.TVector3::PseudoRapidity() but
    using np.where instead of try/except to handle division by zero.
    Source: https://root.cern.ch/doc/master/classTVector3.html#aedc6fc6f5f6f3f3d4e4f4e4f4f4f4f4
    Code: https://root.cern.ch/doc/master/TVector3_8cxx_source.html#l00345
    """
    px, py, pz = np.asarray(px, float), np.asarray(py, float), np.asarray(pz, float)
    p = np.sqrt(px**2 + py**2 + pz**2)
    cos_theta = np.where(p != 0, pz / np.where(p != 0, p, 1.0), 1.0)  # avoid div by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        log_eta = -0.5 * np.log((1.0 - cos_theta) / (1.0 + cos_theta))
    edge_eta = np.where(pz == 0, 0.0, np.where(pz > 0, 1e10, -1e10))
    return np.where(cos_theta**2 < 1, log_eta, edge_eta)


def passes_reqs(pid, px, py, pz):
    """Vectorized fiducial requirement check for one or more particles.

    Fiducial requirements:
    - Pseudorapidity (eta) in [2, 4.5]
    - Muon pT > 500 MeV
    - Muon P > 3 GeV
    - Photon pT > 500 MeV
    """
    pid = np.asarray(pid)
    px, py, pz = np.asarray(px, float), np.asarray(py, float), np.asarray(pz, float)
    p = np.sqrt(px**2 + py**2 + pz**2)  # momentum
    pt = np.sqrt(px**2 + py**2)  # transverse momentum
    eta = pseudorapidity(px, py, pz)

    is_muon = np.abs(pid) == 13
    is_photon = pid == 22
    is_eta = pid == 221

    muon_ok = (2.0 < eta) & (eta < 4.5) & (pt > 500) & (p > 3000)
    photon_ok = (pt > 500) & (2.0 < eta) & (eta < 4.5)

    ok = np.select([is_muon, is_photon, is_eta], [muon_ok, photon_ok, True], default=False)
    return ok & (p != 0)  # prevent division by zero (wouldn't pass cuts anyway)


# -------------------------------------------------------------------------------


def event_passes(pids, px, py, pz, expected_pids):
    """
    Group one event's flat generator-level daughter arrays into candidates of
    len(expected_pids) particles each, and return True if the last candidate
    matching `expected_pids` has all daughters pass the fiducial requirements.

    Preserves the original script's per-event semantics: an event with more
    than one matching candidate (very rare, generator-level MC normally has
    exactly one) keeps only the last matching candidate's pass/fail status
    rather than requiring every candidate to pass -- see the similar
    single-candidate caveat in bkg.py's classify_candidate.
    """
    stride = len(expected_pids)
    passed = True
    for start in range(0, len(pids) - stride + 1, stride):
        group = pids[start:start + stride]
        if group != expected_pids:
            continue
        passed = bool(np.all(passes_reqs(
            group,
            px[start:start + stride],
            py[start:start + stride],
            pz[start:start + stride],
        )))
    return passed


def apply_fiducial_reqs(events, gamma_flag):
    """
    Apply fiducial cuts to generator-level particles for every event.

    `events` is an awkward record array with (at least) the mc_pid, mc_px,
    mc_py, mc_pz fields. Returns a boolean numpy array, one entry per event,
    True where the event passes.
    """
    expected = [221, -13, 13, 22] if gamma_flag else [221, -13, 13]

    # Pull the needed branches into plain Python lists up front -- the
    # per-event candidate grouping below is inherently ragged (variable
    # number of candidates per event), so it stays a Python loop, but
    # operating on pre-loaded lists instead of PyROOT's per-entry TTree
    # reads. The fiducial math itself (pseudorapidity/passes_reqs) runs
    # vectorized over each matching candidate's 3-4 particles at once.
    mc_pid = ak.to_list(events['mc_pid'])
    mc_px = ak.to_list(events['mc_px'])
    mc_py = ak.to_list(events['mc_py'])
    mc_pz = ak.to_list(events['mc_pz'])

    n_events = len(mc_pid)
    mask = np.zeros(n_events, dtype=bool)
    logging.info(f'Entries: {n_events}')
    for i in range(n_events):
        if i % 100000 == 0 and i > 0:
            logging.info(f'Processed {i:,d} events, kept {int(mask[:i].sum()):,d}...')
        pids = [int(pid) for pid in mc_pid[i]]
        mask[i] = event_passes(pids, mc_px[i], mc_py[i], mc_pz[i], expected)
    return mask


# -------------------------------------------------------------------------------

# Set up logging to stderr only
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
logging.getLogger().addHandler(logging.StreamHandler())

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('infile', help='Input ROOT file')
# eta -> mu+ mu- (gamma) flag
parser.add_argument(
    '-g', '--gamma',
    action='store_true',
    help='eta -> mu+ mu- gamma flag',
)
args = parser.parse_args()

infile = args.infile
# Put outfile in same directory as infile with fixed name
outfile = '/'.join(infile.split('/')[:-1]) + '/fiducial_requirements.root'
logging.info(f'Reading from {infile}, writing to {outfile}.')

with uproot.open(infile) as fin:
    events = fin['tree'].arrays(library='ak')

mask = apply_fiducial_reqs(events, args.gamma)
kept = events[mask]

logging.info(f'Total kept entries: {len(kept)}')

with uproot.recreate(outfile) as fout:
    fout['tree'] = kept

logging.info(f'Done: wrote reduced tree with fiducial requirements to {outfile}.')
