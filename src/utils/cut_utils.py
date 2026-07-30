"""
cut_utils.py

ROOT I/O and candidate-level selection. This is the only module in the fitting
pipeline that knows about files, trees, or jagged structure -- everything
downstream (fit_utils, fit_*, plot_*) works in flat numpy.

How the ntuple is laid out
--------------------------
An event contains zero or more candidates, and each candidate has 2+ daughters
(depending on the decay). The file does NOT store it that way. Taking a real
event with three candidates:

    tag_dtf_m    = [1080.5, 551.5, 1022.0]      3 candidates
    prt_pnn_mu   = [0.42, 0.947, -100,          9 daughters, ONE FLAT LIST
                    0.42, 0.947, -100,          (all candidates run together)
                    0.42, 0.947, -100]
    prt_idx_mom  = [0,0,0, 1,1,1, 2,2,2]        parent candidate of each daughter

`tag_*` is one value per candidate; `prt_*` is one value per daughter, with every
candidate's daughters concatenated into a single per-event list. There is no
candidate level in the stored `prt_*` -- `prt_idx_mom` is the only thing that says
where one candidate's daughters end and the next begin. `_group_daughters()`
uses it to insert that missing level, which is what lets a per-daughter quantity
(muon PID) become a per-candidate decision.

Selection convention
--------------------
The cut functions return **per-candidate boolean masks**, shaped like `tag_dtf_m`,
so a caller composes them with plain `&`. They deliberately do not filter the tree
themselves: indexing a whole tree with a per-candidate mask silently truncates the
`prt_*` fields (a 3-candidate/9-daughter event comes back with 3 candidates and 3
daughters, no error). Mask individual `tag_*` branches instead.
"""

import copy
import sys

import awkward as ak
import numpy as np

# Diagnostics accumulated as a run progresses -- events dropped and why, which
# trigger lines were actually available. These are printed to stderr as they
# happen, and collected here so a driver can persist them next to its results
# instead of leaving them in terminal scrollback. One process = one run, so a
# module-level record is enough; call reset_report() to start a fresh one.
_REPORT = {}


def report():
    """Diagnostics collected so far, safe to embed in a results file."""
    return copy.deepcopy(_REPORT)


def reset_report():
    """Discard collected diagnostics (only needed if a process does two runs)."""
    _REPORT.clear()

# Trigger-on-signal branches OR'd together by `tos_mask`. These are stored as
# float64 0.0/1.0 rather than bool, hence the `> 0.5` comparison there.
TOS_TRIGGERS = ['tag_l0_tos0', 'tag_l0_tos1', 'tag_hlt1_tos0',
                'tag_hlt1_tos1', 'tag_hlt2_tos0', 'tag_hlt2_tos1']


def read_branches(root_path, branches, tree='tree', entry_stop=None):
    """
    Read `branches` as awkward arrays. The single place library= is chosen.

    entry_stop : optional cap on events read, for quick iteration on large files.
    """
    import uproot
    with uproot.open(root_path) as f:
        return f[tree].arrays(list(branches), library='ak', entry_stop=entry_stop)


def drop_malformed(tree, daughter_fields=(), candidate_fields=(),
                   tag_ref='tag_dtf_m', idx_branch='prt_idx_mom'):
    """
    Drop events that cannot be interpreted, so everything downstream can trust
    the candidate structure. Three kinds go, each counted separately on stderr:

      - events with no candidates: nothing to select, and they would turn into
        None entries at the best-candidate step
      - events whose `daughter_fields` are inconsistent with `prt_idx_mom` --
        a differing length, or an index pointing past the last candidate
      - events where a `candidate_fields` branch holds a partial number of
        entries, i.e. neither one per candidate nor none at all. Some trigger
        branches are written per event only when the line was evaluated, so an
        empty list is meaningful (`tos_mask` reads it as "did not fire"), but a
        partial list gives no way to tell which candidate each entry belongs
        to, so those events are dropped rather than guessed at.

    Reports counts rather than failing: these are upstream production
    artefacts, and how many there are varies a lot by channel.
    """
    n_before = len(tree)
    n_tag = ak.num(tree[tag_ref])

    has_candidates = n_tag > 0
    well_formed_daughters = ak.all(tree[idx_branch] < n_tag, axis=1)
    for field in daughter_fields:
        well_formed_daughters = well_formed_daughters & (
            ak.num(tree[field]) == ak.num(tree[idx_branch]))
    whole_candidate_fields = ak.ones_like(n_tag, dtype=bool)
    partial_by_field = {}
    for field in candidate_fields:
        n_field = ak.num(tree[field])
        whole = (n_field == n_tag) | (n_field == 0)
        partial_by_field[field] = int(ak.sum(has_candidates & ~whole))
        whole_candidate_fields = whole_candidate_fields & whole

    keep = has_candidates & well_formed_daughters & whole_candidate_fields
    n_kept = int(ak.sum(keep))
    _REPORT['events'] = {
        'read': n_before,
        'kept': n_kept,
        'dropped_total': n_before - n_kept,
        'dropped_no_candidates': int(ak.sum(~has_candidates)),
        'dropped_malformed_daughters': int(ak.sum(has_candidates & ~well_formed_daughters)),
        'dropped_partial_candidate_branch': int(ak.sum(has_candidates & ~whole_candidate_fields)),
        'partial_by_branch': {f: n for f, n in partial_by_field.items() if n},
    }
    print(f'[INFO] {n_before - n_kept} events dropped, {n_kept} kept '
          f'(no candidates: {int(ak.sum(~has_candidates))}, '
          f'malformed {idx_branch}: {int(ak.sum(has_candidates & ~well_formed_daughters))}, '
          f'partial candidate branch: {int(ak.sum(has_candidates & ~whole_candidate_fields))})',
          file=sys.stderr)
    # Name the offending branches: a branch that is chronically partial is a
    # data-quality problem worth chasing upstream, not just a few lost events.
    for field, n_partial in sorted(partial_by_field.items(), key=lambda kv: -kv[1]):
        if n_partial:
            print(f'[WARN] {field}: partial in {n_partial} events '
                  f'({100.0 * n_partial / n_before:.1f}%) -- those events dropped, '
                  f'entries cannot be matched to candidates', file=sys.stderr)
    return tree[keep]


def _align_candidates(tree, field, tag_ref='tag_dtf_m'):
    """
    Expand a per-candidate branch to exactly one value per candidate, filling 0
    for events where the branch is absent (stored as an empty list).

    Assumes `drop_malformed()` has already run, so each event's branch holds
    either one entry per candidate or none at all.
    """
    n_tag = ak.num(tree[tag_ref])
    # Per-candidate flag: does this candidate's event carry the branch at all?
    present = ak.num(tree[field]) == n_tag
    present_per_candidate = ak.to_numpy(
        ak.flatten(ak.broadcast_arrays(present, tree[tag_ref])[0]))

    values = np.zeros(int(ak.sum(n_tag)))
    values[present_per_candidate] = ak.to_numpy(ak.flatten(tree[field]))
    return ak.unflatten(values, n_tag)


def _group_daughters(tree, field, tag_ref='tag_dtf_m', idx_branch='prt_idx_mom'):
    """
    Insert the missing candidate level into a per-daughter `prt_*` field, so
    [event][daughter] becomes [event][candidate][daughter]. See the module
    docstring for why this level is absent in the first place.

    Daughters are grouped by consecutive runs of `prt_idx_mom`. That index
    restarts at 0 in every event, so it is first offset by the number of
    candidates in all preceding events -- otherwise the last candidate of one
    event and the first of the next would merge into a single run.

    Assumes `drop_malformed()` has already run.
    """
    n_tag = ak.num(tree[tag_ref])
    if len(tree) == 0: return tree[field]

    counts_per_event = ak.to_numpy(n_tag)
    # Candidates in all events before this one; makes each candidate id global.
    offset = np.concatenate([[0], np.cumsum(counts_per_event)[:-1]])
    counts = ak.run_lengths(ak.flatten(tree[idx_branch] + offset))

    # Two unflattens: daughters -> candidates, then candidates -> events.
    per_candidate = ak.unflatten(ak.flatten(tree[field]), counts)
    return ak.unflatten(per_candidate, n_tag)


def tos_mask(tree, triggers=TOS_TRIGGERS, tag_ref='tag_dtf_m'):
    """
    Per-candidate mask: True where the candidate is TOS for at least one trigger.

    Trigger branches hold 0.0/1.0 as float64, so `> 0.5` recovers the boolean.
    They are also aligned first: an HLT2 line that was not evaluated for an
    event is stored as an empty list rather than a row of zeros, and counts as
    not fired.
    """
    n_tag = ak.num(tree[tag_ref])
    n_events = len(tree)
    not_evaluated, ignored = {}, []
    mask = None
    for trigger in triggers:
        # An absent branch contributes nothing to the OR, so a line that is
        # absent everywhere is silently not part of the selection at all.
        # Report it: that is a trigger being ignored, and it should be chased
        # rather than assumed intentional.
        n_absent = int(ak.sum(ak.num(tree[trigger]) != n_tag))
        if n_absent:
            not_evaluated[trigger] = n_absent
            level = 'WARN' if n_absent == n_events else 'INFO'
            if n_absent == n_events: ignored.append(trigger)
            print(f'[{level}] {trigger}: not evaluated in {n_absent}/{n_events} events '
                  f'({100.0 * n_absent / n_events:.1f}%), counted as not fired'
                  + (' -- this trigger contributes nothing here'
                     if n_absent == n_events else ''), file=sys.stderr)
        fired = _align_candidates(tree, trigger, tag_ref) > 0.5
        mask = fired if mask is None else (mask | fired)

    _REPORT['triggers'] = {
        'requested': list(triggers),
        # Lines fully available everywhere; these are the ones actually driving
        # the OR in every event.
        'always_evaluated': [t for t in triggers if t not in not_evaluated],
        # branch -> number of events where it was absent and read as not fired
        'not_evaluated_in_events': not_evaluated,
        # Absent in every event, so contributing nothing to the selection
        'ignored': ignored,
    }
    return mask


def pnn_mask(tree, threshold=0.95, pnn_field='prt_pnn_mu', pid_field='prt_pid'):
    """
    Per-candidate mask: True where *every* muon daughter has pnn_mu > threshold.

    Only muons (|pid| == 13) are tested -- photons carry a -100 sentinel in
    prt_pnn_mu and would fail any real threshold.

    `ak.all(..., axis=-1)` collapses the innermost (daughter) level, turning a
    per-daughter comparison into one bool per candidate, shaped like tag_dtf_m.
    """
    pid = _group_daughters(tree, pid_field)
    pnn = _group_daughters(tree, pnn_field)
    is_muon = abs(pid) == 13
    return ak.all(pnn[is_muon] > threshold, axis=-1)


def selection_mask(tree, use_tos=True, pnn_mu_min=0.95, tag_ref='tag_dtf_m'):
    """
    Per-candidate mask combining the standard cuts, so the fit drivers and the
    diagnostic plotters select identically.

    use_tos     : require the candidate to be TOS for at least one trigger
    pnn_mu_min  : require every muon daughter above this; None drops the cut

    With both disabled every candidate passes, which is the "no selection"
    baseline in a selection comparison.
    """
    mask = ak.ones_like(tree[tag_ref], dtype=bool)
    stages = [('all candidates', int(ak.sum(ak.num(mask))))]
    if use_tos:
        mask = mask & tos_mask(tree, tag_ref=tag_ref)
        stages.append(('TOS', int(ak.sum(mask))))
    if pnn_mu_min is not None:
        mask = mask & pnn_mask(tree, pnn_mu_min)
        stages.append((f'pnn_mu > {pnn_mu_min}', int(ak.sum(mask))))

    _REPORT['cuts'] = {
        'tos': use_tos,
        'pnn_mu_min': pnn_mu_min,
        # Candidates surviving after each cut is applied, in order
        'cutflow_candidates': [{'stage': s, 'candidates': n} for s, n in stages],
    }
    return mask


def find_best_candidates(metrics, method='min', target=None):
    """
    Reduce each event to its single best candidate, by `metrics`.

    metrics : per-candidate awkward array to rank on (e.g. tag_dtf_chi2),
              already filtered to the candidates that passed selection
    method  : 'min' / 'max' -- smallest / largest metric wins
              'closest'     -- metric nearest `target` wins
    target  : required for method='closest'

    Returns a per-event index into `metrics` with `keepdims=True`, i.e. a
    one-element list per event rather than a bare number. That keeps the array
    jagged when you index with it, so the result is still one-candidate-per-event
    and every other per-candidate array can be indexed the same way.

    Every event in `metrics` must have at least one candidate -- `ak.argmin` on
    an empty list yields None, which becomes a silently masked numpy entry later.
    Filter with `ak.num(metrics) > 0` first.
    """
    if method == 'min':
        return ak.argmin(metrics, axis=1, keepdims=True)
    if method == 'max':
        return ak.argmax(metrics, axis=1, keepdims=True)
    if method == 'closest':
        if target is None:
            raise ValueError("method='closest' requires a target value")
        return ak.argmin(abs(metrics - target), axis=1, keepdims=True)
    raise ValueError(f"Invalid method {method!r}. Must be 'min', 'max', or 'closest'.")


def transverse_momentum(px, py):
    """pT = sqrt(px^2 + py^2), elementwise over any-shaped px/py fields."""
    return np.sqrt(px**2 + py**2)
