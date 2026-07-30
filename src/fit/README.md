# README

## Pipeline

```
driver_*.py   read ROOT → apply selection → reduce to 1 candidate/event → numpy masses
   ↓ numpy array
fit_*.py      fit the numbers. Blind to files, cuts, candidates.
   ↓ result dict
plot_*.py     render
```

Awkward arrays own everything jagged (ROOT I/O, event/candidate structure) and
live only in `src/utils/cut_utils.py`. Numpy owns everything flat (histogram,
PDFs, curves) and lives in `src/utils/fit_utils.py` and the `fit_*` modules.
The boundary is one `ak.to_numpy(ak.flatten(...))` call in the driver.

The drivers are the only entry point -- `fit_dcb.py` / `fit_gauss.py` are
library modules with no CLI, since they take an array rather than a file.

## Usage

One-shot fit + plot. Every flag shown explicitly, with its default value:

For double crystal ball:

```
python src/fit/driver_dcb.py data.root \
    --outdir out \
    --max-events 500000 \
    --xmin 480 --xmax 620 --nbins 80 \
    --mean-init 547.9 --sigma-init 10.0 \
    --alphaL-init 1.5 --nL-init 2.0 \
    --alphaR-init 1.5 --nR-init 2.0
```

For gaussian:

```
python src/fit/driver_gauss.py data.root \
    --outdir out \
    --max-events 500000 \
    --xmin 480 --xmax 620 --nbins 80 \
    --mean-init 547.9 --sigma-init 10.0
```

Each driver runs the fit, saves `results_dcb.json` / `results_gauss.json`,
and writes `dcb_mass_fit.png` / `gauss_mass_fit.png` into `--outdir`.
`--max-events` caps how many events are read, for quick iteration on large
files; omit it to read everything.

The drivers only expose fit hyperparameters as flags -- plot styling (title,
font, stats box, etc.) isn't driver-configurable; it comes from each plot
script's `DEFAULT_CONFIG` (see "Overriding plot style" below).

## Changing the selection

Selection lives at the top of each driver as plain module-level constants, not
CLI flags -- editing them is the intended way to change what gets fit:

```python
MASS_BRANCH     = 'tag_dtf_m'      # the quantity being fit
METRIC_BRANCH   = 'tag_dtf_chi2'   # ranks candidates within an event
DAUGHTER_FIELDS = ['prt_pid', 'prt_pnn_mu']
PNN_MU_MIN      = 0.95
```

The cut helpers in `cut_utils` return **per-candidate boolean masks** shaped
like `tag_dtf_m`, so adding a cut is one line:

```python
keep = cut_utils.tos_mask(tree) & cut_utils.pnn_mask(tree, PNN_MU_MIN)
```

Two things to know before editing that block:

- **Never index the whole tree with a per-candidate mask.** `tree[mask]`
  silently truncates the `prt_*` fields -- a 3-candidate / 9-daughter event
  comes back with 3 candidates and 3 daughters, no error. Mask individual
  `tag_*` branches instead, which is what the drivers do.
- **Drop empty events before the best-candidate step.** `ak.argmin` on an
  event with no surviving candidates returns `None`, which becomes a silently
  masked numpy entry. The drivers filter on `ak.num(metric) > 0` first.

Cuts are applied *before* the best-candidate reduction, so the survivor is the
best among candidates that passed rather than an event being discarded because
its arbitrarily-chosen candidate failed. Reordering is a matter of moving two
lines in the driver.

`drop_malformed()` removes events with no candidates and events whose `prt_*`
branches are inconsistent with `prt_idx_mom` (~0.1% of current ntuples), and
reports the count to stderr.

## Known ntuple quirks

Open items, surfaced by the drivers on stderr every run. Each is worked around
rather than solved, and each is worth chasing upstream in the ntuple production.

**Ragged HLT2 trigger branches.** `tag_hlt2_tos0` / `tag_hlt2_tos1` are not
always written with one entry per candidate. Two distinct cases:

- *Absent* (empty list): the line was not evaluated for that event. Treated as
  "did not fire" and zero-filled by `_align_candidates()`. Reported as
  `[INFO] <branch>: not evaluated in N/M events`. If a branch is ever absent in
  *all* events it contributes nothing to the TOS OR -- i.e. that trigger is
  silently not part of the selection -- so it is escalated to `[WARN]`.
- *Partial* (fewer entries than candidates): there is no way to tell which
  candidate each entry belongs to, so `drop_malformed()` drops those events
  rather than assuming positional correspondence. Reported as
  `[WARN] <branch>: partial in N events`.

Measured on `ntuple/lhcb/.../20260716/combined_files.root`:

| channel | partial (events dropped) | `tag_hlt2_tos0` absent | `tag_hlt2_tos1` absent |
|---|---|---|---|
| `eta2mumu` | none | none | none |
| `eta2mumugamma` | 16.8% | 88.3% | 52.1% |

So `eta2mumu` is clean, while `eta2mumugamma` loses ~17% of events to the
ambiguity and runs with HLT2 contributing to the TOS OR only where it was
actually evaluated. **Unresolved:** why the HLT2 branches are written per event
rather than per candidate, and what a partial list is meant to signify.

**Inconsistent `prt_*` lengths.** In a small fraction of events (~0.1% in MC)
the `prt_*` branches disagree in length with `prt_idx_mom`, or an index points
past the last candidate. Those events are dropped and counted.

> **Caveat on the best-candidate metric.** In current ntuples `tag_dtf_chi2` is
> identical across all candidates in ~61% of multi-candidate events, with a
> median mass spread of ~300 MeV within those events -- so `method='min'` is
> effectively picking candidate 0 arbitrarily most of the time.
> `find_best_candidates()` takes `metric=` / `method=` so a better discriminant
> can be substituted.

## Overriding plot style

Plot as a standalone CLI (only takes the JSON path -- uses `plot_dcb.py`'s
`DEFAULT_CONFIG` as-is):

```
python src/plot/plot_dcb.py out/results_dcb.json
```

To override any plot style value, either edit `DEFAULT_CONFIG` in
`plot_dcb.py` directly, or call `plot_dcb()` as a function with every value
defined explicitly:

```python
import sys
sys.path.insert(0, "src/plot")
sys.path.insert(0, "src/utils")
import plot_utils
from plot_dcb import plot_dcb

data = plot_utils.load_data("out/results_dcb.json")
plot_dcb(data, {
    "output":     "out/dcb_mass_fit.png",
    "title":      "Eta Mass",
    "font":       "sans",          # "serif" (STIXGeneral) or "sans" (Inter)
    "x_lim":      (420, 680),      # or None to auto-range from the data
    "pull_ylim":  (-4, 4),         # or None to auto-range from the pulls
    "text_loc":   "upper left",    # "auto", or any plot_style.POSITIONS key
                                    # (see "Config keys" below for the full list)
    "show_stats": True,            # False to omit the stats box entirely
    "legend_loc": "upper right",   # any matplotlib legend loc string
})
```

`plot_gauss` works identically, with `plot_gauss.py`'s own `DEFAULT_CONFIG`.

Each fit type's config stays entirely in that fit type's own file -- nothing is shared or
copied between DCB and Gauss, so their titles/ranges/inits can diverge independently.

## Fitting an arbitrary array

Because the fit modules take a numpy array, they can fit anything -- not just
what a driver produced:

```python
import sys
sys.path.insert(0, "src/fit")
from fit_dcb import fit_dcb

result = fit_dcb(masses, xmin=480, xmax=620, nbins=80)   # masses: 1-D ndarray
```

`result['meta']` records `n_input` (how many values were fit) but no filename;
the driver adds `input_file` and `selection` afterwards for provenance.

### Config keys

These are the keys `plot_fit_data()` reads out of a config dict (whether that
dict is `DEFAULT_CONFIG` as-is, or `DEFAULT_CONFIG` with overrides layered on
top via `plot_dcb(data, {...})` / `plot_gauss(data, {...})`):

- `output`: path to save the PNG to.
- `title`: plot title.
- `x_lim` (optional): `(xmin, xmax)` tuple; omit or set `None` to auto-range
  from the data.
- `pull_ylim` (optional): `(ymin, ymax)` tuple for the pull panel; omit or
  set `None` to auto-range from the pulls.
- `font` (optional): `"serif"` (default, STIXGeneral) or `"sans"` (Inter,
  falls back to a generic sans-serif if Inter isn't installed).
- `show_stats` (optional, default `True`): set `False` to omit the stats box
  entirely.
- `text_loc` (optional, default `"auto"`): where the stats box goes.
  `"auto"` picks the emptiest of the upper positions by default; or pass any
  `plot_style.POSITIONS` key directly for manual placement: `"upper left"`,
  `"upper center"`, `"upper right"`, `"middle left"`, `"middle right"`,
  `"lower left"`, `"lower center"`, `"lower right"`.
- `legend_loc` (optional): any matplotlib legend `loc` string, default `"best"`.
- `build_text`: function `(params, chi2_per_ndof) -> (plot_text, term_text)`
  that builds the stats box content -- see `build_text()` in `plot_dcb.py` /
  `plot_gauss.py` for how it's assembled via `plot_style.build_stats_text()`.

### How config flows

Each layer only knows about the layer directly below it, and only passes down what it
needs to override -- no script reaches into another script's config object:

- `driver_dcb.py` / `driver_gauss.py` own the fit hyperparameters (as CLI flags) and call
  `fit_dcb()` / `fit_gauss()` directly. They then call `plot_dcb(result, {"output": ...})`
  / `plot_gauss(result, {"output": ...})`, passing only the one override they need.
- `plot_dcb.py` / `plot_gauss.py` own `DEFAULT_CONFIG` (styling: title, font, stats-box
  placement, `build_text`) and expose a `plot_dcb(data, config=None)` / `plot_gauss(data,
  config=None)` function that layers any caller overrides on top of `DEFAULT_CONFIG` before
  calling into the shared engine. Running `python plot_dcb.py results.json` standalone (no
  driver) just uses `DEFAULT_CONFIG` as-is.
- `plot_mass_fit.py` is the shared rendering engine (`plot_fit_data(data, config)`) -- it owns
  no defaults of its own, it just renders whatever config dict it's handed. It builds on two
  shared helper modules: `src/utils/plot_utils.py` (generic JSON I/O, font setup, CLI arg
  parsing) and `src/utils/plot_style.py` (the stats-box/placement/save-and-report "house
  style," shared by any current or future fit-plot type).
