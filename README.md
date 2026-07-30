# anaroot

LHCb analysis of η → μμ(γ): background classification, mass fits, and Dalitz
analysis. Fit methodology and the motivation behind the choice of signal/
background functions are documented separately (LaTeX writeup); this repo is
inline-documented (each script has a `"""usage docstring"""` at the top) and
isn't meant to duplicate that.

## Layout

```
src/root/   file prep: fiducial cuts, concatenation (PyROOT for TFileMerger; uproot elsewhere)
src/ana/    MC-truth background classification and diagnostics
src/fit/    mass fits (Gaussian, DCB, Dalitz bin-by-bin) -- see src/fit/README.md
src/plot/   plotting scripts
src/utils/  shared library code (fit math, plot styling, cuts)
ntuple/     input ROOT files
out/        generated plots, JSON fit results
```

ROOT I/O is uproot+awkward everywhere except `src/root/concat_roots.py`, which
uses PyROOT's `TFileMerger` for whole-file merging (no uproot equivalent).

## Environment

```
uv venv
source .venv/bin/activate
uv pip install -r requirements.txt
```

PyROOT itself isn't pip-installable -- it comes from the LCG/`lb-conda` stack,
needed only for `concat_roots.py`.

## Tests

```
pytest tests/
```

## Linting

Rule set/config lives in `ruff.toml` -- only quote-style is enabled by
default (single quotes for strings, double for docstrings); other rules are
listed there commented-out with what each one does. Uncomment a rule code to
turn it on.

```
uvx ruff check src/          # scan only, no changes
uvx ruff check src/ --fix    # apply auto-fixes for whatever's enabled
```

(No `ruff` install needed -- `uvx` fetches and caches it on first use. If you
don't have `uv`, `python3 -m ruff check src/` works too, as long as `ruff` is
installed in your active environment.)

This runs the *linter* (`ruff check`), not the full-file `ruff format`
formatter -- see the note at the top of `ruff.toml` before turning that on,
it reformats far more than quotes.
