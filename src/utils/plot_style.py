"""
plot_style.py

Repeatable plot style for figures: text-box stat summaries, text placement (auto
or manual), and summary statistics both to terminal and matplotlib.
"""

import numpy as np

# Axes-fraction (x, y, ha, va) for each named placement. "auto" placement
# (see choose_best_position) picks among a caller-supplied subset of these.
POSITIONS = {
    'upper left':   (0.025, 0.95, 'left',   'top'),
    'upper center': (0.5,   0.95, 'center', 'top'),
    'upper right':  (0.975, 0.95, 'right',  'top'),
    'middle left':  (0.025, 0.5,  'left',   'center'),
    'middle right': (0.975, 0.5,  'right',  'center'),
    'lower left':   (0.025, 0.05, 'left',   'bottom'),
    'lower center': (0.5,   0.05, 'center', 'bottom'),
    'lower right':  (0.975, 0.05, 'right',  'bottom'),
}

def choose_best_position(ax, plot_text, curve_x, curve_y, loc_candidates=None):
    """
    Pick the position (from `loc_candidates`, default: all of POSITIONS)
    whose rendered stats box overlaps the curve the least.

    This actually draws `plot_text` at each candidate position, measures its
    real rendered bounding box via the canvas renderer (so it automatically
    accounts for font size, number of lines, and box padding), and counts how
    many (curve_x, curve_y) points fall inside it.

    Ties resolve to whichever position is listed first in
    `loc_candidates` / POSITIONS, so upper positions are preferred over
    middle/lower ones when multiple are equally clear.
    """
    loc_candidates = loc_candidates or list(POSITIONS.keys())
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    xf, yf = ax.transAxes.inverted().transform(
        ax.transData.transform(np.column_stack([curve_x, curve_y]))
    ).T
    scores = {}
    for pos in loc_candidates:
        px, py, ha, va = POSITIONS[pos]
        artist = ax.text(px, py, plot_text, transform=ax.transAxes,
                          horizontalalignment=ha, verticalalignment=va,
                          bbox=dict(boxstyle='round,pad=0.4'))
        bbox_axes = artist.get_window_extent(renderer=renderer).transformed(
            ax.transAxes.inverted())
        inside = ((xf >= bbox_axes.x0) & (xf <= bbox_axes.x1) &
                  (yf >= bbox_axes.y0) & (yf <= bbox_axes.y1))
        scores[pos] = int(inside.sum())
        artist.remove()
    return min(scores, key=scores.get)

def build_stats_text(stats, extra=None):
    """
    Build the matplotlib (LaTeX) and terminal (plain Unicode) renderings of a
    fit's stats box from one shared spec, so each symbol is written once.

    stats: list of "rows"; each row is a list of one or more
        (latex_symbol, plain_symbol, value, error, fmt, unit) tuples. Tuples
        within a row are joined with ", " onto one line (so e.g. alphaL and
        nL can share a row, matching a compact stats box layout).
    extra: optional (plot_line, term_line) tuple appended as a final row,
        for values that aren't a value/error pair (e.g. chi2/ndof).

    Returns (plot_text, term_text).
    """
    def render_row(row, latex):
        """
        Render a single row of the stats box for latex or terminal output.
        """
        parts = []
        for sym_latex, sym_plain, value, error, fmt, unit in row:
            sym = f'${sym_latex}$' if latex else sym_plain
            unit_str = f' {unit}' if unit else ''
            parts.append(f'{sym} = {value:{fmt}} ± {error:{fmt}}{unit_str}')
        return ', '.join(parts)

    plot_lines = [render_row(row, latex=True) for row in stats]
    term_lines = [render_row(row, latex=False) for row in stats]
    if extra is not None:
        # Append the extra row the output exactly as given (no val/err format).
        plot_lines.append(extra[0])
        term_lines.append(extra[1])
    return '\n'.join(plot_lines), '\n'.join(term_lines)

def add_stats_box(ax, plot_text, loc='auto', show=True, x=None, curve=None,
                   loc_candidates=None):
    """
    Draw a boxed text annotation on `ax` at `loc` (a POSITIONS key, or "auto"
    to resolve via choose_best_position(ax, plot_text, x, curve,
    loc_candidates)). Does nothing if show=False or plot_text is empty --
    this is the "statistics / no statistics" toggle.
    """
    if not show or not plot_text:
        return
    if loc == 'auto':
        loc = choose_best_position(ax, plot_text, x, curve,
                                    loc_candidates=loc_candidates)
    tx, ty, ha, va = POSITIONS[loc]
    ax.text(tx, ty, plot_text,
            transform=ax.transAxes,
            horizontalalignment=ha,
            verticalalignment=va,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor='0.3', alpha=0.85))

def print_summary(text, width=48, title='Fit summary'):
    """Print a summary to the terminal."""
    print()
    print('─' * width)
    print(f'{title:^{width-4}}')
    print('─' * width)
    print(text)
    print('─' * width)

def save_plot(fig, output, dpi=300):
    """Save the plot to disk."""
    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    print(f'[DONE] Plot saved to {output}.')
