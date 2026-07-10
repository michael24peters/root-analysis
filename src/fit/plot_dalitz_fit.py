"""
plot_dalitz_fit.py

Plot the bin-by-bin Dalitz fit yield map produced by fit_dalitz.py.

Same visual style as src/plot/plot_dalitz.py (viridis pcolormesh, same axis
labels), but colored by fitted signal yield (n_sig) per bin instead of raw
candidate counts, and bins where the mass fit didn't run or didn't converge
are masked white -- same idea as plot_dalitz.py masking low-count bins, but
here it's "no reliable yield" rather than "too few raw counts".

Usage:
    python plot_dalitz_fit.py results_dalitz.json [--output dalitz_fit.png]
"""

import argparse
import json
import os
import numpy as np
import matplotlib.pyplot as plt

# PDG masses [GeV], used only to draw the physical Dalitz-boundary reference
# lines below -- real eta -> mu mu gamma signal is kinematically confined to:
#   m12 = m^2(gamma, mu+) in [m_mu^2, (m_eta - m_mu)^2]
#   m23 = m^2(mu+, mu-)   in [(2 m_mu)^2, m_eta^2]
# (background candidates can and do land outside this box; that's why
# fit_dalitz.py's default grid is padded a bit beyond it rather than cut
# exactly at it -- see its docstring)
M_MU = 0.105658
M_ETA = 0.547862
M12_PHYS = (M_MU**2, (M_ETA - M_MU)**2)
M23_PHYS = ((2 * M_MU)**2, M_ETA**2)


def parse_args():
    parser = argparse.ArgumentParser(description="Plot a bin-by-bin Dalitz fit yield map.")
    parser.add_argument("input", help="JSON produced by fit_dalitz.py")
    parser.add_argument("--output", default=None,
                        help="Output PNG path (default: dalitz_fit.png next to the input JSON)")
    return parser.parse_args()


args = parse_args()

outfile = args.output or os.path.join(
    os.path.dirname(os.path.abspath(args.input)), "dalitz_fit.png")
print(f"[INFO] Reading from {args.input} and writing to {outfile}.")

with open(args.input) as f:
    result = json.load(f)

meta = result["meta"]
bins = result["bins"]
n12, n23 = meta["n12"], meta["n23"]

# fit_dalitz.py writes bins in row-major (i over m12, j over m23) order, so
# the edges can be read straight off the grid without re-matching floats.
m12_edges = [bins[i * n23]["m12_lo"] for i in range(n12)] + [bins[(n12 - 1) * n23]["m12_hi"]]
m23_edges = [bins[j]["m23_lo"] for j in range(n23)] + [bins[n23 - 1]["m23_hi"]]

# NaN for bins that weren't fit (too few entries) or didn't converge, so
# pcolormesh can mask them white like plot_dalitz.py does for low counts.
n_sig = np.array([b["n_sig"] if b["valid"] else np.nan for b in bins]).reshape(n12, n23)
n_sig_masked = np.ma.masked_invalid(n_sig)

fig, ax = plt.subplots(figsize=(7, 6))
mesh = ax.pcolormesh(m12_edges, m23_edges, n_sig_masked.T, cmap="viridis")
fig.colorbar(mesh, ax=ax, label="Fitted signal yield")

# Mark the physical signal boundary (see M12_PHYS/M23_PHYS above).
ax.axvline(M12_PHYS[0], color="red", ls="--", lw=1, label="physical signal window")
ax.axvline(M12_PHYS[1], color="red", ls="--", lw=1)
ax.axhline(M23_PHYS[0], color="red", ls="--", lw=1)
ax.axhline(M23_PHYS[1], color="red", ls="--", lw=1)
ax.legend(loc="upper right", fontsize=8)

ax.set_xlabel(r"$m^2_{\gamma,\mu^+}$ [GeV$^2$]")
ax.set_ylabel(r"$m^2_{\mu^+,\mu^-}$ [GeV$^2$]")
ax.set_title(r"$\eta \to \mu\mu\gamma$ Dalitz plot (bin-by-bin fitted yield)")

fig.tight_layout()
fig.savefig(outfile, dpi=150)
print(f"Saved: {outfile}")
