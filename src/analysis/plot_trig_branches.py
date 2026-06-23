# Quick diagnostic: plot all trigger branches as histograms to PNG.
# Usage: python src/analysis/plot_trig_branches.py input.root

import sys
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

branches = [
    "tag_l0_tos0", "tag_l0_tos1",
    "tag_hlt1_tos0", "tag_hlt1_tos1",
    "tag_hlt2_tos0", "tag_hlt2_tos1", "tag_hlt2_tos2", "tag_hlt2_tis",
]

with uproot.open(sys.argv[1]) as f:
    tree = f["tree"]
    for name in branches:
        flat = ak.to_numpy(ak.flatten(tree[name].array(library="ak")))
        fig, ax = plt.subplots()
        ax.hist(flat, bins=50)
        ax.set_title(name)
        ax.set_xlabel("value")
        ax.set_ylabel("candidates")
        fig.savefig(f"{name}.png", dpi=100)
        plt.close(fig)
        print(f"{name}: min={flat.min():.3f} max={flat.max():.3f} mean={flat.mean():.4f} nonzero={np.sum(flat != 0)}")
