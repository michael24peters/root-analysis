"""
eta_peak_fit.py
───────────────
Reads a ROOT file, extracts the `tag_m` branch (reconstructed eta candidate
masses), fits a peak model (Gaussian signal + Chebyshev polynomial background)
using PyROOT / RooFit, and saves the result as a PNG.

Usage
─────
    python eta_peak_fit.py input.root [output.png]

Positional arguments (all optional after the ROOT file):
    output     : output image file name                (default: "eta_fit.png")

Optional arguments:
    --xmin      : lower mass range for fit/plot [MeV] (default: 480)
    --xmax      : upper mass range for fit/plot [MeV] (default: 620)
    --nbins     : number of histogram bins            (default: 80)
"""

import sys
import os
import argparse

# ── argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Fit the eta peak in the tag_m branch of a ROOT file."
)
parser.add_argument("rootfile",           help="Path to the input ROOT file")
parser.add_argument("output", nargs="?", default="eta_fit.png",
                    help="Output PNG file (default: eta_fit.png)")
parser.add_argument("--xmin",  type=float, default=480.0,
                    help="Lower mass range for fit / plot [MeV] (default: 480)")
parser.add_argument("--xmax",  type=float, default=620.0,
                    help="Upper mass range for fit / plot [MeV] (default: 620)")
parser.add_argument("--nbins", type=int,   default=80,
                    help="Number of histogram bins (default: 80)")
args = parser.parse_args()

# ── ROOT import ───────────────────────────────────────────────────────────────
try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
except ImportError:
    sys.exit("PyROOT is not available. Please install ROOT with Python bindings.")

# ── load data ─────────────────────────────────────────────────────────────────
if not os.path.exists(args.rootfile):
    sys.exit(f"File not found: {args.rootfile}")

rfile = ROOT.TFile.Open(args.rootfile, "READ")
if not rfile or rfile.IsZombie():
    sys.exit(f"Could not open ROOT file: {args.rootfile}")

tree = rfile.Get("tree")
if not tree:
    keys = [k.GetName() for k in rfile.GetListOfKeys()]
    sys.exit(f"TTree 'tree' not found.  Available keys: {keys}")

# ── fill histogram from branch ────────────────────────────────────────────────
xmin, xmax = args.xmin, args.xmax
h = ROOT.TH1D("h_tagm", "tag_m;m_{#eta} (MeV);Counts", args.nbins, xmin, xmax)
h.Sumw2()  # sum of squares of weights for error bars
tree.Draw("tag_m>>h_tagm", "", "goff")   # "goff" = no graphics
n_total = int(h.GetEntries())
print(f"[info] Entries in [{xmin}, {xmax}] MeV : {n_total}")
if n_total == 0:
    sys.exit("Histogram is empty – check branch name and mass range.")

# ── RooFit model setup ────────────────────────────────────────────────────────
m = ROOT.RooRealVar("m", "m_{#eta} (MeV)", xmin, xmax)

# — signal: Gaussian —
# name, title, init, min, max
mean  = ROOT.RooRealVar("mean",  "mean",   548.0,  520.0,  580.0)   # eta PDG ≈ 547.9 MeV
sigma = ROOT.RooRealVar("sigma", "sigma",  10.0,   1.0,    40.0)
gauss = ROOT.RooGaussian("gauss", "Gaussian signal", m, mean, sigma)

# — background: 2nd-order Chebyshev polynomial —
# Chebyshev polynomial used for creating a smooth background PDF.  Mass is
# normalized to [-1, 1] range. The first three Chebyshev polynomials are:
#
# T0(x) = 1, T1(x) = x, T2(x) = 2x^2 - 1.
#
# The background PDF is then modeled by a linear combination of Chebyshev
# polynomials:
#
# B(x) = 1 + c0 * T1(x) + c1 * T2(x), 
# 
# where c0, c1 are fit params that build the background PDF, x is the normalized
# mass variable. c0 controls slope (linear) and c1 controls curvature 
# (quadratic).
c0 = ROOT.RooRealVar("c0", "c0",  0.0, -5.0, 5.0)
c1 = ROOT.RooRealVar("c1", "c1",  0.0, -5.0, 5.0)
bkg = ROOT.RooChebychev("bkg", "Chebyshev background", m, ROOT.RooArgList(c0, c1))

# — signal + background yields —
n_sig = ROOT.RooRealVar("n_sig", "signal yield",     n_total * 0.5, 0, n_total * 2)
n_bkg = ROOT.RooRealVar("n_bkg", "background yield", n_total * 0.5, 0, n_total * 2)
model = ROOT.RooAddPdf("model", "Gaussian + Chebyshev",
                       ROOT.RooArgList(gauss, bkg),
                       ROOT.RooArgList(n_sig, n_bkg))

# ── import histogram into RooFit ──────────────────────────────────────────────
data = ROOT.RooDataHist("data", "data", ROOT.RooArgList(m), h)

# ── fit ───────────────────────────────────────────────────────────────────────
# Maximum-likelihood fit for mean, sigma, n_sig, n_bkg, c0, c1
fit_result = model.fitTo(
    data,
    ROOT.RooFit.Save(True),
    ROOT.RooFit.PrintLevel(-1),
    ROOT.RooFit.Warnings(False),
)

print("\n[fit result]")
print(f"  mean  = {mean.getVal():.2f} ± {mean.getError():.2f} MeV")
print(f"  sigma = {sigma.getVal():.2f} ± {sigma.getError():.2f} MeV")
print(f"  n_sig = {n_sig.getVal():.0f} ± {n_sig.getError():.0f}")
print(f"  n_bkg = {n_bkg.getVal():.0f} ± {n_bkg.getError():.0f}")
'''
fit status codes from RooFitResult::status():
Status	Meaning
0	    Converged successfully — the fit found a reliable minimum
1	    Covariance matrix was made positive-definite (approximate errors)
2	    Hessian is invalid (errors unreliable)
3	    EDM (estimated distance to minimum) above tolerance
4	    Reached call limit before converging
'''
print(f"  fit status : {fit_result.status()}  (0 = converged)")


# ── pull plot ─────────────────────────────────────────────────────────────────

frame_diag = m.frame(ROOT.RooFit.Bins(args.nbins))
data.plotOn(frame_diag,
            ROOT.RooFit.Name("h_data"),
            ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))  # plot data on frame
model.plotOn(frame_diag, ROOT.RooFit.Name("fit_curve"))  # plot model on frame

fit_curve = frame_diag.getCurve("fit_curve")
bin_ctrs = []  # bin centre positions [MeV]
pulls = []

for i in range(1, args.nbins + 1):
    bc = h.GetBinCenter(i)  # bin center position [MeV]
    d = h.GetBinContent(i)  # observed data count
    err = h.GetBinError(i)  # sum of squares of weights error
    f = fit_curve.Eval(bc)  # model prediction at bin center

    bin_ctrs.append(bc)  # store bin centre for pull plot x-axis

    if err > 0: pulls.append((d - f) / err)
    else: pulls.append(0.0)

# ── canvas layout──────────────────────────────────────────────────────────────
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

canvas = ROOT.TCanvas("c", "Gaussian mass fit", 900, 950)

# Pad y-boundaries (with only two pads: main ~60%, pull ~40%)
Y_MAIN_BOT  = 0.40
Y_P2_BOT    = 0.00

L_MARGIN = 0.13
R_MARGIN = 0.05

pad_main  = ROOT.TPad("pad_main",  "", 0.0, Y_MAIN_BOT, 1.0, 1.00)
pad_pull = ROOT.TPad("pad_pull", "", 0.0, Y_P2_BOT,   1.0, Y_MAIN_BOT)

for pad in (pad_main, pad_pull):
    pad.SetLeftMargin(L_MARGIN)
    pad.SetRightMargin(R_MARGIN)
    pad.Draw()

pad_main.SetTopMargin(0.08)
pad_main.SetBottomMargin(0.02)    # no x-axis labels drawn here

pad_pull.SetTopMargin(0.04)
pad_pull.SetBottomMargin(0.35)   # enlarged bottom → room for x-axis title

# ─── main fit pad ─────────────────────────────────────────────────────────────
pad_main.cd()

frame = m.frame(ROOT.RooFit.Bins(args.nbins))

data.plotOn(
    frame,
    ROOT.RooFit.Name("data"),
    ROOT.RooFit.MarkerStyle(20),
    ROOT.RooFit.MarkerSize(0.85),
    ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
)

model.plotOn(
    frame,
    ROOT.RooFit.Name("total"),
    ROOT.RooFit.LineColor(ROOT.kBlue + 1),
    ROOT.RooFit.LineWidth(2),
)

model.plotOn(
    frame,
    ROOT.RooFit.Name("sig_only"),
    ROOT.RooFit.Components("gauss"),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineColor(ROOT.kGreen + 2),
    ROOT.RooFit.LineWidth(2),
    # ROOT.RooFit.Precision(1e-5),
)

model.plotOn(
    frame,
    ROOT.RooFit.Name("bkg_only"),
    ROOT.RooFit.Components("bkg"),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineColor(ROOT.kRed + 1),
    ROOT.RooFit.LineWidth(2),
)

frame.GetXaxis().SetLabelSize(0)
frame.GetXaxis().SetTitleSize(0)
frame.GetYaxis().SetTitle("Counts")
frame.GetYaxis().SetTitleSize(0.055)
frame.GetYaxis().SetTitleOffset(0.95)
frame.GetYaxis().SetLabelSize(0.047)
frame.Draw()

# — legend —
leg = ROOT.TLegend(0.65, 0.62, 0.96, 0.88)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.032)
leg.AddEntry(frame.findObject("data"),     "Data",               "PE")
leg.AddEntry(frame.findObject("total"),    "Total fit",          "L")
leg.AddEntry(frame.findObject("sig_only"), "Signal (Gauss)",     "L")
leg.AddEntry(frame.findObject("bkg_only"), "Background (Cheb.)", "L")
leg.Draw()

# — fit parameter box (TLatex) —
lat = ROOT.TLatex()
lat.SetNDC()
lat.SetTextSize(0.040)
lat.SetTextAlign(12)

y0, dy = 0.85, 0.06  # starting y and line spacing
lat.DrawLatex(0.16, y0,        f"#mu  = {mean.getVal():.2f} #pm {mean.getError():.2f}  MeV")
lat.DrawLatex(0.16, y0 - dy,   f"#sigma = {sigma.getVal():.2f} #pm {sigma.getError():.2f}  MeV")
lat.DrawLatex(0.16, y0 - 2*dy, f"N_{{sig}} = {n_sig.getVal():.0f} #pm {n_sig.getError():.0f}")

pad_main.Update()

# ─────────────────────────────────────────────────────────────────────────────
#  Pull pad helper
# ─────────────────────────────────────────────────────────────────────────────
def draw_pull_pad(pad, pull_vals, bin_ctrs, xmin, xmax, nbins,
                  ylabel, fill_color, show_x_axis=False):
    """
    Draw a pull histogram on 'pad'.

    Parameters
    ----------
    pad         : TPad — target pad (already Draw()'n on the canvas)
    pull_vals   : list of float — pull value per bin
    bin_ctrs    : list of float — bin centre positions [MeV]
    xmin, xmax  : float — x range
    nbins       : int   — number of bins
    ylabel      : str   — y-axis label string (ROOT TLatex syntax)
    fill_color  : int   — ROOT colour index for the histogram fill
    show_x_axis : bool  — draw x-axis title & labels (only on the bottom pad)

    Returns
    -------
    (hp, ref_lines) — ROOT objects that must stay alive (Python GC protection)
    """
    pad.cd()

    # Build a TH1D filled with pull values
    hp = ROOT.TH1D(f"hp_{id(pad)}", "", nbins, xmin, xmax)
    for i, pv in enumerate(pull_vals, start=1):
        hp.SetBinContent(i, pv)

    # Symmetric y range: at least ±3σ, or ±1.4 × max(|pull|) if larger
    pmax   = max((abs(p) for p in pull_vals), default=3.0)
    yrange = max(pmax * 1.4, 3.0)
    hp.SetMaximum( yrange)
    hp.SetMinimum(-yrange)

    hp.SetLineColor(fill_color)
    hp.SetFillColorAlpha(fill_color, 0.42)
    hp.SetLineWidth(1)

    # Axis size: slightly larger than main pad (60%/40% ≈ 1.5× scaling)
    # Pull pad's 40% of canvas height requires proportionally larger labels

    hp.GetYaxis().SetTitle(ylabel)
    hp.GetYaxis().SetTitleSize(0.080)
    hp.GetYaxis().SetLabelSize(0.070)
    hp.GetYaxis().SetTitleOffset(0.50)
    hp.GetYaxis().SetNdivisions(504)

    if show_x_axis:
        hp.GetXaxis().SetTitle("m_{#eta}  [MeV]")
        hp.GetXaxis().SetTitleSize(0.080)
        hp.GetXaxis().SetLabelSize(0.070)
        hp.GetXaxis().SetTitleOffset(0.85)
    else:
        hp.GetXaxis().SetLabelSize(0)
        hp.GetXaxis().SetTitleSize(0)

    hp.Draw("HIST")

    # Reference lines:  solid black at 0,  dashed grey at ±2σ
    ref_lines = []
    for yval, color, style in [
        ( 0.0, ROOT.kBlack,    ROOT.kSolid),
        (-2.0, ROOT.kGray + 1, ROOT.kDashed),
        ( 2.0, ROOT.kGray + 1, ROOT.kDashed),
    ]:
        ln = ROOT.TLine(xmin, yval, xmax, yval)
        ln.SetLineColor(color)
        ln.SetLineStyle(style)
        ln.SetLineWidth(1)
        ln.Draw()
        ref_lines.append(ln)   # keep reference to prevent GC

    pad.Update()
    return hp, ref_lines


# ── Draw pull pad:  (d − f) / σ_d ────────────────────────────────────────────
pull_hist, ref_lines = draw_pull_pad(
    pad_pull, pulls, bin_ctrs, xmin, xmax, args.nbins,
    ylabel     = "(y_{data}#minus y_{fit}) / #sigma_{data}",
    fill_color = ROOT.kOrange + 7,
    show_x_axis = True,           # bottom pad: show x-axis with title
)

canvas.Update()


# ─── Save plot ────────────────────────────────────────────────────────────────
os.makedirs("out", exist_ok=True)
outpath = os.path.join("out", args.output)
canvas.SaveAs(outpath)
print(f"[done] Plot saved to: {outpath}")

rfile.Close()
