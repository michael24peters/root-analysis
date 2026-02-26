"""
eta_peak_fit.py
───────────────
Reads a ROOT file, extracts the `tag_m` branch (reconstructed eta candidate
masses), fits a peak model (Gaussian signal + Chebyshev polynomial background)
using PyROOT / RooFit, and saves the result as a PNG.

Usage
─────
    python eta_peak_fit.py input.root [tree_name] [output.png]

Positional arguments (all optional after the ROOT file):
    tree_name  : name of the TTree containing `tag_m`  (default: "tree")
    output     : output image file name                (default: "eta_fit.png")

Optional arguments:
    --xmin      : lower mass range for fit/plot [MeV] (default: 480)
    --xmax      : upper mass range for fit/plot [MeV] (default: 620)
    --nbins     : number of histogram bins            (default: 80)
    TODO: implement sideband constrined fit (args are not yet implemented)
    --sideband-lo LO1 LO2 : lower sideband region [LO1, LO2] MeV
    --sideband-hi HI1 HI2 : upper sideband region [HI1, HI2] MeV
"""

import sys
import os
import argparse

# ── argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Fit the eta peak in the tag_m branch of a ROOT file."
)
parser.add_argument("rootfile",           help="Path to the input ROOT file")
parser.add_argument("tree",  nargs="?", default="tree",
                    help="TTree name inside the ROOT file (default: tree)")
parser.add_argument("output", nargs="?", default="eta_fit.png",
                    help="Output PNG file (default: eta_fit.png)")
parser.add_argument("--xmin",  type=float, default=480.0,
                    help="Lower mass range for fit / plot [MeV] (default: 480)")
parser.add_argument("--xmax",  type=float, default=620.0,
                    help="Upper mass range for fit / plot [MeV] (default: 620)")
parser.add_argument("--nbins", type=int,   default=80,
                    help="Number of histogram bins (default: 80)")
parser.add_argument("--sideband-lo", type=float, nargs=2, metavar=("LO1","LO2"),
                    default=[480.0, 520.0],
                    help="Lower sideband region excluded from signal fit")
parser.add_argument("--sideband-hi", type=float, nargs=2, metavar=("HI1","HI2"),
                    default=[580.0, 620.0],
                    help="Upper sideband region excluded from signal fit")
args = parser.parse_args()

# ── ROOT imports ──────────────────────────────────────────────────────────────
try:
    import ROOT
    ROOT.gROOT.SetBatch(True)          # suppress interactive canvas pop-ups
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
except ImportError:
    sys.exit("PyROOT is not available. Please install ROOT with Python bindings.")

# ── load data ─────────────────────────────────────────────────────────────────
if not os.path.exists(args.rootfile):
    sys.exit(f"File not found: {args.rootfile}")

rfile = ROOT.TFile.Open(args.rootfile, "READ")
if not rfile or rfile.IsZombie():
    sys.exit(f"Could not open ROOT file: {args.rootfile}")

tree = rfile.Get(args.tree)
if not tree:
    sys.exit(
        f"TTree '{args.tree}' not found in {args.rootfile}.\n"
        f"Available keys: {[k.GetName() for k in rfile.GetListOfKeys()]}"
    )

# ── fill histogram from branch ────────────────────────────────────────────────
xmin, xmax = args.xmin, args.xmax
h = ROOT.TH1D("h_tagm", "tag_m;m_{#eta} (MeV);Counts",
              args.nbins, xmin, xmax)
h.Sumw2()
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
'''
Chebyshev polynomial used for creating a smooth background PDF.  Mass is
normalized to [-1, 1] range. The first three Chebyshev polynomials are:

T0(x) = 1, T1(x) = x, T2(x) = 2x^2 - 1.

The background PDF is then modeled by a linear combination of Chebyshev
polynomials:

B(x) = 1 + c0 * T1(x) + c1 * T2(x), 

where c0, c1 are fit params that build the background PDF, x is the normalized
mass variable. c0 controls slope (linear) and c1 controls curvature (quadratic)
'''
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

# ── build plot ────────────────────────────────────────────────────────────────
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

c = ROOT.TCanvas("c", "eta peak fit", 900, 700)
c.SetLeftMargin(0.13)
c.SetBottomMargin(0.13)
c.SetRightMargin(0.05)

frame = m.frame(ROOT.RooFit.Title("eta Peak Fit"))

data.plotOn(
    frame,
    ROOT.RooFit.Name("data"),
    ROOT.RooFit.MarkerStyle(20),
    ROOT.RooFit.MarkerSize(0.8),
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
    ROOT.RooFit.Name("bkg_only"),
    ROOT.RooFit.Components("bkg"),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineColor(ROOT.kRed + 1),
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

frame.GetXaxis().SetTitle("m_{#eta}  (MeV)")
frame.GetYaxis().SetTitle("Counts")
frame.GetXaxis().SetTitleSize(0.035)
frame.GetYaxis().SetTitleSize(0.035)
frame.GetXaxis().SetLabelSize(0.03)
frame.GetYaxis().SetLabelSize(0.03)

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
lat.SetTextSize(0.033)
lat.SetTextAlign(12)

y0, dy = 0.85, 0.06   # starting y and line spacing
lat.DrawLatex(0.16, y0,        f"#mu  = {mean.getVal():.2f} #pm {mean.getError():.2f}  MeV")
lat.DrawLatex(0.16, y0 - dy,   f"#sigma = {sigma.getVal():.2f} #pm {sigma.getError():.2f}  MeV")
lat.DrawLatex(0.16, y0 - 2*dy, f"N_{{sig}} = {n_sig.getVal():.0f} #pm {n_sig.getError():.0f}")

# — vertical line at fitted mean —
# vline = ROOT.TLine(mean.getVal(), frame.GetMinimum(), mean.getVal(), frame.GetMaximum())
# vline.SetLineColor(ROOT.kBlue + 1)
# vline.SetLineStyle(ROOT.kDashed)
# vline.SetLineWidth(1)
# vline.Draw()

c.Update()
c.SaveAs('out/' + args.output)
print(f"\n[done] Plot saved to: {'out/' + args.output}")

rfile.Close()