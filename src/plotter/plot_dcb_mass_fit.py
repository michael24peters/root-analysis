#!/usr/bin/env python3
"""
dcb_mass_fit.py  –  Double Crystal Ball + 2nd-order Chebyshev fit
                    for the eta → μ+μ−γ reconstructed mass spectrum
──────────────────────────────────────────────────────────────────────────────
Reads a ROOT file whose TTree is always "tree" and whose eta candidate mass
branch is always "tag_m".  Fits a Double Crystal Ball (DCB) signal model on
top of a 2nd-order Chebyshev polynomial background using PyROOT / RooFit,
then produces:
  • a multi-panel PNG  — main fit  +  two pull panels
  • a JSON results file — parameters, errors, χ²/ndf, fit status

─────────────────────────────────────────────
  The Double Crystal Ball PDF
─────────────────────────────────────────────
A Gaussian core with independent power-law tails on each side of the peak:

          ╔  A_L · (B_L − t)^{−n_L}    if   t < −α_L       (left  tail)
  f(t) =  ║  exp(−t²/2)                if  −α_L ≤ t ≤ α_R  (Gaussian core)
          ╚  A_R · (B_R + t)^{−n_R}    if   t >  α_R        (right tail)

  where   t = (x − μ) / σ
  and the continuity + differentiability constants are
          A_{L,R} = (n / |α|)^n  · exp(−α²/2)
          B_{L,R} =  n / |α|    −  |α|

α_{L,R} > 0 : distance from the mean (in σ units) where the PDF switches
              from Gaussian to power-law.
n_{L,R} > 1 : power-law exponent; larger → harder / more abrupt tail falloff.

─────────────────────────────────────────────
  Usage
─────────────────────────────────────────────
  python dcb_mass_fit.py input.root [output.png] [options]

  # Asymmetric DCB (default)
  python dcb_mass_fit.py data.root fit.png --xmin 490 --xmax 610

  # Symmetric DCB  (α_R ≡ α_L, n_R ≡ n_L)
  python dcb_mass_fit.py data.root fit.png --symmetric

─────────────────────────────────────────────
  CLI reference
─────────────────────────────────────────────
Positional:
  rootfile                Input ROOT file (TTree="tree", branch="tag_m")
  output          (opt.)  Output PNG file                (default: dcb_fit.png)

Range / binning:
  --xmin  FLOAT           Lower mass bound [MeV]         (default: 480)
  --xmax  FLOAT           Upper mass bound [MeV]         (default: 620)
  --nbins INT             Number of histogram bins       (default: 80)

DCB initial values (starting points for the minimiser):
  --mean-init   FLOAT     Initial mean   [MeV]           (default: 548.0)
  --sigma-init  FLOAT     Initial σ      [MeV]           (default:  10.0)
  --alphaL-init FLOAT     Initial α_L threshold          (default:   1.5)
  --nL-init     FLOAT     Initial n_L  power             (default:   5.0)
  --alphaR-init FLOAT     Initial α_R threshold          (default:   1.5)
  --nR-init     FLOAT     Initial n_R  power             (default:   5.0)

DCB mode:
  --symmetric             Constrain α_R = α_L  and  n_R = n_L

Output:
  --output-json PATH      JSON results file  (default: dcb_fit_results.json)
  --dump-pulls            Save per-bin pull values to out/pulls.npy
"""

import sys
import os
import json
import argparse
import math
import numpy as np

# ─── CLI ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Double Crystal Ball + Chebyshev fit for an eta mass peak.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

# IO
parser.add_argument("rootfile",
                    help="Input ROOT file")
parser.add_argument("output", nargs="?", default="dcb_fit.png",
                    help="Output PNG file (default: dcb_fit.png)")

# Range / binning
g_range = parser.add_argument_group("range / binning")
g_range.add_argument("--xmin",  type=float, default=480.0,
                     help="Lower mass bound [MeV] (default: 480)")
g_range.add_argument("--xmax",  type=float, default=620.0,
                     help="Upper mass bound [MeV] (default: 620)")
g_range.add_argument("--nbins", type=int,   default=80,
                     help="Number of histogram bins (default: 80)")

# DCB initial values
g_dcb = parser.add_argument_group("DCB initial values")
g_dcb.add_argument("--mean-init",   type=float, default=548.0,
                   help="Initial mean [MeV] (default: 548.0)")
g_dcb.add_argument("--sigma-init",  type=float, default=10.0,
                   help="Initial sigma [MeV] (default: 10.0)")
g_dcb.add_argument("--alphaL-init", type=float, default=1.5,
                   help="Initial alpha_L (default: 1.5)")
g_dcb.add_argument("--nL-init",     type=float, default=5.0,
                   help="Initial n_L (default: 5.0)")
g_dcb.add_argument("--alphaR-init", type=float, default=1.5,
                   help="Initial alpha_R (default: 1.5)")
g_dcb.add_argument("--nR-init",     type=float, default=5.0,
                   help="Initial n_R (default: 5.0)")
g_dcb.add_argument("--symmetric",   action="store_true",
                   help="Constrain alpha_R = alpha_L and n_R = n_L (symmetric DCB)")

# Output
g_out = parser.add_argument_group("output")
g_out.add_argument("--output-json", default="dcb_fit_results.json",
                   help="JSON results file name/path written under out/ (default: dcb_fit_results.json)")
g_out.add_argument("--dump-pulls", action="store_true", default=False,
                   help="Save per-bin pull values to out/pulls.npy")

args = parser.parse_args()
xmin, xmax = args.xmin, args.xmax


# ─── ROOT import ──────────────────────────────────────────────────────────────
try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
except ImportError:
    sys.exit("PyROOT is not available. Please install ROOT with Python bindings.")


# ─── Load data ────────────────────────────────────────────────────────────────
if not os.path.exists(args.rootfile):
    sys.exit(f"File not found: {args.rootfile}")

rfile = ROOT.TFile.Open(args.rootfile, "READ")
if not rfile or rfile.IsZombie():
    sys.exit(f"Cannot open ROOT file: {args.rootfile}")

tree = rfile.Get("tree")
if not tree:
    keys = [k.GetName() for k in rfile.GetListOfKeys()]
    sys.exit(f"TTree 'tree' not found.  Available keys: {keys}")


# ─── Fill histogram from "tag_m" branch ───────────────────────────────────────
h = ROOT.TH1D("h_tagm", "tag_m;m_{#eta} (MeV);Counts", args.nbins, xmin, xmax)
h.Sumw2()  # sum of squares of weights for error bars
tree.Draw("tag_m>>h_tagm", "", "goff")   # "goff" = no graphics output
n_total = int(h.GetEntries())
print(f"[info] Entries in [{xmin}, {xmax}] MeV : {n_total}")
if n_total == 0:
    sys.exit("Histogram is empty — check branch name and mass range.")


# ─── RooFit observable ────────────────────────────────────────────────────────
m = ROOT.RooRealVar("m", "m_{#eta}  [MeV]", xmin, xmax)


# TODO: possibly share x0 but have independent sigma
# possibly should get fit params for each crystal ball fnc, then propagate to
# each of those; should look into HOW I'm doing this and figure out how to do
# phil's request
# ─── Signal model: Double Crystal Ball (RooCrystalBall) ───────────────────────
#
# RooCrystalBall(name, title, x, x0, sigma, alphaL, nL, alphaR, nR)
#
#   x      — the observable (reconstructed eta mass)
#   x0     — peak position  μ  (mean of the Gaussian core)
#   sigma  — width σ of the Gaussian core (shared by both sides)
#   alphaL — threshold α_L > 0:  left  tail starts at  x = μ − α_L·σ
#   nL     — power-law index n_L > 1 for the left  tail
#   alphaR — threshold α_R > 0:  right tail starts at  x = μ + α_R·σ
#   nR     — power-law index n_R > 1 for the right tail
#
# Physical interpretation for η→μ+μ−γ:
#   • The left  tail models energy losses from bremsstrahlung / FSR
#     and detector material interactions (particles radiating before detection).
#   • The right tail models mis-reconstructed events with extra energy deposits,
#     or overlapping tracks contributing to the invariant mass.
#   • σ captures the core detector mass resolution.
#
mean   = ROOT.RooRealVar("mean",   "mean",   args.mean_init,   510.0, 590.0)
sigma  = ROOT.RooRealVar("sigma",  "sigma",  args.sigma_init,    1.0,  50.0)
alphaL = ROOT.RooRealVar("alphaL", "alphaL", args.alphaL_init,  0.10,  10.0)
nL     = ROOT.RooRealVar("nL",     "nL",     args.nL_init,      1.01, 100.0)
alphaR = ROOT.RooRealVar("alphaR", "alphaR", args.alphaR_init,  0.10,  10.0)
nR     = ROOT.RooRealVar("nR",     "nR",     args.nR_init,      1.01, 100.0)

if args.symmetric:
    # Symmetric DCB: tie right-tail parameters to left-tail
    # RooFormulaVar evaluates "@0" = alphaL at every function call, so
    # alphaR_eff automatically tracks alphaL as the minimiser moves it.
    # Fitting only α_L and n_L reduces the free parameter count by 2.
    alphaR_eff = ROOT.RooFormulaVar(
        "alphaR_eff", "alphaR (= alphaL)", "@0", ROOT.RooArgList(alphaL)
    )
    nR_eff = ROOT.RooFormulaVar(
        "nR_eff", "nR (= nL)", "@0", ROOT.RooArgList(nL)
    )
    print("[info] Symmetric DCB: α_R ≡ α_L ,  n_R ≡ n_L")
else:
    # Asymmetric DCB: left and right tail parameters float independently
    alphaR_eff = alphaR
    nR_eff     = nR

# Build the DCB PDF using ROOT's native implementation
dcb = ROOT.RooCrystalBall(
    "dcb", "Double Crystal Ball",
    m, mean, sigma,
    alphaL, nL,  # left  tail
    alphaR_eff, nR_eff,  # right tail (may be tied in symmetric mode)
)


# ─── Background model: Chebyshev polynomial ───────────────────────────────────
#
# RooChebychev normalises x to x̃ ∈ [−1, 1] across the fit range, then builds:
#
#   B(x̃) = 1  +  c0·T₁(x̃)  +  c1·T₂(x̃)
#
# where the first three Chebyshev polynomials are:
#   T₀(x̃) = 1         (constant — absorbed into normalisation)
#   T₁(x̃) = x̃         (c0 controls the linear tilt / slope)
#   T₂(x̃) = 2x̃² − 1   (c1 controls quadratic curvature)
#
# This is a flexible, well-behaved background parametrisation that:
#   • is guaranteed positive over the fit range for small |c0|, |c1|
#   • has orthogonal basis functions → c0, c1 are nearly uncorrelated
#
c0 = ROOT.RooRealVar("c0", "c0",  0.0, -5.0, 5.0)
c1 = ROOT.RooRealVar("c1", "c1",  0.0, -5.0, 5.0)
bkg = ROOT.RooChebychev("bkg", "Chebyshev",
                         m, ROOT.RooArgList(c0, c1))


# ─── Extended combined model:  N_sig · DCB  +  N_bkg · Chebyshev ──────────────
#
# RooAddPdf in extended mode floats the total number of events as part of the
# likelihood, rather than treating it as fixed.
#
n_sig = ROOT.RooRealVar("n_sig", "signal yield",
                         n_total * 0.5, 0.0, n_total * 2.0)
n_bkg = ROOT.RooRealVar("n_bkg", "background yield",
                         n_total * 0.5, 0.0, n_total * 2.0)

model = ROOT.RooAddPdf(
    "model", "DCB + Chebyshev",
    ROOT.RooArgList(dcb, bkg),
    ROOT.RooArgList(n_sig, n_bkg),
)


# ─── Import histogram into RooFit and run the fit ─────────────────────────────
data = ROOT.RooDataHist("data", "data", ROOT.RooArgList(m), h)

print("[info] Running extended maximum-likelihood fit …")
fit_result = model.fitTo(
    data,
    ROOT.RooFit.Save(True),          # return RooFitResult (status, covariance)
    ROOT.RooFit.Extended(True),      # extended ML: floats N_total
    ROOT.RooFit.PrintLevel(-1),      # suppress MINUIT output to terminal
    ROOT.RooFit.Warnings(False),
)


# ─── χ²/ndf and per-bin fit predictions ───────────────────────────────────────
#
# Variables:
#  d_i: observed data count in bin i
#  f_i: model prediction for bin i (from fit)
#  chi2: how far data is from model prediction
#  Δᵢ: normalized residual for bin i, i.e., Δᵢ = (dᵢ − fᵢ) / √fᵢ
# 
# Strategy:
#   1. Plot data + model onto a dedicated RooFit frame.
#   2. Use frame.chiSquare(nFreeParams) → χ²/ndf  (Poisson normalisation,
#      i.e. denominator per bin = model prediction f_i).
#   3. Extract the fitted RooCurve from the frame; evaluate it at each bin
#      centre to get f_i for the pull calculation.

frame_diag = m.frame(ROOT.RooFit.Bins(args.nbins))
data.plotOn(frame_diag,
            ROOT.RooFit.Name("h_data"),
            ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))  # plot data on frame
model.plotOn(frame_diag, ROOT.RooFit.Name("fit_curve"))  # plot model on frame

# Number of freely floating parameters in this fit
n_float = fit_result.floatParsFinal().getSize()

# frame.chiSquare(curveName, histName, nPar) returns χ²/ndf
# where ndf = (number of occupied bins) − nPar(ams)
# This uses the Poisson definition for std: Δᵢ = (dᵢ − fᵢ) / √fᵢ
chi2_per_ndf = frame_diag.chiSquare("fit_curve", "h_data", n_float)

# Count occupied bins to recover the absolute χ²
n_occupied = sum(1 for i in range(1, args.nbins + 1)
                 if h.GetBinContent(i) > 0)
ndf      = n_occupied - n_float
chi2_abs = chi2_per_ndf * ndf


# ─── Per-bin values for the pull histograms ───────────────────────────────────
#
# Pull definition: data-based residual  (d − f) / σ_d 
#
# We normalize by the DATA uncertainty (from SumW2), not model prediction.
#
# Physics rationale (from advisor): Data is the ground truth with measured/real
# uncertainty; the fit is just one proposed explanation. So deviations should be
# scaled by what the data actually measured, not by what an imperfect model predicts.
#
# For unweighted events: σ_d = sqrt(d), so pull = (d − f) / sqrt(d).
#
# Note: RooFit's frame.chiSquare() uses Poisson-style (d − f) / sqrt(f) for its
# global χ²/ndf summary, but our per-bin pull plot uses data-based (d − f) / σ_d,
# which is the appropriate residual for visual diagnostics on unweighted data.
#
# Distributed as N(0,1) if model is correct and bin statistics are sufficient.

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

# Save per-bin pull values to a numpy file
if args.dump_pulls:
    np.save(os.path.join("out", "pulls.npy"), np.array(pulls))
    print("[info] Per-bin pull values saved to: out/pulls.npy")


# ─── Print fit summary to terminal ────────────────────────────────────────────
STATUS_CODES = {
    0: "Converged successfully",
    1: "Covariance forced positive-definite (approximate errors)",
    2: "Hessian invalid (errors unreliable)",
    3: "EDM above tolerance (fit may not be at minimum)",
    4: "Call limit reached before convergence",
}
status_code = fit_result.status()

# Resolve right-tail parameter values for printing
# (in symmetric mode these equal the left-tail values)
aR_val = alphaL.getVal()   if args.symmetric else alphaR.getVal()
aR_err = alphaL.getError() if args.symmetric else alphaR.getError()
nR_val = nL.getVal()       if args.symmetric else nR.getVal()
nR_err = nL.getError()     if args.symmetric else nR.getError()

W = 62  # table width
print()
print("═" * W)
print(f"{'DCB + CHEBYSHEV  FIT  RESULT':^{W}}")
print("═" * W)
print(f"  {'Parameter':<24} {'Value':>10}    {'Error':<10}")
print("─" * W)
print(f"  {'mean   [MeV]':<24} {mean.getVal():>10.4f}  ±  {mean.getError():<10.4f}")
print(f"  {'sigma  [MeV]':<24} {sigma.getVal():>10.4f}  ±  {sigma.getError():<10.4f}")
print(f"  {'α_L':<24} {alphaL.getVal():>10.4f}  ±  {alphaL.getError():<10.4f}")
print(f"  {'n_L':<24} {nL.getVal():>10.4f}  ±  {nL.getError():<10.4f}")
if args.symmetric:
    print(f"  {'α_R  [tied to α_L]':<24} {aR_val:>10.4f}     (no independent error)")
    print(f"  {'n_R  [tied to n_L]':<24} {nR_val:>10.4f}     (no independent error)")
else:
    print(f"  {'α_R':<24} {alphaR.getVal():>10.4f}  ±  {alphaR.getError():<10.4f}")
    print(f"  {'n_R':<24} {nR.getVal():>10.4f}  ±  {nR.getError():<10.4f}")
print(f"  {'c0  (Cheb. slope)':<24} {c0.getVal():>10.4f}  ±  {c0.getError():<10.4f}")
print(f"  {'c1  (Cheb. curv.)':<24} {c1.getVal():>10.4f}  ±  {c1.getError():<10.4f}")
print(f"  {'N_sig':<24} {n_sig.getVal():>10.1f}  ±  {n_sig.getError():<10.1f}")
print(f"  {'N_bkg':<24} {n_bkg.getVal():>10.1f}  ±  {n_bkg.getError():<10.1f}")
print("─" * W)
print(f"  χ²          = {chi2_abs:.2f}")
print(f"  ndf         = {ndf}  "
      f"(occupied bins = {n_occupied},  free params = {n_float})")
print(f"  χ²/ndf      = {chi2_per_ndf:.4f}")
print(f"  fit status  = {status_code}  →  {STATUS_CODES.get(status_code, 'unknown')}")
print("═" * W)
print()


# ─── Canvas layout ────────────────────────────────────────────────────────────

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

canvas = ROOT.TCanvas("c", "DCB mass fit", 900, 950)

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


# ─── Main fit pad ─────────────────────────────────────────────────────────────
pad_main.cd()

frame = m.frame(ROOT.RooFit.Bins(args.nbins))

# Data points with SumW2 error bars
data.plotOn(
    frame,
    ROOT.RooFit.Name("data"),
    ROOT.RooFit.MarkerStyle(20),
    ROOT.RooFit.MarkerSize(0.85),
    ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
)

# Total model (DCB + Chebyshev), solid blue
model.plotOn(
    frame,
    ROOT.RooFit.Name("total"),
    ROOT.RooFit.LineColor(ROOT.kBlue + 1),
    ROOT.RooFit.LineWidth(2),
)

# Signal component only (DCB), dashed green
model.plotOn(
    frame,
    ROOT.RooFit.Name("sig_only"),
    ROOT.RooFit.Components("dcb"),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineColor(ROOT.kGreen + 2),
    ROOT.RooFit.LineWidth(2),
)

# Background component only (Chebyshev), dashed red
model.plotOn(
    frame,
    ROOT.RooFit.Name("bkg_only"),
    ROOT.RooFit.Components("bkg"),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineColor(ROOT.kRed + 1),
    ROOT.RooFit.LineWidth(2),
)

# Hide x-axis on the main pad (shared visual axis is on the bottom pull pad)
frame.GetXaxis().SetLabelSize(0)
frame.GetXaxis().SetTitleSize(0)
frame.GetYaxis().SetTitle("Counts")
frame.GetYaxis().SetTitleSize(0.055)
frame.GetYaxis().SetTitleOffset(0.95)
frame.GetYaxis().SetLabelSize(0.047)
frame.Draw()

# --- Legend 1: curve / component identification ---
leg1 = ROOT.TLegend(0.65, 0.62, 0.96, 0.88)
leg1.SetBorderSize(0)
leg1.SetFillStyle(0)
# leg1.SetTextFont(42)
leg1.SetTextSize(0.032)
leg1.AddEntry(frame.findObject("data"),     "Data",                     "PE")
leg1.AddEntry(frame.findObject("total"),    "Total fit",                 "L")
leg1.AddEntry(frame.findObject("sig_only"), "Signal (DCB)",              "L")
leg1.AddEntry(frame.findObject("bkg_only"), "Background (Cheb.)",        "L")
leg1.Draw()

# --- Legend 2: fit parameters (TLatex) ---
#
# Drawn on the main pad in NDC coordinates.
# Each row corresponds to one fit parameter or summary statistic.
#
lat = ROOT.TLatex()
lat.SetNDC()
lat.SetTextFont(42)
lat.SetTextSize(0.040)
lat.SetTextAlign(12)    # left-align horizontally, centre vertically

x0_lat = 0.16          # left edge of the annotation block
y0_lat = 0.85          # top y position
dy_lat = 0.06         # vertical spacing between rows

param_rows = [
    f"#mu = {mean.getVal():.2f} #pm {mean.getError():.2f}  MeV",
    f"#sigma = {sigma.getVal():.2f} #pm {sigma.getError():.2f}  MeV",
    f"#alpha_{{L}} = {alphaL.getVal():.3f} #pm {alphaL.getError():.3f}",
    f"n_{{L}} = {nL.getVal():.2f} #pm {nL.getError():.2f}",
    (f"#alpha_{{R}} = {aR_val:.3f} #pm {aR_err:.3f}"
     + ("  [sym]" if args.symmetric else "")),
    (f"n_{{R}} = {nR_val:.2f} #pm {nR_err:.2f}"
     + ("  [sym]" if args.symmetric else "")),
    f"N_{{sig}} = {n_sig.getVal():.0f} #pm {n_sig.getError():.0f}",
    f"#chi^{{2}}/ndf = {chi2_per_ndf:.3f}",
]
for k, row in enumerate(param_rows):
    lat.DrawLatex(x0_lat, y0_lat - k * dy_lat, row)

pad_main.Update()


# ─── Pull pad helper ──────────────────────────────────────────────────────────
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


# ─── Save JSON results ────────────────────────────────────────────────────────
results = {
    "meta": {
        "input_file":    args.rootfile,
        "fit_model":     "Double Crystal Ball + Chebyshev",
        "symmetric_dcb": args.symmetric,
        "mass_range":    [xmin, xmax],
        "nbins":         args.nbins,
        "n_entries":     n_total,
    },
    "fit_status": {
        "code":    status_code,
        "meaning": STATUS_CODES.get(status_code, "unknown"),
    },
    "chi2": {
        "chi2":            round(chi2_abs,     4),
        "ndf":             ndf,
        "chi2_per_ndf":    round(chi2_per_ndf, 6),
        "n_occupied_bins": n_occupied,
        "n_free_params":   n_float,
        "definition":      "Poisson: sum((d_i - f_i)^2 / f_i) over occupied bins",
    },
    "parameters": {
        "mean": {
            "value": round(mean.getVal(),   6),
            "error": round(mean.getError(), 6),
            "unit":  "MeV",
        },
        "sigma": {
            "value": round(sigma.getVal(),   6),
            "error": round(sigma.getError(), 6),
            "unit":  "MeV",
        },
        "alphaL": {
            "value": round(alphaL.getVal(),   6),
            "error": round(alphaL.getError(), 6),
        },
        "nL": {
            "value": round(nL.getVal(),   6),
            "error": round(nL.getError(), 6),
        },
        "alphaR": {
            "value": round(aR_val, 6),
            "error": round(aR_err, 6),
            "note":  "tied to alphaL (symmetric mode)" if args.symmetric else "free",
        },
        "nR": {
            "value": round(nR_val, 6),
            "error": round(nR_err, 6),
            "note":  "tied to nL (symmetric mode)" if args.symmetric else "free",
        },
        "c0": {
            "value": round(c0.getVal(),   6),
            "error": round(c0.getError(), 6),
            "note":  "Chebyshev linear coefficient (slope)",
        },
        "c1": {
            "value": round(c1.getVal(),   6),
            "error": round(c1.getError(), 6),
            "note":  "Chebyshev quadratic coefficient (curvature)",
        },
        "n_sig": {
            "value": round(n_sig.getVal(),   2),
            "error": round(n_sig.getError(), 2),
        },
        "n_bkg": {
            "value": round(n_bkg.getVal(),   2),
            "error": round(n_bkg.getError(), 2),
        },
    },
}

# Write json output, make sure it goes to out/
os.makedirs("out", exist_ok=True)
json_outpath = args.output_json
if not os.path.isabs(json_outpath):
    norm_json = os.path.normpath(json_outpath)
    out_prefix = "out" + os.sep
    if norm_json != "out" and not norm_json.startswith(out_prefix):
        json_outpath = os.path.join("out", json_outpath)

with open(json_outpath, "w") as jf:
    json.dump(results, jf, indent=2)
print(f"[done] JSON results saved to: {json_outpath}")

# Close ROOT file
rfile.Close()