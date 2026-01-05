###############################################################################
# Step 4 of 4                                                                 #
# Script to make mass plots from histograms and save to png files             #
# Author: Michael Peters                                                      #
###############################################################################
# TODO: consider pull plot. See tutorial; either do difference / bin error OR 
# ratio.

import ROOT
import sys

# Optional command line arguments: include legend, include stats box.
# By default, no legend or stats box.
include_legend = False
include_stats = False
sig_file = False
if len(sys.argv) > 1:
    if 'legend' in sys.argv[1:]:
        include_legend = True
    if 'stats' in sys.argv[1:]:
        include_stats = True
    if 'sig' in sys.argv[1:]:
        sig_file = True

if sig_file:
    infile = 'hist/sig_hist_m.root'
    fileheader = 'figs/sig/tag_m'
else:
    infile = 'hist/hist_m.root'
    fileheader = 'figs/minbias/tag_m'

print(f'Reading from {infile} and writing to {fileheader}_*.png')

tfile = ROOT.TFile.Open(infile, 'READ')

# Get histograms from TFile
hsig = tfile.Get('sig')
hbkg = tfile.Get('bkg')

# Configure histogram
for hist in (hbkg, hsig):
    # hist.Sumw2()  # statistical uncertainties by sum of weights squared
    # Keep histogram in memory (not remove when files are closed)
    hist.SetDirectory(0)
    hist.SetStats(0)
    hist.GetXaxis().SetTitle("Mass [MeV]")
    hist.GetYaxis().SetTitle("Events")

# Histogram styling
# Signal styling - black points with error bars
hsig.SetMarkerStyle(20)  # filled circle
hsig.SetMarkerSize(.5)
hsig.SetLineWidth(1)
hsig.SetLineColor(ROOT.kBlack)
hsig.SetMarkerColor(ROOT.kBlack)

# Background styling - gray fill
hbkg.SetFillStyle(1001)
hbkg.SetFillColor(ROOT.kGray+1)
hbkg.SetLineColor(ROOT.kGray+1)

# Create canvas
canvas = ROOT.TCanvas('canvas')
canvas.cd()


# =============================================================================

# Configure legend
def draw_legend(hsig=None, hbkg=None):
    print("Drawing legend...")  # debug
    # Place legend in top-right corner partially inside plot area
    # (xlow, ylow, xup, yup)
    leg = ROOT.TLegend(0.79, 0.78, 0.98, 0.88)
    leg.SetTextSize(0.03)
    leg.SetBorderSize(0)  # 0: remove border
    # leg.SetFillStyle(0)  # 0: transparent background
    if hsig: leg.AddEntry(hsig, "Signal", "l")  # l: line
    if hbkg: leg.AddEntry(hbkg, "Background", "f")  # f: fill
    leg.Draw()
    return leg

# =============================================================================
# Draw each histogram separately

# signal
hsig.SetTitle('signal tag mass')
hsig.Draw('pe1x0')  # p e 1 x0 : points; error bars; error bar lines; no x bars
if include_stats: hsig.SetStats(1)
canvas.Print(f'{fileheader}_sig.png')

# background
hbkg.SetTitle('background tag mass')
hbkg.Draw('h')  # h: histogram
if include_stats: hbkg.SetStats(1)
canvas.Print(f'{fileheader}_bkg.png')

# Clear custom configs
for hist in (hbkg, hsig): 
    hist.SetStats(0)
    hist.SetTitle('')
# Clear the canvas
canvas.Clear()

# =============================================================================
# Draw all histograms on one canvas

# Title on combined plot
hbkg.SetTitle('combined tag mass')

# Find the largest y value among all histograms
def get_max_with_error(hist, err=False, margin=1):
    max_val = -float('inf')
    for i in range(1, hist.GetNbinsX()+1):
        val = hist.GetBinContent(i)
        bin_err = hist.GetBinError(i) if err else 0
        if val + bin_err > max_val:
            max_val = val + bin_err
    # Add a small margin to avoid clipping error bars
    return max_val + margin

# Find maximum y value considering error bars for signal
ymax = max(get_max_with_error(hsig, err=True), get_max_with_error(hbkg))

# Set the maximum for all histograms to ensure consistent y range
for hist in (hsig, hbkg):
    hist.SetMaximum(ymax)
    hist.SetMinimum(0)

hbkg.Draw('h')
hsig.Draw('pe1x0, same')
if include_legend: leg = draw_legend(hsig=hsig, hbkg=hbkg)
# Save canvas to file
canvas.Print(f'{fileheader}_combined.png')

# Clear the canvas
# Reset formatting
for hist in (hsig, hbkg):
    hist.SetMaximum()
    hist.SetMinimum()
canvas.Clear()

# =============================================================================
# Draw histograms on split-panel canvas

# Format histograms
# Top plot should have y axis ticks to only be integers, y axis label, no x axis label, 

# Title on top plot
hsig.SetTitle('signal vs background tag mass')
ROOT.gStyle.SetTitleFontSize(0.07)
# No title on bottom plot
hbkg.SetTitle('')

# Top pad (70% of canvas)
# Set axes labels and value sizes
pad1_label_size = 0.04
pad1_title_size = 0.04
hsig.GetXaxis().SetTitleSize(pad1_title_size)
hsig.GetYaxis().SetTitleSize(pad1_title_size)
hsig.GetXaxis().SetLabelSize(pad1_label_size)
hsig.GetYaxis().SetLabelSize(pad1_label_size)

# Bottom pad (30% of canvas)
# Set axes labels and value sizes
pad2_label_size = 0.12
pad2_title_size = 0.12
hbkg.GetXaxis().SetTitleSize(pad2_title_size)
hbkg.GetYaxis().SetTitleSize(pad2_title_size)
hbkg.GetXaxis().SetLabelSize(pad2_label_size)
hbkg.GetYaxis().SetLabelSize(pad2_label_size)

# No x axis label on top plot
hsig.GetXaxis().SetTitle('')
# No y axis label on bottom plot
hbkg.GetYaxis().SetTitle('')

# Set y axis ticks to only be integers
hsig.GetYaxis().SetNdivisions(5)
hbkg.GetYaxis().SetNdivisions(4)

# Set up first pad - takes up top 70% of canvas
# Args: name, title, xlow, ylow, xup, yup
pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.3, 1, 1)
pad1.Draw()
pad1.cd()
# Top margin larger to accommodate title
pad1.SetTopMargin(0.15)
# Bottom margin 0 to connect to bottom pad
pad1.SetBottomMargin(0)
# Draw signal to top pad
hsig.Draw('pe1x0')
if include_legend: leg = draw_legend(hsig=hsig, hbkg=hbkg)

# Remove 0 value label on y-axis top plot 
pad1.Update()
hsig.GetYaxis().ChangeLabel(1, -1, -1, -1, -1, -1, " ")

# Return to canvas level
canvas.cd()

# Set up second pad - takes up bottom 30% of canvas
# Args: name, title, xlow, ylow, xup, yup
pad2 = ROOT.TPad('pad2', 'pad2', 0, 0.05, 1, 0.3)
pad2.Draw()
pad2.cd()
# Formatting - top margin 0 to connect to top pad
pad2.SetTopMargin(0)
# Bottom margin larger to accommodate x-axis labels
pad2.SetBottomMargin(0.25)
# Draw background to bottom pad
hbkg.Draw('h')
canvas.Print(f'{fileheader}_split.png')

# Clear the canvas
canvas.Clear()

# Reset formatting
hbkg.SetFillColor(ROOT.kGray+1)
hbkg.SetLineColor(ROOT.kGray+1)
# Clear the canvas
canvas.Clear()

print(f'Done: wrote plots to {fileheader}_*.png')
