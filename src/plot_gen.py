###############################################################################
# Step 4 of 4                                                                 #
# Script to make mass plots from histograms and save to png files             #
# Author: Michael Peters                                                      #
###############################################################################

# consider:
# pull plot: see example, either do difference / bin error OR ratio

import ROOT
import sys

# Optional command line argument: include stats box (none by default)
include_stats = False
sig_file = False
if len(sys.argv) > 1:
    if 'stats' in sys.argv[1:]:
        include_stats = True
    if 'sig' in sys.argv[1:]:
        sig_file = True

if sig_file:
    infile = 'hist/sig_hist_gen.root'
    fileheader = 'figs/sig/'
else: 
    infile = 'hist/hist_gen.root'
    fileheader = 'figs/minbias/'

print(f'Reading from {infile} and writing to {fileheader}*.png')

tfile = ROOT.TFile.Open(infile, 'READ')

# Get histograms from TFile
hpid = tfile.Get('mc_pid')
hm = tfile.Get('mc_m')
hp = tfile.Get('mc_p')
hpt = tfile.Get('mc_pt')
hpz = tfile.Get('mc_pz')

# Configure histogram
for hist in (hpid, hm, hp, hpt, hpz):
    # hist.Sumw2()  # statistical uncertainties by sum of weights squared
    # Keep histogram in memory (not remove when files are closed)
    hist.SetDirectory(0)
    hist.SetStats(0)
    hist.GetYaxis().SetTitle("Events")
    hist.GetXaxis().SetTitle("Momentum [MeV/c]")
    hist.SetFillStyle(0)  # no fill style
    hist.SetLineColor(ROOT.kBlue)

# Create canvas
canvas = ROOT.TCanvas('canvas')
canvas.cd()

# Draw each histogram separately
hists = [hpid, hp,
         hm,
         hpt, hpz]
titles = ['mc pid', 'mc momentum',
          'mc mass',
          'mc transverse momentum', 'mc pz']
names = ['mc_pid', 'mc_p', 
         'mc_m', 
         'mc_pt', 'mc_pz']

for hist, name, title in zip(hists, names, titles):
    hist.SetTitle(title)
    # Include special formatting for specific histograms to look nice
    if name == 'mc_m': hist.GetXaxis().SetTitle("Mass [MeV]")
    elif name in ['mc_pid']: hist.GetXaxis().SetTitle("Particle ID")
    if name == 'mc_pid': hist.GetXaxis().SetRangeUser(-14.5, 223.5)
    if name == 'mc_m': hist.GetXaxis().SetRangeUser(546.5, 549.5)
    hist.Draw('h')
    if include_stats:
        hist.SetStats(1)
        canvas.Update()  # Ensure stats box is created
        if name == 'mc_pid':
            # Shift to top middle third of plot
            stat = hist.GetListOfFunctions().FindObject("stats")
            stat.SetX1NDC(0.28)
            stat.SetX2NDC(0.48)
            stat.SetY1NDC(0.70)
            stat.SetY2NDC(0.85)
        else: hist.SetStats(1)
    canvas.Print(f'{fileheader}{name}.png')

# Clear the canvas
canvas.Clear()

print(f'Done: wrote plots to {fileheader}gen*.png')
