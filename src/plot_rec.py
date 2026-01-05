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
    infile = 'hist/sig_hist_rec.root'
    fileheader = 'figs/sig/rec'
else:
    infile = 'hist/hist_rec.root'
    fileheader = 'figs/minbias/rec'

print(f'Reading from {infile} and writing to {fileheader}_*.png')

tfile = ROOT.TFile.Open(infile, 'READ')

# Get histograms from TFile
htpid = tfile.Get('tag_pid')
htm = tfile.Get('tag_m')
htp = tfile.Get('tag_p')
htpt = tfile.Get('tag_pt')
htpz = tfile.Get('tag_pz')
hppid = tfile.Get('prt_pid')
hpp = tfile.Get('prt_p')
hppz = tfile.Get('prt_pz')
hppt = tfile.Get('prt_pt')

# Configure histogram
for hist in (htpid, htm, htp, htpt, htpz, hppid, hpp, hppz, hppt):
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
hists = [htpid, htp, 
         htpt, htpz,
         htm,
         hppid, hpp, 
         hppt, hppz]
titles = ['reconstructed tag pid', 'reconstructed tag momentum', 
          'reconstructed tag transverse momentum', 'reconstructed tag pz',
          'reconstructed tag mass',
          'reconstructed daughter pid', 'reconstructed daughter momentum', 
          'reconstructed daughter transverse momentum', 'reconstructed daughter pz']
names = ['tag_pid', 'tag_p', 
         'tag_pt', 'tag_pz', 
         'tag_m',
         'prt_pid', 'prt_p', 
         'prt_pt', 'prt_pz']

for hist, name, title in zip(hists, names, titles):
    hist.SetTitle(title)
    if name == 'tag_m': hist.GetXaxis().SetTitle("Mass [MeV]")  # special case
    elif name in ['prt_pid', 'tag_pid']: hist.GetXaxis().SetTitle("Particle ID")  # special case
    hist.Draw('h')
    if include_stats:
        hist.SetStats(1)
        canvas.Update()  # Ensure stats box is created
        if name == 'prt_pid':
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

print(f'Done: wrote plots to {fileheader}*.png')
