import ROOT
from array import array

# TODO: Consider only filling histograms and returning array of histograms,
# instead of writing to file here. This would make the function more flexible.
def create_histograms(outfile, binwidths, arrays, names):
    '''Create histograms from arrays of data and save to ROOT TFile.
    Args:
        outfile (str): output ROOT file name
        binwidths (list<float>): list of bin widths for each histogram
        arrays (list<list<float>>): list of arrays of data for each histogram
        names (list<str>): list of histogram names
    '''
    # Check if binwidths, arrays, names have same length
    if not (len(binwidths) == len(arrays) == len(names)):
        raise ValueError("binwidths, arrays, and names must have the same length.")
    
    histfile = ROOT.TFile.Open(outfile, "RECREATE")
    # Switch scope to being in output TFile
    histfile.cd()

    # Loop over histogram variables, create histograms
    for i, (arr, name) in enumerate(zip(arrays, names)):
        binwidth = binwidths[i]
        xmin = round(min(arr) - binwidth)
        # Since the hist is shifted to left, need to add another binwidth to the
        # right side (max).
        xmax = round(max(arr) + 2*binwidth)
        nbins = int((xmax - xmin) / binwidth)
        bins = []
        # nbins + 1 to give margin on right edge
        for i in range(nbins+1):
            # Shifted by -0.5 to get proper binning
            bins.append(-0.5 + xmin + i * binwidth)
        # Args: source, title;x-axis label;y-axis label, bins, xmin, xmax
        # Create histogram with manual binning
        hist = ROOT.TH1D(name, name, nbins, array('d', bins))
        # Fill histogram
        for val in arr: hist.Fill(val)
        
        # Write histogram to TFile
        hist.Write()

    # Close histogram TFile
    histfile.Close()
    return
