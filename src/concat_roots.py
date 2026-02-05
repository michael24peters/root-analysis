################################################################################
# Script to concatenate multiple ROOT files into a single ROOT file.           #
# Author: Michael Peters                                                       #
################################################################################

import ROOT
import os

# User-inputted path to nested search ROOT files to concatenate
print('Provide the base path containing the ROOT files to concatenate.')
print('The base path must contain magup and magdown subdirectories.')
print('Example: ntuple/mc_minbias_20260124/')
input_path = input('Path: ').strip()
if not input_path.endswith('/'): input_path += '/'
# Now find all ROOT files in this path and its subdirectories
infiles = []
# Search the magup and magdown subdirectories within the provided path
for subdir in ['magup', 'magdown']:
    dir = os.path.join(input_path, subdir)
    if not os.path.exists(dir):
        print(f'Warning: Subdirectory {dir} does not exist. Skipping.')
        continue
    for root, _, files in os.walk(dir):
        for file in files:
            if file.endswith('.root'):
                infiles.append(os.path.join(root, file))
# Write outfile to base path
outfile = os.path.join(input_path, 'combined_files.root')

print(f'Concatenating files:')
for infile in infiles:
    print(f' - {infile}')
print(f'Into single output file: {outfile}')

# Create a TFileMerger to merge the files
merger = ROOT.TFileMerger()
merger.OutputFile(outfile)
# Open each input file and add it to the list
for infile in infiles:
    tfile = ROOT.TFile.Open(infile, 'READ')
    # Add TFile to merger
    merger.AddFile(tfile)
# Use ROOT's TFileMerger to merge the files
merger.Merge()
print(f'Successfully created {outfile}.')
