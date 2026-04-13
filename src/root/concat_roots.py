################################################################################
# Script to concatenate multiple ROOT files into a single ROOT file.           #
# Author: Michael Peters                                                       #
# To run: run src/concat_roots.py from anaroot/ directory.                     #
################################################################################

import ROOT
import os
import logging

# Set up logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    filename=args.outfile + '.log'
)
# Also print to stderr
logging.getLogger().addHandler(logging.StreamHandler())

# User-inputted path to nested search ROOT files to concatenate
print('Provide the base path containing the ROOT files to concatenate.')
print('The base path must contain magup and magdown subdirectories.')
print('Example: ntuple/mc_minbias_20260124/')
input_path = input('Path: ').strip()
if not input_path.endswith('/'): input_path += '/'
# Now find all ROOT files in this path and its subdirectories
infiles = []
# Search the magup and magdown subdirectories within the provided path
logging.info(f'Searching for ROOT files...')
for subdir in ['magup', 'magdown']:
    dir = os.path.join(input_path, subdir)
    if not os.path.exists(dir):
        logging.warning(f'Subdirectory {dir} does not exist. Skipping.')
        continue
    for root, _, files in os.walk(dir):
        for file in files:
            if file.endswith('.root'):
                infiles.append(os.path.join(root, file))
# Write outfile to base path
outfile = os.path.join(input_path, 'combined_files.root')

logging.info(f'Concatenating files:')
for infile in infiles:
    logging.info(f' - {infile}')
logging.info(f'Into single output file: {outfile}')

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
logging.info(f'Successfully created {outfile}.')
