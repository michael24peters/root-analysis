################################################################################
# Create input and output files and directories for ROOT hist/plot scripts.    #
# Author: Michael Peters                                                       #
################################################################################

import argparse
import os

def parse_args(subdir, name):
    '''
    Parse command line arguments for input and output files. Creates output
    directory if it doesn't exist.
    
    Args:
        subdir (str): subdirectory name (e.g. 'hist', 'plot'). Also effects file
                      name (e.g. 'hist_gen.root', 'plot_gen.png', etc.)
        name (str): name to include in output file name (e.g. 'gen', 'mass',
                    'rec')
    
    Returns:
        infile (str): input ROOT file name
        outfile (str): output ROOT file name
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--infile',
        help='Input ROOT file'
    )
    parser.add_argument(
        '-s', '--signal',
        action='store_true',
        help='Use signal file'
    )
    parser.add_argument(
        '-d', '--decay',
        help='Decay mode (e.g. eta2mumu, eta2mumugamma, etc.)'
    )
    args = parser.parse_args()
    
    # Check for valid decay argument
    valid_decays = ['eta2mumu', 'eta2mumugamma']
    if args.decay not in valid_decays:
        raise ValueError(f'Invalid decay mode. Must be one of: {valid_decays}')

    infile = args.infile
    decay = args.decay

    # Create directory if it doesn't exist
    if args.decay: outdir = f'out/{subdir}/{decay}/'
    else: raise ValueError('Decay mode must be specified with -d flag.')
    if args.signal: outdir += 'signal/'
    else: outdir += 'minbias/'
    os.makedirs(outdir, exist_ok=True)
    
    # Define outfile
    outfile = f'{outdir}{subdir}_{name}.root'
    
    return infile, outfile, decay