################################################################################
# Create input and output files and directories for ROOT hist/plot scripts.    #
# Author: Michael Peters                                                       #
################################################################################

import argparse
import os

def parse_hist_args(name):
    '''
    Parse command line arguments for histogram scripts. Creates output
    directory if it doesn't exist.

    Args:
        name (str): name to include in output file name (e.g. 'gen', 'mass',
                    'rec')

    Returns:
        infile (str): input ROOT file name
        outfile (str): output ROOT file name
        decay (str): decay mode
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
    outdir = f'out/hist/{decay}/'
    if args.signal: outdir += 'signal/'
    else: outdir += 'minbias/'
    os.makedirs(outdir, exist_ok=True)

    outfile = f'{outdir}hist_{name}.root'

    return infile, outfile, decay


def parse_plot_args(name):
    '''
    Parse command line arguments for plot scripts. Creates output directory
    if it doesn't exist.

    Args:
        name (str): name to include in output file prefix (e.g. 'gen', 'mass',
                    'rec')

    Returns:
        infile (str): input ROOT file name
        fileheader (str): output file prefix for PNG files
        decay (str): decay mode
        include_stats (bool): whether to include stats box
        include_legend (bool): whether to include legend
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
    parser.add_argument(
        '--stats',
        action='store_true',
        help='Include stats box'
    )
    parser.add_argument(
        '--legend',
        action='store_true',
        help='Include legend'
    )
    args = parser.parse_args()

    # Check for valid decay argument
    valid_decays = ['eta2mumu', 'eta2mumugamma']
    if args.decay not in valid_decays:
        raise ValueError(f'Invalid decay mode. Must be one of: {valid_decays}')

    infile = args.infile
    decay = args.decay

    # Create directory if it doesn't exist
    outdir = f'out/plot/{decay}/'
    if args.signal: outdir += 'signal/'
    else: outdir += 'minbias/'
    os.makedirs(outdir, exist_ok=True)

    fileheader = f'{outdir}plot_{name}'

    return infile, fileheader, decay, args.stats, args.legend