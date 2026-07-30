"""
pnnmu_cuts.py

TODO: Outdated script to run bkg.py with a range of PROBNNmu cuts and save the
output to text files. Will be updated with a more flexible and robust generic
cutting script.
"""

import subprocess
import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    '-i', '--infile',
    help='Input ROOT file'
)
parser.add_argument(
    '-o', '--outfile',
    help='Output text file'
)
args = parser.parse_args()

infile = args.infile
outfile = args.outfile

PY = ['lb-conda', 'default', 'python3']

commands = []
pnnmu_cuts = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975]
outfile_ext = outfile.split('.')[-1]
outfile_name = '.'.join(outfile.split('.')[:-1])

for cut in pnnmu_cuts:
    # Put cut value between file extension and name
    # e.g. out/bkg_ana_0.400.txt
    outfile_cut = f'{outfile_name}_{cut:.3f}.{outfile_ext}'
    commands.append(['src/ana/bkg.py', infile, '-o', outfile_cut,
                     '--pnnmu_cut', str(cut)])

for args in commands:
    subprocess.run(PY + args, check=True)
