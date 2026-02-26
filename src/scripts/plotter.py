import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '--apply-cuts',
    action='store_true',
    help='Apply fiducial cuts before histogramming'
)
parser.add_argument(
    '-s', '--signal',
    action='store_true',
    help='Process signal data'
)
args = parser.parse_args()

# User inputted header for input files
print('Enter the header path for input files (e.g. ntuple/mc_minbias_20260124/):')
header = input('Path: ').strip()
if not header.endswith('/'): header += '/'

PY = ["lb-conda", "default", "python3"]

commands = []
if args.apply_cuts: 
    infile = header + 'combined_files.root'
    commands.append(['src/fid_reqs.py', '-i', infile])
if args.signal:
    infile = header + ('fiducial_requirements.root')
    commands += [
        ['src/hist_rec.py', '-i', infile, '-s'],
        ['src/plot_rec.py', 'stats', 'legend', 'signal'],
        ['src/hist_gen.py', '-i', infile, '-s'],
        ['src/plot_gen.py', 'stats', 'legend', 'signal'],
        ['src/hist_mass.py', '-i', infile, '-s'],
        ['src/plot_mass.py', 'stats', 'legend', 'signal'],
    ]
else:
    infile = header + ('fiducial_requirements.root')
    commands += [
        ['src/hist_rec.py', '-i', infile],
        ['src/plot_rec.py', 'stats', 'legend'],
        ['src/hist_gen.py', '-i', infile],
        ['src/plot_gen.py', 'stats', 'legend'],
        ['src/hist_mass.py', '-i', infile],
        ['src/plot_mass.py', 'stats', 'legend'],
    ]

for args in commands:
    subprocess.run(PY + args, check=True)
