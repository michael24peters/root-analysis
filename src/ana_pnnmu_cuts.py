import subprocess
import argparse
from pathlib import Path
import re

parser = argparse.ArgumentParser()
parser.add_argument(
    '-i', '--infile',
    help='Input ROOT file'
)
parser.add_argument(
    '-o', '--outfile',
    help='Output text file'
)
parser.add_argument(
    '-d', '--dir',
    help='Directory with text outputs to parse'
)
args = parser.parse_args()

infile = args.infile
outfile = args.outfile
indir = args.dir

if indir:
    indir_path = Path(indir)
    if not indir_path.is_dir():
        raise ValueError(f'Input directory not found: {indir}')

    txt_files = sorted(indir_path.glob('*.txt'))
    if not txt_files:
        raise ValueError(f'No .txt files found in directory: {indir}')

    candidates_re = re.compile(r'^Total candidates processed:\s*(\d+)\s*$')
    signal_re = re.compile(r'^Total signal candidates:\s*(\d+)\s*$')

    results = {}
    for txt_file in txt_files:
        candidates = None
        signal = None
        with txt_file.open('r', encoding='utf-8') as handle:
            for line in handle:
                line = line.strip()
                if candidates is None:
                    match = candidates_re.match(line)
                    if match:
                        candidates = int(match.group(1))
                        continue
                if signal is None:
                    match = signal_re.match(line)
                    if match:
                        signal = int(match.group(1))
                        continue

        results[txt_file.name] = {
            'candidates': candidates,
            'signal': signal,
        }

    print('Parsed candidate/signal counts:')
    for name, counts in results.items():
        print(f'{name}: candidates={counts["candidates"]}, signal={counts["signal"]}')
    raise SystemExit(0)

PY = ["lb-conda", "default", "python3"]

commands = []
pnnmu_cuts = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975]
outfile_ext = outfile.split('.')[-1]
outfile_name = '.'.join(outfile.split('.')[:-1])

for cut in pnnmu_cuts:
    # Put cut value between file extension and name
    # e.g. out/bkg_ana_0.400.txt
    outfile_cut = f'{outfile_name}_{cut:.3f}.{outfile_ext}'
    commands.append(['src/bkg_ana.py', '-i', infile, '-o', outfile_cut, 
                     '--pnnmu_cut', str(cut)])

for args in commands:
    subprocess.run(PY + args, check=True)
