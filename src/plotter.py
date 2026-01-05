import subprocess

PY = ["lb-conda", "default", "python3"]

commands = [
    # Signal
    # ['src/rec_hist.py', 'sig'],
    # ['src/rec_plot.py', 'sig', 'stats', 'legend'],
    # ['src/gen_hist.py', 'sig'],
    # ['src/gen_plot.py', 'sig', 'stats', 'legend'],
    # ['src/mass_hist.py', 'sig'],
    # ['src/mass_plot.py', 'sig', 'stats', 'legend'],
    # Minbias
    ['src/hist_rec.py'],
    ['src/plot_rec.py', 'stats', 'legend'],
    ['src/hist_gen.py'],
    ['src/plot_gen.py', 'stats', 'legend'],
    ['src/hist_mass.py'],
    ['src/plot_mass.py', 'stats', 'legend'],
]

for args in commands:
    subprocess.run(PY + args, check=True)
