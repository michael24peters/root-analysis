"""
plot_utils.py

Generic plotting helpers: font setup, JSON I/O, param-dict unpacking, and the
argparse for a JSON input file.

For "house style" presentation logic (stats-box text, text placement,
save/report boilerplate), see plot_style.py.
"""

import argparse
import json

import matplotlib.pyplot as plt

FONTS = {
    'serif': {'family': 'STIXGeneral', 'mathtext.fontset': 'stix'},
    'sans':  {'family': 'Inter', 'mathtext.fontset': 'dejavusans'},
}

def apply_font(font):
    """Set font.family/mathtext.fontset for the given style ('serif' or 'sans')."""
    style = dict(FONTS.get(font, FONTS['serif']))
    if style['family'] == 'Inter':
        try:
            import matplotlib.font_manager as fm
            fm.findfont('Inter', fallback_to_default=False)
        except Exception:
            print('[INFO] Inter font not found, falling back to default sans-serif.')
            style['family'] = 'sans-serif'
    plt.rcParams['font.family'] = style['family']
    plt.rcParams['mathtext.fontset'] = style['mathtext.fontset']

def parse_args():
    parser = argparse.ArgumentParser(description='Plot fit results from JSON file')
    parser.add_argument('input', help='Input JSON file')
    return parser.parse_args()

def load_data(infile):
    with open(infile) as f: return json.load(f)

def val_err(param):
    return param.get('value', float('nan')), param.get('error', float('nan'))
