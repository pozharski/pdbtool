#! /usr/bin/env python3

headerhelp = \
''' 
    PDB file analysis plots.  Actions to perform on the input file are 
    defined by -p option. Allowed actions are:
    
    bvalues             Plots average per-residue Bfactors, 
                        including all atoms, backbone and side chain 
                        columns.
'''
# There should be a way to do this with relative imports
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('inpath',
                    help='The input PDB file.')
parser.add_argument('-p', '--plot', action='append',
                    choices = ['bvalues'],
                    default = [],
                    metavar='', help='Plotting options.')
args = parser.parse_args()

from pdbtool import ReadPDBfile
model = ReadPDBfile(args.inpath)

from pdbwindows import BWindow
from matplotlib.pyplot import figure, show

for whatoprint in args.plot:
    if whatoprint == 'bvalues':
        bv_fig = figure(FigureClass=BWindow)
        bv_fig.set_model(model)
        bv_fig.plot()
        show()
