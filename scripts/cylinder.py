#! /usr/bin/env python3

headerhelp = \
''' 
    Calculates a cylinder to cover the protein.
'''
# There should be a way to do this with relative imports
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('inpath',
                    help='The input PDB file.')
parser.add_argument('--Dmin', type=float, default=3.2,
                    help='Size of cushion, defaults to 3.2A.')
parser.add_argument('--ranges',
                    help='Residue range selection, chain, start, end. \
                          Comma-separated, no spaces. Separate chains with \
                          forward slash. Example: \
                          A,50-55,72-80/B,50-55 will select residues \
                          50 to 55 and 72 to 80 in chain A, and also \
                          residues 50 to 55 in chain B.')
args = parser.parse_args()

from scipy import sqrt, array, pi

from pdbtool import ReadPDBfile
model = ReadPDBfile(args.inpath)

if args.ranges:
    from helper import parse_ranges
    model = model.extract_range(parse_ranges(args.ranges, model.GetChains()))

tin = model.GetInertiaTensor()
model.shift(-tin.comass)
model.transform(tin.vr)
x,y,z = model.GetCoordinateArray().T
Rc = sqrt(array([max(y**2+z**2), max(x**2+z**2), max(x**2+y**2)])) + args.Dmin
Hc = array([x.ptp(), y.ptp(), z.ptp()]) + 2*args.Dmin
Vc = pi*Rc**2*Hc
Ri = array([x if x>0 else 0.0 for x in sqrt(array([min(y**2+z**2), min(x**2+z**2), min(x**2+y**2)])) - args.Dmin])
print("Cylinders:\n%9s %9s %9s %9s" % ("Diameter","Height","Volume","Inner diameter"))
print("%9.1f %9.1f %9.2f %9.2f\n"*3 % tuple(array([2*Rc,Hc,Vc/1e6,2*Ri]).T[Vc.argsort()].flatten()))
