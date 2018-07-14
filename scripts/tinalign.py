#! /usr/bin/env python3

headerhelp = \
'''
TINALIGN aligns two structures by aligning their tensors of inertia.
'''

# There should be a way to do this with relative imports
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('fixmodel',
                    help='Fixed model PDB file.')
parser.add_argument('movmodel',
                    help='Moving model PDB file.')
parser.add_argument('outmodel',
                    help='Output model PDB file.')
parser.add_argument('--chids',
                    help='Selection of chids')
parser.add_argument('--spin', choices = ['x','y','z'],
                    help='Spin around selected axis.')
parser.add_argument('--swap', action='store_true',
                    help='Swap similar axis')

args = parser.parse_args()

from pdbtool import ReadPDBfile, pdbmolecule
from scipy import array

# Import models
fixmodel = ReadPDBfile(args.fixmodel)
movmodel = ReadPDBfile(args.movmodel)
# Get inertia tensors
if args.chids is not None:
    fixtin = fixmodel.GetInertiaTensor('chids', chids=args.chids)
    movtin = movmodel.GetInertiaTensor('chids', chids=args.chids)
else:
    fixtin = fixmodel.GetInertiaTensor()
    movtin = movmodel.GetInertiaTensor()
# Center the moving models
movmodel.shift(-movtin.comass)
if args.swap:
    if movtin.abc[0]+movtin.abc[2] < 2*movtin.abc[1]:
        vr = movtin.vr[:,[1,0,2]]
    else:
        vr = movtin.vr[:,[0,2,1]]
else:
    vr = movtin.vr
movmodel.transform(fixtin.vr.T.dot(vr))
if args.spin == 'x':
    movmodel.transform(array([[1,0,0],[0,-1,0],[0,0,-1]]))
elif args.spin == 'y':
    movmodel.transform(array([[-1,0,0],[0,1,0],[0,0,-1]]))
elif args.spin == 'z':
    movmodel.transform(array([[-1,0,0],[0,-1,0],[0,0,1]]))
movmodel.shift(fixtin.comass)

movmodel.writePDB(args.outmodel)

