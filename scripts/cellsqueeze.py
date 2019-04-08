#! /usr/bin/env python3

headerhelp = \
''' 
    Reduces the unit cell to match the protein shape.  Currently
    only works in P1, which is consistent with the script's primary
    purpose of removing extra space from cryoEM models.
'''
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('inpath',
                    help='Input PDB file.')
parser.add_argument('outpath',
                    help='Output PDB file.')
parser.add_argument('--mrc', 
                    help='Corresponding map file')
parser.add_argument('--output-mrc',
                    help='Output map file')
parser.add_argument('--cushion', type=float, default=2.0,
                    help='Cushion size around the molecule.')
args = parser.parse_args()

from pdbtool import ReadPDBfile, pdbmolecule, cell_and_center
from scipy import floor, ceil, array

model = ReadPDBfile(args.inpath)

if args.mrc is not None:
    import mrc
    emmap = mrc.mrc(args.mrc)
    hxyz = array([emmap.ucell.a/emmap.mx, emmap.ucell.b/emmap.my, emmap.ucell.c/emmap.mz])
else:
    emmap = None

r = model.GetCoordinateArray()

lfb = r.min(0)-args.cushion
if emmap is not None:
    map_lfb = floor(lfb/hxyz).astype(int)
    lfb = map_lfb*hxyz

rbt = r.max(0)+args.cushion
if emmap is not None:
    map_rbt = ceil(rbt/hxyz).astype(int)
    rbt = map_rbt*hxyz

newabc = rbt-lfb
model.shift(-lfb)
celline = 'CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1                     \n' % tuple(newabc)
model.writePDB(args.outpath, header=celline)

if emmap is not None:
    emmap.nx, emmap.ny, emmap.nz = map_rbt-map_lfb
    emmap.mx, emmap.my, emmap.mz = map_rbt-map_lfb
    emmap.ucell.a, emmap.ucell.b, emmap.ucell.c = array([emmap.mx,emmap.my,emmap.mz])*hxyz
    emmap.data = emmap.data[map_lfb[2]:map_rbt[2],map_lfb[1]:map_rbt[1],map_lfb[0]:map_rbt[0]]
    if args.output_mrc is None:
        args.output_mrc = args.mrc+'.out'
    emmap.write(args.output_mrc)
    
    
