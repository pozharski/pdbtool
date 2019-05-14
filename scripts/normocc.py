#! /usr/bin/env python3

headerhelp = \
''' 
    Normalize occupancies for alternate conformers that exceed total
    occupancy of 1.    
'''
# There should be a way to do this with relative imports
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('inpath',
                    help='The input PDB file.')
parser.add_argument('outpath', nargs='?',
                    help='The output PDB file.  If not specified, input model is simply checked for occupancy issues.')
parser.add_argument('-f', '--force-overwrite', 
                    action='store_true',
                    help='Overwrite existing files.')
args = parser.parse_args()

from pdbtool import ReadPDBfile
model = ReadPDBfile(args.inpath)

over_occupied = [x for x in model.GetAtoms() if float(x.GetOccupancy())>1]

if len(over_occupied):
    print("Found %d atoms with occupancy>1, corrected" % len(over_occupied))
    for atom in over_occupied:
        atom.SetOccupancy(1.0)

acatoms = model.atom_getter('altconf')

print("Found %d atoms with non-trivial alternate conformer IDs" % len(acatoms))

acgroups = []
for atom in acatoms:
    if atom not in sum(acgroups,[]):
        acgroups.append(model.atom_getter('altgroup', atom=atom))

print("Found %d groups of alternate conformer atoms" % len(acgroups))

over_occupied = [g for g in acgroups if sum([a.GetOccupancy() for a in g])>1]

for g in over_occupied:
    total_occupancy = sum([a.GetOccupancy() for a in g])
    for a in g:
        a.SetOccupancy(round(a.GetOccupancy()/total_occupancy,2))

print("Found %d atom groups with occupancy>1" % len(over_occupied))

resids = set([x.resid() for x in sum(over_occupied,[])])
print("Problematic residues:")
print('\n'.join(resids))

if args.outpath:
    if not os.access(args.outpath, os.F_OK) or args.force_overwrite:
        model.writePDB(args.outpath)
        print("Corrected model saved as %s" % args.outpath)
    else:
        print("File %s already exists, skipping output" % args.outpath)





