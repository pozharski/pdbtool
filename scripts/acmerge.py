#! /usr/bin/env python3

headerhelp = \
''' 
    Merges two models as alternative conformations
'''
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('inpath1',
                    help='The input PDB file of the first model.')
parser.add_argument('inpath2',
                    help='The input PDB file of the second model.')
parser.add_argument('outpath',
                    help='The output PDB file.')
parser.add_argument('--acletters', default='AB',
                    help='Letters used to designate multiple conformers.')
args = parser.parse_args()

altloc1, altloc2 = args.acletters[:2]

from pdbtool import ReadPDBfile, pdbmolecule

model1 = ReadPDBfile(args.inpath1)
model2 = ReadPDBfile(args.inpath2)

for atom in model1.atoms:
    atom.SetAltLoc(altloc1)
for atom in model2.atoms:
    atom.SetAltLoc(altloc2)

model = pdbmolecule()

while model1.GetAtomNumber():
    atom = model1.PopAtom(0)
    model.AppendAtom(atom)
    for i in model2.findalt(atom)[::-1]:
        model.AppendAtom(model2.PopAtom(i))
    print(model.GetAtomNumber(),model1.GetAtomNumber(),model2.GetAtomNumber(),"               ",end='\r')
model.AppendAtom(model2.GetAtoms())
print(model.GetAtomNumber(),model1.GetAtomNumber(),model2.GetAtomNumber())
model.SerialReset()
print("Writing output model...                 ", end="\r")
model.writePDB(args.outpath)
print("                        done.")


