#! /usr/bin/env python3

headerhelp = \
''' 
    Calculate buried surface area upon complex formation.
    Components of the complex are defined via first/second parameters
    with syntax of chain,start,end (comma-separated, no spaces). 
    Separate chains with forward slash. 
    Example: A,50-55,72-80/B,50-55 
    will select residues 50 to 55 and 72 to 80 in chain A, and also 
                residues 50 to 55              in chain B.
    If second selection is omitted, second component will comprise the
    rest of the structure.
'''
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import freesasa
import numpy as np

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('inpath',
                    help='Input PDB file for complex.')
parser.add_argument('inpath2', nargs='?',
                    help='Input PDB file for first component (optional).')
parser.add_argument('inpath3', nargs='?',
                    help='Input PDB file for second component (optional).')
parser.add_argument('--first',
                    help='First component definition, required unless components are provided explicitly as PDB files.')
parser.add_argument('--second', 
                    help='Second component definition.  If omitted, assumes that it is the rest of the complex.')
parser.add_argument('--csv', 
                    help='Output CSV file')
parser.add_argument('--bpdb', 
                    help='Output PDB file name with BSA written as B-factors.')
parser.add_argument('--verbosity', '-v', action='count', default=0,
                    help='Verbosity level.')


args = parser.parse_args()

from pdbtool import ReadPDBfile
from tempfile import mktemp

mol = ReadPDBfile(args.inpath)
if args.inpath2:
    mol1 = mol
    mol2 = ReadPDBfile(args.inpath2)
    fout1, fout2 = args.inpath, args.inpath2
else:
    range1 = mol.parse_range(args.first)
    mol1 = mol.extract_range(range1)
    fout1 = mktemp(os.extsep+'pdb')
    mol1.writePDB(fout1)
    if args.second is None:
        range2 = mol.parse_range(args.first)
        mol2 = mol.extract_range(range1, inverse=True)
    else:
        range2 = mol.parse_range(args.second)
        mol2 = mol.extract_range(range2)
    fout2 = mktemp(os.extsep+'pdb')
    mol2.writePDB(fout2)
if args.inpath3:
    fout3 = args.inpath3
    mol3 = ReadPDBfile(args.inpath3)
else:
    mol3 = mol1.copy()
    mol3.AppendMolecule(mol2)
    fout3 = mktemp(os.extsep+'pdb')
    mol3.writePDB(fout3)

#print(fout1, fout2,fout3)

asa1 = freesasa.calc(freesasa.Structure(fout1))
asa2 = freesasa.calc(freesasa.Structure(fout2))
asa3 = freesasa.calc(freesasa.Structure(fout3))

#print(dir(list(list(asa1.residueAreas().items())[0][1].items())[0][1]))

d1 = dict([(k,dict(v.items())) for k,v in asa1.residueAreas().items()])
d2 = dict([(k,dict(v.items())) for k,v in asa2.residueAreas().items()])

fnan = float('nan')

bsa = {}

for chid,chv in asa3.residueAreas().items():
    bsa[chid] = {}
    for resi, resv in chv.items():
        resv1 = d1[chid].get(resi,None) if chid in d1 else None
        resv2 = d2[chid].get(resi,None) if chid in d2 else None
        resvC = resv1 if resv1 else resv2
        bsa[chid][resi] = [resvC.__getattribute__(k) if resvC else fnan for k in ['total','apolar','polar','mainChain','sideChain']]+[1 if resv1 else 2]
        if args.verbosity>1:
            print("%1s %-5d %5.1f %5.1f %5.1f %5.1f %5.1f %2d" % tuple([chid,int(resi)]+bsa[chid][resi]))

comp1, comp2  = {}, {}
for chid,chv in bsa.items():
    comp1[chid] = [sum(x) for x in zip(*[resv[:-1] for resi, resv in chv.items() if resv[-1]==1])]
    comp2[chid] = [sum(x) for x in zip(*[resv[:-1] for resi, resv in chv.items() if resv[-1]==2])]

print("Component    Total     Apolar     Polar    Backbone  Siode chain")

print("    1      %9.1f %9.1f %9.1f %9.1f %9.1f" % tuple([sum(x) for x in zip(*[v for k,v in comp1.items() if v])]))
print("    2      %9.1f %9.1f %9.1f %9.1f %9.1f" % tuple([sum(x) for x in zip(*[v for k,v in comp2.items() if v])]))


if args.csv:
    import csv
    with open(args.csv,"w") as fcsv:
        csv_writer = csv.writer(fcsv)
        csv_writer.writerow(['Chain ID','Residue number','total','apolar','polar','mainChain','sideChain','Component'])
        for chid,chv in bsa.items():
            for resi, resv in chv.items():
                csv_writer.writerow([chid,resi]+resv)

sys.exit()

print("--- Total area ---")
print("%.0f  ---> %.0f" % (asa1.totalArea(), asa2.totalArea()))
print("--- Initial conformation ---")
print("ChainID      Apolar     Polar       Main      Side       Total")
print('\n'.join(['%s       %10.0f %10.0f %10.0f %10.0f %10.0f ' % t for t in [(k,sum([r.apolar for r in v.values()]),
                                                             sum([r.polar for r in v.values()]),
                                                             sum([r.mainChain for r in v.values()]),
                                                             sum([r.sideChain for r in v.values()]),
                                                             sum([r.total for r in v.values()])) for k, v in asa1.residueAreas().items()]]))
print("--- Final conformation ---")
print("ChainID      Apolar     Polar       Main      Side       Total")
print('\n'.join(['%s       %10.0f %10.0f %10.0f %10.0f %10.0f ' % t for t in [(k,sum([r.apolar for r in v.values()]),
                                                             sum([r.polar for r in v.values()]),
                                                             sum([r.mainChain for r in v.values()]),
                                                             sum([r.sideChain for r in v.values()]),
                                                             sum([r.total for r in v.values()])) for k, v in asa2.residueAreas().items()]]))
print("--- Buried area ---")
if args.chids:
    chid1,chid2 = args.chids.split(',')
    print("Apolar   : %10.0f A^2" % (sum([r.apolar for r in asa1.residueAreas()[chid1].values()]) - sum([r.apolar for r in asa2.residueAreas()[chid2].values()])))
    print("Polar    : %10.0f A^2" % (sum([r.polar for r in asa1.residueAreas()[chid1].values()]) - sum([r.polar for r in asa2.residueAreas()[chid2].values()])))
    print("Mainchain: %10.0f A^2" % (sum([r.mainChain for r in asa1.residueAreas()[chid1].values()]) - sum([r.mainChain for r in asa2.residueAreas()[chid2].values()])))
    print("Sidechain: %10.0f A^2" % (sum([r.sideChain for r in asa1.residueAreas()[chid1].values()]) - sum([r.sideChain for r in asa2.residueAreas()[chid2].values()])))
    print("Total    : %10.0f A^2" % (sum([r.total for r in asa1.residueAreas()[chid1].values()]) - sum([r.total for r in asa2.residueAreas()[chid2].values()])))
else:
    print("Apolar   : %10.0f A^2" % (sum([sum([r.apolar for r in v.values()]) for k, v in asa1.residueAreas().items()])-sum([sum([r.apolar for r in v.values()]) for k, v in asa2.residueAreas().items()])))
    print("Polar    : %10.0f A^2" % (sum([sum([r.polar for r in v.values()]) for k, v in asa1.residueAreas().items()])-sum([sum([r.polar for r in v.values()]) for k, v in asa2.residueAreas().items()])))
    print("Mainchain: %10.0f A^2" % (sum([sum([r.mainChain for r in v.values()]) for k, v in asa1.residueAreas().items()])-sum([sum([r.mainChain for r in v.values()]) for k, v in asa2.residueAreas().items()])))
    print("Sidechain: %10.0f A^2" % (sum([sum([r.sideChain for r in v.values()]) for k, v in asa1.residueAreas().items()])-sum([sum([r.sideChain for r in v.values()]) for k, v in asa2.residueAreas().items()])))
    print("Total    : %10.0f A^2" % (sum([sum([r.total for r in v.values()]) for k, v in asa1.residueAreas().items()])-sum([sum([r.total for r in v.values()]) for k, v in asa2.residueAreas().items()])))

if args.bpdb:
    from pdbtool import ReadPDBfile
    model1 = ReadPDBfile(args.inpath1)
    model1.SetBfactorValues([asa1.atomArea(i)-asa2.atomArea(i) for i in range(model1.GetAtomNumber())])
    model1.writePDB(args.inpath1+os.extsep+'asadiff'+os.extsep+'pdb')
    model2 = ReadPDBfile(args.inpath2)
    model2.SetBfactorValues([asa1.atomArea(i)-asa2.atomArea(i) for i in range(model1.GetAtomNumber())])
    model2.writePDB(args.inpath2+os.extsep+'asadiff'+os.extsep+'pdb')

