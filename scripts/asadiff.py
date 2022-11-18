#! /usr/bin/env python3

headerhelp = \
''' 
    Buried surface area calculator.
'''
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import freesasa

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('inpath1',
                    help='Initial conformation PDB file.')
parser.add_argument('inpath2',
                    help='Final conformation PDB file.')
parser.add_argument('--chids',
                    help='Chain ID match.')
parser.add_argument('--bpdb', 
                    action='store_true',
                    help='Output PDB files with ASA written as B-factors.')

args = parser.parse_args()


asa1 = freesasa.calc(freesasa.Structure(args.inpath1))
asa2 = freesasa.calc(freesasa.Structure(args.inpath2))

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

