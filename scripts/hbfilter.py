#! /usr/bin/env python3

headerhelp = \
''' 
    HBFILTER selects a subset of contacts from a database
'''
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('--sqlpath', default='hbpisa.sqlite',
                    help='HB database file.  Defaults to hbpisa.sqlite.')
parser.add_argument('--sqlpisa', default='pisa.sqlite',
                    help='PISA database file.  Defaults to pisa.sqlite.')
parser.add_argument('--sqlout', default='hbfiltered.sqlite',
                    help='Output HB database file.  Defaults to hbfiltered.sqlite.')
parser.add_argument('--mmsize',
                    help='Oligomerization number.  Could be a comma separated list.')
parser.add_argument('--same-residue', action='store_true',
                    help='Filter for identical residues forming a bond')
parser.add_argument('--do-checks', action='store_true',
                    help='Check the PISA database integrity.')
args = parser.parse_args()

from pdbminer import pisa
from aconts import hbond_pdbase

pisabase = pisa.pisa_dbsres_pdbase(args.sqlpisa)
hpdbase = hbond_pdbase(args.sqlpath)
hpdbout = hbond_pdbase(args.sqlout)

pdbs = pisabase.get_pdbs()
hbs = hpdbase.get_hbonds()

if args.mmsize is not None:
    mmsizes = [int(x) for x in args.mmsize.split(',')]
    print("Only include oligomers of the following sizes: %s" % args.mmsize)
    print("Start with %d PDB codes" % len(pdbs))
    pdbs = [x for x in pdbs if int(x[1].mmsize) in mmsizes]
    print("Reduced to %d PDB codes" % len(pdbs))
    codes = [x[0] for x in pdbs]
    print("Filtering %d bond entries" % len(hbs))
    hbs=[x for x in hbs if x[0] in codes]
    print("Reduced to %d bond entries" % len(hbs))

if args.same_residue:
    print("Filtering %d bond entries" % len(hbs))
    hbs = [x for x in hbs if x[1].resid1[1:] == x[1].resid2[1:]]
    print("Reduced to %d bond entries" % len(hbs))

for code in set([x[0] for x in hbs]):
    hpdbout.insert_new_code(code, 1)

for hb in hbs:
    hpdbout.insert_hb(hb[0],hb[1])

hpdbout.close()
