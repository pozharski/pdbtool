#! /usr/bin/env python3

headerhelp = \
''' 
    HBFILTER selects a subset of contacts from a database
'''
import os, sys, shutil
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
parser.add_argument('--in-place', action='store_true',
                    help='Modify the database file in place, use with caution')
parser.add_argument('--mmsize',
                    help='Oligomerization number.  Could be a comma separated list.')
parser.add_argument('--res2',
                    help='Filter by second residue type')
parser.add_argument('--same-residue', action='store_true',
                    help='Filter for identical residues forming a bond')
parser.add_argument('--not-same-residue', action='store_true',
                    help='Filter for identical residues forming a bond')
parser.add_argument('--do-checks', action='store_true',
                    help='Check the PISA database integrity.')
parser.add_argument('--pvalue', type=float,
                    help='pvalue cutoff.  This requires that the type of hydrogen is specified.')
parser.add_argument('--hbtype', '-b',
                    help="Hydrogen bond type.")
parser.add_argument('--symmetric',
                    action='store_true',
                    help='Bonds are chemically symmetric so both directions should be included.')
parser.add_argument('--no-metals', action='store_true',
                    help='Filter out false bonds that result from metal binding sites')

args = parser.parse_args()

from pdbminer import pisa
from aconts import hbond_pdbase

pisabase = pisa.pisa_dbsres_pdbase(args.sqlpisa)
if args.do_checks:
    pisabase.download_check()
if args.in_place:
    hpdbout = hbond_pdbase(args.sqlpath)
else:
    shutil.copyfile(args.sqlpath, args.sqlout)
    hpdbout = hbond_pdbase(args.sqlout)

#hbs = hpdbase.get_hbonds()

if args.no_metals:
    for code in hpdbase.fetch_processed_codes():
        print("Processing %s..." % code)
        molpath = pisabase.get_path(code)
        mol = ReadPDBfile(molpath[0])
        metals = mol.atom_getter('metals')
        if len(metals):
            pass
        else:
            print("No metals found, skip")

if args.mmsize is not None:
    mmsizes = [int(x) for x in args.mmsize.split(',')]
    print("Only include oligomers of the following sizes: %s" % args.mmsize)
    pdbs = pisabase.get_pdbs()
    print("Start with %d PDB codes" % len(pdbs))
    codes2remove = [x[0] for x in pdbs if int(x[1].mmsize) not in mmsizes]
    num2remove = len(codes2remove)
    print("Remove %d PDB codes" % num2remove)
    print("Filtering %d bond entries" % hpdbout.get_hbond_number())
    for (i, code) in enumerate(codes2remove):
        hpdbout.delete_by_code(code)
        hpdbout.commit()
        print("%8d codes still to remove" % (num2remove-i-1), end='\r')
    print("Reduced to %d bond entries" % hpdbout.get_hbond_number())

if args.res2 is not None:
    print("Only include %s as the second residue" % (args.res2))
    print("Start with %d hydrogen bonds" % hpdbout.get_hbond_number())
    hpdbout.filterby(res2=args.res2)
    print("%d hydrogen bonds left after filtering" % hpdbout.get_hbond_number())

if args.same_residue:
    hpdbout.filter_same_residue()
if args.not_same_residue:
    hpdbout.filter_same_residue(False)

if args.pvalue:
    if args.hbtype:
        hpdbout.filter_pvalues(args.hbtype, args.pvalue, args.symmetric)
    else:
        print("Hydrogen bond type not specified, skipping p-value filtering.")
        

hpdbout.close()
