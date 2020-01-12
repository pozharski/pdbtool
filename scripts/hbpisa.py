#! /usr/bin/env python3

headerhelp = \
''' 
    Builds HB database for dataset derived from PISA
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
parser.add_argument('--hbtype', '-b',
                    help="Hydrogen bond type (mandatory).")
parser.add_argument('--do-checks', action='store_true',
                    help='Check the PISA database integrity.')
args = parser.parse_args()

from pdbminer import pisa
from pdbtool import ReadPDBfile
from aconts import hbond_pdbase, HBreader, DAKERNEL_PARAMS

if args.hbtype is None:
    sys.exit('Hydrogen bond type unknown.  Are you sure you understand how this works?')
try:
    exec('from aconts import ' + args.hbtype + ' as HydroBonds')
except ImportError:
    sys.exit('Hydrogen bond type '+args.hbtype+' is not defined.  Please choose from the following list:\n'+'\n'.join(sorted(DAKERNEL_PARAMS.keys())))
    
pisabase = pisa.pisa_dbsres_pdbase(args.sqlpisa)
if args.do_checks:
    pisabase.download_check()
hpdbase = hbond_pdbase(args.sqlpath)
for code in hpdbase.filter_codes(pisabase.fetch_codes()):
    print("Processing %s..." % code)
    molpath = pisabase.get_path(code)
    mol = ReadPDBfile(molpath[0])
    if mol.is_multi_model() or mol.GetAtomNumber()==0:
        if mol.GetAtomNumber()==0:
            print(code + ' is an empty structure.  Failed download?',file=sys.stderr)
        else:
            print(code + ' appears to be multi-model, skip', file=sys.stderr)
    else:
        hb = HydroBonds(mol)
        hpdbase.insert_new_code(code)
        for item in hb.get_readers():
            hpdbase.insert_hb(code, item)
        hpdbase.code_lock(code)
        hpdbase.commit()
