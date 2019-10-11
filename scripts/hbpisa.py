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
                    help='HB database file')
parser.add_argument('--sqlpisa', default='pisa.sqlite',
                    help='PISA database file')
parser.add_argument('--do-checks', action='store_true',
                    help='Check the PISA database integrity.')
args = parser.parse_args()

from pdbminer import pisa
from pdbtool import ReadPDBfile

pisabase = pisa.pisa_dbsres_pdbase(args.sqlpisa)
if args.do_checks:
    pisabase.download_check()
code = pisabase.fetch_codes()[0]
print(pisabase.get_path(code))
mol = ReadPDBfile(pisabase.get_path(code))

