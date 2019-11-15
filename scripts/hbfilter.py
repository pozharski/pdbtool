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
parser.add_argument('--mmsize',
                    help='Oligomerization number.  Could be a comma separated list.')
parser.add_argument('--do-checks', action='store_true',
                    help='Check the PISA database integrity.')
args = parser.parse_args()

from pdbminer import pisa
from aconts import hbond_pdbase

pisabase = pisa.pisa_dbsres_pdbase(args.sqlpisa)
hpdbase = hbond_pdbase(args.sqlpath)

for code in 
