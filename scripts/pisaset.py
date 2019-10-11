#! /usr/bin/env python3

headerhelp = \
''' 
    Updates PISA datasets.
'''
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('--inpath', '-i',
                    help='The input files listing PISA entries.')
parser.add_argument('--sqlpath', default='pisa.sqlite',
                    help='Database file')
parser.add_argument('--no-checks', action='store_true',
                    help='Do not check the database integrity.')
args = parser.parse_args()

from pdbminer import pisa

dbase = pisa.pisa_dbsres_pdbase(args.sqlpath)
if not args.no_checks:
    dbase.download_check()
if args.inpath is not None:
    dbase.process_codes(args.inpath)
