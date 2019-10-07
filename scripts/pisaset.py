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
parser.add_argument('inpath',
                    help='The input files listing PISA entries.')
parser.add_argument('--sqlpath', default='pisa.sqlite',
                    help='Database file')
args = parser.parse_args()

from pdbminer import pisa

dbase = pisa.pisa_dbsres_pdbase(args.sqlpath)
dbase.process_codes(args.inpath)
