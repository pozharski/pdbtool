'''
DATREFINE reads a hb database that contains distance/angle/torsion 
columns and refines the DAT parameters.
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('--sqlpath', default='hbpisa.sqlite',
                    help='HB database file.  Defaults to hbpisa.sqlite.')

parser.add_argument('--torshift', type=float,
                    help='Shift the torsion angle domain')
parser.add_argument('--kmeans',
                    type=int, default=2,
                    help='Number of clusters to use to perform K-means clustering analysis.')

args = parser.parse_args()

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from aconts import hbond_pdbase
