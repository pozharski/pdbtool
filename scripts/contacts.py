#! /usr/bin/env python3

headerhelp = \
''' 
    Calculates all the close contacts.
'''
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('inpath',
                    help='The input PDB file.')
parser.add_argument('--Dmax', type=float, default=4.0,
                    help='Maximum contact distance, defaults to 4.0A.')
args = parser.parse_args()

from pdbtool import ReadPDBfile
model = ReadPDBfile(args.inpath)

print('\n'.join([' : '.join([model.GetAtomTitle(i),model.GetAtomTitle(j),"%.2f" % d]) for i,j,d in model.AllContacts(rmax=args.Dmax)]))
