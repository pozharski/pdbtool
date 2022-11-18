#! /usr/bin/env python3

# There should be a way to do this with relative imports
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import epyactions

def main():
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    headerhelp = \
'''
PDB file comparisons.  Actions to perform on two input files are 
defined by -a option. Allowed actions are:

    bdifference         Calculate and assign B-factor differences


------------------------------------------------------------------------

'''
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                            description=headerhelp)
    parser.add_argument('inpath1',
                        help='First PDB file.')
    parser.add_argument('inpath2',
                        help='Second PDB file.')
    parser.add_argument('outpath1', nargs='?',
                        default='pdbtrickout1.pdb',
                        help='The first output PDB file.')
    parser.add_argument('outpath2', nargs='?',
                        default='pdbtrickout2.pdb',
                        help='The second output PDB file.')
    parser.add_argument('-a', '--action', action='append',
                        choices = [	'bdifference',
                                    ],
                        default = [],
                        metavar = '', help='Action to perform')
    args = parser.parse_args()

    from pdbtool import ReadPDBfile as read_pdb_file

    model1 = read_pdb_file(args.inpath1)
    if model1 is None:
        sys.exit('Failed to read model from '+args.inpath1)
    model2 = read_pdb_file(args.inpath2)
    if model1 is None:
        sys.exit('Failed to read model from '+args.inpath2)

    for action in args.action:
        epyactions.__getattribute__(action.replace('-','_').lower())(args,  [model1,model2])

if __name__ == "__main__":
    main()
