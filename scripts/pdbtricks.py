#! /usr/bin/env python3

# There should be a way to do this with relative imports
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import epyactions

def main():
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    headerhelp = \
'''PDB file manipulations.  Actions to perform on the input file are 
defined by -a option. Allowed actions are:

    extract-chains      Extract only chains specified by --chids option
    extract-ranges      Extracts atoms that belong to the list of 
                        residue ranges
    rename-chains       Rename chains using the pattern defined by
                        --chids option.  For example, "AB,CD" will
                        rename A to B and C to D.
    rjust-resid         Makes sure residue names are right-justified
    tinertia-ranges     Outputs tensor of inertia information for 
                        specified ranges
    tinertia-slider     Outputs tensor of inertia information along the
                        sequence
    center-orient       Centers the molecule at the center of mass and
                        orients it along inertia axes
    set-b-per-chain     Sets B-factors to chain average

Program will also print various information extracted from the input 
PDB file. Output is defined by -p option.  Currently supported choices 
are:

    bvalue              Prints the list of average per-residue Bfactors, 
                        including all atoms, backbone and side chain 
                        columns.
    chains              Prints the list of chains with number of atoms 
                        in each and average B factor.
    phipsi              Prints the list of backbone torsions.
    bcontrast           Calculate the average B-factor of a particular
                        residue range (must be specified with --ranges
                        option) and that of its immediate environment,
                        defined as all the non-water atoms within 
                        cutoff distance
    hbonds              Prints the list of hydrogen bonds
    range_bs            Prints the average B-factor for a selected range
    rgyration           Prints radius of gyration for a selected range

--------------------------------------------------------------------------------

'''
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                            description=headerhelp)
    parser.add_argument('inpath',
                        help='The input PDB file.')
    parser.add_argument('outpath', nargs='?',
                        default='pdbtrickout.pdb',
                        help='The output PDB file.')
    parser.add_argument('-a', '--action', action='append',
                        choices = [	'extract-chains',
                                    'rename-chains',
                                    'rjust-resid', 
                                    'extract-ranges',
                                    'tinertia-ranges',
                                    'tinertia-slider',
                                    'center-orient',
                                    'set-b-per-chain',
                                    'strip-waters'],
                        default = [],
                        metavar = '', help='Action to perform')
    parser.add_argument('-p', '--outprint', action='append',
                        choices = [ 'bvalue', 
                                    'chains', 
                                    'phipsi', 
                                    'bcontrast', 
                                    'resgem', 
                                    'range_bs',
                                    'rgyration'],
                        default = [],
                        metavar='', help='Information to print out.')
    parser.add_argument('--resid', 
                        help='Residue id parameter for various commands.')
    parser.add_argument('--bvalue-print',
                        action='store_true',
                        help='Print per-residue B-balues.')
    parser.add_argument('--chids',
                        help='Chain IDs for various selections.')
    parser.add_argument('--extract-chains',
                        action='store_true',
                        help='Extract specified chains.')
    parser.add_argument('--rjust-resid', 
                        action='store_true',
                        help='Make sure resids are right-justified.')
    parser.add_argument('--ranges',
                        dest='ranges',
                        help='Residue range selection, chain, start, end. \
                              Comma-separated, no spaces. Separate chains with \
                              forward slash. Example: \
                              A,50-55,72-80/B,50-55 will select residues \
                              50 to 55 and 72 to 80 in chain A, and also \
                              residues 50 to 55 in chain B.')
    parser.add_argument('--window-size', type=int, default=5,
                        help='Sliding sequence window size.')
    parser.add_argument('--rcutoff', type=float, default=4.0,
                        help='Distance cutoff, defaults to 4A')
    args = parser.parse_args()

    from pdbtool import ReadPDBfile as read_pdb_file
    from helper import range_check, parse_ranges
    from scipy import array

    model = read_pdb_file(args.inpath)
    if model is None:
        sys.exit('Failed to read model from '+args.inpath)

    args.ranges = parse_ranges(args.ranges, model.GetChains()) if args.ranges else None


    for action in args.action:
        epyactions.__getattribute__(action.replace('-','_').lower())(args,  model)
    for whatoprint in args.outprint:
        epyactions.__getattribute__('print_'+whatoprint.replace('-','_').lower())(args,  model)

if __name__ == "__main__":
    main()
