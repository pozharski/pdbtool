#! /usr/bin/env python3

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

--------------------------------------------------------------------------------

'''

import sys

from argparse import ArgumentParser, RawDescriptionHelpFormatter
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
                                'center-orient',],
                    default = [],
                    metavar = '', help='Action to perform')
parser.add_argument('-p', '--outprint', action='append',
                    choices = ['bvalue', 'chains', 'phipsi', 'bcontrast', 'resgem'],
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

ranges = parse_ranges(args.ranges, model.GetChains()) if args.ranges else None

for whatoprint in args.outprint:
    if whatoprint == 'bvalue':
        resnames = model.GetResidueNames()
        resids, b0, b1, b2 = model.GetResidueBvectorByChain()
        for chid in resids:
            if not args.chids or chid in args.chids:
                for (i, resnum) in enumerate(resids[chid]):
                    print('%3s %5s %6.2f %6.2f %6.2f' % (resnames[chid+resnum],\
                                                        chid+resnum.strip(),\
                                                        b0[chid][i],\
                                                        b1[chid][i],\
                                                        b2[chid][i]))
    elif whatoprint == 'chains':
        chanums = model.GetChains()
        bavs = model.GetChainAverageBfactor()
        print('Chain   Atoms   <B>')
        for key in chanums:
            print(key + str(chanums[key]).rjust(12) + " %6.2f" % bavs[key])
    elif whatoprint == 'phipsi':
        phi, psi = model.PhiPsiList()
        for resid in model.resid_lister('p_bb'):
            print("%6s %8.2f %8.2f" % (resid, phi.get(resid, float('nan')), psi .get(resid, float('nan'))))
    elif whatoprint == 'bcontrast':
        resids = [y for y in set([x.resid() for x in model.atoms]) if array([range_check(int(y[1:-1]),x) for x in ranges.get(y[0],[])]).any()]
        listik = model.atom_lister('resids', resids=resids)
        print('%d atoms selected for analysis\n<B>sel = %.2f' % (len(listik), model.GetAverageBfactor(listik=listik)))
        vicatoms = model.atom_lister('vicinity', model.atom_lister('notwater'), corelist=listik, rcutoff=args.rcutoff)
        print('%d non-water atoms found within %.2f Angstroms\n<B>vic = %.2f' % (len(vicatoms), args.rcutoff, model.GetAverageBfactor(listik=vicatoms)))
    elif whatoprint == 'hbonds':
        pass
    elif whatoprint == 'resgem':
        residue = model.get_residues('resid',resid=args.resid)[args.resid]
        b,a,t,m = residue.BondsAnglesTorsions(printout=True)
        print('Residue %s (%s)' % (args.resid, residue.get_res_name()))
        print('----------- Bonds ----------- ')
        for x in b:
            print(x)
        print('----------- Angles ----------- ')
        for x in a:
            print(x)
        print('----------- Torsions ----------- ')
        for x in t:
            print(x)
        print('----------- Impropers ----------- ')
        for x in m:
            print(x)
        
    
if len(args.action) > 1:
    sys.exit('Multiple action not yet supported.  One at a time, please.')

for whatodo in args.action:
    if whatodo == 'extract-chains':
        if args.chids:
            model.writePDBchains(args.outpath, args.chids)
        else:
            print('This does not compute - extract chains but no chains listed?')
    elif whatodo == 'rjust-resid':
        model.rjust_res_names()
        model.writePDB(args.outpath)
    elif whatodo == 'extract-ranges':
        if ranges:
            model.extract_range(ranges).writePDB(args.outpath)
        else:
            print('This does not compute - extract but no ranges listed?')
    elif whatodo == 'tinertia-ranges':
        if ranges:
            ramodel = model.extract_range(ranges)
            tin = ramodel.GetInertiaTensor(ramodel.atom_lister('backbone'))
        else:
            tin = model.GetInertiaTensor(model.atom_lister('backbone'))
        print(tin.report())
        print(tin.frame())
    elif whatodo == 'tinertia-slider':
        pass
    elif whatodo == 'center-orient':
        remodel = model.copy()
        remodel.shift(-remodel.GetCoM())
        tin_remodel = remodel.GetInertiaTensor()
        remodel.transform(tin_remodel.vr)
        remodel.writePDB(args.outpath)
    elif whatodo == 'rename-chains':
        if args.chids:
            for chid1,chid2 in [(x[0],x[1]) for x in args.chids.split(',')]:
                model.rename_chain(chid1,chid2)
            model.writePDB(args.outpath)
        else:
            print('This does not compute - rename chains but no chain pairs listed?')
