#! /usr/bin/env python3

headerhelp = \
'''
HBNCOMP analyzes changes in hydrogen bonding patterns.
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('model1',
                    help='First model PDB file.')
parser.add_argument('model2',
                    help='Second model PDB file.')
parser.add_argument('--ohmax', type=float, default=2.5,
					help='O...H distance cutoff for main chain hydrogen bonds')
parser.add_argument('--nomax', type=float, default=3.2,
					help='N...O distance cutoff for main chain hydrogen bonds')
parser.add_argument('--dhamin', type=float, default=90.0,
					help='Minimum donor...hydrogen...acceptor angle in hydrogen bonds')
parser.add_argument('--pcutoff', type=float, default=0.01,
					help='p-value cutoff for including hydrogen bonds')
parser.add_argument('--bondtype', default='all',
                                        help='Comma-separated list of bond types you wish to show.  Defaults to "all"')
args = parser.parse_args()
from .pdbtool import ReadPDBfile as read_pdb_file
from .pdbtool import pdbmolecule
import aconts, sys
model1 = read_pdb_file(args.model1)
if model1 is None:
    sys.exit('Failed to read the model coordinates from '+args.model1)
model2 = read_pdb_file(args.model2)
if model2 is None:
    sys.exit('Failed to read the model coordinates from '+args.model2)
if args.bondtype.lower() == 'all':
    bb1 = model1.backbone()
    bb2 = model2.backbone()
    hb1 = bb1.MChbonds(args.ohmax, args.nomax, args.dhamin)
    hb2 = bb2.MChbonds(args.ohmax, args.nomax, args.dhamin)
    print("--------------------------------------------------------------------------------")
    print("Conserved main chain hydrogen bonds")
    for hbond in sorted(list(set(hb1).intersection(hb2))):
        ticks = int(10*(hb1[hbond][1]-hb2[hbond][1]))
        print("%7s%7s %6.2f %6.2f %10s|%-10s" % (hbond[:6], hbond[6:], hb1[hbond][1], hb2[hbond][1], '*'*min(ticks,10), '*'*min(-ticks,10)))
    print("Broken mainchain hydrogen bonds")
    for hbond in sorted(list(set(hb1).difference(hb2))):
        if hbond[:6] not in bb2.protons or hbond[6:] not in bb2.residues:
            print("%7s%7s %6.2f ??????" % (hbond[:6], hbond[6:], hb1[hbond][1]))
        else:
            otherbond = bb2.get_mchbond(hbond[:6], hbond[6:])
            ticks = int(10*(hb1[hbond][1]-otherbond[1]))
            print("%7s%7s %6.2f %6.2f %10s|%-10s" % (hbond[:6], hbond[6:], hb1[hbond][1], otherbond[1], '*'*min(ticks,10), '*'*min(-ticks,10)))
    print("Newly formed mainchain hydrogen bonds")
    for hbond in sorted(list(set(hb2).difference(hb1))):
        if hbond[:6] not in bb1.protons or hbond[6:] not in bb1.residues:
            print("%7s%7s ?????? %6.2f" % (hbond[:6], hbond[6:], hb2[hbond][1]))
        else:
            otherbond = bb1.get_mchbond(hbond[:6], hbond[6:])
            ticks = int(10*(otherbond[1]-hb2[hbond][1]))
            print("%7s%7s %6.2f %6.2f %10s|%-10s" % (hbond[:6], hbond[6:], otherbond[1], hb2[hbond][1], '*'*min(ticks,10), '*'*min(-ticks,10)))
    print("--------------------------------------------------------------------------------")
print("Side chain hydrogen bonds listed by type")
if args.bondtype.lower() == 'all':
    keys = [
    'ThrOG1toBackbone','SerOGtoBackbone','TyrOHtoBackbone',
    'TrpNE1toBackbone',
    'LystoBackbone',
    'HistoBackbone', 'CarboxamideNtoBackbone',
    'TyrOHtoAnionic',
    'ThrOG1toAnionic',
    'SerOGtoAnionic',
    'TyrOHtoCarboxamideO',
    'SerOGtoCarboxamideO',
    'ThrOG1toCarboxamideO',
    'CarboxamideNtoCarboxamideO',
    'TrpNE1toAnionic',
    'TrpNE1toCarboxamideO',
    'SerOGtoSerOG',
    'CarboxamideNToAnionic',
    'ThrOG1toThrOG1',
    'TyrOHtoTyrOH',
    'SerOGtoThrOG1', 
    'ThrOG1toSerOG',
    'TyrOHtoThrOG1',
    'TrpNE1toTyrOH',
    'TrpNE1toSerOG',
        ]
else:
    keys = args.bondtype.split(',')
for key in keys:
    print("--------------------------------------------------------------------------------")
    HydroBonds = eval('aconts.'+key)
    hb1 = HydroBonds(model1, rcutoff=min(3.2,aconts.pcutoff(0.01, key)))
    hb1.pfilter(args.pcutoff)
    hb2 = HydroBonds(model2, rcutoff=min(3.2,aconts.pcutoff(0.01, key)))
    hb2.pfilter(args.pcutoff)
    print(HydroBonds.__doc__)
    hb1.report_diffs(hb2)
