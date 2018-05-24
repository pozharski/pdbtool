#! /usr/bin/env python3
headerhelp = \
'''
HBLIST lists the hydrogen bonds in a structure.
'''
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('model',
                    help='Model PDB file.')
parser.add_argument('--ohmax', type=float, default=2.5,
					help='O...H distance cutoff for main chain hydrogen bonds')
parser.add_argument('--nomax', type=float, default=3.2,
					help='N...O distance cutoff for main chain hydrogen bonds')
parser.add_argument('--dhamin', type=float, default=90.0,
					help='Minimum donor...hydrogen...acceptor angle in hydrogen bonds')
parser.add_argument('--pcutoff', type=float, default=0.01,
					help='p-value cutoff for including hydrogen bonds')
parser.add_argument('-b', '--bondtype', default='all',
                    help='Comma-separated list of bond types you wish to show.  Defaults to "all"')
parser.add_argument('--ranges',
                    dest='ranges',
                    help='Residue range selection, chain, start, end. \
                          Comma-separated, no spaces. Separate chains with \
                          forward slash. Example: \
                          A,50-55,72-80/B,50-55 will select residues \
                          50 to 55 and 72 to 80 in chain A, and also \
                          residues 50 to 55 in chain B.')
args = parser.parse_args()

from pdbtool import ReadPDBfile as read_pdb_file
import aconts

model = read_pdb_file(args.model)
if model is None:
    sys.exit('Failed to read the model coordinates from '+args.model)

if args.ranges:
    ranges = dict([tuple([r[0], [tuple([int(x) for x in s.split('-')]) for s in r.split(',')[1:]]]) for r in args.ranges.split('/')])
    try:
        uniranges = ranges.pop('*')
        for chid in model.GetChains():
            ranges[chid] = ranges.get(chid, [])
            ranges[chid].extend(uniranges)
    except KeyError:
        pass
else:
    ranges = None

if ranges is not None:
    resids = model.range_residues(ranges)

if args.bondtype.lower() == 'all' or args.bondtype.lower() == 'backbone':
    residnames = model.get_residnames()
    bb = model.backbone()
    hb = bb.MChbonds(args.ohmax, args.nomax, args.dhamin)
    if ranges is not None:
        hb = dict([(k,v) for k,v in hb.items() if k[:6] in resids or k[6:] in resids])
    print("--------------------------------------------------------------------------------")
    print("Main chain hydrogen bonds (%d)" % len(hb))
    print("Residue1     Residue2       D(O-H)  D(N-O)  DHA")
    for k,v in hb.items():
        print("%3s%7s N %3s%7s O %6.2f  %6.2f  %6.1f " % tuple([residnames.get(k[:6]), k[:6], residnames.get(k[6:]), k[6:]]+v))

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
    if key in dir(aconts):
        print("--------------------------------------------------------------------------------")
        HydroBonds = eval('aconts.'+key)
        hb = HydroBonds(model, rcutoff=min(3.2,aconts.pcutoff(0.01, key)))
        hb.pfilter(args.pcutoff)
        if ranges is not None:
            hb.resid_filter(resids)
        print(HydroBonds.__doc__)
        hb.report()
