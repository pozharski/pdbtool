#! /usr/bin/env python3

headerhelp = \
'''
TINCOMP analyzes structural changes using tensor of inertia changes
throughout protein structure.
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('model1',
                    help='First model PDB file.')
parser.add_argument('model2',
                    help='Second model PDB file.')
parser.add_argument('segfile',
                    help='Input file with structural segments.')
parser.add_argument('--simple-output',
                    action='store_true',
                    help='Formatted output for easy parsing.')
args = parser.parse_args()
from .pdbtool import ReadPDBfile as read_pdb_file
from .pdbtool import pdbmolecule
model1 = read_pdb_file(args.model1)
model2 = read_pdb_file(args.model2)
with open(args.segfile) as segf:
    segrs = {}
    codes = []
    for line in segf:
        try:
            if line.strip()[0] not in '#;!%':
                key = line.split()[0].lower()
                if key == 'range':
                    code = line.split()[1]
                    segrs[code]={}
                    for segr in line.split()[2:]:
                        chid = segr[0]
                        if chid in segrs[code]:
                            segrs[code][chid].append(tuple(map(int,segr[1:].split('-'))))
                        else:
                            segrs[code][chid] = [tuple(map(int,segr[1:].split('-')))]
                    codes.append(code)
        except:
            if line.strip():
                print("Warning: the following line in the segment file was uninterpretable:")
                print('|'+line.strip()+'|')
model1bb = pdbmolecule(atoms=model1.atom_getter('backbone'))
model2bb = pdbmolecule(atoms=model2.atom_getter('backbone'))
tin1 = {}
tin2 = {}
for key in codes:
    value = segrs[key]
    tin1[key] = model1bb.extract_range(value).GetInertiaTensor()
    tin2[key] = model2bb.extract_range(value).GetInertiaTensor()
    print('SEGMENT '+key)
    if args.simple_output:
        print(tin1[key].numcomp(tin2[key], key))
    else:
        print(tin1[key].reportcomp(tin2[key]))
print('Segment pairwise shifts/angles')
for (i1,key1) in enumerate(codes):
    for key2 in codes[i1+1:]:
            print(' %10s - %10s : %6.1f %6.1f vs %6.1f %6.1f' % tuple([key1, key2, tin1[key1].shift(tin1[key2]), tin1[key1].angles(tin1[key2])[0], tin2[key1].shift(tin2[key2]), tin2[key1].angles(tin2[key2])[0]]))
