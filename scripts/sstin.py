#! /usr/bin/env python
headerhelp = \
'''
SSTIN runs the survey of inertia ellipsoids of secondary structure
        elements.
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('-p', '---pdb', 
    help='Input pdb file.  This parameter is either a PDB (not case-sensitive) code or path to a PDB file.')
parser.add_argument('-l', '--list-file',
    help="List of PDB codes to process.  Results are stored in SQLite database.")
parser.add_argument('-d', '--sqlite-file', default='sstin.sqlite',
    help="SQLite database to store/retrieve processed results.")
args = parser.parse_args()

import os, sys, sqlite3

from pdbtool import get_model, get_ssrecords
from sstin_sqlite import sstin_dbase

if args.pdb:
    model = get_model(args.pdb).Extract('p_bb')
    if model:
        ssx = get_ssrecords(args.pdb, sanitize=True)
        helices = ssx.get_helices(model)
        htins = map(lambda x : x.GetInertiaTensor(), helices)
        strands = ssx.get_strands(model)
        blens = map(lambda x : x.GetProteinResidueNumber(), strands)
        btins = map(lambda x : x.GetInertiaTensor(), strands)
        for h, t in filter(lambda x : x[1].valid, zip(*[ssx.helices, htins])):
            print "Helix #%-3d: %5.2f %5.2f %5.2f %4d residues long, %s" % tuple([h.serNum] + t.eccent + [h.length, h.helix_class()])
        for b, length, t in filter(lambda x : x[2].valid, zip(*[ssx.strands, blens, btins])):
            print "Sheet %4s Strand #%-3d: %5.2f %5.2f %5.2f %4d residues long, %s" % tuple([b.sheetID, b.strand] + t.eccent + [length, b.strand_sense()])
if args.list_file:
    if os.access(args.list_file, os. R_OK):
        with open(args.list_file) as flist:
            codes = map(lambda x : x.split()[0], flist.readlines())
        dbase = sstin_dbase(args.sqlite_file)
        codes = dbase.filter_codes(codes)
        for code in codes:
            sys.stdout.write("Processing "+code+"... ")
            dbase.insert_new_code(code)
            model = get_model(code).Extract('p_bb')
            if model:
                ssx = get_ssrecords(code, sanitize=True)
                helices = ssx.get_helices(model)
                htins = map(lambda x : x.GetInertiaTensor(), helices)
                strands = ssx.get_strands(model)
                blens = map(lambda x : x.GetProteinResidueNumber(), strands)
                btins = map(lambda x : x.GetInertiaTensor(), strands)
                for h, t in filter(lambda x : x[1].valid, zip(*[ssx.helices, htins])):
                    dbase.insert_helix(code, h, t)
                for b, length, t in filter(lambda x : x[2].valid, zip(*[ssx.strands, blens, btins])):
                    dbase.insert_strand(code, b, t, length)
                dbase.code_lock(code)
                dbase.commit()
            sys.stdout.write("done\n")
        dbase.close()
