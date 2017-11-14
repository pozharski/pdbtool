from download import pathcode_check
import pdb_parser
from pdbtool import ReadPDBfile

def get_model(pathcode):
    fpath = pathcode_check(pathcode)
    if fpath:
        return ReadPDBfile(fpath)
    else:
        return None

def get_ssrecords(pathcode, sanitize=False):
    fpath = pathcode_check(pathcode)
    if fpath:
        ssx = pdb_parser.get_ssrecords(fpath)
        if sanitize:
            ssx.sanitize()
        return ssx
    else:
        return None
