from .download import pathcode_check
from . import pdb_parser
from .pdbtool import ReadPDBfile

def get_model(pathcode, softcheck=True, active=True):
    '''
        Return pdbmolecule object.  If pathcode is an existing PDB file,
        it will be imported.  Otherwise, pathcode will be expected to
        be a PDB ID.  Local folder ./pdb-download will be checked and
        if file is not found, download will be attempted from the 
        Protein Data Bank.
    '''
    fpath = pathcode_check(pathcode, softcheck, active, return_stream=True)
    if fpath:
        model = ReadPDBfile(fpath)
        fpath.close()
        return model
    else:
        return None

def get_ssrecords(pathcode, sanitize=False):
    '''
        Return ssrecords object.  Rules for finding/retrieving the 
        data file are the same as in get_model.  
    '''
    fpath = pathcode_check(pathcode)
    if fpath:
        ssx = pdb_parser.get_ssrecords(fpath)
        if sanitize:
            ssx.sanitize()
        return ssx
    else:
        return None
