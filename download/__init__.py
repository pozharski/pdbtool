download_folder = 'pdb-download'

import urllib, os, sys, gzip

def pathcode_check(value):
    if os.access(value, os.R_OK):
        return value
    if not os.path.isdir(download_folder):
        os.mkdir(download_folder)
        sys.stderr('Standard download folder created at %s\n', os.path.abspath(download_folder))
    fpath = os.path.join(download_folder, value.lower()+os.extsep+'pdb')
    if os.access(fpath, os.R_OK):
        return fpath
    pdbfile = urllib.urlretrieve('http://www.rcsb.org/pdb/files/'+value.lower()+'.pdb.gz')
    if pdbfile[1].type == 'text/html':
        sys.stderr.write('Failed to retrieve PDB entry '+value+'.  You need to rethink your life.\n')
        return None
    with open(fpath, 'w') as fout:
        fin = gzip.open(pdbfile[0])
        for line in fin:
            fout.write(line)
        fin.close()
    return fpath
