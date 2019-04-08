import urllib.request, os, sys, gzip

mirror_folder = os.environ.get('PDBMIRROR','pdb-mirror')
download_folder = os.environ.get('PDBDOWNLOAD','pdb-download')

def combo_open(fpath, return_stream=False, compressed=False):
    if return_stream:
        return gzip.open(fpath,'r') if compressed else open(fpath)
    else:
        return fpath

def pathcode_check(value, softcheck=True, active=True, return_stream=False):
    if os.access(value, os.R_OK):
        return combo_open(value, return_stream=return_stream)
    if os.path.isdir(mirror_folder):
        mirror_found = True
    else:
        mirror_found = False
        if not os.path.isdir(download_folder):
            os.mkdir(download_folder)
            sys.stderr('Standard download folder created at %s\n', os.path.abspath(download_folder))
    if mirror_found:
        fpath = os.path.join(mirror_folder,value[1:3],'pdb'+value+os.extsep+'ent'+os.extsep+'gz')
    else:
        fpath = os.path.join(download_folder, value.lower()+os.extsep+'pdb')
    if os.access(fpath, os.R_OK):
        return combo_open(fpath, return_stream=return_stream, compressed=mirror_found)
    elif mirror_found:
        msg = "File "+fpath+" not found.\nTry re-syncing the mirror."
        if softcheck:
            sys.stderr.write(msg)
            return None
        else:
            sys.exit(msg)
    else:
        if active:
            pdbfile = urllib.request.urlretrieve('http://www.rcsb.org/pdb/files/'+value.lower()+'.pdb.gz')
            if pdbfile[1].type == 'text/html':
                sys.stderr.write('Failed to retrieve PDB entry '+value+'.  You need to rethink your life.\n')
                return None
            with open(fpath, 'w') as fout:
                fin = gzip.open(pdbfile[0])
                for line in fin:
                    fout.write(line)
                fin.close()
            return combo_open(fpath, return_stream=return_stream)
        else:
            msg = "File "+fpath+" not found.\nTry re-syncing the mirror."
        if softcheck:
            sys.stderr.write(msg)
            return None
        else:
            sys.exit(msg)
