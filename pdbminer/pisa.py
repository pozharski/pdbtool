#! /usr/bin/env python3

# There should be a way to do this with relative imports
import os, sys, re
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pdbminer

dbsres_ptrn = re.compile("\s*\d*\s*(\S{4})\s*(\d*)\s*(\d*)(.{19})\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)")

class pisa_dbsres_reader(pdbminer.general_reader):
    ''' Reads the single value line produced by PBDe PISA search '''
    def readline(self, line=False):
        if not line:
            chunks = [None, None, '', None, None, None]
        else:
            chunks = dbsres_ptrn.match(line).groups()[1:]
            chunks = [int(x) for x in chunks[:2]]+[chunks[2].strip()]+[float(x) for x in chunks[3:7]]
        self.mmsize, self.symnum, self.spacegroup, self.asa, self.bsa, self.deltaG = chunks
    def report(self):
        return "%6d%4d%19s%11.1f%11.1f%11.1f" % (self.mmsize, self.symnum, self.spacegroup, self.asa, self.bsa, self.deltaG)

class pisa_dbsres_pdbase(pdbminer.pdbase):
    tables = [('pisa_dbsres', pisa_dbsres_reader(), '')]
    def process_codes(self, fpath):
        with open(fpath) as fin:
            precodes = self.fetch_processed_codes()
            for line in fin.readlines():
                m = dbsres_ptrn.match(line)
                if m:
                    code = m.groups()[0]
                    if code not in precodes:
                        self.insert_new_code(code)
                        reader = pisa_dbsres_reader(line)
                        pisafolder = os.path.join(os.environ.get('PISA_PATH',os.path.join(os.getcwd(),'pisa_download')),str(reader.mmsize))
                        if not os.access(pisafolder, os.F_OK):
                            os.makedirs(pisafolder)
                        pdbpath = os.path.join(pisafolder,'pisa_'+code+'.pdb')
                        self.insert_new_item('pisa_dbsres', code, reader)
                        if os.access(pdbpath, os.F_OK):
                            print("Entry %s already downloaded" % (code))
                            retval = 1
                        else:
                            print("Downloading entry %s... " % (code), end='')
                            retval = int(self.pisa_download(code,pdbpath))
                        if retval:
                            print("success")
                        else:
                            print("failed")
                        self.code_lock(retval)
        self.commit()

    def download_check(self, fDownload=False):
        ''' Checks the database for downloaded structures. Needs standard
            PISA database output file as fpath.  Use fDownload flag to
            specify whether missing data should be donwladed. '''
        for code in self.fetch_codes():
            for item in self.get_items('pisa_dbsres', code):
                pisafolder = os.path.join(os.environ.get('PISA_PATH',os.path.join(os.getcwd(),'pisa_download')),str(item.mmsize))
                if not os.access(pisafolder, os.F_OK):
                    os.makedirs(pisafolder)
                pdbpath = os.path.join(pisafolder,'pisa_'+code+'.pdb')
                if os.access(pdbpath, os.F_OK):
                    print("Entry %s already downloaded, update status" % (code))
                    self.code_lock(code)
                else:
                    if fDownload:
                        print("Downloading entry %s... " % (code), end='')
                        retval = int(self.pisa_download(code,pdbpath))
                        if retval:
                            print("success")
                            self.code_lock(code)
                        else:
                            print("failed")
                            self.code_unlock(code)
        self.commit()

    def pisa_download(self, code, path):
        ''' Downloads a model from PISA server. '''
        from urllib.request import urlopen
        try:
            with urlopen('http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/multimer.pdb?'+code+':1,1') as source:
                with open(path, 'wb') as destination:
                    destination.write(source.read())
            return True
        except:
            sys.stderr.write("Download failed with %s exception: %s\n" % (str(sys.exc_info[0]),str(sys.exc_info[1])))
            return False

    def get_path(self, code):
        pdbpath = []
        for item in self.get_items('pisa_dbsres', code):
            pisafolder = os.path.join(os.environ.get('PISA_PATH',os.path.join(os.getcwd(),'pisa_download')),str(item.mmsize))
            pdbpath.append(os.path.join(pisafolder,'pisa_'+code+'.pdb'))
        return pdbpath
