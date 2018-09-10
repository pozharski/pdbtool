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
            chunks = [int(x) for x in chunks[:2]]+[chunks[2].strip()]+[float(x) for x in chunks[4:7]]
        self.mmsize, self.symnum, self.spacegroup, self.asa, self.bsa, self.deltaG = chunks
    def report(self):
        return "%6d%4d%19s%11.1f%11.1f%11.1f" % (self.mmsize, self.symnum, self.spacegroup, self.asa, self.bsa, self.deltaG)

class pisa_dbsres_pdbase(pdbminer.pdbase):
    tables = [('pisa_dbsres', pisa_dbsres_reader(), '')]
    def process_codes(self, fpath):
        with open(fpath) as fin:
            for line in fin.readlines():
                m = dbsres_ptrn.match(line)
                if m:
                    code = m.groups()[0]
                    self.insert_new_code(code)
                    self.insert_new_item('pisa_dbsres', code, pisa_dbsres_reader(line))
                    self.code_lock(code)
