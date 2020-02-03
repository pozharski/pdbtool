'''
CIF-file reading module.
'''

import re
from collections import OrderedDict

class text_reader:
    '''
    General text file reader class.  Uses regular expressions to read
    selected lines or blocks of lines.
    '''
    def __init__(self, src):
        '''
        Input must be a text stream opened with open method in text mode.
        '''
        self.src = src
        back_pos = src.tell()
        self.flen = src.seek(0,2)
        src.seek(back_pos)
    def eof(self):
        return self.src.tell() == self.flen
    def close(self):
        src.close()
    def reset(self):
        self.src.seek(0)
    def find_and_read(self, ptrn):
        mptrn = re.compile(ptrn)
        while True:
            line = self.src.readline()
            retval = mptrn.match(line)
            if retval:
                return retval.groups() if mptrn.groups else True
            elif not line:
                return False
    def read_until(self, ptrn, skip_back=False):
        mptrn = re.compile(ptrn)
        retlines = []
        while True:
            if skip_back:
                back_pos = self.src.tell()
            line = self.src.readline()
            if mptrn.match(line):
                if skip_back:
                    sel.src.seek(back_pos)
                return retlines
            if line:
                retlines.append(line)
            else:
                return False
    def read_ahead(self, n=1):
        back_pos = self.src.tell()
        retval = [self.src.readline() for i in range(n)]
        self.src.seek(back_pos)
        return retval
    def skip_until(self, ptrn, skip_back=False):
        mptrn = re.compile(ptrn)
        while True:
            if skip_back:
                back_pos = self.src.tell()
            line = self.src.readline()
            if mptrn.match(line):
                if skip_back:
                    sel.src.seek(back_pos)
                return True
            if not line:
                return False
    
class cif_reader(text_reader):
    '''
    CIF file reader class.  Reads and processes data blocks.
    '''
    def loop_check(self):
        return self.read_ahead()[0].strip() == 'loop_'
    def read_block(self):
        lines = self.read_until('#')
        if not lines:
            return lines 
        mptrn=re.compile("_([^\.]*)\.([^\s]*)\s*(.*)")
        ret_dict = OrderedDict()
        for line in lines:
            readout = mptrn.match(line)
            if readout:
                bname, key, value = readout.groups()
                ret_dict[key] = value.strip().strip("'")
            else:
                ret_dict[key] += line.strip().strip("'")
        return bname, ret_dict
    def read_loop(self):
        self.skip_until('loop_')
        lines = self.read_until('#')
        if lines is False:
            return False
        fptrn=re.compile("_([^\.]*)\.([^\s]*)")
        kptrn = re.compile("('[^']*')+?")
        keys, values = [], []
        for line in lines:
            field_read = fptrn.match(line)
            if field_read:
                bname, key = field_read.groups()
                keys.append(key)
            else:
                values.append(sum([[x[1:-1]] if x[0]=="'" else x.split() for x in kptrn.split(line) if x],[]))
        values, ret_dict = list(zip(*values)), OrderedDict()
        for (i, key) in enumerate(keys):
            ret_dict[key] = values[i]
        return bname, ret_dict

class pdbentry:
    '''
    PDB entry class.  Reads a CIF file and stores data as dictionary 
        of dictionaries.
    '''
    def __init__(self, src):
        reader = cif_reader(src)
        self.code = reader.find_and_read("data_(.*)")[0]
        self.datablocks, self.arrays = OrderedDict(), OrderedDict()
        while True:
            if reader.eof():
                break
            if reader.loop_check():
                key, value = reader.read_loop()
                self.arrays[key] = value
                #print('----------')
                #print(key)
                #print('\n'.join([str(k)+": "+'|'.join(v) for k,v in self.arrays[key].items()]))
            else:
                readout = reader.read_block()
                if readout:
                    self.datablocks[readout[0]] = readout[1]
                    #print('----------')
                    #print(readout[0])
                    #print('\n'.join(["%s: %s" % (k,v) for k,v in self.datablocks[readout[0]].items()]))
        
if __name__ == '__main__':
    import sys
    with open(sys.argv[1]) as fin:
        x = pdbentry(fin)
        print(x.datablocks.keys())
