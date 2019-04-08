''' 
    MRC module implements I/O, manipulation and analysis of the 
    MRC map file format.
'''
from binarrays import read_byte, read_integer, read_float, read_char, read_shortint, read_cints, read_fints, read_shortuint
from binarrays import write_char, write_integer, write_float, write_byte, write_shortint, write_cint, write_fint, write_shortuint
from scipy import array

mode_read = {
                0   : read_char,
                1   : read_shortint,
                2   : read_float,
                3   : read_cints,
                4   : read_fints,
                6   : read_shortuint,
            }

mode_write = {
                0   : write_char,
                1   : write_shortint,
                2   : write_float,
                3   : write_cint,
                4   : write_fint,
                6   : write_shortuint,
            }

class mrc(object):
    def __init__(self, fname=None):
        if fname is not None:
            self.read(fname)

    def read(self, fname):
        with open(fname,'rb') as fin:
            self.nx, self.ny, self.nz = read_integer(fin, 3)
            self.mode = read_integer(fin)
            self.nxstart, self.nystart, self.nzstart = read_integer(fin, 3)
            self.mx, self.my, self.mz = read_integer(fin, 3)
            self.ucell = unit_cell(cell=read_float(fin, 6))
            self.mapc, self.mapr, self.maps = read_integer(fin, 3)
            self.dmin, self.dmax, self.dmean = read_float(fin, 3)
            self.ispg, self.nsymbt = read_integer(fin, 2)
            self.extra = read_byte(fin, 100)
            self.origin = read_integer(fin, 3)
            self.kwrd = read_byte(fin, 4)
            self.machst = read_byte(fin, 4) # Not handled yet
            self.rms = read_float(fin)
            self.nlabl = read_integer(fin)
            labels = read_byte(fin, 800)
            self.labels = [labels[80*i:80*(i+1)] for i in range(self.nlabl)]
            symops = read_byte(fin, self.nsymbt)
            self.symops = [symops[80*i:80*(i+1)] for i in range(int(self.nsymbt/80))]
            self.data = array([[mode_read[self.mode](fin,self.nx) for row in range(self.ny)] for section in range(self.nz)])

    def write(self, fname):
        with open(fname,'wb') as fout:
            write_integer(fout, (self.nx, self.ny, self.nz))
            write_integer(fout, self.mode)
            write_integer(fout, (self.nxstart, self.nystart, self.nzstart))
            write_integer(fout, (self.mx, self.my, self.mz))
            write_float(fout, self.ucell.tolist())
            write_integer(fout, (self.mapc, self.mapr, self.maps))
            self.update_stats()
            write_float(fout, (self.dmin, self.dmax, self.dmean))
            write_integer(fout, (self.ispg, self.nsymbt))
            write_byte(fout, self.extra)
            write_integer(fout, self.origin)
            write_byte(fout, self.kwrd)
            write_byte(fout, self.machst)
            write_float(fout,  self.rms)
            write_integer(fout, self.nlabl)
            for label in self.labels:
                write_byte(fout, label)
            write_byte(fout, [0]*80*(10-self.nlabl))
            for symop in self.symops:
                write_byte(fout, symop)
            mode_write[self.mode](fout, self.data.flat)
            
    def update_stats(self):
        self.dmin = self.data.min()
        self.dmax = self.data.max()
        self.dmean = self.data.mean()
        self.rms = self.data.std()

class unit_cell(object):
    def __init__(self, *args, **kwds):
        if 'cell' in kwds:
            (self.a, self.b, self.c,
            self.alpha, self.beta, self.gamma) = kwds['cell']
    def tolist(self):
        return [self.a, self.b, self.c, self.alpha, self.beta, self.gamma]
