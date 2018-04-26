'''
    Secondary structure record parser.  Recognizes HELIX and SHEET
    records.
'''

helix_class = {
                1 :    "Right-handed alpha",
                2 :    "Right-handed omega",
                3 :    "Right-handed pi",
                4 :    "Right-handed gamma",
                5 :    "Right-handed 3 - 10",
                6 :    "Left-handed alpha",
                7 :    "Left-handed omega",
                8 :    "Left-handed gamma",
                9 :    "2 - 7 ribbon/helix",
                10:    "Polyproline"
            }
strand_sense = {
                -1 : "anti-parallel",
                 0 : "first",
                 1 : "parallel"
                }

import re

class pdbr_ss:
    def __init__(self, line):
        if re.match('HELIX ', line):
            pdbr_helix.__init__(self, line)
        elif re.match('SHEET ', line):
            pdbr_strand.__init__(self, line)
        else:
            raise ValueError('Not a HELIX/SHEET record:\n'+line)
    def is_one_chain(self):
        return self.initChainID == self.endChainID
    def no_icodes(self):
        return not bool((self.initICode+self.endICode).strip())
    def is_incremental(self):
        return self.endSeqNum>self.initSeqNum
    def can_range(self):
        return self.is_one_chain() and self.no_icodes() and self.is_incremental()
    def residue_range(self):
        if self.can_range():
            return {self.initChainID: [(self.initSeqNum,self.endSeqNum)]}
    def extract(self, model):
        if self.can_range():
            return model.extract_range(self.residue_range())

class pdbr_helix(pdbr_ss):
    def __init__(self, line):
        if re.match('HELIX ', line):
            try:
                self.serNum = int(line[7:10])
                self.helixID = line[11:14]
                self.initResName = line[15:18]
                self.initChainID = line[19]
                self.initSeqNum = int(line[21:25])
                self.initICode = line[25]
                self.endResName = line[27:30]
                self.endChainID = line[31]
                self.endSeqNum = int(line[33:37])
                self.endICode = line[37]
                self.helixClass = int(line[38:40])
                self.comment = line[40:70]
                self.length = int(line[71:76])
            except:
                print('\nFault: '+line.strip())
                raise
        else:
            raise ValueError('Not a HELIX record:\n'+line)
    def helix_class(self):
        return helix_class.get(self.helixClass)

class pdbr_strand(pdbr_ss):
    def __init__(self, line):
        if re.match('SHEET ', line):
            try:
                self.strand = int(line[7:10])
                self.sheetID = line[11:14]
                self.numStrands = int(line[14:16])
                self.initResName = line[17:20]
                self.initChainID = line[21]
                self.initSeqNum = int(line[22:26])
                self.initICode = line[26]
                self.endResName = line[28:31]
                self.endChainID = line[32]
                self.endSeqNum = int(line[33:37])
                self.endICode = line[37]
                try:
                    self.sense = int(line[38:40])
                except ValueError:
                    self.sense = 0
                if self.sense and line[41:].strip():
                    self.curAtom = line[41:45]
                    self.curResName = line[45:48]
                    self.curChainId = line[49]
                    self.curResSeq = int(line[50:54])
                    self.curICode = line[54]
                    self.prevAtom = line[56:60]
                    self.prevResName = line[60:63]
                    self.prevChainId = line[64]
                    self.prevResSeq = int(line[65:69])
                    self.prevICode = line[69]
            except:
                print('\nFault: '+line.strip())
                raise
        else:
            raise ValueError('Not a SHEET record:\n'+line)
    def strand_sense(self):
        return strand_sense.get(self.sense)
    def pdbrecord(self):
        return 'SHEET %3d %4s%2d %3s %1s %4d%1s  %3s %1s %4d%1s %2d' % (self.strand,self.sheetID,self.numStrands,self.initResName,self.initChainID,self.initSeqNum,self.initICode,self.endResName,self.endChainID,self.endSeqNum,self.endICode,self.sense)

class ssrecords:
    def __init__(self, path):
        with open(path) as f:
            lines = [x for x in f.readlines() if re.match('(HELIX )|(SHEET )', x)]
        self.helices = [pdbr_helix(x) for x in [x for x in lines if re.match('HELIX ', x)]]
        self.strands = [pdbr_strand(x) for x in [x for x in lines if re.match('SHEET ', x)]]
        self.quarantine = []
    def get_helices(self, model):
        return [x.extract(model) for x in self.helices]
    def get_strands(self, model):
        return [x.extract(model) for x in self.strands]
    def sanitize(self):
        self.quarantine += [x for x in self.helices if not x.can_range()]
        self.quarantine += [x for x in self.strands if not x.can_range()]
        self.helices = [x for x in self.helices if x.can_range()]
        self.strands = [x for x in self.strands if x.can_range()]
        
