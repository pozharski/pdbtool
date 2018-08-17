#! /usr/bin/env python3

class RamaReader(object):
    ''' Reads the lines produced by phenix.ramalyze '''
    def __init__(self, line=False, *args, **kwds):
        self.readline(line)
    def readline(self, line=False):
        if line:
            chunks = line.split(':')
            self.chain = chunks[0][1]
            self.resn = int(chunks[0][2:6])
            self.ac = chunks[0][6]
            self.resi = chunks[0][-3:]
            self.score, self.phi, self.psi = [float(x) for x in chunks[1:4]]
            self.evaluation, self.type = [x for x in chunks[4:]]
        else:
            chunks = None
            self.chain, self.ac, self.resi, self.evaluation, self.type = ['']*5
            self.resn, self.score, self.phi, self.psi = [None]*4
        self._extra_reads(chunks)
    def _extra_reads(self, chunks):
        pass
    def report(self):
        return "%2s%4d%s%s:%.2f:%.2f:%.2f:%s:%s" % (self.chain,self.resn,self.ac,self.resi,self.score,self.phi,self.psi,self.evaluation,self.type)
    def _extra_report(self):
        return ''

from ..pdbminer import pdbase
class rama_pdbase(pdbase):
    tables = [('rama', RamaReader(), '')]
