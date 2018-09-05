import os, sys, subprocess
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class TotalClashScoreReader(object):
    ''' Reads the single value line produced by "phenix.clashscore verbose=False" '''
    def __init__(self, line=False, *args, **kwds):
        self.readline(line)
    def readline(self, line=False):
        if line:
            self.clashscore = float(line)
        else:
            self.clashscore = None
    def report(self):
        return "%.2f" % (self.clashscore)
    def _extra_report(self):
        return ''

from pdbminer import pdbase
class totalclashscore_pdbase(pdbase):
    tables = [('clashscore', TotalClashScoreReader(), '')]
    def process_code(self, code, fpath):
        phenix_cmd = "phenix.clashscore verbose=False "+fpath
        exitcode, output = subprocess.getstatusoutput(full_cmd)
        if not exitcode:
            self.insert_new_code(code)
            self.insert_new_item('clashscore',code,TotalClashScoreReader(output))
            self.code_lock(code)
        else:
            print("ERROR: pnenix.clashscore failed with the following output:")
            print(output)

