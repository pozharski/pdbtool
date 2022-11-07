''' This module includes various useful classes that don't seem to belong
    anywhere else. '''

import sys, os
import array as binarray

class sock:
    '''
    Class replaces standard output with a pipe.  There is no 
    functionality to read from it (other than directly accessing
    sock.r).  To restore standard output, use sock.close().
    '''
    def __init__(self):
        '''
        Create a pipe, replace standard output with it and store the 
        previous stdout to restore upon closing the sock.
        '''
        self.r, self.w = os.pipe()
        self.restorout = sys.stdout
        sys.stdout = self

    def close(self):
        '''
        Close the sock and restore standard output.
        '''
        os.close(self.r)
        os.close(self.w)
        sys.stdout = self.restorout

    def write(self, line):
        '''
            Output a line into the sock.
        '''
        os.write(self.w, line)

class progressbar:
    ''' The text progress bar.  It can be used in lieu of the
            wx.Gauge as long as only SetRange() and SetValue() methods
            are in use. GetParent() method returns None.  The purpose
            of this is to make methods work both from command line
            and in wx-based applications. '''

    def __init__(self, maxrange=10, loud=False):
        self.range = maxrange
        if loud:
            self.SetRange(maxrange)

    def SetRange(self, maxrange):
        ''' Sets the maximum range of the progress bar. '''
        sys.stderr.write('\r'+' '*(self.range+2))
        self.range = maxrange
        sys.stderr.write('\r|'+' '*(maxrange)+'|')

    def SetValue(self, value):
        ''' Sets the value of the progress bar. '''
        sys.stderr.write('\r|'+'-'*value+' '*(self.range-value)+'|')
        if value in (0, self.range):
            sys.stderr.write('\n')

    def GetParent(self):
        ''' This method is present for the sake of compatibility and 
            returns None. '''
        return None

def readarray(fin, astype='I', num=1):
    '''
        Reads an element from a binary stream.
    '''
    var = binarray.array(astype)
    var.fromfile(fin, num)
    return var

def list_to_float(items):
    '''
        This method converts a list of elements into the list of floats.
    '''
    return [float(x) for x in items]

def list_to_int(items):
    '''
        This method converts a list of elements 
        into the list of integers.
    '''
    return [int(x) for x in items]

def list_to_strip(items):
    '''
        This method returns the list of the strings stripped.
    '''
    return [x.strip() for x in items]

def range_check(value, limits):
    return value >= limits[0] and value <= limits[1]

def parse_ranges(ranges, chids):
    ranges = dict([tuple([r[0], [tuple([int(x) for x in s.split('-')]) for s in r.split(',')[1:]]]) for r in ranges.split('/')])
    try:
        uniranges = ranges.pop('*',False)
        if uniranges:
            for chid in chids:
                ranges[chid] = ranges.get(chid, [])
                ranges[chid].extend(uniranges)
    except KeyError:
        pass
    return ranges
