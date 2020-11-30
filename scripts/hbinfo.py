#! /usr/bin/env python3

headerhelp = \
''' 
    Returns list of hydrogen bond types defined in aconts module.
    It still prints some spurious class info, just ignore it.
'''
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import aconts

print('\n'.join(["%s\n\t%s" % (x.__name__,x.__doc__) for x in aconts.HydrogenBond.__subclasses__()]))
