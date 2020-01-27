#!/usr/bin/env python3
headerhelp = \
'''
DATREFINE reads a hb database that contains distance/angle/torsion 
columns and refines the DAT parameters.
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('--sqlpath', default='hbpisa.sqlite',
                    help='HB database file.  Defaults to hbpisa.sqlite.')

parser.add_argument('--torshift', type=float,
                    help='Shift the torsion angle domain')
parser.add_argument('--pvalue', type=float, default=0.05,
                    help='p-value cutoff to use for refinement.')
parser.add_argument('--hbtype', '-b',
                    help="Hydrogen bond type.")
parser.add_argument('--title', default='DAT clusters',
                    help='Title of the figure.')

args = parser.parse_args()

import os, sys, shutil
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from aconts import hbond_pdbase, DATKERNEL_PARAMS, sigmap
from scipy import array, ones, bincount, isfinite, ravel
from matplotlib.pyplot import figure, title, subplot, show, grid, xlabel, ylabel, plot
from scipy.stats import iqr
from tempfile import mkstemp

tf,tpath = mkstemp()
os.close(tf)
shutil.copyfile(args.sqlpath, tpath)
hpdb = hbond_pdbase(tpath)
hpdb.filter_pvalues(args.hbtype, args.pvalue)

distance,angle,torsion=array([(x.d,x.angle1,x.tor1) for x in [x[1] for x in hpdb.get_hbonds()]]).T.astype(float)

ind = isfinite(distance+angle+torsion)
distance, angle, torsion = distance[ind], angle[ind], torsion[ind]

if args.torshift!=None:
    torsion[torsion<args.torshift] += 360.0

figure(args.title)

from scipy.cluster.vq import kmeans2, whiten

d0,a0,t0=array(DATKERNEL_PARAMS[args.hbtype]).T[[0,2,4]].tolist()
Ndat = len(d0)
z = whiten(array([d0+distance.tolist(),a0+angle.tolist(),t0+torsion.tolist()]).T)
centroid, label = kmeans2(z[Ndat:], z[:Ndat], minit='matrix')

freqs = 100*bincount(label)/len(label)
c_d = [distance[label==i].mean() for i in range(Ndat)]
c_a = [angle[label==i].mean() for i in range(Ndat)]
c_t = [torsion[label==i].mean() for i in range(Ndat)]
s_d = sigmap([iqr(distance[label==i],scale='normal') for i in range(Ndat)], args.pvalue)
s_a = sigmap([iqr(angle[label==i],scale='normal') for i in range(Ndat)], args.pvalue)
s_t = sigmap([iqr(torsion[label==i],scale='normal') for i in range(Ndat)], args.pvalue)
for (i, c) in enumerate(array([[c_d,s_d],[c_a,s_a],[c_t,s_t]]).T):
    print("Cluster #%-2d: " % (i+1) + "%8.3f %8.3f "*len(c.T) % tuple(ravel(c.T)) + "%5.1f%%" % (freqs[i]))
    datvalues = ravel(array([c[0],2*c[1]**2]).T)
    print("               (" + ','.join(['%.3f']*len(datvalues)) % tuple(datvalues) + ")")

subplot(311)
title('Distance-Angle', position=(0.05,0.95), ha='left', va='top', size='xx-large')
for i in range(Ndat):
    plot(distance[label==i],angle[label==i],'rgbcmy'[i%6]+',',c_d[i],c_a[i],'rgbcmy'[i%6]+'*')
xlabel(r'D, $\AA$', size='xx-large')
ylabel(r'$\alpha$, $\degree$', position=(0.9,0), size='xx-large', va='center')
grid()

subplot(312)
title('Distance-Torsion', position=(0.05,0.95), ha='left', va='top', size='xx-large')
for i in range(Ndat):
    plot(distance[label==i],torsion[label==i],'rgbcmy'[i%6]+',',c_d[i],c_t[i],'rgbcmy'[i%6]+'*')
xlabel(r'D, $\AA$', size='xx-large')
ylabel(r'$\phi$, $\degree$', position=(0.9,0), size='xx-large', va='center')
grid()

subplot(313)
title('Angle-Torsion', position=(0.05,0.95), ha='left', va='top', size='xx-large')
for i in range(Ndat):
    plot(angle[label==i],torsion[label==i],'rgbcmy'[i%6]+',',c_a[i],c_t[i],'rgbcmy'[i%6]+'*')
xlabel(r'$\alpha$, $\degree$', size='xx-large')
ylabel(r'$\phi$, $\degree$', position=(0.9,0), size='xx-large', va='center')
grid()

show()

