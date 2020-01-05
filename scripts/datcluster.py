#!/usr/bin/env python3
headerhelp = \
'''
DATCONTUR reads a hb database that contains distance/angle/torsion columns
and presents the results as distributions.
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('--sqlpath', default='hbpisa.sqlite',
                    help='HB database file.  Defaults to hbpisa.sqlite.')
parser.add_argument('--symmetric',
                    action='store_true',
                    help='Bonds are chemically symmetric so both directions should be included.')
parser.add_argument('--d1', type=float,
                    help='Distance left boundary.')
parser.add_argument('--d2', type=float,
                    help='Distance right boundary.')
parser.add_argument('--a1', type=float,
                    help='Angle left boundary.')
parser.add_argument('--a2', type=float,
                    help='Angle right boundary.')
parser.add_argument('--t1', type=float,
                    help='Torsion left boundary.')
parser.add_argument('--t2', type=float,
                    help='Torsion right boundary.')
parser.add_argument('--torshift', type=float,
                    help='Shift the torsion angle domain')
parser.add_argument('--kmeans',
                    type=int, default=2,
                    help='Number of clusters to use to perform K-means clustering analysis.')
parser.add_argument('--weighted', action='store_true',
                    help='Weight distribution by redundancies (input file should be processed by SEQFILTER')
parser.add_argument('--title', default='DAT clusters',
                    help='Title of the figure.')
args = parser.parse_args()

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from aconts import hbond_pdbase

from scipy import array, ones, bincount, isfinite, ravel
from matplotlib.pyplot import figure, title, subplot, show, grid, xlabel, ylabel, plot
from scipy.stats import iqr

'''
def f2gauss(xy, p):
    A, xo, yo, sx, sy = p
    return (A*exp(-(xy[0]-xo)**2/sx-(xy[1]-yo)**2/sy))
def lsq2gauss(p, x, y, z, n):
    return sum(((f2gauss((x,y),p)-z)**n).ravel())
def lsqwgauss(p, x, y, z):
    dfv = f2gauss((x,y),p) - z
    w = z
    return sum((w*(dfv**2)).ravel())
def da_pvalue(kp, d, a):
    xo, sx, yo, sy = kp
    return exp(-sqrt((d-xo)**2/sx+(a-yo)**2/sy))
'''

hpdb = hbond_pdbase(args.sqlpath)
if args.symmetric:
    distance,angle,torsion=array([(x.d,x.angle1,x.tor1) for x in [x[1] for x in hpdb.get_hbonds()]] +
                                 [(x.d,x.angle2,x.tor2) for x in [x[1] for x in hpdb.get_hbonds()]]).T.astype(float)
else:
    distance,angle,torsion=array([(x.d,x.angle1,x.tor1) for x in [x[1] for x in hpdb.get_hbonds()]]).T.astype(float)
redcount = ones(len(distance))

ind = isfinite(distance+angle+torsion)
distance, angle, torsion, redcount = distance[ind], angle[ind], torsion[ind], redcount[ind]

if args.weighted:
    pass # This needs completion

if args.torshift!=None:
    torsion[torsion<args.torshift] += 360.0

if args.d1==None: 
    args.d1=distance.min()
else:
    ind = distance>args.d1
    distance, angle, torsion, redcount = distance[ind], angle[ind], torsion[ind], redcount[ind]
if args.d2==None: 
    args.d2=distance.max()
else:
    ind = distance<args.d2
    distance, angle, torsion, redcount = distance[ind], angle[ind], torsion[ind], redcount[ind]
if args.a1==None: 
    args.a1=angle.min()
else:
    ind = angle>args.a1
    distance, angle, torsion, redcount = distance[ind], angle[ind], torsion[ind], redcount[ind]
if args.a2==None: 
    args.a2=angle.max()
else:
    ind = angle<args.a2
    distance, angle, torsion, redcount = distance[ind], angle[ind], torsion[ind], redcount[ind]
if args.t1==None: 
    args.t1=torsion.min()
else:
    ind = torsion>args.t1
    distance, angle, torsion, redcount = distance[ind], angle[ind], torsion[ind], redcount[ind]
if args.t2==None: 
    args.t2=torsion.max()
else:
    ind = torsion<args.t2
    distance, angle, torsion, redcount = distance[ind], angle[ind], torsion[ind], redcount[ind]


figure(args.title)

from scipy.cluster.vq import kmeans2, whiten
z = whiten(array([distance,angle,torsion]).T)
centroid, label = kmeans2(z, args.kmeans, minit='points')
freqs = 100*bincount(label)/len(label)
c_d = [distance[label==i].mean() for i in range(args.kmeans)]
c_a = [angle[label==i].mean() for i in range(args.kmeans)]
c_t = [torsion[label==i].mean() for i in range(args.kmeans)]
s_d = [iqr(distance[label==i],scale='normal') for i in range(args.kmeans)]
s_a = [iqr(angle[label==i],scale='normal') for i in range(args.kmeans)]
s_t = [iqr(torsion[label==i],scale='normal') for i in range(args.kmeans)]
for (i, c) in enumerate(array([[c_d,s_d],[c_a,s_a],[c_t,s_t]]).T):
    print("Cluster #%-2d: " % (i+1) + "%8.3f %8.3f "*len(c.T) % tuple(ravel(c.T)) + "%5.1f%%" % (freqs[i]))
    datvalues = ravel(array([c[0],2*c[1]**2]).T)
    print("               (" + ','.join("%.3f"*len(datvalues)) % (datvalues) + ")")

subplot(311)
title('Distance-Angle', position=(0.05,0.95), ha='left', va='top', size='xx-large')
for i in range(args.kmeans):
    plot(distance[label==i],angle[label==i],'rgbcmy'[i%6]+',',c_d[i],c_a[i],'rgbcmy'[i%6]+'*')
xlabel(r'D, $\AA$', size='xx-large')
ylabel(r'$\alpha$, $\degree$', position=(0.9,0), size='xx-large', va='center')
grid()

subplot(312)
title('Distance-Torsion', position=(0.05,0.95), ha='left', va='top', size='xx-large')
for i in range(args.kmeans):
    plot(distance[label==i],torsion[label==i],'rgbcmy'[i%6]+',',c_d[i],c_t[i],'rgbcmy'[i%6]+'*')
xlabel(r'D, $\AA$', size='xx-large')
ylabel(r'$\phi$, $\degree$', position=(0.9,0), size='xx-large', va='center')
grid()

subplot(313)
title('Angle-Torsion', position=(0.05,0.95), ha='left', va='top', size='xx-large')
for i in range(args.kmeans):
    plot(angle[label==i],torsion[label==i],'rgbcmy'[i%6]+',',c_a[i],c_t[i],'rgbcmy'[i%6]+'*')
xlabel(r'$\alpha$, $\degree$', size='xx-large')
ylabel(r'$\phi$, $\degree$', position=(0.9,0), size='xx-large', va='center')
grid()

show()
