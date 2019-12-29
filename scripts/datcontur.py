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
parser.add_argument('--dstep', type=float, default=0.02,
                    help='Distance step size.')
parser.add_argument('--astep', type=float, default=1.0,
                    help='Angle step size.')
parser.add_argument('--tstep', type=float, default=10.0,
                    help='Torsion step size.')
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
parser.add_argument('--normden',
                    action='store_true',
                    help='Use density option when building histograms.')
parser.add_argument('--dakernel',
                    action='store_true',
                    help='Calculate parameters of the DA kernel (slow).')
parser.add_argument('--weighted', action='store_true',
                    help='Weight distribution by redundancies (input file should be processed by SEQFILTER')
parser.add_argument('--title', default='DAT distributions',
                    help='Title of the figure.')
args = parser.parse_args()

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from aconts import hbond_pdbase

from scipy import loadtxt, array, arange, histogram, unravel_index, isfinite, meshgrid, argmax, exp, sqrt, pi, ones
from scipy.optimize import fmin
from matplotlib.pyplot import subplot, contourf, contour, show, grid, figure, title, hist, xlabel, ylabel, gca

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

hpdb = hbond_pdbase(args.sqlpath)
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

td = arange(args.d1,args.d2,args.dstep)
ta = arange(args.a1,args.a2,args.astep)
tt=arange(args.t1,args.t2,args.tstep)

figure(args.title)

subplot(321)
title('Distance distribution', position=(0.05,0.95), ha='left', va='top', size='xx-large')
hval=hist(distance, td, weights=redcount)
xlabel(r'D, $\AA$', position=(0.9,0), size='xx-large', va='center')
gca().yaxis.set_ticklabels([])
grid()

subplot(322)
title('Distance vs angle', position=(0.05,0.95), color='white', ha='left', va='top', size='xx-large')
zda=array([histogram(angle[abs(distance-x)<=0.5*args.dstep], bins=ta, weights=redcount[abs(distance-x)<=0.5*args.dstep], density=args.normden)[0] for x in td])
contourf(0.5*(ta[1:]+ta[:-1]),td,array(zda),50)
xlabel(r'$\alpha$, $\degree$', position=(0.9,0), size='xx-large', va='center')
ylabel(r'D, $\AA$', size='xx-large')
grid()

if args.dakernel:
    xv,yv=meshgrid(0.5*(ta[1:]+ta[:-1]), td)
    ti, tj = unravel_index(argmax(zda),zda.shape)
    xo = xv[ti,tj]
    yo = yv[ti,tj]
    sx = 2*sum(sum(zda>zda.max()*exp(-1)))/pi*args.astep**2
    sy = 2*sum(sum(zda>zda.max()*exp(-1)))/pi*args.dstep**2
    Anorm=zda.sum()*sqrt(pi*(sx+sy))*args.dstep*args.astep
    p = [Anorm, xo, yo, sx, sy]
    print("Round 1:")
    print("Initial parameters: D/sigmaD = %.2f / %.2f, A/sigmaA = %.1f / %.1f" % (p[2], sqrt(p[4]/2), p[1], sqrt(p[3]/2)))
    pp = fmin(lsqwgauss, p, (xv, yv, zda))
    print("Refined parameters: D/sigmaD = %.2f / %.2f, A/sigmaA = %.1f / %.1f" % (pp[2], sqrt(pp[4]/2), pp[1], sqrt(pp[3]/2)))
    print("Round 2:")
    zv = f2gauss((xv,yv),pp)
    vv = zv.max()*array([exp(-8), exp(-4.5), exp(-2), exp(-0.5)])
    contour(xv,yv,zv,vv,colors='white')
    ppp = fmin(lsqwgauss, pp, (xv, yv, zda))
    print("Refined parameters: D/sigmaD = %.2f / %.2f, A/sigmaA = %.1f / %.1f" % (ppp[2], sqrt(ppp[4]/2), ppp[1], sqrt(ppp[3]/2)))
    zv = f2gauss((xv,yv),ppp)
    vv = zv.max()*array([exp(-8), exp(-4.5), exp(-2), exp(-0.5)])
    contour(xv,yv,zv,vv,colors='yellow')
    print("Parameters of the kernel estimator (X - distance, Y - angle):")
    print("A = %.4f\nX0 = %.4f\nSx = %.4f\nY0 = %.4f\nSy = %.4f" % (1./pi/sqrt(ppp[3]*ppp[4]), ppp[2], ppp[4], ppp[1], ppp[3]))
    kps = (ppp[2], ppp[4], ppp[1], ppp[3])
    print("ACONTS line: (%.4f, %.4f, %.4f, %.4f)" % kps)
    pvs = da_pvalue(kps, distance, angle)
    Neff = sum(redcount[pvs>0.05])
    print("Lower estimate of the effective number of measurements: %d" % int(Neff))
    print("Average distance: %.4f +- %.4f" % (ppp[2], sqrt(ppp[4]/2/Neff)))
    print("Average angle: %.3f +- %.3f" % (ppp[1], sqrt(ppp[3]/2/Neff)))

subplot(323)
title('Angle distribution', position=(0.05,0.95), ha='left', va='top', size='xx-large')
hval=hist(angle, ta, weights=redcount)
gca().yaxis.set_ticklabels([])
xlabel(r'$\alpha$, $\degree$', position=(0.9,0), size='xx-large', va='center')
grid()

subplot(324)
title('Distance vs torsion', position=(0.05,0.95), color='white', ha='left', va='top', size='xx-large')
zdt=array([histogram(torsion[abs(distance-x)<=0.5*args.dstep], bins=tt, weights=redcount[abs(distance-x)<=0.5*args.dstep], density=args.normden)[0] for x in td])
contourf(0.5*(tt[1:]+tt[:-1]),td,array(zdt),50)
xlabel(r'$\phi$, $\degree$', position=(0.9,0), size='xx-large', va='center')
ylabel(r'D, $\AA$', size='xx-large')
grid()

subplot(325)
title('Torsion distribution', position=(0.05,0.95), ha='left', va='top', size='xx-large')
hval=hist(torsion, tt, weights=redcount)
gca().yaxis.set_ticklabels([])
xlabel(r'$\phi$, $\degree$', position=(0.9,0), size='xx-large', va='center')
grid()

subplot(326)
title('Angle vs torsion', position=(0.05,0.95), color='white', ha='left', va='top', size='xx-large')
zat=array([histogram(torsion[abs(angle-x)<=0.5*args.astep], bins=tt, weights=redcount[abs(angle-x)<=0.5*args.astep], density=args.normden)[0] for x in ta])
contourf(0.5*(tt[1:]+tt[:-1]),ta,array(zat),50)
xlabel(r'$\phi$, $\degree$', position=(0.9,0), size='xx-large', va='center')
ylabel(r'$\alpha$, $\degree$', size='xx-large')
grid()

show()
