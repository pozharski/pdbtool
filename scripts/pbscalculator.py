#! /usr/bin/env python3

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

version = '0.2'

defopts = {
'xyzin1'    :   '',
'xyzin2'    :   '',
'xyzout1'   :   '__NONE__',
'xyzout2'   :   '__NONE__',
'mode'      :   'ALL',
'fitm'      :   '',
'chai'      :   '_',
'nbin'      :   '50',
'nmax'      :   '3',
'wcut'      :   '95',
'outp'      :   'ALL',
'refc'      :   '0',
'bfac'      :   'PROBABILITY',
'bave'      :   'NO',
'chrn'      :   'NO'}

import lsq
from math import sqrt

print ' ------------------- '
print '|PBSCALC version 0.2|'
print ' ------------------- '


for (key, value) in zip(*[sys.argv[1::2],sys.argv[2::2]]):
    defopts[key.lower()] = value

comments = ''
ignored = ''
#for line in sys.stdin:
while True:
    try:
        line = raw_input()
    except EOFError:
        break
    if line.strip():
        if line[:4].strip().upper() == 'END':
            break
        try:
            (key, value) = line.strip().split()[:2]
            kword = key[:4].lower()
            if line[0] == '#':
                comments += line[1:]
            elif kword in defopts:
                defopts[kword] = value
            else:
                ignored += '# - ' + line
        except:
            print 'Failed to interpret this line in the input script:\n"' + line.strip() + '"'
if comments:
    print 'Comments from the input file:\n\n'+comments+'\n'
if ignored:
    print "The following commands were ignored (no such keyword):\n\n"+ignored

if not defopts['fitm']:
    defopts['fitm'] = defopts['mode']
for key in ['mode','fitm','outp','bfac','bave','chrn']:
    defopts[key] = defopts[key].upper()
exitFlag = False
for (key, value) in defopts.iteritems():
    if not value:
        print 'Missing mandatory parameter '+key.upper()
        exitFlag = True
if exitFlag:
    sys.exit()

(ch1,ch2)=list(defopts['chai'].split('_'))
chm = (list(ch1),list(ch2))
algn_dict=None
if defopts['mode'] == 'SEQUENCE_ALIGNMENT' or defopts['fitm'] == 'SEQUENCE_ALIGNMENT':
    algn_dict=lsq.SSMSmatch(defopts['xyzin1'], defopts['xyzin2'])
resCheck = (defopts['chrn'] == 'YES')
al = lsq.Aligner(defopts['xyzin1'], defopts['xyzin2'], fListBuild=False)
al.BuildList(mode=eval('lsq.MATCH_'+defopts['fitm']), chains2match=chm, algn_dict=algn_dict, resCheck=resCheck)
print 'Kabsch algorithm applied to align molecules'
print str(al.GetMatchLength()) + ' atoms were selected for alignment using "'+defopts['fitm']+'" mode'
print 'Residue name checking is ' + ('ON' if resCheck else 'OFF')
al.FitKabsch()
for i in range(int(defopts['refc'])):
    al.FitPBS()
al.BuildList(mode=eval('lsq.MATCH_'+defopts['mode']), chains2match=chm, algn_dict=algn_dict, resCheck=resCheck)
print str(al.GetMatchLength()) + ' atoms matched using "'+defopts['mode']+'" mode'
print 'Residue name checking is ' + ('ON' if resCheck else 'OFF')
if al.GetMatchLength():
    outp = defopts['outp'].lower()
    Nx = int(defopts['nbin'])
    if outp == 'all' or outp == 'rmsd' or outp == 'rmsdpbs':
        print "\n=================================="
        print "RMSD between two models = %6.3f" % al.rmsd()
    if outp == 'all' or outp == 'pbs' or outp == 'rmsdpbs':
        print "PBS  between two models = %6.3f" % al.pbs()
        print "==================================\n"
    if outp == 'all' or outp == 'quantile':
        print "\n=================================="
        print 'PERCENTILE    P.B.S., A'
        for z in  zip(*al.QuantileSpread()):
            print '%6d     %8.3f' % z
        print "==================================\n"
    if outp == 'all' or outp == 'distr':
        (s,w,tfft) = al.FitHistogram(Nx=Nx)
        print "==The distribution of distances==="
        print " D, A      p.d.f.     fit"
        for z in  zip(*tfft):
            print '%6.3f   %8.4f  %8.4f' % z
        print "==================================\n"
        print 'Fit of the distance distribution to Maxwell-Boltzmann\nresults in estimated spread of\n     D=%6.3f A\nwhich accounts for \n      w=%4.1f %%\nof atoms in the match.\n' % (sqrt(1.5)*s,100.0*w)
        print "==================================\n"
    (s,w,(t,f,ft)) = al.FitHistogramUpToN(N=int(defopts['nmax']),Nx=Nx,wcut=.01*float(defopts['wcut']))
    if outp == 'all' or outp == 'mult':
        print "Now fitting the distance distribution to up to "+defopts['nmax']+" groups of atoms"
        if len(s) > 1:
            print "Group #      D,A         w"
            for (i,(ss,ww)) in enumerate(zip(*[s,w])):
                print "  %3d      %6.3f      %5.3f" % (i+1,ss,ww)
            linout = " D, A      p.d.f.     fit"
            fft = zip(*ft)
            for i in range(len(ft)):
                 linout += "  %8d" % (i+1)
            print "\n=================================="
            print linout
            for (i,tt) in enumerate(t):
                linout = '%6.3f   %8.4f  %8.4f' % (tt, f[i], sum(fft[i]))
                for ff in fft[i]:
                    linout += '  %8.4f' % ff
                print linout
        else:
            print "Unimodal distribution is sufficient to account for over "+defopts['wcut']+"% of atoms."
    if defopts['bfac'] == 'PROBABILITY':
        al.SetBprob(s,w)
    elif defopts['bfac'] == 'ZSCORE':
        al.SetBZscore()
    if defopts['bave'] == 'YES':
        al.SetB_PeResidue()
    if defopts['xyzout1'] == '__NONE__':
        defopts['xyzout1'] = ''
    if defopts['xyzout2'] == '__NONE__':
        defopts['xyzout2'] = ''
    al.WriteModels(defopts['xyzout1'],defopts['xyzout2'])
    if outp == 'all' or outp == 'residues':
        al.SetBZscore()
        ((resn1,(bf1,bb1,bs1)),(resn2,(bf2,bb2,bs2))) = al.GetBLists()
        print "==The Z-scores of individual residues==="
        print "============= Molecule #1 =============="
        print "("+defopts['xyzin1']+")"
        print 'Residue  Total Backbone Sidechain' 
        for (i,name1) in enumerate(resn1):
            print name1+'  %6.3f %6.3f   %6.3f' % (bf1[i],bb1[i],bs1[i])
        print "============= Molecule #2 =============="
        print "("+defopts['xyzin2']+")"
        print 'Residue  Total Backbone Sidechain' 
        for (i,name2) in enumerate(resn2):
            print name2+'  %6.3f %6.3f   %6.3f' % (bf2[i],bb2[i],bs2[i])
        print "========================================"
        
else:
    print "Nothing to match!  Selection is too restrictive."
