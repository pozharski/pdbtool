#! /usr/bin/env python3
headerhelp = \
'''
BBCOMP runs the comparison of backbone angles (phi/psi/omega) between
    two conformations.
'''


from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('model1',
                    help='First model PDB file.')
parser.add_argument('model2',
                    help='Second model PDB file.')
parser.add_argument('-l', '--list', action='store_true',
                    help="List the residues over the cutoff values")
parser.add_argument('--delta-phi', type=float, default=10.0,
					help='Phi angle change cutoff.')
parser.add_argument('--delta-psi', type=float, default=10.0,
					help='Psi angle change cutoff.')
parser.add_argument('--delta-omega', type=float, default=10.0,
					help='Omega angle change cutoff.')
parser.add_argument('--bbout',
                    help='Path to the output PDB file showing backbone differences as B-factors.')
args = parser.parse_args()

from pdbtool import ReadPDBfile as read_pdb_file
from scipy import degrees, radians, sin, arcsin, sqrt
from scipy.stats import scoreatpercentile

model1 = read_pdb_file(args.model1)
if model1 is None:
    sys.exit('Failed to read the model coordinates from '+args.model1)
model2 = read_pdb_file(args.model2)
if model2 is None:
    sys.exit('Failed to read the model coordinates from '+args.model2)

phi1,psi1,omega1 = model1.BackboneTorsions()
phi2,psi2,omega2 = model2.BackboneTorsions()

phikeys = set(phi1).union(phi2)
psikeys = set(psi1).union(psi2)
omegakeys = set(omega1).union(omega2)
reskeys = phikeys.union(psikeys).union(omegakeys)

phis = dict([(key,[phi1.get(key),phi2.get(key)]) for key in phikeys])
psis = dict([(key,[psi1.get(key),psi2.get(key)]) for key in psikeys])
omegas = dict([(key,[omega1.get(key),omega2.get(key)]) for key in omegakeys])

phiv = dict([(k,degrees(arcsin(sin(radians(v[0]-v[1]))))) for k,v in phis.items() if not None in v])
psiv = dict([(k,degrees(arcsin(sin(radians(v[0]-v[1]))))) for k,v in psis.items() if not None in v])
omegav = dict([(k,degrees(arcsin(sin(radians(v[0]-v[1]))))) for k,v in omegas.items() if not None in v])

phiscores = [abs(x) for x in phiv.values() if x is not None]
psiscores = [abs(x) for x in psiv.values() if x is not None]
omegascores = [abs(x) for x in omegav.values() if x is not None]
print("Percentiles    50%    90%    95%    99%   N")
print("phi         %6.1f %6.1f %6.1f %6.1f" % tuple(scoreatpercentile(phiscores,[50,90,95,99])),str(len(phiscores)))
print("psi         %6.1f %6.1f %6.1f %6.1f" % tuple(scoreatpercentile(psiscores,[50,90,95,99])),str(len(psiscores)))
print("omega       %6.1f %6.1f %6.1f %6.1f" % tuple(scoreatpercentile(omegascores,[50,90,95,99])),str(len(omegascores)))

if args.list:
    for key in sorted(reskeys):
        if key in phis:
            if None in phis[key]:
                print("%7s %s %s --- phi missing" % tuple([key]+['None  ' if v is None else '%6.1f' % (v) for v in phis[key]]))
            elif abs(phiv[key])>args.delta_phi:
                print("%7s %6.1f %6.1f %6.1f --- phi" % tuple([key,phiv[key]]+phis[key]))
        if key in psis:
            if None in psis[key]:
                print("%7s %s %s --- psi missing" % tuple([key]+['None  ' if v is None else '%6.1f' % (v) for v in psis[key]]))
            elif abs(psiv[key])>args.delta_psi:
                print("%7s %6.1f %6.1f %6.1f --- psi" % tuple([key,psiv[key]]+psis[key]))
        if key in omegas:
            if None in omegas[key]:
                print("%7s %s %s --- omega missing" % tuple([key]+['None  ' if v is None else '%6.1f' % (v) for v in omegas[key]]))
            elif abs(omegav[key])>args.delta_omega:
                print("%7s %6.1f %6.1f %6.1f --- omega" % tuple([key,omegav[key]]+omegas[key]))

if args.bbout is not None:
    bvalues = dict([(key, sqrt(sum([x**2 for x in [phiv.get(key),psiv.get(key),omegav.get(key)] if x]))) for key in reskeys])
    for atom in model1.atom_getter():
        atom.SetB(bvalues.get(atom.resid(),0.0))
    model1.writePDB(args.bbout)
    
        
