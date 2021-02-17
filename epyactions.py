from helper import range_check
import numpy as np

def strip_waters(args, model):
    model.WriteAtomList(args.outpath, model.atom_lister('notwater'), 'cell')

def extract_chains(args, model):
    if args.chids:
        model.writePDBchains(args.outpath, args.chids)
    else:
        print('This does not compute - extract chains but no chains listed?')

def set_b_per_chain(args, model):
    bavs = model.GetChainAverageBfactor()
    remodel = model.copy()
    for key,value in bavs.items():
        remodel.SetBfactorValues(float(value), what='chid', chid=key)
    remodel.writePDB(args.outpath)

def rename_chains(args, model):
    if args.chids:
        for chid1,chid2 in [(x[0],x[1]) for x in args.chids.split(',')]:
            model.rename_chain(chid1,chid2)
        model.writePDB(args.outpath)
    else:
        print('This does not compute - rename chains but no chain pairs listed?')

def rjust_resid(args, model):
    model.rjust_res_names()
    model.writePDB(args.outpath)

def extract_ranges(args, model):
    if args.ranges:
        model.extract_range(args.ranges).writePDB(args.outpath)
    else:
        print('This does not compute - extract but no ranges listed?')

def tinertia_ranges(args, model):
    if args.ranges:
        ramodel = model.extract_range(args.ranges)
        tin = ramodel.GetInertiaTensor('backbone')
    else:
        tin = model.GetInertiaTensor('backbone')
    print(tin.report())
    print(tin.frame())
    print(tin.anisou())

def tinertia_slider(args, model):
    pass

def center_orient(args, model):
    remodel = model.copy()
    remodel.shift(-remodel.GetCoM())
    tin_remodel = remodel.GetInertiaTensor()
    remodel.transform(tin_remodel.vr)
    remodel.writePDB(args.outpath)

def print_bvalue(args, model):
    resnames = model.GetResidueNames()
    resids, b0, b1, b2 = model.GetResidueBvectorByChain()
    for chid in resids:
        if not args.chids or chid in args.chids:
            for (i, resnum) in enumerate(resids[chid]):
                print('%3s %5s %6.2f %6.2f %6.2f' % (resnames[chid+resnum],\
                                                    chid+resnum.strip(),\
                                                    b0[chid][i],\
                                                    b1[chid][i],\
                                                    b2[chid][i]))

def print_chains(args, model):
    chanums = model.GetChains()
    bavs = model.GetChainAverageBfactor()
    print('Chain   Atoms   <B>')
    for key in chanums:
        print(key + str(chanums[key]).rjust(12) + " %6.2f" % bavs[key])

def print_phipsi(args, model):
    phi, psi = model.PhiPsiList()
    for resid in model.resid_lister('p_bb'):
        print("%6s %8.2f %8.2f" % (resid, phi.get(resid, float('nan')), psi .get(resid, float('nan'))))

def print_bcontrast(args, model):
    resids = [y for y in set([x.resid() for x in model.atoms]) if np.array([range_check(int(y[1:-1]),x) for x in args.ranges.get(y[0],[])]).any()]
    listik = model.atom_lister('resids', resids=resids)
    print('%d atoms selected for analysis\n<B>sel = %.2f' % (len(listik), model.GetAverageBfactor(listik=listik)))
    vicatoms = model.atom_lister('vicinity', model.atom_lister('notwater'), corelist=listik, rcutoff=args.rcutoff)
    print('%d non-water atoms found within %.2f Angstroms\n<B>vic = %.2f' % (len(vicatoms), args.rcutoff, model.GetAverageBfactor(listik=vicatoms)))

def print_hbonds(args, model):
    pass

def print_resgem(args, model):
    residue = model.get_residues('resid',resid=args.resid)[args.resid]
    b,a,t,m = residue.BondsAnglesTorsions(printout=True)
    print('Residue %s (%s)' % (args.resid, residue.get_res_name()))
    print('----------- Bonds ----------- ')
    for x in b:
        print(x)
    print('----------- Angles ----------- ')
    for x in a:
        print(x)
    print('----------- Torsions ----------- ')
    for x in t:
        print(x)
    print('----------- Impropers ----------- ')
    for x in m:
        print(x)

def print_range_bs(args, model):
    if args.ranges:
        print("%6.2f" % model.extract_range(args.ranges).GetAverageBfactor())
    elif args.chids:
        print("%6.2f" % model.extract_chains(args.chids).GetAverageBfactor())
    else:
        print("%6.2f" % model.GetAverageBfactor())

def print_rgyration(args, model):
    if args.ranges:
        print("%6.2f" % model.extract_range(args.ranges).Rgyration())
    elif args.chids:
        print("%6.2f" % model.extract_chains(args.chids).Rgyration())
    else:
        print("%6.2f" % model.Rgyration())
    
