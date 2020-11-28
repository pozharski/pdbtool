from pdbtool import pdbmolecule, ReadPDBfile
from scipy import exp, sqrt, log, pi, fmod, array
from scipy.special import erfc
from scipy.optimize import fsolve

TWOBYSQRTPI = 2 / sqrt(pi)
FOURBYTHREESQRTPI = 4 / 3 / sqrt(pi)

DAKERNEL_PARAMS = {
    'TyrOHtoAnionic'            : (2.6164, 0.0257, 114.8657,  76.3537),
    'TyrOHtoAspOD'              : (2.6149, 0.0233, 115.0366,  75.3831),
    'TyrOHtoGluOE'              : (2.6187, 0.0294, 114.6956,  78.0677),
    'TyrOHtoCarboxamideO'       : (2.6584, 0.0251, 114.6189,  71.8621),
    'TyrOHtoAsnOD1'             : (2.6544, 0.0198, 114.4375,  67.8227),
    'TyrOHtoGlnOE1'             : (2.6681, 0.0383, 114.7680,  91.6975),
    'TyrOHtoTyrOH'              : (2.7120, 0.0587, 115.5739, 156.4656),
    'TyrOHtoThrOG1'             : (2.7089, 0.0483, 115.8939, 145.4301),
    'TyrOHtoHis'                : (2.7048, 0.0242, 116.9498,  95.7449),
    'TyrOHtoSerOG'              : (2.7157, 0.0398, 115.3121, 131.9901),
    'TyrOHtoBackbone'           : (2.6650, 0.0224, 113.3358,  75.7802),
    'TrpNE1toAnionic'           : (2.8718, 0.0362, 127.1243, 281.2524),
    'TrpNE1toCarboxamideO'      : (2.8638, 0.0368, 129.2766, 416.3360),
    'TrpNE1toAsnOD1'            : (2.8777, 0.0573, 129.3067, 933.8468),
    'TrpNE1toGlnOE1'            : (2.8583, 0.0388, 129.1742, 340.5477),
    'ThrOG1toAnionic'           : (2.6803, 0.0404, 110.6122, 168.8166),
    'TrpNE1toAspOD'             : (2.8746, 0.0401, 126.3796, 342.5100),
    'TrpNE1toGluOE'             : (2.8713, 0.0460, 127.2237, 300.2707),
    'TrpNE1toTyrOH'             : (2.9806, 0.0240, 117.9992, 384.3989),
    'TrpNE1toSerOG'             : (2.9090, 0.0327, 125.5139, 218.3441),
    'TrpNE1toThrOG1'            : (2.9552, 0.0962, 125.2985, 686.3387),
    'TrpNE1toHydroxyl'          : (2.9360, 0.0460, 124.1150, 395.5564),
    'TrpNE1toHis'               : (2.9858, 0.1686, 124.7804, 469.5073), 
    'TrpNE1toBackbone'          : (2.8971, 0.0377, 129.0206, 643.4862),
    'SerOGtoAnionic'            : (2.6576, 0.0350, 107.3850, 232.5428),
    'SerOGtoTyrOH'              : (2.7176, 0.0647, 113.6966, 584.9499),
    'AnionicToAnionic'          : (2.5441, 0.0374, 117.1554,  32.7188),
    'SerOGtoCarboxamideO'       : (2.6960, 0.0451, 108.5755, 321.7805),
    'ThrOG1toCarboxamideO'      : (2.7380, 0.0492, 112.2677, 192.7388),
    'CarboxamideNtoCarboxamideO': (2.9280, 0.0557, 116.3558, 253.4728),
    'SerOGtoSerOG'              : (2.7441, 0.0440, 109.3314, 239.4919),
    'CarboxamideNToAnionic'     : (2.9438, 0.1145, 120.0690, 298.5215),
    'ThrOG1toThrOG1'            : (2.7691, 0.0514, 117.6728, 244.1296),
    'SerOGtoThrOG1'             : (2.7527, 0.0471, 113.7098, 578.1031), 
    'ThrOG1toSerOG'             : (2.7594, 0.0569, 118.8537, 496.2742),
    'CarboxamideNToTyrOH'       : (3.0485, 0.1052, 116.2582, 344.1567),
    'CarboxamideNToThrOG1'      : (2.9576, 0.0623, 118.4897, 303.5630),
    'CarboxamideNToSerOG'       : (2.9532, 0.1051, 118.2686, 357.3830),
    'CarboxamideNToHydroxyl'    : (2.9770, 0.0726, 118.3963, 304.6859),
    'LysNZtoAnionic'            : (2.8104, 0.0799, 105.4847, 433.1885),
    'LysNZtoCarboxamideO'       : (2.8323, 0.0878, 103.5184, 542.9631),
    'LysNZtoHydroxyl'           : (2.9322, 0.1182, 101.3103, 1106.9116),
    'LysNZtoSerOG'              : (2.9112, 0.1343, 103.6097, 968.3556),
    'LysNZtoThrOG1'             : (2.8751, 0.0871, 103.7193, 943.1082),
    'ArgNEtoAnionic'            : (2.8553, 0.0454, 116.6745, 103.5459),
    'LystoBackbone'             : (2.8393, 0.0791, 103.1729, 665.7838),
    'AsnNtoBackbone'            : (2.9194, 0.0497, 116.7913, 257.5720),
    'GlnNtoBackbone'            : (2.9480, 0.0640, 118.1744, 297.2644),
    'CarboxamideNtoBackbone'    : (2.9311, 0.0559, 117.3667, 277.5849),
    'ThrOG1toBackbone'          : (2.7694, 0.0745, 105.6189, 239.0208),
    'SerOGtoBackbone'           : (2.7438, 0.0769, 105.6017, 464.1153),
    'ThrOG1toTyrOH'             : (2.7195, 0.0430, 119.4327, 505.4977),
    'BackboneNtoBackboneO'      : (2.9359, 0.0767, 119.2287, 102.2321),
    'BackboneNtoAspOD'          : (2.9290, 0.0702, 124.9471, 274.0948),
    'BackboneNtoHis'            : (3.0180, 0.0393, 125.7693, 143.6583),
    'BackboneNtoHisNE2'         : (2.9845, 0.0477, 118.3224,  78.9542),
    'BackboneNtoHisND1'         : (3.0180, 0.0350, 126.8302, 101.5115),
    'BackboneNtoCarboxamideO'   : (2.9109, 0.0556, 120.0171, 321.0081),
    'BackboneNtoAsnOD1'         : (2.9266, 0.0591, 123.6668, 345.3618),
    'BackboneNtoGlnOE1'         : (2.8906, 0.0419, 115.5596, 120.6308),
    'BackboneNtoSerOG'          : (3.0921, 0.0684, 127.4116, 59.7068),
    'BackboneNtoThrOG1'         : (3.1113, 0.0694, 127.3112, 73.5292),
    'BackboneNtoTyrOH'          : (2.9480, 0.0559, 117.7464, 156.2317),
    'HistoBackbone'             : (2.8433, 0.0799, 120.3890, 648.0583),
    'HistoHis'                  : (2.8847, 0.0810, 123.7450, 320.4533),
}

DATKERNEL_PARAMS = {
    'AnionicToAnionic'          : [(2.559,0.041,118.383,309.864,179.909,1746.392),
                                   (2.571,0.043,123.385,302.686,361.596,2410.010)],
    'TyrOHtoAnionic'            : [(2.613,0.030,116.558,103.868,361.762,781.720),
                                   (2.620,0.031,116.565,97.802,180.024,809.980)],
    'TyrOHtoAspOD'              : [(2.616,0.029,116.343,87.408,180.229,751.473),
                                   (2.614,0.028,116.103,100.761,361.988,758.558)],
    'TyrOHtoGlnOE1'             : [(2.612,0.031,116.976,109.722,180.032,871.674),
                                   (2.620,0.034,117.089,114.210,361.910,832.923)],
    'TyrOHtoCarboxamideO'       : [(2.683,0.021,114.899,64.139,181.838,514.651),
                                   (2.679,0.023,115.001,64.627,362.605,572.587)],
    'TyrOHtoTyrOH'              : [(2.756,0.076,116.737,202.242,178.847,2523.879),
                                   (2.751,0.077,116.970,221.608,357.934,3877.188)],
}

def pcutoff(pv, hbtp):
    xo, sx = DAKERNEL_PARAMS[hbtp][:2]
    return xo - sqrt(sx)*log(pv)
def da_pvalue(hbtp, d, a):
    xo, sx, yo, sy = DAKERNEL_PARAMS[hbtp]
    return exp(-sqrt((d-xo)**2/sx+(a-yo)**2/sy))
def erf3(x):
    return erfc(x) + TWOBYSQRTPI*x*exp(-x**2)
def dat_pvalues(hbtp, d, a, t):
    return [erf3(sqrt((d-xo)**2/sx+(a-yo)**2/sy+(fmod(abs(t-zo-180),360)-180)**2/sz)) for (xo,sx,yo,sy,zo,sz) in DATKERNEL_PARAMS[hbtp]]
def erfinv3(p):
    def f3(x,p):
        return erf3(x)-p
    return fsolve(f3,1,p)
def sigmap(s,p):
    x = erfinv3(p)
    return s / (1-p-FOURBYTHREESQRTPI*x**3*exp(-x**2))


class AtomContact:
    '''
    Atom contacts extraction, storage, manipulation and analysis.
    '''
    def __init__(self, model, *args, **kwds):
        self._model = model
        self._aconts = {}
    def store_rats(self, rats, rcutoff=4.0, sameres=False, userkey=False):
        '''
        Extract and store contacts for certain combinations of residue
        and atom names.  Argument rats is a list of such combinations to
        include.  For example, ['LYSNZ', 'TYROH'] will include NZ atom
        of every lysine and OH atom for every tyrosine when they are as
        close as rcutoff.
        '''
        if userkey:
            key = userkey
        else:
            'rats_'+'_'.join(sorted(rats))
        if key not in self._aconts:
            listik = self._model.atom_lister('rat', rats=rats)
            submol = pdbmolecule(atoms=self._model.atom_getter('all', listik))
            self._aconts[key] = [(listik[x[0]], listik[x[1]], x[2]) for x in submol.AllContacts(rcutoff)]
            if not sameres:
                self._aconts[key] = [x for x in self._aconts[key] if not self._model.same_residue(x[0], x[1])]
    def store_lists(self, listik1, listik2, key, rcutoff=4.0, sameres=False):
        if key not in self._aconts:
            self._aconts[key] = list(self._model.IndexContacts(listik1, listik2, rmax=rcutoff))
            if not sameres:
                self._aconts[key] = [x for x in self._aconts[key] if not self._model.same_residue(x[0], x[1])]
    def report(self):
        for key, value in self._aconts.items():
            print("%s %5d" % (key, len(value)))
            for i,j,r in value[:3]:
                print("%-16s %-16s %.2f " % (self._model.GetAtomTitle(i), self._model.GetAtomTitle(j), r))
    def key_report(self, key):
        retval = ''
        value = self._aconts.get(key)
        if value is not None:
            for i,j,r in value[:3]:
                retval += "%-16s %-16s %.2f \n" % (self._model.GetAtomTitle(i), self._model.GetAtomTitle(j), r)
        return retval

class ACreader(object):
    ''' Reads the report lines produced by contact classes '''
    def __init__(self, line=False, *args, **kwds):
        self.readline(line)
    def readline(self, line=False):
        if type(line) is str:
            chunks = line.split()
        else:
            chunks = line
        if line:
            self.resn1, self.resid1, self.atom1 = chunks[:3]
            self.resn2, self.resid2, self.atom2 = chunks[3:6]
            self.d = float(chunks[6])
        else:
            chunks = None
            self.resn1, self.resid1, self.atom1, self.resn2, self.resid2, self.atom2 = ['']*6
            self.d = float('nan')
        self._extra_reads(chunks)
    def _extra_reads(self, chunks):
        pass
    def swap(self):
        self.resn2, self.resid2, self.atom2, self.resn1, self.resid1, self.atom1 = self.resn1, self.resid1, self.atom1, self.resn2, self.resid2, self.atom2
        self._extra_swap()
    def _extra_swap(self):
        pass
    def report(self):
        return "%3s %5s %-5s %3s %5s %-5s %7.2f " % (self.resn1, self.resid1,self.atom1, self.resn2, self.resid2, self.atom2, self.d) + self._extra_report()
    def _extra_report(self):
        return ''
    def items(self):
        return [self.resn1, self.resid1,self.atom1, self.resn2, self.resid2, self.atom2, self.d] + self._extra_items()
    def _extra_items(self):
        return []

from pdbminer import pdbase
class atom_contact_pdbase(pdbase):
    tables = [('atom_contacts', ACreader(), '')]
    def process_code_rats(self, code, fpath, rats, rcutoff=4.0):
        ac = AtomContact(ReadPDBfile(fpath))
        ac.store_rats(rats, rcutoff=rcutoff, userkey='whatever')
        self.insert_new_code(code)
        things = [ACreader(x) for x in ac.key_report('whatever').strip('\n').split('\n') if len(x)]
        if len(things):
            for thing in things:
                self.insert_new_item('atom_contacts',code,thing)
            print("Processed "+code+"... (got %d)" % (len(things)))
        else:
            print("Processed "+code+"...")
        self.code_lock(code)

class HBreader(ACreader):
    def _extra_reads(self, chunks):
        if chunks:
            self.angle1, self.tor1 = [float(x) for x in chunks[7:9]]
            self.angle2, self.tor2 = [float(x) for x in chunks[9:11]]
            self.atomi1, self.atomi2 = [int(x) for x in chunks[11:13]]
        else:
            self.angle1, self.tor1, self.angle2, self.tor2 = [float('nan')]*4
            self.atomi1, self.atomi2 = [-1,-1]
    def _extra_swap(self):
        self.atomi1, self.angle1, self.tor1, self.atomi2, self.angle2, self.tor2 = self.atomi2, self.angle2, self.tor2, self.atomi1, self.angle1, self.tor1
    def _extra_report(self):
        return "%10.2f %10.2f %10.2f %10.2f" % (self.angle1, self.tor1, self.angle2, self.tor2)
    def _extra_items(self):
        return [self.angle1, self.tor1, self.angle2, self.tor2]
    def bonafide(self,d=None,a=None,t=None,sd=0.01,sa=1.0,st=10.0):
        if t is None:
            t = self.tor1
        if a is None:
            a = self.angle1
        if d is None:
            d = self.d
        return abs(d-self.d)<sd and abs(a-self.angle1)<sa and abs(t-self.tor1)<st
    def get_dat(self):
        return [self.d, self.angle1, self.tor1]
    def get_reverse_dat(self):
        return [self.d, self.angle2, self.tor2]

class hbond_pdbase(pdbase):
    tables = [('hydrogen_bonds', HBreader(), '')]
    def insert_hb(self, code, item):
        self.insert_new_item('hydrogen_bonds', code, item)
    def get_hbonds(self, pdbcode=None):
        return self.get_items('hydrogen_bonds', pdbcode)
    def get_hbond_codetree(self):
        return self.get_items_codetree('hydrogen_bonds')
    def get_hbond_number(self):
        self.commit()
        return self.get_item_number('hydrogen_bonds')
    def delete_by_code(self, code):
        self.delete_items('hydrogen_bonds', code)
    def filter_same_residue(self, sameres=True):
        if sameres:
            self.execute('DELETE FROM hydrogen_bonds WHERE substr(resid1,2)!=substr(resid2,2)')
        else:
            self.execute('DELETE FROM hydrogen_bonds WHERE substr(resid1,2)=substr(resid2,2)')
    def filterby(self, *args, **kwds):
        if 'res2' in kwds:
            self.execute('DELETE FROM hydrogen_bonds WHERE resn2!=?',tuple([kwds['res2']]))
    def print_pvalues(self, hbtp, pvalue=0.05, fSym=False):
        for (code,hb) in self.get_hbonds():
            if hb.d and hb.angle1 and hb.tor1:
                pv = dat_pvalues(hbtp,hb.d, hb.angle1, hb.tor1)
                if (array(pv)>pvalue).any():
                    print("%5.2f %7.1f %7.1f " % (hb.d, hb.angle1, hb.tor1) + "%10.3g "*len(pv) % tuple(pv))
                if fSym:
                    if hb.angle2 and hb.tor2:
                        pv = dat_pvalues(hbtp,hb.d, hb.angle2, hb.tor2)
                        if (array(pv)>pvalue).any():
                            print("%5.2f %7.1f %7.1f " % (hb.d, hb.angle2, hb.tor2) + "%10.3g "*len(pv) % tuple(pv))
    def filter_pvalues(self, hbtp, pvalue=0.05, fSym=False):
        hbs = self.get_hbonds()
        print("Filter by p>%.3f" % (pvalue))
        print("Database contains %d hydrogen bonds" % self.get_hbond_number())
        self.delete_all('hydrogen_bonds')
        for (code,hb) in hbs:
            if hb.d and hb.angle1 and hb.tor1:
                pv = dat_pvalues(hbtp, hb.d, hb.angle1, hb.tor1)
                if (array(pv)>pvalue).any():
                    self.insert_hb(code,hb)
                    print("%5.2f %7.1f %7.1f " % (hb.d, hb.angle1, hb.tor1) + "%10.3g "*len(pv) % tuple(pv), end='\r')
                elif fSym:
                    if hb.angle2 and hb.tor2:
                        pv = dat_pvalues(hbtp,hb.d, hb.angle2, hb.tor2)
                        if (array(pv)>pvalue).any():
                            hb.swap()
                            self.insert_hb(code,hb)
                            print("%5.2f %7.1f %7.1f " % (hb.d, hb.angle1, hb.tor1) + "%10.3g "*len(pv) % tuple(pv), end='\r')
        print("%d hydrogen bonds passed the filter                    " % self.get_hbond_number())
        self.commit()
    def remove_metals(self, pisabase):
        hbs = self.get_hbond_codetree()
        print("Remove false bonds due to metal sites")
        print("Database contains %d hydrogen bonds" % self.get_hbond_number())
        self.delete_all('hydrogen_bonds')
        for code, items in hbs.items():
            print("Processing %s..." % code)
            molpath = pisabase.get_path(code)
            mol = ReadPDBfile(molpath[0])
            metals = mol.atom_lister('metal')
            print('Metals: ',len(metals))
            hbis = set(sum([[x.atomi1,x.atomi2] for x in items], []))
            nearmetals = set(mol.atom_lister('vicinity', listik=hbis, corelist=metals))
            print('Coordinated atoms: ',len(nearmetals))
            for item in items:
                if not len(nearmetals.intersection([item.atomi1,item.atomi2])):
                    self.insert_hb(code,item)
                else:
                    print("Removed %s %31s  " % (code, item.report()))
        print("%d hydrogen bonds passed the filter                    " % self.get_hbond_number())
        self.commit()

def read_contacts(fname, contact_label='OO_CONTACT', resreader=ACreader):
    with open(fname) as fin:
        lines = fin.readlines()
    contacts, params = [], {}
    for line in lines:
        try:
            key, value = line.split(':')[:2]
        except:
            sys.stderr.write(sys.argv[1]+':'+line+'\n')
            raise
        if key.strip() == contact_label:
            contacts.append(resreader(value))
        else:
            params[key.strip()] = value.strip()
    return contacts, params

class CysSGtoCysSG(AtomContact):
    ''' Disulfide bonds ''' 
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        AtomContact.__init__(self, model, *args, **kwds)
        self.store_rats(['CYSSG'], rcutoff=rcutoff, userkey='SSbonds')


class HydrogenBond(AtomContact):
    def report(self):
        for h,g,i,j,k,l,r,angle1,tor1,angle2,tor2 in self._aconts[self.__class__.__name__]:
            print("%-16s %-16s %8.2f %8.2f %8.2f %8.2f %8.2f" % (self._model.GetAtomTitle(i), self._model.GetAtomTitle(j), r, angle1,tor1,angle2,tor2))
    def bond_data(self):
        return  self._aconts[self.__class__.__name__]
    def evalbond(self, atomi, atomj):
        resindi = self._model.atom_lister('resid', resid=atomi.GetResID())
        resindj = self._model.atom_lister('resid', resid=atomj.GetResID())
        i = self._model.find(atomi, resindi)
        if i is None:
            return [float('nan'), float('nan'), float('nan'), float('nan')]
        j = self._model.find(atomj, resindj)
        if j is None:
            return [float('nan'), float('nan'), float('nan'), float('nan')]
        d = self._model.distance(i,j)
        g,h = self.get_gh(i)
        if h is None:
            return [d, float('nan'), float('nan'), float('nan')]
        elif g is None:
            a = self._model.angle(h,i,j)
            return [d, a, float('nan'), self.da_pvalue(d,a)]
        else:
            a = self._model.angle(h,i,j)
            return [d, a, self._model.torsion(g,h,i,j), self.da_pvalue(d,a)]
    def store_lists(self, listik1, listik2, key, rcutoff=4.0, sameres=False):
        AtomContact.store_lists(self, listik1, listik2, key, rcutoff=rcutoff, sameres=sameres)
        hbonds = []
        for (i,j,r) in self._aconts[self.__class__.__name__]:
            g, h = self.get_gh(i)
            if h is None:
                angle1, tor1 = float('nan'), float('nan')
            elif g is None:
                angle1, tor1 = self._model.angle(h,i,j),float('nan')
            else:
                angle1, tor1 = self._model.angle(h,i,j),self._model.torsion(g,h,i,j)
            k, l = self.get_kl(j)
            if k is None:
                angle2, tor2 = float('nan'), float('nan')
            elif l is None:
                angle2, tor2 = self._model.angle(k,j,i),float('nan')
            else:
                angle2, tor2 = self._model.angle(k,j,i),self._model.torsion(l,k,j,i)
            hbonds.append((h,g,i,j,k,l,r,angle1,tor1,angle2,tor2))
        self._aconts[self.__class__.__name__] = hbonds
    def get_gh(self, i):
        dummy=self._model.GetAtom(i).copy()
        resind=self._model.atom_lister('resid', resid=dummy.GetResID())
        nameh, nameg = self.get_namegh(dummy.GetName())
        dummy.set_name(nameh)
        h = self._model.find(dummy, resind)
        dummy.set_name(nameg)
        g = self._model.find(dummy, resind)
        return (g, h)
    def get_kl(self, j):
        dummy=self._model.GetAtom(j).copy()
        resind=self._model.atom_lister('resid', resid=dummy.GetResID())
        namek, namel = self.get_namekl(dummy.GetName())
        dummy.set_name(namek)
        k = self._model.find(dummy, resind)
        dummy.set_name(namel)
        l = self._model.find(dummy, resind)
        return (k, l)
    def get_namegh(self, namei):
        pass
    def get_namekl(self, namej):
        pass
    def da_pvalue(self, d, a):
        return da_pvalue(self.__class__.__name__, d, a)
    def pfilter(self, p):
        self._aconts[self.__class__.__name__] = [x for x in self._aconts[self.__class__.__name__] if self.da_pvalue(x[6],x[7])>p]
    def resid_filter(self, resids):
        self._aconts[self.__class__.__name__] = [x for x in self._aconts[self.__class__.__name__] if self._model.GetAtom(x[2]).resid() in resids or self._model.GetAtom(x[3]).resid() in resids]
    def report_pval(self):
        for h,g,i,j,k,l,r,angle1,tor1,angle2,tor2 in self._aconts[self.__class__.__name__]:
            print("%-16s %-16s %8.2f %8.2f %8.2f %8.2f" % (self._model.GetAtomTitle(i), self._model.GetAtomTitle(j), r, angle1,tor1,self.da_pvalue(r,angle1)))
    def get_hbonds(self):
        hbonds = {}
        for h,g,i,j,k,l,r,angle1,tor1,angle2,tor2 in self._aconts[self.__class__.__name__]:
            hbonds[self._model.GetAtomTitle(i)+'-'+self._model.GetAtomTitle(j)] = (i,j,r,angle1,tor1,self.da_pvalue(r,angle1))
        return hbonds
    def get_item_sets(self):
        itemsets = []
        for h,g,i,j,k,l,r,angle1,tor1,angle2,tor2 in self._aconts[self.__class__.__name__]:
            itemsets.append([self._model.GetAtomResidueName(i), self._model.GetAtomResID(i), self._model.GetAtomName(i),
                             self._model.GetAtomResidueName(j), self._model.GetAtomResID(j), self._model.GetAtomName(j),
                             r, angle1, tor1, angle2, tor2, i, j])
        return itemsets
    def get_readers(self, reader=HBreader):
        return [reader(x) for x in self.get_item_sets()]
    def get_hbids(self):
        return [self._model.GetAtomTitle(x[2])+'-'+self._model.GetAtomTitle(x[3]) for x in self._aconts[self.__class__.__name__]]
    def conserved_hbonds(self, other):
        return sorted(list(set(self.get_hbids()).intersection(other.get_hbids())))
    def broken_hbonds(self, other):
        return sorted(list(set(self.get_hbids()).difference(other.get_hbids())))
    def new_hbonds(self, other):
        return sorted(list(set(other.get_hbids()).difference(self.get_hbids())))
    def report_diffs(self, other):
        selfbonds = self.get_hbonds()
        otherbonds = other.get_hbonds()
        hbs = self.conserved_hbonds(other)
        print("Conserved hydrogen bonds (%d)" % (len(hbs)))
        if len(hbs):
            for name in hbs:
                print("%-20s %-20s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" % tuple(name.split('-') + list(selfbonds[name][2:]) + list(otherbonds[name][2:])))
        hbs = self.broken_hbonds(other)
        print("Broken hydrogen bonds (%d)" % (len(hbs)))
        if len(hbs):
            for name in hbs:
                atomi = self._model.GetAtom(selfbonds[name][0])
                atomj = self._model.GetAtom(selfbonds[name][1])
                print("%-20s %-20s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" % tuple(name.split('-') + list(selfbonds[name][2:]) + other.evalbond(atomi, atomj)))
        hbs = self.new_hbonds(other)
        print("Newly formed hydrogen bonds (%d)" % (len(hbs)))
        if len(hbs):
            for name in hbs:
                atomi = other._model.GetAtom(otherbonds[name][0])
                atomj = other._model.GetAtom(otherbonds[name][1])
                print("%-20s %-20s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" % tuple(name.split('-') + list(otherbonds[name][2:]) + self.evalbond(atomi, atomj)))

class X2YHBond(HydrogenBond):
    '''  Generic hydrogen bond class defined via residue types.
        Following parameters have to be provided:
            rats1 - donor atom types
            rats2 - acceptor atom types
            sameres - whether to exclude intra-residue bonds 
            namegh - get_namegh method reference 
            namekl - get_namekl method reference
    '''
    def __init__(self, model, rats1, rats2, sameres, namegh, namekl, rcutoff=3.2, *args, **kwds):
        self.namegh = namegh
        self.namekl = namekl
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=rats1)
        listik2 = model.atom_lister('rat', rats=rats2)
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=sameres)
    def get_namegh(self, namei):
        return self.namegh(namei)
    def get_namekl(self, namej):
        return self.namekl(namej)

class X2bbHBond(HydrogenBond):
    '''  Generic hydrogen bond class with backbone oxygen acceptor defined via donor residue type.
        Following parameters must be provided:
            rats - donor atom types
            sameres - whether to exclude intra-residue bonds 
            namegh - get_namegh method reference
    '''
    def __init__(self, model, rats, sameres, namegh, rcutoff=3.2, *args, **kwds):
        self.namegh = namegh
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=rats)
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=sameres)
    def get_namegh(self, namei):
        return self.namegh(namei)
    def get_namekl(self, namej):
        return get_bbo_names()

class bb2YHBond(HydrogenBond):
    '''  Generic hydrogen bond class with backbone nitrogen donor defined via acceptor residue type.
        Following parameters have to be provided:
            rats - acceptor atom types
            sameres - whether to exclude intra-residue bonds 
            namekl - get_namekl method reference
    '''
    def __init__(self, model, rats, sameres, namekl, rcutoff=3.2, *args, **kwds):
        self.namekl = namekl
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('name_N', listik=model.atom_lister('p_bb'))
        listik2 = model.atom_lister('rat', rats=rats)
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=sameres)
    def get_gh(self, i):
        c_list = self._model.atom_lister('vicinity', listik=self._model.atom_lister('name_C', listik=self._model.atom_lister('p_bb')), corelist=[i], rcutoff=2.0)
        if len(c_list):
            h = c_list[0]
            dummy=self._model.GetAtom(h).copy()
            dummy.set_name(' CA ')
            g = self._model.find(dummy, self._model.atom_lister('resid', resid=dummy.GetResID()))
            return (g, h)
        return (None, None)
    def get_namekl(self, namej):
        return self.namekl(namej)



# Tyrosine donor classes

class TyrOHtoAnionic(X2YHBond):
    ''' Tyr OH donating hydrogen to a Glu/Asp side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'], sameres=True, namegh=get_tyr_names, namekl=get_de_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoAspOD(X2YHBond):
    ''' Tyr OH donating hydrogen to an Asp side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['ASPOD1', 'ASPOD2'], sameres=True, namegh=get_tyr_names, namekl=get_asp_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoGluOE(X2YHBond):
    ''' Tyr OH donating hydrogen to a Glu side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['GLUOE1', 'GLUOE2'], sameres=True, namegh=get_tyr_names, namekl=get_glu_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoAsnOD1(X2YHBond):
    ''' Tyr OG1 donating hydrogen to a Asn side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['ASNOD1'], sameres=True, namegh=get_tyr_names, namekl=get_asn_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoGlnOE1(X2YHBond):
    ''' Tyr OG1 donating hydrogen to a Gln side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['GLNOE1'], sameres=True, namegh=get_tyr_names, namekl=get_gln_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoCarboxamideO(X2YHBond):
    ''' Tyr OG1 donating hydrogen to a carboxamide (N/Q) side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['GLNOE1', 'ASNOD1'], sameres=True, namegh=get_tyr_names, namekl=get_nqo_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoTyrOH(X2YHBond):
    ''' Tyrosine OG1 donating hydrogen to a Tyrosine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['TYROH'], sameres=False, namegh=get_tyr_names, namekl=get_tyr_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoSerOG(X2YHBond):
    ''' Tyrosine OG1 donating hydrogen to a Serine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['SEROG'], sameres=True, namegh=get_tyr_names, namekl=get_ser_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoThrOG1(X2YHBond):
    ''' Tyrosine OG1 donating hydrogen to a Threonine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['THROG1'], sameres=True, namegh=get_tyr_names, namekl=get_thr_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoHis(X2YHBond):
    ''' Tyrosine OG1 donating hydrogen to a Histidine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['HISND1','HISNE2'], sameres=True, namegh=get_tyr_names, namekl=get_his_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoCysSG(X2YHBond):
    ''' Tyrosine OG1 donating hydrogen to a Cysteine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TYROH'], rats2=['CYSSG'], sameres=True, namegh=get_tyr_names, namekl=get_cys_names, rcutoff=rcutoff, *args, **kwds)
class TyrOHtoBackbone(X2bbHBond):
    ''' Tyr hydroxyl to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2bbHBond.__init__(self, model, rats=['TYROH'], sameres=True, namegh=get_tyr_names, rcutoff=rcutoff, *args, **kwds)

# Tryptophan donor classes

class TrpNE1toAnionic(X2YHBond):
    ''' Tryptophan donating hydrogen to an anionic (Glu, Asp) side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'], sameres=True, namegh=get_trp_names, namekl=get_de_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toAspOD(X2YHBond):
    ''' Tryptophan donating hydrogen to an Aspartate side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['ASPOD1', 'ASPOD2'], sameres=True, namegh=get_trp_names, namekl=get_asp_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toGluOE(X2YHBond):
    ''' Tryptophan donating hydrogen to an Aspartate side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['GLUOE1', 'GLUOE2'], sameres=True, namegh=get_trp_names, namekl=get_glu_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toCarboxamideO(X2YHBond):
    ''' Tryptophan donating hydrogen to a carboxamide (N-H-O) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['GLNOE1', 'ASNOD1'], sameres=True, namegh=get_trp_names, namekl=get_nqo_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toAsnOD1(X2YHBond):
    ''' Tryptophan donating hydrogen to an asparagine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['ASNOD1'], sameres=True, namegh=get_trp_names, namekl=get_nqo_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toGlnOE1(X2YHBond):
    ''' Tryptophan donating hydrogen to a glutamine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['GLNOE1'], sameres=True, namegh=get_trp_names, namekl=get_nqo_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toSerOG(X2YHBond):
    ''' Tryptophan donating hydrogen to a Serine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['SEROG'], sameres=True, namegh=get_trp_names, namekl=get_ser_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toThrOG1(X2YHBond):
    ''' Tryptophan donating hydrogen to a Threonine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['THROG1'], sameres=True, namegh=get_trp_names, namekl=get_thr_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toTyrOH(X2YHBond):
    ''' Tryptophan donating hydrogen to a Tyrosine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['TYROH'], sameres=True, namegh=get_trp_names, namekl=get_tyr_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toHydroxyl(X2YHBond):
    ''' Tryptophan donating hydrogen to a Tyr/Ser/Thr side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['TYROH','THROG1','SEROG'], sameres=True, namegh=get_trp_names, namekl=get_sty_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toHis(X2YHBond):
    ''' Tryptophan donating hydrogen to a histidine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['TRPNE1'], rats2=['HISND1','HISNE2'], sameres=True, namegh=get_trp_names, namekl=get_his_names, rcutoff=rcutoff, *args, **kwds)
class TrpNE1toBackbone(X2bbHBond):
    ''' Tryptophan nitrogen (donor) to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2bbHBond.__init__(self, model, rats=['TRPNE1'], sameres=True, namegh=get_trp_names, rcutoff=rcutoff, *args, **kwds)

# Backbone donor classes

class BackboneNtoCysSG(bb2YHBond):
    ''' Backbone nitrogen (donor) to Cys sulfur acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['CYSSG'], sameres=False, namekl=get_cys_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoAnionic(bb2YHBond):
    ''' Backbone nitrogen (donor) to Asp/Glu oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'], sameres=False, namekl=get_de_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoAspOD(bb2YHBond):
    ''' Backbone nitrogen (donor) to Asp oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['ASPOD1', 'ASPOD2'], sameres=False, namekl=get_asp_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoGluOE(bb2YHBond):
    ''' Backbone nitrogen (donor) to Glu oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['GLUOE1', 'GLUOE2'], sameres=False, namekl=get_glu_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoCarboxamideO(bb2YHBond):
    ''' Backbone nitrogen (donor) to carboxamide (Asn/Gln) oxygen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['GLNOE1', 'ASNOD1'], sameres=False, namekl=get_nqo_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoAsnOD1(bb2YHBond):
    ''' Backbone nitrogen (donor) to carboxamide Asn oxygen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['ASNOD1'], sameres=False, namekl=get_asn_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoGlnOE1(bb2YHBond):
    ''' Backbone nitrogen (donor) to carboxamide Gln oxygen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['GLNOE1'], sameres=False, namekl=get_gln_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoHydroxyl(bb2YHBond):
    ''' Backbone nitrogen (donor) to Tyr/Ser/Thr side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['TYROH', 'SEROG', 'THROG1'], sameres=False, namekl=get_sty_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoSerOG(bb2YHBond):
    ''' Backbone nitrogen (donor) to Ser side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['SEROG'], sameres=False, namekl=get_ser_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoThrOG1(bb2YHBond):
    ''' Backbone nitrogen (donor) to Thr side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['THROG1'], sameres=False, namekl=get_thr_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoTyrOH(bb2YHBond):
    ''' Backbone nitrogen (donor) to Tyr side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['TYROH'], sameres=False, namekl=get_tyr_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoHis(bb2YHBond):
    ''' Backbone nitrogen (donor) to histidine nitrogen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['HISND1','HISNE2'], sameres=False, namekl=get_his_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoHisNE2(bb2YHBond):
    ''' Backbone nitrogen (donor) to histidine NE2 acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['HISNE2'], sameres=False, namekl=get_his_names, rcutoff=rcutoff, *args, **kwds)
class BackboneNtoHisND1(bb2YHBond):
    ''' Backbone nitrogen (donor) to histidine ND1 acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        bb2YHBond.__init__(self, model, rats=['HISND1'], sameres=False, namekl=get_his_names, rcutoff=rcutoff, *args, **kwds)

class BackboneNHydrogenBond(HydrogenBond):
    def get_gh(self, i):
        c_list = self._model.atom_lister('vicinity', listik=self._model.atom_lister('name_C', listik=self._model.atom_lister('p_bb')), corelist=[i], rcutoff=2.0)
        if len(c_list):
            h = c_list[0]
            dummy=self._model.GetAtom(h).copy()
            dummy.set_name(' CA ')
            g = self._model.find(dummy, self._model.atom_lister('resid', resid=dummy.GetResID()))
            return (g, h)
        return (None, None)

class BackboneNtoBackboneO(BackboneNHydrogenBond):
    ''' Backbone nitrogen (donor) to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        bblistik = model.atom_lister('p_bb')
        listik1 = model.atom_lister('name_N', listik=bblistik)
        listik2 = model.atom_lister('name_O', listik=bblistik)
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff)
        self._aconts[self.__class__.__name__] = [x for x in self._aconts[self.__class__.__name__] if self._model.distance(x[2],x[4])>2.0]
    def get_namekl(self, namej):
        return get_bbo_names()

# Histidine donor classes

class HistoHis(X2YHBond):
    ''' Histidine nitrogen (donor) to histidine nitrogen (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['HISND1','HISNE2'], rats2=['HISND1','HISNE2'], sameres=False, namegh=get_his_names, namekl=get_his_names, rcutoff=rcutoff, *args, **kwds)
class HistoAnionic(X2YHBond):
    ''' Histidine nitrogen (donor) to Asp/Glu oxygen (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['HISND1','HISNE2'], rats2=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'], sameres=True, namegh=get_his_names, namekl=get_de_names, rcutoff=rcutoff, *args, **kwds)
class HistoCarboxamideO(X2YHBond):
    ''' Histidine nitrogen (donor) to Gln/Asn oxygen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['HISND1','HISNE2'], rats2=['GLNOE1', 'ASNOD1'], sameres=True, namegh=get_his_names, namekl=get_nqo_names, rcutoff=rcutoff, *args, **kwds)
class HistoHydroxyl(X2YHBond):
    ''' Histidine nitrogen (donor) to Ser/Thr/Tyr hydroxyl '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['HISND1','HISNE2'], rats2=['TYROH', 'SEROG', 'THROG1'], sameres=True, namegh=get_his_names, namekl=get_sty_names, rcutoff=rcutoff, *args, **kwds)
class HistoCysSG(X2YHBond):
    ''' Histidine nitrogen (donor) to Cysteine sulfur (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2YHBond.__init__(self, model, rats1=['HISND1','HISNE2'], rats2=['CYSSG'], sameres=True, namegh=get_his_names, namekl=get_cys_names, rcutoff=rcutoff, *args, **kwds)
class HistoBackbone(X2bbHBond):
    ''' Histidine nitrogen (donor) to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        X2bbHBond.__init__(self, model, rats=['HISND1','HISNE2'], sameres=False, namegh=get_his_names, rcutoff=rcutoff, *args, **kwds)

class SerOGtoAspOD(AtomContact):
    ''' Ser OG donating hydrogen to an Asp side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        AtomContact.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)

class SerOGtoGluOE(AtomContact):
    ''' Ser OG donating hydrogen to a Glu side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        AtomContact.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('rat', rats=['GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)

class SerOGtoAnionic(HydrogenBond):
    ''' Ser OG donating hydrogen to a Glu/Asp side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        if namej in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namej in ['OE1','OE2']:
            return ' CD ', ' CG '

class AnionictoSerOG(HydrogenBond):
    ''' Glu/Asp side chain donating hydrogen to a Ser OG '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        listik2 = model.atom_lister('rat', rats=['SEROG'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        if namei in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namei in ['OE1','OE2']:
            return ' CD ', ' CG '
    def get_namekl(self, namej):
        return ' CB ', ' CA '
    
class ThrOG1toAspOD(AtomContact):
    ''' Thr OG1 donating hydrogen to an Asp side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        AtomContact.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)

class ThrOG1toGluOE(AtomContact):
    ''' Thr OG1 donating hydrogen to a Glu side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        AtomContact.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('rat', rats=['GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)

class ThrOG1toAnionic(HydrogenBond):
    ''' Thr OG1 donating hydrogen to a Glu/Asp side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        if namej in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namej in ['OE1','OE2']:
            return ' CD ', ' CG '

class AnionictoThrOG1(HydrogenBond):
    ''' Glu/Asp side chain donating hydrogen to a Thr OG1 '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        listik2 = model.atom_lister('rat', rats=['THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        if namei in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namei in ['OE1','OE2']:
            return ' CD ', ' CG '
    def get_namekl(self, namej):
        return ' CB ', ' CA '
    
class AnionictoTyrOH(HydrogenBond):
    ''' Glu/Asp side chain donating hydrogen to a Tyr OH '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        listik2 = model.atom_lister('rat', rats=['TYROH'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        if namei in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namei in ['OE1','OE2']:
            return ' CD ', ' CG '
    def get_namekl(self, namej):
        return ' CZ ', ' CE1'
    
class HydroxylToAnionic(HydrogenBond):
    ''' Hydroxyl (Tyr, Ser, Thr) donating hydrogen to an anionic (Glu, Asp) side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        if namei == 'OH':
            return ' CZ ', ' CE1'
        if namei in ['OG', 'OG1']:
            return ' CB ', ' CA '
    def get_namekl(self, namej):
        if namej in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namej in ['OE1','OE2']:
            return ' CD ', ' CG '

class AnionicToAnionic(HydrogenBond):
    ''' Carboxylic acid (protonated Glu, Asp) donating hydrogen to a carboxylate (Glu, Asp) (Rare!) ''' 
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff)
    def get_namegh(self, namei):
        return get_de_names(namei)
    def get_namekl(self, namej):
        return get_de_names(namej)

class AnionicToHydroxyl(HydrogenBond):
    ''' Anionic (Glu, Asp) side chain donating hydrogen to a hydroxyl (Tyr, Ser, Thr) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        listik2 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        if namei in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namei in ['OE1','OE2']:
            return ' CD ', ' CG '
    def get_namekl(self, namej):
        if namej == 'OH':
            return ' CZ ', ' CE1'
        if namej in ['OG', 'OG1']:
            return ' CB ', ' CA '

class SerOGtoAsnOD1(HydrogenBond):
    ''' Ser OG donating hydrogen to a Asn side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('rat', rats=['ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CG ', ' CB '

class SerOGtoGlnOE1(HydrogenBond):
    ''' Ser OG donating hydrogen to a Gln side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CD ', ' CG '

class SerOGtoCarboxamideO(HydrogenBond):
    ''' Ser OG donating hydrogen to a carboxamide (N/Q) side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1', 'ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        if namej in ['OD1']:
            return ' CG ', ' CB '
        if namej in ['OE1']:
            return ' CD ', ' CG '

class ThrOG1toAsnOD1(HydrogenBond):
    ''' Thr OG1 donating hydrogen to a Asn side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('rat', rats=['ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CG ', ' CB '

class ThrOG1toGlnOE1(HydrogenBond):
    ''' Thr OG1 donating hydrogen to a Gln side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CD ', ' CG '

class ThrOG1toCarboxamideO(HydrogenBond):
    ''' Thr OG1 donating hydrogen to a carboxamide (N/Q) side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1', 'ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        if namej in ['OD1']:
            return ' CG ', ' CB '
        if namej in ['OE1']:
            return ' CD ', ' CG '

class HydroxylToCarboxamideO(HydrogenBond):
    ''' Hydroxyl (Tyr, Ser, Thr) donating hydrogen to a carboxamide (N/Q) side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1', 'ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        if namei == 'OH':
            return ' CZ ', ' CE1'
        if namei in ['OG', 'OG1']:
            return ' CB ', ' CA '
    def get_namekl(self, namej):
        if namej in ['OD1']:
            return ' CG ', ' CB '
        if namej in ['OE1']:
            return ' CD ', ' CG '

class CarboxamideNtoCarboxamideO(HydrogenBond):
    ''' Carboxamide donating hydrogen to a carboxamide (N-H-O) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2', 'ASNND2'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1', 'ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff)
    def get_namegh(self, namei):
        if namei == 'ND2':
            return ' CG ', ' CB '
        if namei == 'NE2':
            return ' CD ', ' CG '
    def get_namekl(self, namej):
        if namej in ['OD1']:
            return ' CG ', ' CB '
        if namej in ['OE1']:
            return ' CD ', ' CG '

class CarboxamideNToHis(HydrogenBond):
    '''  Asn/Gln carboxamide nitrogen (donor) to histidine nitrogen (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2', 'ASNND2'])
        listik2 = model.atom_lister('rat', rats=['HISND1','HISNE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_nqn_names(namei)
    def get_namekl(self, namej):
        return get_his_names(namej)

class CarboxamideNToAnionic(HydrogenBond):
    '''  Carboxamide donating hydrogen to an anionic (Glu, Asp) side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2', 'ASNND2'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        if namei == 'ND2':
            return ' CG ', ' CB '
        if namei == 'NE2':
            return ' CD ', ' CG '
    def get_namekl(self, namej):
        if namej in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namej in ['OE1','OE2']:
            return ' CD ', ' CG '

class CarboxamideNToHydroxyl(HydrogenBond):
    '''  Asn/Gln carboxamide nitrogen (donor) to Ser/Thr/Tyr hydroxyl oxygen (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2', 'ASNND2'])
        listik2 = model.atom_lister('rat', rats=['SEROG', 'THROG1', 'TYROH'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        if namei == 'ND2':
            return ' CG ', ' CB '
        if namei == 'NE2':
            return ' CD ', ' CG '
    def get_namekl(self, namej):
        if namej == 'OH':
            return ' CZ ', ' CE1'
        if namej in ['OG', 'OG1']:
            return ' CB ', ' CA '

class CarboxamideNToTyrOH(HydrogenBond):
    '''  Asn/Gln carboxamide nitrogen (donor) to Tyr hydroxyl oxygen (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2', 'ASNND2'])
        listik2 = model.atom_lister('rat', rats=['TYROH'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_nqn_names(namei)
    def get_namekl(self, namej):
        return get_tyr_names(namej)

class CarboxamideNToThrOG1(HydrogenBond):
    '''  Asn/Gln carboxamide nitrogen (donor) to Thr hydroxyl oxygen (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2', 'ASNND2'])
        listik2 = model.atom_lister('rat', rats=['THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_nqn_names(namei)
    def get_namekl(self, namej):
        return get_thr_names(namej)

class CarboxamideNToSerOG(HydrogenBond):
    '''  Asn/Gln carboxamide nitrogen (donor) to Ser hydroxyl oxygen (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2', 'ASNND2'])
        listik2 = model.atom_lister('rat', rats=['SEROG'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_nqn_names(namei)
    def get_namekl(self, namej):
        return get_ser_names(namej)

class CarboxamideNtoCysSG(HydrogenBond):
    '''  Asn/Gln carboxamide nitrogen (donor)to Cysteine sulfur (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2', 'ASNND2'])
        listik2 = model.atom_lister('rat', rats=['CYSSG'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_nqn_names(namei)
    def get_namekl(self, namej):
        return get_cys_names(namej)

class SerOGtoSerOG(HydrogenBond):
    ''' Ser OG donating hydrogen to a Serine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('rat', rats=['SEROG'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CB ', ' CA '

class ThrOG1toThrOG1(HydrogenBond):
    ''' Thr OG1 donating hydrogen to a Threonine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('rat', rats=['THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CB ', ' CA '

class SerOGtoThrOG1(HydrogenBond):
    ''' Ser OG donating hydrogen to a Threonine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('rat', rats=['THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CB ', ' CA '

class ThrOG1toSerOG(HydrogenBond):
    ''' Thr oG1 donating hydrogen to a Serine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('rat', rats=['SEROG'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CB ', ' CA '

class SerOGtoTyrOH(HydrogenBond):
    ''' Ser OG donating hydrogen to a Tyrosine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('rat', rats=['TYROH'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CZ ', ' CE1'

class ThrOG1toTyrOH(HydrogenBond):
    ''' Thr OG1 donating hydrogen to a Tyrosine side chain '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('rat', rats=['TYROH'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CB ', ' CA '
    def get_namekl(self, namej):
        return ' CZ ', ' CE1'

class LysNZtoAnionic(HydrogenBond):
    ''' Lysine amine nitrogen (donor) to Glu/Asp hydroxyl acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['LYSNZ'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CE ', ' CD '
    def get_namekl(self, namej):
        if namej in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namej in ['OE1','OE2']:
            return ' CD ', ' CG '

class LysNZtoCarboxamideO(HydrogenBond):
    ''' Lysine amine nitrogen (donor) to Gln/Asn oxygen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['LYSNZ'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1', 'ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CE ', ' CD '
    def get_namekl(self, namej):
        if namej == 'OD1':
            return ' CG ', ' CB '
        if namej == 'OE1':
            return ' CD ', ' CG '

class LysNZtoHydroxyl(HydrogenBond):
    ''' Lysine amine nitrogen (donor) to Ser/Thr/Tyr hydroxyl '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['LYSNZ'])
        listik2 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_lys_names()
    def get_namekl(self, namej):
        return get_sty_names(namej)

class LysNZtoSerOG(HydrogenBond):
    ''' Lysine amine nitrogen (donor) to Ser hydroxyl '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['LYSNZ'])
        listik2 = model.atom_lister('rat', rats=['SEROG'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_lys_names()
    def get_namekl(self, namej):
        return get_ser_names(namej)

class LysNZtoThrOG1(HydrogenBond):
    ''' Lysine amine nitrogen (donor) to Thr hydroxyl '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['LYSNZ'])
        listik2 = model.atom_lister('rat', rats=['THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_lys_names()
    def get_namekl(self, namej):
        return get_thr_names(namej)

class LysNZtoCysSG(HydrogenBond):
    ''' Lysine amine nitrogen (donor) to Cysteine sulfur (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['LYSNZ'])
        listik2 = model.atom_lister('rat', rats=['CYSSG'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_lys_names()
    def get_namekl(self, namej):
        return get_cys_names(namej)

class LysNZtoHis(HydrogenBond):
    ''' Lysine amine nitrogen (donor) to Histidine '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['LYSNZ'])
        listik2 = model.atom_lister('rat', rats=['HISND1','HISNE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_lys_names()
    def get_namekl(self, namej):
        return get_his_names(namej)

class ArgNEtoAnionic(HydrogenBond):
    ''' Arginine amide (donor) to Glu/Asp hydroxyl acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNE'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CD ', ' CG '
    def get_namekl(self, namej):
        if namej in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namej in ['OE1','OE2']:
            return ' CD ', ' CG '

class ArgNEtoCarboxamideO(HydrogenBond):
    ''' Arginine amide (donor) to Gln/Asn oxygen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNE'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1', 'ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CD ', ' CG '
    def get_namekl(self, namej):
        if namej == 'OD1':
            return ' CG ', ' CB '
        if namej == 'OE1':
            return ' CD ', ' CG '

class ArgNEtoHydroxyl(HydrogenBond):
    ''' Arginine amide (donor) to Ser/Thr/Tyr hydroxyl '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNE'])
        listik2 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CD ', ' CG '
    def get_namekl(self, namej):
        if namej == 'OH':
            return ' CZ ', ' CE1'
        if namej in ['OG', 'OG1']:
            return ' CB ', ' CA '

class ArgAminetoAnionic(HydrogenBond):
    ''' Arginine amine (donor) to Glu/Asp hydroxyl acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNH1','ARGNH2'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CZ ', ' NE '
    def get_namekl(self, namej):
        if namej in ['OD1','OD2']:
            return ' CG ', ' CB '
        if namej in ['OE1','OE2']:
            return ' CD ', ' CG '

class ArgAminetoCarboxamideO(HydrogenBond):
    ''' Arginine amine (donor) to Gln/Asn oxygen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNH1','ARGNH2'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1', 'ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CZ ', ' NE '
    def get_namekl(self, namej):
        if namej == 'OD1':
            return ' CG ', ' CB '
        if namej == 'OE1':
            return ' CD ', ' CG '

class ArgAminetoHydroxyl(HydrogenBond):
    ''' Arginine amine (donor) to Ser/Thr/Tyr hydroxyl '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNH1','ARGNH2'])
        listik2 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return ' CZ ', ' NE '
    def get_namekl(self, namej):
        if namej == 'OH':
            return ' CZ ', ' CE1'
        if namej in ['OG', 'OG1']:
            return ' CB ', ' CA '

class ArgtoAnionic(HydrogenBond):
    ''' Arginine guanidinium (donor) to Glu/Asp oxygen (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNE','ARGNH1','ARGNH2'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_arg_names(namei)
    def get_namekl(self, namej):
        return get_de_names(namej)

class ArgtoCarboxamideO(HydrogenBond):
    ''' Arginine guanidinium (donor) to Gln/Asn oxygen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNE','ARGNH1','ARGNH2'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1', 'ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_arg_names(namei)
    def get_namekl(self, namej):
        if namej == 'OD1':
            return ' CG ', ' CB '
        if namej == 'OE1':
            return ' CD ', ' CG '

class ArgtoHydroxyl(HydrogenBond):
    ''' Arginine guanidinium (donor) to Ser/Thr/Tyr hydroxyl '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNE', 'ARGNH1','ARGNH2'])
        listik2 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_arg_names(namei)
    def get_namekl(self, namej):
        if namej == 'OH':
            return ' CZ ', ' CE1'
        if namej in ['OG', 'OG1']:
            return ' CB ', ' CA '

class ArgtoBackbone(HydrogenBond):
    ''' Arginine guanidinium (donor) to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNE', 'ARGNH1','ARGNH2'])
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_arg_names(namei)
    def get_namekl(self, namej):
        return get_bbo_names()

class ArgtoCysSG(HydrogenBond):
    ''' Arginine guanidinium (donor) to Cysteine sulfur (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNE', 'ARGNH1','ARGNH2'])
        listik2 = model.atom_lister('rat', rats=['CYSSG'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_arg_names(namei)
    def get_namekl(self, namej):
        return get_cys_names(namej)

class ArgToHis(HydrogenBond):
    ''' Arginine guanidinium (donor) to histidine nitrogen (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ARGNE', 'ARGNH1','ARGNH2'])
        listik2 = model.atom_lister('rat', rats=['HISND1','HISNE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_arg_names(namei)
    def get_namekl(self, namej):
        return get_his_names(namej)


class LystoBackbone(HydrogenBond):
    ''' Lysine nitrogen (donor) to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['LYSNZ'])
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_lys_names()
    def get_namekl(self, namej):
        return get_bbo_names()

class HydroxylToBackbone(HydrogenBond):
    ''' Ser/Thr/Tyr hydroxyl to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_sty_names(namei)
    def get_namekl(self, namej):
        return get_bbo_names()

class HydroxylToCysSG(HydrogenBond):
    '''  Ser/Thr/Tyr hydroxyl (donor) to Cysteine sulfur (acceptor) '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        listik2 = model.atom_lister('rat', rats=['CYSSG'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_sty_names(namei)
    def get_namekl(self, namej):
        return get_cys_names(namej)

class SerOGtoBackbone(HydrogenBond):
    ''' Ser hydroxyl to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['SEROG'])
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_ser_names(namei)
    def get_namekl(self, namej):
        return get_bbo_names()

class ThrOG1toBackbone(HydrogenBond):
    ''' Thr hydroxyl to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['THROG1'])
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_thr_names(namei)
    def get_namekl(self, namej):
        return get_bbo_names()

class CarboxamideNtoBackbone(HydrogenBond):
    ''' Asn/Gln nitrogen (donor) to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2', 'ASNND2'])
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_nqn_names(namei)
    def get_namekl(self, namej):
        return get_bbo_names()

class AsnNtoBackbone(HydrogenBond):
    ''' Asn nitrogen (donor) to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['ASNND2'])
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_asn_names(namei)
    def get_namekl(self, namej):
        return get_bbo_names()

class GlnNtoBackbone(HydrogenBond):
    ''' Gln nitrogen (donor) to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['GLNNE2'])
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_gln_names(namei)
    def get_namekl(self, namej):
        return get_bbo_names()

class CysSGtoBackbone(HydrogenBond):
    ''' Cysteine sulfur (donor) to backbone oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['CYSSG'])
        listik2 = model.atom_lister('name_O', listik=model.atom_lister('p_bb'))
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff)
    def get_namegh(self, namei):
        return get_cys_names(namei)
    def get_namekl(self, namej):
        return get_bbo_names()

class CysSGtoAnionic(HydrogenBond):
    ''' Cysteine sulfur (donor) to Asp/Glu oxygen '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['CYSSG'])
        listik2 = model.atom_lister('rat', rats=['ASPOD1', 'ASPOD2', 'GLUOE1', 'GLUOE2'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_cys_names(namei)
    def get_namekl(self, namej):
        return get_de_names(namej)

class CysSGCarboxamideO(HydrogenBond):
    ''' Cysteine sulfur (donor) to Gln/Asn oxygen acceptor '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds)
        listik1 = model.atom_lister('rat', rats=['CYSSG'])
        listik2 = model.atom_lister('rat', rats=['GLNOE1', 'ASNOD1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_cys_names(namei)
    def get_namekl(self, namej):
        return get_nqo_names(namej)

class CysSGtoHydroxyl(HydrogenBond):
    ''' Cysteine sulfur (donor) to Ser/Thr/Tyr hydroxyl '''
    def __init__(self, model, rcutoff=3.2, *args, **kwds):
        HydrogenBond.__init__(self, model, *args, **kwds) 
        listik1 = model.atom_lister('rat', rats=['CYSSG'])
        listik2 = model.atom_lister('rat', rats=['TYROH', 'SEROG', 'THROG1'])
        self.store_lists(listik1, listik2, self.__class__.__name__, rcutoff=rcutoff, sameres=True)
    def get_namegh(self, namei):
        return get_cys_names(namei)
    def get_namekl(self, namej):
        return get_sty_names(namej)

def get_cys_names(name):
    return ' CB ', ' CA '
def get_bbo_names():
    return ' C  ', ' CA '
def get_bbn_names():
    pass
def get_lys_names():
    return ' CE ', ' CD '
def get_his_names(name):
    if name == 'ND1':
        return ' CG ', ' CB '
    if name == 'NE2':
        return ' CD2', ' CG '
def get_tyr_names(name):
    return ' CZ ', ' CE1'
def get_trp_names(name):
    return ' CD1', ' CG '
def get_thr_names(name):
    return ' CB ', ' CA '
def get_ser_names(name):
    return ' CB ', ' CA '
def get_sty_names(name):
    if name == 'OH':
        return ' CZ ', ' CE1'
    if name in ['OG', 'OG1']:
        return ' CB ', ' CA '
def get_nqo_names(name):
    if name == 'OD1':
        return ' CG ', ' CB '
    if name == 'OE1':
        return ' CD ', ' CG '
def get_nqn_names(name):
    if name == 'ND2':
        return ' CG ', ' CB '
    if name == 'NE2':
        return ' CD ', ' CG '
def get_asn_names(name):
    return  ' CG ', ' CB '
def get_gln_names(name):
    return ' CD ', ' CG '
def get_de_names(name):
    if name in ['OD1','OD2']:
        return ' CG ', ' CB '
    if name in ['OE1','OE2']:
        return ' CD ', ' CG '
def get_asp_names(name):
    return ' CG ', ' CB '
def get_glu_names(name):
    return ' CD ', ' CG '
def get_arg_names(name):
    if name == 'NE':
        return ' CD ', ' CG '
    if name[:2] == 'NH':
        return ' CZ ', ' NE '
