'''
Module for reading PDB-files.  
Includes pdbatom and pdbmolecule classes
'''

import gzip, urllib.request, os, random, math, sys, re, copy, logging, time

import pdbnames, SpaceGroups
from helper import progressbar
from rotate import transform_list
from tinertia import TInertia
from scipy.linalg import eigh
from scipy import   array, cos, sin, pi, radians, sqrt, dot, cross, \
                    randn, zeros, matrix, ones, floor, nonzero, \
                    degrees, arccos, arctan2
from collections import Counter

def read_multi_model_pdb(pdbin, remark_parser=None):
    if type(pdbin) == file:
        source = pdbin
    elif type(pdbin) == str:
        try:
            source = open(pdbin)
        except IOError:
            return None
    else:
        return None
    cell = None
    models = []
    if remark_parser is not None:
        remarks = [remark_parser()]
    for line in source:
        if line[:5] == 'MODEL':
            atoms, anisous = [], []
        elif line[:6] == 'ENDMDL':
            models.append(pdbmolecule(atoms=atoms, cell=cell, anisous=anisous))
            if remark_parser is not None:
                remarks.append(remark_parser())
        elif line[:6] == 'CRYST1':
            cell = pdbcell(line)
        elif line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
            atoms.append(pdbatom(line))
        elif line[:6] == 'ANISOU':
            anisou = pdbanisou(line)
            if anisou == atoms[-1]:
                atoms[-1].SetUij(anisou.GetUij())
            else:
                anisous.append(anisou)
        elif line[:6] == 'REMARK':
            if remark_parser is not None:
                remarks[-1].parse_line(line)
    source.close()
    if remark_parser is not None:
        return models, remarks
    else:
        return models

def ReadPDBfile(pdbin, readcell=True):
    atoms, anisous = [], []
    modelN = 0
    if type(pdbin) == str:
        try:
            source = open(pdbin)
        except IOError:
            return None
    else:
        source = pdbin
    cell = None
    for line in [x.decode() if type(x) is not str else x for x in source ]:
        if line[:5] == 'MODEL':
            modelN += 1
        if line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
            atoms.append(pdbatom(line))
#            xyz.append(atoms[-1].GetR())
        elif line[:6] == 'CRYST1':
            if readcell:
                cell = pdbcell(line)
        elif line[:6] == 'ANISOU':
            anisou = pdbanisou(line)
            if anisou == atoms[-1]:
                atoms[-1].SetUij(anisou.GetUij())
            else:
                anisous.append(anisou)
    if type(pdbin) == str:
        source.close()
    mol = pdbmolecule(atoms=atoms, cell=cell, anisous=anisous)
    if modelN > 0:
        sys.stderr.write('Warning: Input of ReadPDBfile appears to be multi-model (N=%d).  Use read_multi_model_pdb instead.\n' % modelN)
        mol.modelN = modelN
    else:
        mol.modelN = 1
    return mol

def ReadPDBremarks(pdbin):
    if type(pdbin) == file:
        source = pdbin
    elif type(pdbin) == str:
        try:
            source = open(pdbin)
        except IOError:
            return None
    else:
        return None
    remarks = RemarkParser(source)
    if type(pdbin) == str:
        source.close()
    return remarks

def ReadPDBCell(pdbin):
    if type(pdbin) == file:
        source = pdbin
    elif type(pdbin) == str:
        try:
            source = open(pdbin)
        except IOError:
            return None
    else:
        return None
    for line in source:
        if line[:6] == 'CRYST1':
            if type(pdbin) == str:
                source.close()
            return pdbcell(line)
    return None

def ReadPDBTLS(pdbin):
    if type(pdbin) == file:
        source = pdbin
    elif type(pdbin) == str:
        try:
            source = open(pdbin)
        except IOError:
            return None
    else:
        return None
    tls = TLSparser(source)
    if type(pdbin) == str:
        source.close()
    return tls

def WritePDBTLS(tls, tlspath):
    if tls.GetGroupNumber():
        ftls = open(tlspath,'w')
    for (i,group) in enumerate(tls.GetGroups()):
        ftls.write('TLS '+str(i)+'\n')
        for rng in group.GetRanges():
            ftls.write('RANGE '+rng+'\n')
    if tls.GetGroupNumber():
        ftls.close()

def ReadPDBNCS(pdbin):
    if type(pdbin) == file:
        source = pdbin
    elif type(pdbin) == str:
        try:
            source = open(pdbin)
        except IOError:
            return None
    else:
        return None
    ncs = NCSparser(source)
    if type(pdbin) == str:
        source.close()
    return ncs

def WriteNCScode(ncs, ncspath):
    ''' Deduces the NCS definitions from the ncs object and writes them into 
        the ncspath. Returns True on success and False otherwise (e.g. if
        the ncs definitions are corrupted and the resulting output file
        has zero length; the zero length output file is also deleted). '''
    if ncs.GetGroupNumber():
        fncs = open(ncspath,'w')
        for (i,group) in enumerate(ncs.GetGroups()):
            ncs_line =  group.GetCommand()
            if ncs_line:
                fncs.write(ncs_line + '\n')
        ncs_bytes = fncs.tell()
        fncs.close()
        if ncs_bytes:
            return True
        else:
            os.remove(ncspath)
    return False

def WriteNMRstyle(models, name, header=None):
    fout = open(name, 'w')
    tyhe = type(header)
    if tyhe is str:
        fout.write(header)
    elif tyhe is list or tyhe is tuple:
        for line in header:
            fout.write(line)
    for (ix,model) in enumerate(models):
        fout.write('MODEL     %4d\n' % (ix+1))
        for atom in model.GetAtoms():
            fout.write(atom.GetAtomRecord())
        fout.write('ENDMDL\n')
    fout.close()

def cell_and_center(molecule, scale=1.0, cushion=0.0):
    ''' The method returns the CRYST1 line for the P1 cell that could hold
        the entire molecule and the copy of the molecule with shifted
        coordinates that place it in the center. '''
    r = molecule.GetCoordinateArray()
    gabarit = r.ptp(0) + 2*cushion
    abc = scale * gabarit
    retmol = molecule.copy()
    retmol.shift(0.5 * scale * gabarit - r.mean(0))
    celline =  'CRYST1'
    celline += '%9.3f%9.3f%9.3f' % tuple(abc)
    celline += '  90.00  90.00  90.00 P 1                     \n'
    return (celline, retmol)

def getatomid(atom):
    return atom.atomid()

class pdbcell:

    def __init__(self, line):
        self.cryst=line
        self.sg_HM = line[55:66].strip()
        self.cell=line.split()[1:7]
        self.sg = SpaceGroups.xHMN.get(self.sg_HM)
        ang = radians(array([float(self.cell[3]), float(self.cell[4]), float(self.cell[5])]))
        a,b,c = float(self.cell[0]),float(self.cell[1]),float(self.cell[2])
        cang = cos(ang)
        sing = sin(ang[-1])
        v = math.sqrt(1-(cang**2).sum()+(cang**2).prod())
        self.Mcf = array([ [1.0/a,                                  0.0,                                    0.0],
                           [-cang[2]/a/sing,                        1.0/b/sing,                             0.0],
                           [(cang[0]*cang[2]-cang[1])/(a*v*sing),   (cang[1]*cang[2]-cang[0])/(b*v*sing),   sing/c/v]]).T
        self.Mfc = array([ [a,          0.0,                                0.0     ],
                           [b*cang[2],  b*sing,                             0.0     ],
                           [c*cang[1],  c*(cang[0]-cang[1]*cang[2])/sing,   c*v/sing]]).T

    def GetLine(self):
        return self.cryst

    def GetSGHM(self):
        return self.sg_HM

    def GetSG(self):
        return self.sg

    def GetMcf(self):
        return self.Mcf

    def GetMfc(self):
        return self.Mfc

    def GetCellParameters(self):
        return self.cell

    def GetA(self):
        return self.cell[0]

    def GetB(self):
        return self.cell[1]

    def GetC(self):
        return self.cell[2]

    def GetAlpha(self):
        return self.cell[3]

    def GetBeta(self):
        return self.cell[4]

    def GetGamma(self):
        return self.cell[5]

class pdbanisou:
    ''' Defines ANISOU record. '''
    def __init__(self, line):
        self.record = pdbrecord(line)
        self.uij = pdbuij(line[28:70])
        
    def __eq__(self, other):
        return (self.resid()==other.resid()) and (self.name()==other.name()) and (self.altLoc()==other.altLoc())

    def __ne__(self, other):
        return not self==other

# --- Elements of the ATOM line

    def charge(self):
        return self.record.charge()

    def element(self):
        return self.record.element()

    def altLoc(self):
        return self.record.altLoc()

    def segid(self):
        return self.record.segid()
   
    def atomid(self):
        return self.record.atomid()

    def serial(self):
        return self.record.serial()

    def name(self):
        return self.record.name()

    def resid(self):
        return self.record.resid()

# ---
    def GetAnisouRecord(self):
        return ('ANISOU' + '%5d ' % self.serial() + self.atomid() + ' ' + self.uij.GetString() + self.segid().rjust(6) + self.element() + self.charge()).rstrip('\n') + '\n'

    def GetUij(self):
        return self.uij

class pdbuij:
    ''' Defines anisotropic thermal parameters, Uij.  Input values must be 
        the tuple in the following order:
            uij = (u11, u22, u33, u12, u13, u23)
        or  the string conforming to the Uij part of the PDB formatted
        ANISOU record, i.e. 6 integers 7 symbols wide each.'''
    def __init__(self, uij):
        if type(uij) == str:
            self.uij = (int(uij[:7]),
                        int(uij[7:14]),
                        int(uij[14:21]),
                        int(uij[21:28]),
                        int(uij[28:35]),
                        int(uij[35:]))
        elif type(uij) == tuple:
            self.uij = uij
        else:
            raise ValueError('Tuple or string expected to define Uij, got "'+str(type(uij))+'" instead') 

    def GetString(self):
        return '%7d%7d%7d%7d%7d%7d' % self.uij

    def GetAnisotropy(self):
        u = self.uij
        Uij = [[u[0],u[3],u[4]],[u[3],u[1],u[5]],[u[4],u[5],u[2]]]
        v = eigh(Uij)[0]
        return min(v)/max(v)

    def GetValues(self):
        return self.uij

    def SetValues(self, uij):
        if type(uij) == str:
            self.uij = (int(uij[:7]),
                        int(uij[7:14]),
                        int(uij[14:21]),
                        int(uij[21:28]),
                        int(uij[28:35]),
                        int(uij[35:]))
        elif type(uij) == tuple:
            self.uij = uij

class pdbrecord:

    def __init__(self, line):
        line = line.strip().ljust(80)+'\n'
        if line[76:78] == '  ':
            if line[12].isdigit():
                self.line = line[:76] + ' ' + line[13] + line[78:]
            else:
                self.line = line[:76] +   line[12:14]  + line[78:]
        else:
            self.line = line

    def iCode(self):
        return self.line[26]

    def occupancy(self):
        return self.line[54:60]

    def element(self):
        return self.line[76:78]

    def charge(self):
        return self.line[78:80]

    def altLoc(self):
        return self.line[16]

    def segid(self):
        return self.line[72:76].strip()

    def atomid(self):
        return self.line[12:27]

    def chainID(self):
        return self.line[21]

    def serial(self):
        return int(self.line[6:11])

    def name(self):
        return self.line[12:16].strip()

    def resName(self):
        return self.line[17:20].strip()

    def rat(self):
        return self.line[17:20]+self.line[12:16].strip()

    def resSeq(self):
        return int(self.line[22:26])

    def resid(self):
        return self.line[21:27]

    def tempFactor(self):
        return float(self.line[60:66])

    def SetOccupancy(self, o):
        self.line = self.line[:54] + '%6.2f' % o + self.line[60:]

    def SetAltLoc(self, value):
        self.line = self.line[:16]+value+self.line[17:]

    def SetSegid(self, value):
        self.line = self.line[:72] + value.ljust(4)[:4] + self.line[76:]

    def SetChain(self, value):
        self.line = self.line[:21]+value[0]+self.line[22:]

    def SetSerial(self, value):
        self.line = self.line[:6] + '%5d' % value + self.line[11:]

    def SetB(self, value):
        self.line = self.line[:60] + '%6.2f' % value + self.line[66:]

    def SetResID(self, resid):
        self.line = self.line[:21]+resid+self.line[27:]

    def set_res_name(self, name):
        self.line = self.line[:17]+name+self.line[20:]

    def set_name(self, name):
        self.line = self.line[:12]+name+self.line[16:]

class pdbatom:

    def __init__(self, line):
        self.record = pdbrecord(line)
        self.xyz = array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        self.uij = None

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        return (self.resid()==other.resid()) and (self.name()==other.name()) and (self.altLoc()==other.altLoc())

    def __ne__(self, other):
        return not self==other

    def alt(self, other):
        return (self.resid()==other.resid()) and (self.name()==other.name())

    def same_alt(self, other):
        return self.altLoc()==' ' or other.altLoc()==' ' or self.altLoc()==other.altLoc()

    def same_residue(self, other):
        return (self.resid()==other.resid())

    def same_chain(self, other):
        return self.chainID() == other.chainID()

    def test(self,name=None,altloc=None,resn=None,resid=None):
        if name != None and self.name() != name:
            return False
        elif altloc != None and self.altLoc() != altloc:
            return False
        elif resn != None and self.resName() != resn:
            return False
        elif resid != None and self.resid() != resid:
            return False
        else:
            return True

    def test_resid(self, resid):
        return self.resid() == resid

# --- Elements of the ATOM line

    def iCode(self):
        return self.record.iCode()

    def occupancy(self):
        return self.record.occupancy()

    def charge(self):
        return self.record.charge()

    def element(self):
        return self.record.element()

    def altLoc(self):
        return self.record.altLoc()

    def segid(self):
        return self.record.segid()

    def atomid(self):
        return self.record.atomid()

    def chainID(self):
        return self.record.chainID()

    def serial(self):
        return self.record.serial()

    def name(self):
        return self.record.name()

    def mass(self):
        return pdbnames.GetMass(self.record.element().strip())

    def resName(self):
        return self.record.resName()

    def resSeq(self):
        return self.record.resSeq()

    def resid(self):
        return self.record.resid()

    def tempFactor(self):
        return self.record.tempFactor()

# ---

    def SetUij(self, uij):
        self.uij = uij

    def GetUij(self):
        return self.uij

    def GetUijValues(self):
        return self.uij.GetValues()

    def SetUijValues(self, uij):
        self.uij.SetValues(uij)

    def prime_uij(self, overwrite=False):
        ''' Initializes the anisotropic ADPs for the atom using its isotropic
            B-factor.  Use the overwrite flag to force the anisotropic ADP
            reset when they are alreay present. '''
        if not self.IsAnisotropic() or overwrite:
            ueq = 126.65148*self.GetB()
            self.uij = pdbuij((ueq,ueq,ueq,0.0,0.0,0.0))

    def GetAnisotropy(self):
        if self.uij:
            return self.uij.GetAnisotropy()
        else:
            return 1.0

    def GetChain(self):
        ''' Returns atom's chain ID. '''
        return self.chainID()

    def SetChain(self, value):
        ''' Sets atom chain ID to value (must be a single character). '''
        self.record.SetChain(value)

    def GetAtomID(self):
        return self.record.atomid()

    def GetResID(self):
        return self.resid()

    def SetResID(self, resid):
        ''' Sets the atom resid (must be in correct format, i.e. "A 156B"). '''
        self.record.SetResID(resid)

    def get_alt_resid(self):
        return (self.resid() + '-' + self.altLoc()).strip().strip('-')

    def GetElement(self):
        return self.element().strip() 

    def GetName(self):
        return self.record.name()

    def GetAtomRecord(self):
        return self.record.line[:30]+'%8.3f%8.3f%8.3f' % tuple(self.xyz)+self.record.line[54:]

    def IsAnisotropic(self):
        return bool(self.uij)

    def GetAnisouRecord(self):
        if self.uij:
            return ('ANISOU' + '%5d ' % self.serial() + self.atomid() + ' ' + self.uij.GetString() + self.segid().rjust(6) + self.element() + self.charge()).rstrip('\n') + '\n'
        else:
            return ''

    def GetAtomTitle(self):
        return self.resName() + ' ' + self.chainID() + ('%d' % self.resSeq() + self.iCode()).ljust(6) + ' ' + self.name()

    def GetAltLoc(self):
        return self.altLoc()

    def SetAltLoc(self, value):
        self.record.SetAltLoc(value)

    def GetSegid(self):
        return self.segid()

    def SetSegid(self,value):
        self.record.SetSegid(value)

    def HasAltConf(self):
        return bool(self.altLoc().strip())

    def HasCNSAltConf(self):
        return bool(re.search('AC[0-9]+',self.segid()))

    def GetCNSAltLoc(self):
        return re.search('AC[0-9]+',self.segid()).group()[-1]

    def get_res_name(self):
        return self.record.resName()

    def set_res_name(self, name):
        self.record.set_res_name(name)

    def rjust_res_name(self):
        self.record.set_res_name(self.record.resName().strip().rjust(3))

    def GetResTitle(self):
        return self.resName() + ' ' + self.chainID() + '%d' % self.resSeq() + self.iCode()

    def GetB(self):
        ''' Return the atomic ADP. '''
        return self.record.tempFactor()

    def SetB(self, B):
        ''' Set the atomic ADP. '''
        self.record.SetB(B)

    def GetOccupancy(self):
        return float(self.occupancy())

    def SetOccupancy(self, o):
        self.record.SetOccupancy(o)

    def IsWater(self):
        return self.resName() == 'HOH'

    def IsPolar(self):
        return pdbnames.IsPolar(self.element().strip())

    def IsProteinBackbone(self):
        ''' True if atom belongs to protein backbone, False otherwise. '''
        return pdbnames.Is3Amino(self.resName()) and pdbnames.MaybeBackbone(self.name())

    def IsProtein(self):
        ''' True if atom belongs to a protein, False otherwise. '''
        return pdbnames.Is3Amino(self.resName())

    def IsHetero(self):
        ''' True if atom is a heteroatom (i.e. not protein or water). '''
        return pdbnames.IsHetero(self.resName())

    def IsMetal(self):
        ''' True if atom is a metal (as defined in pdbnames). '''
        return pdbnames.IsMetal(self.name())

    def IsBackbone(self):
        ''' True if atom belongs to protein/DNA backbone, False otherwise. '''
        return pdbnames.MaybeBackbone(self.name()) and pdbnames.NotWater(self.resName())

    def NotBackbone(self):
        ''' True if atom does not belong protein/DNA backbone, False otherwise. '''
        return pdbnames.NotBackbone(self.name())

    def IsProteinSidechain(self):
        ''' True if atom belongs to protein side chain, False otherwise. '''
        return self.IsProtein() and pdbnames.NotBackbone(self.name())

    def shift(self, xyz):
        self.xyz += xyz

    def transform(self, M):
        self.xyz = array(dot(M,self.xyz.T))[0]

    def GetR(self): #convert to GetXYZarray
        return self.xyz

    def SetR(self, vector): #convert to SetXYZarray
        self.xyz = vector

    def GetXYZarray(self): #convert to GetR
        return self.xyz

    def SetXYZarray(self, xyz): #convert to SetR
        self.xyz = xyz

    def GetSerial(self):
        return self.serial()

    def SetSerial(self, value):
        self.record.SetSerial(value)

    def ApplyOperator(self, op):
        self.xyz = array(op(self.xyz[0],self.xyz[1],self.xyz[2]))

    def copy(self):
        return copy.deepcopy(self)

    def get_vdw_radius(self):
        return pdbnames.GetVDWRadius(self.element().strip())

    def set_name(self, name):
        self.record.set_name(name)

    def rat(self):
        return self.record.rat()

    def dx(self, other):
        return other.xyz[0]-self.xyz[0]
    def dy(self, other):
        return other.xyz[1]-self.xyz[1]
    def dz(self, other):
        return other.xyz[2]-self.xyz[2]
    def dr(self, other):
        return other.xyz-self.xyz
    def distance(self, other):
        return sqrt((self.dr(other)**2).sum())
    def R2(self, other):
        return (self.dr(other)**2).sum()
    def noise(self, xnoise=0.1, bnoise=0.1, occnoise=0.0):
        self.xyz += xnoise*randn(3)
        self.SetB(math.fabs(self.GetB()*random.gauss(1,bnoise)))
        if self.IsAnisotropic():
            self.SetUijValues(tuple((array(self.GetUijValues())*abs(1.0+bnoise*randn(6))).astype(int)))
        atomocc = float(atom.GetOccupancy())
        if atomocc<1.0 and atomocc>0.0:
            newocc = random.gauss(atomocc,occnoise)
            while newocc>1.0 or newocc<0.0:
                newocc = random.gauss(atomocc,occnoise)
            self.SetOccupancy(newocc)

class single_pdbatom(pdbatom):

    def __init__(self, line):
        self.record = pdbrecord(line)
        self.xyz = array([float(line[30:38]), float(line[38:46]), float(line[46:54])],dtype=float)
        self.uij = None

class pdbmolecule:

    def __init__(self, code=None, atoms=None, cell=None, anisous=None):
        self.modelN = 0
        if atoms:
            self.atoms = atoms
        else:
            self.atoms = []
        self.cell = cell
        if code != None:
            self.pdbid = code.lower()
            self.pdbname = self.pdbid+'.pdb'
            if not os.access('pdb-download',os.R_OK):
                os.mkdir('pdb-download')
            if not os.access('pdb-download/'+self.pdbname,os.R_OK):
                try:
                    pdbfile=urllib.request.urlretrieve('http://www.rcsb.org/pdb/files/'+self.pdbid+'.pdb.gz')
                    fout = open('pdb-download/'+self.pdbname, 'w')
                    fin = gzip.open(pdbfile[0])
                    for line in fin:
                        fout.write(line)
                    fin.close()
                    fout.close()
                except IOError:
                    fin.close()
                    fout.close()
                    sys.stderr.write("PDBMOLECULE: I/O error.  File "+self.pdbname+" can't be opened or retrieved from PDB\n")
            fin = open('pdb-download/'+self.pdbname)
            for line in fin:
                if line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
                    self.atoms.append(pdbatom(line))
                if line[:5] == 'MODEL':
                    self.modelN += 1
            fin.close()
        if anisous:
            for anisou in anisous:
                i = self.find(anisou)
                if i>=0:
                    self.atoms[i].SetUij(anisou.GetUij())
        if self.modelN == 0:
            self.modelN = 1
        self.cartesian = True

    def __len__(self):
        return self.GetAtomNumber()

    def __iter__(self):
        return self.atoms.__iter__()

    def backbone(self):
        return backbone(self)

    def is_multi_model(self):
        return self.modelN > 1

    def rjust_res_names(self):
        for atom in self.atoms:
            atom.rjust_res_name()

    def prime_uij(self, overwrite=False, what='all', listik=False, *args, **kwargs):
        ''' Initializes anisotropic ADPs.  The overwrite flag defines if
            the ANISOU record will be overwritten if already persent.
            Selection syntax the same as in atom_lister method.'''
        for atom in self.atom_getter(what, listik, *args, **kwargs):
            atom.prime_uij(overwrite=overwrite)

    def GetSpaceGroup(self):
        if self.cell:
            return self.cell.GetSG()

    def SerialReset(self):
        ''' Reset serial numbers of atoms.  Use this after atom insertion/removal/rearrangement.'''
        for (i,atom) in enumerate(self.atoms):
            atom.SetSerial(i+1)

    def BfactorReset(self, b=20.0):
        ''' Reset B-factors of all atoms to the supplied value. '''
        for atom in self.atoms:
            atom.SetB(b)

    def SetBfactorValues(self, bvalues, what='all', listik=False, *args, **kwargs):
        ''' Sets the atomic B-factors using the provided values mapped 
            into the list of atoms. If bvalues is a single number, all 
            the atoms in the list will have the same B-factor.  
            Selection syntax the same as in atom_lister method. '''
        if type(bvalues) == float:
            for atom in self.atom_getter(what, listik, *args, **kwargs):
                atom.SetB(bvalues)
        else:
            atoms = self.atom_getter(what, listik, *args, **kwargs)
            assert len(bvalues)==len(atoms), 'Shape mismatch between selection and Bvalues vector.'
            for (i, atom) in enumerate(atoms):
                atom.SetB(bvalues[i])

    def set_occupancies(self, values, what='all', listik=False, *args, **kwargs):
        ''' Sets the atomic occupancies using the provided values mapped
            into  the list of atoms.  If values is a single number, all 
            the atoms in the list will have the same occupancy.  
            Selection syntax the same as in atom_lister method. '''
        if type(values) == float:
            for atom in self.atom_getter(what, listik, *args, **kwargs):
                atom.SetOccupancy(values)
        else:
            atoms = self.atom_getter(what, listik, *args, **kwargs)
            assert len(values)==len(atoms), 'Shape mismatch between selection and values vector.'
            for (i, atom) in enumerate(atoms):
                atom.SetOccupancy(values[i])

    def PopAtom(self, atomi):
        ''' Return the atom and delete it from the molecule.  Use with caution. '''
        return self.atoms.pop(atomi)

    def DeleteAtoms(self, listik):
        ''' Delete atoms from the molecule.  Use with caution. Atoms
            with indices in listik are removed. '''
        for i in sorted(listik, reverse=True):
            del self.atoms[i]

    def InsertAtom(self, atoms, atomi=0):
        ''' Insert atoms before position atomi in the atom list.  Notice that an atom
            cannot be inserted at the end, use AppendAtom() for that.  By default,
            atom is inserted at the beginning of the atom list.  Notice that syntax
            provides different order of parameters as compared to insert method of
            a regular python list.  Insert either a single atom object or list of atoms.'''
        if type(atoms) == list:
            self.atoms = self.atoms[:atomi]+atoms+self.atoms[atomi:]
        else:
            self.atoms.insert(atomi, atoms)

    def AppendAtom(self, atoms):
        ''' Append atoms at the end of the atom list. '''
        if type(atoms) == list:
            self.atoms += atoms
        else:
            self.atoms.append(atoms)

    def AppendMolecule(self, other):
        ''' Append atoms from another molecule. No sanity checks. '''
        self.AppendAtom(other.GetAtoms())

    def matchAlts(self, other):
        listik1 = self.resid_lister(listik=self.merge_listers(['altconf','protein']))
        listik2 = other.resid_lister(listik=other.merge_listers(['altconf','protein']))
        listik = list(set(listik1).union(listik2))
        residues1 = self.get_residues('resids', resids=listik)
        residues2 = other.get_residues('resids', resids=listik)
        listik = list(set(residues1.keys()).intersection(residues2.keys()))
        for resid in listik:
            n1, n2 = residues1[resid].GetAltNum(), residues2[resid].GetAltNum()
            if n1<2:
                if n2==0:
                    residues1[resid].AssignAlts(' ')
                else:
                    residues1[resid].renameSingle2Alt(residues2[resid])
            else:
                if n1==n2:
                    residues1[resid].renameAlt2Alt(residues2[resid])
                elif n2<2:
                    residues1[resid].renameAlt2Single(residues2[resid])
                else:
                    sys.stdout.write(resid+": different number of alt confs not yet supported\n")

    def GetChangedAlts(self, other, fProteinOnly=True):
        ''' Returns the tuple containing two lists.  First lists the residues that 
            have alternate conformers that are missing in the other model.  The second
            lists the single conformer residues that expand into multiple conformations
            in the other model. '''
        listik1 = self.ListAltResidues(fProteinOnly=fProteinOnly)
        listik2 = other.ListAltResidues(fProteinOnly=fProteinOnly)
        return (list(set(listik1).difference(listik2)), list(set(listik2).difference(listik1)))

    def acCNS2PDB(self):
        ''' Looks for CNS-styled alternate conformers and converts them to PDB format. '''
        altlist = []
        for (i,atom) in enumerate(self.atoms):
            if atom.HasCNSAltConf():
                altlist.append(i)
                atom.SetAltLoc(atom.GetCNSAltLoc())
        while altlist:
            lead = altlist.pop(0)
            acgroup = self.findalt(self.GetAtom(lead),altlist)
            for (i,atomi) in enumerate(acgroup):
                self.InsertAtom(self.PopAtom(atomi),lead+1+i)
                altlist.remove(atomi)
                for (j,atomj) in enumerate(altlist):
                    if atomj < atomi:
                        altlist[j] += 1
        self.SerialReset()

    def acPDB2CNS(self):
        ''' Converts alternate conformers into CNS style. '''
        altlist = self.ListAltConf()
        while altlist:
            lead = altlist.pop(0)
            self.GetAtom(lead).SetSegid('AC1 ')
            acgroup = self.findalt(self.GetAtom(lead),altlist)
            for (i,atomi) in enumerate(acgroup):
                self.AppendAtom(self.PopAtom(atomi))
                self.GetAtom(-1).SetSegid(('AC'+str(i+2)).ljust(4))
                altlist.remove(atomi)
                for j in range(len(altlist)):
                    altlist[j] -= 1
        self.SerialReset()

    def find(self, atom2find, subset=None):
        ''' Finds the atom in the molecule matching atom2find and returns its index.  If no such atom exists, returns None.
            Atoms are matched by residue ID (chain ID + residue number), atom name and alternate conformer ID. 
            Searching within part of the molecule will speed things up and can be accomplished by passing subset 
            parameter, which is a list of atom indices to include in the search.  Such lists can be generated by
            some class methods, including atom_lister().
        '''
        if subset is None:
            for (i,atom) in enumerate(self.atoms):
                if atom == atom2find:
                    return i
        else:
            for i in subset:
                if self.atoms[i] == atom2find:
                    return i
        return None

    def find_atom(self, atom2find, listik=False):
        ''' Finds the atom in the molecule matching atom2find and returns the 
            atom object.  If no such atom exists, returns None. Atoms are matched 
            by residue ID (chain ID + residue number), atom name and alternate 
            conformer ID. Searching within part of the molecule will speed 
            things up and can be accomplished by passing subset parameter, 
            which is a list of atom indices to include in the search.  Such 
            lists can be generated by atom_lister() method. 
        '''
        atoms = self.GetListedAtoms(self.__ensure_listik_(listik))
        if atoms.count(atom2find):
            return atoms[atoms.index(atom2find)]
        else:
            return None

    def findalt(self, atom2find, subset=None):
        ''' Similar to find() method, but returns the list of alternate conformer atoms. '''
        altatoms = []
        if not subset:
            for (i,atom) in enumerate(self.atoms):
                if atom.alt(atom2find):
                    altatoms.append(i)
        else:
            for i in subset:
                if self.atoms[i].alt(atom2find):
                    altatoms.append(i)
        return altatoms

    def shift(self, T, listik=False):
        listik = self.__ensure_listik_(listik)
        T = array(T)
        for atom in self.GetListedAtoms(listik):
            atom.shift(T)

    def transform(self, U, listik=False):
        listik = self.__ensure_listik_(listik)
        for atom in self.GetListedAtoms(listik):
            atom.transform(matrix(U))

    def ExtractSelection(self, listik):
        return pdbmolecule(atoms=self.GetListedAtoms(listik), cell=self.cell)

    def Extract(self, what, listik=False, *args, **kwargs):
        '''
            Make a copy of the molecule selecting only atoms of the specified
            type.  Atom type is selected by what argument which allows the
            following options.  It is a case-insensitive string.
            'all'               - all atoms
            'backbone'          - backbone atoms
            'sidechain'         - side chain atoms
            'water'             - water molecules
            'element_X'         - chemical element (X = H,C,O,N,P,CA, etc)
            'name_X'            - atom name (X = atom name)
            'protein'           - protein atoms
            'p_bb' or 
            'proteinbackbone'   - protein backbone atoms     
            'resid'             - specific residue ID, pass the value 
                                    as extra resid parameter      
        '''
        return pdbmolecule(atoms=self.atom_getter(what, listik, *args, **kwargs), cell=self.cell)

    def get_elements(self, what='all', listik=False, *args, **kwargs):
        return sorted(set([a.element() for a in self.atom_getter(what, listik, *args, **kwargs)]))

    def get_atom_types(self, what='all', listik=False, *args, **kwargs):
        return sorted(set([a.name() for a in self.atom_getter(what, listik, *args, **kwargs)]))

# --- "GetAtomSomething" methods ---

    def GetAtomNumber(self):
        ''' Returns the total number of atoms in the molecule. '''
        return len(self.atoms)

    def GetAtomNumberByType(self):
        ''' Returns the tuple with the number of atoms in the following categories:
            all, protein, protein backbone, protein side chains, waters,
            heteroatoms. '''
        Nall, Nprot, Nbb, Nsc, Nwat, Nhetero = len(self.atoms), 0, 0, 0, 0, 0
        for atom in self.atoms:
            if atom.IsProtein():
                Nprot += 1
                if atom.IsProteinBackbone():
                    Nbb += 1
                else:
                    Nsc += 1
            elif atom.IsWater():
                Nwat += 1
            else:
                Nhetero += 1
        return (Nall, Nprot, Nbb, Nsc, Nwat, Nhetero)

    def GetAtom(self,i):
        ''' Returns the ith atom object. '''
        return self.atoms[i]

    def GetAtoms(self):
        ''' Returns the list of atom objects comprising the molecule. '''
        return self.atoms

    def GetListedAtoms(self, listik):
        ''' Returns the list of atoms from selection. '''
        return array(self.atoms)[listik].tolist()

    def GetAtomResID(self, i):
        ''' Returns the resid of the residue atom belongs to (e.g. "A 234B" or 
            "L  16 "). '''
        return self.atoms[i].GetResID()

    def get_atom_alt_resid(self, i):
        ''' Returns the resid of the residue atom belongs to, including the
            alternate conformer identifier, if non-empty (e.g. B47-A, notice 
            the dash separator). '''
        return self.atoms[i].get_alt_resid()

    def GetAtomResidueName(self, i):
        ''' Returns the 3-letter code of the residue name. '''
        return self.GetAtom(i).get_res_name()

    def GetAtomTitle(self, i):
        ''' Returns the formatted title of the atom (e.g. "ALA L27B C", "TYR L32  CG"). '''
        return self.atoms[i].GetAtomTitle()

    def GetAtomXYZarray(self,i): #convert to GetR
        ''' Returns the vector of atom coordinates as [x,y,z]. '''
        return self.atoms[i].GetXYZarray()

    def GetAtomR(self, i):
        ''' Returns the vector of atom coordinates as array([x,y,z]). '''
        return self.atoms[i].GetR()

    def GetAtomName(self, i):
        return self.atoms[i].name()

    def GetAtomB(self, i):
        ''' Returns the atom temperature factor. '''
        return self.atoms[i].GetB()

    def GetAtomAnisotropy(self, i):
        ''' Returns the ratio of longest to shortest thermal ellipsoid axis.  
            If there is no ANISOU record associated with the atom, defaults to 1.0. '''
        return self.atoms[i].GetAnisotropy()

    def get_atom_vdw_radius(self, i):
        return self.atoms[i].get_vdw_radius()

# --- Get various  vectors

    def GetBvector(self, what='all', listik=False, *args, **kwargs):
        ''' Returns an array of B-factor values for a list of atoms.
            Defaults to all atoms in the molecule. '''
        return array([a.tempFactor() for a in self.atom_getter(what, listik, *args, **kwargs)])

    def GetResidueBvector(self, mode='lin'):
        bs = []
        residues = self.get_residues()
        resids = sorted(residues.keys())
        bs = [residues[resid].GetGroupedBs(mode) for resid in resids]
        return (resids, bs)

    def GetResidueBvectorByChain(self):
        b0, b1, b2, resids = {}, {}, {}, {}
        residues = self.get_residues()
        for resid in sorted(residues.keys()):
            chid = resid[0]
            ball, bmain, bside = residues[resid].GetGroupedBs()
            if chid in resids:
                resids[chid].append(resid[1:])
                b0[chid].append(ball)
                b1[chid].append(bmain)
                b2[chid].append(bside)
            else:
                resids[chid] = [resid[1:]]
                b0[chid] = [ball]
                b1[chid] = [bmain]
                b2[chid] = [bside]
        return (resids, b0, b1, b2)

    def GetOccupancyVector(self, listik=False):
        ''' Returns an array of occupancy values for a list of atoms.
            Defaults to all atoms in the molecule. '''
        return array([a.GetOccupancy() for a in self.atom_getter(listik=listik)])

    def GetMassVector(self, listik=False):
        ''' Returns an array of mass values for a list of atoms.
            Defaults to all atoms in the molecule. '''
        return array([a.mass() for a in self.atom_getter(listik=listik)])

    def GetMOVector(self, what='all', listik=False, *args, **kwargs):
        ''' Returns an array of occupancy weighted masses or a list of 
            atoms. Defaults to all atoms in the molecule. '''
        return array([a.GetOccupancy()*a.mass() for a in self.atom_getter(what, listik, *args, **kwargs)])

# --- Methods that return atom lists


    def atom_getter(self, what='all', listik=False, *args, **kwargs):
        ''' Returns the list of atoms based on the string defining the 
            group. Recognized values are (any mix of lower and upper case 
            characters):
            'all'               - all atoms
            'backbone'          - backbone atoms
            'sidechain'         - side chain atoms
            'water'             - water molecules
            'notwater'          - other than water molecules
            'polar'             - polar atoms
            'hetero'            - heteroatoms
            'metal'             - metals
            'element_X'         - chemical element (X = H,C,O,N,P,CA, etc)
            'name_X'            - atom name (X = atom name)
            'names'             - atom from a list of names , pass the
                                    list of names as names parameter
            'protein'           - protein atoms
            'p_bb' or 
            'proteinbackbone'   - protein backbone atoms
            'p_sc' or
            'proteinsidechain'  - protein sidec chain atoms
            'resid'             - specific residue ID, pass the value 
                                    as extra resid parameter
            'resids'            - atoms from a set of residues, pass
                                    the list of residue IDs as resids
                                    parameter
            'chid'              - specific chain ID, pass the value 
                                    as extra chid parameter    
            'chids'             - atoms from a set of chains, pass
                                    the list of chain IDs as chids
                                    parameter
            'resnames'          - atoms from a set of residue types, pass
                                    the list of residue types as resnames
                                    parameter
            'vicinity'          - atoms in the vicinity of a group of 
                                    "core atoms".  Must pass the list of core
                                    atoms as coreatoms argument and optional
                                    rcutoff (defaults to 4A).  Notice that the
                                    atom list may become unsorted
            'sphere'            - similar to vicinity, but includes core
                                    atoms as well
            'rat'               - specific residue/atom name combination.
                                    For example, ASPOD1 designates OD1
                                    atoms in aspartic acids.  Must provide 
                                    the list of combos to include as rats 
                                    parameter
            'altconf'           - atoms that have alternate conformer
                                    identifier.  If acvalue parameter is
                                    provided, only atoms that match it
                                    are selected (this could be an
                                    inclusive list/string).
            'altgroup'          - atoms that match the atom specified 
                                    in atom parameter.  Will include the
                                    seed atom itself.
        '''
        whatlow = what.lower()
        atoms = self.GetListedAtoms(self.__ensure_listik_(listik))
        if whatlow == 'all':
            return atoms
        elif whatlow == 'backbone':
            return [x for x in atoms if x.IsBackbone()]
        elif whatlow == 'sidechain':
            return [x for x in atoms if x.NotBackbone()]
        elif whatlow == 'water':
            return [x for x in atoms if x.IsWater()]
        elif whatlow == 'notwater':
            return [x for x in atoms if not x.IsWater()]
        elif whatlow == 'polar':
            return [x for x in atoms if x.IsPolar()]
        elif whatlow == 'hetero':
            return [x for x in atoms if x.IsHetero()]
        elif whatlow == 'metal':
            return [x for x in atoms if x.IsMetal()]
        elif whatlow[:8] == 'element_':
            return [x for x in atoms if x.GetElement() == what[8:]]
        elif whatlow[:5] == 'name_':
            return [x for x in atoms if x.GetName() == what[5:]]
        elif whatlow == 'names':
            return [x for x in atoms if x.GetName() in kwargs['names']]
        elif whatlow == 'protein':
            return [x for x in atoms if x.IsProtein()]
        elif whatlow == 'proteinbackbone' or whatlow == 'p_bb':
            return [x for x in atoms if x.IsProteinBackbone()]
        elif whatlow == 'proteinsidechain' or whatlow == 'p_sc':
            return [x for x in atoms if x.IsProteinSidechain()]
        elif whatlow == 'resid':
            return [x for x in atoms if x.GetResID() == kwargs['resid']]
        elif whatlow == 'resids':
            return [x for x in atoms if x.GetResID() in kwargs['resids']]
        elif whatlow == 'chid':
            return [x for x in atoms if x.chainID() == kwargs['chid']]
        elif whatlow == 'chids':
            if kwargs.get('chids') is not None:
                return [x for x in atoms if x.chainID() in kwargs['chids']]
            return atoms
        elif whatlow == 'resnames':
            return [x for x in atoms if x.get_res_name() in kwargs['resnames']]
        elif whatlow == 'vicinity':
            r2cutoff = kwargs.get('rcutoff', 4.0)**2
            coreXYZ = array([x.GetR() for x in kwargs['coreatoms']])
            return sorted([x for x in set(atoms).difference(kwargs['coreatoms']) if (((coreXYZ-x.GetR())**2).sum(1)<=r2cutoff).any()], key=getatomid)
        elif whatlow == 'sphere':
            r2cutoff = kwargs.get('rcutoff', 4.0)**2
            coreXYZ = array([x.GetR() for x in kwargs['coreatoms']])
            return [x for x in set(atoms).difference(kwargs['coreatoms']) if (((coreXYZ-x.GetR())**2).sum(1)<=r2cutoff).any()]+kwargs['coreatoms']
        elif whatlow == 'rat':
            return [x for x in atoms if x.rat() in kwargs['rats']]
        elif whatlow == 'altconf':
            if kwargs.get('acvalue') is None:
                return [x for x in atoms if x.HasAltConf()]
            else:
                return [x for x in atoms if x.altLoc() in kwargs['acvalue']]
        elif whatlow == 'altgroup':
            return [x for x in atoms if kwargs['atom'].alt(x)]

    def atom_lister(self, what='all', listik=False, *args, **kwargs):
        ''' Returns the list of atom indices based on the string defining the 
            group. Recognized values are (any mix of lower and upper case 
            characters):
            'all'               - all atoms
            'backbone'          - backbone atoms
            'sidechain'         - side chain atoms
            'water'             - water molecules
            'notwater'          - other than water molecules
            'polar'             - polar atoms
            'hetero'            - heteroatoms
            'metal'             - metals
            'element_X'         - chemical element (X = H,C,O,N,P,CA, etc)
            'name_X'            - atom name (X = atom name)
            'names'             - atom from a list of names, pass the
                                    list of names as names parameter
            'protein'           - protein atoms
            'p_bb' or 
            'proteinbackbone'   - protein backbone atoms   
            'p_sc' or
            'proteinsidechain'  - protein sidec chain atoms
            'resid'             - specific residue ID, pass the value 
                                    as extra resid parameter    
            'resids'            - atoms from a set of residues, pass
                                    the list of residue IDs as resids
                                    parameter
            'chid'              - specific chain ID, pass the value 
                                    as extra chid parameter    
            'chids'             - atoms from a set of chains, pass
                                    the list of chain IDs as chids
                                    parameter
            'resnames'          - atoms from a set of residue types, pass
                                    the list of residue types as resnames
                                    parameter
            'vicinity'          - atoms in the vicinity of a group of 
                                    "core atoms".  Must pass the list of core
                                    atoms as corelist argument and optional
                                    rcutoff (defaults to 4A)
            'sphere'            - similar to vicinity, but includes core
                                    atoms as well
            'rat'               - specific residue/atom name combination.
                                    For example, ASPOD1 designates OD1
                                    atoms in aspartic acids.  Must provide 
                                    the list of combos to include as rats 
                                    parameter
            'altconf'           - atoms that have alternate conformer
                                    identifier.  If acvalue parameter is
                                    provided, only atoms that match it
                                    are selected(this could be an
                                    inclusive list/string).
            'altgroup'          - atoms that match the atom specified 
                                    in atomid parameter.  Will include 
                                    the seed atom itself.
        '''
        whatlow = what.lower()
        listik = self.__ensure_listik_(listik)
        if whatlow == 'all':
            return listik
        elif whatlow == 'backbone':
            return [i for i in listik if self.atoms[i].IsBackbone()]
        elif whatlow == 'sidechain':
            return [i for i in listik if self.atoms[i].NotBackbone()]
        elif whatlow == 'water':
            return [i for i in listik if self.atoms[i].IsWater()]
        elif whatlow == 'notwater':
            return [i for i in listik if not self.atoms[i].IsWater()]
        elif whatlow == 'polar':
            return [i for i in listik if self.atoms[i].IsPolar()]
        elif whatlow == 'hetero':
            return [i for i in listik if self.atoms[i].IsHetero()]
        elif whatlow == 'metal':
            return [i for i in listik if self.atoms[i].IsMetal()]
        elif whatlow[:8] == 'element_':
            return [i for i in listik if self.atoms[i].GetElement() == what[8:]]
        elif whatlow[:5] == 'name_':
            return [i for i in listik if self.atoms[i].GetName() == what[5:]]
        elif whatlow == 'names':
            return [i for i in listik if self.atoms[i].GetName() in kwargs['names']]
        elif whatlow == 'protein':
            return [i for i in listik if self.atoms[i].IsProtein()]
        elif whatlow == 'proteinbackbone' or whatlow == 'p_bb':
            return [i for i in listik if self.atoms[i].IsProteinBackbone()]
        elif whatlow == 'proteinsidechain' or whatlow == 'p_sc':
            return [i for i in listik if self.atoms[i].IsProteinSidechain()]
        elif whatlow == 'resid':
            return [i for i in listik if self.atoms[i].GetResID() == kwargs['resid']]
        elif whatlow == 'resnames':
            return [i for i in listik if self.atoms[i].get_res_name() in kwargs['resnames']]
        elif whatlow == 'resids':
            return [i for i in listik if self.atoms[i].GetResID() in kwargs['resids']]
        elif whatlow == 'chid':
            return [i for i in listik if self.atoms[i].chainID() == kwargs['chid']]
        elif whatlow == 'chids':
            if kwargs.get('chids') is not None:
                return [i for i in listik if self.atoms[i].chainID() in kwargs['chids']]
            return listik
        elif whatlow == 'vicinity':
            r2cutoff = kwargs.get('rcutoff', 4.0)**2
            coreXYZ = self.GetCoordinateArray(listik=kwargs['corelist'])
            return sorted([i for i in set(listik).difference(kwargs['corelist']) if (((coreXYZ-self.GetAtomR(i))**2).sum(1)<=r2cutoff).any()])
        elif whatlow == 'sphere':
            r2cutoff = kwargs.get('rcutoff', 4.0)**2
            coreXYZ = self.GetCoordinateArray(listik=kwargs['corelist'])
            return sorted([i for i in set(listik).difference(kwargs['corelist']) if (((coreXYZ-self.GetAtomR(i))**2).sum(1)<=r2cutoff).any()]+kwargs['corelist'])
        elif whatlow == 'rat':
            return [i for i in listik if self.atoms[i].rat() in kwargs['rats']]
        elif whatlow == 'altconf':
            if kwargs.get('acvalue') is None:
                return [i for i in listik if self.atoms[i].HasAltConf()]
            else:
                return [i for i in listik if self.atoms[i].altLoc() in kwargs['acvalue']]
        elif whatlow == 'altgroup':
            return [i for i in listik if self.atoms[i].alt(kwargs['atomid'])]

    def merge_listers(self, whats, listik=False, *args, **kwargs):
        '''
            Merge multiple atom listers.
        '''
        for what in whats:
            listik = self.atom_lister(what, listik, *args, **kwargs)
        return listik

    def merge_getters(self, whats, listik=False, *args, **kwargs):
        '''
            Merge multiple atom getters.
        '''
        return [self.atoms[i] for i in self.merge_listers(whats, listik, *args, **kwargs)]

    def resid_lister(self, what='all', listik=False, *args, **kwargs):
        ''' Returs the list of residue IDs from an atom selection.
            Selection parameters are the same as in atom_getter method.
        '''
        return list(set([a.resid() for a in self.atom_getter(what, listik, *args, **kwargs)]))

    def list_bcutoff(self, value, listik=False):
        ''' Returns the list of indices for atoms with B-factors less or 
            equal to the  cutoff value provided. '''
        return [i for i in self.__ensure_listik_(listik) if self.atoms[i].GetB()<=value]

    def list_high_occupancy(self, value, listik=False):
        ''' Returns the list of indices for atoms with occupancy above or
            equal to the cutoff value provided. '''
        return [i for i in self.__ensure_listik_(listik) if self.atoms[i].GetOccupancy()>=value]

    def list_low_occupancy(self, value, listik=False):
        ''' Returns the list of indices for atoms with occupancy below the
            cutoff value provided. '''
        return [i for i in self.__ensure_listik_(listik) if self.atoms[i].GetOccupancy()<value]

    def __enumerate_atoms_(self, listik=False):
        if not listik:
            return enumerate(self.atoms)
        else:
            return [(i,self.GetAtom(i)) for i in listik]

    def __ensure_listik_(self, listik=False):
        if listik is False:
            return range(len(self.atoms))
        else:
            return listik

    def __ensure_atoms_(self, listik=False):
        if listik is False:
            return self.atoms
        else:
            return [self.atoms[i] for i in listik]

    def get_residnames(self, what='all', listik=False, *args, **kwargs):
        ''' Returns a distionary with residue IDs as keys and residue
            names as values. Selection parameters are the same as in 
            atom_getter method. '''
        return dict(set([(a.resid(),a.resName()) for a in self.atom_getter(what, listik, *args, **kwargs)]))
        

    def ListByAltConf(self, ac, listik=False):
        '''
        Return atom indices (from the listik if provided) with the selected
        alternate conformer code.
        '''
        return [i for i,atom in self.__enumerate_atoms_(listik) if atom.GetAltLoc() == ac]

    def ListAltConf(self, listik=False):
        '''
        Return atom indices (from the listik if provided) that are
        labeled as alternate conformers.
        '''
        return [i for i,atom in self.__enumerate_atoms_(listik) if atom.HasAltConf()]

    def ListChainSplit(self, listik=False):
        ''' Returns the dictionary of list of atoms from individual chains.'''
        chainsplit = {}
        for (i,atom) in self.__enumerate_atoms_(listik):
            chid = atom.GetChain()
            try:
                chainsplit[chid].append(i)
            except KeyError:
                chainsplit[chid] = [i]
        return chainsplit
        
    def ListCompleteResidues(self, listik):
        return self.atom_lister('resids', resids=self.resid_lister(listik=listik))

    def GetProteinResidueNumber(self):
        return len(self.resid_lister('protein'))

    def ListAltResidues(self, fProteinOnly=False):
        ''' Returns the list of residues that contain alternate conformers.
            This is slower but more robust version of ListResidues(hasAltConf=True),
            as it selects residues via pdbresidue.HasAltConf() method,
            which returns False whenever there is only one altLoc in the residue.
            The ListResidues methods goes through the list of all the atoms and 
            uses the pdbatom.HasAltConf() method, which returns False only if 
            altLoc is empty.  Two methods will therefore return different lists
            whenever there are residues that have non-empty altLoc, even if the other
            conformation is not present.  While this is somewhat of an unusual situation,
            matchAlts() method will produce it when only one of the two matched structures
            has an alternate conformation in a residue. '''
        residues = [x for x in self.get_residues() if x.HasAltConf()]
        if fProteinOnly:
            residues = [x for x in residues if x.IsAminoAcid()]
        return list(residues.keys())

# --- Methods returning average values

    def GetAverageBfactor(self, what='all', listik=False, *args, **kwargs):
        ''' Returns the average B-factor of a group of atoms.  
            Selection parameters are the same as in atom_lister method.
        '''
        return self.GetListAverageB(self.atom_lister(what, listik, *args, **kwargs))

    def GetChainAverageBfactor(self, what='all', listik=False, *args, **kwargs):
        ''' Returns a dictionary of average B-factors for each chain.
            Selection parameters are the same as in atom_lister method.
        '''
        return dict([(chid,self.GetListAverageB(chind)) for (chid,chind) in self.ListChainSplit(self.atom_lister(what, listik, *args, **kwargs)).items()])

    def GetListAverageB(self, listik):
        bvals = self.GetBvector(listik=listik)
        ovals = self.GetOccupancyVector(listik)
        return sum(bvals*ovals)/sum(ovals)

    def get_listed_altconf(self, listik):
        return list(set([a.GetAltLoc() for a in self.atom_getter(listik=listik)]))

    def GetCrystLine(self):
        ''' Return the CRYST1 pdb record for the molecule. '''
        return '' if self.cell == None else self.cell.GetLine()

    def GetSymmOps(self):
        return SpaceGroups.symops[self.cell.GetSG()] if self.cell else None

    def GetSymmNum(self):
        return len(SpaceGroups.symops[self.cell.GetSG()]) if self.cell else 0

    def GetFractionalCoordinates(self, r):
        return dot(self.cell.GetMcf(),r.T).T if self.cell else None

    def Fractionalize(self):
        if self.cartesian:
            if self.cell:
                M = self.cell.GetMcf()
                for atom in self.atoms:
                    atom.transform(M)
                self.cartesian = False

    def Cartesianize(self):
        if not self.cartesian:
            if self.cell:
                M = self.cell.GetMfc()
                for atom in self.atoms:
                    atom.transform(M)
                self.cartesian = True

    def ApplySymmetryOperator(self, i):
        if i>0 and i < self.GetSymmNum():
            symop = self.GetSymmOps()[i]
            self.Fractionalize()
            for atom in self.atoms:
                atom.ApplyOperator(symop)
            self.Cartesianize()

    def CreateSymate(self, i, cell_shift=array([0.0,0.0,0.0])):
        if i==0:
            return self.copy()
        elif i>0 and i < self.GetSymmNum():
            x = self.copy()
            x.ApplySymmetryOperator(i)
            if array(cell_shift).any():
                x.CellShift(cell_shift)
            return x
        return None

    def GetCoordinateArray(self, what='all', listik=False, *args, **kwargs):
        ''' Return the Nx3 array of individual atom coordinates. '''
        return array([x.GetR() for x in self.atom_getter(what, listik, *args, **kwargs)])

    def SymateDistance(self, i, cell_shift=array([0.0,0.0,0.0])):
        symate = self.CreateSymate(i, cell_shift)
        return self.ShortestApproach(symate) if symate else None

    def ShortestApproach(self, other):
        xyz0 = self.GetCoordinateArray()
        xyz1 = other.GetCoordinateArray()
        N = self.GetAtomNumber()
        return sqrt(array((matrix(xyz1.T[0]).T*ones(N)).T-matrix(xyz0.T[0]).T*ones(N))**2+array((matrix(xyz1.T[1]).T*ones(N)).T-matrix(xyz0.T[1]).T*ones(N))**2+array((matrix(xyz1.T[2]).T*ones(N)).T-matrix(xyz0.T[2]).T*ones(N))**2).min()

    def ShortestDistances(self, other):
        xyz0 = self.GetCoordinateArray()
        xyz1 = other.GetCoordinateArray()
        N0, N1 = self.GetAtomNumber(), other.GetAtomNumber()
        return sqrt(array((matrix(xyz1.T[0]).T*ones(N0)).T-matrix(xyz0.T[0]).T*ones(N1))**2 + array((matrix(xyz1.T[1]).T*ones(N0)).T-matrix(xyz0.T[1]).T*ones(N1))**2 + array((matrix(xyz1.T[2]).T*ones(N0)).T-matrix(xyz0.T[2]).T*ones(N1))**2).min(1)

    def Distances2Bfactors(self, other):
        self.SetBfactorValues(bvalues=self.ShortestDistances(other))

    def SymNeighbors(self, cutoff=4.0, celldepth=1, verbose=False):
        self.__CellCheck_()
        symlist, symates = [], []
        for i in range(self.GetSymmNum()):
            symates.append(self.CreateSymate(i))
        for symate in symates:
            symate.MoveCenter2UnitCell()
        origin = symates[0].copy()
        back_shift = floor(self.GetCenter(True))-floor(origin.GetCenter(True))
        for a in range(-celldepth,celldepth+1):
            for b in range(-celldepth,celldepth+1):
                for c in range(-celldepth,celldepth+1):
                    for (i,symate) in enumerate(symates):
                        symate.CellShift([a,b,c])
                        D = origin.ShortestApproach(symate)
                        if D <= cutoff:
                            symcopy = symate.copy()
                            symcopy.CellShift(back_shift)
                            symlist.append([i,[a,b,c],D,symcopy])
                        if verbose:
                            if D <= cutoff:
                                print('#%2d + [%3d %3d %3d] -- %6.1f <-- include' % (i,a,b,c,D))
                            else:
                                print('#%2d + [%3d %3d %3d] -- %6.1f ' % (i,a,b,c,D))
                        symate.CellShift([-a,-b,-c])
        return symlist

    def WriteSymNeighbors(self, pdbfile, cutoff=4.0, celldepth=1, verbose=False, nmrstyle=False):
        symlist = self.SymNeighbors(cutoff=cutoff, celldepth=celldepth, verbose=verbose)
        if nmrstyle:
            WriteNMRstyle(zip(*symlist)[3],pdbfile,self.GetCrystLine())
        else:
            for (i,model) in enumerate(zip(*symlist)[3]):
                model.writePDB(pdbfile+str(i)+os.extsep+'pdb',self.GetCrystLine())

    def GetCenter(self, fractional=False):
        if fractional and self.cell:
            self.Fractionalize()
            retval = self.GetCoordinateArray().mean(0)
            self.Cartesianize()
        else:
            retval = self.GetCoordinateArray().mean(0)
        return retval      

    def WriteSymate(self, i, pdbfile):
        x = self.CreateSymate(i)
        if x:
            x.writePDB(pdbfile)

    def WriteSymates(self, pdbfile, move2cell=True):
        self.__CellCheck_()
        for i in range(self.GetSymmNum()):
            symate = self.CreateSymate(i)
            if move2cell:
                symate.MoveCenter2UnitCell()
            symate.writePDB(pdbfile+str(i)+os.extsep+'pdb')

    def __CellCheck_(self, line='Attempt to apply symmetry failed - no unit cell specified.'):
        assert self.cell, line

    def Move2UnitCell(self):
        self.__CellCheck_()
        self.Fractionalize()
        for atom in self.atoms:
            atom.xyz %= 1
        self.Cartesianize()

    def MoveCenter2UnitCell(self):
        self.__CellCheck_()
        vec = -floor(self.GetCenter(True))
        self.CellShift(vec)
        return vec

    def CellShift(self, vec=array([0.0,0.0,0.0])):
        self.__CellCheck_()
        if array(vec).any():
            self.Fractionalize()
            for atom in self.atoms:
                atom.shift(array(vec))
            self.Cartesianize()

    def distance(self,i,j):
        return sqrt(((self.atoms[i].xyz-self.atoms[j].xyz)**2).sum())

    def select(self,name=None,altloc=None,resn=None,resid=None):
        return [i for i,atom in self.__enumerate_atoms_() if atom.test(name,altloc,resn,resid)]

    def idselect(self, atid):
        return self.select(name=atid.GetAtomName(), resid=atid.ResID(), altloc=atid.GetAltLoc())[0]

    def angle(self,i,j,k):
        r1 = self.r(j,i)
        r2 = self.r(j,k)
        return degrees(arccos((r1*r2).sum()/sqrt((r1**2).sum()*(r2**2).sum())))

    def r(self,i,j):
        return self.atoms[j].xyz - self.atoms[i].xyz

    def dx(self,i,j):
        return self.atoms[j].xyz[0] - self.atoms[i].xyz[0]
    def dy(self,i,j):
        return self.atoms[j].xyz[1] - self.atoms[i].xyz[1]
    def dz(self,i,j):
        return self.atoms[j].xyz[2] - self.atoms[i].xyz[2]

    def torsion(self,i,j,k,l):
        b1 = self.r(i,j)
        b2 = self.r(j,k)
        b3 = self.r(k,l)
        b2xb1 = cross(b2,b1)
        x = (b3*cross(b2,b2xb1)).sum()
        y = sqrt((b2**2).sum())*(b3*b2xb1).sum()
        return degrees(arctan2(-y,x))

    def pick(self,atomid):
        ind = None
        for (i,atom) in enumerate(self.atoms):
            if atom.atomid() == atomid:
                ind = i
                break
        return ind

    def environment(self, atomi, rmax=3.2, polaronly=True):
        ''' Returns the environment of the atom.  List of the contact atoms is
            generated using maximum distance of rmax and only include polar atoms
            by default (use polaronly=False to include carbons).
            Method returns the dictionary of bonds.  Keys are atom IDs and values
            are bond distances.
            Intra-residue contacts are included.  Use first_shell() method to
            obtain the inter-residue contacts only.'''
        atom = self.atoms[atomi]
        listik = self.atom_lister('polar') if polaronly else self.atom_lister()
        listik = [i for (i,a) in self.__enumerate_atoms_(listik) if abs(atom.dx(a))<rmax]
        listik = [i for (i,a) in self.__enumerate_atoms_(listik) if abs(atom.dy(a))<rmax]
        listik = [i for (i,a) in self.__enumerate_atoms_(listik) if abs(atom.dz(a))<rmax]
        D = [(i,atom.distance(a)) for (i,a) in self.__enumerate_atoms_(listik)]
        bonds = dict([(i,x) for (i,x) in D if x < rmax])
        bonds.pop(atomi,0)
        return bonds

    def same_residue(self, atomi, atomj):
        ''' Returns True if atoms belong to the same residue (including
            alternate conformers), False otherwise. '''
        return self.atoms[atomi].same_residue(self.atoms[atomj])

    def first_shell(self, atomi, rmax=3.2, polaronly=True):
        ''' This is similar to environment() method, but excludes the intra-residue
            contacts (e.g. covalent bonds). '''
        return dict([(i,x) for (i,x) in self.environment(atomi, rmax, polaronly).items() if not self.same_residue(atomi,i)])

    def GetCoM(self, listik=False):
        mo = self.GetMOVector(listik=listik)
        r = self.GetCoordinateArray(listik=listik)
        return (mo*r.T).sum(1)/sum(mo)

    def Rgyration(self, listik=False):
        mo = self.GetMOVector(listik=listik).T
        r = self.GetCoordinateArray(listik=listik).T
        Rcenter = (mo*r).sum(1)/mo.sum()
        return sqrt(sum(mo*((r.T-Rcenter)**2).sum(1))/sum(mo))

    def GetInertiaTensor(self, what='all', listik=False, *args, **kwargs):
        listik = self.atom_lister(what, listik, *args, **kwargs)
        mo = self.GetMOVector(listik=listik)
        x, y, z = self.GetCoordinateArray(listik=listik).T
        return TInertia(mo, x, y, z)

    def WriteNeighborResidues(self, pdbFile, resid, rmax=4.0, self_exclude=False, header=None):
        residlist = self.resid_lister(listik=self.atom_lister('vicinity',rcutoff=rmax,corelist=self.atom_lister('resid',resid=resid)))
        if not self_exclude:
            residlist.append(resid)
        self.WriteResidueList(pdbFile=pdbFile, residlist=residlist, header=header)

    def WriteResidueList(self, pdbFile, residlist, header=None):
        with open(pdbFile, 'w') as fout:
            if header == 'cell':
                self.__headwrite_(fout, self.GetCrystLine())
            else:
                self.__headwrite_(fout, header)
            for atom in self.atom_getter('resids',resids=residlist):
                fout.write(atom.GetAtomRecord())
            fout.write('END   \n')

    def WriteAtomList(self, pdbFile, atomlist, header=None):
        with open(pdbFile, 'w') as fout:
            if header == 'cell':
                self.__headwrite_(fout, self.GetCrystLine())
            else:
                self.__headwrite_(fout, header)
            for atomi in sorted(atomlist):
                fout.write(self.GetAtom(atomi).GetAtomRecord())
            fout.write('END   \n')

    def get_residues(self, what='all', listik=False, *args, **kwargs):
        ''' Returns the dictionary of residues.  Atom selection is the
            same as in atom_lister method. Notice that residues will be 
            completed in case your atom selection is partial for a
            particular residue.
        '''
        atoms = self.atom_getter(listik=self.ListCompleteResidues(self.atom_lister(what, listik, *args, **kwargs)))
        residues = {}
        for atom in atoms:
            resid = atom.resid()
            if resid not in residues.keys():
                residues[resid] = [atom]
            else:
                residues[resid].append(atom)
        return dict([(resid,pdbresidue(ratoms)) for resid,ratoms in residues.items()])

    def GetResidueNames(self, listik=False):
        return dict([(a.resid(),a.get_res_name()) for a in self.__ensure_atoms_(listik)])

    def GetChains(self):
        ''' Returns the dictionary of chains with chain IDs as keys and
            number of atoms in each chain as values. '''
        return Counter([a.chainID() for a in self.atom_getter()])

    def GetChainIDs(self):
        ''' Returns the sorted list of chain IDs present.'''
        return sorted(set([a.chainID() for a in self.atom_getter()]))

    def rename_chain(self, chid1, chid2, forced=False, what='all', listik=False, *args, **kwargs):
        ''' Changes chain ID for atoms in in chid1 to chid2.  Asserts
            that chid2 doesn't exist unless forced.'''
        if chid1 != chid2:
            if not forced:
                assert chid2 not in self.GetChainIDs(), 'Chain '+chid2+' already present, exiting.'
            listik = self.atom_lister(what, listik, *args, **kwargs)
            for atom in self.atom_getter('chid', listik, chid=chid1):
                atom.SetChain(chid2)

    def __headwrite_(self, fout, header=None):
        if header is None:
            fout.write(self.GetCrystLine())
        else:
            tyhe = type(header)
            if tyhe is str:
                fout.write(header)
            elif tyhe is list or tyhe is tuple:
                for line in header:
                    fout.write(line)

    def writePDB(self, pdbFile, header=None, mode='w'):
        ''' Write the PDB file.  pdbFile is the path to the output file,
            header may contain single line, list or tuple  of lines 
            which will be placed at the beginning of the file.
        ''' 
        with open(pdbFile, mode) as fout:
            self.__headwrite_(fout, header)
            for atom in self.atoms:
                fout.write(atom.GetAtomRecord())
                fout.write(atom.GetAnisouRecord())
            fout.write('END   \n')

    def writePDBchains(self, pdbFile, chains, header=None):
        ''' Write the PDB file that only include the specified chains.
            pdbFile is the path to the output file. chains could be either
            a list or a string with selection chain IDs.'''
        with open(pdbFile, 'w') as fout:
            self.__headwrite_(fout, header)
            for atom in self.atom_getter('chids', chids=chains):
                fout.write(atom.GetAtomRecord())
                fout.write(atom.GetAnisouRecord())
            fout.write('END   \n')

    def extract_range(self, ranges):
        ''' Returns the copy of the molecule that only contains atoms from
            the supplied dictionary of ranges.The selection is defined by ranges
            dictionary, which has chain IDs as keys and list of tuples for 
            inclusive ranges of residue numbers.  For example, 
                ranges={'A': [(2,20), (27, 33), (45,45)],
                        'B': []}
            will select A2-A20, A27-A33, A45 and nothing from chain B (which
            is the same as omitting it).'''
        extracted_atoms = []
        for atom in self.atoms:
            chid, resn = atom.GetChain(), atom.resSeq()
            if atom.GetChain() in ranges:
                for limits in ranges[atom.GetChain()]:
                    if self.__range_checker_(resn, limits):
                        extracted_atoms.append(atom.copy())
                        break
        return pdbmolecule(atoms=extracted_atoms, cell=self.cell)

    def range_residues(self, ranges):
        resids = set([a.resid() for a in self.atom_getter() if a.chainID() in ranges])
        resids = [x for x in resids if sum([self.__range_checker_(int(x[1:]),y) for y in ranges.get(x[0])])]
        return resids

    def __range_checker_(self, value, limits):
        return value >= limits[0] and value <= limits[1]

    def extract_chains(self, chains):
        ''' Returns the copy of the molecule that only contains atoms from
            the supplied list of chains (the latter could be either list or
            string of symbols. '''
        extracted_atoms = [a.copy() for a in self.atom_getter('chids', chids=chains)]
        return pdbmolecule(atoms=extracted_atoms, cell=self.cell)

    def rename_chains(self, chains):
        ''' Rename chain IDs.  Supply the dictionary for renaming. '''
        for atom in self.atoms:
            chid = atom.GetChain()
            if chid in chains:
                atom.SetChain(chains[chid])

    def noise(self, xnoise=0.1, bnoise=0.1, occnoise=0.0):
        for atom in self.atoms:
            atom.noise(xnoise, bnoise, occnoise)

    def copy(self):
        return copy.deepcopy(self)

    def PhiPsiList(self):
        return backbone(self).PhiPsiList()

    def BackboneTorsions(self):
        return backbone(self).BackboneTorsions()

    def W2Wdistance(self, water1, water2):
        return water1.GetAtom('O').distance(water2.GetAtom('O'))

    def PolarContacts(self, rmax=3.2, include_interchain=True):
        ''' 
            Return the list of polar contacts. Excludes atoms from different 
            alternate conformers and same residues.  May exclude interchain
            contacts (say if you are not intersted in crystal contacts of a
            confirmed monomer) if include_interchain is set to false.
        '''
        return self.AllContacts(rmax, self.atom_lister('polar')) 

    def AllContacts(self, rmax=4.0, listik=False):
        '''
        Return the list of atoms and distances that are within rmax. 
        '''
        contacts = []
        listik = array(self.__ensure_listik_(listik))
        xyz = self.xyz(listik)
        rmax2 = rmax**2
        for i in range(len(listik)):
            r2 = ((xyz[i+1:]-xyz[i])**2).sum(1)
            ind = r2<=rmax2
            N = sum(ind)
            if sum(ind)>0:
                contacts += zip(*[(ones(N)*i).astype(int),listik[nonzero(ind)[0]+i+1],sqrt(r2[ind])])
        return contacts

    def IndexContacts(self, ind1, ind2, rmax=4.0):
        '''
        Return the list of atoms and distances that are within rmax.
        First atom is taken fron ind1 selection, second from ind2.
        '''
        contacts = []
        if len(ind2):
            xyz1, xyz2 = self.xyz(ind1), self.xyz(ind2)
            rmax2 = rmax**2
            for i in range(len(ind1)):
                r2 = ((xyz2-xyz1[i])**2).sum(1)
                ind = r2<=rmax2
                N = sum(ind)
                if sum(ind)>0:
                    contacts += zip(*[(ones(N)*ind1[i]).astype(int),array(ind2)[ind],sqrt(r2[ind])])
        return contacts

    def NeighborDatabase(self):
        return backbone(self).NeighborDatabase()

    def IndexDistances(self, indA, indB):
        d = []
        for atomi in indA:
            d.append([])
            for atomj in indB:
                d[-1].append(self.distance(atomi,atomj))
        return d

    def ChiList(self):
        chis = {}
        for residue in self.get_residues('protein').values():
            chichi = residue.GetChis()
            for chi in list(chichi.keys()):
                if not chichi[chi]:
                    chichi.pop(chi)
            if chichi:
                chis[residue.GetResID()+'|'+residue.get_res_name()] = chichi
            if residue.HasAltConf():
                for ac in residue.GetAltCodes()[1:]:
                    chichi = residue.GetAltChis(ac)
                    for chi in list(chichi.keys()):
                        if not chichi[chi]:
                            chichi.pop(chi)
                if chichi:
                    chis[residue.GetResID()+ac+'|'+residue.get_res_name()] = chichi
        return chis

    def bperesidue(self):
        for residue in self.get_residues().values():
            residue.Baverage()

    def GetBlimits(self):
        b = [a.GetB() for a in self.atoms]
        return (min(b),max(b))

    def get_occupancy_estimates(self):
        return dict([(resid,self.get_residues('resid',resid=resid)[resid].esitimate_ac_occupancies()) for resid in self.ListAltResidues()])

    def print_occupancy_estimates(self):
        occs = self.get_occupancy_estimates()
        for resid in sorted(occs):
            print(resid + ' %5.2f'*len(occs[resid]) % tuple(occs[resid]))

    def xyz(self, listik=False):
        return array([a.GetR() for a in self.atom_getter(listik=listik)])

class vpdbmolecule(pdbmolecule):
    '''
    Extension of pdbmolecule class that adds a lot of verbose methods.
    These are methods that can be replaced by a single line of code,
    but are provided for convenience.
    '''
    def ListAltConfTypes(self, listik=False):
        '''
        Returns the list of alternate conformer labels present in the
        molecule.
        '''
        return sorted(set([x.GetAltLoc() for x in self.atom_getter('altconf', listik=listik)]))

class backbone:
    '''
    Protein backbone atoms manipulation.
    '''
    def __init__(self, molecule):
        residues = molecule.get_residues('protein')
        self.segments = []
        self.residues = residues
        self.protons = {}
        for resid in residues:
            newFlag, tailFlag, joinFlag = True, True, False
            for (i,segment) in enumerate(self.segments):
                if self.join(residues[resid], segment[0]):
                    segment.insert(0, residues[resid])
                    newFlag = False
                    break
                elif self.join(segment[-1], residues[resid]):
                    segment.append(residues[resid])
                    newFlag, tailFlag = False, False
                    break
            if newFlag:
                self.segments.append([residues[resid]])
            else:
                if tailFlag:
                    for (j, segment) in enumerate(self.segments):
                        if self.join(segment[-1], residues[resid]):
                            joinFlag = True
                            break
                else:
                    for (j, segment) in enumerate(self.segments):
                        if self.join(residues[resid], segment[0]):
                            joinFlag = True
                            break
                if joinFlag:
                    if tailFlag:
                        self.segments[j].extend(self.segments[i])
                        self.segments.pop(i)
                    else:
                        self.segments[i].extend(self.segments[j])
                        self.segments.pop(j)
        for segment in self.segments:
            for (i,residue) in enumerate(segment[1:]):
                atomC  = segment[i].GetAtom('C')
                atomCA = residue.GetAtom('CA')
                atomN  = residue.GetAtom('N')
                if atomC and atomCA and atomN:
                    self.protons[residue.GetResID()] =  0.7792*(atomN.GetXYZarray()-atomC.GetXYZarray()) + 0.7029*(atomN.GetXYZarray()-atomCA.GetXYZarray()) + atomN.GetXYZarray()
        self.resegment = {}
        for (i,segment) in enumerate(self.segments):
            for (j,residue) in enumerate(segment):
                self.resegment[residue.GetResID()] = (i,j)
                
    def PredictAmideProton(self, resid):
        return self.protons[resid]

    def GetResegment(self):
        return self.resegment

    def MChbonds(self, OHmax = 2.5, ONmax = 3.2, DHAmin=90.0):
        hbonds = {}
        for resN, proton in self.protons.items():
            for resO, oxyresidue in self.residues.items():
                oxygen = oxyresidue.GetAtom('O').GetR()
                d1 = distance2(oxygen,proton)
                if d1 < OHmax:
                    nitrogen = self.residues[resN].GetAtom('N').GetR()
                    d2 = distance2(oxygen, nitrogen)
                    if d2 < ONmax:
                        if resN != resO:
                            DHA = angle2(nitrogen,proton,oxygen)
                            if DHA > DHAmin:
                                hbonds[resN+resO] = [d1, d2, DHA]
                    break
        return hbonds

    def get_mchbond(self, res1, res2):
        '''
            Returns the array of the parameters of mainchain hydrogen
            bond from res1 (N donor) to res2 (O acceptor).  Array
            includes predicted O...H distance, O-N distance, and 
            predicted N-H...O angle.
        '''
        proton = self.protons[res1]
        nitrogen = self.residues[res1].GetAtom('N').GetR()
        oxygen = self.residues[res2].GetAtom('O').GetR()
        return array([distance2(oxygen,proton), distance2(oxygen, nitrogen), angle2(nitrogen,proton,oxygen)])

    def Twist(self):
        twists = []
        for segment in self.segments:
            twists.append([])
            for (i,residue) in enumerate(segment[1:]):
                twist = torsion(segment[i].GetAtom('O'), segment[i].GetAtom('C'), residue.GetAtom('C'), residue.GetAtom('O'))
                if twist <= 0:
                    twists[-1].append(180+twist)
                else:
                    twists[-1].append(twist-180)
        return twists

    def SegmentTwist(self, resid1, resid2):
        return sum(self.BetaTwist(resid1, resid2))

    def BetaTwist(self, resid1, resid2):
        (seg1,res1) = self.resegment[resid1]
        (seg2,res2) = self.resegment[resid2]
        if seg1 != seg2:
            return None
        twists, segment = [], self.segments[seg1]
        for i in range(res1,res2):
            twist = torsion(segment[i].GetAtom('O'), segment[i].GetAtom('C'), segment[i+1].GetAtom('C'), segment[i+1].GetAtom('O'))
            if twist <= 0:
                twists.append(180+twist)
            else:
                twists.append(twist-180)
        return twists

    def GetSegment(self, Segment):
        (resid1, resid2) = Segment
        (seg1,res1) = self.resegment[resid1]
        (seg2,res2) = self.resegment[resid2]
        residues = []
        for segment in range(seg1,seg2+1):
            if segment == seg1:
                resA = res1
            else:
                resA = 0
            if segment == seg2:
                resB = res2
            else:
                resB = len(self.segments[segment])-1
            residues.append(self.segments[segment][resA])
            for res in range(resA,resB):
                residues.append(self.segments[segment][res+1])
        return residues

    def Residue2SegmentDistance(self, resid1, resid2, resid3, fullReturn=False):
        r0 = self.residues[resid1].GetCA().GetXYZarray()
        (seg1,res1) = self.resegment[resid2]
        (seg2,res2) = self.resegment[resid3]
        d, names = [], []
        for segment in range(seg1,seg2+1):
            if segment == seg1:
                resA = res1
            else:
                resA = 0
            if segment == seg2:
                resB = res2
            else:
                resB = len(self.segments[segment])-1
            cres = self.segments[segment][resA]
            r1 = cres.GetCA().GetXYZarray()
            names.append(cres.get_res_name()+' '+cres.GetResID())
            for res in range(resA,resB):
                cres = self.segments[segment][res+1]
                r2 = cres.GetCA().GetXYZarray()
                d.append(point2segment(r1,r2,r0))
                names.append(cres.get_res_name()+' '+cres.GetResID())
                r1=r2
        if fullReturn:
            return (d, names)
        else:
            return [min(d)]

    def Segment2SegmentDistance(self, Segment1, Segment2, returnType=0):
        (seg11, res11), (seg12, res12) = self.resegment[Segment1[0]], self.resegment[Segment1[1]]
        (seg21, res21), (seg22, res22) = self.resegment[Segment2[0]], self.resegment[Segment2[1]]
        d = []
        for seg1 in range(seg11, seg12+1):
            if seg1 == seg11:
                start1 = res11
            else:
                start1 = 0
            if seg1 == seg12:
                finis1 = res12
            else:
                finis1 = len(self.segments[seg1])-1
            r11 = self.segments[seg1][start1].GetCA().GetXYZarray()
            for r1 in range(start1, finis1):
                d.append([])
                r12 = self.segments[seg1][r1+1].GetCA().GetXYZarray()
                for seg2 in range(seg21, seg22+1):
                    if seg2 == seg21:
                        start2 = res21
                    else:
                        start2 = 0
                    if seg2 == seg22:
                        finis2 = res22
                    else:
                        finis2 = len(self.segments[seg2])-1
                    r21 = self.segments[seg2][start2].GetCA().GetXYZarray()
                    for r2 in range(start2, finis2):
                        r22 = self.segments[seg2][r2+1].GetCA().GetXYZarray()
                        d[-1].append(segment2segment((r11,r12),(r21,r22)))
                        r21 = r22
                r11 = r12
        if returnType==1:
            dd = []
            for i in range(len(d)):
                dd.append(min(d[i]))
            return (dd, self.GetSegment(Segment1))
        elif returnType==2:
            dd = []
            dzip = list(zip(*d))
            for i in range(len(dzip)):
                dd.append(min(dzip[i]))
            return (dd,self.GetSegment(Segment2))
        else:
            return (d, (self.GetSegment(Segment1), self.GetSegment(Segment2)))
        

    def join(self, tail, head):
        c, n = tail.GetAtom('C'), head.GetAtom('N')
        if c and n:
            return distance(c, n) < 2.5
        else:
            return False

    def PhiPsiList(self):
        phi, psi = {}, {}
        for segment in self.segments:
            if len(segment) > 1:
                residue = segment[0]
                n, ca, c, n_ = (residue.GetAtom('N'), residue.GetAtom('CA'), residue.GetAtom('C'), segment[1].GetAtom('N'))
                if n and ca and c and n_:
                    psi[residue.GetResID()] = torsion(n, ca, c, n_)
                else:
                    logging.info('Missing atoms in \''+residue.GetResID()+'\' - no psi angle calculated')
                for (i,residue) in enumerate(segment[1:-1]):
                    resid = residue.GetResID()
                    (c_, n, ca, c, n_) =  (segment[i].GetAtom('C'), residue.GetAtom('N'), residue.GetAtom('CA'), residue.GetAtom('C'), segment[i+2].GetAtom('N'))
                    if c_ and n and ca and c:
                        phi[resid] = torsion(c_, n, ca, c)
                    else:
                        logging.info('Missing atoms in \''+resid+'\' - no phi angle calculated')
                    if n and ca and c and n_:
                        psi[resid] = torsion(n, ca, c, n_)
                    else:
                        logging.info('Missing atoms in \''+resid+'\' - no psi angle calculated')
                residue = segment[-1]
                (c_, n, ca, c) =  (segment[-2].GetAtom('C'), residue.GetAtom('N'), residue.GetAtom('CA'), residue.GetAtom('C'))
                if c_ and n and ca and c:
                    phi[residue.GetResID()] = torsion(c_, n, ca, c)
                else:
                    logging.info('Missing atoms in \''+residue.GetResID()+'\' - no phi angle calculated')
        return (phi, psi)

    def BackboneTorsions(self):
        phi, psi, omega = {}, {}, {}
        for segment in self.segments:
            if len(segment) > 1:
                residue = segment[0]
                n, ca, c, o, n_, ca_ = (residue.GetAtom('N'), residue.GetAtom('CA'), residue.GetAtom('C'), residue.GetAtom('O'), segment[1].GetAtom('N'), segment[1].GetAtom('CA'))
                if n and ca and c and n_:
                    psi[residue.GetResID()] = torsion(n, ca, c, n_)
                else:
                    logging.info('Missing atoms in \''+residue.GetResID()+'\' - no psi angle calculated')
                if o and c and n_ and ca_:
                    omega[residue.GetResID()] = torsion(o, c, n_, ca_)
                for (i,residue) in enumerate(segment[1:-1]):
                    resid = residue.GetResID()
                    (c_, n, ca, c, o, n_, ca_) =  (segment[i].GetAtom('C'), residue.GetAtom('N'), residue.GetAtom('CA'), residue.GetAtom('C'), residue.GetAtom('O'), segment[i+2].GetAtom('N'), segment[i+2].GetAtom('CA'))
                    if c_ and n and ca and c:
                        phi[resid] = torsion(c_, n, ca, c)
                    else:
                        logging.info('Missing atoms in \''+resid+'\' - no phi angle calculated')
                    if n and ca and c and n_:
                        psi[resid] = torsion(n, ca, c, n_)
                    else:
                        logging.info('Missing atoms in \''+resid+'\' - no psi angle calculated')
                    if o and c and n_ and ca_:
                        omega[resid] = torsion(o, c, n_, ca_)
                    else:
                        logging.info('Missing atoms in \''+resid+'\' - no omega angle calculated')
                residue = segment[-1]
                (c_, n, ca, c) =  (segment[-2].GetAtom('C'), residue.GetAtom('N'), residue.GetAtom('CA'), residue.GetAtom('C'))
                if c_ and n and ca and c:
                    phi[residue.GetResID()] = torsion(c_, n, ca, c)
                else:
                    logging.info('Missing atoms in \''+residue.GetResID()+'\' - no phi angle calculated')
        return (phi, psi, omega)

    def NeighborDatabase(self):
        ndb = {}
        for segment in self.segments:
            resids = []
            for residue in segment:
                resids.append(residue.GetResID())
            ndb[resids[0]] = [resids[1]]
            for (i, resid) in enumerate(resids[1:-1]):
                ndb[resid] = [resids[i],resids[i+2]]
            ndb[resids[-1]] = resids[-2]
        return ndb

#--- check below the line
class pdbresidue:

    def __init__(self, atoms):
        self.origatoms = atoms
        self.atoms, self.altcodes, self.altatoms = {}, [], []
        for atom in atoms:
            if atom.HasAltConf():
                ac = atom.GetAltLoc()
                if ac not in self.altcodes:
                    self.altcodes.append(ac)
                    self.altatoms.append({})
        if not len(self.altcodes):
            self.altatoms.append({})
        for atom in atoms:
            if atom.HasAltConf():
                ac = atom.GetAltLoc()
                self.altatoms[self.altcodes.index(ac)][atom.GetName()] = atom
            else:
                for aconf in self.altatoms:
                    aconf[atom.GetName()] = atom
        self.atoms = self.altatoms[0]
        self.name = atom.get_res_name()
        self.resid = atom.GetResID()
    def r(self,i,j):
        return self.origatoms[j].xyz - self.origatoms[i].xyz
    def distance(self,i,j):
        return math.sqrt(((self.origatoms[i].xyz-self.origatoms[j].xyz)**2).sum())
    def angle(self,i,j,k):
        r1 = self.r(j,i)
        r2 = self.r(j,k)
        return math.degrees(arccos((r1*r2).sum()/math.sqrt((r1**2).sum()*(r2**2).sum())))
    def torsion(self,i,j,k,l):
        b1 = self.r(i,j)
        b2 = self.r(j,k)
        b3 = self.r(k,l)
        b2xb1 = cross(b2,b1)
        x = (b3*cross(b2,b2xb1)).sum()
        y = math.sqrt((b2**2).sum())*(b3*b2xb1).sum()
        return math.degrees(math.atan2(-y,x))

    def rmsdAlt(self, other):
        r = []
        for altgroup1 in self.altatoms:
            r.append([])
            for altgroup2 in other.altatoms:
                r2 = []
                for key in altgroup1:
                    if key in altgroup2:
                        r2.append(altgroup1[key].GetXYZarray()-altgroup2[key].GetXYZarray())
                r[-1].append(sqrt(3*(array(r2)**2).mean()))
        return array(r)

    def matchAlt2Alt(self, other):
        a = self.rmsdAlt(other)
        ind1 = list(range(len(self.altcodes)))
        ind2 = list(range(len(other.altcodes)))
        matchlist = []
        while ind1 or ind2:
            i = ind1[a[ind1,ind2].squeeze().min(-1).argmin()]
            j = ind2[a[i].squeeze()[ind2].argmin()]
            matchlist.append((self.altcodes[i], other.altcodes[j]))
            ind1.remove(i)
            ind2.remove(j)
        return matchlist

    def matchSingle2Alt(self, other):
        return other.altcodes[self.rmsdAlt(other).argmin()]

    def matchAlt2Single(self, other):
        return self.altcodes[other.rmsdAlt(self).argmin()]

    def renameAlts(self, codes):
        codict=dict(codes)
        for atom in self.origatoms:
            if atom.altLoc() in codict:
                atom.SetAltLoc(codict[atom.altLoc()])
        self.__init__(self.origatoms)

    def renameAlt2Alt(self, other):
        self.renameAlts(self.matchAlt2Alt(other))

    def renameAlt2Single(self, other):
        other.AssignAlts(self.matchAlt2Single(other))

    def AssignAlts(self, ac):
        for atom in self.origatoms:
            atom.SetAltLoc(ac)
        self.__init__(self.origatoms)

    def renameSingle2Alt(self, other):
        self.AssignAlts(self.matchSingle2Alt(other))

    def get_res_name(self):
        return self.name

    def GetResID(self):
        return self.resid

    def GetTitle(self):
        return self.name+self.resid[0]+self.resid[1:].strip()

    def HasAltConf(self):
        return bool(len(self.altatoms)-1)

    def GetChis(self, purge=False):
        chis = {}
        # chi1     side chain torsion angle for atoms N,CA,CB,*G.         
        chis['chi1']  = self.GetChi(('N$','CA$','CB$','.G$'))
        # chi11     side chain torsion angle for atoms N,CA,CB,*G1. 
        chis['chi11'] = self.GetChi(('N$','CA$','CB$','.G1$'))
        # chi12     side chain torsion angle for atoms N,CA,CB,*G2.
        chis['chi12'] = self.GetChi(('N$','CA$','CB$','.G2$'))
        # chi2     side chain torsion angle for atoms CA,CB,*G,*D.
        # (includes exception for ILE CG1)
        chis['chi2'] = self.GetChi(('CA$','CB$','.G$|CG1$','.D$'))
        # chi21     side chain torsion angle for atoms CA,CB,*G,*D1. 
        chis['chi21'] = self.GetChi(('CA$','CB$','.G$','.D1$'))
        # chi22     side chain torsion angle for atoms CA,CB,*G,*D2. 
        chis['chi22'] = self.GetChi(('CA$','CB$','.G$','.D2$'))
        # chi3     side chain torsion angle for atoms CB,*G,*D,*E.
        chis['chi3'] = self.GetChi(('CB$','.G$','.D$','.E$'))
        # chi31     side chain torsion angle for atoms CB,*G,*D,*E1.
        chis['chi31'] = self.GetChi(('CB$','.G$|OG1$','.D$','.E1$'))
        # chi32     side chain torsion angle for atoms CB,*G,*D,*E2.
        chis['chi32'] = self.GetChi(('CB$','.G$|OG1$','.D$','.E2$'))
        # chi33    defined for TPO O3P
        chis['chi33'] = self.GetChi(('CB$','OG1$','^P$','O3P$'))
        # chi4     side chain torsion angle for atoms *G,*D,*E,*Z.
        chis['chi4'] = self.GetChi(('.G$','.D$','.E$','.Z$'))
        # chi51     side chain torsion angle for atoms *D,*E,*Z, NH1.
        chis['chi51'] = self.GetChi(('.D$','.E$','.Z$','NH1$'))
        # chi52     side chain torsion angle for atoms *D,*E,*Z, NH2. 
        chis['chi52'] = self.GetChi(('.D$','.E$','.Z$','NH2$'))
        # Phosphothreonine
        if self.name == 'TPO':
            chis['chi2']  = self.GetChi(('CA$','CB$','OG1$','^P$'))
            chis['chi31'] = self.GetChi(('CB$','OG1$','^P$','O1P$'))
            chis['chi32'] = self.GetChi(('CB$','OG1$','^P$','O2P$'))
            chis['chi33'] = self.GetChi(('CB$','OG1$','^P$','O3P$'))
        if self.name == 'PTR':
            chis['chi61'] = self.GetChi(('CE1$','CZ$','OH$','^P$'))
            chis['chi62'] = self.GetChi(('CE2$','CZ$','OH$','^P$'))
            chis['chi71'] = self.GetChi(('CZ$','OH$','^P$','O1P$'))
            chis['chi72'] = self.GetChi(('CZ$','OH$','^P$','O2P$'))
            chis['chi73'] = self.GetChi(('CZ$','OH$','^P$','O3P$'))
        if purge:
            keys = list(chis.keys())
            for key in keys:
                if not chis[key]:
                    del chis[key]
        return chis

    def GetAltChis(self, code, purge=False):
        chis = {}
        chis['chi1']  = self.GetAltChi(('N$','CA$','CB$','.G$'),code)
        chis['chi11'] = self.GetAltChi(('N$','CA$','CB$','.G1$'),code)
        chis['chi12'] = self.GetAltChi(('N$','CA$','CB$','.G2$'),code)
        chis['chi2'] = self.GetAltChi(('CA$','CB$','.G$|CG1$','.D$'),code)
        chis['chi21'] = self.GetAltChi(('CA$','CB$','.G$','.D1$'),code)
        chis['chi22'] = self.GetAltChi(('CA$','CB$','.G$','.D2$'),code)
        chis['chi3'] = self.GetAltChi(('CB$','.G$','.D$','.E$'),code)
        chis['chi31'] = self.GetAltChi(('CB$','.G$','.D$','.E1$'),code)
        chis['chi32'] = self.GetAltChi(('CB$','.G$','.D$','.E2$'),code)
        chis['chi4'] = self.GetAltChi(('.G$','.D$','.E$','.Z$'),code)
        chis['chi51'] = self.GetAltChi(('.D$','.E$','.Z$','NH1$'),code)
        chis['chi52'] = self.GetAltChi(('.D$','.E$','.Z$','NH2$'),code)
        if self.name == 'TPO':
            chis['chi2']  = self.GetAltChi(('CA$','CB$','OG1$','^P$'),code)
            chis['chi31'] = self.GetAltChi(('CB$','OG1$','^P$','O1P$'),code)
            chis['chi32'] = self.GetAltChi(('CB$','OG1$','^P$','O2P$'),code)
            chis['chi33'] = self.GetAltChi(('CB$','OG1$','^P$','O3P$'),code)
        if self.name == 'PTR':
            chis['chi61'] = self.GetAltChi(('CE1$','CZ$','OH$','^P$'),code)
            chis['chi62'] = self.GetAltChi(('CE2$','CZ$','OH$','^P$'),code)
            chis['chi71'] = self.GetAltChi(('CZ$','OH$','^P$','O1P$'),code)
            chis['chi72'] = self.GetAltChi(('CZ$','OH$','^P$','O2P$'),code)
            chis['chi73'] = self.GetAltChi(('CZ$','OH$','^P$','O3P$'),code)
        if purge:
            keys = chis.keys()
            for key in keys:
                if not chis[key]:
                    del chis[key]
        return chis

    def GetChi(self,masks):
        chi = []
        atomgroups = (self.GetMaskAtoms(masks[0]), self.GetMaskAtoms(masks[1]), self.GetMaskAtoms(masks[2]), self.GetMaskAtoms(masks[3]))
        for i in atomgroups[0]:
            iName = i.GetName()
            for j in atomgroups[1]:
                jName = j.GetName()
                for k in atomgroups[2]:
                    kName = k.GetName()
                    for l in atomgroups[3]:
                        lName = l.GetName()
                        chi.append((iName+'-'+jName+'-'+kName+'-'+lName, torsion(i,j,k,l)))
        return chi

    def GetAltChi(self,masks,code):
        chi = []
        atomgroups = (self.GetMaskAltAtoms(masks[0],code), self.GetMaskAltAtoms(masks[1],code), self.GetMaskAltAtoms(masks[2],code), self.GetMaskAltAtoms(masks[3],code))
        for i in atomgroups[0]:
            iName = i.GetName()
            for j in atomgroups[1]:
                jName = j.GetName()
                for k in atomgroups[2]:
                    kName = k.GetName()
                    for l in atomgroups[3]:
                        lName = l.GetName()
                        chi.append((iName+'-'+jName+'-'+kName+'-'+lName, torsion(i,j,k,l)))
        return chi

    def GetMaskAtoms(self, mask):
        return [self.atoms[name] for name in self.atoms.keys() if re.match(mask,name)]

    def GetMaskAltAtoms(self, mask, code):
        ac = self.altcodes.index(code)
        atoms = [self.altatoms[ac][name] for name in self.altatoms[ac].keys() if re.match(mask,name)]
        return atoms if atoms else self.GetMaskAtoms(mask)

    def GetAtom(self, name):
        try:
            return self.atoms[name]
        except KeyError:
            return None

    def GetAtoms(self):
        return list(self.atoms.values())

    def get_atom_names(self):
        ''' Return the list of atom names in the residue, sorted in 
            random order. '''
        return list(self.atoms.keys())

    def get_atom_names_sorted(self):
        ''' Return the list of atom names in the residue, sorted alphabetically.
        '''
        return sorted(self.atoms.keys())

    def get_atom_number(self):
        return len(self.atoms)

    def IsAminoAcid(self):
        return pdbnames.Is3Amino(self.name)

    def GetCA(self):
        return self.atoms['CA']

    def GetAtomR(self, name):
        return self.atoms[name].GetXYZarray()

    def SetAtomR(self, name, vector):
        self.atoms[name].SetXYZarray(vector)

    def GetAltCodes(self):
        return self.altcodes

    def GetAltNum(self):
        return len(self.altcodes)

    def get_b_vectors(self):
        ''' Returns the array that contains vectors of B-factors for all
            the conformers (only one vector when the residue has only
            one conformer).  Atoms are sorted by name, use 
            get_atom_names_sorted() method  to get the list of atoms.'''  
        bvec = []
        for atmgrp in self.altatoms:
            bvec.append([])
            for name in sorted(atmgrp):
                bvec[-1].append(atmgrp[name].GetB())
        return array(bvec)

    def get_r_array(self):
        ''' Returns the array of atomic coordinates for each atom.  
            Atoms are sorted by name, use get_atom_names_sorted() method 
            to get the list of atoms.'''  
        rvec = []
        for atmgrp in self.altatoms:
            rvec.append([])
            for name in sorted(atmgrp):
                rvec[-1].append(atmgrp[name].GetR())
        return array(rvec)

    def get_alt_distances(self):
        ''' Returns the vectors of distances among the corresponding atoms
            in alternate conformers.  Every vector contains distances for the
            name-sorted atoms (use get_atom_names_sorted() method to obtain
            the list of atom names).  If more than two conformers are present,
            the possible pairs are iterated over.  For example, if three
            conformers are present, A/B/C, the method returns pairwise 
            distances in the following order - A-B, A-C, B-C.  Use 
            GetAltCodes() method to get the list of conformer codes. 
            Note that the return array is squeezed, i.e. when only two
            conformers are present, the return value is a simple vector. '''
        num = self.GetAltNum()
        rvecs = self.get_r_array()
        dvecs = []
        for i in range(num-1):
            for j in range (i+1, num):
                dvecs.append(sqrt(((rvecs[i]-rvecs[j])**2).sum(1)))
        return array(dvecs).squeeze()

    def get_b_ratios(self):
        ''' Returns the vectors of b-factor ratios for the corresponding atoms
            in alternate conformers.  Every vector contains ratios for the
            name-sorted atoms (use get_atom_names_sorted() method to obtain
            the list of atom names).  If more than two conformers are present,
            the possible pairs are iterated over.  For example, if three
            conformers are present, A/B/C, the method returns pairwise 
            ratios in the following order - A-B, A-C, B-C.  Use 
            GetAltCodes() method to get the list of conformer codes. 
            Note that the return array is squeezed, i.e. when only two
            conformers are present, the return value is a simple vector. '''
        num = self.GetAltNum()
        bvecs = self.get_b_vectors()
        dvecs = []
        for i in range(num-1):
            for j in range (i+1, num):
                dvecs.append(bvecs[i]/bvecs[j])
        return array(dvecs).squeeze()

    def esitimate_ac_occupancies(self, cutoff=1.0):
        ''' Estimate the occupancies of the alternate conformers.  The 
            algorithm assumes that the occupancy correction at the peak is
            dure to broadening of the electron density due to increased
            B-factor.  If the occupancies during refinement were not equal,
            then the returned values should be treated as relative correction
            factors, i.e. they should be applied after dividing by the number
            of the conformers. '''
        bvecs = self.get_b_ratios()
        if self.get_atom_number() == 1:
            fracs = array([1.0, sqrt(bvecs)])
        else:
            dvecs = self.get_alt_distances()
            num = self.GetAltNum() - 1
            ind = dvecs>cutoff
            if num > 1:
                ind = ind.all(0)
                fracs = array([1.0]+(sqrt(bvecs[:num].compress(ind,1))).mean(1).tolist())
            else:
                fracs=array([1.0,sqrt(bvecs[ind]).mean()])
        return fracs/sum(fracs)

    def Baverage(self):
        b, ball, oall = [], [], []
        of = 1.0 / float(len(self.altatoms))
        for atom in self.origatoms:
            ball.append(atom.GetB())
            oall.append(atom.GetOccupancy())
        bmean = (array(ball)*array(oall)).sum()/sum(oall)
        for ac in self.altatoms:
            bs, os = [], []
            for atom in ac.values():
                bs.append(atom.GetB())
                if atom.HasAltConf():
                    os.append(atom.GetOccupancy())
                else:
                    os.append(of*atom.GetOccupancy())
            b.append((array(bs)*array(os)).sum()/sum(os))
        for atom in self.origatoms:
            if atom.HasAltConf():
                atom.SetB(b[self.altcodes.index(atom.GetAltLoc())])
            else:
                atom.SetB(bmean)

    def GetGroupedBs(self, mode='lin'):
        ''' Returns average B-factor overall and for backbone and non-backbone 
            atoms. Default mode is linear average, but 'rms' will return root 
            mean square, which can be useful if the B-factor column is used to
            carry some other parameter that needs to be averaged this way.
            Averaging is occupancy-weighted.'''
        b = array([x.GetB() for x in self.origatoms])
        o = array([x.GetOccupancy() for x in self.origatoms])
        f = array([x.IsBackbone() for x in self.origatoms]).astype(int)
        n_all, n_bb, n_sc = sum(o), sum(f*o), sum((1-f)*o)
        if mode == 'lin':
            return [sum(b*o)/n_all if n_all else float('nan'),
                    sum(b*f*o)/n_bb if n_bb else float('nan'),
                    sum(b*(1-f)*o)/n_sc if n_sc else float('nan')]
        elif mode == 'rms':
            return [sqrt(sum(o*b**2)/n_all) if n_all else float('nan'),
                    sqrt(sum(f*o*b**2)/n_bb) if n_bb else float('nan'),
                    sqrt(sum((1-f)*o*b**2)/n_sc) if n_sc else float('nan')]
        else:
            raise KeyError("Incorrect b-factor averaging mode '"+mode+"'")

    def GetAtomsByElement(self, element):
        atoms = []
        for atom in self.atoms:
            if atom.GetElement() == element:
                atoms.append(atom)
        return atoms

    def GetNamesByElement(self, element):
        names = []
        for atom in self.GetAtoms():
            if atom.GetElement() == element:
                names.append(atom.GetName())
        return names

    def GetCoM(self):
        a, b = zeros(3), 0.0
        for atom in self.origatoms:
            m, o = pdbnames.GetMass(atom.GetElement().upper()), atom.GetOccupancy()
            a += atom.GetR() * m * o
            b += m * o
        return a/b

    def NonHydrogens(self):
        ''' Returns the list of non-hydrogen atoms. '''
        names = []
        for atom in self.GetAtoms():
            if atom.GetElement() != 'H':
                names.append(atom.GetName())
        return names

    def Polars(self):
        ''' Return the list of polar atoms.  Currently, this means everything but
            hydrogens and carbons.'''
        names = []
        for atom in self.GetAtoms():
            if atom.GetElement() not in 'CH':
                names.append(atom.GetName())
        return names

    def Ri(self, name, names=None):
        ''' Returns the dictionary of distances from an atom to other atoms in the
            residue.  If the list of names is not supplied, all atoms in the
            residue will be included in the list.'''
        if not names:
            names = self.atoms.keys()
        R, atomi = {}, self.GetAtom(name)
        if atomi:
            for n in names:
                atomj = self.GetAtom(n)
                if atomj:
                    R[n] = distance(atomi, atomj)
            return R
        return None

    def Rij(self, names1=None, names2=None):
        ''' Returns the dictionary of dictionaries with interatomic distances.
            Two list of names are needed, if not supplied, will include all the atoms.
            No checking is done for redundancy, which normally should not be a
            problem since single residue does not contain many atoms.  For massive
            calculations, use pbdmolecule.IndexDistances() method instead. '''
        if not names1:
            names1 = list(self.atoms.keys())
        if not names2:
            names2 = list(self.atoms.keys())
        R = {}
        for name in names1:
            Ri = self.Ri(name,names2)
            if Ri:
                R[name] = Ri
        return R

    def CovalentBonds(self, names1=None, names2=None, cutoff=2.0):
        ''' Returns the dictionary of dictionaries with covalent bonds formed between
            atoms.  It is not currently too smart since it simply finds atom pairs
            closer than cutoff distance, with no respect for atom types.  Thus, it
            will by default miss disulfide bonds and if residue contains hydrogens
            there will be false positives. '''
        rij, cb = self.Rij(names1, names2), {}
        for name1 in rij:
            bonds = {}
            for name2 in rij[name1]:
                if rij[name1][name2] < cutoff:
                    bonds[name2]=rij[name1][name2]
            if name1 in bonds:
                del bonds[name1]
            if bonds:
                cb[name1]=bonds
        return cb

    def FindCarboxylates(self):
        ''' Identifies carboxylates in the residue.  The following algorithm is used:
            - find all the carbons making two covalent bonds to different oxygens
            - check if oxygens are "terminal" '''
        carbons = self.GetNamesByElement('C')
        oxygens = self.GetNamesByElement('O')
        heavies = self.NonHydrogens()
        cb = self.CovalentBonds(carbons,oxygens)
        ob = self.CovalentBonds(oxygens,heavies)
        carboxylates = []
        for key in cb:
            if len(cb[key]) == 2:
                oxys = list(cb[key].keys())
                if len(ob[oxys[0]]) == 1 and len(ob[oxys[1]]) == 1:
                    carboxylates.append((key, oxys[0], oxys[1]))
        return carboxylates

    def Flip(self):
        ''' Flips the side chain around in (quasi)symmetrical residues.  Only
        operates on Asp, Glu, Phe, His, Asn, Gln, Tyr. '''
        if self.name == 'ASP':
            od1 = self.GetAtomR('OD1')
            od2 = self.GetAtomR('OD2')
            self.SetAtomR('OD1', od2)
            self.SetAtomR('OD2', od1)
        elif self.name == 'GLU':
            oe1 = self.GetAtomR('OE1')
            oe2 = self.GetAtomR('OE2')
            self.SetAtomR('OE1', oe2)
            self.SetAtomR('OE2', oe1)
        elif self.name == 'ASN':
            od1 = self.GetAtomR('OD1')
            nd2 = self.GetAtomR('ND2')
            cg = self.GetAtomR('CG')
            self.SetAtomR('OD1', 0.927*nd2 + 0.073*cg)
            self.SetAtomR('ND2', 1.079*od1 - 0.079*cg)
        elif self.name == 'GLN':
            oe1 = self.GetAtomR('OE1')
            ne2 = self.GetAtomR('NE2')
            cd = self.GetAtomR('CD')
            self.SetAtomR('OE1', 0.927*ne2 + 0.073*cd)
            self.SetAtomR('NE2', 1.079*oe1 - 0.079*cd)
        elif self.name == 'HIS':
            # this will produce slightly distorted geometry
            cd2 = self.GetAtomR('CD2')
            nd1 = self.GetAtomR('ND1')
            ce1 = self.GetAtomR('CE1')
            ne2 = self.GetAtomR('NE2')
            self.SetAtomR('CD2', nd1)
            self.SetAtomR('ND1', cd2)
            self.SetAtomR('CE1', ne2)
            self.SetAtomR('NE2', ce1)
    def __pairlist_(self, N):
        x = [[(i,j) for j in range(i+1,N)] for i in range(N-1)]
        return array([item for sublist in x for item in sublist])
    def BondsAnglesTorsions(self, Rcutoff=2.0, printout=False):
        Nats = len(self.origatoms)
        ps = self.__pairlist_(Nats)
        ds = array([self.distance(x[0],x[1]) for x in ps])
        bonds = ps[ds<Rcutoff]
        Nbs = len(bonds)
        ps = self.__pairlist_(Nbs)
        angles = array([array(x)[array([x.count(y) for y in x]).argsort()[array([0,2,1])]].tolist() for x in [x for x in [bonds[x[0]].tolist()+bonds[x[1]].tolist() for x in ps] if len(set(x))==3]])
        Nas = len(angles)
        ps = self.__pairlist_(Nas)
        torsions = array([x[:4] if x.tolist().count(x[-1])==2 else x[:3].tolist()+[x[-1]] for x in array([x[:3][::-1].tolist()+x[3:].tolist() if x.tolist().count(x[0])==2 else x for x in [x for x in array([x for x in [angles[x[0]].tolist()+angles[x[1]].tolist() for x in ps] if len(set(x))==4]) if x[1]!=x[4]]])])
        impropers = array(list(set([(x[1],x[0],x[2],x[3]) for x in [[x[1]]+sorted(set(x)-set([x[1]])) for x in [x for x in [x for x in [angles[x[0]].tolist()+angles[x[1]].tolist() for x in ps] if len(set(x))==4] if x[1]==x[4]]]])))
        if printout:
            bonds = ['%5s --- %5s : %10.3f' % (self.origatoms[x[0]].name(), self.origatoms[x[1]].name(), self.distance(x[0],x[1])) for x in bonds]
            angles = ['%5s - %5s - %5s : %10.2f' % (self.origatoms[x[0]].name(), self.origatoms[x[1]].name(), self.origatoms[x[2]].name(), self.angle(x[0],x[1],x[2])) for x in angles]
            torsions = ['%5s - %5s - %5s - %5s : %10.2f' % (self.origatoms[x[0]].name(), self.origatoms[x[1]].name(), self.origatoms[x[2]].name(), self.origatoms[x[3]].name(), self.torsion(x[0],x[1],x[2],x[3])) for x in torsions]
            impropers = ['%5s - %5s - %5s - %5s : %10.2f' % (self.origatoms[x[0]].name(), self.origatoms[x[1]].name(), self.origatoms[x[2]].name(), self.origatoms[x[3]].name(), self.torsion(x[0],x[1],x[2],x[3])) for x in impropers]
        return bonds, angles, torsions, impropers
        
class RemarkParser:

    def __init__(self, source):
        self.remarks = {}
        for line in source:
            if re.match('REMARK', line) != None:
                try:
                    (key,value)=line[6:].split(':')
                    if value.split():
                        try:
                            if ';' in value:
                                fv = float(value.split(';')[1])
                            else:
                                fv = float(value.split()[0])
                            if key.strip() in REFMAC_REMARK_FLOATS:
                                self.remarks[REFMAC_REMARK_FLOATS[key.strip()]] = fv
                            elif key.strip() in PHENIX_REMARK_FLOATS:
                                self.remarks[PHENIX_REMARK_FLOATS[key.strip()]] = fv
                        except KeyError:
                            pass
                except ValueError:
                    pass

    def GetKeys(self):
        return list(self.remarks.keys())

    def GetValue(self, key):
        try:
            return self.remarks[key]
        except KeyError:
            return None

    def GetRs(self):
        return (self.GetValue('R'), self.GetValue('RWORK'),self.GetValue('RFREE'))

    def GetR(self):
        return self.GetValue('R')

    def GetRwork(self):
        return self.GetValue('RWORK')

    def GetRfree(self):
        return self.GetValue('RFREE')

    def GetESUs(self):
        return (self.GetValue('ESU_R'), self.GetValue('ESU_RFREE'), self.GetValue('ESU_ML'), self.GetValue('ESU_B_ML'))

    def GetResolution(self):
        return self.GetValue('RESOLUTION')

    def GetNref(self):
        return self.GetValue('NREF')

    def GetNatm(self):
        return self.GetValue('NATOMS')

REFMAC_REMARK_FLOATS = {
'3   R VALUE     (WORKING + TEST SET)'                      :   'R',
'3   R VALUE            (WORKING SET)'                      :   'RWORK',
'3   FREE R VALUE'                                          :   'RFREE',
'3   ESU BASED ON R VALUE                            (A)'   :   'ESU_R',
'3   ESU BASED ON FREE R VALUE                       (A)'   :   'ESU_RFREE',
'3   ESU BASED ON MAXIMUM LIKELIHOOD                 (A)'   :   'ESU_ML',
'3   ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD (A**2)'   :   'ESU_B_ML',
'3   BOND LENGTHS REFINED ATOMS        (A)'                 :   'RMSD_BONDS',
'3   RESOLUTION RANGE HIGH (ANGSTROMS)'                     :   'RESOLUTION',
'3   NUMBER OF REFLECTIONS'                                 :   'NREF',
'3   ALL ATOMS'                                             :   'NATOMS'
}

PHENIX_REMARK_FLOATS = {
'3   COORDINATE ERROR (MAXIMUM-LIKELIHOOD BASED)'           :   'ESU_ML',
'3   BOND'                                                  :   'RMSD_BONDS',
'3    ALL'                                                  :   'NATOMS'
}

class TLSparser:

    def __init__(self, source):
        self.groups = []
        tlsFlag, tensorFlag = True, 0
        for line in source:
            if tlsFlag:
                if line[:24] == 'REMARK   3   TLS GROUP :':
                    tlsFlag = False
                    self.groups.append(TLSgroup())
            else:
                if tensorFlag:
                    self.groups[-1].SetTLS(line.split()[2][:3], float(line.split()[3]))
                    self.groups[-1].SetTLS(line.split()[4][:3], float(line.split()[5]))
                    if line.split()[2][0] == 'S':
                        self.groups[-1].SetTLS(line.split()[6][:3], float(line.split()[7]))
                    tensorFlag -= 1
                else:
                    if line[:29] == 'REMARK   3    RESIDUE RANGE :':
                        self.groups[-1].AddRange(line[29:].strip())
                    elif line[:35] == 'REMARK   3    ORIGIN FOR THE GROUP ':
                        (x,y,z) = split_nequal(line, 39, 3, 9)
                        self.groups[-1].SetOrigin((float(x),float(y),float(z)))
                    elif re.match('REMARK *3 *[TLS] TENSOR',line):
                        tensorFlag = 3
                    elif line[:24] == 'REMARK   3   TLS GROUP :':
                        self.groups.append(TLSgroup())

    def GetGroupNumber(self):
        return len(self.groups)

    def GetGroups(self):
        return self.groups

    def GetGroup(self,i):
        return self.groups[i]

class TLSgroup:

    def __init__(self, ranges=None, origin=None, T=[[0.0]*3]*3, L=[[0.0]*3]*3, S=[[0.0]*3]*3):
        if ranges == None:
            self.ranges = []
        else:
            self.ranges = ranges
        self.origin = origin
        self.TLS = {'T':T, 'L':L, 'S':S}

    def AddRange(self, line):
        (chain_start,residue_start,chain_end,residue_end) = line.split()
        self.ranges.append("'"+chain_start[0]+residue_start.rjust(4)+"' '"+chain_end[0]+residue_end.rjust(4) +"' ")

    def SetOrigin(self, xyz):
        self.origin = xyz

    def SetTLS(self,what,value):
        (matr, si, sj) = list(what)
        i, j = int(si)-1, int(sj)-1
        self.TLS[matr][i][j] = value
        if i!=j and matr!='S':
            self.TLS[matr][j][i] = value

    def GetOrigin(self):
        return self.origin

    def GetRanges(self):
        return self.ranges

    def GetRangeNumber(self):
        return len(self.ranges)

    def GetRange(self, i):
        return self.ranges[i]

    def GetTLS(self,matr):
        return self.TLS[matr]

class NCSparser:

    def __init__(self, source):
        self.groups = []
        ncsFlag, spanFlag = False, False
        for line in source:
            if ncsFlag:
                if spanFlag:
                    if line[:70] == 'REMARK   3                   GROUP CHAIN        COUNT   RMS     WEIGHT':
                        spanFlag, ncsFlag = False, False
                    else:
                        items = line.split()[2:]
                        span, chain, res1, res2 = int(items[0]), items[1], items[2], items[4]
                        if items[5] == 'NULL':
                            ncode = 1
                        else:
                            ncode = int(items[5])
                        self.groups[-1].AddChainSpan(span, chain, res1, res2, ncode)
                else:
                    if line[:47] == 'REMARK   3     CHAIN NAMES                    :':
                        self.groups[-1].SetChains(line.split(':')[1].split())
                    elif line[:47] == 'REMARK   3     NUMBER OF COMPONENTS NCS GROUP :':
                        Nsp = int(line.split(':')[1])
                    elif line[:59] == 'REMARK   3       COMPONENT C  SSSEQI  TO  C   SSSEQI   CODE':
                        spanFlag = True
                self.__groupCatcher(line)
            else:
                ncsFlag = self.__groupCatcher(line)

    def __groupCatcher(self, line):
        if line[:44] == 'REMARK   3  NCS GROUP NUMBER               :':
            self.groups.append(NCSgroup())
            return True

    def GetGroupNumber(self):
        return len(self.groups)

    def GetGroups(self):
        return self.groups

    def GetGroup(self,i):
        return self.groups[i]

class NCSgroup:

    def __init__(self, chains=None, spans=None):
        if chains:
            self.chains = chains
        else:
            self.chains = []
        if spans:
            self.spans = spans
        else:
            self.spans = {}

    def GetChains(self):
        return self.chains

    def SetChains(self, chains):
        self.chains = chains

    def AppendChain(self, chain):
        self.chains.append(chain)

    def DeleteChain(self, chain):
        if chain in self.chains:
            self.chains.remove(chain)

    def GetChainNumber(self):
        return len(self.chains)

    def GetSpanNumber(self):
        return len(self.spans)

    def AddChainSpan(self, span, chain, res1, res2, ncode):
        if span in self.spans:
            if res1 == self.spans[span][1] and res2 == self.spans[span][2]:
                self.spans[span][0].append(chain)
        else:
            self.spans[span] = ([chain], res1, res2, ncode)

    def GetCommand(self, refware='refmac'):
        if refware == 'refmac':
            Nch = self.GetChainNumber()
            if Nch<2 or Nch>40:
                return None
            line = 'ncsr nchains %d chains ' % len(self.chains)
            for chain in self.chains:
                line += chain + ' '
            line += 'nspans ' + str(len(self.spans)) + ' '
            for span in sorted(self.spans):
                line += self.spans[span][1] + ' ' + self.spans[span][2] + ' ' + str(self.spans[span][3]) + ' '
            if 'NULL' in line:
                return None
            return line
        return None

class fakeatom:

    def __init__(self, xyz):
        self.xyz = array(xyz)

def distance2(r1, r2):
    return math.sqrt(((r1-r2)**2).sum())

def r(atom1, atom2):
    return atom2.xyz-atom1.xyz

def torsion(i,j,k,l):
    b1 = r(i,j)
    b2 = r(j,k)
    b3 = r(k,l)
    b2xb1 = cross(b2,b1)
    x = (b3*cross(b2,b2xb1)).sum()
    y = math.sqrt((b2**2).sum())*(b3*b2xb1).sum()
    return math.degrees(math.atan2(-y,x))

def angle(i,j,k):
    r1 = r(i,j)
    r2 = r(k,j)
    return math.degrees(arccos((r1*r2).sum()/math.sqrt((r1**2).sum()*(r2**2).sum())))

def angle2(i,j,k):
    r1 = i-j
    r2 = k-j
    return math.degrees(arccos((r1*r2).sum()/math.sqrt((r1**2).sum()*(r2**2).sum())))

def distance(i,j):
    return math.sqrt(((i.xyz-j.xyz)**2).sum())

def point2segment(a1, a2, b):
    p = a2 - a1
    t = ((b-a1)*p).sum() / (p**2).sum()
    if t<=0:
        return math.sqrt(((b-a1)**2).sum())
    elif t >=1:
        return math.sqrt(((b-a2)**2).sum())
    else:
        return math.sqrt(((b-a1-t*p)**2).sum())

def segment2segment(seg1, seg2):
    (a1,b1), (a2,b2) = seg1, seg2
    p1 = b1 - a1
    p2 = b2 - a2
    D = (p1**2).sum()*(p2**2).sum() - ((p1*p2).sum())**2
    if not D:
        return point2segment(a1, b1, a2)
    D1 = (p2**2).sum() * (p1*(a2-a1)).sum() + (p1*p2).sum() * (p2*(a1-a2)).sum()
    D2 = (p1**2).sum() * (p2*(a1-a2)).sum() + (p1*p2).sum() * (p1*(a2-a1)).sum()
    t1, t2 = D1/D, D2/D
    if t1<=0:
        return point2segment(a2,b2,a1)
    elif t1>=1:
        return point2segment(a2,b2,b1)
    elif t2<=0:
        return point2segment(a1,b1,a2)
    elif t2>=1:
        return point2segment(a1,b1,b2)
    else:
        return math.sqrt(((a1-a2+p1*t1-p2*t2)**2).sum())

def MakeResID(chid, resn, icode=' '):
    return chid + '%4s' % resn + icode

def MakeAtomID(name, chid, resn, altloc=' ', icode=' '):
    return name+altloc+chid + '%4s' % resn + icode

def SplitResID(resid):
    return (resid[0], int(resid[1:5]), resid[5])

def split_nequal(line, start=0, number=1, width=1):
    ''' Returns the list of substrings taken from the input line.  Each
        substring is width long, and number chunks are returned.  It is
        also possible to start from start position.  This method can be
        useful to, for instance, extract numbers that are formatted to 
        certain length without regard for the number of significant 
        digits (which could fail if plain split() string mthod is used). '''
    return [line[i:i+width] for i in range(start, start+number*width, width)]
