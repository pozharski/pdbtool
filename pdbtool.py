'''
Module for reading PDB-files.  
Includes pdbatom and pdbmolecule classes
'''

import gzip, urllib.request, os, random, math, sys, re, copy, logging, time

from . import pdbnames, SpaceGroups
from .helper import progressbar
from .rotate import transform_list
from tinertia import TInertia
from scipy.linalg import eigh
from scipy import   array, cos, sin, pi, radians, sqrt, dot, cross, \
                    randn, zeros, matrix, ones, floor, nonzero

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
    for line in source:
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

def cell_and_center(molecule, scale=2.0):
    ''' The method returns the CRYST1 line for the P1 cell that could hold
        the entire molecule and the copy of the molecule with shifted
        coordinates that place it in the center. '''
    r = molecule.GetCoordinateArray()
    gabarit = r.ptp(0)
    abc = scale * gabarit
    retmol = molecule.copy()
    retmol.shift(0.5 * scale * gabarit - r.mean(0))
    celline =  'CRYST1'
    celline += '%9.3f%9.3f%9.3f' % tuple(abc)
    celline += '  90.00  90.00  90.00 P 1                     \n'
    return (celline, retmol)

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
        return pdbnames.IsPolar(self.name())

    def IsProteinBackbone(self):
        ''' True if atom belongs to protein backbone, False otherwise. '''
        return pdbnames.Is3Amino(self.resName()) and pdbnames.MaybeBackbone(self.name())

    def IsProtein(self):
        ''' True if atom belongs to a protein, False otherwise. '''
        return pdbnames.Is3Amino(self.resName())

    def IsHetero(self):
        return pdbnames.IsHetero(self.resName)

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

    def backbone(self):
        return backbone(self)

    def is_multi_model(self):
        return self.modelN > 1

    def rjust_res_names(self):
        for atom in self.atoms:
            atom.rjust_res_name()

    def prime_uij(self, overwrite=False, listik=None):
        ''' Initializes the anisotropik ADP for the atoms in listik (default
            is all the atoms in the molecule).  The overwrite flag defines if
            the ANISOU record will be overwritten if already persent.'''
        for atomi in listik if listik else range(self.GetAtomNumber()):
            self.atoms[atomi].prime_uij(overwrite=overwrite)

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

    def SetBfactorValues(self, bvalues, listik=None):
        ''' Sets the atomic B-factors using the provided values mapped into the
            list of atoms. If bvalues is a single number, all the atoms
            in the list will have the same B-factor.  If listik is not
            supplied, all atoms will have their B-factors reset. '''
        if not listik:
            listik = range(self.GetAtomNumber())
        if type(bvalues) == float:
            for atomi in listik:
                self.atoms[atomi].SetB(bvalues)
        else:
            assert len(bvalues)==len(listik), 'Shape mismatch between selection and Bvalues vector.'
            for (i, atomi) in enumerate(listik):
                self.atoms[atomi].SetB(bvalues[i])

    def set_occupancies(self, values, listik=None):
        ''' Sets the atomic occupancies using the provided values mapped into 
            the list of atoms.  If values is a single number, all the atoms
            in the list will have the same occupancy.  If listik is not
            supplied, all atoms will have their occupancies reset. '''
        if not listik:
            listik = range(self.GetAtomNumber())
        if type(values) == float:
            for atomi in listik:
                self.atoms[atomi].SetOccupancy(values)
        else:
            assert len(values)==len(listik), 'Shape mismatch between selection and values vector.'
            for (i, atomi) in enumerate(listik):
                self.atoms[atomi].SetOccupancy(values[i])

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
        listik1 = self.ListResidues(hasAltConf=True, fProteinOnly=True)
        listik2 = other.ListResidues(hasAltConf=True, fProteinOnly=True)
        listik = list(set(listik1).union(listik2))
        residues1 = self.GetResiduesFromList(listik)
        residues2 = other.GetResiduesFromList(listik)
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

    def get_elements(self, listik=None):
        if listik is None:
            listik = self.atom_lister()
        elif type(listik) == str:
            listik = self.atom_lister(listik)
        return sorted(set([a.GetName() for a in self.GetListedAtoms(listik)]))

    def get_atom_types(self, listik=None):
        if listik is None:
            listik = self.atom_lister()
        elif type(listik) == str:
            listik = self.atom_lister(listik)
        return sorted(set([a.GetName() for a in self.GetListedAtoms(listik)]))

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
        return list(array(self.atoms)[listik])

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

    def GetBvector(self, listik=None):
        ''' Returns an array of B-factor values for a list of atoms.
            Defaults to all atoms in the molecule. '''
        if not listik:
            return array([x.GetB() for x in self.atoms])
        else:
            return array([x.GetB() for x in array(self.atoms)[listik]])

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

    def GetOccupancyVector(self, listik=None):
        ''' Returns an array of occupancy values for a list of atoms.
            Defaults to all atoms in the molecule. '''
        return array([a.GetOccupancy() for a in self.atom_getter(listik=listik)])

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
            'rat'               - specific residue/atom name combination.
                                    For example, ASPOD1 designates OD1
                                    atoms in aspartic acids.  Must provide 
                                    the list of combos to include as rats 
                                    parameter
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
            return [x for x in set(atoms).difference(kwargs['coreatoms']) if (((coreXYZ-x.GetR())**2).sum(1)<=r2cutoff).any()]
        elif whatlow == 'rat':
            return [x for x in atoms if x.rat in kwargs['rats']]

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
            'rat'               - specific residue/atom name combination.
                                    For example, ASPOD1 designates OD1
                                    atoms in aspartic acids.  Must provide 
                                    the list of combos to include as rats 
                                    parameter
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
            return [i for i in atoms if self.atoms[i].get_res_name() in kwargs['resnames']]
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
            coreXYZ = self.GetListCoordinateArray(kwargs['corelist'])
            return sorted([i for i in set(listik).difference(kwargs['corelist']) if (((coreXYZ-self.GetAtomR(i))**2).sum(1)<=r2cutoff).any()])
        elif whatlow == 'rat':
            return [i for i in listik if self.atoms[i].rat() in kwargs['rats']]

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

    def ListAltConfTypes(self, listik=False):
        '''
        Returns the list of alternate conformer labels present in the
        molecule.
        '''
        return sorted(set([self.atoms[i].GetAltLoc() for i in self.__ensure_listik_(listik) if self.atoms[i].HasAltConf()]))

    def ListChainSplit(self):
        ''' Returns the dictionary of list of atoms from individual chains.'''
        chainsplit = {}
        for (i,atom) in enumerate(self.atoms):
            chid = atom.GetChain()
            try:
                chainsplit[chid].append(i)
            except KeyError:
                chainsplit[chid] = [i]
        return chainsplit
        
    def ListBackbone(self, listik=False):
        return filter(lambda i : self.atoms[i].IsBackbone(), self.__ensure_listik_(listik))

    def ListSideChains(self, listik=False):
        return filter(lambda i : self.atoms[i].NotBackbone(), self.__ensure_listik_(listik))

    def ListCompleteResidues(self, source):
        residues = []
        for atomi in source:
            resid = self.atoms[atomi].GetResID()
            if resid not in residues:
                residues.append(resid)
        ind = []
        for (i,atom) in enumerate(self.atoms):
            if atom.GetResID() in residues:
                ind.append(i)
        return ind

# ---

    def ListResidues(self, doSort=True, byChains=False, hasAltConf=False, fProteinOnly=False, resnames=False):
        ''' Returns the list of residues in the molecule.  If doSort (default),
            the list will be sorted using resid_compare() method.  Otherwise, the 
            list will be in the same order in whic atoms appear in the pdb file.
            If byChains, the method will return the disctionary with chain IDs
            as keys and list of residues in the corresponding chains as values.
            If hasAltConf, only the residues that have alternate conformers will be
            included (otherwise a complete list is returned).
            If resnames list of residue names is provided, only the
            residues with matching names will be selected.'''
        residues = []
        if hasAltConf:
            source = self.ListAltConf()
        else:
            source = range(self.GetAtomNumber())
        if fProteinOnly:
            source = list(set(source).intersection(self.atom_lister('protein')))
        if resnames:
            source = filter(lambda i : self.atoms[i].get_res_name() in resnames,
                                                                        source)
        for atomi in source:
            resid = self.atoms[atomi].GetResID()
            if resid not in residues:
                residues.append(resid)
        if byChains:
            dicres = {}
            for resid in residues:
                chid = resid[0]
                if chid not in dicres:
                    dicres[chid] = [resid]
                else:
                    dicres[chid].append(resid)
            if doSort:
                for chid in dicres:
                    dicres[chid].sort()
            return dicres
        else:
            if doSort:
                residues.sort()
            return residues

    def ListProteinResidues(self):
        residues = []
        for atom in self.atoms:
            resid, resname = atom.GetResID(), atom.get_res_name()
            if pdbnames.Is3Amino(resname):
                if resid not in residues:
                    residues.append(resid)
        return residues

    def GetProteinResidueNumber(self):
        return len(self.ListProteinResidues())

    def ListHeteroResidues(self):
        ''' Returns the list of residues to be considered heteroatoms.  Currently,
            the method select residues other than twenty canonical amino acids and
            waters. '''
        return list(set([a.resid() for a in self.atom_getter('hetero')]))

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
        listik = []
        residues = self.get_residues()
        for resid in residues:
            if residues[resid].HasAltConf():
                listik.append(resid)
                if fProteinOnly:
                    if not residues[resid].IsAminoAcid():
                        listik.pop()
        return listik

# --- Methods returning average values

    def GetAverageBfactor(self, group='all', chids=None, selection=[], *args, **kwargs):
        ''' Returns the average B-factor of a group of atoms.  If chids is specified,
            selection is restricted to certain chains.  The type of selection is
            determined by the group parameter which may have the following values:
                'all'       - all atoms in the molecule
                'protein'   - only protein atoms
                'proteinbb' - only protein backbone atoms
                'proteinsc' - only protein side chain atoms
                'water'     - only water atoms (residue name HOH)
                'hetero'    - heteroatoms (e.g. not protein and not waters)
                'resname'   - listed residue names (provide resids as parameter)
                'list'      - index of atoms passed as selection
                'chains'    - individual chains (returns a dictionary) '''
        if group.lower()=='all':
            ind = self.atom_lister('chids', chids=chids)
        elif group.lower()=='protein':
            if chids is not None:
                ind = self.merge_listers(['chids','protein'], chids=chids)
            else:
                ind = self.atom_lister('protein')
        elif group.lower()=='proteinbb':
            ind = self.merge_listers(['proteinbackbone','chids'], chids=chids)
        elif group.lower()=='proteinsc':
            ind = self.merge_listers(['proteinsidechain','chids'], chids=chids)
        elif group.lower()=='water':
            ind = self.merge_listers(['water','chids'], chids=chids)
        elif group.lower()=='hetero':
            ind = self.merge_listers(['hetero','chids'], chids=chids)
        elif group.lower()=='resname':
            ind = self.atom_lister('resids', resids=resids)
            if chids is not None:
                ind = self.atom_lister('chids', ind, chids=chids)
        elif group.lower()=='list':
            ind = selection
        elif group.lower()=='chains':
            bvs = {}
            for (chid,chind) in self.ListChainSplit().items():
                bvs[chid] = self.GetListAverageB(chind)
            return bvs
        else:
            return float('nan')
        return self.GetListAverageB(ind)

    def GetListAverageB(self, listik):
        bvals = self.GetBvector(listik)
        ovals = self.GetOccupancyVector(listik)
        return sum(bvals*ovals)/sum(ovals)

# ---

    def get_listed_resids(self, listik, doSort=True):
        residues = []
        for atomi in listik:
            resid = self.atoms[atomi].GetResID()
            if resid not in residues:
                residues.append(resid)
        if doSort:
            residues.sort()
        return residues

    def get_listed_altconf(self, listik, doSort=True):
        altconf = []
        for atomi in listik:
            ac = self.atoms[atomi].GetAltLoc()
            if ac not in altconf:
                altconf.append(ac)
        if doSort:
            altconf.sort()
        return altconf

    def GetCrystLine(self):
        ''' Return the CRYST1 pdb record for the molecule. '''
        if self.cell == None:
            return ''
        else:
            return self.cell.GetLine()

    def GetSymmOps(self):
        if self.cell:
            return SpaceGroups.symops[self.cell.GetSG()]
        else:
            return None

    def GetSymmNum(self):
        if self.cell:
            return len(SpaceGroups.symops[self.cell.GetSG()])
        else:
            return 0

    def GetFractionalCoordinates(self, r):
        if self.cell:
            M = self.cell.GetMcf()
            return dot(M,r.T).T
        else:
            return None

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
        if i>0:
            if i < self.GetSymmNum():
                symop = self.GetSymmOps()[i]
                self.Fractionalize()
                for atom in self.atoms:
                    atom.ApplyOperator(symop)
                self.Cartesianize()

    def CreateSymate(self, i, cell_shift=array([0.0,0.0,0.0])):
        if i==0:
            return self.copy()
        elif i>0:
            if i < self.GetSymmNum():
                x = self.copy()
                x.ApplySymmetryOperator(i)
                if array(cell_shift).any():
                    x.CellShift(cell_shift)
                return x
        return None

    def GetListCoordinateArray(self, listik=None):
        if listik:
            xyz = []
            for i in listik:
                xyz.append(self.GetAtom(i).GetR())
            return array(xyz)
        else:
            return self.GetCoordinateArray()

    def GetCoordinateArray(self):
        ''' Return the Nx3 array of individual atom coordinates. '''
        return array(map(lambda x : x.GetR(), self.atoms))

    def SymateDistance(self, i, cell_shift=array([0.0,0.0,0.0])):
        symate = self.CreateSymate(i, cell_shift)
        if symate:
            return self.ShortestApproach(symate)
#            xyz0 = self.GetCoordinateArray()
#            xyz1 = symate.GetCoordinateArray()
#            N = self.GetAtomNumber()
#            return sqrt(array((matrix(xyz1.T[0]).T*ones(N)).T-matrix(xyz0.T[0]).T*ones(N))**2+array((matrix(xyz1.T[1]).T*ones(N)).T-matrix(xyz0.T[1]).T*ones(N))**2+array((matrix(xyz1.T[2]).T*ones(N)).T-matrix(xyz0.T[2]).T*ones(N))**2).min()
        return None

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

    def Distances2Bfactors(self, other, listik=None):
        self.SetBfactorValues(bvalues=self.ShortestDistances(other), listik=listik)

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
            x.writePDBwithCell(pdbfile)

    def WriteSymates(self, pdbfile, move2cell=True):
        self.__CellCheck_()
        for i in range(self.GetSymmNum()):
            symate = self.CreateSymate(i)
            if move2cell:
                symate.MoveCenter2UnitCell()
            symate.writePDBwithCell(pdbfile+str(i)+os.extsep+'pdb')

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
        return math.sqrt(((self.atoms[i].xyz-self.atoms[j].xyz)**2).sum())
#        return math.sqrt((self.atoms[i].x-self.atoms[j].x)**2 + (self.atoms[i].y-self.atoms[j].y)**2 + (self.atoms[i].z-self.atoms[j].z)**2)

    def select(self,name=None,altloc=None,resn=None,resid=None):
        ind = []
        for (i,atom) in enumerate(self.atoms):
            if atom.test(name,altloc,resn,resid):
                ind.append(i)
        return ind

    def idselect(self, atid):
        return self.select(name=atid.GetAtomName(), resid=atid.ResID(), altloc=atid.GetAltLoc())[0]

    def angle(self,i,j,k):
        r1 = self.r(j,i)
        r2 = self.r(j,k)
        return math.degrees(math.acos((r1*r2).sum()/math.sqrt((r1**2).sum()*(r2**2).sum())))
#        return math.degrees(math.acos(self.dotprod(r1,r2)/(self.vecnorm(r1)*self.vecnorm(r2))))

    def r(self,i,j):
        return self.atoms[j].xyz - self.atoms[i].xyz
#        return [self.atoms[j].x-self.atoms[i].x, self.atoms[j].y-self.atoms[i].y, self.atoms[j].z-self.atoms[i].z]

    def torsion(self,i,j,k,l):
        b1 = self.r(i,j)
        b2 = self.r(j,k)
        b3 = self.r(k,l)
        b2xb1 = cross(b2,b1)
        x = (b3*cross(b2,b2xb1)).sum()
        y = math.sqrt((b2**2).sum())*(b3*b2xb1).sum()
        return math.degrees(math.atan2(-y,x))

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
        bonds = {}
        for (atomj,atom) in enumerate(self.atoms):
            if not (polaronly and atom.GetElement()=='C'):
                D = self.distance(atomi,atomj)
                if D < rmax:
                    bonds[atomj] = D
        bonds.pop(atomi,0)
        return bonds

    def same_residue(self, atomi, atomj):
        ''' Returns True if atoms belong to the same residue (including
            alternate conformers), False otherwise. '''
        return self.atoms[atomi].same_residue(self.atoms[atomj])

    def first_shell(self, atomi, rmax=3.2, polaronly=True):
        ''' This is similar to environment() method, but excludes the intra-residue
            contacts (e.g. covalent bonds). '''
        bonds = {}
        for (atomj,atom) in enumerate(self.atoms):
            if not (polaronly and atom.GetElement()=='C'):
                D = self.distance(atomi,atomj)
                if D < rmax:
                    if not self.same_residue(atomi,atomj):
                        bonds[atomj] = D
        bonds.pop(atomi,0)
        return bonds

    def GetCoM(self, listik=[]):
        if not listik:
            listik = range(self.GetAtomNumber())
        a, b = zeros(3), 0.0
        for i in listik:
            m, o = pdbnames.GetMass(self.atoms[i].GetElement().upper()), self.atoms[i].GetOccupancy()
            a += self.atoms[i].GetR() * m * o
            b += m * o
        return a/b

    def Rgyration(self, listik=None):
        if listik is None:
            listik = range(self.GetAtomNumber())
        mo = array([pdbnames.GetMass(x.element().strip())*x.GetOccupancy() for x in array(self.atoms)[listik]]).T
        r = array([x.GetR() for x in array(self.atoms)[listik]]).T
        Rcenter = (mo*r).sum(1)/mo.sum()
        return sqrt(sum(mo*((r.T-Rcenter)**2).sum(1))/sum(mo))

    def GetInertiaTensor(self, listik=[]):
        if not listik:
            listik = range(self.GetAtomNumber())
        atoms = self.GetListedAtoms(listik)
        m = array([pdbnames.GetMass(x.element().strip()) for x in atoms])
        o = array([x.GetOccupancy() for x in atoms])
        x, y, z = array([x.GetR() for x in atoms]).T
        return TInertia(m*o, x, y, z)

    def SelectionNeighborhood(self, selection, rmax=4.0, listik=None):
        ''' Returns the list of atoms that are within rmax from an atom in the 
            selection. Search can be narrowed to atoms in listik, which 
            defaults to all atoms.'''
        if listik is None:
            listik = range(self.GetAtomNumber())
        if len(listik) == 0:
            return array([])
        listik = list(set(listik).difference(selection))
        xyz0 = self.GetListCoordinateArray(selection)
        xyz1 = self.GetListCoordinateArray(listik)
        N0, N1 = len(xyz0), len(xyz1)
        return array(listik)[nonzero(sqrt(array((matrix(xyz1.T[0]).T*ones(N0)).T-matrix(xyz0.T[0]).T*ones(N1))**2 + array((matrix(xyz1.T[1]).T*ones(N0)).T-matrix(xyz0.T[1]).T*ones(N1))**2 + array((matrix(xyz1.T[2]).T*ones(N0)).T-matrix(xyz0.T[2]).T*ones(N1))**2).min(0)<rmax)]

    def NeighborAtoms(self, i, rmax=4.0, listik=None):
        ''' Return the list of atoms within rmax from atom i (could be
            either index or pdbatom object. Search can be 
            narrowed to atoms in listik, which defaults to all atoms.'''
        nelist = []
        if not listik:
            listik = range(self.GetAtomNumber())
        if type(i) is int:
            for j in listik:
                if self.distance(i,j) < rmax:
                    nelist.append(j)
        else:
            for j in listik:
                if distance(i,self.atoms[j]) < rmax:
                    nelist.append(j)
        return nelist

    def NeighborResidues(self, resid, rmax=4.0):
        residlist = []
        for i in self.atom_lister('resid', resid=resid):
            for j in self.NeighborAtoms(i=i, rmax=rmax):
                this_resid = self.atoms[j].GetResID()
                if this_resid not in residlist:
                    residlist.append(this_resid)
        return sorted(residlist)

    def WriteNeighborResidues(self, pdbFile, resid, rmax=4.0, self_exclude=False, header=None):
        residlist = self.NeighborResidues(resid=resid, rmax=rmax)
        if self_exclude:
            residlist.remove(resid)
        self.WriteResidueList(pdbFile=pdbFile, residlist=residlist, header=header)

    def WriteResidueList(self, pdbFile, residlist, header=None):
        fout = open(pdbFile, 'w')
        if header == 'cell':
            self.__headwrite_(fout, self.GetCrystLine())
        else:
            self.__headwrite_(fout, header)
        for atom in self.atoms:
            if atom.GetResID() in residlist:
                fout.write(atom.GetAtomRecord())
        fout.close()

    def WriteAtomList(self, pdbFile, atomlist, header=None):
        fout = open(pdbFile, 'w')
        if header == 'cell':
            self.__headwrite_(fout, self.GetCrystLine())
        else:
            self.__headwrite_(fout, header)
        for atomi in sorted(atomlist):
            fout.write(self.GetAtom(atomi).GetAtomRecord())
        fout.write('END   \n')
        fout.close()

    def GetResTitle(self, i):
        return self.atoms[i].GetResTitle()

    def get_residues(self, resname=False):
        ''' Returns the dictionary of residues.
            If resname is specified, only the residues of that type are 
            included, e.g. get_residues(resname='HIS') will return the
            dictionary of histidines. '''
        residues = {}
        for atom in self.atoms:
            if resname:
                if atom.get_res_name()!=resname:
                    continue
            resid = atom.GetResID()
            if resid not in residues.keys():
                residues[resid] = [atom]
            else:
                residues[resid].append(atom)
        for resid in residues.keys():
            residues[resid] = pdbresidue(residues[resid])
        return residues

    def GetResiduesFromList(self, listik):
        residues = {}
        for atom in self.atoms:
            resid = atom.GetResID()
            if resid in listik:
                if resid not in residues.keys():
                    residues[resid] = [atom]
                else:
                    residues[resid].append(atom)
        for resid in residues.keys():
            residues[resid] = pdbresidue(residues[resid])
        return residues

    def GetProteinResidues(self):
        residues = {}
        for atom in self.atoms:
            resid, resname = atom.GetResID(), atom.get_res_name()
            if pdbnames.Is3Amino(resname):
                if resid not in residues.keys():
                    residues[resid] = [atom]
                else:
                    residues[resid].append(atom)
        for resid in residues.keys():
            residues[resid] = pdbresidue(residues[resid])
        return residues

    def GetResidue(self, resid):
        ''' Returns residue class for the group of atoms with resid.  If resid is
            malformed (e.g. such residue is not found in the molecule), returns None. '''
        atoms = []
        for atom in self.atoms:
            if atom.GetResID() == resid:
                atoms.append(atom)
        if atoms:
            return pdbresidue(atoms)
        else:
            return None

    def GetResidueNames(self, listik=False):
        return sorted(set(map(lambda i : self.atoms[i].get_res_name(), self.__ensure_listik_(listik))))

    def GetChains(self):
        ''' Returns the dictionary of chains with chain IDs as keys and
            number of atoms in each chain as values. '''
        chains = {}
        for atom in self.atoms:
            chid = atom.GetChain()
            try:
                chains[chid] += 1
            except KeyError:
                chains[chid] = 1
        return chains

    def __headwrite_(self, fout, header):
        tyhe = type(header)
        if tyhe is str:
            fout.write(header)
        elif tyhe is list or tyhe is tuple:
            for line in header:
                fout.write(line)
            if not header:
                fout.write(self.GetCrystLine())

    def writePDB(self, pdbFile, header=None, mode='w'):
        ''' Write the PDB file.  pdbFile is the path to the output file, header may contain
            single line, list or tuple  of lines which will be placed at the beginning of the file.''' 
        fout = open(pdbFile, mode)
        if header == 'cell':
            self.__headwrite_(fout, self.GetCrystLine())
        else:
            self.__headwrite_(fout, header)
        for atom in self.atoms:
            fout.write(atom.GetAtomRecord())
            fout.write(atom.GetAnisouRecord())
        fout.close()

    def writePDBwithCell(self, pdbFile):
        ''' Write the PDB file and include the CRYST1 record of the unit cell parameters.
            pdbFile is the path to the output file.'''
        self.writePDB(pdbFile, self.GetCrystLine())

    def writePDBchains(self, pdbFile, chains, header=None):
        ''' Write the PDB file that only include the specified chains.
            pdbFile is the path to the output file. chains could be either
            a list or a string with selection chain IDs.'''
        fout = open(pdbFile, 'w')
        if header == 'cell':
            self.__headwrite_(fout, self.GetCrystLine())
        else:
            self.__headwrite_(fout, header)
        for atom in self.atoms:
            if atom.GetChain() in chains:
                fout.write(atom.GetAtomRecord())
        fout.close()

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

    def __range_checker_(self, value, limits):
        return value >= limits[0] and value <= limits[1]

    def extract_chains(self, chains):
        ''' Returns the copy of the molecule that only contains atoms from
            the supplied list of chains (the latter could be either list or
            string of symbols. '''
        extracted_atoms = []
        for atom in self.atoms:
            if atom.GetChain() in chains:
                extracted_atoms.append(atom.copy())
        return pdbmolecule(atoms=extracted_atoms, cell=self.cell)

    def rename_chains(self, chains):
        ''' Rename chain IDs.  Supply the dictionary for renaming. '''
        for atom in self.atoms:
            chid = atom.GetChain()
            if chid in chains:
                atom.SetChain(chains[chid])

    def noise(self, xnoise=0.1, bnoise=0.1, occnoise=0.0):
        for atom in self.atoms:
            atom.xyz = atom.xyz+xnoise*randn(3)
#            (atom.x, atom.y, atom.z) = atom.xyz #remove
            atom.SetB(math.fabs(atom.GetB()*random.gauss(1,bnoise)))
            if atom.IsAnisotropic():
                atom.SetUijValues(tuple((array(atom.GetUijValues())*abs(1.0+bnoise*randn(6))).astype(int)))
            atomocc = float(atom.GetOccupancy())
            if atomocc<1.0 and atomocc>0.0:
                newocc = random.gauss(atomocc,occnoise)
                while newocc>1.0 or newocc<0.0:
                    newocc = random.gauss(atomocc,occnoise)
                atom.SetOccupancy(newocc)

    def copy(self):
        return copy.deepcopy(self)

    def PhiPsiList(self):
        return backbone(self).PhiPsiList()
        reslist = self.ListProteinResidues()
        residues = self.GetProteinResidues()
        phi, psi = {}, {}
        for (i,resid) in enumerate(reslist):
            try:
                if i:
                    previd = reslist[i-1]
                    if distance(residues[resid].GetAtom('N'),residues[previd].GetAtom('C')) < 2.5:
                        phi[resid] = torsion(residues[previd].GetAtom('C'), residues[resid].GetAtom('N'), residues[resid].GetAtom('CA'), residues[resid].GetAtom('C'))
                nextid = reslist[i+1]
                psi[resid] = torsion(residues[resid].GetAtom('N'), residues[resid].GetAtom('CA'), residues[resid].GetAtom('C'), residues[nextid].GetAtom('N'))
            except (KeyError, IndexError):
                logging.info('No peptide bond atoms for ' + resid)
            finally:
                pass
        return (phi, psi)

    def BackboneTorsions(self):
        reslist = self.ListProteinResidues()
        residues = self.GetProteinResidues()
        phi, psi, omega = {}, {}, {}
        for (i,resid) in enumerate(reslist):
            try:
                if i:
                    previd = reslist[i-1]
                    if distance(residues[resid].GetAtom('N'),residues[previd].GetAtom('C')) < 2.5:
                        phi[resid] = torsion(residues[previd].GetAtom('C'), residues[resid].GetAtom('N'), residues[resid].GetAtom('CA'), residues[resid].GetAtom('C'))
                nextid = reslist[i+1]
                psi[resid] = torsion(residues[resid].GetAtom('N'), residues[resid].GetAtom('CA'), residues[resid].GetAtom('C'), residues[nextid].GetAtom('N'))
                if distance(residues[resid].GetAtom('C'),residues[nextid].GetAtom('N')) < 2.5:
                    omega[resid] = torsion(residues[resid].GetAtom('O'), residues[resid].GetAtom('C'), residues[nextid].GetAtom('N'), residues[nextid].GetAtom('CA'))
            except (KeyError, IndexError):
                logging.info('No peptide bond atoms for ' + resid)
            finally:
                pass
        return (phi, psi, omega)

    def GetWaters(self):
        ''' Return the dictionary of waters as residues.'''
        waters = {}
        residues = self.get_residues()
        for resid in residues:
            if residues[resid].get_res_name() == 'HOH':
                waters[resid] = residues[resid]
        return waters

    def W2Wdistance(self, water1, water2):
        return distance(water1.GetAtom('O'),water2.GetAtom('O'))

    def PolarContacts(self, rmax=3.2, include_interchain=True):
        ''' 
            Return the list of polar contacts. Excludes atoms from different 
            alternate conformers and same residues.  May exclude interchain
            contacts (say if you are not intersted in crystal contacts of a
            confirmed monomer) if include_interchain is set to false.
        '''
        indp = self.atom_lister('polar')
        contacts = []
        for (i,atomi) in enumerate(indp):
            resid = self.GetAtomResID(atomi)
            d = self.IndexDistances([atomi],indp[(i+1):])
            for (j,atomj) in enumerate(indp[(i+1):]):
                if d[0][j] < rmax:
                    if self.GetAtomResID(atomj) != resid:
                        if self.GetAtom(atomi).same_alt(self.GetAtom(atomj)):
                            if include_interchain or self.GetAtom(atomi).same_chain(self.GetAtom(atomj)):
                                contacts.append((atomi, atomj, d[0][j]))
        return contacts

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

    def HydrogenBonds(self, rmax=3.2, list_info=False, ndb=None):
        if not ndb:
            ndb = backbone(self).NeighborDatabase()
        bonds = self.PolarContacts(rmax)
        mchb = self.BackboneHbonds()
        poplist = []
        for (i, (i1, i2, d)) in enumerate(bonds):
            flag = False
            resid1, resid2 = self.GetAtomResID(i1), self.GetAtomResID(i2)
            name1, name2 = self.GetAtom(i1).GetName(), self.GetAtom(i2).GetName()
            try:
                if resid1 in ndb[resid2]:
                    flag = True
                elif name1=='N' and name2=='O':
                    if resid1 in mchb:
                        if mchb[resid1][0]!=resid2:
                            flag = True
                elif name2=='N' and name1=='O':
                    if resid2 in mchb:
                        if mchb[resid2][0]!=resid1:
                            flag = True
            except KeyError:
                pass
            if not flag:
                res1, res2 = self.GetAtom(i1).get_res_name(), self.GetAtom(i2).get_res_name()
                flag = True
                try:
                    if pdbnames.MayHBond(res1,name1,res2,name2):
                        flag = False
                except KeyError:
                    flag = False
            if flag:
                poplist.append(i)
        for i in poplist[::-1]:
            del bonds[i]
        if list_info:
            infobonds = []
            for (a1, a2, d) in bonds:
                atom1, atom2 = self.GetAtomTitle(a1), self.GetAtomTitle(a2)
                infobonds.append((atom1[:8].strip(),atom1[8:].strip(),atom2[:8].strip(),atom2[8:].strip()))
            return infobonds
        else:
            return bonds

    def CACAdistances(self):
        ind = self.atom_lister('name_CA')
        resids, d = [], []
        for (i,atomi) in enumerate(ind):
            d.append([])
            resids.append(self.get_atom_alt_resid(atomi).replace(' ',''))
            for (j,atomj) in enumerate(ind):
                if j<i:
                    d[i].append(d[j][i])
                elif j==i:
                    d[i].append(0.0)
                else:
                    d[i].append(self.distance(atomi,atomj))
        return (resids, d)

    def ABdistances(self, nameA, nameB):
        indA = self.atom_lister('name_'+nameA)
        indB = self.atom_lister('name_'+nameB)
        resids, d = [[],[]], []
        for (i,atomi) in enumerate(indA):
            d.append([])
            resids[0].append(self.GetAtomResID(atomi))
            for (j,atomj) in enumerate(indB):
                d[i].append(self.distance(atomi,atomj))
        for atom in indB:
            resids[1].append(self.GetAtomResID(atom))
        return (resids, d)

    def IndexDistances(self, indA, indB):
        d = []
        for atomi in indA:
            d.append([])
            for atomj in indB:
                d[-1].append(self.distance(atomi,atomj))
        return d

    def IsProteinAtom(self,i):
        return pdbnames.Is3Amino(self.GetAtomResidueName(i))

    def IsProteinBackboneAtom(self, i):
        return self.atoms[i].IsProteinBackbone()

    def IsProteinSideChainAtom(self, i):
        return self.atoms[i].IsProteinSidechain()

    def IsBackboneAtom(self, i):
        return self.atoms[i].IsBackbone()

    def IsSideChainAtom(self, i):
        return self.atoms[i].IsSidechain()

    def IsHeteroAtom(self,i):
        return pdbnames.IsHetero(self.GetAtomResidueName(i))

    def IsWater(self, i):
        return pdbnames.IsWater(self.GetAtomResidueName(i))

    def BackboneHbonds(self):
        bb=backbone(self)
        return bb.MChbonds()

    def ChiList(self):
        chis = {}
        for residue in self.GetProteinResidues().values():
            chichi = residue.GetChis()
            for chi in chichi.keys():
                if not chichi[chi]:
                    chichi.pop(chi)
            if chichi:
                chis[residue.GetResID()+'|'+residue.get_res_name()] = chichi
            if residue.HasAltConf():
                for ac in residue.GetAltCodes()[1:]:
                    chichi = residue.GetAltChis(ac)
                    for chi in chichi.keys():
                        if not chichi[chi]:
                            chichi.pop(chi)
                if chichi:
                    chis[residue.GetResID()+ac+'|'+residue.get_res_name()] = chichi
        return chis

    def bperesidue(self):
        for residue in self.get_residues().values():
            residue.Baverage()

    def GetBlimits(self):
        b = []
        for atom in self.atoms:
            b.append(atom.GetB())
        return (min(b),max(b))

    def resid_compare(self, resid1, resid2):
        chid = cmp(resid1[0], resid2[0])
        if chid:
            return chid
        try:
            resn1 = int(resid1[1:])
        except ValueError:
            resn1 = int(resid1[1:-1])
        try:
            resn2 = int(resid2[1:])
        except ValueError:
            resn2 = int(resid2[1:-1])
        resnum = cmp(resn1,resn2)
        if resnum:
            return resnum
        return cmp(resid1,resid2)

    def get_occupancy_estimates(self):
        occs = {}
        for resid in self.ListAltResidues():
            occs[resid] = self.GetResidue(resid).esitimate_ac_occupancies()
        return occs

    def print_occupancy_estimates(self):
        occs = self.get_occupancy_estimates()
        for resid in sorted(occs):
            print(resid + ' %5.2f'*len(occs[resid]) % tuple(occs[resid]))

    def xyz(self, listik=False):
        listik = self.__ensure_listik_(listik)
        return array(map(lambda x : x.GetR(), self.atom_getter('all', listik)))

class backbone:
    '''
    Protein backbone atoms manipulation.
    '''
    def __init__(self, molecule):
        residues = molecule.GetProteinResidues()
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
            dzip = zip(*d)
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
#        reslist = self.ListProteinResidues()
#        residues = self.GetProteinResidues()
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
        return math.degrees(math.acos((r1*r2).sum()/math.sqrt((r1**2).sum()*(r2**2).sum())))
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
        ind1 = range(len(self.altcodes))
        ind2 = range(len(other.altcodes))
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
            keys = chis.keys()
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
        atoms = []
        for name in self.atoms.keys():
            if re.match(mask,name):
                atoms.append(self.atoms[name])
        return atoms

    def GetMaskAltAtoms(self, mask, code):
        atoms, ac = [], self.altcodes.index(code)
        for name in self.altatoms[ac].keys():
            if re.match(mask,name):
                atoms.append(self.altatoms[ac][name])
        if atoms:
            return atoms
        else:
            return self.GetMaskAtoms(mask)

    def GetAtom(self, name):
        try:
            return self.atoms[name]
        except KeyError:
            return None

    def GetAtoms(self):
        return self.atoms.values()

    def get_atom_names(self):
        ''' Return the list of atom names in the residue, sorted in 
            random order. '''
        return self.atoms.keys()

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
        b = array(map(lambda x : x.GetB(), self.origatoms))
        o = array(map(lambda x : x.GetOccupancy(), self.origatoms))
        f = array(map(lambda x : x.IsBackbone(), self.origatoms)).astype(int)
        if mode == 'lin':
            return [sum(b*o)/sum(o),sum(b*f*o)/sum(f*o),sum(b*(1-f)*o)/sum((1-f)*o)]
        elif mode == 'rms':
            return [sqrt(sum(o*b**2)/sum(o)),sqrt(sum(f*o*b**2)/sum(f*o)),sqrt(sum((1-f)*o*b**2)/sum((1-f)*o))]
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
            names1 = self.atoms.keys()
        if not names2:
            names2 = self.atoms.keys()
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
                oxys = cb[key].keys()
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
    def BondsAnglesTorsions(self, Rcutoff=2.0, printout=False):
        Nats = len(self.origatoms)
        ps = array(reduce(lambda x,y: x+y, [[(k,i) for i in range(k+1,Nats)] for k in range(Nats)]))
        ds = array(map(lambda x : self.distance(x[0],x[1]), ps))
        bonds = ps[ds<Rcutoff]
        Nbs = len(bonds)
        ps = array(reduce(lambda x,y: x+y, [[(k,i) for i in range(k+1,Nbs)] for k in range(Nbs)]))
        angles = array(map(lambda x : array(x)[array([x.count(y) for y in x]).argsort()[array([0,2,1])]].tolist(), filter(lambda x : len(set(x))==3, map(lambda x : bonds[x[0]].tolist()+bonds[x[1]].tolist(), ps))))
        Nas = len(angles)
        ps = array(reduce(lambda x,y: x+y, [[(k,i) for i in range(k+1,Nas)] for k in range(Nas)]))
        torsions = array(map(lambda x : x[:4] if x.tolist().count(x[-1])==2 else x[:3].tolist()+[x[-1]], array(map(lambda x : x[:3][::-1].tolist()+x[3:].tolist() if x.tolist().count(x[0])==2 else x,filter(lambda x : x[1]!=x[4], array(filter(lambda x : len(set(x))==4, map(lambda x : angles[x[0]].tolist()+angles[x[1]].tolist(), ps))))))))
        impropers = array(list(set(map(lambda x : (x[1],x[0],x[2],x[3]), map(lambda x : [x[1]]+sorted(set(x)-set([x[1]])), filter(lambda x : x[1]==x[4], filter(lambda x : len(set(x))==4, map(lambda x : angles[x[0]].tolist()+angles[x[1]].tolist(), ps))))))))
        if printout:
            bonds = map(lambda x : '%5s --- %5s : %10.3f' % (self.origatoms[x[0]].name(), self.origatoms[x[1]].name(), self.distance(x[0],x[1])), bonds)
            angles = map(lambda x : '%5s - %5s - %5s : %10.2f' % (self.origatoms[x[0]].name(), self.origatoms[x[1]].name(), self.origatoms[x[2]].name(), self.angle(x[0],x[1],x[2])), angles)
            torsions = map(lambda x : '%5s - %5s - %5s - %5s : %10.2f' % (self.origatoms[x[0]].name(), self.origatoms[x[1]].name(), self.origatoms[x[2]].name(), self.origatoms[x[3]].name(), self.torsion(x[0],x[1],x[2],x[3])), torsions)
            impropers = map(lambda x : '%5s - %5s - %5s - %5s : %10.2f' % (self.origatoms[x[0]].name(), self.origatoms[x[1]].name(), self.origatoms[x[2]].name(), self.origatoms[x[3]].name(), self.torsion(x[0],x[1],x[2],x[3])), impropers)
        return bonds, angles, torsions, impropers
#class idatom:

#    def __init__(self, chainid = ' ', resid = 1, atomn = 'CA', altloc = ' ', iCode = ' '):
#        self.chid = chainid[0]
#        self.resi = int(resid)
#        self.name = atomn.strip()
#        self.altc = altloc
#        self.icod = iCode

#    def ResID(self):
#        return self.chid + '%4d' % self.resi + self.icod

#    def GetAtomName(self):
#        return self.name

#    def GetAltLoc(self):
#        return self.altc
        
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
        return self.remarks.keys()

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
    return math.degrees(math.acos((r1*r2).sum()/math.sqrt((r1**2).sum()*(r2**2).sum())))

def angle2(i,j,k):
    r1 = i-j
    r2 = k-j
    return math.degrees(math.acos((r1*r2).sum()/math.sqrt((r1**2).sum()*(r2**2).sum())))

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
