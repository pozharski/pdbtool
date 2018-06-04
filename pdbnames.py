AMINO_ACIDS = {
'ALA':['C','CA','CB','N','O'],
'CYS':['C','CA','CB','N','O','SG'],
'ASP':['C','CA','CB','CG','N','O','OD1','OD2'],
'GLU':['C','CA','CB','CD','CG','N','O','OE1','OE2'],
'PHE':['C','CA','CB','CD1','CD2','CE1','CE2','CG','CZ','N','O'],
'GLY':['C','CA','N','O'],
'HIS':['C','CA','CB','CD2','CE1','CG','N','ND1','NE2','O'],
'ILE':['C','CA','CB','CD1','CG1','CG2','N','O'],
'LYS':['C','CA','CB','CD','CE','CG','N','NZ','O'],
'LEU':['C','CA','CB','CD1','CD2','CG','N','O'],
'MET':['C','CA','CB','CE','CG','N','O','SD'],
'ASN':['C','CA','CB','CG','N','ND2','O','OD1'],
'PRO':['C','CA','CB','CD','CG','N','O'],
'GLN':['C','CA','CB','CD','CG','N','NE2','O','OE1'],
'ARG':['C','CA','CB','CD','CG','CZ','N','NE','NH1','NH2','O'],
'SER':['C','CA','CB','N','O','OG'],
'THR':['C','CA','CB','CG2','N','O','OG1'],
'VAL':['C','CA','CB','CG1','CG2','N','O'],
'TRP':['C','CA','CB','CD1','CD2','CE2','CE3','CG','CH2','CZ2','CZ3','N','NE1','O'],
'TYR':['C','CA','CB','CD1','CD2','CE1','CE2','CG','CZ','N','O','OH'],
'CME':['C','CA','CB','CE','CZ','N','O','OH','SD','SG'],
'TPO':['C','CA','CB','CG2','N','O','OG1','P','O1P','O2P','O3P'],
'PTR':['C','CA','CB','CD1','CD2','CE1','CE2','CG','CZ','N','O','OH','P','O1P','O2P','O3P'],
'UNK':['C','CA','CB','N','O'],
'MSE':['C','CA','CB','CE','CG','N','O','SE'],
'DMK':['C','CA','CB','CG1','CG2','CG3','N','O','OD1','OD2'],
'DMH':['C','CA','CB','CG','CE1','CE2','N','ND2','O','OD1'],
}

DONORS = {
'ALA':['N'],
'CYS':['N','SG'],
'ASP':['N'],
'GLU':['N'],
'PHE':['N'],
'GLY':['N'],
'HIS':['N','ND1','NE2'],
'ILE':['N'],
'LYS':['N','NZ'],
'LEU':['N'],
'MET':['N'],
'ASN':['N','ND2'],
'PRO':['N'],
'GLN':['N','NE2'],
'ARG':['N','NE','NH1','NH2'],
'SER':['N','OG'],
'THR':['N','OG1'],
'VAL':['N'],
'TRP':['N','NE1'],
'TYR':['N','OH'],
'CME':['N','OH'],
'TPO':['N','OG1','O1P','O2P','O3P'],
'PTR':['N','OH','O1P','O2P','O3P'],
'HOH':['O']
}

ACCEPTORS = {
'ALA':['O'],
'CYS':['O','SG'],
'ASP':['O','OD1','OD2'],
'GLU':['O','OE1','OE2'],
'PHE':['O'],
'GLY':['O'],
'HIS':['ND1','O'],
'ILE':['O'],
'LYS':['O'],
'LEU':['O'],
'MET':['O','SD'],
'ASN':['O','OD1'],
'PRO':['O'],
'GLN':['O','OE1'],
'ARG':['O'],
'SER':['O','OG'],
'THR':['O','OG1'],
'VAL':['O'],
'TRP':['O'],
'TYR':['O','OH'],
'CME':['O','OH', 'SD', 'SG'],
'HOH':['O'],
'SO4':['O1','O2','O3','O4'],
'TPO':['O','OG1','O1P','O2P','O3P'],
'PTR':['O','OH','O1P','O2P','O3P']
}

EXTRA_DONORS_ACCEPTORS = {
'ASP':['OD1','OD2'],
'GLU':['OE1','OE2'],
'HIS':['CD2','CE1'],
'ASN':['ND2','OD1'],
'GLN':['NE2','OE1'],
}

BBLIST = ['CA','C','O','N','P','O1P','O2P','O3*','O4*','O5*','C1*','C2*','C3*','C4*','C5*',"O3'","O4'","O5'","C1'","C2'","C3'","C4'","C5*"]

CHIS = {
'CYS': ['chi1'], 
'ASP': ['chi1', 'chi21', 'chi22'], 
'SER': ['chi1'], 
'GLN': ['chi1', 'chi2', 'chi31', 'chi32'], 
'LYS': ['chi1', 'chi2', 'chi3', 'chi4'], 
'PRO': ['chi1', 'chi2'], 
'THR': ['chi11', 'chi12'], 
'PHE': ['chi1', 'chi21', 'chi22'], 
'ALA': [], 
'HIS': ['chi1', 'chi21', 'chi22'], 
'GLY': [], 
'ILE': ['chi11', 'chi12', 'chi2'], 
'GLU': ['chi1', 'chi2', 'chi31', 'chi32'], 
'LEU': ['chi1', 'chi21', 'chi22'], 
'ARG': ['chi1', 'chi2', 'chi3', 'chi4', 'chi51', 'chi52'], 
'TRP': ['chi1', 'chi21', 'chi22'], 
'VAL': ['chi11', 'chi12'], 
'ASN': ['chi1', 'chi21', 'chi22'], 
'TYR': ['chi1', 'chi21', 'chi22'], 
'MET': ['chi1', 'chi2', 'chi3'],
'TPO': ['chi11', 'chi12'],
'PTR': ['chi1', 'chi21', 'chi22']
}

CHIMIN = {
'CYS': ['chi1'], 
'ASP': ['chi1', 'chi21'], 
'SER': ['chi1'], 
'GLN': ['chi1', 'chi2', 'chi31'], 
'LYS': ['chi1', 'chi2', 'chi3', 'chi4'], 
'PRO': ['chi1', 'chi2'], 
'THR': ['chi11'], 
'PHE': ['chi1', 'chi21'], 
'ALA': [], 
'HIS': ['chi1', 'chi21'], 
'GLY': [], 
'ILE': ['chi11', 'chi2'], 
'GLU': ['chi1', 'chi2', 'chi31'], 
'LEU': ['chi1', 'chi21'], 
'ARG': ['chi1', 'chi2', 'chi3', 'chi4'], 
'TRP': ['chi1', 'chi21'], 
'VAL': ['chi11'], 
'ASN': ['chi1', 'chi21'], 
'TYR': ['chi1', 'chi21'], 
'MET': ['chi1', 'chi2', 'chi3'],
'TPO': ['chi11'],
'PTR': ['chi1', 'chi21', 'chi62', 'chi73']
}

AA_SIDE = {
'ALA':['CB'],
'CYS':['CB','SG'],
'ASP':['CB','CG','OD1','OD2'],
'GLU':['CB','CD','CG','OE1','OE2'],
'PHE':['CB','CD1','CD2','CE1','CE2','CG','CZ'],
'GLY':[],
'HIS':['CB','CD2','CE1','CG','N','ND1','NE2','O'],
'ILE':['CB','CD1','CG1','CG2'],
'LYS':['CB','CD','CE','CG','N','NZ','O'],
'LEU':['CB','CD1','CD2','CG'],
'MET':['CB','CE','CG','SD'],
'ASN':['CB','CG','N','ND2','O','OD1'],
'PRO':['CB','CD','CG'],
'GLN':['CB','CD','CG','N','NE2','O','OE1'],
'ARG':['CB','CD','CG','CZ','N','NE','NH1','NH2','O'],
'SER':['CB','OG'],
'THR':['CB','CG2','OG1'],
'VAL':['CB','CG1','CG2'],
'TRP':['CB','CD1','CD2','CE2','CE3','CG','CH2','CZ2','CZ3','N','NE1','O'],
'TYR':['CB','CD1','CD2','CE1','CE2','CG','CZ','OH'],
'TPO':['CB','CG2','OG1','P','O1P','O2P','O3P'],
'PTR':['CB','CD1','CD2','CE1','CE2','CG','CZ','OH','P','O1P','O2P','O3P']

}

MASS = {
'H' : 1.008,
'HE': 4.003,
'LI': 6.941,
'C' : 12.011,
'N' : 14.007,
'O' : 15.999,
'F' : 18.998,
'NA': 22.990,
'MG': 24.305,
'P' : 30.974,
'S' : 32.066,
'CL': 34.453,
'K' : 39.098,
'CA': 40.078,
'Ca': 40.078,
'MN': 54.938,
'FE': 55.845,
'CO': 58.933,
'NI': 58.693,
'CU': 63.546,
'ZN': 65.380,
'SE': 78.960,
'HG': 200.59,
'PB': 207.2,
'BR': 79.904,
'I' : 126.9045
}

VDWRADIUS = {
'H' : 1.20,
'HE': 1.40,
'LI': 1.82,
'C' : 1.70,
'N' : 1.55,
'O' : 1.52,
'F' : 1.47,
'NE': 1.54,
'NA': 2.27,
'MG': 1.73,
'SI': 2.10,
'P' : 1.80,
'S' : 1.80,
'CL': 1.75,
'K' : 2.75,
'CA': 2.40,
'Ca': 2.40,
'MN': 2.00,
'FE': 2.00,
'CO': 2.00,
'NI': 1.63,
'CU': 1.40,
'ZN': 1.39,
'SE': 1.90,
'HG': 1.55,
'PB': 2.02
}

POLARS = ['F', 'SE', 'CL',  'N', 'O', 'P', 'S']
METALS = ['MN', 'MG', 'HG', 'CA', 'Ca', 'CL', 'K', 'FE', 'CO', 'NA', 'LI', 'ZN', 'PB', 'CU', 'NI']

def Is3Amino(name):
    return name in AMINO_ACIDS
def IsWater(name):
    return name == 'HOH'
def NotWater(name):
    return name != 'HOH'
def IsHetero(name):
    return not (Is3Amino(name) or IsWater(name))
def MaybeBackbone(name):
    return name in BBLIST
def NotBackbone(name):
    return name not in BBLIST
def IsDonor(resn, atom):
    return atom in DONORS[resn]
def IsAcceptor(resn, atom):
    return atom in ACCEPTORS[resn]
def MayHBond(resn1,atom1,resn2,atom2):
    return (IsDonor(resn1,atom1) and IsAcceptor(resn2,atom2)) or (IsDonor(resn2,atom2) and IsAcceptor(resn1,atom1))
def GetMass(name):
    return MASS[name]
def GetVDWRadius(name):
    return VDWRADIUS[name]
def IsPolar(name):
    return name in POLARS
def IsMetal(name):
    return name in METALS
