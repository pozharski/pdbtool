from numpy import array, matrix, sqrt, arccos, degrees, pi, diag, matmul
from numpy.linalg import eig, inv

class TInertia:
    '''
    Inertia tensor calculations for a system of point masses.
    '''
    def __init__(self, m=None, x=None, y=None, z=None, xyz=None, mxyz=None):
        '''
        Class requires an array of masses and three arrayswith coordinates.
        Alternatively, Nx3 or 3xN array xyz can be supplied to define
        coordinates.  If mxyz is provided, it should be 4xN or Nx4 array.
        In special cases of 3x3 of 4x4 arrays, parameters are expected to be 
        arranged in rows.
        '''
        if mxyz is not None:
            if mxyz.shape.index(4) == 0:
                m,x,y,z = mxyz
            else:
                m,x,y,z = mxyz.T
        if xyz is not None:
            if xyz.shape.index(3) == 0:
                x,y,z = xyz
            else:
                x,y,z = mxyz.T
        self.mass = sum(m)
        if self.mass > 0.0:
            self.comass = array([   sum(m*x)/self.mass,
                                    sum(m*y)/self.mass,
                                    sum(m*z)/self.mass ])
            self.m = m
            self.x = x - self.comass[0]
            self.y = y - self.comass[1]
            self.z = z - self.comass[2]
            Ixx = sum(m*(self.y**2+self.z**2))
            Iyy = sum(m*(self.x**2+self.z**2))
            Izz = sum(m*(self.x**2+self.y**2))
            Ixy = - sum(m*self.x*self.y)
            Iyz = - sum(m*self.y*self.z)
            Ixz = - sum(m*self.x*self.z)
            self.Ixyz = matrix(array([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]]))
            self.w, self.vr = eig(self.Ixyz)
            self.abc = sqrt(((2.5*self.w/self.mass)*array([[-1,1,1],[1,-1,1],[1,1,-1]])).sum(1)).real
            self.vr = self.vr.T[self.abc.argsort()[::-1]]
            self.w = self.w[self.abc.argsort()[::-1]]
            self.abc = self.abc[self.abc.argsort()[::-1]]
            self.eccent = [self.abc[0], self.abc[0]/sqrt(self.abc[1]**2+self.abc[2]**2), self.abc[1]/self.abc[2]]
            self.valid = True
        else:
            self.valid = False

    def angles(self, other):
        return degrees(arccos(abs((self.vr*other.vr).sum(1))))

    def shift3(self, other):
        return other.comass-self.comass

    def shift(self, other):
        return sqrt(sum(self.shift3(other)**2))

    def report(self):
        retline  = 'Total mass = %.1f\n' % self.mass
        retline += 'Center of mass coordinates:\n'
        retline += ' Xcom = %.1f\n Ycom = %.1f\n Zcom = %.1f\n' % tuple(self.comass)
        retline += 'Tensor of inertia:\n'
        Ixyz = self.Ixyz.getA()
        retline += '  %15.2f %15.2f %15.2f\n' % (Ixyz[0][0],Ixyz[0][1],Ixyz[0][2])
        retline += '  %15.2f %15.2f %15.2f\n' % (Ixyz[1][0],Ixyz[1][1],Ixyz[1][2])
        retline += '  %15.2f %15.2f %15.2f\n' % (Ixyz[2][0],Ixyz[2][1],Ixyz[2][2])
        retline += 'Inertia ellipsoid:\n'
        retline += ' Axis 1 (%10.3f) : %10.3f %10.3f %10.3f\n' % tuple([self.abc[0]]+self.vr[0].tolist())
        retline += ' Axis 2 (%10.3f) : %10.3f %10.3f %10.3f\n' % tuple([self.abc[1]]+self.vr[1].tolist())
        retline += ' Axis 3 (%10.3f) : %10.3f %10.3f %10.3f\n' % tuple([self.abc[2]]+self.vr[2].tolist())
        return retline

    def numcomp(self, other, prefix=''):
        retline  = prefix+'_MASS %.2f %.2f\n' % (self.mass, other.mass)
        retline += prefix+'_COM %.2f %.2f %.2f\n' % tuple(self.shift3(other))
        retline += prefix+'_ABC %.2f %.2f %.2f %.2f %.2f %.2f\n' % tuple(self.abc.tolist()+other.abc.tolist())
        return retline

    def reportcomp(self, other):
        retline  = 'Total mass change: %.1f vs %.1f\n' % (self.mass, other.mass)
        if self.mass != other.mass:
            retline += 'WARNING: segment mass have changed. Results will reflect both changes in segment orientation and composition.\n'
        retline += 'Center of mass translation:\n %.1f %.1f %.1f\n' % tuple(self.shift3(other))
        retline += 'Inertia ellipsoid dimensions: (%.1f %.1f %.1f) -> (%.1f %.1f %.1f)\n' % tuple(self.abc.tolist()+other.abc.tolist())
        return retline

    def frame(self):
        retline  = 'ATOM      1  X   UNK M   1    %8.3f%8.3f%8.3f  1.00 20.00              \n' % tuple(self.comass)
        retline += 'ATOM      2  X   UNK A  -1    %8.3f%8.3f%8.3f  1.00 20.00              \n' % tuple(self.comass-self.abc[0]*self.vr[0])
        retline += 'ATOM      3  X   UNK A   1    %8.3f%8.3f%8.3f  1.00 20.00              \n' % tuple(self.comass+self.abc[0]*self.vr[0])
        retline += 'ATOM      4  X   UNK B  -1    %8.3f%8.3f%8.3f  1.00 20.00              \n' % tuple(self.comass-self.abc[1]*self.vr[1])
        retline += 'ATOM      5  X   UNK B   1    %8.3f%8.3f%8.3f  1.00 20.00              \n' % tuple(self.comass+self.abc[1]*self.vr[1])
        retline += 'ATOM      6  X   UNK C  -1    %8.3f%8.3f%8.3f  1.00 20.00              \n' % tuple(self.comass-self.abc[2]*self.vr[2])
        retline += 'ATOM      7  X   UNK C   1    %8.3f%8.3f%8.3f  1.00 20.00              \n' % tuple(self.comass+self.abc[2]*self.vr[2])
        return retline

    def anisou(self, serid=1, aname='X', altloc=' ', rname='UNK', chid='X',rid=1,icode=' ',element='X',charge=''):
        Exyz = matmul(self.vr.T,matmul(diag(self.abc**2),self.vr))
        retline  = 'ATOM  %5d %4s%1s%3s ' % (serid, aname,altloc,rname)
        retline += '%1s%4d%1s' % (chid,rid,icode)
        retline += '   %8.3f%8.3f%8.3f  1.00' % tuple(self.comass)
        retline += '%6.2f' % (8*pi**2*(sum(self.abc)))
        retline += '          %2s%2s\n' % (element,charge)
        retline += 'ANISOU%5d %4s%1s%3s ' % (serid, aname,altloc,rname)
        retline += '%1s%4d%1s' % (chid,rid,icode)
        retline += ' %7d%7d%7d' % (Exyz[0][0],Exyz[1][1],Exyz[2][2])
        retline += '%7d%7d%7d  ' % (Exyz[0][1],Exyz[0][2],Exyz[1][2])
        retline += '    %2s%2s\n' % (element,charge)
        return retline
