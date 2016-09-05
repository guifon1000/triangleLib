from math import hypot
import numpy as np
#import sympy

class Quaternion:
    def __init__(self,a,b):
        self.a=a
        self.b=b

    def __str__(self):
        aff='('
        aff+=str(self.a.real)+')+('
        aff+=str(self.a.imag)+')i+('
        aff+=str(self.b.real)+')j+('
        aff+=str(self.b.imag)+')k'
        return aff

    def __neg__(self):
        return Quaternion(-self.a,-self.b)

    def __add__(self,other):
        return Quaternion(self.a+other.a,self.b+other.b)

    def __sub__(self,other):
        return Quaternion(self.a-other.a,self.b-other.b)

    def __mul__(self,other):
        c=self.a*other.a-self.b*other.b.conjugate()
        d=self.a*other.b+self.b*other.a.conjugate()
        return Quaternion(c,d)

    def __rmul__(self,k):
        return Quaternion(self.a*k,self.b*k)

    def __abs__(self):
        return hypot(abs(self.a),abs(self.b))

    def conjugate(self):
        return Quaternion(self.a.conjugate(),-self.b)

    def __div__(self,other):
        return self*(1./abs(other)**2*other.conjugate())

    def __pow__(self,n):
        r=1
        for i in range(n):
            r=r*self
        return r
    def quaternion2matrix(self):
        mat=np.zeros((3,3),dtype=np.float)
        x=self.a.real
        y=self.a.imag
        z=self.b.real
        w=self.b.imag
        mat[0,0]=1.-2*y**2.-2.*z**2. 
        mat[0,1]=2.*x*y-2.*z*w
        mat[0,2]=2.*x*z+2.*y*w
        mat[1,0]=2.*x*y+2.*z*w
        mat[1,1]=1.-2*x**2.-2.*z**2.
        mat[1,2]=2.*y*z-2.*x*w
        mat[2,0]=2.*x*z-2.*y*w
        mat[2,1]=2.*y*z+2.*x*w
        mat[2,2]=1.-2*x**2.-2.*y**2.
        return mat


def matrix2quaternion(m):
    diag=np.diag(m)
    np.append(diag,1.)
    
    tr= np.trace(m)+1.
    if tr>0.:
        s=0.5/np.sqrt(tr)
        x=(m[2,1]-m[1,2])*s
        y=(m[0,2]-m[2,0])*s
        z=(m[1,0]-m[0,1])*s
        w=0.25/s
        return Quaternion(np.complex(w,x),np.complex(y,z))
         

if __name__=='__main__':
    a=np.complex(0.5,-0.5)
    b=np.complex(-0.5,-0.5)
    q=Quaternion(a,b)
    print q
    print '-------------------------------'
    m=q.quaternion2matrix()
    print m
    print '-------------------------------'
    print matrix2quaternion(m).quaternion2matrix()-m
