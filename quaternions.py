from math import hypot

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
