from sympy.abc import x
import sympy as sp
from sympy.utilities.lambdify import lambdify, implemented_function
from sympy.utilities.lambdify import lambdastr
import numpy as np
import matplotlib.pyplot as plt
from sympy import Function




class function_0(object):
    def __init__(self,v):
        self.f = lambdify( x , v)
        self.v = v
    def plot(self,t,step=0.1):
        plt.clf()
        ordo = []
        absc=np.arange(t[0],t[1],step)
        for xi in absc:
            ordo.append(self.f(xi))
        plt.plot(absc,ordo)
        plt.show()
    def __str__(self):
        return str(self.v)
    def __add__(self,other):
        return function_0(self.v + other.v)
    def __mul__(self,other):
        return function_0(self.v * other.v)
    def __rmul__(self,other):
        return function_0(self.v * other.v)
    def __call__(self,v):
        return self.v(v)


class basisFunction(function_0):
    def __init__(self,i):
        super(basisFunction,self).__init__(sp.cos(float(i)*x))

class Basis:
    def __init__(self,N):
        ba = []
        for i in range(N):
            ba.append(basisFunction(i))
        self.ba = ba



def family(N):
    fam = [None]*N
    for i in range(N):
        fam[i] = basisFunction(i)
    return fam


if __name__=='__main__':
    Omega=(0,1)
    N=10
    f0=function_0(0)
    alpha=[]
    basis = Basis(N)
    for i in range(N):
        alpha.append(1./float(i+1))
    
    
    mat = np.zeros((N,N))
    for i in range(N):
        fi = basis.ba[i]
        for j in range(N):
            fj = basis.ba[j]
            mat[i,j]=sp.integrate(fi(x, Omega[0], Omega[1])*fj((x, Omega[0], Omega[1])))


    for i,f in enumerate(basis.ba):
      f0+=function_0(alpha[i])*f
    g=function_0(f0)
    g.plot((0,100),0.001)

    
