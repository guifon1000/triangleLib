from sympy.abc import x
import sympy as sp
from sympy.utilities.lambdify import lambdify, implemented_function
from sympy.utilities.lambdify import lambdastr
import numpy as np
import matplotlib.pyplot as plt
from sympy import Function
from sympy.solvers.solvers import solve_linear_system




class function_0(object):
    def __init__(self,v):
        self.f = lambdify( x , v)
        self.v = v
    def plot(self,t,step=0.1):
        ordo = []
        absc=np.arange(t[0],t[1],step)
        for xi in absc:
            ordo.append(self.f(xi))
        plt.plot(absc,ordo)
    def __str__(self):
        return str(self.v)
    def __add__(self,other):
        return function_0(self.v + other.v)
    def __mul__(self,other):
        return function_0(self.v * other.v)
    def __rmul__(self,other):
        if hasattr(other,'v'):return function_0(self.v * other.v)
        else:return function_0(self.v * other)
    def __call__(self,v):
        return self.v(v)


class basisFunction(function_0):
    def __init__(self,i):
        if np.mod(i,2)==0:super(basisFunction,self).__init__(sp.cos(float(i)*x))
        else:super(basisFunction,self).__init__(sp.sin(float(i)*x))

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




if __name__=='__maiiin__':
    Omega=(0,10)
    targ = function_0(2*x+sp.cos(x)+0.001*x**2)
    f0 = function_0(x+sp.sin(x))
    plt.clf()
    targ.plot(Omega)
    f0.plot(Omega)
    v = sp.integrate(f0*f0,(x,Omega[0],Omega[1]))
    u = sp.integrate(f0*targ, (x, Omega[0], Omega[1]))
    print float(u/v)
    test = float(u/v)*f0
    test.plot(Omega)


    plt.show()

if __name__=='__main__':
    Omega=(-2,2)
    N=15
    f0=function_0(0)
    targ = function_0(x*x*x+5*x**2+9*x+18)
    #alpha=[]
    basis = Basis(N)
    
    
    mat = sp.zeros(N,N)
    rht = sp.zeros(N,1)
    for i in range(N):
        fi = basis.ba[i]
        for j in range(i,N):
            fj = basis.ba[j]
            mat[i,j]= sp.integrate(fi*fj,(x, Omega[0], Omega[1]))
            print i,j
            mat[j,i]=mat[i,j]
        ff = fi*targ
        rht[i,0] = sp.integrate(ff, (x, Omega[0], Omega[1]))
    print np.array(mat)
    print '-----------------'
    print rht
    print '-----------------'
    alpha=mat.LUsolve(rht)    
    #print sp.solvers.solve_linear_system_LU(mat,rht)
    print '####################################'
    plt.clf()
    for i,f in enumerate(basis.ba):
        f0+=function_0(alpha[i])*f
    g=function_0(f0)
    g.plot(Omega)
    targ.plot(Omega)
    plt.show()

    
