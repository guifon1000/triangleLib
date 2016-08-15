from sympy.abc import x
import sympy as sp
from sympy.utilities.lambdify import lambdify, implemented_function
from sympy.utilities.lambdify import lambdastr
import numpy as np
import matplotlib.pyplot as plt



def func(x):
    return np.cos(x)


def family(N):
    fam = [None]*N
    for i in range(N):
        fam[i] = lambdify(x, sp.cos(float(i)*x)   )
    return fam


if __name__=='__main__':
    fam = family(5)
    absc= np.arange(0.,2.*np.pi,0.01)
    plt.clf() 
    for f in fam:
        ordo = []
        for xx in absc: 
            ordo.append(f(xx))
        plt.plot(absc,ordo)
    plt.show()
    
