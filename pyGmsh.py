import numpy as np

def read4CornersSurfMSH(name):
    f=open(name,'r').readlines()
    points=[]
    faces=[]
    for i,l in enumerate(f):
        if '$MeshFormat' in l:
            vers=f[i+1]
        if '$Nodes' in l:
            nb = int(f[i+1])
            for j in range(i+2,i+2+nb):
                points.append(f[j])
        if '$Elements' in l:
            nb = int(f[i+1])
            for j in range(i+2,i+2+nb):
                if f[j].split()[1]=='3':faces.append(f[j])
    print faces


if __name__=='__main__':
    read4CornersSurfMSH('testGMSH.msh')

