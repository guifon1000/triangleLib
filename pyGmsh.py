import numpy as np
import triangleLib as tl
import operator



def returnIndex(tab,d):
    j=None
    for i,el in enumerate(tab):
        if el.index == d:
            j=i
    return j

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
                el=f[j].split()
                ind = int(el[0])
                x = float(el[1])
                y = float(el[2])
                z = float(el[3])
                p=tl.Point(x,y,z,ind)
                points.append(p)
    points = sorted(points,key = operator.attrgetter('index'))



    for i,l in enumerate(f):
        if '$Elements' in l:
            nb = int(f[i+1])
            for j in range(i+2,i+2+nb):
                lis = f[j].split()
                last = len(lis)-1
                if lis[1]=='3':
                    pts = []
                    for i in range(4):
                        ploc = int(lis[last-i])
                        posi = returnIndex(points,ploc)
                        p = points[posi]
                        pts.append(p)
                    pan = tl.Panel(pts[0],pts[1],pts[2],pts[3]) 
                    faces.append(pan)
    return points,faces

if __name__=='__main__':
    read4CornersSurfMSH('testGMSH.msh')

