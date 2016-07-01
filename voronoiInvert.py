import pyevtk.hl as pv
import numpy as np
import triangles_to_VTK as tv
import triangleLib as tl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time


tri=[]

k=[]
iStart=-1
iEnd=-1

for i,l in enumerate(open("VOILES.TEC","r")):
    if "ZONE" in l:
        print i
        if iStart==-1:
            iStart=i+1
        if iStart!=-1 and i > iStart:
            iEnd=i

print iStart, iEnd
first=[]
faces=[]
f_voiles=open("VOILES.TEC","r").readlines()[iStart:iEnd]

for l in f_voiles:
    if len(l.split())>3:
        vs=l.split()
        v=[]
        for s in vs:
            v.append(float(s))
        first.append(v)
    if len(l.split())==3:
        faces.append([int(l.split()[0])-1,int(l.split()[1])-1,int(l.split()[2])-1])

maxipt=-1
for v in faces:maxipt=max(np.max(v),maxipt)
print maxipt
print len(first)

x=[]
y=[]
z=[]

for e in first:
    x.append(e[0])
    y.append(e[1])
    z.append(e[2])
xs=np.array(x)
ys=np.array(y)
zs=np.array(z)

faces=np.array(faces)



def distance(p0,p1):
    d2=(p1[0]-p0[0])**2.+(p1[1]-p0[1])**2.+(p1[2]-p0[2])**2.
    return np.sqrt(d2)



PrincipalStressX=[]
PrincipalStressY=[]



#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

k=[]


for i,l in enumerate(open('CHAMPS_PLIS.TEC',"r").readlines()):
    if "ZONE" in l and iStart==-1:iStart=i
    if "ZONE" in l and iStart!=-1:iEnd=i


	
for l in open('CHAMPS_PLIS.TEC',"r"):
    if len(l.split())>0:
        if "VARIABLES" in l.split()[0]:
           for v in l.split()[1:]:
               if "=" not in v:k.append(v.replace("\"","").replace(',',''))
d={}
N=len(k)
xm=[]
ym=[]
zm=[]
for n in k:d[n]=[]
for l in open('CHAMPS_PLIS.TEC',"r").readlines()[iStart:iEnd]:
    if ("ZONE" not in l) and ("VARIABLES" not in l) and (len(l.split())==N) :
        val=l.split()
        for i,v in enumerate(k):
            d[k[i]].append(float(val[i]))

        xm.append(float(val[0]))
        ym.append(float(val[1]))
        zm.append(float(val[2]))



lpli=open('CHAMPS_PLIS.TEC',"r").readlines()

princStress=[]

for i,l in enumerate(lpli):
    if "ZONE" in l and 'Contraintes_Principales' in l:
        print l
        for j in range(i+1,len(lpli)):
            l2=lpli[j]
            if "ZONE" in l2 :
                break
            else:
                l3=[float(ll) for ll in l2.split()]
                princStress.append(l3)
while [] in princStress:princStress.remove([])

        





princStress=np.array(princStress)

sfValX=[]
sfValY=[]


N0=len(lpli)
xg=[]
yg=[]
zg=[]
chrono=[]

for f in faces:
    
    #print '------------ triangle -------------'
    t=tl.Triangle()
    p1=tl.Point()
    p2=tl.Point()
    p3=tl.Point()
    i1=f[0]
    p1.setPos(first[i1][0],first[i1][1],first[i1][2])
    i1=f[1]
    p2.setPos(first[i1][0],first[i1][1],first[i1][2])    
    i1=f[2]
    p3.setPos(first[i1][0],first[i1][1],first[i1][2])
    t.setPoints(p1,p2,p3)
    #vor = t.circumCenter()
    vor = t.gravityCenter()
    xg.append(vor.x)
    yg.append(vor.y)    
    zg.append(vor.z)
    Nf = np.size(lpli)
    #lpli=[]
    toRem=[]
    dist=np.ones(Nf)*1000000000.
    ind=[]
    found=False

 
    t0 = time.clock()
    for lm in range(len(princStress)):
        l=princStress[lm]
        pref=tl.Point()
        pref.setPos(float(l[0]),float(l[1]),float(l[2]))
        dist[lm]=tl.distance(pref,vor)
        ind.append(lm)

        
    mini=dist[np.argmin(dist)]
    imin=np.argmin(dist)

    for kk in range(np.size(dist)):
        d=dist[kk]
        
        if d==mini:
            toRem.append(princStress[ind[kk]])
    for e in toRem:
        print e
        np.delete(princStress,e)
    print len(princStress)
    print toRem
    t1 = time.clock()
    chrono.append(t1-t0)

plt.plot(range(len(chrono)),chrono)
plt.show()

pv.pointsToVTK("./pointsGF", np.array(xg), np.array(yg), np.array(zg),data = None)

pv.pointsToVTK("./pointsMD", np.array(xm), np.array(ym), np.array(zm),data = None)



print '---------------------------------------------'


tv.triangle_faces_to_VTK("mesh",
                      x=xs, y=ys, z=zs,
                      faces=faces,
                      point_data=None,
                      cell_data=None)

    
