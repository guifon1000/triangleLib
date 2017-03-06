import bidLib
import triangleLib as tl
import matplotlib.pyplot as plt
import json


class corner(tl.Point):
    def __init__(self,x,y,name):
        super(corner,self).__init__(x,y,0.,index = None,name=name)

class Wall(list):
    def __init__(self,pStart,pEnd,name=None):
        self.p0 = pStart
        self.p1 = pEnd
        self.length = tl.distance(pStart,pEnd)
        if name == None:
            self.name = self.p0.name+self.p1.name
        else : 
            self.name = name
    def plot(self,plt):
        plt.plot([self.p0[0],self.p1[0]],[self.p0[1],self.p1[1]],'-k')
        plt.text(0.5*(self.p0[0]+self.p1[0]),\
                0.5*(self.p0[1]+self.p1[1]),\
                "%.2f" % self.length,color='red')


def points(walls):
    pts = []
    for w in walls:
        if w.p0 not in pts : pts.append(w.p0)
        if w.p1 not in pts : pts.append(w.p1)
    return pts


def scatter(points,plt):
    for p in points :
        plt.scatter(p[0],p[1])
        if p.name is not None:
            plt.text(p[0]+0.2,p[1]+0.2,p.name)#+'('+str(p[0])+';'+str(p[1])+')')


#
#          
#



if __name__=='__main__' :
    xA = 0.
    yA = 4.37+3.52
    A = corner(xA,yA,'A')
    xB = A[0]+19.23
    yB = yA
    B = corner(xB,yB,'B')
    xC = xB
    yC = yB - 7.14
    C = corner(xC,yC,'C')
    xD = xC - 7.65
    yD = yC
    D = corner(xD, yD,'D')
    xE = xD
    yE = yD - 0.8               #pas sur 0.8...
    E = corner(xE, yE,'E')
    xF = xE - 9.0
    yF = yE
    F = corner(xF,yF,'F')
    xG = xF
    yG = yF + 4.37
    G = corner(xG,yG,'G')
    xH = xA
    yH = yG
    H = corner(xH,yH,'H') 
    walls = []
    # external walls
    print 'polyline'
    pl = tl.Polyline(B,C,E,F,G,H,A)
    extPolyline = pl
    walls.append(Wall(B,C))
    walls.append(Wall(C,E))
    walls.append(Wall(E,F))
    walls.append(Wall(F,G))
    walls.append(Wall(G,H))
    walls.append(Wall(H,A))
    
    # internal walls

    A1 = corner(xA+4.3,yA,'A1')
    H1 = corner(xA+4.3,yH,'H1')
    walls.append(Wall(A1,H1))
    A2 = corner(A1[0]+1.80,yA,'A2')
    K1 = corner(A2[0],A2[1]-2.,'K1')
    walls.append(Wall(A2,K1))
    A3 = corner(A2[0]+2.10,yA,'A3')
    K2 = corner(A3[0],A3[1]-2.,'K2')
    K3 = corner(K2[0],K2[1]-1.,'K3')
    walls.append(Wall(A3,K2))
    walls.append(Wall(K2,K3))
    walls.append(Wall(K1,K2))
    K4 = corner(K3[0]+0.84,K3[1],'K4')
    walls.append(Wall(K3,K4))

    K5 = corner(K4[0],K4[1]+1.66,'K5')
    walls.append(Wall(K4,K5))
    K6 = corner(K3[0],K5[1],'K6')

    walls.append(Wall(K5,K6))


    A4 = corner(A3[0] + 1.87,yA,'A4')

    A5 = corner(A4[0] + 3.2,yA,'A5')
    K7 = corner(A4[0],K4[1],'K7')
    K8 = corner(A5[0],K7[1],'K8')
    C1 = corner(K8[0],C[1],'C1')
    walls.append(Wall(A4,K7))
    walls.append(Wall(K4,K7))
    walls.append(Wall(A5,K8))
    walls.append(Wall(K7,K8))
    walls.append(Wall(K8,C1))
    L1 = corner(D[0],D[1]+2.2,'L1')
    L2 = corner(L1[0]-5.,L1[1],'L2')

    M1 = corner(L2[0],H1[1],'M1')
    walls.append(Wall(L2,M1))
    E1 = corner(L2[0],E[1],'E1')
    walls.append(Wall(D,L1))
    walls.append(Wall(L1,L2))
    walls.append(Wall(L2,E1))
    walls.append(Wall(G,H1))
    walls.append(Wall(H1,M1))
    walls.append(Wall(A,A1))
    walls.append(Wall(A1,A2))
    walls.append(Wall(A2,A3))
    walls.append(Wall(A3,A4))
    walls.append(Wall(A4,A5))
    walls.append(Wall(A5,B))

    plt.clf()

    for i,l in enumerate(extPolyline[1:]):
         xl = [extPolyline[i-1][0],extPolyline[i][0]]
         yl = [extPolyline[i-1][1],extPolyline[i][1]]
         plt.plot(xl,yl,'k',linewidth=2)
        

    #for w in walls:w.plot(plt)
    pts = points(walls)
    out = {}
    out['points'] = {}
    out['internal_walls'] = {}
    out['external_walls'] = {}
    for p in pts:
        out['points'][p.name] = p
    for w in walls:
        out['internal_walls'][w.name]=[w.p0.name, w.p1.name]
    with open('romarine.json', 'w') as fp:
        json.dump(out, fp, indent = 4)
    #scatter(pts,plt)
    plt.axis('equal')
    plt.show()

