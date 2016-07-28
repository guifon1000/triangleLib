import numpy as np
import triangleLib as tl
import time
import sys




def createDict(name, addTriangulation=None,values=None):
    if addTriangulation:print 'triangulation added'
    f=open(name,'r')
    fl=f.readlines()
    d={}
    for l in fl:
        if "ZONE" in l:
            na=l.replace("ZONE","").replace("T","").replace("=","").replace('\n','').replace("\r","").replace(" ","")
            if values :
                if na in values:
                    if not d.has_key(na):d[na]=[]
    
    


    


    if addTriangulation and values :
        print 'experimental   : IN !'
        vtkd={}
 
        for k in d.keys():
            appendedPoints=[]
            triMat=[]
            vtkd[k+'_X']=np.zeros(len(addTriangulation))
            vtkd[k+'_Y']=np.zeros(len(addTriangulation))
            vtkd[k+'_Z']=np.zeros(len(addTriangulation))
            for i in range(len(fl)):
                itot=-1
                if str(k) in fl[i]:
                    for j in range(i+1,len(fl)):
                        if len(fl[j].split())!=0:
                            if "ZONE" in fl[j]:
                                break
                            
                            arPoint=[float(x) for x in fl[j].split()[0:3]]
                            values=[abs(float(x)) for x in fl[j].split()[3:6]]
                            pp=tl.Point()
                            pp.setPos(arPoint[0],arPoint[1],arPoint[2])
                            dd=-1.
                            intr=False
                            if len(appendedPoints)==0:intr=False
                            for l,ap in enumerate(appendedPoints[::-1]):
                                pa=tl.Point()
                                pa.setPos(ap[0],ap[1],ap[2])
                                itr=len(appendedPoints)-l
                                if tl.distance(pa,pp)==0.:
                                    intr=True
                                    itr=len(appendedPoints)-l
                                    break
                            
                            if not intr:
                                v=arPoint
                                v.append(abs(values[0]))
                                v.append(abs(values[1]))
                                #vs=tl.Vector()
                                #vs.set(values[0],values[1],values[2])
                                #vs.projectOnCoordinatesSystem(self.localSys[0],self.localSys[1],self.localSys[2])
                                minDis=1.e12
                                itt=None
                                for kt,t in enumerate(addTriangulation):
                                    ptl=tl.Point()
                                    ptl.setPos(arPoint[0],arPoint[1],arPoint[2])
                                    if tl.distance(t.cg,ptl)<minDis:
                                        minDis=tl.distance(t.cg,ptl)
                                        itt=kt

                                    
                                vtkd[k+'_X'][itt]=values[0]
                                vtkd[k+'_Y'][itt]=values[1]
                                vtkd[k+'_Z'][itt]=values[2]
                                appendedPoints.append(arPoint)
                                itot+=1
                                print (float(itot)/len(addTriangulation)*100.)," % done           \r",
                            else: 
                                vtkd[k+'_X'][itr]+=values[0]
                                vtkd[k+'_Y'][itr]+=values[1]
                                vtkd[k+'_Z'][itr]+=values[2]
        di2={}                      
        for k in d.keys():
            di2[k]=[[0.,0.,0.]]*len(addTriangulation)
            cvb=[]
            for j in range(len(addTriangulation)):
                t=addTriangulation[j]
                v=tl.Vector()
                v.set(vtkd[k+'_X'][j],vtkd[k+'_Y'][j],vtkd[k+'_Z'][j])
                vc=v.projectOnCoordinatesSystem(t.localSys[0],t.localSys[1],t.localSys[2])
                di2[k][j]=np.abs([vc.x,vc.y,vc.z])
        di3={}
        for k in di2.keys():
            di3[k+'_X']=np.array([di2[k][jj][0] for jj in range(len(addTriangulation))])
            di3[k+'_Y']=np.array([di2[k][jj][1] for jj in range(len(addTriangulation))])
            di3[k+'_Z']=np.array([di2[k][jj][2] for jj in range(len(addTriangulation))])
        print di3
        return di3
            
        
                    

def expressVectorOnTriangles(vtkd,addTriangulation):
    vtk2={}
    for i,t in enumerate(addTriangulation):
        for k in vtkd.keys():
            vec=tl.Vector()
            
            
            
        

    
def VTKdict(d,triangulation,values=None):
    vtkd={}
    print " filling the dictionary for VTK  "
    for k in d.keys():
        assign=False
        if values :
            if k in values:assign=True
        if 'Contraintes_Principales' in k and assign:
            t1=time.clock()
            kx='principalStress_X'
            ky='principalStress_Y'
        if 'Contraintes_Secondaires' in k and assign:
            t1=time.clock()
            kx='secondaryStress_X'
            ky='secondaryStress_Y'
        if 'Deformations_Principales' in k and assign:
            t1=time.clock()
            kx='principalStrain_X'
            ky='principalStrain_Y'
        if 'Deformations_Secondaires' in k and assign:
            t1=time.clock()
            kx='secondaryStrain_X'
            ky='secondaryStrain_Y'
        vx=[]
        vy=[]
        if assign:
            print 'VTK face data for '+str(k)+' !!' 
	    for iit in range(len(triangulation)):
		t=triangulation[iit]
		vinter=[]
		for vec in d[k]:
		    if vec[-1]==iit:
			vinter.append(vec)
		if len(vinter)==0:
		    vv=[0.]*7
		    vinter.append(vv)
		lastV=vinter[0]
		for l in vinter : 
		    lastV[3]=max(lastV[3],l[3])
		    lastV[4]=max(lastV[4],l[4])      
		vx.append(lastV[3])         
		vy.append(lastV[4])  
	    vtkd[kx]=np.array(vx)       
	    vtkd[ky]=np.array(vy)      
	    t2=time.clock()
	    print 'took '+str(t2-t1)+" seconds"
    return vtkd   

if __name__=='__main__':
    d=createDict("CHAMPS_PLIS.TEC",addTriangulation=1.)
