import triangleLib as tl
import json


def read_msh_file(name,**kwargs):
    d = {}
    lmsh = open(name,'r').readlines()
    for i,l in enumerate(lmsh) :
        if '$PhysicalNames' in l:
            d['physical']=[]
            N = int(lmsh[i+1])
            for j in range(N):
                il = i+j+2
                lpn = lmsh[il].split()
                number = int(lpn[1])
                name = str(lpn[2].replace('\"',''))
                d['physical'].append([number,name])
            i = i+N+2
            break
    i0 = i

    for i,l in enumerate(lmsh) :
        if '$Nodes' in l:
            d['nodes']=[]
            N = int(lmsh[i+1])
            for j in range(N):
                il = i+j+2
                lpn = lmsh[il].split()
                print('point index :'+str(int(lpn[0])))
                pt = str((float(lpn[1]),float(lpn[2]),float(lpn[3]))).replace(',',' ')
                d['nodes'].append(pt)
            i = i+N+2
            break
    i0 = i

    for i,l in enumerate(lmsh) :
        if '$Elements' in l:
            d['triangles']=[]
            N = int(lmsh[i+1])
            for j in range(N):
                il = i+j+2
                lpn = lmsh[il].split()
                if int(lpn[1])==2:
                    print '---------- LPN --------------'
                    print lpn
                    print '-----------------------------'
                    print('triangle :'+str(int(lpn[0])))
                    tags = int(lpn[2])
                    phys = int(lpn[1+tags])
                    print ('physical :'+str(phys))
                    p0 = int(lpn[2+tags+1])
                    p1 = int(lpn[2+tags+2])
                    p2 = int(lpn[2+tags+3])
                    print p0,p1,p2
                    for k in d['physical']:
                        if k[0]==phys:
                            print 'k ='+str(k)
                            print('triangle belongs to '+k[1])
                            pt = ((p0-1,p1-1,p2-1),k[0]-1)
                            d['triangles'].append(str(pt).replace(',',''))
                            break
            i = i+N+2
            break
    i0 = i
    print lmsh[i0]
    return d


def write_fms_file(name,**kwargs):
    f = open(name+('.fms'),'w')
    d=kwargs
    print '------------------------------------------'
    print '------------------ FMS PART --------------'
    print '------------------------------------------'
    print d.keys()
    f.write('// patch names and types\n')
    f.write(str(len(d['physical']))+'\n')
    f.write('(\n')
    for i,k in enumerate(d['physical']):
        f.write(k[1]+'\n')
        if k[1]=='cylinder':
            f.write('wall\n')
        else:
            f.write('patch\n')
        f.write('\n')
    f.write(')\n\n')
    f.write('// coordinates of surface points\n')
    f.write(str(len(d['nodes']))+'\n')
    f.write('(\n')
    for p in d['nodes']:
        f.write(str(p)+'\n')
    f.write(')\n\n')
    f.write('// list of triangles\n')
    f.write(str(len(d['triangles']))+'\n')
    f.write('(\n')
    for p in d['triangles']:
        f.write(str(p)+'\n')
    f.write(')\n\n')

        


if __name__=='__main__':
    name = 'cylinder_zero.msh'
    d = read_msh_file(name)
    write_fms_file(name.replace('.msh',''),**d)
