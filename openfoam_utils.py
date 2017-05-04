import mesh
import subprocess
import sys
sys.path.append('./profiles/')
from splineProfileMultiParam import Profile
import readMSH
import os


def parse_dict(name):
    f=open(name,'r').readlines()
    header = []
    headerApp = False
    dico = ''
    i0 = 0
    headerClosed = True
    for i,l in enumerate(f):
        if l.startswith('/*'):
            headerClosed = False
            headerApp = True 
            for j in range(i,len(f)):
                if f[j].replace('\n','').endswith('*/'):
                    i0 = i+ j
                    header.append(f[j])
                    headerApp = False
                    break
                if headerApp:
                    header.append(f[j]) 
                else :
                    headerClosed = True
                    break

    for i in range(i0,len(f)):
        l = f[i]
        if l.startswith('/') or l.startswith('\\') or l.startswith('|'):
            pass
        else:
            dico+=l.replace('\n','')
    for l in header:print l.replace('\n','')
    print('//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')
    print(dico)
    
    print('//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')

    inDict = False
    istart = None
    iend = None
    separ = ['{',';','}',':']
    newDico=''
    for i,c in enumerate(dico):
        if c in separ:
            istart = i
            for j in range(i+1,len(dico)-1):
                if dico[j] in separ:
                    iend = j
                    substring = dico[istart:iend+1]
                    for s in separ:
                        substring = substring.replace(s,'')
                    if len(substring.split()) == 2:
                        sssplit = substring.split()
                        if (sssplit[0] not in separ and sssplit[1] not in separ):
                            dico = dico.replace(substring,sssplit[0]+': '+sssplit[1])
                    break
    dico = dico.replace(';}','},').replace(';',',').replace('{',':{').replace('\"','\'')
    dico = dico.replace(',}','}')
    for s in separ+[',']:
        dico = dico.replace(s,' '+s+' ')
    import json
    import re
    displit = dico.split()
    noquotes = []
    for s in displit:
        ire = re.search(r'^[-+]?[0-9]+$',s)
        if ire and (s not in noquotes) :
            noquotes.append(s)
        rre = re.search(r'[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?',s)
        if rre and (s not in noquotes) :
            noquotes.append(s)
    nd = ''
    for s in displit :
        if (s in noquotes) or (s in separ+[',']) :
            ns = (' '+s+' ')
        else:
            ns = (' "'+s+'" ')
        nd += ns
    nd = nd.replace(' ','')
    

    nd = nd.replace(',}','}')
    nd = nd.replace('}\"','},\"')


    print '----------------------------------------------'
    nd = '{'+nd+'}'
    print nd
    print '----------------------------------------------'
    print nd[260:]
    fjs = open(name+'.json','w')
    djs = json.loads(nd)
    
    json.dump(djs, fjs,indent=4)
    print djs



def create_2D_project(pl,**kwargs):
    subprocess.call(['mkdir','openfon_projects'])
    os.chdir('openfon_projects')
    subprocess.call(['mkdir',kwargs['name']])
    os.chdir(kwargs['name'])
    pl.box_2d(**kwargs) 
    subprocess.call(['mkdir','0'])
    subprocess.call(['mkdir','constant'])
    subprocess.call(['mkdir','system'])
    subprocess.call(['/home/fon/gmsh-2.16.0-Linux/bin/gmsh',kwargs['name']+'.geo','-2','-o',kwargs['name']+'.msh','>',kwargs['name']+'.mshlog'])
    print('MSH file successfully written')
    d = readMSH.read_msh_file(kwargs['name'])
    readMSH.write_fms_file(kwargs['name'],**d)
    print('FMS file successfully written')
    os.chdir('../../')
if __name__ == '__main__' :
    f = Profile(typ = 'fon',par = [0.82,0.21,0.13,0.04,0.029],npt = 100)
    pl = f.polyline()
    create_2D_project(pl,name = 'anais',facdom=10.)
    parse_dict('test_2D_cfmesh/system/meshDict')
