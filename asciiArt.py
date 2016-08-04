def openFile(name):
    fl=open(name+'.asciiart').readlines()
    for l in fl:print l.replace('\n','')
