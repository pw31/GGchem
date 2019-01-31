import numpy as np
import matplotlib.pyplot as plt
import sys

#==========================================================================
# read Sharp & Huebner data
#==========================================================================
SHnam   = np.loadtxt('SharpHuebner.dat',skiprows=0,usecols=[0],dtype='str')
SHcoeff = np.loadtxt('SharpHuebner.dat',skiprows=0,usecols=[1,2,3,4,5])

#==========================================================================
# read GGchem data
#==========================================================================
f = open('../../data/DustChem.dat','r')
lines = f.readlines()[:]
f.close
Nline = len(lines)
Ncond = len(SHnam)
Nfound = 0
f = open('DustChem_SH90.dat','w')
f.write('dust species\n')
f.write('============\n')
f.write('%d\n' %(Ncond))
for i in range(0,Ncond):
  cond = SHnam[i]+'[s]'
  lenc = len(cond)
  found = 0
  for j in range(0,Nline):
    line = lines[j]
    if (line[0:lenc]==cond):
      found = 1
      Nfound = Nfound+1
      print cond,j
      f.write('\n')
      correct = 0
      while True:
        line = lines[j]
        if (line=='\n'): break
        if (line[0]=='#'): correct=1
        if (correct==1): line = '#'+line[1:]
        if (correct==0): f.write(line)
        j = j+1
      c = SHcoeff[i,0:5]
      f.write('# Sharp & Huebner (1990):\n')
      f.write('  1 %12.5e %12.5e %12.5e %12.5e %12.5e\n' %(c[0],c[1],c[2],c[3],c[4]))

  if (found==0):
    print cond," not found."
    f.write('\n')
    f.write(cond+'\n')
    f.write('# Sharp & Huebner (1990):\n')
    f.write('  1 %12.5e %12.5e %12.5e %12.5e %12.5e\n' %(c[0],c[1],c[2],c[3],c[4]))
      
f.close
print Nfound," condensates found."
print Ncond-Nfound," condensates not found."
