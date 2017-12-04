import numpy as np
import matplotlib.pyplot as plt

f = open('fort.99','r')
header = f.readline()
lines = f.readlines()[:]
f.close
keyword = np.array(header.split())
Nel  = keyword.size-3
dmax = np.zeros(Nel,dtype='double')
dmin = dmax + 1.E+99
iline = 1
for line in lines:
  data = line.split()
  if (data[0]=='Tg'): continue
  dat = np.array(data,dtype='double')
  iline = iline+1
  for i in range(0,Nel):
    val = dat[i+3]
    #print dmin[i],val,dmax[i]
    if (val>dmax[i]): dmax[i]=val
    if (val<dmin[i]): dmin[i]=val
                     
for i in range(0,Nel):
  print keyword[i+3],dmin[i],dmax[i]
  
   
