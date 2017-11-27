import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter
plt.rcParams['axes.linewidth'] = 1.5

single_figures = 0

# units are cgs
bar = 1.E+6     # 1 bar in dyn/cm2
atm = 1.013E+6  # 1 atm in dyn/cm2
R   = 8.3144598 # J/mol/K
Rcal= 1.987     # cal/mol

f = open('../data/dispol_BarklemCollet.dat','r')
header = f.readline()
lines = f.readlines()[:]
f.close
NBC = int(header.split()[0])
BCname = np.empty(NBC, dtype="S12")
BCdat = np.zeros([NBC,5])
for i in range(0,NBC*2,2):
  #print lines[i]
  BCname[i/2] = lines[i].split()[0]
  tmp = lines[i+1].split()[1:6]
  for j in range(0,5):
    BCdat[i/2,j] = float(tmp[j])
print "BC species",BCname

f = open('../data/dispol_GGchem.dat','r')
header = f.readline()
lines = f.readlines()[:]
f.close
NGG = int(header)
GGname = np.empty(NGG, dtype="S12")
GGdat  = np.zeros([NGG,5])
GGcode = np.zeros(NGG,dtype='int')
for i in range(0,NGG*2,2):
  #print lines[i]
  GGname[i/2] = lines[i].split()[0]
  GGcode[i/2] = int(lines[i+1].split()[0])
  dat = lines[i+1].split()[1:6]
  for j in range(0,5):
    GGdat[i/2,j] = float(dat[j])
print "GG species",GGname

f = open('../data/dispol_StockKitzmann.dat','r')
header = f.readline()
lines = f.readlines()[:]
f.close
NSK = int(header)
SKname = np.empty(NSK, dtype="S12")
SKdat = np.zeros([NSK,5])
for i in range(0,NSK*2,2):
  #print lines[i]
  SKname[i/2] = lines[i].split()[0]
  dat = lines[i+1].split()[1:6]
  for j in range(0,5):
    SKdat[i/2,j] = float(dat[j])
print "SK species",SKname

f = open('../data/dispol_Tsuji.dat','r')
header = f.readline()
lines = f.readlines()[:]
f.close
NT = int(header)
Tname = np.empty(NT, dtype="S12")
Tdat = np.zeros([NT,5])
for i in range(0,NT*2,2):
  #print lines[i]
  Tname[i/2] = lines[i].split()[0]
  dat = lines[i+1].split()[1:6]
  for j in range(0,5):
    Tdat[i/2,j] = float(dat[j])
print "Tsuji species",Tname

f = open('../data/dispol_SharpHuebner.dat','r')
header = f.readline()
lines = f.readlines()[:]
f.close
NSH = int(header)
SHname = np.empty(NSH, dtype="S12")
SHdat = np.zeros([NSH,5])
SHcode = np.zeros(NSH,dtype='int')
for i in range(0,NSH*2,2):
  #print lines[i]
  SHname[i/2] = lines[i].split()[0]
  dat = lines[i+1].split()[1:6]
  SHcode[i/2] = lines[i+1].split()[0]
  #print i,SHname[i/2],dat
  for j in range(0,5):
    SHdat[i/2,j] = float(dat[j])
print "SH species",SHname
print
name = raw_input("which species? ")
Natom = int(raw_input("Natom=? "))

Tmin = 100.0
Tmax = 6500.0
T  = np.arange(Tmin,Tmax,1.0)
T2 = np.arange(400.0,Tmax,1.0)

#================== Fig 1 ===========================
if (single_figures==0): pp=PdfPages('kp.pdf')
if (single_figures==1): pp=PdfPages('kp_'+name+'_kp.pdf')
plt.figure(figsize=(4,4))

kpmean = 0*T
Nmean = 0
ln10 = np.log(10)
ind1 = np.where(name.upper()==GGname)[0]
if (len(ind1)>0):
  if (GGcode[ind1[0]]==1):
    a     = GGdat[ind1[0],0:5]
    Th    = 5040.0/T
    lnkp1 = a[0] + a[1]*Th + a[2]*Th**2 + a[3]*Th**3 + a[4]*Th**4
    #print "GGchem",GGname[ind1[0]],a
    plt.plot(T,lnkp1/ln10,c='black',lw=2.5,label='old GGchem')
    print "GGchem",GGname[ind1[0]],a
    kpmean = kpmean+lnkp1
    Nmean = Nmean+1
  else:
    ind1=''  

ind2 = np.where(name==SKname)[0]
if (len(ind2)>0):
  a     = SKdat[ind2[0],0:5]
  dGRT  = a[0]/T + a[1]*np.log(T) + a[2] + a[3]*T + a[4]*T**2
  lnkp2 = dGRT + (1-Natom)*np.log(bar)
  print "S&K",SKname[ind2[0]],a
  plt.plot(T,lnkp2/ln10,c='green',lw=3.5,label='Stock')
  kpmean = kpmean+lnkp2
  Nmean = Nmean+1

ind4 = np.where(name==SHname)[0]
if (len(ind4)>0):
  if (SHcode[ind4[0]]==3):
    a     = SHdat[ind4[0],0:5]
    dG    = a[0]/T + a[1] + a[2]*T + a[3]*T**2 + a[4]*T**3
    lnkp4 = -dG/(Rcal*T) + (1-Natom)*np.log(atm)
    print 'S&H',SHname[ind4[0]],a
    plt.plot(T,lnkp4/ln10,c='orange',lw=1.5,label='Sharp & Huebner')
    kpmean = kpmean+lnkp4
    Nmean = Nmean+1
  else:
    ind4 = ''
    
ind5 = np.where(name==Tname)[0]
if (len(ind5)>0):
  a     = Tdat[ind5[0],0:5]
  Th    = 5040.0/T2
  lnkp5 = -(a[0] + a[1]*Th + a[2]*Th**2 + a[3]*Th**3 + a[4]*Th**4)*np.log(10)
  Th    = 5040.0/T
  lnkp51= -(a[0] + a[1]*Th + a[2]*Th**2 + a[3]*Th**3 + a[4]*Th**4)*np.log(10)
  print 'Tsuji',Tname[ind5[0]],a
  plt.plot(T2,lnkp5/ln10,c='magenta',lw=1.5,label='Tsuji')

ind3 = np.where(name==BCname)[0]
if (len(ind3)>0):
  a     = BCdat[ind3[0],0:5]
  dGRT  = a[0]/T + a[1]*np.log(T) + a[2] + a[3]*T + a[4]*T**2
  lnkp3 = dGRT + (1-Natom)*np.log(bar)
  print 'BC',BCname[ind3[0]],a
  plt.plot(T,lnkp3/ln10,c='blue',ls='--',lw=1.5,label='Barklem & Collet')
  kpmean = kpmean+lnkp3
  Nmean = Nmean+1

kpmean = kpmean/Nmean
plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=15)
plt.ylabel(r'$\log_{10} k_p \mathrm{[cgs]}$',fontsize=15)
plt.title(name)
plt.xlim(Tmin,Tmax)
plt.xscale('log')
plt.subplots_adjust(left=0.21, right=0.94, top=0.94, bottom=0.14)
plt.tick_params(axis='both', which='major', length=6,width=1.5)
plt.tick_params(axis='both', which='minor', length=4,width=1)
plt.legend(loc='upper right',fontsize=10)
plt.savefig(pp,format='pdf')

#================== Fig 2 ===========================
if (single_figures==0): plt.clf()
if (single_figures==1): pp.close()
if (single_figures==1): pp=PdfPages('kp_'+name+'_kperr.pdf')
plt.figure(figsize=(4,4))

ymin = 0.0
ymax = 0.0
iT = np.where(T>300)[0]
iT1= np.where(T>600)[0]
iT2= np.where(T>1000)[0]
print iT
if (len(ind1)>0):
  plt.plot(T,(lnkp1-kpmean)/ln10,c='black',lw=2.5,label='old GGchem')
  ymin = np.min([ymin,np.min(lnkp1[iT]-kpmean[iT])])
  ymax = np.max([ymax,np.max(lnkp1[iT]-kpmean[iT])])
  print "GG",ymin,ymax

if (len(ind2)>0):
  plt.plot(T,(lnkp2-kpmean)/ln10,c='green',lw=3.5,label='Stock')
  ymin = np.min([ymin,np.min(lnkp2[iT]-kpmean[iT])])
  ymax = np.max([ymax,np.max(lnkp2[iT]-kpmean[iT])])
  print "SK",ymin,ymax

if (len(ind4)>0):
  plt.plot(T,(lnkp4-kpmean)/ln10,c='orange',lw=1.5,label='Sharp & Huebner')
  ymin = np.min([ymin,np.min(lnkp4[iT1]-kpmean[iT1])])
  ymax = np.max([ymax,np.max(lnkp4[iT1]-kpmean[iT1])])
  print "SH",ymin,ymax

if (len(ind5)>0):
  plt.plot(T,(lnkp51-kpmean)/ln10,c='magenta',lw=1.5,label='Tsuji')
  ymin = np.min([ymin,np.min(lnkp51[iT2]-kpmean[iT2])])
  ymax = np.max([ymax,np.max(lnkp51[iT2]-kpmean[iT2])])
  print "Tsu",ymin,ymax

if (len(ind3)>0):
  plt.plot(T,(lnkp3-kpmean)/ln10,c='blue',ls='--',lw=1.5,label='Barklem & Collet')
  ymin = np.min([ymin,np.min(lnkp3[iT]-kpmean[iT])])
  ymax = np.max([ymax,np.max(lnkp3[iT]-kpmean[iT])])
  print "BC",ymin,ymax

plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=15)
plt.ylabel(r'$\log_{10} k_p - \langle\log_{10} k_p\rangle \mathrm{[cgs]}$',fontsize=15)
plt.title(name)
plt.xlim(Tmin,Tmax)
#ymax = ymax+0.5*(ymax-ymin)
#dy = ymax-ymin
#plt.ylim((ymin-dy)/ln10,(ymax+dy)/ln10)
plt.ylim(-4,+4)
plt.xscale('log')
plt.subplots_adjust(left=0.19, right=0.94, top=0.94, bottom=0.14)
plt.tick_params(axis='both', which='major', length=6,width=1.5)
plt.tick_params(axis='both', which='minor', length=4,width=1)
#plt.legend(loc='upper right',fontsize=10)
plt.savefig(pp,format='pdf')

#================== Fig 3 ===========================
if (single_figures==0): plt.clf()
if (single_figures==1): pp.close()
if (single_figures==1): pp=PdfPages('kp_'+name+'_dG.pdf')

Gmean = 0*T
Nmean = 0
fig,ax = plt.subplots(figsize=(4,4))
if (len(ind1)>0):
  dG  = ((1-Natom)*np.log(bar) - lnkp1)*R*T/1000
  Gmean = Gmean+dG
  Nmean = Nmean+1
  plt.plot(T,dG,c='black',lw=2.5,label='old GGchem')

if (len(ind2)>0):
  dG  = ((1-Natom)*np.log(bar) - lnkp2)*R*T/1000
  Gmean = Gmean+dG
  Nmean = Nmean+1
  plt.plot(T,dG,c='green',lw=3.5,label='Stock')

if (len(ind4)>0):
  dG = ((1-Natom)*np.log(bar) - lnkp4)*R*T/1000
  #Gmean = Gmean+dG
  #Nmean = Nmean+1
  plt.plot(T,dG,c='orange',lw=1.5,label='Sharp & Huebner')

if (len(ind5)>0):
  dG2= ((1-Natom)*np.log(bar) - lnkp5)*R*T2/1000
  #Gmean = Gmean+dG2
  #Nmean = Nmean+1
  plt.plot(T2,dG2,c='magenta',lw=1.5,label='Tsuji')

if (len(ind3)>0):
  dG  = ((1-Natom)*np.log(bar) - lnkp3)*R*T/1000
  Gmean = Gmean+dG
  Nmean = Nmean+1
  plt.plot(T,dG,c='blue',ls='--',lw=1.5,label='Barklem & Collet')

Gmean = Gmean/Nmean
plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=15)
plt.ylabel(r'$\Delta G_{\rm f}^\theta \mathrm{[kJ/mol]}$',fontsize=15)
plt.title(name)
plt.subplots_adjust(left=0.22, right=0.94, top=0.94, bottom=0.14)
plt.legend(loc='upper left',fontsize=10)
plt.tick_params(axis='both', which='major', length=6,width=1.5)
plt.tick_params(axis='both', which='minor', length=4,width=1)
plt.xlim(0,Tmax)
minorLocator = MultipleLocator(500.0)
ax.xaxis.set_minor_locator(minorLocator)
#minorLocator = MultipleLocator(50.0)
#ax.yaxis.set_minor_locator(minorLocator)
plt.savefig(pp,format='pdf')

#================== Fig 4 ===========================
if (single_figures==0): plt.clf()
if (single_figures==1): pp.close()
if (single_figures==1): pp=PdfPages('kp_'+name+'_dGerr.pdf')
fig,ax = plt.subplots(figsize=(4,4))

ymin = 0.0
ymax = 0.0
if (len(ind1)>0):
  dG = ((1-Natom)*np.log(bar) - lnkp1)*R*T/1000
  plt.plot(T,dG-Gmean,c='black',lw=2.5,label='old GGchem')
  ymin = np.min([ymin,np.min(dG[iT]-Gmean[iT])])
  ymax = np.max([ymax,np.max(dG[iT]-Gmean[iT])])
  print "GG",ymin,ymax

if (len(ind2)>0):
  dG = ((1-Natom)*np.log(bar) - lnkp2)*R*T/1000
  plt.plot(T,dG-Gmean,c='green',lw=3.5,label='Stock')
  ymin = np.min([ymin,np.min(dG[iT]-Gmean[iT])])
  ymax = np.max([ymax,np.max(dG[iT]-Gmean[iT])])
  print "SK",ymin,ymax

if (len(ind4)>0):
  dG = ((1-Natom)*np.log(bar) - lnkp4)*R*T/1000
  plt.plot(T,dG-Gmean,c='orange',lw=1.5,label='Sharp & Huebner')
  ymin = np.min([ymin,np.min(dG[iT1]-Gmean[iT1])])
  ymax = np.max([ymax,np.max(dG[iT1]-Gmean[iT1])])
  print "SH",ymin,ymax

if (len(ind5)>0):
  dG = ((1-Natom)*np.log(bar) - lnkp51)*R*T/1000
  plt.plot(T,dG-Gmean,c='magenta',lw=1.5,label='Tsuji')
  ymin = np.min([ymin,np.min(dG[iT2]-Gmean[iT2])])
  ymax = np.max([ymax,np.max(dG[iT2]-Gmean[iT2])])
  print "Tsu",ymin,ymax

if (len(ind3)>0):
  dG = ((1-Natom)*np.log(bar) - lnkp3)*R*T/1000
  plt.plot(T,dG-Gmean,c='blue',ls='--',lw=1.5,label='Barklem & Collet')
  ymin = np.min([ymin,np.min(dG[iT]-Gmean[iT])])
  ymax = np.max([ymax,np.max(dG[iT]-Gmean[iT])])
  print "BK",ymin,ymax

plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=15)
plt.ylabel(r'$\Delta G_{\rm f}^\theta - \langle\Delta G_{\rm f}^\theta\rangle \mathrm{[kJ/mol]}$',fontsize=15)
plt.title(name)
plt.xlim(0,Tmax)
#plt.xscale('log')
#ymax = ymax+0.5*(ymax-ymin)
#dy = ymax-ymin
#plt.ylim(ymin-dy,ymax+dy)
plt.ylim(-30,+30)
plt.subplots_adjust(left=0.19, right=0.94, top=0.94, bottom=0.14)
#plt.legend(loc='upper left',fontsize=10)
plt.tick_params(axis='both', which='major', length=6,width=1.5)
plt.tick_params(axis='both', which='minor', length=4,width=1)
minorLocator = MultipleLocator(500.0)
ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(5.0)
ax.yaxis.set_minor_locator(minorLocator)
plt.savefig(pp,format='pdf')

pp.close()
print ' '
print 'written output to kp.pdf.'

  
   
