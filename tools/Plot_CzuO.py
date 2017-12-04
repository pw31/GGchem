import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 2.0
pp = PdfPages('ggchem.pdf')

file   = 'Static_Conc_CzuO_nocond.dat'
file   = 'Static_Conc.dat'
data   = open(file)
dummy  = data.readline()
dimens = data.readline()
dimens = np.array(dimens.split())
NELEM1 = int(dimens[0])
NMOLE1 = int(dimens[1])
NDUST1 = int(dimens[2])
NPOINT = int(dimens[3])
header = data.readline()
data.close()
dat1 = np.loadtxt(file,skiprows=3)
keyword = np.array(header.split())
iC = np.where('epsC'==keyword)[0][0]
iO = np.where('epsO'==keyword)[0][0]
CzuO1 = 10**dat1[:,iC]/10**dat1[:,iO]
ntot1   = 0.0*CzuO1
for i in range(3,4+NELEM1+NMOLE1):  # electrons, all atoms, ions and cations
  ntot1 = ntot1 + 10**dat1[:,i]
lntot1 = np.log10(ntot1)

file    = 'Static_Conc_CzuO_eqcond.dat'
data    = open(file)
dummy   = data.readline()
dimens  = data.readline()
dimens  = np.array(dimens.split())
NELEM2  = int(dimens[0])
NMOLE2  = int(dimens[1])
NDUST2  = int(dimens[2])
NPOINT  = int(dimens[3])
header  = data.readline()
data.close()
dat2    = np.loadtxt(file,skiprows=3)
keyword = np.array(header.split())
nHtot2  = dat2[:,1]                   # n<H> [cm-3]
lnHtot2 = np.log10(nHtot2)          
CzuO2   = 10**dat2[:,iC]/10**dat2[:,iO]
for i in range(0,NPOINT):
  CzuO2[i] = 0.3+1.1*i/(NPOINT-1)
ntot2   = 0.0*CzuO2
for i in range(3,4+NELEM1+NMOLE1):  # electrons, all atoms, ions and cations
  ntot2 = ntot2 + 10**dat2[:,i]
lntot2 = np.log10(ntot1)
iC = np.where('epsC'==keyword)[0][0]
iO = np.where('epsO'==keyword)[0][0]

bar  = 1.E+6                       # 1 bar in dyn/cm2 
col  = ['black','blue','silver','red','gold','darkorange','darkorchid','aqua','cadetblue','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod','darkkhaki','pink','moccasin']
Ncolor = len(col)
col  = col*10
styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
widt = [2]*Ncolor*10

#================== some important molecules ====================
fig,ax = plt.subplots(figsize=(8,5))
mols  = ['CO','H2O','HCN','CH4','C2H2']
mols  = np.array(mols)
count = 0
for mol in mols:
  i = np.where(mol==keyword)[0][0]
  yy = dat1[:,i]-lntot1            # log10 nmol/ntot
  plt.plot(CzuO1,yy,c=col[count],lw=3,label=mol)
  yy = dat2[:,i]-lntot2   
  plt.plot(CzuO2,yy,c=col[count],lw=3,ls='--')
  count = count + 1
plt.xlabel(r'$\mathrm{C/O}$',fontsize=23)
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=23)
plt.xlim(0.3,1.4)
plt.ylim(-13,-2.5)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=8, width=2.0, which='major')
plt.tick_params('both', length=4, width=1.5, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
#minorLocator = MultipleLocator(1)
#ax.yaxis.set_minor_locator(minorLocator)
leg = plt.legend(loc='lower right',fontsize=14,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

H2O = dat1[:,np.where('H2O'==keyword)[0][0]]
CH4 = dat1[:,np.where('CH4'==keyword)[0][0]]
index = np.where(np.abs(H2O-CH4)==np.min(np.abs(H2O-CH4)))[0][0]
print "transition C/O=",CzuO1[index]

#================ the gas phase element abundances ===================
fig,ax = plt.subplots(figsize=(8,6))
count = 0
for i in range(4+NELEM1+NMOLE1+2*NDUST1,4+NELEM1+NMOLE1+2*NDUST1+NELEM1,1):
  elm = keyword[i]
  element = elm[3:]
  yy = dat2[:,i]               # log10 eps
  if (np.max(yy)>-20):
    plt.plot(CzuO2,yy,c=col[count],ls=styl[count],lw=widt[count],label=element)
    count = count+1
plt.xlabel(r'$\mathrm{C/O}$',fontsize=20)
plt.ylabel(r'$\log\,\epsilon_{\rm gas}$',fontsize=20)
plt.xlim(0.3,1.4)
plt.ylim(-13,1)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
sz = np.min([11,1+195.0/count])
plt.legend(loc='lower right',fontsize=sz,fancybox=True)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== solid particle densities ===================
solids = []
smean = []
nmax = float(-100)
ymin = -8.0  #-7.6
ymax = -3.5
for i in range(4+NELEM2+NMOLE2,4+NELEM2+NMOLE2+NDUST2,1):
  solid = keyword[i]
  solids.append(solid[1:])
  smean.append(np.mean(dat2[:,i])) 
  ind = np.where(keyword == 'n'+solid[1:])[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  yy = dat2[:,ind]                 # log10 nsolid/n<H>
  nmax = np.max([nmax,np.max(yy)])
if (nmax>-99):
  print solids
  fig,ax = plt.subplots(figsize=(8,5))
  indices = np.argsort(smean)
  count = 0
  for isolid in reversed(indices):
    solid = solids[isolid]
    ind = np.where(keyword == 'n'+solid)[0]
    if (np.size(ind) == 0): continue
    ind = ind[0]
    yy = dat2[:,ind]               # log10 nsolid/n<H>
    yy = yy - lntot2 + lnHtot2
    nmax = np.max([nmax,np.max(yy)])
    if (np.max(yy)>ymin):
      plt.plot(CzuO2,yy,c=col[count],ls=styl[count],lw=3,label=solid+'[s]')
      count = count + 1
  plt.xlabel(r'$\mathrm{C/O}$',fontsize=23)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{tot}$',fontsize=23)
  plt.xlim(0.3,1.4)
  plt.ylim(ymin,ymax)
  #minorLocator = MultipleLocator(0.5)
  #ax.yaxis.set_minor_locator(minorLocator)
  minorLocator = MultipleLocator(1.0)
  ax.yaxis.set_major_locator(minorLocator)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=8, width=2.0, which='major')
  plt.tick_params('both', length=4, width=1.5, which='minor')
  leg = plt.legend(loc='lower right',fontsize=11,fancybox=True)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

  #for iT in range(0,NPOINT):
  #  iact = 0
  #  outp = ' '
  #  for i in range(4+NELEM+NMOLE+NDUST,4+NELEM+NMOLE+2*NDUST,1):
  #    #print keyword[i],dat[iT,i]
  #    if (dat[iT,i]>-200): 
  #      iact=iact+1
  #      outp=outp+' '+keyword[i][1:]
  #  print Tg[iT],iact,outp  

pp.close()
print '... written output to ggchem.pdf.'


