import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
pp = PdfPages('ggchem.pdf')

file    = 'Static_fast.dat'
data    = open(file)
dummy   = data.readline()
dimens  = data.readline()
dimens  = np.array(dimens.split())
NELEM1  = int(dimens[0])
NMOLE1  = int(dimens[1])
NDUST1  = int(dimens[2])
header  = data.readline()
data.close()
dat1 = np.loadtxt(file,skiprows=3)
keyword1 = np.array(header.split())
Tg1     = dat1[:,0]                 # T [K]
ntot    = 0.0*Tg1
for i in range(3,4+NELEM1+NMOLE1):  # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat1[:,i]
lntot1 = np.log10(ntot)

file    = 'Static_gas.dat'
data    = open(file)
dummy   = data.readline()
dimens  = data.readline()
dimens  = np.array(dimens.split())
NELEM2  = int(dimens[0])
NMOLE2  = int(dimens[1])
NDUST2  = int(dimens[2])
header  = data.readline()
data.close()
dat2 = np.loadtxt(file,skiprows=3)
keyword2 = np.array(header.split())
Tg2     = dat2[:,0]                 # T [K]
ntot    = 0.0*Tg2
for i in range(3,4+NELEM2+NMOLE2):  # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat2[:,i]
lntot2 = np.log10(ntot)

file    = 'Static_cond.dat'
data    = open(file)
dummy   = data.readline()
dimens  = data.readline()
dimens  = np.array(dimens.split())
NELEM3  = int(dimens[0])
NMOLE3  = int(dimens[1])
NDUST3  = int(dimens[2])
header  = data.readline()
data.close()
dat3 = np.loadtxt(file,skiprows=3)
keyword3 = np.array(header.split())
Tg3     = dat3[:,0]                 # T [K]
ntot    = 0.0*Tg3
for i in range(3,4+NELEM3+NMOLE3):  # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat3[:,i]
lntot3 = np.log10(ntot)

bar   = 1.E+6                       # 1 bar in dyn/cm2 
Tmin  = 500
Tmax  = 3000
sep   = 100
col = ['darkgoldenrod','darkgray','darkgreen','darkmagenta','red','darkorange','darkorchid','aqua','cadetblue']
col2 = ['aquamarine','beige','darkolivegreen','bisque','burlywood','chartreuse','chocolate','coral','cornflowerblue','crimson','darkcyan','darkkhaki']


#================== some important molecules ====================
fig,ax = plt.subplots()
mols  = ['CO','CO2','CH4','N2','NH3','HCN','C2H2','C2H4','H2O']
mols  = np.array(mols)
count = 0
for mol in mols:
  i = np.where(mol==keyword1)[0][0]
  yy = dat1[:,i]-lntot1            # log10 nmol/ntot
  plt.plot(Tg1,yy,c=col[count],lw=3,label=mol)
  i = np.where(mol==keyword2)[0][0]
  yy = dat2[:,i]-lntot2            # log10 nmol/ntot
  plt.plot(Tg2,yy,c=col[count],lw=2,ls='--')
  i = np.where(mol==keyword3)[0][0]
  yy = dat3[:,i]-lntot3            # log10 nmol/ntot
  plt.plot(Tg3,yy,c=col[count],lw=2,ls=':')
  count = count + 1
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=20)
plt.xlim(Tmin,Tmax)
plt.ylim(-15,-2)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
minorLocator = MultipleLocator(sep)
ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(1)
ax.yaxis.set_minor_locator(minorLocator)
plt.legend(loc='lower right',fontsize=11,fancybox=True)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

pp.close()
print '... written output to ggchem.pdf.'


