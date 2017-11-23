import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 2.0
pp = PdfPages('ggchem.pdf')

file    = 'Static_Conc.dat'
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
keyword = np.array(header.split())
Tg1     = dat1[:,0]                 # T [K]
ntot    = 0.0*Tg1
for i in range(3,4+NELEM1+NMOLE1):  # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat1[:,i]
lntot1 = np.log10(ntot)
iC = np.where('epsC'==keyword)[0][0]
iO = np.where('epsO'==keyword)[0][0]
CzuO = 10**dat1[:,iC]/10**dat1[:,iO]

bar   = 1.E+6                       # 1 bar in dyn/cm2 
col = ['black','blue','chartreuse','red','pink','darkgoldenrod','darkgray','darkmagenta','red','darkorange','darkorchid','aqua','cadetblue']
col2 = ['aquamarine','beige','darkolivegreen','bisque','burlywood','chartreuse','chocolate','coral','cornflowerblue','crimson','darkcyan','darkkhaki']


#================== some important molecules ====================
fig,ax = plt.subplots(figsize=(8,5))
mols  = ['CO','H2O','HCN','CH4','C2H2']
mols  = np.array(mols)
count = 0
for mol in mols:
  i = np.where(mol==keyword)[0][0]
  yy = dat1[:,i]-lntot1            # log10 nmol/ntot
  plt.plot(CzuO,yy,c=col[count],lw=3,label=mol)
  count = count + 1
plt.xlabel(r'$\mathrm{C/O}$',fontsize=23)
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=23)
plt.xlim(0.3,1.4)
plt.ylim(-12.5,-2.5)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=8, width=2.0, which='major')
plt.tick_params('both', length=4, width=1.5, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
#minorLocator = MultipleLocator(1)
#ax.yaxis.set_minor_locator(minorLocator)
plt.legend(loc='lower right',fontsize=14,fancybox=True)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

pp.close()
print '... written output to ggchem.pdf.'


