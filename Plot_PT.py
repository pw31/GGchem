import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
pp = PdfPages('ggchem.pdf')

file   = 'Static_Conc.dat'
data   = open(file)
dummy  = data.readline()
dimens = data.readline()
dimens = np.array(dimens.split())
NELEM  = int(dimens[0])
NMOLE  = int(dimens[1])
NDUST  = int(dimens[2])
NPOINT = int(dimens[3])
header = data.readline()
data.close()
dat = np.loadtxt(file,skiprows=3)
keyword = np.array(header.split())

bar   = 1.E+6                    # 1 bar in dyn/cm2 
Tg    = dat[:,0]                 # T [K]
nHtot = dat[:,1]                 # n<H> [cm-3]
lognH = np.log10(nHtot)          
press = dat[:,2]                 # p [dyn/cm2]
pmin  = np.min(press)/bar
pmax  = np.max(press)/bar
pmin  = pmin*0.9
pmax  = pmax*1.1
nHmin = np.min(nHtot)
nHmax = np.max(nHtot)
nHmin = nHmin*0.9
nHmax = nHmax*1.1
Tmin  = np.min(Tg)
Tmax  = np.max(Tg)
#if (Tmax>4*Tmin): Tmax=4*Tmin
#if (Tmin<Tmax/3): Tmin=Tmax/3
sep = 20
if (Tmax-Tmin>1500): sep=100
if (Tmax-Tmin<500): sep=10
Tmin  = Tmin*0.95
Tmax  = Tmax*1.1
styl  = ['-','-','-','-','-','-','-','--','--','--','--','--','--','--',':',':',':',':',':',':',':','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.']
widt  = [ 2 , 2 , 2 , 2 , 2 , 2 , 2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  ]

#================== temperature-pressure structure ====================
fig,ax = plt.subplots()
plt.plot(Tg,press/bar,lw=4)
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$p\ \mathrm{[bar]}$',fontsize=20)
plt.xlim(Tmin,Tmax)
plt.ylim(pmin,pmax)
plt.yscale('log')
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
minorLocator = MultipleLocator(sep)
ax.xaxis.set_minor_locator(minorLocator)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== setup ===================
iii = np.where((Tg>Tmin) & (Tg<Tmax))[0]
solids = []
smean = []
nmax = float(-100)
for i in range(4+NELEM+NMOLE,4+NELEM+NMOLE+NDUST,1):
  solid = keyword[i]
  solids.append(solid[1:])
  smean.append(np.mean(dat[iii,i])) 
  ind = np.where(keyword == 'n'+solid[1:])[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  yy = dat[:,ind]               # log10 nsolid/n<H>
  nmax = np.max([nmax,np.max(yy[iii])])
  #print solid[1:],ind,np.max(yy[iii])
if (nmax>-99):
  print solids
  fig,ax = plt.subplots()
  indices = np.argsort(smean)

#================== phase diagrams ===================
fig,ax = plt.subplots()
for isolid in reversed(indices):
  solid = solids[isolid]
  ind = np.where(keyword == 'S'+solid)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  S = 10**dat[:,ind]              # S
  Tphase = np.array(Tg)
  Pphase = np.array(press)/bar
  Pphase = np.log10(Pphase)
  Sphase = np.array(S) #.reshape(300,300)
#  plt.tricontourf(Tphase, Pphase, Sphase, 30)
#  plt.contourf(Sphase, 30)
  plt.tricontour(Tphase, Pphase, Sphase, [1.0],color='blue')
  for iliquid in reversed(indices):
    liquid = solids[iliquid]
    if (liquid != solid+'[l]'): continue       #include liquid
    ind = np.where(keyword == 'S'+liquid)[0]
    if (np.size(ind) == 0): continue
    ind = ind[0]
    S = 10**dat[:,ind]              # S
    Tphase = np.array(Tg)
    Pphase = np.array(press)/bar
    Pphase = np.log10(Pphase)
    Sphase = np.array(S) #.reshape(100,100)
    plt.tricontour(Tphase, Pphase, Sphase, [1.0],color='red')
  plt.title(solid,fontsize=20)
  plt.xlabel(r'$T\mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ P\mathrm{[bar]}$',fontsize=20)
  xmin = np.amin(Tphase)
  xmax = np.amax(Tphase)
  ymin = np.amin(Pphase)
  ymax = np.amax(Pphase)
  ax.set_xlim([xmin,xmax])
  ax.set_ylim([ymin,ymax])
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  minorLocator = MultipleLocator(sep)
  ax.xaxis.set_minor_locator(minorLocator)
#  cbr = plt.colorbar()
#  cbr.set_label(r'$supersaturation ratio$',fontsize=16)
#  plt.show()
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

pp.close()
print '... written output to ggchem.pdf.'


