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
NELEM  = 21
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
pmin  = pmin/1.2
pmax  = pmax*1.2
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

#================== solid particle densities ===================
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
  count = 0
  for isolid in reversed(indices):
    solid = solids[isolid]
    ind = np.where(keyword == 'n'+solid)[0]
    if (np.size(ind) == 0): continue
    ind = ind[0]
    yy = dat[:,ind]               # log10 nsolid/n<H>
    nmax = np.max([nmax,np.max(yy[iii])])
    if (np.max(yy[iii])>nmax-20):
      plt.plot(Tg,yy,ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  plt.title('condensates',fontsize=20)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$',fontsize=20)
  #plt.xscale('log')
  plt.xlim(Tmin,Tmax)
  plt.ylim(nmax-9,nmax+0.5)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  minorLocator = MultipleLocator(sep)
  ax.xaxis.set_minor_locator(minorLocator)
  plt.legend(loc='lower right',fontsize=9,fancybox=True)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

#================== supersaturation ratios ===================
  fig,ax = plt.subplots()
  count = 0
  for isolid in reversed(indices):
    solid = solids[isolid]
    ind = np.where(keyword == 'S'+solid)[0]
    if (np.size(ind) == 0): continue
    ind = ind[0]
    #print solid,ind
    logS = dat[:,ind]              # log10 S
    if (np.max(logS[iii])>-6):
      plt.plot(Tg,logS,ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  plt.title('supersaturation ratios',fontsize=20)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ S$',fontsize=20)
  #plt.xscale('log')
  plt.xlim(Tmin,Tmax)
  plt.ylim(-7,0.5)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  minorLocator = MultipleLocator(sep)
  ax.xaxis.set_minor_locator(minorLocator)
  plt.legend(loc='lower right',fontsize=6,fancybox=True)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

  fig,ax = plt.subplots()
  count = 0
  for isolid in reversed(indices):
    solid = solids[isolid]
    ind = np.where(keyword == 'S'+solid)[0]
    if (np.size(ind) == 0): continue
    ind = ind[0]
    #print solid,ind
    S = 10**dat[:,ind]              # S
    if (np.max(S[iii])>0.2):
      plt.plot(Tg,S,ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  plt.title('supersaturation ratios',fontsize=20)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$S$',fontsize=20)
  #plt.xscale('log')
  plt.xlim(Tmin,Tmax)
  plt.ylim(0,1.05)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  plt.legend(loc='lower right',fontsize=7,fancybox=True)
  minorLocator = MultipleLocator(sep)
  ax.xaxis.set_minor_locator(minorLocator)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

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
  Sphase = np.array(S)
#  plt.tricontourf(Tphase, Pphase, Sphase, 30)
  plt.tricontour(Tphase, Pphase, Sphase, [0.9],color='blue')
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
    Sphase = np.array(S)
    plt.tricontour(Tphase, Pphase, Sphase, [0.9],color='red')
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

#================== some important molecules ====================
fig,ax = plt.subplots()
mols  = ['H2','H','N2','H2O','O2','CO','CO2','CH4','NH3','C2H2','O3','He','el']
nmax  = np.float(0)
for mol in range(4,4+NELEM+NMOLE,1):
  yy = dat[:,mol]                # log10 nmol [cm-3]
  yy = yy - lognH                # log10 nmol/n<H>
  nmax = np.max([nmax,np.max(yy[iii])])
count = 0
for mol in mols:
  ind = np.where(keyword == mol)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  print mol,ind
  yy = dat[:,ind]                # log10 nmol [cm-3]
  yy = yy - lognH                # log10 nmol/n<H>
  if (np.max(yy)>nmax-6):
    plt.plot(Tg,yy,ls=styl[count],lw=widt[count],label=mol)
    count = count + 1
plt.title('important molecules',fontsize=20)
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{\langle H\rangle}$',fontsize=20)
plt.xscale('log')
plt.xlim(Tmin,Tmax)
plt.ylim(nmax-8,nmax+0.5)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
plt.legend(loc='lower right',fontsize=10,fancybox=True)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== where are the elements? ================
ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','el']
allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li','+']
exlist = [' He ',' Cl CL Ca CA Cr ',' ',' Na NA Ni ',' ',' Si SI ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','charge carriers']
limits = [2,5,3.5,6,6,5,6,4,7,8,6,6,6,6,7,5]   
for i in range(0,16):
  fig,ax = plt.subplots()
  el = ellist[i]
  al = allist[i]
  ex = exlist[i]
  limit = limits[i]
  titel = titels[i]
  print titel+" ..."
  nmax = np.float(-100)
  nmin = np.float(0)
  mollist = []
  abulist = []
  maxy = 0.0*dat[:,0]
  for mol in range(3,4+NELEM+NMOLE,1):
    molname = keyword[mol]
    ind = str.find(molname,el)
    if (ind < 0): 
      ind = str.find(molname,al)
    if (ind < 0 and el=='el'): 
      ind = str.find(molname,'-')
    if (ind >= 0):
      next = molname[ind:ind+2]
      #print mol,keyword[mol],next,str.find(ex,next),len(next)
      if (len(next)==1 or str.find(ex,next)==-1 or molname=='SIS'):
        yy = dat[:,mol]                # log10 nmol [cm-3]
        yy = yy - lognH                # log10 nmol/n<H>
        nmax = np.max([nmax,np.max(yy[iii])])
        maxy = maxy + 10**yy
        if (molname=='el'): nmin = np.min([nmin,np.min(yy[iii])])
        mollist.append(mol)   
        abulist.append(np.mean(yy))
  indices = np.argsort(abulist)
  count = 0
  maxy = np.log10(maxy)
  nmin = np.min([nmin,np.min(maxy[iii])-limit,nmax-12])
  for ind in reversed(indices):
    mol = mollist[ind]
    abu = abulist[ind]
    molname = keyword[mol]
    yy = dat[:,mol]                # log10 nmol [cm-3]
    yy = yy - lognH                # log10 nmol/n<H>
    if (np.max(yy[iii]-maxy[iii])>-limit or molname=='el'):
      print molname,abu
      plt.plot(Tg,yy,ls=styl[count],lw=widt[count],label=molname)
      count = count + 1
  plt.title(titel,fontsize=20)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{\langle H\rangle}$',fontsize=20)
  plt.xscale('log')
  plt.xlim(Tmin,Tmax)
  plt.ylim(nmin,nmax+1)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  #minorLocator = MultipleLocator(sep)
  #ax.xaxis.set_minor_locator(minorLocator)
  plt.legend(loc='lower right',fontsize=10,fancybox=True)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()


pp.close()
print '... written output to ggchem.pdf.'


