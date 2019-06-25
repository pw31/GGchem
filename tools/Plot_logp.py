import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
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
NPOINT = len(dat[0:])

bar   = 1.E+6                    # 1 bar in dyn/cm2 
Tg    = dat[:,0]                 # T [K]
nHtot = dat[:,1]                 # n<H> [cm-3]
lognH = np.log10(nHtot)          
press = dat[:,2]/bar             # p [bar]
lp    = np.log10(press)
pmin  = np.min(lp)
pmax  = np.max(lp)
Narg  = len(sys.argv)
if (Narg>1): pmin=float(sys.argv[1])
if (Narg>2): pmax=float(sys.argv[2])
iii   = np.where((lp>pmin) & (lp<pmax))[0]
Tmin  = np.min(Tg[iii])
Tmax  = np.max(Tg[iii])
Tmin  = Tmin*0.9
Tmax  = Tmax*1.1
nHmin = np.min(nHtot[iii])
nHmax = np.max(nHtot[iii])
nHmin = nHmin*0.9
nHmax = nHmax*1.1
if (nHmax>nHmin*5): 
  nHmin = nHmin/2.0
  nHmax = nHmax*2.0
#sep = 20
#if (Tmax-Tmin>1500): sep=100
#if (Tmax-Tmin>1000): sep=50
#if (Tmax-Tmin<600): sep=20
#if (Tmax-Tmin<400): sep=10
colo = ['blue','black','silver','red','darkorange','gold','darkorchid','aqua','cadetblue','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod','darkkhaki','pink','moccasin']
#'darkolivegreen','darkmagenta','aquamarine','coral','burlywood',
#'beige','darkorange','crimson','darkcyan','bisque'
Ncolor = len(colo)
colo = colo*10
styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
widt = [2]*Ncolor*10

#================== temperature-pressure structure ====================
fig,ax = plt.subplots()
plt.plot(lp,Tg,lw=4)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=16)
plt.ylabel(r'$T\ \mathrm{[K]}$',fontsize=16)
plt.xlim(pmin,pmax)
plt.ylim(Tmin,Tmax)
plt.tick_params(axis='both', labelsize=13)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== temperature-density structure ====================
fig,ax = plt.subplots()
plt.plot(lp,nHtot,lw=4)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=16)
plt.ylabel(r'$n_\mathrm{\langle H\rangle}\ \mathrm{[cm^{-3}]}$',fontsize=16)
plt.xlim(pmin,pmax)
plt.ylim(nHmin,nHmax)
if (nHmax>nHmin*5): plt.yscale('log')
plt.tick_params(axis='both', labelsize=13)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== the dust/gas mass ratio ====================
fig,ax = plt.subplots()
ind = np.where(keyword=='dust/gas')[0][0]
log10_dust_gas = dat[:,ind]
ymax = np.max(log10_dust_gas)
if (ymax>-10):
  plt.plot(lp,10**log10_dust_gas,lw=4)
  plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{dust/gas}$',fontsize=20)
  plt.xlim(pmin,pmax)
  plt.ylim(10**(ymax-8),10**ymax*3)
  plt.yscale('log')
  plt.tick_params(axis='both', labelsize=15)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  ax.yaxis.set_minor_locator(LogLocator(subs=[2,3,4,5,6,7,8,9]))
  #minorLocator = MultipleLocator(sep)
  #ax.xaxis.set_minor_locator(minorLocator)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

#================ the gas phase element abundances ===================
fig,ax = plt.subplots()
count = 0
ymax = -100.0
for i in range(4+NELEM+NMOLE+2*NDUST,4+NELEM+NMOLE+2*NDUST+NELEM,1):
  elm = keyword[i]
  element = elm[3:]
  yy = dat[:,i]               # log10 eps
  ymax=np.max([ymax,np.max(yy)])            
  if (np.max(yy)>-20):
    plt.plot(lp,yy,c=colo[count],ls=styl[count],lw=widt[count],label=element)
    count = count+1
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\log\,\epsilon_{\rm gas}$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(ymax-12,ymax+0.3)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
#minorLocator = MultipleLocator(1)
#ax.yaxis.set_minor_locator(minorLocator)
sz = np.min([11,1+195.0/count])
leg = plt.legend(loc='lower right',fontsize=sz,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== solid particle densities ===================
solids = []
smean = []
ymax = -100.0
for i in range(4+NELEM+NMOLE,4+NELEM+NMOLE+NDUST,1):
  solid = keyword[i]
  solids.append(solid[1:])
  smean.append(np.mean(dat[iii,i])) 
  ind = np.where(keyword == 'n'+solid[1:])[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  yy = dat[:,ind]               # log10 nsolid/n<H>
  ymax = np.max([ymax,np.max(yy[iii])])
  #print solid[1:],ind,np.max(yy[iii])
indices = np.argsort(smean)
if (ymax>-99):
  ymin = ymax-8
  print solids
  fig,ax = plt.subplots()
  count = 0
  for isolid in reversed(indices):
    solid = solids[isolid]
    ind = np.where(keyword == 'n'+solid)[0]
    if (np.size(ind) == 0): continue
    ind = ind[0]
    yy = dat[:,ind]               # log10 nsolid/n<H>
    ymax = np.max([ymax,np.max(yy[iii])])
    if (np.max(yy[iii])>ymin):
      plt.plot(lp[iii],yy[iii],c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  plt.title('condensates',fontsize=20)
  plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$',fontsize=20)
  #plt.xscale('log')
  plt.xlim(pmin,pmax)
  plt.ylim(ymin,ymax+0.3)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  #minorLocator = MultipleLocator(sep)
  #ax.xaxis.set_minor_locator(minorLocator)
  #minorLocator = MultipleLocator(1.0)
  #ax.yaxis.set_minor_locator(minorLocator)
  sz = np.min([11,1+195.0/count])
  leg = plt.legend(loc='lower right',fontsize=11,fancybox=True,
             handlelength=2.5,prop={'size':sz})
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

  for iT in range(0,NPOINT):
    iact = 0
    outp = ' '
    for i in range(4+NELEM+NMOLE+NDUST,4+NELEM+NMOLE+2*NDUST,1):
      #print keyword[i],dat[iT,i]
      if (dat[iT,i]>-200): 
        iact=iact+1
        outp=outp+' '+keyword[i][1:]
    print Tg[iT],iact,outp  

#==============================================================
fig,ax = plt.subplots()
count = 0
for isolid in reversed(indices):
  solid = solids[isolid]
  ind = np.where(keyword == 'S'+solid)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  #print solid,ind
  S = 10**dat[:,ind]              # S
  if (np.max(S[iii])>0.7):
    plt.plot(lp,S,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
    count = count + 1
#plt.title('supersaturation ratios',fontsize=20)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$S$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(0,1.05)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
sz = np.min([13,1+195.0/count])
leg = plt.legend(loc='lower right',fontsize=10,fancybox=True,
           handlelength=3,prop={'size':sz})
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== supersaturation ratios ===================
fig,ax = plt.subplots()
count = 0
for isolid in reversed(indices):
  solid = solids[isolid]
  #solid = solid[1:]
  ind = np.where(keyword == 'S'+solid)[0]
  if (np.size(ind) == 0): continue
  if ('[l]' in solid): continue
  ind = ind[0]
  logS = dat[:,ind]              # log10 S
  if (np.max(logS[iii])>-0.2):
    print solid
    plt.plot(lp,logS,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
    count = count + 1
plt.title('supersaturation ratios',fontsize=20)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ S$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(-10,10)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
sz = np.min([13,1+195.0/count])
col = 1
if (count>30): 
  sz = np.min([13,1+90.0/count*2])
  col = 2
leg = plt.legend(loc='upper right',fontsize=sz,fancybox=True,
                 handlelength=3,prop={'size':sz},ncol=col)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== some important molecules ====================
fig,ax = plt.subplots()
mols  = ['H2','H','N2','H2O','O2','CO','CO2','CH4','NH3','C2H2','el']
mols  = ['CO2','N2','SO2','COS','H2S','S2','CO','H2O','H2','HCL','HF','Ar','Ne','H2SO4']
mols  = np.array(mols)
ntot  = 0.0*nHtot
for i in range(3,4+NELEM+NMOLE): # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat[:,i]
lntot = np.log10(ntot)
count = 0
for i in range(3,4+NELEM+NMOLE): 
  mol = keyword[i]
  yy = dat[:,i]-lntot            # log10 nmol/ntot
  crit = -5
  ind = np.where(mols == mol)[0]
  if (np.size(ind)>0): crit=-7
  #print i,mol,ind,np.size(ind)
  if (np.max(yy[iii])>crit):
    plt.plot(lp,yy,c=colo[count],ls=styl[count],lw=widt[count],label=mol)
    count = count + 1
plt.title('important molecules',fontsize=20)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(-8,0.2)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
if (Tmax/Tmin>10):
  plt.xscale('log')
  #else:  
  #minorLocator = MultipleLocator(sep)
  #ax.xaxis.set_minor_locator(minorLocator)

#minorLocator = MultipleLocator(0.2)
#ax.yaxis.set_minor_locator(minorLocator)
leg = plt.legend(loc='lower left',fontsize=11,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== where are the elements? ================
ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','F','P','NI','MN','CR','ZN','ZR','RB','CU','B','BR','V','SR','W','el']
allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li',' ',' ','Ni','Mn','Cr','Zn','Zr','Rb','Cu',' ','Br',' ','Sr',' ','+']
exlist = [' He ',' Cl CL Ca CA Cr CR Co Cu CU ',' ',' Na NA Ni NI Ne NE',' ',' Si SI Sr SR ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' Fe FE ',' ',' ',' ',' ',' ',' ',' ',' ',' Br BR ',' ',' ',' ',' ',' ']
titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','fluorine','phosphorus','nickel','manganese','chromium','zinc','zirconium','rubidium','copper','boron','bromine','vanadium','strontium','tungston','charge carriers']
limits = [2,5,2.5,6,6,5,6,4,7,8,6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5]   
for i in range(0,30):
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
      next1 = molname[ind:ind+2]
      next2 = molname[ind-1:ind+1]
      #print keyword[mol],next1,str.find(ex,next1),len(next1)
      if (len(next1)==1 or str.find(ex,next1)==-1 or molname=='SIS'):
        if (next2!='MN' and next2!='ZN'):
          yy = dat[:,mol]                # log10 nmol [cm-3]
          yy = yy - lognH                # log10 nmol/n<H>
          nmax = np.max([nmax,np.max(yy[iii])])
          maxy = maxy + 10**yy
          if (molname=='el'): nmin = np.min([nmin,np.min(yy[iii])])
          mollist.append(mol)   
          abulist.append(np.mean(yy))
  if (nmax==-100): continue
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
      plt.plot(lp,yy,c=colo[count],ls=styl[count],lw=widt[count],label=molname)
      count = count + 1
  plt.title(titel,fontsize=20)
  plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{\langle H\rangle}$',fontsize=20)
  plt.xlim(pmin,pmax)
  plt.ylim(nmin,nmax+1)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  if (Tmax/Tmin>10):
    plt.xscale('log')
    #else:  
    #minorLocator = MultipleLocator(sep)
    #ax.xaxis.set_minor_locator(minorLocator)

  #minorLocator = MultipleLocator(1.0)
  #if (nmax-nmin>50): minorLocator = MultipleLocator(2.0)
  #if (nmax-nmin>100): minorLocator = MultipleLocator(5.0)
  #if (nmax-nmin>200): minorLocator = MultipleLocator(10.0)
  #ax.yaxis.set_minor_locator(minorLocator)
  sz = np.min([11,1+195.0/count])
  col = 1
  if (count>30): 
    sz = np.min([9,1+195.0/count*2])
    col = 2
  leg = plt.legend(loc='lower right',fontsize=10,fancybox=True,
             handlelength=3,prop={'size':sz},ncol=col)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()


pp.close()
print '... written output to ggchem.pdf.'


