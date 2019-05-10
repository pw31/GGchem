import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
pp = PdfPages('ggchem.pdf')
import sys

mmHg = 1.3328E+3   # 1 mmHg in dyn/cm2
def yaws(T,A,B,C,D,E):
    #print A,B,C,D,E
    #print T
    val = 10**(A + B/T + C*np.log10(T) + D*T + E*T**2)
    return val*mmHg

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
press = dat[:,2]                 # p [dyn/cm2]
Tmin  = np.min(Tg)
Tmax  = np.max(Tg)
#if (Tmax>4*Tmin): Tmax=4*Tmin
#if (Tmin<Tmax/3): Tmin=Tmax/3
#Tmax  = 2800
#Tmin  = 2000
Narg  = len(sys.argv)
if (Narg>1): Tmin=float(sys.argv[1])
if (Narg>2): Tmax=float(sys.argv[2])
iii   = np.where((Tg>Tmin) & (Tg<Tmax))[0]
pmin  = np.min(press[iii])/bar
pmax  = np.max(press[iii])/bar
pmin  = pmin*0.9
pmax  = pmax*1.1
if (pmax>pmin*5): 
  pmin = pmin/2.0
  pmax = pmax*2.0
nHmin = np.min(nHtot[iii])
nHmax = np.max(nHtot[iii])
nHmin = nHmin*0.9
nHmax = nHmax*1.1
if (nHmax>nHmin*5): 
  nHmin = nHmin/2.0
  nHmax = nHmax*2.0
sep = 20
if (Tmax-Tmin>1500): sep=100
if (Tmax-Tmin>1000): sep=50
if (Tmax-Tmin<600): sep=20
if (Tmax-Tmin<400): sep=10
#Tmin  = Tmin*0.95
#Tmax  = Tmax*1.1
colo = ['blue','black','silver','red','darkorange','gold','darkorchid','aqua','cadetblue','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod','darkkhaki','pink','moccasin']
#'darkolivegreen','darkmagenta','aquamarine','coral','burlywood',
#'beige','darkorange','crimson','darkcyan','bisque'
Ncolor = len(colo)
colo = colo*10
styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
widt = [2]*Ncolor*10

#================== temperature-pressure structure ====================
fig,ax = plt.subplots()
plt.plot(Tg,press/bar,lw=4)
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$p\ \mathrm{[bar]}$',fontsize=20)
plt.xlim(Tmin,Tmax)
plt.ylim(pmin,pmax)
if (pmax>pmin*5): plt.yscale('log')
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
minorLocator = MultipleLocator(sep)
ax.xaxis.set_minor_locator(minorLocator)
cAl2O3 = [  14.1611, -2.8238E+04, -7.3843E-01, -3.7413E-07,  2.2086E-11, 2421, 3253 ]
cSiO2  = [-378.5210,  6.5473E+03,  1.3150E+02, -3.5774E-02,  3.4220E-06, 1883, 2503 ]
cH2O   = [  29.8605, -3.1522E+03, -7.3037E+00,  2.4247E-09,  1.8090E-06,  273,  647 ]
cFe    = [  11.5549, -1.9538E+04, -6.2549E-01, -2.7182E-09,  1.9086E-13, 1808, 3008 ]
cNH3   = [  37.1575, -2.0277E+03, -1.1601E+01,  7.4625E-03, -9.5811E-12,  195,  406 ]
cCO2   = [  35.0187, -1.5119E+03, -1.1335E+01,  9.3383E-03,  7.7626E-10,  216,  305 ]
#cc = cCO2
#Tyaws1 = cc[5]
#Tyaws2 = cc[6]
#Tyaws  = np.arange(Tyaws1, Tyaws2, 0.1)
#pvap   = yaws(Tyaws,*cc[0:5])
#Tvap = np.arange(Tmin, Tmax, 0.1)
#pvap = 10.0**(7.8974 - 2409.4/Tvap)
#plt.plot(Tvap,pvap,c='red',lw=2,ls='--')
#fmt=ScalarFormatter(useOffset=False)
#fmt.set_scientific(False)
#ax.yaxis.set_major_formatter(fmt)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== temperature-density structure ====================
fig,ax = plt.subplots()
plt.plot(Tg,nHtot,lw=4)
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$n_\mathrm{\langle H\rangle}\ \mathrm{[cm^{-3}]}$',fontsize=20)
plt.xlim(Tmin,Tmax)
plt.ylim(nHmin,nHmax)
if (nHmax>nHmin*5): plt.yscale('log')
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
minorLocator = MultipleLocator(sep)
ax.xaxis.set_minor_locator(minorLocator)
#fmt=ScalarFormatter(useOffset=False)
#fmt.set_scientific(False)
#ax.yaxis.set_major_formatter(fmt)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== the dust/gas mass ratio ====================
fig,ax = plt.subplots()
ind = np.where(keyword=='dust/gas')[0][0]
log10_dust_gas = dat[:,ind]
ymax = np.max(log10_dust_gas)
if (ymax>-10):
  plt.plot(Tg,10**log10_dust_gas,lw=4)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{dust/gas}$',fontsize=20)
  plt.xlim(Tmin,Tmax)
  plt.ylim(10**(ymax-8),10**ymax*3)
  plt.yscale('log')
  plt.tick_params(axis='both', labelsize=15)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  ax.yaxis.set_minor_locator(LogLocator(subs=[2,3,4,5,6,7,8,9]))
  minorLocator = MultipleLocator(sep)
  ax.xaxis.set_minor_locator(minorLocator)
  #fmt=ScalarFormatter(useOffset=False)
  #fmt.set_scientific(False)
  #ax.yaxis.set_major_formatter(fmt)
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
    plt.plot(Tg,yy,c=colo[count],ls=styl[count],lw=widt[count],label=element)
    count = count+1
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$\log\,\epsilon_{\rm gas}$',fontsize=20)
plt.xlim(Tmin,Tmax)
plt.ylim(ymax-12,ymax+0.3)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
minorLocator = MultipleLocator(sep)
ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(1)
ax.yaxis.set_minor_locator(minorLocator)
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
  ymin = -20
if (ymax>-99):
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
    ymax = np.max([ymax,np.max(yy[iii])])
    if (np.max(yy[iii])>-99): print solid,ind,np.max(yy[iii])
    if (np.max(yy[iii])>ymin):
      plt.plot(Tg[iii],yy[iii],c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  #plt.title('condensates',fontsize=20)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$',fontsize=20)
  #plt.xscale('log')
  plt.xlim(Tmin,Tmax)
  plt.ylim(ymax-10,ymax+0.3)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  minorLocator = MultipleLocator(sep)
  ax.xaxis.set_minor_locator(minorLocator)
  minorLocator = MultipleLocator(1.0)
  ax.yaxis.set_minor_locator(minorLocator)
  sz = np.min([11,1+195.0/count])
  #sz = sz+1
  leg = plt.legend(loc='lower right',fontsize=11,fancybox=True,
             handlelength=2.5,prop={'size':sz})
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  #plt.show()
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
      plt.plot(Tg,logS,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
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
  sz = np.min([13,1+195.0/count])
  col = 1
  if (count>30): 
    sz = np.min([13,1+195.0/count*2])
    col = 2
  leg = plt.legend(loc='lower right',fontsize=10,fancybox=True,
             handlelength=3,prop={'size':sz},ncol=col)
  leg.get_frame().set_alpha(0.7)
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
    if (np.max(S[iii])>0.7):
      plt.plot(Tg,S,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  #plt.title('supersaturation ratios',fontsize=20)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$S$',fontsize=20)
  #plt.xscale('log')
  plt.xlim(Tmin,Tmax)
  plt.ylim(0,1.05)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  sz = np.min([13,1+195.0/count])
  #sz = sz+3
  leg = plt.legend(loc='lower right',fontsize=10,fancybox=True,
             handlelength=3,prop={'size':sz})
  leg.get_frame().set_alpha(0.7)
  minorLocator = MultipleLocator(sep)
  ax.xaxis.set_minor_locator(minorLocator)
  plt.tight_layout()
  #plt.show()
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
  crit = -2.5
  crit = -5
  ind = np.where(mols == mol)[0]
  if (np.size(ind)>0): crit=-5
  if (np.size(ind)>0): crit=-7
  #print i,mol,ind,np.size(ind)
  if (np.max(yy[iii])>crit):
    plt.plot(Tg,yy,c=colo[count],ls=styl[count],lw=widt[count],label=mol)
    count = count + 1
plt.title('important molecules',fontsize=20)
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=20)
plt.xlim(Tmin,Tmax)
plt.ylim(-6.2,0.2)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
if (Tmax/Tmin>10):
  plt.xscale('log')
else:  
  minorLocator = MultipleLocator(sep)
  ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(0.2)
ax.yaxis.set_minor_locator(minorLocator)
leg = plt.legend(loc='lower right',fontsize=11,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== where are the elements? ================
ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','F','P','NI','MN','CR','ZN','ZR','RB','CU','B','BR','V','SR','W','el']
allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li',' ',' ','Ni','Mn','Cr','Zn','Zr','Rb','Cu',' ','Br',' ','Sr',' ','+']
exlist = [' He ',' Cl CL Ca CA Cr CR Co Cu CU ',' ',' Ne NE Na NA Ni NI ',' ',' Si SI Sr SR ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' Fe FE ',' ',' ',' ',' ',' ',' ',' ',' ',' Br BR ',' ',' ',' ',' ',' ']
titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','fluorine','phosphorus','nickel','manganese','chromium','zinc','zirconium','rubidium','copper','boron','bromine','vanadium','strontium','tungston','charge carriers']
limits = [2,5,4,6,6,5,6,4,7,8,6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5]   
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
      plt.plot(Tg,yy,c=colo[count],ls=styl[count],lw=widt[count],label=molname)
      count = count + 1
  plt.title(titel,fontsize=20)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{\langle H\rangle}$',fontsize=20)
  plt.xlim(Tmin,Tmax)
  plt.ylim(nmin,nmax+1)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  if (Tmax/Tmin>10):
    plt.xscale('log')
  else:  
    minorLocator = MultipleLocator(sep)
    ax.xaxis.set_minor_locator(minorLocator)
  minorLocator = MultipleLocator(1.0)
  if (nmax-nmin>50): minorLocator = MultipleLocator(2.0)
  if (nmax-nmin>100): minorLocator = MultipleLocator(5.0)
  if (nmax-nmin>200): minorLocator = MultipleLocator(10.0)
  ax.yaxis.set_minor_locator(minorLocator)
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


