import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['axes.linewidth'] = 1.3
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.major.size'] = plt.rcParams['ytick.major.size'] = 7
plt.rcParams['xtick.minor.size'] = plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = plt.rcParams['ytick.major.width'] = 1.6
plt.rcParams['xtick.labelsize'] = plt.rcParams['ytick.labelsize'] = 15
pp = PdfPages('structure.pdf')
import sys

file  = 'Structure.dat'
dat   = np.loadtxt(file,skiprows=1)
zz    = dat[:,0]                 # z [km]
rho   = dat[:,1]                 # rho[g/cm3]
press = dat[:,2]                 # pressure [dyn/cm2]
Tg    = dat[:,3]                 # T [K]
nHtot = dat[:,4]                 # n<H> [cm-3]
mu    = dat[:,5]                 # mean molecular weight [amu]
gg    = dat[:,6]                 # gravity [m/s2]
Hp    = dat[:,7]                 # scale height [km]

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
zmin  = np.min(zz)
zmax  = np.max(zz)
Narg  = len(sys.argv)
if (Narg>1): zmin=float(sys.argv[1])
if (Narg>2): zmax=float(sys.argv[2])
iii   = np.where((zz>zmin) & (zz<zmax))[0]
pmin  = np.min(press[iii])/bar
pmax  = np.max(press[iii])/bar
pmin  = pmin*0.9
pmax  = pmax*1.1
Tmin  = np.min(Tg[iii])
Tmax  = np.max(Tg[iii])
Tmin  = Tmin*0.9
Tmax  = Tmax*1.1
ntot  = 0.0*nHtot
for i in range(3,4+NELEM+NMOLE): # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat[:,i]
logntot = np.log10(ntot)
colo = ['blue','black','silver','red','darkorange','gold','darkorchid','aqua','cadetblue','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod','darkkhaki','pink','moccasin']
#'darkolivegreen','darkmagenta','aquamarine','coral','burlywood',
#'beige','darkorange','crimson','darkcyan','bisque'
Ncolor = len(colo)
colo = colo*10
styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
widt = [2]*Ncolor*10
locmin = LogLocator(base=10.0,subs=(0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
locmaj = np.array([100,200,500,1000,2000,5000,10000])
locmaj = locmaj[np.where(locmaj<=Tmax)[0]]
locmaj = locmaj[np.where(locmaj>=Tmin)[0]]

#================== pressure structure ====================
fig,ax = plt.subplots(figsize=(5.5,7.5))
plt.plot(press/bar,zz,lw=4)
plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
plt.xlabel(r'$p\ \mathrm{[bar]}$',fontsize=20)
plt.xscale('log')
plt.xlim(pmin,pmax)
plt.ylim(zmin,zmax)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== temperature structure ====================
fig,ax = plt.subplots(figsize=(5.5,7.5))
plt.plot(Tg,zz,lw=4)
plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.xlim(Tmin,Tmax)
plt.ylim(zmin,zmax)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== mean molecular weight ====================
fig,ax = plt.subplots(figsize=(6.5,5.5))
plt.plot(mu,zz,lw=4)
plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
plt.xlabel(r'$\mu\ \mathrm{[amu]}$',fontsize=20)
plt.ylim(zmin,zmax)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== the dust/gas mass ratio ====================
fig,ax = plt.subplots(figsize=(6.5,5.5))
ind = np.where(keyword=='dust/gas')[0][0]
log10_dust_gas = dat[:,ind]
ymax = np.max(log10_dust_gas)
if (ymax>-10):
  plt.plot(10**log10_dust_gas,zz,lw=4)
  plt.xlabel(r'$\mathrm{dust/gas}$',fontsize=20)
  plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
  plt.ylim(zmin,zmax)
  #plt.xlim(10**(ymax-8),10**(ymax+1))
  plt.xlim(1.E-20,1.E-2)
  plt.xscale('log')
  #ax.yaxis.set_minor_locator(LogLocator(subs=[2,3,4,5,6,7,8,9]))
  #ax.xaxis.set_minor_locator(locmin)
  #ax.set_xticks(locmaj)
  #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

#================ the gas phase element abundances ===================
fig,ax = plt.subplots(figsize=(6.5,5.5))
count = 0
ymax = -100.0
for i in range(4+NELEM+NMOLE+2*NDUST,4+NELEM+NMOLE+2*NDUST+NELEM,1):
  yy = dat[:,i]               # log10 eps
  ymax=np.max([ymax,np.max(yy)])
ymin = ymax-10.0  
for i in range(4+NELEM+NMOLE+2*NDUST,4+NELEM+NMOLE+2*NDUST+NELEM,1):
  elm = keyword[i]
  element = elm[3:]
  yy = dat[:,i]               # log10 eps
  if (np.max(yy[iii])>ymin):
    plt.plot(yy,zz,c=colo[count],ls=styl[count],lw=widt[count],label=element)
    count = count+1
plt.title('gas phase element abundances',fontsize=16)
plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
plt.xlabel(r'$\log\,\epsilon_{\rm gas}$',fontsize=20)
plt.ylim(zmin,zmax)
plt.xlim(ymin,ymax+0.3)
#plt.yscale('log')
#ax.yaxis.set_minor_locator(locmin)
#ax.set_yticks(locmaj)
#ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
sz = np.min([11,1+195.0/count])
leg = plt.legend(loc='best',fontsize=sz,fancybox=True)
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
  yy = dat[:,ind]                # log10 nsolid/n<H>
  yy = yy + lognH - logntot      # log10 nsolid/ntot
  ymax = np.max([ymax,np.max(yy[iii])])
ymin = ymax-10
indices = np.argsort(smean)
if (ymax>-99):
  print solids
  fig,ax = plt.subplots(figsize=(5.5,7.5))
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
      plt.plot(yy[iii],zz[iii],c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  plt.title('condensates',fontsize=16)
  plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=18)
  plt.xlabel(r'$\mathrm{log}_{10}\ n_\mathrm{cond}/n_\mathrm{tot}$',fontsize=18)
  plt.ylim(zmin,zmax)
  plt.xlim(ymin,ymax+0.3)
  ax.xaxis.set_major_locator(MultipleLocator(2.0))
  ax.xaxis.set_minor_locator(MultipleLocator(1.0))
  #ax.xaxis.set_minor_locator(locmin)
  #ax.set_xticks(locmaj)
  #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  leg = plt.legend(loc='best',fontsize=10,fancybox=True,handlelength=2.5)
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

#================== log10 supersaturation ratios ===================
count = 0
for isolid in reversed(indices):
  solid = solids[isolid]
  ind = np.where(keyword == 'S'+solid)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  logS = dat[:,ind]              # log10 S
  if (np.max(logS[iii])>-6):
    count = count + 1
if (count>0):
  fig,ax = plt.subplots(figsize=(6.5,5.5))
  count = 0
  for isolid in reversed(indices):
    solid = solids[isolid]
    ind = np.where(keyword == 'S'+solid)[0]
    print solid,ind
    if (np.size(ind) == 0): continue
    if (count>=50): break
    ind = ind[0]
    logS = dat[:,ind]              # log10 S
    if (np.max(logS[iii])>-2):
      plt.plot(logS,zz,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count+1
  plt.title('supersaturation ratios',fontsize=16)
  plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
  plt.xlabel(r'$\mathrm{log}_{10}\ S$',fontsize=20)
  plt.ylim(zmin,zmax)
  plt.xlim(-7,0.5)
  #ax.xaxis.set_minor_locator(locmin)
  #ax.set_xticks(locmaj)
  #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  sz = np.min([10,1+160.0/count])
  col = 1
  if (count>30): 
    sz = np.min([10,1+160.0/count*2])
    col = 2
  leg = plt.legend(loc='best',fontsize=10,fancybox=True,
                   handlelength=3,prop={'size':sz},ncol=col)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

#================== supersaturation ratios ===================
count = 0
for isolid in reversed(indices):
  solid = solids[isolid]
  ind = np.where(keyword == 'S'+solid)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  #print solid,ind
  S = 10**dat[:,ind]              # S
  if (np.max(S[iii])>0.7):
    count = count + 1
if (count>0):
  fig,ax = plt.subplots(figsize=(6.5,5.5))
  count = 0
  for isolid in reversed(indices):
    solid = solids[isolid]
    ind = np.where(keyword == 'S'+solid)[0]
    if (np.size(ind) == 0): continue
    if (count>=50): break
    ind = ind[0]
    S = 10**dat[:,ind]              # S
    if (np.max(S[iii])>0.7):
      plt.plot(S,zz,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  plt.title('supersaturation ratios',fontsize=15)
  plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
  plt.xlabel(r'$S$',fontsize=20)
  plt.ylim(zmin,zmax)
  plt.xlim(0,1.05)
  #ax.xaxis.set_minor_locator(locmin)
  #ax.set_xticks(locmaj)
  #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  sz = np.min([10,1+160.0/count])
  col = 1
  if (count>30): 
    sz = np.min([10,1+160.0/count*2])
    col = 2
  leg = plt.legend(loc='best',fontsize=10,fancybox=True,
                   handlelength=3,prop={'size':sz},ncol=col)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

#================== some important molecules ====================
fig,ax = plt.subplots(figsize=(6.5,5.5))
mols  = ['H2','N2','H2O','O2','CO','COS','CO2','CH4','NH3','SO2','H2S','el']
mols  = np.array(mols)
count = 0
for i in range(3,4+NELEM+NMOLE): 
  mol = keyword[i]
  yy = dat[:,i]-logntot            # log10 nmol/ntot
  crit = -5.0
  ind = np.where(mols == mol)[0]
  if (np.size(ind)>0): crit=-7.0
  #print i,mol,ind,np.size(ind)
  if (np.max(yy[iii])>crit):
    plt.plot(yy,zz,c=colo[count],ls=styl[count],lw=widt[count],label=mol)
    count = count + 1
plt.title('important molecules',fontsize=15)
plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
plt.xlabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=20)
plt.ylim(zmin,zmax)
plt.xlim(-7.2,0.2)
#ax.xaxis.set_minor_locator(locmin)
#ax.set_xticks(locmaj)
#ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
leg = plt.legend(loc='best',fontsize=8,ncol=3,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== where are the elements? ================
ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','F','P','NI','MN','CR','ZN','ZR','RB','CU','B','BR','V','SR','W','el']
allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li',' ',' ','Ni','Mn','Cr','Zn','Zr','Rb','Cu',' ','Br',' ','Sr',' ','+']
exlist = [' He ',' Cl CL Ca CA Cr CR Co Cu CU ',' ',' Na NA Ni NI Ne NE ',' ',' Si SI Sr SR ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' Fe FE ',' ',' ',' ',' ',' ',' ',' ',' ',' Br BR ',' ',' ',' ',' ',' ']
titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','fluorine','phosphorus','nickel','manganese','chromium','zinc','zirconium','rubidium','copper','boron','bromine','vanadium','strontium','tungston','charge carriers']
limits = [5,5,5,5,6,5,6,4,7,8,6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,1.2]
condensates = indices
for i in range(0,30):
  fig,ax = plt.subplots()
  el = ellist[i]
  al = allist[i]
  ex = exlist[i]
  limit = limits[i]
  titel = titels[i]
  print titel+" ..."
  nmax = np.float(-200)
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
          yy = yy - logntot              # log10 nmol/ntot
          nmax = np.max([nmax,np.max(yy[iii])])
          maxy = maxy + 10**yy
          if (molname=='el'): nmin = np.min([nmin,np.min(yy[iii])])
          mollist.append(mol)   
          abulist.append(np.mean(yy))
  for isolid in condensates:
    solid = solids[isolid]
    isol = np.where(keyword == 'n'+solid)[0]
    if (np.size(isol) == 0): continue
    isol = isol[0]
    search = el
    if (len(el)==2): search=al
    ind = str.find(solid,search)
    found = 0
    while (ind>=0):
      #print solid,ind
      if (len(search)==2): found=1
      if (found==0):  
        if (ind==len(solid)-1): found=2
      if (found==0):  
        next1 = solid[ind+1]
        #print solid,search,next1
        if (next1.isupper() or next1.isdigit() or next1=='['): found=3
      if (found>0): break  
      ind = solid.find(search,ind+1)
      if (ind<0): break
      #print 'try again with rest ',ind,len(solid),solid[ind:]
    if (found>0):
      yy = dat[:,isol]               # log10 nsolid/n<H>
      yy = yy + lognH - logntot      # log10 nmol/ntot
      nmax = np.max([nmax,np.max(yy[iii])])
      maxy = maxy + 10**yy
      #print found,isol,keyword[isol],np.max(yy[iii])
      mollist.append(isol)   
      abulist.append(np.max(yy[iii]))
  if (nmax==-200): continue
  count = 0
  indices = np.argsort(abulist)
  maxy = np.log10(maxy)
  nmin = np.min([nmin,np.min(maxy[iii])-limit,nmax-6])
  for ind in reversed(indices):
    mol = mollist[ind]
    abu = abulist[ind]
    molname = keyword[mol]
    yy = dat[:,mol]                # log10 nmol [cm-3]
    if (mol<=4+NELEM+NMOLE):
      yy = yy - logntot            # log10 nmol/ntot
    else:
      yy = yy + lognH - logntot
      molname = molname[1:]
      if (str.find(molname,'[l]')<0):
        molname = molname+'[s]'
    #print mol,molname,abu,np.max(yy[iii])
    if (np.max(yy[iii]-maxy[iii])>-limit or molname=='el'):
      plt.plot(yy,zz,c=colo[count],ls=styl[count],lw=widt[count],label=molname)
      count = count + 1
  plt.title(titel,fontsize=18)
  plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=18)
  plt.xlabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=18)
  plt.ylim(zmin,zmax)
  plt.xlim(nmin,nmax+1)
  #ax.xaxis.set_minor_locator(locmin)
  #ax.set_xticks(locmaj)
  #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  sz = np.min([11,1+195.0/count])
  col = 1
  if (count>30): 
    sz = np.min([9,1+195.0/count*2])
    col = 2
  leg = plt.legend(loc='best',fontsize=8,fancybox=True,
             handlelength=3,prop={'size':sz},ncol=col)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

pp.close()
print '... written output to structure.pdf.'


