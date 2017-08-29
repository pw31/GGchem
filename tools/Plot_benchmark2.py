import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 2
pp = PdfPages('ggchem.pdf')

single_figures = 0

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
bk    = bk=1.380662E-16
Tg    = dat[:,0]                 # T [K]
nHtot = dat[:,1]                 # n<H> [cm-3]
ntot  = 0.0*nHtot
for i in range(4,4+NELEM+NMOLE): # without el, but including ions and cations
  ntot = ntot + 10**dat[:,i]
  #print keyword[i]
lognn = np.log10(ntot)
press = dat[:,2]/bar             # p [dyn/cm2] -> [bar]
pmin  = np.min(press)
pmax  = np.max(press)
pmin  = pmin/1.2
pmax  = pmax*1.2
nHmin = np.min(nHtot)
nHmax = np.max(nHtot)
nHmin = nHmin*0.9
nHmax = nHmax*1.1
Tmin  = np.min(Tg)
Tmax  = np.max(Tg)
Tmin  = 100.0
Tmax  = 6000.0
TEAmin = 410.0
#if (Tmax>4*Tmin): Tmax=4*Tmin
#if (Tmin<Tmax/3): Tmin=Tmax/3
sep = 20
if (Tmax-Tmin>1500): sep=100
Tmin  = Tmin*0.85
Tmax  = Tmax*1.05

file   = 'TEAoutBench2/results/TEAoutBench2.tea'
data   = open(file)
dum    = data.readline()
dum    = data.readline()
dum    = data.readline()
dum    = data.readline()
dum    = data.readline()
dum    = data.readline()
dum    = data.readline()
header = data.readline()
lines  = data.readlines()
data.close()
sp_tea = np.array(header.split())
print "TEA has ",sp_tea
Npoint = len(lines)-1
Nsp    = len(sp_tea)
#dat2  = np.loadtxt(file,skiprows=8)
dat2   = np.zeros((Npoint,Nsp),dtype='float')
for i in range(0,Npoint):
  lval = lines[i].split()
  for isp in range(0,Nsp):
    dat2[i,isp] = np.float(lval[isp])
p_tea     = dat2[:,0]            # TEA's pressure [bar]
T_tea     = dat2[:,1]            # TEA's temperature
ntot_tea  = p_tea*bar/bk/T_tea   # TEA's ntot [cm-3]
logn_tea  = np.log10(ntot)          
nH_tea    = dat2[:,np.where(sp_tea == 'H_g')[0]]
nHe_tea   = dat2[:,np.where(sp_tea == 'He_ref')[0]]
nH2_tea   = dat2[:,np.where(sp_tea == 'H2_ref')[0]]
nH2O_tea  = dat2[:,np.where(sp_tea == 'H2O_g')[0]]
nO2_tea   = dat2[:,np.where(sp_tea == 'O2_ref')[0]]
nCO_tea   = dat2[:,np.where(sp_tea == 'CO_g')[0]]
nCO2_tea  = dat2[:,np.where(sp_tea == 'CO2_g')[0]]
nCH4_tea  = dat2[:,np.where(sp_tea == 'CH4_g')[0]]
nC2H2_tea = dat2[:,np.where(sp_tea == 'C2H2_g')[0]]
nC2H4_tea = dat2[:,np.where(sp_tea == 'C2H4_g')[0]]
nN2_tea   = dat2[:,np.where(sp_tea == 'N2_ref')[0]]
nNH3_tea  = dat2[:,np.where(sp_tea == 'NH3_g')[0]]
nOH_tea   = dat2[:,np.where(sp_tea == 'OH_g')[0]]

Tind1 = np.where((Tg<Tmax)    & (Tg>Tmin))[0]
Tind2 = np.where((T_tea<Tmax) & (T_tea>TEAmin))[0]

styl = ['-','-','-','-','-','-','-','--','--','--','--','--','--','--',':',':',':',':',':',':',':','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.','-.']
widt  = [ 2 , 2 , 2 , 2 , 2 , 2 , 2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  , 2  ]

#================== temperature-pressure structure ====================
fig,ax = plt.subplots()
plt.plot(Tg,press,lw=4)
plt.plot(T_tea,p_tea,c='lightgray',lw=1)
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$p\ \mathrm{[bar]}$',fontsize=20)
plt.xlim(Tmin,Tmax)
plt.ylim(pmin,pmax)
#plt.xscale('log')
#plt.yscale('log')
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
minorLocator = MultipleLocator(sep)
ax.xaxis.set_minor_locator(minorLocator)
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

#================== some important molecules ====================
fig,ax = plt.subplots()
mols  = ['H2','H','N2','H2O','O2','CO','CO2','CH4','NH3','He','SI(CH3)4']
nmax  = np.float(0)
for mol in range(4,4+NELEM+NMOLE):
  yy = dat[:,mol]                # log10 nmol [cm-3]
  yy = yy - lognn                # log10 nmol/ntot
  nmax = np.max([nmax,np.max(yy)])
count = 0
for mol in mols:
  ind = np.where(keyword == mol)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  print mol,ind
  yy = dat[:,ind]                # log10 nmol [cm-3]
  yy = yy - lognn                # log10 nmol/ntot
  #if (np.max(yy)>nmax-6):
  plt.plot(Tg,yy,ls=styl[count],lw=4,label=mol)
  count = count + 1
plt.plot(T_tea[Tind2],np.log10(nH_tea[Tind2])   ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nH2_tea[Tind2])  ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nHe_tea[Tind2])  ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nCO_tea[Tind2])  ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nH2O_tea[Tind2]) ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nO2_tea[Tind2])  ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nCO2_tea[Tind2]) ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nCH4_tea[Tind2]) ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nN2_tea[Tind2])  ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nNH3_tea[Tind2]) ,c='lightgray',lw=1)
plt.plot(T_tea[Tind2],np.log10(nC2H2_tea[Tind2]),c='lightgray',lw=1)  
plt.title('important molecules',fontsize=20)
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=20)
plt.xlim(Tmin,Tmax)
plt.ylim(nmax-8,nmax+0.5)
plt.xscale('log')
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
plt.legend(loc='lower left',fontsize=10,fancybox=True)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()
if (single_figures==1): pp.close()

#================== where are the elements? ================
ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','F','P','NI','MN','CR','ZN','ZR','RB','CU','B','BR','V','SR','W','el']
allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li',' ',' ','Ni','Mn','Cr','Zn','Zr','Rb','Cu',' ','Br',' ','Sr',' ','+']
exlist = [' He ',' Cl CL Ca CA Cr CR Co Cu CU ',' ',' Na NA Ni NI ',' ',' Si SI Sr SR ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' Fe FE ',' ',' ',' ',' ',' ',' ',' ',' ',' Br BR ',' ',' ',' ',' ',' ']
titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','fluorine','phosphorus','nickel','manganese','chromium','zinc','zirconium','rubidium','copper','boron','bromine','vanadium','strontium','tungsten','charge carriers']
#         H  C   O  N SiS Na ClCa Ti
limits = [2,4.5,3.5,4,7,5,10,3,10,5,10,6,6,12,6,3.5,5,10,10,10,10,10,10,10,10,10,10,10,10,5]    #Orich
#limits = [2,5,3.5,6,6,5,6,4,6,6,6,6,6,6,6.5,5]  #Crich
for i in range(0,30):
  el = ellist[i]
  al = allist[i]
  ex = exlist[i]
  limit = limits[i]
  titel = titels[i]
  print
  print titel+" ..."
  print 'ggchem ...'
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
      plotit = 0
      if (len(next)==1): plotit=1
      if (str.find(ex,next)==-1 or molname=='SIS'): plotit=1
      if (el=='N' and molname=='MNH'): plotit=0
      if (plotit==1):
        yy = dat[:,mol]                # log10 nmol [cm-3]
        yy = yy - lognn                # log10 nmol/ntot
        nmax = np.max([nmax,np.max(yy[Tind1])])
        maxy = maxy + 10**yy
        mollist.append(mol)   
        abulist.append(np.mean(yy[Tind1]))
  if (len(abulist)<=0): continue
  if (single_figures==1):  pp = PdfPages('benchmark_'+titel+'.pdf')
  if (single_figures==1):  fig,ax = plt.subplots(figsize=(7,6))
  indices = np.argsort(abulist)
  count = 0
  maxy = np.log10(maxy)
  nmin = np.min([nmin,np.min(maxy)-limit,nmax-14])
  if (el=='el'): nmin=-30
  for ind in reversed(indices):
    mol = mollist[ind]
    abu = abulist[ind]
    molname = keyword[mol]
    yy = dat[:,mol]                    # log10 nmol [cm-3]
    yy = yy - lognn                    # log10 nmol/ntot
    if (np.max(yy[Tind1]-maxy[Tind1])>-limit):
      print molname,np.max(yy[Tind1]-maxy[Tind1])
      plt.plot(Tg,yy,ls=styl[count],lw=4,label=molname)
      count = count + 1

  if (al<>' '): el=al
  print 'TEA ...'
  NTEA = len(sp_tea)
  maxy = 0.0*dat2[:,0]
  for mol in range(2,NTEA):
    molname = sp_tea[mol]
    ind = str.find(molname,el)
    if (ind >= 0):
      next = molname[ind:ind+2]
      if (len(next)==1 or str.find(ex,next)==-1 or molname=='SiS_g'):
        yy = np.log10(dat2[:,mol])     # log10 nmol/ntot 
        maxy = maxy + 10**yy
  maxy = np.log10(maxy)
  for mol in range(2,NTEA):
    molname = sp_tea[mol]
    ind = str.find(molname,el)
    if (ind >= 0):
      next = molname[ind:ind+2]
      if (len(next)==1 or str.find(ex,next)==-1 or molname=='SiS_g'):
        yy = np.log10(dat2[:,mol])     # log10 nmol/ntot
        plotit = 0
        if (np.max(yy[Tind2]-maxy[Tind2])>-limit): plotit=1
        if (molname.find('ZrF4_g')>=0): plotit=1
        if (molname.find('ZrCl4_g')>=0): plotit=1
        if (molname.find('TiCl3_g')>=0): plotit=1
        if ((el<>'O') and (molname.find('TiOCl2_g')>=0)): plotit=1
        #print el,molname,plotit
        if (plotit==1):
          print sp_tea[mol],np.max(yy[Tind2]-maxy[Tind2])
          plt.plot(T_tea[Tind2],yy[Tind2],c='lightgray',lw=1.5)

  plt.title(titel,fontsize=20)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=20)
  plt.xscale('log')
  plt.xlim(100,Tmax)
  plt.ylim(nmin,nmax+1)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=7, width=1.5, which='major')
  plt.tick_params('both', length=5, width=1, which='minor')
  #minorLocator = MultipleLocator(sep)
  #ax.xaxis.set_minor_locator(minorLocator)
  plt.legend(loc='lower left',fontsize=10,fancybox=True)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  if (single_figures==0): plt.clf()
  if (single_figures==1): pp.close()

if (single_figures==0): pp.close()
print '... written output to ggchem.pdf.'


