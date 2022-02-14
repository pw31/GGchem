import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True 
plt.rcParams['xtick.labelsize'] = plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['xtick.major.size'] = plt.rcParams['ytick.major.size'] = 7
plt.rcParams['xtick.minor.size'] = plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = plt.rcParams['ytick.major.width'] = 1.6
plt.rcParams['font.size'] = 16
pp = PdfPages('ggchem.pdf')
import sys

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

#================== temperature-pressure structure ====================
fig,ax = plt.subplots()
plt.plot(Tg,press/bar,lw=4)
plt.xlabel(r'$T\ \mathrm{[K]}$')
plt.ylabel(r'$p\ \mathrm{[bar]}$')
plt.xlim(Tmin,Tmax)
plt.ylim(pmin,pmax)
ax.set_xscale('log')
if (pmax>pmin*5): ax.set_yscale('log')
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.set_xticks(locmaj)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== temperature-density structure ====================
fig,ax = plt.subplots()
plt.plot(Tg,nHtot,lw=4)
plt.xlabel(r'$T\ \mathrm{[K]}$')
plt.ylabel(r'$n_\mathrm{\langle H\rangle}\ \mathrm{[cm^{-3}]}$')
plt.xlim(Tmin,Tmax)
plt.ylim(nHmin,nHmax)
plt.xscale('log')
if (nHmax>nHmin*5): plt.yscale('log')
ax.set_xticks(locmaj)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
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
  plt.xlabel(r'$T\ \mathrm{[K]}$')
  plt.ylabel(r'$\mathrm{cond/gas}$')
  plt.xlim(Tmin,Tmax)
  plt.ylim(10**(ymax-8),10**ymax*3)
  plt.xscale('log')
  plt.yscale('log')
  ax.set_xticks(locmaj)
  ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
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
plt.title('gas phase element abundances')
plt.xlabel(r'$T\ \mathrm{[K]}$')
plt.ylabel(r'$\log\,\epsilon_{\rm gas}$')
plt.xlim(Tmin,Tmax)
plt.ylim(ymax-11.5,ymax+0.3)
plt.xscale('log')
ax.set_xticks(locmaj)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
sz = np.min([11,1+140.0/count])
leg = plt.legend(loc='best',fontsize=sz,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================ the gas phase element abundances ===================
fig,ax = plt.subplots()
count = 0
ymax = -100.0
sumeps = 0.0*Tg
for i in range(4+NELEM+NMOLE+2*NDUST,4+NELEM+NMOLE+2*NDUST+NELEM,1):
  yy = dat[:,i]               # log10 eps
  sumeps = sumeps + 10**yy
for i in range(4+NELEM+NMOLE+2*NDUST,4+NELEM+NMOLE+2*NDUST+NELEM,1):
  elm = keyword[i]
  element = elm[3:]
  yy = dat[:,i]               # log10 eps
  yy = np.log10(10**yy/sumeps)
  print elm,np.max(yy)
  ymax=np.max([ymax,np.max(yy)])            
  if (np.max(yy)>-20):
    plt.plot(Tg,yy,c=colo[count],ls=styl[count],lw=widt[count],label=element)
    count = count+1
plt.title('relative gas phase element abundances')
plt.xlabel(r'$T\ \mathrm{[K]}$')
plt.ylabel(r'$\log\,\epsilon_{\rm gas}/\epsilon_{\rm gas}^{\rm tot}$')
plt.xlim(Tmin,Tmax)
plt.ylim(ymax-6.0,ymax+0.3)
plt.xscale('log')
ax.set_xticks(locmaj)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
sz = np.min([11,1+140.0/count])
leg = plt.legend(loc='best',fontsize=sz,fancybox=True)
leg.get_frame().set_alpha(0.7)
count = 0
for text in leg.get_texts():
  text.set_color(colo[count])
  count = count+1
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
  ymin = -15
indices = np.argsort(smean)
if (ymax>-99):
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
    if (np.max(yy[iii])>-99): print solid,ind,np.max(yy[iii])
    if (np.max(yy[iii])>ymin):
      plt.plot(Tg[iii],yy[iii],c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  plt.title('condensates')
  plt.xlabel(r'$T\ \mathrm{[K]}$')
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$')
  plt.xscale('log')
  plt.xlim(Tmin,Tmax)
  plt.ylim(ymax-8,ymax+0.3)
  ax.set_xticks(locmaj)
  ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  col = 3
  sz = np.min([11,1+70.0/(count/col)])
  leg = plt.legend(loc='best',fontsize=11,fancybox=True,
                   handlelength=2.5,prop={'size':sz},ncol=3)
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

#================== condensation groups =============================
  family = ['silicates','feldspar','Ca-Al-Ti','Fe-Ni-Cr','graphite','iron-oxides','metal-oxides',
            'sulphide','phyllosilicates','carbonates','halide','P-compounds','ices','other']
  Nfam   = len(family)
  member = []
  member.append(['SiO','SiO2','SiO2[l]','MgSiO3','MgSiO3[l]','Mg2SiO4','Mg2SiO4[l]',
                 'Mn2SiO4','NaAlSiO4','Fe2SiO4','Na2SiO3','CaMgSi2O6','KAlSiO4','CaSiO3',
                 'Ca3Fe2Si3O12','Mn3Al2Si3O12','Na2SiO3[l]','Ca2SiO4','KAlSi2O6',
                 'NaAlSi2O6','MnSiO3','NaFeSi2O6','NaCrSi2O6'])
  member.append(['KAlSi3O8','CaAl2Si2O8','NaAlSi3O8'])
  member.append(['Al2O3','Al2O3[l]','Ca3Al2Si3O12','MnTiO3','CaTiSiO5','Ca2Al2SiO7','Ti4O7',
                 'FeAl2O4','Ca2MgSi2O7','Ca2MgSi2O7','TiO2','CaTiO3','MgAl2O4','FeTiO3',
                 'Mg2TiO4[l]','MgAl2O4[l]','MgTi2O5','MgTi2O5[l]','Fe2TiO4','CaO[l]',
                 'Ti3O5'])
  member.append(['Fe[l]','Fe','Ni[l]','Ni','Cr'])
  member.append(['C'])
  member.append(['Fe2O3','Fe3O4','FeO','FeO[l]'])
  member.append(['Mn2O3','MgCr2O4','Cr2O3','MgFe2O4','MnO','VO','V2O3'])
  member.append(['FeS','FeS[l]','FeS2','MnS','CaSO4','Ni3S2'])
  member.append(['Mg3Si2O9H4','NaMg3AlSi3O12H2','MnAl2SiO7H2','Mg3Si4O12H2','FeAl2SiO7H2',
                 'KFe3AlSi3O12H2','CaAl2Si2O10H4','Fe3Si2O9H4','KMg3AlSi3O12H2','CaAl2Si4O16H8',
                 'Al2Si2O9H4','Ca2FeAl2Si3O13H','Ca2MnAl2Si3O13H','NaAl3Si3O12H2'])
  member.append(['CaMgC2O6','CaCO3'])
  member.append(['NaCl','KCl','MgF2','CaF2'])
  member.append(['Ca5P3O13H','Ca5P3O12F','Mg3P2O8'])
  member.append(['H2O','H2O[l]','NH3'])
  member.append([' '])
  file  = 'data/DustChem.dat'
  data  = open(file)
  lines = data.readlines()
  data.close
  Nline = len(lines)
  elements = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl',
              'Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',
              'Br','Kr','Rb','Sr','Y ','Zr','W ']
  emass    = [1.008 ,4.0026,6.94  ,9.0122,10.81 ,12.011,14.007,15.999,18.998,20.180,22.990,24.305,
              26.982,28.085,30.974,32.06 ,35.45 ,39.948,39.098,40.078,44.956,47.867,50.942,51.996,
              54.938,55.845,58.933,58.693,63.546,65.38 ,69.723,72.63 ,74.922,78.96 ,79.904,83.798,
              85.468,87.62 ,88.906,91.224,183.84]
  Nelem = len(elements)
  dust_gas = 10**log10_dust_gas
  rho_fam = np.zeros([Nfam,len(dust_gas)],dtype='float')
  rho_all = np.zeros(len(dust_gas),dtype='float')
  for isolid in indices:
    solid = solids[isolid]
    ind = np.where(keyword == 'n'+solid)[0]
    if (np.size(ind) == 0): continue
    ind = ind[0]
    yy = 10**dat[:,ind]               # nsolid/n<H>
    if (np.max(yy)<1.E-20): continue
    found = False
    ifam = Nfam-1
    for fam in range(0,Nfam):
      for mem in member[fam]:
        if (solid == mem):
          found = True
          ifam = fam
          break
      if (found): break
    rhom = 0.0  
    mass = 0.0
    for l in range(0,Nline):
      if (len(lines[l])<6): continue  
      test = lines[l].split()[0]
      search = solid
      if (str.find(search,'[l]')<0): search = search+'[s]'
      if (search==test.strip()):
        #print lines[l].strip()
        #print lines[l+1]
        name = lines[l].split()[1] 
        rhod = float(lines[l+1].split()[0])
        Nel = int(lines[l+2].split()[0])
        for el in range(0,Nel):
          stoich = int(lines[l+3+el].split()[0])
          elem   = lines[l+3+el].split()[1].strip()
          for i in range(0,Nelem):
            if (elem==elements[i].strip()):
              mass = mass + stoich*emass[i]
              break
          #print stoich,elem,emass[i]
        break
    if (mass==0.0):
      print "*** condensate not found: ",search      
      stop
    print "%16s %20s %5.3f %7.2f -> %16s" %(search,name,rhod,mass,family[ifam])
    yfam = yy*mass                                          # rho_family/n<H>
    rho_fam[ifam] = rho_fam[ifam] + yfam
    rho_all       = rho_all       + yfam
  
  fig,ax = plt.subplots()
  #plt.plot(Tg,dust_gas,c='black',lw=4)
  count = 0
  for ifam in range(0,Nfam):
    mass_ratio = rho_fam[ifam]
    mass_ratio = mass_ratio[iii]/rho_all[iii]
    if (np.max(mass_ratio)>1.E-5):
      plt.plot(Tg[iii],dust_gas[iii]*mass_ratio,lw=2,c=colo[count],label=family[ifam])
      count = count + 1
  plt.xlabel(r'$T\ \mathrm{[K]}$')
  plt.ylabel(r'$\mathrm{cond/gas}$')
  ymax = np.max(np.log10(dust_gas))
  plt.xlim(Tmin,Tmax)
  plt.ylim(10**(ymax-9),10**ymax*4)
  plt.xscale('log')
  plt.yscale('log')
  #ax.yaxis.set_minor_locator(LogLocator(subs=[2,3,4,5,6,7,8,9]))
  ax.set_xticks(locmaj)
  ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  #ax.set_xticks(locmaj)
  #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  leg = plt.legend(loc='upper center',fontsize=9,fancybox=True,ncol=4)
  leg.get_frame().set_alpha(0.7)
  count = 0
  for text in leg.get_texts():
    text.set_color(colo[count])
    count = count+1
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()
    
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
  fig,ax = plt.subplots()
  count = 0
  for isolid in reversed(indices):
    solid = solids[isolid]
    ind = np.where(keyword == 'S'+solid)[0]
    #print solid,ind
    if (np.size(ind) == 0): continue
    if (count>=50): break
    ind = ind[0]
    logS = dat[:,ind]              # log10 S
    if (np.max(logS[iii])>-6):
      plt.plot(Tg,logS,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count+1
  plt.title('supersaturation ratios')
  plt.xlabel(r'$T\ \mathrm{[K]}$')
  plt.ylabel(r'$\mathrm{log}_{10}\ S$')
  plt.xscale('log')
  plt.xlim(Tmin,Tmax)
  plt.ylim(-7,0.5)
  ax.set_xticks(locmaj)
  ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  sz = np.min([11,1+140.0/count])
  col = 1
  if (count>20): 
    sz = np.min([11,1+200.0/count])
    col = 2
  if (count>40): 
    sz = np.min([11,1+250.0/count])
    col = 3
  leg = plt.legend(loc='best',fontsize=sz,fancybox=True,ncol=col)
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
  fig,ax = plt.subplots()
  count = 0
  for isolid in reversed(indices):
    solid = solids[isolid]
    ind = np.where(keyword == 'S'+solid)[0]
    if (np.size(ind) == 0): continue
    if (count>=50): break
    ind = ind[0]
    S = 10**dat[:,ind]              # S
    if (np.max(S[iii])>0.7):
      plt.plot(Tg,S,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  plt.title('supersaturation ratios')
  plt.xlabel(r'$T\ \mathrm{[K]}$')
  plt.ylabel(r'$S$')
  plt.xscale('log')
  plt.xlim(Tmin,Tmax)
  plt.ylim(0,1.05)
  ax.set_xticks(locmaj)
  ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  sz = np.min([11,1+140.0/count])
  col = 1
  if (count>20): 
    sz = np.min([11,1+200.0/count])
    col = 2
  if (count>40): 
    sz = np.min([11,1+250.0/count])
    col = 3
  leg = plt.legend(loc='best',fontsize=sz,fancybox=True,ncol=col)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

#================== some important molecules ====================
fig,ax = plt.subplots()
mols  = ['H2','H','N2','H2O','O2','CO','CO2','CH4','NH3','C2H2','el']
mols  = np.array(mols)
ntot  = 0.0*nHtot
for i in range(3,4+NELEM+NMOLE): # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat[:,i]
lntot = np.log10(ntot)
count = 0
for i in range(3,4+NELEM+NMOLE): 
  mol = keyword[i]
  yy = dat[:,i]-lntot            # log10 nmol/ntot
  crit = -1.0
  ind = np.where(mols == mol)[0]
  if (np.size(ind)>0): crit=-1.0
  #print i,mol,ind,np.size(ind)
  if (np.max(yy[iii])>crit):
    plt.plot(Tg,yy,c=colo[count],ls=styl[count],lw=widt[count],label=mol)
    count = count + 1
plt.title('important molecules')
plt.xlabel(r'$T\ \mathrm{[K]}$')
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$')
plt.xlim(Tmin,Tmax)
plt.ylim(-6.2,0.2)
plt.xscale('log')
ax.set_xticks(locmaj)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
leg = plt.legend(loc='best',fontsize=11,ncol=2,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== where are the elements? ================
ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','F','P','NI','MN','CR','ZN','ZR','RB','CU','B','BR','V','SR','W','el']
allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li',' ',' ','Ni','Mn','Cr','Zn','Zr','Rb','Cu',' ','Br',' ','Sr',' ','+']
exlist = [' He ',' Cl CL Ca CA Cr CR Co Cu CU ',' ',' Na NA Ni NI Ne NE ',' ',' Si SI Sr SR ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' Fe FE ',' ',' ',' ',' ',' ',' ',' ',' ',' Br BR ',' ',' ',' ',' ',' ']
titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','fluorine','phosphorus','nickel','manganese','chromium','zinc','zirconium','rubidium','copper','boron','bromine','vanadium','strontium','tungston','charge carriers']
limits = [2,2,2,2,6,5,6,4,7,8,6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,1.2]
condensates = indices
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
      nmax = np.max([nmax,np.max(yy[iii])])
      maxy = maxy + 10**yy
      #print found,isol,keyword[isol],np.max(yy[iii])
      mollist.append(isol)   
      abulist.append(np.max(yy[iii]))
  if (nmax==-100): continue
  count = 0
  indices = np.argsort(abulist)
  maxy = np.log10(maxy)
  nmin = np.min([nmin,np.min(maxy[iii])-limit,nmax-12])
  for ind in reversed(indices):
    mol = mollist[ind]
    abu = abulist[ind]
    molname = keyword[mol]
    yy = dat[:,mol]                # log10 nmol [cm-3]
    if (mol<=4+NELEM+NMOLE):
      yy = yy - lognH              # log10 nmol/n<H>
    else:
      molname = molname[1:]
      if (str.find(molname,'[l]')<0):
        molname = molname+'[s]'
    #print mol,molname,abu,np.max(yy[iii])
    if (np.max(yy[iii]-maxy[iii])>-limit or molname=='el'):
      plt.plot(Tg,yy,c=colo[count],ls=styl[count],lw=widt[count],label=molname)
      count = count + 1
  plt.title(titel)
  plt.xlabel(r'$T\ \mathrm{[K]}$')
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{\langle H\rangle}$')
  plt.xlim(Tmin,Tmax)
  plt.ylim(nmin,nmax+1)
  plt.xscale('log')
  ax.set_xticks(locmaj)
  ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  sz = np.min([9,1+120.0/count])
  col = 1
  if (count>20): 
    sz = np.min([9,1+200.0/count])
    col = 2
  if (count>40): 
    sz = np.min([9,1+250.0/count])
    col = 3
  leg = plt.legend(loc='best',fontsize=sz,ncol=col,fancybox=True)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

pp.close()
print '... written output to ggchem.pdf.'


