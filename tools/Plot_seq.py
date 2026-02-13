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
plt.rcParams['font.size'] = 14
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
xx = np.arange(1,NPOINT+1,1)

bar   = 1.E+6                    # 1 bar in dyn/cm2 
nHtot = dat[:,1]                 # n<H> [cm-3]
lognH = np.log10(nHtot)          
press = dat[:,2]                 # p [dyn/cm2]
xmin  = np.min(xx)-0.3
xmax  = np.max(xx)+0.3
Narg  = len(sys.argv)
if (Narg>1): xmin=float(sys.argv[1])
if (Narg>2): xmax=float(sys.argv[2])
iii   = np.where((xx>=xmin) & (xx<=xmax))[0]
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
locmaj = locmaj[np.where(locmaj<=xmax)[0]]
locmaj = locmaj[np.where(locmaj>=xmin)[0]]

#==============================================================
# planet data
#==============================================================
#https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
#Titan: Niemann et al (2005), Nature 438, 779
#1ppm = 1.E-6 = 1.E-4%
planet = ['Earth' ,'Venus', 'Mars','Jupiter','Titan']
N2     = [78.084  ,   3.5 ,   2.7 ,      0 ,   95.0 ]
O2     = [20.946  ,   0.0 ,  0.13 ,      0 ,      0 ]
CO2    = [0.0407  ,  96.5 , 95.32 ,      0 ,      0 ]
CH4    = [0.00018 ,     0 ,     0 , 3000E-4,    4.9 ]
H2O    = [0.80    , 20.E-4,210.E-4,   4.E-4,      0 ]
H2     = [0.000055,     0 ,     0 ,   89.8 ,   0.15 ]
CO     = [      0 , 17.E-4,  0.08 ,      0 ,      0 ]
NH3    = [      0 ,     0 ,     0 ,  260E-4,      0 ]
C2H6   = [      0 ,     0 ,     0 ,  5.8E-4,      0 ]
He     = [      0 ,     0 ,     0 ,   10.2 ,      0 ]
ipl = 0
eepsH = []
eepsC = []
eepsO = []
eepsN = []
eepsHe= []
xpl = []
ypl = []
for pl in planet:
  Htot = 2*H2[ipl] + 4*CH4[ipl] + 2*H2O[ipl] + 3*NH3[ipl] + 6*C2H6[ipl]
  Ctot = CO2[ipl] + CH4[ipl] + 2*C2H6[ipl] + CO[ipl]
  Otot = 2*O2[ipl] + 2*CO2[ipl] + H2O[ipl] + CO[ipl] 
  Ntot = 2*N2[ipl] + NH3[ipl]
  eepsH.append(Htot)
  eepsC.append(Ctot)
  eepsO.append(Otot)
  eepsN.append(Ntot)
  eepsHe.append(He[ipl])
  #xpl.append((Ctot-Otot)/(Ctot+Otot))
  #ypl.append(Htot/(Htot+Ctot+Otot))
  #xpl.append((Ctot-Htot)/(Ctot+Htot))
  #ypl.append(Otot/(Htot+Ctot+Otot))
  ypl.append((Otot-Htot)/(Otot+Htot))
  xpl.append(Ctot/(Htot+Ctot+Otot))
  print '%8s: epsH=%16.8f , epsO=%16.8f , epsC=%16.8f , epsN=%16.8f' \
        % (pl,Htot,Otot,Ctot,Ntot)
  ipl = ipl+1
ypl = np.array(ypl)
xpl = np.array(xpl)

def annotate_planets():
  plt.scatter(xpl,ypl,marker='o',s=60,c='blue',clip_on=False,zorder=5)
  ax.annotate(r'$\rm{Earth}$',xy=(xpl[0],ypl[0]), xycoords='data',zorder=5,
              xytext=( 30,-15), textcoords='offset points',size=12,
              arrowprops=dict(facecolor='blue', shrink=0.15, 
                              headwidth=7, headlength=6, width=2),
              horizontalalignment='left', verticalalignment='center')
  ax.annotate(r'$\rm{Mars, Venus}$',xy=(xpl[1],ypl[1]),xycoords='data',zorder=5,
              xytext=(-40,-20), textcoords='offset points',size=12,
              arrowprops=dict(facecolor='blue', shrink=0.15, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='right', verticalalignment='center')
  ax.annotate(r'$\rm{Jupiter}$',xy=(xpl[3],ypl[3]), xycoords='data',zorder=5,
              xytext=( 30, 15), textcoords='offset points',size=12,
              arrowprops=dict(facecolor='blue', shrink=0.15, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='left', verticalalignment='center')
  ax.annotate(r'$\rm{Titan}$',xy=(xpl[4],ypl[4]), xycoords='data',zorder=5,
              xytext=(-15, 40), textcoords='offset points',size=12,
              arrowprops=dict(facecolor='blue', shrink=0.15, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='right', verticalalignment='center')

#================== temperature-pressure structure ====================
fig,ax = plt.subplots()
plt.plot(xx,press/bar,lw=4)
plt.xlabel(r'$x$')
plt.ylabel(r'$p\ \mathrm{[bar]}$')
plt.xlim(xmin,xmax)
plt.ylim(pmin,pmax)
if (pmax>pmin*5): plt.yscale('log')
#ax.xaxis.set_minor_locator(locmin)
#ax.set_xticks(locmaj)
#ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== temperature-density structure ====================
fig,ax = plt.subplots()
plt.plot(xx,nHtot,lw=4)
plt.xlabel(r'$x$')
plt.ylabel(r'$n_\mathrm{\langle H\rangle}\ \mathrm{[cm^{-3}]}$')
plt.xlim(xmin,xmax)
plt.ylim(nHmin,nHmax)
if (nHmax>nHmin*5): plt.yscale('log')
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== the dust/gas mass ratio ====================
fig,ax = plt.subplots()
ind = np.where(keyword=='dust/gas')[0][0]
log10_dust_gas = dat[:,ind]
ymax = np.max(log10_dust_gas)
if (ymax>-10):
  plt.plot(xx,10**log10_dust_gas,lw=4)
  plt.xlabel(r'$x$')
  plt.ylabel(r'$\mathrm{dust/gas}$')
  plt.xlim(xmin,xmax)
  plt.ylim(10**(ymax-8),10**ymax*3)
  plt.yscale('log')
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
    plt.plot(xx,yy,c=colo[count],ls=styl[count],lw=widt[count],label=element)
    count = count+1
plt.title('gas phase element abundances')
plt.xlabel(r'$x$')
plt.ylabel(r'$\log\,\epsilon_{\rm gas}$')
plt.xlim(xmin,xmax)
plt.ylim(ymax-11.5,ymax+0.3)
sz = np.min([10,1+140.0/count])
leg = plt.legend(loc='lower right',fontsize=sz,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

fig,ax = plt.subplots(figsize=(6,5.5))
ind = np.where(keyword=='epsH')[0][0]
HH = 10**dat[:,ind]
ind = np.where(keyword=='epsC')[0][0]
CC = 10**dat[:,ind]
ind = np.where(keyword=='epsO')[0][0]
OO = 10**dat[:,ind]
ind = np.where(keyword=='epsN')[0][0]
NN = 10**dat[:,ind]
xxx = CC/(HH+CC+OO)
yyy = (OO-HH)/(OO+HH)
xxmin = 0.0
xxmax = 0.34
yymin =-1.0
yymax = 1.0
# the cut
plt.plot(xxx,yyy,c='red')
plt.text(xxx[0], yyy[0], "x=0",color='red',size=12)
plt.text(xxx[NPOINT-1], yyy[NPOINT-1], "x="+str(NPOINT),color='red',size=12)
plt.text(0.12, 0.55,"B",size=20)
plt.text(0.18,-0.25,"C",size=20)
plt.text(0.03,-0.70,"A",size=20)
# lines
epsC1=np.arange(0.0,0.25+1.E-8,0.25/100)
epsO1=0.5-2*epsC1
epsH1=0*epsC1+1.0
xl1 = epsC1/(epsH1+epsC1+epsO1)
yl1 = (epsO1-epsH1)/(epsO1+epsH1)
epsC2=np.arange(0.0,100.0+1.E-8,100.0/1000)
epsO2=0.5+2*epsC2
epsH2=0*epsC2+1.0
xl2 = epsC2/(epsH2+epsC2+epsO2)
yl2 = (epsO2-epsH2)/(epsO2+epsH2)
epsC3=np.arange(0.25,100.0+1.E-8,(100.0-0.25)/1000)
epsO3=2*(epsC3-0.25)
epsH3=0*epsC3+1.0
xl3 = epsC3/(epsH3+epsC3+epsO3)
yl3 = (epsO3-epsH3)/(epsO3+epsH3)
epsC4=np.arange(0.25,0.5+1.E-8,(0.5-0.25)/100)
epsO4=6*(epsC4-0.25)
epsH4=0*epsC4+1.0
xl4 = epsC4/(epsH4+epsC4+epsO4)
yl4 = (epsO4-epsH4)/(epsO4+epsH4)
epsC5=np.arange(0.0,0.5+1.E-8,(0.5)/100)
epsO5=0.5
epsH5=0*epsC5+1.0
xl5 = epsC5/(epsH5+epsC5+epsO5)
yl5 = (epsO5-epsH5)/(epsO5+epsH5)
epsC6=np.arange(1.0/6.0,100+1.E-8,(100.0-1.0/6.0)/1000)
epsO6=2*epsC6-1.0/6.0
epsH6=0*epsC6+1.0
xl6 = epsC6/(epsH6+epsC6+epsO6)
yl6 = (epsO6-epsH6)/(epsO6+epsH6)
plt.plot(xl1,yl1,c='grey',lw=1.5)
plt.plot(xl2,yl2,c='grey',lw=1.5)
plt.plot(xl3,yl3,c='grey',lw=1.5)
annotate_planets()
plt.xlabel(r'$\rm C/(H+O+C)$')
plt.ylabel(r'$\rm (O-H)/(O+H)$')
plt.xlim(xxmin,xxmax)
plt.ylim(yymin,yymax)
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
      plt.plot(xx[iii],yy[iii],c=colo[count],ls=styl[count],lw=widt[count],label=solid)
      count = count + 1
  plt.title('condensates')
  plt.xlabel(r'$x$')
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$')
  plt.xlim(xmin,xmax)
  plt.ylim(ymax-12,ymax+0.3)
  sz = np.min([11,1+140.0/count])
  leg = plt.legend(loc='best',fontsize=sz,fancybox=True,handlelength=2.5)
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
    print xx[iT],iact,outp  

#================== log10 supersaturation ratios ===================
fig,ax = plt.subplots()
count = 0
for isolid in reversed(indices):
  solid = solids[isolid]
  ind = np.where(keyword == 'S'+solid)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  logS = dat[:,ind]              # log10 S
  if (np.max(logS[iii])>-4):
    plt.plot(xx,logS,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
    count = count + 1
plt.plot([xmin,xmax],[0,0],ls="--",c='black',lw=1.0)
plt.title('supersaturation ratios')
plt.xlabel(r'$x$')
plt.ylabel(r'$\mathrm{log}_{10}\ S$')
plt.xlim(xmin,xmax)
plt.ylim(-4,+8)
sz = np.min([11,1+50.0/count])
col = 1
if (count>10): 
  sz = np.min([11,1+50.0/count*2])
  col = 2
leg = plt.legend(loc='best',fontsize=sz,fancybox=True,ncol=col,handlelength=2)
leg.get_frame().set_alpha(0.7)
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
  S = 10**dat[:,ind]
  if (np.max(S[iii])>0.7):
    plt.plot(xx,S,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
    count = count + 1
    if (count>100): break 
if (count>0):
  plt.title('supersaturation ratios')
  plt.xlabel(r'$x$')
  plt.ylabel(r'$S$')
  plt.xlim(xmin,xmax)
  plt.ylim(0,1.05)
  sz = np.min([11,1+50.0/count])
  col = 1
  if (count>10): 
    sz = np.min([11,1+50.0/count*2])
    col = 2
  leg = plt.legend(loc='best',fontsize=sz,fancybox=True,ncol=col,handlelength=2)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

#================== some important molecules ====================
fig,ax = plt.subplots()
mols  = ['CO2','N2','SO2','COS','H2S','S8','CO','H2O','H2','HCL','HF','Ar','Ne','H2SO4']
mols =  ['CO2','N2','SO2','H2S','S8','CO','H2O','H2','HCL','HF','H2SO4','FHO3S','Cl2']
mols  = np.array(mols)
ntot  = 0.0*nHtot
for i in range(3,4+NELEM+NMOLE): # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat[:,i]
lntot = np.log10(ntot)
count = 0
for i in range(3,4+NELEM+NMOLE): 
  mol = keyword[i]
  yy = dat[:,i]-lntot            # log10 nmol/ntot
  crit = -7
  ind = np.where(mols == mol)[0]
  if (np.size(ind)>0): crit=-9
  #print i,mol,ind,np.size(ind)
  if (np.max(yy[iii])>crit):
    plt.plot(xx,yy,c=colo[count],ls=styl[count],lw=widt[count],label=mol)
    count = count + 1
plt.title('important molecules')
plt.xlabel(r'$x$')
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$')
plt.xlim(0,xmax)
plt.ylim(-9.0,0.5)
leg = plt.legend(loc='upper right',fontsize=9,ncol=2,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== atomic pressures ====================
fig,ax = plt.subplots()
mols =  ['H','C','N','O','S','Cl','F']
mols  = np.array(mols)
ntot  = 0.0*nHtot
for i in range(3,4+NELEM+NMOLE): # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat[:,i]
lntot = np.log10(ntot)
count = 0
ymax = -999.0
for i in range(3,4+NELEM+NMOLE): 
  mol = keyword[i]
  yy = dat[:,i]-lntot            # log10 nmol/ntot
  crit = +1
  ind = np.where(mols == mol)[0]
  if (np.size(ind)>0): crit=-999
  #print i,mol,ind,np.size(ind)
  if (np.max(yy[iii])>crit):
    ymax = np.max([ymax,np.max(yy[iii])])
    plt.plot(xx,yy,c=colo[count],ls=styl[count],lw=widt[count],label=mol)
    count = count + 1
plt.title('free atoms')
plt.xlabel(r'$x$')
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{at}/n_\mathrm{tot}$')
ymax = ymax+5
ymin = ymax-130
plt.xlim(0,xmax)
#plt.ylim(ymin,ymax)
leg = plt.legend(loc='best',fontsize=9,ncol=2,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== linear changes ================
count = 0
fig,ax = plt.subplots()
maxchange = 0.0
thresh = 1.E-2
yshow  = 1.E-2
for i in range(3,4+NELEM+NMOLE): 
  mol = keyword[i]
  yy = 10.0**(dat[:,i]-lognH)            # nmol/n<H>
  #yy = 10.0**(dat[:,i]-lntot)           # nmol/ntot
  maxchange = np.max([maxchange,np.max(yy[iii])-np.min(yy[iii])])
for i in range(3,4+NELEM+NMOLE): 
  mol = keyword[i]
  yy = 10.0**(dat[:,i]-lognH)            # nmol/n<H>
  #yy = 10.0**(dat[:,i]-lntot)           # nmol/ntot
  change = np.max(yy[iii])-np.min(yy[iii])
  yy = yy - np.mean(yy[iii])
  if (change>maxchange*thresh):
    print mol,yy[np.max(iii)]
    plt.plot(xx[iii],yy[iii]*1.E+2,c=colo[count],lw=widt[count],label=mol)
    count = count + 1
for isolid in indices:
  solid = solids[isolid]
  ind = np.where(keyword == 'n'+solid)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  yy = 10.0**dat[:,ind]                 # nsolid/n<H>
  #yy = 10.0**(dat[:,ind]+lognH-lntot)   # nsolid/ntot
  change = np.max(yy[iii])-np.min(yy[iii])
  yy = yy - np.mean(yy[iii])
  if (change>maxchange*thresh):
    if (str.find(solid,'[l]')<0):
      solid = solid+'[s]'
    print solid,yy[np.max(iii)]
    plt.plot(xx[iii],yy[iii]*1.E+2,c=colo[count],ls='--',lw=widt[count],label=solid)
    count = count + 1
plt.xlabel(r'$x$')
plt.ylabel(r'$\frac{n-\langle n\rangle}{n_\mathrm{\langle H\rangle}}\,\mathrm{[\%]}$')
plt.xlim(xmin,xmax)
plt.ylim(-yshow,yshow)
leg = plt.legend(loc='best',fontsize=9,ncol=2,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== where are the elements? ================
ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','F','P','NI','MN','CR','ZN','ZR','RB','CU','B','BR','V','SR','W','el']
allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li',' ',' ','Ni','Mn','Cr','Zn','Zr','Rb','Cu',' ','Br',' ','Sr',' ','+']
exlist = [' He ',' Cl CL Ca CA Cr CR Co Cu CU ',' ',' Na NA Ni NI Ne NE ',' ',' Si SI Sr SR ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' Fe FE ',' ',' ',' ',' ',' ',' ',' ',' ',' Br BR ',' ',' ',' ',' ',' ']
titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','fluorine','phosphorus','nickel','manganese','chromium','zinc','zirconium','rubidium','copper','boron','bromine','vanadium','strontium','tungston','charge carriers']
limits = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
condensates = indices
for i in range(0,30):
  fig,ax = plt.subplots()
  el = ellist[i]
  al = allist[i]
  ex = exlist[i]
  limit = limits[i]
  titel = titels[i]
  print titel+" ..."
  nmax = 0.0
  nmin = 1.E+99
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
          #yy = yy - lognH               # log10 nmol/n<H>
          yy = yy - lntot                # log10 nmol/ntot
          nmax = np.max([nmax,np.max(yy[iii])])
          maxy = maxy + 10**yy
          if (molname=='el'): nmin = np.min([nmin,np.min(yy[iii])])
          mollist.append(mol)   
          abulist.append(np.mean(yy))
  nmaxgas = nmax
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
      yy = yy + lognH - lntot        # log10 nsolid/ntot
      nmax = np.max([nmax,np.max(yy[iii])])
      maxy = maxy + 10**yy
      #print found,isol,keyword[isol],np.max(yy[iii])
      mollist.append(isol)   
      abulist.append(np.max(yy[iii]))
  if (nmax==0.0): continue
  count = 0
  indices = np.argsort(abulist)
  maxy = np.log10(maxy)
  nmin = np.min([nmin,np.min(maxy[iii])-limit,nmax-7])
  for ind in reversed(indices):
    mol = mollist[ind]
    abu = abulist[ind]
    molname = keyword[mol]
    yy = dat[:,mol]                # log10 nmol [cm-3]
    if (mol<=4+NELEM+NMOLE):
      lim = limit+(nmax-nmaxgas)
      #yy = yy - lognH             # log10 nmol/n<H>
      yy = yy - lntot              # log10 nmol/ntot
    else:
      yy = yy + lognH - lntot      # log10 nsolid/ntot
      lim = limit
      molname = molname[1:]
      if (str.find(molname,'[l]')<0):
        molname = molname+'[s]'
    #print mol,molname,abu,np.max(yy[iii])
    if (np.max(yy[iii]-maxy[iii])>-lim or molname=='el'):
      yy = 10**yy
      plt.plot(xx,yy,c=colo[count],ls=styl[count],lw=widt[count],label=molname)
      count = count + 1
  plt.title(titel)
  plt.xlabel(r'$x$')
  #plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{\langle H\rangle}$')
  #plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$')
  plt.ylabel(r'$n_\mathrm{mol}/n_\mathrm{tot}$')
  plt.xlim(xmin,xmax)
  #plt.ylim(nmin,nmax+1)
  sz = np.min([11,1+140.0/count])
  col = 1
  if (count>30): 
    sz = np.min([9,1+195.0/count*2])
    col = 2
  leg = plt.legend(loc='best',fontsize=sz,fancybox=True,handlelength=3)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()


pp.close()
print '... written output to ggchem.pdf.'


