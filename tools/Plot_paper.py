import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5

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
delT  = (Tmax-Tmin)*0.30  #0.47
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

#================== the dust/gas mass ratio ====================
pp = PdfPages('dustgas.pdf')
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
pp.close()

#================ the gas phase element abundances ===================
pp = PdfPages('eps.pdf')
fig,ax = plt.subplots()
count = 0
ymax = -100.0
colors = []
for i in range(4+NELEM+NMOLE+2*NDUST,4+NELEM+NMOLE+2*NDUST+NELEM,1):
  elm = keyword[i]
  element = elm[3:]
  yy = dat[:,i]               # log10 eps
  ymax=np.max([ymax,np.max(yy)])            
  if (np.max(yy)>-20):
    plt.plot(Tg,yy,c=colo[count],ls=styl[count],lw=widt[count],label=element)
    colors.append(colo[count])
    count = count+1
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$\log\,\epsilon_{\rm gas}$',fontsize=20)
plt.xlim(100,2500)
plt.ylim(ymax-12.2,0.5)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
minorLocator = MultipleLocator(sep)
ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(1)
ax.yaxis.set_minor_locator(minorLocator)
sz = np.min([11,1+195.0/count])
leg = plt.legend(loc='lower right',fontsize=sz,fancybox=True)
for color,text in zip(colors,leg.get_texts()):
  text.set_color(color)
  text.set_size(11)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()
pp.close()

#================== solid particle densities ===================
T1 = [100,460,1000]
T2 = [460,1000,2250]
y1 = [-7.4,-7.6,-12.2]
y2 = [-3.4,-4.3,-4]
csize = [17,18,18]
filename = ['cond_lowT.pdf','cond_medT.pdf','cond_highT.pdf']
for iplot in range(0,3):
  Tmin = T1[iplot]
  Tmax = T2[iplot]
  ymin = y1[iplot]
  ymax = y2[iplot]
  iii = np.where((Tg>Tmin) & (Tg<Tmax))[0]
  solids = []
  smean = []
  for i in range(4+NELEM+NMOLE,4+NELEM+NMOLE+NDUST,1):
    solid = keyword[i]
    solids.append(solid[1:])
    smean.append(np.mean(dat[iii,i])) 
    ind = np.where(keyword == 'n'+solid[1:])[0]
    if (np.size(ind) == 0): continue
    ind = ind[0]
    yy = dat[:,ind]               # log10 nsolid/n<H>

  if (ymax>-99):
    pp = PdfPages(filename[iplot])
    print solids
    fig,ax = plt.subplots()
    indices = np.argsort(smean)
    colors = []
    count = 0
    for isolid in reversed(indices):
      solid = solids[isolid]
      ind = np.where(keyword == 'n'+solid)[0]
      if (np.size(ind) == 0): continue
      ind = ind[0]
      yy = dat[:,ind]               # log10 nsolid/n<H>
      if (np.max(yy[iii])>ymin):
        plt.plot(Tg[iii],yy[iii],c=colo[count],ls=styl[count],lw=widt[count],label=solid)
        colors.append(colo[count])
        count = count + 1
    #plt.title('condensates',fontsize=20)
    plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
    plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$',fontsize=20)
    plt.xlim(Tmin,Tmax)
    plt.ylim(ymin,ymax)
    plt.tick_params(axis='both', labelsize=14)
    plt.tick_params('both', length=6, width=1.5, which='major')
    plt.tick_params('both', length=3, width=1, which='minor')
    minorLocator = MultipleLocator(sep)
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator = MultipleLocator(1.0)
    ax.yaxis.set_minor_locator(minorLocator)
    leg = plt.legend(loc='lower right',fontsize=11,fancybox=True,
             handlelength=2.5,prop={'size':csize[iplot]},ncol=2)
    for color,text in zip(colors,leg.get_texts()):
      text.set_color(color)
      text.set_size(13)
    plt.tight_layout()
    plt.savefig(pp,format='pdf')
    plt.clf()
    pp.close()  

for iT in range(0,NPOINT):
  iact = 0
  outp = ' '
  for i in range(4+NELEM+NMOLE+NDUST,4+NELEM+NMOLE+2*NDUST,1):
    if (dat[iT,i]>-200): 
      iact=iact+1
      outp=outp+' '+keyword[i][1:]
  print Tg[iT],iact,outp  

print '... written output 3 .pdf files.'


