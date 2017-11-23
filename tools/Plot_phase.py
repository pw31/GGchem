import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5

single_figures = 0
if (single_figures==0): pp=PdfPages('ggchem.pdf')

file   = 'Static_Conc_2D.dat'
#file   = 'results/Static_Conc_2D_nocond.dat'
#file   = 'results/Static_Conc_2D_eqcond.dat'
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
keyword = np.array(header.split())
dat = np.loadtxt(file,skiprows=3)

bar   = 1.E+6                    # 1 bar in dyn/cm2 
Tg    = dat[:,0]                 # T [K]
nHtot = dat[:,1]                 # n<H> [cm-3]
lognH = np.log10(nHtot)          
press = dat[:,2]                 # p [dyn/cm2]
logp  = np.log10(press/bar)
Tmin  = np.min(Tg)
Tmax  = np.max(Tg)
pmin  = np.min(logp)
pmax  = np.max(logp)
#Tmin  = 100
#Tmax  = 460
iii   = np.where((Tg>=Tmin) & (Tg<=Tmax) & (logp>=pmin) & (logp<=pmax))[0]
sep   = 20
if (Tmax-Tmin>1500): sep=100
if (Tmax-Tmin>1000): sep=50
if (Tmax-Tmin<600): sep=20
if (Tmax-Tmin<400): sep=10

for sp in range(3,4+NELEM+NMOLE):
  dat[:,sp] = dat[:,sp]-lognH[:]    # log(nmol) -> log(nmol/n<H>) for molecules 
dat = np.array(dat)

colo = ['gold','cadetblue','coral','blue','beige','chartreuse','darkgreen','red','darkorchid','aqua','burlywood','chocolate','black','darkkhaki','pink','moccasin','cornflowerblue','darkgray']
#'aquamarine','darkgoldenrod','darkorange','crimson','darkcyan','bisque','darkmagenta','darkolivegreen'
Ncolor = len(colo)

#================== where are the elements? ================
ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','F','P','NI','MN','CR','ZN','ZR','RB','CU','B','BR','V','SR','W','el']
allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li',' ',' ','Ni','Mn','Cr','Zn','Zr','Rb','Cu',' ','Br',' ','Sr',' ','+']
exlist = [' He HE ',' Cl CL Ca CA Cr CR Co Cu CU ',' ',' Na NA Ni NI ',' ',' Si SI Sr SR ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' Fe FE ',' ',' ',' ',' ',' ',' ',' ',' ',' Br BR ',' ',' ',' ',' ',' ']
titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','fluorine','phosphorus','nickel','manganese','chromium','zinc','zirconium','rubidium','copper','boron','bromine','vanadium','strontium','tungsten','charge carriers']

xl = np.zeros([NPOINT],dtype=np.float)
xr = np.zeros([NPOINT],dtype=np.float)
yl = np.zeros([NPOINT],dtype=np.float)
yr = np.zeros([NPOINT],dtype=np.float)
for ix in range(0,NPOINT):
  xl[ix] = np.exp(np.log(Tmax)+np.log(Tmin/Tmax)*np.max([0,ix-0.5])/(NPOINT-1))
  xr[ix] = np.exp(np.log(Tmax)+np.log(Tmin/Tmax)*np.min([NPOINT-1,ix+0.5])/(NPOINT-1))
for iy in range(0,NPOINT):
  yl[iy] = pmax+(pmin-pmax)*np.max([0,iy-0.5])/(NPOINT-1)
  yr[iy] = pmax+(pmin-pmax)*np.min([NPOINT-1,iy+0.5])/(NPOINT-1)

for i in range(0,30):
  el = ellist[i]
  al = allist[i]
  ex = exlist[i]
  titel = titels[i]
  print titel+" ..."
  fig,ax = plt.subplots(figsize=(7,6))
  nmax = np.float(-100)
  splist = []
  dat2 = np.zeros([len(dat[:,0]),len(dat[0,:])])
  for sp in range(3,4+NELEM+NMOLE)+range(4+NELEM+NMOLE+NDUST,4+NELEM+NMOLE+2*NDUST):
    spname = keyword[sp]
    ind = str.find(spname,el)
    if (ind < 0): 
      ind = str.find(spname,al)
    if (ind < 0 and el=='el'): 
      ind = str.find(spname,'-')
    if (ind >= 0):
      next1 = spname[ind:ind+2]
      next2 = spname[ind-1:ind+1]
      #print keyword[sp],next1,str.find(ex,next1),len(next1)
      if (len(next1)==1 or str.find(ex,next1)==-1 or spname=='SIS'):
        if (next2!='MN' and next2!='ZN'):
          text1 = spname[ind+len(el):]
          text2 = ''
          for c in text1:
            if c.isdigit(): 
              text2=text2+c
            else:
              break
          stoich = 1.0
          if (len(text2)>0): stoich = float(text2)
          j1=str.find(spname,'(')
          j2=str.find(spname,')')
          if (j1>=0 and j1<ind and ind<j2):
            text2=spname[j2+1:j2+2]
            stoich=stoich*float(text2)
          print el,spname,stoich
          lstoich = np.log10(stoich)
          dat2[:,sp] = dat[:,sp] + lstoich
          ymax = np.max(dat2[:,sp])
          nmax = np.max([nmax,ymax])
          if (ymax>-100): splist.append(sp)   
  if (nmax==-100): continue
  splist = np.array(splist)
  print keyword[splist]
  implist = []
  Nc = 0
  ii = 0
  for iy in range(0,NPOINT):
    for ix in range(0,NPOINT):
      cmax  = np.max(dat2[ii,splist])
      spmax = np.where(dat2[ii,splist]==cmax)[0][0]
      #print keyword[splist[0:3]]
      #print dat[ii,splist[0:3]]
      #print dat2[ii,splist[0:3]]
      if (spmax not in implist): 
        lab = keyword[splist[spmax]]
        if ('[l]' in lab): 
          lab=lab[1:]
        else: 
          if (lab[0]=='n'): lab=lab[1:]+'[s]'
        implist.append(spmax)
        cc = implist.index(spmax)  
        rect = patches.Rectangle((xl[ix],yl[iy]),xr[ix]-xl[ix],yr[iy]-yl[iy],
               linewidth=0.02,edgecolor=colo[cc],facecolor=colo[cc],label=lab)
        ax.add_patch(rect)
      else:  
        cc = implist.index(spmax)  
        rect = patches.Rectangle((xl[ix],yl[iy]),xr[ix]-xl[ix],yr[iy]-yl[iy],
               linewidth=0.02,edgecolor=colo[cc],facecolor=colo[cc])
        ax.add_patch(rect)
      #print ix,iy,Tg[ii],logp[ii],keyword[splist[spmax]],cc
      ii = ii+1
  implist = np.array(keyword[splist[implist]])
  print implist

  plt.title(titel,fontsize=20)
  plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=22)
  plt.ylabel(r'$\mathrm{log}_{10}\ p\ \mathrm{[bar]}$',fontsize=22)
  plt.xlim(Tmin,Tmax)
  plt.ylim(pmax,pmin)
  if (Tmax/Tmin>10):
    plt.xscale('log')
  else:  
    minorLocator = MultipleLocator(sep)
    ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  plt.tick_params(axis='both', labelsize=18)
  plt.tick_params('both', length=11, width=2, which='major')
  plt.tick_params('both', length=8, width=1.5, which='minor')
  #minorLocator = MultipleLocator(1.0)
  #if (nmax-nmin>50): minorLocator = MultipleLocator(2.0)
  #if (nmax-nmin>100): minorLocator = MultipleLocator(5.0)
  #if (nmax-nmin>200): minorLocator = MultipleLocator(10.0)
  #ax.yaxis.set_minor_locator(minorLocator)
  leg = plt.legend(loc='lower left',fontsize=13,fancybox=True)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  if (single_figures==0): plt.savefig(pp,format='pdf')
  if (single_figures==1): fig.savefig('phase_'+titel+'.png')
  plt.clf()

if (single_figures==0): pp.close()
if (single_figures==0): print '... written output to ggchem.pdf.'


