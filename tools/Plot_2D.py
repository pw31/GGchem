import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.tri as tri
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
pp=PdfPages('ggchem.pdf')

single_figures = 0    # 0=no, 1 for pdf, 2 for png

def plot_points():
  plt.scatter(xx,yy,marker='.',s=0.1,c='black')

def annotate_planets():
  plt.scatter(xpl,ypl,marker='o',s=60,c='blue',clip_on=False,zorder=5)
  ax.annotate(r'$\rm{Earth}$',xy=(xpl[0],ypl[0]), xycoords='data',zorder=5,
              xytext=( 15,-30), textcoords='offset points',
              arrowprops=dict(facecolor='black', shrink=0.15, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='left', verticalalignment='center')
  ax.annotate(r'$\rm{Mars, Venus}$',xy=(xpl[1],ypl[1]), xycoords='data',zorder=5,
              xytext=(-40,-20), textcoords='offset points',
              arrowprops=dict(facecolor='black', shrink=0.15, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='right', verticalalignment='center')
  #ax.annotate(r'$\rm{Mars}$',xy=(xpl[2],ypl[2]), xycoords='data',zorder=5,
  #            xytext=(-20,-50), textcoords='offset points',
  #            arrowprops=dict(facecolor='black', shrink=0.15, 
  #                            headwidth=7, headlength=6,width=2),
  #            horizontalalignment='center', verticalalignment='center')
  ax.annotate(r'$\rm{Jupiter}$',xy=(xpl[3],ypl[3]), xycoords='data',zorder=5,
              xytext=( 15, 30), textcoords='offset points',
              arrowprops=dict(facecolor='black', shrink=0.15, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='left', verticalalignment='center')
  ax.annotate(r'$\rm{Titan}$',xy=(xpl[4],ypl[4]), xycoords='data',zorder=5,
              xytext=(-15, 30), textcoords='offset points',
              arrowprops=dict(facecolor='black', shrink=0.15, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='right', verticalalignment='center')
  imo = 0
  for mo in Omodel:
    plt.scatter(xOl[imo],yOl[imo],marker='x',s=20,c=colo[imo], \
                label=mo,clip_on=False,zorder=100)
    imo = imo+1

file   = 'Static_Conc_2D.dat'
zusatz = ''
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
press = dat[:,2]                 # p [dyn/cm2]
ntot  = 0.0*nHtot
for i in range(3,4+NELEM+NMOLE):        # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat[:,i]
ind  = np.where(keyword=='N2')[0][0]    # subtract N2
ntot = ntot - 10**dat[:,ind]
ind  = np.where(keyword=='NH3')[0][0]   # subtract NH3
ntot = ntot - 10**dat[:,ind]
lntot = np.log10(ntot)

ind  = np.where(keyword=='epsH')[0][0]
print keyword[ind]
epsH = 10**dat[:,ind]
ind  = np.where(keyword=='epsC')[0][0]
print keyword[ind]
epsC = 10**dat[:,ind]
ind  = np.where(keyword=='epsO')[0][0]
print keyword[ind]
epsO = 10**dat[:,ind]
ind  = np.where(keyword=='epsN')[0][0]
print keyword[ind]
epsN = 10**dat[:,ind]
ind  = np.where(keyword=='H2O')[0][0]
print keyword[ind]
cH2O = 10**(dat[:,ind]-lntot)
ind  = np.where(keyword=='CO2')[0][0]
print keyword[ind]
cCO2 = 10**(dat[:,ind]-lntot)
ind  = np.where(keyword=='CH4')[0][0]
print keyword[ind]
cCH4 = 10**(dat[:,ind]-lntot)
ind  = np.where(keyword=='H2')[0][0]
print keyword[ind]
cH2  = 10**(dat[:,ind]-lntot)
ind  = np.where(keyword=='O2')[0][0]
print keyword[ind]
cO2  = 10**(dat[:,ind]-lntot)
ind  = np.where(keyword=='CO')[0][0]
print keyword[ind]
cCO  = 10**(dat[:,ind]-lntot)
ind  = np.where(keyword=='SC')[0][0]
print keyword[ind]
SatC = 10**dat[:,ind]
ind  = np.where(keyword=='SH2O')[0][0]
print keyword[ind]
SatH2O = 10**dat[:,ind]
ind  = np.where(keyword=='SH2O[l]')[0][0]
print keyword[ind]
SatH2O_l = 10**dat[:,ind]

xmin = 0.0
xmax = 0.34
ymin =-1.0
ymax = 1.0

#xx = (epsC-epsO)/(epsO+epsC)
#yy = epsH/(epsO+epsC+epsH)
#xlab='(C-O)/(C+O)'
#ylab='H/(H+O+C)'
#xx = (epsC-epsH)/(epsH+epsC)
#yy = epsO/(epsO+epsC+epsH)
#xlab='(C-H)/(C+H)'
#ylab='O/(H+O+C)'
yy = (epsO-epsH)/(epsO+epsH)
xx = epsC/(epsO+epsC+epsH)
#xm = (xmin+xmax)/2
#xx = xm + (epsC/(epsO+epsC+epsH)-xm)*(1.0-yy)/2.0
ylab='(O-H)/(O+H)'
xlab='C/(H+O+C)'

print "temperatures[K]:",np.min(Tg),np.max(Tg)
print " pressures[bar]:",np.min(press)/bar,np.max(press)/bar
print "           epsC:",np.min(epsC),np.max(epsC)
print "           epsO:",np.min(epsO),np.max(epsO)
print "ln(nH2O)",np.max(cH2O)
print "ln(nCO2)",np.max(cCO2)
print "ln(nCH4)",np.max(cCH4)
print "Sat(C[s])",np.max(SatC)
print "Sat(H2O[s])",np.max(SatH2O)
print "Sat(H2O[l])",np.max(SatH2O_l)

SatH2O = np.maximum(SatH2O,SatH2O_l)
SatMax = np.maximum(SatH2O,SatC)

#triang = tri.Triangulation(epsC, epsO)
#Smid = tri.LinearTriInterpolator(triang, SatMax)
#mask = np.where(Smid>1.0,1,0)
#triang.set_mask(mask)

print "S(H2O[s/l]):",np.min(SatH2O),np.max(SatH2O)
print "S(C[s])    :",np.min(SatC),np.max(SatC)
print "Smax       :",np.min(SatMax),np.max(SatMax)

#--- contour levels and tickmarks ---
zmax = 1.0    #0.0
zmin = 1.E-2  #-4.0
zticks1 = []
zticks2 = []
zz = zmax
ic = 0
dz = 0.002
imod = 100
while zz>zmin:
  zticks1.append(zz)
  if ic%imod == 0:
    zticks2.append(zz)
  ic = ic+1
  zz = zz-dz
zticks1.append(zmin)
zticks2.append(zmin)
zticks1 = np.array(zticks1[::-1])
zticks2 = np.array(zticks2[::-1])

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

#=================================
# Oliver model data (400 K, 1 bar)
#=================================
Omodel = ['CI', 'BSE', 'CC', 'PWD-ARE1', \
          'BSE15', 'MORB', 'PWD-CI', 'PWD-BSE', \
          'PWD-CI+4H2O', 'PWD-CI+5H2O']
epsHdat = [0.8080, 2.406e-08, 1.553e-06, 1.0, \
           0.2509, 4.950e-08, 0.8284, 0.1365, \
           1.0, 0.6759386811258722]
epsCdat = [0.06012, 0.1860, 70.41, 3.162e-05, \
           0.001199, 0.1299, 0.06170, 122.9, \
           0.04171, 0.05899]
epsOdat = [0.4769, 0.01500, 141.3, 2.772e-05, \
           0.4890, 0.3236, 0.02230, 246.3, \
           0.2040, 0.04026]
epsNdat = [0.01333, 43852.0, 6176.9, 0.001, \
           1.508e-05, 3476.2, 0.01300, 0.007729, \
           0.008792, 0.01243]
### BURCAT Fe3O4[s] and [l]
Omodel = ['CI','MORB','CC','CC(no_phyl)','BSE','PWD']
epsHdat = [8.0806E-01, 1.9450E-07, 1.4974E-06, 1.0000E+00, 1.6638E-07, 1.0000E+00]
epsCdat = [5.8576E-02, 4.0470E-08, 1.0616E-04, 1.0607E-02, 3.7317E-08, 3.1623E-05]
epsNdat = [1.0776E-02, 1.7209E-04, 9.5952E-03, 9.5952E-03, 1.0555E-03, 1.0000E-03]
epsOdat = [4.0204E-01, 1.3471E-08, 2.1307E-04, 5.2120E-01, 3.8676E-10, 2.7535E-05]

### Fe3O4[s] from SUPCRTBL, no [l] 
Omodel = ['CI','MORB','CC','CC(no_phyl)','BSE','PWD']
epsHdat = [8.0806E-01, 4.8862E-08, 1.4976E-06, 1.0000E+00, 2.5151E-08, 1.0000E+00]
epsCdat = [4.8323E-02, 6.3510E-09, 1.0616E-04, 1.0619E-02, 4.6703E-09, 3.1623E-05]
epsNdat = [1.0776E-02, 1.7209E-04, 9.5952E-03, 9.5952E-03, 1.0555E-03, 1.0000E-03]
epsOdat = [3.8496E-01, 1.5757E-08, 2.1307E-04, 5.2111E-01, 3.8888E-10, 2.7535E-05]

xOl = []
yOl = []
imo = 0
for mo in Omodel:
  Htot = epsHdat[imo]
  Ctot = epsCdat[imo]
  Otot = epsOdat[imo]
  Ntot = epsNdat[imo]
  yOl.append((Otot-Htot)/(Otot+Htot))
  xOl.append(Ctot/(Htot+Ctot+Otot))
  print '%12s: (O-H)/(O+H)=%16.8f , C/(H+C+O))=%16.8f' \
        % (mo,(Otot-Htot)/(Otot+Htot),Ctot/(Htot+Ctot+Otot))
  imo = imo+1
yOl = np.array(yOl)
xOl = np.array(xOl)

colo = ['blue','beige','coral','cadetblue','chartreuse','red','darkorchid','aqua','burlywood','chocolate','black','darkkhaki','pink','moccasin','cornflowerblue','darkgray']
#'darkgreen','gold','aquamarine','darkgoldenrod','darkorange','crimson','darkcyan','bisque','darkmagenta','darkolivegreen'

#====================================
# lines
#====================================
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


#=======================================================================================
mols = ['CH4','H2O','CO2','O2','H2','CO','NH3','N2','H2S','S2','S8','SO2','COS','H2SO4']
#=======================================================================================
for mol in mols:
  ind  = np.where(keyword==mol)[0][0]
  print keyword[ind]+zusatz
  cmol = 10**(dat[:,ind]-lntot)
  if (single_figures==1):  pp = PdfPages(mol+'.pdf')
  fig,ax = plt.subplots()
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])
  plt.xlabel(xlab, fontsize=17)
  plt.ylabel(ylab, fontsize=17)
  plt.subplots_adjust(left=0.09,right=1.03,top=1.0,bottom=0.10)
  if (single_figures==0):  plt.title(mol)
  map = plt.tricontourf(xx,yy,np.sqrt(cmol),zticks1,cmap=plt.cm.rainbow)

  #--- add colorbar ---
  cbr = plt.colorbar(map,pad=0.03,shrink=0.98)
  cbr.ax.get_yaxis().set_ticks([])
  for j, lab in enumerate([0.0001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]):
    zlin = np.sqrt(lab)
    zpos = (zlin-zmin)/(zmax-zmin)
    cbr.ax.text(1.2,zpos,lab,ha='left',va='center',fontsize=11)
    #print j,lab,zlin,zpos
  cbr.set_label(r'$n\,/\,(n_{\rm tot}-n_{N_2})$',fontsize=16)
  cbr.ax.get_yaxis().labelpad = 25

  #plt.tricontour(xx,yy,cH2-3*cCO2,[0.0], \
  #               norm=colors.Normalize(vmin=-1,vmax=1),cmap=cm.binary)
  #plt.tricontour(xx,yy,cO2-3*cCH4,[0.0], \
  #               norm=colors.Normalize(vmin=-1,vmax=1),cmap=cm.binary)
  #plt.tricontour(xx,yy,cCO-0.1*cH2O,[0.0], \
  #               norm=colors.Normalize(vmin=-1,vmax=1),cmap=cm.binary)
  plt.plot(xl1,yl1,c='grey',lw=1.5)
  plt.plot(xl2,yl2,c='grey',lw=1.5)
  plt.plot(xl3,yl3,c='grey',lw=1.5)
  #plt.tricontour(xx,yy,cCO2-cCH4,[0.0],linestyles=['dashed'],linewidths=[0.6])
  #plt.tricontour(xx,yy,cH2O-cCH4,[0.0],linestyles=['dashed'],linewidths=[0.6])
  #plt.tricontour(xx,yy,cH2O-cCO2,[0.0],linestyles=['dashed'],linewidths=[0.6])
  plt.plot(xl4,yl4,c='grey',ls='dashed',lw=1)
  plt.plot(xl5,yl5,c='grey',ls='dashdot',lw=1)
  plt.plot(xl6,yl6,c='grey',ls='dotted',lw=1)
  #plt.tricontour(xx,yy,cH2-3*cCO2,[0.0], \
  #               norm=colors.Normalize(vmin=-1,vmax=1),cmap=cm.binary)
  #plt.tricontour(xx,yy,cO2-3*cCH4,[0.0], \
  #               norm=colors.Normalize(vmin=-1,vmax=1),cmap=cm.binary)
  #plt.tricontour(xx,yy,cCO-0.1*cH2O,[0.0], \
  #               norm=colors.Normalize(vmin=-1,vmax=1),cmap=cm.binary)

  #plt.tricontour(xx,yy,SatC,[0.99], \
  #               norm=colors.Normalize(vmin=0,vmax=2),cmap=cm.Oranges)
  #plt.tricontour(xx,yy,SatH2O,[0.99], \
  #               norm=colors.Normalize(vmin=0,vmax=3),cmap=cm.coolwarm)
  #if (mol=='N2'): plot_points()
  #annotate_planets()
  #plt.title('%s , T = %d K , p = %.3f bar' % (mol,np.max(Tg),np.max(press)/bar))
  #leg = plt.legend(bbox_to_anchor=(0.15,1.0),loc='upper left',fontsize=7,fancybox=True,handlelength=1.5)
  #leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  if (single_figures<2):  plt.savefig(pp,format='pdf')
  if (single_figures==2): plt.savefig(mol+zusatz+'.png')
  if (single_figures==0): plt.clf()
  if (single_figures==1): pp.close()

if (single_figures==0): pp.close()

print '... written output to ggchem.pdf.'


