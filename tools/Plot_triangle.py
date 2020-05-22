import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.tri as tri
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
pp=PdfPages('ggchem.pdf')

def plot_points():
  plt.scatter(xx,yy,marker='.',s=0.1,c='black')

def annotate_planets():
  ax.annotate(r'$\rm{Earth}$',xy=(xpl[0],ypl[0]), xycoords='data',
              xytext=( 15,-30), textcoords='offset points',
              arrowprops=dict(facecolor='black', shrink=0.1, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='left', verticalalignment='center')
  ax.annotate(r'$\rm{Mars, Venus}$',xy=(xpl[1],ypl[1]), xycoords='data',
              xytext=(-40,-20), textcoords='offset points',
              arrowprops=dict(facecolor='black', shrink=0.1, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='right', verticalalignment='center')
  #ax.annotate(r'$\rm{Mars}$',xy=(xpl[2],ypl[2]), xycoords='data',
  #            xytext=(-20,-50), textcoords='offset points',
  #            arrowprops=dict(facecolor='black', shrink=0.1, 
  #                            headwidth=7, headlength=6,width=2),
  #            horizontalalignment='center', verticalalignment='center')
  ax.annotate(r'$\rm{Jupiter}$',xy=(xpl[3],ypl[3]), xycoords='data',
              xytext=( 15, 30), textcoords='offset points',
              arrowprops=dict(facecolor='black', shrink=0.1, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='left', verticalalignment='center')
  ax.annotate(r'$\rm{Titan}$',xy=(xpl[4],ypl[4]), xycoords='data',
              xytext=(-15, 30), textcoords='offset points',
              arrowprops=dict(facecolor='black', shrink=0.1, 
                              headwidth=7, headlength=6,width=2),
              horizontalalignment='right', verticalalignment='center')

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
press = dat[:,2]                 # p [dyn/cm2]
ntot  = 0.0*nHtot
for i in range(3,4+NELEM+NMOLE): # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat[:,i]
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
xmax = 0.4
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
zmax = 0.0
zmin = -4.0
zticks1 = [zmin]
zticks2 = [zmin]
zz = zmin
ic = 0
dz = 0.01
imod = 50
while zz<=zmax:
  zz = zz+dz
  ic = ic+1
  zticks1.append(zz)
  if ic%imod == 0:
    zticks2.append(zz)
zticks1 = np.array(zticks1)
zticks2 = np.array(zticks2)

colo = ['gold','cadetblue','coral','blue','beige','chartreuse','darkgreen','red','darkorchid','aqua','burlywood','chocolate','black','darkkhaki','pink','moccasin','cornflowerblue','darkgray']
#'aquamarine','darkgoldenrod','darkorange','crimson','darkcyan','bisque','darkmagenta','darkolivegreen'
Ncolor = len(colo)

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
  print '%8s: epsO=%16.8f , epsC=%16.8f , epsN=%16.8f , epsHe=%16.8f' \
        % (pl,eepsO[ipl],eepsC[ipl],eepsN[ipl],eepsHe[ipl])
  ipl = ipl+1
ypl = np.array(ypl)
xpl = np.array(xpl)
#xpl = xm + (xpl-xm)*(1.0-ypl)/2.0

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
###
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


#===================================================
mols = ['CH4','H2O','CO2','O2','H2','CO','NH3','N2']
#=================================================== 
for mol in mols:
  ind  = np.where(keyword==mol)[0][0]
  print keyword[ind]
  cmol = 10**(dat[:,ind]-lntot)
  fig,ax = plt.subplots()
  plt.xlim([xmin,xmax])
  plt.ylim([ymin,ymax])
  plt.xlabel(xlab, fontsize=18)
  plt.ylabel(ylab, fontsize=18)
  plt.subplots_adjust(left=0.09,right=1.03,top=1.0,bottom=0.10)
  #plt.tripcolor(triang,np.log10(cmol),zticks1,cmap=plt.cm.jet,shading='gouraud')
  plt.tricontourf(xx,yy,np.log10(cmol),zticks1,cmap=plt.cm.jet)

  #--- add colorbar ---
  cbr = plt.colorbar(ticks=zticks2)
  cbr.set_label(r'$\ln\ n_{\rm %s}\ /\ n_{\rm tot}$' % (mol),fontsize=16)

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
  plt.plot(xl4,yl4,c='grey',ls='--',lw=1)
  plt.plot(xl5,yl5,c='grey',ls='--',lw=1)
  plt.plot(xl6,yl6,c='grey',ls='--',lw=1)
  plt.tricontour(xx,yy,SatC,[0.99], \
                 norm=colors.Normalize(vmin=0,vmax=2),cmap=cm.Oranges)
  plt.tricontour(xx,yy,SatH2O,[0.99], \
                 norm=colors.Normalize(vmin=0,vmax=3),cmap=cm.Oranges)
  #plt.tricontour(xx,yy,cH2-3*cCO2,[0.0], \
  #               norm=colors.Normalize(vmin=-1,vmax=1),cmap=cm.binary)
  #plt.tricontour(xx,yy,cO2-3*cCH4,[0.0], \
  #               norm=colors.Normalize(vmin=-1,vmax=1),cmap=cm.binary)
  #plt.tricontour(xx,yy,cCO-0.1*cH2O,[0.0], \
  #               norm=colors.Normalize(vmin=-1,vmax=1),cmap=cm.binary)
  if (mol=='N2'): plot_points()
  annotate_planets()
  plt.title('%s , T = %d K , p = %.3f bar' % (mol,np.max(Tg),np.max(press)/bar))
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

ticks = np.log10([0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,1.E+300])
fig,ax = plt.subplots()
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xlabel(xlab, fontsize=18)
plt.ylabel(ylab, fontsize=18)
plt.subplots_adjust(left=0.09,right=1.03,top=0.9,bottom=0.10)
plt.tricontourf(xx,yy,np.log10(SatH2O),ticks, \
                norm=colors.Normalize(vmin=-2,vmax=2),cmap=plt.cm.jet)
cbr = plt.colorbar()
cbr.set_label(r'$\log\ S_{\rm H_2O}$',fontsize=16)
plt.tricontour(xx,yy,SatH2O,[0.99], \
                norm=colors.Normalize(vmin=0.99,vmax=2),cmap=cm.gray)
plt.tricontour(xx,yy,SatC,[0.99], \
               norm=colors.Normalize(vmin=0,vmax=2),cmap=cm.Oranges)
plot_points()
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

fig,ax = plt.subplots()
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xlabel(xlab,fontsize=18)
plt.ylabel(ylab,fontsize=18)
plt.subplots_adjust(left=0.09,right=1.03,top=0.98,bottom=0.10)
plt.tricontourf(xx,yy,np.log10(SatC),ticks, \
                norm=colors.Normalize(vmin=-2,vmax=2),cmap=plt.cm.jet)
cbr = plt.colorbar()
cbr.set_label(r'$\log\ S_{\rm C}$',fontsize=16)
plt.tricontour(xx,yy,SatC,[0.99], \
                norm=colors.Normalize(vmin=0.99,vmax=2),cmap=cm.gray)
plt.tricontour(xx,yy,SatH2O,[0.99], \
               norm=colors.Normalize(vmin=0,vmax=2),cmap=cm.Oranges)
plot_points()
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

ticks = np.log10([1.0,2.0,5.0,10.0,20.0,50.0,1.E+300])
fig,ax = plt.subplots()
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.xlabel(xlab, fontsize=18)
plt.ylabel(ylab, fontsize=18)
plt.subplots_adjust(left=0.09,right=1.03,top=0.98,bottom=0.10)
plt.tricontourf(xx,yy,np.log10(SatMax),ticks, \
                norm=colors.Normalize(vmin=-2,vmax=2),cmap=plt.cm.jet)
cbr = plt.colorbar()
cbr.set_label(r'$\log\ \max \{S_{\rm H_2O},S_{\rm C}\}$',fontsize=16)
plt.tricontour(xx,yy,SatMax,[0.99], \
               norm=colors.Normalize(vmin=0,vmax=2),cmap=cm.Oranges)
annotate_planets()
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

pp.close()
print '... written output to ggchem.pdf.'


