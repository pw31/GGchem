import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
from scipy.optimize import curve_fit
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
pp = PdfPages('compare.pdf')

# nature constants and unit conversions 
Pa  = 10                # 1 Pa in dyn/cm2
bar = 1.E+6             # 1 bar in dyn/cm2
atm = 1.013E+6          # 1 atm in dyn/cm2
R   = 8.3144598         # gas constant in J/K/mol
cal = 4.184             # 1 cal in J   
bk  = 1.38065812E-16    # Boltzman constant erg/K

# set fit and plot ranges
Tref = 298.15     # reference temperature
Tmin = 50         # plot range
Tmax = 3500
N    = 200
TT   = np.zeros(N,dtype=float)
for i in range(0,N):
  TT[i] = np.exp(np.log(Tmin) + np.log(Tmax/Tmin)*i/(N-1))

#==========================================================================
# read converted GGchem data
#==========================================================================
print
print "reading BURCAT_generated.dat ..."
f = open('BURCAT_generated.dat','r')
lines = f.readlines()[:]
f.close
i = 0
GGnam = []
GGfit = []
for iline in range(0,99999):
  if (i>=len(lines)): break
  line = lines[i]
  i = i+1
  molnam = line.split()[0]
  #print molnam
  line = lines[i]
  i = i+1
  fit = np.array(line.split(),dtype='double')
  GGnam.append(molnam)
  GGfit.append(fit)
allnam = np.array(sorted(GGnam),dtype='str')
GGnam  = np.array(GGnam)

#==========================================================================
# read original BURCAT data
#==========================================================================
print
print "reading BURCAT.THR ..."
f = open('../BURCAT.THR','r')
lines = f.readlines()[:]
f.close
i = 0
BUnam = []
BUnamU= []
BUfit = []
width = 15
num_fields = 5
for iline in range(0,999999):
  if (i>=len(lines)): break
  line = lines[i]
  i+=1
  if (line[79:80]=='1'):
    if (line.find('G  ')==-1): continue
    molnam = line.split()[0]
    #print line
    fit = []
    line = lines[i]
    i+=1
    for j in range(num_fields): fit.append(float(line[width*j:width*(j+1)])) 
    line = lines[i]
    i+=1
    for j in range(num_fields): fit.append(float(line[width*j:width*(j+1)])) 
    line = lines[i]
    i+=1
    for j in range(num_fields): fit.append(float(line[width*j:width*(j+1)])) 
    fit = np.array(fit,dtype='double')
    #print molnam,fit[14]
    #--- check for same molecule (isomer)
    ind = np.where(np.array(BUnam)==molnam)[0]
    if (len(ind)>0):
      j = ind[0]
      dH298 = BUfit[j][14]
      print "there is another isomere of ",molnam
      print "previous: ",BUnam[j],dH298
      print "    this: ",molnam,fit[14]
      if (fit[14]<dH298):
        print "    the new isomer is more stable! replacing ..."
        BUfit[j] = fit
      else:
        print "    which is less stable => OK"
    else:
      BUnam.append(molnam)
      BUnamU.append(molnam.upper())
      BUfit.append(fit)
BUnamU = np.array(BUnamU)

iplot = 0
for mol in allnam:
  ind = np.where(GGnam==mol)[0]
  i = ind[0]
  mol1 = GGnam[i]
  fit1 = GGfit[i]
  print mol1
  ind = np.where(BUnamU==mol)[0]
  if (len(ind)>0):
    i = ind[0]
    mol2 = BUnam[i]
    fit2 = BUfit[i]
    print "-->",i,mol2
    dG_GG = 0.0*TT
    dG_BU = 0.0*TT
    a0 = fit1[0]
    a1 = fit1[1]
    a2 = fit1[2]
    a3 = fit1[3]
    a4 = fit1[4]
    a5 = fit1[5]
    a6 = fit1[6]
    b0 = fit1[7]
    b1 = fit1[8]
    b2 = fit1[9]
    b3 = fit1[10]
    b4 = fit1[11]
    b5 = fit1[12]
    b6 = fit1[13]
    i = 0
    for T in TT:
      if (T>1000.0):
        dGRT = a0*(1-np.log(T)) -a1*T/2 -a2*T**2/6 -a3*T**3/12 - a4*T**4/20 \
             + a5/T - a6  
      else:
        dGRT = b0*(1-np.log(T)) -b1*T/2 -b2*T**2/6 -b3*T**3/12 - b4*T**4/20 \
             + b5/T - b6  
      dG_GG[i] = dGRT*R*T
      i += 1
    a0 = fit2[0]
    a1 = fit2[1]
    a2 = fit2[2]
    a3 = fit2[3]
    a4 = fit2[4]
    a5 = fit2[5]
    a6 = fit2[6]
    b0 = fit2[7]
    b1 = fit2[8]
    b2 = fit2[9]
    b3 = fit2[10]
    b4 = fit2[11]
    b5 = fit2[12]
    b6 = fit2[13]
    i = 0
    for T in TT:
      if (T>1000.0):
        dGRT = a0*(1-np.log(T)) -a1*T/2 -a2*T**2/6 -a3*T**3/12 - a4*T**4/20 \
             + a5/T - a6  
      else:
        dGRT = b0*(1-np.log(T)) -b1*T/2 -b2*T**2/6 -b3*T**3/12 - b4*T**4/20 \
             + b5/T - b6  
      dG_BU[i] = dGRT*R*T
      i+=1
    iplot+=1  
    plt.figure(figsize=(7,8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.8,1])
    plt.subplot(gs[0])
    plt.plot(TT,dG_GG/1000,c='orange',lw=4.5,label='GGchem')
    plt.plot(TT,dG_BU/1000,c='dodgerblue',lw=2.0,ls='--',label='BURCAT')
    #plt.xscale('log') 
    #plt.xscale('log')
    plt.xlim(0,1.02*Tmax)
    #plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=16)
    plt.ylabel(r'$\Delta_\mathrm{f} G^0\ \mathrm{[kJ/mol]}$',fontsize=16)
    plt.title(mol2)
    plt.legend(loc='best')
    plt.subplot(gs[1])
    plt.plot(TT,(dG_GG-dG_BU)/1000,c='black',lw=2,label='GGchem-BURCAT')
    ymin = np.min(dG_GG-dG_BU)/1000
    ymax = np.max(dG_GG-dG_BU)/1000
    ymin = np.min([ymin-2,-5.0])
    ymax = np.max([ymax+2,+5.0])
    plt.xlim(0,1.02*Tmax)
    plt.ylim(ymin,ymax)
    plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=16)
    plt.ylabel(r'$\Delta_\mathrm{f} G^0\ \mathrm{[kJ/mol]}$',fontsize=16)
    plt.subplots_adjust(left=0.15,right=0.98,bottom=0.07,top=0.95,hspace=0.09)
    plt.legend(loc='best')
    plt.savefig(pp,format='pdf')
    plt.clf()

  #if (mol1=='HCO+'):
  #  pp.close()
  #  stop

pp.close()
print ' '
print ' '
print 'written output to compare.pdf'
