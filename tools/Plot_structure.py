import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
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

bar   = 1.E+6                    # 1 bar in dyn/cm2
zmax  = np.max(zz)
colo = ['blue','black','silver','red','darkorange','gold','darkorchid','aqua','cadetblue','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod','darkkhaki','pink','moccasin']
#'darkolivegreen','darkmagenta','aquamarine','coral','burlywood',
#'beige','darkorange','crimson','darkcyan','bisque'
Ncolor = len(colo)
colo = colo*10
styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
widt = [2]*Ncolor*10

#================== pressure structure ====================
fig,ax = plt.subplots()
plt.plot(press/bar,zz,lw=4)
plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
plt.xlabel(r'$p\ \mathrm{[bar]}$',fontsize=20)
plt.xscale('log')
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== temperature structure ====================
fig,ax = plt.subplots()
plt.plot(Tg,zz,lw=4)
plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== temperature structure ====================
fig,ax = plt.subplots()
plt.plot(mu,zz,lw=4)
plt.ylabel(r'$z\ \mathrm{[km]}$',fontsize=20)
plt.xlabel(r'$\mu\ \mathrm{[amu]}$',fontsize=20)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

pp.close()
print '... written output to structure.pdf.'


