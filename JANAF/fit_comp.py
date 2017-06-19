import matplotlib.pylab as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
pp = PdfPages('fit_comp.pdf')

# Open the file and read everything as a string into 'f'
f = open('fit_set.dat','r')
#Separate the lines
header = f.readline()
lines = f.readlines()[0:]
#Close the file
f.close()

Rgas = 8.1344598   # J/mol/K
cal  = 4.184       # 1 cal in erg
bar  = 1.E+6       # 1 bar in dyn/cm2
mmHg = 1.3328E+3   # 1 mmHg in dyn/cm2


#Define functions used by each fit
def poly(T,A,B,C,D,E):
    val = np.exp(A/T + B + C*T + D*T**2 + E*T**3)
    return val

def yaws(T,A,B,C,D,E):
    val = 10**(A + B/T + C*np.log10(T) + D*T + E*T**2)
    return val*mmHg

def newf(T,A,B,C):
    val = A + B/(T + C)
    return val

def stock(T,A,B,C,D,E):
    val = - Rgas*T*(A/T + B*np.log(T) + C + D*T + E*T**2)
    return val

tmin = 100
tmax = 6000
temp = np.arange(tmin,tmax,1)

line1 = lines[10]
line2 = lines[15]
print line1
print line2
data1 = line1.split()
data2 = line2.split()
for i in range(3,8):
    data1[i] = float(data1[i])
    data2[i] = float(data2[i])
fit1 = poly(temp,*data1[3:8])/bar       #pvap [bar]
fit2 = poly(temp,*data2[3:8])/bar       #pvap [bar]
plt.plot(temp,fit1,label = data1[0])
plt.plot(temp,fit2,label = data2[0])
specie = data1[0][:-3]
value = data1[1]
if (value == 'pvap'): unit = '[bar]'
if (value == 'dg'): unit = '[kJ/mol]'
plt.title(specie)
plt.xlabel('T [K]')
plt.ylabel(value +' '+ unit)
plt.xlim(tmin,tmax)
plt.xscale('log')
plt.yscale('log')
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(frameon=False)
#plt.show()
plt.savefig(pp,format='pdf')
plt.clf()

plt.plot(temp,fit1/fit2)
plt.plot(temp,0*temp+1,c='black',ls='--')
plt.title(specie)
plt.xlabel('T [K]')
plt.ylabel(data1[0]+' / '+data2[0])
plt.xlim(tmin,tmax)
plt.ylim(0,2)
plt.xscale('log')
plt.legend(frameon=False)
plt.savefig(pp,format='pdf')
pp.close()
print '... written output to fit_comp.pdf.'

#for line1 in lines:
#    data1 = line1.split()
#    for line2 in lines:
#        data2 = line2.split()
#        if (data1[0:2] == data2[0:2]):
#            for i in range(3,8):
#                data2[i] = float(data2[i])
#            if (data2[2] == 'Yaws'):
#                plt.plot(temp,yaws(temp,*data2[3:8])*1.3328E-03,label = data2[2])
#            elif (data2[2] == 'S&H'):
#                plt.plot(temp,poly(temp,*data2[3:8])/0.0041868,label = data2[2])
#            elif (data2[2] == 'poly'):
#                plt.plot(temp,poly(temp,*data2[3:8]),label = data2[2])
#            elif (data2[2] == 'Woitke?'):
#                plt.plot(temp,newf(temp,*data2[3:6]),label = data2[2])
#            elif (data2[2] == 'Stock'):
#                plt.plot(temp,stock(temp,*data2[3:8]),label = data2[2])
#    specie = data1[0]
#    value = data1[1]
#    if (value == 'pvap'): unit = '[bar]'
#    if (value == 'dg'): unit = '[kJ/mol]'
#    plt.title(specie)
#    plt.xlabel('T [K]')
#    plt.ylabel(value +' '+ unit)
#    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#    plt.legend(frameon=False)
#    plt.show()
