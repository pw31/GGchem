import matplotlib.pylab as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
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
cal  = 4.184       # 1 cal in J
mmHg = 1.3328E+3   # 1 mmHg in dyn/cm2


#Define functions used by each fit
def poly(T,A,B,C,D,E):
    val = A/T + B + C*T + D*T**2 + E*T**3
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

tmin = 50
tmax = 6000
temp = np.arange(tmin,tmax,1)

line1 = lines[142]
line2 = lines[147]
print line1
print line2
data1 = line1.split()
data2 = line2.split()
for i in range(3,8):
    data1[i] = float(data1[i])
    data2[i] = float(data2[i])
if (data1[2] == 'poly'):
    if (data1[1] == 'dg'):
        fit1 = poly(temp,*data1[3:8])
    else:
        fit1 = np.exp(poly(temp,*data1[3:8]))
elif (data1[2] == 'Woitke?'):
    fit1 = np.exp(newf(temp,*data1[3:6]))
elif (data1[2] == 'Yaws'):
    fit1 = yaws(temp,*data1[3:8])
elif (data1[2] == 'S&H'):
    fit1 = poly(temp,*data1[3:8])*cal
elif (data1[2] == 'Stock'):
    fit1 = stock(temp,*data1[3:8])
if (data2[2] == 'poly'):
    if (data2[1] == 'dg'):
        fit2 = poly(temp,*data2[3:8])
    else:
        fit2 = np.exp(poly(temp,*data2[3:8]))
elif (data2[2] == 'Woitke?'):
    fit2 = np.exp(newf(temp,*data2[3:6]))
elif (data2[2] == 'Yaws'):
    fit2 = yaws(temp,*data2[3:8])
elif (data2[2] == 'S&H'):
    fit2 = poly(temp,*data2[3:8])*cal
elif (data2[2] == 'Stock'):
    fit2 = stock(temp,*data2[3:8])
fig,ax = plt.subplots()
pmax = np.max([fit1,fit2])
pmin = np.min([fit1,fit2])
if (data1[1] == 'pvap'):
    plt.yscale('log')
    pmin = np.max([pmin,pmax*1.E-25])
minorLocator = MultipleLocator(100)
ax.xaxis.set_minor_locator(minorLocator)
plt.plot(temp,fit1,label = data1[0])
plt.plot(temp,fit2,label = data2[0])
specie = data1[0][:-3]
value = data1[1]
if (value == 'pvap'): unit = '[dyn/cm2]'
if (value == 'dg'): unit = '[J/mol]'
plt.title(specie)
plt.xlabel('T [K]')
plt.ylabel(value +' '+ unit)
plt.xlim(tmin,tmax)
plt.ylim(pmin,pmax)
#plt.xscale('log')
plt.legend(frameon=False)
plt.show()
#plt.savefig(pp,format='pdf')
#plt.clf()

fig,ax = plt.subplots()
minorLocator = MultipleLocator(100)
ax.xaxis.set_minor_locator(minorLocator)
plt.plot(temp,fit1/fit2)
plt.plot(temp,0*temp+1,c='black',ls='--')
plt.title(specie)
plt.xlabel('T [K]')
plt.ylabel(data1[0]+' / '+data2[0])
plt.xlim(tmin,tmax)
plt.ylim(0,2)
Nt = temp.size
ibest = 0
qbest = 1.e+99
for i in range(0,Nt):
  q = fit1[i]/fit2[i]
  if (abs(q-1)<qbest):
    ibest=i
    qbest=abs(q-1)
tmelt = temp[ibest]
print "melting temperature: ",tmelt,"K"
print "pvap(melting temp.): ",fit1[ibest],"bar"
plt.plot([tmelt,tmelt],[0,1],c='black',ls='--')
plt.show()
#plt.savefig(pp,format='pdf')
#plt.clf()


#for line1 in lines:
#    data1 = line1.split()
#    for line2 in lines:
#        data2 = line2.split()
#        if (data1[0:2] == data2[0:2]):
#            for i in range(3,8):
#                data2[i] = float(data2[i])
#            if (data2[2] == 'Yaws'):
#                fit = yaws(temp,*data2[3:8])
#                #plt.yscale('log')
#            elif (data2[2] == 'S&H'):
#                fit = poly(temp,*data2[3:8])*cal
#            elif (data2[2] == 'poly'):
#                if (data2[1] == 'dg'):
#                    fit = poly(temp,*data2[3:8])
#                else:
#                    fit = np.exp(poly(temp,*data2[3:8]))
#                    #plt.yscale('log')
#            elif (data2[2] == 'Woitke?'):
#                fit = np.exp(newf(temp,*data2[3:6]))
#                #plt.yscale('log')
#            elif (data2[2] == 'Stock'):
#                fit = stock(temp,*data2[3:8])
#            plt.plot(temp,fit,label = data2[2])
#    specie = data1[0]
#    value = data1[1]
#    if (value == 'pvap'): unit = '[dyn/cm2]'
#    if (value == 'dg'): unit = '[J/mol]'
#    plt.title(specie)
#    plt.xlabel('T [K]')
#    plt.ylabel(value +' '+ unit)
#    plt.legend(frameon=False)
#    #plt.show()
#    plt.savefig(pp,format='pdf')
#    plt.clf()
#pp.close()
#print '... written output to fit_comp.pdf.'
