import matplotlib.pylab as plt
import numpy as np
import re
import os.path
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
bar  = 1.00E+6     # 1 bar in dyn/cm2

chemnam = 'Fe'     #choose specie

latexnm = re.sub(r'(\d+)',r'_\1',chemnam)
for phase in ['_cr','_l']:
  if os.path.isfile(chemnam+phase+'.txt'):
    f = open(chemnam+phase+'.txt','r')
    header = f.readline().split()
    f.close()
    molecule = header[-1]
    print molecule
    atoms = [] 
    stoich = []  
    s = re.findall(r'\D+',molecule)
    for i in s:
        atoms.append(i)  #atoms present
    s = re.findall(r'\d+',molecule)
    for i in s:
        stoich.append(float(i))  #stoich coeff
    print atoms
    print stoich
    data = np.genfromtxt(chemnam+phase+'.txt',skip_header=3,usecols=(0,6),dtype=float,invalid_raise=False)
    data = data[~np.any(np.isnan(data), axis=1)]
    tdat = data[:,0]
    gdat = data[:,1]
    if os.path.isfile(chemnam+'.txt'):
        gas = np.genfromtxt(chemnam+'.txt',skip_header=3,usecols=(0,6),dtype=float,invalid_raise=False)
        gas = gas[~np.any(np.isnan(gas), axis=1)]
        length = len(gas)
        i=0
        while not i==length:
            #print i,length
            if gas[i,0] not in tdat:
                #print data[i]
                gas = np.delete(gas,i,0)
                length = len(gas)
                continue
            i=i+1
        pdat = bar*np.exp(1000*(gdat-gas[:,1])/(Rgas*tdat))   #now in dyn/cm2
        if (phase=='_cr'): spdat = pdat; stdat = tdat
        if (phase=='_l'):  lpdat = pdat; ltdat = tdat
    j=0
    for atom in atoms[:-1]:
        data = np.genfromtxt(atom+'.txt',skip_header=3,usecols=(0,6),dtype=float,invalid_raise=False)
        data = data[~np.any(np.isnan(data), axis=1)]
        length = len(data)
        i=0
        while not i==length:
            #print i,length
            if data[i,0] not in tdat:
                #print data[i]
                data = np.delete(data,i,0)
                length = len(data)
                continue
            i=i+1
        gdat = gdat - stoich[j]*data[:,1]
        j=j+1
    #gdat = 1000*gdat   #now in J/mol
    if (phase=='_cr'): sgdat = gdat; stdat = tdat
    if (phase=='_l'):  lgdat = gdat; ltdat = tdat

#Define functions used by each fit
def poly(T,A,B,C,D,E):
    val = A/T + B + C*T + D*T**2 + E*T**3
    return val

def yaws(T,A,B,C,D,E):
    val = 10**(A + B/T + C*np.log10(T) + D*T + E*T**2)
    return val*mmHg

def newf(T,A,B,C):
    val = np.exp(A + B/(T + C))
    return val

def stock(T,A,B,C,D,E):
    val = - Rgas*T*(A/T + B*np.log(T) + C + D*T + E*T**2)
    return val

tmin = 100
tmax = 6000
temp = np.arange(tmin,tmax,1)

for line1 in lines:
    data1 = line1.split()
    if (data1[0] == chemnam+'[s]'):     #select solid
        for line2 in lines:
            data2 = line2.split()
            if (data2[0] == chemnam+'[l]'):    #select liquid
                if (data1[1:3]==data2[1:3]):
                    #line1 = lines[10]         #manual selection
                    #line2 = lines[13]
                    #data1 = line1.split()
                    #data2 = line2.split()
                    for i in range(3,8):
                        data1[i] = float(data1[i])
                        data2[i] = float(data2[i])
                    if (data1[2] == 'poly'):
                        if (data1[1] == 'dg'):
                            fit1 = poly(temp,*data1[3:8])
                        else:
                            fit1 = np.exp(poly(temp,*data1[3:8]))
                    elif (data1[2] == 'Woitke?'):
                        fit1 = newf(temp,*data1[3:6])
                    #elif (data1[2] == 'Yaws'):
                        #fit1 = yaws(temp,*data1[3:8])
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
                        fit2 = newf(temp,*data2[3:6])
                    #elif (data2[2] == 'Yaws'):
                        #fit2 = yaws(temp,*data2[3:8])
                    elif (data2[2] == 'S&H'):
                        fit2 = poly(temp,*data2[3:8])*cal
                    elif (data2[2] == 'Stock'):
                        fit2 = stock(temp,*data2[3:8])
                    fig,ax = plt.subplots()
                    if (data1[1] == 'dg'):
                        fit1 = fit1/1000 #kJ/mol
                    if (data2[1] == 'dg'):
                        fit2 = fit2/1000 #kJ/mol
                    pmax = np.max([fit1,fit2])
                    pmin = np.min([fit1,fit2])
                    minorLocator = MultipleLocator(100)
                    ax.xaxis.set_minor_locator(minorLocator)
                    plt.plot(temp,fit1,label = data1[0],color='orange')
                    plt.plot(temp,fit2,label = data2[0],color='deepskyblue')
                    specie = data1[0][:-3]
                    value = data1[1]
                    if (value == 'pvap'):
                        unit = '[dyn/cm2]'
                        plt.yscale('log')
                        pmin = np.max([pmin,pmax*1.E-25])
                        plt.scatter(ltdat,lpdat,color='dodgerblue',marker='o',label='data(l)')
                        plt.scatter(stdat,spdat,color='darkorange',marker='+',label='data(s)')
                        plt.ylabel(r'$ p^{\mathrm{vap}} \mathrm{[dyn/cm^2]}$',fontsize=16)
                    if (value == 'dg'):
                        unit = '[kJ/mol]'
                        plt.scatter(ltdat,lgdat,color='dodgerblue',marker='o',label='data(l)')
                        plt.scatter(stdat,sgdat,color='darkorange',marker='+',label='data(s)')
                        plt.ylabel(r'$\Delta_f G^\theta ({{{}}}) \quad \mathrm{{{}}}$'.format(latexnm,unit))
                    #plt.title(specie+' '+data1[2])
                    plt.title(specie)
                    plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=16)
                    plt.xlim(tmin,tmax)
                    plt.ylim(pmin,pmax)
                    #plt.xscale('log')
                    plt.legend(frameon=False,loc='upper left')
                    plt.subplots_adjust(left=0.14, right=0.96, top=0.94, bottom=0.14)
                    plt.show()
                    #plt.savefig(pp,format='pdf')
                    #plt.clf()

                    fig,ax = plt.subplots()
                    minorLocator = MultipleLocator(100)
                    ax.xaxis.set_minor_locator(minorLocator)
                    if (value == 'pvap'):
                      plt.plot(temp,fit1/fit2)
                      plt.plot(temp,0*temp+1,c='black',ls='--')
                      plt.ylabel(r'$ P_{{vap}} (\mathrm{{{}[s]}}) \div P_{{vap}} (\mathrm{{{}[l]}})$'.format(latexnm,latexnm))
                    if (value == 'dg'):
                      plt.plot(temp,fit1-fit2)
                      plt.plot(temp,0*temp,c='black',ls='--')
                      unit = '[kJ/mol]'
                      plt.ylabel(r'$\Delta_f G^\theta (\mathrm{{{}[s]}}) - \Delta_f G^\theta (\mathrm{{{}[l]}})  \quad \mathrm{{{}}}$'.format(latexnm,latexnm,unit))
                    #plt.title(specie+' '+data1[2])
                    plt.title(specie)
                    plt.xlabel(r'$T\ \mathrm{[K]}$')
                    plt.xlim(tmin,tmax)
                    Nt = temp.size
                    ibest = 0
                    qbest = 1.e+99
                    qmax  = 0
                    for i in range(0,Nt):
                      if (value == 'pvap'):
                        q = fit1[i]/fit2[i]
                        if (abs(q-1)<qbest):
                          ibest=i
                          qbest=abs(q-1)
                      if (value == 'dg'):
                        q = fit1[i]-fit2[i]
                        if (abs(q)<qbest):
                          ibest=i
                          qbest=abs(q)
                        if (abs(q)>qmax):
                          qmax = abs(q)
                    tmelt = temp[ibest]
                    print "Fit: ",data1[1:3]
                    print "melting temperature: ",tmelt,"K"
                    print "pvap(melting temp.): ",fit1[ibest],"bar"
                    plt.plot([tmelt,tmelt],[-qmax,qmax+1],c='black',ls='--')
                    plt.ylim(-qmax,qmax+2)
                    plt.show()
                    #plt.savefig(pp,format='pdf')
                    #plt.clf()


for line1 in lines:
    data1 = line1.split()
    if (data1[0][:-3] == chemnam):    #Choose condensate
        for line2 in lines:
            data2 = line2.split()
            if (data1[0:2] == data2[0:2]):
                for i in range(3,8):
                    data2[i] = float(data2[i])
                if (data2[2] == 'Yaws'):
                    fit = yaws(temp,*data2[3:8])
                    plt.yscale('log')
                    lab = 'Yaws(1999)'
                elif (data2[2] == 'S&H'):
                    fit = poly(temp,*data2[3:8])*cal
                    lab = 'S&H(1990)'
                elif (data2[2] == 'poly'):
                    if (data2[1] == 'dg'):
                        fit = poly(temp,*data2[3:8])
                        lab = 'Equation(2)'
                    else:
                        fit = np.exp(poly(temp,*data2[3:8]))
                        plt.yscale('log')
                        lab = 'Equation(6)'
                elif (data2[2] == 'Woitke?'):
                    fit = newf(temp,*data2[3:6])
                    plt.yscale('log')
                    lab = 'Equation(7)'
                elif (data2[2] == 'Stock'):
                    fit = stock(temp,*data2[3:8])
                    lab = 'Equation(3)'
                if (data2[1] == 'dg'):
                        fit = fit/1000 #kJ/mol
                plt.plot(temp,fit,label = lab)
        specie = data1[0]
        value = data1[1]
        if (value == 'pvap'):
            unit = '[dyn/cm2]'
            plt.ylabel(r'$ P_{{vap}} ({{{}}}) \quad \mathrm{{{}}}$'.format(latexnm,unit))
            if (data1[2] != 'Yaws'):
              if (data1[0][-2]=='s'): plt.scatter(stdat,spdat,marker='+',label='data',c='red')
              if (data1[0][-2]=='l'): plt.scatter(ltdat,lpdat,marker='+',label='data',c='red')
        if (value == 'dg'):
            unit = '[kJ/mol]'
            plt.ylabel(r'$\Delta_f G^\theta ({{{}}}) \quad \mathrm{{{}}}$'.format(latexnm,unit))
            if (data2[2] != 'S&H'):
              if (data1[0][-2]=='s'): plt.scatter(stdat,sgdat,marker='+',label='data',c='red')
              if (data1[0][-2]=='l'): plt.scatter(ltdat,lgdat,marker='+',label='data',c='red')
        plt.xlim(tmin,tmax)
        plt.title(specie)
        plt.xlabel(r'$T\ \mathrm{[K]}$')
        plt.legend(frameon=False)
        plt.show()
        #plt.savefig(pp,format='pdf')
        #plt.clf()
#pp.close()
#print '... written output to fit_comp.pdf.'
