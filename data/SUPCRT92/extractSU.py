import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 2
pp = PdfPages('dG.pdf')

# nature constants and unit conversions 
Pa  = 10                # 1 Pa in dyn/cm2
bar = 1.E+6             # 1 bar in dyn/cm2
atm = 1.013E+6          # 1 atm in dyn/cm2
R   = 8.3144598         # gas constant in J/K/mol
cal = 4.184             # 1 cal in J   

# set fit and plot ranges
Tref  = 298.15     # reference temperature
Tmin1 = 100        # fit range
Tmax1 = 2500
Tmin2 = 50         # plot range
Tmax2 = 3500
N     = 100
T1 = np.zeros(N,dtype=float)
T2 = np.zeros(N,dtype=float)
for i in range(0,N):
  T1[i] = np.exp(np.log(Tmin1) + np.log(Tmax1/Tmin1)*i/(N-1))
  T2[i] = Tmin2 + (Tmax2-Tmin2)*i/(N-1)

#==========================================================================
# Fit coefficients for delta_f Gibbs free energy [J/mol] of atoms from NIST
#==========================================================================
atdat = np.loadtxt('NISTatom.dat',skiprows=0,usecols=[1,2,3,4,5])
atnam = np.loadtxt('NISTatom.dat',skiprows=0,usecols=[6],dtype='str')
atdat = np.array(atdat,dtype='float')
atnam = np.array(atnam,dtype='str')
#==========================================================================
def Stock(Temp,a):
  if (np.isscalar(Temp)):
    dG = Stock1(Temp,*a) 
  else:
    N  = len(Temp)
    dG = np.zeros(N,dtype=float)
    i = 0
    for T in Temp:
      dG[i] = Stock1(T,*a) 
      i = i+1
  return(dG)
def Stock1(T,a,b,c,d,e):
  dG1 = -(a/T + b*np.log(T) + c + d*T + e*T**2)*R*T
  return(dG1)

#==========================================================================
# read Sharp & Huebner data
#==========================================================================
SHnam   = np.loadtxt('SharpHuebner.dat',skiprows=0,usecols=[0],dtype='str')
SHcoeff = np.loadtxt('SharpHuebner.dat',skiprows=0,usecols=[1,2,3,4,5])

#==========================================================================
# read GGchem data
#==========================================================================
f = open('../../data/DustChem.dat','r')
lines = f.readlines()[:]
f.close
i = 0
lnam = []
lfit = []
for iline in range(0,999):
  if (i>=len(lines)): break
  line = lines[i]
  i = i+1
  if (line.find('[s]')>0):
    form = line.split()[0]
    name = line.split()[1]
    line = lines[i+1]
    Nel  = int(line.split()[0])
    i = i+Nel+2
    found = 0
    cfit = [0.0,0.0,0.0,0.0,0.0]
    print form,name
    for j in range(0,99):
      info = lines[i].split()
      i = i+1
      #print info
      if (len(info)<1): break
      if (info[0]=='2'): tmp=info[1:6]
      if (info[1]=='2'): tmp=info[2:7]
      if (info[0]=='2' or info[1]=='2'): found=1
    if (found==1):  
      for j in range(0,5):
        cfit[j] = float(tmp[j])
      print cfit
      lnam.append(form)
      lfit.append(cfit)
GGnam   = np.array(lnam,dtype='str')
GGcoeff = np.array(lfit)
print GGnam
print SHnam
#sys.exit()

#==========================================================================
# read SLOP16 database
#==========================================================================
elnam  = np.loadtxt('../../data/Abundances.dat',skiprows=5,usecols=[2],dtype='str')
elmass = np.loadtxt('../../data/Abundances.dat',skiprows=5,usecols=[3])
elnam  = np.array(elnam,dtype='str')
elmass = np.array(elmass,dtype='float')
amu = 1.66055E-24              # atomar mass unit [g]
NA  = 6.022140857E+23          # Avogadro's number [1/mol] 
f = open('slop16.dat','r')
lines = f.readlines()[355:]
f.close()
SLOPname = []
SLOPstoich = []
SLOPrho = []
Nline = len(lines)
for iline in range(0,Nline):
  line = lines[iline]
  if (line.find('gases')>=0): break
  if (line.find('ref:')<0): continue
  name   = lines[iline-2].split()[0]
  stoich = lines[iline-1].split()[1]
  vol    = lines[iline+1].split()[3]   #[cm3/mol]
  mass = 0
  warn = 0
  st   = stoich
  #print name,stoich
  for iel in range(0,9):
    ind1 = st.find('(')
    ind2 = st.find(')')
    if (ind1<0): break
    #print st,ind1,ind2
    elem = st[0:ind1]
    numb = st[ind1+1:ind2]
    st   = st[ind2+1:]
    if (elem==''): warn=1
    if (elem.find('Ru')>=0): warn=1
    if (elem.find('+')>=0): warn=1
    if (numb.find('N')>=0): warn=1
    if (warn==0):
      ind  = np.where(elnam==elem)[0][0]
      #print float(numb),elem,ind,elmass[ind]
      mass = mass + float(numb)*elmass[ind]
  if (warn==0):
    rho = mass*amu*NA/float(vol)  
    print name,stoich,"rho=",rho 
    SLOPname.append(name)
    SLOPstoich.append(stoich)
    SLOPrho.append(rho)
SLOPname   = np.array(SLOPname,dtype='str')
SLOPstoich = np.array(SLOPstoich,dtype='str')
SLOPrho    = np.array(SLOPrho,dtype='float')
#sys.exit()

#==========================================================================
# read SPRONSBL database
#==========================================================================
f = open('spronsbl.dat','r')
lines = f.readlines()[:]
f.close()
leer  = '                      '
iline = 0
Ncond = 0
Ncond2= 0
phase = 1
lform = []
lname = []
lNel  = []
lel   = []
lnum  = []
lq1   = []
lq2   = []
lcoef = []
lGdat = []
lrho  = []
all   = ''
for icond in range(0,9999):
  line1 = lines[iline]
  iline = iline+1
  if (line1.find('abandoned')>0): break
  if (line1.find('Landau')>0): phase=2
  if (line1[0]=='*'): continue
  line2 = lines[iline+0]
  line3 = lines[iline+1]
  line4 = lines[iline+2]
  line5 = lines[iline+3]
  line6 = lines[iline+4]
  line7 = lines[iline+5]
  iline = iline+6
  if (phase==2):
    line8 = lines[iline]
    iline = iline+1
  name  = line1.split()[0]
  form  = line1.split()[1]
  st    = line2.split()[1]
  stoich= st
  #print form,name
  Nel = 0
  el  = []
  num = []
  form2 = ''
  warn  = 0
  durch2 = 1
  durch3 = 1
  durch4 = 1
  durch5 = 1
  Natom  = 0
  for iel in range(0,9):
    ind1 = st.find('(')
    ind2 = st.find(')')
    if (ind1<0): break
    #print st,ind1,ind2
    elem = st[0:ind1]
    numb = st[ind1+1:ind2]
    el.append(elem)
    num.append(numb)
    if (int(float(numb))%2<>0): durch2=0
    if (int(float(numb))%3<>0): durch3=0
    if (int(float(numb))%4<>0): durch4=0
    if (int(float(numb))%5<>0): durch5=0
    st = st[ind2+1:]
    form2 = form2+elem
    if (numb<>'1'): form2 = form2+numb
    if (float(numb)<>int(float(numb))): warn=1  #broken stoichiometry
    if (float(numb)>16): warn=1                 #stoichiometic factor > 16
    if (elem == 'As'): warn=1
    if (elem == 'Ga'): warn=1
    Natom = Natom+int(float(numb))
    Nel = Nel+1

  #print Ncond,warn
  if (warn==1): continue
  Ncond = Ncond+1
  mult = 1
  if (durch2): mult=2
  if (durch3): mult=3
  if (durch4): mult=4
  if (durch5): mult=5
  if (mult>1):
    Natom = Natom/mult
    print form,name,"multi",mult
    el  = []
    num = []
    form2 = ''
    st = stoich
    for iel in range(0,Nel):
      ind1 = st.find('(')
      ind2 = st.find(')')
      if (ind1<0): break
      elem = st[0:ind1]
      numb = st[ind1+1:ind2]
      numb = str(int(float(numb)/mult))
      el.append(elem)
      num.append(numb)
      st = st[ind2+1:]
      form2 = form2+elem
      if (numb<>'1'): form2 = form2+numb
    print form2,el,num  

  # search for density in SLOP16.dat
  rho = 3.0  # default
  ind = np.where(SLOPstoich==stoich)[0]
  if (len(ind)>0): 
    rho = SLOPrho[ind[0]]
    print "density found in SLOP",name,SLOPname[ind[0]],rho
  else:  
    print "density not found",name
  if (all.find(' '+form2+' ')<0): Ncond2=Ncond2+1
  all = all+' '+form2+' '
  lrho.append(rho)
  lform.append(form2)
  lname.append(name)
  lNel.append(Nel)
  lel.append(el)
  lnum.append(num)

  #-------------------------------------
  # name of condensate and stoichiometry
  #-------------------------------------
  print form2,leer[len(form2):],name
  print rho
  el  = np.array(el)
  num = np.array(num,dtype='str')
  print Nel
  for iel in range(0,Nel):
    print("%2d %s" % (int(num[iel]),el[iel]))

  #--------------------------------------------------------
  # delta Gibbs free energy = G(cond,T)-G(ref,Tref) [J/mol]
  #--------------------------------------------------------
  dGref  = np.float(line4.split()[0])
  deltaH = np.float(line4.split()[1])   
  Sref   = np.float(line4.split()[2]) 
  Vref   = np.float(line4.split()[3]) 
  aa     = np.float(line5.split()[0]) 
  bb     = np.float(line5.split()[1])*1.E-5
  cc     = np.float(line5.split()[2]) 
  dd     = np.float(line5.split()[3]) 
  print dGref,deltaH,Sref,Vref
  print aa,bb,cc,dd
  Cp     = aa + bb*T1 + cc/T1**2 + dd/T1**0.5                  # [kJ/mol/K]
  CpdT   = aa*(T1 - Tref)                       \
         + bb/2.0*(T1**2 - Tref**2)             \
         - cc*(1.0/T1 - 1.0/Tref)               \
         + dd/0.5*(T1**0.5 - Tref**0.5)                        # [kJ/mol]
  CpdlnT = aa*np.log(T1/Tref)                   \
         + bb*(T1 - Tref)                       \
         - cc/2.0*(1.0/T1**2 - 1.0/Tref**2)     \
         - dd/0.5*(1.0/T1**0.5 - 1.0/Tref**0.5)                # [kJ/mol]
  dG_su  = dGref - Sref/1000*(T1-Tref) + CpdT - T1*CpdlnT      # [kJ/mol]
  dG_su  = dG_su*1000/mult                                     # [J/mol]
  ind = np.where(np.abs(T1-Tref)==np.min(np.abs(T1-Tref)))[0][0]
  malus = 0
  if (name=='CLINO-ENSTATITE'): malus=50*1000
  if (name=='LOW-TROILITE'): malus=50*1000
  lq1.append(dG_su[ind]/Natom+malus)

  #------------------------------------------------------------------------
  # subtract atomic delta Gibbs free energy = G(atom,T)-G(ref,Tref) [J/mol]
  #------------------------------------------------------------------------
  for iel in range(0,Nel):
    elname = el[iel]
    stoich = int(num[iel])
    ind    = np.where(atnam == elname)[0][0]
    coeff  = atdat[ind,:]
    #print elname,stoich,coeff
    dGatm  = Stock(T1,coeff)
    dG_su  = dG_su - stoich*dGatm

  #-------------------------------
  # fit as function of temperature
  #-------------------------------
  pfit,pcov = curve_fit(Stock1, T1, dG_su)
  print("fit coeffs: %12.5e %12.5e %12.5e %12.5e %12.5e" % 
      (pfit[0],pfit[1],pfit[2],pfit[3],pfit[4]))
  dG_su_fit = Stock(T1,pfit)
  N = len(T1)
  qual = 0.0
  for i in range(0,N):
    qual = qual + (dG_su[i]-dG_su_fit[i])**2
  qual = np.sqrt(qual/N)
  print "fit quality [kJ/mol]",qual/1000

  lq2.append(qual/1000)
  lcoef.append(pfit)
  lGdat.append(dG_su)
  
print
print Ncond2," condensates"
#sys.exit()
sys.stdout.write("now plotting ")

lform = np.array(lform,dtype='str')
lname = np.array(lname,dtype='str')
lNel  = np.array(lNel,dtype='int')
lq1   = np.array(lq1,dtype='float')
lq2   = np.array(lq2,dtype='float')
lcoef = np.array(lcoef,dtype='float')
lGdat = np.array(lGdat,dtype='float')
lrho  = np.array(lrho,dtype='float')
file  = open('DustChemSUPCRTBL.dat','w')
file.write("dust species\n")
file.write("============\n")
file.write("%i\n" % Ncond2) 
file2 = open('DustChemSUPCRTBL.tex','w')
index = np.argsort(lq1)
all   = ''
iplot = 0
Nout  = 0
SUnam = []
for icond in range(0,Ncond): 
  i = index[icond]
  name = lname[i]
  form = lform[i]
  Nel  = lNel[i]
  el   = lel[i]
  num  = lnum[i]
  sortq = lq1[i]
  qual = lq2[i]
  coeff= lcoef[i]
  dG_su= lGdat[i]
  rho  = lrho[i]
  if (all.find(' '+form+' ')>=0): 
    print form,name,"omitted"
    continue
  Nout = Nout+1
  all = all+' '+form+' '
  file.write("\n")
  file.write("%s[s]%s%s\n" %(form,leer[len(form):],name))
  if (rho==3.0):
    file.write("%5.2f     (estimated)\n" % rho)
  else:
    file.write("%5.2f\n" % rho)
  file.write("%i\n" % Nel)
  for iel in range(0,Nel):
    file.write("%2d %s\n"% (int(num[iel]),el[iel]))
  file.write("# Stock-fit to SUPCRTBL database (%i-%i)K, +/-%4.2f kJ/mol\n" \
        % (Tmin1,Tmax1,qual))
  file.write("  5 %15.8e %15.8e %15.8e %15.8e %15.8e \n" \
        % (coeff[0],coeff[1],coeff[2],coeff[3],coeff[4]))
  file2.write(" %3i & %15s & %20s & %5.1f & 5 & %13.6e & %13.6e & %13.6e & %13.6e & %13.6e & $\pm$%4.2f%s \n" \
              % (Nout,form,name,sortq/1000,coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],qual,'\\\\'))

  sys.stdout.write(".")
  sys.stdout.flush()

  iplot = iplot+1
  plt.figure(figsize=(7,8))
  plt.subplot(211)
  # the SUPCRTBL data
  dG_su_fit = Stock(T2,coeff)
  dGmean = dG_su_fit
  SUnam.append(form)
  
  # search for Sharp&Huebner data
  has_SH=0
  ind = np.where(SHnam==form)[0]
  if (len(ind)>0):
    has_SH=1
    iSH = ind[0]
    a0 = SHcoeff[iSH,0]
    a1 = SHcoeff[iSH,1]
    a2 = SHcoeff[iSH,2]
    a3 = SHcoeff[iSH,3]
    a4 = SHcoeff[iSH,4]
    Natom = 0
    for iel in range(0,Nel):
      Natom = Natom+int(num[iel])
    print form,"has Sharp & Huebner (1990) data",Natom
    dG_sh = (a0/T2 + a1 + a2*T2 + a3*T2**2 + a4*T2**3) * cal    # J/mol @ 1atm
    dG_sh = dG_sh + Natom*np.log(atm/bar)*R*T2                  # J/mol @ 1bar
    dGmean = dGmean + dG_sh
    plt.plot(T2,dG_sh/1000,c='orange',lw=4.5,label='Sharp & Huebner 1990')
    
  # search for GGchem data
  has_GG=0
  ind = np.where(GGnam==form+'[s]')[0]
  if (len(ind)>0):
    has_GG = 1
    iGG = ind[0]
    a0 = GGcoeff[iGG,0]
    a1 = GGcoeff[iGG,1]
    a2 = GGcoeff[iGG,2]
    a3 = GGcoeff[iGG,3]
    a4 = GGcoeff[iGG,4]
    print form,"has GGchem data"
    dG_gg = (a0/T2 + a1 + a2*T2 + a3*T2**2 + a4*T2**3)          # J/mol @ 1bar
    dGmean = dGmean + dG_gg
    plt.plot(T2,dG_gg/1000,c='green',lw=1.5,label='old GGchem')
    
  plt.plot(T2,dG_su_fit/1000  ,c='blue',lw=2.0,label='SUPCRTBL fit')
  plt.plot(T1,dG_su/1000 ,c='black',ls='--',lw=3.5,label='SUPCRTBL data')
  plt.xlim(0.0,1.05*Tmax2)
  plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=16)
  plt.ylabel(r'$\Delta_\mathrm{f} G^{\mathrm{1bar}}\mathrm{[kJ/mol]}$',fontsize=16)
  plt.title(form+"  -  "+name)
  plt.subplots_adjust(left=0.13, right=0.98, top=0.94, bottom=0.13)
  plt.legend(loc='upper left')

  if (has_GG+has_SH>0):
    plt.subplot(212)
    dGmean = dGmean/(1+has_GG+has_SH)
    plt.plot(T2,(dG_su_fit-dGmean)/1000,c='blue',lw=2.0,label='SUPCRTBL - MEAN')
    ymin = np.min((dG_su_fit-dGmean)/1000)
    ymax = np.max((dG_su_fit-dGmean)/1000)
    if (has_SH):
      plt.plot(T2,(dG_sh-dGmean)/1000,c='orange',lw=4.5,label='SH90 - MEAN')
      ymin = np.min([ymin,np.min((dG_sh-dGmean)/1000)])
      ymax = np.max([ymax,np.max((dG_sh-dGmean)/1000)])
    if (has_GG):
      plt.plot(T2,(dG_gg-dGmean)/1000,c='green',lw=1.5,label='GGCHEM - MEAN')
      ymin = np.min([ymin,np.min((dG_gg-dGmean)/1000)])
      ymax = np.max([ymax,np.max((dG_gg-dGmean)/1000)])
    ymin = np.min([ymin,-ymax])  
    ymax = np.max([ymax,-ymin])  
    plt.xlim(0.0,1.05*Tmax2)
    plt.ylim(2*ymin,2*ymax)
    plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=16)
    plt.ylabel(r'$\Delta_\mathrm{f} G^{\mathrm{1bar}}\mathrm{[kJ/mol]}$',fontsize=16)
    plt.legend(loc='lower center')

  plt.subplots_adjust(left=0.16,right=0.99,bottom=0.1,top=0.96,hspace=0.15)
  plt.savefig(pp,format='pdf')
  plt.clf()

SUnam = np.array(SUnam,dtype='str')
print SHnam
print SUnam
NN = len(GGnam)
for ii in range(0,NN):
  form = GGnam[ii]
  form = form.strip('[s]')
  print "searching for ",form
  a0 = GGcoeff[ii,0]
  a1 = GGcoeff[ii,1]
  a2 = GGcoeff[ii,2]
  a3 = GGcoeff[ii,3]
  a4 = GGcoeff[ii,4]
  dG_gg = (a0/T2 + a1 + a2*T2 + a3*T2**2 + a4*T2**3)          # J/mol @ 1bar
  dGmean = dG_gg
  has_SH=0
  ind1 = np.where(SHnam==form)[0]
  ind2 = np.where(SUnam==form)[0]
  if (len(ind1)>0 and len(ind2)==0):
    has_SH=1
    iSH = ind1[0]
    a0 = SHcoeff[iSH,0]
    a1 = SHcoeff[iSH,1]
    a2 = SHcoeff[iSH,2]
    a3 = SHcoeff[iSH,3]
    a4 = SHcoeff[iSH,4]
    Natom = 0
    for iel in range(0,Nel):
      Natom = Natom+int(num[iel])
    print form,"has Sharp & Huebner (1990) data",Natom
    dG_sh = (a0/T2 + a1 + a2*T2 + a3*T2**2 + a4*T2**3) * cal    # J/mol @ 1atm
    dG_sh = dG_sh + Natom*np.log(atm/bar)*R*T2                  # J/mol @ 1bar
    dGmean = dGmean + dG_sh
    iplot = iplot+1
    plt.figure(figsize=(7,8))
    plt.subplot(211)
    plt.plot(T2,dG_gg/1000,c='green',lw=1.5,label='old GGchem')
    plt.plot(T2,dG_sh/1000,c='orange',lw=4.5,label='Sharp & Huebner 1990')
    plt.xlim(0.0,1.05*Tmax2)
    plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=16)
    plt.ylabel(r'$\Delta_\mathrm{f} G^{\mathrm{1bar}}\mathrm{[kJ/mol]}$',fontsize=16)
    plt.title(form)
    plt.subplots_adjust(left=0.13, right=0.98, top=0.94, bottom=0.13)
    plt.legend(loc='upper left')

    plt.subplot(212)
    dGmean = dGmean/2
    plt.plot(T2,(dG_sh-dGmean)/1000,c='orange',lw=4.5,label='SH90 - MEAN')
    ymin = np.min((dG_sh-dGmean)/1000)
    ymax = np.max((dG_sh-dGmean)/1000)
    plt.plot(T2,(dG_gg-dGmean)/1000,c='green',lw=1.5,label='GGCHEM - MEAN')
    ymin = np.min([ymin,np.min((dG_gg-dGmean)/1000)])
    ymax = np.max([ymax,np.max((dG_gg-dGmean)/1000)])
    ymin = np.min([ymin,-ymax])  
    ymax = np.max([ymax,-ymin])  
    plt.xlim(0.0,1.05*Tmax2)
    plt.ylim(2*ymin,2*ymax)
    plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=16)
    plt.ylabel(r'$\Delta_\mathrm{f} G^{\mathrm{1bar}}\mathrm{[kJ/mol]}$',fontsize=16)
    plt.legend(loc='lower center')
    plt.savefig(pp,format='pdf')
    plt.clf()

file.close
file2.close

pp.close()
print ' '
print ' '
print 'written output to dG.pdf.'
