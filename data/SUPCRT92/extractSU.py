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
bk  = 1.38065812E-16    # Boltzman constant erg/K

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
SHatom = []
f = open('../../data/DustChem_SH90.dat','r')
lines = f.readlines()[:]
f.close
i=-1
j=0
SHnam = []
for line in lines:
  i=i+1
  if ('[s]' in line):
    Nel = int(float(lines[i+2].split()[0]))
    Nat = 0
    for el in range(0,Nel):
      Nat = Nat + int(float(lines[i+3+el].split()[0]))
    #print line.split()[0],Nel,Nat
    SHnam.append(line.split()[0])
    SHatom.append(Nat)
    j=j+1
SHnam = np.array(SHnam,dtype='str')

#==========================================================================
# read GGchem data
#==========================================================================
f = open('../../data/DustChem.dat','r')
lines = f.readlines()[:]
f.close
i = 0
lnam1 = []
lfit1 = []
lnam2 = []
lfit2 = []
GGall = []
for iline in range(0,999):
  if (i>=len(lines)): break
  line = lines[i]
  i = i+1
  if (line.find('[s]')>0 or line.find('[l]')>0):
    form = line.split()[0]
    name = line.split()[1]
    GGall.append(form)
    line = lines[i+1]
    Nel  = int(line.split()[0])
    i = i+Nel+2
    print form,name
    for j in range(0,99):
      found1 = 0
      found2 = 0
      info = lines[i].split()
      i = i+1
      #print info
      if (len(info)<1): break
      if (info[0]=='2'): tmp=info[1:6]
      if (info[1]=='2'): tmp=info[2:7]
      if (info[0]=='2' or info[1]=='2'): found1=1
      if (info[0]=='5'): tmp=info[1:6]
      if (info[1]=='5'): tmp=info[2:7]
      if (info[0]=='5' or info[1]=='5'):
        if (lines[i-2].find('NIST')>0): found2=1
      if (found1==1):  
        cfit = [0.0,0.0,0.0,0.0,0.0]
        for j in range(0,5):
          cfit[j] = float(tmp[j])
        print cfit
        lnam1.append(form)
        lfit1.append(cfit)
      if (found2==1):  
        cfit = [0.0,0.0,0.0,0.0,0.0]
        for j in range(0,5):
          cfit[j] = float(tmp[j])
        print cfit
        lnam2.append(form)
        lfit2.append(cfit)
GGnam1   = np.array(lnam1,dtype='str')
GGcoeff1 = np.array(lfit1)
GGnam2   = np.array(lnam2,dtype='str')
GGcoeff2 = np.array(lfit2)

#==========================================================================
# read Kitzmann data
#==========================================================================
f = open('../../data/Kitzmann2023/logK_condensates.dat','r')
lines = f.readlines()[:]
f.close
i = 0
lfor = []
lnam = []
lfit = []
lNat = []
lstoich = []
lTmelt = []
for iline in range(0,9999):
  if (i>=len(lines)): break
  line = lines[i]
  i = i+1
  if (line.find('(s)')>0 or line.find('(s,l)')>0):
    line = line.split()
    form = line[0]
    name = line[1]
    for j in range(2,9999):
      if (line[j]==":"): break
      name = name+line[j]
    if (form.find('.')>0): continue  
    Nat = 0
    stoich = []
    for k in range(1,9999):
      j = j+1
      if (line[j]=="#"): break
      j = j+1
      #Nat = Nat + int(line[j])
      Nat = Nat + 1
      stoich.append([line[j],line[j-1]])
    print form,name,Nat
    line = lines[i].strip()
    if (line=='s'):
      i=i+2
      tmp = lines[i].split()
      form = form[:-3]+"[s]"
      cfit = [0.0,0.0,0.0,0.0,0.0]
      for j in range(0,5):
        cfit[j] = float(tmp[j])
      lfor.append(form)
      lnam.append(name)
      lfit.append(cfit)
      lNat.append(Nat)
      lstoich.append(stoich)
      lTmelt.append(float(0.0))
    elif (line=='sl'):
      i=i+1
      Tmelt = float(lines[i].split()[0])
      i=i+1
      tmp = lines[i].split()
      form = form[:-5]+"[s]"
      cfit = [0.0,0.0,0.0,0.0,0.0]
      for j in range(0,5):
        cfit[j] = float(tmp[j])
      lfor.append(form)
      lnam.append(name)
      lfit.append(cfit)
      lNat.append(Nat)
      lstoich.append(stoich)
      lTmelt.append(Tmelt)
      i=i+1
      tmp = lines[i].split()
      form = form[:-3]+"[l]"
      cfit = [0.0,0.0,0.0,0.0,0.0]
      for j in range(0,5):
        cfit[j] = float(tmp[j])
      lfor.append(form)
      lnam.append(name+"(liquid)")
      lfit.append(cfit)
      lNat.append(Nat)
      lstoich.append(stoich)
      lTmelt.append(Tmelt)
    else:
      print "unknown list of condensates"
      print line
      stop
KZnam   = np.array(lfor,dtype='str')
KZname  = np.array(lnam,dtype='str')
KZstoich= lstoich
KZnel   = np.array(lNat)
KZcoeff = np.array(lfit)
KZtmelt = np.array(lTmelt)

#==========================================================================
# read BURCAT data
#==========================================================================
f = open('../../data/DustChem_BURCAT.dat','r')
lines = f.readlines()[:]
f.close
i = 0
lnam  = []
lNfit = []
lTfit = []
lcfit = []
for iline in range(0,9999):
  if (i>=len(lines)): break
  line = lines[i]
  i = i+1
  if (line.find('[s]')>0 or line.find('[l]')>0):
    form = line.split()[0]
    name = line.split()[1]
    line = lines[i+1]
    Nel  = int(line.split()[0])
    i = i+Nel+2
    found = 0
    #print form,name
    for j in range(0,99):
      info = lines[i].split()
      i = i+1
      if (len(info)<1): break
      #print info
      if (info[0]=='7'):
        ifit=int(float(info[1]))
        Tfit=np.zeros(ifit+1)
        for j in range(0,ifit+1):
          Tfit[j]=float(info[2+j])
        #print Tfit
        found=1
        break
    if (found==1):  
      print ifit
      cfit = np.zeros([ifit,14])
      for k in range(0,ifit):
        tmp = lines[i].split()
        #print tmp
        for j in range(0,14):
          cfit[k,j] = float(tmp[j])
        #print cfit[k,:]
        i=i+1
      lnam.append(form)
      lNfit.append(ifit)
      lTfit.append(Tfit)
      lcfit.append(cfit)
BUnam   = np.array(lnam,dtype='str')
BUNfit  = lNfit
BUTfit  = lTfit
BUcoeff = lcfit
#print BUnam
#print BUNfit
#print BUTfit
#print BUcoeff
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

#--------------------------------------
###   write DustChem_SUPCRTBL.dat   ###
#--------------------------------------
lname = np.array(lname,dtype='str')
lNel  = np.array(lNel,dtype='int')
lq1   = np.array(lq1,dtype='float')
lq2   = np.array(lq2,dtype='float')
lcoef = np.array(lcoef,dtype='float')
lGdat = np.array(lGdat,dtype='float')
lrho  = np.array(lrho,dtype='float')
file  = open('DustChem_SUPCRTBL.dat','w')
file.write("dust species\n")
file.write("============\n")
file.write("%i\n" % Ncond2) 
file2 = open('DustChem_SUPCRTBL.tex','w')
index = np.argsort(lq1)
all   = ''
Nout  = 0
SUnam = []
SUname = []
SUcoeff = []
SUdG = []
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
  SUnam.append(form+'[s]')
  SUname.append(name)
  SUcoeff.append(coeff)
  SUdG.append(dG_su)

file.close
file2.close

#--------------------------------------
###   write DustChem_Kitzmann.dat   ###
#--------------------------------------
NKZ = len(KZnam)
file  = open('../Kitzmann2023/DustChem_Kitzmann.dat','w')
file.write("dust species\n")
file.write("============\n")
file.write("%i\n" % NKZ)
for i in range(0,NKZ):
  file.write("\n")
  if (KZname[i].find('liquid')>0):
    file.write("%s   %s   %7.1f\n" %(KZnam[i],KZname[i],KZtmelt[i]))
  else:
    file.write("%s   %s\n" %(KZnam[i],KZname[i]))
  file.write(" 3.00     (estimated)\n")  
  file.write("%i\n" %(KZnel[i]))
  stoich = KZstoich[i]
  for j in range(0,KZnel[i]):
    st = stoich[j]
    file.write("%2i %s\n" %(int(st[0]),st[1]))
  bem = '  !!! NOT IN GGCHEM'
  search = np.str.upper(KZnam[i])
  for GGcond in GGall:
    cc = np.str.upper(GGcond)
    if (cc==search): bem=' '
  file.write("# dG-Stock-fit Kitzmann+2024  %s\n" %(bem))
  coeff = KZcoeff[i]
  file.write(" 5 %15.8e %15.8e %15.8e %15.8e %15.8e\n"
             %(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4]))
file.close()
stop           
           
SUnam = np.array(SUnam,dtype='str')
allcond = []
GGnamU1 = []
GGnamU2 = []
KZnamU = []
BUnamU = []
SUnamU = []
SHnamU = []
for cond in GGnam1:
  cc = np.str.upper(cond)
  GGnamU1.append(cc)
  if (cc not in allcond): allcond.append(cc)
for cond in GGnam2:
  cc = np.str.upper(cond)
  GGnamU2.append(cc)
  if (cc not in allcond): allcond.append(cc)
for cond in KZnam:
  cc = np.str.upper(cond)
  KZnamU.append(cc)
  if (cc not in allcond): allcond.append(cc)
for cond in BUnam: 
  cc = np.str.upper(cond)
  BUnamU.append(cc)
  if (cc not in allcond): allcond.append(cc)
for cond in SUnam: 
  cc = np.str.upper(cond)
  SUnamU.append(cc)
  if (cc not in allcond): allcond.append(cc)
for cond in SHnam: 
  cc = np.str.upper(cond)
  SHnamU.append(cc)
  if (cc not in allcond): allcond.append(cc)
Ncond = len(allcond)
allcond = np.array(sorted(allcond),dtype='str')
GGnamU1 = np.array(GGnamU1,dtype='str')
GGnamU2 = np.array(GGnamU2,dtype='str')
KZnamU  = np.array(KZnamU,dtype='str')
BUnamU  = np.array(BUnamU,dtype='str')
SUnamU  = np.array(SUnamU,dtype='str')
SHnamU  = np.array(SHnamU,dtype='str')
print "GGnam1=",GGnam1
print "GGnam2=",GGnam2
print "KZnam=",KZnam
print "BUnam=",BUnam
print "SUnam=",SUnam
print "SHnam=",SHnam
print "allcond=",allcond

for icond in range(0,Ncond): 

  name = ''
  form = allcond[icond]
  #if (form <> "AL2O3[S]"): continue
  pform = form
  dGmean = 0.0*T2
  ymin =  1.E+99
  ymax = -1.E+99

  # search for Sharp&Huebner data
  has_SH=0
  ind = np.where(SHnamU==form)[0]
  if (len(ind)>0):
    has_SH=1
    iSH = ind[0]
    pform = SHnam[iSH]
    a0 = SHcoeff[iSH,0]
    a1 = SHcoeff[iSH,1]
    a2 = SHcoeff[iSH,2]
    a3 = SHcoeff[iSH,3]
    a4 = SHcoeff[iSH,4]
    Natom = SHatom[iSH]
    print pform,"has Sharp & Huebner (1990) data",Natom
    dG_sh = (a0/T2 + a1 + a2*T2 + a3*T2**2 + a4*T2**3) * cal    # J/mol @ 1atm
    dG_sh = dG_sh + Natom*np.log(atm/bar)*R*T2                  # J/mol @ 1bar
    dGmean = dGmean + dG_sh
    ymin = np.min([ymin,np.min(dG_sh/1000)])
    ymax = np.max([ymax,np.max(dG_sh/1000)])
    
  # search for GGchem data
  has_GG1=0
  ind = np.where(GGnamU1==form)[0]
  if (len(ind)>0):
    has_GG1 = 1
    iGG1 = ind[0]
    pform = GGnam1[iGG1]
    a0 = GGcoeff1[iGG1,0]
    a1 = GGcoeff1[iGG1,1]
    a2 = GGcoeff1[iGG1,2]
    a3 = GGcoeff1[iGG1,3]
    a4 = GGcoeff1[iGG1,4]
    print pform,"has GGchem1 data"
    dG_gg1 = (a0/T2 + a1 + a2*T2 + a3*T2**2 + a4*T2**3)          # J/mol @ 1bar
    dGmean = dGmean + dG_gg1
    ymin = np.min([ymin,np.min(dG_gg1/1000)])
    ymax = np.max([ymax,np.max(dG_gg1/1000)])
  has_GG2=0
  ind = np.where(GGnamU2==form)[0]
  if (len(ind)>0):
    has_GG2 = 1
    iGG2 = ind[0]
    pform = GGnam2[iGG2]
    coeff = GGcoeff2[iGG2,:]
    print pform,"has GGchem2 data"
    dG_gg2 = Stock(T2,coeff)                                     # J/mol @ 1bar
    dGmean = dGmean + dG_gg2
    ymin = np.min([ymin,np.min(dG_gg2/1000)])
    ymax = np.max([ymax,np.max(dG_gg2/1000)])

  # search for SUPCRTBL data
  has_SU=0
  ind = np.where(SUnamU==form)[0]
  if (len(ind)>0):
    has_SU=1
    iSU = ind[0]
    pform = SUnam[iSU]
    name = SUname[iSU]
    print pform,"has SUPCRTBL data   ",name
    coeff = SUcoeff[iSU]
    dG_su = SUdG[iSU]
    dG_su_fit = Stock(T2,coeff)
    dGmean = dGmean + dG_su_fit
    ymin = np.min([ymin,np.min(dG_su_fit/1000)])
    ymax = np.max([ymax,np.max(dG_su_fit/1000)])

  # search for Kitzmann data
  has_KZ=0
  ind = np.where(KZnamU==form)[0]
  if (len(ind)>0):
    has_KZ = 1
    iKZ = ind[0]
    pform = KZnam[iKZ]
    if (name==''): name=KZname[iKZ]
    print pform,"has Kitzmann data"
    coeff = KZcoeff[iKZ,:]
    dG_KZ = Stock(T2,coeff)
    #logKs = a0/T2 + a1*np.log(T2) + a2 + a3*T2 + a4*T2**2
    #dG_KZ = -logKs*R*T2                                        # J/mol @ 1bar
    dGmean = dGmean + dG_KZ
    ymin = np.min([ymin,np.min(dG_KZ/1000)])
    ymax = np.max([ymax,np.max(dG_KZ/1000)])
    
  # search for BURCAT data
  has_BU=0
  ind = np.where(BUnamU==form)[0]
  if (len(ind)>0):
    has_BU = 1
    iBU = ind[0]
    pform = BUnam[iBU]
    Nfit = BUNfit[iBU]
    Tfit = BUTfit[iBU]
    coeff = BUcoeff[iBU]
    print pform,"has BURCAT data",Nfit,Tfit
    dG_bu = 0.0*T2
    i=0
    ifit=0
    for T in T2:
      if (T>Tfit[ifit+1] and ifit<Nfit-1): ifit=ifit+1
      #print T,ifit
      a0 = coeff[ifit,0]
      a1 = coeff[ifit,1]
      a2 = coeff[ifit,2]
      a3 = coeff[ifit,3]
      a4 = coeff[ifit,4]
      a5 = coeff[ifit,5]
      a6 = coeff[ifit,6]
      b0 = coeff[ifit,7]
      b1 = coeff[ifit,8]
      b2 = coeff[ifit,9]
      b3 = coeff[ifit,10]
      b4 = coeff[ifit,11]
      b5 = coeff[ifit,12]
      b6 = coeff[ifit,13]
      if (T>1000.0):
        dGRT = a0*(1-np.log(T)) -a1*T/2 -a2*T**2/6 -a3*T**3/12 - a4*T**4/20 \
             + a5/T - a6  
      else:
        dGRT = b0*(1-np.log(T)) -b1*T/2 -b2*T**2/6 -b3*T**3/12 - b4*T**4/20 \
             + b5/T - b6  
      dG_bu[i] = dGRT*(R*T)
      i += 1
    #dGmean = dGmean + dG_bu
    BUind = np.where((T2>0.9*Tfit[0]) & (T2<1.1*Tfit[Nfit]))
    ymin = np.min([ymin,np.min(dG_bu[BUind]/1000)])
    ymax = np.max([ymax,np.max(dG_bu[BUind]/1000)])

  if (has_SU+has_GG1+has_GG2+has_SH+has_KZ>0):
    dely = 0.05*(ymax-ymin)  
    ymin = ymin-dely
    ymax = ymax+dely
    plt.figure(figsize=(7,8))
    plt.subplot(211)
    if (has_SH): plt.plot(T2,dG_sh/1000,c='orange',lw=4.5,label='Sharp & Huebner 1990')
    if (has_GG1): plt.plot(T2,dG_gg1/1000,c='green',lw=1.5,label='old GGchem')
    if (has_GG2): plt.plot(T2,dG_gg2/1000,c='darkgreen',ls=':',lw=3,label='new GGchem')
    if (has_KZ): plt.plot(T2,dG_KZ/1000,c='magenta',lw=1.5,label='Kitzmann')
    if (has_SU): plt.plot(T2,dG_su_fit/1000,c='blue',lw=2.0,label='SUPCRTBL fit')
    if (has_SU): plt.plot(T1,dG_su/1000,c='black',ls='--',lw=3.5,label='SUPCRTBL data')
    if (has_BU): plt.plot(T2,dG_bu/1000,c='red',ls=':',lw=1.5)
    if (has_BU): plt.plot(T2[BUind],dG_bu[BUind]/1000,c='red',lw=2,label='BURCAT')
    plt.xlim(0.0,1.05*Tmax2)
    plt.ylim(ymin,ymax)
    plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=16)
    plt.ylabel(r'$\Delta_\mathrm{f} G^{\mathrm{1bar}}\mathrm{[kJ/mol]}$',fontsize=16)
    plt.title(pform+"  -  "+name)
    plt.subplots_adjust(left=0.13, right=0.98, top=0.94, bottom=0.13)
    plt.legend(loc='best')

  if (has_SU+has_GG1+has_GG2+has_SH+has_BU+has_KZ>1):
    plt.subplot(212)
    dGmean = dGmean/(has_SU+has_GG1+has_GG2+has_KZ+has_SH)    #+has_BU)
    ymin =  1.E+99
    ymax = -1.E+99
    if (has_SU):
      plt.plot(T2,(dG_su_fit-dGmean)/1000,c='blue',lw=2.0,label='SUPCRTBL - MEAN')
      ymin = np.min([ymin,np.min((dG_su_fit-dGmean)/1000)])
      ymax = np.max([ymax,np.max((dG_su_fit-dGmean)/1000)])
    if (has_SH):
      plt.plot(T2,(dG_sh-dGmean)/1000,c='orange',lw=4.5,label='SH90 - MEAN')
      ymin = np.min([ymin,np.min((dG_sh-dGmean)/1000)])
      ymax = np.max([ymax,np.max((dG_sh-dGmean)/1000)])
    if (has_GG1):
      plt.plot(T2,(dG_gg1-dGmean)/1000,c='green',lw=1.5,label='old GGCHEM - MEAN')
      ymin = np.min([ymin,np.min((dG_gg1-dGmean)/1000)])
      ymax = np.max([ymax,np.max((dG_gg1-dGmean)/1000)])
    if (has_GG2):
      plt.plot(T2,(dG_gg2-dGmean)/1000,c='darkgreen',ls=':',lw=3,label='new GGCHEM - MEAN')
      ymin = np.min([ymin,np.min((dG_gg2-dGmean)/1000)])
      ymax = np.max([ymax,np.max((dG_gg2-dGmean)/1000)])
    if (has_KZ):
      plt.plot(T2,(dG_KZ-dGmean)/1000,c='magenta',lw=1.5,label='Kitzmann - MEAN')
      ymin = np.min([ymin,np.min((dG_KZ-dGmean)/1000)])
      ymax = np.max([ymax,np.max((dG_KZ-dGmean)/1000)])
    if (has_BU):
      plt.plot(T2,(dG_bu-dGmean)/1000,c='red',lw=1.5,ls=':')
      plt.plot(T2[BUind],(dG_bu[BUind]-dGmean[BUind])/1000,c='red',lw=2,label='BURCAT - MEAN')
      #ymin = np.min([ymin,np.min((dG_bu[BUind]-dGmean[BUind])/1000)])
      #ymax = np.max([ymax,np.max((dG_bu[BUind]-dGmean[BUind])/1000)])
    ymin = np.min([ymin,-ymax,-1])  
    ymax = np.max([ymax,-ymin,+1])  
    plt.xlim(0.0,1.05*Tmax2)
    plt.ylim(2*ymin,2*ymax)
    plt.xlabel(r'$T\,\mathrm{[K]}$',fontsize=16)
    plt.ylabel(r'$\Delta_\mathrm{f} G^{\mathrm{1bar}}\mathrm{[kJ/mol]}$',fontsize=16)
    plt.legend(loc='best')

  if (has_SU+has_GG1+has_GG2+has_SH+has_KZ>0):
    plt.subplots_adjust(left=0.16,right=0.99,bottom=0.1,top=0.96,hspace=0.15)
    plt.savefig(pp,format='pdf')
    plt.clf()

  #if (icond==6):
  #  pp.close()
  #  stop

pp.close()
print ' '
print ' '
print 'written output to dG.pdf.'
