import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
#ppp = PdfPages('VenusTstruc.pdf')
from scipy import interpolate

bar = 1.E+6       
atm = 1.013*bar
km  = 1.E+5

#-------- gas species ----------
nCO2   = 0* 96.5E-2   # 94  - 98
nN2    = 0*  3.5E-2   # 1.5 - 5.5
nSO2   = 0*  60.E-6   # 23 - 250
nH2O   = 0*  35.E-6   # 4 - 53
nOCS   = 0*  10.E-6   # 0.25 - 40
nCO    = 0*  26.E-6   # 1 - 36
nHe    = 0*  12.E-6   # 6 - 36 
nHF    = 0* 500.E-9   # 5 - 5000
nHCl   = 0* 420.E-9   # 350 - 470
nH2S   = 0*   3.E-6   # 3 - 80 uncertain!   
nH2    = 0*   0.E-6   # 0 - 10 
nS3    = 0*   5.E-9   # 0.04 - 21
nS4    = 0*   5.E-9   # 1 - 9
nH2SO4 = 0*   0.E-6   # 0 - 3
nNO    = 0*  5.5E-9   # 4 - 7
nO2    = 0* (16.E-6 + 20.E-6)   # 8 - 24
nSF6   = 0* 200.E-9   # 100 - 300

#-------- condensates ---------
nSiO2   = 0*100.0
nMg2SiO4= 0*0.0
nMgAl2O4= 0*0.0
nCaSiO3 = 0*0.0         # test CaCO3[s] + SiO2[s] <-> CaSiO3 + CO2 hypothesis
nCaCO3  = 0*0.0         
nFeS2   = 0*10.0        # test 3 FeS2[s] + 16 CO2 <-> Fe3O4[s] + 6 SO2 + 16 CO
nCaSO4  = 0*10.0
nFe2O3  = 0*0.0         # test 2 Fe3O4[s] + CO2 <-> 3 Fe2O3 + CO
nFe3O4  = 0*10.0
nFeS    = 0*0.0
nNaCl   = 0*0.0

nSiO2       = 0*110.0        # eq.cond. model 1
nMgSiO3     = 0*100.0
nCaMgSi2O6  = 0*24.0
nFeS2       = 0*23.0
nCaSO4      = 0*22.0
nFe3O4      = 0*21.0
nCaAl2Si2O8 = 0*20.0
nNaAlSi3O8  = 0*19.0
nKAlSi3O8   = 0*18.0
nCaTiSiO5   = 0*17.0

nSiO2       = 0*110.0        # eq.cond. model #2
nMgSiO3     = 0*100.0
nAl2O3      = 0*24.0
nCaSO4      = 0*23.0
nFeS2       = 0*22.0
nFe3O4      = 0*21.0
nCaAl2Si2O8 = 0*20.0
nNaAlSi3O8  = 0*19.0
nKAlSi3O8   = 0*18.0
nTiO2       = 0*17.0

#------------ oxide mass fractions ------------
massH  =  1.008
massHe =  4.003
massC  = 12.011
massN  = 14.007
massO  = 15.999
massF  = 18.998
massNe = 20.180
massNa = 22.990
massMg = 24.305
massAl = 26.982
massSi = 28.085
massS  = 32.060
massCl = 35.450
massAr = 39.948
massK  = 39.098
massCa = 40.078
massTi = 47.867
massMn = 54.938
massFe = 55.845
amount = 1000.0

nSiO2  = 0.0
nTiO2  = 0.0
nAl2O3 = 0.0
nFeO   = 0.0
nMnO   = 0.0
nMgO   = 0.0
nCaO   = 0.0
nNa2O  = 0.0
nK2O   = 0.0
nSO3   = 0.0

# nSiO2  = 45.1 *amount/(1*massSi + 2*massO)     # Venera 13
# nTiO2  = 1.59 *amount/(1*massTi + 2*massO)
# nAl2O3 = 15.8 *amount/(2*massAl + 3*massO)
# nFeO   = 9.3  *amount/(1*massFe + 1*massO)
# nMnO   = 0.2  *amount/(1*massMn + 1*massO)
# nMgO   = 11.4 *amount/(1*massMg + 1*massO)
# nCaO   = 7.1  *amount/(1*massCa + 1*massO)
# nNa2O  = 2.0  *amount/(2*massNa + 1*massO)
# nK2O   = 4.0  *amount/(2*massK  + 1*massO)
# nSO3   = 1.62 *amount/(1*massS  + 3*massO)

# nSiO2  = 48.7 *amount/(1*massSi + 2*massO)     # Venera 14
# nTiO2  = 1.25 *amount/(1*massTi + 2*massO)
# nAl2O3 = 17.9 *amount/(2*massAl + 3*massO)
# nFeO   = 8.8  *amount/(1*massFe + 1*massO)
# nMnO   = 0.16 *amount/(1*massMn + 1*massO)
# nMgO   = 8.1  *amount/(1*massMg + 1*massO)
# nCaO   = 10.3 *amount/(1*massCa + 1*massO)
# nNa2O  = 2.4  *amount/(2*massNa + 1*massO)
# nK2O   = 0.2  *amount/(2*massK  + 1*massO)
# nSO3   = 0.88 *amount/(1*massS  + 3*massO)

nH2O  = 0.00000001
nSiO2 = 40.0  *amount/(1*massSi + 2*massO) 
nMgO  = 40.0  *amount/(1*massMg + 1*massO)
nMgF2 = 20.0  *amount/(1*massMg + 2*massF)

#==========================================================================
epsH  = 2*nH2O + 2*nH2S + nHCl + nHF + 2*nH2 + 2*nH2SO4
epsHe = nHe
#epsNe = nNe
#epsAr = nAr
epsC  = nCO2 + nCO + nCaCO3 + nOCS
epsO  = 2*nCO2 + 2*nSO2 + nCO + nH2O + nNO + 4*nH2SO4 + nOCS + 3*nFe2O3 + 4*nFe3O4 + 4*nMg2SiO4 + 2*nO2 + 2*nSiO2 + 3*nCaSiO3 + 3*nCaCO3 + 4*nCaSO4 + 6*nCaMgSi2O6 + 3*nMgSiO3 + 8*nCaAl2Si2O8 + 8*nNaAlSi3O8 + 8*nKAlSi3O8 + 3*nAl2O3 + 4*nMgAl2O4 + nMgO + nCaO + nFeO + nMnO + 3*nSO3 + 2*nTiO2 + 5*nCaTiSiO5 + nK2O + nNa2O
epsN  = 2*nN2 + nNO
epsS  = nSO2 + nH2S + nH2SO4 + nOCS + 3*nS3 + 4*nS4 + nSF6 + nFeS + 2*nFeS2 + nCaSO4 + nSO3
epsCl = nHCl + nNaCl
epsF  = nHF + 6*nSF6 + 2*nMgF2
epsFe = nFeS + nFeS2 + 2*nFe2O3 + 3*nFe3O4 + nFeO
epsSi = nMg2SiO4 + nSiO2 + nCaSiO3 + 2*nCaMgSi2O6 + nMgSiO3 + 2*nCaAl2Si2O8 + 3*nNaAlSi3O8 + 3*nKAlSi3O8 + nCaTiSiO5
epsMg = 2*nMg2SiO4 + nCaMgSi2O6 + nMgSiO3 + nMgAl2O4 + nMgO + nMgF2
epsCa = nCaSiO3 + nCaCO3 + nCaSO4 + nCaMgSi2O6 + nCaAl2Si2O8 + nCaO + nCaTiSiO5
epsAl = 2*nCaAl2Si2O8 + nNaAlSi3O8 + nKAlSi3O8 + 2*nAl2O3 + 2*nMgAl2O4
epsNa = nNaAlSi3O8 + nNaCl + 2*nNa2O
epsK  = nKAlSi3O8 + 2*nK2O
epsTi = nTiO2 + nCaTiSiO5
epsMn = nMnO

print "H  ",12+np.log10(epsH /epsH)
print "He ",12+np.log10(epsHe/epsH)
print "C  ",12+np.log10(epsC /epsH)
print "N  ",12+np.log10(epsN /epsH)
print "O  ",12+np.log10(epsO /epsH)
print "F  ",12+np.log10(epsF /epsH)
print "S  ",12+np.log10(epsS /epsH)
print "Cl ",12+np.log10(epsCl/epsH)
print "Fe ",12+np.log10(epsFe/epsH)
print "Mn ",12+np.log10(epsMn/epsH)
print "Si ",12+np.log10(epsSi/epsH)
print "Mg ",12+np.log10(epsMg/epsH)
print "Ca ",12+np.log10(epsCa/epsH)
print "Al ",12+np.log10(epsAl/epsH)
print "Na ",12+np.log10(epsNa/epsH)
print "K  ",12+np.log10(epsK /epsH)
print "Ti ",12+np.log10(epsTi/epsH)

f = open('abund.in','w')
if (epsH>0):
  f.write("H  %.9f \n"%(12+np.log10(epsH /epsH)))
if (epsHe>0):
  f.write("He %.9f \n"%(12+np.log10(epsHe/epsH)))
if (epsC>0):
  f.write("C  %.9f \n"%(12+np.log10(epsC /epsH)))
if (epsN>0):
  f.write("N  %.9f \n"%(12+np.log10(epsN /epsH)))
if (epsO>0):
  f.write("O  %.9f \n"%(12+np.log10(epsO /epsH)))
if (epsF>0):
  f.write("F  %.9f \n"%(12+np.log10(epsF /epsH)))
if (epsS>0):
  f.write("S  %.9f \n"%(12+np.log10(epsS /epsH)))
if (epsCl>0):
  f.write("Cl %.9f \n"%(12+np.log10(epsCl/epsH)))
if (epsFe>0):
  f.write("Fe %.9f \n"%(12+np.log10(epsFe/epsH)))
if (epsMn>0):
  f.write("Mn %.9f \n"%(12+np.log10(epsMn/epsH)))
if (epsSi>0):
  f.write("Si %.9f \n"%(12+np.log10(epsSi/epsH)))
if (epsMg>0):
  f.write("Mg %.9f \n"%(12+np.log10(epsMg/epsH)))
if (epsCa>0):
  f.write("Ca %.9f \n"%(12+np.log10(epsCa/epsH)))
if (epsAl>0):
  f.write("Al %.9f \n"%(12+np.log10(epsAl/epsH)))
if (epsNa>0):
  f.write("Na %.9f \n"%(12+np.log10(epsNa/epsH)))
if (epsK>0):
  f.write("K  %.9f \n"%(12+np.log10(epsK /epsH)))
if (epsTi>0):
  f.write("Ti %.9f \n"%(12+np.log10(epsTi/epsH)))
f.close()  
