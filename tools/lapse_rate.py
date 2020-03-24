import numpy as np

obj       = ['Earth','Venus','Mars','Jupiter','Titan']
LapseRate = [    6.5,    7.6,   4.3,      2.1,   1.0 ]  # K/km 
mu        = [   28.0,   44.0,  44.0,      2.3,  28.0 ]  # amu
g         = [   9.81,   8.87,  3.72,     24.8,  1.35 ]  # m/s2

amu       = 1.66055E-24              # g
km        = 1.E+5                    # cm
bk        = 1.380662E-16             # erg/K
g         = np.array(g)*100          # cm/s2
mu        = np.array(mu)*amu         # g
LapseRate = np.array(LapseRate)/km   # K/cm

print
print "object    kappa    gamma"
for i in range(0,5):
  kappa = LapseRate[i] * bk/mu[i] / g[i]
  gamma = 1/(1-kappa)
  print obj[i],kappa,gamma
