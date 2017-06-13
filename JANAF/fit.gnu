set xlabel 'T[K]'             # for dG(T) 
set ylabel 'dG[kJ/mol]'
#set xlabel '5040/T[K]'       # for ln kp(5040/T)
#set ylabel 'ln ( gk[cgs] )'
set y2label ''
set xtics
set nox2tics
cal = 4.184

#---------------------------------------------------
### altes ln(gk) fuer CO zum Vergleich  ###   
a0 = -3.12249E+01  
a1 = +2.64659E+01 
a2 = -1.03320E-01  
a3 = +6.78494E-03 
a4 = -1.67181E-04
### altes ln(gk) fuer COS zum Vergleich  ###
a0 = -6.08816E+01  
a1 = +3.33362E+01 
a2 = -1.57416E-02  
a3 = +6.84627E-05  
a4 = -5.92278E-06
### altes ln(gk) fuer Mg+ zum Vergleich  ###
a0 = +2.32315E+01 
a1 = -1.96065E+01  
a2 = +2.84363E-01 
a3 = -1.87388E-02  
a4 = +4.53171E-04
### altes ln(gk) fuer H2SO4 zum Vergleich  ###
a0 = -1.79473E+02  
a1 =  5.90917E+01  
a2 =  3.41623E-02 
a3 = -5.42144E-03  
a4 =  1.57910E-04
### altes ln(gk) fuer C4N2 zum Vergleich  ###
a0 = -1.49717E+02  
a1 =  7.80084E+01  
a2 =  8.53628E-02 
a3 = -7.77629E-03  
a4 =  2.06956E-04
### altes ln(gk) fuer CH4 zum Vergleich  ###
a0 = -1.17097E+02  
a1 =  4.15062E+01 
a2 = -1.23514E-01  
a3 =  3.91271E-03 
a4 = -5.07497E-05
### altes ln(gk) fuer N+ zum Vergleich  ###
a0 =  2.33261E+01 
a1 = -3.54954E+01  
a2 =  2.26902E-01 
a3 = -1.23704E-02  
a4 =  2.46770E-04
#---------------------------------------------------
f(x) = a0 + a1*x**1 + a2*x**2  + a3*x**3 + a4*x**4 
#---------------------------------------------------

### altes dG fuer Fe[s] zum Vergleich  ###
b0 = +7.37828E+5
b1 = -4.22183E+5 
b2 = +1.71919E+2
b3 = -1.76037E-2
b4 = +2.31459E-6
#---------------------------------------------------
g(x) = (b0/x + b1 + b2*x + b3*x**2 + b4*x**3)/1000.0
#---------------------------------------------------

### altes dG fuer FeS[s] (Sharp & Huebner) ###
c0 = -6.83787E+4
c1 = -1.88620E+5
c2 = +6.70590E+1 
c3 = -2.318938E-3
c4 =  0.0
### altes dG fuer Na2SiO3[s] (Sharp & Huebner) ###
c0 =  4.63483E+4
c1 = -7.12298E+5
c2 = +2.06928E+2
c3 = -4.88925E-3
#-------------------------------------------------------
h(x) = (c0/x + c1 + c2*x + c3*x**2 + c4*x**3)/1000.0*cal
#-------------------------------------------------------

#plot 'fitout.dat' us 1:2  title "Daten"  with points, \
#     'fitout.dat' us 1:3  title "Fit"  with lines, \
#      f(x) title "old"  with lines
#pause mouse

plot 'fitout.dat'  us 1:2  title "Daten"  with points, \
     'fitout2.dat' us 1:2  title "Fit"  with lines
     


