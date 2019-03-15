;---------------------------------------------------------------
pro polyfit, T, a, f, df	; Function + partial derivatives
;---------------------------------------------------------------
  print,a
  f = a(0)/T + a(1) + a(2)*T + a(3)*T^2 + a(4)*T^3
  if N_PARAMS() ge 4 then begin	
    df = [[1.d0/T],[1.d0+0.d0*T],[T],[T^2],[T^3]]
  endif  
end

line1 = ''
line2 = ''
file1 = 'MgClF.txt'
file2 = 'Mgref.txt'
file3 = 'Cl2.txt'
file4 = 'F2.txt'
norm2   = 1.0
norm3   = 0.5
norm4   = 0.5
stoich2 = 1.0
stoich3 = 1.0
stoich4 = 1.0
;file1 = 'SO2.txt'
;file2 = 'Sref.txt'
;file3 = 'O2.txt'
;file4 = 'F2.txt'         ; not used
;norm2 = 1.0
;norm3 = 0.5              ; reference state is 1/2 O2
;norm4 = 0.5
;stoich2 = 1.0
;stoich3 = 2.0
;stoich4 = 0.0            ; not used
R = 8.31447                     ;J/K/mol
openr,1,file1
readf,1,line1
readf,1,line2
close,1
readcol,file1,T1,Cp1,S1,GHoverT1,HH1,dH1,dG1,logkf1,skipline=3
GHoverT1 = -GHoverT1
HH1 = HH1*1000.0
dH1 = dH1*1000.0
dG1 = dG1*1000.0

readcol,file2,T2,Cp2,S2,GHoverT2,HH2,dH2,dG2,logkf2,skipline=3
GHoverT2 = -GHoverT2*norm2
HH2 = HH2*1000.0*norm2
dH2 = dH2*1000.0*norm2
dG2 = dG2*1000.0*norm2
GHoverT2 = INTERPOL(GHoverT2,T2,T1)   ; interpolate to T1-temperatures
HH2 = INTERPOL(HH2,T2,T1)
dH2 = INTERPOL(dH2,T2,T1)
dG2 = INTERPOL(dG2,T2,T1)

readcol,file3,T3,Cp3,S3,GHoverT3,HH3,dH3,dG3,logkf3,skipline=3
GHoverT3 = -GHoverT3*norm3
HH3 = HH3*1000.0*norm3                ; reference state is 1/2 H2
dH3 = dH3*1000.0*norm3
dG3 = dG3*1000.0*norm3
GHoverT3 = INTERPOL(GHoverT3,T3,T1)   ; interpolate to T1-temperatures
HH3 = INTERPOL(HH3,T3,T1)
dH3 = INTERPOL(dH3,T3,T1)
dG3 = INTERPOL(dG3,T3,T1)

readcol,file4,T4,Cp4,S4,GHoverT4,HH4,dH4,dG4,logkf4,skipline=3
GHoverT4 = -GHoverT4*norm4
HH4 = HH4*1000.0*norm4                ; reference state is 1/2 H2
dH4 = dH4*1000.0*norm4
dG4 = dG4*1000.0*norm4
GHoverT4 = INTERPOL(GHoverT4,T3,T1)   ; interpolate to T1-temperatures
HH4 = INTERPOL(HH4,T3,T1)
dH4 = INTERPOL(dH4,T3,T1)
dG4 = INTERPOL(dG4,T3,T1)

N=N_ELEMENTS(T1)
for i=0,N-1 do begin
  if (ABS(T1(i)-298.15) LT 1.E-3) then begin
    Href1=dH1(i) 
    Href2=dH2(i) 
    Href3=dH3(i) 
    Href4=dH4(i) 
  endif  
endfor

S1test  = 1.0/T1*HH1 - GHoverT1                       ; works
dG1test = -R*T1*ALOG(10.0)*logkf1                     ; works
dgef    = GHoverT1 - stoich2*GHoverT2 - stoich3*GHoverT3 - stoich4*GHoverT4
dHref   = Href1    - stoich2*Href2    - stoich3*Href3    - stoich4*Href4
dG1test = dHref + T1*dgef                             ; works
dH1test = dHref + HH1 - stoich2*HH2 - stoich3*HH3 - stoich3*HH4     ; works
logkf1test = -dG1test/(R*T1)/ALOG(10.0)

openw,1,"test.txt"
printf,1,line1
printf,1,line2
N=N_ELEMENTS(T1)
for i=0,N-1 do begin
  ;print,T1(i),S1(i),S1test(i)
  print,T1(i),dG1(i)/1000,dG1test(i)/1000
  ;print,T1(i),dH1(i)/1000,dH1test(i)/1000
  printf,1,format='(F7.2,F8.3,F8.3,F9.3,F14.3,F10.3,F16.3,F14.3)',$
         T1(i),Cp1(i),S1(i),GHoverT1(i),HH1(i)/1000,$
         dH1test(i)/1000,dG1test(i)/1000,logkf1test(i)
endfor
close,1

;xx = T1*1.d0
;yy = dG1test/1000.d0/(R*T1)
;ww = 1.d0/T1^0.5
;pp = [1.d0,yy[0],1.d-2,1.d-6,1.d-10]
;sigma = 0.d0*pp
;dGfit = CURVEFIT(xx,yy,ww,pp,function_name='POLYFIT',$
;        /DOUBLE,itmax=1000,TOL=1.E-40,sigma,status=status)
;POLYFIT, T1,pp,dGfit
;print, '     status: ',status
;print, ' parameters: ',pp
;print, '      sigma: ',sigma
;for i=0,N-1 do begin
;  print, T1(i),dG1test(i)/1000,dGfit(i)*(R*T1(i))
;endfor

stop
end
