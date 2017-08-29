file1 = 'Cl.txt'
file2 = 'Cl2ref.txt'
stoich = 0.5
;file1 = 'Ni.txt'
;file2 = 'Niref.txt'
;stoich = 1.0
line1 = ''
line2 = ''
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
S2  = S2*stoich
GHoverT2 = -GHoverT2*stoich
HH2 = HH2*1000.0*stoich
dH2 = dH2*1000.0*stoich
dG2 = dG2*1000.0*stoich
S2  = INTERPOL(S2,T2,T1)
GHoverT2 = INTERPOL(GHoverT2,T2,T1)   ; interpolate to T1-temperatures
HH2 = INTERPOL(HH2,T2,T1)
dH2 = INTERPOL(dH2,T2,T1)
dG2 = INTERPOL(dG2,T2,T1)

Tref = 298.15
N=N_ELEMENTS(T1)
for i=0,N-1 do begin
  if (ABS(T1(i)-Tref) LT 1.E-3) then begin
    Href1=dH1(i) 
    Href2=dH2(i) 
    iref=i
  endif  
endfor

R = 8.31447                              ;J/K/mol
S1test  = 1.0/T1*HH1 - GHoverT1          ; works
dG1test = -R*T1*ALOG(10.0)*logkf1        ; works
dHref   = Href1    - Href2    
dgef    = GHoverT1 - GHoverT2 
dG1test = dHref + T1*dgef                ; works
dH1test = dHref + HH1 - HH2              ; works

x1 = -HH2 + T1*S2 - Tref*S2(iref)
x2 = -HH1 + T1*S1 - Tref*S1(iref)
dGatom1  = dG1 - dG2 - x1
dGatom2  = dG1(iref) - dG2(iref) - x2

openw,1,"atom.txt"
printf,1,line1
printf,1,line2
for i=0,N-1 do begin
  ;print,T1(i),S1(i),S1test(i)
  print,T1(i),dG1(i)/1000,dG1test(i)/1000,$
        dGatom1(i)/1000,dGatom2(i)/1000
  ;print,T1(i),dH1(i)/1000,dH1test(i)/1000
  printf,1,format='(F8.2,F6.1,F6.1,F10.3,F15.3,F10.3,F16.3,F14.3)',$
         T1(i),0.0,0.0,0.0,0.0,0.0,dGatom1(i)/1000,0.0
endfor
close,1
stop
end
