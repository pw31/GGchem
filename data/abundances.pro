readcol,"EarthCrust.dat",format='(I,A,A,F,F,F,F,D)',SKIPLINE=5, $
        nr,elname,el,f1,f2,f3,f4,mass

N=N_ELEMENTS(el)
nfrac=DBLARR(4,N)
for source=1,4 do begin
  if (source EQ 1) then f=f1
  if (source EQ 2) then f=f2
  if (source EQ 3) then f=f3
  if (source EQ 4) then f=f4
  nsum = 0.d0
  for i=0,N-1 do begin
    nsum = nsum + f[i]/mass[i]
  endfor
  for i=0,N-1 do begin
    nfrac[source-1,i] = f[i]/mass[i]/nsum
  endfor  
endfor

openw,1,"Abundances.dat"
for i=0,N-1 do begin
  printf,1,format='(I3,A14,A4,F10.5,4(E11.3))',$
         nr[i],elname[i],el[i],mass[i],nfrac[*,i]
endfor
close,1

end
