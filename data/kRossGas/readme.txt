Dear Peter,

I computed the table you requested.

Solar Composition: Asplund et al. 2009

The table is arranged as follows

1st col: log10(T/K) Range:4-1.7, step=-0.01
2-142 col: log10(Kross) (only gas), corresponding to log10(nH)=6,20,step=0.1

I also attach a fORTRAN 90 routine that performes a bilinear interpolation in log10(T) and log10(nH), and returns interpolated log10(Kross)

ifort -o interp blinear.f90

