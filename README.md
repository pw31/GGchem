# GGchem

(c) Peter Woitke & Christiane Helling 2017

Fast thermo-chemical equilibrium code with or without equilibrium
condensation down to 100K.

Please cite our A&A paper "Equilibrium chemistry down to 100 K. 
Impact of silicates and phyllosilicates on carbon/oxygen ratio"; 
P. Woitke, Ch. Helling, G. H. Hunter, J. D. Millard, 
G. E. Turner, M. Worters, J. Blecic, J. W. Stock; 
2018; Astronomy & Astrophysics 614, 1; 
see GGchemPaper.pdf in this folder. 

We would be interested to hear from you about what kind of applications
you would like to use ggchem for, please let us know via email
Peter Woitke (pw31@st-and.ac.uk) 
as well as if you have any questions or problems. 

If your research results in any publications, please cite the above
article and consider to give us co-author-ship.

### To checkout the git repository and compile the code, use 

> git clone https://github.com/pw31/GGchem  
> cd GGchem/src16  
> cp makefile.prodimo makefile  
> make  

The makefile.prodimo is for ifort compiler, adjust your own makefile if
you want to compile e.g. with gfortran.

### To run the code, type 

> cd ..  
> ./ggchem input/default.in

It will create the output file "Static_Conc.dat", which contains all
computed molecular, atom and ion particle densities, the electron
density, solid and liquid particle densities, and supersaturation
ratios:

  * Tg: gas temperature [K],  
  * nHtot: total hydrogen nuclei particle density [cm-3],  
  * pges: total gas pressure [dyn/cm2]  
  * el ... W: atomic particle densities log10(natom)[cm-3]  
  * {mol}: molecular particle densities log10(nmol)[cm-3]  
  * S{cond}: supersaturation ratios log10(S) [-] of condensates  
  * n{cond}: concentration of condensed units per H nuclues log10(ncond/nHtot) [-]  
  * eps{el}: remainig element abundances in the gas phase [-]  
  * dust/gas: dust to gas mass ratio log10(rho_dust/rho_gas) [-]  
  * dustVol/H: dust volume per H nucleus log10(dust_volume) [cm3]  

S{cond} and n{cond} are used in the header to distinguish between 
supersaturation ratio and concentration of condensed units, whereas
{mol} (without the leading "n") is a molecular particle density.

### To visualise the results, use e.g.

> python tools/Plot_T.py  
> evince ggchem.pdf &  

### Customise your own model

To create your own model, make a copy of default.in and customize it to 
tell GGchem what it should do. You can also look at some of the other *.in 
files to lean from examples. Select or deselect elements by modifying 
the first line, default choice is

 H He C N O Na Mg Si Fe Al Ca Ti S Cl K Li F P V Cr Mn Ni Zr W el

where "el" means to include atomic and molecular ions, and the
electron density as well, assuming charge equilibrium.  Molecules are
included if they are made of the selected elements, otherwise they
will be ignored.

Choose element abundances with parameter abund_pick. The default
choice is abund_pick=3 for solar abundances from Asplund et
al.(2009). There are additional pre-installed options to use data from
"Abundances.dat", including "EarthCrust" (abund_pick=1), "Ocean"
(abund_pick=2) and "Meteorites" (abund_pick=4) as listed in
"Abundances.dat".  If you want any other element abundances, use
(abund_pick=0) followed by a name of a custom file with abundances,
see, e.g. input/model_Crich.in.

Choose sources for equilibrium constants kp(T)-data, default choice is 
dispol_new.dat. There are 6 different fit-formulas implemented, see
details in src16/smchem16.f (function gk). Data files having kp-data
are in folder data:

dispol_StockKitzmann.dat               : 2008, Diplomarbeit TU Berlin  
dispol_StockKitzmann_withoutTsuji.dat  : same, without Tsuji refits  
dispol_BarklemCollet.dat               : 2016, A&A 588, A96  
dispol_SharpHuebner.dat                : 1990, ApJSS 72, 417  
dispol_Tsuji.dat                       : 1973, A&A 23, 411  
dispol_GGchem.dat                      : old NIST-Janaf fits  
dispol_fast.dat                        : 9-molecules from Heng&Tsai 2016  

You can use combinations by setting dispol_file, dispol_file2, 
dispol_file3, dispol_file4 in your MyModel.in file, in which case 
the latter have preference over the former, and will overwrite 
previous data.

Choose whether you want to constrain the pressure (model_pconst=.true.)
or the mass density (model_pconst=.false.).

You can run single point model (model_dim=0), linear track
(model_dim=1) or 2D coverage (model_dim=2). Set parameters Tmin, Tmax
and then pmax, pmin or nHmax, nHmin for model_pconst=.true. or
.false., respectively. In the default model_dim=1 mode, ggchem will
make a linear track in (logp, logT) parameter space with Npoints
points.

If you want to switch on equilibrium condensation, set
model_eqcond=.true. In that mode, the code will be much slower, and
also possibly unstable. Always start from large T and then lower T
SLOWLY with successive calls. The code will create and expand
"database.dat" automatically from the results of every successful
call, such that once you have filled in the (p,T)-plane with many
points, the results will be faster and more reliable. The Gibbs-free
energy data files are in folder data:

DustChem_GGchem.dat      : old GGchem NIST-Janaf fits  
DustChem_SUPCRTBL.dat    : dG-fits from the SUPCRTBL database  
      (Zimmer et al. 2016, Computers and Geosciences, 90, 97)  
DustChem.dat             : currently used collection from both  

The pure gas phase chemistry needs about 0.4 ms per call for T > 1000 K
(real*8 version) and about 3 ms per call for T < 1000 K (real*16
version). These time measurements are for 16 elements + charge. Time
requirement roughly scale as N^3, if N is the number of elements.
The equilibrium condensation code requires many calls of the gas-phase
equilibrium chemistry routine, and takes about 0.02-0.09 sec per call,
depending on how much useful information is found in database.dat.
 
