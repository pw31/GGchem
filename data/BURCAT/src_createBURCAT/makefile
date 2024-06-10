F90          = ifx
F90LINKER    = ifx
#F90          = gfortran
#F90LINKER    = gfortran

# for IFORT-Compiler
#===================
DEFS     = -D IFORT -diag-disable 8290,8291 
###FFLAGS   = -r8 -i4 -O0 -g -traceback -check all,noarg_temp_created,nouninit -gen-interfaces -warn interfaces -fpe1 -mcmodel medium -shared-intel -qopenmp
FFLAGS  = -r8 -i4 -O3 -xHOST -prec-div -fp-model source -mcmodel medium -shared-intel -qopenmp
INCLUDES = 
LFLAGS   = $(FFLAGS)

# for gfortran-Compiler
#======================
#DEFS      = 
#FFLAGS    = -fdefault-real-8 -fdefault-double-8 -g -fbounds-check -fbacktrace 
#FFLAGS    = -fdefault-real-8 -fdefault-double-8 -g -O5
#INCLUDES  = 
#LFLAGS    = $(FFLAGS)

EXECS     = create_BURCAT

OBJECTS  = datamod.o \
           read_parameter.o \
	   is_nan.o \
	   init.o \
           nasa_polynomial.o \
           init_chemistry.o \
	   stindex.o \
           upper.o \
	   massfrac.o \
	   dispol.o \
           main.o

.SUFFIXES : .o .f .F .f90

default: main

main:    $(OBJECTS)
	 $(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS) 

clean:
	rm -f *.o 

realclean:
	rm -f *.o *~ *.mod *_genmod.f90 $(EXECS)

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

