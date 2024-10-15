OBJ = module.o sparse.o pn_aux.o math_lib.o mgmres.o symmetry.o \
		  readHmnR.o inverse.o proteus.o \
       eigen.o ham_qlayer2qlayer.o psi.o rand.o \
  		 ham_slab.o ham_bulk.o ek_slab.o ek_bulk_polar.o ek_bulk.o \
       readinput.o fermisurface.o surfgreen.o surfstat.o \
  		 mat_mul.o ham_ribbon.o ek_ribbon.o \
       fermiarc.o berrycurvature.o \
  		 wanniercenter.o dos.o  orbital_momenta.o \
		 berry.o wanniercenter_adaptive.o \
  		 effective_mass.o findnodes.o \
	main.o
	

# compiler
F90  = gfortran -cpp 
CC = gcc -E 

INCLUDE = #-I${MKLROOT}/include 
WFLAG = #-nogen-interface
#OFLAG = -O3 -static-intel
OFLAG = -O3 -ffree-line-length-512  #-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
FFLAG = $(OFLAG) $(WFLAG)
LFLAG = $(OFLAG)

# ARPACK LIBRARY
ARPACK= #-L/usr/local/lib -larpack

# need to specify if -DSUPERLU is provided in the F90 
SUPERLULIB=#/usr/local/lib/libsuperlu.a


# dynamic linking
LIBS = ${ARPACK} -llapack -lblas
 
main : $(OBJ)
	$(F90) $(LFLAG) $(OBJ) -o pn.x $(LIBS)
	cp -f pn.x ../bin

.SUFFIXES: .o .f90 .c

.f90.o :
	$(F90) $(FFLAG) $(INCLUDE) -c $*.f90
.c.o:
	$(CC) $(CFLAGS) -I$(INCLUDE) -c $< $(VERBOSE)

clean :
	rm -f *.o *.mod *~ pn.x
