# source ./modules to set up module environment for build

ifndef GSL_HOME
 GSL_HOME = /p/transpusers/jbreslau/gsl/gsl-2.7
endif

EXE = condense mcgen fitjac testEsplines particle_sort

# GNU compiler options
FC = gfortran
CC = gcc
CFLAGS    = -Wall -O2
FFLAGS    = -Wall -O2

ifdef PARTICLE_MPI
PSCC = mpicc
TEFC = mpifort
PSFLAGS = -DPARTICLE_MPI
else
PSCC = $(CC)
TEFC = $(FC)
PSFLAGS =
endif

# Includes, F90 Use files
F90MODS = -I./ -I$(NTCC_HOME)/mod
INCLUDES = -I./ -I$(NETCDF_FORTRAN_HOME)/include

# Libraries
L_LAPACK = -L$(LAPACKHOME)/lib -llapack -lblas -lsmlib
L_TRANSP = -L$(NTCC_HOME)/lib \
	-lplasma_state -lps_xplasma2 -lplasma_state_kernel \
	-lxplasma2 -lgeqdsk_mds -lmdstransp -lvaxonly -lnscrunch \
	-lfluxav -lr8bloat -L$(PSPLINE_HOME)/lib -lpspline -lezcdf \
	-llsode -llsode_linpack -lcomput -lportlib
L_NETCDF  = -L$(NETCDF_FORTRAN_HOME)/lib -lnetcdf -lnetcdff
L_MDSPLUS = -L$(MDSPLUS_DIR)/lib -lMdsLib

fitjac: fitjac.c particleIO.o spline_interface.o
	$(CC) $(CFLAGS) -o $@ $< particleIO.o spline_interface.o \
	-L${GSL_HOME} -lgsl -L${NETCDF_C_HOME}/lib -lnetcdf -lm

condense: condense.c
	$(CC) $(CFLAGS) -o $@ $< -L${NETCDF_C_HOME}/lib -lnetcdf -lm

testEsplines: testEsplines.f90 spline_interface.o particleSplines.o
	$(TEFC) $(FFLAGS) testEsplines.f90 -o $@ \
	spline_interface.o particleSplines.o \
	-L${GSL_HOME} -lgsl -lm

particle_sort: particle_sort.c particleIO.o
	$(CC) $(CFLAGS) -o $@ $< particleIO.o -L${NETCDF_C_HOME}/lib -lnetcdf

spline_interface.o: spline_interface.c spline_interface.h
	$(CC) -c $(CFLAGS) -I$(GSL_HOME) -o $@ $<

particleSplines.o: particleSplines.c particleSplines.h spline_interface.h
	$(PSCC) -c $(CFLAGS) $(PSFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $@ $<

%.o: %.f90
	$(FC) -c -fopenmp $(FFLAGS) $(F90MODS) $(INCLUDES) -o $@ $<

all:  $(EXE)

mcgen: mcgen.o
	$(FC) -o $@ -fopenmp $< $(L_TRANSP) $(L_NETCDF) $(L_LAPACK) $(L_MDSPLUS)


clean:
	rm -f *.o $(EXE) *.mod

