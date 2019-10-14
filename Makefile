# source ./modules to set up module environment for build
EXE = condense mcgen fitjac testEsplines

# GNU compiler options
FC = gfortran
CC = gcc
CFLAGS    = -O2
FFLAGS    = -O2

# Includes, F90 Use files
F90MODS = -I./ -I$(NTCCHOME)/mod
INCLUDES = -I./ -I$(NETCDFHOME)/include

# Libraries
L_LAPACK = -L$(LAPACKHOME)/lib -llapack -lblas -lsmlib
L_TRANSP = -L$(NTCCHOME)/lib \
	-lplasma_state -lps_xplasma2 -lplasma_state_kernel \
	-lxplasma2 -lgeqdsk_mds -lmdstransp -lvaxonly -lnscrunch \
	-lfluxav -lr8bloat -lpspline -lezcdf \
	-llsode -llsode_linpack -lcomput -lportlib
L_NETCDF  = -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff
L_MDSPLUS = -L$(MDSPLUS_DIR)/lib -lMdsLib

fitjac: fitjac.c spline_interface.o
	$(CC) $(CFLAGS) -o $@ $< spline_interface.o -L${GSL_HOME}/lib -lgsl -lgslcblas -L${NETCDFC_HOME}/lib -lnetcdf -lm

condense: condense.c
	$(CC) $(CFLAGS) -o $@ $< -L${NETCDFC_HOME}/lib -lnetcdf

testEsplines: testEsplines.f90 spline_interface.o particleSplines.o
	$(FC) $(FFLAGS) testEsplines.f90 -o $@ \
	spline_interface.o particleSplines.o \
	-L${GSL_HOME}/lib -lgsl -lgslcblas -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $@ $< 

%.o: %.f90
	$(FC) -c $(FFLAGS) $(F90MODS) $(INCLUDES) -o $@ $<

all:  $(EXE)

mcgen: mcgen.o
	$(FC) -o $@ $< $(L_TRANSP) $(L_NETCDF) $(L_LAPACK) $(L_MDSPLUS)


clean:
	rm -f *.o $(EXE) *.mod

