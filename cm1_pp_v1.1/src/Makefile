#################################################################################
#                         MAKEFILE FOR CM1_PP CODE                              #
#                                                                               #
#-------------------------------------------------------------------------------#
# Ryan Hastings, 7 Sep 2011                                                     #
#   Modified 23 Jan 2012                                                        #
# v0.8 28 Mar 2013                                                              #
#-------------------------------------------------------------------------------#
# v1.0, Ryan Hastings, 19 Jan 2017                                              #
#################################################################################
#################################################################################
# CONFIGURATION VARIABLES                                                       #

F90=gfortran # compiler to use
NCINC=-I/opt/local/include # location of include files for netcdf
NCLIB=-L/opt/local/lib # location of libraries for netcdf
FFTWINC=-I/usr/local/include # location of include files for fftw
FFTWLIB=-L/usr/local/lib # location of libraries for fftw

#################################################################################
# YOU PROBABLY DON'T WANT TO CHANGE ANYTHING BELOW THIS                         #
#################################################################################

ALL_PROGRAMS = dyn traj mass_flux circuit

CM1PP_code: $(ALL_PROGRAMS)

.PHONY : CM1PP_code

OBJD = globals.o nc_inout.o dyn.o
OBJT = globals.o nc_inout.o traj_utils.o traj.o
OBJM = globals.o nc_inout.o mass_flux.o
OBJC = globals.o nc_inout.o traj_utils.o circuit.o
OUTPUTOPT = -DNETCDF
LINKOPTS  = -lnetcdff

#################################################################################

dyn :	$(OBJD)
	$(F90) -o ../bin/dyn $(OBJD) -fcheck=bounds $(FFTWINC) $(FFTWLIB)-lfftw3 $(NCLIB) -lm $(LINKOPTS)

traj :	$(OBJT)
	$(F90) -o ../bin/traj $(OBJT) -fcheck=bounds $(NCLIB) -lm $(LINKOPTS)

circuit :	$(OBJC)
		$(F90) -o ../bin/circuit $(OBJC) -fcheck=bounds $(NCLIB) -lm $(LINKOPTS)

mass_flux :	$(OBJM)
	$(F90) -o ../bin/mass_flux $(OBJM) -fcheck=bounds $(NCLIB) -lm $(LINKOPTS)

dyn.o :		dyn.f90
		$(F90) -c dyn.f90 -fcheck=bounds $(FFTWINC) $(FFTWLIB) -lfftw3 -lm $(OPTS)

traj.o :	traj.f90
		$(F90) -c traj.f90 -fcheck=bounds $(OPTS)

circuit.o :	circuit.f90
		$(F90) -c circuit.f90 -fcheck=bounds $(OPTS)

mass_flux.o :	mass_flux.f90
	$(F90) -c mass_flux.f90 -fcheck=bounds $(OPTS)

nc_inout.o :	nc_inout.f90
		$(F90) -c nc_inout.f90 -fcheck=bounds $(OPTS) $(NCINC) $(NCLIB)

traj_utils.o :	traj_utils.f90
		$(F90) -c traj_utils.f90 -fcheck=bounds $(OPTS)

globals.o :	globals.f90
		$(F90) -c globals.f90 -fcheck=bounds $(OPTS)

clean:
		rm -rf *o *mod
