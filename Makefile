
# Before making, you need to copy all healix headers in the current folder.
# I don't understand why this should be needed, but it is needed.
# cp /global/homes/e/eschaan/local/Healpix_3.31/include_f90/* ./

#default settings for ifort
F90C     = mpif90
#F90C    = ifort

healpix = $(HEALPIX)
LAPACKL = -mkl=sequential -lmpi -lhealpix -qopenmp

#remove -xHost if cluster is not homogeneous
#add -DHEALPIXI4B if using older healpix and get errors about arguments not matching
FFLAGS = -O3 -xHost -ip -fpp -error-limit 500 -DMPIPIX -DMPI -heap-arrays -g -traceback
#g and traceback allows for easier error handling
#cfitsio = /usr/local/Cluster-Apps/cfitsio/intel/3.300
#cfitsio = $(CFITSIO)
cfitsio = /global/common/sw/cray/cnl6/haswell/cfitsio/3.410/intel/17.0.2.174/xhct5xe

F90FLAGS = $(FFLAGS) -I$(INCLUDE) -I$(healpix)/include_f90 -L$(cfitsio)/lib -L$(healpix)/lib_f90 $(LAPACKL) -lcfitsio

OBJFILES= toms760.o inifile.o utils.o spin_alm_tools.o \
   HealpixObj.o HealpixVis.o

LENSPIX = $(OBJFILES) SimLens.o


default: simlens
#all: simlens recon

spin_alm_tools.o:  utils.o toms760.o
HealpixObj.o: spin_alm_tools.o
HealpixVis.o: HealpixObj.o
SimLens.o: HealpixVis.o inifile.o

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90


simlens: $(LENSPIX)
	$(F90C) -o simlens $(LENSPIX) $(F90FLAGS) $(LINKFLAGS)

recon: $(OBJFILES) LensReconExample.o
	$(F90C) -o recon $(OBJFILES) LensReconExample.o $(F90FLAGS) $(LINKFLAGS)

clean:
	rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il
	cp $(HEALPIX)/include_f90/* ./
