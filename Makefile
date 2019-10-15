# © Shahram @ 2019


CC       = mpicc        # sets the C-compiler (default)
CXX       = mpiCC       # sets the C++-compiler (default)

FC 	 = mpif90

OPTIMIZE = -02 -Wall  -g   # optimization and warning flags (default)

MPICHLIB = -lmpich

CC       =  /share/apps/intel/impi/4.0.1.007/intel64/bin/mpicc #mpicc  # /usr/local/x86_64/gnu/openmpi-1.4.5/bin/mpicc   # sets the C-compiler
CXX      =  mpiCC
#ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
#CC       =  mpiCC
#OPT     +=  -DX86FIX
#endif
OPTIMIZE =   -g -Wall        # -g -Wall   # -fopenmp 
# GSL_INCL =  -I/usr/common/pdsoft/include
# # GSL_LIBS =  -L/usr/common/pdsoft/lib
 GSL_INCL =  -I/share/apps/gsl/include    #-I/usr/local/x86_64/gnu/gsl-1.9/include      #-I/afs/mpa/home/volker/Libs/include
 GSL_LIBS =  -L/share/apps/gsl/lib -Wl,"-R /share/apps/gsl/lib" #-L/usr/local/x86_64/gnu/gsl-1.9/lib   #-L/afs/mpa/home/volker/Libs/lib
 FFTW_INCL=  -I/stuperm/stalei/fftw/include  #-I/lustre/projects/p004_swin/shared/my_fftw_install/include    #-I/afs/mpa/home/volker/Libs/include
 FFTW_LIBS=  -L/stuperm/stalei/fftw/lib -Wl,"-R /stuperm/stalei/fftw/lib"  #-L/lustre/projects/p004_swin/shared/my_fftw_install/lib     #-L/afs/mpa/home/volker/Libs/lib -Xlinker -R -Xlinker /afs/mpa/home/volker/Libs/lib
 MPICHLIB =
 #HDF5INCL = -I/stuperm/stalei/run/HDF5MPI/include -Wall 
HDF5INCL = -I/share/apps/hdf5/hdf5-1.8.20-intel/include -Wall -g
# -DH5_USE_18_API  # -I/share/apps/include  -Wall -DH5_USE_16_API  #-I/usr/local/x86_64/gnu/hdf5-1.8.9-openmpi-psm/include   #-I/afs/mpa/home/volker/Libs/include
 #HDF5LIB  = -L/stuperm/stalei/run/HDF5MPI/lib  -lhdf5_hl  -lhdf5 -lz -lm  
 HDF5LIB  = -L/share/apps/hdf5/hdf5-1.8.20-intel/lib  -lhdf5_hl  -lhdf5 -lz -lm
# -L/share/apps/lib -lhdf5_hl -lhdf5 -lz -lm  #-L/usr/local/x86_64/gnu/hdf5-1.8.9-openmpi-psm/lib -lhdf5 -lz  #-L/afs/mpa/home/volker/Libs/lib -lhdf5 -lz
# #IPM_INTEL = -L/share/apps/ipm/intel/lib -lipm  
 OPT     +=  -DOLD_HDF5
 OPT     +=  -DUNEQUALSOFTENINGS

OPT      += -D DoParallel
# configuration and options
OPT      += -D TestConfig
OPT      += -D PERIODIC
OPT      += -D PMGRID
OPT      += -D PLACEHIGHRESREGION
OPT      += -D ENLARGEREGION
OPT      += -D MULTIPLEDOMAINS
OPT      += -D PEANOHILBERT
OPT      += -D WALLCLOCK
OPT      += -D MYSORT
OPT      += -D DOUBLEPRECISION
OPT      += -D DOUBLEPRECISION_FFTW
OPT      += -D PARTICLE_TAGGING
OPT      += -D FOF
OPT      += -D FOF_PRIMARY_LINK_TYPES 
OPT      += -D FOF_SECONDARY_LINK_TYPES
OPT      += -D FOF_GROUP_MIN_LEN
OPT      += -D SUBFIND
OPT      += -D MAX_NGB_CHECK
OPT      += -D SUBFIND_SAVE_PARTICLELISTS
OPT      += -D ORDER_SNAPSHOTS_BY_ID
OPT      += -D LINKLENGTH
OPT      += -D NO_GAS_CLOUDS
OPT      += -D NOTEST_FOR_IDUNIQUENESS
OPT      += -D NO_ISEND_IRECV_IN_DOMAIN
OPT      += -D FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
OPT      += -D OUTPUTACCELERATION
OPT      += -D OUTPUTTIMESTEP



#
OPTIONS = $(OPTIMIZE) $(OPT)

FOPTIONS = $(OPTIMIZE) $(FOPT)

EXEC   = PtagPP
OBJS = GlobalVars.o main.o EndRun.o InitConfig.o ReadSage.o ReadTag.o PaintStars.o WriteTag.o WriteTaggedSnap.o HashTag.o AnalyzeTable.o
INCL =Makefile
INCL +=GlobalVars.h proto.h
CFLAGS = $(OPTIONS) $(GSL_INCL) $(HDF5INCL)
LIBS = -g $(HDF5LIB)

$(EXEC): $(OBJS) $(FOBJS)
	$(CC) $(OPTIMIZE) $(OBJS) $(FOBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)


$(FOBJS): $(FINCL)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC)





