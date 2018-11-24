############################################################################
#  
#   $Id: Makefile 4663 2012-12-05 23:06:55Z dneilsen1 $
#
#   Master Makefile
#
#   David Neilsen, Louisiana State University
#
#   To compile:
#     (1) The following environment variables must be externall defined:
#         LSV       = -lsv   (if installed, else blank)
#         IMPI      = -I/usr/local/vmi/mpich/include   (on titan.ncsa.)
#                     -I/usr/local/include  (on most systems)
#         LMPI      = -lmpich (Mac,Linux) or -lmpi (sgi)
#         MPIF90    = mpif90  (Linux)
#                     or /usr/local/mpich/bin/mpif90 (Mac)
#                     or f90 (sgi)
#         LOCALIB   = -L/usr/local/lib  or -L/usr/local/lib32
#         F90FLAGS  = -r8 -O0 -g -check all -fpe0 -traceback -fp-model strict
#          F90FLAGS  = -r8 -O3 -xP -no-prec-div -ip -static-libgcc -static-intel
            
#         CCFLAGS   =
#
#         For Intel 8:
#
#         F90_FREEFORM_FLAG   "-free"
#         F90_FIXEDFORM_FLAG  "-fixed"
#
#
#     (2) The following variables inside this makefile must be defined.
#         
#         EQS       = name of directory in ./src that contains the equations,
#                     e.g., hyperCurvedWave.
#         EXECNAME  = name for the executable, e.g., hcwad
#         all       = define what should be compiled on typing 'make',
#                     e.g., all: dirs hcw
#
#         It may be easiest to set EQS and EXECNAME in your environment
#         at compile time.  Then you don't need to modify the Makefile
#         each time a different set of equations are used.
#
############################################################################

SHELL = /bin/sh

HYPER_FLUID=false
HYPER_BSSN=false

# Check to see if EQS is defined in the environment.  
#If not, choose something here.
ifndef EQS
EQS = hyperQuadG
endif

# Check to see if EXECNAME is defined in the environment.  
# If not, choose something here.
ifndef EXECNAME
EXECNAME = hquadg
endif

# If you want/need the associated rk2 and rk3 fields generated
# then uncomment the following:
#     (this would be required for certain projects w/ an elliptic problem)
#export HAD_HIMEM="on"

#
# For more general equations of state, uncomment:
#
#export GEN_EOS="UWM_cand_1"
#export GEN_EOS="UWM_cand_2"
#export GEN_EOS="UWM_cand_3"

# Projects that do not have generalized harmonic (GH)
# variables will not be able to compile some portions
# of the various horizon finder interfaces.
# Adding the macro -DGH indicates that GH variables
# are present.
ifeq ($(EQS),hyperGHMHD)
  EXTRA_MACRO = $(CPPMACRO)-DGH 
endif
ifeq ($(EQS),hyperGHFFMHD)
  EXTRA_MACRO = $(CPPMACRO)-DGH
endif


ifeq ($(EQS),hyperMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperMHD2)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperHD2)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperXMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperGHXMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperGHMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperGHFFMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperBSSN_MHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperBSSN_MHD2)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperBSSN_XMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperRMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperBSSN_RMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperAMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperMHDe)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperGHMHDe)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperBDMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperBDGHMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperCCZ4_XMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperCCZ4SF_RMHD)
  HYPER_FLUID=true
endif
ifeq ($(EQS),hyperCCZ4SF_XMHD)
  HYPER_FLUID=true
endif


ifeq ($(EQS),hyperBSSN_EM)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperBSSN_MHD)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperBSSN_MHD2)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperBSSN_XMHD)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperBSSN_RMHD)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperBSSN)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperDilaton)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperEMDA)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperBSSN_MHDe)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperBSSN_BHEM)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperBSSN_FF)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperBSSN_FFSF)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperCCZ4)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperCCZ4_FF)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperCCZ4SF)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperCCZ4_FFSF)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperCCZ4_XMHD)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperCCZ4SF_RMHD)
  HYPER_BSSN=true
endif
ifeq ($(EQS),hyperCCZ4SF_XMHD)
  HYPER_BSSN=true
endif


ifeq ($(CCTK_F77_QUAD_TYPE),real16)
  HYPER_HAVE_QUAD=true
else
  HYPER_HAVE_QUAD=false
endif

export HYPER_FLUID HYPER_BSSN

ifdef HAD_HIMEM
HYPERRK_FILES = rk3_himem.o rk3tvd_himem.o rk3tvd_fv_himem.o hyperrk_himem.o
else
HYPERRK_FILES = rk3.o rk3tvd.o rk3tvd_fv.o imex_rk3.o hyperrk.o
endif



# Check to see if module flag is set in the environment
ifndef MODULE_INC_DIR
MODULE_INC_DIR = -module 
endif

ifndef CPPFORTSUFFIX
CPPFORTSUFFIX = f
else
CPP_FLAG = -traditional-cpp
endif

ifndef RANLIB_FLAGS
RANLIB_FLAGS = 
endif
#
#
all: dirs hcw

export EQS EXECNAME CPPFORTSUFFIX CPP_FLAG HYPERRK_FILES EXTRA_MACRO 
export HYPER_HAVE_QUAD
export MODULE_INC_DIR

#
#  Directory definitions
#
HEADDIR   = $(PWD)
AMRDIR    = $(PWD)/src/amr
HYPERDIR  = $(PWD)/src/hyper
LIBDIR    = $(PWD)/lib
BINDIR    = $(PWD)/bin
EQDIR=$(HEADDIR)/src/$(EQS)
export HEADDIR HYPERDIR AMRDIR LIBDIR BINDIR EQDIR

#
#  Tools
#
AR        = ar r
RANLIB    = ranlib
#SETUPSCRIPT = $(HEADDIR)/sbin/hypersetup
ifndef HSCRIPT
SETUPSCRIPT = $(HEADDIR)/sbin/hsetup2
else
SETUPSCRIPT = $(HEADDIR)/sbin/$(HSCRIPT)
endif

ifdef GHMHD_INDEP_RES
SETUPFILE = setup_indep_res
else
SETUPFILE = setup
endif

ifeq ($(EQS),hyperBSSN_EM)
ifdef BSSN_EM_INDEP_RES
SETUPFILE = setup_indep_res
else
SETUPFILE = setup
endif
endif

ifdef HDF5_HOME
HDF5_INCLUDE = -I$(HDF5_HOME)/include
HDF5_LIB_DIR = $(HDF5_HOME)/lib
HDF5_LIBS    = -L$(HDF5_LIB_DIR) -lhdf5_fortran -lhdf5 -lsz
export HDF5_INCLUDE HDF5_LIB_DIR HDF5_LIBS
endif


ifdef GSL_DIR
GSL_FLAG = -DGSL
GSL = $(GSL_DIR)
GSL_LIBS = -L$(GSL)/lib -lgsl -lgslcblas
GSL_INCLUDE = -I$(GSL)/include
#GSL_LDFLAGS = -Wl,--rpath -Wl,$(GSL_DIR)/lib
GSL_LDFLAGS = -L$(GSL_DIR)/lib
endif

ifdef OSX_GSL_DIR
GSL_FLAG = -DGSL
GSL = $(GSL_DIR)
GSL_LIBS = -L$(GSL)/lib -lgsl -lgslcblas
GSL_INCLUDE = -I$(GSL)/include
endif

ifdef HOME_LORENE
LORENE_CXX_FLAGS = -I$(HOME_LORENE)/Export/C++/Include -I$(HOME_LORENE)/C++/Include -I$(HOME_LORENE)/C++/Include_extra -I$(EQDIR)/include -I$(HEADDIR)/include $(LOCALINCLUDE)
LORENE_LIBS =-L$(HOME_LORENE)/Lib -llorene_export -llorene -llorenef77 $(LFFTW) 
LORENE_MACRO = -DLORENE
endif

ifdef NEED_LORENE_NAMESPACE
LORENE_CXX_FLAGS += -DNEED_LORENE_NAMESPACE
endif

CXX_FLAGS = $(CXXFLAGS) $(LORENE_CXX_FLAGS)

ifeq ($(CXX_COMMAND),icpc)
  CXX_LINKER_FLAG = -cxxlib
endif

ifeq ($(CXX_COMMAND),xlC)
  CXX_LINKER_FLAG = -lC
  GSL_LDFLAGS =
endif

ifeq ($(CXX_COMMAND),mpixlcxx)
  CXX_LINKER_FLAG = -cxxlib
  GSL_LDFLAGS =
endif

ifdef GNU_LORENE
GNU_LORENE_LIBS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lgfortran
endif

ifdef CUDA
CUDA_LIBS = -L/usr/local/cuda/lib -L/home/matt/NVIDIA_CUDA_SDK/lib -L/home/matt/NVIDIA_CUDA_SDK/common/lib -lcudart -lGL -lGLU -lglut   -L/usr/local/cuda/lib -lcutil
endif


export AR RANLIB SETUPSCRIPT SETUPFILE GSL_FLAG RANLIB_FLAGS

#
#  Libraries and includes
#
#
# Environment variables (define in .cshrc):
#
#     LSV       = -lsv   (if installed, else blank)
#     IMPI      = -I/usr/local/vmi/mpich/include   (on titan.ncsa.)
#                 -I/usr/local/include  (on most systems)
#     LMPI      = -lmpich (Mac,Linux) or -lmpi (sgi)
#     MPIF90    = mpif90  (Linux)
#                 or /usr/local/mpich/bin/mpif90 (Mac) 
#                 or f90 (sgi)
#     LOCALIB   = -L/usr/local/lib  or -L/usr/local/lib32
#     F90FLAGS  =
#     CCFLAGS   =
#

MPIOPTIONS =
#MPIOPTIONS = -mpilog
#MPIOPTIONS = -mpitrace

EQLIB     = -l$(EQS)

ifeq ($(EQS),hyperMHDe)
EQLIB     = -l$(EQS) -lnuc_eos
endif
ifeq ($(EQS),hyperGHMHDe)
EQLIB     = -l$(EQS) -lnuc_eos
endif
ifeq ($(EQS),hyperXMHD)
EQLIB     = -l$(EQS) -lnuc_eos
endif
ifeq ($(EQS),hyperGHXMHD)
EQLIB     = -l$(EQS) -lnuc_eos
endif
ifeq ($(EQS),hyperBSSN_XMHD)
EQLIB     = -l$(EQS) -lnuc_eos
endif
ifeq ($(EQS),hyperCCZ4_XMHD)
EQLIB     = -l$(EQS) -lnuc_eos
endif
ifeq ($(EQS),hyperCCZ4SF_XMHD)
EQLIB     = -l$(EQS) -lnuc_eos
endif



HYPERLIB = -L/$(LIBDIR) -lHyper $(EQLIB) -lHyper 
INCDIRS  = $(IMPI) -I$(EQDIR)/include -I$(HEADDIR)/include \
            $(LOCALINCLUDE) $(HDF5_INCLUDE) \
           -I/usr/local/include $(GSL_INCLUDE)

INCFLAGS = $(INCDIRS) $(MODULE_INC_DIR)$(HEADDIR)/include 
ALLFLAGS = $(F90FLAGS) $(INCFLAGS)

#
# These latter three libraries are for the use of Matt's interpolator.
# The "-lsvml" likely isn't needed for newer versions of the intel compilers.
#
FLIBS     = $(HYPERLIB) $(LMPI) -lbbhutil $(LSV) $(LORENE_LIBS) $(LFFTW) $(GSL_LIBS) \
 $(HDF5_LIBS) -llapack -lblas -lv3util -lvutil $(GNU_LORENE_LIBS) $(CUDA_LIBS)

#FLIBS     = $(HYPERLIB) $(LMPI) -lbbhutil $(LSV) $(LORENE_LIBS) $(GSL_LIBS) -llapack -lblas -lv3util -lvutil $(GNU_LORENE_LIBS) $(CUDA_LIBS)
#FLIBS     = $(HYPERLIB) $(LMPI) -lbbhutil $(LSV) 
CC_       = $(CC)  $(CCFLAGS) $(LORENE_MACRO) $(INCDIRS)
CPP_COMMAND   = cpp
CPP       = $(CPP_COMMAND) $(INCDIRS) $(LORENE_MACRO)

F90_      = $(MPIF90) $(ALLFLAGS) $(F90_FREEFORM_FLAG)
F90_FIXED = $(MPIF90) $(ALLFLAGS) $(F90_FIXEDFORM_FLAG)

F90_LOAD  = $(MPIF90) $(CXX_LINKER_FLAG) $(MPIOPTIONS) $(ALLFLAGS)  $(HYPERLIB) $(LOCALIB) $(GSL_LDFLAGS)

CXX_LOAD  = $(CXX_COMMAND) $(CXXFLAGS) $(MPIOPTIONS) $(HYPERLIB) $(LOCALIB) $(GSL_LDFLAGS)

ifdef NEED_CXX_LOAD
LOADER = $(CXX_LOAD)
else
LOADER = $(F90_LOAD)
endif

export FLIBS CC_ F90_ LOADER F90_FIXED CPP CXX_FLAGS CXX_COMMAND




#
#  Compiling
#

dirs: $(LIBDIR) $(BINDIR)
$(LIBDIR):
	if test ! -d $(LIBDIR); then\
		mkdir $(LIBDIR);\
	fi
$(BINDIR):
	if test ! -d $(BINDIR); then\
		mkdir $(BINDIR);\
	fi

hcw:
	cd $(EQDIR) && $(MAKE) params
	cd $(HYPERDIR) && $(MAKE) derivs
	cd $(EQDIR) && $(MAKE)
	cd $(HYPERDIR) && $(MAKE)
	cd $(AMRDIR) && $(MAKE)

clean:
	cd $(EQDIR) && $(MAKE) clean
	cd $(HYPERDIR) && $(MAKE) clean
	cd $(AMRDIR) && $(MAKE) clean
	cd $(HEADDIR)/include && $(MAKE) clean
	cd $(BINDIR) && rm -f $(EXECNAME)

realclean:
	$(MAKE) clean
	rmdir $(LIBDIR)
