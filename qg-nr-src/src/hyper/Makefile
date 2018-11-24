######################################################################
#
#   $Id$
#   
#   Makefile for hyper
#   David Neilsen
#
######################################################################

SHELL = /bin/sh

.IGNORE =
.SUFFIXES : .F90 .f90

all: libhyper

ifeq ($(HYPER_FLUID),true)
  FLUID_MACRO=-DHYPER_FLUID
else
  FLUID_MACRO=
endif

ifeq ($(HYPER_BSSN),true)
  BSSN_MACRO=-DHYPER_BSSN
else
  BSSN_MACRO=
endif

OBJECTS = mask.o dissipation21.o hyperdissipation.o \
          hypercoords.o HyperInit.o \
          derivs43.o derivs43_pt.o derivs21.o derivs21_pt.o \
          derivs42.o derivs42_pt.o hyperderivs.o \
          hyperauxvars.o hyperfilereader.o charboundary.o \
          hyperboundary.o hyperrhs.o $(HYPERRK_FILES) \
          HyperAnalysis.o hyperfilereader1d.o hypercfl.o

#copy_auxvars_file:
#	if test ! -h $(HYPER_DIR)/hyperauxvars.o; then\
#		ln -s $(GRFLOWER_DIR)/global_auxvars.o $(HYPER_DIR)/hyperauxvars.o;\
#	fi


DERIVATIVE_FILES = mask.o dissipation21.o hypercoords.o derivs43.o \
                   derivs43_pt.o derivs42.o derivs42_pt.o \
                   derivs21.o derivs21_pt.o hyperderivs.o

derivs: $(DERIVATIVE_FILES) 
	echo "Compiling derivatives"

.f90.o:
	$(CPP) $(CPP_FLAG) $(MACRO) $(FLUID_MACRO) $(BSSN_MACRO) $*.f90 $*.f90.$(CPPFORTSUFFIX)
	$(F90_) -c $*.f90.$(CPPFORTSUFFIX)
	mv $*.f90.o $*.o
	rm -f $*.f90.$(CPPFORTSUFFIX)

.F90.o:
	$(F90_FIXED) $(F90SUFFIX) $(FLUID_MACRO) $(BSSN_MACRO) -c $*.F90


HYPERLIB_FULLNAME = libHyper.a
export HYPERLIB_FULLNAME

libhyper: $(OBJECTS)
	$(AR) $(HYPERLIB_FULLNAME) $(OBJECTS)
	mv $(HYPERLIB_FULLNAME) $(LIBDIR)


clean:
	/bin/rm -f *.o *.stb
	/bin/rm -f $(LIBDIR)/$(HYPERLIB_FULLNAME)
	if test -h hyperboundary.f90; then\
		/bin/rm hyperboundary.f90;\
	fi
	if test -h hyperauxvars.f90; then\
		/bin/rm hyperauxvars.f90;\
	fi
