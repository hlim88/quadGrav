.IGNORE:
.SUFFIXES : .F90 .f90 .F .f .cu .cu_o

all: $(EXECNAME)


OBJECTS    = hyperio.o main.o grid.o grid3.o updates.o \
             util.o mem.o problem.o  horizon.o \
             wenointerpolation.o \
             comm_mpi.o level.o level_lib.o grid_lib.o \
             comp_amr_error.o util2.o grid2.o level2.o grid_ell.o \
             horizon_apparent.o horizon_area2.o horizon_bstiso.o \
             horizon_bstksgen.o horizon_deriv.o \
             horizon_grid.o horizon_ilucg.o horizon_inter3d.o \
             horizon_interp.o horizon_io.o horizon_jacobian.o \
             horizon_metric.o horizon_hdf.o \
             horizon_metricdata.o horizon_molecule.o \
             horizon_solver.o horizon_sphere.o horizon_main.o \
             horizon_mijan.o gpu_neutrino.o gather.o cfl.o tracers.o

ifdef CUDA
CUFILES = Neutrino.cu
OBJECTS += $(patsubst %.cu,%.cu_o,$(notdir $(CUFILES)))
CUDA_INCLUDE = -I. -I/usr/local/cuda/include -I/home/matt/NVIDIA_CUDA_SDK/common/inc
CUDA_FLAGS = --compiler-options -fno-strict-aliasing -DUNIX -O3
CUDA_MACRO = -DCUDA
endif  

.cu.cu_o:
	nvcc -o $*.cu_o -c $*.cu $(CUDA_FLAGS) $(CUDA_INCLUDE)

.f.o:
	$(F90_FIXED) $(CUDA_MACRO) -c $*.f
.F.o:
	$(F90_FIXED) $(CPPMACRO)$(MACRO) $(CUDA_MACRO) $(EXTRA_MACRO) -c $*.F
.F90.o:
	$(F90_FIXED) $(F90SUFFIX) -c  $*.F90
.f90.o:
	$(F90_) -c $*.f90
.c.o:
	$(CC_)  -c $*.c


$(EXECNAME): $(OBJECTS) $(LIBDIR)/$(HYPERLIB_FULLNAME) $(LIBDIR)/$(EQLIB_FULLNAME)
	$(LOADER) $(OBJECTS) $(FLIBS) -o $(EXECNAME) 
	#$(F90_LOAD) $(OBJECTS) $(FLIBS) -o $(EXECNAME) -lgsl -lgslcblas
	mv $(EXECNAME) $(BINDIR)

Neutrino.cu_o: Neutrino.cu MersenneTwister.h MersenneTwister_kernel.cu Neutrino_kernel.cu
	nvcc -o Neutrino.cu_o -c Neutrino.cu $(CUDA_FLAGS) $(CUDA_INCLUDE)

actclean: 
	touch main.f comm_mpi.f grid.f

memclean: 
	touch problem.f mem.f main.f grid.F level.f grid3.F

clean: 
	 rm -f mysearch.timings junk? junk core.* *.stb *.o *.cu_o > /dev/null 2>&1
	/bin/rm -f *.mod
	/bin/rm -f grid_ell.f hyperio.F90 comp_amr_error.f

testchkpt: testcheckpt.f
	$(F77) -c testcheckpt.f 
	#$(F77) testcheckpt.o util2.o $(FLIBS) -o testcheckpt
	mpif90 -cxxlib  -O3 -mcmodel=large -i_dynamic -choicemod -r8 -I/usr/local/intel/include -I/home/steve/projects/hadSemiL2/src/hyperGHMHD/include -I/home/steve/projects/hadSemiL2/include -I/usr/local/gnu/include  -I/usr/local/include -I/usr/lib/include -module /home/steve/projects/hadSemiL2/include   -L//home/steve/projects/hadSemiL2/lib -lHyper -lhyperGHMHD -lHyper  -L/usr/local/gnu/lib -L/usr/lib/lib testcheckpt.o comm_mpi.o wenointerpolation.o hyperio.o main.o grid.o grid3.o updates.o util.o mem.o problem.o  horizon.o wenointerpolation.o comm_mpi.o level.o level_lib.o grid_lib.o comp_amr_error.o util2.o grid2.o level2.o grid_ell.o horizon_apparent.o horizon_area2.o horizon_bstiso.o horizon_bstksgen.o horizon_deriv.o horizon_grid.o horizon_ilucg.o horizon_inter3d.o horizon_interp.o horizon_io.o horizon_jacobian.o horizon_metric.o horizon_hdf.o horizon_metricdata.o horizon_molecule.o horizon_solver.o horizon_sphere.o horizon_main.o horizon_mijan.o gpu_neutrino.o gather.o cfl.o -L//home/steve/projects/hadSemiL2/lib -lHyper -lhyperGHMHD -lHyper   -lbbhutil  -L/home/install/Lorene/Lib -llorene_export -llorene -llorenef77 -L/usr/lib/lib -lgsl -lgslcblas  -llapack -lblas -lv3util -lvutil   -o testchkpt

