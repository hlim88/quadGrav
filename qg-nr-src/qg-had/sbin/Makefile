all: files

#FOR_FILE_PROC = f_file_processor.pl -free_format
F77_TO_FREE_FORM = ./f77tofreeform
FOR_FILE_PROC = ./f_file_processor.pl
FIXNUMS = ./f77_fix_numbers
FPP           = cpp


files:
	cat $(HEADDIR)/include/mem.inc | $(FOR_FILE_PROC) | $(FIXNUMS) > xmem.inc
	$(F77_TO_FREE_FORM) xmem.inc > $(HEADDIR)/include/f90mem.inc
	cat $(EQDIR)/$(SETUPFILE) | $(SETUPSCRIPT)
	$(FPP) util.x | $(FOR_FILE_PROC) > hyperutil.F90 
	$(FPP) params.x | $(FOR_FILE_PROC) > hyperparams.F90
	$(FPP) type.x | $(FOR_FILE_PROC) > hypertype.F90
	$(FPP) hyperio.x | $(FOR_FILE_PROC) > hyperio.F90
	mv param.inc .param.inc ; $(FPP) .param.inc | $(FOR_FILE_PROC) > param.inc ; rm .param.inc
	#/bin/rm util.x params.x type.x hyperio.x
	cp hyperutil.F90 hyperparams.F90 hypertype.F90 $(EQDIR)
	cp hyperio.F90 $(HEADDIR)/src/amr
	cp cctk.h fields.inc param.inc surfaces.inc $(HEADDIR)/include
	#cp cctk.h fields.inc param.inc $(HEADDIR)/include
	if test -f cpar.h; then\
		cp cpar.h fpar.h cfunc.h $(HEADDIR)/include;\
	fi
	if test -f cpar.c; then\
		cp cpar.c $(EQDIR);\
	fi
#	echo $(HAD_NUM_RK_ARRAYS)
#	cp $(HYPERRK_FILE) $(HEADDIR)/src/hyper/hyperrk.f90



clean:
	/bin/rm -f *.F90 *.inc *.F *.F90 *.f *.f90 *.x 
	/bin/rm -f cctk.h $(HEADDIR)/include/cctk.h
	/bin/rm -f xmem.inc $(HEADDIR)/include/f90mem.inc
	if test -f cpar.h; then\
		/bin/rm -f cpar.h fpar.h cfunc.h cpar.c;\
		/bin/rm -f $(HEADDIR)/include/cpar.h;\
		/bin/rm -f $(HEADDIR)/include/fpar.h;\
		/bin/rm -f $(HEADDIR)/include/cfunc.h;\
		/bin/rm -f $(EQDIR)/cpar.c;\
	fi

