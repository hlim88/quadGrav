
! Variable timestep modification by Dominic Marcello 2009-08-04
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  comput_cfl:                                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comput_cfl
      implicit     none
      include     'mpif.h'
      include     'param.inc'
      include 	   'glob.inc'
      include     'grid_methods.inc'
      logical      trace
      integer      gi
      integer      finest_level
      integer      lev
      integer      ierr
      real         this_cs
      real 	   max_cs

#ifdef DYNAMIC_CFL
#warning 'Compiling with dynamic timestep'
      trace = .false. 



      lev = 0
      max_cs = 0.0
   5  continue
      if (level_return_existence(lev)) then
        gi = level_return_start(lev)
  10    continue
        if (grid_return_existence(gi)) then
          if (grid_is_local(gi)) then
            call grid_comp_cfl(this_cs,gi)
  	    max_cs = max( max_cs, this_cs )
          end if
          gi = grid_return_sibling(gi)
          goto 10
        end if
        lev = lev + 1
      endif
      this_cs = max_cs
      call mpi_allreduce( this_cs, max_cs, 1, 
     *     MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
      if (ierr .ne. MPI_SUCCESS) then
        write(*,*) 'comput_cfl: MPI Error' 
        call abort 
      endif
      if( max_cs .le. 0.0 ) then 
        write(*,*) 'comput_cfl: Maximum signal speed is', max_cs
        call abort
      endif
      lambda = cfl_lambda / max_cs 
      dt = lambda * h
      if( trace .eqv. .true. ) then
        write(*,*) 'comput_cfl: lambda = ', lambda
      endif

#endif

      return 
      end    ! END: comput_cfl
 ! --------------end modification by Dominic Marcello 2009-08-04


