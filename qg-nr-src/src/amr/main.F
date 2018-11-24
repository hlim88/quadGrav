cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc    NLSM: evolve the nonlinear sigma model S^3 --> S^3                      cc
cc                                                                            cc
cc    Copyright Steve Liebling, Southampton College, LIU 2002                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program nlsm

      implicit none

      ! Include required system MPI include file:
      include 'mpif.h'

      ! Include local mpi variables:
      include 'mpi_stuff.inc'

      ! Include user defined parameters
      include 'param.inc'

      ! My globals:
      include 'glob.inc'

      ! My actions defined on levels/grids:
      include 'action.inc'

      ! Grid methods (functions)
      include 'grid_methods.inc'

      ! Surface stuff
      include 'surfaces.inc'

      ! Tracer stuff
      include 'tracers.inc'

      ! Local variables to main:
      integer  i,       ! generic loop counter
     *         level,
     *         lev_i,
     *         lev_j,
     *         finestlev,
     *         fileunit,
     *         proc,
     *         j
      integer  keep_chkpt_readstate

      integer  proc_pick
      external proc_pick

      ! Testing for unexpected truncation:
      real(kind=8) testreal1,testreal2,diffreal

      ! Code termination variables:
      real(kind=8) init_wtime, current_wtime, elapsed_wtime
      integer      rc
      logical      complete

      logical  next_level_exists, time_to_refine, time_to_state

      !
      ! General trace:
      !
      logical     ltrace
      parameter ( ltrace  = .false. )
      !
      ! Trace refinement:
      !
      logical     ltracer
      parameter ( ltracer = .false. )
      !
      ! Trace MPI calls:
      !
      logical     ltracem
      parameter ( ltracem = .false. )

      !
      ! Initialize MPI and associated variables/arrays:
      !
      if(ltrace)write(*,*)'main: Calling proc_init()'
      call proc_init()

      !
      ! Store initial wall-clock time:
      !
      init_wtime = MPI_WTIME()
      if(ltrace)write(*,*)'main: initial wall-clock time: ',init_wtime

      !
      ! Initialize tree and grid memory structures:
      !
      call init_grid_struct()

      !
      ! Truncation test:
      !
      testreal1 = 3.14159111111111111
      testreal2 = 3.14159111111111111d0
      diffreal  = testreal1-testreal2
      if (diffreal.ne.0) then
         if (myid.eq.master) then
            write(*,*) ' '
            write(*,*) 'main: You need to recompile to force the use'
            write(*,*) '      of 8byte doubles for all reals, otherwise'
            write(*,*) '      the results will suffer from incorrect'
            write(*,*) '      truncation effects.'
            write(*,*) ' '
            write(*,*) ' With Intel compilers, use the -r8 flag and for'
            write(*,*) '      AIX use -qautodbl=dbl4'                   
            write(*,*) ' '
            write(*,*) ' Quitting....'
         end if
         call my_exit('Needs recompilation.')
      end if

      !.................................................................
      if ( myid .eq. master ) then !......  MASTER .....
      !.................................................................
         !
         call get_params()
         !
         ! Read in params, pack them together, 
         ! and broadcast to other proc's:
         !
         call params_pack()
         !call params_dump()
         call send_actionAll(BCASTPAR,0,myid,0,0,0,0,0,0,0,0)
         call params_bcast()
         call surface_init()
         !
         ! Open file for run time statistics output:
         !
         fileout = 8
      open(fileout,file='Runtimestat',form='formatted',status='unknown')
         !
         if (chkpt_readstate.eq.0) then
            !
            ! Create coarse grid and refine:
            !
            lev_i = -1
            if(ltrace)write(*,*)'call level_apply refining level:',lev_i
            call level_apply(lev_i,REFINE)
            ! Only needed at initial time:
            call level_apply(lev_i,SYNCBND)
            lev_i = 0
            ! Variable timestep modification by Dominic Marcello 2009-08-03
            if( variable_timestep .eq. 1 ) then
               call level_apply(lev_i,CFL)
            endif
            ! --------------end modification by Dominic Marcello 2009-08-03
 25         if (level_return_existence(lev_i)) then
               if (ltrace) write(*,*) 'main: Refining level:',lev_i
               call level_apply(lev_i,SETMSK)
               if (use_mask.gt.0) call level_apply(lev_i,BHMASK)
               if (elliptic_solve .lt. 0 .or. elliptic_solve .eq.1) then
                  if (ltrace) write(*,*) 'main: Ellsolve on lev:',lev_i
                  call level_apply(lev_i,ELLSOLVE)
               endif
               ! Inject (this for elliptic solves and FV fields, and
               !         should have only roundoff difference otherwise)
               !if (lev_i.gt.0) call level_apply(lev_i,INJECT)
               if (shadow .gt. 0 .and. clusterreadwrite.ne.1) then
                  !
                  ! Only if we're using the shadow hierarchy:
                  ! 
                  if (lev_i .gt. 0) then
                     write(*,*) 'main: Using shadow hierarchy',lev_i
                     ! 
                     !
                     if(ltrace)write(*,*)'main: Comp. TRE on lev:',lev_i
                     ! Step the  *parent*  level  *once*:
                     if(ltrace)write(*,*)'   Take step on lev:',lev_i-1
                     if(ltrace)write(*,*)'  lev_i-1 from time: ',
     *                    grid_return_time(level_return_start(lev_i-1))
                     call level_apply(lev_i-1,STEP)
                     ! Step the  *current* level  *twice*:
                     if(ltrace)write(*,*)'   Take 1 step  on lev:',lev_i
                     call level_apply(lev_i,  STEP)
                     if (elliptic_solve.lt.0.or.elliptic_solve.eq.1)then
                        if(ltrace)write(*,*)'main:Ellsolve on lev',lev_i
                        call level_apply(lev_i,ELLSOLVE)
                     endif
                     if(ltrace)write(*,*)'   Increment lev count:',lev_i
                     call level_apply(lev_i,LEVINC)
                     if(ltrace)write(*,*)'   Take 1 step  on lev:',lev_i
                     call level_apply(lev_i,  STEP)
                     ! Compute difference via the computation of bounds:
                     if(ltrace)write(*,*)'   Compute bounds/error',lev_i
                     call level_apply(lev_i,  BOUNDS)
                     ! Reset times to zero:
                     if(ltrace)write(*,*)'   Reset times on levels. '
                     call level_apply(lev_i-1,RESETTIME)
                     call level_apply(lev_i,  RESETTIME)
                     ! Reinitialize the parent level and current level:
                     if(ltrace)write(*,*)'   Reinit both levels. '
                     call level_apply(lev_i-1,INIT)
                     call level_apply(lev_i,  INIT)
                     !
                     if(ltrace)write(*,*)'   Injecting level...'
                     call level_apply(lev_i,  INJECT)
                  end if
               end if
               if(ltrace)write(*,*)'   Calling refinement routine...'
               call level_apply(lev_i,REFINE)
               call level_apply(lev_i,SYNCBND)
               lev_i = lev_i + 1
               goto 25
            end if
            !
            if (elliptic_solve.lt.0 .or. elliptic_solve.eq.1) then
               if(ltrace)write(*,*)'main:Ellsolve last: ',lev_i-1
               call level_apply(lev_i-1,ELLSOLVE)
               ! Only calling this again to have the Hamiltonian 
               ! constraint evaluated again (be careful not to write
               ! over elliptic field(s)):
               if(ltrace)write(*,*)'main:Calling INIT : ',lev_i-1
               call level_apply(lev_i-1,INIT)
               do lev_j = 0, lev_i-1
                  call level_apply(lev_j,INIT)
               end do
            endif
            call level_apply(0,TRACERINIT)
            !call tracers_init()
            call level_tree_dump()
            call level_write_step( 0, 0.d0, 0)
            !
         else
            !
            ! Reading state from a checkpoint:
            !
            if (ltrace) write(*,*) 'main: Reading state'
            ! The value of this will be rewritten so keep it:
            keep_chkpt_readstate = chkpt_readstate
            call level_apply(0,READ_STATE)
            call level_tree_dump()
            call mem_stat()
            if (ltrace) write(*,*) 'main: Done. Reading state'
            if (ltrace) call params_dump()
            !
            if (keep_chkpt_readstate .gt. 10) then
               ! Redo parameters from par file so that they supercede:
               if (.true.) write(*,*)'main: Re-broadcasting parameters',
     *           ' from the par file'
               call send_actionAll(BCASTPAR,0,myid,0,0,0,0,0,0,0,0)
               if (ltrace) write(*,*)'main: Calling params_bcast()'
               call params_bcast()
               !    this has already been called!:
               !if (ltrace) write(*,*)'main: Calling surface_init()'
               !call surface_init()
               if (tracers_initial.ne. 0)call level_apply(0,TRACERINIT)
            end if
            !
      open(fileout,file='Runtimestat',form='formatted',status='unknown')
            lev_i = level_return_finest(Levelp,maxlev)+1
            !
         end if

            ! Output the initial data
            !   whether from regular setting of it
            !        or from reading state
      !
      ! Dangerous emulation?
      !
      if (nt0.eq.0 .and. emulate_proc.ge.0 .and. numprocs.eq.1) then
         write(*,*) 'main: Code is emulating another processor'
         write(*,*) 'main: Cross your fingers'
         write(*,*) 'main: Originally: myid = ',myid
         myid = emulate_proc
         write(*,*) 'main: And now:    myid = ',myid
      end if


            if (ltrace) write(*,*) 'main: Outputting coarse grid'
            do lev_j = 0, lev_i-1
               call level_apply(lev_j,OUTPUT)
            end do
            
      if (nt0.eq.0 .and. emulate_proc.ge.0 .and. numprocs.eq.1) then
          myid = 0
      end if
   
         if (ltrace) write(*,*) 'main: Beginning main evol loop'
         !----------------------------
         ! Loop:
         !----------------------------
         do steps = 1, nt0

            lev_i = 0

				! Variable timestep modification by Dominic Marcello 2009-08-03
				if( variable_timestep .eq. 1 ) then
	            call level_apply(lev_i,CFL)
				endif
				! --------------end modification by Dominic Marcello 2009-08-03

  5         continue
            if (ltrace) write(*,*) ' ---step      level:   ',lev_i
            call level_apply(lev_i,STEP)
            if (amrbound_prepost.eq.1 .and .amrbound_timealign.eq.1
     *           .and. lev_i.gt.0) then
               if (ltrace) write(*,*) ' ---A: get bounds on lev: ',lev_i
               call level_apply(lev_i,BOUNDS)
            end if
            if (ltrace) write(*,*) ' ---increment counter: ',lev_i
            call level_apply(lev_i,LEVINC)
            if (ltrace) write(*,*) '     l,c:',lev_i, lev_count(lev_i)
            if ( level_return_existence(lev_i+1)) then
               lev_i            = lev_i + 1
               if (amrbound_prepost.eq. 0) then
                  call level_apply(lev_i,BOUNDS)
                  if (ltrace) write(*,*) ' ---get bounds on lev: ',lev_i
               end if
               goto 5
            end if
            if (elliptic_solve .gt. 0 ) then
               ! At this point, we are on the finest level
               if (ltrace)write(*,*) 'calling ellsolve',lev_i
               call level_apply(lev_i,ELLSOLVE)
            end if
            call level_write_step( lev_i,
     *          grid_return_time(level_return_start(lev_i)),
     *                          lev_count(lev_i)            )
            if (ltrace) write(*,*) ' ---output sdf on lev: ',lev_i
            call level_apply(lev_i,OUTPUT)

         !
         ! REFINE
         !
   34       continue

         if (ltrace) write(*,*) 'main: Checking for collapse...'
         call proc_getfromall()
         call check_for_collapse()

         finestlev      = level_return_finest( Levelp, maxlev )
         if (refine_period_ctrl .eq. 0) then
           if(ltrace)write(*,*)'main: Rerefining after: ',refine_period
           time_to_refine = mod(lev_count(finestlev),refine_period).eq.0
         else
           if(ltrace)write(*,*)'main: Rerefining after: ',refine_deltat
           time_to_refine = mod(
     *          grid_return_time(level_return_start(finestlev)),
     *          refine_deltat ) .lt.
     *          lambda*grid_return_h(level_return_start(finestlev))
         end if
         if (time_to_refine) then
            lev_j = finestlev
            !
            ! Refine coarser grids advanced to this time step:
            !
  35        continue
            if (ltracer) write(*,*) 'main: Time to refine level:',lev_j
            if (ltrace)  write(*,*) 'main: Time to refine level:',lev_j
            call level_apply(lev_j,REFINE)
            call level_apply(lev_j+1,SETMSK)
            lev_j = lev_j - 1
            if ( lev_j .ge. 0) then
               if (refine_period_ctrl .eq. 0) then
                 time_to_refine=mod(lev_count(lev_j),refine_period).eq.0
               else
                 time_to_refine = mod(
     *              grid_return_time(level_return_start(lev_j)),
     *              refine_deltat ) .lt.
     *              lambda*grid_return_h(level_return_start(lev_j))
               end if
            else
               time_to_refine = .false.
            end if
            if (time_to_refine) goto 35
            !
            call level_tree_dump
            !
         end if
         !
         ! END: REFINE
         !


   10       if (lev_i .eq. 0) goto 100

            if (amrbound_prepost .eq. 0) then
               if (ltrace) write(*,*) ' ---get bounds on lev: ',lev_i
               call level_apply(lev_i,BOUNDS)
            end if
            if (ltrace) write(*,*) ' ---step      level:   ',lev_i
            call level_apply(lev_i,STEP)
            if (amrbound_prepost.eq.1 .and. amrbound_timealign.eq.0.and.
     *           mod(lev_count(lev_i),refine_factor)+1.eq.refine_factor)
     *          then
               if (ltrace) write(*,*) ' ---get bounds on lev: ',lev_i
               call level_apply(lev_i,BOUNDS)
            end if
            if (ltrace) write(*,*) ' ---increment counter: ',lev_i
            call level_apply(lev_i,LEVINC)
            if (ltrace) write(*,*) lev_i, lev_count(lev_i)
         !
         ! Three possibilities at this point:
         !   (1) this level has evolved to same time as parent level
         !         --- inject level to parent and step parent level
         !   (2) this level has fine child grids
         !         --- tranfer boundary values to fine level and step 
         !             fine level
         !   (3) this level has no finer grids, and we keep evolving
         !
         next_level_exists = level_return_existence(lev_i+1)
         if (     ( mod(lev_count(lev_i),refine_factor).eq.0 )
     *       .and.( .not. next_level_exists )                ) then

            if (ltrace) write(*,*) ' ---IF                 ',lev_i
   40       continue
            if (ltrace) write(*,*) ' ---output sdf on lev: ',lev_i
            ! Delay output until after the elliptic solve below
            !call level_apply(lev_i,OUTPUT)
            if (ltrace) write(*,*) ' ---inject from level: ',lev_i
            call level_apply(lev_i,INJECT)
            lev_i            = lev_i - 1
            if (ltrace) write(*,*) ' ---compute derv. lev: ',lev_i
            call level_apply(lev_i,COMPUT)
            if (ltrace) write(*,*) '     l,c:',lev_i, lev_count(lev_i)
            if (      (mod(lev_count(lev_i),refine_factor).eq.0)
     *          .and. (lev_i .ne. 0)                             ) then
               if (ltrace) write(*,*) ' Moving up a level:    ',lev_i
               goto 40 
            else if (elliptic_solve .gt. 0 ) then
               !
               ! At this point, we are on the finest level
               ! and have evolved such that all levels are
               ! time aligned, so we can do a global, elliptic solve:
               ! do an elliptic solve from finest to the coarse level: lev_i
               !
               if(ltrace)write(*,*) 'calling ellsolve',lev_i
               call level_apply(lev_i,ELLSOLVE)
            end if
            call level_write_step( lev_i,
     *          grid_return_time(level_return_start(lev_i)),
     *                          lev_count(lev_i)            )
            !
            ! Output on all levels after the elliptic solve (if one)
            !
            do j = level_return_finest( Levelp, maxlev ), lev_i, -1
               if (ltrace) write(*,*) ' ---output sdf on lev: ',j
               call level_apply(j,OUTPUT)
            end do
            !if (ltrace) write(*,*) ' ---output sdf on lev: ',lev_i
            !call level_apply(lev_i,OUTPUT)
         else if ( next_level_exists ) then
            if (ltrace) write(*,*) ' ---ELSE IF            ',lev_i
            lev_i = lev_i + 1
            goto 10
         else
            if (ltrace) write(*,*) ' ---ELSE               ',lev_i
            if (elliptic_solve .gt. 0 ) then
               ! At this point, we are on the finest level
               if(ltrace)write(*,*) 'calling ellsolve',lev_i
               call level_apply(lev_i,ELLSOLVE)
            end if
            call level_write_step( lev_i,
     *          grid_return_time(level_return_start(lev_i)),
     *                          lev_count(lev_i)            )
            if (ltrace) write(*,*) ' ---output sdf on lev: ',lev_i
            call level_apply(lev_i,OUTPUT)
         end if
         goto 34




 100        continue

         if (chkpt_period.gt.0) then
            finestlev      = level_return_finest( Levelp, maxlev )
            ! Have period w/r/t coarse level steps, not finest:
            time_to_state  = mod(lev_count(0),chkpt_period).eq.0
            !time_to_state  = mod(lev_count(finestlev),chkpt_period).eq.0
            ! Also save state if finished w/ last time step:
            if (time_to_state .or. steps.eq.nt0) then
            !if (time_to_state) then
               if (ltrace) write(*,*) 'main: Saving state...'
               call level_apply(finestlev,SAVESTATE)
            end if
         end if

         !
         ! Surface integrals:
         !
         if (asf_period.gt.0 .and. .false.) then
            if ( mod( lev_count(0), asf_period) .eq. 0 ) then
               if (ltrace)
     *            write(*,*) 'main: Computing Surface A functions...'
               call send_actionAll(SURFCOMP,0,myid,0,0,0,0,0,0,0,0)
               call surface_compute()
               if (ltrace)
     *            write(*,*) 'main: Outputing Surface A functions...'
               call surface_out(grid_return_time(1))
            end if
         end if

            !
            ! Check if it is time to terminate:
            !  --Check for existence of file HADCOMPLETE
            !  --Check if we've exceeded user-inputed runtime
            !
            complete = .false.
            if (.not.complete .and. run_wtime .gt. 0.d0) then
               current_wtime = MPI_WTIME()
               elapsed_wtime = current_wtime - init_wtime
               write(*,*)'main: elapsed wtime: ',elapsed_wtime,run_wtime
               if (elapsed_wtime.ge.run_wtime) then
                  complete = .true.
                  write(*,*)'main: Runtime reached: ',elapsed_wtime
               end if
            end if
            if (.not.complete) then
               open(unit=66,file='HADCOMPLETE',status='old',iostat=rc)
               if(rc .ne. 0) then
                  if(ltrace)write(*,*)'main: HADCOMPLETE not found'
               else
                  close(66)
                  complete = .true.
                  write(*,*)'main: HADCOMPLETE file found.'
               end if
            end if
            if (complete) then
               if (chkpt_period.gt.0) then
                  finestlev      = level_return_finest( Levelp, maxlev )
                  if (ltrace) write(*,*) 'main: Saving state...'
                  call level_apply(finestlev,SAVESTATE)
               end if
               !
               write(*,*)'main: Exiting...'
               call my_exit('Code completed.')
            end if

            !
            ! End of main evolution loop (over coarse level steps)
            !
         end do

         fileunit = 7
         open(fileunit, file='DISPERSAL',form='formatted')
         write(fileunit,*) amp,nt0
         close(fileunit)

      !.................................................................
      else  ! ............................  SLAVE  ....
      !.................................................................
         call listen_for_actions()
         if (ltracem) write(*,99) myid, ' main: Slave out of listen_for'
   
      end if

 99      format('[',I3,'] ',A,2I5)


      call my_exit('Normal completion')
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  my_exit:                                                                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine my_exit(message)
      implicit     none
      include     'glob.inc'
      include     'mem.inc'
      include     'param.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'action.inc'

      character*(*) message
      character(1)  junk
      integer      fileunit, finest,
     *             grid_return_finest,
     *             level_return_finest_ever,
     *             nextproc
      real(kind=8)       avgwork
      real(kind=8)       grid_return_time

      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,96)myid,'my_exit: Enter w/ msg:',message
      end if

      !
      ! If we're the master, command all processors to stop
      ! Else we tell the mater to quit, and the master will
      !     then tell everyone else
      !     (The issue here is that all processes need to
      !      get the QUIT command at a time when they 
      !      all expect a command, otherwise you can 
      !      get a process expecting data from a processor
      !      that already quit. This hangs the process since
      !      it won't get the QUIT command until it finishes
      !      with this data it will never get).
      !
      if (myid .eq. master) then
         if(ltrace)write(*,90)myid,'my_exit: As master,sendall'
         !write(*,*) 'myid = ',myid, ' sending all Quit'
         !call send_actionAllslow(QUIT,0,myid,0,0,0,0,0,0,0,0)
         call send_actionAll(QUIT,0,myid,0,0,0,0,0,0,0,0)
      else
         write(*,*) 'myid = ',myid, ' sending master Quit',message
         if(ltrace)write(*,96)myid,'my_exit: As master,sendall:',message
         call send_action(master,QUIT,0,myid,0,0,0,0,0,0,0,0)
         if (message .ne. 'Quitting from listen_for_actions.') then
            if(ltrace)write(*,90)myid,'my_exit: Not quitting yet'
            !write(*,*) 'myid = ',myid, ' Not quitting yet'
            return
         end if
      end if

      if (myid .eq. master) then
         write(*,*) '***** my_exit: *****'
         call proc_dump_wkld()
         call level_tree_dump()
         call grid_dump_freenums
         call mem_stat()
         call level_dump_count()

         write(*,98) myid,'              ................'
         write(*,98) myid,'              maxchi_thresh = ',maxchi_thresh
         write(*,98) myid,'              maxchi        = ',maxchi
         write(*,98) myid,'              maxchi_t      = ',maxchi_t
         write(*,99) myid,'              maxchi_l      = ',maxchi_l
         write(*,99) myid,'              maxchi_i,j,k  = ',maxchi_i,
     *                                        maxchi_j, maxchi_k
         write(*,98) myid,'              ................'
 90      format('[',I4,'] ',A)
 91      format('[',I4,'] ',A,I5)
 96      format('[',I4,'] ',2A)
 97      format('[',I4,'] ',A,F16.8)
 98      format('[',I4,'] ',A,3F16.8)
 99      format('[',I4,'] ',A,3I5)

c     if (myid .eq. master) then
         write(*,90)myid,'my_exit: Exiting...'
         if (trace_level .gt. 90) then
            write(*,91)myid,'my_exit: total_work     = ',total_work
            write(*,91)myid,'my_exit: steps          = ',steps
         end if
         if (steps.le.1) then
            avgwork = 0.d0
         else
            avgwork = total_work / (steps-1)
         end if
         write(*,97)myid,'my_exit: Avg work per coarse step  = ',avgwork

         !
         ! Dump exit status:
         !
!         write(*,*) 'my_exit: Dumping exit status to file: Exitstatus'
!         fileunit = 7
!      open(fileunit,file='Exitstatus',form='formatted',status='unknown')
!         junk = "Z"
!  6      format(A1)
!  8      read(fileunit,6,end=10) junk
!         goto 8
! 10      continue
!#ifdef __GFORTRAN__
!c        with gfortran 4.4 and later versions, I get the error 
!c        Fortran runtime error: Sequential READ or WRITE not allowed 
!c        after EOF marker, possibly use REWIND or BACKSPACE. So adding
!c        the backspace command here for gfortran.
!         backspace(fileunit)
!#endif
!         if (junk.eq."Z") then
!            write(fileunit,*) 't         L AvgWork    Amp          ',
!     *                     '       MaxChi  Maxchi_l  Mxc_t Mem  Msg'
!         end if
!         finest = grid_return_finest()
!         if (finest .le. 0) then
!            write(*,*)'my_exit: Problem with finest grid: ',finest
!            ! There should be at least one grid:
!            finest = 1
!         end if
!         write(fileunit,12) grid_return_time(finest),
!     *                   level_return_finest_ever(),
!c    *                   avgwork,amp,maxchi,maxchi_l,maxchi_t,
!c    *                   avgwork,amp_phi,maxchi,maxchi_l,maxchi_t,
!     *                   avgwork,0.d0   ,maxchi,maxchi_l,maxchi_t,
!     *                    (100.d0/arena)*max_qptr, message
! 12   format(F8.5,' ',I3,' ',F7.1,' ',F21.16,' ',G9.2,' ',I8,' ',F8.3,
!     *          ' ',F4.1,' ',A)
!         close(fileunit)

         write(*,96)myid,'my_exit: Msg: ', message
         write(*,90)myid,'***** my_exit  *****'
      else
         !
         ! All other process do this at exit:
         !
         call mem_stat()
      end if

      call gft_close_all()      ! close all sdf output files
      close(fileout)            ! close runtime statististics file
      call MPI_FINALIZE(ierr)
      stop

  9   format(A,F9.1)
      return
      end       ! END: my_exit

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  fint2:                                                                    cc
cc         Handler for interrupt 2 (signal interrupt).                        cc
cc         Gets registed in ctest.c to handle this interrupt.                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine fint2
       implicit   none
       include   'glob.inc'
       integer     THRESHOLD
       parameter ( THRESHOLD = 3 )

       interrupt_set = interrupt_set + 1

       if (trace_level.gt.0) write(*,9)
     * '            Interrupt detected...',interrupt_set,' / ',THRESHOLD

 9     format(A,I2,A,I2)
       if (interrupt_set .eq. 3)
     *    call my_exit('User selected exit.')

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  check_interrupt:                                                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine check_interrupt()
       implicit   none
       include   'param.inc'
       include   'glob.inc'
       include   'grid.inc'
       integer    user_input, user_input2, user_input3, i
       real(kind=8)     user_inputdp

       !
       ! If no interrupt, simply return
       !
       if (interrupt_set .eq. 0) return

       close(fileout)            ! close runtime statististics file

       !
       ! Figure out what user wants:
       !
  8    write(*,*) 'Enter option:'
       write(*,*) '1) Exit'
       write(*,*) '2) Resume'
       write(*,*) '3) Output run-time statistics'
       write(*,*) '4) Select functions for output'
       write(*,*) '5) Select output period by level'
       write(*,*) '6) Select output style'
       write(*,*) '7) Select threshold value for collapse detect'
       write(*,*) '8) Change tracing parameter: trace_level'
       write(*,*) '9) Output parameters'
       read(*,9,err=8) user_input
  9    format(I9)

       !...............=== 1...................
       if (user_input .eq. 1) then
          write(*,*) 'Exiting....'
          call my_exit('User selected exit.')
       !...............=== 2...................
       else  if (user_input .eq. 2) then
          write(*,*) 'Resuming...'
       !...............=== 3...................
       else if (user_input .eq. 3) then
          call level_tree_dump()
          call level_dump_count()
          call mem_stat()
         write(*,*) '                    ................'
         write(*,*) '                    maxchi_thresh = ',maxchi_thresh
         write(*,*) '                    maxchi        = ',maxchi
         write(*,*) '                    maxchi_t      = ',maxchi_t
         write(*,*) '                    maxchi_l      = ',maxchi_l
         write(*,*) '                    maxchi_i,j,k  = ',maxchi_i,
     *                                        maxchi_j, maxchi_k
         write(*,*) '                    ................'
          goto 8
       !...............=== 4...................
       else if (user_input .eq. 4) then
 11       write(*,*) '   Select grid function:'
          do i = 1, num_gfuncs
             write(*,12) i,gfunc_name(i),gfunc_out(i)
          end do
             write(*,13) num_gfuncs+1,'Go back'
 12       format('   ',I3,') ',A,' ',I2)
 13       format('   ',I3,') ',A)
          read(*,9,err=11) user_input2
          if ((user_input2.gt.0).and.(user_input2.le.num_gfuncs)) then
             gfunc_out(user_input2) = 1-gfunc_out(user_input2)
          else if (user_input2.eq.num_gfuncs+1) then
             goto 8
          else
             write(*,*) '  Invalid option.'
          end if
          goto 11
       !...............=== 5...................
       else if (user_input .eq. 5) then
 14       write(*,*) '   Select output level:'
          do i = 0, maxlev-1
             write(*,15) i,out_period(i)
          end do
 15       format('    Level: ',I3,'   output period: ',I3)
             write(*,*) '   ',maxlev,'  Go back'
          read(*,9,err=14) user_input2
          if ((user_input2.ge.0).and.(user_input2.lt.maxlev)) then
 16          write(*,*) '   Select period for Level:',user_input2
             read(*,9,err=16) user_input3
             if (user_input3.ge.0) then
                out_period(user_input2) = user_input3
             else
                write(*,*) 'Invalid option.'
                goto 16
             end if
          else if (user_input2.eq.maxlev) then
             goto 8
          else
             write(*,*) '   Invalid option.'
          end if
          goto 14
       !...............=== 6...................
       else if (user_input .eq. 6) then
 10       write(*,*) '   Select output style:'
          write(*,*) '   Current value output_style = ', output_style
          write(*,*) '   1) levels output separately (e.g. chi3.sdf)'
          write(*,*) '   2) levels output together   (e.g. chi.sdf)'
          write(*,*) '   3) Both (1) and (2)'
          write(*,*) '   4) Output given level and all finer ones'
          write(*,*) '   5) 1D SDF as function of x'
          write(*,*) '   6) Both (4) and (5)'
          write(*,*) '   7) Go back'
          read(*,9,err=10) user_input2
          if ((user_input2.gt.0).and.(user_input2.lt.7)) then
             output_style = user_input2
          else if (user_input2.eq.7) then
             goto 8
          else
             write(*,*) '  Invalid option.'
             goto 10
          end if
          goto 8
       !...............=== 7...................
       else  if (user_input .eq. 7) then
 17       write(*,*) 'Select value for threshold value:'
          write(*,*) '     maxchi_thresh = ',maxchi_thresh
          read(*,20,err=17) user_inputdp
 20       format(F20.8)
          if (user_inputdp.ge.0) then
             maxchi_thresh = user_inputdp
             write(*,*) '   maxchi_thresh = ',maxchi_thresh
          else
             write(*,*) '   Invalid option.'
             goto 17
          end if
          goto 8
       !...............=== 8...................
       else  if (user_input .eq. 8) then
 21       write(*,*) 'Select value for trace_value (0,100):'
          read(*,9,err=21) user_input2
          if (user_input2.ge.0) then
             trace_level = user_input2
             write(*,*) '   trace_level = ',trace_level
          else
             write(*,*) '   Invalid option.'
             goto 21
          end if
          goto 8
       !...............=== 8...................
       else  if (user_input .eq. 9) then
          call params_dump()
          goto 8
       !...............=== ?...................
       else  
          write(*,*) 'Invalid option.'
          goto 8
       end if

       !
       ! Reset interrupt
       !
       interrupt_set = 0

       !
       ! Re-open runtime statistics file:
       !
      open(fileout,file='Runtimestat',form='formatted',status='unknown')

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  init_grid_struct:                                                         cc
cc                    Initialize various quantities on every processor.       cc
cc                    Good place to put anything that needs to be done        cc
cc                    once on each proc before things get going.              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine init_grid_struct()
       implicit   none
       include   'glob.inc'
       include   'grid.inc'
       include   'mask.inc'
       integer    i, myid, tmpi
       ! For time=0, we get a bigger number and so have to use that:
       integer    time
       parameter (time = 0)
       !real(kind=8) time
       !parameter (  time = 1.d0)
       integer    mystringlength
       external   mystringlength
       integer    proc_return_myid
       external   proc_return_myid

      myid = proc_return_myid()

      !
      ! Initialize memory management variables:
      !
      qptr               = 1
      max_qptr           = qptr
      numgrids           = 0
      max_numgrids       = numgrids
      num_freegrids      = 0
      max_num_freegrids  = num_freegrids
      num_freeblocks     = 0
      max_num_freeblocks = num_freeblocks
      num_blockslost     = 0
      size_blockslost    = 0

      total_work         = 0.d0

      do i = 0, maxlev-1
         Levelp(i)     = -1.d0! init level pointers to null value
         lev_count(i)  = 0    ! init counter on each level to zero
         out_period(i) = 1    ! init output per. on lev. to every step
      end do

      num_injfields    = 0    ! Count up the number of injected fields
      num_synfields    = 0    ! Count up the number of synced fields
      !
      do i = 1, num_gfuncs
      !
         gfunc_out(i)  = 0    ! init flag for each grid func no output
         if (gfunc_type(i) .eq. GFUNC_INTEGRAL) then
            num_injfields = num_injfields + 1
         end if
         !
         ! Count number of sync'ed fields:
         !    (this "if" clause should match that in grid_syncbnd_slave)
         !    (I think whether time is 0 or not should produce the same
         !     count, but either way, we just need the bigger number
         !     to know how big to make the memory buffer to hold things):
         !
         if ((gfunc_type(i) .eq.GFUNC_INTEGRAL_ADV)   .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL  .and.(time.eq.0)).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DER)
     *                 .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DEREL)
     *                 .or.
     *       (gfunc_class(i) .eq. GFUNC_W)
     *                 .or.
     *       (gfunc_name(i) .eq. 'mask')
     *                 .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_ADVFV).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRALFV.and.time.eq.0).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DERFV) ) then
            num_synfields = num_synfields + 1
            !write(*,*)myid,'] synfield ',i,num_synfields
         end if
         !
      end do
      !write(*,*)myid,'] num_synfields=',num_synfields

      !
      ! Initialize tag with processor number for
      ! SDF output:
      !
      call int2str(myid,tag)

      !
      ! Initialize name for output file for storing state:
      !
      chkpt_fnum = 1
      !
      chkpt_fname(1) = tag(1:mystringlength(tag))//'hello.sdf'
      chkpt_fname(2) = tag(1:mystringlength(tag))//'goodbye.sdf'

      ! Mask stuff:
      do i = 1, max_num_masks
         bh_true(i) = .false.
      end do
      mask_level     = -1
      mask_iteration =  0
      forcedmerge    =  .false.

       ! Initialize all elements of the arena to zero
       if (.true.) then
          !write(*,*)myid,'] tag: ', tag
          !write(*,*)myid,'] chkpt_fname(1):',chkpt_fname(1)
          !write(*,*)myid,'] chkpt_fname(2):',chkpt_fname(2)
          !write(*,98)myid,' tag: ', tag
          !write(*,98)myid,' chkpt_fname(1):',chkpt_fname(1)
          !write(*,98)myid,' chkpt_fname(2):',chkpt_fname(2)
          !write(*,98)myid,' Zeroing out the arena'
          call load_scal1D(q,0.d0,arena)
          !write(*,*)myid,'] Done zeroing out the arena'
       end if

  98   format('[',I4,'] ',A,A)

       return
       end      ! END: init_grid_struct

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  get_params:                                                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine get_params()
       implicit     none
       include     'glob.inc'
       include     'param.inc'
       include     'largesmall.inc'
       character(64) paramfile
       integer      argc,
     *              iargc
       logical      double_equal
       external     double_equal

         !--------------------
         ! Get parameters:
         !--------------------
         argc = iargc()
         if ( argc.lt.1 ) then
            write(*,*) 'Usage: nlsm <parameter filename>'
            call my_exit('No parameter file specified.')
         else
           call getarg(1, paramfile)
         end if
         call  params_read(paramfile)
         !!!!!!!!!!!!!!!!!!
         !  Initializations:
         !!!!!!!!!!!!!!!!!!
         if (abs(hxyz0) .lt. SMALLNUMBER) then
         !if (hxyz0 .eq. 0) then
            hx     =  ( maxx0 - minx0 ) / (nx0-1.d0)
            hy     =  ( maxy0 - miny0 ) / (ny0-1.d0)
            hz     =  ( maxz0 - minz0 ) / (nz0-1.d0)
            dt     =  lambda * hx
            h      =  hx
            hxyz0  =  hx
         else
            hx     = hxyz0
            hy     = hxyz0
            hz     = hxyz0
            dt     = lambda * hx
            h      = hx
            minx0 = -hxyz0*(nx0-1.d0)/2.d0
            maxx0 = -minx0
            miny0 = -hxyz0*(ny0-1.d0)/2.d0
            maxy0 = -miny0
            minz0 = -hxyz0*(nz0-1.d0)/2.d0
            maxz0 = -minz0
         end if
         if (deriv_order.eq.4) then
            write(*,*)'main:  * * *'
            write(*,*)'main: This option does not work and has been'
            write(*,*)'main: disabled. Instead using 42'
            write(*,*)'main:  * * *'
            deriv_order = 42
         end if
         if (assume_symmetry.eq.3.or.assume_symmetry.eq.6) then
            write(*,*)'main:Because you want reflect.  symm. about z=0,'
            ! If enforcing symmetry at z=0:
            !minz0 = 0.d0
            ! Otherwise, extend boundary past z=0:
            minz0 = -4.d0*hz
            nz0   = (maxz0-minz0)/hz + 1
            if (mod(nz0,2).eq.0) nz0=nz0+1
            write(*,*) 'main: Forcing minz0 = ',minz0
            write(*,*) 'main:         maxz0 = ',minz0+(nz0-1)*hz
            write(*,*) 'main: and resetting nz0=',nz0
         end if
         if (assume_symmetry.eq.2.or.assume_symmetry.eq.6) then
            write(*,*)'main:Because you want reflect.  symm. about y=0,'
            ! Otherwise, extend boundary past z=0:
            miny0 = -4.d0*hy
            ny0   = (maxy0-miny0)/hy + 1
            if (mod(ny0,2).eq.0) ny0=ny0+1
            write(*,*) 'main: Forcing miny0 = ',miny0
            write(*,*) 'main:         maxy0 = ',miny0+(ny0-1)*hy
            write(*,*) 'main: and resetting ny0=',ny0
         end if
         if (assume_symmetry.eq.1.or.assume_symmetry.eq.6) then
            write(*,*)'main:Because you want reflect.  symm. about x=0,'
            ! Otherwise, extend boundary past z=0:
            minx0 = -4.d0*hx
            nx0   = (maxx0-minx0)/hx + 1
            if (mod(nx0,2).eq.0) nx0=nx0+1
            write(*,*) 'main: Forcing minx0 = ',minx0
            write(*,*) 'main:         maxz0 = ',minx0+(nx0-1)*hx
            write(*,*) 'main: and resetting nx0=',nx0
         end if
         if ( mod(refine_factor,2).ne.0 ) then
            write(*,*) 'main: refine_factor must be even'
            call my_exit('refine_factor must be even')
         end if
         if ( allowedl .ge. maxlev) then
            write(*,*) 'main: allowedl >= than maxlev'
            write(*,*) 'main: allowedl  = ',allowedl
            write(*,*) 'main: maxlev    = ',maxlev
            write(*,*) 'main: Setting allowedl equal to maxlev-1:'
            allowedl = maxlev - 1
            write(*,*) 'main: allowedl = ',allowedl
         end if
         ! Some parameter checking here:
         ! Some parameter checking here:
         ! Some parameter checking here:
         ! Some parameter checking here:
         ! Some parameter checking here:
         if (mindim .le. 2*bound_width) then
            write(*,*) 'main: The minimum size for any grid is'
            write(*,*) 'main: too low in comparison to the bound_width.'
            write(*,*) 'main: Re-setting it to twice the bound_width'
            write(*,*) 'main: old mindim = ',mindim
            write(*,*) 'main: bound_width = ',bound_width
            mindim = 2*bound_width
            write(*,*) 'main: new mindim = ',mindim
         end if
         if ( mod(ghostwidth,2) .eq. 0 ) then
            write(*,*) 'main: Ghostwidth needs to be odd.'
            write(*,*) 'main: Adding one to make it odd safely.'
            write(*,*) 'main: old ghostwidth = ',ghostwidth
            ghostwidth = ghostwidth + 1
            write(*,*) 'main: new ghostwidth = ',ghostwidth
         end if
         if ( ghostwidth .le. bound_width) then
            write(*,*) 'main: Ghostwidth generally needs to be'
            write(*,*) 'main: greater than the bound_width.'
            write(*,*) 'main: Re-setting to be twice the bound_width+1'
            write(*,*) 'main: old ghostwidth = ',ghostwidth
            ghostwidth = 2*bound_width
            write(*,*) 'main: new ghostwidth = ',ghostwidth
         end if
         !
         if ( .not.(double_equal(hx,hy).and.double_equal(hx,hz)) )then
         !if ( (hx.ne.hy).or.(hx.ne.hz) ) then
            write(*,*) 'main: unequal grid spacings not supported.'
            write(*,*) 'main: hx/y/z: ',hx,hy,hz
            write(*,*) 'main: hxyz0:  ',hxyz0
            call params_dump()
            call my_exit('unequal grid spacings no supported')
         else
            h = hx
            if (trace_level.gt.0) then
              write(*,*) 'main: **** BEGIN: Parameters ****'
              write(*,*) 'main: num time:   nt0 = ', nt0
              write(*,*) 'main:        allowedl = ', allowedl
              write(*,*) 'main:             nx0 = ', nx0
              write(*,*) 'main:             ny0 = ', ny0
              write(*,*) 'main:             nz0 = ', nz0
              write(*,*) 'main:             amp = ', amp
              write(*,*) 'main:           idata = ', idata
              write(*,*) 'main:   refine_factor = ', refine_factor
              write(*,*) 'main:    chkpt_period = ', chkpt_period
              write(*,*) 'main: chkpt_readstate = ', chkpt_readstate
              call mem_stat_brief()
            end if
            if (trace_level.gt.3) then
              write(*,*) 'main: grid spacing h = ', h
              write(*,*) 'main: time step   dt = ', dt
              write(*,*) 'main:          epsdis= ', epsdis
              write(*,*) 'main:          lambda= ', lambda
              write(*,*) 'main:     trace_level= ', trace_level
            end if
            if (trace_level.gt.4) call params_dump()
            if (trace_level.gt.0) then
              write(*,*) 'main: **** END: Parameters ****'
            end if
         end if

       return
       end     ! END: get_params()


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  writestate:                                                               cc
cc             Save current state of code/processor to file for               cc
cc             checkpointing.                                                 cc
cc         NB: Added control to have processors take turns dumping state:     cc
cc             chkpt_control == 0 all dump at same time [default]             cc
cc                           == 1 take turns                                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine writestate()
       implicit     none
       include     'glob.inc'
       include     'param.inc'
       include     'mpif.h'
       include     'mpi_stuff.inc'
       !  processer whose turn it is to dump state (if taking turns):
       !    (also used a dummy variable for the MPI calls)
       integer     whoseturn
       !
       logical     ltrace
       parameter ( ltrace  = .true. )
       logical     ltrace2
       parameter ( ltrace2 = .false. )

       call level_test_local()

       if (ltrace2) then
          write(*,99) myid, 'writestate: chkpt_control=',chkpt_control
          write(*,98) myid, 'writestate: Saving state to:',
     *             chkpt_fname(chkpt_fnum)(1:20)
       end if

       ! Just go in numerical order and start with 0:
       whoseturn = 0

       if (chkpt_control .eq. 0 .or. whoseturn .eq. myid) then
          !
 100      if(ltrace2) write(*,99) myid, 'writestate: Saving: parameters'
          if(ltrace)write(*,98) myid, 'writestate: Saving state to:',
     *             chkpt_fname(chkpt_fnum)(1:20)
          call params_writestate(chkpt_fname(chkpt_fnum))
          if (ltrace2) write(*,99) myid, 'writestate: Saving: mem'
          call mem_writestate( chkpt_fname(chkpt_fnum) )
          if (ltrace2) write(*,99) myid, 'writestate: Closing SDFs'
#ifdef AIX
          call gft_close(trim(chkpt_fname(chkpt_fnum))/ /CHAR(0))
#else
          call gft_close(chkpt_fname(chkpt_fnum))
#endif
          !
          if (chkpt_control.gt.0 .and. myid+1.lt.numprocs) then
             ! Let next processor know to go ahead:
             if(ltrace2)write(*,99)myid,'writestate: Go ahead:',myid+1
             call MPI_Send(whoseturn, 0,
     *                     MPI_INTEGER, myid+1, TAG_CHKPT,
     *                     MPI_COMM_WORLD, ierr)
          end if
          !
       else 
          !
          ! Not our turn, wait for signal from our predecessor:
          !
          if(ltrace2)write(*,99)myid,'writestate: Waiting...:',myid-1
          call MPI_Recv(whoseturn,
     *                 0, MPI_INTEGER,
     *                 myid-1, TAG_CHKPT, MPI_COMM_WORLD, status,ierr)
          goto 100
       end if

       if (ltrace2) then
          write(*,99) myid, 'writestate: Done saving state'
       end if

       !
       ! Switch filename so we save to a different file next time:
       !
       chkpt_fnum = 3 - chkpt_fnum

  98   format('[',I4,'] ',A,A)
  99   format('[',I4,'] ',A,3I5)
       return
       end     ! END: writestate

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  readstate:                                                                cc
cc             Read current state of code/processor to file for               cc
cc             checkpointing.                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine readstate()
       implicit     none
       include     'glob.inc'
       include     'grid.inc'
       include     'param.inc'
       include     'mpif.h'
       include     'mpi_stuff.inc'
       ! Have the number input by user select which fileset to read:
       integer     num, i
       logical     stateread
       logical     grid_alive
       external    grid_alive
       !
       logical     ltrace
       parameter ( ltrace = .true. )
       logical     ltrace2
       parameter ( ltrace2 = .false. )

       stateread = .true.

       !
       ! The first decimal digit encodes which filename to use
       !   (the second encodes whether to read parameters from
       !    the chkpt files or from the par file)
       !
       num       = mod(chkpt_readstate,10)
       if (ltrace) then
!         write(*,97)myid,'readstate: chkpt_readstate: num:',num
          write(*,98)myid,'readstate: Reading state from:',
     *         chkpt_fname(num)(1:20)
!         write(*,*)myid,'readstate: name: ',
!    *         chkpt_fname(num)(1:20)
       end if

       !call gft_set_single( chkpt_fname(chkpt_fnum) )

       if (ltrace2) write(*,99) myid, 'readstate: Reading: parameters'
       !call params_readstate(chkpt_fname(chkpt_fnum))
       call params_readstate(chkpt_fname(num))
       if (ltrace2.and.myid.eq.master) call params_dump()
       if (ltrace2) write(*,99) myid, 'readstate: Reading: mem'
       call mem_readstate( chkpt_fname(num) )
       if (ltrace2.and.myid.eq.master) call proc_dump_wkld()

       !call gft_set_multi()

       stateread = .false.

       !if (.not. stateread) then
          !write(*,99) myid, 'readstate: Failed to read state.'
          !call my_exit('readstate: Failed to read state')
       !end if

       if (num.eq.1) then
          if (ltrace2) then 
             write(*,99)myid,'readstate: ',chkpt_fnum
             write(*,98)myid,'readstate: Because we read: ',
     *                                       chkpt_fname(chkpt_fnum)
             write(*,99)myid,'readstate: we will set things so that'
          end if
          chkpt_fnum = 2
          if (ltrace2) then 
             write(*,99)myid,'readstate: the first checkpoint will be'
             write(*,98)myid,'readstate: with the files: ',
     *                                       chkpt_fname(chkpt_fnum)
             write(*,99)myid,'readstate: ',chkpt_fnum
          end if
       end if

       call level_test_local()

       if (ltrace2) then
          call mem_stat()
          !do i = 1, max_numgrids
             !if (grid_alive(i)) call grid_dump_info(i)
          !end do
          write(*,99) myid, 'readstate: Done reading state'
       end if

  97   format('[',I4,'] ',A,I5)
  98   format('[',I4,'] ',A,A)
  99   format('[',I4,'] ',A,3I5)
       return
       end     ! END: readstate

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  surface_out:                                                              cc
cc             Output surface functions and/or integrals.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine surface_out(time)
       implicit     none
       real(kind=8) time
       include     'surfaces.inc'
       include     'mem.inc'
       !
       integer     i, gft_rc, index, j, k
       integer     ntheta, nphi
       integer     gft_out_bbox, mystringlength
       external    gft_out_bbox, mystringlength
       real(kind=8) dtheta, theta, integral1, integral2
       real(kind=8) dphi,   phi, sinphidtdp
       real(kind=8) bbox(2)
       data         bbox / 0.d0, 0.d0 /
       character(20) fname
       logical     ltrace
       parameter ( ltrace  = .false. )
       logical     ltrace2
       parameter ( ltrace2 = .false. )

       if (ltrace) write(*,*) 'surface_out: Outputting surfaces'

#ifdef AIX

       !
       ! Surface A:
       !
       do i = 1, num_afuncs
          if (asf_out(i) .eq. 0 .or. asf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out fct:',i,asf_name(i)

             gft_rc = gft_out_bbox(  trim(asf_name(i))/ /CHAR(0),
     *                           %VAL(time),   asf_shape,
     *                           %VAL(2), asf_bbox, q(asf_pointer(i)) )
          end if

          if (asf_out(i) .eq.1 .or. asf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out Int:',i,asf_name(i)
          end if
       end do

       !
       ! Surface B:
       !
       do i = 1, num_bfuncs
          if (bsf_out(i) .eq. 0 .or. bsf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out fct:',i,bsf_name(i)
             gft_rc = gft_out_bbox(  trim(bsf_name(i))/ /CHAR(0),
     *                            %VAL(time),   bsf_shape,
     *                            %VAL(2), bsf_bbox, q(bsf_pointer(i)) )
          end if
          if (bsf_out(i) .eq.1 .or. bsf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out Int:',i,bsf_name(i)
          end if
       end do

       !
       ! Surface C:
       !
       do i = 1, num_cfuncs
          if (csf_out(i) .eq. 0 .or. csf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out fct:',i,csf_name(i)
             gft_rc = gft_out_bbox(  trim(csf_name(i))/ /CHAR(0),
     *                            %VAL(time),   csf_shape,
     *                            %VAL(2), csf_bbox, q(csf_pointer(i)) )
          end if
          if (csf_out(i) .eq.1 .or. csf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out Int:',i,csf_name(i)
          end if
       end do

#else
       !
       ! Surface A:
       !
       do i = 1, num_afuncs
          if (asf_out(i) .eq. 0 .or. asf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out fct:',i,asf_name(i)
             gft_rc = gft_out_bbox(  asf_name(i), time,   asf_shape,
     *                               2, asf_bbox, q(asf_pointer(i)) )
          end if
          if (asf_out(i) .eq.1 .or. asf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out Int:',i,asf_name(i)
             !
             ! Compute Integral;
             !
             integral1 = 0.d0
             integral2 = 0.d0
             ntheta    = asf_shape(1)
             nphi      = asf_shape(2)
             !
             dtheta = ( asf_bbox(2) - asf_bbox(1) ) / (ntheta-1)
             dphi   = ( asf_bbox(4) - asf_bbox(3) ) / (nphi  -1)
             !
             do k = 1, nphi
                phi        = (k-1)* dphi
                sinphidtdp = sin(phi) * dtheta * dphi
                !
                do j = 1, ntheta
                   theta = (j-1)* dtheta
                   index = (k-1)*ntheta + (j-1)
                   !
                   integral1 =  integral1 + sinphidtdp
                   integral2 =  integral2
     *                             + sinphidtdp *q(asf_pointer(i)+index)
                end do
             end do
             !
             fname = asf_name(i)
             k     = mystringlength(fname)
             fname(k+1:) = '_integral'
             gft_rc = gft_out_bbox( fname, time, 1, 1, bbox, integral2 )
             !
             if (ltrace2) then
                write(*,*)'surface_out: fname     = ',fname
                write(*,*)'surface_out: ntheta    = ',ntheta
                write(*,*)'surface_out: nphi      = ',nphi
                write(*,*)'surface_out: dtheta    = ',dtheta
                write(*,*)'surface_out: dphi      = ',dphi  
                write(*,*)'surface_out: integral1 = ',integral1
                write(*,*)'surface_out: integral2 = ',integral2
             end if
          end if
       end do

       !
       ! Surface B:
       !
       do i = 1, num_bfuncs
          if (bsf_out(i) .eq. 0 .or. bsf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out fct:',i,bsf_name(i)
             gft_rc = gft_out_bbox(  bsf_name(i), time,   bsf_shape,
     *                               2, bsf_bbox, q(bsf_pointer(i)) )
          end if
          if (bsf_out(i) .eq.1 .or. bsf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out Int:',i,bsf_name(i)
             !
             ! Compute Integral;
             !
             integral1 = 0.d0
             integral2 = 0.d0
             ntheta    = bsf_shape(1)
             nphi      = bsf_shape(2)
             !
             dtheta = ( bsf_bbox(2) - bsf_bbox(1) ) / (ntheta-1)
             dphi   = ( bsf_bbox(4) - bsf_bbox(3) ) / (nphi  -1)
             !
             do k = 1, nphi
                phi        = (k-1)* dphi
                sinphidtdp = sin(phi) * dtheta * dphi
                !
                do j = 1, ntheta
                   theta = (j-1)* dtheta
                   index = (k-1)*ntheta + (j-1)
                   !
                   integral1 =  integral1 + sinphidtdp
                   integral2 =  integral2
     *                             + sinphidtdp *q(bsf_pointer(i)+index)
                end do
             end do
             !
             fname = bsf_name(i)
             k     = mystringlength(fname)
             fname(k+1:) = '_integral'
             gft_rc = gft_out_bbox( fname, time, 1, 1, bbox, integral2 )
             !
             if (ltrace2) then
                write(*,*)'surface_out: fname     = ',fname
                write(*,*)'surface_out: ntheta    = ',ntheta
                write(*,*)'surface_out: nphi      = ',nphi
                write(*,*)'surface_out: dtheta    = ',dtheta
                write(*,*)'surface_out: dphi      = ',dphi  
                write(*,*)'surface_out: integral1 = ',integral1
                write(*,*)'surface_out: integral2 = ',integral2
             end if
          end if
       end do

       !
       ! Surface C:
       !
       do i = 1, num_cfuncs
          if (csf_out(i) .eq. 0 .or. csf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out fct:',i,csf_name(i)
             gft_rc = gft_out_bbox(  csf_name(i), time,   csf_shape,
     *                               2, csf_bbox, q(csf_pointer(i)) )
          end if
          if (csf_out(i) .eq.1 .or. csf_out(i) .eq.2) then
             if (ltrace2)write(*,*)'surface_out: Out Int:',i,csf_name(i)
             !
             ! Compute Integral;
             !
             integral1 = 0.d0
             integral2 = 0.d0
             ntheta    = csf_shape(1)
             nphi      = csf_shape(2)
             !
             dtheta = ( csf_bbox(2) - csf_bbox(1) ) / (ntheta-1)
             dphi   = ( csf_bbox(4) - csf_bbox(3) ) / (nphi  -1)
             !
             do k = 1, nphi
                phi        = (k-1)* dphi
                sinphidtdp = sin(phi) * dtheta * dphi
                !
                do j = 1, ntheta
                   theta = (j-1)* dtheta
                   index = (k-1)*ntheta + (j-1)
                   !
                   integral1 =  integral1 + sinphidtdp
                   integral2 =  integral2
     *                             + sinphidtdp *q(csf_pointer(i)+index)
                end do
             end do
             !
             fname = csf_name(i)
             k     = mystringlength(fname)
             fname(k+1:) = '_integral'
             gft_rc = gft_out_bbox( fname, time, 1, 1, bbox, integral2 )
             !
             if (ltrace2) then
                write(*,*)'surface_out: fname     = ',fname
                write(*,*)'surface_out: ntheta    = ',ntheta
                write(*,*)'surface_out: nphi      = ',nphi
                write(*,*)'surface_out: dtheta    = ',dtheta
                write(*,*)'surface_out: dphi      = ',dphi  
                write(*,*)'surface_out: integral1 = ',integral1
                write(*,*)'surface_out: integral2 = ',integral2
             end if
          end if
       end do

#endif


       return
       end     ! END: surface_out

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  surface_compute:                                                          cc
cc             Compute surface functions and/or integrals.                    cc
cc             All processors call this routine to compute the surface funcs. cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine surface_compute()
       implicit     none
       include     'param.inc'
       include     'surfaces.inc'
       include     'mem.inc'
       include     'mpif.h'
       include     'mpi_stuff.inc'
       include     'largesmall.inc'
       !
       integer        i, j, ntheta, nphi, index
       real(kind=8)   x, y, z
       real(kind=8)    theta,  phi
       real(kind=8)   dtheta, dphi
       real(kind=8)   value
       integer        li, gi, own, length
       integer        l, k

       !
       logical     ltrace
       parameter ( ltrace  = .false. )
       logical     ltrace2
       parameter ( ltrace2 = .false. )
       logical     ltraceB
       parameter ( ltraceB = .false. )

       if (ltrace) write(*,99) myid, 'surface_compute: Computing...'

       !....................................
       ! Surface A:
       !....................................
       ntheta = asf_shape(1)
       nphi   = asf_shape(2)
       length = ntheta * nphi
       !
       ! Initialize functions to large number
       ! so that we can tell where we have actually
       ! computed values:
       !
       do k = 1, num_afuncs-1
          call load_scal1D(q(asf_pointer(k)),LARGENUMBER, length)
       end do
       !
       dtheta = ( asf_bbox(2) - asf_bbox(1) ) / (ntheta-1)
       dphi   = ( asf_bbox(4) - asf_bbox(3) ) / (nphi  -1)
       !
       do j = 1, nphi
          phi = (j-1)* dphi
          !
          do i = 1, ntheta
             theta = (i-1)* dtheta
             !
             ! Get Cartesian coordinates:
             !
             x = asf_rconst * cos(theta) * sin(phi)
             y = asf_rconst * sin(theta) * sin(phi)
             z = asf_rconst *              cos(phi)
             if(ltrace2)write(*,98)myid,' Pt: ',theta,phi,x,y,z
             !
             ! If point contained at highest resolution on
             ! this processor, then compute function here:
             !
             call level_finest_containing(x,y,z,li,gi,own)
             if (own .eq. myid) then
                if(ltrace2)write(*,98)myid,' We own Pt:',theta,phi,x,y,z
                !
                ! compute functions at this point:
                !
                !
                if (li.lt.asf_level.and. .false.) then
                   write(*,98)myid,'surface_compute:Warning Coarse Lev',
     .                                theta,phi,x,y,z
                   !write(*,98)myid,'surface_compute: occurs on coarse'
                   !write(*,98)myid,'surface_compute: level which may'
                   !write(*,98)myid,'surface_compute: not exist at this'
                   !write(*,98)myid,'surface_compute: time and therefore'
                   !write(*,98)myid,'surface_compute: requires interpola'
                end if
                call grid_asurface_comp(gi,x,y,z,i,j)
                !
             end if
             !
          end do
       end do
       !
       ! At this point, the functions are computed, but
       ! the values are distributed among the different processors.
       ! Collect them all on the master process:
       !
       if (myid .eq. master) then
          !
          do l = 1, numprocs-1
             if(ltrace)write(*,99)myid,'  Recv data from: ',l
             do k = 1, num_afuncs-1
                if(ltrace)write(*,99)myid,'     field: ',k
                call MPI_Recv(q(asf_pointer(num_afuncs)),
     *                 length, MPI_DOUBLE_PRECISION,
     *                 l, TAG_SURFACE, MPI_COMM_WORLD, status,ierr)
                !
                do j = 1, nphi
                do i = 1, ntheta
                   index = (j-1)*ntheta + (i-1)
                   value = q(asf_pointer(num_afuncs)+index)
                   if (value .lt. LARGE_M_SMALL) then
                      q(asf_pointer(k)+index) = value
                   end if
                end do
                end do
                !
             end do
          end do
          k = 1
          if (ltrace) 
     .          call field_dump_info2D(q(asf_pointer(k)),ntheta,nphi,
     .                                   asf_name(k))
       else
          !
          if(ltrace)write(*,98)myid,'  Send data '
          do k = 1, num_afuncs-1
             if(ltrace)write(*,99)myid,'     field: ',k
             call MPI_Send(q(asf_pointer(k)), length,
     *                     MPI_DOUBLE_PRECISION, master, TAG_SURFACE,
     *                     MPI_COMM_WORLD, ierr)
          end do
       end if

       if (ltrace) write(*,99) myid, 'surface_compute: Done w/ A.'

       !....................................
       ! Surface B:
       !....................................
       ntheta = bsf_shape(1)
       nphi   = bsf_shape(2)
       length = ntheta * nphi
       !
       ! Initialize functions to large number
       ! so that we can tell where we have actually
       ! computed values:
       !
       do k = 1, num_bfuncs-1
          call load_scal1D(q(bsf_pointer(k)),LARGENUMBER, length)
       end do
       !
       dtheta = ( bsf_bbox(2) - bsf_bbox(1) ) / (ntheta-1)
       dphi   = ( bsf_bbox(4) - bsf_bbox(3) ) / (nphi  -1)
       !
       do j = 1, nphi
          phi = (j-1)* dphi
          !
          do i = 1, ntheta
             theta = (i-1)* dtheta
             !
             ! Get Cartesian coordinates:
             !
             x = bsf_rconst * cos(theta) * sin(phi)
             y = bsf_rconst * sin(theta) * sin(phi)
             z = bsf_rconst *              cos(phi)
             if(ltrace2)write(*,98)myid,' Pt: ',theta,phi,x,y,z
             !
             ! If point contained at highest resolution on
             ! this processor, then compute function here:
             !
             call level_finest_containing(x,y,z,li,gi,own)
             if (own .eq. myid) then
                if(ltraceB)write(*,97)myid,' We own B Pt:',
     *                                     theta,phi,x,y,z,gi
                !
                ! compute functions at this point:
                !
                !
                call grid_bsurface_comp(gi,x,y,z,i,j)
                !
             end if
             !
          end do
       end do

       !
       ! At this point, the functions are computed, but
       ! the values are distributed among the different processors.
       ! Collect them all on the master process:
       !
       if (myid .eq. master) then
          !
          do l = 1, numprocs-1
             if(ltrace)write(*,99)myid,'  Recv data from: ',l
             do k = 1, num_bfuncs-1
                if(ltrace)write(*,99)myid,'     field: ',k
                call MPI_Recv(q(bsf_pointer(num_bfuncs)),
     *                 length, MPI_DOUBLE_PRECISION,
     *                 l, TAG_SURFACE, MPI_COMM_WORLD, status,ierr)
                !
                do j = 1, nphi
                do i = 1, ntheta
                   index = (j-1)*ntheta + (i-1)
                   value = q(bsf_pointer(num_bfuncs)+index)
                   if (value .lt. LARGE_M_SMALL) then
                 !write(*,*) 'on B, using value from other rpoc: ',value
                      q(bsf_pointer(k)+index) = value
                   end if
                end do
                end do
                !
             end do
          end do
       else
          !
          if(ltrace)write(*,98)myid,'  Send data '
          do k = 1, num_bfuncs-1
             if(ltrace)write(*,99)myid,'     field: ',k
             call MPI_Send(q(bsf_pointer(k)), length,
     *                     MPI_DOUBLE_PRECISION, master, TAG_SURFACE,
     *                     MPI_COMM_WORLD, ierr)
          end do
       end if
        
c     call field_out2d(q(bsf_pointer(1)),0.d0, 'bsf_1',
c    *             bsf_bbox(1),bsf_bbox(2),
c    *             bsf_bbox(3),bsf_bbox(4),
c    *             bsf_shape(1),bsf_shape(2),myid)

       if (ltrace) write(*,99) myid, 'surface_compute: Done w/ B.'

       !....................................
       ! Surface C:
       !....................................
       ntheta = csf_shape(1)
       nphi   = csf_shape(2)
       length = ntheta * nphi
       !
       ! Initialize functions to large number
       ! so that we can tell where we have actually
       ! computed values:
       !
       do k = 1, num_cfuncs-1
          call load_scal1D(q(csf_pointer(k)),LARGENUMBER, length)
       end do
       !
       dtheta = ( csf_bbox(2) - csf_bbox(1) ) / (ntheta-1)
       dphi   = ( csf_bbox(4) - csf_bbox(3) ) / (nphi  -1)
       !
       do j = 1, nphi
          phi = (j-1)* dphi
          !
          do i = 1, ntheta
             theta = (i-1)* dtheta
             !
             ! Get Cartesian coordinates:
             !
             x = csf_rconst * cos(theta) * sin(phi)
             y = csf_rconst * sin(theta) * sin(phi)
             z = csf_rconst *              cos(phi)
             if(ltrace2)write(*,98)myid,' Pt: ',theta,phi,x,y,z
             !
             ! If point contained at highest resolution on
             ! this processor, then compute function here:
             !
             call level_finest_containing(x,y,z,li,gi,own)
             if (own .eq. myid) then
                if(ltrace2)write(*,98)myid,' We own Pt:',theta,phi,x,y,z
                !
                ! compute functions at this point:
                !
                !
                call grid_csurface_comp(gi,x,y,z,i,j)
                !
             end if
             !
          end do
       end do

       !
       ! At this point, the functions are computed, but
       ! the values are distributed among the different processors.
       ! Collect them all on the master process:
       !
       if (myid .eq. master) then
          !
          do l = 1, numprocs-1
             if(ltrace)write(*,99)myid,'  Recv data from: ',l
             do k = 1, num_cfuncs-1
                if(ltrace)write(*,99)myid,'     field: ',k
                call MPI_Recv(q(csf_pointer(num_cfuncs)),
     *                 length, MPI_DOUBLE_PRECISION,
     *                 l, TAG_SURFACE, MPI_COMM_WORLD, status,ierr)
                !
                do j = 1, nphi
                do i = 1, ntheta
                   index = (j-1)*ntheta + (i-1)
                   value = q(csf_pointer(num_cfuncs)+index)
                   if (value .lt. LARGE_M_SMALL) then
                      q(csf_pointer(k)+index) = value
                   end if
                end do
                end do
                !
             end do
          end do
       else
          !
          if(ltrace)write(*,98)myid,'  Send data '
          do k = 1, num_cfuncs-1
             if(ltrace)write(*,99)myid,'     field: ',k
             call MPI_Send(q(csf_pointer(k)), length,
     *                     MPI_DOUBLE_PRECISION, master, TAG_SURFACE,
     *                     MPI_COMM_WORLD, ierr)
          end do
       end if

       if (ltrace) write(*,99) myid, 'surface_compute: Done w/ C.'
   
 97    format('[',I3,'] ',A,5F12.6,I4)
 98    format('[',I3,'] ',A,5F12.6)
 99    format('[',I3,'] ',A,2I5)


       return
       end     ! END: surface_compute


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  surface_init:                                                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine surface_init()
       implicit     none
       include     'glob.inc'
       include     'param.inc'
       include     'surfaces.inc'
       include     'mem.inc'
       include     'mask.inc'
       include     'mpif.h'
       include     'mpi_stuff.inc'
       !
       logical      secondtime
       save         secondtime
       integer      i, tmpi
       integer      mem_alloc
       external     mem_alloc
       !
       logical     ltrace
       parameter ( ltrace  = .false. )

      ! Use a couple tests so as not to init twice
      ! Only the non-master processes call this twice.
      if (.not.secondtime) then
         secondtime = .true.
      else if (asf_shape(1).eq.asf_ntheta) then
         secondtime = .true.
         write(*,99)myid,'surface_init: Second time through, skipping.'
         return
      end if
      !
      ! Surface A:
      !
      asf_bbox(1)  = 0.d0
      asf_bbox(2)  = 2. * cPi
      asf_bbox(3)  = 0.d0
      asf_bbox(4)  =      cPi
      !
      asf_shape(1) = asf_ntheta
      asf_shape(2) = asf_nphi
      !
      tmpi         = asf_ntheta * asf_nphi
      do i = 1, num_afuncs
         asf_pointer(i) = mem_alloc(tmpi)
      end do
      !
      if (ltrace) then
         write(*,99) myid, 'surface_init: tmpi       =',tmpi
         write(*,99) myid, 'surface_init: asf_ntheta =',asf_ntheta
         write(*,99) myid, 'surface_init: asf_nphi   =',asf_nphi
         write(*,98) myid, 'surface_init: asf_rconst =',asf_rconst
         write(*,99) myid, 'surface_init: asf_period =',asf_period
      end if

      !
      ! Surface B:
      !
      bsf_bbox(1)  = 0.d0
      bsf_bbox(2)  = 2. * cPi
      bsf_bbox(3)  = 0.d0
      bsf_bbox(4)  =      cPi
      !
      bsf_shape(1) = bsf_ntheta
      bsf_shape(2) = bsf_nphi
      !
      tmpi         = bsf_ntheta * bsf_nphi
      do i = 1, num_bfuncs
         bsf_pointer(i) = mem_alloc(tmpi)
      end do
      !
      if (ltrace) then
         write(*,99) myid, 'surface_init: tmpi       =',tmpi
         write(*,99) myid, 'surface_init: bsf_ntheta =',bsf_ntheta
         write(*,99) myid, 'surface_init: bsf_nphi   =',bsf_nphi
         write(*,98) myid, 'surface_init: bsf_rconst =',bsf_rconst
         write(*,99) myid, 'surface_init: bsf_period =',bsf_period
      end if

      !
      ! Surface C:
      !
      csf_bbox(1)  = 0.d0
      csf_bbox(2)  = 2. * cPi
      csf_bbox(3)  = 0.d0
      csf_bbox(4)  =      cPi
      !
      csf_shape(1) = csf_ntheta
      csf_shape(2) = csf_nphi
      !
      tmpi         = csf_ntheta * csf_nphi
      do i = 1, num_cfuncs
         csf_pointer(i) = mem_alloc(tmpi)
      end do

      !
      ! Need to go ahead an allocate even if not needed,
      ! because down the road in a future checkpoint,
      ! someone could turn on horizon finding, so we always
      ! just allocate this stuff:
      !
      ! Apparent Horizon stuff:
      !
      ! *NOT* Using as convention that in
      !     http://mathworld.wolfram.com/SphericalCoordinates.html
      ! Instead, using conventionso fo frans: gr-qc/0407110
      !if (abs(findhorizon).eq.2 .and. use_mask.gt.0) then
      if (.true.) then
         if(ltrace)write(*,99)myid,'surface_init: Initing horizon'
         horizon_shape(1) = horizon_ntheta    ! 0-->pi
         horizon_shape(2) = horizon_nphi      ! 0-->2pi
         !
         tmpi = horizon_ntheta * horizon_nphi
         do i = 1, max_num_masks
            horizon_theta(i) = mem_alloc(tmpi)
            horizon_R(i)     = mem_alloc(tmpi)
            !if (HORIZON_USEICN) then
            if (.true.) then
               ! If using iterative CN for flow evolution need
               ! storage for np1 iterative:
               horizon_theta_np1(i) = mem_alloc(tmpi)
               horizon_R_np1(i)     = mem_alloc(tmpi)
            end if
         end do
         !
      end if

      !
      if (ltrace) then
         write(*,99) myid, 'surface_init: tmpi       =',tmpi
         write(*,99) myid, 'surface_init: csf_ntheta =',csf_ntheta
         write(*,99) myid, 'surface_init: csf_nphi   =',csf_nphi
         write(*,98) myid, 'surface_init: csf_rconst =',csf_rconst
         write(*,99) myid, 'surface_init: csf_period =',csf_period
      end if

 98    format('[',I3,'] ',A,5F12.6)
 99    format('[',I3,'] ',A,2I5)

       return
       end     ! END: surface_init

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  masks_return_num:                                                         cc
cc                     returns the current value of num_masks                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       integer function masks_return_num()
       implicit     none
       include     'mask.inc'
       !
       logical     ltrace
       parameter ( ltrace  = .false. )

       masks_return_num = num_masks

       return
       end     ! END: masks_return_num
