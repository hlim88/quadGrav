cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  ell_solve:                                                                cc
cc                Top level routine to do an elliptic solve.                  cc
cc                Solve proceeds from finest level to coarse level "lev"      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ell_solve( lev )
      implicit none
      integer     lev
      include     'grid_methods.inc'
      include     'grid.inc'
      include     'glob.inc'
      include     'action.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'output.inc'
      include     'param.inc'
      real(kind=8) residnorm, glob_residnorm, oldnorm, percentchg
      real(kind=8) oresidnorm
      integer      i, finest, l, cycle
      logical      converging
      ! coarsest grid nx/y/z
      integer      nxc, nyc, nzc       
      ! factors of 2 coarsening for each dim
      integer      nfact_x, nfact_y, nfact_z 
      ! return codes from factoring
      integer      xrc, yrc, zrc
      !
      !
      ! Else solve by multigrid:
      !
      logical      solve_by_relaxation
      parameter (  solve_by_relaxation  = .false. )
      !
      ! How often should we update the residual while doing relaxation:
      !
      integer         relax_outperiod
      parameter     ( relax_outperiod   = 200 )
      ! Continue sweeping until residual is how small?
      real(kind=8)    relax_resid
      parameter     ( relax_resid       = 1.d-8 )
      ! Continue vcycling until resid at least below:
      real(kind=8)    minvcycres
      parameter     ( minvcycres        = 1.d-3 )
      ! Give up after this many vcycles:
      integer         maxnumvcyc
      parameter     ( maxnumvcyc        = 20    )
      !parameter     ( maxnumvcyc        = 5    )
      ! Number of vcycles for coarser levels in full multigrid:
      !    (if set to zero, then no full multigrid):
      integer         FMGMAX    
      parameter     ( FMGMAX            = 1    )
      integer         FMGlev
      ! How many times to try smoothing if nonconverging?
      integer         maxsmoothattempts, smoothattempts
      parameter     ( maxsmoothattempts = 0     )
      !
      ! Multigrid constants:
      integer         lmax,          minnxc,           maxnxc
      parameter     ( lmax = 12,     minnxc = 5,       maxnxc = 9  )
      !
      logical      do_nothing
      parameter (  do_nothing           = .false. )
      !
      integer         checkfacts
      external        checkfacts
      !
      logical      ltrace
      parameter (  ltrace  = .true. )
      logical      ltrace2
      parameter (  ltrace2 = .false. )
      logical      vtrace
      parameter (  vtrace = .true. )

      if (do_nothing) return

      if (.false. .and. lev.ne.0) then
         write(*,99)myid,'ell_solve: Not solving:',finest,lev
         return
      end if

      finest = level_return_finest(Levelp,maxlev)

      if (.false. .and.myid.eq.master.or.ltrace2) then
 91      format('[',I4,'] ',A,I3,t50,1p,e9.2)
 98      format('[',I4,'] ',A,3F15.9)
 99      format('[',I4,'] ',A,3I7)
         write(*,99)myid,'ell_solve: On levels : ',finest,lev
      end if

      !
      ! Compute the RHS for the elliptic problem:
      !
      if (ltrace2.and.myid.eq.master.or.ltrace2)
     *     write(*,99)myid, 'ell_solve: Compute RHS:     ',finest

      if (solve_by_relaxation) then
         !
         ! Only relax finest existing level:
         !
         if(ltrace)write(*,99)myid,'ell_solve: Solve by relax: ',finest
         i = 0
         call level_ellrhs_local(finest)
         residnorm = level_return_resnrm(finest)
         if (myid.eq.master.or.ltrace2)
     *   write(*,91)myid,'ell_solve: i, residnorm= ',i,residnorm
         do i = 1, maxsweeps
            !
            call level_relax_local(finest)
            call level_syncbndi_local(finest)
            !
            if (mod(i,relax_outperiod).eq.0) then
               residnorm = level_return_resnrm(finest)
               if (myid.eq.master)
     *         write(*,91)myid,'ell_solve: i, residnorm= ',i,residnorm
               call MPI_AllReduce(residnorm,glob_residnorm,1,
     *              MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
               if (glob_residnorm .lt. relax_resid) goto 66
            end if
         end do
 66      if(ltrace)write(*,99)myid,'ell_solve: Relaxation completed',i-1
      else
         !
         ! Working only with coarsest grid
         !
         if (ltrace2 .and. myid.eq.master) then
            write(*,*) '******** parameter dump ************'
            write(*,*) '    number of levels            l =', l
            write(*,*) '    number of v-cycles    nvcycle =', nvcycle
            write(*,*) '    number of pre-sweeps   preswp =', preswp
            write(*,*) '    number of post-sweeps  pstswp =', pstswp
            write(*,*) '************ end dump **************'
         end if
         !-------------------------------------
         ! Factorize Grid Handed-in into Sub-Grids
         !        minnxc,  maxnyc  --- minimum / maximum size of coarse grid
         !        nx0,     ny0     --- fine grid resolution as passed into routine
         !        nfact_x, nfact_y --- number of factors of two in (nx0,ny0)
         !-------------------------------------
         if (.false.) then
            ! such a check would only apply if weren't strictly
            ! solving things on the aMR hierarchy itself:
         xrc = checkfacts(nx0 - 1,maxnxc,minnxc,nxc,nfact_x)
         yrc = checkfacts(ny0 - 1,maxnxc,minnxc,nyc,nfact_y)
         zrc = checkfacts(nz0 - 1,maxnxc,minnxc,nzc,nfact_z)
         if( xrc .lt. 0  .or.  yrc .lt. 0 .or. zrc .lt. 0) then
            write(*,*) 'solve_fas: Grid of size: ', nx0,ny0,nz0
            write(*,*) nx0-1,' factors to ',nfact_x,' factors of 2'
            write(*,*) ' with remaining: ', nxc
            write(*,*) ny0-1,' factors to ',nfact_y,' factors of 2'
            write(*,*) ' with remaining: ', nyc
            write(*,*) nz0-1,' factors to ',nfact_z,' factors of 2'
            write(*,*) ' with remaining: ', nzc
            write(*,*) 'solve_fas: not factorizable down to'
            write(*,*) 'solve_fas: maximum coarse size: ',maxnxc
            write(*,*) 'solve_fas: Quitting...'
            call my_exit('Cannot factorize coarse grid sufficiently')
         end if
         !
         !     Set number of levels and coarse grid res:
         !
         if (nfact_x. lt. nfact_y .and. nfact_x .lt. nfact_z) then
               l   = nfact_x + 1
               nxc = 1 + nxc
               nyc = 1 + nyc*(2**(nfact_y-nfact_x))
               nzc = 1 + nzc*(2**(nfact_z-nfact_x))
         else if (nfact_y .le. nfact_x .and. nfact_y .lt. nfact_z) then
               l   = nfact_y + 1
               nxc = 1 + nxc*(2**(nfact_x-nfact_y))
               nyc = 1 + nyc
               nzc = 1 + nzc*(2**(nfact_z-nfact_y))
         else 
               l   = nfact_z + 1
               nxc = 1 + nxc*(2**(nfact_x-nfact_z))
               nyc = 1 + nyc*(2**(nfact_y-nfact_z))
               nzc = 1 + nzc
         end if
         !
         !     Check number of levels:
         !
         if ( l .gt. lmax ) then
             write(*,*) 'ell_solve: too many levels ', l,lmax
             call my_exit('Too many levels')
         end if
         !
         if (ltrace .and. .false.) then
           write(*,*)'ell_solve: l:           ',l
           write(*,*)'ell_solve: nx/y/z0:     ',nx0,ny0,nz0
           write(*,*)'ell_solve: nx/y/zc:     ',nxc,nyc,nzc
           write(*,*)'ell_solve: nfact_x/y/z: ',nfact_x,nfact_y,nfact_z
         end if
            ! cut the above section out since we are just using the amr hierarchy
         end if
         !
         !-------------------------------------
         ! Create grid structure
         !-------------------------------------

         !-------------------------------------
         ! Use FMG
         !      1) Solve coarse grid problem
         !      2) Prolong solution to finer level
         !      3) Repeat
         !-------------------------------------
         do FMGlev = 0, finest
            oldnorm        = 1.d99
            glob_residnorm = 2.*oldnorm
            converging     = .true.
            smoothattempts = 0
            cycle          = 0
            !
            if (FMGlev.gt.0) then
               if (ltrace) write(*,*) 'ellsolve: prolonging:'
               call level_mgprolongto(FMGlev)
            end if
            if (ltrace) write(*,*) 'ellsolve: Level: ',FMGlev
            call level_ellrhs_local(FMGlev)
            if (vtrace) call level_mg_vtrace(FMGlev,0)
            !
 23         if ( converging .and.
     *           cycle .lt. maxnumvcyc .and.
     *           .not.(FMGlev.ne.finest.and.cycle.gt.FMGMAX-1) .and.
     *             (cycle.lt.nvcycle .or. 
     *              glob_residnorm.gt.minvcycres) ) then
               !
               ! Keep Vcycling if
               !    --- it is driving the residual down
               !    --- we have not reached the maxnumvcyc
               !    --- we are on finest level (otherwise just do one vcycle)
               !    --- either:
               !            ...the user has asked for more vcycles
               !            ...the residual is above some threshold
               !
               cycle      = cycle + 1
               oldnorm    = glob_residnorm
               !
               if(myid.eq.master.and.ltrace2)write(*,91)
     *                  myid,'ell_solve: Vcycling:',cycle,glob_residnorm
               !
               call vcycle(FMGlev,lev,glob_residnorm)
               if(ltrace2)write(*,98)myid,'ell_solve:GR:',glob_residnorm
               if (vtrace) call level_mg_vtrace(FMGlev,cycle)
               !
               percentchg = 100.*(oldnorm-abs(glob_residnorm)) / oldnorm
               converging = percentchg .gt. 0.5
               ! Negative residual means that vcycle did not converge
               if (glob_residnorm .lt. 0) converging = .false.
               !if(ltrace2)write(*,99)myid,'ell_solve:converg',converging
               goto 23
               !
            end if
            !
            if (.not.converging .and. FMGlev .ne. 0) then
               ! If solution not converging we have a big problem
               ! So far this seems to indicate nonsmoothness in
               ! the solution on the finest level
               ! Try smoothing that level and restarting:
               if (myid.eq.master) then
                   write(*,92)myid,'ell_solve: Vcycling: ',
     *             cycle, glob_residnorm
                   write(*,98)myid,'ell_solve:Not converging',percentchg
 92                format('[',I4,'] ',A,I3,F15.9)
               end if
               smoothattempts = smoothattempts + 1
               if (smoothattempts .eq. 22 .and.
     *             smoothattempts .le. maxsmoothattempts) then
                  ! Try relaxing first:
                  if (myid.eq.master) write(*,99)myid,
     *                'ell_solve: Trying relaxing: ',smoothattempts
                  residnorm = level_return_resnrm(FMGlev)
                  call MPI_AllReduce(residnorm,oresidnorm,1,
     *              MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
                  do i = 1, maxsweeps
                     !
                     call level_relax_local(FMGlev)
                     call level_syncbndi_local(FMGlev)
                     !
                     if (mod(i,relax_outperiod).eq.0) then
                        residnorm = level_return_resnrm(FMGlev)
                        call MPI_AllReduce(residnorm,glob_residnorm,1,
     *              MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
                        if (myid.eq.master) write(*,91)myid,
     *                      'ell_solve: i, residnorm= ',i,glob_residnorm
                        percentchg = 100.*(oresidnorm-glob_residnorm)
     *                                    /oresidnorm
                        converging = percentchg .gt. 0.5
                        if (glob_residnorm .lt. relax_resid .and.
     *                      converging) goto 67
                     end if
                  end do
 67               if(.false.)write(*,99)myid,'ell_solve:Relax done:',i-1
                  !
               !else if (smoothattempts .le. maxsmoothattempts) then
                else if (smoothattempts .le. maxsmoothattempts .or.
     *                  abs(glob_residnorm) .gt. 1.7  ) then
                  if (myid.eq.master.or.ltrace2) write(*,99)myid,
     *                'ell_solve: Trying smoothing:',smoothattempts
                  ! Otherwise, try explicit smoothing:
                  call level_mgsmooth(FMGlev,smoothattempts)
                  call level_syncbndi_local(FMGlev)
               end if
               !
               ! Go back now and see if vcycles will work:
               !
               !cycle          = 0
               if (smoothattempts .le. maxsmoothattempts) then
c                 if (myid.eq.master) write(*,*)myid,
c    *                'ell_solve: Smoothattempt:',smoothattempts
                  converging = .true.
                  goto 23
               else
                  if (myid.eq.master) write(*,99)myid,
     *            'ell_solve: Giving up... PROBLEM!!!',FMGlev
c    *            'ell_solve: Giving up... PROBLEM!!!',maxsmoothattempts
               end if
            end if
            if (cycle.ge.maxnumvcyc.and.myid.eq.master)
     *       write(*,98)myid,'ell_solve: Reach max num ',percentchg
            !
            ! Prolong solution to next level
            !    but only if we have a good solution, otherwise,
            !    better just to use previous solution
            !
            if (FMGlev.lt.finest .and. converging .and. cycle.gt.0) then
               if(ltrace2)write(*,99)myid,'ellsolve:Prolong to',FMGlev+1
               call level_mgprolongto(FMGlev+1)
               call level_syncbndi_local(FMGlev+1)
            end if
            !
         end do
         !
      end if    ! if-then-else: solve by vcycling

      !
      ! Do any postprocessing necessary:
      !
      if (ltrace2.and.myid.eq.master)
     *         write(*,99)myid,'ell_solve: Calling ellpost.'
      do i = finest, 1, -1
         if (ltrace2.and.myid.eq.master)
     *         write(*,99)myid,'ell_solve: Inject solution up',i
         call level_mginject_local(i-1)
      end do
      do i = 0, finest
         if (ltrace2.and.myid.eq.master)
     *         write(*,99)myid,'ell_solve: ellpost on: ',i
         call level_ellpost_local(i)
      end do

      if (vtrace) call level_mg_vtrace(finest,2*maxnumvcyc)

      if (ltrace.and.myid.eq.master) then
         write(*,93)myid,'ell_solve: Vcycles: ',
     *             cycle, ' w/ resid: ',glob_residnorm
 93                format('[',I4,'] ',A,I3,A,F15.9)
      end if
      if (ltrace2) write(*,99)myid,'ell_solve: Done.'

      !if (lev.eq.1) call my_exit('ell_solve: quitting')
      return
      end      ! END: ell_solve

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  vcycle:                                                                   cc
cc          Right now, just using AMR hierarchy. At some point,               cc
cc          will need to create appropriate memory structures for             cc
cc          "grids" coarser than the AMR hierarchy's coarsest grid.           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vcycle(l,lcoarse,glob_resnrm)
      implicit      none
      integer       l, lcoarse
      real(kind=8)  glob_resnrm
      include      'glob.inc'
      include      'action.inc'
      include      'grid_methods.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      include      'largesmall.inc'
      integer       m, sweep
      real(kind=8)  resnrm, prev_resnrm, change, percentchg, fineresid
      logical       converging
      !
      ! Minimum fractional change in residual norm on coarse grid
      ! before determining the process failed:
      !
      real(kind=8)  changethresh
      parameter   ( changethresh = -0.00000000001d0 )
      !parameter   ( changethresh = 0.05d0 )

      !
      ! Output coarse grid results with this frequency:
      !
      integer       cg_period
      parameter (   cg_period  =40 )

      logical      vtrace
      parameter (  vtrace  = .false. )
      logical      ltrace
      parameter (  ltrace  = .true. )
      logical      ltrace2
      parameter (  ltrace2 = .false. )

      if (ltrace.and.myid.eq.master) then
 98      format('[',I4,'] ',A,3F15.9)
 99      format('[',I4,'] ',A,3I7)
       write(*,*)'vcycle:                        / / / ******* \\ \\ \\'
      end if
      if (ltrace2.and.myid.eq.master) then
         write(*,*) 'vcycle:                       l      = ',l
         write(*,*) 'vcycle:                       lcoarse= ',lcoarse
         write(*,*) 'vcycle:                       preswp = ',preswp
         write(*,*) 'vcycle:                       pstswp = ',pstswp
         write(*,*) 'vcycle:'
      end if
      !-----------------------------------
      !    Coarse grid correction
      !-----------------------------------
      !if (vtrace) call level_mg_vtrace(0,0*preswp)
      do m = l, 1, -1
         if (vtrace) call level_mg_vtrace(l,0*preswp)
         !
         ! Pre-CGC relaxation sweeps:
         !
         if(ltrace2)write(*,99)myid,'vcycle: Presweeping  level: ',m
         prev_resnrm = LARGENUMBER
         converging  = .true.
         do sweep = 1, preswp
            !
            call level_relax_local(m)
            call level_syncbndi_local(m)
            !
            if (ltrace) then
               resnrm = level_return_resnrm(m)
               call MPI_AllReduce(resnrm,glob_resnrm,1,
     *              MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
               if (myid.eq.master)
     *            call output_resid(m,'Pre-CGC',glob_resnrm)
               if (vtrace) call level_mg_vtrace(m,sweep)
               percentchg = 100.*(prev_resnrm-glob_resnrm)
     *                                    /prev_resnrm
               converging = percentchg .gt. changethresh
               prev_resnrm = glob_resnrm
               if (.not. converging) then
                  !
                  ! If not converging quit early
                  ! before things get worse:
                  !
                  ! If we are not on finest level,
                  ! then recover previous residual:
                  if (m.ne.l) glob_resnrm = fineresid
                  !
                  ! Put negative sign on the residual
                  ! so that ell_solve() knows not
                  ! to try again:
                  !
                  glob_resnrm = - glob_resnrm
                  if (myid.eq.master.or.ltrace2) write(*,99)myid,
     *               'vcycle: Not converging on level:',m
                  return
               end if
            end if
            if (m.eq.l) fineresid = glob_resnrm
         end do
         !
         !
         !  (a)      Compute rhs on (m-1) level:
         !
         !              2h      2h 2h h     2h       h    h h
         !           rhs   :=  L  I  u   + I   (  rhs  - L u  )
         !                          h        h
         !
         !  (b)      Restrict u and transfer to the (m-1) level
         !              
         !            2h         2h  h
         !           u       := I   u
         !                       h
         !
         if(ltrace2)write(*,99)myid,'vcycle: Computing    CGC level: ',m
         !
         call level_mg_cgccompute_local(m)

         !if (vtrace) call level_mg_vtrace(m,10*preswp)
         if(ltrace2)write(*,99)myid,'vcycle: Transferring CGC level: ',m

         call level_cgctransfer_local(m)

         if (vtrace) call level_mg_vtrace(m,10*preswp+1)
         !
      end do

      !-----------------------------------
      !    Solve coarse grid problem
      !-----------------------------------
      ! 
      if (vtrace) call level_mg_vtrace(0,10*preswp+2)
      if (ltrace2) write(*,*) 'vcycle: Solving coarse grid problem'
      sweep       = 0
      glob_resnrm = 2.d0*ell_epsilon
      glob_resnrm = LARGENUMBER
      change      = LARGENUMBER
      do while (  glob_resnrm .gt. ell_epsilon .and.
     *            sweep       .lt. maxsweeps   .and.
     *            change      .gt. changethresh      )
            !
            call level_relax_local(0)
            call level_syncbndi_local(0)
            !
            if (mod(sweep,cg_period).eq.0) then
               prev_resnrm = glob_resnrm
               resnrm      = level_return_resnrm(0)
               !if(ltrace2)write(*,98)myid,'vcycle: allreduce:',resnrm
               call MPI_AllReduce(resnrm,glob_resnrm,1,
     *              MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
               if (ltrace.and.myid.eq.master.or.ltrace2)
     *               call output_resid(0,'CG Solve:',glob_resnrm)
               change = ( prev_resnrm - glob_resnrm ) / glob_resnrm
            end if
            sweep = sweep + 1
      end do
      !if (myid.eq.master) write(*,*) prev_resnrm, change
      if (ltrace.and.myid.eq.master.or.ltrace2)
     *    write(*,13) 'Coarse grid solve took: ',sweep,
     *                   ' sweeps w/ resnrm= ',glob_resnrm
 13   format(A,I6,a,1p,E9.2)
      !if (glob_resnrm .gt. ell_epsilon.and.myid.eq.master) then
      if (glob_resnrm .gt. ell_epsilon) then
         if (myid.eq.master)
     *          write(*,*) 'vcycle: Coarse grid solve failed.'
         ! Put negative sign to indicate failure
         glob_resnrm = - glob_resnrm
         !call my_exit('vcycle: Coarse grid solve failed.')
      end if
      if (vtrace) then
          if(ltrace2)write(*,*) myid, 'outputting vtrace at ',20*preswp
          call level_mg_vtrace(0,20*preswp)
      end if

      !-----------------------------------
      !    Post-CGC fix & relaxation
      !-----------------------------------
      do m = 1, l
         !           Linearly interpolate correction and update
         !           level m unknown.
         !
         !            h       h    h    2h     2h  h
         !           u    := u  + I  ( u   -  I   u  )
         !                         2h          h
         !
         if (ltrace2) write(*,*) 'vcycle: Computing Post-CGC level: ',m
         !
         call level_psttransfer_local(m)
         if (ltrace2) write(*,*) 'vcycle: Computing Post-CGC compute'
         call level_mg_pstcompute_local(m)
         if (ltrace2) write(*,*) 'vcycle: Computing Syncing'
         call level_syncbndi_local(m)
         if (vtrace) call level_mg_vtrace(m,100*preswp)
         !
         !
         ! Post-CGC relaxation sweeps:
         !
         if (ltrace2) write(*,*) 'vcycle: Postsweeping level: ',m
         do sweep = 1, pstswp
            !
            call level_relax_local(m)
            call level_syncbndi_local(m)
            !
            if (ltrace) then
               resnrm = level_return_resnrm(m)
               call MPI_AllReduce(resnrm,glob_resnrm,1,
     *              MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
               if (ltrace.and.myid.eq.master.or.ltrace2)
     *              call output_resid(m,'Post-CGC',glob_resnrm)
               if (vtrace) call level_mg_vtrace(m,200*preswp+sweep)
            end if
         end do
         !
      end do


      if (ltrace) then
      !if (ltrace .and. myid.eq.master) then
       call level_get_resid(l)
       if (myid .eq. master) 
     . write(*,*)'vcycle:                        \\ \\ \\ ******* / / /'
      end if
      if (vtrace) call level_mg_vtrace(l,300*preswp+sweep)

      !if (l.eq.1) call my_exit('vcycle:byebye')

      return
      end      ! END: vcycle

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mg_vtrace:   Has each grid on a level output tracing fields for     cc
cc                     multigrid routine.                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mg_vtrace(l,count)
      implicit      none
      integer       l, count
      integer       level_return_start, grid_return_sibling
      external      level_return_start, grid_return_sibling
      integer       myid, gi
      integer       proc_return_myid
      external      proc_return_myid
      logical       grid_return_existence, grid_is_local
      external      grid_return_existence, grid_is_local
      real(kind=8)  grid_return_time
      external      grid_return_time
      !
      logical       ltrace
      parameter (   ltrace  = .false. )
      logical       ltrace2
      parameter (   ltrace2 = .false. )


 97      format('[',I4,'] ',A,3F15.9)
 98      format('[',I4,'] ',A,3F15.9)
 99      format('[',I4,'] ',A,I7,A,F15.9,A,I7)

      myid = proc_return_myid()
      if (ltrace2) then
         write(*,99) myid, 'level_mg_trace: Level: ',l
         write(*,99) myid, 'level_mg_trace: Count: ',count
      end if

      gi = level_return_start(l)
      if (myid.eq.0.and.grid_return_existence(gi).and.ltrace)
     *        write(*,99) myid,'level_mg_trace: Outputing lev: ',l,
     *          'at t= ',grid_return_time(gi),' count=',count
 7    if (grid_return_existence(gi)) then
          if (grid_is_local(gi)) then
             call grid_mg_trace(gi,count)
          end if
          gi = grid_return_sibling(gi)
          goto 7
      end if

      if (ltrace2) then
         write(*,98)myid,'level_mg_trace: Done.'
      end if

      return
      end       ! END: level_mg_trace

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_return_resnrm:                                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function level_return_resnrm(l)
      implicit      none
      integer       l
      integer       gi, numpoints, totalpts
      integer       nx, ny, nz
      real(kind=8)  norm
      double precision  grid_return_resnrm, normDP
      external      grid_return_resnrm
      integer       level_return_start, grid_return_sibling
      external      level_return_start, grid_return_sibling
      integer       myid
      integer       proc_return_myid
      external      proc_return_myid
      logical       grid_return_existence, grid_is_local
      external      grid_return_existence, grid_is_local
      !
      logical       ltrace
      parameter (   ltrace = .false. )


 98      format('[',I4,'] ',A,3G15.9)
 99      format('[',I4,'] ',A,3I7)

      if (ltrace) then
         myid = proc_return_myid()
         write(*,99) myid, 'level_return_resnrm: Level: ',l
      end if

      totalpts            = 0
      level_return_resnrm = 0.d0

      gi = level_return_start(l)
 7    if (grid_return_existence(gi)) then
          if (grid_is_local(gi)) then
             if(ltrace)write(*,99)myid,'level_return_resnrm:gi=',gi
             call grid_get_dims(gi,nx,ny,nz)
             numpoints = nx * ny * nz
             totalpts  = totalpts + numpoints
             !write(*,*) 'level_return_resnrm: About to call:'
             !write(*,*) 'normDP = ',normDP
             !write(*,*) 'norm   = ',norm
             !write(*,*) 'gi     = ',gi
             !write(*,*) 'level_return_resnrm: ',grid_return_resnrm(gi)
             norm      = grid_return_resnrm(gi)
             !write(*,*)myid,'level_return_resnrm:norm=',norm
             if(ltrace)write(*,98)myid,'level_return_resnrm:norm=',norm
             level_return_resnrm = level_return_resnrm+numpoints*norm**2
          end if
          gi = grid_return_sibling(gi)
          goto 7
      end if

      if(ltrace)write(*,98)myid,'level_return_resnr',level_return_resnrm
      if (totalpts .ne. 0) then
         if(ltrace)write(*,99)myid,'level_return_resnrm:total=',totalpts
         level_return_resnrm = sqrt(level_return_resnrm / totalpts)
      end if

      if (ltrace) then
         write(*,98)myid,'level_return_resnrm:Done ',level_return_resnrm
      end if

      return
      end       ! END: level_return_resnrm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_get_resid:  Computes residual and sticks it into the _st1 fields.   cc
cc                    Used just for being able to output the resid (for now). cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_get_resid(l)
      implicit      none
      integer       l
      integer       gi, numpoints, totalpts
      integer       nx, ny, nz
      real(kind=8)  norm
      double precision  grid_return_resnrm
      external      grid_return_resnrm
      integer       level_return_start   , grid_return_sibling
      external      level_return_start   , grid_return_sibling
      logical       grid_return_existence, grid_is_local
      external      grid_return_existence, grid_is_local

      logical       ltrace
      parameter (   ltrace = .false. )


      if (ltrace) write(*,*) 'level_get_resid:     Level: ',l

      gi = level_return_start(l)
 7    if (grid_return_existence(gi)) then
          if (grid_is_local(gi)) then
             call grid_get_resid(gi)
          end if
          if (ltrace) write(*,*) 'level_get_resid:    gi = ',gi        
          gi = grid_return_sibling(gi)
          goto 7
      end if

      return
      end       ! END: level_get_resid

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_bbox_write:    Write to file the clustering information             cc
cc                    NB: Updating to handle cases where symmetry assumed.    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_bbox_write(li,b_minx,b_maxx,b_miny,b_maxy,
     *                           b_minz,b_maxz,mx,my,mz,myh,time,numbox)
      implicit    none
      integer     li,        numbox
      real(kind=8)mx,my,mz, myh, time
      integer     b_minx(*), b_maxx(*),
     *            b_miny(*), b_maxy(*),
     *            b_minz(*), b_maxz(*)
      integer     filenum,   i
      character(80) junk
      include     'param.inc'

      logical     ltrace
      parameter ( ltrace = .false. )

      !if (numbox.le.0) return

      if (ltrace) then
         write(*,*) 'level_bbox_write: Recording for level:',li
         write(*,*) 'level_bbox_write: min x/y/z: ',mx,my,mz
         write(*,*) 'level_bbox_write: myh:       ',myh
      end if

      filenum = 9
      open(filenum,file='.cluster',form='formatted',status='unknown')

      ! Get to the end of the file to append info:
         if(ltrace)write(*,*) 'level_bbox_write: Getting to EOF:'
         junk = "Z"
  6      format(A80)
  8      read(filenum,6,end=10) junk
         if(ltrace)write(*,*) 'level_bbox_write: ',junk
         goto 8 
 10      continue  
         if(ltrace)write(*,*) 'level_bbox_write: Now at EOF:'

      if (ltrace) write(*,*) 'level_bbox:time,li,numbox:',time,li,numbox

      rewind(filenum)
      write(filenum,99) time, li, numbox, nx0,ny0,nz0
      write(filenum,89) mx,my,mz,myh
  99  format(F18.9,5I7)
  89  format(4F18.9)

      do i = 1, numbox
         write(filenum,98) b_minx(i),b_miny(i),b_minz(i),
     *                     b_maxx(i),b_maxy(i),b_maxz(i)
         if (ltrace)
     *         write(*,97) b_minx(i),b_maxx(i),
     *                     b_miny(i),b_maxy(i),
     *                     b_minz(i),b_maxz(i)
      end do
  98  format(6I8)
  97  format(2I4,'  ',2I4,'  ',2I4)

      close(filenum)
      if (ltrace) write(*,*) 'level_bbox_write: done.     for level:',li

      return
      end      ! END: level_bbox_write

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_bbox_read:     Read from file the clustering information            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_bbox_read(li,b_minx,b_maxx,b_miny,b_maxy,
     *                           b_minz,b_maxz,mx,my,mz,myh,time,numbox,
     *                           assume,status)
      implicit    none
      integer     li,        numbox, assume, status
      real(kind=8)mx,my,mz, myh, time
      integer     b_minx(*), b_maxx(*),
     *            b_miny(*), b_maxy(*),
     *            b_minz(*), b_maxz(*)
      real(kind=8)      tmp_time, factor
      integer     tmp_li, refine, num,
     *            tmp_nx0, tmp_ny0, tmp_nz0
      real(kind=8)tmp_mx,tmp_my,tmp_mz, tmp_myh
      integer     filenum,   i, rc, clusterzero, level_zero
      character(1) junk
      logical      sametime
      logical      double_equal
      external     double_equal
      include     'param.inc'
      !include     'largesmall.inc'

      logical     ltrace
      parameter ( ltrace = .true. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )
      ! Return codes:
      !     success--->found an entry and read in boxes
      !     failure--->found an entry but could not read in boxes
      !     notfound-->could not find the appropriate entry
      integer    SUCCESS,    FAILURE,    NOTFOUND
      parameter (SUCCESS=1,  FAILURE=-1, NOTFOUND=0)

      !
      ! If no indications to create a subgrid, then none should exist:
      !
      numbox = 0

      if (ltrace2) then
         write(*,*) 'level_bbox_read: Looking for level: ', li
         write(*,*) 'level_bbox_read:          at time: ', time
         write(*,*) 'level_bbox_read: min x/y/z: ',mx,my,mz
         write(*,*) 'level_bbox_read: h: ', h
      end if

      status  = FAILURE

      filenum = 9
      open(filenum,file='.cluster',status='old',iostat=rc)
      if(rc .ne. 0) then
        write(*,*) 'level_bbox_read: Unable to open .cluster file'
        call my_exit('Unable to open .cluster file.')
      end if

  5   read(filenum,*,end=10) tmp_time,tmp_li,num,tmp_nx0,tmp_ny0,tmp_nz0
      read(filenum,*,end=10) tmp_mx,tmp_my,tmp_mz,tmp_myh
      if (ltrace2) write(*,*) '  An entry:',tmp_time,tmp_li, num
      if (ltrace2) write(*,*) '       and:',tmp_mx,tmp_my,tmp_mz,tmp_myh
      !sametime = double_equal(tmp_time,time)
      ! Need a less stringent tolerance:
      sametime = abs(tmp_time-time) .le. 1.0d-7
      if (sametime .and. tmp_li .eq. li) then
         if(ltrace)write(*,*)' Our Entry:  t/li/num',tmp_time,tmp_li,num
         numbox = num
         do i = 1, numbox
            read(filenum,*,end=10) b_minx(i),b_miny(i),b_minz(i),
     *                             b_maxx(i),b_maxy(i),b_maxz(i)
            if(ltrace2)write(*,*)' Read a box:',i,
     *                             b_minx(i),b_miny(i),b_minz(i),
     *                             b_maxx(i),b_maxy(i),b_maxz(i)
         end do
         status = SUCCESS
         if(ltrace2)write(*,*)' Done reading boxes',numbox
      else
         if (ltrace2) write(*,*) '   Not our entry',
     *           tmp_time - time,
     *           tmp_li   .eq. li
         do i = 1, num
            read(filenum,*,end=10) junk
         end do
         if (ltrace2) write(*,*) '   Are we past our time?'
         if (tmp_time .le. time .or. sametime) then
            go to 5
         else
            status = NOTFOUND
         end if
      end if

 10   continue  

      if (status.eq.FAILURE.and.tmp_time.gt.time) then
      !if (status.eq.FAILURE) then
         write(*,*) 'level_bbox_read: Unable to read in all the boxes'
         call my_exit('Unable to read in all the boxes')
      end if

      if (status.eq.SUCCESS) then
         !
         ! Rationalize different refinement factors
         !   (assumes right now that resolution same for all dimensions)
         !
         if(ltrace2)write(*,*)' Rationalize for different coarse level?'
         if(ltrace2)write(*,*)' nx0 versus tmp_nx0 = ',nx0, tmp_nx0
         if (nx0.ne.tmp_nx0) then
            do i = 1, numbox
               ! X-axis:
               factor    = ( (nx0-1.d0) / (tmp_nx0-1.d0) )
               if(ltrace2)write(*,*)'X Before: ',b_minx(i),b_maxx(i)
               if (assume.eq.1.or.assume.eq.6) then
                  if(ltrace2)write(*,*)' Rationalizing w/symm in x...'
                  ! Get indices of zero points:
                  clusterzero = (0-tmp_mx)/tmp_myh + 1
                  level_zero  = (0-    mx)/    myh + 1
                  if(ltrace2)write(*,*)' level/clusterzero:',
     *                                   level_zero,clusterzero
                  if (b_minx(i) .lt. clusterzero) then
                     ! If box begins less than zero,
                     ! Then rationalized box begins at least coordinate value:
                     b_minx(i) = 1
                  else
                     ! Otherwise
                     b_minx(i) = NINT((b_minx(i)-clusterzero)*factor)
     *                        + level_zero
                  end if
                  b_maxx(i) = NINT((b_maxx(i)-clusterzero)*factor)
     *                     + level_zero
               else
                  b_minx(i) = NINT( (b_minx(i)-1) * factor) + 1
                  b_maxx(i) = NINT( (b_maxx(i)-1) * factor) + 1
               end if
               if(ltrace2)write(*,*)'X After: ',b_minx(i),b_maxx(i)
               ! Y-axis:
               factor    = ( (ny0-1.d0) / (tmp_ny0-1.d0) )
               if(ltrace2)write(*,*)'Y Before: ',b_miny(i),b_maxy(i)
               if (assume.eq.2.or.assume.eq.6) then
                  if(ltrace2)write(*,*)' Rationalizing w/symm in y...'
                  ! Get indices of zero points:
                  clusterzero = (0-tmp_my)/tmp_myh + 1
                  level_zero  = (0-    my)/    myh + 1
                  if(ltrace2)write(*,*)' level/clusterzero:',
     *                                   level_zero,clusterzero
                  if (b_miny(i) .lt. clusterzero) then
                     ! If box begins less than zero,
                     ! Then rationalized box begins at least coordinate value:
                     b_miny(i) = 1
                  else
                     ! Otherwise
                     b_miny(i) = NINT((b_miny(i)-clusterzero)*factor)
     *                        + level_zero
                  end if
                  b_maxy(i) = NINT((b_maxy(i)-clusterzero)*factor)
     *                     + level_zero
               else
                  b_miny(i) = NINT( (b_miny(i)-1) * factor) + 1
                  b_maxy(i) = NINT( (b_maxy(i)-1) * factor) + 1
               end if
               if(ltrace2)write(*,*)'Y After: ',b_miny(i),b_maxy(i)
               ! Z-axis:
               !factor    = ( (nz0-1.d0) / (tmp_nz0-1.d0) )
               factor    = tmp_myh / myh
               if(ltrace2)write(*,*)'factor: ',factor
               if(ltrace2)write(*,*)'Z Before: ',b_minz(i),b_maxz(i)
               if (assume.eq.3.or.assume.eq.6) then
                  if(ltrace2)write(*,*)' Rationalizing w/symm in z...'
                  ! Get indices of zero points:
                  clusterzero = (0-tmp_mz)/tmp_myh + 1
                  level_zero  = (0-    mz)/    myh + 1
                  if(ltrace2)write(*,*)' level/clusterzero:',
     *                                   level_zero,clusterzero
                  if (b_minz(i) .lt. clusterzero) then
                     ! If box begins less than zero,
                     ! Then rationalized box begins two points before its zero:
                     b_minz(i) = max(1,level_zero-2)
                  else
                     ! Otherwise
                     b_minz(i) = NINT((b_minz(i)-clusterzero)*factor)
     *                        + level_zero
                  end if
                  b_maxz(i) = NINT((b_maxz(i)-clusterzero)*factor)
     *                     + level_zero
               else
                  b_minz(i) = NINT( (b_minz(i)-1) * factor) + 1
                  b_maxz(i) = NINT( (b_maxz(i)-1) * factor) + 1
               end if
               if(ltrace2)write(*,*)'Z After: ',b_minz(i),b_maxz(i)
            end do
         end if
      end if

      close(filenum)

      if (ltrace) then
         write(*,*)'level_bbox_read: Final box(es):',numbox
         do i = 1, numbox
               write(*,*) '..........'
               write(*,99) b_minx(i),b_maxx(i),
     *                     b_miny(i),b_maxy(i),
     *                     b_minz(i),b_maxz(i)
               write(*,*) '..........'
         end do
      end if
      if(ltrace)write(*,*)'level_bbox_read: Done.status=',status

 99   format(2I4,'  ',2I4,'  ',2I4)
      return
      end      ! END: level_bbox_read

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mgprolongto:                                                        cc
cc          Take solution for elliptic problem on one level                   cc
cc          and interpolate/prolong to the next finer level.                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mgprolongto(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi, gj, gi_owner, gj_owner
      integer      round, numrounds, proc
      integer      round_robin
      external     round_robin
      logical      gi_local, gj_local
      logical      ltrace
      parameter (  ltrace = .true. )

      if (ltrace) then
         write(*,99) myid, 'level_mgprolongto: level = ',level
      end if

      !
      ! We are processor "myid" and we need to check with all other
      ! processors to see if we need to sync with them. We want to do
      ! so in a way that avoids having a processor wait for the other
      ! processor to finish with someone else. So we use a round robin.
      ! 
      ! For each "round" of the round robin, we determine which processor
      ! with which to sync.

      if (mod(numprocs,2).eq.0) then
         ! EVEN number of processors:
         numrounds = numprocs - 1
      else
         ! ODD  number of processors:
         numrounds = numprocs
      end if
      ! Add a round to sync pairs of grids both owned by this process:
      numrounds = numrounds + 1

      !
      ! Loop over the rounds of the round robin:
      !
      do round = 1, numrounds
        if(ltrace)write(*,99)myid,'level_psttransfer_local:round=',round
        !
        ! Determine which processor with which to sync:
        !
        if (round .ne. numrounds) then
           proc = round_robin(round, myid, numprocs)
        else
           proc = myid
        end if
        !
        ! Make sure we don't have a "bye" (only happens for odd numprocs):
        !
        if (proc .ge. 0) then
           if (ltrace) then
              write(*,99)myid,'level_mgprolongto: w/proc=',proc
           end if
           !
           ! Loop for all grids (gi) on level. For each gi,
           ! consider pairings of gi and gj where occurs in the
           ! hierarchy after gi (this avoids duplicate pairings),
           ! and sync those pairings (gi,gj):
           !
           gi = level_return_start(level)
           !
 10        if (grid_return_existence(gi)) then
              if (ltrace) write(*,99) myid, '      gi = ',gi
              gi_owner = grid_return_owner(gi)
              gi_local = grid_is_local(gi)
              ! If neither us nor processor "proc" own gi, then skip:
              if (.not. gi_local .and. gi_owner .ne. proc) goto 40
              !
              !gj = grid_return_sibling(gi)
              ! Loop over parent level:
              gj = level_return_start(level-1)
 20           if (grid_return_existence(gj)) then
                 !
                 ! Should we skip this pairing?
                 gj_owner = grid_return_owner(gj)
                 gj_local = grid_is_local(gj)
                 if (.not. gj_local .and. gj_owner .ne. proc) goto 30
                 if (gi_local.and.gj_local .and. (myid.ne.proc)) goto 30
                 ! Does pair intersect?
                 if (ltrace) write(*,99) myid, '      w/ gj = ',gj
                 if ( grid_intersect(gi,gj) ) then
                    if (gj_local) then
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi/gj local'
                          call grid_initMG(gi, gj,2)
                          !call grid_initMG(gi, gj,refine_factor)
                       else
                          if (ltrace) write(*,99) myid, '   gj   local'
                          call grid_init_parentMG(gi, gj)
                       end if
                    else
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi   local'
                          !call grid_psttransfer_slave(gi, gj)
                          call grid_init_childMG(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   none local'
                       end if
                    end if
                 else
                          if (ltrace) write(*,99) myid, '   No intersec'
                 end if
                 ! Look to next sibling
  30             gj = grid_return_sibling(gj)
                 goto 20
              end if
              !
              ! No more siblings w/ which to sync gi,
              ! so, repeat process for next grid on level:
              !
  40          gi = grid_return_sibling(gi)
              goto 10
           end if

        end if
      end do

      if (ltrace) write(*,99) myid, 'level_mgprolongto: Finished.'

  99  format('[',I4,'] ',A,3I5)
      return 
      end    ! END: level_mgprolongto
