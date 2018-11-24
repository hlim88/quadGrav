cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  tracers_init:                                                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine tracers_init()
       implicit     none
       include     'glob.inc'
       include     'param.inc'
       include     'tracers.inc'
       !include     'mem.inc'
       !include     'mask.inc'
       ! Only for testing routine call:
       include     'grid.inc'
       include     'mpif.h'
       include     'mpi_stuff.inc'
       !
       integer      i, j, k, mynx,myny
       real*8       mydx, mydy                
       ! values that will get overwritten with MPI call:
       real*8       UNDEFINED
       parameter (  UNDEFINED = 1d99)
       logical     ltrace
       parameter ( ltrace  = .true. )
       logical     ltrace2
       parameter ( ltrace2  = .false. )

      if ((ltrace.and.myid.eq.0).or.ltrace2) then
         write(*,99) myid, 'tracers_init: Enter: ',tracers_number
      end if

      ! Testing that the interp and interpderiv() routines work:
      if (.false.) then
      call test_interproutines(q(gfunc_pointer(152,1)),
     .               q(gfunc_pointer(153,1)),tracers_x,
     .               gr_nx(1),gr_ny(1),gr_nz(1),tracers_number)
      end if

      tracers_t = 0.d0

      mynx = 1
      myny = 1

      do i = 1, tracers_number
         tracers_x(i) = UNDEFINED
         tracers_y(i) = UNDEFINED
         tracers_z(i) = UNDEFINED
      end do

      ! Some options can be accomplished by master process alone:
      if (.false. .and. myid.eq.0) then
         !
         if (ltrace) then
            write(*,99)myid,'tracers_init: Master initting tracers'
         end if
         if (.true.) then
            !
            ! Distribute across mid-plane
            !
            mynx = sqrt(1.d0*tracers_number)
            myny = tracers_number/mynx
            if (myny*mynx>tracers_number)myny=myny-1
            mydx = (maxx0-minx0)/(mynx-1)
            mydy = (maxy0-miny0)/(myny-1)
            ! For right now just distribute across plance, but
            ! in the future, will want to distribute consistent with
            ! the density so that density of particles corresponds to
            ! fluid density. Or, may want to distribute only in low
            ! density regions
            k    = 0
            do j = 1, myny
            do i = 1, mynx
               k = k + 1
               tracers_x(k) = minx0 + (i-1)*mydx
               tracers_y(k) = miny0 + (j-1)*mydy
               tracers_z(k) = 0.d0
            end do
            end do
         else if (.true.) then
            !
            ! Distribute across the x-axis:
            !
            mydx = (maxx0-minx0)/(tracers_number-1)
            do i = 1, tracers_number
               tracers_x(i) = minx0 + (i-1)*mydx
               tracers_y(i) = 0.d0
               tracers_z(i) = 0.d0
            end do
         end if
         !
      else
         !
         ! All processors participate:
         !
         if (ltrace2) then
            write(*,99)myid,'tracers_init: All processes initing'
         end if
         call tracers_init_sub()
         if (ltrace2) then
            write(*,99)myid,'tracers_init: Back from tracers_init_sub'
         end if
      end if

      !
      !
      if (ltrace2) then
         write(*,99) myid, 'tracers_init:        MPI_Allreduce:'
      end if
      call mat_copy1d(tracers_x,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_x,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_x,tracers_number)
      !
      call mat_copy1d(tracers_y,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_y,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_y,tracers_number)
      !
      call mat_copy1d(tracers_z,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_z,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_z,tracers_number)
      !
      call mat_copy1d(tracers_Vx,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_Vx,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_Vx,tracers_number)
      !
      call mat_copy1d(tracers_Vy,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_Vy,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_Vy,tracers_number)
      !
      call mat_copy1d(tracers_Vz,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_Vz,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_Vz,tracers_number)
      if (ltrace2) then
         write(*,99) myid, 'tracers_init: Done w/MPI_Allreduce.'
      end if

      ! In case not all particles were defined, simply reset
      ! un-initialized ones:
      do i = 1, tracers_number
         if (tracers_x(i) .ge. UNDEFINED) tracers_x(i)=0.d0
         if (tracers_y(i) .ge. UNDEFINED) tracers_y(i)=0.d0
         if (tracers_z(i) .ge. UNDEFINED) tracers_z(i)=0.d0
         if (tracers_Vx(i).ge. UNDEFINED) tracers_Vx(i)=0.d0
         if (tracers_Vy(i).ge. UNDEFINED) tracers_Vy(i)=0.d0
         if (tracers_Vz(i).ge. UNDEFINED) tracers_Vz(i)=0.d0
      end do

      ! Initialize n-1 time data:
      call mat_copy1d(tracers_x, tracers_xprev,  tracers_number)
      call mat_copy1d(tracers_y, tracers_yprev,  tracers_number)
      call mat_copy1d(tracers_z, tracers_zprev,  tracers_number)
      !
      call mat_copy1d(tracers_Vx,tracers_Vxprev, tracers_number)
      call mat_copy1d(tracers_Vy,tracers_Vyprev, tracers_number)
      call mat_copy1d(tracers_Vz,tracers_Vzprev, tracers_number)

      if ((ltrace.and.myid.eq.0).or.ltrace2) then
         call field_dump_info1D(tracers_x,     tracers_number,
     .                                         'tracers_x')
         call field_dump_info1D(tracers_xprev, tracers_number,
     .                                         'tracers_xprev')
         call field_dump_info1D(tracers_Vx,    tracers_number,
     .                                         'tracers_Vx')
         call field_dump_info1D(tracers_Vxprev,tracers_number,
     .                                         'tracers_Vxprev')
         write(*,99) myid, 'tracers_init: tracers_initial=',
     .                  tracers_initial
         write(*,99) myid, 'tracers_init: tracers_period =',
     .                  tracers_period 
         write(*,99) myid, 'tracers_init: tracers_scheme =',
     .                  tracers_scheme 
         write(*,99) myid, 'tracers_init: Done.'
         !write(*,98) myid, 'tracers_init: mydx=',mydx
         !write(*,98) myid, 'tracers_init: mydy=',mydy
         !write(*,99) myid, 'tracers_init: mynx=',mynx
         !write(*,99) myid, 'tracers_init: myny=',myny
         !write(*,99) myid, 'tracers_init: k   =',k
      end if

 98    format('[',I4,'] ',A,5F12.6)
 99    format('[',I4,'] ',A,2I5)

       return
       end     ! END: tracers_init

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  tracers_compute:                                                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine tracers_compute(time)
       implicit     none
       real(kind=8) time
       include     'glob.inc'
       include     'param.inc'
       include     'tracers.inc'
       !include     'mem.inc'
       !include     'mask.inc'
       include     'mpif.h'
       include     'mpi_stuff.inc'
       !
       integer      i, j, k, mynx,myny
       integer      li, gi, own
       real*8       mydx, mydy, myx,myy,myz, mydt
       real*8       myVx, myVy, myVz
       real*8       velx,vely,velz
       real*8       myl2norm
       external     myl2norm
       ! values that will get overwritten with MPI call:
       real*8       UNDEFINED
       parameter (  UNDEFINED = 1d99)
       logical     ltrace
       parameter ( ltrace   = .false. )
       logical     ltrace2
       parameter ( ltrace2  = .false. )

      mydt = time - tracers_t

      if(tracers_t.eq.0) then
        ! To ge the leap-frog scheme, take first step with half-timestep
        if(ltrace2)write(*,99)myid,'tracers_compute: Halving mydt'
        mydt = 0.5d0*mydt
      end if

      if (ltrace) then
         write(*,99) myid, 'tracers_compute: Enter: ',tracers_number
         write(*,98) myid, 'tracers_compute: Last t:',tracers_t
         write(*,98) myid, 'tracers_compute: time:  ',time
         write(*,98) myid, 'tracers_compute: mydt:  ',mydt
         call field_dump_info1D(tracers_x,      tracers_number,'X'     )
         call field_dump_info1D(tracers_xprev,  tracers_number,'Xprev' )
         call field_dump_info1D(tracers_Vx,     tracers_number,'Vx'    )
         call field_dump_info1D(tracers_Vxprev, tracers_number,'Vxprev')
         write(*,98) myid, 'tracers_compute: ......................'
      end if

      ! Switch prev and current positions:
      call mat_copy1d(tracers_x,    tracers_tmp,  tracers_number)
      call mat_copy1d(tracers_xprev,tracers_x,    tracers_number)
      call mat_copy1d(tracers_tmp,  tracers_xprev,tracers_number)
      !
      call mat_copy1d(tracers_y,    tracers_tmp,  tracers_number)
      call mat_copy1d(tracers_yprev,tracers_y,    tracers_number)
      call mat_copy1d(tracers_tmp,  tracers_yprev,tracers_number)
      !
      call mat_copy1d(tracers_z,    tracers_tmp,  tracers_number)
      call mat_copy1d(tracers_zprev,tracers_z,    tracers_number)
      call mat_copy1d(tracers_tmp,  tracers_zprev,tracers_number)
      !
      call mat_copy1d(tracers_Vx,    tracers_tmp,  tracers_number)
      call mat_copy1d(tracers_Vxprev,tracers_Vx,   tracers_number)
      call mat_copy1d(tracers_tmp,  tracers_Vxprev,tracers_number)
      !
      call mat_copy1d(tracers_Vy,    tracers_tmp,  tracers_number)
      call mat_copy1d(tracers_Vyprev,tracers_Vy,   tracers_number)
      call mat_copy1d(tracers_tmp,  tracers_Vyprev,tracers_number)
      !
      call mat_copy1d(tracers_Vz,    tracers_tmp,  tracers_number)
      call mat_copy1d(tracers_Vzprev,tracers_Vz,   tracers_number)
      call mat_copy1d(tracers_tmp,  tracers_Vzprev,tracers_number)
      if (ltrace) then
         write(*,98) myid, 'tracers_compute: Switching arrays:'
         write(*,98) myid, 'tracers_compute: ......................'
         call field_dump_info1D(tracers_x,      tracers_number,'X'     )
         call field_dump_info1D(tracers_xprev,  tracers_number,'Xprev' )
         call field_dump_info1D(tracers_Vx,     tracers_number,'Vx'    )
         call field_dump_info1D(tracers_Vxprev, tracers_number,'Vxprev')
         write(*,98) myid, 'tracers_compute: ......................'
      end if

      do i = 1, tracers_number
         ! Get location of tracer particle:
         myx = tracers_xprev(i)
         myy = tracers_yprev(i)
         myz = tracers_zprev(i)
         myVx= tracers_Vxprev(i)
         myVy= tracers_Vyprev(i)
         myVz= tracers_Vzprev(i)
         ! Find finest grid containing this particle:
         call level_finest_containing(myx,myy,myz,li,gi,own)
         if (ltrace .and. own.lt.0) then
         !if (ltrace2) then
            write(*,97)myid,'tracers_compute: ',i,li,gi,own
         end if
         if (own.eq.-1) then
            if(ltrace2)then
               write(*,99)myid,'tracers_compute: Particle not on grid:'
               write(*,98)myid,'tracers_compute: x/y/z:',myx,myy,myz
            end if
            write(*,97)myid,'tracers_compute:off grid:',i,own,gi,li
            write(*,98)myid,'tracers_compute:off x/y/z:',myx,myy,myz
            tracers_x(i) = 0.d0
            tracers_y(i) = 0.d0
            tracers_z(i) = 0.d0
            ! any particles outside of grid should be reset:
         else if (own.eq.myid) then
            if(ltrace2)write(*,97)myid,'tracers_compute:part. local',i
            write(*,97)myid,'tracers_compute:on grid:',i,own,gi,li
            write(*,98)myid,'tracers_compute:on x/y/z:',myx,myy,myz
            ! Because     local, update particle here:
            call grid_tracer_update(gi,i,mydt)
         else
            ! Because not local, set to UNDEFINED 
            !  so that we'll get value from MPI below:
            if(ltrace2)write(*,97)myid,'tracers_compute:not local',i,own
            !write(*,97)myid,'tracers_compute:not local',i,own,gi,li
            tracers_x(i) = UNDEFINED
            tracers_y(i) = UNDEFINED
            tracers_z(i) = UNDEFINED
            tracers_Vx(i)= UNDEFINED
            tracers_Vy(i)= UNDEFINED
            tracers_Vz(i)= UNDEFINED
         end if
      end do

      if (numprocs.gt.1) then
      ! all processors share data (and get rid of undefined values):
      !     have to use different storage...and so use temporary storage
      !     and copy back:
      call mat_copy1d(tracers_x,tracers_tmp,tracers_number)
      if(ltrace2)then
         call field_dump_info1D(tracers_x, tracers_number, 'tracers_x')
      end if
      call MPI_Allreduce(tracers_x,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_x,tracers_number)
      if(ltrace2)then
         call field_dump_info1D(tracers_x, tracers_number, 'tracers_x')
      end if
      !
      call mat_copy1d(tracers_y,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_y,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_y,tracers_number)
      !
      call mat_copy1d(tracers_z,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_z,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_z,tracers_number)
      !
      call mat_copy1d(tracers_Vx,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_Vx,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_Vx,tracers_number)
      !
      call mat_copy1d(tracers_Vy,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_Vy,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_Vy,tracers_number)
      !
      call mat_copy1d(tracers_Vz,tracers_tmp,tracers_number)
      call MPI_Allreduce(tracers_Vz,tracers_tmp,tracers_number,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call mat_copy1d(tracers_tmp,tracers_Vz,tracers_number)
      end if

      ! Update the time for which the particles were last advected
      tracers_t = time
      
      ! In case not all particles were defined, simply reset
      ! un-initialized ones:
      !    SLL: This should not be necessary
      do i = 1, tracers_number
         if (tracers_x(i) .ge. UNDEFINED .or.
     .       tracers_y(i) .ge. UNDEFINED .or.
     .       tracers_z(i) .ge. UNDEFINED .or.
     .       tracers_Vx(i).ge. UNDEFINED .or.
     .       tracers_Vy(i).ge. UNDEFINED .or.
     .       tracers_Vz(i).ge. UNDEFINED     ) then
            tracers_x(i) =0.d0
            tracers_y(i) =0.d0
            tracers_z(i) =0.d0
            tracers_Vx(i)=0.d0
            tracers_Vy(i)=0.d0
            tracers_Vz(i)=0.d0
            tracers_xprev(i) =0.d0
            tracers_yprev(i) =0.d0
            tracers_zprev(i) =0.d0
            tracers_Vxprev(i)=0.d0
            tracers_Vyprev(i)=0.d0
            tracers_Vzprev(i)=0.d0
         end if
      end do

      if (ltrace) then
         call field_dump_info1D(tracers_x,      tracers_number,'X'     )
         call field_dump_info1D(tracers_xprev,  tracers_number,'Xprev' )
         call field_dump_info1D(tracers_Vx,     tracers_number,'Vx'    )
         call field_dump_info1D(tracers_Vxprev, tracers_number,'Vxprev')
         write(*,98) myid, 'tracers_compute: Done: ',tracers_t
         !write(*,99) myid, 'tracers_compute: mynx=',mynx
         !write(*,99) myid, 'tracers_compute: myny=',myny
      end if

 97    format('[',I3,'] ',A,4I5)
 98    format('[',I3,'] ',A,5F12.6)
 99    format('[',I3,'] ',A,2I5)

       return
       end     ! END: tracers_compute

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  tracers_out:                                                              cc
cc             Output tracer particles.                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine tracers_out(time)
       implicit     none
       real(kind=8) time
       include     'tracers.inc'
       include     'mem.inc'
       !
       integer     i, gft_rc, index, j, k
       integer     ntheta, nphi
       integer     gft_out_bbox, mystringlength
       external    gft_out_bbox, mystringlength
       real(kind=8) dtheta, theta, integral1, integral2
       real(kind=8) dphi,   phi, sinphidtdp
       real(kind=8) bbox(2)
       !data         bbox / 0.d0, 0.d0 /
       character(20) fname
       logical     ltrace
       parameter ( ltrace  = .false. )
       logical     ltrace2
       parameter ( ltrace2 = .false. )

       if (ltrace) write(*,*) 'tracers_out: Outputting surfaces'

#ifdef AIX
!      do i = 1, num_afuncs
!         if (asf_out(i) .eq. 0 .or. asf_out(i) .eq.2) then
!            if (ltrace2)write(*,*)'tracers_out: Out fct:',i,asf_name(i)

!            gft_rc = gft_out_bbox(  trim(asf_name(i))/ /CHAR(0),
!    *                           %VAL(time),   asf_shape,
!    *                           %VAL(2), asf_bbox, q(asf_pointer(i)) )
!         end if

!         if (asf_out(i) .eq.1 .or. asf_out(i) .eq.2) then
!            if (ltrace2)write(*,*)'tracers_out: Out Int:',i,asf_name(i)
!         end if
!      end do
#else
             bbox(1) = 1
             bbox(2) = 1.d0*tracers_number
             gft_rc = gft_out_bbox(  'trac_x', time,   tracers_number,
     *                               1, bbox, tracers_x )
             gft_rc = gft_out_bbox(  'trac_y', time,   tracers_number,
     *                               1, bbox, tracers_y )
             gft_rc = gft_out_bbox(  'trac_z', time,   tracers_number,
     *                               1, bbox, tracers_z )
             gft_rc = gft_out_bbox(  'trac_vx', time,   tracers_number,
     *                               1, bbox, tracers_Vx )
             gft_rc = gft_out_bbox(  'trac_vy', time,   tracers_number,
     *                               1, bbox, tracers_Vy )
             gft_rc = gft_out_bbox(  'trac_vz', time,   tracers_number,
     *                               1, bbox, tracers_Vz )
#endif

       if (ltrace) write(*,*) 'tracers_out: done outputing.'
       return
       end      ! end: tracers_out

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  test_interproutines(field,field2,tracers,nx,ny,nz,nt)
      implicit none
      integer  nx,ny,nz, nt
      real*8   field(nx,ny,nz), field2(nx,ny,nz), tracers(nt)
      !
      integer  i,j,k
      integer  ilow, jlow, klow
      real*8   fracx, fracy, fracz
      real*8   x, y, z, hg, hg2, inv_hg
      real*8   interp, interpderiv
      external interp, interpderiv
   
      ! first test: constant
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         field(i,j,k) = 5.d0
      end do
      end do
      end do

      ! Arbitrary:
      fracx = 0.1d0
      fracy = 0.0d0
      fracz = 0.0d0
      ilow  = (nx+1)/2.d0
      jlow  = (ny+1)/2.d0
      klow  = (nz+1)/2.d0
      do i = 1, nt
       tracers(i)  = interp(  field,
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
       ilow = ilow + (nx/2.0)/nt
      end do
      call field_dump_info1D(tracers, nt, 'const')
      do i = 1, nt
       tracers(i)  = interpderiv(field,
     *               fracx, fracy, fracz, 1.d0,
     *               ilow,jlow,klow,
     *               nx,ny,nz,1)
       ilow = ilow + (nx/2.0)/nt
      end do
      call field_dump_info1D(tracers, nt, 'zero')

      ! second test: ramp
      hg  = (1.0d0-0.00d0)/(nx-1)
      inv_hg = 1.d0/hg
      hg2 = (0.5d0-0.25d0)/(nx-1)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         field(i,j,k) = hg*sqrt((i-1.d0)**2+(j-1.d0)**2+(k-1.d0)**2)
      end do
      end do
      end do

      ! Arbitrary:
      fracx = 0.1d0
      fracy = 0.0d0
      fracz = 0.0d0
      ilow  = (nx+1)/2.d0
      jlow  = (ny+1)/2.d0
      klow  = (nz+1)/2.d0
      do i = 1, nt
       tracers(i)  = interp(  field,
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
       ilow = ilow + (nx/2.0)/nt
      end do
      call field_dump_info1D(tracers, nt, 'ramp')
      do i = 1, nt
       tracers(i)  = interpderiv(field,
     *               fracx, fracy, fracz, 1.d0,
     *               ilow,jlow,klow,
     *               nx,ny,nz,1)
       ilow = ilow + (nx/2.0)/nt
      end do
      hg  = (1.0d0-0.00d0)/(nx-1)
      hg2 = (0.5d0-0.25d0)/(nx-1)
      write(*,*) 'hg,hg2 = ',hg,hg2
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
       x = 0.25d0 + hg2*(i-1.0)
       y = 0.25d0 + hg2*(j-1.0)
       z = 0.25d0 + hg2*(k-1.0)
       ! Index of lower left point of interpolation cell:
       ilow  = INT( (x - 0.d0) / hg ) + 1
       jlow  = INT( (y - 0.d0) / hg ) + 1
       klow  = INT( (z - 0.d0) / hg ) + 1
       ! Fractional distance point is from lower left point:
       fracx = ( x - hg*(ilow-1) ) / hg
       fracy = ( y - hg*(jlow-1) ) / hg
       fracz = ( z - hg*(klow-1) ) / hg
       if (fracx.gt.1.or.fracy.gt.1.or.fracz.gt.1)then
          write(*,*) 'x/y/z:',x,y,z
          write(*,*) 'i/j/klow:',ilow,jlow,klow
          write(*,*) 'diff:',x,hg*(ilow-1),x - hg*(ilow-1)
          write(*,*) 'fracx/y/z:',fracx,fracy,fracz
       end if
       !
       field2(i,j,k)= interpderiv(field,
     *               fracx, fracy, fracz, inv_hg,
     *               ilow,jlow,klow,
     *               nx,ny,nz,3)
!      field2(i,j,k)= interp(field,
!    *               fracx, fracy, fracz,
!    *               ilow,jlow,klow,
!    *               nx,ny,nz)
       !field2(i,j,k)= 3.d0*(i-1.0)
      end do
      end do
      end do
      call field_dump_info1D(tracers, nt, 'constant')
      call field_out1d(tracers,0.d0,'slice',0.d0,1.d0,nt,0)
      call field_out3d(field, 0.d0,
     *             'ramp',
     *             0.d0,1.d0,
     *             0.d0,1.d0,
     *             0.d0,1.d0,
     *             nx,ny,nz,0)
      call field_dump_infoB(field2,nx,ny, nz,'field2')
      call field_out3d(field2, 0.d0,
     *             'rampderiv',
     *             0.25d0,0.5d0,
     *             0.25d0,0.5d0,
     *             0.25d0,0.5d0,
     *             nx,ny,nz,0)

      return
      end ! test_interproutines
