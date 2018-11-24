cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_intp_bnd_master:                                                     cc
cc                         Interpolate in time from parent grid               cc
cc                     and pack and send the data to                          cc
cc                     the process which owns the grid, which                 cc
cc                     should be running grid_intp_bnd_slave.                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_intp_bnd_master(gridnum,parent)
      implicit    none
      integer     gridnum, parent
      include     'grid.inc'
      include     'glob.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     nxp, nyp, nzp, nxc, nyc, nzc,
     *            i, length, owner, gwc,
     *            levp, levc, countp, countc, grid_return_level,
     *            nxs, nys, nzs,
     *            bi, ei, bj, ej, bk, ek
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      real(kind=8)      tp_np1,  tp_n, tc_np1, tmp
      logical     grid_return_existence
      character(3) tmps,tmps2
      !
      ! For FV we need our own variables since we do not
      ! want to use an expanded region even if user calls for cubic interp
      !
      integer     min_icFV,max_icFV, min_jcFV,max_jcFV,min_kcFV,max_kcFV
      integer     min_ipFV,max_ipFV, min_jpFV,max_jpFV,min_kpFV,max_kpFV
      integer     nxsFV, nysFV, nzsFV, lengthFV
      integer     c_bnds_C(6), p_bnds_C(6), shp_C(3), length_C, b_C(3)
      integer     c_bnds_L(6), p_bnds_L(6), shp_L(3), length_L, b_L(3)
      integer     c_bnds_Q(6), p_bnds_Q(6), shp_Q(3), length_Q, b_Q(3)

      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 1 )
      logical     ltrace
      parameter ( ltrace = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )
      ! Determine whether it is safe to reset bounds of vfuncs:
      logical     handlevfuncs
      integer     ADV
      
      owner  = gr_own(gridnum)
      
      !ltrace  = .false.  
      !ltrace2 = .false. 
      !!if (gr_level(gridnum).eq.2) then
      !if (gridnum .eq. 19) then
      !if (gr_level(gridnum).eq.4.and. owner.eq.3) then
      !   ltrace  = .true.
      !   ltrace2 = .true.
      !end if


      levp   = gr_level(parent)
      levc   = gr_level(gridnum)

      countp = lev_count(levp)
      countc = lev_count(levc)

      call load_pointers(parent)

      nxp    = gr_nx(parent)
      nyp    = gr_ny(parent)
      nzp    = gr_nz(parent)

      nxc    = gr_nx(gridnum)
      nyc    = gr_ny(gridnum)
      nzc    = gr_nz(gridnum)

      call grid_find_intersection(gridnum, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, gridnum, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      ! Keep these for FV fields:
      min_ipFV = min_ip
      min_jpFV = min_jp
      min_kpFV = min_kp
      max_ipFV = max_ip
      max_jpFV = max_jp
      max_kpFV = max_kp
      !
      min_icFV = min_ic
      min_jcFV = min_jc
      min_kcFV = min_kc
      max_icFV = max_ic
      max_jcFV = max_jc
      max_kcFV = max_kc
      !
      nxsFV    = max_ipFV - min_ipFV + 1
      nysFV    = max_jpFV - min_jpFV + 1
      nzsFV    = max_kpFV - min_kpFV + 1
      !
      lengthFV = nxsFV * nysFV * nzsFV

      gwc  = ghostwidth /refine_factor + 1

      !
      ! For higher order methods, in each direction either:
      !     (1) expand the region of the parent you send to the child
      !  or
      !     (2) shrink the region on the child in which you interpolate
      !         ...this shouldn't worry you because in this case some
      !         other parent must overlap in that region anyway.
      !
      !
      !
      call set_interp_bounds(  min_ip, max_ip, min_jp, max_jp,
     &                         min_kp, max_kp, min_ic, max_ic,
     &                         min_jc, max_jc, min_kc, max_kc,
     &                         nxp, nyp, nzp, length,
     &                         bi, bj, bk, refine_factor,
     &                         p_bnds_L, c_bnds_L, shp_L, length_L, b_L,
     &                         p_bnds_C, c_bnds_C, shp_C, length_C, b_C,
     &                         p_bnds_Q, c_bnds_Q, shp_Q, length_Q, b_Q)

      if (c_bnds_Q(2) .lt. c_bnds_Q(1)) return
      if (c_bnds_Q(4) .lt. c_bnds_Q(3)) return
      if (c_bnds_Q(6) .lt. c_bnds_Q(5)) return


      ! Notes about these times:
      !     parent grid advanced already w/r/t child grid
      !                hence the subtraction to find tp_n
      !     child grid not yet advanced
      !                hence the addition to find tc_np1
      !
      !         ____________________   parent grid time: tp_np1
      !              .......
      !              .......           child grid time:  tc_np1
      !              .......
      !         ____________________   parent grid time: tp_n
      !

      !
      ! Use level counters instead:
      !    NB: doesn't actually compute the time to which the
      !        grid is evolved. Instead, just computes relative
      !        time differences which is all that's needed.
      !
      tp_n   = 0.d0
      tp_np1 = lambda*gr_h(parent)  
      tmp    = mod(countc,refine_factor) + 1
      tc_np1 = tmp * tp_np1 / refine_factor

      !
      ! People want the quantities computed in analysis
      ! to look "right"...so we can set their boundaries
      ! only if we do so when things are time aligned.
      ! Such quantities do not have two time levels on the
      ! parent so we have to be careful about getting the
      ! data from the right field.
      !
      handlevfuncs = .false.
      if (amrbound_prepost.eq.1.and.amrbound_timealign.eq.0) then
         handlevfuncs = .true.
      end if
      !  Force not to do the vfuncs
      !handlevfuncs = .false.


      !
      ! Checks in routine: int_field_bounds_time()
      ! with accurate length, so bypassing here:
      !
      if (length_C.gt.nxp*nyp*nzp) then
         write(*,*) 'grid_intp_bnd_master: gridnum = ', gridnum
         write(*,*) 'grid_intp_bnd_master: gwc     = ', gwc
         write(*,*) 'grid_intp_bnd_master: nx/y/zc = ', nxc,nyc,nzc
         write(*,*) 'grid_intp_bnd_master: nx/y/zp = ', nxp,nyp,nzp
         write(*,*) 'grid_intp_bnd_master: countc  = ', countc
         write(*,*) 'grid_intp_bnd_master: countp  = ', countp
         write(*,*) 'grid_intp_bnd_master: owner   = ', owner
         write(*,*) 'grid_intp_bnd_master: parent  = ', parent
         write(*,*) 'grid_intp_bnd_master: powner  = ', gr_own(parent)
         write(*,*) 'grid_intp_bnd_master: length  = ', length
         write(*,*) 'grid_intp_bnd_master:nxp*nyp*nzp',nxp*nyp*nzp
        write(*,*)'grid_intp_bnd_master:min_i/j/kp',min_ip,min_jp,min_kp
        write(*,*)'grid_intp_bnd_master:max_i/j/kp',max_ip,max_jp,max_kp
         write(*,*) 'Problem of space in grid_intp_bnd_master'
      end if
      if (ltrace) then
        call int2str(gridnum,tmps)
        call int2str(parent,tmps2)
         write(*,99)myid,'grid_intp_bnd_master: Outto: '//tmps//tmps2
         write(*,99)myid,'grid_intp_bnd_master: gridnum = ', gridnum
         write(*,99)myid,'grid_intp_bnd_master: owner   = ', owner
         write(*,99)myid,'grid_intp_bnd_master: parent  = ', parent
         write(*,99)myid,'grid_intp_bnd_master: powner =',gr_own(parent)
         write(*,99)myid,'grid_intp_bnd_master: length  = ', length
         write(*,98)myid,'grid_intp_bnd_master: time    = ',gr_t(parent)
         write(*,99)myid,'grid_intp_bnd_master: gwc     = ', gwc
         write(*,99)myid,'grid_intp_bnd_master: nx/y/zc = ', nxc,nyc,nzc
         write(*,99)myid,'grid_intp_bnd_master: nx/y/zp = ', nxp,nyp,nzp
         write(*,99)myid,'grid_intp_bnd_master: countc  = ', countc
         write(*,99)myid,'grid_intp_bnd_master: countp  = ', countp
         write(*,99)myid,'grid_intp_bnd_master: owner   = ', owner
         write(*,99)myid,'grid_intp_bnd_master: powner =',gr_own(parent)
         write(*,99)myid,'grid_intp_bnd_master: length  = ', length
         write(*,99)myid,'grid_intp_bnd_master:nxp*nyp*nzp',nxp*nyp*nzp
        write(*,99)myid,'grid_intp_bnd_master:min_i/j/kp',min_ip,min_jp,
     *                              min_kp
        write(*,99)myid,'grid_intp_bnd_master:max_i/j/kp',max_ip,max_jp,
     *                              max_kp
      end if

      !
      ! Assumes that advanced gfunc is next after integral one:
      !
      ! Notes about these pointers:
      !         the parent has already been "STEP"ed and hence
      !         the advanced data was SWAPed into the _n spot
      !             while the other data is now stored in the _np1 spot
      !         this is to explain the order of the two parent pointers
      !             in the call to int_field_bounds below.
      !
      do i = 1, num_gfuncs-1
         ! Set the variables that define the interpolation parameters
         ! for each type of interpolation.
         if (gfunc_interptype(i) .eq. GFUNC_INTERP_LINEAR) then
           call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxs, nys, nzs, length,
     &                              bi, bj, bk,
     &                              p_bnds_L, c_bnds_L, shp_L,
     &                              length_L, b_L)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_CUBIC) then
           call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxs, nys, nzs, length,
     &                              bi, bj, bk,
     &                              p_bnds_C, c_bnds_C, shp_C,
     &                              length_C, b_C)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_QUINTIC .or.
     &            gfunc_interptype(i) .eq. GFUNC_INTERP_WENO) then
           call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxs, nys, nzs, length,
     &                              bi, bj, bk,
     &                              p_bnds_Q, c_bnds_Q, shp_Q,
     &                              length_Q, b_Q)
         end if

         if (  gfunc_type(i) .eq. GFUNC_INTEGRAL .or.
     *        (gfunc_type(i) .eq. GFUNC_INTEGRAL_DER .and. handlevfuncs)
     *       ) then
            call int_bounds_time(q(gr_tmp),
     *                           q(gfunc_pointer(i,parent)),
     *                           q(gfunc_pointer(i+1,parent)),
     *                           tc_np1, tp_np1, tp_n,
     *                           min_ip, min_jp, min_kp,
     *                           max_ip, max_jp, max_kp,
     *                           nxs,nys,nzs,
     *                           nxp,nyp,nzp )
            if (i.eq.outfnum.and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,parent)),gr_t(parent),
     *             'parmasn'//tmps//tmps2,
     *             gr_minx(parent),gr_maxx(parent),
     *             gr_miny(parent),gr_maxy(parent),
     *             gr_minz(parent),gr_maxz(parent),
     *             nxp,nyp,nzp,myid)
!           if (i.eq.1.and. ltrace)
!    *      call field_out3d(q(gfunc_pointer(i+1,parent)),gr_t(parent),
!    *             'parmas'//tmps//tmps2,
!    *             gr_minx(parent),gr_maxx(parent),
!    *             gr_miny(parent),gr_maxy(parent),
!    *             gr_minz(parent),gr_maxz(parent),
!    *             nxp,nyp,nzp,myid)
!           if (i.eq.1.and. ltrace)
!    *      call field_out3d(q(gr_tmp),gr_t(parent),
!    *                     'tmpmas'//tmps//tmps2,
!    *       q(gr_x(parent)+(min_ip-1)),q(gr_x(parent)+(max_ip-1)),
!    *       q(gr_y(parent)+(min_jp-1)),q(gr_y(parent)+(max_jp-1)),
!    *       q(gr_z(parent)+(min_kp-1)),q(gr_z(parent)+(max_kp-1)),
!    *             nxs,nys,nzs,myid)
            if (i.eq.outfnum .and. ltrace) then
               write(*,98)myid,'grid_intp_bnd_master: xmin/max:',
     *     q(gr_x(parent)),q(gr_x(parent+nxp-1))
               write(*,98)myid,'grid_intp_bnd_master: ymin/max:',
     *     q(gr_y(parent)),q(gr_y(parent+nyp-1))
               write(*,98)myid,'grid_intp_bnd_master: zmin/max:',
     *     q(gr_z(parent)),q(gr_z(parent+nzp-1))
               write(*,98)myid,'grid_intp_bnd_master: x:',
     *     q(gr_x(parent)+min_ip-1),q(gr_x(parent)+max_ip-1)
               write(*,98)myid,'grid_intp_bnd_master: y:',
     *     q(gr_y(parent)+min_jp-1),q(gr_y(parent)+max_jp-1)
               write(*,98)myid,'grid_intp_bnd_master: z:',
     *     q(gr_z(parent)+min_kp-1),q(gr_z(parent)+max_kp-1)
            end if
            if (ltrace2) write(*,99)myid,'    Sending data for ',i
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_BND_DATA, MPI_COMM_WORLD, ierr)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV ) then
            !
            ! Finite Volume fields
            !
            ! 1) Interpolate in time on parent       --> gr_tmp
            ! 2) Initialize whole region for FV      --> gr_flag
            ! 3) Paste points into boundary of field --> q()
            call int_bounds_time(q(gr_tmp),
     *                              q(gfunc_pointer(i,parent)),
     *                              q(gfunc_pointer(i+1,parent)),
     *                              tc_np1, tp_np1, tp_n,
     *                              min_ipFV, min_jpFV, min_kpFV,
     *                              max_ipFV, max_jpFV, max_kpFV,
     *                              nxsFV,nysFV,nzsFV,
     *                              nxp,nyp,nzp )
            call MPI_Send(q(gr_tmp), lengthFV, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_BND_DATA, MPI_COMM_WORLD, ierr)
         end if
      end do
 
      if (ltrace) then
         write(*,99)myid,'grid_intp_bnd_master: Done.***'
      end if

 97   format('[',I3,'] ',A,I5,3F10.5)
 98   format('[',I3,'] ',A,3F10.5)
 99   format('[',I3,'] ',A,3I5)

      return
      end    ! END: grid_intp_bnd_master

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_intp_bnd:  For a given grid, we interpolate from its                 cc
cc                  parent values for the "inner" boundaries.                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_intp_bnd(gridnum,parent)
      implicit    none
      integer     gridnum, parent
      include     'grid.inc'
      include     'glob.inc'
      include     'grid_methods.inc'
      integer     nxp, nyp, nzp, nxc, nyc, nzc,
     *            ax,  ay,  az,  i, length,
     *            countc, gwc, nxs, nys, nzs,
     *            bi, ei, bj, ej, bk, ek
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      real(kind=8)      tp_np1,  tp_n, tc_np1, tmp
      integer     myid, proc_return_myid
      character(3) tmps,tmps2
      real(kind=8)      myl2norm3d, errornorm
      external          myl2norm3d
      !
      ! For FV we need our own variables since we do not
      ! want to use an expanded region even if user calls for cubic interp
      !
      integer     min_icFV,max_icFV, min_jcFV,max_jcFV,min_kcFV,max_kcFV
      integer     min_ipFV,max_ipFV, min_jpFV,max_jpFV,min_kpFV,max_kpFV
      integer     nxsFV, nysFV, nzsFV, lengthFV
      integer     p_bnds_L(6), c_bnds_L(6), shp_L(3), length_L, b_L(3)
      integer     p_bnds_C(6), c_bnds_C(6), shp_C(3), length_C, b_C(3)
      integer     p_bnds_Q(6), c_bnds_Q(6), shp_Q(3), length_Q, b_Q(3)

      ! Determine whether it is safe to reset bounds of vfuncs:
      logical     handlevfuncs
      integer     ADV
      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 1 )

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      myid = proc_return_myid()

      !ltrace  = .false.  
      !ltrace2 = .false. 
      !if (gr_level(gridnum).eq.4.and. myid.eq.3.and.gridnum.eq.15) then
      !if (gridnum.eq.179)then
         !ltrace  = .true.
         !ltrace2 = .true.
      !end if


      if (ltrace) write(*,*) 'grid_intp_bnd: gridnum = ', gridnum

      call grid_find_intersection(gridnum, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, gridnum, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      ! Keep these for FV fields:
      min_ipFV = min_ip
      min_jpFV = min_jp
      min_kpFV = min_kp
      max_ipFV = max_ip
      max_jpFV = max_jp
      max_kpFV = max_kp
      !
      min_icFV = min_ic
      min_jcFV = min_jc
      min_kcFV = min_kc
      max_icFV = max_ic
      max_jcFV = max_jc
      max_kcFV = max_kc
      !
      nxsFV    = max_ipFV - min_ipFV + 1
      nysFV    = max_jpFV - min_jpFV + 1
      nzsFV    = max_kpFV - min_kpFV + 1
      !
      lengthFV = nxsFV * nysFV * nzsFV

      nxp    = gr_nx(parent)
      nyp    = gr_ny(parent)
      nzp    = gr_nz(parent)

      nxc    = gr_nx(gridnum)
      nyc    = gr_ny(gridnum)
      nzc    = gr_nz(gridnum)

      !
      ! How many coarse grid points do we need?
      !  If r=2:
      !     gw=1 --> gwc=1
      !     gw=2 --> gwc=2
      !     gw=3 --> gwc=2
      !     gw=4 --> gwc=3
      !     gw=5 --> gwc=3
      !     gw=6 --> gwc=4
      !
      gwc  = ghostwidth /refine_factor + 1

      !
      ! For higher order methods, in each direction either:
      !     (1) expand the region of the parent you send to the child
      !  or
      !     (2) shrink the region on the child in which you interpolate
      !         ...this shouldn't worry you because in this case some
      !         other parent must overlap in that region anyway.
      !
      !
      !
      call set_interp_bounds( min_ip, max_ip, min_jp, max_jp,
     &                        min_kp, max_kp, min_ic, max_ic,
     &                        min_jc, max_jc, min_kc, max_kc,
     &                        nxp, nyp, nzp, length,
     &                        bi, bj, bk, refine_factor,
     &                        p_bnds_L, c_bnds_L, shp_L, length_L, b_L,
     &                        p_bnds_C, c_bnds_C, shp_C, length_C, b_C,
     &                        p_bnds_Q, c_bnds_Q, shp_Q, length_Q, b_Q)


      if (c_bnds_Q(2) .lt. c_bnds_Q(1)) return
      if (c_bnds_Q(4) .lt. c_bnds_Q(3)) return
      if (c_bnds_Q(6) .lt. c_bnds_Q(5)) return

      countc = lev_count(grid_return_level(gridnum))


      !
      ! Checks in routine: int_field_bounds_time()
      ! with accurate length, so bypassing here:
      !
      if (length_C.gt.nxp*nyp*nzp) then
         write(*,*) 'Problem of space in grid_intp_bnd'
         write(*,*) 'length = ',length
         write(*,*) 'nxp*nyp*nzp = ',nxp*nyp*nzp
         call my_exit('Not enough space in grid_intp_bnd')
      end if

      !
      ! Notes about these times:
      !
      !         ____________________   parent grid time: tp_np1
      !              .......
      !              .......           child grid time:  tc_np1
      !              .......
      !         ____________________   parent grid time: tp_n
      !
      !
      ! Use level counters:
      !
      tp_n   = 0.d0
      tp_np1 = lambda*gr_h(parent)  
      tmp    = mod(countc,refine_factor) + 1
      tc_np1 = tmp * tp_np1 / refine_factor

      !
      ! People want the quantities computed in analysis
      ! to look "right"...so we can set their boundaries
      ! only if we do so when things are time aligned.
      ! Such quantities do not have two time levels on the
      ! parent so we have to be careful about getting the
      ! data from the right field.
      !
      handlevfuncs = .false.
      if (amrbound_prepost.eq.1.and.amrbound_timealign.eq.0) then
         handlevfuncs = .true.
      end if
      !  Force not to do the vfuncs
      !handlevfuncs = .false.

      if (ltrace) then
        call int2str(gridnum,tmps)
        call int2str(parent,tmps2)
         write(*,99)myid,'grid_intp_bnd:  Doing bounds on grid:',gridnum
         write(*,99)myid,'grid_intp_bnd: gridnum = ', gridnum
         write(*,99)myid,'grid_intp_bnd: parent  = ', parent
         write(*,99)myid,'grid_intp_bnd: nx/y/zc = ', nxc,nyc,nzc
         write(*,99)myid,'grid_intp_bnd: nx/y/zp = ', nxp,nyp,nzp
         !write(*,99)myid,'grid_intp_bnd: nx/y/zs = ', nxs,nys,nzs
         write(*,99)myid,'grid_intp_bnd:mini/j/kc ',min_ic,min_jc,min_kc
         write(*,99)myid,'grid_intp_bnd:maxi/j/kc ',max_ic,max_jc,max_kc
         write(*,99)myid,'grid_intp_bnd:mini/j/kp ',min_ip,min_jp,min_kp
         write(*,99)myid,'grid_intp_bnd:maxi/j/kp ',max_ip,max_jp,max_kp
         write(*,99)myid,'grid_intp_bnd: gw      = ', ghostwidth
         write(*,99)myid,'grid_intp_bnd: gwc     = ', gwc
         write(*,99)myid,'grid_intp_bnd: tmp     = ', tmp
         write(*,99)myid,'grid_intp_bnd: length  = ', length
         write(*,99)myid,'grid_intp_bnd: nxp*nyp*nzp',nxp*nyp*nzp
         write(*,99)myid,'grid_intp_bnd: amrbound_prepost   = ',
     *                                           amrbound_prepost
         write(*,99)myid,'grid_intp_bnd: amrbound_timealign = ',
     *                                           amrbound_timealign
            write(*,98)myid,'grid_intp_bnd: time    = ', gr_t(gridnum)
         write(*,98)myid,'grid_intp_bnd: tp_n/tp_np1 =',tp_n,tp_np1
         write(*,98)myid,'grid_intp_bnd: tc_np1      =',tc_np1
         write(*,98)myid,'grid_intp_bnd: alpha_n/np1 =',
     *                (tp_np1-tc_np1)/(tp_np1-tp_n),
     *               1.d0-(tp_np1-tc_np1)/(tp_np1-tp_n)
         if (ltrace2) then
            ADV=0
            call field_out3d(q(gfunc_pointer(outfnum+ADV,gridnum)),
     *             gr_t(gridnum),'pret'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
            ADV=1
            call field_out3d(q(gfunc_pointer(outfnum+ADV,gridnum)),
     *             gr_t(gridnum),'pretADV'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
         end if
      end if

      if (amrbound_prepost.eq.1) call grid_swap_n_np1(gridnum)
      call load_pointers(gridnum)

      !
      ! Assumes that advanced gfunc is next after integral one:
      !
      ! Notes about these pointers:
      !         the parent has already been "STEP"ed and hence
      !         the advanced data was SWAPed into the _n spot
      !             while the other data is now stored in the _np1 spot
      !         this is to explain the order of the two parent pointers
      !             in the call to int_field_bounds below.
      !
      do i = 1, num_gfuncs-1
         if (gfunc_interptype(i) .eq. GFUNC_INTERP_LINEAR) then
            call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxs, nys, nzs, length,
     &                              bi, bj, bk,
     &                              p_bnds_L, c_bnds_L, shp_L, length_L,
     &                              b_L)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_CUBIC) then
            call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxs, nys, nzs, length,
     &                              bi, bj, bk,
     &                              p_bnds_C, c_bnds_C, shp_C, length_C,
     &                              b_C)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_QUINTIC .or. 
     &            gfunc_interptype(i) .eq. GFUNC_INTERP_WENO) then
            call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxs, nys, nzs, length,
     &                              bi, bj, bk,
     &                              p_bnds_Q, c_bnds_Q, shp_Q, length_Q,
     &                              b_Q)
         else
            write(*,99)myid,'grid_intp_bnd: Unknown INTERPTYPE'
            call my_exit('Problem in grid_intp_bnd')
         end if
         !
         if (  gfunc_type(i) .eq. GFUNC_INTEGRAL .or.
     *        (gfunc_type(i) .eq. GFUNC_INTEGRAL_DER .and. handlevfuncs)
     *       ) then
            !
            if (gfunc_type(i).eq.GFUNC_INTEGRAL_DER) then
               ADV = 0
            else
               !
               ! For the u-funcs, put data into the advanced field:
               !
               ADV = 1
            end if
            if(ltrace2) write(*,99)myid,'    ADV = ', ADV
            !
            !    1) Copy or interpolate in time whole of overlap
            !       from parent.
            !    2) Interpolate in space whole of overlap into work space
            !    3) Paste data from boundaries into field
            !
            call int_bounds_time(q(gr_tmp),
     *                           q(gfunc_pointer(i,parent)),
     *                           q(gfunc_pointer(i+1,parent)),
     *                           tc_np1, tp_np1, tp_n,
     *                           min_ip, min_jp, min_kp,
     *                           max_ip, max_jp, max_kp,
     *                           nxs,nys,nzs,
     *                           nxp,nyp,nzp )
            !
            ! Spatially interpolate everywhere into temp space (gr_flag):
            !call load_scal1D(q(gr_tmp),1.d0*gridnum,nxs*nys*nzs)
            !call load_scal_mult3D(q(gr_tmp),q(gr_tmp),99.d0,nxs,nys,nzs)
            if (i.eq.outfnum .and. ltrace) then
         write(*,99)myid,'grid_intp_bnd: nx/y/zc = ', nxc,nyc,nzc
         write(*,99)myid,'grid_intp_bnd: nx/y/zp = ', nxp,nyp,nzp
      write(*,99)myid,'grid_intp_bnd: nx/y/zs = ', nxs,nys,nzs
      write(*,99)myid,'grid_intp_bnd:mini/j/kc',min_ic,min_jc,min_kc
      write(*,99)myid,'grid_intp_bnd:maxi/j/kc',max_ic,max_jc,max_kc
      write(*,99)myid,'grid_intp_bnd:mini/j/kp',min_ip,min_jp,min_kp
      write(*,99)myid,'grid_intp_bnd:maxi/j/kp',max_ip,max_jp,max_kp
               write(*,98)myid,'grid_intp_bnd: x:',
     *     q(gr_x(gridnum)+min_ic-1),q(gr_x(gridnum)+max_ic-1)
               write(*,98)myid,'grid_intp_bnd: y:',
     *     q(gr_y(gridnum)+min_jc-1),q(gr_y(gridnum)+max_jc-1)
               write(*,98)myid,'grid_intp_bnd: z:',
     *     q(gr_z(gridnum)+min_kc-1),q(gr_z(gridnum)+max_kc-1)
            end if
            if (i.eq.outfnum.and. ltrace2)
     *      call field_out3d(q(gfunc_pointer(i,parent)),gr_t(gridnum),
     *             'parn'//tmps//tmps2,
     *             gr_minx(parent),gr_maxx(parent),
     *             gr_miny(parent),gr_maxy(parent),
     *             gr_minz(parent),gr_maxz(parent),
     *             nxp,nyp,nzp,myid)
!           if (i.eq.1.and. ltrace)
!    *      call field_out3d(q(gfunc_pointer(i+1,parent)),gr_t(gridnum),
!    *             'parnp1'//tmps//tmps2,
!    *             gr_minx(parent),gr_maxx(parent),
!    *             gr_miny(parent),gr_maxy(parent),
!    *             gr_minz(parent),gr_maxz(parent),
!    *             nxp,nyp,nzp,myid)
            if (i.eq.outfnum.and. ltrace2)
     *      call field_out3d(q(gr_tmp),gr_t(gridnum),
     *                     'tmp'//tmps//tmps2,
     *       q(gr_x(parent)+(min_ip-1)),q(gr_x(parent)+(max_ip-1)),
     *       q(gr_y(parent)+(min_jp-1)),q(gr_y(parent)+(max_jp-1)),
     *       q(gr_z(parent)+(min_kp-1)),q(gr_z(parent)+(max_kp-1)),
     *             nxs,nys,nzs,myid)
            if (i.eq.outfnum .and. ltrace) then
         write(*,99)myid,'grid_intp_bnd: nx/y/zc = ', nxc,nyc,nzc
         write(*,99)myid,'grid_intp_bnd: nx/y/zp = ', nxp,nyp,nzp
      write(*,99)myid,'grid_intp_bnd: nx/y/zs = ', nxs,nys,nzs
      write(*,99)myid,'grid_intp_bnd:mini/j/kc',min_ic,min_jc,min_kc
      write(*,99)myid,'grid_intp_bnd:maxi/j/kc',max_ic,max_jc,max_kc
      write(*,99)myid,'grid_intp_bnd:mini/j/kp',min_ip,min_jp,min_kp
      write(*,99)myid,'grid_intp_bnd:maxi/j/kp',max_ip,max_jp,max_kp
               write(*,98)myid,'grid_intp_bnd: x:',
     *     q(gr_x(gridnum)+min_ic-1),q(gr_x(gridnum)+max_ic-1)
               write(*,98)myid,'grid_intp_bnd: y:',
     *     q(gr_y(gridnum)+min_jc-1),q(gr_y(gridnum)+max_jc-1)
               write(*,98)myid,'grid_intp_bnd: z:',
     *     q(gr_z(gridnum)+min_kc-1),q(gr_z(gridnum)+max_kc-1)
            end if
            if (gfunc_interptype(i) .eq. GFUNC_INTERP_LINEAR) then
               ! Linear interpolation:
               call init_fieldGWa( q(gr_tmp),
     *                             q(gr_flag),
     *                             nxs,nys,nzs,
     *                             min_ic,min_jc,min_kc,
     *                             max_ic,max_jc,max_kc,
     *                             nxc,nyc,nzc,
     *                             refine_factor,ghostwidth)
            else if (gfunc_interptype(i) .eq. GFUNC_INTERP_CUBIC) then
               ! Cubic interpolation:
               call interp_from_parentB( q(gr_tmp),
     *                          q(gr_flag),
     *                          q(gr_chr), q(gr_mask),
     *                          nxs,nys,nzs,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            else if (gfunc_interptype(i) .eq. GFUNC_INTERP_QUINTIC) then
               ! 5th order interpolation:
               call interp_from_parentC( q(gr_tmp),
     *                          q(gr_flag),
     *                          q(gr_chr), q(gr_mask),
     *                          nxs,nys,nzs,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            else if (gfunc_interptype(i) .eq. GFUNC_INTERP_WENO) then
               ! WENO interpolation:
!              call interp_from_parentW( q(gr_tmp),
!    *                          q(gr_flag),
!    *                          q(gr_chr), q(gr_mask),
!    *                          nxs,nys,nzs,
!    *                          min_ic,min_jc,min_kc,
!    *                          max_ic,max_jc,max_kc,
!    *                          bi,    bj,    bk,
!    *                          nxc,nyc,nzc, refine_factor,ghostwidth)
 
            else
               write(*,*)'Unknown interpolation type'
            end if
            if ( shadow .ne. 0 .and. ADV.eq.1) then
               !
               ! Subtract parent data from data on this grid
               ! to affect the self-shadow hierarchy
               !
               call tre_subtract(     q(gr_flag),
     *                                q(gfunc_pointer(i+1,gridnum)),
     *                                q(gr_error),
     *                                q(gr_chr),
     *                                gfunc_shadow(i),
     *                                min_ic, min_jc, min_kc,
     *                                max_ic, max_jc, max_kc,
     *                                nxc,nyc,nzc, num_injfields,
     *                                refine_factor, ghostwidth)
               if(ltrace)write(*,97)myid,'grid_intp_bnd:compute error:',
     *             i,myl2norm3d(q(gr_error),nxc,nyc,nzc)
            end if
            !
            ! Paste data into boundary regions:
            !
            if (i.eq.outfnum.and. ltrace2)
     *      call field_out3d(q(gr_flag),gr_t(gridnum),
     *                     'flag'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
            if(ltrace2) write(*,99)myid,'.   ADV = ', ADV
            if (i.eq.outfnum .and. ltrace2)
     *      call field_out3d(q(gfunc_pointer(i+ADV,gridnum)),
     *             gr_t(gridnum),'prei'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
            if (i.eq.outfnum .and. ltrace2)
     *      call field_out3d(q(gfunc_pointer(gr_chr,gridnum)),
     *             gr_t(gridnum),'intpchr'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
            call paste_boundaries(    q(gr_flag),
     *                                q(gfunc_pointer(i+ADV,gridnum)),
     *                                q(gr_chr),
     *                                min_ic, min_jc, min_kc,
     *                                max_ic, max_jc, max_kc,
     *                                nxc,nyc,nzc,
     *                                refine_factor, ghostwidth)
            if (i.eq.outfnum .and. ltrace2)
     *      call field_out3d(q(gfunc_pointer(i,gridnum)),
     *             gr_t(gridnum),'post'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
            if (i.eq.outfnum .and. ltrace2)
     *      call field_out3d(q(gfunc_pointer(i+ADV,gridnum)),
     *             gr_t(gridnum),'postADV'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV ) then
            !
            ! Finite Volume fields
            !
            ! 1) Interpolate in time on parent       --> gr_tmp
            ! 2) Initialize whole region for FV      --> gr_flag
            ! 3) Paste points into boundary of field --> q()
            call int_bounds_time(q(gr_tmp),
     *                              q(gfunc_pointer(i,parent)),
     *                              q(gfunc_pointer(i+1,parent)),
     *                              tc_np1, tp_np1, tp_n,
     *                              min_ipFV, min_jpFV, min_kpFV,
     *                              max_ipFV, max_jpFV, max_kpFV,
     *                              nxsFV,nysFV,nzsFV,
     *                              nxp,nyp,nzp )
            call init_fieldFV( q(gr_tmp),
     *                            q(gr_flag),
     *                            nxsFV,nysFV,nzsFV,
     *                            min_icFV,min_jcFV,min_kcFV,
     *                            max_icFV,max_jcFV,max_kcFV,
     *                            nxc,nyc,nzc, refine_factor,ghostwidth)
            call paste_boundariesFV(q(gr_flag),
     *                                q(gfunc_pointer(i+1,gridnum)),
     *                                q(gr_chr),
     *                                min_icFV, min_jcFV, min_kcFV,
     *                                max_icFV, max_jcFV, max_kcFV,
     *                                nxc,nyc,nzc,
     *                                refine_factor, ghostwidth)
         end if
      end do
 
      if (shadow.gt.0) then
         if(ltrace)write(*,99)myid,'grid_intp_bnd: Increment chr array'
         call inc_mask_overlapshadow(q(gr_chr),
     *                           q(gr_flag),
     *                           q(gr_error),
     *                           min_ic, min_jc, min_kc,
     *                           max_ic, max_jc, max_kc,
     *                           nxc,nyc,nzc )
      end if

      if (shadow.eq.1) then
         errornorm = myl2norm3d(q(gr_error),nxc,nyc,nzc)
         if(ltrace)write(*,99)myid,'grid_intp_bnd: Errornorm:',errornorm
         if(ltrace)call field_dump_infob(q(gr_error),nxc,nyc,nzc,'erro')
         if(ltrace)  call field_out3d(q(gr_error),
     *             gr_t(gridnum),'intrperror'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
         if (errornorm.gt.1.d8) then
            write(*,99)myid,'grid_intp_bnd: The TRE is HUGE!'
            call field_dump_info(q(gr_error),nxc,nyc,nzc)
            write(*,99)myid,'grid_intp_bnd: gridnum = ', gridnum
            call grid_dump_info(gridnum)
            write(*,99)myid,'grid_intp_bnd: parent  = ', parent
            call grid_dump_info(parent)
            write(*,98)myid,'grid_intp_bnd: time    = ', gr_t(gridnum)
         end if
         if(ltrace2.or.errornorm.gt.1.d8)  call field_out3d(q(gr_error),
     *             gr_t(gridnum),'intrperror'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
      end if

      if (amrbound_prepost.eq.1) then
         if(ltrace)write(*,99)myid,'grid_intp_bnd: Re-swapping pointers'
         call grid_swap_n_np1(gridnum)
      end if

      !if(ltrace)call grid_test(gridnum)
      !write(*,99)myid,'grid_intp_bnd: call grid_test ',gridnum
      !call grid_test(gridnum)
      if(ltrace)write(*,99)myid,'grid_intp_bnd: Done.**********'

 97   format('[',I3,'] ',A,I5,3F10.5)
 98   format('[',I3,'] ',A,3F10.5)
 99   format('[',I3,'] ',A,3I5)


      return
      end    ! END: grid_intp_bnd

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_copy_n_np1:                                                          cc
cc                   Copy fields to advanced time level.                      cc
cc                   Used/called right after initialization so                cc
cc                   that both levels have something there.                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_copy_n_np1( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      integer     i

      !
      ! Assumes that advanced time level exists immediately
      ! following the level n data:                           
      !
      do i = 1, num_gfuncs - 1
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL ) then
            call mat_copy3d( q(gfunc_pointer(i,   gridnum)),
     *                       q(gfunc_pointer(i+1, gridnum)),
     *                       gr_nx(gridnum),
     *                       gr_ny(gridnum),
     *                       gr_nz(gridnum)                   )
         end if
      end do

      return
      end    ! END: grid_copy_n_np1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_h:                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function grid_return_h(gridnum)
      implicit     none
      integer      gridnum
      include     'grid.inc'

      grid_return_h = gr_h(gridnum)

      return
      end    ! END: grid_return_h

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_find_bounds:                                                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_find_bounds(gridnum,minx,maxx,miny,maxy,minz,maxz)
      implicit     none
      integer      gridnum
      real(kind=8) minx,maxx, miny,maxy, minz,maxz
      include     'grid.inc'
      integer      temp, i

      minx = gr_minx(gridnum)
      miny = gr_miny(gridnum)
      minz = gr_minz(gridnum)

      maxx = gr_maxx(gridnum)
      maxy = gr_maxy(gridnum)
      maxz = gr_maxz(gridnum)

      return
      end    ! END: grid_find_bounds

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_swap_n_np1                                                           cc
cc                  Requires that a given field, say "u", is followed         cc
cc                  immediately in order (in fields.inc) by its               cc
cc                  corresponding advanced field, "u_np1".                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_swap_n_np1( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      integer     temp, i

      do i = 1, num_gfuncs - 1
         !if ( gfunc_type(i) .eq. GFUNC_INTEGRAL ) then
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRALFV    ) then
            temp                       = gfunc_pointer(i,  gridnum)
            gfunc_pointer(i,  gridnum) = gfunc_pointer(i+1,gridnum)
            gfunc_pointer(i+1,gridnum) = temp
         end if
      end do

      ! Reload the pointers
      call load_pointers(gridnum)

      return
      end    ! END: grid_swap_n_np1


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_free:                                                                cc
cc               Grid no longer in use.                                       cc
cc               Therefore:    (1) add grid number to list of free grid #s    cc
cc                             (2) free memory (if owner)                     cc
cc                                                                            cc
cc               NB: updated for MPI. Freeing memory only done                cc
cc               if this is the owner of that grid...otherwise, there's       cc
cc               no memory allocated on this processor to be freed.           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_free( grid )
      implicit    none
      integer     grid
      include    'glob.inc'
      include    'grid.inc'
      include    'mpif.h'
      include    'mpi_stuff.inc'
      include    'grid_methods.inc'
      integer     size,    i

      logical     ltrace
      parameter ( ltrace = .false. )

 99   format('[',I3,'] ',A,4I4)

      if (ltrace) then
         write(*,99) myid, 'grid_free: Freeing grid ',grid
      end if

      if (grid.le.0) then
         write(*,*) 'grid_free: WARNING!!!!'
         write(*,*) 'grid_free: Grid Does Not exist: ',grid
         write(*,*) 'grid_free: Returning...'
         return
      end if

      !
      ! Reclaim free'ed grid number:
      !
      if (ltrace) write(*,99)myid,'grid_free:    numgrids = ',numgrids
      if (grid .eq. numgrids) then
         !
         ! Just decrement numgrids...
         ! and check if other grid numbers can be
         !
         numgrids = numgrids - 1
         if (ltrace) write(*,99)
     *      myid, 'grid_free: decrementing numgrids',numgrids
      else if (num_freegrids .lt. maxnumchildren) then
         !
         ! Otherwise keep track of free grid numbers
         !
         num_freegrids                    = num_freegrids + 1
         free_grid_numbers(num_freegrids) = grid
         if (num_freegrids.gt.max_num_freegrids)
     *      max_num_freegrids = num_freegrids
         if (ltrace) 
     *      write(*,99) myid, 'grid_free: adding grid to stored list'
      end if

      if (grid_return_owner(grid).eq.myid) then
         size = gr_nx(grid)*gr_ny(grid)*gr_nz(grid)
         !
         ! Explicitly free the memory:
         !    Note: most efficient if you mem_dealloc in order
         !          that grid_create() allocates them
         !
         call mem_dealloc( gr_x(grid),        gr_nx(grid) )
         call mem_dealloc( gr_y(grid),        gr_ny(grid) )
         call mem_dealloc( gr_z(grid),        gr_nz(grid) )
         do i = 1, num_gfuncs
            call mem_dealloc( gfunc_pointer(i,grid), size )
         end do
      end if

      if (ltrace) then
c        write(*,*) 'grid_free:   grid size     = ',size
         write(*,99) myid,'grid_free:   num_freegrids = ', num_freegrids
c        write(*,*) 'grid_free: Done.'
      end if

      return
      end     ! END: grid_free

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_step2:                                                               cc
cc              This routine evolved from grid_step() which is no longer used.cc
cc              After grid_iter() is called some number of times, this        cc
cc              routine completes the step process by swapping time levels    cc
cc              updating appropriate counters, post step smoothing, and       cc
cc              computing the derived quantities.                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_step2( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'glob.inc'
      include    'chr.inc'
      include    'mpif.h'
      include    'mpi_stuff.inc'
      integer    nx, ny, nz,       i,
     *           minx, miny, minz
      real(kind=8)     h,  dt
      real(kind=8)     myepsdis
      parameter     (  myepsdis = 0.1d0 )
      !real(kind=8)     myl2norm3d
      !external         myl2norm3d
      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 1 )
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      !ltrace  = .false.
      !ltrace2 = .false.
      !if (gridnum.eq.19) then
      !ltrace  = .true.
      !ltrace2 = .true.
      !end if

      if (ltrace) then
         write(*,99)myid,'grid_step2: gridnum = ',gridnum
         write(*,99)myid,'grid_step2: Calling grid_test before',gridnum
         call grid_test(gridnum)
      end if
      call load_pointers(gridnum)

      h    = gr_h(gridnum)
      dt   = lambda*h
      nx   = gr_nx(gridnum)
      ny   = gr_ny(gridnum)
      nz   = gr_nz(gridnum)
      minx = gr_minx(gridnum)
      miny = gr_miny(gridnum)
      minz = gr_minz(gridnum)

      if (ltrace2) then
            call field_out3d(q(gfunc_pointer(outfnum+1,gridnum)),
     *             gr_t(gridnum),'stepADV',
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nx,ny,nz,myid)
            call field_out3d(q(gfunc_pointer(outfnum,gridnum)),
     *             gr_t(gridnum),'step',
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nx,ny,nz,myid)
      end if
      if(ltrace)write(*,99)myid,'grid_step2: gridnum = ',gridnum
      if(ltrace2)call grid_dump_info(gridnum)

       !write(*,99)myid,'step2:Enter gridnum =',gridnum,gr_level(gridnum)
       !call field_dump_info(q(gr_psi4R),nx,ny,nz)

      !if (.false.) then
       ! Try some extra smoothing to control B-field growth in atmosphere:
       !call apply_diss(  q(gr_psi),    q(gr_tmp), myepsdis,    nx,ny,nz)
       !call apply_dissnh(q(gr_psi),    q(gr_tmp), myepsdis, 2, nx,ny,nz)
       !call apply_diss(  q(gr_bx_np1), q(gr_tmp), myepsdis,    nx,ny,nz)
       !call apply_diss(  q(gr_by_np1), q(gr_tmp), myepsdis,    nx,ny,nz)
       !call apply_diss(  q(gr_bz_np1), q(gr_tmp), myepsdis,    nx,ny,nz)
      !end if

c     if (.true.) then
      if (.false.) then
         if (gr_parent(gridnum).ge.0) then       ! if not coarse grid:
            if(ltrace)write(*,99)myid,'grid_step2: Smoothing'
            do i = 1, num_gfuncs
               if (gfunc_type(i) .eq. GFUNC_INTEGRAL_ADV ) then
c                 call smooth_near(q(gfunc_pointer(i,gridnum)),
c    *                             q(gr_chr),nx,ny,nz)
                  call smooth_nearGW(q(gfunc_pointer(i,gridnum)),
     *                             q(gr_chr),
     *                             nx,ny,nz,ghostwidth+1)
                  call smooth_nearGW(q(gfunc_pointer(i,gridnum)),
     *                             q(gr_chr),
     *                             nx,ny,nz,ghostwidth)
               end if
            end do
         end if
      else
         if(ltrace)write(*,99)myid,'grid_step2: Not Smoothing'
      end if

      !
      ! Swap advanced (_np1) to current (_n):
      !    All things outside grid_step() assume (_n) current.
      !
      if (ltrace) then
         write(*,99)myid,'grid_step2: Calling grid_test GG',gridnum
         call grid_test(gridnum)
      end if
      if(ltrace)write(*,99)myid,'grid_step2: Swapping pointers',gridnum
      call grid_swap_n_np1(gridnum)
      call load_pointers(gridnum)
      if (ltrace) then
         write(*,99)myid,'grid_step2: Calling grid_test HH',gridnum
         call grid_test(gridnum)
      end if

       !write(*,99)myid,'step2: afterswap gridnum  =',gridnum
       !call field_dump_info(q(gr_psi4R),nx,ny,nz)
      if(ltrace)write(*,99)myid,'grid_step2: Computing derived qts'
      call grid_comp_deriv(gridnum)    ! acts on XXX_n fields

      if (ltrace2) then
            call field_out3d(q(gfunc_pointer(outfnum+1,gridnum)),
     *             gr_t(gridnum),'stepafterADV',
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nx,ny,nz,myid)
            call field_out3d(q(gfunc_pointer(outfnum,gridnum)),
     *             gr_t(gridnum),'stepafter',
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nx,ny,nz,myid)
      end if

      !
      ! Step time at which grid exists:
      !
      if(ltrace)write(*,99)myid,'grid_step2: advancing time'
      call grid_advance_time(gridnum)
      gr_count(gridnum) = gr_count(gridnum) + 1

      !
      ! Keep track of total "work" done:
      !
       total_work = total_work + gr_work(gridnum)
     *              /refine_factor**gr_level(gridnum)


      if (ltrace2)call grid_dump_info(gridnum)
      if (ltrace) then
         write(*,98)myid,'grid_step2: total_work = ',total_work
         write(*,98)myid,'grid_step2: gr_time = ',gr_t(gridnum)
         write(*,99)myid,'grid_step2: gr_count= ',gr_count(gridnum)
         write(*,99)myid,'grid_step2: Done.'
         write(*,99)myid,'grid_step2: Calling grid_test on',gridnum
         call grid_test(gridnum)
         write(*,99)myid,'grid_step2: Done    = ',gridnum
      end if


      return
 98   format('[',I4,'] ',A,3G15.5)
 99   format('[',I4,'] ',A,3I5)
      end    ! END: grid_step2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_2grid_init:                                                          cc
cc                         Init a newly created grid using data from a        cc
cc                        local old sibling grid.                             cc
cc                        NB: extended to deal with siblings w/ different     cc
cc                            parents.                                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_2grid_init(grid, sibling)
      implicit    none
      integer     grid, sibling
      include     'grid.inc'
      include     'glob.inc'
      include     'grid_methods.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     parent
      integer     nxg, nyg, nzg, nxs, nys, nzs, nxt,nyt,nzt,
     *            i, length, owner,
     *            sowner,
     *            oi, oj, ok,
     *            min_i, max_i,
     *            min_j, max_j,
     *            min_k, max_k,
     *            levp, levc, countp, countc
      character(2) tmps,tmps2
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(sibling)

      if (grid_return_level(grid) .ne. grid_return_level(sibling)) then
         write(*,*) 'grid_2grid_init: Grids not on same level'
         call my_exit('grid_2grid_init: Grids not on same level')
      end if

      nxg    = gr_nx(grid)
      nyg    = gr_ny(grid)
      nzg    = gr_nz(grid)

      nxs    = gr_nx(sibling)
      nys    = gr_ny(sibling)
      nzs    = gr_nz(sibling)

      call grid_find_intersection(sibling,grid, length,
     *            min_i, max_i, min_j, max_j, min_k, max_k)

      nxt = max_i - min_i + 1
      nyt = max_j - min_j + 1
      nzt = max_k - min_k + 1

      if (length.le.0) then
         if (ltrace) write(*,*) 'grid_2grid_init: No intersection'
         return
      end if

      !
      ! Determine offsets between grids:
      !
      if (gr_minx(grid).ge.gr_minx(sibling)) then
         oi = 0
      else
         oi = NINT( (gr_minx(sibling)-gr_minx(grid))/gr_h(grid) )
      end if
      if (gr_miny(grid).ge.gr_miny(sibling)) then
         oj = 0
      else
         oj = NINT( (gr_miny(sibling)-gr_miny(grid))/gr_h(grid) )
      end if
      if (gr_minz(grid).ge.gr_minz(sibling)) then
         ok = 0
      else
         ok = NINT( (gr_minz(sibling)-gr_minz(grid))/gr_h(grid) )
      end if

      if (ltrace) then
         call int2str(grid,tmps)
         call int2str(sibling,tmps2)
         call level_tree_dump()
         call grid_dump_info(grid)
         call grid_dump_info(sibling)
         write(*,*) 'grid_2grid_init:  grid   = ', grid
         write(*,*) 'grid_2grid_init: sibling = ', sibling
         write(*,*) 'grid_2grid_init:  owner  = ', owner
         write(*,*) 'grid_2grid_init: length  = ', length
         write(*,*) 'grid_2grid_init: oi      = ', oi
         write(*,*) 'grid_2grid_init: oj      = ', oj
         write(*,*) 'grid_2grid_init: ok      = ', ok
         write(*,*) 'grid_2grid_init: nxg     = ', nxg
         write(*,*) 'grid_2grid_init: nyg     = ', nyg
         write(*,*) 'grid_2grid_init: nzg     = ', nzg
         write(*,*) 'grid_2grid_init: nxt     = ', nxt
         write(*,*) 'grid_2grid_init: nyt     = ', nyt
         write(*,*) 'grid_2grid_init: nzt     = ', nzt
         write(*,*) 'grid_2grid_init: nxs     = ', nxs
         write(*,*) 'grid_2grid_init: nys     = ', nys
         write(*,*) 'grid_2grid_init: nzs     = ', nzs
         write(*,*) 'grid_2grid_init: min_i   = ', min_i
         write(*,*) 'grid_2grid_init: max_i   = ', max_i
         write(*,*) 'grid_2grid_init: min_j   = ', min_j
         write(*,*) 'grid_2grid_init: max_j   = ', max_j
         write(*,*) 'grid_2grid_init: min_k   = ', min_k
         write(*,*) 'grid_2grid_init: max_k   = ', max_k
           call field_out3d(q(gfunc_pointer(1,grid)),gr_t(grid),
     *             'g2gbeforeg'//tmps//tmps2,
     *             gr_minx(grid),gr_maxx(grid),
     *             gr_miny(grid),gr_maxy(grid),
     *             gr_minz(grid),gr_maxz(grid),
     *             gr_nx(grid),gr_ny(grid),gr_nz(grid),myid)
           call field_out3d(q(gfunc_pointer(1,sibling)),gr_t(sibling),
     *             'g2gbefores'//tmps//tmps2,
     *             gr_minx(sibling),gr_maxx(sibling),
     *             gr_miny(sibling),gr_maxy(sibling),
     *             gr_minz(sibling),gr_maxz(sibling),
     *             gr_nx(sibling),gr_ny(sibling),gr_nz(sibling),myid)
      end if

      !
      ! Copies points of the sibling that overlap with grid to
      ! temporary storage
      !
      !
      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL       .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRALFV     .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL .or.
     *        gfunc_name(i) .eq. 'mask'               .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DER       ) then
            !
            ! Copy grid funcation data from sibling
            ! to gr_tmp storage:
            !
            call copy_overlap(
     *                       q(gr_tmp),
     *                       q(gfunc_pointer(i,sibling)),
     *                       min_i,min_j,min_k,
     *                       nxt,nyt,nzt,
     *                       nxs,nys,nzs )
            !
            ! Copy from storage into the grid to be initted:
            !
            call copy_interp_section( 
     *                       q(gfunc_pointer(i,grid)),
     *                       q(gr_tmp),
     *                       oi, oj, ok,
     *                       nxt,nyt,nzt,
     *                       nxg,nyg,nzg, refine_factor)
         end if
      end do
 
      if (ltrace) then
           call field_out3d(q(gfunc_pointer(1,grid)),gr_t(grid),
     *             'g2gafterg'//tmps//tmps2,
     *             gr_minx(grid),gr_maxx(grid),
     *             gr_miny(grid),gr_maxy(grid),
     *             gr_minz(grid),gr_maxz(grid),
     *             gr_nx(grid),gr_ny(grid),gr_nz(grid),myid)
         write(*,*) 'grid_2grid_init: Done.'
      end if

      return
      end    ! END: grid_2grid_init

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_initialize:                                                          cc
cc                                                                            cc
cc      High-level initialization of grids (meaning it calls                  cc
cc          grid_init_parent, grid_2grid_init etc. )                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_initialize(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      include 'glob.inc'
      include 'mpif.h'
      include 'mpi_stuff.inc'
      include 'action.inc'
      include 'grid_methods.inc'
      integer  level, gi_tmp
      integer  nx, ny, nz
      real(kind=8) t, minx, maxx, miny,maxy,minz,maxz
      character(2) tmps,tmps2
      ! Are grid, parent, and sibling local to this processor?
      logical     local, plocal, slocal
      logical     double_equal
      external    double_equal
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )
      ! To track bugs down via SDF output:
      logical     ltraceS
      parameter ( ltraceS = .false. )
      ! Which field to output if above is true? 
      !    (look in had/include/fields.inc for number)
      integer     fieldnum
      parameter ( fieldnum  = 17 )

      level = grid_return_level(gi)

      local = grid_is_local(gi)
      nx    = gr_nx(gi)
      ny    = gr_ny(gi)
      nz    = gr_nz(gi)
      t     = gr_t(gi)
      minx  = gr_minx(gi)
      maxx  = gr_maxx(gi)
      miny  = gr_miny(gi)
      maxy  = gr_maxy(gi)
      minz  = gr_minz(gi)
      maxz  = gr_maxz(gi)

      if (ltrace2 .or. (ltrace.and.local) ) then
         write(*,99) myid, 'grid_initialize: Initting grid: ',gi
         if (local) then
            write(*,99) myid, 'grid_initialize: Level = ',level
            write(*,98) myid, 'grid_initialize: time  = ',gr_t(gi)
         end if
      end if
      if (local.and.ltrace) call level_tree_dump()


      !
      ! Init grid either by:
      !            -- if (t=0), 'exact' initial data
      !            -- interpolation from coarse grid parent
      !               & transfer values from an old fine grid
      !if (gr_t(gi).eq.0.d0 .and. level.eq.0) then
      !if (gr_t(gi).eq.0.d0) then
      if (double_equal(gr_t(gi),0.d0)) then
         if(ltrace.and.local)write(*,99) myid, ' Init with init data',gi
         if (local) call grid_init_t0(gi)
         if (local) call grid_copy_n_np1(gi)
         if(ltrace.and.local)write(*,99)myid,'DONEInit w/ init data',gi
      else                                
         !
         ! Init from all parents:
         !
         if (ltrace.and.local) write(*,99) myid, ' Init from parents'
         ! No parents if on level 0:
         if (level .eq. 0 ) goto 68
         if (local .and. ltraceS) then
            write(*,99) myid, 'Zeroing (-1.3) out first field'
            call load_scal1d(q(gfunc_pointer(1,gi)),
     *                    -1.3d0,gr_nx(gi)*gr_ny(gi)*gr_nz(gi))
         end if
         !
         ! When interpolating from parents on a newly created
         ! level (so no siblings), we need to know where the
         ! mask is even before we interpolate. We also need
         ! to call this afterward to mask out the fields.
         !
         if (local .and.
     *       level .eq. level_return_finest(Levelp,maxlev)) then
            !
            if (ltrace) write(*,99)myid,'   Init mask w/ bhmask on ',gi
            call grid_bhmask_gfuncs(gi)
         end if
         gi_tmp = level_return_start(level-1)
 47      if ( grid_return_existence(gi_tmp) ) then
            if (.not. grid_intersect(gi,gi_tmp)) goto 48
            plocal = grid_is_local(gi_tmp)
            if (local .and. plocal) then
               if (ltrace)write(*,99)myid,'locally from parent:', gi_tmp
               call grid_init(gi, gi_tmp, refine_factor)
            else if (local) then
               if (ltrace)write(*,99) myid, '    grid_init_child',gi_tmp
               call grid_init_child(gi,gi_tmp)
            else if (plocal) then
               if (ltrace)write(*,99) myid,'    grid_init_parent',gi_tmp
               call grid_init_parent(gi,gi_tmp)
            else 
               if(ltrace.and.local)write(*,99)myid,'nether local',gi_tmp
            end if
            !if (ltrace .and. local) call grid_output_sdf(gi)
            if (local .and. ltraceS) then
              call int2str(gi_tmp,tmps)
              call int2str(gi,tmps2)
              write(*,99)myid,'   Output First'//tmps2//tmps
              write(*,98)myid,'   time = ',t
              write(*,99)myid,'   nx,ny,nz  = ',nx,ny,nz
              call field_out3d(q(gfunc_pointer(fieldnum,gi)),t,
     *             'First'//tmps2//tmps,
     *             minx,maxx,miny, maxy,minz,maxz, nx,ny,nz,myid)
            end if
 48         gi_tmp = grid_return_sibling(gi_tmp)
            if (grid_return_existence(gi_tmp)) goto 47
         end if
         if (local .and. ltraceS) then
           write(*,99)myid,'   Outputting First field_afterpar....'
           call field_out3d(q(gfunc_pointer(fieldnum,gi)),gr_t(gi),
     *             'D_afterpar',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
         end if
         !
         ! Need to deal with the mask:
         ! If we are on the finest level (or creating it), 
         !    then we just call grid_bhmask()
         ! Otherwise, siblings should be able to initialize it:
         !
         if (local .and.
     *       level .eq. level_return_finest(Levelp,maxlev)) then
            !
            if (ltrace) write(*,99)myid,'   Init mask w/ bhmask on ',gi
            call grid_bhmask_gfuncs(gi)
         end if
         !if (local.and.ltrace) call grid_dump_info(gi)
         if (local .and. ltrace) then
           write(*,99)myid,'   Outputting First field_afterparM....'
           call field_out3d(q(gfunc_pointer(fieldnum,gi)),gr_t(gi),
     *             'D_afterparM',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
           call field_out3d(q(gr_mask),gr_t(gi),
     *             'ginit_maskA',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
         end if
         !
         ! Init from siblings:
         !
 68      if (local.and.ltrace) call level_tree_dump()
      !if (.false.) then
         if (local.and.ltrace) write(*,99) myid,'  Init from siblings'
         gi_tmp = level_return_start(level)
 69      if ( grid_return_existence(gi_tmp) ) then
            if (      gi_tmp .eq. gi
     *           .or. .not.grid_intersect(gi,gi_tmp)) goto 70
            if(local.and.ltrace)write(*,99)myid,'Init fm old sib',gi_tmp
            slocal = grid_is_local(gi_tmp)
            if (local .and. slocal) then
               if (ltrace)write(*,99) myid, ' both    local',gi_tmp
               call grid_2grid_init(gi, gi_tmp)
            else if (local) then
               if (ltrace)write(*,99) myid, ' gi      local',gi_tmp
               call grid_2grid_init_slave(gi, gi_tmp)
            else if (slocal) then
               if (ltrace)write(*,99) myid, ' gi_tmp  local',gi_tmp
               call grid_2grid_init_master(gi,gi_tmp)
            else
               if (ltrace)write(*,99) myid, ' neither local',gi,gi_tmp
            end if
 70         gi_tmp = grid_return_sibling(gi_tmp)
            goto 69
         end if
      !end if
         if (ltrace) write(*,99) myid, 'No more siblings: ',gi_tmp
         ! copy just initted field to advanced levels
         if (local) call grid_copy_n_np1(gi)
      end if
      if (local .and. ltraceS) then
           write(*,99)myid,'   Outputting First field_aftersib....'
           call field_out3d(q(gfunc_pointer(fieldnum,gi)),gr_t(gi),
     *             'D_aftersib',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
      end if
      if (local .and. ltrace) then
           write(*,97)myid,' Outputting for grid: ',gi,gr_t(gi)
           write(*,99)myid,'   Outputting First field_aftersibM....'
           if (ltraceS) then
           call field_out3d(q(gfunc_pointer(fieldnum,gi)),gr_t(gi),
     *             'D_aftersibM',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
           call field_out3d(q(gr_mask),gr_t(gi),
     *             'ginit_maskB',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
           end if
      end if
      !
      ! Now compute derived quanitities no matter how
      ! grid was initialized
      !
      if (local) then
         if(ltrace)write(*,99)myid,' calling grid_comp_deriv on',gi
         call grid_comp_deriv(gi)
         !write(*,*) 'grid_initialize: calling grid_symm on: ',gi
         if (assume_symmetry.ne.0) call grid_symm(gi)
         !if(ltrace) call grid_test(gi)
      end if
      if (local.and.ltrace) call grid_output_sdf(gi)
      if (local .and. ltraceS) then
        call grid_output_sdf(gi)
        call field_out3d(q(gr_mask),gr_t(gi), 'mask_init',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
      end if
      if (local.and.ltrace) call grid_dump_info(gi)

      if (ltrace) then
         write(*,99)myid, 'grid_initialize: Done initting grid: ***',gi
         call level_test_local()
      end if

 97   format('[',I3,'] ',A,I5,3F10.5)
 98   format('[',I3,'] ',A,3F10.5)
 99   format('[',I3,'] ',A,3I5)

      return
      end      ! END: grid_initialize

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_bbox_intersection:                                                   cc
cc              Return some measure of overlap between a bounding box         cc
cc              and grid, used for determining best processor to own          cc
cc              a new grid about to be created.                               cc
cc              Dervied from grid_find_intersection.                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function  grid_bbox_intersection(gi, minx,maxx, miny,maxy,
     *                  minz, maxz)
      implicit    none
      integer     gi
      real(kind=8) minx,maxx, miny,maxy, minz,maxz
      include    'grid.inc'
      include    'glob.inc'
      integer     length, si,fi, sj,fj, sk, fk

      logical     ltrace
      parameter ( ltrace = .false.)

      length = 0
      si     = 0
      fi     = 0
      sj     = 0
      fj     = 0
      sk     = 0
      fk     = 0

      ! 
      ! Find bounding box of intersection (if any):
      !
      if ( gr_minx(gi) .ge. minx ) then
         if ( gr_minx(gi) .gt. maxx ) goto 100! no intersection
         si = 1
      else
         if ( gr_maxx(gi) .lt. minx ) goto 100! no intersection
         si = 1 + NINT( (minx - gr_minx(gi)) / gr_h(gi) )
      end if
      if ( gr_maxx(gi) .le. maxx ) then
         fi = gr_nx(gi)
      else
         fi = 1 + NINT( (maxx - gr_minx(gi)) / gr_h(gi) )
      end if
      !
      if ( gr_miny(gi) .ge. miny ) then
         if ( gr_miny(gi) .gt. maxy ) goto 100! no intersection
         sj = 1
      else
         if ( gr_maxy(gi) .lt. miny ) goto 100! no intersection
         sj = 1 + NINT( (miny - gr_miny(gi)) / gr_h(gi) )
      end if
      if ( gr_maxy(gi) .le. maxy ) then
         fj = gr_ny(gi)
      else
         fj = 1 + NINT( (maxy - gr_miny(gi)) / gr_h(gi) )
      end if
      !
      if ( gr_minz(gi) .ge. minz ) then
         if ( gr_minz(gi) .gt. maxz ) goto 100! no intersection
         sk = 1
      else
         if ( gr_maxz(gi) .lt. minz ) goto 100! no intersection
         sk = 1 + NINT( (minz - gr_minz(gi)) / gr_h(gi) )
      end if
      if ( gr_maxz(gi) .le. maxz ) then
         fk = gr_nz(gi)
      else
         fk = 1 + NINT( (maxz - gr_minz(gi)) / gr_h(gi) )
      end if

      length = (fk-sk+1)*(fj-sj+1)*(fi-si+1)

 100  continue

      grid_bbox_intersection = length

      return
      end       ! END: grid_bbox_intersection

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_mask_error:                                                          cc
cc                  Computing the error is project-dependent, but after       cc
cc                  it is computed, two things need to be done to dependent   cc
cc                  on the mask (gr_chr):                                     cc
cc                    (0)  Make sure computed error is strictly >= 0          cc
cc                    (1)  for all AMR boundaries, the error within the       cc
cc                         ghost region needs to be zeroed out there so that  cc
cc                         the clusterer does not try to put a subgrid there. cc
cc                    (2)  for all DD  boundaries, the error needs to be      cc
cc                            zeroed out so that the error from overlapping   cc
cc                            grids is used (else we run into the problem     cc
cc                            where the error gets erroneously high).         cc
cc                    (3)  for all AMR boundaries, the error needs to be      cc
cc                         set to a negative value in region gw/r.            cc
cc                                                                            cc
cc                   10/29/16 SLL: Fixed so that no point is set to 0         cc
cc                   if it has already been set negative.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_mask_error(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      include 'grid_methods.inc'
      include 'glob.inc'
      include 'chr.inc'
      integer  index, index_tmp, myid
      integer  i,   j,   k
      integer  ii,  jj,  kk
      integer  nx,  ny,  nz

      integer     proc_return_myid
      external    proc_return_myid
      real(kind=8)myl2norm3d, errornorm
      external    myl2norm3d

      !
      !  The width of points at Domain Deco boundaries
      !  to mask off.                                  
      !
      integer     ddw
      parameter ( ddw = 1 )

      logical     ltrace
      parameter ( ltrace = .false. )

      myid = proc_return_myid()
      call load_pointers(gi)

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

      errornorm = myl2norm3d(q(gr_error),nx,ny,nz)
      if (errornorm.gt.1.0d7) then
         write(*,99)myid,'grid_mask_error: Possible problem on gi:',gi
         write(*,98)myid,'grid_mask_error: error norm=',errornorm
         call grid_dump_info(gi)
      end if

      if (ltrace) then
         !call load_scal3d( q(gr_error), -10.d0, nx, ny, nz)
 98      format('[',I3,'] ',A,3G15.7)
 99      format('[',I3,'] ',A,3I5)
         write(*,99)myid,'grid_mask_error: Masking error on:',gi
         write(*,99)myid,'grid_mask_error: nx/y/z = ',nx,ny,nz
         write(*,99)myid,'grid_mask_error: Outputting mask_errorA'
         call field_out3d(q(gr_error),gr_t(gi), 'mask_errorA',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
         call field_out3d(q(gr_chr),gr_t(gi), 'mask_chrA',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if

      !
      ! Ensure error starts nowhere negative
      !
      call vector_abs1D(q(gr_error),nx*ny*nz)

      !
      ! X-boundaries
      !
      i = 1
      do k = 1, nz
      do j = 1, ny
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_amr_bdy) ) then
            if(ltrace.and.j.eq.ny/2.and.k.eq.nz/2) then
               write(*,99)myid,' low x AMR'
               write(*,99)myid,'grid_mask_error: i/j/k: ',i,j,k
               write(*,99)myid,'grid_mask_error: chr: ',
     *                                 NINT(q(gr_chr+index))
            end if
            do ii = 1, ghostwidth
               index_tmp = (k-1)*ny*nx + (j-1)*nx + (ii-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
            do ii = 1, ghostwidth/refine_factor
               index_tmp = (k-1)*ny*nx + (j-1)*nx + (ii-1)
               q(gr_error+index_tmp) = -1.d0
            end do
         else if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_deco_bdy) ) then
            if(ltrace.and.j.eq.ny/2.and.k.eq.nz/2) then
               write(*,99)myid,' low x DD'
               write(*,99)myid,'grid_mask_error: i/j/k: ',i,j,k
               write(*,99)myid,'grid_mask_error: chr: ',
     *                                 NINT(q(gr_chr+index))
            end if
            do ii = 1, ddw
               index_tmp = (k-1)*ny*nx + (j-1)*nx + (ii-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
         end if
      end do
      end do
      if (ltrace) then
         call field_out3d(q(gr_error),gr_t(gi), 'mask_erroPXX',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if
      !
      i = nx
      do k = 1, nz
      do j = 1, ny
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_amr_bdy)) then
            if(ltrace.and.j.eq.ny/2.and.k.eq.nz/2) then
               write(*,99)myid,' hi  x AMR'
               write(*,99)myid,'grid_mask_error: i/j/k: ',i,j,k
               write(*,99)myid,'grid_mask_error: chr: ',
     *                                 NINT(q(gr_chr+index))
            end if
            do ii = nx, nx-ghostwidth+1, -1
               index_tmp = (k-1)*ny*nx + (j-1)*nx + (ii-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
            do ii = nx, nx-ghostwidth/refine_factor+1, -1
               index_tmp = (k-1)*ny*nx + (j-1)*nx + (ii-1)
               q(gr_error+index_tmp) = -1.d0
            end do
         else if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_deco_bdy) ) then
            if(ltrace.and.j.eq.ny/2.and.k.eq.nz/2) then
               write(*,99)myid,' hi  x DD'
               write(*,99)myid,'grid_mask_error: i/j/k: ',i,j,k
               write(*,99)myid,'grid_mask_error: chr: ',
     *                                 NINT(q(gr_chr+index))
            end if
            do ii = nx, nx-ddw+1, -1
               index_tmp = (k-1)*ny*nx + (j-1)*nx + (ii-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
         end if
      end do
      end do
      if (ltrace) then
         call field_out3d(q(gr_error),gr_t(gi), 'mask_erroPX',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if
      !
      ! Y-boundaries
      !
      j = 1
      do k = 1, nz
      do i = 1, nx
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_amr_bdy) ) then
            do jj = 1, ghostwidth
               index_tmp = (k-1)*ny*nx + (jj-1)*nx + (i-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
            do jj = 1, ghostwidth/refine_factor
               index_tmp = (k-1)*ny*nx + (jj-1)*nx + (i-1)
               q(gr_error+index_tmp) = -1.d0
            end do
         else if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_deco_bdy) ) then
            do jj = 1, ddw
               index_tmp = (k-1)*ny*nx + (jj-1)*nx + (i-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
         end if
      end do
      end do
      if (ltrace) then
         call field_out3d(q(gr_error),gr_t(gi), 'mask_erroPYY',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if
      !
      j = ny
      do k = 1, nz
      do i = 1, nx
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_amr_bdy) ) then
            do jj = ny, ny-ghostwidth+1, -1
               index_tmp = (k-1)*ny*nx + (jj-1)*nx + (i-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
            do jj = ny, ny-ghostwidth/refine_factor+1, -1
               index_tmp = (k-1)*ny*nx + (jj-1)*nx + (i-1)
               q(gr_error+index_tmp) = -1.d0
            end do
         else if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_deco_bdy) ) then
            do jj = ny, ny-ddw+1, -1
               index_tmp = (k-1)*ny*nx + (jj-1)*nx + (i-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
         end if
      end do
      end do
      if (ltrace) then
         call field_out3d(q(gr_error),gr_t(gi), 'mask_erroPY',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if
      !
      ! Z-boundaries
      !
      k = 1
      do j = 1, ny
      do i = 1, nx
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_amr_bdy) ) then
            do kk = 1, ghostwidth
               index_tmp = (kk-1)*ny*nx + (j-1)*nx + (i-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
            do kk = 1, ghostwidth/refine_factor
               index_tmp = (kk-1)*ny*nx + (j-1)*nx + (i-1)
               q(gr_error+index_tmp) = -1.d0
            end do
         else if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_deco_bdy) ) then
            do kk = 1, ddw
               index_tmp = (kk-1)*ny*nx + (j-1)*nx + (i-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
         end if
      end do
      end do
      if (ltrace) then
         call field_out3d(q(gr_error),gr_t(gi), 'mask_erroPZZ',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if
      !
      k = nz
      do j = 1, ny
      do i = 1, nx
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_amr_bdy) ) then
            do kk = nz, nz-ghostwidth+1, -1
               index_tmp = (kk-1)*ny*nx + (j-1)*nx + (i-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
            do kk = nz, nz-ghostwidth/refine_factor+1, -1
               index_tmp = (kk-1)*ny*nx + (j-1)*nx + (i-1)
               q(gr_error+index_tmp) = -1.d0
            end do
         else if ( NINT(q(gr_chr+index)) .eq. NINT(CHR_deco_bdy) ) then
            do kk = nz, nz-ddw+1, -1
               index_tmp = (kk-1)*ny*nx + (j-1)*nx + (i-1)
               if( q(gr_error+index_tmp).ge.0)
     *         q(gr_error+index_tmp) = 0.d0
            end do
         end if
      end do
      end do

      if (ltrace) then
         write(*,99)myid,'grid_mask_error: Outputting mask_errorB'
         call field_out3d(q(gr_error),gr_t(gi), 'mask_errorB',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
         write(*,99)myid,'grid_mask_error: Done.'
      end if

      return
      end       ! END: grid_mask_error


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_bhmask:                                                              cc
cc                High level masking routine which handles:                   cc
cc                (1) Setting the mask and masking the fields by calling      cc
cc                      grid_bhmask_mask()                                    cc
cc                (2) Determining which points are to be repopulated          cc
cc                (3) Calling grid_repopulate()                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_bhmask(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      include 'mask.inc'
      real(kind=8) r1_2, r2_2
      integer  nx,ny,nz
      integer  i,j,k
      integer  l, gf
      integer  index, indexip1, indexim1
      integer         indexjp1, indexjm1
      integer         indexkp1, indexkm1
      integer  numasked, numrepop, myid
      logical  inmaskedregion
      logical  anyneighborunmasked
      integer  proc_return_myid
      external proc_return_myid
      
      logical     ltrace
      parameter ( ltrace = .false.)

      if (ltrace) then
         myid = proc_return_myid()
         write(*,99)myid,'grid_bhmask: Check mask initially, label 0',gi
         call grid_check_mask(gi,0)
      end if
      !
      ! Nothing to do if no masked regions defined
      !
      if (num_masks .le. 0) return

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

      !
      ! Reset mask to normalized (-1,+1)...
      ! the hyper-side masking changes the values based on
      ! geometry but it is supposed to change it back
      ! by calling this routine at the end of the rk process.
      ! I have not located how it's not fixed back, so
      ! we will call again here.
      !
      if (ltrace) call field_dump_info(q(gr_mask),nx,ny,nz)
      call load_pointers(gi)
      call unfixup_mask(q(gr_mask),nx,ny,nz)

      if (ltrace) then
        write(*,98)myid,'grid_bhmask: time  = ',gr_t(gi)
        write(*,98)myid,'grid_bhmask: h     = ',gr_h(gi)
        do k = 1, max_num_masks
           do i = 1, 6
           write(*,98)myid,'grid_bhmask: mask_coords:',mask_coords(i,k)
           end do
        end do
        call field_out3d(q(gr_mask),gr_t(gi), 'bhmask_pre',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if

      numasked = 0
      numrepop = 0

      call grid_bhmask_mask(gi)

      call load_pointers(gi)

      !
      ! Second, determine which points can be unmasked.
      !    (it is important to separate these two loops
      !     because otherwise the masked region is not 
      !     guaranteed to be rectangular because near  
      !     the corners some points will be able to    
      !     be repopulated that otherwise shouldn't    
      !     simply because the loop to mask them off   
      !     has not gotten to them yet):
      !
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         !
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         !
         inmaskedregion = .false.
         do l = 1, max_num_masks
            if (bh_true(l)) then
               !
               inmaskedregion = inmaskedregion .or. (
     *           ( q(gr_x(gi)+i-1)  .gt. (mask_coords(1,l)+FUZZ) ) .and.
     *           ( q(gr_x(gi)+i-1)  .lt. (mask_coords(2,l)-FUZZ) ) .and.
     *           ( q(gr_y(gi)+j-1)  .gt. (mask_coords(3,l)+FUZZ) ) .and.
     *           ( q(gr_y(gi)+j-1)  .lt. (mask_coords(4,l)-FUZZ) ) .and.
     *           ( q(gr_z(gi)+k-1)  .gt. (mask_coords(5,l)+FUZZ) ) .and.
     *           ( q(gr_z(gi)+k-1)  .lt. (mask_coords(6,l)-FUZZ) ) )
               !
            end if
         end do
         !
         if (NINT(q(gr_mask+index)) .lt. 0 .and.
     *       .not. inmaskedregion ) then
            ! This point should *no* longer be masked:
            indexip1 = (k-1)*ny*nx + (j-1)*nx + (i  )
            indexim1 = (k-1)*ny*nx + (j-1)*nx + (i-2)
            !
            indexjp1 = (k-1)*ny*nx + (j  )*nx + (i-1)
            indexjm1 = (k-1)*ny*nx + (j-2)*nx + (i-1)
            !
            indexkm1 = (k-2)*ny*nx + (j-1)*nx + (i-1)
            indexkp1 = (k  )*ny*nx + (j-1)*nx + (i-1)
            !
            anyneighborunmasked = .false.
            if (i.gt.1 ) anyneighborunmasked = anyneighborunmasked
     *         .or. ( NINT(q(gr_mask+indexim1)).eq.BHMASK_UNMASKED )
            if (i.lt.nx) anyneighborunmasked = anyneighborunmasked
     *         .or. ( NINT(q(gr_mask+indexip1)).eq.BHMASK_UNMASKED )
            !
            if (j.gt.1 ) anyneighborunmasked = anyneighborunmasked
     *         .or. ( NINT(q(gr_mask+indexjm1)).eq.BHMASK_UNMASKED )
            if (j.lt.ny) anyneighborunmasked = anyneighborunmasked
     *         .or. ( NINT(q(gr_mask+indexjp1)).eq.BHMASK_UNMASKED )
            !
            if (k.gt.1 ) anyneighborunmasked = anyneighborunmasked
     *         .or. ( NINT(q(gr_mask+indexkm1)).eq.BHMASK_UNMASKED )
            if (k.lt.nz) anyneighborunmasked = anyneighborunmasked
     *         .or. ( NINT(q(gr_mask+indexkp1)).eq.BHMASK_UNMASKED )
            !
            if (anyneighborunmasked) then
               ! Only unmask if next to unmasked region
               !   (this limits the motion of the mask to a single
               !    grid point in any direction):
               q(gr_mask+index) = 1.d0 * BHMASK_REPOPULATE
               numrepop = numrepop + 1
            end if
         end if
         !
      end do
      end do
      end do

      if (ltrace) then
        write(*,99)myid,'grid_bhmask: Outputting bhmask:'
        write(*,99)myid,'grid_bhmask: numasked  = ',numasked
        write(*,99)myid,'grid_bhmask: numrepop  = ',numrepop
        write(*,99)myid,'grid_bhmask: .....'
        call field_out3d(q(gr_mask),gr_t(gi), 'bhmask',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if

      call grid_repopulate(gi)

      if (ltrace) then
        call field_out3d(q(gr_mask),gr_t(gi), 'bhmaskp',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
         write(*,99)myid,'grid_bhmask: Check mask at the end,label 1',gi
         call grid_check_mask(gi,1)
        write(*,99)myid,'grid_bhmask: Done.'
      end if


 98   format('[',I3,'] ',A,3F15.9)
 99   format('[',I3,'] ',A,3I7)
      return
      end       ! END: grid_bhmask

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_bhmask_mask:                                                         cc
cc                    Lower level routine to set mask and mask fields         cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_bhmask_mask(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      include 'mask.inc'
      integer  nx,ny,nz
      integer  i,j,k
      integer  l, gf
      integer  index
      logical  inmaskedregion
      integer  numasked, myid
      integer  proc_return_myid
      external proc_return_myid
      
      logical     ltrace
      parameter ( ltrace = .false.)

      myid = proc_return_myid()

      !
      ! Nothing to do if no masked regions defined
      !
      if (num_masks .le. 0) then
         if (ltrace) write(*,98)myid,'grid_bhmask_mask: No masks,return'
         return
      end if

      call load_pointers(gi)

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

      if (ltrace) then
        write(*,99)myid,'grid_bhmask_mask: gi    = ',gi
        write(*,98)myid,'grid_bhmask_mask: time  = ',gr_t(gi)
        write(*,98)myid,'grid_bhmask_mask: h     = ',gr_h(gi)
        do k = 1, max_num_masks
           do i = 1, 6
           write(*,98)myid,'grid_bhmask_mask:msk_crds:',mask_coords(i,k)
           end do
        end do
        call field_out3d(q(gr_mask),gr_t(gi), 'bhmask_mask',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if

      numasked = 0
      !
      !
      ! First, mask out points which may not have been masked
      !
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         !
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         !
         inmaskedregion = .false.
         do l = 1, max_num_masks
            if (bh_true(l)) then
               !
               inmaskedregion = inmaskedregion .or. (
     *           ( q(gr_x(gi)+i-1)  .gt. (mask_coords(1,l)+FUZZ) ) .and.
     *           ( q(gr_x(gi)+i-1)  .lt. (mask_coords(2,l)-FUZZ) ) .and.
     *           ( q(gr_y(gi)+j-1)  .gt. (mask_coords(3,l)+FUZZ) ) .and.
     *           ( q(gr_y(gi)+j-1)  .lt. (mask_coords(4,l)-FUZZ) ) .and.
     *           ( q(gr_z(gi)+k-1)  .gt. (mask_coords(5,l)+FUZZ) ) .and.
     *           ( q(gr_z(gi)+k-1)  .lt. (mask_coords(6,l)-FUZZ) ) )
               !
            end if
         end do
         !
         if (inmaskedregion) then
            ! This point should       be masked:
            q(gr_mask+index) = 1.d0 * BHMASK_MASKED
            numasked = numasked + 1
            do gf = 1, num_gfuncs
               if ( gfunc_name(gf) .ne. 'mask'   .and.
     *              gfunc_name(gf) .ne. 'cctk_x' .and.
     *              gfunc_name(gf) .ne. 'cctk_y' .and.
     *              gfunc_name(gf) .ne. 'cctk_z' .and.
     *              gfunc_name(gf) .ne. 'xphys'  .and.
     *              gfunc_name(gf) .ne. 'yphys'  .and.
     *              gfunc_name(gf) .ne. 'zphys'  .and.
     *              gfunc_name(gf) .ne. 'r'      .and.
     *              gfunc_name(gf) .ne. 'chr'        ) then
                  q(gfunc_pointer(gf,gi)+index) = 0.d0
               end if
            end do
         end if
         !
      end do
      end do
      end do


      if (ltrace) then
        write(*,99)myid,'grid_bhmask_mask: Outputting bhmask_maskP'
        call field_out3d(q(gr_mask),gr_t(gi), 'bhmask_maskP',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
        write(*,99)myid,'grid_bhmask_mask: Mask info:'
        call field_dump_info(q(gr_mask),nx,ny,nz)
        write(*,99)myid,'grid_bhmask_mask: Done.'
      end if

 98   format('[',I3,'] ',A,3F15.9)
 99   format('[',I3,'] ',A,3I7)
      return
      end       ! END: grid_bhmask_mask


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_bhmask_gfuncs:                                                       cc
cc                    Lowest level, just masking fields for which mask        cc
cc                    indicates the mask exist....in other words, it does     cc
cc                    not actually set the mask itself, just looks at it.     cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_bhmask_gfuncs(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      include 'mask.inc'
      integer  nx,ny,nz
      integer  i,j,k
      integer  l, gf
      integer  index
      logical  inmaskedregion
      integer  numasked, myid
      integer  proc_return_myid
      external proc_return_myid
      
      logical     ltrace
      parameter ( ltrace = .false.)

      myid = proc_return_myid()

      !
      ! Nothing to do if no masked regions defined
      !
      if (num_masks .le. 0) then
         if(ltrace)write(*,98)myid,'grid_bhmask_gfuncs: No masks,return'
         return
      end if

      call load_pointers(gi)

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

      if (ltrace) then
        write(*,99)myid,'grid_bhmask_gfuncs: gi    = ',gi
        write(*,98)myid,'grid_bhmask_gfuncs: time  = ',gr_t(gi)
        write(*,98)myid,'grid_bhmask_gfuncs: h     = ',gr_h(gi)
        call field_out3d(q(gr_mask),gr_t(gi), 'bhmask_gfuncs',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if

      !
      !
      ! First, mask out points which may not have been masked
      !
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         !
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         !
         inmaskedregion = q(gr_mask+index) .lt. 0.d0
         !
         if (inmaskedregion) then
            ! This point should       be masked:
            do gf = 1, num_gfuncs
               if ( gfunc_name(gf) .ne. 'mask'   .and.
     *              gfunc_name(gf) .ne. 'cctk_x' .and.
     *              gfunc_name(gf) .ne. 'cctk_y' .and.
     *              gfunc_name(gf) .ne. 'cctk_z' .and.
     *              gfunc_name(gf) .ne. 'xphys'  .and.
     *              gfunc_name(gf) .ne. 'yphys'  .and.
     *              gfunc_name(gf) .ne. 'zphys'  .and.
     *              gfunc_name(gf) .ne. 'r'      .and.
     *              gfunc_name(gf) .ne. 'chr'        ) then
                  q(gfunc_pointer(gf,gi)+index) = 0.d0
               end if
            end do
         end if
         !
      end do
      end do
      end do


      if (ltrace) then
        write(*,99)myid,'grid_bhmask_gfuncs: Outputting bhmask_gfuncsB'
        call field_out3d(q(gr_mask),gr_t(gi), 'bhmask_gfuncsB',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
        write(*,99)myid,'grid_bhmask_gfuncs: Mask info:'
        call field_dump_info(q(gr_mask),nx,ny,nz)
        write(*,99)myid,'grid_bhmask_gfuncs: Done.'
      end if

 98   format('[',I3,'] ',A,3F15.9)
 99   format('[',I3,'] ',A,3I7)
      return
      end       ! END: grid_bhmask_gfuncs

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_repopulate:                                                          cc
cc                                                                            cc
cc      Where points come out from being masked, do the following:            cc
cc          (1) reset the mask so that it shows points not masked anymore     cc
cc          (2) repopulate those points on all fields via extrapolation       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_repopulate(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      include 'mask.inc'
      integer  nx,ny,nz, i,j,k, index, mask, gf
      !
      ! Extrapolation variables:
      !
      real(kind=8)  coef_ip1, coef_ip2,coef_ip3, coef_ip4,
     *              coef_im1, coef_im2,coef_im3, coef_im4,
     *              coef_jp1, coef_jp2,coef_jp3, coef_jp4,
     *              coef_jm1, coef_jm2,coef_jm3, coef_jm4,
     *              coef_kp1, coef_kp2,coef_kp3, coef_kp4,
     *              coef_km1, coef_km2,coef_km3, coef_km4,
     *              sum_coef
      integer  index_ip1, index_ip2,index_ip3, index_ip4,
     *         index_im1, index_im2,index_im3, index_im4,
     *         index_jp1, index_jp2,index_jp3, index_jp4,
     *         index_jm1, index_jm2,index_jm3, index_jm4,
     *         index_kp1, index_kp2,index_kp3, index_kp4,
     *         index_km1, index_km2,index_km3, index_km4
      integer  num_unrepop, fptr
      integer  numasked, numrepop, myid
      integer  proc_return_myid
      external proc_return_myid

      !
      ! Use linear extrapolation or cubic:
      !
      logical     linear
      parameter ( linear  = .false.)
      integer     int_order
      parameter ( int_order = 1 )

      logical     ltrace
      parameter ( ltrace  = .false.)
      logical     ltrace2
      parameter ( ltrace2 = .false.)

      call load_pointers(gi)

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

      ! Keep count of how many points are set to be repopulated,
      ! but that cannot be:
      num_unrepop = 0
      numasked    = 0
      numrepop    = 0

      if (ltrace) then
        myid = proc_return_myid()
        write(*,99)myid,'grid_repopulate: gi     = ',gi
        write(*,98)myid,'grid_repopulate: time   = ',gr_t(gi)
        write(*,99)myid,'grid_repopulate: nx/y/z = ',nx,ny,nz
      end if

      !
      ! These two sets of loops need to be separate
      ! so that we don't try to extrapolate using
      ! points which themselves were just extrapolated
      !
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         mask  = NINT( q( gr_mask+index ) )
         if (mask .eq. BHMASK_REPOPULATE) then
            ! This point was masked with the old mask,
            ! but is no longer masked so we need to
            ! repopulate all the fields at those points:
            if (ltrace2)write(*,*)'grid_repopulate:     unmasked:',i,j,k
            index_im4 = (k-1)*ny*nx + (j-1)*nx + (i-5)
            index_im3 = (k-1)*ny*nx + (j-1)*nx + (i-4)
            index_im2 = (k-1)*ny*nx + (j-1)*nx + (i-3)
            index_im1 = (k-1)*ny*nx + (j-1)*nx + (i-2)
            index_ip1 = (k-1)*ny*nx + (j-1)*nx + (i  )
            index_ip2 = (k-1)*ny*nx + (j-1)*nx + (i+1)
            index_ip3 = (k-1)*ny*nx + (j-1)*nx + (i+2)
            index_ip4 = (k-1)*ny*nx + (j-1)*nx + (i+3)
            !
            index_jm4 = (k-1)*ny*nx + (j-5)*nx + (i-1)
            index_jm3 = (k-1)*ny*nx + (j-4)*nx + (i-1)
            index_jm2 = (k-1)*ny*nx + (j-3)*nx + (i-1)
            index_jm1 = (k-1)*ny*nx + (j-2)*nx + (i-1)
            index_jp1 = (k-1)*ny*nx + (j  )*nx + (i-1)
            index_jp2 = (k-1)*ny*nx + (j+1)*nx + (i-1)
            index_jp3 = (k-1)*ny*nx + (j+2)*nx + (i-1)
            index_jp4 = (k-1)*ny*nx + (j+3)*nx + (i-1)
            !
            index_km4 = (k-5)*ny*nx + (j-1)*nx + (i-1)
            index_km3 = (k-4)*ny*nx + (j-1)*nx + (i-1)
            index_km2 = (k-3)*ny*nx + (j-1)*nx + (i-1)
            index_km1 = (k-2)*ny*nx + (j-1)*nx + (i-1)
            index_kp1 = (k  )*ny*nx + (j-1)*nx + (i-1)
            index_kp2 = (k+1)*ny*nx + (j-1)*nx + (i-1)
            index_kp3 = (k+2)*ny*nx + (j-1)*nx + (i-1)
            index_kp4 = (k+3)*ny*nx + (j-1)*nx + (i-1)
            !
            coef_ip1 =  0.d0
            coef_ip2 =  0.d0
            coef_ip3 =  0.d0
            coef_ip4 =  0.d0
            if (i .lt. nx-1) then
               if (NINT(q(gr_mask+index_ip1)).eq.BHMASK_UNMASKED .and.
     *             NINT(q(gr_mask+index_ip2)).eq.BHMASK_UNMASKED ) then
                  if (int_order.eq.1) then
                     coef_ip1 =  1.d0
                  else if (int_order.eq.2) then
                     coef_ip1 =  2.d0
                     coef_ip2 = -1.d0
                  else if (int_order.eq.3) then
                     coef_ip1 =  3.d0
                     coef_ip2 = -3.d0
                     coef_ip3 =  1.d0
	         end if
               end if
            end if
            !
            coef_im1 =  0.d0
            coef_im2 =  0.d0
            coef_im3 =  0.d0
            coef_im4 =  0.d0
            if (i .gt. 2) then
               if (NINT(q(gr_mask+index_im1)).eq.BHMASK_UNMASKED .and.
     *             NINT(q(gr_mask+index_im2)).eq.BHMASK_UNMASKED ) then
                  if (int_order.eq.1) then
                     coef_im1 =  1.d0
                  else if (int_order.eq.2) then
                     coef_im1 =  2.d0
                     coef_im2 = -1.d0
                  else if (int_order.eq.3) then
                     coef_im1 =  3.d0
                     coef_im2 = -3.d0
                     coef_im3 =  1.d0
		 end if
	       end if
            end if
            !
            coef_jp1 =  0.d0
            coef_jp2 =  0.d0
            coef_jp3 =  0.d0
            coef_jp4 =  0.d0
            if (j .lt. ny-1) then
               if (NINT(q(gr_mask+index_jp1)).eq.BHMASK_UNMASKED .and.
     *             NINT(q(gr_mask+index_jp2)).eq.BHMASK_UNMASKED ) then
                  if (int_order.eq.1) then
                     coef_jp1 =  1.d0
                  else if (int_order.eq.2) then
                     coef_jp1 =  2.d0
                     coef_jp2 = -1.d0
                  else if (int_order.eq.3) then
                     coef_jp1 =  3.d0
                     coef_jp2 = -3.d0
                     coef_jp3 =  1.d0
                  end if
               end if
            end if
            !
            coef_jm1 =  0.d0
            coef_jm2 =  0.d0
            coef_jm3 =  0.d0
            coef_jm4 =  0.d0
            if (j .gt. 2) then
               if (NINT(q(gr_mask+index_jm1)).eq.BHMASK_UNMASKED .and.
     *             NINT(q(gr_mask+index_jm2)).eq.BHMASK_UNMASKED ) then
                  if (int_order.eq.1) then
                     coef_jm1 =  1.d0
                  else if (int_order.eq.2) then
                     coef_jm1 =  2.d0
                     coef_jm2 = -1.d0
                  else if (int_order.eq.3) then
                     coef_jm1 =  3.d0
                     coef_jm2 = -3.d0
                     coef_jm3 =  1.d0
                  end if
               end if
            end if
            !
            coef_kp1 =  0.d0
            coef_kp2 =  0.d0
            coef_kp3 =  0.d0
            coef_kp4 =  0.d0
            if ( k .lt. nz-1) then
               if (NINT(q(gr_mask+index_kp1)).eq.BHMASK_UNMASKED .and.
     *             NINT(q(gr_mask+index_kp2)).eq.BHMASK_UNMASKED ) then
                  if (int_order.eq.1) then
                     coef_kp1 =  1.d0
                  else if (int_order.eq.2) then
                     coef_kp1 =  2.d0
                     coef_kp2 = -1.d0
                  else if (int_order.eq.3) then
                     coef_kp1 =  3.d0
                     coef_kp2 = -3.d0
                     coef_kp3 =  1.d0
                  end if
               end if
            end if
            !
            coef_km1 =  0.d0
            coef_km2 =  0.d0
            coef_km3 =  0.d0
            coef_km4 =  0.d0
            if (k .gt. 2) then
               if (NINT(q(gr_mask+index_km1)).eq.BHMASK_UNMASKED .and.
     *             NINT(q(gr_mask+index_km2)).eq.BHMASK_UNMASKED ) then
                  if (int_order.eq.1) then
                     coef_km1 =  1.d0
                  else if (int_order.eq.2) then
                     coef_km1 =  2.d0
                     coef_km2 = -1.d0
                  else if (int_order.eq.3) then
                     coef_km1 =  3.d0
                     coef_km2 = -3.d0
                     coef_km3 =  1.d0
                  end if
               end if
            end if
            !
            sum_coef = coef_im1+coef_im2+coef_ip1+coef_ip2
     *		      +coef_im3+coef_im4+coef_ip3+coef_ip4
     *                +coef_jm1+coef_jm2+coef_jp1+coef_jp2
     *                +coef_jm3+coef_jm4+coef_jp3+coef_jp4
     *                +coef_km1+coef_km2+coef_kp1+coef_kp2
     *                +coef_km3+coef_km4+coef_kp3+coef_kp4

            !
            if (NINT(sum_coef).eq.0) then
               if (ltrace2)
     *         write(*,*) 'grid_repopulate: Cannot extrapolate:',i,j,k
               num_unrepop = num_unrepop
               ! Let us just keep this point as masked:
               q(gr_mask+index) = 1.d0 * BHMASK_MASKED
               do gf = 1, num_gfuncs
                  if ( gfunc_name(gf) .ne. 'mask'   .and.
     *                 gfunc_name(gf) .ne. 'cctk_x' .and.
     *                 gfunc_name(gf) .ne. 'cctk_y' .and.
     *                 gfunc_name(gf) .ne. 'cctk_z' .and.
     *                 gfunc_name(gf) .ne. 'xphys'  .and.
     *                 gfunc_name(gf) .ne. 'yphys'  .and.
     *                 gfunc_name(gf) .ne. 'zphys'  .and.
     *                 gfunc_name(gf) .ne. 'r'      .and.
     *                 gfunc_name(gf) .ne. 'chr'        ) then
                     q(gfunc_pointer(gf,gi)+index) = 0.d0
                  end if
               end do
            else
               numrepop = numrepop + 1
               do gf = 1, num_gfuncs
                  if ( gfunc_name(gf) .ne. 'mask' .and.
     *                 gfunc_name(gf) .ne. 'cctk_x' .and.
     *                 gfunc_name(gf) .ne. 'cctk_y' .and.
     *                 gfunc_name(gf) .ne. 'cctk_z' .and.
     *                 gfunc_name(gf) .ne. 'xphys'  .and.
     *                 gfunc_name(gf) .ne. 'yphys'  .and.
     *                 gfunc_name(gf) .ne. 'zphys'  .and.
     *                 gfunc_name(gf) .ne. 'r'      .and.
     *                 gfunc_name(gf) .ne. 'chr'        ) then
                     !
                     ! Do not reference memory that does not exist:
                     !
                     fptr            = gfunc_pointer(gf,gi)
                     q(fptr + index) = 0.d0
                     !
                     if ( NINT(coef_ip1) .ne. 0)
     *                  q(fptr+index) = q(fptr+index) + 
     *       coef_ip1*q(fptr + index_ip1) + coef_ip2*q(fptr + index_ip2)
     *     + coef_ip3*q(fptr + index_ip3) 

                     !
                     if ( NINT(coef_im1) .ne. 0)
     *                  q(fptr+index) = q(fptr+index) + 
     *       coef_im1*q(fptr + index_im1) + coef_im2*q(fptr + index_im2)
     *     + coef_im3*q(fptr + index_im3) 

                     !
                     if ( NINT(coef_jp1) .ne. 0)
     *                  q(fptr+index) = q(fptr+index) + 
     *       coef_jp1*q(fptr + index_jp1) + coef_jp2*q(fptr + index_jp2)
     *     + coef_jp3*q(fptr + index_jp3) 

                     !
                     if ( NINT(coef_jm1) .ne. 0)
     *                  q(fptr+index) = q(fptr+index) + 
     *       coef_jm1*q(fptr + index_jm1) + coef_jm2*q(fptr + index_jm2)
     *     + coef_jm3*q(fptr + index_jm3) 

                     !
                     if ( NINT(coef_kp1) .ne. 0)
     *                  q(fptr+index) = q(fptr+index) + 
     *       coef_kp1*q(fptr + index_kp1) + coef_kp2*q(fptr + index_kp2)
     *     + coef_kp3*q(fptr + index_kp3) 

                     !
                     if ( NINT(coef_km1) .ne. 0)
     *                  q(fptr+index) = q(fptr+index) + 
     *       coef_km1*q(fptr + index_km1) + coef_km2*q(fptr + index_km2)
     *     + coef_km3*q(fptr + index_km3) 

                     !
                     ! Normalize properly:
                     !
                     q(fptr + index) = q(fptr + index) / sum_coef
                     !
                  end if
               end do
            end if
         end if
      end do
      end do
      end do

      if (ltrace) then
        write(*,*)'grid_repopulate: Done w/ repop, now masking'
        if (num_unrepop.gt.0) then
           write(*,*)'grid_repopulate: There were points that could not'
           write(*,*)'grid_repopulate: be repopulated = ',num_unrepop
        end if
        write(*,99)myid,'grid_repopulate: numasked = ',numasked
        write(*,99)myid,'grid_repopulate: numrepop = ',numrepop
      end if

      numasked = 0
      numrepop = 0
      !
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         index = (k-1)*ny*nx + (j-1)*nx + (i-1)
         mask  = NINT( q( gr_mask+index ) )
         if (ltrace2) write(*,*) i,j,k, mask
         !
         if (mask .eq. BHMASK_REPOPULATE) then
            ! This point was masked with the old mask,
            ! but is no longer masked so we need to
            ! to set the point as no longer being masked:
            if (ltrace2)write(*,*)'grid_repopulate:     unmasked:',i,j,k
            q(gr_mask+index) = 1.d0 * BHMASK_UNMASKED
            numrepop = numrepop + 1
         else if (mask .lt. 0) then
            !  This point is masked with the new mask,
            !  and so let us now set the mask to -1
            !  and mask out the fields in case they have not been already
            if (ltrace2)write(*,*)'grid_repopulate: Newly masked:',i,j,k
            do gf = 1, num_gfuncs
               if ( gfunc_name(gf) .ne. 'mask' .and.
     *                 gfunc_name(gf) .ne. 'cctk_x' .and.
     *                 gfunc_name(gf) .ne. 'cctk_y' .and.
     *                 gfunc_name(gf) .ne. 'cctk_z' .and.
     *                 gfunc_name(gf) .ne. 'xphys'  .and.
     *                 gfunc_name(gf) .ne. 'yphys'  .and.
     *                 gfunc_name(gf) .ne. 'zphys'  .and.
     *                 gfunc_name(gf) .ne. 'r'      .and.
     *              gfunc_name(gf) .ne. 'chr'        ) then
                  q(gfunc_pointer(gf,gi)+index) = 0.d0
               end if
            end do
            numasked = numasked + 1
         !else
            !  This point is not masked and was not...
            !  nothing need be done
         end if
      end do
      end do
      end do

      if (ltrace) then
        write(*,99)myid,'grid_repopulate: numasked = ',numasked
        write(*,99)myid,'grid_repopulate: numrepop = ',numrepop
        write(*,99)myid,'grid_repopulate: Done.'
      end if

 98   format('[',I3,'] ',A,3F15.9)
 99   format('[',I3,'] ',A,3I7)
      return
      end       ! END: grid_repopulate

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_output_sdf:                                                          cc
cc                   Output grid functions to SDF files.                      cc
cc                Output dependent on "output_style":                         cc
cc                     1 --- grids on level X output to    chiX.sdf           cc
cc                     4 --- output all grids at and finer than level to chiX cc
cc                     7 --- output grid X to file chiX                       cc
cc                                                                            cc
cc      NB: For MPI version, only allow option 4                              cc
cc      NB: Adding back in option 1                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_output_sdf(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'glob.inc'
      include    'output.inc'
      include    'grid_methods.inc'
      include    'largesmall.inc'
      integer     ii, jj, kk, index, indexn
      integer     level,      i,  i_start,  i_finish, j,
     *            taglength, namelength, mystringlength,
     *            my_level,  grid_count, lev_diff, count
      logical     no_output, this_level_off
      real(kind=8)      grid_time, level_time, time_diff
      real(kind=8)      bbox(6), bboxFV(6)
      integer           gft_shapeFV(3)
      ! Needed for FV output:
      real(kind=8)      half_hi
      real(kind=8)      y_mid, z_mid
      integer           i_mid, j_mid, k_mid
      integer           ilb, jlb, klb, iub,jub,kub
      integer           nx, ny, nz, nxn,nyn,nzn, rank
      character(2) lev
      character(5) gfunc_name_5char
      integer      num_families, fam, fam_type, fam_level,
     *             fam_period, fam_digit, myid
      real(kind=8) fam_lb1, fam_lb2, fam_lb3,
     *             fam_ub1, fam_ub2, fam_ub3
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )
      logical     ltraceF                ! for "family" output only
      parameter ( ltraceF = .false. )
      real(kind=8) myl2norm3d
      external     myl2norm3d
      integer      proc_return_myid
      external     proc_return_myid
      integer      find_index, find_indexN
      external     find_index, find_indexN
      


 97   format('[',I3,'] ',A,A)
 98   format('[',I3,'] ',A,3F15.9)
 99   format('[',I3,'] ',A,3I7)
      myid = proc_return_myid()

      if (gridnum.le.0) then
         write(*,*) 'grid_output_sdf: nonexistent grid: ',gridnum
         write(*,*) 'grid_output_sdf: returning'
         return
      end if

      level      = grid_return_level(gridnum)
      grid_time  = grid_return_time(gridnum)
      grid_count = grid_return_count(gridnum)

      !
      ! This is the level of the grid
      !
      my_level  = level
      count     = lev_count(my_level)

      if (ltrace) then
         write(*,*) '    %%%%%%%%%%%%%%%%%%%%%%%%%%%'
         write(*,*) '    grid_output_sdf, gridnum = ',gridnum
         write(*,*) '                  my_level   = ',my_level
         write(*,*) '                  grid_time  = ',grid_time
         write(*,*) '                output_style = ',output_style
      end if

      call load_pointers(gridnum)

      gft_shape(1) = gr_nx(gridnum)
      gft_shape(2) = gr_ny(gridnum)
      gft_shape(3) = gr_nz(gridnum)

      nx = gft_shape(1)
      ny = gft_shape(2)
      nz = gft_shape(3)

      bbox(1) = gr_minx(gridnum)
      bbox(2) = gr_maxx(gridnum)
      bbox(3) = gr_miny(gridnum)
      bbox(4) = gr_maxy(gridnum)
      bbox(5) = gr_minz(gridnum)
      bbox(6) = gr_maxz(gridnum)

      !
      ! For FV fields:
      !     (1) offset the bbox to place points at the cells:
      !     (2) Only (1,nx-1)(1,ny-1)(1,nz-1) cells are meainingful
      !
      half_hi = 0.5d0 * gr_h(gridnum)
      !write(*,*) 'gridnum = ',gridnum
      !write(*,*) 'half_hi = ',half_hi
      do i = 1, 6
         bboxFV(i) = bbox(i) + half_hi
         ! Upper bounds decrease by half the grid spacing:
         half_hi = -1.d0*half_hi
         !write(*,*) 'bboxFV(i) = ',i,bboxFV(i)
      end do
      gft_shapeFV(1) = gr_nx(gridnum) - 1
      gft_shapeFV(2) = gr_ny(gridnum) - 1
      gft_shapeFV(3) = gr_nz(gridnum) - 1

      taglength               = mystringlength(tag)
      gft_name(1:taglength)   = tag(1:taglength)
      !write(*,*) 'gft_name = ',gft_name

      if (ltrace2) then
         do i = 1, 6
            write(*,*) 'bbox(i) = ',i,bbox(i)
         end do
         write(*,*) 'gft_name = ',gft_name
      end if

      !
      ! For cases: outstyle == 4 or 6
      ! Then we may need to output the grid even if output
      ! for this level is turned off (e.g. if output for level 0
      ! is turned on, then every fine grid is output when it
      ! is evolved to the time of the coarse grid). In fact,
      ! the grid is likely to have to be output many times
      ! for each (coarser) grid that is being output.
      !
      !    NB: A problem for the MPI version is that
      !    only the master processor keeps track of the lev_count.
      !    Also, the grid counter only gets updated on the owner proc.
      !
      !
  11     continue
         !
         ! Output on level "level"?
         !
         !   (1) Has user selected output on this level?
         !
         no_output = out_period(level).eq.0
         !
         !   (2) Is it time to output at this level?
         !
         if (.not.no_output) no_output =
     *           mod(lev_count(level), out_period(level)) .ne.0

         !
         !   (3) Is this subgrid advanced to the same
         !       "time" as this level?                  
         !
         !
         if (output_style.eq.4 .and. .not.no_output)
     *      no_output = .not.level_time_aligned(my_level,level)


         if (.not.no_output .and. ltrace2) then
         write(*,*)'                     level = ',level
         write(*,*)'                     count = ',count
         write(*,*)'                level_time = ',level_time
         write(*,*)'                 time_diff = ',time_diff
         write(*,*)'         out_period(level) = ',out_period(level)
         write(*,*)'          lev_count(level) = ',lev_count(level)
         end if

         if (no_output) then
            level = level - 1
            if (level .ge. 0 .and. output_style.eq.4) then
               goto 11
            else 
               if (ltrace) write(*,*) '  returning'
               !return
               goto 88
            end if
         end if

      !
      ! Convert integer: level
      !       to string: lev   
      !
      if (output_style .eq. 7) then
         call int2str(gridnum,lev)
      else
         call int2str(level,lev)
      end if

      if (ltrace) then
         write(*,*) 'Level at which output is turned on:',level
         write(*,*) '                  outputting sdfs'
         write(*,*) '                     at time ',gr_t(Levelp(level))
      end if
      if (ltrace2) then
         write(*,*) '    output: grid          = ',gridnum
         write(*,*) '    output: level of grid = ',my_level, lev
         write(*,*) '                     time = ',grid_time
      end if

      !
      ! For dimensionally reduced output:
      !     data is contiguous in x, so pick start (ii=1)
      !     for others, find data corresponding to the
      !     location of the middle of the coarse grid
      !
      y_mid  = 0.5d0*(miny0 + maxy0)
      z_mid  = 0.5d0*(minz0 + maxz0)
      j_mid  = find_index(q(gr_y(gridnum)), y_mid, gr_ny(gridnum))
      k_mid  = find_index(q(gr_z(gridnum)), z_mid, gr_nz(gridnum))

      do j = 1, num_gfuncs
         if (gfunc_out(j).gt.0) then
               gft_name(taglength+1:)  = gfunc_name(j)
               namelength              = mystringlength(gft_name)
               gft_name(namelength+1:) = lev
               if (output_dim .ge. 3 .and. output_dim .le. 6) then
                   !
                   ! Output 3D:
                   !
                   if ( field_is_FV(gfunc_type(j)) ) then
                      ! FV field:
                      call copyFV(q(gr_tmp),q(gfunc_pointer(j,gridnum)),
     *                     gr_nx(gridnum),gr_ny(gridnum),gr_nz(gridnum))
#ifdef AIX
                      gft_rc       =
     *                gft_out_bbox(  trim(gft_name)/ /CHAR(0), 
     *                               %VAL(grid_time),   gft_shapeFV,
     *                               %VAL(3), bboxFV,
     *                               q(gr_tmp) )
#else
                      gft_rc       =
     *                gft_out_bbox(  gft_name, grid_time,   gft_shapeFV,
     *                               3, bboxFV,
     *                               q(gr_tmp) )
#endif
                   else
                      ! FD field:
#ifdef AIX
                      gft_rc       =
     *                gft_out_bbox(trim(gft_name)/ /CHAR(0),
     *                             %VAL(grid_time),
     *                             gft_shape,
     *                             %VAL(3), bbox,
     *                             q(gfunc_pointer(j,gridnum)))
#else
                      gft_rc       =
     *                gft_out_bbox(  gft_name, grid_time, gft_shape,
     *                               3, bbox,
     *                               q(gfunc_pointer(j,gridnum)))
#endif
                  end if
               end if
               if (     (output_dim .eq. 2 .or. output_dim .eq. 5 .or.
     *                   output_dim .eq. 6 .or. output_dim .eq. 7 )
     *             .and. k_mid .gt. 0) then
                   !
                   ! Output 2D:
                   !
                   ii    = 1
                   jj    = 1
                   gft_name(namelength+1:) = '2D'//lev
                   if ( field_is_FV(gfunc_type(j)) ) then
                      ! FV field:
                      index =
     *                   (k_mid-1)*(gr_ny(gridnum)-1)*(gr_nx(gridnum)-1)
     *                 + (jj-1)*(gr_nx(gridnum)-1)
     *                 + (ii-1)
                      call copyFV(q(gr_tmp),q(gfunc_pointer(j,gridnum)),
     *                     gr_nx(gridnum),gr_ny(gridnum),gr_nz(gridnum))
#ifdef AIX
                      gft_rc       =
     *                gft_out_bbox(  trim(gft_name)/ /CHAR(0), 
     *                               %VAL(grid_time),
     *                               gft_shapeFV,
     *                               %VAL(2), bboxFV,
     *                               q(gr_tmp+index) )
#else
                      gft_rc       =
     *                gft_out_bbox(  gft_name, grid_time,   gft_shapeFV,
     *                               2, bboxFV,
     *                               q(gr_tmp+index) )
#endif
                   else
                      ! FD field:
                      index =    (k_mid-1)*gr_ny(gridnum)*gr_nx(gridnum)
     *                         + (jj-1)*gr_nx(gridnum)
     *                         + (ii-1)
#ifdef AIX
                      gft_rc       =
     *                gft_out_bbox(trim(gft_name)/ /CHAR(0),
     *                               %VAL(grid_time),
     *                               gft_shape,%VAL(2), bbox,
     *                               q(gfunc_pointer(j,gridnum)+index))
#else
                      gft_rc       =
     *                gft_out_bbox(  gft_name, grid_time,     gft_shape,
     *                               2, bbox,
     *                               q(gfunc_pointer(j,gridnum)+index))
#endif
                  end if
               end if
               if (     (output_dim .eq. 1 .or. output_dim .eq. 4 .or.
     *                   output_dim .eq. 6 .or. output_dim .eq. 7 )
     *             .and. k_mid .gt. 0 .and. j_mid .gt. 0) then
                   !
                   ! Output 1D:
                   !
                   ii    = 1
                   index =    (k_mid-1)*gr_ny(gridnum)*gr_nx(gridnum)
     *                      + (j_mid-1)*gr_nx(gridnum)
     *                      + (ii-1)
                   gft_name(namelength+1:) = '1D'//lev
                   !write(*,*) 'NAME:   ', gft_name
#ifdef AIX
                   gft_rc       =
     *             gft_out_bbox(trim(gft_name)/ /CHAR(0), 
     *                          %VAL(grid_time), gft_shape,
     *                          %VAL(1), bbox,
     *                            q(gfunc_pointer(j,gridnum)+index))
#else
                   gft_rc       =
     *             gft_out_bbox(  gft_name, grid_time,     gft_shape,
     *                            1, bbox,
     *                            q(gfunc_pointer(j,gridnum)+index))
#endif
               end if
         end if
      end do

      !
      ! See if time to output coarser grids
      !
      if (output_style.eq.4) then
         level = level - 1
         if (level.ge.0) goto 11
      end if

 88   continue

      ! --- BEGIN: Family output  ----
      ! This is just the number families hard-coded in the par setup file:
      num_families = 4
      do fam = 1, num_families
         if (fam.eq.1) then
            fam_type    = output_f1_type
            fam_level   = output_f1_level
            fam_period  = output_f1_period
            fam_lb1     = output_f1_lb1
            fam_lb2     = output_f1_lb2
            fam_lb3     = output_f1_lb3
            fam_ub1     = output_f1_ub1
            fam_ub2     = output_f1_ub2
            fam_ub3     = output_f1_ub3
            lev         = 'F1'
            fam_digit   = 10
         else if (fam.eq.2) then
            fam_type    = output_f2_type
            fam_level   = output_f2_level
            fam_period  = output_f2_period
            fam_lb1     = output_f2_lb1
            fam_lb2     = output_f2_lb2
            fam_lb3     = output_f2_lb3
            fam_ub1     = output_f2_ub1
            fam_ub2     = output_f2_ub2
            fam_ub3     = output_f2_ub3
            lev         = 'F2'
            fam_digit   = 100
         else if (fam.eq.3) then
            fam_type    = output_f3_type
            fam_level   = output_f3_level
            fam_period  = output_f3_period
            fam_lb1     = output_f3_lb1
            fam_lb2     = output_f3_lb2
            fam_lb3     = output_f3_lb3
            fam_ub1     = output_f3_ub1
            fam_ub2     = output_f3_ub2
            fam_ub3     = output_f3_ub3
            lev         = 'F3'
            fam_digit   = 1000
         else if (fam.eq.4) then
            fam_type    = output_f4_type
            fam_level   = output_f4_level
            fam_period  = output_f4_period
            fam_lb1     = output_f4_lb1
            fam_lb2     = output_f4_lb2
            fam_lb3     = output_f4_lb3
            fam_ub1     = output_f4_ub1
            fam_ub2     = output_f4_ub2
            fam_ub3     = output_f4_ub3
            lev         = 'F4'
            fam_digit   = 10000
         end if
         if (ltraceF) then
            write(*,99)myid,'grid_output_sdf: Family: ',fam
            write(*,99)myid,'grid_output_sdf: my_level',my_level
            write(*,99)myid,'grid_output_sdf: type:   ',fam_type
            write(*,99)myid,'grid_output_sdf: flevel: ',fam_level
            write(*,99)myid,'grid_output_sdf: period: ',fam_period
            write(*,99)myid,'grid_output_sdf: count:  ',
     *           lev_count(fam_level)
            write(*,98)myid,'grid_output_sdf: lb1/2/3: ',fam_lb1,
     *           fam_lb2, fam_lb3
            write(*,98)myid,'grid_output_sdf: ub1/2/3: ',fam_ub1,
     *           fam_ub2, fam_ub3
         end if
         !
         no_output = my_level .lt. fam_level
         if(ltraceF)write(*,99)myid,'grid_output_sdfno_output',no_output
         !
         if (.not. no_output) no_output = fam_period .le. 0
         if(ltraceF)write(*,99)myid,'grid_output_sdfno_output',no_output
         !
         if (.not. no_output) no_output = 
     *        mod(lev_count(fam_level),fam_period) .ne. 0
         if(ltraceF)write(*,99)myid,'grid_output_sdfno_output',no_output
         !
         if (.not. no_output) no_output = 
     *       .not.level_time_aligned(my_level,fam_level)
         if(ltraceF)write(*,99)myid,'grid_output_sdfno_output',no_output
         !
         if (.not. no_output) then
            if(ltraceF)write(*,99)myid,'grid_output_sdf: Time to output'
            if (fam_type .eq. 3) then
               !
               ! 3D output over the region (lb1,lb2,lb3)-->(ub1,ub2,ub3)
               !
               ilb  = find_indexN(q(gr_x(gridnum)), fam_lb1,nx)
               jlb  = find_indexN(q(gr_y(gridnum)), fam_lb2,ny)
               klb  = find_indexN(q(gr_z(gridnum)), fam_lb3,nz)
               iub  = find_indexN(q(gr_x(gridnum)), fam_ub1,nx)
               jub  = find_indexN(q(gr_y(gridnum)), fam_ub2,ny)
               kub  = find_indexN(q(gr_z(gridnum)), fam_ub3,nz)
               if (ltraceF) then
                  write(*,99)myid,'grid_output_sdf:i/j/klb=',ilb,jlb,klb
                  write(*,99)myid,'grid_output_sdf:i/j/kub=',iub,jub,kub
               end if
               if(ilb.lt. 1) ilb=1 
               if(jlb.lt. 1) jlb=1 
               if(klb.lt. 1) klb=1 
               if(iub.gt.nx) iub=nx
               if(jub.gt.ny) jub=ny
               if(kub.gt.nz) kub=nz
               if (ltraceF) then
                  write(*,99)myid,'grid_output_sdf: 3D Subregion'
                  write(*,98)myid,'grid_:x(il)=',q(gr_x(gridnum)+ilb-1)
                  write(*,98)myid,'grid_:y(jl)=',q(gr_y(gridnum)+jlb-1)
                  write(*,98)myid,'grid_:z(kl)=',q(gr_z(gridnum)+klb-1)
                  write(*,98)myid,'grid_:x(iu)=',q(gr_x(gridnum)+iub-1)
                  write(*,98)myid,'grid_:y(ju)=',q(gr_y(gridnum)+jub-1)
                  write(*,98)myid,'grid_:z(ku)=',q(gr_z(gridnum)+kub-1)
               end if
               nxn          = iub - ilb + 1
               nyn          = jub - jlb + 1
               nzn          = kub - klb + 1
               if (ltraceF) then
                  write(*,99)myid,'grid_output_sdf:nx/y/zn=',nxn,nyn,nzn
               end if
               rank         = 3
               if (nxn.le.1.or.nyn.le.1.or.nzn.le.1) goto 89
               gft_shape(1) = nxn
               gft_shape(2) = nyn
               gft_shape(3) = nzn
               bbox(1)      = q(gr_x(gridnum)+ilb-1)
               bbox(2)      = q(gr_x(gridnum)+iub-1)
               bbox(3)      = q(gr_y(gridnum)+jlb-1)
               bbox(4)      = q(gr_y(gridnum)+jub-1)
               bbox(5)      = q(gr_z(gridnum)+klb-1)
               bbox(6)      = q(gr_z(gridnum)+kub-1)
            else if (fam_type .eq. 11) then
               !
               ! 1D output along x with constant y=lb2 and z=lb3
               !
               if (gr_miny(gridnum) .gt. fam_lb2+SMALLNUMBER .or.
     *             gr_maxy(gridnum) .lt. fam_lb2-SMALLNUMBER ) goto 89
               if (gr_minz(gridnum) .gt. fam_lb3+SMALLNUMBER .or.
     *             gr_maxz(gridnum) .lt. fam_lb3-SMALLNUMBER ) goto 89
               jlb  = find_index(q(gr_y(gridnum)), fam_lb2,ny)
               klb  = find_index(q(gr_z(gridnum)), fam_lb3,nz)
               if (jlb.le.0 .or. klb.le.0) goto 89
               !jlb  = find_indexN(q(gr_y(gridnum)), fam_lb2,ny)
               !klb  = find_indexN(q(gr_z(gridnum)), fam_lb3,nz)
               if (ltraceF) then
                  write(*,98)myid,'grid_output_sdf: Constant Y&Z output'
                  write(*,98)myid,'grid_:y(i)= ',q(gr_y(gridnum)+jlb-1)
                  write(*,98)myid,'grid_:z(i)= ',q(gr_z(gridnum)+klb-1)
               end if
               nxn          = nx
               gft_shape(1) = nxn
               rank         = 1
               bbox(1)      = gr_minx(gridnum)
               bbox(2)      = gr_maxx(gridnum)
            else if (fam_type .eq. 12) then
               !
               ! 1D output along x with constant x=lb1 and z=lb3
               !
               if (gr_minx(gridnum) .gt. fam_lb1+SMALLNUMBER .or.
     *             gr_maxx(gridnum) .lt. fam_lb1-SMALLNUMBER ) goto 89
               if (gr_minz(gridnum) .gt. fam_lb3+SMALLNUMBER .or.
     *             gr_maxz(gridnum) .lt. fam_lb3-SMALLNUMBER ) goto 89
               ilb  = find_index(q(gr_x(gridnum)), fam_lb1,nx)
               klb  = find_index(q(gr_z(gridnum)), fam_lb3,nz)
               if (ilb.le.0 .or. klb.le.0) goto 89
               if (ltraceF) then
                  write(*,98)myid,'grid_output_sdf: Constant X&Z output'
                  write(*,98)myid,'grid_:x(i)= ',q(gr_x(gridnum)+ilb-1)
                  write(*,98)myid,'grid_:z(i)= ',q(gr_z(gridnum)+klb-1)
               end if
               nyn          = ny
               gft_shape(1) = nyn
               rank         = 1
               bbox(1)      = gr_miny(gridnum)
               bbox(2)      = gr_maxy(gridnum)
            else if (fam_type .eq. 13) then
               !
               ! 1D output along x with constant x=lb1 and y=lb2
               !
               if (gr_minx(gridnum) .gt. fam_lb1+SMALLNUMBER .or.
     *             gr_maxx(gridnum) .lt. fam_lb1-SMALLNUMBER ) goto 89
               if (gr_miny(gridnum) .gt. fam_lb2+SMALLNUMBER .or.
     *             gr_maxy(gridnum) .lt. fam_lb2-SMALLNUMBER ) goto 89
               ilb  = find_index(q(gr_x(gridnum)), fam_lb1,nx)
               jlb  = find_index(q(gr_y(gridnum)), fam_lb2,ny)
               if (ilb.le.0 .or. jlb.le.0) goto 89
               if (ltraceF) then
                  write(*,98)myid,'grid_output_sdf: Constant X&Y output'
                  write(*,98)myid,'grid_:x(i)= ',q(gr_x(gridnum)+ilb-1)
                  write(*,98)myid,'grid_:y(i)= ',q(gr_y(gridnum)+jlb-1)
               end if
               nzn          = nz
               gft_shape(1) = nzn
               rank         = 1
               bbox(1)      = gr_minz(gridnum)
               bbox(2)      = gr_maxz(gridnum)
            else if (fam_type .eq. 21) then
               !
               ! 2D output over the plane defined by x=lb1
               !
               if (gr_minx(gridnum) .gt. fam_lb1+SMALLNUMBER .or.
     *             gr_maxx(gridnum) .lt. fam_lb1-SMALLNUMBER ) goto 89
               ilb  = find_index(q(gr_x(gridnum)), fam_lb1,nx)
               if (ilb.le.0) goto 89
               if (ltraceF) then
                  write(*,98)myid,'grid_output_sdf: Constant X output'
                  write(*,98)myid,'grid_:x(i)= ',q(gr_x(gridnum)+ilb-1)
               end if
               nyn          = ny
               nzn          = nz
               gft_shape(1) = nyn
               gft_shape(2) = nzn
               rank         = 2
               bbox(1)      = gr_miny(gridnum)
               bbox(2)      = gr_maxy(gridnum)
               bbox(3)      = gr_minz(gridnum)
               bbox(4)      = gr_maxz(gridnum)
            else if (fam_type .eq. 22) then
               !
               ! 2D output over the plane defined by y=lb2
               !
               if (gr_miny(gridnum) .gt. fam_lb2+SMALLNUMBER .or.
     *             gr_maxy(gridnum) .lt. fam_lb2-SMALLNUMBER ) goto 89
               jlb  = find_index(q(gr_y(gridnum)), fam_lb2,ny)
               if (jlb.le.0) goto 89
               if (ltraceF) then
                  write(*,98)myid,'grid_output_sdf: Constant Y output'
                  write(*,98)myid,'grid_:y(i)= ',q(gr_y(gridnum)+jlb-1)
               end if
               nxn          = nx
               nzn          = nz
               gft_shape(1) = nxn
               gft_shape(2) = nzn
               rank         = 2
               bbox(1)      = gr_minx(gridnum)
               bbox(2)      = gr_maxx(gridnum)
               bbox(3)      = gr_minz(gridnum)
               bbox(4)      = gr_maxz(gridnum)
            else if (fam_type .eq. 23) then
               !
               ! 2D output over the plane defined by z=lb3
               !
               if (gr_minz(gridnum) .gt. fam_lb3+SMALLNUMBER .or.
     *             gr_maxz(gridnum) .lt. fam_lb3-SMALLNUMBER ) goto 89
               klb  = find_index(q(gr_z(gridnum)), fam_lb3,nz)
               if (klb.le.0) goto 89
               if (ltraceF) then
                  write(*,98)myid,'grid_output_sdf: Constant Z output'
                  write(*,98)myid,'grid_:z(i)= ',q(gr_z(gridnum)+klb-1)
               end if
               nxn          = nx
               nyn          = ny
               gft_shape(1) = nxn
               gft_shape(2) = nyn
               rank         = 2
               bbox(1)      = gr_minx(gridnum)
               bbox(2)      = gr_maxx(gridnum)
               bbox(3)      = gr_miny(gridnum)
               bbox(4)      = gr_maxy(gridnum)
            else
               write(*,99)myid,'grid_output_sdf: Unrecognized or',
     *                 ' unimplemented output type: ',fam_type
               goto 89
            end if
            !
            ! Loop over fields:
            !
            do j = 1, num_gfuncs
            !if (mod( int(gfunc_out(j)/10), fam_digit) .gt. 0) then
            !write(*,*) 'gfunc_name: = ',gfunc_name(j)
            !write(*,*) 'gfunc_out:  = ',gfunc_out(j)
            !write(*,*) 'fam_digit   = ',fam_digit
            !write(*,*) 'digit: ',mod( int(gfunc_out(j)/fam_digit), 10)
            if(ltraceF)write(*,99)myid,'grid_output_sdf:fam/j: ',fam,j
            if (mod( int(gfunc_out(j)/fam_digit), 10) .gt. 0) then
			      gft_name(taglength+1:)  = gfunc_name(j)
			      namelength              = mystringlength(gft_name)
			      gft_name(namelength+1:) = lev
               if(ltraceF)write(*,97)myid,'grid_:output: ',gft_name(1:7)
               !
               if (fam_type .eq. 3) then
                  do kk = klb, kub
                  do jj = jlb, jub
                  do ii = ilb, iub
                     index  = (kk-1  )*nx *ny  + (jj-1  )*nx  + (ii-1  )
                     indexn = (kk-klb)*nxn*nyn + (jj-jlb)*nxn + (ii-ilb)
                     q(gr_tmp+indexn) =q(gfunc_pointer(j,gridnum)+index)
                  end do
                  end do
                  end do
               else if (fam_type .eq. 11) then
                  jj    = jlb
                  kk    = klb
                  do ii = 1, nxn
                     index  = (kk-1  )*nx *ny  + (jj-1  )*nx  + (ii-1  )
                     indexn = (ii-1)
                     q(gr_tmp+indexn) =q(gfunc_pointer(j,gridnum)+index)
                  end do
               else if (fam_type .eq. 12) then
                  ii    = ilb
                  kk    = klb
                  do jj = 1, nyn
                     index  = (kk-1  )*nx *ny  + (jj-1  )*nx  + (ii-1  )
                     indexn = (jj-1)
                     q(gr_tmp+indexn) =q(gfunc_pointer(j,gridnum)+index)
                  end do
               else if (fam_type .eq. 13) then
                  ii    = ilb
                  jj    = jlb
                  do kk = 1, nzn
                     index  = (kk-1  )*nx *ny  + (jj-1  )*nx  + (ii-1  )
                     indexn = (kk-1)
                     q(gr_tmp+indexn) =q(gfunc_pointer(j,gridnum)+index)
                  end do
               else if (fam_type .eq. 21) then
                  ii    = ilb
                  do kk = 1, nzn
                  do jj = 1, nyn
                     index  = (kk-1  )*nx *ny  + (jj-1  )*nx  + (ii-1  )
                     indexn = (kk-1)*nyn + (jj-1)
                     q(gr_tmp+indexn) =q(gfunc_pointer(j,gridnum)+index)
                  end do
                  end do
               else if (fam_type .eq. 22) then
                  jj    = jlb
                  do kk = 1, nzn
                  do ii = 1, nxn
                     index  = (kk-1  )*nx *ny  + (jj-1  )*nx  + (ii-1  )
                     indexn = (kk-1)*nxn + (ii-1)
                     q(gr_tmp+indexn) =q(gfunc_pointer(j,gridnum)+index)
                  end do
                  end do
               else if (fam_type .eq. 23) then
                  kk    = klb
                  do jj = 1, nyn
                  do ii = 1, nxn
                     index  = (kk-1  )*nx *ny  + (jj-1  )*nx  + (ii-1  )
                     indexn = (jj-1)*nxn + (ii-1)
                     q(gr_tmp+indexn) =q(gfunc_pointer(j,gridnum)+index)
                  end do
                  end do
               end if
               !
               ! Do the output:
               !
               if(ltrace2)write(*,99)myid,'grid_output_sdf:call gft_out'
#ifdef AIX
               gft_rc = gft_out_bbox(trim(gft_name)/ /CHAR(0), 
     *                          %VAL(grid_time), gft_shape,
     *                          %VAL(rank), bbox, q(gr_tmp) )
#else
               gft_rc = gft_out_bbox(  gft_name, grid_time, gft_shape,
     *                               rank, bbox, q(gr_tmp) )
#endif
               if(ltrace2)write(*,99)myid,'grid_output_sdf: output done'
               !
            end if
            end do
         else
            if(ltraceF)write(*,99)myid,'grid_output_sdf: Do not output.'
         end if
 89      continue
         if(ltraceF)write(*,99)myid,'grid_output_sdf: Next family...'
      end do   ! loop over families
      ! --- END:   Family output  ----

      if (ltrace) then
         write(*,99)myid,'grid_output_sdf: Done. Returning'
      end if

      return
      end     ! END: subroutine grid_output_sdf

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  field_is_FV:                                                              cc
cc              Returns .TRUE. if field type (from fields.inc) is any of      cc
cc              the FV types. Else, returns .FALSE.                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function field_is_FV(type)
      implicit    none
      integer     type
      include    'grid.inc'

      if (  
     *     type .eq. GFUNC_INTEGRALFV     .or.
     *     type .eq. GFUNC_INTEGRAL_ADVFV .or.
     *     type .eq. GFUNC_INTEGRAL_DERFV  ) then
          field_is_FV = .true.
      else
          field_is_FV = .false.
      end if

      return
      end        ! function: field_is_FV

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_setbndslarge:                                                        cc
cc                                                                            cc
cc      Debugging routine which sets points on the AMR boundary to            cc
cc      LARGENUMBER before setting them by interpolation to test              cc
cc      that all the points are getting set properly.                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_setbndslarge(gi)
      implicit none
      integer  gi
      include    'grid.inc'
      include    'glob.inc'
      include    'largesmall.inc'
      include    'chr.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer  nx,ny,nz, i,j,k,l, index, mychr
      integer  ii,jj,kk, indexii
      real(kind=8)      myl2norm3d
      external          myl2norm3d
      !integer  myid,    proc_return_myid
      !external          myl2norm3d

      logical     ltrace
      parameter ( ltrace = .false.)

      call load_pointers(gi)

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

      if (ltrace) then
        write(*,*)'grid_setbndslarge: gi     = ',gi
        write(*,*)'grid_setbndslarge: nx/y/z = ',nx,ny,nz
      end if

      ! z-boundaries first:
      if(ltrace)write(*,*)'grid_setbndslarge: Z bounds:'
      do k = 1, nz, nz-1
         do j = 1, ny
         do i = 1, nx
            index = (k-1)*ny*nx+(j-1)*nx+(i-1)
            mychr = NINT(q(gr_chr+index))
            if (mychr .eq. NINT(CHR_AMR_BDY)) then
               ! These pts need to be re-set
               !write(*,*) max(1,k-ghostwidth+1),min(nz,k+ghostwidth-1)
               do kk = max(1,k-ghostwidth+1),min(nz,k+ghostwidth-1)
                  indexii = (kk-1)*ny*nx+(j-1)*nx+(i-1)
                  mychr   = NINT(q(gr_chr+indexii))
                  if (mychr .ne. NINT(CHR_DECO_BDY)) then
                     do l = 1, num_gfuncs
                        if (gfunc_type(l) .eq. GFUNC_INTEGRAL)
                        !if (gfunc_type(l) .eq. GFUNC_INTEGRAL_ADV)
     .                     q(gfunc_pointer(l,gi)+indexii)=LARGENUMBER
                     end do
                  end if
               end do
               !
            end if
         end do
         end do
      end do
      ! y-boundaries first:
      if(ltrace)write(*,*)'grid_setbndslarge: Y bounds:'
      do j = 1, ny, ny-1
         do k = 1, nz
         do i = 1, nx
            index = (k-1)*ny*nx+(j-1)*nx+(i-1)
            mychr = NINT(q(gr_chr+index))
            if (mychr .eq. NINT(CHR_AMR_BDY)) then
               ! These pts need to be re-set
               !write(*,*) max(1,j-ghostwidth+1),min(ny,j+ghostwidth-1)
               do jj = max(1,j-ghostwidth+1),min(ny,j+ghostwidth-1)
                  indexii = (k-1)*ny*nx+(jj-1)*nx+(i-1)
                  mychr   = NINT(q(gr_chr+indexii))
                  if (mychr .ne. NINT(CHR_DECO_BDY)) then
                     do l = 1, num_gfuncs
                        !if (gfunc_type(l) .eq. GFUNC_INTEGRAL_ADV)
                        if (gfunc_type(l) .eq. GFUNC_INTEGRAL)
     .                     q(gfunc_pointer(l,gi)+indexii)=LARGENUMBER
                     end do
                  end if
               end do
               !
            end if
         end do
         end do
      end do
      ! x-boundaries first:
      if(ltrace)write(*,*)'grid_setbndslarge: X bounds:'
      do i = 1, nx, nx-1
         do k = 1, nz
         do j = 1, ny
            index = (k-1)*ny*nx+(j-1)*nx+(i-1)
            mychr = NINT(q(gr_chr+index))
            if (mychr .eq. NINT(CHR_AMR_BDY)) then
               ! These pts need to be re-set
               !write(*,*) max(1,i-ghostwidth+1),min(nx,i+ghostwidth-1)
               do ii = max(1,i-ghostwidth+1),min(nx,i+ghostwidth-1)
                  indexii = (k-1)*ny*nx+(j-1)*nx+(ii-1)
                  mychr   = NINT(q(gr_chr+indexii))
                  if (mychr .ne. NINT(CHR_DECO_BDY)) then
                     do l = 1, num_gfuncs
                        !if (gfunc_type(l) .eq. GFUNC_INTEGRAL_ADV)
                        if (gfunc_type(l) .eq. GFUNC_INTEGRAL)
     .                     q(gfunc_pointer(l,gi)+indexii)=LARGENUMBER
                     end do
                  end if
               end do
               !
            end if
         end do
         end do
      end do

      if (ltrace) then
         call field_out2d(q(gfunc_pointer(1,gi)
     .                      +((nz+1)/2-1)*ny*nx),gr_t(gi),
     *                    'setbndslargeN',gr_minx(gi),gr_maxx(gi),
     *                    gr_miny(gi),gr_maxy(gi),nx,ny,myid)
         call field_out2d(q(gfunc_pointer(2,gi)
     .                      +((nz+1)/2-1)*ny*nx),gr_t(gi),
     *                    'setbndslargeNP1',gr_minx(gi),gr_maxx(gi),
     *                    gr_miny(gi),gr_maxy(gi),nx,ny,myid)

        call grid_dump_info(gi)
        write(*,*)'grid_setbndslarge: Done.'
      end if

      return
      end       ! END: grid_setbndslarge

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_zeroerror:                                                           cc
cc                                                                            cc
cc      When using a shadow hierarchy for TRE, we need to zero out the        cc
cc      the error before looping over fields to compute the error.            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_zeroerror(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      integer  nx,ny,nz, i,j,k

      logical     ltrace
      parameter ( ltrace = .false.)

      call load_pointers(gi)

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

      if (ltrace) then
        write(*,*)'grid_zeroerror: gi     = ',gi
        write(*,*)'grid_zeroerror: nx/y/z = ',nx,ny,nz
      end if

      call load_scal1D(q(gr_error),0.d0,nx*ny*nz)

      if (ltrace) then
        write(*,*)'grid_zeroerror: Done.'
      end if

      return
      end       ! END: grid_zeroerror

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_extfind:                                                             cc
cc                                                                            cc
cc      Find the extremum and its coordinates on a given grid but do so       cc
cc      by considering any finer subregions occuring on the grid.             cc
cc      NB: the maximam comes into this routine, so do not re-initialize.     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_extfind(gi,field, extreme, extcoords, OPTION)
      implicit none
      integer  gi, field, OPTION
      real(kind=8) extreme, extcoords(3)
      include 'grid.inc'
      integer      MIN,     MAX,     ABSMIN,     ABSMAX
      parameter (  MIN = 0, MAX = 1, ABSMIN = 2, ABSMAX = 3 )

      integer  nx,ny,nz, i,j,k, index, gi_tmp, lev
      real(kind=8) abstmp, tmp, h
      integer      UNMASKED,      MASKED 
      parameter (  UNMASKED  = 1, MASKED = 0 )
      integer  index_i_j_k,     index_i_jp1_k,
     *         index_i_j_kp1,   index_i_jp1_kp1
      integer  index_ip1_j_k,   index_ip1_jp1_k,
     *         index_ip1_j_kp1, index_ip1_jp1_kp1
      integer  mini1, maxi1, minj1, maxj1, mink1, maxk1
      integer  mini2, maxi2, minj2, maxj2, mink2, maxk2
      integer  length, myid
      logical  mask
      !
      integer  proc_return_myid
      external proc_return_myid
      integer  level_return_start, grid_return_sibling
      external level_return_start, grid_return_sibling
      logical  grid_intersect,     grid_return_existence 
      external grid_intersect,     grid_return_existence

      logical     ltrace
      parameter ( ltrace = .false.)


      if (ltrace) then
         write(*,*) 'grid_extfind: gi/field/OPTION: ',gi,field,OPTION
      end if

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      h    = gr_h(gi)
      lev  = gr_level(gi)
      myid = proc_return_myid()

      if (ltrace)
     *call field_out3d(q(gr_tmp),gr_t(gi), 'extmaskB',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)

      !
      ! Compute the extremum now subject to the mask:
      !
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         index       = (k-1)*ny*nx + (j-1)*nx + (i-1)
         !
         index_i_j_k       = gfunc_pointer(field,gi) +
     *              (k-1)*ny*nx + (j-1)*nx + (i-1)
         !
         if ( NINT(q(gr_tmp+index)) .eq. UNMASKED ) then
            tmp          = q(index_i_j_k)
            abstmp       = abs(tmp)
            if ( (OPTION.EQ.MAX    .and.    tmp .gt. extreme) .or.
     *           (OPTION.eq.ABSMAX .and. abstmp .gt. extreme) .or.
     *           (OPTION.eq.ABSMIN .and. abstmp .lt. extreme) .or.
     *           (OPTION.eq.MIN    .and.    tmp .lt. extreme) ) then
               extreme      = tmp
               if (OPTION.eq.ABSMAX.or.OPTION.eq.ABSMIN)
     .            extreme   = abstmp
               extcoords(1) = q(gr_x(gi) + i-1)
               extcoords(2) = q(gr_y(gi) + j-1)
               extcoords(3) = q(gr_z(gi) + k-1)
            end if
         end if
      end do
      end do
      end do

      if (ltrace) then
         write(*,*) 'grid_extfind: Done. '
      end if

      return
      end       ! END: grid_extfind


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_set_anamask:                                                         cc
cc                                                                            cc
cc      Set the analysis mask for use by the analysis functions (e.g. integra)cc
cc          (1) Mask out points not to be counted:                            cc
cc              a) Mask out where other grids on level are                    cc
cc                 To deal with overlap, to ensure points are only            cc
cc                 counted once, we need to.....                              cc
cc              b) Mask out where finer level exists                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_set_anamask(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      integer  nx,ny,nz, i,j,k, index, gi_tmp, lev
      integer      UNMASKED,      MASKED 
      parameter (  UNMASKED  = 1, MASKED = 0 )
      real(kind=8) tmp, h
      integer  index_i_j_k,     index_i_jp1_k,
     *         index_i_j_kp1,   index_i_jp1_kp1
      integer  index_ip1_j_k,   index_ip1_jp1_k,
     *         index_ip1_j_kp1, index_ip1_jp1_kp1
      integer  mini1, maxi1, minj1, maxj1, mink1, maxk1
      integer  mini2, maxi2, minj2, maxj2, mink2, maxk2
      integer  length, myid
      logical  mask
      !
      integer  proc_return_myid
      external proc_return_myid
      integer  level_return_start, grid_return_sibling
      external level_return_start, grid_return_sibling
      logical  grid_intersect,     grid_return_existence 
      external grid_intersect,     grid_return_existence

      logical     ltrace
      parameter ( ltrace = .false.)


      if (ltrace) then
         write(*,*) 'grid_set_anamask: gi: ',gi
      end if

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      h    = gr_h(gi)
      lev  = gr_level(gi)
      myid = proc_return_myid()

      !
      ! Use gr_tmp as mask storage:
      !      ==1   unmasked
      !      ==0     masked
      !

      ! Unmask all to start:
      call load_scal1D(q(gr_tmp),1.d0*UNMASKED,nx*ny*nz)

      ! Loop over grids on this level
      !   For any pair of grids, need to decide which gets masked
      !   so that points get counted only once.
      !   A criterion which would seem to work is to look 
      !   in the x-direction. If one has the intersection starting
      !   away from 1, that one gets masked. If both have it at 1,
      !   look at next coordinate and so on:
      gi_tmp = level_return_start(lev)
  15  if (grid_return_existence(gi_tmp)) then
         if (ltrace) write(*,*) 'grid_set_anamask: gi_tmp = ',gi_tmp
         if (gi_tmp.ne.gi .and. grid_intersect(gi_tmp,gi)) then
            if(ltrace)write(*,*)'grid_set_anamask:Intersects w/',gi_tmp
            call grid_find_intersection(gi, gi_tmp,length,
     *            mini1, maxi1, minj1, maxj1, mink1, maxk1)
            call grid_find_intersection(gi_tmp, gi,length,
     *            mini2, maxi2, minj2, maxj2, mink2, maxk2)
            mask = .false.
            if (     mini1 .ne. 1 .and. mini2 .eq. 1) then
               mask = .true.
            else if (mini1 .eq. 1 .and. mini2 .ne. 1) then
               mask = .false.
            else if (minj1 .ne. 1 .and. minj2 .eq. 1) then
               mask = .true.
            else if (minj1 .eq. 1 .and. minj2 .ne. 1) then
               mask = .false.
            else if (mink1 .ne. 1 .and. mink2 .eq. 1) then
               mask = .true.
            else if (mink1 .eq. 1 .and. mink2 .ne. 1) then
               mask = .false.
            else if (gi .lt. gi_tmp) then
               mask = .true.
            else
               mask = .false.
            end if
            if (mask) then
               if(ltrace)write(*,*)'grid_set_anamask: Masking w/',gi_tmp
               do k = mink1, maxk1
               do j = minj1, maxj1
               do i = mini1, maxi1
                  index = (k-1)*ny*nx + (j-1)*nx + (i-1)
                  q(gr_tmp+index) = 1.d0*MASKED
               end do
               end do
               end do
            end if
         end if
         gi_tmp = grid_return_sibling(gi_tmp)
         goto 15
      end if

      if (ltrace)
     *call field_out3d(q(gr_tmp),gr_t(gi), 'anamaskA',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)

      ! Loop over grids on the next finer level
      gi_tmp = level_return_start(lev+1)
  25  if (grid_return_existence(gi_tmp)) then
         if(ltrace)write(*,*)'grid_set_anamask:Fine lev gi_tmp=',gi_tmp
         if (grid_intersect(gi_tmp,gi)) then
            if(ltrace)write(*,*)'grid_set_anamask:Intersects w/',gi_tmp
            call grid_find_intersection(gi, gi_tmp,length,
     *            mini1, maxi1, minj1, maxj1, mink1, maxk1)
               do k = mink1, maxk1-1
               do j = minj1, maxj1-1
               do i = mini1, maxi1-1
                  index = (k-1)*ny*nx + (j-1)*nx + (i-1)
                  q(gr_tmp+index) = MASKED
               end do
               end do
               end do
         end if
         gi_tmp = grid_return_sibling(gi_tmp)
         goto 25
      end if

      if (ltrace)
     *call field_out3d(q(gr_tmp),gr_t(gi), 'anamaskB',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)

      !

      if (ltrace) then
         write(*,*) 'grid_set_anamask: Done.'
      end if

      return
      end       ! END: grid_set_anamask


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_intfind:                                                             cc
cc                                                                            cc
cc      Find the integral/l2norm on a given grid subject to being masked      cc
cc      out for intersections with other grids:                               cc
cc          (1) Mask out points not to be counted:                            cc
cc              a) Mask out where other grids on level are                    cc
cc                 To deal with overlap, to ensure points are only            cc
cc                 counted once, we need to.....                              cc
cc              b) Mask out where finer level exists                          cc
cc          (2) Integrate where not masked out                                cc
cc      NB: return number of points over which calculation taken              cc
cc          (used for the l2norm normalization in the end)                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function grid_intfind(gi,field,OPTION,numpts)
      implicit none
      integer  gi, field, OPTION, numpts
      include 'grid.inc'
      integer  nx,ny,nz, i,j,k, index, gi_tmp, lev
      real(kind=8) tmp, h
      integer      INTEGRATE,     L2NORM
      parameter (  INTEGRATE = 2, L2NORM = 3 )
      integer      UNMASKED,      MASKED 
      parameter (  UNMASKED  = 1, MASKED = 0 )
      integer  index_i_j_k,     index_i_jp1_k,
     *         index_i_j_kp1,   index_i_jp1_kp1
      integer  index_ip1_j_k,   index_ip1_jp1_k,
     *         index_ip1_j_kp1, index_ip1_jp1_kp1
      integer  mini1, maxi1, minj1, maxj1, mink1, maxk1
      integer  mini2, maxi2, minj2, maxj2, mink2, maxk2
      integer  length, myid
      logical  mask
      !
      integer  proc_return_myid
      external proc_return_myid
      integer  level_return_start, grid_return_sibling
      external level_return_start, grid_return_sibling
      logical  grid_intersect,     grid_return_existence 
      external grid_intersect,     grid_return_existence

      logical     ltrace
      parameter ( ltrace = .false.)


      if (ltrace) then
         write(*,*) 'grid_intfind: gi/field/OPTION: ',gi,field,OPTION
      end if

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      h    = gr_h(gi)
      lev  = gr_level(gi)
      myid = proc_return_myid()

      !
      ! Use gr_tmp as mask storage:
      !      ==1   unmasked
      !      ==0     masked
      !

      if (ltrace)
     *call field_out3d(q(gr_tmp),gr_t(gi), 'intmaskB',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)

      !
      ! Compute the integral now subject to the mask:
      !
      numpts       = 0
      grid_intfind = 0.d0
      if (OPTION .eq. INTEGRATE) then
         do k = 1, nz-1
         do j = 1, ny-1
         do i = 1, nx-1
            index       = (k-1)*ny*nx + (j-1)*nx + (i-1)
            !
            index_i_j_k       = gfunc_pointer(field,gi) +
     *                 (k-1)*ny*nx + (j-1)*nx + (i-1)
            index_i_jp1_k     = gfunc_pointer(field,gi) +
     *                 (k-1)*ny*nx + (j  )*nx + (i-1)
            index_i_j_kp1     = gfunc_pointer(field,gi) +
     *                 (k  )*ny*nx + (j-1)*nx + (i-1)
            index_i_jp1_kp1   = gfunc_pointer(field,gi) +
     *                 (k  )*ny*nx + (j  )*nx + (i-1)
            !
            index_ip1_j_k     = gfunc_pointer(field,gi) +
     *                 (k-1)*ny*nx + (j-1)*nx + (i  )
            index_ip1_jp1_k   = gfunc_pointer(field,gi) +
     *                 (k-1)*ny*nx + (j  )*nx + (i  )
            index_ip1_j_kp1   = gfunc_pointer(field,gi) +
     *                 (k  )*ny*nx + (j-1)*nx + (i  )
            index_ip1_jp1_kp1 = gfunc_pointer(field,gi) +
     *                 (k  )*ny*nx + (j  )*nx + (i  )
            !
            if ( NINT(q(gr_tmp+index)) .eq. UNMASKED ) then
               tmp          = h**3 * 0.125d0 * (
     *                        q(index_i_j_k      )
     *                      + q(index_i_jp1_k    )
     *                      + q(index_i_j_kp1    )
     *                      + q(index_i_jp1_kp1  )
     *                      + q(index_ip1_j_k    )
     *                      + q(index_ip1_jp1_k  )
     *                      + q(index_ip1_j_kp1  )
     *                      + q(index_ip1_jp1_kp1)
     *                                       )
               grid_intfind = grid_intfind + tmp
               numpts       = numpts + 1
            end if
         end do
         end do
         end do
      else if (OPTION .eq. L2NORM) then
         do k = 1, nz
         do j = 1, ny
         do i = 1, nx
            index       = (k-1)*ny*nx + (j-1)*nx + (i-1)
            index_i_j_k = gfunc_pointer(field,gi) +
     *                 (k-1)*ny*nx + (j-1)*nx + (i-1)
            if ( NINT(q(gr_tmp+index)) .eq. 1 ) then
               tmp          = q(index_i_j_k )**2
               grid_intfind = grid_intfind + tmp
               numpts       = numpts + 1
            end if
         end do
         end do
         end do
         grid_intfind = grid_intfind + tmp
         !write(*,*) 'grid_intfind: L2norm not implemented correctly yet'
      else
         write(*,*) 'grid_intfind: Uknown option: ',OPTION
      end if

      if (ltrace) then
         write(*,*) 'grid_intfind: numpts',numpts
         write(*,*) 'grid_intfind: Done: ',grid_intfind
      end if

      return
      end       ! END: grid_intfind

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_momentfind:                                                          cc
cc                                                                            cc
cc      Compute the moments of a given field D, in particular,                cc
cc      the first:                                                            cc
cc                I^j  = \int D x^j d^3x                                      cc
cc      and second moments:                                                   cc
cc                I^jk = \int D x^j x^k d^3x                                  cc
cc      Moments assumed to be packed into a single array:                     cc
cc                     (  I^xx    I^xy   I^xz  )                              cc
cc                     (  I^x     I^yy   I^yz  )                              cc
cc                     (  I^y     I^z    I^zz  )                              cc
cc      Integrals added to existing values since integrals accumulate         cc
cc      over the various grids.                                               cc
cc      NB: adapted from grid_intfind()                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_momentfind(moms, gi,field)
      implicit none
      real(kind=8) moms(3,3)
      integer  gi, field
      include 'grid.inc'
      integer  nx,ny,nz, i,j,k, index, gi_tmp, lev
      integer  ii,jj
      real(kind=8) tmp, h, dV
      real(kind=8) x_i, x_ip1, y_j,y_jp1, z_k,z_kp1
      integer      UNMASKED,      MASKED 
      parameter (  UNMASKED  = 1, MASKED = 0 )
      integer  index_i_j_k,     index_i_jp1_k,
     *         index_i_j_kp1,   index_i_jp1_kp1
      integer  index_ip1_j_k,   index_ip1_jp1_k,
     *         index_ip1_j_kp1, index_ip1_jp1_kp1
      integer  mini1, maxi1, minj1, maxj1, mink1, maxk1
      integer  mini2, maxi2, minj2, maxj2, mink2, maxk2
      integer  length, myid
      logical  mask
      !
      integer  proc_return_myid
      external proc_return_myid
      integer  level_return_start, grid_return_sibling
      external level_return_start, grid_return_sibling
      logical  grid_intersect,     grid_return_existence 
      external grid_intersect,     grid_return_existence

      logical     ltrace
      parameter ( ltrace = .false.)

      if (ltrace) then
         write(*,*) 'grid_momentfind: gi/field: ',gi,field
      end if

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      h    = gr_h(gi)
      lev  = gr_level(gi)
      myid = proc_return_myid()

      !
      ! Use gr_tmp as mask storage:
      !      ==1   unmasked
      !      ==0     masked
      !

      if (ltrace)
     *call field_out3d(q(gr_tmp),gr_t(gi), 'intmaskB',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)

      !
      ! Compute the integral now subject to the mask:
      !
         do k = 1, nz-1
         do j = 1, ny-1
         do i = 1, nx-1
            index       = (k-1)*ny*nx + (j-1)*nx + (i-1)
            if ( NINT(q(gr_tmp+index)) .eq. UNMASKED ) then
               !
               index_i_j_k       = gfunc_pointer(field,gi) +
     *                 (k-1)*ny*nx + (j-1)*nx + (i-1)
               index_i_jp1_k     = gfunc_pointer(field,gi) +
     *                 (k-1)*ny*nx + (j  )*nx + (i-1)
               index_i_j_kp1     = gfunc_pointer(field,gi) +
     *                 (k  )*ny*nx + (j-1)*nx + (i-1)
               index_i_jp1_kp1   = gfunc_pointer(field,gi) +
     *                 (k  )*ny*nx + (j  )*nx + (i-1)
               !
               index_ip1_j_k     = gfunc_pointer(field,gi) +
     *                 (k-1)*ny*nx + (j-1)*nx + (i  )
               index_ip1_jp1_k   = gfunc_pointer(field,gi) +
     *                 (k-1)*ny*nx + (j  )*nx + (i  )
               index_ip1_j_kp1   = gfunc_pointer(field,gi) +
     *                 (k  )*ny*nx + (j-1)*nx + (i  )
               index_ip1_jp1_kp1 = gfunc_pointer(field,gi) +
     *                 (k  )*ny*nx + (j  )*nx + (i  )
               !
               x_i   = q(gr_x(gi)+i  -1)
               x_ip1 = q(gr_x(gi)+i+1-1)
               y_j   = q(gr_y(gi)+j  -1)
               y_jp1 = q(gr_y(gi)+j+1-1)
               z_k   = q(gr_z(gi)+k  -1)
               z_kp1 = q(gr_z(gi)+k+1-1)
               do jj=1,3
               do ii=1,3
                  dV        = h**3 * 0.125d0 
                  if (ii.eq.1.and.jj.eq.1) then
                     tmp    = dV *(
     *                        q(index_i_j_k      )*x_i**2
     *                      + q(index_i_jp1_k    )*x_i**2
     *                      + q(index_i_j_kp1    )*x_i**2
     *                      + q(index_i_jp1_kp1  )*x_i**2
     *                      + q(index_ip1_j_k    )*x_ip1**2
     *                      + q(index_ip1_jp1_k  )*x_ip1**2
     *                      + q(index_ip1_j_kp1  )*x_ip1**2
     *                      + q(index_ip1_jp1_kp1)*x_ip1**2 )
                  else if (ii.eq.2.and.jj.eq.2) then
                     tmp    = dV *(
     *                        q(index_i_j_k      )*y_j**2
     *                      + q(index_i_jp1_k    )*y_j**2
     *                      + q(index_i_j_kp1    )*y_j**2
     *                      + q(index_i_jp1_kp1  )*y_j**2
     *                      + q(index_ip1_j_k    )*y_jp1**2
     *                      + q(index_ip1_jp1_k  )*y_jp1**2
     *                      + q(index_ip1_j_kp1  )*y_jp1**2
     *                      + q(index_ip1_jp1_kp1)*y_jp1**2 )
                  else if (ii.eq.3.and.jj.eq.3) then
                     tmp    = dV *(
     *                        q(index_i_j_k      )*z_k**2
     *                      + q(index_i_jp1_k    )*z_k**2
     *                      + q(index_i_j_kp1    )*z_k**2
     *                      + q(index_i_jp1_kp1  )*z_k**2
     *                      + q(index_ip1_j_k    )*z_kp1**2
     *                      + q(index_ip1_jp1_k  )*z_kp1**2
     *                      + q(index_ip1_j_kp1  )*z_kp1**2
     *                      + q(index_ip1_jp1_kp1)*z_kp1**2 )
                  else if (ii.eq.1.and.jj.eq.2) then
                     tmp    = dV *(
     *                        q(index_i_j_k      )*x_i*y_j
     *                      + q(index_i_jp1_k    )*x_i*y_j
     *                      + q(index_i_j_kp1    )*x_i*y_j
     *                      + q(index_i_jp1_kp1  )*x_i*y_j
     *                      + q(index_ip1_j_k    )*x_ip1*y_jp1
     *                      + q(index_ip1_jp1_k  )*x_ip1*y_jp1
     *                      + q(index_ip1_j_kp1  )*x_ip1*y_jp1
     *                      + q(index_ip1_jp1_kp1)*x_ip1*y_jp1 )
                  else if (ii.eq.1.and.jj.eq.3) then
                     tmp    = dV *(
     *                        q(index_i_j_k      )*x_i*z_k
     *                      + q(index_i_jp1_k    )*x_i*z_k
     *                      + q(index_i_j_kp1    )*x_i*z_k
     *                      + q(index_i_jp1_kp1  )*x_i*z_k
     *                      + q(index_ip1_j_k    )*x_ip1*z_kp1
     *                      + q(index_ip1_jp1_k  )*x_ip1*z_kp1
     *                      + q(index_ip1_j_kp1  )*x_ip1*z_kp1
     *                      + q(index_ip1_jp1_kp1)*x_ip1*z_kp1 )
                  else if (ii.eq.2.and.jj.eq.3) then
                     tmp    = dV *(
     *                        q(index_i_j_k      )*z_k*y_j
     *                      + q(index_i_jp1_k    )*z_k*y_j
     *                      + q(index_i_j_kp1    )*z_k*y_j
     *                      + q(index_i_jp1_kp1  )*z_k*y_j
     *                      + q(index_ip1_j_k    )*z_kp1*y_jp1
     *                      + q(index_ip1_jp1_k  )*z_kp1*y_jp1
     *                      + q(index_ip1_j_kp1  )*z_kp1*y_jp1
     *                      + q(index_ip1_jp1_kp1)*z_kp1*y_jp1 )
                  else if (ii.eq.2.and.jj.eq.1) then
                     ! Now, compute the first moments:
                     tmp    = dV *(
     *                        q(index_i_j_k      )*x_i
     *                      + q(index_i_jp1_k    )*x_i
     *                      + q(index_i_j_kp1    )*x_i
     *                      + q(index_i_jp1_kp1  )*x_i
     *                      + q(index_ip1_j_k    )*x_ip1
     *                      + q(index_ip1_jp1_k  )*x_ip1
     *                      + q(index_ip1_j_kp1  )*x_ip1
     *                      + q(index_ip1_jp1_kp1)*x_ip1 )
                  else if (ii.eq.3.and.jj.eq.1) then
                     tmp    = dV *(
     *                        q(index_i_j_k      )*y_j
     *                      + q(index_i_jp1_k    )*y_j
     *                      + q(index_i_j_kp1    )*y_j
     *                      + q(index_i_jp1_kp1  )*y_j
     *                      + q(index_ip1_j_k    )*y_jp1
     *                      + q(index_ip1_jp1_k  )*y_jp1
     *                      + q(index_ip1_j_kp1  )*y_jp1
     *                      + q(index_ip1_jp1_kp1)*y_jp1 )
                  else if (ii.eq.3.and.jj.eq.2) then
                     tmp    = dV *(
     *                        q(index_i_j_k      )*z_k
     *                      + q(index_i_jp1_k    )*z_k
     *                      + q(index_i_j_kp1    )*z_k
     *                      + q(index_i_jp1_kp1  )*z_k
     *                      + q(index_ip1_j_k    )*z_kp1
     *                      + q(index_ip1_jp1_k  )*z_kp1
     *                      + q(index_ip1_j_kp1  )*z_kp1
     *                      + q(index_ip1_jp1_kp1)*z_kp1 )
                  end if
                  moms(ii,jj) = moms(ii,jj) + tmp
                  !grid_intfind = grid_intfind + tmp
               end do
               end do
            end if
         end do
         end do
         end do

      if (ltrace) then
         write(*,*) 'grid_momentfind: Done.'
         !write(*,*) 'grid_momentfind: numpts',numpts
      end if

      return
      end       ! END: grid_momentfind

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_get_fname:                                                           cc
cc                  Provide the field name associate with a given field numbercc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_get_fname(fieldnum,fieldname)
      implicit       none
      integer        fieldnum
      character(*) fieldname
      include 'grid.inc'
      !include 'grid_methods.inc'
      !include 'glob.inc'

      fieldname = gfunc_name(fieldnum)

      return
      end       ! END: grid_get_fname

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_syncbndi_slave:                                                      cc
cc                 Sync the boundaries only of ELLIPTIC field.                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_syncbndi_slave(gi1, gi2)
      implicit    none
      integer     gi1, gi2
      include     'grid.inc'
      include     'grid_methods.inc'
      include     'glob.inc'
      include     'chr.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     parent, nxs,nys,nzs, length
      integer     lx, rx, ly, ry, lz, rz
      integer     nx1, ny1, nz1, nx2, ny2, nz2, i,j,k,
     *            min_i1, max_i1,
     *            min_j1, max_j1,
     *            min_k1, max_k1,
     *            min_i2, max_i2,
     *            min_j2, max_j2,
     *            min_k2, max_k2,
     *            properlength, maxlength,
     *            length1, length2, owner, index
      integer     request_recv, request_send, gr_chr_gi1
      logical     first
      logical     fvflag
      logical     double_equal
      external    double_equal
      integer     calc_bnd_lenb
      external    calc_bnd_lenb
      logical     blockingsend
      parameter ( blockingsend  = .false. )
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )


 98   format('[',I3,'] ',A,3F15.9)
 99   format('[',I3,'] ',A,3I7)

      if (ltrace) then
         write(*,99) myid,'grid_syncbndi_slave: gi, gj: ', gi1,gi2
      end if

      !
      ! Get pointers to temporary array:
      !
      call load_pointers(gi1)
      gr_chr_gi1 = gr_chr
      owner = grid_return_owner(gi2)

      nx1 = gr_nx(gi1)
      ny1 = gr_ny(gi1)
      nz1 = gr_nz(gi1)

      nx2 = gr_nx(gi2)
      ny2 = gr_ny(gi2)
      nz2 = gr_nz(gi2)

      call grid_find_intersection(gi1, gi2,length1,
     *            min_i1, max_i1, min_j1, max_j1, min_k1, max_k1)

      !
      ! If no intersection:
      !
      if (length1 .eq. 0) return
      !
      ! If grids do not overlap at their ends
      !         more than twice the bound_width:
      !
      if (    ((max_i1-min_i1).lt.bound_width-1 .and.
     *                (max_i1.eq.nx1.or.min_i1.eq.1))
     *    .or.((max_j1-min_j1).lt.bound_width-1 .and.
     *                (max_j1.eq.ny1.or.min_j1.eq.1))
     *    .or.((max_k1-min_k1).lt.bound_width-1 .and.
     *                (max_k1.eq.nz1.or.min_k1.eq.1)) ) then
         if(ltrace)write(*,99)myid,'grid_syncbndi_slv:not syncg',gi1,gi2
         return
      end if


      call grid_find_intersection(gi2, gi1,length2,
     *            min_i2, max_i2, min_j2, max_j2, min_k2, max_k2)


      !
      !
      ! Not needed anymore:
      !
      lx = 1
      rx = 1
      ly = 1
      ry = 1
      lz = 1
      rz = 1

      !
      ! Above lengths don't work for just the boundary:
      !
      length1 = calc_bnd_lenb(min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx1,ny1,nz1, nx2,ny2,nz2,bound_width)
      length2 = calc_bnd_lenb(min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx2,ny2,nz2, nx1,ny1,nz1,bound_width)


      maxlength = length1
      if (length2 .gt. maxlength) maxlength = length2
      !
      ! If no "effective" intersection:
      !
      if (length1 .eq. 0) return
      !
      ! If overlap is just one point wide, do not sync
      ! since those points can't be interior
      !
      if (  (min_i1.eq.max_i1 .and .max_i1.eq.nx1) .or.
     *      (min_j1.eq.max_j1 .and .max_j1.eq.ny1) .or.
     *      (min_k1.eq.max_k1 .and .max_k1.eq.nz1) .or.
     *      (min_i2.eq.max_i2 .and .max_i2.eq.nx2) .or.
     *      (min_j2.eq.max_j2 .and .max_j2.eq.ny2) .or.
     *      (min_k2.eq.max_k2 .and .max_k2.eq.nz2) ) then
         if(ltrace)write(*,99)myid,'grid_syncbndi_slv:Notsyncng',gi1,gi2
         return
      end if

      if (ltrace) then
      write(*,99)myid,'grid_syncbndi_slave: gi1       = ', gi1
      write(*,99)myid,'grid_syncbndi_slave: gi2       = ', gi2
      call field_dump_info(q(gr_chr),nx1,ny1,nz1)
      write(*,99)myid,'grid_syncbndi_slave: length1&2 =',length1,length2
      write(*,99)myid,'grid_syncbndi_slave: maxlength = ',maxlength
      write(*,99)myid,'grid_syncbndi_slave: min_i,j,k1= ',min_i1,min_j1,
     * min_k1
      write(*,99)myid,'grid_syncbndi_slave: max_i,j,k1= ',max_i1,max_j1,
     * max_k1
      write(*,99)myid,'grid_syncbndi_slave: nx,ny,nz1 = ',nx1,ny1,nz1
      write(*,99)myid,'grid_syncbndi_slave: min_i,j,k2= ',min_i2,min_j2,
     * min_k2
      write(*,99)myid,'grid_syncbndi_slave: max_i,j,k2= ',max_i2,max_j2,
     * max_k2
      write(*,99)myid,'grid_syncbndi_slave: nx,ny,nz2 = ',nx2,ny2,nz2
      write(*,99)myid,'grid_syncbndi_slave: lx, rx    = ',lx,rx
      write(*,99)myid,'grid_syncbndi_slave: ly, ry    = ',ly,ry
      write(*,99)myid,'grid_syncbndi_slave: lz, rz    = ',lz,rz
      write(*,98)myid,'grid_syncbndi_slave: gr_t(gi1) = ',gr_t(gi1)
      end if

      ! All elliptic-type fields currently are finite difference:
      fvflag = .false.

      first = .true.
      do i = 1, num_gfuncs
         !
         ! Keep in mind that this boundary syncing 
         ! takes place *during* the iterative crank-nicholson
         ! and hence we apply the sync'ing to the advanced time level
         !
         if ( gfunc_type(i) .eq.GFUNC_INTEGRAL_DEREL ) then
            !
            if (ltrace2) write(*,99)myid,'  doing field i= ',i
            !
            ! Non-blocking recv first:
            !
            if (ltrace2) write(*,99)myid,'  doing IRecv    i= ',i
            call MPI_IRecv(q(gr_flag), length2, MPI_DOUBLE_PRECISION,
     *           owner, TAG_SYNC_DATA, MPI_COMM_WORLD,request_recv,ierr)
            if (ltrace2) write(*,99)myid,'  IRecv   ierr=',ierr
            if (ltrace2) write(*,99)myid,'  request_recv=',request_recv
            !
            ! Make sure tmp space has been sent, so we can write over it:
            !
            if (.not.first) then
               if (.not.blockingsend) then
                  if(ltrace2)
     *               write(*,99)myid,' Waiting.. send',request_send
                  call MPI_Wait(request_send, status, ierr)
                  if (ltrace2) write(*,99)myid,'  Wait on send done.'
               end if
               if(ltrace2)write(*,99)myid,'  Wait ierr: ',ierr
            else
               first = .false.
            end if
            !
            ! Copy data into tmp now that it's safe:
            !
            if (ltrace2) write(*,99)myid,'  Copying data from ',gi1
            call copy_int_bndd(   q(gr_tmp),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           length1, properlength,
     *            nx1,ny1,nz1, nx2,ny2,nz2, bound_width,fvflag)
            !
            ! Non-blocking send:
            !
            if (.not.blockingsend) then
               if(ltrace2)write(*,99)myid,'  Doing ISend'
               call MPI_ISend(q(gr_tmp), length1, MPI_DOUBLE_PRECISION,
     *           owner, TAG_SYNC_DATA,MPI_COMM_WORLD, request_send,ierr)
               if (ltrace2) then
                   write(*,99)myid,'  ISend ierr: ',ierr
                   write(*,99)myid,'  ISend requestsd:',request_send
               end if
            else
               if (ltrace2) write(*,99)myid,'  Blocking send'
               call MPI_Send(q(gr_tmp), length1, MPI_DOUBLE_PRECISION,
     *           owner, TAG_SYNC_DATA,MPI_COMM_WORLD, ierr)
               if (ltrace2) write(*,99)myid,'  MPI_Send ierr=',ierr
            end if
            !
            if (ltrace2) then
               write(*,99)myid,' properlength =',properlength
               write(*,99)myid,' Sending data to  proc: ',owner
               write(*,99)myid,'             of length: ',length1
               write(*,99)myid,'                w/ tag: ',TAG_SYNC_DATA
               write(*,99)myid,'                w/ com: ',MPI_COMM_WORLD
            end if

            !
            ! Wait until data has been received:
            !
            if(ltrace2)write(*,99)myid,' Waiting.. recv',request_recv
            call MPI_Wait(request_recv, status, ierr)
            if(ltrace2)write(*,99)myid,'  Wait ierr: ',ierr
            if (ltrace2) write(*,99)myid,'  Wait on recv done.'
            !
            if (ltrace2) write(*,99) myid,'  Pasting data from ',owner
            !call load_scal3d(q(gr_flag), -gi2*1.d0, nx1,ny1,nz1 )
            call paste_int_bndd(  q(gr_flag),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           nx2,ny2,nz2,
     *                           length2,
     *                           nx1,ny1,nz1,bound_width,fvflag )
         end if
      end do

      if (.not.first .and. .not.blockingsend) then
         if(ltrace2)write(*,99)myid,' Waiting.. send',request_send
         call MPI_Wait(request_send, status, ierr)
         if(ltrace2)write(*,99)myid,'  Wait ierr: ',ierr
      end if

      if (ltrace.or.ltrace2) then
       write(*,99) myid,'grid_syncbndi_slave: Done. ',gi1,gi2
      end if

      return
      end    ! END: grid_syncbndi_slave

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_syncbndi:                                                            cc
cc                 Sync the boundaries for just *ELLIPTIC* fields.            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_syncbndi(gi1, gi2)
      implicit    none
      integer     gi1, gi2
      include     'grid.inc'
      include     'grid_methods.inc'
      include     'glob.inc'
      include     'chr.inc'
      integer     parent, nxs,nys,nzs, length1, length2
      integer     lx, rx, ly, ry, lz, rz
      integer     nx1, ny1, nz1, nx2, ny2, nz2, i,j,k,
     *            min_i1, max_i1,
     *            min_j1, max_j1,
     *            min_k1, max_k1,
     *            min_i2, max_i2,
     *            min_j2, max_j2,
     *            min_k2, max_k2,   index,
     *            tmpptr1, tmpptr2, properlength
      integer     gr_chr_gi1, gr_chr_gi2
      logical     fvflag
      logical     double_equal
      external    double_equal
      integer     calc_bnd_lenb
      external    calc_bnd_lenb
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      !
      ! Get pointers to temporary array:
      !
      call load_pointers(gi2)
      gr_chr_gi2 = gr_chr
      call load_pointers(gi1)
      gr_chr_gi1 = gr_chr


      nx1 = gr_nx(gi1)
      ny1 = gr_ny(gi1)
      nz1 = gr_nz(gi1)

      nx2 = gr_nx(gi2)
      ny2 = gr_ny(gi2)
      nz2 = gr_nz(gi2)

      call grid_find_intersection(gi1, gi2, length1,
     *            min_i1, max_i1, min_j1, max_j1, min_k1, max_k1)

      !
      ! If no intersection:
      !
      if (length1 .eq. 0) return
      !
      ! If grids do not overlap at their ends
      !         more than twice the bound_width:
      !
      if (    ((max_i1-min_i1).lt.bound_width-1 .and.
     *                (max_i1.eq.nx1.or.min_i1.eq.1))
     *    .or.((max_j1-min_j1).lt.bound_width-1 .and.
     *                (max_j1.eq.ny1.or.min_j1.eq.1))
     *    .or.((max_k1-min_k1).lt.bound_width-1 .and.
     *                (max_k1.eq.nz1.or.min_k1.eq.1)) ) then
         if(ltrace)write(*,*)     'grid_syncbndi:not syncg',gi1,gi2
         return
      end if

      call grid_find_intersection(gi2, gi1, length2,
     *            min_i2, max_i2, min_j2, max_j2, min_k2, max_k2)

      !
      !
      ! Not needed anymore:
      !
      lx = 1
      rx = 1
      ly = 1
      ry = 1
      lz = 1
      rz = 1

      !
      ! The above returned values of the length
      ! are for a vector holding the *entire* 
      ! intersection. We only need to store part of the boundary:
      !
      length1 = calc_bnd_lenb(min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx1,ny1,nz1, nx2,ny2,nz2,bound_width)
      length2 = calc_bnd_lenb(min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx2,ny2,nz2, nx1,ny1,nz1,bound_width)


      if (length1.gt. nx1*ny1*nz1) then
         write(*,*) 'grid_syncbndi: Not enough space for g1'
         call my_exit('grid_syncbndi: Not enough space for g1')
      end if
      if (length2.gt. nx2*ny2*nz2) then
         write(*,*) 'grid_syncbndi: Not enough space for g2'
         call my_exit('grid_syncbndi: Not enough space for g2')
      end if
      !
      ! If no "effective" intersection:
      !
      if (length1 .eq. 0) return
      !
      ! If overlap is just one point wide, do not sync
      ! since those points can't be interior
      !
      if (  (min_i1.eq.max_i1 .and .max_i1.eq.nx1) .or.
     *      (min_j1.eq.max_j1 .and .max_j1.eq.ny1) .or.
     *      (min_k1.eq.max_k1 .and .max_k1.eq.nz1) .or.
     *      (min_i2.eq.max_i2 .and .max_i2.eq.nx2) .or.
     *      (min_j2.eq.max_j2 .and .max_j2.eq.ny2) .or.
     *      (min_k2.eq.max_k2 .and .max_k2.eq.nz2) ) then
         if(ltrace)write(*,*)'grid_syncbndi: Not syncing ',gi1,gi2
         return
      end if

      if (ltrace2) then
         write(*,*) 'grid_syncbndi: gi1, gi2  = ', gi1, gi2
         write(*,*) 'grid_syncbndi: gi1    = ', gi1
         write(*,*) 'grid_syncbndi: gi2    = ', gi2
         write(*,*) 'grid_syncbndi: min_i1 = ', min_i1
         write(*,*) 'grid_syncbndi: max_i1 = ', max_i1
         write(*,*) 'grid_syncbndi: min_j1 = ', min_j1
         write(*,*) 'grid_syncbndi: max_j1 = ', max_j1
         write(*,*) 'grid_syncbndi: min_k1 = ', min_k1
         write(*,*) 'grid_syncbndi: max_k1 = ', max_k1
         write(*,*) 'grid_syncbndi: min_i2 = ', min_i2
         write(*,*) 'grid_syncbndi: max_i2 = ', max_i2
         write(*,*) 'grid_syncbndi: min_j2 = ', min_j2
         write(*,*) 'grid_syncbndi: max_j2 = ', max_j2
         write(*,*) 'grid_syncbnd: min_k2 = ', min_k2
         write(*,*) 'grid_syncbndi: max_k2 = ', max_k2
         write(*,*) 'grid_syncbndi: nX1    = ', nx1,ny1,nz1
         write(*,*) 'grid_syncbndi: nX2    = ', nx2,ny2,nz2
         write(*,*) 'grid_syncbndi: length1= ', length1
         write(*,*) 'grid_syncbndi: length2= ', length2
      end if
      if (ltrace) then
       call level_tree_dump()
       write(*,*)'grid_syncbndi: gi1       = ', gi1
       call grid_dump_info(gi1)
       write(*,*)'grid_syncbndi: gi2       = ', gi2
       call grid_dump_info(gi2)
       write(*,*)'grid_syncbndi: length1= ', length1
       write(*,*)'grid_syncbndi: length2= ', length2
       write(*,*)'grid_syncbndi: min_i,j,k1= ',min_i1,min_j1,min_k1
       write(*,*)'grid_syncbndi: max_i,j,k1= ',max_i1,max_j1,max_k1
       write(*,*)'grid_syncbndi: nx,ny,nz1 = ',nx1,ny1,nz1
       write(*,*)'grid_syncbndi: min_i,j,k2= ',min_i2,min_j2,min_k2
       write(*,*)'grid_syncbndi: max_i,j,k2= ',max_i2,max_j2,max_k2
       write(*,*)'grid_syncbndi: nx,ny,nz2 = ',nx2,ny2,nz2
      end if

      ! All elliptic-type fields currently are finite difference:
      fvflag = .false.

      do i = 1, num_gfuncs
         !
         ! Keep in mind that this boundary syncing 
         ! takes place *during* the iterative crank-nicholson
         ! and hence we apply the sync'ing to the advanced time level
         !
         if ( gfunc_type(i) .eq.GFUNC_INTEGRAL_DEREL) then
            !
            if (ltrace2)
     *       write(*,*) 'grid_syncbndi: gfunc_name:',gfunc_name(i)(1:10)
            !
            call copy_int_bndd(   q(gr_tmp),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           length1, properlength,
     *            nx1,ny1,nz1, nx2,ny2,nz2, bound_width,fvflag)
            if(ltrace2)write(*,*)'grid_syncbndi:proplength',properlength
            call paste_int_bndd(  q(gr_tmp),
     *                           q(gfunc_pointer(i,gi2)),
     *                           q(gr_chr_gi2),
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                           lx,rx, ly,ry, lz,rz,
     *                           nx1,ny1,nz1, properlength,
     *                           nx2,ny2,nz2,bound_width,fvflag )
            !
            !
            call copy_int_bndd(   q(gr_tmp),
     *                           q(gfunc_pointer(i,gi2)),
     *                           q(gr_chr_gi2),
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                           lx,rx, ly,ry, lz,rz,
     *                           length2, properlength,
     *            nx2,ny2,nz2, nx1,ny1,nz1,bound_width,fvflag)
            if(ltrace2)write(*,*)'grid_syncbndi:proplength',properlength
            call paste_int_bndd(  q(gr_tmp),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           nx2,ny2,nz2, properlength,
     *                           nx1,ny1,nz1,bound_width,fvflag )
            !
         end if
         !
      end do

      return
      end    ! END: grid_syncbndi

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_check_mask:                                                          cc
cc                    Purely for debugging, check if mask is purely           cc
cc                    rectangular in that two points define the corners       cc
cc                    of a 3D region which contains *only* masked points      cc
cc                    and *no* masked points appear outside this region.      cc
cc                                                                            cc
cc                NB: with domain decomposition, it would seem that any       cc
cc                    rectangular masked region would be decomposed into      cc
cc                    purely rectangular sub-regions.                         cc
cc                                                                            cc
cc                NB: Assumes, if mask>0 point is unmasked.                   cc
cc                NB: intlabel is so that one can label the sdfdump           cc
cc                    based on where this check is occuring in the code.      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_check_mask(gi,intlabel)
      implicit none
      integer  gi,intlabel
      include 'grid.inc'
      include 'mask.inc'
      integer  nx,ny,nz
      integer  i,j,k, l, myid
      real(kind=8) x,y,z
      integer  num_masked_outside, num_unmasked_inside
      integer  num_nonunitary, num_zero
      integer  maskpt
      integer  proc_return_myid
      external proc_return_myid
      logical  quickreturn, inmaskedregion, problem
      character(2) tmps

      myid = proc_return_myid()

      call load_pointers(gi)

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

      ! Check if no masks active:
      quickreturn = .true.
      do l = 1, max_num_masks
         quickreturn = quickreturn .and. .not.bh_true(l)
      end do
      if (quickreturn) then
          write(*,99)myid,'grid_check_mask: No active masks to check',
     *             gi,intlabel
         return
      !else
          !write(*,99)myid,'grid_check_mask: Checking grid: ',gi
      end if

      ! Count various bad points:
      num_zero            = 0
      num_nonunitary      = 0
      num_masked_outside  = 0
      num_unmasked_inside = 0
      problem             = .false.
      !
      ! Loop over points:
      !
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         !
         maskpt = NINT(q(gr_mask + (k-1)*ny*nx + (j-1)*nx + (i-1)))
         x      = q(gr_x(gi)+i-1)
         y      = q(gr_y(gi)+j-1)
         z      = q(gr_z(gi)+k-1)
         ! Check if point is +1 or -1:
         if (ABS(maskpt).ne.1)  num_nonunitary = num_nonunitary+1
         ! Check if in any active mask regions:
         inmaskedregion = .false.
         do l = 1, max_num_masks
            if (bh_true(l)) then
               inmaskedregion = inmaskedregion .or. (
     *           ( x  .gt. (mask_coords(1,l)+FUZZ) ) .and.
     *           ( x  .lt. (mask_coords(2,l)-FUZZ) ) .and.
     *           ( y  .gt. (mask_coords(3,l)+FUZZ) ) .and.
     *           ( y  .lt. (mask_coords(4,l)-FUZZ) ) .and.
     *           ( z  .gt. (mask_coords(5,l)+FUZZ) ) .and.
     *           ( z  .lt. (mask_coords(6,l)-FUZZ) ) )
            end if
         end do
         !
         if (maskpt .eq. 0) then
            num_zero = num_zero + 1
            write(*,99)myid,'  zero: i/j/k: ',i,j,k
            write(*,98)myid,'  zero: x/y/z: ',x,y,z
         else if (inmaskedregion .and. maskpt .gt. 0) then
            num_unmasked_inside = num_unmasked_inside + 1
            write(*,99)myid,'  unmasked-inside: i/j/k: ',i,j,k
            write(*,98)myid,'  unmasked-inside: x/y/z: ',x,y,z
         else if (.not. inmaskedregion .and. maskpt .lt. 0) then
            num_masked_outside  = num_masked_outside  + 1
            write(*,99)myid,'  masked-outside: i/j/k: ',i,j,k
            write(*,98)myid,'  masked-outside: x/y/z: ',x,y,z
         end if
         !
      end do
      end do
      end do


      if (num_zero       .gt. 0) then
          write(*,99)myid,'grid_check_mask: Mask contains zeros',
     *                 nx,ny,nz
          write(*,99) myid, 'grid_check_mask: num_zero: ',
     *                       num_zero
          problem = .true.
      end if
!     if (num_nonunitary .gt. 0) then
!         write(*,99)myid,'grid_check_mask: Mask is not unitary!',
!    *                 nx,ny,nz
!         write(*,99) myid, 'grid_check_mask: num_nonunitary: ',
!    *                       num_nonunitary
!     end if
      if (num_unmasked_inside .gt. 0) then
          write(*,99)myid,'grid_check_mask: Mask is not rectangular!',
     *                 nx,ny,nz
          write(*,99)myid,'grid_check_mask: num_unmasked_inside: ',
     *                       num_unmasked_inside
          problem = .true.
      end if
      if (num_masked_outside .gt. 0) then
          write(*,99) myid, 'grid_check_mask: Mask is not rectangular!',
     *                 nx,ny,nz
          write(*,99) myid, 'grid_check_mask: num_masked_outside: ',
     *                       num_masked_outside
          problem = .true.
      end if
      if (num_masked_outside .le. 0 .and. num_unmasked_inside.le.0) then
          write(*,99)myid,'grid_check_mask: Mask is rectangular gi:',gi,
     *            intlabel
      end if

      if (problem) then
         write(*,99) myid, 'grid_check_mask: Dumping info  and mask',gi,
     *         intlabel
         write(*,98) myid, 'grid_check_mask: at t = ',gr_t(gi)
         call field_dump_info(q(gr_mask),nx,ny,nz)
         !call grid_dump_info(gi)
         call int2str(intlabel,tmps)
         call field_out3d(q(gr_mask),gr_t(gi), 'check_mask'//tmps,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
      end if

 98   format('[',I3,'] ',A,3F15.9)
 99   format('[',I3,'] ',A,3I7)
      return
      end       ! END: grid_check_mask

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_symm:                                                                cc
cc                    Populate extended boundary according to assumed         cc
cc                    symmetries.                                             cc
cc                    Grid must be local to this proc.                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_symm(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      include 'glob.inc'
      include 'mpif.h'
      include 'mpi_stuff.inc'
      integer  i
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      if(ltrace)write(*,99)myid,'grid_symm: Enforcing symmetry gi: ',gi

      if (gr_minz(gi) .le. 0 .and.
     *    (assume_symmetry .eq. 3 .or. assume_symmetry.eq.6)) then
         do i = 1, num_gfuncs
            if (gfunc_name(i) .ne. 'chr'    .and.
     *          gfunc_name(i) .ne. 'mask'   .and.
     *          gfunc_name(i) .ne. 'r'      .and.
     *          gfunc_name(i) .ne. 'xphys'  .and.
     *          gfunc_name(i) .ne. 'yphys'  .and.
     *          gfunc_name(i) .ne. 'zphys'  .and.
     *          gfunc_name(i) .ne. 'cctk_x' .and.
     *          gfunc_name(i) .ne. 'cctk_y' .and.
     *          gfunc_name(i) .ne. 'cctk_z'      ) then
               if(ltrace2)
     *         write(*,*)'grid_symm: field/i/zsym:',gfunc_name(i),i,
     *                            gfunc_zsym(i)
               call apply_symmetryz(q(gfunc_pointer(i,gi)), 
     *                           q(gr_z(gi)),
     *                           gfunc_zsym(i), 
     *                           gr_nx(gi),gr_ny(gi),gr_nz(gi))
            end if
         end do
      end if

      if (gr_miny(gi) .le. 0 .and.
     *    (assume_symmetry .eq. 2 .or. assume_symmetry.eq.6)) then
         do i = 1, num_gfuncs
            if (gfunc_name(i) .ne. 'chr'    .and.
     *          gfunc_name(i) .ne. 'mask'   .and.
     *          gfunc_name(i) .ne. 'r'      .and.
     *          gfunc_name(i) .ne. 'xphys'  .and.
     *          gfunc_name(i) .ne. 'yphys'  .and.
     *          gfunc_name(i) .ne. 'zphys'  .and.
     *          gfunc_name(i) .ne. 'cctk_x' .and.
     *          gfunc_name(i) .ne. 'cctk_y' .and.
     *          gfunc_name(i) .ne. 'cctk_z'      ) then
               if (ltrace2)
     *         write(*,*) 'grid_symm: field/i/ysym:',gfunc_name(i),i,
     *                            gfunc_ysym(i)
               call apply_symmetryy(q(gfunc_pointer(i,gi)), 
     *                           q(gr_y(gi)),
     *                           gfunc_ysym(i), 
     *                           gr_nx(gi),gr_ny(gi),gr_nz(gi))
            end if
         end do
      end if

      if (gr_minx(gi) .le. 0 .and.
     *    (assume_symmetry .eq. 1 .or. assume_symmetry.eq.6)) then
         do i = 1, num_gfuncs
            if (gfunc_name(i) .ne. 'chr'    .and.
     *          gfunc_name(i) .ne. 'mask'   .and.
     *          gfunc_name(i) .ne. 'r'      .and.
     *          gfunc_name(i) .ne. 'xphys'  .and.
     *          gfunc_name(i) .ne. 'yphys'  .and.
     *          gfunc_name(i) .ne. 'zphys'  .and.
     *          gfunc_name(i) .ne. 'cctk_x' .and.
     *          gfunc_name(i) .ne. 'cctk_y' .and.
     *          gfunc_name(i) .ne. 'cctk_z'      ) then
               if (ltrace2)
     *         write(*,*) 'grid_symm: field/i/xsym:',gfunc_name(i),i,
     *                            gfunc_xsym(i)
               call apply_symmetryx(q(gfunc_pointer(i,gi)), 
     *                           q(gr_x(gi)),
     *                           gfunc_xsym(i), 
     *                           gr_nx(gi),gr_ny(gi),gr_nz(gi))
            end if
         end do
      end if

 98   format('[',I3,'] ',A,3F15.9)
 99   format('[',I3,'] ',A,3I7)
      return
      end       ! END: grid_check_mask


      !-----------------------------------------------------------------
      !
      !
      !
      !-----------------------------------------------------------------
      subroutine set_interp_bounds( 
     &                 min_ip, max_ip, min_jp, max_jp,
     &                 min_kp, max_kp, min_ic, max_ic,
     &                 min_jc, max_jc, min_kc, max_kc,
     &                 nxp, nyp, nzp, length,
     &                 bi, bj, bk, refine_factor,
     &                 p_bnds_L, c_bnds_L, shp_L, length_L, b_L,
     &                 p_bnds_C, c_bnds_C, shp_C, length_C, b_C,
     &                 p_bnds_Q, c_bnds_Q, shp_Q, length_Q, b_Q)

        implicit none
        integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
        integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
        integer     nxp, nyp, nzp, length, bi, bj, bk
        integer     refine_factor
 
        integer     p_bnds_Q(6), c_bnds_Q(6), shp_Q(3), b_Q(3), length_Q
        integer     p_bnds_C(6), c_bnds_C(6), shp_C(3), b_C(3), length_C
        integer     p_bnds_L(6), c_bnds_L(6), shp_L(3), b_L(3), length_L

        ! local vars
        integer     i


  
        ! Keep these for linearly interpolated fields:
        p_bnds_L(1) = min_ip
        p_bnds_L(2) = max_ip
        p_bnds_L(3) = min_jp
        p_bnds_L(4) = max_jp
        p_bnds_L(5) = min_kp
        p_bnds_L(6) = max_kp

        c_bnds_L(1) = min_ic
        c_bnds_L(2) = max_ic
        c_bnds_L(3) = min_jc
        c_bnds_L(4) = max_jc
        c_bnds_L(5) = min_kc
        c_bnds_L(6) = max_kc

        shp_L(1) = p_bnds_L(2) - p_bnds_L(1) + 1
        shp_L(2) = p_bnds_L(4) - p_bnds_L(3) + 1
        shp_L(3) = p_bnds_L(6) - p_bnds_L(5) + 1
        !
        length_L = shp_L(1) * shp_L(2) * shp_L(3)
        b_L(1)   = 0
        b_L(2)   = 0
        b_L(3)   = 0

        ! Set initial values for cubic interpolation
        do i = 1, 6
          p_bnds_C(i) = p_bnds_L(i)
          c_bnds_C(i) = c_bnds_L(i)
        end do
  
        if (p_bnds_C(1) .gt.  1 ) then
           p_bnds_C(1) = p_bnds_C(1) - 1
        else
           c_bnds_C(1) = c_bnds_C(1) + refine_factor
        end if
        if (p_bnds_C(3) .gt.  1 ) then
           p_bnds_C(3) = p_bnds_C(3) - 1
        else
           c_bnds_C(3) = c_bnds_C(3) + refine_factor
        end if
        if (p_bnds_C(5) .gt.  1 ) then
           p_bnds_C(5) = p_bnds_C(5) - 1
        else
           c_bnds_C(5) = c_bnds_C(5) + refine_factor
        end if
        !
        if (p_bnds_C(2) .lt. nxp) then
           p_bnds_C(2) = p_bnds_C(2) + 1
        else
           c_bnds_C(2) = c_bnds_C(2) - refine_factor
        end if
        if (p_bnds_C(4) .lt. nyp)then
           p_bnds_C(4) = p_bnds_C(4) + 1
        else
           c_bnds_C(4) = c_bnds_C(4) - refine_factor
        end if
        if (p_bnds_C(6) .lt. nzp)then
           p_bnds_C(6) = p_bnds_C(6) + 1
        else
           c_bnds_C(6) = c_bnds_C(6) - refine_factor
        end if
  
        !
        ! Maximum amount that would be needed to store parent boundary
        ! data:
        !
        shp_C(1) = p_bnds_C(2) - p_bnds_C(1) + 1
        shp_C(2) = p_bnds_C(4) - p_bnds_C(3) + 1
        shp_C(3) = p_bnds_C(6) - p_bnds_C(5) + 1
        length_C =   shp_C(1) * shp_C(2) * shp_C(3)
        b_C(1)   = 1
        b_C(2)   = 1
        b_C(3)   = 1

        ! Set initial values for quintic interpolation
        ! The quintic stencil is 1 point larger in each direction than
        ! the cubic interpolation stencil.  Set the quintic limits
        ! first to the cubic limits, and then see if they have to be
        ! adjusted by one point.
        do i = 1, 6
          p_bnds_Q(i) = p_bnds_C(i)
          c_bnds_Q(i) = c_bnds_C(i)
        end do
  
        if (p_bnds_Q(1) .gt.  1 ) then
           p_bnds_Q(1) = p_bnds_Q(1) - 1
        else
           c_bnds_Q(1) = c_bnds_Q(1) + refine_factor
        end if
        if (p_bnds_Q(3) .gt.  1 ) then
           p_bnds_Q(3) = p_bnds_Q(3) - 1
        else
           c_bnds_Q(3) = c_bnds_Q(3) + refine_factor
        end if
        if (p_bnds_Q(5) .gt.  1 ) then
           p_bnds_Q(5) = p_bnds_Q(5) - 1
        else
           c_bnds_Q(5) = c_bnds_Q(5) + refine_factor
        end if
        !
        if (p_bnds_Q(2) .lt. nxp) then
           p_bnds_Q(2) = p_bnds_Q(2) + 1
        else
           c_bnds_Q(2) = c_bnds_Q(2) - refine_factor
        end if
        if (p_bnds_Q(4) .lt. nyp)then
           p_bnds_Q(4) = p_bnds_Q(4) + 1
        else
           c_bnds_Q(4) = c_bnds_Q(4) - refine_factor
        end if
        if (p_bnds_Q(6) .lt. nzp)then
           p_bnds_Q(6) = p_bnds_Q(6) + 1
        else
           c_bnds_Q(6) = c_bnds_Q(6) - refine_factor
        end if
  
        !
        ! Maximum amount that would be needed to store parent boundary
        ! data:
        !
        shp_Q(1) = p_bnds_Q(2) - p_bnds_Q(1) + 1
        shp_Q(2) = p_bnds_Q(4) - p_bnds_Q(3) + 1
        shp_Q(3) = p_bnds_Q(6) - p_bnds_Q(5) + 1
        length_Q =   shp_Q(1) * shp_Q(2) * shp_Q(3)
        b_Q(1)   = 2
        b_Q(2)   = 2
        b_Q(3)   = 2
 
        return
      end subroutine

      !-----------------------------------------------------------------
      !
      !
      !
      !-----------------------------------------------------------------
      subroutine load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nx, ny, nz, length,
     &                              bi, bj, bk,
     &                              p_bnds_C, c_bnds_C, shp_C, length_C,
     &                              b_C)

        implicit none
        integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
        integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
        integer     nx, ny, nz, length, bi, bj, bk
   
        integer     p_bnds_C(6), c_bnds_C(6), shp_C(3), length_C, b_C(3)

        ! local vars
        logical, parameter  ::  ltrace = .false.

        if (ltrace) print *,'*** Begin load_interp_const'

        min_ip = p_bnds_C(1)
        min_jp = p_bnds_C(3)
        min_kp = p_bnds_C(5)
        max_ip = p_bnds_C(2)
        max_jp = p_bnds_C(4)
        max_kp = p_bnds_C(6)
        !
        min_ic = c_bnds_C(1)
        min_jc = c_bnds_C(3)
        min_kc = c_bnds_C(5)
        max_ic = c_bnds_C(2)
        max_jc = c_bnds_C(4)
        max_kc = c_bnds_C(6)
        !
        nx    = shp_C(1)
        ny    = shp_C(2)
        nz    = shp_C(3)
        length = length_C
        !
        !
        bi     = b_C(1)
        bj     = b_C(2)
        bk     = b_C(3)

        if (ltrace) print *,'*** END   load_interp_const'

      return
      end

