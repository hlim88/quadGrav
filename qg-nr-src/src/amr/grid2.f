cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_mgsmooth:                                                            cc
cc            Smooth grid.                                                    cc
cc            Try not mess with values near AMR boundaries.                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_mgsmooth(gi,numattempts)
      implicit      none
      integer       gi, numattempts
      include      'grid.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer       nx,  ny,  nz
      integer       nxc, nyc, nzc
      integer       i, u_ptr, u_st1_ptr, u_st2_ptr, j
      real(kind=8)  mix, max, miy, may, miz, maz, time
      real(kind=8)  epsdis
      parameter   ( epsdis = 1.0d0 )
      integer       numdis  
      parameter   ( numdis = 3     )

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     vtrace
      parameter ( vtrace = .false. )


      if (ltrace) then
         write(*,*) 'grid_mgsmooth: grid: ',gi, numattempts
      end if

      call load_pointers(gi)

      nx  = gr_nx(gi)
      ny  = gr_ny(gi)
      nz  = gr_nz(gi)

      time= gr_t(gi)

      mix = gr_minx(gi)
      max = gr_maxx(gi)
      miy = gr_miny(gi)
      may = gr_maxy(gi)
      miz = gr_minz(gi)
      maz = gr_maxz(gi)

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            u_ptr     = gfunc_pointer(i,  gi)
            !
            if (ltrace ) write(*,*) 'grid_mgsmooth: ',gfunc_name(i)(1:7)
            if (vtrace) call field_out3d(q(u_ptr),time,
     *            'SMO_u_old', mix,max,miy,may,miz,maz, nx,ny,nz,myid)
            !
c           call load_scal1D(q(u_ptr),0.d0, nx*ny*nz)
c           !
c           if (numattempts .eq. 1) then
c              do j = 1, numdis
                  call apply_diss(q(u_ptr), q(gr_tmp), epsdis, nx,ny,nz)
c              end do
c           else
c              call field_smooth(q(u_ptr),q(gr_tmp),q(gr_chr),
c    *                 numattempts,nx,ny,nz)
c           end if
            !
            if (vtrace) call field_out3d(q(u_ptr),time,
     *            'SMO_u_new', mix,max,miy,may,miz,maz, nx,ny,nz,myid)
            !
         end if
      end do

      if (ltrace) then
         write(*,*) 'grid_mgsmooth: Finished.'
      end if

      return
      end      ! END: grid_mgsmooth

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_pstcompute:                                                          cc
cc            Compute post-cgc correction.                                    cc
cc            Transferred from the next coarser level should                  cc
cc            be u(2h) residing in storage at u_st1                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_pstcompute(gi)
      implicit      none
      integer       gi
      include      'grid.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer       nx,  ny,  nz
      integer       nxc, nyc, nzc
      integer       i, u_ptr, u_st1_ptr, u_st2_ptr
      real(kind=8)  mix, max, miy, may, miz, maz, time

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     vtrace
      parameter ( vtrace = .false. )


      if (ltrace) then
         write(*,*) 'grid_pstcompute: grid: ',gi
      end if

      call load_pointers(gi)

      nx  = gr_nx(gi)
      ny  = gr_ny(gi)
      nz  = gr_nz(gi)

      time= gr_t(gi)

      nxc = NINT( 0.5d0*(nx+1) )
      nyc = NINT( 0.5d0*(ny+1) )
      nzc = NINT( 0.5d0*(nz+1) )

      mix = gr_minx(gi)
      max = gr_maxx(gi)
      miy = gr_miny(gi)
      may = gr_maxy(gi)
      miz = gr_minz(gi)
      maz = gr_maxz(gi)

         !            h       h    h    2h     2h  h
         !           u    := u  + I  ( u   -  I   u  )
         !                         2h          h
         !
         !                       2h
         !          u_st1   <== u       (from grid_psttranfer)
         !

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            u_ptr     = gfunc_pointer(i,  gi)
            u_st1_ptr = gfunc_pointer(i+1,gi)
            u_st2_ptr = gfunc_pointer(i+2,gi)
            !
            if (vtrace)call field_out3d(q(u_ptr),time,'PST_u_2orig',
     *                    mix,max,miy,may,miz,maz, nx,ny,nz,myid)
            !
            call restrict(q(u_st2_ptr), q(u_ptr),nxc,nyc,nzc)
            if (vtrace) call field_out3d(q(u_st2_ptr),time,'PST_u_2H',
     *                    mix,max,miy,may,miz,maz, nxc,nyc,nzc,myid)
            !
            if (vtrace) call field_out3d(q(u_st1_ptr),time,'PSTu2Hcomp',
     *                     mix,max,miy,may,miz,maz, nxc,nyc,nzc,myid)
            call mat_mat_sub3d( q(u_st2_ptr), q(u_st1_ptr),q(u_st2_ptr),
     *                                                    nxc,nyc,nzc)
            if (vtrace) call field_out3d(q(u_st2_ptr),time,'PSTcorr2H',
     *                     mix,max,miy,may,miz,maz, nxc,nyc,nzc,myid)
            !
            ! Interpolate correction to fine grid:
            !
            call prolong( q(u_st1_ptr), q(u_st2_ptr), nxc,nyc,nzc)
            if (vtrace) call field_out3d(q(u_st1_ptr),time,'PSTcorrH',
     *                        mix,max,miy,may,miz,maz, nx,ny,nz,myid)
            !
            ! Zero out effects of domain boundary:
            !    (this step not actually necessary because these
            !     values get written over once the fields get synced)
            !  SLL (03/16/17): This appears to cause some small spikes
            !     at the corners of interior fine grids and so removed.
            !
            !call zerobounddeco( q(u_st1_ptr), q(gr_chr), 2, nx, ny, nz)
            if (vtrace) call field_out3d(q(u_st1_ptr),time,'PSTcorrrH',
     *                        mix,max,miy,may,miz,maz, nx,ny,nz,myid)
            !
            ! Add correction to field:
            !
            if (vtrace) call field_out3d(q(u_ptr),time,
     *            'PST_u_old', mix,max,miy,may,miz,maz, nx,ny,nz,myid)
             call mat_mat_add3d(q(u_ptr),q(u_ptr),q(u_st1_ptr),nx,ny,nz)
            if (vtrace) call field_out3d(q(u_ptr),time,
     *            'PST_u_new', mix,max,miy,may,miz,maz, nx,ny,nz,myid)
            !
         end if
      end do

      if (ltrace) then
         write(*,*) 'grid_pstcompute: Finished.'
      end if

      return
      end      ! END: grid_pstcompute

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_psttransfer_master:                                                  cc
cc                    This routine owns the parent grid.                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_psttransfer_master(child,parent)
      implicit      none
      integer       child, parent
      include      'grid.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i
      integer     nxt, nyt, nzt, length, nx,ny,nz
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     min_icc, min_jcc, min_kcc, owner
      integer       u_ptr, u_st1_ptr, u_st2_ptr
      real(kind=8)  mix, max, miy, may, miz, maz

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     vtrace
      parameter ( vtrace = .false. )


      call load_pointers(parent)

      owner = gr_own(child)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)
      
      !
      ! Dimensions of child:
      !
      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)
      
      mix = gr_minx(child)
      max = gr_maxx(child)
      miy = gr_miny(child)
      may = gr_maxy(child)
      miz = gr_minz(child)
      maz = gr_maxz(child)

      !
      ! Dimensions of virtual coarse grid (the u(2h) grid);
      !
      nx  = NINT( 0.5d0*(nxc+1) )
      ny  = NINT( 0.5d0*(nyc+1) )
      nz  = NINT( 0.5d0*(nzc+1) )

      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)
      
      ! 
      ! We're creating a 2h field on this grid, so we need
      ! the location of the intersection of the parent with
      ! this 2h child grid:
      !        (NB: interesction should be odd because
      !             child has to intersect parent at parent's grid points)
      ! 
      min_icc = NINT(0.5d0*(min_ic+1))
      min_jcc = NINT(0.5d0*(min_jc+1))
      min_kcc = NINT(0.5d0*(min_kc+1))

      ! 
      ! Dimensions of part to be copied:
      !
      nxt    = max_ip - min_ip + 1
      nyt    = max_jp - min_jp + 1
      nzt    = max_kp - min_kp + 1
      
      length = nxt * nyt * nzt
      
      if (length.gt.nxp*nyp*nzp) then 
         write(*,*) 'Problem of space in grid_psttransfer_master'
         call my_exit('Problem of space in grid_psttransfer_master')
      end if

      if (ltrace) then
       write(*,99)myid,'grid_psttransfer_master: grids:   ',child,parent
       write(*,99)myid,'grid_psttransfer_master: nx/y/z : ',nx, ny, nz
       write(*,99)myid,'grid_psttransfer_master: nx/y/zt: ',nxt,nyt,nzt
       write(*,99)myid,'grid_psttransfer_master: nx/y/zc: ',nxc,nyc,nzc
       write(*,99)myid,'grid_psttransfer_master: nx/y/zp: ',nxp,nyp,nzp
       write(*,99)myid,'grid_psttransfer_master: min_i/j/kp:',min_ip,
     *                                              min_jp,min_kp
       write(*,99)myid,'grid_psttransfer_master: min_i/j/kc:',min_ic,
     *                                              min_jc,min_kc
       write(*,99)myid,'grid_psttransfer_master: min_i/j/kcc:',min_icc,
     *                                              min_jcc, min_kcc
       write(*,99)myid,'grid_psttransfer_master: length:     ',length
      end if

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            u_ptr     = gfunc_pointer(i,  parent)
            u_st1_ptr = gfunc_pointer(i+1,parent)
            u_st2_ptr = gfunc_pointer(i+2,parent)
            !
            ! Copy data to tmp first:
            !
            call copy_overlap( q(gr_tmp),
     *                         q(u_ptr),
     *                         min_ip, min_jp, min_kp,
     *                         nxt,nyt,nzt,
     *                         nxp,nyp,nzp )
            if (vtrace)
     *         call field_out3d(q(gr_tmp),1.d0*parent, 'PSTT_u_par',
     *                   0.d0,1.d0,0.d0,1.d0,0.d0,1.d0,nxt,nyt,nzt,myid)
            !
            if (ltrace) then
               write(*,99) myid,'       master: sending correction'
               write(*,99) myid,'       master: to owner:   ', owner
               write(*,99) myid,'       master: length:     ',length
            end if
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_PST_DATA, MPI_COMM_WORLD, ierr)
            !
       end if
      end do

      if (ltrace) then
         write(*,99) myid,'grid_psttransfer_master: Finished.'
      end if

 99   format('[',I3,'] ',A,3I7)

      return
      end      ! END: grid_psttransfer_master

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_psttransfer_slave:                                                   cc
cc                    This routine owns the child grid.                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_psttransfer_slave(child,parent)
      implicit      none
      integer       child, parent
      include      'grid.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i, phi_parent
      integer     nxt, nyt, nzt, length, nx,ny,nz
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     min_icc, min_jcc, min_kcc, owner
      integer       u_ptr, u_st1_ptr, u_st2_ptr
      real(kind=8)  mix, max, miy, may, miz, maz

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     vtrace
      parameter ( vtrace = .false. )


      call load_pointers(child)

      owner = gr_own(parent)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)
      
      !
      ! Dimensions of child:
      !
      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)
      
      mix = gr_minx(child)
      max = gr_maxx(child)
      miy = gr_miny(child)
      may = gr_maxy(child)
      miz = gr_minz(child)
      maz = gr_maxz(child)

      !
      ! Dimensions of virtual coarse grid (the u(2h) grid);
      !
      nx  = NINT( 0.5d0*(nxc+1) )
      ny  = NINT( 0.5d0*(nyc+1) )
      nz  = NINT( 0.5d0*(nzc+1) )

      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)
      
      ! 
      ! We're creating a 2h field on this grid, so we need
      ! the location of the intersection of the parent with
      ! this 2h child grid:
      !        (NB: interesction should be odd because
      !             child has to intersect parent at parent's grid points)
      ! 
      min_icc = NINT(0.5d0*(min_ic+1))
      min_jcc = NINT(0.5d0*(min_jc+1))
      min_kcc = NINT(0.5d0*(min_kc+1))

      ! 
      ! Dimensions of part to be copied:
      !
      nxt    = max_ip - min_ip + 1
      nyt    = max_jp - min_jp + 1
      nzt    = max_kp - min_kp + 1
      
      length = nxt * nyt * nzt
      
      if (length.gt.nxc*nyc*nzc) then 
         write(*,*) 'Problem of space in grid_psttransfer_slave'
         call my_exit('Problem of space in grid_psttransfer_slave')
      end if

      if (ltrace) then
       write(*,99)myid,'grid_psttransfer_slave: grids:    ',child,parent
       write(*,99)myid,'grid_psttransfer_slave: nx/y/z :  ',nx, ny, nz
       write(*,99)myid,'grid_psttransfer_slave: nx/y/zt:  ',nxt,nyt,nzt
       write(*,99)myid,'grid_psttransfer_slave: nx/y/zc:  ',nxc,nyc,nzc
       write(*,99)myid,'grid_psttransfer_slave: nx/y/zp:  ',nxp,nyp,nzp
       write(*,99)myid,'grid_psttransfer_slave: min_i/j/kc:',min_ic,
     *                                              min_jc,  min_kc
       write(*,99)myid,'grid_psttransfer_slave: min_i/j/kcc:',min_icc,
     *                                              min_jcc, min_kcc
      end if

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            u_ptr     = gfunc_pointer(i,  child)
            u_st1_ptr = gfunc_pointer(i+1,child)
            u_st2_ptr = gfunc_pointer(i+2,child)
            if (ltrace) then
               write(*,99) myid,'       slave: receiving correction'
               write(*,99) myid,'       slave: from owner: ', owner
               write(*,99) myid,'       slave: length:     ',length
            end if
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_PST_DATA, MPI_COMM_WORLD, status,ierr)
  
            !
            ! Transfer into appropriate place in phi_rk1:
            !
            call copy_toplace( q(u_st1_ptr),
     *                         q(gr_tmp),
     *                         min_icc, min_jcc, min_kcc,
     *                         nx, ny, nz,
     *                         nxt,nyt,nzt )
            !
            if (vtrace)
     *         call field_out3d(q(u_st1_ptr),0.d0, 'PSTT_u_2H',
     *                    mix,max,miy,may,miz,maz,nx,ny,nz,myid)
            !
         end if
      end do


      if (ltrace) then
         write(*,99) myid,'grid_psttransfer_slave: Finished.'
      end if

 99   format('[',I3,'] ',A,3I7)

      return
      end      ! END: grid_psttransfer_slave

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_psttransfer:                                                         cc
cc                    This routine gets called for the different              cc
cc                    pairings of child and parent which overlap.             cc
cc                    The goal is to transfer down to the child               cc
cc                    the field(s) u_2H. Once all such pairings finish,       cc
cc                    then grid_pstcompute() processes the complete u_2H      cc
cc                    field to correct u_H.                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_psttransfer(child,parent)
      implicit      none
      integer       child, parent
      include      'grid.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i, phi_parent
      integer     nxt, nyt, nzt, length, nx,ny,nz
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     min_icc, min_jcc, min_kcc
      integer       u_ptr, u_st1_ptr, u_st2_ptr
      real(kind=8)  mix, max, miy, may, miz, maz

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     vtrace
      parameter ( vtrace = .false. )


      call load_pointers(parent)
      !phi_parent = gr_phi

      call load_pointers(child)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)
      
      !
      ! Dimensions of child:
      !
      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)
      
      mix = gr_minx(child)
      max = gr_maxx(child)
      miy = gr_miny(child)
      may = gr_maxy(child)
      miz = gr_minz(child)
      maz = gr_maxz(child)

      !
      ! Dimensions of virtual coarse grid (the u(2h) grid);
      !
      nx  = NINT( 0.5d0*(nxc+1) )
      ny  = NINT( 0.5d0*(nyc+1) )
      nz  = NINT( 0.5d0*(nzc+1) )

      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)
      
      ! 
      ! We're creating a 2h field on this grid, so we need
      ! the location of the intersection of the parent with
      ! this 2h child grid:
      !        (NB: interesction should be odd because
      !             child has to intersect parent at parent's grid points)
      ! 
      min_icc = NINT(0.5d0*(min_ic+1))
      min_jcc = NINT(0.5d0*(min_jc+1))
      min_kcc = NINT(0.5d0*(min_kc+1))

      ! 
      ! Dimensions of part to be copied:
      !
      nxt    = max_ip - min_ip + 1
      nyt    = max_jp - min_jp + 1
      nzt    = max_kp - min_kp + 1
      
      length = nxt * nyt * nzt
      
      if (length.gt.nxc*nyc*nzc) then 
         write(*,*) 'Problem of space in grid_psttransfer'
         call my_exit('Problem of space in grid_psttransfer')
      end if

      if (ltrace) then
         write(*,*)'grid_psttransfer: grids:      ',child, parent
         write(*,*)'grid_psttransfer: nx/y/z :    ',nx, ny, nz
         write(*,*)'grid_psttransfer: nx/y/zt:    ',nxt,nyt,nzt
         write(*,*)'grid_psttransfer: nx/y/zc:    ',nxc,nyc,nzc
         write(*,*)'grid_psttransfer: nx/y/zp:    ',nxp,nyp,nzp
         write(*,*)'grid_psttransfer: min_i/j/kp:', min_ip,
     *                                              min_jp,min_kp
         write(*,*)'grid_psttransfer: min_i/j/kc: ',min_ic,min_jc,min_kc
         write(*,*)'grid_psttransfer: min_i/j/kcc:',min_icc,min_jcc,
     *                                              min_kcc
      end if

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            u_ptr     = gfunc_pointer(i,  child)
            u_st1_ptr = gfunc_pointer(i+1,child)
            !u_st2_ptr = gfunc_pointer(i+2,child)
            !
            ! Copy data to tmp first:
            !
            call copy_overlap( q(gr_tmp),
     *                         q(gfunc_pointer(i,parent)),
     *                         min_ip, min_jp, min_kp,
     *                         nxt,nyt,nzt,
     *                         nxp,nyp,nzp )
            if (vtrace)
     *      call field_out3d(q(gfunc_pointer(i,parent)),gr_t(parent),
     .                       'PSTT_u_Ppar',
     *                gr_minx(parent),gr_maxx(parent),
     *                gr_miny(parent),gr_maxy(parent),
     *                gr_minz(parent),gr_maxz(parent),
     *                nxp,nyp,nzp,myid)
            if (vtrace)
     *      call field_out3d(q(gr_tmp),gr_t(parent), 'PSTT_u_par',
     *                q(gr_x(parent)+min_ip-1),
     *                q(gr_x(parent)+max_ip-1),
     *                q(gr_y(parent)+min_jp-1),
     *                q(gr_y(parent)+max_jp-1),
     *                q(gr_z(parent)+min_kp-1),
     *                q(gr_z(parent)+max_kp-1),
!    *                gr_minx(parent),gr_maxx(parent),
!    *                gr_miny(parent),gr_maxy(parent),
!    *                gr_minz(parent),gr_maxz(parent),
     *                nxt,nyt,nzt,myid)
            !
            ! Transfer into appropriate place in phi_rk1:
            !
            if (vtrace)
     *      call field_out3d(q(u_st1_ptr),gr_t(child),'PSTT_u_2Hprev',
     *                    mix,max,miy,may,miz,maz,nx,ny,nz,myid)
            call copy_toplace(
     *                       q(u_st1_ptr),
     *                       q(gr_tmp),
     *                       min_icc, min_jcc, min_kcc,
     *                       nx, ny, nz,
     *                       nxt,nyt,nzt )

            if (vtrace)
     *      call field_out3d(q(u_st1_ptr),gr_t(child),'PSTT_u_2H',
     *                    mix,max,miy,may,miz,maz,nx,ny,nz,myid)
            !
         end if
      end do


      if (ltrace) then
         write(*,*) 'grid_psttransfer: Finished.'
      end if

      return
      end      ! END: grid_psttransfer

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_cgctransfer_slave:                                                   cc
cc                 We own the child but not the parent                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_cgctransfer_slave(child, parent)
      implicit      none
      integer       child, parent
      include      'grid.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer     nxs,nys,nzs
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i, nx,ny,nz
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     min_icc, min_jcc, min_kcc
      integer       u_ptr, u_st1_ptr, u_st2_ptr, u_rhs_ptr
      integer     length, owner
      real(kind=8)  mixc, maxc, miyc, mayc, mizc, mazc
      real(kind=8)  mixp, maxp, miyp, mayp, mizp, mazp
      real(kind=8)  time

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     vtrace
      parameter ( vtrace = .false. )


      if (ltrace) then
         write(*,99)myid,'grid_cgctransfer_slavechild/prnt',child,parent
      end if

      owner = gr_own(parent)

      call load_pointers(child)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent) 
      nzp = gr_nz(parent)

      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      time= gr_t(child)

      !
      ! Dimensions of virtual coarse grid (the u(2h) grid);
      !
      nx  = NINT( 0.5d0*(nxc+1) )
      ny  = NINT( 0.5d0*(nyc+1) )
      nz  = NINT( 0.5d0*(nzc+1) )

      mixc = gr_minx(child)
      maxc = gr_maxx(child)
      miyc = gr_miny(child)
      mayc = gr_maxy(child)
      mizc = gr_minz(child)
      mazc = gr_maxz(child)

      mixp = gr_minx(parent)
      maxp = gr_maxx(parent)
      miyp = gr_miny(parent)
      mayp = gr_maxy(parent)
      mizp = gr_minz(parent)
      mazp = gr_maxz(parent)

      !
      ! Find indices which bound intersection region:
      !
      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      !
      ! Overlap region on the parent
      !    but with resolution of the parent grid:
      !
      nxs = max_ip - min_ip + 1
      nys = max_jp - min_jp + 1
      nzs = max_kp - min_kp + 1

      ! 
      ! We're creating a 2h field on this grid, so we need
      ! the location of the intersection of the parent with
      ! this 2h child grid:
      !        (NB: interesction should be odd because
      !             child has to intersect parent at parent's grid points)
      ! 
      min_icc = NINT(0.5d0*(min_ic+1))
      min_jcc = NINT(0.5d0*(min_jc+1))
      min_kcc = NINT(0.5d0*(min_kc+1))
      !
      !length = nx * ny * nz
      length = nxs * nys * nzs
      !
      if (ltrace) then
         write(*,99)myid,'grid_cgctransfer_slave:  length=',length
         write(*,99)myid,'grid_cgctransfer_slave:  nx/y/z =',nx,ny,nz
         write(*,99)myid,'grid_cgctransfer_slave:  nx/y/zp=',nxp,nyp,nzp
         write(*,99)myid,'grid_cgctransfer_slave:  nx/y/zs=',nxp,nys,nzs
      end if

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            u_ptr     = gfunc_pointer(i,  child)
            u_st1_ptr = gfunc_pointer(i+1,child)
            u_st2_ptr = gfunc_pointer(i+2,child)
            !
            ! Copy data from child to tmp first:
            !
            call copy_overlap( q(gr_tmp),
     *                         q(u_st1_ptr),
     *                         min_icc, min_jcc, min_kcc,
     *                         nxs,nys,nzs,
     *                         nx ,ny ,nz  )
            call copy_overlap( q(gr_flag),
     *                         q(u_st2_ptr),
     *                         min_icc, min_jcc, min_kcc,
     *                         nxs,nys,nzs,
     *                         nx ,ny ,nz  )
            !
            ! Send to owner of parent:
            !
            if (ltrace) then
               write(*,99) myid,' slave:  sending correction'
               write(*,99) myid,' slave:   to owner: ', owner
            end if
            !
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_CGC_DATA, MPI_COMM_WORLD, ierr)
            call MPI_Send(q(gr_flag), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_CGC_DATA, MPI_COMM_WORLD, ierr)
            !
         else if ( .false. ) then
!        else if ( i .eq.  99 .or. i .eq. 102 .or. i .eq. 105 .or.
!    *             i .eq. 108 .or. i .eq. 111 .or. i .eq. 114 .or.
!    *             i .eq. 135 .or. i .eq. 138 .or. i .eq. 141 .or.
!    *             i .eq. 144 .or. i .eq. 147 .or. i .eq. 150 .or.
!    *             i .eq. 153 .or. i .eq. 156 .or. i .eq. 159 .or.
!    *             i .eq. 162 .or. i .eq. 165 .or. i .eq. 168 .or.
!    *             i .eq. 171 .or. i .eq. 174 .or. i .eq. 177 .or.
!    *             i .eq. 180 .or. i .eq. 183 .or. i .eq. 186     )then
            !
            !u_ptr     = gr_g11_rk3
            u_ptr     = gfunc_pointer(i+2,  child)
           call copy_overlap(q(gr_tmp),q(u_ptr),min_icc,min_jcc,min_kcc,
     *                         nxs,nys,nzs, nx ,ny ,nz  )
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_CGC_DATA, MPI_COMM_WORLD, ierr)
            !
         end if
      end do

 99   format('[',I3,'] ',A,3I7)
      if (ltrace) then
         write(*,99) myid,'grid_cgctransfer_slave: Finished.'
      end if

      return
      end      ! END: grid_cgctransfer_slave

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_cgctransfer_master:                                                  cc
cc                               We own parent.                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_cgctransfer_master(child, parent)
      implicit      none
      integer       child, parent
      include      'grid.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer     nxs,nys,nzs
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i, nx,ny,nz
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     min_icc, min_jcc, min_kcc
      integer       u_ptr, u_st1_ptr, u_st2_ptr, u_rhs_ptr
      integer     length, owner
      real(kind=8)  mixc, maxc, miyc, mayc, mizc, mazc
      real(kind=8)  mixp, maxp, miyp, mayp, mizp, mazp
      real(kind=8)  time


      logical     ltrace
      parameter ( ltrace = .false. )
      logical     vtrace
      parameter ( vtrace = .false. )


      if (ltrace) then
        write(*,99)myid,'grid_cgctransfer_master:ch/par:',child,parent
      end if

      owner = gr_own(child)

      time  = gr_t(child)

      call load_pointers(parent)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent) 
      nzp = gr_nz(parent)

      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      !
      ! Dimensions of virtual coarse grid (the u(2h) grid);
      !
      nx  = NINT( 0.5d0*(nxc+1) )
      ny  = NINT( 0.5d0*(nyc+1) )
      nz  = NINT( 0.5d0*(nzc+1) )

      mixc = gr_minx(child)
      maxc = gr_maxx(child)
      miyc = gr_miny(child)
      mayc = gr_maxy(child)
      mizc = gr_minz(child)
      mazc = gr_maxz(child)

      mixp = gr_minx(parent)
      maxp = gr_maxx(parent)
      miyp = gr_miny(parent)
      mayp = gr_maxy(parent)
      mizp = gr_minz(parent)
      mazp = gr_maxz(parent)

      !
      ! Find indices which bound intersection region:
      !
      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      !
      ! Overlap region on the parent
      !    but with resolution of the parent grid:
      !
      nxs = max_ip - min_ip + 1
      nys = max_jp - min_jp + 1
      nzs = max_kp - min_kp + 1

      ! 
      ! We're creating a 2h field on this grid, so we need
      ! the location of the intersection of the parent with
      ! this 2h child grid:
      !        (NB: interesction should be odd because
      !             child has to intersect parent at parent's grid points)
      ! 
      min_icc = NINT(0.5d0*(min_ic+1))
      min_jcc = NINT(0.5d0*(min_jc+1))
      min_kcc = NINT(0.5d0*(min_kc+1))

      !length = nx * ny * nz
      length = nxs * nys * nzs
      if (length .gt. nxp*nyp*nzp) then
         write(*,*) 'grid_cgctransfer_master: Not enough mem'
         write(*,*) 'grid_cgctransfer_master: length=',length
         write(*,*) 'grid_cgctransfer_master: nx/y/z =',nx,ny,nz
         write(*,*) 'grid_cgctransfer_master: nx/y/zp=',nxp,nyp,nzp
         call my_exit('grid_cgctransfer_master: Not enough mem')
      end if
      if (ltrace) then
         write(*,99)myid,'grid_cgctransfer_master: length=',length
         write(*,99)myid,'grid_cgctransfer_master: nx/y/z =',nx,ny,nz
         write(*,99)myid,'grid_cgctransfer_master: nx/y/zp=',nxp,nyp,nzp
         write(*,99)myid,'grid_cgctransfer_master: nx/y/zs=',nxs,nys,nzs
      end if

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            u_ptr     = gfunc_pointer(i,  parent)
            u_st1_ptr = gfunc_pointer(i+1,parent)
            u_st2_ptr = gfunc_pointer(i+2,parent)
            u_rhs_ptr = gfunc_pointer(i+3,parent)
            !
            if (ltrace) then
               write(*,99) myid,' master: receiving correction'
               write(*,99) myid,' master: from owner: ', owner
            end if
            !
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_CGC_DATA, MPI_COMM_WORLD, status,ierr)
            call MPI_Recv(q(gr_flag), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_CGC_DATA, MPI_COMM_WORLD, status,ierr)
            !
            ! Transfer into appropriate place on parent:
            !
            if (vtrace)
     *         call field_out3d(q(u_ptr),time, 'CGCT_old_par',
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            call copy_toplace( q(u_ptr),
     *                         q(gr_tmp),
     *                         min_ip, min_jp, min_kp,
     *                         nxp,nyp,nzp,
     *                         nxs,nys,nzs )
            if (vtrace)
     *         call field_out3d(q(u_ptr),time, 'CGCT_new_par',
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            !
            if (vtrace)
     *         call field_out3d(q(u_rhs_ptr),time, 'CGCT_old_rho',
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            call copy_toplace( q(u_rhs_ptr),
     *                         q(gr_flag),
     *                         min_ip, min_jp, min_kp,
     *                         nxp,nyp,nzp,
     *                         nxs,nys,nzs )
            if (vtrace)
     *         call field_out3d(q(u_rhs_ptr),time, 'CGCT_new_rho',
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            !
         else if ( .false. ) then
!        else if ( i .eq.  99 .or. i .eq. 102 .or. i .eq. 105 .or.
!    *             i .eq. 108 .or. i .eq. 111 .or. i .eq. 114 .or.
!    *             i .eq. 135 .or. i .eq. 138 .or. i .eq. 141 .or.
!    *             i .eq. 144 .or. i .eq. 147 .or. i .eq. 150 .or.
!    *             i .eq. 153 .or. i .eq. 156 .or. i .eq. 159 .or.
!    *             i .eq. 162 .or. i .eq. 165 .or. i .eq. 168 .or.
!    *             i .eq. 171 .or. i .eq. 174 .or. i .eq. 177 .or.
!    *             i .eq. 180 .or. i .eq. 183 .or. i .eq. 186     )then
            !
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_CGC_DATA, MPI_COMM_WORLD, status,ierr)
            u_ptr     = gfunc_pointer(i,    parent)
            call copy_toplace( q(u_ptr),q(gr_tmp),min_ip, min_jp,min_kp,
     *                         nxp,nyp,nzp, nxs,nys,nzs )
            if (vtrace) call
     *         field_out3d(q(u_ptr),gr_t(parent),'CGCT_'//gfunc_name(i),
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            !
         end if
      end do

      if (ltrace) then
         write(*,99)myid, 'grid_cgctransfer_master: Finished.'
 99      format('[',I3,'] ',A,3I7)
      end if

      return
      end      ! END: grid_cgctransfer_master

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_cgctransfer:                                                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_cgctransfer(child, parent)
      implicit      none
      integer       child, parent
      include      'grid.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      integer     nxs,nys,nzs
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i, nx,ny,nz
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     min_icc, min_jcc, min_kcc
      integer       u_ptr, u_st1_ptr, u_st2_ptr, u_rhs_ptr, u_tmp
      integer       length, parent_chr
      real(kind=8)  mixc, maxc, miyc, mayc, mizc, mazc
      real(kind=8)  mixp, maxp, miyp, mayp, mizp, mazp
      real(kind=8)  time, hp


      logical     ltrace
      parameter ( ltrace = .false. )
      logical     vtrace
      parameter ( vtrace = .false. )
      logical     vtrace2
      parameter ( vtrace2 = .false. )


      if (ltrace) then
         write(*,99)myid,'grid_cgctransfer:child, parent: ',child,parent
      end if

      call load_pointers(parent)
      parent_chr = gr_chr
      call load_pointers(child)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent) 
      nzp = gr_nz(parent)
      hp  = gr_h(parent)

      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      time= gr_t(child)

      !
      ! Dimensions of virtual coarse grid (the u(2h) grid);
      !
      nx  = NINT( 0.5d0*(nxc+1) )
      ny  = NINT( 0.5d0*(nyc+1) )
      nz  = NINT( 0.5d0*(nzc+1) )

      mixc = gr_minx(child)
      maxc = gr_maxx(child)
      miyc = gr_miny(child)
      mayc = gr_maxy(child)
      mizc = gr_minz(child)
      mazc = gr_maxz(child)

      mixp = gr_minx(parent)
      maxp = gr_maxx(parent)
      miyp = gr_miny(parent)
      mayp = gr_maxy(parent)
      mizp = gr_minz(parent)
      mazp = gr_maxz(parent)

      !
      ! Find indices which bound intersection region:
      !
      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      !
      ! Overlap region on the parent
      !    but with resolution of the parent grid:
      !
      nxs = max_ip - min_ip + 1
      nys = max_jp - min_jp + 1
      nzs = max_kp - min_kp + 1

      ! 
      ! We're creating a 2h field on this grid, so we need
      ! the location of the intersection of the parent with
      ! this 2h child grid:
      !        (NB: interesction should be odd because
      !             child has to intersect parent at parent's grid points)
      ! 
      min_icc = NINT(0.5d0*(min_ic+1))
      min_jcc = NINT(0.5d0*(min_jc+1))
      min_kcc = NINT(0.5d0*(min_kc+1))

      if (ltrace) then
         write(*,99)myid,'grid_cgctransfer: length=',length
         write(*,99)myid,'grid_cgctransfer: nx/y/z =',nx,ny,nz
         write(*,99)myid,'grid_cgctransfer: nx/y/zp=',nxp,nyp,nzp
         write(*,99)myid,'grid_cgctransfer: nx/y/zc=',nxc,nyc,nzc
         write(*,99)myid,'grid_cgctransfer: nx/y/zs=',nxs,nys,nzs
         write(*,99)myid,'grid_cgctransfer: min_i/j/kcc=',min_icc,
     .                  min_jcc, min_kcc
         write(*,99)myid,'grid_cgctransfer: min_i/j/kp =',min_ip,
     .                  min_jp, min_kp
      end if

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            u_ptr     = gfunc_pointer(i,  child)
            u_st1_ptr = gfunc_pointer(i+1,child)
            u_st2_ptr = gfunc_pointer(i+2,child)
            u_rhs_ptr = gfunc_pointer(i+3,child)
            !
            !
            ! Copy data from child to tmp first:
            !
            call copy_overlap( q(gr_tmp),
     *                         q(u_st1_ptr),
     *                         min_icc, min_jcc, min_kcc,
     *                         nxs,nys,nzs,
     *                         nx ,ny ,nz  )
            call copy_overlap( q(gr_flag),
     *                         q(u_st2_ptr),
     *                         min_icc, min_jcc, min_kcc,
     *                         nxs,nys,nzs,
     *                         nx ,ny ,nz  )
            !
            ! Transfer into appropriate place on parent:
            !
            !
            u_ptr     = gfunc_pointer(i,  parent)
            u_st1_ptr = gfunc_pointer(i+1,parent)
            u_st2_ptr = gfunc_pointer(i+2,parent)
            u_rhs_ptr = gfunc_pointer(i+3,parent)
            if (vtrace)
     *         call field_out3d(q(u_ptr),time, 'CGCT_old_par',
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            call copy_toplace( q(u_ptr),
     *                         q(gr_tmp),
     *                         min_ip, min_jp, min_kp,
     *                         nxp,nyp,nzp,
     *                         nxs,nys,nzs )
            if (vtrace)
     *         call field_out3d(q(u_ptr),time, 'CGCT_new_par',
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            !
            if (vtrace)
     *         call field_out3d(q(u_rhs_ptr),time, 'CGCT_old_rho',
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            if (vtrace)
     *         call field_out3d(q(gr_flag),  time, 'CGCT_flag_rho',
     *                   mixp+hp*(min_ip-1),mixp+hp*(max_ip-1),
     *                   miyp+hp*(min_jp-1),miyp+hp*(max_jp-1),
     *                   mizp+hp*(min_kp-1),mizp+hp*(max_kp-1),
     *                   nxs,nys,nzs,myid)
        if(ltrace)call field_dump_infoB(q(u_rhs_ptr),nxp,nyp,nzp,'rhs' )
        if(ltrace)call field_dump_infoB(q(gr_flag),  nxs,nys,nzs,'flag')
            call copy_toplace( q(u_rhs_ptr),
     *                         q(gr_flag),
     *                         min_ip, min_jp, min_kp,
     *                         nxp,nyp,nzp,
     *                         nxs,nys,nzs )
         if(ltrace)call field_dump_infoB(q(u_rhs_ptr),nxp,nyp,nzp,'rhs')
            if (vtrace)
     *         call field_out3d(q(u_rhs_ptr),time, 'CGCT_new_rho',
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
c          Not needed (soon to be deleted):
c           call zerobounddeco( q(u_rhs_ptr), q(parent_chr),
c    *                          1, nxp, nyp, nzp)
c           if (vtrace)
c    *         call field_out3d(q(u_rhs_ptr),time, 'CGCT_Nnew_rho',
c    *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            !
         else if ( .false. ) then
!        else if ( i .eq.  99 .or. i .eq. 102 .or. i .eq. 105 .or.
!    *             i .eq. 108 .or. i .eq. 111 .or. i .eq. 114 .or.
!    *             i .eq. 135 .or. i .eq. 138 .or. i .eq. 141 .or.
!    *             i .eq. 144 .or. i .eq. 147 .or. i .eq. 150 .or.
!    *             i .eq. 153 .or. i .eq. 156 .or. i .eq. 159 .or.
!    *             i .eq. 162 .or. i .eq. 165 .or. i .eq. 168 .or.
!    *             i .eq. 171 .or. i .eq. 174 .or. i .eq. 177 .or.
!    *             i .eq. 180 .or. i .eq. 183 .or. i .eq. 186     )then
            !
            u_ptr     = gfunc_pointer(i+2,  child)
           call copy_overlap(q(gr_tmp),q(u_ptr),min_icc,min_jcc,min_kcc,
     *                         nxs,nys,nzs, nx ,ny ,nz  )
            if (vtrace2) call
     *        field_out3d(q(u_ptr),gr_t(parent),'CGCT_R'//gfunc_name(i),
     *                   mixc,maxc,miyc,mayc,mizc,mazc,nx, ny, nz, myid)
            u_ptr     = gfunc_pointer(i,    parent)
            call copy_toplace( q(u_ptr),q(gr_tmp),min_ip, min_jp,min_kp,
     *                         nxp,nyp,nzp, nxs,nys,nzs )
            if (vtrace2) call
     *         field_out3d(q(u_ptr),gr_t(parent),'CGCT_'//gfunc_name(i),
     *                   mixp,maxp,miyp,mayp,mizp,mazp,nxp,nyp,nzp,myid)
            !
         end if
      end do

      ! Re-calculate the residual
      !    ...this is not needed for the solve, but otherwise
      !       the st1 arrays have coarse grid data and if looking
      !       at the vtracing, there looks to be a bug
      call grid_get_resid(child)

      if (ltrace) then
         write(*,99) myid, 'grid_cgctransfer: Finished.'
 99      format('[',I4,'] ',A,3I7)
      end if

      return
      end      ! END: grid_cgctransfer

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_cgccompute:                                                          cc
cc            Compute coarse grid correction.                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_cgccompute(gi)
      implicit      none
      integer       gi
      include      'grid.inc'
      include      'glob.inc'
      include      'mpif.h'
      include      'mpi_stuff.inc'
      real(kind=8)  hf, hc
      integer       nx,  ny,  nz
      integer       nxc, nyc, nzc
      real(kind=8)  mix, max, miy, may, miz, maz, time
      integer       i, u_ptr, u_st1_ptr, u_st2_ptr, u_rhs_ptr
      integer       u_rk1_ptr, u_rk3_ptr
      character*4   tmpchar

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     vtrace
      parameter ( vtrace = .false. )


      if (ltrace) then
         write(*,99)myid,'grid_cgccompute: grid: ',gi
      end if
      if (vtrace) call int2str(gi,tmpchar)
      !if (.true.) call int2str(gi,tmpchar)

      call load_pointers(gi)

      nx  = gr_nx(gi)
      ny  = gr_ny(gi)
      nz  = gr_nz(gi)

      hf  = gr_h(gi)
      hc  = 2.d0*hf

      time= gr_t(gi)

      mix = gr_minx(gi)
      max = gr_maxx(gi)
      miy = gr_miny(gi)
      may = gr_maxy(gi)
      miz = gr_minz(gi)
      maz = gr_maxz(gi)

      nxc = NINT( 0.5d0*(nx+1) )
      nyc = NINT( 0.5d0*(ny+1) )
      nzc = NINT( 0.5d0*(nz+1) )

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

      ! Use the associated storage arrays as temporary memory:
      !
      !                     2h
      ! u_st1       <===   u
      !
      !                       2h
      ! u_st2       <===   rhs
      !
            !
            ! (1) Compute residual of u^h  into the "u_st1" fields
            ! (2) Restrict            u^h    to the "u_st2" fields
            ! (3) Compute lop      of u^2h into the "u_st2" fields
            !
      if(ltrace)write(*,99)myid,'grid_cgccompute: call grid_get_resid()'
      call grid_get_resid(gi)
      !
      if(ltrace)write(*,99)myid,'grid_cgccompute:restrict fields/resids'
      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            u_ptr     = gfunc_pointer(i,  gi)
            u_st1_ptr = gfunc_pointer(i+1,gi)
            u_st2_ptr = gfunc_pointer(i+2,gi)
            !
            call restrict(  q(u_st2_ptr), q(u_st1_ptr), nxc,nyc,nzc)
            call restrict(  q(u_st1_ptr), q(u_ptr),     nxc,nyc,nzc)
            !
            if (vtrace) then
            !if (.true.) then
               ! this code section seems to need to be called for
               ! some strange bug:
               call field_out3d(q(u_st2_ptr),time,'CGC_st2'//tmpchar,
     *                       mix,max,miy,may,miz,maz,nxc,nyc,nzc,myid)
               call field_out3d(q(u_st1_ptr),time, 'CGC_st1'//tmpchar,
     *                       mix,max,miy,may,miz,maz,nxc,nyc,nzc,myid)
            end if
         end if
      end do
      !
      if (ltrace) write(*,99)myid,'grid_cgccompute: call grid_lop()'
      call grid_lop(gi)
      !
      if (ltrace) write(*,99)myid,'grid_cgccompute: computing new rhs'
      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            u_ptr     = gfunc_pointer(i,  gi)
            u_st1_ptr = gfunc_pointer(i+1,gi)
            u_st2_ptr = gfunc_pointer(i+2,gi)
            !
            ! Mask out DD boundaries where data is on a 2h domain
            ! and chr is on a h domain:
            !
            if (vtrace) then
               call field_out3d(q(u_st2_ptr),time,'CGC_rhs_2H'//tmpchar,
     *                       mix,max,miy,may,miz,maz,nxc,nyc,nzc,myid)
               call field_out3d(q(u_st1_ptr),time, 'CGC_u_2H'//tmpchar,
     *                       mix,max,miy,may,miz,maz,nxc,nyc,nzc,myid)
            end if
            !
            if (ltrace) write(*,99)myid,'grid_cgccompute: masking'
            call maskLargeh2h(q(u_st2_ptr),q(gr_chr),nxc,nyc,nzc,
     *                                                      bound_width)
            if (ltrace) write(*,99)myid,'grid_cgccompute: masking'
            call maskLargeh2h(q(u_st1_ptr),q(gr_chr),nxc,nyc,nzc,
     *                                                      bound_width)
            !
            if (vtrace) then
c              call int2str(i,tmpchar)
               call field_out3d(q(u_st2_ptr),time,'CGC_rhs_2H'//tmpchar,
     *                       mix,max,miy,may,miz,maz,nxc,nyc,nzc,myid)
               call field_out3d(q(u_st1_ptr),time, 'CGC_u_2H'//tmpchar,
     *                       mix,max,miy,may,miz,maz,nxc,nyc,nzc,myid)
            end if
            !
         else if ( .false. ) then
!        else if ( i .eq.  99 .or. i .eq. 102 .or. i .eq. 105 .or.
!    *             i .eq. 108 .or. i .eq. 111 .or. i .eq. 114 .or.
!    *             i .eq. 135 .or. i .eq. 138 .or. i .eq. 141 .or.
!    *             i .eq. 144 .or. i .eq. 147 .or. i .eq. 150 .or.
!    *             i .eq. 153 .or. i .eq. 156 .or. i .eq. 159 .or.
!    *             i .eq. 162 .or. i .eq. 165 .or. i .eq. 168 .or.
!    *             i .eq. 171 .or. i .eq. 174 .or. i .eq. 177 .or.
!    *             i .eq. 180 .or. i .eq. 183 .or. i .eq. 186     )then
            !
            u_rk1_ptr     = gfunc_pointer(i,   gi)
            u_rk3_ptr     = gfunc_pointer(i+2, gi)
            !
            ! Need to restrict auxiliary fields...those fields which show
            ! up on the lhs (i.e. in the L operator) but are not solved for
            ! with the elliptic solve:
            !     restrict _rk1 (h) fields -->_rk3 (2h)
            !    (these will then get transferred to the next coarser level)
            !
            call restrict(  q(u_rk3_ptr),  q(u_rk1_ptr),  nxc,nyc,nzc)
            !
            call maskLargeh2h(q(u_rk3_ptr),q(gr_chr),nxc,nyc,nzc,
     *                                                      bound_width)
            !
            if (vtrace) then
             call field_out3d(q(u_rk1_ptr),time,'CGC_'//gfunc_name(i),
     *                   mix,max,miy,may,miz,maz,nx, ny, nz, myid)
             call field_out3d(q(u_rk3_ptr),time,'CGC_'//gfunc_name(i+2),
     *                   mix,max,miy,may,miz,maz,nxc,nyc,nzc,myid)
            end if
            !
         end if
      end do


      if (ltrace) then
         write(*,99)myid,'grid_cgccompute: Finished.'
 99      format('[',I4,'] ',A,3I7)
      end if

      return
      end      ! END: grid_cgcompute

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_initMG:                                                              cc
cc                                                                            cc
cc             Initialize a grid from its parent upon creation.               cc
cc             This "MG" version just works on elliptic variables.            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_initMG( child, parent, r)
      implicit    none
      integer     child, parent, r
      include    'glob.inc'
      include    'grid.inc'
      include    'mask.inc'
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i
      integer     nxt, nyt, nzt, length
      integer     bi, bj, bk
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     parent_mask_ptr
      integer     proc_return_myid
      !
      ! For FV we need our own variables since we do not
      ! want to use an expanded region even if user calls for cubic interp
      !
      integer     min_icFV,max_icFV, min_jcFV,max_jcFV,min_kcFV,max_kcFV
      integer     min_ipFV,max_ipFV, min_jpFV,max_jpFV,min_kpFV,max_kpFV
      integer     nxtFV, nytFV, nztFV, lengthFV

      logical     ltrace
      parameter ( ltrace = .false. )


      call load_pointers(parent)
      parent_mask_ptr = gr_mask

      call load_pointers(child)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)

      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
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
      nxtFV    = max_ipFV - min_ipFV + 1
      nytFV    = max_jpFV - min_jpFV + 1
      nztFV    = max_kpFV - min_kpFV + 1
      !
      lengthFV = nxtFV * nytFV * nztFV

      if (linearbounds .eq. 0) then
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
         if (min_ip .gt.  1 ) then
            min_ip = min_ip - 1
            bi     = 1
         else
            min_ic = min_ic + refine_factor
            bi     = 1
         end if
         if (min_jp .gt.  1 ) then
            min_jp = min_jp - 1
            bj     = 1
         else
            min_jc = min_jc + refine_factor
            bj     = 1
         end if
         if (min_kp .gt.  1 ) then
            min_kp = min_kp - 1
            bk     = 1
         else
            min_kc = min_kc + refine_factor
            bk     = 1
         end if
         !
         if (max_ip .lt. nxp) then
            max_ip = max_ip + 1
         else
            max_ic = max_ic - refine_factor
         end if
         if (max_jp .lt. nyp)then
            max_jp = max_jp + 1
         else
            max_jc = max_jc - refine_factor
         end if
         if (max_kp .lt. nzp)then
            max_kp = max_kp + 1
         else
            max_kc = max_kc - refine_factor
         end if
      else
         bi = 0
         bj = 0
         bk = 0
      end if

      if (max_ic .lt. min_ic) return
      if (max_jc .lt. min_jc) return
      if (max_kc .lt. min_kc) return

      !
      ! Dimensions of part to be copied:
      !
      nxt    = max_ip - min_ip + 1
      nyt    = max_jp - min_jp + 1
      nzt    = max_kp - min_kp + 1

      length = nxt * nyt * nzt

      if (length.gt.nxc*nyc*nzc) then
         write(*,*) 'Problem of space in grid_init'
         call my_exit('Problem of space in grid_init')
      end if

      if (ltrace) then
         write(*,*) 'grid_init: linearbounds:  ',linearbounds
         write(*,*) 'grid_init: child, parent: ',child,parent
         write(*,*) 'grid_init: bi/j/k:        ',bi,bj,bk
         write(*,*) 'grid_init: nx/y/zp:       ',nxp,nyp,nzp
         write(*,*) 'grid_init: nx/y/zc:       ',nxc,nyc,nzc
         write(*,*) 'grid_init: nx/y/zt:       ',nxt,nyt,nzt
         write(*,*) 'grid_init: '
      end if

      !
      ! Initialize mask from parent first
      ! since we need it below:
      !
      if (num_masks .gt. 0) then
         ! First need to transfer the mask over:
         call copy_overlap( q(gr_tmp),
     *                       q(parent_mask_ptr),
     *                       min_ip, min_jp, min_kp,
     *                       nxt,nyt,nzt,
     *                       nxp,nyp,nzp )
         call paste_mask( q(gr_tmp),
     *                       q(gr_mask),
     *                       nxt,nyt,nzt,
     *                       min_ic,min_jc,min_kc,
     *                       max_ic,max_jc,max_kc,
     *                       bi,    bj,    bk,
     *                       nxc,nyc,nzc,refine_factor)
         if (ltrace) 
     *   call field_out3d(q(gr_mask),gr_t(child),'msk',
     *       gr_minx(child),gr_maxx(child),
     *       gr_miny(child),gr_maxy(child),
     *       gr_minz(child),gr_maxz(child),
     *             nxc,nyc,nzc,proc_return_myid())
      end if

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            ! Copy overlap area (perhaps a bit more):
            !
            call copy_overlapb(
     *                       q(gr_tmp),
     *                       q(gfunc_pointer(i,parent)),
     *                       q(parent_mask_ptr),
     *                       min_ip, min_jp, min_kp,
     *                       nxt,nyt,nzt,
     *                       nxp,nyp,nzp )
            if (linearbounds.gt.0) then
               !
               ! Linear interpolation:
               !
               if ( weno_interp .eq. 0 ) then
                 call init_fieldGWa( q(gr_tmp),
     *                            q(gfunc_pointer(i,child)),
     *                            nxt,nyt,nzt,
     *                            min_ic,min_jc,min_kc,
     *                            max_ic,max_jc,max_kc,
     *                            nxc,nyc,nzc, refine_factor,ghostwidth)
               else
                 call init_fieldGWa_weno( q(gr_tmp),
     *                            q(gfunc_pointer(i,child)),
     *                            nxt,nyt,nzt,
     *                            min_ic,min_jc,min_kc,
     *                            max_ic,max_jc,max_kc,
     *                            nxc,nyc,nzc, refine_factor,ghostwidth)
               end if               
            else
               !
               ! Cubic interpolation:
               !
               call interp_from_parentB( q(gr_tmp),
     *                          q(gfunc_pointer(i,child)),
     *                          q(gr_chr), q(gr_mask),
     *                          nxt,nyt,nzt,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            end if
         end if
      end do

      return
      end    ! END: grid_init


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_init_parentMG:                                                       cc
cc                         Copy grid function data within the region          cc
cc                     covered by child to tmp space and send to child.       cc
cc                     Input "grid" is the child grid to be initialized.      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_init_parentMG(grid,parent)
      implicit    none
      integer     grid, parent
      include     'grid.inc'
      include     'glob.inc'
      include     'mask.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     nxp, nyp, nzp, nxc, nyc, nzc, nxt,nyt,nzt,
     *            ax,  ay,  az,  i, length, owner,
     *            levp, levc
      integer     bi, bj, bk
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      !
      ! For FV we need our own variables since we do not
      ! want to use an expanded region even if user calls for cubic interp
      !
      integer     min_icFV,max_icFV, min_jcFV,max_jcFV,min_kcFV,max_kcFV
      integer     min_ipFV,max_ipFV, min_jpFV,max_jpFV,min_kpFV,max_kpFV
      integer     nxtFV, nytFV, nztFV, lengthFV
      logical     ltrace
      parameter ( ltrace = .false. )


      owner  = gr_own(grid)

      call load_pointers(parent)

      nxp    = gr_nx(parent)
      nyp    = gr_ny(parent)
      nzp    = gr_nz(parent)

      nxc    = gr_nx(grid)
      nyc    = gr_ny(grid)
      nzc    = gr_nz(grid)

      call grid_find_intersection(grid, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, grid, length,
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
      nxtFV    = max_ipFV - min_ipFV + 1
      nytFV    = max_jpFV - min_jpFV + 1
      nztFV    = max_kpFV - min_kpFV + 1
      !
      lengthFV = nxtFV * nytFV * nztFV

      if (linearbounds .eq. 0) then
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
         if (min_ip .gt.  1 ) then
            min_ip = min_ip - 1
            bi     = 1
         else
            min_ic = min_ic + refine_factor
            bi     = 1
         end if
         if (min_jp .gt.  1 ) then
            min_jp = min_jp - 1
            bj     = 1
         else
            min_jc = min_jc + refine_factor
            bj     = 1
         end if
         if (min_kp .gt.  1 ) then
            min_kp = min_kp - 1
            bk     = 1
         else
            min_kc = min_kc + refine_factor
            bk     = 1
         end if
         !
         if (max_ip .lt. nxp) then
            max_ip = max_ip + 1
         else
            max_ic = max_ic - refine_factor
         end if
         if (max_jp .lt. nyp)then
            max_jp = max_jp + 1
         else
            max_jc = max_jc - refine_factor
         end if
         if (max_kp .lt. nzp)then
            max_kp = max_kp + 1
         else
            max_kc = max_kc - refine_factor
         end if
      else
         bi = 0
         bj = 0
         bk = 0
      end if

      if (max_ic .lt. min_ic) return
      if (max_jc .lt. min_jc) return
      if (max_kc .lt. min_kc) return

      !
      ! Dimensions of part to be copied:
      !
      nxt    = max_ip - min_ip + 1
      nyt    = max_jp - min_jp + 1
      nzt    = max_kp - min_kp + 1

      length = nxt * nyt * nzt

      if (length.gt.nxp*nyp*nzp) then
         write(*,*) 'Problem of space in grid_init_parent'
         call my_exit('Problem of space in grid_init_parent')
      end if

      if (ltrace) then
         write(*,*) 'grid_init_parent:  grid   = ', grid
         write(*,*) 'grid_init_parent:  owner  = ', owner
         write(*,*) 'grid_init_parent: parent  = ', parent
         write(*,*) 'grid_init_parent: powner  = ', gr_own(parent)
         write(*,*) 'grid_init_parent: length  = ', length
         write(*,*) 'grid_init_parent: linearbounds:  ',linearbounds
         write(*,*) 'grid_init_parent: bi/j/k:        ',bi,bj,bk
         write(*,*) 'grid_init_parent: nx/y/zp:       ',nxp,nyp,nzp
         write(*,*) 'grid_init_parent: nx/y/zc:       ',nxc,nyc,nzc
         write(*,*) 'grid_init_parent: nx/y/zt:       ',nxt,nyt,nzt
         write(*,*) 'grid_init_parent: '
         !call grid_dump_info(grid)
         !call grid_dump_info(parent)
      end if

      ! Initialize mask from parent first
      ! since we need it below:
      !
      if (num_masks .gt. 0) then
         ! First need to transfer the mask over:
         call copy_overlap( q(gr_tmp),
     *                       q(gr_mask),
     *                       min_ip, min_jp, min_kp,
     *                       nxt,nyt,nzt,
     *                       nxp,nyp,nzp )
          call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INIT_DATA, MPI_COMM_WORLD, ierr)
      end if

      !
      ! Copies points of the parent that overlap with child to
      ! temporary storage
      !
      !
      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL) then
            call copy_overlapb(
     *                       q(gr_tmp),
     *                       q(gfunc_pointer(i,parent)),
     *                       q(gr_mask),
     *                       min_ip, min_jp, min_kp,
     *                       nxt,nyt,nzt,
     *                       nxp,nyp,nzp )
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INIT_DATA, MPI_COMM_WORLD, ierr)
         end if
      end do
 
      if (ltrace) then
         write(*,*) 'grid_init_parent: Done.'
      end if

      return
      end    ! END: grid_init_parentMG

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_init_childMG:                                                        cc
cc                         Copy grid funcation data received from parent      cc
cc                     to grid and interpolate information in order to        cc
cc                     initialize a grid.                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_init_childMG(grid, parent)
      implicit    none
      integer     grid, parent
      include     'grid.inc'
      include     'glob.inc'
      include     'mask.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     nxc, nyc, nzc, nxt,nyt,nzt, nxp,nyp,nzp,
     *            i, length, powner,
     *            levp, levc, countp, countc, grid_return_level
      integer     bi, bj, bk
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      logical     grid_return_existence
      !
      ! For FV we need our own variables since we do not
      ! want to use an expanded region even if user calls for cubic interp
      !
      integer     min_icFV,max_icFV, min_jcFV,max_jcFV,min_kcFV,max_kcFV
      integer     min_ipFV,max_ipFV, min_jpFV,max_jpFV,min_kpFV,max_kpFV
      integer     nxtFV, nytFV, nztFV, lengthFV
      logical     ltrace
      parameter ( ltrace = .false. )


      powner = gr_own(parent)

      call load_pointers(grid)

      nxc    = gr_nx(grid)
      nyc    = gr_ny(grid)
      nzc    = gr_nz(grid)

      nxp    = gr_nx(parent)
      nyp    = gr_ny(parent)
      nzp    = gr_nz(parent)

      call grid_find_intersection(grid, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, grid, length,
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
      nxtFV    = max_ipFV - min_ipFV + 1
      nytFV    = max_jpFV - min_jpFV + 1
      nztFV    = max_kpFV - min_kpFV + 1
      !
      lengthFV = nxtFV * nytFV * nztFV

      if (linearbounds .eq. 0) then
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
         if (min_ip .gt.  1 ) then
            min_ip = min_ip - 1
            bi     = 1
         else
            min_ic = min_ic + refine_factor
            bi     = 1
         end if
         if (min_jp .gt.  1 ) then
            min_jp = min_jp - 1
            bj     = 1
         else
            min_jc = min_jc + refine_factor
            bj     = 1
         end if
         if (min_kp .gt.  1 ) then
            min_kp = min_kp - 1
            bk     = 1
         else
            min_kc = min_kc + refine_factor
            bk     = 1
         end if
         !
         if (max_ip .lt. nxp) then
            max_ip = max_ip + 1
         else
            max_ic = max_ic - refine_factor
         end if
         if (max_jp .lt. nyp)then
            max_jp = max_jp + 1
         else
            max_jc = max_jc - refine_factor
         end if
         if (max_kp .lt. nzp)then
            max_kp = max_kp + 1
         else
            max_kc = max_kc - refine_factor
         end if
      else
         bi = 0
         bj = 0
         bk = 0
      end if

      if (max_ic .lt. min_ic) return
      if (max_jc .lt. min_jc) return
      if (max_kc .lt. min_kc) return

      !
      ! Dimensions of part to be copied:
      !
      nxt    = max_ip - min_ip + 1
      nyt    = max_jp - min_jp + 1
      nzt    = max_kp - min_kp + 1

      length =  nxt * nyt * nzt

      if (ltrace) then
         write(*,*) 'grid_init_child:  grid   = ', grid
         write(*,*) 'grid_init_child: parent  = ', parent
         write(*,*) 'grid_init_child: powner  = ', powner
         write(*,*) 'grid_init_child: length  = ', length
         write(*,*) 'grid_init_child: linearbounds:  ',linearbounds
         write(*,*) 'grid_init_child: bi/j/k:        ',bi,bj,bk
         write(*,*) 'grid_init_child: nx/y/zp:       ',nxp,nyp,nzp
         write(*,*) 'grid_init_child: nx/y/zc:       ',nxc,nyc,nzc
         write(*,*) 'grid_init_child: nx/y/zt:       ',nxt,nyt,nzt
         write(*,*) 'grid_init_child: '
         !call grid_dump_info(grid)
         !call grid_dump_info(parent)
      end if

      !
      ! Initialize mask from parent first
      ! since we need it below:
      !
      if (num_masks .gt. 0) then
         call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *               powner, TAG_INIT_DATA, MPI_COMM_WORLD, status,ierr)
         call paste_mask( q(gr_tmp),
     *                       q(gr_mask),
     *                       nxt,nyt,nzt,
     *                       min_ic,min_jc,min_kc,
     *                       max_ic,max_jc,max_kc,
     *                       bi,    bj,    bk,
     *                       nxc,nyc,nzc,refine_factor)
         if (ltrace)
     *   call field_out3d(q(gr_mask),gr_t(grid),'msk',
     *       gr_minx(grid),gr_maxx(grid),
     *       gr_miny(grid),gr_maxy(grid),
     *       gr_minz(grid),gr_maxz(grid),
     *             nxc,nyc,nzc,myid)
      end if

      !
      ! Copies points of the parent that overlap with child to
      ! temporary storage
      !
      !
      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL ) then
            !
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *               powner, TAG_INIT_DATA, MPI_COMM_WORLD, status,ierr)
            !
            if (linearbounds.gt.0) then
               !
               ! Linear interpolation:
               !
               call init_fieldGWa( q(gr_tmp),
     *                          q(gfunc_pointer(i,grid)),
     *                          nxt,nyt,nzt,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            else
               !
               ! Cubic interpolation:
               !
               call interp_from_parentb( q(gr_tmp),
     *                          q(gfunc_pointer(i,grid)),
     *                          q(gr_chr), q(gr_mask),
     *                          nxt,nyt,nzt,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            end if
         end if
      end do
 
      if (ltrace) then
         write(*,*) 'grid_init_child: Done.'
      end if

      return
      end    ! END: grid_init_childMG

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_mginject_master:                                                     cc
cc                 For MPI version for remotely injecting fields.             cc
cc                 Receives correction for processors which own children      cc
cc                 grids and applies those corrections to the respective      cc
cc                 fields. Complements grid_inject_slave.                     cc
cc             NB: Updated for overlapping grid using CHR mask for averaging. cc
cc         NB2:    These mginject routines differ from the normal injection   cc
cc                 in two ways: (i) they only inject the elliptic fields      cc
cc                 and (ii) they do not use ghostwidth, injecting all but     cc
cc                 the final point (so gw==1 effectively)                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_mginject_master(child,parent)
      implicit    none
      integer     child, parent
      include    'grid.inc'
      include    'glob.inc'
      include    'mpif.h'
      include    'mpi_stuff.inc'
      integer     length, owner
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i, nxs,nys,nzs
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     tmp_ghostwidth
      logical     double_equal
      external    double_equal
      character(3) tmps,tmps2
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )


      !call grid_check_mask(parent,20)

      owner  = gr_own(child)

      !
      ! This is the grid this process owns
      !
      call load_pointers(parent)

      !
      ! Parent's dimensions:
      !
      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)

      !
      ! Child's dimensions:
      !
      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      !
      ! After an elliptic solve initially, we inject but want
      ! to inject with ghostwidth==1 no matter what the real
      ! value is:
      !
      tmp_ghostwidth = 1

      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      !
      ! Overlap region on the parent
      !    but with resolution of the parent grid:
      !
      nxs = max_ip - min_ip + 1
      nys = max_jp - min_jp + 1
      nzs = max_kp - min_kp + 1

      !
      ! (linear) Size of grid to be transmitted:
      !
      length = nxs*nys*nzs

      if (ltrace) then
        call int2str(child,tmps)
        call int2str(parent,tmps2)
         write(*,99) myid,'grid_mginject_master: child   = ', child
         write(*,99) myid,'grid_mginject_master: parent  = ', parent
         write(*,99) myid,'grid_mginject_master: length  = ', length
         write(*,99) myid,'grid_mginject_master: owner   = ', owner
         write(*,99) myid,'grid_mginject_master: nx/y/zp = ',nxp,nyp,nzp
         write(*,99) myid,'grid_mginject_master: nx/y/zp = ',nxs,nys,nzs
         write(*,99) myid,'grid_mginject_master: nx/y/zc = ',nxc,nyc,nzc
         write(*,99) myid,'grid_mginject_master: mini/j/kp=',
     *                                            min_ip,min_jp,min_kp
         write(*,99) myid,'grid_mginject_master: mini/j/kc=',
     *                                            min_ic,min_jc,min_kc
         write(*,99) myid,'grid_mginject_master: maxi/j/kc=',
     *                                            max_ic,max_jc,max_kc
      end if

      do i = 1, num_gfuncs
         if (gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL)then
            !
            ! Receive correction for grid owner:
            !
            if (ltrace2) then
               write(*,99) myid,' master: receiving correction'
               write(*,99) myid,' master:     field: ',i
               write(*,99) myid,' master: from owner: ', owner
            end if
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INJ_DATA, MPI_COMM_WORLD, status,ierr)
            if (i.eq.1 .and. ltrace)
     *      call field_out3d(q(gr_tmp),gr_t(parent),
     *           'InjectMTmp'//tmps//tmps2,
     *       q(gr_x(parent)+(min_ip-1)),q(gr_x(parent)+(max_ip-1)),
     *       q(gr_y(parent)+(min_jp-1)),q(gr_y(parent)+(max_jp-1)),
     *       q(gr_z(parent)+(min_kp-1)),q(gr_z(parent)+(max_kp-1)),
     *             nxs,nys,nzs,myid)
            !
            ! Replace values from tmp storage:
            !
            call replace_part_field_avgGWc( 
     *                           q(gfunc_pointer(i,parent)),
     *                           q(gr_chr),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth,
     *                           bound_width )
            if (ltrace2)
     *      call field_dump_info(q(gfunc_pointer(i,parent)),nxs,nys,nzs)
            if (i.eq.1 .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,parent)),gr_t(parent),
     *                   'InjectMPar'//tmps//tmps2,
     *               gr_minx(parent),gr_maxx(parent),
     *               gr_miny(parent),gr_maxy(parent),
     *               gr_minz(parent),gr_maxz(parent),
     *               nxp,nyp,nzp,myid)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *             gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INJ_DATA, MPI_COMM_WORLD, status,ierr)
            !call load_scal1d(q(gr_tmp),-1.d0*child,nxs*nys*nzs)
            call replace_part_field_avgFV(
     *                           q(gfunc_pointer(i,parent)),
     *                           q(gr_chr),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth,
     *                           bound_width )
         end if
      end do

      !
      ! this isn't field specific, but we can use the last field
      ! updated above to get info about which points are being
      ! restricted:
      !
      call inc_mask_overlapB(    q(gr_chr),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth )

      ! Enforce symmetry conditions after injection:
      if (assume_symmetry.ne.0) call grid_symm(parent)

      if (ltrace)
     *      call field_out3d(q(gfunc_pointer(1,parent)),gr_t(parent),
     *                   'InjectMPOst'//tmps//tmps2,
     *               gr_minx(parent),gr_maxx(parent),
     *               gr_miny(parent),gr_maxy(parent),
     *               gr_minz(parent),gr_maxz(parent),
     *               nxp,nyp,nzp,myid)

      !call grid_check_mask(parent,21)

      if (ltrace) then
         call grid_dump_info(parent)
         !call grid_dump_info(child)
         write(*,99)myid,'grid_mginject_master: Done.',parent, child
      end if
 99   format('[',I3,'] ',A,3I7)

      return
      end    ! END: grid_mginject_master

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_mginject_slave:                                                      cc
cc                 For MPI version for remotely injecting fields.             cc
cc                 Restricts fields into tmp storage, and then passes that    cc
cc                 to the parent of this grid so that "grid_inject_master"    cc
cc                 does the rest of the work in terms of updating the fields. cc
cc             NB: Apparently no change to allow for overlapping grids.       cc
cc         NB2:    These mginject routines differ from the normal injection   cc
cc                 in two ways: (i) they only inject the elliptic fields      cc
cc                 and (ii) they do not use ghostwidth, injecting all but     cc
cc                 the final point (so gw==1 effectively)                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_mginject_slave(child,parent)
      implicit    none
      integer     child, parent
      include    'grid.inc'
      include    'glob.inc'
      include    'mpif.h'
      include    'mpi_stuff.inc'
      integer     length, owner
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i, nxs,nys,nzs
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     tmp_ghostwidth
      character(3) tmps,tmps2
      logical     double_equal
      external    double_equal
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      owner  = gr_own(parent)

      call load_pointers(child)

      !
      ! After an elliptic solve initially, we inject but want
      ! to inject with ghostwidth==1 no matter what the real
      ! value is:
      !
      tmp_ghostwidth = 1

      !
      ! Parent's dimensions:
      !
      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)

      !
      ! Child's dimensions:
      !
      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      !
      ! Overlap region on the parent
      !    but with resolution of the parent grid:
      !
      nxs = max_ip - min_ip + 1
      nys = max_jp - min_jp + 1
      nzs = max_kp - min_kp + 1

      !
      ! (linear) Size of grid to be transmitted:
      !
      length = nxs*nys*nzs

      if (ltrace) then
        call int2str(child,tmps)
        call int2str(parent,tmps2)
         write(*,99) myid,'grid_mginject_slave: child   = ', child
         write(*,99) myid,'grid_mginject_slave: parent  = ', parent
         write(*,99) myid,'grid_mginject_slave: length  = ', length
         write(*,99) myid,'grid_mginject_slave: owner   = ', owner
         write(*,99) myid,'grid_mginject_slave: nx/y/zp = ',nxp,nyp,nzp
         write(*,99) myid,'grid_mginject_slave: nx/y/zs = ',nxs,nys,nzs
         write(*,99) myid,'grid_mginject_slave: nx/y/zc = ',nxc,nyc,nzc
         write(*,99) myid,'grid_mginject_slave: mini/j/kp=',
     *                                            min_ip,min_jp,min_kp
         write(*,99) myid,'grid_mginject_slave: mini/j/kc=',
     *                                            min_ic,min_jc,min_kc
         write(*,99) myid,'grid_mginject_slave: maxi/j/kc=',
     *                                            max_ic,max_jc,max_kc
      end if

      do i = 1, num_gfuncs
         if (gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL)then
            !
            ! Restrict to tmp storage:
            !
            if (ltrace2)call field_dump_info(q(gfunc_pointer(i,child)),
     *                                         nxc,nyc,nzc)
            if (i.eq.1 .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,child)),gr_t(child),
     *                   'InjectSChi'//tmps//tmps2,
     *               gr_minx(child),gr_maxx(child),
     *               gr_miny(child),gr_maxy(child),
     *               gr_minz(child),gr_maxz(child),
     *               nxc,nyc,nzc,myid)
            if (i.eq.1 .and. ltrace)
     *      call field_out3d(q(gr_chr), gr_t(child),
     *                   'InjectSChr'//tmps//tmps2,
     *               gr_minx(child),gr_maxx(child),
     *               gr_miny(child),gr_maxy(child),
     *               gr_minz(child),gr_maxz(child),
     *               nxc,nyc,nzc,myid)
            call restrict_fieldGWb(q(gr_tmp),
     *                           q(gfunc_pointer(i,child)),
     *                           q(gr_chr),
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor,tmp_ghostwidth)
            if (ltrace2)call field_dump_info(q(gr_tmp),nxs,nys,nzs)
            if (i.eq.1 .and. ltrace)
     *      call field_out3d(q(gr_tmp),gr_t(child),
     *           'InjectSTmp'//tmps//tmps2,
     *       q(gr_x(child)+(min_ic-1)),q(gr_x(child)+(max_ic-1)),
     *       q(gr_y(child)+(min_jc-1)),q(gr_y(child)+(max_jc-1)),
     *       q(gr_z(child)+(min_kc-1)),q(gr_z(child)+(max_kc-1)),
     *             nxs,nys,nzs,myid)
            !
            ! Send to owner of parent:
            !
            if (ltrace2) then
               write(*,99) myid,' slave:  sending correction'
               write(*,99) myid,' slave:      field: ',i
               write(*,99) myid,' slave:   to owner: ', owner
            end if
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INJ_DATA, MPI_COMM_WORLD, ierr)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *             gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
            call restrict_fieldFV(q(gr_tmp),
     *                           q(gfunc_pointer(i,child)),
     *                           q(gr_chr),
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor,tmp_ghostwidth)
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INJ_DATA, MPI_COMM_WORLD, ierr)
         end if
      end do

      if (ltrace) then
         !call grid_dump_info(parent)
         call grid_dump_info(child)
         write(*,99)myid,'grid_mginject_slave: Done.',parent, child
      end if
 99   format('[',I3,'] ',A,3I7)

      return
      end    ! END: grid_mginject_slave

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_mginject:                                                            cc
cc                 Update parent grid with values from                        cc
cc                 fine grids. Only injects evolved quantities.               cc
cc                 Allow   derived fields to be computed on parent grid.      cc
cc         NB2:    These mginject routines differ from the normal injection   cc
cc                 in two ways: (i) they only inject the elliptic fields      cc
cc                 and (ii) they do not use ghostwidth, injecting all but     cc
cc                 the final point (so gw==1 effectively)                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_mginject(child,parent)
      implicit    none
      integer     child, parent
      include     'grid.inc'
      include     'glob.inc'
      integer     nxs,nys,nzs, parent_mask_pointer
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     tmp_ghostwidth
      logical     double_equal
      external    double_equal
      integer     length
      integer     myid, proc_return_myid
      character(3) tmps,tmps2
      logical     ltrace
      parameter ( ltrace = .false. )

      !call grid_check_mask(parent,10)

      myid = proc_return_myid()
      !
      ! Store the address of the parent's mask function:
      !
      call load_pointers(parent)
      parent_mask_pointer = gr_chr

      call load_pointers(child)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)

      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      !
      ! After an elliptic solve initially, we inject but want
      ! to inject with ghostwidth==1 no matter what the real
      ! value is:
      !
      tmp_ghostwidth = 1

      !
      ! Find indices which bound intersection region:
      !
      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      !
      ! Overlap region on the parent
      !    but with resolution of the parent grid:
      !
      nxs = max_ip - min_ip + 1
      nys = max_jp - min_jp + 1
      nzs = max_kp - min_kp + 1

      if (nxs*nys*nzs .gt. nxc*nyc*nzc) then
         write(*,*) 'grid_mginject: Problem. Not enough storage'
         write(*,*) 'grid_mginject: nx/y/zs   = ',nxs,nys,nzs
         write(*,*) 'grid_mginject: nx/y/zc   = ',nxc,nyc,nzc
         write(*,*) 'grid_mginject: nx*y*zs   = ',nxs*nys*nzs
         write(*,*) 'grid_mginject: nx*y*zc   = ',nxc*nyc*nzc
         call my_exit('Not enough storage in grid_mginject')
      end if

      if (ltrace) then
        call int2str(child,tmps)
        call int2str(parent,tmps2)
         write(*,99)myid,'grid_mginject: parent    = ', parent
         write(*,99)myid,'grid_mginject: child     = ', child
         write(*,99)myid,'grid_mginject: length    = ', length
         write(*,99)myid,'grid_mginject: nx/y/zs   = ',nxs,nys,nzs
         write(*,99)myid,'grid_mginject: nx/y/zc   = ',nxc,nyc,nzc
         write(*,99)myid,'grid_mginject: nx/y/zp   = ',nxp,nyp,nzp
         write(*,99)myid,'grid_mginject: min/max_ic= ',min_ic,max_ic
         write(*,99)myid,'grid_mginject: min/max_jc= ',min_jc,max_jc
         write(*,99)myid,'grid_mginject: min/max_kc= ',min_kc,max_kc
         write(*,99)myid,'grid_mginject: min/max_ip= ',min_ip,max_ip
         write(*,99)myid,'grid_mginject: min/max_jp= ',min_jp,max_jp
         write(*,99)myid,'grid_mginject: min/max_kp= ',min_kp,max_kp
         write(*,98)myid,'grid_mginject: time(parent): ',gr_t(parent)
         write(*,98)myid,'grid_mginject: time(child) : ',gr_t(child)
         call field_dump_info(q(gfunc_pointer(1,parent)),nxp,nyp,nzp)
         call field_dump_info(q(gfunc_pointer(1,child)),nxc,nyc,nzc)
      end if

      do i = 1, num_gfuncs
         if (gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL)then
            if (i.eq.1 .and. ltrace)
     *      call field_out3d(q(gr_chr), gr_t(child),
     *                   'InjectChr'//tmps//tmps2,
     *               gr_minx(child),gr_maxx(child),
     *               gr_miny(child),gr_maxy(child),
     *               gr_minz(child),gr_maxz(child),
     *               nxc,nyc,nzc,myid)
            if (i.eq.1 .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,child)),gr_t(child),
     *                   'InjectChi'//tmps//tmps2,
     *               gr_minx(child),gr_maxx(child),
     *               gr_miny(child),gr_maxy(child),
     *               gr_minz(child),gr_maxz(child),
     *               nxc,nyc,nzc,myid)
            call restrict_fieldGWb(q(gr_tmp),
     *                           q(gfunc_pointer(i,child)),
     *                           q(gr_chr),
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor,tmp_ghostwidth)
            !call load_scal1d(q(gr_tmp),-2.d0,nxs*nys*nzs)
            if (i.eq.1 .and. ltrace)
     *      call field_out3d(q(gr_tmp),gr_t(parent),
     *           'InjectTmp'//tmps//tmps2,
     *       q(gr_x(parent)+(min_ip-1)),q(gr_x(parent)+(max_ip-1)),
     *       q(gr_y(parent)+(min_jp-1)),q(gr_y(parent)+(max_jp-1)),
     *       q(gr_z(parent)+(min_kp-1)),q(gr_z(parent)+(max_kp-1)),
     *             nxs,nys,nzs,myid)
            call replace_part_field_avgGWc( 
     *                           q(gfunc_pointer(i,parent)),
     *                           q(parent_mask_pointer),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth,
     *                           bound_width )
            if (i.eq.1 .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,parent)),gr_t(parent),
     *                   'InjectPar'//tmps//tmps2,
     *               gr_minx(parent),gr_maxx(parent),
     *               gr_miny(parent),gr_maxy(parent),
     *               gr_minz(parent),gr_maxz(parent),
     *               nxp,nyp,nzp,myid)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *             gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
            call restrict_fieldFV(q(gr_tmp),
     *                           q(gfunc_pointer(i,child)),
     *                           q(gr_chr),
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor,tmp_ghostwidth)
            !call load_scal1d(q(gr_tmp),-1.d0*child,nxs*nys*nzs)
            call replace_part_field_avgFV(
     *                           q(gfunc_pointer(i,parent)),
     *                           q(parent_mask_pointer),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth,
     *                           bound_width )
         !
         end if
      end do

      ! Enforce symmetry conditions after injection:
      if (assume_symmetry.ne.0) call grid_symm(parent)

      if (ltrace)
     *      call field_out3d(q(gfunc_pointer(1,parent)),gr_t(parent),
     *                   'InjectPOst'//tmps//tmps2,
     *               gr_minx(parent),gr_maxx(parent),
     *               gr_miny(parent),gr_maxy(parent),
     *               gr_minz(parent),gr_maxz(parent),
     *               nxp,nyp,nzp,myid)

      call inc_mask_overlapB(    q(parent_mask_pointer),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth )

      if (ltrace)
     *     call field_out3d(   q(parent_mask_pointer),
     *             gr_t(parent)*child,'inject_mask',
     *             gr_minx(parent),gr_maxx(parent),
     *             gr_miny(parent),gr_maxy(parent),
     *             gr_minz(parent),gr_maxz(parent),
     *             nxp,nyp,nzp,myid)

      !call grid_check_mask(parent,11)

      if (ltrace) then
         call grid_dump_info(parent)
         call grid_dump_info(child)
         write(*,99)myid,'grid_mginject: Done.',parent, child
      end if

 97   format('[',I3,'] ',A,I5,3F10.5)
 98   format('[',I3,'] ',A,3F10.5)
 99   format('[',I3,'] ',A,3I5)

      return
      end    ! END: grid_mginject
