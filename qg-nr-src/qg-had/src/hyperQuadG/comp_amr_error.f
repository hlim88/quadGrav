cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  comp_amr_error:                                                           cc
cc                   Needs to set the field gr_error to have                  cc
cc                   the appropriate error function....                       cc
cc                   therefore needs to be non-negative function              cc
cc                   of the fields.                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comp_amr_error( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'grid_methods.inc'
      include    'glob.inc'
      include    'param.inc'
      integer     nx, ny, nz, numfields
      real*8      hg,  tmp, tmp_maxchi
      real*8      distance_wanted, distance_actual
      real*8      mask_distance_movedx
      external    mask_distance_movedx
      integer     tmp_maxchi_i, tmp_maxchi_j,tmp_maxchi_k
      integer     level
      real*8      xmin, ymin, zmin, x, y, z
      logical     ltrace
      parameter ( ltrace = .false. )
      integer    i,j,k,index, l
      integer    nd, d_type
      integer    dx_alpha, dy_alpha, dz_alpha
      integer    dx_gt11, dy_gt11, dz_gt11
      integer    dx_gt12, dy_gt12, dz_gt12
      integer    dx_gt13, dy_gt13, dz_gt13
      integer    dx_gt22, dy_gt22, dz_gt22
      integer    dx_gt23, dy_gt23, dz_gt23
      integer    dx_gt33, dy_gt33, dz_gt33
      real*8     local_dx, local_dy, local_dz

      integer     mem_alloc
      external    mem_alloc

      call load_pointers(gridnum)

      if (ltrace) write(*,*) 'comp_amr_error: On grid: ',gridnum

      nx   = gr_nx(gridnum)
      ny   = gr_ny(gridnum)
      nz   = gr_nz(gridnum)
      hg   = gr_h(gridnum)

      level = grid_return_level(gridnum)

      xmin = q(gr_x(gridnum))
      ymin = q(gr_y(gridnum))
      zmin = q(gr_z(gridnum))

!     take derivatives of alpha for the refinement
      nd = nx * ny * nz

!      dx_alpha = mem_alloc(nd)
!      dy_alpha = mem_alloc(nd)
!      dz_alpha = mem_alloc(nd)

!      dx_gt11 = mem_alloc(nd)

      local_dx = hg
      local_dy = hg
      local_dz = hg
      d_type = deriv_order

!      call deriv_x(q(dx_alpha), q(gr_alpha), local_dx, 
!     *             nx, ny, nz, d_type)
!      call deriv_y(q(dy_alpha), q(gr_alpha), local_dy, 
!     *             nx, ny, nz, d_type)
!      call deriv_z(q(dz_alpha), q(gr_alpha), local_dz, 
!     *             nx, ny, nz, d_type)

!      call deriv_x(q(dx_gt11), q(gr_gt11), local_dx, 
!     *             nx, ny, nz, d_type)
!      call deriv_y(q(dy_gt11), q(gr_gt11), local_dy, 
!     *             nx, ny, nz, d_type)
!      call deriv_z(q(dz_gt11), q(gr_gt11), local_dz, 
!     *             nx, ny, nz, d_type)

      if (shadow .eq. 0) then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              x = xmin + hg*(i-1)
              y = ymin + hg*(j-1)
              z = zmin + hg*(k-1)
              index = (k-1)*ny*nx + (j-1)*nx + i

              q(gr_error + index - 1) = hg**2 * (
     *              + q(gr_phiR+index-1)**2/hg**2
     *              + q(gr_phiI+index-1)**2/hg**2
     *                )

            end do
          end do
        end do
      end if

      ! deallocate storage
!      call mem_dealloc(dx_alpha, nd)
!      call mem_dealloc(dy_alpha, nd)
!      call mem_dealloc(dz_alpha, nd)

!      call mem_dealloc(dx_gt11, nd)
!      call mem_dealloc(dy_gt11, nd)
!      call mem_dealloc(dz_gt11, nd)

      return
      end    ! END: comp_amr_error
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_ahfind:                                                              cc
cc                                                                            cc
cc         Project specific routine to find where to excise and hence         cc
cc         should be located in comp_amr_error.f.                             cc
cc         This routines requires centers for up to masks to be specified.    cc
cc         All this is subject to change, but for now this routine is         cc
cc         responsible for determining on this grid whether the quantities    cc
cc             local_left(1..2,1..3)                                          cc
cc             local_right(1..2,1..3)                                         cc
cc         need to be adjusted. These are local to this processor, and        cc
cc         global bounds will be found in level_bhmask_local().               cc
cc                                                                            cc
cc         To see the data structures for this routine, look in:              cc
cc                  had/include/mask.inc                                      cc
cc                                                                            cc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_ahfind(gi)
      implicit none
      integer  gi
      include 'mpif.h'
      include 'mpi_stuff.inc'
      include 'grid.inc'
      include 'mask.inc'
      include 'param.inc'      
      integer  nx,ny,nz, i,j,k, kk
      real(kind=8) my_left(3,max_num_masks),
     *             my_right(3,max_num_masks),
     *             bh_mass, time
      logical     bh_found
      character*3  tmpstr
      character*12 myname

      logical     ltrace
      parameter ( ltrace  = .false.)
      logical     ltrace2
      parameter ( ltrace2 = .false.)
      logical     ltrace3
      parameter ( ltrace3 = .false.)
      logical     ltrace4
      parameter ( ltrace4 = .false.)

      call load_pointers(gi)

      nx = gr_nx(gi)
      ny = gr_ny(gi)
      nz = gr_nz(gi)

c      time = par(P_TIME)

      num_masks = 0

      !
      ! Do the real computation for this grid:
      !
      do kk = 1, max_num_masks
        !
        if (ltrace3) then
           write(*,99) myid, 'grid_ahfind: Center:',kk
           do i = 1, 3
              write(*,98) myid, '    center coord',mask_center(i,kk)
              write(*,98) myid, '    left/right: ',mask_left(i,kk),
     *                                            mask_right(i,kk)
           end do
           write(*,98) myid, '    radius:     ',mask_radius(1,kk)
        end if
        !
        bh_found    = .false.
 
        if (kk .EQ. 1) bh_mass = bh1_mass
        if (kk .EQ. 2) bh_mass = bh2_mass
        if(forcedmerge) bh_mass = bh1_mass + bh2_mass

        if (bh1_mass.lt.0.0) bh_mass = abs(bh1_mass)

        if (local_time .LT. 1E-10) then
        else
        end if
        
	if ( bh_found ) then
            bh_true(kk) = .true.
            num_masks = num_masks + 1
            if(ltrace)write(*,99)myid, 'grid_ahfind: BH found: ',kk
            do i = 1, 3
               if (ltrace2)write(*,98)myid, ' l/right:',
     *                              my_left(i,kk),my_right(i,kk)
               if (my_left(i,kk)  .lt. local_left(i,kk)  )
     *             local_left(i,kk)  = my_left(i,kk)
               if (my_right(i,kk) .gt. local_right(i,kk) )
     *             local_right(i,kk) = my_right(i,kk)
               if(ltrace)write(*,98)myid,'      local_left/right(i,k)=',
     *            local_left(i,kk), local_right(i,kk)
            end do
        end if
         !
         !
        if (ltrace4) then
        end if 
         !
      end do

      if (ltrace2) then
         write(*,99) myid, 'grid_ahfind: gi     = ',gi
         write(*,99) myid, 'grid_ahfind: nx/y/z = ',nx,ny,nz
         write(*,99) myid, 'grid_ahfind: Done.'
      end if

      if (local_time .GT. 1E-10) then
        if(forcedmerge) then
          bh1_x0 = 0.5*(mask_coords(1,1)+mask_coords(2,1))
          bh1_y0 = 0.5*(mask_coords(3,1)+mask_coords(4,1))
          bh1_z0 = 0.5*(mask_coords(5,1)+mask_coords(6,1))
          bh2_x0 = bh1_x0
          bh2_y0 = bh1_y0
          bh2_z0 = bh1_z0
        else
          bh1_x0 = 0.5*(mask_coords(1,1)+mask_coords(2,1))
          bh1_y0 = 0.5*(mask_coords(3,1)+mask_coords(4,1))
          bh1_z0 = 0.5*(mask_coords(5,1)+mask_coords(6,1))
          bh2_x0 = 0.5*(mask_coords(1,2)+mask_coords(2,2))
          bh2_y0 = 0.5*(mask_coords(3,2)+mask_coords(4,2))
          bh2_z0 = 0.5*(mask_coords(5,2)+mask_coords(6,2))
        end if
       
        bh1_exc_rad = mask_radius(1,1)
        bh2_exc_rad = mask_radius(1,2)
      end if



  98  format('[',I3,'] ',A,3F18.6)
  99  format('[',I3,'] ',A,3I5)
      return
      end       ! END: grid_ahfind

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  compute_boundary_bh:                                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_boundary_bh( )
      implicit   none

      return
      end subroutine compute_boundary_bh

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_asurface_comp:                                                        cc
cc                   Takes as input the coordinate of a point located         cc
cc                   on the grid. It needs to compute the various             cc
cc                   surface functions at that point.                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_asurface_comp(gi, x,y,z, itheta, jphi)
      implicit     none
      integer      gi, itheta, jphi
      real(kind=8) x,y,z
      include     'grid.inc'
      include     'surfaces.inc'
      integer      nx, ny, nz
      integer      i,j,k,index, l
      integer      ilow, jlow, klow
      real(kind=8) hg, rpsi4, ipsi4, fracx, fracy, fracz, R2scalar
      real(kind=8) rphi0, iphi0, rphi2, iphi2
      real(kind=8) interp, massADM, lin_momx, lin_momy, ang_momz
      real(kind=8) gthethe, gthephi, gphiphi, gur, tau
      external     interp
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      ! Index of lower left point of interpolation cell:
      ilow  = INT( (x - gr_minx(gi)) / hg ) + 1
      jlow  = INT( (y - gr_miny(gi)) / hg ) + 1
      klow  = INT( (z - gr_minz(gi)) / hg ) + 1

      ! Fractional distance point is from lower left point:
      fracx = ( x - q( gr_x(gi) + (ilow-1) ) ) / hg
      fracy = ( y - q( gr_y(gi) + (jlow-1) ) ) / hg
      fracz = ( z - q( gr_z(gi) + (klow-1) ) ) / hg

      if (ilow .lt. 1) then
         write(*,*) 'grid_surface_comp: x too small ',ilow
         call my_exit('grid_surface_comp: X too small')
      end if
      if (fracx.lt.-0.01 .or.fracy.lt.-0.01.or.fracz.lt.-0.0001) then
         write(*,*) 'grid_surface_comp: fraction cannot be negative'
         write(*,*) 'grid_surface_comp: fracx/y/z =',fracx,fracy,fracz
         call my_exit('grid_surface_comp: fraction cannot be negative')
      end if

      rphi0 = interp(  q(gr_phi0R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      iphi0  = interp(  q(gr_phi0I),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      rphi2 = interp(  q(gr_phi2R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      iphi2  = interp(  q(gr_phi2I),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      rpsi4 = interp(  q(gr_psi4R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      ipsi4  = interp(  q(gr_psi4I),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)

      index = (jphi-1)*asf_shape(1) + (itheta-1)
      q(asf_rap0+index) = rphi0
      q(asf_iap0+index) = iphi0
      q(asf_rap2+index) = rphi2
      q(asf_iap2+index) = iphi2
      q(asf_rap4+index) = rpsi4
      q(asf_iap4+index) = ipsi4


      return
      end    ! END: grid_surface_comp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_bsurface_comp:                                                        cc
cc                   Takes as input the coordinate of a point located         cc
cc                   on the grid. It needs to compute the various             cc
cc                   surface functions at that point.                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_bsurface_comp(gi, x,y,z, itheta, jphi)
      implicit     none
      integer      gi, itheta, jphi
      real(kind=8) x,y,z
      include     'grid.inc'
      include     'surfaces.inc'
      integer      nx, ny, nz
      integer      i,j,k,index, l
      integer      ilow, jlow, klow
      real(kind=8) hg, rpsi4, ipsi4, fracx, fracy, fracz, R2scalar
      real(kind=8) rphi0, iphi0, rphi2, iphi2
      real(kind=8) interp, massADM, lin_momx, lin_momy, ang_momz
      real(kind=8) gthethe, gthephi, gphiphi, gur, tau
      external     interp
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      ! Index of lower left point of interpolation cell:
      ilow  = INT( (x - gr_minx(gi)) / hg ) + 1
      jlow  = INT( (y - gr_miny(gi)) / hg ) + 1
      klow  = INT( (z - gr_minz(gi)) / hg ) + 1

      ! Fractional distance point is from lower left point:
      fracx = ( x - q( gr_x(gi) + (ilow-1) ) ) / hg
      fracy = ( y - q( gr_y(gi) + (jlow-1) ) ) / hg
      fracz = ( z - q( gr_z(gi) + (klow-1) ) ) / hg

      if (ilow .lt. 1) then
         write(*,*) 'grid_surface_comp: x too small ',ilow
         call my_exit('grid_surface_comp: X too small')
      end if
      if (fracx.lt.-0.01 .or.fracy.lt.-0.01.or.fracz.lt.-0.0001) then
         write(*,*) 'grid_surface_comp: fraction cannot be negative'
         write(*,*) 'grid_surface_comp: fracx/y/z =',fracx,fracy,fracz
         call my_exit('grid_surface_comp: fraction cannot be negative')
      end if

      rphi0 = interp(  q(gr_phi0R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      iphi0  = interp(  q(gr_phi0I),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      rphi2 = interp(  q(gr_phi2R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      iphi2  = interp(  q(gr_phi2I),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      rpsi4 = interp(  q(gr_psi4R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      ipsi4  = interp(  q(gr_psi4I),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)

      index = (jphi-1)*bsf_shape(1) + (itheta-1)
      q(bsf_rbp0+index) = rphi0
      q(bsf_ibp0+index) = iphi0
      q(bsf_rbp2+index) = rphi2
      q(bsf_ibp2+index) = iphi2
      q(bsf_rbp4+index) = rpsi4
      q(bsf_ibp4+index) = ipsi4
      return
      end    ! END: grid_surface_comp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_csurface_comp:                                                        cc
cc                   Takes as input the coordinate of a point located         cc
cc                   on the grid. It needs to compute the various             cc
cc                   surface functions at that point.                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_csurface_comp(gi, x,y,z, itheta, jphi)
      implicit     none
      integer      gi, itheta, jphi
      real(kind=8) x,y,z
      include     'grid.inc'
      include     'surfaces.inc'
      integer      nx, ny, nz
      integer      i,j,k,index, l
      integer      ilow, jlow, klow
      real(kind=8) hg, rpsi4, ipsi4, fracx, fracy, fracz, R2scalar
      real(kind=8) rphi0, iphi0, rphi2, iphi2
      real(kind=8) interp, massADM, lin_momx, lin_momy, ang_momz
      real(kind=8) gthethe, gthephi, gphiphi, gur, tau
      external     interp
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      ! Index of lower left point of interpolation cell:
      ilow  = INT( (x - gr_minx(gi)) / hg ) + 1
      jlow  = INT( (y - gr_miny(gi)) / hg ) + 1
      klow  = INT( (z - gr_minz(gi)) / hg ) + 1

      ! Fractional distance point is from lower left point:
      fracx = ( x - q( gr_x(gi) + (ilow-1) ) ) / hg
      fracy = ( y - q( gr_y(gi) + (jlow-1) ) ) / hg
      fracz = ( z - q( gr_z(gi) + (klow-1) ) ) / hg

      if (ilow .lt. 1) then
         write(*,*) 'grid_surface_comp: x too small ',ilow
         call my_exit('grid_surface_comp: X too small')
      end if
      if (fracx.lt.-0.01 .or.fracy.lt.-0.01.or.fracz.lt.-0.0001) then
         write(*,*) 'grid_surface_comp: fraction cannot be negative'
         write(*,*) 'grid_surface_comp: fracx/y/z =',fracx,fracy,fracz
         call my_exit('grid_surface_comp: fraction cannot be negative')
      end if

      rphi0 = interp(  q(gr_phi0R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow, 
     *               nx,ny,nz)
      iphi0  = interp(  q(gr_phi0I), 
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      rphi2 = interp(  q(gr_phi2R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow, 
     *               nx,ny,nz)
      iphi2  = interp(  q(gr_phi2I), 
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      rpsi4 = interp(  q(gr_psi4R),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow, 
     *               nx,ny,nz)
      ipsi4  = interp(  q(gr_psi4I), 
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)

      index = (jphi-1)*csf_shape(1) + (itheta-1)
      q(csf_rcp0+index) = rphi0
      q(csf_icp0+index) = iphi0
      q(csf_rcp2+index) = rphi2
      q(csf_icp2+index) = iphi2
      q(csf_rcp4+index) = rpsi4
      q(csf_icp4+index) = ipsi4


      return
      end    ! END: grid_surface_comp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_modify_data:                                                         cc
cc                   Here one need not do anything. However, this             cc
cc                   gets called after reading in a checkpoint file           cc
cc                   and also after global analysis quantities are called     cc
cc                   and so one can use those to set data.                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_modify_data( gridnum )
      implicit    none
      integer     gridnum
      include    'mpif.h'
      include    'mpi_stuff.inc'
      include    'grid.inc'
      include    'grid_methods.inc'
      include    'glob.inc'
      include    'param.inc'
      integer     nx, ny, nz, numfields
      real*8      hg,  tmp, tmp_maxchi
      real*8      distance_wanted, distance_actual
      real*8      mask_distance_movedx
      external    mask_distance_movedx
      integer     tmp_maxchi_i, tmp_maxchi_j,tmp_maxchi_k
      integer     level
      real*8      xmin, ymin, zmin, x, y, z, xc, yc, zc
      logical     ltrace
      parameter ( ltrace = .true. )
      integer    i,j,k,index, l
      integer    dx_alpha, dy_alpha, dz_alpha, nd, d_type
      real*8     local_dx, local_dy, local_dz
      real*8     xs1, ys1, zs1, xs2, ys2, zs2, dist, t1, t2
      real*8     myP, myrho
      real*8     dxAx,dyAx,dzAx
      real*8     dxAy,dyAy,dzAy
      real*8     dxAz,dyAz,dzAz
      real*8     sdetg_pt, r1sq, r2sq, invsdetg_pt
      real*8     myl2norm3D
      external   myl2norm3D
      !
      ! Factor by which rho must be greater than vacuum to have nonzero Bfield:
      !
      real*8     RHOTHRESH
      parameter (RHOTHRESH = 1d3)
      logical    DONOTHING 
      parameter (DONOTHING = .true.) 
      logical    reset_magnetic_field
      parameter (reset_magnetic_field = .false.) 

      integer     mem_alloc
      external    mem_alloc


      if(ltrace)write(*,99)myid,'grid_modify_data: Done grid: ',gridnum

  98  format('[',I3,'] ',A,3F18.6)
  99  format('[',I3,'] ',A,3I5)
      return
      end    ! END: grid_modify_data

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_test:                                                                cc
cc                   Called at the end of grid_initialize() which             cc
cc                   sets the data on a newly created amr grid.               cc
cc                   This routine is just for testing to see if anything      cc
cc                   is wrong with how it is set (such as a bad metric).      cc
cc                                                                            cc
cc   DWN: This was stub was copied from hyperBSSN_XMHD. It needs to be
cc        modified specifically for XMHD.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_test( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'grid_methods.inc'
      include    'glob.inc'
      include    'param.inc'



      return
      end    ! END: grid_test


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_tracer_update:                                                        cc
cc                   Takes as input the coordinate of a point located         cc
cc                   on the grid. It updatesthis point via convecting         cc
cc                   the particle with a "fluid".                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_tracer_update(gi,particlenum,dt)
      implicit     none
      integer      gi,particlenum
      real(kind=8) dt
      include     'grid.inc'
      include     'tracers.inc'
      integer      nx, ny, nz
      integer      i,j,k,index, l
      integer      ilow, jlow, klow
      real(kind=8) hg, x, y, z
      real(kind=8) velx, vely, velz
      real(kind=8) fracx, fracy, fracz
      real(kind=8) interp
      external     interp
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      ! Coordinates of particle:
      x = tracers_x(particlenum)
      y = tracers_y(particlenum)
      z = tracers_z(particlenum)

      ! Index of lower left point of interpolation cell:
      ilow  = INT( (x - gr_minx(gi)) / hg ) + 1
      jlow  = INT( (y - gr_miny(gi)) / hg ) + 1
      klow  = INT( (z - gr_minz(gi)) / hg ) + 1

      ! Fractional distance point is from lower left point:
      fracx = ( x - q( gr_x(gi) + (ilow-1) ) ) / hg
      fracy = ( y - q( gr_y(gi) + (jlow-1) ) ) / hg
      fracz = ( z - q( gr_z(gi) + (klow-1) ) ) / hg

      if (ltrace) then
         write(*,*) 'grid_tracer_update: Enter...dt=',dt
      end if
      if (ilow .lt. 1) then
         write(*,*) 'grid_tracer_update: x too small ',ilow
         call my_exit('grid_tracer_update: X too small')
      end if
      if (fracx.lt.-0.01 .or.fracy.lt.-0.01.or.fracz.lt.-0.0001) then
         write(*,*) 'grid_tracer_update: fraction cannot be negative'
         write(*,*) 'grid_tracer_update: fracx/y/z =',fracx,fracy,fracz
         call my_exit('grid_tracer_update: fraction cannot be negative')
      end if

      ! Project-dependent: Choose what your particle velocity is:
      velx = interp(  q(gr_shift1),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      vely = interp(  q(gr_shift2),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)
      velz = interp(  q(gr_shift3),
     *               fracx, fracy, fracz,
     *               ilow,jlow,klow,
     *               nx,ny,nz)

      tracers_x(particlenum) = tracers_x(particlenum) + dt*velx
      tracers_y(particlenum) = tracers_y(particlenum) + dt*vely
      tracers_z(particlenum) = tracers_z(particlenum) + dt*velz

      return
      end    ! END: grid_tracer_update

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  tracers_init_sub                                                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine tracers_init_sub()
       implicit     none
       include    'grid.inc'
       include    'grid_methods.inc'
       include    'glob.inc'
       include     'param.inc'
       include     'tracers.inc'
       include     'mpif.h'
       include     'mpi_stuff.inc'
       include     'chr.inc'
       !
       integer      i, j, k, l, mynx,myny,mynz, index
       integer      coarselevel, gi, numcells, offset,part
       integer      numparticles(0:numprocs-1), mynumparticles
       integer      myarray(0:numprocs-1)
       !integer      numparticles(0:10), mynumparticles
       !integer      myarray(0:10)
       integer      SEED
       !parameter (  SEED = 78237 )
       real*8       myx,myy,myz
       real*8       mydx, mydy, allsumrho, sumrho
       real*8       myrho
       !  The range of densities in which  to put particles:
       real*8       MINRHO,        MAXRHO
       parameter (  MINRHO = 5d-8, MAXRHO = 3d-1)
       ! values that will get overwritten with MPI call:
       real*8       UNDEFINED
       parameter (  UNDEFINED = 1d99)
       logical     ltrace
       parameter ( ltrace  = .true. )
       logical     ltrace2
       parameter ( ltrace2  = .false. )

      SEED = 78237

      if (ltrace) then
         write(*,99) myid, 'tracers_initsub:Enter: ',tracers_number
      end if

      if (.false.) then
         !
         ! something super simple just for testing:
         !
            mydx = (maxx0-minx0)/(tracers_number-1)
            do i = 1, tracers_number
               tracers_x(i) = minx0 + (i-1)*mydx
               tracers_y(i) = 0
               tracers_z(i) = 0.d0
            end do
         !
      else if (.true.) then
         !
         ! Distribute based on a range of density:
         !      Consider each (coarse level) gridpoint as the center of a cell of width h
         !      and plase n_i particles in that cell randomly so that n_i propto rho_i (the density)
         !      Since we want \Sigma n_i = N, that means we want:
         !            n_i = (N/\Sigma rho_i) rho_i
         !
         !      First compute Sigma rho_i on coarse level:
         !
         if (ltrace) then
            write(*,98)myid,'tracers_initsub: Low density regions'
            write(*,95)myid,'tracers_initsub: MINRHO=',MINRHO
            write(*,95)myid,'tracers_initsub: MAXRHO=',MAXRHO
         end if
         sumrho      = 0.d0
         numcells    = 0
         coarselevel = 0
         gi = level_return_start(coarselevel)
         if (.not.grid_return_existence(gi)) then
            write(*,99)myid,'tracers_init_sub:coarse level doesnt exist'
            return
         end if
 10      continue
         if (grid_return_existence(gi)) then
            if(ltrace)write(*,99)myid,'tracers_initsub:gi',gi
            if (grid_is_local(gi)) then
               !
               ! Add up rho values
               !
               mynx = gr_nx(gi)
               myny = gr_ny(gi)
               mynz = gr_nz(gi)
               call load_pointers(gi)
               do k = 1, mynz
               do j = 1, myny
               do i = 1, mynx
                  index = (k-1)*mynx*myny+(j-1)*mynx+(i-1)
                  myrho = q(gr_chi+index)
                  if ( (NINT(q(gr_chr+index)).eq.CHR_interior) .and.
     .                 (myrho.gt.MINRHO) .and.  (myrho.lt.MAXRHO) ) then
                     sumrho   = sumrho + q(gr_chi+index)
                     numcells = numcells + 1
                  end if
               end do
               end do
               end do
               if(ltrace)then 
                  write(*,96)myid,'tracers_initsub:gi/sumrho:',
     .                  gi,sumrho
                  write(*,96)myid,'tracers_initsub:numcells/sumrho:',
     .                  numcells,sumrho
               end if
            end if
            gi = grid_return_sibling(gi)
            goto 10
         end if
         !
         !      Second: get from other processors:
         !         --the total sum of rho
         !         --the number of cells each proc has
         !
         call MPI_Allreduce(sumrho,allsumrho,1,
     *          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
         if (ltrace) then
            write(*,98)myid,'tracers_initsub:allsumrho:',allsumrho
         end if
         if (abs(allsumrho).lt.1d-10)allsumrho=1.0

         !
         ! Compute mynumparticles
         !
         mynumparticles = 0
         gi = level_return_start(coarselevel)
         if (.not.grid_return_existence(gi)) then
            write(*,99)myid,'tracers_init_sub:coarse level doesnt exist'
            return
         end if
 15      continue
         if(ltrace)write(*,99)myid,'tracers_initsub:gi',gi
         if (grid_return_existence(gi)) then
            if (grid_is_local(gi)) then
               !
               ! Add Place particles in the desired cells
               !
               mynx = gr_nx(gi)
               myny = gr_ny(gi)
               mynz = gr_nz(gi)
               call load_pointers(gi)
               do k = 1, mynz
               do j = 1, myny
               do i = 1, mynx
                  index = (k-1)*mynx*myny+(j-1)*mynx+(i-1)
                  myrho = q(gr_chi+index)
                  if ( (NINT(q(gr_chr+index)).eq.CHR_interior) .and.
     .                 (myrho.gt.MINRHO) .and.  (myrho.lt.MAXRHO) ) then
                     mynumparticles = mynumparticles 
     .                            + NINT(tracers_number*myrho/allsumrho)
                  end if
               end do
               end do
               end do
               if(ltrace)write(*,99)myid,'tracers_inits:mynumparticles',
     .                              mynumparticles,gi
            end if
            gi = grid_return_sibling(gi)
            goto 15
         end if
         if(ltrace)write(*,99)myid,'tracers_inits:mynumparticles',
     .                              mynumparticles
            do i = 0,numprocs-1
               numparticles(i)=0
            end do
         call MPI_ALLGATHER(mynumparticles,  1,        MPI_INTEGER,
     .                      numparticles(0), 1,        MPI_INTEGER,
     .                      MPI_COMM_WORLD, ierr)
         if (ltrace) then
            write(*,98)myid,'tracers_initsub:Numparticles:'
            do i = 0,numprocs-1
               write(*,99)myid,'  numparticles:',i,numparticles(i)
            end do
         end if
         offset = 0
         do i = 0,myid-1
            offset = offset+numparticles(i)
         end do
         if (ltrace) then
            write(*,99)myid,'tracers_initsub:offset=',offset
         end if
         !
         !      Third: place the particles
         !
         part = 0
         gi = level_return_start(coarselevel)
         if (.not.grid_return_existence(gi)) then
            write(*,99)myid,'tracers_init_sub:coarse level doesnt exist'
            return
         end if
 20      continue
         if(ltrace)write(*,99)myid,'tracers_initsub:gi',gi
         if (grid_return_existence(gi)) then
            if(ltrace)write(*,99)myid,'tracers_initsub:gi exists',gi
            if (grid_is_local(gi)) then
               !
               ! Add Place particles in the desired cells
               !
               mynx = gr_nx(gi)
               myny = gr_ny(gi)
               mynz = gr_nz(gi)
               if(ltrace)write(*,99)myid,'tracers_initsub:gi local',mynx
     .                   ,myny, mynz
               call load_pointers(gi)
               do k = 1, mynz
               do j = 1, myny
               do i = 1, mynx
                  !
                  index = (k-1)*mynx*myny+(j-1)*mynx+(i-1)
                  myrho = q(gr_chi+index)
                  if ( (NINT(q(gr_chr+index)).eq.CHR_interior) .and.
     .                 (myrho.gt.MINRHO) .and.  (myrho.lt.MAXRHO) ) then
                     if(ltrace2)write(*,99)myid,'tracers_initsub: num:',
     .                   NINT(tracers_number*myrho/allsumrho),i,j,k
                     do l = 1, NINT(tracers_number*myrho/allsumrho)
                        if(ltrace2)write(*,99)myid,'tracers_initsub:Loo'
     .                                  ,l,gi
                        if (l+offset.gt.tracers_number) then
                         write(*,99)myid,'tracers_initsub:PROBLEM'
                        else
                           if(ltrace2)write(*,99)myid,'tracers_initsubl'
     .                                           ,l,gi,part
                           myx = gr_minx(gi)+gr_h(gi)*(i-1)
                           myy = gr_miny(gi)+gr_h(gi)*(j-1)
                           myz = gr_minz(gi)+gr_h(gi)*(k-1)
                           tracers_x(part+l+offset) = 
     .                                  myx+(ran(SEED)-0.5)*gr_h(gi)
!    .                                  myx+(1.0-0.5*ran(SEED))*gr_h(gi)
                           tracers_y(part+l+offset) = 
     .                                  myy+(ran(SEED)-0.5)*gr_h(gi)
                           tracers_z(part+l+offset) = 
     .                                  myz+(ran(SEED)-0.5)*gr_h(gi)
                           if(ltrace )write(*,94)myid,'tracers_initsub',
     .                                  part+l,
     .                                  tracers_x(part+l+offset),
     .                                  tracers_y(part+l+offset),
     .                                  tracers_z(part+l+offset)
                        end if
                     end do
                     part = part + NINT(tracers_number*myrho/allsumrho)
                  end if
                  !
               end do
               end do
               end do
               if(ltrace)write(*,99)myid,'tracers_initsub:gi/part:',
     .                      gi,part
            end if
            gi = grid_return_sibling(gi)
            goto 20
         end if
         !
      end if

      if (ltrace) then
         call field_dump_info1D(tracers_x, tracers_number, 'tracers_x')
         call field_dump_info1D(tracers_y, tracers_number, 'tracers_y')
         call field_dump_info1D(tracers_z, tracers_number, 'tracers_z')
         write(*,99) myid, 'tracers_init_sub: Done.', part
      end if

 94    format('[',I3,'] ',A,I7,'   ',5E12.5)
 95    format('[',I3,'] ',A,5E12.5)
 96    format('[',I3,'] ',A,I7,F12.6)
 98    format('[',I3,'] ',A,5F12.6)
 99    format('[',I3,'] ',A,5I7)

       return
       end     ! END: tracers_init_sub
