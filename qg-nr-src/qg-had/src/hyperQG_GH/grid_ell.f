!
!
! grid_ell.f:  These are a collection of project-specific routines
!              which address the elliptic problem. If the project
!              has no elliptic equation to solve, these should just
!              be dummy routines (declared, but no body necessary).
!
!              The idea is that the call sequences to the "normal" 
!              lop(), get_resid(), and relax() routines need to be 
!              hidden from the infrastructure and kept secured in  
!              the project directories.                             
!
!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_mg_trace:                                                            cc
cc            Output fields relevenat to the multigrid solver for tracing     cc
cc            and debugging.                                                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_mg_trace( gi, count )
      implicit    none
      integer     gi, count
      include    'grid.inc'
      include    'glob.inc'
      integer     nx,  ny,  nz, myint
      real*8      hg, time
      integer      proc_return_myid, myid, lev
      external     proc_return_myid
      logical     ltrace
      parameter ( ltrace = .false. )
      integer    i,j,k,index
      integer    nxc,nyc,nzc, kc,indexc
      character(2) tmps

      myid = proc_return_myid()

      call load_pointers(gi)

      if (ltrace) write(*,*) 'grid_mg_trace: On grid: ',gi, count

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)
      lev  = gr_level(gi)
      call int2str(lev,tmps)
      !time = 1.d0 * count
      ! Offset time w/r/t real time:
      myint = NINT(1000*gr_t(gi))
      time = 1.d0*count + 100000*myint

      !
      ! Output 2D field in the middle of the z-axis:
      !
      k    = (nz+1)/2
      index= (k-1)*nx*ny

      ! For the restricted fields:
      nxc    = NINT( 0.5d0*(nx+1) )
      nyc    = NINT( 0.5d0*(ny+1) )
      nzc    = NINT( 0.5d0*(nz+1) )
      kc     = (nzc+1)/2
      indexc = (kc-1)*nxc*nyc

c     call field_out2d( q(gr_myalpha+index),
c    *                  time,
c    *                  'MGmyalpha'//tmps,
c    *                  gr_minx(gi), gr_maxx(gi),
c    *                  gr_miny(gi), gr_maxy(gi),
c    *                  nx,ny, myid)

      return
      end    ! END: grid_mg_trace
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_get_resid:                                                           cc
cc            Compute the residual of the current solution and stick          cc
cc            it into the u_st1 arrays.                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_get_resid( gi )
      implicit    none
      integer     gi
      include    'grid.inc'
      integer     nx,  ny,  nz
      integer     nxc, nyc, nzc
      real*8      hg, hc
      logical     ltrace
      parameter ( ltrace = .false. )
      integer    i,j,k,index

      call load_pointers(gi)

      if (ltrace) write(*,*) 'grid_get_resid: On grid: ',gi

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      nxc = NINT( 0.5d0*(nx+1) )
      nyc = NINT( 0.5d0*(ny+1) )
      nzc = NINT( 0.5d0*(nz+1) )
      hc  = 2.d0*hg

c     call get_resid(q(gr_psi_elliptic_st1),
c    *               q(gr_psi_elliptic),
c    *               q(gr_psi_elliptic_rhs),hg,nx,ny,nz)

      call get_resid(q(gr_psi_el_st1), q(gr_psi_el), q(gr_psi_el_rhs),
     *           q(gr_betax_el_st1), q(gr_betax_el), q(gr_betax_el_rhs),
     *           q(gr_betay_el_st1), q(gr_betay_el), q(gr_betay_el_rhs),
     *           q(gr_betaz_el_st1), q(gr_betaz_el), q(gr_betaz_el_rhs),
     *           q(gr_chr), hg,nx,ny,nz)

      return
      end    ! END: grid_get_resid

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_lop:                                                                 cc
cc            Compute the results of the differential operator associated     cc
cc            with the elliptic problem acting on the set of restricted       cc
cc            fields residing in the u_st1 fields and *add* the result        cc
cc            to the u_st2 arrays.                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_lop( gi )
      implicit    none
      integer     gi
      include    'grid.inc'
      integer     nx,  ny,  nz
      integer     nxc, nyc, nzc
      real*8      hg, hc
      logical     ltrace
      parameter ( ltrace = .false. )
      integer    i,j,k,index

      call load_pointers(gi)

      if (ltrace) write(*,*) 'grid_lop:    On grid: ',gi

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      nxc = NINT( 0.5d0*(nx+1) )
      nyc = NINT( 0.5d0*(ny+1) )
      nzc = NINT( 0.5d0*(nz+1) )
      hc  = 2.d0*hg

      ! Restrict an auxiliary fields
      ! None in this project
      !call restrict(  q(gr_d333_rk2), q(gr_d333_rk1), nxc,nyc,nzc)

c     call lop( q(gr_psi_elliptic_st2),
c    *          q(gr_psi_elliptic_st1),hc,nxc,nyc,nzc)

      call lop(  q(gr_psi_el_st2),   q(gr_psi_el_st1),
     *           q(gr_betax_el_st2), q(gr_betax_el_st1),
     *           q(gr_betay_el_st2), q(gr_betay_el_st1),
     *           q(gr_betaz_el_st2), q(gr_betaz_el_st1),
     *           hc,nxc,nyc,nzc)

      return
      end    ! END: grid_lop

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_ellrhs:                                                              cc
cc               Compute the RHS for the elliptic solve.                      cc
cc               If no elliptic solver in project, just leave empty.          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_ellrhs( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      integer     nx, ny, nz
      real*8      hg, jx, jy, jz
      logical     ltrace
      parameter ( ltrace = .false. )
      integer     i,j,k,index
      real*8      cpi
      parameter ( cpi  =  3.14159 26535 89793 23846 26433 83279 d0)
      real*8      h11, h12, h13, h22, h23, h33, deth,
     *            T00, T01, T02, T03, T11, T12, T13, T22, T23, T33,
     *            Tdu01, Tdu11, Tdu21, Tdu31, Tdu02, Tdu12,
     *            Tdu22, Tdu32, Tdu03, Tdu13, Tdu23, Tdu33,
     *            huu11, huu12, huu13, huu22, huu23, huu33,
     *            nd0, nd1, nd2, nd3,
     *            nu0, nu1, nu2, nu3,
     *            b1, b2, b3, alp, alp2,
     *            g00, g01, g02, g03,
     *            guu00, guu01, guu02, guu03, guu11, guu12, guu13,
     *            guu22, guu23, guu33

      call load_pointers(gridnum)

      if (ltrace) write(*,*) 'grid_ellrhs:    On grid: ',gridnum

      nx   = gr_nx(gridnum)
      ny   = gr_ny(gridnum)
      nz   = gr_nz(gridnum)
      hg   = gr_h(gridnum)

      ! Initialize guesses for elliptic fields:
      call load_scal3d( q(gr_psi_el),   1.d0, nx,ny,nz )
      if (.false.) then
         ! Either start with zero shift....
         call load_scal3d( q(gr_betax_el), 0.d0, nx,ny,nz )
         call load_scal3d( q(gr_betay_el), 0.d0, nx,ny,nz )
         call load_scal3d( q(gr_betaz_el), 0.d0, nx,ny,nz )
      else
         ! Or use what we get from the superposition
         call mat_copy3d( q(gr_shift1), q(gr_betax_el), nx,ny,nz)
         call mat_copy3d( q(gr_shift2), q(gr_betay_el), nx,ny,nz)
         call mat_copy3d( q(gr_shift3), q(gr_betaz_el), nx,ny,nz)
      end if

      !
      ! Compute interior:
      !
      do k = 2, nz-1
         do j = 2, ny-1
            do i = 2, nx-1
               index    = (k-1)*ny*nx + (j-1)*nx + (i-1)
               !
               g00   = q(gr_g00+index)
               g01   = q(gr_g01+index)
               g02   = q(gr_g02+index)
               g03   = q(gr_g03+index)
               !
               h11   = q(gr_g11+index)
               h12   = q(gr_g12+index)
               h13   = q(gr_g13+index)
               h22   = q(gr_g22+index)
               h23   = q(gr_g23+index)
               h33   = q(gr_g33+index)
               !
               ! Stress-energy Tensor (down-down componenents):
               !
               T00   = q(gr_Tmunu00+index)
               T01   = q(gr_Tmunu01+index)
               T02   = q(gr_Tmunu02+index)
               T03   = q(gr_Tmunu03+index)
               T11   = q(gr_Tmunu11+index)
               T12   = q(gr_Tmunu12+index)
               T13   = q(gr_Tmunu13+index)
               T22   = q(gr_Tmunu22+index)
               T23   = q(gr_Tmunu23+index)
               T33   = q(gr_Tmunu33+index)
               !
               huu11 = -h23**2 + h22*h33
               huu12 = h13*h23 - h12*h33
               huu13 = -(h13*h22) + h12*h23
               huu22 = -h13**2 + h11*h33
               huu23 = h12*h13 - h11*h23
               huu33 = -h12**2 + h11*h22
               deth  = h11*huu11 + h12*huu12 + h13*huu13
               huu11 = huu11/deth
               huu12 = huu12/deth
               huu13 = huu13/deth
               huu22 = huu22/deth
               huu23 = huu23/deth
               huu33 = huu33/deth
               !
               b1    = g01*huu11 + g02*huu12 + g03*huu13
               b2    = g01*huu12 + g02*huu22 + g03*huu23
               b3    = g01*huu13 + g02*huu23 + g03*huu33
               alp2  = -g00 + b1**2*h11 + b2**2*h22
     *                 + 2*(b1*(b2*h12 + b3*h13)
     *                 + b2*b3*h23) + b3**2*h33
               alp   = sqrt(alp2)
               !
               guu00 = -(1.0d0/alp2)
               guu01 = b1/alp2
               guu02 = b2/alp2
               guu03 = b3/alp2
               guu11 = -(b1**2/alp2) + huu11
               guu12 = -((b1*b2)/alp2) + huu12
               guu13 = -((b1*b3)/alp2) + huu13
               guu22 = -(b2**2/alp2) + huu22
               guu23 = -((b2*b3)/alp2) + huu23
               guu33 = -(b3**2/alp2) + huu33  
               !
               ! Stress-energy Tensor (down-up componenents):
               !
               Tdu01   = guu01*T00 + guu11*T01 + guu12*T02 + guu13*T03
               Tdu11   = guu01*T01 + guu11*T11 + guu12*T12 + guu13*T13
               Tdu21   = guu01*T02 + guu11*T12 + guu12*T22 + guu13*T23
               Tdu31   = guu01*T03 + guu11*T13 + guu12*T23 + guu13*T33
               Tdu02   = guu02*T00 + guu12*T01 + guu22*T02 + guu23*T03
               Tdu12   = guu02*T01 + guu12*T11 + guu22*T12 + guu23*T13
               Tdu22   = guu02*T02 + guu12*T12 + guu22*T22 + guu23*T23
               Tdu32   = guu02*T03 + guu12*T13 + guu22*T23 + guu23*T33
               Tdu03   = guu03*T00 + guu13*T01 + guu23*T02 + guu33*T03
               Tdu13   = guu03*T01 + guu13*T11 + guu23*T12 + guu33*T13
               Tdu23   = guu03*T02 + guu13*T12 + guu23*T22 + guu33*T23
               Tdu33   = guu03*T03 + guu13*T13 + guu23*T23 + guu33*T33
               !
               ! Unit normal to hypersurface:
               !
               nd0 = -1.d0
               nd1 =  0.d0
               nd2 =  0.d0
               nd3 =  0.d0
               !
               nu0 =  1.d0/1.d0
               nu1 = -b1/1.d0
               nu2 = -b2/1.d0
               nu3 = -b3/1.d0
               !
               ! These current components to be set for
               !  psi  --> - 2 Pi rho
               !  betax--> +16 Pi j^x
               !  betay--> +16 Pi j^y
               !  betaz--> +16 Pi j^z
               !
!              jx    =  16.d0*cPi*(
!    *                     nd0*Tuu01 + nd1*Tuu11 + nd2*Tuu12 + nd3*Tuu13)
!              jy    =  16.d0*cPi*(
!    *                     nd0*Tuu02 + nd1*Tuu12 + nd2*Tuu22 + nd3*Tuu13)
!              jz    =  16.d0*cPi*(
!    *                     nd0*Tuu03 + nd1*Tuu13 + nd2*Tuu23 + nd3*Tuu33)
               !
               jx    =  16.d0*cPi*(
     *                    nu0*Tdu01 + nu1*Tdu11 + nu2*Tdu21 + nu3*Tdu31)
               jy    =  16.d0*cPi*(
     *                    nu0*Tdu02 + nu1*Tdu12 + nu2*Tdu22 + nu3*Tdu32)
               jz    =  16.d0*cPi*(
     *                    nu0*Tdu03 + nu1*Tdu13 + nu2*Tdu23 + nu3*Tdu33)
               !
               ! I am worried that this is rho times sqrt(g) not rho:
               q(gr_psi_el_rhs  +index) = 
     *                              -2.d0*cPi*q(gr_energy_dens+index)
               q(gr_betax_el_rhs+index) = 16.d0 * cPi * jx
               q(gr_betay_el_rhs+index) = 16.d0 * cPi * jy
               q(gr_betaz_el_rhs+index) = 16.d0 * cPi * jz
               !
            end do
         end do
      end do

      !
      ! For debugging only set to a Gaussian:
      !
c     call gaussian3d( q(gr_psi_elliptic_rhs), hg, gr_minx(gridnum),
c    *     gr_miny(gridnum), gr_minz(gridnum),nx,ny,nz,
c    *     1.d0,1.d0,1.d0,0.d0,0.1d0,0.d0,0.d0, 0.d0)

c     if (ltrace) call field_dump_info(q(gr_psi_elliptic_rhs),nx,ny,nz)

      return
      end    ! END: grid_ellrhs

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_relax:                                                               cc
cc            Apply relaxation routine to field(s) on a given grid.           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_relax(gi)
      implicit      none
      integer       gi
      include      'grid.inc'
      include      'glob.inc'

      real(kind=16)   relax
      external       relax
      real(kind=8)   resid

      logical     ltrace
      parameter ( ltrace  = .false. )


      if (ltrace) then
         write(*,*) 'grid_relax: Relaxing grid: ',gi
      end if

      call load_pointers(gi)

c     resid = relax(q(gr_psi_elliptic),q(gr_psi_elliptic_rhs),q(gr_chr),
c    *                gr_h(gi),gr_nx(gi),gr_ny(gi),gr_nz(gi))

      resid = relax(q(gr_psi_el),   q(gr_psi_el_rhs),
     *              q(gr_betax_el), q(gr_betax_el_rhs),
     *              q(gr_betay_el), q(gr_betay_el_rhs),
     *              q(gr_betaz_el), q(gr_betaz_el_rhs),
     *              q(gr_chr), gr_h(gi),gr_nx(gi),gr_ny(gi),gr_nz(gi))

      !
      ! For output and debugging purposes only:
      !
!     call get_resid(q(gr_phi_x),q(gr_phi),q(gr_rho),
!    *                gr_h(gi),gr_nx(gi),gr_ny(gi),gr_nz(gi))

      if (ltrace) then
         write(*,*) 'grid_relax: Finished with resid: ',resid
      end if

      return
      end      ! END: grid_relax

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_resnrm:  Get norm of residual on a given grid (for debugging):cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function grid_return_resnrm(gi)
      implicit      none
      integer       gi
      include      'grid.inc'
      integer       nx, ny, nz

      real(kind=8)  h
      real(kind=16)  get_resnrm
      external      get_resnrm
      real(kind=8)  grid_return_h
      external      grid_return_h
      logical       grid_is_local
      external      grid_is_local

      logical     ltrace
      parameter ( ltrace = .false. )


      call grid_get_dims(gi,nx,ny,nz)
      h = grid_return_h(gi)

      call load_pointers(gi)

      if (.not. grid_is_local(gi)) then
         write(*,*) 'grid_return_resnrm: Only works on local grids'
      end if

c     grid_return_resnrm = 0.d0
c     grid_return_resnrm = get_resnrm(q(gr_psi_elliptic),
c    *                          q(gr_psi_elliptic_rhs),h,nx,ny,nz)

      grid_return_resnrm = get_resnrm(
     *                        q(gr_psi_el),   q(gr_psi_el_rhs),
     *                        q(gr_betax_el), q(gr_betax_el_rhs),
     *                        q(gr_betay_el), q(gr_betay_el_rhs),
     *                        q(gr_betaz_el), q(gr_betaz_el_rhs),
     *                        h,nx,ny,nz)

      return
      end      ! END: grid_return_resnrm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_ellpost:                                                             cc
cc               Things to do after an elliptic solve.                        cc
cc               For MHD project, this is where we clean the B field.         cc
cc               For Maximal Slicing  project, this is where copy into alpha. cc
cc               For CTS ID, this is where we reconstruct the initial metric. cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_ellpost( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'param.inc'
      integer     nx, ny, nz
      real*8      hg
      logical     ltrace
      parameter ( ltrace = .true. )
      integer     i,j,k,index
      real(kind=8)  tmp

      call load_pointers(gridnum)

      if (ltrace) write(*,*) 'grid_ellpost:   On grid: ',gridnum

      nx   = gr_nx(gridnum)
      ny   = gr_ny(gridnum)
      nz   = gr_nz(gridnum)
      hg   = gr_h(gridnum)

      ! Spatial metric:
      call load_scal3d( q(gr_g12), 0.d0, nx,ny,nz )
      call load_scal3d( q(gr_g13), 0.d0, nx,ny,nz )
      call load_scal3d( q(gr_g23), 0.d0, nx,ny,nz )

      ! Not sure if this is necessary:
      call load_scal3d( q(gr_alpha), 1.d0, nx,ny,nz )

      ! Set shift vector for evolution part with elliptic solution:
      call mat_copy3d( q(gr_betax_el), q(gr_shift1), nx,ny,nz)
      call mat_copy3d( q(gr_betay_el), q(gr_shift2), nx,ny,nz)
      call mat_copy3d( q(gr_betaz_el), q(gr_shift3), nx,ny,nz)

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         index           = (k-1)*ny*nx + (j-1)*nx + (i-1)
         tmp             = q(gr_psi_el+index)**4
         ! Diagonal of spatial metric:
         q(gr_g11+index) = q(gr_g11+index) * tmp
         q(gr_g22+index) = q(gr_g22+index) * tmp
         q(gr_g33+index) = q(gr_g33+index) * tmp
         ! 
         q(gr_g00+index) = -1.d0
     *                     + q(gr_shift1+index)**2
     *                     + q(gr_shift2+index)**2
     *                     + q(gr_shift3+index)**2
         q(gr_g01+index) =   q(gr_g11+index) * q(gr_shift1+index)
     *                     + q(gr_g12+index) * q(gr_shift2+index)
     *                     + q(gr_g13+index) * q(gr_shift3+index)
         q(gr_g02+index) =   q(gr_g12+index) * q(gr_shift1+index)
     *                     + q(gr_g22+index) * q(gr_shift2+index)
     *                     + q(gr_g23+index) * q(gr_shift3+index)
         q(gr_g03+index) =   q(gr_g13+index) * q(gr_shift1+index)
     *                     + q(gr_g23+index) * q(gr_shift2+index)
     *                     + q(gr_g33+index) * q(gr_shift3+index)
      end do
      end do
      end do

      ! Need to set Q_ij and rest of metric terms

      ! Do metric derivatives need to be computed?

      return
      end    ! END: grid_ellpost

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  get_resid:                                                                cc
cc             Compute:   res = rhs - L u                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_resid(res_psi, psi, rhs_psi,
     *                     res_bx,  bx,  rhs_bx,
     *                     res_by,  by,  rhs_by,
     *                     res_bz,  bz,  rhs_bz,
     *                     chr, h, nx, ny, nz)
      implicit      none
      integer       nx, ny, nz
      real(kind=8)  res_psi(nx,ny,nz), psi(nx,ny,nz), rhs_psi(nx,ny,nz)
      real(kind=8)  res_bx(nx,ny,nz),  bx(nx,ny,nz),  rhs_bx(nx,ny,nz)
      real(kind=8)  res_by(nx,ny,nz),  by(nx,ny,nz),  rhs_by(nx,ny,nz)
      real(kind=8)  res_bz(nx,ny,nz),  bz(nx,ny,nz),  rhs_bz(nx,ny,nz)
      real(kind=8)  chr(nx,ny,nz),     h
      include      'chr.inc'
      integer       i,j,k

      !
      !  Compute: L u
      !
      call load_scal1D(res_psi,0.d0,nx*ny*nz)
      call load_scal1D(res_bx, 0.d0,nx*ny*nz)
      call load_scal1D(res_by, 0.d0,nx*ny*nz)
      call load_scal1D(res_bz, 0.d0,nx*ny*nz)

      !
      ! Compute: Lu
      !
      call lop(  res_psi, psi, 
     *           res_bx,  bx,
     *           res_by,  by,
     *           res_bz,  bz,
     *          h, nx,ny,nz)

      !
      !  Compute: rhs - L u
      !
      call mat_mat_sub3d( res_psi, rhs_psi, res_psi, nx,ny,nz)
      call mat_mat_sub3d( res_bx,  rhs_bx,  res_bx,  nx,ny,nz)
      call mat_mat_sub3d( res_by,  rhs_by,  res_by,  nx,ny,nz)
      call mat_mat_sub3d( res_bz,  rhs_bz,  res_bz,  nx,ny,nz)

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx, nx-1
         if (NINT(chr(i,j,k)).eq.CHR_deco_bdy) then
            res_bx(i,j,k)  = 0.d0
            res_by(i,j,k)  = 0.d0
            res_bz(i,j,k)  = 0.d0
         end if
      end do
      end do
      end do
      do k = 1, nz
      do i = 1, nx
      do j = 1, ny, ny-1
         if (NINT(chr(i,j,k)).eq.CHR_deco_bdy) then
            res_psi(i,j,k) = 0.d0
            res_bx(i,j,k)  = 0.d0
            res_by(i,j,k)  = 0.d0
            res_bz(i,j,k)  = 0.d0
         end if
      end do
      end do
      end do
      do j = 1, ny
      do i = 1, nx
      do k = 1, nz, nz-1
         if (NINT(chr(i,j,k)).eq.CHR_deco_bdy) then
            res_psi(i,j,k) = 0.d0
            res_bx(i,j,k)  = 0.d0
            res_by(i,j,k)  = 0.d0
            res_bz(i,j,k)  = 0.d0
         end if
      end do
      end do
      end do

      return
      end      ! END: get_resid

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  get_resnrm:                                                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function get_resnrm(psi,rhs_psi,
     *                                     bx, rhs_bx,
     *                                     by, rhs_by,
     *                                     bz, rhs_bz,
     *                                     h,nx,ny,nz)
      implicit      none
      integer       nx, ny, nz
      real(kind=8)  psi(nx,ny,nz),  rhs_psi(nx,ny,nz), h,
     *              bx(nx,ny,nz),   rhs_bx(nx,ny,nz),
     *              by(nx,ny,nz),   rhs_by(nx,ny,nz),
     *              bz(nx,ny,nz),   rhs_bz(nx,ny,nz)
      !
      integer       i,j,k
      !
      real(kind=8)  lhs_psi, lhs_bx, lhs_by, lhs_bz
      real(kind=8)  d2psi_dx2, d2psi_dy2, d2psi_dz2,
     *              d2bx_dx2,  d2bx_dy2,  d2bx_dz2,
     *              d2by_dx2,  d2by_dy2,  d2by_dz2,
     *              d2bz_dx2,  d2bz_dy2,  d2bz_dz2,
     *              Ad11, Ad12, Ad13, Au11, Au12, Au13,
     *                    Ad22, Ad23,       Au22, Au23,
     *                          Ad33,             Au33,
     *              d2bx_dxdy,            d2bx_dxdz,
     *              d2by_dxdy, d2by_dydz,
     *                         d2bz_dydz, d2bz_dxdz,
     *              sum_dbk_dk
      !

      get_resnrm = 0.d0

      do k = 2, nz-1
         do j = 2, ny-1
            do i = 2, nx-1
               !
                  d2psi_dx2 =
     *              (psi(i+1,j,k) - 2.d0*psi(i,j,k) + psi(i-1,j,k))/h**2
                  d2psi_dy2 =
     *              (psi(i,j+1,k) - 2.d0*psi(i,j,k) + psi(i,j-1,k))/h**2
                  d2psi_dz2 =
     *              (psi(i,j,k+1) - 2.d0*psi(i,j,k) + psi(i,j,k-1))/h**2
                  !
                  d2bx_dx2  =
     *              (bx(i+1,j,k) - 2.d0*bx(i,j,k) + bx(i-1,j,k))/h**2
                  d2bx_dy2  =
     *              (bx(i,j+1,k) - 2.d0*bx(i,j,k) + bx(i,j-1,k))/h**2
                  d2bx_dz2  =
     *              (bx(i,j,k+1) - 2.d0*bx(i,j,k) + bx(i,j,k-1))/h**2
                  !
                  d2by_dx2  =
     *              (by(i+1,j,k) - 2.d0*by(i,j,k) + by(i-1,j,k))/h**2
                  d2by_dy2  =
     *              (by(i,j+1,k) - 2.d0*by(i,j,k) + by(i,j-1,k))/h**2
                  d2by_dz2  =
     *              (by(i,j,k+1) - 2.d0*by(i,j,k) + by(i,j,k-1))/h**2
                  !
                  d2bz_dx2  =
     *              (bz(i+1,j,k) - 2.d0*bz(i,j,k) + bz(i-1,j,k))/h**2
                  d2bz_dy2  =
     *              (bz(i,j+1,k) - 2.d0*bz(i,j,k) + bz(i,j-1,k))/h**2
                  d2bz_dz2  =
     *              (bz(i,j,k+1) - 2.d0*bz(i,j,k) + bz(i,j,k-1))/h**2
                  !
                  d2bx_dxdy =  (  (bx(i+1,j+1,k) - bx(i+1,j-1,k))
     *                           -(bx(i-1,j+1,k) - bx(i-1,j-1,k))
     *                           ) / (4.d0*h**2)
                  d2bx_dxdz =  (  (bx(i+1,j,k+1) - bx(i+1,j,k-1))
     *                           -(bx(i-1,j,k+1) - bx(i-1,j,k-1))
     *                           ) / (4.d0*h**2)
                  !
                  d2by_dxdy =  (  (by(i+1,j+1,k) - by(i+1,j-1,k))
     *                           -(by(i-1,j+1,k) - by(i-1,j-1,k))
     *                           ) / (4.d0*h**2)
                  d2by_dydz =  (  (by(i,j+1,k+1) - by(i,j+1,k-1))
     *                           -(by(i,j-1,k+1) - by(i,j-1,k-1))
     *                           ) / (4.d0*h**2)
                  !
                  d2bz_dxdz =  (  (bz(i+1,j,k+1) - bz(i+1,j,k-1))
     *                           -(bz(i-1,j,k+1) - bz(i-1,j,k-1))
     *                           ) / (4.d0*h**2)
                  d2bz_dydz =  (  (bz(i,j+1,k+1) - bz(i,j+1,k-1))
     *                           -(bz(i,j-1,k+1) - bz(i,j-1,k-1))
     *                           ) / (4.d0*h**2)
                  !
                  sum_dbk_dk = (   bx(i+1,j,k)-bx(i-1,j,k)
     *                           + by(i,j+1,k)-by(i,j-1,k)
     *                           + by(i,j,k+1)-by(i,j,k-1)
     *                           ) / (2.d0*h)
                  !
                  Ad11 = 0.5d0*( 
     *                2.d0*(   psi(i+1,j,k)**4*bx(i+1,j,k)
     *                       - psi(i-1,j,k)**4*bx(i-1,j,k) )/(2.d0*h)
     *               -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Ad22 = 0.5d0*( 
     *                2.d0*(   psi(i,j+1,k)**4*by(i,j+1,k)
     *                       - psi(i,j-1,k)**4*by(i,j-1,k) )/(2.d0*h)
     *               -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Ad33 = 0.5d0*( 
     *                2.d0*(   psi(i,j,k+1)**4*bz(i,j,k+1)
     *                       - psi(i,j,k-1)**4*bz(i,j,k-1) )/(2.d0*h)
     *               -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  !
                  Ad12 = 0.5d0*( 
     *                     (   psi(i+1,j,k)**4*by(i+1,j,k)
     *                       - psi(i-1,j,k)**4*by(i-1,j,k) )/(2.d0*h)
     *                   + (   psi(i,j+1,k)**4*bx(i,j+1,k)
     *                       - psi(i,j-1,k)**4*bx(i,j-1,k) )/(2.d0*h)
     *                          )
                  Ad13 = 0.5d0*( 
     *                     (   psi(i+1,j,k)**4*bz(i+1,j,k)
     *                       - psi(i-1,j,k)**4*bz(i-1,j,k) )/(2.d0*h)
     *                   + (   psi(i,j,k+1)**4*bx(i,j,k+1)
     *                       - psi(i,j,k-1)**4*bx(i,j,k-1) )/(2.d0*h)
     *                          )
                  Ad23 = 0.5d0*( 
     *                     (   psi(i,j+1,k)**4*bz(i,j+1,k)
     *                       - psi(i,j-1,k)**4*bz(i,j-1,k) )/(2.d0*h)
     *                   + (   psi(i,j,k+1)**4*by(i,j,k+1)
     *                       - psi(i,j,k-1)**4*by(i,j,k-1) )/(2.d0*h)
     *                          )
                  !
                  Au11 = 0.5d0*(  
     *                        2.d0*(bx(i+1,j,k)-bx(i-1,j,k))/(2.d0*h)
     *                             / psi(i,j,k)**4
     *                      -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Au22 = 0.5d0*(  
     *                        2.d0*(by(i,j+1,k)-by(i,j-1,k))/(2.d0*h)
     *                             / psi(i,j,k)**4
     *                      -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Au33 = 0.5d0*(  
     *                        2.d0*(bz(i,j,k+1)-bz(i,j,k-1))/(2.d0*h)
     *                             / psi(i,j,k)**4
     *                      -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  !
                  Au12 = 0.5d0*(   (by(i+1,j,k)-by(i-1,j,k))/(2.d0*h)
     *                           + (bx(i,j+1,k)-bx(i,j-1,k))/(2.d0*h)
     *                          ) / psi(i,j,k)**4
                  Au13 = 0.5d0*(   (bz(i+1,j,k)-bz(i-1,j,k))/(2.d0*h)
     *                           + (bx(i,j,k+1)-bx(i,j,k-1))/(2.d0*h)
     *                          ) / psi(i,j,k)**4
                  Au23 = 0.5d0*(   (bz(i,j+1,k)-bz(i,j-1,k))/(2.d0*h)
     *                           + (by(i,j,k+1)-by(i,j,k-1))/(2.d0*h)
     *                          ) / psi(i,j,k)**4
                  !
                  !---------------------------------
                  !
                  !    Psi Equation:
                  !
                  !
                  lhs_psi  = (    d2psi_dx2 + d2psi_dy2 + d2psi_dz2
     *                         +  (1.d0 / 8.d0 / psi(i,j,k)**7)
     *                           *(   Ad11*Au11 + Ad12*Au12 + Ad13*Au13
     *                              + Ad12*Au12 + Ad22*Au22 + Ad23*Au23
     *                              + Ad13*Au13 + Ad23*Au23 + Ad33*Au33)
     *                         ) / psi(i,j,k)**5
                  !
                  !---------------------------------
                  !
                  !    Bx Equation:
                  !
                  lhs_bx  = (     d2bx_dx2 + d2bx_dy2 + d2bx_dz2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *(d2bx_dx2 + d2by_dxdy + d2bz_dxdz )
     *                         ) / psi(i,j,k)**(10)
                  !
                  !---------------------------------
                  !
                  !    By Equation:
                  !
                  lhs_by  = (     d2by_dx2 + d2by_dy2 + d2by_dz2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *(d2bx_dxdy + d2by_dy2 + d2bz_dydz )
     *                         ) / psi(i,j,k)**(10)
                  !
                  !---------------------------------
                  !
                  !    Bz Equation:
                  !
                  lhs_bz  = (     d2bz_dx2 + d2bz_dy2 + d2bz_dz2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *(d2bx_dxdz + d2by_dydz + d2bz_dz2 )
     *                         ) / psi(i,j,k)**(10)
               !
               get_resnrm = get_resnrm
     *                              + (lhs_psi - rhs_psi(i,j,k))**2
     *                              + (lhs_bx  - rhs_bx(i,j,k))**2
     *                              + (lhs_by  - rhs_by(i,j,k))**2
     *                              + (lhs_bz  - rhs_bz(i,j,k))**2
            end do
         end do
      end do

      get_resnrm = sqrt(get_resnrm / nx /ny /nz / 4.d0)

      return
      end       ! END: get_resnrm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  lop:  Needs to parallel the finite differencing in relax()                cc
cc        NB: Needs to add the calculation to what resides in                 cc
cc            storage already!                                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lop( result_psi, psi,
     *                result_bx,  bx,
     *                result_by,  by,
     *                result_bz,  bz,
     *                h,nx,ny,nz,gw)
      implicit      none
      integer       nx, ny, nz, gw
      real(kind=8)  h
      real(kind=8)  result_psi(nx,ny,nz), psi(nx,ny,nz), 
     *              result_bx(nx,ny,nz),  bx(nx,ny,nz), 
     *              result_by(nx,ny,nz),  by(nx,ny,nz), 
     *              result_bz(nx,ny,nz),  bz(nx,ny,nz)
      !
      integer       i, j, k
      !
      real(kind=8)  lhs_psi, lhs_bx, lhs_by, lhs_bz
      real(kind=8)  d2psi_dx2, d2psi_dy2, d2psi_dz2,
     *              d2bx_dx2,  d2bx_dy2,  d2bx_dz2,
     *              d2by_dx2,  d2by_dy2,  d2by_dz2,
     *              d2bz_dx2,  d2bz_dy2,  d2bz_dz2,
     *              Ad11, Ad12, Ad13, Au11, Au12, Au13,
     *                    Ad22, Ad23,       Au22, Au23,
     *                          Ad33,             Au33,
     *              d2bx_dxdy,            d2bx_dxdz,
     *              d2by_dxdy, d2by_dydz,
     *                         d2bz_dydz, d2bz_dxdz,
     *              sum_dbk_dk
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'lop:'
      end if

      do k = 2, nz-1
         do j = 2, ny-1
            do i = 2, nx-1
               !
                  d2psi_dx2 =
     *              (psi(i+1,j,k) - 2.d0*psi(i,j,k) + psi(i-1,j,k))/h**2
                  d2psi_dy2 =
     *              (psi(i,j+1,k) - 2.d0*psi(i,j,k) + psi(i,j-1,k))/h**2
                  d2psi_dz2 =
     *              (psi(i,j,k+1) - 2.d0*psi(i,j,k) + psi(i,j,k-1))/h**2
                  !
                  d2bx_dx2  =
     *              (bx(i+1,j,k) - 2.d0*bx(i,j,k) + bx(i-1,j,k))/h**2
                  d2bx_dy2  =
     *              (bx(i,j+1,k) - 2.d0*bx(i,j,k) + bx(i,j-1,k))/h**2
                  d2bx_dz2  =
     *              (bx(i,j,k+1) - 2.d0*bx(i,j,k) + bx(i,j,k-1))/h**2
                  !
                  d2by_dx2  =
     *              (by(i+1,j,k) - 2.d0*by(i,j,k) + by(i-1,j,k))/h**2
                  d2by_dy2  =
     *              (by(i,j+1,k) - 2.d0*by(i,j,k) + by(i,j-1,k))/h**2
                  d2by_dz2  =
     *              (by(i,j,k+1) - 2.d0*by(i,j,k) + by(i,j,k-1))/h**2
                  !
                  d2bz_dx2  =
     *              (bz(i+1,j,k) - 2.d0*bz(i,j,k) + bz(i-1,j,k))/h**2
                  d2bz_dy2  =
     *              (bz(i,j+1,k) - 2.d0*bz(i,j,k) + bz(i,j-1,k))/h**2
                  d2bz_dz2  =
     *              (bz(i,j,k+1) - 2.d0*bz(i,j,k) + bz(i,j,k-1))/h**2
                  !
                  d2bx_dxdy =  (  (bx(i+1,j+1,k) - bx(i+1,j-1,k))
     *                           -(bx(i-1,j+1,k) - bx(i-1,j-1,k))
     *                           ) / (4.d0*h**2)
                  d2bx_dxdz =  (  (bx(i+1,j,k+1) - bx(i+1,j,k-1))
     *                           -(bx(i-1,j,k+1) - bx(i-1,j,k-1))
     *                           ) / (4.d0*h**2)
                  !
                  d2by_dxdy =  (  (by(i+1,j+1,k) - by(i+1,j-1,k))
     *                           -(by(i-1,j+1,k) - by(i-1,j-1,k))
     *                           ) / (4.d0*h**2)
                  d2by_dydz =  (  (by(i,j+1,k+1) - by(i,j+1,k-1))
     *                           -(by(i,j-1,k+1) - by(i,j-1,k-1))
     *                           ) / (4.d0*h**2)
                  !
                  d2bz_dxdz =  (  (bz(i+1,j,k+1) - bz(i+1,j,k-1))
     *                           -(bz(i-1,j,k+1) - bz(i-1,j,k-1))
     *                           ) / (4.d0*h**2)
                  d2bz_dydz =  (  (bz(i,j+1,k+1) - bz(i,j+1,k-1))
     *                           -(bz(i,j-1,k+1) - bz(i,j-1,k-1))
     *                           ) / (4.d0*h**2)
                  !
                  sum_dbk_dk = (   bx(i+1,j,k)-bx(i-1,j,k)
     *                           + by(i,j+1,k)-by(i,j-1,k)
     *                           + by(i,j,k+1)-by(i,j,k-1)
     *                           ) / (2.d0*h)
                  !
                  Ad11 = 0.5d0*( 
     *                2.d0*(   psi(i+1,j,k)**4*bx(i+1,j,k)
     *                       - psi(i-1,j,k)**4*bx(i-1,j,k) )/(2.d0*h)
     *               -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Ad22 = 0.5d0*( 
     *                2.d0*(   psi(i,j+1,k)**4*by(i,j+1,k)
     *                       - psi(i,j-1,k)**4*by(i,j-1,k) )/(2.d0*h)
     *               -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Ad33 = 0.5d0*( 
     *                2.d0*(   psi(i,j,k+1)**4*bz(i,j,k+1)
     *                       - psi(i,j,k-1)**4*bz(i,j,k-1) )/(2.d0*h)
     *               -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  !
                  Ad12 = 0.5d0*( 
     *                     (   psi(i+1,j,k)**4*by(i+1,j,k)
     *                       - psi(i-1,j,k)**4*by(i-1,j,k) )/(2.d0*h)
     *                   + (   psi(i,j+1,k)**4*bx(i,j+1,k)
     *                       - psi(i,j-1,k)**4*bx(i,j-1,k) )/(2.d0*h)
     *                          )
                  Ad13 = 0.5d0*( 
     *                     (   psi(i+1,j,k)**4*bz(i+1,j,k)
     *                       - psi(i-1,j,k)**4*bz(i-1,j,k) )/(2.d0*h)
     *                   + (   psi(i,j,k+1)**4*bx(i,j,k+1)
     *                       - psi(i,j,k-1)**4*bx(i,j,k-1) )/(2.d0*h)
     *                          )
                  Ad23 = 0.5d0*( 
     *                     (   psi(i,j+1,k)**4*bz(i,j+1,k)
     *                       - psi(i,j-1,k)**4*bz(i,j-1,k) )/(2.d0*h)
     *                   + (   psi(i,j,k+1)**4*by(i,j,k+1)
     *                       - psi(i,j,k-1)**4*by(i,j,k-1) )/(2.d0*h)
     *                          )
                  !
                  Au11 = 0.5d0*(  
     *                        2.d0*(bx(i+1,j,k)-bx(i-1,j,k))/(2.d0*h)
     *                             / psi(i,j,k)**4
     *                      -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Au22 = 0.5d0*(  
     *                        2.d0*(by(i,j+1,k)-by(i,j-1,k))/(2.d0*h)
     *                             / psi(i,j,k)**4
     *                      -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Au33 = 0.5d0*(  
     *                        2.d0*(bz(i,j,k+1)-bz(i,j,k-1))/(2.d0*h)
     *                             / psi(i,j,k)**4
     *                      -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  !
                  Au12 = 0.5d0*(   (by(i+1,j,k)-by(i-1,j,k))/(2.d0*h)
     *                           + (bx(i,j+1,k)-bx(i,j-1,k))/(2.d0*h)
     *                          ) / psi(i,j,k)**4
                  Au13 = 0.5d0*(   (bz(i+1,j,k)-bz(i-1,j,k))/(2.d0*h)
     *                           + (bx(i,j,k+1)-bx(i,j,k-1))/(2.d0*h)
     *                          ) / psi(i,j,k)**4
                  Au23 = 0.5d0*(   (bz(i,j+1,k)-bz(i,j-1,k))/(2.d0*h)
     *                           + (by(i,j,k+1)-by(i,j,k-1))/(2.d0*h)
     *                          ) / psi(i,j,k)**4
                  !
                  !---------------------------------
                  !
                  !    Psi Equation:
                  !
                  !
                  lhs_psi  = (    d2psi_dx2 + d2psi_dy2 + d2psi_dz2
     *                         +  (1.d0 / 8.d0 / psi(i,j,k)**7)
     *                           *(   Ad11*Au11 + Ad12*Au12 + Ad13*Au13
     *                              + Ad12*Au12 + Ad22*Au22 + Ad23*Au23
     *                              + Ad13*Au13 + Ad23*Au23 + Ad33*Au33)
     *                         ) / psi(i,j,k)**5
                  !
                  !---------------------------------
                  !
                  !    Bx Equation:
                  !
                  lhs_bx  = (     d2bx_dx2 + d2bx_dy2 + d2bx_dz2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *(d2bx_dx2 + d2by_dxdy + d2bz_dxdz )
     *                         ) / psi(i,j,k)**(10)
                  !
                  !---------------------------------
                  !
                  !    By Equation:
                  !
                  lhs_by  = (     d2by_dx2 + d2by_dy2 + d2by_dz2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *(d2bx_dxdy + d2by_dy2 + d2bz_dydz )
     *                         ) / psi(i,j,k)**(10)
                  !
                  !---------------------------------
                  !
                  !    Bz Equation:
                  !
                  lhs_bz  = (     d2bz_dx2 + d2bz_dy2 + d2bz_dz2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *(d2bx_dxdz + d2by_dydz + d2bz_dz2 )
     *                         ) / psi(i,j,k)**(10)
               !
               result_psi(i,j,k) = result_psi(i,j,k) + lhs_psi
               result_bx(i,j,k)  = result_bx(i,j,k)  + lhs_bx
               result_by(i,j,k)  = result_by(i,j,k)  + lhs_by
               result_bz(i,j,k)  = result_bz(i,j,k)  + lhs_bz
               !
            end do
         end do
      end do

      !
      ! Assuming dirichlet for now:
      !
      call load_scal_bound3d( result_psi, 0.d0, nx,ny,nz)
      call load_scal_bound3d( result_bx,  0.d0, nx,ny,nz)
      call load_scal_bound3d( result_by,  0.d0, nx,ny,nz)
      call load_scal_bound3d( result_bz,  0.d0, nx,ny,nz)

      if (ltrace) then
         write(*,*) 'lop: '
      end if

      return
      end      ! END: lop

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  relax:     See Carlos' notes                                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function relax( psi, rhs_psi,
     *           bx, rhs_bx, by, rhs_by, bz, rhs_bz,
     *           chr, h,nx,ny,nz)
      implicit      none
      integer       nx, ny, nz
      real(kind=8)  h
      real(kind=8)  psi(nx,ny,nz),   rhs_psi(nx,ny,nz),
     *              bx(nx,ny,nz),    rhs_bx(nx,ny,nz),
     *              by(nx,ny,nz),    rhs_by(nx,ny,nz),
     *              bz(nx,ny,nz),    rhs_bz(nx,ny,nz),
     *              chr(nx,ny,nz)
      include      'chr.inc'
      include      'largesmall.inc'
      integer       pass, ib, kswitch
      integer       i, j, k, m
      integer       myid, proc_return_myid
      external            proc_return_myid
      real(kind=8)  lhs_psi, lhs_bx, lhs_by, lhs_bz
      real(kind=8)  jacobian(4,4), res(4), resnrms(4)
      !
      real(kind=8)  d2psi_dx2, d2psi_dy2, d2psi_dz2,
     *              d2bx_dx2,  d2bx_dy2,  d2bx_dz2,
     *              d2by_dx2,  d2by_dy2,  d2by_dz2,
     *              d2bz_dx2,  d2bz_dy2,  d2bz_dz2,
     *              Ad11, Ad12, Ad13, Au11, Au12, Au13,
     *                    Ad22, Ad23,       Au22, Au23,
     *                          Ad33,             Au33,
     *              d2bx_dxdy,            d2bx_dxdz,
     *              d2by_dxdy, d2by_dydz,
     *                         d2bz_dydz, d2bz_dxdz,
     *              sum_dbk_dk
      !
      logical     ltrace
      parameter ( ltrace = .false. )


      if (ltrace) then
         write(*,*) 'relax: Relaxing'
      end if

      if (h .eq. 0) then
         write(*,*) 'Problem with h = ',h
         stop
      end if

      ! Initialize norms to zero:
      do m = 1, 4
         resnrms(m) = 0.d0
      end do

      !
      ! Red-Black ordering of the interior only:
      !
      do pass = 1, 2
         kswitch = pass
         do k = 2, nz-1
            ib = kswitch + 1
            do j = 2, ny-1
               do i = ib, nx-1, 2
                  !
                  d2psi_dx2 =
     *              (psi(i+1,j,k) - 2.d0*psi(i,j,k) + psi(i-1,j,k))/h**2
                  d2psi_dy2 =
     *              (psi(i,j+1,k) - 2.d0*psi(i,j,k) + psi(i,j-1,k))/h**2
                  d2psi_dz2 =
     *              (psi(i,j,k+1) - 2.d0*psi(i,j,k) + psi(i,j,k-1))/h**2
                  !
                  d2bx_dx2  =
     *              (bx(i+1,j,k) - 2.d0*bx(i,j,k) + bx(i-1,j,k))/h**2
                  d2bx_dy2  =
     *              (bx(i,j+1,k) - 2.d0*bx(i,j,k) + bx(i,j-1,k))/h**2
                  d2bx_dz2  =
     *              (bx(i,j,k+1) - 2.d0*bx(i,j,k) + bx(i,j,k-1))/h**2
                  !
                  d2by_dx2  =
     *              (by(i+1,j,k) - 2.d0*by(i,j,k) + by(i-1,j,k))/h**2
                  d2by_dy2  =
     *              (by(i,j+1,k) - 2.d0*by(i,j,k) + by(i,j-1,k))/h**2
                  d2by_dz2  =
     *              (by(i,j,k+1) - 2.d0*by(i,j,k) + by(i,j,k-1))/h**2
                  !
                  d2bz_dx2  =
     *              (bz(i+1,j,k) - 2.d0*bz(i,j,k) + bz(i-1,j,k))/h**2
                  d2bz_dy2  =
     *              (bz(i,j+1,k) - 2.d0*bz(i,j,k) + bz(i,j-1,k))/h**2
                  d2bz_dz2  =
     *              (bz(i,j,k+1) - 2.d0*bz(i,j,k) + bz(i,j,k-1))/h**2
                  !
                  d2bx_dxdy =  (  (bx(i+1,j+1,k) - bx(i+1,j-1,k))
     *                           -(bx(i-1,j+1,k) - bx(i-1,j-1,k))
     *                           ) / (4.d0*h**2)
                  d2bx_dxdz =  (  (bx(i+1,j,k+1) - bx(i+1,j,k-1))
     *                           -(bx(i-1,j,k+1) - bx(i-1,j,k-1))
     *                           ) / (4.d0*h**2)
                  !
                  d2by_dxdy =  (  (by(i+1,j+1,k) - by(i+1,j-1,k))
     *                           -(by(i-1,j+1,k) - by(i-1,j-1,k))
     *                           ) / (4.d0*h**2)
                  d2by_dydz =  (  (by(i,j+1,k+1) - by(i,j+1,k-1))
     *                           -(by(i,j-1,k+1) - by(i,j-1,k-1))
     *                           ) / (4.d0*h**2)
                  !
                  d2bz_dxdz =  (  (bz(i+1,j,k+1) - bz(i+1,j,k-1))
     *                           -(bz(i-1,j,k+1) - bz(i-1,j,k-1))
     *                           ) / (4.d0*h**2)
                  d2bz_dydz =  (  (bz(i,j+1,k+1) - bz(i,j+1,k-1))
     *                           -(bz(i,j-1,k+1) - bz(i,j-1,k-1))
     *                           ) / (4.d0*h**2)
                  !
                  sum_dbk_dk = (   bx(i+1,j,k)-bx(i-1,j,k)
     *                           + by(i,j+1,k)-by(i,j-1,k)
     *                           + by(i,j,k+1)-by(i,j,k-1)
     *                           ) / (2.d0*h)
                  !
                  Ad11 = 0.5d0*( 
     *                2.d0*(   psi(i+1,j,k)**4*bx(i+1,j,k)
     *                       - psi(i-1,j,k)**4*bx(i-1,j,k) )/(2.d0*h)
     *               -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Ad22 = 0.5d0*( 
     *                2.d0*(   psi(i,j+1,k)**4*by(i,j+1,k)
     *                       - psi(i,j-1,k)**4*by(i,j-1,k) )/(2.d0*h)
     *               -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Ad33 = 0.5d0*( 
     *                2.d0*(   psi(i,j,k+1)**4*bz(i,j,k+1)
     *                       - psi(i,j,k-1)**4*bz(i,j,k-1) )/(2.d0*h)
     *               -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  !
                  Ad12 = 0.5d0*( 
     *                     (   psi(i+1,j,k)**4*by(i+1,j,k)
     *                       - psi(i-1,j,k)**4*by(i-1,j,k) )/(2.d0*h)
     *                   + (   psi(i,j+1,k)**4*bx(i,j+1,k)
     *                       - psi(i,j-1,k)**4*bx(i,j-1,k) )/(2.d0*h)
     *                          )
                  Ad13 = 0.5d0*( 
     *                     (   psi(i+1,j,k)**4*bz(i+1,j,k)
     *                       - psi(i-1,j,k)**4*bz(i-1,j,k) )/(2.d0*h)
     *                   + (   psi(i,j,k+1)**4*bx(i,j,k+1)
     *                       - psi(i,j,k-1)**4*bx(i,j,k-1) )/(2.d0*h)
     *                          )
                  Ad23 = 0.5d0*( 
     *                     (   psi(i,j+1,k)**4*bz(i,j+1,k)
     *                       - psi(i,j-1,k)**4*bz(i,j-1,k) )/(2.d0*h)
     *                   + (   psi(i,j,k+1)**4*by(i,j,k+1)
     *                       - psi(i,j,k-1)**4*by(i,j,k-1) )/(2.d0*h)
     *                          )
                  !
                  Au11 = 0.5d0*(  
     *                        2.d0*(bx(i+1,j,k)-bx(i-1,j,k))/(2.d0*h)
     *                             / psi(i,j,k)**4
     *                      -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Au22 = 0.5d0*(  
     *                        2.d0*(by(i,j+1,k)-by(i,j-1,k))/(2.d0*h)
     *                             / psi(i,j,k)**4
     *                      -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  Au33 = 0.5d0*(  
     *                        2.d0*(bz(i,j,k+1)-bz(i,j,k-1))/(2.d0*h)
     *                             / psi(i,j,k)**4
     *                      -(2.d0/3.d0)* sum_dbk_dk
     *                          )
                  !
                  Au12 = 0.5d0*(   (by(i+1,j,k)-by(i-1,j,k))/(2.d0*h)
     *                           + (bx(i,j+1,k)-bx(i,j-1,k))/(2.d0*h)
     *                          ) / psi(i,j,k)**4
                  Au13 = 0.5d0*(   (bz(i+1,j,k)-bz(i-1,j,k))/(2.d0*h)
     *                           + (bx(i,j,k+1)-bx(i,j,k-1))/(2.d0*h)
     *                          ) / psi(i,j,k)**4
                  Au23 = 0.5d0*(   (bz(i,j+1,k)-bz(i,j-1,k))/(2.d0*h)
     *                           + (by(i,j,k+1)-by(i,j,k-1))/(2.d0*h)
     *                          ) / psi(i,j,k)**4
                  !
                  !---------------------------------
                  !
                  !    Psi Equation:
                  !
                  !
                  lhs_psi  = (    d2psi_dx2 + d2psi_dy2 + d2psi_dz2
     *                         +  (1.d0 / 8.d0 / psi(i,j,k)**7)
     *                           *(   Ad11*Au11 + Ad12*Au12 + Ad13*Au13
     *                              + Ad12*Au12 + Ad22*Au22 + Ad23*Au23
     *                              + Ad13*Au13 + Ad23*Au23 + Ad33*Au33)
     *                         ) / psi(i,j,k)**5
                  res(1)   = lhs_psi - rhs_psi(i,j,k)
                  !
                  ! partial Res_psi / partial psi:
                  jacobian(1,1) = -(5.d0 / psi(i,j,k)) * lhs_psi
     *                       +(  -6.d0/h**2
     *                           -(7.d0 / 8.d0 / psi(i,j,k)**8)
     *                             *( Ad11*Au11 + Ad12*Au12 + Ad13*Au13
     *                              + Ad12*Au12 + Ad22*Au22 + Ad23*Au23
     *                              + Ad13*Au13 + Ad23*Au23 + Ad33*Au33)
     *                           -(4.d0 / 8.d0 / psi(i,j,k)**8)
     *                           *(   Ad11*(Au11+(1.d0/3.d0)*sum_dbk_dk)
     *                                          + Ad12*Au12 + Ad13*Au13
     *                              + Ad22*(Au22+(1.d0/3.d0)*sum_dbk_dk)
     *                                          + Ad12*Au12 + Ad23*Au23
     *                              + Ad33*(Au33+(1.d0/3.d0)*sum_dbk_dk)
     *                                          + Ad13*Au13 + Ad23*Au23)
     *                         ) / psi(i,j,k)**5
                  ! partial Res_psi / partial bx:
                  jacobian(1,2) = 0.d0
                  ! partial Res_psi / partial by:
                  jacobian(1,3) = 0.d0
                  ! partial Res_psi / partial bz:
                  jacobian(1,4) = 0.d0
                  !
                  !---------------------------------
                  !
                  !    Bx Equation:
                  !
                  lhs_bx  = (     d2bx_dx2 + d2bx_dy2 + d2bx_dz2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *(d2bx_dx2 + d2by_dxdy + d2bz_dxdz )
     *                         ) / psi(i,j,k)**(10)
                  res(2)  = lhs_bx - rhs_bx(i,j,k)
                  !
                  ! partial Res_bx / partial psi:
                  jacobian(2,1) = -(10.d0 / psi(i,j,k)) * lhs_bx
     *                           - (4.d0 / 3.d0 / psi(i,j,k)**5)
     *                           *(d2bx_dx2 + d2by_dxdy + d2bz_dxdz )
     *                             / psi(i,j,k)**(10)
                  ! partial Res_bx / partial bx:
                  jacobian(2,2) = ( -6.d0/h**2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *( -2.d0/h**2 )
     *                         ) / psi(i,j,k)**(10)
                  ! partial Res_bx / partial by:
                  jacobian(2,3) = 0.d0
                  ! partial Res_bx / partial bz:
                  jacobian(2,4) = 0.d0
                  !
                  !---------------------------------
                  !
                  !    By Equation:
                  !
                  lhs_by  = (     d2by_dx2 + d2by_dy2 + d2by_dz2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *(d2bx_dxdy + d2by_dy2 + d2bz_dydz )
     *                         ) / psi(i,j,k)**(10)
                  res(3)  = lhs_by - rhs_by(i,j,k)
                  !
                  ! partial Res_by / partial psi:
                  jacobian(3,1) = -(10.d0 / psi(i,j,k)) * lhs_by
     *                           - (4.d0 / 3.d0 / psi(i,j,k)**5)
     *                           *(d2bx_dxdy + d2by_dy2 + d2bz_dydz )
     *                             / psi(i,j,k)**(10)
                  ! partial Res_by / partial bx:
                  jacobian(3,2) = 0.d0
                  ! partial Res_by / partial by:
                  jacobian(3,3) = ( -6.d0/h**2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *( -2.d0/h**2 )
     *                         ) / psi(i,j,k)**(10)
                  ! partial Res_by / partial bz:
                  jacobian(3,4) = 0.d0
                  !
                  !---------------------------------
                  !
                  !    Bz Equation:
                  !
                  lhs_bz  = (     d2bz_dx2 + d2bz_dy2 + d2bz_dz2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *(d2bx_dxdz + d2by_dydz + d2bz_dz2 )
     *                         ) / psi(i,j,k)**(10)
                  res(4)  = lhs_bz - rhs_bz(i,j,k)
                  !
                  ! partial Res_bz / partial psi:
                  jacobian(4,1) = -(10.d0 / psi(i,j,k)) * lhs_bz
     *                           - (4.d0 / 3.d0 / psi(i,j,k)**5)
     *                           *(d2bx_dxdz + d2by_dydz + d2bz_dz2 )
     *                             / psi(i,j,k)**(10)
                  ! partial Res_bz / partial bx:
                  jacobian(4,2) = 0.d0
                  ! partial Res_bz / partial by:
                  jacobian(4,3) = 0.d0
                  ! partial Res_bz / partial bz:
                  jacobian(4,4) = ( -6.d0/h**2
     *                         +  (1.d0 / 3.d0 / psi(i,j,k)**4)
     *                           *( -2.d0/h**2 )
     *                         ) / psi(i,j,k)**(10)
                  !
                  !---------------------------------
                  !
                  !   Solve:
                  !
                  call solver4by4_nopivot(jacobian,res)
                  !
                  !   Update fields:
                  !
                  psi(i,j,k) = psi(i,j,k) - res(1)
                  bx(i,j,k)  = bx(i,j,k)  - res(2)
                  by(i,j,k)  = by(i,j,k)  - res(3)
                  bz(i,j,k)  = bz(i,j,k)  - res(4)
                  !
                  do m = 1, 4
                     resnrms(m) = resnrms(m) + res(m)**2
                  end do
                  !
                  !
                  !
               end do
               ib = 5 - ib
            end do
            kswitch = 3 - kswitch
         end do
      end do

      do m = 1, 4
         relax = relax + resnrms(m)
      end do
      relax = sqrt( relax /nx /ny /nz / 4.d0 )

      if (ltrace) then
         write(*,*) 'relax: ',relax
      end if

  95  format('[',I3,'] ',A,2I5,3F18.10)
  96  format('[',I3,'] ',A,I5,3F18.10)
  98  format('[',I3,'] ',A,3F18.10)
  99  format('[',I3,'] ',A,3I5)
      return
      end      ! END: relax

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ........................................................
c
c     Inline 4x4 solver: (NO PIVOTING)
c
c     ........................................................
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine solver4by4_nopivot(a,b)
      implicit none
      real*8  a(4,4),b(4)

      real*8      dum
      real*8      TINY
      parameter ( TINY   = 1.0d-20 )

      if (a(1,1).eq.0) then
          write(*,*) 'solve4by4: a(1,1) is small=', a(1,1)
          a(1,1) = TINY
      end if
      dum     = 1.d0 / a(1,1)
      a(2,1)  = a(2,1) * dum
      a(3,1)  = a(3,1) * dum
      a(4,1)  = a(4,1) * dum
c
c     j=2:
c
c             i=2
      a(2,2)  = a(2,2) - a(2,1)*a(1,2)
c             i=3
      a(3,2)  = a(3,2) - a(3,1)*a(1,2)
c             i=4
      a(4,2)  = a(4,2) - a(4,1)*a(1,2)
      if (a(2,2).eq.0) then
          write(*,*) 'solve4by4: a(2,2) is small=', a(2,2)
          a(2,2) = TINY
      end if
      dum     = 1.d0 / a(2,2)
      a(3,2)  = a(3,2) * dum
      a(4,2)  = a(4,2) * dum
c
c     j=3
c
      a(2,3)  = a(2,3) - a(2,1)*a(1,3)
c             i=3
      a(3,3)  = a(3,3) - a(3,1)*a(1,3) - a(3,2)*a(2,3)
c             i=4
      a(4,3)  = a(4,3) - a(4,1)*a(1,3) - a(4,2)*a(2,3)
      if (a(3,3).eq.0) then
          write(*,*) 'solve4by4: a(3,3) is small=', a(3,3)
          if (a(3,3).gt.0) write(*,*) '>0'
          if (a(3,3).le.0) write(*,*) '<=0'
          a(3,3) = TINY
      end if
      a(4,3)  = a(4,3) / a(3,3)
c
c     j=4:
c
      a(2,4) = a(2,4) - a(2,1)*a(1,4)
      a(3,4) = a(3,4) - a(3,1)*a(1,4) - a(3,2)*a(2,4)
c             i=4
      a(4,4)  = a(4,4) - a(4,1)*a(1,4) - a(4,2)*a(2,4) -a(4,3)*a(3,4)
      if (a(4,4).eq.0) then
          write(*,*) 'solve4by4: a(4,4) is small=', a(4,4)
          a(4,4) = TINY
      end if



c
c     Forward substitution:
c
c          i=1:
c     b(1)  = b(1)
c          i=2:
      b(2) = b(2) - a(2,1)*b(1)
c          i=3:
      b(3) = b(3) - a(3,1)*b(1) - a(3,2)*b(2)
c          i=4:
      b(4) = b(4) - a(4,1)*b(1) - a(4,2)*b(2) - a(4,3)*b(3)

c
c     Backward substitution:
c
c     i=4:
      b(4) = (b(4) ) / a(4,4)
c     i=3:
      b(3) = (b(3) - a(3,4)*b(4) ) / a(3,3)
c     i=2:
      b(2) = (b(2) - a(2,3)*b(3) - a(2,4)*b(4) ) / a(2,2)
c     i=1:
      b(1) = (b(1) - a(1,2)*b(2) - a(1,3)*b(3) - a(1,4)*b(4) ) / a(1,1)

      return
      end
