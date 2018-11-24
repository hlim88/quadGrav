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
!   This solves for the correction u() to the conformal factor psi(x,y,z)
!   for BBH initial data, using the elliptic Hamiltonian constraint  
!   equation.
!
!      As described in http://arxiv.org/abs/gr-qc/9703066                        Eq.  9
!                  and http://relativity.livingreviews.org/Articles/lrr-2000-5/  Eq. 77
!
!                     \nabla^2 u + beta (1+alpha u)^(-7) = 0
!
!      e^(4*phi) = psi^4 and chi = exp(-4 phi ), then
!
!      chi = [e^(4*phi)]^-1 = [psi^4]^-1 = 1/(psi^4)
!
!    ! Need to calculate: tildeA_ij tildeA^ij
!    !
!    ! Using:
!    !            g^{ik} = psi^{-4} delta^{ij}
!    !   and:        chi = psi^(-4)
!    !   and:  tildeA_ij = psi^2*A_ij 
!    !   and:       A_ij = stored grid functions in the code
!    !
!    ! Calculate: tildeA_ij tildeA^ij =
!    !  = g^{ik}g^{jl} tilde A_{ij} tilde A_{kl} 
!    !  = psi^{-8} tilde[A11**2+2*A12**2+2*A13**2+A22**2++2*A23**2+A33**2]
!    !  = psi^{-4}      [A11**2+2*A12**2+2*A13**2+A22**2++2*A23**2+A33**2]
!    !  = chi           [A11**2+2*A12**2+2*A13**2+A22**2++2*A23**2+A33**2]
!    ! These are the non-tilde quantities
!    ! So ultimately we have:
!    ! tilde A_{ij} tilde A^{ij} =chi  [A11**2+A22**2+A33**2+2*(A12**2+A13**2+A23**2)]
!    !       A_{ij}       A^{ij} =chi^2[A11**2+A22**2+A33**2+2*(A12**2+A13**2+A23**2)]
!
!   Note:  First implementation used "u" as defined in the papers.
!          However, if "useuminusone" is set to true in various places
!          in this file, then the problem is solved in terms of epsilon=1-u
!          You also have to change how this correction is set in initial.f90
!
!   Note:  If you change anything in relax(), lop(), or get_resnrm(),
!          you *must* make identical to changes to the other two routines
!          of this triplet or else things will not converge.
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
      !
      integer     nx,ny,nz, lev,myid,myint
      integer     i
      integer     proc_return_myid
      external    proc_return_myid
      real*8      hg, time
      character(2) tmps
      logical     ltrace
      parameter ( ltrace = .false. )


      myid = proc_return_myid()
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

      if(ltrace)write(*,94)myid,'grid_mg_trace:Output fields:',tmps,time
      do i = 1, num_gfuncs
         if (gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL) then
            ! This is an elliptic field:
            call field_out3d(q(gfunc_pointer(i,gi)),time,
     *             'mg_'//gfunc_name(i)//tmps,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
            ! This should be its RHS:
            call field_out3d(q(gfunc_pointer(i+3,gi)),time,
     *             'mg_'//gfunc_name(i+3)//tmps,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
            ! This should be its residual:
            call field_out3d(q(gfunc_pointer(i+1,gi)),time,
     *             'mg_'//gfunc_name(i+1)//tmps,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             nx,ny,nz,myid)
         end if
      end do

  94  format('[',I4,'] ',A,A8,3F18.10)
  95  format('[',I4,'] ',A,2I5,3F18.10)
  96  format('[',I4,'] ',A,I5,3F18.10)
  98  format('[',I4,'] ',A,3F18.10)
  99  format('[',I4,'] ',A,3I5)

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
      integer    myid, proc_return_myid
      external         proc_return_myid

      myid = proc_return_myid()
      if (ltrace) write(*,99)myid,'grid_get_resid: On grid: ',gi

      call load_pointers(gi)

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      call get_resid(q(gr_uell_st1),q(gr_uell), q(gr_uell_rhs),
     .           q(gr_chr),
     .           q(gr_A11), q(gr_A12), q(gr_A13),
     .           q(gr_A22), q(gr_A23), q(gr_A33),
     .           q(gr_chi),
     .           q(gr_x(gi)),q(gr_y(gi)),q(gr_z(gi)),
     .           hg, nx, ny, nz)
c     call get_resid(q(gr_psi_elliptic_st1),
c    *               q(gr_psi_elliptic),
c    *               q(gr_psi_elliptic_rhs),hg,nx,ny,nz)

      if (ltrace) then
         call field_out3d(q(gr_uell_st1), 0.d0,'grid_grst1',
     *             q(gr_x(gi)),   q(gr_x(gi)+nx-1),
     *             q(gr_y(gi)),   q(gr_y(gi)+ny-1),
     *             q(gr_z(gi)),   q(gr_z(gi)+nz-1),
     *             nx,   ny,    nz, myid)
      end if

  94  format('[',I4,'] ',A,A8,3F18.10)
  95  format('[',I4,'] ',A,2I5,3F18.10)
  96  format('[',I4,'] ',A,I5,3F18.10)
  98  format('[',I4,'] ',A,3F18.10)
  99  format('[',I4,'] ',A,3I5)

      return
      end    ! END: grid_get_resid

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_lop:                                                                 cc
cc            Compute the results of the differential operator associated     cc
cc            with the elliptic problem acting on the set of restricted       cc
cc            fields residing in the u_st1 fields and *add* the result        cc
cc            to the u_st2 arrays. This operation occurs on coarsened data    cc
cc            as part of the CGC (coarse grid correction).                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_lop( gi )
      implicit    none
      integer     gi
      include    'grid.inc'
      include    'param.inc'
      integer     nx,  ny,  nz
      integer     nxc, nyc, nzc
      real*8      hg, hc
      integer    i,j,k,index
      integer    myid, proc_return_myid
      external         proc_return_myid
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      myid = proc_return_myid()
      call load_pointers(gi)

  95  format('[',I4,'] ',A,2I5,3F18.10)
  96  format('[',I4,'] ',A,I5,3F18.10)
  98  format('[',I4,'] ',A,5F18.10)
  99  format('[',I4,'] ',A,3I5)
      if(ltrace)write(*,99)myid,'grid_lop:    On grid: ',gi

      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)
      hg   = gr_h(gi)

      nxc = NINT( 0.5d0*(nx+1) )
      nyc = NINT( 0.5d0*(ny+1) )
      nzc = NINT( 0.5d0*(nz+1) )
      hc  = 2.d0*hg

      ! Need coarsened version of coordinates:
      if(ltrace)write(*,99)myid,'grid_lop:    Creating coarsened coords'
      if(ltrace)write(*,99)myid,'grid_lop: nx/y/z: ',nxc,nyc,nzc
      do i = 1,nxc
         q(gr_tmp+i-1)          = q(gr_x(gi)+2*(i-1))
      end do
      do j = 1,nyc
         q(gr_tmp+nxc+j-1)      = q(gr_y(gi)+2*(j-1))
      end do
      do k = 1,nzc
         q(gr_tmp+nxc+nyc+k-1)  = q(gr_z(gi)+2*(k-1))
      end do
      ! Need coarsened versions of A?? fields:
      call restrict( q(gr_A11_rk1),  q(gr_A11),  nxc,nyc,nzc)
      call restrict( q(gr_A12_rk1),  q(gr_A12),  nxc,nyc,nzc)
      call restrict( q(gr_A13_rk1),  q(gr_A13),  nxc,nyc,nzc)
      call restrict( q(gr_A22_rk1),  q(gr_A22),  nxc,nyc,nzc)
      call restrict( q(gr_A23_rk1),  q(gr_A23),  nxc,nyc,nzc)
      call restrict( q(gr_A33_rk1),  q(gr_A33),  nxc,nyc,nzc)
      call restrict( q(gr_chi_rk1),  q(gr_chi),  nxc,nyc,nzc)

      if (ltrace) then
         write(*,99)myid,'grid_lop:   Outputting coarsened coords'
         call field_out1d(q(gr_tmp),        0.d0,'lopx',gr_minx(gi),
     .         gr_maxx(gi), nxc, myid)
         call field_out1d(q(gr_tmp+nxc),    0.d0,'lopy',gr_miny(gi),
     .         gr_maxy(gi), nyc, myid)
         call field_out1d(q(gr_tmp+nxc+nyc),0.d0,'lopz',gr_minz(gi),
     .         gr_maxz(gi), nzc, myid)
         call field_out3d(q(gr_uell_st2), 0.d0,'grid_lopST2pre',
     *             q(gr_tmp),         q(gr_tmp+nxc-1),
     *             q(gr_tmp+nxc),     q(gr_tmp+nxc+nyc-1),
     *             q(gr_tmp+nxc+nyc), q(gr_tmp+nxc+nyc+nzc-1),
     *             nxc,   nyc,    nzc, myid)
      end if

      !
      ! Computing on a coarsened version of the fields here
      ! so be careful to use coarsened values of h and nx, etc:
      !
      call lop(q(gr_uell_st2),q(gr_uell_st1),
!    .              q(gr_A11), q(gr_A12),     q(gr_A13),
!    .              q(gr_A22), q(gr_A23),     q(gr_A33),
!    .              q(gr_chi),
     .              q(gr_A11_rk1), q(gr_A12_rk1), q(gr_A13_rk1),
     .              q(gr_A22_rk1), q(gr_A23_rk1), q(gr_A33_rk1),
     .              q(gr_chi_rk1),
     .              q(gr_tmp),q(gr_tmp+nxc),q(gr_tmp+nxc+nyc),hc,
     .              bh_1_mass, bh_1_x, bh_1_y, bh_1_z,
     .              bh_2_mass, bh_2_x, bh_2_y, bh_2_z,
     .              bh_n, nxc,nyc,nzc)
c     call lop( q(gr_psi_elliptic_st2),
c    *          q(gr_psi_elliptic_st1),hc,nxc,nyc,nzc)

      if (ltrace) then
         write(*,99)myid,'grid_lop:   Done.'
         call field_out3d(q(gr_uell_st1), 0.d0,'grid_lopST1',
     *             q(gr_tmp),         q(gr_tmp+nxc-1),
     *             q(gr_tmp+nxc),     q(gr_tmp+nxc+nyc-1),
     *             q(gr_tmp+nxc+nyc), q(gr_tmp+nxc+nyc+nzc-1),
     *             nxc,   nyc,    nzc, myid)
         call field_out3d(q(gr_uell_st2), 0.d0,'grid_lopST2',
     *             q(gr_tmp),         q(gr_tmp+nxc-1),
     *             q(gr_tmp+nxc),     q(gr_tmp+nxc+nyc-1),
     *             q(gr_tmp+nxc+nyc), q(gr_tmp+nxc+nyc+nzc-1),
     *             nxc,   nyc,    nzc, myid)
      end if

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
      include    'glob.inc'
      include    'param.inc'
      integer     nx, ny, nz
      real*8      hg
      !
      ! The test one can do is to simply choose what uell() should be,
      ! in this case a gaussian pulse, and then set the RHS of the equation
      ! to whatever function it would have to have so that a Gaussian is
      ! indeed a solution. If the elliptic solver finds the that same Gaussian,
      ! then the solver is working correctly.
      !
      logical     dotest
      parameter ( dotest = .false. )
      real(kind=8)  minradius
      parameter (   minradius  = 0.001d0)
      logical     ltrace
      parameter ( ltrace = .false. )
      integer    i,j,k,index
      integer    indexip1, indexim1
      integer    indexjp1, indexjm1
      integer    indexkp1, indexkm1
      real*8     invalpha, alpha,   beta, u_pt, chi_pt
      real*8     z_pt,z1,z2
      real*8     y_pt,y1,y2
      real*8     x_pt,x1,x2
      real*8     invh,invhsq,rv1,rv2, kabkab
      real*8     du_dxx, du_dyy, du_dzz


      call load_pointers(gridnum)

      if (ltrace) write(*,*) 'grid_ellrhs:    On grid: ',gridnum

      nx   = gr_nx(gridnum)
      ny   = gr_ny(gridnum)
      nz   = gr_nz(gridnum)
      hg   = gr_h(gridnum)

      !
      ! Compute interior:
      !
      if (dotest) then
         write(*,*) 'grid_ellrhs: TESTING: choosing Gaussian for uell'
         if (.true.) then
            call gaussian3d(q(gr_tmp), hg,  q(gr_x(gridnum)),
     .                   q(gr_y(gridnum)),q(gr_z(gridnum)),
     .                   nx,ny,nz,
!                        amp, ellipticity_x, ellipticity_y,...
     .                   5.3d0, 1.4d0, 0.6d0,0.d0,2.1d0,0.d0,0.d0,0.d0)
         else
            do k=1,nz
            do j=1,ny
            do i=1,nx
                  index    = (k-1)*ny*nx + (j-1)*nx + (i-1)
                  x_pt     = q(gr_x(gridnum)+i-1)
                  y_pt     = q(gr_y(gridnum)+j-1)
                  z_pt     = q(gr_z(gridnum)+k-1)
                  rv1      = sqrt(x_pt**2+y_pt**2+z_pt**2)
                  if (rv1.lt.minradius) rv1 = minradius
                  !q(gr_tmp+index) = 1.d0 / rv1
                  q(gr_tmp+index) = minradius / rv1
            end do
            end do
            end do
         end if
         write(*,*) 'grid_ellrhs: TESTING: output as dotest_uell.sdf'
         call field_out3d(q(gr_tmp), 0.d0, 'dotest_uell',
     *             q(gr_x(gridnum)), q(gr_x(gridnum)+nx-1),
     *             q(gr_y(gridnum)), q(gr_y(gridnum)+ny-1),
     *             q(gr_z(gridnum)), q(gr_z(gridnum)+nz-1),
     *             nx,   ny,    nz, 0d0)
         write(*,*) 'grid_ellrhs: TESTING: compared uell to this'
      end if
      !
      do k = 1, nz
         z_pt     = q(gr_z(gridnum)+k-1)
         z1       = z_pt - bh_1_z
         z2       = z_pt - bh_2_z
         do j = 1, ny
            y_pt     = q(gr_y(gridnum)+j-1)
            y1       = y_pt - bh_1_y
            y2       = y_pt - bh_2_y
            do i = 1, nx
               !
               !
               invh     = 1.d0 / hg
               invhsq   = invh / hg
               x_pt     = q(gr_x(gridnum)+i-1)
               x1       = x_pt - bh_1_x
               x2       = x_pt - bh_2_x
               !
               rv1      = sqrt(x1**2+y1**2+z1**2)
               rv2      = sqrt(x2**2+y2**2+z2**2)
               !
               indexip1 = (k-1)*ny*nx + (j-1)*nx + (i  )
               indexim1 = (k-1)*ny*nx + (j-1)*nx + (i-2)
               !
               indexjp1 = (k-1)*ny*nx + (j  )*nx + (i-1)
               indexjm1 = (k-1)*ny*nx + (j-2)*nx + (i-1)
               !
               indexkp1 = (k  )*ny*nx + (j-1)*nx + (i-1)
               indexkm1 = (k-2)*ny*nx + (j-1)*nx + (i-1)
               !
               index    = (k-1)*ny*nx + (j-1)*nx + (i-1)
               !
               if (.not. dotest) then
                  q(gr_uell_rhs+index) = 0.d0
               else
                  !
                  ! Substitute the tmp array into the
                  ! equation to be solved with all terms
                  ! on the left-hand side:
                  !
                  u_pt     = q(gr_tmp+index)
                  !
                  du_dxx   = (    q(gr_tmp+indexip1)
     *                          - 2.d0*u_pt
     *                          + q(gr_tmp+indexim1))*invhsq
                  du_dyy   = (    q(gr_tmp+indexjp1)
     *                          - 2.d0*u_pt
     *                          + q(gr_tmp+indexjm1))*invhsq
                  du_dzz   = (    q(gr_tmp+indexkp1)
     *                          - 2.d0*u_pt
     *                          + q(gr_tmp+indexkm1))*invhsq
                  if (rv1.lt.minradius) then
                     invalpha =  0.5d0 * bh_1_mass / minradius
                  else
                     invalpha =  0.5d0 * bh_1_mass / rv1
                  end if
                  if (bh_n .eq. 2) then
                     if (rv2.lt.minradius) then
                        invalpha =  invalpha + 0.5d0*bh_2_mass/minradius
                     else
                        invalpha =  invalpha + 0.5d0*bh_2_mass/rv2
                     end if
                  end if
                  alpha    = 1.d0 / invalpha
                  chi_pt   = 1.d0 + q(gr_chi+index)
                  kabkab   = (
     .                                     q(gr_A11+index)**2
     .                               +2.d0*q(gr_A12+index)**2
     .                               +2.d0*q(gr_A13+index)**2
     .                                    +q(gr_A22+index)**2
     .                               +2.d0*q(gr_A23+index)**2
     .                                    +q(gr_A33+index)**2
     .                                     )*chi_pt
                  beta     = (1.d0/8.d0)*alpha**7*kabkab
                  beta     = 0.d0
                  !
                  q(gr_uell_rhs+index) =
     .               du_dxx + du_dyy + du_dzz
     .             + beta * (1.d0+alpha*(1.d0+u_pt))**(-7)
               end if
            end do
         end do
      end do



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
      include      'param.inc'
      integer       nx, ny, nz, myid
      integer     proc_return_myid
      external    proc_return_myid
      real(kind=16)   relax
      external       relax
      real(kind=8)   resid
      integer     count
      save        count

      logical     ltrace
      parameter ( ltrace  = .false. )

      count = count + 1

      call load_pointers(gi)

      myid = proc_return_myid()
      nx   = gr_nx(gi)
      ny   = gr_ny(gi)
      nz   = gr_nz(gi)

      if (ltrace) then
       write(*,98)myid,'grid_relax: Relaxing grid: ',gi
       write(*,98)myid,'grid_relax: BH1:',bh_1_mass,bh_1_x,bh_1_y,bh_1_z
       if (bh_n.eq.2)
     . write(*,98)myid,'grid_relax: BH2:',bh_2_mass,bh_2_x,bh_2_y,bh_2_z
       write(*,99)myid,'grid_relax: n',bh_n, nx, ny, nz
       call field_dump_infoB(q(gr_uell),nx,ny, nz, 'uell')
       call field_dump_infoB(q(gr_uell_rhs),nx,ny, nz, 'uell_rhs')
      end if

c     resid = relax(q(gr_psi_elliptic),q(gr_psi_elliptic_rhs),q(gr_chr),
c    *                gr_h(gi),gr_nx(gi),gr_ny(gi),gr_nz(gi))
      resid = relax(q(gr_uell),q(gr_uell_rhs),q(gr_chr),
     .              q(gr_A11), q(gr_A12),     q(gr_A13),
     .              q(gr_A22), q(gr_A23),     q(gr_A33),
     .              q(gr_chi), q(gr_tmp),
     .              q(gr_x(gi)),q(gr_y(gi)),  q(gr_z(gi)), gr_h(gi),
     .              bh_1_mass, bh_1_x, bh_1_y, bh_1_z,
     .              bh_2_mass, bh_2_x, bh_2_y, bh_2_z,
     .              bh_n, nx, ny, nz)
      !
      ! For output and debugging purposes only:
      !
      if (.false.) then
      call get_resid(q(gr_uell_st1),q(gr_uell), q(gr_uell_rhs),
     .           q(gr_chr),
     .           q(gr_A11), q(gr_A12), q(gr_A13),
     .           q(gr_A22), q(gr_A23), q(gr_A33),
     .           q(gr_chi),
     .           q(gr_x(gi)),q(gr_y(gi)),q(gr_z(gi)),
     .           gr_h(gi), nx, ny, nz)
!     call get_resid(q(gr_phi_x),q(gr_phi),q(gr_rho),
!    *                gr_h(gi),gr_nx(gi),gr_ny(gi),gr_nz(gi))
         call field_out3d(q(gr_uell_st1), 1.d0*count, 'grid_relaxRES',
     *             q(gr_x(gi)), q(gr_x(gi)+nx-1),
     *             q(gr_y(gi)), q(gr_y(gi)+ny-1),
     *             q(gr_z(gi)), q(gr_z(gi)+nz-1),
     *             nx,   ny,    nz, myid)
      end if

      if (ltrace) then
         write(*,98)myid,'grid_relax: Finished with resid: ',resid
         call field_dump_infoB(q(gr_uell),nx,ny, nz, 'uell')
         call field_dump_infoB(q(gr_uell_rhs),nx,ny, nz, 'uell_rhs')
         call field_dump_infoB(q(gr_uell_st1),nx,ny, nz, 'uell_st1')
      end if

  95  format('[',I4,'] ',A,2I5,3F18.10)
  96  format('[',I4,'] ',A,I5,3F18.10)
  98  format('[',I4,'] ',A,5F15.6)
  99  format('[',I4,'] ',A,5I5)

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
      include      'param.inc'
      integer       nx, ny, nz

      !real(kind=8)  h
      real(kind=16)  get_resnrm
      external      get_resnrm
      real(kind=8)  grid_return_h
      external      grid_return_h
      logical       grid_is_local
      external      grid_is_local

      logical     ltrace
      parameter ( ltrace = .false. )

!     grid_return_resnrm = 9999.0d0

      call grid_get_dims(gi,nx,ny,nz)
      !h = grid_return_h(gi)

      call load_pointers(gi)

      if (.not. grid_is_local(gi)) then
         write(*,*) 'grid_return_resnrm: Only works on local grids'
      end if

      grid_return_resnrm =  get_resnrm(q(gr_uell),q(gr_uell_rhs),
     .              q(gr_A11), q(gr_A12),     q(gr_A13),
     .              q(gr_A22), q(gr_A23),     q(gr_A33),
     .              q(gr_chi),
     .              q(gr_x(gi)),q(gr_y(gi)),  q(gr_z(gi)), gr_h(gi),
     .              bh_1_mass, bh_1_x, bh_1_y, bh_1_z,
     .              bh_2_mass, bh_2_x, bh_2_y, bh_2_z,
     .                      bh_n,   nx,ny,nz)
c     grid_return_resnrm = 0.d0
c     grid_return_resnrm = get_resnrm(q(gr_psi_elliptic),
c    *                          q(gr_psi_elliptic_rhs),h,nx,ny,nz)


      return
      end      ! END: grid_return_resnrm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_ellpost:                                                             cc
cc               Things to do after an elliptic solve.                        cc
cc               For MHD project, this is where we clean the B field.         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_ellpost( gridnum )
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'param.inc'
      integer     nx, ny, nz
      real*8      hg
      ! Only set true for debugging
      logical     donothing
      parameter ( donothing = .false. )
      logical     useuminusone
      parameter ( useuminusone = .true. )
      logical     ltrace
      parameter ( ltrace = .true. )
      integer    i,j,k,index
      integer    indexip1, indexim1
      integer    indexjp1, indexjm1
      integer    indexkp1, indexkm1
      real*8      grad_x_psi, grad_y_psi, grad_z_psi
      real*8      l2change, chi_prev, sqrtchi, psi_prev, change
      integer     proc_return_myid, myid
      external    proc_return_myid
      real(kind=8)  alpha, invalpha, rv1,rv2
      real(kind=8)  x1,x2,y1,y2,z1,z2, x_pt,y_pt,z_pt
      real(kind=8)  minradius
      parameter (   minradius  = 0.001d0)

      call load_pointers(gridnum)

      myid = proc_return_myid()
      if(donothing) then 
         write(*,98)myid,'grid_ellpost: Doing nothing!'
         return
      end if

      nx   = gr_nx(gridnum)
      ny   = gr_ny(gridnum)
      nz   = gr_nz(gridnum)
      hg   = gr_h(gridnum)

      if (ltrace) then
         write(*,*) 'grid_ellpost:   On grid: ',gridnum
         call field_out3d(q(gr_ham), gr_t(gridnum), 'orig_ham',
     *             q(gr_x(gridnum)), q(gr_x(gridnum)+nx-1),
     *             q(gr_y(gridnum)), q(gr_y(gridnum)+ny-1),
     *             q(gr_z(gridnum)), q(gr_z(gridnum)+nz-1),
     *             nx,   ny,    nz, myid)
         call field_out3d(q(gr_chi), gr_t(gridnum), 'orig_chi',
     *             q(gr_x(gridnum)), q(gr_x(gridnum)+nx-1),
     *             q(gr_y(gridnum)), q(gr_y(gridnum)+ny-1),
     *             q(gr_z(gridnum)), q(gr_z(gridnum)+nz-1),
     *             nx,   ny,    nz, myid)
      end if

      !
      ! psi = 1/ alpha + u
      !
      l2change = 0.d0
      do k = 1, nz
         z_pt     = q(gr_z(gridnum)+k-1)
         z1       = z_pt - bh_1_z
         z2       = z_pt - bh_2_z
         do j = 1, ny
            y_pt     = q(gr_y(gridnum)+j-1)
            y1       = y_pt - bh_1_y
            y2       = y_pt - bh_2_y
            do i = 1, nx
               x_pt     = q(gr_x(gridnum)+i-1)
               x1       = x_pt - bh_1_x
               x2       = x_pt - bh_2_x
               !
               rv1      = sqrt(x1**2+y1**2+z1**2)
               rv2      = sqrt(x2**2+y2**2+z2**2)
               !
               index    = (k-1)*ny*nx + (j-1)*nx + (i-1)
               !
               ! chi = 1/(psi^4)
               !     = 1/((psi_prev+u)^4)
               ! psi_prev = 1.0/sqrt(sqrt(chi_prev))
               if (.false.) then
                  ! Can correct the existing chi if doing
                  ! an iterative procedure (not sure if this will work):
                  chi_prev        = q(gr_chi+index)
                  sqrtchi         = sqrt(chi_prev)
                  ! need to use just difference from 1 since 
                  ! u involves a perturbation about 1
                  psi_prev        = 1.d0 / sqrt(sqrtchi) - 1
               else
                  ! Or can use the analytic result
                  ! which will not be iterative:
                  if (rv1.lt.minradius) then
                     invalpha =  0.5d0 * bh_1_mass / minradius
                  else
                     invalpha =  0.5d0 * bh_1_mass / rv1
                  end if
                  if (bh_n .eq. 2) then
                     if (rv2.lt.minradius) then
                        invalpha =  invalpha + 0.5d0*bh_2_mass/minradius
                     else
                        invalpha =  invalpha + 0.5d0*bh_2_mass/rv2
                     end if
                  end if
                  alpha    =  1.d0 / invalpha
                  psi_prev =  invalpha
                  chi_prev =  1.d0 / (1.d0+psi_prev)**4
               end if
               ! Just for diagnostics:
               q(gr_tmp+index)     = psi_prev
               q(gr_flag+index)    = psi_prev+q(gr_uell+index)
               q(gr_uell_st2+index)= chi_prev
               !
               if (useuminusone) then
                  q(gr_chi+index) = 1.d0 /
     .                          (psi_prev+1.d0+q(gr_uell+index))**4
               else
                  q(gr_chi+index) = 1.d0/(psi_prev+q(gr_uell+index))**4
               end if
               change          = q(gr_chi+index) - chi_prev
               l2change        = l2change + change**2
               !
            end do
         end do
      end do

      l2change = sqrt(l2change / (nx) / (ny) / (nz) )

      if(ltrace)then
         write(*,*)'grid_ellpost:   l2change = ',l2change
         call field_out3d(q(gr_uell_st2), gr_t(gridnum), 'orig_chiprev',
     *             q(gr_x(gridnum)), q(gr_x(gridnum)+nx-1),
     *             q(gr_y(gridnum)), q(gr_y(gridnum)+ny-1),
     *             q(gr_z(gridnum)), q(gr_z(gridnum)+nz-1),
     *             nx,   ny,    nz, myid)
         call field_out3d(q(gr_uell), gr_t(gridnum), 'orig_uell',
     *             q(gr_x(gridnum)), q(gr_x(gridnum)+nx-1),
     *             q(gr_y(gridnum)), q(gr_y(gridnum)+ny-1),
     *             q(gr_z(gridnum)), q(gr_z(gridnum)+nz-1),
     *             nx,   ny,    nz, myid)
         call field_out3d(q(gr_tmp), gr_t(gridnum), 'orig_psi',
     *             q(gr_x(gridnum)), q(gr_x(gridnum)+nx-1),
     *             q(gr_y(gridnum)), q(gr_y(gridnum)+ny-1),
     *             q(gr_z(gridnum)), q(gr_z(gridnum)+nz-1),
     *             nx,   ny,    nz, myid)
         call field_out3d(q(gr_flag), gr_t(gridnum), 'new_psi',
     *             q(gr_x(gridnum)), q(gr_x(gridnum)+nx-1),
     *             q(gr_y(gridnum)), q(gr_y(gridnum)+ny-1),
     *             q(gr_z(gridnum)), q(gr_z(gridnum)+nz-1),
     *             nx,   ny,    nz, myid)
         call field_out3d(q(gr_chi),  gr_t(gridnum), 'new_chi',
     *             q(gr_x(gridnum)), q(gr_x(gridnum)+nx-1),
     *             q(gr_y(gridnum)), q(gr_y(gridnum)+ny-1),
     *             q(gr_z(gridnum)), q(gr_z(gridnum)+nz-1),
     *             nx,   ny,    nz, myid)
      end if

  95  format('[',I4,'] ',A,2I5,3F18.10)
  96  format('[',I4,'] ',A,I5,3F18.10)
  98  format('[',I4,'] ',A,5F18.10)
  99  format('[',I4,'] ',A,3I5)

      return
      end    ! END: grid_ellpost

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  get_resid:                                                                cc
cc             Compute:   res = rhs - L u                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_resid(res, u, rhs, chr,
     .           A11, A12, A13, A22, A23, A33, chi, x, y, z,
     .           myh, nx, ny, nz)
      implicit      none
      integer       nx, ny, nz
      real(kind=8)  chr(nx,ny,nz)
      real(kind=8)  res(nx,ny,nz), u(nx,ny,nz), rhs(nx,ny,nz), myh
      real(kind=8)  A11(nx,ny,nz), A12(nx,ny,nz), A13(nx,ny,nz),
     .              A22(nx,ny,nz), A23(nx,ny,nz), A33(nx,ny,nz),
     .              chi(nx,ny,nz), x(nx), y(ny), z(nz)
      include      'param.inc'
      include      'chr.inc'
      integer       i,j,k
      integer    myid, proc_return_myid
      external         proc_return_myid

      myid = proc_return_myid()
      !
      !  Compute: L u
      !
      call load_scal1D(res,0.d0,nx*ny*nz)
      call lop(res, u,
     .         A11, A12, A13,
     .         A22, A23, A33,
     .         chi,
     .         x, y, z, myh,
     .         bh_1_mass, bh_1_x, bh_1_y, bh_1_z,
     .         bh_2_mass, bh_2_x, bh_2_y, bh_2_z,
     .         bh_n, nx, ny, nz )

      !
      !  Compute: rhs - L u
      !
      call mat_mat_sub3d( res,rhs,res,nx,ny,nz)

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx, nx-1
!        if (NINT(chr(i,j,k)).eq.NINT(CHR_deco_bdy)) then
            res(i,j,k) = 0.d0
!        end if
      end do
      end do
      end do
      do k = 1, nz
      do i = 1, nx
      do j = 1, ny, ny-1
!        if (NINT(chr(i,j,k)).eq.NINT(CHR_deco_bdy)) then
            res(i,j,k) = 0.d0
!        end if
      end do
      end do
      end do
      do j = 1, ny
      do i = 1, nx
      do k = 1, nz, nz-1
!        if (NINT(chr(i,j,k)).eq.NINT(CHR_deco_bdy)) then
            res(i,j,k) = 0.d0
!        end if
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
      double precision function get_resnrm(u,rhs,
     .               A11,A12,A13,A22,A23,A33,chi,
     .               x,y,z,h,
     .               bh1m,bh1x,bh1y,bh1z,
     .               bh2m,bh2x,bh2y,bh2z,
     .                                 bh_n, nx,ny,nz)
      implicit      none
      integer       nx, ny, nz, bh_n
      real(kind=8)  u(nx,ny,nz),   rhs(nx,ny,nz)
      real(kind=8)  A11(nx,ny,nz), A12(nx,ny,nz), A13(nx,ny,nz)
      real(kind=8)  A22(nx,ny,nz), A23(nx,ny,nz), A33(nx,ny,nz)
      real(kind=8)  chi(nx,ny,nz)
      real(kind=8)  x(nx),         y(ny),         z(nz), h
      real(kind=8)  bh1m, bh1x, bh1y, bh1z
      real(kind=8)  bh2m, bh2x, bh2y, bh2z

      real(kind=8)  du_dxx, du_dyy, du_dzz, kabkab, tmp
      real(kind=8)  res, lhs_u, invh, invhsq, u_pt, chi_pt
      real(kind=8)  alpha, invalpha, rv1,rv2, beta
      real(kind=8)  x1,x2,y1,y2,z1,z2, x_pt,y_pt,z_pt
      integer       i, j, k
      logical     useuminusone
      parameter ( useuminusone = .true. )
      logical     usemaxnorm
      parameter ( usemaxnorm   = .false. )
      real(kind=8)  minradius
      parameter (   minradius  = 0.001d0)


      get_resnrm = 0.d0

      do k = 2, nz-1
            z_pt     = z(k)
            z1       = z_pt - bh1z
            z2       = z_pt - bh2z
            do j = 2, ny-1
               y_pt     = y(j)
               y1       = y_pt - bh1y
               y2       = y_pt - bh2y
               do i = 2, nx-1
                  !
                  !
                  invh     = 1.d0 / h
                  invhsq   = invh / h
                  x_pt     = x(i)
                  x1       = x_pt - bh1x
                  x2       = x_pt - bh2x
                  !
                  rv1      = sqrt(x1**2+y1**2+z1**2)
                  rv2      = sqrt(x2**2+y2**2+z2**2)
                  !
                  u_pt     = u(i,j,k)
                  !
                  du_dxx   = (    u(i+1,j,k)
     *                          - 2.d0*u_pt
     *                          + u(i-1,j,k)  )*invhsq
                  du_dyy   = (    u(i,j+1,k)
     *                          - 2.d0*u_pt
     *                          + u(i,j-1,k)  )*invhsq
                  du_dzz   = (    u(i,j,k+1)
     *                          - 2.d0*u_pt
     *                          + u(i,j,k-1)  )*invhsq
                  !
                  if (rv1.lt.minradius) then
                     invalpha =  0.5d0 * bh1m / minradius
                  else
                     invalpha =  0.5d0 * bh1m / rv1
                  end if
                  if (bh_n .eq. 2) then
                     if (rv2.lt.minradius) then
                        invalpha =  invalpha + 0.5d0 * bh2m / minradius
                     else
                        invalpha =  invalpha + 0.5d0 * bh2m / rv2
                     end if
                  end if
                  alpha    =  1.d0 / invalpha
                  !
                  chi_pt   = chi(i,j,k)
                  ! Try using this in terms of u:
                  chi_pt   = 1.d0 / (invalpha+1.d0+u_pt)**4
                  !
                  beta     =  (1.d0/8.d0)*alpha**7*(
     .                                     A11(i,j,k)**2
     .                               +2.d0*A12(i,j,k)**2
     .                               +2.d0*A13(i,j,k)**2
     .                                    +A22(i,j,k)**2
     .                               +2.d0*A23(i,j,k)**2
     .                                    +A33(i,j,k)**2
     .                                     )*chi_pt
                  kabkab   =  (1.d0/8.d0)*(
     .                                     A11(i,j,k)**2
     .                               +2.d0*A12(i,j,k)**2
     .                               +2.d0*A13(i,j,k)**2
     .                                    +A22(i,j,k)**2
     .                               +2.d0*A23(i,j,k)**2
     .                                    +A33(i,j,k)**2
     .                                     )*chi_pt
!    .                                     )*chi_pt**2
!    .                                     )
                  !kabkab      = 0.d0
                  tmp      =  alpha / (1.d0+alpha*u_pt+alpha)
                  !
                  if (useuminusone) then
!                    lhs_u    =   du_dxx + du_dyy + du_dzz
!    .                          + beta/(1.d0+alpha*(u_pt+1.d0))**7
                     lhs_u    =   du_dxx + du_dyy + du_dzz
     .                          + kabkab*tmp**7
                  else
                     lhs_u    =   du_dxx + du_dyy + du_dzz
     .                          + beta/(1.d0+alpha*u_pt)**7
                  end if
                  res      = lhs_u - rhs(i,j,k)
                  !
                  if (usemaxnorm) then
                     get_resnrm = max( get_resnrm, abs(res) )
                  else
                     get_resnrm = get_resnrm + res**2
                  end if
                  !
                  !
                  !
            end do
         end do
      end do

      if (usemaxnorm) then
         !
      else
         get_resnrm = sqrt(get_resnrm / nx /ny /nz)
      end if


      return
      end       ! END: get_resnrm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  lop:  Needs to parallel the finite differencing in relax()                cc
cc        NB: Needs to add the calculation to what resides in                 cc
cc            storage already!                                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lop(result,u,
     .               A11,A12,A13,A22,A23,A33,chi,
     .               x,y,z,h,
     .               bh1m,bh1x,bh1y,bh1z,
     .               bh2m,bh2x,bh2y,bh2z,
     .               bh_n, nx,ny,nz)
      implicit      none
      integer       nx, ny, nz, bh_n
      real(kind=8)  result(nx,ny,nz)
      real(kind=8)  u(nx,ny,nz)
      real(kind=8)  A11(nx,ny,nz), A12(nx,ny,nz), A13(nx,ny,nz)
      real(kind=8)  A22(nx,ny,nz), A23(nx,ny,nz), A33(nx,ny,nz)
      real(kind=8)  chi(nx,ny,nz)
      real(kind=8)  x(nx),         y(ny),         z(nz), h
      real(kind=8)  bh1m, bh1x, bh1y, bh1z
      real(kind=8)  bh2m, bh2x, bh2y, bh2z

      real(kind=8)  du_dxx, du_dyy, du_dzz, kabkab, tmp
      real(kind=8)  res, lhs_u, invh, invhsq, u_pt, chi_pt
      real(kind=8)  alpha, invalpha, rv1,rv2, beta
      real(kind=8)  x1,x2,y1,y2,z1,z2, x_pt,y_pt,z_pt
      integer       i, j, k
      integer    myid, proc_return_myid
      external         proc_return_myid
      logical     useuminusone
      parameter ( useuminusone = .true. )
      real(kind=8)  minradius
      parameter (   minradius  = 0.001d0)

      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         myid = proc_return_myid()
         write(*,99) myid,'lop: Enter:',nx,ny,nz
         write(*,99) myid,' bh_n: ',bh_n
         write(*,98) myid,' BH1: ',bh1m,bh1x,bh1y,bh1z
         write(*,98) myid,' BH2: ',bh2m,bh2x,bh2y,bh2z
         call field_out3d(A11, 0.d0, 'lopA11',
     *             x(1), x(nx),
     *             y(1), y(ny),
     *             z(1), z(nz),
     *             nx,   ny,    nz, 0d0)
         call field_out3d(A23, 0.d0, 'lopA23',
     *             x(1), x(nx),
     *             y(1), y(ny),
     *             z(1), z(nz),
     *             nx,   ny,    nz, 0d0)
      end if
  95  format('[',I4,'] ',A,2I5,3F18.10)
  96  format('[',I4,'] ',A,I5,3F18.10)
  98  format('[',I4,'] ',A,5F18.10)
  99  format('[',I4,'] ',A,3I5)

      do k = 2, nz-1
            z_pt     = z(k)
            z1       = z_pt - bh1z
            z2       = z_pt - bh2z
            do j = 2, ny-1
               y_pt     = y(j)
               y1       = y_pt - bh1y
               y2       = y_pt - bh2y
               do i = 2, nx-1
                  !
                  !
                  invh     = 1.d0 / h
                  invhsq   = invh / h
                  x_pt     = x(i)
                  x1       = x_pt - bh1x
                  x2       = x_pt - bh2x
                  !
                  rv1      = sqrt(x1**2+y1**2+z1**2)
                  rv2      = sqrt(x2**2+y2**2+z2**2)
                  !
                  u_pt     = u(i,j,k)
                  !
                  du_dxx   = (    u(i+1,j,k)
     *                          - 2.d0*u_pt
     *                          + u(i-1,j,k)  )*invhsq
                  du_dyy   = (    u(i,j+1,k)
     *                          - 2.d0*u_pt
     *                          + u(i,j-1,k)  )*invhsq
                  du_dzz   = (    u(i,j,k+1)
     *                          - 2.d0*u_pt
     *                          + u(i,j,k-1)  )*invhsq
                  !
                  if (rv1.lt.minradius) then
                     invalpha =  0.5d0 * bh1m / minradius
                  else
                     invalpha =  0.5d0 * bh1m / rv1
                  end if
                  if (bh_n .eq. 2) then
                     if (rv2.lt.minradius) then
                        invalpha =  invalpha + 0.5d0 * bh2m / minradius
                     else
                        invalpha =  invalpha + 0.5d0 * bh2m / rv2
                     end if
                  end if
                  alpha    =  1.d0 / invalpha
                  !
                  chi_pt   = chi(i,j,k)
                  ! Try using this in terms of u:
                  chi_pt   = 1.d0 / (invalpha+1.d0+u_pt)**4
                  !
                  beta     =  (1.d0/8.d0)*alpha**7*(
     .                                     A11(i,j,k)**2
     .                               +2.d0*A12(i,j,k)**2
     .                               +2.d0*A13(i,j,k)**2
     .                                    +A22(i,j,k)**2
     .                               +2.d0*A23(i,j,k)**2
     .                                    +A33(i,j,k)**2
     .                                     )*chi_pt
                  kabkab   =  (1.d0/8.d0)*(
     .                                     A11(i,j,k)**2
     .                               +2.d0*A12(i,j,k)**2
     .                               +2.d0*A13(i,j,k)**2
     .                                    +A22(i,j,k)**2
     .                               +2.d0*A23(i,j,k)**2
     .                                    +A33(i,j,k)**2
     .                                     )*chi_pt
!    .                                     )*chi_pt**2
!    .                                     )
                  !kabkab      = 0.d0
                  tmp      =  alpha / (1.d0+alpha*u_pt+alpha)
                  !
                  if (useuminusone) then
!                    lhs_u    =   du_dxx + du_dyy + du_dzz
!    .                          + beta/(1.d0+alpha*(u_pt+1.d0))**7
                     lhs_u    =   du_dxx + du_dyy + du_dzz
     .                          + kabkab*tmp**7
                  else
                     lhs_u    =   du_dxx + du_dyy + du_dzz
     .                          + beta/(1.d0+alpha*u_pt)**7
                  end if
                  result(i,j,k) = result(i,j,k) + lhs_u
                  !
                  !
            end do
         end do
      end do

      !
      ! Assuming dirichlet for now:
      !
      call load_scal_bound3d( result,0.d0, nx,ny,nz)

      if (ltrace) then
         write(*,*) 'lop: Done.'
         call field_dump_infoB(result,nx,ny, nz, 'result')
      end if

      return
      end      ! END: lop

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  relax:                                                                    cc
cc               Solving Hamiltonian constraint for u as per                  cc
cc               Brandt & Bergman                                             cc
cc               http://arxiv.org/abs/gr-qc/9703066                           cc
cc               nabla^2 u = -beta(1+alpha u)^{-7}                            cc
cc               1/alpha = Sigma_i  m_i / (2|r-r_i|)                          cc
cc               beta    = (1/8) alpha^6 K^{ab} K_{ab}                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function relax(u,rhs,chr,
     .                                A11,A12,A13,A22,A23,A33,chi,
     .                                work, x,y,z,h,
     .                                bh1m,bh1x,bh1y,bh1z,
     .                                bh2m,bh2x,bh2y,bh2z,
     .                                bh_n, nx,ny,nz)
      implicit      none
      integer       nx, ny, nz, bh_n
      real(kind=8)  u(nx,ny,nz),   rhs(nx,ny,nz), chr(nx,ny,nz)
      real(kind=8)  A11(nx,ny,nz), A12(nx,ny,nz), A13(nx,ny,nz)
      real(kind=8)  A22(nx,ny,nz), A23(nx,ny,nz), A33(nx,ny,nz)
      real(kind=8)  chi(nx,ny,nz), work(nx,ny,nz)
      real(kind=8)  x(nx),         y(ny),         z(nz), h
      real(kind=8)  bh1m, bh1x, bh1y, bh1z
      real(kind=8)  bh2m, bh2x, bh2y, bh2z
      include      'chr.inc'
      !
      real(kind=8)  residual, jacobian, dpsi
      real(kind=8)  du_dxx, du_dyy, du_dzz, kabkab, tmp
      real(kind=8)  res, lhs_u, invh, invhsq, u_pt, chi_pt
      real(kind=8)  alpha, invalpha, rv1,rv2, beta
      real(kind=8)  x1,x2,y1,y2,z1,z2, x_pt,y_pt,z_pt
      integer       pass, ib, kswitch
      integer       i, j, k
      integer     proc_return_myid, myid
      external    proc_return_myid
      logical     ltrace
      parameter ( ltrace = .false. )
      logical     useuminusone
      parameter ( useuminusone = .true. )
      logical     usemaxnorm
      parameter ( usemaxnorm   = .false. )
      real(kind=8)  minradius
      parameter (   minradius  = 0.001d0)

!     relax = 0.999999999d1

      !relax = sqrt( relax /nx /ny /nz )

      myid = proc_return_myid()
      if (ltrace) then
         !write(*,*)myid,'relax: ',relax
         write(*,98)myid,'relax: bh1m/x/y/z: ',bh1m,bh1x,bh1y,bh1z
         if (bh_n.gt.1)
     .   write(*,98)myid,'relax: bh2m/x/y/z: ',bh2m,bh2x,bh2y,bh2z
         !write(*,99)myid,'relax: bh_n: ',bh_n
         !call field_out1d(x,0.d0,'xrelax',x(1),x(nx),nx,myid)
         !call field_out1d(y,0.d0,'yrelax',y(1),y(ny),ny,myid)
         !call field_out1d(z,0.d0,'zrelax',z(1),z(nz),nz,myid)
      end if

      !
      ! Red-Black ordering of the interior only:
      !
      do pass = 1, 2
         kswitch = pass
         do k = 2, nz-1
            ib = kswitch + 1
            z_pt     = z(k)
            z1       = z_pt - bh1z
            z2       = z_pt - bh2z
            do j = 2, ny-1
               y_pt     = y(j)
               y1       = y_pt - bh1y
               y2       = y_pt - bh2y
               do i = ib, nx-1, 2
                  !
                  !
                  invh     = 1.d0 / h
                  invhsq   = invh / h
                  x_pt     = x(i)
                  x1       = x_pt - bh1x
                  x2       = x_pt - bh2x
                  !
                  rv1      = sqrt(x1**2+y1**2+z1**2)
                  rv2      = sqrt(x2**2+y2**2+z2**2)
                  !
                  u_pt     = u(i,j,k)
                  !
                  du_dxx   = (    u(i+1,j,k)
     *                          - 2.d0*u_pt
     *                          + u(i-1,j,k)  )*invhsq
                  du_dyy   = (    u(i,j+1,k)
     *                          - 2.d0*u_pt
     *                          + u(i,j-1,k)  )*invhsq
                  du_dzz   = (    u(i,j,k+1)
     *                          - 2.d0*u_pt
     *                          + u(i,j,k-1)  )*invhsq
                  !
                  if (rv1.lt.minradius) then
                     invalpha =  0.5d0 * bh1m / minradius
                  else
                     invalpha =  0.5d0 * bh1m / rv1
                  end if
                  if (bh_n .eq. 2) then
                     if (rv2.lt.minradius) then
                        invalpha =  invalpha + 0.5d0 * bh2m / minradius
                     else
                        invalpha =  invalpha + 0.5d0 * bh2m / rv2
                     end if
                  end if
                  alpha    =  1.d0 / invalpha
                  !
                  chi_pt   = chi(i,j,k)
                  ! Try using this in terms of u:
                  chi_pt   = 1.d0 / (invalpha+1.d0+u_pt)**4
                  !
     .            ! Need to calculate: tildeA_ij tildeA^ij
     .            !
     .            ! Using:
     .            !            g^{ik} = psi^{-4} delta^{ij}
     .            !   and:        chi = psi^(-4)
     .            !   and:  tildeA_ij = psi^2*A_ij 
     .            !   and:       A_ij = stored grid functions in the code
     .            !
     .            ! Calculate: tildeA_ij tildeA^ij =
     .            !  = g^{ik}g^{jl} tilde A_{ij} tilde A_{kl} 
     .            !  = psi^{-8} tilde[A11**2+2*A12**2+2*A13**2+A22**2++2*A23**2+A33**2]
     .            !  = psi^{-4}      [A11**2+2*A12**2+2*A13**2+A22**2++2*A23**2+A33**2]
     .            !  = chi           [A11**2+2*A12**2+2*A13**2+A22**2++2*A23**2+A33**2]
     .            ! These are the non-tilde quantities
     .            ! So ultimately we have:
     .            ! tilde A_{ij} tilde A^{ij} =A11**2+A22**2+A33**2
                  kabkab   =  (1.d0/8.d0)*(
     .                                     A11(i,j,k)**2
     .                               +2.d0*A12(i,j,k)**2
     .                               +2.d0*A13(i,j,k)**2
     .                                    +A22(i,j,k)**2
     .                               +2.d0*A23(i,j,k)**2
     .                                    +A33(i,j,k)**2
     .                                     )*chi_pt
!    .                                     )*chi_pt**2
!    .                                     )
                  tmp      =  alpha / (1.d0+alpha*u_pt+alpha)
                  !work(i,j,k) = kabkab
                  !work(i,j,k) = beta
                  !work(i,j,k) = alpha
                  !work(i,j,k) = tmp
                  work(i,j,k) = kabkab*tmp**7
                  !kabkab      = 0.d0
                  !
                  if (useuminusone) then
                     lhs_u    =   du_dxx + du_dyy + du_dzz
     .                          + kabkab*tmp**7
!    .                          + beta/(1.d0+alpha*(u_pt+1.d0))**7
                     jacobian = -6.d0*invhsq 
     .                          -7.d0*kabkab*tmp**8
!    .                   -7.d0*beta*alpha/(1.d0+alpha*(u_pt+1.d0))**8
                  else
                     lhs_u    =   du_dxx + du_dyy + du_dzz
     .                          + beta/(1.d0+alpha*u_pt)**7
                     jacobian = -6.d0*invhsq 
     .                          -7.d0*beta*alpha/(1.d0+alpha*u_pt)**8
                  end if
                  res      = lhs_u - rhs(i,j,k)
                  u(i,j,k) = u_pt  - res / jacobian
                  !
                  if (usemaxnorm) then
                     relax = max( relax, abs(res) )
                  else
                     relax = relax + res**2
                  end if
                  !
                  !
                  !
               end do
               ib = 5 - ib
            end do
            kswitch = 3 - kswitch
         end do
      end do

      if (usemaxnorm) then
         !relax = max( relax, abs(res) )
      else
         relax = sqrt( relax/nx/ny/nz )
      end if

      if (ltrace) then
         write(*,94)myid,'relax: Norm: ',relax
         call field_out3d(work, 0.d0, 'work_kabkab7',
         !call field_out3d(work, 0.d0, 'work_AA',
         !call field_out3d(work, 0.d0, 'work_alpha',
         !call field_out3d(work, 0.d0, 'work_beta',
     *             x(1), x(nx), y(1), y(ny), z(1), z(nz),
     *             nx,   ny,    nz, myid)
         call field_out3d(rhs, 0.d0, 'relax_rhs',
     *             x(1), x(nx), y(1), y(ny), z(1), z(nz),
     *             nx,   ny,    nz, myid)
      end if

  94  format('[',I4,'] ',A,5G13.6)
  95  format('[',I4,'] ',A,2I5,3F18.10)
  96  format('[',I4,'] ',A,I5,3F18.10)
  98  format('[',I4,'] ',A,5F13.6)
  99  format('[',I4,'] ',A,3I5)

      return
      end      ! END: relax

