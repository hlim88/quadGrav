#include "cctk.h"


!----------------------------------------------------------------
!
! 
!_________________________________________________________________
#ifdef GR_MHD
module bssn_initial
#else
module initial
#endif

use params
implicit none

CONTAINS
#ifdef GR_MHD
  subroutine bssn_initialdata(u, u2, v, w, par)
#else
  subroutine initialdata(u, u2, v, w, imask, par)
#endif
    use GF

#if defined GR_MHD
  use bssn_auxvars
#else
  use auxvars
#endif
    implicit none


    type(gridfunction), dimension(NU_G)   :: u, u2
    type(gridfunction), dimension(NV_G)   :: v
    type(gridfunction), dimension(NW)     :: w
#ifdef GR_MHD
    ! imask is now passed by HyperInit
#else
    CCTK_INT, dimension(:,:,:)            :: imask
#endif
    CCTK_REAL, dimension(NPAR)               :: par
    integer myrank,ierr

!    include 'bssnfunctions_fortran.inc'
    include 'mpif.h'

    !! local vars
    CCTK_INT, dimension(3)  :: shp
    logical, parameter      :: ltrace2 = .false.
    logical                 :: setchi
      logical     useuminusone
      parameter ( useuminusone = .true. )


    CCTK_REAL myl2norm3d

! travis adding new stuff here:
    integer :: i_idatatype,i_btype
    CCTK_REAL, dimension(:,:,:),pointer  :: xx, yy, zz
    integer                 :: i,j,k,m, idtype, nx, ny, nz

! EM Variables:
    CCTK_REAL :: v_bamp, v_e1amp, v_e2amp, v_sfamp, v_sfamp_axn

! black hole variables:
    CCTK_INT  :: n_bhs
    CCTK_REAL :: mass1,mass2,mass3,mass4
    CCTK_REAL :: bh1x,bh1y,bh1z
    CCTK_REAL :: bh2x,bh2y,bh2z
    CCTK_REAL :: bh3x,bh3y,bh3z
    CCTK_REAL :: bh4x,bh4y,bh4z
    CCTK_REAL :: rv1,rv2,rv3,rv4,vpsibl
    CCTK_REAL :: x1,y1,z1
    CCTK_REAL :: x2,y2,z2
    CCTK_REAL :: x3,y3,z3
    CCTK_REAL :: x4,y4,z4
    CCTK_REAL :: vn1(1:3),vn2(1:3),vn3(1:3),vn4(1:3)
    CCTK_REAL :: vp1(1:3),vp2(1:3),vp3(1:3),vp4(1:3)
    CCTK_REAL :: vs1(1:3),vs2(1:3),vs3(1:3),vs4(1:3)
    CCTK_REAL :: deltaij(1:3,1:3)
    CCTK_REAL :: epijk(1:3,1:3,1:3)
    CCTK_REAL :: spin1,spin2,spin3,spin4
    CCTK_REAL :: spin1_th,spin1_phi,spin2_th,spin2_phi
    CCTK_REAL :: spin3_th,spin3_phi,spin4_th,spin4_phi

! some temporary variables for building things...
    CCTK_INT   :: i1,i2,i3,i4,i5,i6
    CCTK_REAL  :: v1,v2,v3,v4
    CCTK_REAL  :: vt1,vt2,vt3,vt4

    CCTK_REAL  :: Gamt(3), Atd(3,3), gtd(3,3), chi, alpha, Betau(3), trK

! for the noise stability test
    CCTK_REAL, DIMENSION(50) :: harvest
    CCTK_REAL  :: rb_amp

! variables to construct a better u function:
    CCTK_INT :: i_better_u
    CCTK_REAL :: vp1tot,vp2tot,vp3tot,vp4tot
    CCTK_REAL :: vpsibl_u,vpsibl_u2,v_u_corr
    CCTK_REAL :: v_u_j1,v_u_j2,v_u_j3,v_u_j4
    CCTK_REAL :: v_u_p1,v_u_p2,v_u_p3,v_u_p4
    CCTK_REAL :: v_u_c1,v_u_c2,v_u_c3,v_u_c4
    CCTK_REAL :: l_r,amp_capr,amp_capj,amp_capp
    CCTK_REAL :: mu_j,p2_mu_j,u0_j,u2_j
    CCTK_REAL :: mu_p,p2_mu_p,u0_p,u2_p

! variables for the boson star
    CCTK_INT  :: lev, interp, istat
    CCTK_INT  :: dnx, dny, dnz, dnr
    CCTK_REAL :: minx0,miny0,minz0, hi, vx, vy, vz
    CCTK_REAL :: minx, miny, minz, omega1, omega2, xcenter, ycenter, zcenter
    CCTK_REAL, allocatable, dimension(:)   :: fdata1d      
    CCTK_REAL, allocatable, dimension(:)   :: x1d, y1d, z1d, coords
    CCTK_REAL, allocatable, dimension(:)   :: coords1d

! variables specific to Einstein-Maxwell-Dilaton  
    CCTK_INT  :: emd_bh_type 
    CCTK_REAL :: phi_0, dil_alpha, charge_magnetic, charge_electric, charge1
    CCTK_REAL :: rbar, mass1_sq, charge1_sq, dil_alpha_sq, discriminant
    CCTK_REAL :: sqrt_disc, r_plus, r_minus, rbar_1, rbar_2, rbar_H, beta_0 
    CCTK_REAL :: r_mi_rH, r_pl_rH, r_pl_r1, r_pl_r2, alpha_sq 
    CCTK_REAL :: tmp_dil, tmp_Erbar, tmp_Brbar 
    CCTK_REAL :: emd_gaussian_amp, emd_gaussian_cent, emd_gaussian_width
    CCTK_REAL :: emd_gaussian_pert
! variables specific to the rotating BH solution in EMD 
    CCTK_REAL :: xbar, ybar, zbar, rhobar
    CCTK_REAL :: cos_theta, sin_theta, cos_phi, sin_phi
    CCTK_REAL :: r_BL, rho_sq, Delta, BB, a_spin, mass, v_boost, inv_chi
    CCTK_REAL :: CC, CC_tothird, dCC_dr, dCC_dtheta
    CCTK_REAL :: Sigma, dSigma_dr, dSigma_dtheta, K_rphi, K_thetaphi
    CCTK_REAL :: gamtilde_tmp, K_ij_tmp1, K_ij_tmp2, betau_phi
    CCTK_REAL :: F_tr, F_ttheta, F_rphi, F_thetaphi
    CCTK_REAL :: E_up_tmp1, E_up_tmp2, B_up_tmp1, B_up_tmp2
! variables specific to the Einstein-Maxwell-Dilaton-Axion, Kerr-Sen (HL)
    CCTK_REAL :: axn_alpha, axn_mass, axn_infty, sen_alpha
    CCTK_REAL :: rbar_plus, rbar_minus, a_ang_par
! variables specific to the Rasheed-Larsen
    CCTK_REAL :: H_comm1, H_comm2, H_1_RL, H_2_RL, H_3_RL, B_RL 
    CCTK_REAL :: p_mag, q_elec, Sigma_pq, c_0_RL, c_1_RL
    CCTK_REAL :: dH_1_RL_dr, dH_1_RL_dtheta, dH_2_RL_dr, dH_2_RL_dtheta
    CCTK_REAL :: dSigma_pq_dr, dSigma_pq_dtheta, Qt_elec, Pt_mag

    interface
      subroutine read_sdf_file(filename, f, xf, yf, zf, nx, ny, nz, &
                            level, fdata, coords, dnx, dny, dnz,    &
                            interp, cctk_llb, cctk_gshp)
          implicit none
          CCTK_INT        nx, ny, nz, level
          CCTK_INT        dnx, dny, dnz
          CCTK_INT        interp, cctk_llb(3), cctk_gshp(3)
          CCTK_REAL       f(nx,ny,nz)
          CCTK_REAL       xf(nx), yf(ny), zf(nz)
          CCTK_REAL       fdata(dnx,dny,dnz)
          CCTK_REAL       coords(dnx+dny+dnz)
          character*(*)   filename
      end subroutine read_sdf_file
	
      subroutine read_sdf_file1d(filename, f, xf, yf, zf, nx, ny, nz,    &
                            level, fdata1d, coords1d, dnr, interp, xc, yc, zc)
          implicit none
          CCTK_INT        nx, ny, nz, level
          CCTK_INT        dnr
          CCTK_INT        interp
          CCTK_REAL       f(nx,ny,nz)
          CCTK_REAL       xf(nx), yf(ny), zf(nz)
          CCTK_REAL       fdata1d(dnr)
          CCTK_REAL       coords1d(dnr)
          character*(*)   filename
          CCTK_REAL        xc, yc, zc
      end subroutine read_sdf_file1d		
		
      subroutine sdf_file_mem(filename, level, shp)
          implicit none
          character*(*)   filename
          CCTK_INT         ::  shp(3), level
      end subroutine sdf_file_mem
    end interface


! get the coordinates
    xx => w(H_XPHYS)%d
    yy => w(H_YPHYS)%d
    zz => w(H_ZPHYS)%d

    nx = nint(par(P_NX))
    ny = nint(par(P_NY))
    nz = nint(par(P_NZ))

! read in the black hole parameters:
    n_bhs = par(P_BH_N)
! not using this currently - could later...
    i_idatatype = par(P_BH_TYPE)
! EM fields initial data type:
    i_btype = par(P_INITIAL_B)
    v_bamp = par(P_BAMP)
    v_e1amp = par(P_E1_AMP)
    v_e2amp = par(P_E2_AMP)
    v_sfamp = par(P_SF_AMP)
    v_sfamp_axn = par(P_SF_AMP_AXN)

! black hole 1:
    mass1 = par(P_BH_1_MASS)
    bh1x  = par(P_BH_1_X)
    bh1y  = par(P_BH_1_Y)
    bh1z  = par(P_BH_1_Z)
! we give the black holes momentum here, not velocity - 
! thus needs to be scaled if the mass changes...
    vp1(1) = par(P_BH_1_PX)
    vp1(2) = par(P_BH_1_PY)
    vp1(3) = par(P_BH_1_PZ)
    vp1tot = sqrt(vp1(1)**2+vp1(2)**2+vp1(3)**2)
    spin1     = par(P_BH_1_SPIN)
    spin1_th  = par(P_BH_1_STH)
    spin1_phi = par(P_BH_1_SPHI)
    vs1(1) = spin1*sin(spin1_th)*cos(spin1_phi)
    vs1(2) = spin1*sin(spin1_th)*sin(spin1_phi)
    vs1(3) = spin1*cos(spin1_th)

! black hole 2 (need to set BH_N = 2 to use...)
    mass2 = par(P_BH_2_MASS)
    bh2x  = par(P_BH_2_X)
    bh2y  = par(P_BH_2_Y)
    bh2z  = par(P_BH_2_Z)
    vp2(1) = par(P_BH_2_PX)
    vp2(2) = par(P_BH_2_PY)
    vp2(3) = par(P_BH_2_PZ)
    vp2tot = sqrt(vp2(1)**2+vp2(2)**2+vp2(3)**2)
    spin2     = par(P_BH_2_SPIN)
    spin2_th  = par(P_BH_2_STH)
    spin2_phi = par(P_BH_2_SPHI)
    vs2(1) = spin2*sin(spin2_th)*cos(spin2_phi)
    vs2(2) = spin2*sin(spin2_th)*sin(spin2_phi)
    vs2(3) = spin2*cos(spin2_th)

! black hole 3:
    mass3 = par(P_BH_3_MASS)
    bh3x  = par(P_BH_3_X)
    bh3y  = par(P_BH_3_Y)
    bh3z  = par(P_BH_3_Z)
    vp3(1) = par(P_BH_3_PX)
    vp3(2) = par(P_BH_3_PY)
    vp3(3) = par(P_BH_3_PZ)
    vp3tot = sqrt(vp3(1)**2+vp3(2)**2+vp3(3)**2)
    spin3     = par(P_BH_3_SPIN)
    spin3_th  = par(P_BH_3_STH)
    spin3_phi = par(P_BH_3_SPHI)
    vs3(1) = spin3*sin(spin3_th)*cos(spin3_phi)
    vs3(2) = spin3*sin(spin3_th)*sin(spin3_phi)
    vs3(3) = spin3*cos(spin3_th)

! black hole 4 (need to set BH_N = 4 to use...)
    mass4 = par(P_BH_4_MASS)
    bh4x  = par(P_BH_4_X)
    bh4y  = par(P_BH_4_Y)
    bh4z  = par(P_BH_4_Z)
    vp4(1) = par(P_BH_4_PX)
    vp4(2) = par(P_BH_4_PY)
    vp4(3) = par(P_BH_4_PZ)
    vp4tot = sqrt(vp4(1)**2+vp4(2)**2+vp4(3)**2)
    spin4     = par(P_BH_4_SPIN)
    spin4_th  = par(P_BH_4_STH)
    spin4_phi = par(P_BH_4_SPHI)
    vs4(1) = spin4*sin(spin4_th)*cos(spin4_phi)
    vs4(2) = spin4*sin(spin4_th)*sin(spin4_phi)
    vs4(3) = spin4*cos(spin4_th)

! define parameters specific to EMD BHs
    emd_bh_type = par(P_EMD_BH_TYPE)
    !emd_gaussian_amp   = par(P_EMD_GAUSSIAN_AMP) 
    !emd_gaussian_cent  = par(P_EMD_GAUSSIAN_CENT) 
    !emd_gaussian_width = par(P_EMD_GAUSSIAN_WIDTH) 

! define the antisymmetric tensor, and the identity 3-metric:
    epijk(:,:,:) = 0.d0
    epijk(1,2,3) = 1.d0;epijk(2,3,1) = 1.d0;epijk(3,1,2) = 1.d0
    epijk(1,3,2) = -1.d0;epijk(3,2,1) = -1.d0;epijk(2,1,3) = -1.d0

    deltaij(:,:) = 0.d0
    deltaij(1,1) = 1.d0;deltaij(2,2) = 1.d0;deltaij(3,3) = 1.d0

    write(0,*) '|Chi|: ',myl2norm3d(u(H_CHI)%d,nx,ny,nz),xx(2,1,1)-xx(1,1,1)
    if (myl2norm3d(u(H_CHI)%d,nx,ny,nz).gt.0 .and. .false.) then
      write(0,*)'...Not setting chi, assuming elliptic solver did so'
      setchi = .false.
    else
      write(0,*)'...Chi does not appear set, so setting in initial.f90'
      setchi = .true.
      if (useuminusone) then
         call load_scal1D(v(H_UELL)%d,0.0, nx*ny*nz)
      else
         call load_scal1D(v(H_UELL)%d,1.0, nx*ny*nz)
      end if
    end if

    if (n_bhs.le.0) then
       print *,'No Holes!'
    else if (n_bhs.gt.4) then
       print *,'Too Many Holes!'
    endif

    if (ltrace2) then
      write(0,*)'...Begin BSSN Initial'
    end if

    shp = shape(u(1)%d)
    if (ltrace2) write(*,*)'begin shape = ',shp(1), shp(2), shp(3)

    if (shp(1) .ne. nint(par(P_NX))) then
      write(0,*)'*** initial.f90:  Problem with inconsistent', & 
     &              ' array sizes'
      write(0,*) 'shp(1)    = ',shp(1)
      write(0,*) 'par(P_NX) = ',nint(par(P_NX))
  end if
    if (shp(2) .ne. nint(par(P_NY))) then
      write(0,*)'*** initial.f90:  Problem with inconsistent', &
   &              ' array sizes'
      write(0,*) 'shp(2)    = ',shp(2)
      write(0,*) 'par(P_NY) = ',nint(par(P_NY))
    end if
    if (shp(3) .ne. nint(par(P_NZ))) then
      write(0,*)'*** initial.f90:  Problem with inconsistent', &
   &              ' array sizes'
      write(0,*) 'shp(3)    = ',shp(3)
      write(0,*) 'par(P_NZ) = ',nint(par(P_NZ))
    end if


    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    ! -----------------------   INITIAL DATA  -----------------------
    if ( myrank .eq. 0 ) then
      write(0,*)'~~~~BSSN INITIAL DATA idtype = ',nint(par(P_IDTYPE))
    end if

    ! --- Check to see if we need to initialize the geometry here,
    !     or if it comes in already set.  A fixed geometry is defined
    !     here.
!    write(0,*) 'GEOMETRY TYPE = ',nint(par(P_GEOMETRY_TYPE))
    if (nint(par(P_GEOMETRY_TYPE)) .eq. 0) then
      ! assume that geometry is set elsewhere, as in hyperGR.  Do nothing.
!       print *,'travis here, lets start setting things up - need coordinates'
!     carlos: formulas from eq (66) for psi = 1 + m1/(2*r1) + m2/(2*r2)
!     and (70) for the extrinsic curvature Atilde_ij = psi^{-6} * Abar_ij
!     in  http://relativity.livingreviews.org/Articles/lrr-2000-5/

       do k=1,nz
          do j=1,ny
             do i=1,nx
                
                x1 = xx(i,j,k) - bh1x
                y1 = yy(i,j,k) - bh1y
                z1 = zz(i,j,k) - bh1z
                
                rv1 = sqrt(x1**2+y1**2+z1**2)
                if (rv1.lt.0.0001d0) then
                   rv1=0.0001d0
                endif
                vn1(1) = x1/rv1
                vn1(2) = y1/rv1
                vn1(3) = z1/rv1
                
                if (n_bhs.eq.2) then
                   x2 = xx(i,j,k) - bh2x
                   y2 = yy(i,j,k) - bh2y
                   z2 = zz(i,j,k) - bh2z
                   
                   rv2 = sqrt(x2**2+y2**2+z2**2)
                   if (rv2.lt.0.0001d0) then
                      rv2=0.0001d0
                   endif
                   vn2(1) = x2/rv2
                   vn2(2) = y2/rv2
                   vn2(3) = z2/rv2
                endif


                ! Brill-Lindquist conformal factor:
                vpsibl = 1.d0 + mass1/(2.d0*rv1)
                if (n_bhs.eq.2) then
                   vpsibl = vpsibl + mass2/(2.d0*rv2)
                endif

! hard-code in for the time being:
                i_better_u = 1
                if (i_better_u.eq.1) then
                   v_u_corr = 0.d0
                   
                   if (abs(spin1).gt.0.000001) then
                      amp_capj = 4.d0*spin1/mass1**2
                      amp_capr = 2.d0*rv1/mass1
                      l_r = 1.d0/(1.d0 + amp_capr)
                      u0_j = (l_r + l_r**2 + l_r**3 - 4.d0*l_r**4 + 2.d0*l_r**5)/40.d0
                      u2_j = -l_r**5/20.d0
                      mu_j =         vn1(1)*vs1(1)
                      mu_j =  mu_j + vn1(2)*vs1(2)
                      mu_j = (mu_j + vn1(3)*vs1(3))/abs(spin1)
                      p2_mu_j = (3.d0*mu_j**2 - 1.d0)/2.d0
                      v_u_j1 = amp_capj**2*(u0_j+u2_j*amp_capr**2*p2_mu_j)
                      v_u_corr = v_u_corr + v_u_j1
                   endif
                   
                   if (vp1tot.gt.0.000001) then
                      amp_capp = 2.d0*vp1tot/mass1
                      amp_capr = 2.d0*rv1/mass1
                      l_r = 1.d0/(1.d0 + amp_capr)
                      u0_p = l_r - 2.d0*l_r**2 + 2.d0*l_r**3 
                      u0_p = (u0_p - l_r**4 + 0.2d0*l_r**5)*(5.d0/32.d0)
                      u2_p = 15.d0*l_r + 132.d0*l_r**2 + 53.d0*l_r**3
                      u2_p = u2_p + 96.d0*l_r**4 + 82.d0*l_r**5
                      u2_p = u2_p + (84.d0/amp_capr)*(l_r**5+log(l_r)/amp_capr)
                      u2_p = (u2_p)/(80.d0*amp_capr)
                      mu_p =        vn1(1)*vp1(1)/vp1tot
                      mu_p = mu_p + vn1(2)*vp1(2)/vp1tot
                      mu_p = mu_p + vn1(3)*vp1(3)/vp1tot
                      p2_mu_p = (3.d0*mu_p**2 - 1.d0)/2.d0
                      v_u_p1 = amp_capp**2*(u0_p+u2_p*p2_mu_p)
                      v_u_corr = v_u_corr + v_u_p1
                   endif
                   
                   if ((vp1tot.gt.0.000001).and.(abs(spin1).gt.0.000001)) then
                      v1 =      (vp1(2)*vs1(3)-vp1(3)*vs1(2))*vn1(1)
                      v1 = v1 + (vp1(3)*vs1(1)-vp1(1)*vs1(3))*vn1(2)
                      v1 = v1 + (vp1(1)*vs1(2)-vp1(2)*vs1(1))*vn1(3)
                      v1 = v1*(16.d0/mass1**4)*rv1
                      
                      amp_capr = 2.d0*rv1/mass1
                      l_r = 1.d0/(1.d0 + amp_capr)
                      
                      v2 = 1.d0 + 5.d0*amp_capr + 10.d0*amp_capr**2
                      
                      v_u_c1 = (v1*v2*l_r**5)/80.d0
                      v_u_corr = v_u_corr + v_u_c1
                   endif
                   
                   if (n_bhs.eq.2) then
                      
                      if (abs(spin2).gt.0.000001) then
                         amp_capj = 4.d0*spin2/mass2**2
                         amp_capr = 2.d0*rv2/mass2
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_j = (l_r + l_r**2 + l_r**3 - 4.d0*l_r**4 + 2.d0*l_r**5)/40.d0
                         u2_j = -l_r**5/20.d0
                         mu_j =         vn2(1)*vs2(1)
                         mu_j =  mu_j + vn2(2)*vs2(2)
                         mu_j = (mu_j + vn2(3)*vs2(3))/abs(spin2)
                         p2_mu_j = (3.d0*mu_j**2 - 1.d0)/2.d0
                         v_u_j2 = amp_capj**2*(u0_j+u2_j*amp_capr**2*p2_mu_j)
                         v_u_corr = v_u_corr + v_u_j2
                      endif
                      
                      if (vp2tot.gt.0.000001) then
                         amp_capp = 2.d0*vp2tot/mass2
                         amp_capr = 2.d0*rv2/mass2
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_p = l_r - 2.d0*l_r**2 + 2.d0*l_r**3 
                         u0_p = (u0_p - l_r**4 + 0.2d0*l_r**5)*(5.d0/32.d0)
                         u2_p = 15.d0*l_r + 132.d0*l_r**2 + 53.d0*l_r**3
                         u2_p = u2_p + 96.d0*l_r**4 + 82.d0*l_r**5
                         u2_p = u2_p + (84.d0/amp_capr)*(l_r**5+log(l_r)/amp_capr)
                         u2_p = (u2_p)/(80.d0*amp_capr)
                         mu_p =        vn2(1)*vp2(1)/vp2tot
                         mu_p = mu_p + vn2(2)*vp2(2)/vp2tot
                         mu_p = mu_p + vn2(3)*vp2(3)/vp2tot
                         p2_mu_p = (3.d0*mu_p**2 - 1.d0)/2.d0
                         v_u_p2 = amp_capp**2*(u0_p+u2_p*p2_mu_p)
                         v_u_corr = v_u_corr + v_u_p2
                      endif
                      
                      if ((vp2tot.gt.0.000001).and.(abs(spin2).gt.0.000001)) then
                         v1 =      (vp2(2)*vs2(3)-vp2(3)*vs2(2))*vn2(1)
                         v1 = v1 + (vp2(3)*vs2(1)-vp2(1)*vs2(3))*vn2(2)
                         v1 = v1 + (vp2(1)*vs2(2)-vp2(2)*vs2(1))*vn2(3)
                         v1 = v1*(16.d0/mass2**4)*rv2
                         
                         amp_capr = 2.d0*rv2/mass2
                         l_r = 1.d0/(1.d0 + amp_capr)
                         
                         v2 = 1.d0 + 5.d0*amp_capr + 10.d0*amp_capr**2
                         
                         v_u_c2 = (v1*v2*l_r**5)/80.d0
                         v_u_corr = v_u_corr + v_u_c2
                      endif

                   endif
                  
! vpsibl_u will be used for the conformal factor, 
                   vpsibl_u  = vpsibl + v_u_corr
! vpsibl_u2 is for the Aij terms...
! since the corrections are first order...
! adding half of the correction seems to give the best results...
! update - do a fit for spin = 0.6...
!                   vpsibl_u2 = vpsibl
!                   vpsibl_u2 = vpsibl + v_u_corr/2.d0
!                   vpsibl_u2 = vpsibl + 0.313d0*v_u_corr
                   vpsibl_u2 = vpsibl + v_u_corr
                else
                   vpsibl_u  = vpsibl
                   vpsibl_u2 = vpsibl
                endif
                
                u(H_ALPHA)%d(i,j,k) = 1.d0/vpsibl_u**2
                v2 =  1.d0/vpsibl_u**4

                if (setchi) then
                   u(H_CHI)%d(i,j,k) = v2
                end if
                u(H_TRK)%d(i,j,k) =  0.0d0
                
                u(H_SHIFT1)%d(i,j,k) =  0.0d0
                u(H_SHIFT2)%d(i,j,k) =  0.0d0
                u(H_SHIFT3)%d(i,j,k) =  0.0d0
                   
                u(H_GAM1)%d(i,j,k) =  0.0d0
                u(H_GAM2)%d(i,j,k) =  0.0d0
                u(H_GAM3)%d(i,j,k) =  0.0d0
                
                u(H_GB1)%d(i,j,k) =  0.0d0
                u(H_GB2)%d(i,j,k) =  0.0d0
                u(H_GB3)%d(i,j,k) =  0.0d0
                
                u(H_GT11)%d(i,j,k) =  1.0d0
                u(H_GT12)%d(i,j,k) =  0.0d0
                u(H_GT13)%d(i,j,k) =  0.0d0
                u(H_GT22)%d(i,j,k) =  1.0d0
                u(H_GT23)%d(i,j,k) =  0.0d0
                u(H_GT33)%d(i,j,k) =  1.0d0
                

                do i1 = 1,3
                   do i2 = i1,3
                      
                      v2 = 0.d0
                      do i3 = 1,3
                         do i4 = 1,3
                            vt1 = epijk(i1,i3,i4)*vs1(i3)*vn1(i4)*vn1(i2)
                            vt2 = epijk(i2,i3,i4)*vs1(i3)*vn1(i4)*vn1(i1)
                            v2 = v2 + vt1 + vt2
                         enddo
                      enddo
                      
                      v3 = vp1(i1)*vn1(i2) + vp1(i2)*vn1(i1)
                      vt1 = 0.d0
                      do i3 = 1,3
                         vt1 = vt1 + vp1(i3)*vn1(i3)
                      enddo
                      vt1 = vt1*(vn1(i1)*vn1(i2) - deltaij(i1,i2))
                      v3 = v3 + vt1
                      
                      v1 = 3.d0/(vpsibl_u2**6*rv1**3)
                      v4 = v1*(v2+(rv1/2.d0)*v3)
                      
                      if (n_bhs.eq.2) then
                         v2 = 0.d0
                         do i3 = 1,3
                            do i4 = 1,3
                               vt1 = epijk(i1,i3,i4)*vs2(i3)*vn2(i4)*vn2(i2)
                               vt2 = epijk(i2,i3,i4)*vs2(i3)*vn2(i4)*vn2(i1)
                               v2 = v2 + vt1 + vt2
                            enddo
                         enddo
                         
                         v3 = vp2(i1)*vn2(i2) + vp2(i2)*vn2(i1)
                         vt1 = 0.d0
                         do i3 = 1,3
                            vt1 = vt1 + vp2(i3)*vn2(i3)
                         enddo
                         vt1 = vt1*(vn2(i1)*vn2(i2) - deltaij(i1,i2))
                         v3 = v3 + vt1
                         
                         v1 = 3.d0/(vpsibl_u2**6*rv2**3)
                         v4 = v4 + v1*(v2+(rv2/2.d0)*v3)
                      endif
                      
                      if ((i1.eq.1).and.(i2.eq.1)) then                            
                         u(H_A11)%d(i,j,k) = v4
                      else if ((i1.eq.1).and.(i2.eq.2)) then
                         u(H_A12)%d(i,j,k) = v4
                      else if ((i1.eq.1).and.(i2.eq.3)) then
                         u(H_A13)%d(i,j,k) = v4
                      else if ((i1.eq.2).and.(i2.eq.2)) then
                         u(H_A22)%d(i,j,k) = v4
                      else if ((i1.eq.2).and.(i2.eq.3)) then
                         u(H_A23)%d(i,j,k) = v4
                      else if ((i1.eq.3).and.(i2.eq.3)) then
                         u(H_A33)%d(i,j,k) = v4
                      endif
                   enddo
                enddo
               
! and add the EM initial data...

                if (i_btype.eq.1) then
                   u(H_EX)%d(i,j,k) = 0.d0
                   u(H_EY)%d(i,j,k) = 0.d0 
                   u(H_EZ)%d(i,j,k) = 0.d0
                   u(H_BX)%d(i,j,k) = 0.d0 
                   u(H_BY)%d(i,j,k) = 0.d0
                   u(H_BZ)%d(i,j,k) = v_bamp/((vpsibl)**6)+v_bamp*1.d-6
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 
                elseif (i_btype.eq.2) then
                   ! Electrically charged BH, exact for single non-boosted 
                   u(H_BX)%d(i,j,k) = 0.0
                   u(H_BY)%d(i,j,k) = 0.0
                   u(H_BZ)%d(i,j,k) = 0.0
                   ! radial electric field E^r = Q/(r^2*psi^6)
                   ! eq (3.16) from http://arxiv.org/abs/0907.1151
                   u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 

                   if (n_bhs.eq.2) then
                     u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                   end if
                elseif (i_btype.eq.3) then
                   ! Magnetic monopole BH, exact for single non-boosted 
                   u(H_EX)%d(i,j,k) = 0.0
                   u(H_EY)%d(i,j,k) = 0.0
                   u(H_EZ)%d(i,j,k) = 0.0
                   ! radial electric field E^r = Q/(r^2*psi^6)
                   ! eq (3.16) from http://arxiv.org/abs/0907.1151
                   u(H_BX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_BY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_BZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 

                   if (n_bhs.eq.2) then
                     u(H_BX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_BY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_BZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                   end if

                elseif (i_btype.eq.4) then
                   ! DYonic BH:
                   ! Magnetic monopole BH, exact for single non-boosted 
                   u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   ! radial electric field E^r = Q/(r^2*psi^6)
                   ! eq (3.16) from http://arxiv.org/abs/0907.1151
                   u(H_BX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_BY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_BZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 

                   if (n_bhs.eq.2) then
                     u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                     u(H_BX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_BY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_BZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                   end if
                elseif (i_btype.eq.5) then
                   ! electric BH with magnetic BH:
                   ! Electrically charged BH, exact for single non-boosted 
                   u(H_BX)%d(i,j,k) = 0.0
                   u(H_BY)%d(i,j,k) = 0.0
                   u(H_BZ)%d(i,j,k) = 0.0
                   ! radial electric field E^r = Q/(r^2*psi^6)
                   ! eq (3.16) from http://arxiv.org/abs/0907.1151
                   u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHIR)%d(i,j,k) = v_sfamp-v_e1amp*exp(-rv1**2)

                   if (n_bhs.eq.2) then
                     u(H_BX)%d(i,j,k) = v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_BY)%d(i,j,k) = v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_BZ)%d(i,j,k) = v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                     u(H_PHIR)%d(i,j,k) = u(H_PHIR)%d(i,j,k)+v_e2amp*exp(-rv2**2)
                   end if

                endif

! and add the SF initial data...

                if (i_btype.ne.5) then
                   u(H_PHIR)%d(i,j,k) = v_sfamp
                end if
                u(H_PHII)%d(i,j,k) = 0.0d0
                u(H_PIR)%d(i,j,k) = 0.0d0
                u(H_PII)%d(i,j,k) = 0.0d0

                
             enddo
          enddo
       enddo

    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 1) then
      !-----------------------------------------------------------
      ! Kerr-Schild data
      !-----------------------------------------------------------
        do k = 1, nz
          do j = 1, ny
             do i = 1, nx

                x1 = xx(i,j,k) - bh1x
                y1 = yy(i,j,k) - bh1y
                z1 = zz(i,j,k) - bh1z


                call ks_initial_data(gtd, Atd, chi, trK, Gamt, alpha, Betau, &
           &                 x1, y1, z1)

                u(H_ALPHA)%d(i,j,k)  =  alpha
                u(H_SHIFT1)%d(i,j,k) =  Betau(1)
                u(H_SHIFT2)%d(i,j,k) =  Betau(2)
                u(H_SHIFT3)%d(i,j,k) =  Betau(3)

                u(H_GB1)%d(i,j,k) =  0.0d0
                u(H_GB2)%d(i,j,k) =  0.0d0
                u(H_GB3)%d(i,j,k) =  0.0d0

                u(H_GAM1)%d(i,j,k) =  Gamt(1)
                u(H_GAM2)%d(i,j,k) =  Gamt(2)
                u(H_GAM3)%d(i,j,k) =  Gamt(3)

                u(H_GT11)%d(i,j,k) =  gtd(1,1)
                u(H_GT12)%d(i,j,k) =  gtd(1,2)
                u(H_GT13)%d(i,j,k) =  gtd(1,3)
                u(H_GT22)%d(i,j,k) =  gtd(2,2)
                u(H_GT23)%d(i,j,k) =  gtd(2,3)
                u(H_GT33)%d(i,j,k) =  gtd(3,3)

                u(H_A11)%d(i,j,k) =  Atd(1,1)
                u(H_A12)%d(i,j,k) =  Atd(1,2)
                u(H_A13)%d(i,j,k) =  Atd(1,3)
                u(H_A22)%d(i,j,k) =  Atd(2,2)
                u(H_A23)%d(i,j,k) =  Atd(2,3)
                u(H_A33)%d(i,j,k) =  Atd(3,3)

                u(H_CHI)%d(i,j,k) = chi
                u(H_TRK)%d(i,j,k) = trK

                u(H_EX)%d(i,j,k) = cos(x1)
                u(H_EY)%d(i,j,k) = sin(y1)
                u(H_EZ)%d(i,j,k) = cos(z1)
                u(H_BX)%d(i,j,k) = sin(x1)
                u(H_BY)%d(i,j,k) = cos(y1)
                u(H_BZ)%d(i,j,k) = sin(z1)
                u(H_PSI_EM)%d(i,j,k) = 0.d0 
                u(H_PHI_EM)%d(i,j,k) = 0.d0 

!               if (i .eq. 11 .and. j .eq. 11 .and. k .eq. 11) then
!                 print *,'I x1d = ',x1
!                 print *,'I y1d = ',y1
!                 print *,'I z1d = ',z1
!                 print *,'I gtd(11)  = ',gtd(1,1)
!                 print *,'I gtd(22)  = ',gtd(2,2)
!                 print *,'I Atd(11)  = ',Atd(1,1)
!                 print *,'I Atd(22)  = ',Atd(2,2)
!                 print *,'I chi      = ',chi
!                 print *,'I trK      = ',trK
!               end if


             enddo
          enddo
       enddo

    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 2) then
      !-----------------------------------------------------------
      ! FAKE DATA
      !-----------------------------------------------------------

       do k = 1, nz
          do j = 1, ny
             do i = 1, nx

                x1 = xx(i,j,k)
                y1 = yy(i,j,k)
                z1 = zz(i,j,k)

                u(H_ALPHA)%d(i,j,k)  =  1.0d0 - 0.25d0*sin(x1)
                u(H_SHIFT1)%d(i,j,k) =  0.0d0
                u(H_SHIFT2)%d(i,j,k) =  0.0d0
                u(H_SHIFT3)%d(i,j,k) =  0.0d0

      Gamt(1) = 5.D0*cos(x1)/(10.D0*sin(x1)+26.D0-1.D0*cos(x1)**2)
      Gamt(2) = -5.D0*sin(y1)/(25.D0+10.D0*cos(y1)+cos(y1)**2)
      Gamt(3) = 0.d0

                u(H_GAM1)%d(i,j,k) =  Gamt(1)
                u(H_GAM2)%d(i,j,k) =  Gamt(2)
                u(H_GAM3)%d(i,j,k) =  Gamt(3)

                u(H_GB1)%d(i,j,k) =  0.0d0
                u(H_GB2)%d(i,j,k) =  0.0d0
                u(H_GB3)%d(i,j,k) =  0.0d0

                u(H_GT11)%d(i,j,k) =  1.0d0+0.2d0*sin(x1)
                u(H_GT12)%d(i,j,k) =  0.0d0
                u(H_GT13)%d(i,j,k) =  0.0d0
                u(H_GT22)%d(i,j,k) =  1.0d0+0.2d0*cos(y1)
                u(H_GT23)%d(i,j,k) =  0.0d0
                u(H_GT33)%d(i,j,k) =  1.0d0/(u(H_GT11)%d(i,j,k) &
     &                    *u(H_GT22)%d(i,j,k));

      Atd(1,1) = exp(-0.4D1*cos(x1)*sin(y1))*(cos(x1)-0.3333333333D0*ex&
     &p(4*cos(x1)*sin(y1))*(1+0.2D0*sin(x1))*(5.D0*exp(-4.D0*cos(x1)*sin&
     &(y1))/(5.D0+sin(x1))*cos(x1)+5.D0*exp(-4.D0*cos(x1)*sin(y1))/(5.D0&
     &+cos(y1))*cos(y1)+0.4D-1*(25.D0+5.D0*cos(y1)+5.D0*sin(x1)+sin(x1)*&
     &cos(y1))*exp(-4.D0*cos(x1)*sin(y1))*cos(z1)))
      Atd(1,2) = 0.D0
      Atd(1,3) = 0.D0
      Atd(2,1) = 0.D0
      Atd(2,2) = exp(-0.4D1*cos(x1)*sin(y1))*(cos(y1)-0.3333333333D0*ex&
     &p(4*cos(x1)*sin(y1))*(1+0.2D0*cos(y1))*(5.D0*exp(-4.D0*cos(x1)*sin&
     &(y1))/(5.D0+sin(x1))*cos(x1)+5.D0*exp(-4.D0*cos(x1)*sin(y1))/(5.D0&
     &+cos(y1))*cos(y1)+0.4D-1*(25.D0+5.D0*cos(y1)+5.D0*sin(x1)+sin(x1)*&
     &cos(y1))*exp(-4.D0*cos(x1)*sin(y1))*cos(z1)))
      Atd(2,3) = 0.D0
      Atd(3,1) = 0.D0
      Atd(3,2) = 0.D0
      Atd(3,3) = exp(-0.4D1*cos(x1)*sin(y1))*(cos(z1)-0.3333333333D0*ex&
     &p(4*cos(x1)*sin(y1))/(1+0.2D0*sin(x1))/(1+0.2D0*cos(y1))*(5.D0*exp&
     &(-4.D0*cos(x1)*sin(y1))/(5.D0+sin(x1))*cos(x1)+5.D0*exp(-4.D0*cos(&
     &x1)*sin(y1))/(5.D0+cos(y1))*cos(y1)+0.4D-1*(25.D0+5.D0*cos(y1)+5.D&
     &0*sin(x1)+sin(x1)*cos(y1))*exp(-4.D0*cos(x1)*sin(y1))*cos(z1)))



                u(H_A11)%d(i,j,k) = Atd(1,1)
                u(H_A12)%d(i,j,k) = Atd(1,2)
                u(H_A13)%d(i,j,k) = Atd(1,3)
                u(H_A22)%d(i,j,k) = Atd(2,2)
                u(H_A23)%d(i,j,k) = Atd(2,3)
                u(H_A33)%d(i,j,k) = Atd(3,3)

                u(H_CHI)%d(i,j,k) = exp(-4.0d0*cos(x1)*sin(y1))
                u(H_TRK)%d(i,j,k) = 5.D0*exp(-4.D0*cos(x1)*sin(y1))&
     &/(5.D0+sin(x1))*cos(x1)+5.D0*e&
     &xp(-4.D0*cos(x1)*sin(y1))/(5.D0+cos(y1))*cos(y1)+0.4D-1*(25.D0+5.D&
     &0*cos(y1)+5.D0*sin(x1)+sin(x1)*cos(y1))*exp(-4.D0*cos(x1)*sin(y1))&
     &*cos(z1)

                if (i_btype.eq.1) then
                   u(H_EX)%d(i,j,k) = cos(x1)
                   u(H_EY)%d(i,j,k) = sin(y1)
                   u(H_EZ)%d(i,j,k) = cos(z1)
                   u(H_BX)%d(i,j,k) = sin(x1)
                   u(H_BY)%d(i,j,k) = cos(y1)
                   u(H_BZ)%d(i,j,k) = sin(z1)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 
                endif
 

             enddo
          enddo
       enddo


    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 3) then
      !-----------------------------------------------------------
      ! noise data
      !-----------------------------------------------------------
        rb_amp = par(P_MR_AMP)
        do k = 1, nz
          do j = 1, ny
             do i = 1, nx

                call random_number(harvest)
                harvest = 2.0*(harvest - 0.5)
                
                u(H_ALPHA)%d(i,j,k)  = 1.0d0 + rb_amp*harvest(1)
                u(H_SHIFT1)%d(i,j,k) = rb_amp*harvest(31) 
                u(H_SHIFT2)%d(i,j,k) = rb_amp*harvest(32)
                u(H_SHIFT3)%d(i,j,k) = rb_amp*harvest(33)

                u(H_GB1)%d(i,j,k) = rb_amp*harvest(2)
                u(H_GB2)%d(i,j,k) = rb_amp*harvest(3)
                u(H_GB3)%d(i,j,k) = rb_amp*harvest(4)

                u(H_GAM1)%d(i,j,k) = rb_amp*harvest(5) 
                u(H_GAM2)%d(i,j,k) = rb_amp*harvest(6)
                u(H_GAM3)%d(i,j,k) = rb_amp*harvest(7)

                u(H_GT11)%d(i,j,k) =  1.0d0 + rb_amp*harvest(8)
                u(H_GT12)%d(i,j,k) =  0.0d0 + rb_amp*harvest(9)
                u(H_GT13)%d(i,j,k) =  0.0d0 + rb_amp*harvest(10)
                u(H_GT22)%d(i,j,k) =  1.0d0 + rb_amp*harvest(11)
                u(H_GT23)%d(i,j,k) =  0.0d0 + rb_amp*harvest(12)
                u(H_GT33)%d(i,j,k) =  1.0d0 + rb_amp*harvest(13)

                u(H_A11)%d(i,j,k) =  rb_amp*harvest(14)
                u(H_A12)%d(i,j,k) =  rb_amp*harvest(15)
                u(H_A13)%d(i,j,k) =  rb_amp*harvest(16)
                u(H_A22)%d(i,j,k) =  rb_amp*harvest(17)
                u(H_A23)%d(i,j,k) =  rb_amp*harvest(18)
                u(H_A33)%d(i,j,k) =  rb_amp*harvest(19)

                u(H_CHI)%d(i,j,k) = 1.0d0 + rb_amp*harvest(20)
                u(H_TRK)%d(i,j,k) = rb_amp*harvest(21)

                u(H_EX)%d(i,j,k) = rb_amp*harvest(22)
                u(H_EY)%d(i,j,k) = rb_amp*harvest(23)
                u(H_EZ)%d(i,j,k) = rb_amp*harvest(24)
                u(H_BX)%d(i,j,k) = rb_amp*harvest(25)
                u(H_BY)%d(i,j,k) = rb_amp*harvest(26)
                u(H_BZ)%d(i,j,k) = rb_amp*harvest(27)
                u(H_PSI_EM)%d(i,j,k) = rb_amp*harvest(28)
                u(H_PHI_EM)%d(i,j,k) = rb_amp*harvest(29)

             enddo
          enddo
       enddo

    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 4) then
      !-----------------------------------------------------------
      ! Kerr-Schild Kerr black hole data
      !-----------------------------------------------------------
        do k = 1, nz
          do j = 1, ny
             do i = 1, nx

                x1 = xx(i,j,k) - bh1x
                y1 = yy(i,j,k) - bh1y
                z1 = zz(i,j,k) - bh1z


                call kskerr_initial_data(gtd, Atd, chi, trK, Gamt, &
                     & alpha, Betau, x1, y1, z1, mass1, spin1)

                u(H_ALPHA)%d(i,j,k)  =  alpha
                u(H_SHIFT1)%d(i,j,k) =  Betau(1)
                u(H_SHIFT2)%d(i,j,k) =  Betau(2)
                u(H_SHIFT3)%d(i,j,k) =  Betau(3)

                u(H_GB1)%d(i,j,k) =  0.0d0
                u(H_GB2)%d(i,j,k) =  0.0d0
                u(H_GB3)%d(i,j,k) =  0.0d0

                u(H_GAM1)%d(i,j,k) =  Gamt(1)
                u(H_GAM2)%d(i,j,k) =  Gamt(2)
                u(H_GAM3)%d(i,j,k) =  Gamt(3)

                u(H_GT11)%d(i,j,k) =  gtd(1,1)
                u(H_GT12)%d(i,j,k) =  gtd(1,2)
                u(H_GT13)%d(i,j,k) =  gtd(1,3)
                u(H_GT22)%d(i,j,k) =  gtd(2,2)
                u(H_GT23)%d(i,j,k) =  gtd(2,3)
                u(H_GT33)%d(i,j,k) =  gtd(3,3)

                u(H_A11)%d(i,j,k) =  Atd(1,1)
                u(H_A12)%d(i,j,k) =  Atd(1,2)
                u(H_A13)%d(i,j,k) =  Atd(1,3)
                u(H_A22)%d(i,j,k) =  Atd(2,2)
                u(H_A23)%d(i,j,k) =  Atd(2,3)
                u(H_A33)%d(i,j,k) =  Atd(3,3)

                u(H_CHI)%d(i,j,k) = chi
                u(H_TRK)%d(i,j,k) = trK

                u(H_EX)%d(i,j,k) = 0.d0
                u(H_EY)%d(i,j,k) = 0.d0
                u(H_EZ)%d(i,j,k) = 0.d0
                u(H_BX)%d(i,j,k) = 0.d0
                u(H_BY)%d(i,j,k) = 0.d0
                u(H_BZ)%d(i,j,k) = 0.d0
                u(H_PSI_EM)%d(i,j,k) = 0.d0 
                u(H_PHI_EM)%d(i,j,k) = 0.d0 


             enddo
          enddo
       enddo


    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 5) then
      !-----------------------------------------------------------
      ! 4 black hole data
      !-----------------------------------------------------------
      ! assume that geometry is set elsewhere, as in hyperGR.  Do nothing.
!       print *,'travis here, lets start setting things up - need coordinates'
!     carlos: formulas from eq (66) for psi = 1 + m1/(2*r1) + m2/(2*r2)
!     and (70) for the extrinsic curvature Atilde_ij = psi^{-6} * Abar_ij
!     in  http://relativity.livingreviews.org/Articles/lrr-2000-5/

       do k=1,nz
          do j=1,ny
             do i=1,nx
                
                x1 = xx(i,j,k) - bh1x
                y1 = yy(i,j,k) - bh1y
                z1 = zz(i,j,k) - bh1z
                
                rv1 = sqrt(x1**2+y1**2+z1**2)
                if (rv1.lt.0.0001d0) then
                   rv1=0.0001d0
                endif
                vn1(1) = x1/rv1
                vn1(2) = y1/rv1
                vn1(3) = z1/rv1
                
                if (n_bhs.ge.2) then
                   x2 = xx(i,j,k) - bh2x
                   y2 = yy(i,j,k) - bh2y
                   z2 = zz(i,j,k) - bh2z
                   
                   rv2 = sqrt(x2**2+y2**2+z2**2)
                   if (rv2.lt.0.0001d0) then
                      rv2=0.0001d0
                   endif
                   vn2(1) = x2/rv2
                   vn2(2) = y2/rv2
                   vn2(3) = z2/rv2
                endif

                if (n_bhs.ge.3) then
                   x3 = xx(i,j,k) - bh3x
                   y3 = yy(i,j,k) - bh3y
                   z3 = zz(i,j,k) - bh3z
                   
                   rv3 = sqrt(x3**2+y3**2+z3**2)
                   if (rv3.lt.0.0001d0) then
                      rv3=0.0001d0
                   endif
                   vn3(1) = x3/rv3
                   vn3(2) = y3/rv3
                   vn3(3) = z3/rv3
                endif

                if (n_bhs.eq.4) then
                   x4 = xx(i,j,k) - bh4x
                   y4 = yy(i,j,k) - bh4y
                   z4 = zz(i,j,k) - bh4z
                   
                   rv4 = sqrt(x4**2+y4**2+z4**2)
                   if (rv4.lt.0.0001d0) then
                      rv4=0.0001d0
                   endif
                   vn4(1) = x4/rv4
                   vn4(2) = y4/rv4
                   vn4(3) = z4/rv4
                endif


                ! Brill-Lindquist conformal factor:
                vpsibl = 1.d0 + mass1/(2.d0*rv1)
                if (n_bhs.ge.2) then
                   vpsibl = vpsibl + mass2/(2.d0*rv2)
                endif
                if (n_bhs.ge.3) then
                   vpsibl = vpsibl + mass3/(2.d0*rv3)
                endif
                if (n_bhs.eq.4) then
                   vpsibl = vpsibl + mass4/(2.d0*rv4)
                endif

! hard-code in for the time being:
                i_better_u = 1
                if (i_better_u.eq.1) then
                   v_u_corr = 0.d0
                   
                   if (abs(spin1).gt.0.000001) then
                      amp_capj = 4.d0*spin1/mass1**2
                      amp_capr = 2.d0*rv1/mass1
                      l_r = 1.d0/(1.d0 + amp_capr)
                      u0_j = (l_r + l_r**2 + l_r**3 - 4.d0*l_r**4 + 2.d0*l_r**5)/40.d0
                      u2_j = -l_r**5/20.d0
                      mu_j =         vn1(1)*vs1(1)
                      mu_j =  mu_j + vn1(2)*vs1(2)
                      mu_j = (mu_j + vn1(3)*vs1(3))/abs(spin1)
                      p2_mu_j = (3.d0*mu_j**2 - 1.d0)/2.d0
                      v_u_j1 = amp_capj**2*(u0_j+u2_j*amp_capr**2*p2_mu_j)
                      v_u_corr = v_u_corr + v_u_j1
                   endif
                   
                   if (vp1tot.gt.0.000001) then
                      amp_capp = 2.d0*vp1tot/mass1
                      amp_capr = 2.d0*rv1/mass1
                      l_r = 1.d0/(1.d0 + amp_capr)
                      u0_p = l_r - 2.d0*l_r**2 + 2.d0*l_r**3 
                      u0_p = (u0_p - l_r**4 + 0.2d0*l_r**5)*(5.d0/32.d0)
                      u2_p = 15.d0*l_r + 132.d0*l_r**2 + 53.d0*l_r**3
                      u2_p = u2_p + 96.d0*l_r**4 + 82.d0*l_r**5
                      u2_p = u2_p + (84.d0/amp_capr)*(l_r**5+log(l_r)/amp_capr)
                      u2_p = (u2_p)/(80.d0*amp_capr)
                      mu_p =        vn1(1)*vp1(1)/vp1tot
                      mu_p = mu_p + vn1(2)*vp1(2)/vp1tot
                      mu_p = mu_p + vn1(3)*vp1(3)/vp1tot
                      p2_mu_p = (3.d0*mu_p**2 - 1.d0)/2.d0
                      v_u_p1 = amp_capp**2*(u0_p+u2_p*p2_mu_p)
                      v_u_corr = v_u_corr + v_u_p1
                   endif
                   
                   if ((vp1tot.gt.0.000001).and.(abs(spin1).gt.0.000001)) then
                      v1 =      (vp1(2)*vs1(3)-vp1(3)*vs1(2))*vn1(1)
                      v1 = v1 + (vp1(3)*vs1(1)-vp1(1)*vs1(3))*vn1(2)
                      v1 = v1 + (vp1(1)*vs1(2)-vp1(2)*vs1(1))*vn1(3)
                      v1 = v1*(16.d0/mass1**4)*rv1
                      
                      amp_capr = 2.d0*rv1/mass1
                      l_r = 1.d0/(1.d0 + amp_capr)
                      
                      v2 = 1.d0 + 5.d0*amp_capr + 10.d0*amp_capr**2
                      
                      v_u_c1 = (v1*v2*l_r**5)/80.d0
                      v_u_corr = v_u_corr + v_u_c1
                   endif
                   
                   if (n_bhs.ge.2) then
                      
                      if (abs(spin2).gt.0.000001) then
                         amp_capj = 4.d0*spin2/mass2**2
                         amp_capr = 2.d0*rv2/mass2
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_j = (l_r + l_r**2 + l_r**3 - 4.d0*l_r**4 + 2.d0*l_r**5)/40.d0
                         u2_j = -l_r**5/20.d0
                         mu_j =         vn2(1)*vs2(1)
                         mu_j =  mu_j + vn2(2)*vs2(2)
                         mu_j = (mu_j + vn2(3)*vs2(3))/abs(spin2)
                         p2_mu_j = (3.d0*mu_j**2 - 1.d0)/2.d0
                         v_u_j2 = amp_capj**2*(u0_j+u2_j*amp_capr**2*p2_mu_j)
                         v_u_corr = v_u_corr + v_u_j2
                      endif
                      
                      if (vp2tot.gt.0.000001) then
                         amp_capp = 2.d0*vp2tot/mass2
                         amp_capr = 2.d0*rv2/mass2
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_p = l_r - 2.d0*l_r**2 + 2.d0*l_r**3 
                         u0_p = (u0_p - l_r**4 + 0.2d0*l_r**5)*(5.d0/32.d0)
                         u2_p = 15.d0*l_r + 132.d0*l_r**2 + 53.d0*l_r**3
                         u2_p = u2_p + 96.d0*l_r**4 + 82.d0*l_r**5
                         u2_p = u2_p + (84.d0/amp_capr)*(l_r**5+log(l_r)/amp_capr)
                         u2_p = (u2_p)/(80.d0*amp_capr)
                         mu_p =        vn2(1)*vp2(1)/vp2tot
                         mu_p = mu_p + vn2(2)*vp2(2)/vp2tot
                         mu_p = mu_p + vn2(3)*vp2(3)/vp2tot
                         p2_mu_p = (3.d0*mu_p**2 - 1.d0)/2.d0
                         v_u_p2 = amp_capp**2*(u0_p+u2_p*p2_mu_p)
                         v_u_corr = v_u_corr + v_u_p2
                      endif
                      
                      if ((vp2tot.gt.0.000001).and.(abs(spin2).gt.0.000001)) then
                         v1 =      (vp2(2)*vs2(3)-vp2(3)*vs2(2))*vn2(1)
                         v1 = v1 + (vp2(3)*vs2(1)-vp2(1)*vs2(3))*vn2(2)
                         v1 = v1 + (vp2(1)*vs2(2)-vp2(2)*vs2(1))*vn2(3)
                         v1 = v1*(16.d0/mass2**4)*rv2
                         
                         amp_capr = 2.d0*rv2/mass2
                         l_r = 1.d0/(1.d0 + amp_capr)
                         
                         v2 = 1.d0 + 5.d0*amp_capr + 10.d0*amp_capr**2
                         
                         v_u_c2 = (v1*v2*l_r**5)/80.d0
                         v_u_corr = v_u_corr + v_u_c2
                      endif

                   endif

                   if (n_bhs.ge.3) then
                      
                      if (abs(spin3).gt.0.000001) then
                         amp_capj = 4.d0*spin3/mass3**2
                         amp_capr = 2.d0*rv3/mass3
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_j = (l_r + l_r**2 + l_r**3 - 4.d0*l_r**4 + 2.d0*l_r**5)/40.d0
                         u2_j = -l_r**5/20.d0
                         mu_j =         vn3(1)*vs3(1)
                         mu_j =  mu_j + vn3(2)*vs3(2)
                         mu_j = (mu_j + vn3(3)*vs3(3))/abs(spin3)
                         p2_mu_j = (3.d0*mu_j**2 - 1.d0)/2.d0
                         v_u_j3 = amp_capj**2*(u0_j+u2_j*amp_capr**2*p2_mu_j)
                         v_u_corr = v_u_corr + v_u_j3
                      endif
                      
                      if (vp3tot.gt.0.000001) then
                         amp_capp = 2.d0*vp3tot/mass3
                         amp_capr = 2.d0*rv3/mass3
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_p = l_r - 2.d0*l_r**2 + 2.d0*l_r**3 
                         u0_p = (u0_p - l_r**4 + 0.2d0*l_r**5)*(5.d0/32.d0)
                         u2_p = 15.d0*l_r + 132.d0*l_r**2 + 53.d0*l_r**3
                         u2_p = u2_p + 96.d0*l_r**4 + 82.d0*l_r**5
                         u2_p = u2_p + (84.d0/amp_capr)*(l_r**5+log(l_r)/amp_capr)
                         u2_p = (u2_p)/(80.d0*amp_capr)
                         mu_p =        vn3(1)*vp3(1)/vp3tot
                         mu_p = mu_p + vn3(2)*vp3(2)/vp3tot
                         mu_p = mu_p + vn3(3)*vp3(3)/vp3tot
                         p2_mu_p = (3.d0*mu_p**2 - 1.d0)/2.d0
                         v_u_p3 = amp_capp**2*(u0_p+u2_p*p2_mu_p)
                         v_u_corr = v_u_corr + v_u_p3
                      endif
                      
                      if ((vp3tot.gt.0.000001).and.(abs(spin3).gt.0.000001)) then
                         v1 =      (vp3(2)*vs3(3)-vp3(3)*vs3(2))*vn3(1)
                         v1 = v1 + (vp3(3)*vs3(1)-vp3(1)*vs3(3))*vn3(2)
                         v1 = v1 + (vp3(1)*vs3(2)-vp3(2)*vs3(1))*vn3(3)
                         v1 = v1*(16.d0/mass3**4)*rv3
                         
                         amp_capr = 2.d0*rv3/mass3
                         l_r = 1.d0/(1.d0 + amp_capr)
                         
                         v2 = 1.d0 + 5.d0*amp_capr + 10.d0*amp_capr**2
                         
                         v_u_c3 = (v1*v2*l_r**5)/80.d0
                         v_u_corr = v_u_corr + v_u_c3
                      endif

                   endif

                   if (n_bhs.eq.4) then
                      
                      if (abs(spin4).gt.0.000001) then
                         amp_capj = 4.d0*spin4/mass4**2
                         amp_capr = 2.d0*rv4/mass4
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_j = (l_r + l_r**2 + l_r**3 - 4.d0*l_r**4 + 2.d0*l_r**5)/40.d0
                         u2_j = -l_r**5/20.d0
                         mu_j =         vn4(1)*vs4(1)
                         mu_j =  mu_j + vn4(2)*vs4(2)
                         mu_j = (mu_j + vn4(3)*vs4(3))/abs(spin4)
                         p2_mu_j = (3.d0*mu_j**2 - 1.d0)/2.d0
                         v_u_j4 = amp_capj**2*(u0_j+u2_j*amp_capr**2*p2_mu_j)
                         v_u_corr = v_u_corr + v_u_j4
                      endif
                      
                      if (vp4tot.gt.0.000001) then
                         amp_capp = 2.d0*vp4tot/mass4
                         amp_capr = 2.d0*rv4/mass4
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_p = l_r - 2.d0*l_r**2 + 2.d0*l_r**3 
                         u0_p = (u0_p - l_r**4 + 0.2d0*l_r**5)*(5.d0/32.d0)
                         u2_p = 15.d0*l_r + 132.d0*l_r**2 + 53.d0*l_r**3
                         u2_p = u2_p + 96.d0*l_r**4 + 82.d0*l_r**5
                         u2_p = u2_p + (84.d0/amp_capr)*(l_r**5+log(l_r)/amp_capr)
                         u2_p = (u2_p)/(80.d0*amp_capr)
                         mu_p =        vn4(1)*vp4(1)/vp4tot
                         mu_p = mu_p + vn4(2)*vp4(2)/vp4tot
                         mu_p = mu_p + vn4(3)*vp4(3)/vp4tot
                         p2_mu_p = (3.d0*mu_p**2 - 1.d0)/2.d0
                         v_u_p4 = amp_capp**2*(u0_p+u2_p*p2_mu_p)
                         v_u_corr = v_u_corr + v_u_p4
                      endif
                      
                      if ((vp4tot.gt.0.000001).and.(abs(spin4).gt.0.000001)) then
                         v1 =      (vp4(2)*vs4(3)-vp4(3)*vs4(2))*vn4(1)
                         v1 = v1 + (vp4(3)*vs4(1)-vp4(1)*vs4(3))*vn4(2)
                         v1 = v1 + (vp4(1)*vs4(2)-vp4(2)*vs4(1))*vn4(3)
                         v1 = v1*(16.d0/mass4**4)*rv4
                         
                         amp_capr = 2.d0*rv4/mass4
                         l_r = 1.d0/(1.d0 + amp_capr)
                         
                         v2 = 1.d0 + 5.d0*amp_capr + 10.d0*amp_capr**2
                         
                         v_u_c4 = (v1*v2*l_r**5)/80.d0
                         v_u_corr = v_u_corr + v_u_c4
                      endif

                   endif
                  
! vpsibl_u will be used for the conformal factor, 
                   vpsibl_u  = vpsibl + v_u_corr
! vpsibl_u2 is for the Aij terms...
! since the corrections are first order...
! adding half of the correction seems to give the best results...
! update - do a fit for spin = 0.6...
!                   vpsibl_u2 = vpsibl
!                   vpsibl_u2 = vpsibl + v_u_corr/2.d0
!                   vpsibl_u2 = vpsibl + 0.313d0*v_u_corr
                   vpsibl_u2 = vpsibl + v_u_corr
                else
                   vpsibl_u  = vpsibl
                   vpsibl_u2 = vpsibl
                endif
                
!                u(H_ALPHA)%d(i,j,k) = 1.d0/vpsibl**2
!                v2 =  1.d0/vpsibl**4
                u(H_ALPHA)%d(i,j,k) = 1.d0/vpsibl_u**2
                v2 =  1.d0/vpsibl_u**4
!                u(H_PHI)%d(i,j,k) = log( vpsibl )
!                u(H_CHI)%d(i,j,k) = exp( -4.d0 * log( vpsibl ) )
                if (setchi) then
                   u(H_CHI)%d(i,j,k) = v2
                end if
                u(H_TRK)%d(i,j,k) =  0.0d0
                
                u(H_SHIFT1)%d(i,j,k) =  0.0d0
                u(H_SHIFT2)%d(i,j,k) =  0.0d0
                u(H_SHIFT3)%d(i,j,k) =  0.0d0
                   
                u(H_GAM1)%d(i,j,k) =  0.0d0
                u(H_GAM2)%d(i,j,k) =  0.0d0
                u(H_GAM3)%d(i,j,k) =  0.0d0
                
                u(H_GB1)%d(i,j,k) =  0.0d0
                u(H_GB2)%d(i,j,k) =  0.0d0
                u(H_GB3)%d(i,j,k) =  0.0d0
                
                u(H_GT11)%d(i,j,k) =  1.0d0
                u(H_GT12)%d(i,j,k) =  0.0d0
                u(H_GT13)%d(i,j,k) =  0.0d0
                u(H_GT22)%d(i,j,k) =  1.0d0
                u(H_GT23)%d(i,j,k) =  0.0d0
                u(H_GT33)%d(i,j,k) =  1.0d0
                

                do i1 = 1,3
                   do i2 = i1,3
                      
                      v2 = 0.d0
                      do i3 = 1,3
                         do i4 = 1,3
                            vt1 = epijk(i1,i3,i4)*vs1(i3)*vn1(i4)*vn1(i2)
                            vt2 = epijk(i2,i3,i4)*vs1(i3)*vn1(i4)*vn1(i1)
                            v2 = v2 + vt1 + vt2
                         enddo
                      enddo
                      
                      v3 = vp1(i1)*vn1(i2) + vp1(i2)*vn1(i1)
                      vt1 = 0.d0
                      do i3 = 1,3
                         vt1 = vt1 + vp1(i3)*vn1(i3)
                      enddo
                      ! travis, please check the following,
                      ! there is a sign wrong
                      !vt1 = vt1*(deltaij(i1,i2)-vn1(i1)*vn1(i2))
                      vt1 = vt1*(vn1(i1)*vn1(i2) - deltaij(i1,i2))
                      v3 = v3 + vt1
                      
!                      v1 = 3.d0/(vpsibl**6*rv1**3)
                      v1 = 3.d0/(vpsibl_u2**6*rv1**3)
                      v4 = v1*(v2+(rv1/2.d0)*v3)
                      
                      if (n_bhs.ge.2) then
                         v2 = 0.d0
                         do i3 = 1,3
                            do i4 = 1,3
                               vt1 = epijk(i1,i3,i4)*vs2(i3)*vn2(i4)*vn2(i2)
                               vt2 = epijk(i2,i3,i4)*vs2(i3)*vn2(i4)*vn2(i1)
                               v2 = v2 + vt1 + vt2
                            enddo
                         enddo
                         
                         v3 = vp2(i1)*vn2(i2) + vp2(i2)*vn2(i1)
                         vt1 = 0.d0
                         do i3 = 1,3
                            vt1 = vt1 + vp2(i3)*vn2(i3)
                         enddo
                         ! the same wrong sign
                         !vt1 = vt1*(deltaij(i1,i2)-vn2(i1)*vn2(i2))
                         vt1 = vt1*(vn2(i1)*vn2(i2) - deltaij(i1,i2))
                         v3 = v3 + vt1
                         
!                         v1 = 3.d0/(vpsibl**6*rv2**3)
                         v1 = 3.d0/(vpsibl_u2**6*rv2**3)
                         v4 = v4 + v1*(v2+(rv2/2.d0)*v3)
                      endif

                      if (n_bhs.ge.3) then
                         v2 = 0.d0
                         do i3 = 1,3
                            do i4 = 1,3
                               vt1 = epijk(i1,i3,i4)*vs3(i3)*vn3(i4)*vn3(i2)
                               vt2 = epijk(i2,i3,i4)*vs3(i3)*vn3(i4)*vn3(i1)
                               v2 = v2 + vt1 + vt2
                            enddo
                         enddo
                         
                         v3 = vp3(i1)*vn3(i2) + vp3(i2)*vn3(i1)
                         vt1 = 0.d0
                         do i3 = 1,3
                            vt1 = vt1 + vp3(i3)*vn3(i3)
                         enddo
                         ! the same wrong sign
                         !vt1 = vt1*(deltaij(i1,i2)-vn2(i1)*vn2(i2))
                         vt1 = vt1*(vn3(i1)*vn3(i2) - deltaij(i1,i2))
                         v3 = v3 + vt1
                         
!                         v1 = 3.d0/(vpsibl**6*rv2**3)
                         v1 = 3.d0/(vpsibl_u2**6*rv3**3)
                         v4 = v4 + v1*(v2+(rv3/2.d0)*v3)
                      endif

                      if (n_bhs.eq.4) then
                         v2 = 0.d0
                         do i3 = 1,3
                            do i4 = 1,3
                               vt1 = epijk(i1,i3,i4)*vs4(i3)*vn4(i4)*vn4(i2)
                               vt2 = epijk(i2,i3,i4)*vs4(i3)*vn4(i4)*vn4(i1)
                               v2 = v2 + vt1 + vt2
                            enddo
                         enddo
                         
                         v3 = vp4(i1)*vn4(i2) + vp4(i2)*vn4(i1)
                         vt1 = 0.d0
                         do i3 = 1,3
                            vt1 = vt1 + vp4(i3)*vn4(i3)
                         enddo
                         ! the same wrong sign
                         !vt1 = vt1*(deltaij(i1,i2)-vn2(i1)*vn2(i2))
                         vt1 = vt1*(vn4(i1)*vn4(i2) - deltaij(i1,i2))
                         v3 = v3 + vt1
                         
!                         v1 = 3.d0/(vpsibl**6*rv2**3)
                         v1 = 3.d0/(vpsibl_u2**6*rv4**3)
                         v4 = v4 + v1*(v2+(rv4/2.d0)*v3)
                      endif
                      
                      if ((i1.eq.1).and.(i2.eq.1)) then                            
                         u(H_A11)%d(i,j,k) = v4
                      else if ((i1.eq.1).and.(i2.eq.2)) then
                         u(H_A12)%d(i,j,k) = v4
                      else if ((i1.eq.1).and.(i2.eq.3)) then
                         u(H_A13)%d(i,j,k) = v4
                      else if ((i1.eq.2).and.(i2.eq.2)) then
                         u(H_A22)%d(i,j,k) = v4
                      else if ((i1.eq.2).and.(i2.eq.3)) then
                         u(H_A23)%d(i,j,k) = v4
                      else if ((i1.eq.3).and.(i2.eq.3)) then
                         u(H_A33)%d(i,j,k) = v4
                      endif
                   enddo
                enddo
               
! and add the EM initial data...
                if (i_btype.eq.1) then
                   u(H_EX)%d(i,j,k) = 0.d0
                   !u(H_EY)%d(i,j,k) = - 2.d0 * v_bamp/((vpsibl)**6)
                   u(H_EY)%d(i,j,k) = 0.d0 
                   u(H_EZ)%d(i,j,k) = 0.d0
                   !u(H_BX)%d(i,j,k) = v_bamp * exp( - rv1**2 )  
                   u(H_BX)%d(i,j,k) = 0.d0 
                   u(H_BY)%d(i,j,k) = 0.d0
                   u(H_BZ)%d(i,j,k) = v_bamp/((vpsibl)**6)+v_bamp*1.d-6
                   !u(H_PSI_EM)%d(i,j,k) = v_bamp/((vpsibl)**6)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   !u(H_PHI_EM)%d(i,j,k) = v_bamp/((vpsibl)**6)
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 
                endif
                
             enddo
          enddo
       enddo

    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 10) then
      !-----------------------------------------------------------
      ! 1 boson star
      !-----------------------------------------------------------

        xcenter = par(P_BS_1_X0)
        ycenter = par(P_BS_1_Y0)
        zcenter = par(P_BS_1_Z0)
        omega1  = par(P_BS_1_OMEGA)       
        vx      = par(P_BS_1_VX)
        vy      = par(P_BS_1_VY)
        vz      = par(P_BS_1_VZ)
                  
        lev    = nint(par(P_READ_DATA_LEVEL))
        interp = nint(par(P_INTERP_ID))
     
        call sdf_file_mem("psi1d.sdf", lev, shp)

        allocate(fdata1d(shp(1)),coords1d(shp(1)),STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for arrays to read sdf files'
          stop
        end if

        ! it is most useful to have 1-d coordinate arrays.
        allocate (x1d(nx), y1d(ny), z1d(nz), STAT=istat)
        if (istat .ne. 0) then
          write(*,*)'DEF_INITIAL:  Could not allocate memory'
          write(*,*)'for 1d coordinate arrays'
          stop
        end if
        do i = 1, nx
          x1d(i) = xx(i,1,1)
        end do
        do i = 1, ny
          y1d(i) = yy(1,i,1)
        end do
        do i = 1, nz
          z1d(i) = zz(1,1,i)
        end do

       minx0 = par(P_MINX0)
       miny0 = par(P_MINY0)
       minz0 = par(P_MINZ0)

       hi   = x1d(2) - x1d(1)
       !
       ! Need minimums for this grid alone:
       !
       minx = x1d(1)
       miny = y1d(1)
       minz = z1d(1)
       dnr  = shp(1)
      
       call read_sdf_file1d("psi1d.sdf", u(H_TRK)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       u(H_CHI)%d  = (u(H_TRK)%d)**(-4.0)
       u(H_GT11)%d =  1.0d0
       u(H_GT12)%d =  0.0d0
       u(H_GT13)%d =  0.0d0
       u(H_GT22)%d =  1.0d0
       u(H_GT23)%d =  0.0d0
       u(H_GT33)%d =  1.0d0

      ! HOW TO COMPUTE THE GAM1? is it 0??
       u(H_GAM1)%d =  0.0d0
       u(H_GAM2)%d =  0.0d0
       u(H_GAM3)%d =  0.0d0

       u(H_TRK)%d = 0.0d0
       u(H_A11)%d = 0.0d0
       u(H_A12)%d = 0.0d0
       u(H_A13)%d=  0.0d0
       u(H_A22)%d = 0.0d0
       u(H_A23)%d = 0.0d0
       u(H_A33)%d = 0.0d0
      
       call read_sdf_file1d("alpha1d.sdf", u(H_ALPHA)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)

       u(H_SHIFT1)%d =  0.0d0
       u(H_SHIFT2)%d =  0.0d0
       u(H_SHIFT3)%d =  0.0d0
       u(H_GB1)%d =  0.0d0
       u(H_GB2)%d =  0.0d0
       u(H_GB3)%d =  0.0d0


       call read_sdf_file1d("phir1d.sdf", u(H_PHIR)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
                   coords1d, dnr, interp, xcenter, ycenter, zcenter)

!       u(H_PHIR)%d = (u(H_PHIR)%d)/sqrt(2.0)
       u(H_PHII)%d = 0.0d0
       u(H_PII)%d  = -(omega1*u(H_PHIR)%d)/(u(H_ALPHA)%d)
       u(H_PIR)%d  = 0.0d0

! and add the EM initial data...

       u(H_EX)%d = 0.d0
       u(H_EY)%d = 0.d0 
       u(H_EZ)%d = 0.d0
       u(H_BX)%d = 0.d0 
       u(H_BY)%d = 0.d0
       u(H_BZ)%d = 0.0d0
       u(H_PSI_EM)%d = 0.d0 
       u(H_PHI_EM)%d = 0.d0 


      !-----------------------------------------------------------
      ! 2 boson star
      !-----------------------------------------------------------
       if (n_bhs.eq.2) then

         xcenter = par(P_BS_2_X0)
         ycenter = par(P_BS_2_Y0)
         zcenter = par(P_BS_2_Z0)
         omega2  = par(P_BS_2_OMEGA)       
         vx      = par(P_BS_2_VX)
         vy      = par(P_BS_2_VY)
         vz      = par(P_BS_2_VZ)

         call read_sdf_file1d("psi1d.sdf", u(H_TRK)%d,  &
                     x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
	  	     coords1d, dnr, interp, xcenter, ycenter, zcenter)

         u(H_CHI)%d  = u(H_CHI)%d + (u(H_TRK)%d)**(-4.0) - 1.0d0
         u(H_GT11)%d = 1.0d0
         u(H_GT12)%d = 0.0d0
         u(H_GT13)%d = 0.0d0
         u(H_GT22)%d = 1.0d0
         u(H_GT23)%d = 0.0d0
         u(H_GT33)%d = 1.0d0

         ! HOW TO COMPUTE THE GAM1? is it 0??
         u(H_GAM1)%d = 0.0d0
         u(H_GAM2)%d = 0.0d0
         u(H_GAM3)%d = 0.0d0

         u(H_TRK)%d = 0.0d0
         u(H_A11)%d = 0.0d0
         u(H_A12)%d = 0.0d0
         u(H_A13)%d = 0.0d0
         u(H_A22)%d = 0.0d0
         u(H_A23)%d = 0.0d0
         u(H_A33)%d = 0.0d0

         call read_sdf_file1d("alpha1d.sdf", u(H_TRK)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

         u(H_ALPHA)%d  = u(H_ALPHA)%d + u(H_TRK)%d - 1.0d0
         u(H_SHIFT1)%d = 0.0d0
         u(H_SHIFT2)%d = 0.0d0
         u(H_SHIFT3)%d = 0.0d0
         u(H_GB1)%d = 0.0d0
         u(H_GB2)%d = 0.0d0
         u(H_GB3)%d = 0.0d0

         call read_sdf_file1d("phir1d.sdf", u(H_PIR)%d,  &
                   x1d, y1d, z1d, nx, ny, nz, lev, fdata1d,  &
		   coords1d, dnr, interp, xcenter, ycenter, zcenter)

         u(H_PHIR)%d = u(H_PHIR)%d + u(H_PIR)%d 
         u(H_PHII)%d = 0.0d0
         u(H_PII)%d  = u(H_PII)%d -(omega2*u(H_PIR)%d)/(u(H_TRK)%d)  
         u(H_PIR)%d  = 0.0d0
         u(H_TRK)%d  = 0.0d0

       end if         

! and add the EM initial data...

       u(H_EX)%d = 0.d0
       u(H_EY)%d = 0.d0 
       u(H_EZ)%d = 0.d0
       u(H_BX)%d = 0.d0 
       u(H_BY)%d = 0.d0
       u(H_BZ)%d = 0.0d0
       u(H_PSI_EM)%d = 0.d0 
       u(H_PHI_EM)%d = 0.d0 

! Below is the IDs for EMD/EMDA theories
! GHS solution
    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 11) then
      !-----------------------------------------------------------
      ! single Einstein-Maxwell-Dilaton (EMD) black hole 
      ! (the GHS solution)  
      !-----------------------------------------------------------
      write(0,*)' initial data:  single Einstein-Maxwell-Dilaton (EMD) black hole ' 


      ! Input values should include alpha_0 (dil_alpha), the BH 
      ! mass (M), the BH charge (Q), and the asymptotic value of  
      ! the dilaton (phi_0).  The charge can be either magnetic or 
      ! electric.  
      phi_0     = par(P_DIL_INFTY)
      dil_alpha = par(P_DIL_ALPHA)

      if ( emd_bh_type .eq. 0 ) then ! magnetic  
         charge_magnetic = v_e1amp 
         charge_electric = 0.d0 
         charge1 = charge_magnetic * exp( - dil_alpha * phi_0 ) 
      else if ( emd_bh_type .eq. 1 ) then !electric 
         charge_magnetic = 0.d0 
         charge_electric = v_e1amp 
         charge1 = charge_electric * exp( dil_alpha * phi_0 ) 
      else 
         write(0,*)' initial data:  unknown EMD BH type ' 
         stop 
      endif 

      mass1_sq     = mass1 * mass1
      charge1_sq   = charge1 * charge1
      dil_alpha_sq = dil_alpha * dil_alpha

      discriminant = 1.d0 + (1.d0 - dil_alpha_sq)*charge1_sq/mass1_sq
      if ( discriminant .lt. 0.d0 ) then  
         write(0,*) ' initial data:  discriminant is < 0 ' 
      endif
      sqrt_disc = sqrt( discriminant )

      r_plus  = mass1 * ( 1.d0 + sqrt_disc ) 
      r_minus = charge1_sq/mass1*(1.d0 + dil_alpha_sq)/(1.d0 + sqrt_disc) 
     
      rbar_1 = 0.5d0 * ( sqrt( r_plus ) - sqrt( r_minus) ) 
      rbar_1 = rbar_1 * rbar_1
      rbar_2 = 0.5d0 * ( sqrt( r_plus ) + sqrt( r_minus) ) 
      rbar_2 = rbar_2 * rbar_2 
      rbar_H = 0.25d0 * ( r_plus - r_minus ) 
      beta_0 =  2.d0 * dil_alpha_sq / ( 1.d0 + dil_alpha_sq ) 
      
      do k=1,nz
         do j=1,ny
            do i=1,nx

               x1 = xx(i,j,k) - bh1x
               y1 = yy(i,j,k) - bh1y
               z1 = zz(i,j,k) - bh1z

               rbar = sqrt(x1*x1 + y1*y1 + z1*z1)
               if ( rbar .lt. 0.0001d0 ) then
                  rbar = 0.0001d0
               endif

               vn1(1) = x1 / rbar
               vn1(2) = y1 / rbar
               vn1(3) = z1 / rbar

               r_mi_rH = 1.d0 - ( rbar_H / rbar ) 
               r_pl_rH = 1.d0 + ( rbar_H / rbar ) 
               r_pl_r1 = 1.d0 + ( rbar_1 / rbar ) 
               r_pl_r2 = 1.d0 + ( rbar_2 / rbar ) 
              
               alpha_sq =   r_mi_rH**2 * r_pl_rH**(2.d0 - 2.d0 * beta_0 ) & 
     &                    / ( r_pl_r1 * r_pl_r2 )**(2.d0 - beta_0)  
               
               ! this can be thought of as inv_chi         
               vpsibl =   (r_pl_r1*r_pl_r2)**(2.d0-beta_0)     & 
     &                  * r_pl_rH**(2.d0*beta_0) 


               u(H_CHI)%d(i,j,k)  = 1.0 / vpsibl  
               u(H_GT11)%d(i,j,k) = 1.0d0
               u(H_GT12)%d(i,j,k) = 0.0d0
               u(H_GT13)%d(i,j,k) = 0.0d0
               u(H_GT22)%d(i,j,k) = 1.0d0
               u(H_GT23)%d(i,j,k) = 0.0d0
               u(H_GT33)%d(i,j,k) = 1.0d0

               u(H_GAM1)%d(i,j,k) = 0.0d0
               u(H_GAM2)%d(i,j,k) = 0.0d0
               u(H_GAM3)%d(i,j,k) = 0.0d0
        
               u(H_TRK)%d(i,j,k) = 0.0d0
               u(H_A11)%d(i,j,k) = 0.0d0
               u(H_A12)%d(i,j,k) = 0.0d0
               u(H_A13)%d(i,j,k)=  0.0d0
               u(H_A22)%d(i,j,k) = 0.0d0
               u(H_A33)%d(i,j,k) = 0.0d0

               u(H_ALPHA)%d(i,j,k)  = sqrt( alpha_sq )
               u(H_SHIFT1)%d(i,j,k) = 0.0d0
               u(H_SHIFT2)%d(i,j,k) = 0.0d0
               u(H_SHIFT3)%d(i,j,k) = 0.0d0

               u(H_GB1)%d(i,j,k) =  0.0d0
               u(H_GB2)%d(i,j,k) =  0.0d0
               u(H_GB3)%d(i,j,k) =  0.0d0

               tmp_dil = 0.5d0*beta_0/dil_alpha*log(r_pl_rH*r_pl_rH/r_pl_r1/r_pl_r2) 
               if ( emd_bh_type .eq. 0 ) then ! magnetic  
                  u(H_PHIR)%d(i,j,k) = phi_0 - tmp_dil  
               else if ( emd_bh_type .eq. 1 ) then ! electric  
                  u(H_PHIR)%d(i,j,k) = phi_0 + tmp_dil  
               else 
                  write(0,*)' initial data:  unknown EMD BH type ' 
                  stop 
               endif 
               u(H_PIR)%d(i,j,k)  = 0.0d0

               if ( emd_bh_type .eq. 0 ) then ! magnetic  
                  tmp_Brbar = charge_magnetic/vpsibl/sqrt(vpsibl)/rbar**2  
                  u(H_EX)%d(i,j,k) = 0.d0
                  u(H_EY)%d(i,j,k) = 0.d0
                  u(H_EZ)%d(i,j,k) = 0.d0
                  u(H_BX)%d(i,j,k) = vn1(1) * tmp_Brbar 
                  u(H_BY)%d(i,j,k) = vn1(2) * tmp_Brbar 
                  u(H_BZ)%d(i,j,k) = vn1(3) * tmp_Brbar 
                  u(H_PSI_EM)%d(i,j,k) = 0.d0
                  u(H_PHI_EM)%d(i,j,k) = 0.d0
               else if ( emd_bh_type .eq. 1 ) then ! electric  
                  tmp_Erbar = - charge_electric/(rbar*r_pl_r1*r_pl_r2)**2/sqrt(vpsibl) 
                  u(H_EX)%d(i,j,k) = vn1(1) * tmp_Erbar 
                  u(H_EY)%d(i,j,k) = vn1(2) * tmp_Erbar 
                  u(H_EZ)%d(i,j,k) = vn1(3) * tmp_Erbar 
                  u(H_BX)%d(i,j,k) = 0.d0
                  u(H_BY)%d(i,j,k) = 0.d0
                  u(H_BZ)%d(i,j,k) = 0.d0
                  u(H_PSI_EM)%d(i,j,k) = 0.d0
                  u(H_PHI_EM)%d(i,j,k) = 0.d0
               else 
                  write(0,*)' initial data:  unknown EMD BH type ' 
                  stop 
               endif 


               !emd_gaussian_pert = 1.0 + emd_gaussian_amp * exp(-(rbar-emd_gaussian_cent)**2 / emd_gaussian_width**2 )
               !u(H_EX)%d(i,j,k) = u(H_EX)%d(i,j,k) * emd_gaussian_pert
               !u(H_EY)%d(i,j,k) = u(H_EY)%d(i,j,k) * emd_gaussian_pert
               !u(H_EZ)%d(i,j,k) = u(H_EZ)%d(i,j,k) * emd_gaussian_pert
!
               !u(H_BX)%d(i,j,k) = u(H_BX)%d(i,j,k) * emd_gaussian_pert
               !u(H_BY)%d(i,j,k) = u(H_BY)%d(i,j,k) * emd_gaussian_pert
               !u(H_BZ)%d(i,j,k) = u(H_BZ)%d(i,j,k) * emd_gaussian_pert
!
               !u(H_PHIR)%d(i,j,k) = u(H_PHIR)%d(i,j,k) * emd_gaussian_pert
               !!!u(H_PIR)%/d(i,j,k)  = u(H_PIR)%d(i,j,k) + emd_gaussian_pert  


             enddo
          enddo
       enddo


! HHKK solution
    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 12) then
      !-------------------------------------------------------------
      ! Single, rotating Einstein-Maxwell-Dilaton (EMD) black hole. 
      ! This is the Horne-Horowitz-Kaluza-Klein (HHKK) solution for 
      ! alpha_0 = sqrt(3) .  
      !-------------------------------------------------------------
      write(0,*)' initial data:  single, rotating Einstein-Maxwell-Dilaton (EMD) black hole for alpha_0 = sqrt(3)'


      ! The alpha_0 (dil_alpha) value should be sqrt(3) for the exact 
      ! solution.  Free parameters in this case include the mass (M), 
      ! the BH charge (Q_e or Q_m) and the spin parameter (J).  The 
      ! charge can be either magnetic or electric.  The asymptotic value of 
      ! the dilaton is taken to be 0.  
      phi_0     = par(P_DIL_INFTY)
      dil_alpha = par(P_DIL_ALPHA)

      if ( phi_0 .ne. 0.0 ) then
         write(0,*)' initial data:  The asymptotic value of the dilaton '
         write(0,*)' initial data:  has been set to something other than '
         write(0,*)' initial data:  zero, namely phi_0 = ', phi_0
         write(0,*)' initial data:  However, the exact solution used and '
         write(0,*)' initial data:  hard coded here is zero. '
      endif

      if ( dil_alpha .ne. sqrt(3.0) ) then
         write(0,*)' initial data:  The exact solution used here assumes '
         write(0,*)' initial data:  that alpha_0 is sqrt(3) for Kaluza-Klein.'
         write(0,*)' initial data:  But alpha_0 has been set to ', dil_alpha
         write(0,*)' initial data:  You can do that, but this will not be '
         write(0,*)' initial data:  the exact solution. '
      endif

      if ( emd_bh_type .eq. 0 ) then ! magnetic  
         charge_magnetic = v_e1amp
         charge_electric = 0.d0
         charge1 = charge_magnetic * exp( - dil_alpha * phi_0 )
      else if ( emd_bh_type .eq. 1 ) then !electric 
         charge_magnetic = 0.d0
         charge_electric = v_e1amp
         charge1 = charge_electric * exp( dil_alpha * phi_0 )
      else
         write(0,*)' initial data:  unknown EMD BH type '
         stop
      endif

      mass1_sq     = mass1 * mass1
      charge1_sq   = charge1 * charge1
      dil_alpha_sq = dil_alpha * dil_alpha

      discriminant = 1.d0 + (1.d0 - dil_alpha_sq)*charge1_sq/mass1_sq
      if ( discriminant .lt. 0.d0 ) then
         write(0,*) ' initial data:  discriminant is < 0 '
      endif
      sqrt_disc = sqrt( discriminant )

      r_plus  = mass1 * ( 1.d0 + sqrt_disc )
      r_minus = charge1_sq/mass1*(1.d0 + dil_alpha_sq)/(1.d0 + sqrt_disc)

      rbar_1 = 0.5d0 * ( sqrt( r_plus ) - sqrt( r_minus) )
     rbar_1 = rbar_1 * rbar_1
      rbar_2 = 0.5d0 * ( sqrt( r_plus ) + sqrt( r_minus) )
      rbar_2 = rbar_2 * rbar_2
      rbar_H = 0.25d0 * ( r_plus - r_minus )
      beta_0 =  2.d0 * dil_alpha_sq / ( 1.d0 + dil_alpha_sq )

      do k=1,nz
         do j=1,ny
            do i=1,nx

               xbar = xx(i,j,k) - bh1x
               ybar = yy(i,j,k) - bh1y
               zbar = zz(i,j,k) - bh1z

               rhobar = sqrt(xbar*xbar + ybar*ybar)
               rbar   = sqrt(xbar*xbar + ybar*ybar + zbar*zbar)

               if ( rhobar .lt. 1.d-6 ) then
                  rhobar = 1.d-6
               endif
               if ( rbar .lt. 1.d-6 ) then
                  rbar = 1.d-6
               endif

               cos_theta =   zbar / rbar
               sin_theta = rhobar / rbar

               cos_phi = xbar / rhobar
               sin_phi = ybar / rhobar

               !vn1(1) = x1 / rbar
               !vn1(2) = y1 / rbar
               !vn1(3) = z1 / rbar

               ! we need the Boyer-Lindquist type radial coordinate as it 
               ! shows up in our various metric and EM functions 
               r_BL = ( (rbar + 0.5*mass)**2 - (0.5*a_spin)**2 ) / rbar

               ! build quantities used in the metric 
               rho_sq = r_BL**2 + (a_spin * cos_theta)**2
               Delta  = r_BL * (r_BL - 2.d0*mass) + a_spin**2
               BB     =   1.d0                                         &
     &                  + 2.d0*mass*r_BL*v_boost**2/(1.d0-v_boost**2)  &
     &                        / rho_sq
               BB     = sqrt(BB)
               Sigma  =   BB * rho_sq * ( r_BL**2 + a_spin**2 )   &
     &                  + 2.d0 * mass * r_BL * (a_spin*sin_theta)**2

               r_mi_rH = 1.d0 - ( rbar_H / rbar )
               r_pl_rH = 1.d0 + ( rbar_H / rbar )
               r_pl_r1 = 1.d0 + ( rbar_1 / rbar )
               r_pl_r2 = 1.d0 + ( rbar_2 / rbar )

               alpha_sq =   Delta * rho_sq * BB / Sigma
               ! inv_chi is what we called vpsibl in the static EMD BH case 
               inv_chi  = rbar**2 / (BB * rho_sq * Sigma)
               inv_chi  = inv_chi**(1.d0/3.d0)
               chi = 1.d0 / inv_chi
               CC = ( BB * rho_sq )**2 / Sigma
               CC_tothird     = CC**(1.0d0/3.d0)

               !u(H_CHI)%d(i,j,k)  = 1.0 / inv_chi 
               u(H_CHI)%d(i,j,k)  = chi

               u(H_GT11)%d(i,j,k) =   ( CC * cos_phi**2 + sin_phi**2 )  &
     &                              / CC_tothird**2
               u(H_GT12)%d(i,j,k) =   cos_phi * sin_phi * ( CC - 1.d0 ) &
     &                              / CC_tothird**2
               u(H_GT13)%d(i,j,k) = 0.0d0
               u(H_GT22)%d(i,j,k) =   ( CC * sin_phi**2 + cos_phi**2 )  &
     &                              / CC_tothird**2
               u(H_GT23)%d(i,j,k) = 0.0d0
               u(H_GT33)%d(i,j,k) = CC_tothird

               dSigma_dr =     ( r_BL + mass*v_boost**2/(1.d0 - v_boost**2) ) &
     &                       * ( r_BL**2 + a_spin**2 )                        &
     &                       / BB                                             &
     &                     + BB * r_BL                                        &
     &                          * ( r_BL**2 + a_spin**2 + 2.d0 * rho_sq )     &
     &                     + 2.d0 * mass * (a_spin*sin_theta**2)

               dSigma_dtheta = - 2.d0 * Delta                 &
     &                         -   ( r_BL**2 + a_spin**2 )    &
     &                           * ( 1.d0 - BB )**2           &
     &                           / BB


               dCC_dr = (   2.d0 * (   r_BL                          &
     &                               + mass * v_boost**2             &
     &                                      / ( 1.d0 - v_boost**2 )  &
     &                             )                                 &
     &                           * rho_sq                            &
     &                           * Sigma                             &
     &                    + 2.d0 * r_BL * BB**2 * rho_sq * Sigma     &
     &                    + (BB*rho_sq)**2 * dSigma_dr               &
     &                  ) / Sigma**2

               dCC_dtheta = ( - 2.d0 * a_spin**2                  &
     &                               * sin_theta                  &
     &                               * cos_theta                  &
     &                               * Sigma                      &
     &                               * rho_sq                     &
     &                               * ( 1.d0 + BB**2 )           &
     &                        + (BB*rho_sq)**2 * dSigma_dtheta    &
     &                      ) / Sigma**2

               gamtilde_tmp =   (     sin_theta                            &
     &                              * (   sin_theta * sqrt(Delta) * dCC_dr &
     &                                  + cos_theta * dCC_dtheta           &
     &                                )                                    &
     &                            - 3.d0 * CC * ( 1.d0 - CC )              &
     &                          )                                          &
     &                            / (   3.d0                               &
     &                                * rhobar                             &
     &                                * CC_tothird**4                      &
     &                              )


               u(H_GAM1)%d(i,j,k) = cos_phi * gamtilde_tmp
               u(H_GAM2)%d(i,j,k) = sin_phi * gamtilde_tmp
               u(H_GAM3)%d(i,j,k) =   (   cos_theta * sqrt(Delta) * dCC_dr &
     &                                  - sin_theta * dCC_dtheta           &
     &                                )                                    &
     &                                  / (   3.d0                         &
     &                                      * rbar                         &
     &                                      * CC_tothird**4                &
     &                                    )

               u(H_TRK)%d(i,j,k) = 0.0d0

               K_rphi =   a_spin                           &
     &                  * mass                             &
     &                  * sin_theta**2                     &
     &                  * ( r_BL * dSigma_dr - Sigma )     &
     &                  / sqrt( Delta * rho_sq * Sigma )

               K_thetaphi = - 2.d0 * a_spin**3                      &
     &                             * mass                           &
     &                             * r_BL                           &
     &                             * sin_theta**3                   &
     &                             * cos_theta                      &
     &                             / rho_sq                         &
     &                             * sqrt( Delta / rho_sq / Sigma )

               K_ij_tmp1 =   (   sin_theta * sqrt(Delta) * K_rphi  &
     &                         + cos_theta * K_thetaphi            &
     &                       )                                     &
     &                         / rbar                              &
     &                         / rhobar

               K_ij_tmp2 =   (   cos_theta * sqrt(Delta) * K_rphi  &
     &                         - sin_theta * K_thetaphi            &
     &                       )                                     &
     &                         / rbar                              &
     &                         / rhobar


               u(H_A11)%d(i,j,k) = - 2.d0 * chi*sin_phi*cos_phi*K_ij_tmp1
               u(H_A22)%d(i,j,k) = - u(H_A11)%d(i,j,k)
               u(H_A12)%d(i,j,k) =   chi*(cos_phi**2 - sin_phi**2)*K_ij_tmp1
               u(H_A13)%d(i,j,k) = - chi*sin_phi * K_ij_tmp2
               u(H_A23)%d(i,j,k) =   chi*cos_phi * K_ij_tmp2
               u(H_A33)%d(i,j,k) = 0.0d0

               u(H_ALPHA)%d(i,j,k)  = sqrt( alpha_sq )

               betau_phi = - 2.d0 * a_spin                     &
     &                            * mass                       &
     &                            * r_BL                       &
     &                            / sqrt( 1.d0 - v_boost**2 )  &
     &                            / Sigma


               u(H_SHIFT1)%d(i,j,k) = - ybar * betau_phi
               u(H_SHIFT2)%d(i,j,k) =   xbar * betau_phi
               u(H_SHIFT3)%d(i,j,k) = 0.0d0

               u(H_GB1)%d(i,j,k) =  0.0d0
               u(H_GB2)%d(i,j,k) =  0.0d0
               u(H_GB3)%d(i,j,k) =  0.0d0

               tmp_dil = 0.5d0 * sqrt(3.d0) * log(BB)

               if ( emd_bh_type .eq. 0 ) then ! magnetic  
                  u(H_PHIR)%d(i,j,k) =   tmp_dil
               else if ( emd_bh_type .eq. 1 ) then ! electric  
                  u(H_PHIR)%d(i,j,k) = - tmp_dil
               else
                  write(0,*)' initial data:  unknown EMD BH type '
                  stop
               endif
               u(H_PIR)%d(i,j,k)  = 0.0d0


               if ( emd_bh_type .eq. 0 ) then ! magnetic 
                  F_tr       = - 2.d0 * a_spin                     &
     &                                * mass                       &
     &                                * r_BL                       &
     &                                * v_boost                    &
     &                                / sqrt( 1.d0 - v_boost**2 )  &
     &                                * cos_theta                  &
     &                                / rho_sq**2
                  F_ttheta =   a_spin                              &
     &                       * mass                                &
     &                       * v_boost                             &
     &                       / sqrt( 1.d0 - v_boost**2 )           &
     &                       * ( 2.d0 * r_BL**2 - rho_sq )         &
     &                       / rho_sq**2                           &
     &                       * sin_theta

                  F_rphi = - 2.d0 * mass                           &
     &                            * r_BL                           &
     &                            * v_boost                        &
     &                            / ( 1.d0 - v_boost**2 )          &
     &                            / rho_sq**2                      &
     &                            * a_spin**2                      &
     &                            * sin_theta**2                   &
     &                            * cos_theta

                  F_thetaphi = -   mass                            &
     &                           * v_boost                         &
     &                           / ( 1.d0 - v_boost**2 )           &
     &                           * ( 2.d0 * r_BL**2 - rho_sq )     &
     &                           / rho_sq**2                       &
     &                           * ( r_BL**2 + a_spin**2 )         &
     &                           * sin_theta
               else if ( emd_bh_type .eq. 1 ) then ! electric  
                  F_tr       =   mass                              &
     &                         * v_boost                           &
     &                         / ( 1.d0 - v_boost**2 )             &
     &                         * ( 2.d0 * r_BL**2 - rho_sq )       &
     &                         / ( BB**2 * rho_sq )**2
                  F_ttheta   = - 2.d0 * mass                       &
     &                                * a_spin**2                  &
     &                                * r_BL                       &
     &                                * v_boost                    &
     &                                / ( 1.d0 - v_boost**2 )      &
     &                                / ( BB**2 * rho_sq )**2      &
     &                                * sin_theta**2               &
     &                                * cos_theta
                  F_rphi     =   a_spin                            &
     &                         * mass                              &
     &                         * v_boost                           &
     &                         / sqrt( 1.d0 - v_boost**2 )         &
     &                         * ( 2.d0 * r_BL**2 - rho_sq )       &
     &                         / ( BB**2 * rho_sq )**2             &
     &                         * sin_theta**2

                  F_thetaphi = - 2.d0 * a_spin                     &
     &                                * mass                       &
     &                                * r_BL                       &
     &                                * v_boost                    &
     &                                / ( BB**2 * rho_sq )**2      &
     &                                / sqrt( 1.d0 - v_boost**2 )  &
     &                                * (   BB**2 * rho_sq         &
     &                                    + ( a_spin * sin_theta )**2  &
     &                                  )                          &
     &                                * sin_theta                  &
     &                                * cos_theta
               else
                  write(0,*)' initial data:  unknown EMD BH type '
                  stop
               endif


               E_up_tmp1 = -   chi                                          &
     &                       / alpha                                        &
     &                       / CC_tothird                                   &
     &                       / rbar                                         &
     &                       * (     sin_theta                              &
     &                             * sqrt(Delta)                            &
     &                             * ( F_tr + betau_phi * F_rphi )          &
     &                           +   cos_theta                              &
     &                             * ( F_ttheta + betau_phi * F_thetaphi )  &
     &                         )

               E_up_tmp2 = -   chi                                          &
     &                       / alpha                                        &
     &                       / CC_tothird                                   &
     &                       / rbar                                         &
     &                       * (     cos_theta                              &
     &                             * sqrt(Delta)                            &
     &                             * ( F_tr + betau_phi * F_rphi )          &
     &                           -   sin_theta                              &
     &                             * ( F_ttheta + betau_phi * F_thetaphi )  &
     &                         )

               u(H_EX)%d(i,j,k) = cos_phi * E_up_tmp1
               u(H_EY)%d(i,j,k) = sin_phi * E_up_tmp1
               u(H_EZ)%d(i,j,k) =           E_up_tmp2

               B_up_tmp1 =   chi**(1.5d0)                           &
     &                     / rbar                                   &
     &                     / rhobar                                 &
     &                     * (   sin_theta * sqrt(Delta) * F_rphi   &
     &                         + cos_theta * F_thetaphi             &
     &                       )

               B_up_tmp2 = -   chi**(1.5d0)                           &
     &                       / rbar                                   &
     &                       / rhobar                                 &
     &                       * (   cos_theta * sqrt(Delta) * F_rphi   &
     &                           - sin_theta * F_thetaphi             &
     &                         )

               u(H_BX)%d(i,j,k) = cos_phi * B_up_tmp2
               u(H_BY)%d(i,j,k) = sin_phi * B_up_tmp2
               u(H_BZ)%d(i,j,k) =           B_up_tmp1


               u(H_PSI_EM)%d(i,j,k) = 0.d0
               u(H_PHI_EM)%d(i,j,k) = 0.d0

               !emd_gaussian_pert = 1.0 + emd_gaussian_amp * exp(-(rbar-emd_gaussian_cent)**2 / emd_gaussian_width**2 )  
               !u(H_EX)%d(i,j,k) = u(H_EX)%d(i,j,k) * emd_gaussian_pert    
               !u(H_EY)%d(i,j,k) = u(H_EY)%d(i,j,k) * emd_gaussian_pert  
               !u(H_EZ)%d(i,j,k) = u(H_EZ)%d(i,j,k) * emd_gaussian_pert  
              ! 
               !u(H_BX)%d(i,j,k) = u(H_BX)%d(i,j,k) * emd_gaussian_pert  
               !u(H_BY)%d(i,j,k) = u(H_BY)%d(i,j,k) * emd_gaussian_pert  
               !u(H_BZ)%d(i,j,k) = u(H_BZ)%d(i,j,k) * emd_gaussian_pert  
              ! 
               !u(H_PHIR)%d(i,j,k) = u(H_PHIR)%d(i,j,k) * emd_gaussian_pert  
               !u(H_PIR)%d(i,j,k)  = u(H_PIR)%d(i,j,k) + emd_gaussian_pert  


             enddo
          enddo
       enddo

! KS solution
    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 13) then
      !-------------------------------------------------------------
      ! Single, Einstein-Maxwell-Dilaton-Axion (EMDA) black hole. 
      ! This is comes from Kerr-Sen BH
      !-------------------------------------------------------------
      write(0,*)' initial data: Single, EMDA black hole from Kerr-Sen'

      phi_0     = par(P_DIL_INFTY)
      dil_alpha = par(P_DIL_ALPHA)
      axn_alpha = par(P_AXN_ALPHA)
      sen_alpha = par(P_SEN_ALPHA) !Alpha value for Kerr-Sen data

      if ( phi_0 .ne. 0.0 ) then
         write(0,*)' initial data:  The asymptotic value of the dilaton '
         write(0,*)' initial data:  has been set to something other than '
         write(0,*)' initial data:  zero, namely phi_0 = ', phi_0
         write(0,*)' initial data:  However, the exact solution used and '
         write(0,*)' initial data:  hard coded here is zero. '
      endif

      if ( dil_alpha .ne. sqrt(3.0) ) then
         write(0,*)' initial data:  The exact solution used here assumes '
         write(0,*)' initial data:  that alpha_0 is sqrt(3) for Kaluza-Klein.'
         write(0,*)' initial data:  But alpha_0 has been set to ', dil_alpha
         write(0,*)' initial data:  You can do that, but this will not be '
         write(0,*)' initial data:  the exact solution. '
      endif

      if ( sen_alpha == 0 ) then
        write(0,*)' Your sen_alpha = ', sen_alpha
        write(0,*)' Please use non-zero sen_alpha.'
      endif

      discriminant = 1.d0 + (1.d0 - dil_alpha_sq)*charge1_sq/mass1_sq
      if ( discriminant .lt. 0.d0 ) then
         write(0,*) ' initial data:  discriminant is < 0 '
      endif
      sqrt_disc = sqrt( discriminant )

      ! Some variables for quasi-isotropic system for Kerr-Sen metric
      rbar_1 = mass1 + a_ang_par
      rbar_2 = mass1 - a_ang_par
      rbar_plus = 0.5d0 * (r_BL - mass1 + sqrt( (r_BL - mass1) * (r_BL - mass1) &
      &           - rbar_1 * rbar_2)) 
      rbar_minus = 0.5d0 * (r_BL - mass1 - sqrt( (r_BL - mass1) * (r_BL - mass1) &
      &            -  rbar_1 * rbar_2)) 
    
    !  EMD stuffs, keep for while (HL)    
    !  r_plus  = mass1 * ( 1.d0 + sqrt_disc )
    !  r_minus = charge1_sq/mass1*(1.d0 + dil_alpha_sq)/(1.d0 + sqrt_disc)

    !  rbar_H = 0.25d0 * ( r_plus - r_minus )
    !

    !  beta_0 =  2.d0 * dil_alpha_sq / ( 1.d0 + dil_alpha_sq )

      do k=1,nz
         do j=1,ny
            do i=1,nx

               xbar = xx(i,j,k) - bh1x
               ybar = yy(i,j,k) - bh1y
               zbar = zz(i,j,k) - bh1z

               rhobar = sqrt(xbar*xbar + ybar*ybar)
               rbar   = sqrt(xbar*xbar + ybar*ybar + zbar*zbar)

               if ( rhobar .lt. 1.d-6 ) then
                  rhobar = 1.d-6
               endif
               if ( rbar .lt. 1.d-6 ) then
                  rbar = 1.d-6
               endif

               cos_theta =   zbar / rbar
               sin_theta = rhobar / rbar

               cos_phi = xbar / rhobar
               sin_phi = ybar / rhobar

               !vn1(1) = x1 / rbar
               !vn1(2) = y1 / rbar
               !vn1(3) = z1 / rbar

               ! we need the Boyer-Lindquist type radial coordinate as it 
               ! shows up in our various metric and EM functions 
               ! This comes from Kerr-Sen

               r_BL = ( (rbar + rbar_1 ) * (rbar + rbar_2) ) / rbar

               ! build quantities used in the metric 
               rho_sq = r_BL**2 + (a_ang_par * cos_theta)**2 + 2 * mass * r_BL &
     &                  * (sinh(sen_alpha))**2
               Delta  = r_BL * (r_BL - 2.d0*mass1) + a_ang_par**2
               Sigma  = (rho_sq + (a_ang_par * sin_theta)**2)**2 - Delta       &
     &                  * (a_ang_par * sin_theta)**2  

              ! --------------------------------------------------------------
              ! This is EMD stuffs, keep these for while (HL)
              ! r_mi_rH = 1.d0 - ( rbar_H / rbar )
              ! r_pl_rH = 1.d0 + ( rbar_H / rbar )
              ! r_pl_r1 = 1.d0 + ( rbar_1 / rbar )
              ! r_pl_r2 = 1.d0 + ( rbar_2 / rbar )
  
               alpha_sq =   Delta * rho_sq / Sigma
              ! inv_chi is what we called vpsibl in the static EMD BH case 
               inv_chi  = rho_sq*Sigma/rbar**6
               inv_chi  = inv_chi**(1.d0/3.d0)
               chi = 1.d0 / inv_chi
              ! CC = ( BB * rho_sq )**2 / Sigma
              ! CC_tothird     = CC**(1.0d0/3.d0)

               !u(H_CHI)%d(i,j,k)  = 1.0 / inv_chi 
               !u(H_CHI)%d(i,j,k)  = chi
              ! ---------------------------------------------------------------

              CC = rho_sq**2 / Sigma
              CC_tothird = CC**(1.0d0/3.0d0)


               u(H_GT11)%d(i,j,k) =   ( CC * cos_phi**2 + sin_phi**2 )  &
     &                              / CC_tothird**2
               u(H_GT12)%d(i,j,k) =   cos_phi * sin_phi * ( CC - 1.d0 ) &
     &                              / CC_tothird**2
               u(H_GT13)%d(i,j,k) = 0.0d0
               u(H_GT22)%d(i,j,k) =   ( CC * sin_phi**2 + cos_phi**2 )  &
     &                              / CC_tothird**2
               u(H_GT23)%d(i,j,k) = 0.0d0
               u(H_GT33)%d(i,j,k) = CC_tothird

               dSigma_dr = 4.0d0 * rho_sq * (r_BL + mass1               &  
     &                     * (sinh(sen_alpha)**2))                      &
     &                     + 2.0d0*((a_ang_par*sin_theta)**2)           &
     &                     *(r_BL + mass1 * cosh(2.0d0*sen_alpha))
     
    
              dSigma_dtheta = - 2.d0 * Delta * a_ang_par**2             &
     &                        * sin_theta * cos_theta     


               dCC_dr = (   2.d0 * (   r_BL                          &
     &                               + mass * v_boost**2             &
     &                                      / ( 1.d0 - v_boost**2 )  &
     &                             )                                 &
     &                           * rho_sq                            &
     &                           * Sigma                             &
     &                    + 2.d0 * r_BL * BB**2 * rho_sq * Sigma     &
     &                    + (BB*rho_sq)**2 * dSigma_dr               &
     &                  ) / Sigma**2

               dCC_dtheta = ( - 2.d0 * a_spin**2                  &
     &                               * sin_theta                  &
     &                               * cos_theta                  &
     &                               * Sigma                      &
     &                               * rho_sq                     &
     &                               * ( 1.d0 + BB**2 )           &
     &                        + (BB*rho_sq)**2 * dSigma_dtheta    &
     &                      ) / Sigma**2

               gamtilde_tmp =   (     sin_theta                            &
     &                              * (   sin_theta * sqrt(Delta) * dCC_dr &
     &                                  + cos_theta * dCC_dtheta           &
     &                                )                                    &
     &                            - 3.d0 * CC * ( 1.d0 - CC )              &
     &                          )                                          &
     &                            / (   3.d0                               &
     &                                * rhobar                             &
     &                                * CC_tothird**4                      &
     &                              )


               u(H_GAM1)%d(i,j,k) = cos_phi * gamtilde_tmp
               u(H_GAM2)%d(i,j,k) = sin_phi * gamtilde_tmp
               u(H_GAM3)%d(i,j,k) =   (   cos_theta * sqrt(Delta) * dCC_dr &
     &                                  - sin_theta * dCC_dtheta           &
     &                                )                                    &
     &                                  / (   3.d0                         &
     &                                      * rbar                         &
     &                                      * CC_tothird**4                &
     &                                    )

               u(H_TRK)%d(i,j,k) = 0.0d0

               K_rphi = a_ang_par                                   &
     &                  * mass1                                     &
     &                  * sin_theta**2                              &
     &                  * sinh(sen_alpha)**2                        &
     &                  * ( r_BL * dSigma_dr - Sigma )              &
     &                  / (rho_sq * sqrt( Delta * rho_sq * Sigma ))

               K_thetaphi = a_ang_par                               &
     &                      * mass1                                 &
     &                      * r_BL                                  &
     &                      * sin_theta**2                          &
     &                      * cosh(sen_alpha)**2                    &
     &                      * dSigma_dtheta                         &
     &                      /(rho_sq * sqrt( Delta * rho_sq * Sigma ))

               K_ij_tmp1 =   (   sin_theta * sqrt(Delta) * K_rphi  &
     &                         + cos_theta * K_thetaphi            &
     &                       )                                     &
     &                         / rbar                              &
     &                         / rhobar

               K_ij_tmp2 =   (   cos_theta * sqrt(Delta) * K_rphi  &
     &                         - sin_theta * K_thetaphi            &
     &                       )                                     &
     &                         / rbar                              &
     &                         / rhobar


               u(H_A11)%d(i,j,k) = - 2.d0 * chi*sin_phi*cos_phi*K_ij_tmp1
               u(H_A22)%d(i,j,k) = - u(H_A11)%d(i,j,k)
               u(H_A12)%d(i,j,k) =   chi*(cos_phi**2 - sin_phi**2)*K_ij_tmp1
               u(H_A13)%d(i,j,k) = - chi*sin_phi * K_ij_tmp2
               u(H_A23)%d(i,j,k) =   chi*cos_phi * K_ij_tmp2
               u(H_A33)%d(i,j,k) = 0.0d0

               u(H_ALPHA)%d(i,j,k)  = sqrt( alpha_sq )

               betau_phi = - 2.d0 * a_ang_par                  &
     &                            * mass1                      &
     &                            * r_BL                       &
     &                            * cosh(sen_alpha)**2         &
     &                            / Sigma


               u(H_SHIFT1)%d(i,j,k) = - ybar * betau_phi
               u(H_SHIFT2)%d(i,j,k) =   xbar * betau_phi
               u(H_SHIFT3)%d(i,j,k) = 0.0d0

               u(H_GB1)%d(i,j,k) =  0.0d0
               u(H_GB2)%d(i,j,k) =  0.0d0
               u(H_GB3)%d(i,j,k) =  0.0d0

               tmp_dil = 0.5d0 * sqrt(3.d0) * log(BB)

               u(H_PIR)%d(i,j,k)  = 0.0d0


                  F_tr     = mass1                                   &
     &                       * sinh(2.0d0 * sen_alpha)               &
     &                       * (r_BL**2 - (a_ang_par * cos_theta)**2)&
     &                       / (sqrt(2.0d0) * rho_sq**2)             

                  F_ttheta = - sqrt(2.0d0) * mass1                 &
     &                       * sinh(2.0d0 * sen_alpha)             &
     &                       * r_BL                                &
     &                       * a_ang_par**2                        &
     &                       * sin_theta * cos_theta               &
     &                       / rho_sq**2                           

                  F_rphi = a_ang_par                               &
     &                     * mass1                                 &
     &                     * sinh(2.0d0 *sen_alpha)                &
     &                     * (r_BL**2 - (a_ang_par * cos_theta)**2)&
     &                     * sin_theta**2                          &
                           / (sqrt(2.0d0) * rho_sq**2)

                  F_thetaphi = - sqrt(2.0d0) * mass1                    &
     &                         * a_ang_par                              &
     &                         * mass1                                  &
     &                         * sinh(2.0d0 * sen_alpha)                &
     &                         * r_BL                                   &
     &                         * (r_BL**2 + a_ang_par**2 + 2.0d0        &
     &                           * mass1 * r_BL * (sinh(sen_alpha)**2)) &
     &                         * sin_theta * cos_theta                  &
     &                         /rho_sq**2    

               E_up_tmp1 = -   chi                                          &
     &                       / alpha                                        &
     &                       / CC_tothird                                   &
     &                       / rbar                                         &
     &                       * (     sin_theta                              &
     &                             * sqrt(Delta)                            &
     &                             * ( F_tr + betau_phi * F_rphi )          &
     &                           +   cos_theta                              &
     &                             * ( F_ttheta + betau_phi * F_thetaphi )  &
     &                         )

               E_up_tmp2 = -   chi                                          &
     &                       / alpha                                        &
     &                       / CC_tothird                                   &
     &                       / rbar                                         &
     &                       * (     cos_theta                              &
     &                             * sqrt(Delta)                            &
     &                             * ( F_tr + betau_phi * F_rphi )          &
     &                           -   sin_theta                              &
     &                             * ( F_ttheta + betau_phi * F_thetaphi )  &
     &                         )

               u(H_EX)%d(i,j,k) = cos_phi * E_up_tmp1
               u(H_EY)%d(i,j,k) = sin_phi * E_up_tmp1
               u(H_EZ)%d(i,j,k) =           E_up_tmp2

               B_up_tmp1 =   chi**(1.5d0)                           &
     &                     / rbar                                   &
     &                     / rhobar                                 &
     &                     * (   sin_theta * sqrt(Delta) * F_rphi   &
     &                         + cos_theta * F_thetaphi             &
     &                       )

               B_up_tmp2 = -   chi**(1.5d0)                           &
     &                       / rbar                                   &
     &                       / rhobar                                 &
     &                       * (   cos_theta * sqrt(Delta) * F_rphi   &
     &                           - sin_theta * F_thetaphi             &
     &                         )

               u(H_BX)%d(i,j,k) = cos_phi * B_up_tmp2
               u(H_BY)%d(i,j,k) = sin_phi * B_up_tmp2
               u(H_BZ)%d(i,j,k) =           B_up_tmp1


               u(H_PSI_EM)%d(i,j,k) = 0.d0
               u(H_PHI_EM)%d(i,j,k) = 0.d0

             enddo
          enddo
       enddo !End EMDA, Kerr-Sen

! RL solution
    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 14) then
      !-------------------------------------------------------------
      ! Single, Rasheed-Larsen BH for rotating EMD solution
      !-------------------------------------------------------------
      write(0,*)' initial data: Single, Rasheed-Larsen BH'

      phi_0     = par(P_DIL_INFTY)
      dil_alpha = par(P_DIL_ALPHA)
      axn_alpha = par(P_AXN_ALPHA)

      if ( phi_0 .ne. 0.0 ) then
         write(0,*)' initial data:  The asymptotic value of the dilaton '
         write(0,*)' initial data:  has been set to something other than '
         write(0,*)' initial data:  zero, namely phi_0 = ', phi_0
         write(0,*)' initial data:  However, the exact solution used and '
         write(0,*)' initial data:  hard coded here is zero. '
      endif

      if ( dil_alpha .ne. sqrt(3.0) ) then
         write(0,*)' initial data:  The exact solution used here assumes '
         write(0,*)' initial data:  that alpha_0 is sqrt(3) for Kaluza-Klein.'
         write(0,*)' initial data:  But alpha_0 has been set to ', dil_alpha
         write(0,*)' initial data:  You can do that, but this will not be '
         write(0,*)' initial data:  the exact solution. '
      endif

      if ( sen_alpha == 0 ) then
        write(0,*)' Your sen_alpha = ', sen_alpha
        write(0,*)' Please use non-zero sen_alpha.'
      endif

      discriminant = 1.d0 + (1.d0 - dil_alpha_sq)*charge1_sq/mass1_sq
      if ( discriminant .lt. 0.d0 ) then
         write(0,*) ' initial data:  discriminant is < 0 '
      endif
      sqrt_disc = sqrt( discriminant )

      ! Some variables for quasi-isotropic system for Rasheed-Larsen metric
      rbar_plus = 0.5d0 * (r_BL - mass1 + sqrt( (r_BL**2 + a_ang_par**2   &
    &           -2.d0 * mass1* r_BL ) ) )
      rbar_minus = 0.5d0 * (r_BL - mass1 - sqrt( (r_BL**2 + a_ang_par**2   &
    &           -2.d0 * mass1* r_BL ) ) )
    !  EMD stuffs, keep for while (HL)    
    !  r_plus  = mass1 * ( 1.d0 + sqrt_disc )
    !  r_minus = charge1_sq/mass1*(1.d0 + dil_alpha_sq)/(1.d0 + sqrt_disc)

    !  rbar_H = 0.25d0 * ( r_plus - r_minus )
    !


      beta_0 =  2.d0 * dil_alpha_sq / ( 1.d0 + dil_alpha_sq )

      do k=1,nz
         do j=1,ny
            do i=1,nx

               xbar = xx(i,j,k) - bh1x
               ybar = yy(i,j,k) - bh1y
               zbar = zz(i,j,k) - bh1z

               rhobar = sqrt(xbar*xbar + ybar*ybar)
               rbar   = sqrt(xbar*xbar + ybar*ybar + zbar*zbar)

               if ( rhobar .lt. 1.d-6 ) then
                  rhobar = 1.d-6
               endif
               if ( rbar .lt. 1.d-6 ) then
                  rbar = 1.d-6
               endif

               cos_theta =   zbar / rbar
               sin_theta = rhobar / rbar

               cos_phi = xbar / rhobar
               sin_phi = ybar / rhobar

               !vn1(1) = x1 / rbar
               !vn1(2) = y1 / rbar
               !vn1(3) = z1 / rbar

               ! we need the Boyer-Lindquist type radial coordinate as it 
               ! shows up in our various metric and EM functions 

               r_BL = ( (rbar + mass1/2.d0 )**2 - a_ang_par**2/4.d0 ) / rbar

               ! Define electric and magnetic charge from constants p and q
               ! (or Define P and Q from p and q naviely..)
        
               Pt_mag = (p_mag * (p_mag**2 - 4.d0 * mass1**2))                  &
     &                 /(4.d0 * (p_mag + q_elec))

               Qt_elec = (q_elec * (q_elec**2 - 4.d0 * mass1**2))               &
     &                  /(4.d0 * (p_mag + q_elec))

               ! build quantities used in the metric
               H_3_RL = r_BL * (r_BL - 2.d0*mass1) + (a_ang_par*cos_theta)**2
               H_comm1 = ( (p_mag - 2.d0 * mass1) * (q_elec - 2.d0 * mass1) )  &
     &                  / ( 2.d0 * (p_mag + q_elec) )
               H_comm2 = sqrt((q_elec**2 -4.d0 * mass1**2) * (p_mag**2 - 4.d0 &
     &                   * mass1**2)) * (a_ang_par * cos_theta)**2             &
     &                   / (2.d0 * mass1 * (p_mag + q_elec) )
               H_1_RL = H_3_RL + r_BL*p_mag + p_mag*H_comm1 - p_mag*H_comm2
               H_2_RL = H_3_RL + r_BL*q_elec + q_elec*H_comm1 - q_elec*H_comm2
               B_RL = sqrt(p_mag*q_elec) * ( (p_mag * q_elec + 4.d0 * mass1**2) * &
     &                mass1 * (p_mag - 2 * mass1) * (q_elec - 2.d0 * mass1) *     &
     &                (a_ang_par * sin_theta)**2 ) / (2 * mass1 *              &
     &                (p_mag + q_elec) * H_3_RL)
               rho_sq = r_BL**2 + (a_ang_par * cos_theta)**2 
               Delta  = r_BL * (r_BL - 2.d0*mass1) + a_ang_par**2
               c_0_RL = (p_mag * q_elec + 4.d0 * mass1)                        & 
     &                  /(2.d0 * mass1 * (p_mag + q_elec) - 1.d0)
               c_1_RL = sqrt(c_0_RL * (c_0_RL + 2.d0)) 
               Sigma_pq  = (r_BL**2 + a_ang_par**2)                             &
     &                     * (rho_sq + r_BL * (q_elec-2.d0*mass1)               &
     &                     + r_BL * (p_mag-2.d0*mass1)                          &
     &                     + 0.5d0 * (p_mag-2.d0*mass1) * (q_elec-2.d0*mass1))  &
     &                     + 2.d0 * mass1 * r_BL * (a_ang_par*sin_theta)**2     &
     &                     + Delta*c_1_RL*(q_elec-p_mag)*a_ang_par*cos_theta**2 &
     &                     + (p_mag-2.d0*mass1)*(q_elec-2.d0*mass1)             &
     &                     * (r_BL**2 + r_BL/(2.d0*(p_mag+q_elec)               &
     &                     * (q_elec * (p_mag-2.d0*mass1)                       &
     &                     + p_mag * (q_elec-2.d0*mass1))))                     &
     &                     + p_mag * q_elec                                     &
     &                     * ((c_0_RL*mass1)**2 - (c_1_RL*a_ang_par)**2)
              ! --------------------------------------------------------------
              ! This is EMD stuffs, keep these for while (HL)
              ! r_mi_rH = 1.d0 - ( rbar_H / rbar )
              ! r_pl_rH = 1.d0 + ( rbar_H / rbar )
              ! r_pl_r1 = 1.d0 + ( rbar_1 / rbar )
              ! r_pl_r2 = 1.d0 + ( rbar_2 / rbar )
  
               alpha_sq =   Delta * rho_sq * BB / Sigma
              ! inv_chi is what we called vpsibl in the static EMD BH case 
               inv_chi  = rbar**2 / (BB * rho_sq * Sigma)
               inv_chi  = inv_chi**(1.d0/3.d0)
               chi = 1.d0 / inv_chi
              ! CC = ( BB * rho_sq )**2 / Sigma
              ! CC_tothird     = CC**(1.0d0/3.d0)

               !u(H_CHI)%d(i,j,k)  = 1.0 / inv_chi 
               !u(H_CHI)%d(i,j,k)  = chi
              ! ---------------------------------------------------------------

              CC = (H_1_RL * H_2_RL) / Sigma_pq
              CC_tothird = CC**(1.0d0/3.0d0)

              dH_1_RL_dr = 2.d0 * r_BL + p_mag - 2.d0 * mass1
              dH_2_RL_dr = 2.d0 * r_BL + q_elec - 2.d0 * mass1
              dH_1_RL_dtheta = -2.d0 * sin_theta * cos_theta * a_ang_par**2   &
     &                         + p_mag * c_1_RL * a_ang_par * sin_theta 
              dH_2_RL_dtheta = -2.d0 * sin_theta * cos_theta * a_ang_par**2   &
     &                         + q_elec * c_1_RL * a_ang_par * sin_theta 

              dSigma_pq_dr = 2.d0 * r_BL * (rho_sq + r_BL * (q_elec-2.d0*mass1)&
     &                       + r_BL * (p_mag-2.d0*mass1)                       &
     &                       + 0.5d0 * (p_mag-2.d0*mass1)*(q_elec-2.d0*mass1)) &
     &                       + (r_BL**2 + a_ang_par**2)                        &
     &                       * (2.d0*r_BL+q_elec-2.d0*mass1+p_mag-2.d0*mass1)  &
     &                       + 2.d0*mass1*(a_ang_par*sin_theta)**2             &
     &                       + 2.d0 * (r_BL-mass1) * c_1_RL * (q_elec-p_mag)   &
     &                       * cos_theta                                       &
     &                       + (p_mag-2.d0*mass1) * (q_elec-2.d0*mass1)        &
     &                       * (2.d0*r_BL + (q_elec * (p_mag - 2.d0*mass1)     &
     &                       + p_mag * (q_elec - 2.d0*mass1))                  &
     &                       /(2.d0*(p_mag+q_elec)))   
              dSigma_pq_dtheta = -2.d0 * (r_BL**2 + a_ang_par**2)              &
     &                           * a_ang_par**2 * sin_theta * cos_theta        &
     &                           + 4.d0 * mass1 * r_BL * a_ang_par**2          &
     &                           * sin_theta * cos_theta                       &
     &                           - c_1_RL*(q_elec - p_mag)*Delta*sin_theta

               u(H_GT11)%d(i,j,k) =   ( CC * cos_phi**2 + sin_phi**2 )  &
     &                              / CC_tothird**2
               u(H_GT12)%d(i,j,k) =   cos_phi * sin_phi * ( CC - 1.d0 ) &
     &                              / CC_tothird**2
               u(H_GT13)%d(i,j,k) = 0.0d0
               u(H_GT22)%d(i,j,k) =   ( CC * sin_phi**2 + cos_phi**2 )  &
     &                              / CC_tothird**2
               u(H_GT23)%d(i,j,k) = 0.0d0
               u(H_GT33)%d(i,j,k) = CC_tothird

   ! --------------------------------------------------------------------
   ! Kerr-Sen Sigma, differs with Rasheed-Larsen
   !            dSigma_dr = 4.0d0 * rho_sq * (r_BL + mass1               &  
   !  &                     * (sinh(sen_alpha)**2))                      &
   !  &                     + 2.0d0*((a_ang_par*sin_theta)**2)           &
   !  &                     *(r_BL + mass1 * cosh(2.0d0*sen_alpha))
     
    
   !           dSigma_dtheta = - 2.d0 * Delta * a_ang_par**2             &
   !  &                        * sin_theta * cos_theta     

   ! ---------------------------------------------------------------------   
               dCC_dr = (   2.d0 * (   r_BL                          &
     &                               + mass * v_boost**2             &
     &                                      / ( 1.d0 - v_boost**2 )  &
     &                             )                                 &
     &                           * rho_sq                            &
     &                           * Sigma                             &
     &                    + 2.d0 * r_BL * BB**2 * rho_sq * Sigma     &
     &                    + (BB*rho_sq)**2 * dSigma_dr               &
     &                  ) / Sigma**2

               dCC_dtheta = ( - 2.d0 * a_spin**2                  &
     &                               * sin_theta                  &
     &                               * cos_theta                  &
     &                               * Sigma                      &
     &                               * rho_sq                     &
     &                               * ( 1.d0 + BB**2 )           &
     &                        + (BB*rho_sq)**2 * dSigma_dtheta    &
     &                      ) / Sigma**2

               gamtilde_tmp =   (     sin_theta                            &
     &                              * (   sin_theta * sqrt(Delta) * dCC_dr &
     &                                  + cos_theta * dCC_dtheta           &
     &                                )                                    &
     &                            - 3.d0 * CC * ( 1.d0 - CC )              &
     &                          )                                          &
     &                            / (   3.d0                               &
     &                                * rhobar                             &
     &                                * CC_tothird**4                      &
     &                              )


               u(H_GAM1)%d(i,j,k) = cos_phi * gamtilde_tmp
               u(H_GAM2)%d(i,j,k) = sin_phi * gamtilde_tmp
               u(H_GAM3)%d(i,j,k) =   (   cos_theta * sqrt(Delta) * dCC_dr &
     &                                  - sin_theta * dCC_dtheta           &
     &                                )                                    &
     &                                  / (   3.d0                         &
     &                                      * rbar                         &
     &                                      * CC_tothird**4                &
     &                                    )

               u(H_TRK)%d(i,j,k) = 0.0d0

               K_rphi = a_ang_par                                   &
     &                  * mass1                                     &
     &                  * sin_theta**2                              &
     &                  * sinh(sen_alpha)**2                        &
     &                  * ( r_BL * dSigma_dr - Sigma )              &
     &                  / (rho_sq * sqrt( Delta * rho_sq * Sigma ))

               K_thetaphi = a_ang_par                               &
     &                      * mass1                                 &
     &                      * r_BL                                  &
     &                      * sin_theta**2                          &
     &                      * cosh(sen_alpha)**2                    &
     &                      * dSigma_dtheta                         &
     &                      /(rho_sq * sqrt( Delta * rho_sq * Sigma ))

               K_ij_tmp1 =   (   sin_theta * sqrt(Delta) * K_rphi  &
     &                         + cos_theta * K_thetaphi            &
     &                       )                                     &
     &                         / rbar                              &
     &                         / rhobar

               K_ij_tmp2 =   (   cos_theta * sqrt(Delta) * K_rphi  &
     &                         - sin_theta * K_thetaphi            &
     &                       )                                     &
     &                         / rbar                              &
     &                         / rhobar


               u(H_A11)%d(i,j,k) = - 2.d0 * chi*sin_phi*cos_phi*K_ij_tmp1
               u(H_A22)%d(i,j,k) = - u(H_A11)%d(i,j,k)
               u(H_A12)%d(i,j,k) =   chi*(cos_phi**2 - sin_phi**2)*K_ij_tmp1
               u(H_A13)%d(i,j,k) = - chi*sin_phi * K_ij_tmp2
               u(H_A23)%d(i,j,k) =   chi*cos_phi * K_ij_tmp2
               u(H_A33)%d(i,j,k) = 0.0d0

               u(H_ALPHA)%d(i,j,k)  = sqrt( alpha_sq )

               betau_phi = - 2.d0 * a_ang_par                  &
     &                            * mass1                      &
     &                            * r_BL                       &
     &                            * cosh(sen_alpha)**2         &
     &                            / Sigma


               u(H_SHIFT1)%d(i,j,k) = - ybar * betau_phi
               u(H_SHIFT2)%d(i,j,k) =   xbar * betau_phi
               u(H_SHIFT3)%d(i,j,k) = 0.0d0

               u(H_GB1)%d(i,j,k) =  0.0d0
               u(H_GB2)%d(i,j,k) =  0.0d0
               u(H_GB3)%d(i,j,k) =  0.0d0

               tmp_dil = 0.5d0 * sqrt(3.d0) * log(BB)

               u(H_PIR)%d(i,j,k)  = 0.0d0


                  F_tr     = mass1                                   &
     &                       * sinh(2.0d0 * sen_alpha)               &
     &                       * (r_BL**2 - (a_ang_par * cos_theta)**2)&
     &                       / (sqrt(2.0d0) * rho_sq**2)             

                  F_ttheta = - sqrt(2.0d0) * mass1                 &
     &                       * sinh(2.0d0 * sen_alpha)             &
     &                       * r_BL                                &
     &                       * a_ang_par**2                        &
     &                       * sin_theta * cos_theta               &
     &                       / rho_sq**2                           

                  F_rphi = a_ang_par                               &
     &                     * mass1                                 &
     &                     * sinh(2.0d0 *sen_alpha)                &
     &                     * (r_BL**2 - (a_ang_par * cos_theta)**2)&
     &                     * sin_theta**2                          &
                           / (sqrt(2.0d0) * rho_sq**2)

                  F_thetaphi = - sqrt(2.0d0) * mass1                    &
     &                         * a_ang_par                              &
     &                         * mass1                                  &
     &                         * sinh(2.0d0 * sen_alpha)                &
     &                         * r_BL                                   &
     &                         * (r_BL**2 + a_ang_par**2 + 2.0d0        &
     &                           * mass1 * r_BL * (sinh(sen_alpha)**2)) &
     &                         * sin_theta * cos_theta                  &
     &                         /rho_sq**2    

               E_up_tmp1 = -   chi                                          &
     &                       / alpha                                        &
     &                       / CC_tothird                                   &
     &                       / rbar                                         &
     &                       * (     sin_theta                              &
     &                             * sqrt(Delta)                            &
     &                             * ( F_tr + betau_phi * F_rphi )          &
     &                           +   cos_theta                              &
     &                             * ( F_ttheta + betau_phi * F_thetaphi )  &
     &                         )

               E_up_tmp2 = -   chi                                          &
     &                       / alpha                                        &
     &                       / CC_tothird                                   &
     &                       / rbar                                         &
     &                       * (     cos_theta                              &
     &                             * sqrt(Delta)                            &
     &                             * ( F_tr + betau_phi * F_rphi )          &
     &                           -   sin_theta                              &
     &                             * ( F_ttheta + betau_phi * F_thetaphi )  &
     &                         )

               u(H_EX)%d(i,j,k) = cos_phi * E_up_tmp1
               u(H_EY)%d(i,j,k) = sin_phi * E_up_tmp1
               u(H_EZ)%d(i,j,k) =           E_up_tmp2

               B_up_tmp1 =   chi**(1.5d0)                           &
     &                     / rbar                                   &
     &                     / rhobar                                 &
     &                     * (   sin_theta * sqrt(Delta) * F_rphi   &
     &                         + cos_theta * F_thetaphi             &
     &                       )

               B_up_tmp2 = -   chi**(1.5d0)                           &
     &                       / rbar                                   &
     &                       / rhobar                                 &
     &                       * (   cos_theta * sqrt(Delta) * F_rphi   &
     &                           - sin_theta * F_thetaphi             &
     &                         )

               u(H_BX)%d(i,j,k) = cos_phi * B_up_tmp2
               u(H_BY)%d(i,j,k) = sin_phi * B_up_tmp2
               u(H_BZ)%d(i,j,k) =           B_up_tmp1


               u(H_PSI_EM)%d(i,j,k) = 0.d0
               u(H_PHI_EM)%d(i,j,k) = 0.d0

             enddo
          enddo
       enddo !End Rasheed-Larsen

    else if (nint(par(P_GEOMETRY_TYPE)) .eq. 15) then
    ! Simple binary blackhole in EMDA. We add additional field into
    ! moving puncture form

    write(0,*)' initial data: Simple binary black hole in EMDA'

       do k=1,nz
          do j=1,ny
             do i=1,nx
                
                x1 = xx(i,j,k) - bh1x
                y1 = yy(i,j,k) - bh1y
                z1 = zz(i,j,k) - bh1z
                
                rv1 = sqrt(x1**2+y1**2+z1**2)
                if (rv1.lt.0.0001d0) then
                   rv1=0.0001d0
                endif
                vn1(1) = x1/rv1
                vn1(2) = y1/rv1
                vn1(3) = z1/rv1
                
                if (n_bhs.eq.2) then
                   x2 = xx(i,j,k) - bh2x
                   y2 = yy(i,j,k) - bh2y
                   z2 = zz(i,j,k) - bh2z
                   
                   rv2 = sqrt(x2**2+y2**2+z2**2)
                   if (rv2.lt.0.0001d0) then
                      rv2=0.0001d0
                   endif
                   vn2(1) = x2/rv2
                   vn2(2) = y2/rv2
                   vn2(3) = z2/rv2
                endif


                ! Brill-Lindquist conformal factor:
                vpsibl = 1.d0 + mass1/(2.d0*rv1)
                if (n_bhs.eq.2) then
                   vpsibl = vpsibl + mass2/(2.d0*rv2)
                endif

! hard-code in for the time being:
                i_better_u = 1
                if (i_better_u.eq.1) then
                   v_u_corr = 0.d0
                   
                   if (abs(spin1).gt.0.000001) then
                      amp_capj = 4.d0*spin1/mass1**2
                      amp_capr = 2.d0*rv1/mass1
                      l_r = 1.d0/(1.d0 + amp_capr)
                      u0_j = (l_r + l_r**2 + l_r**3 - 4.d0*l_r**4 + 2.d0*l_r**5)/40.d0
                      u2_j = -l_r**5/20.d0
                      mu_j =         vn1(1)*vs1(1)
                      mu_j =  mu_j + vn1(2)*vs1(2)
                      mu_j = (mu_j + vn1(3)*vs1(3))/abs(spin1)
                      p2_mu_j = (3.d0*mu_j**2 - 1.d0)/2.d0
                      v_u_j1 = amp_capj**2*(u0_j+u2_j*amp_capr**2*p2_mu_j)
                      v_u_corr = v_u_corr + v_u_j1
                   endif
                   
                   if (vp1tot.gt.0.000001) then
                      amp_capp = 2.d0*vp1tot/mass1
                      amp_capr = 2.d0*rv1/mass1
                      l_r = 1.d0/(1.d0 + amp_capr)
                      u0_p = l_r - 2.d0*l_r**2 + 2.d0*l_r**3 
                      u0_p = (u0_p - l_r**4 + 0.2d0*l_r**5)*(5.d0/32.d0)
                      u2_p = 15.d0*l_r + 132.d0*l_r**2 + 53.d0*l_r**3
                      u2_p = u2_p + 96.d0*l_r**4 + 82.d0*l_r**5
                      u2_p = u2_p + (84.d0/amp_capr)*(l_r**5+log(l_r)/amp_capr)
                      u2_p = (u2_p)/(80.d0*amp_capr)
                      mu_p =        vn1(1)*vp1(1)/vp1tot
                      mu_p = mu_p + vn1(2)*vp1(2)/vp1tot
                      mu_p = mu_p + vn1(3)*vp1(3)/vp1tot
                      p2_mu_p = (3.d0*mu_p**2 - 1.d0)/2.d0
                      v_u_p1 = amp_capp**2*(u0_p+u2_p*p2_mu_p)
                      v_u_corr = v_u_corr + v_u_p1
                   endif
                   
                   if ((vp1tot.gt.0.000001).and.(abs(spin1).gt.0.000001)) then
                      v1 =      (vp1(2)*vs1(3)-vp1(3)*vs1(2))*vn1(1)
                      v1 = v1 + (vp1(3)*vs1(1)-vp1(1)*vs1(3))*vn1(2)
                      v1 = v1 + (vp1(1)*vs1(2)-vp1(2)*vs1(1))*vn1(3)
                      v1 = v1*(16.d0/mass1**4)*rv1
                      
                      amp_capr = 2.d0*rv1/mass1
                      l_r = 1.d0/(1.d0 + amp_capr)
                      
                      v2 = 1.d0 + 5.d0*amp_capr + 10.d0*amp_capr**2
                      
                      v_u_c1 = (v1*v2*l_r**5)/80.d0
                      v_u_corr = v_u_corr + v_u_c1
                   endif
                   
                   if (n_bhs.eq.2) then
                      
                      if (abs(spin2).gt.0.000001) then
                         amp_capj = 4.d0*spin2/mass2**2
                         amp_capr = 2.d0*rv2/mass2
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_j = (l_r + l_r**2 + l_r**3 - 4.d0*l_r**4 + 2.d0*l_r**5)/40.d0
                         u2_j = -l_r**5/20.d0
                         mu_j =         vn2(1)*vs2(1)
                         mu_j =  mu_j + vn2(2)*vs2(2)
                         mu_j = (mu_j + vn2(3)*vs2(3))/abs(spin2)
                         p2_mu_j = (3.d0*mu_j**2 - 1.d0)/2.d0
                         v_u_j2 = amp_capj**2*(u0_j+u2_j*amp_capr**2*p2_mu_j)
                         v_u_corr = v_u_corr + v_u_j2
                      endif
                      
                      if (vp2tot.gt.0.000001) then
                         amp_capp = 2.d0*vp2tot/mass2
                         amp_capr = 2.d0*rv2/mass2
                         l_r = 1.d0/(1.d0 + amp_capr)
                         u0_p = l_r - 2.d0*l_r**2 + 2.d0*l_r**3 
                         u0_p = (u0_p - l_r**4 + 0.2d0*l_r**5)*(5.d0/32.d0)
                         u2_p = 15.d0*l_r + 132.d0*l_r**2 + 53.d0*l_r**3
                         u2_p = u2_p + 96.d0*l_r**4 + 82.d0*l_r**5
                         u2_p = u2_p + (84.d0/amp_capr)*(l_r**5+log(l_r)/amp_capr)
                         u2_p = (u2_p)/(80.d0*amp_capr)
                         mu_p =        vn2(1)*vp2(1)/vp2tot
                         mu_p = mu_p + vn2(2)*vp2(2)/vp2tot
                         mu_p = mu_p + vn2(3)*vp2(3)/vp2tot
                         p2_mu_p = (3.d0*mu_p**2 - 1.d0)/2.d0
                         v_u_p2 = amp_capp**2*(u0_p+u2_p*p2_mu_p)
                         v_u_corr = v_u_corr + v_u_p2
                      endif
                      
                      if ((vp2tot.gt.0.000001).and.(abs(spin2).gt.0.000001)) then
                         v1 =      (vp2(2)*vs2(3)-vp2(3)*vs2(2))*vn2(1)
                         v1 = v1 + (vp2(3)*vs2(1)-vp2(1)*vs2(3))*vn2(2)
                         v1 = v1 + (vp2(1)*vs2(2)-vp2(2)*vs2(1))*vn2(3)
                         v1 = v1*(16.d0/mass2**4)*rv2
                         
                         amp_capr = 2.d0*rv2/mass2
                         l_r = 1.d0/(1.d0 + amp_capr)
                         
                         v2 = 1.d0 + 5.d0*amp_capr + 10.d0*amp_capr**2
                         
                         v_u_c2 = (v1*v2*l_r**5)/80.d0
                         v_u_corr = v_u_corr + v_u_c2
                      endif

                   endif
                  
! vpsibl_u will be used for the conformal factor, 
                   vpsibl_u  = vpsibl + v_u_corr
! vpsibl_u2 is for the Aij terms...
! since the corrections are first order...
! adding half of the correction seems to give the best results...
! update - do a fit for spin = 0.6...
!                   vpsibl_u2 = vpsibl
!                   vpsibl_u2 = vpsibl + v_u_corr/2.d0
!                   vpsibl_u2 = vpsibl + 0.313d0*v_u_corr
                   vpsibl_u2 = vpsibl + v_u_corr
                else
                   vpsibl_u  = vpsibl
                   vpsibl_u2 = vpsibl
                endif
                
                u(H_ALPHA)%d(i,j,k) = 1.d0/vpsibl_u**2
                v2 =  1.d0/vpsibl_u**4

                if (setchi) then
                   u(H_CHI)%d(i,j,k) = v2
                end if
                u(H_TRK)%d(i,j,k) =  0.0d0
                
                u(H_SHIFT1)%d(i,j,k) =  0.0d0
                u(H_SHIFT2)%d(i,j,k) =  0.0d0
                u(H_SHIFT3)%d(i,j,k) =  0.0d0
                   
                u(H_GAM1)%d(i,j,k) =  0.0d0
                u(H_GAM2)%d(i,j,k) =  0.0d0
                u(H_GAM3)%d(i,j,k) =  0.0d0
                
                u(H_GB1)%d(i,j,k) =  0.0d0
                u(H_GB2)%d(i,j,k) =  0.0d0
                u(H_GB3)%d(i,j,k) =  0.0d0
                
                u(H_GT11)%d(i,j,k) =  1.0d0
                u(H_GT12)%d(i,j,k) =  0.0d0
                u(H_GT13)%d(i,j,k) =  0.0d0
                u(H_GT22)%d(i,j,k) =  1.0d0
                u(H_GT23)%d(i,j,k) =  0.0d0
                u(H_GT33)%d(i,j,k) =  1.0d0
                

                do i1 = 1,3
                   do i2 = i1,3
                      
                      v2 = 0.d0
                      do i3 = 1,3
                         do i4 = 1,3
                            vt1 = epijk(i1,i3,i4)*vs1(i3)*vn1(i4)*vn1(i2)
                            vt2 = epijk(i2,i3,i4)*vs1(i3)*vn1(i4)*vn1(i1)
                            v2 = v2 + vt1 + vt2
                         enddo
                      enddo
                      
                      v3 = vp1(i1)*vn1(i2) + vp1(i2)*vn1(i1)
                      vt1 = 0.d0
                      do i3 = 1,3
                         vt1 = vt1 + vp1(i3)*vn1(i3)
                      enddo
                      vt1 = vt1*(vn1(i1)*vn1(i2) - deltaij(i1,i2))
                      v3 = v3 + vt1
                      
                      v1 = 3.d0/(vpsibl_u2**6*rv1**3)
                      v4 = v1*(v2+(rv1/2.d0)*v3)
                      
                      if (n_bhs.eq.2) then
                         v2 = 0.d0
                         do i3 = 1,3
                            do i4 = 1,3
                               vt1 = epijk(i1,i3,i4)*vs2(i3)*vn2(i4)*vn2(i2)
                               vt2 = epijk(i2,i3,i4)*vs2(i3)*vn2(i4)*vn2(i1)
                               v2 = v2 + vt1 + vt2
                            enddo
                         enddo
                         
                         v3 = vp2(i1)*vn2(i2) + vp2(i2)*vn2(i1)
                         vt1 = 0.d0
                         do i3 = 1,3
                            vt1 = vt1 + vp2(i3)*vn2(i3)
                         enddo
                         vt1 = vt1*(vn2(i1)*vn2(i2) - deltaij(i1,i2))
                         v3 = v3 + vt1
                         
                         v1 = 3.d0/(vpsibl_u2**6*rv2**3)
                         v4 = v4 + v1*(v2+(rv2/2.d0)*v3)
                      endif
                      
                      if ((i1.eq.1).and.(i2.eq.1)) then                            
                         u(H_A11)%d(i,j,k) = v4
                      else if ((i1.eq.1).and.(i2.eq.2)) then
                         u(H_A12)%d(i,j,k) = v4
                      else if ((i1.eq.1).and.(i2.eq.3)) then
                         u(H_A13)%d(i,j,k) = v4
                      else if ((i1.eq.2).and.(i2.eq.2)) then
                         u(H_A22)%d(i,j,k) = v4
                      else if ((i1.eq.2).and.(i2.eq.3)) then
                         u(H_A23)%d(i,j,k) = v4
                      else if ((i1.eq.3).and.(i2.eq.3)) then
                         u(H_A33)%d(i,j,k) = v4
                      endif
                   enddo
                enddo
               
! and add the EM initial data...

                if (i_btype.eq.1) then
                   u(H_EX)%d(i,j,k) = 0.d0
                   u(H_EY)%d(i,j,k) = 0.d0 
                   u(H_EZ)%d(i,j,k) = 0.d0
                   u(H_BX)%d(i,j,k) = 0.d0 
                   u(H_BY)%d(i,j,k) = 0.d0
                   u(H_BZ)%d(i,j,k) = v_bamp/((vpsibl)**6)+v_bamp*1.d-6
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 
                elseif (i_btype.eq.2) then
                   ! Electrically charged BH, exact for single non-boosted 
                   u(H_BX)%d(i,j,k) = 0.0
                   u(H_BY)%d(i,j,k) = 0.0
                   u(H_BZ)%d(i,j,k) = 0.0
                   ! radial electric field E^r = Q/(r^2*psi^6)
                   ! eq (3.16) from http://arxiv.org/abs/0907.1151
                   u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 

                   if (n_bhs.eq.2) then
                     u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                   end if
                elseif (i_btype.eq.3) then
                   ! Magnetic monopole BH, exact for single non-boosted 
                   u(H_EX)%d(i,j,k) = 0.0
                   u(H_EY)%d(i,j,k) = 0.0
                   u(H_EZ)%d(i,j,k) = 0.0
                   ! radial electric field E^r = Q/(r^2*psi^6)
                   ! eq (3.16) from http://arxiv.org/abs/0907.1151
                   u(H_BX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_BY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_BZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 

                   if (n_bhs.eq.2) then
                     u(H_BX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_BY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_BZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                   end if

                elseif (i_btype.eq.4) then
                   ! DYonic BH:
                   ! Magnetic monopole BH, exact for single non-boosted 
                   u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   ! radial electric field E^r = Q/(r^2*psi^6)
                   ! eq (3.16) from http://arxiv.org/abs/0907.1151
                   u(H_BX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_BY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_BZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 

                   if (n_bhs.eq.2) then
                     u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                     u(H_BX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_BY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_BZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2) &
     &                                + v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                   end if
                elseif (i_btype.eq.5) then
                   ! electric BH with magnetic BH:
                   ! Electrically charged BH, exact for single non-boosted 
                   u(H_BX)%d(i,j,k) = 0.0
                   u(H_BY)%d(i,j,k) = 0.0
                   u(H_BZ)%d(i,j,k) = 0.0
                   ! radial electric field E^r = Q/(r^2*psi^6)
                   ! eq (3.16) from http://arxiv.org/abs/0907.1151
                   u(H_EX)%d(i,j,k) = v_e1amp*vn1(1)/(vpsibl_u**6*rv1**2)
                   u(H_EY)%d(i,j,k) = v_e1amp*vn1(2)/(vpsibl_u**6*rv1**2)
                   u(H_EZ)%d(i,j,k) = v_e1amp*vn1(3)/(vpsibl_u**6*rv1**2)
                   u(H_PSI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHI_EM)%d(i,j,k) = 0.d0 
                   u(H_PHIR)%d(i,j,k) = v_sfamp-v_e1amp*exp(-rv1**2)

                   if (n_bhs.eq.2) then
                     u(H_BX)%d(i,j,k) = v_e2amp*vn2(1)/(vpsibl_u**6*rv2**2)
                     u(H_BY)%d(i,j,k) = v_e2amp*vn2(2)/(vpsibl_u**6*rv2**2)
                     u(H_BZ)%d(i,j,k) = v_e2amp*vn2(3)/(vpsibl_u**6*rv2**2) 
                     u(H_PHIR)%d(i,j,k) = u(H_PHIR)%d(i,j,k)+v_e2amp*exp(-rv2**2)
                   end if

                endif

! and add the SF initial data...
! Axion field is added like dilaton field

                if (i_btype.ne.5) then
                   u(H_PHIR)%d(i,j,k) = v_sfamp
                   u(H_PHII)%d(i,j,k) = v_sfamp_axn 
                end if
                u(H_PIR)%d(i,j,k) = 0.0d0
                u(H_PII)%d(i,j,k) = 0.0d0

                
             enddo
          enddo
       enddo !End BBH

    else
       write(0,*)'BSSN initial: unkown geometry_type ',&
                &   nint(par(P_GEOMETRY_TYPE))
       write(0,*)'BSSN+MHD requires geometry_type = 0'
      stop
    end if 

    call EnforceBSSN(u, v, w, par)


    return
  end subroutine

!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
  subroutine ks_initial_data(gtd, Atd, Chi, TrK, Gamt, alpha, Betau, &
           &                 x, y, z)
    implicit none
    CCTK_REAL, dimension(3,3) :: gtd, Atd
    CCTK_REAL, dimension(3)   :: Gamt, Betau
    CCTK_REAL                 :: Chi, TrK, alpha, x, y, z

    ! local
    CCTK_REAL, parameter      :: M = 1.0d0
      real*8  t1, t2, t3, t4, t5
      real*8  t6, t7, t8, t16, t17
      real*8  t20, t22, t23, t26, t29
      real*8  t38, t10, t12, t13, t15
      real*8  t18, t21, t34, t37, t39
      real*8  t44, t45, t49, t52, t53
      real*8  t55, t60, t65, t69, t92
      real*8  t11, t27, t28, t31, t33
      real*8  t36



      t1 = x ** 2
      t2 = y ** 2
      t3 = z ** 2
      t4 = t1 + t2 + t3
      t5 = sqrt(t4)
      t6 = t1 * t5
      t7 = t2 * t5
      t8 = t3 * t5
      t16 = (0.1D1 / t5 * (t5 + 0.2D1 * M)) ** (0.1D1 / 0.3D1)
      t17 = 0.1D1 / t16
      t20 = 0.1D1 / t5 / t4
      t22 = M * t17
      t23 = x * t20
      t26 = 0.2D1 * y * t23 * t22
      t29 = 0.2D1 * z * t23 * t22
      t38 = 0.2D1 * z * y * t20 * t22
      gtd(1,1) = t20 * t17 * (t6 + t7 + t8 + 0.2D1 * t1 * M)
      gtd(1,2) = t26
      gtd(1,3) = t29
      gtd(2,1) = t26
      gtd(2,2) = t20 * t17 * (t6 + t7 + t8 + 0.2D1 * t2 * M)
      gtd(2,3) = t38
      gtd(3,1) = t29
      gtd(3,2) = t38
      gtd(3,3) = t20 * t17 * (t6 + t7 + t8 + 0.2D1 * t3 * M)
      t1 = x ** 2
      t2 = y ** 2
      t3 = z ** 2
      t4 = t1 + t2 + t3
      t5 = sqrt(t4)
      t7 = 0.1D1 / t5 * M
      t10 = log(0.1D1 + 0.2D1 * t7)
      t12 = exp(-0.3333333333333333D0 * t10)
      t13 = t4 ** 2
      t15 = 0.1D1 / t13 * M
      t18 = 0.1D1 / (t5 + 0.2D1 * M)
      t20 = sqrt(t18 * t5)
      t21 = t7 + 0.2D1
      t29 = 0.1D1 / t5 / t4 * M
      t34 = 0.1D1 / t4
      t37 = t5 + 0.3D1 * M
      t39 = t18 * t37 * t20
      t44 = t20 * t15
      t45 = x * t21
      t49 = M ** 2
      t52 = 0.1D1 / t5 / t13 * t49
      t53 = x * t52
      t55 = t18 * t37
      t60 = (-0.2D1 * y * t45 * t44 - 0.1333333333333333D1 * t55 * t20 *&
     & y * t53) * t12
      t65 = t55 * t20 * z
      t69 = (-0.2D1 * z * t45 * t44 - 0.1333333333333333D1 * t65 * t53) &
     &* t12
      t92 = (-0.2D1 * z * y * t21 * t44 - 0.1333333333333333D1 * t65 * y&
     & * t52) * t12
      Atd(1,1) = (-0.2D1 * (t1 * t21 - t1 - t2 - t3) * t20 * t15 - 0.666&
     &6666666666666D0 * t39 * t34 * M * (0.1D1 + 0.2D1 * t1 * t29)) * t1&
     &2
      Atd(1,2) = t60
      Atd(1,3) = t69
      Atd(2,1) = t60
      Atd(2,2) = (-0.2D1 * (t2 * t21 - t1 - t2 - t3) * t20 * t15 - 0.666&
     &6666666666666D0 * t39 * t34 * M * (0.1D1 + 0.2D1 * t2 * t29)) * t1&
     &2
      Atd(2,3) = t92
      Atd(3,1) = t69
      Atd(3,2) = t92
      Atd(3,3) = (-0.2D1 * (t3 * t21 - t1 - t2 - t3) * t20 * t15 - 0.666&
     &6666666666666D0 * t39 * t34 * M * (0.1D1 + 0.2D1 * t3 * t29)) * t1&
     &2
      t1 = x ** 2
      t2 = y ** 2
      t3 = z ** 2
      t5 = sqrt(t1 + t2 + t3)
      t10 = (0.1D1 + 0.2D1 / t5 * M) ** (0.1D1 / 0.3D1)
      Chi = 0.1D1 / t10
      t1 = x ** 2
      t2 = y ** 2
      t3 = z ** 2
      t4 = t1 + t2 + t3
      t7 = sqrt(t4)
      t10 = 0.1D1 / (t7 + 0.2D1 * M)
      t12 = sqrt(t7 * t10)
      TrK = 0.2D1 * M / t4 * t12 * (t7 + 0.3D1 * M) * t10
      t1 = x ** 2
      t4 = y ** 2
      t5 = z ** 2
      t6 = t1 + t4 + t5
      t7 = sqrt(t6)
      t11 = M ** 2
      t18 = 0.5D1 * t1 * M + t1 * t7 + t4 * t7 + t5 * t7 + 0.6D1 * t7 * &
     &t11 + 0.5D1 * t5 * M + 0.5D1 * t4 * M
      t22 = t7 + 0.2D1 * M
      t23 = t22 ** 2
      t27 = (0.1D1 / t7 * t22) ** (0.1D1 / 0.3D1)
      t28 = t27 ** 2
      t31 = t6 ** 2
      t33 = 0.1D1 / t31 / t28 / t23
      t36 = M * t18
      Gamt(1) = 0.8D1 / 0.3D1 * t33 * M * x * t18
      Gamt(2) = 0.8D1 / 0.3D1 * t33 * y * t36
      Gamt(3) = 0.8D1 / 0.3D1 * t33 * z * t36
      t1 = x ** 2
      t2 = y ** 2
      t3 = z ** 2
      t5 = sqrt(t1 + t2 + t3)
      alpha = sqrt(t5 / (t5 + 0.2D1 * M))
      t1 = x ** 2
      t2 = y ** 2
      t3 = z ** 2
      t5 = sqrt(t1 + t2 + t3)
      t7 = 0.1D1 / t5 * M
      t10 = 0.1D1 / (t5 + 0.2D1 * M)
      Betau(1) = 0.2D1 * t7 * t10 * x
      Betau(2) = 0.2D1 * t10 * y * t7
      Betau(3) = 0.2D1 * t10 * z * t7

    return
  end subroutine


  ! Single Kerr black hole in Kerr-Schild coordinates.  Spin along z axis.
  subroutine kskerr_initial_data(gtd, Atd, Chi, TrK, Gamt, alpha, Betau, &
       &                 x, y, z, m, ain)
    implicit none
    CCTK_REAL, dimension(3,3) :: gtd, Atd
    CCTK_REAL, dimension(3)   :: Gamt, Betau
    CCTK_REAL                 :: Chi, TrK, alpha, x, y, z, M, ain

    real*8 :: r, rbar, r1, r2, detgamma, traceA, a, af, re, rplus, r1e, r2e
    real*8 :: interpconst, r3e, ractual, reactual

    a = ain
    r =         Sqrt(-a**2 + x**2 + y**2 + z**2 + &
     &    Sqrt(4*a**2*z**2 + (-a**2 + x**2 + y**2 + z**2)**2))/Sqrt(2.D0)
    re = Sqrt(x**2+y**2+z**2)
    rplus = (M + sqrt(M**2 - a**2))

    ! Store the actual values of r, re
    ractual = r
    reactual = re

    ! For r < r2 stuff the BH by modifying all the r/s that appear in
    ! the denominators of the initial data functions so that no
    ! singularities occur (call it rbar).  r1 is the smallest value
    ! rbar is allowed to take.
    r1 = 0.0D0 * rplus
    r2 = 0.7D0 * rplus

    r1e = Sqrt(r1**2+a**2)*r1*re / Sqrt(r1**2*(x**2+y**2) + (r1**2+a**2)*z**2 )
    r2e = Sqrt(r2**2+a**2)*r2*re / Sqrt(r2**2*(x**2+y**2) + (r2**2+a**2)*z**2 )

    if (r < r2 .AND. re > r1) then
       a = (ain*(re - r1)**6*((r1 - r2e)**5 + &
            &      (re - r2e)*(6*(r1 - r2e)**4 + &
            &         7*(re - r2e)*(36*re**3 + 3*r1**3 + 18*re**2*(r1 - 7*r2e) - &
            &            17*r1**2*r2e + 43*r1*r2e**2 - 65*r2e**3 + &
            &            4*re*(2*r1**2 - 13*r1*r2e + 38*r2e**2)))))/(r1 - r2e)**11
    else if (re < r1) then
       a = 0.
    end if

    !r =         Sqrt(-a**2 + x**2 + y**2 + z**2 + &
    !     &    Sqrt(4*a**2*z**2 + (-a**2 + x**2 + y**2 + z**2)**2))/Sqrt(2.D0)

    !r1 = 0.1D0 * (M + sqrt(M**2.D0 - a**2.D0))
    !r2 = 0.5D0 * (M + sqrt(M**2.D0 - a**2.D0))


    r1 = 0.4D0 * rplus
    r2 = 0.7D0 * rplus

    r1e = Sqrt(r1**2+a**2)*r1*re / Sqrt(r1**2*(x**2+y**2) + (r1**2+a**2)*z**2 )
    r2e = Sqrt(r2**2+a**2)*r2*re / Sqrt(r2**2*(x**2+y**2) + (r2**2+a**2)*z**2 )

    interpconst = 0.15
    r3e = r1e - 0.2*M
    if (r < r2) then
       ! First interpolate a between 0 and ain.
       !a = (ain*r**6*(-252*r**5 + 1386*r**4*r2 - 3080*r**3*r2**2 + 3465*r**2*r2**3 - &
       !     &      1980*r*r2**4 + 462*r2**5))/r2**11
       !if (r > r1) then
       !  a = (ain*(r - r1)**6*((r1 - r2)**5 + &
       !     &      (r - r2)*(6*(r1 - r2)**4 + &
       !     &         7*(r - r2)*(36*r**3 + 3*r1**3 + 18*r**2*(r1 - 7*r2) - &
       !     &            17*r1**2*r2 + 43*r1*r2**2 - 65*r2**3 + &
       !     &            4*r*(2*r1**2 - 13*r1*r2 + 38*r2**2)))))/(r1 - r2)**11
       !else
       !  a = 0.
       !end if 
       ! Re-define r to match that new a.
       !r =         Sqrt(-a**2 + x**2 + y**2 + z**2 + &
       !     &    Sqrt(4*a**2*z**2 + (-a**2 + x**2 + y**2 + z**2)**2))/Sqrt(2.D0)
       ! Rescale r.  This might increase rs value from previously, 
       ! but I do not think it will be a problem.
       ! C5 version
       rbar = (126*r**11*(2*r1 - r2) + 35*r**9*(88*r1 - 45*r2)*r2**2 + &
            &    30*r**7*(66*r1 - 35*r2)*r2**4 + r1*r2**11 + &
            &    42*r**6*r2**5*(-11*r1 + 6*r2) + 45*r**8*r2**3*(-77*r1 + 40*r2) + &
            &    14*r**10*r2*(-99*r1 + 50*r2))/r2**11
       ! add a bit to rbar to smooth it out
       if (re > r3e) then
          rbar = rbar + interpconst - (interpconst*(r3e - re)**5* &
               &     (-56*r2e**3 + r3e**3 + 5*r3e**2*re + 15*r3e*re**2 + 35*re**3 + &
               &       28*r2e**2*(r3e + 5*re) - 8*r2e*(r3e**2 + 5*r3e*re + 15*re**2)))/ &
               &   (r2e - r3e)**8
       else
          rbar = rbar + interpconst
       end if
    else
       rbar = r
    end if

    ! Try this:
    r = rbar

    gtd(1,1) =         (r**3*(r*((a**2 + r**2)**2 + 2*M*r*x**2) + 4*a*M*r*x*y + &
     &       2*a**2*M*y**2) + a**2*(a**2 + r**2)**2*z**2)/&
     &  ((a**2 + rbar**2)**2*(rbar**4 + a**2*z**2)**0.6666666666666666*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333)
    gtd(1,2) =         (2*M*r**3*(r*x + a*y)*(-(a*x) + r*y))/&
     &  ((a**2 + rbar**2)**2*(rbar**4 + a**2*z**2)**0.6666666666666666*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333)
    gtd(1,3) =         (2*M*r**2*(r*x + a*y)*z)/&
     &  ((a**2 + rbar**2)*(rbar**4 + a**2*z**2)**0.6666666666666666*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333)
    gtd(2,1) =         (2*M*r**3*(r*x + a*y)*(-(a*x) + r*y))/&
     &  ((a**2 + rbar**2)**2*(rbar**4 + a**2*z**2)**0.6666666666666666*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333)
    gtd(2,2) =         (r**3*(a**4*r + r**5 + 2*a**2*(r**3 + M*x**2) - 4*a*M*r*x*y + &
     &       2*M*r**2*y**2) + a**2*(a**2 + r**2)**2*z**2)/&
     &  ((a**2 + rbar**2)**2*(rbar**4 + a**2*z**2)**0.6666666666666666*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333)
    gtd(2,3) =         (2*M*r**2*(-(a*x) + r*y)*z)/&
     &  ((a**2 + rbar**2)*(rbar**4 + a**2*z**2)**0.6666666666666666*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333)
    gtd(3,1) =         (2*M*r**2*(r*x + a*y)*z)/&
     &  ((a**2 + rbar**2)*(rbar**4 + a**2*z**2)**0.6666666666666666*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333)
    gtd(3,2) =         (2*M*r**2*(-(a*x) + r*y)*z)/&
     &  ((a**2 + rbar**2)*(rbar**4 + a**2*z**2)**0.6666666666666666*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333)
    gtd(3,3) =         (r**4 + (a**2 + 2*M*r)*z**2)/&
     &  ((rbar**4 + a**2*z**2)**0.6666666666666666*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333)

    Atd(1,1) =         (2*M*r**2*(r**10*(-(a**4*r*(3*M + r)) + &
     &         r**2*(3*M + 2*r)*(r**3 - (4*M + 3*r)*x**2) - &
     &         3*a**3*(2*M + r)*x*y - &
     &         a*r*(24*M**2 + 28*M*r + 9*r**2)*x*y + &
     &         a**2*(r**4 - M*(12*M + 5*r)*y**2)) + &
     &      a**2*r**6*(r*(4*M + 3*r)*(-a**4 + a**2*r**2 + 2*r**4) + &
     &         4*r*(3*a**2*(2*M + r) + M*r*(M + 2*r))*x**2 + &
     &         2*a*(3*a**2*(3*M + r) - r*(-4*M**2 + M*r + 3*r**2))*x*&
     &          y + 4*a**2*M*(M - r)*y**2)*z**2 + &
     &      a**4*r**3*(-(a**4*(M + 3*r)) + r**4*(5*M + 6*r) + &
     &         r**2*(M + 6*r)*x**2 + 9*a**3*x*y + a*r*(2*M + 3*r)*x*y + &
     &         a**2*(3*(r**3 + 4*r*x**2) + M*(4*r**2 + y**2)))*z**4 + &
     &      a**6*(-a**4 + a**2*r**2 + 2*r**4)*z**6))/&
     &  (3.*(a**2 + rbar**2)**2*&
     &    (rbar**4 + a**2*z**2)**2.1666666666666665*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.8333333333333333)
    Atd(1,2) =         (M*r**5*(-2*r**9*(3*M + 2*r)*(4*M + 3*r)*x*y + &
     &      a*r**8*(24*M**2 + 28*M*r + 9*r**2)*(x - y)*(x + y) - &
     &      2*a**6*(M - 12*r)*x*y*z**4 + 9*a**7*(-x**2 + y**2)*z**4 + &
     &      2*a**2*M*r**5*x*y*(r**2*(12*M + 5*r) + 4*(M + 2*r)*z**2) - &
     &      a**5*r*(x - y)*(x + y)*z**2*&
     &       (6*r**2*(3*M + r) + (2*M + 3*r)*z**2) + &
     &      2*a**4*r**2*x*y*z**2*&
     &       (4*r*(-M**2 + 7*M*r + 3*r**2) + (M + 6*r)*z**2) + &
     &      a**3*r**4*(x - y)*(x + y)*&
     &       (3*r**3*(2*M + r) + 2*(-4*M**2 + M*r + 3*r**2)*z**2)))/&
     &  (3.*(a**2 + rbar**2)**2*&
     &    (rbar**4 + a**2*z**2)**2.1666666666666665*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.8333333333333333)
    Atd(1,3) =         (M*r**3*z*(r**7*(-2*r*&
     &          (6*a**2*(2*M + r) + r*(3*M + 2*r)*(4*M + 3*r))*x - &
     &         a*(9*a**2*(2*M + r) + r*(24*M**2 + 28*M*r + 9*r**2))*y)&
     &       + 2*a**2*r**3*(4*M*r*(3*a**2 + r*(M + 2*r))*x + &
     &         a*(M - r)*(3*a**2 + r*(4*M + 3*r))*y)*z**2 + &
     &      a**4*(2*r*(6*a**2 + r*(M + 6*r))*x + &
     &         a*(3*a**2 + r*(2*M + 3*r))*y)*z**4))/&
     &  (3.*(a**2 + rbar**2)*(rbar**4 + a**2*z**2)**2.1666666666666665*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.8333333333333333)
    Atd(2,1) =         (M*r**5*(-2*r**9*(3*M + 2*r)*(4*M + 3*r)*x*y + &
     &      a*r**8*(24*M**2 + 28*M*r + 9*r**2)*(x - y)*(x + y) - &
     &      2*a**6*(M - 12*r)*x*y*z**4 + 9*a**7*(-x**2 + y**2)*z**4 + &
     &      2*a**2*M*r**5*x*y*(r**2*(12*M + 5*r) + 4*(M + 2*r)*z**2) - &
     &      a**5*r*(x - y)*(x + y)*z**2*&
     &       (6*r**2*(3*M + r) + (2*M + 3*r)*z**2) + &
     &      2*a**4*r**2*x*y*z**2*&
     &       (4*r*(-M**2 + 7*M*r + 3*r**2) + (M + 6*r)*z**2) + &
     &      a**3*r**4*(x - y)*(x + y)*&
     &       (3*r**3*(2*M + r) + 2*(-4*M**2 + M*r + 3*r**2)*z**2)))/&
     &  (3.*(a**2 + rbar**2)**2*&
     &    (rbar**4 + a**2*z**2)**2.1666666666666665*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.8333333333333333)
    Atd(2,2) =         (2*M*r**2*(r**10*(-(a**4*r*(3*M + r)) + &
     &         a**2*(r**4 - M*(12*M + 5*r)*x**2) + &
     &         3*a**3*(2*M + r)*x*y + &
     &         a*r*(24*M**2 + 28*M*r + 9*r**2)*x*y + &
     &         r**2*(3*M + 2*r)*(r**3 - (4*M + 3*r)*y**2)) + &
     &      a**2*r**6*(-(a**4*r*(4*M + 3*r)) + 2*r**5*(4*M + 3*r) - &
     &         6*a**3*(3*M + r)*x*y + &
     &         2*a*r*(-4*M**2 + M*r + 3*r**2)*x*y + &
     &         4*M*r**2*(M + 2*r)*y**2 + &
     &         a**2*(4*M*r**3 + 3*r**4 + 4*M**2*x**2 - 4*M*r*x**2 + &
     &            12*r*(2*M + r)*y**2))*z**2 + &
     &      a**4*r**3*(-(a**4*(M + 3*r)) + r**4*(5*M + 6*r) - &
     &         9*a**3*x*y - a*r*(2*M + 3*r)*x*y + r**2*(M + 6*r)*y**2 + &
     &         a**2*(M*(4*r**2 + x**2) + 3*(r**3 + 4*r*y**2)))*z**4 + &
     &      a**6*(-a**4 + a**2*r**2 + 2*r**4)*z**6))/&
     &  (3.*(a**2 + rbar**2)**2*&
     &    (rbar**4 + a**2*z**2)**2.1666666666666665*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.8333333333333333)
    Atd(2,3) =         (M*r**3*z*(r**7*(a*(9*a**2*(2*M + r) + &
     &            r*(24*M**2 + 28*M*r + 9*r**2))*x - &
     &         2*r*(6*a**2*(2*M + r) + r*(3*M + 2*r)*(4*M + 3*r))*y) + &
     &      2*a**2*r**3*(-(a*(M - r)*(3*a**2 + r*(4*M + 3*r))*x) + &
     &         4*M*r*(3*a**2 + r*(M + 2*r))*y)*z**2 + &
     &      a**4*(-(a*(3*a**2 + r*(2*M + 3*r))*x) + &
     &         2*r*(6*a**2 + r*(M + 6*r))*y)*z**4))/&
     &  (3.*(a**2 + rbar**2)*(rbar**4 + a**2*z**2)**2.1666666666666665*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.8333333333333333)
    Atd(3,1) =         (M*r**3*z*(r**7*(-2*r*&
     &          (6*a**2*(2*M + r) + r*(3*M + 2*r)*(4*M + 3*r))*x - &
     &         a*(9*a**2*(2*M + r) + r*(24*M**2 + 28*M*r + 9*r**2))*y)&
     &       + 2*a**2*r**3*(4*M*r*(3*a**2 + r*(M + 2*r))*x + &
     &         a*(M - r)*(3*a**2 + r*(4*M + 3*r))*y)*z**2 + &
     &      a**4*(2*r*(6*a**2 + r*(M + 6*r))*x + &
     &         a*(3*a**2 + r*(2*M + 3*r))*y)*z**4))/&
     &  (3.*(a**2 + rbar**2)*(rbar**4 + a**2*z**2)**2.1666666666666665*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.8333333333333333)
    Atd(3,2) =         (M*r**3*z*(r**7*(a*(9*a**2*(2*M + r) + &
     &            r*(24*M**2 + 28*M*r + 9*r**2))*x - &
     &         2*r*(6*a**2*(2*M + r) + r*(3*M + 2*r)*(4*M + 3*r))*y) + &
     &      2*a**2*r**3*(-(a*(M - r)*(3*a**2 + r*(4*M + 3*r))*x) + &
     &         4*M*r*(3*a**2 + r*(M + 2*r))*y)*z**2 + &
     &      a**4*(-(a*(3*a**2 + r*(2*M + 3*r))*x) + &
     &         2*r*(6*a**2 + r*(M + 6*r))*y)*z**4))/&
     &  (3.*(a**2 + rbar**2)*(rbar**4 + a**2*z**2)**2.1666666666666665*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.8333333333333333)
    Atd(3,3) =         (2*M*r**2*(r**11*(3*M + 2*r) - &
     &      r**7*(r*(3*M + 2*r)*(4*M + 3*r) + 2*a**2*(8*M + 3*r))*&
     &       z**2 + a**2*r**3*(a**2*(5*M - 6*r) + 4*M*r*(M + 2*r))*&
     &       z**4 + a**4*(2*a**2 + r*(M + 6*r))*z**6))/&
     &  (3.*(rbar**4 + a**2*z**2)**2.1666666666666665*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.8333333333333333)

    ! For Chi, also set r -> rbar in the numerator
    Chi =         (rbar**4 + a**2*z**2)**0.3333333333333333/&
     &  (2*M*rbar**3 + rbar**4 + a**2*z**2)**0.3333333333333333

    TrK =         (2*M*r**2*((r**4 + a**2*z**2)**2 + &
     &      M*(3*r**7 + a**2*r**3*z**2)))/&
     &  ((rbar**4 + a**2*z**2)**1.5*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.5)

    Gamt(1) =         (2*M*r**2*(r**10*(a**2*(2*M + r)*x + 4*r**2*(3*M + r)*x + &
     &         a*r*(10*M + 3*r)*y) + &
     &      a**2*r**6*(-2*a**2*(3*M + r)*x + r**2*(6*M + 7*r)*x + &
     &         3*a*r*(4*M + 3*r)*y)*z**2 + &
     &      a**4*r**3*(-3*a**2*x + 2*r*(M + 3*r)*x + a*(2*M + 9*r)*y)*&
     &       z**4 + 3*a**6*(r*x + a*y)*z**6))/&
     &  (3.*(a**2 + rbar**2)*(rbar**4 + a**2*z**2)**2.3333333333333335*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.6666666666666667)
    Gamt(2) =         (2*M*r**2*(r**10*(-(a*r*(10*M + 3*r)*x) + a**2*(2*M + r)*y + &
     &         4*r**2*(3*M + r)*y) - &
     &      a**2*r**6*(3*a*r*(4*M + 3*r)*x + 2*a**2*(3*M + r)*y - &
     &         r**2*(6*M + 7*r)*y)*z**2 + &
     &      a**4*r**3*(-(a*(2*M + 9*r)*x) - 3*a**2*y + 2*r*(M + 3*r)*y)*&
     &       z**4 + 3*a**6*(-(a*x) + r*y)*z**6))/&
     &  (3.*(a**2 + rbar**2)*(rbar**4 + a**2*z**2)**2.3333333333333335*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.6666666666666667)
    Gamt(3) =         (2*M*r*z*(4*r**11*(3*M + r) + &
     &      2*a**4*r**3*z**2*(r**2*(-M + r) + (M + 3*r)*z**2) + &
     &      a**6*(-(r**2*z**4) + 3*z**6) + &
     &      a**2*r**7*(3*r**3 + 7*r*z**2 + 6*M*(r**2 + z**2))))/&
     &  (3.*(rbar**4 + a**2*z**2)**2.3333333333333335*&
     &    (2*M*rbar**3 + rbar**4 + a**2*z**2)**1.6666666666666667)

    alpha = Sqrt(r**4 + a**2*z**2)/Sqrt(2*M*rbar**3 + rbar**4 + a**2*z**2)

    Betau(1) =         (2*M*r**3*(r*x + a*y))/&
     &  ((a**2 + rbar**2)*(2*M*rbar**3 + rbar**4 + a**2*z**2))
    Betau(2) =         (2*M*r**3*(-(a*x) + r*y))/&
     &  ((a**2 + rbar**2)*(2*M*rbar**3 + rbar**4 + a**2*z**2))
    Betau(3) = (2*M*r**2*z)/(2*M*rbar**3 + rbar**4 + a**2*z**2)

    ! If we have stuffed the BH, it may be the case that gtd no longer
    ! has unit determinant, or Atd is no longer traceless.  This will
    ! be fixed by calling the subroutine EnforceBSSN by subroutine
    ! initialdata.
       
    ! Set lapse=1 and shift=0
    alpha = 1
    Betau = 0

    ! Switch back to actual values of r, re
    r = ractual
    re = reactual

    ! Smooth out TrK -- reduce its value in the stuffed region
!!$    r1 = 0.3D0 * rplus
!!$    r2 = 0.7D0 * rplus
!!$
!!$    r1e = Sqrt(r1**2+a**2)*r1*re / Sqrt(r1**2*(x**2+y**2) + (r1**2+a**2)*z**2 )
!!$    r2e = Sqrt(r2**2+a**2)*r2*re / Sqrt(r2**2*(x**2+y**2) + (r2**2+a**2)*z**2 )
!!$
!!$    if (r < r2 .AND. re > r1) then
!!$       TrK = (TrK*(re - r1)**6*((r1 - r2e)**5 + &
!!$            &      (re - r2e)*(6*(r1 - r2e)**4 + &
!!$            &         7*(re - r2e)*(36*re**3 + 3*r1**3 + 18*re**2*(r1 - 7*r2e) - &
!!$            &            17*r1**2*r2e + 43*r1*r2e**2 - 65*r2e**3 + &
!!$            &            4*re*(2*r1**2 - 13*r1*r2e + 38*r2e**2)))))/(r1 - r2e)**11
!!$    else if (re < r1) then
!!$       TrK = 0.
!!$    end if

    r1 = 0.4D0 * rplus
    r2 = 0.8D0 * rplus

    r1e = Sqrt(r1**2+ain**2)*r1*re / Sqrt(r1**2*(x**2+y**2) + (r1**2+ain**2)*z**2 )
    r2e = Sqrt(r2**2+ain**2)*r2*re / Sqrt(r2**2*(x**2+y**2) + (r2**2+ain**2)*z**2 )

    interpconst = 0.1
    r3e = r2e - 0.1*M
    if (r < r2) then
       if (re > r3e) then
          TrK=TrK*(interpconst + ((-1 + interpconst)*(r3e - re)**5* &
     &     (126*r2e**4 + r3e**4 + 5*r3e**3*re + 15*r3e**2*re**2 + 35*r3e*re**3 + &
     &       70*re**4 - 84*r2e**3*(r3e + 5*re) + &
     &       36*r2e**2*(r3e**2 + 5*r3e*re + 15*re**2) - &
     &       9*r2e*(r3e**3 + 5*r3e**2*re + 15*r3e*re**2 + 35*re**3)))/&
     &   (r2e - r3e)**9)
       else
          TrK=TrK*interpconst
       end if
    end if


    return
  end subroutine


 end module

