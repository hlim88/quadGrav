!------------------------------------------------------------------------
!
!   $Id: def_initial.f90,v 1.28 2010-03-13 00:05:31 megevand Exp $
!
!------------------------------------------------------------------------
#include "cctk.h"
#include "basic.h"

#if defined GH_MHD || defined GH_FLOWER
  module initial_gr
#else
  module initial
#endif
  use params
#if defined GH_MHD || defined GH_FLOWER
  use gh_exact_quantities
#else
  use exact_quantities
#endif
  use m_gr_id_file

  implicit none
CONTAINS
#if defined GH_MHD || defined GH_FLOWER
  subroutine initialdatagr(u0,u2,v,w,par)
#else
  subroutine initialdata(u0,u2,v,w,imask,par)
#endif
    use GF

    type(gridfunction), dimension(NU_G) :: u0, u2
    type(gridfunction), dimension(NV_G) :: v
    type(gridfunction), dimension(NW) :: w
#if defined GH_MHD || defined GH_FLOWER
    ! imask is now passed by HyperInit
#else
    CCTK_INT, dimension(:,:,:)        :: imask
#endif
    CCTK_REAL, dimension(:) :: par

    !! local vars
    CCTK_REAL, dimension(:,:,:), pointer :: xx, yy, zz
    CCTK_REAL               :: r,r0,delta_r,amp_phi,amp_phi_dot
    CCTK_REAL		    :: factor,factor2,mass
    CCTK_INT                :: i,j,k, idtype, nx, ny, nz
    CCTK_INT                :: istat, rhs_type, ii, jj, kk, iround
    logical, parameter      :: ltrace = .false.
    CCTK_REAL	            :: x,y,z,rad_exc,dx
    logical :: masking
    
    CCTK_REAL                         :: phir_i, phic_i, phim_i, phin_i 
    CCTK_REAL, dimension(0:3)         :: H_i,dphir_i,dphic_i,dphim_i,dphin_i,dH_i
    CCTK_REAL, dimension(0:3,0:3)     :: g_i, K_i
    CCTK_REAL, dimension(1:3,0:3,0:3) :: dg_i
    
    CCTK_REAL, allocatable, dimension(:)  :: x1d, y1d, z1d
    CCTK_REAL :: myl2norm3d
    ! String stuff:
    CCTK_REAL :: hi, tmp
    CCTK_REAL :: xcR, zcR,xcL, zcL, xxL, zzL, ffL, rrL
    CCTK_REAL ::                    xxR, zzR, ffR, rrR
    CCTK_REAL :: chargeL, chargeR, charge, stringr, blendwidth, string_eta
    CCTK_REAL :: chitilde, ztilde, rtilde, ff
    logical   :: ONLEFT
    integer   :: myid, proc_return_myid

    interface
      subroutine sdf_file_mem(filename, shp)
         implicit none
         character*(*)   filename
         CCTK_INT         ::  shp(3)
      end subroutine sdf_file_mem
        subroutine D43NoMask_simple(Du, u, direction, par)
          use params
          implicit none
          CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
          CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
          CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
          CCTK_INT :: direction
        end subroutine D43NoMask_simple

    end interface



    myid = proc_return_myid()
    if (ltrace) write(*,*) myid, "-------------------Entering initial---------------------"

    ! SLL: initialize dh_i since the exact() routines just add to it
    !      and depend on it getting initiailized explicitly
    dH_i=0.d0

    xx => w(H_XPHYS)%d
    yy => w(H_YPHYS)%d
    zz => w(H_ZPHYS)%d

    if(par(P_USE_MASK).eq.1) then
      masking=.true.
    else
      masking =.false.
    end if


    do k=1,NU_G
       u0(k)%d = 0.0
    end do

    idtype = nint(par(P_IDTYPE))

    dx = par(P_DX)
    mass = par(P_BH1_MASS)
    nx = nint(par(P_NX))
    ny = nint(par(P_NY))
    nz = nint(par(P_NZ))

    if (idtype .EQ. -50) then
      print*, "calling Lorene for the initial data"
      call read_lorene_bbh(u0, v, w, par)
!!MM
      print*, "computing the numerical derivatives of the metric"
      call D43NoMask_simple(u0(H_D1G00)%d, u0(H_G00)%d, 1, par)
      call D43NoMask_simple(u0(H_D1G01)%d, u0(H_G01)%d, 1, par)
      call D43NoMask_simple(u0(H_D1G02)%d, u0(H_G02)%d, 1, par)
      call D43NoMask_simple(u0(H_D1G03)%d, u0(H_G03)%d, 1, par)
      call D43NoMask_simple(u0(H_D1G11)%d, u0(H_G11)%d, 1, par)
      call D43NoMask_simple(u0(H_D1G12)%d, u0(H_G12)%d, 1, par)
      call D43NoMask_simple(u0(H_D1G13)%d, u0(H_G13)%d, 1, par)
      call D43NoMask_simple(u0(H_D1G22)%d, u0(H_G22)%d, 1, par)
      call D43NoMask_simple(u0(H_D1G23)%d, u0(H_G23)%d, 1, par)
      call D43NoMask_simple(u0(H_D1G33)%d, u0(H_G33)%d, 1, par)

      call D43NoMask_simple(u0(H_D2G00)%d, u0(H_G00)%d, 2, par)
      call D43NoMask_simple(u0(H_D2G01)%d, u0(H_G01)%d, 2, par)
      call D43NoMask_simple(u0(H_D2G02)%d, u0(H_G02)%d, 2, par)
      call D43NoMask_simple(u0(H_D2G03)%d, u0(H_G03)%d, 2, par)
      call D43NoMask_simple(u0(H_D2G11)%d, u0(H_G11)%d, 2, par)
      call D43NoMask_simple(u0(H_D2G12)%d, u0(H_G12)%d, 2, par)
      call D43NoMask_simple(u0(H_D2G13)%d, u0(H_G13)%d, 2, par)
      call D43NoMask_simple(u0(H_D2G22)%d, u0(H_G22)%d, 2, par)
      call D43NoMask_simple(u0(H_D2G23)%d, u0(H_G23)%d, 2, par)
      call D43NoMask_simple(u0(H_D2G33)%d, u0(H_G33)%d, 2, par)

      call D43NoMask_simple(u0(H_D3G00)%d, u0(H_G00)%d, 3, par)
      call D43NoMask_simple(u0(H_D3G01)%d, u0(H_G01)%d, 3, par)
      call D43NoMask_simple(u0(H_D3G02)%d, u0(H_G02)%d, 3, par)
      call D43NoMask_simple(u0(H_D3G03)%d, u0(H_G03)%d, 3, par)
      call D43NoMask_simple(u0(H_D3G11)%d, u0(H_G11)%d, 3, par)
      call D43NoMask_simple(u0(H_D3G12)%d, u0(H_G12)%d, 3, par)
      call D43NoMask_simple(u0(H_D3G13)%d, u0(H_G13)%d, 3, par)
      call D43NoMask_simple(u0(H_D3G22)%d, u0(H_G22)%d, 3, par)
      call D43NoMask_simple(u0(H_D3G23)%d, u0(H_G23)%d, 3, par)	
      call D43NoMask_simple(u0(H_D3G33)%d, u0(H_G33)%d, 3, par)
!!MM

    end if
    if (idtype .EQ. -60) then
      call read_cook_bbh(u0, v, w, par)
    end if

    if ( (idtype .LE. -1) .AND. (idtype .GE. -49) ) then
      ! read initial data from files
      if (ltrace) then
        write(*,*)'@@@ Read data from files'
      end if      
      call file_data_gr(u0, u2, w, xx, yy, zz, nx, ny, nz, par)
    end if

    ! Setting flat space
    if (idtype .EQ. 0) then
      call gh_flatspace(u0, par)
      return
    else if (idtype .EQ. 200) then
      !
      ! Setting cosmic strings:
      !    Need to edit the following so that m^2 can be negative:
      !       had/src/hyperGHSF/rhs.f90           -- remove "**2" from occurences of "sf_m"
      !       had/src/hyperGHSF/auxvars.f90       -- remove "**2" from occurences of "sf_m"
      !
      !    Setting the potential:
      !       sf_lambda =  1.0
      !       sf_mass   = -2*lambda*string*eta**2
      !       
      !    string_eta   = scale of symmetry breaking
      !       
      !    String Loop:
      !           R0  --- radius of loop
      !       
      !    Parallel Strings: (left and right)
      !           center of string in x-z plane (xcL,zcL)
      !           charge of string: charge = +/- 1
      !           radius of string: stringr
      !
      call gh_flatspace(u0, par)
      hi = w(H_XPHYS)%d(2,1,1)-w(H_XPHYS)%d(1,1,1)
      write(*,*) 'cosmic string, hi=', hi
      blendwidth = 0.5
      stringr    = 0.5
      string_eta = 0.01
      !
      zcL     =  0.0
      xcL     = -1.5
      chargeL = -1.0
      !
      zcR     =  0.0
      xcR     =  1.5
      chargeR = +1.0
      !
      R0      =  10.0
      !
      if (.true.) then
         !
         ! String loop:
         !
         write(*,*) 'cosmic string, loop'
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  chitilde = sqrt(xx(i,j,k)**2+yy(i,j,k)**2) - R0
                  ztilde   = zz(i,j,k)
                  rtilde   = sqrt( chitilde**2 + ztilde**2 )
                  if (rtilde.lt.stringr) then
                     ff = rtilde/stringr * string_eta
                  else
                     ff = string_eta
                  end if
                  !
                  if (rtilde.ne.0) then
                     u0(H_PHIR)%d(i,j,k) = ff*chitilde / rtilde
                     u0(H_PHIC)%d(i,j,k) = ff*  ztilde / rtilde
                  else
                     u0(H_PHIR)%d(i,j,k) = 0.0
                     u0(H_PHIC)%d(i,j,k) = 0.0
                  end if
                  !
                  u0(H_PIR)%d(i,j,k)  = 0.0
                  u0(H_PIC)%d(i,j,k)  = 0.0
                  !
               end do
            end do
         end do
         ! END: String loop
      else
         !
         ! Two parallel straight strings:
         !
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  !
                  ONLEFT = .false.
                  if (xx(i,1,1) .lt. 0) ONLEFT = .true.
                  !
                  xxL    = xx(i,1,1) - xcL
                  zzL    = zz(1,1,k) - zcL
                  rrL    = sqrt( xxL**2 + zzL**2)
                  if (rrL.lt.stringr) then
                     ffL = rrL/stringr * string_eta
                  else
                     ffL = string_eta
                  end if
                  !
                  xxR    = xx(i,1,1) - xcR
                  zzR    = zz(1,1,k) - zcR
                  rrR    = sqrt( xxR**2 + zzR**2)
                  if (rrR.lt.stringr) then
                     ffR = rrR/stringr * string_eta
                  else
                     ffR = string_eta
                  end if
                  !write(*,*) 'cosmic string:',xxx,zzz,rr,ff
                  !
                  if ( abs(xx(i,1,1)) .lt. blendwidth ) then
                     tmp = ( blendwidth - xx(i,1,1) ) / (2.0*blendwidth)
                     u0(H_PHIR)%d(i,j,k) =   tmp      *chargeL*ffL*xxL/rrL &
                                           + (1.0-tmp)*chargeR*ffR*xxR/rrR
                     u0(H_PHIC)%d(i,j,k) =   tmp      *chargeL*ffL*zzL/rrL &
                                           + (1.0-tmp)*chargeR*ffR*zzR/rrR
                     if (rrL.eq.0.or.rrR.eq.0) then
                        u0(H_PHIR)%d(i,j,k) = 0.0
                        u0(H_PHIC)%d(i,j,k) = 0.0
                     end if
                  else if (ONLEFT) then
                     u0(H_PHIR)%d(i,j,k) = chargeL * ffL * xxL / rrL
                     u0(H_PHIC)%d(i,j,k) = chargeL * ffL * zzL / rrL
                     if (rrL.eq.0) then
                        u0(H_PHIR)%d(i,j,k) = 0.0
                        u0(H_PHIC)%d(i,j,k) = 0.0
                     end if
                  else 
                     u0(H_PHIR)%d(i,j,k) = chargeR * ffR * xxR / rrR
                     u0(H_PHIC)%d(i,j,k) = chargeR * ffR * zzR / rrR
                     if (rrR.eq.0) then
                        u0(H_PHIR)%d(i,j,k) = 0.0
                        u0(H_PHIC)%d(i,j,k) = 0.0
                     end if
                  end if
                  !
                  u0(H_PIR)%d(i,j,k)  = 0.0
                  u0(H_PIC)%d(i,j,k)  = 0.0
               end do
            end do
         end do 
         ! END: Parallel strings
      end if
      call deriv2_x(   u0(H_D1PHIR)%d, u0(H_PHIR)%d, hi, nx, ny, nz )
      call deriv2_y(   u0(H_D2PHIR)%d, u0(H_PHIR)%d, hi, nx, ny, nz )
      call deriv2_z(   u0(H_D3PHIR)%d, u0(H_PHIR)%d, hi, nx, ny, nz )
      call deriv2_x(   u0(H_D1PHIC)%d, u0(H_PHIC)%d, hi, nx, ny, nz )
      call deriv2_y(   u0(H_D2PHIC)%d, u0(H_PHIC)%d, hi, nx, ny, nz )
      call deriv2_z(   u0(H_D3PHIC)%d, u0(H_PHIC)%d, hi, nx, ny, nz )
      return
      !
    end if


    iround = 1
           
    do k=1,nz
      do j=1,ny
        do i=1,nx
            
	  x = xx(i,j,k); y=yy(i,j,k); z=zz(i,j,k)

          g_i(0,0) = u0(H_G00)%d(i,j,k)
          g_i(0,1) = u0(H_G01)%d(i,j,k)
          g_i(0,2) = u0(H_G02)%d(i,j,k)
          g_i(0,3) = u0(H_G03)%d(i,j,k)	  	           
          g_i(1,1) = u0(H_G11)%d(i,j,k)
          g_i(1,2) = u0(H_G12)%d(i,j,k)
          g_i(1,3) = u0(H_G13)%d(i,j,k)
          g_i(2,2) = u0(H_G22)%d(i,j,k)
          g_i(2,3) = u0(H_G23)%d(i,j,k)
          g_i(3,3) = u0(H_G33)%d(i,j,k)

	  K_i(0,0) = u0(H_K00)%d(i,j,k)
	  K_i(0,1) = u0(H_K01)%d(i,j,k)
	  K_i(0,2) = u0(H_K02)%d(i,j,k)
	  K_i(0,3) = u0(H_K03)%d(i,j,k)	  	  	  
	  K_i(1,1) = u0(H_K11)%d(i,j,k)
          K_i(1,2) = u0(H_K12)%d(i,j,k)
          K_i(1,3) = u0(H_K13)%d(i,j,k)
          K_i(2,2) = u0(H_K22)%d(i,j,k)
          K_i(2,3) = u0(H_K23)%d(i,j,k)
          K_i(3,3) = u0(H_K33)%d(i,j,k)	  

          dg_i(1,0,0) = u0(H_D1G00)%d(i,j,k)
          dg_i(1,0,1) = u0(H_D1G01)%d(i,j,k)
          dg_i(1,0,2) = u0(H_D1G02)%d(i,j,k)
          dg_i(1,0,3) = u0(H_D1G03)%d(i,j,k)
          dg_i(1,1,1) = u0(H_D1G11)%d(i,j,k)
          dg_i(1,1,2) = u0(H_D1G12)%d(i,j,k)
          dg_i(1,1,3) = u0(H_D1G13)%d(i,j,k)
          dg_i(1,2,2) = u0(H_D1G22)%d(i,j,k)
          dg_i(1,2,3) = u0(H_D1G23)%d(i,j,k)
          dg_i(1,3,3) = u0(H_D1G33)%d(i,j,k)

          dg_i(2,0,0) = u0(H_D2G00)%d(i,j,k)
          dg_i(2,0,1) = u0(H_D2G01)%d(i,j,k)
          dg_i(2,0,2) = u0(H_D2G02)%d(i,j,k)
          dg_i(2,0,3) = u0(H_D2G03)%d(i,j,k)
	  dg_i(2,1,1) = u0(H_D2G11)%d(i,j,k)
          dg_i(2,1,2) = u0(H_D2G12)%d(i,j,k)
          dg_i(2,1,3) = u0(H_D2G13)%d(i,j,k)
          dg_i(2,2,2) = u0(H_D2G22)%d(i,j,k)
          dg_i(2,2,3) = u0(H_D2G23)%d(i,j,k)
          dg_i(2,3,3) = u0(H_D2G33)%d(i,j,k)
	  
          dg_i(3,0,0) = u0(H_D3G00)%d(i,j,k)
          dg_i(3,0,1) = u0(H_D3G01)%d(i,j,k)
          dg_i(3,0,2) = u0(H_D3G02)%d(i,j,k)
          dg_i(3,0,3) = u0(H_D3G03)%d(i,j,k)
          dg_i(3,1,1) = u0(H_D3G11)%d(i,j,k)
          dg_i(3,1,2) = u0(H_D3G12)%d(i,j,k)
          dg_i(3,1,3) = u0(H_D3G13)%d(i,j,k)
          dg_i(3,2,2) = u0(H_D3G22)%d(i,j,k)
          dg_i(3,2,3) = u0(H_D3G23)%d(i,j,k)
          dg_i(3,3,3) = u0(H_D3G33)%d(i,j,k)

	  phim_i     = u0(H_PHIM)%d(i,j,k)
	  dphim_i(0) = u0(H_PIM)%d(i,j,k)
	  dphim_i(1) = u0(H_D1PHIM)%d(i,j,k)
          dphim_i(2) = u0(H_D2PHIM)%d(i,j,k)
          dphim_i(3) = u0(H_D3PHIM)%d(i,j,k)

	  phin_i     = u0(H_PHIN)%d(i,j,k)
	  dphin_i(0) = u0(H_PIN)%d(i,j,k)
	  dphin_i(1) = u0(H_D1PHIN)%d(i,j,k)
          dphin_i(2) = u0(H_D2PHIN)%d(i,j,k)
          dphin_i(3) = u0(H_D3PHIN)%d(i,j,k)

	  phir_i     = u0(H_PHIR)%d(i,j,k)
	  dphir_i(0) = u0(H_PIR)%d(i,j,k)
	  dphir_i(1) = u0(H_D1PHIR)%d(i,j,k)
          dphir_i(2) = u0(H_D2PHIR)%d(i,j,k)
          dphir_i(3) = u0(H_D3PHIR)%d(i,j,k)

	  phic_i     = u0(H_PHIC)%d(i,j,k)
	  dphic_i(0) = u0(H_PIC)%d(i,j,k)
	  dphic_i(1) = u0(H_D1PHIC)%d(i,j,k)
          dphic_i(2) = u0(H_D2PHIC)%d(i,j,k)
          dphic_i(3) = u0(H_D3PHIC)%d(i,j,k)


! if there is a mask do not compute those points, except when there is the 
! addition of two solutions; the black hole is always the last one, but 
! the first non-singular solution has to be defined everywhere for the derivatives
	  if ((masking .EQV. .false.) .OR. &
	 &    ((masking .EQV. .true.) .AND. (w(H_MASK)%d(i,j,k).gt.0)) .OR. &
!!MM
!!	 &    (par(P_idtype).ge.100) ) then
         &    (par(P_idtype).ge.100) .OR. (par(P_idtype).eq.-50) ) then
!!MM
#if defined GH_MHD || defined GH_FLOWER
    	    call exact_gr(iround,par,x,y,z,g_i,K_i,dg_i,H_i, dH_i, &
         &  phir_i,dphir_i,phic_i,dphic_i,phim_i,dphim_i,phin_i,dphin_i)
#else
    	    call exact(iround,par,x,y,z,g_i,K_i,dg_i,H_i, dH_i, &
         &  phir_i,dphir_i,phic_i,dphic_i,phim_i,dphim_i,phin_i,dphin_i)
#endif
	  end if

	  ! the metric components                
          u0(H_G00)%d(i,j,k) = g_i(0,0)
          u0(H_G01)%d(i,j,k) = g_i(0,1)
          u0(H_G02)%d(i,j,k) = g_i(0,2)
          u0(H_G03)%d(i,j,k) = g_i(0,3)
          u0(H_G11)%d(i,j,k) = g_i(1,1)
          u0(H_G12)%d(i,j,k) = g_i(1,2)
          u0(H_G13)%d(i,j,k) = g_i(1,3)
          u0(H_G22)%d(i,j,k) = g_i(2,2)
          u0(H_G23)%d(i,j,k) = g_i(2,3)
          u0(H_G33)%d(i,j,k) = g_i(3,3)

	  u0(H_K00)%d(i,j,k) = K_i(0,0)
	  u0(H_K01)%d(i,j,k) = K_i(0,1)
	  u0(H_K02)%d(i,j,k) = K_i(0,2)
	  u0(H_K03)%d(i,j,k) = K_i(0,3)
	  u0(H_K11)%d(i,j,k) = K_i(1,1)
          u0(H_K12)%d(i,j,k) = K_i(1,2)
          u0(H_K13)%d(i,j,k) = K_i(1,3)
          u0(H_K22)%d(i,j,k) = K_i(2,2)
          u0(H_K23)%d(i,j,k) = K_i(2,3)
          u0(H_K33)%d(i,j,k) = K_i(3,3)
	  
	  
          ! their space derivatives
          u0(H_D1G00)%d(i,j,k) = dg_i(1,0,0)
          u0(H_D1G01)%d(i,j,k) = dg_i(1,0,1)
          u0(H_D1G02)%d(i,j,k) = dg_i(1,0,2)
          u0(H_D1G03)%d(i,j,k) = dg_i(1,0,3)
          u0(H_D1G11)%d(i,j,k) = dg_i(1,1,1)
          u0(H_D1G12)%d(i,j,k) = dg_i(1,1,2)
          u0(H_D1G13)%d(i,j,k) = dg_i(1,1,3)
          u0(H_D1G22)%d(i,j,k) = dg_i(1,2,2)
          u0(H_D1G23)%d(i,j,k) = dg_i(1,2,3)
          u0(H_D1G33)%d(i,j,k) = dg_i(1,3,3)

          u0(H_D2G00)%d(i,j,k) = dg_i(2,0,0)
          u0(H_D2G01)%d(i,j,k) = dg_i(2,0,1)
          u0(H_D2G02)%d(i,j,k) = dg_i(2,0,2)
          u0(H_D2G03)%d(i,j,k) = dg_i(2,0,3)          
	  u0(H_D2G11)%d(i,j,k) = dg_i(2,1,1)
          u0(H_D2G12)%d(i,j,k) = dg_i(2,1,2)
          u0(H_D2G13)%d(i,j,k) = dg_i(2,1,3)
          u0(H_D2G22)%d(i,j,k) = dg_i(2,2,2)
          u0(H_D2G23)%d(i,j,k) = dg_i(2,2,3)
          u0(H_D2G33)%d(i,j,k) = dg_i(2,3,3)

          u0(H_D3G00)%d(i,j,k) = dg_i(3,0,0)
          u0(H_D3G01)%d(i,j,k) = dg_i(3,0,1)
          u0(H_D3G02)%d(i,j,k) = dg_i(3,0,2)
          u0(H_D3G03)%d(i,j,k) = dg_i(3,0,3)
          u0(H_D3G11)%d(i,j,k) = dg_i(3,1,1)
          u0(H_D3G12)%d(i,j,k) = dg_i(3,1,2)
          u0(H_D3G13)%d(i,j,k) = dg_i(3,1,3)
          u0(H_D3G22)%d(i,j,k) = dg_i(3,2,2)
          u0(H_D3G23)%d(i,j,k) = dg_i(3,2,3)
          u0(H_D3G33)%d(i,j,k) = dg_i(3,3,3)

          ! their time derivatives 
          u0(H_H0)%d(i,j,k) = H_i(0)         
	  u0(H_H1)%d(i,j,k) = H_i(1)
          u0(H_H2)%d(i,j,k) = H_i(2)
          u0(H_H3)%d(i,j,k) = H_i(3)

          u0(H_G0)%d(i,j,k)   = dH_i(0)
          u0(H_D1H0)%d(i,j,k) = dH_i(1)
          u0(H_D2H0)%d(i,j,k) = dH_i(2)
          u0(H_D3H0)%d(i,j,k) = dH_i(3)
  
          u0(H_PHIR)%d(i,j,k)   = phir_i
          u0(H_PIR)%d(i,j,k)    = dphir_i(0)
          u0(H_D1PHIR)%d(i,j,k) = dphir_i(1)
          u0(H_D2PHIR)%d(i,j,k) = dphir_i(2)
          u0(H_D3PHIR)%d(i,j,k) = dphir_i(3)

	  u0(H_PHIC)%d(i,j,k)   = phic_i
	  u0(H_PIC)%d(i,j,k)    = dphic_i(0)
	  u0(H_D1PHIC)%d(i,j,k) = dphic_i(1)
          u0(H_D2PHIC)%d(i,j,k) = dphic_i(2)
          u0(H_D3PHIC)%d(i,j,k) = dphic_i(3)

	  u0(H_PHIM)%d(i,j,k)   = phim_i
	  u0(H_PIM)%d(i,j,k)    = dphim_i(0)	    	    
	  u0(H_D1PHIM)%d(i,j,k) = dphim_i(1)
          u0(H_D2PHIM)%d(i,j,k) = dphim_i(2)
          u0(H_D3PHIM)%d(i,j,k) = dphim_i(3)

	  u0(H_PHIN)%d(i,j,k)   = phin_i
	  u0(H_PIN)%d(i,j,k)    = dphin_i(0)	    	    
	  u0(H_D1PHIN)%d(i,j,k) = dphin_i(1)
          u0(H_D2PHIN)%d(i,j,k) = dphin_i(2)
          u0(H_D3PHIN)%d(i,j,k) = dphin_i(3)
	  	    
        end do
      end do
    end do

    if ((par(P_idtype).eq.10) .OR. (par(P_idtype).eq.11) &
    &    .OR.(par(P_idtype).eq.-50) .OR. (par(P_idtype).ge.100) ) then
        call D43NoMask_simple(u0(H_D1G00)%d, u0(H_G00)%d, 1, par)
        call D43NoMask_simple(u0(H_D1G01)%d, u0(H_G01)%d, 1, par)
        call D43NoMask_simple(u0(H_D1G02)%d, u0(H_G02)%d, 1, par)
        call D43NoMask_simple(u0(H_D1G03)%d, u0(H_G03)%d, 1, par)
        call D43NoMask_simple(u0(H_D1G11)%d, u0(H_G11)%d, 1, par)
        call D43NoMask_simple(u0(H_D1G12)%d, u0(H_G12)%d, 1, par)
        call D43NoMask_simple(u0(H_D1G13)%d, u0(H_G13)%d, 1, par)
        call D43NoMask_simple(u0(H_D1G22)%d, u0(H_G22)%d, 1, par)
        call D43NoMask_simple(u0(H_D1G23)%d, u0(H_G23)%d, 1, par)        
	call D43NoMask_simple(u0(H_D1G33)%d, u0(H_G33)%d, 1, par)

        call D43NoMask_simple(u0(H_D2G00)%d, u0(H_G00)%d, 2, par)
        call D43NoMask_simple(u0(H_D2G01)%d, u0(H_G01)%d, 2, par)
        call D43NoMask_simple(u0(H_D2G02)%d, u0(H_G02)%d, 2, par)
        call D43NoMask_simple(u0(H_D2G03)%d, u0(H_G03)%d, 2, par)
        call D43NoMask_simple(u0(H_D2G11)%d, u0(H_G11)%d, 2, par)
        call D43NoMask_simple(u0(H_D2G12)%d, u0(H_G12)%d, 2, par)
        call D43NoMask_simple(u0(H_D2G13)%d, u0(H_G13)%d, 2, par)
        call D43NoMask_simple(u0(H_D2G22)%d, u0(H_G22)%d, 2, par)
        call D43NoMask_simple(u0(H_D2G23)%d, u0(H_G23)%d, 2, par)
        call D43NoMask_simple(u0(H_D2G33)%d, u0(H_G33)%d, 2, par)

        call D43NoMask_simple(u0(H_D3G00)%d, u0(H_G00)%d, 3, par)
        call D43NoMask_simple(u0(H_D3G01)%d, u0(H_G01)%d, 3, par)
        call D43NoMask_simple(u0(H_D3G02)%d, u0(H_G02)%d, 3, par)
        call D43NoMask_simple(u0(H_D3G03)%d, u0(H_G03)%d, 3, par)
        call D43NoMask_simple(u0(H_D3G11)%d, u0(H_G11)%d, 3, par)
        call D43NoMask_simple(u0(H_D3G12)%d, u0(H_G12)%d, 3, par)
        call D43NoMask_simple(u0(H_D3G13)%d, u0(H_G13)%d, 3, par)
        call D43NoMask_simple(u0(H_D3G22)%d, u0(H_G22)%d, 3, par)
        call D43NoMask_simple(u0(H_D3G23)%d, u0(H_G23)%d, 3, par)	
        call D43NoMask_simple(u0(H_D3G33)%d, u0(H_G33)%d, 3, par)

        call D43NoMask_simple(u0(H_D1H0)%d, u0(H_H0)%d, 1, par)
        call D43NoMask_simple(u0(H_D2H0)%d, u0(H_H0)%d, 2, par)
        call D43NoMask_simple(u0(H_D3H0)%d, u0(H_H0)%d, 3, par)

    end if


    if (ltrace) then
      write(*,*) myid, '@@@ L2 norms from files after round 1'
      write(*,*) '@@@ :  ||g00|| = ',myl2norm3d(u0(H_G00)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g01|| = ',myl2norm3d(u0(H_G01)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g02|| = ',myl2norm3d(u0(H_G02)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g03|| = ',myl2norm3d(u0(H_G03)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g11|| = ',myl2norm3d(u0(H_G11)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g12|| = ',myl2norm3d(u0(H_G12)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g13|| = ',myl2norm3d(u0(H_G13)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g22|| = ',myl2norm3d(u0(H_G22)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g23|| = ',myl2norm3d(u0(H_G23)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g33|| = ',myl2norm3d(u0(H_G33)%d, nx, ny, nz)

      write(*,*) '@@@ :  ||d1g00|| = ',myl2norm3d(u0(H_D1G00)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g01|| = ',myl2norm3d(u0(H_D1G01)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g02|| = ',myl2norm3d(u0(H_D1G02)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g03|| = ',myl2norm3d(u0(H_D1G03)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g11|| = ',myl2norm3d(u0(H_D1G11)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g12|| = ',myl2norm3d(u0(H_D1G12)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g13|| = ',myl2norm3d(u0(H_D1G13)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g22|| = ',myl2norm3d(u0(H_D1G22)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g23|| = ',myl2norm3d(u0(H_D1G23)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g33|| = ',myl2norm3d(u0(H_D1G33)%d, nx, ny, nz)
            
      write(*,*) '@@@ :  ||k00|| = ',myl2norm3d(u0(H_K00)%d, nx, ny, nz)      
      write(*,*) '@@@ :  ||k01|| = ',myl2norm3d(u0(H_K01)%d, nx, ny, nz)      
      write(*,*) '@@@ :  ||k02|| = ',myl2norm3d(u0(H_K02)%d, nx, ny, nz)      
      write(*,*) '@@@ :  ||k03|| = ',myl2norm3d(u0(H_K03)%d, nx, ny, nz)                        
      write(*,*) '@@@ :  ||k11|| = ',myl2norm3d(u0(H_K11)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k12|| = ',myl2norm3d(u0(H_K12)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k13|| = ',myl2norm3d(u0(H_K13)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k22|| = ',myl2norm3d(u0(H_K22)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k23|| = ',myl2norm3d(u0(H_K23)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k33|| = ',myl2norm3d(u0(H_K33)%d, nx, ny, nz)

    
      write(*,*) '@@@ :  ||phir||  = ',myl2norm3d(u0(H_PHIR)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||pir||   = ',myl2norm3d(u0(H_PIR)%d, nx, ny, nz)	
      write(*,*) '@@@ :  ||d1phir||= ',myl2norm3d(u0(H_D1PHIR)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d2phir||= ',myl2norm3d(u0(H_D2PHIR)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d3phir||= ',myl2norm3d(u0(H_D3PHIR)%d, nx, ny, nz)
        
      write(*,*) '@@@ :  ||phic||  = ',myl2norm3d(u0(H_PHIC)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||pic||   = ',myl2norm3d(u0(H_PIC)%d, nx, ny, nz)	
      write(*,*) '@@@ :  ||d1phic||= ',myl2norm3d(u0(H_D1PHIC)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d2phic||= ',myl2norm3d(u0(H_D2PHIC)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d3phic||= ',myl2norm3d(u0(H_D3PHIC)%d, nx, ny, nz)
      
      write(*,*) '@@@ :  ||phim||  = ',myl2norm3d(u0(H_PHIM)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||pim||   = ',myl2norm3d(u0(H_PIM)%d, nx, ny, nz)	
      write(*,*) '@@@ :  ||d1phim||= ',myl2norm3d(u0(H_D1PHIM)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d2phim||= ',myl2norm3d(u0(H_D2PHIM)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d3phim||= ',myl2norm3d(u0(H_D3PHIM)%d, nx, ny, nz)
    
    end if
    
    
    
!   special cases where there is addition of two solutions, like Teuk+Kerr        
    if ((par(P_idtype).ge.100)) then
      iround = 2
           
      do k=1,nz
      do j=1,ny
      do i=1,nx
            
	  x = xx(i,j,k); y=yy(i,j,k); z=zz(i,j,k)

          g_i(0,0) = u0(H_G00)%d(i,j,k)
          g_i(0,1) = u0(H_G01)%d(i,j,k)
          g_i(0,2) = u0(H_G02)%d(i,j,k)
          g_i(0,3) = u0(H_G03)%d(i,j,k)	  	           
          g_i(1,1) = u0(H_G11)%d(i,j,k)
          g_i(1,2) = u0(H_G12)%d(i,j,k)
          g_i(1,3) = u0(H_G13)%d(i,j,k)
          g_i(2,2) = u0(H_G22)%d(i,j,k)
          g_i(2,3) = u0(H_G23)%d(i,j,k)
          g_i(3,3) = u0(H_G33)%d(i,j,k)

	  K_i(0,0) = u0(H_K00)%d(i,j,k)
	  K_i(0,1) = u0(H_K01)%d(i,j,k)
	  K_i(0,2) = u0(H_K02)%d(i,j,k)
	  K_i(0,3) = u0(H_K03)%d(i,j,k)	  	  	  
	  K_i(1,1) = u0(H_K11)%d(i,j,k)
          K_i(1,2) = u0(H_K12)%d(i,j,k)
          K_i(1,3) = u0(H_K13)%d(i,j,k)
          K_i(2,2) = u0(H_K22)%d(i,j,k)
          K_i(2,3) = u0(H_K23)%d(i,j,k)
          K_i(3,3) = u0(H_K33)%d(i,j,k)	  

          dg_i(1,0,0) = u0(H_D1G00)%d(i,j,k)
          dg_i(1,0,1) = u0(H_D1G01)%d(i,j,k)
          dg_i(1,0,2) = u0(H_D1G02)%d(i,j,k)
          dg_i(1,0,3) = u0(H_D1G03)%d(i,j,k)
          dg_i(1,1,1) = u0(H_D1G11)%d(i,j,k)
          dg_i(1,1,2) = u0(H_D1G12)%d(i,j,k)
          dg_i(1,1,3) = u0(H_D1G13)%d(i,j,k)
          dg_i(1,2,2) = u0(H_D1G22)%d(i,j,k)
          dg_i(1,2,3) = u0(H_D1G23)%d(i,j,k)
          dg_i(1,3,3) = u0(H_D1G33)%d(i,j,k)

          dg_i(2,0,0) = u0(H_D2G00)%d(i,j,k)
          dg_i(2,0,1) = u0(H_D2G01)%d(i,j,k)
          dg_i(2,0,2) = u0(H_D2G02)%d(i,j,k)
          dg_i(2,0,3) = u0(H_D2G03)%d(i,j,k)          
	  dg_i(2,1,1) = u0(H_D2G11)%d(i,j,k)
          dg_i(2,1,2) = u0(H_D2G12)%d(i,j,k)
          dg_i(2,1,3) = u0(H_D2G13)%d(i,j,k)
          dg_i(2,2,2) = u0(H_D2G22)%d(i,j,k)
          dg_i(2,2,3) = u0(H_D2G23)%d(i,j,k)
          dg_i(2,3,3) = u0(H_D2G33)%d(i,j,k)
	  
          dg_i(3,0,0) = u0(H_D3G00)%d(i,j,k)
          dg_i(3,0,1) = u0(H_D3G01)%d(i,j,k)
          dg_i(3,0,2) = u0(H_D3G02)%d(i,j,k)
          dg_i(3,0,3) = u0(H_D3G03)%d(i,j,k)
          dg_i(3,1,1) = u0(H_D3G11)%d(i,j,k)
          dg_i(3,1,2) = u0(H_D3G12)%d(i,j,k)
          dg_i(3,1,3) = u0(H_D3G13)%d(i,j,k)
          dg_i(3,2,2) = u0(H_D3G22)%d(i,j,k)
          dg_i(3,2,3) = u0(H_D3G23)%d(i,j,k)
          dg_i(3,3,3) = u0(H_D3G33)%d(i,j,k)

	  phim_i     = u0(H_PHIM)%d(i,j,k)
	  dphim_i(0) = u0(H_PIM)%d(i,j,k)
	  dphim_i(1) = u0(H_D1PHIM)%d(i,j,k)
          dphim_i(2) = u0(H_D2PHIM)%d(i,j,k)
          dphim_i(3) = u0(H_D3PHIM)%d(i,j,k)

	  phir_i     = u0(H_PHIR)%d(i,j,k)
	  dphir_i(0) = u0(H_PIR)%d(i,j,k)
	  dphir_i(1) = u0(H_D1PHIR)%d(i,j,k)
          dphir_i(2) = u0(H_D2PHIR)%d(i,j,k)
          dphir_i(3) = u0(H_D3PHIR)%d(i,j,k)

	  phic_i     = u0(H_PHIC)%d(i,j,k)
	  dphic_i(0) = u0(H_PIC)%d(i,j,k)
	  dphic_i(1) = u0(H_D1PHIC)%d(i,j,k)
          dphic_i(2) = u0(H_D2PHIC)%d(i,j,k)
          dphic_i(3) = u0(H_D3PHIC)%d(i,j,k)

	  	  	  	  	      
	  if ((masking .EQV. .false.) .OR. &
	 &    ((masking .EQV. .true.) .AND. (w(H_MASK)%d(i,j,k).gt.0)) ) then
#if defined GH_MHD || defined GH_FLOWER
    	    call exact_gr(iround,par,x,y,z,g_i,K_i,dg_i,H_i, dH_i, &
	  &  phir_i,dphir_i,phic_i,dphic_i,phim_i,dphim_i,phin_i,dphin_i)
#else
    	    call exact(iround,par,x,y,z,g_i,K_i,dg_i,H_i, dH_i, &
	  &  phir_i,dphir_i,phic_i,dphic_i,phim_i,dphim_i,phin_i,dphin_i)
#endif
	  end if

	  ! the metric components                
          u0(H_G00)%d(i,j,k) = g_i(0,0)
          u0(H_G01)%d(i,j,k) = g_i(0,1)
          u0(H_G02)%d(i,j,k) = g_i(0,2)
          u0(H_G03)%d(i,j,k) = g_i(0,3)
          u0(H_G11)%d(i,j,k) = g_i(1,1)
          u0(H_G12)%d(i,j,k) = g_i(1,2)
          u0(H_G13)%d(i,j,k) = g_i(1,3)
          u0(H_G22)%d(i,j,k) = g_i(2,2)
          u0(H_G23)%d(i,j,k) = g_i(2,3)
          u0(H_G33)%d(i,j,k) = g_i(3,3)

	  u0(H_K00)%d(i,j,k) = K_i(0,0)
	  u0(H_K01)%d(i,j,k) = K_i(0,1)
	  u0(H_K02)%d(i,j,k) = K_i(0,2)
	  u0(H_K03)%d(i,j,k) = K_i(0,3)
	  u0(H_K11)%d(i,j,k) = K_i(1,1)
          u0(H_K12)%d(i,j,k) = K_i(1,2)
          u0(H_K13)%d(i,j,k) = K_i(1,3)
          u0(H_K22)%d(i,j,k) = K_i(2,2)
          u0(H_K23)%d(i,j,k) = K_i(2,3)
          u0(H_K33)%d(i,j,k) = K_i(3,3)
	  
	  
          ! their space derivatives
          u0(H_D1G00)%d(i,j,k) = dg_i(1,0,0)
          u0(H_D1G01)%d(i,j,k) = dg_i(1,0,1)
          u0(H_D1G02)%d(i,j,k) = dg_i(1,0,2)
          u0(H_D1G03)%d(i,j,k) = dg_i(1,0,3)
          u0(H_D1G11)%d(i,j,k) = dg_i(1,1,1)
          u0(H_D1G12)%d(i,j,k) = dg_i(1,1,2)
          u0(H_D1G13)%d(i,j,k) = dg_i(1,1,3)
          u0(H_D1G22)%d(i,j,k) = dg_i(1,2,2)
          u0(H_D1G23)%d(i,j,k) = dg_i(1,2,3)
          u0(H_D1G33)%d(i,j,k) = dg_i(1,3,3)

          u0(H_D2G00)%d(i,j,k) = dg_i(2,0,0)
          u0(H_D2G01)%d(i,j,k) = dg_i(2,0,1)
          u0(H_D2G02)%d(i,j,k) = dg_i(2,0,2)
          u0(H_D2G03)%d(i,j,k) = dg_i(2,0,3)          
	  u0(H_D2G11)%d(i,j,k) = dg_i(2,1,1)
          u0(H_D2G12)%d(i,j,k) = dg_i(2,1,2)
          u0(H_D2G13)%d(i,j,k) = dg_i(2,1,3)
          u0(H_D2G22)%d(i,j,k) = dg_i(2,2,2)
          u0(H_D2G23)%d(i,j,k) = dg_i(2,2,3)
          u0(H_D2G33)%d(i,j,k) = dg_i(2,3,3)

          u0(H_D3G00)%d(i,j,k) = dg_i(3,0,0)
          u0(H_D3G01)%d(i,j,k) = dg_i(3,0,1)
          u0(H_D3G02)%d(i,j,k) = dg_i(3,0,2)
          u0(H_D3G03)%d(i,j,k) = dg_i(3,0,3)
          u0(H_D3G11)%d(i,j,k) = dg_i(3,1,1)
          u0(H_D3G12)%d(i,j,k) = dg_i(3,1,2)
          u0(H_D3G13)%d(i,j,k) = dg_i(3,1,3)
          u0(H_D3G22)%d(i,j,k) = dg_i(3,2,2)
          u0(H_D3G23)%d(i,j,k) = dg_i(3,2,3)
          u0(H_D3G33)%d(i,j,k) = dg_i(3,3,3)

          ! their time derivatives 
          u0(H_H0)%d(i,j,k) = H_i(0)         
          u0(H_H1)%d(i,j,k) = H_i(1)
          u0(H_H2)%d(i,j,k) = H_i(2)
          u0(H_H3)%d(i,j,k) = H_i(3)

	  u0(H_G0)%d(i,j,k)   = dH_i(0)
	  u0(H_D1H0)%d(i,j,k) = dH_i(1)
          u0(H_D2H0)%d(i,j,k) = dH_i(2)
          u0(H_D3H0)%d(i,j,k) = dH_i(3)
	  
	  u0(H_PHIR)%d(i,j,k)   = phir_i
	  u0(H_PIR)%d(i,j,k)    = dphir_i(0)
	  u0(H_D1PHIR)%d(i,j,k) = dphir_i(1)
          u0(H_D2PHIR)%d(i,j,k) = dphir_i(2)
          u0(H_D3PHIR)%d(i,j,k) = dphir_i(3)

	  u0(H_PHIC)%d(i,j,k)   = phic_i
	  u0(H_PIC)%d(i,j,k)    = dphic_i(0)
	  u0(H_D1PHIC)%d(i,j,k) = dphic_i(1)
          u0(H_D2PHIC)%d(i,j,k) = dphic_i(2)
          u0(H_D3PHIC)%d(i,j,k) = dphic_i(3)

	  u0(H_PHIM)%d(i,j,k)   = phim_i
	  u0(H_PIM)%d(i,j,k)    = dphim_i(0)
	  u0(H_D1PHIM)%d(i,j,k) = dphim_i(1)
          u0(H_D2PHIM)%d(i,j,k) = dphim_i(2)
          u0(H_D3PHIM)%d(i,j,k) = dphim_i(3)
	  	    
      end do
      end do
      end do
    
    end if  
    
	            
    if (ltrace) then
      write(*,*) '@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'    
      write(*,*) myid, '@@@ L2 norms from files after round 2'
      write(*,*) '@@@ :  ||g00|| = ',myl2norm3d(u0(H_G00)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g01|| = ',myl2norm3d(u0(H_G01)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g02|| = ',myl2norm3d(u0(H_G02)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g03|| = ',myl2norm3d(u0(H_G03)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g11|| = ',myl2norm3d(u0(H_G11)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g12|| = ',myl2norm3d(u0(H_G12)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g13|| = ',myl2norm3d(u0(H_G13)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g22|| = ',myl2norm3d(u0(H_G22)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g23|| = ',myl2norm3d(u0(H_G23)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||g33|| = ',myl2norm3d(u0(H_G33)%d, nx, ny, nz)

      write(*,*) '@@@ :  ||d1g00|| = ',myl2norm3d(u0(H_D1G00)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g01|| = ',myl2norm3d(u0(H_D1G01)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g02|| = ',myl2norm3d(u0(H_D1G02)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g03|| = ',myl2norm3d(u0(H_D1G03)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g11|| = ',myl2norm3d(u0(H_D1G11)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g12|| = ',myl2norm3d(u0(H_D1G12)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g13|| = ',myl2norm3d(u0(H_D1G13)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g22|| = ',myl2norm3d(u0(H_D1G22)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g23|| = ',myl2norm3d(u0(H_D1G23)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d1g33|| = ',myl2norm3d(u0(H_D1G33)%d, nx, ny, nz)
            
      write(*,*) '@@@ :  ||k00|| = ',myl2norm3d(u0(H_K00)%d, nx, ny, nz)      
      write(*,*) '@@@ :  ||k01|| = ',myl2norm3d(u0(H_K01)%d, nx, ny, nz)      
      write(*,*) '@@@ :  ||k02|| = ',myl2norm3d(u0(H_K02)%d, nx, ny, nz)      
      write(*,*) '@@@ :  ||k03|| = ',myl2norm3d(u0(H_K03)%d, nx, ny, nz)                        
      write(*,*) '@@@ :  ||k11|| = ',myl2norm3d(u0(H_K11)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k12|| = ',myl2norm3d(u0(H_K12)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k13|| = ',myl2norm3d(u0(H_K13)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k22|| = ',myl2norm3d(u0(H_K22)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k23|| = ',myl2norm3d(u0(H_K23)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||k33|| = ',myl2norm3d(u0(H_K33)%d, nx, ny, nz)

      write(*,*) '@@@ :  ||phir||  = ',myl2norm3d(u0(H_PHIR)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||pir||   = ',myl2norm3d(u0(H_PIR)%d, nx, ny, nz)	
      write(*,*) '@@@ :  ||d1phir||= ',myl2norm3d(u0(H_D1PHIR)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d2phir||= ',myl2norm3d(u0(H_D2PHIR)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d3phir||= ',myl2norm3d(u0(H_D3PHIR)%d, nx, ny, nz)
        
      write(*,*) '@@@ :  ||phic||  = ',myl2norm3d(u0(H_PHIC)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||pic||   = ',myl2norm3d(u0(H_PIC)%d, nx, ny, nz)	
      write(*,*) '@@@ :  ||d1phic||= ',myl2norm3d(u0(H_D1PHIC)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d2phic||= ',myl2norm3d(u0(H_D2PHIC)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d3phic||= ',myl2norm3d(u0(H_D3PHIC)%d, nx, ny, nz)
      
      write(*,*) '@@@ :  ||phim||  = ',myl2norm3d(u0(H_PHIM)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||pim||   = ',myl2norm3d(u0(H_PIM)%d, nx, ny, nz)	
      write(*,*) '@@@ :  ||d1phim||= ',myl2norm3d(u0(H_D1PHIM)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d2phim||= ',myl2norm3d(u0(H_D2PHIM)%d, nx, ny, nz)
      write(*,*) '@@@ :  ||d3phim||= ',myl2norm3d(u0(H_D3PHIM)%d, nx, ny, nz)
    
    end if
    
    if (ltrace) write(*,*) myid, "--------------------- exiting initial -----------------------"
    
  end subroutine

  !---------------------------------------------------------------------
  !
  !
  !
  !---------------------------------------------------------------------
  subroutine read_lorene_bbh(u0, v, w, par)
    implicit none

    type(gridfunction), dimension(NU_G)               :: u0
    type(gridfunction), dimension(NV_G)               :: v
    type(gridfunction), dimension(NW)                 :: w
    CCTK_REAL, dimension(NPAR)                        :: par

    ! local vars
    CCTK_INT                             :: nx, ny, nz
    CCTK_REAL, dimension(:,:,:), pointer :: g11, g12, g13, g22, g23, g33
    CCTK_REAL, dimension(:,:,:), pointer :: K11, K12, K13, K22, K23, K33
    CCTK_REAL, dimension(:,:,:), pointer :: lapse, shift1, shift2, shift3
    CCTK_REAL, dimension(:,:,:), pointer :: x, y, z

!!MM: We'll also need these:
    CCTK_INT         :: i, j, k
    CCTK_REAL        :: alp, b1, b2, b3
    CCTK_REAL        :: huu11,huu22,huu33,huu12,huu13,huu23,deth
!!MM

    x    => w(G_XPHYS)%d
    y    => w(G_YPHYS)%d
    z    => w(G_ZPHYS)%d

    ! The gauge variables will be store in the g0* components
    lapse  => u0(G_G00)%d
    shift1 => u0(G_G01)%d
    shift2 => u0(G_G02)%d
    shift3 => u0(G_G03)%d

    g11 => u0(G_G11)%d
    g12 => u0(G_G12)%d
    g13 => u0(G_G13)%d
    g22 => u0(G_G22)%d
    g23 => u0(G_G23)%d
    g33 => u0(G_G33)%d

    K11 => u0(G_K11)%d
    K12 => u0(G_K12)%d
    K13 => u0(G_K13)%d
    K22 => u0(G_K22)%d
    K23 => u0(G_K23)%d
    K33 => u0(G_K33)%d

    nx = nint(par(P_NX))
    ny = nint(par(P_NY))
    nz = nint(par(P_NZ))

#ifdef LORENE 
    call cook_bbh(   g11, g12, g13,    &
         &           g22, g23, g33,    &
         &           K11, K12, K13, &
         &           K22, K23, K33, &
         &           lapse, shift1, shift2, shift3, &
         &           x, y, z, nx, ny, nz)
#else
    write(0,*)'ERROR: HAD was compiled without the ability to read LORENE data'
    write(0,*)'       specify the location of the LORENE libraries and '
    write(0,*)'       recompile'
    call my_exit('Lorene not available for intial data')
#endif

!!MM: Copying this from hyperGHEM. I'll comment out the computation
!!    of g0a done in subroutine exact
    ! we read the lapse and shift in g0a, now we will compute
    ! the real g0a and save them in the right place
    print*, "computing the metric components"
    do i=1,nx
    do j=1,ny
    do k=1,nz
      shift1(i,j,k) = 0.0d0
      shift2(i,j,k) = 0.0d0
      shift3(i,j,k) = 0.0d0


!      xx = x(i,j,k) - par(P_BH1_X0)
!      yy = y(i,j,k) - par(P_BH1_Y0)
!      zz = z(i,j,k) - par(P_BH1_Z0)
!      rr = dsqrt(xx**2 + yy**2 + zz**2)
!      alpA   = 1.0d0/dsqrt(1.0d0 + 2.0d0*par(P_BH1_MASS)/rr)

!      if (rr .LT. 2) then
!         lapse(i,j,k)= 0.5*(lapse(i,j,k) + alpA)
!      end if
      huu11 = -g23(i,j,k)**2 + g22(i,j,k)*g33(i,j,k)
      huu12 = g13(i,j,k)*g23(i,j,k) - g12(i,j,k)*g33(i,j,k)
      huu13 = -(g13(i,j,k)*g22(i,j,k)) + g12(i,j,k)*g23(i,j,k)
      huu22 = -g13(i,j,k)**2 + g11(i,j,k)*g33(i,j,k)
      huu23 = g12(i,j,k)*g13(i,j,k) - g11(i,j,k)*g23(i,j,k)
      huu33 = -g12(i,j,k)**2 + g11(i,j,k)*g22(i,j,k)
      deth = g11(i,j,k)*huu11 + g12(i,j,k)*huu12 + g13(i,j,k)*huu13

      lapse(i,j,k)= 1.0d0/(deth)**(0.05)
      if (lapse(i,j,k) .LT. 0.2) lapse(i,j,k)=0.2

      alp = lapse(i,j,k)
      b1  = shift1(i,j,k)
      b2  = shift2(i,j,k)
      b3  = shift3(i,j,k)

      lapse(i,j,k) = -alp*alp +  b1*b1*g11(i,j,k) + b2*b2*g22(i,j,k) &
   &               + b3*b3*g33(i,j,k) + 2.0*(b1*(b2*g12(i,j,k) + b3*g13(i,j,k)) &
   &               + b2*b3*g23(i,j,k))
      shift1(i,j,k) = g11(i,j,k)*b1 + g12(i,j,k)*b2 + g13(i,j,k)*b3
      shift2(i,j,k) = g12(i,j,k)*b1 + g22(i,j,k)*b2 + g23(i,j,k)*b3
      shift3(i,j,k) = g13(i,j,k)*b1 + g23(i,j,k)*b2 + g33(i,j,k)*b3

!      lapse(i,j,k) = -1.0d0

    end do
    end do
    end do
!!MM

    return
  end subroutine

  !---------------------------------------------------------------------
  !
  !
  !
  !---------------------------------------------------------------------
  subroutine read_cook_bbh(u0, v, w, par)
    implicit none

    type(gridfunction), dimension(NU_G)               :: u0
    type(gridfunction), dimension(NV_G)               :: v
    type(gridfunction), dimension(NW)                 :: w
    CCTK_REAL, dimension(NPAR)                        :: par

    ! local vars
    CCTK_INT                             :: nx, ny, nz
    CCTK_REAL, dimension(:,:,:), pointer :: g11, g12, g13, g22, g23, g33
    CCTK_REAL, dimension(:,:,:), pointer :: K11, K12, K13, K22, K23, K33
    CCTK_REAL, dimension(:,:,:), pointer :: lapse, shift1, shift2, shift3
    CCTK_REAL, dimension(:,:,:), pointer :: x, y, z

    x    => w(G_XPHYS)%d
    y    => w(G_YPHYS)%d
    z    => w(G_ZPHYS)%d

    ! The gauge variables will be store in the g0* components
    lapse  => u0(G_G00)%d
    shift1 => u0(G_G01)%d
    shift2 => u0(G_G02)%d
    shift3 => u0(G_G03)%d

    g11 => u0(G_G11)%d
    g12 => u0(G_G12)%d
    g13 => u0(G_G13)%d
    g22 => u0(G_G22)%d
    g23 => u0(G_G23)%d
    g33 => u0(G_G33)%d

    K11 => u0(G_K11)%d
    K12 => u0(G_K12)%d
    K13 => u0(G_K13)%d
    K22 => u0(G_K22)%d
    K23 => u0(G_K23)%d
    K33 => u0(G_K33)%d

    nx = nint(par(P_NX))
    ny = nint(par(P_NY))
    nz = nint(par(P_NZ))

#if defined COOK
    call cook_bbh(   g11, g12, g13,    &
         &           g22, g23, g33,    &
         &           K11, K12, K13, &
         &           K22, K23, K33, &
         &           lapse, shift1, shift2, shift3, &
         &           x, y, z, nx, ny, nz)
#else
    write(*,*) ' You must define a c++ compiler in order to use this routine '
#endif

    return
  end subroutine


  !---------------------------------------------------------------------
    
 end module
