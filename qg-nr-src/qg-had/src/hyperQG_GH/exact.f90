#include "cctk.h"

!$Id: exact.f90,v 1.74 2011-09-29 01:52:08 carlos Exp $

#if defined GH_MHD || defined GH_FLOWER
  module gh_exact_quantities
#else
  module exact_quantities
#endif

use params
implicit none

contains

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !
  ! This subroutine computes the exact quantities, needed for initial data, boundary conditions, etc.
  ! As in the whole hyperGR, in a derivative the first index(es) denote 
  ! direction of differentiation, e.g. diff_beta(i,j) = partial_i beta_j
  ! Also, diff is used for derivatives of quantities that are never evolved, 
  ! and D for finite differencing (i.e. things that do evolve). Here everything
  ! is analytic, so it would nt make any difference, but the same notation is 
  ! used for consistency.
  !
  ! IMPORTANT: dd_a is the analytica value of d, while D_a is the derivative of a
  !
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#if defined GH_MHD || defined GH_FLOWER
  subroutine exact_gr(round, par, xold, yold, zold, g_e, K_e, dg_e, &
          &  H_e, dH_e, phir_e, dphir_e, phic_e, dphic_e, &
          &  phim_e, dphim_e, phin_e, dphin_e)
#else
  subroutine exact(round, par, xold, yold, zold, g_e, K_e, dg_e, &
          &  H_e, dH_e, phir_e, dphir_e, phic_e, dphic_e, &
          &  phim_e, dphim_e, phin_e, dphin_e)
#endif

!    use exact_kerrschild
    implicit none

    CCTK_REAL, dimension(:), intent(in) :: par  
    CCTK_REAL, intent(in) :: xold,yold,zold
    
    CCTK_REAL                         :: phir_e,phic_e,phim_e,phin_e 
    CCTK_REAL, dimension(0:3)         :: H_e,dphir_e,dphic_e,dphim_e,dphin_e,dH_e
    CCTK_REAL, dimension(0:3,0:3)     :: g_e, K_e, X_e
    CCTK_REAL, dimension(1:3,0:3,0:3) :: dg_e

    ! local vars
    CCTK_INT :: i,j,kk, li,lj,lk, lm,idtype, nx, matching, initial_H, round
    CCTK_REAL :: factor, t, pi, dx, Lx, epsilon, eps
    CCTK_REAL, dimension(3) :: xx
    CCTK_REAL, dimension(3,3) :: delta
    !for the mr
    CCTK_REAL ::  mr_amp, mr_shift, vx, vy, v2, cgamma
    CCTK_REAL, DIMENSION(100) :: harvest
    !for the twisted
    CCTK_REAL :: a,b,c,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15
    CCTK_REAL :: xc, yc, zc, xo, ro, nfac, psi, drpsi, lap, omega, width
    CCTK_REAL :: x, y, z, r, theta, phi, M, t0
    CCTK_REAL :: drdx, drdy, drdz, amp, lambda
    CCTK_REAL :: dthetadx, dthetady, dthetadz, dphidx, dphidy, dphidz
    CCTK_REAL :: psi1,psi2,psiT,drpsi1,drpsi2,drpsiT
    CCTK_REAL :: frp, ftt1, ftt2, ftp, fpp1, fpp2
    CCTK_REAL :: dAdt, dBdt, dCdt
        
    CCTK_REAL :: alp,alp2,deth,b1,b2,b3,alp_b,alpA,alpB
    CCTK_REAL :: h11,h12,h13,h22,h23,h33            
    CCTK_REAL :: huu11,huu12,huu13,huu22,huu23,huu33    
    CCTK_REAL :: guu00,guu01,guu02,guu03,guu11,guu12,guu13,guu22,guu23,guu33        
    CCTK_REAL :: g00,g01,g02,g03,g11,g12,g13,g22,g23,g33    
    CCTK_REAL :: K00,K01,K02,K03,K11,K12,K13,K22,K23,K33        
    CCTK_REAL :: D000,D001,D002,D003,D011,D012,D013,D022,D023,D033    
    CCTK_REAL :: D100,D101,D102,D103,D111,D112,D113,D122,D123,D133    
    CCTK_REAL :: D200,D201,D202,D203,D211,D212,D213,D222,D223,D233    
    CCTK_REAL :: D300,D301,D302,D303,D311,D312,D313,D322,D323,D333
    CCTK_REAL :: G000,G001,G002,G003,G011,G012,G013,G022,G023,G033    
    CCTK_REAL :: G100,G101,G102,G103,G111,G112,G113,G122,G123,G133    
    CCTK_REAL :: G200,G201,G202,G203,G211,G212,G213,G222,G223,G233    
    CCTK_REAL :: G300,G301,G302,G303,G311,G312,G313,G322,G323,G333                
    CCTK_REAL :: eK11,eK12,eK13,eK22,eK23,eK33        
    CCTK_REAL :: AKud11,AKud12,AKud13,AKud21,AKud22,AKud23,AKud31,AKud32,AKud33
    CCTK_REAL :: BKud11,BKud12,BKud13,BKud21,BKud22,BKud23,BKud31,BKud32,BKud33    
    CCTK_REAL :: eKud11,eKud12,eKud13,eKud21,eKud22,eKud23,eKud31,eKud32,eKud33        
    CCTK_REAL :: dxalpA,dyalpA,dzalpA,dxalpB,dyalpB,dzalpB,dxalp,dyalp,dzalp
    CCTK_REAL :: trG0,trG1,trG2,trG3,Z0,Z1,Z2,Z3
    CCTK_REAL :: H0,H1,H2,H3,d0H0,d1H0,d2H0,d3H0
    CCTK_REAL :: d1alp,d2alp,d3alp
    
    CCTK_REAL :: r2, dtr, dxr, dyr, dzr, H, nult, nulx, nuly, nulz
    CCTK_REAL :: dtH, dtnult, dtnulx, dtnuly, dtnulz
    CCTK_REAL :: dxH, dxnult, dxnulx, dxnuly, dxnulz
    CCTK_REAL :: dyH, dynult, dynulx, dynuly, dynulz
    CCTK_REAL :: dzH, dznult, dznulx, dznuly, dznulz
                        

!added by luis for the Teukolsky data II
    CCTK_REAL :: gout(4,4), gD(4,4,4), gC(4,4), DgC(4,4,4), Ja(4,4), Dja(4,4,4)
    CCTK_REAL :: A0,B0,C0,dr_A,dr_B,dr_C,dt_A,dt_B,dt_C,lam,p
    CCTK_REAL :: dx_r,dy_r,dz_r,dx_theta,dy_theta,dz_theta,dx_phi,dy_phi,dz_phi
    CCTK_REAL :: dxx_r,dxy_r,dxz_r,dxx_theta,dxy_theta,dxz_theta,dxx_phi,dxy_phi,dxz_phi
    CCTK_REAL :: dyx_r,dyy_r,dyz_r,dyx_theta,dyy_theta,dyz_theta,dyx_phi,dyy_phi,dyz_phi
    CCTK_REAL :: dzx_r,dzy_r,dzz_r,dzx_theta,dzy_theta,dzz_theta,dzx_phi,dzy_phi,dzz_phi
    CCTK_REAL :: costeta,cosphi,sinteta,sinphi, Temp1, Temp2, outgoing
    CCTK_REAL :: frr, frt, frf, f1tt, f2tt, ftf, f1ff, f2ff, del
    CCTK_INT  :: ll,lp,mult

    !added for the boosted BS
    CCTK_REAL :: gold(0:3,0:3), Dgold(0:3,0:3,0:3), gnew(0:3,0:3), Dgnew(0:3,0:3,0:3)
    CCTK_REAL :: phirold(0:3), phicold(0:3), Dphirold(0:3), Dphicold(0:3)    
    CCTK_REAL :: phirnew(0:3), phicnew(0:3), Dphirnew(0:3), Dphicnew(0:3)        
    CCTK_REAL :: Jud(0:3,0:3), DJud(0:3,0:3,0:3)
    CCTK_REAL :: rad, sigma, dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, cvx, cvy
    
    logical, parameter :: ltrace = .false.

    !________________________________________________________________________________!

    idtype    = nint(par(P_IDTYPE))
    initial_H = nint(par(P_INITIAL_H))
    
    if (ltrace) then
       write(*,*)'@@@ beginning exact.  idtype = ',idtype
    end if

    pi = 3.1415926535897932d0
    t = par(P_TIME)
    mr_amp = par(P_MR_AMP)    

    x = xold; y = yold; z = zold    
    r = sqrt(x**2+y**2+z**2)

    ! Initialize everything to Minkowski if the data is not given    
    if ((idtype .ge. 0) .AND. (round .eq. 1)) then
      g_e(0,0) =-1.0
      g_e(0,1) = 0.0 
      g_e(0,2) = 0.0
      g_e(0,3) = 0.0
      g_e(1,1) = 1.0
      g_e(1,2) = 0.0
      g_e(1,3) = 0.0
      g_e(2,2) = 1.0 
      g_e(2,3) = 0.0
      g_e(3,3) = 1.0

      K_e  = 0.0
      dg_e = 0.0
      H_e  = 0.0
      dH_e = 0.0

      phir_e  = 0.0
      dphir_e = 0.0
      phic_e  = 0.0
      dphic_e = 0.0
      phim_e  = 0.0
      dphim_e = 0.0
      phin_e  = 0.0
      dphin_e = 0.0

    end if

    !############################################################################
    !###### INITIAL DATA FOR LBNS + BOSON from idfile ###########################
    !############################################################################
    if (idtype.eq.-5) then

#include "GH_initialK.inc"


    !############################################################################
    !###### INITIAL DATA FOR 1-2-3 BOOSTED BOSON STAR (1D) from idfile ##########
    !############################################################################
    else if ((idtype.eq.-21) .OR. (idtype.eq.-22) .OR. (idtype.eq.-23) .OR. (idtype.eq.-24) .OR. &
   &    (idtype.eq.-6) .OR. (idtype.eq.-7) .OR. (idtype.eq.-8) .OR. (idtype.eq.-9) .OR. &
   &    (idtype.eq.-11) .OR. (idtype.eq.-12) .OR. (idtype.eq.-31) .OR. (idtype.eq.-32) ) then

     !--- computing the extrinsic curvature from the time derivative of ---
     !--- the metric that was saved in K_ab and the other quantities ------
     !---------------------------------------------------------------------

#include "GH_initialK.inc"

     !-----------------------------------------------------------------------
     !---computing the pi with the full definition---------------------------
     !-----------------------------------------------------------------------
      dphir_e(0) = -( dphir_e(0) - (b1*dphir_e(1) + b2*dphir_e(2) + b3*dphir_e(3)) )/alp
      dphic_e(0) = -( dphic_e(0) - (b1*dphic_e(1) + b2*dphic_e(2) + b3*dphic_e(3)) )/alp

!      dphim_e(0) = -( dphim_e(0) - (b1*dphim_e(1) + b2*dphim_e(2) + b3*dphim_e(3)) )/alp
!      dphin_e(0) = -( dphin_e(0) - (b1*dphin_e(1) + b2*dphin_e(2) + b3*dphin_e(3)) )/alp


    !########################################################################
    !###### INITIAL DATA FOR 1(+1) BOSON STAR (1D) + 1 BH ISO ###############
    !########################################################################
    else if (idtype.eq.-10) then

     !-------------------------------------------------------------------
     ! compute the spacetime of the isotropic black hole-----------------
     !-------------------------------------------------------------------
      M        = par(P_BH1_MASS)
      matching = par(P_BH1_MATCHING)
      ro       = M/2.0d0

      x = xold - par(P_BH1_X0)
      y = yold - par(P_BH1_Y0)
      z = zold - par(P_BH1_Z0)
      r = sqrt( x**2 + y**2 + z**2 )
      if (r.lt.1E-3) r=1E-3
      
      psi1 = (1.0d0 + M/(2.0*r))
      drpsi1 = - M/(2.0*r**2)      
          
     !--Dressing (2 derivatives continuous)---------------------------------
      if ( (matching .EQ. 1) .AND. (r .LT. ro) )  then
	psi1 = (1.0 + M/ro) - (M/ro**3)*(r**2) + (M/(2*ro**4))*(r**3) 
	drpsi1 = -(2.0*M/(ro**3))*r + ((3.0*M)/(2.0*ro**4))*(r**2)
      end if

      
     !-------------------------------------------------------------------           
     !----the addition of the black hole with the previous BS------------
     !-------------------------------------------------------------------      
          
      g_e(1,1) = g_e(1,1) + psi1**4 - 1.0d0
      g_e(2,2) = g_e(1,1) + psi1**4 - 1.0d0
      g_e(3,3) = g_e(1,1) + psi1**4 - 1.0d0

      dg_e(1,1,1) = dg_e(1,1,1) + 4.0*(psi1**3)*(x/r)*drpsi1
      dg_e(2,1,1) = dg_e(2,1,1) + 4.0*(psi1**3)*(y/r)*drpsi1
      dg_e(3,1,1) = dg_e(3,1,1) + 4.0*(psi1**3)*(z/r)*drpsi1

      dg_e(1,2,2) = dg_e(1,2,2) + 4.0*(psi1**3)*(x/r)*drpsi1
      dg_e(1,3,3) = dg_e(1,3,3) + 4.0*(psi1**3)*(x/r)*drpsi1
      
      dg_e(2,2,2) = dg_e(2,2,2) + 4.0*(psi1**3)*(y/r)*drpsi1
      dg_e(2,3,3) = dg_e(2,3,3) + 4.0*(psi1**3)*(y/r)*drpsi1
      
      dg_e(3,2,2) = dg_e(3,2,2) + 4.0*(psi1**3)*(z/r)*drpsi1
      dg_e(3,3,3) = dg_e(3,3,3) + 4.0*(psi1**3)*(z/r)*drpsi1


!#############################################################################                    
!########## INITIAL DATA FOR 1 BOSON STAR (1D) + 1 KERR ######################
!#############################################################################  

    elseif (idtype.eq.-11) then 
    
     !-------------------------------------------------------------------      
     ! compute the spacetime of the kerr black hole----------------------
     !-------------------------------------------------------------------      
      M = par(P_BH1_MASS)            
      a = par(P_BH1_SPIN)                  

      x = xold - par(P_BH1_X0)
      y = yold - par(P_BH1_Y0)
      z = zold - par(P_BH1_Z0)      
      
      r = sqrt(x**2 + y**2 + z**2)
      if (r.lt.1E-4) r=1E-4
                              
#include "spinKS.inc"     
 
#include "KS_GH.inc"             
          
     !-------------------------------------------------------------------           
     !----the addition of the black hole with the previous BS------------
     !-------------------------------------------------------------------      
     
     !! -- as a first approximation, we will assume than the lapse and shift-
     !! ---are given by the black hole --------------------------------------

!      g_e(0,0) = g_e(0,0) + g00 + 1.0d0
!      g_e(0,1) = g01 + g_e(0,1)
!      g_e(0,2) = g02 + g_e(0,2)
!      g_e(0,3) = g03 + g_e(0,3)                                           

      g_e(0,0) = g00
      g_e(0,1) = g01
      g_e(0,2) = g02
      g_e(0,3) = g03            
      g_e(1,1) = g_e(1,1) + g11 - 1.0d0
      g_e(1,2) = g_e(1,2) + g12
      g_e(1,3) = g_e(1,3) + g13
      g_e(2,2) = g_e(2,2) + g22 - 1.0d0
      g_e(2,3) = g_e(2,3) + g23
      g_e(3,3) = g_e(3,3) + g33 - 1.0d0

!      dg_e(1,0,0) = dg_e(1,0,0) + D100      
!      dg_e(1,0,1) = dg_e(1,0,1) + D101
!      dg_e(1,0,2) = dg_e(1,0,2) + D102
!      dg_e(1,0,3) = dg_e(1,0,3) + D103                  
      
      dg_e(1,0,0) = D100      
      dg_e(1,0,1) = D101
      dg_e(1,0,2) = D102
      dg_e(1,0,3) = D103                        
      dg_e(1,1,1) = dg_e(1,1,1) + D111
      dg_e(1,1,2) = dg_e(1,1,2) + D112
      dg_e(1,1,3) = dg_e(1,1,3) + D113
      dg_e(1,2,2) = dg_e(1,2,2) + D122
      dg_e(1,2,3) = dg_e(1,2,3) + D123
      dg_e(1,3,3) = dg_e(1,3,3) + D133

      
      
!      dg_e(2,0,0) = dg_e(2,0,0) + D200      
!      dg_e(2,0,1) = dg_e(2,0,1) + D201
!      dg_e(2,0,2) = dg_e(2,0,2) + D202
!      dg_e(2,0,3) = dg_e(2,0,3) + D203                                    
      
      dg_e(2,0,0) = D200      
      dg_e(2,0,1) = D201
      dg_e(2,0,2) = D202
      dg_e(2,0,3) = D203                              
      dg_e(2,1,1) = dg_e(2,1,1) + D211
      dg_e(2,1,2) = dg_e(2,1,2) + D212
      dg_e(2,1,3) = dg_e(2,1,3) + D213
      dg_e(2,2,2) = dg_e(2,2,2) + D222
      dg_e(2,2,3) = dg_e(2,2,3) + D223
      dg_e(2,3,3) = dg_e(2,3,3) + D233      

!      dg_e(3,0,0) = dg_e(3,0,0) + D300      
!      dg_e(3,0,1) = dg_e(3,0,1) + D301
!      dg_e(3,0,2) = dg_e(3,0,2) + D302
!      dg_e(3,0,3) = dg_e(3,0,3) + D303                              

      dg_e(3,0,0) = D300      
      dg_e(3,0,1) = D301
      dg_e(3,0,2) = D302
      dg_e(3,0,3) = D303                        
      dg_e(3,1,1) = dg_e(3,1,1) + D311
      dg_e(3,1,2) = dg_e(3,1,2) + D312
      dg_e(3,1,3) = dg_e(3,1,3) + D313
      dg_e(3,2,2) = dg_e(3,2,2) + D322
      dg_e(3,2,3) = dg_e(3,2,3) + D323
      dg_e(3,3,3) = dg_e(3,3,3) + D333

      K_e(1,1) = K_e(1,1) + K11
      K_e(1,2) = K_e(1,2) + K12
      K_e(1,3) = K_e(1,3) + K13
      K_e(2,2) = K_e(2,2) + K22
      K_e(2,3) = K_e(2,3) + K23
      K_e(3,3) = K_e(3,3) + K33
      
      
!#############################################################################                    
!########## INITIAL DATA FOR 1 BOSON STAR (1D) + 1 BOOSTED KERR ##############
!#############################################################################  

    elseif (idtype.eq.-12) then 
    
     !-------------------------------------------------------------------      
     ! compute the spacetime of the kerr black hole----------------------
     !-------------------------------------------------------------------      
      M = par(P_BH1_MASS)            
      vx = par(P_BH1_VX)
      vy = par(P_BH1_VY)                  

      x = xold - par(P_BH1_X0)
      y = yold - par(P_BH1_Y0)
      z = zold - par(P_BH1_Z0)      
      
      r = sqrt(x**2 + y**2 + z**2)
      if (r.lt.1E-2) r=1E-2      
            
#include "boostKS.inc"             

#include "KS_GH.inc"             
                      
     !-------------------------------------------------------------------           
     !----the addition of the black hole with the previous BS------------
     !-------------------------------------------------------------------      
     
     !! -- as a first approximation, we will assume than the lapse and shift-
     !! ---are given by the black hole --------------------------------------

!      g_e(0,0) = g_e(0,0) + g00 + 1.0d0
!      g_e(0,1) = g01 + g_e(0,1)
!      g_e(0,2) = g02 + g_e(0,2)
!      g_e(0,3) = g03 + g_e(0,3)                                           

      g_e(0,0) = g00
      g_e(0,1) = g01
      g_e(0,2) = g02
      g_e(0,3) = g03            
      g_e(1,1) = g_e(1,1) + g11 - 1.0d0
      g_e(1,2) = g_e(1,2) + g12
      g_e(1,3) = g_e(1,3) + g13
      g_e(2,2) = g_e(2,2) + g22 - 1.0d0
      g_e(2,3) = g_e(2,3) + g23
      g_e(3,3) = g_e(3,3) + g33 - 1.0d0

!      dg_e(1,0,0) = dg_e(1,0,0) + D100      
!      dg_e(1,0,1) = dg_e(1,0,1) + D101
!      dg_e(1,0,2) = dg_e(1,0,2) + D102
!      dg_e(1,0,3) = dg_e(1,0,3) + D103                  
      
      dg_e(1,0,0) = D100      
      dg_e(1,0,1) = D101
      dg_e(1,0,2) = D102
      dg_e(1,0,3) = D103                        
      dg_e(1,1,1) = dg_e(1,1,1) + D111
      dg_e(1,1,2) = dg_e(1,1,2) + D112
      dg_e(1,1,3) = dg_e(1,1,3) + D113
      dg_e(1,2,2) = dg_e(1,2,2) + D122
      dg_e(1,2,3) = dg_e(1,2,3) + D123
      dg_e(1,3,3) = dg_e(1,3,3) + D133

      
      
!      dg_e(2,0,0) = dg_e(2,0,0) + D200      
!      dg_e(2,0,1) = dg_e(2,0,1) + D201
!      dg_e(2,0,2) = dg_e(2,0,2) + D202
!      dg_e(2,0,3) = dg_e(2,0,3) + D203                                    
      
      dg_e(2,0,0) = D200      
      dg_e(2,0,1) = D201
      dg_e(2,0,2) = D202
      dg_e(2,0,3) = D203                              
      dg_e(2,1,1) = dg_e(2,1,1) + D211
      dg_e(2,1,2) = dg_e(2,1,2) + D212
      dg_e(2,1,3) = dg_e(2,1,3) + D213
      dg_e(2,2,2) = dg_e(2,2,2) + D222
      dg_e(2,2,3) = dg_e(2,2,3) + D223
      dg_e(2,3,3) = dg_e(2,3,3) + D233      

!      dg_e(3,0,0) = dg_e(3,0,0) + D300      
!      dg_e(3,0,1) = dg_e(3,0,1) + D301
!      dg_e(3,0,2) = dg_e(3,0,2) + D302
!      dg_e(3,0,3) = dg_e(3,0,3) + D303                              

      dg_e(3,0,0) = D300      
      dg_e(3,0,1) = D301
      dg_e(3,0,2) = D302
      dg_e(3,0,3) = D303                        
      dg_e(3,1,1) = dg_e(3,1,1) + D311
      dg_e(3,1,2) = dg_e(3,1,2) + D312
      dg_e(3,1,3) = dg_e(3,1,3) + D313
      dg_e(3,2,2) = dg_e(3,2,2) + D322
      dg_e(3,2,3) = dg_e(3,2,3) + D323
      dg_e(3,3,3) = dg_e(3,3,3) + D333

      K_e(1,1) = K_e(1,1) + K11
      K_e(1,2) = K_e(1,2) + K12
      K_e(1,3) = K_e(1,3) + K13
      K_e(2,2) = K_e(2,2) + K22
      K_e(2,3) = K_e(2,3) + K23
      K_e(3,3) = K_e(3,3) + K33
            
    !#################### INITIAL DATA FOR 1 BOSON STAR (3D)#############  
    elseif (idtype.eq.-2) then 

      alp = g_e(0,0)
      d1alp = dg_e(1,0,0)
      d2alp = dg_e(2,0,0)
      d3alp = dg_e(3,0,0)            

      g_e(0,0) = - alp*alp
      dg_e(1,0,0) = - 2.0*alp*d1alp
      dg_e(2,0,0) = - 2.0*alp*d2alp
      dg_e(3,0,0) = - 2.0*alp*d3alp            

      alp_b      = K_e(0,0)         
      K_e(0,0)   = 0.0d0
      dphic_e(0) = (omega/alp_b)*dphic_e(0)
          
    !#################### INITIAL DATA FOR 1 MASSLESS SF & 1 COMPLEX######  
    elseif (idtype.eq.-1) then 
     
      alp = g_e(0,0)
      
      g_e(0,0) = - alp*alp
                
    !#################### RANDOM INITIAL DATA IN TOP OF MINKOWSKI ############
    elseif (idtype.eq.1) then 
      
      if (ltrace) write(*,*)'@@@ runing idype 1 -----'
    
      mr_amp = par(P_MR_AMP)
      mr_shift = par(P_MR_SHIFT)      
      call random_number(harvest)
      harvest = 2.0*(harvest - 0.5)      

      g_e(0,0) = g_e(0,0) + mr_amp*harvest(1)
      g_e(0,1) = mr_shift + mr_amp*harvest(2)
      g_e(0,2) = mr_shift + mr_amp*harvest(3)
      g_e(0,3) = mr_shift + mr_amp*harvest(4)
      g_e(1,1) = g_e(1,1) + mr_amp*harvest(5)
      g_e(1,2) = g_e(1,2) + mr_amp*harvest(6)
      g_e(1,3) = g_e(1,3) + mr_amp*harvest(7)
      g_e(2,2) = g_e(2,2) + mr_amp*harvest(8)
      g_e(2,3) = g_e(2,3) + mr_amp*harvest(9)
      g_e(3,3) = g_e(3,3) + mr_amp*harvest(10)                                    

      dg_e(1,0,0) = dg_e(1,0,0) + mr_amp*harvest(11)
      dg_e(2,0,0) = dg_e(2,0,0) + mr_amp*harvest(12)
      dg_e(3,0,0) = dg_e(3,0,0) + mr_amp*harvest(13)
      
      dg_e(1,0,1) = dg_e(1,0,1) + mr_amp*harvest(14)            
      dg_e(1,0,2) = dg_e(1,0,2) + mr_amp*harvest(15)            
      dg_e(1,0,3) = dg_e(1,0,3) + mr_amp*harvest(16)                        
      dg_e(2,0,1) = dg_e(2,0,1) + mr_amp*harvest(17)            
      dg_e(2,0,2) = dg_e(2,0,2) + mr_amp*harvest(18)            
      dg_e(2,0,3) = dg_e(2,0,3) + mr_amp*harvest(19)                        
      dg_e(3,0,1) = dg_e(3,0,1) + mr_amp*harvest(20)            
      dg_e(3,0,2) = dg_e(3,0,2) + mr_amp*harvest(21)            
      dg_e(3,0,3) = dg_e(3,0,3) + mr_amp*harvest(22)                        
      
      dg_e(1,1,1) = dg_e(1,1,1) + mr_amp*harvest(23)
      dg_e(1,1,2) = dg_e(1,1,2) + mr_amp*harvest(24)
      dg_e(1,1,3) = dg_e(1,1,3) + mr_amp*harvest(25)
      dg_e(1,2,2) = dg_e(1,2,2) + mr_amp*harvest(26)
      dg_e(1,2,3) = dg_e(1,2,3) + mr_amp*harvest(27)
      dg_e(1,3,3) = dg_e(1,3,3) + mr_amp*harvest(28)                                    
      dg_e(2,1,1) = dg_e(2,1,1) + mr_amp*harvest(29)
      dg_e(2,1,2) = dg_e(2,1,2) + mr_amp*harvest(30)
      dg_e(2,1,3) = dg_e(2,1,3) + mr_amp*harvest(31)
      dg_e(2,2,2) = dg_e(2,2,2) + mr_amp*harvest(32)
      dg_e(2,2,3) = dg_e(2,2,3) + mr_amp*harvest(33)
      dg_e(2,3,3) = dg_e(2,3,3) + mr_amp*harvest(34)                                    
      dg_e(3,1,1) = dg_e(3,1,1) + mr_amp*harvest(35)
      dg_e(3,1,2) = dg_e(3,1,2) + mr_amp*harvest(36)
      dg_e(3,1,3) = dg_e(3,1,3) + mr_amp*harvest(37)
      dg_e(3,2,2) = dg_e(3,2,2) + mr_amp*harvest(38)
      dg_e(3,2,3) = dg_e(3,2,3) + mr_amp*harvest(39)
      dg_e(3,3,3) = dg_e(3,3,3) + mr_amp*harvest(40)                                    
      
      H_e(0) = H_e(0) + mr_amp*harvest(41)     
      H_e(1) = H_e(1) + mr_amp*harvest(42)     
      H_e(2) = H_e(2) + mr_amp*harvest(43)     
      H_e(3) = H_e(3) + mr_amp*harvest(44)                 
      
      K_e(0,0) = K_e(0,0) + mr_amp*harvest(45)
      K_e(0,1) = K_e(0,1) + mr_amp*harvest(46)
      K_e(0,2) = K_e(0,2) + mr_amp*harvest(47)
      K_e(0,3) = K_e(0,3) + mr_amp*harvest(48)                        
      K_e(1,1) = K_e(1,1) + mr_amp*harvest(49)
      K_e(1,2) = K_e(1,2) + mr_amp*harvest(50)
      K_e(1,3) = K_e(1,3) + mr_amp*harvest(51)
      K_e(2,2) = K_e(2,2) + mr_amp*harvest(52)
      K_e(2,3) = K_e(2,3) + mr_amp*harvest(53)
      K_e(3,3) = K_e(3,3) + mr_amp*harvest(54)                                    
      
      phim_e = phim_e + mr_amp*harvest(55)
      phir_e = phir_e + mr_amp*harvest(56)
      phic_e = phic_e + mr_amp*harvest(57)            

      dphim_e(0) = dphim_e(0) + mr_amp*harvest(58)
      dphim_e(1) = dphim_e(1) + mr_amp*harvest(59)
      dphim_e(2) = dphim_e(2) + mr_amp*harvest(60)
      dphim_e(3) = dphim_e(3) + mr_amp*harvest(61)
            
      dphir_e(0) = dphir_e(0) + mr_amp*harvest(62)
      dphir_e(1) = dphir_e(1) + mr_amp*harvest(63)
      dphir_e(2) = dphir_e(2) + mr_amp*harvest(64)
      dphir_e(3) = dphir_e(3) + mr_amp*harvest(65)

      dphic_e(0) = dphic_e(0) + mr_amp*harvest(66)
      dphic_e(1) = dphic_e(1) + mr_amp*harvest(67)
      dphic_e(2) = dphic_e(2) + mr_amp*harvest(68)
      dphic_e(3) = dphic_e(3) + mr_amp*harvest(69)

      dH_e(0) = dH_e(0) + mr_amp*harvest(70)
      dH_e(1) = dH_e(1) + mr_amp*harvest(71)
      dH_e(2) = dH_e(2) + mr_amp*harvest(72)
      dH_e(3) = dH_e(3) + mr_amp*harvest(73)
                  
                                                    
!#################### FLAT IN TWISTED COORDINATES #######################
   elseif (idtype.eq.2) then 
      a = par(P_TW_AMP)
      b = par(P_TW_WIDTH)
      nx = par(P_NX)
      dx = par(P_DX)
      Lx = (nx- 2)*dx
      c = 2*pi/Lx 
        
!#include "messytwisted_3D.inc"             
!#include "messytwisted_2D.inc"     

      g_e(1,1) = g11
      g_e(1,2) = g12
      g_e(1,3) = g13
      g_e(2,2) = g22
      g_e(2,3) = g23
      g_e(3,3) = g33

      dg_e(1,1,1) = D111
      dg_e(1,1,2) = D112
      dg_e(1,1,3) = D113
      dg_e(1,2,2) = D122
      dg_e(1,2,3) = D123
      dg_e(1,3,3) = D133      
      dg_e(2,1,1) = D211
      dg_e(2,1,2) = D212
      dg_e(2,1,3) = D213
      dg_e(2,2,2) = D222
      dg_e(2,2,3) = D223
      dg_e(2,3,3) = D233      
      dg_e(3,1,1) = D311
      dg_e(3,1,2) = D312
      dg_e(3,1,3) = D313
      dg_e(3,2,2) = D322
      dg_e(3,2,3) = D323
      dg_e(3,3,3) = D333

!#################### KERR-SCHILD BH WITHOUT SPIN #######################            
    elseif (idtype.eq.3) then 

      M = par(P_BH1_MASS)            
      
      x = xold - par(P_BH1_X0)
      y = yold - par(P_BH1_Y0)
      z = zold - par(P_BH1_Z0)      
      
      r = sqrt(x**2 + y**2 + z**2)
      if (r.lt.1E-2) r=1E-2
            
      g_e(0,0) = -1.0 + 2.0*M/r
      g_e(0,1) = 2.0*M*x/r**2
      g_e(0,2) = 2.0*M*y/r**2
      g_e(0,3) = 2.0*M*z/r**2
      g_e(1,1) = 1.0 + 2.0*M*x**2/r**3
      g_e(1,2) = 2.0*M*x*y/r**3
      g_e(1,3) = 2.0*M*x*z/r**3
      g_e(2,2) = 1.0 + 2.0*M*y**2/r**3
      g_e(2,3) = 2.0*M*y*z/r**3
      g_e(3,3) = 1.0 + 2.0*M*z**2/r**3

      dg_e(1,0,0) = -2.0*M*x/r**3
      dg_e(2,0,0) = -2.0*M*y/r**3
      dg_e(3,0,0) = -2.0*M*z/r**3		

      dg_e(1,0,1) = 2.0*M*(r*r - 2.0*x*x)/r**4
      dg_e(1,0,2) = -4.0*M*x*y/r**4
      dg_e(1,0,3) = -4.0*M*x*z/r**4
      dg_e(2,0,1) = -4.0*M*y*x/r**4	
      dg_e(2,0,2) = 2.0*M*(r*r - 2.0*y*y)/r**4
      dg_e(2,0,3) = -4.0*M*y*z/r**4
      dg_e(3,0,1) = -4.0*M*z*x/r**4
      dg_e(3,0,2) = -4.0*M*z*y/r**4		
      dg_e(3,0,3) = 2.0*M*(r*r - 2.0*z*z)/r**4	

      dg_e(1,1,1) = 2.0*M*(x*r**3 + x*r**3 - 3.0*r*x*x*x)/r**6
      dg_e(1,1,2) = 2.0*M*(y*r**3 - 3.0*r*x*x*y)/r**6
      dg_e(1,1,3) = 2.0*M*(z*r**3 - 3.0*r*x*x*z)/r**6
      dg_e(1,2,2) = 2.0*M*(- 3.0*r*x*y*y)/r**6
      dg_e(1,2,3) = 2.0*M*(- 3.0*r*x*y*z)/r**6
      dg_e(1,3,3) = 2.0*M*(- 3.0*r*x*z*z)/r**6
      dg_e(2,1,1) = 2.0*M*(- 3.0*r*y*x*x)/r**6
      dg_e(2,1,2) = 2.0*M*(x*r**3 - 3.0*r*y*x*y)/r**6
      dg_e(2,1,3) = 2.0*M*(- 3.0*r*y*x*z)/r**6
      dg_e(2,2,2) = 2.0*M*(y*r**3 + y*r**3 - 3.0*r*y*y*y)/r**6
      dg_e(2,2,3) = 2.0*M*(z*r**3 - 3.0*r*y*y*z)/r**6
      dg_e(2,3,3) = 2.0*M*(- 3.0*r*y*z*z)/r**6
      dg_e(3,1,1) = 2.0*M*(- 3.0*r*z*x*x)/r**6
      dg_e(3,1,2) = 2.0*M*(- 3.0*r*z*x*y)/r**6
      dg_e(3,1,3) = 2.0*M*(x*r**3 - 3.0*r*z*x*z)/r**6
      dg_e(3,2,2) = 2.0*M*(- 3.0*r*z*y*y)/r**6
      dg_e(3,2,3) = 2.0*M*(y*r**3 - 3.0*r*z*y*z)/r**6
      dg_e(3,3,3) = 2.0*M*(z*r**3 + z*r**3 - 3.0*r*z*z*z)/r**6
	
!-----the standard extrinsic curvature---------------------
      K_e(1,1) = -2.D0*((M/r+2)*x**2-x**2-y**2-z**2)*M/(x*&
     &*2+y**2+z**2)**2/sqrt(1.D0+2.D0*M/r)
      K_e(1,2) = -2.D0*(M/r+2)*x*y*M/(r**2)**2/s&
     &qrt(1.D0+2.D0*M/r)
      K_e(1,3) = -2.D0*(M/r+2)*x*z*M/(r**2)**2/s&
     &qrt(1.D0+2.D0*M/r)
      K_e(2,2) = -2.D0*((M/r+2)*y**2-x**2-y**2-z**2)*M/(x*&
     &*2+y**2+z**2)**2/sqrt(1.D0+2.D0*M/r)
      K_e(2,3) = -2.D0*(M/r+2)*y*z*M/(r**2)**2/s&
     &qrt(1.D0+2.D0*M/r)
      K_e(3,3) = -2.D0*((M/r+2)*z**2-x**2-y**2-z**2)*M/(x*&
     &*2+y**2+z**2)**2/sqrt(1.D0+2.D0*M/r)

     
!#################### BH IN ISOTROPICS  #########################
    elseif (idtype.eq.4) then 

      M = par(P_BH1_MASS)
      matching = par(P_BH1_MATCHING)
      ro = M/2.0d0

      if (r.lt.1E-4) r=1E-4
      
      psi = (1.0d0 + M/(2.0*r))
      drpsi = - M/(2.0*r**2)      
                  
!     Dressing (2 derivatives continuous)            
      if ( (matching .EQ. 1) .AND. (r .LT. ro) )  then
	psi = (1.0 + M/ro) - (M/ro**3)*(r**2) + (M/(2*ro**4))*(r**3) 
	drpsi = -(2.0*M/(ro**3))*r + ((3.0*M)/(2.0*ro**4))*(r**2)
      end if
                                    
      g_e(1,1) = psi**4
      g_e(2,2) = g_e(1,1)
      g_e(3,3) = g_e(1,1)

      dg_e(1,1,1) = 4.0*(psi**3)*(x/r)*drpsi
      dg_e(1,2,2) = dg_e(1,1,1)
      dg_e(1,3,3) = dg_e(1,1,1)
      dg_e(2,1,1) = 4.0*(psi**3)*(y/r)*drpsi
      dg_e(2,2,2) = dg_e(2,1,1)
      dg_e(2,3,3) = dg_e(2,1,1)
      dg_e(3,1,1) = 4.0*(psi**3)*(z/r)*drpsi
      dg_e(3,2,2) = dg_e(3,1,1)
      dg_e(3,3,3) = dg_e(3,1,1)      
     
     
!#################### KERR-SCHILD BH WITH SPIN #######################                          
    elseif (idtype.eq.5) then 
      
      M     = par(P_BH1_MASS)            
      a     = par(P_BH1_SPIN)     
      theta = par(P_BH1_SPIN_THETA)             

      x = xold - par(P_BH1_X0)
      y = yold - par(P_BH1_Y0)
      z = zold - par(P_BH1_Z0)      
       
      if ((abs(x).lt.1E-2) .AND. (abs(y).lt.1E-2)&
     &    .AND. (abs(z).lt.1E-2)) then
      
        x=1E-2
        y=1E-2
        z=1E-2
      
      end if

      r = sqrt(x**2 + y**2 + z**2)
      if (r.lt.1E-2) r=1E-2
#include "spinKS.inc" 
                        
!#include "spinorientedKS.inc"                                           
            
#include "KS_GH.inc"             

      g_e(0,0) = g00
      g_e(0,1) = g01
      g_e(0,2) = g02
      g_e(0,3) = g03                                           
      g_e(1,1) = g11
      g_e(1,2) = g12
      g_e(1,3) = g13
      g_e(2,2) = g22
      g_e(2,3) = g23
      g_e(3,3) = g33

      dg_e(1,0,0) = D100
      dg_e(1,0,1) = D101
      dg_e(1,0,2) = D102
      dg_e(1,0,3) = D103                  
      dg_e(1,1,1) = D111
      dg_e(1,1,2) = D112
      dg_e(1,1,3) = D113
      dg_e(1,2,2) = D122
      dg_e(1,2,3) = D123
      dg_e(1,3,3) = D133

      dg_e(2,0,0) = D200
      dg_e(2,0,1) = D201
      dg_e(2,0,2) = D202
      dg_e(2,0,3) = D203                                    
      dg_e(2,1,1) = D211
      dg_e(2,1,2) = D212
      dg_e(2,1,3) = D213
      dg_e(2,2,2) = D222
      dg_e(2,2,3) = D223
      dg_e(2,3,3) = D233      

      dg_e(3,0,0) = D300
      dg_e(3,0,1) = D301
      dg_e(3,0,2) = D302
      dg_e(3,0,3) = D303                              
      dg_e(3,1,1) = D311
      dg_e(3,1,2) = D312
      dg_e(3,1,3) = D313
      dg_e(3,2,2) = D322
      dg_e(3,2,3) = D323
      dg_e(3,3,3) = D333

      K_e(1,1) = K11
      K_e(1,2) = K12
      K_e(1,3) = K13
      K_e(2,2) = K22
      K_e(2,3) = K23
      K_e(3,3) = K33
      
 !#################### BOOSTED BH ##################################           
     
    elseif (idtype.eq.6) then 
      
      M = par(P_BH1_MASS)            
      vx = par(P_BH1_VX)
      vy = par(P_BH1_VY)                  

      x = xold - par(P_BH1_X0)
      y = yold - par(P_BH1_Y0)
      z = zold - par(P_BH1_Z0)      
      
      r = sqrt(x**2 + y**2 + z**2)
      if (r.lt.1E-2) r=1E-2      
            
#include "boostKS.inc"             

#include "KS_GH.inc"             

      g_e(0,0) = g00
      g_e(0,1) = g01
      g_e(0,2) = g02
      g_e(0,3) = g03                                           
      g_e(1,1) = g11
      g_e(1,2) = g12
      g_e(1,3) = g13
      g_e(2,2) = g22
      g_e(2,3) = g23
      g_e(3,3) = g33

      dg_e(1,0,0) = D100
      dg_e(1,0,1) = D101
      dg_e(1,0,2) = D102
      dg_e(1,0,3) = D103                  
      dg_e(1,1,1) = D111
      dg_e(1,1,2) = D112
      dg_e(1,1,3) = D113
      dg_e(1,2,2) = D122
      dg_e(1,2,3) = D123
      dg_e(1,3,3) = D133

      dg_e(2,0,0) = D200
      dg_e(2,0,1) = D201
      dg_e(2,0,2) = D202
      dg_e(2,0,3) = D203                                    
      dg_e(2,1,1) = D211
      dg_e(2,1,2) = D212
      dg_e(2,1,3) = D213
      dg_e(2,2,2) = D222
      dg_e(2,2,3) = D223
      dg_e(2,3,3) = D233      

      dg_e(3,0,0) = D300
      dg_e(3,0,1) = D301
      dg_e(3,0,2) = D302
      dg_e(3,0,3) = D303                              
      dg_e(3,1,1) = D311
      dg_e(3,1,2) = D312
      dg_e(3,1,3) = D313
      dg_e(3,2,2) = D322
      dg_e(3,2,3) = D323
      dg_e(3,3,3) = D333

      K_e(1,1) = K11
      K_e(1,2) = K12
      K_e(1,3) = K13
      K_e(2,2) = K22
      K_e(2,3) = K23
      K_e(3,3) = K33

!#################### 1+1 KERR BLACK HOLES #######################              
    elseif (idtype.eq.7.or.idtype.eq.17) then 
    
!     the first black hole      
            
      M = par(P_BH1_MASS)            
      a = par(P_BH1_SPIN)                  
      vx = par(P_BH1_VX)
      vy = par(P_BH1_VY)

      x = xold - par(P_BH1_X0)
      y = yold - par(P_BH1_Y0)
      z = zold - par(P_BH1_Z0)      
      
      r = sqrt(x**2 + y**2 + z**2)
      if (r.lt.1E-2) r=1E-2
     
if(idtype.gt.10) then 
#include "boostKS.inc"
else                      
#include "spinKS.inc"             
end if

#include "KS_GH.inc"             

      g_e(0,0) = g00
      g_e(0,1) = g01
      g_e(0,2) = g02
      g_e(0,3) = g03                                           
      g_e(1,1) = g11
      g_e(1,2) = g12
      g_e(1,3) = g13
      g_e(2,2) = g22
      g_e(2,3) = g23
      g_e(3,3) = g33

      dg_e(1,0,0) = D100
      dg_e(1,0,1) = D101
      dg_e(1,0,2) = D102
      dg_e(1,0,3) = D103                  
      dg_e(1,1,1) = D111
      dg_e(1,1,2) = D112
      dg_e(1,1,3) = D113
      dg_e(1,2,2) = D122
      dg_e(1,2,3) = D123
      dg_e(1,3,3) = D133

      dg_e(2,0,0) = D200
      dg_e(2,0,1) = D201
      dg_e(2,0,2) = D202
      dg_e(2,0,3) = D203                                    
      dg_e(2,1,1) = D211
      dg_e(2,1,2) = D212
      dg_e(2,1,3) = D213
      dg_e(2,2,2) = D222
      dg_e(2,2,3) = D223
      dg_e(2,3,3) = D233      

      dg_e(3,0,0) = D300
      dg_e(3,0,1) = D301
      dg_e(3,0,2) = D302
      dg_e(3,0,3) = D303                              
      dg_e(3,1,1) = D311
      dg_e(3,1,2) = D312
      dg_e(3,1,3) = D313
      dg_e(3,2,2) = D322
      dg_e(3,2,3) = D323
      dg_e(3,3,3) = D333

      K_e(1,1) = K11
      K_e(1,2) = K12
      K_e(1,3) = K13
      K_e(2,2) = K22
      K_e(2,3) = K23
      K_e(3,3) = K33

      ! compute the lapseA and the AKud            
      h11 = g11
      h12 = g12
      h13 = g13
      h22 = g22
      h23 = g23
      h33 = g33
                                          
      huu11 = -h23**2 + h22*h33
      huu12 = h13*h23 - h12*h33
      huu13 = -(h13*h22) + h12*h23
      huu22 = -h13**2 + h11*h33
      huu23 = h12*h13 - h11*h23
      huu33 = -h12**2 + h11*h22
      deth = h11*huu11 + h12*huu12 + h13*huu13
      huu11 = huu11/deth
      huu12 = huu12/deth
      huu13 = huu13/deth
      huu22 = huu22/deth
      huu23 = huu23/deth
      huu33 = huu33/deth
           
      AKud11 = huu11*K11 + huu12*K12 + huu13*K13
      AKud12 = huu11*K12 + huu12*K22 + huu13*K23
      AKud13 = huu11*K13 + huu12*K23 + huu13*K33
      AKud21 = huu12*K11 + huu22*K12 + huu23*K13
      AKud22 = huu12*K12 + huu22*K22 + huu23*K23
      AKud23 = huu12*K13 + huu22*K23 + huu23*K33
      AKud31 = huu13*K11 + huu23*K12 + huu33*K13
      AKud32 = huu13*K12 + huu23*K22 + huu33*K23
      AKud33 = huu13*K13 + huu23*K23 + huu33*K33
      
      alpA   = 1.0d0/sqrt(1.0d0 + 2.0d0*M/r) 
      dxalpA = M*(alpA**3)*x/(r**3)
      dyalpA = M*(alpA**3)*y/(r**3)      
      dzalpA = M*(alpA**3)*z/(r**3)       

            
!     the second black hole      
      
      M = par(P_BH2_MASS)            
      a = par(P_BH2_SPIN)                  
      vx = -par(P_BH1_VX)
      vy = -par(P_BH1_VY)

      x = xold - par(P_BH2_X0)
      y = yold - par(P_BH2_Y0)
      z = zold - par(P_BH2_Z0)      
      
      r = sqrt(x**2 + y**2 + z**2)
      if (r.lt.1E-2) r=1E-2

if(idtype.gt.10) then
#include "boostKS.inc"
else
#include "spinKS.inc"
end if                              
        
#include "KS_GH.inc"             

      g_e(0,1) = g_e(0,1) + g01
      g_e(0,2) = g_e(0,2) + g02
      g_e(0,3) = g_e(0,3) + g03                                           
      g_e(1,1) = g_e(1,1) + g11 - 1.0d0
      g_e(1,2) = g_e(1,2) + g12
      g_e(1,3) = g_e(1,3) + g13
      g_e(2,2) = g_e(2,2) + g22 - 1.0d0
      g_e(2,3) = g_e(2,3) + g23
      g_e(3,3) = g_e(3,3) + g33 - 1.0d0

      dg_e(1,0,1) = dg_e(1,0,1) + D101
      dg_e(1,0,2) = dg_e(1,0,2) + D102
      dg_e(1,0,3) = dg_e(1,0,3) + D103                  
      dg_e(1,1,1) = dg_e(1,1,1) + D111
      dg_e(1,1,2) = dg_e(1,1,2) + D112
      dg_e(1,1,3) = dg_e(1,1,3) + D113
      dg_e(1,2,2) = dg_e(1,2,2) + D122
      dg_e(1,2,3) = dg_e(1,2,3) + D123
      dg_e(1,3,3) = dg_e(1,3,3) + D133

      dg_e(2,0,1) = dg_e(2,0,1) + D201
      dg_e(2,0,2) = dg_e(2,0,2) + D202
      dg_e(2,0,3) = dg_e(2,0,3) + D203                                    
      dg_e(2,1,1) = dg_e(2,1,1) + D211
      dg_e(2,1,2) = dg_e(2,1,2) + D212
      dg_e(2,1,3) = dg_e(2,1,3) + D213
      dg_e(2,2,2) = dg_e(2,2,2) + D222
      dg_e(2,2,3) = dg_e(2,2,3) + D223
      dg_e(2,3,3) = dg_e(2,3,3) + D233      

      dg_e(3,0,1) = dg_e(3,0,1) + D301
      dg_e(3,0,2) = dg_e(3,0,2) + D302
      dg_e(3,0,3) = dg_e(3,0,3) + D303                              
      dg_e(3,1,1) = dg_e(3,1,1) + D311
      dg_e(3,1,2) = dg_e(3,1,2) + D312
      dg_e(3,1,3) = dg_e(3,1,3) + D313
      dg_e(3,2,2) = dg_e(3,2,2) + D322
      dg_e(3,2,3) = dg_e(3,2,3) + D323
      dg_e(3,3,3) = dg_e(3,3,3) + D333

      K_e(1,1) = K_e(1,1) + K11
      K_e(1,2) = K_e(1,2) + K12
      K_e(1,3) = K_e(1,3) + K13
      K_e(2,2) = K_e(2,2) + K22
      K_e(2,3) = K_e(2,3) + K23
      K_e(3,3) = K_e(3,3) + K33

      ! compute the lapseB and the BKud            
      h11 = g11
      h12 = g12
      h13 = g13
      h22 = g22
      h23 = g23
      h33 = g33
                                          
      huu11 = -h23**2 + h22*h33
      huu12 = h13*h23 - h12*h33
      huu13 = -(h13*h22) + h12*h23
      huu22 = -h13**2 + h11*h33
      huu23 = h12*h13 - h11*h23
      huu33 = -h12**2 + h11*h22
      deth = h11*huu11 + h12*huu12 + h13*huu13
      huu11 = huu11/deth
      huu12 = huu12/deth
      huu13 = huu13/deth
      huu22 = huu22/deth
      huu23 = huu23/deth
      huu33 = huu33/deth
           
      BKud11 = huu11*K11 + huu12*K12 + huu13*K13
      BKud12 = huu11*K12 + huu12*K22 + huu13*K23
      BKud13 = huu11*K13 + huu12*K23 + huu13*K33
      BKud21 = huu12*K11 + huu22*K12 + huu23*K13
      BKud22 = huu12*K12 + huu22*K22 + huu23*K23
      BKud23 = huu12*K13 + huu22*K23 + huu23*K33
      BKud31 = huu13*K11 + huu23*K12 + huu33*K13
      BKud32 = huu13*K12 + huu23*K22 + huu33*K23
      BKud33 = huu13*K13 + huu23*K23 + huu33*K33
            
      alpB   = 1.0d0/sqrt(1.0d0 + 2.0d0*M/r) 
      dxalpB = M*(alpB**3)*x/(r**3)
      dyalpB = M*(alpB**3)*y/(r**3)      
      dzalpB = M*(alpB**3)*z/(r**3)       

                  
 !    now the lapse addition
      alp2 = 1.0d0/( alpA**(-2) + alpB**(-2) - 1.0d0) 
      alp = sqrt(alp2)
      
      h11 = g_e(1,1)
      h12 = g_e(1,2)
      h13 = g_e(1,3)
      h22 = g_e(2,2)
      h23 = g_e(2,3)
      h33 = g_e(3,3)
                                          
      huu11 = -h23**2 + h22*h33
      huu12 = h13*h23 - h12*h33
      huu13 = -(h13*h22) + h12*h23
      huu22 = -h13**2 + h11*h33
      huu23 = h12*h13 - h11*h23
      huu33 = -h12**2 + h11*h22
      deth = h11*huu11 + h12*huu12 + h13*huu13
      huu11 = huu11/deth
      huu12 = huu12/deth
      huu13 = huu13/deth
      huu22 = huu22/deth
      huu23 = huu23/deth
      huu33 = huu33/deth
       
      b1 = g_e(0,1)*huu11 + g_e(0,2)*huu12 + g_e(0,3)*huu13
      b2 = g_e(0,1)*huu12 + g_e(0,2)*huu22 + g_e(0,3)*huu23
      b3 = g_e(0,1)*huu13 + g_e(0,2)*huu23 + g_e(0,3)*huu33
      g_e(0,0) = -alp2 + b1**2*h11 + b2**2*h22 + 2*(b1*(b2*h12 + b3*h13) &
     &      + b2*b3*h23) + b3**2*h33
      
!     now the extrinsic curvature addition
      eKud11 = AKud11 + BKud11
      eKud12 = AKud12 + BKud12
      eKud13 = AKud13 + BKud13
      eKud21 = AKud21 + BKud21
      eKud22 = AKud22 + BKud22
      eKud23 = AKud23 + BKud23
      eKud31 = AKud31 + BKud31
      eKud32 = AKud32 + BKud32
      eKud33 = AKud33 + BKud33
                                                
      K_e(1,1) = h11*eKud11 + h12*eKud21 + h13*eKud31 
      K_e(2,2) = h12*eKud12 + h22*eKud22 + h23*eKud32 
      K_e(3,3) = h13*eKud13 + h23*eKud23 + h33*eKud33 
      
      K_e(1,2) = h11*eKud12 + h12*eKud22 + h13*eKud32 
      K_e(1,3) = h11*eKud13 + h12*eKud23 + h13*eKud33 
      K_e(2,3) = h12*eKud13 + h22*eKud23 + h23*eKud33                               

!     compute now the space derivatives of g00
      dxalp = (alp**3)*(dxalpA/(alpA**3) + dxalpB/(alpB**3))
      dyalp = (alp**3)*(dyalpA/(alpA**3) + dyalpB/(alpB**3))
      dzalp = (alp**3)*(dzalpA/(alpA**3) + dzalpB/(alpB**3))                 

      dg_e(1,0,0) = -2.0*alp*dxalp  &
     &  + 2.0*( b1*dg_e(1,0,1) + b2*dg_e(1,0,2) + b3*dg_e(1,0,3) ) &
     &  - ( b1*b1*dg_e(1,1,1) + b2*b2*dg_e(1,2,2) + b3*b3*dg_e(1,3,3) &
     &    + 2.0d0*(b1*b2*dg_e(1,1,2) + b1*b3*dg_e(1,1,3) + b2*b3*dg_e(1,2,3)) )
      
      dg_e(2,0,0) = -2.0*alp*dyalp  &
     &  + 2.0*( b1*dg_e(2,0,1) + b2*dg_e(2,0,2) + b3*dg_e(2,0,3) ) &
     &  - ( b1*b1*dg_e(2,1,1) + b2*b2*dg_e(2,2,2) + b3*b3*dg_e(2,3,3) &
     &    + 2.0d0*(b1*b2*dg_e(2,1,2) + b1*b3*dg_e(2,1,3) + b2*b3*dg_e(2,2,3)) )

      dg_e(3,0,0) = -2.0*alp*dzalp  &
     &  + 2.0*( b1*dg_e(3,0,1) + b2*dg_e(3,0,2) + b3*dg_e(3,0,3) ) &
     &  - ( b1*b1*dg_e(3,1,1) + b2*b2*dg_e(3,2,2) + b3*b3*dg_e(3,3,3) &
     &    + 2.0d0*(b1*b2*dg_e(3,1,2) + b1*b3*dg_e(3,1,3) + b2*b3*dg_e(3,2,3)) )
                           
!######################################################################
!#################### 1+1 ISOTROPIC BLACK HOLES #######################              
!######################################################################

    elseif (idtype.eq.8) then 
    
!     the first black hole --------------     

      M = par(P_BH1_MASS)
      matching = par(P_BH1_MATCHING)
      ro = M/2.0d0

      x = xold - par(P_BH1_X0)
      y = yold - par(P_BH1_Y0)
      z = zold - par(P_BH1_Z0)
      r = sqrt( x**2 + y**2 + z**2 )

      if (r.lt.1E-4) r=1E-4
      
      psi1 = (1.0d0 + M/(2.0*r))
      drpsi1 = - M/(2.0*r**2)      

!     Dressing (2 derivatives continuous)            
      if ( (matching .EQ. 1) .AND. (r .LT. ro) )  then
	psi1 = (1.0 + M/ro) - (M/ro**3)*(r**2) + (M/(2*ro**4))*(r**3) 
	drpsi1 = -(2.0*M/(ro**3))*r + ((3.0*M)/(2.0*ro**4))*(r**2)
      end if

      dg_e(1,1,1) = 4.0*(psi1**3)*(x/r)*drpsi1
      dg_e(2,1,1) = 4.0*(psi1**3)*(y/r)*drpsi1
      dg_e(3,1,1) = 4.0*(psi1**3)*(z/r)*drpsi1


!     the second black hole -------------     

      M = par(P_BH2_MASS)
      matching = par(P_BH2_MATCHING)
      ro = M/2.0d0

      x = xold - par(P_BH2_X0)
      y = yold - par(P_BH2_Y0)
      z = zold - par(P_BH2_Z0)
      r = sqrt( x**2 + y**2 + z**2 )
      if (r.lt.1E-4) r=1E-4
      
      psi2 = (1.0d0 + M/(2.0*r))
      drpsi2 = - M/(2.0*r**2)      
          
!     Dressing (2 derivatives continuous)            
      if ( (matching .EQ. 1) .AND. (r .LT. ro) )  then
	psi2 = (1.0 + M/ro) - (M/ro**3)*(r**2) + (M/(2*ro**4))*(r**3) 
	drpsi2 = -(2.0*M/(ro**3))*r + ((3.0*M)/(2.0*ro**4))*(r**2)
      end if
                                    
!     the addition of the black holes -------------     

      psiT   = psi1 + psi2 - 1.0d0
      drpsiT = drpsi1 + drpsi2

      g_e(1,1) = psiT**4
      g_e(2,2) = g_e(1,1)
      g_e(3,3) = g_e(1,1)

      dg_e(1,1,1) = dg_e(1,1,1) + 4.0*(psi2**3)*(x/r)*drpsi2
      dg_e(2,1,1) = dg_e(2,1,1) + 4.0*(psi2**3)*(y/r)*drpsi2
      dg_e(3,1,1) = dg_e(3,1,1) + 4.0*(psi2**3)*(z/r)*drpsi2

      dg_e(1,2,2) = dg_e(1,1,1)
      dg_e(1,3,3) = dg_e(1,1,1)
      dg_e(2,2,2) = dg_e(2,1,1)
      dg_e(2,3,3) = dg_e(2,1,1)
      dg_e(3,2,2) = dg_e(3,1,1)
      dg_e(3,3,3) = dg_e(3,1,1)

!#######################################################################       
!#################### BINARY BLACK HOLE FROM LORENE ####################       
!#######################################################################

    else if (idtype.eq.-50) then 
!!MM: g0a are computed (differently) in subroutine read_lorene_bbh
!
!       ! just compute the g0a from the lapse and the shift
!       ! which are saved in g0a precisely

!       alp = g_e(0,0)
!       b1  = g_e(0,1)
!       b2  = g_e(0,2)
!       b3  = g_e(0,3)

!       g_e(0,0) = -alp*alp +  b1*b1*g_e(1,1) + b2*b2*g_e(2,2) + b3*b3*g_e(3,3) &
!    &           + 2.0*(b1*(b2*g_e(1,2) + b3*g_e(1,3))+ b2*b3*g_e(2,3))
!       g_e(0,1) = g_e(1,1)*b1 + g_e(1,2)*b2 + g_e(1,3)*b3
!       g_e(0,2) = g_e(1,2)*b1 + g_e(2,2)*b2 + g_e(2,3)*b3
!       g_e(0,3) = g_e(1,3)*b1 + g_e(2,3)*b2 + g_e(3,3)*b3
!!MM


!#######################################################################
!#################### TEUKOLSKY WAVES ##################################
!#######################################################################

    else if ((idtype.eq.100) .OR. (idtype.eq.101)) then 

      if ((round .eq. 1) .AND. (par(P_TW_AMP) .gt. 1E-12))  then
      
!------ teukolsky wave -------------------------------

        p = 1.0
        amp  = par(P_TW_AMP)      
        lam  = par(P_TW_WIDTH)
        mult = int(par(P_TW_MULTIPOLE))
        a    = par(P_TW_TIME_SHIFT)
            
        eps = -1.0d0
        if (idtype.eq.101) then
	  eps = 1.0d0
        end if
        a = -eps*a
	
        x = xold 
        y = yold 
        z = zold 
     
        r = sqrt(x**2 + y**2 + z**2)
        
        if (r.gt.(0.5*par(P_DX))) then	   
        
          if (sqrt(x**2+y**2).gt.(0.5*par(P_DX))) then

            theta = acos(z/r)
            phi = atan2(y,x)

            sinteta = sqrt(x**2+y**2)/r
            costeta = z/r

            dx_r = x/r
            dy_r = y/r
            dz_r = z/r

            sinphi = y / sqrt(x**2+y**2)
            cosphi = x / sqrt(x**2+y**2)
            dx_phi = -y/( x*x + y*y )
            dy_phi =  x/( x*x + y*y )
            dz_phi = 0.0
     
            dxx_phi = 2.*y*x/(x**2+y**2)**2
  	    dxy_phi = -(x**2-y**2)/(x**2+y**2)**2
	    dxz_phi = 0.0
	    dyy_phi = -2*y*x/(x**2+y**2)**2
	    dyz_phi = 0.0
	    dzz_phi = 0.0

            dx_theta = x*z/( r**2*sqrt( x**2 + y**2) )
            dy_theta = y*z/( r**2*sqrt( x**2 + y**2) )
            dz_theta = -sqrt( x*x + y*y )/(r*r)

            dxx_theta = -z*(2.*x**4+x**2*y**2-y**4-z**2*y**2)/r**4/(x**2+y**2)**1.5
            dxy_theta = -x*y*z*(3.*x**2+3.*y**2+z**2)/r**4/(x**2+y**2)**1.5
            dyx_theta = dxy_theta
            dxz_theta = x*(x**2+y**2-z**2)/r**4/(x**2+y**2)**0.5
            dzx_theta = dxz_theta
            dyy_theta = z*(x**4-x**2*y**2-2.*y**4+x**2*z**2)/r**4/(x**2+y**2)**1.5
            dzy_theta = y*(x**2+y**2-z**2)/r**4/(x**2+y**2)**0.5
            dyz_theta = dzy_theta

	
	    dyx_phi = dxy_phi
	    dzx_phi = dxz_phi
	    dzy_phi = dyz_phi

	    dxx_r = (y**2+z**2)/r**3
	    dyy_r = (x**2+z**2)/r**3
	    dzz_r = (x**2+y**2)/r**3
	    dxy_r = -(y*x)/r**3
	    dyx_r = -(y*x)/r**3
	    dxz_r = -(z*x)/r**3
	    dzx_r = -(z*x)/r**3
	    dyz_r = -(y*z)/r**3
	    dzy_r = -(y*z)/r**3

	    dzz_theta = 2*z*(x**2+y**2)**0.5/r**4
 
          else
	
	    r = abs(z)
	    p = 0.0

	  end if
	                                      
#include "TEUK.in"            

          g_e(0,0) = gC(4,4)
          g_e(0,1) = gC(1,4) 
          g_e(0,2) = gc(2,4)
          g_e(0,3) = gc(3,4)                                         
          g_e(1,1) = gC(1,1)
          g_e(1,2) = gC(1,2)
          g_e(1,3) = gC(1,3)
          g_e(2,2) = gC(2,2)
          g_e(2,3) = gC(2,3)
          g_e(3,3) = gC(3,3)
      
          dg_e(1:3,0,0) = DgC(1:3,4,4)
          dg_e(1:3,0,1) = DgC(1:3,1,4)
          dg_e(1:3,0,2) = Dgc(1:3,2,4)
          dg_e(1:3,0,3) = Dgc(1:3,3,4)                                         
          dg_e(1:3,1:3,1:3) = DgC(1:3,1:3,1:3)
                    
          K_e(1:3,1:3) = - 0.5d0*DgC(4,1:3,1:3)
          K_e(0,1:3) = - 0.5d0*DgC(4,4,1:3)

        end if 
	
      else if ((round .eq. 2) .AND. (par(P_BH1_MASS) .gt. 0.0))  then
           
!---- spining black hole--------------------------------------

        M = par(P_BH1_MASS)            
        a = par(P_BH1_SPIN)                  
        theta = par(P_BH1_SPIN_THETA)      

        x = xold - par(P_BH1_X0)
        y = yold - par(P_BH1_Y0)
        z = zold - par(P_BH1_Z0)      

        if ((abs(x).lt.1E-3) .AND. (abs(y).lt.1E-3)&
     &    .AND. (abs(z).lt.1E-3)) then
      
          x=1E-3
          y=1E-3
          z=1E-3
      
        end if
                    	      
!        r = dsqrt(x**2 + y**2 + z**2)
!        if (r.lt.1E-2) r=1E-2            
!#include "spinKS.inc"             

#include "spinorientedKS.inc"                                           

#include "KS_GH.inc"             

        g_e(0,0) = g_e(0,0) + g00 + 1.0d0
        g_e(0,1) = g_e(0,1) + g01
        g_e(0,2) = g_e(0,2) + g02
        g_e(0,3) = g_e(0,3) + g03                                           
        g_e(1,1) = g_e(1,1) + g11 - 1.0d0
        g_e(1,2) = g_e(1,2) + g12
        g_e(1,3) = g_e(1,3) + g13
        g_e(2,2) = g_e(2,2) + g22 - 1.0d0
        g_e(2,3) = g_e(2,3) + g23
        g_e(3,3) = g_e(3,3) + g33 - 1.0d0

        dg_e(1,0,0) = dg_e(1,0,0) + D100	
        dg_e(1,0,1) = dg_e(1,0,1) + D101
        dg_e(1,0,2) = dg_e(1,0,2) + D102
        dg_e(1,0,3) = dg_e(1,0,3) + D103                  
        dg_e(1,1,1) = dg_e(1,1,1) + D111
        dg_e(1,1,2) = dg_e(1,1,2) + D112
        dg_e(1,1,3) = dg_e(1,1,3) + D113
        dg_e(1,2,2) = dg_e(1,2,2) + D122
        dg_e(1,2,3) = dg_e(1,2,3) + D123
        dg_e(1,3,3) = dg_e(1,3,3) + D133

        dg_e(2,0,0) = dg_e(2,0,0) + D200	        
	dg_e(2,0,1) = dg_e(2,0,1) + D201
        dg_e(2,0,2) = dg_e(2,0,2) + D202
        dg_e(2,0,3) = dg_e(2,0,3) + D203                                    
        dg_e(2,1,1) = dg_e(2,1,1) + D211
        dg_e(2,1,2) = dg_e(2,1,2) + D212
        dg_e(2,1,3) = dg_e(2,1,3) + D213
        dg_e(2,2,2) = dg_e(2,2,2) + D222
        dg_e(2,2,3) = dg_e(2,2,3) + D223
        dg_e(2,3,3) = dg_e(2,3,3) + D233      

        dg_e(3,0,0) = dg_e(3,0,0) + D300        
	dg_e(3,0,1) = dg_e(3,0,1) + D301
        dg_e(3,0,2) = dg_e(3,0,2) + D302
        dg_e(3,0,3) = dg_e(3,0,3) + D303                              
        dg_e(3,1,1) = dg_e(3,1,1) + D311
        dg_e(3,1,2) = dg_e(3,1,2) + D312
        dg_e(3,1,3) = dg_e(3,1,3) + D313
        dg_e(3,2,2) = dg_e(3,2,2) + D322
        dg_e(3,2,3) = dg_e(3,2,3) + D323
        dg_e(3,3,3) = dg_e(3,3,3) + D333

        K_e(1,1) = K_e(1,1) + K11
        K_e(1,2) = K_e(1,2) + K12
        K_e(1,3) = K_e(1,3) + K13
        K_e(2,2) = K_e(2,2) + K22
        K_e(2,3) = K_e(2,3) + K23
        K_e(3,3) = K_e(3,3) + K33
			
      end if	

!#######################################################################
!#################### Z- WAVES #########################################
!#######################################################################

    ! Initialize everything to Minkowski plus a gaussian in Z0
    else if (idtype.eq.88) then 
      g_e(0,0) =-1.0
      g_e(0,1) = 0.0 
      g_e(0,2) = 0.0
      g_e(0,3) = 0.0
      g_e(1,1) = 1.0
      g_e(1,2) = 0.0
      g_e(1,3) = 0.0
      g_e(2,2) = 1.0 
      g_e(2,3) = 0.0
      g_e(3,3) = 1.0

      K_e  = 0.0
      dg_e = 0.0
      H_e  = 0.0
      dH_e = 0.0

      phir_e  = 0.0
      dphir_e = 0.0
      phic_e  = 0.0
      dphic_e = 0.0
      phim_e  = 0.0
      dphim_e = 0.0

      ! setting a gaussian for the Z0
      amp   = par(P_TW_AMP)
      width = par(P_TW_WIDTH)
      ro    = par(P_TW_TIME_SHIFT)

      x = xold
      y = yold
      z = zold
!      r = dsqrt( x**2 + y**2 + z**2 )
      r = sqrt( x**2 + y**2 )
      if (r.lt.1E-3) r=1E-3

      Z0 = amp*exp(-(r-ro)**2/(width**2))
      Z1 = 0.0; Z2 =0.0; Z3 =0.0
!      Z1 = Z0*(-2.0d0*(r-ro)/(width**2))*x/r
!      Z2 = phim_e*(-2.0d0*(r-ro)/(width**2))*y/r
!      Z3 = phim_e*(-2.0d0*(r-ro)/(width**2))*z/r
!      Z0 =+(b1*dphim_e(1) + b2*dphim_e(2) + b3*dphim_e(3))/alp


!#################### sino pues a tomar poc culo #######################
    else 

      write(*,*)'this initial data is not defined!!'
      STOP
    end if

    if (ltrace)  write(*,*)'@@@ now for the initial Ps'

!---------------------------------------------------------------------
!--compute now the initial values of the gauge variables -------------
!---------------------------------------------------------------------
    
#include "GH_initialgauge.inc"

      H_e(0) = H0
      H_e(1) = H1
      H_e(2) = H2
      H_e(3) = H3

      K_e(0,0) = K00
      K_e(0,1) = K01
      K_e(0,2) = K02
      K_e(0,3) = K03
      K_e(1,1) = K11
      K_e(1,2) = K12
      K_e(1,3) = K13
      K_e(2,2) = K22
      K_e(2,3) = K23
      K_e(3,3) = K33    
  

    if (ltrace)  write(*,*)'@@@ finish the Ps succesfully'

  end subroutine

  subroutine gh_flatspace(u0, par)
    use GF
    implicit none
    type(gridfunction), dimension(NU_G)   :: u0
    CCTK_REAL, dimension(:)               :: par

    u0(H_G00)%d = -1.0d0
    u0(H_G01)%d = 0.0d0
    u0(H_G02)%d = 0.0d0
    u0(H_G03)%d = 0.0d0
    u0(H_G11)%d = 1.0d0
    u0(H_G12)%d = 0.0d0
    u0(H_G13)%d = 0.0d0
    u0(H_G22)%d = 1.0d0
    u0(H_G23)%d = 0.0d0
    u0(H_G33)%d = 1.0d0

    u0(H_K00)%d = 0.0d0
    u0(H_K01)%d = 0.0d0
    u0(H_K02)%d = 0.0d0
    u0(H_K03)%d = 0.0d0
    u0(H_K11)%d = 0.0d0
    u0(H_K12)%d = 0.0d0
    u0(H_K13)%d = 0.0d0
    u0(H_K22)%d = 0.0d0
    u0(H_K23)%d = 0.0d0
    u0(H_K33)%d = 0.0d0

    u0(H_D1G00)%d = 0.0d0
    u0(H_D1G01)%d = 0.0d0
    u0(H_D1G02)%d = 0.0d0
    u0(H_D1G03)%d = 0.0d0
    u0(H_D1G11)%d = 0.0d0
    u0(H_D1G12)%d = 0.0d0
    u0(H_D1G13)%d = 0.0d0
    u0(H_D1G22)%d = 0.0d0
    u0(H_D1G23)%d = 0.0d0
    u0(H_D1G33)%d = 0.0d0

    u0(H_D2G00)%d = 0.0d0
    u0(H_D2G01)%d = 0.0d0
    u0(H_D2G02)%d = 0.0d0
    u0(H_D2G03)%d = 0.0d0
    u0(H_D2G11)%d = 0.0d0
    u0(H_D2G12)%d = 0.0d0
    u0(H_D2G13)%d = 0.0d0
    u0(H_D2G22)%d = 0.0d0
    u0(H_D2G23)%d = 0.0d0
    u0(H_D2G33)%d = 0.0d0

    u0(H_D3G00)%d = 0.0d0
    u0(H_D3G01)%d = 0.0d0
    u0(H_D3G02)%d = 0.0d0
    u0(H_D3G03)%d = 0.0d0
    u0(H_D3G11)%d = 0.0d0
    u0(H_D3G12)%d = 0.0d0
    u0(H_D3G13)%d = 0.0d0
    u0(H_D3G22)%d = 0.0d0
    u0(H_D3G23)%d = 0.0d0
    u0(H_D3G33)%d = 0.0d0

    u0(H_PHIM)%d = 0.0d0
    u0(H_PIM)%d = 0.0d0
    u0(H_D1PHIM)%d = 0.0d0
    u0(H_D2PHIM)%d = 0.0d0
    u0(H_D3PHIM)%d = 0.0d0

    u0(H_PHIN)%d = 0.0d0
    u0(H_PIN)%d = 0.0d0
    u0(H_D1PHIN)%d = 0.0d0
    u0(H_D2PHIN)%d = 0.0d0
    u0(H_D3PHIN)%d = 0.0d0

    u0(H_PHIR)%d = 0.0d0
    u0(H_PIR)%d = 0.0d0
    u0(H_D1PHIR)%d = 0.0d0
    u0(H_D2PHIR)%d = 0.0d0
    u0(H_D3PHIR)%d = 0.0d0

    u0(H_PHIC)%d = 0.0d0
    u0(H_PIC)%d = 0.0d0
    u0(H_D1PHIC)%d = 0.0d0
    u0(H_D2PHIC)%d = 0.0d0
    u0(H_D3PHIC)%d = 0.0d0

    return
  end subroutine
end module 



