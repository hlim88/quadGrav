!--------------------------------------------------------------------
!
!  $Id: analysis.f90,v 1.59 2009-09-11 18:45:11 steve Exp $
!
!--------------------------------------------------------------------

#include "cctk.h"

#include "basic.h"
#include "GH_rhs.h"
#include "Dus.h"
#include "GH_analysis.h"

#if defined GH_MHD || defined GH_FLOWER
module mod_gh_analysis
#else
module mod_analysis
#endif
  use params
  use UTILEQS
  use GF
  use MPSI4
  use MRHS
    
implicit none
CONTAINS

  !--------------------------------------------------------------------
  !
  ! GF analysis
  !
  !--------------------------------------------------------------------

#if defined GH_MHD || defined GH_FLOWER
  subroutine Xgh_analysis(u2G, u0G, dxuG, dyuG, dzuG, v, dxv, dyv, dzv, w, par)
#else
  subroutine analysis(u2G, u0G, dxuG, dyuG, dzuG, v, dxv, dyv, dzv, w, par)
#endif
    implicit none
    type(gridfunction), dimension(NU_G)   :: u2G,u0G
    type(gridfunction), dimension(NV_G)   :: v, dxv, dyv, dzv
    type(gridfunction), dimension(NW)     :: w
    type(gridfunction), dimension(NU_G)   :: dxuG,dyuG,dzuG
    CCTK_REAL, dimension(:)               :: par


    call my_exit('analysis : the point_wise_analysis has to be set to 1')    
        
    return
  end subroutine
  
  
  
  
  
  !--------------------------------------------------------------------
  !
  ! Point wise analysis, there is a parameter in hyper called
  !     point_wise_analysis in order to call this routine
  !
  !--------------------------------------------------------------------
  
#if defined GH_MHD || defined GH_FLOWER
  subroutine gh_analysis_pt(u_pt, u0_pt, urk_pt, dxu_pt, dyu_pt, dzu_pt, &
                         v_pt , dxv_pt, dyv_pt, dzv_pt, w_pt, par)
#else
  subroutine analysis_pt(u_pt, u0_pt, urk_pt, dxu_pt, dyu_pt, dzu_pt, &
                         v_pt , dxv_pt, dyv_pt, dzv_pt, w_pt, par)
#endif
    implicit none
    CCTK_REAL, dimension(NU_G), intent(in)     :: dxu_pt,dyu_pt,dzu_pt
    CCTK_REAL, dimension(NV_G), intent(in)     :: dxv_pt,dyv_pt,dzv_pt
    CCTK_REAL, dimension(NU_G), intent(in)     :: u0_pt,u_pt
    CCTK_REAL, dimension(NU_G), intent(inout)  :: urk_pt    
    CCTK_REAL, dimension(NV_G), intent(inout)  :: v_pt
    CCTK_REAL, dimension(NW), intent(inout)    :: w_pt
    CCTK_REAL, dimension(NPAR),intent(inout)   :: par
    CCTK_REAL :: x,y,z,r,rcyl,rho,pi,rcoor

    CCTK_REAL                             :: alp, alp2, deth, strK, sKK, Rnn
    CCTK_REAL                             :: Sdds11, Sdds22, Sdds23, Sdds33
    CCTK_REAL, dimension(I_1:I_3)         :: beta, trD, trE, trDu, trEu
    CCTK_REAL, dimension(I_1:I_3)         :: Dnn
    CCTK_REAL, dimension(I_0:I_3)         :: tu, trG, Zd, ttD, xiu
    CCTK_REAL, dimension(I_0:I_3)         :: Duphim, Duphin, Duphir, Duphic
    CCTK_REAL, dimension(IS3_11:IS3_33)   :: h, huu, eK, Suu, Gndd, sK, Sdd, Ri
    CCTK_REAL, dimension(IN3_11:IN3_33)   :: eKud, Sud, Dn, sKud, Cr, Ci
    CCTK_REAL, dimension(IN4_00:IN4_33)   :: tD
    CCTK_REAL, dimension(IS4_00:IS4_33)   :: guu, D0
    CCTK_REAL, dimension(IS3_111:IS3_333) :: Dduu, Dudd
    CCTK_REAL, dimension(IN3_111:IN3_333) :: Dddu
    CCTK_REAL, dimension(IS4_000:IS4_333) :: Gudd, Gddd 

    ! for the momentum constraint
    CCTK_REAL, dimension(I_1:I_3)         :: trGu
    CCTK_REAL, dimension(IS3_11:IS3_33)   :: Xt
    CCTK_REAL, dimension(IS3_111:IS3_333) :: dXt, deK, dceK
    CCTK_REAL, dimension(IN4_00:IN4_33)   :: dbeta
    CCTK_REAL  :: S1, S2, S3, momx_const, momy_const, momz_const

    CCTK_REAL :: etrK, trQ, Q, trB, tau, atau, aS1, aS2, aS3
    CCTK_REAL :: sm, energy_density, noether_density
    CCTK_REAL :: dphidphir, dphidphic, dphidphim, dphidphin
    CCTK_REAL :: dphim0, dphin0, dphir0, dphic0    
    CCTK_REAL :: H0e, H1e, H2e, H3e    

    CCTK_REAL :: energy_const, allconst_pt     
    CCTK_REAL :: order1_const_pt, order2_const_pt, physZ_const_pt
    CCTK_REAL :: rad_speed_pt, rad_exp_pt, ra, rb, rc
    CCTK_REAL :: nr1, nr2, nr3, nrnorm, nru1, nru2, nru3
    CCTK_REAL :: dx_nr1,dx_nr2,dx_nr3,dy_nr1,dy_nr2,dy_nr3,dz_nr1,dz_nr2,dz_nr3
    CCTK_REAL :: term1, term2, term3, term4

    CCTK_REAL :: z4rhs0_pt, z4rhs1_pt, z4rhs2_pt, z4rhs3_pt, trT
    
    CCTK_INT  :: gauge_type, psi4_analysis, li, lj, lk, lm, temp
    CCTK_INT  :: la, lb, lc, ld, le, lf, lg, lh    
    
    ! for the psi4
    CCTK_REAL, dimension(3) 		   :: dalp, M
    CCTK_REAL, dimension(3,3)              :: KtrK, K, h_uu, Ricci3, JA
    CCTK_REAL, dimension(4) 		   :: Tl,Tn,Tm_R,Tm_I    
    CCTK_REAL, dimension(4,4)              :: g, g_uu
    CCTK_REAL, dimension(4,4,4)            :: dg, Chd
    CCTK_REAL, dimension(4,4,4,4)          :: D_dg, R4
    CCTK_REAL, dimension(3,3,3,3)          :: R3
    CCTK_REAL                              :: Psi4R, Psi4I, rt, R3scalar, R2scalar
    CCTK_REAL                              :: ADMmass,BONDImass,curvature,xinorm,gur
    CCTK_REAL                              :: lin_momx,lin_momy,lin_momz
    CCTK_REAL                              :: ang_momx,ang_momy,ang_momz,kang_momz
    CCTK_REAL                              :: costheta,sintheta,cosphi,sinphi,Y2p2,Y20
    CCTK_REAL                              :: bhdist, bhsize, Exx, Eyy, Ezz,absm1,absm2
    ! Save code from lapse getting too small
    !    (not wanted when doing a critical search):
    logical, parameter :: rescue   = .true.
   
    integer nummask
    external masks_return_num
    integer masks_return_num
 
!---------------------------------------------------------    
                    
    x = w_pt(H_XPHYS)
    y = w_pt(H_YPHYS)
    z = w_pt(H_ZPHYS)
        
    gauge_type = par(P_GAUGE_TYPE)    
    psi4_analysis = par(P_PSI4_ANALYSIS)
    sm = par(P_SF_MASS)        
    
    pi = 3.141592653589793d0
!    r = sqrt(x*x + y*y + z*z - a*a)
    r = sqrt(x*x + y*y + z*z)
    if (r .LT. 1E-4) r=1E-4

    rcyl = sqrt(x*x + y*y)
    if (rcyl .LT. 1E-4) rcyl=1E-4

!-----the first order constraints------------------------------------------

    order1_const_pt = &
          + (D100 - dx_g00)**2 + (D200 - dy_g00)**2 + (D300 - dz_g00)**2 & 
          + (D101 - dx_g01)**2 + (D201 - dy_g01)**2 + (D301 - dz_g01)**2 & 
          + (D102 - dx_g02)**2 + (D202 - dy_g02)**2 + (D302 - dz_g02)**2 & 
          + (D103 - dx_g03)**2 + (D203 - dy_g03)**2 + (D303 - dz_g03)**2 & 
          + (D111 - dx_g11)**2 + (D211 - dy_g11)**2 + (D311 - dz_g11)**2 & 
          + (D113 - dx_g13)**2 + (D213 - dy_g13)**2 + (D313 - dz_g13)**2 & 
          + (D122 - dx_g22)**2 + (D222 - dy_g22)**2 + (D322 - dz_g22)**2 & 
          + (D123 - dx_g23)**2 + (D223 - dy_g23)**2 + (D323 - dz_g23)**2 & 
          + (D133 - dx_g33)**2 + (D233 - dy_g33)**2 + (D333 - dz_g33)**2     
    order1_const_pt = sqrt(order1_const_pt/30.0)
  
    order2_const_pt = &
          + (dx_D200 - dy_D100)**2 + (dx_D300 - dz_D100)**2 &
          + (dy_D300 - dz_D200)**2 &
          + (dx_D201 - dy_D101)**2 + (dx_D301 - dz_D101)**2 &
          + (dy_D301 - dz_D201)**2 &
          + (dx_D202 - dy_D102)**2 + (dx_D302 - dz_D102)**2 &
          + (dy_D302 - dz_D202)**2 &
          + (dx_D203 - dy_D103)**2 + (dx_D303 - dz_D103)**2 &
          + (dy_D303 - dz_D203)**2 &
          + (dx_D211 - dy_D111)**2 + (dx_D311 - dz_D111)**2 &
          + (dy_D311 - dz_D211)**2 &
          + (dx_D212 - dy_D112)**2 + (dx_D312 - dz_D112)**2 &
          + (dy_D312 - dz_D212)**2 &
          + (dx_D213 - dy_D113)**2 + (dx_D313 - dz_D113)**2 &
          + (dy_D313 - dz_D213)**2 &
          + (dx_D222 - dy_D122)**2 + (dx_D322 - dz_D122)**2 &
          + (dy_D322 - dz_D222)**2 &
          + (dx_D223 - dy_D123)**2 + (dx_D323 - dz_D123)**2 &
          + (dy_D323 - dz_D223)**2 &
          + (dx_D233 - dy_D133)**2 + (dx_D333 - dz_D133)**2 &
          + (dy_D333 - dz_D233)**2 
    order2_const_pt = sqrt(order2_const_pt/30.0)


!-----the 3+1 quantities----------------------------------------------------
       
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
       if(deth.lt.0.1) deth = 0.1
       huu11 = huu11/deth
       huu12 = huu12/deth
       huu13 = huu13/deth
       huu22 = huu22/deth
       huu23 = huu23/deth
       huu33 = huu33/deth
                     
       b1 = g01*huu11 + g02*huu12 + g03*huu13
       b2 = g01*huu12 + g02*huu22 + g03*huu23
       b3 = g01*huu13 + g02*huu23 + g03*huu33
       alp2 = -g00 + b1**2*h11 + b2**2*h22 + 2*(b1*(b2*h12 + b3*h13) &
     &      + b2*b3*h23) + b3**2*h33
       if(alp2.lt.1E-4) alp2=1E-4
       alp = sqrt(alp2)
      
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
       
       eK11 = -(-alp*K11 - D101 - D101 &
      &    +b1*(D111 + D111) + b2*(D112 + D112) &
      &    + b3*(D113 + D113))/(2.0*alp)
       eK12 = -(-alp*K12 - D102 - D201 &
      &    +b1*(D112 + D211) + b2*(D122 + D212) &
      &    + b3*(D123 + D213))/(2.0*alp)
       eK13 = -(-alp*K13 - D103 - D301 &
      &    +b1*(D113 + D311) + b2*(D123 + D312) &
      &    + b3*(D133 + D313))/(2.0*alp)
       eK22 = -(-alp*K22 - D202 - D202 &
      &    +b1*(D212 + D212) + b2*(D222 + D222) &
      &    + b3*(D223 + D223))/(2.0*alp)
       eK23 = -(-alp*K23 - D203 - D302 &
      &    +b1*(D213 + D312) + b2*(D223 + D322) &
      &    + b3*(D233 + D323))/(2.0*alp)
       eK33 = -(-alp*K33 - D303 - D303 &
      &    +b1*(D313 + D313) + b2*(D323 + D323) &
      &    + b3*(D333 + D333))/(2.0*alp)


!-----a lot of stuff for the energy constraint---------------------------
      
       Dddu111 = huu11*D111 + huu12*D112 + huu13*D113
       Dddu112 = huu12*D111 + huu22*D112 + huu23*D113
       Dddu113 = huu13*D111 + huu23*D112 + huu33*D113
       Dddu121 = huu11*D112 + huu12*D122 + huu13*D123
       Dddu122 = huu12*D112 + huu22*D122 + huu23*D123
       Dddu123 = huu13*D112 + huu23*D122 + huu33*D123
       Dddu131 = huu11*D113 + huu12*D123 + huu13*D133
       Dddu132 = huu12*D113 + huu22*D123 + huu23*D133
       Dddu133 = huu13*D113 + huu23*D123 + huu33*D133
       Dddu211 = huu11*D211 + huu12*D212 + huu13*D213
       Dddu212 = huu12*D211 + huu22*D212 + huu23*D213
       Dddu213 = huu13*D211 + huu23*D212 + huu33*D213
       Dddu221 = huu11*D212 + huu12*D222 + huu13*D223
       Dddu222 = huu12*D212 + huu22*D222 + huu23*D223
       Dddu223 = huu13*D212 + huu23*D222 + huu33*D223
       Dddu231 = huu11*D213 + huu12*D223 + huu13*D233
       Dddu232 = huu12*D213 + huu22*D223 + huu23*D233
       Dddu233 = huu13*D213 + huu23*D223 + huu33*D233
       Dddu311 = huu11*D311 + huu12*D312 + huu13*D313
       Dddu312 = huu12*D311 + huu22*D312 + huu23*D313
       Dddu313 = huu13*D311 + huu23*D312 + huu33*D313
       Dddu321 = huu11*D312 + huu12*D322 + huu13*D323
       Dddu322 = huu12*D312 + huu22*D322 + huu23*D323
       Dddu323 = huu13*D312 + huu23*D322 + huu33*D323
       Dddu331 = huu11*D313 + huu12*D323 + huu13*D333
       Dddu332 = huu12*D313 + huu22*D323 + huu23*D333
       Dddu333 = huu13*D313 + huu23*D323 + huu33*D333
       
       Dduu111 = huu11**2*D111 + huu12**2*D122 + 2*(huu11*(huu12*D112 +&
     & huu13*D113) + huu12*huu13*D123) + huu13**2*D133
       Dduu112 = huu12**2*D112 + huu11* (huu12*D111 + huu22*D112 +&
     & huu23*D113) + huu12*(huu13*D113 + huu22*D122 + huu23*D123) +&
     & huu13* (huu22*D123 + huu23*D133)
       Dduu113 = huu13**2*D113 + huu11* (huu23*D112 + huu33*D113) +&
     & huu12*(huu23*D122 + huu33*D123) + huu13*(huu11*D111 + huu12*D112&
     & + huu23*D123 + huu33*D133)
       Dduu122 = huu12**2*D111 + huu22**2*D122 + 2*(huu12*(huu22*D112 +&
     & huu23*D113) + huu22*huu23*D123) + huu23**2*D133
       Dduu123 = huu13*(huu12*D111 + huu22*D112) + huu23**2*D123 +&
     & huu33* (huu12*D113 + huu22*D123) + huu23*(huu12*D112 + huu13*D113&
     & + huu22*D122 + huu33*D133)
       Dduu133 = huu13**2*D111 + huu23**2*D122 + 2*(huu13*(huu23*D112 +&
     & huu33*D113) + huu23*huu33*D123) + huu33**2*D133
       Dduu211 = huu11**2*D211 + huu12**2*D222 + 2*(huu11*(huu12*D212 +&
     & huu13*D213) + huu12*huu13*D223) + huu13**2*D233
       Dduu212 = huu12**2*D212 + huu11* (huu12*D211 + huu22*D212 +&
     & huu23*D213) + huu12*(huu13*D213 + huu22*D222 + huu23*D223) +&
     & huu13* (huu22*D223 + huu23*D233)
       Dduu213 = huu13**2*D213 + huu11* (huu23*D212 + huu33*D213) +&
     & huu12*(huu23*D222 + huu33*D223) + huu13*(huu11*D211 + huu12*D212&
     & + huu23*D223 + huu33*D233)
       Dduu222 = huu12**2*D211 + huu22**2*D222 + 2*(huu12*(huu22*D212 +&
     & huu23*D213) + huu22*huu23*D223) + huu23**2*D233
       Dduu223 = huu13*(huu12*D211 + huu22*D212) + huu23**2*D223 +&
     & huu33* (huu12*D213 + huu22*D223) + huu23*(huu12*D212 + huu13*D213&
     & + huu22*D222 + huu33*D233)
       Dduu233 = huu13**2*D211 + huu23**2*D222 + 2*(huu13*(huu23*D212 +&
     & huu33*D213) + huu23*huu33*D223) + huu33**2*D233
       Dduu311 = huu11**2*D311 + huu12**2*D322 + 2*(huu11*(huu12*D312 +&
     & huu13*D313) + huu12*huu13*D323) + huu13**2*D333
       Dduu312 = huu12**2*D312 + huu11* (huu12*D311 + huu22*D312 +&
     & huu23*D313) + huu12*(huu13*D313 + huu22*D322 + huu23*D323) +&
     & huu13* (huu22*D323 + huu23*D333)
       Dduu313 = huu13**2*D313 + huu11* (huu23*D312 + huu33*D313) +&
     & huu12*(huu23*D322 + huu33*D323) + huu13*(huu11*D311 + huu12*D312&
     & + huu23*D323 + huu33*D333)
       Dduu322 = huu12**2*D311 + huu22**2*D322 + 2*(huu12*(huu22*D312 +&
     & huu23*D313) + huu22*huu23*D323) + huu23**2*D333
       Dduu323 = huu13*(huu12*D311 + huu22*D312) + huu23**2*D323 +&
     & huu33* (huu12*D313 + huu22*D323) + huu23*(huu12*D312 + huu13*D313&
     & + huu22*D322 + huu33*D333)
       Dduu333 = huu13**2*D311 + huu23**2*D322 + 2*(huu13*(huu23*D312 +&
     & huu33*D313) + huu23*huu33*D323) + huu33**2*D333
       Dudd111 = huu11*D111 + huu12*D211 + huu13*D311
       Dudd112 = huu11*D112 + huu12*D212 + huu13*D312
       Dudd113 = huu11*D113 + huu12*D213 + huu13*D313
       Dudd122 = huu11*D122 + huu12*D222 + huu13*D322
       Dudd123 = huu11*D123 + huu12*D223 + huu13*D323
       Dudd133 = huu11*D133 + huu12*D233 + huu13*D333
       Dudd211 = huu12*D111 + huu22*D211 + huu23*D311
       Dudd212 = huu12*D112 + huu22*D212 + huu23*D312
       Dudd213 = huu12*D113 + huu22*D213 + huu23*D313
       Dudd222 = huu12*D122 + huu22*D222 + huu23*D322
       Dudd223 = huu12*D123 + huu22*D223 + huu23*D323
       Dudd233 = huu12*D133 + huu22*D233 + huu23*D333
       Dudd311 = huu13*D111 + huu23*D211 + huu33*D311
       Dudd312 = huu13*D112 + huu23*D212 + huu33*D312
       Dudd313 = huu13*D113 + huu23*D213 + huu33*D313
       Dudd322 = huu13*D122 + huu23*D222 + huu33*D322
       Dudd323 = huu13*D123 + huu23*D223 + huu33*D323
       Dudd333 = huu13*D133 + huu23*D233 + huu33*D333
       trD1 = Dddu111 + Dddu122 + Dddu133
       trD2 = Dddu211 + Dddu222 + Dddu233
       trD3 = Dddu311 + Dddu322 + Dddu333
       trE1 = Dudd111 + Dudd212 + Dudd313
       trE2 = Dudd112 + Dudd222 + Dudd323
       trE3 = Dudd113 + Dudd223 + Dudd333
       trDu1 = huu11*trD1 + huu12*trD2 + huu13*trD3
       trDu2 = huu12*trD1 + huu22*trD2 + huu23*trD3
       trDu3 = huu13*trD1 + huu23*trD2 + huu33*trD3
       trEu1 = huu11*trE1 + huu12*trE2 + huu13*trE3
       trEu2 = huu12*trE1 + huu22*trE2 + huu23*trE3
       trEu3 = huu13*trE1 + huu23*trE2 + huu33*trE3

       !with the first index up                                             
       G111 = Dddu111 - Dudd111/2.
       G112 = (Dddu121 + Dddu211 - Dudd112)/2.
       G113 = (Dddu131 + Dddu311 - Dudd113)/2.
       G122 = Dddu221 - Dudd122/2.
       G123 = (Dddu231 + Dddu321 - Dudd123)/2.
       G133 = Dddu331 - Dudd133/2.
       G211 = Dddu112 - Dudd211/2.
       G212 = (Dddu122 + Dddu212 - Dudd212)/2.
       G213 = (Dddu132 + Dddu312 - Dudd213)/2.
       G222 = Dddu222 - Dudd222/2.
       G223 = (Dddu232 + Dddu322 - Dudd223)/2.
       G233 = Dddu332 - Dudd233/2.
       G311 = Dddu113 - Dudd311/2.
       G312 = (Dddu123 + Dddu213 - Dudd312)/2.
       G313 = (Dddu133 + Dddu313 - Dudd313)/2.
       G322 = Dddu223 - Dudd322/2.
       G323 = (Dddu233 + Dddu323 - Dudd323)/2.
       G333 = Dddu333 - Dudd333/2.
                    
       eKud11 = huu11*eK11 + huu12*eK12 + huu13*eK13
       eKud12 = huu11*eK12 + huu12*eK22 + huu13*eK23
       eKud13 = huu11*eK13 + huu12*eK23 + huu13*eK33
       eKud21 = huu12*eK11 + huu22*eK12 + huu23*eK13
       eKud22 = huu12*eK12 + huu22*eK22 + huu23*eK23
       eKud23 = huu12*eK13 + huu22*eK23 + huu23*eK33
       eKud31 = huu13*eK11 + huu23*eK12 + huu33*eK13
       eKud32 = huu13*eK12 + huu23*eK22 + huu33*eK23
       eKud33 = huu13*eK13 + huu23*eK23 + huu33*eK33
       etrK = eKud11 + eKud22 + eKud33       
       
!       energy_const = (dx_D122 - dx_D212 - dy_D112 + dy_D211)*huu12**2 -&
!     & dz_D113*huu13**2 + (dx_D133 - dx_D313 + dz_D311)*huu13**2 +&
!     & ((-2*dx_D123 + dz_D112)*huu11 - dz_D211*huu11 + (-dx_D223 +&
!     & 2*dy_D213 + dz_D122)*huu12 - dz_D212*huu12 + (-dx_D323 + dy_D133&
!     & + 2*dz_D312)*huu13)*huu23 - dz_D223*huu23**2 + (dy_D233 -&
!     & dy_D323 + dz_D322)*huu23**2 + huu13*(-(dx_D312*huu12) -&
!     & 2*dy_D213*huu22 + (dx_D223 - dx_D322 + dy_D123 + dy_D312 +&
!     & dz_D212)*huu22 + (dx_D233 - dy_D313)*huu23 - dz_D123*huu23) -&
!     & huu13*(dz_D112*huu12 + dz_D122*huu22 + dz_D213*huu23) -&
!     & dz_D322*huu22*huu33 + ((-dx_D133 + dz_D113)*huu11 + (dy_D313 -&
!     & 2*dz_D312)*huu12 + (-dy_D233 + dy_D323 + dz_D223)*huu22)*huu33 +&
!     & huu11*((-dx_D122 + dx_D212 + dy_D112)*huu22 + (dx_D213 + dx_D312&
!     & + dy_D113)*huu23 + dx_D313*huu33) + huu12*(-(dy_D113*huu13) +&
!     & (2*dx_D123 - dx_D213 + dy_D311 + dz_D211)* huu13 + (dx_D322 -&
!     & dy_D123)*huu23 - dy_D312*huu23 - dy_D133*huu33 + (-dx_D233 +&
!     & dx_D323 + dz_D123 + dz_D213)*huu33) - huu11*(dy_D211*huu22 +&
!     & dy_D311*huu23 + dz_D311*huu33) - eKud11**2 - eKud22**2 -&
!     & 2*(eKud12*eKud21 + eKud13*eKud31 + eKud23*eKud32) - eKud33**2 -&
!     & trE1*trEu1 - trE2*trEu2 - trE3*trEu3 + etrK**2 + trEu1*trD1 -&
!     & (trD1*trDu1)/4. + trD2*(trEu2 - trDu2/4.) + trD3*(trEu3 -&
!     & trDu3/4.) - Dddu111*Dduu111 - Dddu121*Dduu112 - Dddu131*Dduu113 -&
!     & Dddu221*Dduu122 - Dddu231*Dduu123 - Dddu331*Dduu133 -&
!     & Dddu112*Dduu211 - Dddu122*Dduu212 - Dddu132*Dduu213 -&
!     & Dddu222*Dduu222 - Dddu232*Dduu223 - Dddu332*Dduu233 -&
!     & Dddu113*Dduu311 - Dddu123*Dduu312 - Dddu133*Dduu313 -&
!     & Dddu223*Dduu322 - Dddu233*Dduu323 - Dddu333*Dduu333 +&
!     & (G111*Dduu111 + G122*Dduu122 + G133*Dduu133 + G211*Dduu211 +&
!     & G222*Dduu222 + G233*Dduu233 + G311*Dduu311 + G322*Dduu322 +&
!     & G333*Dduu333)/ 2. + Dduu111*Dudd111 + Dduu112*(G112 - Dddu211 +&
!     & 2*Dudd112) + Dduu113*(G113 - Dddu311 + 2*Dudd113) +&
!     & Dduu122*Dudd122 + Dduu123* (G123 - Dddu321 + 2*Dudd123) +&
!     & Dduu133*Dudd133 + Dduu211*Dudd211 + Dduu212*(G212 - Dddu212 +&
!     & 2*Dudd212) + Dduu213*(G213 - Dddu312 + 2*Dudd213) +&
!     & Dduu222*Dudd222 + Dduu223* (G223 - Dddu322 + 2*Dudd223) +&
!     & Dduu233*Dudd233 + Dduu311*Dudd311 + Dduu312*(G312 - Dddu213 +&
!     & 2*Dudd312) + Dduu313*(G313 - Dddu313 + 2*Dudd313) +&
!     & Dduu322*Dudd322 + Dduu323* (G323 - Dddu323 + 2*Dudd323) +&
!     & Dduu333*Dudd333

       Ri11 = -G212**2 - 2*G213*G312 - G313**2 + G111*(G212 + G313) +&
     & G211*(-G112 + G222 + G323) + G311*(-G113 + G223 + G333) +&
     & (-dx_D122/2. + dy_D112)*huu22 + (-dx_D123 + dy_D113 +&
     & dz_D112)*huu23 + (-dx_D133/2. + dz_D113)*huu33 +&
     & (-(dx_D311*huu13) - dy_D211*huu22 - (dy_D311 + dz_D211)*huu23 -&
     & dz_D311*huu33)/2. - trE1*Dddu111 - trE2*Dddu112 - trE3*Dddu113 +&
     & Dddu112*Dddu121 + Dddu113*Dddu131 + Dddu123*Dddu132 + ((-dx_D211&
     & + dy_D111)*huu12 + dz_D111*huu13 + Dddu111**2 + Dddu122**2 +&
     & Dddu133**2 + trE1*Dudd111 + trE2*Dudd211 + trE3*Dudd311)/2.+0.d0
       Ri12 = -(G122*G211) - G123*G311 + G112*(G212 + G313) - G213*G322&
     & + (G212 - G313)*G323 + G312*G333 + ((dx_D211 - dy_D111)*huu11 +&
     & dy_D122*huu22)/4. - (dy_D133/4. + dz_D312/2.)*huu33 +&
     & (-(dx_D222*huu22) - dx_D233*huu33)/4. + (-(dy_D112*huu12) -&
     & (dx_D312 + dy_D113)*huu13 - (dx_D223 + dy_D312)*huu23 -&
     & trE1*Dddu121 - trE2*Dddu212 - trE3*(Dddu123 + Dddu213))/2. +&
     & ((dx_D122 - dx_D212 + dy_D211)*huu12 + (dx_D123 + dz_D211)*huu13&
     & + (dy_D213 + dz_D122)*huu23 + (dz_D123 + dz_D213)*huu33 -&
     & trE2*Dddu122 - trE1*Dddu211 + Dddu111*Dddu211 + Dddu121*Dddu212 +&
     & Dddu131*Dddu213 + Dddu112*Dddu221 + Dddu122*Dddu222 +&
     & Dddu132*Dddu223 + Dddu113*Dddu231 + Dddu123*Dddu232 +&
     & Dddu133*Dddu233 + trE1*Dudd112 + trE2*Dudd212 + trE3*Dudd312)/ 2.+0.d0
       Ri13 = -(G123*G211) + G213*G222 - G133*G311 - G233*G312 +&
     & G223*(-G212 + G313) + G113*(G212 + G313) - (dy_D213/2. + (dx_D322&
     & + dz_D122)/4.)*huu22 - (dx_D333*huu33)/4. + ((dx_D311 -&
     & dz_D111)*huu11 + dz_D133*huu33)/4. + (-(dz_D112*huu12) - (dx_D313&
     & + dz_D113)*huu13 - (dx_D323 + dz_D213)*huu23 - trE1*Dddu131 -&
     & trE2*Dddu312 - trE3*(Dddu133 + Dddu313))/2. + ((dx_D123 -&
     & dx_D213 + dy_D311)*huu12 + (dx_D133 + dz_D311)*huu13 + (dy_D123 +&
     & dy_D312)*huu22 + (dy_D133 + dz_D312)*huu23 - trE2*Dddu132 -&
     & trE1*Dddu311 + Dddu111*Dddu311 + Dddu121*Dddu312 +&
     & Dddu131*Dddu313 + Dddu112*Dddu321 + Dddu122*Dddu322 +&
     & Dddu132*Dddu323 + Dddu113*Dddu331 + Dddu123*Dddu332 +&
     & Dddu133*Dddu333 + trE1*Dudd113 + trE2*Dudd213 + trE3*Dudd313)/ 2.+0.d0
       Ri22 = -G112**2 - 2*G123*G312 + G122*(G111 - G212 + G313) -&
     & G323**2 + G222*(G112 + G323) + G322*(G113 - G223 + G333) +&
     & (-dx_D122/2. + dx_D212)*huu11 - dy_D213*huu13 + (dx_D223 -&
     & dx_D322/2. + dz_D212)*huu13 + (-dy_D233/2. + dz_D223)*huu33 +&
     & (-(dy_D122*huu12) - dz_D122*huu13 - dy_D322*huu23 -&
     & dz_D322*huu33)/2. - trE1*Dddu221 + Dddu212*Dddu221 - trE2*Dddu222&
     & - trE3*Dddu223 + Dddu213*Dddu231 + Dddu223*Dddu232 +&
     & (-(dy_D211*huu11) + dx_D222*huu12 + dz_D222*huu23 + Dddu211**2 +&
     & Dddu222**2 + Dddu233**2 + trE1*Dudd122 + trE2*Dudd222 +&
     & trE3*Dudd322)/ 2.+0.d0
       Ri23 = G111*G123 - G122*G213 - G133*G312 - G233*G322 +&
     & G113*(-G112 + G323) + G223*(G112 + G323) + (-(dz_D211*huu11) -&
     & dz_D222*huu22 - dy_D333*huu33)/4. + (-(dy_D311*huu11) +&
     & dy_D322*huu22 + dz_D233*huu33)/4. + (-((dy_D123 + dz_D212)*huu12)&
     & - (dy_D313 + dz_D123)*huu13 - (dy_D323 + dz_D223)*huu23 -&
     & trE1*(Dddu231 + Dddu321) - trE2*Dddu322 - trE3*(Dddu233 +&
     & Dddu323))/2. + ((-dx_D123 + dx_D213 + dx_D312)*huu11 + (dx_D322 +&
     & dy_D213)*huu12 + (dx_D233 + dz_D312)*huu13 + (dy_D233 +&
     & dz_D322)*huu23 + Dddu211*Dddu311 + Dddu221*Dddu312 +&
     & Dddu231*Dddu313 + Dddu212*Dddu321 + Dddu222*Dddu322 +&
     & Dddu232*(-trE2 + Dddu323) + Dddu213*Dddu331 + Dddu223*Dddu332 +&
     & Dddu233*Dddu333 + trE1*Dudd123 + trE2*Dudd223 + trE3*Dudd323)/2.+0.d0
       Ri33 = -G113**2 - 2*G123*G213 - G223**2 + G133*(G111 + G212 -&
     & G313) + G233*(G112 + G222 - G323) + (G113 + G223)*G333 +&
     & (-dx_D133/2. + dx_D313)*huu11 + (-dx_D233/2. + dx_D323 +&
     & dy_D313)*huu12 - (dy_D133/2. + dz_D312)*huu12 + (-dy_D233/2. +&
     & dy_D323)*huu22 + (-(dz_D133*huu13) - dz_D322*huu22 -&
     & dz_D233*huu23)/2. + Dddu312*Dddu321 - trE1*Dddu331 +&
     & Dddu313*Dddu331 - trE2*Dddu332 + Dddu323*Dddu332 - trE3*Dddu333 +&
     & (-(dz_D311*huu11) + dx_D333*huu13 + dy_D333*huu23 + Dddu311**2 +&
     & Dddu322**2 + Dddu333**2 + trE1*Dudd133 + trE2*Dudd233 +&
     & trE3*Dudd333)/2.+0.d0

       energy_const = -eKud11**2 - eKud22**2 - 2*(eKud12*eKud21 +&
     & eKud13*eKud31 + eKud23*eKud32) - eKud33**2 + huu11*Ri11 +&
     & huu22*Ri22 + 2*(huu12*Ri12 + huu13*Ri13 + huu23*Ri23) +&
     & huu33*Ri33 + etrK**2+0.d0

!-----the energy density with the sqrt(h) from the volume factor----------     

       tu0 = 1.0/alp; tu1 = -b1/alp; tu2 = -b2/alp; tu3 = -b3/alp 
       tau = tu0*tu0*T00 + tu1*tu1*T11 + tu2*tu2*T22 + tu3*tu3*T33 & 
     &      + 2.d0*(tu0*tu1*T01 + tu0*tu2*T02 + tu0*tu3*T03 &
     &            + tu1*tu2*T12 + tu1*tu3*T13 + tu2*tu3*T23)

       ! comment this two lines for testing
       energy_const = energy_const - 16.d0*pi*tau
       energy_density = tau*sqrt(deth)

       ! and uncomment this ones
!       energy_density = 16.d0*pi*tau

!       if (r .LT. 1E-2) then
!          print*, "in constraints"
!          print*, "x,y,z",x,y,z
!          print*, "tau=", tau
!          print*, "T11,T12,T13=", T01
!          print*, "energy_const=", energy_const     
!          print*, "16*pi*rhoADM", 16.d0*pi*tau
!       end if

!----- or the momentum constraint ---------------------------------------     

       D000 = -alp*K00 + b1*D100 + b2*D200 + b3*D300 
       D001 = -alp*K01 + b1*D101 + b2*D201 + b3*D301
       D002 = -alp*K02 + b1*D102 + b2*D202 + b3*D302 
       D003 = -alp*K03 + b1*D103 + b2*D203 + b3*D303 
       D011 = -alp*K11 + b1*D111 + b2*D211 + b3*D311 
       D012 = -alp*K12 + b1*D112 + b2*D212 + b3*D312 
       D013 = -alp*K13 + b1*D113 + b2*D213 + b3*D313 
       D022 = -alp*K22 + b1*D122 + b2*D222 + b3*D322 
       D023 = -alp*K23 + b1*D123 + b2*D223 + b3*D323
       D033 = -alp*K33 + b1*D133 + b2*D233 + b3*D333

       tD00 = tu0*D000 + tu1*D001 + tu2*D002 + tu3*D003+0.d0
       tD01 = tu0*D001 + tu1*D011 + tu2*D012 + tu3*D013+0.d0
       tD02 = tu0*D002 + tu1*D012 + tu2*D022 + tu3*D023+0.d0
       tD03 = tu0*D003 + tu1*D013 + tu2*D023 + tu3*D033+0.d0
       tD10 = tu0*D100 + tu1*D101 + tu2*D102 + tu3*D103+0.d0
       tD11 = tu0*D101 + tu1*D111 + tu2*D112 + tu3*D113+0.d0
       tD12 = tu0*D102 + tu1*D112 + tu2*D122 + tu3*D123+0.d0
       tD13 = tu0*D103 + tu1*D113 + tu2*D123 + tu3*D133+0.d0
       tD20 = tu0*D200 + tu1*D201 + tu2*D202 + tu3*D203+0.d0
       tD21 = tu0*D201 + tu1*D211 + tu2*D212 + tu3*D213+0.d0
       tD22 = tu0*D202 + tu1*D212 + tu2*D222 + tu3*D223+0.d0
       tD23 = tu0*D203 + tu1*D213 + tu2*D223 + tu3*D233+0.d0
       tD30 = tu0*D300 + tu1*D301 + tu2*D302 + tu3*D303+0.d0
       tD31 = tu0*D301 + tu1*D311 + tu2*D312 + tu3*D313+0.d0
       tD32 = tu0*D302 + tu1*D312 + tu2*D322 + tu3*D323+0.d0
       tD33 = tu0*D303 + tu1*D313 + tu2*D323 + tu3*D333+0.d0
       ttD0 = tD00*tu0 + tD01*tu1 + tD02*tu2 + tD03*tu3+0.d0
       ttD1 = tD10*tu0 + tD11*tu1 + tD12*tu2 + tD13*tu3+0.d0
       ttD2 = tD20*tu0 + tD21*tu1 + tD22*tu2 + tD23*tu3+0.d0
       ttD3 = tD30*tu0 + tD31*tu1 + tD32*tu2 + tD33*tu3+0.d0

       d1alp = -(alp*ttD1)/2.+0.d0
       d2alp = -(alp*ttD2)/2.+0.d0
       d3alp = -(alp*ttD3)/2.+0.d0

       d1b1 = alp*(huu11*tD11 + huu12*tD12 + huu13*tD13)+0.d0
       d1b2 = alp*(huu12*tD11 + huu22*tD12 + huu23*tD13)+0.d0
       d1b3 = alp*(huu13*tD11 + huu23*tD12 + huu33*tD13)+0.d0
       d2b1 = alp*(huu11*tD21 + huu12*tD22 + huu13*tD23)+0.d0
       d2b2 = alp*(huu12*tD21 + huu22*tD22 + huu23*tD23)+0.d0
       d2b3 = alp*(huu13*tD21 + huu23*tD22 + huu33*tD23)+0.d0
       d3b1 = alp*(huu11*tD31 + huu12*tD32 + huu13*tD33)+0.d0
       d3b2 = alp*(huu12*tD31 + huu22*tD32 + huu23*tD33)+0.d0
       d3b3 = alp*(huu13*tD31 + huu23*tD32 + huu33*tD33)+0.d0

       !--- some auxiliary fields ---------------------------------          
       X11 = 2*D101 - 2*(b1*D111 + b2*D112 + b3*D113)+0.d0
       X12 = D102 + D201 - b1*D112 - b2*D122 &
     &     - b3*D123 - b1*D211 &
     &     - b2*D212 - b3*D213+0.d0
       X13 = D103 + D301 - b1*D113 &
     &     - b2*D123 - b3*D133 - b1*D311 -&
     & b2*D312 - b3*D313+0.d0
       X22 = 2*D202 - 2*(b1*D212 + b2*D222 + b3*D223)+0.d0
       X23 = D203 + D302 - b1*D213 - b2*D223 &
     & - b3*D233 - b1*D312 -&
     & b2*D322 - b3*D323+0.d0
       X33 = 2*D303 - 2*(b1*D313 + b2*D323 + b3*D333)+0.d0

       !--- the contracted G with the first index up------------ 
       trGu1 = G111*huu11 + G122*huu22 + 2*(G112*huu12 +&
     & G113*huu13 + G123*huu23) + G133*huu33+0.d0
       trGu2 = G211*huu11 + G222*huu22 + 2*(G212*huu12 +&
     & G213*huu13 + G223*huu23) + G233*huu33+0.d0
       trGu3 = G311*huu11 + G322*huu22 + 2*(G312*huu12 +&
     & G313*huu13 + G323*huu23) + G333*huu33+0.d0

      ! the derivative of the auxiliary fields --------------------
      dx_X12 = dx_D102 + dx_D201 - b1*dx_D112 &
     &       - b2*dx_D122 -&
     &       b3*dx_D123 - b1*dx_D211 - b2*dx_D212 &
     &       - b3*dx_D213 - d1b1*D112 -&
     &       d1b2*D122 - d1b3*D123 - d1b1*D211 &
     &       - d1b2*D212 - d1b3*D213+0.d0
      dx_X13 = dx_D103 + dx_D301 - b1*dx_D113 &
     &       - b2*dx_D123 -&
     &       b3*dx_D133 - b1*dx_D311 &
     &       - b2*dx_D312 - b3*dx_D313 &
     &       - d1b1*D113 -&
     &       d1b2*D123 - d1b3*D133 - d1b1*D311 &
     &       - d1b2*D312 - d1b3*D313+0.d0
      dx_X22 = 2*dx_D202 - 2*(b1*dx_D212 + b2*dx_D222 &
     &       + b3*dx_D223 +&
     &       d1b1*D212 + d1b2*D222 + d1b3*D223)+0.d0
      dx_X23 = dx_D203 + dx_D302 - b1*dx_D213 &
     &       - b2*dx_D223 -&
     &       b3*dx_D233 - b1*dx_D312 - b2*dx_D322 &
     &       - b3*dx_D323 - d1b1*D213 -&
     &       d1b2*D223 - d1b3*D233 - d1b1*D312 &
     &       - d1b2*D322 - d1b3*D323+0.d0
      dx_X33 = 2*dx_D303 - 2*(b1*dx_D313 + b2*dx_D323 + b3*dx_D333 +&
     & d1b1*D313 + d1b2*D323 + d1b3*D333)+0.d0
      dy_X11 = 2*dy_D101 - 2*(b1*dy_D111 + b2*dy_D112 + b3*dy_D113 +&
     & d2b1*D111 + d2b2*D112 + d2b3*D113)+0.d0
      dy_X12 = dy_D102 + dy_D201 - b1*dy_D112 &
     &       - b2*dy_D122 -&
     &         b3*dy_D123 - b1*dy_D211 - b2*dy_D212 &
     &       - b3*dy_D213 - d2b1*D112 -&
     &         d2b2*D122 - d2b3*D123 - d2b1*D211 &
     &       - d2b2*D212 - d2b3*D213+0.d0
      dy_X13 = dy_D103 + dy_D301 - b1*dy_D113 &
     &       - b2*dy_D123 -&
     &       b3*dy_D133 - b1*dy_D311 - b2*dy_D312 &
     &       - b3*dy_D313 - d2b1*D113 -&
     &         d2b2*D123 - d2b3*D133 - d2b1*D311 &
     &       - d2b2*D312 - d2b3*D313+0.d0
      dy_X23 = dy_D203 + dy_D302 - b1*dy_D213 &
     &       - b2*dy_D223 -&
     &         b3*dy_D233 - b1*dy_D312 - b2*dy_D322 &
     &       - b3*dy_D323 - d2b1*D213 -&
     &         d2b2*D223 - d2b3*D233 - d2b1*D312 &
     &       - d2b2*D322 - d2b3*D323+0.d0
      dy_X33 = 2*dy_D303 - 2*(b1*dy_D313 + b2*dy_D323 + b3*dy_D333 +&
     & d2b1*D313 + d2b2*D323 + d2b3*D333)+0.d0
      dz_X11 = 2*dz_D101 - 2*(b1*dz_D111 + b2*dz_D112 + b3*dz_D113 +&
     & d3b1*D111 + d3b2*D112 + d3b3*D113)+0.d0
      dz_X12 = dz_D102 + dz_D201 - b1*dz_D112 &
     &       - b2*dz_D122 -&
     &         b3*dz_D123 - b1*dz_D211 - b2*dz_D212 &
     &       - b3*dz_D213 - d3b1*D112 -&
     &         d3b2*D122 - d3b3*D123 - d3b1*D211 &
     &       - d3b2*D212 - d3b3*D213+0.d0
      dz_X13 = dz_D103 + dz_D301 - b1*dz_D113 &
     &       - b2*dz_D123 -&
     &         b3*dz_D133 - b1*dz_D311 - b2*dz_D312 &
     &       - b3*dz_D313 - d3b1*D113 -&
     &         d3b2*D123 - d3b3*D133 - d3b1*D311 &
     &       - d3b2*D312 - d3b3*D313+0.d0
      dz_X22 = 2*dz_D202 - 2*(b1*dz_D212 + b2*dz_D222 &
     &       + b3*dz_D223 +&
     &         d3b1*D212 + d3b2*D222 + d3b3*D223)+0.d0
      dz_X23 = dz_D203 + dz_D302 - b1*dz_D213 &
     &       - b2*dz_D223 -&
     &         b3*dz_D233 - b1*dz_D312 - b2*dz_D322 &
     &       - b3*dz_D323 - d3b1*D213 -&
     &         d3b2*D223 - d3b3*D233 - d3b1*D312 &
     &       - d3b2*D322 - d3b3*D323+0.d0

      ! the derivative of the real extrinsic curvature------------
      dx_eK12 = 0.5*(dx_K12 + dx_X12/alp) - (0.5*d1alp*X12)/alp**2+0.d0
      dx_eK13 = 0.5*(dx_K13 + dx_X13/alp) - (0.5*d1alp*X13)/alp**2+0.d0
      dx_eK22 = 0.5*(dx_K22 + dx_X22/alp) - (0.5*d1alp*X22)/alp**2+0.d0
      dx_eK23 = 0.5*(dx_K23 + dx_X23/alp) - (0.5*d1alp*X23)/alp**2+0.d0
      dx_eK33 = 0.5*(dx_K33 + dx_X33/alp) - (0.5*d1alp*X33)/alp**2+0.d0
      dy_eK11 = 0.5*(dy_K11 + dy_X11/alp) - (0.5*d2alp*X11)/alp**2+0.d0
      dy_eK12 = 0.5*(dy_K12 + dy_X12/alp) - (0.5*d2alp*X12)/alp**2+0.d0
      dy_eK13 = 0.5*(dy_K13 + dy_X13/alp) - (0.5*d2alp*X13)/alp**2+0.d0
      dy_eK23 = 0.5*(dy_K23 + dy_X23/alp) - (0.5*d2alp*X23)/alp**2+0.d0
      dy_eK33 = 0.5*(dy_K33 + dy_X33/alp) - (0.5*d2alp*X33)/alp**2+0.d0
      dz_eK11 = 0.5*(dz_K11 + dz_X11/alp) - (0.5*d3alp*X11)/alp**2+0.d0
      dz_eK12 = 0.5*(dz_K12 + dz_X12/alp) - (0.5*d3alp*X12)/alp**2+0.d0
      dz_eK13 = 0.5*(dz_K13 + dz_X13/alp) - (0.5*d3alp*X13)/alp**2+0.d0
      dz_eK22 = 0.5*(dz_K22 + dz_X22/alp) - (0.5*d3alp*X22)/alp**2+0.d0
      dz_eK23 = 0.5*(dz_K23 + dz_X23/alp) - (0.5*d3alp*X23)/alp**2+0.d0

      S1 = -(tu0*T01 + tu1*T11 + tu2*T12 + tu3*T13)
      S2 = -(tu0*T02 + tu1*T12 + tu2*T22 + tu3*T23)
      S3 = -(tu0*T03 + tu1*T13 + tu2*T23 + tu3*T33)

      ! the momentum constraint---------------------------------------
      momx_const = (-dx_eK12 + dy_eK11)*huu12 + (-dx_eK13 + dz_eK11)*huu13 +&
     & (-dx_eK22 + dy_eK12)*huu22 + (-2*dx_eK23 + dy_eK13 +&
     & dz_eK12)*huu23 + (-dx_eK33 + dz_eK13)*huu33 - 8*pi*S1 -&
     & eK11*trGu1 - eK12*trGu2 - eK13*trGu3 + eK11*Dduu111 +&
     & 2*eK12*Dduu112 + 2*eK13*Dduu113 + eK22*Dduu122 + 2*eK23*Dduu123 +&
     & eK33*Dduu133+0.d0
      momy_const = (dx_eK12 - dy_eK11)*huu11 + (dx_eK22 - dy_eK12)*huu12 +&
     & (dx_eK23 - 2*dy_eK13 + dz_eK12)*huu13 + (-dy_eK23 +&
     & dz_eK22)*huu23 + (-dy_eK33 + dz_eK23)*huu33 - 8*pi*S2 -&
     & eK12*trGu1 - eK22*trGu2 - eK23*trGu3 + eK11*Dduu211 +&
     & 2*eK12*Dduu212 + 2*eK13*Dduu213 + eK22*Dduu222 + 2*eK23*Dduu223 +&
     & eK33*Dduu233+0.d0
      momz_const = (dx_eK13 - dz_eK11)*huu11 + (dx_eK23 + dy_eK13 -&
     & 2*dz_eK12)*huu12 + (dx_eK33 - dz_eK13)*huu13 + (dy_eK23 -&
     & dz_eK22)*huu22 + (dy_eK33 - dz_eK23)*huu23 - 8*pi*S3 - eK13*trGu1&
     & - eK23*trGu2 - eK33*trGu3 + eK11*Dduu311 + 2*eK12*Dduu312 +&
     & 2*eK13*Dduu313 + eK22*Dduu322 + 2*eK23*Dduu323 + eK33*Dduu333+0.d0

!---- the Wyel tensor----------------------------------------------------
       dcx_eK12 = dx_eK12 - eK11*G112 - eK22*G211 - eK12*(G111 + G212) -&
     & eK23*G311 - eK13*G312+0.d0
       dcx_eK13 = dx_eK13 - eK11*G113 - eK23*G211 - eK12*G213 -&
     & eK33*G311 - eK13*(G111 + G313)+0.d0
       dcx_eK22 = dx_eK22 - 2*(eK12*G112 + eK22*G212 + eK23*G312)+0.d0
       dcx_eK23 = dx_eK23 - eK13*G112 - eK12*G113 - eK22*G213 -&
     & eK33*G312 - eK23*(G212 + G313)+0.d0
       dcx_eK33 = dx_eK33 - 2*(eK13*G113 + eK23*G213 + eK33*G313)+0.d0
       dcy_eK11 = dy_eK11 - 2*(eK11*G112 + eK12*G212 + eK13*G312)+0.d0
       dcy_eK12 = dy_eK12 - eK11*G122 - eK22*G212 - eK12*(G112 + G222) -&
     & eK23*G312 - eK13*G322+0.d0
       dcy_eK13 = dy_eK13 - eK11*G123 - eK23*G212 - eK12*G223 -&
     & eK33*G312 - eK13*(G112 + G323)+0.d0
       dcy_eK23 = dy_eK23 - eK13*G122 - eK12*G123 - eK22*G223 -&
     & eK33*G322 - eK23*(G222 + G323)+0.d0
       dcy_eK33 = dy_eK33 - 2*(eK13*G123 + eK23*G223 + eK33*G323)+0.d0
       dcz_eK11 = dz_eK11 - 2*(eK11*G113 + eK12*G213 + eK13*G313)+0.d0
       dcz_eK12 = dz_eK12 - eK11*G123 - eK22*G213 - eK12*(G113 + G223) -&
     & eK23*G313 - eK13*G323+0.d0
       dcz_eK13 = dz_eK13 - eK11*G133 - eK23*G213 - eK12*G233 -&
     & eK33*G313 - eK13*(G113 + G333)+0.d0
       dcz_eK22 = dz_eK22 - 2*(eK12*G123 + eK22*G223 + eK23*G323)+0.d0
       dcz_eK23 = dz_eK23 - eK13*G123 - eK12*G133 - eK22*G233 -&
     & eK33*G323 - eK23*(G223 + G333)+0.d0
      
       Cr11 = eK12*eKud21 + eK13*eKud31 + eK11*(eKud11 - etrK) + Ri11+0.d0
       Cr12 = eK22*eKud21 + eK23*eKud31 + eK12*(eKud11 - etrK) + Ri12+0.d0
       Cr13 = eK23*eKud21 + eK33*eKud31 + eK13*(eKud11 - etrK) + Ri13+0.d0
       Cr21 = eK11*eKud12 + eK13*eKud32 + eK12*(eKud22 - etrK) + Ri12+0.d0
       Cr22 = eK12*eKud12 + eK23*eKud32 + eK22*(eKud22 - etrK) + Ri22+0.d0
       Cr23 = eK13*eKud12 + eK33*eKud32 + eK23*(eKud22 - etrK) + Ri23+0.d0
       Cr31 = eK11*eKud13 + eK12*eKud23 + eK13*(eKud33 - etrK) + Ri13+0.d0
       Cr32 = eK12*eKud13 + eK22*eKud23 + eK23*(eKud33 - etrK) + Ri23+0.d0
       Cr33 = eK13*eKud13 + eK23*eKud23 + eK33*(eKud33 - etrK) + Ri33+0.d0

       Ci11 = dcy_eK13 - dcz_eK12+0.d0
       Ci12 = dcy_eK23 - dcz_eK22+0.d0
       Ci13 = dcy_eK33 - dcz_eK23+0.d0
       Ci21 = -dcx_eK13 + dcz_eK11+0.d0
       Ci22 = -dcx_eK23 + dcz_eK12+0.d0
       Ci23 = -dcx_eK33 + dcz_eK13+0.d0
       Ci31 = dcx_eK12 - dcy_eK11+0.d0
       Ci32 = dcx_eK22 - dcy_eK12+0.d0
       Ci33 = dcx_eK23 - dcy_eK13+0.d0

       ! the frame
       Tm_R(1)= (x*z)/(r*rcyl*sqrt(2.0d0))
       Tm_R(2)= (y*z)/(r*rcyl*sqrt(2.0d0))
       Tm_R(3)=-(x*x+y*y)/(r*rcyl*sqrt(2.0d0))
       Tm_R(4)= 0.0d0
       Tm_I(1)=-y/(rcyl*sqrt(2.0d0))
       Tm_I(2)= x/(rcyl*sqrt(2.0d0))
       Tm_I(3)= 0.0d0
       Tm_I(4)= 0.0d0
       Tl(1)= -x/(r*sqrt(2.0d0))
       Tl(2)= -y/(r*sqrt(2.0d0))
       Tl(3)= -z/(r*sqrt(2.0d0))
       Tl(4)= 1.0/sqrt(2.0d0)
       Tn(1)= x/(r*sqrt(2.0d0))
       Tn(2)= y/(r*sqrt(2.0d0))
       Tn(3)= z/(r*sqrt(2.0d0))
       Tn(4)= 1.0/sqrt(2.0d0)

       ! the psi4
       Psi4R = -(Cr11*Tm_I(1)**2) - Cr22*Tm_I(2)**2 - (Cr31*Tm_I(1) +&
     & (Cr23 + Cr32)*Tm_I(2))*Tm_I(3) - Cr33*Tm_I(3)**2 - Tm_I(1)*&
     & ((Cr12 + Cr21)*Tm_I(2) + Cr13*Tm_I(3)) + Cr11*Tm_R(1)**2 +&
     & Cr22*Tm_R(2)**2 + (Ci31*Tm_I(1) + (Ci23 + Ci32)*Tm_I(2) +&
     & 2*Ci33*Tm_I(3) + Cr31*Tm_R(1) + Cr32*Tm_R(2))* Tm_R(3) +&
     & Cr33*Tm_R(3)**2 + Tm_I(1)*((Ci12 + Ci21)*Tm_R(2) + Ci13*Tm_R(3))&
     & + Tm_R(1)*(2*Ci11*Tm_I(1) + (Ci12 + Ci21)*Tm_I(2) + (Ci13 +&
     & Ci31)*Tm_I(3) + Cr12*Tm_R(2) + Cr13*Tm_R(3)) + Tm_R(2)*&
     & (2*Ci22*Tm_I(2) + (Ci23 + Ci32)*Tm_I(3) + Cr21*Tm_R(1) +&
     & Cr23*Tm_R(3))+0.d0
       Psi4I = -(Ci11*Tm_I(1)**2) - Ci22*Tm_I(2)**2 - Ci33*Tm_I(3)**2 +&
     & Ci11*Tm_R(1)**2 - (Cr21*Tm_I(1) + Cr32*Tm_I(3))*Tm_R(2) +&
     & Ci22*Tm_R(2)**2 - (Cr31*Tm_I(1) + Cr32*Tm_I(2))* Tm_R(3) +&
     & (Ci31*Tm_R(1) + (Ci23 + Ci32)*Tm_R(2))* Tm_R(3) + Ci33*Tm_R(3)**2&
     & + Tm_R(1)*((Ci12 + Ci21)*Tm_R(2) + Ci13*Tm_R(3)) - Tm_I(1)*((Ci12&
     & + Ci21)*Tm_I(2) + Ci13*Tm_I(3) + 2*Cr11*Tm_R(1) + Cr12*Tm_R(2) +&
     & Cr13*Tm_R(3)) - Tm_I(2)*((Cr12 + Cr21)*Tm_R(1) + 2*Cr22*Tm_R(2) +&
     & Cr23*Tm_R(3)) - Tm_I(3)* (Ci31*Tm_I(1) + (Ci23 + Ci32)*Tm_I(2) +&
     & (Cr13 + Cr31)*Tm_R(1) + Cr23*Tm_R(2) + 2*Cr33*Tm_R(3))+0.d0

       Psi4R = Psi4R * r
       Psi4I = Psi4I * r

!-----the Noether density with the sqrt(h) from the volume factor----------

!      noether_density = 2.0d0*deth*(phir*pic - phic*pir)
      noether_density = sqrt(deth)*(phir*pic - phic*pir)&
     &                + sqrt(deth)*(phim*pin - phin*pim)

!-----the Zs quantities----------------------------------------------------
                  
       
       Gd000 = D000/2.
       Gd001 = D100/2.
       Gd002 = D200/2.
       Gd003 = D300/2.
       Gd011 = -D011/2. + D101
       Gd012 = (-D012 + D102 + D201)/2.
       Gd013 = (-D013 + D103 + D301)/2.
       Gd022 = -D022/2. + D202
       Gd023 = (-D023 + D203 + D302)/2.
       Gd033 = -D033/2. + D303
       Gd100 = D001 - D100/2.
       Gd101 = D011/2.
       Gd102 = (D012 - D102 + D201)/2.
       Gd103 = (D013 - D103 + D301)/2.
       Gd111 = D111/2.
       Gd112 = D211/2.
       Gd113 = D311/2.
       Gd122 = -D122/2. + D212
       Gd123 = (-D123 + D213 + D312)/2.
       Gd133 = -D133/2. + D313
       Gd200 = D002 - D200/2.
       Gd201 = (D012 + D102 - D201)/2.
       Gd202 = D022/2.
       Gd203 = (D023 - D203 + D302)/2.
       Gd211 = D112 - D211/2.
       Gd212 = D122/2.
       Gd213 = (D123 - D213 + D312)/2.
       Gd222 = D222/2.
       Gd223 = D322/2.
       Gd233 = -D233/2. + D323
       Gd300 = D003 - D300/2.
       Gd301 = (D013 + D103 - D301)/2.
       Gd302 = (D023 + D203 - D302)/2.
       Gd303 = D033/2.
       Gd311 = D113 - D311/2.
       Gd312 = (D123 + D213 - D312)/2.
       Gd313 = D133/2.
       Gd322 = D223 - D322/2.
       Gd323 = D233/2.
       Gd333 = D333/2.

       trG0 = guu00*Gd000 + guu11*Gd011 + guu22*Gd022 + guu33*Gd033 &
      &   +2.0*(guu01*Gd001 + guu02*Gd002 + guu03*Gd003 +&
      &         guu12*Gd012 + guu13*Gd013 + guu23*Gd023) 
      
       trG1 = guu00*Gd100 + guu11*Gd111 + guu22*Gd122 + guu33*Gd133 &
      &   +2.0*(guu01*Gd101 + guu02*Gd102 + guu03*Gd103 +&
      &         guu12*Gd112 + guu13*Gd113 + guu23*Gd123) 

       trG2 = guu00*Gd200 + guu11*Gd211 + guu22*Gd222 + guu33*Gd233 &
      &   +2.0*(guu01*Gd201 + guu02*Gd202 + guu03*Gd203 +&
      &         guu12*Gd212 + guu13*Gd213 + guu23*Gd223) 
 
       trG3 = guu00*Gd300 + guu11*Gd311 + guu22*Gd322 + guu33*Gd333 &
      &   +2.0*(guu01*Gd301 + guu02*Gd302 + guu03*Gd303 +&
      &         guu12*Gd312 + guu13*Gd313 + guu23*Gd323)       
                        
       Z0 = (-H0 - trG0)/2.
       Z1 = (-H1 - trG1)/2.
       Z2 = (-H2 - trG2)/2.
       Z3 = (-H3 - trG3)/2.
             
       physZ_const_pt = Z0**2 + Z1**2 + Z2**2 + Z3**2
       physZ_const_pt = sqrt(physZ_const_pt/4.0) 

      !------------------------------------------------------
      ! Calculate the components of the source term
      ! obeyed by the physZ constraints:
      !      \Box Z^\alpha = - R^\alpha_\nu Z^\nu
      ! (Eq. 2.8 of Choptuik&Sorkin http://arxiv.org/abs/0908.2500)
      ! Rewrite in terms of stress-energy tensor:
      !      \Box Z^\alpha = - Z^nu g^{alpha mu} [T_{mu nu}-(1/2)Tg_{mu nu}]
      !------------------------------------------------------
      trT  =   guu00*T00 + guu01*T01 + guu02*T02 + guu03*T03 &
     &       + guu01*T01 + guu11*T11 + guu12*T12 + guu13*T13 &
     &       + guu02*T02 + guu12*T12 + guu22*T22 + guu23*T23 &
     &       + guu03*T03 + guu13*T13 + guu23*T23 + guu33*T33 

!     z4rhs0_pt = -Z0*( guu00*( T00-0.5d0*trT*g00)  &
!    &                 +guu01*( T01-0.5d0*trT*g01)  &
!    &                 +guu02*( T02-0.5d0*trT*g02)  &
!    &                 +guu03*( T03-0.5d0*trT*g03) )
!     z4rhs1_pt = -Z1*( guu01*( T01-0.5d0*trT*g01)  &
!    &                 +guu11*( T11-0.5d0*trT*g11)  &
!    &                 +guu12*( T12-0.5d0*trT*g12)  &
!    &                 +guu13*( T13-0.5d0*trT*g13) )
!     z4rhs2_pt = -Z2*( guu02*( T02-0.5d0*trT*g02)  &
!    &                 +guu12*( T12-0.5d0*trT*g12)  &
!    &                 +guu22*( T22-0.5d0*trT*g22)  &
!    &                 +guu23*( T23-0.5d0*trT*g23) )
!     z4rhs3_pt = -Z3*( guu03*( T03-0.5d0*trT*g03)  &
!    &                 +guu13*( T13-0.5d0*trT*g13)  &
!    &                 +guu23*( T23-0.5d0*trT*g23)  &
!    &                 +guu33*( T33-0.5d0*trT*g33) )
      z4rhs0_pt = -     guu00*( T00-0.5d0*trT*g00)
      z4rhs1_pt = -     guu11*( T11-0.5d0*trT*g11)
      z4rhs2_pt = -     guu22*( T22-0.5d0*trT*g22)
      z4rhs3_pt = -     guu33*( T33-0.5d0*trT*g33)


       
      !-- compute the angular momentum from a killing----------

       xiu0 = 0.0d0
       xiu1 = -y
       xiu2 = x 
       xiu3 = 0.0d0
 
       xinorm = g00*xiu0*xiu0 + g11*xiu1*xiu1 + g22*xiu2*xiu2 + g33*xiu3*xiu3 &
     &      + 2.0d0*( g01*xiu0*xiu1 + g02*xiu0*xiu2 + g03*xiu0*xiu3 &
     &              + g12*xiu1*xiu2 + g13*xiu1*xiu3 + g23*xiu2*xiu3)

       if (abs(xinorm) .LT. 1E-8) xinorm=1E-8
          
       xiu0 = xiu0/xinorm; xiu1 = xiu1/xinorm; xiu2 = xiu2/xinorm; xiu3 = xiu3/xinorm
      
       kang_momz = - sqrt(deth)*&
     &	   (xiu0*tu0*T00 + xiu0*tu1*T01 + xiu0*tu2*T02 &
     &      + xiu0*tu3*T03 + &
     &      xiu1*tu0*T01 + xiu1*tu1*T11 + xiu1*tu2*T12 &
     &      + xiu1*tu3*T13 + &
     &      xiu2*tu0*T02 + xiu2*tu1*T12 + xiu2*tu2*T22 &
     &      + xiu2*tu3*T23 + &
     &      xiu3*tu0*T03 + xiu3*tu1*T13 + xiu3*tu2*T23 &
     &      + xiu3*tu3*T33 )
      
!--------------------------------------------------------------------------------
!-----computing the frame, the Riemann and the Psi4 -----------------------------       
!--------------------------------------------------------------------------------
       if (psi4_analysis .GE. 1) then

         !---------the 4-metric----------------------------------
         g(1,1) = g11; g(1,2) = g12; g(1,3) = g13; g(1,4) = g01
         g(2,1) = g12; g(2,2) = g22; g(2,3) = g23; g(2,4) = g02
         g(3,1) = g13; g(3,2) = g23; g(3,3) = g33; g(3,4) = g03
         g(4,1) = g01; g(4,2) = g02; g(4,3) = g03; g(4,4) = g00

         !---------the inverse of the 4-metric--------------------
         g_uu(1,1) = guu11; g_uu(1,2) = guu12; g_uu(1,3) = guu13; g_uu(1,4) = guu01
         g_uu(2,1) = guu12; g_uu(2,2) = guu22; g_uu(2,3) = guu23; g_uu(2,4) = guu02
         g_uu(3,1) = guu13; g_uu(3,2) = guu23; g_uu(3,3) = guu33; g_uu(3,4) = guu03
         g_uu(4,1) = guu01; g_uu(4,2) = guu02; g_uu(4,3) = guu03; g_uu(4,4) = guu00

         !---------the space derivatives of the 4-metric---------
         dg(1,1,1) = u_pt(H_D1G11); dg(1,1,2) = u_pt(H_D1G12); dg(1,1,3) = u_pt(H_D1G13); dg(1,1,4) = u_pt(H_D1G01)
         dg(1,2,1) = u_pt(H_D1G12); dg(1,2,2) = u_pt(H_D1G22); dg(1,2,3) = u_pt(H_D1G23); dg(1,2,4) = u_pt(H_D1G02)
         dg(1,3,1) = u_pt(H_D1G13); dg(1,3,2) = u_pt(H_D1G23); dg(1,3,3) = u_pt(H_D1G33); dg(1,3,4) = u_pt(H_D1G03)
         dg(1,4,1) = u_pt(H_D1G01); dg(1,4,2) = u_pt(H_D1G02); dg(1,4,3) = u_pt(H_D1G03); dg(1,4,4) = u_pt(H_D1G00)

         dg(2,1,1) = u_pt(H_D2G11); dg(2,1,2) = u_pt(H_D2G12); dg(2,1,3) = u_pt(H_D2G13); dg(2,1,4) = u_pt(H_D2G01)
         dg(2,2,1) = u_pt(H_D2G12); dg(2,2,2) = u_pt(H_D2G22); dg(2,2,3) = u_pt(H_D2G23); dg(2,2,4) = u_pt(H_D2G02)
         dg(2,3,1) = u_pt(H_D2G13); dg(2,3,2) = u_pt(H_D2G23); dg(2,3,3) = u_pt(H_D2G33); dg(2,3,4) = u_pt(H_D2G03)
         dg(2,4,1) = u_pt(H_D2G01); dg(2,4,2) = u_pt(H_D2G02); dg(2,4,3) = u_pt(H_D2G03); dg(2,4,4) = u_pt(H_D2G00)
        
         dg(3,1,1) = u_pt(H_D3G11); dg(3,1,2) = u_pt(H_D3G12); dg(3,1,3) = u_pt(H_D3G13); dg(3,1,4) = u_pt(H_D3G01)
         dg(3,2,1) = u_pt(H_D3G12); dg(3,2,2) = u_pt(H_D3G22); dg(3,2,3) = u_pt(H_D3G23); dg(3,2,4) = u_pt(H_D3G02)
         dg(3,3,1) = u_pt(H_D3G13); dg(3,3,2) = u_pt(H_D3G23); dg(3,3,3) = u_pt(H_D3G33); dg(3,3,4) = u_pt(H_D3G03)
         dg(3,4,1) = u_pt(H_D3G01); dg(3,4,2) = u_pt(H_D3G02); dg(3,4,3) = u_pt(H_D3G03); dg(3,4,4) = u_pt(H_D3G00)

         !---------the "time" derivatives of the 4-metric---------
         dg(4,1,1) = D011; dg(4,1,2) = D012; dg(4,1,3) = D013; dg(4,1,4) = D001
         dg(4,2,1) = D012; dg(4,2,2) = D022; dg(4,2,3) = D023; dg(4,2,4) = D002
         dg(4,3,1) = D013; dg(4,3,2) = D023; dg(4,3,3) = D033; dg(4,3,4) = D003
         dg(4,4,1) = D001; dg(4,4,2) = D002; dg(4,4,3) = D003; dg(4,4,4) = D000

	 !---------the space derivates of the 4-metric------------
	 D_dg = 0.0d0
         !---x direction -------------	
	 D_dg(1,1,1,1) = dxu_pt(H_D1G11); D_dg(1,1,1,2) = dxu_pt(H_D1G12)
	 D_dg(1,1,1,3) = dxu_pt(H_D1G13); D_dg(1,1,1,4) = dxu_pt(H_D1G01)
         D_dg(1,1,2,1) = dxu_pt(H_D1G12); D_dg(1,1,2,2) = dxu_pt(H_D1G22)
	 D_dg(1,1,2,3) = dxu_pt(H_D1G23); D_dg(1,1,2,4) = dxu_pt(H_D1G02)
         D_dg(1,1,3,1) = dxu_pt(H_D1G13); D_dg(1,1,3,2) = dxu_pt(H_D1G23)
	 D_dg(1,1,3,3) = dxu_pt(H_D1G33); D_dg(1,1,3,4) = dxu_pt(H_D1G03)
         D_dg(1,1,4,1) = dxu_pt(H_D1G01); D_dg(1,1,4,2) = dxu_pt(H_D1G02)
	 D_dg(1,1,4,3) = dxu_pt(H_D1G03); D_dg(1,1,4,4) = dxu_pt(H_D1G00)

         D_dg(1,2,1,1) = dxu_pt(H_D2G11); D_dg(1,2,1,2) = dxu_pt(H_D2G12)
	 D_dg(1,2,1,3) = dxu_pt(H_D2G13); D_dg(1,2,1,4) = dxu_pt(H_D2G01)
         D_dg(1,2,2,1) = dxu_pt(H_D2G12); D_dg(1,2,2,2) = dxu_pt(H_D2G22)
	 D_dg(1,2,2,3) = dxu_pt(H_D2G23); D_dg(1,2,2,4) = dxu_pt(H_D2G02)
         D_dg(1,2,3,1) = dxu_pt(H_D2G13); D_dg(1,2,3,2) = dxu_pt(H_D2G23)
	 D_dg(1,2,3,3) = dxu_pt(H_D2G33); D_dg(1,2,3,4) = dxu_pt(H_D2G03)
         D_dg(1,2,4,1) = dxu_pt(H_D2G01); D_dg(1,2,4,2) = dxu_pt(H_D2G02)
	 D_dg(1,2,4,3) = dxu_pt(H_D2G03); D_dg(1,2,4,4) = dxu_pt(H_D2G00)

         D_dg(1,3,1,1) = dxu_pt(H_D3G11); D_dg(1,3,1,2) = dxu_pt(H_D3G12)
         D_dg(1,3,1,3) = dxu_pt(H_D3G13); D_dg(1,3,1,4) = dxu_pt(H_D3G01)
         D_dg(1,3,2,1) = dxu_pt(H_D3G12); D_dg(1,3,2,2) = dxu_pt(H_D3G22)
         D_dg(1,3,2,3) = dxu_pt(H_D3G23); D_dg(1,3,2,4) = dxu_pt(H_D3G02)
         D_dg(1,3,3,1) = dxu_pt(H_D3G13); D_dg(1,3,3,2) = dxu_pt(H_D3G23)
         D_dg(1,3,3,3) = dxu_pt(H_D3G33); D_dg(1,3,3,4) = dxu_pt(H_D3G03)
         D_dg(1,3,4,1) = dxu_pt(H_D3G01); D_dg(1,3,4,2) = dxu_pt(H_D3G02)
         D_dg(1,3,4,3) = dxu_pt(H_D3G03); D_dg(1,3,4,4) = dxu_pt(H_D3G00)

         !---y direction -------------
         D_dg(2,1,1,1) = dyu_pt(H_D1G11); D_dg(2,1,1,2) = dyu_pt(H_D1G12)
         D_dg(2,1,1,3) = dyu_pt(H_D1G13); D_dg(2,1,1,4) = dyu_pt(H_D1G01)
         D_dg(2,1,2,1) = dyu_pt(H_D1G12); D_dg(2,1,2,2) = dyu_pt(H_D1G22)
         D_dg(2,1,2,3) = dyu_pt(H_D1G23); D_dg(2,1,2,4) = dyu_pt(H_D1G02)
         D_dg(2,1,3,1) = dyu_pt(H_D1G13); D_dg(2,1,3,2) = dyu_pt(H_D1G23)
         D_dg(2,1,3,3) = dyu_pt(H_D1G33); D_dg(2,1,3,4) = dyu_pt(H_D1G03)
         D_dg(2,1,4,1) = dyu_pt(H_D1G01); D_dg(2,1,4,2) = dyu_pt(H_D1G02)
         D_dg(2,1,4,3) = dyu_pt(H_D1G03); D_dg(2,1,4,4) = dyu_pt(H_D1G00)

         D_dg(2,2,1,1) = dyu_pt(H_D2G11); D_dg(2,2,1,2) = dyu_pt(H_D2G12)
         D_dg(2,2,1,3) = dyu_pt(H_D2G13); D_dg(2,2,1,4) = dyu_pt(H_D2G01)
         D_dg(2,2,2,1) = dyu_pt(H_D2G12); D_dg(2,2,2,2) = dyu_pt(H_D2G22)
         D_dg(2,2,2,3) = dyu_pt(H_D2G23); D_dg(2,2,2,4) = dyu_pt(H_D2G02)
         D_dg(2,2,3,1) = dyu_pt(H_D2G13); D_dg(2,2,3,2) = dyu_pt(H_D2G23)
         D_dg(2,2,3,3) = dyu_pt(H_D2G33); D_dg(2,2,3,4) = dyu_pt(H_D2G03)
         D_dg(2,2,4,1) = dyu_pt(H_D2G01); D_dg(2,2,4,2) = dyu_pt(H_D2G02)
         D_dg(2,2,4,3) = dyu_pt(H_D2G03); D_dg(2,2,4,4) = dyu_pt(H_D2G00)

         D_dg(2,3,1,1) = dyu_pt(H_D3G11); D_dg(2,3,1,2) = dyu_pt(H_D3G12)
         D_dg(2,3,1,3) = dyu_pt(H_D3G13); D_dg(2,3,1,4) = dyu_pt(H_D3G01)
         D_dg(2,3,2,1) = dyu_pt(H_D3G12); D_dg(2,3,2,2) = dyu_pt(H_D3G22)
         D_dg(2,3,2,3) = dyu_pt(H_D3G23); D_dg(2,3,2,4) = dyu_pt(H_D3G02)
         D_dg(2,3,3,1) = dyu_pt(H_D3G13); D_dg(2,3,3,2) = dyu_pt(H_D3G23)
         D_dg(2,3,3,3) = dyu_pt(H_D3G33); D_dg(2,3,3,4) = dyu_pt(H_D3G03)
         D_dg(2,3,4,1) = dyu_pt(H_D3G01); D_dg(2,3,4,2) = dyu_pt(H_D3G02)
         D_dg(2,3,4,3) = dyu_pt(H_D3G03); D_dg(2,3,4,4) = dyu_pt(H_D3G00)

         !---z direction -------------
         D_dg(3,1,1,1) = dzu_pt(H_D1G11); D_dg(3,1,1,2) = dzu_pt(H_D1G12)
         D_dg(3,1,1,3) = dzu_pt(H_D1G13); D_dg(3,1,1,4) = dzu_pt(H_D1G01)
         D_dg(3,1,2,1) = dzu_pt(H_D1G12); D_dg(3,1,2,2) = dzu_pt(H_D1G22)
         D_dg(3,1,2,3) = dzu_pt(H_D1G23); D_dg(3,1,2,4) = dzu_pt(H_D1G02)
         D_dg(3,1,3,1) = dzu_pt(H_D1G13); D_dg(3,1,3,2) = dzu_pt(H_D1G23)
         D_dg(3,1,3,3) = dzu_pt(H_D1G33); D_dg(3,1,3,4) = dzu_pt(H_D1G03)
         D_dg(3,1,4,1) = dzu_pt(H_D1G01); D_dg(3,1,4,2) = dzu_pt(H_D1G02)
         D_dg(3,1,4,3) = dzu_pt(H_D1G03); D_dg(3,1,4,4) = dzu_pt(H_D1G00)

         D_dg(3,2,1,1) = dzu_pt(H_D2G11); D_dg(3,2,1,2) = dzu_pt(H_D2G12)
         D_dg(3,2,1,3) = dzu_pt(H_D2G13); D_dg(3,2,1,4) = dzu_pt(H_D2G01)
         D_dg(3,2,2,1) = dzu_pt(H_D2G12); D_dg(3,2,2,2) = dzu_pt(H_D2G22)
         D_dg(3,2,2,3) = dzu_pt(H_D2G23); D_dg(3,2,2,4) = dzu_pt(H_D2G02)
         D_dg(3,2,3,1) = dzu_pt(H_D2G13); D_dg(3,2,3,2) = dzu_pt(H_D2G23)
         D_dg(3,2,3,3) = dzu_pt(H_D2G33); D_dg(3,2,3,4) = dzu_pt(H_D2G03)
         D_dg(3,2,4,1) = dzu_pt(H_D2G01); D_dg(3,2,4,2) = dzu_pt(H_D2G02)
         D_dg(3,2,4,3) = dzu_pt(H_D2G03); D_dg(3,2,4,4) = dzu_pt(H_D2G00)

         D_dg(3,3,1,1) = dzu_pt(H_D3G11); D_dg(3,3,1,2) = dzu_pt(H_D3G12)
         D_dg(3,3,1,3) = dzu_pt(H_D3G13); D_dg(3,3,1,4) = dzu_pt(H_D3G01)
         D_dg(3,3,2,1) = dzu_pt(H_D3G12); D_dg(3,3,2,2) = dzu_pt(H_D3G22)
         D_dg(3,3,2,3) = dzu_pt(H_D3G23); D_dg(3,3,2,4) = dzu_pt(H_D3G02)
         D_dg(3,3,3,1) = dzu_pt(H_D3G13); D_dg(3,3,3,2) = dzu_pt(H_D3G23)
         D_dg(3,3,3,3) = dzu_pt(H_D3G33); D_dg(3,3,3,4) = dzu_pt(H_D3G03) 
         D_dg(3,3,4,1) = dzu_pt(H_D3G01); D_dg(3,3,4,2) = dzu_pt(H_D3G02)
         D_dg(3,3,4,3) = dzu_pt(H_D3G03); D_dg(3,3,4,4) = dzu_pt(H_D3G00) 

         !---t direction -------------
         call compute_rhs_analysis(urk_pt, u_pt, v_pt, dxu_pt, dyu_pt, dzu_pt, &
                                   dxv_pt, dyv_pt, dzv_pt, w_pt, par , D_dg)
              
         !--compute the frame--------- 
         call frame(Tl, Tn, Tm_R, Tm_I, x, y, z, alp, beta, g, g_uu, deth)

         !--compute the 4D Rieman-----
         call rieman(R4, Chd, g_uu, g, dg, D_dg, 4)
 
         !-- compute the psi 4--------
         call calc_psi4(Psi4R, Psi4I, R4, Tn, Tm_R, Tm_I)
    
         Psi4R = Psi4R * r
         Psi4I = Psi4I * r

         !---the Bondi gur -----------------------------------
         !----------------------------------------------------
         !--------- r_a = (l_a - n_a)/sqrt(2)-----------------
 
         gur = 0.0d0
         do la=1,4
           do lb=1,4
             gur = gur + g(la,lb)*Tn(la)*((Tl(lb)-Tn(lb))/sqrt(2.0d0))
           end do
         end do

         !-- compute the ADM mass ----------------------------------
         ADMmass = 0.0         
         do li=1,3
           do lm=1,3
             ADMmass = ADMmass + &
     &          (dg(li,li,lm)-dg(lm,li,li))*(Tl(lm)-Tn(lm))*r**2
           end do
         end do
         ADMmass = ADMmass/(16*pi*sqrt(2.0d0)) 
 
         !-- compute the linear and angular momentum-----------------
         !--the M quantities------------
         KtrK(1,1) = eK11 - g11*etrK
         KtrK(1,2) = eK12 - g12*etrK
         KtrK(1,3) = eK13 - g13*etrK
         KtrK(2,2) = eK22 - g22*etrK
         KtrK(2,3) = eK23 - g23*etrK
         KtrK(3,3) = eK33 - g33*etrK

         KtrK(2,1) = KtrK(1,2)
         KtrK(3,1) = KtrK(1,3)
         KtrK(3,2) = KtrK(2,3)

         M = 0.0d0
         do li=1,3; do lm=1,3
            M(li) = M(li) + KtrK(li,lm)*(Tl(lm)-Tn(lm))*r*r
         end do; end do  

         lin_momx = M(1)/(8*pi*sqrt(2.0d0)) 
         lin_momy = M(2)/(8*pi*sqrt(2.0d0)) 
         lin_momz = M(3)/(8*pi*sqrt(2.0d0))  	

         ang_momx = (y*lin_momz - z*lin_momy)
         ang_momy = (z*lin_momx - x*lin_momz)
         ang_momz = (x*lin_momy - y*lin_momx)
 
       end if

       
!--------------------------------------------------------------------------       
!----- compute the curvature ----------------------------------------------
!--------------------------------------------------------------------------       

       if ((psi4_analysis .EQ. 2) .OR. (psi4_analysis .EQ. 4)) then         
              
         call cal_curvature(curvature, R4, g_uu)
 
       end if       
              
!--------------------------------------------------------------------------       
!----- compute the 2-curvature on the spheres normal to n------------------
!--------------------------------------------------------------------------       

       if (psi4_analysis .GE. 3) then         

         !-- compute the 3 curvature -------------------------------  
         K(1,1) = eK11; K(1,2) = eK12; K(1,3) = eK13 
         K(2,1) = eK12; K(2,2) = eK22; K(2,3) = eK23 
         K(3,1) = eK13; K(3,2) = eK23; K(3,3) = eK33 

         h_uu(1,1) = huu11; h_uu(1,2) = huu12; h_uu(1,3) = huu13
         h_uu(2,1) = huu12; h_uu(2,2) = huu22; h_uu(2,3) = huu23
         h_uu(3,1) = huu13; h_uu(3,2) = huu23; h_uu(3,3) = huu33
  
         R3 = 0.0d0	 
         do li=1,3; do lj=1,3; do lk=1,3; do lm=1,3 
             R3(li,lj,lk,lm) = R4(li,lj,lk,lm) - K(li,lk)*K(lj,lm) + K(li,lm)*K(lk,lj)
         end do; end do; end do; end do
  
         Ricci3 = 0.0d0
         do li=1,3; do lj=1,3; do lk=1,3; do lm=1,3 
            Ricci3(li,lj) = Ricci3(li,lj) + R3(lk,li,lm,lj)*h_uu(lk,lm)
         end do; end do; end do; end do
          
         R3scalar = 0.0d0
         do li=1,3; do lj=1,3
            R3scalar = R3scalar + Ricci3(li,lj)*h_uu(li,lj)
         end do; end do
 
        !--other auxiliary variables -------------------	 
         nr1 = x; nr2 = y; nr3 = z             
         nrnorm = sqrt(huu11*nr1*nr1 + huu22*nr2*nr2 + huu33*nr3*nr3 +& 
     &            2.0d0*(huu12*nr1*nr2 + huu13*nr1*nr3 + huu23*nr2*nr3))
         if (nrnorm .LT. 1E-8) nrnorm=1E-8 
         nr1 = nr1/nrnorm
         nr2 = nr2/nrnorm
         nr3 = nr3/nrnorm      
                                                 
         nru1 = huu11*nr1 + huu12*nr2 + huu13*nr3
         nru2 = huu12*nr1 + huu22*nr2 + huu23*nr3
         nru3 = huu13*nr1 + huu23*nr2 + huu33*nr3
                            
         Rnn = Ricci3(1,1)*nru1*nru1 + Ricci3(2,2)*nru2*nru2 + Ricci3(3,3)*nru3*nru3 &
     &    + 2.0d0*(Ricci3(1,2)*nru1*nru2 + Ricci3(1,3)*nru1*nru3 + Ricci3(2,3)*nru2*nru3)
 
         Gndd11 = nr1*G111 + nr2*G211 + nr3*G311
         Gndd12 = nr1*G112 + nr2*G212 + nr3*G312
         Gndd13 = nr1*G113 + nr2*G213 + nr3*G313
         Gndd22 = nr1*G122 + nr2*G222 + nr3*G322
         Gndd23 = nr1*G123 + nr2*G223 + nr3*G323
         Gndd33 = nr1*G133 + nr2*G233 + nr3*G333
              
         Dnn1 = D111*nru1*nru1 + D122*nru2*nru2 + D133*nru3*nru3 + &
     &    2.0d0*(D112*nru1*nru2 + D113*nru1*nru3 + D123*nru2*nru3) 
         Dnn2 = D211*nru1*nru1 + D222*nru2*nru2 + D233*nru3*nru3 + &
     &    2.0d0*(D212*nru1*nru2 + D213*nru1*nru3 + D223*nru2*nru3) 
         Dnn3 = D311*nru1*nru1 + D322*nru2*nru2 + D333*nru3*nru3 + &
     &    2.0d0*(D312*nru1*nru2 + D313*nru1*nru3 + D323*nru2*nru3) 
       
         D1n1 = (1.0d0 - nr1*nru1)/nrnorm + 0.50d0*Dnn1*nr1 - Gndd11
         D2n1 = (0.0d0 - nr1*nru2)/nrnorm + 0.50d0*Dnn2*nr1 - Gndd12
         D3n1 = (0.0d0 - nr1*nru3)/nrnorm + 0.50d0*Dnn3*nr1 - Gndd13                     
         D1n2 = (0.0d0 - nr2*nru1)/nrnorm + 0.50d0*Dnn1*nr2 - Gndd12
         D2n2 = (1.0d0 - nr2*nru2)/nrnorm + 0.50d0*Dnn2*nr2 - Gndd22
         D3n2 = (0.0d0 - nr2*nru3)/nrnorm + 0.50d0*Dnn3*nr2 - Gndd23              
         D1n3 = (0.0d0 - nr3*nru1)/nrnorm + 0.50d0*Dnn1*nr3 - Gndd13
         D2n3 = (0.0d0 - nr3*nru2)/nrnorm + 0.50d0*Dnn2*nr3 - Gndd23
         D3n3 = (1.0d0 - nr3*nru3)/nrnorm + 0.50d0*Dnn3*nr3 - Gndd33                     

         !-define the 2 projector over the hs----
         Sud11 = 1.0d0 - nru1*nr1
         Sud12 = - nru1*nr2
         Sud13 = - nru1*nr3	
         Sud21 = - nru2*nr1
         Sud22 = 1.0d0 - nru2*nr2
         Sud23 = - nru2*nr3        
         Sud31 = - nru3*nr1
         Sud32 = - nru3*nr2
         Sud33 = 1.0d0 - nru3*nr3

         Suu11 = huu11 - nru1*nru1
         Suu12 = huu12 - nru1*nru2
         Suu13 = huu13 - nru1*nru3	
         Suu22 = huu22 - nru2*nru2
         Suu23 = huu23 - nru2*nru3        
         Suu33 = huu33 - nru3*nru3
         
	 sK11 = D1n1*Sud11**2 + D2n2*Sud21**2 + (D3n1*Sud11 + (D2n3 +&
     &   D3n2)*Sud21)*Sud31 + D3n3*Sud31**2 + Sud11* ((D1n2 + D2n1)*Sud21&
     &   + D1n3*Sud31)
         sK12 = Sud12*(D1n1*Sud11 + D1n2*Sud21 + D1n3*Sud31) +&
     &   Sud22*(D2n1*Sud11 + D2n2*Sud21 + D2n3*Sud31) + (D3n1*Sud11 +&
     &   D3n2*Sud21 + D3n3*Sud31)*Sud32
         sK13 = Sud13*(D1n1*Sud11 + D1n2*Sud21 + D1n3*Sud31) +&
     &   Sud23*(D2n1*Sud11 + D2n2*Sud21 + D2n3*Sud31) + (D3n1*Sud11 +&
     &   D3n2*Sud21 + D3n3*Sud31)*Sud33
         sK22 = D1n1*Sud12**2 + D2n2*Sud22**2 + (D3n1*Sud12 + (D2n3 +&
     &   D3n2)*Sud22)*Sud32 + D3n3*Sud32**2 + Sud12* ((D1n2 + D2n1)*Sud22&
     &   + D1n3*Sud32)
         sK23 = Sud13*(D1n1*Sud12 + D1n2*Sud22 + D1n3*Sud32) +&
     &   Sud23*(D2n1*Sud12 + D2n2*Sud22 + D2n3*Sud32) + (D3n1*Sud12 +&
     &   D3n2*Sud22 + D3n3*Sud32)*Sud33
         sK33 = D1n1*Sud13**2 + D2n2*Sud23**2 + (D3n1*Sud13 + (D2n3 +&
     &   D3n2)*Sud23)*Sud33 + D3n3*Sud33**2 + Sud13* ((D1n2 + D2n1)*Sud23&
     &   + D1n3*Sud33)
           
         sKud11 = Suu11*sK11 + Suu12*sK12 + Suu13*sK13
         sKud12 = Suu11*sK12 + Suu12*sK22 + Suu13*sK23
         sKud13 = Suu11*sK13 + Suu12*sK23 + Suu13*sK33
         sKud21 = Suu12*sK11 + Suu22*sK12 + Suu23*sK13
         sKud22 = Suu12*sK12 + Suu22*sK22 + Suu23*sK23
         sKud23 = Suu12*sK13 + Suu22*sK23 + Suu23*sK33
         sKud31 = Suu13*sK11 + Suu23*sK12 + Suu33*sK13
         sKud32 = Suu13*sK12 + Suu23*sK22 + Suu33*sK23
         sKud33 = Suu13*sK13 + Suu23*sK23 + Suu33*sK33
         strK = sKud11 + sKud22 + sKud33       

         sKK = sKud11*sKud11 + sKud22*sKud22 + sKud33*sKud33 &
     &     + 2.0d0*(sKud12*sKud21 + sKud13*sKud31 + sKud23*sKud32)	 
       
         !---finally, the curvature in the 2-surface ----------
         R2scalar = R3scalar - 2*Rnn + strK*strK - sKK
         R2scalar = R2scalar*r*r
  
         !---the intrinsic metric induced on the 2-sphere ----
         !--- in spherical coordinates -----------------------
          
         Sdd11 = (h11 - nr1*nr1)
         Sdd12 = (h12 - nr1*nr2)
         Sdd13 = (h13 - nr1*nr3)
         Sdd22 = (h22 - nr2*nr2)
         Sdd23 = (h23 - nr2*nr3)        
         Sdd33 = (h33 - nr3*nr3)

         rho = sqrt(x*x + y*y)
         if ( rho .LT. 1E-6 ) rho = 1E-6  
         costheta = z/r
         sintheta = rho/r
         cosphi   = x/rho
         sinphi   = y/rho
           
         JA(1,1) = cosphi*sintheta
         JA(2,1) = sinphi*sintheta
         JA(3,1) = cosphi

         JA(1,2) = r*cosphi*costheta
         JA(2,2) = r*sinphi*costheta
         JA(3,2) = -r*sintheta
  
         JA(1,3) = -r*sinphi*sintheta
         JA(2,3) = r*cosphi*sintheta
         JA(3,3) = 0.0
     
!         Sdds11 = JA(1,1)*JA(1,1)*Sdd11 + JA(2,1)*JA(2,1)*Sdd22 + JA(3,1)*JA(3,1)*Sdd33 &
!     &     + 2.0*(JA(1,1)*JA(2,1)*Sdd12 + JA(1,1)*JA(3,1)*Sdd13 + JA(2,1)*JA(3,1)*Sdd23) 
         Sdds22 = JA(1,2)*JA(1,2)*Sdd11 + JA(2,2)*JA(2,2)*Sdd22 + JA(3,2)*JA(3,2)*Sdd33 &
     &     + 2.0*(JA(1,2)*JA(2,2)*Sdd12 + JA(1,2)*JA(3,2)*Sdd13 + JA(2,2)*JA(3,2)*Sdd23) 
         Sdds33 = JA(1,3)*JA(1,3)*Sdd11 + JA(2,3)*JA(2,3)*Sdd22 + JA(3,3)*JA(3,3)*Sdd33 &
     &     + 2.0*(JA(1,3)*JA(2,3)*Sdd12 + JA(1,3)*JA(3,3)*Sdd13 + JA(2,3)*JA(3,3)*Sdd23) 
         Sdds23 = JA(1,2)*JA(1,3)*Sdd11 + JA(2,2)*JA(2,3)*Sdd22 + JA(3,2)*JA(3,3)*Sdd33 &
     &     + 2.0*(JA(1,2)*JA(2,3)*Sdd12 + JA(1,2)*JA(3,3)*Sdd13 + JA(2,2)*JA(3,3)*Sdd23) 

         Sdds22 = Sdds22/(r*r)
         Sdds23 = Sdds23/(r*r)               
         Sdds33 = Sdds33/(r*r)
       
       end if
    

    v_pt(H_ORDER1_CONST)    = order1_const_pt
    v_pt(H_ORDER2_CONST)    = order2_const_pt
    v_pt(H_PHYSZ_CONST)     = physZ_const_pt
    v_pt(H_ENERGY_DENS)     = energy_density
    v_pt(H_ENERGY_CONST)    = energy_const
    v_pt(H_MOMX_CONST)      = momx_const
    v_pt(H_MOMY_CONST)      = momy_const
    v_pt(H_MOMZ_CONST)      = momz_const
    v_pt(H_NOETHER_DENS)    = noether_density    
    v_pt(H_ALPHA)           = alp    
    v_pt(H_SQDETG)          = sqrt(deth)
    
    v_pt(H_PSI4R)           = Psi4R 
    v_pt(H_PSI4I)           = Psi4I
    v_pt(H_massADM)         = ADMmass
    v_pt(H_CURVATURE)       = curvature
    v_pt(H_R2SCALAR)        = R2scalar
    v_pt(H_GTHETHE)         = Sdds22
    v_pt(H_GTHEPHI)         = Sdds23
    v_pt(H_GPHIPHI)         = Sdds33        
    v_pt(H_GUR)             = gur
!    v_pt(H_RAD_EXP)        = rad_exp_pt 
!    v_pt(H_RAD_SPEED)      = rad_speed_pt
    v_pt(H_LIN_MOMX)        = lin_momx
    v_pt(H_LIN_MOMY)        = lin_momy    
    v_pt(H_ANG_MOMZ)        = ang_momz    
    v_pt(H_KANG_MOMZ)       = kang_momz

    w_pt(H_Z4RHS0)          = z4rhs0_pt
    w_pt(H_Z4RHS1)          = z4rhs1_pt
    w_pt(H_Z4RHS2)          = z4rhs2_pt
    w_pt(H_Z4RHS3)          = z4rhs3_pt
                   
    ! -- extra dissipation around the BHs--------------

    bhdist = sqrt( (par(P_BH1_X0)-par(P_BH2_X0))**2 &
    &             +(par(P_BH1_Y0)-par(P_BH2_Y0))**2 &
    &             +(par(P_BH1_Z0)-par(P_BH2_Z0))**2 )

     absm1 = abs(par(P_BH1_mass))
     absm2 = abs(par(P_BH2_mass))
 
    bhsize = absm1 + absm2 

    !
    Exx = x
    Eyy = y
    Ezz = z

    W_pt(H_WDISS) = 1.0d0

    rcoor = sqrt(x**2+y**2+z**2)
    if(rcoor.gt.50.*bhsize.and.(rcoor).lt.60*bhsize) then
       w_pt(H_WDISS) = 1. + par(P_EXTRADISSOUT)*&
                     &(rcoor-50.*bhsize)/(10.*bhsize)
    else if (rcoor.ge.60.*bhsize) then
       w_pt(H_WDISS) = 1. + par(P_EXTRADISSOUT)
    end if


    nummask = masks_return_num()


    if(par(P_NBHOLES).eq.1.or.int(nummask).eq.1) then
!              write(*,*) 'pero esto', par(P_BH1_X0),par(P_BH2_X0)
!              write(*,*) 'pero esto', par(P_NBHOLES), nummask

          if(par(P_EXTRADISS).gt.0) then
                 rt = sqrt((Exx-par(P_BH1_X0))**2+(Eyy-par(P_BH1_Y0))**2+&
               &           (Ezz-par(P_BH1_Z0))**2)

                 if( rt.le.3.*absm1) then
                  W_pt(H_WDISS) = 1. + par(P_EXTRADISS)
                 else if (rt.gt.3.*absm1.and.rt.le.5.5*absm1) then
                  W_pt(H_WDISS) = 1. + par(P_EXTRADISS)* &
               &                    (5.5*absm1-rt)/(2.5*absm1)
                 end if
           end if
    else if (par(P_NBHOLES).eq.2.or.int(nummask).eq.2) then

           IF(bhdist .gt. 5.5*bhsize .AND. par(P_TIME).gt.1e-3) THEN
!             write(*,*) 'trato esto', par(P_BH1_X0),par(P_BH2_X0)
              if(par(P_EXTRADISS).gt.0) then
                 rt = sqrt((Exx-par(P_BH1_X0))**2+(Eyy-par(P_BH1_Y0))**2+&
               &           (Ezz-par(P_BH1_Z0))**2)

                 if( rt.le.3.*absm1) then
                  W_pt(H_WDISS) = 1. + par(P_EXTRADISS)
                 else if (rt.gt.3.*absm1.and.rt.le.5*absm1) then
                  W_pt(H_WDISS) = 1. + par(P_EXTRADISS)* &
               &                    (5*absm1-rt)/(2.*absm1)
                 end if

                 rt = sqrt((Exx-par(P_BH2_X0))**2+(Eyy-par(P_BH2_Y0))**2+&
               &           (Ezz-par(P_BH2_Z0))**2)

                 if( rt.le.3.*absm2) then
                  W_pt(H_WDISS) = 1. + par(P_EXTRADISS)
                 else if (rt.gt.3.*absm2.and.rt.le.5*absm2) then
                  W_pt(H_WDISS) = 1. + par(P_EXTRADISS)* &
               &                    (5*absm2-rt)/(2.*absm2)
                 end if
              end if

            ELSE IF(bhdist .le.5.5*bhsize .AND. par(P_TIME).gt.1e-3) THEN

!           write(*,*) 'hago esto', absm1*par(P_BH1_X0)+absm2*par(P_BH2_X0)
            Exx = x - (absm1*par(P_BH1_X0)+absm2*par(P_BH2_X0))/bhsize
            Eyy = y - (absm1*par(P_BH1_y0)+absm2*par(P_BH2_y0))/bhsize
            Ezz = z - (absm1*par(P_BH1_z0)+absm2*par(P_BH2_z0))/bhsize

             if(par(P_EXTRADISS).gt.0) then
                 rt = sqrt(Exx**2+Eyy**2+Ezz**2)

                 if( rt.le.3.*bhsize) then
                  W_pt(H_WDISS) = 1. + par(P_EXTRADISS)
                 else if (rt.gt.3.*bhsize.and.rt.le.5.5*bhsize) then
                  W_pt(H_WDISS) = 1. + par(P_EXTRADISS)* &
               &                    (5.5*bhsize-rt)/(2.5*bhsize)
                 end if

              end if

            END IF

      end if
 
  
    return
  end subroutine


      
 end module
