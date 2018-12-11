!-------------------------------------------------------------------------
!
!  $Id: rhs.f90,v 1.50 2010-05-04 18:30:35 tgarrett Exp $
!
!-------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_DefineThorn.h"

#include "basic.h"
#include "Dus.h"
#include "GH_rhs.h"

#if defined GH_MHD || defined GH_FLOWER
  module gh_rhs_mine
#else
  module rhs_mine
#endif

  use params
#if defined GH_MHD || defined GH_FLOWER
  use gh_auxvars
#else
  use auxvars
#endif

  implicit none

contains

  !---------------------------------------------------------------------
  !  
  ! SUBROUTINE CALCRHS
  !    Computes rhs of the 34 variables. Done with maple
  !
  !---------------------------------------------------------------------

#if defined GH_MHD || defined GH_FLOWER
  subroutine gh_calcrhs (rhs, u_pt, v_pt, dxu_pt, dyu_pt, dzu_pt, &
                         dxv_pt, dyv_pt, dzv_pt, w_pt, par)
#else
  subroutine calcrhs (rhs, u_pt, v_pt, dxu_pt, dyu_pt, dzu_pt, &
                      dxv_pt, dyv_pt, dzv_pt, w_pt, par)
#endif
    implicit none

    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    CCTK_REAL, dimension(NU_G), intent(inout) :: rhs
    CCTK_REAL, dimension(NU_G), intent(in)    :: dxu_pt,dyu_pt,dzu_pt
    CCTK_REAL, dimension(NU_G), intent(in)    :: u_pt
    CCTK_REAL, dimension(NV_G), intent(in)    :: v_pt,dxv_pt,dyv_pt,dzv_pt

    CCTK_REAL, dimension(NW), intent(in)     :: w_pt
    CCTK_REAL, dimension(NPAR),intent(inout) :: par

    CCTK_INT  :: i, iter
    logical   :: do_evolution

    ! Function arrays and time derivatives
    CCTK_REAL                             :: phir_dot, phic_dot, phim_dot, phin_dot
    CCTK_REAL                             :: pir_dot, pic_dot, pim_dot, pin_dot, G0_dot
    CCTK_REAL                             :: Fpir, Fpic, Fpim, Fpin, FG0
    CCTK_REAL, dimension(I_1:I_3)         :: Fdphim, Fdphin, Fdphir, Fdphic, FdH0
    CCTK_REAL, dimension(I_1:I_3)         :: dphim_dot, dphin_dot, dphir_dot, dphic_dot, dH0_dot
    CCTK_REAL, dimension(I_0:I_3)         :: H_dot, Duphim, Duphin, Duphir, Duphic
    CCTK_REAL, dimension(IS4_00:IS4_33)   :: g_dot, K_dot, FK
    CCTK_REAL, dimension(IS4_100:IS4_333) :: dg_dot, FD
        
    CCTK_REAL                             :: alp, alp2, deth, sqdetg
    CCTK_REAL                             :: phi2rc, phi2mn, trT, ttK, tZ, d0H0
    CCTK_REAL, dimension(I_0:I_3)         :: beta, tu, t, trG, Zd, Hu, ttD, tK
    CCTK_REAL, dimension(I_0:I_3)         :: DuH0
    CCTK_REAL, dimension(IS3_11:IS3_33)   :: h, huu
    CCTK_REAL, dimension(IN4_00:IN4_33)   :: tD
    CCTK_REAL, dimension(IS4_00:IS4_33)   :: guu, D0
    CCTK_REAL, dimension(IS4_000:IS4_333) :: Gudd, Du
    CCTK_REAL :: dx_H0e, dy_H0e, dz_H0e
    CCTK_REAL :: dphim0, dphin0, dphir0, dphic0
            
    ! other stuff
    CCTK_REAL          :: pi, sigma0, sigma1, sigma2, sf_m, sf_lam, sigma0_slope,sigma0_mu
    CCTK_REAL          :: chi1, chi2, chin, lchi1, lchi2, ltime, ltime1, ltime0, frac
    CCTK_INT           :: debug, rhs_type, gauge_type, damping_matter
    CCTK_INT           :: evolve_geometry, evolve_scalarfield ,psi4_analysis 
    logical, parameter :: ltrace   = .true.
    logical, parameter :: nancheck = .false.
    CCTK_INT           :: myid, proc_return_myid, m
    logical            :: isanan, nanfound
    ! Save code from lapse getting too small
    !    (not wanted when doing a critical search):
    logical, parameter :: rescue   = .true.

    CCTK_REAL	       :: fact, x, y, z, r, damping_factor, damping_rmax, decay,amp_phi
    CCTK_REAL          :: absm1, absm2, masstot, caltf, caltg, alphaadjust_decay

! travis here, adding rho to further damp tZ:
    CCTK_REAL          :: rho
    CCTK_INT           :: i_use_rho

    x = w_pt(H_XPHYS)
    y = w_pt(H_YPHYS)
    z = w_pt(H_ZPHYS)
    r = sqrt(x*x + y*y + z*z)
    !call mask_return_closestdist(x,y,z,r)

    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    myid     = proc_return_myid()
    nanfound = .false.
    if (ltrace) then
       write(0,*)myid,']  @@@ rhs_maple: begin (after initial includes)'
    end if
    if (nancheck) then
       !
       ! Checking for nans on input
       !
       do m = 1, NU_G
          if ( isanan(u_pt(m)) .or. isanan(dxu_pt(m)).or.isanan(dyu_pt(m)).or.isanan(dzu_pt(m)) ) then
             nanfound = .true.
          end if
       end do
       do m = 1, NV_G
          if ( isanan(v_pt(m)) .or. isanan(dxv_pt(m)).or.isanan(dyv_pt(m)).or.isanan(dzv_pt(m)) ) then
             nanfound = .true.
          end if
       end do
       do m = 1, NW
          if ( isanan(w_pt(m)) ) then
             nanfound = .true.
          end if
       end do
       if (nanfound) then
          write(0,*)myid,']  calc_rhs: Nan found on input:'
          do m = 1, NU_G
             write(0,*)myid,']  calc_rhs:      RHS(m):    ',m,RHS(m)
             write(0,*)myid,']  calc_rhs:      u_pt(m):   ',m,u_pt(m)
             write(0,*)myid,']  calc_rhs:      dxu_pt(m): ',m,dxu_pt(m)
             write(0,*)myid,']  calc_rhs:      dyu_pt(m): ',m,dyu_pt(m)
             write(0,*)myid,']  calc_rhs:      dzu_pt(m): ',m,dzu_pt(m)
          end do
          do m = 1, NV_G
             write(0,*)myid,']  calc_rhs:      v_pt(m):   ',m,v_pt(m)
             write(0,*)myid,']  calc_rhs:      dxv_pt(m): ',m,dxv_pt(m)
             write(0,*)myid,']  calc_rhs:      dyv_pt(m): ',m,dyv_pt(m)
             write(0,*)myid,']  calc_rhs:      dzv_pt(m): ',m,dzv_pt(m)
          end do
          do m = 1, NW
             write(0,*)myid,']  calc_rhs:      w_pt(m):   ',m,w_pt(m)
          end do
          call my_exit('rhs: Nan found on input')
       end if
    end if

    evolve_geometry    = par(P_EVOLVE_GEOMETRY)
    evolve_scalarfield = par(P_EVOLVE_SCALARFIELD)
    psi4_analysis      = par(P_PSI4_ANALYSIS)
    rhs_type   = par(P_RHS_TYPE)
    gauge_type = par(P_GAUGE_TYPE)
    damping_matter = par(P_DAMPING_MATTER)
    sigma0     = par(P_SIGMA0)
    sigma1     = par(P_SIGMA1)
    sigma2     = par(P_SIGMA2)
    sigma0_slope = par(P_SIGMA0_SLOPE)
    lchi1      = par(P_CHI1)
    lchi2      = par(P_CHI2)
    ltime      = par(P_TIME)
    ltime1     = par(P_CHIT1)
    ltime0     = par(P_CHIT0)
    chin       = par(P_CHIN)
    sf_m       = par(P_SF_MASS)
    sf_lam     = par(P_SF_LAMBDA)
    pi         = 3.1415926535897932d0
    iter       = nint(par(P_RK_ITER))
    damping_rmax   = par(P_DAMPING_RMAX)
    damping_factor = par(P_DAMPING_FACTOR)
    amp_phi = par(P_AMP_PHI)
    absm1 = par(P_BH1_mass)
    absm2 = par(P_BH2_mass)
    masstot = absm1+absm2
    alphaadjust_decay = par(P_ALPHAADJUST_DECAY)

    if ((psi4_analysis .GE. 1) .AND. (iter .EQ. 1)) then 
       do_evolution = .false.
    else
       do_evolution = .true.
    end if            

!--comment this line in order to reuse the last step of analysis    
   do_evolution = .true.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !----OPTION TO AVOID EVOLVING THE FIRST ITER ----------!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print*, "QUE PASA AQUI"
    print*, "QUE PASA AQUI"
    print*, "QUE PASA AQUI"
    print*, "QUE PASA AQUI"
    print*, "QUE PASA AQUI"
    print*, "QUE PASA AQUI"
              
    if (do_evolution) then

    ! smooth transition from t=0 to t=t1 
    chi1 = 0.d0
    chi2 = 0.d0

    if ((gauge_type.eq.3) .AND. (ltime.gt.ltime0)) then
      frac = (ltime -ltime0)
      caltg = (2.-exp(-frac/(17.5*masstot)))*(1.-exp(-frac**2/(15.*masstot)**2))
      caltf = caltg
      chi1 = lchi1*caltg
      chi2 = lchi2!*caltf
    end if    

    ! smooth transition from damping to factor*damping
    decay = 1.0d0
    if (r .GT. abs(damping_rmax)*masstot) then
      if(damping_rmax.gt.0) decay = (damping_rmax*masstot/r)**2.5
      if(damping_rmax.lt.0) decay = (exp(-(r+damping_rmax*masstot)**2 &
     &                                                       /(15.*masstot)**2)) 

      sigma0     = sigma0*decay
      sigma1     = sigma1*decay
      sigma2     = sigma2*decay
    end if 

      chi1   = chi1*exp(-r**2/(40.*masstot)**2)
      chi2   = chi2*exp(-r**2/(40.*masstot)**2)

    !lower the value of chi1, if time is sufficiently long
     if(amp_phi.lt.0.and.ltime.ge.abs(amp_phi)) then
!      chi1 = chi1 * exp((abs(amp_phi)-ltime)/5.)
      chi1 = chi1 * max(exp(0.25*(abs(amp_phi)-ltime)/(absm1+absm2)),0.001/(absm1+absm2))
     end if

     ! damping the constraints in the presence of matter
     ! 1 is the simplest option, there are more space for different alternatives
     sigma0_mu = 0.0d0
     if (int(damping_matter) .EQ. 1) then
        sigma0_mu = sigma0_slope*mod(ltime,2.)
     end if
     sigma0 = sigma0*(1.0d0 + sigma0_mu)


    ! initialize RHS to zero
    g_dot = 0.0d0;  K_dot = 0.0d0; dg_dot = 0.0d0
    H_dot = 0.0d0;  dH0_dot = 0.0d0; G0_dot = 0.0d0
    FK = 0.0d0; FD = 0.0d0
    pir_dot = 0.0d0;  pic_dot = 0.0d0;  pim_dot = 0.0d0
    phir_dot = 0.0d0; phic_dot = 0.0d0; phim_dot = 0.0d0

    RHS = 0.0d0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---- compute the 3+ metric components ----------------
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
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
    if(rescue .and. deth.lt.0.1) deth = 0.1
    huu11 = huu11/deth
    huu12 = huu12/deth
    huu13 = huu13/deth
    huu22 = huu22/deth
    huu23 = huu23/deth
    huu33 = huu33/deth
      
    fact = 1.0d0 
    b1 = g01*huu11 + g02*huu12 + g03*huu13
    b2 = g01*huu12 + g02*huu22 + g03*huu23
    b3 = g01*huu13 + g02*huu23 + g03*huu33
    alp2 = -g00 + b1**2*h11 + b2**2*h22 + 2*(b1*(b2*h12 + b3*h13) &
   &      + b2*b3*h23) + b3**2*h33
    if (rescue .and. alp2.le.0.01) then
!          write(0,*)myid,']  calc_rhs: alp2 <=0, Fixing it: alp2,dx',alp2,par(p_dx)
          alp2 = 0.01d0
          fact = 0.0d0
    end if
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

   !---normal to the 3D hypersurfaces---------------------
      
    tu0 = 1.0/alp; tu1 = -b1/alp; tu2 = -b2/alp; tu3 = -b3/alp 
    t0 = -alp; t1 = 0.0; t2 = 0.0; t3 = 0.0 

!we might need to adjust the decay with alpha to
!be aware that a BH might form
   if(alphaadjust_decay.gt.0) then
   sigma0 = sigma0/alp**(alphaadjust_decay)
   sigma1 = sigma1/alp**(alphaadjust_decay)
   sigma2 = sigma2/alp**(alphaadjust_decay)
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---- compute the auxiliary variables ----------------
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
        
   !---the matter terms -------------------------------    
                             
    trT = guu00*T00 + guu11*T11 + guu22*T22 &
     & + 2.0d0*(guu01*T01 +&
     & guu02*T02 + guu03*T03 + guu12*T12 &
     & + guu13*T13 + guu23*T23) +&
     & guu33*T33

   !---the auxiliary variables ---------------------------   
    Du100 = huu11*D100 + huu12*D200 + huu13*D300+0.d0
    Du101 = huu11*D101 + huu12*D201 + huu13*D301+0.d0
    Du102 = huu11*D102 + huu12*D202 + huu13*D302+0.d0
    Du103 = huu11*D103 + huu12*D203 + huu13*D303+0.d0
    Du111 = huu11*D111 + huu12*D211 + huu13*D311+0.d0
    Du112 = huu11*D112 + huu12*D212 + huu13*D312+0.d0
    Du113 = huu11*D113 + huu12*D213 + huu13*D313+0.d0
    Du122 = huu11*D122 + huu12*D222 + huu13*D322+0.d0
    Du123 = huu11*D123 + huu12*D223 + huu13*D323+0.d0
    Du133 = huu11*D133 + huu12*D233 + huu13*D333+0.d0
    Du200 = huu12*D100 + huu22*D200 + huu23*D300+0.d0
    Du201 = huu12*D101 + huu22*D201 + huu23*D301+0.d0
    Du202 = huu12*D102 + huu22*D202 + huu23*D302+0.d0
    Du203 = huu12*D103 + huu22*D203 + huu23*D303+0.d0
    Du211 = huu12*D111 + huu22*D211 + huu23*D311+0.d0
    Du212 = huu12*D112 + huu22*D212 + huu23*D312+0.d0
    Du213 = huu12*D113 + huu22*D213 + huu23*D313+0.d0
    Du222 = huu12*D122 + huu22*D222 + huu23*D322+0.d0
    Du223 = huu12*D123 + huu22*D223 + huu23*D323+0.d0
    Du233 = huu12*D133 + huu22*D233 + huu23*D333+0.d0
    Du300 = huu13*D100 + huu23*D200 + huu33*D300+0.d0
    Du301 = huu13*D101 + huu23*D201 + huu33*D301+0.d0
    Du302 = huu13*D102 + huu23*D202 + huu33*D302+0.d0
    Du303 = huu13*D103 + huu23*D203 + huu33*D303+0.d0
    Du311 = huu13*D111 + huu23*D211 + huu33*D311+0.d0
    Du312 = huu13*D112 + huu23*D212 + huu33*D312+0.d0
    Du313 = huu13*D113 + huu23*D213 + huu33*D313+0.d0
    Du322 = huu13*D122 + huu23*D222 + huu33*D322+0.d0
    Du323 = huu13*D123 + huu23*D223 + huu33*D323+0.d0
    Du333 = huu13*D133 + huu23*D233 + huu33*D333+0.d0      
   
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
   
    G000 = D000/2.
    G001 = D100/2.
    G002 = D200/2.
    G003 = D300/2.
    G011 = -D011/2. + D101
    G012 = (-D012 + D102 + D201)/2.
    G013 = (-D013 + D103 + D301)/2.
    G022 = -D022/2. + D202
    G023 = (-D023 + D203 + D302)/2.
    G033 = -D033/2. + D303
    G100 = D001 - D100/2.
    G101 = D011/2.
    G102 = (D012 - D102 + D201)/2.
    G103 = (D013 - D103 + D301)/2.
    G111 = D111/2.
    G112 = D211/2.
    G113 = D311/2.
    G122 = -D122/2. + D212
    G123 = (-D123 + D213 + D312)/2.
    G133 = -D133/2. + D313
    G200 = D002 - D200/2.
    G201 = (D012 + D102 - D201)/2.
    G202 = D022/2.
    G203 = (D023 - D203 + D302)/2.
    G211 = D112 - D211/2.
    G212 = D122/2.
    G213 = (D123 - D213 + D312)/2.
    G222 = D222/2.
    G223 = D322/2.
    G233 = -D233/2. + D323
    G300 = D003 - D300/2.
    G301 = (D013 + D103 - D301)/2.
    G302 = (D023 + D203 - D302)/2.
    G303 = D033/2.
    G311 = D113 - D311/2.
    G312 = (D123 + D213 - D312)/2.
    G313 = D133/2.
    G322 = D223 - D322/2.
    G323 = D233/2.
    G333 = D333/2.
    
    trG0 = G000*guu00 + G011*guu11 + G022*guu22 + 2* (G001*guu01 +&
     & G002*guu02 + G003*guu03 + G012*guu12 + G013*guu13 + G023*guu23) +&
     & G033*guu33
    trG1 = G100*guu00 + G111*guu11 + G122*guu22 + 2* (G101*guu01 +&
     & G102*guu02 + G103*guu03 + G112*guu12 + G113*guu13 + G123*guu23) +&
     & G133*guu33
    trG2 = G200*guu00 + G211*guu11 + G222*guu22 + 2* (G201*guu01 +&
     & G202*guu02 + G203*guu03 + G212*guu12 + G213*guu13 + G223*guu23) +&
     & G233*guu33
    trG3 = G300*guu00 + G311*guu11 + G322*guu22 + 2* (G301*guu01 +&
     & G302*guu02 + G303*guu03 + G312*guu12 + G313*guu13 + G323*guu23) +&
     & G333*guu33
     
    Z0 = (-H0 - trG0)/2.
    Z1 = (-H1 - trG1)/2.
    Z2 = (-H2 - trG2)/2.
    Z3 = (-H3 - trG3)/2.
    
    Hu0 = guu00*H0 + guu01*H1 + guu02*H2 + guu03*H3
    Hu1 = guu01*H0 + guu11*H1 + guu12*H2 + guu13*H3
    Hu2 = guu02*H0 + guu12*H1 + guu22*H2 + guu23*H3
    Hu3 = guu03*H0 + guu13*H1 + guu23*H2 + guu33*H3
    
    !contractions with tu
    tD00 = tu0*D000 + tu1*D001 + tu2*D002 + tu3*D003
    tD01 = tu0*D001 + tu1*D011 + tu2*D012 + tu3*D013
    tD02 = tu0*D002 + tu1*D012 + tu2*D022 + tu3*D023
    tD03 = tu0*D003 + tu1*D013 + tu2*D023 + tu3*D033
    tD10 = tu0*D100 + tu1*D101 + tu2*D102 + tu3*D103
    tD11 = tu0*D101 + tu1*D111 + tu2*D112 + tu3*D113
    tD12 = tu0*D102 + tu1*D112 + tu2*D122 + tu3*D123
    tD13 = tu0*D103 + tu1*D113 + tu2*D123 + tu3*D133
    tD20 = tu0*D200 + tu1*D201 + tu2*D202 + tu3*D203
    tD21 = tu0*D201 + tu1*D211 + tu2*D212 + tu3*D213
    tD22 = tu0*D202 + tu1*D212 + tu2*D222 + tu3*D223
    tD23 = tu0*D203 + tu1*D213 + tu2*D223 + tu3*D233
    tD30 = tu0*D300 + tu1*D301 + tu2*D302 + tu3*D303
    tD31 = tu0*D301 + tu1*D311 + tu2*D312 + tu3*D313
    tD32 = tu0*D302 + tu1*D312 + tu2*D322 + tu3*D323
    tD33 = tu0*D303 + tu1*D313 + tu2*D323 + tu3*D333
    
    ttD0 = tD00*tu0 + tD01*tu1 + tD02*tu2 + tD03*tu3
    ttD1 = tD10*tu0 + tD11*tu1 + tD12*tu2 + tD13*tu3
    ttD2 = tD20*tu0 + tD21*tu1 + tD22*tu2 + tD23*tu3
    ttD3 = tD30*tu0 + tD31*tu1 + tD32*tu2 + tD33*tu3
    
    tK0 = K00*tu0 + K01*tu1 + K02*tu2 + K03*tu3
    tK1 = K01*tu0 + K11*tu1 + K12*tu2 + K13*tu3
    tK2 = K02*tu0 + K12*tu1 + K22*tu2 + K23*tu3
    tK3 = K03*tu0 + K13*tu1 + K23*tu2 + K33*tu3
    
    ttK = tK0*tu0 + tK1*tu1 + tK2*tu2 + tK3*tu3
    tZ = tu0*Z0 + tu1*Z1 + tu2*Z2 + tu3*Z3

    rho = par(P_SIGMA0_RHO)
    tZ = tZ*(1.0d0+rho)


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !----OPTION TO AVOID EVOLVING THE SPACETIME  ----------!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   if (evolve_geometry .EQ. 1) then    
           
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---- the evolution of the gauge  --------------------
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! H are either sources or independent quantities with no evolution
   if ((gauge_type .EQ. 0) .OR. (gauge_type .EQ. 1)) then

      H_dot(I_0)= 0.0d0; H_dot(I_1)= 0.0d0; H_dot(I_2)= 0.0d0; H_dot(I_3)= 0.0d0     
      dx_H0e = dx_H0; dy_H0e = dy_H0; dz_H0e = dz_H0

   ! H are independent quantities, the HAKE gauge     
   elseif (gauge_type .EQ. 2) then
                 
      write(0,*)myid,' the HAKE gauge does not work '
      write(0,*)myid,' change gauge_type to something different from 2 !!'      
      call my_exit('rhs: gauge_type not allowed')       
      
   elseif (gauge_type .EQ. 3) then
     
      H_dot(I_0) = b1*d1H0 + b2*d2H0 + b3*d3H0 - alp*G0 
      H_dot(I_1) = 0.0d0; H_dot(I_2) = 0.0d0; H_dot(I_3) = 0.0d0   
      dx_H0e = d1H0; dy_H0e = d2H0; dz_H0e = d3H0         
                             
      d0H0 = H_dot(I_0)
      du0H0 = d0H0*guu00 + d1H0*guu01 + d2H0*guu02 + d3H0*guu03
      du1H0 = d0H0*guu01 + d1H0*guu11 + d2H0*guu12 + d3H0*guu13
      du2H0 = d0H0*guu02 + d1H0*guu12 + d2H0*guu22 + d3H0*guu23
      du3H0 = d0H0*guu03 + d1H0*guu13 + d2H0*guu23 + d3H0*guu33

      FG0 = -(b1*dx_G0 + b2*dy_G0 + b3*dz_G0) &
   &       + alp*(huu11*dx_d1H0 + huu12*dx_d2H0 + huu13*dx_d3H0 &
   &            + huu12*dy_d1H0 + huu22*dy_d2H0 + huu23*dy_d3H0 &
   &            + huu13*dz_d1H0 + huu23*dz_d2H0 + huu33*dz_d3H0)&
   &       + sigma2*(b1*dx_H0 + b2*dy_H0 + b3*dz_H0) 
    
      FdH0(I_1) = - (b1*dx_d1H0 + b2*dy_d1H0 + b3*dz_d1H0) & 
   &              + alp*dx_G0 - alp*sigma2*dx_H0
      FdH0(I_2) = - (b1*dx_d2H0 + b2*dy_d2H0 + b3*dz_d2H0) & 
   &              + alp*dy_G0 - alp*sigma2*dy_H0
      FdH0(I_3) = - (b1*dx_d3H0 + b2*dy_d3H0 + b3*dz_d3H0) & 
   &              + alp*dz_G0 - alp*sigma2*dz_H0
      

      G0_dot = (b1*d1H0 + b2*d2H0 + b3*d3H0)*sigma2 + &
     & alp*(-((d1H0*huu11 + d2H0*huu12 + d3H0*huu13)*tK1) - (d1H0*huu12 &
     & + d2H0*huu22 + d3H0*huu23)* tK2 - (d1H0*huu13 + d2H0*huu23 + &
     & d3H0*huu33)*tK3 + du0H0*trG0 + du1H0*trG1 + du2H0*trG2 + &
     & du3H0*trG3 - G0*(chi2 + ttK/2.))  - chi1*(alp -1.0)*(alp**(1.0-chin))
      dH0_dot(I_1) = alp*((d2H0*huu12 + d3H0*huu13)*tD11 + (d2H0*huu22 &
     & + d3H0*huu23)*tD12 + (d2H0*huu23 + d3H0*huu33)*tD13 + &
     & d1H0*(-sigma2 + huu11*tD11 + huu12*tD12 + huu13*tD13) + &
     & (G0*ttD1)/2.)
      dH0_dot(I_2) = alp*((d1H0*huu11 + d3H0*huu13)*tD21 + (d1H0*huu12 &
     & + d3H0*huu23)*tD22 + (d1H0*huu13 + d3H0*huu33)*tD23 + &
     & d2H0*(-sigma2 + huu12*tD21 + huu22*tD22 + huu23*tD23) + &
     & (G0*ttD2)/2.)
      dH0_dot(I_3) = alp*((d1H0*huu11 + d2H0*huu12)*tD31 + (d1H0*huu12 &
     & + d2H0*huu22)*tD32 + (d1H0*huu13 + d2H0*huu23)*tD33 + &
     & d3H0*(-sigma2 + huu13*tD31 + huu23*tD32 + huu33*tD33) + &
     & (G0*ttD3)/2.)
                                                  
   end if   ! --- gauge_type ---------------------------------

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---- the principal part of the Ks --------------------
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      
    FK(IS4_00) = - (b1*dx_K00 + b2*dy_K00 + b3*dz_K00) & 
   &       + sigma2*(b1*dx_g00 + b2*dy_g00 + b3*dz_g00) &
   &       + alp*(huu11*dx_D100 + huu12*dx_D200 + huu13*dx_D300 &
   &            + huu12*dy_D100 + huu22*dy_D200 + huu23*dy_D300 &
   &            + huu13*dz_D100 + huu23*dz_D200 + huu33*dz_D300) &
   &       + alp*H_dot(I_0) + alp*H_dot(I_0)
     
    FK(IS4_01) = - (b1*dx_K01 + b2*dy_K01 + b3*dz_K01) & 
   &       + sigma2*(b1*dx_g01 + b2*dy_g01 + b3*dz_g01) &
   &       + alp*(huu11*dx_D101 + huu12*dx_D201 + huu13*dx_D301 &
   &            + huu12*dy_D101 + huu22*dy_D201 + huu23*dy_D301 &
   &            + huu13*dz_D101 + huu23*dz_D201 + huu33*dz_D301)&
   &       + alp*dx_H0e + alp*H_dot(I_1)

    FK(IS4_02) = - (b1*dx_K02 + b2*dy_K02 + b3*dz_K02) & 
   &       + sigma2*(b1*dx_g02 + b2*dy_g02 + b3*dz_g02) &
   &       + alp*(huu11*dx_D102 + huu12*dx_D202 + huu13*dx_D302 &
   &            + huu12*dy_D102 + huu22*dy_D202 + huu23*dy_D302 &
   &            + huu13*dz_D102 + huu23*dz_D202 + huu33*dz_D302)&
   &       + alp*dy_H0e + alp*H_dot(I_2) 

    FK(IS4_03) = - (b1*dx_K03 + b2*dy_K03 + b3*dz_K03) & 
   &       + sigma2*(b1*dx_g03 + b2*dy_g03 + b3*dz_g03) &
   &       + alp*(huu11*dx_D103 + huu12*dx_D203 + huu13*dx_D303 &
   &            + huu12*dy_D103 + huu22*dy_D203 + huu23*dy_D303 &
   &            + huu13*dz_D103 + huu23*dz_D203 + huu33*dz_D303)&
   &       + alp*dz_H0e + alp*H_dot(I_3)

    FK(IS4_11) = - (b1*dx_K11 + b2*dy_K11 + b3*dz_K11) & 
   &       + sigma2*(b1*dx_g11 + b2*dy_g11 + b3*dz_g11) &
   &       + alp*(huu11*dx_D111 + huu12*dx_D211 + huu13*dx_D311 &
   &            + huu12*dy_D111 + huu22*dy_D211 + huu23*dy_D311 &
   &            + huu13*dz_D111 + huu23*dz_D211 + huu33*dz_D311)&
   &       + alp*dx_H1 + alp*dx_H1

    FK(IS4_12) = - (b1*dx_K12 + b2*dy_K12 + b3*dz_K12) & 
   &       + sigma2*(b1*dx_g12 + b2*dy_g12 + b3*dz_g12) &
   &       + alp*(huu11*dx_D112 + huu12*dx_D212 + huu13*dx_D312 &
   &            + huu12*dy_D112 + huu22*dy_D212 + huu23*dy_D312 &
   &            + huu13*dz_D112 + huu23*dz_D212 + huu33*dz_D312)&
   &       + alp*dx_H2 + alp*dy_H1

    FK(IS4_13) = - (b1*dx_K13 + b2*dy_K13 + b3*dz_K13) & 
   &       + sigma2*(b1*dx_g13 + b2*dy_g13 + b3*dz_g13) &
   &       + alp*(huu11*dx_D113 + huu12*dx_D213 + huu13*dx_D313 &
   &            + huu12*dy_D113 + huu22*dy_D213 + huu23*dy_D313 &
   &            + huu13*dz_D113 + huu23*dz_D213 + huu33*dz_D313)&
   &       + alp*dz_H1 + alp*dx_H3

    FK(IS4_22) = - (b1*dx_K22 + b2*dy_K22 + b3*dz_K22) & 
   &       + sigma2*(b1*dx_g22 + b2*dy_g22 + b3*dz_g22) &
   &       + alp*(huu11*dx_D122 + huu12*dx_D222 + huu13*dx_D322 &
   &            + huu12*dy_D122 + huu22*dy_D222 + huu23*dy_D322 &
   &            + huu13*dz_D122 + huu23*dz_D222 + huu33*dz_D322)&
   &       + alp*dy_H2 + alp*dy_H2

    FK(IS4_23) = - (b1*dx_K23 + b2*dy_K23 + b3*dz_K23) & 
   &       + sigma2*(b1*dx_g23 + b2*dy_g23 + b3*dz_g23) &
   &       + alp*(huu11*dx_D123 + huu12*dx_D223 + huu13*dx_D323 &
   &            + huu12*dy_D123 + huu22*dy_D223 + huu23*dy_D323 &
   &            + huu13*dz_D123 + huu23*dz_D223 + huu33*dz_D323)&
   &       + alp*dy_H3 + alp*dz_H2
                                        
    FK(IS4_33) = - (b1*dx_K33 + b2*dy_K33 + b3*dz_K33) & 
   &       + sigma2*(b1*dx_g33 + b2*dy_g33 + b3*dz_g33) &
   &       + alp*(huu11*dx_D133 + huu12*dx_D233 + huu13*dx_D333 &
   &            + huu12*dy_D133 + huu22*dy_D233 + huu23*dy_D333 &
   &            + huu13*dz_D133 + huu23*dz_D233 + huu33*dz_D333)&
   &       + alp*dz_H3 + alp*dz_H3
     
              
!   the principal part of the Diab
    FD(IS4_100) = - (b1*dx_D100 + b2*dy_D100 + b3*dz_D100) & 
   &              + alp*dx_K00 - alp*sigma2*dx_g00
    FD(IS4_101) = - (b1*dx_D101 + b2*dy_D101 + b3*dz_D101) & 
   &              + alp*dx_K01 - alp*sigma2*dx_g01
    FD(IS4_102) = - (b1*dx_D102 + b2*dy_D102 + b3*dz_D102) & 
   &              + alp*dx_K02 - alp*sigma2*dx_g02
    FD(IS4_103) = - (b1*dx_D103 + b2*dy_D103 + b3*dz_D103) & 
   &              + alp*dx_K03 - alp*sigma2*dx_g03
    FD(IS4_111) = - (b1*dx_D111 + b2*dy_D111 + b3*dz_D111) & 
   &              + alp*dx_K11 - alp*sigma2*dx_g11
    FD(IS4_112) = - (b1*dx_D112 + b2*dy_D112 + b3*dz_D112) & 
   &              + alp*dx_K12 - alp*sigma2*dx_g12
    FD(IS4_113) = - (b1*dx_D113 + b2*dy_D113 + b3*dz_D113) & 
   &              + alp*dx_K13 - alp*sigma2*dx_g13
    FD(IS4_122) = - (b1*dx_D122 + b2*dy_D122 + b3*dz_D122) & 
   &              + alp*dx_K22 - alp*sigma2*dx_g22
    FD(IS4_123) = - (b1*dx_D123 + b2*dy_D123 + b3*dz_D123) & 
   &              + alp*dx_K23 - alp*sigma2*dx_g23
    FD(IS4_133) = - (b1*dx_D133 + b2*dy_D133 + b3*dz_D133) & 
   &              + alp*dx_K33 - alp*sigma2*dx_g33
 
    FD(IS4_200) = - (b1*dx_D200 + b2*dy_D200 + b3*dz_D200) & 
   &              + alp*dy_K00 - alp*sigma2*dy_g00
    FD(IS4_201) = - (b1*dx_D201 + b2*dy_D201 + b3*dz_D201) & 
   &              + alp*dy_K01 - alp*sigma2*dy_g01
    FD(IS4_202) = - (b1*dx_D202 + b2*dy_D202 + b3*dz_D202) & 
   &              + alp*dy_K02 - alp*sigma2*dy_g02
    FD(IS4_203) = - (b1*dx_D203 + b2*dy_D203 + b3*dz_D203) & 
   &              + alp*dy_K03 - alp*sigma2*dy_g03
    FD(IS4_211) = - (b1*dx_D211 + b2*dy_D211 + b3*dz_D211) & 
   &              + alp*dy_K11 - alp*sigma2*dy_g11
    FD(IS4_212) = - (b1*dx_D212 + b2*dy_D212 + b3*dz_D212) & 
   &              + alp*dy_K12 - alp*sigma2*dy_g12
    FD(IS4_213) = - (b1*dx_D213 + b2*dy_D213 + b3*dz_D213) & 
   &              + alp*dy_K13 - alp*sigma2*dy_g13
    FD(IS4_222) = - (b1*dx_D222 + b2*dy_D222 + b3*dz_D222) & 
   &              + alp*dy_K22 - alp*sigma2*dy_g22
    FD(IS4_223) = - (b1*dx_D223 + b2*dy_D223 + b3*dz_D223) & 
   &              + alp*dy_K23 - alp*sigma2*dy_g23
    FD(IS4_233) = - (b1*dx_D233 + b2*dy_D233 + b3*dz_D233) & 
   &              + alp*dy_K33 - alp*sigma2*dy_g33
  
    FD(IS4_300) = - (b1*dx_D300 + b2*dy_D300 + b3*dz_D300) & 
   &              + alp*dz_K00 - alp*sigma2*dz_g00
    FD(IS4_301) = - (b1*dx_D301 + b2*dy_D301 + b3*dz_D301) & 
   &              + alp*dz_K01 - alp*sigma2*dz_g01
    FD(IS4_302) = - (b1*dx_D302 + b2*dy_D302 + b3*dz_D302) & 
   &              + alp*dz_K02 - alp*sigma2*dz_g02
    FD(IS4_303) = - (b1*dx_D303 + b2*dy_D303 + b3*dz_D303) & 
   &              + alp*dz_K03 - alp*sigma2*dz_g03
    FD(IS4_311) = - (b1*dx_D311 + b2*dy_D311 + b3*dz_D311) & 
   &              + alp*dz_K11 - alp*sigma2*dz_g11
    FD(IS4_312) = - (b1*dx_D312 + b2*dy_D312 + b3*dz_D312) & 
   &              + alp*dz_K12 - alp*sigma2*dz_g12
    FD(IS4_313) = - (b1*dx_D313 + b2*dy_D313 + b3*dz_D313) & 
   &              + alp*dz_K13 - alp*sigma2*dz_g13
    FD(IS4_322) = - (b1*dx_D322 + b2*dy_D322 + b3*dz_D322) & 
   &              + alp*dz_K22 - alp*sigma2*dz_g22
    FD(IS4_323) = - (b1*dx_D323 + b2*dy_D323 + b3*dz_D323) & 
   &              + alp*dz_K23 - alp*sigma2*dz_g23
    FD(IS4_333) = - (b1*dx_D333 + b2*dy_D333 + b3*dz_D333) & 
   &              + alp*dz_K33 - alp*sigma2*dz_g33

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---- the sources terms -------------------------------
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   
    g_dot(IS4_00) = -(alp*K00) + b1*D100 + b2*D200 + b3*D300
    g_dot(IS4_01) = -(alp*K01) + b1*D101 + b2*D201 + b3*D301
    g_dot(IS4_02) = -(alp*K02) + b1*D102 + b2*D202 + b3*D302
    g_dot(IS4_03) = -(alp*K03) + b1*D103 + b2*D203 + b3*D303
    g_dot(IS4_11) = -(alp*K11) + b1*D111 + b2*D211 + b3*D311
    g_dot(IS4_12) = -(alp*K12) + b1*D112 + b2*D212 + b3*D312
    g_dot(IS4_13) = -(alp*K13) + b1*D113 + b2*D213 + b3*D313
    g_dot(IS4_22) = -(alp*K22) + b1*D122 + b2*D222 + b3*D322
    g_dot(IS4_23) = -(alp*K23) + b1*D123 + b2*D223 + b3*D323
    g_dot(IS4_33) = -(alp*K33) + b1*D133 + b2*D233 + b3*D333

    K_dot(IS4_00) = sigma2*(b1*D100 + b2*D200 + b3*D300) +&
     & alp*(-(G011*(8*G003*guu01*guu13 + 4*(G022*guu12**2 +&
     & G033*guu13**2))) - G033*(4*G022*guu23**2 + 8*G003*guu03*guu33) -&
     & 8*(G000*guu00* (G001*guu01 + G002*guu02 + G003*guu03) +&
     & G001*((G002*guu00 + G023*guu03)*guu12 + (G003*guu00 +&
     & G033*guu03)* guu13 + guu02*(G012*guu11 + G022*guu12 +&
     & G023*guu13)) + (G033*(G002*guu03 + G012*guu13) + G023*&
     & (G012*guu12 + G013*guu13))* guu23 + guu03* (G002*(G003*guu02 +&
     & G013*guu12) + G013* (G000*guu01 + G001*guu11 + G003*guu13) +&
     & G003*G023*guu23) + G002*(G022*guu02*guu22 + G003*guu00*guu23) +&
     & G012*((G001*guu01 + G002*guu02 + G011*guu11)* guu12 +&
     & guu13*(G003*guu02 + G023*guu22) + G013*guu11*guu23) +&
     & guu02*((G003*G022 + G002*G023)*guu23 + G003*G023*guu33) +&
     & guu12*(G022* (G012*guu22 + G013*guu23) + G023*(G011*guu13 +&
     & G013*guu33)) + guu01*((G001*G002 + G000*G012)*guu02 + G001*&
     & (G003*guu03 + G011*guu11 + G013*guu13) + G002* (G011*guu12 +&
     & G012*guu22) + G003*(G012*guu23 + G013*guu33)) +&
     & G013*(G002*guu01*guu23 + guu13*(G011*guu11 + G012*guu12 +&
     & G033*guu33))) + G000*(-8*G023*guu02*guu03 - 4*(G022*guu02**2 +&
     & G033*guu03**2) + 2*Hu0) - 2*(G000**2*guu00**2 + G011**2*guu11**2&
     & + G022**2*guu22**2 + G033**2*guu33**2 + guu00*K00**2 +&
     & guu11*K01**2) - guu22*(8*G023* (G002*guu03 + G022*guu23) +&
     & 4*(G012**2*guu11 + G023**2*guu33) + 2*K02**2) -&
     & guu33*(4*G013**2*guu11 + 8*G023*G033*guu23 + 2*K03**2) +&
     & pi*(-16*T00 + 8*g00*trT) - (K00*ttK)/2. - tK1*Du100 - tK2*Du200 -&
     & tK3*Du300 + 2*(G100*Hu1 + G200*Hu2 + G300*Hu3 + g00*sigma0*tZ +&
     & (guu00*D100 + guu01*D101 + guu02*D102 + guu03*D103)*Du100 +&
     & (guu01*D100 + guu11*D101 + guu12*D102 + guu13*D103)*Du101 +&
     & (guu02*D100 + guu12*D101 + guu22*D102 + guu23*D103)*Du102 +&
     & (guu03*D100 + guu13*D101 + guu23*D102 + guu33*D103)*Du103 +&
     & (guu00*D200 + guu01*D201 + guu02*D202 + guu03*D203)*Du200 +&
     & (guu01*D200 + guu11*D201 + guu12*D202 + guu13*D203)*Du201 +&
     & (guu02*D200 + guu12*D201 + guu22*D202 + guu23*D203)*Du202 +&
     & (guu03*D200 + guu13*D201 + guu23*D202 + guu33*D203)*Du203 +&
     & (guu00*D300 + guu01*D301 + guu02*D302 + guu03*D303)*Du300 +&
     & (guu01*D300 + guu11*D301 + guu12*D302 + guu13*D303)*Du301 +&
     & (guu02*D300 + guu12*D301 + guu22*D302 + guu23*D303)*Du302 +&
     & (guu03*D300 + guu13*D301 + guu23*D302 + guu33*D303)*Du303) -&
     & 4*((G001**2 + G000*G011)* guu01**2 + G001**2*guu00*guu11 +&
     & G012**2*guu12**2 + G013**2*guu13**2 + G002**2*(guu02**2 +&
     & guu00*guu22) + G023**2*guu23**2 + G003**2*(guu03**2 +&
     & guu00*guu33) + guu23*K02*K03 + K00*(guu01*K01 + guu02*K02 +&
     & guu03*K03) + K01*(guu12*K02 + guu13*K03) + sigma0*t0*Z0))+0.d0
    
    K_dot(IS4_01) = sigma2*(b1*D101 + b2*D201 + b3*D301) +&
     & alp*(-((4*G002*G102 + 2*G000*G122)* guu02**2) - G100*(4*guu00*&
     & (G001*guu01 + G002*guu02) + 2*(G000*guu00**2 + G033*guu03**2)) -&
     & (4*G012*G112 + 2*G011*G122)* guu12**2 - (4*G013*G113 +&
     & 2*G011*G133)*guu13**2 - (4*G023*G123 + 2*G022*G133)* guu23**2 -&
     & 4* (G101*(G001*guu01**2 + guu01*(G000*guu00 + G002*guu02 +&
     & G003*guu03) + guu03*(G013*guu11 + G023*guu12)) +&
     & G001*((G111*guu01 + G112*guu02 + G113*guu03)* guu11 +&
     & G102*guu00*guu12 + G113*guu01*guu13) + G002*((G112*guu01 +&
     & G122*guu02)*guu22 + G113*guu01*guu23) + (G013*(G113*guu11 +&
     & G123*guu12) + G133* (G013*guu13 + G023*guu23))* guu33 + guu01*&
     & ((G001*G102 + G000*G112)* guu02 + (G013*G100 + G001*G103 +&
     & G000*G113)*guu03 + (G011*G102 + G002*G111)* guu12 + (G011*G103 +&
     & G003*G111)* guu13 + G101* (G011*guu11 + G012*guu12 + G013*guu13)&
     & + (G012*G103 + G003*G112)* guu23 + G102* (G012*guu22 +&
     & G013*guu23) + G013*G103*guu33) + guu02*((G002*G103 +&
     & G000*G123)*guu03 + G100* (G012*guu01 + G023*guu03) + (G002*G112 +&
     & G001*G122)* guu12 + (G012*G103 + G003*G112)* guu13 + G101*&
     & (G012*guu11 + G022*guu12 + G023*guu13) + (G022*G103 + G003*G122)*&
     & guu23 + G102* (G000*guu00 + G003*guu03 + G022*guu22 + G023*guu23)&
     & + G023*G103*guu33) + guu03*((G002*G113 + G001*G123)*guu12 +&
     & (G033*G101 + G003*G113 + G001*G133)*guu13 + (G003*G123 +&
     & G002*G133)* guu23 + G102* (G023*guu22 + G033*guu23) +&
     & G033*G103*guu33) + guu12*(G102* (G012*guu02 + G013*guu03) +&
     & (G012*G113 + G011*G123)* guu13 + G111* (G012*guu11 + G023*guu13)&
     & + (G013*G122 + G012*G123)* guu23 + G112* (G001*guu01 + G011*guu11&
     & + G022*guu22 + G023*guu23) + G023*G113*guu33) +&
     & G003*(G103*guu03**2 + guu00*(G101*guu13 + G102*guu23) +&
     & (G113*guu01 + G123*guu02 + G133*guu03)*guu33) + guu00*((G003*G100&
     & + G000*G103)*guu03 + G101* (G001*guu11 + G002*guu12) +&
     & G002*G102*guu22 + G103*(G001*guu13 + G002*guu23 + G003*guu33)) +&
     & guu22*(G012*(G112*guu11 + G122*guu12) + G023*G122*guu23 +&
     & G123*(G002*guu03 + G012*guu13 + G023*guu33)) +&
     & guu13*(G013*(G103*guu03 + G111*guu11 + G112*guu12) +&
     & G012*G133*guu23 + G123* (G001*guu02 + G013*guu23) + G112*&
     & (G023*guu22 + G033*guu23) + G113*(G011*guu11 + G033*guu33)) +&
     & guu23*((G013*G112 + G012*G113)*guu11 + G022*G113*guu12 +&
     & G023*(G103*guu03 + G113*guu13) + G123*(G002*guu02 + G022*guu22 +&
     & G033*guu33))) + pi*(-16*T01 + 8*g01*trT) - (K01*ttK)/2. -&
     & tK1*Du101 - tK2*Du201 - tK3*Du301 - 2*((G011*G100 + G000*G111)*&
     & guu01**2 + G111*(G011*guu11**2 + G033*guu13**2) +&
     & G022*(G100*guu02**2 + G111*guu12**2 + G122*guu22**2) +&
     & G033*G122*guu23**2 + G133*(G000*guu03**2 + G033*guu33**2) +&
     & (guu12*K02 + guu13*K03)*K11 + guu01*(K01**2 + K00*K11) +&
     & (guu02*K00 + guu22*K02 + guu23*K03)* K12 + (guu03*K00 + guu23*K02&
     & + guu33*K03)*K13 + K01*(guu00*K00 + guu02*K02 + guu03*K03 +&
     & guu11*K11 + guu12*K12 + guu13*K13) + sigma0*t1*Z0) + 2*(G001*Hu0&
     & + G101*Hu1 + G201*Hu2 + G301*Hu3 + (guu00*D100 + guu01*D101 +&
     & guu02*D102 + guu03*D103)*Du101 + (guu01*D100 + guu11*D101 +&
     & guu12*D102 + guu13*D103)*Du111 + (guu02*D100 + guu12*D101 +&
     & guu22*D102 + guu23*D103)*Du112 + (guu03*D100 + guu13*D101 +&
     & guu23*D102 + guu33*D103)*Du113 + (guu00*D200 + guu01*D201 +&
     & guu02*D202 + guu03*D203)*Du201 + (guu01*D200 + guu11*D201 +&
     & guu12*D202 + guu13*D203)*Du211 + (guu02*D200 + guu12*D201 +&
     & guu22*D202 + guu23*D203)*Du212 + (guu03*D200 + guu13*D201 +&
     & guu23*D202 + guu33*D203)*Du213 + (guu00*D300 + guu01*D301 +&
     & guu02*D302 + guu03*D303)*Du301 + (guu01*D300 + guu11*D301 +&
     & guu12*D302 + guu13*D303)*Du311 + (guu02*D300 + guu12*D301 +&
     & guu22*D302 + guu23*D303)*Du312 + (guu03*D300 + guu13*D301 +&
     & guu23*D302 + guu33*D303)*Du313 + sigma0*(g01*tZ - t0*Z1)))+0.d0
    
    K_dot(IS4_02) = sigma2*(b1*D102 + b2*D202 + b3*D302) +&
     & alp*(-((4*G002*G202 + 2*G000*G222)* guu02**2) - G200*(4*guu00*&
     & (G001*guu01 + G002*guu02) + 2*(G000*guu00**2 + G033*guu03**2)) -&
     & (4*G012*G212 + 2*G011*G222)* guu12**2 - (4*G013*G213 +&
     & 2*G011*G233)*guu13**2 - (4*G023*G223 + 2*G022*G233)* guu23**2 -&
     & 4* (G201*(G001*guu01**2 + guu01*(G000*guu00 + G002*guu02 +&
     & G003*guu03) + guu03*(G013*guu11 + G023*guu12)) +&
     & G001*((G211*guu01 + G212*guu02 + G213*guu03)* guu11 +&
     & G202*guu00*guu12 + G213*guu01*guu13) + G002*((G212*guu01 +&
     & G222*guu02)*guu22 + G213*guu01*guu23) + (G013*(G213*guu11 +&
     & G223*guu12) + G233* (G013*guu13 + G023*guu23))* guu33 + guu01*&
     & ((G001*G202 + G000*G212)* guu02 + (G013*G200 + G001*G203 +&
     & G000*G213)*guu03 + (G011*G202 + G002*G211)* guu12 + (G011*G203 +&
     & G003*G211)* guu13 + G201* (G011*guu11 + G012*guu12 + G013*guu13)&
     & + (G012*G203 + G003*G212)* guu23 + G202* (G012*guu22 +&
     & G013*guu23) + G013*G203*guu33) + guu02*((G002*G203 +&
     & G000*G223)*guu03 + G200* (G012*guu01 + G023*guu03) + (G002*G212 +&
     & G001*G222)* guu12 + (G012*G203 + G003*G212)* guu13 + G201*&
     & (G012*guu11 + G022*guu12 + G023*guu13) + (G022*G203 + G003*G222)*&
     & guu23 + G202* (G000*guu00 + G003*guu03 + G022*guu22 + G023*guu23)&
     & + G023*G203*guu33) + guu03*((G002*G213 + G001*G223)*guu12 +&
     & (G033*G201 + G003*G213 + G001*G233)*guu13 + (G003*G223 +&
     & G002*G233)* guu23 + G202* (G023*guu22 + G033*guu23) +&
     & G033*G203*guu33) + guu12*(G202* (G012*guu02 + G013*guu03) +&
     & (G012*G213 + G011*G223)* guu13 + G211* (G012*guu11 + G023*guu13)&
     & + (G013*G222 + G012*G223)* guu23 + G212* (G001*guu01 + G011*guu11&
     & + G022*guu22 + G023*guu23) + G023*G213*guu33) +&
     & G003*(G203*guu03**2 + guu00*(G201*guu13 + G202*guu23) +&
     & (G213*guu01 + G223*guu02 + G233*guu03)*guu33) + guu00*((G003*G200&
     & + G000*G203)*guu03 + G201* (G001*guu11 + G002*guu12) +&
     & G002*G202*guu22 + G203*(G001*guu13 + G002*guu23 + G003*guu33)) +&
     & guu22*(G012*(G212*guu11 + G222*guu12) + G023*G222*guu23 +&
     & G223*(G002*guu03 + G012*guu13 + G023*guu33)) +&
     & guu13*(G013*(G203*guu03 + G211*guu11 + G212*guu12) +&
     & G012*G233*guu23 + G223* (G001*guu02 + G013*guu23) + G212*&
     & (G023*guu22 + G033*guu23) + G213*(G011*guu11 + G033*guu33)) +&
     & guu23*((G013*G212 + G012*G213)*guu11 + G022*G213*guu12 +&
     & G023*(G203*guu03 + G213*guu13) + G223*(G002*guu02 + G022*guu22 +&
     & G033*guu33))) + pi*(-16*T02 + 8*g02*trT) - (K02*ttK)/2. -&
     & tK1*Du102 - tK2*Du202 - tK3*Du302 - 2*((G011*G200 + G000*G211)*&
     & guu01**2 + G211*(G011*guu11**2 + G033*guu13**2) +&
     & G022*(G200*guu02**2 + G211*guu12**2 + G222*guu22**2) +&
     & G033*G222*guu23**2 + G233*(G000*guu03**2 + G033*guu33**2) +&
     & (guu01*K00 + guu11*K01 + guu13*K03)* K12 + (guu12*K01 +&
     & guu23*K03)*K22 + guu02*(K02**2 + K00*K22) + (guu03*K00 +&
     & guu13*K01 + guu33*K03)* K23 + K02*(guu00*K00 + guu01*K01 +&
     & guu03*K03 + guu12*K12 + guu22*K22 + guu23*K23) + sigma0*t2*Z0) +&
     & 2*(G002*Hu0 + G102*Hu1 + G202*Hu2 + G302*Hu3 + (guu00*D100 +&
     & guu01*D101 + guu02*D102 + guu03*D103)*Du102 + (guu01*D100 +&
     & guu11*D101 + guu12*D102 + guu13*D103)*Du112 + (guu02*D100 +&
     & guu12*D101 + guu22*D102 + guu23*D103)*Du122 + (guu03*D100 +&
     & guu13*D101 + guu23*D102 + guu33*D103)*Du123 + (guu00*D200 +&
     & guu01*D201 + guu02*D202 + guu03*D203)*Du202 + (guu01*D200 +&
     & guu11*D201 + guu12*D202 + guu13*D203)*Du212 + (guu02*D200 +&
     & guu12*D201 + guu22*D202 + guu23*D203)*Du222 + (guu03*D200 +&
     & guu13*D201 + guu23*D202 + guu33*D203)*Du223 + (guu00*D300 +&
     & guu01*D301 + guu02*D302 + guu03*D303)*Du302 + (guu01*D300 +&
     & guu11*D301 + guu12*D302 + guu13*D303)*Du312 + (guu02*D300 +&
     & guu12*D301 + guu22*D302 + guu23*D303)*Du322 + (guu03*D300 +&
     & guu13*D301 + guu23*D302 + guu33*D303)*Du323 + sigma0*(g02*tZ -&
     & t0*Z2)))+0.d0
    
    K_dot(IS4_03) = sigma2*(b1*D103 + b2*D203 + b3*D303) +&
     & alp*(-((4*G002*G302 + 2*G000*G322)* guu02**2) - G300*(4*guu00*&
     & (G001*guu01 + G002*guu02) + 2*(G000*guu00**2 + G033*guu03**2)) -&
     & (4*G012*G312 + 2*G011*G322)* guu12**2 - (4*G013*G313 +&
     & 2*G011*G333)*guu13**2 - (4*G023*G323 + 2*G022*G333)* guu23**2 -&
     & 4* (G301*(G001*guu01**2 + guu01*(G000*guu00 + G002*guu02 +&
     & G003*guu03) + guu03*(G013*guu11 + G023*guu12)) +&
     & G001*((G311*guu01 + G312*guu02 + G313*guu03)* guu11 +&
     & G302*guu00*guu12 + G313*guu01*guu13) + G002*((G312*guu01 +&
     & G322*guu02)*guu22 + G313*guu01*guu23) + (G013*(G313*guu11 +&
     & G323*guu12) + G333* (G013*guu13 + G023*guu23))* guu33 + guu01*&
     & ((G001*G302 + G000*G312)* guu02 + (G013*G300 + G001*G303 +&
     & G000*G313)*guu03 + (G011*G302 + G002*G311)* guu12 + (G011*G303 +&
     & G003*G311)* guu13 + G301* (G011*guu11 + G012*guu12 + G013*guu13)&
     & + (G012*G303 + G003*G312)* guu23 + G302* (G012*guu22 +&
     & G013*guu23) + G013*G303*guu33) + guu02*((G002*G303 +&
     & G000*G323)*guu03 + G300* (G012*guu01 + G023*guu03) + (G002*G312 +&
     & G001*G322)* guu12 + (G012*G303 + G003*G312)* guu13 + G301*&
     & (G012*guu11 + G022*guu12 + G023*guu13) + (G022*G303 + G003*G322)*&
     & guu23 + G302* (G000*guu00 + G003*guu03 + G022*guu22 + G023*guu23)&
     & + G023*G303*guu33) + guu03*((G002*G313 + G001*G323)*guu12 +&
     & (G033*G301 + G003*G313 + G001*G333)*guu13 + (G003*G323 +&
     & G002*G333)* guu23 + G302* (G023*guu22 + G033*guu23) +&
     & G033*G303*guu33) + guu12*(G302* (G012*guu02 + G013*guu03) +&
     & (G012*G313 + G011*G323)* guu13 + G311* (G012*guu11 + G023*guu13)&
     & + (G013*G322 + G012*G323)* guu23 + G312* (G001*guu01 + G011*guu11&
     & + G022*guu22 + G023*guu23) + G023*G313*guu33) +&
     & G003*(G303*guu03**2 + guu00*(G301*guu13 + G302*guu23) +&
     & (G313*guu01 + G323*guu02 + G333*guu03)*guu33) + guu00*((G003*G300&
     & + G000*G303)*guu03 + G301* (G001*guu11 + G002*guu12) +&
     & G002*G302*guu22 + G303*(G001*guu13 + G002*guu23 + G003*guu33)) +&
     & guu22*(G012*(G312*guu11 + G322*guu12) + G023*G322*guu23 +&
     & G323*(G002*guu03 + G012*guu13 + G023*guu33)) +&
     & guu13*(G013*(G303*guu03 + G311*guu11 + G312*guu12) +&
     & G012*G333*guu23 + G323* (G001*guu02 + G013*guu23) + G312*&
     & (G023*guu22 + G033*guu23) + G313*(G011*guu11 + G033*guu33)) +&
     & guu23*((G013*G312 + G012*G313)*guu11 + G022*G313*guu12 +&
     & G023*(G303*guu03 + G313*guu13) + G323*(G002*guu02 + G022*guu22 +&
     & G033*guu33))) + pi*(-16*T03 + 8*g03*trT) - (K03*ttK)/2. -&
     & tK1*Du103 - tK2*Du203 - tK3*Du303 - 2*((G011*G300 + G000*G311)*&
     & guu01**2 + G311*(G011*guu11**2 + G033*guu13**2) +&
     & G022*(G300*guu02**2 + G311*guu12**2 + G322*guu22**2) +&
     & G033*G322*guu23**2 + G333*(G000*guu03**2 + G033*guu33**2) +&
     & (guu01*K00 + guu11*K01 + guu12*K02)* K13 + (guu02*K00 + guu12*K01&
     & + guu22*K02)*K23 + (guu13*K01 + guu23*K02)*K33 + K03*(guu00*K00 &
     & + guu01*K01 + guu02*K02 + guu13*K13 + guu23*K23 + guu33*K33) &
     & + guu03*(K03**2 + K00*K33) + sigma0*t3*Z0) + 2*(G003*Hu0 +&
     & G103*Hu1 + G203*Hu2 + G303*Hu3 + (guu00*D100 + guu01*D101 +&
     & guu02*D102 + guu03*D103)*Du103 + (guu01*D100 + guu11*D101 +&
     & guu12*D102 + guu13*D103)*Du113 + (guu02*D100 + guu12*D101 +&
     & guu22*D102 + guu23*D103)*Du123 + (guu03*D100 + guu13*D101 +&
     & guu23*D102 + guu33*D103)*Du133 + (guu00*D200 + guu01*D201 +&
     & guu02*D202 + guu03*D203)*Du203 + (guu01*D200 + guu11*D201 +&
     & guu12*D202 + guu13*D203)*Du213 + (guu02*D200 + guu12*D201 +&
     & guu22*D202 + guu23*D203)*Du223 + (guu03*D200 + guu13*D201 +&
     & guu23*D202 + guu33*D203)*Du233 + (guu00*D300 + guu01*D301 +&
     & guu02*D302 + guu03*D303)*Du303 + (guu01*D300 + guu11*D301 +&
     & guu12*D302 + guu13*D303)*Du313 + (guu02*D300 + guu12*D301 +&
     & guu22*D302 + guu23*D303)*Du323 + (guu03*D300 + guu13*D301 +&
     & guu23*D302 + guu33*D303)*Du333 + sigma0*(g03*tZ - t0*Z3)))+0.d0
    
    K_dot(IS4_11) = sigma2*(b1*D111 + b2*D211 + b3*D311) +&
     & alp*(-(G100*(8*G123*guu02*guu03 + 4*(G122*guu02**2 +&
     & G133*guu03**2))) - G133*(4*G122*guu23**2 + 8*G103*guu03*guu33) -&
     & 8*(G100*guu00* (G101*guu01 + G102*guu02 + G103*guu03) +&
     & G101*((G102*guu00 + G123*guu03)*guu12 + (G103*guu00 +&
     & G133*guu03)* guu13 + guu02*(G112*guu11 + G122*guu12 +&
     & G123*guu13)) + (G133*(G102*guu03 + G112*guu13) + G123*&
     & (G112*guu12 + G113*guu13))* guu23 + guu03* (G102*(G103*guu02 +&
     & G113*guu12) + G113* (G100*guu01 + G101*guu11 + G103*guu13) +&
     & G103*G123*guu23) + G102*(G122*guu02*guu22 + G103*guu00*guu23) +&
     & G112*((G101*guu01 + G102*guu02 + G111*guu11)* guu12 +&
     & guu13*(G103*guu02 + G123*guu22) + G113*guu11*guu23) +&
     & guu02*((G103*G122 + G102*G123)*guu23 + G103*G123*guu33) +&
     & guu12*(G122* (G112*guu22 + G113*guu23) + G123*(G111*guu13 +&
     & G113*guu33)) + guu01*((G101*G102 + G100*G112)*guu02 + G101*&
     & (G103*guu03 + G111*guu11 + G113*guu13) + G102* (G111*guu12 +&
     & G112*guu22) + G103*(G112*guu23 + G113*guu33)) +&
     & G113*(G102*guu01*guu23 + guu13*(G111*guu11 + G112*guu12 +&
     & G133*guu33))) + G111*(-8*G103*guu01*guu13 - 4*(G122*guu12**2 +&
     & G133*guu13**2) + 2*Hu1) - 2*(G100**2*guu00**2 + G111**2*guu11**2&
     & + G122**2*guu22**2 + G133**2*guu33**2 + guu00*K01**2 +&
     & guu11*K11**2) - guu22*(8*G123* (G102*guu03 + G122*guu23) +&
     & 4*(G112**2*guu11 + G123**2*guu33) + 2*K12**2) -&
     & guu33*(4*G113**2*guu11 + 8*G123*G133*guu23 + 2*K13**2) +&
     & pi*(-16*T11 + 8*g11*trT) - (K11*ttK)/2. - tK1*Du111 - tK2*Du211 -&
     & tK3*Du311 + 2*(G011*Hu0 + G211*Hu2 + G311*Hu3 + g11*sigma0*tZ +&
     & (guu00*D101 + guu01*D111 + guu02*D112 + guu03*D113)*Du101 +&
     & (guu01*D101 + guu11*D111 + guu12*D112 + guu13*D113)*Du111 +&
     & (guu02*D101 + guu12*D111 + guu22*D112 + guu23*D113)*Du112 +&
     & (guu03*D101 + guu13*D111 + guu23*D112 + guu33*D113)*Du113 +&
     & (guu00*D201 + guu01*D211 + guu02*D212 + guu03*D213)*Du201 +&
     & (guu01*D201 + guu11*D211 + guu12*D212 + guu13*D213)*Du211 +&
     & (guu02*D201 + guu12*D211 + guu22*D212 + guu23*D213)*Du212 +&
     & (guu03*D201 + guu13*D211 + guu23*D212 + guu33*D213)*Du213 +&
     & (guu00*D301 + guu01*D311 + guu02*D312 + guu03*D313)*Du301 +&
     & (guu01*D301 + guu11*D311 + guu12*D312 + guu13*D313)*Du311 +&
     & (guu02*D301 + guu12*D311 + guu22*D312 + guu23*D313)*Du312 +&
     & (guu03*D301 + guu13*D311 + guu23*D312 + guu33*D313)*Du313) -&
     & 4*((G101**2 + G100*G111)* guu01**2 + G101**2*guu00*guu11 +&
     & G112**2*guu12**2 + G113**2*guu13**2 + G102**2*(guu02**2 +&
     & guu00*guu22) + G123**2*guu23**2 + G103**2*(guu03**2 +&
     & guu00*guu33) + guu23*K12*K13 + K01*(guu01*K11 + guu02*K12 +&
     & guu03*K13) + K11*(guu12*K12 + guu13*K13) + sigma0*t1*Z1))+0.d0
    
    K_dot(IS4_12) = sigma2*(b1*D112 + b2*D212 + b3*D312) +&
     & alp*(-((4*G102*G202 + 2*G100*G222)* guu02**2) - G200*(4*guu00*&
     & (G101*guu01 + G102*guu02) + 2*(G100*guu00**2 + G133*guu03**2)) -&
     & (4*G112*G212 + 2*G111*G222)* guu12**2 - (4*G113*G213 +&
     & 2*G111*G233)*guu13**2 - (4*G123*G223 + 2*G122*G233)* guu23**2 -&
     & 4* (G201*(G101*guu01**2 + guu01*(G100*guu00 + G102*guu02 +&
     & G103*guu03) + guu03*(G113*guu11 + G123*guu12)) +&
     & G101*((G211*guu01 + G212*guu02 + G213*guu03)* guu11 +&
     & G202*guu00*guu12 + G213*guu01*guu13) + G102*((G212*guu01 +&
     & G222*guu02)*guu22 + G213*guu01*guu23) + (G113*(G213*guu11 +&
     & G223*guu12) + G233* (G113*guu13 + G123*guu23))* guu33 + guu01*&
     & ((G101*G202 + G100*G212)* guu02 + (G113*G200 + G101*G203 +&
     & G100*G213)*guu03 + (G111*G202 + G102*G211)* guu12 + (G111*G203 +&
     & G103*G211)* guu13 + G201* (G111*guu11 + G112*guu12 + G113*guu13)&
     & + (G112*G203 + G103*G212)* guu23 + G202* (G112*guu22 +&
     & G113*guu23) + G113*G203*guu33) + guu02*((G102*G203 +&
     & G100*G223)*guu03 + G200* (G112*guu01 + G123*guu03) + (G102*G212 +&
     & G101*G222)* guu12 + (G112*G203 + G103*G212)* guu13 + G201*&
     & (G112*guu11 + G122*guu12 + G123*guu13) + (G122*G203 + G103*G222)*&
     & guu23 + G202* (G100*guu00 + G103*guu03 + G122*guu22 + G123*guu23)&
     & + G123*G203*guu33) + guu03*((G102*G213 + G101*G223)*guu12 +&
     & (G133*G201 + G103*G213 + G101*G233)*guu13 + (G103*G223 +&
     & G102*G233)* guu23 + G202* (G123*guu22 + G133*guu23) +&
     & G133*G203*guu33) + guu12*(G202* (G112*guu02 + G113*guu03) +&
     & (G112*G213 + G111*G223)* guu13 + G211* (G112*guu11 + G123*guu13)&
     & + (G113*G222 + G112*G223)* guu23 + G212* (G101*guu01 + G111*guu11&
     & + G122*guu22 + G123*guu23) + G123*G213*guu33) +&
     & G103*(G203*guu03**2 + guu00*(G201*guu13 + G202*guu23) +&
     & (G213*guu01 + G223*guu02 + G233*guu03)*guu33) + guu00*((G103*G200&
     & + G100*G203)*guu03 + G201* (G101*guu11 + G102*guu12) +&
     & G102*G202*guu22 + G203*(G101*guu13 + G102*guu23 + G103*guu33)) +&
     & guu22*(G112*(G212*guu11 + G222*guu12) + G123*G222*guu23 +&
     & G223*(G102*guu03 + G112*guu13 + G123*guu33)) +&
     & guu13*(G113*(G203*guu03 + G211*guu11 + G212*guu12) +&
     & G112*G233*guu23 + G223* (G101*guu02 + G113*guu23) + G212*&
     & (G123*guu22 + G133*guu23) + G213*(G111*guu11 + G133*guu33)) +&
     & guu23*((G113*G212 + G112*G213)*guu11 + G122*G213*guu12 +&
     & G123*(G203*guu03 + G213*guu13) + G223*(G102*guu02 + G122*guu22 +&
     & G133*guu33))) + pi*(-16*T12 + 8*g12*trT) - (K12*ttK)/2. -&
     & tK1*Du112 - tK2*Du212 - tK3*Du312 - 2*((G111*G200 + G100*G211)*&
     & guu01**2 + G211*(G111*guu11**2 + G133*guu13**2) +&
     & G122*(G200*guu02**2 + G211*guu12**2 + G222*guu22**2) +&
     & G133*G222*guu23**2 + G233*(G100*guu03**2 + G133*guu33**2) +&
     & K02*(guu00*K01 + guu01*K11 + guu02*K12 + guu03*K13) + (guu02*K01&
     & + guu23*K13)*K22 + guu12*(K12**2 + K11*K22) + (guu03*K01 +&
     & guu13*K11 + guu33*K13)* K23 + K12*(guu01*K01 + guu11*K11 +&
     & guu13*K13 + guu22*K22 + guu23*K23) + sigma0*t2*Z1) + 2*&
     & (G012*Hu0 + G112*Hu1 + G212*Hu2 + G312*Hu3 + (guu00*D101 +&
     & guu01*D111 + guu02*D112 + guu03*D113)*Du102 + (guu01*D101 +&
     & guu11*D111 + guu12*D112 + guu13*D113)*Du112 + (guu02*D101 +&
     & guu12*D111 + guu22*D112 + guu23*D113)*Du122 + (guu03*D101 +&
     & guu13*D111 + guu23*D112 + guu33*D113)*Du123 + (guu00*D201 +&
     & guu01*D211 + guu02*D212 + guu03*D213)*Du202 + (guu01*D201 +&
     & guu11*D211 + guu12*D212 + guu13*D213)*Du212 + (guu02*D201 +&
     & guu12*D211 + guu22*D212 + guu23*D213)*Du222 + (guu03*D201 +&
     & guu13*D211 + guu23*D212 + guu33*D213)*Du223 + (guu00*D301 +&
     & guu01*D311 + guu02*D312 + guu03*D313)*Du302 + (guu01*D301 +&
     & guu11*D311 + guu12*D312 + guu13*D313)*Du312 + (guu02*D301 +&
     & guu12*D311 + guu22*D312 + guu23*D313)*Du322 + (guu03*D301 +&
     & guu13*D311 + guu23*D312 + guu33*D313)*Du323 + sigma0*(g12*tZ -&
     & t1*Z2)))+0.d0
    
    K_dot(IS4_13) = sigma2*(b1*D113 + b2*D213 + b3*D313) +&
     & alp*(-((4*G102*G302 + 2*G100*G322)* guu02**2) - G300*(4*guu00*&
     & (G101*guu01 + G102*guu02) + 2*(G100*guu00**2 + G133*guu03**2)) -&
     & (4*G112*G312 + 2*G111*G322)* guu12**2 - (4*G113*G313 +&
     & 2*G111*G333)*guu13**2 - (4*G123*G323 + 2*G122*G333)* guu23**2 -&
     & 4* (G301*(G101*guu01**2 + guu01*(G100*guu00 + G102*guu02 +&
     & G103*guu03) + guu03*(G113*guu11 + G123*guu12)) +&
     & G101*((G311*guu01 + G312*guu02 + G313*guu03)* guu11 +&
     & G302*guu00*guu12 + G313*guu01*guu13) + G102*((G312*guu01 +&
     & G322*guu02)*guu22 + G313*guu01*guu23) + (G113*(G313*guu11 +&
     & G323*guu12) + G333* (G113*guu13 + G123*guu23))* guu33 + guu01*&
     & ((G101*G302 + G100*G312)* guu02 + (G113*G300 + G101*G303 +&
     & G100*G313)*guu03 + (G111*G302 + G102*G311)* guu12 + (G111*G303 +&
     & G103*G311)* guu13 + G301* (G111*guu11 + G112*guu12 + G113*guu13)&
     & + (G112*G303 + G103*G312)* guu23 + G302* (G112*guu22 +&
     & G113*guu23) + G113*G303*guu33) + guu02*((G102*G303 +&
     & G100*G323)*guu03 + G300* (G112*guu01 + G123*guu03) + (G102*G312 +&
     & G101*G322)* guu12 + (G112*G303 + G103*G312)* guu13 + G301*&
     & (G112*guu11 + G122*guu12 + G123*guu13) + (G122*G303 + G103*G322)*&
     & guu23 + G302* (G100*guu00 + G103*guu03 + G122*guu22 + G123*guu23)&
     & + G123*G303*guu33) + guu03*((G102*G313 + G101*G323)*guu12 +&
     & (G133*G301 + G103*G313 + G101*G333)*guu13 + (G103*G323 +&
     & G102*G333)* guu23 + G302* (G123*guu22 + G133*guu23) +&
     & G133*G303*guu33) + guu12*(G302* (G112*guu02 + G113*guu03) +&
     & (G112*G313 + G111*G323)* guu13 + G311* (G112*guu11 + G123*guu13)&
     & + (G113*G322 + G112*G323)* guu23 + G312* (G101*guu01 + G111*guu11&
     & + G122*guu22 + G123*guu23) + G123*G313*guu33) +&
     & G103*(G303*guu03**2 + guu00*(G301*guu13 + G302*guu23) +&
     & (G313*guu01 + G323*guu02 + G333*guu03)*guu33) + guu00*((G103*G300&
     & + G100*G303)*guu03 + G301* (G101*guu11 + G102*guu12) +&
     & G102*G302*guu22 + G303*(G101*guu13 + G102*guu23 + G103*guu33)) +&
     & guu22*(G112*(G312*guu11 + G322*guu12) + G123*G322*guu23 +&
     & G323*(G102*guu03 + G112*guu13 + G123*guu33)) +&
     & guu13*(G113*(G303*guu03 + G311*guu11 + G312*guu12) +&
     & G112*G333*guu23 + G323* (G101*guu02 + G113*guu23) + G312*&
     & (G123*guu22 + G133*guu23) + G313*(G111*guu11 + G133*guu33)) +&
     & guu23*((G113*G312 + G112*G313)*guu11 + G122*G313*guu12 +&
     & G123*(G303*guu03 + G313*guu13) + G323*(G102*guu02 + G122*guu22 +&
     & G133*guu33))) + pi*(-16*T13 + 8*g13*trT) - (K13*ttK)/2. -&
     & tK1*Du113 - tK2*Du213 - tK3*Du313 - 2*((G111*G300 + G100*G311)*&
     & guu01**2 + G311*(G111*guu11**2 + G133*guu13**2) +&
     & G122*(G300*guu02**2 + G311*guu12**2 + G322*guu22**2) +&
     & G133*G322*guu23**2 + G333*(G100*guu03**2 + G133*guu33**2) +&
     & K03*(guu00*K01 + guu01*K11 + guu02*K12 + guu03*K13) + (guu02*K01&
     & + guu12*K11 + guu22*K12)* K23 + (guu03*K01 + guu23*K12)*K33 +&
     & K13*(guu01*K01 + guu11*K11 + guu12*K12 + guu23*K23 + guu33*K33)&
     & + guu13*(K13**2 + K11*K33) + sigma0*t3*Z1) + 2*(G013*Hu0 +&
     & G113*Hu1 + G213*Hu2 + G313*Hu3 + (guu00*D101 + guu01*D111 +&
     & guu02*D112 + guu03*D113)*Du103 + (guu01*D101 + guu11*D111 +&
     & guu12*D112 + guu13*D113)*Du113 + (guu02*D101 + guu12*D111 +&
     & guu22*D112 + guu23*D113)*Du123 + (guu03*D101 + guu13*D111 +&
     & guu23*D112 + guu33*D113)*Du133 + (guu00*D201 + guu01*D211 +&
     & guu02*D212 + guu03*D213)*Du203 + (guu01*D201 + guu11*D211 +&
     & guu12*D212 + guu13*D213)*Du213 + (guu02*D201 + guu12*D211 +&
     & guu22*D212 + guu23*D213)*Du223 + (guu03*D201 + guu13*D211 +&
     & guu23*D212 + guu33*D213)*Du233 + (guu00*D301 + guu01*D311 +&
     & guu02*D312 + guu03*D313)*Du303 + (guu01*D301 + guu11*D311 +&
     & guu12*D312 + guu13*D313)*Du313 + (guu02*D301 + guu12*D311 +&
     & guu22*D312 + guu23*D313)*Du323 + (guu03*D301 + guu13*D311 +&
     & guu23*D312 + guu33*D313)*Du333 + sigma0*(g13*tZ - t1*Z3)))+0.d0
    
    K_dot(IS4_22) = sigma2*(b1*D122 + b2*D222 + b3*D322) +&
     & alp*(-(G200*(8*G223*guu02*guu03 + 4*(G222*guu02**2 +&
     & G233*guu03**2))) - G211*(8*G203*guu01*guu13 + 4*(G222*guu12**2 +&
     & G233*guu13**2)) - G233*(4*G222*guu23**2 + 8*G203*guu03*guu33) -&
     & 8*(G200*guu00* (G201*guu01 + G202*guu02 + G203*guu03) +&
     & G201*((G202*guu00 + G223*guu03)*guu12 + (G203*guu00 +&
     & G233*guu03)* guu13 + guu02*(G212*guu11 + G222*guu12 +&
     & G223*guu13)) + (G233*(G202*guu03 + G212*guu13) + G223*&
     & (G212*guu12 + G213*guu13))* guu23 + guu03* (G202*(G203*guu02 +&
     & G213*guu12) + G213* (G200*guu01 + G201*guu11 + G203*guu13) +&
     & G203*G223*guu23) + G202*(G222*guu02*guu22 + G203*guu00*guu23) +&
     & G212*((G201*guu01 + G202*guu02 + G211*guu11)* guu12 +&
     & guu13*(G203*guu02 + G223*guu22) + G213*guu11*guu23) +&
     & guu02*((G203*G222 + G202*G223)*guu23 + G203*G223*guu33) +&
     & guu12*(G222* (G212*guu22 + G213*guu23) + G223*(G211*guu13 +&
     & G213*guu33)) + guu01*((G201*G202 + G200*G212)*guu02 + G201*&
     & (G203*guu03 + G211*guu11 + G213*guu13) + G202* (G211*guu12 +&
     & G212*guu22) + G203*(G212*guu23 + G213*guu33)) +&
     & G213*(G202*guu01*guu23 + guu13*(G211*guu11 + G212*guu12 +&
     & G233*guu33))) - 2*(G200**2*guu00**2 + G211**2*guu11**2 +&
     & G222**2*guu22**2 + G233**2*guu33**2 + guu00*K02**2 +&
     & guu11*K12**2) - guu22*(8*G223* (G202*guu03 + G222*guu23) +&
     & 4*(G212**2*guu11 + G223**2*guu33) + 2*K22**2) -&
     & guu33*(4*G213**2*guu11 + 8*G223*G233*guu23 + 2*K23**2) +&
     & pi*(-16*T22 + 8*g22*trT) - (K22*ttK)/2. - tK1*Du122 - tK2*Du222 -&
     & tK3*Du322 + 2*(G022*Hu0 + G122*Hu1 + G222*Hu2 + G322*Hu3 +&
     & g22*sigma0*tZ + (guu00*D102 + guu01*D112 + guu02*D122 +&
     & guu03*D123)*Du102 + (guu01*D102 + guu11*D112 + guu12*D122 +&
     & guu13*D123)*Du112 + (guu02*D102 + guu12*D112 + guu22*D122 +&
     & guu23*D123)*Du122 + (guu03*D102 + guu13*D112 + guu23*D122 +&
     & guu33*D123)*Du123 + (guu00*D202 + guu01*D212 + guu02*D222 +&
     & guu03*D223)*Du202 + (guu01*D202 + guu11*D212 + guu12*D222 +&
     & guu13*D223)*Du212 + (guu02*D202 + guu12*D212 + guu22*D222 +&
     & guu23*D223)*Du222 + (guu03*D202 + guu13*D212 + guu23*D222 +&
     & guu33*D223)*Du223 + (guu00*D302 + guu01*D312 + guu02*D322 +&
     & guu03*D323)*Du302 + (guu01*D302 + guu11*D312 + guu12*D322 +&
     & guu13*D323)*Du312 + (guu02*D302 + guu12*D312 + guu22*D322 +&
     & guu23*D323)*Du322 + (guu03*D302 + guu13*D312 + guu23*D322 +&
     & guu33*D323)*Du323) - 4*((G201**2 + G200*G211)* guu01**2 +&
     & G201**2*guu00*guu11 + G212**2*guu12**2 + G213**2*guu13**2 +&
     & G202**2*(guu02**2 + guu00*guu22) + G223**2*guu23**2 +&
     & G203**2*(guu03**2 + guu00*guu33) + guu23*K22*K23 + K02*(guu01*K12&
     & + guu02*K22 + guu03*K23) + K12*(guu12*K22 + guu13*K23) +&
     & sigma0*t2*Z2))+0.d0
    
    K_dot(IS4_23) = sigma2*(b1*D123 + b2*D223 + b3*D323) +&
     & alp*(-((4*G202*G302 + 2*G200*G322)* guu02**2) - G300*(4*guu00*&
     & (G201*guu01 + G202*guu02) + 2*(G200*guu00**2 + G233*guu03**2)) -&
     & (4*G212*G312 + 2*G211*G322)* guu12**2 - (4*G213*G313 +&
     & 2*G211*G333)*guu13**2 - (4*G223*G323 + 2*G222*G333)* guu23**2 -&
     & 4* (G301*(G201*guu01**2 + guu01*(G200*guu00 + G202*guu02 +&
     & G203*guu03) + guu03*(G213*guu11 + G223*guu12)) +&
     & G201*((G311*guu01 + G312*guu02 + G313*guu03)* guu11 +&
     & G302*guu00*guu12 + G313*guu01*guu13) + G202*((G312*guu01 +&
     & G322*guu02)*guu22 + G313*guu01*guu23) + (G213*(G313*guu11 +&
     & G323*guu12) + G333* (G213*guu13 + G223*guu23))* guu33 + guu01*&
     & ((G201*G302 + G200*G312)* guu02 + (G213*G300 + G201*G303 +&
     & G200*G313)*guu03 + (G211*G302 + G202*G311)* guu12 + (G211*G303 +&
     & G203*G311)* guu13 + G301* (G211*guu11 + G212*guu12 + G213*guu13)&
     & + (G212*G303 + G203*G312)* guu23 + G302* (G212*guu22 +&
     & G213*guu23) + G213*G303*guu33) + guu02*((G202*G303 +&
     & G200*G323)*guu03 + G300* (G212*guu01 + G223*guu03) + (G202*G312 +&
     & G201*G322)* guu12 + (G212*G303 + G203*G312)* guu13 + G301*&
     & (G212*guu11 + G222*guu12 + G223*guu13) + (G222*G303 + G203*G322)*&
     & guu23 + G302* (G200*guu00 + G203*guu03 + G222*guu22 + G223*guu23)&
     & + G223*G303*guu33) + guu03*((G202*G313 + G201*G323)*guu12 +&
     & (G233*G301 + G203*G313 + G201*G333)*guu13 + (G203*G323 +&
     & G202*G333)* guu23 + G302* (G223*guu22 + G233*guu23) +&
     & G233*G303*guu33) + guu12*(G302* (G212*guu02 + G213*guu03) +&
     & (G212*G313 + G211*G323)* guu13 + G311* (G212*guu11 + G223*guu13)&
     & + (G213*G322 + G212*G323)* guu23 + G312* (G201*guu01 + G211*guu11&
     & + G222*guu22 + G223*guu23) + G223*G313*guu33) +&
     & G203*(G303*guu03**2 + guu00*(G301*guu13 + G302*guu23) +&
     & (G313*guu01 + G323*guu02 + G333*guu03)*guu33) + guu00*((G203*G300&
     & + G200*G303)*guu03 + G301* (G201*guu11 + G202*guu12) +&
     & G202*G302*guu22 + G303*(G201*guu13 + G202*guu23 + G203*guu33)) +&
     & guu22*(G212*(G312*guu11 + G322*guu12) + G223*G322*guu23 +&
     & G323*(G202*guu03 + G212*guu13 + G223*guu33)) +&
     & guu13*(G213*(G303*guu03 + G311*guu11 + G312*guu12) +&
     & G212*G333*guu23 + G323* (G201*guu02 + G213*guu23) + G312*&
     & (G223*guu22 + G233*guu23) + G313*(G211*guu11 + G233*guu33)) +&
     & guu23*((G213*G312 + G212*G313)*guu11 + G222*G313*guu12 +&
     & G223*(G303*guu03 + G313*guu13) + G323*(G202*guu02 + G222*guu22 +&
     & G233*guu33))) + pi*(-16*T23 + 8*g23*trT) - (K23*ttK)/2. -&
     & tK1*Du123 - tK2*Du223 - tK3*Du323 - 2*((G211*G300 + G200*G311)*&
     & guu01**2 + G311*(G211*guu11**2 + G233*guu13**2) +&
     & G222*(G300*guu02**2 + G311*guu12**2 + G322*guu22**2) +&
     & G233*G322*guu23**2 + G333*(G200*guu03**2 + G233*guu33**2) +&
     & K03*(guu00*K02 + guu01*K12 + guu02*K22 + guu03*K23) +&
     & K13*(guu01*K02 + guu11*K12 + guu12*K22 + guu13*K23) + (guu03*K02&
     & + guu13*K12)*K33 + K23*(guu02*K02 + guu12*K12 + guu22*K22 +&
     & guu33*K33) + guu23*(K23**2 + K22*K33) + sigma0*t3*Z2) +&
     & 2*(G023*Hu0 + G123*Hu1 + G223*Hu2 + G323*Hu3 + (guu00*D102 +&
     & guu01*D112 + guu02*D122 + guu03*D123)*Du103 + (guu01*D102 +&
     & guu11*D112 + guu12*D122 + guu13*D123)*Du113 + (guu02*D102 +&
     & guu12*D112 + guu22*D122 + guu23*D123)*Du123 + (guu03*D102 +&
     & guu13*D112 + guu23*D122 + guu33*D123)*Du133 + (guu00*D202 +&
     & guu01*D212 + guu02*D222 + guu03*D223)*Du203 + (guu01*D202 +&
     & guu11*D212 + guu12*D222 + guu13*D223)*Du213 + (guu02*D202 +&
     & guu12*D212 + guu22*D222 + guu23*D223)*Du223 + (guu03*D202 +&
     & guu13*D212 + guu23*D222 + guu33*D223)*Du233 + (guu00*D302 +&
     & guu01*D312 + guu02*D322 + guu03*D323)*Du303 + (guu01*D302 +&
     & guu11*D312 + guu12*D322 + guu13*D323)*Du313 + (guu02*D302 +&
     & guu12*D312 + guu22*D322 + guu23*D323)*Du323 + (guu03*D302 +&
     & guu13*D312 + guu23*D322 + guu33*D323)*Du333 + sigma0*(g23*tZ -&
     & t2*Z3)))+0.d0
    
    K_dot(IS4_33) = sigma2*(b1*D133 + b2*D233 + b3*D333) +&
     & alp*(-(G300*(8*G323*guu02*guu03 + 4*(G322*guu02**2 +&
     & G333*guu03**2))) - G311*(8*G303*guu01*guu13 + 4*(G322*guu12**2 +&
     & G333*guu13**2)) - 8*(G300*guu00* (G301*guu01 + G302*guu02 +&
     & G303*guu03) + G301*((G302*guu00 + G323*guu03)*guu12 + (G303*guu00&
     & + G333*guu03)* guu13 + guu02*(G312*guu11 + G322*guu12 +&
     & G323*guu13)) + (G333*(G302*guu03 + G312*guu13) + G323*&
     & (G312*guu12 + G313*guu13))* guu23 + guu03* (G302*(G303*guu02 +&
     & G313*guu12) + G313* (G300*guu01 + G301*guu11 + G303*guu13) +&
     & G303*G323*guu23) + G302*(G322*guu02*guu22 + G303*guu00*guu23) +&
     & G312*((G301*guu01 + G302*guu02 + G311*guu11)* guu12 +&
     & guu13*(G303*guu02 + G323*guu22) + G313*guu11*guu23) +&
     & guu02*((G303*G322 + G302*G323)*guu23 + G303*G323*guu33) +&
     & guu12*(G322* (G312*guu22 + G313*guu23) + G323*(G311*guu13 +&
     & G313*guu33)) + guu01*((G301*G302 + G300*G312)*guu02 + G301*&
     & (G303*guu03 + G311*guu11 + G313*guu13) + G302* (G311*guu12 +&
     & G312*guu22) + G303*(G312*guu23 + G313*guu33)) +&
     & G313*(G302*guu01*guu23 + guu13*(G311*guu11 + G312*guu12 +&
     & G333*guu33))) + G333*(-4*G322*guu23**2 - 8*G303*guu03*guu33 +&
     & 2*Hu3) - 2*(G300**2*guu00**2 + G311**2*guu11**2 +&
     & G322**2*guu22**2 + G333**2*guu33**2 + guu00*K03**2 +&
     & guu11*K13**2) - guu22*(8*G323* (G302*guu03 + G322*guu23) +&
     & 4*(G312**2*guu11 + G323**2*guu33) + 2*K23**2) + pi*(-16*T33 +&
     & 8*g33*trT) - (K33*ttK)/2. - tK1*Du133 - tK2*Du233 - tK3*Du333 +&
     & 2*(G033*Hu0 + G133*Hu1 + G233*Hu2 + g33*sigma0*tZ + (guu00*D103 +&
     & guu01*D113 + guu02*D123 + guu03*D133)*Du103 + (guu01*D103 +&
     & guu11*D113 + guu12*D123 + guu13*D133)*Du113 + (guu02*D103 +&
     & guu12*D113 + guu22*D123 + guu23*D133)*Du123 + (guu03*D103 +&
     & guu13*D113 + guu23*D123)*Du133 + (guu00*D203 + guu01*D213 +&
     & guu02*D223 + guu03*D233)*Du203 + (guu01*D203 + guu11*D213 +&
     & guu12*D223 + guu13*D233)*Du213 + (guu02*D203 + guu12*D213 +&
     & guu22*D223 + guu23*D233)*Du223 + (guu03*D203 + guu13*D213 +&
     & guu23*D223)*Du233 + (guu00*D303 + guu01*D313 + guu02*D323 +&
     & guu03*D333)*Du303 + (guu01*D303 + guu11*D313 + guu12*D323 +&
     & guu13*D333)*Du313 + (guu02*D303 + guu12*D313 + guu22*D323 +&
     & guu23*D333)*Du323 + (guu03*D303 + guu13*D313 + guu23*D323)*Du333)&
     & + guu33*(-4*G313**2*guu11 - 8*G323*G333*guu23 - 2*K33**2 +&
     & 2*(D133*Du133 + D233*Du233 + D333*Du333)) - 4*((G301**2 +&
     & G300*G311)*guu01**2 + G301**2*guu00*guu11 + G312**2*guu12**2 +&
     & G313**2*guu13**2 + G302**2*(guu02**2 + guu00*guu22) +&
     & G323**2*guu23**2 + G303**2*(guu03**2 + guu00*guu33) +&
     & guu23*K23*K33 + K03*(guu01*K13 + guu02*K23 + guu03*K33) +&
     & K13*(guu12*K23 + guu13*K33) + sigma0*t3*Z3))+0.d0
    
    dg_dot(IS4_100) = alp*((K00*ttD1)/2. - sigma2*D100 + tD11*Du100 +&
     & tD12*Du200 + tD13*Du300)+0.d0
    dg_dot(IS4_101) = alp*((K01*ttD1)/2. - sigma2*D101 + tD11*Du101 +&
     & tD12*Du201 + tD13*Du301)+0.d0
    dg_dot(IS4_102) = alp*((K02*ttD1)/2. - sigma2*D102 + tD11*Du102 +&
     & tD12*Du202 + tD13*Du302)+0.d0
    dg_dot(IS4_103) = alp*((K03*ttD1)/2. - sigma2*D103 + tD11*Du103 +&
     & tD12*Du203 + tD13*Du303)+0.d0
    dg_dot(IS4_111) = alp*((K11*ttD1)/2. - sigma2*D111 + tD11*Du111 +&
     & tD12*Du211 + tD13*Du311)+0.d0
    dg_dot(IS4_112) = alp*((K12*ttD1)/2. - sigma2*D112 + tD11*Du112 +&
     & tD12*Du212 + tD13*Du312)+0.d0
    dg_dot(IS4_113) = alp*((K13*ttD1)/2. - sigma2*D113 + tD11*Du113 +&
     & tD12*Du213 + tD13*Du313)+0.d0
    dg_dot(IS4_122) = alp*((K22*ttD1)/2. - sigma2*D122 + tD11*Du122 +&
     & tD12*Du222 + tD13*Du322)+0.d0
    dg_dot(IS4_123) = alp*((K23*ttD1)/2. - sigma2*D123 + tD11*Du123 +&
     & tD12*Du223 + tD13*Du323)+0.d0
    dg_dot(IS4_133) = alp*((K33*ttD1)/2. - sigma2*D133 + tD11*Du133 +&
     & tD12*Du233 + tD13*Du333)+0.d0
    
    dg_dot(IS4_200) = alp*((K00*ttD2)/2. - sigma2*D200 + tD21*Du100 +&
     & tD22*Du200 + tD23*Du300)+0.d0
    dg_dot(IS4_201) = alp*((K01*ttD2)/2. - sigma2*D201 + tD21*Du101 +&
     & tD22*Du201 + tD23*Du301)+0.d0
    dg_dot(IS4_202) = alp*((K02*ttD2)/2. - sigma2*D202 + tD21*Du102 +&
     & tD22*Du202 + tD23*Du302)+0.d0
    dg_dot(IS4_203) = alp*((K03*ttD2)/2. - sigma2*D203 + tD21*Du103 +&
     & tD22*Du203 + tD23*Du303)+0.d0
    dg_dot(IS4_211) = alp*((K11*ttD2)/2. - sigma2*D211 + tD21*Du111 +&
     & tD22*Du211 + tD23*Du311)+0.d0
    dg_dot(IS4_212) = alp*((K12*ttD2)/2. - sigma2*D212 + tD21*Du112 +&
     & tD22*Du212 + tD23*Du312)+0.d0
    dg_dot(IS4_213) = alp*((K13*ttD2)/2. - sigma2*D213 + tD21*Du113 +&
     & tD22*Du213 + tD23*Du313)+0.d0
    dg_dot(IS4_222) = alp*((K22*ttD2)/2. - sigma2*D222 + tD21*Du122 +&
     & tD22*Du222 + tD23*Du322)+0.d0
    dg_dot(IS4_223) = alp*((K23*ttD2)/2. - sigma2*D223 + tD21*Du123 +&
     & tD22*Du223 + tD23*Du323)+0.d0
    dg_dot(IS4_233) = alp*((K33*ttD2)/2. - sigma2*D233 + tD21*Du133 +&
     & tD22*Du233 + tD23*Du333)+0.d0
    
    dg_dot(IS4_300) = alp*((K00*ttD3)/2. - sigma2*D300 + tD31*Du100 +&
     & tD32*Du200 + tD33*Du300)+0.d0
    dg_dot(IS4_301) = alp*((K01*ttD3)/2. - sigma2*D301 + tD31*Du101 +&
     & tD32*Du201 + tD33*Du301)+0.d0
    dg_dot(IS4_302) = alp*((K02*ttD3)/2. - sigma2*D302 + tD31*Du102 +&
     & tD32*Du202 + tD33*Du302)+0.d0
    dg_dot(IS4_303) = alp*((K03*ttD3)/2. - sigma2*D303 + tD31*Du103 +&
     & tD32*Du203 + tD33*Du303)+0.d0
    dg_dot(IS4_311) = alp*((K11*ttD3)/2. - sigma2*D311 + tD31*Du111 +&
     & tD32*Du211 + tD33*Du311)+0.d0
    dg_dot(IS4_312) = alp*((K12*ttD3)/2. - sigma2*D312 + tD31*Du112 +&
     & tD32*Du212 + tD33*Du312)+0.d0
    dg_dot(IS4_313) = alp*((K13*ttD3)/2. - sigma2*D313 + tD31*Du113 +&
     & tD32*Du213 + tD33*Du313)+0.d0
    dg_dot(IS4_322) = alp*((K22*ttD3)/2. - sigma2*D322 + tD31*Du122 +&
     & tD32*Du222 + tD33*Du322)+0.d0
    dg_dot(IS4_323) = alp*((K23*ttD3)/2. - sigma2*D323 + tD31*Du123 +&
     & tD32*Du223 + tD33*Du323)+0.d0
    dg_dot(IS4_333) = alp*((K33*ttD3)/2. - sigma2*D333 + tD31*Du133 +&
     & tD32*Du233 + tD33*Du333)+0.d0
      
    
    end if !--evolve_geometry-------------------------------------------
      
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !----EVOLVE ONLY THE SCALAR FIELD ---------------------!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if (evolve_scalarfield .EQ. 1) then   
   
   !---auxiliary fields ------------------------------
      
    phi2rc = phic**2 + phir**2    
    phi2mn = phin**2 + phim**2
    
    dphim0 = b1*dphim1 + b2*dphim2 + b3*dphim3 - alp*pim
    dphin0 = b1*dphin1 + b2*dphin2 + b3*dphin3 - alp*pin
    dphir0 = b1*dphir1 + b2*dphir2 + b3*dphir3 - alp*pir
    dphic0 = b1*dphic1 + b2*dphic2 + b3*dphic3 - alp*pic
            
    duphim0 = dphim0*guu00 + dphim1*guu01 + dphim2*guu02 + dphim3*guu03
    duphim1 = dphim0*guu01 + dphim1*guu11 + dphim2*guu12 + dphim3*guu13
    duphim2 = dphim0*guu02 + dphim1*guu12 + dphim2*guu22 + dphim3*guu23
    duphim3 = dphim0*guu03 + dphim1*guu13 + dphim2*guu23 + dphim3*guu33

    duphin0 = dphin0*guu00 + dphin1*guu01 + dphin2*guu02 + dphin3*guu03
    duphin1 = dphin0*guu01 + dphin1*guu11 + dphin2*guu12 + dphin3*guu13
    duphin2 = dphin0*guu02 + dphin1*guu12 + dphin2*guu22 + dphin3*guu23
    duphin3 = dphin0*guu03 + dphin1*guu13 + dphin2*guu23 + dphin3*guu33

    duphir0 = dphir0*guu00 + dphir1*guu01 + dphir2*guu02 + dphir3*guu03
    duphir1 = dphir0*guu01 + dphir1*guu11 + dphir2*guu12 + dphir3*guu13
    duphir2 = dphir0*guu02 + dphir1*guu12 + dphir2*guu22 + dphir3*guu23
    duphir3 = dphir0*guu03 + dphir1*guu13 + dphir2*guu23 + dphir3*guu33

    duphic0 = dphic0*guu00 + dphic1*guu01 + dphic2*guu02 + dphic3*guu03
    duphic1 = dphic0*guu01 + dphic1*guu11 + dphic2*guu12 + dphic3*guu13
    duphic2 = dphic0*guu02 + dphic1*guu12 + dphic2*guu22 + dphic3*guu23
    duphic3 = dphic0*guu03 + dphic1*guu13 + dphic2*guu23 + dphic3*guu33
            
   !---the principal part ---------------------------
      
    Fpim = -(b1*dx_pim + b2*dy_pim + b3*dz_pim) &
   &       + alp*(huu11*dx_d1phim + huu12*dx_d2phim + huu13*dx_d3phim &
   &            + huu12*dy_d1phim + huu22*dy_d2phim + huu23*dy_d3phim &
   &            + huu13*dz_d1phim + huu23*dz_d2phim + huu33*dz_d3phim)&
   &       + sigma2*(b1*dx_phim + b2*dy_phim + b3*dz_phim) 

    Fpin = -(b1*dx_pin + b2*dy_pin + b3*dz_pin) &
   &       + alp*(huu11*dx_d1phin + huu12*dx_d2phin + huu13*dx_d3phin &
   &            + huu12*dy_d1phin + huu22*dy_d2phin + huu23*dy_d3phin &
   &            + huu13*dz_d1phin + huu23*dz_d2phin + huu33*dz_d3phin)&
   &       + sigma2*(b1*dx_phin + b2*dy_phin + b3*dz_phin) 
    
    Fpir = -(b1*dx_pir + b2*dy_pir + b3*dz_pir) &
   &       + alp*(huu11*dx_d1phir + huu12*dx_d2phir + huu13*dx_d3phir &
   &            + huu12*dy_d1phir + huu22*dy_d2phir + huu23*dy_d3phir &
   &            + huu13*dz_d1phir + huu23*dz_d2phir + huu33*dz_d3phir)&
   &       + sigma2*(b1*dx_phir + b2*dy_phir + b3*dz_phir) 

    Fpic = -(b1*dx_pic + b2*dy_pic + b3*dz_pic) &
   &       + alp*(huu11*dx_d1phic + huu12*dx_d2phic + huu13*dx_d3phic &
   &            + huu12*dy_d1phic + huu22*dy_d2phic + huu23*dy_d3phic &
   &            + huu13*dz_d1phic + huu23*dz_d2phic + huu33*dz_d3phic)&
   &       + sigma2*(b1*dx_phic + b2*dy_phic + b3*dz_phic) 

   
    Fdphim(I_1) = - (b1*dx_d1phim + b2*dy_d1phim + b3*dz_d1phim) & 
   &              + alp*dx_pim - alp*sigma2*dx_phim
    Fdphim(I_2) = - (b1*dx_d2phim + b2*dy_d2phim + b3*dz_d2phim) & 
   &              + alp*dy_pim - alp*sigma2*dy_phim
    Fdphim(I_3) = - (b1*dx_d3phim + b2*dy_d3phim + b3*dz_d3phim) & 
   &              + alp*dz_pim - alp*sigma2*dz_phim

    Fdphin(I_1) = - (b1*dx_d1phin + b2*dy_d1phin + b3*dz_d1phin) & 
   &              + alp*dx_pin - alp*sigma2*dx_phin
    Fdphin(I_2) = - (b1*dx_d2phin + b2*dy_d2phin + b3*dz_d2phin) & 
   &              + alp*dy_pin - alp*sigma2*dy_phin
    Fdphin(I_3) = - (b1*dx_d3phin + b2*dy_d3phin + b3*dz_d3phin) & 
   &              + alp*dz_pin - alp*sigma2*dz_phin
      
    Fdphir(I_1) = - (b1*dx_d1phir + b2*dy_d1phir + b3*dz_d1phir) & 
   &              + alp*dx_pir - alp*sigma2*dx_phir
    Fdphir(I_2) = - (b1*dx_d2phir + b2*dy_d2phir + b3*dz_d2phir) & 
   &              + alp*dy_pir - alp*sigma2*dy_phir
    Fdphir(I_3) = - (b1*dx_d3phir + b2*dy_d3phir + b3*dz_d3phir) & 
   &              + alp*dz_pir - alp*sigma2*dz_phir

    Fdphic(I_1) = - (b1*dx_d1phic + b2*dy_d1phic + b3*dz_d1phic) & 
   &              + alp*dx_pic - alp*sigma2*dx_phic
    Fdphic(I_2) = - (b1*dx_d2phic + b2*dy_d2phic + b3*dz_d2phic) & 
   &              + alp*dy_pic - alp*sigma2*dy_phic
    Fdphic(I_3) = - (b1*dx_d3phic + b2*dy_d3phic + b3*dz_d3phic) & 
   &              + alp*dz_pic - alp*sigma2*dz_phic
       
   !---the source terms ---------------------------
     
    phim_dot = b1*dphim1 + b2*dphim2 + b3*dphim3 - alp*pim+0.d0
    pim_dot = (b1*dphim1 + b2*dphim2 + b3*dphim3)* sigma2 +&
     & alp*(phim*(phi2mn*sf_lam + sf_m**2) - (dphim1*huu11 + dphim2*huu12&
     & + dphim3*huu13)*tK1 - (dphim1*huu12 + dphim2*huu22 +&
     & dphim3*huu23)*tK2 - (dphim1*huu13 + dphim2*huu23 +&
     & dphim3*huu33)*tK3 + duphim0*trG0 + duphim1*trG1 + duphim2*trG2 +&
     & duphim3*trG3 - (pim*ttK)/2.)+0.d0    
    dphim_dot(I_1) = alp*((dphim2*huu12 + dphim3*huu13)*tD11 +&
     & (dphim2*huu22 + dphim3*huu23)*tD12 + (dphim2*huu23 +&
     & dphim3*huu33)*tD13 + dphim1*(-sigma2 + huu11*tD11 + huu12*tD12 +&
     & huu13*tD13) + (pim*ttD1)/2.)+0.d0
    dphim_dot(I_2) = alp*((dphim1*huu11 + dphim3*huu13)*tD21 +&
     & (dphim1*huu12 + dphim3*huu23)*tD22 + (dphim1*huu13 +&
     & dphim3*huu33)*tD23 + dphim2*(-sigma2 + huu12*tD21 + huu22*tD22 +&
     & huu23*tD23) + (pim*ttD2)/2.)+0.d0
    dphim_dot(I_3) = alp*((dphim1*huu11 + dphim2*huu12)*tD31 +&
     & (dphim1*huu12 + dphim2*huu22)*tD32 + (dphim1*huu13 +&
     & dphim2*huu23)*tD33 + dphim3*(-sigma2 + huu13*tD31 + huu23*tD32 +&
     & huu33*tD33) + (pim*ttD3)/2.)+0.d0
    
    phin_dot = b1*dphin1 + b2*dphin2 + b3*dphin3 - alp*pin+0.d0
    pin_dot = (b1*dphin1 + b2*dphin2 + b3*dphin3)* sigma2 +&
     & alp*(phin*(phi2mn*sf_lam + sf_m**2) - (dphin1*huu11 + dphin2*huu12&
     & + dphin3*huu13)*tK1 - (dphin1*huu12 + dphin2*huu22 +&
     & dphin3*huu23)*tK2 - (dphin1*huu13 + dphin2*huu23 +&
     & dphin3*huu33)*tK3 + duphin0*trG0 + duphin1*trG1 + duphin2*trG2 +&
     & duphin3*trG3 - (pin*ttK)/2.)+0.d0    
    dphin_dot(I_1) = alp*((dphin2*huu12 + dphin3*huu13)*tD11 +&
     & (dphin2*huu22 + dphin3*huu23)*tD12 + (dphin2*huu23 +&
     & dphin3*huu33)*tD13 + dphin1*(-sigma2 + huu11*tD11 + huu12*tD12 +&
     & huu13*tD13) + (pin*ttD1)/2.)+0.d0
    dphin_dot(I_2) = alp*((dphin1*huu11 + dphin3*huu13)*tD21 +&
     & (dphin1*huu12 + dphin3*huu23)*tD22 + (dphin1*huu13 +&
     & dphin3*huu33)*tD23 + dphin2*(-sigma2 + huu12*tD21 + huu22*tD22 +&
     & huu23*tD23) + (pin*ttD2)/2.)+0.d0
    dphin_dot(I_3) = alp*((dphin1*huu11 + dphin2*huu12)*tD31 +&
     & (dphin1*huu12 + dphin2*huu22)*tD32 + (dphin1*huu13 +&
     & dphin2*huu23)*tD33 + dphin3*(-sigma2 + huu13*tD31 + huu23*tD32 +&
     & huu33*tD33) + (pin*ttD3)/2.)+0.d0
 
    phir_dot = b1*dphir1 + b2*dphir2 + b3*dphir3 - alp*pir+0.d0
    pir_dot = (b1*dphir1 + b2*dphir2 + b3*dphir3)* sigma2 +&
     & alp*(phir*(phi2rc*sf_lam + sf_m**2) - (dphir1*huu11 + dphir2*huu12&
     & + dphir3*huu13)*tK1 - (dphir1*huu12 + dphir2*huu22 +&
     & dphir3*huu23)*tK2 - (dphir1*huu13 + dphir2*huu23 +&
     & dphir3*huu33)*tK3 + duphir0*trG0 + duphir1*trG1 + duphir2*trG2 +&
     & duphir3*trG3 - (pir*ttK)/2.)+0.d0    
    dphir_dot(I_1) = alp*((dphir2*huu12 + dphir3*huu13)*tD11 +&
     & (dphir2*huu22 + dphir3*huu23)*tD12 + (dphir2*huu23 +&
     & dphir3*huu33)*tD13 + dphir1*(-sigma2 + huu11*tD11 + huu12*tD12 +&
     & huu13*tD13) + (pir*ttD1)/2.)+0.d0
    dphir_dot(I_2) = alp*((dphir1*huu11 + dphir3*huu13)*tD21 +&
     & (dphir1*huu12 + dphir3*huu23)*tD22 + (dphir1*huu13 +&
     & dphir3*huu33)*tD23 + dphir2*(-sigma2 + huu12*tD21 + huu22*tD22 +&
     & huu23*tD23) + (pir*ttD2)/2.)+0.d0
    dphir_dot(I_3) = alp*((dphir1*huu11 + dphir2*huu12)*tD31 +&
     & (dphir1*huu12 + dphir2*huu22)*tD32 + (dphir1*huu13 +&
     & dphir2*huu23)*tD33 + dphir3*(-sigma2 + huu13*tD31 + huu23*tD32 +&
     & huu33*tD33) + (pir*ttD3)/2.)+0.d0
    
    phic_dot = b1*dphic1 + b2*dphic2 + b3*dphic3 - alp*pic+0.d0
    pic_dot = (b1*dphic1 + b2*dphic2 + b3*dphic3)* sigma2 +&
     & alp*(phic*(phi2rc*sf_lam + sf_m**2) - (dphic1*huu11 + dphic2*huu12&
     & + dphic3*huu13)*tK1 - (dphic1*huu12 + dphic2*huu22 +&
     & dphic3*huu23)*tK2 - (dphic1*huu13 + dphic2*huu23 +&
     & dphic3*huu33)*tK3 + duphic0*trG0 + duphic1*trG1 + duphic2*trG2 +&
     & duphic3*trG3 - (pic*ttK)/2.)+0.d0    
    dphic_dot(I_1) = alp*((dphic2*huu12 + dphic3*huu13)*tD11 +&
     & (dphic2*huu22 + dphic3*huu23)*tD12 + (dphic2*huu23 +&
     & dphic3*huu33)*tD13 + dphic1*(-sigma2 + huu11*tD11 + huu12*tD12 +&
     & huu13*tD13) + (pic*ttD1)/2.)+0.d0
    dphic_dot(I_2) = alp*((dphic1*huu11 + dphic3*huu13)*tD21 +&
     & (dphic1*huu12 + dphic3*huu23)*tD22 + (dphic1*huu13 +&
     & dphic3*huu33)*tD23 + dphic2*(-sigma2 + huu12*tD21 + huu22*tD22 +&
     & huu23*tD23) + (pic*ttD2)/2.)+0.d0
    dphic_dot(I_3) = alp*((dphic1*huu11 + dphic2*huu12)*tD31 +&
     & (dphic1*huu12 + dphic2*huu22)*tD32 + (dphic1*huu13 +&
     & dphic2*huu23)*tD33 + dphic3*(-sigma2 + huu13*tD31 + huu23*tD32 +&
     & huu33*tD33) + (pic*ttD3)/2.)+0.d0
  
    end if  !---if evolve_scalarfield-----------------------------------------

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!  THE RHS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
   
    rhs(H_H0) = H_dot(I_0)
    rhs(H_H1) = H_dot(I_1)
    rhs(H_H2) = H_dot(I_2)
    rhs(H_H3) = H_dot(I_3)                    

    rhs(H_G0)   = - FG0 + G0_dot
    rhs(H_D1H0) = - FdH0(I_1) + dH0_dot(I_1)
    rhs(H_D2H0) = - FdH0(I_2) + dH0_dot(I_2)
    rhs(H_D3H0) = - FdH0(I_3) + dH0_dot(I_3)       
                
    rhs(H_G00) = g_dot(IS4_00)
    rhs(H_G01) = g_dot(IS4_01)
    rhs(H_G02) = g_dot(IS4_02)
    rhs(H_G03) = g_dot(IS4_03)
    rhs(H_G11) = g_dot(IS4_11)
    rhs(H_G12) = g_dot(IS4_12)
    rhs(H_G13) = g_dot(IS4_13)
    rhs(H_G22) = g_dot(IS4_22)
    rhs(H_G23) = g_dot(IS4_23)
    rhs(H_G33) = g_dot(IS4_33)
    
    rhs(H_K00) = - FK(IS4_00) + K_dot(IS4_00)
    rhs(H_K01) = - FK(IS4_01) + K_dot(IS4_01)
    rhs(H_K02) = - FK(IS4_02) + K_dot(IS4_02)
    rhs(H_K03) = - FK(IS4_03) + K_dot(IS4_03)
    rhs(H_K11) = - FK(IS4_11) + K_dot(IS4_11)
    rhs(H_K12) = - FK(IS4_12) + K_dot(IS4_12)
    rhs(H_K13) = - FK(IS4_13) + K_dot(IS4_13)
    rhs(H_K22) = - FK(IS4_22) + K_dot(IS4_22)
    rhs(H_K23) = - FK(IS4_23) + K_dot(IS4_23)
    rhs(H_K33) = - FK(IS4_33) + K_dot(IS4_33)

    rhs(H_D1G00) = - FD(IS4_100) + dg_dot(IS4_100)
    rhs(H_D1G01) = - FD(IS4_101) + dg_dot(IS4_101)
    rhs(H_D1G02) = - FD(IS4_102) + dg_dot(IS4_102)
    rhs(H_D1G03) = - FD(IS4_103) + dg_dot(IS4_103)
    rhs(H_D1G11) = - FD(IS4_111) + dg_dot(IS4_111)
    rhs(H_D1G12) = - FD(IS4_112) + dg_dot(IS4_112)
    rhs(H_D1G13) = - FD(IS4_113) + dg_dot(IS4_113)
    rhs(H_D1G22) = - FD(IS4_122) + dg_dot(IS4_122)
    rhs(H_D1G23) = - FD(IS4_123) + dg_dot(IS4_123)
    rhs(H_D1G33) = - FD(IS4_133) + dg_dot(IS4_133)
    
    rhs(H_D2G00) = - FD(IS4_200) + dg_dot(IS4_200)
    rhs(H_D2G01) = - FD(IS4_201) + dg_dot(IS4_201)
    rhs(H_D2G02) = - FD(IS4_202) + dg_dot(IS4_202)
    rhs(H_D2G03) = - FD(IS4_203) + dg_dot(IS4_203)
    rhs(H_D2G11) = - FD(IS4_211) + dg_dot(IS4_211)
    rhs(H_D2G12) = - FD(IS4_212) + dg_dot(IS4_212)
    rhs(H_D2G13) = - FD(IS4_213) + dg_dot(IS4_213)
    rhs(H_D2G22) = - FD(IS4_222) + dg_dot(IS4_222)
    rhs(H_D2G23) = - FD(IS4_223) + dg_dot(IS4_223)
    rhs(H_D2G33) = - FD(IS4_233) + dg_dot(IS4_233)
      
    rhs(H_D3G00) = - FD(IS4_300) + dg_dot(IS4_300)
    rhs(H_D3G01) = - FD(IS4_301) + dg_dot(IS4_301)
    rhs(H_D3G02) = - FD(IS4_302) + dg_dot(IS4_302)
    rhs(H_D3G03) = - FD(IS4_303) + dg_dot(IS4_303)
    rhs(H_D3G11) = - FD(IS4_311) + dg_dot(IS4_311)
    rhs(H_D3G12) = - FD(IS4_312) + dg_dot(IS4_312)
    rhs(H_D3G13) = - FD(IS4_313) + dg_dot(IS4_313)
    rhs(H_D3G22) = - FD(IS4_322) + dg_dot(IS4_322)
    rhs(H_D3G23) = - FD(IS4_323) + dg_dot(IS4_323)
    rhs(H_D3G33) = - FD(IS4_333) + dg_dot(IS4_333)              

   if (evolve_scalarfield .EQ. 1) then   
!   the scalar fields
    rhs(H_PHIM)   = phim_dot
    rhs(H_PIM)    = - Fpim + pim_dot
    rhs(H_D1PHIM) = - Fdphim(I_1) + dphim_dot(I_1)
    rhs(H_D2PHIM) = - Fdphim(I_2) + dphim_dot(I_2)
    rhs(H_D3PHIM) = - Fdphim(I_3) + dphim_dot(I_3)    

    rhs(H_PHIN)   = phin_dot
    rhs(H_PIN)    = - Fpin + pin_dot
    rhs(H_D1PHIN) = - Fdphin(I_1) + dphin_dot(I_1)
    rhs(H_D2PHIN) = - Fdphin(I_2) + dphin_dot(I_2)
    rhs(H_D3PHIN) = - Fdphin(I_3) + dphin_dot(I_3)    

    rhs(H_PHIR)   = phir_dot
    rhs(H_PIR)    = - Fpir + pir_dot
    rhs(H_D1PHIR) = - Fdphir(I_1) + dphir_dot(I_1)
    rhs(H_D2PHIR) = - Fdphir(I_2) + dphir_dot(I_2)
    rhs(H_D3PHIR) = - Fdphir(I_3) + dphir_dot(I_3)    

    rhs(H_PHIC)   = phic_dot
    rhs(H_PIC)    = - Fpic + pic_dot
    rhs(H_D1PHIC) = - Fdphic(I_1) + dphic_dot(I_1)
    rhs(H_D2PHIC) = - Fdphic(I_2) + dphic_dot(I_2)
    rhs(H_D3PHIC) = - Fdphic(I_3) + dphic_dot(I_3)    

!    if (r .LT. 1E-2) then
      print*, "***************************************************"
      print*, "in rhs"
      print*, "***************************************************"
      print*, "x,y,z",x,y,z

      print*, "rhs phiR = ", phir_dot
      print*, "rhs phiC = ", phic_dot
      print*, "rhs piR = ", pir_dot
      print*, "rhs piC = ", pic_dot
!    end if

   else
    rhs(H_PHIM)   = 0.d0
    rhs(H_PIM)    = 0.d0
    rhs(H_D1PHIM) = 0.d0
    rhs(H_D2PHIM) = 0.d0
    rhs(H_D3PHIM) = 0.d0
    rhs(H_PHIR)   = 0.d0
    rhs(H_PIR)    = 0.d0
    rhs(H_D1PHIR) = 0.d0
    rhs(H_D2PHIR) = 0.d0
    rhs(H_D3PHIR) = 0.d0
    rhs(H_PHIC)   = 0.d0
    rhs(H_PIC)    = 0.d0
    rhs(H_D1PHIC) = 0.d0
    rhs(H_D2PHIC) = 0.d0
    rhs(H_D3PHIC) = 0.d0
   end if

    do m=1,NU_G
      rhs(m) = rhs(m)*fact
    end do

    end if   !--- if do_evolution ------------------------------------------   

    
           
    if (nancheck) then
       !
       ! Checking for nans, but probably do not
       ! want to bother unless problems arise:
       !
       if (alp2.le.0) then
          write(0,*)myid,']  calc_rhs: alp2<=0!!!!!!!'
       end if
       if (deth.eq.0) then
          write(0,*)myid,']  calc_rhs: Deth=0!!!!!!!'
       end if
       if (isanan(guu00)) then
          write(0,*)myid,']  calc_rhs: guu00 is a NAN!!!'
          nanfound = .true.
       end if
       if (isanan(huu11)) then
          write(0,*)myid,']  calc_rhs: huu11 is a NAN!!!'
          nanfound = .true.
       end if
       do m = 1, NU_G
          if ( isanan(rhs(m)) ) then
             write(0,*)myid,']  calc_rhs: Found Nan: ',m,rhs(m)
             nanfound = .true.
          end if
       end do
       if (nanfound) then
          write(0,*)myid,']  calc_rhs: Nan found on output:'
          do m = 1, NU_G
             write(0,*)myid,']  calc_rhs:      RHS(m):    ',m,RHS(m)
             write(0,*)myid,']  calc_rhs:      u_pt(m):   ',m,u_pt(m)
             write(0,*)myid,']  calc_rhs:      dxu_pt(m): ',m,dxu_pt(m)
             write(0,*)myid,']  calc_rhs:      dyu_pt(m): ',m,dyu_pt(m)
             write(0,*)myid,']  calc_rhs:      dzu_pt(m): ',m,dzu_pt(m)
          end do
          do m = 1, NV_G
             write(0,*)myid,']  calc_rhs:      v_pt(m):   ',m,v_pt(m)
             write(0,*)myid,']  calc_rhs:      dxv_pt(m): ',m,dxv_pt(m)
             write(0,*)myid,']  calc_rhs:      dyv_pt(m): ',m,dyv_pt(m)
             write(0,*)myid,']  calc_rhs:      dzv_pt(m): ',m,dzv_pt(m)
          end do
          do m = 1, NW
             write(0,*)myid,']  calc_rhs:      w_pt(m):   ',m,w_pt(m)
          end do
          call my_exit('rhs: Nan found on output')
       end if
    end if
                                  
#if defined GH_MHD || defined GH_FLOWER
  end subroutine gh_calcrhs
#else
  end subroutine calcrhs
#endif


  !------------------------------------------------------------


  
#ifndef GH_MHD
#ifndef GH_FLOWER
!----------------------------------------------------------------------
!  This routine is empty
!----------------------------------------------------------------------
   !-----------------------------------------------------------------
   !
   !-----------------------------------------------------------------
    subroutine calcrhs2(u_dot, u, v, w, dxu, dyu, dzu, dxv, dyv, dzv, &
         imask, par)

      use GF
      implicit none

      type(gridfunction), dimension(NU_G)   :: dxu, dyu, dzu
      type(gridfunction), dimension(NV_G)   :: dxv, dyv, dzv
      type(gridfunction), dimension(NU_G)   :: u, u_dot
      type(gridfunction), dimension(NW  )   :: w
      type(gridfunction), dimension(NV_G)   :: v
      CCTK_REAL, dimension(:)             :: par
      CCTK_INT, dimension(:,:,:)          :: imask
      
      return
    end subroutine

#endif
#endif

  !#######################################################################!

#if defined GH_MHD || defined GH_FLOWER
   end module gh_rhs_mine
#else
   end module rhs_mine
#endif




