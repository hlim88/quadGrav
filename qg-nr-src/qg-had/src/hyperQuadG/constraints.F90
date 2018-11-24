!-----------------------------------------------------------------------
!
!   $Id: constraints.F90,v 1.7 2010-09-04 19:20:58 carlos Exp $
!
!-----------------------------------------------------------------------
#include "cctk.h"

!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
!      subroutine cal_bssn_analysis(u, v, w,
      subroutine calc_constraints(u, v, w,
#include "CR_con_sub_args.h"
     &           x1d, y1d, z1d, nx, ny, nz, par)

        use GF
        use params
        implicit none

        type(gridfunction), dimension(NU_G)               :: u
        type(gridfunction), dimension(NV_G)               :: v
        type(gridfunction), dimension(NW)                 :: w
        CCTK_REAL, dimension(NPAR),intent(in)             :: par
        CCTK_INT                                          :: nx, ny, nz
        CCTK_REAL, dimension(nx)                          :: x1d
        CCTK_REAL, dimension(ny)                          :: y1d
        CCTK_REAL, dimension(nz)                          :: z1d
#include "CR_con_decl_sub_args.h"

        ! local vars
        CCTK_REAL, dimension(:,:,:), pointer :: Alpha3D
        CCTK_REAL, dimension(:,:,:), pointer :: shiftx
        CCTK_REAL, dimension(:,:,:), pointer :: shifty
        CCTK_REAL, dimension(:,:,:), pointer :: shiftz
        CCTK_REAL, dimension(:,:,:), pointer :: gtxx
        CCTK_REAL, dimension(:,:,:), pointer :: gtxy
        CCTK_REAL, dimension(:,:,:), pointer :: gtxz
        CCTK_REAL, dimension(:,:,:), pointer :: gtyy
        CCTK_REAL, dimension(:,:,:), pointer :: gtyz
        CCTK_REAL, dimension(:,:,:), pointer :: gtzz
        CCTK_REAL, dimension(:,:,:), pointer :: Atxx
        CCTK_REAL, dimension(:,:,:), pointer :: Atxy
        CCTK_REAL, dimension(:,:,:), pointer :: Atxz
        CCTK_REAL, dimension(:,:,:), pointer :: Atyy
        CCTK_REAL, dimension(:,:,:), pointer :: Atyz
        CCTK_REAL, dimension(:,:,:), pointer :: Atzz
        CCTK_REAL, dimension(:,:,:), pointer :: chi3D
        CCTK_REAL, dimension(:,:,:), pointer :: trK3D
        CCTK_REAL, dimension(:,:,:), pointer :: Gamtx
        CCTK_REAL, dimension(:,:,:), pointer :: Gamty
        CCTK_REAL, dimension(:,:,:), pointer :: Gamtz
        CCTK_REAL, dimension(:,:,:), pointer :: gxx
        CCTK_REAL, dimension(:,:,:), pointer :: gxy
        CCTK_REAL, dimension(:,:,:), pointer :: gxz
        CCTK_REAL, dimension(:,:,:), pointer :: gyy
        CCTK_REAL, dimension(:,:,:), pointer :: gyz
        CCTK_REAL, dimension(:,:,:), pointer :: gzz
        CCTK_REAL, dimension(:,:,:), pointer :: div_e
        CCTK_REAL, dimension(:,:,:), pointer :: div_b
        CCTK_REAL, dimension(:,:,:), pointer :: Ex
        CCTK_REAL, dimension(:,:,:), pointer :: Ey
        CCTK_REAL, dimension(:,:,:), pointer :: Ez
        CCTK_REAL, dimension(:,:,:), pointer :: Bx
        CCTK_REAL, dimension(:,:,:), pointer :: By
        CCTK_REAL, dimension(:,:,:), pointer :: Bz
        CCTK_REAL, dimension(:,:,:), pointer :: phiR_sf
        CCTK_REAL, dimension(:,:,:), pointer :: phiI_sf 
        CCTK_REAL, dimension(:,:,:), pointer :: piR_sf
        CCTK_REAL, dimension(:,:,:), pointer :: piI_sf 
        CCTK_REAL, dimension(:,:,:), pointer :: psi4R3d
        CCTK_REAL, dimension(:,:,:), pointer :: psi4I3d
        CCTK_REAL, dimension(:,:,:), pointer :: phi2R3d
        CCTK_REAL, dimension(:,:,:), pointer :: phi2I3d
        CCTK_REAL, dimension(:,:,:), pointer :: phi0R3d
        CCTK_REAL, dimension(:,:,:), pointer :: phi0I3d
        CCTK_REAL, dimension(:,:,:), pointer :: ham
        CCTK_REAL, dimension(:,:,:), pointer :: momx
        CCTK_REAL, dimension(:,:,:), pointer :: momy
        CCTK_REAL, dimension(:,:,:), pointer :: momz
        CCTK_REAL, dimension(:,:,:), pointer :: tr_A
        CCTK_REAL, dimension(:,:,:), pointer :: detgt_m1
        CCTK_REAL, dimension(:,:,:), pointer :: gamtx_con
        CCTK_REAL, dimension(:,:,:), pointer :: gamty_con
        CCTK_REAL, dimension(:,:,:), pointer :: gamtz_con
        CCTK_REAL, dimension(:,:,:), pointer :: calgamtx_con
        CCTK_REAL, dimension(:,:,:), pointer :: calgamty_con
        CCTK_REAL, dimension(:,:,:), pointer :: calgamtz_con

        CCTK_REAL                      :: GNewton
        CCTK_REAL                      :: trK, Alpha, chi, em4chi
        CCTK_REAL                      :: rho_ADM, trPsi1, trS, chi_rhs
!        CCTK_REAL                      :: Alpha_rhs, trK_rhs
        CCTK_REAL                      :: eta, trK0
        CCTK_REAL                      :: inv_chi
        CCTK_REAL,  dimension(3)       :: Betau !, Bu_rhs, Betau_rhs
        CCTK_REAL,  dimension(3)       :: Jtd_ADM, CalGamt, Gamt
!        CCTK_REAL,  dimension(3)       :: d_Alpha, Del_Alpha
        CCTK_REAL,  dimension(3)       :: d_Alpha
!        CCTK_REAL,  dimension(3)       :: d_chi, Del_chi, d_trK
        CCTK_REAL,  dimension(3)       :: d_chi, d_trK
!        CCTK_REAL,  dimension(3)       :: Gamtu_rhs
        CCTK_REAL,  dimension(3)       :: Bu
        CCTK_REAL,  dimension(3,3)     :: gd, gu, gtd, gtu
        CCTK_REAL,  dimension(3,3)     :: EWeyl, BWeyl
        CCTK_REAL,  dimension(3,3)     :: Atd, Atu, Atud, dd_chi, Kd
!        CCTK_REAL,  dimension(3,3)     :: DeltDelt_chi, Rpd, Rtd, R
        CCTK_REAL,  dimension(3,3)     :: Rpd, Rtd, R
!        CCTK_REAL,  dimension(3,3)     :: DelDel_Alpha, DtDt_chi
        CCTK_REAL,  dimension(3,3)     :: Sd, dd_Alpha, d_Gamt
        CCTK_REAL,  dimension(3,3)     :: d_Betau, Psi1, Psi1TF
!        CCTK_REAL,  dimension(3,3)     :: Atd_rhs, gtd_rhs
!        CCTK_REAL,  dimension(3,3)     :: DelDel_chi
        CCTK_REAL,  dimension(3,3,3)   :: dd_Betau, d_gtd, d_Atd
        CCTK_REAL,  dimension(3,3,3)   :: Del_Kd, d_Kd
        CCTK_REAL,  dimension(3,3,3,3) :: dd_gtd
        CCTK_REAL,  dimension(0:3,0:3) :: g4d, g4u, Tu, Tsfu, Temu
        CCTK_REAL,  dimension(0:3,0:3) :: Femd, Femu, Femud
        CCTK_REAL                      :: detgd, idetgd
        CCTK_REAL                      :: detgtd, idetgtd
        CCTK_REAL,  dimension(3,3,3)   :: Ctd, Ct, Chr
        CCTK_REAL                      :: four_Pi_G_exp4chi
        CCTK_REAL                      :: third_trPsi1
        CCTK_REAL                      :: G_exp4chi  ! set?
        CCTK_REAL                      :: expm4chi, exp4chi    ! set?


        CCTK_REAL                      :: HamCon, Rscalar
!        CCTK_REAL,  dimension(3)       :: MomCon, MomCon_v2, MomCon_v3
        CCTK_REAL,  dimension(3)       :: MomCon

        CCTK_REAL                      :: trA, detgtm1
        CCTK_REAL,  dimension(3)       :: gamt_con, calgamt_con
        CCTK_REAL,  dimension(3)       :: cconfun, calcconfun

        CCTK_REAL                      :: dphi4sq, phi2, sfV, dVdphi2
        CCTK_REAL                      :: dphiR4sq, dphiI4sq, dphiIR4sq
        CCTK_REAL                      :: dphiRI4sq, sfVR, sfVI
        CCTK_REAL                      :: dVdphiR2, dVdphiI2
!       maple vars start (need to modify)
        CCTK_REAL                      :: t65, t75, t76, t79, t95, t99
        CCTK_REAL                      :: t40
!       maple vars end
        CCTK_REAL                      :: phiR,phiI,piR,piI
        CCTK_REAL, dimension(3)        :: d_phiR, d_phiI, d_piR, d_piI
        CCTK_REAL,  dimension(0:3)     :: d_phiR4d, d_phiI4d
        CCTK_REAL,  dimension(0:3)     :: d_phiR4u, d_phiI4u
        CCTK_REAL                      :: dil_alpha, dil_mass
        CCTK_REAL                      :: axn_alpha, axn_mass

        CCTK_REAL,  dimension(3)       :: emEu, emBu
        CCTK_REAL,  dimension(3,3)     :: d_emEu, d_emBu
        CCTK_REAL,  dimension(0:3,0:3) :: dualFemu
        CCTK_REAL                      :: divE, divB, Fsq 

        CCTK_REAL                      :: eight_pi_G, sixteen_pi_G
        CCTK_REAL                      :: four_pi_G

        CCTK_REAL                      :: xpsi4R, xpsi4I
        
        CCTK_REAL                      :: temp0, temp1, temp2, temp3, temp4
        CCTK_REAL                      :: temp1Ru, temp2Ru, temp3Ru 
        CCTK_REAL                      :: temp1Iu, temp2Iu, temp3Iu 

        CCTK_REAL, parameter      :: fourth = 0.25d0
        CCTK_REAL, parameter      :: threehalves = 1.5d0
        CCTK_REAL, parameter      :: third = 0.3333333333333333d0
        CCTK_REAL, parameter      :: half = 0.5d0
        CCTK_REAL, parameter      :: twothirds = 0.6666666666666667d0
        CCTK_REAL, parameter      :: one       = 1.0d0
        CCTK_REAL, parameter      :: two       = 2.0d0
        CCTK_REAL, parameter      :: six       = 6.0d0

! temporary variables:
        CCTK_REAL :: v1,v2,v3,v4,v5,v6
        CCTK_INT  :: i1,i2,i3,i4,i5,i6, la, lb, lc
! point-wise variables:
! note - the time index goes to 4, not 0
        CCTK_REAL :: y_ij(3,3),yuu_ij(3,3)
        CCTK_REAL :: g(4,4),g_uu(4,4),bui(3),bdi(3)
! phi2 variables: 
        CCTK_REAL :: Phi2R,Phi2I,Phi0R,Phi0I
! frame variables: 
        CCTK_REAL, dimension(4) :: Tl,Tn,Tm_R,Tm_I
! Faraday tensor:
        CCTK_REAL, dimension(4,4) :: F_uu,F_du,F_dd
!  various EM fields:
        CCTK_REAL :: Bu0,Bu1,Bu2,Bu3,Bd0,Bd1,Bd2,Bd3
        CCTK_REAL :: Eu0,Eu1,Eu2,Eu3,Ed0,Ed1,Ed2,Ed3
        CCTK_REAL :: g00,g01,g02,g03,g10,g11,g12,g13
        CCTK_REAL :: g20,g21,g22,g23,g30,g31,g32,g33
        CCTK_REAL :: tu0,tu1,tu2,tu3
        CCTK_REAL :: alp,deth,eps
        CCTK_REAL :: h11,h12,h13,h22,h23,h33
        CCTK_REAL :: huu11,huu12,huu13,huu22,huu23,huu33
! local values of x,y,z
        CCTK_REAL :: x2,y2,z2, rad, rcyl
        CCTK_REAL :: f_tr, f_rphi,f_tphi,f_thephi,vb1,vb2

        logical, parameter             :: ltrace   = .false.
        CCTK_INT    ::  i, j, k
#include "constraints_temp_vars.h"


        ! G is Newton"s gravitation constant
        GNewton = 1.0d0

        four_pi_G  = 12.5663706143591729d0 * GNewton 
        eight_pi_G = 25.13274123d0 * GNewton
        sixteen_pi_G = 50.26548246d0 * GNewton

        dil_alpha = par(P_DIL_ALPHA)
        dil_mass  = par(P_DIL_MASS)

        axn_alpha = par(P_AXN_ALPHA)
        axn_mass  = par(P_AXN_MASS)

        Alpha3D     => u(G_ALPHA)%d
        shiftx      => u(G_SHIFT1)%d
        shifty      => u(G_SHIFT2)%d
        shiftz      => u(G_SHIFT3)%d
        chi3D       => u(G_CHI)%d
        trK3D       => u(G_TRK)%d
        gtxx        => u(G_GT11)%d
        gtxy        => u(G_GT12)%d
        gtxz        => u(G_GT13)%d
        gtyy        => u(G_GT22)%d
        gtyz        => u(G_GT23)%d
        gtzz        => u(G_GT33)%d
        Atxx        => u(G_A11)%d
        Atxy        => u(G_A12)%d
        Atxz        => u(G_A13)%d
        Atyy        => u(G_A22)%d
        Atyz        => u(G_A23)%d
        Atzz        => u(G_A33)%d
        Gamtx       => u(G_GAM1)%d
        Gamty       => u(G_GAM2)%d
        Gamtz       => u(G_GAM3)%d

        Ex          => u(G_EX)%d
        Ey          => u(G_EY)%d
        Ez          => u(G_EZ)%d
        Bx          => u(G_BX)%d
        By          => u(G_BY)%d
        Bz          => u(G_BZ)%d

        phiR_sf      =>  u(G_PHIR)%d 
        phiI_sf      =>  u(G_PHII)%d 
        piR_sf       =>  u(G_PIR)%d 
        piI_sf       =>  u(G_PII)%d 

        gxx         => v(G_G11)%d
        gxy         => v(G_G12)%d
        gxz         => v(G_G13)%d
        gyy         => v(G_G22)%d
        gyz         => v(G_G23)%d
        gzz         => v(G_G33)%d

!        Tmn00       => v(G_TUU00)%d
!        Tmn01       => v(G_TUU01)%d
!        Tmn02       => v(G_TUU02)%d
!        Tmn03       => v(G_TUU03)%d
!        Tmn11       => v(G_TUU11)%d
!        Tmn12       => v(G_TUU12)%d
!        Tmn13       => v(G_TUU13)%d
!        Tmn22       => v(G_TUU22)%d
!        Tmn23       => v(G_TUU23)%d
!        Tmn33       => v(G_TUU33)%d

        ham         => w(H_HAM)%d
        momx        => w(H_MOMX)%d
        momy        => w(H_MOMY)%d
        momz        => w(H_MOMZ)%d

        div_e        => w(H_DIV_E)%d
        div_b        => w(H_DIV_B)%d

        psi4R3d      => v(G_PSI4R)%d
        psi4I3d      => v(G_PSI4I)%d
        phi0R3d      => v(G_PHI0R)%d
        phi0I3d      => v(G_PHI0I)%d
        phi2R3d      => v(G_PHI2R)%d
        phi2I3d      => v(G_PHI2I)%d

        tr_A          => w(H_TR_A)%d
        detgt_m1      => w(H_DETGT_M1)%d

        gamtx_con     => w(H_GAMTX_CON)%d
        gamty_con     => w(H_GAMTY_CON)%d
        gamtz_con     => w(H_GAMTZ_CON)%d

        calgamtx_con  => w(H_CALGAMTX_CON)%d
        calgamty_con  => w(H_CALGAMTY_CON)%d
        calgamtz_con  => w(H_CALGAMTZ_CON)%d


        do k = 3, nz-2
        do j = 3, ny-2
        do i = 3, nx-2

          gd(1,1) = gxx(i,j,k)
          gd(1,2) = gxy(i,j,k)
          gd(1,3) = gxz(i,j,k)
          gd(2,1) = gxy(i,j,k)
          gd(2,2) = gyy(i,j,k)
          gd(2,3) = gyz(i,j,k)
          gd(3,1) = gxz(i,j,k)
          gd(3,2) = gyz(i,j,k)
          gd(3,3) = gzz(i,j,k)

          gtd(1,1) = gtxx(i,j,k)
          gtd(1,2) = gtxy(i,j,k)
          gtd(1,3) = gtxz(i,j,k)
          gtd(2,1) = gtxy(i,j,k)
          gtd(2,2) = gtyy(i,j,k)
          gtd(2,3) = gtyz(i,j,k)
          gtd(3,1) = gtxz(i,j,k)
          gtd(3,2) = gtyz(i,j,k)
          gtd(3,3) = gtzz(i,j,k)

          Atd(1,1) = Atxx(i,j,k)
          Atd(1,2) = Atxy(i,j,k)
          Atd(1,3) = Atxz(i,j,k)
          Atd(2,1) = Atxy(i,j,k)
          Atd(2,2) = Atyy(i,j,k)
          Atd(2,3) = Atyz(i,j,k)
          Atd(3,1) = Atxz(i,j,k)
          Atd(3,2) = Atyz(i,j,k)
          Atd(3,3) = Atzz(i,j,k)

          trK      = trK3D(i,j,k)
          chi      = chi3D(i,j,k)

          Gamt(1) = Gamtx(i,j,k)
          Gamt(2) = Gamty(i,j,k)
          Gamt(3) = Gamtz(i,j,k)

          Alpha    = Alpha3D(i,j,k)
          Betau(1) = shiftx(i,j,k)
          Betau(2) = shifty(i,j,k)
          Betau(3) = shiftz(i,j,k)

!          Tu(0,0)  = Tmn00(i,j,k)
!          Tu(0,1)  = Tmn01(i,j,k)
!          Tu(0,2)  = Tmn02(i,j,k)
!          Tu(0,3)  = Tmn03(i,j,k)
!
!          Tu(1,0)  = Tu(0,1)
!          Tu(1,1)  = Tmn11(i,j,k)
!          Tu(1,2)  = Tmn12(i,j,k)
!          Tu(1,3)  = Tmn13(i,j,k)
!
!          Tu(2,0)  = Tu(0,2)
!          Tu(2,1)  = Tu(1,2)
!          Tu(2,2)  = Tmn22(i,j,k)
!          Tu(2,3)  = Tmn23(i,j,k)
!
!          Tu(3,0)  = Tu(0,3)
!          Tu(3,1)  = Tu(1,3)
!          Tu(3,2)  = Tu(2,3)
!          Tu(3,3)  = Tmn33(i,j,k)

          Tu(0,0)  = 0.d0 
          Tu(0,1)  = 0.d0 
          Tu(0,2)  = 0.d0 
          Tu(0,3)  = 0.d0 

          Tu(1,0)  = 0.d0 
          Tu(1,1)  = 0.d0 
          Tu(1,2)  = 0.d0 
          Tu(1,3)  = 0.d0 

          Tu(2,0)  = 0.d0 
          Tu(2,1)  = 0.d0 
          Tu(2,2)  = 0.d0 
          Tu(2,3)  = 0.d0 

          Tu(3,0)  = 0.d0 
          Tu(3,1)  = 0.d0
          Tu(3,2)  = 0.d0 
          Tu(3,3)  = 0.d0 

          d_chi(1)      = dx_chi(i,j,k)
          d_chi(2)      = dy_chi(i,j,k)
          d_chi(3)      = dz_chi(i,j,k)

          dd_chi(1,1)   = dxx_chi(i,j,k)
          dd_chi(1,2)   = dxy_chi(i,j,k)
          dd_chi(1,3)   = dxz_chi(i,j,k)
          dd_chi(2,2)   = dyy_chi(i,j,k)
          dd_chi(2,3)   = dyz_chi(i,j,k)
          dd_chi(3,3)   = dzz_chi(i,j,k)

          d_trK(1)      = dx_trK(i,j,k)
          d_trK(2)      = dy_trK(i,j,k)
          d_trK(3)      = dz_trK(i,j,k)

          d_Gamt(1,1)   = dx_Gamtx(i,j,k)
          d_Gamt(2,1)   = dy_Gamtx(i,j,k)
          d_Gamt(3,1)   = dz_Gamtx(i,j,k)
          d_Gamt(1,2)   = dx_Gamty(i,j,k)
          d_Gamt(2,2)   = dy_Gamty(i,j,k)
          d_Gamt(3,2)   = dz_Gamty(i,j,k)
          d_Gamt(1,3)   = dx_Gamtz(i,j,k)
          d_Gamt(2,3)   = dy_Gamtz(i,j,k)
          d_Gamt(3,3)   = dz_Gamtz(i,j,k)

          d_gtd(1,1,1)      = dx_gtxx(i,j,k)
          d_gtd(2,1,1)      = dy_gtxx(i,j,k)
          d_gtd(3,1,1)      = dz_gtxx(i,j,k)

          d_gtd(1,1,2)      = dx_gtxy(i,j,k)
          d_gtd(2,1,2)      = dy_gtxy(i,j,k)
          d_gtd(3,1,2)      = dz_gtxy(i,j,k)

          d_gtd(1,1,3)      = dx_gtxz(i,j,k)
          d_gtd(2,1,3)      = dy_gtxz(i,j,k)
          d_gtd(3,1,3)      = dz_gtxz(i,j,k)

          d_gtd(1,2,1)      = d_gtd(1,1,2)
          d_gtd(2,2,1)      = d_gtd(2,1,2)
          d_gtd(3,2,1)      = d_gtd(3,1,2)

          d_gtd(1,2,2)      = dx_gtyy(i,j,k)
          d_gtd(2,2,2)      = dy_gtyy(i,j,k)
          d_gtd(3,2,2)      = dz_gtyy(i,j,k)

          d_gtd(1,2,3)      = dx_gtyz(i,j,k)
          d_gtd(2,2,3)      = dy_gtyz(i,j,k)
          d_gtd(3,2,3)      = dz_gtyz(i,j,k)

          d_gtd(1,3,1)      = d_gtd(1,1,3)
          d_gtd(2,3,1)      = d_gtd(2,1,3)
          d_gtd(3,3,1)      = d_gtd(3,1,3)

          d_gtd(1,3,2)      = d_gtd(1,2,3)
          d_gtd(2,3,2)      = d_gtd(2,2,3)
          d_gtd(3,3,2)      = d_gtd(3,2,3)

          d_gtd(1,3,3)      = dx_gtzz(i,j,k)
          d_gtd(2,3,3)      = dy_gtzz(i,j,k)
          d_gtd(3,3,3)      = dz_gtzz(i,j,k)

          d_Atd(1,1,1)      = dx_Atxx(i,j,k)
          d_Atd(2,1,1)      = dy_Atxx(i,j,k)
          d_Atd(3,1,1)      = dz_Atxx(i,j,k)

          d_Atd(1,1,2)      = dx_Atxy(i,j,k)
          d_Atd(2,1,2)      = dy_Atxy(i,j,k)
          d_Atd(3,1,2)      = dz_Atxy(i,j,k)

          d_Atd(1,1,3)      = dx_Atxz(i,j,k)
          d_Atd(2,1,3)      = dy_Atxz(i,j,k)
          d_Atd(3,1,3)      = dz_Atxz(i,j,k)

          d_Atd(1,2,1)      = d_Atd(1,1,2)
          d_Atd(2,2,1)      = d_Atd(2,1,2)
          d_Atd(3,2,1)      = d_Atd(3,1,2)

          d_Atd(1,2,2)      = dx_Atyy(i,j,k)
          d_Atd(2,2,2)      = dy_Atyy(i,j,k)
          d_Atd(3,2,2)      = dz_Atyy(i,j,k)

          d_Atd(1,2,3)      = dx_Atyz(i,j,k)
          d_Atd(2,2,3)      = dy_Atyz(i,j,k)
          d_Atd(3,2,3)      = dz_Atyz(i,j,k)

          d_Atd(1,3,1)      = d_Atd(1,1,3)
          d_Atd(2,3,1)      = d_Atd(2,1,3)
          d_Atd(3,3,1)      = d_Atd(3,1,3)

          d_Atd(1,3,2)      = d_Atd(1,2,3)
          d_Atd(2,3,2)      = d_Atd(2,2,3)
          d_Atd(3,3,2)      = d_Atd(3,2,3)

          d_Atd(1,3,3)      = dx_Atzz(i,j,k)
          d_Atd(2,3,3)      = dy_Atzz(i,j,k)
          d_Atd(3,3,3)      = dz_Atzz(i,j,k)

          dd_gtd(1,1,1,1)      = dxx_gtxx(i,j,k)
          dd_gtd(1,2,1,1)      = dxy_gtxx(i,j,k)
          dd_gtd(1,3,1,1)      = dxz_gtxx(i,j,k)
          dd_gtd(2,2,1,1)      = dyy_gtxx(i,j,k)
          dd_gtd(2,3,1,1)      = dyz_gtxx(i,j,k)
          dd_gtd(3,3,1,1)      = dzz_gtxx(i,j,k)

          dd_gtd(1,1,1,2)      = dxx_gtxy(i,j,k)
          dd_gtd(1,2,1,2)      = dxy_gtxy(i,j,k)
          dd_gtd(1,3,1,2)      = dxz_gtxy(i,j,k)
          dd_gtd(2,2,1,2)      = dyy_gtxy(i,j,k)
          dd_gtd(2,3,1,2)      = dyz_gtxy(i,j,k)
          dd_gtd(3,3,1,2)      = dzz_gtxy(i,j,k)

          dd_gtd(1,1,1,3)      = dxx_gtxz(i,j,k)
          dd_gtd(1,2,1,3)      = dxy_gtxz(i,j,k)
          dd_gtd(1,3,1,3)      = dxz_gtxz(i,j,k)
          dd_gtd(2,2,1,3)      = dyy_gtxz(i,j,k)
          dd_gtd(2,3,1,3)      = dyz_gtxz(i,j,k)
          dd_gtd(3,3,1,3)      = dzz_gtxz(i,j,k)

          dd_gtd(1,1,2,2)      = dxx_gtyy(i,j,k)
          dd_gtd(1,2,2,2)      = dxy_gtyy(i,j,k)
          dd_gtd(1,3,2,2)      = dxz_gtyy(i,j,k)
          dd_gtd(2,2,2,2)      = dyy_gtyy(i,j,k)
          dd_gtd(2,3,2,2)      = dyz_gtyy(i,j,k)
          dd_gtd(3,3,2,2)      = dzz_gtyy(i,j,k)

          dd_gtd(1,1,2,3)      = dxx_gtyz(i,j,k)
          dd_gtd(1,2,2,3)      = dxy_gtyz(i,j,k)
          dd_gtd(1,3,2,3)      = dxz_gtyz(i,j,k)
          dd_gtd(2,2,2,3)      = dyy_gtyz(i,j,k)
          dd_gtd(2,3,2,3)      = dyz_gtyz(i,j,k)
          dd_gtd(3,3,2,3)      = dzz_gtyz(i,j,k)

          dd_gtd(1,1,3,3)      = dxx_gtzz(i,j,k)
          dd_gtd(1,2,3,3)      = dxy_gtzz(i,j,k)
          dd_gtd(1,3,3,3)      = dxz_gtzz(i,j,k)
          dd_gtd(2,2,3,3)      = dyy_gtzz(i,j,k)
          dd_gtd(2,3,3,3)      = dyz_gtzz(i,j,k)
          dd_gtd(3,3,3,3)      = dzz_gtzz(i,j,k)

          emEu(1)              = Ex(i,j,k) 
          emEu(2)              = Ey(i,j,k) 
          emEu(3)              = Ez(i,j,k) 

          emBu(1)              = Bx(i,j,k) 
          emBu(2)              = By(i,j,k) 
          emBu(3)              = Bz(i,j,k) 

          d_emEu(1,1)          = dx_Ex(i,j,k) 
          d_emEu(1,2)          = dx_Ey(i,j,k) 
          d_emEu(1,3)          = dx_Ez(i,j,k) 
          d_emEu(2,1)          = dy_Ex(i,j,k) 
          d_emEu(2,2)          = dy_Ey(i,j,k) 
          d_emEu(2,3)          = dy_Ez(i,j,k) 
          d_emEu(3,1)          = dz_Ex(i,j,k) 
          d_emEu(3,2)          = dz_Ey(i,j,k) 
          d_emEu(3,3)          = dz_Ez(i,j,k) 

          d_emBu(1,1)          = dx_Bx(i,j,k) 
          d_emBu(1,2)          = dx_By(i,j,k) 
          d_emBu(1,3)          = dx_Bz(i,j,k) 
          d_emBu(2,1)          = dy_Bx(i,j,k) 
          d_emBu(2,2)          = dy_By(i,j,k) 
          d_emBu(2,3)          = dy_Bz(i,j,k) 
          d_emBu(3,1)          = dz_Bx(i,j,k) 
          d_emBu(3,2)          = dz_By(i,j,k) 
          d_emBu(3,3)          = dz_Bz(i,j,k) 

          phiR               = phiR_Sf(i,j,k)   
          phiI               = phiI_sf(i,j,k)   
          piR                = piR_sf(i,j,k)   
          piI                = piI_sf(i,j,k)   

          d_phiR(1)          = dx_phiR(i,j,k)
          d_phiR(2)          = dy_phiR(i,j,k)
          d_phiR(3)          = dz_phiR(i,j,k)

          d_phiI(1)          = dx_phiI(i,j,k)
          d_phiI(2)          = dy_phiI(i,j,k)
          d_phiI(3)          = dz_phiI(i,j,k)

          d_piR(1)           = dx_piR(i,j,k)
          d_piR(2)           = dy_piR(i,j,k)
          d_piR(3)           = dz_piR(i,j,k)

          d_piI(1)           = dx_piI(i,j,k)
          d_piI(2)           = dy_piI(i,j,k)
          d_piI(3)           = dz_piI(i,j,k)


!------------------------------------------------
         ! compute first the tetrad and the phi2

         ! make the lowered indice 4-metric g4dd
         ! get the lowered physical metric first:

          v1 = 1.d0/chi
          y_ij(1,1) = v1*gd(1,1)
          y_ij(1,2) = v1*gd(1,2)
          y_ij(1,3) = v1*gd(1,3)
          y_ij(2,2) = v1*gd(2,2)
          y_ij(2,3) = v1*gd(2,3)
          y_ij(3,3) = v1*gd(3,3)
          y_ij(2,1) = y_ij(1,2)
          y_ij(3,1) = y_ij(1,3)
          y_ij(3,2) = y_ij(2,3)
          deth = v1**3   

          ! create the up-index spatial metric:
          yuu_ij(1,1)=(-y_ij(2,3)**2 + y_ij(2,2)*y_ij(3,3))/deth
          yuu_ij(1,2)=(y_ij(1,3)*y_ij(2,3) - y_ij(1,2)*y_ij(3,3))/deth
          yuu_ij(1,3)=(-(y_ij(1,3)*y_ij(2,2))+y_ij(1,2)*y_ij(2,3))/deth
          yuu_ij(2,2)=(-y_ij(1,3)**2 + y_ij(1,1)*y_ij(3,3))/deth
          yuu_ij(2,3)=(y_ij(1,2)*y_ij(1,3) - y_ij(1,1)*y_ij(2,3))/deth
          yuu_ij(3,3)=(-y_ij(1,2)**2 + y_ij(1,1)*y_ij(2,2))/deth
          yuu_ij(2,1)=yuu_ij(1,2)
          yuu_ij(3,1)=yuu_ij(1,3)
          yuu_ij(3,2)=yuu_ij(2,3)

! plug in the raised index shift:
          bui(1) = Betau(1)
          bui(2) = Betau(2)
          bui(3) = Betau(3)

! now get the lowered shift vectors:
          do i1 = 1,3
             bdi(i1) = 0.d0
             do i2 = 1,3
                bdi(i1) = bdi(i1) + y_ij(i1,i2)*bui(i2)
             enddo
          enddo

! now create the down-index 4-metric:
          v1 = -(Alpha)**2
          do i1 = 1,3
             v1 = v1 + bdi(i1)*bui(i1)
          enddo
          g(4,4) = v1
          do i1 = 1,3
             g(i1,4) = bdi(i1)
             g(4,i1) = bdi(i1)
          enddo
          do i1 = 1,3
             do i2 = 1,3
                g(i1,i2) = y_ij(i1,i2)
             enddo
          enddo

! and the up-index 4-metric:
          v1 = 1.d0/(Alpha)**2
          g_uu(4,4) = -v1
          do i1 = 1,3
             g_uu(i1,4) = v1*bui(i1)
             g_uu(4,i1) = v1*bui(i1)
          enddo
          do i1 = 1,3
             do i2 = 1,3
                v2 = yuu_ij(i1,i2)
                v2 = v2 - v1*bui(i1)*bui(i2)
                g_uu(i1,i2) = v2
             enddo
          enddo

         ! get the Faraday tensor:
          Bu0 = 0.d0
          Bu1 = emBu(1)
          Bu2 = emBu(2)
          Bu3 = emBu(3)

          Bd0 = g(4,4)*Bu0 + g(4,1)*Bu1 + g(4,2)*Bu2 + g(4,3)*Bu3
          Bd1 = g(4,1)*Bu0 + g(1,1)*Bu1 + g(1,2)*Bu2 + g(1,3)*Bu3
          Bd2 = g(4,2)*Bu0 + g(1,2)*Bu1 + g(2,2)*Bu2 + g(2,3)*Bu3
          Bd3 = g(4,3)*Bu0 + g(1,3)*Bu1 + g(2,3)*Bu2 + g(3,3)*Bu3

          Eu0 = 0.d0
          Eu1 = emEu(1)
          Eu2 = emEu(2)
          Eu3 = emEu(3)
  
          alp = Alpha           
          tu0 = 1.0/alp; tu1 = -bui(1)/alp; 
          tu2 = -bui(2)/alp; tu3 = -bui(3)/alp
          eps = 1.0d0/sqrt(deth)

          x2 = x1d(i)
          y2 = y1d(j)
          z2 = z1d(k)
          rad = sqrt(x2**2 + y2**2 + z2**2)

          !--compute the frame--------- 
          call frame(Tl, Tn, Tm_R, Tm_I, x2, y2, z2, alp, 
     &                     bui, g, g_uu, deth)

          F_uu(1,1) = 0.0d0
          F_uu(1,2) = tu1*Eu2 - tu2*Eu1 + eps*Bd3
          F_uu(1,3) = tu1*Eu3 - tu3*Eu1 - eps*Bd2
          F_uu(1,4) = tu1*Eu0 - tu0*Eu1
          F_uu(2,1) = tu2*Eu1 - tu1*Eu2 - eps*Bd3
          F_uu(2,2) = 0.0d0
          F_uu(2,3) = tu2*Eu3 - tu3*Eu2 + eps*Bd1
          F_uu(2,4) = tu2*Eu0 - tu0*Eu2
          F_uu(3,1) = tu3*Eu1 - tu1*Eu3 + eps*Bd2
          F_uu(3,2) = tu3*Eu2 - tu2*Eu3 - eps*Bd1
          F_uu(3,3) = 0.0d0
          F_uu(3,4) = tu3*Eu0 - tu0*Eu3
          F_uu(4,1) = tu0*Eu1 - tu1*Eu0
          F_uu(4,2) = tu0*Eu2 - tu2*Eu0
          F_uu(4,3) = tu0*Eu3 - tu3*Eu0
          F_uu(4,4) = 0.0d0

!	  print *,'about to enter calc_phi2'	
          call calc_phi2(Phi2R,Phi2I,Phi0R,Phi0I,
     &             g,F_uu,Tl,Tn,Tm_R,Tm_I,x2,y2,z2)

!-------------------------------------------------


#include "constraints.h"


          ham(i,j,k)           = HamCon
          momx(i,j,k)          = MomCon(1)
          momy(i,j,k)          = MomCon(2)
          momz(i,j,k)          = MomCon(3)


!          div_e(i,j,k)         = sixteen_pi_G*rho_ADM

          div_e(i,j,k)         = divE 
          div_b(i,j,k)         = divB 

          psi4R3d(i,j,k)       = rad * xpsi4R
          psi4I3d(i,j,k)       = rad * xpsi4I

          phi0R3d(i,j,k)       = rad * Phi0R
          phi0I3d(i,j,k)       = rad * Phi0I
          phi2R3d(i,j,k)       = rad * Phi2R
          phi2I3d(i,j,k)       = rad * Phi2I

          tr_A(i,j,k)          = trA
          detgt_m1(i,j,k)      = detgtm1

          gamtx_con(i,j,k)     = gamt_con(1)
          gamty_con(i,j,k)     = gamt_con(2)
          gamtz_con(i,j,k)     = gamt_con(3)

          calgamtx_con(i,j,k)  = calgamt_con(1)
          calgamty_con(i,j,k)  = calgamt_con(2)
          calgamtz_con(i,j,k)  = calgamt_con(3)

         if (ltrace) then
         if (i .eq. 31 .and. j .eq. 31 .and. k .eq. 31) then
           print *,'x1d = ',x1d(31)
           print *,'y1d = ',y1d(31)
           print *,'z1d = ',z1d(31)
           print*, "guu00, guu01=", g4u(0,0), g4u(0,1)
           print*, "guu02, guu03=", g4u(0,2), g4u(0,3)
           print*, "guu11, guu12, guu13=", g4u(1,1), g4u(1,2), g4u(1,3)
           print*, "guu22, guu23, guu33=", g4u(2,2), g4u(2,3), g4u(3,3)

           print *,'alpha       = ',alpha
           print *,'piI         = ',piI
           print *,'phiR        = ',phiR
           print *,'d_phiR      = ',d_phiR     

           print*, "dphir=",d_phiR4d(0),d_phiR4d(1)
           print*, d_phiR4d(2),d_phiR4d(3)
           print*, "dphic=",d_phiI4d(0),d_phiI4d(1)
           print*, d_phiI4d(2),d_phiI4d(3)
           print*, "dphiru=",d_phiR4u(0),d_phiR4u(1)
           print*, d_phiR4u(2),d_phiR4u(3)
           print*, "dphicu=",d_phiI4u(0),d_phiI4u(1)
           print*, d_phiI4u(2),d_phiI4u(3)

           print*, "emEu=",emEu(1),emEu(2),emEu(3)
           print*, "emBu=",emBu(1),emBu(2),emBu(3)

           print*, "Tsfu=",Tsfu
           print*, "Temu=",Temu
           print*, "Tu=",Tu

           print *,'rho_AMD  = ',sixteen_pi_G*rho_ADM
           print *,'ham      = ',HamCon
           print *,'momx     = ',MomCon(1)
           print *,'momy     = ',MomCon(2)
           print *,'momz     = ',MomCon(3)
!           print *,' '
!           print *,'gtd(11)  = ',gtd(1,1)
!           print *,'gtd(22)  = ',gtd(2,2)
!           print *,'Atd(11)  = ',Atd(1,1)
!           print *,'Atd(22)  = ',Atd(2,2)
!           print *,'chi      = ',chi
!           print *,'trK      = ',trK
!           print *,'Rt(11)   = ',Rtd(1,1)
!           print *,'Rt(22)   = ',Rtd(2,2)
!           print *,'Rt(33)   = ',Rtd(3,3)
!           print *,'Rphi(11) = ',inv_chi*Rpd(1,1)
!           print *,'Rphi(22) = ',inv_chi*Rpd(2,2)
!           print *,'Rphi(33) = ',inv_chi*Rpd(3,3)
!           print *,'gtu(11)  = ',gtu(1,1)
!           print *,'gtu(22)  = ',gtu(2,2)
!           stop
!
         end if
         end if

        end do
        end do
        end do

        ! - Set the constraints to zero on the boundaries.  The boundary
        ! - conditions are not valid in the grid interior.
        do k = 1, nz
        do j = 1, ny
          do i = 1, 3
            ham(i,j,k)           = 0.0d0
            momx(i,j,k)          = 0.0d0
            momy(i,j,k)          = 0.0d0
            momz(i,j,k)          = 0.0d0
          end do
          do i = nx-2, nx
            ham(i,j,k)           = 0.0d0
            momx(i,j,k)          = 0.0d0
            momy(i,j,k)          = 0.0d0
            momz(i,j,k)          = 0.0d0
          end do
        end do
        end do

        do j = 1, ny
        do i = 1, nx
          do k = 1, 3
            ham(i,j,k)           = 0.0d0
            momx(i,j,k)          = 0.0d0
            momy(i,j,k)          = 0.0d0
            momz(i,j,k)          = 0.0d0
          end do
          do k = nz-2, nz
            ham(i,j,k)           = 0.0d0
            momx(i,j,k)          = 0.0d0
            momy(i,j,k)          = 0.0d0
            momz(i,j,k)          = 0.0d0
          end do
        end do
        end do

        do k = 1, nz
        do i = 1, nx
          do j = 1, 3
            ham(i,j,k)           = 0.0d0
            momx(i,j,k)          = 0.0d0
            momy(i,j,k)          = 0.0d0
            momz(i,j,k)          = 0.0d0
          end do
          do j = ny-2, ny
            ham(i,j,k)           = 0.0d0
            momx(i,j,k)          = 0.0d0
            momy(i,j,k)          = 0.0d0
            momz(i,j,k)          = 0.0d0
          end do
        end do
        end do

        return
      end subroutine
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
