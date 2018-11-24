#include "cctk.h"
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
      subroutine cal_bssn_rhs(dtu, u, v, w, imask,
#include "BR_bssnrhs_sub_args.h"
     &                    x1d, y1d, z1d, nx, ny, nz, par)

        use GF
        use params
!       use ieee_arithmetic 
!     &                   isnan=>IEEE_IS_NAN, isfinite=>ieee_is_finite
        implicit none

        type(gridfunction), dimension(NU_G)               :: dtu
        type(gridfunction), dimension(NU_G)               :: u
        type(gridfunction), dimension(NV_G)               :: v
        type(gridfunction), dimension(NW)                 :: w
        CCTK_INT, dimension(nx,ny,nz)                     :: imask
        CCTK_REAL, dimension(NPAR),intent(in)             :: par
        CCTK_INT                                          :: nx, ny, nz
        CCTK_REAL, dimension(nx)                          :: x1d(nx)
        CCTK_REAL, dimension(ny)                          :: y1d(ny)
        CCTK_REAL, dimension(nz)                          :: z1d(nz)
#include "BR_bssnrhs_decl_sub_args.h"

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
!        CCTK_REAL, dimension(:,:,:), pointer :: Tmn00, Tmn01, Tmn02
!        CCTK_REAL, dimension(:,:,:), pointer :: Tmn03, Tmn11, Tmn12
!        CCTK_REAL, dimension(:,:,:), pointer :: Tmn13, Tmn22, Tmn23
!        CCTK_REAL, dimension(:,:,:), pointer :: Tmn33
        CCTK_REAL, dimension(:,:,:), pointer :: gxx 
        CCTK_REAL, dimension(:,:,:), pointer :: gxy 
        CCTK_REAL, dimension(:,:,:), pointer :: gxz
        CCTK_REAL, dimension(:,:,:), pointer :: gyy 
        CCTK_REAL, dimension(:,:,:), pointer :: gyz 
        CCTK_REAL, dimension(:,:,:), pointer :: gzz
        CCTK_REAL, dimension(:,:,:), pointer :: gbx 
        CCTK_REAL, dimension(:,:,:), pointer :: gby 
        CCTK_REAL, dimension(:,:,:), pointer :: gbz

        CCTK_REAL, dimension(:,:,:), pointer :: dt_Alpha3D
        CCTK_REAL, dimension(:,:,:), pointer :: dt_shiftx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_shifty 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_shiftz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtxx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtxy
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtxz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtyy
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtyz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtzz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atxx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atxy
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atxz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atyy
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atyz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atzz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_chi3D
        CCTK_REAL, dimension(:,:,:), pointer :: dt_trK3D
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Gamtx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Gamty
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Gamtz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gbx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gby
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gbz

        CCTK_REAL, dimension(:,:,:), pointer :: Ex
        CCTK_REAL, dimension(:,:,:), pointer :: Ey 
        CCTK_REAL, dimension(:,:,:), pointer :: Ez
        CCTK_REAL, dimension(:,:,:), pointer :: Bx
        CCTK_REAL, dimension(:,:,:), pointer :: By
        CCTK_REAL, dimension(:,:,:), pointer :: Bz
        CCTK_REAL, dimension(:,:,:), pointer :: Psi_em
        CCTK_REAL, dimension(:,:,:), pointer :: Phi_em 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Ex
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Ey
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Ez 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Bx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_By
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Bz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Psi_em
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Phi_em 

        CCTK_REAL, dimension(:,:,:), pointer :: phiR_sf
        CCTK_REAL, dimension(:,:,:), pointer :: phiI_sf 
        CCTK_REAL, dimension(:,:,:), pointer :: piR_sf
        CCTK_REAL, dimension(:,:,:), pointer :: piI_sf 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_phiR
        CCTK_REAL, dimension(:,:,:), pointer :: dt_phiI 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_piR
        CCTK_REAL, dimension(:,:,:), pointer :: dt_piI 

        CCTK_REAL, dimension(:,:,:), pointer :: x3d
        CCTK_REAL, dimension(:,:,:), pointer :: y3d
        CCTK_REAL, dimension(:,:,:), pointer :: z3d

	CCTK_REAL, dimension(:,:,:), pointer :: sdetg

        CCTK_INT                       :: i, j, k, ii
        CCTK_INT                       :: rc, err

        CCTK_REAL                      :: dx, dy, dz
        CCTK_REAL                      :: trK, Alpha, chi, G
        CCTK_REAL                      :: rho_ADM, trPsi1, tr_pT
        CCTK_REAL                      :: Alpha_rhs, trK_rhs, chi_rhs
        CCTK_REAL,  dimension(3)       :: Betau, Bu_rhs, Betau_rhs
        CCTK_REAL,  dimension(3)       :: Jtd_ADM, CalGamt, Gamt 
        CCTK_REAL,  dimension(3)       :: d_Alpha
        CCTK_REAL,  dimension(3)       :: d_chi, d_trK
        CCTK_REAL,  dimension(3)       :: Gamt_rhs
        CCTK_REAL,  dimension(3)       :: Bu
        CCTK_REAL,  dimension(3,3)     :: gd, gu, gtd, gtu
        CCTK_REAL,  dimension(3,3)     :: Atd, Atu, Atud, dd_chi
        CCTK_REAL,  dimension(3,3)     :: Rpd, Rtd, Rpd_1
!        CCTK_REAL,  dimension(3,3)     :: DelDel_Alpha, pTtd_ADM
        CCTK_REAL,  dimension(3,3)     :: pTtd_ADM
        CCTK_REAL,  dimension(3,3)     :: dd_Alpha, d_Gamt 
        CCTK_REAL,  dimension(3,3)     :: d_Betau, Psi1, Psi1TF
        CCTK_REAL,  dimension(3,3)     :: Atd_rhs, gtd_rhs
        CCTK_REAL,  dimension(3,3)     :: Atd_rhs_xxx, gtd_rhs_xxx !maple vars
        CCTK_REAL,  dimension(3,3,3)   :: dd_Betau, d_gtd, d_Atd
        CCTK_REAL,  dimension(3,3,3,3) :: dd_gtd 
        CCTK_REAL,  dimension(0:3,0:3) :: g4d, g4u, Tu, Tsfu, Temu
        CCTK_REAL                      :: detgd, idetgd
        CCTK_REAL                      :: detgtd, idetgtd
        CCTK_REAL,  dimension(3,3,3)   :: Ctd, Ct, Chr
        CCTK_REAL                      :: twelve_Pi_G
        CCTK_REAL                      :: third_trPsi1
!        CCTK_REAL                      :: expphi,  exp2phi,  exp4phi  ! set?
!        CCTK_REAL                      :: expmphi, expm2phi, expm4phi ! set?
        CCTK_REAL                      :: inv_chi, sqrt_chi 
!        CCTK_REAL                      :: inv_sqrt_chi  

        CCTK_REAL,  dimension(3)       ::   emEu,       emBu 
        CCTK_REAL,  dimension(3,3)     :: d_emEu,     d_emBu 
        CCTK_REAL,  dimension(3)       ::   emEu_rhs,   emBu_rhs
        CCTK_REAL                      :: emPsi, emPhi   
        CCTK_REAL                      :: emPsi_rhs, emPhi_rhs   
        CCTK_REAL, dimension(3)        :: d_emPhi, d_emPsi
        CCTK_REAL                      :: divE, divB, Bsq, inv_Bsq 
        CCTK_REAL                      :: div_Beta, Fsq, qE  
        CCTK_REAL                      :: dphiR4sq, dphiI4sq, sfVR, sfVI 
        CCTK_REAL                      :: dVdphiR2, dVdphiI2
        CCTK_REAL                      :: dphiIR4sq, dphiRI4sq
        CCTK_REAL, dimension(3)        :: d_div_Beta, Ju, Jd
        CCTK_REAL,  dimension(0:3,0:3) :: Femd, Femud, Femu, dualFemu


        CCTK_REAL                      :: dphi4sq, phi2, sfV, dVdphi2
        CCTK_REAL                      :: phiR_rhs,phiI_rhs
        CCTK_REAL                      :: piR_rhs,piI_rhs
        CCTK_REAL                      :: phiR,phiI,piR,piI
        CCTK_REAL, dimension(3)        :: d_phiR, d_phiI, d_piR, d_piI
        CCTK_REAL,  dimension(3,3)     :: dd_phiR, dd_phiI
        CCTK_REAL,  dimension(0:3)     :: d_phiR4d, d_phiI4d
        CCTK_REAL,  dimension(0:3)     :: d_phiR4u, d_phiI4u

 
        CCTK_REAL, dimension(3)        :: xpt 
        CCTK_REAL                      :: r_coord
 
        CCTK_REAL                      :: det_gtd, idet_gtd
        CCTK_REAL                      :: det_gtd_to_third
        CCTK_REAL                      :: det_gtd_to_neg_third
        CCTK_REAL                      :: trace_Atd
        CCTK_REAL                      :: neg_one_third_trace_Atd

        CCTK_INT                       :: lambda_1,  lambda_2
        CCTK_INT                       :: lambda_3,  lambda_4
        CCTK_INT                       :: lambda_f0, lambda_f1
        CCTK_INT                       :: lambda_f2, lambda_f3

        CCTK_INT                       :: eta_damping, eta_damping_exp
        CCTK_REAL                      :: R_0
        CCTK_REAL                      :: eta, trK0, feta
        CCTK_REAL                      :: dil_alpha, dil_mass
        CCTK_REAL                      :: axn_alpha, axn_mass
!       maple vars : start (Need to remove)
        CCTK_REAL                      :: t65, t99, t40, t74, t97
        CCTK_REAL                      :: t122, t124, t126, t163, t164
        CCTK_REAL                      :: t165, t169, t175, t198
!       maple vars :end
        CCTK_REAL                      :: kappa_1, kappa_2 !! for the EM damping
        CCTK_REAL                      :: chi_floor 
	CCTK_REAL	               :: chi_min, dissipSF, extradiss
	CCTK_REAL		       :: dissfact

        CCTK_INT                       :: force_free 

        CCTK_INT                       :: evolve_geometry
        CCTK_INT                       :: evolve_em_field
        CCTK_INT                       :: evolve_scalar_field

        CCTK_REAL, parameter    :: twelve = 12.0d0
        CCTK_REAL, parameter    :: six    = 6.0d0
        CCTK_REAL, parameter    :: four   = 4.0d0
        CCTK_REAL, parameter    :: three  = 3.0d0
        CCTK_REAL, parameter    :: two    = 2.0d0
        CCTK_REAL, parameter    :: one    = 1.0d0
        CCTK_REAL, parameter    :: threehalves  = 1.5d0
        CCTK_REAL, parameter    :: threefourths = 0.75d0
        CCTK_REAL, parameter    :: twothirds = 0.6666666666666667d0
        CCTK_REAL, parameter    :: half   = 0.5d0
        CCTK_REAL, parameter    :: third  = 0.3333333333333333d0
        CCTK_REAL, parameter    :: fourth = 0.25d0 
        CCTK_REAL, parameter    :: sixth  = 0.1666666666666667d0
        CCTK_REAL, parameter    :: four_pi_G  = 12.5663706143591729d0 
        CCTK_REAL, parameter    :: eight_pi_G = 25.1327412287183459d0 

        character(len=*), parameter    :: FMT3='(3E12.3)'
        character(len=*), parameter    :: FMT6='(6E12.3)'

        logical, parameter             :: ltrace      = .false.
        logical, parameter             :: ltrace2     = .false.
        logical, parameter             :: error_check = .false.

!       This variable should be:  turn_off = .false.
        logical, parameter             :: turn_off    = .false.

        ! include declarations for Maple temp variables
#include "bssn_emtest_temp_vars.h"

        if (turn_off) then
          return
        end if
 
!       print *,'begin: ',w(H_XPHYS)%d(31,23,31), 
!    &                    w(H_YPHYS)%d(31,23,31), 
!    &                    w(H_ZPHYS)%d(31,23,31)
!       print *,'begin: dx_alpha  ',dx_alpha(31,23,31)
!       print *,'begin: dy_alpha  ',dy_alpha(31,23,31)
!       print *,'begin: dz_alpha  ',dz_alpha(31,23,31)

        ! G is Newton"s gravitation constant
        G = 1.0d0
        lambda_1  = NINT(par(P_BSSN_LAMBDA_1))  
        lambda_2  = NINT(par(P_BSSN_LAMBDA_2))  
        lambda_3  = NINT(par(P_BSSN_LAMBDA_3))  
        lambda_4  = NINT(par(P_BSSN_LAMBDA_4))  

        lambda_f0 = par(P_BSSN_LAMBDA_F0)  
        lambda_f1 = par(P_BSSN_LAMBDA_F1)  
        lambda_f2 = par(P_BSSN_LAMBDA_F2)  
        lambda_f3 = par(P_BSSN_LAMBDA_F3)  

        eta_damping     = NINT(par(P_BSSN_ETA_DAMPING))
        eta_damping_exp = NINT(par(P_BSSN_ETA_DAMPING_EXP))

        R_0  = par(P_BSSN_R_0)
        eta  = par(P_BSSN_ETA)
        trK0 = par(P_BSSN_TRK0)

        kappa_1 = par(P_BSSN_KAPPA1)
        kappa_2 = par(P_BSSN_KAPPA2)
	chi_min = par(P_KAPPA_MAX)
	extradiss = PAR(P_EXTRADISS)
	dissfact = par(P_MAGCASE)

        chi_floor = par(P_BSSN_CHI_FLOOR)

        dil_alpha = par(P_DIL_ALPHA)
        dil_mass  = par(P_DIL_MASS)

        axn_alpha = par(P_AXN_ALPHA)
        axn_mass  = par(P_AXN_MASS)

        dx = par(P_DX)
        dy = par(P_DY)
        dz = par(P_DZ)

        evolve_geometry     = nint(par(P_EVOLVE_GEOMETRY))
        evolve_em_field     = nint(par(P_EVOLVE_EM_FIELD))
        evolve_scalar_field = nint(par(P_EVOLVE_SCALAR_FIELD))

        if ( ltrace ) then
            write(0,*)' cal_bssn_rhs:  Checking parameters '
            write(0,*)'           These should be integers: '
            write(0,*)'                lambda_1  = ', lambda_1
            write(0,*)'                lambda_2  = ', lambda_2
            write(0,*)'                lambda_3  = ', lambda_3
            write(0,*)'                lambda_4  = ', lambda_4
            write(0,*)'                lambda_f0 = ', lambda_f0
            write(0,*)'                lambda_f1 = ', lambda_f1
            write(0,*)'           These should be reals: '
            write(0,*)'                eta       = ', eta
            write(0,*)'                trK0      = ', trK0
            write(0,*)'                kappa_1   = ', kappa_1
            write(0,*)'                kappa_2   = ', kappa_1
            write(0,*)'                chi_floor = ', chi_floor
            write(0,*)'    These should be integers: '
            write(0,*)'    evolve_geometry = ', evolve_geometry
            write(0,*)'    evolve_em_field = ', evolve_em_field
            write(0,*)'    evolve_scalar_field = ', evolve_scalar_field
        endif

        rho_ADM = 0.0d0
        tr_pT   = 0.0d0
        Jtd_ADM(1)   = 0.0d0
        Jtd_ADM(2)   = 0.0d0
        Jtd_ADM(3)   = 0.0d0
        pTtd_ADM(1,1)   = 0.0d0
        pTtd_ADM(2,1)   = 0.0d0
        pTtd_ADM(3,1)   = 0.0d0
        pTtd_ADM(1,2)   = 0.0d0
        pTtd_ADM(2,2)   = 0.0d0
        pTtd_ADM(3,2)   = 0.0d0
        pTtd_ADM(1,3)   = 0.0d0
        pTtd_ADM(2,3)   = 0.0d0
        pTtd_ADM(3,3)   = 0.0d0

       Ju(1)   = 0.0d0
       Ju(2)   = 0.0d0
       Ju(3)   = 0.0d0
       Jd(1)   = 0.0d0
       Jd(2)   = 0.0d0
       Jd(3)   = 0.0d0
!       Sd(1,1)   = 0.0d0
!       Sd(2,1)   = 0.0d0
!       Sd(3,1)   = 0.0d0
!       Sd(1,2)   = 0.0d0
!       Sd(2,2)   = 0.0d0
!       Sd(3,2)   = 0.0d0
!       Sd(1,3)   = 0.0d0
!       Sd(2,3)   = 0.0d0
!       Sd(3,3)   = 0.0d0

        x3d         => w(G_XPHYS)%d 
        y3d         => w(G_YPHYS)%d 
        z3d         => w(G_ZPHYS)%d 

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
        gbx         => u(G_GB1)%d
        gby         => u(G_GB2)%d
        gbz         => u(G_GB3)%d

        dt_Alpha3D     => dtu(G_ALPHA)%d
        dt_shiftx      => dtu(G_SHIFT1)%d
        dt_shifty      => dtu(G_SHIFT2)%d
        dt_shiftz      => dtu(G_SHIFT3)%d
        dt_chi3D       => dtu(G_CHI)%d
        dt_trK3D       => dtu(G_TRK)%d
        dt_gtxx        => dtu(G_GT11)%d
        dt_gtxy        => dtu(G_GT12)%d
        dt_gtxz        => dtu(G_GT13)%d
        dt_gtyy        => dtu(G_GT22)%d
        dt_gtyz        => dtu(G_GT23)%d
        dt_gtzz        => dtu(G_GT33)%d
        dt_Atxx        => dtu(G_A11)%d
        dt_Atxy        => dtu(G_A12)%d
        dt_Atxz        => dtu(G_A13)%d
        dt_Atyy        => dtu(G_A22)%d
        dt_Atyz        => dtu(G_A23)%d
        dt_Atzz        => dtu(G_A33)%d
        dt_Gamtx       => dtu(G_GAM1)%d
        dt_Gamty       => dtu(G_GAM2)%d
        dt_Gamtz       => dtu(G_GAM3)%d
        dt_gbx         => dtu(G_GB1)%d
        dt_gby         => dtu(G_GB2)%d
        dt_gbz         => dtu(G_GB3)%d

        gxx         => v(G_G11)%d
        gxy         => v(G_G12)%d
        gxz         => v(G_G13)%d
        gyy         => v(G_G22)%d
        gyz         => v(G_G23)%d
        gzz         => v(G_G33)%d
	sdetg	    => v(G_SDETG)%d

        ! EM fields 
        Ex          =>  u(G_EX)%d 
        Ey          =>  u(G_EY)%d 
        Ez          =>  u(G_EZ)%d 
        Bx          =>  u(G_BX)%d 
        By          =>  u(G_BY)%d 
        Bz          =>  u(G_BZ)%d 
        Psi_em      =>  u(G_PSI_EM)%d 
        Phi_em      =>  u(G_PHI_EM)%d 
        
        dt_Ex          =>  dtu(G_EX)%d 
        dt_Ey          =>  dtu(G_EY)%d 
        dt_Ez          =>  dtu(G_EZ)%d 
        dt_Bx          =>  dtu(G_BX)%d 
        dt_By          =>  dtu(G_BY)%d 
        dt_Bz          =>  dtu(G_BZ)%d 
        dt_Psi_em      =>  dtu(G_PSI_EM)%d 
        dt_Phi_em      =>  dtu(G_PHI_EM)%d 

        ! Scalar fields 
        phiR_sf      =>  u(G_PHIR)%d 
        phiI_sf      =>  u(G_PHII)%d 
        piR_sf       =>  u(G_PIR)%d 
        piI_sf       =>  u(G_PII)%d 

        dt_phiR     =>  dtu(G_PHIR)%d 
        dt_phiI     =>  dtu(G_PHII)%d 
        dt_piR      =>  dtu(G_PIR)%d 
        dt_piI      =>  dtu(G_PII)%d 

        
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx

          xpt(1) = x3d(i,j,k)
          xpt(2) = y3d(i,j,k)
          xpt(3) = z3d(i,j,k)

        ! Construct a radial coordinate from the Cartesian coordinates 
        ! in the usual way.  We need this if we damp the eta factor in 
        ! the gauge conditions.  (See arxiv:1003.0859) 
          r_coord = sqrt(    xpt(1)*xpt(1)
     &                     + xpt(2)*xpt(2)
     &                     + xpt(3)*xpt(3)  )

! define somehow the value of feta, the damped eta
	  feta = eta
	  if(r_coord.ge.R_0) then
             feta=eta*(R_0/r_coord)**eta_damping_exp
             kappa_1 = par(P_BSSN_KAPPA1)*(R_0/r_coord)**eta_damping_exp
             kappa_2 = par(P_BSSN_KAPPA2)*(R_0/r_coord)**eta_damping_exp
          end if

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



!***********************************************************************
!          ! Enforce the BSSN constraints, i.e. the unit determinant of 
!          ! gtd and tracelessness of Atd.  
!            det_gtd =   gtd(1,1)*( gtd(2,2)*gtd(3,3) - gtd(2,3)**2)
!     &              - gtd(1,2)**2*gtd(3,3)
!     &              + 2.d0*gtd(1,2)*gtd(1,3)*gtd(2,3)
!     &              - gtd(1,3)**2*gtd(2,2)
!          
!          det_gtd_to_third = det_gtd**(third)
!          det_gtd_to_neg_third = 1.d0 / det_gtd_to_third
!          
!          gtd(1,1) = det_gtd_to_neg_third * gtd(1,1)
!          gtd(1,2) = det_gtd_to_neg_third * gtd(1,2)
!          gtd(1,3) = det_gtd_to_neg_third * gtd(1,3)
!          gtd(2,2) = det_gtd_to_neg_third * gtd(2,2)
!          gtd(2,3) = det_gtd_to_neg_third * gtd(2,3)
!          gtd(3,3) = det_gtd_to_neg_third * gtd(3,3)
!          
!          idet_gtd = 0.1D1/det_gtd
!          gtu(1,1) = idet_gtd*(gtd(2,2)*gtd(3,3)-gtd(2,3)**2)
!          gtu(1,2) = idet_gtd*(-gtd(1,2)*gtd(3,3)+gtd(1,3)*gtd(2,3))
!          gtu(1,3) = idet_gtd*(gtd(1,2)*gtd(2,3)-gtd(1,3)*gtd(2,2))
!          gtu(2,1) = gtu(1,2)
!          gtu(2,2) = idet_gtd*(gtd(1,1)*gtd(3,3)-gtd(1,3)**2)
!          gtu(2,3) = idet_gtd*(-gtd(1,1)*gtd(2,3)+gtd(1,2)*gtd(1,3))
!          gtu(3,1) = gtu(1,3)
!          gtu(3,2) = gtu(2,3)
!          gtu(3,3) = idet_gtd*(gtd(1,1)*gtd(2,2)-gtd(1,2)**2)
!          
!          trace_Atd =   Atd(1,1)*gtu(1,1)
!     &                + Atd(2,2)*gtu(2,2)
!     &                + Atd(3,3)*gtu(3,3)
!     &                + 2.d0 * (   Atd(1,2)*gtu(1,2)
!     &                           + Atd(1,3)*gtu(1,3)
!     &                           + Atd(2,3)*gtu(2,3)  )
!          
!          neg_one_third_trace_Atd = - third * trace_Atd
!            
!          Atd(1,1) = Atd(1,1) + neg_one_third_trace_Atd * gtd(1,1)
!          Atd(1,2) = Atd(1,2) + neg_one_third_trace_Atd * gtd(1,2)
!          Atd(1,3) = Atd(1,3) + neg_one_third_trace_Atd * gtd(1,3)
!          Atd(2,1) = Atd(1,2)
!          Atd(2,2) = Atd(2,2) + neg_one_third_trace_Atd * gtd(2,2)
!          Atd(2,3) = Atd(2,3) + neg_one_third_trace_Atd * gtd(2,3)
!          Atd(3,1) = Atd(1,3)
!          Atd(3,2) = Atd(2,3) 
!          Atd(3,3) = Atd(3,3) + neg_one_third_trace_Atd * gtd(3,3)
!          
!         detgt(i,j,k) = gtd(1,1)*(gtd(2,2)*gtd(3,3) - gtd(2,3)**2) 
!     &                + gtd(1,2)*(gtd(1,3)*gtd(2,3) - gtd(1,2)*gtd(3,3))
!     &                + gtd(1,3)*(gtd(1,2)*gtd(2,3) - gtd(1,3)*gtd(2,2))
!
!         trA(i,j,k) =   Atd(1,1)*gtu(1,1)
!     &                + Atd(2,2)*gtu(2,2)
!     &                + Atd(3,3)*gtu(3,3)
!     &                + 2.d0 * (   Atd(1,2)*gtu(1,2)
!     &                           + Atd(1,3)*gtu(1,3)
!     &                           + Atd(2,3)*gtu(2,3)  )
!
!
!
!***********************************************************************




          trK      = trK3D(i,j,k)
          chi      = chi3D(i,j,k)

! enforce a floor on chi 
          chi = max( chi3D(i,j,k) , chi_floor )

! and take the square root of chi (this is needed in the E&M equations)  
          if ( ltrace ) then
            write(0,*)' cal_bssn_rhs:  Checking that chi is nonnegative'
            if ( chi .lt. 0.d0 ) then
               write(0,*)' chi = ', chi
               write(0,*)' and the square root of chi will be undefined'
            endif
          endif

          sqrt_chi = sqrt( chi )

          Gamt(1) = Gamtx(i,j,k)
          Gamt(2) = Gamty(i,j,k)
          Gamt(3) = Gamtz(i,j,k)

          Alpha    = Alpha3D(i,j,k)
          Betau(1) = shiftx(i,j,k)
          Betau(2) = shifty(i,j,k)
          Betau(3) = shiftz(i,j,k)

          d_Alpha(1)    = dx_alpha(i,j,k)
          d_Alpha(2)    = dy_alpha(i,j,k)
          d_Alpha(3)    = dz_alpha(i,j,k)

          dd_Alpha(1,1) = dxx_alpha(i,j,k)
          dd_Alpha(1,2) = dxy_alpha(i,j,k)
          dd_Alpha(1,3) = dxz_alpha(i,j,k)
          dd_Alpha(2,2) = dyy_alpha(i,j,k)
          dd_Alpha(2,3) = dyz_alpha(i,j,k)
          dd_Alpha(3,3) = dzz_alpha(i,j,k)

          d_chi(1)      = dx_chi(i,j,k)
          d_chi(2)      = dy_chi(i,j,k)
          d_chi(3)      = dz_chi(i,j,k)

          dd_chi(1,1)   = dxx_chi(i,j,k)
          dd_chi(1,2)   = dxy_chi(i,j,k)
          dd_chi(1,3)   = dxz_chi(i,j,k)
          dd_chi(2,2)   = dyy_chi(i,j,k)
          dd_chi(2,3)   = dyz_chi(i,j,k)
          dd_chi(3,3)   = dzz_chi(i,j,k)

          d_trK(1)        = dx_trK(i,j,k)
          d_trK(2)        = dy_trK(i,j,k)
          d_trK(3)        = dz_trK(i,j,k)

          d_Gamt(1,1)   = dx_Gamtx(i,j,k)
          d_Gamt(2,1)   = dy_Gamtx(i,j,k)
          d_Gamt(3,1)   = dz_Gamtx(i,j,k)
          d_Gamt(1,2)   = dx_Gamty(i,j,k)
          d_Gamt(2,2)   = dy_Gamty(i,j,k)
          d_Gamt(3,2)   = dz_Gamty(i,j,k)
          d_Gamt(1,3)   = dx_Gamtz(i,j,k)
          d_Gamt(2,3)   = dy_Gamtz(i,j,k)
          d_Gamt(3,3)   = dz_Gamtz(i,j,k)

          d_Betau(1,1)   = dx_shiftx(i,j,k)
          d_Betau(2,1)   = dy_shiftx(i,j,k)
          d_Betau(3,1)   = dz_shiftx(i,j,k)
          d_Betau(1,2)   = dx_shifty(i,j,k)
          d_Betau(2,2)   = dy_shifty(i,j,k)
          d_Betau(3,2)   = dz_shifty(i,j,k)
          d_Betau(1,3)   = dx_shiftz(i,j,k)
          d_Betau(2,3)   = dy_shiftz(i,j,k)
          d_Betau(3,3)   = dz_shiftz(i,j,k)

          dd_Betau(1,1,1)   = dxx_shiftx(i,j,k)
          dd_Betau(1,2,1)   = dxy_shiftx(i,j,k)
          dd_Betau(1,3,1)   = dxz_shiftx(i,j,k)
          dd_Betau(2,1,1)   = dd_Betau(1,2,1)
          dd_Betau(2,2,1)   = dyy_shiftx(i,j,k)
          dd_Betau(2,3,1)   = dyz_shiftx(i,j,k)
          dd_Betau(3,1,1)   = dd_Betau(1,3,1)
          dd_Betau(3,2,1)   = dd_Betau(2,3,1)
          dd_Betau(3,3,1)   = dzz_shiftx(i,j,k)

          dd_Betau(1,1,2)   = dxx_shifty(i,j,k)
          dd_Betau(1,2,2)   = dxy_shifty(i,j,k)
          dd_Betau(1,3,2)   = dxz_shifty(i,j,k)
          dd_Betau(2,1,2)   = dd_Betau(1,2,2)
          dd_Betau(2,2,2)   = dyy_shifty(i,j,k)
          dd_Betau(2,3,2)   = dyz_shifty(i,j,k)
          dd_Betau(3,1,2)   = dd_Betau(1,3,2)
          dd_Betau(3,2,2)   = dd_Betau(2,3,2)
          dd_Betau(3,3,2)   = dzz_shifty(i,j,k)

          dd_Betau(1,1,3)   = dxx_shiftz(i,j,k)
          dd_Betau(1,2,3)   = dxy_shiftz(i,j,k)
          dd_Betau(1,3,3)   = dxz_shiftz(i,j,k)
          dd_Betau(2,1,3)   = dd_Betau(1,2,3)
          dd_Betau(2,2,3)   = dyy_shiftz(i,j,k)
          dd_Betau(2,3,3)   = dyz_shiftz(i,j,k)
          dd_Betau(3,1,3)   = dd_Betau(1,3,3)
          dd_Betau(3,2,3)   = dd_Betau(2,3,3)
          dd_Betau(3,3,3)   = dzz_shiftz(i,j,k)

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

          Bu(1)                = gbx(i,j,k)
          Bu(2)                = gby(i,j,k)
          Bu(3)                = gbz(i,j,k)
 
          emEu(1)              = Ex(i,j,k) 
          emEu(2)              = Ey(i,j,k) 
          emEu(3)              = Ez(i,j,k) 
        
          emBu(1)              = Bx(i,j,k) 
          emBu(2)              = By(i,j,k) 
          emBu(3)              = Bz(i,j,k) 
        
          emPsi                = Psi_em(i,j,k)   
          emPhi                = Phi_em(i,j,k)   
        
          d_emEu(1,1)             = dx_Ex(i,j,k)  
          d_emEu(1,2)             = dx_Ey(i,j,k)  
          d_emEu(1,3)             = dx_Ez(i,j,k)  
          d_emEu(2,1)             = dy_Ex(i,j,k)  
          d_emEu(2,2)             = dy_Ey(i,j,k)  
          d_emEu(2,3)             = dy_Ez(i,j,k)  
          d_emEu(3,1)             = dz_Ex(i,j,k)  
          d_emEu(3,2)             = dz_Ey(i,j,k)  
          d_emEu(3,3)             = dz_Ez(i,j,k)  
            
          d_emBu(1,1)             = dx_Bx(i,j,k)  
          d_emBu(1,2)             = dx_By(i,j,k)  
          d_emBu(1,3)             = dx_Bz(i,j,k)  
          d_emBu(2,1)             = dy_Bx(i,j,k)  
          d_emBu(2,2)             = dy_By(i,j,k)  
          d_emBu(2,3)             = dy_Bz(i,j,k)  
          d_emBu(3,1)             = dz_Bx(i,j,k)  
          d_emBu(3,2)             = dz_By(i,j,k)  
          d_emBu(3,3)             = dz_Bz(i,j,k)  

          d_emPsi(1)              = dx_Psi_em(i,j,k)
          d_emPsi(2)              = dy_Psi_em(i,j,k)
          d_emPsi(3)              = dz_Psi_em(i,j,k)

          d_emPhi(1)              = dx_Phi_em(i,j,k)
          d_emPhi(2)              = dy_Phi_em(i,j,k)
          d_emPhi(3)              = dz_Phi_em(i,j,k)

          phiR                = phiR_sf(i,j,k)   
          phiI                = phiI_sf(i,j,k)   
          piR                 = piR_sf(i,j,k)   
          piI                 = piI_sf(i,j,k)   

          d_phiR(1)           = dx_phiR(i,j,k)
          d_phiR(2)           = dy_phiR(i,j,k)
          d_phiR(3)           = dz_phiR(i,j,k)

          d_phiI(1)           = dx_phiI(i,j,k)
          d_phiI(2)           = dy_phiI(i,j,k)
          d_phiI(3)           = dz_phiI(i,j,k)

          d_piR(1)           = dx_piR(i,j,k)
          d_piR(2)           = dy_piR(i,j,k)
          d_piR(3)           = dz_piR(i,j,k)

          d_piI(1)           = dx_piI(i,j,k)
          d_piI(2)           = dy_piI(i,j,k)
          d_piI(3)           = dz_piI(i,j,k)

          dd_phiR(1,1)       = dxx_phiR(i,j,k)
          dd_phiR(1,2)       = dxy_phiR(i,j,k)
          dd_phiR(1,3)       = dxz_phiR(i,j,k)
          dd_phiR(2,2)       = dyy_phiR(i,j,k)
          dd_phiR(2,3)       = dyz_phiR(i,j,k)
          dd_phiR(3,3)       = dzz_phiR(i,j,k)

          dd_phiI(1,1)       = dxx_phiI(i,j,k)
          dd_phiI(1,2)       = dxy_phiI(i,j,k)
          dd_phiI(1,3)       = dxz_phiI(i,j,k)
          dd_phiI(2,2)       = dyy_phiI(i,j,k)
          dd_phiI(2,3)       = dyz_phiI(i,j,k)
          dd_phiI(3,3)       = dzz_phiI(i,j,k)
!dissipation for SF
	  dissipSF = dissfact*dx* (dxx_phiR(i,j,k) 
     &           + dyy_phiR(i,j,k)+dzz_phiR(i,j,k)) 
            

#include "bssn_emtest.h" 

          if (evolve_geometry .eq. 1) then
            dt_Alpha3D(i,j,k)   = Alpha_rhs
            dt_shiftx(i,j,k)    = Betau_rhs(1)
            dt_shifty(i,j,k)    = Betau_rhs(2)
            dt_shiftz(i,j,k)    = Betau_rhs(3)

            dt_chi3D(i,j,k)     = chi_rhs
            dt_trK3D(i,j,k)     = trK_rhs

            dt_gtxx(i,j,k)      = gtd_rhs(1,1)
            dt_gtxy(i,j,k)      = gtd_rhs(1,2)
            dt_gtxz(i,j,k)      = gtd_rhs(1,3)
            dt_gtyy(i,j,k)      = gtd_rhs(2,2)
            dt_gtyz(i,j,k)      = gtd_rhs(2,3)
            dt_gtzz(i,j,k)      = gtd_rhs(3,3)

            dt_Atxx(i,j,k)      = Atd_rhs(1,1)
            dt_Atxy(i,j,k)      = Atd_rhs(1,2)
            dt_Atxz(i,j,k)      = Atd_rhs(1,3)
            dt_Atyy(i,j,k)      = Atd_rhs(2,2)
            dt_Atyz(i,j,k)      = Atd_rhs(2,3)
            dt_Atzz(i,j,k)      = Atd_rhs(3,3)

            dt_Gamtx(i,j,k)     = Gamt_rhs(1)
            dt_Gamty(i,j,k)     = Gamt_rhs(2)
            dt_Gamtz(i,j,k)     = Gamt_rhs(3)
  
            dt_gbx(i,j,k)       = Bu_rhs(1)
            dt_gby(i,j,k)       = Bu_rhs(2)
            dt_gbz(i,j,k)       = Bu_rhs(3)
          else
            dt_gtxx(i,j,k) = 0.0d0
            dt_gtxy(i,j,k) = 0.0d0
            dt_gtxz(i,j,k) = 0.0d0
            dt_gtyy(i,j,k) = 0.0d0
            dt_gtyz(i,j,k) = 0.0d0
            dt_gtzz(i,j,k) = 0.0d0
  
            dt_Atxx(i,j,k) = 0.0d0
            dt_Atxy(i,j,k) = 0.0d0
            dt_Atxz(i,j,k) = 0.0d0
            dt_Atyy(i,j,k) = 0.0d0
            dt_Atyz(i,j,k) = 0.0d0
            dt_Atzz(i,j,k) = 0.0d0
  
            dt_trK3D(i,j,k) = 0.0d0
            dt_chi3D(i,j,k) = 0.0d0
  
            dt_Gamtx(i,j,k) = 0.0d0
            dt_Gamty(i,j,k) = 0.0d0
            dt_Gamtz(i,j,k) = 0.0d0
            
            dt_Alpha3D(i,j,k) = 0.0d0
            dt_shiftx(i,j,k)  = 0.0d0
            dt_shifty(i,j,k)  = 0.0d0
            dt_shiftz(i,j,k)  = 0.0d0
            
            dt_gbx(i,j,k) = 0.0d0
            dt_gby(i,j,k) = 0.0d0
            dt_gbz(i,j,k) = 0.0d0
          end if



          if (evolve_em_field .eq. 1) then
            !load charge and current
	    v(G_J1)%d(i,j,k) = Ju(1)
            v(G_J2)%d(i,j,k) = Ju(2)
            v(G_J3)%d(i,j,k) = Ju(3)
            v(G_charge)%d(i,j,k) = qE

            dt_Ex(i,j,k)        = emEu_rhs(1) 
            dt_Ey(i,j,k)        = emEu_rhs(2) 
            dt_Ez(i,j,k)        = emEu_rhs(3) 
  
            dt_Bx(i,j,k)        = emBu_rhs(1) 
            dt_By(i,j,k)        = emBu_rhs(2) 
            dt_Bz(i,j,k)        = emBu_rhs(3) 

            dt_Psi_em(i,j,k)    = emPsi_rhs 
            dt_Phi_em(i,j,k)    = emPhi_rhs 
            if (.false. ) then
               ! avoid divergence cleaning at the center of BHs:
               !dt_Psi_em(i,j,k)    = emPsi_rhs * alpha**2
               !dt_Phi_em(i,j,k)    = emPhi_rhs * alpha**2
               dt_Psi_em(i,j,k)    = emPsi_rhs * chi**2
               dt_Phi_em(i,j,k)    = emPhi_rhs * chi**2
            end if

          else
            dt_Ex(i,j,k)        = 0.0d0
            dt_Ey(i,j,k)        = 0.0d0
            dt_Ez(i,j,k)        = 0.0d0
  
            dt_Bx(i,j,k)        = 0.0d0
            dt_By(i,j,k)        = 0.0d0
            dt_Bz(i,j,k)        = 0.0d0

            dt_Psi_em(i,j,k)    = 0.0d0
            dt_Phi_em(i,j,k)    = 0.0d0
          end if



          if (evolve_scalar_field .eq. 1) then
!add dissipation a la KO 1st order to phi
            if(extradiss.gt.0.5.and.chi.lt.chi_min) then
	     phiR_rhs = phiR_rhs+dissipSF*(1.0d0-chi/chi_min) 
!	     phiI_rhs = phiI_rhs*(sdetg_max/sdetg(i,j,k)) 
	    end if

            dt_phiR(i,j,k)    = phiR_rhs
            dt_phiI(i,j,k)    = phiI_rhs 
            dt_piR(i,j,k)     = piR_rhs 
            dt_piI(i,j,k)     = piI_rhs 
          else
            dt_phiR(i,j,k)    = 0.0d0
            dt_phiI(i,j,k)    = 0.0d0 
            dt_piR(i,j,k)     = 0.0d0 
            dt_piI(i,j,k)     = 0.0d0 
          end if

         if ( ltrace ) then
         if (k .eq. 31 .and. i .eq. 31 .and. j .eq. 31) then
           print *,'000000000000000000000000000000000000000000'
           print *,'000000000000000000000000000000000000000000'
           print *,'000000000000000000000000000000000000000000'
           print *,'000000000000000000000000000000000000000000'

           print *,'! x = ',w(H_XPHYS)%d(i,j,k)
           print *,'! y = ',w(H_YPHYS)%d(i,j,k)
           print *,'! z = ',w(H_ZPHYS)%d(i,j,k)

           print *,'! dt_phiR = ',dt_phiR(i,j,k)
           print *,'! dt_phiI  = ',dt_phiI(i,j,k)
           print *,'! dt_piR = ', dt_piR(i,j,k)
           print *,'! dt_piI  = ',dt_piI(i,j,k)

           print*, "! Ju=",Ju

           print *,'! dt_Ex = ',dt_Ex(i,j,k)
           print *,'! dt_Bx  = ',dt_Bx(i,j,k)
           print *,'! Alpha = ',Alpha
           print *,'! d3_Alpha = ',d_Alpha(3)
           print *,'! trK      = ',trK
           print *,'! Ju(1)    = ',Ju(1)

           print *, "d_phiR(1) = ", d_phiR(1)
           print *, "d_phiR(2) = ", d_phiR(2)
           print *, "d_phiR(3) = ", d_phiR(3)

           print *, "d_phiI(1) = ", d_phiI(1)
           print *, "d_phiI(2) = ", d_phiI(2)
           print *, "d_phiI(3) = ", d_phiI(3)

           print *,'d_phiR4d(0) = ', d_phiR4d(0)
           print *,'d_phiR4d(1) = ', d_phiR4d(1)
           print *,'d_phiR4d(2) = ', d_phiR4d(2)
           print *,'d_phiR4d(3) = ', d_phiR4d(3)


           print *,'d_phiR4u(0) = ', d_phiR4u(0)
           print *,'d_phiR4u(1) = ', d_phiR4u(1)
           print *,'d_phiR4u(2) = ', d_phiR4u(2)
           print *,'d_phiR4u(3) = ', d_phiR4u(3)

           print *,'d_phiI4u(0) = ', d_phiI4u(0)
           print *,'d_phiI4u(1) = ', d_phiI4u(1)
           print *,'d_phiI4u(2) = ', d_phiI4u(2)
           print *,'d_phiI4u(3) = ', d_phiI4u(3)


           print *,'guu(00) = ', g4u(0,0)
           print *,'guu(01) = ', g4u(0,1)
           print *,'guu(02) = ', g4u(0,2)
           print *,'guu(03) = ', g4u(0,3)
           print *,'guu(10) = ', g4u(1,0)
           print *,'guu(11) = ', g4u(1,1)
           print *,'guu(12) = ', g4u(1,2)
           print *,'guu(13) = ', g4u(1,3)
           print *,'guu(20) = ', g4u(2,0)
           print *,'guu(21) = ', g4u(2,1)
           print *,'guu(22) = ', g4u(2,2)
           print *,'guu(23) = ', g4u(2,3)
           print *,'guu(30) = ', g4u(3,0)
           print *,'guu(31) = ', g4u(3,1)
           print *,'guu(32) = ', g4u(3,2)
           print *,'guu(33) = ', g4u(3,3)


           print *,'Tsf(00) = ', Tsfu(0,0)
           print *,'Tsf(01) = ', Tsfu(0,1)
           print *,'Tsf(02) = ', Tsfu(0,2)
           print *,'Tsf(03) = ', Tsfu(0,3)
           print *,'Tsf(10) = ', Tsfu(1,0)
           print *,'Tsf(11) = ', Tsfu(1,1)
           print *,'Tsf(12) = ', Tsfu(1,2)
           print *,'Tsf(13) = ', Tsfu(1,3)
           print *,'Tsf(20) = ', Tsfu(2,0)
           print *,'Tsf(21) = ', Tsfu(2,1)
           print *,'Tsf(22) = ', Tsfu(2,2)
           print *,'Tsf(23) = ', Tsfu(2,3)
           print *,'Tsf(30) = ', Tsfu(3,0)
           print *,'Tsf(31) = ', Tsfu(3,1)
           print *,'Tsf(32) = ', Tsfu(3,2)
           print *,'Tsf(33) = ', Tsfu(3,3)

           print *,'Tu(00) = ', Tu(0,0)
           print *,'Tu(01) = ', Tu(0,1)
           print *,'Tu(02) = ', Tu(0,2)
           print *,'Tu(03) = ', Tu(0,3)
           print *,'Tu(10) = ', Tu(1,0)
           print *,'Tu(11) = ', Tu(1,1)
           print *,'Tu(12) = ', Tu(1,2)
           print *,'Tu(13) = ', Tu(1,3)
           print *,'Tu(20) = ', Tu(2,0)
           print *,'Tu(21) = ', Tu(2,1)
           print *,'Tu(22) = ', Tu(2,2)
           print *,'Tu(23) = ', Tu(2,3)
           print *,'Tu(30) = ', Tu(3,0)
           print *,'Tu(31) = ', Tu(3,1)
           print *,'Tu(32) = ', Tu(3,2)
           print *,'Tu(33) = ', Tu(3,3)

           print *,'alp_rhs  = ', Alpha_rhs
           print *,'beta_rhs(1)  = ',Betau_rhs(1)
           print *,'beta_rhs(2)  = ',Betau_rhs(2)
           print *,'beta_rhs(3)  = ',Betau_rhs(3)

           print *,'chi_rhs(1)  = ',chi_rhs
           print *,'trK_rhs(1)  = ',trK_rhs

           print *,'Atd_rhs(11)  = ',Atd_rhs(1,1)
           print *,'Atd_rhs(12)  = ',Atd_rhs(1,2)
           print *,'Atd_rhs(13)  = ',Atd_rhs(1,3)
           print *,'Atd_rhs(22)  = ',Atd_rhs(2,2)
           print *,'Atd_rhs(23)  = ',Atd_rhs(2,3)
           print *,'Atd_rhs(33)  = ',Atd_rhs(3,3)
           print *,'dtg_rhs(11)  = ',gtd_rhs(1,1)
           print *,'dtg_rhs(12)  = ',gtd_rhs(1,2)
           print *,'dtg_rhs(13)  = ',gtd_rhs(1,3)
           print *,'dtg_rhs(22)  = ',gtd_rhs(2,2)
           print *,'dtg_rhs(23)  = ',gtd_rhs(2,3)
           print *,'dtg_rhs(33)  = ',gtd_rhs(3,3)

           print *,'Gamt_rhs(1)  = ',Gamt_rhs(1)
           print *,'Gamt_rhs(2)  = ',Gamt_rhs(2)
           print *,'Gamt_rhs(3)  = ',Gamt_rhs(3)

           print *,'Bu_rhs(1)  = ',Bu_rhs(1)
           print *,'Bu_rhs(2)  = ',Bu_rhs(2)
           print *,'Bu_rhs(3)  = ',Bu_rhs(3)

         end if
         end if

          if (error_check) then
            err = 1
            call check_isfinite(dt_Alpha3D(i,j,k), 'dt_Alpha3d', rc)
            err = err*rc
            call check_isfinite(dt_shiftx(i,j,k), 'dt_shiftx', rc)
            err = err*rc
            call check_isfinite(dt_shifty(i,j,k), 'dt_shifty', rc)
            err = err*rc
            call check_isfinite(dt_shiftz(i,j,k), 'dt_shiftz', rc)
            err = err*rc
  
            call check_isfinite(dt_chi3D(i,j,k), 'dt_chi3D', rc)
            err = err*rc
            call check_isfinite(dt_trK3D(i,j,k), 'dt_trK3D', rc)
            err = err*rc
  
            call check_isfinite(dt_gtxx(i,j,k), 'dt_gtxx', rc)
            err = err*rc
            call check_isfinite(dt_gtxy(i,j,k), 'dt_gtxy', rc)
            err = err*rc
            call check_isfinite(dt_gtxz(i,j,k), 'dt_gtxz', rc)
            err = err*rc
            call check_isfinite(dt_gtyy(i,j,k), 'dt_gtyy', rc)
            err = err*rc
            call check_isfinite(dt_gtyz(i,j,k), 'dt_gtyz', rc)
            err = err*rc
            call check_isfinite(dt_gtzz(i,j,k), 'dt_gtzz', rc)
            err = err*rc
  
            call check_isfinite(dt_Atxx(i,j,k), 'dt_Atxx', rc)
            err = err*rc
            call check_isfinite(dt_Atxy(i,j,k), 'dt_Atxy', rc)
            err = err*rc
            call check_isfinite(dt_Atxz(i,j,k), 'dt_Atxz', rc)
            err = err*rc
            call check_isfinite(dt_Atyy(i,j,k), 'dt_Atyy', rc)
            err = err*rc
            call check_isfinite(dt_Atyz(i,j,k), 'dt_Atyz', rc)
            err = err*rc
            call check_isfinite(dt_Atzz(i,j,k), 'dt_Atzz', rc)
            err = err*rc
  
            call check_isfinite(dt_Gamtx(i,j,k), 'dt_Gamtx', rc)
            err = err*rc
            call check_isfinite(dt_Gamty(i,j,k), 'dt_Gamty', rc)
            err = err*rc
            call check_isfinite(dt_Gamtz(i,j,k), 'dt_Gamtz', rc)
            err = err*rc
  
            call check_isfinite(dt_gbx(i,j,k), 'dt_gbx', rc)
            err = err*rc
            call check_isfinite(dt_gby(i,j,k), 'dt_gby', rc)
            err = err*rc
            call check_isfinite(dt_gbz(i,j,k), 'dt_gbz', rc)
            err = err*rc
  
            call check_isfinite(dt_Ex(i,j,k), 'dt_Ex', rc)
            err = err*rc
            call check_isfinite(dt_Ey(i,j,k), 'dt_Ey', rc)
            err = err*rc
            call check_isfinite(dt_Ez(i,j,k), 'dt_Ez', rc)
            err = err*rc
  
            call check_isfinite(dt_Bx(i,j,k), 'dt_Bx', rc)
            err = err*rc
            call check_isfinite(dt_By(i,j,k), 'dt_By', rc)
            err = err*rc
            call check_isfinite(dt_Bz(i,j,k), 'dt_Bz', rc)
            err = err*rc
  
            call check_isfinite(dt_Psi_em(i,j,k), 'dt_Psi_em', rc)
            err = err*rc
            call check_isfinite(dt_Phi_em(i,j,k), 'dt_Phi_em', rc)
            err = err*rc
  
            if (err .eq. 0) then
              write(*,*)'### Problem with update'
              write(*,*)' location i,j,k: ',i,j,k
              write(*,*)' location nx,ny,nz: ',nx,ny,nz
              write(*,*)' metric row 1', gd(1,1), gd(1,2), gd(1,3)
              write(*,*)' metric row 2', gd(2,1), gd(2,2), gd(2,3)
              write(*,*)' metric row 3', gd(3,1), gd(3,2), gd(3,3)
  
              write(*,*)' gtd ', gtd(1,1), gtd(1,2), gtd(1,3)
              write(*,*)' gtd ', gtd(2,1), gtd(2,2), gtd(2,3)
              write(*,*)' gtd ', gtd(3,1), gtd(3,2), gtd(3,3)
  
              write(*,*)' Atd ', Atd(1,1), Atd(1,2), Atd(1,3)
              write(*,*)' Atd ', Atd(2,1), Atd(2,2), Atd(2,3)
              write(*,*)' Atd ', Atd(3,1), Atd(3,2), Atd(3,3)
  
              write(*,*)' trK,chi', trK, chi
  
              write(*,*)' Gamt ', Gamt(1), Gamt(2), Gamt(3)
  
              write(*,*)' Alpha ', Alpha
              write(*,*)' Betau ', Betau(1), Betau(2), Betau(3)
  
              write(*,*)' Tu ', Tu(0,0), Tu(0,1), Tu(0,2)
              write(*,*)' Tu ', Tu(0,3), Tu(1,0), Tu(1,1)
              write(*,*)' Tu ', Tu(1,2), Tu(1,3), Tu(2,0)
              write(*,*)' Tu ', Tu(2,1), Tu(2,2), Tu(2,3)
              write(*,*)' Tu ', Tu(3,0), Tu(3,1), Tu(3,2)
              write(*,*)' Tu ', Tu(3,3)
  
       
              write(*,*) ' d_Alpha ', d_Alpha(1), d_Alpha(2), d_Alpha(3)
  
              write(*,*) ' dd_Alpha ', dd_Alpha(1,1), dd_Alpha(1,2), 
     &                                 dd_Alpha(1,3)
              write(*,*)' dd_Alpha ', dd_Alpha(2,2), dd_Alpha(2,3),
     &                                dd_Alpha(3,3)
  
              write(*,*)' d_chi ', d_chi(1), d_chi(2), d_chi(3)
  
              write(*,*)' dd_chi ', dd_chi(1,1),dd_chi(1,2), dd_chi(1,3)
              write(*,*)' dd_chi ', dd_chi(2,2),dd_chi(2,3), dd_chi(3,3)

              write(*,*)' d_trK ', d_trK(1), d_trK(2), d_trK(3)

              write(*,*)' d_Gamt ', d_Gamt(1,1),d_Gamt(2,1), d_Gamt(3,1)
              write(*,*)' d_Gamt ', d_Gamt(1,2),d_Gamt(2,2), d_Gamt(3,2)
              write(*,*)' d_Gamt ', d_Gamt(1,3),d_Gamt(2,3), d_Gamt(3,3)

              write(*,*)' d_Betau 1'
              write(*,FMT3) d_Betau(1,1), d_Betau(2,1), d_Betau(3,1)
              write(*,*)' d_Betau 2'
              write(*,FMT3) d_Betau(1,2), d_Betau(2,2), d_Betau(3,2)
              write(*,*)' d_Betau 3'
              write(*,FMT3) d_Betau(1,3), d_Betau(2,3), d_Betau(3,3)

              write(*,*)' dd_Betau 1'
              write(*,FMT6) dd_Betau(1,1,1), dd_Betau(1,2,1),
     &                    dd_Betau(1,3,1), dd_Betau(2,1,1),
     &                    dd_Betau(2,2,1), dd_Betau(2,3,1)

              write(*,FMT3) dd_Betau(3,1,1), dd_Betau(3,2,1),
     &                    dd_Betau(3,3,1)

              write(*,*)' dd_Betau 2'
              write(*,FMT6) dd_Betau(1,1,2), dd_Betau(1,2,2),
     &                    dd_Betau(1,3,2), dd_Betau(2,1,2),
     &                    dd_Betau(2,2,2), dd_Betau(2,3,2)
              write(*,FMT3) dd_Betau(3,1,2), dd_Betau(3,2,2),
     &                    dd_Betau(3,3,2)

              write(*,*)' dd_Betau 2'
              write(*,FMT6) dd_Betau(1,1,3), dd_Betau(1,2,3),
     &                    dd_Betau(1,3,3), dd_Betau(2,1,3),
     &                    dd_Betau(2,2,3), dd_Betau(2,3,3)
              write(*,FMT3) dd_Betau(3,1,3), dd_Betau(3,2,3),
     &                    dd_Betau(3,3,3)

              write(*,*)' d_gtd 11'
              write(*,FMT3) d_gtd(1,1,1), d_gtd(2,1,1), d_gtd(3,1,1) 
              write(*,*)' d_gtd 12'
              write(*,FMT3) d_gtd(1,1,2), d_gtd(2,1,2), d_gtd(3,1,2)
              write(*,*)' d_gtd 13'
              write(*,FMT3) d_gtd(1,1,3), d_gtd(2,1,3), d_gtd(3,1,3)
              write(*,*)' d_gtd 21'
              write(*,FMT3) d_gtd(1,2,1), d_gtd(2,2,1), d_gtd(3,2,1)
              write(*,*)' d_gtd 22'
              write(*,FMT3) d_gtd(1,2,2), d_gtd(2,2,2), d_gtd(3,2,2)
              write(*,*)' d_gtd 23'
              write(*,FMT3) d_gtd(1,2,3), d_gtd(2,2,3), d_gtd(3,2,3)
              write(*,*)' d_gtd 31'
              write(*,FMT3) d_gtd(1,3,1), d_gtd(2,3,1), d_gtd(3,3,1)
              write(*,*)' d_gtd 32'
              write(*,FMT3) d_gtd(1,3,2), d_gtd(2,3,2), d_gtd(3,3,2)
              write(*,*)' d_gtd 33'
              write(*,FMT3) d_gtd(1,3,3), d_gtd(2,3,3), d_gtd(3,3,3)

              write(*,*)' d_Atd 11'
              write(*,FMT3) d_Atd(1,1,1), d_Atd(2,1,1), d_Atd(3,1,1) 
              write(*,*)' d_Atd 12'
              write(*,FMT3) d_Atd(1,1,2), d_Atd(2,1,2), d_Atd(3,1,2)
              write(*,*)' d_Atd 13'
              write(*,FMT3) d_Atd(1,1,3), d_Atd(2,1,3), d_Atd(3,1,3)
              write(*,*)' d_Atd 21'
              write(*,FMT3) d_Atd(1,2,1), d_Atd(2,2,1), d_Atd(3,2,1)
              write(*,*)' d_Atd 22'
              write(*,FMT3) d_Atd(1,2,2), d_Atd(2,2,2), d_Atd(3,2,2)
              write(*,*)' d_Atd 23'
              write(*,FMT3) d_Atd(1,2,3), d_Atd(2,2,3), d_Atd(3,2,3)
              write(*,*)' d_Atd 31'
              write(*,FMT3) d_Atd(1,3,1), d_Atd(2,3,1), d_Atd(3,3,1)
              write(*,*)' d_Atd 32'
              write(*,FMT3) d_Atd(1,3,2), d_Atd(2,3,2), d_Atd(3,3,2)
              write(*,*)' d_Atd 33'
              write(*,FMT3) d_Atd(1,3,3), d_Atd(2,3,3), d_Atd(3,3,3)

              write(*,*)' dd_gtd 11'
              write(*,FMT6) dd_gtd(1,1,1,1), dd_gtd(1,2,1,1),
     &                    dd_gtd(1,3,1,1), dd_gtd(2,2,1,1),
     &                    dd_gtd(2,3,1,1), dd_gtd(3,3,1,1)

              write(*,*)' dd_gtd 12'
              write(*,FMT6) dd_gtd(1,1,1,2), dd_gtd(1,2,1,2),
     &                    dd_gtd(1,3,1,2), dd_gtd(2,2,1,2),
     &                    dd_gtd(2,3,1,2), dd_gtd(3,3,1,2)

              write(*,*)' dd_gtd 13'
              write(*,FMT6) dd_gtd(1,1,1,3), dd_gtd(1,2,1,3),
     &                    dd_gtd(1,3,1,3), dd_gtd(2,2,1,3),
     &                    dd_gtd(2,3,1,3), dd_gtd(3,3,1,3)

              write(*,*)' dd_gtd 22'
              write(*,FMT6) dd_gtd(1,1,2,2), dd_gtd(1,2,2,2),
     &                    dd_gtd(1,3,2,2), dd_gtd(2,2,2,2),
     &                    dd_gtd(2,3,2,2), dd_gtd(3,3,2,2)

              write(*,*)' dd_gtd 23'
              write(*,FMT6) dd_gtd(1,1,2,3), dd_gtd(1,2,2,3),
     &                    dd_gtd(1,3,2,3), dd_gtd(2,2,2,3),
     &                    dd_gtd(2,3,2,3), dd_gtd(3,3,2,3)

              write(*,*)' dd_gtd 33'
              write(*,FMT6) dd_gtd(1,1,3,3), dd_gtd(1,2,3,3),
     &                    dd_gtd(1,3,3,3), dd_gtd(2,2,3,3),
     &                    dd_gtd(2,3,3,3), dd_gtd(3,3,3,3)

              write(*,*)' Bu'
              write(*,FMT3) Bu(1), Bu(2), Bu(3)

            end if
          end if

          if (ltrace2) then
          if (i .eq. 11 .and. j .eq. 11 .and. k .eq. 11) then
            print *,'gtu(11)  = ',gtu(1,1)
            print *,'gtu(22)  = ',gtu(2,2)
            print *,'d_Gamt(11)   = ',d_Gamt(1,1)
            print *,'d_Gamt(21)   = ',d_Gamt(2,1)
            print *,'d_Gamt(22)   = ',d_Gamt(2,2)
            print *,'d_Gamt(33)   = ',d_Gamt(3,3)
            print *,'dd_gtd(1111) = ',dd_gtd(1,1,1,1)
            print *,'dd_gtd(1122) = ',dd_gtd(1,1,2,2)
            print *,'dd_gtd(3333) = ',dd_gtd(3,3,3,3)
            print *,'Ct(111)      = ',0.5d0*Ct(1,1,1)
            print *,'Ct(112)      = ',0.5d0*Ct(1,1,2)
            print *,'Ct(211)      = ',0.5d0*Ct(2,1,1)
            print *,'Ct(222)      = ',0.5d0*Ct(2,2,2)
            print *,'Ctd(222)     = ',0.5d0*Ctd(2,2,2)
            print *,'rho_AMD      = ',rho_ADM
            print *,'Atu(11)      = ',Atu(1,1)
            print *,'Atu(22)      = ',Atu(2,2)
            print *,'Atu(33)      = ',Atu(3,3)
            print *,'Atd_rhs(11)  = ',Atd_rhs(1,1)
            print *,'Atd_rhs(12)  = ',Atd_rhs(1,2)
            print *,'Atd_rhs(22)  = ',Atd_rhs(2,2)
            print *,'Atd_rhs(33)  = ',Atd_rhs(3,3)
            print *,'Gamt_rhs(1)  = ',Gamt_rhs(1)
            print *,'Gamt_rhs(2)  = ',Gamt_rhs(2)
            print *,'Gamt_rhs(3)  = ',Gamt_rhs(3)
            do ii = 1, NU_G
              print *,'RHS:  ',ii , dtu(ii)%d(i,j,k)
            end do
!           stop
          end if
          end if


        end do
        end do
        end do



        return
      end  subroutine

!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
      subroutine cal_bssn_rhs_lie(dtu, u, v, w, imask,
#include "BR_bssnlie_sub_args.h"
     &                    nx, ny, nz, par)
        use GF
        use params
        implicit none

        type(gridfunction), dimension(NU_G)               :: dtu
        type(gridfunction), dimension(NU_G)               :: u
        type(gridfunction), dimension(NV_G)               :: v
        type(gridfunction), dimension(NW)                 :: w
        CCTK_INT, dimension(nx,ny,nz)                     :: imask
        CCTK_REAL, dimension(NPAR),intent(in)             :: par
        CCTK_INT                                          :: nx, ny, nz

#include "BR_bssnlie_decl_sub_args.h"

        ! local vars
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Alpha3D
        CCTK_REAL, dimension(:,:,:), pointer :: dt_shiftx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_shifty
        CCTK_REAL, dimension(:,:,:), pointer :: dt_shiftz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtxx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtxy
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtxz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtyy
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtyz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gtzz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atxx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atxy
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atxz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atyy
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atyz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Atzz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_chi3D
        CCTK_REAL, dimension(:,:,:), pointer :: dt_trK3D
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Gamtx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Gamty
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Gamtz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gbx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gby
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gbz

        CCTK_REAL, dimension(:,:,:), pointer :: shiftx
        CCTK_REAL, dimension(:,:,:), pointer :: shifty
        CCTK_REAL, dimension(:,:,:), pointer :: shiftz
 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Ex
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Ey
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Ez
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Bx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_By
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Bz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Psi_em
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Phi_em

        CCTK_REAL, dimension(:,:,:), pointer :: dt_phiR
        CCTK_REAL, dimension(:,:,:), pointer :: dt_phiI
        CCTK_REAL, dimension(:,:,:), pointer :: dt_piR
        CCTK_REAL, dimension(:,:,:), pointer :: dt_piI

        CCTK_REAL,  dimension(3)       :: Betau
        CCTK_REAL,  dimension(3,3)     :: gtd_rhs, Atd_rhs 
        CCTK_REAL                      :: chi_rhs, trK_rhs, Alpha_rhs
        CCTK_REAL,  dimension(3)       :: Bu_rhs, Betau_rhs, Gamt_rhs
        CCTK_REAL,  dimension(3)       :: adv_d_chi, adv_d_trK
        CCTK_REAL,  dimension(3)       :: adv_d_Alpha
        CCTK_REAL,  dimension(3,3)     :: adv_d_Betau, adv_d_Bu
        CCTK_REAL,  dimension(3,3)     :: adv_d_Gamt  
        CCTK_REAL,  dimension(3,3,3)   :: adv_d_gtd, adv_d_Atd

        CCTK_REAL,  dimension(3,3)     :: adv_d_emEu , adv_d_emBu 
        CCTK_REAL,  dimension(3)       :: adv_d_emPsi, adv_d_emPhi  
        CCTK_REAL,  dimension(3)       :: emEu_rhs, emBu_rhs
        CCTK_REAL                      :: emPhi_rhs, emPsi_rhs

        CCTK_REAL                      :: phiR_rhs, phiI_rhs
        CCTK_REAL                      :: piR_rhs,piI_rhs
!        CCTK_REAL                      :: phiR,phiI,piR,piI
        CCTK_REAL, dimension(3)        :: adv_d_phiR, adv_d_phiI
        CCTK_REAL, dimension(3)        :: adv_d_piR, adv_d_piI


        CCTK_INT                       :: lambda_1,  lambda_2
        CCTK_INT                       :: lambda_3,  lambda_4
        CCTK_INT                       :: evolve_geometry
        CCTK_INT                       :: evolve_em_field
        CCTK_INT                       :: evolve_scalar_field
!       maple work vars : start
        CCTK_REAL                      :: t1, t2, t4, t7  
!       maple work vars : end
        CCTK_INT                       ::  i, j, k, ii
        CCTK_REAL                      ::  dx, dy, dz
        CCTK_INT                       :: rc, err
        logical, parameter             :: error_check = .false.
!       This variable should be:  turn_off = .false.
        logical, parameter             :: turn_off    = .false.
        logical, parameter             :: ltrace      = .false.
        logical, parameter             :: ltrace2     = .false.

        if (turn_off) then
          return
        end if

        if (ltrace2) then
          write(0,*) '### Begin cal_bssn_rhs_lie'
        end if

        dx = par(P_DX)
        dy = par(P_DY)
        dz = par(P_DZ)

        lambda_1 = NINT(par(P_BSSN_LAMBDA_1))
        lambda_2 = NINT(par(P_BSSN_LAMBDA_2))
        lambda_3 = NINT(par(P_BSSN_LAMBDA_3))
        lambda_4 = NINT(par(P_BSSN_LAMBDA_4))

        evolve_geometry     = nint(par(P_EVOLVE_GEOMETRY))
        evolve_em_field     = nint(par(P_EVOLVE_EM_FIELD))
        evolve_scalar_field = nint(par(P_EVOLVE_SCALAR_FIELD))

        shiftx      => u(G_SHIFT1)%d
        shifty      => u(G_SHIFT2)%d
        shiftz      => u(G_SHIFT3)%d

        dt_Alpha3D     => dtu(G_ALPHA)%d
        dt_shiftx      => dtu(G_SHIFT1)%d
        dt_shifty      => dtu(G_SHIFT2)%d
        dt_shiftz      => dtu(G_SHIFT3)%d
        dt_chi3D       => dtu(G_CHI)%d
        dt_trK3D       => dtu(G_TRK)%d
        dt_gtxx        => dtu(G_GT11)%d
        dt_gtxy        => dtu(G_GT12)%d
        dt_gtxz        => dtu(G_GT13)%d
        dt_gtyy        => dtu(G_GT22)%d
        dt_gtyz        => dtu(G_GT23)%d
        dt_gtzz        => dtu(G_GT33)%d
        dt_Atxx        => dtu(G_A11)%d
        dt_Atxy        => dtu(G_A12)%d
        dt_Atxz        => dtu(G_A13)%d
        dt_Atyy        => dtu(G_A22)%d
        dt_Atyz        => dtu(G_A23)%d
        dt_Atzz        => dtu(G_A33)%d
        dt_Gamtx       => dtu(G_GAM1)%d
        dt_Gamty       => dtu(G_GAM2)%d
        dt_Gamtz       => dtu(G_GAM3)%d
        dt_gbx         => dtu(G_GB1)%d
        dt_gby         => dtu(G_GB2)%d
        dt_gbz         => dtu(G_GB3)%d
        dt_Ex          => dtu(G_EX)%d
        dt_Ey          => dtu(G_EY)%d
        dt_Ez          => dtu(G_EZ)%d
        dt_Bx          => dtu(G_BX)%d
        dt_By          => dtu(G_BY)%d
        dt_Bz          => dtu(G_BZ)%d
        dt_Psi_em      => dtu(G_PSI_EM)%d
        dt_Phi_em      => dtu(G_PHI_EM)%d
        dt_phiR        => dtu(G_PHIR)%d
        dt_phiI        => dtu(G_PHII)%d
        dt_piR         => dtu(G_PIR)%d
        dt_piI         => dtu(G_PII)%d

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          !print *,'i,j,k', i, j, k
          Betau(1) = shiftx(i,j,k)
          Betau(2) = shifty(i,j,k)
          Betau(3) = shiftz(i,j,k)

          !print *,'gtd: i,j,k', i, j, k
          gtd_rhs(1,1) = dt_gtxx(i,j,k)
          gtd_rhs(1,2) = dt_gtxy(i,j,k)
          gtd_rhs(1,3) = dt_gtxz(i,j,k)
          gtd_rhs(2,1) = dt_gtxy(i,j,k)
          gtd_rhs(2,2) = dt_gtyy(i,j,k)
          gtd_rhs(2,3) = dt_gtyz(i,j,k)
          gtd_rhs(3,1) = dt_gtxz(i,j,k)
          gtd_rhs(3,2) = dt_gtyz(i,j,k)
          gtd_rhs(3,3) = dt_gtzz(i,j,k)

          Atd_rhs(1,1) = dt_Atxx(i,j,k)
          Atd_rhs(1,2) = dt_Atxy(i,j,k)
          Atd_rhs(1,3) = dt_Atxz(i,j,k)
          Atd_rhs(2,1) = dt_Atxy(i,j,k)
          Atd_rhs(2,2) = dt_Atyy(i,j,k)
          Atd_rhs(2,3) = dt_Atyz(i,j,k)
          Atd_rhs(3,1) = dt_Atxz(i,j,k)
          Atd_rhs(3,2) = dt_Atyz(i,j,k)
          Atd_rhs(3,3) = dt_Atzz(i,j,k)

          trK_rhs      = dt_trK3D(i,j,k)
          chi_rhs      = dt_chi3D(i,j,k)

          Gamt_rhs(1) = dt_Gamtx(i,j,k)
          Gamt_rhs(2) = dt_Gamty(i,j,k)
          Gamt_rhs(3) = dt_Gamtz(i,j,k)

          Alpha_rhs    = dt_Alpha3D(i,j,k)
          Betau_rhs(1) = dt_shiftx(i,j,k)
          Betau_rhs(2) = dt_shifty(i,j,k)
          Betau_rhs(3) = dt_shiftz(i,j,k)

          Bu_rhs(1)    = dt_gbx(i,j,k)
          Bu_rhs(2)    = dt_gby(i,j,k)
          Bu_rhs(3)    = dt_gbz(i,j,k)

          emEu_rhs(1)    = dt_Ex(i,j,k)
          emEu_rhs(2)    = dt_Ey(i,j,k)
          emEu_rhs(3)    = dt_Ez(i,j,k)
          emBu_rhs(1)    = dt_Bx(i,j,k)
          emBu_rhs(2)    = dt_By(i,j,k)
          emBu_rhs(3)    = dt_Bz(i,j,k)
          emPhi_rhs      = dt_Phi_em(i,j,k)
          emPsi_rhs      = dt_Psi_em(i,j,k)

          phiR_rhs      = dt_phiR(i,j,k)
          phiI_rhs      = dt_phiI(i,j,k)
          piR_rhs       = dt_piR(i,j,k)
          piI_rhs       = dt_piI(i,j,k)

          !print *,'Gamt: i,j,k', i, j, k
          adv_d_Gamt(1,1) = adx_Gamtx(i,j,k)
          adv_d_Gamt(1,2) = adx_Gamty(i,j,k)
          adv_d_Gamt(1,3) = adx_Gamtz(i,j,k)
          adv_d_Gamt(2,1) = ady_Gamtx(i,j,k)
          adv_d_Gamt(2,2) = ady_Gamty(i,j,k)
          adv_d_Gamt(2,3) = ady_Gamtz(i,j,k)
          adv_d_Gamt(3,1) = adz_Gamtx(i,j,k)
          adv_d_Gamt(3,2) = adz_Gamty(i,j,k)
          adv_d_Gamt(3,3) = adz_Gamtz(i,j,k)

          adv_d_Atd(1,1,1) = adx_Atxx(i,j,k)
          adv_d_Atd(2,1,1) = ady_Atxx(i,j,k)
          adv_d_Atd(3,1,1) = adz_Atxx(i,j,k)
          adv_d_Atd(1,1,2) = adx_Atxy(i,j,k)
          adv_d_Atd(2,1,2) = ady_Atxy(i,j,k)
          adv_d_Atd(3,1,2) = adz_Atxy(i,j,k)
          adv_d_Atd(1,1,3) = adx_Atxz(i,j,k)
          adv_d_Atd(2,1,3) = ady_Atxz(i,j,k)
          adv_d_Atd(3,1,3) = adz_Atxz(i,j,k)
          adv_d_Atd(1,2,1) = adv_d_Atd(1,1,2)
          adv_d_Atd(2,2,1) = adv_d_Atd(2,1,2)
          adv_d_Atd(3,2,1) = adv_d_Atd(3,1,2)
          adv_d_Atd(1,2,2) = adx_Atyy(i,j,k)
          adv_d_Atd(2,2,2) = ady_Atyy(i,j,k)
          adv_d_Atd(3,2,2) = adz_Atyy(i,j,k)
          adv_d_Atd(1,2,3) = adx_Atyz(i,j,k)
          adv_d_Atd(2,2,3) = ady_Atyz(i,j,k)
          adv_d_Atd(3,2,3) = adz_Atyz(i,j,k)
          adv_d_Atd(1,3,1) = adv_d_Atd(1,1,3)
          adv_d_Atd(2,3,1) = adv_d_Atd(2,1,3)
          adv_d_Atd(3,3,1) = adv_d_Atd(3,1,3)
          adv_d_Atd(1,3,2) = adv_d_Atd(1,2,3)
          adv_d_Atd(2,3,2) = adv_d_Atd(2,2,3)
          adv_d_Atd(3,3,2) = adv_d_Atd(3,2,3)
          adv_d_Atd(1,3,3) = adx_Atzz(i,j,k)
          adv_d_Atd(2,3,3) = ady_Atzz(i,j,k)
          adv_d_Atd(3,3,3) = adz_Atzz(i,j,k)

          adv_d_gtd(1,1,1) = adx_gtxx(i,j,k)
          adv_d_gtd(2,1,1) = ady_gtxx(i,j,k)
          adv_d_gtd(3,1,1) = adz_gtxx(i,j,k)
          adv_d_gtd(1,1,2) = adx_gtxy(i,j,k)
          adv_d_gtd(2,1,2) = ady_gtxy(i,j,k)
          adv_d_gtd(3,1,2) = adz_gtxy(i,j,k)
          adv_d_gtd(1,1,3) = adx_gtxz(i,j,k)
          adv_d_gtd(2,1,3) = ady_gtxz(i,j,k)
          adv_d_gtd(3,1,3) = adz_gtxz(i,j,k)
          adv_d_gtd(1,2,1) = adv_d_gtd(1,1,2)
          adv_d_gtd(2,2,1) = adv_d_gtd(2,1,2)
          adv_d_gtd(3,2,1) = adv_d_gtd(3,1,2)
          adv_d_gtd(1,2,2) = adx_gtyy(i,j,k)
          adv_d_gtd(2,2,2) = ady_gtyy(i,j,k)
          adv_d_gtd(3,2,2) = adz_gtyy(i,j,k)
          adv_d_gtd(1,2,3) = adx_gtyz(i,j,k)
          adv_d_gtd(2,2,3) = ady_gtyz(i,j,k)
          adv_d_gtd(3,2,3) = adz_gtyz(i,j,k)
          adv_d_gtd(1,3,1) = adv_d_gtd(1,1,3)
          adv_d_gtd(2,3,1) = adv_d_gtd(2,1,3)
          adv_d_gtd(3,3,1) = adv_d_gtd(3,1,3)
          adv_d_gtd(1,3,2) = adv_d_gtd(1,2,3)
          adv_d_gtd(2,3,2) = adv_d_gtd(2,2,3)
          adv_d_gtd(3,3,2) = adv_d_gtd(3,2,3)
          adv_d_gtd(1,3,3) = adx_gtzz(i,j,k)
          adv_d_gtd(2,3,3) = ady_gtzz(i,j,k)
          adv_d_gtd(3,3,3) = adz_gtzz(i,j,k)

          adv_d_chi(1)     = adx_chi(i,j,k)
          adv_d_chi(2)     = ady_chi(i,j,k)
          adv_d_chi(3)     = adz_chi(i,j,k)

          adv_d_trK(1)     = adx_trK(i,j,k)
          adv_d_trK(2)     = ady_trK(i,j,k)
          adv_d_trK(3)     = adz_trK(i,j,k)

          adv_d_Alpha(1)   = adx_alpha(i,j,k)
          adv_d_Alpha(2)   = ady_alpha(i,j,k)
          adv_d_Alpha(3)   = adz_alpha(i,j,k)

          !print *,'Betau: i,j,k', i, j, k
          adv_d_Betau(1,1) = adx_shiftx(i,j,k)
          adv_d_Betau(2,1) = ady_shiftx(i,j,k)
          adv_d_Betau(3,1) = adz_shiftx(i,j,k)
          adv_d_Betau(1,2) = adx_shifty(i,j,k)
          adv_d_Betau(2,2) = ady_shifty(i,j,k)
          adv_d_Betau(3,2) = adz_shifty(i,j,k)
          adv_d_Betau(1,3) = adx_shiftz(i,j,k)
          adv_d_Betau(2,3) = ady_shiftz(i,j,k)
          adv_d_Betau(3,3) = adz_shiftz(i,j,k)

          !print *,'Bu: i,j,k', i, j, k
          adv_d_Bu(1,1)    = adx_gbx(i,j,k)
          adv_d_Bu(2,1)    = ady_gbx(i,j,k)
          adv_d_Bu(3,1)    = adz_gbx(i,j,k)
          adv_d_Bu(1,2)    = adx_gby(i,j,k)
          adv_d_Bu(2,2)    = ady_gby(i,j,k)
          adv_d_Bu(3,2)    = adz_gby(i,j,k)
          adv_d_Bu(1,3)    = adx_gbz(i,j,k)
          adv_d_Bu(2,3)    = ady_gbz(i,j,k)
          adv_d_Bu(3,3)    = adz_gbz(i,j,k)

          !print *,'emEu: i,j,k', i, j, k
          adv_d_emEu(1,1)     = adx_Ex(i,j,k)
          adv_d_emEu(1,2)     = adx_Ey(i,j,k)
          adv_d_emEu(1,3)     = adx_Ez(i,j,k)
          adv_d_emEu(2,1)     = ady_Ex(i,j,k)
          adv_d_emEu(2,2)     = ady_Ey(i,j,k)
          adv_d_emEu(2,3)     = ady_Ez(i,j,k)
          adv_d_emEu(3,1)     = adz_Ex(i,j,k)
          adv_d_emEu(3,2)     = adz_Ey(i,j,k)
          adv_d_emEu(3,3)     = adz_Ez(i,j,k)

          !print *,'emBu: i,j,k', i, j, k
          adv_d_emBu(1,1)     = adx_Bx(i,j,k)
          adv_d_emBu(1,2)     = adx_By(i,j,k)
          adv_d_emBu(1,3)     = adx_Bz(i,j,k)
          adv_d_emBu(2,1)     = ady_Bx(i,j,k)
          adv_d_emBu(2,2)     = ady_By(i,j,k)
          adv_d_emBu(2,3)     = ady_Bz(i,j,k)
          adv_d_emBu(3,1)     = adz_Bx(i,j,k)
          adv_d_emBu(3,2)     = adz_By(i,j,k)
          adv_d_emBu(3,3)     = adz_Bz(i,j,k)

          adv_d_emPsi(1)   = adx_Psi_em(i,j,k)
          adv_d_emPsi(2)   = ady_Psi_em(i,j,k)
          adv_d_emPsi(3)   = adz_Psi_em(i,j,k)

          !print *,'emPhi: i,j,k', i, j, k
          adv_d_emPhi(1)   = adx_Phi_em(i,j,k)
          adv_d_emPhi(2)   = ady_Phi_em(i,j,k)
          adv_d_emPhi(3)   = adz_Phi_em(i,j,k)

          ! for the scalar fields
          adv_d_phiR(1)   = adx_phiR(i,j,k)
          adv_d_phiR(2)   = ady_phiR(i,j,k)
          adv_d_phiR(3)   = adz_phiR(i,j,k)

          adv_d_phiI(1)   = adx_phiI(i,j,k)
          adv_d_phiI(2)   = ady_phiI(i,j,k)
          adv_d_phiI(3)   = adz_phiI(i,j,k)

          adv_d_piR(1)   = adx_piR(i,j,k)
          adv_d_piR(2)   = ady_piR(i,j,k)
          adv_d_piR(3)   = adz_piR(i,j,k)

          adv_d_piI(1)   = adx_piI(i,j,k)
          adv_d_piI(2)   = ady_piI(i,j,k)
          adv_d_piI(3)   = adz_piI(i,j,k)


#include "bssn_emtest_adv.h"

          if (evolve_geometry .eq. 1) then
            !print *,'update rhs gt: i,j,k', i, j, k
            dt_gtxx(i,j,k) = gtd_rhs(1,1)
            dt_gtxy(i,j,k) = gtd_rhs(1,2)
            dt_gtxz(i,j,k) = gtd_rhs(1,3)
            dt_gtyy(i,j,k) = gtd_rhs(2,2)
            dt_gtyz(i,j,k) = gtd_rhs(2,3)
            dt_gtzz(i,j,k) = gtd_rhs(3,3)
  
            dt_Atxx(i,j,k) = Atd_rhs(1,1)
            dt_Atxy(i,j,k) = Atd_rhs(1,2)
            dt_Atxz(i,j,k) = Atd_rhs(1,3)
            dt_Atyy(i,j,k) = Atd_rhs(2,2)
            dt_Atyz(i,j,k) = Atd_rhs(2,3)
            dt_Atzz(i,j,k) = Atd_rhs(3,3)
  
            dt_trK3D(i,j,k) = trK_rhs
            dt_chi3D(i,j,k) = chi_rhs
  
            dt_Gamtx(i,j,k) = Gamt_rhs(1)
            dt_Gamty(i,j,k) = Gamt_rhs(2)
            dt_Gamtz(i,j,k) = Gamt_rhs(3)
            
            dt_Alpha3D(i,j,k) = Alpha_rhs
            dt_shiftx(i,j,k)  = Betau_rhs(1)
            dt_shifty(i,j,k)  = Betau_rhs(2)
            dt_shiftz(i,j,k)  = Betau_rhs(3)
            
            !print *,'update rhs gb: i,j,k', i, j, k
            dt_gbx(i,j,k) = Bu_rhs(1)
            dt_gby(i,j,k) = Bu_rhs(2)
            dt_gbz(i,j,k) = Bu_rhs(3)
          else
            dt_gtxx(i,j,k) = 0.0d0
            dt_gtxy(i,j,k) = 0.0d0
            dt_gtxz(i,j,k) = 0.0d0
            dt_gtyy(i,j,k) = 0.0d0
            dt_gtyz(i,j,k) = 0.0d0
            dt_gtzz(i,j,k) = 0.0d0
  
            dt_Atxx(i,j,k) = 0.0d0
            dt_Atxy(i,j,k) = 0.0d0
            dt_Atxz(i,j,k) = 0.0d0
            dt_Atyy(i,j,k) = 0.0d0
            dt_Atyz(i,j,k) = 0.0d0
            dt_Atzz(i,j,k) = 0.0d0
  
            dt_trK3D(i,j,k) = 0.0d0
            dt_chi3D(i,j,k) = 0.0d0
  
            dt_Gamtx(i,j,k) = 0.0d0
            dt_Gamty(i,j,k) = 0.0d0
            dt_Gamtz(i,j,k) = 0.0d0
            
            dt_Alpha3D(i,j,k) = 0.0d0
            dt_shiftx(i,j,k)  = 0.0d0
            dt_shifty(i,j,k)  = 0.0d0
            dt_shiftz(i,j,k)  = 0.0d0
            
            dt_gbx(i,j,k) = 0.0d0
            dt_gby(i,j,k) = 0.0d0
            dt_gbz(i,j,k) = 0.0d0
 
          end if


          if (evolve_em_field .eq. 1) then
            !print *,'update rhs E: i,j,k', i, j, k
            dt_Ex(i,j,k) = emEu_rhs(1)
            dt_Ey(i,j,k) = emEu_rhs(2)
            dt_Ez(i,j,k) = emEu_rhs(3)

            !print *,'update rhs B: i,j,k', i, j, k
            dt_Bx(i,j,k) = emBu_rhs(1)
            dt_By(i,j,k) = emBu_rhs(2)
            dt_Bz(i,j,k) = emBu_rhs(3)

            dt_Phi_em(i,j,k) = emPhi_rhs
            dt_Psi_em(i,j,k) = emPsi_rhs

          else 
            !print *,'update rhs E: i,j,k', i, j, k
            dt_Ex(i,j,k) = 0.0d0
            dt_Ey(i,j,k) = 0.0d0
            dt_Ez(i,j,k) = 0.0d0

            !print *,'update rhs B: i,j,k', i, j, k
            dt_Bx(i,j,k) = 0.0d0
            dt_By(i,j,k) = 0.0d0
            dt_Bz(i,j,k) = 0.0d0

            dt_Phi_em(i,j,k) = 0.0d0
            dt_Psi_em(i,j,k) = 0.0d0

          end if



          if (evolve_scalar_field .eq. 1) then
            !print *,'update rhs phir: i,j,k', i, j, k
            dt_phiR(i,j,k)    = phiR_rhs 
            dt_phiI(i,j,k)    = phiI_rhs 
            dt_piR(i,j,k)     = piR_rhs 
            dt_piI(i,j,k)     = piI_rhs 

          else 
            dt_phiR(i,j,k)    = 0.0d0
            dt_phiI(i,j,k)    = 0.0d0
            dt_piR(i,j,k)     = 0.0d0 
            dt_piI(i,j,k)     = 0.0d0 

          end if



          if (error_check) then
            err = 1
            call check_isfinite(dt_Alpha3D(i,j,k), 'adv dt_Alpha3d', rc)
            err = err*rc
            call check_isfinite(dt_shiftx(i,j,k), 'adv dt_shiftx', rc)
            err = err*rc
            call check_isfinite(dt_shifty(i,j,k), 'adv dt_shifty', rc)
            err = err*rc
            call check_isfinite(dt_shiftz(i,j,k), 'adv dt_shiftz', rc)
            err = err*rc
  
            call check_isfinite(dt_chi3D(i,j,k), 'adv dt_chi3D', rc)
            err = err*rc
            call check_isfinite(dt_trK3D(i,j,k), 'adv dt_trK3D', rc)
            err = err*rc
  
            call check_isfinite(dt_gtxx(i,j,k), 'adv dt_gtxx', rc)
            err = err*rc
            call check_isfinite(dt_gtxy(i,j,k), 'adv dt_gtxy', rc)
            err = err*rc
            call check_isfinite(dt_gtxz(i,j,k), 'adv dt_gtxz', rc)
            err = err*rc
            call check_isfinite(dt_gtyy(i,j,k), 'adv dt_gtyy', rc)
            err = err*rc
            call check_isfinite(dt_gtyz(i,j,k), 'adv dt_gtyz', rc)
            err = err*rc
            call check_isfinite(dt_gtzz(i,j,k), 'adv dt_gtzz', rc)
            err = err*rc
  
            call check_isfinite(dt_Atxx(i,j,k), 'adv dt_Atxx', rc)
            err = err*rc
            call check_isfinite(dt_Atxy(i,j,k), 'adv dt_Atxy', rc)
            err = err*rc
            call check_isfinite(dt_Atxz(i,j,k), 'adv dt_Atxz', rc)
            err = err*rc
            call check_isfinite(dt_Atyy(i,j,k), 'adv dt_Atyy', rc)
            err = err*rc
            call check_isfinite(dt_Atyz(i,j,k), 'adv dt_Atyz', rc)
            err = err*rc
            call check_isfinite(dt_Atzz(i,j,k), 'adv dt_Atzz', rc)
            err = err*rc
  
            call check_isfinite(dt_Gamtx(i,j,k), 'adv dt_Gamtx', rc)
            err = err*rc
            call check_isfinite(dt_Gamty(i,j,k), 'adv dt_Gamty', rc)
            err = err*rc
            call check_isfinite(dt_Gamtz(i,j,k), 'adv dt_Gamtz', rc)
            err = err*rc
  
            call check_isfinite(dt_gbx(i,j,k), 'adv dt_gbx', rc)
            err = err*rc
            call check_isfinite(dt_gby(i,j,k), 'adv dt_gby', rc)
            err = err*rc
            call check_isfinite(dt_gbz(i,j,k), 'adv dt_gbz', rc)
            err = err*rc
  
            if (err .eq. 0) then
              write(*,*)'Problem with update.  Values not finite.'
              write(*,*)' location i,j,k: ',i,j,k
              call my_exit('NAN found in ADVECTIVE bssnrhs.  Die here.')
            end if
          end if

!         if (i .eq. 11 .and. j .eq. 11 .and. k .eq. 11) then
!           print *,'gtu(11)  = ',gtu(1,1)
!           print *,'gtu(22)  = ',gtu(2,2)
!           print *,'d_Gamt(11)   = ',d_Gamt(1,1)
!           print *,'d_Gamt(21)   = ',d_Gamt(2,1)
!           print *,'d_Gamt(22)   = ',d_Gamt(2,2)
!           print *,'d_Gamt(33)   = ',d_Gamt(3,3)
!           print *,'dd_gtd(1111) = ',dd_gtd(1,1,1,1)
!           print *,'dd_gtd(1122) = ',dd_gtd(1,1,2,2)
!           print *,'dd_gtd(3333) = ',dd_gtd(3,3,3,3)
!           print *,'Ct(111)      = ',0.5d0*Ct(1,1,1)
!           print *,'Ct(112)      = ',0.5d0*Ct(1,1,2)
!           print *,'Ct(211)      = ',0.5d0*Ct(2,1,1)
!           print *,'Ct(222)      = ',0.5d0*Ct(2,2,2)
!           print *,'Ctd(222)     = ',0.5d0*Ctd(2,2,2)
!           print *,'rho_AMD      = ',rho_ADM
!           print *,'Atu(11)      = ',Atu(1,1)
!           print *,'Atu(22)      = ',Atu(2,2)
!           print *,'Atu(33)      = ',Atu(3,3)
!           print *,'Atd_rhs(11)  = ',Atd_rhs(1,1)
!           print *,'Atd_rhs(12)  = ',Atd_rhs(1,2)
!           print *,'Atd_rhs(22)  = ',Atd_rhs(2,2)
!           print *,'Atd_rhs(33)  = ',Atd_rhs(3,3)
!           print *,'Gamt_rhs(1)  = ',Gamt_rhs(1)
!           print *,'Gamt_rhs(2)  = ',Gamt_rhs(2)
!           print *,'Gamt_rhs(3)  = ',Gamt_rhs(3)
!           do ii = 1, NU_G
!             print *,'Final RHS:  ',ii , dtu(ii)%d(i,j,k)
!           end do
!           stop
!         end if




        end do
        end do
        end do


        return
      end  subroutine

!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
