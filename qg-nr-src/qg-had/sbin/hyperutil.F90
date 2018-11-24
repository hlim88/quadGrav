      MODULE UTILEQS
      use params
      implicit none
      CONTAINS
      subroutine assign_ptrs_init(gt11,gt12,gt13,gt22,gt23,gt33,A11,A12,
     &A13,A22,A23,A33,chi,trK,Gam1,Gam2,Gam3,alpha,shift1,shift2,shift3,
     &gb1,gb2,gb3,Ex,Ey,Ez,Bx,By,Bz,Phi_em,Psi_em,phiR,phiI,piR,piI,g11,
     &g12,g13,g22,g23,g33,sdetg,rad_exp,psi4R,psi4I,massADM,massBONDI,cu
     &rvature,phi2R,phi2I,phi0R,phi0I,vbr_omega,vbth_omega,charge,J1,J2,
     &J3,poyntingx_dens,poyntingy_dens,poyntingz_dens,uell,ham,momx,momy
     &,momz,div_B,div_E,tr_A,detgt_m1,gamtx_con,gamty_con,gamtz_con,calg
     &amtx_con,calgamty_con,calgamtz_con,cctk_x,cctk_y,cctk_z,r,xphys,yp
     &hys,zphys,mask,wdiss,chr,error,flag,u0, v, w,nx,ny,nz)
      use GF
      implicit none
      integer :: nx,ny,nz
      type(gridfunction), dimension(NU) :: u0
      type(gridfunction), dimension(NV) :: v
      type(gridfunction), dimension(NW) :: w
      real(kind=8) gt11(nx,ny,nz)
      real(kind=8) gt12(nx,ny,nz)
      real(kind=8) gt13(nx,ny,nz)
      real(kind=8) gt22(nx,ny,nz)
      real(kind=8) gt23(nx,ny,nz)
      real(kind=8) gt33(nx,ny,nz)
      real(kind=8) A11(nx,ny,nz)
      real(kind=8) A12(nx,ny,nz)
      real(kind=8) A13(nx,ny,nz)
      real(kind=8) A22(nx,ny,nz)
      real(kind=8) A23(nx,ny,nz)
      real(kind=8) A33(nx,ny,nz)
      real(kind=8) chi(nx,ny,nz)
      real(kind=8) trK(nx,ny,nz)
      real(kind=8) Gam1(nx,ny,nz)
      real(kind=8) Gam2(nx,ny,nz)
      real(kind=8) Gam3(nx,ny,nz)
      real(kind=8) alpha(nx,ny,nz)
      real(kind=8) shift1(nx,ny,nz)
      real(kind=8) shift2(nx,ny,nz)
      real(kind=8) shift3(nx,ny,nz)
      real(kind=8) gb1(nx,ny,nz)
      real(kind=8) gb2(nx,ny,nz)
      real(kind=8) gb3(nx,ny,nz)
      real(kind=8) Ex(nx,ny,nz)
      real(kind=8) Ey(nx,ny,nz)
      real(kind=8) Ez(nx,ny,nz)
      real(kind=8) Bx(nx,ny,nz)
      real(kind=8) By(nx,ny,nz)
      real(kind=8) Bz(nx,ny,nz)
      real(kind=8) Phi_em(nx,ny,nz)
      real(kind=8) Psi_em(nx,ny,nz)
      real(kind=8) phiR(nx,ny,nz)
      real(kind=8) phiI(nx,ny,nz)
      real(kind=8) piR(nx,ny,nz)
      real(kind=8) piI(nx,ny,nz)
      real(kind=8) g11(nx,ny,nz)
      real(kind=8) g12(nx,ny,nz)
      real(kind=8) g13(nx,ny,nz)
      real(kind=8) g22(nx,ny,nz)
      real(kind=8) g23(nx,ny,nz)
      real(kind=8) g33(nx,ny,nz)
      real(kind=8) sdetg(nx,ny,nz)
      real(kind=8) rad_exp(nx,ny,nz)
      real(kind=8) psi4R(nx,ny,nz)
      real(kind=8) psi4I(nx,ny,nz)
      real(kind=8) massADM(nx,ny,nz)
      real(kind=8) massBONDI(nx,ny,nz)
      real(kind=8) curvature(nx,ny,nz)
      real(kind=8) phi2R(nx,ny,nz)
      real(kind=8) phi2I(nx,ny,nz)
      real(kind=8) phi0R(nx,ny,nz)
      real(kind=8) phi0I(nx,ny,nz)
      real(kind=8) vbr_omega(nx,ny,nz)
      real(kind=8) vbth_omega(nx,ny,nz)
      real(kind=8) charge(nx,ny,nz)
      real(kind=8) J1(nx,ny,nz)
      real(kind=8) J2(nx,ny,nz)
      real(kind=8) J3(nx,ny,nz)
      real(kind=8) poyntingx_dens(nx,ny,nz)
      real(kind=8) poyntingy_dens(nx,ny,nz)
      real(kind=8) poyntingz_dens(nx,ny,nz)
      real(kind=8) uell(nx,ny,nz)
      real(kind=8) ham(nx,ny,nz)
      real(kind=8) momx(nx,ny,nz)
      real(kind=8) momy(nx,ny,nz)
      real(kind=8) momz(nx,ny,nz)
      real(kind=8) div_B(nx,ny,nz)
      real(kind=8) div_E(nx,ny,nz)
      real(kind=8) tr_A(nx,ny,nz)
      real(kind=8) detgt_m1(nx,ny,nz)
      real(kind=8) gamtx_con(nx,ny,nz)
      real(kind=8) gamty_con(nx,ny,nz)
      real(kind=8) gamtz_con(nx,ny,nz)
      real(kind=8) calgamtx_con(nx,ny,nz)
      real(kind=8) calgamty_con(nx,ny,nz)
      real(kind=8) calgamtz_con(nx,ny,nz)
      real(kind=8) cctk_x(nx,ny,nz)
      real(kind=8) cctk_y(nx,ny,nz)
      real(kind=8) cctk_z(nx,ny,nz)
      real(kind=8) r(nx,ny,nz)
      real(kind=8) xphys(nx,ny,nz)
      real(kind=8) yphys(nx,ny,nz)
      real(kind=8) zphys(nx,ny,nz)
      real(kind=8) mask(nx,ny,nz)
      real(kind=8) wdiss(nx,ny,nz)
      real(kind=8) chr(nx,ny,nz)
      real(kind=8) error(nx,ny,nz)
      real(kind=8) flag(nx,ny,nz)
      target :: gt11
      target :: gt12
      target :: gt13
      target :: gt22
      target :: gt23
      target :: gt33
      target :: A11
      target :: A12
      target :: A13
      target :: A22
      target :: A23
      target :: A33
      target :: chi
      target :: trK
      target :: Gam1
      target :: Gam2
      target :: Gam3
      target :: alpha
      target :: shift1
      target :: shift2
      target :: shift3
      target :: gb1
      target :: gb2
      target :: gb3
      target :: Ex
      target :: Ey
      target :: Ez
      target :: Bx
      target :: By
      target :: Bz
      target :: Phi_em
      target :: Psi_em
      target :: phiR
      target :: phiI
      target :: piR
      target :: piI
      target :: g11
      target :: g12
      target :: g13
      target :: g22
      target :: g23
      target :: g33
      target :: sdetg
      target :: rad_exp
      target :: psi4R
      target :: psi4I
      target :: massADM
      target :: massBONDI
      target :: curvature
      target :: phi2R
      target :: phi2I
      target :: phi0R
      target :: phi0I
      target :: vbr_omega
      target :: vbth_omega
      target :: charge
      target :: J1
      target :: J2
      target :: J3
      target :: poyntingx_dens
      target :: poyntingy_dens
      target :: poyntingz_dens
      target :: uell
      target :: ham
      target :: momx
      target :: momy
      target :: momz
      target :: div_B
      target :: div_E
      target :: tr_A
      target :: detgt_m1
      target :: gamtx_con
      target :: gamty_con
      target :: gamtz_con
      target :: calgamtx_con
      target :: calgamty_con
      target :: calgamtz_con
      target :: cctk_x
      target :: cctk_y
      target :: cctk_z
      target :: r
      target :: mask
      target :: wdiss
      target :: chr
      target :: error
      target :: flag
      u0(H_GT11)%d => gt11
      u0(H_GT12)%d => gt12
      u0(H_GT13)%d => gt13
      u0(H_GT22)%d => gt22
      u0(H_GT23)%d => gt23
      u0(H_GT33)%d => gt33
      u0(H_A11)%d => A11
      u0(H_A12)%d => A12
      u0(H_A13)%d => A13
      u0(H_A22)%d => A22
      u0(H_A23)%d => A23
      u0(H_A33)%d => A33
      u0(H_CHI)%d => chi
      u0(H_TRK)%d => trK
      u0(H_GAM1)%d => Gam1
      u0(H_GAM2)%d => Gam2
      u0(H_GAM3)%d => Gam3
      u0(H_ALPHA)%d => alpha
      u0(H_SHIFT1)%d => shift1
      u0(H_SHIFT2)%d => shift2
      u0(H_SHIFT3)%d => shift3
      u0(H_GB1)%d => gb1
      u0(H_GB2)%d => gb2
      u0(H_GB3)%d => gb3
      u0(H_EX)%d => Ex
      u0(H_EY)%d => Ey
      u0(H_EZ)%d => Ez
      u0(H_BX)%d => Bx
      u0(H_BY)%d => By
      u0(H_BZ)%d => Bz
      u0(H_PHI_EM)%d => Phi_em
      u0(H_PSI_EM)%d => Psi_em
      u0(H_PHIR)%d => phiR
      u0(H_PHII)%d => phiI
      u0(H_PIR)%d => piR
      u0(H_PII)%d => piI
      u0(H_GT11)%excise_val = 0
      u0(H_GT11)%dissipation = 1
      u0(H_GT11)%diss_factor = 1
      u0(H_GT11)%zsym = 1
      u0(H_GT11)%ysym = 1
      u0(H_GT11)%xsym = 1
      u0(H_GT11)%stiff = 0
      u0(H_GT11)%take_dx = 0
      u0(H_GT11)%take_dy = 0
      u0(H_GT11)%take_dz = 0
      u0(H_GT12)%excise_val = 0
      u0(H_GT12)%dissipation = 1
      u0(H_GT12)%diss_factor = 1
      u0(H_GT12)%zsym = 1
      u0(H_GT12)%ysym = 1
      u0(H_GT12)%xsym = 1
      u0(H_GT12)%stiff = 0
      u0(H_GT12)%take_dx = 0
      u0(H_GT12)%take_dy = 0
      u0(H_GT12)%take_dz = 0
      u0(H_GT13)%excise_val = 0
      u0(H_GT13)%dissipation = 1
      u0(H_GT13)%diss_factor = 1
      u0(H_GT13)%zsym = 1
      u0(H_GT13)%ysym = 1
      u0(H_GT13)%xsym = 1
      u0(H_GT13)%stiff = 0
      u0(H_GT13)%take_dx = 0
      u0(H_GT13)%take_dy = 0
      u0(H_GT13)%take_dz = 0
      u0(H_GT22)%excise_val = 0
      u0(H_GT22)%dissipation = 1
      u0(H_GT22)%diss_factor = 1
      u0(H_GT22)%zsym = 1
      u0(H_GT22)%ysym = 1
      u0(H_GT22)%xsym = 1
      u0(H_GT22)%stiff = 0
      u0(H_GT22)%take_dx = 0
      u0(H_GT22)%take_dy = 0
      u0(H_GT22)%take_dz = 0
      u0(H_GT23)%excise_val = 0
      u0(H_GT23)%dissipation = 1
      u0(H_GT23)%diss_factor = 1
      u0(H_GT23)%zsym = 1
      u0(H_GT23)%ysym = 1
      u0(H_GT23)%xsym = 1
      u0(H_GT23)%stiff = 0
      u0(H_GT23)%take_dx = 0
      u0(H_GT23)%take_dy = 0
      u0(H_GT23)%take_dz = 0
      u0(H_GT33)%excise_val = 0
      u0(H_GT33)%dissipation = 1
      u0(H_GT33)%diss_factor = 1
      u0(H_GT33)%zsym = 1
      u0(H_GT33)%ysym = 1
      u0(H_GT33)%xsym = 1
      u0(H_GT33)%stiff = 0
      u0(H_GT33)%take_dx = 0
      u0(H_GT33)%take_dy = 0
      u0(H_GT33)%take_dz = 0
      u0(H_A11)%excise_val = 0
      u0(H_A11)%dissipation = 1
      u0(H_A11)%diss_factor = 1
      u0(H_A11)%zsym = 1
      u0(H_A11)%ysym = 1
      u0(H_A11)%xsym = 1
      u0(H_A11)%stiff = 0
      u0(H_A11)%take_dx = 0
      u0(H_A11)%take_dy = 0
      u0(H_A11)%take_dz = 0
      u0(H_A12)%excise_val = 0
      u0(H_A12)%dissipation = 1
      u0(H_A12)%diss_factor = 1
      u0(H_A12)%zsym = 1
      u0(H_A12)%ysym = 1
      u0(H_A12)%xsym = 1
      u0(H_A12)%stiff = 0
      u0(H_A12)%take_dx = 0
      u0(H_A12)%take_dy = 0
      u0(H_A12)%take_dz = 0
      u0(H_A13)%excise_val = 0
      u0(H_A13)%dissipation = 1
      u0(H_A13)%diss_factor = 1
      u0(H_A13)%zsym = 1
      u0(H_A13)%ysym = 1
      u0(H_A13)%xsym = 1
      u0(H_A13)%stiff = 0
      u0(H_A13)%take_dx = 0
      u0(H_A13)%take_dy = 0
      u0(H_A13)%take_dz = 0
      u0(H_A22)%excise_val = 0
      u0(H_A22)%dissipation = 1
      u0(H_A22)%diss_factor = 1
      u0(H_A22)%zsym = 1
      u0(H_A22)%ysym = 1
      u0(H_A22)%xsym = 1
      u0(H_A22)%stiff = 0
      u0(H_A22)%take_dx = 0
      u0(H_A22)%take_dy = 0
      u0(H_A22)%take_dz = 0
      u0(H_A23)%excise_val = 0
      u0(H_A23)%dissipation = 1
      u0(H_A23)%diss_factor = 1
      u0(H_A23)%zsym = 1
      u0(H_A23)%ysym = 1
      u0(H_A23)%xsym = 1
      u0(H_A23)%stiff = 0
      u0(H_A23)%take_dx = 0
      u0(H_A23)%take_dy = 0
      u0(H_A23)%take_dz = 0
      u0(H_A33)%excise_val = 0
      u0(H_A33)%dissipation = 1
      u0(H_A33)%diss_factor = 1
      u0(H_A33)%zsym = 1
      u0(H_A33)%ysym = 1
      u0(H_A33)%xsym = 1
      u0(H_A33)%stiff = 0
      u0(H_A33)%take_dx = 0
      u0(H_A33)%take_dy = 0
      u0(H_A33)%take_dz = 0
      u0(H_CHI)%excise_val = 0
      u0(H_CHI)%dissipation = 1
      u0(H_CHI)%diss_factor = 1
      u0(H_CHI)%zsym = 1
      u0(H_CHI)%ysym = 1
      u0(H_CHI)%xsym = 1
      u0(H_CHI)%stiff = 0
      u0(H_CHI)%take_dx = 0
      u0(H_CHI)%take_dy = 0
      u0(H_CHI)%take_dz = 0
      u0(H_TRK)%excise_val = 0
      u0(H_TRK)%dissipation = 1
      u0(H_TRK)%diss_factor = 1
      u0(H_TRK)%zsym = 1
      u0(H_TRK)%ysym = 1
      u0(H_TRK)%xsym = 1
      u0(H_TRK)%stiff = 0
      u0(H_TRK)%take_dx = 0
      u0(H_TRK)%take_dy = 0
      u0(H_TRK)%take_dz = 0
      u0(H_GAM1)%excise_val = 0
      u0(H_GAM1)%dissipation = 1
      u0(H_GAM1)%diss_factor = 1
      u0(H_GAM1)%zsym = 1
      u0(H_GAM1)%ysym = 1
      u0(H_GAM1)%xsym = 1
      u0(H_GAM1)%stiff = 0
      u0(H_GAM1)%take_dx = 0
      u0(H_GAM1)%take_dy = 0
      u0(H_GAM1)%take_dz = 0
      u0(H_GAM2)%excise_val = 0
      u0(H_GAM2)%dissipation = 1
      u0(H_GAM2)%diss_factor = 1
      u0(H_GAM2)%zsym = 1
      u0(H_GAM2)%ysym = 1
      u0(H_GAM2)%xsym = 1
      u0(H_GAM2)%stiff = 0
      u0(H_GAM2)%take_dx = 0
      u0(H_GAM2)%take_dy = 0
      u0(H_GAM2)%take_dz = 0
      u0(H_GAM3)%excise_val = 0
      u0(H_GAM3)%dissipation = 1
      u0(H_GAM3)%diss_factor = 1
      u0(H_GAM3)%zsym = 1
      u0(H_GAM3)%ysym = 1
      u0(H_GAM3)%xsym = 1
      u0(H_GAM3)%stiff = 0
      u0(H_GAM3)%take_dx = 0
      u0(H_GAM3)%take_dy = 0
      u0(H_GAM3)%take_dz = 0
      u0(H_ALPHA)%excise_val = 0
      u0(H_ALPHA)%dissipation = 1
      u0(H_ALPHA)%diss_factor = 1
      u0(H_ALPHA)%zsym = 1
      u0(H_ALPHA)%ysym = 1
      u0(H_ALPHA)%xsym = 1
      u0(H_ALPHA)%stiff = 0
      u0(H_ALPHA)%take_dx = 0
      u0(H_ALPHA)%take_dy = 0
      u0(H_ALPHA)%take_dz = 0
      u0(H_SHIFT1)%excise_val = 0
      u0(H_SHIFT1)%dissipation = 1
      u0(H_SHIFT1)%diss_factor = 1
      u0(H_SHIFT1)%zsym = 1
      u0(H_SHIFT1)%ysym = 1
      u0(H_SHIFT1)%xsym = 1
      u0(H_SHIFT1)%stiff = 0
      u0(H_SHIFT1)%take_dx = 0
      u0(H_SHIFT1)%take_dy = 0
      u0(H_SHIFT1)%take_dz = 0
      u0(H_SHIFT2)%excise_val = 0
      u0(H_SHIFT2)%dissipation = 1
      u0(H_SHIFT2)%diss_factor = 1
      u0(H_SHIFT2)%zsym = 1
      u0(H_SHIFT2)%ysym = 1
      u0(H_SHIFT2)%xsym = 1
      u0(H_SHIFT2)%stiff = 0
      u0(H_SHIFT2)%take_dx = 0
      u0(H_SHIFT2)%take_dy = 0
      u0(H_SHIFT2)%take_dz = 0
      u0(H_SHIFT3)%excise_val = 0
      u0(H_SHIFT3)%dissipation = 1
      u0(H_SHIFT3)%diss_factor = 1
      u0(H_SHIFT3)%zsym = 1
      u0(H_SHIFT3)%ysym = 1
      u0(H_SHIFT3)%xsym = 1
      u0(H_SHIFT3)%stiff = 0
      u0(H_SHIFT3)%take_dx = 0
      u0(H_SHIFT3)%take_dy = 0
      u0(H_SHIFT3)%take_dz = 0
      u0(H_GB1)%excise_val = 0
      u0(H_GB1)%dissipation = 1
      u0(H_GB1)%diss_factor = 1
      u0(H_GB1)%zsym = 1
      u0(H_GB1)%ysym = 1
      u0(H_GB1)%xsym = 1
      u0(H_GB1)%stiff = 0
      u0(H_GB1)%take_dx = 0
      u0(H_GB1)%take_dy = 0
      u0(H_GB1)%take_dz = 0
      u0(H_GB2)%excise_val = 0
      u0(H_GB2)%dissipation = 1
      u0(H_GB2)%diss_factor = 1
      u0(H_GB2)%zsym = 1
      u0(H_GB2)%ysym = 1
      u0(H_GB2)%xsym = 1
      u0(H_GB2)%stiff = 0
      u0(H_GB2)%take_dx = 0
      u0(H_GB2)%take_dy = 0
      u0(H_GB2)%take_dz = 0
      u0(H_GB3)%excise_val = 0
      u0(H_GB3)%dissipation = 1
      u0(H_GB3)%diss_factor = 1
      u0(H_GB3)%zsym = 1
      u0(H_GB3)%ysym = 1
      u0(H_GB3)%xsym = 1
      u0(H_GB3)%stiff = 0
      u0(H_GB3)%take_dx = 0
      u0(H_GB3)%take_dy = 0
      u0(H_GB3)%take_dz = 0
      u0(H_EX)%excise_val = 0
      u0(H_EX)%dissipation = 1
      u0(H_EX)%diss_factor = 1
      u0(H_EX)%zsym = 1
      u0(H_EX)%ysym = 1
      u0(H_EX)%xsym = 1
      u0(H_EX)%stiff = 0
      u0(H_EX)%take_dx = 0
      u0(H_EX)%take_dy = 0
      u0(H_EX)%take_dz = 0
      u0(H_EY)%excise_val = 0
      u0(H_EY)%dissipation = 1
      u0(H_EY)%diss_factor = 1
      u0(H_EY)%zsym = 1
      u0(H_EY)%ysym = 1
      u0(H_EY)%xsym = 1
      u0(H_EY)%stiff = 0
      u0(H_EY)%take_dx = 0
      u0(H_EY)%take_dy = 0
      u0(H_EY)%take_dz = 0
      u0(H_EZ)%excise_val = 0
      u0(H_EZ)%dissipation = 1
      u0(H_EZ)%diss_factor = 1
      u0(H_EZ)%zsym = 1
      u0(H_EZ)%ysym = 1
      u0(H_EZ)%xsym = 1
      u0(H_EZ)%stiff = 0
      u0(H_EZ)%take_dx = 0
      u0(H_EZ)%take_dy = 0
      u0(H_EZ)%take_dz = 0
      u0(H_BX)%excise_val = 0
      u0(H_BX)%dissipation = 1
      u0(H_BX)%diss_factor = 1
      u0(H_BX)%zsym = 1
      u0(H_BX)%ysym = 1
      u0(H_BX)%xsym = 1
      u0(H_BX)%stiff = 0
      u0(H_BX)%take_dx = 0
      u0(H_BX)%take_dy = 0
      u0(H_BX)%take_dz = 0
      u0(H_BY)%excise_val = 0
      u0(H_BY)%dissipation = 1
      u0(H_BY)%diss_factor = 1
      u0(H_BY)%zsym = 1
      u0(H_BY)%ysym = 1
      u0(H_BY)%xsym = 1
      u0(H_BY)%stiff = 0
      u0(H_BY)%take_dx = 0
      u0(H_BY)%take_dy = 0
      u0(H_BY)%take_dz = 0
      u0(H_BZ)%excise_val = 0
      u0(H_BZ)%dissipation = 1
      u0(H_BZ)%diss_factor = 1
      u0(H_BZ)%zsym = 1
      u0(H_BZ)%ysym = 1
      u0(H_BZ)%xsym = 1
      u0(H_BZ)%stiff = 0
      u0(H_BZ)%take_dx = 0
      u0(H_BZ)%take_dy = 0
      u0(H_BZ)%take_dz = 0
      u0(H_PHI_EM)%excise_val = 0
      u0(H_PHI_EM)%dissipation = 1
      u0(H_PHI_EM)%diss_factor = 1
      u0(H_PHI_EM)%zsym = 1
      u0(H_PHI_EM)%ysym = 1
      u0(H_PHI_EM)%xsym = 1
      u0(H_PHI_EM)%stiff = 0
      u0(H_PHI_EM)%take_dx = 0
      u0(H_PHI_EM)%take_dy = 0
      u0(H_PHI_EM)%take_dz = 0
      u0(H_PSI_EM)%excise_val = 0
      u0(H_PSI_EM)%dissipation = 1
      u0(H_PSI_EM)%diss_factor = 1
      u0(H_PSI_EM)%zsym = 1
      u0(H_PSI_EM)%ysym = 1
      u0(H_PSI_EM)%xsym = 1
      u0(H_PSI_EM)%stiff = 0
      u0(H_PSI_EM)%take_dx = 0
      u0(H_PSI_EM)%take_dy = 0
      u0(H_PSI_EM)%take_dz = 0
      u0(H_PHIR)%excise_val = 0
      u0(H_PHIR)%dissipation = 1
      u0(H_PHIR)%diss_factor = 1
      u0(H_PHIR)%zsym = 1
      u0(H_PHIR)%ysym = 1
      u0(H_PHIR)%xsym = 1
      u0(H_PHIR)%stiff = 0
      u0(H_PHIR)%take_dx = 0
      u0(H_PHIR)%take_dy = 0
      u0(H_PHIR)%take_dz = 0
      u0(H_PHII)%excise_val = 0
      u0(H_PHII)%dissipation = 1
      u0(H_PHII)%diss_factor = 1
      u0(H_PHII)%zsym = 1
      u0(H_PHII)%ysym = 1
      u0(H_PHII)%xsym = 1
      u0(H_PHII)%stiff = 0
      u0(H_PHII)%take_dx = 0
      u0(H_PHII)%take_dy = 0
      u0(H_PHII)%take_dz = 0
      u0(H_PIR)%excise_val = 0
      u0(H_PIR)%dissipation = 1
      u0(H_PIR)%diss_factor = 1
      u0(H_PIR)%zsym = 1
      u0(H_PIR)%ysym = 1
      u0(H_PIR)%xsym = 1
      u0(H_PIR)%stiff = 0
      u0(H_PIR)%take_dx = 0
      u0(H_PIR)%take_dy = 0
      u0(H_PIR)%take_dz = 0
      u0(H_PII)%excise_val = 0
      u0(H_PII)%dissipation = 1
      u0(H_PII)%diss_factor = 1
      u0(H_PII)%zsym = 1
      u0(H_PII)%ysym = 1
      u0(H_PII)%xsym = 1
      u0(H_PII)%stiff = 0
      u0(H_PII)%take_dx = 0
      u0(H_PII)%take_dy = 0
      u0(H_PII)%take_dz = 0
      v(H_G11)%d => g11
      v(H_G12)%d => g12
      v(H_G13)%d => g13
      v(H_G22)%d => g22
      v(H_G23)%d => g23
      v(H_G33)%d => g33
      v(H_SDETG)%d => sdetg
      v(H_RAD_EXP)%d => rad_exp
      v(H_PSI4R)%d => psi4R
      v(H_PSI4I)%d => psi4I
      v(H_MASSADM)%d => massADM
      v(H_MASSBONDI)%d => massBONDI
      v(H_CURVATURE)%d => curvature
      v(H_PHI2R)%d => phi2R
      v(H_PHI2I)%d => phi2I
      v(H_PHI0R)%d => phi0R
      v(H_PHI0I)%d => phi0I
      v(H_VBR_OMEGA)%d => vbr_omega
      v(H_VBTH_OMEGA)%d => vbth_omega
      v(H_CHARGE)%d => charge
      v(H_J1)%d => J1
      v(H_J2)%d => J2
      v(H_J3)%d => J3
      v(H_POYNTINGX_DENS)%d => poyntingx_dens
      v(H_POYNTINGY_DENS)%d => poyntingy_dens
      v(H_POYNTINGZ_DENS)%d => poyntingz_dens
      v(H_UELL)%d => uell
      w(H_HAM)%d => ham
      w(H_MOMX)%d => momx
      w(H_MOMY)%d => momy
      w(H_MOMZ)%d => momz
      w(H_DIV_B)%d => div_B
      w(H_DIV_E)%d => div_E
      w(H_TR_A)%d => tr_A
      w(H_DETGT_M1)%d => detgt_m1
      w(H_GAMTX_CON)%d => gamtx_con
      w(H_GAMTY_CON)%d => gamty_con
      w(H_GAMTZ_CON)%d => gamtz_con
      w(H_CALGAMTX_CON)%d => calgamtx_con
      w(H_CALGAMTY_CON)%d => calgamty_con
      w(H_CALGAMTZ_CON)%d => calgamtz_con
      w(H_X)%d => cctk_x
      w(H_Y)%d => cctk_y
      w(H_Z)%d => cctk_z
      w(H_R)%d => r
      w(H_XPHYS)%d => cctk_x
      w(H_YPHYS)%d => cctk_y
      w(H_ZPHYS)%d => cctk_z
      w(H_MASK)%d => mask
      w(H_WDISS)%d => wdiss
      w(H_CHR)%d => chr
      w(H_ERROR)%d => error
      w(H_FLAG)%d => flag
      return
      end subroutine assign_ptrs_init
      subroutine assign_ptrs_fields(gt11, gt11_np1,gt12, gt12_np1,gt13, 
     &gt13_np1,gt22, gt22_np1,gt23, gt23_np1,gt33, gt33_np1,A11, A11_np1
     &,A12, A12_np1,A13, A13_np1,A22, A22_np1,A23, A23_np1,A33, A33_np1,
     &chi, chi_np1,trK, trK_np1,Gam1, Gam1_np1,Gam2, Gam2_np1,Gam3, Gam3
     &_np1,alpha, alpha_np1,shift1, shift1_np1,shift2, shift2_np1,shift3
     &, shift3_np1,gb1, gb1_np1,gb2, gb2_np1,gb3, gb3_np1,Ex, Ex_np1,Ey,
     & Ey_np1,Ez, Ez_np1,Bx, Bx_np1,By, By_np1,Bz, Bz_np1,Phi_em, Phi_em
     &_np1,Psi_em, Psi_em_np1,phiR, phiR_np1,phiI, phiI_np1,piR, piR_np1
     &,piI, piI_np1,g11,g12,g13,g22,g23,g33,sdetg,rad_exp,psi4R,psi4I,ma
     &ssADM,massBONDI,curvature,phi2R,phi2I,phi0R,phi0I,vbr_omega,vbth_o
     &mega,charge,J1,J2,J3,poyntingx_dens,poyntingy_dens,poyntingz_dens,
     &uell,ham,momx,momy,momz,div_B,div_E,tr_A,detgt_m1,gamtx_con,gamty_
     &con,gamtz_con,calgamtx_con,calgamty_con,calgamtz_con,cctk_x,cctk_y
     &,cctk_z,r,xphys,yphys,zphys,mask,wdiss,chr,error,flag,u2, u0, v, w
     &, nx, ny, nz)
      use GF
      implicit none
      integer :: nx, ny, nz
      type(gridfunction), dimension(NU) :: u2, u0
      type(gridfunction), dimension(NW) :: w
      type(gridfunction), dimension(NV) :: v
      real(kind=8) gt11(nx,ny,nz),gt11_np1(nx,ny,nz)
      real(kind=8) gt12(nx,ny,nz),gt12_np1(nx,ny,nz)
      real(kind=8) gt13(nx,ny,nz),gt13_np1(nx,ny,nz)
      real(kind=8) gt22(nx,ny,nz),gt22_np1(nx,ny,nz)
      real(kind=8) gt23(nx,ny,nz),gt23_np1(nx,ny,nz)
      real(kind=8) gt33(nx,ny,nz),gt33_np1(nx,ny,nz)
      real(kind=8) A11(nx,ny,nz),A11_np1(nx,ny,nz)
      real(kind=8) A12(nx,ny,nz),A12_np1(nx,ny,nz)
      real(kind=8) A13(nx,ny,nz),A13_np1(nx,ny,nz)
      real(kind=8) A22(nx,ny,nz),A22_np1(nx,ny,nz)
      real(kind=8) A23(nx,ny,nz),A23_np1(nx,ny,nz)
      real(kind=8) A33(nx,ny,nz),A33_np1(nx,ny,nz)
      real(kind=8) chi(nx,ny,nz),chi_np1(nx,ny,nz)
      real(kind=8) trK(nx,ny,nz),trK_np1(nx,ny,nz)
      real(kind=8) Gam1(nx,ny,nz),Gam1_np1(nx,ny,nz)
      real(kind=8) Gam2(nx,ny,nz),Gam2_np1(nx,ny,nz)
      real(kind=8) Gam3(nx,ny,nz),Gam3_np1(nx,ny,nz)
      real(kind=8) alpha(nx,ny,nz),alpha_np1(nx,ny,nz)
      real(kind=8) shift1(nx,ny,nz),shift1_np1(nx,ny,nz)
      real(kind=8) shift2(nx,ny,nz),shift2_np1(nx,ny,nz)
      real(kind=8) shift3(nx,ny,nz),shift3_np1(nx,ny,nz)
      real(kind=8) gb1(nx,ny,nz),gb1_np1(nx,ny,nz)
      real(kind=8) gb2(nx,ny,nz),gb2_np1(nx,ny,nz)
      real(kind=8) gb3(nx,ny,nz),gb3_np1(nx,ny,nz)
      real(kind=8) Ex(nx,ny,nz),Ex_np1(nx,ny,nz)
      real(kind=8) Ey(nx,ny,nz),Ey_np1(nx,ny,nz)
      real(kind=8) Ez(nx,ny,nz),Ez_np1(nx,ny,nz)
      real(kind=8) Bx(nx,ny,nz),Bx_np1(nx,ny,nz)
      real(kind=8) By(nx,ny,nz),By_np1(nx,ny,nz)
      real(kind=8) Bz(nx,ny,nz),Bz_np1(nx,ny,nz)
      real(kind=8) Phi_em(nx,ny,nz),Phi_em_np1(nx,ny,nz)
      real(kind=8) Psi_em(nx,ny,nz),Psi_em_np1(nx,ny,nz)
      real(kind=8) phiR(nx,ny,nz),phiR_np1(nx,ny,nz)
      real(kind=8) phiI(nx,ny,nz),phiI_np1(nx,ny,nz)
      real(kind=8) piR(nx,ny,nz),piR_np1(nx,ny,nz)
      real(kind=8) piI(nx,ny,nz),piI_np1(nx,ny,nz)
      real(kind=8) g11(nx,ny,nz)
      real(kind=8) g12(nx,ny,nz)
      real(kind=8) g13(nx,ny,nz)
      real(kind=8) g22(nx,ny,nz)
      real(kind=8) g23(nx,ny,nz)
      real(kind=8) g33(nx,ny,nz)
      real(kind=8) sdetg(nx,ny,nz)
      real(kind=8) rad_exp(nx,ny,nz)
      real(kind=8) psi4R(nx,ny,nz)
      real(kind=8) psi4I(nx,ny,nz)
      real(kind=8) massADM(nx,ny,nz)
      real(kind=8) massBONDI(nx,ny,nz)
      real(kind=8) curvature(nx,ny,nz)
      real(kind=8) phi2R(nx,ny,nz)
      real(kind=8) phi2I(nx,ny,nz)
      real(kind=8) phi0R(nx,ny,nz)
      real(kind=8) phi0I(nx,ny,nz)
      real(kind=8) vbr_omega(nx,ny,nz)
      real(kind=8) vbth_omega(nx,ny,nz)
      real(kind=8) charge(nx,ny,nz)
      real(kind=8) J1(nx,ny,nz)
      real(kind=8) J2(nx,ny,nz)
      real(kind=8) J3(nx,ny,nz)
      real(kind=8) poyntingx_dens(nx,ny,nz)
      real(kind=8) poyntingy_dens(nx,ny,nz)
      real(kind=8) poyntingz_dens(nx,ny,nz)
      real(kind=8) uell(nx,ny,nz)
      real(kind=8) ham(nx,ny,nz)
      real(kind=8) momx(nx,ny,nz)
      real(kind=8) momy(nx,ny,nz)
      real(kind=8) momz(nx,ny,nz)
      real(kind=8) div_B(nx,ny,nz)
      real(kind=8) div_E(nx,ny,nz)
      real(kind=8) tr_A(nx,ny,nz)
      real(kind=8) detgt_m1(nx,ny,nz)
      real(kind=8) gamtx_con(nx,ny,nz)
      real(kind=8) gamty_con(nx,ny,nz)
      real(kind=8) gamtz_con(nx,ny,nz)
      real(kind=8) calgamtx_con(nx,ny,nz)
      real(kind=8) calgamty_con(nx,ny,nz)
      real(kind=8) calgamtz_con(nx,ny,nz)
      real(kind=8) cctk_x(nx,ny,nz)
      real(kind=8) cctk_y(nx,ny,nz)
      real(kind=8) cctk_z(nx,ny,nz)
      real(kind=8) r(nx,ny,nz)
      real(kind=8) xphys(nx,ny,nz)
      real(kind=8) yphys(nx,ny,nz)
      real(kind=8) zphys(nx,ny,nz)
      real(kind=8) mask(nx,ny,nz)
      real(kind=8) wdiss(nx,ny,nz)
      real(kind=8) chr(nx,ny,nz)
      real(kind=8) error(nx,ny,nz)
      real(kind=8) flag(nx,ny,nz)
      target :: gt11, gt11_np1
      target :: gt12, gt12_np1
      target :: gt13, gt13_np1
      target :: gt22, gt22_np1
      target :: gt23, gt23_np1
      target :: gt33, gt33_np1
      target :: A11, A11_np1
      target :: A12, A12_np1
      target :: A13, A13_np1
      target :: A22, A22_np1
      target :: A23, A23_np1
      target :: A33, A33_np1
      target :: chi, chi_np1
      target :: trK, trK_np1
      target :: Gam1, Gam1_np1
      target :: Gam2, Gam2_np1
      target :: Gam3, Gam3_np1
      target :: alpha, alpha_np1
      target :: shift1, shift1_np1
      target :: shift2, shift2_np1
      target :: shift3, shift3_np1
      target :: gb1, gb1_np1
      target :: gb2, gb2_np1
      target :: gb3, gb3_np1
      target :: Ex, Ex_np1
      target :: Ey, Ey_np1
      target :: Ez, Ez_np1
      target :: Bx, Bx_np1
      target :: By, By_np1
      target :: Bz, Bz_np1
      target :: Phi_em, Phi_em_np1
      target :: Psi_em, Psi_em_np1
      target :: phiR, phiR_np1
      target :: phiI, phiI_np1
      target :: piR, piR_np1
      target :: piI, piI_np1
      target :: ham
      target :: momx
      target :: momy
      target :: momz
      target :: div_B
      target :: div_E
      target :: tr_A
      target :: detgt_m1
      target :: gamtx_con
      target :: gamty_con
      target :: gamtz_con
      target :: calgamtx_con
      target :: calgamty_con
      target :: calgamtz_con
      target :: cctk_x
      target :: cctk_y
      target :: cctk_z
      target :: r
      target :: mask
      target :: wdiss
      target :: chr
      target :: error
      target :: flag
      target :: g11
      target :: g12
      target :: g13
      target :: g22
      target :: g23
      target :: g33
      target :: sdetg
      target :: rad_exp
      target :: psi4R
      target :: psi4I
      target :: massADM
      target :: massBONDI
      target :: curvature
      target :: phi2R
      target :: phi2I
      target :: phi0R
      target :: phi0I
      target :: vbr_omega
      target :: vbth_omega
      target :: charge
      target :: J1
      target :: J2
      target :: J3
      target :: poyntingx_dens
      target :: poyntingy_dens
      target :: poyntingz_dens
      target :: uell
      u0(H_GT11)%d => gt11
      u2(H_GT11)%d => gt11_np1
      u0(H_GT12)%d => gt12
      u2(H_GT12)%d => gt12_np1
      u0(H_GT13)%d => gt13
      u2(H_GT13)%d => gt13_np1
      u0(H_GT22)%d => gt22
      u2(H_GT22)%d => gt22_np1
      u0(H_GT23)%d => gt23
      u2(H_GT23)%d => gt23_np1
      u0(H_GT33)%d => gt33
      u2(H_GT33)%d => gt33_np1
      u0(H_A11)%d => A11
      u2(H_A11)%d => A11_np1
      u0(H_A12)%d => A12
      u2(H_A12)%d => A12_np1
      u0(H_A13)%d => A13
      u2(H_A13)%d => A13_np1
      u0(H_A22)%d => A22
      u2(H_A22)%d => A22_np1
      u0(H_A23)%d => A23
      u2(H_A23)%d => A23_np1
      u0(H_A33)%d => A33
      u2(H_A33)%d => A33_np1
      u0(H_CHI)%d => chi
      u2(H_CHI)%d => chi_np1
      u0(H_TRK)%d => trK
      u2(H_TRK)%d => trK_np1
      u0(H_GAM1)%d => Gam1
      u2(H_GAM1)%d => Gam1_np1
      u0(H_GAM2)%d => Gam2
      u2(H_GAM2)%d => Gam2_np1
      u0(H_GAM3)%d => Gam3
      u2(H_GAM3)%d => Gam3_np1
      u0(H_ALPHA)%d => alpha
      u2(H_ALPHA)%d => alpha_np1
      u0(H_SHIFT1)%d => shift1
      u2(H_SHIFT1)%d => shift1_np1
      u0(H_SHIFT2)%d => shift2
      u2(H_SHIFT2)%d => shift2_np1
      u0(H_SHIFT3)%d => shift3
      u2(H_SHIFT3)%d => shift3_np1
      u0(H_GB1)%d => gb1
      u2(H_GB1)%d => gb1_np1
      u0(H_GB2)%d => gb2
      u2(H_GB2)%d => gb2_np1
      u0(H_GB3)%d => gb3
      u2(H_GB3)%d => gb3_np1
      u0(H_EX)%d => Ex
      u2(H_EX)%d => Ex_np1
      u0(H_EY)%d => Ey
      u2(H_EY)%d => Ey_np1
      u0(H_EZ)%d => Ez
      u2(H_EZ)%d => Ez_np1
      u0(H_BX)%d => Bx
      u2(H_BX)%d => Bx_np1
      u0(H_BY)%d => By
      u2(H_BY)%d => By_np1
      u0(H_BZ)%d => Bz
      u2(H_BZ)%d => Bz_np1
      u0(H_PHI_EM)%d => Phi_em
      u2(H_PHI_EM)%d => Phi_em_np1
      u0(H_PSI_EM)%d => Psi_em
      u2(H_PSI_EM)%d => Psi_em_np1
      u0(H_PHIR)%d => phiR
      u2(H_PHIR)%d => phiR_np1
      u0(H_PHII)%d => phiI
      u2(H_PHII)%d => phiI_np1
      u0(H_PIR)%d => piR
      u2(H_PIR)%d => piR_np1
      u0(H_PII)%d => piI
      u2(H_PII)%d => piI_np1
      v(H_G11)%d => g11
      v(H_G12)%d => g12
      v(H_G13)%d => g13
      v(H_G22)%d => g22
      v(H_G23)%d => g23
      v(H_G33)%d => g33
      v(H_SDETG)%d => sdetg
      v(H_RAD_EXP)%d => rad_exp
      v(H_PSI4R)%d => psi4R
      v(H_PSI4I)%d => psi4I
      v(H_MASSADM)%d => massADM
      v(H_MASSBONDI)%d => massBONDI
      v(H_CURVATURE)%d => curvature
      v(H_PHI2R)%d => phi2R
      v(H_PHI2I)%d => phi2I
      v(H_PHI0R)%d => phi0R
      v(H_PHI0I)%d => phi0I
      v(H_VBR_OMEGA)%d => vbr_omega
      v(H_VBTH_OMEGA)%d => vbth_omega
      v(H_CHARGE)%d => charge
      v(H_J1)%d => J1
      v(H_J2)%d => J2
      v(H_J3)%d => J3
      v(H_POYNTINGX_DENS)%d => poyntingx_dens
      v(H_POYNTINGY_DENS)%d => poyntingy_dens
      v(H_POYNTINGZ_DENS)%d => poyntingz_dens
      v(H_UELL)%d => uell
      w(H_HAM)%d => ham
      w(H_MOMX)%d => momx
      w(H_MOMY)%d => momy
      w(H_MOMZ)%d => momz
      w(H_DIV_B)%d => div_B
      w(H_DIV_E)%d => div_E
      w(H_TR_A)%d => tr_A
      w(H_DETGT_M1)%d => detgt_m1
      w(H_GAMTX_CON)%d => gamtx_con
      w(H_GAMTY_CON)%d => gamty_con
      w(H_GAMTZ_CON)%d => gamtz_con
      w(H_CALGAMTX_CON)%d => calgamtx_con
      w(H_CALGAMTY_CON)%d => calgamty_con
      w(H_CALGAMTZ_CON)%d => calgamtz_con
      w(H_X)%d => cctk_x
      w(H_Y)%d => cctk_y
      w(H_Z)%d => cctk_z
      w(H_R)%d => r
      w(H_XPHYS)%d => cctk_x
      w(H_YPHYS)%d => cctk_y
      w(H_ZPHYS)%d => cctk_z
      w(H_MASK)%d => mask
      w(H_WDISS)%d => wdiss
      w(H_CHR)%d => chr
      w(H_ERROR)%d => error
      w(H_FLAG)%d => flag
      u0(H_GT11)%excise_val = 0
      u2(H_GT11)%excise_val = 0
      u0(H_GT11)%dissipation = 1
      u2(H_GT11)%dissipation = 1
      u0(H_GT11)%diss_factor = 1
      u2(H_GT11)%diss_factor = 1
      u0(H_GT11)%zsym = 1
      u0(H_GT11)%ysym = 1
      u0(H_GT11)%xsym = 1
      u2(H_GT11)%zsym = 1
      u2(H_GT11)%ysym = 1
      u2(H_GT11)%xsym = 1
      u0(H_GT11)%stiff = 0
      u2(H_GT11)%stiff = 0
      u2(H_GT11)%take_dx = 0
      u0(H_GT11)%take_dx = 0
      u2(H_GT11)%take_dy = 0
      u0(H_GT11)%take_dy = 0
      u2(H_GT11)%take_dz = 0
      u0(H_GT11)%take_dz = 0
      u0(H_GT12)%excise_val = 0
      u2(H_GT12)%excise_val = 0
      u0(H_GT12)%dissipation = 1
      u2(H_GT12)%dissipation = 1
      u0(H_GT12)%diss_factor = 1
      u2(H_GT12)%diss_factor = 1
      u0(H_GT12)%zsym = 1
      u0(H_GT12)%ysym = 1
      u0(H_GT12)%xsym = 1
      u2(H_GT12)%zsym = 1
      u2(H_GT12)%ysym = 1
      u2(H_GT12)%xsym = 1
      u0(H_GT12)%stiff = 0
      u2(H_GT12)%stiff = 0
      u2(H_GT12)%take_dx = 0
      u0(H_GT12)%take_dx = 0
      u2(H_GT12)%take_dy = 0
      u0(H_GT12)%take_dy = 0
      u2(H_GT12)%take_dz = 0
      u0(H_GT12)%take_dz = 0
      u0(H_GT13)%excise_val = 0
      u2(H_GT13)%excise_val = 0
      u0(H_GT13)%dissipation = 1
      u2(H_GT13)%dissipation = 1
      u0(H_GT13)%diss_factor = 1
      u2(H_GT13)%diss_factor = 1
      u0(H_GT13)%zsym = 1
      u0(H_GT13)%ysym = 1
      u0(H_GT13)%xsym = 1
      u2(H_GT13)%zsym = 1
      u2(H_GT13)%ysym = 1
      u2(H_GT13)%xsym = 1
      u0(H_GT13)%stiff = 0
      u2(H_GT13)%stiff = 0
      u2(H_GT13)%take_dx = 0
      u0(H_GT13)%take_dx = 0
      u2(H_GT13)%take_dy = 0
      u0(H_GT13)%take_dy = 0
      u2(H_GT13)%take_dz = 0
      u0(H_GT13)%take_dz = 0
      u0(H_GT22)%excise_val = 0
      u2(H_GT22)%excise_val = 0
      u0(H_GT22)%dissipation = 1
      u2(H_GT22)%dissipation = 1
      u0(H_GT22)%diss_factor = 1
      u2(H_GT22)%diss_factor = 1
      u0(H_GT22)%zsym = 1
      u0(H_GT22)%ysym = 1
      u0(H_GT22)%xsym = 1
      u2(H_GT22)%zsym = 1
      u2(H_GT22)%ysym = 1
      u2(H_GT22)%xsym = 1
      u0(H_GT22)%stiff = 0
      u2(H_GT22)%stiff = 0
      u2(H_GT22)%take_dx = 0
      u0(H_GT22)%take_dx = 0
      u2(H_GT22)%take_dy = 0
      u0(H_GT22)%take_dy = 0
      u2(H_GT22)%take_dz = 0
      u0(H_GT22)%take_dz = 0
      u0(H_GT23)%excise_val = 0
      u2(H_GT23)%excise_val = 0
      u0(H_GT23)%dissipation = 1
      u2(H_GT23)%dissipation = 1
      u0(H_GT23)%diss_factor = 1
      u2(H_GT23)%diss_factor = 1
      u0(H_GT23)%zsym = 1
      u0(H_GT23)%ysym = 1
      u0(H_GT23)%xsym = 1
      u2(H_GT23)%zsym = 1
      u2(H_GT23)%ysym = 1
      u2(H_GT23)%xsym = 1
      u0(H_GT23)%stiff = 0
      u2(H_GT23)%stiff = 0
      u2(H_GT23)%take_dx = 0
      u0(H_GT23)%take_dx = 0
      u2(H_GT23)%take_dy = 0
      u0(H_GT23)%take_dy = 0
      u2(H_GT23)%take_dz = 0
      u0(H_GT23)%take_dz = 0
      u0(H_GT33)%excise_val = 0
      u2(H_GT33)%excise_val = 0
      u0(H_GT33)%dissipation = 1
      u2(H_GT33)%dissipation = 1
      u0(H_GT33)%diss_factor = 1
      u2(H_GT33)%diss_factor = 1
      u0(H_GT33)%zsym = 1
      u0(H_GT33)%ysym = 1
      u0(H_GT33)%xsym = 1
      u2(H_GT33)%zsym = 1
      u2(H_GT33)%ysym = 1
      u2(H_GT33)%xsym = 1
      u0(H_GT33)%stiff = 0
      u2(H_GT33)%stiff = 0
      u2(H_GT33)%take_dx = 0
      u0(H_GT33)%take_dx = 0
      u2(H_GT33)%take_dy = 0
      u0(H_GT33)%take_dy = 0
      u2(H_GT33)%take_dz = 0
      u0(H_GT33)%take_dz = 0
      u0(H_A11)%excise_val = 0
      u2(H_A11)%excise_val = 0
      u0(H_A11)%dissipation = 1
      u2(H_A11)%dissipation = 1
      u0(H_A11)%diss_factor = 1
      u2(H_A11)%diss_factor = 1
      u0(H_A11)%zsym = 1
      u0(H_A11)%ysym = 1
      u0(H_A11)%xsym = 1
      u2(H_A11)%zsym = 1
      u2(H_A11)%ysym = 1
      u2(H_A11)%xsym = 1
      u0(H_A11)%stiff = 0
      u2(H_A11)%stiff = 0
      u2(H_A11)%take_dx = 0
      u0(H_A11)%take_dx = 0
      u2(H_A11)%take_dy = 0
      u0(H_A11)%take_dy = 0
      u2(H_A11)%take_dz = 0
      u0(H_A11)%take_dz = 0
      u0(H_A12)%excise_val = 0
      u2(H_A12)%excise_val = 0
      u0(H_A12)%dissipation = 1
      u2(H_A12)%dissipation = 1
      u0(H_A12)%diss_factor = 1
      u2(H_A12)%diss_factor = 1
      u0(H_A12)%zsym = 1
      u0(H_A12)%ysym = 1
      u0(H_A12)%xsym = 1
      u2(H_A12)%zsym = 1
      u2(H_A12)%ysym = 1
      u2(H_A12)%xsym = 1
      u0(H_A12)%stiff = 0
      u2(H_A12)%stiff = 0
      u2(H_A12)%take_dx = 0
      u0(H_A12)%take_dx = 0
      u2(H_A12)%take_dy = 0
      u0(H_A12)%take_dy = 0
      u2(H_A12)%take_dz = 0
      u0(H_A12)%take_dz = 0
      u0(H_A13)%excise_val = 0
      u2(H_A13)%excise_val = 0
      u0(H_A13)%dissipation = 1
      u2(H_A13)%dissipation = 1
      u0(H_A13)%diss_factor = 1
      u2(H_A13)%diss_factor = 1
      u0(H_A13)%zsym = 1
      u0(H_A13)%ysym = 1
      u0(H_A13)%xsym = 1
      u2(H_A13)%zsym = 1
      u2(H_A13)%ysym = 1
      u2(H_A13)%xsym = 1
      u0(H_A13)%stiff = 0
      u2(H_A13)%stiff = 0
      u2(H_A13)%take_dx = 0
      u0(H_A13)%take_dx = 0
      u2(H_A13)%take_dy = 0
      u0(H_A13)%take_dy = 0
      u2(H_A13)%take_dz = 0
      u0(H_A13)%take_dz = 0
      u0(H_A22)%excise_val = 0
      u2(H_A22)%excise_val = 0
      u0(H_A22)%dissipation = 1
      u2(H_A22)%dissipation = 1
      u0(H_A22)%diss_factor = 1
      u2(H_A22)%diss_factor = 1
      u0(H_A22)%zsym = 1
      u0(H_A22)%ysym = 1
      u0(H_A22)%xsym = 1
      u2(H_A22)%zsym = 1
      u2(H_A22)%ysym = 1
      u2(H_A22)%xsym = 1
      u0(H_A22)%stiff = 0
      u2(H_A22)%stiff = 0
      u2(H_A22)%take_dx = 0
      u0(H_A22)%take_dx = 0
      u2(H_A22)%take_dy = 0
      u0(H_A22)%take_dy = 0
      u2(H_A22)%take_dz = 0
      u0(H_A22)%take_dz = 0
      u0(H_A23)%excise_val = 0
      u2(H_A23)%excise_val = 0
      u0(H_A23)%dissipation = 1
      u2(H_A23)%dissipation = 1
      u0(H_A23)%diss_factor = 1
      u2(H_A23)%diss_factor = 1
      u0(H_A23)%zsym = 1
      u0(H_A23)%ysym = 1
      u0(H_A23)%xsym = 1
      u2(H_A23)%zsym = 1
      u2(H_A23)%ysym = 1
      u2(H_A23)%xsym = 1
      u0(H_A23)%stiff = 0
      u2(H_A23)%stiff = 0
      u2(H_A23)%take_dx = 0
      u0(H_A23)%take_dx = 0
      u2(H_A23)%take_dy = 0
      u0(H_A23)%take_dy = 0
      u2(H_A23)%take_dz = 0
      u0(H_A23)%take_dz = 0
      u0(H_A33)%excise_val = 0
      u2(H_A33)%excise_val = 0
      u0(H_A33)%dissipation = 1
      u2(H_A33)%dissipation = 1
      u0(H_A33)%diss_factor = 1
      u2(H_A33)%diss_factor = 1
      u0(H_A33)%zsym = 1
      u0(H_A33)%ysym = 1
      u0(H_A33)%xsym = 1
      u2(H_A33)%zsym = 1
      u2(H_A33)%ysym = 1
      u2(H_A33)%xsym = 1
      u0(H_A33)%stiff = 0
      u2(H_A33)%stiff = 0
      u2(H_A33)%take_dx = 0
      u0(H_A33)%take_dx = 0
      u2(H_A33)%take_dy = 0
      u0(H_A33)%take_dy = 0
      u2(H_A33)%take_dz = 0
      u0(H_A33)%take_dz = 0
      u0(H_CHI)%excise_val = 0
      u2(H_CHI)%excise_val = 0
      u0(H_CHI)%dissipation = 1
      u2(H_CHI)%dissipation = 1
      u0(H_CHI)%diss_factor = 1
      u2(H_CHI)%diss_factor = 1
      u0(H_CHI)%zsym = 1
      u0(H_CHI)%ysym = 1
      u0(H_CHI)%xsym = 1
      u2(H_CHI)%zsym = 1
      u2(H_CHI)%ysym = 1
      u2(H_CHI)%xsym = 1
      u0(H_CHI)%stiff = 0
      u2(H_CHI)%stiff = 0
      u2(H_CHI)%take_dx = 0
      u0(H_CHI)%take_dx = 0
      u2(H_CHI)%take_dy = 0
      u0(H_CHI)%take_dy = 0
      u2(H_CHI)%take_dz = 0
      u0(H_CHI)%take_dz = 0
      u0(H_TRK)%excise_val = 0
      u2(H_TRK)%excise_val = 0
      u0(H_TRK)%dissipation = 1
      u2(H_TRK)%dissipation = 1
      u0(H_TRK)%diss_factor = 1
      u2(H_TRK)%diss_factor = 1
      u0(H_TRK)%zsym = 1
      u0(H_TRK)%ysym = 1
      u0(H_TRK)%xsym = 1
      u2(H_TRK)%zsym = 1
      u2(H_TRK)%ysym = 1
      u2(H_TRK)%xsym = 1
      u0(H_TRK)%stiff = 0
      u2(H_TRK)%stiff = 0
      u2(H_TRK)%take_dx = 0
      u0(H_TRK)%take_dx = 0
      u2(H_TRK)%take_dy = 0
      u0(H_TRK)%take_dy = 0
      u2(H_TRK)%take_dz = 0
      u0(H_TRK)%take_dz = 0
      u0(H_GAM1)%excise_val = 0
      u2(H_GAM1)%excise_val = 0
      u0(H_GAM1)%dissipation = 1
      u2(H_GAM1)%dissipation = 1
      u0(H_GAM1)%diss_factor = 1
      u2(H_GAM1)%diss_factor = 1
      u0(H_GAM1)%zsym = 1
      u0(H_GAM1)%ysym = 1
      u0(H_GAM1)%xsym = 1
      u2(H_GAM1)%zsym = 1
      u2(H_GAM1)%ysym = 1
      u2(H_GAM1)%xsym = 1
      u0(H_GAM1)%stiff = 0
      u2(H_GAM1)%stiff = 0
      u2(H_GAM1)%take_dx = 0
      u0(H_GAM1)%take_dx = 0
      u2(H_GAM1)%take_dy = 0
      u0(H_GAM1)%take_dy = 0
      u2(H_GAM1)%take_dz = 0
      u0(H_GAM1)%take_dz = 0
      u0(H_GAM2)%excise_val = 0
      u2(H_GAM2)%excise_val = 0
      u0(H_GAM2)%dissipation = 1
      u2(H_GAM2)%dissipation = 1
      u0(H_GAM2)%diss_factor = 1
      u2(H_GAM2)%diss_factor = 1
      u0(H_GAM2)%zsym = 1
      u0(H_GAM2)%ysym = 1
      u0(H_GAM2)%xsym = 1
      u2(H_GAM2)%zsym = 1
      u2(H_GAM2)%ysym = 1
      u2(H_GAM2)%xsym = 1
      u0(H_GAM2)%stiff = 0
      u2(H_GAM2)%stiff = 0
      u2(H_GAM2)%take_dx = 0
      u0(H_GAM2)%take_dx = 0
      u2(H_GAM2)%take_dy = 0
      u0(H_GAM2)%take_dy = 0
      u2(H_GAM2)%take_dz = 0
      u0(H_GAM2)%take_dz = 0
      u0(H_GAM3)%excise_val = 0
      u2(H_GAM3)%excise_val = 0
      u0(H_GAM3)%dissipation = 1
      u2(H_GAM3)%dissipation = 1
      u0(H_GAM3)%diss_factor = 1
      u2(H_GAM3)%diss_factor = 1
      u0(H_GAM3)%zsym = 1
      u0(H_GAM3)%ysym = 1
      u0(H_GAM3)%xsym = 1
      u2(H_GAM3)%zsym = 1
      u2(H_GAM3)%ysym = 1
      u2(H_GAM3)%xsym = 1
      u0(H_GAM3)%stiff = 0
      u2(H_GAM3)%stiff = 0
      u2(H_GAM3)%take_dx = 0
      u0(H_GAM3)%take_dx = 0
      u2(H_GAM3)%take_dy = 0
      u0(H_GAM3)%take_dy = 0
      u2(H_GAM3)%take_dz = 0
      u0(H_GAM3)%take_dz = 0
      u0(H_ALPHA)%excise_val = 0
      u2(H_ALPHA)%excise_val = 0
      u0(H_ALPHA)%dissipation = 1
      u2(H_ALPHA)%dissipation = 1
      u0(H_ALPHA)%diss_factor = 1
      u2(H_ALPHA)%diss_factor = 1
      u0(H_ALPHA)%zsym = 1
      u0(H_ALPHA)%ysym = 1
      u0(H_ALPHA)%xsym = 1
      u2(H_ALPHA)%zsym = 1
      u2(H_ALPHA)%ysym = 1
      u2(H_ALPHA)%xsym = 1
      u0(H_ALPHA)%stiff = 0
      u2(H_ALPHA)%stiff = 0
      u2(H_ALPHA)%take_dx = 0
      u0(H_ALPHA)%take_dx = 0
      u2(H_ALPHA)%take_dy = 0
      u0(H_ALPHA)%take_dy = 0
      u2(H_ALPHA)%take_dz = 0
      u0(H_ALPHA)%take_dz = 0
      u0(H_SHIFT1)%excise_val = 0
      u2(H_SHIFT1)%excise_val = 0
      u0(H_SHIFT1)%dissipation = 1
      u2(H_SHIFT1)%dissipation = 1
      u0(H_SHIFT1)%diss_factor = 1
      u2(H_SHIFT1)%diss_factor = 1
      u0(H_SHIFT1)%zsym = 1
      u0(H_SHIFT1)%ysym = 1
      u0(H_SHIFT1)%xsym = 1
      u2(H_SHIFT1)%zsym = 1
      u2(H_SHIFT1)%ysym = 1
      u2(H_SHIFT1)%xsym = 1
      u0(H_SHIFT1)%stiff = 0
      u2(H_SHIFT1)%stiff = 0
      u2(H_SHIFT1)%take_dx = 0
      u0(H_SHIFT1)%take_dx = 0
      u2(H_SHIFT1)%take_dy = 0
      u0(H_SHIFT1)%take_dy = 0
      u2(H_SHIFT1)%take_dz = 0
      u0(H_SHIFT1)%take_dz = 0
      u0(H_SHIFT2)%excise_val = 0
      u2(H_SHIFT2)%excise_val = 0
      u0(H_SHIFT2)%dissipation = 1
      u2(H_SHIFT2)%dissipation = 1
      u0(H_SHIFT2)%diss_factor = 1
      u2(H_SHIFT2)%diss_factor = 1
      u0(H_SHIFT2)%zsym = 1
      u0(H_SHIFT2)%ysym = 1
      u0(H_SHIFT2)%xsym = 1
      u2(H_SHIFT2)%zsym = 1
      u2(H_SHIFT2)%ysym = 1
      u2(H_SHIFT2)%xsym = 1
      u0(H_SHIFT2)%stiff = 0
      u2(H_SHIFT2)%stiff = 0
      u2(H_SHIFT2)%take_dx = 0
      u0(H_SHIFT2)%take_dx = 0
      u2(H_SHIFT2)%take_dy = 0
      u0(H_SHIFT2)%take_dy = 0
      u2(H_SHIFT2)%take_dz = 0
      u0(H_SHIFT2)%take_dz = 0
      u0(H_SHIFT3)%excise_val = 0
      u2(H_SHIFT3)%excise_val = 0
      u0(H_SHIFT3)%dissipation = 1
      u2(H_SHIFT3)%dissipation = 1
      u0(H_SHIFT3)%diss_factor = 1
      u2(H_SHIFT3)%diss_factor = 1
      u0(H_SHIFT3)%zsym = 1
      u0(H_SHIFT3)%ysym = 1
      u0(H_SHIFT3)%xsym = 1
      u2(H_SHIFT3)%zsym = 1
      u2(H_SHIFT3)%ysym = 1
      u2(H_SHIFT3)%xsym = 1
      u0(H_SHIFT3)%stiff = 0
      u2(H_SHIFT3)%stiff = 0
      u2(H_SHIFT3)%take_dx = 0
      u0(H_SHIFT3)%take_dx = 0
      u2(H_SHIFT3)%take_dy = 0
      u0(H_SHIFT3)%take_dy = 0
      u2(H_SHIFT3)%take_dz = 0
      u0(H_SHIFT3)%take_dz = 0
      u0(H_GB1)%excise_val = 0
      u2(H_GB1)%excise_val = 0
      u0(H_GB1)%dissipation = 1
      u2(H_GB1)%dissipation = 1
      u0(H_GB1)%diss_factor = 1
      u2(H_GB1)%diss_factor = 1
      u0(H_GB1)%zsym = 1
      u0(H_GB1)%ysym = 1
      u0(H_GB1)%xsym = 1
      u2(H_GB1)%zsym = 1
      u2(H_GB1)%ysym = 1
      u2(H_GB1)%xsym = 1
      u0(H_GB1)%stiff = 0
      u2(H_GB1)%stiff = 0
      u2(H_GB1)%take_dx = 0
      u0(H_GB1)%take_dx = 0
      u2(H_GB1)%take_dy = 0
      u0(H_GB1)%take_dy = 0
      u2(H_GB1)%take_dz = 0
      u0(H_GB1)%take_dz = 0
      u0(H_GB2)%excise_val = 0
      u2(H_GB2)%excise_val = 0
      u0(H_GB2)%dissipation = 1
      u2(H_GB2)%dissipation = 1
      u0(H_GB2)%diss_factor = 1
      u2(H_GB2)%diss_factor = 1
      u0(H_GB2)%zsym = 1
      u0(H_GB2)%ysym = 1
      u0(H_GB2)%xsym = 1
      u2(H_GB2)%zsym = 1
      u2(H_GB2)%ysym = 1
      u2(H_GB2)%xsym = 1
      u0(H_GB2)%stiff = 0
      u2(H_GB2)%stiff = 0
      u2(H_GB2)%take_dx = 0
      u0(H_GB2)%take_dx = 0
      u2(H_GB2)%take_dy = 0
      u0(H_GB2)%take_dy = 0
      u2(H_GB2)%take_dz = 0
      u0(H_GB2)%take_dz = 0
      u0(H_GB3)%excise_val = 0
      u2(H_GB3)%excise_val = 0
      u0(H_GB3)%dissipation = 1
      u2(H_GB3)%dissipation = 1
      u0(H_GB3)%diss_factor = 1
      u2(H_GB3)%diss_factor = 1
      u0(H_GB3)%zsym = 1
      u0(H_GB3)%ysym = 1
      u0(H_GB3)%xsym = 1
      u2(H_GB3)%zsym = 1
      u2(H_GB3)%ysym = 1
      u2(H_GB3)%xsym = 1
      u0(H_GB3)%stiff = 0
      u2(H_GB3)%stiff = 0
      u2(H_GB3)%take_dx = 0
      u0(H_GB3)%take_dx = 0
      u2(H_GB3)%take_dy = 0
      u0(H_GB3)%take_dy = 0
      u2(H_GB3)%take_dz = 0
      u0(H_GB3)%take_dz = 0
      u0(H_EX)%excise_val = 0
      u2(H_EX)%excise_val = 0
      u0(H_EX)%dissipation = 1
      u2(H_EX)%dissipation = 1
      u0(H_EX)%diss_factor = 1
      u2(H_EX)%diss_factor = 1
      u0(H_EX)%zsym = 1
      u0(H_EX)%ysym = 1
      u0(H_EX)%xsym = 1
      u2(H_EX)%zsym = 1
      u2(H_EX)%ysym = 1
      u2(H_EX)%xsym = 1
      u0(H_EX)%stiff = 0
      u2(H_EX)%stiff = 0
      u2(H_EX)%take_dx = 0
      u0(H_EX)%take_dx = 0
      u2(H_EX)%take_dy = 0
      u0(H_EX)%take_dy = 0
      u2(H_EX)%take_dz = 0
      u0(H_EX)%take_dz = 0
      u0(H_EY)%excise_val = 0
      u2(H_EY)%excise_val = 0
      u0(H_EY)%dissipation = 1
      u2(H_EY)%dissipation = 1
      u0(H_EY)%diss_factor = 1
      u2(H_EY)%diss_factor = 1
      u0(H_EY)%zsym = 1
      u0(H_EY)%ysym = 1
      u0(H_EY)%xsym = 1
      u2(H_EY)%zsym = 1
      u2(H_EY)%ysym = 1
      u2(H_EY)%xsym = 1
      u0(H_EY)%stiff = 0
      u2(H_EY)%stiff = 0
      u2(H_EY)%take_dx = 0
      u0(H_EY)%take_dx = 0
      u2(H_EY)%take_dy = 0
      u0(H_EY)%take_dy = 0
      u2(H_EY)%take_dz = 0
      u0(H_EY)%take_dz = 0
      u0(H_EZ)%excise_val = 0
      u2(H_EZ)%excise_val = 0
      u0(H_EZ)%dissipation = 1
      u2(H_EZ)%dissipation = 1
      u0(H_EZ)%diss_factor = 1
      u2(H_EZ)%diss_factor = 1
      u0(H_EZ)%zsym = 1
      u0(H_EZ)%ysym = 1
      u0(H_EZ)%xsym = 1
      u2(H_EZ)%zsym = 1
      u2(H_EZ)%ysym = 1
      u2(H_EZ)%xsym = 1
      u0(H_EZ)%stiff = 0
      u2(H_EZ)%stiff = 0
      u2(H_EZ)%take_dx = 0
      u0(H_EZ)%take_dx = 0
      u2(H_EZ)%take_dy = 0
      u0(H_EZ)%take_dy = 0
      u2(H_EZ)%take_dz = 0
      u0(H_EZ)%take_dz = 0
      u0(H_BX)%excise_val = 0
      u2(H_BX)%excise_val = 0
      u0(H_BX)%dissipation = 1
      u2(H_BX)%dissipation = 1
      u0(H_BX)%diss_factor = 1
      u2(H_BX)%diss_factor = 1
      u0(H_BX)%zsym = 1
      u0(H_BX)%ysym = 1
      u0(H_BX)%xsym = 1
      u2(H_BX)%zsym = 1
      u2(H_BX)%ysym = 1
      u2(H_BX)%xsym = 1
      u0(H_BX)%stiff = 0
      u2(H_BX)%stiff = 0
      u2(H_BX)%take_dx = 0
      u0(H_BX)%take_dx = 0
      u2(H_BX)%take_dy = 0
      u0(H_BX)%take_dy = 0
      u2(H_BX)%take_dz = 0
      u0(H_BX)%take_dz = 0
      u0(H_BY)%excise_val = 0
      u2(H_BY)%excise_val = 0
      u0(H_BY)%dissipation = 1
      u2(H_BY)%dissipation = 1
      u0(H_BY)%diss_factor = 1
      u2(H_BY)%diss_factor = 1
      u0(H_BY)%zsym = 1
      u0(H_BY)%ysym = 1
      u0(H_BY)%xsym = 1
      u2(H_BY)%zsym = 1
      u2(H_BY)%ysym = 1
      u2(H_BY)%xsym = 1
      u0(H_BY)%stiff = 0
      u2(H_BY)%stiff = 0
      u2(H_BY)%take_dx = 0
      u0(H_BY)%take_dx = 0
      u2(H_BY)%take_dy = 0
      u0(H_BY)%take_dy = 0
      u2(H_BY)%take_dz = 0
      u0(H_BY)%take_dz = 0
      u0(H_BZ)%excise_val = 0
      u2(H_BZ)%excise_val = 0
      u0(H_BZ)%dissipation = 1
      u2(H_BZ)%dissipation = 1
      u0(H_BZ)%diss_factor = 1
      u2(H_BZ)%diss_factor = 1
      u0(H_BZ)%zsym = 1
      u0(H_BZ)%ysym = 1
      u0(H_BZ)%xsym = 1
      u2(H_BZ)%zsym = 1
      u2(H_BZ)%ysym = 1
      u2(H_BZ)%xsym = 1
      u0(H_BZ)%stiff = 0
      u2(H_BZ)%stiff = 0
      u2(H_BZ)%take_dx = 0
      u0(H_BZ)%take_dx = 0
      u2(H_BZ)%take_dy = 0
      u0(H_BZ)%take_dy = 0
      u2(H_BZ)%take_dz = 0
      u0(H_BZ)%take_dz = 0
      u0(H_PHI_EM)%excise_val = 0
      u2(H_PHI_EM)%excise_val = 0
      u0(H_PHI_EM)%dissipation = 1
      u2(H_PHI_EM)%dissipation = 1
      u0(H_PHI_EM)%diss_factor = 1
      u2(H_PHI_EM)%diss_factor = 1
      u0(H_PHI_EM)%zsym = 1
      u0(H_PHI_EM)%ysym = 1
      u0(H_PHI_EM)%xsym = 1
      u2(H_PHI_EM)%zsym = 1
      u2(H_PHI_EM)%ysym = 1
      u2(H_PHI_EM)%xsym = 1
      u0(H_PHI_EM)%stiff = 0
      u2(H_PHI_EM)%stiff = 0
      u2(H_PHI_EM)%take_dx = 0
      u0(H_PHI_EM)%take_dx = 0
      u2(H_PHI_EM)%take_dy = 0
      u0(H_PHI_EM)%take_dy = 0
      u2(H_PHI_EM)%take_dz = 0
      u0(H_PHI_EM)%take_dz = 0
      u0(H_PSI_EM)%excise_val = 0
      u2(H_PSI_EM)%excise_val = 0
      u0(H_PSI_EM)%dissipation = 1
      u2(H_PSI_EM)%dissipation = 1
      u0(H_PSI_EM)%diss_factor = 1
      u2(H_PSI_EM)%diss_factor = 1
      u0(H_PSI_EM)%zsym = 1
      u0(H_PSI_EM)%ysym = 1
      u0(H_PSI_EM)%xsym = 1
      u2(H_PSI_EM)%zsym = 1
      u2(H_PSI_EM)%ysym = 1
      u2(H_PSI_EM)%xsym = 1
      u0(H_PSI_EM)%stiff = 0
      u2(H_PSI_EM)%stiff = 0
      u2(H_PSI_EM)%take_dx = 0
      u0(H_PSI_EM)%take_dx = 0
      u2(H_PSI_EM)%take_dy = 0
      u0(H_PSI_EM)%take_dy = 0
      u2(H_PSI_EM)%take_dz = 0
      u0(H_PSI_EM)%take_dz = 0
      u0(H_PHIR)%excise_val = 0
      u2(H_PHIR)%excise_val = 0
      u0(H_PHIR)%dissipation = 1
      u2(H_PHIR)%dissipation = 1
      u0(H_PHIR)%diss_factor = 1
      u2(H_PHIR)%diss_factor = 1
      u0(H_PHIR)%zsym = 1
      u0(H_PHIR)%ysym = 1
      u0(H_PHIR)%xsym = 1
      u2(H_PHIR)%zsym = 1
      u2(H_PHIR)%ysym = 1
      u2(H_PHIR)%xsym = 1
      u0(H_PHIR)%stiff = 0
      u2(H_PHIR)%stiff = 0
      u2(H_PHIR)%take_dx = 0
      u0(H_PHIR)%take_dx = 0
      u2(H_PHIR)%take_dy = 0
      u0(H_PHIR)%take_dy = 0
      u2(H_PHIR)%take_dz = 0
      u0(H_PHIR)%take_dz = 0
      u0(H_PHII)%excise_val = 0
      u2(H_PHII)%excise_val = 0
      u0(H_PHII)%dissipation = 1
      u2(H_PHII)%dissipation = 1
      u0(H_PHII)%diss_factor = 1
      u2(H_PHII)%diss_factor = 1
      u0(H_PHII)%zsym = 1
      u0(H_PHII)%ysym = 1
      u0(H_PHII)%xsym = 1
      u2(H_PHII)%zsym = 1
      u2(H_PHII)%ysym = 1
      u2(H_PHII)%xsym = 1
      u0(H_PHII)%stiff = 0
      u2(H_PHII)%stiff = 0
      u2(H_PHII)%take_dx = 0
      u0(H_PHII)%take_dx = 0
      u2(H_PHII)%take_dy = 0
      u0(H_PHII)%take_dy = 0
      u2(H_PHII)%take_dz = 0
      u0(H_PHII)%take_dz = 0
      u0(H_PIR)%excise_val = 0
      u2(H_PIR)%excise_val = 0
      u0(H_PIR)%dissipation = 1
      u2(H_PIR)%dissipation = 1
      u0(H_PIR)%diss_factor = 1
      u2(H_PIR)%diss_factor = 1
      u0(H_PIR)%zsym = 1
      u0(H_PIR)%ysym = 1
      u0(H_PIR)%xsym = 1
      u2(H_PIR)%zsym = 1
      u2(H_PIR)%ysym = 1
      u2(H_PIR)%xsym = 1
      u0(H_PIR)%stiff = 0
      u2(H_PIR)%stiff = 0
      u2(H_PIR)%take_dx = 0
      u0(H_PIR)%take_dx = 0
      u2(H_PIR)%take_dy = 0
      u0(H_PIR)%take_dy = 0
      u2(H_PIR)%take_dz = 0
      u0(H_PIR)%take_dz = 0
      u0(H_PII)%excise_val = 0
      u2(H_PII)%excise_val = 0
      u0(H_PII)%dissipation = 1
      u2(H_PII)%dissipation = 1
      u0(H_PII)%diss_factor = 1
      u2(H_PII)%diss_factor = 1
      u0(H_PII)%zsym = 1
      u0(H_PII)%ysym = 1
      u0(H_PII)%xsym = 1
      u2(H_PII)%zsym = 1
      u2(H_PII)%ysym = 1
      u2(H_PII)%xsym = 1
      u0(H_PII)%stiff = 0
      u2(H_PII)%stiff = 0
      u2(H_PII)%take_dx = 0
      u0(H_PII)%take_dx = 0
      u2(H_PII)%take_dy = 0
      u0(H_PII)%take_dy = 0
      u2(H_PII)%take_dz = 0
      u0(H_PII)%take_dz = 0
      v(H_G11)%take_dx = 0
      v(H_G11)%take_dy = 0
      v(H_G11)%take_dz = 0
      v(H_G12)%take_dx = 0
      v(H_G12)%take_dy = 0
      v(H_G12)%take_dz = 0
      v(H_G13)%take_dx = 0
      v(H_G13)%take_dy = 0
      v(H_G13)%take_dz = 0
      v(H_G22)%take_dx = 0
      v(H_G22)%take_dy = 0
      v(H_G22)%take_dz = 0
      v(H_G23)%take_dx = 0
      v(H_G23)%take_dy = 0
      v(H_G23)%take_dz = 0
      v(H_G33)%take_dx = 0
      v(H_G33)%take_dy = 0
      v(H_G33)%take_dz = 0
      v(H_SDETG)%take_dx = 0
      v(H_SDETG)%take_dy = 0
      v(H_SDETG)%take_dz = 0
      v(H_RAD_EXP)%take_dx = 0
      v(H_RAD_EXP)%take_dy = 0
      v(H_RAD_EXP)%take_dz = 0
      v(H_PSI4R)%take_dx = 0
      v(H_PSI4R)%take_dy = 0
      v(H_PSI4R)%take_dz = 0
      v(H_PSI4I)%take_dx = 0
      v(H_PSI4I)%take_dy = 0
      v(H_PSI4I)%take_dz = 0
      v(H_MASSADM)%take_dx = 0
      v(H_MASSADM)%take_dy = 0
      v(H_MASSADM)%take_dz = 0
      v(H_MASSBONDI)%take_dx = 0
      v(H_MASSBONDI)%take_dy = 0
      v(H_MASSBONDI)%take_dz = 0
      v(H_CURVATURE)%take_dx = 0
      v(H_CURVATURE)%take_dy = 0
      v(H_CURVATURE)%take_dz = 0
      v(H_PHI2R)%take_dx = 0
      v(H_PHI2R)%take_dy = 0
      v(H_PHI2R)%take_dz = 0
      v(H_PHI2I)%take_dx = 0
      v(H_PHI2I)%take_dy = 0
      v(H_PHI2I)%take_dz = 0
      v(H_PHI0R)%take_dx = 0
      v(H_PHI0R)%take_dy = 0
      v(H_PHI0R)%take_dz = 0
      v(H_PHI0I)%take_dx = 0
      v(H_PHI0I)%take_dy = 0
      v(H_PHI0I)%take_dz = 0
      v(H_VBR_OMEGA)%take_dx = 0
      v(H_VBR_OMEGA)%take_dy = 0
      v(H_VBR_OMEGA)%take_dz = 0
      v(H_VBTH_OMEGA)%take_dx = 0
      v(H_VBTH_OMEGA)%take_dy = 0
      v(H_VBTH_OMEGA)%take_dz = 0
      v(H_CHARGE)%take_dx = 0
      v(H_CHARGE)%take_dy = 0
      v(H_CHARGE)%take_dz = 0
      v(H_J1)%take_dx = 0
      v(H_J1)%take_dy = 0
      v(H_J1)%take_dz = 0
      v(H_J2)%take_dx = 0
      v(H_J2)%take_dy = 0
      v(H_J2)%take_dz = 0
      v(H_J3)%take_dx = 0
      v(H_J3)%take_dy = 0
      v(H_J3)%take_dz = 0
      v(H_POYNTINGX_DENS)%take_dx = 0
      v(H_POYNTINGX_DENS)%take_dy = 0
      v(H_POYNTINGX_DENS)%take_dz = 0
      v(H_POYNTINGY_DENS)%take_dx = 0
      v(H_POYNTINGY_DENS)%take_dy = 0
      v(H_POYNTINGY_DENS)%take_dz = 0
      v(H_POYNTINGZ_DENS)%take_dx = 0
      v(H_POYNTINGZ_DENS)%take_dy = 0
      v(H_POYNTINGZ_DENS)%take_dz = 0
      v(H_UELL)%take_dx = 0
      v(H_UELL)%take_dy = 0
      v(H_UELL)%take_dz = 0
      return
      end subroutine assign_ptrs_fields
      subroutine assign_ptrs_rks(gt11_rk1,gt12_rk1,gt13_rk1,gt22_rk1,gt2
     &3_rk1,gt33_rk1,A11_rk1,A12_rk1,A13_rk1,A22_rk1,A23_rk1,A33_rk1,chi
     &_rk1,trK_rk1,Gam1_rk1,Gam2_rk1,Gam3_rk1,alpha_rk1,shift1_rk1,shift
     &2_rk1,shift3_rk1,gb1_rk1,gb2_rk1,gb3_rk1,Ex_rk1,Ey_rk1,Ez_rk1,Bx_r
     &k1,By_rk1,Bz_rk1,Phi_em_rk1,Psi_em_rk1,phiR_rk1,phiI_rk1,piR_rk1,p
     &iI_rk1,urk1,nx,ny,nz)
      use GF
      implicit none
      integer nx, ny, nz
      type(gridfunction), dimension(NU) :: urk1
      real(kind=8) gt11_rk1(nx,ny,nz)
      real(kind=8) gt12_rk1(nx,ny,nz)
      real(kind=8) gt13_rk1(nx,ny,nz)
      real(kind=8) gt22_rk1(nx,ny,nz)
      real(kind=8) gt23_rk1(nx,ny,nz)
      real(kind=8) gt33_rk1(nx,ny,nz)
      real(kind=8) A11_rk1(nx,ny,nz)
      real(kind=8) A12_rk1(nx,ny,nz)
      real(kind=8) A13_rk1(nx,ny,nz)
      real(kind=8) A22_rk1(nx,ny,nz)
      real(kind=8) A23_rk1(nx,ny,nz)
      real(kind=8) A33_rk1(nx,ny,nz)
      real(kind=8) chi_rk1(nx,ny,nz)
      real(kind=8) trK_rk1(nx,ny,nz)
      real(kind=8) Gam1_rk1(nx,ny,nz)
      real(kind=8) Gam2_rk1(nx,ny,nz)
      real(kind=8) Gam3_rk1(nx,ny,nz)
      real(kind=8) alpha_rk1(nx,ny,nz)
      real(kind=8) shift1_rk1(nx,ny,nz)
      real(kind=8) shift2_rk1(nx,ny,nz)
      real(kind=8) shift3_rk1(nx,ny,nz)
      real(kind=8) gb1_rk1(nx,ny,nz)
      real(kind=8) gb2_rk1(nx,ny,nz)
      real(kind=8) gb3_rk1(nx,ny,nz)
      real(kind=8) Ex_rk1(nx,ny,nz)
      real(kind=8) Ey_rk1(nx,ny,nz)
      real(kind=8) Ez_rk1(nx,ny,nz)
      real(kind=8) Bx_rk1(nx,ny,nz)
      real(kind=8) By_rk1(nx,ny,nz)
      real(kind=8) Bz_rk1(nx,ny,nz)
      real(kind=8) Phi_em_rk1(nx,ny,nz)
      real(kind=8) Psi_em_rk1(nx,ny,nz)
      real(kind=8) phiR_rk1(nx,ny,nz)
      real(kind=8) phiI_rk1(nx,ny,nz)
      real(kind=8) piR_rk1(nx,ny,nz)
      real(kind=8) piI_rk1(nx,ny,nz)
        target :: gt11_rk1
        target :: gt12_rk1
        target :: gt13_rk1
        target :: gt22_rk1
        target :: gt23_rk1
        target :: gt33_rk1
        target :: A11_rk1
        target :: A12_rk1
        target :: A13_rk1
        target :: A22_rk1
        target :: A23_rk1
        target :: A33_rk1
        target :: chi_rk1
        target :: trK_rk1
        target :: Gam1_rk1
        target :: Gam2_rk1
        target :: Gam3_rk1
        target :: alpha_rk1
        target :: shift1_rk1
        target :: shift2_rk1
        target :: shift3_rk1
        target :: gb1_rk1
        target :: gb2_rk1
        target :: gb3_rk1
        target :: Ex_rk1
        target :: Ey_rk1
        target :: Ez_rk1
        target :: Bx_rk1
        target :: By_rk1
        target :: Bz_rk1
        target :: Phi_em_rk1
        target :: Psi_em_rk1
        target :: phiR_rk1
        target :: phiI_rk1
        target :: piR_rk1
        target :: piI_rk1
        urk1(H_GT11)%d => gt11_rk1
        urk1(H_GT12)%d => gt12_rk1
        urk1(H_GT13)%d => gt13_rk1
        urk1(H_GT22)%d => gt22_rk1
        urk1(H_GT23)%d => gt23_rk1
        urk1(H_GT33)%d => gt33_rk1
        urk1(H_A11)%d => A11_rk1
        urk1(H_A12)%d => A12_rk1
        urk1(H_A13)%d => A13_rk1
        urk1(H_A22)%d => A22_rk1
        urk1(H_A23)%d => A23_rk1
        urk1(H_A33)%d => A33_rk1
        urk1(H_CHI)%d => chi_rk1
        urk1(H_TRK)%d => trK_rk1
        urk1(H_GAM1)%d => Gam1_rk1
        urk1(H_GAM2)%d => Gam2_rk1
        urk1(H_GAM3)%d => Gam3_rk1
        urk1(H_ALPHA)%d => alpha_rk1
        urk1(H_SHIFT1)%d => shift1_rk1
        urk1(H_SHIFT2)%d => shift2_rk1
        urk1(H_SHIFT3)%d => shift3_rk1
        urk1(H_GB1)%d => gb1_rk1
        urk1(H_GB2)%d => gb2_rk1
        urk1(H_GB3)%d => gb3_rk1
        urk1(H_EX)%d => Ex_rk1
        urk1(H_EY)%d => Ey_rk1
        urk1(H_EZ)%d => Ez_rk1
        urk1(H_BX)%d => Bx_rk1
        urk1(H_BY)%d => By_rk1
        urk1(H_BZ)%d => Bz_rk1
        urk1(H_PHI_EM)%d => Phi_em_rk1
        urk1(H_PSI_EM)%d => Psi_em_rk1
        urk1(H_PHIR)%d => phiR_rk1
        urk1(H_PHII)%d => phiI_rk1
        urk1(H_PIR)%d => piR_rk1
        urk1(H_PII)%d => piI_rk1
      return
      end subroutine assign_ptrs_rks
      subroutine assign_ptrs_derivs(dxu, dyu, dzu, dxv, dyv, dzv, nx, ny
     &, nz)
      use GF
      implicit none
      integer :: nx, ny, nz
      type(gridfunction), dimension(NU) :: dxu, dyu, dzu
      type(gridfunction), dimension(NV) :: dxv, dyv, dzv
      nullify(dxu(H_GT11)%d)
      nullify(dyu(H_GT11)%d)
      nullify(dzu(H_GT11)%d)
      nullify(dxu(H_GT12)%d)
      nullify(dyu(H_GT12)%d)
      nullify(dzu(H_GT12)%d)
      nullify(dxu(H_GT13)%d)
      nullify(dyu(H_GT13)%d)
      nullify(dzu(H_GT13)%d)
      nullify(dxu(H_GT22)%d)
      nullify(dyu(H_GT22)%d)
      nullify(dzu(H_GT22)%d)
      nullify(dxu(H_GT23)%d)
      nullify(dyu(H_GT23)%d)
      nullify(dzu(H_GT23)%d)
      nullify(dxu(H_GT33)%d)
      nullify(dyu(H_GT33)%d)
      nullify(dzu(H_GT33)%d)
      nullify(dxu(H_A11)%d)
      nullify(dyu(H_A11)%d)
      nullify(dzu(H_A11)%d)
      nullify(dxu(H_A12)%d)
      nullify(dyu(H_A12)%d)
      nullify(dzu(H_A12)%d)
      nullify(dxu(H_A13)%d)
      nullify(dyu(H_A13)%d)
      nullify(dzu(H_A13)%d)
      nullify(dxu(H_A22)%d)
      nullify(dyu(H_A22)%d)
      nullify(dzu(H_A22)%d)
      nullify(dxu(H_A23)%d)
      nullify(dyu(H_A23)%d)
      nullify(dzu(H_A23)%d)
      nullify(dxu(H_A33)%d)
      nullify(dyu(H_A33)%d)
      nullify(dzu(H_A33)%d)
      nullify(dxu(H_CHI)%d)
      nullify(dyu(H_CHI)%d)
      nullify(dzu(H_CHI)%d)
      nullify(dxu(H_TRK)%d)
      nullify(dyu(H_TRK)%d)
      nullify(dzu(H_TRK)%d)
      nullify(dxu(H_GAM1)%d)
      nullify(dyu(H_GAM1)%d)
      nullify(dzu(H_GAM1)%d)
      nullify(dxu(H_GAM2)%d)
      nullify(dyu(H_GAM2)%d)
      nullify(dzu(H_GAM2)%d)
      nullify(dxu(H_GAM3)%d)
      nullify(dyu(H_GAM3)%d)
      nullify(dzu(H_GAM3)%d)
      nullify(dxu(H_ALPHA)%d)
      nullify(dyu(H_ALPHA)%d)
      nullify(dzu(H_ALPHA)%d)
      nullify(dxu(H_SHIFT1)%d)
      nullify(dyu(H_SHIFT1)%d)
      nullify(dzu(H_SHIFT1)%d)
      nullify(dxu(H_SHIFT2)%d)
      nullify(dyu(H_SHIFT2)%d)
      nullify(dzu(H_SHIFT2)%d)
      nullify(dxu(H_SHIFT3)%d)
      nullify(dyu(H_SHIFT3)%d)
      nullify(dzu(H_SHIFT3)%d)
      nullify(dxu(H_GB1)%d)
      nullify(dyu(H_GB1)%d)
      nullify(dzu(H_GB1)%d)
      nullify(dxu(H_GB2)%d)
      nullify(dyu(H_GB2)%d)
      nullify(dzu(H_GB2)%d)
      nullify(dxu(H_GB3)%d)
      nullify(dyu(H_GB3)%d)
      nullify(dzu(H_GB3)%d)
      nullify(dxu(H_EX)%d)
      nullify(dyu(H_EX)%d)
      nullify(dzu(H_EX)%d)
      nullify(dxu(H_EY)%d)
      nullify(dyu(H_EY)%d)
      nullify(dzu(H_EY)%d)
      nullify(dxu(H_EZ)%d)
      nullify(dyu(H_EZ)%d)
      nullify(dzu(H_EZ)%d)
      nullify(dxu(H_BX)%d)
      nullify(dyu(H_BX)%d)
      nullify(dzu(H_BX)%d)
      nullify(dxu(H_BY)%d)
      nullify(dyu(H_BY)%d)
      nullify(dzu(H_BY)%d)
      nullify(dxu(H_BZ)%d)
      nullify(dyu(H_BZ)%d)
      nullify(dzu(H_BZ)%d)
      nullify(dxu(H_PHI_EM)%d)
      nullify(dyu(H_PHI_EM)%d)
      nullify(dzu(H_PHI_EM)%d)
      nullify(dxu(H_PSI_EM)%d)
      nullify(dyu(H_PSI_EM)%d)
      nullify(dzu(H_PSI_EM)%d)
      nullify(dxu(H_PHIR)%d)
      nullify(dyu(H_PHIR)%d)
      nullify(dzu(H_PHIR)%d)
      nullify(dxu(H_PHII)%d)
      nullify(dyu(H_PHII)%d)
      nullify(dzu(H_PHII)%d)
      nullify(dxu(H_PIR)%d)
      nullify(dyu(H_PIR)%d)
      nullify(dzu(H_PIR)%d)
      nullify(dxu(H_PII)%d)
      nullify(dyu(H_PII)%d)
      nullify(dzu(H_PII)%d)
      nullify(dxv(H_G11)%d)
      nullify(dyv(H_G11)%d)
      nullify(dzv(H_G11)%d)
      nullify(dxv(H_G12)%d)
      nullify(dyv(H_G12)%d)
      nullify(dzv(H_G12)%d)
      nullify(dxv(H_G13)%d)
      nullify(dyv(H_G13)%d)
      nullify(dzv(H_G13)%d)
      nullify(dxv(H_G22)%d)
      nullify(dyv(H_G22)%d)
      nullify(dzv(H_G22)%d)
      nullify(dxv(H_G23)%d)
      nullify(dyv(H_G23)%d)
      nullify(dzv(H_G23)%d)
      nullify(dxv(H_G33)%d)
      nullify(dyv(H_G33)%d)
      nullify(dzv(H_G33)%d)
      nullify(dxv(H_SDETG)%d)
      nullify(dyv(H_SDETG)%d)
      nullify(dzv(H_SDETG)%d)
      nullify(dxv(H_RAD_EXP)%d)
      nullify(dyv(H_RAD_EXP)%d)
      nullify(dzv(H_RAD_EXP)%d)
      nullify(dxv(H_PSI4R)%d)
      nullify(dyv(H_PSI4R)%d)
      nullify(dzv(H_PSI4R)%d)
      nullify(dxv(H_PSI4I)%d)
      nullify(dyv(H_PSI4I)%d)
      nullify(dzv(H_PSI4I)%d)
      nullify(dxv(H_MASSADM)%d)
      nullify(dyv(H_MASSADM)%d)
      nullify(dzv(H_MASSADM)%d)
      nullify(dxv(H_MASSBONDI)%d)
      nullify(dyv(H_MASSBONDI)%d)
      nullify(dzv(H_MASSBONDI)%d)
      nullify(dxv(H_CURVATURE)%d)
      nullify(dyv(H_CURVATURE)%d)
      nullify(dzv(H_CURVATURE)%d)
      nullify(dxv(H_PHI2R)%d)
      nullify(dyv(H_PHI2R)%d)
      nullify(dzv(H_PHI2R)%d)
      nullify(dxv(H_PHI2I)%d)
      nullify(dyv(H_PHI2I)%d)
      nullify(dzv(H_PHI2I)%d)
      nullify(dxv(H_PHI0R)%d)
      nullify(dyv(H_PHI0R)%d)
      nullify(dzv(H_PHI0R)%d)
      nullify(dxv(H_PHI0I)%d)
      nullify(dyv(H_PHI0I)%d)
      nullify(dzv(H_PHI0I)%d)
      nullify(dxv(H_VBR_OMEGA)%d)
      nullify(dyv(H_VBR_OMEGA)%d)
      nullify(dzv(H_VBR_OMEGA)%d)
      nullify(dxv(H_VBTH_OMEGA)%d)
      nullify(dyv(H_VBTH_OMEGA)%d)
      nullify(dzv(H_VBTH_OMEGA)%d)
      nullify(dxv(H_CHARGE)%d)
      nullify(dyv(H_CHARGE)%d)
      nullify(dzv(H_CHARGE)%d)
      nullify(dxv(H_J1)%d)
      nullify(dyv(H_J1)%d)
      nullify(dzv(H_J1)%d)
      nullify(dxv(H_J2)%d)
      nullify(dyv(H_J2)%d)
      nullify(dzv(H_J2)%d)
      nullify(dxv(H_J3)%d)
      nullify(dyv(H_J3)%d)
      nullify(dzv(H_J3)%d)
      nullify(dxv(H_POYNTINGX_DENS)%d)
      nullify(dyv(H_POYNTINGX_DENS)%d)
      nullify(dzv(H_POYNTINGX_DENS)%d)
      nullify(dxv(H_POYNTINGY_DENS)%d)
      nullify(dyv(H_POYNTINGY_DENS)%d)
      nullify(dzv(H_POYNTINGY_DENS)%d)
      nullify(dxv(H_POYNTINGZ_DENS)%d)
      nullify(dyv(H_POYNTINGZ_DENS)%d)
      nullify(dzv(H_POYNTINGZ_DENS)%d)
      nullify(dxv(H_UELL)%d)
      nullify(dyv(H_UELL)%d)
      nullify(dzv(H_UELL)%d)
      return
      end subroutine assign_ptrs_derivs
      subroutine assign_params(gauge_type,idtype,evolve_geometry,evolve_
     &em_field,evolve_scalar_field,temperature_type,read_file_alp,read_f
     &ile_b,read_file_g,read_file_k,read_file_psi,read_file_phir,read_fi
     &le_phic,read_file_phim,read_data_level,id_d_diff_order,mr_amp,mr_s
     &hift,bamp,e1_amp,e2_amp,initial_b,sf_amp,initial_sf,sf_amp_axn,dil
     &_mass,dil_alpha,dil_infty,emd_bh_type,axn_mass,axn_alpha,axn_infty
     &,sen_alpha,a_ang_par,q_elec,p_mag,bs_1_x0,bs_1_y0,bs_1_z0,bs_1_vx,
     &bs_1_vy,bs_1_vz,bs_1_omega,bs_2_x0,bs_2_y0,bs_2_z0,bs_2_vx,bs_2_vy
     &,bs_2_vz,bs_2_omega,bh_n,bh_type,bh_1_mass,bh_2_mass,bh_3_mass,bh_
     &4_mass,bh_1_x,bh_1_y,bh_1_z,bh_2_x,bh_2_y,bh_2_z,bh_3_x,bh_3_y,bh_
     &3_z,bh_4_x,bh_4_y,bh_4_z,bh_1_px,bh_1_py,bh_1_pz,bh_2_px,bh_2_py,b
     &h_2_pz,bh_3_px,bh_3_py,bh_3_pz,bh_4_px,bh_4_py,bh_4_pz,bh_1_spin,b
     &h_1_sth,bh_1_sphi,bh_2_spin,bh_2_sth,bh_2_sphi,bh_3_spin,bh_3_sth,
     &bh_3_sphi,bh_4_spin,bh_4_sth,bh_4_sphi,outer_boundary,gr_bound_con
     &d,fluid_bound_cond,bssn_lambda_1,bssn_lambda_2,bssn_lambda_3,bssn_
     &lambda_4,bssn_lambda_f0,bssn_lambda_f1,bssn_lambda_f2,bssn_lambda_
     &f3,bssn_eta_damping,bssn_R_0,bssn_eta_damping_exp,bssn_eta,bssn_tr
     &k0,bssn_kappa1,bssn_kappa2,bssn_chi_floor,bssn_adv_derivs,geometry
     &_type,gr_idtype,mhd_idtype,detgwarn,force_free,q1,vinj,Kappajan,dj
     &ump,rpeak,rout,stab,rho_real,magcase,kappa_max,calcDivB,constraint
     &s_analysis,G_scale_factor,B_scale_factor,anti_aligned,bfield_vacuu
     &mfreeze,bfield_taperrange,project_div_B,psi_ch,psi_cr,damp,id_bond
     &i_sonic_r,id_bondi_rhoc,id_bondi_rstart,id_bondi_fixedIBC,id_tov_d
     &epletion,id_tovbh_vx,id_tovbh_vy,id_tovbh_vz,id_tov_average,id_tov
     &_p_pres,id_tov_p_r,id_tov_p_sigma,id_tov_magnetic,id_tov_Asize,id_
     &drns_r_e,id_center_x1,id_center_y1,id_center_x2,id_center_y2,id_vx
     &1,id_vy1,id_vx2,id_vy2,id_perturb_m,id_perturb_p,id_perturb_rho,id
     &_disk_width,id_disk_rin,id_disk_rout,id_disk_B,id_disk_type,id_dis
     &k_angmom,id_disk_potential,id_disk_atmosphere,id_disk_kappa,id_dis
     &k_Asize,id_disk_rhocut,id_disk_magnetic,id_disk_rho0,id_disk_decay
     &,id_disk_kickvel,id_disk_kicktheta,id_disk_kickphi,id_disk_GammaB,
     &id_disk_cB,id_disk_lbound,id_disk_bhmass,nx0,ny0,nz0,nt0,h,hx,hy,h
     &z,hxyz0,run_wtime,amp,idata,maxchi_thresh,maxchi_minctime,local_nx
     &,local_ny,local_nz,global_nx,global_ny,global_nz,local_lower_bnd_x
     &,local_lower_bnd_y,local_lower_bnd_z,bbox1,bbox2,bbox3,bbox4,bbox5
     &,bbox6,nghostzones_x,nghostzones_y,nghostzones_z,dt,dx,dy,dz,local
     &_time,alt_coord_type,pc_coord_trans_rad,pc_coord_trans_width,pc_co
     &ord_scale,bc_type,inner_bound_data,deriv_order,dissipation,sigma_d
     &iss,extradiss,extradissOUT,nbholes,use_mask,mask_type,initial_anal
     &ysis,bh1_mass,bh1_spin,bh1_spin_phi,bh1_spin_theta,bh1_x0,bh1_y0,b
     &h1_z0,bh1_vx,bh1_vy,bh1_exc_rad,bh1_velx,bh2_mass,bh2_spin,bh2_spi
     &n_phi,bh2_spin_theta,bh2_x0,bh2_y0,bh2_z0,bh2_exc_rad,bh2_velx,bou
     &ndary_conditions,penalty,PP,QQ,t0,sigma_t,amp_boundary,amp_random_
     &bc,sigma_rho,interp_id,runge_kutta_type,runge_kutta_bound,point_wi
     &se_analysis,psi4_analysis,rk_iter,asf_ntheta,asf_nphi,asf_period,a
     &sf_level,asf_rconst,bsf_ntheta,bsf_nphi,bsf_period,bsf_level,bsf_r
     &const,csf_ntheta,csf_nphi,csf_period,csf_level,csf_rconst,lambda,r
     &efine_factor,refine_period_ctrl,refine_period,refine_deltat,simple
     &FMR,clusterDD,clusterstyle,diss_afterinj,allowedl,linearbounds,sha
     &dow,minx0,miny0,minz0,maxx0,maxy0,maxz0,trace_level,ethreshold,buf
     &fer,mindim,window,minefficiency,output_style,output_dim,output_f1_
     &type,output_f1_level,output_f1_period,output_f1_lb1,output_f1_lb2,
     &output_f1_lb3,output_f1_ub1,output_f1_ub2,output_f1_ub3,output_f2_
     &type,output_f2_level,output_f2_period,output_f2_lb1,output_f2_lb2,
     &output_f2_lb3,output_f2_ub1,output_f2_ub2,output_f2_ub3,output_f3_
     &type,output_f3_level,output_f3_period,output_f3_lb1,output_f3_lb2,
     &output_f3_lb3,output_f3_ub1,output_f3_ub2,output_f3_ub3,output_f4_
     &type,output_f4_level,output_f4_period,output_f4_lb1,output_f4_lb2,
     &output_f4_lb3,output_f4_ub1,output_f4_ub2,output_f4_ub3,clusterrea
     &dwrite,ghostwidth,update_scheme,amrbound_prepost,amrbound_timealig
     &n,elliptic_solve,nvcycle,preswp,pstswp,maxsweeps,ell_epsilon,num_e
     &vol_iters,chkpt_period,chkpt_readstate,chkpt_control,bound_width,w
     &eno_interp,findhorizon,mask_period,mask_usemin,mask_minfield,horiz
     &on_ntheta,horizon_nphi,horizon_thresh,horizon_nholes,horizon_recen
     &terp,horizon_growth,emulate_proc,assume_symmetry,flush_period,trac
     &ers_period,tracers_initial,tracers_scheme,variable_timestep,cfl_la
     &mbda,periodicBC,hg,dtg,minx,miny,minz,maxx,maxy,maxz,time,par)
      implicit none
      real(kind=8), dimension(NPAR) :: par
      real(kind=8) hg,dtg,minx,miny,minz,maxx,maxy,maxz,time
       integer gauge_type
       integer idtype
       integer evolve_geometry
       integer evolve_em_field
       integer evolve_scalar_field
       integer temperature_type
       integer read_file_alp
       integer read_file_b
       integer read_file_g
       integer read_file_k
       integer read_file_psi
       integer read_file_phir
       integer read_file_phic
       integer read_file_phim
       integer read_data_level
       integer id_d_diff_order
       real(kind=8) mr_amp
       real(kind=8) mr_shift
       real(kind=8) bamp
       real(kind=8) e1_amp
       real(kind=8) e2_amp
       integer initial_b
       real(kind=8) sf_amp
       integer initial_sf
       real(kind=8) sf_amp_axn
       real(kind=8) dil_mass
       real(kind=8) dil_alpha
       real(kind=8) dil_infty
       integer emd_bh_type
       real(kind=8) axn_mass
       real(kind=8) axn_alpha
       real(kind=8) axn_infty
       real(kind=8) sen_alpha
       real(kind=8) a_ang_par
       real(kind=8) q_elec
       real(kind=8) p_mag
       real(kind=8) bs_1_x0
       real(kind=8) bs_1_y0
       real(kind=8) bs_1_z0
       real(kind=8) bs_1_vx
       real(kind=8) bs_1_vy
       real(kind=8) bs_1_vz
       real(kind=8) bs_1_omega
       real(kind=8) bs_2_x0
       real(kind=8) bs_2_y0
       real(kind=8) bs_2_z0
       real(kind=8) bs_2_vx
       real(kind=8) bs_2_vy
       real(kind=8) bs_2_vz
       real(kind=8) bs_2_omega
       integer bh_n
       integer bh_type
       real(kind=8) bh_1_mass
       real(kind=8) bh_2_mass
       real(kind=8) bh_3_mass
       real(kind=8) bh_4_mass
       real(kind=8) bh_1_x
       real(kind=8) bh_1_y
       real(kind=8) bh_1_z
       real(kind=8) bh_2_x
       real(kind=8) bh_2_y
       real(kind=8) bh_2_z
       real(kind=8) bh_3_x
       real(kind=8) bh_3_y
       real(kind=8) bh_3_z
       real(kind=8) bh_4_x
       real(kind=8) bh_4_y
       real(kind=8) bh_4_z
       real(kind=8) bh_1_px
       real(kind=8) bh_1_py
       real(kind=8) bh_1_pz
       real(kind=8) bh_2_px
       real(kind=8) bh_2_py
       real(kind=8) bh_2_pz
       real(kind=8) bh_3_px
       real(kind=8) bh_3_py
       real(kind=8) bh_3_pz
       real(kind=8) bh_4_px
       real(kind=8) bh_4_py
       real(kind=8) bh_4_pz
       real(kind=8) bh_1_spin
       real(kind=8) bh_1_sth
       real(kind=8) bh_1_sphi
       real(kind=8) bh_2_spin
       real(kind=8) bh_2_sth
       real(kind=8) bh_2_sphi
       real(kind=8) bh_3_spin
       real(kind=8) bh_3_sth
       real(kind=8) bh_3_sphi
       real(kind=8) bh_4_spin
       real(kind=8) bh_4_sth
       real(kind=8) bh_4_sphi
       real(kind=8) outer_boundary
       real(kind=8) gr_bound_cond
       real(kind=8) fluid_bound_cond
       integer bssn_lambda_1
       integer bssn_lambda_2
       integer bssn_lambda_3
       integer bssn_lambda_4
       real(kind=8) bssn_lambda_f0
       real(kind=8) bssn_lambda_f1
       real(kind=8) bssn_lambda_f2
       real(kind=8) bssn_lambda_f3
       real(kind=8) bssn_eta_damping
       real(kind=8) bssn_R_0
       real(kind=8) bssn_eta_damping_exp
       real(kind=8) bssn_eta
       real(kind=8) bssn_trk0
       real(kind=8) bssn_kappa1
       real(kind=8) bssn_kappa2
       real(kind=8) bssn_chi_floor
       integer bssn_adv_derivs
       integer geometry_type
       integer gr_idtype
       integer mhd_idtype
       integer detgwarn
       integer force_free
       real(kind=8) q1
       real(kind=8) vinj
       real(kind=8) Kappajan
       real(kind=8) djump
       real(kind=8) rpeak
       real(kind=8) rout
       real(kind=8) stab
       real(kind=8) rho_real
       real(kind=8) magcase
       real(kind=8) kappa_max
       integer calcDivB
       integer constraints_analysis
       real(kind=8) G_scale_factor
       real(kind=8) B_scale_factor
       integer anti_aligned
       integer bfield_vacuumfreeze
       real(kind=8) bfield_taperrange
       integer project_div_B
       real(kind=8) psi_ch
       real(kind=8) psi_cr
       integer damp
       real(kind=8) id_bondi_sonic_r
       real(kind=8) id_bondi_rhoc
       real(kind=8) id_bondi_rstart
       integer id_bondi_fixedIBC
       real(kind=8) id_tov_depletion
       real(kind=8) id_tovbh_vx
       real(kind=8) id_tovbh_vy
       real(kind=8) id_tovbh_vz
       integer id_tov_average
       real(kind=8) id_tov_p_pres
       real(kind=8) id_tov_p_r
       real(kind=8) id_tov_p_sigma
       integer id_tov_magnetic
       real(kind=8) id_tov_Asize
       real(kind=8) id_drns_r_e
       real(kind=8) id_center_x1
       real(kind=8) id_center_y1
       real(kind=8) id_center_x2
       real(kind=8) id_center_y2
       real(kind=8) id_vx1
       real(kind=8) id_vy1
       real(kind=8) id_vx2
       real(kind=8) id_vy2
       integer id_perturb_m
       real(kind=8) id_perturb_p
       real(kind=8) id_perturb_rho
       real(kind=8) id_disk_width
       real(kind=8) id_disk_rin
       real(kind=8) id_disk_rout
       real(kind=8) id_disk_B
       integer id_disk_type
       real(kind=8) id_disk_angmom
       real(kind=8) id_disk_potential
       real(kind=8) id_disk_atmosphere
       real(kind=8) id_disk_kappa
       real(kind=8) id_disk_Asize
       real(kind=8) id_disk_rhocut
       integer id_disk_magnetic
       real(kind=8) id_disk_rho0
       real(kind=8) id_disk_decay
       real(kind=8) id_disk_kickvel
       real(kind=8) id_disk_kicktheta
       real(kind=8) id_disk_kickphi
       real(kind=8) id_disk_GammaB
       real(kind=8) id_disk_cB
       real(kind=8) id_disk_lbound
       real(kind=8) id_disk_bhmass
       integer nx0
       integer ny0
       integer nz0
       integer nt0
       real(kind=8) h
       real(kind=8) hx
       real(kind=8) hy
       real(kind=8) hz
       real(kind=8) hxyz0
       real(kind=8) run_wtime
       real(kind=8) amp
       integer idata
       real(kind=8) maxchi_thresh
       real(kind=8) maxchi_minctime
       integer local_nx
       integer local_ny
       integer local_nz
       integer global_nx
       integer global_ny
       integer global_nz
       integer local_lower_bnd_x
       integer local_lower_bnd_y
       integer local_lower_bnd_z
       integer bbox1
       integer bbox2
       integer bbox3
       integer bbox4
       integer bbox5
       integer bbox6
       integer nghostzones_x
       integer nghostzones_y
       integer nghostzones_z
       real(kind=8) dt
       real(kind=8) dx
       real(kind=8) dy
       real(kind=8) dz
       real(kind=8) local_time
       integer alt_coord_type
       real(kind=8) pc_coord_trans_rad
       real(kind=8) pc_coord_trans_width
       real(kind=8) pc_coord_scale
       integer bc_type
       integer inner_bound_data
       integer deriv_order
       integer dissipation
       real(kind=8) sigma_diss
       real(kind=8) extradiss
       real(kind=8) extradissOUT
       integer nbholes
       integer use_mask
       integer mask_type
       integer initial_analysis
       real(kind=8) bh1_mass
       real(kind=8) bh1_spin
       real(kind=8) bh1_spin_phi
       real(kind=8) bh1_spin_theta
       real(kind=8) bh1_x0
       real(kind=8) bh1_y0
       real(kind=8) bh1_z0
       real(kind=8) bh1_vx
       real(kind=8) bh1_vy
       real(kind=8) bh1_exc_rad
       real(kind=8) bh1_velx
       real(kind=8) bh2_mass
       real(kind=8) bh2_spin
       real(kind=8) bh2_spin_phi
       real(kind=8) bh2_spin_theta
       real(kind=8) bh2_x0
       real(kind=8) bh2_y0
       real(kind=8) bh2_z0
       real(kind=8) bh2_exc_rad
       real(kind=8) bh2_velx
       integer boundary_conditions
       real(kind=8) penalty
       integer PP
       integer QQ
       real(kind=8) t0
       real(kind=8) sigma_t
       real(kind=8) amp_boundary
       real(kind=8) amp_random_bc
       real(kind=8) sigma_rho
       integer interp_id
       integer runge_kutta_type
       integer runge_kutta_bound
       integer point_wise_analysis
       integer psi4_analysis
       integer rk_iter
       integer asf_ntheta
       integer asf_nphi
       integer asf_period
       integer asf_level
       real(kind=8) asf_rconst
       integer bsf_ntheta
       integer bsf_nphi
       integer bsf_period
       integer bsf_level
       real(kind=8) bsf_rconst
       integer csf_ntheta
       integer csf_nphi
       integer csf_period
       integer csf_level
       real(kind=8) csf_rconst
       real(kind=8) lambda
       integer refine_factor
       integer refine_period_ctrl
       integer refine_period
       real(kind=8) refine_deltat
       integer simpleFMR
       integer clusterDD
       integer clusterstyle
       integer diss_afterinj
       integer allowedl
       integer linearbounds
       integer shadow
       real(kind=8) minx0
       real(kind=8) miny0
       real(kind=8) minz0
       real(kind=8) maxx0
       real(kind=8) maxy0
       real(kind=8) maxz0
       integer trace_level
       real(kind=8) ethreshold
       integer buffer
       integer mindim
       integer window
       real(kind=8) minefficiency
       integer output_style
       integer output_dim
       integer output_f1_type
       integer output_f1_level
       integer output_f1_period
       real(kind=8) output_f1_lb1
       real(kind=8) output_f1_lb2
       real(kind=8) output_f1_lb3
       real(kind=8) output_f1_ub1
       real(kind=8) output_f1_ub2
       real(kind=8) output_f1_ub3
       integer output_f2_type
       integer output_f2_level
       integer output_f2_period
       real(kind=8) output_f2_lb1
       real(kind=8) output_f2_lb2
       real(kind=8) output_f2_lb3
       real(kind=8) output_f2_ub1
       real(kind=8) output_f2_ub2
       real(kind=8) output_f2_ub3
       integer output_f3_type
       integer output_f3_level
       integer output_f3_period
       real(kind=8) output_f3_lb1
       real(kind=8) output_f3_lb2
       real(kind=8) output_f3_lb3
       real(kind=8) output_f3_ub1
       real(kind=8) output_f3_ub2
       real(kind=8) output_f3_ub3
       integer output_f4_type
       integer output_f4_level
       integer output_f4_period
       real(kind=8) output_f4_lb1
       real(kind=8) output_f4_lb2
       real(kind=8) output_f4_lb3
       real(kind=8) output_f4_ub1
       real(kind=8) output_f4_ub2
       real(kind=8) output_f4_ub3
       integer clusterreadwrite
       integer ghostwidth
       integer update_scheme
       integer amrbound_prepost
       integer amrbound_timealign
       integer elliptic_solve
       integer nvcycle
       integer preswp
       integer pstswp
       integer maxsweeps
       real(kind=8) ell_epsilon
       integer num_evol_iters
       integer chkpt_period
       integer chkpt_readstate
       integer chkpt_control
       integer bound_width
       integer weno_interp
       integer findhorizon
       integer mask_period
       integer mask_usemin
       integer mask_minfield
       integer horizon_ntheta
       integer horizon_nphi
       real(kind=8) horizon_thresh
       integer horizon_nholes
       integer horizon_recenterp
       real(kind=8) horizon_growth
       integer emulate_proc
       integer assume_symmetry
       integer flush_period
       integer tracers_period
       integer tracers_initial
       integer tracers_scheme
       integer variable_timestep
       real(kind=8) cfl_lambda
       integer periodicBC
      
      
      
      
      
      global_nx = nx0
      global_ny = ny0
      global_nz = nz0
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      dt = dt
      dx = h
      dy = h
      dz = h
      nghostzones_x = bound_width
      nghostzones_y = bound_width
      nghostzones_z = bound_width
      local_time = time
      par(P_GAUGE_TYPE) = gauge_type
      par(P_IDTYPE) = idtype
      par(P_EVOLVE_GEOMETRY) = evolve_geometry
      par(P_EVOLVE_EM_FIELD) = evolve_em_field
      par(P_EVOLVE_SCALAR_FIELD) = evolve_scalar_field
      par(P_TEMPERATURE_TYPE) = temperature_type
      par(P_READ_FILE_ALP) = read_file_alp
      par(P_READ_FILE_B) = read_file_b
      par(P_READ_FILE_G) = read_file_g
      par(P_READ_FILE_K) = read_file_k
      par(P_READ_FILE_PSI) = read_file_psi
      par(P_READ_FILE_PHIR) = read_file_phir
      par(P_READ_FILE_PHIC) = read_file_phic
      par(P_READ_FILE_PHIM) = read_file_phim
      par(P_READ_DATA_LEVEL) = read_data_level
      par(P_ID_D_DIFF_ORDER) = id_d_diff_order
      par(P_MR_AMP) = mr_amp
      par(P_MR_SHIFT) = mr_shift
      par(P_BAMP) = bamp
      par(P_E1_AMP) = e1_amp
      par(P_E2_AMP) = e2_amp
      par(P_INITIAL_B) = initial_b
      par(P_SF_AMP) = sf_amp
      par(P_INITIAL_SF) = initial_sf
      par(P_SF_AMP_AXN) = sf_amp_axn
      par(P_DIL_MASS) = dil_mass
      par(P_DIL_ALPHA) = dil_alpha
      par(P_DIL_INFTY) = dil_infty
      par(P_EMD_BH_TYPE) = emd_bh_type
      par(P_AXN_MASS) = axn_mass
      par(P_AXN_ALPHA) = axn_alpha
      par(P_AXN_INFTY) = axn_infty
      par(P_SEN_ALPHA) = sen_alpha
      par(P_A_ANG_PAR) = a_ang_par
      par(P_Q_ELEC) = q_elec
      par(P_P_MAG) = p_mag
      par(P_BS_1_X0) = bs_1_x0
      par(P_BS_1_Y0) = bs_1_y0
      par(P_BS_1_Z0) = bs_1_z0
      par(P_BS_1_VX) = bs_1_vx
      par(P_BS_1_VY) = bs_1_vy
      par(P_BS_1_VZ) = bs_1_vz
      par(P_BS_1_OMEGA) = bs_1_omega
      par(P_BS_2_X0) = bs_2_x0
      par(P_BS_2_Y0) = bs_2_y0
      par(P_BS_2_Z0) = bs_2_z0
      par(P_BS_2_VX) = bs_2_vx
      par(P_BS_2_VY) = bs_2_vy
      par(P_BS_2_VZ) = bs_2_vz
      par(P_BS_2_OMEGA) = bs_2_omega
      par(P_BH_N) = bh_n
      par(P_BH_TYPE) = bh_type
      par(P_BH_1_MASS) = bh_1_mass
      par(P_BH_2_MASS) = bh_2_mass
      par(P_BH_3_MASS) = bh_3_mass
      par(P_BH_4_MASS) = bh_4_mass
      par(P_BH_1_X) = bh_1_x
      par(P_BH_1_Y) = bh_1_y
      par(P_BH_1_Z) = bh_1_z
      par(P_BH_2_X) = bh_2_x
      par(P_BH_2_Y) = bh_2_y
      par(P_BH_2_Z) = bh_2_z
      par(P_BH_3_X) = bh_3_x
      par(P_BH_3_Y) = bh_3_y
      par(P_BH_3_Z) = bh_3_z
      par(P_BH_4_X) = bh_4_x
      par(P_BH_4_Y) = bh_4_y
      par(P_BH_4_Z) = bh_4_z
      par(P_BH_1_PX) = bh_1_px
      par(P_BH_1_PY) = bh_1_py
      par(P_BH_1_PZ) = bh_1_pz
      par(P_BH_2_PX) = bh_2_px
      par(P_BH_2_PY) = bh_2_py
      par(P_BH_2_PZ) = bh_2_pz
      par(P_BH_3_PX) = bh_3_px
      par(P_BH_3_PY) = bh_3_py
      par(P_BH_3_PZ) = bh_3_pz
      par(P_BH_4_PX) = bh_4_px
      par(P_BH_4_PY) = bh_4_py
      par(P_BH_4_PZ) = bh_4_pz
      par(P_BH_1_SPIN) = bh_1_spin
      par(P_BH_1_STH) = bh_1_sth
      par(P_BH_1_SPHI) = bh_1_sphi
      par(P_BH_2_SPIN) = bh_2_spin
      par(P_BH_2_STH) = bh_2_sth
      par(P_BH_2_SPHI) = bh_2_sphi
      par(P_BH_3_SPIN) = bh_3_spin
      par(P_BH_3_STH) = bh_3_sth
      par(P_BH_3_SPHI) = bh_3_sphi
      par(P_BH_4_SPIN) = bh_4_spin
      par(P_BH_4_STH) = bh_4_sth
      par(P_BH_4_SPHI) = bh_4_sphi
      par(P_OUTER_BOUNDARY) = outer_boundary
      par(P_GR_BOUND_COND) = gr_bound_cond
      par(P_FLUID_BOUND_COND) = fluid_bound_cond
      par(P_BSSN_LAMBDA_1) = bssn_lambda_1
      par(P_BSSN_LAMBDA_2) = bssn_lambda_2
      par(P_BSSN_LAMBDA_3) = bssn_lambda_3
      par(P_BSSN_LAMBDA_4) = bssn_lambda_4
      par(P_BSSN_LAMBDA_F0) = bssn_lambda_f0
      par(P_BSSN_LAMBDA_F1) = bssn_lambda_f1
      par(P_BSSN_LAMBDA_F2) = bssn_lambda_f2
      par(P_BSSN_LAMBDA_F3) = bssn_lambda_f3
      par(P_BSSN_ETA_DAMPING) = bssn_eta_damping
      par(P_BSSN_R_0) = bssn_R_0
      par(P_BSSN_ETA_DAMPING_EXP) = bssn_eta_damping_exp
      par(P_BSSN_ETA) = bssn_eta
      par(P_BSSN_TRK0) = bssn_trk0
      par(P_BSSN_KAPPA1) = bssn_kappa1
      par(P_BSSN_KAPPA2) = bssn_kappa2
      par(P_BSSN_CHI_FLOOR) = bssn_chi_floor
      par(P_BSSN_ADV_DERIVS) = bssn_adv_derivs
      par(P_GEOMETRY_TYPE) = geometry_type
      par(P_GR_IDTYPE) = gr_idtype
      par(P_MHD_IDTYPE) = mhd_idtype
      par(P_DETGWARN) = detgwarn
      par(P_FORCE_FREE) = force_free
      par(P_Q1) = q1
      par(P_VINJ) = vinj
      par(P_KAPPAJAN) = Kappajan
      par(P_DJUMP) = djump
      par(P_RPEAK) = rpeak
      par(P_ROUT) = rout
      par(P_STAB) = stab
      par(P_RHO_REAL) = rho_real
      par(P_MAGCASE) = magcase
      par(P_KAPPA_MAX) = kappa_max
      par(P_CALCDIVB) = calcDivB
      par(P_CONSTRAINTS_ANALYSIS) = constraints_analysis
      par(P_G_SCALE_FACTOR) = G_scale_factor
      par(P_B_SCALE_FACTOR) = B_scale_factor
      par(P_ANTI_ALIGNED) = anti_aligned
      par(P_BFIELD_VACUUMFREEZE) = bfield_vacuumfreeze
      par(P_BFIELD_TAPERRANGE) = bfield_taperrange
      par(P_PROJECT_DIV_B) = project_div_B
      par(P_PSI_CH) = psi_ch
      par(P_PSI_CR) = psi_cr
      par(P_DAMP) = damp
      par(P_ID_BONDI_SONIC_R) = id_bondi_sonic_r
      par(P_ID_BONDI_RHOC) = id_bondi_rhoc
      par(P_ID_BONDI_RSTART) = id_bondi_rstart
      par(P_ID_BONDI_FIXEDIBC) = id_bondi_fixedIBC
      par(P_ID_TOV_DEPLETION) = id_tov_depletion
      par(P_ID_TOVBH_VX) = id_tovbh_vx
      par(P_ID_TOVBH_VY) = id_tovbh_vy
      par(P_ID_TOVBH_VZ) = id_tovbh_vz
      par(P_ID_TOV_AVERAGE) = id_tov_average
      par(P_ID_TOV_P_PRES) = id_tov_p_pres
      par(P_ID_TOV_P_R) = id_tov_p_r
      par(P_ID_TOV_P_SIGMA) = id_tov_p_sigma
      par(P_ID_TOV_MAGNETIC) = id_tov_magnetic
      par(P_ID_TOV_ASIZE) = id_tov_Asize
      par(P_ID_DRNS_R_E) = id_drns_r_e
      par(P_ID_CENTER_X1) = id_center_x1
      par(P_ID_CENTER_Y1) = id_center_y1
      par(P_ID_CENTER_X2) = id_center_x2
      par(P_ID_CENTER_Y2) = id_center_y2
      par(P_ID_VX1) = id_vx1
      par(P_ID_VY1) = id_vy1
      par(P_ID_VX2) = id_vx2
      par(P_ID_VY2) = id_vy2
      par(P_ID_PERTURB_M) = id_perturb_m
      par(P_ID_PERTURB_P) = id_perturb_p
      par(P_ID_PERTURB_RHO) = id_perturb_rho
      par(P_ID_DISK_WIDTH) = id_disk_width
      par(P_ID_DISK_RIN) = id_disk_rin
      par(P_ID_DISK_ROUT) = id_disk_rout
      par(P_ID_DISK_B) = id_disk_B
      par(P_ID_DISK_TYPE) = id_disk_type
      par(P_ID_DISK_ANGMOM) = id_disk_angmom
      par(P_ID_DISK_POTENTIAL) = id_disk_potential
      par(P_ID_DISK_ATMOSPHERE) = id_disk_atmosphere
      par(P_ID_DISK_KAPPA) = id_disk_kappa
      par(P_ID_DISK_ASIZE) = id_disk_Asize
      par(P_ID_DISK_RHOCUT) = id_disk_rhocut
      par(P_ID_DISK_MAGNETIC) = id_disk_magnetic
      par(P_ID_DISK_RHO0) = id_disk_rho0
      par(P_ID_DISK_DECAY) = id_disk_decay
      par(P_ID_DISK_KICKVEL) = id_disk_kickvel
      par(P_ID_DISK_KICKTHETA) = id_disk_kicktheta
      par(P_ID_DISK_KICKPHI) = id_disk_kickphi
      par(P_ID_DISK_GAMMAB) = id_disk_GammaB
      par(P_ID_DISK_CB) = id_disk_cB
      par(P_ID_DISK_LBOUND) = id_disk_lbound
      par(P_ID_DISK_BHMASS) = id_disk_bhmass
      par(P_NX0) = nx0
      par(P_NY0) = ny0
      par(P_NZ0) = nz0
      par(P_NT0) = nt0
      par(P_H) = h
      par(P_HX) = hx
      par(P_HY) = hy
      par(P_HZ) = hz
      par(P_HXYZ0) = hxyz0
      par(P_RUN_WTIME) = run_wtime
      par(P_AMP) = amp
      par(P_IDATA) = idata
      par(P_MAXCHI_THRESH) = maxchi_thresh
      par(P_MAXCHI_MINCTIME) = maxchi_minctime
      par(P_NX) = local_nx
      par(P_NY) = local_ny
      par(P_NZ) = local_nz
      par(P_GLOBAL_NX) = global_nx
      par(P_GLOBAL_NY) = global_ny
      par(P_GLOBAL_NZ) = global_nz
      par(P_LOWER_BND_X) = local_lower_bnd_x
      par(P_LOWER_BND_Y) = local_lower_bnd_y
      par(P_LOWER_BND_Z) = local_lower_bnd_z
      par(P_BBOX1) = bbox1
      par(P_BBOX2) = bbox2
      par(P_BBOX3) = bbox3
      par(P_BBOX4) = bbox4
      par(P_BBOX5) = bbox5
      par(P_BBOX6) = bbox6
      par(P_NGHOSTZONES_X) = nghostzones_x
      par(P_NGHOSTZONES_Y) = nghostzones_y
      par(P_NGHOSTZONES_Z) = nghostzones_z
      par(P_DT) = dt
      par(P_DX) = dx
      par(P_DY) = dy
      par(P_DZ) = dz
      par(P_TIME) = local_time
      par(P_ALT_COORD_TYPE) = alt_coord_type
      par(P_PC_COORD_TRANS_RAD) = pc_coord_trans_rad
      par(P_PC_COORD_TRANS_WIDTH) = pc_coord_trans_width
      par(P_PC_COORD_SCALE) = pc_coord_scale
      par(P_BC_TYPE) = bc_type
      par(P_INNER_BOUND_DATA) = inner_bound_data
      par(P_DERIV_ORDER) = deriv_order
      par(P_DISSIPATION) = dissipation
      par(P_SIGMA_DISS) = sigma_diss
      par(P_EXTRADISS) = extradiss
      par(P_EXTRADISSOUT) = extradissOUT
      par(P_NBHOLES) = nbholes
      par(P_USE_MASK) = use_mask
      par(P_MASK_TYPE) = mask_type
      par(P_INITIAL_ANALYSIS) = initial_analysis
      par(P_BH1_MASS) = bh1_mass
      par(P_BH1_SPIN) = bh1_spin
      par(P_BH1_SPIN_PHI) = bh1_spin_phi
      par(P_BH1_SPIN_THETA) = bh1_spin_theta
      par(P_BH1_X0) = bh1_x0
      par(P_BH1_Y0) = bh1_y0
      par(P_BH1_Z0) = bh1_z0
      par(P_BH1_VX) = bh1_vx
      par(P_BH1_VY) = bh1_vy
      par(P_BH1_EXC_RAD) = bh1_exc_rad
      par(P_BH1_VELX) = bh1_velx
      par(P_BH2_MASS) = bh2_mass
      par(P_BH2_SPIN) = bh2_spin
      par(P_BH2_SPIN_PHI) = bh2_spin_phi
      par(P_BH2_SPIN_THETA) = bh2_spin_theta
      par(P_BH2_X0) = bh2_x0
      par(P_BH2_Y0) = bh2_y0
      par(P_BH2_Z0) = bh2_z0
      par(P_BH2_EXC_RAD) = bh2_exc_rad
      par(P_BH2_VELX) = bh2_velx
      par(P_BOUNDARY_CONDITIONS) = boundary_conditions
      par(P_PENALTY) = penalty
      par(P_PP) = PP
      par(P_QQ) = QQ
      par(P_T0) = t0
      par(P_SIGMA_T) = sigma_t
      par(P_AMP_BOUNDARY) = amp_boundary
      par(P_AMP_RANDOM_BC) = amp_random_bc
      par(P_SIGMA_RHO) = sigma_rho
      par(P_INTERP_ID) = interp_id
      par(P_RUNGE_KUTTA_TYPE) = runge_kutta_type
      par(P_RUNGE_KUTTA_BOUND) = runge_kutta_bound
      par(P_POINT_WISE_ANALYSIS) = point_wise_analysis
      par(P_PSI4_ANALYSIS) = psi4_analysis
      par(P_RK_ITER) = rk_iter
      par(P_ASF_NTHETA) = asf_ntheta
      par(P_ASF_NPHI) = asf_nphi
      par(P_ASF_PERIOD) = asf_period
      par(P_ASF_LEVEL) = asf_level
      par(P_ASF_RCONST) = asf_rconst
      par(P_BSF_NTHETA) = bsf_ntheta
      par(P_BSF_NPHI) = bsf_nphi
      par(P_BSF_PERIOD) = bsf_period
      par(P_BSF_LEVEL) = bsf_level
      par(P_BSF_RCONST) = bsf_rconst
      par(P_CSF_NTHETA) = csf_ntheta
      par(P_CSF_NPHI) = csf_nphi
      par(P_CSF_PERIOD) = csf_period
      par(P_CSF_LEVEL) = csf_level
      par(P_CSF_RCONST) = csf_rconst
      par(P_LAMBDA) = lambda
      par(P_REFINE_FACTOR) = refine_factor
      par(P_REFINE_PERIOD_CTRL) = refine_period_ctrl
      par(P_REFINE_PERIOD) = refine_period
      par(P_REFINE_DELTAT) = refine_deltat
      par(P_SIMPLEFMR) = simpleFMR
      par(P_CLUSTERDD) = clusterDD
      par(P_CLUSTERSTYLE) = clusterstyle
      par(P_DISS_AFTERINJ) = diss_afterinj
      par(P_ALLOWEDL) = allowedl
      par(P_LINEARBOUNDS) = linearbounds
      par(P_SHADOW) = shadow
      par(P_MINX0) = minx0
      par(P_MINY0) = miny0
      par(P_MINZ0) = minz0
      par(P_MAXX0) = maxx0
      par(P_MAXY0) = maxy0
      par(P_MAXZ0) = maxz0
      par(P_TRACE_LEVEL) = trace_level
      par(P_ETHRESHOLD) = ethreshold
      par(P_BUFFER) = buffer
      par(P_MINDIM) = mindim
      par(P_WINDOW) = window
      par(P_MINEFFICIENCY) = minefficiency
      par(P_OUTPUT_STYLE) = output_style
      par(P_OUTPUT_DIM) = output_dim
      par(P_OUTPUT_F1_TYPE) = output_f1_type
      par(P_OUTPUT_F1_LEVEL) = output_f1_level
      par(P_OUTPUT_F1_PERIOD) = output_f1_period
      par(P_OUTPUT_F1_LB1) = output_f1_lb1
      par(P_OUTPUT_F1_LB2) = output_f1_lb2
      par(P_OUTPUT_F1_LB3) = output_f1_lb3
      par(P_OUTPUT_F1_UB1) = output_f1_ub1
      par(P_OUTPUT_F1_UB2) = output_f1_ub2
      par(P_OUTPUT_F1_UB3) = output_f1_ub3
      par(P_OUTPUT_F2_TYPE) = output_f2_type
      par(P_OUTPUT_F2_LEVEL) = output_f2_level
      par(P_OUTPUT_F2_PERIOD) = output_f2_period
      par(P_OUTPUT_F2_LB1) = output_f2_lb1
      par(P_OUTPUT_F2_LB2) = output_f2_lb2
      par(P_OUTPUT_F2_LB3) = output_f2_lb3
      par(P_OUTPUT_F2_UB1) = output_f2_ub1
      par(P_OUTPUT_F2_UB2) = output_f2_ub2
      par(P_OUTPUT_F2_UB3) = output_f2_ub3
      par(P_OUTPUT_F3_TYPE) = output_f3_type
      par(P_OUTPUT_F3_LEVEL) = output_f3_level
      par(P_OUTPUT_F3_PERIOD) = output_f3_period
      par(P_OUTPUT_F3_LB1) = output_f3_lb1
      par(P_OUTPUT_F3_LB2) = output_f3_lb2
      par(P_OUTPUT_F3_LB3) = output_f3_lb3
      par(P_OUTPUT_F3_UB1) = output_f3_ub1
      par(P_OUTPUT_F3_UB2) = output_f3_ub2
      par(P_OUTPUT_F3_UB3) = output_f3_ub3
      par(P_OUTPUT_F4_TYPE) = output_f4_type
      par(P_OUTPUT_F4_LEVEL) = output_f4_level
      par(P_OUTPUT_F4_PERIOD) = output_f4_period
      par(P_OUTPUT_F4_LB1) = output_f4_lb1
      par(P_OUTPUT_F4_LB2) = output_f4_lb2
      par(P_OUTPUT_F4_LB3) = output_f4_lb3
      par(P_OUTPUT_F4_UB1) = output_f4_ub1
      par(P_OUTPUT_F4_UB2) = output_f4_ub2
      par(P_OUTPUT_F4_UB3) = output_f4_ub3
      par(P_CLUSTERREADWRITE) = clusterreadwrite
      par(P_GHOSTWIDTH) = ghostwidth
      par(P_UPDATE_SCHEME) = update_scheme
      par(P_AMRBOUND_PREPOST) = amrbound_prepost
      par(P_AMRBOUND_TIMEALIGN) = amrbound_timealign
      par(P_ELLIPTIC_SOLVE) = elliptic_solve
      par(P_NVCYCLE) = nvcycle
      par(P_PRESWP) = preswp
      par(P_PSTSWP) = pstswp
      par(P_MAXSWEEPS) = maxsweeps
      par(P_ELL_EPSILON) = ell_epsilon
      par(P_NUM_EVOL_ITERS) = num_evol_iters
      par(P_CHKPT_PERIOD) = chkpt_period
      par(P_CHKPT_READSTATE) = chkpt_readstate
      par(P_CHKPT_CONTROL) = chkpt_control
      par(P_BOUND_WIDTH) = bound_width
      par(P_WENO_INTERP) = weno_interp
      par(P_FINDHORIZON) = findhorizon
      par(P_MASK_PERIOD) = mask_period
      par(P_MASK_USEMIN) = mask_usemin
      par(P_MASK_MINFIELD) = mask_minfield
      par(P_HORIZON_NTHETA) = horizon_ntheta
      par(P_HORIZON_NPHI) = horizon_nphi
      par(P_HORIZON_THRESH) = horizon_thresh
      par(P_HORIZON_NHOLES) = horizon_nholes
      par(P_HORIZON_RECENTERP) = horizon_recenterp
      par(P_HORIZON_GROWTH) = horizon_growth
      par(P_EMULATE_PROC) = emulate_proc
      par(P_ASSUME_SYMMETRY) = assume_symmetry
      par(P_FLUSH_PERIOD) = flush_period
      par(P_TRACERS_PERIOD) = tracers_period
      par(P_TRACERS_INITIAL) = tracers_initial
      par(P_TRACERS_SCHEME) = tracers_scheme
      par(P_VARIABLE_TIMESTEP) = variable_timestep
      par(P_CFL_LAMBDA) = cfl_lambda
      par(P_PERIODICBC) = periodicBC
      par(P_H) = hg
      par(P_DX) = hg
      par(P_DY) = hg
      par(P_DZ) = hg
      par(P_DT) = dtg
      return
      end subroutine assign_params
      
      
      
      
      
      
      subroutine set_bbox(bbox1, bbox2, bbox3, bbox4, bbox5, bbox6,
     & chr, nx, ny, nz)
      implicit none
      integer bbox1, bbox2, bbox3, bbox4, bbox5, bbox6
      integer nx, ny, nz
      real(kind=8) chr(nx,ny,nz)
      
      integer nxm, nym, nzm
      nxm = nx/2
      nym = ny/2
      nzm = nz/2
      
      if (nint(chr(1, nym, nzm)) .eq. 2) then
         bbox1 = 1
      else
         bbox1 = 0
      end if
      
      if (nint(chr(nx, nym, nzm)) .eq. 3) then
         bbox2 = 1
      else
         bbox2 = 0
      end if
      
      if (nint(chr(nxm, 1, nzm)) .eq. 4) then
         bbox3 = 1
      else
         bbox3 = 0
      end if
      
      if (nint(chr(nxm, ny, nzm)) .eq. 5) then
         bbox4 = 1
      else
         bbox4 = 0
      end if
      
      if (nint(chr(nxm, nym, 1)) .eq. 6) then
         bbox5 = 1
      else
         bbox5 = 0
      end if
      
      if (nint(chr(nxm, nym, nz)) .eq. 7) then
         bbox6 = 1
      else
         bbox6 = 0
      end if
      return
      end subroutine set_bbox
      END MODULE UTILEQS
