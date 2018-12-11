      MODULE UTILEQS
      use params
      implicit none
      CONTAINS
      subroutine assign_ptrs_init(g00,g01,g02,g03,g11,g12,g13,g22,g23,g3
     &3,K00,K01,K02,K03,K11,K12,K13,K22,K23,K33,d1g00,d1g01,d1g02,d1g03,
     &d1g11,d1g12,d1g13,d1g22,d1g23,d1g33,d2g00,d2g01,d2g02,d2g03,d2g11,
     &d2g12,d2g13,d2g22,d2g23,d2g33,d3g00,d3g01,d3g02,d3g03,d3g11,d3g12,
     &d3g13,d3g22,d3g23,d3g33,H0,H1,H2,H3,G0,d1H0,d2H0,d3H0,phim,pim,d1p
     &him,d2phim,d3phim,phin,pin,d1phin,d2phin,d3phin,phir,pir,d1phir,d2
     &phir,d3phir,phic,pic,d1phic,d2phic,d3phic,alpha,shift1,shift2,shif
     &t3,sqdetg,Tmunu00,Tmunu01,Tmunu02,Tmunu03,Tmunu11,Tmunu12,Tmunu13,
     &Tmunu22,Tmunu23,Tmunu33,order1_const,order2_const,physZ_const,rad_
     &exp,rad_speed,energy_dens,noether_dens,energy_const,momx_const,mom
     &y_const,momz_const,psi4R,psi4I,massADM,curvature,lin_momx,lin_momy
     &,ang_momz,kang_momz,R2scalar,gthethe,gthephi,gphiphi,gur,psi_el,be
     &tax_el,betay_el,betaz_el,ADM_k11,ADM_k12,ADM_k13,ADM_k22,ADM_k23,A
     &DM_k33,z4rhs0,z4rhs1,z4rhs2,z4rhs3,cctk_x,cctk_y,cctk_z,r,xphys,yp
     &hys,zphys,mask,wdiss,chr,error,flag,u0, v, w,nx,ny,nz)
      use GF
      implicit none
      integer :: nx,ny,nz
      type(gridfunction), dimension(NU) :: u0
      type(gridfunction), dimension(NV) :: v
      type(gridfunction), dimension(NW) :: w
      real(kind=8) g00(nx,ny,nz)
      real(kind=8) g01(nx,ny,nz)
      real(kind=8) g02(nx,ny,nz)
      real(kind=8) g03(nx,ny,nz)
      real(kind=8) g11(nx,ny,nz)
      real(kind=8) g12(nx,ny,nz)
      real(kind=8) g13(nx,ny,nz)
      real(kind=8) g22(nx,ny,nz)
      real(kind=8) g23(nx,ny,nz)
      real(kind=8) g33(nx,ny,nz)
      real(kind=8) K00(nx,ny,nz)
      real(kind=8) K01(nx,ny,nz)
      real(kind=8) K02(nx,ny,nz)
      real(kind=8) K03(nx,ny,nz)
      real(kind=8) K11(nx,ny,nz)
      real(kind=8) K12(nx,ny,nz)
      real(kind=8) K13(nx,ny,nz)
      real(kind=8) K22(nx,ny,nz)
      real(kind=8) K23(nx,ny,nz)
      real(kind=8) K33(nx,ny,nz)
      real(kind=8) d1g00(nx,ny,nz)
      real(kind=8) d1g01(nx,ny,nz)
      real(kind=8) d1g02(nx,ny,nz)
      real(kind=8) d1g03(nx,ny,nz)
      real(kind=8) d1g11(nx,ny,nz)
      real(kind=8) d1g12(nx,ny,nz)
      real(kind=8) d1g13(nx,ny,nz)
      real(kind=8) d1g22(nx,ny,nz)
      real(kind=8) d1g23(nx,ny,nz)
      real(kind=8) d1g33(nx,ny,nz)
      real(kind=8) d2g00(nx,ny,nz)
      real(kind=8) d2g01(nx,ny,nz)
      real(kind=8) d2g02(nx,ny,nz)
      real(kind=8) d2g03(nx,ny,nz)
      real(kind=8) d2g11(nx,ny,nz)
      real(kind=8) d2g12(nx,ny,nz)
      real(kind=8) d2g13(nx,ny,nz)
      real(kind=8) d2g22(nx,ny,nz)
      real(kind=8) d2g23(nx,ny,nz)
      real(kind=8) d2g33(nx,ny,nz)
      real(kind=8) d3g00(nx,ny,nz)
      real(kind=8) d3g01(nx,ny,nz)
      real(kind=8) d3g02(nx,ny,nz)
      real(kind=8) d3g03(nx,ny,nz)
      real(kind=8) d3g11(nx,ny,nz)
      real(kind=8) d3g12(nx,ny,nz)
      real(kind=8) d3g13(nx,ny,nz)
      real(kind=8) d3g22(nx,ny,nz)
      real(kind=8) d3g23(nx,ny,nz)
      real(kind=8) d3g33(nx,ny,nz)
      real(kind=8) H0(nx,ny,nz)
      real(kind=8) H1(nx,ny,nz)
      real(kind=8) H2(nx,ny,nz)
      real(kind=8) H3(nx,ny,nz)
      real(kind=8) G0(nx,ny,nz)
      real(kind=8) d1H0(nx,ny,nz)
      real(kind=8) d2H0(nx,ny,nz)
      real(kind=8) d3H0(nx,ny,nz)
      real(kind=8) phim(nx,ny,nz)
      real(kind=8) pim(nx,ny,nz)
      real(kind=8) d1phim(nx,ny,nz)
      real(kind=8) d2phim(nx,ny,nz)
      real(kind=8) d3phim(nx,ny,nz)
      real(kind=8) phin(nx,ny,nz)
      real(kind=8) pin(nx,ny,nz)
      real(kind=8) d1phin(nx,ny,nz)
      real(kind=8) d2phin(nx,ny,nz)
      real(kind=8) d3phin(nx,ny,nz)
      real(kind=8) phir(nx,ny,nz)
      real(kind=8) pir(nx,ny,nz)
      real(kind=8) d1phir(nx,ny,nz)
      real(kind=8) d2phir(nx,ny,nz)
      real(kind=8) d3phir(nx,ny,nz)
      real(kind=8) phic(nx,ny,nz)
      real(kind=8) pic(nx,ny,nz)
      real(kind=8) d1phic(nx,ny,nz)
      real(kind=8) d2phic(nx,ny,nz)
      real(kind=8) d3phic(nx,ny,nz)
      real(kind=8) alpha(nx,ny,nz)
      real(kind=8) shift1(nx,ny,nz)
      real(kind=8) shift2(nx,ny,nz)
      real(kind=8) shift3(nx,ny,nz)
      real(kind=8) sqdetg(nx,ny,nz)
      real(kind=8) Tmunu00(nx,ny,nz)
      real(kind=8) Tmunu01(nx,ny,nz)
      real(kind=8) Tmunu02(nx,ny,nz)
      real(kind=8) Tmunu03(nx,ny,nz)
      real(kind=8) Tmunu11(nx,ny,nz)
      real(kind=8) Tmunu12(nx,ny,nz)
      real(kind=8) Tmunu13(nx,ny,nz)
      real(kind=8) Tmunu22(nx,ny,nz)
      real(kind=8) Tmunu23(nx,ny,nz)
      real(kind=8) Tmunu33(nx,ny,nz)
      real(kind=8) order1_const(nx,ny,nz)
      real(kind=8) order2_const(nx,ny,nz)
      real(kind=8) physZ_const(nx,ny,nz)
      real(kind=8) rad_exp(nx,ny,nz)
      real(kind=8) rad_speed(nx,ny,nz)
      real(kind=8) energy_dens(nx,ny,nz)
      real(kind=8) noether_dens(nx,ny,nz)
      real(kind=8) energy_const(nx,ny,nz)
      real(kind=8) momx_const(nx,ny,nz)
      real(kind=8) momy_const(nx,ny,nz)
      real(kind=8) momz_const(nx,ny,nz)
      real(kind=8) psi4R(nx,ny,nz)
      real(kind=8) psi4I(nx,ny,nz)
      real(kind=8) massADM(nx,ny,nz)
      real(kind=8) curvature(nx,ny,nz)
      real(kind=8) lin_momx(nx,ny,nz)
      real(kind=8) lin_momy(nx,ny,nz)
      real(kind=8) ang_momz(nx,ny,nz)
      real(kind=8) kang_momz(nx,ny,nz)
      real(kind=8) R2scalar(nx,ny,nz)
      real(kind=8) gthethe(nx,ny,nz)
      real(kind=8) gthephi(nx,ny,nz)
      real(kind=8) gphiphi(nx,ny,nz)
      real(kind=8) gur(nx,ny,nz)
      real(kind=8) psi_el(nx,ny,nz)
      real(kind=8) betax_el(nx,ny,nz)
      real(kind=8) betay_el(nx,ny,nz)
      real(kind=8) betaz_el(nx,ny,nz)
      real(kind=8) ADM_k11(nx,ny,nz)
      real(kind=8) ADM_k12(nx,ny,nz)
      real(kind=8) ADM_k13(nx,ny,nz)
      real(kind=8) ADM_k22(nx,ny,nz)
      real(kind=8) ADM_k23(nx,ny,nz)
      real(kind=8) ADM_k33(nx,ny,nz)
      real(kind=8) z4rhs0(nx,ny,nz)
      real(kind=8) z4rhs1(nx,ny,nz)
      real(kind=8) z4rhs2(nx,ny,nz)
      real(kind=8) z4rhs3(nx,ny,nz)
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
      target :: g00
      target :: g01
      target :: g02
      target :: g03
      target :: g11
      target :: g12
      target :: g13
      target :: g22
      target :: g23
      target :: g33
      target :: K00
      target :: K01
      target :: K02
      target :: K03
      target :: K11
      target :: K12
      target :: K13
      target :: K22
      target :: K23
      target :: K33
      target :: d1g00
      target :: d1g01
      target :: d1g02
      target :: d1g03
      target :: d1g11
      target :: d1g12
      target :: d1g13
      target :: d1g22
      target :: d1g23
      target :: d1g33
      target :: d2g00
      target :: d2g01
      target :: d2g02
      target :: d2g03
      target :: d2g11
      target :: d2g12
      target :: d2g13
      target :: d2g22
      target :: d2g23
      target :: d2g33
      target :: d3g00
      target :: d3g01
      target :: d3g02
      target :: d3g03
      target :: d3g11
      target :: d3g12
      target :: d3g13
      target :: d3g22
      target :: d3g23
      target :: d3g33
      target :: H0
      target :: H1
      target :: H2
      target :: H3
      target :: G0
      target :: d1H0
      target :: d2H0
      target :: d3H0
      target :: phim
      target :: pim
      target :: d1phim
      target :: d2phim
      target :: d3phim
      target :: phin
      target :: pin
      target :: d1phin
      target :: d2phin
      target :: d3phin
      target :: phir
      target :: pir
      target :: d1phir
      target :: d2phir
      target :: d3phir
      target :: phic
      target :: pic
      target :: d1phic
      target :: d2phic
      target :: d3phic
      target :: alpha
      target :: shift1
      target :: shift2
      target :: shift3
      target :: sqdetg
      target :: Tmunu00
      target :: Tmunu01
      target :: Tmunu02
      target :: Tmunu03
      target :: Tmunu11
      target :: Tmunu12
      target :: Tmunu13
      target :: Tmunu22
      target :: Tmunu23
      target :: Tmunu33
      target :: order1_const
      target :: order2_const
      target :: physZ_const
      target :: rad_exp
      target :: rad_speed
      target :: energy_dens
      target :: noether_dens
      target :: energy_const
      target :: momx_const
      target :: momy_const
      target :: momz_const
      target :: psi4R
      target :: psi4I
      target :: massADM
      target :: curvature
      target :: lin_momx
      target :: lin_momy
      target :: ang_momz
      target :: kang_momz
      target :: R2scalar
      target :: gthethe
      target :: gthephi
      target :: gphiphi
      target :: gur
      target :: psi_el
      target :: betax_el
      target :: betay_el
      target :: betaz_el
      target :: ADM_k11
      target :: ADM_k12
      target :: ADM_k13
      target :: ADM_k22
      target :: ADM_k23
      target :: ADM_k33
      target :: z4rhs0
      target :: z4rhs1
      target :: z4rhs2
      target :: z4rhs3
      target :: cctk_x
      target :: cctk_y
      target :: cctk_z
      target :: r
      target :: mask
      target :: wdiss
      target :: chr
      target :: error
      target :: flag
      u0(H_G00)%d => g00
      u0(H_G01)%d => g01
      u0(H_G02)%d => g02
      u0(H_G03)%d => g03
      u0(H_G11)%d => g11
      u0(H_G12)%d => g12
      u0(H_G13)%d => g13
      u0(H_G22)%d => g22
      u0(H_G23)%d => g23
      u0(H_G33)%d => g33
      u0(H_K00)%d => K00
      u0(H_K01)%d => K01
      u0(H_K02)%d => K02
      u0(H_K03)%d => K03
      u0(H_K11)%d => K11
      u0(H_K12)%d => K12
      u0(H_K13)%d => K13
      u0(H_K22)%d => K22
      u0(H_K23)%d => K23
      u0(H_K33)%d => K33
      u0(H_D1G00)%d => d1g00
      u0(H_D1G01)%d => d1g01
      u0(H_D1G02)%d => d1g02
      u0(H_D1G03)%d => d1g03
      u0(H_D1G11)%d => d1g11
      u0(H_D1G12)%d => d1g12
      u0(H_D1G13)%d => d1g13
      u0(H_D1G22)%d => d1g22
      u0(H_D1G23)%d => d1g23
      u0(H_D1G33)%d => d1g33
      u0(H_D2G00)%d => d2g00
      u0(H_D2G01)%d => d2g01
      u0(H_D2G02)%d => d2g02
      u0(H_D2G03)%d => d2g03
      u0(H_D2G11)%d => d2g11
      u0(H_D2G12)%d => d2g12
      u0(H_D2G13)%d => d2g13
      u0(H_D2G22)%d => d2g22
      u0(H_D2G23)%d => d2g23
      u0(H_D2G33)%d => d2g33
      u0(H_D3G00)%d => d3g00
      u0(H_D3G01)%d => d3g01
      u0(H_D3G02)%d => d3g02
      u0(H_D3G03)%d => d3g03
      u0(H_D3G11)%d => d3g11
      u0(H_D3G12)%d => d3g12
      u0(H_D3G13)%d => d3g13
      u0(H_D3G22)%d => d3g22
      u0(H_D3G23)%d => d3g23
      u0(H_D3G33)%d => d3g33
      u0(H_H0)%d => H0
      u0(H_H1)%d => H1
      u0(H_H2)%d => H2
      u0(H_H3)%d => H3
      u0(H_G0)%d => G0
      u0(H_D1H0)%d => d1H0
      u0(H_D2H0)%d => d2H0
      u0(H_D3H0)%d => d3H0
      u0(H_PHIM)%d => phim
      u0(H_PIM)%d => pim
      u0(H_D1PHIM)%d => d1phim
      u0(H_D2PHIM)%d => d2phim
      u0(H_D3PHIM)%d => d3phim
      u0(H_PHIN)%d => phin
      u0(H_PIN)%d => pin
      u0(H_D1PHIN)%d => d1phin
      u0(H_D2PHIN)%d => d2phin
      u0(H_D3PHIN)%d => d3phin
      u0(H_PHIR)%d => phir
      u0(H_PIR)%d => pir
      u0(H_D1PHIR)%d => d1phir
      u0(H_D2PHIR)%d => d2phir
      u0(H_D3PHIR)%d => d3phir
      u0(H_PHIC)%d => phic
      u0(H_PIC)%d => pic
      u0(H_D1PHIC)%d => d1phic
      u0(H_D2PHIC)%d => d2phic
      u0(H_D3PHIC)%d => d3phic
      u0(H_G00)%excise_val = 0
      u0(H_G00)%dissipation = 1
      u0(H_G00)%diss_factor = 1
      u0(H_G00)%zsym = 1
      u0(H_G00)%ysym = 1
      u0(H_G00)%xsym = 1
      u0(H_G00)%stiff = 0
      u0(H_G00)%take_dx = 1
      u0(H_G00)%take_dy = 1
      u0(H_G00)%take_dz = 1
      u0(H_G01)%excise_val = 0
      u0(H_G01)%dissipation = 1
      u0(H_G01)%diss_factor = 1
      u0(H_G01)%zsym = 1
      u0(H_G01)%ysym = 1
      u0(H_G01)%xsym = -1
      u0(H_G01)%stiff = 0
      u0(H_G01)%take_dx = 1
      u0(H_G01)%take_dy = 1
      u0(H_G01)%take_dz = 1
      u0(H_G02)%excise_val = 0
      u0(H_G02)%dissipation = 1
      u0(H_G02)%diss_factor = 1
      u0(H_G02)%zsym = 1
      u0(H_G02)%ysym = -1
      u0(H_G02)%xsym = 1
      u0(H_G02)%stiff = 0
      u0(H_G02)%take_dx = 1
      u0(H_G02)%take_dy = 1
      u0(H_G02)%take_dz = 1
      u0(H_G03)%excise_val = 0
      u0(H_G03)%dissipation = 1
      u0(H_G03)%diss_factor = 1
      u0(H_G03)%zsym = -1
      u0(H_G03)%ysym = 1
      u0(H_G03)%xsym = 1
      u0(H_G03)%stiff = 0
      u0(H_G03)%take_dx = 1
      u0(H_G03)%take_dy = 1
      u0(H_G03)%take_dz = 1
      u0(H_G11)%excise_val = 0
      u0(H_G11)%dissipation = 1
      u0(H_G11)%diss_factor = 1
      u0(H_G11)%zsym = 1
      u0(H_G11)%ysym = 1
      u0(H_G11)%xsym = 1
      u0(H_G11)%stiff = 0
      u0(H_G11)%take_dx = 1
      u0(H_G11)%take_dy = 1
      u0(H_G11)%take_dz = 1
      u0(H_G12)%excise_val = 0
      u0(H_G12)%dissipation = 1
      u0(H_G12)%diss_factor = 1
      u0(H_G12)%zsym = 1
      u0(H_G12)%ysym = -1
      u0(H_G12)%xsym = -1
      u0(H_G12)%stiff = 0
      u0(H_G12)%take_dx = 1
      u0(H_G12)%take_dy = 1
      u0(H_G12)%take_dz = 1
      u0(H_G13)%excise_val = 0
      u0(H_G13)%dissipation = 1
      u0(H_G13)%diss_factor = 1
      u0(H_G13)%zsym = -1
      u0(H_G13)%ysym = 1
      u0(H_G13)%xsym = -1
      u0(H_G13)%stiff = 0
      u0(H_G13)%take_dx = 1
      u0(H_G13)%take_dy = 1
      u0(H_G13)%take_dz = 1
      u0(H_G22)%excise_val = 0
      u0(H_G22)%dissipation = 1
      u0(H_G22)%diss_factor = 1
      u0(H_G22)%zsym = 1
      u0(H_G22)%ysym = 1
      u0(H_G22)%xsym = 1
      u0(H_G22)%stiff = 0
      u0(H_G22)%take_dx = 1
      u0(H_G22)%take_dy = 1
      u0(H_G22)%take_dz = 1
      u0(H_G23)%excise_val = 0
      u0(H_G23)%dissipation = 1
      u0(H_G23)%diss_factor = 1
      u0(H_G23)%zsym = -1
      u0(H_G23)%ysym = -1
      u0(H_G23)%xsym = 1
      u0(H_G23)%stiff = 0
      u0(H_G23)%take_dx = 1
      u0(H_G23)%take_dy = 1
      u0(H_G23)%take_dz = 1
      u0(H_G33)%excise_val = 0
      u0(H_G33)%dissipation = 1
      u0(H_G33)%diss_factor = 1
      u0(H_G33)%zsym = 1
      u0(H_G33)%ysym = 1
      u0(H_G33)%xsym = 1
      u0(H_G33)%stiff = 0
      u0(H_G33)%take_dx = 1
      u0(H_G33)%take_dy = 1
      u0(H_G33)%take_dz = 1
      u0(H_K00)%excise_val = 0
      u0(H_K00)%dissipation = 1
      u0(H_K00)%diss_factor = 1
      u0(H_K00)%zsym = 1
      u0(H_K00)%ysym = 1
      u0(H_K00)%xsym = 1
      u0(H_K00)%stiff = 0
      u0(H_K00)%take_dx = 1
      u0(H_K00)%take_dy = 1
      u0(H_K00)%take_dz = 1
      u0(H_K01)%excise_val = 0
      u0(H_K01)%dissipation = 1
      u0(H_K01)%diss_factor = 1
      u0(H_K01)%zsym = 1
      u0(H_K01)%ysym = 1
      u0(H_K01)%xsym = -1
      u0(H_K01)%stiff = 0
      u0(H_K01)%take_dx = 1
      u0(H_K01)%take_dy = 1
      u0(H_K01)%take_dz = 1
      u0(H_K02)%excise_val = 0
      u0(H_K02)%dissipation = 1
      u0(H_K02)%diss_factor = 1
      u0(H_K02)%zsym = 1
      u0(H_K02)%ysym = -1
      u0(H_K02)%xsym = 1
      u0(H_K02)%stiff = 0
      u0(H_K02)%take_dx = 1
      u0(H_K02)%take_dy = 1
      u0(H_K02)%take_dz = 1
      u0(H_K03)%excise_val = 0
      u0(H_K03)%dissipation = 1
      u0(H_K03)%diss_factor = 1
      u0(H_K03)%zsym = -1
      u0(H_K03)%ysym = 1
      u0(H_K03)%xsym = 1
      u0(H_K03)%stiff = 0
      u0(H_K03)%take_dx = 1
      u0(H_K03)%take_dy = 1
      u0(H_K03)%take_dz = 1
      u0(H_K11)%excise_val = 0
      u0(H_K11)%dissipation = 1
      u0(H_K11)%diss_factor = 1
      u0(H_K11)%zsym = 1
      u0(H_K11)%ysym = 1
      u0(H_K11)%xsym = 1
      u0(H_K11)%stiff = 0
      u0(H_K11)%take_dx = 1
      u0(H_K11)%take_dy = 1
      u0(H_K11)%take_dz = 1
      u0(H_K12)%excise_val = 0
      u0(H_K12)%dissipation = 1
      u0(H_K12)%diss_factor = 1
      u0(H_K12)%zsym = 1
      u0(H_K12)%ysym = -1
      u0(H_K12)%xsym = -1
      u0(H_K12)%stiff = 0
      u0(H_K12)%take_dx = 1
      u0(H_K12)%take_dy = 1
      u0(H_K12)%take_dz = 1
      u0(H_K13)%excise_val = 0
      u0(H_K13)%dissipation = 1
      u0(H_K13)%diss_factor = 1
      u0(H_K13)%zsym = -1
      u0(H_K13)%ysym = 1
      u0(H_K13)%xsym = -1
      u0(H_K13)%stiff = 0
      u0(H_K13)%take_dx = 1
      u0(H_K13)%take_dy = 1
      u0(H_K13)%take_dz = 1
      u0(H_K22)%excise_val = 0
      u0(H_K22)%dissipation = 1
      u0(H_K22)%diss_factor = 1
      u0(H_K22)%zsym = 1
      u0(H_K22)%ysym = 1
      u0(H_K22)%xsym = 1
      u0(H_K22)%stiff = 0
      u0(H_K22)%take_dx = 1
      u0(H_K22)%take_dy = 1
      u0(H_K22)%take_dz = 1
      u0(H_K23)%excise_val = 0
      u0(H_K23)%dissipation = 1
      u0(H_K23)%diss_factor = 1
      u0(H_K23)%zsym = -1
      u0(H_K23)%ysym = -1
      u0(H_K23)%xsym = 1
      u0(H_K23)%stiff = 0
      u0(H_K23)%take_dx = 1
      u0(H_K23)%take_dy = 1
      u0(H_K23)%take_dz = 1
      u0(H_K33)%excise_val = 0
      u0(H_K33)%dissipation = 1
      u0(H_K33)%diss_factor = 1
      u0(H_K33)%zsym = 1
      u0(H_K33)%ysym = 1
      u0(H_K33)%xsym = 1
      u0(H_K33)%stiff = 0
      u0(H_K33)%take_dx = 1
      u0(H_K33)%take_dy = 1
      u0(H_K33)%take_dz = 1
      u0(H_D1G00)%excise_val = 0
      u0(H_D1G00)%dissipation = 1
      u0(H_D1G00)%diss_factor = 1
      u0(H_D1G00)%zsym = 1
      u0(H_D1G00)%ysym = 1
      u0(H_D1G00)%xsym = -1
      u0(H_D1G00)%stiff = 0
      u0(H_D1G00)%take_dx = 1
      u0(H_D1G00)%take_dy = 1
      u0(H_D1G00)%take_dz = 1
      u0(H_D1G01)%excise_val = 0
      u0(H_D1G01)%dissipation = 1
      u0(H_D1G01)%diss_factor = 1
      u0(H_D1G01)%zsym = 1
      u0(H_D1G01)%ysym = 1
      u0(H_D1G01)%xsym = 1
      u0(H_D1G01)%stiff = 0
      u0(H_D1G01)%take_dx = 1
      u0(H_D1G01)%take_dy = 1
      u0(H_D1G01)%take_dz = 1
      u0(H_D1G02)%excise_val = 0
      u0(H_D1G02)%dissipation = 1
      u0(H_D1G02)%diss_factor = 1
      u0(H_D1G02)%zsym = 1
      u0(H_D1G02)%ysym = -1
      u0(H_D1G02)%xsym = -1
      u0(H_D1G02)%stiff = 0
      u0(H_D1G02)%take_dx = 1
      u0(H_D1G02)%take_dy = 1
      u0(H_D1G02)%take_dz = 1
      u0(H_D1G03)%excise_val = 0
      u0(H_D1G03)%dissipation = 1
      u0(H_D1G03)%diss_factor = 1
      u0(H_D1G03)%zsym = -1
      u0(H_D1G03)%ysym = 1
      u0(H_D1G03)%xsym = -1
      u0(H_D1G03)%stiff = 0
      u0(H_D1G03)%take_dx = 1
      u0(H_D1G03)%take_dy = 1
      u0(H_D1G03)%take_dz = 1
      u0(H_D1G11)%excise_val = 0
      u0(H_D1G11)%dissipation = 1
      u0(H_D1G11)%diss_factor = 1
      u0(H_D1G11)%zsym = 1
      u0(H_D1G11)%ysym = 1
      u0(H_D1G11)%xsym = -1
      u0(H_D1G11)%stiff = 0
      u0(H_D1G11)%take_dx = 1
      u0(H_D1G11)%take_dy = 1
      u0(H_D1G11)%take_dz = 1
      u0(H_D1G12)%excise_val = 0
      u0(H_D1G12)%dissipation = 1
      u0(H_D1G12)%diss_factor = 1
      u0(H_D1G12)%zsym = 1
      u0(H_D1G12)%ysym = -1
      u0(H_D1G12)%xsym = 1
      u0(H_D1G12)%stiff = 0
      u0(H_D1G12)%take_dx = 1
      u0(H_D1G12)%take_dy = 1
      u0(H_D1G12)%take_dz = 1
      u0(H_D1G13)%excise_val = 0
      u0(H_D1G13)%dissipation = 1
      u0(H_D1G13)%diss_factor = 1
      u0(H_D1G13)%zsym = -1
      u0(H_D1G13)%ysym = 1
      u0(H_D1G13)%xsym = 1
      u0(H_D1G13)%stiff = 0
      u0(H_D1G13)%take_dx = 1
      u0(H_D1G13)%take_dy = 1
      u0(H_D1G13)%take_dz = 1
      u0(H_D1G22)%excise_val = 0
      u0(H_D1G22)%dissipation = 1
      u0(H_D1G22)%diss_factor = 1
      u0(H_D1G22)%zsym = 1
      u0(H_D1G22)%ysym = 1
      u0(H_D1G22)%xsym = -1
      u0(H_D1G22)%stiff = 0
      u0(H_D1G22)%take_dx = 1
      u0(H_D1G22)%take_dy = 1
      u0(H_D1G22)%take_dz = 1
      u0(H_D1G23)%excise_val = 0
      u0(H_D1G23)%dissipation = 1
      u0(H_D1G23)%diss_factor = 1
      u0(H_D1G23)%zsym = -1
      u0(H_D1G23)%ysym = -1
      u0(H_D1G23)%xsym = -1
      u0(H_D1G23)%stiff = 0
      u0(H_D1G23)%take_dx = 1
      u0(H_D1G23)%take_dy = 1
      u0(H_D1G23)%take_dz = 1
      u0(H_D1G33)%excise_val = 0
      u0(H_D1G33)%dissipation = 1
      u0(H_D1G33)%diss_factor = 1
      u0(H_D1G33)%zsym = 1
      u0(H_D1G33)%ysym = 1
      u0(H_D1G33)%xsym = -1
      u0(H_D1G33)%stiff = 0
      u0(H_D1G33)%take_dx = 1
      u0(H_D1G33)%take_dy = 1
      u0(H_D1G33)%take_dz = 1
      u0(H_D2G00)%excise_val = 0
      u0(H_D2G00)%dissipation = 1
      u0(H_D2G00)%diss_factor = 1
      u0(H_D2G00)%zsym = 1
      u0(H_D2G00)%ysym = -1
      u0(H_D2G00)%xsym = 1
      u0(H_D2G00)%stiff = 0
      u0(H_D2G00)%take_dx = 1
      u0(H_D2G00)%take_dy = 1
      u0(H_D2G00)%take_dz = 1
      u0(H_D2G01)%excise_val = 0
      u0(H_D2G01)%dissipation = 1
      u0(H_D2G01)%diss_factor = 1
      u0(H_D2G01)%zsym = 1
      u0(H_D2G01)%ysym = -1
      u0(H_D2G01)%xsym = -1
      u0(H_D2G01)%stiff = 0
      u0(H_D2G01)%take_dx = 1
      u0(H_D2G01)%take_dy = 1
      u0(H_D2G01)%take_dz = 1
      u0(H_D2G02)%excise_val = 0
      u0(H_D2G02)%dissipation = 1
      u0(H_D2G02)%diss_factor = 1
      u0(H_D2G02)%zsym = 1
      u0(H_D2G02)%ysym = 1
      u0(H_D2G02)%xsym = 1
      u0(H_D2G02)%stiff = 0
      u0(H_D2G02)%take_dx = 1
      u0(H_D2G02)%take_dy = 1
      u0(H_D2G02)%take_dz = 1
      u0(H_D2G03)%excise_val = 0
      u0(H_D2G03)%dissipation = 1
      u0(H_D2G03)%diss_factor = 1
      u0(H_D2G03)%zsym = -1
      u0(H_D2G03)%ysym = -1
      u0(H_D2G03)%xsym = 1
      u0(H_D2G03)%stiff = 0
      u0(H_D2G03)%take_dx = 1
      u0(H_D2G03)%take_dy = 1
      u0(H_D2G03)%take_dz = 1
      u0(H_D2G11)%excise_val = 0
      u0(H_D2G11)%dissipation = 1
      u0(H_D2G11)%diss_factor = 1
      u0(H_D2G11)%zsym = 1
      u0(H_D2G11)%ysym = -1
      u0(H_D2G11)%xsym = 1
      u0(H_D2G11)%stiff = 0
      u0(H_D2G11)%take_dx = 1
      u0(H_D2G11)%take_dy = 1
      u0(H_D2G11)%take_dz = 1
      u0(H_D2G12)%excise_val = 0
      u0(H_D2G12)%dissipation = 1
      u0(H_D2G12)%diss_factor = 1
      u0(H_D2G12)%zsym = 1
      u0(H_D2G12)%ysym = 1
      u0(H_D2G12)%xsym = -1
      u0(H_D2G12)%stiff = 0
      u0(H_D2G12)%take_dx = 1
      u0(H_D2G12)%take_dy = 1
      u0(H_D2G12)%take_dz = 1
      u0(H_D2G13)%excise_val = 0
      u0(H_D2G13)%dissipation = 1
      u0(H_D2G13)%diss_factor = 1
      u0(H_D2G13)%zsym = -1
      u0(H_D2G13)%ysym = -1
      u0(H_D2G13)%xsym = -1
      u0(H_D2G13)%stiff = 0
      u0(H_D2G13)%take_dx = 1
      u0(H_D2G13)%take_dy = 1
      u0(H_D2G13)%take_dz = 1
      u0(H_D2G22)%excise_val = 0
      u0(H_D2G22)%dissipation = 1
      u0(H_D2G22)%diss_factor = 1
      u0(H_D2G22)%zsym = 1
      u0(H_D2G22)%ysym = -1
      u0(H_D2G22)%xsym = 1
      u0(H_D2G22)%stiff = 0
      u0(H_D2G22)%take_dx = 1
      u0(H_D2G22)%take_dy = 1
      u0(H_D2G22)%take_dz = 1
      u0(H_D2G23)%excise_val = 0
      u0(H_D2G23)%dissipation = 1
      u0(H_D2G23)%diss_factor = 1
      u0(H_D2G23)%zsym = -1
      u0(H_D2G23)%ysym = 1
      u0(H_D2G23)%xsym = 1
      u0(H_D2G23)%stiff = 0
      u0(H_D2G23)%take_dx = 1
      u0(H_D2G23)%take_dy = 1
      u0(H_D2G23)%take_dz = 1
      u0(H_D2G33)%excise_val = 0
      u0(H_D2G33)%dissipation = 1
      u0(H_D2G33)%diss_factor = 1
      u0(H_D2G33)%zsym = 1
      u0(H_D2G33)%ysym = -1
      u0(H_D2G33)%xsym = 1
      u0(H_D2G33)%stiff = 0
      u0(H_D2G33)%take_dx = 1
      u0(H_D2G33)%take_dy = 1
      u0(H_D2G33)%take_dz = 1
      u0(H_D3G00)%excise_val = 0
      u0(H_D3G00)%dissipation = 1
      u0(H_D3G00)%diss_factor = 1
      u0(H_D3G00)%zsym = -1
      u0(H_D3G00)%ysym = 1
      u0(H_D3G00)%xsym = 1
      u0(H_D3G00)%stiff = 0
      u0(H_D3G00)%take_dx = 1
      u0(H_D3G00)%take_dy = 1
      u0(H_D3G00)%take_dz = 1
      u0(H_D3G01)%excise_val = 0
      u0(H_D3G01)%dissipation = 1
      u0(H_D3G01)%diss_factor = 1
      u0(H_D3G01)%zsym = -1
      u0(H_D3G01)%ysym = 1
      u0(H_D3G01)%xsym = -1
      u0(H_D3G01)%stiff = 0
      u0(H_D3G01)%take_dx = 1
      u0(H_D3G01)%take_dy = 1
      u0(H_D3G01)%take_dz = 1
      u0(H_D3G02)%excise_val = 0
      u0(H_D3G02)%dissipation = 1
      u0(H_D3G02)%diss_factor = 1
      u0(H_D3G02)%zsym = -1
      u0(H_D3G02)%ysym = -1
      u0(H_D3G02)%xsym = 1
      u0(H_D3G02)%stiff = 0
      u0(H_D3G02)%take_dx = 1
      u0(H_D3G02)%take_dy = 1
      u0(H_D3G02)%take_dz = 1
      u0(H_D3G03)%excise_val = 0
      u0(H_D3G03)%dissipation = 1
      u0(H_D3G03)%diss_factor = 1
      u0(H_D3G03)%zsym = 1
      u0(H_D3G03)%ysym = 1
      u0(H_D3G03)%xsym = 1
      u0(H_D3G03)%stiff = 0
      u0(H_D3G03)%take_dx = 1
      u0(H_D3G03)%take_dy = 1
      u0(H_D3G03)%take_dz = 1
      u0(H_D3G11)%excise_val = 0
      u0(H_D3G11)%dissipation = 1
      u0(H_D3G11)%diss_factor = 1
      u0(H_D3G11)%zsym = -1
      u0(H_D3G11)%ysym = 1
      u0(H_D3G11)%xsym = 1
      u0(H_D3G11)%stiff = 0
      u0(H_D3G11)%take_dx = 1
      u0(H_D3G11)%take_dy = 1
      u0(H_D3G11)%take_dz = 1
      u0(H_D3G12)%excise_val = 0
      u0(H_D3G12)%dissipation = 1
      u0(H_D3G12)%diss_factor = 1
      u0(H_D3G12)%zsym = -1
      u0(H_D3G12)%ysym = -1
      u0(H_D3G12)%xsym = -1
      u0(H_D3G12)%stiff = 0
      u0(H_D3G12)%take_dx = 1
      u0(H_D3G12)%take_dy = 1
      u0(H_D3G12)%take_dz = 1
      u0(H_D3G13)%excise_val = 0
      u0(H_D3G13)%dissipation = 1
      u0(H_D3G13)%diss_factor = 1
      u0(H_D3G13)%zsym = 1
      u0(H_D3G13)%ysym = 1
      u0(H_D3G13)%xsym = -1
      u0(H_D3G13)%stiff = 0
      u0(H_D3G13)%take_dx = 1
      u0(H_D3G13)%take_dy = 1
      u0(H_D3G13)%take_dz = 1
      u0(H_D3G22)%excise_val = 0
      u0(H_D3G22)%dissipation = 1
      u0(H_D3G22)%diss_factor = 1
      u0(H_D3G22)%zsym = -1
      u0(H_D3G22)%ysym = 1
      u0(H_D3G22)%xsym = 1
      u0(H_D3G22)%stiff = 0
      u0(H_D3G22)%take_dx = 1
      u0(H_D3G22)%take_dy = 1
      u0(H_D3G22)%take_dz = 1
      u0(H_D3G23)%excise_val = 0
      u0(H_D3G23)%dissipation = 1
      u0(H_D3G23)%diss_factor = 1
      u0(H_D3G23)%zsym = 1
      u0(H_D3G23)%ysym = -1
      u0(H_D3G23)%xsym = 1
      u0(H_D3G23)%stiff = 0
      u0(H_D3G23)%take_dx = 1
      u0(H_D3G23)%take_dy = 1
      u0(H_D3G23)%take_dz = 1
      u0(H_D3G33)%excise_val = 0
      u0(H_D3G33)%dissipation = 1
      u0(H_D3G33)%diss_factor = 1
      u0(H_D3G33)%zsym = -1
      u0(H_D3G33)%ysym = 1
      u0(H_D3G33)%xsym = 1
      u0(H_D3G33)%stiff = 0
      u0(H_D3G33)%take_dx = 1
      u0(H_D3G33)%take_dy = 1
      u0(H_D3G33)%take_dz = 1
      u0(H_H0)%excise_val = 0
      u0(H_H0)%dissipation = 1
      u0(H_H0)%diss_factor = 1
      u0(H_H0)%zsym = 1
      u0(H_H0)%ysym = 1
      u0(H_H0)%xsym = 1
      u0(H_H0)%stiff = 0
      u0(H_H0)%take_dx = 1
      u0(H_H0)%take_dy = 1
      u0(H_H0)%take_dz = 1
      u0(H_H1)%excise_val = 0
      u0(H_H1)%dissipation = 1
      u0(H_H1)%diss_factor = 1
      u0(H_H1)%zsym = 1
      u0(H_H1)%ysym = 1
      u0(H_H1)%xsym = 1
      u0(H_H1)%stiff = 0
      u0(H_H1)%take_dx = 1
      u0(H_H1)%take_dy = 1
      u0(H_H1)%take_dz = 1
      u0(H_H2)%excise_val = 0
      u0(H_H2)%dissipation = 1
      u0(H_H2)%diss_factor = 1
      u0(H_H2)%zsym = 1
      u0(H_H2)%ysym = 1
      u0(H_H2)%xsym = 1
      u0(H_H2)%stiff = 0
      u0(H_H2)%take_dx = 1
      u0(H_H2)%take_dy = 1
      u0(H_H2)%take_dz = 1
      u0(H_H3)%excise_val = 0
      u0(H_H3)%dissipation = 1
      u0(H_H3)%diss_factor = 1
      u0(H_H3)%zsym = 1
      u0(H_H3)%ysym = 1
      u0(H_H3)%xsym = 1
      u0(H_H3)%stiff = 0
      u0(H_H3)%take_dx = 1
      u0(H_H3)%take_dy = 1
      u0(H_H3)%take_dz = 1
      u0(H_G0)%excise_val = 0
      u0(H_G0)%dissipation = 1
      u0(H_G0)%diss_factor = 1
      u0(H_G0)%zsym = 1
      u0(H_G0)%ysym = 1
      u0(H_G0)%xsym = 1
      u0(H_G0)%stiff = 0
      u0(H_G0)%take_dx = 1
      u0(H_G0)%take_dy = 1
      u0(H_G0)%take_dz = 1
      u0(H_D1H0)%excise_val = 0
      u0(H_D1H0)%dissipation = 1
      u0(H_D1H0)%diss_factor = 1
      u0(H_D1H0)%zsym = 1
      u0(H_D1H0)%ysym = 1
      u0(H_D1H0)%xsym = -1
      u0(H_D1H0)%stiff = 0
      u0(H_D1H0)%take_dx = 1
      u0(H_D1H0)%take_dy = 1
      u0(H_D1H0)%take_dz = 1
      u0(H_D2H0)%excise_val = 0
      u0(H_D2H0)%dissipation = 1
      u0(H_D2H0)%diss_factor = 1
      u0(H_D2H0)%zsym = 1
      u0(H_D2H0)%ysym = -1
      u0(H_D2H0)%xsym = 1
      u0(H_D2H0)%stiff = 0
      u0(H_D2H0)%take_dx = 1
      u0(H_D2H0)%take_dy = 1
      u0(H_D2H0)%take_dz = 1
      u0(H_D3H0)%excise_val = 0
      u0(H_D3H0)%dissipation = 1
      u0(H_D3H0)%diss_factor = 1
      u0(H_D3H0)%zsym = -1
      u0(H_D3H0)%ysym = 1
      u0(H_D3H0)%xsym = 1
      u0(H_D3H0)%stiff = 0
      u0(H_D3H0)%take_dx = 1
      u0(H_D3H0)%take_dy = 1
      u0(H_D3H0)%take_dz = 1
      u0(H_PHIM)%excise_val = 0
      u0(H_PHIM)%dissipation = 1
      u0(H_PHIM)%diss_factor = 1
      u0(H_PHIM)%zsym = 1
      u0(H_PHIM)%ysym = 1
      u0(H_PHIM)%xsym = 1
      u0(H_PHIM)%stiff = 0
      u0(H_PHIM)%take_dx = 1
      u0(H_PHIM)%take_dy = 1
      u0(H_PHIM)%take_dz = 1
      u0(H_PIM)%excise_val = 0
      u0(H_PIM)%dissipation = 1
      u0(H_PIM)%diss_factor = 1
      u0(H_PIM)%zsym = 1
      u0(H_PIM)%ysym = 1
      u0(H_PIM)%xsym = 1
      u0(H_PIM)%stiff = 0
      u0(H_PIM)%take_dx = 1
      u0(H_PIM)%take_dy = 1
      u0(H_PIM)%take_dz = 1
      u0(H_D1PHIM)%excise_val = 0
      u0(H_D1PHIM)%dissipation = 1
      u0(H_D1PHIM)%diss_factor = 1
      u0(H_D1PHIM)%zsym = 1
      u0(H_D1PHIM)%ysym = 1
      u0(H_D1PHIM)%xsym = -1
      u0(H_D1PHIM)%stiff = 0
      u0(H_D1PHIM)%take_dx = 1
      u0(H_D1PHIM)%take_dy = 1
      u0(H_D1PHIM)%take_dz = 1
      u0(H_D2PHIM)%excise_val = 0
      u0(H_D2PHIM)%dissipation = 1
      u0(H_D2PHIM)%diss_factor = 1
      u0(H_D2PHIM)%zsym = 1
      u0(H_D2PHIM)%ysym = -1
      u0(H_D2PHIM)%xsym = 1
      u0(H_D2PHIM)%stiff = 0
      u0(H_D2PHIM)%take_dx = 1
      u0(H_D2PHIM)%take_dy = 1
      u0(H_D2PHIM)%take_dz = 1
      u0(H_D3PHIM)%excise_val = 0
      u0(H_D3PHIM)%dissipation = 1
      u0(H_D3PHIM)%diss_factor = 1
      u0(H_D3PHIM)%zsym = -1
      u0(H_D3PHIM)%ysym = 1
      u0(H_D3PHIM)%xsym = 1
      u0(H_D3PHIM)%stiff = 0
      u0(H_D3PHIM)%take_dx = 1
      u0(H_D3PHIM)%take_dy = 1
      u0(H_D3PHIM)%take_dz = 1
      u0(H_PHIN)%excise_val = 0
      u0(H_PHIN)%dissipation = 1
      u0(H_PHIN)%diss_factor = 1
      u0(H_PHIN)%zsym = 1
      u0(H_PHIN)%ysym = 1
      u0(H_PHIN)%xsym = 1
      u0(H_PHIN)%stiff = 0
      u0(H_PHIN)%take_dx = 1
      u0(H_PHIN)%take_dy = 1
      u0(H_PHIN)%take_dz = 1
      u0(H_PIN)%excise_val = 0
      u0(H_PIN)%dissipation = 1
      u0(H_PIN)%diss_factor = 1
      u0(H_PIN)%zsym = 1
      u0(H_PIN)%ysym = 1
      u0(H_PIN)%xsym = 1
      u0(H_PIN)%stiff = 0
      u0(H_PIN)%take_dx = 1
      u0(H_PIN)%take_dy = 1
      u0(H_PIN)%take_dz = 1
      u0(H_D1PHIN)%excise_val = 0
      u0(H_D1PHIN)%dissipation = 1
      u0(H_D1PHIN)%diss_factor = 1
      u0(H_D1PHIN)%zsym = 1
      u0(H_D1PHIN)%ysym = 1
      u0(H_D1PHIN)%xsym = -1
      u0(H_D1PHIN)%stiff = 0
      u0(H_D1PHIN)%take_dx = 1
      u0(H_D1PHIN)%take_dy = 1
      u0(H_D1PHIN)%take_dz = 1
      u0(H_D2PHIN)%excise_val = 0
      u0(H_D2PHIN)%dissipation = 1
      u0(H_D2PHIN)%diss_factor = 1
      u0(H_D2PHIN)%zsym = 1
      u0(H_D2PHIN)%ysym = -1
      u0(H_D2PHIN)%xsym = 1
      u0(H_D2PHIN)%stiff = 0
      u0(H_D2PHIN)%take_dx = 1
      u0(H_D2PHIN)%take_dy = 1
      u0(H_D2PHIN)%take_dz = 1
      u0(H_D3PHIN)%excise_val = 0
      u0(H_D3PHIN)%dissipation = 1
      u0(H_D3PHIN)%diss_factor = 1
      u0(H_D3PHIN)%zsym = -1
      u0(H_D3PHIN)%ysym = 1
      u0(H_D3PHIN)%xsym = 1
      u0(H_D3PHIN)%stiff = 0
      u0(H_D3PHIN)%take_dx = 1
      u0(H_D3PHIN)%take_dy = 1
      u0(H_D3PHIN)%take_dz = 1
      u0(H_PHIR)%excise_val = 0
      u0(H_PHIR)%dissipation = 1
      u0(H_PHIR)%diss_factor = 1
      u0(H_PHIR)%zsym = 1
      u0(H_PHIR)%ysym = 1
      u0(H_PHIR)%xsym = 1
      u0(H_PHIR)%stiff = 0
      u0(H_PHIR)%take_dx = 1
      u0(H_PHIR)%take_dy = 1
      u0(H_PHIR)%take_dz = 1
      u0(H_PIR)%excise_val = 0
      u0(H_PIR)%dissipation = 1
      u0(H_PIR)%diss_factor = 1
      u0(H_PIR)%zsym = 1
      u0(H_PIR)%ysym = 1
      u0(H_PIR)%xsym = 1
      u0(H_PIR)%stiff = 0
      u0(H_PIR)%take_dx = 1
      u0(H_PIR)%take_dy = 1
      u0(H_PIR)%take_dz = 1
      u0(H_D1PHIR)%excise_val = 0
      u0(H_D1PHIR)%dissipation = 1
      u0(H_D1PHIR)%diss_factor = 1
      u0(H_D1PHIR)%zsym = 1
      u0(H_D1PHIR)%ysym = 1
      u0(H_D1PHIR)%xsym = -1
      u0(H_D1PHIR)%stiff = 0
      u0(H_D1PHIR)%take_dx = 1
      u0(H_D1PHIR)%take_dy = 1
      u0(H_D1PHIR)%take_dz = 1
      u0(H_D2PHIR)%excise_val = 0
      u0(H_D2PHIR)%dissipation = 1
      u0(H_D2PHIR)%diss_factor = 1
      u0(H_D2PHIR)%zsym = 1
      u0(H_D2PHIR)%ysym = -1
      u0(H_D2PHIR)%xsym = 1
      u0(H_D2PHIR)%stiff = 0
      u0(H_D2PHIR)%take_dx = 1
      u0(H_D2PHIR)%take_dy = 1
      u0(H_D2PHIR)%take_dz = 1
      u0(H_D3PHIR)%excise_val = 0
      u0(H_D3PHIR)%dissipation = 1
      u0(H_D3PHIR)%diss_factor = 1
      u0(H_D3PHIR)%zsym = -1
      u0(H_D3PHIR)%ysym = 1
      u0(H_D3PHIR)%xsym = 1
      u0(H_D3PHIR)%stiff = 0
      u0(H_D3PHIR)%take_dx = 1
      u0(H_D3PHIR)%take_dy = 1
      u0(H_D3PHIR)%take_dz = 1
      u0(H_PHIC)%excise_val = 0
      u0(H_PHIC)%dissipation = 1
      u0(H_PHIC)%diss_factor = 1
      u0(H_PHIC)%zsym = 1
      u0(H_PHIC)%ysym = 1
      u0(H_PHIC)%xsym = 1
      u0(H_PHIC)%stiff = 0
      u0(H_PHIC)%take_dx = 1
      u0(H_PHIC)%take_dy = 1
      u0(H_PHIC)%take_dz = 1
      u0(H_PIC)%excise_val = 0
      u0(H_PIC)%dissipation = 1
      u0(H_PIC)%diss_factor = 1
      u0(H_PIC)%zsym = 1
      u0(H_PIC)%ysym = 1
      u0(H_PIC)%xsym = 1
      u0(H_PIC)%stiff = 0
      u0(H_PIC)%take_dx = 1
      u0(H_PIC)%take_dy = 1
      u0(H_PIC)%take_dz = 1
      u0(H_D1PHIC)%excise_val = 0
      u0(H_D1PHIC)%dissipation = 1
      u0(H_D1PHIC)%diss_factor = 1
      u0(H_D1PHIC)%zsym = 1
      u0(H_D1PHIC)%ysym = 1
      u0(H_D1PHIC)%xsym = -1
      u0(H_D1PHIC)%stiff = 0
      u0(H_D1PHIC)%take_dx = 1
      u0(H_D1PHIC)%take_dy = 1
      u0(H_D1PHIC)%take_dz = 1
      u0(H_D2PHIC)%excise_val = 0
      u0(H_D2PHIC)%dissipation = 1
      u0(H_D2PHIC)%diss_factor = 1
      u0(H_D2PHIC)%zsym = 1
      u0(H_D2PHIC)%ysym = -1
      u0(H_D2PHIC)%xsym = 1
      u0(H_D2PHIC)%stiff = 0
      u0(H_D2PHIC)%take_dx = 1
      u0(H_D2PHIC)%take_dy = 1
      u0(H_D2PHIC)%take_dz = 1
      u0(H_D3PHIC)%excise_val = 0
      u0(H_D3PHIC)%dissipation = 1
      u0(H_D3PHIC)%diss_factor = 1
      u0(H_D3PHIC)%zsym = -1
      u0(H_D3PHIC)%ysym = 1
      u0(H_D3PHIC)%xsym = 1
      u0(H_D3PHIC)%stiff = 0
      u0(H_D3PHIC)%take_dx = 1
      u0(H_D3PHIC)%take_dy = 1
      u0(H_D3PHIC)%take_dz = 1
      v(H_ALPHA)%d => alpha
      v(H_SHIFT1)%d => shift1
      v(H_SHIFT2)%d => shift2
      v(H_SHIFT3)%d => shift3
      v(H_SQDETG)%d => sqdetg
      v(H_TMUNU00)%d => Tmunu00
      v(H_TMUNU01)%d => Tmunu01
      v(H_TMUNU02)%d => Tmunu02
      v(H_TMUNU03)%d => Tmunu03
      v(H_TMUNU11)%d => Tmunu11
      v(H_TMUNU12)%d => Tmunu12
      v(H_TMUNU13)%d => Tmunu13
      v(H_TMUNU22)%d => Tmunu22
      v(H_TMUNU23)%d => Tmunu23
      v(H_TMUNU33)%d => Tmunu33
      v(H_ORDER1_CONST)%d => order1_const
      v(H_ORDER2_CONST)%d => order2_const
      v(H_PHYSZ_CONST)%d => physZ_const
      v(H_RAD_EXP)%d => rad_exp
      v(H_RAD_SPEED)%d => rad_speed
      v(H_ENERGY_DENS)%d => energy_dens
      v(H_NOETHER_DENS)%d => noether_dens
      v(H_ENERGY_CONST)%d => energy_const
      v(H_MOMX_CONST)%d => momx_const
      v(H_MOMY_CONST)%d => momy_const
      v(H_MOMZ_CONST)%d => momz_const
      v(H_PSI4R)%d => psi4R
      v(H_PSI4I)%d => psi4I
      v(H_MASSADM)%d => massADM
      v(H_CURVATURE)%d => curvature
      v(H_LIN_MOMX)%d => lin_momx
      v(H_LIN_MOMY)%d => lin_momy
      v(H_ANG_MOMZ)%d => ang_momz
      v(H_KANG_MOMZ)%d => kang_momz
      v(H_R2SCALAR)%d => R2scalar
      v(H_GTHETHE)%d => gthethe
      v(H_GTHEPHI)%d => gthephi
      v(H_GPHIPHI)%d => gphiphi
      v(H_GUR)%d => gur
      v(H_PSI_EL)%d => psi_el
      v(H_BETAX_EL)%d => betax_el
      v(H_BETAY_EL)%d => betay_el
      v(H_BETAZ_EL)%d => betaz_el
      w(H_ADM_K11)%d => ADM_k11
      w(H_ADM_K12)%d => ADM_k12
      w(H_ADM_K13)%d => ADM_k13
      w(H_ADM_K22)%d => ADM_k22
      w(H_ADM_K23)%d => ADM_k23
      w(H_ADM_K33)%d => ADM_k33
      w(H_Z4RHS0)%d => z4rhs0
      w(H_Z4RHS1)%d => z4rhs1
      w(H_Z4RHS2)%d => z4rhs2
      w(H_Z4RHS3)%d => z4rhs3
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
      subroutine assign_ptrs_fields(g00, g00_np1,g01, g01_np1,g02, g02_n
     &p1,g03, g03_np1,g11, g11_np1,g12, g12_np1,g13, g13_np1,g22, g22_np
     &1,g23, g23_np1,g33, g33_np1,K00, K00_np1,K01, K01_np1,K02, K02_np1
     &,K03, K03_np1,K11, K11_np1,K12, K12_np1,K13, K13_np1,K22, K22_np1,
     &K23, K23_np1,K33, K33_np1,d1g00, d1g00_np1,d1g01, d1g01_np1,d1g02,
     & d1g02_np1,d1g03, d1g03_np1,d1g11, d1g11_np1,d1g12, d1g12_np1,d1g1
     &3, d1g13_np1,d1g22, d1g22_np1,d1g23, d1g23_np1,d1g33, d1g33_np1,d2
     &g00, d2g00_np1,d2g01, d2g01_np1,d2g02, d2g02_np1,d2g03, d2g03_np1,
     &d2g11, d2g11_np1,d2g12, d2g12_np1,d2g13, d2g13_np1,d2g22, d2g22_np
     &1,d2g23, d2g23_np1,d2g33, d2g33_np1,d3g00, d3g00_np1,d3g01, d3g01_
     &np1,d3g02, d3g02_np1,d3g03, d3g03_np1,d3g11, d3g11_np1,d3g12, d3g1
     &2_np1,d3g13, d3g13_np1,d3g22, d3g22_np1,d3g23, d3g23_np1,d3g33, d3
     &g33_np1,H0, H0_np1,H1, H1_np1,H2, H2_np1,H3, H3_np1,G0, G0_np1,d1H
     &0, d1H0_np1,d2H0, d2H0_np1,d3H0, d3H0_np1,phim, phim_np1,pim, pim_
     &np1,d1phim, d1phim_np1,d2phim, d2phim_np1,d3phim, d3phim_np1,phin,
     & phin_np1,pin, pin_np1,d1phin, d1phin_np1,d2phin, d2phin_np1,d3phi
     &n, d3phin_np1,phir, phir_np1,pir, pir_np1,d1phir, d1phir_np1,d2phi
     &r, d2phir_np1,d3phir, d3phir_np1,phic, phic_np1,pic, pic_np1,d1phi
     &c, d1phic_np1,d2phic, d2phic_np1,d3phic, d3phic_np1,alpha,shift1,s
     &hift2,shift3,sqdetg,Tmunu00,Tmunu01,Tmunu02,Tmunu03,Tmunu11,Tmunu1
     &2,Tmunu13,Tmunu22,Tmunu23,Tmunu33,order1_const,order2_const,physZ_
     &const,rad_exp,rad_speed,energy_dens,noether_dens,energy_const,momx
     &_const,momy_const,momz_const,psi4R,psi4I,massADM,curvature,lin_mom
     &x,lin_momy,ang_momz,kang_momz,R2scalar,gthethe,gthephi,gphiphi,gur
     &,psi_el,betax_el,betay_el,betaz_el,ADM_k11,ADM_k12,ADM_k13,ADM_k22
     &,ADM_k23,ADM_k33,z4rhs0,z4rhs1,z4rhs2,z4rhs3,cctk_x,cctk_y,cctk_z,
     &r,xphys,yphys,zphys,mask,wdiss,chr,error,flag,u2, u0, v, w, nx, ny
     &, nz)
      use GF
      implicit none
      integer :: nx, ny, nz
      type(gridfunction), dimension(NU) :: u2, u0
      type(gridfunction), dimension(NW) :: w
      type(gridfunction), dimension(NV) :: v
      real(kind=8) g00(nx,ny,nz),g00_np1(nx,ny,nz)
      real(kind=8) g01(nx,ny,nz),g01_np1(nx,ny,nz)
      real(kind=8) g02(nx,ny,nz),g02_np1(nx,ny,nz)
      real(kind=8) g03(nx,ny,nz),g03_np1(nx,ny,nz)
      real(kind=8) g11(nx,ny,nz),g11_np1(nx,ny,nz)
      real(kind=8) g12(nx,ny,nz),g12_np1(nx,ny,nz)
      real(kind=8) g13(nx,ny,nz),g13_np1(nx,ny,nz)
      real(kind=8) g22(nx,ny,nz),g22_np1(nx,ny,nz)
      real(kind=8) g23(nx,ny,nz),g23_np1(nx,ny,nz)
      real(kind=8) g33(nx,ny,nz),g33_np1(nx,ny,nz)
      real(kind=8) K00(nx,ny,nz),K00_np1(nx,ny,nz)
      real(kind=8) K01(nx,ny,nz),K01_np1(nx,ny,nz)
      real(kind=8) K02(nx,ny,nz),K02_np1(nx,ny,nz)
      real(kind=8) K03(nx,ny,nz),K03_np1(nx,ny,nz)
      real(kind=8) K11(nx,ny,nz),K11_np1(nx,ny,nz)
      real(kind=8) K12(nx,ny,nz),K12_np1(nx,ny,nz)
      real(kind=8) K13(nx,ny,nz),K13_np1(nx,ny,nz)
      real(kind=8) K22(nx,ny,nz),K22_np1(nx,ny,nz)
      real(kind=8) K23(nx,ny,nz),K23_np1(nx,ny,nz)
      real(kind=8) K33(nx,ny,nz),K33_np1(nx,ny,nz)
      real(kind=8) d1g00(nx,ny,nz),d1g00_np1(nx,ny,nz)
      real(kind=8) d1g01(nx,ny,nz),d1g01_np1(nx,ny,nz)
      real(kind=8) d1g02(nx,ny,nz),d1g02_np1(nx,ny,nz)
      real(kind=8) d1g03(nx,ny,nz),d1g03_np1(nx,ny,nz)
      real(kind=8) d1g11(nx,ny,nz),d1g11_np1(nx,ny,nz)
      real(kind=8) d1g12(nx,ny,nz),d1g12_np1(nx,ny,nz)
      real(kind=8) d1g13(nx,ny,nz),d1g13_np1(nx,ny,nz)
      real(kind=8) d1g22(nx,ny,nz),d1g22_np1(nx,ny,nz)
      real(kind=8) d1g23(nx,ny,nz),d1g23_np1(nx,ny,nz)
      real(kind=8) d1g33(nx,ny,nz),d1g33_np1(nx,ny,nz)
      real(kind=8) d2g00(nx,ny,nz),d2g00_np1(nx,ny,nz)
      real(kind=8) d2g01(nx,ny,nz),d2g01_np1(nx,ny,nz)
      real(kind=8) d2g02(nx,ny,nz),d2g02_np1(nx,ny,nz)
      real(kind=8) d2g03(nx,ny,nz),d2g03_np1(nx,ny,nz)
      real(kind=8) d2g11(nx,ny,nz),d2g11_np1(nx,ny,nz)
      real(kind=8) d2g12(nx,ny,nz),d2g12_np1(nx,ny,nz)
      real(kind=8) d2g13(nx,ny,nz),d2g13_np1(nx,ny,nz)
      real(kind=8) d2g22(nx,ny,nz),d2g22_np1(nx,ny,nz)
      real(kind=8) d2g23(nx,ny,nz),d2g23_np1(nx,ny,nz)
      real(kind=8) d2g33(nx,ny,nz),d2g33_np1(nx,ny,nz)
      real(kind=8) d3g00(nx,ny,nz),d3g00_np1(nx,ny,nz)
      real(kind=8) d3g01(nx,ny,nz),d3g01_np1(nx,ny,nz)
      real(kind=8) d3g02(nx,ny,nz),d3g02_np1(nx,ny,nz)
      real(kind=8) d3g03(nx,ny,nz),d3g03_np1(nx,ny,nz)
      real(kind=8) d3g11(nx,ny,nz),d3g11_np1(nx,ny,nz)
      real(kind=8) d3g12(nx,ny,nz),d3g12_np1(nx,ny,nz)
      real(kind=8) d3g13(nx,ny,nz),d3g13_np1(nx,ny,nz)
      real(kind=8) d3g22(nx,ny,nz),d3g22_np1(nx,ny,nz)
      real(kind=8) d3g23(nx,ny,nz),d3g23_np1(nx,ny,nz)
      real(kind=8) d3g33(nx,ny,nz),d3g33_np1(nx,ny,nz)
      real(kind=8) H0(nx,ny,nz),H0_np1(nx,ny,nz)
      real(kind=8) H1(nx,ny,nz),H1_np1(nx,ny,nz)
      real(kind=8) H2(nx,ny,nz),H2_np1(nx,ny,nz)
      real(kind=8) H3(nx,ny,nz),H3_np1(nx,ny,nz)
      real(kind=8) G0(nx,ny,nz),G0_np1(nx,ny,nz)
      real(kind=8) d1H0(nx,ny,nz),d1H0_np1(nx,ny,nz)
      real(kind=8) d2H0(nx,ny,nz),d2H0_np1(nx,ny,nz)
      real(kind=8) d3H0(nx,ny,nz),d3H0_np1(nx,ny,nz)
      real(kind=8) phim(nx,ny,nz),phim_np1(nx,ny,nz)
      real(kind=8) pim(nx,ny,nz),pim_np1(nx,ny,nz)
      real(kind=8) d1phim(nx,ny,nz),d1phim_np1(nx,ny,nz)
      real(kind=8) d2phim(nx,ny,nz),d2phim_np1(nx,ny,nz)
      real(kind=8) d3phim(nx,ny,nz),d3phim_np1(nx,ny,nz)
      real(kind=8) phin(nx,ny,nz),phin_np1(nx,ny,nz)
      real(kind=8) pin(nx,ny,nz),pin_np1(nx,ny,nz)
      real(kind=8) d1phin(nx,ny,nz),d1phin_np1(nx,ny,nz)
      real(kind=8) d2phin(nx,ny,nz),d2phin_np1(nx,ny,nz)
      real(kind=8) d3phin(nx,ny,nz),d3phin_np1(nx,ny,nz)
      real(kind=8) phir(nx,ny,nz),phir_np1(nx,ny,nz)
      real(kind=8) pir(nx,ny,nz),pir_np1(nx,ny,nz)
      real(kind=8) d1phir(nx,ny,nz),d1phir_np1(nx,ny,nz)
      real(kind=8) d2phir(nx,ny,nz),d2phir_np1(nx,ny,nz)
      real(kind=8) d3phir(nx,ny,nz),d3phir_np1(nx,ny,nz)
      real(kind=8) phic(nx,ny,nz),phic_np1(nx,ny,nz)
      real(kind=8) pic(nx,ny,nz),pic_np1(nx,ny,nz)
      real(kind=8) d1phic(nx,ny,nz),d1phic_np1(nx,ny,nz)
      real(kind=8) d2phic(nx,ny,nz),d2phic_np1(nx,ny,nz)
      real(kind=8) d3phic(nx,ny,nz),d3phic_np1(nx,ny,nz)
      real(kind=8) alpha(nx,ny,nz)
      real(kind=8) shift1(nx,ny,nz)
      real(kind=8) shift2(nx,ny,nz)
      real(kind=8) shift3(nx,ny,nz)
      real(kind=8) sqdetg(nx,ny,nz)
      real(kind=8) Tmunu00(nx,ny,nz)
      real(kind=8) Tmunu01(nx,ny,nz)
      real(kind=8) Tmunu02(nx,ny,nz)
      real(kind=8) Tmunu03(nx,ny,nz)
      real(kind=8) Tmunu11(nx,ny,nz)
      real(kind=8) Tmunu12(nx,ny,nz)
      real(kind=8) Tmunu13(nx,ny,nz)
      real(kind=8) Tmunu22(nx,ny,nz)
      real(kind=8) Tmunu23(nx,ny,nz)
      real(kind=8) Tmunu33(nx,ny,nz)
      real(kind=8) order1_const(nx,ny,nz)
      real(kind=8) order2_const(nx,ny,nz)
      real(kind=8) physZ_const(nx,ny,nz)
      real(kind=8) rad_exp(nx,ny,nz)
      real(kind=8) rad_speed(nx,ny,nz)
      real(kind=8) energy_dens(nx,ny,nz)
      real(kind=8) noether_dens(nx,ny,nz)
      real(kind=8) energy_const(nx,ny,nz)
      real(kind=8) momx_const(nx,ny,nz)
      real(kind=8) momy_const(nx,ny,nz)
      real(kind=8) momz_const(nx,ny,nz)
      real(kind=8) psi4R(nx,ny,nz)
      real(kind=8) psi4I(nx,ny,nz)
      real(kind=8) massADM(nx,ny,nz)
      real(kind=8) curvature(nx,ny,nz)
      real(kind=8) lin_momx(nx,ny,nz)
      real(kind=8) lin_momy(nx,ny,nz)
      real(kind=8) ang_momz(nx,ny,nz)
      real(kind=8) kang_momz(nx,ny,nz)
      real(kind=8) R2scalar(nx,ny,nz)
      real(kind=8) gthethe(nx,ny,nz)
      real(kind=8) gthephi(nx,ny,nz)
      real(kind=8) gphiphi(nx,ny,nz)
      real(kind=8) gur(nx,ny,nz)
      real(kind=8) psi_el(nx,ny,nz)
      real(kind=8) betax_el(nx,ny,nz)
      real(kind=8) betay_el(nx,ny,nz)
      real(kind=8) betaz_el(nx,ny,nz)
      real(kind=8) ADM_k11(nx,ny,nz)
      real(kind=8) ADM_k12(nx,ny,nz)
      real(kind=8) ADM_k13(nx,ny,nz)
      real(kind=8) ADM_k22(nx,ny,nz)
      real(kind=8) ADM_k23(nx,ny,nz)
      real(kind=8) ADM_k33(nx,ny,nz)
      real(kind=8) z4rhs0(nx,ny,nz)
      real(kind=8) z4rhs1(nx,ny,nz)
      real(kind=8) z4rhs2(nx,ny,nz)
      real(kind=8) z4rhs3(nx,ny,nz)
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
      target :: g00, g00_np1
      target :: g01, g01_np1
      target :: g02, g02_np1
      target :: g03, g03_np1
      target :: g11, g11_np1
      target :: g12, g12_np1
      target :: g13, g13_np1
      target :: g22, g22_np1
      target :: g23, g23_np1
      target :: g33, g33_np1
      target :: K00, K00_np1
      target :: K01, K01_np1
      target :: K02, K02_np1
      target :: K03, K03_np1
      target :: K11, K11_np1
      target :: K12, K12_np1
      target :: K13, K13_np1
      target :: K22, K22_np1
      target :: K23, K23_np1
      target :: K33, K33_np1
      target :: d1g00, d1g00_np1
      target :: d1g01, d1g01_np1
      target :: d1g02, d1g02_np1
      target :: d1g03, d1g03_np1
      target :: d1g11, d1g11_np1
      target :: d1g12, d1g12_np1
      target :: d1g13, d1g13_np1
      target :: d1g22, d1g22_np1
      target :: d1g23, d1g23_np1
      target :: d1g33, d1g33_np1
      target :: d2g00, d2g00_np1
      target :: d2g01, d2g01_np1
      target :: d2g02, d2g02_np1
      target :: d2g03, d2g03_np1
      target :: d2g11, d2g11_np1
      target :: d2g12, d2g12_np1
      target :: d2g13, d2g13_np1
      target :: d2g22, d2g22_np1
      target :: d2g23, d2g23_np1
      target :: d2g33, d2g33_np1
      target :: d3g00, d3g00_np1
      target :: d3g01, d3g01_np1
      target :: d3g02, d3g02_np1
      target :: d3g03, d3g03_np1
      target :: d3g11, d3g11_np1
      target :: d3g12, d3g12_np1
      target :: d3g13, d3g13_np1
      target :: d3g22, d3g22_np1
      target :: d3g23, d3g23_np1
      target :: d3g33, d3g33_np1
      target :: H0, H0_np1
      target :: H1, H1_np1
      target :: H2, H2_np1
      target :: H3, H3_np1
      target :: G0, G0_np1
      target :: d1H0, d1H0_np1
      target :: d2H0, d2H0_np1
      target :: d3H0, d3H0_np1
      target :: phim, phim_np1
      target :: pim, pim_np1
      target :: d1phim, d1phim_np1
      target :: d2phim, d2phim_np1
      target :: d3phim, d3phim_np1
      target :: phin, phin_np1
      target :: pin, pin_np1
      target :: d1phin, d1phin_np1
      target :: d2phin, d2phin_np1
      target :: d3phin, d3phin_np1
      target :: phir, phir_np1
      target :: pir, pir_np1
      target :: d1phir, d1phir_np1
      target :: d2phir, d2phir_np1
      target :: d3phir, d3phir_np1
      target :: phic, phic_np1
      target :: pic, pic_np1
      target :: d1phic, d1phic_np1
      target :: d2phic, d2phic_np1
      target :: d3phic, d3phic_np1
      target :: ADM_k11
      target :: ADM_k12
      target :: ADM_k13
      target :: ADM_k22
      target :: ADM_k23
      target :: ADM_k33
      target :: z4rhs0
      target :: z4rhs1
      target :: z4rhs2
      target :: z4rhs3
      target :: cctk_x
      target :: cctk_y
      target :: cctk_z
      target :: r
      target :: mask
      target :: wdiss
      target :: chr
      target :: error
      target :: flag
      target :: alpha
      target :: shift1
      target :: shift2
      target :: shift3
      target :: sqdetg
      target :: Tmunu00
      target :: Tmunu01
      target :: Tmunu02
      target :: Tmunu03
      target :: Tmunu11
      target :: Tmunu12
      target :: Tmunu13
      target :: Tmunu22
      target :: Tmunu23
      target :: Tmunu33
      target :: order1_const
      target :: order2_const
      target :: physZ_const
      target :: rad_exp
      target :: rad_speed
      target :: energy_dens
      target :: noether_dens
      target :: energy_const
      target :: momx_const
      target :: momy_const
      target :: momz_const
      target :: psi4R
      target :: psi4I
      target :: massADM
      target :: curvature
      target :: lin_momx
      target :: lin_momy
      target :: ang_momz
      target :: kang_momz
      target :: R2scalar
      target :: gthethe
      target :: gthephi
      target :: gphiphi
      target :: gur
      target :: psi_el
      target :: betax_el
      target :: betay_el
      target :: betaz_el
      u0(H_G00)%d => g00
      u2(H_G00)%d => g00_np1
      u0(H_G01)%d => g01
      u2(H_G01)%d => g01_np1
      u0(H_G02)%d => g02
      u2(H_G02)%d => g02_np1
      u0(H_G03)%d => g03
      u2(H_G03)%d => g03_np1
      u0(H_G11)%d => g11
      u2(H_G11)%d => g11_np1
      u0(H_G12)%d => g12
      u2(H_G12)%d => g12_np1
      u0(H_G13)%d => g13
      u2(H_G13)%d => g13_np1
      u0(H_G22)%d => g22
      u2(H_G22)%d => g22_np1
      u0(H_G23)%d => g23
      u2(H_G23)%d => g23_np1
      u0(H_G33)%d => g33
      u2(H_G33)%d => g33_np1
      u0(H_K00)%d => K00
      u2(H_K00)%d => K00_np1
      u0(H_K01)%d => K01
      u2(H_K01)%d => K01_np1
      u0(H_K02)%d => K02
      u2(H_K02)%d => K02_np1
      u0(H_K03)%d => K03
      u2(H_K03)%d => K03_np1
      u0(H_K11)%d => K11
      u2(H_K11)%d => K11_np1
      u0(H_K12)%d => K12
      u2(H_K12)%d => K12_np1
      u0(H_K13)%d => K13
      u2(H_K13)%d => K13_np1
      u0(H_K22)%d => K22
      u2(H_K22)%d => K22_np1
      u0(H_K23)%d => K23
      u2(H_K23)%d => K23_np1
      u0(H_K33)%d => K33
      u2(H_K33)%d => K33_np1
      u0(H_D1G00)%d => d1g00
      u2(H_D1G00)%d => d1g00_np1
      u0(H_D1G01)%d => d1g01
      u2(H_D1G01)%d => d1g01_np1
      u0(H_D1G02)%d => d1g02
      u2(H_D1G02)%d => d1g02_np1
      u0(H_D1G03)%d => d1g03
      u2(H_D1G03)%d => d1g03_np1
      u0(H_D1G11)%d => d1g11
      u2(H_D1G11)%d => d1g11_np1
      u0(H_D1G12)%d => d1g12
      u2(H_D1G12)%d => d1g12_np1
      u0(H_D1G13)%d => d1g13
      u2(H_D1G13)%d => d1g13_np1
      u0(H_D1G22)%d => d1g22
      u2(H_D1G22)%d => d1g22_np1
      u0(H_D1G23)%d => d1g23
      u2(H_D1G23)%d => d1g23_np1
      u0(H_D1G33)%d => d1g33
      u2(H_D1G33)%d => d1g33_np1
      u0(H_D2G00)%d => d2g00
      u2(H_D2G00)%d => d2g00_np1
      u0(H_D2G01)%d => d2g01
      u2(H_D2G01)%d => d2g01_np1
      u0(H_D2G02)%d => d2g02
      u2(H_D2G02)%d => d2g02_np1
      u0(H_D2G03)%d => d2g03
      u2(H_D2G03)%d => d2g03_np1
      u0(H_D2G11)%d => d2g11
      u2(H_D2G11)%d => d2g11_np1
      u0(H_D2G12)%d => d2g12
      u2(H_D2G12)%d => d2g12_np1
      u0(H_D2G13)%d => d2g13
      u2(H_D2G13)%d => d2g13_np1
      u0(H_D2G22)%d => d2g22
      u2(H_D2G22)%d => d2g22_np1
      u0(H_D2G23)%d => d2g23
      u2(H_D2G23)%d => d2g23_np1
      u0(H_D2G33)%d => d2g33
      u2(H_D2G33)%d => d2g33_np1
      u0(H_D3G00)%d => d3g00
      u2(H_D3G00)%d => d3g00_np1
      u0(H_D3G01)%d => d3g01
      u2(H_D3G01)%d => d3g01_np1
      u0(H_D3G02)%d => d3g02
      u2(H_D3G02)%d => d3g02_np1
      u0(H_D3G03)%d => d3g03
      u2(H_D3G03)%d => d3g03_np1
      u0(H_D3G11)%d => d3g11
      u2(H_D3G11)%d => d3g11_np1
      u0(H_D3G12)%d => d3g12
      u2(H_D3G12)%d => d3g12_np1
      u0(H_D3G13)%d => d3g13
      u2(H_D3G13)%d => d3g13_np1
      u0(H_D3G22)%d => d3g22
      u2(H_D3G22)%d => d3g22_np1
      u0(H_D3G23)%d => d3g23
      u2(H_D3G23)%d => d3g23_np1
      u0(H_D3G33)%d => d3g33
      u2(H_D3G33)%d => d3g33_np1
      u0(H_H0)%d => H0
      u2(H_H0)%d => H0_np1
      u0(H_H1)%d => H1
      u2(H_H1)%d => H1_np1
      u0(H_H2)%d => H2
      u2(H_H2)%d => H2_np1
      u0(H_H3)%d => H3
      u2(H_H3)%d => H3_np1
      u0(H_G0)%d => G0
      u2(H_G0)%d => G0_np1
      u0(H_D1H0)%d => d1H0
      u2(H_D1H0)%d => d1H0_np1
      u0(H_D2H0)%d => d2H0
      u2(H_D2H0)%d => d2H0_np1
      u0(H_D3H0)%d => d3H0
      u2(H_D3H0)%d => d3H0_np1
      u0(H_PHIM)%d => phim
      u2(H_PHIM)%d => phim_np1
      u0(H_PIM)%d => pim
      u2(H_PIM)%d => pim_np1
      u0(H_D1PHIM)%d => d1phim
      u2(H_D1PHIM)%d => d1phim_np1
      u0(H_D2PHIM)%d => d2phim
      u2(H_D2PHIM)%d => d2phim_np1
      u0(H_D3PHIM)%d => d3phim
      u2(H_D3PHIM)%d => d3phim_np1
      u0(H_PHIN)%d => phin
      u2(H_PHIN)%d => phin_np1
      u0(H_PIN)%d => pin
      u2(H_PIN)%d => pin_np1
      u0(H_D1PHIN)%d => d1phin
      u2(H_D1PHIN)%d => d1phin_np1
      u0(H_D2PHIN)%d => d2phin
      u2(H_D2PHIN)%d => d2phin_np1
      u0(H_D3PHIN)%d => d3phin
      u2(H_D3PHIN)%d => d3phin_np1
      u0(H_PHIR)%d => phir
      u2(H_PHIR)%d => phir_np1
      u0(H_PIR)%d => pir
      u2(H_PIR)%d => pir_np1
      u0(H_D1PHIR)%d => d1phir
      u2(H_D1PHIR)%d => d1phir_np1
      u0(H_D2PHIR)%d => d2phir
      u2(H_D2PHIR)%d => d2phir_np1
      u0(H_D3PHIR)%d => d3phir
      u2(H_D3PHIR)%d => d3phir_np1
      u0(H_PHIC)%d => phic
      u2(H_PHIC)%d => phic_np1
      u0(H_PIC)%d => pic
      u2(H_PIC)%d => pic_np1
      u0(H_D1PHIC)%d => d1phic
      u2(H_D1PHIC)%d => d1phic_np1
      u0(H_D2PHIC)%d => d2phic
      u2(H_D2PHIC)%d => d2phic_np1
      u0(H_D3PHIC)%d => d3phic
      u2(H_D3PHIC)%d => d3phic_np1
      v(H_ALPHA)%d => alpha
      v(H_SHIFT1)%d => shift1
      v(H_SHIFT2)%d => shift2
      v(H_SHIFT3)%d => shift3
      v(H_SQDETG)%d => sqdetg
      v(H_TMUNU00)%d => Tmunu00
      v(H_TMUNU01)%d => Tmunu01
      v(H_TMUNU02)%d => Tmunu02
      v(H_TMUNU03)%d => Tmunu03
      v(H_TMUNU11)%d => Tmunu11
      v(H_TMUNU12)%d => Tmunu12
      v(H_TMUNU13)%d => Tmunu13
      v(H_TMUNU22)%d => Tmunu22
      v(H_TMUNU23)%d => Tmunu23
      v(H_TMUNU33)%d => Tmunu33
      v(H_ORDER1_CONST)%d => order1_const
      v(H_ORDER2_CONST)%d => order2_const
      v(H_PHYSZ_CONST)%d => physZ_const
      v(H_RAD_EXP)%d => rad_exp
      v(H_RAD_SPEED)%d => rad_speed
      v(H_ENERGY_DENS)%d => energy_dens
      v(H_NOETHER_DENS)%d => noether_dens
      v(H_ENERGY_CONST)%d => energy_const
      v(H_MOMX_CONST)%d => momx_const
      v(H_MOMY_CONST)%d => momy_const
      v(H_MOMZ_CONST)%d => momz_const
      v(H_PSI4R)%d => psi4R
      v(H_PSI4I)%d => psi4I
      v(H_MASSADM)%d => massADM
      v(H_CURVATURE)%d => curvature
      v(H_LIN_MOMX)%d => lin_momx
      v(H_LIN_MOMY)%d => lin_momy
      v(H_ANG_MOMZ)%d => ang_momz
      v(H_KANG_MOMZ)%d => kang_momz
      v(H_R2SCALAR)%d => R2scalar
      v(H_GTHETHE)%d => gthethe
      v(H_GTHEPHI)%d => gthephi
      v(H_GPHIPHI)%d => gphiphi
      v(H_GUR)%d => gur
      v(H_PSI_EL)%d => psi_el
      v(H_BETAX_EL)%d => betax_el
      v(H_BETAY_EL)%d => betay_el
      v(H_BETAZ_EL)%d => betaz_el
      w(H_ADM_K11)%d => ADM_k11
      w(H_ADM_K12)%d => ADM_k12
      w(H_ADM_K13)%d => ADM_k13
      w(H_ADM_K22)%d => ADM_k22
      w(H_ADM_K23)%d => ADM_k23
      w(H_ADM_K33)%d => ADM_k33
      w(H_Z4RHS0)%d => z4rhs0
      w(H_Z4RHS1)%d => z4rhs1
      w(H_Z4RHS2)%d => z4rhs2
      w(H_Z4RHS3)%d => z4rhs3
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
      u0(H_G00)%excise_val = 0
      u2(H_G00)%excise_val = 0
      u0(H_G00)%dissipation = 1
      u2(H_G00)%dissipation = 1
      u0(H_G00)%diss_factor = 1
      u2(H_G00)%diss_factor = 1
      u0(H_G00)%zsym = 1
      u0(H_G00)%ysym = 1
      u0(H_G00)%xsym = 1
      u2(H_G00)%zsym = 1
      u2(H_G00)%ysym = 1
      u2(H_G00)%xsym = 1
      u0(H_G00)%stiff = 0
      u2(H_G00)%stiff = 0
      u2(H_G00)%take_dx = 1
      u0(H_G00)%take_dx = 1
      u2(H_G00)%take_dy = 1
      u0(H_G00)%take_dy = 1
      u2(H_G00)%take_dz = 1
      u0(H_G00)%take_dz = 1
      u0(H_G01)%excise_val = 0
      u2(H_G01)%excise_val = 0
      u0(H_G01)%dissipation = 1
      u2(H_G01)%dissipation = 1
      u0(H_G01)%diss_factor = 1
      u2(H_G01)%diss_factor = 1
      u0(H_G01)%zsym = 1
      u0(H_G01)%ysym = 1
      u0(H_G01)%xsym = -1
      u2(H_G01)%zsym = 1
      u2(H_G01)%ysym = 1
      u2(H_G01)%xsym = -1
      u0(H_G01)%stiff = 0
      u2(H_G01)%stiff = 0
      u2(H_G01)%take_dx = 1
      u0(H_G01)%take_dx = 1
      u2(H_G01)%take_dy = 1
      u0(H_G01)%take_dy = 1
      u2(H_G01)%take_dz = 1
      u0(H_G01)%take_dz = 1
      u0(H_G02)%excise_val = 0
      u2(H_G02)%excise_val = 0
      u0(H_G02)%dissipation = 1
      u2(H_G02)%dissipation = 1
      u0(H_G02)%diss_factor = 1
      u2(H_G02)%diss_factor = 1
      u0(H_G02)%zsym = 1
      u0(H_G02)%ysym = -1
      u0(H_G02)%xsym = 1
      u2(H_G02)%zsym = 1
      u2(H_G02)%ysym = -1
      u2(H_G02)%xsym = 1
      u0(H_G02)%stiff = 0
      u2(H_G02)%stiff = 0
      u2(H_G02)%take_dx = 1
      u0(H_G02)%take_dx = 1
      u2(H_G02)%take_dy = 1
      u0(H_G02)%take_dy = 1
      u2(H_G02)%take_dz = 1
      u0(H_G02)%take_dz = 1
      u0(H_G03)%excise_val = 0
      u2(H_G03)%excise_val = 0
      u0(H_G03)%dissipation = 1
      u2(H_G03)%dissipation = 1
      u0(H_G03)%diss_factor = 1
      u2(H_G03)%diss_factor = 1
      u0(H_G03)%zsym = -1
      u0(H_G03)%ysym = 1
      u0(H_G03)%xsym = 1
      u2(H_G03)%zsym = -1
      u2(H_G03)%ysym = 1
      u2(H_G03)%xsym = 1
      u0(H_G03)%stiff = 0
      u2(H_G03)%stiff = 0
      u2(H_G03)%take_dx = 1
      u0(H_G03)%take_dx = 1
      u2(H_G03)%take_dy = 1
      u0(H_G03)%take_dy = 1
      u2(H_G03)%take_dz = 1
      u0(H_G03)%take_dz = 1
      u0(H_G11)%excise_val = 0
      u2(H_G11)%excise_val = 0
      u0(H_G11)%dissipation = 1
      u2(H_G11)%dissipation = 1
      u0(H_G11)%diss_factor = 1
      u2(H_G11)%diss_factor = 1
      u0(H_G11)%zsym = 1
      u0(H_G11)%ysym = 1
      u0(H_G11)%xsym = 1
      u2(H_G11)%zsym = 1
      u2(H_G11)%ysym = 1
      u2(H_G11)%xsym = 1
      u0(H_G11)%stiff = 0
      u2(H_G11)%stiff = 0
      u2(H_G11)%take_dx = 1
      u0(H_G11)%take_dx = 1
      u2(H_G11)%take_dy = 1
      u0(H_G11)%take_dy = 1
      u2(H_G11)%take_dz = 1
      u0(H_G11)%take_dz = 1
      u0(H_G12)%excise_val = 0
      u2(H_G12)%excise_val = 0
      u0(H_G12)%dissipation = 1
      u2(H_G12)%dissipation = 1
      u0(H_G12)%diss_factor = 1
      u2(H_G12)%diss_factor = 1
      u0(H_G12)%zsym = 1
      u0(H_G12)%ysym = -1
      u0(H_G12)%xsym = -1
      u2(H_G12)%zsym = 1
      u2(H_G12)%ysym = -1
      u2(H_G12)%xsym = -1
      u0(H_G12)%stiff = 0
      u2(H_G12)%stiff = 0
      u2(H_G12)%take_dx = 1
      u0(H_G12)%take_dx = 1
      u2(H_G12)%take_dy = 1
      u0(H_G12)%take_dy = 1
      u2(H_G12)%take_dz = 1
      u0(H_G12)%take_dz = 1
      u0(H_G13)%excise_val = 0
      u2(H_G13)%excise_val = 0
      u0(H_G13)%dissipation = 1
      u2(H_G13)%dissipation = 1
      u0(H_G13)%diss_factor = 1
      u2(H_G13)%diss_factor = 1
      u0(H_G13)%zsym = -1
      u0(H_G13)%ysym = 1
      u0(H_G13)%xsym = -1
      u2(H_G13)%zsym = -1
      u2(H_G13)%ysym = 1
      u2(H_G13)%xsym = -1
      u0(H_G13)%stiff = 0
      u2(H_G13)%stiff = 0
      u2(H_G13)%take_dx = 1
      u0(H_G13)%take_dx = 1
      u2(H_G13)%take_dy = 1
      u0(H_G13)%take_dy = 1
      u2(H_G13)%take_dz = 1
      u0(H_G13)%take_dz = 1
      u0(H_G22)%excise_val = 0
      u2(H_G22)%excise_val = 0
      u0(H_G22)%dissipation = 1
      u2(H_G22)%dissipation = 1
      u0(H_G22)%diss_factor = 1
      u2(H_G22)%diss_factor = 1
      u0(H_G22)%zsym = 1
      u0(H_G22)%ysym = 1
      u0(H_G22)%xsym = 1
      u2(H_G22)%zsym = 1
      u2(H_G22)%ysym = 1
      u2(H_G22)%xsym = 1
      u0(H_G22)%stiff = 0
      u2(H_G22)%stiff = 0
      u2(H_G22)%take_dx = 1
      u0(H_G22)%take_dx = 1
      u2(H_G22)%take_dy = 1
      u0(H_G22)%take_dy = 1
      u2(H_G22)%take_dz = 1
      u0(H_G22)%take_dz = 1
      u0(H_G23)%excise_val = 0
      u2(H_G23)%excise_val = 0
      u0(H_G23)%dissipation = 1
      u2(H_G23)%dissipation = 1
      u0(H_G23)%diss_factor = 1
      u2(H_G23)%diss_factor = 1
      u0(H_G23)%zsym = -1
      u0(H_G23)%ysym = -1
      u0(H_G23)%xsym = 1
      u2(H_G23)%zsym = -1
      u2(H_G23)%ysym = -1
      u2(H_G23)%xsym = 1
      u0(H_G23)%stiff = 0
      u2(H_G23)%stiff = 0
      u2(H_G23)%take_dx = 1
      u0(H_G23)%take_dx = 1
      u2(H_G23)%take_dy = 1
      u0(H_G23)%take_dy = 1
      u2(H_G23)%take_dz = 1
      u0(H_G23)%take_dz = 1
      u0(H_G33)%excise_val = 0
      u2(H_G33)%excise_val = 0
      u0(H_G33)%dissipation = 1
      u2(H_G33)%dissipation = 1
      u0(H_G33)%diss_factor = 1
      u2(H_G33)%diss_factor = 1
      u0(H_G33)%zsym = 1
      u0(H_G33)%ysym = 1
      u0(H_G33)%xsym = 1
      u2(H_G33)%zsym = 1
      u2(H_G33)%ysym = 1
      u2(H_G33)%xsym = 1
      u0(H_G33)%stiff = 0
      u2(H_G33)%stiff = 0
      u2(H_G33)%take_dx = 1
      u0(H_G33)%take_dx = 1
      u2(H_G33)%take_dy = 1
      u0(H_G33)%take_dy = 1
      u2(H_G33)%take_dz = 1
      u0(H_G33)%take_dz = 1
      u0(H_K00)%excise_val = 0
      u2(H_K00)%excise_val = 0
      u0(H_K00)%dissipation = 1
      u2(H_K00)%dissipation = 1
      u0(H_K00)%diss_factor = 1
      u2(H_K00)%diss_factor = 1
      u0(H_K00)%zsym = 1
      u0(H_K00)%ysym = 1
      u0(H_K00)%xsym = 1
      u2(H_K00)%zsym = 1
      u2(H_K00)%ysym = 1
      u2(H_K00)%xsym = 1
      u0(H_K00)%stiff = 0
      u2(H_K00)%stiff = 0
      u2(H_K00)%take_dx = 1
      u0(H_K00)%take_dx = 1
      u2(H_K00)%take_dy = 1
      u0(H_K00)%take_dy = 1
      u2(H_K00)%take_dz = 1
      u0(H_K00)%take_dz = 1
      u0(H_K01)%excise_val = 0
      u2(H_K01)%excise_val = 0
      u0(H_K01)%dissipation = 1
      u2(H_K01)%dissipation = 1
      u0(H_K01)%diss_factor = 1
      u2(H_K01)%diss_factor = 1
      u0(H_K01)%zsym = 1
      u0(H_K01)%ysym = 1
      u0(H_K01)%xsym = -1
      u2(H_K01)%zsym = 1
      u2(H_K01)%ysym = 1
      u2(H_K01)%xsym = -1
      u0(H_K01)%stiff = 0
      u2(H_K01)%stiff = 0
      u2(H_K01)%take_dx = 1
      u0(H_K01)%take_dx = 1
      u2(H_K01)%take_dy = 1
      u0(H_K01)%take_dy = 1
      u2(H_K01)%take_dz = 1
      u0(H_K01)%take_dz = 1
      u0(H_K02)%excise_val = 0
      u2(H_K02)%excise_val = 0
      u0(H_K02)%dissipation = 1
      u2(H_K02)%dissipation = 1
      u0(H_K02)%diss_factor = 1
      u2(H_K02)%diss_factor = 1
      u0(H_K02)%zsym = 1
      u0(H_K02)%ysym = -1
      u0(H_K02)%xsym = 1
      u2(H_K02)%zsym = 1
      u2(H_K02)%ysym = -1
      u2(H_K02)%xsym = 1
      u0(H_K02)%stiff = 0
      u2(H_K02)%stiff = 0
      u2(H_K02)%take_dx = 1
      u0(H_K02)%take_dx = 1
      u2(H_K02)%take_dy = 1
      u0(H_K02)%take_dy = 1
      u2(H_K02)%take_dz = 1
      u0(H_K02)%take_dz = 1
      u0(H_K03)%excise_val = 0
      u2(H_K03)%excise_val = 0
      u0(H_K03)%dissipation = 1
      u2(H_K03)%dissipation = 1
      u0(H_K03)%diss_factor = 1
      u2(H_K03)%diss_factor = 1
      u0(H_K03)%zsym = -1
      u0(H_K03)%ysym = 1
      u0(H_K03)%xsym = 1
      u2(H_K03)%zsym = -1
      u2(H_K03)%ysym = 1
      u2(H_K03)%xsym = 1
      u0(H_K03)%stiff = 0
      u2(H_K03)%stiff = 0
      u2(H_K03)%take_dx = 1
      u0(H_K03)%take_dx = 1
      u2(H_K03)%take_dy = 1
      u0(H_K03)%take_dy = 1
      u2(H_K03)%take_dz = 1
      u0(H_K03)%take_dz = 1
      u0(H_K11)%excise_val = 0
      u2(H_K11)%excise_val = 0
      u0(H_K11)%dissipation = 1
      u2(H_K11)%dissipation = 1
      u0(H_K11)%diss_factor = 1
      u2(H_K11)%diss_factor = 1
      u0(H_K11)%zsym = 1
      u0(H_K11)%ysym = 1
      u0(H_K11)%xsym = 1
      u2(H_K11)%zsym = 1
      u2(H_K11)%ysym = 1
      u2(H_K11)%xsym = 1
      u0(H_K11)%stiff = 0
      u2(H_K11)%stiff = 0
      u2(H_K11)%take_dx = 1
      u0(H_K11)%take_dx = 1
      u2(H_K11)%take_dy = 1
      u0(H_K11)%take_dy = 1
      u2(H_K11)%take_dz = 1
      u0(H_K11)%take_dz = 1
      u0(H_K12)%excise_val = 0
      u2(H_K12)%excise_val = 0
      u0(H_K12)%dissipation = 1
      u2(H_K12)%dissipation = 1
      u0(H_K12)%diss_factor = 1
      u2(H_K12)%diss_factor = 1
      u0(H_K12)%zsym = 1
      u0(H_K12)%ysym = -1
      u0(H_K12)%xsym = -1
      u2(H_K12)%zsym = 1
      u2(H_K12)%ysym = -1
      u2(H_K12)%xsym = -1
      u0(H_K12)%stiff = 0
      u2(H_K12)%stiff = 0
      u2(H_K12)%take_dx = 1
      u0(H_K12)%take_dx = 1
      u2(H_K12)%take_dy = 1
      u0(H_K12)%take_dy = 1
      u2(H_K12)%take_dz = 1
      u0(H_K12)%take_dz = 1
      u0(H_K13)%excise_val = 0
      u2(H_K13)%excise_val = 0
      u0(H_K13)%dissipation = 1
      u2(H_K13)%dissipation = 1
      u0(H_K13)%diss_factor = 1
      u2(H_K13)%diss_factor = 1
      u0(H_K13)%zsym = -1
      u0(H_K13)%ysym = 1
      u0(H_K13)%xsym = -1
      u2(H_K13)%zsym = -1
      u2(H_K13)%ysym = 1
      u2(H_K13)%xsym = -1
      u0(H_K13)%stiff = 0
      u2(H_K13)%stiff = 0
      u2(H_K13)%take_dx = 1
      u0(H_K13)%take_dx = 1
      u2(H_K13)%take_dy = 1
      u0(H_K13)%take_dy = 1
      u2(H_K13)%take_dz = 1
      u0(H_K13)%take_dz = 1
      u0(H_K22)%excise_val = 0
      u2(H_K22)%excise_val = 0
      u0(H_K22)%dissipation = 1
      u2(H_K22)%dissipation = 1
      u0(H_K22)%diss_factor = 1
      u2(H_K22)%diss_factor = 1
      u0(H_K22)%zsym = 1
      u0(H_K22)%ysym = 1
      u0(H_K22)%xsym = 1
      u2(H_K22)%zsym = 1
      u2(H_K22)%ysym = 1
      u2(H_K22)%xsym = 1
      u0(H_K22)%stiff = 0
      u2(H_K22)%stiff = 0
      u2(H_K22)%take_dx = 1
      u0(H_K22)%take_dx = 1
      u2(H_K22)%take_dy = 1
      u0(H_K22)%take_dy = 1
      u2(H_K22)%take_dz = 1
      u0(H_K22)%take_dz = 1
      u0(H_K23)%excise_val = 0
      u2(H_K23)%excise_val = 0
      u0(H_K23)%dissipation = 1
      u2(H_K23)%dissipation = 1
      u0(H_K23)%diss_factor = 1
      u2(H_K23)%diss_factor = 1
      u0(H_K23)%zsym = -1
      u0(H_K23)%ysym = -1
      u0(H_K23)%xsym = 1
      u2(H_K23)%zsym = -1
      u2(H_K23)%ysym = -1
      u2(H_K23)%xsym = 1
      u0(H_K23)%stiff = 0
      u2(H_K23)%stiff = 0
      u2(H_K23)%take_dx = 1
      u0(H_K23)%take_dx = 1
      u2(H_K23)%take_dy = 1
      u0(H_K23)%take_dy = 1
      u2(H_K23)%take_dz = 1
      u0(H_K23)%take_dz = 1
      u0(H_K33)%excise_val = 0
      u2(H_K33)%excise_val = 0
      u0(H_K33)%dissipation = 1
      u2(H_K33)%dissipation = 1
      u0(H_K33)%diss_factor = 1
      u2(H_K33)%diss_factor = 1
      u0(H_K33)%zsym = 1
      u0(H_K33)%ysym = 1
      u0(H_K33)%xsym = 1
      u2(H_K33)%zsym = 1
      u2(H_K33)%ysym = 1
      u2(H_K33)%xsym = 1
      u0(H_K33)%stiff = 0
      u2(H_K33)%stiff = 0
      u2(H_K33)%take_dx = 1
      u0(H_K33)%take_dx = 1
      u2(H_K33)%take_dy = 1
      u0(H_K33)%take_dy = 1
      u2(H_K33)%take_dz = 1
      u0(H_K33)%take_dz = 1
      u0(H_D1G00)%excise_val = 0
      u2(H_D1G00)%excise_val = 0
      u0(H_D1G00)%dissipation = 1
      u2(H_D1G00)%dissipation = 1
      u0(H_D1G00)%diss_factor = 1
      u2(H_D1G00)%diss_factor = 1
      u0(H_D1G00)%zsym = 1
      u0(H_D1G00)%ysym = 1
      u0(H_D1G00)%xsym = -1
      u2(H_D1G00)%zsym = 1
      u2(H_D1G00)%ysym = 1
      u2(H_D1G00)%xsym = -1
      u0(H_D1G00)%stiff = 0
      u2(H_D1G00)%stiff = 0
      u2(H_D1G00)%take_dx = 1
      u0(H_D1G00)%take_dx = 1
      u2(H_D1G00)%take_dy = 1
      u0(H_D1G00)%take_dy = 1
      u2(H_D1G00)%take_dz = 1
      u0(H_D1G00)%take_dz = 1
      u0(H_D1G01)%excise_val = 0
      u2(H_D1G01)%excise_val = 0
      u0(H_D1G01)%dissipation = 1
      u2(H_D1G01)%dissipation = 1
      u0(H_D1G01)%diss_factor = 1
      u2(H_D1G01)%diss_factor = 1
      u0(H_D1G01)%zsym = 1
      u0(H_D1G01)%ysym = 1
      u0(H_D1G01)%xsym = 1
      u2(H_D1G01)%zsym = 1
      u2(H_D1G01)%ysym = 1
      u2(H_D1G01)%xsym = 1
      u0(H_D1G01)%stiff = 0
      u2(H_D1G01)%stiff = 0
      u2(H_D1G01)%take_dx = 1
      u0(H_D1G01)%take_dx = 1
      u2(H_D1G01)%take_dy = 1
      u0(H_D1G01)%take_dy = 1
      u2(H_D1G01)%take_dz = 1
      u0(H_D1G01)%take_dz = 1
      u0(H_D1G02)%excise_val = 0
      u2(H_D1G02)%excise_val = 0
      u0(H_D1G02)%dissipation = 1
      u2(H_D1G02)%dissipation = 1
      u0(H_D1G02)%diss_factor = 1
      u2(H_D1G02)%diss_factor = 1
      u0(H_D1G02)%zsym = 1
      u0(H_D1G02)%ysym = -1
      u0(H_D1G02)%xsym = -1
      u2(H_D1G02)%zsym = 1
      u2(H_D1G02)%ysym = -1
      u2(H_D1G02)%xsym = -1
      u0(H_D1G02)%stiff = 0
      u2(H_D1G02)%stiff = 0
      u2(H_D1G02)%take_dx = 1
      u0(H_D1G02)%take_dx = 1
      u2(H_D1G02)%take_dy = 1
      u0(H_D1G02)%take_dy = 1
      u2(H_D1G02)%take_dz = 1
      u0(H_D1G02)%take_dz = 1
      u0(H_D1G03)%excise_val = 0
      u2(H_D1G03)%excise_val = 0
      u0(H_D1G03)%dissipation = 1
      u2(H_D1G03)%dissipation = 1
      u0(H_D1G03)%diss_factor = 1
      u2(H_D1G03)%diss_factor = 1
      u0(H_D1G03)%zsym = -1
      u0(H_D1G03)%ysym = 1
      u0(H_D1G03)%xsym = -1
      u2(H_D1G03)%zsym = -1
      u2(H_D1G03)%ysym = 1
      u2(H_D1G03)%xsym = -1
      u0(H_D1G03)%stiff = 0
      u2(H_D1G03)%stiff = 0
      u2(H_D1G03)%take_dx = 1
      u0(H_D1G03)%take_dx = 1
      u2(H_D1G03)%take_dy = 1
      u0(H_D1G03)%take_dy = 1
      u2(H_D1G03)%take_dz = 1
      u0(H_D1G03)%take_dz = 1
      u0(H_D1G11)%excise_val = 0
      u2(H_D1G11)%excise_val = 0
      u0(H_D1G11)%dissipation = 1
      u2(H_D1G11)%dissipation = 1
      u0(H_D1G11)%diss_factor = 1
      u2(H_D1G11)%diss_factor = 1
      u0(H_D1G11)%zsym = 1
      u0(H_D1G11)%ysym = 1
      u0(H_D1G11)%xsym = -1
      u2(H_D1G11)%zsym = 1
      u2(H_D1G11)%ysym = 1
      u2(H_D1G11)%xsym = -1
      u0(H_D1G11)%stiff = 0
      u2(H_D1G11)%stiff = 0
      u2(H_D1G11)%take_dx = 1
      u0(H_D1G11)%take_dx = 1
      u2(H_D1G11)%take_dy = 1
      u0(H_D1G11)%take_dy = 1
      u2(H_D1G11)%take_dz = 1
      u0(H_D1G11)%take_dz = 1
      u0(H_D1G12)%excise_val = 0
      u2(H_D1G12)%excise_val = 0
      u0(H_D1G12)%dissipation = 1
      u2(H_D1G12)%dissipation = 1
      u0(H_D1G12)%diss_factor = 1
      u2(H_D1G12)%diss_factor = 1
      u0(H_D1G12)%zsym = 1
      u0(H_D1G12)%ysym = -1
      u0(H_D1G12)%xsym = 1
      u2(H_D1G12)%zsym = 1
      u2(H_D1G12)%ysym = -1
      u2(H_D1G12)%xsym = 1
      u0(H_D1G12)%stiff = 0
      u2(H_D1G12)%stiff = 0
      u2(H_D1G12)%take_dx = 1
      u0(H_D1G12)%take_dx = 1
      u2(H_D1G12)%take_dy = 1
      u0(H_D1G12)%take_dy = 1
      u2(H_D1G12)%take_dz = 1
      u0(H_D1G12)%take_dz = 1
      u0(H_D1G13)%excise_val = 0
      u2(H_D1G13)%excise_val = 0
      u0(H_D1G13)%dissipation = 1
      u2(H_D1G13)%dissipation = 1
      u0(H_D1G13)%diss_factor = 1
      u2(H_D1G13)%diss_factor = 1
      u0(H_D1G13)%zsym = -1
      u0(H_D1G13)%ysym = 1
      u0(H_D1G13)%xsym = 1
      u2(H_D1G13)%zsym = -1
      u2(H_D1G13)%ysym = 1
      u2(H_D1G13)%xsym = 1
      u0(H_D1G13)%stiff = 0
      u2(H_D1G13)%stiff = 0
      u2(H_D1G13)%take_dx = 1
      u0(H_D1G13)%take_dx = 1
      u2(H_D1G13)%take_dy = 1
      u0(H_D1G13)%take_dy = 1
      u2(H_D1G13)%take_dz = 1
      u0(H_D1G13)%take_dz = 1
      u0(H_D1G22)%excise_val = 0
      u2(H_D1G22)%excise_val = 0
      u0(H_D1G22)%dissipation = 1
      u2(H_D1G22)%dissipation = 1
      u0(H_D1G22)%diss_factor = 1
      u2(H_D1G22)%diss_factor = 1
      u0(H_D1G22)%zsym = 1
      u0(H_D1G22)%ysym = 1
      u0(H_D1G22)%xsym = -1
      u2(H_D1G22)%zsym = 1
      u2(H_D1G22)%ysym = 1
      u2(H_D1G22)%xsym = -1
      u0(H_D1G22)%stiff = 0
      u2(H_D1G22)%stiff = 0
      u2(H_D1G22)%take_dx = 1
      u0(H_D1G22)%take_dx = 1
      u2(H_D1G22)%take_dy = 1
      u0(H_D1G22)%take_dy = 1
      u2(H_D1G22)%take_dz = 1
      u0(H_D1G22)%take_dz = 1
      u0(H_D1G23)%excise_val = 0
      u2(H_D1G23)%excise_val = 0
      u0(H_D1G23)%dissipation = 1
      u2(H_D1G23)%dissipation = 1
      u0(H_D1G23)%diss_factor = 1
      u2(H_D1G23)%diss_factor = 1
      u0(H_D1G23)%zsym = -1
      u0(H_D1G23)%ysym = -1
      u0(H_D1G23)%xsym = -1
      u2(H_D1G23)%zsym = -1
      u2(H_D1G23)%ysym = -1
      u2(H_D1G23)%xsym = -1
      u0(H_D1G23)%stiff = 0
      u2(H_D1G23)%stiff = 0
      u2(H_D1G23)%take_dx = 1
      u0(H_D1G23)%take_dx = 1
      u2(H_D1G23)%take_dy = 1
      u0(H_D1G23)%take_dy = 1
      u2(H_D1G23)%take_dz = 1
      u0(H_D1G23)%take_dz = 1
      u0(H_D1G33)%excise_val = 0
      u2(H_D1G33)%excise_val = 0
      u0(H_D1G33)%dissipation = 1
      u2(H_D1G33)%dissipation = 1
      u0(H_D1G33)%diss_factor = 1
      u2(H_D1G33)%diss_factor = 1
      u0(H_D1G33)%zsym = 1
      u0(H_D1G33)%ysym = 1
      u0(H_D1G33)%xsym = -1
      u2(H_D1G33)%zsym = 1
      u2(H_D1G33)%ysym = 1
      u2(H_D1G33)%xsym = -1
      u0(H_D1G33)%stiff = 0
      u2(H_D1G33)%stiff = 0
      u2(H_D1G33)%take_dx = 1
      u0(H_D1G33)%take_dx = 1
      u2(H_D1G33)%take_dy = 1
      u0(H_D1G33)%take_dy = 1
      u2(H_D1G33)%take_dz = 1
      u0(H_D1G33)%take_dz = 1
      u0(H_D2G00)%excise_val = 0
      u2(H_D2G00)%excise_val = 0
      u0(H_D2G00)%dissipation = 1
      u2(H_D2G00)%dissipation = 1
      u0(H_D2G00)%diss_factor = 1
      u2(H_D2G00)%diss_factor = 1
      u0(H_D2G00)%zsym = 1
      u0(H_D2G00)%ysym = -1
      u0(H_D2G00)%xsym = 1
      u2(H_D2G00)%zsym = 1
      u2(H_D2G00)%ysym = -1
      u2(H_D2G00)%xsym = 1
      u0(H_D2G00)%stiff = 0
      u2(H_D2G00)%stiff = 0
      u2(H_D2G00)%take_dx = 1
      u0(H_D2G00)%take_dx = 1
      u2(H_D2G00)%take_dy = 1
      u0(H_D2G00)%take_dy = 1
      u2(H_D2G00)%take_dz = 1
      u0(H_D2G00)%take_dz = 1
      u0(H_D2G01)%excise_val = 0
      u2(H_D2G01)%excise_val = 0
      u0(H_D2G01)%dissipation = 1
      u2(H_D2G01)%dissipation = 1
      u0(H_D2G01)%diss_factor = 1
      u2(H_D2G01)%diss_factor = 1
      u0(H_D2G01)%zsym = 1
      u0(H_D2G01)%ysym = -1
      u0(H_D2G01)%xsym = -1
      u2(H_D2G01)%zsym = 1
      u2(H_D2G01)%ysym = -1
      u2(H_D2G01)%xsym = -1
      u0(H_D2G01)%stiff = 0
      u2(H_D2G01)%stiff = 0
      u2(H_D2G01)%take_dx = 1
      u0(H_D2G01)%take_dx = 1
      u2(H_D2G01)%take_dy = 1
      u0(H_D2G01)%take_dy = 1
      u2(H_D2G01)%take_dz = 1
      u0(H_D2G01)%take_dz = 1
      u0(H_D2G02)%excise_val = 0
      u2(H_D2G02)%excise_val = 0
      u0(H_D2G02)%dissipation = 1
      u2(H_D2G02)%dissipation = 1
      u0(H_D2G02)%diss_factor = 1
      u2(H_D2G02)%diss_factor = 1
      u0(H_D2G02)%zsym = 1
      u0(H_D2G02)%ysym = 1
      u0(H_D2G02)%xsym = 1
      u2(H_D2G02)%zsym = 1
      u2(H_D2G02)%ysym = 1
      u2(H_D2G02)%xsym = 1
      u0(H_D2G02)%stiff = 0
      u2(H_D2G02)%stiff = 0
      u2(H_D2G02)%take_dx = 1
      u0(H_D2G02)%take_dx = 1
      u2(H_D2G02)%take_dy = 1
      u0(H_D2G02)%take_dy = 1
      u2(H_D2G02)%take_dz = 1
      u0(H_D2G02)%take_dz = 1
      u0(H_D2G03)%excise_val = 0
      u2(H_D2G03)%excise_val = 0
      u0(H_D2G03)%dissipation = 1
      u2(H_D2G03)%dissipation = 1
      u0(H_D2G03)%diss_factor = 1
      u2(H_D2G03)%diss_factor = 1
      u0(H_D2G03)%zsym = -1
      u0(H_D2G03)%ysym = -1
      u0(H_D2G03)%xsym = 1
      u2(H_D2G03)%zsym = -1
      u2(H_D2G03)%ysym = -1
      u2(H_D2G03)%xsym = 1
      u0(H_D2G03)%stiff = 0
      u2(H_D2G03)%stiff = 0
      u2(H_D2G03)%take_dx = 1
      u0(H_D2G03)%take_dx = 1
      u2(H_D2G03)%take_dy = 1
      u0(H_D2G03)%take_dy = 1
      u2(H_D2G03)%take_dz = 1
      u0(H_D2G03)%take_dz = 1
      u0(H_D2G11)%excise_val = 0
      u2(H_D2G11)%excise_val = 0
      u0(H_D2G11)%dissipation = 1
      u2(H_D2G11)%dissipation = 1
      u0(H_D2G11)%diss_factor = 1
      u2(H_D2G11)%diss_factor = 1
      u0(H_D2G11)%zsym = 1
      u0(H_D2G11)%ysym = -1
      u0(H_D2G11)%xsym = 1
      u2(H_D2G11)%zsym = 1
      u2(H_D2G11)%ysym = -1
      u2(H_D2G11)%xsym = 1
      u0(H_D2G11)%stiff = 0
      u2(H_D2G11)%stiff = 0
      u2(H_D2G11)%take_dx = 1
      u0(H_D2G11)%take_dx = 1
      u2(H_D2G11)%take_dy = 1
      u0(H_D2G11)%take_dy = 1
      u2(H_D2G11)%take_dz = 1
      u0(H_D2G11)%take_dz = 1
      u0(H_D2G12)%excise_val = 0
      u2(H_D2G12)%excise_val = 0
      u0(H_D2G12)%dissipation = 1
      u2(H_D2G12)%dissipation = 1
      u0(H_D2G12)%diss_factor = 1
      u2(H_D2G12)%diss_factor = 1
      u0(H_D2G12)%zsym = 1
      u0(H_D2G12)%ysym = 1
      u0(H_D2G12)%xsym = -1
      u2(H_D2G12)%zsym = 1
      u2(H_D2G12)%ysym = 1
      u2(H_D2G12)%xsym = -1
      u0(H_D2G12)%stiff = 0
      u2(H_D2G12)%stiff = 0
      u2(H_D2G12)%take_dx = 1
      u0(H_D2G12)%take_dx = 1
      u2(H_D2G12)%take_dy = 1
      u0(H_D2G12)%take_dy = 1
      u2(H_D2G12)%take_dz = 1
      u0(H_D2G12)%take_dz = 1
      u0(H_D2G13)%excise_val = 0
      u2(H_D2G13)%excise_val = 0
      u0(H_D2G13)%dissipation = 1
      u2(H_D2G13)%dissipation = 1
      u0(H_D2G13)%diss_factor = 1
      u2(H_D2G13)%diss_factor = 1
      u0(H_D2G13)%zsym = -1
      u0(H_D2G13)%ysym = -1
      u0(H_D2G13)%xsym = -1
      u2(H_D2G13)%zsym = -1
      u2(H_D2G13)%ysym = -1
      u2(H_D2G13)%xsym = -1
      u0(H_D2G13)%stiff = 0
      u2(H_D2G13)%stiff = 0
      u2(H_D2G13)%take_dx = 1
      u0(H_D2G13)%take_dx = 1
      u2(H_D2G13)%take_dy = 1
      u0(H_D2G13)%take_dy = 1
      u2(H_D2G13)%take_dz = 1
      u0(H_D2G13)%take_dz = 1
      u0(H_D2G22)%excise_val = 0
      u2(H_D2G22)%excise_val = 0
      u0(H_D2G22)%dissipation = 1
      u2(H_D2G22)%dissipation = 1
      u0(H_D2G22)%diss_factor = 1
      u2(H_D2G22)%diss_factor = 1
      u0(H_D2G22)%zsym = 1
      u0(H_D2G22)%ysym = -1
      u0(H_D2G22)%xsym = 1
      u2(H_D2G22)%zsym = 1
      u2(H_D2G22)%ysym = -1
      u2(H_D2G22)%xsym = 1
      u0(H_D2G22)%stiff = 0
      u2(H_D2G22)%stiff = 0
      u2(H_D2G22)%take_dx = 1
      u0(H_D2G22)%take_dx = 1
      u2(H_D2G22)%take_dy = 1
      u0(H_D2G22)%take_dy = 1
      u2(H_D2G22)%take_dz = 1
      u0(H_D2G22)%take_dz = 1
      u0(H_D2G23)%excise_val = 0
      u2(H_D2G23)%excise_val = 0
      u0(H_D2G23)%dissipation = 1
      u2(H_D2G23)%dissipation = 1
      u0(H_D2G23)%diss_factor = 1
      u2(H_D2G23)%diss_factor = 1
      u0(H_D2G23)%zsym = -1
      u0(H_D2G23)%ysym = 1
      u0(H_D2G23)%xsym = 1
      u2(H_D2G23)%zsym = -1
      u2(H_D2G23)%ysym = 1
      u2(H_D2G23)%xsym = 1
      u0(H_D2G23)%stiff = 0
      u2(H_D2G23)%stiff = 0
      u2(H_D2G23)%take_dx = 1
      u0(H_D2G23)%take_dx = 1
      u2(H_D2G23)%take_dy = 1
      u0(H_D2G23)%take_dy = 1
      u2(H_D2G23)%take_dz = 1
      u0(H_D2G23)%take_dz = 1
      u0(H_D2G33)%excise_val = 0
      u2(H_D2G33)%excise_val = 0
      u0(H_D2G33)%dissipation = 1
      u2(H_D2G33)%dissipation = 1
      u0(H_D2G33)%diss_factor = 1
      u2(H_D2G33)%diss_factor = 1
      u0(H_D2G33)%zsym = 1
      u0(H_D2G33)%ysym = -1
      u0(H_D2G33)%xsym = 1
      u2(H_D2G33)%zsym = 1
      u2(H_D2G33)%ysym = -1
      u2(H_D2G33)%xsym = 1
      u0(H_D2G33)%stiff = 0
      u2(H_D2G33)%stiff = 0
      u2(H_D2G33)%take_dx = 1
      u0(H_D2G33)%take_dx = 1
      u2(H_D2G33)%take_dy = 1
      u0(H_D2G33)%take_dy = 1
      u2(H_D2G33)%take_dz = 1
      u0(H_D2G33)%take_dz = 1
      u0(H_D3G00)%excise_val = 0
      u2(H_D3G00)%excise_val = 0
      u0(H_D3G00)%dissipation = 1
      u2(H_D3G00)%dissipation = 1
      u0(H_D3G00)%diss_factor = 1
      u2(H_D3G00)%diss_factor = 1
      u0(H_D3G00)%zsym = -1
      u0(H_D3G00)%ysym = 1
      u0(H_D3G00)%xsym = 1
      u2(H_D3G00)%zsym = -1
      u2(H_D3G00)%ysym = 1
      u2(H_D3G00)%xsym = 1
      u0(H_D3G00)%stiff = 0
      u2(H_D3G00)%stiff = 0
      u2(H_D3G00)%take_dx = 1
      u0(H_D3G00)%take_dx = 1
      u2(H_D3G00)%take_dy = 1
      u0(H_D3G00)%take_dy = 1
      u2(H_D3G00)%take_dz = 1
      u0(H_D3G00)%take_dz = 1
      u0(H_D3G01)%excise_val = 0
      u2(H_D3G01)%excise_val = 0
      u0(H_D3G01)%dissipation = 1
      u2(H_D3G01)%dissipation = 1
      u0(H_D3G01)%diss_factor = 1
      u2(H_D3G01)%diss_factor = 1
      u0(H_D3G01)%zsym = -1
      u0(H_D3G01)%ysym = 1
      u0(H_D3G01)%xsym = -1
      u2(H_D3G01)%zsym = -1
      u2(H_D3G01)%ysym = 1
      u2(H_D3G01)%xsym = -1
      u0(H_D3G01)%stiff = 0
      u2(H_D3G01)%stiff = 0
      u2(H_D3G01)%take_dx = 1
      u0(H_D3G01)%take_dx = 1
      u2(H_D3G01)%take_dy = 1
      u0(H_D3G01)%take_dy = 1
      u2(H_D3G01)%take_dz = 1
      u0(H_D3G01)%take_dz = 1
      u0(H_D3G02)%excise_val = 0
      u2(H_D3G02)%excise_val = 0
      u0(H_D3G02)%dissipation = 1
      u2(H_D3G02)%dissipation = 1
      u0(H_D3G02)%diss_factor = 1
      u2(H_D3G02)%diss_factor = 1
      u0(H_D3G02)%zsym = -1
      u0(H_D3G02)%ysym = -1
      u0(H_D3G02)%xsym = 1
      u2(H_D3G02)%zsym = -1
      u2(H_D3G02)%ysym = -1
      u2(H_D3G02)%xsym = 1
      u0(H_D3G02)%stiff = 0
      u2(H_D3G02)%stiff = 0
      u2(H_D3G02)%take_dx = 1
      u0(H_D3G02)%take_dx = 1
      u2(H_D3G02)%take_dy = 1
      u0(H_D3G02)%take_dy = 1
      u2(H_D3G02)%take_dz = 1
      u0(H_D3G02)%take_dz = 1
      u0(H_D3G03)%excise_val = 0
      u2(H_D3G03)%excise_val = 0
      u0(H_D3G03)%dissipation = 1
      u2(H_D3G03)%dissipation = 1
      u0(H_D3G03)%diss_factor = 1
      u2(H_D3G03)%diss_factor = 1
      u0(H_D3G03)%zsym = 1
      u0(H_D3G03)%ysym = 1
      u0(H_D3G03)%xsym = 1
      u2(H_D3G03)%zsym = 1
      u2(H_D3G03)%ysym = 1
      u2(H_D3G03)%xsym = 1
      u0(H_D3G03)%stiff = 0
      u2(H_D3G03)%stiff = 0
      u2(H_D3G03)%take_dx = 1
      u0(H_D3G03)%take_dx = 1
      u2(H_D3G03)%take_dy = 1
      u0(H_D3G03)%take_dy = 1
      u2(H_D3G03)%take_dz = 1
      u0(H_D3G03)%take_dz = 1
      u0(H_D3G11)%excise_val = 0
      u2(H_D3G11)%excise_val = 0
      u0(H_D3G11)%dissipation = 1
      u2(H_D3G11)%dissipation = 1
      u0(H_D3G11)%diss_factor = 1
      u2(H_D3G11)%diss_factor = 1
      u0(H_D3G11)%zsym = -1
      u0(H_D3G11)%ysym = 1
      u0(H_D3G11)%xsym = 1
      u2(H_D3G11)%zsym = -1
      u2(H_D3G11)%ysym = 1
      u2(H_D3G11)%xsym = 1
      u0(H_D3G11)%stiff = 0
      u2(H_D3G11)%stiff = 0
      u2(H_D3G11)%take_dx = 1
      u0(H_D3G11)%take_dx = 1
      u2(H_D3G11)%take_dy = 1
      u0(H_D3G11)%take_dy = 1
      u2(H_D3G11)%take_dz = 1
      u0(H_D3G11)%take_dz = 1
      u0(H_D3G12)%excise_val = 0
      u2(H_D3G12)%excise_val = 0
      u0(H_D3G12)%dissipation = 1
      u2(H_D3G12)%dissipation = 1
      u0(H_D3G12)%diss_factor = 1
      u2(H_D3G12)%diss_factor = 1
      u0(H_D3G12)%zsym = -1
      u0(H_D3G12)%ysym = -1
      u0(H_D3G12)%xsym = -1
      u2(H_D3G12)%zsym = -1
      u2(H_D3G12)%ysym = -1
      u2(H_D3G12)%xsym = -1
      u0(H_D3G12)%stiff = 0
      u2(H_D3G12)%stiff = 0
      u2(H_D3G12)%take_dx = 1
      u0(H_D3G12)%take_dx = 1
      u2(H_D3G12)%take_dy = 1
      u0(H_D3G12)%take_dy = 1
      u2(H_D3G12)%take_dz = 1
      u0(H_D3G12)%take_dz = 1
      u0(H_D3G13)%excise_val = 0
      u2(H_D3G13)%excise_val = 0
      u0(H_D3G13)%dissipation = 1
      u2(H_D3G13)%dissipation = 1
      u0(H_D3G13)%diss_factor = 1
      u2(H_D3G13)%diss_factor = 1
      u0(H_D3G13)%zsym = 1
      u0(H_D3G13)%ysym = 1
      u0(H_D3G13)%xsym = -1
      u2(H_D3G13)%zsym = 1
      u2(H_D3G13)%ysym = 1
      u2(H_D3G13)%xsym = -1
      u0(H_D3G13)%stiff = 0
      u2(H_D3G13)%stiff = 0
      u2(H_D3G13)%take_dx = 1
      u0(H_D3G13)%take_dx = 1
      u2(H_D3G13)%take_dy = 1
      u0(H_D3G13)%take_dy = 1
      u2(H_D3G13)%take_dz = 1
      u0(H_D3G13)%take_dz = 1
      u0(H_D3G22)%excise_val = 0
      u2(H_D3G22)%excise_val = 0
      u0(H_D3G22)%dissipation = 1
      u2(H_D3G22)%dissipation = 1
      u0(H_D3G22)%diss_factor = 1
      u2(H_D3G22)%diss_factor = 1
      u0(H_D3G22)%zsym = -1
      u0(H_D3G22)%ysym = 1
      u0(H_D3G22)%xsym = 1
      u2(H_D3G22)%zsym = -1
      u2(H_D3G22)%ysym = 1
      u2(H_D3G22)%xsym = 1
      u0(H_D3G22)%stiff = 0
      u2(H_D3G22)%stiff = 0
      u2(H_D3G22)%take_dx = 1
      u0(H_D3G22)%take_dx = 1
      u2(H_D3G22)%take_dy = 1
      u0(H_D3G22)%take_dy = 1
      u2(H_D3G22)%take_dz = 1
      u0(H_D3G22)%take_dz = 1
      u0(H_D3G23)%excise_val = 0
      u2(H_D3G23)%excise_val = 0
      u0(H_D3G23)%dissipation = 1
      u2(H_D3G23)%dissipation = 1
      u0(H_D3G23)%diss_factor = 1
      u2(H_D3G23)%diss_factor = 1
      u0(H_D3G23)%zsym = 1
      u0(H_D3G23)%ysym = -1
      u0(H_D3G23)%xsym = 1
      u2(H_D3G23)%zsym = 1
      u2(H_D3G23)%ysym = -1
      u2(H_D3G23)%xsym = 1
      u0(H_D3G23)%stiff = 0
      u2(H_D3G23)%stiff = 0
      u2(H_D3G23)%take_dx = 1
      u0(H_D3G23)%take_dx = 1
      u2(H_D3G23)%take_dy = 1
      u0(H_D3G23)%take_dy = 1
      u2(H_D3G23)%take_dz = 1
      u0(H_D3G23)%take_dz = 1
      u0(H_D3G33)%excise_val = 0
      u2(H_D3G33)%excise_val = 0
      u0(H_D3G33)%dissipation = 1
      u2(H_D3G33)%dissipation = 1
      u0(H_D3G33)%diss_factor = 1
      u2(H_D3G33)%diss_factor = 1
      u0(H_D3G33)%zsym = -1
      u0(H_D3G33)%ysym = 1
      u0(H_D3G33)%xsym = 1
      u2(H_D3G33)%zsym = -1
      u2(H_D3G33)%ysym = 1
      u2(H_D3G33)%xsym = 1
      u0(H_D3G33)%stiff = 0
      u2(H_D3G33)%stiff = 0
      u2(H_D3G33)%take_dx = 1
      u0(H_D3G33)%take_dx = 1
      u2(H_D3G33)%take_dy = 1
      u0(H_D3G33)%take_dy = 1
      u2(H_D3G33)%take_dz = 1
      u0(H_D3G33)%take_dz = 1
      u0(H_H0)%excise_val = 0
      u2(H_H0)%excise_val = 0
      u0(H_H0)%dissipation = 1
      u2(H_H0)%dissipation = 1
      u0(H_H0)%diss_factor = 1
      u2(H_H0)%diss_factor = 1
      u0(H_H0)%zsym = 1
      u0(H_H0)%ysym = 1
      u0(H_H0)%xsym = 1
      u2(H_H0)%zsym = 1
      u2(H_H0)%ysym = 1
      u2(H_H0)%xsym = 1
      u0(H_H0)%stiff = 0
      u2(H_H0)%stiff = 0
      u2(H_H0)%take_dx = 1
      u0(H_H0)%take_dx = 1
      u2(H_H0)%take_dy = 1
      u0(H_H0)%take_dy = 1
      u2(H_H0)%take_dz = 1
      u0(H_H0)%take_dz = 1
      u0(H_H1)%excise_val = 0
      u2(H_H1)%excise_val = 0
      u0(H_H1)%dissipation = 1
      u2(H_H1)%dissipation = 1
      u0(H_H1)%diss_factor = 1
      u2(H_H1)%diss_factor = 1
      u0(H_H1)%zsym = 1
      u0(H_H1)%ysym = 1
      u0(H_H1)%xsym = 1
      u2(H_H1)%zsym = 1
      u2(H_H1)%ysym = 1
      u2(H_H1)%xsym = 1
      u0(H_H1)%stiff = 0
      u2(H_H1)%stiff = 0
      u2(H_H1)%take_dx = 1
      u0(H_H1)%take_dx = 1
      u2(H_H1)%take_dy = 1
      u0(H_H1)%take_dy = 1
      u2(H_H1)%take_dz = 1
      u0(H_H1)%take_dz = 1
      u0(H_H2)%excise_val = 0
      u2(H_H2)%excise_val = 0
      u0(H_H2)%dissipation = 1
      u2(H_H2)%dissipation = 1
      u0(H_H2)%diss_factor = 1
      u2(H_H2)%diss_factor = 1
      u0(H_H2)%zsym = 1
      u0(H_H2)%ysym = 1
      u0(H_H2)%xsym = 1
      u2(H_H2)%zsym = 1
      u2(H_H2)%ysym = 1
      u2(H_H2)%xsym = 1
      u0(H_H2)%stiff = 0
      u2(H_H2)%stiff = 0
      u2(H_H2)%take_dx = 1
      u0(H_H2)%take_dx = 1
      u2(H_H2)%take_dy = 1
      u0(H_H2)%take_dy = 1
      u2(H_H2)%take_dz = 1
      u0(H_H2)%take_dz = 1
      u0(H_H3)%excise_val = 0
      u2(H_H3)%excise_val = 0
      u0(H_H3)%dissipation = 1
      u2(H_H3)%dissipation = 1
      u0(H_H3)%diss_factor = 1
      u2(H_H3)%diss_factor = 1
      u0(H_H3)%zsym = 1
      u0(H_H3)%ysym = 1
      u0(H_H3)%xsym = 1
      u2(H_H3)%zsym = 1
      u2(H_H3)%ysym = 1
      u2(H_H3)%xsym = 1
      u0(H_H3)%stiff = 0
      u2(H_H3)%stiff = 0
      u2(H_H3)%take_dx = 1
      u0(H_H3)%take_dx = 1
      u2(H_H3)%take_dy = 1
      u0(H_H3)%take_dy = 1
      u2(H_H3)%take_dz = 1
      u0(H_H3)%take_dz = 1
      u0(H_G0)%excise_val = 0
      u2(H_G0)%excise_val = 0
      u0(H_G0)%dissipation = 1
      u2(H_G0)%dissipation = 1
      u0(H_G0)%diss_factor = 1
      u2(H_G0)%diss_factor = 1
      u0(H_G0)%zsym = 1
      u0(H_G0)%ysym = 1
      u0(H_G0)%xsym = 1
      u2(H_G0)%zsym = 1
      u2(H_G0)%ysym = 1
      u2(H_G0)%xsym = 1
      u0(H_G0)%stiff = 0
      u2(H_G0)%stiff = 0
      u2(H_G0)%take_dx = 1
      u0(H_G0)%take_dx = 1
      u2(H_G0)%take_dy = 1
      u0(H_G0)%take_dy = 1
      u2(H_G0)%take_dz = 1
      u0(H_G0)%take_dz = 1
      u0(H_D1H0)%excise_val = 0
      u2(H_D1H0)%excise_val = 0
      u0(H_D1H0)%dissipation = 1
      u2(H_D1H0)%dissipation = 1
      u0(H_D1H0)%diss_factor = 1
      u2(H_D1H0)%diss_factor = 1
      u0(H_D1H0)%zsym = 1
      u0(H_D1H0)%ysym = 1
      u0(H_D1H0)%xsym = -1
      u2(H_D1H0)%zsym = 1
      u2(H_D1H0)%ysym = 1
      u2(H_D1H0)%xsym = -1
      u0(H_D1H0)%stiff = 0
      u2(H_D1H0)%stiff = 0
      u2(H_D1H0)%take_dx = 1
      u0(H_D1H0)%take_dx = 1
      u2(H_D1H0)%take_dy = 1
      u0(H_D1H0)%take_dy = 1
      u2(H_D1H0)%take_dz = 1
      u0(H_D1H0)%take_dz = 1
      u0(H_D2H0)%excise_val = 0
      u2(H_D2H0)%excise_val = 0
      u0(H_D2H0)%dissipation = 1
      u2(H_D2H0)%dissipation = 1
      u0(H_D2H0)%diss_factor = 1
      u2(H_D2H0)%diss_factor = 1
      u0(H_D2H0)%zsym = 1
      u0(H_D2H0)%ysym = -1
      u0(H_D2H0)%xsym = 1
      u2(H_D2H0)%zsym = 1
      u2(H_D2H0)%ysym = -1
      u2(H_D2H0)%xsym = 1
      u0(H_D2H0)%stiff = 0
      u2(H_D2H0)%stiff = 0
      u2(H_D2H0)%take_dx = 1
      u0(H_D2H0)%take_dx = 1
      u2(H_D2H0)%take_dy = 1
      u0(H_D2H0)%take_dy = 1
      u2(H_D2H0)%take_dz = 1
      u0(H_D2H0)%take_dz = 1
      u0(H_D3H0)%excise_val = 0
      u2(H_D3H0)%excise_val = 0
      u0(H_D3H0)%dissipation = 1
      u2(H_D3H0)%dissipation = 1
      u0(H_D3H0)%diss_factor = 1
      u2(H_D3H0)%diss_factor = 1
      u0(H_D3H0)%zsym = -1
      u0(H_D3H0)%ysym = 1
      u0(H_D3H0)%xsym = 1
      u2(H_D3H0)%zsym = -1
      u2(H_D3H0)%ysym = 1
      u2(H_D3H0)%xsym = 1
      u0(H_D3H0)%stiff = 0
      u2(H_D3H0)%stiff = 0
      u2(H_D3H0)%take_dx = 1
      u0(H_D3H0)%take_dx = 1
      u2(H_D3H0)%take_dy = 1
      u0(H_D3H0)%take_dy = 1
      u2(H_D3H0)%take_dz = 1
      u0(H_D3H0)%take_dz = 1
      u0(H_PHIM)%excise_val = 0
      u2(H_PHIM)%excise_val = 0
      u0(H_PHIM)%dissipation = 1
      u2(H_PHIM)%dissipation = 1
      u0(H_PHIM)%diss_factor = 1
      u2(H_PHIM)%diss_factor = 1
      u0(H_PHIM)%zsym = 1
      u0(H_PHIM)%ysym = 1
      u0(H_PHIM)%xsym = 1
      u2(H_PHIM)%zsym = 1
      u2(H_PHIM)%ysym = 1
      u2(H_PHIM)%xsym = 1
      u0(H_PHIM)%stiff = 0
      u2(H_PHIM)%stiff = 0
      u2(H_PHIM)%take_dx = 1
      u0(H_PHIM)%take_dx = 1
      u2(H_PHIM)%take_dy = 1
      u0(H_PHIM)%take_dy = 1
      u2(H_PHIM)%take_dz = 1
      u0(H_PHIM)%take_dz = 1
      u0(H_PIM)%excise_val = 0
      u2(H_PIM)%excise_val = 0
      u0(H_PIM)%dissipation = 1
      u2(H_PIM)%dissipation = 1
      u0(H_PIM)%diss_factor = 1
      u2(H_PIM)%diss_factor = 1
      u0(H_PIM)%zsym = 1
      u0(H_PIM)%ysym = 1
      u0(H_PIM)%xsym = 1
      u2(H_PIM)%zsym = 1
      u2(H_PIM)%ysym = 1
      u2(H_PIM)%xsym = 1
      u0(H_PIM)%stiff = 0
      u2(H_PIM)%stiff = 0
      u2(H_PIM)%take_dx = 1
      u0(H_PIM)%take_dx = 1
      u2(H_PIM)%take_dy = 1
      u0(H_PIM)%take_dy = 1
      u2(H_PIM)%take_dz = 1
      u0(H_PIM)%take_dz = 1
      u0(H_D1PHIM)%excise_val = 0
      u2(H_D1PHIM)%excise_val = 0
      u0(H_D1PHIM)%dissipation = 1
      u2(H_D1PHIM)%dissipation = 1
      u0(H_D1PHIM)%diss_factor = 1
      u2(H_D1PHIM)%diss_factor = 1
      u0(H_D1PHIM)%zsym = 1
      u0(H_D1PHIM)%ysym = 1
      u0(H_D1PHIM)%xsym = -1
      u2(H_D1PHIM)%zsym = 1
      u2(H_D1PHIM)%ysym = 1
      u2(H_D1PHIM)%xsym = -1
      u0(H_D1PHIM)%stiff = 0
      u2(H_D1PHIM)%stiff = 0
      u2(H_D1PHIM)%take_dx = 1
      u0(H_D1PHIM)%take_dx = 1
      u2(H_D1PHIM)%take_dy = 1
      u0(H_D1PHIM)%take_dy = 1
      u2(H_D1PHIM)%take_dz = 1
      u0(H_D1PHIM)%take_dz = 1
      u0(H_D2PHIM)%excise_val = 0
      u2(H_D2PHIM)%excise_val = 0
      u0(H_D2PHIM)%dissipation = 1
      u2(H_D2PHIM)%dissipation = 1
      u0(H_D2PHIM)%diss_factor = 1
      u2(H_D2PHIM)%diss_factor = 1
      u0(H_D2PHIM)%zsym = 1
      u0(H_D2PHIM)%ysym = -1
      u0(H_D2PHIM)%xsym = 1
      u2(H_D2PHIM)%zsym = 1
      u2(H_D2PHIM)%ysym = -1
      u2(H_D2PHIM)%xsym = 1
      u0(H_D2PHIM)%stiff = 0
      u2(H_D2PHIM)%stiff = 0
      u2(H_D2PHIM)%take_dx = 1
      u0(H_D2PHIM)%take_dx = 1
      u2(H_D2PHIM)%take_dy = 1
      u0(H_D2PHIM)%take_dy = 1
      u2(H_D2PHIM)%take_dz = 1
      u0(H_D2PHIM)%take_dz = 1
      u0(H_D3PHIM)%excise_val = 0
      u2(H_D3PHIM)%excise_val = 0
      u0(H_D3PHIM)%dissipation = 1
      u2(H_D3PHIM)%dissipation = 1
      u0(H_D3PHIM)%diss_factor = 1
      u2(H_D3PHIM)%diss_factor = 1
      u0(H_D3PHIM)%zsym = -1
      u0(H_D3PHIM)%ysym = 1
      u0(H_D3PHIM)%xsym = 1
      u2(H_D3PHIM)%zsym = -1
      u2(H_D3PHIM)%ysym = 1
      u2(H_D3PHIM)%xsym = 1
      u0(H_D3PHIM)%stiff = 0
      u2(H_D3PHIM)%stiff = 0
      u2(H_D3PHIM)%take_dx = 1
      u0(H_D3PHIM)%take_dx = 1
      u2(H_D3PHIM)%take_dy = 1
      u0(H_D3PHIM)%take_dy = 1
      u2(H_D3PHIM)%take_dz = 1
      u0(H_D3PHIM)%take_dz = 1
      u0(H_PHIN)%excise_val = 0
      u2(H_PHIN)%excise_val = 0
      u0(H_PHIN)%dissipation = 1
      u2(H_PHIN)%dissipation = 1
      u0(H_PHIN)%diss_factor = 1
      u2(H_PHIN)%diss_factor = 1
      u0(H_PHIN)%zsym = 1
      u0(H_PHIN)%ysym = 1
      u0(H_PHIN)%xsym = 1
      u2(H_PHIN)%zsym = 1
      u2(H_PHIN)%ysym = 1
      u2(H_PHIN)%xsym = 1
      u0(H_PHIN)%stiff = 0
      u2(H_PHIN)%stiff = 0
      u2(H_PHIN)%take_dx = 1
      u0(H_PHIN)%take_dx = 1
      u2(H_PHIN)%take_dy = 1
      u0(H_PHIN)%take_dy = 1
      u2(H_PHIN)%take_dz = 1
      u0(H_PHIN)%take_dz = 1
      u0(H_PIN)%excise_val = 0
      u2(H_PIN)%excise_val = 0
      u0(H_PIN)%dissipation = 1
      u2(H_PIN)%dissipation = 1
      u0(H_PIN)%diss_factor = 1
      u2(H_PIN)%diss_factor = 1
      u0(H_PIN)%zsym = 1
      u0(H_PIN)%ysym = 1
      u0(H_PIN)%xsym = 1
      u2(H_PIN)%zsym = 1
      u2(H_PIN)%ysym = 1
      u2(H_PIN)%xsym = 1
      u0(H_PIN)%stiff = 0
      u2(H_PIN)%stiff = 0
      u2(H_PIN)%take_dx = 1
      u0(H_PIN)%take_dx = 1
      u2(H_PIN)%take_dy = 1
      u0(H_PIN)%take_dy = 1
      u2(H_PIN)%take_dz = 1
      u0(H_PIN)%take_dz = 1
      u0(H_D1PHIN)%excise_val = 0
      u2(H_D1PHIN)%excise_val = 0
      u0(H_D1PHIN)%dissipation = 1
      u2(H_D1PHIN)%dissipation = 1
      u0(H_D1PHIN)%diss_factor = 1
      u2(H_D1PHIN)%diss_factor = 1
      u0(H_D1PHIN)%zsym = 1
      u0(H_D1PHIN)%ysym = 1
      u0(H_D1PHIN)%xsym = -1
      u2(H_D1PHIN)%zsym = 1
      u2(H_D1PHIN)%ysym = 1
      u2(H_D1PHIN)%xsym = -1
      u0(H_D1PHIN)%stiff = 0
      u2(H_D1PHIN)%stiff = 0
      u2(H_D1PHIN)%take_dx = 1
      u0(H_D1PHIN)%take_dx = 1
      u2(H_D1PHIN)%take_dy = 1
      u0(H_D1PHIN)%take_dy = 1
      u2(H_D1PHIN)%take_dz = 1
      u0(H_D1PHIN)%take_dz = 1
      u0(H_D2PHIN)%excise_val = 0
      u2(H_D2PHIN)%excise_val = 0
      u0(H_D2PHIN)%dissipation = 1
      u2(H_D2PHIN)%dissipation = 1
      u0(H_D2PHIN)%diss_factor = 1
      u2(H_D2PHIN)%diss_factor = 1
      u0(H_D2PHIN)%zsym = 1
      u0(H_D2PHIN)%ysym = -1
      u0(H_D2PHIN)%xsym = 1
      u2(H_D2PHIN)%zsym = 1
      u2(H_D2PHIN)%ysym = -1
      u2(H_D2PHIN)%xsym = 1
      u0(H_D2PHIN)%stiff = 0
      u2(H_D2PHIN)%stiff = 0
      u2(H_D2PHIN)%take_dx = 1
      u0(H_D2PHIN)%take_dx = 1
      u2(H_D2PHIN)%take_dy = 1
      u0(H_D2PHIN)%take_dy = 1
      u2(H_D2PHIN)%take_dz = 1
      u0(H_D2PHIN)%take_dz = 1
      u0(H_D3PHIN)%excise_val = 0
      u2(H_D3PHIN)%excise_val = 0
      u0(H_D3PHIN)%dissipation = 1
      u2(H_D3PHIN)%dissipation = 1
      u0(H_D3PHIN)%diss_factor = 1
      u2(H_D3PHIN)%diss_factor = 1
      u0(H_D3PHIN)%zsym = -1
      u0(H_D3PHIN)%ysym = 1
      u0(H_D3PHIN)%xsym = 1
      u2(H_D3PHIN)%zsym = -1
      u2(H_D3PHIN)%ysym = 1
      u2(H_D3PHIN)%xsym = 1
      u0(H_D3PHIN)%stiff = 0
      u2(H_D3PHIN)%stiff = 0
      u2(H_D3PHIN)%take_dx = 1
      u0(H_D3PHIN)%take_dx = 1
      u2(H_D3PHIN)%take_dy = 1
      u0(H_D3PHIN)%take_dy = 1
      u2(H_D3PHIN)%take_dz = 1
      u0(H_D3PHIN)%take_dz = 1
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
      u2(H_PHIR)%take_dx = 1
      u0(H_PHIR)%take_dx = 1
      u2(H_PHIR)%take_dy = 1
      u0(H_PHIR)%take_dy = 1
      u2(H_PHIR)%take_dz = 1
      u0(H_PHIR)%take_dz = 1
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
      u2(H_PIR)%take_dx = 1
      u0(H_PIR)%take_dx = 1
      u2(H_PIR)%take_dy = 1
      u0(H_PIR)%take_dy = 1
      u2(H_PIR)%take_dz = 1
      u0(H_PIR)%take_dz = 1
      u0(H_D1PHIR)%excise_val = 0
      u2(H_D1PHIR)%excise_val = 0
      u0(H_D1PHIR)%dissipation = 1
      u2(H_D1PHIR)%dissipation = 1
      u0(H_D1PHIR)%diss_factor = 1
      u2(H_D1PHIR)%diss_factor = 1
      u0(H_D1PHIR)%zsym = 1
      u0(H_D1PHIR)%ysym = 1
      u0(H_D1PHIR)%xsym = -1
      u2(H_D1PHIR)%zsym = 1
      u2(H_D1PHIR)%ysym = 1
      u2(H_D1PHIR)%xsym = -1
      u0(H_D1PHIR)%stiff = 0
      u2(H_D1PHIR)%stiff = 0
      u2(H_D1PHIR)%take_dx = 1
      u0(H_D1PHIR)%take_dx = 1
      u2(H_D1PHIR)%take_dy = 1
      u0(H_D1PHIR)%take_dy = 1
      u2(H_D1PHIR)%take_dz = 1
      u0(H_D1PHIR)%take_dz = 1
      u0(H_D2PHIR)%excise_val = 0
      u2(H_D2PHIR)%excise_val = 0
      u0(H_D2PHIR)%dissipation = 1
      u2(H_D2PHIR)%dissipation = 1
      u0(H_D2PHIR)%diss_factor = 1
      u2(H_D2PHIR)%diss_factor = 1
      u0(H_D2PHIR)%zsym = 1
      u0(H_D2PHIR)%ysym = -1
      u0(H_D2PHIR)%xsym = 1
      u2(H_D2PHIR)%zsym = 1
      u2(H_D2PHIR)%ysym = -1
      u2(H_D2PHIR)%xsym = 1
      u0(H_D2PHIR)%stiff = 0
      u2(H_D2PHIR)%stiff = 0
      u2(H_D2PHIR)%take_dx = 1
      u0(H_D2PHIR)%take_dx = 1
      u2(H_D2PHIR)%take_dy = 1
      u0(H_D2PHIR)%take_dy = 1
      u2(H_D2PHIR)%take_dz = 1
      u0(H_D2PHIR)%take_dz = 1
      u0(H_D3PHIR)%excise_val = 0
      u2(H_D3PHIR)%excise_val = 0
      u0(H_D3PHIR)%dissipation = 1
      u2(H_D3PHIR)%dissipation = 1
      u0(H_D3PHIR)%diss_factor = 1
      u2(H_D3PHIR)%diss_factor = 1
      u0(H_D3PHIR)%zsym = -1
      u0(H_D3PHIR)%ysym = 1
      u0(H_D3PHIR)%xsym = 1
      u2(H_D3PHIR)%zsym = -1
      u2(H_D3PHIR)%ysym = 1
      u2(H_D3PHIR)%xsym = 1
      u0(H_D3PHIR)%stiff = 0
      u2(H_D3PHIR)%stiff = 0
      u2(H_D3PHIR)%take_dx = 1
      u0(H_D3PHIR)%take_dx = 1
      u2(H_D3PHIR)%take_dy = 1
      u0(H_D3PHIR)%take_dy = 1
      u2(H_D3PHIR)%take_dz = 1
      u0(H_D3PHIR)%take_dz = 1
      u0(H_PHIC)%excise_val = 0
      u2(H_PHIC)%excise_val = 0
      u0(H_PHIC)%dissipation = 1
      u2(H_PHIC)%dissipation = 1
      u0(H_PHIC)%diss_factor = 1
      u2(H_PHIC)%diss_factor = 1
      u0(H_PHIC)%zsym = 1
      u0(H_PHIC)%ysym = 1
      u0(H_PHIC)%xsym = 1
      u2(H_PHIC)%zsym = 1
      u2(H_PHIC)%ysym = 1
      u2(H_PHIC)%xsym = 1
      u0(H_PHIC)%stiff = 0
      u2(H_PHIC)%stiff = 0
      u2(H_PHIC)%take_dx = 1
      u0(H_PHIC)%take_dx = 1
      u2(H_PHIC)%take_dy = 1
      u0(H_PHIC)%take_dy = 1
      u2(H_PHIC)%take_dz = 1
      u0(H_PHIC)%take_dz = 1
      u0(H_PIC)%excise_val = 0
      u2(H_PIC)%excise_val = 0
      u0(H_PIC)%dissipation = 1
      u2(H_PIC)%dissipation = 1
      u0(H_PIC)%diss_factor = 1
      u2(H_PIC)%diss_factor = 1
      u0(H_PIC)%zsym = 1
      u0(H_PIC)%ysym = 1
      u0(H_PIC)%xsym = 1
      u2(H_PIC)%zsym = 1
      u2(H_PIC)%ysym = 1
      u2(H_PIC)%xsym = 1
      u0(H_PIC)%stiff = 0
      u2(H_PIC)%stiff = 0
      u2(H_PIC)%take_dx = 1
      u0(H_PIC)%take_dx = 1
      u2(H_PIC)%take_dy = 1
      u0(H_PIC)%take_dy = 1
      u2(H_PIC)%take_dz = 1
      u0(H_PIC)%take_dz = 1
      u0(H_D1PHIC)%excise_val = 0
      u2(H_D1PHIC)%excise_val = 0
      u0(H_D1PHIC)%dissipation = 1
      u2(H_D1PHIC)%dissipation = 1
      u0(H_D1PHIC)%diss_factor = 1
      u2(H_D1PHIC)%diss_factor = 1
      u0(H_D1PHIC)%zsym = 1
      u0(H_D1PHIC)%ysym = 1
      u0(H_D1PHIC)%xsym = -1
      u2(H_D1PHIC)%zsym = 1
      u2(H_D1PHIC)%ysym = 1
      u2(H_D1PHIC)%xsym = -1
      u0(H_D1PHIC)%stiff = 0
      u2(H_D1PHIC)%stiff = 0
      u2(H_D1PHIC)%take_dx = 1
      u0(H_D1PHIC)%take_dx = 1
      u2(H_D1PHIC)%take_dy = 1
      u0(H_D1PHIC)%take_dy = 1
      u2(H_D1PHIC)%take_dz = 1
      u0(H_D1PHIC)%take_dz = 1
      u0(H_D2PHIC)%excise_val = 0
      u2(H_D2PHIC)%excise_val = 0
      u0(H_D2PHIC)%dissipation = 1
      u2(H_D2PHIC)%dissipation = 1
      u0(H_D2PHIC)%diss_factor = 1
      u2(H_D2PHIC)%diss_factor = 1
      u0(H_D2PHIC)%zsym = 1
      u0(H_D2PHIC)%ysym = -1
      u0(H_D2PHIC)%xsym = 1
      u2(H_D2PHIC)%zsym = 1
      u2(H_D2PHIC)%ysym = -1
      u2(H_D2PHIC)%xsym = 1
      u0(H_D2PHIC)%stiff = 0
      u2(H_D2PHIC)%stiff = 0
      u2(H_D2PHIC)%take_dx = 1
      u0(H_D2PHIC)%take_dx = 1
      u2(H_D2PHIC)%take_dy = 1
      u0(H_D2PHIC)%take_dy = 1
      u2(H_D2PHIC)%take_dz = 1
      u0(H_D2PHIC)%take_dz = 1
      u0(H_D3PHIC)%excise_val = 0
      u2(H_D3PHIC)%excise_val = 0
      u0(H_D3PHIC)%dissipation = 1
      u2(H_D3PHIC)%dissipation = 1
      u0(H_D3PHIC)%diss_factor = 1
      u2(H_D3PHIC)%diss_factor = 1
      u0(H_D3PHIC)%zsym = -1
      u0(H_D3PHIC)%ysym = 1
      u0(H_D3PHIC)%xsym = 1
      u2(H_D3PHIC)%zsym = -1
      u2(H_D3PHIC)%ysym = 1
      u2(H_D3PHIC)%xsym = 1
      u0(H_D3PHIC)%stiff = 0
      u2(H_D3PHIC)%stiff = 0
      u2(H_D3PHIC)%take_dx = 1
      u0(H_D3PHIC)%take_dx = 1
      u2(H_D3PHIC)%take_dy = 1
      u0(H_D3PHIC)%take_dy = 1
      u2(H_D3PHIC)%take_dz = 1
      u0(H_D3PHIC)%take_dz = 1
      v(H_ALPHA)%take_dx = 0
      v(H_ALPHA)%take_dy = 0
      v(H_ALPHA)%take_dz = 0
      v(H_SHIFT1)%take_dx = 0
      v(H_SHIFT1)%take_dy = 0
      v(H_SHIFT1)%take_dz = 0
      v(H_SHIFT2)%take_dx = 0
      v(H_SHIFT2)%take_dy = 0
      v(H_SHIFT2)%take_dz = 0
      v(H_SHIFT3)%take_dx = 0
      v(H_SHIFT3)%take_dy = 0
      v(H_SHIFT3)%take_dz = 0
      v(H_SQDETG)%take_dx = 0
      v(H_SQDETG)%take_dy = 0
      v(H_SQDETG)%take_dz = 0
      v(H_TMUNU00)%take_dx = 0
      v(H_TMUNU00)%take_dy = 0
      v(H_TMUNU00)%take_dz = 0
      v(H_TMUNU01)%take_dx = 0
      v(H_TMUNU01)%take_dy = 0
      v(H_TMUNU01)%take_dz = 0
      v(H_TMUNU02)%take_dx = 0
      v(H_TMUNU02)%take_dy = 0
      v(H_TMUNU02)%take_dz = 0
      v(H_TMUNU03)%take_dx = 0
      v(H_TMUNU03)%take_dy = 0
      v(H_TMUNU03)%take_dz = 0
      v(H_TMUNU11)%take_dx = 0
      v(H_TMUNU11)%take_dy = 0
      v(H_TMUNU11)%take_dz = 0
      v(H_TMUNU12)%take_dx = 0
      v(H_TMUNU12)%take_dy = 0
      v(H_TMUNU12)%take_dz = 0
      v(H_TMUNU13)%take_dx = 0
      v(H_TMUNU13)%take_dy = 0
      v(H_TMUNU13)%take_dz = 0
      v(H_TMUNU22)%take_dx = 0
      v(H_TMUNU22)%take_dy = 0
      v(H_TMUNU22)%take_dz = 0
      v(H_TMUNU23)%take_dx = 0
      v(H_TMUNU23)%take_dy = 0
      v(H_TMUNU23)%take_dz = 0
      v(H_TMUNU33)%take_dx = 0
      v(H_TMUNU33)%take_dy = 0
      v(H_TMUNU33)%take_dz = 0
      v(H_ORDER1_CONST)%take_dx = 0
      v(H_ORDER1_CONST)%take_dy = 0
      v(H_ORDER1_CONST)%take_dz = 0
      v(H_ORDER2_CONST)%take_dx = 0
      v(H_ORDER2_CONST)%take_dy = 0
      v(H_ORDER2_CONST)%take_dz = 0
      v(H_PHYSZ_CONST)%take_dx = 0
      v(H_PHYSZ_CONST)%take_dy = 0
      v(H_PHYSZ_CONST)%take_dz = 0
      v(H_RAD_EXP)%take_dx = 0
      v(H_RAD_EXP)%take_dy = 0
      v(H_RAD_EXP)%take_dz = 0
      v(H_RAD_SPEED)%take_dx = 0
      v(H_RAD_SPEED)%take_dy = 0
      v(H_RAD_SPEED)%take_dz = 0
      v(H_ENERGY_DENS)%take_dx = 0
      v(H_ENERGY_DENS)%take_dy = 0
      v(H_ENERGY_DENS)%take_dz = 0
      v(H_NOETHER_DENS)%take_dx = 0
      v(H_NOETHER_DENS)%take_dy = 0
      v(H_NOETHER_DENS)%take_dz = 0
      v(H_ENERGY_CONST)%take_dx = 0
      v(H_ENERGY_CONST)%take_dy = 0
      v(H_ENERGY_CONST)%take_dz = 0
      v(H_MOMX_CONST)%take_dx = 0
      v(H_MOMX_CONST)%take_dy = 0
      v(H_MOMX_CONST)%take_dz = 0
      v(H_MOMY_CONST)%take_dx = 0
      v(H_MOMY_CONST)%take_dy = 0
      v(H_MOMY_CONST)%take_dz = 0
      v(H_MOMZ_CONST)%take_dx = 0
      v(H_MOMZ_CONST)%take_dy = 0
      v(H_MOMZ_CONST)%take_dz = 0
      v(H_PSI4R)%take_dx = 0
      v(H_PSI4R)%take_dy = 0
      v(H_PSI4R)%take_dz = 0
      v(H_PSI4I)%take_dx = 0
      v(H_PSI4I)%take_dy = 0
      v(H_PSI4I)%take_dz = 0
      v(H_MASSADM)%take_dx = 0
      v(H_MASSADM)%take_dy = 0
      v(H_MASSADM)%take_dz = 0
      v(H_CURVATURE)%take_dx = 0
      v(H_CURVATURE)%take_dy = 0
      v(H_CURVATURE)%take_dz = 0
      v(H_LIN_MOMX)%take_dx = 0
      v(H_LIN_MOMX)%take_dy = 0
      v(H_LIN_MOMX)%take_dz = 0
      v(H_LIN_MOMY)%take_dx = 0
      v(H_LIN_MOMY)%take_dy = 0
      v(H_LIN_MOMY)%take_dz = 0
      v(H_ANG_MOMZ)%take_dx = 0
      v(H_ANG_MOMZ)%take_dy = 0
      v(H_ANG_MOMZ)%take_dz = 0
      v(H_KANG_MOMZ)%take_dx = 0
      v(H_KANG_MOMZ)%take_dy = 0
      v(H_KANG_MOMZ)%take_dz = 0
      v(H_R2SCALAR)%take_dx = 0
      v(H_R2SCALAR)%take_dy = 0
      v(H_R2SCALAR)%take_dz = 0
      v(H_GTHETHE)%take_dx = 0
      v(H_GTHETHE)%take_dy = 0
      v(H_GTHETHE)%take_dz = 0
      v(H_GTHEPHI)%take_dx = 0
      v(H_GTHEPHI)%take_dy = 0
      v(H_GTHEPHI)%take_dz = 0
      v(H_GPHIPHI)%take_dx = 0
      v(H_GPHIPHI)%take_dy = 0
      v(H_GPHIPHI)%take_dz = 0
      v(H_GUR)%take_dx = 0
      v(H_GUR)%take_dy = 0
      v(H_GUR)%take_dz = 0
      v(H_PSI_EL)%take_dx = 0
      v(H_PSI_EL)%take_dy = 0
      v(H_PSI_EL)%take_dz = 0
      v(H_BETAX_EL)%take_dx = 0
      v(H_BETAX_EL)%take_dy = 0
      v(H_BETAX_EL)%take_dz = 0
      v(H_BETAY_EL)%take_dx = 0
      v(H_BETAY_EL)%take_dy = 0
      v(H_BETAY_EL)%take_dz = 0
      v(H_BETAZ_EL)%take_dx = 0
      v(H_BETAZ_EL)%take_dy = 0
      v(H_BETAZ_EL)%take_dz = 0
      return
      end subroutine assign_ptrs_fields
      subroutine assign_ptrs_rks(g00_rk1,g01_rk1,g02_rk1,g03_rk1,g11_rk1
     &,g12_rk1,g13_rk1,g22_rk1,g23_rk1,g33_rk1,K00_rk1,K01_rk1,K02_rk1,K
     &03_rk1,K11_rk1,K12_rk1,K13_rk1,K22_rk1,K23_rk1,K33_rk1,d1g00_rk1,d
     &1g01_rk1,d1g02_rk1,d1g03_rk1,d1g11_rk1,d1g12_rk1,d1g13_rk1,d1g22_r
     &k1,d1g23_rk1,d1g33_rk1,d2g00_rk1,d2g01_rk1,d2g02_rk1,d2g03_rk1,d2g
     &11_rk1,d2g12_rk1,d2g13_rk1,d2g22_rk1,d2g23_rk1,d2g33_rk1,d3g00_rk1
     &,d3g01_rk1,d3g02_rk1,d3g03_rk1,d3g11_rk1,d3g12_rk1,d3g13_rk1,d3g22
     &_rk1,d3g23_rk1,d3g33_rk1,H0_rk1,H1_rk1,H2_rk1,H3_rk1,G0_rk1,d1H0_r
     &k1,d2H0_rk1,d3H0_rk1,phim_rk1,pim_rk1,d1phim_rk1,d2phim_rk1,d3phim
     &_rk1,phin_rk1,pin_rk1,d1phin_rk1,d2phin_rk1,d3phin_rk1,phir_rk1,pi
     &r_rk1,d1phir_rk1,d2phir_rk1,d3phir_rk1,phic_rk1,pic_rk1,d1phic_rk1
     &,d2phic_rk1,d3phic_rk1,urk1,nx,ny,nz)
      use GF
      implicit none
      integer nx, ny, nz
      type(gridfunction), dimension(NU) :: urk1
      real(kind=8) g00_rk1(nx,ny,nz)
      real(kind=8) g01_rk1(nx,ny,nz)
      real(kind=8) g02_rk1(nx,ny,nz)
      real(kind=8) g03_rk1(nx,ny,nz)
      real(kind=8) g11_rk1(nx,ny,nz)
      real(kind=8) g12_rk1(nx,ny,nz)
      real(kind=8) g13_rk1(nx,ny,nz)
      real(kind=8) g22_rk1(nx,ny,nz)
      real(kind=8) g23_rk1(nx,ny,nz)
      real(kind=8) g33_rk1(nx,ny,nz)
      real(kind=8) K00_rk1(nx,ny,nz)
      real(kind=8) K01_rk1(nx,ny,nz)
      real(kind=8) K02_rk1(nx,ny,nz)
      real(kind=8) K03_rk1(nx,ny,nz)
      real(kind=8) K11_rk1(nx,ny,nz)
      real(kind=8) K12_rk1(nx,ny,nz)
      real(kind=8) K13_rk1(nx,ny,nz)
      real(kind=8) K22_rk1(nx,ny,nz)
      real(kind=8) K23_rk1(nx,ny,nz)
      real(kind=8) K33_rk1(nx,ny,nz)
      real(kind=8) d1g00_rk1(nx,ny,nz)
      real(kind=8) d1g01_rk1(nx,ny,nz)
      real(kind=8) d1g02_rk1(nx,ny,nz)
      real(kind=8) d1g03_rk1(nx,ny,nz)
      real(kind=8) d1g11_rk1(nx,ny,nz)
      real(kind=8) d1g12_rk1(nx,ny,nz)
      real(kind=8) d1g13_rk1(nx,ny,nz)
      real(kind=8) d1g22_rk1(nx,ny,nz)
      real(kind=8) d1g23_rk1(nx,ny,nz)
      real(kind=8) d1g33_rk1(nx,ny,nz)
      real(kind=8) d2g00_rk1(nx,ny,nz)
      real(kind=8) d2g01_rk1(nx,ny,nz)
      real(kind=8) d2g02_rk1(nx,ny,nz)
      real(kind=8) d2g03_rk1(nx,ny,nz)
      real(kind=8) d2g11_rk1(nx,ny,nz)
      real(kind=8) d2g12_rk1(nx,ny,nz)
      real(kind=8) d2g13_rk1(nx,ny,nz)
      real(kind=8) d2g22_rk1(nx,ny,nz)
      real(kind=8) d2g23_rk1(nx,ny,nz)
      real(kind=8) d2g33_rk1(nx,ny,nz)
      real(kind=8) d3g00_rk1(nx,ny,nz)
      real(kind=8) d3g01_rk1(nx,ny,nz)
      real(kind=8) d3g02_rk1(nx,ny,nz)
      real(kind=8) d3g03_rk1(nx,ny,nz)
      real(kind=8) d3g11_rk1(nx,ny,nz)
      real(kind=8) d3g12_rk1(nx,ny,nz)
      real(kind=8) d3g13_rk1(nx,ny,nz)
      real(kind=8) d3g22_rk1(nx,ny,nz)
      real(kind=8) d3g23_rk1(nx,ny,nz)
      real(kind=8) d3g33_rk1(nx,ny,nz)
      real(kind=8) H0_rk1(nx,ny,nz)
      real(kind=8) H1_rk1(nx,ny,nz)
      real(kind=8) H2_rk1(nx,ny,nz)
      real(kind=8) H3_rk1(nx,ny,nz)
      real(kind=8) G0_rk1(nx,ny,nz)
      real(kind=8) d1H0_rk1(nx,ny,nz)
      real(kind=8) d2H0_rk1(nx,ny,nz)
      real(kind=8) d3H0_rk1(nx,ny,nz)
      real(kind=8) phim_rk1(nx,ny,nz)
      real(kind=8) pim_rk1(nx,ny,nz)
      real(kind=8) d1phim_rk1(nx,ny,nz)
      real(kind=8) d2phim_rk1(nx,ny,nz)
      real(kind=8) d3phim_rk1(nx,ny,nz)
      real(kind=8) phin_rk1(nx,ny,nz)
      real(kind=8) pin_rk1(nx,ny,nz)
      real(kind=8) d1phin_rk1(nx,ny,nz)
      real(kind=8) d2phin_rk1(nx,ny,nz)
      real(kind=8) d3phin_rk1(nx,ny,nz)
      real(kind=8) phir_rk1(nx,ny,nz)
      real(kind=8) pir_rk1(nx,ny,nz)
      real(kind=8) d1phir_rk1(nx,ny,nz)
      real(kind=8) d2phir_rk1(nx,ny,nz)
      real(kind=8) d3phir_rk1(nx,ny,nz)
      real(kind=8) phic_rk1(nx,ny,nz)
      real(kind=8) pic_rk1(nx,ny,nz)
      real(kind=8) d1phic_rk1(nx,ny,nz)
      real(kind=8) d2phic_rk1(nx,ny,nz)
      real(kind=8) d3phic_rk1(nx,ny,nz)
        target :: g00_rk1
        target :: g01_rk1
        target :: g02_rk1
        target :: g03_rk1
        target :: g11_rk1
        target :: g12_rk1
        target :: g13_rk1
        target :: g22_rk1
        target :: g23_rk1
        target :: g33_rk1
        target :: K00_rk1
        target :: K01_rk1
        target :: K02_rk1
        target :: K03_rk1
        target :: K11_rk1
        target :: K12_rk1
        target :: K13_rk1
        target :: K22_rk1
        target :: K23_rk1
        target :: K33_rk1
        target :: d1g00_rk1
        target :: d1g01_rk1
        target :: d1g02_rk1
        target :: d1g03_rk1
        target :: d1g11_rk1
        target :: d1g12_rk1
        target :: d1g13_rk1
        target :: d1g22_rk1
        target :: d1g23_rk1
        target :: d1g33_rk1
        target :: d2g00_rk1
        target :: d2g01_rk1
        target :: d2g02_rk1
        target :: d2g03_rk1
        target :: d2g11_rk1
        target :: d2g12_rk1
        target :: d2g13_rk1
        target :: d2g22_rk1
        target :: d2g23_rk1
        target :: d2g33_rk1
        target :: d3g00_rk1
        target :: d3g01_rk1
        target :: d3g02_rk1
        target :: d3g03_rk1
        target :: d3g11_rk1
        target :: d3g12_rk1
        target :: d3g13_rk1
        target :: d3g22_rk1
        target :: d3g23_rk1
        target :: d3g33_rk1
        target :: H0_rk1
        target :: H1_rk1
        target :: H2_rk1
        target :: H3_rk1
        target :: G0_rk1
        target :: d1H0_rk1
        target :: d2H0_rk1
        target :: d3H0_rk1
        target :: phim_rk1
        target :: pim_rk1
        target :: d1phim_rk1
        target :: d2phim_rk1
        target :: d3phim_rk1
        target :: phin_rk1
        target :: pin_rk1
        target :: d1phin_rk1
        target :: d2phin_rk1
        target :: d3phin_rk1
        target :: phir_rk1
        target :: pir_rk1
        target :: d1phir_rk1
        target :: d2phir_rk1
        target :: d3phir_rk1
        target :: phic_rk1
        target :: pic_rk1
        target :: d1phic_rk1
        target :: d2phic_rk1
        target :: d3phic_rk1
        urk1(H_G00)%d => g00_rk1
        urk1(H_G01)%d => g01_rk1
        urk1(H_G02)%d => g02_rk1
        urk1(H_G03)%d => g03_rk1
        urk1(H_G11)%d => g11_rk1
        urk1(H_G12)%d => g12_rk1
        urk1(H_G13)%d => g13_rk1
        urk1(H_G22)%d => g22_rk1
        urk1(H_G23)%d => g23_rk1
        urk1(H_G33)%d => g33_rk1
        urk1(H_K00)%d => K00_rk1
        urk1(H_K01)%d => K01_rk1
        urk1(H_K02)%d => K02_rk1
        urk1(H_K03)%d => K03_rk1
        urk1(H_K11)%d => K11_rk1
        urk1(H_K12)%d => K12_rk1
        urk1(H_K13)%d => K13_rk1
        urk1(H_K22)%d => K22_rk1
        urk1(H_K23)%d => K23_rk1
        urk1(H_K33)%d => K33_rk1
        urk1(H_D1G00)%d => d1g00_rk1
        urk1(H_D1G01)%d => d1g01_rk1
        urk1(H_D1G02)%d => d1g02_rk1
        urk1(H_D1G03)%d => d1g03_rk1
        urk1(H_D1G11)%d => d1g11_rk1
        urk1(H_D1G12)%d => d1g12_rk1
        urk1(H_D1G13)%d => d1g13_rk1
        urk1(H_D1G22)%d => d1g22_rk1
        urk1(H_D1G23)%d => d1g23_rk1
        urk1(H_D1G33)%d => d1g33_rk1
        urk1(H_D2G00)%d => d2g00_rk1
        urk1(H_D2G01)%d => d2g01_rk1
        urk1(H_D2G02)%d => d2g02_rk1
        urk1(H_D2G03)%d => d2g03_rk1
        urk1(H_D2G11)%d => d2g11_rk1
        urk1(H_D2G12)%d => d2g12_rk1
        urk1(H_D2G13)%d => d2g13_rk1
        urk1(H_D2G22)%d => d2g22_rk1
        urk1(H_D2G23)%d => d2g23_rk1
        urk1(H_D2G33)%d => d2g33_rk1
        urk1(H_D3G00)%d => d3g00_rk1
        urk1(H_D3G01)%d => d3g01_rk1
        urk1(H_D3G02)%d => d3g02_rk1
        urk1(H_D3G03)%d => d3g03_rk1
        urk1(H_D3G11)%d => d3g11_rk1
        urk1(H_D3G12)%d => d3g12_rk1
        urk1(H_D3G13)%d => d3g13_rk1
        urk1(H_D3G22)%d => d3g22_rk1
        urk1(H_D3G23)%d => d3g23_rk1
        urk1(H_D3G33)%d => d3g33_rk1
        urk1(H_H0)%d => H0_rk1
        urk1(H_H1)%d => H1_rk1
        urk1(H_H2)%d => H2_rk1
        urk1(H_H3)%d => H3_rk1
        urk1(H_G0)%d => G0_rk1
        urk1(H_D1H0)%d => d1H0_rk1
        urk1(H_D2H0)%d => d2H0_rk1
        urk1(H_D3H0)%d => d3H0_rk1
        urk1(H_PHIM)%d => phim_rk1
        urk1(H_PIM)%d => pim_rk1
        urk1(H_D1PHIM)%d => d1phim_rk1
        urk1(H_D2PHIM)%d => d2phim_rk1
        urk1(H_D3PHIM)%d => d3phim_rk1
        urk1(H_PHIN)%d => phin_rk1
        urk1(H_PIN)%d => pin_rk1
        urk1(H_D1PHIN)%d => d1phin_rk1
        urk1(H_D2PHIN)%d => d2phin_rk1
        urk1(H_D3PHIN)%d => d3phin_rk1
        urk1(H_PHIR)%d => phir_rk1
        urk1(H_PIR)%d => pir_rk1
        urk1(H_D1PHIR)%d => d1phir_rk1
        urk1(H_D2PHIR)%d => d2phir_rk1
        urk1(H_D3PHIR)%d => d3phir_rk1
        urk1(H_PHIC)%d => phic_rk1
        urk1(H_PIC)%d => pic_rk1
        urk1(H_D1PHIC)%d => d1phic_rk1
        urk1(H_D2PHIC)%d => d2phic_rk1
        urk1(H_D3PHIC)%d => d3phic_rk1
      return
      end subroutine assign_ptrs_rks
      subroutine assign_ptrs_derivs(dx_g00,dx_g01,dx_g02,dx_g03,dx_g11,d
     &x_g12,dx_g13,dx_g22,dx_g23,dx_g33,dx_K00,dx_K01,dx_K02,dx_K03,dx_K
     &11,dx_K12,dx_K13,dx_K22,dx_K23,dx_K33,dx_d1g00,dx_d1g01,dx_d1g02,d
     &x_d1g03,dx_d1g11,dx_d1g12,dx_d1g13,dx_d1g22,dx_d1g23,dx_d1g33,dx_d
     &2g00,dx_d2g01,dx_d2g02,dx_d2g03,dx_d2g11,dx_d2g12,dx_d2g13,dx_d2g2
     &2,dx_d2g23,dx_d2g33,dx_d3g00,dx_d3g01,dx_d3g02,dx_d3g03,dx_d3g11,d
     &x_d3g12,dx_d3g13,dx_d3g22,dx_d3g23,dx_d3g33,dx_H0,dx_H1,dx_H2,dx_H
     &3,dx_G0,dx_d1H0,dx_d2H0,dx_d3H0,dy_g00,dy_g01,dy_g02,dy_g03,dy_g11
     &,dy_g12,dy_g13,dy_g22,dy_g23,dy_g33,dy_K00,dy_K01,dy_K02,dy_K03,dy
     &_K11,dy_K12,dy_K13,dy_K22,dy_K23,dy_K33,dy_d1g00,dy_d1g01,dy_d1g02
     &,dy_d1g03,dy_d1g11,dy_d1g12,dy_d1g13,dy_d1g22,dy_d1g23,dy_d1g33,dy
     &_d2g00,dy_d2g01,dy_d2g02,dy_d2g03,dy_d2g11,dy_d2g12,dy_d2g13,dy_d2
     &g22,dy_d2g23,dy_d2g33,dy_d3g00,dy_d3g01,dy_d3g02,dy_d3g03,dy_d3g11
     &,dy_d3g12,dy_d3g13,dy_d3g22,dy_d3g23,dy_d3g33,dy_H0,dy_H1,dy_H2,dy
     &_H3,dy_G0,dy_d1H0,dy_d2H0,dy_d3H0,dz_g00,dz_g01,dz_g02,dz_g03,dz_g
     &11,dz_g12,dz_g13,dz_g22,dz_g23,dz_g33,dz_K00,dz_K01,dz_K02,dz_K03,
     &dz_K11,dz_K12,dz_K13,dz_K22,dz_K23,dz_K33,dz_d1g00,dz_d1g01,dz_d1g
     &02,dz_d1g03,dz_d1g11,dz_d1g12,dz_d1g13,dz_d1g22,dz_d1g23,dz_d1g33,
     &dz_d2g00,dz_d2g01,dz_d2g02,dz_d2g03,dz_d2g11,dz_d2g12,dz_d2g13,dz_
     &d2g22,dz_d2g23,dz_d2g33,dz_d3g00,dz_d3g01,dz_d3g02,dz_d3g03,dz_d3g
     &11,dz_d3g12,dz_d3g13,dz_d3g22,dz_d3g23,dz_d3g33,dz_H0,dz_H1,dz_H2,
     &dz_H3,dz_G0,dz_d1H0,dz_d2H0,dz_d3H0,dx_phim,dx_phin,dx_phir,dx_phi
     &c,dx_pim,dx_pin,dx_pir,dx_pic,dx_d1phim,dx_d1phin,dx_d1phir,dx_d1p
     &hic,dx_d2phim,dx_d2phin,dx_d2phir,dx_d2phic,dx_d3phim,dx_d3phin,dx
     &_d3phir,dx_d3phic,dy_phim,dy_phin,dy_phir,dy_phic,dy_pim,dy_pin,dy
     &_pir,dy_pic,dy_d1phim,dy_d1phin,dy_d1phir,dy_d1phic,dy_d2phim,dy_d
     &2phin,dy_d2phir,dy_d2phic,dy_d3phim,dy_d3phin,dy_d3phir,dy_d3phic,
     &dz_phim,dz_phin,dz_phir,dz_phic,dz_pim,dz_pin,dz_pir,dz_pic,dz_d1p
     &him,dz_d1phin,dz_d1phir,dz_d1phic,dz_d2phim,dz_d2phin,dz_d2phir,dz
     &_d2phic,dz_d3phim,dz_d3phin,dz_d3phir,dz_d3phic,dxu, dyu, dzu, dxv
     &, dyv, dzv, nx, ny, nz)
      use GF
      implicit none
      integer :: nx, ny, nz
      type(gridfunction), dimension(NU) :: dxu, dyu, dzu
      type(gridfunction), dimension(NV) :: dxv, dyv, dzv
      real(kind=8) dx_g00(nx,ny,nz)
      real(kind=8) dx_g01(nx,ny,nz)
      real(kind=8) dx_g02(nx,ny,nz)
      real(kind=8) dx_g03(nx,ny,nz)
      real(kind=8) dx_g11(nx,ny,nz)
      real(kind=8) dx_g12(nx,ny,nz)
      real(kind=8) dx_g13(nx,ny,nz)
      real(kind=8) dx_g22(nx,ny,nz)
      real(kind=8) dx_g23(nx,ny,nz)
      real(kind=8) dx_g33(nx,ny,nz)
      real(kind=8) dx_K00(nx,ny,nz)
      real(kind=8) dx_K01(nx,ny,nz)
      real(kind=8) dx_K02(nx,ny,nz)
      real(kind=8) dx_K03(nx,ny,nz)
      real(kind=8) dx_K11(nx,ny,nz)
      real(kind=8) dx_K12(nx,ny,nz)
      real(kind=8) dx_K13(nx,ny,nz)
      real(kind=8) dx_K22(nx,ny,nz)
      real(kind=8) dx_K23(nx,ny,nz)
      real(kind=8) dx_K33(nx,ny,nz)
      real(kind=8) dx_d1g00(nx,ny,nz)
      real(kind=8) dx_d1g01(nx,ny,nz)
      real(kind=8) dx_d1g02(nx,ny,nz)
      real(kind=8) dx_d1g03(nx,ny,nz)
      real(kind=8) dx_d1g11(nx,ny,nz)
      real(kind=8) dx_d1g12(nx,ny,nz)
      real(kind=8) dx_d1g13(nx,ny,nz)
      real(kind=8) dx_d1g22(nx,ny,nz)
      real(kind=8) dx_d1g23(nx,ny,nz)
      real(kind=8) dx_d1g33(nx,ny,nz)
      real(kind=8) dx_d2g00(nx,ny,nz)
      real(kind=8) dx_d2g01(nx,ny,nz)
      real(kind=8) dx_d2g02(nx,ny,nz)
      real(kind=8) dx_d2g03(nx,ny,nz)
      real(kind=8) dx_d2g11(nx,ny,nz)
      real(kind=8) dx_d2g12(nx,ny,nz)
      real(kind=8) dx_d2g13(nx,ny,nz)
      real(kind=8) dx_d2g22(nx,ny,nz)
      real(kind=8) dx_d2g23(nx,ny,nz)
      real(kind=8) dx_d2g33(nx,ny,nz)
      real(kind=8) dx_d3g00(nx,ny,nz)
      real(kind=8) dx_d3g01(nx,ny,nz)
      real(kind=8) dx_d3g02(nx,ny,nz)
      real(kind=8) dx_d3g03(nx,ny,nz)
      real(kind=8) dx_d3g11(nx,ny,nz)
      real(kind=8) dx_d3g12(nx,ny,nz)
      real(kind=8) dx_d3g13(nx,ny,nz)
      real(kind=8) dx_d3g22(nx,ny,nz)
      real(kind=8) dx_d3g23(nx,ny,nz)
      real(kind=8) dx_d3g33(nx,ny,nz)
      real(kind=8) dx_H0(nx,ny,nz)
      real(kind=8) dx_H1(nx,ny,nz)
      real(kind=8) dx_H2(nx,ny,nz)
      real(kind=8) dx_H3(nx,ny,nz)
      real(kind=8) dx_G0(nx,ny,nz)
      real(kind=8) dx_d1H0(nx,ny,nz)
      real(kind=8) dx_d2H0(nx,ny,nz)
      real(kind=8) dx_d3H0(nx,ny,nz)
      real(kind=8) dy_g00(nx,ny,nz)
      real(kind=8) dy_g01(nx,ny,nz)
      real(kind=8) dy_g02(nx,ny,nz)
      real(kind=8) dy_g03(nx,ny,nz)
      real(kind=8) dy_g11(nx,ny,nz)
      real(kind=8) dy_g12(nx,ny,nz)
      real(kind=8) dy_g13(nx,ny,nz)
      real(kind=8) dy_g22(nx,ny,nz)
      real(kind=8) dy_g23(nx,ny,nz)
      real(kind=8) dy_g33(nx,ny,nz)
      real(kind=8) dy_K00(nx,ny,nz)
      real(kind=8) dy_K01(nx,ny,nz)
      real(kind=8) dy_K02(nx,ny,nz)
      real(kind=8) dy_K03(nx,ny,nz)
      real(kind=8) dy_K11(nx,ny,nz)
      real(kind=8) dy_K12(nx,ny,nz)
      real(kind=8) dy_K13(nx,ny,nz)
      real(kind=8) dy_K22(nx,ny,nz)
      real(kind=8) dy_K23(nx,ny,nz)
      real(kind=8) dy_K33(nx,ny,nz)
      real(kind=8) dy_d1g00(nx,ny,nz)
      real(kind=8) dy_d1g01(nx,ny,nz)
      real(kind=8) dy_d1g02(nx,ny,nz)
      real(kind=8) dy_d1g03(nx,ny,nz)
      real(kind=8) dy_d1g11(nx,ny,nz)
      real(kind=8) dy_d1g12(nx,ny,nz)
      real(kind=8) dy_d1g13(nx,ny,nz)
      real(kind=8) dy_d1g22(nx,ny,nz)
      real(kind=8) dy_d1g23(nx,ny,nz)
      real(kind=8) dy_d1g33(nx,ny,nz)
      real(kind=8) dy_d2g00(nx,ny,nz)
      real(kind=8) dy_d2g01(nx,ny,nz)
      real(kind=8) dy_d2g02(nx,ny,nz)
      real(kind=8) dy_d2g03(nx,ny,nz)
      real(kind=8) dy_d2g11(nx,ny,nz)
      real(kind=8) dy_d2g12(nx,ny,nz)
      real(kind=8) dy_d2g13(nx,ny,nz)
      real(kind=8) dy_d2g22(nx,ny,nz)
      real(kind=8) dy_d2g23(nx,ny,nz)
      real(kind=8) dy_d2g33(nx,ny,nz)
      real(kind=8) dy_d3g00(nx,ny,nz)
      real(kind=8) dy_d3g01(nx,ny,nz)
      real(kind=8) dy_d3g02(nx,ny,nz)
      real(kind=8) dy_d3g03(nx,ny,nz)
      real(kind=8) dy_d3g11(nx,ny,nz)
      real(kind=8) dy_d3g12(nx,ny,nz)
      real(kind=8) dy_d3g13(nx,ny,nz)
      real(kind=8) dy_d3g22(nx,ny,nz)
      real(kind=8) dy_d3g23(nx,ny,nz)
      real(kind=8) dy_d3g33(nx,ny,nz)
      real(kind=8) dy_H0(nx,ny,nz)
      real(kind=8) dy_H1(nx,ny,nz)
      real(kind=8) dy_H2(nx,ny,nz)
      real(kind=8) dy_H3(nx,ny,nz)
      real(kind=8) dy_G0(nx,ny,nz)
      real(kind=8) dy_d1H0(nx,ny,nz)
      real(kind=8) dy_d2H0(nx,ny,nz)
      real(kind=8) dy_d3H0(nx,ny,nz)
      real(kind=8) dz_g00(nx,ny,nz)
      real(kind=8) dz_g01(nx,ny,nz)
      real(kind=8) dz_g02(nx,ny,nz)
      real(kind=8) dz_g03(nx,ny,nz)
      real(kind=8) dz_g11(nx,ny,nz)
      real(kind=8) dz_g12(nx,ny,nz)
      real(kind=8) dz_g13(nx,ny,nz)
      real(kind=8) dz_g22(nx,ny,nz)
      real(kind=8) dz_g23(nx,ny,nz)
      real(kind=8) dz_g33(nx,ny,nz)
      real(kind=8) dz_K00(nx,ny,nz)
      real(kind=8) dz_K01(nx,ny,nz)
      real(kind=8) dz_K02(nx,ny,nz)
      real(kind=8) dz_K03(nx,ny,nz)
      real(kind=8) dz_K11(nx,ny,nz)
      real(kind=8) dz_K12(nx,ny,nz)
      real(kind=8) dz_K13(nx,ny,nz)
      real(kind=8) dz_K22(nx,ny,nz)
      real(kind=8) dz_K23(nx,ny,nz)
      real(kind=8) dz_K33(nx,ny,nz)
      real(kind=8) dz_d1g00(nx,ny,nz)
      real(kind=8) dz_d1g01(nx,ny,nz)
      real(kind=8) dz_d1g02(nx,ny,nz)
      real(kind=8) dz_d1g03(nx,ny,nz)
      real(kind=8) dz_d1g11(nx,ny,nz)
      real(kind=8) dz_d1g12(nx,ny,nz)
      real(kind=8) dz_d1g13(nx,ny,nz)
      real(kind=8) dz_d1g22(nx,ny,nz)
      real(kind=8) dz_d1g23(nx,ny,nz)
      real(kind=8) dz_d1g33(nx,ny,nz)
      real(kind=8) dz_d2g00(nx,ny,nz)
      real(kind=8) dz_d2g01(nx,ny,nz)
      real(kind=8) dz_d2g02(nx,ny,nz)
      real(kind=8) dz_d2g03(nx,ny,nz)
      real(kind=8) dz_d2g11(nx,ny,nz)
      real(kind=8) dz_d2g12(nx,ny,nz)
      real(kind=8) dz_d2g13(nx,ny,nz)
      real(kind=8) dz_d2g22(nx,ny,nz)
      real(kind=8) dz_d2g23(nx,ny,nz)
      real(kind=8) dz_d2g33(nx,ny,nz)
      real(kind=8) dz_d3g00(nx,ny,nz)
      real(kind=8) dz_d3g01(nx,ny,nz)
      real(kind=8) dz_d3g02(nx,ny,nz)
      real(kind=8) dz_d3g03(nx,ny,nz)
      real(kind=8) dz_d3g11(nx,ny,nz)
      real(kind=8) dz_d3g12(nx,ny,nz)
      real(kind=8) dz_d3g13(nx,ny,nz)
      real(kind=8) dz_d3g22(nx,ny,nz)
      real(kind=8) dz_d3g23(nx,ny,nz)
      real(kind=8) dz_d3g33(nx,ny,nz)
      real(kind=8) dz_H0(nx,ny,nz)
      real(kind=8) dz_H1(nx,ny,nz)
      real(kind=8) dz_H2(nx,ny,nz)
      real(kind=8) dz_H3(nx,ny,nz)
      real(kind=8) dz_G0(nx,ny,nz)
      real(kind=8) dz_d1H0(nx,ny,nz)
      real(kind=8) dz_d2H0(nx,ny,nz)
      real(kind=8) dz_d3H0(nx,ny,nz)
      real(kind=8) dx_phim(nx,ny,nz)
      real(kind=8) dx_phin(nx,ny,nz)
      real(kind=8) dx_phir(nx,ny,nz)
      real(kind=8) dx_phic(nx,ny,nz)
      real(kind=8) dx_pim(nx,ny,nz)
      real(kind=8) dx_pin(nx,ny,nz)
      real(kind=8) dx_pir(nx,ny,nz)
      real(kind=8) dx_pic(nx,ny,nz)
      real(kind=8) dx_d1phim(nx,ny,nz)
      real(kind=8) dx_d1phin(nx,ny,nz)
      real(kind=8) dx_d1phir(nx,ny,nz)
      real(kind=8) dx_d1phic(nx,ny,nz)
      real(kind=8) dx_d2phim(nx,ny,nz)
      real(kind=8) dx_d2phin(nx,ny,nz)
      real(kind=8) dx_d2phir(nx,ny,nz)
      real(kind=8) dx_d2phic(nx,ny,nz)
      real(kind=8) dx_d3phim(nx,ny,nz)
      real(kind=8) dx_d3phin(nx,ny,nz)
      real(kind=8) dx_d3phir(nx,ny,nz)
      real(kind=8) dx_d3phic(nx,ny,nz)
      real(kind=8) dy_phim(nx,ny,nz)
      real(kind=8) dy_phin(nx,ny,nz)
      real(kind=8) dy_phir(nx,ny,nz)
      real(kind=8) dy_phic(nx,ny,nz)
      real(kind=8) dy_pim(nx,ny,nz)
      real(kind=8) dy_pin(nx,ny,nz)
      real(kind=8) dy_pir(nx,ny,nz)
      real(kind=8) dy_pic(nx,ny,nz)
      real(kind=8) dy_d1phim(nx,ny,nz)
      real(kind=8) dy_d1phin(nx,ny,nz)
      real(kind=8) dy_d1phir(nx,ny,nz)
      real(kind=8) dy_d1phic(nx,ny,nz)
      real(kind=8) dy_d2phim(nx,ny,nz)
      real(kind=8) dy_d2phin(nx,ny,nz)
      real(kind=8) dy_d2phir(nx,ny,nz)
      real(kind=8) dy_d2phic(nx,ny,nz)
      real(kind=8) dy_d3phim(nx,ny,nz)
      real(kind=8) dy_d3phin(nx,ny,nz)
      real(kind=8) dy_d3phir(nx,ny,nz)
      real(kind=8) dy_d3phic(nx,ny,nz)
      real(kind=8) dz_phim(nx,ny,nz)
      real(kind=8) dz_phin(nx,ny,nz)
      real(kind=8) dz_phir(nx,ny,nz)
      real(kind=8) dz_phic(nx,ny,nz)
      real(kind=8) dz_pim(nx,ny,nz)
      real(kind=8) dz_pin(nx,ny,nz)
      real(kind=8) dz_pir(nx,ny,nz)
      real(kind=8) dz_pic(nx,ny,nz)
      real(kind=8) dz_d1phim(nx,ny,nz)
      real(kind=8) dz_d1phin(nx,ny,nz)
      real(kind=8) dz_d1phir(nx,ny,nz)
      real(kind=8) dz_d1phic(nx,ny,nz)
      real(kind=8) dz_d2phim(nx,ny,nz)
      real(kind=8) dz_d2phin(nx,ny,nz)
      real(kind=8) dz_d2phir(nx,ny,nz)
      real(kind=8) dz_d2phic(nx,ny,nz)
      real(kind=8) dz_d3phim(nx,ny,nz)
      real(kind=8) dz_d3phin(nx,ny,nz)
      real(kind=8) dz_d3phir(nx,ny,nz)
      real(kind=8) dz_d3phic(nx,ny,nz)
      target :: dx_g00
      target :: dx_g01
      target :: dx_g02
      target :: dx_g03
      target :: dx_g11
      target :: dx_g12
      target :: dx_g13
      target :: dx_g22
      target :: dx_g23
      target :: dx_g33
      target :: dx_K00
      target :: dx_K01
      target :: dx_K02
      target :: dx_K03
      target :: dx_K11
      target :: dx_K12
      target :: dx_K13
      target :: dx_K22
      target :: dx_K23
      target :: dx_K33
      target :: dx_d1g00
      target :: dx_d1g01
      target :: dx_d1g02
      target :: dx_d1g03
      target :: dx_d1g11
      target :: dx_d1g12
      target :: dx_d1g13
      target :: dx_d1g22
      target :: dx_d1g23
      target :: dx_d1g33
      target :: dx_d2g00
      target :: dx_d2g01
      target :: dx_d2g02
      target :: dx_d2g03
      target :: dx_d2g11
      target :: dx_d2g12
      target :: dx_d2g13
      target :: dx_d2g22
      target :: dx_d2g23
      target :: dx_d2g33
      target :: dx_d3g00
      target :: dx_d3g01
      target :: dx_d3g02
      target :: dx_d3g03
      target :: dx_d3g11
      target :: dx_d3g12
      target :: dx_d3g13
      target :: dx_d3g22
      target :: dx_d3g23
      target :: dx_d3g33
      target :: dx_H0
      target :: dx_H1
      target :: dx_H2
      target :: dx_H3
      target :: dx_G0
      target :: dx_d1H0
      target :: dx_d2H0
      target :: dx_d3H0
      target :: dy_g00
      target :: dy_g01
      target :: dy_g02
      target :: dy_g03
      target :: dy_g11
      target :: dy_g12
      target :: dy_g13
      target :: dy_g22
      target :: dy_g23
      target :: dy_g33
      target :: dy_K00
      target :: dy_K01
      target :: dy_K02
      target :: dy_K03
      target :: dy_K11
      target :: dy_K12
      target :: dy_K13
      target :: dy_K22
      target :: dy_K23
      target :: dy_K33
      target :: dy_d1g00
      target :: dy_d1g01
      target :: dy_d1g02
      target :: dy_d1g03
      target :: dy_d1g11
      target :: dy_d1g12
      target :: dy_d1g13
      target :: dy_d1g22
      target :: dy_d1g23
      target :: dy_d1g33
      target :: dy_d2g00
      target :: dy_d2g01
      target :: dy_d2g02
      target :: dy_d2g03
      target :: dy_d2g11
      target :: dy_d2g12
      target :: dy_d2g13
      target :: dy_d2g22
      target :: dy_d2g23
      target :: dy_d2g33
      target :: dy_d3g00
      target :: dy_d3g01
      target :: dy_d3g02
      target :: dy_d3g03
      target :: dy_d3g11
      target :: dy_d3g12
      target :: dy_d3g13
      target :: dy_d3g22
      target :: dy_d3g23
      target :: dy_d3g33
      target :: dy_H0
      target :: dy_H1
      target :: dy_H2
      target :: dy_H3
      target :: dy_G0
      target :: dy_d1H0
      target :: dy_d2H0
      target :: dy_d3H0
      target :: dz_g00
      target :: dz_g01
      target :: dz_g02
      target :: dz_g03
      target :: dz_g11
      target :: dz_g12
      target :: dz_g13
      target :: dz_g22
      target :: dz_g23
      target :: dz_g33
      target :: dz_K00
      target :: dz_K01
      target :: dz_K02
      target :: dz_K03
      target :: dz_K11
      target :: dz_K12
      target :: dz_K13
      target :: dz_K22
      target :: dz_K23
      target :: dz_K33
      target :: dz_d1g00
      target :: dz_d1g01
      target :: dz_d1g02
      target :: dz_d1g03
      target :: dz_d1g11
      target :: dz_d1g12
      target :: dz_d1g13
      target :: dz_d1g22
      target :: dz_d1g23
      target :: dz_d1g33
      target :: dz_d2g00
      target :: dz_d2g01
      target :: dz_d2g02
      target :: dz_d2g03
      target :: dz_d2g11
      target :: dz_d2g12
      target :: dz_d2g13
      target :: dz_d2g22
      target :: dz_d2g23
      target :: dz_d2g33
      target :: dz_d3g00
      target :: dz_d3g01
      target :: dz_d3g02
      target :: dz_d3g03
      target :: dz_d3g11
      target :: dz_d3g12
      target :: dz_d3g13
      target :: dz_d3g22
      target :: dz_d3g23
      target :: dz_d3g33
      target :: dz_H0
      target :: dz_H1
      target :: dz_H2
      target :: dz_H3
      target :: dz_G0
      target :: dz_d1H0
      target :: dz_d2H0
      target :: dz_d3H0
      target :: dx_phim
      target :: dx_phin
      target :: dx_phir
      target :: dx_phic
      target :: dx_pim
      target :: dx_pin
      target :: dx_pir
      target :: dx_pic
      target :: dx_d1phim
      target :: dx_d1phin
      target :: dx_d1phir
      target :: dx_d1phic
      target :: dx_d2phim
      target :: dx_d2phin
      target :: dx_d2phir
      target :: dx_d2phic
      target :: dx_d3phim
      target :: dx_d3phin
      target :: dx_d3phir
      target :: dx_d3phic
      target :: dy_phim
      target :: dy_phin
      target :: dy_phir
      target :: dy_phic
      target :: dy_pim
      target :: dy_pin
      target :: dy_pir
      target :: dy_pic
      target :: dy_d1phim
      target :: dy_d1phin
      target :: dy_d1phir
      target :: dy_d1phic
      target :: dy_d2phim
      target :: dy_d2phin
      target :: dy_d2phir
      target :: dy_d2phic
      target :: dy_d3phim
      target :: dy_d3phin
      target :: dy_d3phir
      target :: dy_d3phic
      target :: dz_phim
      target :: dz_phin
      target :: dz_phir
      target :: dz_phic
      target :: dz_pim
      target :: dz_pin
      target :: dz_pir
      target :: dz_pic
      target :: dz_d1phim
      target :: dz_d1phin
      target :: dz_d1phir
      target :: dz_d1phic
      target :: dz_d2phim
      target :: dz_d2phin
      target :: dz_d2phir
      target :: dz_d2phic
      target :: dz_d3phim
      target :: dz_d3phin
      target :: dz_d3phir
      target :: dz_d3phic
      dxu(H_G00)%d => dx_g00
      dyu(H_G00)%d => dy_g00
      dzu(H_G00)%d => dz_g00
      dxu(H_G01)%d => dx_g01
      dyu(H_G01)%d => dy_g01
      dzu(H_G01)%d => dz_g01
      dxu(H_G02)%d => dx_g02
      dyu(H_G02)%d => dy_g02
      dzu(H_G02)%d => dz_g02
      dxu(H_G03)%d => dx_g03
      dyu(H_G03)%d => dy_g03
      dzu(H_G03)%d => dz_g03
      dxu(H_G11)%d => dx_g11
      dyu(H_G11)%d => dy_g11
      dzu(H_G11)%d => dz_g11
      dxu(H_G12)%d => dx_g12
      dyu(H_G12)%d => dy_g12
      dzu(H_G12)%d => dz_g12
      dxu(H_G13)%d => dx_g13
      dyu(H_G13)%d => dy_g13
      dzu(H_G13)%d => dz_g13
      dxu(H_G22)%d => dx_g22
      dyu(H_G22)%d => dy_g22
      dzu(H_G22)%d => dz_g22
      dxu(H_G23)%d => dx_g23
      dyu(H_G23)%d => dy_g23
      dzu(H_G23)%d => dz_g23
      dxu(H_G33)%d => dx_g33
      dyu(H_G33)%d => dy_g33
      dzu(H_G33)%d => dz_g33
      dxu(H_K00)%d => dx_K00
      dyu(H_K00)%d => dy_K00
      dzu(H_K00)%d => dz_K00
      dxu(H_K01)%d => dx_K01
      dyu(H_K01)%d => dy_K01
      dzu(H_K01)%d => dz_K01
      dxu(H_K02)%d => dx_K02
      dyu(H_K02)%d => dy_K02
      dzu(H_K02)%d => dz_K02
      dxu(H_K03)%d => dx_K03
      dyu(H_K03)%d => dy_K03
      dzu(H_K03)%d => dz_K03
      dxu(H_K11)%d => dx_K11
      dyu(H_K11)%d => dy_K11
      dzu(H_K11)%d => dz_K11
      dxu(H_K12)%d => dx_K12
      dyu(H_K12)%d => dy_K12
      dzu(H_K12)%d => dz_K12
      dxu(H_K13)%d => dx_K13
      dyu(H_K13)%d => dy_K13
      dzu(H_K13)%d => dz_K13
      dxu(H_K22)%d => dx_K22
      dyu(H_K22)%d => dy_K22
      dzu(H_K22)%d => dz_K22
      dxu(H_K23)%d => dx_K23
      dyu(H_K23)%d => dy_K23
      dzu(H_K23)%d => dz_K23
      dxu(H_K33)%d => dx_K33
      dyu(H_K33)%d => dy_K33
      dzu(H_K33)%d => dz_K33
      dxu(H_D1G00)%d => dx_d1g00
      dyu(H_D1G00)%d => dy_d1g00
      dzu(H_D1G00)%d => dz_d1g00
      dxu(H_D1G01)%d => dx_d1g01
      dyu(H_D1G01)%d => dy_d1g01
      dzu(H_D1G01)%d => dz_d1g01
      dxu(H_D1G02)%d => dx_d1g02
      dyu(H_D1G02)%d => dy_d1g02
      dzu(H_D1G02)%d => dz_d1g02
      dxu(H_D1G03)%d => dx_d1g03
      dyu(H_D1G03)%d => dy_d1g03
      dzu(H_D1G03)%d => dz_d1g03
      dxu(H_D1G11)%d => dx_d1g11
      dyu(H_D1G11)%d => dy_d1g11
      dzu(H_D1G11)%d => dz_d1g11
      dxu(H_D1G12)%d => dx_d1g12
      dyu(H_D1G12)%d => dy_d1g12
      dzu(H_D1G12)%d => dz_d1g12
      dxu(H_D1G13)%d => dx_d1g13
      dyu(H_D1G13)%d => dy_d1g13
      dzu(H_D1G13)%d => dz_d1g13
      dxu(H_D1G22)%d => dx_d1g22
      dyu(H_D1G22)%d => dy_d1g22
      dzu(H_D1G22)%d => dz_d1g22
      dxu(H_D1G23)%d => dx_d1g23
      dyu(H_D1G23)%d => dy_d1g23
      dzu(H_D1G23)%d => dz_d1g23
      dxu(H_D1G33)%d => dx_d1g33
      dyu(H_D1G33)%d => dy_d1g33
      dzu(H_D1G33)%d => dz_d1g33
      dxu(H_D2G00)%d => dx_d2g00
      dyu(H_D2G00)%d => dy_d2g00
      dzu(H_D2G00)%d => dz_d2g00
      dxu(H_D2G01)%d => dx_d2g01
      dyu(H_D2G01)%d => dy_d2g01
      dzu(H_D2G01)%d => dz_d2g01
      dxu(H_D2G02)%d => dx_d2g02
      dyu(H_D2G02)%d => dy_d2g02
      dzu(H_D2G02)%d => dz_d2g02
      dxu(H_D2G03)%d => dx_d2g03
      dyu(H_D2G03)%d => dy_d2g03
      dzu(H_D2G03)%d => dz_d2g03
      dxu(H_D2G11)%d => dx_d2g11
      dyu(H_D2G11)%d => dy_d2g11
      dzu(H_D2G11)%d => dz_d2g11
      dxu(H_D2G12)%d => dx_d2g12
      dyu(H_D2G12)%d => dy_d2g12
      dzu(H_D2G12)%d => dz_d2g12
      dxu(H_D2G13)%d => dx_d2g13
      dyu(H_D2G13)%d => dy_d2g13
      dzu(H_D2G13)%d => dz_d2g13
      dxu(H_D2G22)%d => dx_d2g22
      dyu(H_D2G22)%d => dy_d2g22
      dzu(H_D2G22)%d => dz_d2g22
      dxu(H_D2G23)%d => dx_d2g23
      dyu(H_D2G23)%d => dy_d2g23
      dzu(H_D2G23)%d => dz_d2g23
      dxu(H_D2G33)%d => dx_d2g33
      dyu(H_D2G33)%d => dy_d2g33
      dzu(H_D2G33)%d => dz_d2g33
      dxu(H_D3G00)%d => dx_d3g00
      dyu(H_D3G00)%d => dy_d3g00
      dzu(H_D3G00)%d => dz_d3g00
      dxu(H_D3G01)%d => dx_d3g01
      dyu(H_D3G01)%d => dy_d3g01
      dzu(H_D3G01)%d => dz_d3g01
      dxu(H_D3G02)%d => dx_d3g02
      dyu(H_D3G02)%d => dy_d3g02
      dzu(H_D3G02)%d => dz_d3g02
      dxu(H_D3G03)%d => dx_d3g03
      dyu(H_D3G03)%d => dy_d3g03
      dzu(H_D3G03)%d => dz_d3g03
      dxu(H_D3G11)%d => dx_d3g11
      dyu(H_D3G11)%d => dy_d3g11
      dzu(H_D3G11)%d => dz_d3g11
      dxu(H_D3G12)%d => dx_d3g12
      dyu(H_D3G12)%d => dy_d3g12
      dzu(H_D3G12)%d => dz_d3g12
      dxu(H_D3G13)%d => dx_d3g13
      dyu(H_D3G13)%d => dy_d3g13
      dzu(H_D3G13)%d => dz_d3g13
      dxu(H_D3G22)%d => dx_d3g22
      dyu(H_D3G22)%d => dy_d3g22
      dzu(H_D3G22)%d => dz_d3g22
      dxu(H_D3G23)%d => dx_d3g23
      dyu(H_D3G23)%d => dy_d3g23
      dzu(H_D3G23)%d => dz_d3g23
      dxu(H_D3G33)%d => dx_d3g33
      dyu(H_D3G33)%d => dy_d3g33
      dzu(H_D3G33)%d => dz_d3g33
      dxu(H_H0)%d => dx_H0
      dyu(H_H0)%d => dy_H0
      dzu(H_H0)%d => dz_H0
      dxu(H_H1)%d => dx_H1
      dyu(H_H1)%d => dy_H1
      dzu(H_H1)%d => dz_H1
      dxu(H_H2)%d => dx_H2
      dyu(H_H2)%d => dy_H2
      dzu(H_H2)%d => dz_H2
      dxu(H_H3)%d => dx_H3
      dyu(H_H3)%d => dy_H3
      dzu(H_H3)%d => dz_H3
      dxu(H_G0)%d => dx_G0
      dyu(H_G0)%d => dy_G0
      dzu(H_G0)%d => dz_G0
      dxu(H_D1H0)%d => dx_d1H0
      dyu(H_D1H0)%d => dy_d1H0
      dzu(H_D1H0)%d => dz_d1H0
      dxu(H_D2H0)%d => dx_d2H0
      dyu(H_D2H0)%d => dy_d2H0
      dzu(H_D2H0)%d => dz_d2H0
      dxu(H_D3H0)%d => dx_d3H0
      dyu(H_D3H0)%d => dy_d3H0
      dzu(H_D3H0)%d => dz_d3H0
      dxu(H_PHIM)%d => dx_phim
      dyu(H_PHIM)%d => dy_phim
      dzu(H_PHIM)%d => dz_phim
      dxu(H_PIM)%d => dx_pim
      dyu(H_PIM)%d => dy_pim
      dzu(H_PIM)%d => dz_pim
      dxu(H_D1PHIM)%d => dx_d1phim
      dyu(H_D1PHIM)%d => dy_d1phim
      dzu(H_D1PHIM)%d => dz_d1phim
      dxu(H_D2PHIM)%d => dx_d2phim
      dyu(H_D2PHIM)%d => dy_d2phim
      dzu(H_D2PHIM)%d => dz_d2phim
      dxu(H_D3PHIM)%d => dx_d3phim
      dyu(H_D3PHIM)%d => dy_d3phim
      dzu(H_D3PHIM)%d => dz_d3phim
      dxu(H_PHIN)%d => dx_phin
      dyu(H_PHIN)%d => dy_phin
      dzu(H_PHIN)%d => dz_phin
      dxu(H_PIN)%d => dx_pin
      dyu(H_PIN)%d => dy_pin
      dzu(H_PIN)%d => dz_pin
      dxu(H_D1PHIN)%d => dx_d1phin
      dyu(H_D1PHIN)%d => dy_d1phin
      dzu(H_D1PHIN)%d => dz_d1phin
      dxu(H_D2PHIN)%d => dx_d2phin
      dyu(H_D2PHIN)%d => dy_d2phin
      dzu(H_D2PHIN)%d => dz_d2phin
      dxu(H_D3PHIN)%d => dx_d3phin
      dyu(H_D3PHIN)%d => dy_d3phin
      dzu(H_D3PHIN)%d => dz_d3phin
      dxu(H_PHIR)%d => dx_phir
      dyu(H_PHIR)%d => dy_phir
      dzu(H_PHIR)%d => dz_phir
      dxu(H_PIR)%d => dx_pir
      dyu(H_PIR)%d => dy_pir
      dzu(H_PIR)%d => dz_pir
      dxu(H_D1PHIR)%d => dx_d1phir
      dyu(H_D1PHIR)%d => dy_d1phir
      dzu(H_D1PHIR)%d => dz_d1phir
      dxu(H_D2PHIR)%d => dx_d2phir
      dyu(H_D2PHIR)%d => dy_d2phir
      dzu(H_D2PHIR)%d => dz_d2phir
      dxu(H_D3PHIR)%d => dx_d3phir
      dyu(H_D3PHIR)%d => dy_d3phir
      dzu(H_D3PHIR)%d => dz_d3phir
      dxu(H_PHIC)%d => dx_phic
      dyu(H_PHIC)%d => dy_phic
      dzu(H_PHIC)%d => dz_phic
      dxu(H_PIC)%d => dx_pic
      dyu(H_PIC)%d => dy_pic
      dzu(H_PIC)%d => dz_pic
      dxu(H_D1PHIC)%d => dx_d1phic
      dyu(H_D1PHIC)%d => dy_d1phic
      dzu(H_D1PHIC)%d => dz_d1phic
      dxu(H_D2PHIC)%d => dx_d2phic
      dyu(H_D2PHIC)%d => dy_d2phic
      dzu(H_D2PHIC)%d => dz_d2phic
      dxu(H_D3PHIC)%d => dx_d3phic
      dyu(H_D3PHIC)%d => dy_d3phic
      dzu(H_D3PHIC)%d => dz_d3phic
      nullify(dxv(H_ALPHA)%d)
      nullify(dyv(H_ALPHA)%d)
      nullify(dzv(H_ALPHA)%d)
      nullify(dxv(H_SHIFT1)%d)
      nullify(dyv(H_SHIFT1)%d)
      nullify(dzv(H_SHIFT1)%d)
      nullify(dxv(H_SHIFT2)%d)
      nullify(dyv(H_SHIFT2)%d)
      nullify(dzv(H_SHIFT2)%d)
      nullify(dxv(H_SHIFT3)%d)
      nullify(dyv(H_SHIFT3)%d)
      nullify(dzv(H_SHIFT3)%d)
      nullify(dxv(H_SQDETG)%d)
      nullify(dyv(H_SQDETG)%d)
      nullify(dzv(H_SQDETG)%d)
      nullify(dxv(H_TMUNU00)%d)
      nullify(dyv(H_TMUNU00)%d)
      nullify(dzv(H_TMUNU00)%d)
      nullify(dxv(H_TMUNU01)%d)
      nullify(dyv(H_TMUNU01)%d)
      nullify(dzv(H_TMUNU01)%d)
      nullify(dxv(H_TMUNU02)%d)
      nullify(dyv(H_TMUNU02)%d)
      nullify(dzv(H_TMUNU02)%d)
      nullify(dxv(H_TMUNU03)%d)
      nullify(dyv(H_TMUNU03)%d)
      nullify(dzv(H_TMUNU03)%d)
      nullify(dxv(H_TMUNU11)%d)
      nullify(dyv(H_TMUNU11)%d)
      nullify(dzv(H_TMUNU11)%d)
      nullify(dxv(H_TMUNU12)%d)
      nullify(dyv(H_TMUNU12)%d)
      nullify(dzv(H_TMUNU12)%d)
      nullify(dxv(H_TMUNU13)%d)
      nullify(dyv(H_TMUNU13)%d)
      nullify(dzv(H_TMUNU13)%d)
      nullify(dxv(H_TMUNU22)%d)
      nullify(dyv(H_TMUNU22)%d)
      nullify(dzv(H_TMUNU22)%d)
      nullify(dxv(H_TMUNU23)%d)
      nullify(dyv(H_TMUNU23)%d)
      nullify(dzv(H_TMUNU23)%d)
      nullify(dxv(H_TMUNU33)%d)
      nullify(dyv(H_TMUNU33)%d)
      nullify(dzv(H_TMUNU33)%d)
      nullify(dxv(H_ORDER1_CONST)%d)
      nullify(dyv(H_ORDER1_CONST)%d)
      nullify(dzv(H_ORDER1_CONST)%d)
      nullify(dxv(H_ORDER2_CONST)%d)
      nullify(dyv(H_ORDER2_CONST)%d)
      nullify(dzv(H_ORDER2_CONST)%d)
      nullify(dxv(H_PHYSZ_CONST)%d)
      nullify(dyv(H_PHYSZ_CONST)%d)
      nullify(dzv(H_PHYSZ_CONST)%d)
      nullify(dxv(H_RAD_EXP)%d)
      nullify(dyv(H_RAD_EXP)%d)
      nullify(dzv(H_RAD_EXP)%d)
      nullify(dxv(H_RAD_SPEED)%d)
      nullify(dyv(H_RAD_SPEED)%d)
      nullify(dzv(H_RAD_SPEED)%d)
      nullify(dxv(H_ENERGY_DENS)%d)
      nullify(dyv(H_ENERGY_DENS)%d)
      nullify(dzv(H_ENERGY_DENS)%d)
      nullify(dxv(H_NOETHER_DENS)%d)
      nullify(dyv(H_NOETHER_DENS)%d)
      nullify(dzv(H_NOETHER_DENS)%d)
      nullify(dxv(H_ENERGY_CONST)%d)
      nullify(dyv(H_ENERGY_CONST)%d)
      nullify(dzv(H_ENERGY_CONST)%d)
      nullify(dxv(H_MOMX_CONST)%d)
      nullify(dyv(H_MOMX_CONST)%d)
      nullify(dzv(H_MOMX_CONST)%d)
      nullify(dxv(H_MOMY_CONST)%d)
      nullify(dyv(H_MOMY_CONST)%d)
      nullify(dzv(H_MOMY_CONST)%d)
      nullify(dxv(H_MOMZ_CONST)%d)
      nullify(dyv(H_MOMZ_CONST)%d)
      nullify(dzv(H_MOMZ_CONST)%d)
      nullify(dxv(H_PSI4R)%d)
      nullify(dyv(H_PSI4R)%d)
      nullify(dzv(H_PSI4R)%d)
      nullify(dxv(H_PSI4I)%d)
      nullify(dyv(H_PSI4I)%d)
      nullify(dzv(H_PSI4I)%d)
      nullify(dxv(H_MASSADM)%d)
      nullify(dyv(H_MASSADM)%d)
      nullify(dzv(H_MASSADM)%d)
      nullify(dxv(H_CURVATURE)%d)
      nullify(dyv(H_CURVATURE)%d)
      nullify(dzv(H_CURVATURE)%d)
      nullify(dxv(H_LIN_MOMX)%d)
      nullify(dyv(H_LIN_MOMX)%d)
      nullify(dzv(H_LIN_MOMX)%d)
      nullify(dxv(H_LIN_MOMY)%d)
      nullify(dyv(H_LIN_MOMY)%d)
      nullify(dzv(H_LIN_MOMY)%d)
      nullify(dxv(H_ANG_MOMZ)%d)
      nullify(dyv(H_ANG_MOMZ)%d)
      nullify(dzv(H_ANG_MOMZ)%d)
      nullify(dxv(H_KANG_MOMZ)%d)
      nullify(dyv(H_KANG_MOMZ)%d)
      nullify(dzv(H_KANG_MOMZ)%d)
      nullify(dxv(H_R2SCALAR)%d)
      nullify(dyv(H_R2SCALAR)%d)
      nullify(dzv(H_R2SCALAR)%d)
      nullify(dxv(H_GTHETHE)%d)
      nullify(dyv(H_GTHETHE)%d)
      nullify(dzv(H_GTHETHE)%d)
      nullify(dxv(H_GTHEPHI)%d)
      nullify(dyv(H_GTHEPHI)%d)
      nullify(dzv(H_GTHEPHI)%d)
      nullify(dxv(H_GPHIPHI)%d)
      nullify(dyv(H_GPHIPHI)%d)
      nullify(dzv(H_GPHIPHI)%d)
      nullify(dxv(H_GUR)%d)
      nullify(dyv(H_GUR)%d)
      nullify(dzv(H_GUR)%d)
      nullify(dxv(H_PSI_EL)%d)
      nullify(dyv(H_PSI_EL)%d)
      nullify(dzv(H_PSI_EL)%d)
      nullify(dxv(H_BETAX_EL)%d)
      nullify(dyv(H_BETAX_EL)%d)
      nullify(dzv(H_BETAX_EL)%d)
      nullify(dxv(H_BETAY_EL)%d)
      nullify(dyv(H_BETAY_EL)%d)
      nullify(dzv(H_BETAY_EL)%d)
      nullify(dxv(H_BETAZ_EL)%d)
      nullify(dyv(H_BETAZ_EL)%d)
      nullify(dzv(H_BETAZ_EL)%d)
      return
      end subroutine assign_ptrs_derivs
      subroutine assign_params(gauge_type,initial_H,idtype,boost_type,rh
     &s_type,evolve_geometry,evolve_scalarfield,CPBC_type,damping_matter
     &,sigma0,sigma1,sigma2,sigma0_rho,chi1,chi2,chin,chit0,chit1,dampin
     &g_rmax,damping_factor,sigma0_slope,alphaadjust_decay,read_file_alp
     &,read_file_b,read_file_g,read_file_k,read_file_psi,read_file_phir,
     &read_file_phic,read_file_phim,read_file_rho,read_data_level,id_d_d
     &iff_order,mr_amp,mr_shift,tw_amp,tw_width,tw_multipole,tw_time_shi
     &ft,bh1_matching,bh2_matching,sf_ampm,sf_omega1,sf_omega2,sf_omega3
     &,sf_mass,sf_lambda,sf_phase,bs1_x0,bs1_y0,bs1_z0,bs2_x0,bs2_y0,bs2
     &_z0,bs3_x0,bs3_y0,bs3_z0,bs1_vx,bs1_vy,bs1_vz,bs2_vx,bs2_vy,bs2_vz
     &,bs3_vx,bs3_vy,bs3_vz,G_scale_factor,number_surfaces,outer_boundar
     &y,masked_points,dynamical_method,amp_phi,nx0,ny0,nz0,nt0,h,hx,hy,h
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
       integer initial_H
       integer idtype
       integer boost_type
       integer rhs_type
       integer evolve_geometry
       integer evolve_scalarfield
       integer CPBC_type
       integer damping_matter
       real(kind=8) sigma0
       real(kind=8) sigma1
       real(kind=8) sigma2
       real(kind=8) sigma0_rho
       real(kind=8) chi1
       real(kind=8) chi2
       real(kind=8) chin
       real(kind=8) chit0
       real(kind=8) chit1
       real(kind=8) damping_rmax
       real(kind=8) damping_factor
       real(kind=8) sigma0_slope
       real(kind=8) alphaadjust_decay
       integer read_file_alp
       integer read_file_b
       integer read_file_g
       integer read_file_k
       integer read_file_psi
       integer read_file_phir
       integer read_file_phic
       integer read_file_phim
       integer read_file_rho
       integer read_data_level
       integer id_d_diff_order
       real(kind=8) mr_amp
       real(kind=8) mr_shift
       real(kind=8) tw_amp
       real(kind=8) tw_width
       integer tw_multipole
       real(kind=8) tw_time_shift
       integer bh1_matching
       integer bh2_matching
       real(kind=8) sf_ampm
       real(kind=8) sf_omega1
       real(kind=8) sf_omega2
       real(kind=8) sf_omega3
       real(kind=8) sf_mass
       real(kind=8) sf_lambda
       real(kind=8) sf_phase
       real(kind=8) bs1_x0
       real(kind=8) bs1_y0
       real(kind=8) bs1_z0
       real(kind=8) bs2_x0
       real(kind=8) bs2_y0
       real(kind=8) bs2_z0
       real(kind=8) bs3_x0
       real(kind=8) bs3_y0
       real(kind=8) bs3_z0
       real(kind=8) bs1_vx
       real(kind=8) bs1_vy
       real(kind=8) bs1_vz
       real(kind=8) bs2_vx
       real(kind=8) bs2_vy
       real(kind=8) bs2_vz
       real(kind=8) bs3_vx
       real(kind=8) bs3_vy
       real(kind=8) bs3_vz
       real(kind=8) G_scale_factor
       integer number_surfaces
       real(kind=8) outer_boundary
       real(kind=8) masked_points
       integer dynamical_method
       real(kind=8) amp_phi
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
      par(P_INITIAL_H) = initial_H
      par(P_IDTYPE) = idtype
      par(P_BOOST_TYPE) = boost_type
      par(P_RHS_TYPE) = rhs_type
      par(P_EVOLVE_GEOMETRY) = evolve_geometry
      par(P_EVOLVE_SCALARFIELD) = evolve_scalarfield
      par(P_CPBC_TYPE) = CPBC_type
      par(P_DAMPING_MATTER) = damping_matter
      par(P_SIGMA0) = sigma0
      par(P_SIGMA1) = sigma1
      par(P_SIGMA2) = sigma2
      par(P_SIGMA0_RHO) = sigma0_rho
      par(P_CHI1) = chi1
      par(P_CHI2) = chi2
      par(P_CHIN) = chin
      par(P_CHIT0) = chit0
      par(P_CHIT1) = chit1
      par(P_DAMPING_RMAX) = damping_rmax
      par(P_DAMPING_FACTOR) = damping_factor
      par(P_SIGMA0_SLOPE) = sigma0_slope
      par(P_ALPHAADJUST_DECAY) = alphaadjust_decay
      par(P_READ_FILE_ALP) = read_file_alp
      par(P_READ_FILE_B) = read_file_b
      par(P_READ_FILE_G) = read_file_g
      par(P_READ_FILE_K) = read_file_k
      par(P_READ_FILE_PSI) = read_file_psi
      par(P_READ_FILE_PHIR) = read_file_phir
      par(P_READ_FILE_PHIC) = read_file_phic
      par(P_READ_FILE_PHIM) = read_file_phim
      par(P_READ_FILE_RHO) = read_file_rho
      par(P_READ_DATA_LEVEL) = read_data_level
      par(P_ID_D_DIFF_ORDER) = id_d_diff_order
      par(P_MR_AMP) = mr_amp
      par(P_MR_SHIFT) = mr_shift
      par(P_TW_AMP) = tw_amp
      par(P_TW_WIDTH) = tw_width
      par(P_TW_MULTIPOLE) = tw_multipole
      par(P_TW_TIME_SHIFT) = tw_time_shift
      par(P_BH1_MATCHING) = bh1_matching
      par(P_BH2_MATCHING) = bh2_matching
      par(P_SF_AMPM) = sf_ampm
      par(P_SF_OMEGA1) = sf_omega1
      par(P_SF_OMEGA2) = sf_omega2
      par(P_SF_OMEGA3) = sf_omega3
      par(P_SF_MASS) = sf_mass
      par(P_SF_LAMBDA) = sf_lambda
      par(P_SF_PHASE) = sf_phase
      par(P_BS1_X0) = bs1_x0
      par(P_BS1_Y0) = bs1_y0
      par(P_BS1_Z0) = bs1_z0
      par(P_BS2_X0) = bs2_x0
      par(P_BS2_Y0) = bs2_y0
      par(P_BS2_Z0) = bs2_z0
      par(P_BS3_X0) = bs3_x0
      par(P_BS3_Y0) = bs3_y0
      par(P_BS3_Z0) = bs3_z0
      par(P_BS1_VX) = bs1_vx
      par(P_BS1_VY) = bs1_vy
      par(P_BS1_VZ) = bs1_vz
      par(P_BS2_VX) = bs2_vx
      par(P_BS2_VY) = bs2_vy
      par(P_BS2_VZ) = bs2_vz
      par(P_BS3_VX) = bs3_vx
      par(P_BS3_VY) = bs3_vy
      par(P_BS3_VZ) = bs3_vz
      par(P_G_SCALE_FACTOR) = G_scale_factor
      par(P_NUMBER_SURFACES) = number_surfaces
      par(P_OUTER_BOUNDARY) = outer_boundary
      par(P_MASKED_POINTS) = masked_points
      par(P_DYNAMICAL_METHOD) = dynamical_method
      par(P_AMP_PHI) = amp_phi
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
