#include "cctk.h"

      !----------------------------------------------------------------
      !
      !
      !----------------------------------------------------------------
      subroutine calc_bssn_bcs_rhs( dtu, u, v, w, imask,
#include "BC_sub_args.h"
     &                              x, y, z, nx, ny, nz, par) 
        use GF
        use params
        implicit   none

        type(gridfunction), dimension(NU_G)     :: dtu 
        type(gridfunction), dimension(NU_G)     :: u 
        type(gridfunction), dimension(NV_G)     :: v 
        type(gridfunction), dimension(NW)       :: w 
        
        CCTK_INT , dimension(nx,ny,nz)          :: imask  
        CCTK_REAL, dimension(NPAR), intent(in)  :: par 
        CCTK_INT                                :: nx, ny, nz
        CCTK_REAL, dimension(nx), intent(in)    :: x
        CCTK_REAL, dimension(ny), intent(in)    :: y
        CCTK_REAL, dimension(nz), intent(in)    :: z

#include "BC_decl_sub_args.h"

        CCTK_REAL  dx, dy, dz 

        ! local variables 
        CCTK_REAL, dimension(:,:,:), pointer :: Alpha3D 
        CCTK_REAL, dimension(:,:,:), pointer :: chi3D
        CCTK_REAL, dimension(:,:,:), pointer :: trK3D
        CCTK_REAL, dimension(:,:,:), pointer :: shiftx
        CCTK_REAL, dimension(:,:,:), pointer :: shifty
        CCTK_REAL, dimension(:,:,:), pointer :: shiftz
        CCTK_REAL, dimension(:,:,:), pointer :: Gamtx
        CCTK_REAL, dimension(:,:,:), pointer :: Gamty
        CCTK_REAL, dimension(:,:,:), pointer :: Gamtz
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
        CCTK_REAL, dimension(:,:,:), pointer :: gbx
        CCTK_REAL, dimension(:,:,:), pointer :: gby
        CCTK_REAL, dimension(:,:,:), pointer :: gbz

        CCTK_REAL, dimension(:,:,:), pointer :: Ex
        CCTK_REAL, dimension(:,:,:), pointer :: Ey
        CCTK_REAL, dimension(:,:,:), pointer :: Ez
        CCTK_REAL, dimension(:,:,:), pointer :: Bx
        CCTK_REAL, dimension(:,:,:), pointer :: By
        CCTK_REAL, dimension(:,:,:), pointer :: Bz
        CCTK_REAL, dimension(:,:,:), pointer :: Phi_em
        CCTK_REAL, dimension(:,:,:), pointer :: Psi_em

        CCTK_REAL, dimension(:,:,:), pointer :: PhiR
        CCTK_REAL, dimension(:,:,:), pointer :: PhiI
        CCTK_REAL, dimension(:,:,:), pointer :: PiR
        CCTK_REAL, dimension(:,:,:), pointer :: PiI

        CCTK_REAL, dimension(:,:,:), pointer :: dt_Alpha3D 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_chi3D 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_trK3D 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_shiftx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_shifty 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_shiftz 
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Gamtx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Gamty
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Gamtz
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
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gbx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gby
        CCTK_REAL, dimension(:,:,:), pointer :: dt_gbz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Ex
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Ey
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Ez
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Bx
        CCTK_REAL, dimension(:,:,:), pointer :: dt_By
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Bz
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Phi_em
        CCTK_REAL, dimension(:,:,:), pointer :: dt_Psi_em
        CCTK_REAL, dimension(:,:,:), pointer :: dt_PhiR
        CCTK_REAL, dimension(:,:,:), pointer :: dt_PhiI
        CCTK_REAL, dimension(:,:,:), pointer :: dt_PiR
        CCTK_REAL, dimension(:,:,:), pointer :: dt_PiI

        CCTK_REAL, dimension(:,:,:), pointer :: chr

        CCTK_INT    ::  i, j, k 
        CCTK_INT    ::  evolve_geometry, evolve_em_field
        CCTK_INT    ::  evolve_scalar_field
        CCTK_INT    ::  ichr

        CCTK_REAL   :: inv_r , r
        logical  , parameter     ::   ltrace = .false.

        CCTK_REAL, parameter     ::   Alpha3D_falloff = 1.0d0
        CCTK_REAL, parameter     ::     chi3D_falloff = 1.0d0
        CCTK_REAL, parameter     ::     trK3D_falloff = 1.0d0
      !!  CCTK_REAL, parameter     ::     trK3D_falloff = 2.0d0
        CCTK_REAL, parameter     ::    shiftx_falloff = 1.0d0
        CCTK_REAL, parameter     ::    shifty_falloff = 1.0d0
        CCTK_REAL, parameter     ::    shiftz_falloff = 1.0d0
        CCTK_REAL, parameter     ::      gtxx_falloff = 1.0d0
        CCTK_REAL, parameter     ::      gtxy_falloff = 1.0d0
        CCTK_REAL, parameter     ::      gtxz_falloff = 1.0d0
        CCTK_REAL, parameter     ::      gtyy_falloff = 1.0d0
        CCTK_REAL, parameter     ::      gtyz_falloff = 1.0d0
        CCTK_REAL, parameter     ::      gtzz_falloff = 1.0d0
      !!  CCTK_REAL, parameter     ::      Atxx_falloff = 1.0d0
      !!  CCTK_REAL, parameter     ::      Atxy_falloff = 1.0d0
      !!  CCTK_REAL, parameter     ::      Atxz_falloff = 1.0d0
      !!  CCTK_REAL, parameter     ::      Atyy_falloff = 1.0d0
      !!  CCTK_REAL, parameter     ::      Atyz_falloff = 1.0d0
      !!  CCTK_REAL, parameter     ::      Atzz_falloff = 1.0d0
        CCTK_REAL, parameter     ::      Atxx_falloff = 2.0d0
        CCTK_REAL, parameter     ::      Atxy_falloff = 2.0d0
        CCTK_REAL, parameter     ::      Atxz_falloff = 2.0d0
        CCTK_REAL, parameter     ::      Atyy_falloff = 2.0d0
        CCTK_REAL, parameter     ::      Atyz_falloff = 2.0d0
        CCTK_REAL, parameter     ::      Atzz_falloff = 2.0d0
      !!  CCTK_REAL, parameter     ::     Gamtx_falloff = 1.0d0
      !!  CCTK_REAL, parameter     ::     Gamty_falloff = 1.0d0
      !!  CCTK_REAL, parameter     ::     Gamtz_falloff = 1.0d0
        CCTK_REAL, parameter     ::     Gamtx_falloff = 2.0d0
        CCTK_REAL, parameter     ::     Gamty_falloff = 2.0d0
        CCTK_REAL, parameter     ::     Gamtz_falloff = 2.0d0
        CCTK_REAL, parameter     ::       gbx_falloff = 1.0d0
        CCTK_REAL, parameter     ::       gby_falloff = 1.0d0
        CCTK_REAL, parameter     ::       gbz_falloff = 1.0d0

        CCTK_REAL, parameter     ::        Ex_falloff = 2.0d0
        CCTK_REAL, parameter     ::        Ey_falloff = 2.0d0
        CCTK_REAL, parameter     ::        Ez_falloff = 2.0d0
        CCTK_REAL, parameter     ::        Bx_falloff = 2.0d0
        CCTK_REAL, parameter     ::        By_falloff = 2.0d0
        CCTK_REAL, parameter     ::        Bz_falloff = 2.0d0
        CCTK_REAL, parameter     ::    Phi_em_falloff = 3.0d0
        CCTK_REAL, parameter     ::    Psi_em_falloff = 3.0d0

        CCTK_REAL, parameter     ::    PhiR_falloff = 1.0d0
        CCTK_REAL, parameter     ::    PhiI_falloff = 1.0d0
        CCTK_REAL, parameter     ::     PiR_falloff = 2.0d0
        CCTK_REAL, parameter     ::     PiI_falloff = 2.0d0


        CCTK_REAL, parameter     ::   Alpha3D_asymptotic = 1.0d0
        CCTK_REAL, parameter     ::     chi3D_asymptotic = 1.0d0
        CCTK_REAL, parameter     ::     trK3D_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::    shiftx_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::    shifty_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::    shiftz_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::      gtxx_asymptotic = 1.0d0
        CCTK_REAL, parameter     ::      gtxy_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::      gtxz_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::      gtyy_asymptotic = 1.0d0
        CCTK_REAL, parameter     ::      gtyz_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::      gtzz_asymptotic = 1.0d0
        CCTK_REAL, parameter     ::      Atxx_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::      Atxy_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::      Atxz_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::      Atyy_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::      Atyz_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::      Atzz_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::     Gamtx_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::     Gamty_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::     Gamtz_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::       gbx_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::       gby_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::       gbz_asymptotic = 0.0d0

        CCTK_REAL, parameter     ::        Ex_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::        Ey_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::        Ez_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::        Bx_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::        By_asymptotic = 0.0d0
        CCTK_REAL                ::        Bz_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::    Phi_em_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::    Psi_em_asymptotic = 0.0d0

        CCTK_REAL                ::    PhiR_asymptotic 
        CCTK_REAL, parameter     ::    PhiI_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::     PiR_asymptotic = 0.0d0
        CCTK_REAL, parameter     ::     PiI_asymptotic = 0.0d0


        CCTK_INT                 ::  n_bhs
        CCTK_REAL                ::  mass1, mass2
        CCTK_REAL                ::  bh1x, bh1y, bh1z
        CCTK_REAL                ::  bh2x, bh2y, bh2z, vpsibl
        CCTK_REAL                ::  x1, y1, z1, x2, y2, z2, rv1, rv2


        dx = par(P_DX) 
        dy = par(P_DY) 
        dz = par(P_DZ) 

        evolve_geometry     = nint(par(P_EVOLVE_GEOMETRY))
        evolve_em_field     = nint(par(P_EVOLVE_EM_FIELD))
        evolve_scalar_field = nint(par(P_EVOLVE_SCALAR_FIELD))

        PhiR_asymptotic = par(P_DIL_INFTY)

        Alpha3D     =>   u(G_ALPHA)%d
        shiftx      =>   u(G_SHIFT1)%d
        shifty      =>   u(G_SHIFT2)%d
        shiftz      =>   u(G_SHIFT3)%d
        chi3D       =>   u(G_CHI)%d
        trK3D       =>   u(G_TRK)%d
        gtxx        =>   u(G_GT11)%d
        gtxy        =>   u(G_GT12)%d
        gtxz        =>   u(G_GT13)%d
        gtyy        =>   u(G_GT22)%d
        gtyz        =>   u(G_GT23)%d
        gtzz        =>   u(G_GT33)%d
        Atxx        =>   u(G_A11)%d
        Atxy        =>   u(G_A12)%d
        Atxz        =>   u(G_A13)%d
        Atyy        =>   u(G_A22)%d
        Atyz        =>   u(G_A23)%d
        Atzz        =>   u(G_A33)%d
        Gamtx       =>   u(G_GAM1)%d
        Gamty       =>   u(G_GAM2)%d
        Gamtz       =>   u(G_GAM3)%d
        gbx         =>   u(G_GB1)%d
        gby         =>   u(G_GB2)%d
        gbz         =>   u(G_GB3)%d

        Ex          =>   u(G_EX)%d
        Ey          =>   u(G_EY)%d
        Ez          =>   u(G_EZ)%d
        Bx          =>   u(G_BX)%d
        By          =>   u(G_BY)%d
        Bz          =>   u(G_BZ)%d
        Psi_em      =>   u(G_PSI_EM)%d
        Phi_em      =>   u(G_PHI_EM)%d

        PhiR        =>   u(G_PHIR)%d
        PhiI        =>   u(G_PHII)%d
        PiR         =>   u(G_PIR)%d
        PiI         =>   u(G_PII)%d

        dt_Alpha3D  =>   dtu(G_ALPHA)%d
        dt_shiftx   =>   dtu(G_SHIFT1)%d
        dt_shifty   =>   dtu(G_SHIFT2)%d
        dt_shiftz   =>   dtu(G_SHIFT3)%d
        dt_chi3D    =>   dtu(G_CHI)%d
        dt_trK3D    =>   dtu(G_TRK)%d
        dt_gtxx     =>   dtu(G_GT11)%d
        dt_gtxy     =>   dtu(G_GT12)%d
        dt_gtxz     =>   dtu(G_GT13)%d
        dt_gtyy     =>   dtu(G_GT22)%d
        dt_gtyz     =>   dtu(G_GT23)%d
        dt_gtzz     =>   dtu(G_GT33)%d
        dt_Atxx     =>   dtu(G_A11)%d
        dt_Atxy     =>   dtu(G_A12)%d
        dt_Atxz     =>   dtu(G_A13)%d
        dt_Atyy     =>   dtu(G_A22)%d
        dt_Atyz     =>   dtu(G_A23)%d
        dt_Atzz     =>   dtu(G_A33)%d
        dt_Gamtx    =>   dtu(G_GAM1)%d
        dt_Gamty    =>   dtu(G_GAM2)%d
        dt_Gamtz    =>   dtu(G_GAM3)%d
        dt_gbx      =>   dtu(G_GB1)%d
        dt_gby      =>   dtu(G_GB2)%d
        dt_gbz      =>   dtu(G_GB3)%d

        dt_Ex       =>   dtu(G_EX)%d
        dt_Ey       =>   dtu(G_EY)%d
        dt_Ez       =>   dtu(G_EZ)%d
        dt_Bx       =>   dtu(G_BX)%d
        dt_By       =>   dtu(G_BY)%d
        dt_Bz       =>   dtu(G_BZ)%d
        dt_Psi_em   =>   dtu(G_PSI_EM)%d
        dt_Phi_em   =>   dtu(G_PHI_EM)%d

        dt_PhiR     =>   dtu(G_PHIR)%d
        dt_PhiI     =>   dtu(G_PHII)%d
        dt_PiR      =>   dtu(G_PIR)%d
        dt_PiI      =>   dtu(G_PII)%d

        chr         =>   w(H_CHR)%d

!       black hole 1:
!        n_bhs = nint(par(P_BH_N))
!        mass1 = par(P_BH_1_MASS)
!        bh1x  = par(P_BH_1_X)
!        bh1y  = par(P_BH_1_Y)
!        bh1z  = par(P_BH_1_Z)

!       black hole 2:
!        mass2 = par(P_BH_2_MASS)
!        bh2x  = par(P_BH_2_X)
!        bh2y  = par(P_BH_2_Y)
!        bh2z  = par(P_BH_2_Z)


! We will label the faces in the following way: (Some of the equals signs below really denote "is in the interval")  
!  (F1)  x= xmin        y=(ymin,ymax)  z=(zmin,zmax)  
!  (F2)  x=      xmax   y=(ymin,ymax)  z=(zmin,zmax)  
!  (F3)  x=(xmin,xmax)  y= ymin        z=(zmin,zmax)  
!  (F4)  x=(xmin,xmax)  y=      ymax   z=(zmin,zmax)  
!  (F5)  x=(xmin,xmax)  y=(ymin,ymax)  z= zmin 
!  (F6)  x=(xmin,xmax)  y=(ymin,ymax)  z=      zmax 

!  We will also label the edges:
!  (E1)  x= xmin        y= ymin        z=(zmin,zmax)  
!  (E2)  x= xmin        y=      ymax   z=(zmin,zmax) 
!  (E3)  x=      xmax   y= ymin        z=(zmin,zmax) 
!  (E4)  x=      xmax   y=      ymax   z=(zmin,zmax) 
!  (E5)  x= xmin        y=(ymin,ymax)  z=zmin
!  (E6)  x= xmin        y=(ymin,ymax)  z=zmax  
!  (E7)  x=      xmax   y=(ymin,ymax)  z=zmin
!  (E8)  x=      xmax   y=(ymin,ymax)  z=zmax  
!  (E9)  x=(xmin,xmax)  y=ymin         z=zmin
!  (E10) x=(xmin,xmax)  y=ymin         z=zmax  
!  (E11) x=(xmin,xmax)  y=ymax         z=zmin
!  (E12) x=(xmin,xmax)  y=ymax         z=zmax  

! Labeling the corners 
!  (C1)  x= xmin        y= ymin        z= zmin 
!  (C2)  x= xmin        y= ymin        z=      zmax  
!  (C3)  x= xmin        y=      ymax   z= zmin 
!  (C4)  x= xmin        y=      ymax   z=      zmax  
!  (C5)  x=      xmax   y= ymin        z= zmin 
!  (C6)  x=      xmax   y= ymin        z=      zmax  
!  (C7)  x=      xmax   y=      ymax   z= zmin 
!  (C8)  x=      xmax   y=      ymax   z=      zmax  


! F1, F2, E1-E8, C1-C8  
        do k = 1, nz
        do j = 1, ny 
        do i = 1, nx, nx-1  
           ichr = nint(chr(i,j,k))
           if (ichr .ne. 8 .and. ichr .ne. 1) then

           r = sqrt(x(i)**2 + y(j)**2 + z(k)**2) 
           if (r .gt. 1.0d-1) then
             inv_r = 1.0d0 / r
           else 
             inv_r = 1.0d0 
             if (ltrace) 
     &         write(0,*)'calc_bssn_bcs_rhs: rsq < 0.1. chr=',chr(i,j,k)
           end if
                              
           if (evolve_geometry .eq. 1) then
           dt_Alpha3D(i,j,k) =   x(i) * dx_alpha(i,j,k) 
     &                         + y(j) * dy_alpha(i,j,k) 
     &                         + z(k) * dz_alpha(i,j,k)  
     &                         + Alpha3D_falloff 
     &                             * (   Alpha3D(i,j,k)  
     &                                 - Alpha3D_asymptotic ) 
           dt_Alpha3D(i,j,k) = - inv_r * dt_Alpha3D(i,j,k) 

           dt_shiftx(i,j,k) =    x(i) * dx_shiftx(i,j,k) 
     &                         + y(j) * dy_shiftx(i,j,k) 
     &                         + z(k) * dz_shiftx(i,j,k)  
     &                         + shiftx_falloff 
     &                             * (   shiftx(i,j,k)  
     &                                 - shiftx_asymptotic ) 
           dt_shiftx(i,j,k) = - inv_r * dt_shiftx(i,j,k) 

           dt_shifty(i,j,k) =    x(i) * dx_shifty(i,j,k) 
     &                         + y(j) * dy_shifty(i,j,k) 
     &                         + z(k) * dz_shifty(i,j,k)  
     &                         + shifty_falloff 
     &                             * (   shifty(i,j,k)  
     &                                 - shifty_asymptotic ) 
           dt_shifty(i,j,k) = - inv_r * dt_shifty(i,j,k) 

           dt_shiftz(i,j,k) =    x(i) * dx_shiftz(i,j,k) 
     &                         + y(j) * dy_shiftz(i,j,k) 
     &                         + z(k) * dz_shiftz(i,j,k)  
     &                         + shiftz_falloff 
     &                             * (   shiftz(i,j,k) 
     &                                 - shiftz_asymptotic ) 
           dt_shiftz(i,j,k) = - inv_r * dt_shiftz(i,j,k) 

           dt_chi3D(i,j,k) =    x(i) * dx_chi(i,j,k) 
     &                        + y(j) * dy_chi(i,j,k) 
     &                        + z(k) * dz_chi(i,j,k)  
     &                        + chi3D_falloff 
     &                            * (   chi3D(i,j,k) 
     &                                - chi3D_asymptotic ) 
           dt_chi3D(i,j,k) = - inv_r * dt_chi3D(i,j,k) 

           dt_trK3D(i,j,k) =    x(i) * dx_trK(i,j,k) 
     &                        + y(j) * dy_trK(i,j,k) 
     &                        + z(k) * dz_trK(i,j,k)  
     &                        + trK3D_falloff 
     &                            * (   trK3D(i,j,k) 
     &                                - trK3D_asymptotic ) 
           dt_trK3D(i,j,k) = - inv_r * dt_trK3D(i,j,k) 

           dt_gtxx(i,j,k) =    x(i) * dx_gtxx(i,j,k) 
     &                       + y(j) * dy_gtxx(i,j,k) 
     &                       + z(k) * dz_gtxx(i,j,k)  
     &                       + gtxx_falloff 
     &                           * (   gtxx(i,j,k) 
     &                               - gtxx_asymptotic ) 
           dt_gtxx(i,j,k) = - inv_r * dt_gtxx(i,j,k) 

           dt_gtxy(i,j,k) =    x(i) * dx_gtxy(i,j,k) 
     &                       + y(j) * dy_gtxy(i,j,k) 
     &                       + z(k) * dz_gtxy(i,j,k)  
     &                       + gtxy_falloff 
     &                           * (   gtxy(i,j,k) 
     &                               - gtxy_asymptotic ) 
           dt_gtxy(i,j,k) = - inv_r * dt_gtxy(i,j,k) 

           dt_gtxz(i,j,k) =    x(i) * dx_gtxz(i,j,k) 
     &                       + y(j) * dy_gtxz(i,j,k) 
     &                       + z(k) * dz_gtxz(i,j,k)  
     &                       + gtxz_falloff 
     &                           * (   gtxz(i,j,k) 
     &                               - gtxz_asymptotic ) 
           dt_gtxz(i,j,k) = - inv_r * dt_gtxz(i,j,k) 

           dt_gtyy(i,j,k) =    x(i) * dx_gtyy(i,j,k) 
     &                       + y(j) * dy_gtyy(i,j,k) 
     &                       + z(k) * dz_gtyy(i,j,k)  
     &                       + gtyy_falloff 
     &                           * (   gtyy(i,j,k) 
     &                               - gtyy_asymptotic ) 
           dt_gtyy(i,j,k) = - inv_r * dt_gtyy(i,j,k) 

           dt_gtyz(i,j,k) =    x(i) * dx_gtyz(i,j,k) 
     &                       + y(j) * dy_gtyz(i,j,k) 
     &                       + z(k) * dz_gtyz(i,j,k)  
     &                       + gtyz_falloff 
     &                           * (   gtyz(i,j,k) 
     &                               - gtyz_asymptotic ) 
           dt_gtyz(i,j,k) = - inv_r * dt_gtyz(i,j,k) 

           dt_gtzz(i,j,k) =    x(i) * dx_gtzz(i,j,k) 
     &                       + y(j) * dy_gtzz(i,j,k) 
     &                       + z(k) * dz_gtzz(i,j,k)  
     &                       + gtzz_falloff 
     &                           * (   gtzz(i,j,k) 
     &                               - gtzz_asymptotic ) 
           dt_gtzz(i,j,k) = - inv_r * dt_gtzz(i,j,k) 

           dt_Atxx(i,j,k) =    x(i) * dx_Atxx(i,j,k) 
     &                       + y(j) * dy_Atxx(i,j,k) 
     &                       + z(k) * dz_Atxx(i,j,k)  
     &                       + Atxx_falloff 
     &                           * (   Atxx(i,j,k) 
     &                               - Atxx_asymptotic ) 
           dt_Atxx(i,j,k) = - inv_r * dt_Atxx(i,j,k) 

           dt_Atxy(i,j,k) =    x(i) * dx_Atxy(i,j,k) 
     &                       + y(j) * dy_Atxy(i,j,k) 
     &                       + z(k) * dz_Atxy(i,j,k)  
     &                       + Atxy_falloff 
     &                           * (   Atxy(i,j,k) 
     &                               - Atxy_asymptotic ) 
           dt_Atxy(i,j,k) = - inv_r * dt_Atxy(i,j,k) 

           dt_Atxz(i,j,k) =    x(i) * dx_Atxz(i,j,k) 
     &                       + y(j) * dy_Atxz(i,j,k) 
     &                       + z(k) * dz_Atxz(i,j,k)  
     &                       + Atxz_falloff 
     &                           * (   Atxz(i,j,k) 
     &                               - Atxz_asymptotic ) 
           dt_Atxz(i,j,k) = - inv_r * dt_Atxz(i,j,k) 

           dt_Atyy(i,j,k) =    x(i) * dx_Atyy(i,j,k) 
     &                       + y(j) * dy_Atyy(i,j,k) 
     &                       + z(k) * dz_Atyy(i,j,k)  
     &                       + Atyy_falloff 
     &                           * (   Atyy(i,j,k) 
     &                               - Atyy_asymptotic ) 
           dt_Atyy(i,j,k) = - inv_r * dt_Atyy(i,j,k) 

           dt_Atyz(i,j,k) =    x(i) * dx_Atyz(i,j,k) 
     &                       + y(j) * dy_Atyz(i,j,k) 
     &                       + z(k) * dz_Atyz(i,j,k)  
     &                       + Atyz_falloff 
     &                           * (   Atyz(i,j,k) 
     &                               - Atyz_asymptotic ) 
           dt_Atyz(i,j,k) = - inv_r * dt_Atyz(i,j,k) 

           dt_Atzz(i,j,k) =    x(i) * dx_Atzz(i,j,k) 
     &                       + y(j) * dy_Atzz(i,j,k) 
     &                       + z(k) * dz_Atzz(i,j,k)  
     &                       + Atzz_falloff 
     &                           * (   Atzz(i,j,k) 
     &                               - Atzz_asymptotic ) 
           dt_Atzz(i,j,k) = - inv_r * dt_Atzz(i,j,k) 

           dt_Gamtx(i,j,k) =    x(i) * dx_Gamtx(i,j,k) 
     &                        + y(j) * dy_Gamtx(i,j,k) 
     &                        + z(k) * dz_Gamtx(i,j,k)  
     &                        + Gamtx_falloff 
     &                            * (   Gamtx(i,j,k) 
     &                                - Gamtx_asymptotic ) 
           dt_Gamtx(i,j,k) = - inv_r * dt_Gamtx(i,j,k) 

           dt_Gamty(i,j,k) =    x(i) * dx_Gamty(i,j,k) 
     &                        + y(j) * dy_Gamty(i,j,k) 
     &                        + z(k) * dz_Gamty(i,j,k)  
     &                        + Gamty_falloff 
     &                            * (   Gamty(i,j,k) 
     &                                - Gamty_asymptotic ) 
           dt_Gamty(i,j,k) = - inv_r * dt_Gamty(i,j,k) 

           dt_Gamtz(i,j,k) =    x(i) * dx_Gamtz(i,j,k) 
     &                        + y(j) * dy_Gamtz(i,j,k) 
     &                        + z(k) * dz_Gamtz(i,j,k)  
     &                        + Gamtz_falloff 
     &                            * (   Gamtz(i,j,k) 
     &                                - Gamtz_asymptotic ) 
           dt_Gamtz(i,j,k) = - inv_r * dt_Gamtz(i,j,k) 

           dt_gbx(i,j,k) =      x(i) * dx_gbx(i,j,k) 
     &                        + y(j) * dy_gbx(i,j,k) 
     &                        + z(k) * dz_gbx(i,j,k)  
     &                        + gbx_falloff 
     &                            * (   gbx(i,j,k) 
     &                                - gbx_asymptotic ) 
           dt_gbx(i,j,k) = - inv_r * dt_gbx(i,j,k) 

           dt_gby(i,j,k) =      x(i) * dx_gby(i,j,k) 
     &                        + y(j) * dy_gby(i,j,k) 
     &                        + z(k) * dz_gby(i,j,k)  
     &                        + gby_falloff 
     &                            * (   gby(i,j,k) 
     &                                - gby_asymptotic ) 
           dt_gby(i,j,k) = - inv_r * dt_gby(i,j,k) 

           dt_gbz(i,j,k) =      x(i) * dx_gbz(i,j,k) 
     &                        + y(j) * dy_gbz(i,j,k) 
     &                        + z(k) * dz_gbz(i,j,k)  
     &                        + gbz_falloff 
     &                            * (   gbz(i,j,k) 
     &                                - gbz_asymptotic ) 
           dt_gbz(i,j,k) = - inv_r * dt_gbz(i,j,k) 

           end if
           if (evolve_em_field .eq. 1) then
           dt_Ex(i,j,k) =      x(i) * dx_Ex(i,j,k) 
     &                        + y(j) * dy_Ex(i,j,k) 
     &                        + z(k) * dz_Ex(i,j,k)  
     &                        + Ex_falloff 
     &                            * (   Ex(i,j,k) 
     &                                - Ex_asymptotic ) 
           dt_Ex(i,j,k) = - inv_r * dt_Ex(i,j,k) 
           !dt_Ex(i,j,k) = 0.d0 * dt_Ex(i,j,k) 

           dt_Ey(i,j,k) =      x(i) * dx_Ey(i,j,k) 
     &                        + y(j) * dy_Ey(i,j,k) 
     &                        + z(k) * dz_Ey(i,j,k)  
     &                        + Ey_falloff 
     &                            * (   Ey(i,j,k) 
     &                                - Ey_asymptotic ) 
           dt_Ey(i,j,k) = - inv_r * dt_Ey(i,j,k) 
           !dt_Ey(i,j,k) = 0.d0 

           dt_Ez(i,j,k) =      x(i) * dx_Ez(i,j,k) 
     &                        + y(j) * dy_Ez(i,j,k) 
     &                        + z(k) * dz_Ez(i,j,k)  
     &                        + Ez_falloff 
     &                            * (   Ez(i,j,k) 
     &                                - Ez_asymptotic ) 
            dt_Ez(i,j,k) = - inv_r * dt_Ez(i,j,k) 
           !dt_Ez(i,j,k) = 0.d0 

           dt_Bx(i,j,k) =      x(i) * dx_Bx(i,j,k) 
     &                        + y(j) * dy_Bx(i,j,k) 
     &                        + z(k) * dz_Bx(i,j,k)  
     &                        + Bx_falloff 
     &                            * (   Bx(i,j,k) 
     &                                - Bx_asymptotic ) 
            dt_Bx(i,j,k) = - inv_r * dt_Bx(i,j,k) 
           !dt_Bx(i,j,k) = 0.d0 

           dt_By(i,j,k) =      x(i) * dx_By(i,j,k) 
     &                        + y(j) * dy_By(i,j,k) 
     &                        + z(k) * dz_By(i,j,k)  
     &                        + By_falloff 
     &                            * (   By(i,j,k) 
     &                                - By_asymptotic ) 
            dt_By(i,j,k) = - inv_r * dt_By(i,j,k) 
           !dt_By(i,j,k) = 0.d0 

           dt_Bz(i,j,k) =      x(i) * dx_Bz(i,j,k) 
     &                        + y(j) * dy_Bz(i,j,k) 
     &                        + z(k) * dz_Bz(i,j,k)  
     &                        + Bz_falloff 
     &                            * (   Bz(i,j,k) 
     &                                - Bz_asymptotic ) 
            dt_Bz(i,j,k) = - inv_r * dt_Bz(i,j,k) 
           !dt_Bz(i,j,k) = 0.d0 

           dt_Phi_em(i,j,k) =   x(i) * dx_Phi_em(i,j,k) 
     &                        + y(j) * dy_Phi_em(i,j,k) 
     &                        + z(k) * dz_Phi_em(i,j,k)  
     &                        + Phi_em_falloff 
     &                            * (   Phi_em(i,j,k) 
     &                                - Phi_em_asymptotic ) 
            dt_Phi_em(i,j,k) = - inv_r * dt_Phi_em(i,j,k) 
           !dt_Phi_em(i,j,k) = 0.d0 

           dt_Psi_em(i,j,k) =   x(i) * dx_Psi_em(i,j,k) 
     &                        + y(j) * dy_Psi_em(i,j,k) 
     &                        + z(k) * dz_Psi_em(i,j,k)  
     &                        + Psi_em_falloff 
     &                            * (   Psi_em(i,j,k) 
     &                                - Psi_em_asymptotic ) 
            dt_Psi_em(i,j,k) = - inv_r * dt_Psi_em(i,j,k) 
           !dt_Psi_em(i,j,k) = 0.d0 
            end if
            if (evolve_scalar_field .eq. 1) then
            dt_PhiR(i,j,k) =   x(i) * dx_PhiR(i,j,k) 
     &                       + y(j) * dy_PhiR(i,j,k) 
     &                       + z(k) * dz_PhiR(i,j,k)  
     &                       + PhiR_falloff 
     &                            * (   PhiR(i,j,k) 
     &                                - PhiR_asymptotic ) 
            dt_PhiR(i,j,k) = - inv_r * dt_PhiR(i,j,k) 

            dt_PhiI(i,j,k) =   x(i) * dx_PhiI(i,j,k) 
     &                      + y(j) * dy_PhiI(i,j,k) 
     &                      + z(k) * dz_PhiI(i,j,k)  
     &                      + PhiI_falloff 
     &                            * (   PhiI(i,j,k) 
     &                                - PhiI_asymptotic ) 
            dt_PhiI(i,j,k) = - inv_r * dt_PhiI(i,j,k) 

            dt_PiR(i,j,k) =   x(i) * dx_PiR(i,j,k) 
     &                     + y(j) * dy_PiR(i,j,k) 
     &                     + z(k) * dz_PiR(i,j,k)  
     &                     + PiR_falloff 
     &                            * (   PiR(i,j,k) 
     &                                - PiR_asymptotic ) 
            dt_PiR(i,j,k) = - inv_r * dt_PiR(i,j,k) 

            dt_PiI(i,j,k) =   x(i) * dx_PiI(i,j,k) 
     &                     + y(j) * dy_PiI(i,j,k) 
     &                     + z(k) * dz_PiI(i,j,k)  
     &                     + PiI_falloff 
     &                            * (   PiI(i,j,k) 
     &                                - PiI_asymptotic ) 
            dt_PiI(i,j,k) = - inv_r * dt_PiI(i,j,k) 
            end if
            end if

        end do 
        end do 
        end do 


! F3, F4, E9-E12   
!       do k = 1, nz
!       do j = 1, ny   , ny-1 
!       do i = 2, nx-1   
 
        do k = 1, nz
        do j = 1, ny   , ny-1 
        do i = 1, nx
           ichr = nint(chr(i,j,k))
           if (ichr .ne. 8 .and. ichr .ne. 1) then

           r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)
           if (r .gt. 1.0d-1) then
             inv_r = 1.0d0 / r
           else
             inv_r = 1.0d0
             if (ltrace) 
     &         write(0,*)'calc_bssn_bcs_rhs: rsq < 0.1. chr=',chr(i,j,k)
           end if

           if (evolve_geometry .eq. 1) then                   
           dt_Alpha3D(i,j,k) =   x(i) * dx_alpha(i,j,k) 
     &                         + y(j) * dy_alpha(i,j,k) 
     &                         + z(k) * dz_alpha(i,j,k)  
     &                         + Alpha3D_falloff 
     &                             * (   Alpha3D(i,j,k)  
     &                                 - Alpha3D_asymptotic ) 
           dt_Alpha3D(i,j,k) = - inv_r * dt_Alpha3D(i,j,k) 

           dt_shiftx(i,j,k) =    x(i) * dx_shiftx(i,j,k) 
     &                         + y(j) * dy_shiftx(i,j,k) 
     &                         + z(k) * dz_shiftx(i,j,k)  
     &                         + shiftx_falloff 
     &                             * (   shiftx(i,j,k)  
     &                                 - shiftx_asymptotic ) 
           dt_shiftx(i,j,k) = - inv_r * dt_shiftx(i,j,k) 

           dt_shifty(i,j,k) =    x(i) * dx_shifty(i,j,k) 
     &                         + y(j) * dy_shifty(i,j,k) 
     &                         + z(k) * dz_shifty(i,j,k)  
     &                         + shifty_falloff 
     &                             * (   shifty(i,j,k)  
     &                                 - shifty_asymptotic ) 
           dt_shifty(i,j,k) = - inv_r * dt_shifty(i,j,k) 

           dt_shiftz(i,j,k) =    x(i) * dx_shiftz(i,j,k) 
     &                         + y(j) * dy_shiftz(i,j,k) 
     &                         + z(k) * dz_shiftz(i,j,k)  
     &                         + shiftz_falloff 
     &                             * (   shiftz(i,j,k) 
     &                                 - shiftz_asymptotic ) 
           dt_shiftz(i,j,k) = - inv_r * dt_shiftz(i,j,k) 

           dt_chi3D(i,j,k) =    x(i) * dx_chi(i,j,k) 
     &                        + y(j) * dy_chi(i,j,k) 
     &                        + z(k) * dz_chi(i,j,k)  
     &                        + chi3D_falloff 
     &                            * (   chi3D(i,j,k) 
     &                                - chi3D_asymptotic ) 
           dt_chi3D(i,j,k) = - inv_r * dt_chi3D(i,j,k) 

           dt_trK3D(i,j,k) =    x(i) * dx_trK(i,j,k) 
     &                        + y(j) * dy_trK(i,j,k) 
     &                        + z(k) * dz_trK(i,j,k)  
     &                        + trK3D_falloff 
     &                            * (   trK3D(i,j,k) 
     &                                - trK3D_asymptotic ) 
           dt_trK3D(i,j,k) = - inv_r * dt_trK3D(i,j,k) 

           dt_gtxx(i,j,k) =    x(i) * dx_gtxx(i,j,k) 
     &                       + y(j) * dy_gtxx(i,j,k) 
     &                       + z(k) * dz_gtxx(i,j,k)  
     &                       + gtxx_falloff 
     &                           * (   gtxx(i,j,k) 
     &                               - gtxx_asymptotic ) 
           dt_gtxx(i,j,k) = - inv_r * dt_gtxx(i,j,k) 

           dt_gtxy(i,j,k) =    x(i) * dx_gtxy(i,j,k) 
     &                       + y(j) * dy_gtxy(i,j,k) 
     &                       + z(k) * dz_gtxy(i,j,k)  
     &                       + gtxy_falloff 
     &                           * (   gtxy(i,j,k) 
     &                               - gtxy_asymptotic ) 
           dt_gtxy(i,j,k) = - inv_r * dt_gtxy(i,j,k) 

           dt_gtxz(i,j,k) =    x(i) * dx_gtxz(i,j,k) 
     &                       + y(j) * dy_gtxz(i,j,k) 
     &                       + z(k) * dz_gtxz(i,j,k)  
     &                       + gtxz_falloff 
     &                           * (   gtxz(i,j,k) 
     &                               - gtxz_asymptotic ) 
           dt_gtxz(i,j,k) = - inv_r * dt_gtxz(i,j,k) 

           dt_gtyy(i,j,k) =    x(i) * dx_gtyy(i,j,k) 
     &                       + y(j) * dy_gtyy(i,j,k) 
     &                       + z(k) * dz_gtyy(i,j,k)  
     &                       + gtyy_falloff 
     &                           * (   gtyy(i,j,k) 
     &                               - gtyy_asymptotic ) 
           dt_gtyy(i,j,k) = - inv_r * dt_gtyy(i,j,k) 

           dt_gtyz(i,j,k) =    x(i) * dx_gtyz(i,j,k) 
     &                       + y(j) * dy_gtyz(i,j,k) 
     &                       + z(k) * dz_gtyz(i,j,k)  
     &                       + gtyz_falloff 
     &                           * (   gtyz(i,j,k) 
     &                               - gtyz_asymptotic ) 
           dt_gtyz(i,j,k) = - inv_r * dt_gtyz(i,j,k) 

           dt_gtzz(i,j,k) =    x(i) * dx_gtzz(i,j,k) 
     &                       + y(j) * dy_gtzz(i,j,k) 
     &                       + z(k) * dz_gtzz(i,j,k)  
     &                       + gtzz_falloff 
     &                           * (   gtzz(i,j,k) 
     &                               - gtzz_asymptotic ) 
           dt_gtzz(i,j,k) = - inv_r * dt_gtzz(i,j,k) 

           dt_Atxx(i,j,k) =    x(i) * dx_Atxx(i,j,k) 
     &                       + y(j) * dy_Atxx(i,j,k) 
     &                       + z(k) * dz_Atxx(i,j,k)  
     &                       + Atxx_falloff 
     &                           * (   Atxx(i,j,k) 
     &                               - Atxx_asymptotic ) 
           dt_Atxx(i,j,k) = - inv_r * dt_Atxx(i,j,k) 

           dt_Atxy(i,j,k) =    x(i) * dx_Atxy(i,j,k) 
     &                       + y(j) * dy_Atxy(i,j,k) 
     &                       + z(k) * dz_Atxy(i,j,k)  
     &                       + Atxy_falloff 
     &                           * (   Atxy(i,j,k) 
     &                               - Atxy_asymptotic ) 
           dt_Atxy(i,j,k) = - inv_r * dt_Atxy(i,j,k) 

           dt_Atxz(i,j,k) =    x(i) * dx_Atxz(i,j,k) 
     &                       + y(j) * dy_Atxz(i,j,k) 
     &                       + z(k) * dz_Atxz(i,j,k)  
     &                       + Atxz_falloff 
     &                           * (   Atxz(i,j,k) 
     &                               - Atxz_asymptotic ) 
           dt_Atxz(i,j,k) = - inv_r * dt_Atxz(i,j,k) 

           dt_Atyy(i,j,k) =    x(i) * dx_Atyy(i,j,k) 
     &                       + y(j) * dy_Atyy(i,j,k) 
     &                       + z(k) * dz_Atyy(i,j,k)  
     &                       + Atyy_falloff 
     &                           * (   Atyy(i,j,k) 
     &                               - Atyy_asymptotic ) 
           dt_Atyy(i,j,k) = - inv_r * dt_Atyy(i,j,k) 

           dt_Atyz(i,j,k) =    x(i) * dx_Atyz(i,j,k) 
     &                       + y(j) * dy_Atyz(i,j,k) 
     &                       + z(k) * dz_Atyz(i,j,k)  
     &                       + Atyz_falloff 
     &                           * (   Atyz(i,j,k) 
     &                               - Atyz_asymptotic ) 
           dt_Atyz(i,j,k) = - inv_r * dt_Atyz(i,j,k) 

           dt_Atzz(i,j,k) =    x(i) * dx_Atzz(i,j,k) 
     &                       + y(j) * dy_Atzz(i,j,k) 
     &                       + z(k) * dz_Atzz(i,j,k)  
     &                       + Atzz_falloff 
     &                           * (   Atzz(i,j,k) 
     &                               - Atzz_asymptotic ) 
           dt_Atzz(i,j,k) = - inv_r * dt_Atzz(i,j,k) 

           dt_Gamtx(i,j,k) =    x(i) * dx_Gamtx(i,j,k) 
     &                        + y(j) * dy_Gamtx(i,j,k) 
     &                        + z(k) * dz_Gamtx(i,j,k)  
     &                        + Gamtx_falloff 
     &                            * (   Gamtx(i,j,k) 
     &                                - Gamtx_asymptotic ) 
           dt_Gamtx(i,j,k) = - inv_r * dt_Gamtx(i,j,k) 

           dt_Gamty(i,j,k) =    x(i) * dx_Gamty(i,j,k) 
     &                        + y(j) * dy_Gamty(i,j,k) 
     &                        + z(k) * dz_Gamty(i,j,k)  
     &                        + Gamty_falloff 
     &                            * (   Gamty(i,j,k) 
     &                                - Gamty_asymptotic ) 
           dt_Gamty(i,j,k) = - inv_r * dt_Gamty(i,j,k) 

           dt_Gamtz(i,j,k) =    x(i) * dx_Gamtz(i,j,k) 
     &                        + y(j) * dy_Gamtz(i,j,k) 
     &                        + z(k) * dz_Gamtz(i,j,k)  
     &                        + Gamtz_falloff 
     &                            * (   Gamtz(i,j,k) 
     &                                - Gamtz_asymptotic ) 
           dt_Gamtz(i,j,k) = - inv_r * dt_Gamtz(i,j,k) 

           dt_gbx(i,j,k) =      x(i) * dx_gbx(i,j,k) 
     &                        + y(j) * dy_gbx(i,j,k) 
     &                        + z(k) * dz_gbx(i,j,k)  
     &                        + gbx_falloff 
     &                            * (   gbx(i,j,k) 
     &                                - gbx_asymptotic ) 
           dt_gbx(i,j,k) = - inv_r * dt_gbx(i,j,k) 

           dt_gby(i,j,k) =      x(i) * dx_gby(i,j,k) 
     &                        + y(j) * dy_gby(i,j,k) 
     &                        + z(k) * dz_gby(i,j,k)  
     &                        + gby_falloff 
     &                            * (   gby(i,j,k) 
     &                                - gby_asymptotic ) 
           dt_gby(i,j,k) = - inv_r * dt_gby(i,j,k) 

           dt_gbz(i,j,k) =      x(i) * dx_gbz(i,j,k) 
     &                        + y(j) * dy_gbz(i,j,k) 
     &                        + z(k) * dz_gbz(i,j,k)  
     &                        + gbz_falloff 
     &                            * (   gbz(i,j,k) 
     &                                - gbz_asymptotic ) 
           dt_gbz(i,j,k) = - inv_r * dt_gbz(i,j,k) 
           end if
           if (evolve_em_field .eq. 1) then
           dt_Ex(i,j,k) =      x(i) * dx_Ex(i,j,k) 
     &                        + y(j) * dy_Ex(i,j,k) 
     &                        + z(k) * dz_Ex(i,j,k)  
     &                        + Ex_falloff 
     &                            * (   Ex(i,j,k) 
     &                                - Ex_asymptotic ) 
            dt_Ex(i,j,k) = - inv_r * dt_Ex(i,j,k) 
           !dt_Ex(i,j,k) = 0.d0 

           dt_Ey(i,j,k) =      x(i) * dx_Ey(i,j,k) 
     &                        + y(j) * dy_Ey(i,j,k) 
     &                        + z(k) * dz_Ey(i,j,k)  
     &                        + Ey_falloff 
     &                            * (   Ey(i,j,k) 
     &                                - Ey_asymptotic ) 
           dt_Ey(i,j,k) = - inv_r * dt_Ey(i,j,k) 
           !dt_Ey(i,j,k) = 0.d0 

           dt_Ez(i,j,k) =      x(i) * dx_Ez(i,j,k) 
     &                        + y(j) * dy_Ez(i,j,k) 
     &                        + z(k) * dz_Ez(i,j,k)  
     &                        + Ez_falloff 
     &                            * (   Ez(i,j,k) 
     &                                - Ez_asymptotic ) 
           dt_Ez(i,j,k) = - inv_r * dt_Ez(i,j,k) 
           !dt_Ez(i,j,k) = 0.d0 

           dt_Bx(i,j,k) =      x(i) * dx_Bx(i,j,k) 
     &                        + y(j) * dy_Bx(i,j,k) 
     &                        + z(k) * dz_Bx(i,j,k)  
     &                        + Bx_falloff 
     &                            * (   Bx(i,j,k) 
     &                                - Bx_asymptotic ) 
           dt_Bx(i,j,k) = - inv_r * dt_Bx(i,j,k) 
           !dt_Bx(i,j,k) = 0.d0 

           dt_By(i,j,k) =      x(i) * dx_By(i,j,k) 
     &                        + y(j) * dy_By(i,j,k) 
     &                        + z(k) * dz_By(i,j,k)  
     &                        + By_falloff 
     &                            * (   By(i,j,k) 
     &                                - By_asymptotic ) 
           dt_By(i,j,k) = - inv_r * dt_By(i,j,k) 
           !dt_By(i,j,k) = 0.d0

           dt_Bz(i,j,k) =      x(i) * dx_Bz(i,j,k) 
     &                        + y(j) * dy_Bz(i,j,k) 
     &                        + z(k) * dz_Bz(i,j,k)  
     &                        + Bz_falloff 
     &                            * (   Bz(i,j,k) 
     &                                - Bz_asymptotic ) 
           dt_Bz(i,j,k) = - inv_r * dt_Bz(i,j,k) 
           !dt_Bz(i,j,k) = 0.d0 

           dt_Phi_em(i,j,k) =      x(i) * dx_Phi_em(i,j,k) 
     &                        + y(j) * dy_Phi_em(i,j,k) 
     &                        + z(k) * dz_Phi_em(i,j,k)  
     &                        + Phi_em_falloff 
     &                            * (   Phi_em(i,j,k) 
     &                                - Phi_em_asymptotic ) 
           dt_Phi_em(i,j,k) = - inv_r * dt_Phi_em(i,j,k) 
           !dt_Phi_em(i,j,k) = 0.d0

           dt_Psi_em(i,j,k) =      x(i) * dx_Psi_em(i,j,k) 
     &                        + y(j) * dy_Psi_em(i,j,k) 
     &                        + z(k) * dz_Psi_em(i,j,k)  
     &                        + Psi_em_falloff 
     &                            * (   Psi_em(i,j,k) 
     &                                - Psi_em_asymptotic ) 
           dt_Psi_em(i,j,k) = - inv_r * dt_Psi_em(i,j,k) 
           !dt_Psi_em(i,j,k) = 0.d0 
           end if
           if (evolve_scalar_field .eq. 1) then
           dt_PhiR(i,j,k) =   x(i) * dx_PhiR(i,j,k) 
     &                      + y(j) * dy_PhiR(i,j,k) 
     &                      + z(k) * dz_PhiR(i,j,k)  
     &                      + PhiR_falloff 
     &                            * (   PhiR(i,j,k) 
     &                                - PhiR_asymptotic ) 
            dt_PhiR(i,j,k) = - inv_r * dt_PhiR(i,j,k) 

            dt_PhiI(i,j,k) =   x(i) * dx_PhiI(i,j,k) 
     &                      + y(j) * dy_PhiI(i,j,k) 
     &                      + z(k) * dz_PhiI(i,j,k)  
     &                      + PhiI_falloff 
     &                            * (   PhiI(i,j,k) 
     &                                - PhiI_asymptotic ) 
            dt_PhiI(i,j,k) = - inv_r * dt_PhiI(i,j,k) 

            dt_PiR(i,j,k) =   x(i) * dx_PiR(i,j,k) 
     &                     + y(j) * dy_PiR(i,j,k) 
     &                     + z(k) * dz_PiR(i,j,k)  
     &                     + PiR_falloff 
     &                            * (   PiR(i,j,k) 
     &                                - PiR_asymptotic ) 
            dt_PiR(i,j,k) = - inv_r * dt_PiR(i,j,k) 

            dt_PiI(i,j,k) =   x(i) * dx_PiI(i,j,k) 
     &                     + y(j) * dy_PiI(i,j,k) 
     &                     + z(k) * dz_PiI(i,j,k)  
     &                     + PiI_falloff 
     &                            * (   PiI(i,j,k) 
     &                                - PiI_asymptotic ) 
            dt_PiI(i,j,k) = - inv_r * dt_PiI(i,j,k) 
            end if
            end if
        end do 
        end do 
        end do 


! F3 
        do k = 1, nz  ,  nz-1
!        do j = 2, ny-1       
!        do i = 2, nx-1 
        do j = 1, ny
        do i = 1, nx
           ichr = nint(chr(i,j,k))
           if (ichr .ne. 8 .and. ichr .ne. 1) then

           r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)
           if (r .gt. 1.0d-1) then
             inv_r = 1.0d0 / r
           else
             inv_r = 1.0d0
             if (ltrace) 
     &         write(0,*)'calc_bssn_bcs_rhs: rsq < 0.1. chr=',chr(i,j,k)
           end if

           if (evolve_geometry .eq. 1) then

           dt_Alpha3D(i,j,k) =   x(i) * dx_alpha(i,j,k) 
     &                         + y(j) * dy_alpha(i,j,k) 
     &                         + z(k) * dz_alpha(i,j,k)  
     &                         + Alpha3D_falloff 
     &                             * (   Alpha3D(i,j,k)  
     &                                 - Alpha3D_asymptotic ) 
           dt_Alpha3D(i,j,k) = - inv_r * dt_Alpha3D(i,j,k) 

           dt_shiftx(i,j,k) =    x(i) * dx_shiftx(i,j,k) 
     &                         + y(j) * dy_shiftx(i,j,k) 
     &                         + z(k) * dz_shiftx(i,j,k)  
     &                         + shiftx_falloff 
     &                             * (   shiftx(i,j,k)  
     &                                 - shiftx_asymptotic ) 
           dt_shiftx(i,j,k) = - inv_r * dt_shiftx(i,j,k) 

           dt_shifty(i,j,k) =    x(i) * dx_shifty(i,j,k) 
     &                         + y(j) * dy_shifty(i,j,k) 
     &                         + z(k) * dz_shifty(i,j,k)  
     &                         + shifty_falloff 
     &                             * (   shifty(i,j,k)  
     &                                 - shifty_asymptotic ) 
           dt_shifty(i,j,k) = - inv_r * dt_shifty(i,j,k) 

           dt_shiftz(i,j,k) =    x(i) * dx_shiftz(i,j,k) 
     &                         + y(j) * dy_shiftz(i,j,k) 
     &                         + z(k) * dz_shiftz(i,j,k)  
     &                         + shiftz_falloff 
     &                             * (   shiftz(i,j,k) 
     &                                 - shiftz_asymptotic ) 
           dt_shiftz(i,j,k) = - inv_r * dt_shiftz(i,j,k) 

           dt_chi3D(i,j,k) =    x(i) * dx_chi(i,j,k) 
     &                        + y(j) * dy_chi(i,j,k) 
     &                        + z(k) * dz_chi(i,j,k)  
     &                        + chi3D_falloff 
     &                            * (   chi3D(i,j,k) 
     &                                - chi3D_asymptotic ) 
           dt_chi3D(i,j,k) = - inv_r * dt_chi3D(i,j,k) 

           dt_trK3D(i,j,k) =    x(i) * dx_trK(i,j,k) 
     &                        + y(j) * dy_trK(i,j,k) 
     &                        + z(k) * dz_trK(i,j,k)  
     &                        + trK3D_falloff 
     &                            * (   trK3D(i,j,k) 
     &                                - trK3D_asymptotic ) 
           dt_trK3D(i,j,k) = - inv_r * dt_trK3D(i,j,k) 

           dt_gtxx(i,j,k) =    x(i) * dx_gtxx(i,j,k) 
     &                       + y(j) * dy_gtxx(i,j,k) 
     &                       + z(k) * dz_gtxx(i,j,k)  
     &                       + gtxx_falloff 
     &                           * (   gtxx(i,j,k) 
     &                               - gtxx_asymptotic ) 
           dt_gtxx(i,j,k) = - inv_r * dt_gtxx(i,j,k) 

           dt_gtxy(i,j,k) =    x(i) * dx_gtxy(i,j,k) 
     &                       + y(j) * dy_gtxy(i,j,k) 
     &                       + z(k) * dz_gtxy(i,j,k)  
     &                       + gtxy_falloff 
     &                           * (   gtxy(i,j,k) 
     &                               - gtxy_asymptotic ) 
           dt_gtxy(i,j,k) = - inv_r * dt_gtxy(i,j,k) 

           dt_gtxz(i,j,k) =    x(i) * dx_gtxz(i,j,k) 
     &                       + y(j) * dy_gtxz(i,j,k) 
     &                       + z(k) * dz_gtxz(i,j,k)  
     &                       + gtxz_falloff 
     &                           * (   gtxz(i,j,k) 
     &                               - gtxz_asymptotic ) 
           dt_gtxz(i,j,k) = - inv_r * dt_gtxz(i,j,k) 

           dt_gtyy(i,j,k) =    x(i) * dx_gtyy(i,j,k) 
     &                       + y(j) * dy_gtyy(i,j,k) 
     &                       + z(k) * dz_gtyy(i,j,k)  
     &                       + gtyy_falloff 
     &                           * (   gtyy(i,j,k) 
     &                               - gtyy_asymptotic ) 
           dt_gtyy(i,j,k) = - inv_r * dt_gtyy(i,j,k) 

           dt_gtyz(i,j,k) =    x(i) * dx_gtyz(i,j,k) 
     &                       + y(j) * dy_gtyz(i,j,k) 
     &                       + z(k) * dz_gtyz(i,j,k)  
     &                       + gtyz_falloff 
     &                           * (   gtyz(i,j,k) 
     &                               - gtyz_asymptotic ) 
           dt_gtyz(i,j,k) = - inv_r * dt_gtyz(i,j,k) 

           dt_gtzz(i,j,k) =    x(i) * dx_gtzz(i,j,k) 
     &                       + y(j) * dy_gtzz(i,j,k) 
     &                       + z(k) * dz_gtzz(i,j,k)  
     &                       + gtzz_falloff 
     &                           * (   gtzz(i,j,k) 
     &                               - gtzz_asymptotic ) 
           dt_gtzz(i,j,k) = - inv_r * dt_gtzz(i,j,k) 

           dt_Atxx(i,j,k) =    x(i) * dx_Atxx(i,j,k) 
     &                       + y(j) * dy_Atxx(i,j,k) 
     &                       + z(k) * dz_Atxx(i,j,k)  
     &                       + Atxx_falloff 
     &                           * (   Atxx(i,j,k) 
     &                               - Atxx_asymptotic ) 
           dt_Atxx(i,j,k) = - inv_r * dt_Atxx(i,j,k) 

           dt_Atxy(i,j,k) =    x(i) * dx_Atxy(i,j,k) 
     &                       + y(j) * dy_Atxy(i,j,k) 
     &                       + z(k) * dz_Atxy(i,j,k)  
     &                       + Atxy_falloff 
     &                           * (   Atxy(i,j,k) 
     &                               - Atxy_asymptotic ) 
           dt_Atxy(i,j,k) = - inv_r * dt_Atxy(i,j,k) 

           dt_Atxz(i,j,k) =    x(i) * dx_Atxz(i,j,k) 
     &                       + y(j) * dy_Atxz(i,j,k) 
     &                       + z(k) * dz_Atxz(i,j,k)  
     &                       + Atxz_falloff 
     &                           * (   Atxz(i,j,k) 
     &                               - Atxz_asymptotic ) 
           dt_Atxz(i,j,k) = - inv_r * dt_Atxz(i,j,k) 

           dt_Atyy(i,j,k) =    x(i) * dx_Atyy(i,j,k) 
     &                       + y(j) * dy_Atyy(i,j,k) 
     &                       + z(k) * dz_Atyy(i,j,k)  
     &                       + Atyy_falloff 
     &                           * (   Atyy(i,j,k) 
     &                               - Atyy_asymptotic ) 
           dt_Atyy(i,j,k) = - inv_r * dt_Atyy(i,j,k) 

           dt_Atyz(i,j,k) =    x(i) * dx_Atyz(i,j,k) 
     &                       + y(j) * dy_Atyz(i,j,k) 
     &                       + z(k) * dz_Atyz(i,j,k)  
     &                       + Atyz_falloff 
     &                           * (   Atyz(i,j,k) 
     &                               - Atyz_asymptotic ) 
           dt_Atyz(i,j,k) = - inv_r * dt_Atyz(i,j,k) 

           dt_Atzz(i,j,k) =    x(i) * dx_Atzz(i,j,k) 
     &                       + y(j) * dy_Atzz(i,j,k) 
     &                       + z(k) * dz_Atzz(i,j,k)  
     &                       + Atzz_falloff 
     &                           * (   Atzz(i,j,k) 
     &                               - Atzz_asymptotic ) 
           dt_Atzz(i,j,k) = - inv_r * dt_Atzz(i,j,k) 

           dt_Gamtx(i,j,k) =    x(i) * dx_Gamtx(i,j,k) 
     &                        + y(j) * dy_Gamtx(i,j,k) 
     &                        + z(k) * dz_Gamtx(i,j,k)  
     &                        + Gamtx_falloff 
     &                            * (   Gamtx(i,j,k) 
     &                                - Gamtx_asymptotic ) 
           dt_Gamtx(i,j,k) = - inv_r * dt_Gamtx(i,j,k) 

           dt_Gamty(i,j,k) =    x(i) * dx_Gamty(i,j,k) 
     &                        + y(j) * dy_Gamty(i,j,k) 
     &                        + z(k) * dz_Gamty(i,j,k)  
     &                        + Gamty_falloff 
     &                            * (   Gamty(i,j,k) 
     &                                - Gamty_asymptotic ) 
           dt_Gamty(i,j,k) = - inv_r * dt_Gamty(i,j,k) 

           dt_Gamtz(i,j,k) =    x(i) * dx_Gamtz(i,j,k) 
     &                        + y(j) * dy_Gamtz(i,j,k) 
     &                        + z(k) * dz_Gamtz(i,j,k)  
     &                        + Gamtz_falloff 
     &                            * (   Gamtz(i,j,k) 
     &                                - Gamtz_asymptotic ) 
           dt_Gamtz(i,j,k) = - inv_r * dt_Gamtz(i,j,k) 

           dt_gbx(i,j,k) =      x(i) * dx_gbx(i,j,k) 
     &                        + y(j) * dy_gbx(i,j,k) 
     &                        + z(k) * dz_gbx(i,j,k)  
     &                        + gbx_falloff 
     &                            * (   gbx(i,j,k) 
     &                                - gbx_asymptotic ) 
           dt_gbx(i,j,k) = - inv_r * dt_gbx(i,j,k) 

           dt_gby(i,j,k) =      x(i) * dx_gby(i,j,k) 
     &                        + y(j) * dy_gby(i,j,k) 
     &                        + z(k) * dz_gby(i,j,k)  
     &                        + gby_falloff 
     &                            * (   gby(i,j,k) 
     &                                - gby_asymptotic ) 
           dt_gby(i,j,k) = - inv_r * dt_gby(i,j,k) 

           dt_gbz(i,j,k) =      x(i) * dx_gbz(i,j,k) 
     &                        + y(j) * dy_gbz(i,j,k) 
     &                        + z(k) * dz_gbz(i,j,k)  
     &                        + gbz_falloff 
     &                            * (   gbz(i,j,k) 
     &                                - gbz_asymptotic ) 
           dt_gbz(i,j,k) = - inv_r * dt_gbz(i,j,k) 
           end if
           if (evolve_em_field .eq. 1) then
           dt_Ex(i,j,k) =      x(i) * dx_Ex(i,j,k) 
     &                        + y(j) * dy_Ex(i,j,k) 
     &                        + z(k) * dz_Ex(i,j,k)  
     &                        + Ex_falloff 
     &                            * (   Ex(i,j,k) 
     &                                - Ex_asymptotic ) 
           dt_Ex(i,j,k) = - inv_r * dt_Ex(i,j,k) 
           !dt_Ex(i,j,k) = 0.d0

           dt_Ey(i,j,k) =      x(i) * dx_Ey(i,j,k) 
     &                        + y(j) * dy_Ey(i,j,k) 
     &                        + z(k) * dz_Ey(i,j,k)  
     &                        + Ey_falloff 
     &                            * (   Ey(i,j,k) 
     &                                - Ey_asymptotic ) 
           dt_Ey(i,j,k) = - inv_r * dt_Ey(i,j,k) 
           !dt_Ey(i,j,k) = 0.d0

           dt_Ez(i,j,k) =      x(i) * dx_Ez(i,j,k) 
     &                        + y(j) * dy_Ez(i,j,k) 
     &                        + z(k) * dz_Ez(i,j,k)  
     &                        + Ez_falloff 
     &                            * (   Ez(i,j,k) 
     &                                - Ez_asymptotic ) 
           dt_Ez(i,j,k) = - inv_r * dt_Ez(i,j,k) 
          ! dt_Ez(i,j,k) = 0.d0

           dt_Bx(i,j,k) =      x(i) * dx_Bx(i,j,k) 
     &                        + y(j) * dy_Bx(i,j,k) 
     &                        + z(k) * dz_Bx(i,j,k)  
     &                        + Bx_falloff 
     &                            * (   Bx(i,j,k) 
     &                                - Bx_asymptotic ) 
           dt_Bx(i,j,k) = - inv_r * dt_Bx(i,j,k) 
           !dt_Bx(i,j,k) = 0.d0 

           dt_By(i,j,k) =      x(i) * dx_By(i,j,k) 
     &                        + y(j) * dy_By(i,j,k) 
     &                        + z(k) * dz_By(i,j,k)  
     &                        + By_falloff 
     &                            * (   By(i,j,k) 
     &                                - By_asymptotic ) 
           dt_By(i,j,k) = - inv_r * dt_By(i,j,k) 
           !dt_By(i,j,k) = 0.d0 

           dt_Bz(i,j,k) =      x(i) * dx_Bz(i,j,k) 
     &                        + y(j) * dy_Bz(i,j,k) 
     &                        + z(k) * dz_Bz(i,j,k)  
     &                        + Bz_falloff 
     &                            * (   Bz(i,j,k) 
     &                                - Bz_asymptotic ) 
           dt_Bz(i,j,k) = - inv_r * dt_Bz(i,j,k) 
           !dt_Bz(i,j,k) = 0.d0 

           dt_Phi_em(i,j,k) =      x(i) * dx_Phi_em(i,j,k) 
     &                        + y(j) * dy_Phi_em(i,j,k) 
     &                        + z(k) * dz_Phi_em(i,j,k)  
     &                        + Phi_em_falloff 
     &                            * (   Phi_em(i,j,k) 
     &                                - Phi_em_asymptotic ) 
           dt_Phi_em(i,j,k) = - inv_r * dt_Phi_em(i,j,k) 
           !dt_Phi_em(i,j,k) = 0.d0

           dt_Psi_em(i,j,k) =      x(i) * dx_Psi_em(i,j,k) 
     &                        + y(j) * dy_Psi_em(i,j,k) 
     &                        + z(k) * dz_Psi_em(i,j,k)  
     &                        + Psi_em_falloff 
     &                            * (   Psi_em(i,j,k) 
     &                                - Psi_em_asymptotic ) 
           dt_Psi_em(i,j,k) = - inv_r * dt_Psi_em(i,j,k) 
           !dt_Psi_em(i,j,k) = 0.d0 
           end if                  
           if (evolve_scalar_field .eq. 1) then
           dt_PhiR(i,j,k) =   x(i) * dx_PhiR(i,j,k) 
     &                      + y(j) * dy_PhiR(i,j,k) 
     &                      + z(k) * dz_PhiR(i,j,k)  
     &                      + PhiR_falloff 
     &                            * (   PhiR(i,j,k) 
     &                                - PhiR_asymptotic ) 
            dt_PhiR(i,j,k) = - inv_r * dt_PhiR(i,j,k) 

            dt_PhiI(i,j,k) =   x(i) * dx_PhiI(i,j,k) 
     &                      + y(j) * dy_PhiI(i,j,k) 
     &                      + z(k) * dz_PhiI(i,j,k)  
     &                      + PhiI_falloff 
     &                            * (   PhiI(i,j,k) 
     &                                - PhiI_asymptotic ) 
            dt_PhiI(i,j,k) = - inv_r * dt_PhiI(i,j,k) 

            dt_PiR(i,j,k) =   x(i) * dx_PiR(i,j,k) 
     &                     + y(j) * dy_PiR(i,j,k) 
     &                     + z(k) * dz_PiR(i,j,k)  
     &                     + PiR_falloff 
     &                            * (   PiR(i,j,k) 
     &                                - PiR_asymptotic ) 
            dt_PiR(i,j,k) = - inv_r * dt_PiR(i,j,k) 

           dt_PiI(i,j,k) =   x(i) * dx_PiI(i,j,k) 
     &                     + y(j) * dy_PiI(i,j,k) 
     &                     + z(k) * dz_PiI(i,j,k)  
     &                     + PiI_falloff 
     &                            * (   PiI(i,j,k) 
     &                                - PiI_asymptotic ) 
            dt_PiI(i,j,k) = - inv_r * dt_PiI(i,j,k) 
            end if
            end if                  
        end do 
        end do 
        end do 



        return
      end subroutine

#if 1
      !----------------------------------------------------------------
      !
      !
      !----------------------------------------------------------------
      subroutine bssn_evil_bcs( dtu, u, v,w, imask, 
     &                          x, y, z, nx, ny, nz, par ) 
        use GF
        use params
        implicit   none

        type(gridfunction), dimension(NU_G)     :: dtu 
        type(gridfunction), dimension(NU_G)     :: u 
        type(gridfunction), dimension(NV_G)     :: v 
        type(gridfunction), dimension(NW)       :: w 
        
        CCTK_INT , dimension(nx,ny,nz)          :: imask  
        CCTK_REAL, dimension(NPAR), intent(in)  :: par 
        CCTK_INT                                :: nx, ny, nz
        CCTK_REAL, dimension(nx), intent(in)    :: x
        CCTK_REAL, dimension(ny), intent(in)    :: y
        CCTK_REAL, dimension(nz), intent(in)    :: z

        ! local variables 
        CCTK_INT    ::  i, j, k , m, ichr
        CCTK_REAL   ::  inv_r , r
        CCTK_REAL   ::  bc_r_min, val, val1, val2
        CCTK_REAL, pointer, dimension(:,:,:)    :: chr

        logical  , parameter     ::   ltrace = .false.
        logical  , parameter     ::   periodic = .true.

 
        chr  => w(G_CHR)%d

        inv_r = -1.0d0
        bc_r_min = 1.0d-6
        val1 =  0.0d0
!        val1 =  sqrt(inv_r)
        val2 =  sqrt(inv_r)

        ! periodic boundary conditions
        if (periodic) then

           print*, "periodic boundary conditions"

           do m = 1, NU_G
              u(m)%d(1,:,:)    = u(m)%d(nx-7,:,:)
              u(m)%d(2,:,:)    = u(m)%d(nx-6,:,:)
              u(m)%d(3,:,:)    = u(m)%d(nx-5,:,:)
              u(m)%d(4,:,:)    = u(m)%d(nx-4,:,:)
              u(m)%d(nx-3,:,:) = u(m)%d(5,:,:)
              u(m)%d(nx-2,:,:) = u(m)%d(6,:,:)
              u(m)%d(nx-1,:,:) = u(m)%d(7,:,:)
              u(m)%d(nx,:,:)   = u(m)%d(8,:,:)

              u(m)%d(:,1,:)    = u(m)%d(:,ny-7,:)
              u(m)%d(:,2,:)    = u(m)%d(:,ny-6,:)
              u(m)%d(:,3,:)    = u(m)%d(:,ny-5,:)
              u(m)%d(:,4,:)    = u(m)%d(:,ny-4,:)
              u(m)%d(:,ny-3,:) = u(m)%d(:,5,:)
              u(m)%d(:,ny-2,:) = u(m)%d(:,6,:)
              u(m)%d(:,ny-1,:) = u(m)%d(:,7,:)
              u(m)%d(:,ny,:)   = u(m)%d(:,8,:)

              u(m)%d(:,:,1)    = u(m)%d(:,:,nz-7)
              u(m)%d(:,:,2)    = u(m)%d(:,:,nz-6)
              u(m)%d(:,:,3)    = u(m)%d(:,:,nz-5)
              u(m)%d(:,:,4)    = u(m)%d(:,:,nz-4)
              u(m)%d(:,:,nz-3) = u(m)%d(:,:,5)
              u(m)%d(:,:,nz-2) = u(m)%d(:,:,6)
              u(m)%d(:,:,nz-1) = u(m)%d(:,:,7)
              u(m)%d(:,:,nz)   = u(m)%d(:,:,8)

              dtu(m)%d(1,:,:)    = dtu(m)%d(nx-7,:,:)
              dtu(m)%d(2,:,:)    = dtu(m)%d(nx-6,:,:)
              dtu(m)%d(3,:,:)    = dtu(m)%d(nx-5,:,:)
              dtu(m)%d(4,:,:)    = dtu(m)%d(nx-4,:,:)
              dtu(m)%d(nx-3,:,:) = dtu(m)%d(5,:,:)
              dtu(m)%d(nx-2,:,:) = dtu(m)%d(6,:,:)
              dtu(m)%d(nx-1,:,:) = dtu(m)%d(7,:,:)
              dtu(m)%d(nx,:,:)   = dtu(m)%d(8,:,:)

              dtu(m)%d(:,1,:)    = dtu(m)%d(:,ny-7,:)
              dtu(m)%d(:,2,:)    = dtu(m)%d(:,ny-6,:)
              dtu(m)%d(:,3,:)    = dtu(m)%d(:,ny-5,:)
              dtu(m)%d(:,4,:)    = dtu(m)%d(:,ny-4,:)
              dtu(m)%d(:,ny-3,:) = dtu(m)%d(:,5,:)
              dtu(m)%d(:,ny-2,:) = dtu(m)%d(:,6,:)
              dtu(m)%d(:,ny-1,:) = dtu(m)%d(:,7,:)
              dtu(m)%d(:,ny,:)   = dtu(m)%d(:,8,:)

              dtu(m)%d(:,:,1)    = dtu(m)%d(:,:,nz-7)
              dtu(m)%d(:,:,2)    = dtu(m)%d(:,:,nz-6)
              dtu(m)%d(:,:,3)    = dtu(m)%d(:,:,nz-5)
              dtu(m)%d(:,:,4)    = dtu(m)%d(:,:,nz-4)
              dtu(m)%d(:,:,nz-3) = dtu(m)%d(:,:,5)
              dtu(m)%d(:,:,nz-2) = dtu(m)%d(:,:,6)
              dtu(m)%d(:,:,nz-1) = dtu(m)%d(:,:,7)
              dtu(m)%d(:,:,nz)   = dtu(m)%d(:,:,8)

           end do
        


        ! evil ones
        else

        do k = 1, nz
        do j = 1, ny 
        do i = 1, nx, nx-1  
!          r = sqrt(x(i)**2 + y(j)**2 + z(k)**2) 
!          if (r .lt. bc_r_min) then
!            val = val1
!          else
!            val = val2
!          end if

           ichr = nint(chr(i,j,k))
           if (ichr .ne. 0 .and. ichr .ne. 1 .and. ichr .ne. 8) then
             val = 0.d0
           else
             val = val1
           end if

           do m = 1, NU_G
              dtu(m)%d(i,j,k) = val
           end do
        end do
        end do
        end do

        do k = 1, nz
        do j = 1, ny, ny-1 
        do i = 1, nx
!          r = sqrt(x(i)**2 + y(j)**2 + z(k)**2) 
!          if (r .lt. bc_r_min) then
!            val = val1
!          else
!            val = val2
!          end if

           ichr = nint(chr(i,j,k))
           if (ichr .ne. 0 .and. ichr .ne. 1 .and. ichr .ne. 8) then
             val = 0.d0
           else
             val = val1
           end if
           do m = 1, NU_G
              dtu(m)%d(i,j,k) = val
           end do
        end do
        end do
        end do

        do k = 1, nz, nz-1
        do j = 1, ny
        do i = 1, nx
!          r = sqrt(x(i)**2 + y(j)**2 + z(k)**2) 
!          if (r .lt. bc_r_min) then
!            val = val1
!          else
!            val = val2
!          end if

           ichr = nint(chr(i,j,k))
           if (ichr .ne. 0 .and. ichr .ne. 1 .and. ichr .ne. 8) then
             val = 0.d0
           else
             val = val1
           end if

           do m = 1, NU_G
              dtu(m)%d(i,j,k) = val
           end do
        end do
        end do
        end do

        end if ! end of periodic/evil boundary conditions

        return
      end subroutine
#endif 



