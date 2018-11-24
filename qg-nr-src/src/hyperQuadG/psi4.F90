#include "cctk.h"
!-----------------------------------------------------------------------
!       
!   This routine calculates the real and imaginary parts of the Weyl     
!   scalar, Psi_4.  This is done with respect to a tetrad calculated     
!   in calc_tetrad (clever, eh?).    
!       
!   We should point out that an assumption is made that this will be  
!   calculated in a vacuum region.  If not, there are additional matter 
!   terms in the calculation that are being neglected in this.   
!    
!-----------------------------------------------------------------------
      subroutine calc_psi4(u, v, w,
#include "CR_con_sub_args.h"
     &                    nx, ny, nz, par)
        
        use GF
        use params 
        implicit none
        
        type(gridfunction), dimension(NU_G)               :: dtu
        type(gridfunction), dimension(NU_G)               :: u
        type(gridfunction), dimension(NV_G)               :: v
        type(gridfunction), dimension(NW)                 :: w
        CCTK_REAL, dimension(NPAR),intent(in)             :: par
        CCTK_INT                                          :: nx, ny, nz
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
        CCTK_REAL, dimension(:,:,:), pointer :: gxx 
        CCTK_REAL, dimension(:,:,:), pointer :: gxy 
        CCTK_REAL, dimension(:,:,:), pointer :: gxz
        CCTK_REAL, dimension(:,:,:), pointer :: gyy 
        CCTK_REAL, dimension(:,:,:), pointer :: gyz 
        CCTK_REAL, dimension(:,:,:), pointer :: gzz
        CCTK_REAL, dimension(:,:,:), pointer :: gbx
        CCTK_REAL, dimension(:,:,:), pointer :: gby
        CCTK_REAL, dimension(:,:,:), pointer :: gbz
        CCTK_REAL, dimension(:,:,:), pointer :: psi4r
        CCTK_REAL, dimension(:,:,:), pointer :: psi4i
        CCTK_REAL, dimension(:,:,:), pointer :: x3d
        CCTK_REAL, dimension(:,:,:), pointer :: y3d
        CCTK_REAL, dimension(:,:,:), pointer :: z3d
        CCTK_REAL, dimension(3)              :: xpt

!use phi0R and phi0I for Phi22 and flux of scalar field
        CCTK_REAL, dimension(:,:,:), pointer :: phi0R
        CCTK_REAL, dimension(:,:,:), pointer :: phi0I
        CCTK_REAL, dimension(:,:,:), pointer :: phiR
        CCTK_REAL, dimension(:,:,:), pointer :: piR
        
        CCTK_INT                       :: i, j, k

        CCTK_REAL                      :: dx, dy, dz 
        CCTK_REAL                      :: trK, chi, G, inv_chi
        CCTK_REAL,  dimension(3)       :: CalGamt, Gamt
        CCTK_REAL,  dimension(3)       :: d_chi
        CCTK_REAL,  dimension(3,3)     :: gd, gu, gtd, gtu
        CCTK_REAL,  dimension(3,3)     :: Atd, Atu, Atud, dd_chi
        CCTK_REAL,  dimension(3,3)     :: Rpd, Rtd, R
        CCTK_REAL,  dimension(3,3)     :: d_Gamt
        CCTK_REAL,  dimension(3,3)     :: DelDel_chi
        CCTK_REAL,  dimension(3,3,3)   :: d_gtd, d_Atd
        CCTK_REAL,  dimension(3,3,3,3) :: dd_gtd
        CCTK_REAL                      :: detgtd, idetgtd
        CCTK_REAL,  dimension(3,3,3)   :: Ctd, Ct, Chr
!        CCTK_REAL                      :: exp4chi, expm4chi 
        CCTK_REAL                      :: r_coord 
        CCTK_REAL                      :: psi4_real, psi4_imag, 
     &                                    dil_alpha
!       maple vars start (need to modify)
        CCTK_REAL                      :: t31, t64, t67, t70, t143
        CCTK_REAL                      :: t145, t148, t177, t201, t272
        CCTK_REAL                      :: t274, t294, t296, t316, t399
        CCTK_REAL                      :: t402, t405, t429, t498, t501
        CCTK_REAL                      :: t503, t519
!       maple vars end

        CCTK_REAL, dimension(3)        :: r_np, m_np_real, m_np_imag
        CCTK_REAL, dimension(0:3)      :: r_u, m_r_u, m_i_u

        CCTK_REAL, parameter      :: one = 1.0d0
        CCTK_REAL, parameter      :: two = 2.0d0
        CCTK_REAL, parameter      :: threehalves = 1.5d0
        CCTK_REAL, parameter      :: half = 0.5d0
        CCTK_REAL, parameter      :: third = 0.3333333333333333d0
        CCTK_REAL, parameter      :: fourth = 0.25d0

        ! include declarations for Maple temp variables
#include "psi4_temp_vars.h"

        ! G is Newton"s gravitation constant
        G = 1.0d0

        dx = par(P_DX)
        dy = par(P_DY)
        dz = par(P_DZ)
	dil_alpha = par(P_DIL_ALPHA)

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

        gxx         => v(G_G11)%d
        gxy         => v(G_G12)%d
        gxz         => v(G_G13)%d
        gyy         => v(G_G22)%d
        gyz         => v(G_G23)%d
        gzz         => v(G_G33)%d
        
	phiR        => u(G_PHIR)%d
	piR         => u(G_PIR)%d
        psi4R       => v(G_PSI4R)%d
        psi4I       => v(G_PSI4I)%d
        phi0R       => v(G_PHI0R)%d
        phi0I       => v(G_PHI0I)%d
        
        do k = 2, nz-1 
        do j = 2, ny-1 
        do i = 2, nx-1 

          xpt(1) = x3d(i,j,k)
          xpt(2) = y3d(i,j,k)
          xpt(3) = z3d(i,j,k)

          gd(1,1) = gxx(i,j,k)
          gd(1,2) = gxy(i,j,k)
          gd(1,3) = gxz(i,j,k)
          gd(2,1) = gxy(i,j,k)
          gd(2,2) = gyy(i,j,k)
          gd(2,3) = gyz(i,j,k)
          gd(3,1) = gxz(i,j,k)
          gd(3,2) = gyz(i,j,k)
          gd(3,3) = gzz(i,j,k)
       
          call calc_tetrad( r_u, m_r_u, m_i_u, gd, xpt )

          r_np(1)    =  r_u(1)
          r_np(2)    =  r_u(2) 
          r_np(3)    =  r_u(3) 
          m_np_real(1)  =  m_r_u(1) 
          m_np_real(2)  =  m_r_u(2) 
          m_np_real(3)  =  m_r_u(3)
          m_np_imag(1)  =  m_i_u(1) 
          m_np_imag(2)  =  m_i_u(2)
          m_np_imag(3)  =  m_i_u(3) 


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
        
          d_chi(1)      = dx_chi(i,j,k)
          d_chi(2)      = dy_chi(i,j,k)
          d_chi(3)      = dz_chi(i,j,k)

          dd_chi(1,1)   = dxx_chi(i,j,k)
          dd_chi(1,2)   = dxy_chi(i,j,k)
          dd_chi(1,3)   = dxz_chi(i,j,k)
          dd_chi(2,2)   = dyy_chi(i,j,k)
          dd_chi(2,3)   = dyz_chi(i,j,k)
          dd_chi(3,3)   = dzz_chi(i,j,k)

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

#include "psi4.h" 

          r_coord = sqrt(    xpt(1)*xpt(1)
     &                     + xpt(2)*xpt(2)
     &                     + xpt(3)*xpt(3)  )

          psi4R(i,j,k)  =  0.25 * r_coord * psi4_real
          psi4I(i,j,k)  =  0.25 * r_coord * psi4_imag


!I can also get psi22 which i will keep in phi0R
!i will use a hack where i assume christoffels go like 1/r at best
!so i can ignore them for contributing a net 1/r^2. I also
!assume phi(t-r)/r, so phi,t = - phi,r and so phi22 = - phi,tr = -pi,r
	  phi0R(i,j,k) = 4.*dil_alpha*exp(-2.*dil_alpha*phiR(i,j,k)) * 
     &             (r_np(1)*dx_piR(i,j,k) 
     &            + r_np(2)*dy_piR(i,j,k)
     &            + r_np(3)*dz_piR(i,j,k) )* r_coord

	  phi0I(i,j,k) = piR(i,j,k) * 
     &             (r_np(1)*dx_piR(i,j,k) 
     &            + r_np(2)*dy_piR(i,j,k)
     &            + r_np(3)*dz_piR(i,j,k) )* r_coord**2

!last, we also need to correct to the jordan frame psi4
          psi4R(i,j,k)  = psi4R(i,j,k) * exp(-2.*dil_alpha*phiR(i,j,k))
          psi4I(i,j,k)  = psi4I(i,j,k) * exp(-2.*dil_alpha*phiR(i,j,k))

        end do
        end do
        end do
          
          
          
        return
      end  subroutine
          

#include "cctk.h"

      !----------------------------------------------------------------
      !
      !  This calculates the tetrad used in the Psi_4 calculation.  
      !  
      !  In an early version of this, we calculated the full tetrad:  
      !  (l,u,m,mbar) where l=n-r, u=n+r, m=theta+iphi, etc and where 
      !  n is the usual normal to the foliation and (r,theta,phi) are  
      !  an orthonormal set with respect to the 3-metric.  However, we 
      !  now use a calculation of Psi_4 which has been simplified and 
      !  only uses r and mbar.  As a result, we don't need the full 
      !  tetrad.  Comments below highlight those portions no longer 
      !  explicitly needed.   
      !  
      !  Finally, note that the tetrad is *not* actually normalized.   
      !  In the definition here, we have left off the usual factor of 
      !  1/sqrt(2).  This is okay because we put in an extra factor of   
      !  1/4 = ( 1/sqrt(2) )^4 in the definition of Psi_4.   
      !  
      !----------------------------------------------------------------
      subroutine calc_tetrad( r_u, m_r_u, m_i_u, hd, xpt) 
        implicit none

        CCTK_REAL, dimension(0:3) ::  r_u, m_r_u, m_i_u 
        CCTK_REAL, dimension(3,3) ::  hd
        CCTK_REAL, dimension(3)   ::  xpt

        ! local vars
        CCTK_INT                  ::  a, b
        CCTK_REAL, dimension(3)   ::  theta, phi
        CCTK_REAL                 ::  x, y, z
        CCTK_REAL                 ::  inn_prod,  inv_inn_prod
        CCTK_REAL                 ::  inn_prod1, inn_prod2 
       
        CCTK_REAL                 ::  invsqrt2 = 0.7071067811865475244d0 
       
        ! set the initial set of independent spatial vectors 
        x = xpt(1)
        y = xpt(2)
        z = xpt(3)
      
        theta(1) = x * z
        theta(2) = y * z
        theta(3) = - ( x*x + y*y ) 

        phi(1) = - y 
        phi(2) =   x 
        phi(3) =   0.d0 


        if (   (abs(x) .gt. 1.d-5) .or.
     &         (abs(y) .gt. 1.d-5) .or.
     &         (abs(z) .gt. 1.d-5)       ) then



!        if ( (abs(x) .lt. 1.d-7) .and. (abs(y) .lt. 1.d-7) ) then   
!
!          r_u(0) = 0.d0   
!          r_u(1) = 0.d0   
!          r_u(2) = 0.d0   
!          r_u(3) = 0.d0   
!
!          m_r_u(0) = 0.d0   
!          m_r_u(1) = 0.d0   
!          m_r_u(2) = 0.d0   
!          m_r_u(3) = 0.d0   
!
!          m_i_u(0) = 0.d0   
!          m_i_u(1) = 0.d0   
!          m_i_u(2) = 0.d0   
!          m_i_u(3) = 0.d0   
!        end if 
        ! FIXME Question about logic here.  Should this be an if-then/else
        !       loop?

          inn_prod = 0.d0 
          do a = 1, 3
            inn_prod = inn_prod + hd(a,a)*xpt(a)*xpt(a) 
          end do 

          do a = 1  , 2 
          do b = a+1, 3 
            inn_prod = inn_prod + 2*hd(a,b)*xpt(a)*xpt(b) 
          end do
          end do

          inv_inn_prod = 1.d0 / sqrt( inn_prod )  

          do a = 1, 3
            xpt(a) = xpt(a) * inv_inn_prod 
          end do 

          inn_prod  = 0.d0 
          inn_prod1 = 0.d0 
          do a = 1, 3  
          do b = 1, 3  
            inn_prod  = inn_prod  + hd(a,b)*theta(a)*theta(b) 
            inn_prod1 = inn_prod1 + hd(a,b)*  xpt(a)*theta(b) 
          end do 
          end do 

          inv_inn_prod = 1.d0 / sqrt( inn_prod - inn_prod1*inn_prod1 )  

          do a = 1, 3
            theta(a) = (theta(a) - inn_prod1*xpt(a)) * inv_inn_prod 
          end do 

          inn_prod  = 0.d0 
          inn_prod1 = 0.d0 
          inn_prod2 = 0.d0 
          do a = 1, 3 
          do b = 1, 3 
             inn_prod  = inn_prod  + hd(a,b)*  phi(a)*phi(b) 
             inn_prod1 = inn_prod1 + hd(a,b)*  xpt(a)*phi(b) 
             inn_prod2 = inn_prod2 + hd(a,b)*theta(a)*phi(b) 
          end do 
          end do  

          inv_inn_prod = 1.d0 / sqrt(  inn_prod  
     &                               - (   inn_prod1*inn_prod1
     &                                   + inn_prod2*inn_prod2 ) 
     &                              )  

          do a = 1, 3
             phi(a) = (                   phi(a) 
     &                  - (   inn_prod1*  xpt(a) 
     &                      + inn_prod2*theta(a) ) 
     &                ) * inv_inn_prod
          end do 

   
          r_u(0) = 0.d0  
          r_u(1) = xpt(1)  
          r_u(2) = xpt(2)  
          r_u(3) = xpt(3)  
  
          m_r_u(0) = 0.d0  
          m_r_u(1) = theta(1)  
          m_r_u(2) = theta(2)  
          m_r_u(3) = theta(3)  
  
          m_i_u(0) = 0.d0  
          m_i_u(1) = phi(1)  
          m_i_u(2) = phi(2)  
          m_i_u(3) = phi(3)  

        else

          r_u(0) = 0.d0
          r_u(1) = 0.d0
          r_u(2) = 0.d0
          r_u(3) = 0.d0

          m_r_u(0) = 0.d0
          m_r_u(1) = 0.d0
          m_r_u(2) = 0.d0
          m_r_u(3) = 0.d0

          m_i_u(0) = 0.d0
          m_i_u(1) = 0.d0
          m_i_u(2) = 0.d0
          m_i_u(3) = 0.d0
        end if

        return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! need to calculate phi2 - copied from compute_psi4.f90 in hyperGHEM

#include "cctk.h"

        subroutine calc_phi2(PHI2_R, PHI2_I, PHI0_R, PHI0_I, 
     &                       g, F_uu, l, n, m_R, m_I, x, y, z)
	implicit none
	real*8, dimension(4,4) :: g, F_uu, F_du, F_dd
	real*8, dimension(1:4) :: l, n, m_R, m_I
	real*8, dimension(1:4) :: ld, nd, md_R, md_I, r, t, f, u
	real*8 :: PHI2_R, PHI2_I, PHI0_R, PHI0_I, x, y, z
        real*8 :: rsph, rcyl, sq2
        !local vars
	integer:: la,lb,li
	
!	print *,'inside calc_phi2 - going well so far...'

        ! compute the tetrad with the index down
        ld = 0.0d0
        nd = 0.0d0
        md_R = 0.0d0
        md_I = 0.0d0

	do la=1,4
	do lb=1,4
	   ld(la)   = ld(la) + g(la,lb)*l(lb)
	   nd(la)   = nd(la) + g(la,lb)*n(lb)
	   md_R(la) = md_R(la) + g(la,lb)*m_R(lb)
	   md_I(la) = md_I(la) + g(la,lb)*m_I(lb)
	end do
	end do

        ! compute the phi2
	PHI2_R = 0.0d0
	PHI2_I = 0.0d0
	PHI0_R = 0.0d0
	PHI0_I = 0.0d0

	do la=1,4
	do lb=1,4
	   PHI2_R = PHI2_R + F_uu(la,lb)*md_R(la)*nd(lb)
	   PHI2_I = PHI2_I - F_uu(la,lb)*md_I(la)*nd(lb)
	   PHI0_R = PHI0_R + F_uu(la,lb)*ld(la)*md_R(lb)
	   PHI0_I = PHI0_I + F_uu(la,lb)*ld(la)*md_I(lb)
	end do
	end do

	return
	end subroutine calc_phi2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! need to calculate phi2 - copied from compute_psi4.f90 in hyperGHEM

#include "cctk.h"

	subroutine frame(l, n, m_R, m_I, x, y, z, alpha, beta, g, g_uu, deth)
	implicit none
	
	real*8 :: x, y, z, g(4,4), alpha, beta(3)
	real*8 :: l(4), n(4), m_R(4), m_I(4), g_uu(4,4), h_uu(3,3)
	real*8 :: u(4), rr(4), t(4), f(4), sq2
	
!local quantities for intermediate steps
	real*8 :: v(3,3),tw1,tw2,tw3, vn(3,3), fud
	real*8 :: eps(3,3,3),r, temp, deth,w(3,3)
	real*8 :: t2,t4,t7,t9,t12,t15,t20,t24,t30
	integer :: li,lj,lk,ld
	logical, parameter :: ltrace =.false.
	
!define the total antysmmetric tensor
       	eps = 0.0d0
	eps(1,2,3) = 1.0d0; eps(2,3,1) = 1.0d0; eps(3,1,2) = 1.0d0
	eps(1,3,2) = -1.0d0; eps(2,1,3) = -1.0d0; eps(3,2,1) = -1.0d0
	
	r = sqrt(x**2+y**2+z**2)
	fud = 0.0
        if( (abs(y).lt.1.e-5) .and. (abs(x).lt.1.e-5) ) then
	  if (ltrace) print*, "we are on the z axis; tetrad set to zero"
	  l = 0.0; n=0.0d0; m_R = 0.0d0; m_I = 0.0d0
	else

	v(1,1) = -y; v(1,2) = x; v(1,3) = 0.0d0
	v(2,1) = x; v(2,2) = y; v(2,3) = z
	
        ! which one is the right one?
!	do li=1,3	
!        do lj=1,3
!           h_uu(li,lj) = g_uu(li,lj) + beta(li)*beta(lj)/(alpha*alpha)
!        end do
!        end do

	do li=1,3	
	   temp=0.0
	   do lj=1,3
	     do lk=1,3
	       do ld=1,3
	         temp =temp + g_uu(li,lj)*eps(lj,lk,ld)*v(1,lk)*v(2,ld)
!	         temp =temp + h_uu(li,lj)*eps(lj,lk,ld)*v(1,lk)*v(2,ld)
	       end do
	     end do
	   end do	
	   v(3,li)=temp
	end do

!now redefine things
! vector 1

!now get the inner product to orthonormalize wrt the metric

        temp = 0.0d0
	do lk=1,3
	  do ld=1,3
	     temp = temp + g(lk,ld)*v(1,lk)*v(1,ld)	
	  end do
	end do
	w(1,1) = temp 
     	      
	if(w(1,1).lt.1.e-10) then
          l = 0.0; n=0.0d0; m_R = 0.0d0; m_I = 0.0d0
	  if (ltrace) print*, 'trouble with w11 in frame',
     &                         tw1,v(1,1),v(1,2),x,y,z
	  GOTO 100	  
!	  STOP
	  tw1 = 1.0
	else
          tw1 = sqrt(w(1,1))
	  tw1 = 1.0d0/tw1
	end if
	  
	v(1,:) = v(1,:) * tw1
	
! vector 2
	
!now get the inner product to orthonormalize wrt the metric
	do li=1,2
	  temp = 0.0d0
	  do lk=1,3
	     do ld=1,3
	       temp = temp + g(lk,ld)*v(li,lk)*v(2,ld)	
	     end do
	  end do
	  w(li,2) = temp
	end do

!orthogonalize	
	v(2,:) = (v(2,:) - v(1,:)*w(1,2)) 

!normalize correctly
         temp = 0.0d0
          do lk=1,3
             do ld=1,3
               temp = temp + g(lk,ld)*v(2,lk)*v(2,ld)
             end do
         end do
         w(2,2) = temp

        if(w(2,2).lt.1.e-10) then
          l = 0.0; n=0.0d0; m_R = 0.0d0; m_I = 0.0d0
          if (ltrace) print*, 'trouble with w22 in frame',
     &                         tw2,v(2,1),v(2,2)
          GOTO 100
          tw2 = 1.0
        else
          tw2 = sqrt(w(2,2))
          tw2 = 1.0d0/tw2
        end if
        v(2,:) = v(2,:) * tw2

! vector 3	
!now get the inner product to orthonormalize wrt the metric
	do li=1,3
	  temp = 0.0d0
	  do lk=1,3
	    do ld=1,3
	      temp = temp + g(lk,ld)*v(li,lk)*v(3,ld)	
	    end do
	  end do
	  w(li,3) = temp	
	end do

!orthogonalize	
	v(3,:) = (v(3,:) - v(1,:)*w(1,3) - v(2,:)*w(2,3)) 
!normalize correctly
        temp = 0.0d0
        do lk=1,3
           do ld=1,3
             temp = temp + g(lk,ld)*v(3,lk)*v(3,ld)
           end do
        end do
        w(3,3) = temp

        if(w(3,3).lt.1.e-10) then
          l = 0.0; n=0.0d0; m_R = 0.0d0; m_I = 0.0d0
          if (ltrace) print*, 'trouble with w33 in frame', x,y,z,w(3,3)
          GOTO 100
          tw3 = 1.0
        else
          tw3 = sqrt(w(3,3))
          tw3 = 1.0d0/tw3
        end if
	v(3,:) = v(3,:)*tw3
	
	
!wit these constrcut the tetrad	
	rr(1) = v(2,1); rr(2) = v(2,2); rr(3) = v(2,3); rr(4) = 0.0
	t(1) = v(3,1); t(2) = v(3,2); t(3) = v(3,3); t(4) = 0.0
	f(1) = v(1,1); f(2) = v(1,2); f(3) = v(1,3); f(4) = 0.0
	u(1) = -1./alpha*beta(1);u(2)=-1./alpha*beta(2)
	u(3) = -1./alpha*beta(3); u(4) = 1./alpha

	sq2 = 1./sqrt(2.0d0)
	
	do li=1,4
	  l(li) = sq2 * (u(li) + rr(li))
	  n(li) = sq2 * (u(li) - rr(li))
	  m_r(li) = sq2 * t(li)
	  m_i(li) = sq2 * f(li) 
	end do

100 	CONTINUE

	end if
	
	end subroutine frame	
