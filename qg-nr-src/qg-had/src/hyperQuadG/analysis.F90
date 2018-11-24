!--------------------------------------------------------------------
!
!  $Id: analysis.F90,v 1.11 2011-07-06 13:07:37 carlos Exp $
!
!--------------------------------------------------------------------

#include "cctk.h"

#if defined GR_MHD
      module mod_bssn_analysis
#else
      module mod_analysis
#endif

        use params
        use UTILEQS
        use GF
        use HYPER_DERIVS_42

       !  use MPSI4
       !  use MRHS
    
      implicit none
      CONTAINS

      !--------------------------------------------------------------------
      !
      !    GF analysis
      !
      !--------------------------------------------------------------------

#if defined GR_MHD
      subroutine bssn_analysis(u2, u0, dxu, dyu, dzu, v, 
     &                         dxv, dyv, dzv, w, par)
#else
      subroutine analysis(u2, u0, dxu, dyu, dzu, v, 
     &                         dxv, dyv, dzv, w, par)
#endif
        implicit none
        type(gridfunction), dimension(NU_G)   :: u2,u0
        type(gridfunction), dimension(NV_G)   :: v, dxv, dyv, dzv
        type(gridfunction), dimension(NW)     :: w
        type(gridfunction), dimension(NU_G)   :: dxu,dyu,dzu
        CCTK_REAL, dimension(:)               :: par

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
        CCTK_REAL, dimension(:,:,:), pointer :: chi 
        CCTK_REAL, dimension(:,:,:), pointer :: trK
        CCTK_REAL, dimension(:,:,:), pointer :: Gamtx 
        CCTK_REAL, dimension(:,:,:), pointer :: Gamty 
        CCTK_REAL, dimension(:,:,:), pointer :: Gamtz
!       CCTK_REAL, dimension(:,:,:), pointer :: Tmn00, Tmn01, Tmn02
!       CCTK_REAL, dimension(:,:,:), pointer :: Tmn03, Tmn11, Tmn12
!       CCTK_REAL, dimension(:,:,:), pointer :: Tmn13, Tmn22, Tmn23
!       CCTK_REAL, dimension(:,:,:), pointer :: Tmn33
        CCTK_REAL, dimension(:,:,:), pointer :: gxx 
        CCTK_REAL, dimension(:,:,:), pointer :: gxy 
        CCTK_REAL, dimension(:,:,:), pointer :: gxz
        CCTK_REAL, dimension(:,:,:), pointer :: gyy 
        CCTK_REAL, dimension(:,:,:), pointer :: gyz
        CCTK_REAL, dimension(:,:,:), pointer :: gzz
        CCTK_REAL, dimension(:,:,:), pointer :: gbx
        CCTK_REAL, dimension(:,:,:), pointer :: gby
        CCTK_REAL, dimension(:,:,:), pointer :: gbz

        CCTK_REAL, dimension(:,:,:), pointer :: Ex 
        CCTK_REAL, dimension(:,:,:), pointer :: Ey 
        CCTK_REAL, dimension(:,:,:), pointer :: Ez 
        CCTK_REAL, dimension(:,:,:), pointer :: Bx 
        CCTK_REAL, dimension(:,:,:), pointer :: By 
        CCTK_REAL, dimension(:,:,:), pointer :: Bz
        CCTK_REAL, dimension(:,:,:), pointer :: Psi_em 
        CCTK_REAL, dimension(:,:,:), pointer :: Phi_em
        CCTK_REAL, dimension(:,:,:), pointer :: phiR 
        CCTK_REAL, dimension(:,:,:), pointer :: phiI 
        CCTK_REAL, dimension(:,:,:), pointer :: piR 
        CCTK_REAL, dimension(:,:,:), pointer :: piI 

        CCTK_REAL, dimension(:,:,:), pointer :: x 
        CCTK_REAL, dimension(:,:,:), pointer :: y 
        CCTK_REAL, dimension(:,:,:), pointer :: z

        CCTK_INT  :: shp(3), nx, ny, nz, nd, warn, direc
        CCTK_INT  :: d_type, dd_type, rc
        CCTK_INT  :: psi4_analysis 
        CCTK_INT  :: constraints_analysis 
        CCTK_REAL :: dx, dy, dz
        CCTK_REAL :: G_scale_factor, B_scale_factor
        CCTK_REAL, allocatable, dimension(:) :: x1d, y1d, z1d
     
        CCTK_INT     :: i, j, k, la, lb, lc
        character*64 :: string

        logical, parameter        :: ltrace  = .false.
        logical, parameter        :: ltrace2 = .false.

        ! memory pointers
#include "AN_declare_ptrs.h"

        include 'mem.inc'
        integer  mem_alloc
        external mem_alloc

! Travis here, adding more stuff for the phi2 analysis

        CCTK_REAL, dimension(:,:,:), pointer :: alpha
        CCTK_REAL, dimension(:,:,:), pointer :: shift1,shift2,shift3
! temporary variables:
        CCTK_REAL :: v1,v2,v3,v4,v5,v6
        CCTK_INT  :: i1,i2,i3,i4,i5,i6
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
        CCTK_REAL :: x2,y2,z2, r, rcyl
        CCTK_REAL :: f_tr, f_rphi,f_tphi,f_thephi,vb1,vb2

         CCTK_REAL :: myl2norm3d
         external  :: myl2norm3d

!for extradiss out
	CCTK_REAL :: rcoor, absm1, absm2, bhsize

        if (ltrace) write(0,*) '*** BEGIN bssn_analysis'
        if (ltrace) write(0,*) '*** d_type=',par(P_DERIV_ORDER)

        ! allocate memory for first derivs and second derivs
        shp     = shape(u2(1)%d)
        nx      = shp(1)
        ny      = shp(2)
        nz      = shp(3)
        nd      = nx*ny*nz

        if (nx .ne. nint(par(P_NX))) then
          write(0,*)'*** analysis.f90:  Problem with inconsistent',
     &              ' array sizes'
          write(0,*) 'nx        = ',nx
          write(0,*) 'par(P_NX) = ',nint(par(P_NX))
        end if
        if (ny .ne. nint(par(P_NY))) then
          write(0,*)'*** analysis.f90:  Problem with inconsistent',
     &              ' array sizes'
          write(0,*) 'ny        = ',ny
          write(0,*) 'par(P_NY) = ',nint(par(P_NY))
        end if
        if (nz .ne. nint(par(P_NZ))) then
          write(0,*)'*** analysis.f90:  Problem with inconsistent',
     &              ' array sizes'
          write(0,*) 'nz        = ',nz
          write(0,*) 'par(P_NZ) = ',nint(par(P_NZ))
        end if

        allocate(x1d(shp(1)), y1d(shp(2)), z1d(shp(3)))

        ! --- The coords
        x      => w(G_XPHYS)%d
        y      => w(G_YPHYS)%d
        z      => w(G_ZPHYS)%d
        x1d(:) = x(:,1,1)
        y1d(:) = y(1,:,1)
        z1d(:) = z(1,1,:)


        dx      = par(P_DX)
        dy      = par(P_DY)
        dz      = par(P_DZ)
        d_type  = nint(par(P_DERIV_ORDER))
        dd_type = nint(par(P_DERIV_ORDER))
        warn    = 1
        G_scale_factor = par(P_G_SCALE_FACTOR)
        B_scale_factor = par(P_B_SCALE_FACTOR)

        psi4_analysis        = nint(par(P_PSI4_ANALYSIS))
        constraints_analysis = nint(par(P_CONSTRAINTS_ANALYSIS))

        chi         => u2(G_CHI)%d
        trK         => u2(G_TRK)%d
        gtxx        => u2(G_GT11)%d
        gtxy        => u2(G_GT12)%d
        gtxz        => u2(G_GT13)%d
        gtyy        => u2(G_GT22)%d
        gtyz        => u2(G_GT23)%d
        gtzz        => u2(G_GT33)%d
        Atxx        => u2(G_A11)%d
        Atxy        => u2(G_A12)%d
        Atxz        => u2(G_A13)%d
        Atyy        => u2(G_A22)%d
        Atyz        => u2(G_A23)%d
        Atzz        => u2(G_A33)%d
        Gamtx       => u2(G_GAM1)%d
        Gamty       => u2(G_GAM2)%d
        Gamtz       => u2(G_GAM3)%d
 
! get the lapse and shift also
        alpha       => u2(G_ALPHA)%d
        shift1      => u2(G_SHIFT1)%d
        shift2      => u2(G_SHIFT2)%d
        shift3      => u2(G_SHIFT3)%d

        Ex          => u2(G_EX)%d
        Ey          => u2(G_EY)%d
        Ez          => u2(G_EZ)%d
        Bx          => u2(G_BX)%d
        By          => u2(G_BY)%d
        Bz          => u2(G_BZ)%d
        Phi_em      => u2(G_PHI_EM)%d
        Psi_em      => u2(G_PSI_EM)%d

        phiR        => u2(G_PHIR)%d
        phiI        => u2(G_PHII)%d
        piR         => u2(G_PIR)%d
        piI         => u2(G_PII)%d

        call cal_adm_metric(
     &                    v(G_G11)%d,  v(G_G12)%d,  v(G_G13)%d,  
     &                    v(G_G22)%d,  v(G_G23)%d,  v(G_G33)%d,
     &                    gtxx, gtxy, gtxz, gtyy, gtyz, gtzz,
     &                    chi,  nx,  ny,  nz)


#include "AN_alloc.h"

#include "AN_derivs.h"

! calculate the Hamiltonian and momentum constraints 

!        if (do_cal_bssn_analysis .ge. 1) then

        if ( constraints_analysis == 1 ) then
!        call cal_bssn_analysis(u2, v, w,
        call calc_constraints(u2, v, w,
#include "AN_call_analysis.h"
     &             x1d, y1d, z1d, nx, ny, nz, par)
!         write(*,*) 'analysis: |ham|:',
!     &          myl2norm3d(w(H_HAM)%d, nx, ny, nz)
        end if


! calculate psi4 
        if ( psi4_analysis == 1 ) then          
          call calc_psi4(u2, v, w, 
#include "AN_call_analysis.h"
     &                   nx, ny, nz, par)
        end if



#include "AN_dealloc.h"

#if 0
! Adding phi2 stuff here

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx

! make the lowered indice 4-metric g4dd
! get the lowered spatial metric first:
          v1 = 1.d0 / chi(i,j,k)
          y_ij(1,1) = v1*gtxx(i,j,k)
          y_ij(1,2) = v1*gtxy(i,j,k)
          y_ij(1,3) = v1*gtxz(i,j,k)
          y_ij(2,2) = v1*gtyy(i,j,k)
          y_ij(2,3) = v1*gtyz(i,j,k)
          y_ij(3,3) = v1*gtzz(i,j,k)
          y_ij(2,1) = y_ij(1,2)
          y_ij(3,1) = y_ij(1,3)
          y_ij(3,2) = y_ij(2,3)

! and the determinant:

!!          deth = exp(12.d0*phi(i,j,k))   
          deth = v1**3   

! create the up-index spatial metric:

          h11 = y_ij(1,1)
          h12 = y_ij(1,2)
          h13 = y_ij(1,3)
          h22 = y_ij(2,2)
          h23 = y_ij(2,3)
          h33 = y_ij(3,3)

          huu11 = -h23**2 + h22*h33
          huu12 = h13*h23 - h12*h33
          huu13 = -(h13*h22) + h12*h23
          huu22 = -h13**2 + h11*h33
          huu23 = h12*h13 - h11*h23
          huu33 = -h12**2 + h11*h22

          yuu_ij(1,1) = huu11/deth
          yuu_ij(1,2) = huu12/deth
          yuu_ij(1,3) = huu13/deth
          yuu_ij(2,2) = huu22/deth
          yuu_ij(2,3) = huu23/deth
          yuu_ij(3,3) = huu33/deth
          yuu_ij(2,1) = yuu_ij(1,2)
          yuu_ij(3,1) = yuu_ij(1,3)
          yuu_ij(3,2) = yuu_ij(2,3)

! plug in the raised index shift:

          bui(1) = shift1(i,j,k)
          bui(2) = shift2(i,j,k)
          bui(3) = shift3(i,j,k)

! now get the lowered shift vectors:

          do i1 = 1,3
             bdi(i1) = 0.d0
             do i2 = 1,3
                bdi(i1) = bdi(i1) + y_ij(i1,i2)*bui(i2)
             enddo
          enddo

! now create the down-index 4-metric:

          v1 = -(alpha(i,j,k))**2
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

          v1 = 1.d0/(alpha(i,j,k))**2
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

! get the individual pieces of g:

          g00 = g(4,4)
          g01 = g(4,1)
          g02 = g(4,2)
          g03 = g(4,3)
          g10 = g(1,4)
          g11 = g(1,1)
          g12 = g(1,2)
          g13 = g(1,3)
          g20 = g(2,4)
          g21 = g(2,1)
          g22 = g(2,2)
          g23 = g(2,3)
          g30 = g(3,4)
          g31 = g(3,1)
          g32 = g(3,2)
          g33 = g(3,3)

! get the Faraday tensor:

          Bu0 = 0.d0
          Bu1 = Bx(i,j,k)
          Bu2 = By(i,j,k)
          Bu3 = Bz(i,j,k)

          Bd0 = g00*Bu0 + g01*Bu1 + g02*Bu2 + g03*Bu3
          Bd1 = g01*Bu0 + g11*Bu1 + g12*Bu2 + g13*Bu3
          Bd2 = g02*Bu0 + g12*Bu1 + g22*Bu2 + g23*Bu3
          Bd3 = g03*Bu0 + g13*Bu1 + g23*Bu2 + g33*Bu3

          Eu0 = 0.d0
          Eu1 = Ex(i,j,k)
          Eu2 = Ey(i,j,k)
          Eu3 = Ez(i,j,k)
 
          Ed1 =  g11*Eu1 + g12*Eu2 + g13*Eu3
          Ed2 =  g12*Eu1 + g22*Eu2 + g23*Eu3
          Ed3 =  g13*Eu1 + g23*Eu2 + g33*Eu3

          v(G_poyntingx_density)%d(i,j,k) = (Ed2*Bd3 - Ed3*Bd2)/dsqrt(deth)
          v(G_poyntingy_density)%d(i,j,k) = (Ed3*Bd1 - Ed1*Bd3)/dsqrt(deth)
          v(G_poyntingz_density)%d(i,j,k) = (Ed1*Bd2 - Ed2*Bd1)/dsqrt(deth)

 
          alp = alpha(i,j,k)
          ! deth = exp(12.d0*phi(i,j,k))  -- already calculated above 
           
          tu0 = 1.0/alp; tu1 = -bui(1)/alp; 
          tu2 = -bui(2)/alp; tu3 = -bui(3)/alp
          eps = 1.0d0/dsqrt(deth)

          x2 = x1d(i)
          y2 = y1d(j)
          z2 = z1d(k)
          r = sqrt(x2**2 + y2**2 + z2**2)

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
&                  g,F_uu,Tl,Tn,Tm_R,Tm_I,x2,y2,z2)

! ok, found the phis, load them into the v funcs 

          v(G_PHI2R)%d(i,j,k) = r * Phi2R
!          v(G_PHI2R)%d(i,j,k) = Tl(1)
          v(G_PHI2I)%d(i,j,k) = r * Phi2I 
!          v(G_PHI2I)%d(i,j,k) = Tn(1)
          v(G_PHI0R)%d(i,j,k) = r * Phi0R
!          v(G_PHI0R)%d(i,j,k) = Tm_R(1)
          v(G_PHI0I)%d(i,j,k) = r * Phi0I
!          v(G_PHI0I)%d(i,j,k) = Tm_I(1)

      do la=1,4
         do lb=1,4
            F_du(la,lb)=0.d0
            do lc=1,4
               F_du(la,lb)= F_du(la,lb) + g(la,lc)*F_uu(lc,lb)
            enddo
         enddo
      enddo

      do la=1,4
         do lb=1,4
            F_dd(la,lb)=0.d0
            do lc=1,4
               F_dd(la,lb)= F_dd(la,lb) + g(la,lc)*F_du(lb,lc)
            enddo
         enddo
      enddo

      ! using the jacobians omega1 = f_tr/f_rphi
      rcyl = sqrt(x2**2 + y2**2)
      f_tr = (x2/r)*F_dd(4,1) + (y2/r)*F_dd(4,2) + (z2/r)*F_dd(4,3)
      f_rphi = (x2*y2/r)*(-F_dd(1,1) + F_dd(2,2)) 
     &       + (rcyl*rcyl/r)*F_dd(1,2)
     &       - (y2*z2/r)*F_dd(3,1) + (x2*z2/r)*F_dd(3,2)
      if (abs(f_rphi) .gt. 1E-10) then
         vb1 = f_tr/f_rphi
      else
         vb1 = 0.0d0
      end if
      f_tphi = (x2*z2/rcyl)*F_dd(4,1) + (y2*z2/rcyl)*F_dd(4,2) 
     &       - rcyl*F_dd(4,3)
      f_thephi = (x2*y2*z2/rcyl)*(-F_dd(1,1) + F_dd(2,2)) 
     &         + z2*rcyl*F_dd(1,2)  + y2*rcyl*F_dd(3,1) 
     &         - x2*rcyl*F_dd(3,2)
      if (abs(f_thephi) .gt. 1E-10) then
         vb2 = f_tphi/f_thephi
      else
         vb2 = 0.0d0
      end if


      v(G_VBR_OMEGA)%d(i,j,k)       = vb1
      v(G_VBTH_OMEGA)%d(i,j,k)      = vb2

          if (ltrace2) then
            if (i .eq. 31 .and. j .eq. 0 .and. k .eq. 0) then
              print *,'phi2R = ', Phi2R
              print *,'phi2I = ', Phi2I
              print *,'g     = ',g
              print *,'F_uu  = ',F_uu
              print *,'Tl    = ',Tl
              print *,'Tn    = ',Tn
              print *,'Tm_R  = ',Tm_R
              print *,'x2    = ',x2
              print *,'y2    = ',y2
              print *,'z2    = ',z2
            end if
          end if

        enddo
        enddo
        enddo
#endif

!this more/less out
	absm1 = par(P_BH_1_MASS)
	absm2 = par(P_BH_2_MASS)
	bhsize = absm1+absm2

	do i=1,nx
	do j=1,ny
	do k=1,nz
  	  rcoor = sqrt(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2)

 	  if(rcoor.gt.40*bhsize .and. rcoor.lt.60*bhsize) then
             w(H_WDISS)%d(i,j,k)=1.+par(P_EXTRADISSOUT)* 
     &                          (rcoor-40.*bhsize)/(20.*bhsize)
	  else if(rcoor.ge.60.*bhsize) then
	     w(H_WDISS)%d(i,j,k) = 1.+ par(P_EXTRADISSOUT)
	  end if
	end do
	end do
	end do
        if (ltrace) write(0,*) '*** END bssn_analysis'

        return
      end subroutine

      !---------------------------------------------------------
      !
      ! Point wise analysis, there is a parameter in hyper called
      !     point_wise_analysis in order to call this routine
      !
      !---------------------------------------------------------
   
      subroutine analysis_pt(u_pt, u0_pt, urk_pt, dxu_pt, dyu_pt, 
     &                       dzu_pt, v_pt , dxv_pt, dyv_pt, dzv_pt, 
     &                       w_pt, par)
        implicit none
        CCTK_REAL, dimension(NU), intent(inout)  :: u_pt,u0_pt,dxu_pt
        CCTK_REAL, dimension(NU), intent(inout)  :: dyu_pt
        CCTK_REAL, dimension(NU), intent(inout)  :: dzu_pt
        CCTK_REAL, dimension(NU), intent(inout)  :: urk_pt
        CCTK_REAL, dimension(NV), intent(inout)  :: v_pt
        CCTK_REAL, dimension(NV), intent(inout)  :: dxv_pt
        CCTK_REAL, dimension(NV), intent(inout)  :: dyv_pt
        CCTK_REAL, dimension(NV), intent(inout)  :: dzv_pt
        CCTK_REAL, dimension(NW), intent(in)     :: w_pt
        CCTK_REAL, dimension(NPAR),intent(inout) :: par
       
        return
      end subroutine analysis_pt


      end module  
