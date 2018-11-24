#include "cctk.h"
      !----------------------------------------------------------------
      !
      !
      !
      !----------------------------------------------------------------
      MODULE m_enforce_bssn
      use GF
      use params
      implicit none

      CONTAINS

      !----------------------------------------------------------------
      !
      !
      !
      !----------------------------------------------------------------
      subroutine enforce_bssn_constraints (u, w, par)
        implicit   none

        type(gridfunction), dimension(NU_G)     :: u
        type(gridfunction), dimension(NW)       :: w
        CCTK_REAL, dimension(NPAR), intent(in)  :: par
        CCTK_INT                                :: nx, ny, nz

        ! local variables
        CCTK_REAL, dimension(:,:,:), pointer :: gtxx, gtxy, gtxz
        CCTK_REAL, dimension(:,:,:), pointer ::       gtyy, gtyz
        CCTK_REAL, dimension(:,:,:), pointer ::             gtzz
        CCTK_REAL, dimension(:,:,:), pointer :: Atxx, Atxy, Atxz
        CCTK_REAL, dimension(:,:,:), pointer ::       Atyy, Atyz
        CCTK_REAL, dimension(:,:,:), pointer ::             Atzz

        CCTK_REAL, dimension(:,:,:), pointer ::          Alpha3d
        CCTK_REAL, dimension(:,:,:), pointer ::            chi3d

        CCTK_REAL, dimension(:,:,:), pointer ::            tr_A 
        CCTK_REAL, dimension(:,:,:), pointer ::         detgt_m1 

        CCTK_INT    ::  i, j, k, shp(3)
        CCTK_REAL   ::  trace_Atd, det_gtd_to_third
        CCTK_REAL   ::  det_gtd, det_gtd_to_neg_third 
        CCTK_REAL   ::  neg_one_third_trace_Atd 
        CCTK_REAL, dimension(3,3)   ::  gtd, gtu, Atd
        CCTK_REAL   :: idet_gtd
        CCTK_REAL   :: chi_floor 
        CCTK_REAL, parameter ::  one_third = 0.333333333333333333d0 

        shp = shape(u(1)%d)
        nx = shp(1)
        ny = shp(2)
        nz = shp(3)

        if (nx .ne. nint(par(P_NX))) then
          write(0,*)'*** enforce_bssn.F90:  Problem with inconsistent',
     &              ' array sizes'
          write(0,*) 'nx        = ',nx
          write(0,*) 'par(P_NX) = ',nint(par(P_NX))
        end if
        if (ny .ne. nint(par(P_NY))) then
          write(0,*)'*** enforce_bssn.F90:  Problem with inconsistent',
     &              ' array sizes'
          write(0,*) 'ny        = ',ny
          write(0,*) 'par(P_NY) = ',nint(par(P_NY))
        end if
        if (nz .ne. nint(par(P_NZ))) then
          write(0,*)'*** enforce_bssn.F90:  Problem with inconsistent',
     &              ' array sizes'
          write(0,*) 'nz        = ',nz
          write(0,*) 'par(P_NZ) = ',nint(par(P_NZ))
        end if

        chi_floor = par(P_BSSN_CHI_FLOOR) 

        Alpha3d     => u(G_ALPHA)%d
        chi3d       => u(G_CHI)%d

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

        tr_A        => w(G_TR_A)%d
        detgt_m1    => w(G_DETGT_M1)%d

        do k = 1, nz 
        do j = 1, ny 
        do i = 1, nx 

!          ! apply floor to the lapse and chi
!          !Alpha3d(i,j,k) = max(Alpha3d(i,j,k), 1.0d-8)
!          chi3d(i,j,k)   = max(chi3d(i,j,k), 1.0d-4)
          


          ! Require gtd to have unit determinant
          gtd(1,1) = gtxx(i,j,k) 
          gtd(1,2) = gtxy(i,j,k) 
          gtd(1,3) = gtxz(i,j,k) 
          gtd(2,1) = gtd(1,2)
          gtd(2,2) = gtyy(i,j,k) 
          gtd(2,3) = gtyz(i,j,k) 
          gtd(3,1) = gtd(1,3)
          gtd(3,2) = gtd(2,3)
          gtd(3,3) = gtzz(i,j,k) 
     
          det_gtd =   gtd(1,1)*( gtd(2,2)*gtd(3,3) - gtd(2,3)**2) 
     &              - gtd(1,2)**2*gtd(3,3)
     &              + 2*gtd(1,2)*gtd(1,3)*gtd(2,3) 
     &              - gtd(1,3)**2*gtd(2,2) 

          det_gtd_to_third = det_gtd**(one_third)  
          det_gtd_to_neg_third = 1.d0 / det_gtd_to_third

          gtd(1,1) = det_gtd_to_neg_third * gtd(1,1)  
          gtd(1,2) = det_gtd_to_neg_third * gtd(1,2)  
          gtd(1,3) = det_gtd_to_neg_third * gtd(1,3)  
          gtd(2,2) = det_gtd_to_neg_third * gtd(2,2)  
          gtd(2,3) = det_gtd_to_neg_third * gtd(2,3)  
          gtd(3,3) = det_gtd_to_neg_third * gtd(3,3)  

          det_gtd = gtd(1,1)*gtd(2,2)*gtd(3,3)-gtd(1,1)*gtd(2,3)**2
     &            - gtd(1,2)**2*gtd(3,3)+2*gtd(1,2)*gtd(1,3)*gtd(2,3)
     &            - gtd(1,3)**2*gtd(2,2)

!          detgt(i,j,k) = det_gtd
          detgt_m1(i,j,k) = det_gtd - 1.d0 

          if (abs(det_gtd-1.0d0) .gt. 1.0d-6) then
            write(0,*)'enforce_bssn_constraint: det(gtd) != 1.',
     &                 ' det=',det_gtd
            write(0,*)'   gtd(1,1) ',gtd(1,1)
            write(0,*)'   gtd(1,2) ',gtd(1,2)
            write(0,*)'   gtd(1,3) ',gtd(1,3)
            write(0,*)'   gtd(2,2) ',gtd(2,2)
            write(0,*)'   gtd(2,3) ',gtd(2,3)
            write(0,*)'   gtd(3,3) ',gtd(3,3)
          end if

          idet_gtd = 0.1D1/det_gtd
          gtu(1,1) = idet_gtd*(gtd(2,2)*gtd(3,3)-gtd(2,3)**2)
          gtu(1,2) = idet_gtd*(-gtd(1,2)*gtd(3,3)+gtd(1,3)*gtd(2,3))
          gtu(1,3) = idet_gtd*(gtd(1,2)*gtd(2,3)-gtd(1,3)*gtd(2,2))
          gtu(2,1) = gtu(1,2)
          gtu(2,2) = idet_gtd*(gtd(1,1)*gtd(3,3)-gtd(1,3)**2)
          gtu(2,3) = idet_gtd*(-gtd(1,1)*gtd(2,3)+gtd(1,2)*gtd(1,3))
          gtu(3,1) = gtu(1,3)
          gtu(3,2) = gtu(2,3)
          gtu(3,3) = idet_gtd*(gtd(1,1)*gtd(2,2)-gtd(1,2)**2)

!         Require Atd to be traceless.
          Atd(1,1) = Atxx(i,j,k) 
          Atd(1,2) = Atxy(i,j,k) 
          Atd(1,3) = Atxz(i,j,k) 
          Atd(2,1) = Atd(1,2)
          Atd(2,2) = Atyy(i,j,k) 
          Atd(2,3) = Atyz(i,j,k) 
          Atd(3,1) = Atd(1,3)
          Atd(3,2) = Atd(2,3)
          Atd(3,3) = Atzz(i,j,k) 

          trace_Atd =   Atd(1,1)*gtu(1,1) 
     &                + Atd(2,2)*gtu(2,2)   
     &                + Atd(3,3)*gtu(3,3)   
     &                + 2.d0 * (   Atd(1,2)*gtu(1,2) 
     &                           + Atd(1,3)*gtu(1,3)   
     &                           + Atd(2,3)*gtu(2,3)  )    

          neg_one_third_trace_Atd = - one_third * trace_Atd 

          Atd(1,1) = Atd(1,1) + neg_one_third_trace_Atd * gtd(1,1)   
          Atd(1,2) = Atd(1,2) + neg_one_third_trace_Atd * gtd(1,2)   
          Atd(1,3) = Atd(1,3) + neg_one_third_trace_Atd * gtd(1,3)   
          Atd(2,2) = Atd(2,2) + neg_one_third_trace_Atd * gtd(2,2)   
          Atd(2,3) = Atd(2,3) + neg_one_third_trace_Atd * gtd(2,3)   
          Atd(3,3) = Atd(3,3) + neg_one_third_trace_Atd * gtd(3,3)   

          tr_A(i,j,k)  =   Atd(1,1)*gtu(1,1) 
     &                 + Atd(2,2)*gtu(2,2)   
     &                 + Atd(3,3)*gtu(3,3)   
     &                 + 2.d0 * (   Atd(1,2)*gtu(1,2) 
     &                            + Atd(1,3)*gtu(1,3)   
     &                            + Atd(2,3)*gtu(2,3)  )    

          Atxx(i,j,k) = Atd(1,1) 
          Atxy(i,j,k) = Atd(1,2) 
          Atxz(i,j,k) = Atd(1,3) 
          Atyy(i,j,k) = Atd(2,2) 
          Atyz(i,j,k) = Atd(2,3) 
          Atzz(i,j,k) = Atd(3,3) 

          gtxx(i,j,k) = gtd(1,1) 
          gtxy(i,j,k) = gtd(1,2) 
          gtxz(i,j,k) = gtd(1,3) 
          gtyy(i,j,k) = gtd(2,2) 
          gtyz(i,j,k) = gtd(2,3) 
          gtzz(i,j,k) = gtd(3,3) 


          if ( chi3D(i,j,k) .lt. chi_floor ) then
            chi3D(i,j,k) = chi_floor
          endif


        end do 
        end do 
        end do 

        return 
      end subroutine 


      END MODULE 
