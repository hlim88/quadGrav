#include "cctk.h"

!---------------------------------------------------------------------------
!
! $Id:
!
!---------------------------------------------------------------------------

!----------------------------------------------------------------------
!
!     ADM2BSSN
!
!----------------------------------------------------------------------
      subroutine adm2bssn(gtd11, gtd12, gtd13, gtd22, gtd23, gtd33, 
     &                    Atd11, Atd12, Atd13, Atd22, Atd23, Atd33, 
     &                    chi3D, trK3D, Gamt1, Gamt2, Gamt3,           
     &                    gd11,  gd12,  gd13,  gd22,  gd23,  gd33,  
     &                    Kd11,  Kd12,  Kd13,  Kd22,  Kd23,  Kd33,  
     &                    x1d, y1d, z1d, dx, dy, dz,
     &                    d_type, nx,  ny,  nz)
        implicit none

        CCTK_INT                       :: nx, ny, nz, d_type
        CCTK_REAL, dimension(nx,ny,nz) :: gtd11, gtd12, gtd13, gtd22, 
     &                                    gtd23, gtd33
        CCTK_REAL, dimension(nx,ny,nz) :: Atd11, Atd12, Atd13, Atd22, 
     &                                    Atd23, Atd33
        CCTK_REAL, dimension(nx,ny,nz) :: chi3D, trK3D, Gamt1, Gamt2, 
     &                                    Gamt3
        CCTK_REAL, dimension(nx,ny,nz) :: gd11, gd12, gd13, gd22, 
     &                                    gd23, gd33
        CCTK_REAL, dimension(nx,ny,nz) :: Kd11, Kd12, Kd13, Kd22, 
     &                                    Kd23, Kd33
        CCTK_REAL, dimension(nx)       :: x1d
        CCTK_REAL, dimension(ny)       :: y1d
        CCTK_REAL, dimension(nz)       :: z1d
        CCTK_REAL :: dx, dy, dz

        ! local vars
        CCTK_INT  :: i, j, k
        CCTK_INT  :: nd, pp
        CCTK_INT  :: dx_gtd11, dx_gtd12, dx_gtd13
        CCTK_INT  :: dx_gtd22, dx_gtd23, dx_gtd33
        CCTK_INT  :: dy_gtd11, dy_gtd12, dy_gtd13
        CCTK_INT  :: dy_gtd22, dy_gtd23, dy_gtd33
        CCTK_INT  :: dz_gtd11, dz_gtd12, dz_gtd13
        CCTK_INT  :: dz_gtd22, dz_gtd23, dz_gtd33
  
        CCTK_REAL :: chi, trK, detgd, idetgd, detgtd, idetgtd, inv_chi 
        CCTK_REAL, dimension(3,3)   :: gd, gtd, Ad, Atd, gu, gtu, Kd
        CCTK_REAL, dimension(3)     :: Gamt
        CCTK_REAL, dimension(3,3,3) :: Ct, Ctd, Dgtd

        CCTK_REAL, parameter :: two = 2.d0 
        CCTK_REAL, parameter :: half = 0.5d0 
        CCTK_REAL, parameter :: third = 3.3333333333333333d-1 

        logical, parameter :: ltrace   = .false.
        CCTK_INT :: mem_alloc
        include 'mem.inc'
      
        ! Declare Maple temp vars
#include "adm2bssn_maple_temp_vars.h"

        if (ltrace) write(0,*)'*** BEGIN adm2bssn'

        nd = nx*ny*nz
        dx = x1d(2) - x1d(1)
        dy = y1d(2) - y1d(1)
        dz = z1d(2) - z1d(1)
      
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          gd(1,1) = gd11(i,j,k)
          gd(1,2) = gd12(i,j,k)
          gd(1,3) = gd13(i,j,k)
          gd(2,2) = gd22(i,j,k)
          gd(2,3) = gd23(i,j,k)
          gd(3,3) = gd33(i,j,k)
          gd(2,1) = gd(1,2)
          gd(3,1) = gd(1,3)
          gd(3,2) = gd(2,3)
      
          Kd(1,1) = Kd11(i,j,k)
          Kd(1,2) = Kd12(i,j,k)
          Kd(1,3) = Kd13(i,j,k)
          Kd(2,2) = Kd22(i,j,k)
          Kd(2,3) = Kd23(i,j,k)
          Kd(3,3) = Kd33(i,j,k)
          Kd(2,1) = Kd(1,2)
          Kd(3,1) = Kd(1,3)
          Kd(3,2) = Kd(2,3)
     
#include "adm2bssn.h"

          chi3D(i,j,k) = chi
          trK3D(i,j,k) = trK

          gtd11(i,j,k) = gtd(1,1)
          gtd12(i,j,k) = gtd(1,2)
          gtd13(i,j,k) = gtd(1,3)
          gtd22(i,j,k) = gtd(2,2)
          gtd23(i,j,k) = gtd(2,3)
          gtd33(i,j,k) = gtd(3,3)
      
          Atd11(i,j,k) = Atd(1,1)
          Atd12(i,j,k) = Atd(1,2)
          Atd13(i,j,k) = Atd(1,3)
          Atd22(i,j,k) = Atd(2,2)
          Atd23(i,j,k) = Atd(2,3)
          Atd33(i,j,k) = Atd(3,3)

        end do
        end do
        end do

        if (ltrace) then
          write(0,*)'ADM2BSSN:  Allocating memory for derivs',nd
        end if

        dx_gtd11 = mem_alloc(nd)
        dx_gtd12 = mem_alloc(nd)
        dx_gtd13 = mem_alloc(nd)
        dx_gtd22 = mem_alloc(nd)
        dx_gtd23 = mem_alloc(nd)
        dx_gtd33 = mem_alloc(nd)

        dy_gtd11 = mem_alloc(nd)
        dy_gtd12 = mem_alloc(nd)
        dy_gtd13 = mem_alloc(nd)
        dy_gtd22 = mem_alloc(nd)
        dy_gtd23 = mem_alloc(nd)
        dy_gtd33 = mem_alloc(nd)

        dz_gtd11 = mem_alloc(nd)
        dz_gtd12 = mem_alloc(nd)
        dz_gtd13 = mem_alloc(nd)
        dz_gtd22 = mem_alloc(nd)
        dz_gtd23 = mem_alloc(nd)
        dz_gtd33 = mem_alloc(nd)

        call deriv_x(q(dx_gtd11), gtd11, dx, nx, ny, nz, d_type)
        call deriv_x(q(dx_gtd12), gtd12, dx, nx, ny, nz, d_type)
        call deriv_x(q(dx_gtd13), gtd13, dx, nx, ny, nz, d_type)
        call deriv_x(q(dx_gtd22), gtd22, dx, nx, ny, nz, d_type)
        call deriv_x(q(dx_gtd23), gtd23, dx, nx, ny, nz, d_type)
        call deriv_x(q(dx_gtd33), gtd33, dx, nx, ny, nz, d_type)

        call deriv_y(q(dy_gtd11), gtd11, dy, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtd12), gtd12, dy, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtd13), gtd13, dy, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtd22), gtd22, dy, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtd23), gtd23, dy, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtd33), gtd33, dy, nx, ny, nz, d_type)

        call deriv_z(q(dz_gtd11), gtd11, dz, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtd12), gtd12, dz, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtd13), gtd13, dz, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtd22), gtd22, dz, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtd23), gtd23, dz, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtd33), gtd33, dz, nx, ny, nz, d_type)

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          gtd(1,1) = gtd11(i,j,k)
          gtd(1,2) = gtd12(i,j,k)
          gtd(1,3) = gtd13(i,j,k)
          gtd(2,2) = gtd22(i,j,k)
          gtd(2,3) = gtd23(i,j,k)
          gtd(3,3) = gtd33(i,j,k)
          gtd(2,1) = gtd(1,2)
          gtd(3,1) = gtd(1,3)
          gtd(3,2) = gtd(2,3)
      
          pp = (i-1) + nx*(j-1) + nx*ny*(k-1)
      
          Dgtd(1,1,1) = q(dx_gtd11 + pp)
          Dgtd(1,1,2) = q(dx_gtd12 + pp)
          Dgtd(1,1,3) = q(dx_gtd13 + pp)
          Dgtd(1,2,1) = Dgtd(1,1,2)
          Dgtd(1,2,2) = q(dx_gtd22 + pp)
          Dgtd(1,2,3) = q(dx_gtd23 + pp)
          Dgtd(1,3,1) = Dgtd(1,1,3)
          Dgtd(1,3,2) = Dgtd(1,2,3)
          Dgtd(1,3,3) = q(dx_gtd33 + pp)

          Dgtd(2,1,1) = q(dy_gtd11 + pp)
          Dgtd(2,1,2) = q(dy_gtd12 + pp)
          Dgtd(2,1,3) = q(dy_gtd13 + pp)
          Dgtd(2,2,1) = Dgtd(2,1,2)
          Dgtd(2,2,2) = q(dy_gtd22 + pp)
          Dgtd(2,2,3) = q(dy_gtd23 + pp)
          Dgtd(2,3,1) = Dgtd(2,1,3)
          Dgtd(2,3,2) = Dgtd(2,2,3)
          Dgtd(2,3,3) = q(dy_gtd33 + pp)
      
          Dgtd(3,1,1) = q(dz_gtd11 + pp)
          Dgtd(3,1,2) = q(dz_gtd12 + pp)
          Dgtd(3,1,3) = q(dz_gtd13 + pp)
          Dgtd(3,2,1) = Dgtd(3,1,2)
          Dgtd(3,2,2) = q(dz_gtd22 + pp)
          Dgtd(3,2,3) = q(dz_gtd23 + pp)
          Dgtd(3,3,1) = Dgtd(3,1,3)
          Dgtd(3,3,2) = Dgtd(3,2,3)
          Dgtd(3,3,3) = q(dz_gtd33 + pp)

#include "adm2bssn2.h"

          Gamt1(i,j,k) = Gamt(1)
          Gamt2(i,j,k) = Gamt(2)
          Gamt3(i,j,k) = Gamt(3)

        end do
        end do
        end do

        call mem_dealloc(dx_gtd11, nd)
        call mem_dealloc(dx_gtd12, nd)
        call mem_dealloc(dx_gtd13, nd)
        call mem_dealloc(dx_gtd22, nd)
        call mem_dealloc(dx_gtd23, nd)
        call mem_dealloc(dx_gtd33, nd)

        call mem_dealloc(dy_gtd11, nd)
        call mem_dealloc(dy_gtd12, nd)
        call mem_dealloc(dy_gtd13, nd)
        call mem_dealloc(dy_gtd22, nd)
        call mem_dealloc(dy_gtd23, nd)
        call mem_dealloc(dy_gtd33, nd)

        call mem_dealloc(dz_gtd11, nd)
        call mem_dealloc(dz_gtd12, nd)
        call mem_dealloc(dz_gtd13, nd)
        call mem_dealloc(dz_gtd22, nd)
        call mem_dealloc(dz_gtd23, nd)
        call mem_dealloc(dz_gtd33, nd)

        if (ltrace) write(0,*)'*** END adm2bssn'
 
        return
      end subroutine

!----------------------------------------------------------------------
!
!     ADM2BSSN
!
!----------------------------------------------------------------------
      subroutine cal_adm_metric(
     &                    gd11,  gd12,  gd13,  gd22,  gd23,  gd33,  
     &                    gtd11, gtd12, gtd13, gtd22, gtd23, gtd33, 
     &                    chi3D,
     &                    nx,  ny,  nz)
        implicit none
        CCTK_INT                       :: nx, ny, nz
        CCTK_REAL, dimension(nx,ny,nz) :: gtd11, gtd12, gtd13, gtd22, 
     &                                    gtd23, gtd33
        CCTK_REAL, dimension(nx,ny,nz) :: gd11, gd12, gd13, gd22, 
     &                                    gd23, gd33
        CCTK_REAL, dimension(nx,ny,nz) :: chi3D

        ! local vars
        CCTK_INT  :: i, j, k
        CCTK_REAL              :: inv_chi 

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          inv_chi = 1.d0 / chi3D(i,j,k)
          gd11(i,j,k) = inv_chi * gtd11(i,j,k)
          gd12(i,j,k) = inv_chi * gtd12(i,j,k)
          gd13(i,j,k) = inv_chi * gtd13(i,j,k)
          gd22(i,j,k) = inv_chi * gtd22(i,j,k)
          gd23(i,j,k) = inv_chi * gtd23(i,j,k)
          gd33(i,j,k) = inv_chi * gtd33(i,j,k)
        end do
        end do
        end do
       
        return
      end subroutine

!----------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------
      subroutine cal_adm_exccurv(
     &                    Kd11,   Kd12,   Kd13,   Kd22,   Kd23,   Kd33,
     &                    A11,    A12,    A13,    A22,    A23,    A33,
     &                    tr_K,
     &                    gtd11,  gtd12,  gtd13,  gtd22,  gtd23,  gtd33,
     &                    chi3d,  x1d,    y1d,    z1d,
     &                    nx,  ny,  nz)
        implicit none
        CCTK_INT                       :: nx, ny, nz
        CCTK_REAL, dimension(nx,ny,nz) :: Kd11, Kd12, Kd13, Kd22,
     &                                    Kd23, Kd33
        CCTK_REAL, dimension(nx,ny,nz) :: A11, A12, A13, A22,
     &                                    A23, A33, tr_K
        CCTK_REAL, dimension(nx,ny,nz) :: gtd11, gtd12, gtd13, gtd22,
     &                                    gtd23, gtd33
        CCTK_REAL, dimension(nx,ny,nz) :: chi3d
        CCTK_REAL, dimension(nx)       :: x1d
        CCTK_REAL, dimension(ny)       :: y1d
        CCTK_REAL, dimension(nz)       :: z1d

        ! local vars
        CCTK_INT  :: i, j, k
        CCTK_REAL :: third_TrK, inv_chi  
        CCTK_REAL, parameter :: third = 0.33333333333333333d0

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          inv_chi = 1.d0 / chi3d(i,j,k)
          third_TrK   = third * tr_K(i,j,k)

          Kd11(i,j,k) = inv_chi * (A11(i,j,k) + third_TrK*gtd11(i,j,k))
          Kd12(i,j,k) = inv_chi * (A12(i,j,k) + third_TrK*gtd12(i,j,k))
          Kd13(i,j,k) = inv_chi * (A13(i,j,k) + third_TrK*gtd13(i,j,k))
          Kd22(i,j,k) = inv_chi * (A22(i,j,k) + third_TrK*gtd22(i,j,k))
          Kd23(i,j,k) = inv_chi * (A23(i,j,k) + third_TrK*gtd23(i,j,k))
          Kd33(i,j,k) = inv_chi * (A33(i,j,k) + third_TrK*gtd33(i,j,k))

        end do
        end do
        end do

        return
      end subroutine

!----------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------

