!-----------------------------------------------------------------------
!
!    $Id$
!
!-----------------------------------------------------------------------

#include "cctk.h"

MODULE HYPER_DERIVS_21
use params

CONTAINS

!-----------------------------------------------------------------------
!
!    SUBROUTINE D21NoMask:
!
!    implement Strand derivative.  2nd order interior, 1st order boundary
!
!-----------------------------------------------------------------------
  subroutine D21NoMask(Du, u, direction, par)
    implicit none
    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT :: direction


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, j, k,nx, ny, nz, shp(3)
    

    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)
    
    if (direction.eq.1) then 
      Du(1,:,:) = (u(2,:,:)-u(1,:,:))/dx
      do i = 2, nx-1
        Du(i,:,:) = (u(i+1,:,:)-u(i-1,:,:))/(2*dx)
      end do
      Du(nx,:,:) = (u(nx,:,:)-u(nx-1,:,:))/dx
    else if (direction.eq.2) then 
      Du(:,1,:) = (u(:,2,:)-u(:,1,:))/dy
      do j = 2, ny-1
        Du(:,j,:) = (u(:,j+1,:)-u(:,j-1,:))/(2*dy)
      end do
      Du(:,ny,:) = (u(:,ny,:)-u(:,ny-1,:))/dy
    else if (direction==3) then 
      Du(:,:,1) = (u(:,:,2)-u(:,:,1))/dz
      do k = 2, nz-1
        Du(:,:,k) = (u(:,:,k+1)-u(:,:,k-1))/(2*dz)
      end do
      Du(:,:,nz) = (u(:,:,nz)-u(:,:,nz-1))/dz
    else 
      write(*,*) "D21NoMask: Wrong direction in the derivative" 
      STOP
    end if 

    return
  end subroutine D21NoMask


!-----------------------------------------------------------------------
!
!    D21
!
!  At times a derivative operator may be legitimately undefined.
!  This can happen in the ghost zones, if they are very close to
!  a mask boundary.  At such times the derivative is simply set
!  to zero.  For debugging purposes, the ltrace2 parameter can be
!  used to control whether an message is printed for undefined 
!  derivatives.
!
!-----------------------------------------------------------------------
  subroutine D21(Du, u, direction, x, y, z, mask, imask, par)
    implicit none

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: x, y, z
    CCTK_INT :: direction
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: mask
    CCTK_INT, dimension(:,:,:), INTENT(IN)   :: imask


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, j, k,nx, ny, nz, shp(3)
    CCTK_INT                                 :: d
    CCTK_INT  :: gz_xmin, gz_xmax, gz_ymin, gz_ymax, gz_zmin, gz_zmax
    CCTK_INT  :: bbox1, bbox2, bbox3, bbox4, bbox5, bbox6
    logical, parameter :: debug_deriv = .false.

    !  If ltrace2 is TRUE, then a warning message is printed for 
    !  undefined derivatives are flagged.
    !  If ltrace2 is FALSE, then no warning message is printed.
    logical, parameter :: ltrace2 = .false.


    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    !
    ! Check for ghostzones
    !
    gz_xmin = 1
    gz_xmax = nx 
    gz_ymin = 1
    gz_ymax = ny 
    gz_zmin = 1
    gz_zmax = nz 

    bbox1 = nint(par(P_BBOX1))
    bbox2 = nint(par(P_BBOX2))
    bbox3 = nint(par(P_BBOX3))
    bbox4 = nint(par(P_BBOX4))
    bbox5 = nint(par(P_BBOX5))
    bbox6 = nint(par(P_BBOX6))

    if (bbox1 .eq. 0) then
      gz_xmin = gz_xmin + nint(par(P_NGHOSTZONES_X))
    end if
    if (bbox2 .eq. 0) then
      gz_xmax = gz_xmax - nint(par(P_NGHOSTZONES_X))
    end if

    if (bbox3 .eq. 0) then
      gz_ymin = gz_ymin + nint(par(P_NGHOSTZONES_Y))
    end if
    if (bbox4 .eq. 0) then
      gz_ymax = gz_ymax - nint(par(P_NGHOSTZONES_Y))
    end if

    if (bbox5 .eq. 0) then
      gz_zmin = gz_zmin + nint(par(P_NGHOSTZONES_Z))
    end if
    if (bbox6 .eq. 0) then
      gz_zmax = gz_zmax - nint(par(P_NGHOSTZONES_Z))
    end if

    if (direction.eq.1) then
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = imask(i,j,k) / 10000
        if (d .eq. P_STENCIL_CENTER) then
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = (u(i+1,j,k) - u(i-1,j,k)) / (2.0*dx)
        else if (d .eq. P_STENCIL_CENTER_LEFT .or. &
                 d .eq. P_STENCIL_CENTER_RIGHT) then
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = (u(i+1,j,k) - u(i-1,j,k)) / (2.0*dx)
        else if (d .eq. P_STENCIL_LEFT) then
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = (u(i+1,j,k) - u(i,j,k)) / dx
        else if (d .eq. P_STENCIL_RIGHT) then
          if (debug_deriv) then
            if (nint(mask(i,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = (u(i,j,k) - u(i-1,j,k)) / dx
        else if (d .eq. P_STENCIL_MASK) then
          if (debug_deriv) then
            if (nint(mask(i,j,k)) .gt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = (2.0*u(i+1,j,k) - u(i,j,k) - u(i-1,j,k))/(3.0*dx)
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = (-2.0*u(i-1,j,k) + u(i,j,k) + u(i+1,j,k))/(3.0*dx)
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = 4.0/(7.0*dx)&
                      *(u(i+1,j,k) - 0.25*u(i,j,k) - 0.75*u(i-1,j,k))
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = 4.0/(7.0*dx)&
                      *(-u(i-1,j,k) + 0.25*u(i,j,k) + 0.75*u(i+1,j,k))
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
        else if (d .eq. P_STENCIL_UNDEF) then
          if (ltrace2) then
            if ((bbox1==0 .and. i < gz_xmin) .or.  &
                (bbox2==0 .and. i > gz_xmax)) then
              write(0,*)'d21: OK---Undefined x derivative type in ghostzone'
            else
              write(0,*)'d21: ### Undefined x derivative outside of ghost zone!'
              write(0,*)'d21: imask                          ',imask(i,j,k)
              write(0,*)'d21: i, nx                          ',i, nx
              write(0,*)'d21: gz_xmin                        ',gz_xmin
              write(0,*)'d21: gz_xmax                        ',gz_xmax
              write(0,*)'d21: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_X)
              write(0,*)'d21: par(bbox1), par(bbox2)         ',&
                                                par(P_BBOX1), par(P_BBOX2)
              write(0,*)'d21: bbox1, bbox2                   ',&
                                                bbox1,bbox2
              stop
            end if
          end if
          Du(i,j,k) = 0.0
        else
          write(0,*)'d21: Illegal x derivative type ',d
          write(0,*)'d21: imask                     ',imask(i,j,k)
          write(0,*)'     at location ',x(i,j,k),y(i,j,k),z(i,j,k)
          Du(i,j,k) = 0.0
	  stop
        end if
      end do
      end do
      end do

    else if (direction.eq.2) then

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod((imask(i,j,k)/100), 100)
        if (d .eq. P_STENCIL_CENTER) then
          Du(i,j,k) = (u(i,j+1,k) - u(i,j-1,k)) / (2.0*dy)
        else if (d .eq. P_STENCIL_LEFT) then
          Du(i,j,k) = (u(i,j+1,k) - u(i,j,k)) / dy
        else if (d .eq. P_STENCIL_CENTER_LEFT .or. &
                 d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = (u(i,j+1,k) - u(i,j-1,k)) / (2.0*dy)
        else if (d .eq. P_STENCIL_RIGHT) then
          Du(i,j,k) = (u(i,j,k) - u(i,j-1,k)) / dy
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = (2.0*u(i,j+1,k) - u(i,j,k) - u(i,j-1,k))/(3.0*dy)
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = (-2.0*u(i,j-1,k) + u(i,j,k) + u(i,j+1,k))/(3.0*dy)
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = 4.0/(7.0*dy)&
                      *(u(i,j+1,k) - 0.25*u(i,j,k) - 0.75*u(i,j-1,k))
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = 4.0/(7.0*dy)&
                      *(-u(i,j-1,k) + 0.25*u(i,j,k) + 0.75*u(i,j+1,k))
        else if (d .eq. P_STENCIL_UNDEF) then
          if (ltrace2) then
            if ((bbox3==0 .and. j < gz_ymin) .or. &
                (bbox4==0 .and. j > gz_ymax)) then
              write(0,*)'d21: OK---Undefined y derivative type in ghostzone'
            else
              write(0,*)'d21: ### Undefined y derivative outside of ghost zone!'
              write(0,*)'d21: imask                          ',imask(i,j,k)
              write(0,*)'d21: j, ny                          ',j, ny
              write(0,*)'d21: gz_ymin                        ',gz_ymin
              write(0,*)'d21: gz_ymax                        ',gz_ymax
              write(0,*)'d21: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_Y)
              write(0,*)'d21: par(bbox3), par(bbox4)         ',&
                                                par(P_BBOX3), par(P_BBOX4)
              write(0,*)'d21: bbox3, bbox4                   ',&
                                                bbox3,bbox4
              stop
            end if
          end if
          Du(i,j,k) = 0.0
        else
          write(0,*)'d21: Illegal y derivative type   ',d
          write(0,*)'d21: imask                       ',imask(i,j,k)
          write(0,*)'     at location ',x(i,j,k),y(i,j,k),z(i,j,k)
          Du(i,j,k) = 0.0
	  stop
        end if
      end do
      end do
      end do

    else if (direction.eq.3) then
 
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod(imask(i,j,k),100)
        if (d .eq. P_STENCIL_CENTER) then
          Du(i,j,k) = (u(i,j,k+1) - u(i,j,k-1)) / (2.0*dz)
        else if (d .eq. P_STENCIL_LEFT) then
          Du(i,j,k) = (u(i,j,k+1) - u(i,j,k)) / dz
        else if (d .eq. P_STENCIL_CENTER_LEFT .or. &
                 d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = (u(i,j,k+1) - u(i,j,k-1)) / (2.0*dz)
        else if (d .eq. P_STENCIL_RIGHT) then
          Du(i,j,k) = (u(i,j,k) - u(i,j,k-1)) / dz
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = (2.0*u(i,j,k+1) - u(i,j,k) - u(i,j,k-1))/(3.0*dz)
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = (-2.0*u(i,j,k-1) + u(i,j,k) + u(i,j,k+1))/(3.0*dz)
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = 4.0/(7.0*dz)&
                      *(u(i,j,k+1) - 0.25*u(i,j,k) - 0.75*u(i,j,k-1))
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = 4.0/(7.0*dz)&
                      *(-u(i,j,k-1) + 0.25*u(i,j,k) + 0.75*u(i,j,k+1))
        else if (d .eq. P_STENCIL_UNDEF) then
          if (ltrace2) then
            if ((bbox5==0 .and. k < gz_zmin) .or. &
                (bbox6==0 .and. k > gz_zmax)) then
              write(0,*)'d21: OK---Undefined z derivative type in ghostzone'
            else
              write(0,*)'d21: ### Undefined z derivative outside of ghost zone!'
              write(0,*)'d21: imask                          ',imask(i,j,k)
              write(0,*)'d21: k, nz                          ',k, nz
              write(0,*)'d21: gz_zmin                        ',gz_zmin
              write(0,*)'d21: gz_zmax                        ',gz_zmax
              write(0,*)'d21: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_Z)
              write(0,*)'d21: par(bbox5), par(bbox6)         ',&
                                                par(P_BBOX5), par(P_BBOX6)
              write(0,*)'d21: bbox5, bbox6                   ',&
                                                bbox5,bbox6
              !stop
            end if
          end if
          Du(i,j,k) = 0.0
        else
          write(0,*)'d21: Illegal z derivative type   ',d
          write(0,*)'d21: imask                       ',imask(i,j,k)
          write(0,*)'     at location ',x(i,j,k),y(i,j,k),z(i,j,k)
          Du(i,j,k) = 0.0
	  stop
        end if
      end do
      end do
      end do

    else
      write(0,*)'d21:  unknown direction ',direction
    end if

    return
  end subroutine D21

!-----------------------------------------------------------------------
!
! Working on this subroutine.  Not ready yet.
! SUBROUTINE DD2
!
!   calculates the second order derivative.
!
!-----------------------------------------------------------------------
  subroutine DD2(Du, u, direction, x, y, z, mask, imask, par)
    implicit none

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: x, y, z
    CCTK_INT :: direction
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: mask
    CCTK_INT, dimension(:,:,:), INTENT(IN)   :: imask

    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, j, k,nx, ny, nz, shp(3)
    CCTK_INT  :: d
!   CCTK_INT  :: gz_xmin, gz_xmax, gz_ymin, gz_ymax, gz_zmin, gz_zmax
!   CCTK_INT  :: bbox1, bbox2, bbox3, bbox4, bbox5, bbox6
    CCTK_INT  :: gz_xmin, gz_xmax, gz_ymin, gz_ymax, gz_zmin, gz_zmax
    CCTK_INT  :: bbox1, bbox2
    logical, parameter :: debug_deriv = .false.
    !  If ltrace2 is TRUE, then a warning message is printed for
    !  undefined derivatives are flagged.
    !  If ltrace2 is FALSE, then no warning message is printed.
    logical, parameter :: ltrace2 = .false.


    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    if (direction.eq.1) then
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = imask(i,j,k) / 10000
        if (d .eq. P_STENCIL_CENTER) then
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = (u(i+1,j,k) -2.0*u(i,j,k) +  u(i-1,j,k)) / dx**2
        else if (d .eq. P_STENCIL_CENTER_LEFT .or. &
                 d .eq. P_STENCIL_CENTER_RIGHT) then
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = (u(i+1,j,k) -2.0*u(i,j,k) +  u(i-1,j,k)) / dx**2
        else if (d .eq. P_STENCIL_LEFT) then
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = (u(i+2,j,k) -2.0*u(i+1,j,k) +  u(i,j,k)) / dx**2
        else if (d .eq. P_STENCIL_RIGHT) then
          if (debug_deriv) then
            if (nint(mask(i,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = (u(i-2,j,k) -2.0*u(i-1,j,k) +  u(i,j,k)) / dx**2
        else if (d .eq. P_STENCIL_MASK) then
          if (debug_deriv) then
            if (nint(mask(i,j,k)) .gt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = (2.0*u(i+1,j,k) - u(i,j,k) - u(i-1,j,k))/(3.0*dx)
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = (-2.0*u(i-1,j,k) + u(i,j,k) + u(i+1,j,k))/(3.0*dx)
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = 4.0/(7.0*dx)&
                      *(u(i+1,j,k) - 0.25*u(i,j,k) - 0.75*u(i-1,j,k))
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = 4.0/(7.0*dx)&
                      *(-u(i-1,j,k) + 0.25*u(i,j,k) + 0.75*u(i+1,j,k))
          if (debug_deriv) then
            if (nint(mask(i+1,j,k)) .lt. 0 .or. nint(mask(i-1,j,k)) .lt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i+1,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              stop
            end if
          end if
        else if (d .eq. P_STENCIL_UNDEF) then
          if (ltrace2) then
            if ((bbox1==0 .and. i < gz_xmin) .or.  &
                (bbox2==0 .and. i > gz_xmax)) then
              write(0,*)'d21: OK---Undefined x derivative type in ghostzone'
            else
              write(0,*)'d21: ### Undefined x derivative outside of ghost zone!'
              write(0,*)'d21: imask                          ',imask(i,j,k)
              write(0,*)'d21: i, nx                          ',i, nx
              write(0,*)'d21: gz_xmin                        ',gz_xmin
              write(0,*)'d21: gz_xmax                        ',gz_xmax
              write(0,*)'d21: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_X)
              write(0,*)'d21: par(bbox1), par(bbox2)         ',&
                                                par(P_BBOX1), par(P_BBOX2)
              write(0,*)'d21: bbox1, bbox2                   ',&
                                                bbox1,bbox2
              stop
            end if
          end if
          Du(i,j,k) = 0.0
        else
          write(0,*)'d21: Illegal x derivative type ',d
          write(0,*)'d21: imask                     ',imask(i,j,k)
          write(0,*)'     at location ',x(i,j,k),y(i,j,k),z(i,j,k)
          Du(i,j,k) = 0.0
          stop
        end if
      end do
      end do
      end do

!   else if (direction.eq.2) then
    end if




    return
  end subroutine DD2

end module HYPER_DERIVS_21
!#######################################################################
!#######################################################################
! ###
! ###
! ###                 E N D  O F  M O D U L E
! ###
!#######################################################################
!#######################################################################


!-----------------------------------------------------------------------
!
!    SUBROUTINE D2_NoMask
!
!    a simple first order derivative operator.  All points are second
!    order in accuracy.  Not to be used in evolutions.
!
!-----------------------------------------------------------------------
  subroutine D2_NoMask(Du, u, direction, par)
    use params
    implicit none

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT :: direction


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, j, k,nx, ny, nz, shp(3)
    logical, parameter :: ltrace = .false.

    if (ltrace) write(*,*) '...begin D2_NoMask.  direction = ',direction

    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    if (direction.eq.1) then
      Du(1,:,:) = (-3.0*u(1,:,:) + 4.0*u(2,:,:) - u(3,:,:))/(2*dx)
      do i = 2, nx-1
        Du(i,:,:) = (u(i+1,:,:) - u(i-1,:,:))/(2*dx)
      end do
      Du(nx,:,:) = (3.0*u(nx,:,:) - 4.0*u(nx-1,:,:) + u(nx-2,:,:))/(2*dx)
    else if (direction.eq.2) then
      Du(:,1,:) = (-3.0*u(:,1,:) + 4.0*u(:,2,:) - u(:,3,:))/(2*dy)
      do j = 2, ny-1
        Du(:,j,:) = (u(:,j+1,:)-u(:,j-1,:))/(2*dy)
      end do
      Du(:,ny,:) = (3.0*u(:,ny,:) - 4.0*u(:,ny-1,:) + u(:,ny-2,:))/(2*dy)
    else if (direction==3) then
      Du(:,:,1) = (-3.0*u(:,:,1) + 4.0*u(:,:,2) - u(:,:,3))/(2*dz)
      do k = 2, nz-1
        Du(:,:,k) = (u(:,:,k+1)-u(:,:,k-1))/(2*dz)
      end do
      Du(:,:,nz) = (3.0*u(:,:,nz) - 4.0*u(:,:,nz-1) + u(:,:,nz-2))/(2*dz)
    else
      write(*,*) "D2_NoMask: Wrong direction in the derivative"
      STOP
    end if


    return
  end subroutine D2_NoMask

!-----------------------------------------------------------------------
!
!    SUBROUTINE D43NoMask_simple
!
!    A standard fourth order derivative operator, 3rd order at the
!    boundaries.  This is the simple operator, and does not satisfy
!    summation by parts.
!
!-----------------------------------------------------------------------
  subroutine D43NoMask_simple(Du, u, direction, par)
    use params
    implicit none

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT :: direction


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, nx, ny, nz, shp(3)


    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    if (direction.eq.1) then
      !derivatives in x

      Du(1,:,:) = (-25.0d0*u(1,:,:) + 48.0d0*u(2,:,:) - 36.0d0*u(3,:,:) &
                  &  + 16.0d0*u(4,:,:) - 3.0d0*u(5,:,:)) / (12.0d0 * dx)
      Du(2,:,:) = (-3.0d0*u(1,:,:) - 10.0d0*u(2,:,:) + 18.0d0*u(3,:,:)  &
                  &  -6.0d0*u(4,:,:) + u(5,:,:)) / (12.0d0 * dx)
      do i = 3, nx-2
         Du(i,:,:) = (u(i-2,:,:)-u(i+2,:,:)+8.0d0*(u(i+1,:,:)-u(i-1,:,:))) &
                & / (12.0d0 * dx)
      end do
      Du(nx-1,:,:) = (-u(nx-4,:,:) + 6.0d0*u(nx-3,:,:) - 18.0d0*u(nx-2,:,:) &
                &  +10.0d0*u(nx-1,:,:) + 3.0d0*u(nx,:,:)) / (12.0d0*dx)
      Du(nx,:,:)   = (3.0d0*u(nx-4,:,:) - 16.0d0*u(nx-3,:,:) &
                &  + 36.0d0*u(nx-2,:,:) &
                &  -48.0d0*u(nx-1,:,:) + 25.0d0*u(nx,:,:)) / (12.0d0*dx)

    else if (direction.eq.2) then
      !derivatives in y

      Du(:,1,:) = (-25.0d0*u(:,1,:) + 48.0d0*u(:,2,:) - 36.0d0*u(:,3,:) &
             &  + 16.0d0*u(:,4,:) - 3.0d0*u(:,5,:)) / (12.0d0 * dy)
      Du(:,2,:) = (-3.0d0*u(:,1,:) - 10.0d0*u(:,2,:) + 18.0d0*u(:,3,:)  &  
             &  -6.0d0*u(:,4,:) + u(:,5,:)) / (12.0d0 * dy)

      do i = 3, ny-2 
         Du(:,i,:) = (u(:,i-2,:)-u(:,i+2,:)+8.0d0*(u(:,i+1,:)-u(:,i-1,:))) &
                & / (12.0d0 * dy)
      end do

      Du(:,ny-1,:) = (-u(:,ny-4,:) + 6.0d0*u(:,ny-3,:) - 18.0d0*u(:,ny-2,:) &
                &  +10.0d0*u(:,ny-1,:) + 3.0d0*u(:,ny,:)) / (12.0d0*dy)
      Du(:,ny,:)   = (3.0d0*u(:,ny-4,:) - 16.0d0*u(:,ny-3,:) &
                &  + 36.0d0*u(:,ny-2,:) &
                &  - 48.0d0*u(:,ny-1,:) + 25.0d0*u(:,ny,:)) / (12.0d0*dy)

    else if (direction.eq.3) then
      ! derivatives in z

      Du(:,:,1) = (-25.0d0*u(:,:,1) + 48.0d0*u(:,:,2) - 36.0d0*u(:,:,3) &
             &  + 16.0d0*u(:,:,4) - 3.0d0*u(:,:,5)) / (12.0d0 * dz)
      Du(:,:,2) = (-3.0d0*u(:,:,1) - 10.0d0*u(:,:,2) + 18.0d0*u(:,:,3)  &  
             &  -6.0d0*u(:,:,4) + u(:,:,5)) / (12.0d0 * dz)

      do i = 3, nz-2 
         Du(:,:,i) = (u(:,:,i-2)-u(:,:,i+2)+8.0d0*(u(:,:,i+1)-u(:,:,i-1))) &
                & / (12.0d0 * dz)
      end do

      Du(:,:,nz-1) = (-u(:,:,nz-4) + 6.0d0*u(:,:,nz-3) - 18.0d0*u(:,:,nz-2) &
                &  +10.0d0*u(:,:,nz-1) + 3.0d0*u(:,:,nz)) / (12.0d0*dz)
      Du(:,:,nz)   = (3.0d0*u(:,:,nz-4) - 16.0d0*u(:,:,nz-3) &
                &  + 36.0d0*u(:,:,nz-2) &
                &  - 48.0d0*u(:,:,nz-1) + 25.0d0*u(:,:,nz)) / (12.0d0*dz)

    else
      write(*,*)'Error in D43NoMask_simple.  Unknown direction ',&
                 direction
      stop
    end if

    return
  end subroutine D43NoMask_simple


  
  !-----------------------------------------------------------------------
!
!    SUBROUTINE D43NoMask_simple
!
!    A standard fourth order derivative operator, 3rd order at the
!    boundaries.  This is the simple operator, and does not satisfy
!    summation by parts.
!
!-----------------------------------------------------------------------
  subroutine D42NoMask_simple(Du, u, direction, par)
    use params
    implicit none

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT :: direction


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, j, k,nx, ny, nz, shp(3)


    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    if (direction.eq.1) then
      !derivatives in x

      Du(1,:,:) = (-48.d0*u(1,:,:) + 59.d0*u(2,:,:) - 8.d0*u(3,:,:) - 3.0d0*u(4,:,:))/(34.d0*dx)
      Du(2,:,:) = (u(3,:,:) - u(1,:,:))/(2.d0*dx)
      Du(3,:,:) = (8.d0*u(1,:,:) - 59.d0*u(2,:,:) + 59.d0*u(4,:,:) - 8.d0*u(5,:,:))/(86.d0*dx)
      Du(4,:,:) = (3.d0*u(1,:,:) - 59.d0*u(3,:,:) + 64.d0*u(5,:,:) - 8.d0*u(6,:,:))/(98.d0*dx)
      do i = 5, nx-4
        Du(i,:,:) = (-u(i+2,:,:) + 8.d0*u(i+1,:,:) - 8.d0*u(i-1,:,:) + u(i-2,:,:))/(12.d0*dx)
      end do
      Du(nx-3,:,:) = (-3.d0*u(nx,:,:) + 59.d0*u(nx-2,:,:) - 64.d0*u(nx-4,:,:) + 8.d0*u(nx-5,:,:))/(98.d0*dx)
      Du(nx-2,:,:) = (-8.d0*u(nx,:,:) + 59.d0*u(nx-1,:,:) - 59.d0*u(nx-3,:,:) + 8.d0*u(nx-4,:,:))/(86.d0*dx)
      Du(nx-1,:,:) = (u(nx,:,:) - u(nx-2,:,:))/(2.d0*dx)
      Du(nx,:,:)   = (48.d0*u(nx,:,:) - 59.d0*u(nx-1,:,:) + 8.d0*u(nx-2,:,:) + 3.0d0*u(nx-3,:,:))/(34.d0*dx)

    else if (direction.eq.2) then
      !derivatives in y

      Du(:,1,:) = (-48.d0*u(:,1,:) + 59.d0*u(:,2,:) - 8.d0*u(:,3,:) - 3.0d0*u(:,4,:))/(34.d0*dy)
      Du(:,2,:) = (u(:,3,:) - u(:,1,:))/(2.d0*dy)
      Du(:,3,:) = (8.d0*u(:,1,:) - 59.d0*u(:,2,:) + 59.d0*u(:,4,:) - 8.d0*u(:,5,:))/(86.d0*dy)
      Du(:,4,:) = (3.d0*u(:,1,:) - 59.d0*u(:,3,:) + 64.d0*u(:,5,:) - 8.d0*u(:,6,:))/(98.d0*dy)
      do j = 5, ny-4
        Du(:,j,:) = (-u(:,j+2,:) + 8.d0*u(:,j+1,:) - 8.d0*u(:,j-1,:) + u(:,j-2,:))/(12.d0*dy)
      end do
      Du(:,ny-3,:) = (-3.d0*u(:,ny,:) + 59.d0*u(:,ny-2,:) - 64.d0*u(:,ny-4,:) + 8.d0*u(:,ny-5,:))/(98.d0*dy)
      Du(:,ny-2,:) = (-8.d0*u(:,ny,:) + 59.d0*u(:,ny-1,:) - 59.d0*u(:,ny-3,:) + 8.d0*u(:,ny-4,:))/(86.d0*dy)
      Du(:,ny-1,:) = (u(:,ny,:) - u(:,ny-2,:))/(2.d0*dy)
      Du(:,ny,:)   = (48.d0*u(:,ny,:) - 59.d0*u(:,ny-1,:) + 8.d0*u(:,ny-2,:) + 3.0d0*u(:,ny-3,:))/(34.d0*dy)

    else if (direction.eq.3) then
      ! derivatives in z
      
      Du(:,:,1) = (-48.d0*u(:,:,1) + 59.d0*u(:,:,2) - 8.d0*u(:,:,3) - 3.0d0*u(:,:,4))/(34.d0*dz)
      Du(:,:,2) = (u(:,:,3) - u(:,:,1))/(2.d0*dz)
      Du(:,:,3) = (8.d0*u(:,:,1) - 59.d0*u(:,:,2) + 59.d0*u(:,:,4) - 8.d0*u(:,:,5))/(86.d0*dz)
      Du(:,:,4) = (3.d0*u(:,:,1) - 59.d0*u(:,:,3) + 64.d0*u(:,:,5) - 8.d0*u(:,:,6))/(98.d0*dz)
      do k = 5, nz-4
        Du(:,:,k) = (-u(:,:,k+2) + 8.d0*u(:,:,k+1) - 8.d0*u(:,:,k-1) + u(:,:,k-2))/(12.d0*dz)
      end do
      Du(:,:,nz-3) = (-3.d0*u(:,:,nz) + 59.d0*u(:,:,nz-2) - 64.d0*u(:,:,nz-4) + 8.d0*u(:,:,nz-5))/(98.d0*dz)
      Du(:,:,nz-2) = (-8.d0*u(:,:,nz) + 59.d0*u(:,:,nz-1) - 59.d0*u(:,:,nz-3) + 8.d0*u(:,:,nz-4))/(86.d0*dz)
      Du(:,:,nz-1) = (u(:,:,nz) - u(:,:,nz-2))/(2.d0*dz)
      Du(:,:,nz)   = (48.d0*u(:,:,nz) - 59.d0*u(:,:,nz-1) + 8.d0*u(:,:,nz-2) + 3.0d0*u(:,:,nz-3))/(34.d0*dz)

    else
      write(*,*)'Error in D43NoMask_simple.  Unknown direction ',&
                 direction
      stop
    end if

    return
  end subroutine D42NoMask_simple
