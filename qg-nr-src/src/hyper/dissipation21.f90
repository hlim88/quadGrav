
#include "cctk.h"

MODULE HYPER_DIS_21
use params

CONTAINS

!-----------------------------------------------------------------------
!
!    SUBROUTINE Dis21NoMask:
!
!    implements Dissipation.  2nd order interior, 1st order boundary
!    reuses the space allocated for the standard derivatives
!
!-----------------------------------------------------------------------
  subroutine Dis21NoMask(Du, u, direction, par)

    CCTK_REAL, DIMENSION(NPAR), INTENT(IN)      :: par
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
      Du(1,:,:) = 2.*(   u(2,:,:) & 
                       - u(1,:,:)  )/dx
      do i = 2, nx-1
        Du(i,:,:) =       ( u(i+1,:,:) & 
                    - 2. * u(i  ,:,:) & 
                    +      u(i-1,:,:) )/dx
      end do
      Du(nx,:,:) = - 2. * (   u(nx  ,:,:) & 
                            - u(nx-1,:,:)  )/dx
    
    else if (direction.eq.2) then 
      Du(:,1,:) = 2.*(   u(:,2,:) & 
                       - u(:,1,:)  )/dy
      do j = 2, ny-1
        Du(:,j,:) =       ( u(:,j+1,:) & 
                    - 2. * u(:,j  ,:) & 
                    +      u(:,j-1,:) )/dy
      end do
      Du(:,ny,:) = -2. * (   u(:,ny  ,:) & 
                           - u(:,ny-1,:)  )/dy
    
    else if (direction==3) then 
      Du(:,:,1) = 2. * (   u(:,:,2) & 
                         - u(:,:,1) )/dz
      do k = 2, nz-1
        Du(:,:,k) =      (  u(:,:,k+1) &  
                    - 2. * u(:,:,k  ) & 
                    +      u(:,:,k-1) )/dz
      end do
      Du(:,:,nz) = - 2. * (   u(:,:,nz  ) & 
                            - u(:,:,nz-1) )/dz
    
    else 
      write(*,*) "D21NoMask: Wrong direction in the derivative" 
      STOP
    
    end if 

    return
  end subroutine Dis21NoMask


!-----------------------------------------------------------------------
!
!    Dis21
!
!-----------------------------------------------------------------------
  subroutine Dis21(Du, u, direction, imask, par)
    implicit none

    CCTK_REAL, DIMENSION(NPAR), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT                                 :: direction
    CCTK_INT, dimension(:,:,:)               :: imask
    CCTK_INT                                 :: d


    ! local vars
    CCTK_REAL :: dx, dy, dz, c
    CCTK_INT  :: i, j, k,nx, ny, nz, shp(3)

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
      c=1./dx
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = imask(i,j,k) / 10000
        if (d .eq. P_STENCIL_CENTER) then
          Du(i,j,k) = c*(u(i+1,j,k)-2.0*u(i,j,k)+u(i-1,j,k))
        else if (d .eq. P_STENCIL_CENTER_RIGHT .or. &
                 d .eq. P_STENCIL_CENTER_LEFT) then
          Du(i,j,k) = c*(u(i+1,j,k)-2.0*u(i,j,k)+u(i-1,j,k))
        else if (d .eq. P_STENCIL_LEFT) then
          Du(i,j,k) = c*(2.0*(u(i+1,j,k) - u(i,j,k)))
        else if (d .eq. P_STENCIL_RIGHT) then
          Du(i,j,k) = c*(-2.0*(u(i,j,k) - u(i-1,j,k)))
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = c*(u(i+1,j,k)-2.0*u(i,j,k)+u(i-1,j,k))   ! not correct 
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = c*(u(i+1,j,k)-2.0*u(i,j,k)+u(i-1,j,k))   ! not correct
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = c*(u(i+1,j,k)-2.0*u(i,j,k)+u(i-1,j,k))   ! not correct
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = c*(u(i+1,j,k)-2.0*u(i,j,k)+u(i-1,j,k))   ! not correct
        else if (d .eq. P_STENCIL_UNDEF) then
          Du(i,j,k) = 0.0
        else
          if (ltrace2) then
            write(0,*)'diss21: Illegal x derivative type ',d
            write(0,*)'diss21: imask                     ',imask(i,j,k)
            write(0,*)'diss21: at (i,j,k) location: ',i,j,k
          end if
          Du(i,j,k) = 0.0
        end if
      end do
      end do
      end do

    else if (direction.eq.2) then

        c = 1./dy
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod((imask(i,j,k)/100),100)
        if (d .eq. P_STENCIL_CENTER) then
          Du(i,j,k) = c*(u(i,j+1,k) -2.0*u(i,j,k) + u(i,j-1,k))
        else if (d .eq. P_STENCIL_LEFT) then
          Du(i,j,k) = c*(2.0*(u(i,j+1,k) - u(i,j,k)))
        else if (d .eq. P_STENCIL_CENTER_LEFT .or. &
                 d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = c*(u(i,j+1,k) -2.0*u(i,j,k) + u(i,j-1,k))
        else if (d .eq. P_STENCIL_RIGHT) then
          Du(i,j,k) = c*(-2.0*(u(i,j,k) - u(i,j-1,k)))
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = c*(u(i,j+1,k) -2.0*u(i,j,k) + u(i,j-1,k))   ! not correct
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = c*(u(i,j+1,k) -2.0*u(i,j,k) + u(i,j-1,k))   ! not correct
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = c*(u(i,j+1,k) -2.0*u(i,j,k) + u(i,j-1,k))   ! not correct
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = c*(u(i,j+1,k) -2.0*u(i,j,k) + u(i,j-1,k))   ! not correct
        else if (d .eq. P_STENCIL_UNDEF) then
          Du(i,j,k) = 0.0
        else
          if (ltrace2) then
            write(0,*)'diss21: Illegal y derivative type ',d
            write(0,*)'diss21: imask                     ',imask(i,j,k)
            write(0,*)'diss21: at (i,j,k) location ',i,j,k
          end if
          Du(i,j,k) = 0.0
        end if
      end do
      end do
      end do

    else if (direction.eq.3) then

        c = 1./dz
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod(imask(i,j,k),100)
        if (d .eq. P_STENCIL_CENTER) then
          Du(i,j,k) = c*(u(i,j,k+1) - 2.0*u(i,j,k) + u(i,j,k-1))
        else if (d .eq. P_STENCIL_LEFT) then
          Du(i,j,k) = c*(2.0*(u(i,j,k+1) - u(i,j,k)))
        else if (d .eq. P_STENCIL_CENTER_LEFT .or. &
                 d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = c*(u(i,j,k+1) - 2.0*u(i,j,k) + u(i,j,k-1))
        else if (d .eq. P_STENCIL_RIGHT) then
          Du(i,j,k) = c*(-2.0*(u(i,j,k) - u(i,j,k-1)))
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = c*(u(i,j,k+1) - 2.0*u(i,j,k) + u(i,j,k-1))   ! not correct
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = c*(u(i,j,k+1) - 2.0*u(i,j,k) + u(i,j,k-1)  ) ! not correct
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = c*(u(i,j,k+1) - 2.0*u(i,j,k) + u(i,j,k-1) )  ! not correct
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = c*(u(i,j,k+1) - 2.0*u(i,j,k) + u(i,j,k-1))   ! not correct
        else if (d .eq. P_STENCIL_UNDEF) then
          Du(i,j,k) = 0.0
        else
          if (ltrace2) then
            write(0,*)'diss21: Illegal z derivative type ',d
            write(0,*)'diss21: imask                     ',imask(i,j,k)
            write(0,*)'diss21: at (i,j,k) location ',i,j,k
          end if
          Du(i,j,k) = 0.0
        end if
      end do
      end do
      end do

    else
      write(0,*)'d21:  unknown direction ',direction
    end if

    return
  end subroutine Dis21

!-----------------------------------------------------------------------
!
!    SUBROUTINE KODis21NoMask:
!
!    implements Kreiss-Oliger Dissipation.
!    reuses the space allocated for the standard derivatives
!
!-----------------------------------------------------------------------
  subroutine KODis21NoMask(Du, u, direction, par)

    CCTK_REAL, DIMENSION(NPAR), INTENT(IN)      :: par
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
       Du(1,:,:) = - 2.0 * (       u(3,:,:) &  
                             - 2.0*u(2,:,:) &  
                             +     u(1,:,:)  ) / dx

       Du(2,:,:) = - (       u(4,:,:) &  
                       - 4.0*u(3,:,:) &  
                       + 5.0*u(2,:,:) &
                       - 2.0*u(1,:,:)  ) /dx

       do i = 3, nx-2
          Du(i,:,:) = - (       u(i+2,:,:) &  
                          - 4.0*u(i+1,:,:) &  
                          + 6.0*u(i  ,:,:) &
                          - 4.0*u(i-1,:,:) &  
                          +     u(i-2,:,:)  ) / dx
       end do

       Du(nx-1,:,:) = - (       u(nx-3,:,:) &  
                          - 4.0*u(nx-2,:,:) &  
                          + 5.0*u(nx-1,:,:) &
                          - 2.0*u(nx  ,:,:)  ) / dx

       Du(nx,:,:) = - 2.0 * (       u(nx-2,:,:) &  
                              - 2.0*u(nx-1,:,:) &  
                                  + u(nx  ,:,:)  ) / dx

    else if (direction.eq.2) then
       Du(:,1,:) = - 2.0 * (       u(:,3,:) &  
                             - 2.0*u(:,2,:) &  
                                 + u(:,1,:) ) / dy

       Du(:,2,:) = - (       u(:,4,:) &  
                       - 4.0*u(:,3,:) &  
                       + 5.0*u(:,2,:) &
                       - 2.0*u(:,1,:)  ) /dy

       do j = 3, ny-2
          Du(:,j,:) = - (       u(:,j+2,:) &  
                          - 4.0*u(:,j+1,:) &  
                          + 6.0*u(:,j  ,:) &
                          - 4.0*u(:,j-1,:) &  
                          +     u(:,j-2,:)  ) / dy
       end do

       Du(:,ny-1,:) = - (       u(:,ny-3,:) &  
                          - 4.0*u(:,ny-2,:) &  
                          + 5.0*u(:,ny-1,:) &
                          - 2.0*u(:,ny  ,:)  ) / dy

       Du(:,ny,:) = - 2.0 * (       u(:,ny-2,:) &  
                              - 2.0*u(:,ny-1,:) &  
                              +     u(:,ny  ,:)  ) / dy


    else if (direction==3) then

       Du(:,:,1) = - 2.0 * (       u(:,:,3) &  
                             - 2.0*u(:,:,2) &  
                             +     u(:,:,1)  ) / dz

       Du(:,:,2) = - (       u(:,:,4) &  
                       - 4.0*u(:,:,3) &  
                       + 5.0*u(:,:,2) &
                       - 2.0*u(:,:,1)  ) /dz

       do k = 3, nz-2
          Du(:,:,k) = - (       u(:,:,k+2) &  
                          - 4.0*u(:,:,k+1) &  
                          + 6.0*u(:,:,k  ) &
                          - 4.0*u(:,:,k-1) &  
                          +     u(:,:,k-2)  ) / dz
       end do

       Du(:,:,nz-1) = - (       u(:,:,nz-3) &  
                          - 4.0*u(:,:,nz-2) &  
                          + 5.0*u(:,:,nz-1) &
                          - 2.0*u(:,:,nz  )  ) / dz

       Du(:,:,nz) = - 2.0 * (       u(:,:,nz-2) &  
                              - 2.0*u(:,:,nz-1) &  
                              +     u(:,:,nz  )  ) / dz

    else
      write(*,*) "D21NoMask: Wrong direction in the derivative"
      STOP
    end if

    return
  end subroutine KODis21NoMask


!!-----------------------------------------------------------------------
!!
!!    SUBROUTINE KODis42NoMask:
!!
!!    implements 3-6 Dissipation.
!!    reuses the space allocated for the standard derivatives
!!
!!-----------------------------------------------------------------------
!  subroutine KODis42NoMask(Du, u, direction, par)
! 
!    CCTK_REAL, DIMENSION(NPAR) , INTENT(IN)  :: par
!    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
!    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
!    CCTK_INT  :: direction
! 
! 
!    ! local vars
!    CCTK_REAL :: dx, dy, dz
!    CCTK_INT  :: i, j, k, nx, ny, nz, shp(3)
! 
!
!    shp = shape(u)
!    nx = shp(1)
!    ny = shp(2)
!    nz = shp(3)
!
!    dx = par(P_DX)
!    dy = par(P_DY)
!    dz = par(P_DZ)
!
!    if (direction.eq.1) then
!       Du(1,:,:) = -2.0 * (       u(3,:,:) 
!                            - 2.0*u(2,:,:) 
!                            +     u(1,:,:)  ) / dx
!
!       Du(2,:,:) = - (       u(4,:,:) & 
!                       - 4.0*u(3,:,:) & 
!                       + 5.0*u(2,:,:) &
!                       - 2.0*u(1,:,:)   ) /dx
!
!       do i = 4, nx-3
!          Du(i,:,:) = pre_factor_6 
!                      * ( -      u(i-3,:,:) & 
!                          +  6.0*u(i-2,:,:) & 
!                          - 15.0*u(i-1,:,:) & 
!                          + 20.0*u(i  ,:,:) &
!                          - 15.0*u(i+1,:,:) & 
!                          +  6.0*u(i+2,:,:) & 
!                          -      u(i+3,:,:)   ) 
!       end do
!
!       Du(nx-1,:,:) = - (       u(nx-3,:,:) &  
!                          - 4.0*u(nx-2,:,:) & 
!                          + 5.0*u(nx-1,:,:) &
!                          - 2.0*u(nx  ,:,:)   ) / dx
!
!       Du(nx,:,:) = -2.0 * (       u(nx-2,:,:) &  
!                             - 2.0*u(nx-1,:,:) &  
!                             +     u(nx  ,:,:)   ) / dx
!
!    else if (direction.eq.2) then
!       Du(:,1,:) = -2.0 * (       u(:,3,:) &  
!                            - 2.0*u(:,2,:) &  
!                            +     u(:,1,:)   ) / dy
!
!       Du(:,2,:) = - (       u(:,4,:) & 
!                       - 4.0*u(:,3,:) & 
!                       + 5.0*u(:,2,:) &
!                       - 2.0*u(:,1,:)   ) /dy
!
!       do j = 3, ny-2
!          Du(:,j,:) = - (u(:,j+2,:) - 4.0*u(:,j+1,:) + 6.0*u(:,j,:) &
!                             - 4.0*u(:,j-1,:) + u(:,j-2,:)) / dy
!       end do
!
!       Du(:,ny-1,:) = - (u(:,ny-3,:) - 4.0*u(:,ny-2,:) + 5.0*u(:,ny-1,:) &
!                                - 2.0*u(:,ny,:)) / dy
!
!       Du(:,ny,:) = -2.0 * (u(:,ny-2,:) - 2.0*u(:,ny-1,:) + u(:,ny,:)) / dy
!
!
!    else if (direction==3) then
!
!       Du(:,:,1) = -2.0 * (u(:,:,3) - 2.0*u(:,:,2) + u(:,:,1)) / dz
!
!       Du(:,:,2) = - (u(:,:,4) - 4.0*u(:,:,3) + 5.0*u(:,:,2) &
!                                - 2.0*u(:,:,1)) /dz
!
!       do k = 3, nz-2
!          Du(:,:,k) = - (u(:,:,k+2) - 4.0*u(:,:,k+1) + 6.0*u(:,:,k) &
!                             - 4.0*u(:,:,k-1) + u(:,:,k-2)) / dz
!       end do
!
!       Du(:,:,nz-1) = - (u(:,:,nz-3) - 4.0*u(:,:,nz-2) + 5.0*u(:,:,nz-1) &
!                                - 2.0*u(:,:,nz)) / dz
!
!       Du(:,:,nz) = -2.0 * (u(:,:,nz-2) - 2.0*u(:,:,nz-1) + u(:,:,nz)) / dz
!
!    else
!      write(*,*) "D21NoMask: Wrong direction in the derivative"
!      STOP
!    end if
!
!    return
!  end subroutine KODis42NoMask
 

!-----------------------------------------------------------------------
!
!    KODis21
!
!-----------------------------------------------------------------------
  subroutine KODis21(Du, u, direction, imask, par)
    implicit none

    CCTK_REAL, DIMENSION(NPAR), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT                                 :: direction
    CCTK_INT, dimension(:,:,:)               :: imask
    CCTK_INT                                 :: d


    ! local vars
    CCTK_REAL :: dx, dy, dz, c
    CCTK_INT  :: i, j, k,nx, ny, nz, shp(3)

    !  If ltrace2 is TRUE, then a warning message is printed for
    !  undefined derivatives are flagged.
    !  If ltrace2 is FALSE, then no warning message is printed.
    logical, parameter :: ltrace2 = .true.


    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)


    if (direction.eq.1) then
      c = 1.0/dx
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = imask(i,j,k) / 10000
        if (d .eq. P_STENCIL_CENTER) then
          Du(i,j,k) = -c * (u(i+2,j,k) - 4.0*u(i+1,j,k) + 6.0*u(i,j,k) &
                                - 4.0*u(i-1,j,k) + u(i-2,j,k))
        else if (d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = -c * (u(i-2,j,k) - 4.0*u(i-1,j,k) + 5.0*u(i,j,k) &
                                - 2.0*u(i+1,j,k))
        else if (d .eq. P_STENCIL_CENTER_LEFT) then
          Du(i,j,k) = -c * (u(i+2,j,k) - 4.0*u(i+1,j,k) + 5.0*u(i,j,k) &
                                - 2.0*u(i-1,j,k))
        else if (d .eq. P_STENCIL_LEFT) then
          Du(i,j,k) = -2.0 * c * (u(i+2,j,k) - 2.0*u(i+1,j,k) + u(i,j,k))
        else if (d .eq. P_STENCIL_RIGHT) then
          Du(i,j,k) = -2.0 * c * (u(i-2,j,k) - 2.0*u(i-1,j,k) + u(i,j,k))
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = -2.0/3.0*c* (2.0*u(i+2,j,k) - 4.0*u(i+1,j,k) &
                             + 3.0*u(i,j,k) - 2.0*u(i-1,j,k) + u(i-2,j,k))
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = -2.0/3.0*c* (u(i+2,j,k) - 2.0*u(i+1,j,k) &
                             + 3.0*u(i,j,k) - 4.0*u(i-1,j,k) + 2.0*u(i-2,j,k))
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = -2.0/7.0*c * (4.0*u(i+2,j,k) - 8.0*u(i+1,j,k)  &
                             + 7.0*u(i,j,k) - 6.0*u(i-1,j,k) + 3.0*u(i-2,j,k))
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = -2.0/7.0*c * (3.0*u(i+2,j,k) - 6.0*u(i+1,j,k)  &
                             + 7.0*u(i,j,k) - 8.0*u(i-1,j,k) + 4.0*u(i-2,j,k))
        else if (d .eq. P_STENCIL_UNDEF) then
          Du(i,j,k) = 0.0
        else
          if (ltrace2) then
            write(0,*)'KOdis21: invalid x derivative type ',d
            write(0,*)'     at location ',i,j,k
          end if
          Du(i,j,k) = 0.0
        end if
      end do
      end do
      end do
    end if

    if (direction.eq.2) then
      c = 1.0/dy
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod((imask(i,j,k)/100),100)
        if (d .eq. P_STENCIL_CENTER) then
          Du(i,j,k) = -c * (u(i,j+2,k) - 4.0*u(i,j+1,k) + 6.0*u(i,j,k) &
                                - 4.0*u(i,j-1,k) + u(i,j-2,k))
        else if (d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = -c * (u(i,j-2,k) - 4.0*u(i,j-1,k) + 5.0*u(i,j,k) &
                                - 2.0*u(i,j+1,k))
        else if (d .eq. P_STENCIL_CENTER_LEFT) then
          Du(i,j,k) = -c * (u(i,j+2,k) - 4.0*u(i,j+1,k) + 5.0*u(i,j,k) &
                                - 2.0*u(i,j-1,k))
        else if (d .eq. P_STENCIL_LEFT) then
          Du(i,j,k) = -2.0 * c * (u(i,j+2,k) - 2.0*u(i,j+1,k) + u(i,j,k))
        else if (d .eq. P_STENCIL_RIGHT) then
          Du(i,j,k) = -2.0 * c * (u(i,j-2,k) - 2.0*u(i,j-1,k) + u(i,j,k))
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = -2.0/3.0*c* (2.0*u(i,j+2,k) - 4.0*u(i,j+1,k) &
                             + 3.0*u(i,j,k) - 2.0*u(i,j-1,k) + u(i,j-2,k))
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = -2.0/3.0*c* (u(i,j+2,k) - 2.0*u(i,j+1,k) &
                             + 3.0*u(i,j,k) - 4.0*u(i,j-1,k) + 2.0*u(i,j-2,k))
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = -2.0/7.0*c * (4.0*u(i,j+2,k) - 8.0*u(i,j+1,k)  &
                             + 7.0*u(i,j,k) - 6.0*u(i,j-1,k) + 3.0*u(i,j-2,k))
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = -2.0/7.0*c * (3.0*u(i,j+2,k) - 6.0*u(i,j+1,k)  &
                             + 7.0*u(i,j,k) - 8.0*u(i,j-1,k) + 4.0*u(i,j-2,k))
        else if (d .eq. P_STENCIL_UNDEF) then
          Du(i,j,k) = 0.0
        else
          if (ltrace2) then
            write(0,*)'KOdis21: invalid y derivative type ',d
            write(0,*)'     at location ',i,j,k
          end if
          Du(i,j,k) = 0.0
        end if
      end do
      end do
      end do
    end if

    if (direction.eq.3) then
      c = 1.0/dz
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod(imask(i,j,k),100)
        if (d .eq. P_STENCIL_CENTER) then
          Du(i,j,k) = -c * (u(i,j,k+2) - 4.0*u(i,j,k+1) + 6.0*u(i,j,k) &
                                - 4.0*u(i,j,k-1) + u(i,j,k-2))
        else if (d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = -c * (u(i,j,k-2) - 4.0*u(i,j,k-1) + 5.0*u(i,j,k) &
                                - 2.0*u(i,j,k+1))
        else if (d .eq. P_STENCIL_CENTER_LEFT) then
          Du(i,j,k) = -c * (u(i,j,k+2) - 4.0*u(i,j,k+1) + 5.0*u(i,j,k) &
                                - 2.0*u(i,j,k-1))
        else if (d .eq. P_STENCIL_LEFT) then
          Du(i,j,k) = -2.0 * c * (u(i,j,k+2) - 2.0*u(i,j,k+1) + u(i,j,k))
        else if (d .eq. P_STENCIL_RIGHT) then
          Du(i,j,k) = -2.0 * c * (u(i,j,k-2) - 2.0*u(i,j,k-1) + u(i,j,k))
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_EDGE_P) then
          Du(i,j,k) = -2.0/3.0*c* (2.0*u(i,j,k+2) - 4.0*u(i,j,k+1) &
                             + 3.0*u(i,j,k) - 2.0*u(i,j,k-1) + u(i,j,k-2))
        else if (d .eq. P_STENCIL_EDGE_M) then
          Du(i,j,k) = -2.0/3.0*c* (u(i,j,k+2) - 2.0*u(i,j,k+1) &
                             + 3.0*u(i,j,k) - 4.0*u(i,j,k-1) + 2.0*u(i,j,k-2))
        else if (d .eq. P_STENCIL_VERTEX_P) then
          Du(i,j,k) = -2.0/7.0*c * (4.0*u(i,j,k+2) - 8.0*u(i,j,k+1)  &
                             + 7.0*u(i,j,k) - 6.0*u(i,j,k-1) + 3.0*u(i,j,k-2))
        else if (d .eq. P_STENCIL_VERTEX_M) then
          Du(i,j,k) = -2.0/7.0*c * (3.0*u(i,j,k+2) - 6.0*u(i,j,k+1)  &
                             + 7.0*u(i,j,k) - 8.0*u(i,j,k-1) + 4.0*u(i,j,k-2))
        else if (d .eq. P_STENCIL_UNDEF) then
          Du(i,j,k) = 0.0
        else
          if (ltrace2) then
            write(0,*)'KOdis21: invalid z derivative type ',d
            write(0,*)' at location ',i,j,k
	    write(0,*) 'more', nx,ny,nz,imask(i,j,k)
          end if
          Du(i,j,k) = 0.0
        end if
      end do
      end do
      end do
    end if




    return
  end subroutine KODis21

  
!-----------------------------------------------------------------------
!
!    KODis42
!
!-----------------------------------------------------------------------
  subroutine KODis42(Du, u, direction, imask, par)
    implicit none

    CCTK_REAL, DIMENSION(NPAR), INTENT(IN)   :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT                                 :: direction
    CCTK_INT, dimension(:,:,:)               :: imask
    CCTK_INT                                 :: d


    ! local vars
    CCTK_REAL :: dx, dy, dz, c
    CCTK_INT  :: i, j, k, nx, ny, nz, shp(3)
    CCTK_REAL :: smr3,smr2,smr1,spr1,spr2,spr3


    !  If ltrace2 is TRUE, then a warning message is printed for
    !  undefined derivatives are flagged.
    !  If ltrace2 is FALSE, then no warning message is printed.
    logical, parameter :: ltrace2 = .true.


    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    !local weights for the derivatives
    smr3=59.0/48.0;smr2=43.0/48.0;smr1=49.0/48.0
    spr3=smr3;spr2=smr2;spr1=smr1

    if (direction.eq.1) then
      c = 1.0/dx
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = imask(i,j,k) / 10000

        if (d .eq. P_STENCIL_CENTER) then
	  Du(i,j,k) = c*(u(i+3,j,k)+u(i-3,j,k)-6.*(u(i+2,j,k)+u(i-2,j,k))+ 15.*(u(i+1,j,k)+u(i-1,j,k))-20.*u(i,j,k))
        else if (d .eq. P_STENCIL_LEFT) then
	  Du(i,j,k) = c*(0.0)
	else if (d .eq. P_STENCIL_CENTER_LEFT) then
          Du(i,j,k) = c*(u(i+3,j,k)-3.*u(i+2,j,k)+3.*u(i+1,j,k)-u(i,j,k))/smr3
	else if (d .eq. P_STENCIL_CENTER_LEFT1) then
          Du(i,j,k) = c*(u(i+3,j,k)-6.*u(i+2,j,k)+12.*u(i+1,j,k)-10.*u(i,j,k)+3.*u(i-1,j,k))/smr2
	else if (d .eq. P_STENCIL_CENTER_LEFT2) then
          Du(i,j,k) = c*(u(i+3,j,k)-6.*u(i+2,j,k)+15.*u(i+1,j,k)-19.*u(i,j,k)+12.*u(i-1,j,k)-3.*u(i-2,j,k))/smr1
	else if (d .eq. P_STENCIL_CENTER_RIGHT2) then
          Du(i,j,k) = c*(u(i-3,j,k)-6.*u(i-2,j,k)+15.*u(i-1,j,k)-19.*u(i,j,k)+12.*u(i+1,j,k)-3.*u(i+2,j,k))/spr1
	else if (d .eq. P_STENCIL_CENTER_RIGHT1) then
          Du(i,j,k) = c*(u(i-3,j,k)-6.*u(i-2,j,k)+12.*u(i-1,j,k)-10.*u(i,j,k)+3.*u(i+1,j,k))/spr2
	else if (d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = c*(u(i-3,j,k)-3.*u(i-2,j,k)+3.*u(i-1,j,k)-u(i,j,k))/spr3
	else if (d .eq. P_STENCIL_RIGHT) then
	  Du(i,j,k) =  c*(0.0)     
	else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_UNDEF) then
          Du(i,j,k) = 0.0
        else
          if (ltrace2) then
            write(0,*)'KOdis42: invalid x derivative type ',d
            write(0,*)'     at location ',i,j,k
	    write(0,*) 'more', nx,ny,nz,imask(i,j,k)
          end if
          Du(i,j,k) = 0.0
        end if

      end do
      end do
      end do
    end if


    if (direction.eq.2) then
      c = 1.0/dy
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod((imask(i,j,k)/100),100)

        if (d .eq. P_STENCIL_CENTER) then
	  Du(i,j,k) = c*(u(i,j+3,k)+u(i,j-3,k)-6.*(u(i,j+2,k)+u(i,j-2,k))+ 15.*(u(i,j+1,k)+u(i,j-1,k))-20.*u(i,j,k))
        else if (d .eq. P_STENCIL_LEFT) then
	  Du(i,j,k) = c*(0.0)
        else if (d .eq. P_STENCIL_CENTER_LEFT) then
          Du(i,j,k) = c*(u(i,j+3,k)-3.*u(i,j+2,k)+3.*u(i,j+1,k)-u(i,j,k))/smr3
	else if (d .eq. P_STENCIL_CENTER_LEFT1) then
          Du(i,j,k) = c*(u(i,j+3,k)-6.*u(i,j+2,k)+12.*u(i,j+1,k)-10.*u(i,j,k)+3.*u(i,j-1,k))/smr2
	else if (d .eq. P_STENCIL_CENTER_LEFT2) then
          Du(i,j,k) = c*(u(i,j+3,k)-6.*u(i,j+2,k)+15.*u(i,j+1,k)-19.*u(i,j,k)+12.*u(i,j-1,k)-3.*u(i,j-2,k))/smr1
        else if (d .eq. P_STENCIL_CENTER_RIGHT2) then
          Du(i,j,k) = c*(u(i,j-3,k)-6.*u(i,j-2,k)+15.*u(i,j-1,k)-19.*u(i,j,k)+12.*u(i,j+1,k)-3.*u(i,j+2,k))/spr1
        else if (d .eq. P_STENCIL_CENTER_RIGHT1) then
          Du(i,j,k) = c*(u(i,j-3,k)-6.*u(i,j-2,k)+12.*u(i,j-1,k)-10.*u(i,j,k)+3.*u(i,j+1,k))/spr2
        else if (d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = c*(u(i,j-3,k)-3.*u(i,j-2,k)+3.*u(i,j-1,k)-u(i,j,k))/spr3
        else if (d .eq. P_STENCIL_RIGHT) then
	  Du(i,j,k) =  c*(0.0)
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_UNDEF) then
          Du(i,j,k) = 0.0
        else
          if (ltrace2) then
            write(0,*)'KOdis42: invalid y derivative type ',d
            write(0,*)'     at location ',i,j,k
	    write(0,*) 'more', nx,ny,nz,imask(i,j,k)
          end if
          Du(i,j,k) = 0.0
        end if

      end do
      end do
      end do
    end if


    if (direction.eq.3) then
      c = 1.0/dz
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod(imask(i,j,k),100)

        if (d .eq. P_STENCIL_CENTER) then
	  Du(i,j,k) = c*(u(i,j,k+3)+u(i,j,k-3)-6.*(u(i,j,k+2)+u(i,j,k-2))+ 15.*(u(i,j,k+1)+u(i,j,k-1))-20.*u(i,j,k))
        else if (d .eq. P_STENCIL_LEFT) then
	  Du(i,j,k) = c*(0.0)
        else if (d .eq. P_STENCIL_CENTER_LEFT) then
          Du(i,j,k) = c*(u(i,j,k+3)-3.*u(i,j,k+2)+3.*u(i,j,k+1)-u(i,j,k))/smr3
        else if (d .eq. P_STENCIL_CENTER_LEFT1) then
          Du(i,j,k) = c*(u(i,j,k+3)-6.*u(i,j,k+2)+12.*u(i,j,k+1)-10.*u(i,j,k)+3.*u(i,j,k-1))/smr2
        else if (d .eq. P_STENCIL_CENTER_LEFT2) then
          Du(i,j,k) = c*(u(i,j,k+3)-6.*u(i,j,k+2)+15.*u(i,j,k+1)-19.*u(i,j,k)+12.*u(i,j,k-1)-3.*u(i,j,k-2))/smr1
        else if (d .eq. P_STENCIL_CENTER_RIGHT2) then
          Du(i,j,k) = c*(u(i,j,k-3)-6.*u(i,j,k-2)+15.*u(i,j,k-1)-19.*u(i,j,k)+12.*u(i,j,k+1)-3.*u(i,j,k+2))/spr1
        else if (d .eq. P_STENCIL_CENTER_RIGHT1) then
          Du(i,j,k) = c*(u(i,j,k-3)-6.*u(i,j,k-2)+12.*u(i,j,k-1)-10.*u(i,j,k)+3.*u(i,j,k+1))/spr2
	else if (d .eq. P_STENCIL_CENTER_RIGHT) then
          Du(i,j,k) = c*(u(i,j,k-3)-3.*u(i,j,k-2)+3.*u(i,j,k-1)-u(i,j,k))/spr3
	else if (d .eq. P_STENCIL_RIGHT) then
	  Du(i,j,k) =  c*(0.0)
        else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_UNDEF) then
          Du(i,j,k) = 0.0
        else
          if (ltrace2) then
            write(0,*)'KOdis42: invalid z derivative type ',d
            write(0,*)'     at location ',i,j,k
	    write(0,*) 'more', nx,ny,nz,imask(i,j,k)
          end if
          Du(i,j,k) = 0.0
        end if

      end do
      end do
      end do
    end if

    return
  end subroutine KODis42
  
  






!-----------------------------------------------------------------------
!
!    SUBROUTINE KODis42NoMask:
!
!-----------------------------------------------------------------------
  subroutine KODis42NoMask(Du, u, direction, par)

    CCTK_REAL, DIMENSION(NPAR) , INTENT(IN)  :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT  :: direction


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, j, k, nx, ny, nz, shp(3)

    CCTK_REAL :: pre_factor_6_dx, pre_factor_6_dy, pre_factor_6_dz  
    CCTK_REAL :: smr3,smr2,smr1,spr1,spr2,spr3
    CCTK_INT :: lef, rig

    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    pre_factor_6_dx = -1.d0 / 64.d0 / dx  
    pre_factor_6_dy = -1.d0 / 64.d0 / dy  
    pre_factor_6_dz = -1.d0 / 64.d0 / dz  

    smr3=59.d0/48.d0*64*dx;smr2=43.d0/48.d0*64*dx;smr1=49.d0/48.d0*64*dx
    spr3=smr3;spr2=smr2;spr1=smr1


    if (direction.eq.1) then
    lef = 4
    rig = nx-3

       Du(1,:,:) =  (u(lef,:,:)-3.*u(lef-1,:,:)+3.*u(lef-2,:,:)-u(lef-3,:,:))/smr3
       Du(2,:,:) =  (u(lef+1,:,:)-6.*u(lef,:,:)+12.*u(lef-1,:,:)  &
      &                                 -10.*u(lef-2,:,:)+3.*u(lef-3,:,:))/smr2
       Du(3,:,:) =  (u(lef+2,:,:)-6.*u(lef+1,:,:)+15.*u(lef,:,:)  &
      &                                 -19.*u(lef-1,:,:)+12.*u(lef-2,:,:)-3.*u(lef-3,:,:))/smr1

       do i = lef, rig
          Du(i,:,:) = pre_factor_6_dx       &  
                      * ( -      u(i-3,:,:) & 
                          +  6.0*u(i-2,:,:) & 
                          - 15.0*u(i-1,:,:) & 
                          + 20.0*u(i  ,:,:) &
                          - 15.0*u(i+1,:,:) & 
                          +  6.0*u(i+2,:,:) & 
                          -      u(i+3,:,:)   ) 
       end do

       Du(nx-2,:,:) = (u(rig-2,:,:)-6.*u(rig-1,:,:)+15.*u(rig,:,:) &
     &                         -19.*u(rig+1,:,:)+12.*u(rig+2,:,:)-3.*u(rig+3,:,:))/spr1 
       Du(nx-1,:,:) = (u(rig-1,:,:)-6.*u(rig,:,:)+12.*u(rig+1,:,:) &
     &                          -10.*u(rig+2,:,:)+3.*u(rig+3,:,:))/spr2
       Du(nx  ,:,:) = (u(rig,:,:)-3.*u(rig+1,:,:)+3.*u(rig+2,:,:)-u(rig+3,:,:))/spr3

    else if (direction.eq.2) then
    smr3=59.d0/48.d0*64*dy;smr2=43.d0/48.d0*64*dy;smr1=49.d0/48.d0*64*dy
    spr3=smr3;spr2=smr2;spr1=smr1

    lef = 4
    rig = ny-3

       Du(:,1,:) =  (u(:,lef,:)-3.*u(:,lef-1,:)+3.*u(:,lef-2,:)-u(:,lef-3,:))/smr3
       Du(:,2,:) =  (u(:,lef+1,:)-6.*u(:,lef,:)+12.*u(:,lef-1,:)  &
      &                                 -10.*u(:,lef-2,:)+3.*u(:,lef-3,:))/smr2
       Du(:,3,:) =  (u(:,lef+2,:)-6.*u(:,lef+1,:)+15.*u(:,lef,:)  &
      &                                 -19.*u(:,lef-1,:)+12.*u(:,lef-2,:)-3.*u(:,lef-3,:))/smr1

       do j = lef, rig
          Du(:,j,:) = pre_factor_6_dy       &    
                      * ( -      u(:,j-3,:) & 
                          +  6.0*u(:,j-2,:) & 
                          - 15.0*u(:,j-1,:) & 
                          + 20.0*u(:,j  ,:) &
                          - 15.0*u(:,j+1,:) & 
                          +  6.0*u(:,j+2,:) & 
                          -      u(:,j+3,:)   ) 
       end do

       Du(:,ny-2,:) = (u(:,rig-2,:)-6.*u(:,rig-1,:)+15.*u(:,rig,:) &
     &                         -19.*u(:,rig+1,:)+12.*u(:,rig+2,:)-3.*u(:,rig+3,:))/spr1 
       Du(:,ny-1,:) = (u(:,rig-1,:)-6.*u(:,rig,:)+12.*u(:,rig+1,:) &
     &                          -10.*u(:,rig+2,:)+3.*u(:,rig+3,:))/spr2
       Du(:,ny  ,:) = (u(:,rig,:)-3.*u(:,rig+1,:)+3.*u(:,rig+2,:)-u(:,rig+3,:))/spr3


    else if (direction==3) then
    smr3=59.d0/48.d0*64*dz;smr2=43.d0/48.d0*64*dz;smr1=49.d0/48.d0*64*dz
    spr3=smr3;spr2=smr2;spr1=smr1

    lef = 4
    rig = nz-3

       Du(:,:,1) =  (u(:,:,lef)-3.*u(:,:,lef-1)+3.*u(:,:,lef-2)-u(:,:,lef-3))/smr3
       Du(:,:,2) =  (u(:,:,lef+1)-6.*u(:,:,lef)+12.*u(:,:,lef-1)  &
      &                                 -10.*u(:,:,lef-2)+3.*u(:,:,lef-3))/smr2
       Du(:,:,3) =  (u(:,:,lef+2)-6.*u(:,:,lef+1)+15.*u(:,:,lef)  &
      &                                 -19.*u(:,:,lef-1)+12.*u(:,:,lef-2)-3.*u(:,:,lef-3))/smr1

       do k = lef, rig 
          Du(:,:,k) = pre_factor_6_dz       &            
                      * ( -      u(:,:,k-3) & 
                          +  6.0*u(:,:,k-2) & 
                          - 15.0*u(:,:,k-1) & 
                          + 20.0*u(:,:,k  ) &
                          - 15.0*u(:,:,k+1) & 
                          +  6.0*u(:,:,k+2) & 
                          -      u(:,:,k+3)   ) 
       end do      
            
       Du(:,:,nz-2) = (u(:,:,rig-2)-6.*u(:,:,rig-1)+15.*u(:,:,rig) &
     &                         -19.*u(:,:,rig+1)+12.*u(:,:,rig+2)-3.*u(:,:,rig+3))/spr1 
       Du(:,:,nz-1) = (u(:,:,rig-1)-6.*u(:,:,rig)+12.*u(:,:,rig+1) &
     &                          -10.*u(:,:,rig+2)+3.*u(:,:,rig+3))/spr2
       Du(:,:,nz  ) = (u(:,:,rig)-3.*u(:,:,rig+1)+3.*u(:,:,rig+2)-u(:,:,rig+3))/spr3

    else
      write(*,*) "D21NoMask: Wrong direction in the derivative"
      STOP
    end if

    return
  end subroutine KODis42NoMask




!-----------------------------------------------------------------------
!
!    SUBROUTINE KODis642NoMask:
!
!-----------------------------------------------------------------------
  subroutine KODis642NoMask(Du, u, direction, par)

    CCTK_REAL, DIMENSION(NPAR) , INTENT(IN)  :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT  :: direction


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, j, k, nx, ny, nz, shp(3)

    CCTK_REAL :: pre_factor_8_dx, pre_factor_8_dy, pre_factor_8_dz  
    CCTK_REAL :: smr3,smr2,smr1,spr1,spr2,spr3,smr4,spr4 
    CCTK_INT :: lef, rig

    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    pre_factor_8_dx = + 1.d0 / 256.d0 / dx  
    pre_factor_8_dy = + 1.d0 / 256.d0 / dy  
    pre_factor_8_dz = + 1.d0 / 256.d0 / dz  

    smr4=17.d0/48.d0*256*dx;smr3=59.d0/48.d0*256*dx;smr2=43.d0/48.d0*256*dx;smr1=49.d0/48.d0*256*dx
    spr4=smr4;spr3=smr3;spr2=smr2;spr1=smr1

    if (direction.eq.1) then
        lef=5; rig=nx-4

        du(lef-4,:,:) = -(-u(lef,:,:)+4.*u(lef-1,:,:)-6.*u(lef-2,:,:)   &
     &                                  +4.*u(lef-3,:,:)-u(lef-4,:,:))/smr4
        du(lef-3,:,:) = -(2.*u(lef,:,:)-9*u(lef-1,:,:)+15.*u(lef-2,:,:)  &
     &                                  -11.*u(lef-3,:,:)+3.*u(lef-4,:,:))/smr3
        du(lef-2,:,:) = -(-u(lef+1,:,:)+3.*u(lef,:,:)-8.*u(lef-2,:,:)  &
     &                                  +9.*u(lef-3,:,:)-3.*u(lef-4,:,:))/smr2
        du(lef-1,:,:) = -(-u(lef+2,:,:)+6.*u(lef+1,:,:)-14.*u(lef,:,:)+15.*u(lef-1,:,:)   &
     &                                  -6.*u(lef-2,:,:)-u(lef-3,:,:)+u(lef-4,:,:))/smr1


       do i = lef, rig
          Du(i,:,:) = pre_factor_8_dx         & 
                      * (        u(i-4,:,:)   & 
                          -  8.0*u(i-3,:,:)   & 
                          + 28.0*u(i-2,:,:)   & 
                          - 56.0*u(i-1,:,:)   & 
                          + 70.0*u(i  ,:,:)   &
                          - 56.0*u(i+1,:,:)   & 
                          + 28.0*u(i+2,:,:)   & 
                          -  8.0*u(i+3,:,:)   & 
                          +      u(i+4,:,:)     ) 
       end do

        du(rig+4,:,:) = -(-u(rig,:,:)+4.*u(rig+1,:,:)-6.*u(rig+2,:,:)  &
     &                                     +4.*u(rig+3,:,:)-u(rig+4,:,:))/spr4
        du(rig+3,:,:) = -(2.*u(rig,:,:)-9*u(rig+1,:,:)+15.*u(rig+2,:,:)  &
     &                                     -11.*u(rig+3,:,:)+3.*u(rig+4,:,:))/spr3
        du(rig+2,:,:) = -(-u(rig-1,:,:)+3.*u(rig,:,:)-8.*u(rig+2,:,:)   &
     &                                      +9.*u(rig+3,:,:)-3.*u(rig+4,:,:))/spr2
        du(rig+1,:,:) = -(-u(rig-2,:,:)+6.*u(rig-1,:,:)-14.*u(rig,:,:)+15.*u(rig+1,:,:)   &
     &                                      -6.*u(rig+2,:,:)-u(rig+3,:,:)+u(rig+4,:,:))/spr1


    else if (direction.eq.2) then
        lef=5; rig=ny-4
    smr4=17.d0/48.d0*256*dy;smr3=59.d0/48.d0*256*dy;smr2=43.d0/48.d0*256*dy;smr1=49.d0/48.d0*256*dy
    spr4=smr4;spr3=smr3;spr2=smr2;spr1=smr1

        du(:,lef-4,:) = -(-u(:,lef,:)+4.*u(:,lef-1,:)-6.*u(:,lef-2,:)   &
     &                                  +4.*u(:,lef-3,:)-u(:,lef-4,:))/smr4
        du(:,lef-3,:) = -(2.*u(:,lef,:)-9*u(:,lef-1,:)+15.*u(:,lef-2,:)  &
     &                                  -11.*u(:,lef-3,:)+3.*u(:,lef-4,:))/smr3
        du(:,lef-2,:) = -(-u(:,lef+1,:)+3.*u(:,lef,:)-8.*u(:,lef-2,:)  &
     &                                  +9.*u(:,lef-3,:)-3.*u(:,lef-4,:))/smr2
        du(:,lef-1,:) = -(-u(:,lef+2,:)+6.*u(:,lef+1,:)-14.*u(:,lef,:)+15.*u(:,lef-1,:)   &
     &                                  -6.*u(:,lef-2,:)-u(:,lef-3,:)+u(:,lef-4,:))/smr1

       do j = lef, rig 
          Du(:,j,:) = pre_factor_8_dy         & 
                      * (        u(:,j-4,:)   & 
                          -  8.0*u(:,j-3,:)   & 
                          + 28.0*u(:,j-2,:)   & 
                          - 56.0*u(:,j-1,:)   & 
                          + 70.0*u(:,j  ,:)   &
                          - 56.0*u(:,j+1,:)   & 
                          + 28.0*u(:,j+2,:)   & 
                          -  8.0*u(:,j+3,:)   & 
                          +      u(:,j+4,:)     ) 
       end do

        du(:,rig+4,:) = -(-u(:,rig,:)+4.*u(:,rig+1,:)-6.*u(:,rig+2,:)  &
     &                                     +4.*u(:,rig+3,:)-u(:,rig+4,:))/spr4
        du(:,rig+3,:) = -(2.*u(:,rig,:)-9*u(:,rig+1,:)+15.*u(:,rig+2,:)  &
     &                                     -11.*u(:,rig+3,:)+3.*u(:,rig+4,:))/spr3
        du(:,rig+2,:) = -(-u(:,rig-1,:)+3.*u(:,rig,:)-8.*u(:,rig+2,:)   &
     &                                      +9.*u(:,rig+3,:)-3.*u(:,rig+4,:))/spr2
        du(:,rig+1,:) = -(-u(:,rig-2,:)+6.*u(:,rig-1,:)-14.*u(:,rig,:)+15.*u(:,rig+1,:)   &
     &                                      -6.*u(:,rig+2,:)-u(:,rig+3,:)+u(:,rig+4,:))/spr1


    else if (direction==3) then
        lef=5; rig=nz-4
    smr4=17.d0/48.d0*256*dz;smr3=59.d0/48.d0*256*dz;smr2=43.d0/48.d0*256*dz;smr1=49.d0/48.d0*256*dz
    spr4=smr4;spr3=smr3;spr2=smr2;spr1=smr1

        du(:,:,lef-4) = -(-u(:,:,lef)+4.*u(:,:,lef-1)-6.*u(:,:,lef-2)   &
     &                                  +4.*u(:,:,lef-3)-u(:,:,lef-4))/smr4
        du(:,:,lef-3) = -(2.*u(:,:,lef)-9*u(:,:,lef-1)+15.*u(:,:,lef-2)  &
     &                                  -11.*u(:,:,lef-3)+3.*u(:,:,lef-4))/smr3
        du(:,:,lef-2) = -(-u(:,:,lef+1)+3.*u(:,:,lef)-8.*u(:,:,lef-2)  &
     &                                  +9.*u(:,:,lef-3)-3.*u(:,:,lef-4))/smr2
        du(:,:,lef-1) = -(-u(:,:,lef+2)+6.*u(:,:,lef+1)-14.*u(:,:,lef)+15.*u(:,:,lef-1)   &
     &                                  -6.*u(:,:,lef-2)-u(:,:,lef-3)+u(:,:,lef-4))/smr1

       do k = lef, rig 
          Du(:,:,k) = pre_factor_8_dz         & 
                      * (        u(:,:,k-4)   & 
                          -  8.0*u(:,:,k-3)   & 
                          + 28.0*u(:,:,k-2)   & 
                          - 56.0*u(:,:,k-1)   & 
                          + 70.0*u(:,:,k  )   &
                          - 56.0*u(:,:,k+1)   & 
                          + 28.0*u(:,:,k+2)   & 
                          -  8.0*u(:,:,k+3)   & 
                          +      u(:,:,k+4)     ) 
       end do      
            
        du(:,:,rig+4) = -(-u(:,:,rig)+4.*u(:,:,rig+1)-6.*u(:,:,rig+2)  &
     &                                     +4.*u(:,:,rig+3)-u(:,:,rig+4))/spr4
        du(:,:,rig+3) = -(2.*u(:,:,rig)-9*u(:,:,rig+1)+15.*u(:,:,rig+2)  &
     &                                     -11.*u(:,:,rig+3)+3.*u(:,:,rig+4))/spr3
        du(:,:,rig+2) = -(-u(:,:,rig-1)+3.*u(:,:,rig)-8.*u(:,:,rig+2)   &
     &                                      +9.*u(:,:,rig+3)-3.*u(:,:,rig+4))/spr2
        du(:,:,rig+1) = -(-u(:,:,rig-2)+6.*u(:,:,rig-1)-14.*u(:,:,rig)+15.*u(:,:,rig+1)   &
     &                                      -6.*u(:,:,rig+2)-u(:,:,rig+3)+u(:,:,rig+4))/spr1

    else
      write(*,*) "D21NoMask: Wrong direction in the derivative"
      STOP
    end if

    return
  end subroutine KODis642NoMask










END MODULE HYPER_DIS_21
