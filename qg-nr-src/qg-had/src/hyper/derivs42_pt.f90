!-----------------------------------------------------------------------
!
!    $Id$
!
!-----------------------------------------------------------------------

#include "cctk.h"

MODULE HYPER_DERIVS_42_PT
use params
use hypercoords
implicit none

CONTAINS

!-----------------------------------------------------------------------
!
!    SUBROUTINE D42NoMask_pt:
!
!    implement Strand derivative.  4th order interior, 2nd order boundary
!    following equations (B1) and (B2) from gr-qc/0308007
!
!-----------------------------------------------------------------------
  subroutine D42NoMask_pt(Dxu, Dyu, Dzu, u, wpt, i, j, k, par)
    implicit none

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(NW), INTENT(IN)     :: wpt
    CCTK_REAL, INTENT(OUT)                   :: Dxu, Dyu, Dzu
    CCTK_INT :: i, j, k


    ! local vars
    CCTK_REAL :: dx, dy, dz 
    CCTK_REAL :: dxc, dyc, dzc, jac(3,3)
    CCTK_INT  :: nx, ny, nz, shp(3)
 
    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    
    !-------------------------------------------------
    !  x derivs
    !-------------------------------------------------
    if (i .eq. 1) then
      Dxu = (-48.d0*u(1,j,k) + 59.d0*u(2,j,k) - 8.d0*u(3,j,k) - 3.0d0*u(4,j,k))/(34.d0*dx)
    else if (i .eq. 2) then
      Dxu = (u(3,j,k) - u(1,j,k))/(2.d0*dx)
    else if (i .eq. 3) then
      Dxu = (8.d0*u(1,j,k) - 59.d0*u(2,j,k) + 59.d0*u(4,j,k) - 8.d0*u(5,j,k))/(86.d0*dx)
    else if (i .eq. 4) then
      Dxu = (3.d0*u(1,j,k) - 59.d0*u(3,j,k) + 64.d0*u(5,j,k) - 8.d0*u(6,j,k))/(98.d0*dx)
    else if (i .gt. 4 .and. i .lt. nx-3) then
      Dxu = (-u(i+2,j,k) + 8.d0*u(i+1,j,k) - 8.d0*u(i-1,j,k) + u(i-2,j,k))/(12.d0*dx)
    else if (i .eq. nx-3) then
      Dxu = (-3.d0*u(nx,j,k) + 59.d0*u(nx-2,j,k) - 64.d0*u(nx-4,j,k) + 8.d0*u(nx-5,j,k))/(98.d0*dx)
    else if (i .eq. nx-2) then
      Dxu = (-8.d0*u(nx,j,k) + 59.d0*u(nx-1,j,k) - 59.d0*u(nx-3,j,k) + 8.d0*u(nx-4,j,k))/(86.d0*dx)
    else if (i .eq. nx-1) then
      Dxu = (u(nx,j,k) - u(nx-2,j,k))/(2.d0*dx)
    else if (i .eq. nx) then
      Dxu = (48.d0*u(nx,j,k) - 59.d0*u(nx-1,j,k) + 8.d0*u(nx-2,j,k) + 3.0d0*u(nx-3,j,k))/(34.d0*dx)
    else
      write(0,*)'D42NOMASK_PT>> i out of range in x-deriv. i=',i
      STOP
    end if
    
    !-------------------------------------------------
    !  y derivs
    !-------------------------------------------------
    if (j .eq. 1) then
      Dyu = (-48.d0*u(i,1,k) + 59.d0*u(i,2,k) - 8.d0*u(i,3,k) - 3.0d0*u(i,4,k))/(34.d0*dy)
    else if (j .eq. 2) then
      Dyu = (u(i,3,k) - u(i,1,k))/(2.d0*dy)
    else if (j .eq. 3) then
      Dyu = (8.d0*u(i,1,k) - 59.d0*u(i,2,k) + 59.d0*u(i,4,k) - 8.d0*u(i,5,k))/(86.d0*dy)
    else if (j .eq. 4) then
      Dyu = (3.d0*u(i,1,k) - 59.d0*u(i,3,k) + 64.d0*u(i,5,k) - 8.d0*u(i,6,k))/(98.d0*dy)
    else if (j .gt. 4 .and. j .lt. ny-3) then
      Dyu = (-u(i,j+2,k) + 8.d0*u(i,j+1,k) - 8.d0*u(i,j-1,k) + u(i,j-2,k))/(12.d0*dy)
    else if (j .eq. ny-3) then
      Dyu = (-3.d0*u(i,ny,k) + 59.d0*u(i,ny-2,k) - 64.d0*u(i,ny-4,k) + 8.d0*u(i,ny-5,k))/(98.d0*dy)
    else if (j .eq. ny-2) then
      Dyu = (-8.d0*u(i,ny,k) + 59.d0*u(i,ny-1,k) - 59.d0*u(i,ny-3,k) + 8.d0*u(i,ny-4,k))/(86.d0*dy)
    else if (j .eq. ny-1) then
      Dyu = (u(i,ny,k) - u(i,ny-2,k))/(2.d0*dy)
    else if (j .eq. ny) then
      Dyu = (48.d0*u(i,ny,k) - 59.d0*u(i,ny-1,k) + 8.d0*u(i,ny-2,k) + 3.0d0*u(i,ny-3,k))/(34.d0*dy)
    else
      write(0,*)'D42NOMASK_PT>> j out of range in y-deriv. j=',j
      STOP
    end if
    
    !-------------------------------------------------
    !  z derivs
    !-------------------------------------------------
    if (k .eq. 1) then
      Dzu = (-48.d0*u(i,j,1) + 59.d0*u(i,j,2) - 8.d0*u(i,j,3) - 3.0d0*u(i,j,4))/(34.d0*dz)
    else if (k .eq. 2) then
      Dzu = (u(i,j,3) - u(i,j,1))/(2.d0*dz)
    else if (k .eq. 3) then
      Dzu = (8.d0*u(i,j,1) - 59.d0*u(i,j,2) + 59.d0*u(i,j,4) - 8.d0*u(i,j,5))/(86.d0*dz)
    else if (k .eq. 4) then
      Dzu = (3.d0*u(i,j,1) - 59.d0*u(i,j,3) + 64.d0*u(i,j,5) - 8.d0*u(i,j,6))/(98.d0*dz)
    else if (k .gt. 4 .and. k .lt. nz-3) then
      Dzu = (-u(i,j,k+2) + 8.d0*u(i,j,k+1) - 8.d0*u(i,j,k-1) + u(i,j,k-2))/(12.d0*dz)
    else if (k .eq. nz-3) then
      Dzu = (-3.d0*u(i,j,nz) + 59.d0*u(i,j,nz-2) - 64.d0*u(i,j,nz-4) + 8.d0*u(i,j,nz-5))/(98.d0*dz)
    else if (k .eq. nz-2) then
      Dzu = (-8.d0*u(i,j,nz) + 59.d0*u(i,j,nz-1) - 59.d0*u(i,j,nz-3) + 8.d0*u(i,j,nz-4))/(86.d0*dz)
    else if (k .eq. nz-1) then
      Dzu = (u(i,j,nz) - u(i,j,nz-2))/(2.d0*dz)
    else if (k .eq. nz) then
      Dzu = (48.d0*u(i,j,nz) - 59.d0*u(i,j,nz-1) + 8.d0*u(i,j,nz-2) + 3.0d0*u(i,j,nz-3))/(34.d0*dz)
    else
      write(0,*)'D42NOMASK_PT>> k out of range in z-deriv. k=',k
      STOP
    end if

    !
    ! Check to see if we are using fish-eye coords, and transform
    ! the derivatives is necessary.
    !
    if (P_ENABLE_ALT_COORDS .eq. 1 .and. &
    nint(par(P_ALT_COORD_TYPE)) .ne. 0) then
       if (nint(par(P_ALT_COORD_TYPE)) .eq. 1) then
          call  pcjac_pt(jac, wpt, par)
          dxc = Dxu
          dyc = Dyu
          dzc = Dzu
          Dxu = jac(1,1)*dxc + jac(1,2)*dyc + jac(1,3)*dzc
          Dyu = jac(2,1)*dxc + jac(2,2)*dyc + jac(2,3)*dzc
          Dzu = jac(3,1)*dxc + jac(3,2)*dyc + jac(3,3)*dzc
       else 
          write(0,*)'D42NoMask_PT: unknown alt_coord_type ',&
                     par(P_ALT_COORD_TYPE)
          stop
       end if
    end if


    return
  end subroutine D42NoMask_pt


!-----------------------------------------------------------------------
!
!    D42_pt
!
!  This routine has two debugging options:
!       debug_deriv -> if true, then the mask is checked to see that 
!                      all points used to calculate the derivative are
!                      unmasked, i.e., valid computational points.
!       ltrace2     -> Checks that the undefined derivative type occurs
!                      only in ghostzones, where it may occur when using MPI.
!
!  At times a derivative operator may be legitimately undefined.
!  This can happen in the ghost zones, if they are very close to
!  a mask boundary.  At such times the derivative is simply set
!  to zero.  For debugging purposes, the ltrace2 parameter can be
!  used to control whether an message is printed for undefined 
!  derivatives.
!
!-----------------------------------------------------------------------
  subroutine D42_pt(Dxu, Dyu,  Dzu,  u, wpt,  &
                    imask, i, j, k, par)
    implicit none

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, INTENT(OUT)                   :: Dxu, Dyu, Dzu
    CCTK_REAL, DIMENSION(NW), INTENT(IN)     :: wpt
    CCTK_INT                                 :: i, j, k
    CCTK_INT                                 :: imask


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_REAL :: dxc, dyc, dzc, jac(3,3)
    CCTK_INT  :: d, nx, ny, nz, shp(3)
    CCTK_INT  :: gz_xmin, gz_xmax, gz_ymin, gz_ymax, gz_zmin, gz_zmax
    CCTK_INT  :: bbox1, bbox2, bbox3, bbox4, bbox5, bbox6
 

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

    if (ltrace2) then
      !
      ! If ltrace2==.true., then we check that all undefined derivatives
      ! are hidden in ghostzones.  This is redundant, as the mask routine
      ! now also checks this.  But it is good to set ltrace2=.true. for 
      ! debugging now and then.
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
    end if

    !-----------------------------------------------------------------
    ! Sanity check
    !-----------------------------------------------------------------
    if (i .lt. 1 .or. i .gt. nx) then
      write(0,*)'D42_PT >> i out of range: i/nx = ',i,nx
    end if
    if (j .lt. 1 .or. j .gt. ny) then
      write(0,*)'D42_PT >> j out of range: y/ny = ',j,ny
    end if
    if (k .lt. 1 .or. k .gt. nz) then
      write(0,*)'D42_PT >> k out of range: k/nz = ',k,nz
    end if

    !-----------------------------------------------------------------
    ! x-derivatives
    !-----------------------------------------------------------------

    d = imask / 10000
    if (d .eq. P_STENCIL_CENTER) then
      Dxu = (-u(i+2,j,k) + 8.d0*u(i+1,j,k) - 8.d0*u(i-1,j,k) + u(i-2,j,k))/(12.d0*dx)
    else if (d .eq. P_STENCIL_LEFT) then
      Dxu = (-48.d0*u(i,j,k) + 59.d0*u(i+1,j,k) - 8.d0*u(i+2,j,k) - 3.0d0*u(i+3,j,k))/(34.d0*dx)
    else if (d .eq. P_STENCIL_CENTER_LEFT) then
      Dxu = (u(i+1,j,k) - u(i-1,j,k))/(2.d0*dx)
    else if (d .eq. P_STENCIL_CENTER_LEFT1) then
      Dxu = (8.d0*u(i-2,j,k) - 59.d0*u(i-1,j,k) + 59.d0*u(i+1,j,k) - 8.d0*u(i+2,j,k))/(86.d0*dx)
    else if (d .eq. P_STENCIL_CENTER_LEFT2) then
      Dxu = (3.d0*u(i-3,j,k) - 59.d0*u(i-1,j,k) + 64.d0*u(i+1,j,k) - 8.d0*u(i+2,j,k))/(98.d0*dx)
    else if (d .eq. P_STENCIL_CENTER_RIGHT2) then
      Dxu = (-3.d0*u(i+3,j,k) + 59.d0*u(i+1,j,k) - 64.d0*u(i-1,j,k) + 8.d0*u(i-2,j,k))/(98.d0*dx)      
    else if (d .eq. P_STENCIL_CENTER_RIGHT1) then
      Dxu = (-8.d0*u(i+2,j,k) + 59.d0*u(i+1,j,k) - 59.d0*u(i-1,j,k) + 8.d0*u(i-2,j,k))/(86.d0*dx)
    else if (d .eq. P_STENCIL_CENTER_RIGHT) then
      Dxu = (u(i+1,j,k) - u(i-1,j,k))/(2.d0*dx)    
    else if (d .eq. P_STENCIL_RIGHT) then
      Dxu = (48.d0*u(i,j,k) - 59.d0*u(i-1,j,k) + 8.d0*u(i-2,j,k) + 3.0d0*u(i-3,j,k))/(34.d0*dx)
    else if (d .eq. P_STENCIL_MASK) then
      Dxu = 0.0          
    else if (d .eq. P_STENCIL_UNDEF) then
      Dxu = 0.0
    else
      write(0,*)'d42: Illegal x derivative type ',d
      write(0,*)'d42: imask                     ',imask
      write(0,*)'     at computational location ',&
                 wpt(H_X),wpt(H_Y),wpt(H_Z)
      write(0,*)'     at physical location ',&
                 wpt(H_XPHYS),wpt(H_YPHYS),wpt(H_ZPHYS)
      Dxu = 0.0
      STOP
    end if
   
    !-----------------------------------------------------------------
    ! y-derivatives
    !-----------------------------------------------------------------

    d = mod((imask/100), 100)
    if (d .eq. P_STENCIL_CENTER) then
      Dyu = (-u(i,j+2,k) + 8.d0*u(i,j+1,k) - 8.d0*u(i,j-1,k) + u(i,j-2,k))/(12.d0*dy)
    else if (d .eq. P_STENCIL_LEFT) then
      Dyu = (-48.d0*u(i,j,k) + 59.d0*u(i,j+1,k) - 8.d0*u(i,j+2,k) - 3.0d0*u(i,j+3,k))/(34.d0*dy)
    else if (d .eq. P_STENCIL_CENTER_LEFT) then
      Dyu = (u(i,j+1,k) - u(i,j-1,k))/(2.d0*dy)
    else if (d .eq. P_STENCIL_CENTER_LEFT1) then
      Dyu = (8.d0*u(i,j-2,k) - 59.d0*u(i,j-1,k) + 59.d0*u(i,j+1,k) - 8.d0*u(i,j+2,k))/(86.d0*dy)
    else if (d .eq. P_STENCIL_CENTER_LEFT2) then
      Dyu = (3.d0*u(i,j-3,k) - 59.d0*u(i,j-1,k) + 64.d0*u(i,j+1,k) - 8.d0*u(i,j+2,k))/(98.d0*dy)
    else if (d .eq. P_STENCIL_CENTER_RIGHT2) then
      Dyu = (-3.d0*u(i,j+3,k) + 59.d0*u(i,j+1,k) - 64.d0*u(i,j-1,k) + 8.d0*u(i,j-2,k))/(98.d0*dy)      
    else if (d .eq. P_STENCIL_CENTER_RIGHT1) then
      Dyu = (-8.d0*u(i,j+2,k) + 59.d0*u(i,j+1,k) - 59.d0*u(i,j-1,k) + 8.d0*u(i,j-2,k))/(86.d0*dy)
    else if (d .eq. P_STENCIL_CENTER_RIGHT) then
      Dyu = (u(i,j+1,k) - u(i,j-1,k))/(2.d0*dy)    
    else if (d .eq. P_STENCIL_RIGHT) then
      Dyu = (48.d0*u(i,j,k) - 59.d0*u(i,j-1,k) + 8.d0*u(i,j-2,k) + 3.0d0*u(i,j-3,k))/(34.d0*dy)
    else if (d .eq. P_STENCIL_MASK) then
      Dxu = 0.0            
    else if (d .eq. P_STENCIL_UNDEF) then
      Dyu = 0.0
    else
      write(0,*)'d42: Illegal y derivative type ',d
      write(0,*)'d42: imask                     ',imask
      write(0,*)'     at computational location ',&
                 wpt(H_X),wpt(H_Y),wpt(H_Z)
      write(0,*)'     at physical location ',&
                 wpt(H_XPHYS),wpt(H_YPHYS),wpt(H_ZPHYS)
      Dyu = 0.0
      STOP
    end if
    
    !-----------------------------------------------------------------
    ! z-derivatives
    !-----------------------------------------------------------------

    d = mod(imask,100)
    if (d .eq. P_STENCIL_CENTER) then
      Dzu = (-u(i,j,k+2) + 8.d0*u(i,j,k+1) - 8.d0*u(i,j,k-1) + u(i,j,k-2))/(12.d0*dz)
    else if (d .eq. P_STENCIL_LEFT) then
      Dzu = (-48.d0*u(i,j,k) + 59.d0*u(i,j,k+1) - 8.d0*u(i,j,k+2) - 3.0d0*u(i,j,k+3))/(34.d0*dz)
    else if (d .eq. P_STENCIL_CENTER_LEFT) then
      Dzu = (u(i,j,k+1) - u(i,j,k-1))/(2.d0*dz)
    else if (d .eq. P_STENCIL_CENTER_LEFT1) then
      Dzu = (8.d0*u(i,j,k-2) - 59.d0*u(i,j,k-1) + 59.d0*u(i,j,k+1) - 8.d0*u(i,j,k+2))/(86.d0*dz)
    else if (d .eq. P_STENCIL_CENTER_LEFT2) then
      Dzu = (3.d0*u(i,j,k-3) - 59.d0*u(i,j,k-1) + 64.d0*u(i,j,k+1) - 8.d0*u(i,j,k+2))/(98.d0*dz)
    else if (d .eq. P_STENCIL_CENTER_RIGHT2) then
      Dzu = (-3.d0*u(i,j,k+3) + 59.d0*u(i,j,k+1) - 64.d0*u(i,j,k-1) + 8.d0*u(i,j,k-2))/(98.d0*dz)      
    else if (d .eq. P_STENCIL_CENTER_RIGHT1) then
      Dzu = (-8.d0*u(i,j,k+2) + 59.d0*u(i,j,k+1) - 59.d0*u(i,j,k-1) + 8.d0*u(i,j,k-2))/(86.d0*dz)
    else if (d .eq. P_STENCIL_CENTER_RIGHT) then
      Dzu = (u(i,j,k+1) - u(i,j,k-1))/(2.d0*dz)    
    else if (d .eq. P_STENCIL_RIGHT) then
      Dzu = (48.d0*u(i,j,k) - 59.d0*u(i,j,k-1) + 8.d0*u(i,j,k-2) + 3.0d0*u(i,j,k-3))/(34.d0*dz)
    else if (d .eq. P_STENCIL_MASK) then
      Dxu = 0.0            
    else if (d .eq. P_STENCIL_UNDEF) then
      Dzu = 0.0
    else
      write(0,*)'d42: Illegal z derivative type ',d
      write(0,*)'d42: imask                     ',imask
      write(0,*)'     at computational location ',&
                 wpt(H_X),wpt(H_Y),wpt(H_Z)
      write(0,*)'     at physical location ',&
                 wpt(H_XPHYS),wpt(H_YPHYS),wpt(H_ZPHYS)
      Dzu = 0.0
      STOP
    end if
    

    !
    ! Check to see if we are using fish-eye coords, and transform
    ! the derivatives is necessary.
    !
    if (P_ENABLE_ALT_COORDS .eq. 1 .and. &
    nint(par(P_ALT_COORD_TYPE)) .ne. 0) then
       if (nint(par(P_ALT_COORD_TYPE)) .eq. 1) then
          call  pcjac_pt(jac, wpt, par)
          dxc = Dxu
          dyc = Dyu
          dzc = Dzu
          Dxu = jac(1,1)*dxc + jac(1,2)*dyc + jac(1,3)*dzc
          Dyu = jac(2,1)*dxc + jac(2,2)*dyc + jac(2,3)*dzc
          Dzu = jac(3,1)*dxc + jac(3,2)*dyc + jac(3,3)*dzc
       else 
          write(0,*)'D42_PT: unknown alt_coord_type ',par(P_ALT_COORD_TYPE)
          stop
       end if
    end if


    return
  end subroutine D42_pt

END MODULE HYPER_DERIVS_42_PT
