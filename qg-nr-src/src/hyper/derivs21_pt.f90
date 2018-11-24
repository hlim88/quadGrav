!-----------------------------------------------------------------------
!
!    $Id$
!
!-----------------------------------------------------------------------

#include "cctk.h"

MODULE HYPER_DERIVS_21_PT
use params
use hypercoords
implicit none

CONTAINS

!-----------------------------------------------------------------------
!
!    SUBROUTINE D21NoMask_pt:
!
!    implement Strand derivative.  2nd order interior, 1st order boundary
!
!-----------------------------------------------------------------------
  subroutine D21NoMask_pt(Dxu, Dyu, Dzu, u, wpt, i, j, k, par)
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
      Dxu = (u(i+1,j,k)-u(i,j,k))/dx
    else if (i .gt. 1 .and. i .lt. nx) then
      Dxu = (u(i+1,j,k)-u(i-1,j,k))/(2*dx)
    else if (i .eq. nx) then
      Dxu = (u(nx,j,k)-u(nx-1,j,k))/dx
    else
      write(0,*)'D21NOMASK_PT>> i out of range in x-deriv. i=',i
    end if
    !-------------------------------------------------
    !  y derivs
    !-------------------------------------------------
    if (j .eq. 1) then
      Dyu = (u(i,j+1,k)-u(i,j,k))/dy
    else if (j .gt. 1 .and. j .lt. ny) then
      Dyu = (u(i,j+1,k)-u(i,j-1,k))/(2*dy)
    else if (j .eq. ny) then
      Dyu = (u(i,ny,k)-u(i,ny-1,k))/dy
    else
      write(0,*)'D21NOMASK_PT>> j out of range in y-deriv. j=',j
    end if

    !-------------------------------------------------
    !  z derivs
    !-------------------------------------------------
    if (k .eq. 1) then
      Dzu = (u(i,j,2)-u(i,j,1))/dz
    else if (k .gt. 1 .and. k .lt. nz) then
      Dzu = (u(i,j,k+1)-u(i,j,k-1))/(2*dz)
    else if (k .eq. nz) then
      Dzu = (u(i,j,nz)-u(i,j,nz-1))/dz
    else
      write(0,*)'D21NOMASK_PT>> k out of range in z-deriv. k=',k
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
          write(0,*)'D21NoMask_PT: unknown alt_coord_type ',&
                     par(P_ALT_COORD_TYPE)
          stop
       end if
    end if


    return
  end subroutine D21NoMask_pt


!-----------------------------------------------------------------------
!
!    D21_pt
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
  subroutine D21_pt(Dxu, Dyu,  Dzu,  u, wpt,  &
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
      write(0,*)'D21_PT >> i out of range: i/nx = ',i,nx
    end if
    if (j .lt. 1 .or. j .gt. ny) then
      write(0,*)'D21_PT >> j out of range: y/ny = ',j,ny
    end if
    if (k .lt. 1 .or. k .gt. nz) then
      write(0,*)'D21_PT >> k out of range: k/nz = ',k,nz
    end if

    !-----------------------------------------------------------------
    ! x-derivatives
    !-----------------------------------------------------------------

    d = imask / 10000
    if (d .eq. P_STENCIL_CENTER) then
      Dxu = (u(i+1,j,k) - u(i-1,j,k)) / (2.0*dx)
    else if (d .eq. P_STENCIL_CENTER_LEFT .or. &
             d .eq. P_STENCIL_CENTER_RIGHT) then
      Dxu = (u(i+1,j,k) - u(i-1,j,k)) / (2.0*dx)
    else if (d .eq. P_STENCIL_LEFT) then
      Dxu = (u(i+1,j,k) - u(i,j,k)) / dx
    else if (d .eq. P_STENCIL_RIGHT) then
      Dxu = (u(i,j,k) - u(i-1,j,k)) / dx
    else if (d .eq. P_STENCIL_MASK) then
      Dxu = 0.0
   else if (d .eq. P_STENCIL_EDGE_P) then
     Dxu =  (2.0*u(i+1,j,k) - u(i,j,k) - u(i-1,j,k))/(3.0*dx)
   else if (d .eq. P_STENCIL_EDGE_M) then
     Dxu = (-2.0*u(i-1,j,k) + u(i,j,k) + u(i+1,j,k))/(3.0*dx)
   else if (d .eq. P_STENCIL_VERTEX_P) then
     Dxu = 4.0/(7.0*dx)&
                 *(u(i+1,j,k) - 0.25*u(i,j,k) - 0.75*u(i-1,j,k))
   else if (d .eq. P_STENCIL_VERTEX_M) then
     Dxu = 4.0/(7.0*dx)&
                 *(-u(i-1,j,k) + 0.25*u(i,j,k) + 0.75*u(i+1,j,k))
    else if (d .eq. P_STENCIL_UNDEF) then
      Dxu = 0.0
    else
      write(0,*)'d21: Illegal x derivative type ',d
      write(0,*)'d21: imask                     ',imask
      write(0,*)'     at computational location ',&
                 wpt(H_X),wpt(H_Y),wpt(H_Z)
      write(0,*)'     at physical location ',&
                 wpt(H_XPHYS),wpt(H_YPHYS),wpt(H_ZPHYS)
      Dxu = 0.0
      stop
    end if

    !-----------------------------------------------------------------
    ! y-derivatives
    !-----------------------------------------------------------------

    d = mod((imask/100), 100)
    if (d .eq. P_STENCIL_CENTER) then
      Dyu = (u(i,j+1,k) - u(i,j-1,k)) / (2.0*dy)
    else if (d .eq. P_STENCIL_LEFT) then
      Dyu = (u(i,j+1,k) - u(i,j,k)) / dy
    else if (d .eq. P_STENCIL_CENTER_LEFT .or. &
            d .eq. P_STENCIL_CENTER_RIGHT) then
      Dyu = (u(i,j+1,k) - u(i,j-1,k)) / (2.0*dy)
    else if (d .eq. P_STENCIL_RIGHT) then
      Dyu = (u(i,j,k) - u(i,j-1,k)) / dy
    else if (d .eq. P_STENCIL_MASK) then
      Dyu = 0.0
    else if (d .eq. P_STENCIL_EDGE_P) then
      Dyu =  (2.0*u(i,j+1,k) - u(i,j,k) - u(i,j-1,k))/(3.0*dy)
    else if (d .eq. P_STENCIL_EDGE_M) then
      Dyu = (-2.0*u(i,j-1,k) + u(i,j,k) + u(i,j+1,k))/(3.0*dy)
    else if (d .eq. P_STENCIL_VERTEX_P) then
      Dyu = 4.0/(7.0*dy)&
                  *(u(i,j+1,k) - 0.25*u(i,j,k) - 0.75*u(i,j-1,k))
    else if (d .eq. P_STENCIL_VERTEX_M) then
      Dyu = 4.0/(7.0*dy)&
                  *(-u(i,j-1,k) + 0.25*u(i,j,k) + 0.75*u(i,j+1,k))
    else if (d .eq. P_STENCIL_UNDEF) then
      if (ltrace2) then
        if ((bbox3==0 .and. j < gz_ymin) .or. &
            (bbox4==0 .and. j > gz_ymax)) then
          write(0,*)'d21: OK---Undefined y derivative type in ghostzone'
        else
          write(0,*)'d21: ### Undefined y derivative outside of ghost zone!'
          write(0,*)'d21: imask                          ',imask
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
      Dyu = 0.0
    else
      write(0,*)'d21: Illegal y derivative type   ',d
      write(0,*)'d21: imask                       ',imask
      write(0,*)'     at computational location ',&
                 wpt(H_X),wpt(H_Y),wpt(H_Z)
      write(0,*)'     at physical location ',&
                 wpt(H_XPHYS),wpt(H_YPHYS),wpt(H_ZPHYS)
      Dyu = 0.0
      stop
    end if

    !-----------------------------------------------------------------
    ! z-derivatives
    !-----------------------------------------------------------------
 
    d = mod(imask,100)
    if (d .eq. P_STENCIL_CENTER) then
      Dzu = (u(i,j,k+1) - u(i,j,k-1)) / (2.0*dz)
    else if (d .eq. P_STENCIL_LEFT) then
      Dzu = (u(i,j,k+1) - u(i,j,k)) / dz
    else if (d .eq. P_STENCIL_CENTER_LEFT .or. &
             d .eq. P_STENCIL_CENTER_RIGHT) then
      Dzu = (u(i,j,k+1) - u(i,j,k-1)) / (2.0*dz)
    else if (d .eq. P_STENCIL_RIGHT) then
      Dzu = (u(i,j,k) - u(i,j,k-1)) / dz
    else if (d .eq. P_STENCIL_MASK) then
      Dzu = 0.0
    else if (d .eq. P_STENCIL_EDGE_P) then
      Dzu =  (2.0*u(i,j,k+1) - u(i,j,k) - u(i,j,k-1))/(3.0*dz)
    else if (d .eq. P_STENCIL_EDGE_M) then
      Dzu = (-2.0*u(i,j,k-1) + u(i,j,k) + u(i,j,k+1))/(3.0*dz)
    else if (d .eq. P_STENCIL_VERTEX_P) then
      Dzu = 4.0/(7.0*dz)&
                  *(u(i,j,k+1) - 0.25*u(i,j,k) - 0.75*u(i,j,k-1))
    else if (d .eq. P_STENCIL_VERTEX_M) then
      Dzu = 4.0/(7.0*dz)&
                  *(-u(i,j,k-1) + 0.25*u(i,j,k) + 0.75*u(i,j,k+1))
    else if (d .eq. P_STENCIL_UNDEF) then
      if (ltrace2) then
        if ((bbox5==0 .and. k < gz_zmin) .or. &
            (bbox6==0 .and. k > gz_zmax)) then
          write(0,*)'d21: OK---Undefined z derivative type in ghostzone'
        else
          write(0,*)'d21: ### Undefined z derivative outside of ghost zone!'
          write(0,*)'d21: imask                          ',imask
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
      Dzu = 0.0
    else
      write(0,*)'d21: Illegal z derivative type   ',d
      write(0,*)'d21: imask                       ',imask
      write(0,*)'     at computational location ',&
                 wpt(H_X),wpt(H_Y),wpt(H_Z)
      write(0,*)'     at physical location ',&
                 wpt(H_XPHYS),wpt(H_YPHYS),wpt(H_ZPHYS)
      Dzu = 0.0
      stop
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
          write(0,*)'D21_PT: unknown alt_coord_type ',par(P_ALT_COORD_TYPE)
          stop
       end if
    end if


    return
  end subroutine D21_pt

END MODULE HYPER_DERIVS_21_PT
