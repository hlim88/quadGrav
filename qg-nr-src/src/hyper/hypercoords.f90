!-------------------------------------------------------------------------
!
!   $Id$
!
!   This module contains a couple routines required for the
!   pincushion coordinates.
!
!-------------------------------------------------------------------------

#include "cctk.h"

MODULE HYPERCOORDS
  use PARAMS
  use GF
  implicit none

  CONTAINS

  !------------------------------------------------------------
  !
  !  Maps a 1d, uniform coordinate array to a 3d array, as used 
  !  by Cactus.  Added for compatibility with Steves code.
  !
  !------------------------------------------------------------
  subroutine mapCoords1dTo3d(x3d, y3d, z3d, x1d, y1d, z1d, nx, ny, nz)
    implicit none

    CCTK_INT                         ::  nx, ny, nz
    CCTK_REAL, dimension(nx, ny, nz) ::  x3d, y3d, z3d
    CCTK_REAL, dimension(nx)         ::  x1d
    CCTK_REAL, dimension(ny)         ::  y1d
    CCTK_REAL, dimension(nz)         ::  z1d


    ! local vars
    CCTK_INT                         ::  i, j, k
    logical, parameter               ::  ltrace = .false.

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
      x3d(i,j,k) = x1d(i)
    end do
    end do
    end do

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
      y3d(i,j,k) = y1d(j)
    end do
    end do
    end do

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
      z3d(i,j,k) = z1d(k)
    end do
    end do
    if (ltrace) then
      write(0,*) ' z3d(1,1,k), z1d(k) = ',z3d(1,1,k), z1d(k)
    end if
    end do

    return
  end subroutine mapCoords1dTo3d

  !------------------------------------------------------------
  !
  !
  !
  !
  !------------------------------------------------------------
  subroutine def_alt_coords(w, par)
    implicit none

    type(gridfunction), dimension(NW)    :: w
    CCTK_REAL, dimension(NPAR)           :: par

    ! local vars
    CCTK_REAL, dimension(NW)             :: wpt
    CCTK_REAL, dimension(:,:,:), pointer :: xp, yp, zp, xc, yc, zc, rc
    CCTK_INT                             :: i, j, k, nx, ny, nz
    CCTK_REAL                            :: s

    ! Check that coordinate parameters are valid.
    ! Do it here once, instead of nx*ny*nz times down below.
    if (nint(par(P_ALT_COORD_TYPE)) == 1) then
      s = par(P_PC_COORD_TRANS_WIDTH)
      if (s == 0.0) then
        write(0,*)'DEF_ALT_COORDS >>> the coord_trans_width is zero ',s
        stop
      end if
    end if

    xc => w(H_X)%d
    yc => w(H_Y)%d
    zc => w(H_Z)%d
    rc => w(H_R)%d
    xp => w(H_XPHYS)%d
    yp => w(H_YPHYS)%d
    zp => w(H_ZPHYS)%d

    nx = nint(par(P_NX))
    ny = nint(par(P_NY))
    nz = nint(par(P_NZ))

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx

      if (nint(par(P_ALT_COORD_TYPE)) == 0) then
        ! Physical coordinates are the computational (Reg. Cart. grid)
        xp(i,j,k) = xc(i,j,k)
        yp(i,j,k) = yc(i,j,k)
        zp(i,j,k) = zc(i,j,k)
      else if (nint(par(P_ALT_COORD_TYPE)) == 1) then
        ! Pincushion coordinates

        wpt(H_X) = xc(i,j,k)
        wpt(H_Y) = yc(i,j,k)
        wpt(H_Z) = zc(i,j,k)
        wpt(H_R) = rc(i,j,k)

        call pccoord_pt(wpt, par)

        xp(i,j,k) = wpt(H_XPHYS)
        yp(i,j,k) = wpt(H_YPHYS)
        zp(i,j,k) = wpt(H_ZPHYS)
        !print *,'i,j,k ',xc(i,j,k),xp(i,j,k)

      else
        write(0,*)'Unknown alternative coordinate type ', &
                  nint(par(P_ALT_COORD_TYPE))
        stop
      end if

    end do
    end do
    end do
    
    return
  end subroutine def_alt_coords

  !------------------------------------------------------------
  !
  !
  !
  !
  !------------------------------------------------------------
  subroutine pccoord_pt(w,par)
    implicit none
    CCTK_REAL, dimension(NW)   :: w
    CCTK_REAL, dimension(NPAR) :: par

    ! local vars
    CCTK_REAL  :: xc, yc, zc, rc, rp
    CCTK_REAL  :: a, s, r0

    s = par(P_PC_COORD_TRANS_WIDTH)
    r0 = par(P_PC_COORD_TRANS_RAD)
    a = par(P_PC_COORD_SCALE)

    xc = w(H_X)
    yc = w(H_Y)
    zc = w(H_Z)
    rc = w(H_R)

    if (rc < 1.0e-11) then
      ! the origin is special
      w(H_XPHYS) = 0.0
      w(H_YPHYS) = 0.0
      w(H_ZPHYS) = 0.0
    else
      rp = a*rc + (1 - a)*s/(2*tanh(r0/s)) &
                 *( log(cosh((rc + r0)/s)) - log(cosh((rc - r0)/s)))

      w(H_XPHYS) = rp/rc * xc
      w(H_YPHYS) = rp/rc * yc
      w(H_ZPHYS) = rp/rc * zc
    end if


    return
  end subroutine pccoord_pt

  !------------------------------------------------------------
  !
  ! SUBROUTINE PCJAC
  !
  !   Calculate the Jacobian of coord. transformation to correctly
  !   transform derivatives
  !
  !------------------------------------------------------------
  subroutine pcjac_pt(jac, wpt, par)
    implicit none

    CCTK_REAL, dimension(3,3)              :: jac
    CCTK_REAL, dimension(NW)               :: wpt
    CCTK_REAL, dimension(NPAR)             :: par

    ! local vars
    CCTK_REAL   :: detjac, idetjac, ijac(3,3)
    CCTK_REAL   :: xc, yc, zc, r0, s, a

    ! Maple vars
    CCTK_REAL  t2, t3, t4, t5
    CCTK_REAL  t7, t8, t11, t12, t14
    CCTK_REAL  t15, t17, t18, t19, t21
    CCTK_REAL  t22, t23, t24, t25, t28
    CCTK_REAL  t29, t30, t31, t32, t33
    CCTK_REAL  t36, t37, t38, t44, t48
    CCTK_REAL  t52, t60, t64, t72


      xc = wpt(H_X)
      yc = wpt(H_Y)
      zc = wpt(H_Z)
      r0 = par(P_PC_COORD_TRANS_RAD)
      s  = par(P_PC_COORD_TRANS_WIDTH)
      a  = par(P_PC_COORD_SCALE)

      t2 = (1-a)*s
      t3 = xc**2
      t4 = yc**2
      t5 = zc**2
      t7 = sqrt(t3+t4+t5)

      if (t7 .ne. 0.0) then
         t8 = t7**2
         t11 = t2/t8/t7
         t12 = 1/s
         t14 = tanh(r0*t12)
         t15 = 1/t14
         t17 = (t7+r0)*t12
         t18 = cosh(t17)
         t19 = log(t18)
         t21 = (t7-r0)*t12
         t22 = cosh(t21)
         t23 = log(t22)
         t24 = t19-t23
         t25 = t15*t24
         t28 = 1/t7
         t29 = t28*t15
         t30 = sinh(t17)
         t31 = t30*t28
         t32 = xc*t12
         t33 = 1/t18
         t36 = sinh(t21)
         t37 = t36*t28
         t38 = 1/t22
         t44 = -t11*t25*xc+t2*t29*(t31*t32*t33-t37*t32*t38)
         t48 = t2*t29*t24/2
         t52 = yc*t12
         t60 = -t11*t25*yc+t2*t29*(t31*t52*t33-t37*t52*t38)
         t64 = zc*t12
         t72 = -t11*t25*zc+t2*t29*(t31*t64*t33-t37*t64*t38)
         ijac(1,1) = t44*xc/2+a+t48
         ijac(1,2) = t60*xc/2
         ijac(1,3) = t72*xc/2
         ijac(2,1) = t44*yc/2
         ijac(2,2) = t60*yc/2+a+t48
         ijac(2,3) = t72*yc/2
         ijac(3,1) = t44*zc/2
         ijac(3,2) = t60*zc/2
         ijac(3,3) = t72*zc/2+a+t48
         detjac = ijac(1,1)*ijac(2,2)*ijac(3,3)&
                     -ijac(1,1)*ijac(2,3)*ijac(3,2) &
                     -ijac(2,1)*ijac(1,2)*ijac(3,3) &
                     +ijac(2,1)*ijac(1,3)*ijac(3,2) &
                     +ijac(3,1)*ijac(1,2)*ijac(2,3) &
                     -ijac(3,1)*ijac(1,3)*ijac(2,2)
         if (detjac .eq. 0.0d0) then
            write(0,*)'PINCUSHION: detjac = 0'
            stop
         end if
         idetjac = 1.0/detjac
         jac(1,1) = idetjac*(ijac(2,2)*ijac(3,3)-ijac(2,3)*ijac(3,2))
         jac(1,2) = idetjac*(-ijac(1,2)*ijac(3,3)+ijac(1,3)*ijac(3,2))
         jac(1,3) = idetjac*(ijac(1,2)*ijac(2,3)-ijac(1,3)*ijac(2,2))
         jac(2,1) = idetjac*(-ijac(2,1)*ijac(3,3)+ijac(2,3)*ijac(3,1))
         jac(2,2) = idetjac*(ijac(1,1)*ijac(3,3)-ijac(1,3)*ijac(3,1))
         jac(2,3) = idetjac*(-ijac(1,1)*ijac(2,3)+ijac(1,3)*ijac(2,1))
         jac(3,1) = idetjac*(ijac(2,1)*ijac(3,2)-ijac(2,2)*ijac(3,1))
         jac(3,2) = idetjac*(-ijac(1,1)*ijac(3,2)+ijac(1,2)*ijac(3,1))
         jac(3,3) = idetjac*(ijac(1,1)*ijac(2,2)-ijac(1,2)*ijac(2,1))

      else
         jac(1,1) = 1.0
         jac(1,2) = 0.0
         jac(1,3) = 0.0
         jac(2,1) = 0.0
         jac(2,2) = 1.0
         jac(2,3) = 0.0
         jac(3,1) = 0.0
         jac(3,2) = 0.0
         jac(3,3) = 1.0
      end if

    return
  end subroutine pcjac_pt

END MODULE HYPERCOORDS


