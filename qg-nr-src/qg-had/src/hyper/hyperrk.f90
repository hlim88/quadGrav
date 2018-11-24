!-------------------------------------------------------------------------
!
!  $Id$
!
!-------------------------------------------------------------------------

#include "cctk.h"

subroutine HyperEvolveRK(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, &
                         par)
  use params
  use GF
  implicit none

  type(gridfunction), dimension(NU):: u0,u2, urk1, dxu, dyu, dzu
  type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
  type(gridfunction), dimension(NW):: w
  CCTK_REAL, dimension(NPAR)       :: par

  ! local vars
  CCTK_INT  :: rk_type
  logical, parameter :: ltrace = .false.
  CCTK_INT:: proc_return_myid, myid

  CCTK_REAL:: myl2norm3d
  CCTK_INT :: m,nx,ny,nz

  interface
    subroutine RK3(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, par)
      use params
      use UTILEQS
      use GF
      use HYPER_RHS
      use HYPER_DISSIPATION
      use HYPER_BOUNDARY
      use MOD_MASK
      use hypercoords
      implicit none
      type(gridfunction), dimension(NU):: u0,u2, urk1, dxu, dyu, dzu
      type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
      type(gridfunction), dimension(NW):: w
      CCTK_REAL, dimension(NPAR)           :: par
    end subroutine RK3
    !
    subroutine TVDRK3(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, par)
      use params
      use UTILEQS
      use GF
      use HYPER_RHS
      use HYPER_DISSIPATION
      use HYPER_BOUNDARY
      use MOD_MASK
      use hypercoords
      implicit none
      type(gridfunction), dimension(NU):: u0,u2, urk1, dxu, dyu, dzu
      type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
      type(gridfunction), dimension(NW):: w
      CCTK_REAL, dimension(NPAR)           :: par
    end subroutine TVDRK3
    !
    subroutine TVDRK3_FV(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, par)
      use params
      use UTILEQS
      use GF
      use HYPER_RHS
      use HYPER_DISSIPATION
      use HYPER_BOUNDARY
      use MOD_MASK
      use hypercoords
      implicit none
      type(gridfunction), dimension(NU):: u0,u2, urk1, dxu, dyu, dzu
      type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
      type(gridfunction), dimension(NW):: w
      CCTK_REAL, dimension(NPAR)           :: par
    end subroutine TVDRK3_FV
  end interface

  myid = proc_return_myid()

  !! get rid of
  nx = nint(par(P_NX))
  ny = nint(par(P_NY))
  nz = nint(par(P_NZ))

  if (ltrace) then
     write(0,*)myid, '] >>>  Hyperrk begins'
     do m=1,NU
       write(0,*) myid, '] hyperrk: u0 - 1A ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) myid, '] hyperrk: u2 - 1A ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
     end do
  end if

  rk_type = nint(par(P_RUNGE_KUTTA_TYPE))

  if (rk_type .eq. 1) then
    if (ltrace) write(0,*) myid, '] >>>  Hyperrk: calling RK3'
    call  RK3(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, par)
  else if (rk_type .eq. 2) then
    if (ltrace) write(0,*) myid, '] >>>  Hyperrk: calling TVDRK3'
    call  TVDRK3(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, par)
  else if (rk_type .eq. 3) then
    if (ltrace) write(0,*) myid, '] >>>  Hyperrk: calling TVDRK3_FV'
    call  TVDRK3_FV(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, par)
  else
    write(0,*)'HyperEvolveRK:  unknown RK type = ',rk_type
    write(0,*)'    This executable compiled for low memory'
    write(0,*)'    Only rk_type = 1, 2, 3 available'
    write(0,*)'    For rk_type = 11,12,13, recompile with the HAD_HIMEM'
    write(0,*)'    environment variable.'
    stop
  end if


  if (ltrace) then
     do m=1,NU
       write(0,*) myid, '] hyperrk: u0 - 1B ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) myid, '] hyperrk: u2 - 1B ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
     end do
     write(0,*)myid, '] >>>  Hyperrk done.'
  end if

  return
end subroutine HyperEvolveRK

!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
subroutine HyperEvolveStiffRK(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, &
                              urk1, sr1, sr2, sr3, sr4, par)
  use params
  use GF
  implicit none

  type(gridfunction), dimension(NU):: u0,u2, urk1, dxu, dyu, dzu
  type(gridfunction), dimension(NU):: sr1, sr2, sr3, sr4
  type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
  type(gridfunction), dimension(NW):: w
  CCTK_REAL, dimension(NPAR)       :: par

  ! local vars
  CCTK_INT  :: rk_type
  logical, parameter :: ltrace = .false.
  CCTK_REAL:: myl2norm3d
  CCTK_INT :: m,nx,ny,nz

  interface
    subroutine IMEX_RK3(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, &
                        sr1, sr2, sr3, sr4, par)
      use params
      use UTILEQS
      use GF
      use HYPER_RHS
      use HYPER_DISSIPATION
      use HYPER_BOUNDARY
      use MOD_MASK
      use hypercoords
      implicit none
      type(gridfunction), dimension(NU):: u0,u2, urk1, dxu, dyu, dzu
      type(gridfunction), dimension(NU):: sr1, sr2, sr3, sr4
      type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
      type(gridfunction), dimension(NW):: w
      CCTK_REAL, dimension(NPAR)           :: par
    end subroutine IMEX_RK3
  end interface

  call  IMEX_RK3(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, urk1, &
                 sr1, sr2, sr3, sr4, par)

end subroutine HyperEvolveStiffRK
