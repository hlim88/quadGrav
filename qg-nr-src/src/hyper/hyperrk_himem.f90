!-------------------------------------------------------------------------
!
!  $Id$
!
!-------------------------------------------------------------------------

#include "cctk.h"

subroutine HyperEvolveRK(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
     urk2, urk3, par)
  use params
  use GF
  implicit none

  type(gridfunction), dimension(NU):: u0,u2, &
                                      urk1,urk2,urk3,&
                                      dxu,dyu,dzu

  type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
  type(gridfunction), dimension(NW):: w
  CCTK_REAL, dimension(NPAR)           :: par

  ! local vars
  CCTK_INT  :: rk_type

  interface
    subroutine RK3_himem(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
         urk2, urk3, par)
      use params
      use UTILEQS
      use GF
      use HYPER_RHS
      use HYPER_DISSIPATION
      use HYPER_BOUNDARY
      use MOD_MASK
      use hypercoords
      implicit none
      type(gridfunction), dimension(NU):: u0,u2, urk1,urk2,urk3,&
                                          dxu,dyu,dzu
      type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
      type(gridfunction), dimension(NW):: w
      CCTK_REAL, dimension(NPAR)           :: par
    end subroutine RK3_himem
    subroutine TVDRK3_himem(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
         urk2, urk3, par)
      use params
      use UTILEQS
      use GF
      use HYPER_RHS
      use HYPER_DISSIPATION
      use HYPER_BOUNDARY
      use MOD_MASK
      use hypercoords
      implicit none
      type(gridfunction), dimension(NU):: u0,u2, urk1,urk2,urk3,&
                                          dxu,dyu,dzu
      type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
      type(gridfunction), dimension(NW):: w
      CCTK_REAL, dimension(NPAR)           :: par
    end subroutine TVDRK3_himem
    subroutine TVDRK3_FV_himem(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
         urk2, urk3, par)
      use params
      use UTILEQS
      use GF
      use HYPER_RHS
      use HYPER_DISSIPATION
      use HYPER_BOUNDARY
      use MOD_MASK
      use hypercoords
      implicit none
      type(gridfunction), dimension(NU):: u0,u2, urk1,urk2,urk3,&
                                          dxu,dyu,dzu
      type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
      type(gridfunction), dimension(NW):: w
      CCTK_REAL, dimension(NPAR)           :: par
    end subroutine TVDRK3_FV_himem

  end interface

  rk_type = nint(par(P_RUNGE_KUTTA_TYPE))

  if (rk_type .eq. 11) then
    call  RK3_himem(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
         urk2, urk3, par)
  else if (rk_type .eq. 12) then
    call  TVDRK3_himem(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
         urk2, urk3, par)
  else if (rk_type .eq. 13) then
    call  TVDRK3_FV_himem(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
         urk2, urk3, par)
  else
    write(0,*)'HyperEvolveRK:  unknown RK type.  runge_kutta_type = ',rk_type
    write(0,*)'    This executable compiled for high memory'
    write(0,*)'    Only runge_kutta_type = 11, 12, 13 available'
    write(0,*)'    For runge_kutta_type = 1,2,3, recompile without the HAD_HIMEM'
    write(0,*)'    environment variable.'
    stop
  end if

  call HyperAnalysis(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, par)

 

  return
end subroutine HyperEvolveRK
