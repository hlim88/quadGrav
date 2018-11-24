
! Variable timestep modification by Dominic Marcello 2009-08-03
!-------------------------------------------------------------------------
!
!  $Id$
!
!-------------------------------------------------------------------------

#include "cctk.h"
#ifdef DYNAMIC_CFL

subroutine HyperComputeCFL(cs, u, v, w, par)
  use cfl_mine
  use params
  use GF
  implicit none

  type(gridfunction), dimension(NU):: u
  type(gridfunction), dimension(NV):: v
  type(gridfunction), dimension(NW):: w
  CCTK_REAL, intent(inout) :: cs
  CCTK_REAL, dimension(NPAR)       :: par

  ! local vars
  logical, parameter :: ltrace = .false.
  CCTK_INT :: nx,ny,nz

  nx = nint(par(P_NX))
  ny = nint(par(P_NY))
  nz = nint(par(P_NZ))

  call calccfl (cs, u, v, w, par)


  return
end subroutine HyperComputeCFL

! -------------end modification by Dominic Marcello 2009-08-03

#endif
