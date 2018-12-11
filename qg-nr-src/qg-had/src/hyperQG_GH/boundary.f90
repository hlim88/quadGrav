!-------------------------------------------------------------------------
!
! $Id: boundary.f90,v 1.2 2007-03-27 22:25:26 dneilsen Exp $
!
!-------------------------------------------------------------------------

#include "cctk.h"

MODULE HYPER_BOUNDARY
  use params
  use GF
  use m_char_boundary
  implicit none

CONTAINS

!---------------------------------------------------------------------------
!
!  Some comments: what is called u in this subroutine is actually 
!  u or u_dot, depending on whether we are giving boundary conditions 
!  to the rhs or the vars.
!
!____________________________________________________________________________

  subroutine hyperboundary(u, dxu, dyu, dzu, udot, v, w, imask, par)
    type(gridfunction), dimension(NU_G), intent(inout) :: u, dxu,dyu,dzu,udot
    type(gridfunction), dimension(NV_G) :: v
    type(gridfunction), dimension(NW) :: w
    CCTK_REAL, dimension(:)           :: par
    CCTK_INT, dimension(:,:,:)        :: imask
    ! local vars
    CCTK_REAL, dimension(:,:,:), pointer :: xx, yy, zz
    CCTK_REAL, DIMENSION(3) :: n  ! outward unit normal
    CCTK_INT :: imin,imax,jmin,jmax,kmin,kmax, IBBOX1, IBBOX2, IBBOX3,&
         IBBOX4, IBBOX5, IBBOX6, nx, ny, nz, k, i, j, l, NB
    CCTK_INT :: idx, idy, idz, mask_stencil
    logical bitantz
    logical double_equal

    logical, parameter :: ltrace = .false.


    if (ltrace) write(*,*)' @@@@@  hyperboundary '


    call char_boundary(u, dxu, dyu, dzu, udot, v, w, imask, par)

    return
  end subroutine hyperboundary

end module
