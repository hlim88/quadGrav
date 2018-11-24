!--------------------------------------------------------------------------
!$Id: exact.f90,v 1.1 2010-07-17 22:09:42 carlos Exp $
! 
!  This subroutine can be used to supply analytic boundary conditions.
!  It is called by hyperboundary.
!
!--------------------------------------------------------------------------

#include "cctk.h"

module exact_u
  use params
  implicit none
CONTAINS

  subroutine exactu(par,x,y,z,uexact_pt,w_pt)
    use GF

    CCTK_REAL, dimension(:) :: par
    CCTK_REAL, dimension(NU) :: uexact_pt
    CCTK_REAL, dimension(NW) :: w_pt
    CCTK_REAL :: x,y,z

  end subroutine exactu

end module exact_u
