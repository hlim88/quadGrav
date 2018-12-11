#include "cctk.h"

module exact_u
  use params
#if defined GH_MHD || defined GH_FLOWER
  use gh_exact_quantities
#else
  use exact_quantities
#endif
!  use change_vars_mod

  implicit none
CONTAINS

  subroutine exactu(par,x,y,z,uexact_pt,w_pt)
    use GF

    real(kind=8), dimension(:) :: par
    real(kind=8), dimension(NU) :: uexact_pt
    real(kind=8), dimension(NW) :: w_pt
    real(kind=8) :: x,y,z    
    
  end subroutine exactu
     
end module exact_u
