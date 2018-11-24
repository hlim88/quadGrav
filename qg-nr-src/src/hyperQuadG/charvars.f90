#include "cctk.h"
#if defined GR_MHD || defined GH_MHD
module bssn_charvars
#else
module charvars
#endif
  use params
  implicit none
  private
#if defined GR_MHD || defined GH_MHD
  public bssn_prim2char, bssn_char2prim
#else
  public prim2char, char2prim
#endif
  
contains
  
  !#########################################################################!
#if defined GR_MHD || defined GH_MHD
  subroutine bssn_prim2char(par,v,w,direction,u,w_in,w_out,w_0,T,uevolved,uexact,dxu,dyu,dzu,sources) 
#else
  subroutine prim2char(par,v,w,direction,u,w_in,w_out,w_0,T,uevolved,uexact,dxu,dyu,dzu,sources)
#endif
    implicit none

    CCTK_REAL, dimension(:), intent(in) :: par
    CCTK_REAL, dimension(:), intent(in) :: v, w
    character(len=*), intent(in):: direction   
    CCTK_REAL, dimension(:), intent(in) :: u,uevolved,uexact,dxu,dyu,dzu
    CCTK_REAL, dimension(:), intent(in) :: w_in, w_out,sources
    CCTK_REAL, dimension(:), intent(in) :: w_0
    CCTK_REAL, dimension(3,3)            :: T

    return
  end subroutine
  
  !##########################################################################!
#if defined GR_MHD || defined GH_MHD
  subroutine bssn_char2prim(par,v,w,direction,u,w_in,w_out,w_0,T,uevolved,uexact) 
#else
  subroutine char2prim(par,v,w,direction,u,w_in,w_out,w_0,T,uevolved,uexact) 
#endif
    implicit none

    CCTK_REAL, dimension(:), intent(in) :: par
    CCTK_REAL, dimension(:), intent(in) :: w,uevolved,uexact,v
    character(len=*), intent(in):: direction
    CCTK_REAL, dimension(:), intent(in) :: u
    CCTK_REAL, dimension(:), intent(in) :: w_in, w_out
    CCTK_REAL, dimension(:), intent(in) :: w_0
    CCTK_REAL, dimension(3,3)             :: T

    return
  end subroutine

  !#######################################################################!

end module
