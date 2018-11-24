#include "cctk.h"


module HYPER_BOUNDARY
   use params
   use GF
   use m_char_boundary
   implicit none

   contains 

!---------------------------------------------------------------------------
!
!  Some comments: what is called u in this subroutine is actually 
!  u or u_dot, depending on whether we are giving boundary conditions 
!  to the rhs or the vars.
!
!____________________________________________________________________________

  subroutine hyperboundary(u, dxu, dyu, dzu, udot, v, w, imask, par)
    type(gridfunction), dimension(NU), intent(inout) :: u, dxu,dyu,dzu,udot
    type(gridfunction), dimension(NV) :: v
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

    if (nint(par(P_RUNGE_KUTTA_BOUND)) == 2) then
      ! Put post rk3 update boundary conditions here
      ! Certain boundary conditions need to be applied
      ! at certain times, either before or after the rk3 update.
      if (nint(par(P_BOUNDARY_CONDITIONS)) == 8) then
        if (ltrace) write(*,*)' @@@@@ outflow boundary conditions for EM fields '
        ! This calls outflow BCs for the E&M fields
        call outflow_boundary(u, w, imask, par)
      end if

    else if (nint(par(P_RUNGE_KUTTA_BOUND)) == 1) then
      ! Apply pre rk3 update boundary conditions here
    end if


    return
  end subroutine hyperboundary


  !--------------------------------------------------------------------
  !
  !
  !
  !--------------------------------------------------------------------
  subroutine outflow_boundary(u, w, imask, par)
    implicit   none
    type(gridfunction), dimension(NU), intent(inout) :: u
    type(gridfunction), dimension(NW) :: w
    CCTK_REAL, dimension(:)           :: par
    CCTK_INT, dimension(:,:,:)        :: imask

    ! local vars
    CCTK_INT                          :: nx, ny, nz, ng
    
    nx = nint(par(P_NX))  ! local size of arrays
    ny = nint(par(P_NY))
    nz = nint(par(P_NZ))
    ng = nint(par(P_BOUND_WIDTH))

    call apply_outflow(u(G_EX)%d, w(H_chr)%d, nx, ny, nz, ng, par)
    call apply_outflow(u(G_EY)%d, w(H_chr)%d, nx, ny, nz, ng, par)
    call apply_outflow(u(G_EZ)%d, w(H_chr)%d, nx, ny, nz, ng, par)

    call apply_outflow(u(G_BX)%d, w(H_chr)%d, nx, ny, nz, ng, par)
    call apply_outflow(u(G_BY)%d, w(H_chr)%d, nx, ny, nz, ng, par)
    call apply_outflow(u(G_BZ)%d, w(H_chr)%d,nx, ny, nz, ng, par)

    call apply_outflow(u(G_PHI_EM)%d, w(H_chr)%d, nx, ny, nz, ng, par)
    call apply_outflow(u(G_PSI_EM)%d, w(H_chr)%d, nx, ny, nz, ng, par)

    return
  end subroutine outflow_boundary

  !--------------------------------------------------------------------
  !
  !
  !
  !--------------------------------------------------------------------
  subroutine fixed_boundary(udot, w, imask, par)
    implicit   none
    type(gridfunction), dimension(NU), intent(inout) :: udot
    type(gridfunction), dimension(NW) :: w
    CCTK_REAL, dimension(:)           :: par
    CCTK_INT, dimension(:,:,:)        :: imask

    ! local vars
!   CCTK_INT                          :: nx, ny, nz, ng
    
!   nx = nint(par(P_NX))  ! local size of arrays
!   ny = nint(par(P_NY))
!   nz = nint(par(P_NZ))
!   ng = nint(par(P_BOUND_WIDTH))

!   write(0,*)'fixed_boundary: Nothing here.'

    return
  end subroutine fixed_boundary

end module HYPER_BOUNDARY
