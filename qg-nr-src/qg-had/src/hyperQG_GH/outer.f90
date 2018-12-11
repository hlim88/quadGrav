#include "cctk.h"


module boundary_mine
   use params

   use GF
   implicit none

   contains 

  !--------------------------------------------------------------------
  !
  ! Not used...because of namespace issues, look in
  !           had/src/hyper/charboundary.f90
  !
  !--------------------------------------------------------------------
  subroutine outflow_boundary(u, w, imask, par)
    implicit   none
    type(gridfunction), dimension(NU_G), intent(inout) :: u
    type(gridfunction), dimension(NW) :: w
    CCTK_REAL, dimension(:)           :: par
    CCTK_INT, dimension(:,:,:)        :: imask

    ! local vars
    CCTK_INT                          :: nx, ny, nz, ng
    
    nx = nint(par(P_NX))  ! local size of arrays
    ny = nint(par(P_NY))
    nz = nint(par(P_NZ))
    ng = nint(par(P_GHOSTWIDTH))

    !call apply_outflow(u(H_D)%d,   nx, ny, nz, ng)
    !call apply_outflow(u(H_SX)%d,  nx, ny, nz, ng)
    !call apply_outflow(u(H_SY)%d,  nx, ny, nz, ng)
    !call apply_outflow(u(H_SZ)%d,  nx, ny, nz, ng)
    !call apply_outflow(u(H_TAU)%d, nx, ny, nz, ng)
    !call apply_outflow(u(H_BX)%d,  nx, ny, nz, ng)
    !call apply_outflow(u(H_BY)%d,  nx, ny, nz, ng)
    !call apply_outflow(u(H_BZ)%d,  nx, ny, nz, ng)
    !call apply_outflow(u(H_PSI)%d, nx, ny, nz, ng)

    return
  end subroutine outflow_boundary

  !--------------------------------------------------------------------
  !
  ! Not used...because of namespace issues, look in
  !           had/src/hyper/charboundary.f90
  !
  !
  !--------------------------------------------------------------------
  subroutine apply_outflow(f, nx, ny, nz, ng)
    implicit   none
    CCTK_INT                       :: nx, ny, nz, ng
    CCTK_REAL, dimension(nx,ny,nz) :: f

    ! local vars
    CCTK_INT                       :: i, j, k

    ! xmin
    do k = 1, nz
      do j = 1, ny
        do i = 1, ng
          f(i,j,k) = f(ng+1,j,k)
        end do
      end do
    end do

    ! xmax
    do k = 1, nz
      do j = 1, ny
        do i = nx-ng+1, nx
          f(i,j,k) = f(nx-ng,j,k)
        end do
      end do
    end do

    ! ymin
    do k = 1, nz
      do i = 1, nx
        do j = 1, ng
          f(i,j,k) = f(i,ng+1,k)
        end do
      end do
    end do

    ! ymax
    do k = 1, nz
      do i = 1, nx
        do j = ny-ng+1, ny
          f(i,j,k) = f(i,ny-ng,k)
        end do
      end do
    end do

    ! zmin
    do j = 1, ny
      do i = 1, nx
        do k = 1, ng
          f(i,j,k) = f(i,j,ng+1)
        end do
      end do
    end do

    ! zmax
    do j = 1, ny
      do i = 1, nx
        do k = nz-ng+1, nz
          f(i,j,k) = f(i,j,nz-ng)
        end do
      end do
    end do

    return
  end subroutine apply_outflow
  !--------------------------------------------------------------------
  !
  !
  !
  !--------------------------------------------------------------------
  subroutine my_boundary(u, udot, w, imask, par)
    implicit   none
    type(gridfunction), dimension(NU_G), intent(inout) :: u, udot
    type(gridfunction), dimension(NW) :: w
    CCTK_REAL, dimension(:)           :: par
    CCTK_INT, dimension(:,:,:)        :: imask
    return
  end subroutine my_boundary

!-----------------------------------------------------------------------
  ! Not used...because of namespace issues, look in
  !           had/src/hyper/charboundary.f90

 subroutine PWN2D_boundary(u, udot, v, w, imask, par)
    implicit   none
    type(gridfunction), dimension(NU_G), intent(inout) :: u, udot
    type(gridfunction), dimension(NW) :: w
    type(gridfunction), dimension(NV_G) :: v
    CCTK_REAL, dimension(:)           :: par
    CCTK_INT, dimension(:,:,:)        :: imask

    ! local vars

    CCTK_REAL  :: rho
    CCTK_REAL  :: vx
    CCTK_REAL  :: vy
    CCTK_REAL  :: vz
    CCTK_REAL  :: p
    CCTK_REAL  :: Bx
    CCTK_REAL  :: By
    CCTK_REAL  :: Bz


    CCTK_REAL, dimension(:,:,:), pointer :: xx, yy, zz
   
    CCTK_REAL, DIMENSION(3) :: n  ! outward unit normal
    CCTK_INT :: imin,imax,jmin,jmax,kmin,kmax, & 
                   nx, ny, nz, ng, k, i, j, l, NB    
    CCTK_INT :: idx, idy, idz, mask_stencil
    CCTK_INT :: IBBOX1, IBBOX2, IBBOX3,&
                IBBOX4, IBBOX5, IBBOX6  

    CCTK_REAL            :: R0     
    CCTK_REAL            :: xc     
    CCTK_REAL            :: yc     = 0.d0
    CCTK_REAL            :: zc     = 0.d0

    CCTK_REAL  gamma, h, vsq, Wsq, Bv, Bsq, hWsq, vabs, theta 

    logical bitantz
    logical double_equal

    logical, parameter :: ltrace = .false.
    logical, parameter :: debugg = .false.


 end subroutine PWN2D_boundary
!
!
!
  ! Not used...because of namespace issues, look in
  !           had/src/hyper/charboundary.f90

 
subroutine PWN3D_boundary(u, udot, v, w, imask, par)
    implicit   none
    type(gridfunction), dimension(NU_G), intent(inout) :: u, udot
    type(gridfunction), dimension(NW) :: w
    type(gridfunction), dimension(NV_G) :: v
    CCTK_REAL, dimension(:)           :: par
    CCTK_INT, dimension(:,:,:)        :: imask

    ! local vars

    CCTK_REAL  :: rho
    CCTK_REAL  :: vx
    CCTK_REAL  :: vy
    CCTK_REAL  :: vz
    CCTK_REAL  :: p
    CCTK_REAL  :: Bx
    CCTK_REAL  :: By
    CCTK_REAL  :: Bz


    CCTK_REAL, dimension(:,:,:), pointer :: xx, yy, zz
   
    CCTK_REAL, DIMENSION(3) :: n  ! outward unit normal
    CCTK_INT :: imin,imax,jmin,jmax,kmin,kmax, & 
                   nx, ny, nz, ng, k, i, j, l, NB    
    CCTK_INT :: idx, idy, idz, mask_stencil
    CCTK_INT :: IBBOX1, IBBOX2, IBBOX3,&
                IBBOX4, IBBOX5, IBBOX6  

    CCTK_REAL            :: R0     
    CCTK_REAL            :: xc     
    CCTK_REAL            :: yc     = 0.d0
    CCTK_REAL            :: zc     = 0.d0

    CCTK_REAL  :: gamma, h, vsq, Wsq, Bv, Bsq 
    CCTK_REAL  :: hWsq, vabs, theta, phi, r_xy, r_xyz 

    logical bitantz
    logical double_equal

    logical, parameter :: ltrace = .false.
    logical, parameter :: debugg = .false.



  end subroutine PWN3D_boundary
  
  
end module boundary_mine
