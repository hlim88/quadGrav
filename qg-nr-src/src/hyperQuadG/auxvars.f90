!----------------------------------------------------------------------
!
!  $Id: auxvars.f90,v 1.3 2010-07-30 00:06:38 ehirsch Exp $
!
!----------------------------------------------------------------------

#include "cctk.h"


#if defined GR_MHD
  module bssn_auxvars
#else
  module auxvars
#endif


  use params

  implicit none
  CONTAINS

  !--------------------------------------------------------------------
  !
  ! Define the stress energy tensor for 1 complex and
  ! 1 real massless scalar field
  !
  !--------------------------------------------------------------------
  
#if defined GR_MHD
  subroutine bssn_define_aux_vars(v_pt, w_pt, u_pt, par)
#else
  subroutine define_aux_vars(v_pt, w_pt, u_pt, par)
#endif

    CCTK_REAL, dimension(NU_G), intent(IN)    :: u_pt   
    CCTK_REAL, dimension(NV_G), intent(INOUT) :: v_pt
    CCTK_REAL, dimension(NW)                :: w_pt
    CCTK_REAL, dimension(:)                 :: par

    
    !! local vars
    logical, parameter               :: ltrace2 = .false.

    !myid     = proc_return_myid()

    return
  end subroutine

  !--------------------------------------------------------------------
  !
  !
  !
  !--------------------------------------------------------------------
  subroutine EnforceBSSN(u, v, w, par)
    use GF
    use m_enforce_bssn
    implicit none

    type(gridfunction), dimension(NV) :: v
    type(gridfunction), dimension(NU) :: u
    type(gridfunction), dimension(NW) :: w
    CCTK_REAL, dimension(:)           :: par

    ! local vars
    CCTK_INT, dimension(3) :: shp
    CCTK_INT               :: nx, ny, nz
!    CCTK_INT               :: ohmlaw_type
    CCTK_INT               :: force_free 
    CCTK_REAL, dimension(:,:,:), pointer :: chi
    CCTK_REAL, dimension(:,:,:), pointer :: gtxx, gtxy, gtxz, gtyy, gtyz, gtzz

    shp     = shape(u(1)%d)
    nx      = shp(1)
    ny      = shp(2)
    nz      = shp(3)

    if (nx .ne. nint(par(P_NX))) then
      write(0,*)'*** auxvars.f90:  Problem with inconsistent', &
     &              ' array sizes'
      write(0,*) 'nx        = ',nx
      write(0,*) 'par(P_NX) = ',nint(par(P_NX))
    end if
    if (ny .ne. nint(par(P_NY))) then
      write(0,*)'*** auxvars.f90:  Problem with inconsistent', &
     &              ' array sizes'
      write(0,*) 'ny        = ',ny
      write(0,*) 'par(P_NY) = ',nint(par(P_NY))
    end if
    if (nz .ne. nint(par(P_NZ))) then
      write(0,*)'*** auxvars.f90:  Problem with inconsistent', &
     &              ' array sizes'
      write(0,*) 'nz        = ',nz
      write(0,*) 'par(P_NZ) = ',nint(par(P_NZ))
    end if


    chi         => u(H_CHI)%d
    gtxx        => u(H_GT11)%d
    gtxy        => u(H_GT12)%d
    gtxz        => u(H_GT13)%d
    gtyy        => u(H_GT22)%d
    gtyz        => u(H_GT23)%d
    gtzz        => u(H_GT33)%d

!   if (apply_evil_bc) then
!     call bssn_evil_mpi(u, v, w, par)
!   end if

    call enforce_bssn_constraints(u, w, par)

    call cal_adm_metric(                                            &
     &                    v(G_G11)%d,  v(G_G12)%d,  v(G_G13)%d,     &
     &                    v(G_G22)%d,  v(G_G23)%d,  v(G_G33)%d,     &
     &                    gtxx, gtxy, gtxz, gtyy, gtyz, gtzz,       &
     &                    chi,  nx,  ny,  nz)

!     ohmlaw_type = nint(par(P_OHMLAW_TYPE))
     force_free = nint(par(P_FORCE_FREE))

!    if (ohmlaw_type == 1) then
    if (force_free == 1) then
      call enforce_force_free(u, v, w, par)
    end if

    
  end subroutine EnforceBSSN

  !---------------------------------------------------------------------
  !
  !
  !
  !---------------------------------------------------------------------
  subroutine enforce_force_free(u, v, w, par)
    use GF
    implicit none
    type(gridfunction), dimension(NV_G) :: v
    type(gridfunction), dimension(NU_G) :: u
    type(gridfunction), dimension(NW) :: w
    CCTK_REAL, dimension(:)           :: par

    ! local vars
    CCTK_INT, dimension(3)               :: shp
    CCTK_INT                             :: nx, ny, nz, i, j, k
    CCTK_REAL                            :: EB, Esq, Bsq
    CCTK_REAL, dimension(3)              :: Eu, Bu
    CCTK_REAL, dimension(3,3)            :: gd
    CCTK_REAL, dimension(:,:,:), pointer :: gd11, gd12, gd13, gd22, gd23, gd33
    CCTK_REAL, dimension(:,:,:), pointer :: Ex, Ey, Ez, Bx, By, Bz
    logical, parameter                   :: ltrace = .false.

    if (ltrace) then
      write(0,*)'enforce_force_free: begin '
    end if
    shp     = shape(u(1)%d)
    nx      = shp(1)
    ny      = shp(2)
    nz      = shp(3)

    gd11 => v(G_G11)%d
    gd12 => v(G_G12)%d
    gd13 => v(G_G13)%d
    gd22 => v(G_G22)%d
    gd23 => v(G_G23)%d
    gd33 => v(G_G33)%d

    Ex   => u(G_EX)%d
    Ey   => u(G_EY)%d
    Ez   => u(G_EZ)%d
    Bx   => u(G_BX)%d
    By   => u(G_BY)%d
    Bz   => u(G_BZ)%d
    
    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
         gd(1,1) = gd11(i,j,k)
         gd(1,2) = gd12(i,j,k)
         gd(1,3) = gd13(i,j,k)
         gd(2,1) = gd(1,2)
         gd(2,2) = gd22(i,j,k)
         gd(2,3) = gd23(i,j,k)
         gd(3,1) = gd(1,3)
         gd(3,2) = gd(2,3)
         gd(3,3) = gd33(i,j,k)

         Eu(1)   = Ex(i,j,k)
         Eu(2)   = Ey(i,j,k)
         Eu(3)   = Ez(i,j,k)

         Bu(1)   = Bx(i,j,k)
         Bu(2)   = By(i,j,k)
         Bu(3)   = Bz(i,j,k)

         EB = gd(1,1)*Eu(1)*Bu(1)+gd(1,2)*Eu(1)*Bu(2) &
     &      + gd(1,3)*Eu(1)*Bu(3)+gd(1,2)*Eu(2)*Bu(1) &
     &      + gd(2,2)*Eu(2)*Bu(2)+gd(2,3)*Eu(2)*Bu(3) &
     &      + gd(1,3)*Eu(3)*Bu(1)+gd(2,3)*Eu(3)*Bu(2) &
     &      + gd(3,3)*Eu(3)*Bu(3)
 
         Bsq = gd(1,1)*Bu(1)*Bu(1) + 2.0d0*gd(1,2)*Bu(1)*Bu(2)       &
     &          + 2.0d0*gd(1,3)*Bu(1)*Bu(3) + gd(2,2)*Bu(2)*Bu(2)    &
     &          + 2.0d0*gd(2,3)*Bu(2)*Bu(3) + gd(3,3)*Bu(3)*Bu(3) 

         if (Bsq > 1.0d-16) then
           Eu(1) = Eu(1) - EB/Bsq*Bu(1)
           Eu(2) = Eu(2) - EB/Bsq*Bu(2)
           Eu(3) = Eu(3) - EB/Bsq*Bu(3)
         end if

         Esq = gd(1,1)*Eu(1)*Eu(1) + 2.0d0*gd(1,2)*Eu(1)*Eu(2)      &
     &          + 2.0d0*gd(1,3)*Eu(1)*Eu(3) + gd(2,2)*Eu(2)*Eu(2)   &
     &          + 2.0d0*gd(2,3)*Eu(2)*Eu(3) + gd(3,3)*Eu(3)*Eu(3)

         if (Esq > Bsq) then
           if (Esq > 1.0d-16) then
             Eu(1) = Eu(1)*sqrt(Bsq/Esq)
             Eu(2) = Eu(2)*sqrt(Bsq/Esq)
             Eu(3) = Eu(3)*sqrt(Bsq/Esq)
           else 
             Eu(1) = 0.0d0
             Eu(2) = 0.0d0
             Eu(3) = 0.0d0
           end if
         end if

         Ex(i,j,k) = Eu(1)
         Ey(i,j,k) = Eu(2)
         Ez(i,j,k) = Eu(3)
    end do
    end do
    end do

    return
  end subroutine enforce_force_free


end module
