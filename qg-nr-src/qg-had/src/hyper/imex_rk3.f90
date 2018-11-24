!-------------------------------------------------------------------------
!
!  $Id$
!
!-------------------------------------------------------------------------

#include "cctk.h"


#ifndef IMEX_RK
!----------------------------------------------------------------------
!
! This is an empty routine that is used if IMEX_RK is not defined.
! IMEX_RK is an environment variable defined in the hsetup2 script
! when some of the evolved variables are stiff.
! 
!----------------------------------------------------------------------
subroutine IMEX_RK3(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1,  &
                    sr1, sr2, sr3, sr4, par)
  use params
  use UTILEQS
  use GF
  use rhs_mine
  use HYPER_DISSIPATION
  use HYPER_BOUNDARY
  use MOD_MASK
  use hypercoords
  use HYPER_AUXVARS
#ifdef HYPER_FLUID
  use m_defeos
#endif

  implicit none

  CCTK_INT     istat, rc, nt, nx, ny, nz, err, m, use_mask, dissipation
  CCTK_INT     update_scheme, amrbound_prepost
  CCTK_INT     boundary_conditions
  type(gridfunction), dimension(NU):: u0,u2, &
                                      urk1,dxu,dyu,dzu,  &
                                      sr1, sr2, sr3, sr4

  type(gridfunction), dimension(NV)  :: v, dxv, dyv, dzv
  type(gridfunction), dimension(NW)  :: w
  CCTK_REAL, dimension(NPAR)         :: par
  return
end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------


#else

subroutine IMEX_RK3(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1,  &
                    sr1, sr2, sr3, sr4, par)
  use params
  use UTILEQS
  use GF
  use rhs_mine
  use HYPER_DISSIPATION
  use HYPER_BOUNDARY
  use MOD_MASK
  use hypercoords
  use HYPER_AUXVARS
#ifdef HYPER_FLUID
  use m_defeos
#endif


  implicit none

  CCTK_INT     istat, rc, nt, nx, ny, nz, err, m, use_mask, dissipation
  CCTK_INT     update_scheme, amrbound_prepost
  CCTK_INT     boundary_conditions
  type(gridfunction), dimension(NU):: u0,u2, &
                                      urk1,dxu,dyu,dzu,  &
                                      sr1, sr2, sr3, sr4

  type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
  type(gridfunction), dimension(NW):: w
  CCTK_REAL, dimension(NPAR)           :: par

  CCTK_INT :: iter
  CCTK_REAL :: dt, dx, sigma_diss
  CCTK_INT, allocatable, dimension(:,:,:) :: imask, imask2
  CCTK_REAL, allocatable, dimension(:,:,:):: dxu_diss, dyu_diss, dzu_diss
  CCTK_INT CCTK_Equals
  CCTK_INT:: handle, err_red, i, j, k, order
  CCTK_REAL:: norm_rk, max_g11, min_g11
  CCTK_REAL, pointer, dimension(:,:,:) :: chr
  CCTK_REAL, dimension(:,:,:),pointer  :: xx, yy, zz
  CCTK_INT:: bc_type
  logical, parameter                   :: no_bc   = .false.
  ! variables for tracing
  logical, parameter                   :: debug1  = .false.
  logical, parameter                   :: ltrace  = .false.
  logical, parameter                   :: ltrace2 = .false.
  logical, parameter                   :: ltrace3 = .false.
  CCTK_INT  :: tubi, tuei, tvbi, tvei, midx, midy, midz, imex_type, eos_type
  CCTK_REAL :: myl2norm3d
  CCTK_REAL :: one_third,two_third,one_sixth,one_eigth,three_fourth,one_fourth
  CCTK_REAL :: alpha, beta, eta
  CCTK_REAL :: at21,at31,at32,at41,at42,at43,at51,at52,at53,at54
  CCTK_REAL :: a11,a21,a22,a31,a32,a33,a41,a42,a43,a44,a51,a52,a53,a54,a55
  CCTK_REAL :: omt1,omt2,omt3,omt4,omt5,om1,om2,om3,om4,om5
  CCTK_REAL :: aa
  CCTK_REAL :: time
  CCTK_INT, save :: first_call = 1

  one_third    = 1.0d0/3.0d0
  two_third    = 2.0d0/3.0d0
  one_fourth   = 1.0d0/4.0d0
  three_fourth = 3.0d0/4.0d0
  one_sixth    = 1.0d0/6.0d0
  one_eigth    = 1.0d0/8.0d0

  time = par(P_TIME)
  imex_type = nint(par(P_IMEX_TYPE))

  if (imex_type .EQ. 1) then
    ! parameters for the IMEX-SSP3(4,3,3)
    at21 = 0.0d0
    at31 = 0.0d0; at32 = 1.0d0
    at41 = 0.0d0; at42 = 1.0d0/4.0d0; at43 = 1.0d0/4.0d0
    at51 = 0.0d0; at52 = 1.0d0/6.0d0; at53 = 1.0d0/6.0d0; at54 = 4.0d0/6.0d0
    omt1 = 0.0d0; omt2 = 1.0d0/6.0d0; omt3 = 1.0d0/6.0d0; omt4 = 4.0d0/6.0d0; omt5 = 0.0d0

    alpha = 0.24169426078821
    beta  = 0.06042356519705
    eta   = 0.12915286960590

    a11 = alpha
    a21 =-alpha; a22 = alpha
    a31 = 0.0d0; a32 = 1.0d0 - alpha; a33 = alpha
    a41 = beta; a42 = eta; a43 = 0.5d0 - beta - eta - alpha; a44 = alpha
    a51 = 0.0d0; a52 = 1.0d0/6.0d0; a53 = 1.0d0/6.0d0; a54 = 4.0d0/6.0d0; a55 = 0.0d0
    om1 = 0.0d0; om2 = 1.0d0/6.0d0; om3 = 1.0d0/6.0d0; om4 = 4.0d0/6.0d0; om5 = 0.0d0

  else if (imex_type .EQ. 2) then
    at21 = 0.0d0
    at31 = 0.0d0; at32 = 1.0d0
    at41 = 0.0d0; at42 = 1.0d0/4.0d0; at43 = 1.0d0/4.0d0
    at51 = 0.0d0; at52 = 1.0d0/6.0d0; at53 = 1.0d0/6.0d0; at54 = 4.0d0/6.0d0
    omt1 = 0.0d0; omt2 = 1.0d0/6.0d0; omt3 = 1.0d0/6.0d0; omt4 = 4.0d0/6.0d0; omt5 = 0.0d0

    alpha = 4.0d0/6.0d0

    a11 = alpha
    a21 =-alpha; a22 = alpha
    a31 = 0.0d0; a32 = 1.0d0 - alpha; a33 = alpha
    a41 = (+ 2.0d0*alpha*alpha - 1.0d0 + 2.0d0*alpha)/(8.0d0*alpha)
    a42 = (- 4.0d0*alpha*alpha + 1.0d0)/(8.0d0*alpha)
    a43 = (- 3.0d0*alpha + 1.0d0)/4.0d0
    a44 = alpha
    a51 = 0.0d0; a52 = 1.0d0/6.0d0; a53 = 0.0d0; a54 = 4.0d0/6.0d0; a55 = 1.0d0/6.0d0
    om1 = 0.0d0; om2 = 1.0d0/6.0d0; om3 = 0.0d0; om4 = 4.0d0/6.0d0; om5 = 1.0d0/6.0d0

  else if (imex_type .EQ. 3) then
    ! parameters for the IMEX-SSP3(4,3,3), like the 1 
    at21 = 0.0d0
    at31 = 0.0d0; at32 = 1.0d0
    at41 = 0.0d0; at42 = 1.0d0/4.0d0; at43 = 1.0d0/4.0d0
    at51 = 0.0d0; at52 = 1.0d0/6.0d0; at53 = 1.0d0/6.0d0; at54 = 4.0d0/6.0d0
    omt1 = 0.0d0; omt2 = 1.0d0/6.0d0; omt3 = 1.0d0/6.0d0; omt4 = 4.0d0/6.0d0; omt5 = 0.0d0

    alpha = 0.5d0
    beta  = 1.0d0/8.0d0
    eta   = 0.0d0

    a11 = alpha
    a21 =-alpha; a22 = alpha
    a31 = 0.0d0; a32 = 1.0d0 - alpha; a33 = alpha
    a41 = beta; a42 = eta; a43 = 0.5d0 - beta - eta - alpha; a44 = alpha
    a51 = 0.0d0; a52 = 1.0d0/6.0d0; a53 = 1.0d0/6.0d0; a54 = 4.0d0/6.0d0; a55 = 0.0d0
    om1 = 0.0d0; om2 = 1.0d0/6.0d0; om3 = 1.0d0/6.0d0; om4 = 4.0d0/6.0d0; om5 = 0.0d0

  else
    print*, "this imex_type does not exist", imex_type
  end if

  ! use these variables to set output for traciing
  ! tubi - Traciing U Begin Index
  ! tuei - Traciing U End Index
  tubi = 1
  tuei = 9
  tuei = NU
  tvbi = 1
  tvei = NV

  iter = nint(par(P_RK_ITER))
  chr => w(H_CHR)%d
  xx  => w(H_XPHYS)%d
  yy  => w(H_YPHYS)%d
  zz  => w(H_ZPHYS)%d


  if (ltrace) then
     write(0,*)'>>>  IMEX RK3 begins'
  end if

  nx                    = nint(par(P_NX))
  ny                    = nint(par(P_NY))
  nz                    = nint(par(P_NZ))
  update_scheme         = nint(par(P_update_scheme))
  amrbound_prepost      = nint(par(P_amrbound_prepost))
  use_mask              = nint(par(P_USE_MASK))
  order                 = nint(par(P_DERIV_ORDER))
  dissipation           = nint(par(P_DISSIPATION))
  sigma_diss            = par(P_SIGMA_DISS)
  boundary_conditions   = nint(par(P_BOUNDARY_CONDITIONS))
  bc_type               = nint(par(P_BC_TYPE))

  allocate(dzu_diss(nx,ny,nz), STAT=istat)
  if (istat .ne. 0) then
    write(0,*)'>>> IMEX RK3: Could not allocate memory for dzu_diss'
    write(0,*)'     size : ',nx,ny,nz
  end if
  allocate(dyu_diss(nx,ny,nz), STAT=istat)
  if (istat .ne. 0) then
    write(0,*)'>>> IMEX RK3: Could not allocate memory for dyu_diss'
    write(0,*)'     size : ',nx,ny,nz
  end if
  allocate(dxu_diss(nx,ny,nz), STAT=istat)
  if (istat .ne. 0) then
    write(0,*)'>>> IMEX RK3: Could not allocate memory for dxu_diss'
    write(0,*)'     size : ',nx,ny,nz
  end if

  if (P_ENABLE_ALT_COORDS==1) then
    if (ltrace) then
      write(0,*)'>>> IMEX RK3: Calling def_alt_coords'
    end if
    call def_alt_coords(w, par)
  end if

  !
  ! This is not (always) the right value,
  ! it is the global dx:
  !
  dx = par(P_DX)
  dt = par(P_DT)
  midx = nint((nx - 1.)/2.+1)+5
  midy = nint((ny - 1.)/2.+1)+5
  midz = nint((nz - 1.)/2.+1)
!  midx = 1
!  midy = 1
!  midz = 1


  if (use_mask .eq. 1) then
    allocate(imask(nx,ny,nz),STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'IMEX RK3: Could not allocate memory for imask array'
      write(0,*)'     size : ',nx,ny,nz
    end if
    allocate(imask2(nx,ny,nz),STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'IMEX RK3: Could not allocate memory for imask2 array'
      write(0,*)'     size : ',nx,ny,nz
    end if

    call set_mask(w, imask,imask2, par)

  end if

  !###################################
  !######### 1st RK step #############
  !###################################

  if (iter .eq. 1) then

#ifdef HYPER_FLUID
  if (first_call .eq. 1) then
    first_call = 0
    eos_type = nint(par(P_EOS_TYPE))
    if (eos_type .lt. 100) then
      ! Ideal gas EOS or UWM EOS
      call qsetgeneos(par)
    else if (eos_type .eq. 100) then
      ! Ott EOS
      if (have_eos_table .eq. 0) then
        call initialize_nuc_eos(par, rc)
      end if
    else
      ! Error
      call my_exit('RK3TVD:  Unknown eos_type')
    end if
  end if
#endif

  ! copy from u0 to u2
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
     do m=1, NU
        u2(m)%d(i,j,k) = u0(m)%d(i,j,k)
     end do
  end do
  end do
  end do

  ! computing the auxiliary vars (ie, conductivities)
  call hyperauxvars(v, w, u2, par)      

  if (ltrace) then
     write(0,*)'>>> IMEX RK3: begin first step'
     write(0,*)'>>> IMEX RK3: dt,  dx = ',dt, dx
     pause
  end if

  if (debug1) then
     write(0,*)'>>> a point , before anything'
     write(0,*)'>>> mid=',midx,midy,midz
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey0=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q2=',u2(H_Q)%d(midx,midy,midz)
  end if

  if (ltrace2) then
     write(0,*)'>>>IMEX RK3: initial data'
     do m=tubi,tuei
        write(0,*) 'imex rk3: ||u||  ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
     end do
     do m=tvbi,tvei
        write(0,*) 'imex rk3: ||v||  ',m,myl2norm3d(v(m)%d, nx, ny, nz)
     end do
  end if

  if (ltrace) write(0,*)'>>> IMEX RK3: solving the first stiff step'

  ! the first substep is just the implicit part
#ifdef HYPER_FLUID
  ! initialize the fields for the first sub-step of the first RK step
  ! solve at the same time the implicit stiff part
  ! and the transformation from con to prim
  if (ltrace) write(0,*)'>>> IMEX RK3: Calling SolveStiffPrimVars'
  aa = a11 * dt
  call SolveStiffPrimVars(v, w, u2, sr1, aa, par)
#endif

  ! periodic boundary conditions in a single processor 
  ! this is not present in the 1D code, let us comment it by now
!  if (bc_type .eq. 3) then
!    if (ltrace) then
!      write(0,*) '>>> RK3: Calling periodic BC'
!    end if
!    do m=1, NU
!       call copy_bc2(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
!    end do
!  end if

  if (debug1) then
     write(0,*)'>>> after the first stiff step'
     write(0,*) "point x=",xx(midx,midy,midz),yy(midx,midy,midz),zz(midx,midy,midz)
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey0=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> sr1Ex=',sr1(H_EX)%d(midx,midy,midz),"sr1Ey=",sr1(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q=',u2(H_Q)%d(midx,midy,midz)
     pause
  end if

  if (ltrace2) then
     write(0,*)'>>>IMEX RK3: after the first stiff step'
     do m=tubi,tuei
        write(0,*) 'imex rk3: ||u||  ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
     end do
     do m=tvbi,tvei
        write(0,*) 'imex rk3: ||v||  ',m,myl2norm3d(v(m)%d, nx, ny, nz)
     end do
  end if

  ! the first substep of the first update only have implicit part
  ! for Steve: notice that sr1 does not depend on the neighbours, it is only point-wise
  ! computed. Should it be interpolated on the boundaries too?
  if (ltrace) write(0,*)'>>> IMEX RK3: update u1()'
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
!    if ( nint(chr(i,j,k)).eq.0 .and. amrbound_prepost.eq.1 ) then
      !
      do m=1, NU
        if (u2(m)%stiff == 1) then
          u2(m)%d(i,j,k) = u0(m)%d(i,j,k) + a21 * dt * sr1(m)%d(i,j,k)
!        else 
!          u2(m)%d(i,j,k) = u0(m)%d(i,j,k)
        end if
      end do
      !
!    end if
  end do
  end do
  end do

  if (ltrace2) then
     write(0,*)'>>>IMEX RK3: after the first update'
     do m=tubi,tuei
        write(0,*) 'imex rk3: ||u||  ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
     end do
     do m=tvbi,tvei
        write(0,*) 'imex rk3: ||v||  ',m,myl2norm3d(v(m)%d, nx, ny, nz)
     end do
  end if

  ! the second substep is the implicit part and explicit one
  if (ltrace) write(0,*)'>>> IMEX RK3: solving the second stiff step'
#ifdef HYPER_FLUID
  ! initialize the fields for the first sub-step of the first RK step
  aa = a22 * dt
  ! save the primitive fields and the stiff part at sr2
  call SolveStiffPrimVars(v, w, u2, sr2, aa, par)
#endif

  ! periodic boundary conditions in a single processor 
  if (bc_type .eq. 3) then
    do m=1, NU
       call copy_bc2(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
    end do
  end if

  if (debug1) then
     write(0,*)'>>> after the second stiff step'
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey0=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> sr2Ex=',sr2(H_EX)%d(midx,midy,midz),"sr2Ey=",sr2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q=',u2(H_Q)%d(midx,midy,midz)
     pause
  end if

  if (ltrace2) then
     write(0,*)'>>>IMEX RK3: after the second stiff step'
     do m=tubi,tuei
        write(0,*) 'imex rk3: ||u||  ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
     end do
     do m=tvbi,tvei
        write(0,*) 'imex rk3: ||v||  ',m,myl2norm3d(v(m)%d, nx, ny, nz)
     end do
  end if

  ! computing the auxiliary vars (ie, conductivities)
  call hyperauxvars(v, w, u2, par)      

  ! now the explicit part, notice that only the vs (primitive) should be
  ! used here for both the fluxes and the sources
  call calcrhs2( urk1, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, imask2, par)

  if (ltrace2) then
    write(0,*) ' IMEX RK3: After calling the first rhs '
    do m=tubi,tuei    
       write(0,*) ' u0   - A ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) ' urk1 - A ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if

  ! now it makes sense to compute dissipation
  if (dissipation > 0) then
    if (ltrace) write(0,*)'>>> IMEX RK3: hyperdissipation()', sigma_diss
    do m = 1, NU

      if (u2(m)%dissipation .eq. 1) then
        call hyperdissipation(dxu_diss, dyu_diss, dzu_diss, u2(m)%d, w, imask, par)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          urk1(m)%d(i,j,k) = urk1(m)%d(i,j,k)+ &
     &      W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
        end do
        end do
        end do
      end if

    end do
  end if
  par(P_RUNGE_KUTTA_BOUND) = 1.0
  call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)

  if (ltrace) write(0,*)'>>> IMEX RK3: update u2()'
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
!    if ( nint(chr(i,j,k)).eq.0 .and. amrbound_prepost.eq.1 ) then
      !
      ! If setting amr boundaries from parents after the evolution step,
      ! then we can go ahead and write over the boundary points:
      ! Otherwise, as long as the point is not an amr or dd boundary point,
      ! then go ahead and write over:
      !
      do m=1, NU
        if (u2(m)%stiff == 1) then
          u2(m)%d(i,j,k) = u0(m)%d(i,j,k) +  dt * urk1(m)%d(i,j,k) &
                         + dt * ( a31 * sr1(m)%d(i,j,k) + a32 * sr2(m)%d(i,j,k) ) 
        else
          u2(m)%d(i,j,k) = u0(m)%d(i,j,k) +  dt * urk1(m)%d(i,j,k)
        end if
      end do
      !
!    end if
  end do
  end do
  end do

  if (ltrace2) then
    write(0,*) ' After updating the first rk '
    do m=tubi,tuei
       write(0,*) ' u0   - B ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) ' u2   - B ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) ' urk1 - B ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if

  ! now that it has been update we can impose the boundary conditions
  par(P_RUNGE_KUTTA_BOUND) = 2.0
  call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)

  if (ltrace2) then
    write(0,*) ' After applying the first boundary condition '
    do m=tubi,tuei
       write(0,*) ' u0   - C ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) ' u2   - C ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) ' urk1 - C ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    do m=tvbi,tvei
       write(0,*) 'v - C  ',m,myl2norm3d(v(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if

  if (use_mask .eq. 1) then
    if (ltrace) write(0,*)'>>> IMEX RK3: use_mask'
    call mask_gfuncs(u2, w, par)
  end if

  if (debug1) then
     write(0,*)'>>> before the third stiff step'
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey2=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q2=',u2(H_Q)%d(midx,midy,midz)
     write(0,*)'>>> Ex0=',u0(H_EX)%d(midx,midy,midz),"Ey0=",u0(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q0=',u0(H_Q)%d(midx,midy,midz)
     write(0,*)'>>> sr*=',dt*(a31*sr1(H_EX)%d(midx,midy,midz)+a32*sr2(H_EX)%d(midx,midy,midz))
     write(0,*)'>>> D0=',u0(H_D)%d(midx,midy,midz),"tau0=",u0(H_TAU)%d(midx,midy,midz)
     write(0,*)'>>> S0x=',u0(H_SX)%d(midx,midy,midz),"S0y=",u0(H_SY)%d(midx,midy,midz)
     write(0,*)'>>> D2=',u2(H_D)%d(midx,midy,midz),"tau2=",u2(H_TAU)%d(midx,midy,midz)
     write(0,*)'>>> S2x=',u2(H_SX)%d(midx,midy,midz),"S2y=",u2(H_SY)%d(midx,midy,midz)
     write(0,*)'>>> rkD=',urk1(H_D)%d(midx,midy,midz),"rktau=",urk1(H_TAU)%d(midx,midy,midz)
     write(0,*)'>>> rkSx=',urk1(H_SX)%d(midx,midy,midz),"rkSy=",urk1(H_SY)%d(midx,midy,midz)
  end if

#ifdef HYPER_BSSN
  call EnforceBSSN(u2, v, w, par)
#endif

  if (ltrace) write(0,*)'>>> IMEX RK3: solving the third stiff step'
  ! this is the last sub-step of the first update
#ifdef HYPER_FLUID
  aa = a33 * dt
  ! save the primitive fields and the stiff part at sr3
  call SolveStiffPrimVars(v, w, u2, sr3, aa, par)
#endif

 ! periodic boundary conditions in a single processor 
  if (bc_type .eq. 3) then
      do m=1, NU
         call copy_bc2(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      end do
  end if

  if (debug1) then
     write(0,*)'>>> after the third stiff step'
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey0=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> sr3Ex=',sr3(H_EX)%d(midx,midy,midz),"sr3Ey=",sr3(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q=',u2(H_Q)%d(midx,midy,midz)
     pause
  end if

  if (ltrace) then
     write(0,*)'>>> IMEX RK3: Done step 1'
  end if

  !######### 2nd RK step #############
  else if (iter .eq. 2) then

  if (ltrace) then
     write(0,*)'>>> IMEX RK3: begin second step'
  end if

  if (ltrace3) then
     write(0,*)'>>> before calling rhs in the 2nd RK step'
     write(0,*)'>>> mid=',midx,midy,midz
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey2=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q2=',u2(H_Q)%d(midx,midy,midz)
     write(0,*)'>>> Ex0=',u0(H_EX)%d(midx,midy,midz),"Ey0=",u0(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q0=',u0(H_Q)%d(midx,midy,midz)
     write(0,*)'>>> sr=',sr1(H_EX)%d(midx,midy,midz),sr2(H_EX)%d(midx,midy,midz),sr3(H_EX)%d(midx,midy,midz)
  end if

  ! computing the auxiliary vars (ie, conductivities)
  call hyperauxvars(v, w, u2, par)      

  ! computing the rhs
  call calcrhs2(urk1, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, imask2, par)

  if (ltrace2) then
    write(0,*) ' After calling the second rhs '
    do m=tubi,tuei    
       write(0,*) ' urk2 - A ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if
  
  if (dissipation > 0) then
    if (ltrace) write(0,*)'>>>IMEX  RK3: hyperdissipation()', sigma_diss
    do m=1,NU
      if (u2(m)%dissipation .eq. 1) then
        call hyperdissipation(dxu_diss, dyu_diss, dzu_diss, u2(m)%d, w, imask, par)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          urk1(m)%d(i,j,k) = urk1(m)%d(i,j,k)+ &
     &       W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
        end do
        end do
        end do
      end if
    end do
  end if
  par(P_RUNGE_KUTTA_BOUND) = 1.0
  call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)

  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
!    if ( nint(chr(i,j,k)).eq.0 .and. amrbound_prepost.eq.1 ) then
      !
      do m=1, NU
        if (u2(m)%stiff == 1) then
           u2(m)%d(i,j,k) = three_fourth * u0(m)%d(i,j,k) + one_fourth * u2(m)%d(i,j,k) &
                          + one_fourth * dt * urk1(m)%d(i,j,k) &
                          + dt * (  (a41 - one_fourth*a31) * sr1(m)%d(i,j,k) &
                                 +  (a42 - one_fourth*a32) * sr2(m)%d(i,j,k) &
                                 +  (a43 - one_fourth*a33) * sr3(m)%d(i,j,k)  )
        else
           u2(m)%d(i,j,k) = three_fourth * u0(m)%d(i,j,k) + one_fourth * u2(m)%d(i,j,k) &
                          + one_fourth * dt * urk1(m)%d(i,j,k) 
        end if
      end do
      !
!    end if
  end do
  end do
  end do

  if (debug1) then
     write(0,*)'>>> after the fourth stiff step'
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey2=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> sr4Ex=',sr4(H_EX)%d(midx,midy,midz),"sr4Ey=",sr4(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q=',u2(H_Q)%d(midx,midy,midz)
     pause
  end if

  if (ltrace2) then
    write(0,*) ' After updating the second rk '
    do m=tubi,tuei
       write(0,*) ' u2   - B ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) ' urk2 - B ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if

  ! now that it has been update we can impose the boundary conditions
  par(P_RUNGE_KUTTA_BOUND) = 2.0
  call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)

  if (use_mask .eq. 1) then
    call mask_gfuncs(u2, w, par)
  end if

  if (ltrace3) then
     write(0,*)'>>> before the fourth stiff step'
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey2=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q2=',u2(H_Q)%d(midx,midy,midz)
  end if

#ifdef HYPER_BSSN
  call EnforceBSSN(u2, v, w, par)
#endif

  if (ltrace) write(0,*)'>>> IMEX RK3: solving the fourth stiff step'
  ! this is the last sub-step of the second update
#ifdef HYPER_FLUID
  aa = a44 * dt
  ! save the primitive fields and the stiff part at sr1
  call SolveStiffPrimVars(v, w, u2, sr4, aa, par)
#endif

 ! periodic boundary conditions in a single processor 
  if (bc_type .eq. 3) then
      do m=1, NU
         call copy_bc2(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      end do
  end if

  if (ltrace3) then
     write(0,*)'>>> after the fourth stiff step'
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey2=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> sr4Ex=',sr4(H_EX)%d(midx,midy,midz),"sr4Ey=",sr4(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q=',u2(H_Q)%d(midx,midy,midz)
  end if

  !######### 3rd RK step ##############
  else if (iter .eq. 3) then
  if (ltrace) then
     write(0,*)'>>> IMEX RK3: begin third step'
  end if

  ! computing the auxiliary vars (ie, conductivities)
  call hyperauxvars(v, w, u2, par)      

  ! compute the rhs
  call calcrhs2( urk1, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, imask2, par)

  if (debug1) then
     write(0,*) "step 5 just after calrhs"
     write(0,*) "point x=",xx(midx,midy,midz),yy(midx,midy,midz),zz(midx,midy,midz)
     write(0,*)' D=',u2(H_D)%d(midx,midy,midz),"tau=",u2(H_TAU)%d(midx,midy,midz)
     write(0,*)' S=',u2(H_SX)%d(midx,midy,midz),u2(H_SY)%d(midx,midy,midz),u2(H_SZ)%d(midx,midy,midz)
     write(0,*)' B=',u2(H_BX)%d(midx,midy,midz),u2(H_BY)%d(midx,midy,midz),u2(H_BZ)%d(midx,midy,midz)
     write(0,*)' E=',u2(H_EX)%d(midx,midy,midz),u2(H_EY)%d(midx,midy,midz),u2(H_EZ)%d(midx,midy,midz)
     write(0,*)'rkE=',urk1(H_EX)%d(midx,midy,midz),urk1(H_EY)%d(midx,midy,midz),urk1(H_EZ)%d(midx,midy,midz)
     write(0,*)' q=',u2(H_Q)%d(midx,midy,midz),"rkq",urk1(H_Q)%d(midx,midy,midz)
     pause
  end if


  if (ltrace2) then
    write(0,*) ' After calling the third rhs '
    do m=tubi,tuei    
       write(0,*) ' urk3 - A ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if
  
  if (dissipation > 0) then
    if (ltrace) write(0,*)'>>>IMEX RK3: hyperdissipation()', sigma_diss
    do m=1,NU
      if (u2(m)%dissipation .eq. 1) then
        call hyperdissipation(dxu_diss, dyu_diss, dzu_diss, u2(m)%d, w, imask, par)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          urk1(m)%d(i,j,k) = urk1(m)%d(i,j,k)+ &
     &       W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
        end do
        end do
        end do
      end if
    end do
  end if

  par(P_RUNGE_KUTTA_BOUND) = 1.0
  call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)

  if (ltrace2) then
    write(0,*) ' After applying the third boundary condition '
    do m=tubi,tuei
       write(0,*) ' u2   - C ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) ' urk3 - C ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if

  if (debug1) then
     write(0,*)'>>> before 5 before explicit updating'
     write(0,*) "point x=",xx(midx,midy,midz),yy(midx,midy,midz),zz(midx,midy,midz)
     write(0,*)' D=',u2(H_D)%d(midx,midy,midz),"tau=",u2(H_TAU)%d(midx,midy,midz)
     write(0,*)' S=',u2(H_SX)%d(midx,midy,midz),u2(H_SY)%d(midx,midy,midz),u2(H_SZ)%d(midx,midy,midz)
     write(0,*)' B=',u2(H_BX)%d(midx,midy,midz),u2(H_BY)%d(midx,midy,midz),u2(H_BZ)%d(midx,midy,midz)
     write(0,*)' E=',u2(H_EX)%d(midx,midy,midz),u2(H_EY)%d(midx,midy,midz),u2(H_EZ)%d(midx,midy,midz)
     write(0,*)'rkE=',urk1(H_EX)%d(midx,midy,midz),urk1(H_EY)%d(midx,midy,midz),urk1(H_EZ)%d(midx,midy,midz)
     write(0,*)' q=',u2(H_Q)%d(midx,midy,midz),"rkq",urk1(H_Q)%d(midx,midy,midz)
     write(0,*)'>>> urk1Ex=',urk1(H_EX)%d(midx,midy,midz)
     write(0,*)'>>> sr1Ex=',sr1(H_EX)%d(midx,midy,midz)
     write(0,*)'>>> sr2Ex=',sr2(H_EX)%d(midx,midy,midz)
     write(0,*)'>>> sr3Ex=',sr3(H_EX)%d(midx,midy,midz)
     write(0,*)'>>> sr4Ex=',sr4(H_EX)%d(midx,midy,midz)
     pause
  end if


  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
!    if ( nint(chr(i,j,k)).eq.0 .and. amrbound_prepost.eq.1 ) then
      !
      do m=1, NU
        if (u2(m)%stiff == 1) then
           u2(m)%d(i,j,k) = one_third * u0(m)%d(i,j,k) + two_third * u2(m)%d(i,j,k) &
                + two_third * dt * urk1(m)%d(i,j,k) &
                + dt * (  (a51 - two_third * a41) * sr1(m)%d(i,j,k) &
                       +  (a52 - two_third * a42) * sr2(m)%d(i,j,k) &
                       +  (a53 - two_third * a43) * sr3(m)%d(i,j,k) &
                       +  (a54 - two_third * a44) * sr4(m)%d(i,j,k)  )
        else
           u2(m)%d(i,j,k) = one_third * u0(m)%d(i,j,k) + two_third * u2(m)%d(i,j,k) &
                + two_third * dt * urk1(m)%d(i,j,k)
        end if
      end do
      !
!    end if
  end do
  end do
  end do

  if (debug1) then
     write(0,*) "step 5 after explicit updating"
     write(0,*) "point x=",xx(midx,midy,midz),yy(midx,midy,midz),zz(midx,midy,midz)
     write(0,*)' D=',u2(H_D)%d(midx,midy,midz),"tau=",u2(H_TAU)%d(midx,midy,midz)
     write(0,*)' S=',u2(H_SX)%d(midx,midy,midz),u2(H_SY)%d(midx,midy,midz),u2(H_SZ)%d(midx,midy,midz)
     write(0,*)' B=',u2(H_BX)%d(midx,midy,midz),u2(H_BY)%d(midx,midy,midz),u2(H_BZ)%d(midx,midy,midz)
     write(0,*)' E=',u2(H_EX)%d(midx,midy,midz),u2(H_EY)%d(midx,midy,midz),u2(H_EZ)%d(midx,midy,midz)
     pause
  end if

  if (ltrace2) then
    write(0,*) ' After updating the third rk '
    do m=tubi,tuei
       write(0,*) ' u2   - B ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) ' urk3 - B ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if

  if (ltrace3) then
     write(0,*)'>>> the last step of the RK does not contain solving stiff part'
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey2=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q2=',u2(H_Q)%d(midx,midy,midz)
  end if

#ifdef HYPER_BSSN
  call EnforceBSSN(u2, v, w, par)
#endif

  if (ltrace) write(0,*)'>>> IMEX RK3: inverting the last contoprim'
#ifdef HYPER_FLUID
  ! perform a simple conversion from conserved to primitive
  aa = a55 * dt
  call SolveStiffPrimVars(v, w, u2, sr4, aa, par)
#endif

 ! periodic boundary conditions in a single processor 
  if (bc_type .eq. 3) then
      do m=1, NU
         call copy_bc2(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      end do
  end if

  if (ltrace3) then
     write(0,*)'>>> after the last stiff step'
     write(0,*)'>>> Ex2=',u2(H_EX)%d(midx,midy,midz),"Ey2=",u2(H_EY)%d(midx,midy,midz)
     write(0,*)'>>> q2=',u2(H_Q)%d(midx,midy,midz)
  end if

  ! applying the boundary conditions  
  par(P_RUNGE_KUTTA_BOUND) = 2.0
  call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)

  if (use_mask .eq. 1) then
    call mask_gfuncs(u2, w, par)
  end if

  if (bc_type .eq. 3) then
      do m=1, NU
         call copy_bc2(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      end do
  end if

  if (ltrace) then
     !call output3DB( u2(H_Ez)%d,"Ezitr3end",xx(1,1,1),xx(nx,1,1),yy(1,1,1),yy(1,ny,1),zz(1,1,1),zz(1,1,nz),time,0,nx,ny,nz)
     write(0,*) '||Ezitr3end||  ',myl2norm3d(u2(H_Ez)%d, nx, ny, nz)
  end if

  end if !---iter



  if (allocated(dxu_diss) .and. allocated(dyu_diss) .and. allocated(dzu_diss))&
  then
    if (ltrace2) write(0,*)'>>> IMEX RK3: about to deallocate diss vars...'
    deallocate(dzu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'>>> IMEX RK3: Could not deallocate dzu_diss: istat=', istat
    end if
    deallocate(dxu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'>>> IMEX RK3: Could not deallocate dxu_diss: istat=', istat
    end if
    deallocate(dyu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'>>> IMEX RK3: Could not deallocate dyu_diss: istat=', istat
    end if
    if (ltrace2) write(0,*)'>>> RK3: ...diss vars deallocated'
  else
    write(0,*)'>>> IMEX RK3: diss vars not allocated as expected!'
  end if


  if (use_mask .eq. 1) then
    !w(H_MASK)%d = real(imask)
    if (allocated(imask)) then
      deallocate(imask,imask2)
    else
      write(0,*)'>>> IMEX RK3: imask not allocated as expected!'
    end if
  end if

  if (ltrace) then
    do i = 1, nx
    do j = 1, ny
    do k = 1, nz
      if (abs(w(H_MASK)%d(i,j,k)) .lt. 1.0e-6) then
        print *,'IMEX RK3TVD: mask problem, end rk3tvd'
        stop
      end if
    end do
    end do
    end do
  end if


  if (ltrace) then
     write(0,*)'>>>IMEX RK3: return'
     if (ltrace2) then
        do m=tubi,tuei
           write(0,*) ' u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
           write(0,*) ' u0   - ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
           write(0,*) ' urk3 - ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
           write(0,*) ' '
        end do
     end if
  end if

  return
end subroutine IMEX_RK3
#endif
