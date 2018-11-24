!-------------------------------------------------------------------------
!
!  $Id$
!
!-------------------------------------------------------------------------

#include "cctk.h"


subroutine TVDRK3_himem(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
     urk2, urk3, par)
  use params
  use UTILEQS
  use GF
  use rhs_mine
  use HYPER_DISSIPATION
  use HYPER_BOUNDARY
  use MOD_MASK
  use hypercoords
  use HYPER_AUXVARS



  implicit none

  CCTK_INT     istat, rc, nt, nx, ny, nz, err, m, use_mask, dissipation
  CCTK_INT     update_scheme, amrbound_prepost
  CCTK_INT     boundary_conditions
  type(gridfunction), dimension(NU):: u0,u2, &
                                      urk1,urk2,urk3,&
                                      dxu,dyu,dzu

  type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
  type(gridfunction), dimension(NW):: w
  CCTK_REAL, dimension(NPAR)           :: par



  CCTK_INT  :: iter
  CCTK_REAL :: dt, dx, sigma, sigma_diss
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
  logical, parameter                   :: ltrace  = .false.
  logical, parameter                   :: ltrace2 = .false.
  CCTK_INT tubi, tuei, tvbi, tvei
  CCTK_REAL :: myl2norm3d

  ! use these variables to set output for traciing
  ! tubi - Traciing U Begin Index
  ! tuei - Traciing U End Index
  tubi = 1
  tuei = 9
  tvbi = 1
  tvei = 0

  iter = nint(par(P_RK_ITER))
  chr => w(H_CHR)%d
  xx  => w(H_XPHYS)%d
  yy  => w(H_YPHYS)%d
  zz  => w(H_ZPHYS)%d
                                                                                                                        
  !hi   = xx(2,1,1)-xx(1,1,1)
  !
  ! Need minimums for this grid alone:
  !
  !minx = xx(1,1,1)
  !miny = yy(1,1,1)
  !minz = zz(1,1,1)
  sigma =0. ! move this later 

  if (ltrace) then
     write(0,*)'>>>  RK3_himem begins.  iter=',iter
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
    write(0,*)'>>> RK3_himem: Could not allocate memory for dzu_diss'
    write(0,*)'     size : ',nx,ny,nz
  end if
  allocate(dyu_diss(nx,ny,nz), STAT=istat)
  if (istat .ne. 0) then
    write(0,*)'>>> RK3_himem: Could not allocate memory for dyu_diss'
    write(0,*)'     size : ',nx,ny,nz
  end if
  allocate(dxu_diss(nx,ny,nz), STAT=istat)
  if (istat .ne. 0) then
    write(0,*)'>>> RK3_himem: Could not allocate memory for dxu_diss'
    write(0,*)'     size : ',nx,ny,nz
  end if

  if (P_ENABLE_ALT_COORDS==1) then
    if (ltrace) then
      write(0,*)'>>> RK3_himem: Calling def_alt_coords'
    end if
    call def_alt_coords(w, par)
  end if

  !
  ! This is not (always) the right value,
  ! it is the global dx:
  !
  dx = par(P_DX)
  dt = par(P_DT)

  if (use_mask .eq. 1) then
    allocate(imask(nx,ny,nz),STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'RK3_himemTVD: Could not allocate memory for imask array'
      write(0,*)'     size : ',nx,ny,nz
    end if
    allocate(imask2(nx,ny,nz),STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'RK3_himemTVD: Could not allocate memory for imask2 array'
      write(0,*)'     size : ',nx,ny,nz
    end if


    call set_mask(w, imask,imask2, par)

  end if

  !###################################
  !######### 1st RK step #############
  !###################################

  if (iter .eq. 1) then
  if (ltrace) then
     write(0,*)'>>> RK3_himem: begin first step'
     write(0,*)'>>> RK3_himem: dt,  dx = ',dt, dx
  end if

  if (ltrace2) then
     write(0,*)'>>>RK3_himem: initial data'
     do m=tubi,tuei
        write(0,*) '||u||  ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
     end do
     do m=tvbi,tvei
        write(0,*) '||v||  ',m,myl2norm3d(v(m)%d, nx, ny, nz)
     end do
  end if

  if ( par(P_ALLOWEDL) .gt. 0 ) then
    ! Ensure that primitive variables are never interpolated
    ! Problem interpolation will only occur when there are multiple levels
#ifdef HYPER_FLUID
    call SolvePrimVars(v, w, u0, par)
#endif
  end if

  call calcrhs2(   urk1,  u0,   v,    w,  dxu,  dyu,  dzu,  dxv, dyv,&
       & dzv, imask2, par)

  if (ltrace2) then
    do m=tubi,tuei
       write(0,*) ' u0   - A ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) ' urk1 - A ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if

  if (dissipation > 0) then
    if (ltrace) write(0,*)'>>> RK3_himem: hyperdissipation()', sigma_diss
    do m = 1, NU
      if (u0(m)%dissipation .eq. 1) then
        call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,   u0(m)%d&
             &, w,         imask,     par)
!        urk1(m)%d = urk1(m)%d + W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss +  dyu_diss +  dzu_diss)
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
  call hyperboundary(u0, dxu, dyu, dzu, urk1, v, w, imask, par)

  if (ltrace) write(0,*)'>>> RK3_himem: update u2()'
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
!   if (amrbound_prepost.eq.1          &
!       .or. ( nint(chr(i,j,k)) .ne. 1 .and.  &
!              nint(chr(i,j,k)) .ne. 8      )) then
    if (       (.not. no_bc .or. (no_bc .and. nint(chr(i,j,k)).eq.0))        &
         .and. (      amrbound_prepost.eq.1                                  &
                .or.  ( nint(chr(i,j,k)) .ne. 1 .and. nint(chr(i,j,k)).ne.8))&
        ) then
      !
      ! If setting amr boundaries from parents after the evolution step,
      ! then we can go ahead and write over the boundary points:
      ! Otherwise, as long as the point is not an amr or dd boundary point,
      ! then go ahead and write over:
      !
      do m=1, NU
        u2(m)%d(i,j,k) = u0(m)%d(i,j,k) + dt*urk1(m)%d(i,j,k)
      end do
    else if (nint(chr(i,j,k)) .eq. 1 .and. amrbound_prepost.eq.0) then
      ! Need to interpolate on the boundaries so that
      ! when complete step taken, u2(i,j,k) has the values
      ! given from the parents:
      do m=1, NU
         u2(m)%d(i,j,k) = 0.5d0*( u0(m)%d(i,j,k) + u2(m)%d(i,j,k) )
      end do
    end if
  end do
  end do
  end do

  if (ltrace2) then
    do m=tubi,tuei
       write(0,*) ' u0   - B ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) ' u2   - B ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) ' urk1 - B ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if


!
!  boundary_conditions .eq. 9  ---> PWN2D_boundary
!  boundary_conditions .eq. 10 ---> PWN3D_boundary
!
! if ( (boundary_conditions .eq. 8).OR. & 
!      (boundary_conditions .eq. 9).OR. &
!      (boundary_conditions .eq. 10) )then
      par(P_RUNGE_KUTTA_BOUND) = 2.0
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
! end if

  if (ltrace2) then
    do m=tubi,tuei
       write(0,*) ' u0   - C ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) ' u2   - C ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) ' urk1 - C ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if

  if (use_mask .eq. 1) then
    if (ltrace) write(0,*)'>>> RK3_himem: use_mask'
    call mask_gfuncs(u2, w, par)
 end if

  if (ltrace2) then
    do m=tubi,tuei
       write(0,*) ' u0   - D ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) ' u2   - D ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) ' urk1 - D ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) ''
  end if

#ifdef HYPER_FLUID
 call SolvePrimVars(v, w, u2, par)
#endif

 if (ltrace) write(0,*)'>>> RK3_himem: Done step 1'


  !######### 2nd RK step #############
  else if (iter .eq. 2) then
  if (ltrace) then
     write(0,*)'>>> RK3_himem: begin second step'
  end if

  call calcrhs2(   urk2,  u2,   v,    w,   dxu,  dyu,  dzu,  dxv, dyv, dzv, &
                   imask2,     par)

  if (dissipation > 0) then
    if (ltrace) write(0,*)'>>> RK3_himem: hyperdissipation()', sigma_diss
    do m=1,NU
      if (u2(m)%dissipation .eq. 1) then
        call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,   u2(m)%d,    &
                              w,         imask,     par)
!        urk2(m)%d = urk2(m)%d + W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss +  dyu_diss +  dzu_diss)
        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          urk2(m)%d(i,j,k) = urk2(m)%d(i,j,k)+ &
     &       W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
        end do
        end do
        end do
      end if
    end do
  end if
  par(P_RUNGE_KUTTA_BOUND) = 1.0
  call hyperboundary(u2, dxu, dyu, dzu, urk2, v, w, imask, par)

  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
    !if (amrbound_prepost.eq.1          &
    !    .or. ( nint(chr(i,j,k)) .ne. 1 .and.  &
    !           nint(chr(i,j,k)) .ne. 8      )) then
    if (       (.not. no_bc .or. (no_bc .and. nint(chr(i,j,k)).eq.0))        &
         .and. (      amrbound_prepost.eq.1                                  &
                .or.  ( nint(chr(i,j,k)) .ne. 1 .and. nint(chr(i,j,k)).ne.8))&
        ) then
      do m=1, NU
            u2(m)%d(i,j,k) =u0(m)%d(i,j,k) + dt*0.25*( urk2(m)%d(i,j,k) &
                                           + urk1(m)%d(i,j,k) )
      end do
    else if (nint(chr(i,j,k)) .eq. 1 .and. amrbound_prepost.eq.0) then
      ! Need to interpolate on the boundaries so that
      ! when complete step taken, u2(i,j,k) has the values
      ! given from the parents:
      !    Since we ve already overwritten the original u2,
      ! we need to reconstruct it.
      ! We need to stick into u2() the value:
      !             (3/4) u2tilde + (1/4) u0
      ! where u2tilde() is the original value of u2.
      ! Since we overwrote u2() with:
      !         u2 = (1/2) u2tilde + (1/2) u0
      ! so
      !         u2tilde = 2 u2 - u0
      ! so that
      !         u2 = (3/2) u2tilde - (1/2) u0
      !
      do m=1, NU
         u2(m)%d(i,j,k) = 1.5d0*u2(m)%d(i,j,k) - 0.5d0*u0(m)%d(i,j,k)
      end do
    end if
  end do
  end do
  end do


!  boundary_conditions .eq. 9  ---> PWN2D_boundary
!  boundary_conditions .eq. 10 ---> PWN3D_boundary


!   if ( (boundary_conditions .eq. 8).OR. & 
!        (boundary_conditions .eq. 9).OR. &
!        (boundary_conditions .eq. 10) )then
      par(P_RUNGE_KUTTA_BOUND) = 2.0
     call hyperboundary(u2, dxu, dyu, dzu, urk2, v, w, imask, par)
!   end if

  if (use_mask .eq. 1) then
    call mask_gfuncs(u2, w, par)
  end if

#ifdef HYPER_FLUID
  call SolvePrimVars(v, w, u2, par)
#endif

  !######### 3rd RK step ##############
  else if (iter .eq. 3) then
  if (ltrace) then
     write(0,*)'>>> RK3_himem: begin third step'
  end if

  call calcrhs2(   urk3,  u2,   v,     w,    dxu,  dyu,  dzu,  dxv, dyv, dzv, &
                   imask2,      par)

  if (dissipation > 0) then
    if (ltrace) write(0,*)'>>> RK3_himem: hyperdissipation()', sigma_diss
    do m=1,NU
      if (u2(m)%dissipation .eq. 1) then
        call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,   u2(m)%d,    &
                              w,         imask,     par)
!        urk3(m)%d = urk3(m)%d + W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss +  dyu_diss +  dzu_diss)
        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          urk3(m)%d(i,j,k) = urk3(m)%d(i,j,k)+ &
     &       W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
        end do
        end do
        end do
      end if
    end do
  end if

  par(P_RUNGE_KUTTA_BOUND) = 1.0
  call hyperboundary(u2, dxu, dyu, dzu, urk3, v, w, imask, par)

  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
    !if (amrbound_prepost.eq.1          &
    !    .or. ( nint(chr(i,j,k)) .ne. 1 .and.  &
    !           nint(chr(i,j,k)) .ne. 8      )) then
    if (       (.not. no_bc .or. (no_bc .and. nint(chr(i,j,k)).eq.0))        &
         .and. (      amrbound_prepost.eq.1                                  &
                .or.  ( nint(chr(i,j,k)) .ne. 1 .and. nint(chr(i,j,k)).ne.8))&
        ) then
      do m=1, NU
            u2(m)%d(i,j,k) =u0(m)%d(i,j,k) + (dt/6.0) * ( urk1(m)%d(i,j,k)  &
                           + urk2(m)%d(i,j,k) + 4.0*urk3(m)%d(i,j,k)  )
      end do
    else if (nint(chr(i,j,k)) .eq. 1 .and. amrbound_prepost.eq.0) then
      ! Need to interpolate on the boundaries so that
      ! when complete step taken, u2(i,j,k) has the values
      ! given from the parents:
      !    Interpolation not straightforward, see comments above.
      !
      do m=1, NU
         u2(m)%d(i,j,k) = (4.d0*u2(m)%d(i,j,k) - u0(m)%d(i,j,k) ) /3.d0
      end do
    end if
  end do
  end do
  end do

!  boundary_conditions .eq. 9  ---> PWN2D_boundary
!  boundary_conditions .eq. 10 ---> PWN3D_boundary

! if ( (boundary_conditions .eq. 8).OR. & 
!      (boundary_conditions .eq. 9).OR. &
!      (boundary_conditions .eq. 10) )then
     par(P_RUNGE_KUTTA_BOUND) = 2.0
     call hyperboundary(u2, dxu, dyu, dzu, urk3, v, w, imask, par)
! end if

  if (use_mask .eq. 1) then
    call mask_gfuncs(u2, w, par)
  end if

#ifdef HYPER_FLUID
    call SolvePrimVars(v, w, u2, par)
#endif

! call HyperAnalysis(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, &
!    urk2, urk3, par)

  end if !---iter

  if (allocated(dxu_diss) .and. allocated(dyu_diss) .and. allocated(dzu_diss))&
  then
    if (ltrace2) write(0,*)'>>> RK3_himem: about to deallocate diss vars...'
    deallocate(dzu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'>>> RK3_himem: Could not deallocate dzu_diss: istat=', istat
    end if
    deallocate(dxu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'>>> RK3_himem: Could not deallocate dxu_diss: istat=', istat
    end if
    deallocate(dyu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'>>> RK3_himem: Could not deallocate dyu_diss: istat=', istat
    end if
    if (ltrace2) write(0,*)'>>> RK3_himem: ...diss vars deallocated'
  else
    write(0,*)'>>> RK3_himem: diss vars not allocated as expected!'
  end if


  if (use_mask .eq. 1) then
    !w(H_MASK)%d = real(imask)
    if (allocated(imask)) then
      deallocate(imask,imask2)
    else
      write(0,*)'>>> RK3_himem: imask not allocated as expected!'
    end if
  end if

  if (ltrace) then
    do i = 1, nx
    do j = 1, ny
    do k = 1, nz
      if (abs(w(H_MASK)%d(i,j,k)) .lt. 1.0e-6) then
        print *,'RK3_himemTVD: mask problem, end rk3tvd'
        stop
      end if
    end do
    end do
    end do
  end if


  if (ltrace) then
     write(0,*)'>>>RK3_himem: return'
     if (ltrace2) then
        do m=tubi,tuei
           write(0,*) ' u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
           write(0,*) ' u0   - ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
           write(0,*) ' urk1 - ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
           write(0,*) ' urk2 - ',m,myl2norm3d(urk2(m)%d, nx, ny, nz)
           write(0,*) ' urk3 - ',m,myl2norm3d(urk3(m)%d, nx, ny, nz)
           write(0,*) ' '
        end do
     end if
  end if

  return
end subroutine TVDRK3_himem
