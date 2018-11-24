!-------------------------------------------------------------------------
!
!  $Id$
!
!-------------------------------------------------------------------------

#include "cctk.h"


subroutine TVDRK3(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, par)
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
                                      urk1,dxu,dyu,dzu

  type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
  type(gridfunction), dimension(NW):: w
  CCTK_REAL, dimension(NPAR)           :: par

  CCTK_INT :: iter
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
  logical, parameter                   :: do_check_rk_vars = .false.
  character*(64)                       :: string
  CCTK_INT tubi, tuei, tvbi, tvei
  CCTK_REAL :: myl2norm3d
  CCTK_REAL :: one_third,two_third
  CCTK_INT  :: eos_type
  CCTK_INT  :: myid, proc_return_myid
  CCTK_INT, save :: first_call = 1
  logical    level_grid_exists_within
  external   level_grid_exists_within

  myid = proc_return_myid()
  one_third = 1.0d0/3.0d0
  two_third = 2.0d0/3.0d0
  
  
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
                                                                                                                        
  !hi   = xx(2,1,1)-xx(1,1,1)
  !
  ! Need minimums for this grid alone:
  !
  !minx = xx(1,1,1)
  !miny = yy(1,1,1)
  !minz = zz(1,1,1)
  sigma =0. ! move this later 

  if (ltrace) then
     write(0,*)myid,'] >>>  RK3 begins'
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
    write(0,*)'>>> RK3: Could not allocate memory for dzu_diss'
    write(0,*)'     size : ',nx,ny,nz
  end if
  allocate(dyu_diss(nx,ny,nz), STAT=istat)
  if (istat .ne. 0) then
    write(0,*)'>>> RK3: Could not allocate memory for dyu_diss'
    write(0,*)'     size : ',nx,ny,nz
  end if
  allocate(dxu_diss(nx,ny,nz), STAT=istat)
  if (istat .ne. 0) then
    write(0,*)'>>> RK3: Could not allocate memory for dxu_diss'
    write(0,*)'     size : ',nx,ny,nz
  end if

  if (P_ENABLE_ALT_COORDS==1) then
    if (ltrace) then
      write(0,*)myid,'] >>> RK3: Calling def_alt_coords'
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
      write(0,*)'RK3TVD: Could not allocate memory for imask array'
      write(0,*)'     size : ',nx,ny,nz
    end if
    allocate(imask2(nx,ny,nz),STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'RK3TVD: Could not allocate memory for imask2 array'
      write(0,*)'     size : ',nx,ny,nz
    end if


    call set_mask(w, imask,imask2, par)

  else
    allocate(imask(nx,ny,nz),STAT=istat)
    allocate(imask2(nx,ny,nz),STAT=istat)
    imask = 1
    imask2 = 1
  end if

  if (do_check_rk_vars) then
    write(string,'(A)') 'BEGIN OF RK3'
    call check_rk_vars(u0, v, w, par, string)
  end if

  !###################################
  !######### 1st RK step #############
  !###################################

  if (iter .eq. 1) then
  if (ltrace) then
     write(0,*)myid,'] >>> RK3: begin first step'
     write(0,*)myid,'] >>> RK3: dt,  dx = ',dt, dx
  end if

  if (ltrace2) then
     write(0,*)myid,'] >>>RK3: initial data'
     do m=tubi,tuei
        write(0,*) myid,'] rk3: ||u||  ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
     end do
     do m=tvbi,tvei
        write(0,*) myid,'] rk3: ||v||  ',m,myl2norm3d(v(m)%d, nx, ny, nz)
     end do
  end if

#ifdef HYPER_FLUID
  if (first_call .eq. 1) then
    first_call = 0
    eos_type = nint(par(P_EOS_TYPE))
    if (eos_type .lt. 100) then
      ! Ideal gas EOS or UWM EOS
      if (have_eos_table .eq. 0) then
        call qsetgeneos(par)
      end if
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

  if ( par(P_ALLOWEDL) .gt. 0 ) then
    ! Ensure that primitive variables are never interpolated
    ! Problem interpolation will only occur when there are multiple levels
    if (ltrace) write(0,*)myid,'] >>>RK3: Calling SolvePrimVars'
    call SolvePrimVars(v, w, u0, par)
    if (ltrace) write(0,*)myid,'] >>>RK3: Back from SolvePrimVars'
  end if
#endif

  call calcrhs2(   urk1,  u0,   v,    w,  dxu,  dyu,  dzu,  dxv, dyv,&
       & dzv, imask2, par)

  if (ltrace2) then
    write(0,*) myid,']  After calling the first rhs '
    do m=tubi,tuei    
!      write(0,*) ' u0   - A ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
!      write(0,*) ' urk1 - A ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
       write(0,25) myid, m, myl2norm3d(u0(m)%d, nx, ny, nz), &
                  &   myl2norm3d(urk1(m)%d, nx, ny, nz)
 25    format(I3,'] RK3 A:  ',i4,' ||u|| = ',es9.2,' ||dtu|| = ',es9.2)
    end do
    write(0,*) myid,'] '
  end if

  if (dissipation > 0) then
    if (ltrace) write(0,*)myid,'] >>> RK3: hyperdissipation()', sigma_diss
    do m = 1, NU
      if (u0(m)%dissipation .eq. 1) then
        call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,   u0(m)%d&
             &, w,         imask,     par)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          urk1(m)%d(i,j,k) = urk1(m)%d(i,j,k)+                   &
     &      W(H_WDISS)%d(i,j,k)*sigma_diss*u0(m)%diss_factor*    &
     &      (dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
        end do
        end do
        end do
      end if
    end do
  end if

  par(P_RUNGE_KUTTA_BOUND) = 1.0
  call hyperboundary(u0, dxu, dyu, dzu, urk1, v, w, imask, par)

  if (ltrace2) then
    write(0,*) myid,']  After calling hyperboundary'
    do m=tubi,tuei
       write(0,28) myid, m, myl2norm3d(u0(m)%d, nx, ny, nz), &
         & myl2norm3d(urk1(m)%d, nx, ny, nz)
 28    format(I3,'] RK3 B1:  ',i4,' ||u0|| = ',es9.2,' ||urk1| = ',es9.2)
    end do
    write(0,*) myid,'] '
  end if

  if (ltrace) write(0,*)myid,'] >>> RK3: update u2()', NU,nx,ny,nz
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
    if (       (.not. no_bc .or. (no_bc .and. nint(chr(i,j,k)).eq.0))        &
         .and. (      amrbound_prepost.eq.1                                  &
                .or.  ( nint(chr(i,j,k)) .ne. 1 .and. nint(chr(i,j,k)).ne.8))&
        ) then
      do m=1, NU
        u2(m)%d(i,j,k) = u0(m)%d(i,j,k) + dt*urk1(m)%d(i,j,k)
      end do
    end if
  end do
  end do
  end do

  if (ltrace2) then
    write(0,*) myid,']  After updating the first rk '
    do m=tubi,tuei
!      write(0,*) ' u0   - B ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
!      write(0,*) ' u2   - B ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
!      write(0,*) ' urk1 - B ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
       write(0,30) m, myl2norm3d(u0(m)%d, nx, ny, nz), &
         & myl2norm3d(u2(m)%d, nx, ny, nz), myl2norm3d(urk1(m)%d, nx, ny, nz)
 30    format('RK3 B2:  ',i4,' ||u0|| = ',es9.2,' ||u2|| = ',es9.2, &
         &   ' ||urk1|| = ',es9.2)
    end do
    write(0,*) ' '
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
    write(0,*) myid,']  After applying the first boundary condition '
    do m=tubi,tuei
!      write(0,*) ' u0   - C ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
!      write(0,*) ' u2   - C ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
!      write(0,*) ' urk1 - C ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
       write(0,35) m, myl2norm3d(u0(m)%d, nx, ny, nz), &
        &  myl2norm3d(u2(m)%d, nx, ny, nz), myl2norm3d(urk1(m)%d, nx, ny, nz)
 35    format('RK3 C:  ',i4,' ||u0|| = ',es9.2,' ||u2|| = ',es9.2, &
         &   ' ||urk1|| = ',es9.2)
    end do
    write(0,*) ' '
  end if

  if (use_mask .eq. 1) then
    if (ltrace) write(0,*)'>>> RK3: use_mask'
    call mask_gfuncs(u2, w, par)
 end if

  if (ltrace2) then
    write(0,*) myid,']  After applying the first mask '  
    do m=tubi,tuei
!      write(0,*) ' u0   - D ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
!      write(0,*) ' u2   - D ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
!      write(0,*) ' urk1 - D ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
       write(0,40) m, myl2norm3d(u0(m)%d, nx, ny, nz), &
        &  myl2norm3d(u2(m)%d, nx, ny, nz), myl2norm3d(urk1(m)%d, nx, ny, nz)
 40    format('RK3 D:  ',i4,' ||u0|| = ',es9.2,' ||u2|| = ',es9.2, &
         &   ' ||urk1|| = ',es9.2)
    end do
    write(0,*) myid,'] '
  end if


#ifdef HYPER_BSSN
  if (ltrace2) then
    write(0,*)myid,']  Before EnforceBSSN'
    do m = tubi, tuei
       write(0,*) myid,']  u2 D1 - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
    end do
  end if

  call EnforceBSSN(u2, v, w, par)

  if (ltrace2) then
    write(0,*)myid,']  After EnforceBSSN'
    do m = tubi, tuei
       write(0,*) myid,']  u2 D2 - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
    end do
  end if
#endif



#ifdef HYPER_FLUID
 call SolvePrimVars(v, w, u2, par)
#endif

 ! periodic boundary conditions in a single processor 
 if (bc_type .eq. 3) then
    if (ltrace) then
      write(0,*)myid,'] >>> RK3: Calling periodic BC'
    end if
    do m=1, NU
       call copy_bc2(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
    end do
 end if


 if (ltrace2) then
   write(0,*)myid,']  End of Step 1' 
   do m = tubi, tuei
      write(0,*) myid,']  u2 D3 - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
   end do
 end if

 if (ltrace) write(0,*)myid,'] >>> RK3: Done step 1'


  !######### 2nd RK step #############
  else if (iter .eq. 2) then
  if (ltrace) then
     write(0,*)myid,'] >>> RK3: begin second step'
  end if

  call calcrhs2(   urk1,  u2,   v,    w,   dxu,  dyu,  dzu,  dxv, dyv, dzv, &
                   imask2,     par)

  if (ltrace2) then
    write(0,*) myid,']  After calling the second rhs '
    do m=tubi,tuei    
       write(0,*) myid,']  urk2 - E ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) myid,'] '
  end if
  
  if (dissipation > 0) then
    if (ltrace) write(0,*)myid,'] >>> RK3: hyperdissipation()', sigma_diss
    do m=1,NU
      if (u2(m)%dissipation .eq. 1) then
        call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,   u2(m)%d,    &
                              w,         imask,     par)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          urk1(m)%d(i,j,k) = urk1(m)%d(i,j,k)+                  &
     &       W(H_WDISS)%d(i,j,k)*sigma_diss*u2(m)%diss_factor*  &
     &       (dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
        end do
        end do
        end do
      end if
    end do
  end if
  par(P_RUNGE_KUTTA_BOUND) = 1.0
  call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)

  !if (myid.eq.6 .and.level_grid_exists_within(19,2)) then
     !write(*,*) myid,'TVDRK3 L: calling grid_test(161)'
     !call grid_test(161)
  !end if

  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
    if (       (.not. no_bc .or. (no_bc .and. nint(chr(i,j,k)).eq.0))        &
         .and. (      amrbound_prepost.eq.1                                  &
                .or.  ( nint(chr(i,j,k)) .ne. 1 .and. nint(chr(i,j,k)).ne.8))&
        ) then
      do m=1, NU
            u2(m)%d(i,j,k) = 0.75* u0(m)%d(i,j,k) + 0.25 * u2(m)%d(i,j,k) &
                                                  + 0.25 * dt * urk1(m)%d(i,j,k)
      end do
    end if
  end do
  end do
  end do

  if (ltrace2) then
    write(0,*) myid,']  After updating the second rk '
    do m=tubi,tuei
       write(0,*) myid,']  u2   - F ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) myid,']  urk2 - F ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) myid,'] '
  end if

!  boundary_conditions .eq. 9  ---> PWN2D_boundary
!  boundary_conditions .eq. 10 ---> PWN3D_boundary


!   if ( (boundary_conditions .eq. 8).OR. & 
!        (boundary_conditions .eq. 9).OR. &
!        (boundary_conditions .eq. 10) )then
      par(P_RUNGE_KUTTA_BOUND) = 2.0
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
!   end if

  if (use_mask .eq. 1) then
    call mask_gfuncs(u2, w, par)
  end if

  if (ltrace) then
     write(0,*)myid,'] >>>RK3: G'
     if (ltrace2) then
        do m=tubi,tuei
           write(0,*) myid,']  u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
        end do
     end if
  end if
#ifdef HYPER_BSSN
  call EnforceBSSN(u2, v, w, par)
#endif

  if (ltrace) then
     write(0,*)myid,'] >>>RK3: H'
     if (ltrace2) then
        do m=tubi,tuei
           write(0,*) myid,']  u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
        end do
     end if
  end if

#ifdef HYPER_FLUID
  call SolvePrimVars(v, w, u2, par)
#endif

  if (bc_type .eq. 3) then
      do m=1, NU
         call copy_bc2(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      end do
  end if

  if (ltrace) then
     write(0,*)myid,'] >>>RK3: I'
     if (ltrace2) then
        do m=tubi,tuei
           write(0,*) myid,']  u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
        end do
     end if
  end if

  !######### 3rd RK step ##############
  else if (iter .eq. 3) then
  if (ltrace) then
     write(0,*)myid,'] >>> RK3: begin third step'
  end if

  call calcrhs2(   urk1,  u2,   v,     w,    dxu,  dyu,  dzu,  dxv, dyv, dzv, &
                   imask2,      par)

  if (ltrace2) then
    write(0,*) myid,']  After calling the third rhs '
    do m=tubi,tuei    
       write(0,*) myid,']  urk3 - K ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) myid,'] '
  end if
  
  if (dissipation > 0) then
    if (ltrace) write(0,*)myid,'] >>> RK3: hyperdissipation()', sigma_diss
    do m=1,NU
      if (u2(m)%dissipation .eq. 1) then
        call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,   u2(m)%d,    &
                              w,         imask,     par)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          urk1(m)%d(i,j,k) = urk1(m)%d(i,j,k)+                    &
     &       W(H_WDISS)%d(i,j,k)*sigma_diss*u2(m)%diss_factor*    &
     &       (dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
        end do
        end do
        end do
      end if
    end do
  end if

  par(P_RUNGE_KUTTA_BOUND) = 1.0
  call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)

  if (ltrace2) then
    write(0,*) myid,']  After applying the third boundary condition '
    do m=tubi,tuei
       write(0,*) myid,']  u2   - L ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) myid,']  urk3 - L ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) myid,'] '
  end if
  
  
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
    if (       (.not. no_bc .or. (no_bc .and. nint(chr(i,j,k)).eq.0))        &
         .and. (      amrbound_prepost.eq.1                                  &
                .or.  ( nint(chr(i,j,k)) .ne. 1 .and. nint(chr(i,j,k)).ne.8))&
        ) then
      do m=1, NU
            u2(m)%d(i,j,k) = one_third * u0(m)%d(i,j,k) &
                                        + two_third * u2(m)%d(i,j,k) &
                                        + two_third * dt * urk1(m)%d(i,j,k)      
      end do
    end if
  end do
  end do
  end do


  if (ltrace2) then
    write(0,*) myid,']  After updating the third rk '
    do m=tubi,tuei
       write(0,*) myid,']  u2   - M ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) myid,']  urk3 - M ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) myid,'] '
  end if
  
  
!  boundary_conditions .eq. 9  ---> PWN2D_boundary
!  boundary_conditions .eq. 10 ---> PWN3D_boundary

! if ( (boundary_conditions .eq. 8).OR. & 
!      (boundary_conditions .eq. 9).OR. &
!      (boundary_conditions .eq. 10) )then
     par(P_RUNGE_KUTTA_BOUND) = 2.0
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
! end if

  if (use_mask .eq. 1) then
    call mask_gfuncs(u2, w, par)
  end if

  if (ltrace) then
     write(0,*)myid,'] >>>RK3: U'
     if (ltrace2) then
        do m=tubi,tuei
           write(0,*) myid,']  u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
        end do
     end if
  end if


#ifdef HYPER_BSSN
  call EnforceBSSN(u2, v, w, par)
#endif
  if (ltrace) then
     write(0,*)myid,'] >>>RK3: V'
     if (ltrace2) then
        do m=tubi,tuei
           write(0,*) myid,']  u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
        end do
     end if
  end if


#ifdef HYPER_FLUID
    call SolvePrimVars(v, w, u2, par)
#endif
  if (ltrace) then
     write(0,*)myid,'] >>>RK3: W'
     if (ltrace2) then
        do m=tubi,tuei
           write(0,*) myid,']  u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
        end do
     end if
  end if

  if (bc_type .eq. 3) then
      do m=1, NU
         call copy_bc2(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      end do
  end if

  end if !---iter

  if (do_check_rk_vars) then
    write(string,'(A)') 'END OF RK3'
    call check_rk_vars(u2, v, w, par, string)
  end if


  if (allocated(dxu_diss) .and. allocated(dyu_diss) .and. allocated(dzu_diss))&
  then
    if (ltrace2) write(0,*)myid,'] >>> RK3: about to deallocate diss vars...'
    deallocate(dzu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)myid,'] >>> RK3: Could not deallocate dzu_diss: istat=', istat
    end if
    deallocate(dxu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)myid,'] >>> RK3: Could not deallocate dxu_diss: istat=', istat
    end if
    deallocate(dyu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)myid,'] >>> RK3: Could not deallocate dyu_diss: istat=', istat
    end if
    if (ltrace2) write(0,*)myid,'] >>> RK3: ...diss vars deallocated'
  else
    write(0,*)myid,'] >>> RK3: diss vars not allocated as expected!'
  end if


  if (use_mask .eq. 1) then
    !w(H_MASK)%d = real(imask)
    if (allocated(imask)) then
      deallocate(imask,imask2)
    else
      write(0,*)myid,'] >>> RK3: imask not allocated as expected!'
    end if
  else 
      deallocate(imask,imask2)
  end if

  if (ltrace.and. .true.) then
    do i = 1, nx
    do j = 1, ny
    do k = 1, nz
      if (abs(w(H_MASK)%d(i,j,k)) .lt. 1.0e-6) then
        print *,myid,'] RK3TVD: mask problem, end rk3tvd'
        call my_exit('RK3TVD:  mask problem')
        !stop
      end if
    end do
    end do
    end do
  end if

  if (ltrace) then
     write(0,*)myid,'] >>>RK3: return'
     if (ltrace2) then
        do m=tubi,tuei
           write(0,*) myid,']  u2   - ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
           write(0,*) myid,']  u0   - ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
           write(0,*) myid,']  urk3 - ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
           write(0,*) myid,']  '
        end do
     end if
  end if

  return
end subroutine TVDRK3

!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
subroutine copy_bc2(q, chr, numx, numy, numz, xlwb, xupb, &
                              ylwb, yupb, zlwb, zupb )
implicit  none

! subroutine arguments
 CCTK_INT, intent (in) :: numx
 CCTK_INT, intent (in) :: numy
 CCTK_INT, intent (in) :: numz

 CCTK_REAL, dimension(numx,numy,numz), intent (inout) :: q
 CCTK_REAL, dimension(numx,numy,numz), intent (in)    :: chr

 CCTK_INT, intent (in) :: xlwb
 CCTK_INT, intent (in) :: xupb
 CCTK_INT, intent (in) :: ylwb
 CCTK_INT, intent (in) :: yupb
 CCTK_INT, intent (in) :: zlwb
 CCTK_INT, intent (in) :: zupb

!
!  impose boundary conditions
!
!
 CCTK_INT :: I, J, K

logical              :: ltrace = .false.

!return

if (ltrace) then
   write(*,*) 'copy_bc: numx/y/z = ',numx,numy,numz
end if

do K = 1, numz
   do J = 1, numy
      if (nint(chr(1, J,K)).ne.1) &
!       q(1,J,K) = q(numx-5,J,K)      
!       q(2,J,K) = q(numx-4,J,K)      
!       q(3,J,K) = q(numx-3,J,K)      
       q(1,J,K) = q(numx-8,J,K)      
       q(2,J,K) = q(numx-7,J,K)      
       q(3,J,K) = q(numx-6,J,K)      
       q(4,J,K) = q(numx-5,J,K)      
      if (nint(chr(numx,J,K)).ne.1) &
!       q(numx-2,J,K) = q(4,J,K)      
!       q(numx-1,J,K) = q(5,J,K)      
!       q(numx,J,K)   = q(6,J,K)      
       q(numx-4,J,K) = q(5,J,K)      
       q(numx-3,J,K) = q(6,J,K)      
       q(numx-2,J,K) = q(7,J,K)      
       q(numx-1,J,K) = q(8,J,K)      
       q(numx,J,K)   = q(9,J,K)          
   enddo
enddo

do K = 1, numz
   do I = 1, numx
      if (nint(chr(I,1,   K)).ne.1) &
!       q(I,1,K) = q(I,numy-5,K)
!       q(I,2,K) = q(I,numy-4,K)
!       q(I,3,K) = q(I,numy-3,K)
       q(I,1,K) = q(I,numy-8,K)      
       q(I,2,K) = q(I,numy-7,K)      
       q(I,3,K) = q(I,numy-6,K)      
       q(I,4,K) = q(I,numy-5,K)
      if (nint(chr(I,numy,K)).ne.1) &
!       q(I,numy-2,K) = q(I,4,K)
!       q(I,numy-1,K) = q(I,5,K)
!       q(I,numy,K)   = q(I,6,K)
       q(I,numy-4,K) = q(I,5,K)      
       q(I,numy-3,K) = q(I,6,K)      
       q(I,numy-2,K) = q(I,7,K)      
       q(I,numy-1,K) = q(I,8,K)      
       q(I,numy,K)   = q(I,9,K)      
   enddo
enddo

do J = 1, numy
   do I = 1, numx
      if (nint(chr(I,J,1   )).ne.1) &
!       q(I,J,1) = q(I,J,numz-5 )
!       q(I,J,2) = q(I,J,numz-4 )
!       q(I,J,3) = q(I,J,numz-3 )
       q(I,J,1) = q(I,J,numz-8)      
       q(I,J,2) = q(I,J,numz-7)      
       q(I,J,3) = q(I,J,numz-6)      
       q(I,J,4) = q(I,J,numz-5)
      if (nint(chr(I,J,numz)).ne.1) &
!       q(I,J,numz-2) = q(I,J,4)      
!       q(I,J,numz-1) = q(I,J,5)      
!       q(I,J,numz)   = q(I,J,6)      
       q(I,J,numz-4) = q(I,J,5)      
       q(I,J,numz-3) = q(I,J,6)      
       q(I,J,numz-2) = q(I,J,7)      
       q(I,J,numz-1) = q(I,J,8)      
       q(I,J,numz)   = q(I,J,9)      
   enddo
enddo

end subroutine copy_bc2
!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
subroutine check_rk_vars(u, v, w, par, string)
  use params
  use GF
  implicit none
  type(gridfunction), dimension(NU):: u
  type(gridfunction), dimension(NV):: v
  type(gridfunction), dimension(NW):: w
  CCTK_REAL, dimension(NPAR)       :: par
  character*(*)                    :: string

  ! local vars
  CCTK_INT                         :: nx, ny, nz
  CCTK_INT                         :: i, j, k, m, rc

  nx                    = nint(par(P_NX))
  ny                    = nint(par(P_NY))
  nz                    = nint(par(P_NZ))

  do m = 1, NU
    call check_func(u(m)%d, nx, ny, nz, rc)
    if (rc .eq. -1) then
      write(0,*)'Problem in U var, m = ',m,' at ',string
    end if
  end do
  do m = 1, NV
    call check_func(v(m)%d, nx, ny, nz, rc)
    if (rc .eq. -1) then
      write(0,*)'Problem in V var, m = ',m,' at ',string
    end if
  end do

  return
end subroutine check_rk_vars

!------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------
subroutine check_func(f, nx, ny, nz, rc)
  implicit none
  CCTK_INT                         :: nx, ny, nz, rc
  CCTK_REAL, dimension(nx,ny,nz)   :: f

  ! local vars
  CCTK_INT                         :: i, j, k

  rc = 1

  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
    if (abs(f(i,j,k)) .gt. 1.0d6) then
      rc = -1
    end if
  end do
  end do
  end do

  return
end subroutine check_func

