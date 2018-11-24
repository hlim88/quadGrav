!-------------------------------------------------------------------------
!
!  $Id$
!
!
!
!-------------------------------------------------------------------------

#include "cctk.h"


subroutine RK3(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, urk1, par)
  use params
  use UTILEQS
  use GF
  use HYPER_RHS
  use HYPER_DISSIPATION
  use HYPER_BOUNDARY
  use MOD_MASK
  use hypercoords


  implicit none

  CCTK_INT     istat, rc, nt, nx, ny, nz, err, m, use_mask, dissipation
  CCTK_INT     update_scheme, amrbound_prepost
  CCTK_INT     boundary_conditions
  type(gridfunction), dimension(NU):: u0,u2,urk1,&
                                      dxu,dyu,dzu

  type(gridfunction), dimension(NV):: v, dxv, dyv, dzv
  type(gridfunction), dimension(NW):: w
  CCTK_REAL, dimension(NPAR)           :: par

  CCTK_INT ::   iter, imask_pt
  CCTK_REAL :: dt, dx, sigma, sigma_diss, normv
  CCTK_INT, allocatable, dimension(:,:,:) :: imask, imask2
  CCTK_REAL, allocatable, dimension(:,:,:):: dxu_diss, dyu_diss, dzu_diss
  logical, parameter                   :: ltrace  = .false.
  !logical                              :: ltrace  
  logical, parameter                   :: ltrace2 = .false.
  ! Trace boundary treatment:
  logical, parameter                   :: ltraceB = .false.
  ! For debugging, switch to turn off all BC treatment:
  logical, parameter                   :: no_bc   = .false.
  CCTK_INT CCTK_Equals
  CCTK_INT:: handle, err_red, i, j, k, kk, order
  CCTK_REAL:: norm_rk, max_g11, min_g11
  CCTK_REAL, pointer, dimension(:,:,:) :: chr
  CCTK_REAL, dimension(:,:,:),pointer  :: xx, yy, zz
  CCTK_INT:: bc_type, UNDEF, DEF

  CCTK_INT:: proc_return_myid, myid
  CCTK_REAL :: myl2norm3d, time
  CCTK_REAL :: one_third,two_third
  CCTK_INT :: rcc
  CCTK_INT :: assume_symmetry
  logical :: double_equal
  external   double_equal
  

  UNDEF     = P_STENCIL_X_UNDEF + P_STENCIL_Y_UNDEF + P_STENCIL_Z_UNDEF
  DEF       = P_STENCIL_CENTER

  one_third = 1.0d0/3.0d0
  two_third = 2.0d0/3.0d0

  iter = nint(par(P_RK_ITER))
  assume_symmetry = nint(par(P_ASSUME_SYMMETRY))
  chr => w(H_CHR)%d
  xx  => w(H_XPHYS)%d
  yy  => w(H_YPHYS)%d
  zz  => w(H_ZPHYS)%d
  nx                    = nint(par(P_NX))
  ny                    = nint(par(P_NY))
  nz                    = nint(par(P_NZ))
                                                                                                                        
  !ltrace = .false.
  myid = proc_return_myid()
  !if (myid.eq.1.and.nx.eq.39.and.ny.eq.41.and.nz.eq.39.and.double_equal(xx(1,1,1),-5.0).and.double_equal(yy(1,1,1),-1.66666666666667)) ltrace=.true.
  sigma =0. ! move this later 

  time = par(P_TIME)
  if (ltrace.or.ltrace2.or.ltraceB) then
     write(0,*)myid, '] >>>  RK3 begins'
     !call output3DA(u2(H_B1)%d,"h_b1pre",xx,yy,zz,time,myid)
     !write(0,*)myid, '] >>>  0 u2(H_B1)%d(31,32,14)=',u2(H_B1)%d(31,32,14)
  end if
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


  !define coordinates if not using the default cartesian computational coords.
  if (P_ENABLE_ALT_COORDS==1) then
    if (ltrace) then
      write(0,*)myid, '] >>> RK3: Calling def_alt_coords'
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
      write(0,*)'RK3: Could not allocate memory for imask array'
      write(0,*)'     size : ',nx,ny,nz
    end if

    allocate(imask2(nx,ny,nz),STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'RK3: Could not allocate memory for imask2 array'
      write(0,*)'     size : ',nx,ny,nz
    end if

      call Set_mask(w, imask,imask2, par)
!     myid = proc_return_myid()
!      write(0,*) 'about to in rk', myid
!      rcc = check_imask(imask,par)
!      if(rcc.eq.0) write(0,*) 'ops rk3', myid, rcc

   end if
      
  if (ltrace2) then
     write(0,*)myid, '] >>> RK3----params:'
     write(0,*)myid, ']        iter        = ',iter
     write(0,*)myid, ']        nx, ny, nz  = ',nx, ny, nz
     write(0,*)myid, ']        use_mask    = ',use_mask
     write(0,*)myid, ']        dissipation = ',dissipation
     write(0,*)myid, ']        boundary_conditions = ',boundary_conditions
     write(0,*)myid, ']        bc_type = ',bc_type     
     write(0,*)myid, ']        update_scheme = ', update_scheme
     write(0,*)myid, '] >>> RK3----end params:'
  end if

  if (ltrace2) then
     write(0,*)myid, '] >>>RK3: data upon entry:'
     do m=1,NU
        write(0,*)myid, '] RK3 ID: ||u0||  ',iter,m,myl2norm3d(u0(m)%d, nx, ny, nz)
     end do
     do m=1,NU
        write(0,*)myid, '] RK3 ID: ||u2||  ',iter,m,myl2norm3d(u2(m)%d, nx, ny, nz)
     end do
     do m=1,NV
        write(0,*)myid, '] RK3 ID: ||v||  ',iter,m,myl2norm3d(v(m)%d, nx, ny, nz)
     end do
  end if

  !if (ltrace)write(0,*)myid, '] >>>RK3  0.5   u0(H_B1)%d(31,32,14)=',  u0(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.5   u2(H_B1)%d(31,32,14)=',  u2(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.5 urk1(H_B1)%d(31,32,14)=',urk1(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.5   u0(H_E3)%d(31,32,14)=',  u0(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.5   u2(H_E3)%d(31,32,14)=',  u2(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.5 urk1(H_E3)%d(31,32,14)=',urk1(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.5  dxu(H_B2)%d(31,32,14)=', dxu(H_B2)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.5  dxu(H_B3)%d(31,32,14)=', dxu(H_B3)%d(31,32,14)

  !###################################
  !######### 1st RK step #############
  !###################################

  if (iter .eq. 1) then
  if (ltrace) then
     write(0,*)myid,'] >>> RK3: begin first step'
     write(0,*)myid,'] >>> RK3: dt,  dx = ',dt, dx
  end if

  if (boundary_conditions.eq.1.or.boundary_conditions.eq.2.or.boundary_conditions.eq.7 .and. .not.no_bc) then
     if (ltraceB) write(0,*)myid,'] >>> RK3: hyperboundary()'
     call hyperboundary(u0, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

  !if (ltrace)write(0,*)myid, '] >>>RK3  0.75   u0(H_B1)%d(31,32,14)=',  u0(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.75   u2(H_B1)%d(31,32,14)=',  u2(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.75 urk1(H_B1)%d(31,32,14)=',urk1(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.75   u0(H_E3)%d(31,32,14)=',  u0(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.75   u2(H_E3)%d(31,32,14)=',  u2(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.75 urk1(H_E3)%d(31,32,14)=',urk1(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.75  dxu(H_B2)%d(31,32,14)=', dxu(H_B2)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.75  dxu(H_B3)%d(31,32,14)=', dxu(H_B3)%d(31,32,14)

  if (ltrace) write(0,*)myid,'] >>> RK3: hyperrhs()'
  if (ltrace2) then
    do m=1,NU
       write(0,*) myid, '] u0   - 1 ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
    end do
    do m=1,NV
       write(0,*) myid, '] v    - 1 ',m,myl2norm3d(v(m)%d, nx, ny, nz)
    end do
  end if
  
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.79   u0(H_B1)%d(31,32,14)=',  u0(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.79   u2(H_B1)%d(31,32,14)=',  u2(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.79 urk1(H_B1)%d(31,32,14)=',urk1(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.79   u0(H_E3)%d(31,32,14)=',  u0(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.79   u2(H_E3)%d(31,32,14)=',  u2(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.79 urk1(H_E3)%d(31,32,14)=',urk1(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.79  dxu(H_B2)%d(31,32,14)=', dxu(H_B2)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.79  dxu(H_B3)%d(31,32,14)=', dxu(H_B3)%d(31,32,14)

  call hyperrhs(   urk1,  u0,   v,    w,  dxu,  dyu,  dzu,  dxv, dyv,&
       & dzv, imask2, par)
       
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.85   u0(H_B1)%d(31,32,14)=',  u0(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.85   u2(H_B1)%d(31,32,14)=',  u2(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.85 urk1(H_B1)%d(31,32,14)=',urk1(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.85   u0(H_E3)%d(31,32,14)=',  u0(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.85   u2(H_E3)%d(31,32,14)=',  u2(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.85 urk1(H_E3)%d(31,32,14)=',urk1(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.85  dxu(H_B2)%d(31,32,14)=', dxu(H_B2)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.85  dxu(H_B3)%d(31,32,14)=', dxu(H_B3)%d(31,32,14)

  if (ltrace2) then
    do m=1,NU
       write(0,*) myid, ']  u0   - A ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
       write(0,*) myid, ']  urk1 - A ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) myid,'] '
  end if

  if (dissipation > 0) then
    if (ltrace) write(0,*)myid, '] >>> RK3: hyperdissipation()', sigma_diss,dissipation
    do m=1,NU
      if (dissipation.ne.3) then
      call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,   u0(m)%d&
           &, w,         imask,     par)  
      if(ltrace2)write(0,*)myid, '] >>> RK3: back, now applying...'
      if(ltrace2)write(0,*)myid, '] >>> RK3: adding to urk1()...'
      do k = 1,nz
      do j = 1,ny
      do i = 1,nx
         imask_pt = DEF
         if (use_mask .eq. 1) imask_pt = imask(i,j,k)
         if (imask_pt .ne.UNDEF) then
           if(ltrace2)write(0,*)myid, '] >>> RK3: dx/y/zu_diss',dxu_diss(i,j,k),dyu_diss(i,j,k), dzu_diss(i,j,k)
           urk1(m)%d(i,j,k) = urk1(m)%d(i,j,k)+ &
     &   W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
         end if
      end do
      end do
      end do

      else
      call apply_diss_torhs( urk1(m)%d, u0(m)%d, sigma_diss,nx,ny,nz)
      end if
    end do
  end if

  !if (ltrace)write(0,*)myid, '] >>>RK3  0.87   u0(H_B1)%d(31,32,14)=',  u0(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.87   u2(H_B1)%d(31,32,14)=',  u2(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.87 urk1(H_B1)%d(31,32,14)=',urk1(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.87   u0(H_E3)%d(31,32,14)=',  u0(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.87   u2(H_E3)%d(31,32,14)=',  u2(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.87 urk1(H_E3)%d(31,32,14)=',urk1(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.87  dxu(H_B2)%d(31,32,14)=', dxu(H_B2)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.87  dxu(H_B3)%d(31,32,14)=', dxu(H_B3)%d(31,32,14)

  if (ltrace2) then
     do m=1,NU
        write(0,*)myid,']  urk1 - B ',iter,m,myl2norm3d(urk1(m)%d, nx, ny, nz)
     end do
     !write(0,*) ' '
  end if

  if (boundary_conditions .gt. 2 .and. bc_type .ne. 3 .and. &
      boundary_conditions .ne. 8 .and. .not.no_bc) then
     if (ltraceB) write(0,*)myid,'] >>> RK3: hyperboundary'
     if (ltraceB) write(0,*)myid,'] >>> RK3: dxu = ',dxu(1)%d(1,2,3)
     call hyperboundary(u0, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

  if (ltrace) write(0,*)myid,'] >>> RK3: update u2()',iter
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
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
         if (update_scheme .eq. 3) then
            u2(m)%d(i,j,k) = u0(m)%d(i,j,k) + dt * urk1(m)%d(i,j,k)
         else
            write(0,*) "Error, update_scheme not known"            
         end if
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
  
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.89   u0(H_B1)%d(31,32,14)=',  u0(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.89   u2(H_B1)%d(31,32,14)=',  u2(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.89 urk1(H_B1)%d(31,32,14)=',urk1(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.89   u0(H_E3)%d(31,32,14)=',  u0(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.89   u2(H_E3)%d(31,32,14)=',  u2(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.89 urk1(H_E3)%d(31,32,14)=',urk1(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.89  dxu(H_B2)%d(31,32,14)=', dxu(H_B2)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.89  dxu(H_B3)%d(31,32,14)=', dxu(H_B3)%d(31,32,14)

  if (boundary_conditions .eq. 8) then
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

  !if (ltrace)write(0,*)myid, '] >>>RK3  0.90   u0(H_B1)%d(31,32,14)=',  u0(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.90   u2(H_B1)%d(31,32,14)=',  u2(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.90 urk1(H_B1)%d(31,32,14)=',urk1(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.90   u0(H_E3)%d(31,32,14)=',  u0(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.90   u2(H_E3)%d(31,32,14)=',  u2(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.90 urk1(H_E3)%d(31,32,14)=',urk1(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.90  dxu(H_B2)%d(31,32,14)=', dxu(H_B2)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.90  dxu(H_B3)%d(31,32,14)=', dxu(H_B3)%d(31,32,14)

  if (bc_type.eq.3.and. .not.no_bc) then
     if (ltraceB) write(*,*) myid,'] RK3: Using copy boundary condition:', iter
      do m=1, NU
         call copy_bc(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      end do
      !
      ! This will not compile in other projects, because it refers specifically
      ! to a field which does not exist in them:
      !
      !call jet_bc(u2(H_RHO)%d,    chr,xx,yy,zz,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      !call jet_bc(u2(H_RHO_U_X)%d,chr,xx,yy,zz,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
  end if

  !if (ltrace)write(0,*)myid, '] >>>RK3  0.91   u0(H_B1)%d(31,32,14)=',  u0(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.91   u2(H_B1)%d(31,32,14)=',  u2(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.91 urk1(H_B1)%d(31,32,14)=',urk1(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.91   u0(H_E3)%d(31,32,14)=',  u0(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.91   u2(H_E3)%d(31,32,14)=',  u2(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.91 urk1(H_E3)%d(31,32,14)=',urk1(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.91  dxu(H_B2)%d(31,32,14)=', dxu(H_B2)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.91  dxu(H_B3)%d(31,32,14)=', dxu(H_B3)%d(31,32,14)

  if (assume_symmetry .ne. 0 .and. .false.) then
     !
     ! Enforce symmetry conditions
     !
     if (.true.) write(*,*) myid,'] RK3: Enforcing symmetry: ', iter
      do m=1, NU
         call enforce_symmetry(u2(m)%d,chr,assume_symmetry,nx,ny,nz)
      end do
  end if

  if (use_mask .eq. 1) then
    if (ltrace) write(0,*)myid,'] >>> RK3: use_mask'
    call mask_gfuncs(u2, w, par)
  end if
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.95   u0(H_B1)%d(31,32,14)=',  u0(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.95   u2(H_B1)%d(31,32,14)=',  u2(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.95 urk1(H_B1)%d(31,32,14)=',urk1(H_B1)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.95   u0(H_E3)%d(31,32,14)=',  u0(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.95   u2(H_E3)%d(31,32,14)=',  u2(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.95 urk1(H_E3)%d(31,32,14)=',urk1(H_E3)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.95  dxu(H_B2)%d(31,32,14)=', dxu(H_B2)%d(31,32,14)
  !if (ltrace)write(0,*)myid, '] >>>RK3  0.95  dxu(H_B3)%d(31,32,14)=', dxu(H_B3)%d(31,32,14)

 if (ltrace) write(0,*)myid,'] >>> RK3: Done step 1'
 !if (ltrace) call output3DA(u2(H_B1)%d,"h_b1posONE",xx,yy,zz,time,myid)
 !if (ltrace) write(0,*)myid, '] >>>  I1 u2(H_B1)%d(31,32,14)=',u2(H_B1)%d(31,32,14)
 
 
  !######### 2nd RK step #############  
  else if (iter .eq. 2) then
  if (ltrace) then
     write(0,*)myid,'] >>> RK3: begin second step'
  end if

  if (boundary_conditions.eq.1.or.boundary_conditions.eq.2.or.boundary_conditions.eq.7 .and. .not.no_bc) then
     if (ltraceB) write(0,*)myid,'] >>> RK3: hyperboundary()'
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

  call hyperrhs(   urk1,  u2,   v,    w,   dxu,  dyu,  dzu,  dxv, dyv, dzv, &
                   imask2,     par)

  if (ltrace2) then
    do m=1,NU
       write(0,*) myid, ']  u2   - E ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) myid, ']  urk2 - E ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) myid,'] '
  end if
       
  if (dissipation > 0) then
    if (ltrace) write(0,*)myid,'] >>> RK3: hyperdissipation()', sigma_diss
    do m=1,NU
      if (dissipation .ne. 3) then
      call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,   u2(m)%d,    &
                            w,         imask,     par)  

      do k = 1,nz
      do j = 1,ny
      do i = 1,nx
         imask_pt = UNDEF
         imask_pt = DEF
         if (use_mask .eq. 1) imask_pt = imask(i,j,k)
         if (imask_pt .ne.UNDEF) then
         !write(*,*) '] RK3: applying dissipation B:',imask(i,j,k)
         urk1(m)%d(i,j,k) = urk1(m)%d(i,j,k) + &
     &           W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
         end if
      end do
      end do
      end do
      else
        call apply_diss_torhs( urk1(m)%d, u0(m)%d, sigma_diss,nx,ny,nz)
      end if
    end do
  end if

  if (boundary_conditions.gt.2.and.bc_type.ne.3.and.boundary_conditions.ne.8.and. .not.no_bc) then
     if (ltraceB) write(0,*)myid,'] >>> RK3: hyperboundary()'
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

  if (ltrace2) then
    do m=1,NU
       write(0,*) myid, ']  u2   - F ',m,myl2norm3d(u2(m)%d, nx, ny, nz)
       write(0,*) myid, ']  urk2 - F ',m,myl2norm3d(urk1(m)%d, nx, ny, nz)
    end do
    write(0,*) myid,'] '
  end if
       
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
         if (update_scheme .eq. 3) then
            u2(m)%d(i,j,k) = 0.75* u0(m)%d(i,j,k) + 0.25 * u2(m)%d(i,j,k) &
                                                  + 0.25 * dt * urk1(m)%d(i,j,k)
         else
            write(0,*) "Error, update_scheme not known"            
         end if
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

  if (boundary_conditions .eq. 8) then
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

  if (bc_type.eq.3.and. .not.no_bc) then
     if (ltraceB) write(*,*) myid,'] RK3: Using copy boundary condition:',iter
      do m=1, NU
         call copy_bc(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      end do
      !
      ! This will not compile in other projects, because it refers specifically
      ! to a field which does not exist in them:
      !
      !call jet_bc(u2(H_RHO)%d,    chr,xx,yy,zz,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      !call jet_bc(u2(H_RHO_U_X)%d,chr,xx,yy,zz,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
  end if

  if (assume_symmetry .ne. 0 .and. .false.) then
     !
     ! Enforce symmetry conditions
     !
     if (.true.) write(*,*) myid,'] RK3: Enforcing symmetry: ', iter
      do m=1, NU
         call enforce_symmetry(u2(m)%d,chr,assume_symmetry,nx,ny,nz)
      end do
  end if

  if (use_mask .eq. 1) then
    call mask_gfuncs(u2, w, par)
  end if
  
 !if (ltrace) call output3DA(u2(H_B1)%d,"h_b1posTWO",xx,yy,zz,time,myid)
 !if (ltrace) write(0,*)myid, '] >>>  I2 u2(H_B1)%d(31,32,14)=',u2(H_B1)%d(31,32,14)
  
  !######### 3rd RK step ##############
  else if (iter .eq. 3) then
  if (ltrace) then
     write(0,*)myid,'] >>> RK3: begin third step'
  end if

  if (boundary_conditions.eq.1.or.boundary_conditions.eq.2.and. .not.no_bc) then
     if (ltraceB) write(0,*)myid,'] >>> RK3: hyperboundary()'
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

  call hyperrhs(   urk1,  u2,   v,     w,    dxu,  dyu,  dzu,  dxv, dyv, dzv, &
                   imask2,      par) 

  if (dissipation > 0) then
    if (ltrace) write(0,*)'>>> RK3: hyperdissipation()', sigma_diss
    do m=1,NU
      if (dissipation .ne. 3) then
      call hyperdissipation(dxu_diss,  dyu_diss,  dzu_diss,   u2(m)%d,    &
                            w,         imask,     par)

      do k = 1,nz
      do j = 1,ny
      do i = 1,nx
         imask_pt = UNDEF
         imask_pt = DEF
         if (use_mask .eq. 1) imask_pt = imask(i,j,k)
         if (imask_pt .ne.UNDEF) then
         !write(*,*) '] RK3: applying dissipation C:',imask(i,j,k)
         urk1(m)%d(i,j,k) = urk1(m)%d(i,j,k)+&
     &          W(H_WDISS)%d(i,j,k)*sigma_diss*(dxu_diss(i,j,k)+dyu_diss(i,j,k)+dzu_diss(i,j,k))
         end if
      end do
      end do
      end do

      else
        call apply_diss_torhs( urk1(m)%d, u0(m)%d, sigma_diss,nx,ny,nz)
      end if
    end do
  end if

  if (boundary_conditions.gt.2.and.bc_type.ne.3.and.boundary_conditions.ne.8.and. .not.no_bc) then
     if (ltraceB) write(0,*)myid,'] >>> RK3: hyperboundary()'
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

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
         if (update_scheme .eq. 3) then
            u2(m)%d(i,j,k) = one_third * u0(m)%d(i,j,k) &
					+ 2.*one_third * u2(m)%d(i,j,k) &
                                        + 2.*one_third * dt * urk1(m)%d(i,j,k)
         else
            write(0,*) "Error, update_scheme not known"            
         end if
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

  if (bc_type.eq.3.and. .not.no_bc) then
     if (ltraceB) write(*,*) myid,'] RK3: Using copy boundary condition:',iter
      do m=1, NU
         call copy_bc(u2(m)%d,chr,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      end do
      !
      ! This will not compile in other projects, because it refers specifically
      ! to a field which does not exist in them:
      !
      !call jet_bc(u2(H_RHO)%d,    chr,xx,yy,zz,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
      !call jet_bc(u2(H_RHO_U_X)%d,chr,xx,yy,zz,nx,ny,nz,2,nx-1,2,ny-1,2,nz-1)
  end if

  if (assume_symmetry .ne. 0 .and. .false.) then
     !
     ! Enforce symmetry conditions
     !
     if (.true.) write(*,*) myid,'] RK3: Enforcing symmetry: ', iter
      do m=1, NU
         call enforce_symmetry(u2(m)%d,chr,assume_symmetry,nx,ny,nz)
      end do
  end if


  if (use_mask .eq. 1) then
    call mask_gfuncs(u2, w, par)
  end if

  if (boundary_conditions .eq. 8) then
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

  if (boundary_conditions.eq.1.or.boundary_conditions.eq.2.and.bc_type.ne.3.and. .not.no_bc) then
     if (ltraceB) write(0,*)myid,'] >>> RK3: hyperboundary()'
     call hyperboundary(u2, dxu, dyu, dzu, urk1, v, w, imask, par)
  end if

! call HyperAnalysis(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, par)
     
  end if !---iter

  !###################################
  !######### End of RK   #############
  !###################################
    
  
  if (allocated(dxu_diss) .and. allocated(dyu_diss) .and. allocated(dzu_diss))&
  then
    if (ltrace2) write(0,*)'>>> RK3: about to deallocate diss vars...'
    deallocate(dzu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'>>> RK3: Could not deallocate dzu_diss: istat=', istat
    end if
    deallocate(dxu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'>>> RK3: Could not deallocate dxu_diss: istat=', istat
    end if
    deallocate(dyu_diss, STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'>>> RK3: Could not deallocate dyu_diss: istat=', istat
    end if
    if (ltrace2) write(0,*)'>>> RK3: ...diss vars deallocated'
  else
    write(0,*)'>>> RK3: diss vars not allocated as expected!'
  end if


  if (use_mask .eq. 1) then
    call unfixup_mask(w(H_MASK)%d,nx,ny,nz)
    if (allocated(imask)) then
      deallocate(imask,imask2)
    else
      write(0,*)'>>> RK3: imask not allocated as expected!'
    end if
  end if

  if (ltrace) then
     if (ltrace2) then
        do m=1,NU
           write(0,*)myid, ']  RK3: '
           write(0,*)myid, ']  RK3: u2   - ',iter,m,myl2norm3d(u2(m)%d, nx, ny, nz)
           write(0,*)myid, ']  RK3: u0   - ',iter,m,myl2norm3d(u0(m)%d, nx, ny, nz)
           write(0,*)myid, ']  RK3: '
        end do
     end if
     !call output3DA(u2(H_B1)%d,"h_b1posFIN",xx,yy,zz,time,myid)
    !if (ltrace) write(0,*)myid, '] >>> 9 u2(H_B1)%d(31,32,14)=',u2(H_B1)%d(31,32,14)
     write(0,*)myid, '] >>>RK3: return'
  end if




  return
end subroutine RK3


!
!
!
subroutine copy_bc(q, chr, numx, numy, numz, xlwb, xupb, &
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
       q(1,J,K) = q(numx-1,J,K)      
!      q(1,   J,K) = q(2,     J,K)
      if (nint(chr(numx,J,K)).ne.1) &
       q(numx,J,K) = q(2,J,K)      
!      q(numx,J,K) = q(numx-1,J,K)
   enddo
enddo

do K = 1, numz
   do I = 1, numx
      if (nint(chr(I,1,   K)).ne.1) &
       q(I,1,K) = q(I,numy-1,K)
!      q(I,1,   K) = q(I,2,     K)      
      if (nint(chr(I,numy,K)).ne.1) &
      !if (nint(chr(I,numy,K)).ne.1) then
       q(I,numy,K) = q(I,2,K)
!      q(I,numy,K) = q(I,numy-1,K)      
      !else
      !write(*,*) i,k,chr(I,numy,k)
      !end if
   enddo
enddo

do J = 1, numy
   do I = 1, numx
      if (nint(chr(I,J,1   )).ne.1) &
       q(I,J,1) = q(I,J,numz-1 )
!      q(I,J,1   ) = q(I,J,2     )      
      if (nint(chr(I,J,numz)).ne.1) &
       q(I,J,numz) = q(I,J,2)      
!      q(I,J,numz) = q(I,J,numz-1)
   enddo
enddo

end  subroutine copy_bc

!
!
!
subroutine jet_bc(q, chr, x,y,z,numx, numy, numz, xlwb, xupb, &
                              ylwb, yupb, zlwb, zupb )
implicit  none

! subroutine arguments
CCTK_INT, intent (in) :: numx
CCTK_INT, intent (in) :: numy
CCTK_INT, intent (in) :: numz

CCTK_REAL, dimension(numx,numy,numz), intent (inout) :: q
CCTK_REAL, dimension(numx,numy,numz), intent (in)    :: chr
CCTK_REAL, dimension(numx,numy,numz), intent (in)    :: x,y,z

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
CCTK_REAL            :: amp    = 1.d0
CCTK_REAL            :: R0     = 0.d0
CCTK_REAL            :: delta  = 2.d0
CCTK_REAL            :: yc     = 0.d0
CCTK_REAL            :: zc     = 0.d0

return

if (ltrace) then
   write(*,*) 'jet_bc: numx/y/z = ',numx,numy,numz
end if

i = 1
do K = 1, numz
   do J = 1, numy
      if (nint(chr(1,   J,K)).ne.1) &
      q(1,   J,K) = amp * exp( -( (y(i,j,k)-yc)**2+(z(i,j,k)-zc)**2) )
   enddo
enddo

end  subroutine jet_bc

subroutine enforce_symmetry(q, chr, assume_symmetry, numx, numy, numz)
implicit  none

! subroutine arguments
 CCTK_INT, intent (in) :: numx
 CCTK_INT, intent (in) :: numy
 CCTK_INT, intent (in) :: numz

 CCTK_REAL, dimension(numx,numy,numz), intent (inout) :: q
 CCTK_REAL, dimension(numx,numy,numz), intent (in)    :: chr

 CCTK_INT, intent (in) :: assume_symmetry

!
!  impose boundary conditions
!
!
 CCTK_INT :: I, J, K

logical              :: ltrace = .false.


if (ltrace) then
   write(*,*) 'enforce_symmetry: numx/y/z = ',numx,numy,numz
end if

if (assume_symmetry .eq. 3) then
   K = 1
   do J = 1, numy
      do I = 1, numx
         ! see ../../include/chr.inc for value CHR_Refl_bdy
         if (nint(chr(I,J,K)).eq.29) then
            q(I,J,1) = q(I,J,2)
         end if
      enddo
   enddo
else
   if (ltrace) then
      write(*,*) 'enforce_symmetry: Unknown/Unimplemented value for assume_symmetry: ',assume_symmetry
   end if
end if


end  subroutine enforce_symmetry

	subroutine output3DA(var,name,x1,x2,x3,time,myproc)
	implicit none
	CCTK_REAL, dimension(:,:,:) :: var
	CCTK_REAL, dimension(:,:,:) :: x1,x2,x3
	CCTK_REAL :: time
	character(len=*) :: name
	
	CCTK_INT myproc, istat
		
	real*8, allocatable, dimension(:) :: tempcoord
	integer :: nx, ny, nz, ret, gft_out_full
	character(20):: g11f
	character(32):: cnames
	
	nx=size(var(:,1,1))
	ny=size(var(1,:,1))
	nz=size(var(1,1,:))
	
	cnames = 'x|y|z'

	allocate(tempcoord(nx+ny+nz),STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'output3DA >>> can not allocate tempcoord'
          write(0,*)'             size: ',nx,ny,nz
        end if

	
	tempcoord(1:nx) = x1(:,1,1)
	tempcoord(nx+1:nx+ny) = x2(1,:,1)
	tempcoord(nx+ny+1:nx+ny+nz) = x3(1,1,:)
	
        write (g11f,'(A,1i3)') name,myproc
	if(myproc.lt.0) g11f = name
	
	
          write(0,*)'output3DA >>> outputting ',g11f
	ret=gft_out_full(g11f,time,(/nx,ny,nz/),cnames,3, &
                                tempcoord,var)
				
				
        if (allocated(tempcoord)) then
	  deallocate(tempcoord)
        else
          write(0,*)'output3DA >>> tempcoord not allocated as expected!'
        end if
				
	end subroutine output3DA
	
