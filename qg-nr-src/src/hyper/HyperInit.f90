#include "cctk.h"

  subroutine HyperInit(u0, u2, v, w, dxu, dyu, dzu, dxv, dyv, dzv, par)
  
  use params
  use GF
  use UTILEQS
  use initial
  use mod_mask
  use hypercoords


  implicit none
  
  CCTK_INT    nx, ny, nz, nt, ii, istat, use_mask, i, initial_analysis
  type(gridfunction), dimension(NU)    :: u0,u2,dxu,dyu,dzu
  type(gridfunction), dimension(NV)    :: v,dxv, dyv, dzv
  type(gridfunction), dimension(NW)    :: w

  CCTK_REAL, dimension(NPAR)           :: par
  CCTK_INT, dimension(:,:,:), allocatable :: imask,imask2
  logical, parameter                   :: ltrace = .false.
  CCTK_INT :: rcc
  integer  :: myid, proc_return_myid

  ! external function
  CCTK_REAL myl2norm3d
    
  myid = proc_return_myid()
  nx = nint(par(P_NX))
  ny = nint(par(P_NY))
  nz = nint(par(P_NZ))
  use_mask = nint(par(P_USE_MASK))
  initial_analysis = nint(par(P_INITIAL_ANALYSIS))
    
!assign pointers
  if (ltrace) then
    print *,myid,'------------pointers assigned-----------------------'
    print *,myid,'>>> HYPER_INIT: nx,ny,nz ', nx,ny,nz
    print *,myid,'>>> HYPER_INIT: shape of w(mask) ',shape(w(H_MASK)%d)
    print *,myid,'>>> HYPER_INIT: shape of u(1)    ',shape(u0(1)%d)
  end if
  if (ltrace) then
    write(0,*)myid,'>>> HYPER_INIT: x ',w(H_X)%d(1,1,1),w(H_X)%d(2,1,1),w(H_X)%d(3,1,1)
    write(0,*)myid,'>>> HYPER_INIT: y ',w(H_Y)%d(1,1,1),w(H_Y)%d(1,2,1),w(H_Y)%d(1,3,1)
    write(0,*)myid,'>>> HYPER_INIT: z ',w(H_Z)%d(1,1,1),w(H_Z)%d(1,1,2),w(H_Z)%d(1,1,3)
  end if

    
  !define coordinates if not using the default cartesian computational coords.
  if (P_ENABLE_ALT_COORDS==1) then
    call def_alt_coords(w, par)
  end if

  if (use_mask .eq. 1) then
    if (ltrace) print *,' allocating imask'
    allocate(imask2(nx,ny,nz),imask(nx, ny, nz),STAT=istat)
    if (istat .ne. 0) then
      write(0,*)'Could not allocate imask in HyperInitial'
      write(0,*)'   nx = ',nx
      write(0,*)'   ny = ',ny
      write(0,*)'   nz = ',nz
      stop
    end if

    if (ltrace) print *,' stat = ',istat
!    write(0,*) 'about to ini'
    call set_mask(w, imask,imask2, par)
!    rcc = check_imask(imask,par)
!    if(rcc.eq.0) write(0,*) 'ops, inini'
  else 
    allocate(imask2(nx,ny,nz),imask(nx, ny, nz),STAT=istat)
    imask = 1.
  end if

  if (ltrace) print*, myid,'>>> about to enter initial '

  call initialdata(u0,u2,v,w,imask,par)
 
  if (ltrace) print*, myid,'>>> BACK from initial '

  if (use_mask .eq. 1) then
    call mask_gfuncs(u0, w, par)
    if (ltrace) print *,myid,'>>> about to deallocate imask'

    if (allocated(imask)) then
      deallocate(imask,imask2)
    else
      write(0,*)myid,'>>> HYPER_INIT: imask not allocated as expected!'
    end if

    if (ltrace) print *,' end of deallocate imask'
  end if

  if (ltrace) then
    write(0,*) myid,'>>> HyperInit:  Initial data norms', NU
    do i = 1, NU
      write(0,*) '  ',i,myl2norm3d(u0(i)%d, nx, ny, nz)
    end do
  end if

  if (initial_analysis .eq. 1) then
    if (ltrace) write(0,*) myid,'>>> HyperInit: calling to HyperAnalysis at initial time'  
    call HyperAnalysis(u0, u0, v, w, dxu, dyu, dzu, dxv, dyv, dzv, par)
  end if     
     
 if (ltrace) write(0,*) myid,'>>> HyperInit: DONE'

  return
  end subroutine HyperInit
    
  
  
