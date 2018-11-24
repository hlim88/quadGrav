#include "cctk.h"

!--------------- R E A D   M E ------------------------------------------
!
! Note as of Fri Mar 14 10:17:45 CST 2003
!
! Currently this routine is very simple minded.  It does the following:
!   (a) Assigns pointers to u, v, w, as well as the derivatives,
!   (b) Defines the parameter array
!   (c) Calls auxvars to define v
!   (d) Calls hyperderivs to calculate dxu, dyu, dzu, dxv, dyv, dzv
!   (e) Calls a general analysis routine that should be in the equation thorn,
!       which is not a point-wise routine.       
!
! The analysis functions are presumably stored in w, but no assumptions
! at this level are made.  NOTE: Point-wise derivatives are NOT supported.
! In the future a point-wise analysis routine probably should be written.
!
!------------------------------------------------------------------------


  subroutine HyperAnalysis(u0, u2, v, w, dxu,dyu,dzu, dxv, dyv, dzv, par)
  use params
  use GF
  use UTILEQS
  use MOD_ANALYSIS
  use mod_mask
  use hypercoords
  use HYPER_AUXVARS
  use HYPER_DERIVS


  implicit none
  
  type(gridfunction), dimension(NU)    :: u2, u0, dxu, dyu, dzu
  type(gridfunction), dimension(NV)    :: v, dxv, dyv, dzv
  type(gridfunction), dimension(NW)    :: w
  CCTK_REAL, dimension(NPAR)           :: par
  
  CCTK_REAL, dimension(NU)             :: u_pt,u0_pt,urk_pt
  CCTK_REAL, dimension(NU)             :: dxu_pt,dyu_pt,dzu_pt  
  CCTK_REAL, dimension(NV)             :: v_pt,dxv_pt,dyv_pt,dzv_pt
  CCTK_REAL, dimension(NW)             :: w_pt

  CCTK_INT  :: ii, jj ,kk, m, nx, ny, nz, nt, istat, err, use_mask, imask_pt, order   
  CCTK_INT  :: ilo,ihi,jlo,jhi,klo,khi,rank,ierr, ng
   
  CCTK_INT, dimension(:,:,:), allocatable :: imask,imask2
  logical, parameter                   :: ltrace = .false.
  logical, parameter                   :: ltrace2 = .false.  
  CCTK_INT:: proc_return_myid, myid
  CCTK_REAL :: myl2norm3d

  myid = proc_return_myid()
     !write(0,*)myid, '] >>>  HyperAnalysis Skipping, returning'
     !return

  if (ltrace.or.ltrace2) then
     write(0,*)myid, '] >>>  HyperAnalysis begins'
  end if
  
  if (ltrace) then
    write(0,*)myid,' >>> HyperAnalysis'
    write(0,*)myid,' >>> p_point_wise_derivatives = ',P_POINT_WISE_DERIVATIVES
    write(0,*)myid,' >>> p_point_wise_analysis = ',nint(par(P_POINT_WISE_ANALYSIS))
  end if
 
  nx = nint(par(P_NX))
  ny = nint(par(P_NY))
  nz = nint(par(P_NZ))
  if (ltrace.or.ltrace2) write(0,*)myid, '] >>>  HyperAnalysis: nx/y/z: ',nx,ny,nz
  ng = par(P_BOUND_WIDTH)  
  if (ltrace.or.ltrace2) write(0,*)myid, '] >>>  HyperAnalysis: ng: ',ng
  use_mask = nint(par(P_USE_MASK))
  if (ltrace.or.ltrace2) write(0,*)myid, '] >>>  HyperAnalysis: use_mask: ',use_mask
  order = nint(par(P_DERIV_ORDER))
  if (ltrace.or.ltrace2) write(0,*)myid, '] >>>  HyperAnalysis: order: ',order
  
  !define coordinates if not using the default cartesian computational coords.
  if (P_ENABLE_ALT_COORDS==1) then
    if (ltrace.or.ltrace2) write(0,*)myid, '] >>>  HyperAnalysis calling def_alt_coords'
    call def_alt_coords(w, par)
  end if

  if (use_mask .eq. 1) then
    if (ltrace.or.ltrace2) write(0,*)myid, '] >>>  HyperAnalysis allocating imask'
    allocate(imask(nx,ny,nz),STAT=istat)
    if (istat .ne. 0) then
      write(0,*)myid,'Hyperanalysis: Could not allocate memory for imask array'
      write(0,*)myid,'     size : ',nx,ny,nz
    end if

    if (ltrace.or.ltrace2) write(0,*)myid, '] >>>  HyperAnalysis allocating imask2'
    allocate(imask2(nx,ny,nz),STAT=istat)
    if (istat .ne. 0) then
      write(0,*)myid,'HyperAnalysis: Could not allocate memory for imask2 array'
      write(0,*)myid,'     size : ',nx,ny,nz
    end if

    if (ltrace.or.ltrace2) write(0,*)myid, '] >>>  HyperAnalysis set_mask'
    call set_mask(w, imask,imask2, par)
  end if
  
  if (ltrace.or.ltrace2) write(0,*)myid, '] >>>  HyperAnalysis calling hyperauxvars'
  call hyperauxvars(v, w, u2, par)  
  
  if (ltrace2) then
    do m=1,NU
       write(0,*) myid, '] hyperanalysis: u0 - 1A ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
    end do
    do m=1,NV
       write(0,*) myid, '] hyperanalysis: v - 1A ',m,myl2norm3d(v(m)%d, nx, ny, nz)
    end do
    do m=1,NW
       write(0,*) myid, '] hyperanalysis: w - 1A ',m,myl2norm3d(w(m)%d, nx, ny, nz)
    end do
  end if
  
  
  if (nint(par(P_POINT_WISE_ANALYSIS)) .eq. 0) then 
    !
    if (P_POINT_WISE_DERIVATIVES .eq. 0) then
      if (ltrace) then
        write(0,*)myid,'>>> HYPER_ANALYSIS: calling hyperderivs NU/V:',NU,NV
      end if
      call hyperderivs(dxu,   dyu,  dzu,  dxv,  dyv,  dzv, &
                       u2,    v,    w,     imask2,      par)
    end if
    !
    if (ltrace) write(0,*)'>>> HYPER_ANALYSIS: calling analysis'
    !
    call analysis(u2, u0, dxu, dyu, dzu, v, dxv, dyv, dzv, w, par)
    !
  else
    !
    if (ltrace) write(0,*)'>>> HYPER_ANALYSIS: calling hyperderivs and analysis point to point'
    !
     do kk = 1, nz
       do jj = 1, ny
         do ii = 1, nx

          if ((w(H_MASK)%d(ii,jj,kk).gt.0) .OR. (use_mask .eq. 0) ) then
  
            do m = 1, NU
              u_pt(m) = u2(m)%d(ii,jj,kk)
            end do
            do m = 1, NV
              v_pt(m) = v(m)%d(ii,jj,kk)
            end do
            do m = 1, NW
              w_pt(m) = w(m)%d(ii,jj,kk)
            end do
            if (use_mask .eq. 1) then
               ! This applies if imask() was allocated
               imask_pt = imask2(ii,jj,kk)
            else
               ! Just a dummy value:
               imask_pt = 1
            end if

            call hyperderivs_pt(dxu_pt,   dyu_pt,  dzu_pt,                &
                                dxv_pt,   dyv_pt,  dzv_pt,                &
                                u2,       v,       w_pt,    imask_pt,     & 
                                ii,       jj,      kk,      par           )

            call analysis_pt(u_pt, u_pt, urk_pt, dxu_pt, dyu_pt, dzu_pt, &
                             v_pt, dxv_pt, dyv_pt, dzv_pt, w_pt, par)      

            do m=1, NV
              v(m)%d(ii,jj,kk) = v_pt(m)
            end do
            do m = 1, NW
              w(m)%d(ii,jj,kk) = w_pt(m)
            end do
 
         end if
  
    
        end do
      end do
    end do
    
  end if      
  
  if (ltrace) then
     write(0,*)'>>> HYPER_ANALYSIS: return from analysis'
  end if
    
  if (ltrace2) then
    do m=1,NU
       write(0,*) myid, '] hyperanalysis: u0 - 1B ',m,myl2norm3d(u0(m)%d, nx, ny, nz)
    end do
    do m=1,NV
       write(0,*) myid, '] hyperanalysis: v - 1B ',m,myl2norm3d(v(m)%d, nx, ny, nz)
    end do
    do m=1,NW
       write(0,*) myid, '] hyperanalysis: w - 1B ',m,myl2norm3d(w(m)%d, nx, ny, nz)
    end do
  end if
  
  
  if (use_mask .eq. 1) then
    if (allocated(imask)) then
      if (ltrace) print *,' about to deallocate imask'
      deallocate(imask)
    else
      write(0,*)'>>> HYPER_ANALYSIS: imask not allocated as expected!'
    end if
    if (allocated(imask2)) then
      if (ltrace) print *,' about to deallocate imask2'
      deallocate(imask2)
    else
      write(0,*)'>>> HYPER_ANALYSIS: imask2 not allocated as expected!'
    end if
    if (ltrace) print *,' end of deallocate imask'
  end if
  
  if (ltrace) then
     write(0,*)'>>> HYPER_ANALYSIS: exit from HyperAnalysis'
  end if

  return
  end subroutine HyperAnalysis
    
  
  
