!-----------------------------------------------------------------------
!
!    $Id$
!
!-----------------------------------------------------------------------

#include "cctk.h"

MODULE HYPER_DERIVS_42
use params

CONTAINS

!-----------------------------------------------------------------------
!
!    SUBROUTINE D41NoMask:
!
!    implement Strand derivative.  4th order interior, 2nd order boundary
!
!-----------------------------------------------------------------------
  subroutine D42NoMask(Du, u, direction, par)
    implicit none
    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT :: direction
    logical, parameter :: ltrace2 = .false.
    integer :: myid, proc_return_myid


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, j, k,nx, ny, nz, shp(3)
    

    myid = proc_return_myid()
    if(ltrace2)write(0,*)myid,'] D42NoMask: Enter, direction:',direction
    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)
    
    if (direction.eq.1) then 
      if(ltrace2)write(0,*)myid,'] D42NoMask: direction==1'
      Du(1,:,:) = (-48.d0*u(1,:,:) + 59.d0*u(2,:,:) - 8.d0*u(3,:,:) - 3.0d0*u(4,:,:))/(34.d0*dx)
      Du(2,:,:) = (u(3,:,:) - u(1,:,:))/(2.d0*dx)
      Du(3,:,:) = (8.d0*u(1,:,:) - 59.d0*u(2,:,:) + 59.d0*u(4,:,:) - 8.d0*u(5,:,:))/(86.d0*dx)
      Du(4,:,:) = (3.d0*u(1,:,:) - 59.d0*u(3,:,:) + 64.d0*u(5,:,:) - 8.d0*u(6,:,:))/(98.d0*dx)
      do i = 5, nx-4
        Du(i,:,:) = (-u(i+2,:,:) + 8.d0*u(i+1,:,:) - 8.d0*u(i-1,:,:) + u(i-2,:,:))/(12.d0*dx)
      end do
      Du(nx-3,:,:) = (-3.d0*u(nx,:,:) + 59.d0*u(nx-2,:,:) - 64.d0*u(nx-4,:,:) + 8.d0*u(nx-5,:,:))/(98.d0*dx)
      Du(nx-2,:,:) = (-8.d0*u(nx,:,:) + 59.d0*u(nx-1,:,:) - 59.d0*u(nx-3,:,:) + 8.d0*u(nx-4,:,:))/(86.d0*dx)
      Du(nx-1,:,:) = (u(nx,:,:) - u(nx-2,:,:))/(2.d0*dx)
      Du(nx,:,:)   = (48.d0*u(nx,:,:) - 59.d0*u(nx-1,:,:) + 8.d0*u(nx-2,:,:) + 3.0d0*u(nx-3,:,:))/(34.d0*dx)
    else if (direction.eq.2) then 
      if(ltrace2)write(0,*)myid,'] D42NoMask: direction==2'
      Du(:,1,:) = (-48.d0*u(:,1,:) + 59.d0*u(:,2,:) - 8.d0*u(:,3,:) - 3.0d0*u(:,4,:))/(34.d0*dy)
      Du(:,2,:) = (u(:,3,:) - u(:,1,:))/(2.d0*dy)
      Du(:,3,:) = (8.d0*u(:,1,:) - 59.d0*u(:,2,:) + 59.d0*u(:,4,:) - 8.d0*u(:,5,:))/(86.d0*dy)
      Du(:,4,:) = (3.d0*u(:,1,:) - 59.d0*u(:,3,:) + 64.d0*u(:,5,:) - 8.d0*u(:,6,:))/(98.d0*dy)
      do j = 5, ny-4
        Du(:,j,:) = (-u(:,j+2,:) + 8.d0*u(:,j+1,:) - 8.d0*u(:,j-1,:) + u(:,j-2,:))/(12.d0*dy)
      end do
      Du(:,ny-3,:) = (-3.d0*u(:,ny,:) + 59.d0*u(:,ny-2,:) - 64.d0*u(:,ny-4,:) + 8.d0*u(:,ny-5,:))/(98.d0*dy)
      Du(:,ny-2,:) = (-8.d0*u(:,ny,:) + 59.d0*u(:,ny-1,:) - 59.d0*u(:,ny-3,:) + 8.d0*u(:,ny-4,:))/(86.d0*dy)
      Du(:,ny-1,:) = (u(:,ny,:) - u(:,ny-2,:))/(2.d0*dy)
      Du(:,ny,:)   = (48.d0*u(:,ny,:) - 59.d0*u(:,ny-1,:) + 8.d0*u(:,ny-2,:) + 3.0d0*u(:,ny-3,:))/(34.d0*dy)
    else if (direction==3) then 
      if(ltrace2)write(0,*)myid,'] D42NoMask: direction==3'
      Du(:,:,1) = (-48.d0*u(:,:,1) + 59.d0*u(:,:,2) - 8.d0*u(:,:,3) - 3.0d0*u(:,:,4))/(34.d0*dz)
      Du(:,:,2) = (u(:,:,3) - u(:,:,1))/(2.d0*dz)
      Du(:,:,3) = (8.d0*u(:,:,1) - 59.d0*u(:,:,2) + 59.d0*u(:,:,4) - 8.d0*u(:,:,5))/(86.d0*dz)
      Du(:,:,4) = (3.d0*u(:,:,1) - 59.d0*u(:,:,3) + 64.d0*u(:,:,5) - 8.d0*u(:,:,6))/(98.d0*dz)
      do k = 5, nz-4
        Du(:,:,k) = (-u(:,:,k+2) + 8.d0*u(:,:,k+1) - 8.d0*u(:,:,k-1) + u(:,:,k-2))/(12.d0*dz)
      end do
      Du(:,:,nz-3) = (-3.d0*u(:,:,nz) + 59.d0*u(:,:,nz-2) - 64.d0*u(:,:,nz-4) + 8.d0*u(:,:,nz-5))/(98.d0*dz)
      Du(:,:,nz-2) = (-8.d0*u(:,:,nz) + 59.d0*u(:,:,nz-1) - 59.d0*u(:,:,nz-3) + 8.d0*u(:,:,nz-4))/(86.d0*dz)
      Du(:,:,nz-1) = (u(:,:,nz) - u(:,:,nz-2))/(2.d0*dz)
      Du(:,:,nz)   = (48.d0*u(:,:,nz) - 59.d0*u(:,:,nz-1) + 8.d0*u(:,:,nz-2) + 3.0d0*u(:,:,nz-3))/(34.d0*dz)
    else 
      if(ltrace2)write(0,*)myid,'] D42NoMask: Undefined direction: ',direction
      write(*,*) "D42NoMask: Wrong direction in the derivative" 
      !STOP
    end if 
    if(ltrace2)write(0,*)myid,'] D42NoMask: Done. ',direction

    return
  end subroutine D42NoMask


!-----------------------------------------------------------------------
!
!    D42
!
!  At times a derivative operator may be legitimately undefined.
!  This can happen in the ghost zones, if they are very close to
!  a mask boundary.  At such times the derivative is simply set
!  to zero.  For debugging purposes, the ltrace2 parameter can be
!  used to control whether an message is printed for undefined 
!  derivatives.
!
!-----------------------------------------------------------------------
  subroutine D42(Du, u, direction, x, y, z, mask, imask, par)
    implicit none

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: x, y, z
    CCTK_INT :: direction
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: mask
    CCTK_INT, dimension(:,:,:), INTENT(IN)   :: imask


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: i, j, k,nx, ny, nz, shp(3)
    CCTK_INT                                 :: d
    CCTK_INT  :: gz_xmin, gz_xmax, gz_ymin, gz_ymax, gz_zmin, gz_zmax
    CCTK_INT  :: bbox1, bbox2, bbox3, bbox4, bbox5, bbox6
    logical, parameter :: debug_deriv = .false.

    !  If ltrace2 is TRUE, then a warning message is printed for 
    !  undefined derivatives are flagged.
    !  If ltrace2 is FALSE, then no warning message is printed.
    logical, parameter :: ltrace2 = .false.
    integer :: myid, proc_return_myid
    real(kind=8) :: myl2norm3d

    myid = proc_return_myid()
    if (ltrace2) write(0,*) myid,'] D42: Enter, direction:',direction

    shp = shape(u)
    nx  = shp(1)
    ny  = shp(2)
    nz  = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    if (ltrace2) then
       write(0,*)myid,'] d42:I:u,du:',direction,myl2norm3d(Du, nx, ny, nz), myl2norm3d(u, nx,ny,nz)
       write(0,*)myid,'] d42:nghostzones_x:',nint(par(P_NGHOSTZONES_X))
       write(0,*)myid,'] d42:nghostzones_y:',nint(par(P_NGHOSTZONES_Y))
       write(0,*)myid,'] d42:nghostzones_z:',nint(par(P_NGHOSTZONES_Z))
    end if
    !
    ! Check for ghostzones
    !
    gz_xmin = 1
    gz_xmax = nx 
    gz_ymin = 1
    gz_ymax = ny 
    gz_zmin = 1
    gz_zmax = nz 

    bbox1 = nint(par(P_BBOX1))
    bbox2 = nint(par(P_BBOX2))
    bbox3 = nint(par(P_BBOX3))
    bbox4 = nint(par(P_BBOX4))
    bbox5 = nint(par(P_BBOX5))
    bbox6 = nint(par(P_BBOX6))

    if (bbox1 .eq. 0) then
      gz_xmin = gz_xmin + nint(par(P_NGHOSTZONES_X))
    end if
    if (bbox2 .eq. 0) then
      gz_xmax = gz_xmax - nint(par(P_NGHOSTZONES_X))
    end if

    if (bbox3 .eq. 0) then
      gz_ymin = gz_ymin + nint(par(P_NGHOSTZONES_Y))
    end if
    if (bbox4 .eq. 0) then
      gz_ymax = gz_ymax - nint(par(P_NGHOSTZONES_Y))
    end if

    if (bbox5 .eq. 0) then
      gz_zmin = gz_zmin + nint(par(P_NGHOSTZONES_Z))
    end if
    if (bbox6 .eq. 0) then
      gz_zmax = gz_zmax - nint(par(P_NGHOSTZONES_Z))
    end if

    if (direction.eq.1) then
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = imask(i,j,k) / 10000
        if (d .eq. P_STENCIL_CENTER.and.i.gt.2.and.i.le.nx-2) then
	       Du(i,j,k) = (-u(i+2,j,k) + 8.d0*u(i+1,j,k) - 8.d0*u(i-1,j,k) + u(i-2,j,k))/(12.d0*dx)
        else if (d .eq. P_STENCIL_LEFT           .and.i.le.nx-3) then
          Du(i,j,k) = (-48.d0*u(i,j,k) + 59.d0*u(i+1,j,k) - 8.d0*u(i+2,j,k) - 3.0d0*u(i+3,j,k))/(34.d0*dx)
	     else if (d .eq. P_STENCIL_CENTER_LEFT.and.i.gt.1.and.i.le.nx-1) then
          Du(i,j,k) = (u(i+1,j,k) - u(i-1,j,k))/(2.d0*dx)
	     else if (d .eq. P_STENCIL_CENTER_LEFT1.and.i.gt.2.and.i.le.nx-2) then
          Du(i,j,k) = (8.d0*u(i-2,j,k) - 59.d0*u(i-1,j,k) + 59.d0*u(i+1,j,k) - 8.d0*u(i+2,j,k))/(86.d0*dx) 
	     else if (d .eq. P_STENCIL_CENTER_LEFT2.and.i.gt.3.and.i.le.nx-2) then
          Du(i,j,k) = (3.d0*u(i-3,j,k) - 59.d0*u(i-1,j,k) + 64.d0*u(i+1,j,k) - 8.d0*u(i+2,j,k))/(98.d0*dx)
	     else if (d .eq. P_STENCIL_CENTER_RIGHT2.and.i.gt.2.and.i.le.nx-3) then
          Du(i,j,k) = (-3.d0*u(i+3,j,k) + 59.d0*u(i+1,j,k) - 64.d0*u(i-1,j,k) + 8.d0*u(i-2,j,k))/(98.d0*dx)
	     else if (d .eq. P_STENCIL_CENTER_RIGHT1.and.i.gt.2.and.i.le.nx-2) then
          Du(i,j,k) = (-8.d0*u(i+2,j,k) + 59.d0*u(i+1,j,k) - 59.d0*u(i-1,j,k) + 8.d0*u(i-2,j,k))/(86.d0*dx)
	     else if (d .eq. P_STENCIL_CENTER_RIGHT.and.i.gt.1.and.i.le.nx-1) then
          Du(i,j,k) = (u(i+1,j,k) - u(i-1,j,k))/(2.d0*dx)    
	     else if (d .eq. P_STENCIL_RIGHT.and.i.gt.3) then
	       Du(i,j,k) = (48.d0*u(i,j,k) - 59.d0*u(i-1,j,k) + 8.d0*u(i-2,j,k) + 3.0d0*u(i-3,j,k))/(34.d0*dx)      
	     else if (d .eq. P_STENCIL_MASK) then
          if (debug_deriv) then
            if (nint(mask(i,j,k)) .gt. 0) then
              write(0,*)'Attempt to take undefined derivative'
              write(0,*)'  mask  ',mask(i,j,k), mask(i-1,j,k)
              write(0,*)'  imask ',imask(i,j,k)
              !stop
            end if
          end if
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_UNDEF) then
          if (ltrace2) then
            if ((bbox1==0 .and. i < gz_xmin) .or.  &
                (bbox2==0 .and. i > gz_xmax)) then
              write(0,*)myid,'d42: OK---Undefined x derivative type in ghostzone'
            else
              write(0,*)myid,'d42: ### Undefined x derivative outside of ghost zone!'
              write(0,*)myid,'d42: imask                          ',imask(i,j,k)
              write(0,*)myid,'d42: i, nx                          ',i, nx
              write(0,*)myid,'d42: gz_xmin                        ',gz_xmin
              write(0,*)myid,'d42: gz_xmax                        ',gz_xmax
              write(0,*)myid,'d42: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_X)
              write(0,*)myid,'d42: par(bbox1), par(bbox2)         ',&
                                                par(P_BBOX1), par(P_BBOX2)
              write(0,*)myid,'d42: bbox1, bbox2                   ',&
                                                bbox1,bbox2
              !stop
            end if
          end if
          Du(i,j,k) = 0.0
        else
          write(0,*)myid,'d42: PROBLEM PROBLEM   direction:',direction
          write(0,*)myid,'] d42: Illegal x derivative type ',d
          write(0,*)myid,'] d42: imask                     ',imask(i,j,k)
          write(0,*)myid,']      at location ',x(i,j,k),y(i,j,k),z(i,j,k)
          Du(i,j,k) = 0.0
	       !stop
        end if      
      end do
      end do
      end do
      if(ltrace2)write(0,*)myid,'d42:A:u,du:',direction,myl2norm3d(Du, nx, ny, nz), myl2norm3d(u, nx,ny,nz)
      !
    else if (direction.eq.2) then
      !
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod((imask(i,j,k)/100), 100)
        if (d .eq. P_STENCIL_CENTER.and.j.gt.2.and.j.le.ny-2) then
          Du(i,j,k) = (-u(i,j+2,k) + 8.d0*u(i,j+1,k) - 8.d0*u(i,j-1,k) + u(i,j-2,k))/(12.d0*dy)
        else if (d .eq. P_STENCIL_LEFT.and.           j.le.ny-3) then
          Du(i,j,k) = (-48.d0*u(i,j,k) + 59.d0*u(i,j+1,k) - 8.d0*u(i,j+2,k) - 3.0d0*u(i,j+3,k))/(34.d0*dy)
        else if (d .eq. P_STENCIL_CENTER_LEFT.and.j.gt.1.and.j.le.ny-1) then
          Du(i,j,k) = (u(i,j+1,k) - u(i,j-1,k))/(2.d0*dy)
        else if (d .eq. P_STENCIL_CENTER_LEFT1.and.j.gt.2.and.j.le.ny-2) then
          Du(i,j,k) = (8.d0*u(i,j-2,k) - 59.d0*u(i,j-1,k) + 59.d0*u(i,j+1,k) - 8.d0*u(i,j+2,k))/(86.d0*dy)
        else if (d .eq. P_STENCIL_CENTER_LEFT2.and.j.gt.3.and.j.le.ny-2) then
          Du(i,j,k) = (3.d0*u(i,j-3,k) - 59.d0*u(i,j-1,k) + 64.d0*u(i,j+1,k) - 8.d0*u(i,j+2,k))/(98.d0*dy)
        else if (d .eq. P_STENCIL_CENTER_RIGHT2.and.j.gt.2.and.j.le.ny-3) then
          Du(i,j,k) = (-3.d0*u(i,j+3,k) + 59.d0*u(i,j+1,k) - 64.d0*u(i,j-1,k) + 8.d0*u(i,j-2,k))/(98.d0*dy)  
        else if (d .eq. P_STENCIL_CENTER_RIGHT1.and.j.gt.2.and.j.le.ny-2) then
          Du(i,j,k) = (-8.d0*u(i,j+2,k) + 59.d0*u(i,j+1,k) - 59.d0*u(i,j-1,k) + 8.d0*u(i,j-2,k))/(86.d0*dy)
        else if (d .eq. P_STENCIL_CENTER_RIGHT.and.j.gt.1.and.j.le.ny-1) then
          Du(i,j,k) = (u(i,j+1,k) - u(i,j-1,k))/(2.d0*dy)    
        else if (d .eq. P_STENCIL_RIGHT.and.j.gt.3              ) then
          Du(i,j,k) = (48.d0*u(i,j,k) - 59.d0*u(i,j-1,k) + 8.d0*u(i,j-2,k) + 3.0d0*u(i,j-3,k))/(34.d0*dy)
	else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0
        else if (d .eq. P_STENCIL_UNDEF) then
          if (ltrace2) then
            if ((bbox3==0 .and. j < gz_ymin) .or. &
                (bbox4==0 .and. j > gz_ymax)) then
              write(0,*)myid,'] d2: OK---Undefined y derivative type in ghostzone'
            else
              write(0,*)myid,'d42: PROBLEM PROBLEM   direction:',direction
              write(0,*)myid,'] d42: ### Undefined y derivative outside of ghost zone!'
              write(0,*)myid,'] d42: imask                          ',imask(i,j,k)
              write(0,*)myid,'] d42: j, ny                          ',j, ny
              write(0,*)myid,'] d42: gz_ymin                        ',gz_ymin
              write(0,*)myid,'] d42: gz_ymax                        ',gz_ymax
              write(0,*)myid,'] d42: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_Y)
              write(0,*)myid,'] d42: par(bbox3), par(bbox4)         ',&
                                                par(P_BBOX3), par(P_BBOX4)
              write(0,*)myid,'] d42: bbox3, bbox4                   ',&
                                                bbox3,bbox4
              !stop
            end if
          end if
          Du(i,j,k) = 0.0
        else
          write(0,*)myid,'] d42: Illegal y derivative type   ',d
          write(0,*)myid,'] d42: imask                       ',imask(i,j,k)
          write(0,*)myid,']      at location ',x(i,j,k),y(i,j,k),z(i,j,k)
          Du(i,j,k) = 0.0
	  !stop
        end if
      end do
      end do
      end do

      
    else if (direction.eq.3) then
    
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
        d = mod(imask(i,j,k),100)
        if (d .eq. P_STENCIL_CENTER.and.k.gt.2.and.k.le.nz-2) then
          Du(i,j,k) = (-u(i,j,k+2) + 8.d0*u(i,j,k+1) - 8.d0*u(i,j,k-1) + u(i,j,k-2))/(12.d0*dz)
        else if (d .eq. P_STENCIL_LEFT           .and.k.le.nz-3) then
          Du(i,j,k) = (-48.d0*u(i,j,k) + 59.d0*u(i,j,k+1) - 8.d0*u(i,j,k+2) - 3.0d0*u(i,j,k+3))/(34.d0*dz)
        else if (d .eq. P_STENCIL_CENTER_LEFT.and.k.gt.1.and.k.le.nz-1) then
          Du(i,j,k) = (u(i,j,k+1) - u(i,j,k-1))/(2.d0*dz)
        else if (d .eq. P_STENCIL_CENTER_LEFT1.and.k.gt.2.and.k.le.nz-2) then
          Du(i,j,k) = (8.d0*u(i,j,k-2) - 59.d0*u(i,j,k-1) + 59.d0*u(i,j,k+1) - 8.d0*u(i,j,k+2))/(86.d0*dz)
        else if (d .eq. P_STENCIL_CENTER_LEFT2.and.k.gt.3.and.k.le.nz-2) then
          Du(i,j,k) = (3.d0*u(i,j,k-3) - 59.d0*u(i,j,k-1) + 64.d0*u(i,j,k+1) - 8.d0*u(i,j,k+2))/(98.d0*dz)
        else if (d .eq. P_STENCIL_CENTER_RIGHT2.and.k.gt.2.and.k.le.nz-3) then
          Du(i,j,k) = (-3.d0*u(i,j,k+3) + 59.d0*u(i,j,k+1) - 64.d0*u(i,j,k-1) + 8.d0*u(i,j,k-2))/(98.d0*dz)
        else if (d .eq. P_STENCIL_CENTER_RIGHT1.and.k.gt.2.and.k.le.nz-2) then
          Du(i,j,k) = (-8.d0*u(i,j,k+2) + 59.d0*u(i,j,k+1) - 59.d0*u(i,j,k-1) + 8.d0*u(i,j,k-2))/(86.d0*dz)
        else if (d .eq. P_STENCIL_CENTER_RIGHT.and.k.gt.1.and.k.le.nz-1) then
          Du(i,j,k) = (u(i,j,k+1) - u(i,j,k-1))/(2.d0*dz)    
        else if (d .eq. P_STENCIL_RIGHT.and.k.gt.3              ) then
          Du(i,j,k) = (48.d0*u(i,j,k) - 59.d0*u(i,j,k-1) + 8.d0*u(i,j,k-2) + 3.0d0*u(i,j,k-3))/(34.d0*dz)
	else if (d .eq. P_STENCIL_MASK) then
          Du(i,j,k) = 0.0        
	else if (d .eq. P_STENCIL_UNDEF) then
          if (ltrace2) then
            if ((bbox5==0 .and. k < gz_zmin) .or. &
                (bbox6==0 .and. k > gz_zmax)) then
              write(0,*)myid,'] d42: OK---Undefined z derivative type in ghostzone'
            else
              write(0,*)myid,'d42: PROBLEM PROBLEM   direction:',direction
              write(0,*)myid,'] d42: ### Undefined z derivative outside of ghost zone!'
              write(0,*)myid,'] d42: imask                          ',imask(i,j,k)
              write(0,*)myid,'] d42: k, nz                          ',k, nz
              write(0,*)myid,'] d42: gz_zmin                        ',gz_zmin
              write(0,*)myid,'] d42: gz_zmax                        ',gz_zmax
              write(0,*)myid,'] d42: nghostzones                    ',&
                                                          par(P_NGHOSTZONES_Z)
              write(0,*)myid,'] d42: par(bbox5), par(bbox6)         ',&
                                                par(P_BBOX5), par(P_BBOX6)
              write(0,*)myid,'] d42: bbox5, bbox6                   ',&
                                                bbox5,bbox6
              !stop
            end if
          end if
          Du(i,j,k) = 0.0
        else
          write(0,*)myid,'] d42: Illegal z derivative type   ',d
          write(0,*)myid,'] d42: imask                       ',imask(i,j,k)
          write(0,*)myid,'     at location ',x(i,j,k),y(i,j,k),z(i,j,k)
          Du(i,j,k) = 0.0
	  !stop
        end if
      end do
      end do
      end do

    else
      write(0,*)myid,'] d42:  unknown direction ',direction
    end if

    if (ltrace2) then
       write(0,*)myid,'] d42:u,du:',direction,myl2norm3d(Du, nx, ny, nz), myl2norm3d(u, nx,ny,nz)
    end if

    return
  end subroutine D42

end module HYPER_DERIVS_42
!#######################################################################
!#######################################################################
! ###
! ###
! ###                 E N D  O F  M O D U L E
! ###
!#######################################################################
!#######################################################################
