#include "cctk.h"
   Module OUTPUT
        implicit none

        CONTAINS
	
	
	subroutine output3D(var,name,x1,x2,x3,time,myproc)
	implicit none
	CCTK_REAL, dimension(:,:,:) :: var
	CCTK_REAL, dimension(:) :: x1,x2,x3
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
          write(0,*)'output3D >>> can not allocate tempcoord'
          write(0,*)'             size: ',nx,ny,nz
        end if

	
	tempcoord(1:nx) = x1
	tempcoord(nx+1:nx+ny) = x2
	tempcoord(nx+ny+1:nx+ny+nz) = x3
	
        write (g11f,'(A,1i3)') name,myproc
	if(myproc.lt.0) g11f = name
	
	
	ret=gft_out_full(g11f,time,(/nx,ny,nz/),cnames,3, &
                                tempcoord,var)
				
				
        if (allocated(tempcoord)) then
	  deallocate(tempcoord)
        else
          write(0,*)'output3D >>> tempcoord not allocated as expected!'
        end if
				
	end subroutine output3D
	

!!!----------------- 2D ----------------------
	subroutine output2D(var,name,x1,x2,time,myproc)
	implicit none
	CCTK_REAL, dimension(:,:) :: var
	CCTK_REAL, dimension(:) :: x1,x2
	CCTK_REAL :: time
	character(len=*) :: name
	
	CCTK_INT myproc, istat
		
	real*8, allocatable, dimension(:) :: tempcoord
	integer :: nx, ny, nz, ret, gft_out_full
	character(20):: g11f
	character(32):: cnames
	
	nx=size(var(:,1))
	ny=size(var(1,:))
	
	cnames = 'x|y'
	
	allocate(tempcoord(nx+ny),STAT=istat)
        if (istat .ne. 0) then
          write(0,*)'output2D >>> can not allocate tempcoord'
          write(0,*)'             size: ',nx,ny
        end if

	
	tempcoord(1:nx) = x1
	tempcoord(nx+1:nx+ny) = x2
	
        write (g11f,'(A,1i3)') name,myproc
	if(myproc.lt.0) g11f = name
	
	
	ret=gft_out_full(g11f,time,(/nx,ny/),cnames,2, &
                                tempcoord,var)
				
				
        if (allocated(tempcoord)) then
          deallocate(tempcoord)
        else
          write(0,*)'output2D >>> tempcoord not allocated as expected!'
        end if

				
	end subroutine output2D				

!!!----------------- 1D ----------------------
	subroutine output1D(var,name,x1,time,myproc)
	implicit none
	CCTK_REAL, dimension(:) :: var
	CCTK_REAL, dimension(:) :: x1
	CCTK_REAL :: time
	character(len=*) :: name
	
	CCTK_INT myproc
		
	integer :: nx, ny, nz, ret, gft_out_full
	character(20):: g11f
	character(32):: cnames
	
	nx=size(var(:))
	
	cnames = 'x'
		
        write (g11f,'(A,1i3)') name,myproc
	if(myproc.lt.0) g11f = name
	
	ret=gft_out_full(g11f,time,(/nx/),cnames,2, &
                                x1,var)
				
				
				
	end subroutine output1D		

   End Module OUTPUT
				
