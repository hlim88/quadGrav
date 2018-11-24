!-----------------------------------------------------------------------
!
!    $Id$
!
!-----------------------------------------------------------------------

#include "cctk.h"

MODULE HYPER_DERIVS_43
use params

CONTAINS

!-----------------------------------------------------------------------
!
!    SUBROUTINE D43NoMask:
!
!    implement Strand derivative.  4th order interior, 3rd order boundary
!
!-----------------------------------------------------------------------
  subroutine D43NoMask(Du, u, direction, par)

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL, DIMENSION(:,:,:), INTENT(OUT) :: Du
    CCTK_INT :: direction


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: j, k,nx, ny, nz, shp(3)
 

    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)
    
    if (direction.eq.1) then
       do k=1, nz
        do j=1, ny	 
	  call derivs34(du(:,j,k),u(:,j,k),dx,nx)
	end do
       end do
    else if (direction.eq.2) then 
       do k=1, nz
        do j=1, nx    
	  call derivs34(du(j,:,k),u(j,:,k),dy,ny)
	end do
       end do
    else if (direction.eq.3) then 
       do k=1, ny
        do j=1, nx     
	  call derivs34(du(j,k,:),u(j,k,:),dz,nz)
	end do
	end do
    else 
      write(*,*) "D43NoMask: Wrong direction in the derivative" 
      STOP
    end if 

    return
  end subroutine D43NoMask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! BASIC ROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine derivs34(du,u,dx,Nx)
        real*8, dimension(:), intent(in) :: u
        real*8, dimension(:), intent(inout) :: du
        real*8 :: dx, h_1
        integer :: i, nx

        real*8 :: d00,d01,d02,d03, &
        &    d10,d11,d12,d13,d14,d15, &
        &    d20,d21,d22,d23,d24,d25, &
        &    d30,d31,d32,d33,d34,d35,d36, &
        &    d40,d41,d42,d43,d44,d45,d46, &
        &    d1,d2


     parameter( d00 = -1.833333333333333333333333d0 )
     parameter( d01 =  3.0d0 )
     parameter( d02 =  -1.5d0 )
     parameter( d03 =  0.333333333333333333333333d0 )

     parameter( d10 = -0.38942207148531184298d0 )
     parameter( d11 = -0.26953763903486946056d0 )
     parameter( d12 =  0.63903793765926293838d0 )
     parameter( d13 =  0.094332736084546377480d0 )
     parameter( d14 = -0.080518371580844513359d0 )
     parameter( d15 =  0.0061074083572165009295d0 )

     parameter( d20 =  0.11124996667625322721d0 )
     parameter( d21 = -0.78615310943278550936d0 )
     parameter( d22 =  0.19877943763527643222d0 )
     parameter( d23 =  0.50808067692835148792d0 )
     parameter( d24 = -0.024137062412656370601d0 )
     parameter( d25 = -0.0078199093944392672116d0 )

     parameter( d30 =  0.019051206094885019047822d0 )
     parameter( d31 =  0.026931104200732614181666d0 )
     parameter( d32 = -0.633860292039252305642283d0 )
     parameter( d33 =  0.051772670918649366462688d0 )
     parameter( d34 =  0.592764606048964306931634d0 )
     parameter( d35 = -0.054368814269840675877468d0 )
     parameter( d36 = -0.002290480954138351040607d0 )

     parameter( d40 = -0.002498706495423627386248d0 )
     parameter( d41 =  0.005463924453044550084942d0 )
     parameter( d42 =  0.087024805619019315445041d0 )
     parameter( d43 = -0.686097670431383548237962d0 )
     parameter( d44 =  0.018985530480943661987934d0 )
     parameter( d45 =  0.659895344563505072850627d0 )
     parameter( d46 = -0.082773228189705424744336d0 )

     parameter( d1 = 0.6666666666666666666666666d0 )
     parameter( d2 = -0.0833333333333333333333333d0 )

     h_1 = dx


    du(1) = (d00*u(1)+ d01*u(2)+d02*u(3)+d03*u(4))/dx
    du(nx) = -(d00*u(nx)+ d01*u(nx-1)+d02*u(nx-2)+d03*u(nx-3))/dx


    du(2) = (d10*u(1)+d11*u(2)+d12*u(3)+d13*u(4)+d14*u(5)+d15*u(6))/dx
    du(nx-1) = -(d10*u(nx)+d11*u(nx-1)+d12*u(nx-2)+d13*u(nx-3)+&
    &d14*u(nx-4)+d15*u(nx-5))/dx


    du(3) = (d20*u(1)+d21*u(2)+d22*u(3)+d23*u(4)+d24*u(5)+d25*u(6))/dx
    du(nx-2) = -(d20*u(nx)+d21*u(nx-1)+d22*u(nx-2)+d23*u(nx-3)+&
    &d24*u(nx-4)+d25*u(nx-5))/dx


    du(4) = (d30*u(1)+d31*u(2)+d32*u(3)+d33*u(4)+d34*u(5)+d35*u(6)+d36*u(7))/dx
    du(nx-3) = -(d30*u(nx)+d31*u(nx-1)+d32*u(nx-2)+d33*u(nx-3)+&
    &d34*u(nx-4)+d35*u(nx-5)+d36*u(nx-6))/dx


    du(5) = (d40*u(1)+d41*u(2)+d42*u(3)+d43*u(4)+d44*u(5)+d45*u(6)+d46*u(7))/dx
    du(nx-4) = -(d40*u(nx)+d41*u(nx-1)+d42*u(nx-2)+d43*u(nx-3)+&
    &d44*u(nx-4)+d45*u(nx-5)+d46*u(nx-6))/dx


    do i = 6, nx-5
    du(i) = (-d2*u(i-2)-d1*u(i-1)+d1*u(i+1)+d2*u(i+2))/dx
    end do
    end subroutine derivs34
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



END MODULE HYPER_DERIVS_43



  
