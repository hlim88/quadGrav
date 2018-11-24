!-----------------------------------------------------------------------
!
!    $Id$
!
!-----------------------------------------------------------------------

#include "cctk.h"

MODULE HYPER_DERIVS_43_PT
use params

CONTAINS

!-----------------------------------------------------------------------
!
!    SUBROUTINE D43NoMask_pt:
!
!-----------------------------------------------------------------------
  subroutine D43NoMask_pt(Dxu, Dyu, Dzu, u, i, j, k, par)

    CCTK_REAL, DIMENSION(:), INTENT(IN)      :: par
    CCTK_REAL, DIMENSION(:,:,:), INTENT(IN)  :: u
    CCTK_REAL                                :: Dxu, Dyu, Dzu
    CCTK_INT                                 :: i, j, k


    ! local vars
    CCTK_REAL :: dx, dy, dz
    CCTK_INT  :: nx, ny, nz, shp(3)
    
    real*8 :: d00,d01,d02,d03, &
        &    d10,d11,d12,d13,d14,d15, &
        &    d20,d21,d22,d23,d24,d25, &
        &    d30,d31,d32,d33,d34,d35,d36, &
        &    d40,d41,d42,d43,d44,d45,d46, &
        &    d1,d2


     parameter( d00 = -1.833333333333333333333333d0 )
     parameter(  d01 =  3.0d0 )
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

    shp = shape(u)
    nx = shp(1)
    ny = shp(2)
    nz = shp(3)

    dx = par(P_DX)
    dy = par(P_DY)
    dz = par(P_DZ)

    if (i .eq. 1) then
      Dxu = (d00*u(1,j,k)+ d01*u(2,j,k)+d02*u(3,j,k)+d03*u(4,j,k))/dx
    else if (i .eq. 2) then
      Dxu = (d10*u(1,j,k)+d11*u(2,j,k)+d12*u(3,j,k)+d13*u(4,j,k)&
          & +d14*u(5,j,k)+d15*u(6,j,k))/dx
    else if (i .eq. 3) then
      Dxu = (d20*u(1,j,k)+d21*u(2,j,k)+d22*u(3,j,k)+d23*u(4,j,k)&
          & +d24*u(5,j,k)+d25*u(6,j,k))/dx
    else if (i .eq. 4) then
      Dxu = (d30*u(1,j,k)+d31*u(2,j,k)+d32*u(3,j,k)+d33*u(4,j,k)&
          & +d34*u(5,j,k)+d35*u(6,j,k)+d36*u(7,j,k))/dx
    else if (i .eq. 5) then
      Dxu = (d40*u(1,j,k)+d41*u(2,j,k)+d42*u(3,j,k)+d43*u(4,j,k)&
          & +d44*u(5,j,k)+d45*u(6,j,k)+d46*u(7,j,k))/dx
    else if (i .gt. 5 .and. i .lt. nx-4) then
      Dxu = (-d2*u(i-2,j,k)-d1*u(i-1,j,k)+d1*u(i+1,j,k)+d2*u(i+2,j,k))/dx
    else if (i .eq. nx-4) then
      Dxu = -(d40*u(nx,j,k)+d41*u(nx-1,j,k)+d42*u(nx-2,j,k)+d43*u(nx-3,j,k)+&
      &d44*u(nx-4,j,k)+d45*u(nx-5,j,k)+d46*u(nx-6,j,k))/dx
    else if (i .eq. nx-3) then
      Dxu = -(d30*u(nx,j,k)+d31*u(nx-1,j,k)+d32*u(nx-2,j,k)+d33*u(nx-3,j,k)+&
      &d34*u(nx-4,j,k)+d35*u(nx-5,j,k)+d36*u(nx-6,j,k))/dx
    else if (i .eq. nx-2) then
      Dxu = -(d20*u(nx,j,k)+d21*u(nx-1,j,k)+d22*u(nx-2,j,k)+d23*u(nx-3,j,k)+&
      &d24*u(nx-4,j,k)+d25*u(nx-5,j,k))/dx
    else if (i .eq. nx-1) then
      Dxu = -(d10*u(nx,j,k)+d11*u(nx-1,j,k)+d12*u(nx-2,j,k)+d13*u(nx-3,j,k)+&
      &d14*u(nx-4,j,k)+d15*u(nx-5,j,k))/dx
    else if (i .eq. nx) then
      Dxu = -(d00*u(nx,j,k)+d01*u(nx-1,j,k)+d02*u(nx-2,j,k)+d03*u(nx-3,j,k))/dx
    else
      write(0,*) ' i out of range', i
      stop
    end if

    if (j .eq. 1) then
      Dyu = (d00*u(i,1,k)+ d01*u(i,2,k)+d02*u(i,3,k)+d03*u(i,4,k))/dy
    else if (j .eq. 2) then
      Dyu = (d10*u(i,1,k)+d11*u(i,2,k)+d12*u(i,3,k)+d13*u(i,4,k)&
          & +d14*u(i,5,k)+d15*u(i,6,k))/dy
    else if (j .eq. 3) then
      Dyu = (d20*u(i,1,k)+d21*u(i,2,k)+d22*u(i,3,k)+d23*u(i,4,k)&
          & +d24*u(i,5,k)+d25*u(i,6,k))/dy
    else if (j .eq. 4) then
      Dyu = (d30*u(i,1,k)+d31*u(i,2,k)+d32*u(i,3,k)+d33*u(i,4,k)&
          & +d34*u(i,5,k)+d35*u(i,6,k)+d36*u(i,7,k))/dy
    else if (j .eq. 5) then
      Dyu = (d40*u(i,1,k)+d41*u(i,2,k)+d42*u(i,3,k)+d43*u(i,4,k)&
          & +d44*u(i,5,k)+d45*u(i,6,k)+d46*u(i,7,k))/dy
    else if (j .gt. 5 .and. j .lt. ny-4) then
      Dyu = (-d2*u(i,j-2,k)-d1*u(i,j-1,k)+d1*u(i,j+1,k)+d2*u(i,j+2,k))/dy
    else if (j .eq. ny-4) then
      Dyu = -(d40*u(i,ny,k)+d41*u(i,ny-1,k)+d42*u(i,ny-2,k)+d43*u(i,ny-3,k)+&
      &d44*u(i,ny-4,k)+d45*u(i,ny-5,k)+d46*u(i,ny-6,k))/dy
    else if (j .eq. ny-3) then
      Dyu = -(d30*u(i,ny,k)+d31*u(i,ny-1,k)+d32*u(i,ny-2,k)+d33*u(i,ny-3,k)+&
      &d34*u(i,ny-4,k)+d35*u(i,ny-5,k)+d36*u(i,ny-6,k))/dy
    else if (j .eq. ny-2) then
      Dyu = -(d20*u(i,ny,k)+d21*u(i,ny-1,k)+d22*u(i,ny-2,k)+d23*u(i,ny-3,k)+&
      &d24*u(i,ny-4,k)+d25*u(i,ny-5,k))/dy
    else if (j .eq. ny-1) then
      Dyu = -(d10*u(i,ny,k)+d11*u(i,ny-1,k)+d12*u(i,ny-2,k)+d13*u(i,ny-3,k)+&
      &d14*u(i,ny-4,k)+d15*u(i,ny-5,k))/dy
    else if (j .eq. ny) then
      Dyu = -(d00*u(i,ny,k)+d01*u(i,ny-1,k)+d02*u(i,ny-2,k)+d03*u(i,ny-3,k))/dy
    else
      write(0,*) ' j out of range', j
      stop
    end if

    if (k .eq. 1) then
      Dzu = (d00*u(i,j,1)+ d01*u(i,j,2)+d02*u(i,j,3)+d03*u(i,j,4))/dz
    else if (k .eq. 2) then
      Dzu = (d10*u(i,j,1)+d11*u(i,j,2)+d12*u(i,j,3)+d13*u(i,j,4)&
          & +d14*u(i,j,5)+d15*u(i,j,6))/dz
    else if (k .eq. 3) then
      Dzu = (d20*u(i,j,1)+d21*u(i,j,2)+d22*u(i,j,3)+d23*u(i,j,4)&
          & +d24*u(i,j,5)+d25*u(i,j,6))/dz
    else if (k .eq. 4) then
      Dzu = (d30*u(i,j,1)+d31*u(i,j,2)+d32*u(i,j,3)+d33*u(i,j,4)&
          & +d34*u(i,j,5)+d35*u(i,j,6)+d36*u(i,j,7))/dz
    else if (k .eq. 5) then
      Dzu = (d40*u(i,j,1)+d41*u(i,j,2)+d42*u(i,j,3)+d43*u(i,j,4)&
          & +d44*u(i,j,5)+d45*u(i,j,6)+d46*u(i,j,7))/dz
    else if (k .gt. 5 .and. k .lt. nz-4) then
      Dzu = (-d2*u(i,j,k-2)-d1*u(i,j,k-1)+d1*u(i,j,k+1)+d2*u(i,j,k+2))/dz
    else if (k .eq. nz-4) then
      Dzu = -(d40*u(i,j,nz)+d41*u(i,j,nz-1)+d42*u(i,j,nz-2)+d43*u(i,j,nz-3)+&
      &d44*u(i,j,nz-4)+d45*u(i,j,nz-5)+d46*u(i,j,nz-6))/dz
    else if (k .eq. nz-3) then
      Dzu = -(d30*u(i,j,nz)+d31*u(i,j,nz-1)+d32*u(i,j,nz-2)+d33*u(i,j,nz-3)+&
      &d34*u(i,j,nz-4)+d35*u(i,j,nz-5)+d36*u(i,j,nz-6))/dz
    else if (k .eq. nz-2) then
      Dzu = -(d20*u(i,j,nz)+d21*u(i,j,nz-1)+d22*u(i,j,nz-2)+d23*u(i,j,nz-3)+&
      &d24*u(i,j,nz-4)+d25*u(i,j,nz-5))/dz
    else if (k .eq. nz-1) then
      Dzu = -(d10*u(i,j,nz)+d11*u(i,j,nz-1)+d12*u(i,j,nz-2)+d13*u(i,j,nz-3)+&
      &d14*u(i,j,nz-4)+d15*u(i,j,nz-5))/dz
    else if (k .eq. nz) then
      Dzu = -(d00*u(i,j,nz)+d01*u(i,j,nz-1)+d02*u(i,j,nz-2)+d03*u(i,j,nz-3))/dz
    else
      write(0,*) ' k out of range', k
      stop
    end if

    return
    end subroutine D43NoMask_pt
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



END MODULE HYPER_DERIVS_43_PT



  
