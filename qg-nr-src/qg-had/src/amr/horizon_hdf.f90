!-----------------------------------------------------------------------
!   hdf.f
!
!   Mijan Huq 
!   University of Texas at Austin.
!
!   Created Wed Jun 15 14:18:24 CDT 1994 MFH
!
!   Modified Tue Oct 25 10:42:36 CDT 1994 MFH
!      . Changed from reading in psi.hdf to g11.hdf g12.hdf ...
!
!   Modified Fri Mar 10 16:47:42 CST 1995 MFH
!      . Modified to read hdf data with  gf3_read_hdf_f77 rather
!        than the brief reader which does not store coordinates.
!
!   Contains routines that will read in metric data from hdf files.
!
!   Modified Fri May 16 2003 by Matt Anderson
!   to read in hdf5 data files.
!-----------------------------------------------------------------------
      subroutine HDFMetricData(nx, ny, nz, level,grx, gry, grz, &
                g11, g12, g13, g22, g23, g33, &
                k11, k12, k13, k22, k23, k33,msk,hx,hy,hz)

! This module contains all necessary modules for hdf5 interface
!        use hdf5

        implicit none

        integer     nx, ny, nz, level

!       Coordinates to be figured out from grid.dat
        real*8      grx(nx), gry(ny), grz(nz)
       
!       Metric components to be read in from hdf files
        real*8      g11(nx,ny,nz) ,g12(nx,ny,nz), g13(nx,ny,nz)
        real*8      g22(nx,ny,nz) ,g23(nx,ny,nz), g33(nx,ny,nz)

!       Extrinsic curvature components in from hdf files
        real*8      k11(nx,ny,nz) ,k12(nx,ny,nz), k13(nx,ny,nz)
        real*8      k22(nx,ny,nz) ,k23(nx,ny,nz), k33(nx,ny,nz)

!       Mask
        real*8      msk(nx,ny,nz)

!		  Grid spaceings
        real*8      hx, hy, hz, min, max
	real*8 min1, min2, min3, max1, max2, max3

!       Local variables

        integer     i, j, k
        integer     data_dims(3)

!        integer(hid_t) :: file_id     ! file identifier
!        integer(hid_t) :: dset_id     ! dataset identifier
!        integer(hid_t) :: dataspace   ! dataspace identifier
!        integer(hid_t) :: memspace    ! memspace identifier
!        integer :: error ! error flag
!
! debug
        real*8 x,y,z,r
        real*8 kd22,gd23,t1,t5,t7,gd33
        real*8 t3,kd11,gd22,t4,kd13,t12
        real*8 kd33,gd13,kd12,gd11,t10,t2
        real*8 t14,gd12,kd23

        hx = grx(2)-grx(1)
        hy = gry(2)-gry(1)
        hz = grz(2)-grz(1)
!        open(file='grid.dat',unit=11)
!        read(11,*)min1,max1
!        hx = (max1 - min1)/float(NX-1)
!        read(11,*)min2,max2
!        hy = (max2 - min2)/float(NY-1)
!        read(11,*)min3,max3
!        hz = (max3 - min3)/float(NZ-1)
!        close(unit=11)
!        write(6,*)'(hx,hy,hz)=','(',hx,hy,hz,')'
        
!        do i = 1, nx
!          grx(i) = min1 + (i-1)*hx
!        end do
!        do i = 1, ny
!          gry(i) = min2 + (i-1)*hy
!        end do
!        do i = 1, nz
!          grz(i) = min3 + (i-1)*hz
!        end do

        if (.true.) then 
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if (sqrt(grx(i)**2 + gry(j)**2 + grz(k)**2) .lt. &
                    (0.5)) then
                msk(i,j,k) = 0.0
                g11(i,j,k) = 1.0
                g12(i,j,k) = 0.0
                g13(i,j,k) = 0.0
                g22(i,j,k) = 1.0
                g23(i,j,k) = 0.0
                g33(i,j,k) = 1.0
                k11(i,j,k) = 0.0
                k12(i,j,k) = 0.0
                k13(i,j,k) = 0.0
                k22(i,j,k) = 0.0
                k23(i,j,k) = 0.0
                k33(i,j,k) = 0.0
              else
                x = grx(i) 
                y = gry(j) 
                z = grz(k) 
                r = sqrt(x*x+y*y+z*z)
               t1 = x**2
               t2 = r**2
               gd11 = 1.D0+2.D0*t1/t2/r
               t2 = r**2
               gd12 = 2.D0*x*y/t2/r
               t2 = r**2
               gd13 = 2.D0*x*z/t2/r
               t1 = y**2
               t2 = r**2
               gd22 = 1.D0+2.D0*t1/t2/r
               t2 = r**2
               gd23 = 2.D0*y*z/t2/r
               t1 = z**2
               t2 = r**2
               gd33 = 1.D0+2.D0*t1/t2/r
               t1 = x**2
               t4 = r**2
               t7 = t4**2
               t14 = sqrt((r+2.D0)/r)
               kd11 = -2.D0*(t1+2.D0*t1*r-t4*r)/t7/r/t14
               t1 = 1.D0/r
               t4 = r**2
               t5 = t4**2
               t10 = sqrt(1.D0+2.D0*t1)
               kd12 = -2.D0*(t1+2.D0)*x*y/t5/t10
               t1 = 1.D0/r
               t4 = r**2
               t5 = t4**2
               t10 = sqrt(1.D0+2.D0*t1)
               kd13 = -2.D0*(t1+2.D0)*x*z/t5/t10
               t1 = 1.D0/r
               t3 = y**2
               t5 = r**2
               t7 = t5**2
               t12 = sqrt(1.D0+2.D0*t1)
               kd22 = -2.D0*((t1+2.D0)*t3-t5)/t7/t12
               t1 = 1.D0/r
               t4 = r**2
               t5 = t4**2
               t10 = sqrt(1.D0+2.D0*t1)
               kd23 = -2.D0*(t1+2.D0)*y*z/t5/t10
               t1 = 1.D0/r
               t3 = z**2
               t5 = r**2
               t7 = t5**2
               t12 = sqrt(1.D0+2.D0*t1)
               kd33 = -2.D0*((t1+2.D0)*t3-t5)/t7/t12
                msk(i,j,k) = 1.0
                g11(i,j,k) = gd11
                g12(i,j,k) = gd12
                g13(i,j,k) = gd13
                g22(i,j,k) = gd22
                g23(i,j,k) = gd23
                g33(i,j,k) = gd33
                k11(i,j,k) = kd11
                k12(i,j,k) = kd12
                k13(i,j,k) = kd13
                k22(i,j,k) = kd22
                k23(i,j,k) = kd23
                k33(i,j,k) = kd33
              end if
            end do
          end do
        end do        
        end if

        return
      end
