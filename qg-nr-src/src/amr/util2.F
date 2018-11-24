#include "cctk.h"

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  turok_spergel3d:                                                          cc
cc                  Initializes initial data to generalized                   cc
cc                  Turok-Spergel solution in 3D.                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine turok_spergel3d( chi, chit, h,
     *           minx, miny, minz,nx,ny,nz,epsilon, Delta,ex,ey)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     chi(nx,ny,nz), chit(nx,ny,nz),h, ex, ey, r, d,
     *           minx, miny, minz, epsilon, Delta
      real(kind=8)     x,y,z
      integer    i,j,k

      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'turok_spergel3d:   h = ',h
      end if

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               r           = sqrt( ex*x**2 + ey*y**2 + z**2 )
               chi(i,j,k)  = 2.d0 * epsilon * atan( r / Delta )
               chit(i,j,k) = 2.d0 * epsilon * r / ( Delta**2 + r**2 )
            end do
         end do
      end do


      return
      end    ! END: turok_spergel3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  toroid3d:                                                                 cc
cc                   Initializes field to a toroid pulse.                     cc
cc                                                                            cc
cc field(x,y,z) = a*exp( -(sqrt(e_x*x^2+e_y*y^2)-r0)^2/d^2 )*exp(-z^2/d^2)    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine toroid3d( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z
      integer    i,j,k

      logical     ltrace
      parameter ( ltrace = .false. )

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               field(i,j,k) =
     *              a*exp(-(sqrt(ex*(x-xc)**2+ey*(y-yc)**2)-r)**2/d**2 )
     *              *exp( -z**2 / d**2 )
            end do
         end do
      end do

      return
      end    ! END: toroid3d


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3das:                                                             cc
cc                   Initializes field to an antisymmetric Gaussian pulse.    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3das( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z
      integer    i,j,k

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
c              field(i,j,k) = a*(x+y)*(z+y)*(x+z)*exp( -(sqrt(
c              field(i,j,k) = a*x**2*exp( -(sqrt(
               field(i,j,k) = a*x   *exp( -(sqrt(
     *             ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2)-r)**2/d**2 )
            end do
         end do
      end do

      return
      end    ! END: gaussian3das

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  tanhball:                                                                 cc
cc                   Initializes field to a Tanh ball.                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tanhball( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z, rhat
      integer    i,j,k

      if (d.eq.0) then
         write(*,*) 'tanhball: PROBLEM: Gaussian width cannot be zero'
         write(*,*) 'tanhball: setting d=1'
         d = 1.d0
      end if

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x    = minx + h*(i-1.d0)
               rhat = sqrt(ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2)
               field(i,j,k) = a*
     *            *0.5d0*( tanh( ( rhat+R)/d ) + 1.d0 )
     *            *0.5d0*( tanh( (-rhat+R)/d ) + 1.d0 )
            end do
         end do
      end do

      return
      end    ! END: tanhball


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3d:                                                               cc
cc                   Initializes field to a Gaussian pulse.                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3d( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z
      integer    i,j,k

      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'gaussian3d:   h = ',h
         write(*,*) 'gaussian3d: minx = ',minx
         write(*,*) 'gaussian3d: miny = ',miny
         write(*,*) 'gaussian3d: minz = ',minz
         write(*,*) 'gaussian3d:    a = ',a
         write(*,*) 'gaussian3d:   ex = ',ex
         write(*,*) 'gaussian3d:   ey = ',ey
         write(*,*) 'gaussian3d:    r = ',r
         write(*,*) 'gaussian3d:    d = ',d
         write(*,*) 'gaussian3d:   xc = ',xc
         write(*,*) 'gaussian3d:   yc = ',yc
         write(*,*) 'gaussian3d:   zc = ',zc
      end if

      if (d.eq.0) then
         write(*,*) 'gaussian3d: PROBLEM: Gaussian width cannot be zero'
         write(*,*) 'gaussian3d: setting d=1'
         d = 1.d0
      end if

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               field(i,j,k) = a*exp( -(sqrt(
     *             ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2)-r)**2/d**2 )
c    *             ex*(x-xc)**2 + ey*(y-yc)**2            )-r)**2/d**2 )
            end do
         end do
      end do

      if (ltrace) then
         write(*,*) 'gaussian3d: maxx = ',x
         write(*,*) 'gaussian3d: maxy = ',y
         write(*,*) 'gaussian3d: maxz = ',z
      end if

      return
      end    ! END: gaussian3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3d_add_x:    Initializes field to x deriv of Gaussian pulse.      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3d_add_x( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z, P
      integer    i,j,k

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               P = sqrt( ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2 )
               field(i,j,k) = field(i,j,k) +a*exp( -(P-r)**2/d**2 )
     *             * (-2.d0)*((P-r)/d**2)*ex*(x-xc)
            end do
         end do
      end do

      return
      end    ! END: gaussian3d_add_x
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3d_add_y:    Initializes field to x deriv of Gaussian pulse.      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3d_add_y( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z, P
      integer    i,j,k

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               P = sqrt( ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2 )
               field(i,j,k) = field(i,j,k) +a*exp( -(P-r)**2/d**2 )
     *             * (-2.d0)*((P-r)/d**2)*ey*(y-yc)
            end do
         end do
      end do

      return
      end    ! END: gaussian3d_add_y
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3d_add_z:    Initializes field to x deriv of Gaussian pulse.      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3d_add_z( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z, P
      integer    i,j,k

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               P = sqrt( ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2 )
               field(i,j,k) = field(i,j,k) +a*exp( -(P-r)**2/d**2 )
     *             * (-2.d0)*((P-r)/d**2)*(z-zc)
            end do
         end do
      end do

      return
      end    ! END: gaussian3d_add_z
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3d_x:    Initializes field to x deriv of Gaussian pulse.          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3d_x( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z, P
      integer    i,j,k

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               P = sqrt( ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2 )
               field(i,j,k) = a*exp( -(P-r)**2/d**2 )
     *             * (-2.d0)*((P-r)/d**2)*ex*(x-xc)
            end do
         end do
      end do

      return
      end    ! END: gaussian3d_x
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3d_y:    Initializes field to y deriv of Gaussian pulse.          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3d_y( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z, P
      integer    i,j,k

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               P = sqrt( ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2 )
               field(i,j,k) = a*exp( -(P-r)**2/d**2 )
     *             * (-2.d0)*((P-r)/d**2)*ey*(y-yc)
            end do
         end do
      end do

      return
      end    ! END: gaussian3d_y

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3d_z:    Initializes field to x deriv of Gaussian pulse.          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3d_z( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d,xc,yc,zc)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz, xc,yc,zc
      real(kind=8)     x,y,z, P
      integer    i,j,k

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               P = sqrt( ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2 )
               field(i,j,k) = a*exp( -(P-r)**2/d**2 )
     *             * (-2.d0)*((P-r)/d**2)*(z-zc)
            end do
         end do
      end do

      return
      end    ! END: gaussian3d_z

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3d_r:                                                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3d_r( field, h, minx, miny, minz,nx,ny,nz,
     *           a,ex,ey,r,d)
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), h, a, ex, ey, r, d,
     *           minx, miny, minz
      real(kind=8)     x,y,z, tempr
      integer    i,j,k

      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'gaussian3d_r:   h = ',h
         write(*,*) 'gaussian3d_r: minx = ',minx
         write(*,*) 'gaussian3d_r: miny = ',miny
         write(*,*) 'gaussian3d_r: minz = ',minz
         write(*,*) 'gaussian3d_r:    a = ',a
         write(*,*) 'gaussian3d_r:   ex = ',ex
         write(*,*) 'gaussian3d_r:   ey = ',ey
         write(*,*) 'gaussian3d_r:    r = ',r
         write(*,*) 'gaussian3d_r:    d = ',d
      end if

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               tempr = sqrt(ex*x**2 + ey*y**2 + z**2)
               field(i,j,k) = a*exp( -(tempr-r)**2/d**2 )
     *                         *(-2.d0)*(tempr-r)/d**2
            end do
         end do
      end do

      return
      end    ! END: gaussian3d_r

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  vector_abs1D:                                                             cc
cc                  Take the absolute value of all elements in an array.      cc
cc                  Works in place.                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vector_abs1D( field, nx)
      implicit      none
      integer       nx
      real(kind=8)  field(nx)
      integer       i

      do i = 1, nx
         field(i) = abs(field(i))
      end do

      return
      end    ! END: vector_abs1D

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_r:                                                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_r( r, h, minx, miny, minz,nx,ny,nz )
      implicit    none
      integer    nx, ny, nz
      real(kind=8)     r(nx,ny,nz), h,
     *           minx, miny, minz
      real(kind=8)     x,y,z
      integer    i,j,k

      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'load_r: h    = ',h
         write(*,*) 'load_r: minx = ',minx
         write(*,*) 'load_r: miny = ',miny
         write(*,*) 'load_r: minz = ',minz
         write(*,*) 'load_r:   nx = ',nx
         write(*,*) 'load_r:   ny = ',ny
         write(*,*) 'load_r:   nz = ',nz
      end if

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               x        = minx + h*(i-1.d0)
               y        = miny + h*(j-1.d0)
               z        = minz + h*(k-1.d0)
               r(i,j,k) = sqrt(x**2 + y**2 + z**2)
            end do
         end do
      end do

      if (ltrace) then
         write(*,*) 'load_r: r(nx,ny,nz) =',r(nx,ny,nz)
      end if

      return
      end    ! END: load_r


c-------------------------------------------------------------------------------
      real(kind=8) function cube_root(x)
!
!       Copyright 1994, Miles Ellis, Ivor Philips and Tom Lahey
!       Adapted from web: http://rsc.anu.edu.au/~harry/COURSES/MATHMETH/node53.html
!
c-------------------------------------------------------------------------------
      implicit none

      ! Function to calculate the cube root of a positive real number

      REAL*8 x

      REAL*8       epsilon
      PARAMETER  ( epsilon = 1.0D-20 )

c     write(*,*) 27.d0**(1/3)

      IF ( ABS(x) .lt. epsilon ) THEN
         cube_root = 0.d0
      ELSE IF (x .lt. 0) THEN
         cube_root = -EXP(LOG(-x)/3.d00)
      ELSE
         cube_root = EXP(LOG(x)/3.d00)
      END IF

      end ! END: function cube_root()

c-------------------------------------------------------------------------------
      subroutine  field_dump_info(u, nx, ny, nz)
c-------------------------------------------------------------------------------
      implicit none
      integer  nx,        ny,   nz
      real(kind=8)   u(nx,ny,nz)
      integer  i,         j,    k
      real(kind=8)   l2norm, max, min
      include       'largesmall.inc'
      integer        myid
      integer        proc_return_myid
      external       proc_return_myid
      logical        isanan
      external       isanan


      myid     = proc_return_myid()
      max      = -LARGENUMBER
      min      =  LARGENUMBER
      l2norm   = 0.d0

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               l2norm = l2norm + u(i,j,k)**2
               if (u(i,j,k).lt.min) min = u(i,j,k)
               if (u(i,j,k).gt.max) max = u(i,j,k)
            end do
         end do
      end do

      l2norm = sqrt(l2norm/nx/ny/nz)

      write(*,90)myid,'field_dump_info: l2norm, max,min=',l2norm,max,min
!90   format(A,3F18.3)
 90   format( '[',I4,'] ',A,3G18.6)

      if (isanan(l2norm)) then
         write(*,90)myid,'field_dump_info: Field has a NaN.@@@@'
         !write(*,*) 'field_dump_info: Field has a NaN.@@@@'
         !stop
      end if

      return
      end

c-------------------------------------------------------------------------------
      subroutine  field_dump_infoB(u, nx, ny, nz, name)
c-------------------------------------------------------------------------------
      implicit none
      integer  nx,        ny,   nz
      real(kind=8)   u(nx,ny,nz)
      character*(*)  name
      integer  i,         j,    k
      integer        myid, proc_return_myid
      external             proc_return_myid
      real(kind=8)   l2norm, max, min
      include       'largesmall.inc'
      logical        isanan
      external       isanan


      myid     = proc_return_myid()

      max      = -LARGENUMBER
      min      =  LARGENUMBER
      l2norm   = 0.d0

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               l2norm = l2norm + u(i,j,k)**2
               if (u(i,j,k).lt.min) min = u(i,j,k)
               if (u(i,j,k).gt.max) max = u(i,j,k)
            end do
         end do
      end do

      l2norm = sqrt(l2norm/nx/ny/nz)

      write(*,91) myid,name,l2norm,max,min
 91   format('[',I4,'] ',A10,' l2nm/mx/mn ',3G14.6)
 90   format(A,3F18.8)

      if (isanan(l2norm)) then
         write(*,*) 'field_dump_info: Field has a NaN.@@@@'
         !stop
      end if

      return
      end
c-------------------------------------------------------------------------------
      subroutine  field_dump_stats(u, nx, ny, nz, name)
c-------------------------------------------------------------------------------
      implicit none
      integer  nx,        ny,   nz
      real(kind=8)   u(nx,ny,nz)
      character*(*)  name
      integer  i,         j,    k
      integer        myid, proc_return_myid
      external             proc_return_myid
      real(kind=8)   l2norm, max, min
      include       'largesmall.inc'
      logical        isanan
      external       isanan
      integer       numaboveone,numzero,numneg,numpos,numone

      numaboveone = 0
      numone      = 0
      numzero     = 0
      numneg      = 0
      numpos      = 0

      myid     = proc_return_myid()

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (u(i,j,k).gt.1) numaboveone = numaboveone + 1
               if (abs(u(i,j,k)-1.d0).lt.SMALLNUMBER) numone = numone+1
               if (abs(u(i,j,k))     .lt.SMALLNUMBER) numzero=numzero+1
               if (u(i,j,k).lt.0) numneg = numneg + 1
               if (u(i,j,k).gt.0) numpos = numpos + 1
            end do
         end do
      end do


      write(*,91) myid,name,numaboveone,numone,numzero,numneg,numpos
      write(*,92) myid,name,nx,ny,nz
 91   format('[',I4,'] ',A10,' >1  =1  0  <0  >0 ',5I8)
 92   format('[',I4,'] ',A10,' nx/y/z: ',3I8)

      return
      end

c-------------------------------------------------------------------------------
      subroutine  field_dump_info_nj(u, nx, ny, nz, inan,jnan,knan)
c-------------------------------------------------------------------------------
      implicit none
      integer  nx,        ny,   nz
      real(kind=8)   u(nx,ny,nz)
      integer  i,         j,    k, inan,jnan,knan
      real(kind=8)   l2norm, max, min
      include       'largesmall.inc'
      !real(kind=8)   LARGENUMBER
      !parameter    ( LARGENUMBER = 9.d98 )
      logical        isanan
      external       isanan


      max      = -LARGENUMBER
      min      =  LARGENUMBER
      l2norm   = 0.d0

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (isanan(u(i,j,k))) then
                  write (*,*) "nan found at", i,j,k 
                  inan =i
                  jnan = j
                  knan = k
               end if
               l2norm = l2norm + u(i,j,k)**2
               if (u(i,j,k).lt.min) min = u(i,j,k)
               if (u(i,j,k).gt.max) max = u(i,j,k)
            end do
         end do
      end do

      l2norm = sqrt(l2norm/nx/ny/nz)

      write(*,90) 'field_dump_info_nj:l2norm, max, min =',l2norm,max,min
 90   format(A,3F18.3)

      if (isanan(l2norm)) then
         write(*,*) 'field_dump_info_nj: Field has a NaN.@@@@'
         !stop
      end if

      return
      end

c-------------------------------------------------------------------------------
      subroutine  field_dump_info1D(u, nx, name)
c-------------------------------------------------------------------------------
      implicit none
      integer  nx
      real(kind=8)   u(nx)
      character*(*)  name
      integer  i,         j
      integer        myid, proc_return_myid
      external             proc_return_myid
      real(kind=8)   l2norm, max, min
      include       'largesmall.inc'
      logical        isanan
      external       isanan


      myid     = proc_return_myid()

      max      = -LARGENUMBER
      min      =  LARGENUMBER
      l2norm   = 0.d0

            do i = 1, nx
               l2norm = l2norm + u(i)**2
               if (u(i).lt.min) min = u(i)
               if (u(i).gt.max) max = u(i)
            end do

      l2norm = sqrt(l2norm/nx)

      write(*,91) myid,name,l2norm,max,min
 91   format('[',I4,'] ',A10,' l2nm/mx/mn ',3F18.9)
 90   format(A,3F18.8)

      if (isanan(l2norm)) then
         write(*,*) 'field_dump_info1D: Field has a NaN.@@@@'
      end if

      return
      end

c-------------------------------------------------------------------------------
      subroutine  field_dump_info2D(u, nx, ny, name)
c-------------------------------------------------------------------------------
      implicit none
      integer  nx,        ny
      real(kind=8)   u(nx,ny)
      character*(*)  name
      integer  i,         j
      integer        myid, proc_return_myid
      external             proc_return_myid
      real(kind=8)   l2norm, max, min
      include       'largesmall.inc'
      logical        isanan
      external       isanan


      myid     = proc_return_myid()

      max      = -LARGENUMBER
      min      =  LARGENUMBER
      l2norm   = 0.d0

         do j = 1, ny
            do i = 1, nx
               l2norm = l2norm + u(i,j)**2
               if (u(i,j).lt.min) min = u(i,j)
               if (u(i,j).gt.max) max = u(i,j)
            end do
         end do

      l2norm = sqrt(l2norm/nx/ny)

      write(*,91) myid,name,l2norm,max,min
 91   format('[',I4,'] ',A10,' l2nm/mx/mn ',3F18.9)
 90   format(A,3F18.8)

      if (isanan(l2norm)) then
         write(*,*) 'field_dump_info2D: Field has a NaN.@@@@'
         !stop
      end if

      return
      end

c-------------------------------------------------------------------------------
      real(kind=8) function myl2norm(u, nx)
c-------------------------------------------------------------------------------
      implicit none
      integer  nx
      real(kind=8)   u(nx)
      integer  i
      logical        isanan, isainf
      external       isanan, isainf
      logical        ltrace
      parameter  (   ltrace = .false. )


      myl2norm = 0.d0

      if (nx.le.0) then
         write(*,*) 'myl2norm: ERROR'
      end if

            do i = 1, nx
               if (ltrace.and.isanan(u(i))) then
                  write(*,*)'myl2norm: Nan found at: ',i,u(i)
               end if
               if (ltrace.and.isainf(u(i))) then
                  write(*,*)'myl2norm: Inf found at: ',i,u(i)
               end if
               myl2norm = myl2norm + u(i)**2
            end do

      myl2norm = sqrt(myl2norm/nx)
      return
      end

c-------------------------------------------------------------------------------
      real(kind=8) function myl2norm3d(u, nx, ny, nz)
c-------------------------------------------------------------------------------
      implicit none
      integer  nx,        ny,   nz
      real(kind=8)   u(nx,ny,nz)
      integer  i,         j,    k
      logical        isanan, isainf
      external       isanan, isainf
      logical        ltrace
      parameter  (   ltrace = .false. )


      myl2norm3d = 0.d0

      if ((nx.le.0).or.(ny.le.0).or.(nz.le.0)) then
         write(*,*) 'myl2norm3d: ERROR'
      end if

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (ltrace.and.isanan(u(i,j,k))) then
                  write(*,*)'myl2norm3d: Nan found at: ',i,j,k,
     *                        u(i,j,k)
               end if
               if (ltrace.and.isainf(u(i,j,k))) then
                  write(*,*)'myl2norm3d: Inf found at: ',i,j,k,
     *                        u(i,j,k)
               end if
               myl2norm3d = myl2norm3d + u(i,j,k)**2
            end do
         end do
      end do

      myl2norm3d = sqrt(myl2norm3d/nx/ny/nz)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  dump5by5:                                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dump5by5( field,ax,ay,az,nx,ny,nz )
      implicit    none
      integer     ax,ay,az,  nx,ny,nz
      real(kind=8)      field(nx,ny,nz)
      integer     i, j

      write(*,*) '.....'
      do j = 1, 5
         do i = 1, 5
            write(*,10) j,i,field(1+ax-1,j+ay-1,i+az-1),
     *                      field(2+ax-1,j+ay-1,i+az-1),
     *                      field(3+ax-1,j+ay-1,i+az-1),
     *                      field(4+ax-1,j+ay-1,i+az-1),
     *                      field(5+ax-1,j+ay-1,i+az-1)
         end do
      end do
 10   format('x,',I2,',',I2,': ',5F9.4)

      return
      end    ! END: dump5by5

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_min1d:                                                               cc
cc          finds minimum of 1d vector                                        cc
cc          If multiple regions of the same minimum value,                    cc
cc          it finds the one with the largest span.                           cc
cc             min   = minimum value found                                    cc
cc             min_i = index of minimum where first found                     cc
cc             span  = number of successive points at this minimum            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_min1d( u1, min, min_i, span, nx)
      implicit    none
      integer     min_i,   span, nx
      real(kind=8)      u1(nx),  min
      integer     min_i_tmp,   span_tmp
      real(kind=8)      min_tmp
      integer     i
      !real(kind=8)      LARGENUMBER
      !parameter ( LARGENUMBER = 9.d98 )
      include       'largesmall.inc'
      logical     within_span, within_span_tmp
      logical     ltrace
      parameter ( ltrace      = .false. )

      if (ltrace) then
         write(*,*) 'find_min1d: nx     = ',nx
         write(*,*) 'find_min1d: u1(1)  = ',u1(1)
         write(*,*) 'find_min1d: u1(nx) = ',u1(nx)
      end if

      if (nx.le.0) then
            min   = LARGENUMBER
            min_i = 0
            span  = 0
            return
      end if

      min     = LARGENUMBER
      min_tmp = LARGENUMBER + 1

      do i = 1, nx
         if (ltrace) write(*,*) '      i, u1(i)=',i,u1(i)
         if ( u1(i) .lt. min ) then
            min             = u1(i)
            min_i           = i
            span            = 1
            within_span     = .true.
            within_span_tmp = .false.
         else if ( (u1(i).eq.min) .and. within_span) then
            span            = span + 1
         else if (        (u1(i).eq.min)
     *             .and. .not.within_span
     *             .and. .not.within_span_tmp) then
            min_tmp         = min
            min_i_tmp       = i
            span_tmp        = 1
            within_span_tmp = .true.
         else if ( (u1(i).eq.min) .and. within_span_tmp) then
            span_tmp        = span_tmp + 1
         else                !   (u1(i) .ne. min)
            within_span     = .false.
            within_span_tmp = .false.
            if ( (min .eq. min_tmp) .and. (span_tmp.gt.span) ) then
               span  = span_tmp
               min_i = min_i_tmp
            end if
         end if
      end do

      if (ltrace) then
         write(*,*) 'find_min1d: min   = ', min
         write(*,*) 'find_min1d: min_i = ', min_i
         write(*,*) 'find_min1d: span  = ', span
      end if

      return
      end      ! END: find_min1d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_max:                                                                 cc
cc          finds maximum of some vector                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8)  function  find_max( u, num )
      implicit      none
      integer       num
      real(kind=8)  u(num)
      integer       i

      find_max = u(1)
      do i = 2, num
         if ( u(i) .gt. find_max ) find_max = u(i)
      end do

      return
      end      ! END: find_max

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_min:                                                                 cc
cc          finds minimum of some vector                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8)  function  find_min( u, num )
      implicit      none
      integer       num
      real(kind=8)  u(num)
      integer       i

      find_min = u(1)
      do i = 2, num
         if ( u(i) .lt. find_min ) find_min = u(i)
      end do

      return
      end      ! END: find_min

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_min_index1D:                                                         cc
cc          finds index of minimum of a 1D vector                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function  find_min_index1D( u, num )
      implicit      none
      integer       num
      real(kind=8)  u(num)
      integer       i
      real(kind=8)  mintmp
      logical     ltrace
      parameter ( ltrace = .false. )

      mintmp           = u(1)
      find_min_index1D = 1
      
      do i = 2, num
         if ( u(i) .lt. mintmp ) then
            find_min_index1D = i
            mintmp           = u(i)
         end if
      end do

      if (ltrace) then
         do i = 1, num
            write(*,*) 'find_min_index1D: i, u(i) =',i, u(i)
         end do
            write(*,*) 'find_min_index1D: Index of min:',
     *                        find_min_index1D
      end if
      

      return
      end      ! END: find_min_index1d


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  compute_disallowed:                                                       cc
cc             Compute signatures and Laplacians for disallowed points        cc
cc             in a given cluster.                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_disallowed( flagarray,
     *                      sigi,        sigj,       sigk,
     *                      lapi,        lapj,       lapk,
     *                      bbminx,      bbminy,     bbminz,
     *                      bbmaxx,      bbmaxy,     bbmaxz,
     *                      nx,          ny,         nz,
     *                      disallowedlogical                    )
      implicit    none
      integer     nx,         ny,         nz,
     *            bbminx,     bbminy,     bbminz,
     *            bbmaxx,     bbmaxy,     bbmaxz
      real(kind=8)      flagarray(nx,ny,nz), 
     *            sigi(nx),   sigj(ny),   sigk(nz),
     *            lapi(nx),   lapj(ny),   lapk(nz)
      logical     disallowedlogical
      integer     i,          j,          k
      integer     lengthx,    lengthy,    lengthz
      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'compute_disallowed: bbminx/maxx: ',bbminx,bbmaxx
         write(*,*) 'compute_disallowed: bbminy/maxy: ',bbminy,bbmaxy
         write(*,*) 'compute_disallowed: bbminz/maxz: ',bbminz,bbmaxz
         write(*,*) 'compute_disallowed: nx/y/z:      ',nx,ny,nz
      end if

      !
      ! Zero out signature:
      !   (only need to zero where bounding box exists:
      !
      do i = bbminx, bbmaxx
         sigi(i) = 0.d0
      end do
      do j = bbminy, bbmaxy
         sigj(j) = 0.d0
      end do
      do k = bbminz, bbmaxz
         sigk(k) = 0.d0
      end do


      disallowedlogical = .false.

      do k = bbminz, bbmaxz
         do j = bbminy, bbmaxy
            do i = bbminx, bbmaxx
               if ( flagarray(i,j,k) .lt. 0 ) then
                  disallowedlogical = .true.
                  sigi(i)    = sigi(i) + 1.d0
                  sigj(j)    = sigj(j) + 1.d0
                  sigk(k)    = sigk(k) + 1.d0
               end if
            end do
         end do
      end do

      !
      ! Only compute lapi/j/k if disallowed points are here:
      !
      if (disallowedlogical) then
         lengthx = bbmaxx-bbminx+1
         lengthy = bbmaxy-bbminy+1
         lengthz = bbmaxz-bbminz+1
         !write(*,*) lengthx,lengthy, lengthz
         call  compute_lap(lapi(bbminx), sigi(bbminx), lengthx )
         call  compute_lap(lapj(bbminy), sigj(bbminy), lengthy )
         call  compute_lap(lapk(bbminz), sigk(bbminz), lengthz )
      end if

      if (ltrace) then
         write(*,*) 'compute_disallowed: disallowedlogical: ',
     *                                   disallowedlogical
      end if

      return
      end     ! END: compute_disallowed

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_bbox:                                                                cc
cc             Find bounding box for a binary [0,1] array.                    cc
cc             Note: finds bounding box within *existing* box                 cc
cc                   defined by [bbminx,bbmaxx][][].                          cc
cc             Likewise computes:                                             cc
cc                   sigi, sigj, sigk  ---Signature lines                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_bbox( flagarray,
     *                      sigi,        sigj,       sigk,
     *                      bbminx,      bbminy,     bbminz,
     *                      bbmaxx,      bbmaxy,     bbmaxz,
     *                      nx,          ny,         nz,
     *                      efficiency                           )
      implicit    none
      integer     nx,         ny,         nz,
     *            bbminx,     bbminy,     bbminz,
     *            bbmaxx,     bbmaxy,     bbmaxz
      real(kind=8)      flagarray(nx,ny,nz), 
     *            sigi(nx),   sigj(ny),   sigk(nz),
     *            efficiency
      integer     i,          j,          k,
     *            bbminx_tmp, bbminy_tmp, bbminz_tmp,
     *            bbmaxx_tmp, bbmaxy_tmp, bbmaxz_tmp
      logical     ltrace
      parameter ( ltrace = .false. )

c     if (ltrace) then
c        write(*,*) 'find_bbox: bbminx, bbmaxx =',bbminx,bbmaxx
c        write(*,*) 'find_bbox: bbminy, bbmaxy =',bbminy,bbmaxy
c        write(*,*) 'find_bbox: bbminz, bbmaxz =',bbminz,bbmaxz
c        write(*,*) 'find_bbox: nx,ny,nz=',nx,ny,nz
c     end if
      !
      ! Zero out signature:
      !   (only need to zero where maximal bounding box exists:
      !
      do i = bbminx, bbmaxx
         sigi(i) = 0.d0
      end do
      do j = bbminy, bbmaxy
         sigj(j) = 0.d0
      end do
      do k = bbminz, bbmaxz
         sigk(k) = 0.d0
      end do

      !
      ! Define temporary bounding box:
      !
      bbminx_tmp = nx
      bbminy_tmp = ny
      bbminz_tmp = nz
      bbmaxx_tmp = 1
      bbmaxy_tmp = 1
      bbmaxz_tmp = 1

      efficiency = 0.d0

      do k = bbminz, bbmaxz
         do j = bbminy, bbmaxy
            do i = bbminx, bbmaxx
               if ( flagarray(i,j,k) .gt. 0 ) then
                  efficiency = efficiency + 1.d0
                  sigi(i)    = sigi(i) + 1.d0
                  sigj(j)    = sigj(j) + 1.d0
                  sigk(k)    = sigk(k) + 1.d0
                  if (i.gt.bbmaxx_tmp) bbmaxx_tmp = i
                  if (i.lt.bbminx_tmp) bbminx_tmp = i
                  if (j.gt.bbmaxy_tmp) bbmaxy_tmp = j
                  if (j.lt.bbminy_tmp) bbminy_tmp = j
                  if (k.gt.bbmaxz_tmp) bbmaxz_tmp = k
                  if (k.lt.bbminz_tmp) bbminz_tmp = k
               end if
            end do
         end do
c        if (ltrace) write(*,*) k, sigk(k)
      end do

      !
      ! replace input bbox with  the bounding box found here:
      !
      bbminx = bbminx_tmp
      bbminy = bbminy_tmp
      bbminz = bbminz_tmp
      bbmaxx = bbmaxx_tmp
      bbmaxy = bbmaxy_tmp
      bbmaxz = bbmaxz_tmp

c     if (ltrace) write(*,*) 'find_bbox: numpoints = ',efficiency

      efficiency = efficiency / (bbmaxx - bbminx + 1)
     *                        / (bbmaxy - bbminy + 1)
     *                        / (bbmaxz - bbminz + 1)

      if (ltrace) then
         write(*,*) 'find_bbox: efficiency     = ',efficiency
         write(*,*) 'find_bbox: bbminx, bbmaxx =',bbminx,bbmaxx
         write(*,*) 'find_bbox: bbminy, bbmaxy =',bbminy,bbmaxy
         write(*,*) 'find_bbox: bbminz, bbmaxz =',bbminz,bbmaxz
         do i = bbminx, bbmaxx
            write(*,*) i, sigi(i)
         end do
      end if

      return
      end     ! END: find_bbox

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_absmin3d:                                                            cc
cc             Find minimum of absolute value of a 3d field.                  cc
cc             Return absolute value of maximum and index to point.           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_absmin3d( field,   max,
     *                          max_i,   max_j, max_k,
     *                          nx,      ny,    nz     )
      implicit    none
      integer     max_i, max_j, max_k,
     *            nx,    ny,    nz
      real(kind=8)      field(nx,ny,nz), max
      integer     i,     j,     k
      real(kind=8)      tmp
      include       'largesmall.inc'
      !real(kind=8)   LARGENUMBER
      !parameter    ( LARGENUMBER = 9.d98 )

      max = LARGENUMBER
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               tmp = abs(field(i,j,k))
               if ( tmp .lt. max) then
                  max   = tmp
                  max_i = i
                  max_j = j
                  max_k = k
               end if
            end do
         end do
      end do

      return
      end      ! END: find_absmin3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_absmax3d:                                                            cc
cc             Find maximum of absolute value of a 3d field.                  cc
cc             Return absolute value of maximum and index to point.           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_absmax3d( field,   max,
     *                          max_i,   max_j, max_k,
     *                          nx,      ny,    nz     )
      implicit    none
      integer     max_i, max_j, max_k,
     *            nx,    ny,    nz
      real(kind=8)      field(nx,ny,nz), max
      integer     i,     j,     k
      real(kind=8)      tmp

      max = 0.d0
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               tmp = abs(field(i,j,k))
               if ( tmp .gt. max) then
                  max   = tmp
                  max_i = i
                  max_j = j
                  max_k = k
               end if
            end do
         end do
      end do

      return
      end      ! END: find_absmax3d


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interp_from_parentB:                                                      cc
cc                                                                            cc
cc           1) copy data points to corresponding points in grid function     cc
cc           2) interpolate in between these points                           cc
cc                                                                            cc
cc    NB: This version receives a black hole mask and acts accordingly        cc
cc        near the excision surface.                                          cc
cc    NB: bi/j/k: denote the number of points included in u_c "before"        cc
cc                the fine grid data starts in each direction, because        cc
cc                for higher order interpolation, one sends extra points.     cc
cc                Generally, these will be either 0 or 1.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp_from_parentB( u_c, u_f, chr, mask,
     *                         nxc, nyc, nzc,
     *                         mini,minj,mink,
     *                         maxi,maxj,maxk,
     *                         bi,  bj,  bk,
     *                         nxf, nyf, nzf, r, gw)
      implicit      none
      integer       nxc,nyc,nzc, nxf,nyf,nzf, r, gw
      integer       mini,minj,mink, bi,bj,bk
      integer       maxi,maxj,maxk
      real(kind=8)  u_c(nxc,nyc,nzc), u_f(nxf,nyf,nzf),
     *              chr(nxf,nyf,nzf), mask(nxf,nyf,nzf)
      include      'chr.inc'
      integer       ic,jc,kc, if,jf,kf, i
      real(kind=8)  temp, mytmp, mytmp2
      real(kind=8)  field_i_plus2, field_i_plus, field_i,
     *              field_i_minus, field_i_minus2
      real(kind=8)  field_j_plus2, field_j_plus, field_j,
     *              field_j_minus, field_j_minus2
      real(kind=8)  field_k_plus2, field_k_plus, field_k,
     *              field_k_minus, field_k_minus2
      real(kind=8)   interp_cubicb
      external       interp_cubicb
      logical        aligned_i, aligned_j, aligned_k, unmasked
      !
      ! Number which represents data that has been masked off
      ! (from the other grid, in this case, a parent grid):
      !
      include       'largesmall.inc'
      !
      logical        ltrace
      parameter (    ltrace  = .false. )
      logical        ltrace2
      parameter (    ltrace2 = .false. )

 8    format(A,3I3,' | ',3I3)
 9    format(A,3I3,' | ',3I3,' | ',I2)

      if (ltrace) then
         write(*,*) 'interp_from_parent: '
         write(*,*) 'interp_from_parent: bi/j/k:  ',bi,bj,bk
         write(*,*) 'interp_from_parent: nx/y/zf: ',nxf,nyf,nzf
         write(*,*) 'interp_from_parent: nx/y/zc: ',nxc,nyc,nzc
         write(*,*) 'interp_from_parent: mini/j/k:',mini,minj,mink
         write(*,*) 'interp_from_parent: maxi/j/k:',maxi,maxj,maxk
      end if
      if ( maxi.lt.mini .or. maxj.lt.minj .or. maxk .lt.mink) then
         write(*,*) 'interp_from_parent: '
         write(*,*) 'interp_from_parent: bi/j/k:  ',bi,bj,bk
         write(*,*) 'interp_from_parent: nx/y/zf: ',nxf,nyf,nzf
         write(*,*) 'interp_from_parent: nx/y/zc: ',nxc,nyc,nzc
         write(*,*) 'interp_from_parent: mini/j/k:',mini,minj,mink
         write(*,*) 'interp_from_parent: maxi/j/k:',maxi,maxj,maxk
         write(*,*) 'interp_from_parent: max < min, returning'
         return
      end if

      !
      ! Do we have sufficient data to interpolate only symmetrically?
      !
      ic = (mini - mini)/r + 1 + bi
      if ( ic .le. 1 ) then
         write(*,*) 'interp_from_parent: Not sufficient for small x'
      end if
      ic = (maxi - mini)/r + 1 + bi
      if ( ic .ge. nxc ) then
         write(*,*) 'interp_from_parent: Not sufficient for large x',ic
         write(*,*) 'interp_from_parent: bi/j/k:  ',bi,bj,bk
         write(*,*) 'interp_from_parent: nx/y/zf: ',nxf,nyf,nzf
         write(*,*) 'interp_from_parent: nx/y/zc: ',nxc,nyc,nzc
         write(*,*) 'interp_from_parent: mini/j/k:',mini,minj,mink
         write(*,*) 'interp_from_parent: maxi/j/k:',maxi,maxj,maxk
      end if
      jc = (minj - minj)/r + 1 + bj
      if ( jc .le. 1 ) then
         write(*,*) 'interp_from_parent: Not sufficient for small y'
      end if
      jc = (maxj - minj)/r + 1 + bj
      if ( jc .ge. nyc ) then
         write(*,*) 'interp_from_parent: Not sufficient for large y',jc
      end if
      kc = (mink - mink)/r + 1 + bk
      if ( kc .le. 1 ) then
         write(*,*) 'interp_from_parent: Not sufficient for small z'
      end if
      kc = (maxk - mink)/r + 1 + bk
      if ( kc .ge. nzc ) then
         write(*,*) 'interp_from_parent: Not sufficient for large z',kc
      end if

      !call field_dump_info(mask,nxf,nyf,nzf)
      !call field_dump_info(u_c,nxc,nyc,nzc)
      if (ltrace) call field_dump_info(u_f,nxf,nyf,nzf)
      !call field_out3d(mask,0.d0, 'mask',0.d0,1.d0,0.d0,1.d0, 
!    *                       0.d0,1.d0,nxf,nyf,nzf,0)

      do kf = mink, maxk, r
         kc = (kf - mink)/r + 1 + bk
         do jf = minj, maxj, r
            jc = (jf - minj)/r + 1 + bj
            do if = mini, maxi-r, r
               ic = (if - mini)/r + 1 + bi
               !
               ! copy points in common to both grids
               !
               unmasked =  NINT(mask(if,jf,kf)) .ge. 0
               if (unmasked .and. u_c(ic,jc,kc).lt.LARGE_M_SMALL) then
                  u_f(if,jf,kf) = u_c(ic,jc,kc)
               end if
               !
               !.....interpolate across x-direction
               !
               do i = 1, r-1
                  !
                  temp     = (1.d0*i) / r
                  unmasked =  NINT(mask(if+i,jf,kf)) .ge. 0
                  if (unmasked) then
                     if (ic .gt. 1 .and. NINT(mask(max(if-r,1),jf,kf))
     *                  .ge.0.and. u_c(ic-1,jc,kc).lt.LARGE_M_SMALL)then
                        field_i_minus   = u_c(ic-1,jc,kc)
                     else 
                        field_i_minus   = LARGENUMBER
                     end if
                     if ( NINT(mask(if,jf,kf)).ge.0 
     *                       .and. u_c(ic,jc,kc).lt.LARGE_M_SMALL ) then
                        field_i         = u_c(ic,  jc,kc)
                     else 
                        field_i         = LARGENUMBER
                     end if
                     if(ic+1.le.nxc.and.NINT(mask(if+r,jf,kf)).ge.0
     *                      .and. u_c(ic+1,jc,kc).lt.LARGE_M_SMALL )then
                        field_i_plus    = u_c(ic+1,jc,kc)
                     else 
                        field_i_plus    = LARGENUMBER
                     end if
                   if(ic+2.le.nxc .and.NINT(mask(min(if+2*r,nxf),jf,kf))
     *                 .ge.0 .and. u_c(ic+2,jc,kc).lt.LARGE_M_SMALL)then
                        field_i_plus2   = u_c(ic+2,jc,kc)
                     else 
                        field_i_plus2   = LARGENUMBER
                     end if
                     mytmp2          =interp_cubicb(   field_i_plus2,
     *                                                 field_i_plus,
     *                                                 field_i,
     *                                                 field_i_minus,
     *                                                 1.d0+temp,
     *                                              3.d0,2.d0,1.d0,0.d0)
                     if (mytmp2 .lt. 1.d80) u_f(if+i,jf,kf) = mytmp2
                  else
                     if ( u_f(if+i,jf,kf) .ge. 1.d80) then
                       if (ltrace) then
                       write(*,*)'interp_from_parentb: Point is masked'
                       write(*,*)'but has a bad value: ',u_f(if+i,jf,kf)
                       write(*,*)'Setting to zero:'
                       end if
                       u_f(if+i,jf,kf) = 0.d0
                     end if
                  end if
                  if ( u_f(if+i,jf,kf) .ge. 1.d80 .and. ltrace) then
                        write(*,*) 'interp_from_parentb: Bad u_f in X',
     *                   u_f(if+i,jf,kf)
                        write(*,*) '     i/j/kf        = ', if,jf,kf
                        write(*,*) '     i/j/kc        = ', ic,jc,kc
                        write(*,*) '     nx/y/zf       = ',nxf,nyf,nzf
                        write(*,*) '     nx/y/zc       = ',nxc,nyc,nzc
                       write(*,*) '     mini/j/k      = ',mini,minj,mink
                       write(*,*) '     maxi/j/k      = ',maxi,maxj,maxk
                        write(*,*) '     i             = ',i
                        write(*,*) 'field_i_plus2 = ',field_i_plus2
                        write(*,*) 'field_i_plus  = ',field_i_plus 
                        write(*,*) 'field_i       = ',field_i
                        write(*,*) 'field_i_minus = ',field_i_minus
                        write(*,*) 'mask if+2r    = ',mask(if+2*r,jf,kf)
                        write(*,*) 'mask if+r     = ',mask(if+  r,jf,kf)
                        write(*,*) 'mask if       = ',mask(if,    jf,kf)
                        write(*,*) 'mask if-r     = ',mask(if-  r,jf,kf)
                        write(*,*) 'u_c(ic+2)     = ',u_c(ic+2,jc,kc)
                        write(*,*) 'u_c(ic+1)     = ',u_c(ic+1,jc,kc)
                        write(*,*) 'u_c(ic  )     = ',u_c(ic  ,jc,kc)
                        write(*,*) 'u_c(ic-1)     = ',u_c(ic-1,jc,kc)
                        write(*,*) 'u_f(if  )     = ',u_f(if  ,jf,kf)
                        write(*,*) '     mytmp         = ',mytmp
                        write(*,*) '     mytmp2        = ',mytmp2
                        write(*,*) '     unmasked      = ',unmasked
                  end if
               end do
            end do
            if       = maxi
            unmasked = NINT(mask(if,jf,kf)) .ge. 0
            ic       = (if - mini)/r + 1 + bi
            if (unmasked .and. u_c(ic,jc,kc).lt.LARGE_M_SMALL) then
               u_f(if,jf,kf) = u_c(ic,jc,kc)
            end if
         end do
      end do
      if (ltrace) call field_dump_info(u_f,nxf,nyf,nzf)
      if (ltrace) write(*,*) 'interp_from_parent:  Done w/ x'

      !
      !.....interpolate in y-direction
      !
      do kf = mink, maxk, r
         kc = (kf - mink)/r + 1 + bk
         do jf = minj, maxj-r, r
            jc = (jf - minj)/r + 1 + bj
            do if = mini, maxi
               ic = (if - mini)/r + 1 + bi
               aligned_i =  (ic-1-bi)*r+mini .eq. if
               do i = 1, r-1
                  temp     = (1.d0*i) / r
                  unmasked =  NINT(mask(if,jf+i,kf)) .ge. 0
                  if (unmasked) then
                  !
                  ! Establish the j_minus value:
                  !
                  if (jf .gt. minj .and
     *               .NINT(mask(if,max(jf-r,1),kf)).ge.0)then
                     field_j_minus   = u_f(if,jf-r,kf)
                  else if (aligned_i .and. jc.gt.1 .and.
     *                            u_c(ic,jc-1,kc).lt.LARGE_M_SMALL) then
                     field_j_minus   = u_c(ic,jc-1,  kc)
                  else if (ic+2.le.nxc .and. jc.gt.1.and.ic.gt.1) then
                     field_j_minus   = interp_cubicb(
     *                                             u_c(ic+2,jc-1,kc),
     *                                             u_c(ic+1,jc-1,kc),
     *                                             u_c(ic,  jc-1,kc),
     *                                             u_c(ic-1,jc-1,kc),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                  else if (ic.gt.2 .and. jc.gt.1.and.ic+1.le.nxc) then
                     field_j_minus   = interp_cubicb(
     *                                             u_c(ic+1,jc-1,kc),
     *                                             u_c(ic,  jc-1,kc),
     *                                             u_c(ic-1,jc-1,kc),
     *                                             u_c(ic-2,jc-1,kc),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                  else
                      !write(*,8)'intrp_parent: e',if,jf+i,kf,ic,jc,kc
                     field_j_minus   = LARGENUMBER
                  end if
                  !
                  ! Establish the j_plus2 value:
                  !
                  if (jf+2*r .le. maxj .and.
     *                    NINT(mask(if,min(jf+2*r,nyf),kf)).ge.0) then
                     field_j_plus2   = u_f(if,jf+2*r,kf)
                  else if (aligned_i .and. jc+2.le.nyc .and.
     *                            u_c(ic,jc+2,kc).lt.LARGE_M_SMALL) then
                     field_j_plus2   = u_c(ic,jc+2,  kc)
                  else if (ic+2.le.nxc .and. jc+2.le.nyc) then
                     field_j_plus2   = interp_cubicb(
     *                                             u_c(ic+2,jc+2,kc),
     *                                             u_c(ic+1,jc+2,kc),
     *                                             u_c(ic,  jc+2,kc),
     *                                             u_c(ic-1,jc+2,kc),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                  else if (ic.gt.2 .and. jc+2.le.nyc) then
                     field_j_plus2   = interp_cubicb(
     *                                             u_c(ic+1,jc+2,kc),
     *                                             u_c(ic,  jc+2,kc),
     *                                             u_c(ic-1,jc+2,kc),
     *                                             u_c(ic-2,jc+2,kc),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                  else
                   !write(*,8)'intrp_parent: d',if,jf+i,kf,ic,jc,kc
                     field_j_plus2   = LARGENUMBER
                  end if
                  !
                  ! In case one of these points is masked, 
                  ! communicate it to interpolation routine:
                  !
                  if (NINT(mask(if,jf+r,kf)) .ge. 0) then
                     field_j_plus = u_f(if,jf+r,kf)
                  else
                     field_j_plus = LARGENUMBER
                  end if
                  if (NINT(mask(if,jf,  kf)) .ge. 0) then
                     field_j      = u_f(if,jf,  kf)
                  else
                     field_j      = LARGENUMBER
                  end if
                  !
                  ! Do the interpolation:
                  !
                  u_f(if,jf+i,kf) = interp_cubicb(  field_j_plus2,
     *                                              field_j_plus,
     *                                              field_j,
     *                                              field_j_minus,
     *                                              1.d0+temp,
     *                                              3.d0,2.d0,1.d0,0.d0)
                  else
                     ! Point is masked off:
                     if ( u_f(if,jf+i,kf) .ge. 1.d80) then
                       if (ltrace) then
                       write(*,*)'interp_from_parentb:Y Point is masked'
                       write(*,*)'but has a bad value: ',u_f(if,jf+i,kf)
                       write(*,*)'Setting to zero:'
                       end if
                       u_f(if,jf+i,kf) = 0.d0
                     end if
                  end if
                  if ( u_f(if,jf+i,kf) .ge. 1.d80 .and. ltrace) then
                        write(*,*) 'interp_from_parentb: Y: Bad u_f',
     *                   u_f(if,jf+i,kf),if,jf,kf,ic,jc,kc
                       write(*,*)'interp_from_parentb:Y:field_j_plus2:',
     *                   field_j_plus2
                       write(*,*)'interp_from_parentb:Y:field_j_plus: ',
     *                   field_j_plus
                       write(*,*)'interp_from_parentb:Y:field_j:      ',
     *                   field_j
                       write(*,*)'interp_from_parentb:Y:field_j_minus:',
     *                   field_j_minus
                        write(*,*) 'mask jf+2r    = ',mask(if,jf+2*r,kf)
                        write(*,*) 'mask jf+r     = ',mask(if,jf+  r,kf)
                        write(*,*) 'mask jf       = ',mask(if,jf    ,kf)
                        write(*,*) 'mask jf-r     = ',mask(if,jf-  r,kf)
                        write(*,*) 'u_f(jf+2r)    = ', u_f(if,jf+2*r,kf)
                        write(*,*) 'u_f(jf+ r)    = ', u_f(if,jf+  r,kf)
                        write(*,*) 'u_f(jf   )    = ', u_f(if,jf    ,kf)
                        write(*,*) 'u_f(jf- r)    = ', u_f(if,jf-  r,kf)
                        write(*,*) '     unmasked      = ',unmasked
                        write(*,*) 'Setting to zero....'
                       u_f(if,jf+i,kf) = 0.d0
                  end if
               end do
            end do
         end do
      end do
      if (ltrace) call field_dump_info(u_f,nxf,nyf,nzf)
      if (ltrace) write(*,*) 'interp_from_parent:  Done w/ y'

      !
      !.....interpolate in z-direction
      !
      do kf = mink, maxk-r, r
         kc = (kf - mink)/r + 1 + bk
         do jf = minj, maxj
            jc = (jf - minj)/r + 1 + bj
            aligned_j =  (jc-1-bj)*r+minj .eq. jf
            do if = mini, maxi
               ic = (if - mini)/r + 1 + bi
               aligned_i =  (ic-1-bi)*r+mini .eq. if
               do i = 1, r-1
                  temp     = (1.d0*i) / r
                  unmasked =  NINT(mask(if,jf,kf+i)) .ge. 0
                  if (unmasked) then
                  if(ltrace2)write(*,*)'interp: Z: i/j/kf: ',if,jf,kf
                  if(ltrace2)write(*,*)'interp: Z: i/j/kc: ',ic,jc,kc
                  if(ltrace2)write(*,*)'interp: Z: aligned_i/j: ',
     .              aligned_i,aligned_j
                  !
                  ! Establish the k_minus value:
                  !
                  if (kf-r.ge.mink .and
     *               .NINT(mask(if,jf,max(kf-r,1))).ge.0)then
                     field_k_minus   = u_f(if,jf,kf-r)
                     if(ltrace2)write(*,*)'interp: field_k_minus A'
                  else if (aligned_i .and. aligned_j .and. kc.gt.1) then
                     field_k_minus   = u_c(ic,jc,kc-1)
                     if(ltrace2)write(*,*)'interp: field_k_minus B'
                  else if (aligned_j.and.ic+2.le.nxc.and.kc.gt.1
     .                              .and.ic.gt.1) then
                     field_k_minus   = interp_cubicb(
     *                                             u_c(ic+2,jc,kc-1),
     *                                             u_c(ic+1,jc,kc-1),
     *                                             u_c(ic,  jc,kc-1),
     *                                             u_c(ic-1,jc,kc-1),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     if(ltrace2)write(*,*)'interp: field_k_minus C'
                  else if (aligned_j.and.ic+1.le.nxc.and.kc.gt.1
     .                              .and.ic.gt.2) then
                     field_k_minus   = interp_cubicb(
     *                                             u_c(ic+1,jc,kc-1),
     *                                             u_c(ic,  jc,kc-1),
     *                                             u_c(ic-1,jc,kc-1),
     *                                             u_c(ic-2,jc,kc-1),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     if(ltrace2)write(*,*)'interp: field_k_minus D'
                  else if (aligned_i.and.jc+2.le.nyc.and.kc.gt.1
     .                              .and.jc.gt.1) then
                     field_k_minus   = interp_cubicb(
     *                                             u_c(ic,jc+2,kc-1),
     *                                             u_c(ic,jc+1,kc-1),
     *                                             u_c(ic,jc  ,kc-1),
     *                                             u_c(ic,jc-1,kc-1),
     *                                 1.d0+0.5d0*(jf-(jc-1-bj)*r-minj),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     if(ltrace2)write(*,*)'interp: field_k_minus E'
                  else if (aligned_i.and.jc+1.le.nyc.and.kc.gt.1
     .                              .and.jc.gt.2) then
                     field_k_minus   = interp_cubicb(
     *                                             u_c(ic,jc+1,kc-1),
     *                                             u_c(ic,jc  ,kc-1),
     *                                             u_c(ic,jc-1,kc-1),
     *                                             u_c(ic,jc-2,kc-1),
     *                                 2.d0+0.5d0*(jf-(jc-1-bj)*r-minj),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     if(ltrace2)write(*,*)'interp: field_k_minus F'
                  else if (jc+2.le.nyc.and.jc.gt.1) then
                     if(ltrace2)write(*,*)'interp: field_k_minus G'
                       !
                       ! Now we have to interpolate just to get
                       ! the value we need to do the interpolation in z:
                       !
                       if (ic+2 .le. nxc.and.ic.gt.1) then
                           field_j_plus2   = interp_cubicb(
     *                                             u_c(ic+2,jc+2,kc-1),
     *                                             u_c(ic+1,jc+2,kc-1),
     *                                             u_c(ic  ,jc+2,kc-1),
     *                                             u_c(ic-1,jc+2,kc-1),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                           field_j_plus    = interp_cubicb(
     *                                             u_c(ic+2,jc+1,kc-1),
     *                                             u_c(ic+1,jc+1,kc-1),
     *                                             u_c(ic  ,jc+1,kc-1),
     *                                             u_c(ic-1,jc+1,kc-1),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                           field_j         = interp_cubicb(
     *                                             u_c(ic+2,jc,kc-1),
     *                                             u_c(ic+1,jc,kc-1),
     *                                             u_c(ic  ,jc,kc-1),
     *                                             u_c(ic-1,jc,kc-1),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                           field_j_minus   = interp_cubicb(
     *                                             u_c(ic+2,jc-1,kc-1),
     *                                             u_c(ic+1,jc-1,kc-1),
     *                                             u_c(ic  ,jc-1,kc-1),
     *                                             u_c(ic-1,jc-1,kc-1),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                       else if (ic.gt.2) then
                           field_j_plus2   = interp_cubicb(
     *                                             u_c(ic+1,jc+2,kc-1),
     *                                             u_c(ic,  jc+2,kc-1),
     *                                             u_c(ic-1,jc+2,kc-1),
     *                                             u_c(ic-2,jc+2,kc-1),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                           field_j_plus    = interp_cubicb(
     *                                             u_c(ic+1,jc+1,kc-1),
     *                                             u_c(ic,  jc+1,kc-1),
     *                                             u_c(ic-1,jc+1,kc-1),
     *                                             u_c(ic-2,jc+1,kc-1),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                           field_j         = interp_cubicb(
     *                                             u_c(ic+1,jc,kc-1),
     *                                             u_c(ic,  jc,kc-1),
     *                                             u_c(ic-1,jc,kc-1),
     *                                             u_c(ic-2,jc,kc-1),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                           field_j_minus   = interp_cubicb(
     *                                             u_c(ic+1,jc-1,kc-1),
     *                                             u_c(ic,  jc-1,kc-1),
     *                                             u_c(ic-1,jc-1,kc-1),
     *                                             u_c(ic-2,jc-1,kc-1),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                       else 
                         write(*,8)'intrp_parent: P',if,jf,kf+i,ic,jc,kc
                       end if
                       field_k_minus   = interp_cubicb(
     *                                             field_j_plus2,
     *                                             field_j_plus,
     *                                             field_j,
     *                                             field_j_minus,
     *                                 1.d0+0.5d0*(jf-(jc-1-bj)*r-minj),
     *                                             3.d0,2.d0,1.d0,0.d0)
                       !if(ltrace2)write(*,*)'interp:field_k_minux to100'
                       !field_k_minus   = 100.d0
                  else
                      !write(*,8)'intrp_parent: F',if,jf+i,kf,ic,jc,kc
                      if(ltrace2)write(*,*)'interp:no field_k_minus val'
                       field_k_minus   = LARGENUMBER
                  end if
                  !
                  ! Establish the k_plus2 value:
                  !
                  if (kf+2*r .le. maxk .and.
     *                NINT(mask(if,jf,min(kf+2*r,nzf))).ge.0) then
                     field_k_plus2   = u_f(if,jf,kf+2*r)
                     if(ltrace2)write(*,*)'interp: field_k_plus2 A'
                  else if (aligned_i.and.aligned_j.and.kc+2.le.nzc) then
                     field_k_plus2   = u_c(ic,jc,kc+2)
                     if(ltrace2)write(*,*)'interp: field_k_plus2 B'
                  else if(aligned_j.and.ic+2.le.nxc.and.kc+2.le.nzc
     .                                  .and.ic.gt.1)then
                     field_k_plus2   = interp_cubicb(
     *                                             u_c(ic+2,jc,kc+2),
     *                                             u_c(ic+1,jc,kc+2),
     *                                             u_c(ic,  jc,kc+2),
     *                                             u_c(ic-1,jc,kc+2),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     if(ltrace2)write(*,*)'interp: field_k_plus2 C'
                  else if(aligned_j.and.ic+1.le.nxc.and.kc+2.le.nzc
     .                             .and.ic.gt.2)then
                        field_k_plus2   = interp_cubicb(
     *                                             u_c(ic+1,jc,kc+2),
     *                                             u_c(ic  ,jc,kc+2),
     *                                             u_c(ic-1,jc,kc+2),
     *                                             u_c(ic-2,jc,kc+2),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     if(ltrace2)write(*,*)'interp: field_k_plus2 D'
                  else if(aligned_i.and.jc+2.le.nyc.and.kc+2.le.nzc
     .                             .and.jc.gt.1)then
                        field_k_plus2   = interp_cubicb(
     *                                             u_c(ic,jc+2,kc+2),
     *                                             u_c(ic,jc+1,kc+2),
     *                                             u_c(ic,jc  ,kc+2),
     *                                             u_c(ic,jc-1,kc+2),
     *                                 1.d0+0.5d0*(jf-(jc-1-bj)*r-minj),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     if(ltrace2)write(*,*)'interp: field_k_plus2 E'
                  else if(aligned_i.and.jc+1.le.nyc.and.kc+2.le.nzc
     .                             .and.jc.gt.2)then
                        field_k_plus2   = interp_cubicb(
     *                                             u_c(ic,jc+1,kc+2),
     *                                             u_c(ic,jc  ,kc+2),
     *                                             u_c(ic,jc-1,kc+2),
     *                                             u_c(ic,jc-2,kc+2),
     *                                 2.d0+0.5d0*(jf-(jc-1-bj)*r-minj),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     if(ltrace2)write(*,*)'interp: field_k_plus2 F'
                  else if (jc+2.le.nyc.and.kc+2.le.nzc.and.ic.gt.1
     .                                .and.ic+2.le.nxc) then
                     if(ltrace2)write(*,*)'interp: field_k_plus2 G'
                       !
                       ! Now we have to interpolate just to get
                       ! the value we need to do the interpolation in z:
                       !
                     if (ic+2 .le. nxc) then
                        field_j_plus2   = interp_cubicb(
     *                                             u_c(ic+2,jc+2,kc+2),
     *                                             u_c(ic+1,jc+2,kc+2),
     *                                             u_c(ic  ,jc+2,kc+2),
     *                                             u_c(ic-1,jc+2,kc+2),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j_plus    = interp_cubicb(
     *                                             u_c(ic+2,jc+1,kc+2),
     *                                             u_c(ic+1,jc+1,kc+2),
     *                                             u_c(ic  ,jc+1,kc+2),
     *                                             u_c(ic-1,jc+1,kc+2),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j         = interp_cubicb(
     *                                             u_c(ic+2,jc,kc+2),
     *                                             u_c(ic+1,jc,kc+2),
     *                                             u_c(ic  ,jc,kc+2),
     *                                             u_c(ic-1,jc,kc+2),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j_minus   = interp_cubicb(
     *                                             u_c(ic+2,jc-1,kc+2),
     *                                             u_c(ic+1,jc-1,kc+2),
     *                                             u_c(ic  ,jc-1,kc+2),
     *                                             u_c(ic-1,jc-1,kc+2),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     else if (ic.gt.2) then
                        field_j_plus2   = interp_cubicb(
     *                                             u_c(ic+1,jc+2,kc+2),
     *                                             u_c(ic,  jc+2,kc+2),
     *                                             u_c(ic-1,jc+2,kc+2),
     *                                             u_c(ic-2,jc+2,kc+2),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j_plus    = interp_cubicb(
     *                                             u_c(ic+1,jc+1,kc+2),
     *                                             u_c(ic,  jc+1,kc+2),
     *                                             u_c(ic-1,jc+1,kc+2),
     *                                             u_c(ic-2,jc+1,kc+2),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j         = interp_cubicb(
     *                                             u_c(ic+1,jc,kc+2),
     *                                             u_c(ic,  jc,kc+2),
     *                                             u_c(ic-1,jc,kc+2),
     *                                             u_c(ic-2,jc,kc+2),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j_minus   = interp_cubicb(
     *                                             u_c(ic+1,jc-1,kc+2),
     *                                             u_c(ic,  jc-1,kc+2),
     *                                             u_c(ic-1,jc-1,kc+2),
     *                                             u_c(ic-2,jc-1,kc+2),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     else
                      write(*,8)'intrp_parent: Q',if,jf,kf+i,ic,jc,kc
                     end if
                     if(ltrace2)write(*,*)'interp:interp in Y',
     .                   'to get field_k_plus2'
                     field_k_plus2   = interp_cubicb(
     *                                             field_j_plus2,
     *                                             field_j_plus,
     *                                             field_j,
     *                                             field_j_minus,
     *                                 1.d0+0.5d0*(jf-(jc-1-bj)*r-minj),
     *                                             3.d0,2.d0,1.d0,0.d0)
                  else if (jc.gt.2.and.jc+1.le.nyc) then
                     if (ic+2 .le. nxc) then
                        field_j_plus    = interp_cubicb(
     *                                             u_c(ic+2,jc+1,kc+2),
     *                                             u_c(ic+1,jc+1,kc+2),
     *                                             u_c(ic  ,jc+1,kc+2),
     *                                             u_c(ic-1,jc+1,kc+2),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j         = interp_cubicb(
     *                                             u_c(ic+2,jc,kc+2),
     *                                             u_c(ic+1,jc,kc+2),
     *                                             u_c(ic  ,jc,kc+2),
     *                                             u_c(ic-1,jc,kc+2),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j_minus   = interp_cubicb(
     *                                             u_c(ic+2,jc-1,kc+2),
     *                                             u_c(ic+1,jc-1,kc+2),
     *                                             u_c(ic  ,jc-1,kc+2),
     *                                             u_c(ic-1,jc-1,kc+2),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j_minus2  = interp_cubicb(
     *                                             u_c(ic+2,jc-2,kc+2),
     *                                             u_c(ic+1,jc-2,kc+2),
     *                                             u_c(ic  ,jc-2,kc+2),
     *                                             u_c(ic-1,jc-2,kc+2),
     *                                 1.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     else if (ic.gt.2) then
                        field_j_plus    = interp_cubicb(
     *                                             u_c(ic+1,jc+1,kc+2),
     *                                             u_c(ic,  jc+1,kc+2),
     *                                             u_c(ic-1,jc+1,kc+2),
     *                                             u_c(ic-2,jc+1,kc+2),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j         = interp_cubicb(
     *                                             u_c(ic+1,jc,kc+2),
     *                                             u_c(ic,  jc,kc+2),
     *                                             u_c(ic-1,jc,kc+2),
     *                                             u_c(ic-2,jc,kc+2),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j_minus   = interp_cubicb(
     *                                             u_c(ic+1,jc-1,kc+2),
     *                                             u_c(ic,  jc-1,kc+2),
     *                                             u_c(ic-1,jc-1,kc+2),
     *                                             u_c(ic-2,jc-1,kc+2),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                        field_j_minus2  = interp_cubicb(
     *                                             u_c(ic+1,jc-2,kc+2),
     *                                             u_c(ic,  jc-2,kc+2),
     *                                             u_c(ic-1,jc-2,kc+2),
     *                                             u_c(ic-2,jc-2,kc+2),
     *                                 2.d0+0.5d0*(if-(ic-1-bi)*r-mini),
     *                                             3.d0,2.d0,1.d0,0.d0)
                     else
                         write(*,8)'intrp_parent: S',if,jf,kf+i,ic,jc,kc
                     end if
                     if(ltrace2)write(*,*)'interp:B interp in Y',
     .                   'to get field_k_plus2'
                     field_k_plus2   = interp_cubicb(
     *                                             field_j_plus,
     *                                             field_j,
     *                                             field_j_minus,
     *                                             field_j_minus2,
     *                                 2.d0+0.5d0*(jf-(jc-1-bj)*r-minj),
     *                                             3.d0,2.d0,1.d0,0.d0)
                  else
                     write(*,8)'intrp_parent: G',if,jf+i,kf,ic,jc,kc
                  end if
                  !
                  ! In case one of these points is masked, 
                  ! communicate it to interpolation routine:
                  !
                  if (NINT(mask(if,jf,kf+r)).ge. 0) then
                     field_k_plus = u_f(if,jf,kf+r)
                  else
                     if(ltrace2)write(*,*)'interp:no field_k_plus value'
                     field_k_plus = LARGENUMBER
                  end if
                  if (NINT(mask(if,jf,kf  )).ge. 0) then
                     field_k      = u_f(if,jf,kf  )
                  else
                     if(ltrace2)write(*,*)'interp:no field_k      value'
                     field_k      = LARGENUMBER
                  end if
                  !
                  ! Do the interpolation:
                  !
                  u_f(if,jf,kf+i) = interp_cubicb(field_k_plus2,
     *                                           field_k_plus,
     *                                           field_k,
     *                                           field_k_minus,
     *                                           1.d0+temp,
     *                                           3.d0,2.d0,1.d0,0.d0)
                  else
                     ! Point is masked off:
                     if ( u_f(if,jf,kf+i) .ge. 1.d80) then
                       if (ltrace) then
                       write(*,*)'interp_from_parentb:Z Point is masked'
                       write(*,*)'but has a bad value: ',u_f(if,jf,kf+i)
                       write(*,*)'Setting to zero:'
                       end if
                       u_f(if,jf,kf+i) = 0.d0
                     end if
                  end if
                  if ( u_f(if,jf,kf+i) .ge. 1.d80 .and. ltrace) then
                        write(*,*) 'interp_from_parentb: Z: Bad u_f',
     *                   u_f(if,jf,kf+i),if,jf,kf,ic,jc,kc
                       write(*,*)'interp_from_parentb:Z:field_k_plus2:',
     *                   field_k_plus2
                       write(*,*)'interp_from_parentb:Z:field_k_plus: ',
     *                   field_k_plus
                       write(*,*)'interp_from_parentb:Z:field_k:      ',
     *                   field_k
                       write(*,*)'interp_from_parentb:Z:field_k_minus:',
     *                   field_k_minus
                  end if
               end do
            end do
         end do
      end do
      if (ltrace) call field_dump_info(u_f,nxf,nyf,nzf)
      if (ltrace) write(*,*) 'interp_from_parent:  Done w/ z'

      return
      end    ! END: interp_from_parentb

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interp_from_parentB_weno:                                                 cc
cc                                                                            cc
cc           1) copy data points to corresponding points in grid function     cc
cc           2) interpolate in between these points                           cc
cc                                                                            cc
cc    NB: This version receives a black hole mask and acts accordingly        cc
cc        near the excision surface.                                          cc
cc    NB: bi/j/k: denote the number of points included in u_c "before"        cc
cc                the fine grid data starts in each direction, because        cc
cc                for higher order interpolation, one sends extra points.     cc
cc                Generally, these will be either 0 or 1.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp_from_parentB_weno( u_c, u_f, chr, mask,
     *                         nxc, nyc, nzc,
     *                         mini,minj,mink,
     *                         maxi,maxj,maxk,
     *                         bi,  bj,  bk,
     *                         nxf, nyf, nzf, r, gw)
      implicit      none
      integer       nxc,nyc,nzc, nxf,nyf,nzf, r, gw
      integer       mini,minj,mink, bi,bj,bk
      integer       maxi,maxj,maxk
      real(kind=8)  u_c(nxc,nyc,nzc), u_f(nxf,nyf,nzf),
     *              chr(nxf,nyf,nzf), mask(nxf,nyf,nzf)
      include      'chr.inc'
      integer       ic,jc,kc, if,jf,kf, i
      real(kind=8)  temp
      real(kind=8)  field_j_plus3
      real(kind=8)  field_i_plus2, field_i_plus, field_i,
     *              field_i_minus, field_i_minus2
      real(kind=8)  field_j_plus2, field_j_plus, field_j,
     *              field_j_minus, field_j_minus2
      real(kind=8)  field_k_plus2, field_k_plus, field_k,
     *              field_k_minus, field_k_minus2
      real(kind=8) u(5)
      real(kind=8) x(5)
      real(kind=8) xi
      real(kind=8)   weno_interp_1d
      external       weno_interp_1d
      logical        aligned_i, aligned_j, aligned_k, unmasked
      !
      ! Number which represents data that has been masked off
      ! (from the other grid, in this case, a parent grid):
      !
      include       'largesmall.inc'
      !real(kind=8)   LARGENUMBER
      !parameter    ( LARGENUMBER = 9.d98 )
      !
      logical        ltrace
      parameter (    ltrace = .false. )

 8    format(A,3I3,' | ',3I3)

      if (ltrace) then
         write(*,*) 'interp_from_parent: '
         write(*,*) 'interp_from_parent: bi/j/k:  ',bi,bj,bk
         write(*,*) 'interp_from_parent: nx/y/zf: ',nxf,nyf,nzf
         write(*,*) 'interp_from_parent: nx/y/zc: ',nxc,nyc,nzc
         write(*,*) 'interp_from_parent: mini/j/k:',mini,minj,mink
         write(*,*) 'interp_from_parent: maxi/j/k:',maxi,maxj,maxk
      end if

      write(*,*) '*************WARNING**********************'
      write(*,*) 'Routine interp_from_parentB_weno is not finished'
      write(*,*) 'This routine should not be called yet.'
      write(*,*) '*************WARNING**********************'

      !
      ! Do we have sufficient data to interpolate only symmetrically?
      !
      ic = (mini - mini)/r + 1 + bi
      if ( ic .le. 1 ) then
         write(*,*) 'interp_from_parent: Not sufficient for small x'
      end if
      ic = (maxi - mini)/r + 1 + bi
      if ( ic .ge. nxc ) then
         write(*,*) 'interp_from_parent: Not sufficient for large x',ic
         write(*,*) 'interp_from_parent: bi/j/k:  ',bi,bj,bk
         write(*,*) 'interp_from_parent: nx/y/zf: ',nxf,nyf,nzf
         write(*,*) 'interp_from_parent: nx/y/zc: ',nxc,nyc,nzc
         write(*,*) 'interp_from_parent: mini/j/k:',mini,minj,mink
         write(*,*) 'interp_from_parent: maxi/j/k:',maxi,maxj,maxk
      end if
      jc = (minj - minj)/r + 1 + bj
      if ( jc .le. 1 ) then
         write(*,*) 'interp_from_parent: Not sufficient for small y'
      end if
      jc = (maxj - minj)/r + 1 + bj
      if ( jc .ge. nyc ) then
         write(*,*) 'interp_from_parent: Not sufficient for large y',jc
      end if
      kc = (mink - mink)/r + 1 + bk
      if ( kc .le. 1 ) then
         write(*,*) 'interp_from_parent: Not sufficient for small z'
      end if
      kc = (maxk - mink)/r + 1 + bk
      if ( kc .ge. nzc ) then
         write(*,*) 'interp_from_parent: Not sufficient for large z',kc
      end if

      !call field_dump_info(mask,nxf,nyf,nzf)
      !call field_dump_info(u_c,nxc,nyc,nzc)
c     call field_out3d(mask,0.d0, 'mask',0.d0,1.d0,0.d0,1.d0, 
c    *                       0.d0,1.d0,nxf,nyf,nzf,0)

c
c    Coordinates for weno interp
      x(1) = 0.0
      x(2) = 1.0
      x(3) = 2.0
      x(4) = 3.0
      x(5) = 4.0

      do kf = mink, maxk, r
         kc = (kf - mink)/r + 1 + bk
         do jf = minj, maxj, r
            jc = (jf - minj)/r + 1 + bj
            do if = mini, maxi-r, r
               ic = (if - mini)/r + 1 + bi
               !
               ! copy points in common to both grids
               !
               unmasked =  mask(if,jf,kf) .ge. 0
               if (unmasked .and. u_c(ic,jc,kc).lt.LARGE_M_SMALL) then
               !if (unmasked .and. u_c(ic,jc,kc).lt.LARGENUMBER) then
                  u_f(if,jf,kf) = u_c(ic,jc,kc)
               end if
               !
               !.....interpolate across x-direction
               !
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  unmasked =  mask(if+i,jf,kf) .ge. 0
                  !if (.not. unmasked) write(*,*) 'masked at:',if+i,jf,kf
                  if (ic .gt. 1 .and. ic+3 .le. nxc .and. unmasked )then
                     u(1) = u_c(ic-1,jc,kc)
                     u(2) = u_c(ic,jc,kc)
                     u(3) = u_c(ic+1,jc,kc)
                     u(4) = u_c(ic+2,jc,kc)
                     u(5) = u_c(ic+3,jc,kc)
                     u_f(if+i,jf,kf) =weno_interp_1d(1.d0+temp,u,x)
                  else if (ic .gt. 2 .and. ic+2 .le. nxc
     *                      .and. unmasked ) then
                     u(1) = u_c(ic-2,jc,kc)
                     u(2) = u_c(ic-1,jc,kc)
                     u(3) = u_c(ic,jc,kc)
                     u(4) = u_c(ic+1,jc,kc)
                     u(5) = u_c(ic+2,jc,kc)
                     u_f(if+i,jf,kf) =weno_interp_1d(2.d0+temp,u,x)
                     write(*,*) 'non-symmetric at large x'
                  else
                   !write(*,8)'interp_from_parent: a',if+i,jf,kf,ic,jc,kc
                  end if
               end do
            end do
            if = maxi
            unmasked =  mask(if,jf,kf) .ge. 0
            if (unmasked) then
               ic = (if - mini)/r + 1 + bi
               u_f(if,jf,kf) = u_c(ic,jc,kc)
            end if
         end do
      end do

      !
      !.....interpolate in y-direction
      !
      do kf = mink, maxk, r
         kc = (kf - mink)/r + 1 + bk
         do jf = minj, maxj-r, r
            jc = (jf - minj)/r + 1 + bj
            do if = mini, maxi
               ic = (if - mini)/r + 1 + bi
               aligned_i =  (ic-1-bi)*r+mini .eq. if
               do i = 1, r-1
                 temp = (1.d0*i) / r
                 unmasked =  mask(if,jf+i,kf) .ge. 0
                 if (unmasked) then
                 !
                 ! Establish the j_minus value:
                 !
                 if (jf .gt. minj .and. mask(if,jf-r,kf).ge.0) then
                 !if (jf .gt. minj) then
                    field_j_minus   = u_f(if,jf-r,kf)
                 else if (aligned_i .and. jc.gt.1) then
                    field_j_minus   = u_c(ic,jc-1,  kc)
                 else if (ic+3.le.nxc .and. ic .gt. 1 .and. 
     *                jc.gt.1) then
                    u(1) = u_c(ic-1,jc-1,kc)
                    u(2) = u_c(ic,jc-1,kc)
                    u(3) = u_c(ic+1,jc-1,kc)
                    u(4) = u_c(ic+2,jc-1,kc)
                    u(5) = u_c(ic+3,jc-1,kc)
                    field_j_minus   = weno_interp_1d(
     *                    1.d0+0.5d0*(if-(ic-1-bi)*r-mini),u,x)
                 else if (ic.gt.2 .and. ic+2 .le. nxc .and. 
     *                jc.gt.1) then
                    u(1) = u_c(ic-2,jc-1,kc)
                    u(2) = u_c(ic-1,jc-1,kc)
                    u(3) = u_c(ic,jc-1,kc)
                    u(4) = u_c(ic+1,jc-1,kc)
                    u(5) = u_c(ic+2,jc-1,kc)
                    field_j_minus   = weno_interp_1d(
     *                   2.d0+0.5d0*(if-(ic-1-bi)*r-mini),u,x)
                 else
                     write(*,8)'intrp_parent: e',if,jf+i,kf,ic,jc,kc
                 end if
                 !
                 ! Establish the j_plus2 value:
                 !
                 if (jf+2*r .le. maxj .and.
     *                      mask(if,jf+2*r,kf).ge.0) then
                    field_j_plus2   = u_f(if,jf+2*r,kf)
                 else if (aligned_i .and. jc+2.le.nyc) then
                    field_j_plus2   = u_c(ic,jc+2,  kc)
                 else if (ic+3.le.nxc .and.
     *                    ic .gt. 1 .and. jc+2.le.nyc) then
                    u(1) = u_c(ic-1,jc+2,kc)
                    u(2) = u_c(ic,jc+2,kc)
                    u(3) = u_c(ic+1,jc+2,kc)
                    u(4) = u_c(ic+2,jc+2,kc)
                    u(5) = u_c(ic+3,jc+2,kc)
                    field_j_plus2   = weno_interp_1d(
     *                  1.d0+0.5d0*(if-(ic-1-bi)*r-mini),u,x)
                 else if (ic.gt.2 .and. ic+2.le.nxc 
     *                .and. jc+2.le.nyc) then
                    u(1) = u_c(ic-2,jc+2,kc)
                    u(2) = u_c(ic-1,jc+2,kc)
                    u(3) = u_c(ic,jc+2,kc)
                    u(4) = u_c(ic+1,jc+2,kc)
                    u(5) = u_c(ic+2,jc+2,kc)
                    field_j_plus2   = weno_interp_1d(
     *                  2.d0+0.5d0*(if-(ic-1-bi)*r-mini),u,x)
                 else
                  write(*,8)'intrp_parent: d',if,jf+i,kf,ic,jc,kc
                 end if
                 !
                 ! Establish the j_plus3 value:
                 !
                 if (jf+3*r .le. maxj .and.
     *                                    mask(if,jf+3*r,kf).ge.0) then
                    field_j_plus3   = u_f(if,jf+3*r,kf)
                 else if (aligned_i .and. jc+3.le.nyc) then
                    field_j_plus3   = u_c(ic,jc+3,  kc)
                 else if (ic+3.le.nxc .and.
     *                    ic .gt. 1 .and. jc+3.le.nyc) then
                    u(1) = u_c(ic-1,jc+3,kc)
                    u(2) = u_c(ic,jc+3,kc)
                    u(3) = u_c(ic+1,jc+3,kc)
                    u(4) = u_c(ic+2,jc+3,kc)
                    u(5) = u_c(ic+3,jc+3,kc)
                    field_j_plus3   = weno_interp_1d(
     *                  1.d0+0.5d0*(if-(ic-1-bi)*r-mini),u,x)
                 else if (ic.gt.2 .and. ic+2.le.nxc .and. 
     *                               jc+3.le.nyc) then
                    u(1) = u_c(ic-2,jc+3,kc)
                    u(2) = u_c(ic-1,jc+3,kc)
                    u(3) = u_c(ic,jc+3,kc)
                    u(4) = u_c(ic+1,jc+3,kc)
                    u(5) = u_c(ic+2,jc+3,kc)
                    field_j_plus3   = weno_interp_1d(
     *                  2.d0+0.5d0*(if-(ic-1-bi)*r-mini),u,x)
                 else
                  write(*,8)'intrp_parent: d',if,jf+i,kf,ic,jc,kc
                 end if
                 !
                 ! In case one of these points is masked, 
                 ! communicate it to interpolation routine:
                 !
                 if (mask(if,jf+r,kf) .ge. 0) then
                    field_j_plus = u_f(if,jf+r,kf)
                 else
                    field_j_plus = LARGENUMBER
                 end if
                 if (mask(if,jf,  kf) .ge. 0) then
                    field_j      = u_f(if,jf,  kf)
                 else
                    field_j      = LARGENUMBER
                 end if
                 !
                 ! Do the interpolation:
                 !
                 u(1) = field_j_minus
                 u(2) = field_j
                 u(3) = field_j_plus
                 u(4) = field_j_plus2
                 u(5) = field_j_plus3
                 u_f(if,jf+i,kf) = weno_interp_1d(1.d0+temp,u,x)
                 end if
               end do
            end do
         end do
      end do

      return
      end    ! END: interp_from_parentb
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  weno_interp_1d
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     real(kind=8) function weno_interp_1d(xi,u,x)
!     implicit      none

!     real(kind=8) xi
!     real(kind=8) x(5)
!     real(kind=8) u(5)

!     real(kind=8) epsilon 
!     real(kind=8) c5(5),c(3,3),d(3),beta(3)
!     real(kind=8) omega_t(3),omega(3)
!     real(kind=8) xs(3,3),us(3)
!     real(kind=8) so,ux
!     integer m,i,j

!     epsilon = 1.e-6

!     Coefficients for the end points for the 5-point Lagrange polynomial
!     c5(1) = (xi - x(2))*(xi - x(3))*(xi-x(4))*(xi-x(5))/ 
!    *        ((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))*(x(1)-x(5)))
!     c5(5) = (xi - x(1))*(xi - x(2))*(xi-x(3))*(xi-x(4))/ 
!    *        ((x(5)-x(1))*(x(5)-x(2))*(x(5)-x(3))*(x(5)-x(4)))
! 
!     Coefficients for the substencils 
!     do i=1,3
!       xs(1,i) = x(i+2)
!       xs(2,i) = x(i+1)
!       xs(3,i) = x(i)
!     end do
!     do i=1,3
!       c(i,1) = (xi-xs(i,2))*(xi-xs(i,3))
!    *        /((xs(i,1)-xs(i,2))*(xs(i,1)-xs(i,3)))
!       c(i,2) = (xi-xs(i,1))*(xi-xs(i,3))
!    *        /((xs(i,2)-xs(i,1))*(xs(i,2)-xs(i,3)))
!       c(i,3) = (xi-xs(i,1))*(xi-xs(i,2))
!    *        /((xs(i,3)-xs(i,1))*(xs(i,3)-xs(i,2)))
!     end do
!
!      The 5-point Lagrange stencil is expanded in terms of the substencils
!      The expansion coefficients are
!     d(1) = c5(5) / c(1,3)
!     d(3) = c5(1) / c(3,1)
!     d(2) = 1.d0 - d(1) - d(3) 
!
!     Smoothness indicators
!     beta(1) = (10.0d0*u(3)*u(3) - 31.0d0*u(3)*u(4) 
!    *              + 25.0d0*u(4)*u(4) + 11.0d0*u(3)*u(5)
!    *              - 19.0d0*u(4)*u(5) + 4.0d0*u(5)*u(5)) / 3.0d0;
!     beta(2) = (4.0d0*u(2)*u(2) - 13.0d0*u(2)*u(3) 
!    *              + 13.0d0*u(3)*u(3) + 5.0d0*u(2)*u(4)
!    *              -13.0d0*u(3)*u(4) + 4.0d0*u(4)*u(4)) / 3.0d0;
!     beta(3) = (4.0d0*u(1)*u(1) - 19.0d0*u(1)*u(2) 
!    *                  + 25.0d0*u(2)*u(2) + 11.0d0*u(1)*u(3)
!    *                  - 31.0d0*u(2)*u(3) + 10.0d0*u(3)*u(3)) / 3.0d0

!     so = 0.d0
!     do i=1,3
!       omega_t(i) = d(i)/(epsilon + beta(i))**2
!       so = so + omega_t(i)
!     end do

!     do i=1,3
!       omega(i) = omega_t(i)/so
!     end do

!     ux = 0.d0
!     do i=1,3
!       us(i) = 0.d0
!       do j=1,3
!         us(i) = us(i) + c(i,j)*u((3-i)+j)
!       end do
!       ux = ux + omega(i)*us(i)
!     end do

!     weno_interp_1d = ux

!     return
!     end function weno_interp_1d  ! END: weno_interp_1d


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  Create polynomials for weno stencils that extend out to positions that
cc    will allow interorlation up to the boundaries
cc
cc    the polynomials will be the interpolated values for i+1/2
cc    p(1) will be associated with the stencil that has the leftmost point
cc    at i-3.  That will make p(6) be the stencil that starts at i+2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine weno_polynomials(a,b,c,d,e,f,g,h,p)
!       implicit none
!       real*8, intent(IN) :: a,b,c,d,e,f,g,h
!       real*8, intent(OUT) :: p(6)

!       p(1) = (47.0d0*a-67.d0*b+26.d0*c)/6.d0
!       p(2) = (26.0d0*b-31.d0*c+11.d0*d)/6.d0
!       p(3) = (11.0d0*c- 7.d0*d+ 2.d0*e)/6.d0
!       p(4) = ( 2.d0*d+ 5.d0*e- 1.d0*f)/6.d0
!       p(5) = (-1.d0*e+ 5.d0*f+ 2.d0*g)/6.d0
!       p(6) = ( 2.d0*f- 7.d0*g+11.d0*h)/6.d0

!     end subroutine weno_polynomials

!c The smoothness functions associated with the above polynomials
!     subroutine weno_smoothness(a,b,c,d,e,f,g,h,s)
!       implicit none
!       real *8, intent(IN) :: a,b,c,d,e,f,g,h
!       real *8, intent(OUT):: s(6)

!       s(1)=(c-2.*b+a)**2+(a*a+2.*b*b+c*c-2.*b*(a+c))/2.
!       s(2)=(d-2.*c+b)**2+(b*b+2.*c*c+d*d-2.*c*(b+d))/2.
!       s(3)=(e-2.*d+c)**2+(c*c+2.*d*d+e*e-2.*d*(c+e))/2.
!       s(4)=(f-2.*e+d)**2+(d*d+2.*e*e+f*f-2.*e*(d+f))/2.
!       s(5)=(g-2.*f+e)**2+(e*e+2.*f*f+g*g-2.*f*(e+g))/2.
!       s(6)=(h-2.*g+f)**2+(f*f+2.*g*g+h*h-2.*g*(f+h))/2.
!     end subroutine weno_smoothness

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  paste_mask:     Copy from parent mask info from storage into child        cc 
cc                  and then fill in mask (as opposed to interpolate.         cc 
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine paste_mask( u_c, u_f, nxc,nyc,nzc,
     *                         mini,minj,mink,
     *                         maxi,maxj,maxk,
     *                         bi,bj,bk,
     *                         nxf, nyf, nzf, r)
      implicit    none
      integer     nxc,nyc,nzc, nxf,nyf,nzf, r
      integer     bi,bj,bk
      integer     mini,minj,mink
      integer     maxi,maxj,maxk
      real(kind=8)      u_c(nxc,nyc,nzc), u_f(nxf,nyf,nzf)
      include     'mask.inc'
      integer     ic,jc,kc, if,jf,kf, i
      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'paste_mask: nx/y/zc = ',nxc,nyc,nzc
         write(*,*) 'paste_mask: nx/y/zf = ',nxf,nyf,nzf
         write(*,*) 'paste_mask: bi/j/k  = ',bi,bj,bk
         write(*,*) 'paste_mask: mini/j/k= ',mini,minj,mink
         write(*,*) 'paste_mask: maxi/j/k= ',maxi,maxj,maxk
         write(*,*) 'paste_mask: r =       ', r
      end if

c     .....copy points in common to both grids
      do kf = mink, maxk, r
         kc = (kf - mink)/r + 1 + bk
         do jf = minj, maxj, r
            jc = (jf - minj)/r + 1 + bj
            do if = mini, maxi-r, r
               ic = (if - mini)/r + 1 + bi
               u_f(if,jf,kf) = u_c(ic,jc,kc)
               !.....copy across x-direction (if 
               do i = 1, r-1
                  if (u_c(ic,  jc,kc) .lt. 0.d0 .or.
     *                u_c(ic+1,jc,kc) .lt. 0.d0 ) then
                     u_f(if+i,jf,kf) =   1.d0*BHMASK_MASKED
                  else
                     u_f(if+i,jf,kf) =   1.d0*BHMASK_UNMASKED
                  end if
               end do
            end do
            if = maxi
            ic = (if - mini)/r + 1
            u_f(if,jf,kf) = u_c(ic,jc,kc)
         end do
      end do

      if (ltrace) write(*,*) 'paste_mask: go in y direction'

c     .....interpolate in y-direction
      do kf = mink, maxk, r
         do jf = minj, maxj-r, r
            do if = mini, maxi
               do i = 1, r-1
c                 if (NINT(u_f(if,jf,  kf)) .eq. BHMASK_MASKED .or.
c    *                NINT(u_f(if,jf+r,kf)) .eq. BHMASK_MASKED ) then
                  if (u_f(if,jf,  kf) .lt. 0.d0 .or.
     *                u_f(if,jf+r,kf) .lt. 0.d0 ) then
                     u_f(if,jf+i,kf) =   1.d0*BHMASK_MASKED
                  else
                     u_f(if,jf+i,kf) =   1.d0*BHMASK_UNMASKED
                  end if
               end do
            end do
         end do
      end do

      if (ltrace) write(*,*) 'paste_mask: go in z direction'

c     .....interpolate in z-direction
      do kf = mink, maxk-r, r
         do jf = minj, maxj
            do if = mini, maxi
               do i = 1, r-1
                  if (u_f(if,jf,kf  ) .lt. 0.d0 .or.
     *                u_f(if,jf,kf+r) .lt. 0.d0 ) then
c                 if (NINT(u_f(if,jf,kf  )) .eq. BHMASK_MASKED .or.
c    *                NINT(u_f(if,jf,kf+r)) .eq. BHMASK_MASKED ) then
                     u_f(if,jf,kf+i) =   1.d0*BHMASK_MASKED
                  else
                     u_f(if,jf,kf+i) =   1.d0*BHMASK_UNMASKED
                  end if
               end do
            end do
         end do
      end do

      if (ltrace) write(*,*) 'paste_mask: Done.'

      return
      end    ! END: paste_mask


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  init_fieldGWa:                                                            cc
cc                                                                            cc
cc           NB: newest version for non-unique parents.                       cc
cc                                                                            cc
cc           NB: ghostwidth aware. Uses cubic interpolation within            cc
cc           the boundary of the ghostwidth points, otherwise linear.         cc
cc                                                                            cc
cc           Copies from storage into the field() and then interpolates.      cc
cc           Used to recv part of a field over a network via MPI.             cc
cc           Similar to lin_interpolated3d().                                 cc
cc                                                                            cc
cc           1) copy data points to corresponding points in grid function     cc
cc           2) interpolate in between these points                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_fieldGWa( u_c, u_f, nxc,nyc,nzc,
     *                         mini,minj,mink,
     *                         maxi,maxj,maxk,
     *                         nxf, nyf, nzf, r, gw)
      implicit    none
      integer     nxc,nyc,nzc, nxf,nyf,nzf, r, gw
      integer     mini,minj,mink
      integer     maxi,maxj,maxk
      real(kind=8)      u_c(nxc,nyc,nzc), u_f(nxf,nyf,nzf)
      integer     ic,jc,kc, if,jf,kf, i
      real(kind=8)      temp
      real(kind=8)      interp_cubic
      external          interp_cubic
      !
      ! Use linear interpolation everywhere:
      !    (must do this, if you want cubic use another routine)
      !
      logical     forcelinear
      parameter ( forcelinear = .true. )
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'init_fieldGWa:     nx/y/zc: ',nxc,nyc,nzc
         write(*,*) 'init_fieldGWa:     nx/y/zc: ',nxf,nyf,nzf
         write(*,*) 'init_fieldGWa:    mini/j/k: ',mini,minj,mink
         write(*,*) 'init_fieldGWa:    mini/j/k: ',maxi,maxj,maxk
         write(*,*) 'init_fieldGWa:        r/gw: ',r,gw
         write(*,*) 'init_fieldGWa: forcelinear: ',forcelinear
      end if
c     .....copy points in common to both grids
      do kf = mink, maxk, r
         kc = (kf - mink)/r + 1
         do jf = minj, maxj, r
            jc = (jf - minj)/r + 1
            do if = mini, maxi-r, r
               ic = (if - mini)/r + 1
               if(ltrace)write(*,*)'iGWa: ',if,jf,kf,"<--",ic,jc,kc
               u_f(if,jf,kf) = u_c(ic,jc,kc)
               !.....interpolate across x-direction
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if (ic .lt. 2 .or. ic .gt. nxc-2 .or.forcelinear) then
               if(ltrace)write(*,*)'iGWa: ',if+i,jf,kf,"-|-",ic,  jc,kc,
     *                                                       ic+1,jc,kc,
     *                                                       temp
                     u_f(if+i,jf,kf) =   (1.d0-temp) * u_c(ic  ,jc,kc)
     *                                 +  temp       * u_c(ic+1,jc,kc)
                     if (ltrace.and. .not. forcelinear)
     *                       write(*,*) ' linear:',if+i,jf,kf
                  else
                     u_f(if+i,jf,kf) = interp_cubic(u_c(ic+2,jc,kc),
     *                                              u_c(ic+1,jc,kc),
     *                                              u_c(ic  ,jc,kc),
     *                                              u_c(ic-1,jc,kc),
     *                                              1.d0+temp,
     *                                              3.d0,2.d0,1.d0,0.d0)
                  end if
               end do
            end do
            if = maxi
            ic = (if - mini)/r + 1
            if(ltrace)write(*,*)'iGWa: ',if,jf,kf,"<--",ic,jc,kc
            u_f(if,jf,kf) = u_c(ic,jc,kc)
         end do
      end do

c     .....interpolate in y-direction
      do kf = mink, maxk, r
         do jf = minj, maxj-r, r
            do if = mini, maxi
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if (jf .le. r .or. jf .gt. nyf-2*r.or.forcelinear)then
               if(ltrace)write(*,*)'iGWa: ',if,jf+i,kf,"-|-",jf,jf+r,
     *                                                       temp
                     u_f(if,jf+i,kf) =   (1.d0-temp) * u_f(if,jf  ,kf)
     *                                 +  temp       * u_f(if,jf+r,kf)
                  else
                     u_f(if,jf+i,kf) = interp_cubic(u_f(if,jf+2*r,kf),
     *                                              u_f(if,jf+  r,kf),
     *                                              u_f(if,jf    ,kf),
     *                                              u_f(if,jf-  r,kf),
     *                                              1.d0+temp,
     *                                              3.d0,2.d0,1.d0,0.d0)
                  end if
               end do
            end do
         end do
      end do

c     .....interpolate in z-direction
      do kf = mink, maxk-r, r
         do jf = minj, maxj
            do if = mini, maxi
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if (kf .le. r .or. kf .gt. nzf-2*r.or.forcelinear)then
                     u_f(if,jf,kf+i) =   (1.d0-temp) * u_f(if,jf,kf  )
     *                                 +  temp       * u_f(if,jf,kf+r)
                  else
                     u_f(if,jf,kf+i) = interp_cubic(u_f(if,jf,kf+2*r),
     *                                              u_f(if,jf,kf+  r),
     *                                              u_f(if,jf,kf    ),
     *                                              u_f(if,jf,kf-  r),
     *                                              1.d0+temp,
     *                                              3.d0,2.d0,1.d0,0.d0)
                  end if
               end do
            end do
         end do
      end do

      return
      end    ! END: init_fieldGWa

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  init_fieldGWa_weno:                                                       cc
cc                                                                            cc
cc           NB: newest version for non-unique parents.                       cc
cc                                                                            cc
cc           NB: ghostwidth aware. Uses cubic interpolation within            cc
cc           the boundary of the ghostwidth points, otherwise linear.         cc
cc                                                                            cc
cc           Copies from storage into the field() and then interpolates.      cc
cc           Used to recv part of a field over a network via MPI.             cc
cc           Similar to lin_interpolated3d().                                 cc
cc                                                                            cc
cc           1) copy data points to corresponding points in grid function     cc
cc           2) interpolate in between these points                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_fieldGWa_weno( u_c, u_f, nxc,nyc,nzc,
     *                         mini,minj,mink,
     *                         maxi,maxj,maxk,
     *                         nxf, nyf, nzf, r, gw)
      implicit    none
      integer     nxc,nyc,nzc, nxf,nyf,nzf, r, gw
      integer     mini,minj,mink
      integer     maxi,maxj,maxk
      real(kind=8)      u_c(nxc,nyc,nzc), u_f(nxf,nyf,nzf)
      integer     ic,jc,kc, if,jf,kf, i
      real(kind=8)      temp
      real(kind=8)      u(5)
      real(kind=8)      x(5)
      real(kind=8)      xi
      real(kind=8)      weno_interp_1d
      external          weno_interp_1d
      !
      ! Use linear interpolation everywhere:
      !
      logical     forcelinear
      parameter ( forcelinear = .false. )
      !
      logical     ltrace
      parameter ( ltrace = .false. )

c     Coordinates for weno interp
      x(1) = 0.0
      x(2) = 1.0
      x(3) = 2.0
      x(4) = 3.0
      x(5) = 4.0

c     .....copy points in common to both grids
      do kf = mink, maxk, r
         kc = (kf - mink)/r + 1
         do jf = minj, maxj, r
            jc = (jf - minj)/r + 1
            do if = mini, maxi-r, r
               ic = (if - mini)/r + 1
               u_f(if,jf,kf) = u_c(ic,jc,kc)
               !.....interpolate across x-direction
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if (ic .lt. 3 .or. ic .gt. nxc-2 .or.forcelinear) then
                     u_f(if+i,jf,kf) =   (1.d0-temp) * u_c(ic  ,jc,kc)
     *                                 +  temp       * u_c(ic+1,jc,kc)
                     if (ltrace.and. .not. forcelinear)
     *                       write(*,*) ' linear:',if+i,jf,kf
                  else
                     u(1) = u_c(ic-2,jc,kc)
                     u(2) = u_c(ic-1,jc,kc)
                     u(3) = u_c(ic  ,jc,kc)
                     u(4) = u_c(ic+1,jc,kc)
                     u(5) = u_c(ic+2,jc,kc)
                     u_f(if+i,jf,kf) = weno_interp_1d(2.d0+temp,u,x)
                  end if
               end do
            end do
            if = maxi
            ic = (if - mini)/r + 1
            u_f(if,jf,kf) = u_c(ic,jc,kc)
         end do
      end do

c     .....interpolate in y-direction
      do kf = mink, maxk, r
         do jf = minj, maxj-r, r
            do if = mini, maxi
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if (jf .le. 2*r .or. jf .gt. 
     *                 nyf-2*r.or.forcelinear)then
                     u_f(if,jf+i,kf) =   (1.d0-temp) * u_f(if,jf  ,kf)
     *                                 +  temp       * u_f(if,jf+r,kf)
                  else
                     u(1) = u_f(if,jf-2*r,kf)
                     u(2) = u_f(if,jf-  r,kf)
                     u(3) = u_f(if,jf    ,kf)
                     u(4) = u_f(if,jf+  r,kf)
                     u(5) = u_f(if,jf+2*r,kf)
                     u_f(if,jf+i,kf) = weno_interp_1d(2.d0+temp,u,x)
                  end if
               end do
            end do
         end do
      end do

c     .....interpolate in z-direction
      do kf = mink, maxk-r, r
         do jf = minj, maxj
            do if = mini, maxi
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if (kf .le. 2*r .or. kf .gt. nzf-2*r
     *                 .or.forcelinear)then
                     u_f(if,jf,kf+i) =   (1.d0-temp) * u_f(if,jf,kf  )
     *                                 +  temp       * u_f(if,jf,kf+r)
                  else
                     u(1) = u_f(if,jf,kf-2*r)
                     u(2) = u_f(if,jf,kf-  r)
                     u(3) = u_f(if,jf,kf    )
                     u(4) = u_f(if,jf,kf+  r)
                     u(5) = u_f(if,jf,kf+2*r)
                     u_f(if,jf,kf+i) = weno_interp_1d(2.d0+temp,u,x)
                  end if
               end do
            end do
         end do
      end do

      return
      end    ! END: init_fieldGWa

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  copy_overlapb:                                                            cc
cc                                                                            cc
cc           Copies one contiguous block of a field to storage.               cc
cc           Used to send part of a field over a network via MPI.             cc
cc           Partners with subroutine init_field().                           cc
cc           NB: this version 'masks' out data based on the mask              cc
cc               for use with interp_parentb() and black hole masks.          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copy_overlapb(storage, field,mask,ax,ay,az,nxt,nyt,nzt,
     *                         nx, ny, nz )
      implicit    none
      integer     nxt, nyt,nzt, nx,ny,nz, ax,ay,az
      real(kind=8)      storage(nxt,nyt,nzt), field(nx,ny,nz),
     *                  mask(nx,ny,nz)
      integer     i,j,k, it,jt,kt
      include       'largesmall.inc'
      logical        double_equal
      external       double_equal
      logical        ltrace
      parameter    ( ltrace = .false. )

      if (ltrace) then
         call field_dump_info(mask,nx,ny,nz)
         call field_out3d(mask,0.d0, 'pmask',-3.d0,3.d0,-3.d0,3.d0,
     *                       -3.d0,3.d0,nx,ny,nz,0)
      end if

      do kt = 1, nzt
         k  = kt + az - 1
         do jt = 1, nyt
            j  = jt + ay - 1
            do it = 1, nxt
               i  = it + ax - 1
               if (NINT(mask(i,j,k)) .lt. 0) then
                  ! Value does not exist because it is masked
                  if (ltrace) write(*,*) 'copy_overlapb: masked ',i,j,k
                  storage(it,jt,kt) = LARGENUMBER
               else
                  storage(it,jt,kt) = field(i,j,k)
               end if
            end do
         end do
      end do

      return
      end    ! END: copy_overlapb

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  copy_overlap:                                                             cc
cc                                                                            cc
cc           Copies one contiguous block of a field to storage.               cc
cc           Used to send part of a field over a network via MPI.             cc
cc           Partners with subroutine init_field().                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copy_overlap( storage, field, ax,ay,az,nxt,nyt,nzt,
     *                         nx, ny, nz )
      implicit    none
      integer     nxt, nyt,nzt, nx,ny,nz, ax,ay,az
      real(kind=8)      storage(nxt,nyt,nzt), field(nx,ny,nz)
      integer     i,j,k, it,jt,kt
      logical     ltrace
      parameter ( ltrace = .false. )

      do kt = 1, nzt
         k  = kt + az - 1
         do jt = 1, nyt
            j  = jt + ay - 1
            do it = 1, nxt
               i  = it + ax - 1
               storage(it,jt,kt) = field(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: copy_overlap

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  setmsk_int_bndC:                                                          cc
cc                                                                            cc
cc           Sets mask function for boundary points lying in interior         cc
cc           of other grids.                                                  cc
cc               bwidth determines the width of the region which is           cc
cc                      masked off, but...                                    cc
cc               gwidth determines the width of the region which is           cc
cc                      left alone at the boundaries.                         cc
cc           (1) Set all the intersection points on the boundaries to DECO    cc
cc           (2) Set bwidth-points in from those as DECO boundaries           cc
cc           (3) Except those points which are themselves on the boundaries   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setmsk_int_bndc( chr,
     *                   min_i,  max_i,  min_j,  max_j,  min_k,  max_k,
     *                   min_io, max_io, min_jo, max_jo, min_ko, max_ko,
     *                   nx,  ny,  nz,
     *                   nxo, nyo, nzo, gwidth, bwidth)
      implicit    none
      integer     nx,  ny,  nz,
     *            nxo, nyo, nzo,
     *            min_i,  max_i,  min_j,  max_j,  min_k, max_k,
     *            min_io, max_io, min_jo, max_jo, min_ko,max_ko,
     *            gwidth, bwidth
      real(kind=8)      chr(nx,ny,nz)
      include    'chr.inc'
      !
      ! Indices for the various loops "s" for start and "f" for finish
      !
      integer     i,j,k, is,js,ks, if,jf,kf
      integer     ib, ie, jb, je, kb, ke
      integer     chrpt, AMR_BDY
      logical     both
      logical     ltrace
      parameter ( ltrace = .false. )

      AMR_BDY = NINT(CHR_AMR_BDY)

      ib = min_i
      ie = max_i
      jb = min_j
      je = max_j
      kb = min_k
      ke = max_k
      if (ltrace) then
         write(*,*) 'setmsk_int_bndC: *****'
         write(*,*) 'setmsk_int_bndC: chr_deco_bdy   = ',chr_deco_bdy
         write(*,*) 'setmsk_int_bndC: bwidth         = ',bwidth
         write(*,*) 'setmsk_int_bndC: gwidth         = ',gwidth
         write(*,*) 'setmsk_int_bndC: min_i, max_i   = ',min_i,max_i
         write(*,*) 'setmsk_int_bndC: min_j, max_j   = ',min_j,max_j
         write(*,*) 'setmsk_int_bndC: min_k, max_k   = ',min_k,max_k
         write(*,*) 'setmsk_int_bndC: min_io,max_io  = ',min_io,max_io
         write(*,*) 'setmsk_int_bndC: min_jo,max_jo  = ',min_jo,max_jo
         write(*,*) 'setmsk_int_bndC: min_ko,max_ko  = ',min_ko,max_ko
         write(*,*) 'setmsk_int_bndC: nx,  ny,  nz   = ',nx,ny,nz
         write(*,*) 'setmsk_int_bndC: nxo, nyo, nzo  = ',nxo,nyo,nzo
         write(*,*) 'setmsk_int_bndC: ib/e:          = ',ib, ie
         write(*,*) 'setmsk_int_bndC: jb/e:          = ',jb, je
         write(*,*) 'setmsk_int_bndC: kb/e:          = ',kb, ke
         call field_dump_infob(chr,nx,ny,nz,'setmsk_int_bndC:')
      end if

      both = .false.

      !
      !   Z-face
      !
      if (min_ko.eq. 1 .and. max_ko.eq. nzo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'setmsk_int_bnd:  skipping z-faces'
         goto 10
      else if (min_ko.eq. 1 ) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  upper    z-face'
         ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
         kf = max_k
      else if (max_ko.eq. nzo) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  lower    z-face'
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'setmsk_int_bnd:  both     z-faces'
         both = .true.
         ks   = min_k
         kf   = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      end if
   5  if (ltrace) write(*,*) 'setmsk_int_bnd:  ks     =',ks,kf
      do k = ks, kf
         do j = jb, je
            do i = ib, ie
               chrpt = NINT(chr(i,j,k))
               if ( (k.ne.1).and.(k.ne.nz)
     *              .and. 
     *              chrpt.eq.AMR_BDY .and. 
     *              ( i.eq.1.or.i.eq.nx.or.j.eq.1.or.j.eq.ny) ) then
                  ! Do nothing 
                  ! do not change any AMR boundary points that are not
                  ! on this face
                  if(ltrace)write(*,*)'setmsk_int_bnd: No',i,j,k,chrpt
               else if(((i.eq. 1.and.min_i.eq. 1.and.min_io.eq.  1) .or.
     *                  (i.eq.nx.and.max_i.eq.nx.and.max_io.eq.nxo) .or.
     *                  (j.eq. 1.and.min_j.eq. 1.and.min_jo.eq.  1) .or.
     *                  (j.eq.ny.and.max_j.eq.ny.and.max_jo.eq.nyo))
     *                 .and.chrpt.eq.AMR_BDY)then
                  ! If both grids line up in a coordinate, then
                  ! do not overwrite edge values in case it is a
                  ! real AMR boundary...if not, the intersection 
                  ! with some other grid must overlap it.
                  if(ltrace)write(*,*)'setmsk_int: No edge ',i,j,k,chrpt
               else
                  if(ltrace)write(*,*)'setmsk_int_bnd:Yes',i,j,k,chrpt
                  chr(i,j,k) = 1.d0*CHR_deco_bdy
               end if
            end do
         end do
      end do
      if (both) then
         ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
         kf = max_k
         both = .false.
         goto  5
      end if

      !
      !   Y-face
      !
 10   if (min_jo.eq. 1 .and. max_jo.eq. nyo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'setmsk_int_bnd:  skipping y-faces'
         goto 20
      else if (min_jo.eq. 1 ) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  upper    y-face'
         js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
         jf = max_j
      else if (max_jo.eq. nyo) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  lower    y-face'
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'setmsk_int_bnd:  both     y-faces'
         both = .true.
         js   = min_j
         jf   = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      end if
  15  if (ltrace) write(*,*) 'setmsk_int_bnd:  js     =',js,jf
      do j = js, jf
         do k = kb, ke
            do i = ib, ie
               chrpt = NINT(chr(i,j,k))
               if ( (j.ne.1).and.(j.ne.ny)
     *              .and. 
     *              chrpt.eq.AMR_BDY .and. 
     *              ( i.eq.1.or.i.eq.nx.or.k.eq.1.or.k.eq.nz) ) then
                  ! Do nothing
                  ! do not change any boundary points that are not
                  ! on this face
                  if(ltrace)write(*,*)'setmsk_int_bnd: No',i,j,k,chrpt
               else if(((i.eq. 1.and.min_i.eq. 1.and.min_io.eq.  1) .or.
     *                  (i.eq.nx.and.max_i.eq.nx.and.max_io.eq.nxo) .or.
     *                  (k.eq. 1.and.min_k.eq. 1.and.min_ko.eq.  1) .or.
     *                  (k.eq.nz.and.max_k.eq.nz.and.max_ko.eq.nzo))
     *                 .and.chrpt.eq.AMR_BDY)then
                  ! If both grids line up in a coordinate, then
                  ! do not overwrite edge values in case it is a
                  ! real AMR boundary...if not, the intersection 
                  ! with some other grid must overlap it.
                  if(ltrace)write(*,*)'setmsk_int: No edge ',i,j,k,chrpt
               else
                  if(ltrace)write(*,*)'setmsk_int_bnd:Yes',i,j,k,chrpt
                  chr(i,j,k) = 1.d0*CHR_deco_bdy
               end if
            end do
         end do
      end do
      if (both) then
         js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
         jf = max_j
         both = .false.
         goto 15
      end if

      !
      !   X-face
      !
  20  if (min_io.eq. 1 .and. max_io.eq. nxo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'setmsk_int_bnd:  skipping x-faces'
         goto 30
      else if (min_io.eq. 1 ) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  upper    x-face'
         is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
         if = max_i
      else if (max_io.eq. nxo) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  lower    x-face'
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'setmsk_int_bnd:  both     x-faces'
         both = .true.
         is   = min_i
         if   = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      end if
  25  if (ltrace) write(*,*) 'setmsk_int_bnd:  is     =',is,if
      if(ltrace)write(*,*)'setmsk_int_bnd: jb/e,kb/e:',jb,je,kb,ke
      do i = is, if
         do k = kb, ke
            do j = jb, je
               chrpt = NINT(chr(i,j,k))
               if ( (i.ne.1).and.(i.ne.nx)
     *              .and. 
     *              chrpt.eq.AMR_BDY .and. 
     *              ( j.eq.1.or.j.eq.ny.or.k.eq.1.or.k.eq.nz) ) then
                  ! Do nothing
                  ! do not change any boundary points that are not
                  ! on this face
                  if(ltrace)write(*,*)'setmsk_int_bnd: No',i,j,k,chrpt
               else if(((j.eq. 1.and.min_j.eq. 1.and.min_jo.eq.  1) .or.
     *                  (j.eq.ny.and.max_j.eq.ny.and.max_jo.eq.nyo) .or.
     *                  (k.eq. 1.and.min_k.eq. 1.and.min_ko.eq.  1) .or.
     *                  (k.eq.nz.and.max_k.eq.nz.and.max_ko.eq.nzo))
     *                 .and.chrpt.eq.AMR_BDY)then
                  ! If both grids line up in a coordinate, then
                  ! do not overwrite edge values in case it is a
                  ! real AMR boundary...if not, the intersection 
                  ! with some other grid must overlap it.
                  if(ltrace)write(*,*)'setmsk_int: No edge ',i,j,k,chrpt
               else
                  if(ltrace)write(*,*)'setmsk_int_bnd:Yes',i,j,k,chrpt
                  chr(i,j,k) = 1.d0*CHR_deco_bdy
               end if
            end do
         end do
      end do
      if (both) then
         is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
         if = max_i
         both = .false.
         goto  25
      end if

  30  continue

      if(ltrace)call field_dump_infob(chr,nx,ny,nz,'setmsk_int_bndC:')

      return
      end    ! END: setmsk_int_bndC

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  setmsk_int_bndB:                                                          cc
cc                                                                            cc
cc           Sets mask function for boundary points lying in interior         cc
cc           of other grids.                                                  cc
cc           NB: aware of whether boundaries border true outer boundary       cc
cc               if smallx==0    then lower x boundary is interior            cc
cc                  smallx==1    then lower x boundary is outer boundary      cc
cc           NB: this version is aware of the ghostwidth.                     cc
cc               bwidth determines the width of the region which is           cc
cc                      masked off, but...                                    cc
cc               gwidth determines the width of the region which is           cc
cc                      left alone at the boundaries.                         cc
cc           NB: the min_[i,j,k]  variables denote the intersection region    cc
cc               w/r/t the given grid of dimensions [nx,ny,nz]                cc
cc               the min_[i,j,k]o variables denote the intersection region    cc
cc               w/r/t the other grid of dimensions [nxo,nyo,nzo]             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setmsk_int_bndb( chr,
     *                   min_i,  max_i,  min_j,  max_j,  min_k,  max_k,
     *                   min_io, max_io, min_jo, max_jo, min_ko, max_ko,
     *                   smallx, largex,
     *                   smally, largey,
     *                   smallz, largez,
     *                   nx,  ny,  nz,
     *                   nxo, nyo, nzo, gwidth, bwidth)
      implicit    none
      integer     nx,  ny,  nz,
     *            nxo, nyo, nzo,
     *            min_i,  max_i,  min_j,  max_j,  min_k, max_k,
     *            min_io, max_io, min_jo, max_jo, min_ko,max_ko,
     *                   smallx, largex,
     *                   smally, largey,
     *                   smallz, largez,
     *            gwidth, bwidth
      real(kind=8)      chr(nx,ny,nz)
      include    'chr.inc'
      !
      ! Indices for the various loops "s" for start and "f" for finish
      !
      integer     i,j,k, is,js,ks, if,jf,kf
      integer     ib, ie, jb, je, kb, ke
      !
      ! Boundary width appropriate to each dimension
      !
      integer     bwi, bwj, bwk 
      !
      ! Length of overlap regions:
      !
      integer     li,  lj,  lk
      logical     both
      logical     ltrace
      parameter ( ltrace = .false. )

      !
      ! If intersection starts at boundary of other grid,
      ! then only points beyond the boundary width get copied over
      ! unelss our grid borders the true outer boundary:
      !
      if (min_io .eq. 1 .and. smallx .eq. 0) then
         ib = min_i + bwidth
      else
         ib = min_i
      end if
      if (max_io .eq. nxo .and. largex .eq. 0) then
         ie = max_i - bwidth
      else
         ie = max_i
      end if
      !
      if (min_jo .eq. 1 .and. smally .eq. 0) then
         jb = min_j + bwidth
      else
         jb = min_j
      end if
      if (max_jo .eq. nyo .and. largey .eq. 0) then
         je = max_j - bwidth
      else
         je = max_j
      end if
      !
      if (min_ko .eq. 1 .and. smallz .eq. 0) then
         kb = min_k + bwidth
      else
         kb = min_k
      end if
      if (max_ko .eq. nzo .and. largez .eq. 0) then
         ke = max_k - bwidth
      else
         ke = max_k
      end if

      if (ltrace) then
         write(*,*) 'setmsk_int_bndB: *****'
         !write(*,*) 'setmsk_int_bndB: chr(1,1,1)     = ',chr(1,1,1)
         write(*,*) 'setmsk_int_bndB: chr_deco_bdy   = ',chr_deco_bdy
         write(*,*) 'setmsk_int_bndB: bwidth         = ',bwidth
         write(*,*) 'setmsk_int_bndB: gwidth         = ',gwidth
         write(*,*) 'setmsk_int_bndB: min_i, max_i   = ',min_i,max_i
         write(*,*) 'setmsk_int_bndB: min_j, max_j   = ',min_j,max_j
         write(*,*) 'setmsk_int_bndB: min_k, max_k   = ',min_k,max_k
         write(*,*) 'setmsk_int_bndB: min_io,max_io  = ',min_io,max_io
         write(*,*) 'setmsk_int_bndB: min_jo,max_jo  = ',min_jo,max_jo
         write(*,*) 'setmsk_int_bndB: min_ko,max_ko  = ',min_ko,max_ko
         write(*,*) 'setmsk_int_bndB: nx,  ny,  nz   = ',nx,ny,nz
         write(*,*) 'setmsk_int_bndB: nxo, nyo, nzo  = ',nxo,nyo,nzo
         write(*,*) 'setmsk_int_bndB: small/large x  = ',smallx,largex
         write(*,*) 'setmsk_int_bndB: small/large y  = ',smally,largey
         write(*,*) 'setmsk_int_bndB: small/large z  = ',smallz,largez
         write(*,*) 'setmsk_int_bndB: ib/e:          = ',ib, ie
         write(*,*) 'setmsk_int_bndB: jb/e:          = ',jb, je
         write(*,*) 'setmsk_int_bndB: kb/e:          = ',kb, ke
      end if

      both = .false.


      !
      !   Z-face
      !
      !if (ltrace) call field_dump_info(chr,nx,ny,nz)
      if (min_ko.eq. 1 .and. max_ko.eq. nzo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'setmsk_int_bnd:  skipping z-faces'
         goto 10
      else if (min_ko.eq. 1 ) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  upper    z-face'
         ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
         kf = max_k
      else if (max_ko.eq. nzo) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  lower    z-face'
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'setmsk_int_bnd:  both     z-faces'
         both = .true.
         ks   = min_k
         kf   = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      end if
   5  if (ltrace) write(*,*) 'setmsk_int_bnd:  ks     =',ks,kf
      do k = ks, kf
         do j = jb, je
            do i = ib, ie
               chr(i,j,k) = CHR_deco_bdy
            end do
         end do
      end do
      if (both) then
         ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
         kf = max_k
         both = .false.
         goto  5
      end if

      !
      !   Y-face
      !
 10   if (min_jo.eq. 1 .and. max_jo.eq. nyo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'setmsk_int_bnd:  skipping y-faces'
         goto 20
      else if (min_jo.eq. 1 ) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  upper    y-face'
         js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
         jf = max_j
      else if (max_jo.eq. nyo) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  lower    y-face'
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'setmsk_int_bnd:  both     y-faces'
         both = .true.
         js   = min_j
         jf   = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      end if
  15  if (ltrace) write(*,*) 'setmsk_int_bnd:  js     =',js,jf
      do j = js, jf
         do k = kb, ke
            do i = ib, ie
               chr(i,j,k) = CHR_deco_bdy
            end do
         end do
      end do
      if (both) then
         js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
         jf = max_j
         both = .false.
         goto 15
      end if

      !
      !   X-face
      !
  20  if (min_io.eq. 1 .and. max_io.eq. nxo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'setmsk_int_bnd:  skipping x-faces'
         goto 30
      else if (min_io.eq. 1 ) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  upper    x-face'
         is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
         if = max_i
      else if (max_io.eq. nxo) then
         if (ltrace) write(*,*) 'setmsk_int_bnd:  lower    x-face'
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'setmsk_int_bnd:  both     x-faces'
         both = .true.
         is   = min_i
         if   = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      end if
  25  if (ltrace) write(*,*) 'setmsk_int_bnd:  is     =',is,if
      do i = is, if
         do k = kb, ke
            do j = jb, je
               chr(i,j,k) = CHR_deco_bdy
            end do
         end do
      end do
      if (both) then
         is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
         if = max_i
         both = .false.
         goto  25
      end if

  30  continue


      return
      end    ! END: setmsk_int_bndB

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  copy_int_bndD:                                                            cc
cc                                                                            cc
cc           Adapted from previous version for finite volume fields.          cc
cc           The only meaingful difference in this context is that FV         cc
cc           fields do not extend to the last grid point.                     cc
cc           If fvflag==true, then field is treated as FV.                    cc
cc                                                                            cc
cc           Copies boundary data into a storage vector for use in            cc
cc           sync'ing one grid's boundary data with the interior of           cc
cc           an overlapping grid. Used in conjuction with grid_syncbnd.       cc
cc           Notes:                                                           cc
cc                  (1) Passed to this procedure is the intersection with     cc
cc                      another grid. However, only data values on the        cc
cc                      interior of this grid need be copied.                 cc
cc                  (2) The order in which the data is passed in must,        cc
cc                      of course, match how it is unpacked in paste_int_bnd  cc
cc                  (3) Edges are copied in duplicate to avoid making the     cc
cc                      code too complex.                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copy_int_bndD( storage, field, chr,
     *                   min_i, max_i, min_j, max_j, min_k, max_k,
     *                   min_io,max_io,min_jo,max_jo,min_ko,max_ko,
     *                   lx,rx, ly,ry, lz,rz,
     *                   length, properlength, nx, ny, nz,nxo,nyo,nzo,
     *                   bwidth, fvflag)
      implicit    none
      integer     nx, ny, nz, min_i, max_i, min_j, max_j, min_k, max_k
      integer     nxo,nyo,nzo,min_io,max_io,min_jo,max_jo,min_ko,max_ko
      integer     lx,rx, ly,ry, lz,rz
      integer     length, properlength, bwidth, tmpbwidth
      logical     fvflag
      real(kind=8)      storage(length), field(nx,ny,nz), chr(nx,ny,nz)
      include     'largesmall.inc'
      integer     i,j,k, ii, is, js, ks, ib,jb,kb, ie,je,ke
      integer                if, jf, kf
      logical     both
      ! From ../../include/chr.inc
      integer     DECOBDY
      parameter ( DECOBDY = 8 )
      logical     ltrace
      parameter ( ltrace = .false. )

      !ltrace = .false.
      !if (fvflag) ltrace = .true.
      !
      ! Initialize and then see if we need both sides
      !
      both = .false.

      !
      ! Count/index points stored
      !
      ii   = 0

      ! try *not* adjusting for bwidth:
      if (.true.) then
         tmpbwidth = 0
         !write(*,*) 'copy_int_bnd:using no bwidth'
      else
         tmpbwidth = bwidth
      end if

      !
      ! Paste points but need to avoid simply switching
      ! points from one grid to another ... in other words
      ! points shouldn't be copied and pasted twice:
      !    unless grids share a common boundary, offset
      !    the regions which get pasted.
      !
      if (min_i .eq. 1 .and. min_io .eq. 1 .and. lx .eq. 1) then
         ib = min_i
      else if (min_i .eq. 1) then
         ib = min_i + tmpbwidth
      else
         ib = min_i
      end if
      if (max_i .eq. nx  .and. max_io .eq. nxo .and. rx .eq. 1) then
         ie = max_i
      else if (max_i .eq. nx) then
         ie = max_i - tmpbwidth
         if (ie.lt.ib) ie = ib
      else
         ie = max_i
         if (ie.lt.ib) ib = ie
      end if
      !
      if (min_j .eq. 1 .and. min_jo .eq. 1 .and. ly .eq. 1) then
         jb = min_j
      else if (min_j .eq. 1) then
         jb = min_j + tmpbwidth
      else
         jb = min_j
      end if
      if (max_j .eq. ny  .and. max_jo .eq. nyo .and. ry .eq. 1) then
         je = max_j
      else if (max_j .eq. ny) then
         je = max_j - tmpbwidth
         if (je.lt.jb) je = jb
      else
         je = max_j
         if (je.lt.jb) jb = je
      end if
      !
      if (min_k .eq. 1 .and. min_ko .eq. 1 .and. lz .eq. 1) then
         kb = min_k
      else if (min_k .eq. 1 ) then
         kb = min_k + tmpbwidth
      else
         kb = min_k
      end if
      if (max_k .eq. nz  .and. max_ko .eq. nzo .and. rz .eq. 1) then
         ke = max_k
      else if (max_k .eq. nz) then
         ke = max_k - tmpbwidth
         if (ke.lt.kb) ke = kb
      else
         ke = max_k
         if (ke.lt.kb) kb = ke
      end if

      !
      if (ltrace) then
         write(*,*) 'copy_int_bnd:  nx,ny,nz = ',nx,ny,nz
         write(*,*) 'copy_int_bnd:  min_i,max_i = ',min_i,max_i
         write(*,*) 'copy_int_bnd:  min_j,max_j = ',min_j,max_j
         write(*,*) 'copy_int_bnd:  min_k,max_k = ',min_k,max_k
         write(*,*) 'copy_int_bnd:  ib,ie  = ',ib,ie
         write(*,*) 'copy_int_bnd:  jb,je  = ',jb,je
         write(*,*) 'copy_int_bnd:  kb,ke  = ',kb,ke
         write(*,*) 'copy_int_bnd:  both   = ',both
         write(*,*) 'copy_int_bnd:  bwidth = ',bwidth
         write(*,*) 'copy_int_bnd:  fvflag = ',fvflag
         write(*,*) 'copy_int_bnd: ',lx,rx,ly,ry,lz,rz
      end if

      !
      !   Z-face
      !
      if ( (min_k .eq.1.and.max_k .eq.nz ) .or.
     *     (min_ko.eq.1.and.max_ko.eq.nzo)     ) then
!     if (min_k .eq. 1 .and. max_k .eq. nz) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'copy_int_bnd:  skipping z-faces'
         goto 10
      else if (min_k .eq. 1 ) then
         if (ltrace) write(*,*) 'copy_int_bnd:  upper    z-face'
         if (fvflag) then
            ks = max( max_k - (bwidth-1) -1, (max_k+min_k)/2 )
            kf = max_k - 1
         else
            ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
            kf = max_k
         end if
      else if (max_k .eq. nz) then
         if (ltrace) write(*,*) 'copy_int_bnd:  lower    z-face'
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'copy_int_bnd:  both     z-faces'
         both = .true.
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      end if
   5  continue
      if (ltrace) write(*,*) 'copy_int_bndB: ks =',ks,kf
      do k = ks, kf
         do j = jb, je
            do i = ib, ie
               ii          = ii + 1
               if (nint(chr(i,j,k)).ne.DECOBDY) then
                  storage(ii) = field(i,j,k)
               else
                  storage(ii) = LARGENUMBER
               end if
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            ks = max( max_k - (bwidth-1) -1, (max_k+min_k)/2 )
            kf = max_k - 1
         else
            ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
            kf = max_k
         end if
         both = .false.
         goto  5
      end if
      if (ltrace) write(*,*) 'copy_int_bnd:  z-faces done, ii=',ii

      !
      !   Y-face
      !
  10  if ( (min_j .eq.1.and.max_j .eq.ny ) .or.
     *     (min_jo.eq.1.and.max_jo.eq.nyo)     ) then
! 10  if (min_j .eq. 1 .and. max_j .eq. ny) then
         !  y-face boundaries fall on boundaries
         !  (not interiors) so skip to x-faces
         if (ltrace) write(*,*) 'copy_int_bnd:  skipping y-faces'
         goto 20
      else if (min_j .eq. 1 ) then
         if (ltrace) write(*,*) 'copy_int_bnd:  upper    y-face'
         if (fvflag) then
            js = max( max_j - (bwidth-1) -1, (max_j+min_j)/2 )
            jf = max_j - 1
         else
            js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
            jf = max_j
         end if
      else if (max_j .eq. ny) then
         if (ltrace) write(*,*) 'copy_int_bnd:  lower    y-face'
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'copy_int_bnd:  both     y-faces'
         both = .true.
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      end if
  15  continue
      if (ltrace) write(*,*) 'copy_int_bndB: js =',js,jf
      do j = js, jf
         do k = kb, ke
            do i = ib, ie
               ii          = ii + 1
               if (nint(chr(i,j,k)).ne.DECOBDY) then
                  storage(ii) = field(i,j,k)
               else
                  storage(ii) = LARGENUMBER
               end if
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            js = max( max_j - (bwidth-1) -1, (max_j+min_j)/2 )
            jf = max_j - 1
         else
            js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
            jf = max_j
         end if
         both = .false.
         goto 15
      end if
      if (ltrace) write(*,*) 'copy_int_bnd:  y-faces done, ii=',ii

      !
      !   X-face
      !
  20  if ( (min_i .eq.1.and.max_i .eq.nx ) .or.
     *     (min_io.eq.1.and.max_io.eq.nxo)     ) then
! 20  if (min_i .eq. 1 .and. max_i .eq. nx) then
         !  x-face boundaries fall on boundaries
         !  (not interiors) so done
         if (ltrace) write(*,*) 'copy_int_bnd:  skipping x-faces'
         goto 30
      else if (min_i .eq. 1) then
         if (ltrace) write(*,*) 'copy_int_bnd:  upper    x-face'
         if (fvflag) then
            is = max( max_i - (bwidth-1) -1, (max_i+min_i)/2 )
            if = max_i - 1
         else
            is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
            if = max_i
         end if
      else if (max_i .eq. nx) then
         if (ltrace) write(*,*) 'copy_int_bnd:  lower    x-face'
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'copy_int_bnd:  both     x-faces'
         both = .true.
         is   = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      end if
  25  continue
      if (ltrace) write(*,*) 'copy_int_bndB: is = ',is,if
      do i = is, if
         do k = kb, ke
            do j = jb, je
               ii          = ii + 1
               if (nint(chr(i,j,k)).ne.DECOBDY) then
                  storage(ii) = field(i,j,k)
               else
                  storage(ii) = LARGENUMBER
               end if
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            is = max( max_i - (bwidth-1) -1, (max_i+min_i)/2 )
            if = max_i - 1
         else
            is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
            if = max_i
         end if
         both = .false.
         goto 25
      end if


  30  continue

      if (ltrace) write(*,*) 'copy_int_bnd:          Done, ii=',ii

      if (ii.gt.length) then
         write(*,*) 'Problem in copy_int_bnd',ii,length
         write(*,*) 'Need more storage'
         stop
      end if

      properlength = ii

      return
      end    ! END: copy_int_bndD
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  copy_int_bndC:                                                            cc
cc                                                                            cc
cc           Adapted from previous version for finite volume fields.          cc
cc           The only meaingful difference in this context is that FV         cc
cc           fields do not extend to the last grid point.                     cc
cc           If fvflag==true, then field is treated as FV.                    cc
cc                                                                            cc
cc           Copies boundary data into a storage vector for use in            cc
cc           sync'ing one grid's boundary data with the interior of           cc
cc           an overlapping grid. Used in conjuction with grid_syncbnd.       cc
cc           Notes:                                                           cc
cc                  (1) Passed to this procedure is the intersection with     cc
cc                      another grid. However, only data values on the        cc
cc                      interior of this grid need be copied.                 cc
cc                  (2) The order in which the data is passed in must,        cc
cc                      of course, match how it is unpacked in paste_int_bnd  cc
cc                  (3) Edges are copied in duplicate to avoid making the     cc
cc                      code too complex.                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copy_int_bndC( storage, field,
     *                   min_i, max_i, min_j, max_j, min_k, max_k,
     *                   min_io,max_io,min_jo,max_jo,min_ko,max_ko,
     *                   lx,rx, ly,ry, lz,rz,
     *                   length, properlength, nx, ny, nz,nxo,nyo,nzo,
     *                   bwidth, fvflag)
      implicit    none
      integer     nx, ny, nz, min_i, max_i, min_j, max_j, min_k, max_k
      integer     nxo,nyo,nzo,min_io,max_io,min_jo,max_jo,min_ko,max_ko
      integer     lx,rx, ly,ry, lz,rz
      integer     length, properlength, bwidth
      logical     fvflag
      real(kind=8)      storage(length), field(nx,ny,nz)
      integer     i,j,k, ii, is, js, ks, ib,jb,kb, ie,je,ke
      integer                if, jf, kf
      logical     both
      logical     ltrace
      parameter ( ltrace = .false. )

      !ltrace = .false.
      !if (fvflag) ltrace = .true.
      !
      ! Initialize and then see if we need both sides
      !
      both = .false.

      !
      ! Count/index points stored
      !
      ii   = 0

      !
      ! Paste points but need to avoid simply switching
      ! points from one grid to another ... in other words
      ! points shouldn't be copied and pasted twice:
      !    unless grids share a common boundary, offset
      !    the regions which get pasted.
      !
      if (min_i .eq. 1 .and. min_io .eq. 1 .and. lx .eq. 1) then
         ib = min_i
      else if (min_i .eq. 1) then
         ib = min_i + bwidth
      else
         ib = min_i
      end if
      if (max_i .eq. nx  .and. max_io .eq. nxo .and. rx .eq. 1) then
         ie = max_i
      else if (max_i .eq. nx) then
         ie = max_i - bwidth
         if (ie.lt.ib) ie = ib
      else
         ie = max_i
         if (ie.lt.ib) ib = ie
      end if
      !
      if (min_j .eq. 1 .and. min_jo .eq. 1 .and. ly .eq. 1) then
         jb = min_j
      else if (min_j .eq. 1) then
         jb = min_j + bwidth
      else
         jb = min_j
      end if
      if (max_j .eq. ny  .and. max_jo .eq. nyo .and. ry .eq. 1) then
         je = max_j
      else if (max_j .eq. ny) then
         je = max_j - bwidth
         if (je.lt.jb) je = jb
      else
         je = max_j
         if (je.lt.jb) jb = je
      end if
      !
      if (min_k .eq. 1 .and. min_ko .eq. 1 .and. lz .eq. 1) then
         kb = min_k
      else if (min_k .eq. 1 ) then
         kb = min_k + bwidth
      else
         kb = min_k
      end if
      if (max_k .eq. nz  .and. max_ko .eq. nzo .and. rz .eq. 1) then
         ke = max_k
      else if (max_k .eq. nz) then
         ke = max_k - bwidth
         if (ke.lt.kb) ke = kb
      else
         ke = max_k
         if (ke.lt.kb) kb = ke
      end if

      !
      if (ltrace) then
         write(*,*) 'copy_int_bnd:  nx,ny,nz = ',nx,ny,nz
         write(*,*) 'copy_int_bnd:  min_i,max_i = ',min_i,max_i
         write(*,*) 'copy_int_bnd:  min_j,max_j = ',min_j,max_j
         write(*,*) 'copy_int_bnd:  min_k,max_k = ',min_k,max_k
         write(*,*) 'copy_int_bnd:  ib,ie  = ',ib,ie
         write(*,*) 'copy_int_bnd:  jb,je  = ',jb,je
         write(*,*) 'copy_int_bnd:  kb,ke  = ',kb,ke
         write(*,*) 'copy_int_bnd:  both   = ',both
         write(*,*) 'copy_int_bnd:  bwidth = ',bwidth
         write(*,*) 'copy_int_bnd:  fvflag = ',fvflag
         write(*,*) 'copy_int_bnd: ',lx,rx,ly,ry,lz,rz
      end if

      !
      !   Z-face
      !
      if (min_k .eq. 1 .and. max_k .eq. nz) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'copy_int_bnd:  skipping z-faces'
         goto 10
      else if (min_k .eq. 1 ) then
         if (ltrace) write(*,*) 'copy_int_bnd:  upper    z-face'
         if (fvflag) then
            ks = max( max_k - (bwidth-1) -1, (max_k+min_k)/2 )
            kf = max_k - 1
         else
            ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
            kf = max_k
         end if
      else if (max_k .eq. nz) then
         if (ltrace) write(*,*) 'copy_int_bnd:  lower    z-face'
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'copy_int_bnd:  both     z-faces'
         both = .true.
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      end if
   5  continue
      if (ltrace) write(*,*) 'copy_int_bndB: ks =',ks,kf
      do k = ks, kf
         do j = jb, je
            do i = ib, ie
               ii          = ii + 1
               storage(ii) = field(i,j,k)
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            ks = max( max_k - (bwidth-1) -1, (max_k+min_k)/2 )
            kf = max_k - 1
         else
            ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
            kf = max_k
         end if
         both = .false.
         goto  5
      end if
      if (ltrace) write(*,*) 'copy_int_bnd:  z-faces done, ii=',ii

      !
      !   Y-face
      !
  10  if (min_j .eq. 1 .and. max_j .eq. ny) then
         !  y-face boundaries fall on boundaries
         !  (not interiors) so skip to x-faces
         if (ltrace) write(*,*) 'copy_int_bnd:  skipping y-faces'
         goto 20
      else if (min_j .eq. 1 ) then
         if (ltrace) write(*,*) 'copy_int_bnd:  upper    y-face'
         if (fvflag) then
            js = max( max_j - (bwidth-1) -1, (max_j+min_j)/2 )
            jf = max_j - 1
         else
            js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
            jf = max_j
         end if
      else if (max_j .eq. ny) then
         if (ltrace) write(*,*) 'copy_int_bnd:  lower    y-face'
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'copy_int_bnd:  both     y-faces'
         both = .true.
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      end if
  15  continue
      if (ltrace) write(*,*) 'copy_int_bndB: js =',js,jf
      do j = js, jf
         do k = kb, ke
            do i = ib, ie
               ii          = ii + 1
               storage(ii) = field(i,j,k)
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            js = max( max_j - (bwidth-1) -1, (max_j+min_j)/2 )
            jf = max_j - 1
         else
            js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
            jf = max_j
         end if
         both = .false.
         goto 15
      end if
      if (ltrace) write(*,*) 'copy_int_bnd:  y-faces done, ii=',ii

      !
      !   X-face
      !
  20  if (min_i .eq. 1 .and. max_i .eq. nx) then
         !  x-face boundaries fall on boundaries
         !  (not interiors) so done
         if (ltrace) write(*,*) 'copy_int_bnd:  skipping x-faces'
         goto 30
      else if (min_i .eq. 1) then
         if (ltrace) write(*,*) 'copy_int_bnd:  upper    x-face'
         if (fvflag) then
            is = max( max_i - (bwidth-1) -1, (max_i+min_i)/2 )
            if = max_i - 1
         else
            is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
            if = max_i
         end if
      else if (max_i .eq. nx) then
         if (ltrace) write(*,*) 'copy_int_bnd:  lower    x-face'
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'copy_int_bnd:  both     x-faces'
         both = .true.
         is   = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      end if
  25  continue
      if (ltrace) write(*,*) 'copy_int_bndB: is = ',is,if
      do i = is, if
         do k = kb, ke
            do j = jb, je
               ii          = ii + 1
               storage(ii) = field(i,j,k)
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            is = max( max_i - (bwidth-1) -1, (max_i+min_i)/2 )
            if = max_i - 1
         else
            is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
            if = max_i
         end if
         both = .false.
         goto 25
      end if


  30  continue

      if (ltrace) write(*,*) 'copy_int_bnd:          Done, ii=',ii

      if (ii.gt.length) then
         write(*,*) 'Problem in copy_int_bnd',ii,length
         write(*,*) 'Need more storage'
         stop
      end if

      properlength = ii

      return
      end    ! END: copy_int_bndC


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  paste_int_bndD:                                                           cc
cc                                                                            cc
cc           Adapted from previous version for finite volume fields.          cc
cc           The only meaingful difference in this context is that FV         cc
cc           fields do not extend to the last grid point.                     cc
cc           If fvflag==true, then field is treated as FV.                    cc
cc                                                                            cc
cc           Pastes boundary data into a field   vector for use in            cc
cc           sync'ing one grid's boundary data with the interior of           cc
cc           an overlapping grid. Used in conjuction with grid_syncbnd.       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine paste_int_bndD( storage, field,chr,
     *                   min_i,  max_i,  min_j,  max_j,  min_k,  max_k,
     *                   min_io, max_io, min_jo, max_jo, min_ko, max_ko,
     *                   lx,rx, ly,ry, lz,rz,
     *                   nxo, nyo, nzo,
     *                   length, nx, ny, nz, bwidth, fvflag )
      implicit    none
      integer     nx, ny, nz, min_i, max_i, min_j, max_j,  min_k, max_k
      !
      ! info about the *o*ther grid with which we exchange data:
      !
      integer     nxo,nyo,nzo,min_io,max_io,min_jo,max_jo, min_ko,max_ko
      !
      ! Flags for the boundaries:
      !    lx/y/z --- left  boundary ==0 if will get overwritten by siblings
      !    rx/y/z --- right boundary ==1 if not
      !
      integer     lx,rx, ly,ry, lz,rz
      integer     length, bwidth, tmpbwidth
      logical     fvflag
      real(kind=8)      storage(length), field(nx,ny,nz), chr(nx,ny,nz)
      include     'largesmall.inc'
      integer     i,j,k, ii, is, js, ks, ib,jb,kb,ie,je,ke
      integer                if, jf, kf
      ! Just for debugging
      integer     numptscopied
      logical     both
      ! From ../../include/chr.inc
      integer     DECOBDY
      parameter ( DECOBDY = 8 )
      integer     AMRBDY
      parameter ( AMRBDY = 1 )


      logical     ltrace
      parameter ( ltrace = .false. )

      numptscopied = 0
      both = .false.
      ii   = 0


      ! try *not* adjusting for bwidth:
      if (.true.) then
         tmpbwidth = 0
         !write(*,*) 'paste_int_bnd:using no bwidth'
      else
         tmpbwidth = bwidth
      end if


      !
      ! Paste points but need to avoid simply switching
      ! points from one grid to another ... in other words
      ! points shouldn't be copied and pasted twice:
      !    unless grids share a common boundary, offset
      !    the regions which get pasted.
      !
      if (min_io .eq. 1 .and. min_i .eq. 1 .and. lx .eq. 1) then
         ib = min_i
      else if (min_io .eq. 1) then
         ib = min_i + tmpbwidth
      else
         ib = min_i
      end if
      if (max_io .eq. nxo .and. max_i .eq. nx .and. rx .eq. 1) then
         ie = max_i
      else if (max_io .eq. nxo) then
         ie = max_i - tmpbwidth
         if (ie.lt.ib) ie = ib
      else
         ie = max_i
         if (ie.lt.ib) ib = ie
      end if
      !
      if (min_jo .eq. 1 .and. min_j .eq. 1 .and. ly .eq. 1) then
         jb = min_j
      else if (min_jo .eq. 1) then
         jb = min_j + tmpbwidth
      else
         jb = min_j
      end if
      if (max_jo .eq. nyo .and. max_j .eq. ny .and. ry .eq. 1) then
         je = max_j
      else if (max_jo .eq. nyo ) then
         je = max_j - tmpbwidth
         if (je.lt.jb) je = jb
      else
         je = max_j
         if (je.lt.jb) jb = je
      end if
      !
      if (min_ko .eq. 1 .and. min_k .eq. 1 .and. lz .eq. 1) then
         kb = min_k
      else if (min_ko .eq. 1) then
         kb = min_k + tmpbwidth
      else
         kb = min_k
      end if
      if (max_ko .eq. nzo .and. max_k .eq. nz .and. rz .eq. 1) then
         ke = max_k
      else if (max_ko .eq. nzo) then
         ke = max_k - tmpbwidth
         if (ke.lt.kb) ke = kb
      else
         ke = max_k
         if (ke.lt.kb) kb = ke
      end if



      if (ltrace) then
         write(*,*) 'paste_int_bndD: nx, ny, nz  = ',nx, ny, nz
         write(*,*) 'paste_int_bndD: nxo,nyo,nzo = ',nxo,nyo,nzo
         write(*,*) 'paste_int_bndD: min_i, min_j, min_k  =',
     *                                        min_i, min_j, min_k
         write(*,*) 'paste_int_bndD: max_i, max_j, max_k  =',
     *                                        max_i, max_j, max_k
         write(*,*) 'paste_int_bndD: min_io,min_jo,min_ko =',
     *                                        min_io,min_jo,min_ko
         write(*,*) 'paste_int_bndD: max_io,max_jo,max_ko =',
     *                                        max_io,max_jo,max_ko
         write(*,*) 'paste_int_bndD: ib,ie  = ',ib,ie
         write(*,*) 'paste_int_bndD: jb,je  = ',jb,je
         write(*,*) 'paste_int_bndD kb,ke  = ',kb,ke
         write(*,*) 'paste_int_bndD: ',lx,rx,ly,ry,lz,rz
         write(*,*) 'paste_int_bndD: fvflag = ',fvflag
         numptscopied = 0
      end if

      !
      !   Z-face
      !
      if ( (min_k .eq.1.and.max_k .eq.nz ) .or.
     *     (min_ko.eq.1.and.max_ko.eq.nzo)     ) then
!     if (min_ko .eq. 1 .and. max_ko .eq. nzo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'paste_int_bndD: skipping z-faces'
         goto 10
      else if (min_ko .eq. 1 ) then
         if (ltrace) write(*,*) 'paste_int_bndD: upper    z-face'
         if (fvflag) then
            ks = max( max_k - (bwidth-1)-1, (max_k+min_k)/2 )
            kf = max_k - 1
         else
            ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
            kf = max_k
         end if
      else if (max_ko .eq. nzo) then
         if (ltrace) write(*,*) 'paste_int_bndD: lower    z-face'
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'paste_int_bndD: both     z-faces'
         both = .true.
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      end if
  5   continue
      if (ltrace) write(*,*) 'paste_int_bndD: ks = ',ks,kf
      do k = ks, kf
         do j = jb, je
            do i = ib, ie
               ii           = ii + 1
               if (nint(chr(i,j,k)).eq.DECOBDY .and.
!              if ( (nint(chr(i,j,k)).eq.DECOBDY.or.
!    *               nint(chr(i,j,k)).eq.AMRBDY     ) .and.
     *             storage(ii)     .lt.LARGENUMBER ) then
                  field(i,j,k) = storage(ii)
                  numptscopied = numptscopied + 1
               end if
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            ks = max( max_k - (bwidth-1)-1, (max_k+min_k)/2 )
            kf = max_k - 1
         else
            ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
            kf = max_k
         end if
         both = .false.
         goto  5
      end if

      if (ltrace) write(*,*) 'paste_int_bndD: z-faces done, ii=',ii

      !
      !   Y-face
      !
  10  if ( (min_j .eq.1.and.max_j .eq.ny ) .or.
     *     (min_jo.eq.1.and.max_jo.eq.nyo)     ) then
!10   if (min_jo .eq. 1 .and. max_jo .eq. nyo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'paste_int_bndD: skipping y-faces'
         goto 20
      else if (min_jo .eq. 1 ) then
         if (ltrace) write(*,*) 'paste_int_bndD: upper    y-face'
         if (fvflag) then
            js = max( max_j - (bwidth-1)-1, (max_j+min_j)/2 )
            jf = max_j - 1
         else
            js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
            jf = max_j
         end if
      else if (max_jo .eq. nyo) then
         if (ltrace) write(*,*) 'paste_int_bndD: lower    y-face'
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'paste_int_bndD: both     y-faces'
         both = .true.
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      end if
 15   continue
      if (ltrace) write(*,*) 'paste_int_bndD: js = ',js,jf
      do j = js, jf
         do k = kb, ke
            do i = ib, ie
               ii           = ii + 1
               if (nint(chr(i,j,k)).eq.DECOBDY .and.
!              if ( (nint(chr(i,j,k)).eq.DECOBDY.or.
!    *               nint(chr(i,j,k)).eq.AMRBDY     ) .and.
     *             storage(ii)     .lt.LARGENUMBER ) then
                  field(i,j,k) = storage(ii)
                  numptscopied = numptscopied + 1
               end if
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            js = max( max_j - (bwidth-1)-1, (max_j+min_j)/2 )
            jf = max_j - 1
         else
            js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
            jf = max_j
         end if
         both = .false.
         goto 15
      end if

      if (ltrace) write(*,*) 'paste_int_bndD: y-faces done, ii=',ii

      !
      !   X-faces:
      !
  20  if ( (min_i .eq.1.and.max_i .eq.nx ) .or.
     *     (min_io.eq.1.and.max_io.eq.nxo)     ) then
!20   if (min_io .eq. 1 .and. max_io .eq. nxo) then
         !  x-face boundaries fall on boundaries
         !  (not interiors) so skip 
         if (ltrace) write(*,*) 'paste_int_bndD: skipping x-faces'
         goto 30
      else if (min_io .eq. 1 ) then
         if (ltrace) write(*,*) 'paste_int_bndD: upper    x-face'
         if (fvflag) then
            is = max( max_i - (bwidth-1)-1, (max_i+min_i)/2 )
            if = max_i - 1
         else
            is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
            if = max_i
         end if
      else if (max_io .eq. nxo) then
         if (ltrace) write(*,*) 'paste_int_bndD: lower    x-face'
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'paste_int_bndD: both     x-faces'
         both = .true.
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      end if
  25  continue
      if (ltrace) write(*,*) 'paste_int_bndD: is = ',is,if
      do i = is, if
         do k = kb, ke
            do j = jb, je
               ii           = ii + 1
               if (nint(chr(i,j,k)).eq.DECOBDY .and.
!              if ( (nint(chr(i,j,k)).eq.DECOBDY.or.
!    *               nint(chr(i,j,k)).eq.AMRBDY     ) .and.
     *             storage(ii)     .lt.LARGENUMBER ) then
                  field(i,j,k) = storage(ii)
                  numptscopied = numptscopied + 1
               end if
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            is = max( max_i - (bwidth-1)-1, (max_i+min_i)/2 )
            if = max_i - 1
         else
            is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
            if = max_i
         end if
         both = .false.
         goto 25
      end if

  30  continue

      if (ltrace) write(*,*) 'paste_int_bndD: x-faces done, ii=',ii

      if (ii.ne.length) then
      !if (.false.) then
         write(*,*) 'paste_int_bndD: ii     = ', ii
         write(*,*) 'paste_int_bndD: length = ', length
         write(*,*) 'paste_int_bndD: ib,ie  = ',ib,ie
         write(*,*) 'paste_int_bndD: jb,je  = ',jb,je
         write(*,*) 'paste_int_bndD: kb,ke  = ',kb,ke
         write(*,*) 'Problem in paste_int_bndD'
         write(*,*) 'Problem in paste_int_bndD'
         write(*,*) 'Problem in paste_int_bndD'
         write(*,*) 'Problem in paste_int_bndD'
         write(*,*) 'Problem in paste_int_bndD'
         write(*,*) 'Problem in paste_int_bndD'
         write(*,*) 'Problem in paste_int_bndD'
         write(*,*) 'Problem in paste_int_bndD'
         write(*,*) 'Problem in paste_int_bndD'
         write(*,*) 'Problem in paste_int_bndD'
         !stop
      end if

      if(ltrace)write(*,*)'paste_int_bndD: Done. numptscopied:',
     .                     numptscopied

      return
      end    ! END: paste_int_bndD

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  paste_int_bndC:                                                           cc
cc                                                                            cc
cc           Adapted from previous version for finite volume fields.          cc
cc           The only meaingful difference in this context is that FV         cc
cc           fields do not extend to the last grid point.                     cc
cc           If fvflag==true, then field is treated as FV.                    cc
cc                                                                            cc
cc           Pastes boundary data into a field   vector for use in            cc
cc           sync'ing one grid's boundary data with the interior of           cc
cc           an overlapping grid. Used in conjuction with grid_syncbnd.       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine paste_int_bndC( storage, field,
     *                   min_i,  max_i,  min_j,  max_j,  min_k,  max_k,
     *                   min_io, max_io, min_jo, max_jo, min_ko, max_ko,
     *                   lx,rx, ly,ry, lz,rz,
     *                   nxo, nyo, nzo,
     *                   length, nx, ny, nz, bwidth, fvflag )
      implicit    none
      integer     nx, ny, nz, min_i, max_i, min_j, max_j,  min_k, max_k
      !
      ! info about the *o*ther grid with which we exchange data:
      !
      integer     nxo,nyo,nzo,min_io,max_io,min_jo,max_jo, min_ko,max_ko
      !
      ! Flags for the boundaries:
      !    lx/y/z --- left  boundary ==0 if will get overwritten by siblings
      !    rx/y/z --- right boundary ==1 if not
      !
      integer     lx,rx, ly,ry, lz,rz
      integer     length, bwidth
      logical     fvflag
      real(kind=8)      storage(length), field(nx,ny,nz)
      integer     i,j,k, ii, is, js, ks, ib,jb,kb,ie,je,ke
      integer                if, jf, kf
      logical     both

      logical     ltrace
      parameter ( ltrace = .false. )

      both = .false.
      ii   = 0

      !
      ! Paste points but need to avoid simply switching
      ! points from one grid to another ... in other words
      ! points shouldn't be copied and pasted twice:
      !    unless grids share a common boundary, offset
      !    the regions which get pasted.
      !
      if (min_io .eq. 1 .and. min_i .eq. 1 .and. lx .eq. 1) then
         ib = min_i
      else if (min_io .eq. 1) then
         ib = min_i + bwidth
      else
         ib = min_i
      end if
      if (max_io .eq. nxo .and. max_i .eq. nx .and. rx .eq. 1) then
         ie = max_i
      else if (max_io .eq. nxo) then
         ie = max_i - bwidth
         if (ie.lt.ib) ie = ib
      else
         ie = max_i
         if (ie.lt.ib) ib = ie
      end if
      !
      if (min_jo .eq. 1 .and. min_j .eq. 1 .and. ly .eq. 1) then
         jb = min_j
      else if (min_jo .eq. 1) then
         jb = min_j + bwidth
      else
         jb = min_j
      end if
      if (max_jo .eq. nyo .and. max_j .eq. ny .and. ry .eq. 1) then
         je = max_j
      else if (max_jo .eq. nyo ) then
         je = max_j - bwidth
         if (je.lt.jb) je = jb
      else
         je = max_j
         if (je.lt.jb) jb = je
      end if
      !
      if (min_ko .eq. 1 .and. min_k .eq. 1 .and. lz .eq. 1) then
         kb = min_k
      else if (min_ko .eq. 1) then
         kb = min_k + bwidth
      else
         kb = min_k
      end if
      if (max_ko .eq. nzo .and. max_k .eq. nz .and. rz .eq. 1) then
         ke = max_k
      else if (max_ko .eq. nzo) then
         ke = max_k - bwidth
         if (ke.lt.kb) ke = kb
      else
         ke = max_k
         if (ke.lt.kb) kb = ke
      end if



      if (ltrace) then
         write(*,*) 'paste_int_bndFV: nx, ny, nz  = ',nx, ny, nz
         write(*,*) 'paste_int_bndFV: nxo,nyo,nzo = ',nxo,nyo,nzo
         write(*,*) 'paste_int_bndFV: min_i, min_j, min_k  =',
     *                                        min_i, min_j, min_k
         write(*,*) 'paste_int_bndFV: max_i, max_j, max_k  =',
     *                                        max_i, max_j, max_k
         write(*,*) 'paste_int_bndFV: min_io,min_jo,min_ko =',
     *                                        min_io,min_jo,min_ko
         write(*,*) 'paste_int_bndFV: max_io,max_jo,max_ko =',
     *                                        max_io,max_jo,max_ko
         write(*,*) 'paste_int_bndFV: ib,ie  = ',ib,ie
         write(*,*) 'paste_int_bndFV: jb,je  = ',jb,je
         write(*,*) 'paste_int_bndFV: kb,ke  = ',kb,ke
         write(*,*) 'paste_int_bndFV: ',lx,rx,ly,ry,lz,rz
         write(*,*) 'paste_int_bndFV: fvflag = ',fvflag
      end if

      !
      !   Z-face
      !
      if (min_ko .eq. 1 .and. max_ko .eq. nzo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'paste_int_bndFV: skipping z-faces'
         goto 10
      else if (min_ko .eq. 1 ) then
         if (ltrace) write(*,*) 'paste_int_bndFV: upper    z-face'
         if (fvflag) then
            ks = max( max_k - (bwidth-1)-1, (max_k+min_k)/2 )
            kf = max_k - 1
         else
            ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
            kf = max_k
         end if
      else if (max_ko .eq. nzo) then
         if (ltrace) write(*,*) 'paste_int_bndFV: lower    z-face'
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'paste_int_bndFV: both     z-faces'
         both = .true.
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
      end if
  5   continue
      if (ltrace) write(*,*) 'paste_int_bndFV: ks = ',ks,kf
      do k = ks, kf
         do j = jb, je
            do i = ib, ie
               ii           = ii + 1
               field(i,j,k) = storage(ii)
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            ks = max( max_k - (bwidth-1)-1, (max_k+min_k)/2 )
            kf = max_k - 1
         else
            ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
            kf = max_k
         end if
         both = .false.
         goto  5
      end if

      if (ltrace) write(*,*) 'paste_int_bndFV: z-faces done, ii=',ii

      !
      !   Y-face
      !
 10   if (min_jo .eq. 1 .and. max_jo .eq. nyo) then
         !  z-face boundaries fall on boundaries
         !  (not interiors) so skip to y-faces
         if (ltrace) write(*,*) 'paste_int_bndFV: skipping y-faces'
         goto 20
      else if (min_jo .eq. 1 ) then
         if (ltrace) write(*,*) 'paste_int_bndFV: upper    y-face'
         if (fvflag) then
            js = max( max_j - (bwidth-1)-1, (max_j+min_j)/2 )
            jf = max_j - 1
         else
            js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
            jf = max_j
         end if
      else if (max_jo .eq. nyo) then
         if (ltrace) write(*,*) 'paste_int_bndFV: lower    y-face'
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'paste_int_bndFV: both     y-faces'
         both = .true.
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
      end if
 15   continue
      if (ltrace) write(*,*) 'paste_int_bndFV: js = ',js,jf
      do j = js, jf
         do k = kb, ke
            do i = ib, ie
               ii           = ii + 1
               field(i,j,k) = storage(ii)
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            js = max( max_j - (bwidth-1)-1, (max_j+min_j)/2 )
            jf = max_j - 1
         else
            js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
            jf = max_j
         end if
         both = .false.
         goto 15
      end if

      if (ltrace) write(*,*) 'paste_int_bndFV: y-faces done, ii=',ii

      !
      !   X-faces:
      !
 20   if (min_io .eq. 1 .and. max_io .eq. nxo) then
         !  x-face boundaries fall on boundaries
         !  (not interiors) so skip 
         if (ltrace) write(*,*) 'paste_int_bndFV: skipping x-faces'
         goto 30
      else if (min_io .eq. 1 ) then
         if (ltrace) write(*,*) 'paste_int_bndFV: upper    x-face'
         if (fvflag) then
            is = max( max_i - (bwidth-1)-1, (max_i+min_i)/2 )
            if = max_i - 1
         else
            is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
            if = max_i
         end if
      else if (max_io .eq. nxo) then
         if (ltrace) write(*,*) 'paste_int_bndFV: lower    x-face'
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      else 
         if (ltrace) write(*,*) 'paste_int_bndFV: both     x-faces'
         both = .true.
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
      end if
  25  continue
      if (ltrace) write(*,*) 'paste_int_bndFV: is = ',is,if
      do i = is, if
         do k = kb, ke
            do j = jb, je
               ii           = ii + 1
               field(i,j,k) = storage(ii)
            end do
         end do
      end do
      if (both) then
         if (fvflag) then
            is = max( max_i - (bwidth-1)-1, (max_i+min_i)/2 )
            if = max_i - 1
         else
            is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
            if = max_i
         end if
         both = .false.
         goto 25
      end if

  30  continue

      if (ltrace) write(*,*) 'paste_int_bndFV: x-faces done, ii=',ii

      if (ii.ne.length) then
      !if (.false.) then
         write(*,*) 'paste_int_bndFV: ii     = ', ii
         write(*,*) 'paste_int_bndFV: length = ', length
         write(*,*) 'paste_int_bndFV: ib,ie  = ',ib,ie
         write(*,*) 'paste_int_bndFV: jb,je  = ',jb,je
         write(*,*) 'paste_int_bndFV: kb,ke  = ',kb,ke
         write(*,*) 'Problem in paste_int_bndFV'
         write(*,*) 'Problem in paste_int_bndFV'
         write(*,*) 'Problem in paste_int_bndFV'
         write(*,*) 'Problem in paste_int_bndFV'
         write(*,*) 'Problem in paste_int_bndFV'
         write(*,*) 'Problem in paste_int_bndFV'
         write(*,*) 'Problem in paste_int_bndFV'
         write(*,*) 'Problem in paste_int_bndFV'
         write(*,*) 'Problem in paste_int_bndFV'
         write(*,*) 'Problem in paste_int_bndFV'
         stop
      end if

      return
      end    ! END: paste_int_bndC

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_index:                                                               cc
cc                                                                            cc
cc           Finds the index into an array for which that element             cc
cc           matches an inputted value. Assumes array ordered.                cc
cc           If not present, returns -1.                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function find_index(array_in, value, nx )
      implicit     none
      integer      nx
      real(kind=8) array_in(nx), value
      include     'largesmall.inc'
      integer      i, i_lo, i_hi
      logical      double_equal
      external     double_equal
      !
      logical      ltrace
      parameter (  ltrace = .false. )

      !
      ! Initialize to value not present
      !
      find_index = -1

      !
      ! If value less than beginning or greater than end,
      ! then the value is not present:
      !
      if (value .lt. array_in(1) -SMALLNUMBER .or.
     *    value .gt. array_in(nx)+SMALLNUMBER     ) then
      !if (value .lt. array_in(1) .or. value .gt. array_in(nx)) then
         if (ltrace) write(*,*) 'find_index: value not in range'
         if (ltrace) write(*,*) 'find_index: value     = ',value
         if (ltrace) write(*,*) 'find_index: array(1 ) = ',array_in(1)
         if (ltrace) write(*,*) 'find_index: array(nx) = ',array_in(nx)
         return
      end if

      !
      ! Binary search:
      !
      i_lo = 1
      i_hi = nx
 100  if (i_hi-i_lo .le. 1) then
         if (ltrace) write(*,*) 'find_index: value not in array'
         return
      end if
      i      = ( i_lo + i_hi ) / 2
      if (ltrace) write(*,*) 'find_index: i, array_in(i) = ',
     *                        i, array_in(i)
      !if (array_in(i) .eq. value) then
      if (double_equal(array_in(i),value)) then
         find_index = i
         if (ltrace) write(*,*) 'find_index: value found at: ',i
         return
      else if (array_in(i) .gt. value) then
         i_hi = i
      else
         i_lo = i
      end if
      goto 100

      return
      end      ! END: find_index
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_indexN:                                                              cc
cc                                                                            cc
cc           Finds the NEAREST index into an array for which that element     cc
cc           matches an inputted value. Assumes array ordered.                cc
cc           If element is less    than range, returns 1                      cc
cc           If element is greater than range, returns nx                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function find_indexN(array_in, value, nx )
      implicit     none
      integer      nx
      real(kind=8) array_in(nx), value
      integer      i, i_lo, i_hi
      logical      double_equal
      external     double_equal
      !
      logical      ltrace
      parameter (  ltrace = .false. )

      !
      ! Initialize to value not present
      !
      find_indexN = -1

      if (value .le. array_in(1)) then
         find_indexN = 1
         return
      else if (value .ge. array_in(nx)) then
         find_indexN = nx
         return
      end if

      !
      ! Binary search:
      !
      i_lo = 1
      i_hi = nx
 100  if (i_hi-i_lo .le. 1) then
         if (ltrace) write(*,*) 'find_indexN: value in between',
     *    array_in(i_lo),array_in(i_hi), i_lo,i_hi
         find_indexN = i_lo+NINT(  (value         -array_in(i_lo) ) /
     *                             (array_in(i_hi)+array_in(i_lo) )   )
         return
      end if
      i      = ( i_lo + i_hi ) / 2
      if (ltrace) write(*,*) 'find_index: i, array_in(i) = ',
     *                        i, array_in(i)
      if (double_equal(array_in(i),value)) then
         find_indexN = i
         if (ltrace) write(*,*) 'find_indexN: value found at: ',i
         return
      else if (array_in(i) .gt. value) then
         i_hi = i
      else
         i_lo = i
      end if
      goto 100

      return
      end      ! END: find_indexN

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  calc_bnd_lenB:                                                            cc
cc                                                                            cc
cc           Computes length of vector storage necessary to hold              cc
cc           boundary of intersection of two grids.                           cc
cc           Used in conjunction with paste_int_bnd and copy_int_bnd()        cc
cc           to sync grids.                                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function calc_bnd_lenb(
     *                        min_i, max_i, min_j, max_j, min_k, max_k,
     *                        min_io,max_io,min_jo,max_jo,min_ko,max_ko,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx,ny,nz,nxo,nyo,nzo,bwidth)
      implicit    none
      integer     min_i, max_i, min_j, max_j,  min_k, max_k,
     *            min_io,max_io,min_jo,max_jo, min_ko,max_ko,
     *            lx,rx, ly,ry, lz,rz,
     *            nx,ny,nz, nxo,nyo,nzo, bwidth, tmpbwidth
      integer     ib,ie, jb,je, kb,ke
      integer     if, is,
     *            jf, js,
     *            kf, ks,
     *            width
      integer     calc_bnd_len

      logical     ltrace
      parameter ( ltrace = .false. )


      ! try *not* adjusting for bwidth:
      if (.true.) then
         tmpbwidth = 0
         !write(*,*) 'calc_bnd_len:using no bwidth'
      else
         tmpbwidth = bwidth
      end if

      if (ltrace) then
         write(*,*) 'calc_bnd_len: '
         write(*,*) 'calc_bnd_len: min_i,j,k=',min_i,min_j,min_k
         write(*,*) 'calc_bnd_len: max_i,j,k=',max_i,max_j,max_k
         write(*,*) 'calc_bnd_len: nx,ny,nz =',nx,ny,nz
         write(*,*) 'calc_bnd_len: min_i/j/ko=',min_io,min_jo,min_ko
         write(*,*) 'calc_bnd_len: max_i/j/ko=',max_io,max_jo,max_ko
         write(*,*) 'calc_bnd_len: nx/y/zo   =',nxo,nyo,nzo
         write(*,*) 'calc_bnd_len: lx, rx    =',lx,rx
         write(*,*) 'calc_bnd_len: ly, ry    =',ly,ry
         write(*,*) 'calc_bnd_len: lz, rz    =',lz,rz
         write(*,*) 'calc_bnd_len: bwidth    =',bwidth

         write(*,*) 'calc_bnd_len: nx,ny,nz =',nx,ny,nz
      end if

      !
      ! Paste points but need to avoid simply switching
      ! points from one grid to another ... in other words
      ! points shouldn't be copied and pasted twice:
      !    unless grids share a common boundary, offset
      !    the regions which get pasted.
      !
      if (min_i .eq. 1 .and. min_io .eq. 1 .and. lx .eq. 1) then
         ib = min_i
      else if (min_i .eq. 1 ) then
         ib = min_i + tmpbwidth
      else
         ib = min_i
      end if
      if (max_i .eq. nx .and. max_io .eq. nxo .and. rx .eq. 1) then
         ie = max_i
      else if (max_i .eq. nx) then
         ie = max_i - tmpbwidth
         if (ie.lt.ib) ie = ib
      else
         ie = max_i
         if (ie.lt.ib) ib = ie
      end if
      !
      if (min_j .eq. 1 .and. min_jo .eq. 1 .and. ly .eq. 1) then
         jb = min_j
      else if (min_j .eq. 1 ) then
         jb = min_j + tmpbwidth
      else
         jb = min_j
      end if
      if (max_j .eq. ny .and. max_jo .eq. nyo .and. ry .eq. 1) then
         je = max_j
      else if (max_j .eq. ny) then
         je = max_j - tmpbwidth
         if (je.lt.jb) je = jb
      else
         je = max_j
         if (je.lt.jb) jb = je
      end if
      !
      if (min_k .eq. 1 .and. min_ko .eq. 1 .and. lz .eq. 1) then
         kb = min_k
      else if (min_k .eq. 1 ) then
         kb = min_k + tmpbwidth
      else
         kb = min_k
      end if
      if (max_k .eq. nz .and. max_ko .eq. nzo .and. rz .eq. 1) then
         ke = max_k
      else if (max_k .eq. nz) then
         ke = max_k - tmpbwidth
         if (ke.lt.kb) ke = kb
      else
         ke = max_k
         if (ke.lt.kb) kb = ke
      end if

      !
      ! Add up the points for each face:
      !
      calc_bnd_len = 0

      !
      !   Z faces:
      !
      !if (min_k .eq. 1 .and. max_k .eq. nz) then
      if ( (min_k .eq.1.and.max_k .eq.nz ) .or.
     *     (min_ko.eq.1.and.max_ko.eq.nzo)     ) then
         goto 10
      else if (min_k .eq. 1 ) then
         ! upper z-face
         ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
         kf = max_k
         width = max( 0, kf-ks+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (je-jb+1)*(ie-ib+1)
      else if (max_k .eq. nz) then
         ! lower z-face
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
         width = max( 0, kf-ks+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (je-jb+1)*(ie-ib+1)
      else
         ! both z-faces
         !
         ks = max( max_k - (bwidth-1), (max_k+min_k)/2 + 1 )
         kf = max_k
         width = max( 0, kf-ks+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (je-jb+1)*(ie-ib+1)
         !
         ks = min_k
         kf = min( ks + (bwidth-1), (max_k+min_k)/2 - 1 )
         width = max( 0, kf-ks+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (je-jb+1)*(ie-ib+1)
      end if
      !
      !   Y faces:
      !
  10  if ( (min_j .eq.1.and.max_j .eq.ny ) .or.
     *     (min_jo.eq.1.and.max_jo.eq.nyo)     ) then
! 10  if (min_j .eq. 1 .and. max_j .eq. ny) then
         goto 20
      else if (min_j .eq. 1 ) then
         ! upper y-face
         js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
         jf = max_j
         width = max( 0, jf-js+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (ke-kb+1)*(ie-ib+1)
      else if (max_j .eq. ny) then
         ! lower y-face
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
         width = max( 0, jf-js+1 )
         calc_bnd_len = calc_bnd_len 
     *                    +      width * (ke-kb+1)*(ie-ib+1)
      else
         ! both y-faces
         js = max( max_j - (bwidth-1), (max_j+min_j)/2 + 1 )
         jf = max_j
         width = max( 0, jf-js+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (ke-kb+1)*(ie-ib+1)
         !
         js = min_j
         jf = min( js + (bwidth-1), (max_j+min_j)/2 - 1 )
         width = max( 0, jf-js+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (ke-kb+1)*(ie-ib+1)
      end if
      !
      !   X faces:
      !
  20  if ( (min_i .eq.1.and.max_i .eq.nx ) .or.
     *     (min_io.eq.1.and.max_io.eq.nxo)     ) then
! 20  if (min_i .eq. 1 .and. max_i .eq. nx) then
         goto 30
      else if (min_i .eq. 1 ) then
         ! upper x-face
         is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
         if = max_i
         width = max( 0, if-is+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (ke-kb+1)*(je-jb+1)
      else if (max_i .eq. nx) then
         ! lower x-face
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
         width = max( 0, if-is+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (ke-kb+1)*(je-jb+1)
      else
         ! both x-faces
         is = max( max_i - (bwidth-1), (max_i+min_i)/2 + 1 )
         if = max_i
         width = max( 0, if-is+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (ke-kb+1)*(je-jb+1)
         !
         is = min_i
         if = min( is + (bwidth-1), (max_i+min_i)/2 - 1 )
         width = max( 0, if-is+1 )
         calc_bnd_len = calc_bnd_len
     *                    +      width * (ke-kb+1)*(je-jb+1)
      end if

  30  continue

      if (ltrace) then
         write(*,*) 'calc_bnd_len: ib-->ie:   ',ib,ie
         write(*,*) 'calc_bnd_len: jb-->je:   ',jb,je
         write(*,*) 'calc_bnd_len: kb-->ke:   ',kb,ke
         write(*,*) 'calc_bnd_len: is-->if:   ',is,if
         write(*,*) 'calc_bnd_len: js-->jf:   ',js,jf
         write(*,*) 'calc_bnd_len: ks-->kf:   ',ks,kf
         write(*,*) 'calc_bnd_len: Returning: ', calc_bnd_len
         write(*,*) ''
      end if

      calc_bnd_lenb = calc_bnd_len

      return
      end    ! END: calc_bnd_lenb

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  remove_negs:                                                              cc
cc                  Reset all negative values to zero.                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine remove_negs(field, nx,   ny,       nz   )
      implicit    none
      integer     nx, ny, nz
      real(kind=8)      field(nx,ny,nz)
      integer     i,j,k
      real(kind=8)      tmp
      logical     ltrace
      parameter ( ltrace = .false. )

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         tmp = field(i,j,k)
         if (tmp .lt. 0) field(i,j,k) = 0.d0
      end do
      end do
      end do

      return
      end        ! END: subroutine remove_negs


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  integ2_z:                                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integ2_z( field, integrand, h, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz),
     .           integrand(nx,ny,nz),  h
      integer    i,j,k

      !
      ! Assume constant of integration is zero:
      !
      do j = 1, ny
      do i = 1, nx
         field(i,j,1) = 0.d0
      end do
      end do

      do k = 1, nz-1
      do j = 1, ny
      do i = 1, nx
         field(i,j,k+1) =   field(i,j,k) + 0.5d0 * h
     *                * (integrand(i,j,k+1) + integrand(i,j,k))
      end do
      end do
      end do

      return
      end    ! END: integ2_z

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  find_inflect:                                                             cc
cc          finds largest inflection point in 1d vector                       cc
cc          Used for clustering ala Berger and Rigtousis '91.                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_inflect( u1, max, max_i, nx)
      implicit    none
      integer     max_i, nx
      real(kind=8)      max, u1(nx)
      integer     i
      real(kind=8)      tmp, tmp2
      logical     ltrace
      parameter ( ltrace      = .false. )


      max     = 0.d0
      tmp     = u1(1)

      do i = 2, nx
         if ( u1(i)*tmp .lt. 0.d0 ) then
            !
            ! If opposite sign, then an inflection point:
            !
            tmp2            = abs(u1(i)-tmp)
            if (ltrace) write(*,*) ' Inflection pt found: ',i,tmp2
            if (tmp2 .gt. max) then
               if (ltrace) write(*,*) ' New max'
               max   = tmp2
               max_i = i - 1
            end if
         end if
         tmp = u1(i)
      end do


      return
      end      ! END: find_inflect

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  compute_lap:                                                              cc
cc          computes Laplacian terms from signature arrays...                 cc
cc          Used for clustering ala Berger and Rigtousis '91.                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_lap( lap, sig, nx)
      implicit    none
      integer     nx
      real(kind=8)      lap(nx), sig(nx)
      integer     i

      lap(1)  = 0.d0
      lap(nx) = 0.d0
      do i = 2, nx-1
         lap(i) = 2.d0*sig(i) - sig(i-1) - sig(i+1)
      end do


      return
      end      ! END: compute_lap

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mat_buffer:                                                               cc
cc              Takes a flag array, and flags points withing "buffer" pts     cc
cc              NB: if a point is disallowed, it will not get buffered/flaggedcc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mat_buffer( flag, work, buffer, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz, buffer
      real(kind=8)     flag(nx,ny,nz), work(nx,ny,nz)
      integer    i,j,k, ii,jj,kk

      !
      ! Copy flag to work array:
      !
      call mat_copy3d( flag, work, nx,ny,nz)

      !
      ! Remake flag array with buffer:
      !
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if ( NINT(work(i,j,k)) .eq. 1) then
                  do kk = k-buffer, k+buffer
                  do jj = j-buffer, j+buffer
                  do ii = i-buffer, i+buffer
                     if (ii.ge.1 .and. ii.le.nx .and.
     *                   jj.ge.1 .and. jj.le.ny .and.
     *                   kk.ge.1 .and. kk.le.nz ) then
!    *                   NINT(work(ii,jj,kk)).ge.0) then
                         if (NINT(work(ii,jj,kk)).ge.0)
     *                      flag(ii,jj,kk) = work(i,j,k)
                     end if
                  end do
                  end do
                  end do
               end if
            end do
         end do
      end do

      return
      end    ! END: mat_buffer

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  field_smooth:                                                             cc
cc            Trying out different smoothing for multigrid method.            cc
cc            Work is temporary storage (if needed)                           cc
cc            CHR is the characteristic array to stay away from AMR           cc
cc                because we need this level to still sit smoothly w/ parent  cc
cc            num is an integer describing on what attempt we are in case we  cc
cc                want to keep trying with more aggressive forms of smoothing cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine field_smooth( field, work, chr, num, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz, num
      real(kind=8)     field(nx,ny,nz), work(nx,ny,nz), chr(nx,ny,nz)
      include   'chr.inc'
      integer    i,j,k, pass, kswitch, ib

      !
      ! Copy flag to work array:
      !
      call mat_copy3d( field, work, nx,ny,nz)

      !
      !
      if (num .eq. 1 .or. num .eq. 2 .and. .false.) then
         ! Do either 1 or 2 passes of red-black in which we take average:
         do pass = 1, num
         !do pass = 1, 2
            kswitch = pass
            do k = 2, nz-1
               ib = kswitch + 1
               do j = 2, ny-1
                  do i = ib, nx-1, 2
                     field(i,j,k) = 0.5d0 * work(i,j,k)
     *                             +0.5d0*(1.d0/6.d0)*(
     *                                          work(i+1,j  ,k  )
     *                                         +work(i-1,j  ,k  )
     *                                         +work(i,  j+1,k  )
     *                                         +work(i,  j-1,k  )
     *                                         +work(i,  j,  k+1)
     *                                         +work(i,  j,  k-1)
     *                                        )
                  end do
                  ib = 5 - ib
               end do
               kswitch = 3 - kswitch
            end do
         end do
      else
         !
         ! Remake flag array with buffer:
         !
         do pass = 1, 1
            do k = 2, nz-1
               do j = 2, ny-1
                  do i = 2, nx-1
                  field(i,j,k) = 0.5d0 * work(i,j,k)
     *                          +0.5d0*(1.d0/6.d0)*(
     *                                    work(i+1,j  ,k  )
     *                                   +work(i-1,j  ,k  )
     *                                   +work(i,  j+1,k  )
     *                                   +work(i,  j-1,k  )
     *                                   +work(i,  j,  k+1)
     *                                   +work(i,  j,  k-1)
     *                                  )
                  end do
               end do
            end do
            call mat_copy3d( field, work, nx,ny,nz)
         end do
      end if

      return
      end    ! END: field_smooth

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  smooth3D:                                                                 cc
cc            Using storage array, replace each point by the average over     cc
cc            some stencil.                                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smooth3D( field, work, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), work(nx,ny,nz)
      integer    i,j,k

      !
      ! Copy flag to work array:
      !
      call mat_copy3d( field, work, nx,ny,nz)

      !
      ! Remake flag array with buffer:
      !
      do k = 2, nz-1
         do j = 2, ny-1
            do i = 2, nx-1
               field(i,j,k) = 0.5d0 * work(i,j,k)
     *                       +0.5d0*(1.d0/6.d0)*(
     *                                    work(i+1,j  ,k  )
     *                                   +work(i-1,j  ,k  )
     *                                   +work(i,  j+1,k  )
     *                                   +work(i,  j-1,k  )
     *                                   +work(i,  j,  k+1)
     *                                   +work(i,  j,  k-1)
     *                                  )
            end do
         end do
      end do

      return
      end    ! END: mat_buffer

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  double_equal:                                                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function double_equal( first, second )
      implicit    none
      real(kind=8)      first, second
      include       'largesmall.inc'
      !real(kind=8)      SMALLNUMBER
      !parameter       ( SMALLNUMBER = 1.0d-14)
      logical           ltrace
      parameter       ( ltrace      = .false. )

      double_equal = abs(first-second) .lt. SMALLNUMBER

      if (ltrace) then
         write(*,*)'double_equal:first,second',first,second,double_equal
         write(*,*)'double_equal:diff ',abs(first-second),SMALLNUMBER
         write(*,*)'double_equal:log',abs(first-second).lt.SMALLNUMBER
      end if

      return
      end    ! END: double_equal

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  isainf:                                                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function isainf( x )
      implicit none
      real(kind=8)   x
      include       'largesmall.inc'

      if (x.lt.LARGE_M_SMALL) then
         isainf = .false.
      else
         isainf = .true.
      end if
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  isanan:                                                                   cc
cc         This uses the built-in function isnan() which may (neh, probably)  cc
cc         is not portable. So, substitute what works on your system, or      cc
cc         just not use this routine by always returning false.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function isanan( x )
      implicit none
      real(kind=8)   x

      !if (x.le.abs(x)) then
      !if (x.ne.x) then

#if defined AIX || defined __GNUC__ || defined __PGI
      isanan = .false.
#else
      if (isnan(x)) then
         isanan = .true.
      else
         isanan = .false.
      end if
#endif

      return
      end

c     ........................................................
c
c     Pull factors of 2 out of a number:
c     Let:    x = (2^n) * P
c     Until:  min_x0 =< P =< max_x0
c
c
c     Input:     x           --- number to be factored
c                max_x0      --- maximum remainder, P
c                min_x0      --- minimum remainder, P
c
c     Output:    xc          --- remainder, P
c                nx          --- number of factors of 2, n
c                checkfacts = +1    if xc <= max_x0
c                           = -1    if xc >  max_x0
c
c     ........................................................
      integer function checkfacts(x, max_x0, min_x0, xc, nx)
      implicit   none
      integer    x,    max_x0,   min_x0, xc,    nx
      logical    ltrace
      parameter (ltrace = .false.)

      nx         = 0
      xc         = x
      do while ( (mod(xc,2).eq.0) .and. (xc.gt.min_x0) )
         nx      = nx + 1
         xc      = xc / 2
      end do

c
c     Test if grid factorizable down to max size:
c
      if ( xc .le. max_x0 ) then
         checkfacts = +1
      else
         checkfacts = -1
      end if

      if (ltrace) then
         write(*,*) 'checkfacts = ', checkfacts
         write(*,*) '         x = ', x
         write(*,*) '    max_x0 = ', max_x0
         write(*,*) '    min_x0 = ', min_x0
         write(*,*) '        xc = ', xc
         write(*,*) '        nx = ', nx
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  restrict:                                                                 cc
cc             Half-weighted with straight injection on boundary.             cc
cc             Used right now exclusively for Multigrid.                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine restrict( uc, uf, nxc, nyc, nzc)
      implicit     none
      integer      nxc, nyc, nzc
      real(kind=8) uf(2*nxc-1,2*nyc-1,2*nzc-1), uc(nxc,nyc,nzc)
      integer      ic, jc, kc, if, jf, kf
      !
      ! Straight injection is a bit easier in parallel
      ! but apparently does not allow for vcycles to 
      ! drive the residual down
      !
      logical      straightinject
      parameter (  straightinject = .false.)

      logical      ltrace
      parameter (  ltrace = .false. )
      
      if (ltrace) then
         write(*,*) 'restrict: nx/y/zc: ',nxc,nyc,nzc
         write(*,*) 'restrict: nx/y/zf: ',2*nxc-1,2*nyc-1,2*nzc-1
      end if

      do kc = 2, nzc-1
         kf = 2*kc -1
         do jc = 2, nyc-1
            jf = 2*jc -1
            do ic = 2, nxc-1
               if = 2*ic -1
               if (straightinject) then
               uc(ic,jc,kc) =                 uf(if,  jf,  kf  )
               else
               uc(ic,jc,kc) =    0.5d0    *   uf(if,  jf,  kf  )
     *                         + (1.d0/12.d0)*(
     *                                        uf(if+1,jf,  kf  )
     *                                      + uf(if-1,jf,  kf  )
     *                                      + uf(if,  jf+1,kf  )
     *                                      + uf(if,  jf-1,kf  )
     *                                      + uf(if,  jf,  kf+1)
     *                                      + uf(if,  jf,  kf-1)    )
               end if
            end do
         end do
      end do

      !
      ! Straight injection on both x-boundaries:
      !
      do ic = 1, nxc, nxc-1
         if = 2*ic -1
         do kc = 1, nzc
            kf = 2*kc -1
            do jc = 1, nyc
               jf = 2*jc -1
               uc(ic,jc,kc) = uf(if,jf,kf)
            end do
         end do
      end do
      !
      ! Straight injection on both y-boundaries:
      !
      do jc = 1, nyc, nyc-1
         jf = 2*jc -1
         do kc = 1, nzc
            kf = 2*kc -1
            do ic = 1, nxc
               if = 2*ic -1
               uc(ic,jc,kc) = uf(if,jf,kf)
            end do
         end do
      end do
      !
      ! Straight injection on both z-boundaries:
      !
      do kc = 1, nzc, nzc-1
         kf = 2*kc -1
         do jc = 1, nyc
            jf = 2*jc -1
            do ic = 1, nxc
               if = 2*ic -1
               uc(ic,jc,kc) = uf(if,jf,kf)
            end do
         end do
      end do

      if (ltrace) write(*,*) 'restrict: i/j/kf: ',if,jf,kf

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  zerobounddeco:                                                            cc
cc             Zero-out the boundary of a field where the chr array           cc
cc             indicates a domain decomposition boundary.                     cc
cc             Takes as input the width in grid points to zero-out.           cc
cc             Used in the multigrid post-sweep computation.                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zerobounddeco( u, chr, width, nx, ny, nz)
      implicit     none
      integer      nx, ny, nz, width
      real(kind=8) u(nx,ny,nz), chr(nx,ny,nz)
      include     'chr.inc'
      integer      i,j,k, l, m
      
         if (NINT(chr(nx,ny/2,nz/2)).eq.chr_deco_bdy) then
            do k = 1, nz
            do j = 1, ny
               do i = nx-(width-1), nx
                  u(i,j,k) = 0.d0
               end do
            end do
            end do
         end if
         if (NINT(chr(1,ny/2,nz/2)).eq.chr_deco_bdy) then
            do k = 1, nz
            do j = 1, ny
               do i = 1, 1+(width-1)
                  u(i,j,k) = 0.d0
               end do
            end do
            end do
         end if
         !
         if (NINT(chr(nx/2,ny,nz/2)).eq.chr_deco_bdy) then
            do k = 1, nz
            do i = 1, nx
               do j = ny-(width-1), ny
                  u(i,j,k) = 0.d0
               end do
            end do
            end do
         end if
         if (NINT(chr(nx/2,1,nz/2)).eq.chr_deco_bdy) then
            do k = 1, nz
            do i = 1, nx
               do j = 1, 1+(width-1)
                  u(i,j,k) = 0.d0
               end do
            end do
            end do
         end if
         !
         if (NINT(chr(nx/2,ny/2,nz)).eq.chr_deco_bdy) then
            do j = 1, ny
            do i = 1, nx
               do k = nz-(width-1), nz
                  u(i,j,k) = 0.d0
               end do
            end do
            end do
         end if
         if (NINT(chr(nx/2,ny/2,1)).eq.chr_deco_bdy) then
            do j = 1, ny
            do i = 1, nx
               do k = 1, 1+(width-1)
                  u(i,j,k) = 0.d0
               end do
            end do
            end do
         end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  prolong:                                                                  cc
cc             Interpolate for multigrid.                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prolong( uf, uc, nxc, nyc, nzc)
      implicit     none
      integer      nxc, nyc, nzc
      real(kind=8) uf(2*nxc-1,2*nyc-1,2*nzc-1), uc(nxc,nyc,nzc)
      integer      ic, jc, kc, if, jf, kf
      integer      nxf, nyf, nzf, gw
      integer      mini, maxi, minj, maxj, mink, maxk
      integer      r
      parameter (  r = 2)
      
      !
      ! Careful here. Need to make sure the routine below;
      !    (1) is forced to use linear interpolation
      !    (2) pays no attention to "gw"
      !
      nxf  = 2*nxc-1
      nyf  = 2*nyc-1
      nzf  = 2*nzc-1
      gw   = 1
      mini = 1
      maxi = nxf
      minj = 1
      maxj = nyf
      mink = 1
      maxk = nzf

      call init_fieldGWa( uc, uf, nxc,nyc,nzc,
     *                         mini,minj,mink,
     *                         maxi,maxj,maxk,
     *                         nxf, nyf, nzf, r, gw)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  copy_toplace:                                                             cc
cc                                                                            cc
cc           Takes data stored in storage and copies into a grid              cc
cc           function with lower corner position at displacement/anchor       cc
cc           point: (ax,ay,az).                                               cc
cc           For use with multigrid routines.                                 cc
cc           NB: Assumes that points from originating grid located on         cc
cc           domain decomposition boundaries have been "masked" out by        cc
cc           setting them to LARGENUMBER. This saves us from passing around   cc
cc           that grid's CHR() array.                                         cc
cc           NB: Used with copy_overlap().                                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copy_toplace( field, storage, ax,ay,az, nx, ny, nz,
     *                                                   nxt,nyt,nzt)
      implicit    none
      integer     nx,ny,nz, nxt, nyt, nzt, ax,ay,az
      real(kind=8)storage(nxt,nyt,nzt), field(nx,ny,nz)
      integer     i,j,k, it,jt,kt
      include       'largesmall.inc'
      logical     ltrace
      parameter ( ltrace = .false. )

      do kt = 1, nzt
         k  = kt + az - 1
         do jt = 1, nyt
            j  = jt + ay - 1
            do it = 1, nxt
               i  = it + ax - 1
               if (storage(it,jt,kt) .lt. LARGE_M_SMALL) then
                  field(i,j,k) = storage(it,jt,kt)
               end if
            end do
         end do
      end do


      return
      end    ! END: copy_toplace


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  maskLargeh2h:                                                             cc
cc             Mask out (with LARGENUMBER) field values which are             cc
cc             characterized as DD boundaries. Assumes chr()                  cc
cc             inputted array is defined on a grid of resolution 'h'          cc
cc             whereas the field to be masked is on   resolution '2h.'        cc
cc             Ultimately, copy_toplace() knows to then ignore these values.  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine maskLargeh2h( uc, chrf, nxc, nyc, nzc, bwidth)
      implicit     none
      integer      nxc, nyc, nzc, bwidth
      real(kind=8) chrf(2*nxc-1,2*nyc-1,2*nzc-1), uc(nxc,nyc,nzc)
      include     'chr.inc'
      integer      ic, jc, kc, if, jf, kf, nxf,nyf,nzf
      include     'largesmall.inc'

      logical      ltrace
      parameter (  ltrace = .false. )
      
      if (ltrace) then
         write(*,*) 'maskLargeh2h: '
      end if

      nxf = 2*nxc-1
      nyf = 2*nyc-1
      nzf = 2*nzc-1

      !
      ! Do simple thing first:
      !
      if (NINT(chrf(1,nyf/2,nzf/2)) .eq. CHR_deco_bdy) then
         do kc = 1, nzc
            kf = 2*kc -1
            do jc = 1, nyc
               jf = 2*jc -1
               do if = 1, bwidth, 2
                   ic = NINT( 0.5d0*(if+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
               end do
            end do
         end do
      end if
      if (NINT(chrf(nxf,nyf/2,nzf/2)) .eq. CHR_deco_bdy) then
         do kc = 1, nzc
            kf = 2*kc -1
            do jc = 1, nyc
               jf = 2*jc -1
               do if = nxf, nxf-bwidth+1, -2
                   ic = NINT( 0.5d0*(if+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
               end do
            end do
         end do
      end if
      !
      if (NINT(chrf(nxf/2,1,nzf/2)) .eq. CHR_deco_bdy) then
         do kc = 1, nzc
            kf = 2*kc -1
            do ic = 1, nxc
               if = 2*ic -1
               do jf = 1, bwidth, 2
                   jc = NINT( 0.5d0*(jf+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
               end do
            end do
         end do
      end if
      if (NINT(chrf(nxf/2,nyf,nzf/2)) .eq. CHR_deco_bdy) then
         do kc = 1, nzc
            kf = 2*kc -1
            do ic = 1, nxc
               if = 2*ic -1
               do jf = nyf, nyf-bwidth+1, -2
                   jc = NINT( 0.5d0*(jf+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
               end do
            end do
         end do
      end if
      !
      if (NINT(chrf(nxf/2,nyf/2,1)) .eq. CHR_deco_bdy) then
         do jc = 1, nyc
            jf = 2*jc -1
            do ic = 1, nxc
               if = 2*ic -1
               do kf = 1, bwidth, 2
                   kc = NINT( 0.5d0*(kf+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
               end do
            end do
         end do
      end if
      if (NINT(chrf(nxf/2,nyf/2,nzf)) .eq. CHR_deco_bdy) then
         do jc = 1, nyc
            jf = 2*jc -1
            do ic = 1, nxc
               if = 2*ic -1
               do kf = nzf, nzf-bwidth+1, -2
                   kc = NINT( 0.5d0*(kf+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
               end do
            end do
         end do
      end if
      !
      return
      !
      do kc = 1, nzc
         kf = 2*kc -1
         do jc = 1, nyc
            jf = 2*jc -1
            if (NINT(chrf(1,jf,kf)) .eq. CHR_deco_bdy) then
                do if = 1, bwidth, 2
                   ic = NINT( 0.5d0*(if+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
                end do
            end if
            if (NINT(chrf(nxf,jf,kf)) .eq. CHR_deco_bdy) then
                do if = nxf, nxf-bwidth+1, -2
                   ic = NINT( 0.5d0*(if+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
                end do
            end if
         end do
      end do
      !
      do kc = 1, nzc
         kf = 2*kc -1
         do ic = 1, nxc
            if = 2*ic -1
            if (NINT(chrf(if,1,kf)) .eq. CHR_deco_bdy) then
                do jf = 1, bwidth, 2
                   jc = NINT( 0.5d0*(jf+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
                end do
            end if
            if (NINT(chrf(if,nyf,kf)) .eq. CHR_deco_bdy) then
                do jf = nyf, nyf-bwidth+1, -2
                   jc = NINT( 0.5d0*(jf+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
                end do
            end if
         end do
      end do
      !
      do jc = 1, nyc
         jf = 2*jc -1
         do ic = 1, nxc
            if = 2*ic -1
            if (NINT(chrf(if,jf,1)) .eq. CHR_deco_bdy) then
                do kf = 1, bwidth, 2
                   kc = NINT( 0.5d0*(kf+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
                end do
            end if
            if (NINT(chrf(if,jf,nzf)) .eq. CHR_deco_bdy) then
                do kf = nzf, nzf-bwidth+1, -2
                   kc = NINT( 0.5d0*(kf+1) )
                   uc(ic,jc,kc) =    LARGENUMBER
                end do
            end if
         end do
      end do


      return
      end     ! END: maskLargeh2h

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  output_resid:                                                             cc
cc               Only gets residual from first grid in hierarchy for now.     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine output_resid(level,string,resnrm)
      implicit none
      integer       level
      character*(*) string
      real(kind=8)  resnrm

      integer       indent
      parameter  (  indent = 2 )

      character*128 spaces
      data          spaces  / ' ' /

      write(*,100)  spaces(1:1+(level-1)*indent),level,string, resnrm

 100  format(a,'Level ',i2,' ',a,' sweep ',t50,1p,e9.2)

      return
      end       ! END: output_resid

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  field_out1d:                                                              cc
cc                    Quickly output a field to SDF.                          cc
cc                    Used for debugging, else use grid_output_sdf().         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine field_out1d(field,time, name,mix,max,nx,id)
      implicit      none
      integer       nx, id
      real(kind=8)  field(nx),time
      real(kind=8)  mix, max
      character*(*) name
      !
      character(18) prefix, gft_name
      integer       rank
      parameter (   rank = 1)
      integer       gft_rc, gft_shape(rank), lenpre
      real(kind=8)  bbox(2*rank)

      integer       gft_out_bbox, mystringlength
      external      gft_out_bbox, mystringlength

      logical      ltrace
      parameter (  ltrace = .false. )

      gft_shape(1) = nx

      bbox(1) = mix
      bbox(2) = max

      call int2str(id,prefix)

      lenpre              = mystringlength(prefix)
      gft_name(1:lenpre)  = prefix(1:lenpre)
      gft_name(lenpre+1:) = name

      if (ltrace) then
          write(*,*) 'field_out1d: name:     ',name
          write(*,*) 'field_out1d: id:       ',id
          write(*,*) 'field_out1d: time:     ',time
          write(*,*) 'field_out1d: prefix:   ',prefix
          write(*,*) 'field_out1d: lenpre:   ',lenpre
          write(*,*) 'field_out1d: gft_name: ',gft_name
          write(*,*) 'field_out1d: nx:       ',nx
          write(*,*) 'field_out1d: mi/max:   ',mix,max
      end if

#ifdef AIX
      if (ltrace) write(*,*) 'field_out1d: using AIX method:'
      gft_rc       = gft_out_bbox(trim(gft_name)/ /CHAR(0), 
     *                            %VAL(time), gft_shape,
     *                            %VAL(rank), bbox, field)
#else
      if (ltrace) write(*,*) 'field_out1d: using default method:'
      gft_rc       = gft_out_bbox(gft_name, time,     gft_shape,
     *                            rank, bbox, field)
#endif

      if (ltrace) then
          write(*,*) 'field_out1d: Done.     ',gft_rc
      end if

      return
      end        ! END: field_out1d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  field_out2d:                                                              cc
cc                    Quickly output a field to SDF.                          cc
cc                    Used for debugging, else use grid_output_sdf().         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine field_out2d(field,time, name,mix,max,miy,may, 
     *                       nx,ny,id)
      implicit      none
      integer       nx,ny, id
      real(kind=8)  field(nx,ny),time
      real(kind=8)  mix, max, miy, may
      character*(*) name
      !
      character(18) prefix, gft_name
      integer       rank
      parameter (   rank = 2)
      integer       gft_rc, gft_shape(rank), lenpre
      real(kind=8)  bbox(2*rank)

      integer       gft_out_bbox, mystringlength
      external      gft_out_bbox, mystringlength

      logical      ltrace
      parameter (  ltrace = .false. )

      gft_shape(1) = nx
      gft_shape(2) = ny

      bbox(1) = mix
      bbox(2) = max
      bbox(3) = miy
      bbox(4) = may

      call int2str(id,prefix)

      lenpre              = mystringlength(prefix)
      gft_name(1:lenpre)  = prefix(1:lenpre)
      gft_name(lenpre+1:) = name

      if (ltrace) then
          write(*,*) 'field_out2d: name:     ',name
          write(*,*) 'field_out2d: id:       ',id
          write(*,*) 'field_out2d: time:     ',time
          write(*,*) 'field_out2d: prefix:   ',prefix
          write(*,*) 'field_out2d: lenpre:   ',lenpre
          write(*,*) 'field_out2d: gft_name: ',gft_name
      end if

#ifdef AIX
      gft_rc       = gft_out_bbox(trim(gft_name)/ /CHAR(0), 
     *                            %VAL(time), gft_shape,
     *                            %VAL(rank), bbox, field)
#else
      gft_rc       = gft_out_bbox(gft_name, time,     gft_shape,
     *                            rank, bbox, field)
#endif

      return
      end        ! END: field_out2d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  field_out3d:                                                              cc
cc                    Quickly output a field to SDF.                          cc
cc                    Used for debugging, else use grid_output_sdf().         cc
cc                    Does not record coordinate info, just maps to 0..1.     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine field_out3d(field,time, name,mix,max,miy,may, 
     *                       miz,maz,nx,ny,nz,id)
      implicit      none
      integer       nx,ny,nz, id
      real(kind=8)  field(nx,ny,nz),time
      real(kind=8)  mix, max, miy, may, miz, maz
      character*(*) name
      !
      character(18) prefix, gft_name
      integer       rank
      parameter (   rank = 3)
      integer       gft_rc, gft_shape(rank), lenpre
      real(kind=8)  bbox(2*rank)

      integer       gft_out_bbox, mystringlength
      external      gft_out_bbox, mystringlength

      logical      ltrace
      parameter (  ltrace = .false. )

      gft_shape(1) = nx
      gft_shape(2) = ny
      gft_shape(3) = nz

      bbox(1) = 0.d0
      bbox(2) = 1.d0
      bbox(3) = 0.d0
      bbox(4) = 1.d0
      bbox(5) = 0.d0
      bbox(6) = 1.d0
      bbox(1) = mix
      bbox(2) = max
      bbox(3) = miy
      bbox(4) = may
      bbox(5) = miz
      bbox(6) = maz

      call int2str(id,prefix)

      !Does not work for some reason:
      !gft_name            = prefix//name
      !
      !lenpre              = len(prefix)
      lenpre              = mystringlength(prefix)
      gft_name(1:lenpre)  = prefix(1:lenpre)
      gft_name(lenpre+1:) = name

      if (ltrace) then
          write(*,*) 'field_out3d: name:     ',name
          write(*,*) 'field_out3d: id:       ',id
          write(*,*) 'field_out3d: time:     ',time
          write(*,*) 'field_out3d: prefix:   ',prefix
          write(*,*) 'field_out3d: lenpre:   ',lenpre
          write(*,*) 'field_out3d: gft_name: ',gft_name
      end if

#ifdef AIX
      gft_rc       = gft_out_bbox(trim(gft_name)/ /CHAR(0), 
     *                           %VAL(time), gft_shape,
     *                           %VAL(rank), bbox, field)
#else
      gft_rc       = gft_out_bbox(gft_name, time,     gft_shape,
     *                            rank, bbox, field)
#endif

      return
      end        ! END: field_out3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  test_round_robin:                                                         cc
cc                    Used just to test the round_robin() routine by          cc
cc                    generating the appropriate round robin tables.          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine test_round_robin(numprocs)
      implicit   none
      integer    numprocs
      integer    numrounds, round, i
      integer    round_robin
      external   round_robin

      if (mod(numprocs,2).eq.0) then
         ! Even:
         numrounds = numprocs - 1
      else
         ! Odd:
         numrounds = numprocs
      end if

      write(*,*) 'test_round_robin: numprocs = ',numprocs
      do round = 1, numrounds
         write(*,*) 'Round: ',round
         do i = 0,numprocs-1
         !do i = 0,numprocs/2-1
            write(*,*) '    Match: ',i,round_robin(round,i,numprocs)
         end do
         !write(*,*) '        .....ensure 1-to-1.....'
         !do i = numprocs-1, numprocs/2, -1
            !write(*,*) '    Match: ',i,round_robin(round,i,numprocs)
         !end do
      end do
      write(*,*) ''

      return
      end        ! END: test_round_robin

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  round_robin:                                                              cc
cc                    Returns next "opponent" in a round robin.               cc
cc                    Uses N as the fixed point, and rotates a polygon        cc
cc                    with numprocs-1 vertices around it.                     cc
cc                    Assumes "contestants" are numbered 0...numprocs-1       cc
cc                    If numprocs is odd, returns -1 for the bye.             cc
cc                    See for example:                                        cc
cc                    http://www.devenezia.com/downloads/round-robin/         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function round_robin(round, myid, numprocs)
      implicit   none
      integer    round, myid, numprocs
      integer    N        ! highest rank proc
      integer    roundm1  ! round-1
      logical    ltrace
      parameter (ltrace = .false.)

      if (mod(numprocs,2).eq.0) then
         ! Even:
         N = numprocs - 1
      else
         ! Odd:
         N = numprocs
      end if

      if (round.lt.1) then
         write(*,*) 'round_robin: Round must be >0',round
         stop
      end if

      roundm1 = round - 1
      if (myid .le. N - roundm1*2) then
         round_robin = N - myid - roundm1*2
      else if (myid .eq. N - roundm1) then
         round_robin = N
      else if (myid .eq. N) then
         round_robin = N - roundm1
      else 
         round_robin = 2*N - roundm1*2 - myid
 10      if (round_robin .lt. 0) then
            round_robin = round_robin + N
            goto 10
         end if
      end if

      if (round_robin .eq. numprocs) then
         ! In this case, numprocs is odd, and this processor gets a bye:
         round_robin = -1
      end if

      if (ltrace) write(*,99) myid, 'round_robin: Proc: ',round_robin

  99  format('[',I4,'] ',A,3I5)

      return
      end        ! END: round_robin

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  log_baseN:                                                                cc
cc                    Computes log of a number with general base.             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function log_baseN(number,base)
      implicit none
      integer  base
      real(kind=8) number

      if (number.le.0) then
         write(*,*) 'log_baseN: Cannot take log of neg number'
         log_baseN = -1.d0
      else
         log_baseN = log(number) / log(1.d0*base)
      end if

      return
      end        ! END: log_baseN

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  log_baseN_INT:                                                            cc
cc                    Computes log of an integer returning an integer.        cc
cc                    So returns the maximum power of the (integer) base      cc
cc                    which is less than or equal to the given integer.       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function log_baseN_INT(number,base)
      implicit none
      integer  base, number
      real(kind=8) numberfloat

      if (number.le.0) then
         write(*,*) 'log_baseN_INT: Cannot take log of neg number'
         log_baseN_INT = -1
      else
         !
         ! Add 0.5d0 so that taking the log of powers of 2
         ! are safe with respect to all the casting of types.
         !    for example, the log of 8 is 2, but we could
         !    end up with 1.999999 which would then cast down to 1.
         !    However, adding a fraction before the cast, means we'll safely
         !    get 2.???? something and end up with the right answer,
         !    and we'll avoid having the log of 7 return with 3.
         !
         numberfloat   = 1.d0*number + 0.5d0
         log_baseN_INT = INT( log(numberfloat) / log(1.d0*base) )
      end if

      return
      end        ! END: log_baseN_INT

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc set_bfield:                                                                cc
cc             Set the components of the b-field so that                      cc
cc             divergenceless constraint is satisfied.                        cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_bfield(fieldbz, fieldby, fieldbx, x, y, z,
     *                      nx,ny,nz, ampbx,ampby, abx,aby, bbx,bby, cb)
      implicit     none
      integer      nx,ny,nz
      real(kind=8) fieldbz(nx,ny,nz),fieldby(nx,ny,nz),fieldbx(nx,ny,nz)
      real(kind=8) x(nx,ny,nz),      y(nx,ny,nz),      z(nx,ny,nz)
      real(kind=8) ampbx,ampby,abx,aby,bbx,bby, cb
      integer      i, j, k

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         fieldbx(i,j,k) = ampbx * exp( -x(i,j,k)**2/abx**2 )
     *                          * exp( -y(i,j,k)**2/bbx**2 )
     *                          * exp( -z(i,j,k)**2/cb**2  )
     *                          * z(i,j,k)/cb**2
         fieldby(i,j,k) = ampby * exp( -x(i,j,k)**2/aby**2 )
     *                          * exp( -y(i,j,k)**2/bby**2 )
     *                          * exp( -z(i,j,k)**2/cb**2  )
     *                          * z(i,j,k)/cb**2
         fieldbz(i,j,k) = -exp(-z(i,j,k)**2/cb**2)*(
     *                    (x(i,j,k) / abx**2) *
     *                    ampbx * exp( -x(i,j,k)**2/abx**2 )
     *                          * exp( -y(i,j,k)**2/bbx**2 )
     *                  + (y(i,j,k) / bby**2) *
     *                    ampby * exp( -x(i,j,k)**2/aby**2 )
     *                          * exp( -y(i,j,k)**2/bby**2 )
     *                          )
      end do
      end do
      end do

      return
      end       ! END: set_bfield

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc set_efield:                                                                cc
cc             Set the components of the e-field so that                      cc
cc             divergenceless constraint is satisfied.                        cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_efield(fieldez, fieldey, fieldex, fieldphi, x,y,z,
     *               kappa, nx,ny,nz, ampbx,ampby, abx,aby, bbx,bby, cb)
      implicit     none
      integer      nx,ny,nz
      real(kind=8) fieldez(nx,ny,nz),fieldey(nx,ny,nz),fieldex(nx,ny,nz)
      real(kind=8) fieldphi(nx,ny,nz), kappa
      real(kind=8) x(nx,ny,nz),      y(nx,ny,nz),      z(nx,ny,nz)
      real(kind=8) ampbx,ampby,abx,aby,bbx,bby, cb
      integer      i, j, k

      !
      ! First, do exactly the same thing as we would do for
      ! the magnetic field:
      !
      call set_bfield(fieldez, fieldey, fieldex, x,y,z,
     *                nx,ny,nz, ampbx,ampby, abx,aby, bbx,bby, cb)

      !
      ! Now, since each of those fields is really: exp(-2*kappa*phi)* \vec E
      ! we need to multiply by exp(+2*kappa*phi):
      !
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         fieldex(i,j,k) = exp(2.d0*kappa*fieldphi(i,j,k))*fieldex(i,j,k)
         fieldey(i,j,k) = exp(2.d0*kappa*fieldphi(i,j,k))*fieldey(i,j,k)
         fieldez(i,j,k) = exp(2.d0*kappa*fieldphi(i,j,k))*fieldez(i,j,k)
      end do
      end do
      end do

      return
      end       ! END: set_efield

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc mask_cube:                                                                 cc
cc             Sets a mask function to -1 within a specified                  cc
cc             cube.                                                          cc
cc             The convention within hyper for the mask is:                   cc
cc                  -1.d0     ---     masked                                  cc
cc                 > 0.d0     --- not masked                                  cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mask_cube(mask, radius, centerx, centery, centerz,
     *                     x, y, z, nx, ny, nz )
      implicit     none
      integer      nx,ny,nz
      real(kind=8) mask(nx,ny,nz), x(nx), y(ny), z(nz)
      real(kind=8) radius, centerx, centery, centerz
      integer      i, j, k
      real(kind=8) r2

      logical    ltrace
      parameter (ltrace = .false.)

      r2 = radius**2

      if (ltrace) then
         write(*,*) 'mask_cube: radius      = ',radius
         write(*,*) 'mask_cube: r2          = ',r2
         write(*,*) 'mask_cube: centerx/y/z = ',centerx,centery,centerz
      end if

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if ( ( (x(i)-centerx)**2 .lt. r2 ) .and.
     *        ( (y(j)-centery)**2 .lt. r2 ) .and.
     *        ( (z(k)-centerz)**2 .lt. r2 ) ) then
            !if (ltrace) write(*,*) 'mask_cube: masking pt: ',i,j,k
            mask(i,j,k) = -1.d0
!         else   
!           mask(i,j,k) = +1.d0
         end if
      end do
      end do
      end do


      return
      end       ! END: mask_cube 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc mask_distance_movedx:                                                      cc
cc             Determine the distance in X that the mask has moved            cc
cc             since its initial setup.                                       cc
cc             NB: assumes entire masked region is on grid,                   cc
cc                 else, if no distance can be dtermined, returns -1.d0       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function mask_distance_movedx(mask, radius,
     *                         centerx, x, nx, ny, nz )
      implicit     none
      integer      nx,ny,nz
      real(kind=8) mask(nx,ny,nz), x(nx)
      real(kind=8) radius, centerx
      real(kind=8) newcenter
      integer      i, j, k, min_i, max_i

      logical    ltrace
      parameter (ltrace = .false.)


      min_i = nx
      max_i =  1
      do i = 1, nx
         do k = 1, nz
         do j = 1, ny
            if ( mask(i,j,k) .lt. 0) then
               if (i.lt. min_i) min_i = i
               if (i.gt. max_i) max_i = i
            end if
         end do
         end do
      end do

      if (min_i .eq. nx .and. max_i .eq. 1) then
         mask_distance_movedx = -1.d0
      else
         newcenter            = 0.5d0*(x(max_i)+x(min_i))
         mask_distance_movedx = newcenter - centerx
      end if

      if (ltrace) then
       write(*,*)'mask_distance_movedx: min_i    =',min_i
       write(*,*)'mask_distance_movedx: max_i    =',max_i
       write(*,*)'mask_distance_movedx: nx       =',nx
       write(*,*)'mask_distance_movedx: newcenter=',newcenter
       write(*,*)'mask_distance_movedx: centerx  =',centerx
       write(*,*)'mask_distance_movedx: distance =',mask_distance_movedx
      end if

      return
      end       ! END: mask_distance_movedx

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc remask_cube:                                                               cc
cc             Sets a mask function to +1 in a specified region               cc
cc             for the purposes of moving a mask. Hence, should               cc
cc             be called in conjuntion with grid_repopulate() which           cc
cc             then leaves the mask in the state of being all                 cc
cc             -1's and 0's.                                                  cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine remask_cube(mask, radius, centerx, centery, centerz,
     *                     x, y, z, nx, ny, nz )
      implicit     none
      integer      nx,ny,nz
      real(kind=8) mask(nx,ny,nz), x(nx), y(ny), z(nz)
      real(kind=8) radius, centerx, centery, centerz
      integer      i, j, k
      real(kind=8) r2
      real(kind=8) fuzz_scale
      parameter  ( fuzz_scale = 0.1d-4 )

      logical    ltrace
      parameter (ltrace = .false.)

      r2 = radius**2

      if (ltrace) then
         write(*,*)'remask_cube: radius      = ',radius
         write(*,*)'remask_cube: r2          = ',r2
         write(*,*)'remask_cube: centerx/y/z = ',centerx,centery,centerz
      end if

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if ( ( (x(i)-centerx)**2 .lt. r2 ) .and.
     *        ( (y(j)-centery)**2 .lt. r2 ) .and.
     *        ( (z(k)-centerz)**2 .lt. r2 ) ) then
            mask(i,j,k) = +10.d0
         end if
      end do
      end do
      end do


      return
      end       ! END: remask_cube 

!------------------------------------------------------------------------
!
!     Copied from had/src/hyper/mask.f90 (and converted to F77)
!
!------------------------------------------------------------------------
      subroutine myfixup_mask(mask, nx,ny,nz)
      implicit   none
      integer      nx,ny,nz
      real(kind=8) mask(nx,ny,nz)
      integer      i, j, k, im, iml
      logical    ltrace2
      parameter (ltrace2 = .false.)

      if (ltrace2) then
         write(*,*)'myfixup_mask: .......'
         write(*,*)'myfixup_mask: nx/y/z: ',nx,ny,nz
         write(*,*)'myfixup_mask: .......'
      end if

      !--------------------------------------------------------
      ! March along x and check for mask boundaries
      !--------------------------------------------------------
      do k = 1, nz
      do j = 1, ny

         do i = 2, nx
            im  = nint(mask(i,j,k))
            iml = nint(mask(i-1,j,k))

            if (im*iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  mask(i-1,j,k) = mask(i-1,j,k) - 1.0
                  if (ltrace2) then
                     write(0,*)'AFound boundary in x ',mask(i-1,j,k)
                  end if
               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  mask(i,j,k) = mask(i,j,k) - 1.0
                  if (ltrace2) then
                     write(0,*)'BFound boundary in x ',mask(i,j,k)
                  end if
               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'FIXUP_MASK >> mask problem in x sweep',
     *                                    im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do

      !--------------------------------------------------------
      ! March along y and check for mask boundaries
      !--------------------------------------------------------
      do k = 1, nz
      do i = 1, nx

         do j = 2, ny
            im  = nint(mask(i,j,k))
            iml = nint(mask(i,j-1,k))

            if (im*iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  mask(i,j-1,k) = mask(i,j-1,k) - 1.0
                  if (ltrace2) then
                     write(0,*)'AFound boundary in y ',mask(i,j-1,k)
                  end if
               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  mask(i,j,k) = mask(i,j,k) - 1.0
                  if (ltrace2) then
                     write(0,*)'BFound boundary in y ',mask(i,j,k)
                  end if
               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'FIXUP_MASK >> mask problem in y sweep',
     *                                    im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do

      !--------------------------------------------------------
      ! March along y and check for mask boundaries
      !--------------------------------------------------------
      do j = 1, ny
      do i = 1, nx

         do k = 2, nz
            im  = nint(mask(i,j,k))
            iml = nint(mask(i,j,k-1))

            if (im*iml .lt. 0) then
               !
               ! We have a boundary
               !
               if (im .eq. 1 .and. iml .lt. 0) then
                  !
                  ! The mask is to the left
                  !
                  mask(i,j,k-1) = mask(i,j,k-1) - 1.0
                  if (ltrace2) then
                     write(0,*)'AFound boundary in z ',mask(i,j,k-1)
                  end if
               else if (im .lt. 0 .and. iml .eq. 1) then
                  !
                  !  The mask is to the right
                  !
                  mask(i,j,k) = mask(i,j,k) - 1.0
                  if (ltrace2) then
                     write(0,*)'BFound boundary in z ',mask(i,j,k)
                  end if
               else
                  !
                  ! Unexpected values in the mask
                  !
                  write(0,*)'FIXUP_MASK >> mask problem in z sweep',
     *                                    im, iml,mask(i,j,k)
               end if
            end if

         end do

      end do
      end do

      if (ltrace2) write(*,*)'myfixup_mask: Finished. '

      return
      end subroutine myfixup_mask

!------------------------------------------------------------------------
!
!     Un-does what fixup_mask() and restores mask to MASKED
!     and UNMASKED values
!
!------------------------------------------------------------------------
      subroutine unfixup_mask(mask, nx,ny,nz)
      implicit   none
      integer      nx,ny,nz
      real(kind=8) mask(nx,ny,nz)
      include     'mask.inc'
      integer      i, j, k
      logical    ltrace2
      parameter (ltrace2 = .false.)

      if (ltrace2) then
         write(*,*)'unfixup_mask: .......'
         write(*,*)'unfixup_mask: nx/y/z: ',nx,ny,nz
         write(*,*)'unfixup_mask: .......'
      end if

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         !
         if (mask(i,j,k) .lt. 0) then
            mask(i,j,k) =  1.d0 * BHMASK_MASKED
         else if (mask(i,j,k) .gt. 0) then
            mask(i,j,k) =  1.d0 * BHMASK_UNMASKED
         else 
            write(*,*) 'unfixup_mask: Why is mask zero?'
         end if
         !
      end do
      end do
      end do

      if (ltrace2) write(*,*)'unfixup_mask: Finished. '

      return
      end subroutine unfixup_mask

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  getabsmax                                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8)  function getabsmax( u, nx)
      implicit    none
      integer     nx
      real(kind=8)      u(nx)
      integer     i

      getabsmax = abs(u(1))
      do i = 2, nx
         if ( abs(u(i)) .gt. getabsmax) then
            getabsmax = abs(u(i))
         end if
      end do

      return
      end      ! END: getabsmax


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc    grow_bboxall:                                                           cc
cc              Expand bboxes to include a ghostwidth but                     cc
cc              be careful not to:                                            cc
cc                  1) include any disallowed points                          cc
cc                  2) go past the bounds of the grid                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grow_bboxall( flagarray, value, width, bmini,bmaxi,
     *                      bminj,bmaxj,bmink,bmaxk, nx,ny,nz,numb)
      implicit     none
      integer      numb
      integer      nx, ny, nz, width, bmini(numb),bmaxi(numb),
     *             bminj(numb),bmaxj(numb), bmink(numb),bmaxk(numb)
      real(kind=8) flagarray(nx,ny,nz), value
      !
      integer     i, j, k, l, num, m, n
      logical      double_equal
      external     double_equal
      logical      allowed, alreadyoverlap
      logical      ltrace
      parameter (  ltrace = .false. )
      ! It may make sense not to let grids overlap
      ! in which case we can avoid causing it here
      ! although if grids already overlap this routine will not undue it
      logical      DONOTCAUSEOVERLAP
      parameter (  DONOTCAUSEOVERLAP = .true. )

      if (ltrace) then
         write(*,*) 'grow_bboxall: value:     ',value
         write(*,*) 'grow_bboxall: width:     ',width
         write(*,*) 'grow_bboxall: numb:      ',numb
         write(*,*) 'grow_bboxall: nx/y/z:    ',nx,ny,nz
         do l = 1, numb
         write(*,*) 'grow_bboxall: bmin/axi:  ',bmini(l),bmaxi(l)
         write(*,*) 'grow_bboxall: bmin/axj:  ',bminj(l),bmaxj(l)
         write(*,*) 'grow_bboxall: bmin/axk:  ',bmink(l),bmaxk(l)
         end do
      end if

      ! increase each box by 1 each time so that
      ! the boxes are increased as close to symmetrically as possible
      do m = 1, width
         do l = 1, numb
            !
            ! Small x direction
            !
            allowed = .false.
            i       = bmini(l)-1
            if (i.gt.0) then
               allowed = .true.
               do k = bmink(l), bmaxk(l)
                  do j = bminj(l), bmaxj(l)
                     if (double_equal(flagarray(i,j,k),value)) then
                        allowed = .false.
                     end if
                  end do
               end do
               if (allowed .and. DONOTCAUSEOVERLAP) then
                  ! check if this would cause overlap
                  do n = 1, numb
                     if (n.ne.l) then
                        ! only check against other boxes obviously:
                        alreadyoverlap =(bmaxk(n).lt.bmaxk(l) .and. 
     .                                   bmaxk(n).gt.bmink(l) ).or.
     .                                  (bmaxk(l).lt.bmink(n) .and.
     .                                   bmaxk(l).gt.bmink(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxj(n).lt.bmaxj(l) .and. 
     .                                   bmaxj(n).gt.bminj(l) ).or.
     .                                  (bmaxj(l).lt.bmaxj(n) .and.
     .                                   bmaxj(l).gt.bminj(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxi(n).lt.bmaxi(l) .and. 
     .                                   bmaxi(n).gt.bmini(l) ).or.
     .                                  (bmaxi(l).lt.bmaxi(n) .and.
     .                                   bmaxi(l).gt.bmini(n) )
                        ! If they already overlap, then do not worry
                        if (.not.alreadyoverlap) then
                           allowed     = .not. (
     .                                  (bmaxi(n).le.bmaxi(l) .and. 
     .                                   bmaxi(n).ge.bmini(l) ).or.
     .                                  (bmaxi(l).le.bmaxi(n) .and.
     .                                   bmaxi(l).ge.bmini(n) ) )
                        end if
                     end if   ! end if n!=l
                  end do
               end if    ! end if DONOTCAUSEOVERLAP
            end if       ! end if i>0
            if (allowed) then
               bmini(l) = bmini(l) - 1
               if(ltrace)write(*,*)'grow_bbox: Changing bmini:',bmini(l)
            end if
            !
            ! Large x direction
            !
            allowed = .false.
            i       = bmaxi(l)+1
            if (i.lt.nx) then
               allowed = .true.
               do k = bmink(l), bmaxk(l)
                  do j = bminj(l), bmaxj(l)
                     if (double_equal(flagarray(i,j,k),value)) then
                        allowed = .false.
                     end if
                  end do
               end do
               if (allowed .and. DONOTCAUSEOVERLAP) then
                  ! check if this would cause overlap
                  do n = 1, numb
                     if (n.ne.l) then
                        ! only check against other boxes obviously:
                        alreadyoverlap =(bmaxk(n).lt.bmaxk(l) .and. 
     .                                   bmaxk(n).gt.bmink(l) ).or.
     .                                  (bmaxk(l).lt.bmink(n) .and.
     .                                   bmaxk(l).gt.bmink(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxj(n).lt.bmaxj(l) .and. 
     .                                   bmaxj(n).gt.bminj(l) ).or.
     .                                  (bmaxj(l).lt.bmaxj(n) .and.
     .                                   bmaxj(l).gt.bminj(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxi(n).lt.bmaxi(l) .and. 
     .                                   bmaxi(n).gt.bmini(l) ).or.
     .                                  (bmaxi(l).lt.bmaxi(n) .and.
     .                                   bmaxi(l).gt.bmini(n) )
                        ! If they already overlap, then do not worry
                        if (.not.alreadyoverlap) then
                           allowed     = .not. (
     .                                  (bmaxi(n).le.bmaxi(l) .and. 
     .                                   bmaxi(n).ge.bmini(l) ).or.
     .                                  (bmaxi(l).le.bmaxi(n) .and.
     .                                   bmaxi(l).ge.bmini(n) ) )
                        end if
                     end if   ! end if n!=l
                  end do
               end if    ! end if DONOTCAUSEOVERLAP
            end if       ! end if i>0
            if (allowed) then
               bmaxi(l) = bmaxi(l) + 1
               if(ltrace)write(*,*)'grow_bbox: Changing bmaxi:',bmaxi(l)
            end if
            !
            ! Small y direction
            !
            allowed = .false.
            j       = bminj(l)-1
            if (j.gt.0) then
               allowed = .true.
               do k = bmink(l), bmaxk(l)
                  do i = bmini(l), bmaxi(l)
                     if (double_equal(flagarray(i,j,k),value)) then
                        allowed = .false.
                     end if
                  end do
               end do
               if (allowed .and. DONOTCAUSEOVERLAP) then
                  ! check if this would cause overlap
                  do n = 1, numb
                     if (n.ne.l) then
                        ! only check against other boxes obviously:
                        alreadyoverlap =(bmaxk(n).lt.bmaxk(l) .and. 
     .                                   bmaxk(n).gt.bmink(l) ).or.
     .                                  (bmaxk(l).lt.bmink(n) .and.
     .                                   bmaxk(l).gt.bmink(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxj(n).lt.bmaxj(l) .and. 
     .                                   bmaxj(n).gt.bminj(l) ).or.
     .                                  (bmaxj(l).lt.bmaxj(n) .and.
     .                                   bmaxj(l).gt.bminj(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxi(n).lt.bmaxi(l) .and. 
     .                                   bmaxi(n).gt.bmini(l) ).or.
     .                                  (bmaxi(l).lt.bmaxi(n) .and.
     .                                   bmaxi(l).gt.bmini(n) )
                        ! If they already overlap, then do not worry
                        if (.not.alreadyoverlap) then
                           allowed     = .not. (
     .                                  (bmaxj(n).le.bmaxj(l) .and. 
     .                                   bmaxj(n).ge.bminj(l) ).or.
     .                                  (bmaxj(l).le.bmaxj(n) .and.
     .                                   bmaxj(l).ge.bminj(n) ) )
                        end if
                     end if   ! end if n!=l
                  end do
               end if    ! end if DONOTCAUSEOVERLAP
            end if       ! end if i>0
            if (allowed) then
               bminj(l) = bminj(l) - 1
               if(ltrace)write(*,*)'grow_bbox: Changing bminj:',bminj(l)
            end if
            !
            ! Large y direction
            !
            allowed = .false.
            j       = bmaxj(l)+1
            if (j.lt.ny) then
               allowed = .true.
               do k = bmink(l), bmaxk(l)
                  do i = bmini(l), bmaxi(l)
                     if (double_equal(flagarray(i,j,k),value)) then
                        allowed = .false.
                     end if
                  end do
               end do
               if (allowed .and. DONOTCAUSEOVERLAP) then
                  ! check if this would cause overlap
                  do n = 1, numb
                     if (n.ne.l) then
                        ! only check against other boxes obviously:
                        alreadyoverlap =(bmaxk(n).lt.bmaxk(l) .and. 
     .                                   bmaxk(n).gt.bmink(l) ).or.
     .                                  (bmaxk(l).lt.bmink(n) .and.
     .                                   bmaxk(l).gt.bmink(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxj(n).lt.bmaxj(l) .and. 
     .                                   bmaxj(n).gt.bminj(l) ).or.
     .                                  (bmaxj(l).lt.bmaxj(n) .and.
     .                                   bmaxj(l).gt.bminj(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxi(n).lt.bmaxi(l) .and. 
     .                                   bmaxi(n).gt.bmini(l) ).or.
     .                                  (bmaxi(l).lt.bmaxi(n) .and.
     .                                   bmaxi(l).gt.bmini(n) )
                        ! If they already overlap, then do not worry
                        if (.not.alreadyoverlap) then
                           allowed     = .not. (
     .                                  (bmaxj(n).le.bmaxj(l) .and. 
     .                                   bmaxj(n).ge.bminj(l) ).or.
     .                                  (bmaxj(l).le.bmaxj(n) .and.
     .                                   bmaxj(l).ge.bminj(n) ) )
                        end if
                     end if   ! end if n!=l
                  end do
               end if    ! end if DONOTCAUSEOVERLAP
            end if       ! end if i>0
            if (allowed) then
               bmaxj(l) = bmaxj(l) + 1
               if(ltrace)write(*,*)'grow_bbox: Changing bmaxj:',bmaxj(l)
            end if
            !
            ! Small z direction
            !
            allowed = .false.
            k       = bmink(l)-1
            if (k.gt.0) then
               allowed = .true.
               do j = bminj(l), bmaxj(l)
                  do i = bmini(l), bmaxi(l)
                     if (double_equal(flagarray(i,j,k),value)) then
                        allowed = .false.
                     end if
                  end do
               end do
               if (allowed .and. DONOTCAUSEOVERLAP) then
                  ! check if this would cause overlap
                  do n = 1, numb
                     if (n.ne.l) then
                        ! only check against other boxes obviously:
                        alreadyoverlap =(bmaxk(n).lt.bmaxk(l) .and. 
     .                                   bmaxk(n).gt.bmink(l) ).or.
     .                                  (bmaxk(l).lt.bmink(n) .and.
     .                                   bmaxk(l).gt.bmink(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxj(n).lt.bmaxj(l) .and. 
     .                                   bmaxj(n).gt.bminj(l) ).or.
     .                                  (bmaxj(l).lt.bmaxj(n) .and.
     .                                   bmaxj(l).gt.bminj(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxi(n).lt.bmaxi(l) .and. 
     .                                   bmaxi(n).gt.bmini(l) ).or.
     .                                  (bmaxi(l).lt.bmaxi(n) .and.
     .                                   bmaxi(l).gt.bmini(n) )
                        ! If they already overlap, then do not worry
                        if (.not.alreadyoverlap) then
                           allowed     = .not. (
     .                                  (bmaxk(n).le.bmaxk(l) .and. 
     .                                   bmaxk(n).ge.bmink(l) ).or.
     .                                  (bmaxk(l).le.bmaxk(n) .and.
     .                                   bmaxk(l).ge.bmink(n) ) )
                        end if
                     end if   ! end if n!=l
                  end do
               end if    ! end if DONOTCAUSEOVERLAP
            end if       ! end if i>0
            if (allowed) then
               bmink(l) = bmink(l) - 1
               if(ltrace)write(*,*)'grow_bbox: Changing bmink:',bmink(l)
            end if
            !
            ! Large z direction
            !
            allowed = .false.
            k       = bmaxk(l)+1
            if (k.lt.nz) then
               allowed = .true.
               do j = bminj(l), bmaxj(l)
                  do i = bmini(l), bmaxi(l)
                     if (double_equal(flagarray(i,j,k),value)) then
                        allowed = .false.
                     end if
                  end do
               end do
               if (allowed .and. DONOTCAUSEOVERLAP) then
                  ! check if this would cause overlap
                  do n = 1, numb
                     if (n.ne.l) then
                        ! only check against other boxes obviously:
                        alreadyoverlap =(bmaxk(n).lt.bmaxk(l) .and. 
     .                                   bmaxk(n).gt.bmink(l) ).or.
     .                                  (bmaxk(l).lt.bmink(n) .and.
     .                                   bmaxk(l).gt.bmink(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxj(n).lt.bmaxj(l) .and. 
     .                                   bmaxj(n).gt.bminj(l) ).or.
     .                                  (bmaxj(l).lt.bmaxj(n) .and.
     .                                   bmaxj(l).gt.bminj(n) )
                        if (alreadyoverlap)
     .                  alreadyoverlap =(bmaxi(n).lt.bmaxi(l) .and. 
     .                                   bmaxi(n).gt.bmini(l) ).or.
     .                                  (bmaxi(l).lt.bmaxi(n) .and.
     .                                   bmaxi(l).gt.bmini(n) )
                        ! If they already overlap, then do not worry
                        if (.not.alreadyoverlap) then
                           allowed     = .not. (
     .                                  (bmaxk(n).le.bmaxk(l) .and. 
     .                                   bmaxk(n).ge.bmink(l) ).or.
     .                                  (bmaxk(l).le.bmaxk(n) .and.
     .                                   bmaxk(l).ge.bmink(n) ) )
                        end if
                     end if   ! end if n!=l
                  end do
               end if    ! end if DONOTCAUSEOVERLAP
            end if       ! end if i>0
            if (allowed) then
               bmaxk(l) = bmaxk(l) + 1
               if(ltrace)write(*,*)'grow_bbox: Changing bmaxk:',bmaxk(l)
            end if
            !
         end do   ! end loop over numb
      end do   ! end loop over width

      return
      end      ! END: grow_bboxall

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc    grow_bbox:                                                              cc
cc              Expand bboxes to include a ghostwidth but                     cc
cc              be careful not to:                                            cc
cc                  1) include any disallowed points                          cc
cc                  2) go past the bounds of the grid                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grow_bbox( flagarray, value, width, bmini,bmaxi,
     *                      bminj,bmaxj,bmink,bmaxk, nx,ny,nz)
      implicit     none
      integer      nx, ny, nz, width, bmini,bmaxi,
     *             bminj,bmaxj, bmink,bmaxk
      real(kind=8) flagarray(nx,ny,nz), value
      !
      integer     i, j, k, num
      logical      double_equal
      external     double_equal
      logical      allowed
      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace) then
         write(*,*) 'grow_bbox: value:     ',value
         write(*,*) 'grow_bbox: width:     ',width
         write(*,*) 'grow_bbox: bmin/axi:  ',bmini,bmaxi
         write(*,*) 'grow_bbox: bmin/axj:  ',bminj,bmaxj
         write(*,*) 'grow_bbox: bmin/axk:  ',bmink,bmaxk
         write(*,*) 'grow_bbox: nx/y/z:    ',nx,ny,nz
      end if

      if (.false.  ) then
         !
         ! No disallowed points to worry about:
         !
         bmini = bmini - width
         bmaxi = bmaxi + width
         if (bmini .lt.  1) bmini =  1
         if (bmaxi .gt. nx) bmaxi = nx
         !
         bminj = bminj - width
         bmaxj = bmaxj + width
         if (bminj .lt.  1) bminj =  1
         if (bmaxj .gt. ny) bmaxj = ny
         !
         bmink = bmink - width
         bmaxk = bmaxk + width
         if (bmink .lt.  1) bmink =  1
         if (bmaxk .gt. nz) bmaxk = nz
         !
      else
         !
         ! Try to expand out "width" points
         ! but do not include disallowed points
         !
         ! X min:
         !
         num     = 0
   10    i       = bmini - num - 1
         if (i.lt. 1) goto 15
         allowed = .true.
         do k = bmink, bmaxk
            do j = bminj, bmaxj
               if (double_equal(flagarray(i,j,k),value)) then
                  allowed = .false.
               end if
            end do
         end do
         if (allowed) then
            ! Increase the number by which we can grow box:
            num = num + 1
            if (num .lt. width) goto 10
         end if
   15    continue
         bmini = bmini - num
         if (ltrace) write(*,*) 'grow_bbox: Changing bmini: ',bmini,num
         !
         ! X max:
         !
         num     = 0
   20    i       = bmaxi + num + 1
         if (i.gt.nx) goto 25
         allowed = .true.
         do k = bmink, bmaxk
            do j = bminj, bmaxj
               if (double_equal(flagarray(i,j,k),value)) then
                  allowed = .false.
               end if
            end do
         end do
         if (allowed) then
            ! Increase the number by which we can grow box:
            num = num + 1
            if (num .lt. width) goto 20
         end if
   25    continue
         bmaxi = bmaxi + num
         if (ltrace) write(*,*) 'grow_bbox: Changing bmaxi: ',bmaxi,num
         !
         ! Y min:
         !
         num     = 0
   30    j       = bminj - num - 1
         if (j.lt. 1) goto 35
         allowed = .true.
         do k = bmink, bmaxk
            do i = bmini, bmaxi
               if (double_equal(flagarray(i,j,k),value)) then
                  allowed = .false.
               end if
            end do
         end do
         if (allowed) then
            ! Increase the number by which we can grow box:
            num = num + 1
            if (num .lt. width) goto 30
         end if
   35    continue
         bminj = bminj - num
         if (ltrace) write(*,*) 'grow_bbox: Changing bminj: ',bminj,num
         !
         ! Y max:
         !
         num     = 0
   40    j       = bmaxj + num + 1
         if (j.gt.ny) goto 45
         allowed = .true.
         do k = bmink, bmaxk
            do i = bmini, bmaxi
               if (double_equal(flagarray(i,j,k),value)) then
                  allowed = .false.
               end if
            end do
         end do
         if (allowed) then
            ! Increase the number by which we can grow box:
            num = num + 1
            if (num .lt. width) goto 40
         end if
   45    continue
         bmaxj = bmaxj + num
         if (ltrace) write(*,*) 'grow_bbox: Changing bmaxj: ',bmaxj,num
         !
         ! Z min:
         !
         num     = 0
   50    k       = bmink - num - 1
         if (k.lt. 1) goto 55
         allowed = .true.
         do j = bminj, bmaxj
            do i = bmini, bmaxi
               if (double_equal(flagarray(i,j,k),value)) then
                  allowed = .false.
               end if
            end do
         end do
         if (allowed) then
            ! Increase the number by which we can grow box:
            num = num + 1
            if (num .lt. width) goto 50
         end if
   55    continue
         bmink = bmink - num
         if (ltrace) write(*,*) 'grow_bbox: Changing bmink: ',bmink,num
         !
         ! Z max:
         !
         num     = 0
   60    k       = bmaxk + num + 1
         if (k.gt.nz) goto 65
         allowed = .true.
         do j = bminj, bmaxj
            do i = bmini, bmaxi
               if (double_equal(flagarray(i,j,k),value)) then
                  allowed = .false.
               end if
            end do
         end do
         if (allowed) then
            ! Increase the number by which we can grow box:
            num = num + 1
            if (num .lt. width) goto 60
         end if
   65    continue
         bmaxk = bmaxk + num
         if (ltrace) write(*,*) 'grow_bbox: Changing bmaxk: ',bmaxk,num
         !
      end if

      return
      end      ! END: grow_bbox

c-------------------------------------------------------------------------------
      integer function mygft_write_id_int(fname,pname,p,len)
c
c            Write parameters to the SDF file fname, to be read later
c            by the corresponding read command. Cannot use standard lib
c            functions since they require going into GFT_Single mode
c            which closes all the open SDF files.                     
c-------------------------------------------------------------------------------
      implicit none
      character*(*)  fname,  pname
      integer        len
      integer        p(len)
      integer        i, tmplen
      integer        shape(3)
      integer        maxp
      parameter (    maxp = 10 000 )
      real(kind=8)   pdata(maxp), coords(maxp)
      integer        gft_out_full
      external       gft_out_full
      logical        ltrace
      parameter   (  ltrace  = .false. )

      if (ltrace) then
         write(*,*) 'mygft_write_id_int: Write to:'
         write(*,*) 'mygft_write_id_int:      fname: ',fname
         write(*,*) 'mygft_write_id_int:      pname: ',pname
         write(*,*) 'mygft_write_id_int:       p(1): ',p(1)
         write(*,*) 'mygft_write_id_int:        len: ',len
      end if

      if (len .gt. maxp) then
         write(*,*) 'mygft_write_id_int: Need more storage: ',maxp, len
         write(*,*) 'mygft_write_id_int:      fname: ',fname
         write(*,*) 'mygft_write_id_int:      pname: ',pname
         write(*,*) 'mygft_write_id_int:       p(1): ',p(1)
         return
      end if

      !
      ! I am having some strange problem reading in
      ! a field of length one, so when we write it out,
      ! just waste some memory and make it of length 2
      ! and when you read it in, you just reference the first value:
      !
      if (len .eq. 1) then
         tmplen    = 2
         coords(1) = 1.d0
         coords(2) = 2.d0
         pdata(1)  = 1.d0 * p(1)
         pdata(2)  = 1.d0 * p(1)
      else
         tmplen    = len
         do i = 1, tmplen
           coords(i) = 1.d0 * i
           pdata(i)  = 1.d0 * p(i)
         end do
      end if

#ifdef AIX
      shape(1) = tmplen
      mygft_write_id_int =
     *        gft_out_full(trim(fname)/ /CHAR(0),
     *                %VAL(0.d0),shape,
     *                trim(pname)/ /CHAR(0), 
     *                %VAL(1),coords, pdata)
#else
      mygft_write_id_int =
     *        gft_out_full(fname,0.d0,tmplen,pname, 1,coords, pdata)
#endif

      return
      end

c-------------------------------------------------------------------------------
      integer function mygft_write_id_float(fname,pname,p,len)
c-------------------------------------------------------------------------------
      implicit none
      character*(*)  fname,  pname
      integer        len
      integer        tmplen
      integer        shape(3)
      real(kind=8)   p(len)
      real(kind=8)   tmpp(2)
      integer        gft_out_full
      external       gft_out_full

      !
      ! I am having some strange problem reading in
      ! a field of length one, so when we write it out,
      ! just waste some memory and make it of length 2
      ! and when you read it in, you just reference the first value:
      !   NB: SLL 02/09/12 I'm a bit worried what gft_out_full()
      !       might do when trying to access the second element of
      !       p which does not exist. So to be safe declare some
      !       tmpp(2). The int version of this already does so
      !       (probably since it had to convert things to floats
      !        anyway)
      !
      if (len .eq. 1) then
         tmplen    = 2
         tmpp(1)   = p(1)
         tmpp(2)   = tmpp(1)
#ifdef AIX
      shape(1) = tmplen
      mygft_write_id_float =
     *        gft_out_full(trim(fname)/ /CHAR(0),
     *                     %VAL(0.d0),shape,
     *                     trim(pname)/ /CHAR(0),
     *                     %VAL(1),tmpp, tmpp)
#else
      mygft_write_id_float =
     *        gft_out_full(fname,0.d0,tmplen,pname, 1,tmpp, tmpp)
#endif
      else
         tmplen    = len
#ifdef AIX
      shape(1) = tmplen
      mygft_write_id_float =
     *        gft_out_full(trim(fname)/ /CHAR(0),
     *                     %VAL(0.d0),shape,
     *                     trim(pname)/ /CHAR(0),
     *                     %VAL(1),p, p)
#else
      mygft_write_id_float =
     *        gft_out_full(fname,0.d0,tmplen,pname, 1,p, p)
#endif
      end if

      return
      end

c-------------------------------------------------------------------------------
      integer function mygft_read_id_int(fname,pname,p,len)
c-------------------------------------------------------------------------------
      implicit none
      character*(*)  fname,  pname
      integer        len
      integer        p(len)
      integer        i, tmplen
      integer        maxp
      parameter (    maxp = 10 000 )
      real(kind=8)   pdata(maxp), time, coords(maxp)
      integer        rank, shape(3), step
      character(128) cnames
      integer        proc_return_myid, myid
      external       proc_return_myid
      integer        gft_read_full, gft_read_rank, gft_read_shape
      external       gft_read_full, gft_read_rank, gft_read_shape
      integer        gft_read_brief,gft_read_sdf_ntlevs
      external       gft_read_brief,gft_read_sdf_ntlevs
      integer        ntlevs
      integer        last_int
      save           last_int
      logical        startfromlast
      logical        ltrace
      parameter   (  ltrace  = .false. )
      logical        ltrace2
      parameter   (  ltrace2 = .false. )

      myid = proc_return_myid()
      !
      ! I am having some strange problem reading in
      ! a field of length one, so when we write it out,
      ! just waste some memory and make it of length 2
      ! and when you read it in, you just reference the first value:
      !
      if (len .eq. 1) then
         tmplen = 2
      else
         tmplen = len
      end if

#ifdef AIX
      ntlevs = gft_read_sdf_ntlevs(trim(fname)/ /CHAR(0))
#else
      ntlevs            = gft_read_sdf_ntlevs(fname)
#endif

      if (ltrace) then
         write(*,99)myid,'mygft_read_id_int: Reading from: ',fname(1:10)
         write(*,99)myid,'mygft_read_id_int: Looking for:  ',pname
         write(*,98)myid,'mygft_read_id_int: last_int:     ',last_int
         write(*,98)myid,'mygft_read_id_int: ntlevs:       ',ntlevs
      end if

      mygft_read_id_int = 1
      !
      ! start after where last data set was found
      !    (assume reading done in same order as writing)
      !
      if (last_int .ge. 1 .and. last_int .lt. ntlevs) then
         step              = last_int + 1
         startfromlast     = .true.
         if(ltrace)write(*,99)myid,'mygft_read_id_int:Start from last'
      else
         step              = 1
         startfromlast     = .false.
         if(ltrace)write(*,99)myid,'mygft_read_id_int:Start from begin'
      end if
 333  continue
#ifdef AIX
      do while (mygft_read_id_int .eq. 1 .and. 
     *      trim(cnames) .ne. trim(pname)/ /CHAR(0))
#else
      do while (mygft_read_id_int .eq. 1 .and. cnames .ne. pname)
#endif
      cnames='                                                         '
#ifdef AIX
         mygft_read_id_int = gft_read_rank(trim(fname)/ /CHAR(0),
     *                             %VAL(step),rank)
         mygft_read_id_int = gft_read_shape(trim(fname)/ /CHAR(0),
     *                             %VAL(step),shape)
#else
         mygft_read_id_int = gft_read_rank(fname,step,rank)
         mygft_read_id_int = gft_read_shape(fname,step,shape)
#endif
         if (rank .eq. 1 .and. shape(1) .eq. tmplen) then
            !
            if (shape(1) .gt. maxp) then
              write(*,98)myid,'mygft_read_id_int:Not enough storage:',
     *               shape(1),maxp
              mygft_read_id_int = 0
              return
            end if
            if (ltrace2) then
               write(*,99)myid,'mygft_read_id_int:Readfrom:',fname(1:20)
               write(*,98)myid,'mygft_read_id_int:    step:',step
               write(*,98)myid,'mygft_read_id_int:    rank:',rank
               write(*,98)myid,'mygft_read_id_int:   shape:',shape(1)
               write(*,99)myid,'mygft_read_id_int:  cnames:',cnames(1:9)
               write(*,97)myid,'mygft_read_id_int:    time:',time
            end if
#ifdef AIX
         mygft_read_id_int =
     .           gft_read_full(trim(fname)/ /CHAR(0),
     .                           %VAL(step),shape,
     .                           cnames,
     .                           rank,time,
     .                           coords,   pdata           )
#else
         mygft_read_id_int =
     .           gft_read_full(  fname,    step,   shape,
     .                           cnames,   rank,   time,
     .                           coords,   pdata           )
#endif
            if (ltrace2) then
               write(*,99)myid,'mygft_read_id_int:    have:',cnames(1:9)
               write(*,98)myid,'mygft_read_id_int:  at step:',step
               write(*,97)myid,'mygft_read_id_int:    coords:',coords(1)
               write(*,97)myid,'mygft_read_id_int:     pdata:',pdata(1)
            end if
            !
         end if
         step = step + 1
      end do

#ifdef AIX
      if (trim(cnames) .eq. trim(pname)/ /CHAR(0)) then
#else
      if (cnames .eq. pname) then
#endif
         do i = 1, len
           p(i) = NINT(  pdata(i) )
         end do
         last_int = step - 1
         if (ltrace) then
            write(*,99)myid,'mygft_read_id_int: Parameter found: ',
     *                  cnames(1:20),p(1)
         end if
      else
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         write(*,99)myid,'mygft_read_id_int: Parameter not found:',pname
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !stop
         mygft_read_id_int = 0
         if (startfromlast) then
            step              = 1
            startfromlast     = .false.
            if(ltrace)write(*,99)myid,'mygft_read_id_int:ReStart'
            goto 333
         end if
      end if
  97  format('[',I4,'] ',A,3G14.7)
  98  format('[',I4,'] ',A,3I5)
  99  format('[',I4,'] ',A,A)

      return
      end

c-------------------------------------------------------------------------------
      integer function mygft_read_id_float(fname,pname,p,len)
c-------------------------------------------------------------------------------
      implicit none
      character*(*)  fname,  pname
      integer        len
      real(kind=8)   p(len)
      integer        i, tmplen
      integer        maxp
      parameter (    maxp = 10 000 )
      real(kind=8)   pdata(maxp), time, coords(maxp)
      integer        rank, shape(3), step
      character(128) cnames
      integer        proc_return_myid, myid
      external       proc_return_myid
      integer        gft_read_full, gft_read_rank, gft_read_shape
      external       gft_read_full, gft_read_rank, gft_read_shape
      integer        gft_read_brief,gft_read_sdf_ntlevs
      external       gft_read_brief,gft_read_sdf_ntlevs
      integer        ntlevs
      integer        last_float
      save           last_float
      logical        startfromlast
      logical        namesagree
      logical        ltrace
      parameter   (  ltrace  = .false. )
      logical        ltrace2
      parameter   (  ltrace2 = .false. )

      myid = proc_return_myid()
      !
      ! I am having some strange problem reading in
      ! a field of length one, so when we write it out,
      ! just waste some memory and make it of length 2
      ! and when you read it in, you just reference the first value:
      !
      if (len .eq. 1) then
         tmplen = 2
      else
         tmplen = len
      end if

#ifdef AIX
      ntlevs  = gft_read_sdf_ntlevs(trim(fname)/ /CHAR(0))
#else
      ntlevs            = gft_read_sdf_ntlevs(fname)
#endif

      if (ltrace) then
         write(*,99)myid,'mygft_read_id_float: Reading from:',fname(1:9)
         write(*,99)myid,'mygft_read_id_float: Looking for:',pname
         write(*,98)myid,'mygft_read_id_float: ntlevs:      ',ntlevs
         write(*,98)myid,'mygft_read_id_float: last_float:  ',last_float
         write(*,98)myid,'mygft_read_id_float: maxp:        ',maxp
      end if
      mygft_read_id_float = 1
      !
      ! start after where last data set was found
      !    (assume reading done in same order as writing)
      !
      if (last_float .ge. 1 .and. last_float .lt. ntlevs) then
         step              = last_float + 1
         startfromlast     = .true.
         if(ltrace)write(*,99)myid,'mygft_read_id_float:Start from last'
      else
         step              = 1
         startfromlast     = .false.
         if(ltrace)write(*,99)myid,'mygft_read_id_float:Start begin'
      end if
 333  continue
      namesagree = .false.
      do while (mygft_read_id_float .eq. 1 .and. .not.namesagree)
      cnames='                                                         '
#ifdef AIX
         mygft_read_id_float = gft_read_rank(trim(fname)/ /CHAR(0),
     *                          %VAL(step),rank)
         mygft_read_id_float = gft_read_shape(trim(fname)/ /CHAR(0),
     *                          %VAL(step),shape)
#else
         mygft_read_id_float = gft_read_rank(fname,step,rank)
         mygft_read_id_float = gft_read_shape(fname,step,shape)
#endif
         if (rank .eq. 1 .and. shape(1) .eq. tmplen) then
            !
            if (ltrace2) then
             write(*,99)myid,'mygft_read_id_float:Readfrom:',fname(1:20)
             write(*,98)myid,'mygft_read_id_float:    step:',step
             write(*,98)myid,'mygft_read_id_float:    rank:',rank
             write(*,98)myid,'mygft_read_id_float:   shape:',shape(1)
             write(*,99)myid,'mygft_read_id_float:  cnames:',cnames(1:9)
             write(*,97)myid,'mygft_read_id_float:    time:',time
             !write(*,*)'mygft_read_id_float:   namesagree: ',namesagree
            end if
            if (shape(1) .gt. maxp) then
               !
               ! Assuming this is the memory array "q",
               ! and read it in directly:
               !
#ifdef AIX
               mygft_read_id_float =
     .           gft_read_brief(trim(fname)/ /CHAR(0), 
     .                          %VAL(step), p )
#else
               mygft_read_id_float =
     .           gft_read_brief( fname,    step, p )
#endif
               !last_float = step - 1
               last_float = step
               cnames     = pname
               if (ltrace) write(*,99)myid,
     .             'mygft_read_id_float: reading directly ',cnames(1:20)
              return
            end if
#ifdef AIX
         mygft_read_id_float =
     .           gft_read_full(trim(fname)/ /CHAR(0),%VAL(step),
     .                         shape,cnames,rank,time,
     .                           coords,   pdata       )
#else
         mygft_read_id_float =
     .           gft_read_full(  fname,    step,   shape,
     .                           cnames,   rank,   time,
     .                           coords,   pdata           )
#endif
            if (ltrace2) then
               write(*,99)myid,'mygft_read_id_float:   have:',cnames
               write(*,98)myid,'mygft_read_id_float:at step:',step
               write(*,98)myid,'mygft_read_id_float:  coords:',coords(1)
               write(*,99)myid,'mygft_read_id_float:   pdata:',pdata(1)
            end if
            !
         end if
#ifdef AIX
         namesagree = trim(cnames) .eq. trim(pname)/ /CHAR(0)
#else
         namesagree = cnames .eq. pname
#endif
         step = step + 1
      end do

      if (namesagree) then
         do i = 1, len
           p(i) = pdata(i)
         end do
         last_float = step - 1
         if (ltrace) then
          write(*,99)myid,'mygft_read_id_float: Parameter found: ',
     *                  cnames(1:20),p(1)
         end if
      else
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         write(*,99)myid,'mygft_read_id_float:Parameter nt found:',pname
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         !write(*,*) '****************************************'
         if (startfromlast) then
            step              = 1
            startfromlast     = .false.
            if(ltrace)write(*,99)myid,'mygft_read_id_float:ReStart'
            goto 333
         end if
      end if
  97  format('[',I4,'] ',A,3G14.7)
  98  format('[',I4,'] ',A,3I5)
  99  format('[',I4,'] ',A,A)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interp:                                                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function interp( field, fx, fy, fz,
     *                                  i,j,k, nx,ny,nz)
      implicit     none
      integer      i,j,k, nx, ny, nz
      real(kind=8) field(nx,ny,nz),
     *             fx, fy, fz

      logical     ltrace
      parameter ( ltrace = .false. )

      if (i .gt. nx-1 .or. j .gt. ny-1 .or. k .gt. nz-1 .or.
     *    i .lt. 1    .or. j .lt. 1    .or. k .lt. 1         ) return
      !
      interp = 
     *      + (1.d0-fx) * (1.d0-fy) * (1.d0-fz) * field(i,  j,  k  )
      !
     *      + (1.d0-fx) *       fy  *       fz  * field(i,  j+1,k+1)
     *      +       fx  * (1.d0-fy) *       fz  * field(i+1,j,  k+1)
     *      +       fx  *       fy  * (1.d0-fz) * field(i+1,j+1,k  )
      !
     *      +       fx  * (1.d0-fy) * (1.d0-fz) * field(i+1,j,  k  )
     *      + (1.d0-fx) *       fy  * (1.d0-fz) * field(i,  j+1,k  )
     *      + (1.d0-fx) * (1.d0-fy) *       fz  * field(i,  j,  k+1)
      !
     *      +       fx  *       fy  *       fz  * field(i+1,j+1,k+1)

      return
      end    ! END: interp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interpderiv:                                                              cc
cc               similar to interp(), but interpolates the derivitive         cc
cc               in the chosen direction 1-x, 2-y, or 3-z                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function interpderiv( field, fx, fy, fz, inv_dx,
     *                                  i,j,k, nx,ny,nz,direction)
      implicit     none
      integer      i,j,k, nx, ny, nz, direction
      real(kind=8) field(nx,ny,nz), inv_dx,
     *             fx, fy, fz

      real(kind=8) dijk, dip1, djp1, dkp1
      real(kind=8) dip1jp1, dip1kp1, djp1kp1
      real(kind=8) dip1jp1kp1
      logical     ltrace
      parameter ( ltrace = .false. )
      logical     somethingwrong

      somethingwrong = .false.

      if (i .gt. nx-1 .or. j .gt. ny-1 .or. k .gt. nz-1 .or.
     *    i .lt. 1    .or. j .lt. 1    .or. k .lt. 1         ) return
      !
      !...................................................
      if (direction.eq.1) then
         if (fx.le.0.5 .and. i.gt.1 .and. i.lt.nx) then
            dip1       = inv_dx*( field(i+1,j,  k  )-field(i,  j,  k  ))
            dijk       = inv_dx*( field(i,  j,  k  )-field(i-1,j,  k  ))
            !
            dip1jp1    = inv_dx*( field(i+1,j+1,k  )-field(i,  j+1,k  ))
            djp1       = inv_dx*( field(i,  j+1,k  )-field(i-1,j+1,k  ))
            !
            dip1kp1    = inv_dx*( field(i+1,j,  k+1)-field(i,  j,  k+1))
            dkp1       = inv_dx*( field(i,  j,  k+1)-field(i-1,j,  k+1))
            !
            dip1jp1kp1 = inv_dx*( field(i+1,j+1,k+1)-field(i,  j+1,k+1))
            djp1kp1    = inv_dx*( field(i,  j+1,k+1)-field(i-1,j+1,k+1))
            !
            ! Stencil shifts in x-direction because of derivs:
            fx = fx + 0.5d0
         else if (fx.le.0.5 .and. i.eq.1) then
            ! At the boundary and so we could use back/forwards differncing,
            ! but instead just juse derivative a half-step from boundary:
            dip1       = inv_dx*( field(i+1,j,  k  )-field(i,  j,  k  ))
            dijk       = dip1
            !
            dip1jp1    = inv_dx*( field(i+1,j+1,k  )-field(i,  j+1,k  ))
            djp1       = dip1jp1
            !
            dip1kp1    = inv_dx*( field(i+1,j,  k+1)-field(i,  j,  k+1))
            dkp1       = dip1kp1
            !
            dip1jp1kp1 = inv_dx*( field(i+1,j+1,k+1)-field(i,  j+1,k+1))
            djp1kp1    = djp1kp1
            ! Value of fx does not matter now because
         else if (fx.gt. 0.5 .and. i.lt.nx-1) then
            dip1       = inv_dx*( field(i+2,j,  k  )-field(i+1,j,  k  ))
            dijk       = inv_dx*( field(i+1,j,  k  )-field(i,  j,  k  ))
            !
            dip1jp1    = inv_dx*( field(i+2,j+1,k  )-field(i+1,j+1,k  ))
            djp1       = inv_dx*( field(i+1,j+1,k  )-field(i,  j+1,k  ))
            !
            dip1kp1    = inv_dx*( field(i+2,j,  k+1)-field(i+1,j,  k+1))
            dkp1       = inv_dx*( field(i+1,j,  k+1)-field(i,  j,  k+1))
            !
            dip1jp1kp1 = inv_dx*( field(i+2,j+1,k+1)-field(i+1,j+1,k+1))
            djp1kp1    = inv_dx*( field(i+1,j+1,k+1)-field(i,  j+1,k+1))
            !
            ! Stencil shifts in x-direction because of derivs:
            fx = fx - 0.5d0
         else if (fx.gt. 0.5 .and. i.eq.nx-1) then
            ! At the boundary and so we could use back/forwards differncing,
            ! but instead just juse derivative a half-step from boundary:
            dijk       = inv_dx*( field(i+1,j,  k  )-field(i,  j,  k  ))
            dip1       = dijk
            !
            djp1       = inv_dx*( field(i+1,j+1,k  )-field(i,  j+1,k  ))
            dip1jp1    = djp1
            !
            dkp1       = inv_dx*( field(i+1,j,  k+1)-field(i,  j,  k+1))
            dip1kp1    = dip1kp1
            !
            djp1kp1    = inv_dx*( field(i+1,j+1,k+1)-field(i,  j+1,k+1))
            dip1jp1kp1 = dip1kp1
            ! Value of fx does not matter now because
         else 
            write(*,*) 'interpderiv: something very wrong w/ x-deriv'
            write(*,*) 'interpderiv:  i/j/k:',i, j, k
            write(*,*) 'interpderiv: fx/y/z:',fx,fy,fz
            write(*,*) 'interpderiv: nx/y/z:',nx,ny,nz
            somethingwrong = .true.
         end if
      !...................................................
      else if (direction.eq.2) then
         if (fy.le.0.5 .and. j.gt.1 .and. j.lt.ny) then
            djp1       = inv_dx*( field(i,  j+1,k  )-field(i,  j,  k  ))
            dijk       = inv_dx*( field(i,  j,  k  )-field(i,  j-1,k  ))
            !
            dip1jp1    = inv_dx*( field(i+1,j+1,k  )-field(i+1,j,  k  ))
            dip1       = inv_dx*( field(i+1,j,  k  )-field(i+1,j-1,k  ))
            !
            djp1kp1    = inv_dx*( field(i,  j+1,k+1)-field(i,  j,  k+1))
            dkp1       = inv_dx*( field(i,  j,  k+1)-field(i,  j-1,k+1))
            !
            dip1jp1kp1 = inv_dx*( field(i+1,j+1,k+1)-field(i+1,j,  k+1))
            dip1kp1    = inv_dx*( field(i+1,j,  k+1)-field(i+1,j-1,k+1))
            !
            ! Stencil shifts in y-direction because of derivs:
            fy = fy + 0.5d0
         else if (fy.le.0.5 .and. j.eq.1) then
            ! At the boundary and so we could use back/forwards differncing,
            ! but instead just juse derivative a half-step from boundary:
            djp1       = inv_dx*( field(i,  j+1,k  )-field(i,  j,  k  ))
            dijk       = djp1
            !
            dip1jp1    = inv_dx*( field(i+1,j+1,k  )-field(i+1,j,  k  ))
            dip1       = dip1jp1
            !
            djp1kp1    = inv_dx*( field(i,  j+1,k+1)-field(i,  j,  k+1))
            dkp1       = djp1kp1
            !
            dip1jp1kp1 = inv_dx*( field(i+1,j+1,k+1)-field(i+1,j,  k+1))
            dip1kp1    = dip1jp1kp1
            ! Value of fy does not matter now because
         else if (fy.gt. 0.5 .and. j.lt.ny-1) then
            djp1       = inv_dx*( field(i,  j+2,k  )-field(i,  j+1,k  ))
            dijk       = inv_dx*( field(i,  j+1,k  )-field(i,  j,  k  ))
            !
            dip1jp1    = inv_dx*( field(i+1,j+2,k  )-field(i+1,j+1,k  ))
            dip1       = inv_dx*( field(i+1,j+1,k  )-field(i+1,j,  k  ))
            !
            djp1kp1    = inv_dx*( field(i,  j+2,k+1)-field(i,  j+1,k+1))
            dkp1       = inv_dx*( field(i,  j+1,k+1)-field(i,  j,  k+1))
            !
            dip1jp1kp1 = inv_dx*( field(i+1,j+2,k+1)-field(i+1,j+1,k+1))
            dip1kp1    = inv_dx*( field(i+1,j+1,k+1)-field(i+1,j,  k+1))
            !
            ! Stencil shifts in y-direction because of derivs:
            fy = fy - 0.5d0
         else if (fy.gt. 0.5 .and. j.eq.ny-1) then
            ! At the boundary and so we could use back/forwards differncing,
            ! but instead just juse derivative a half-step from boundary:
            dijk       = inv_dx*( field(i,  j+1,k  )-field(i,  j,  k  ))
            djp1       = dijk
            !
            dip1       = inv_dx*( field(i+1,j+1,k  )-field(i+1,j,  k  ))
            dip1jp1    = dip1jp1
            !
            dkp1       = inv_dx*( field(i,  j+1,k+1)-field(i,  j,  k+1))
            djp1kp1    = dkp1
            !
            dip1kp1    = inv_dx*( field(i+1,j+1,k+1)-field(i+1,j,  k+1))
            dip1jp1kp1 = dip1kp1
            ! Value of fy does not matter now because
         else 
            write(*,*) 'interpderiv: something very wrong w/ y-deriv'
            write(*,*) 'interpderiv:  i/j/k:',i, j, k
            write(*,*) 'interpderiv: fx/y/z:',fx,fy,fz
            write(*,*) 'interpderiv: nx/y/z:',nx,ny,nz
            somethingwrong = .true.
         end if
      !...................................................
      else if (direction.eq.3) then
         if (fz.le.0.5 .and. k.gt.1 .and. k.lt.nz) then
            dkp1       = inv_dx*( field(i,  j,  k+1)-field(i,  j,  k  ))
            dijk       = inv_dx*( field(i,  j,  k  )-field(i,  j,  k-1))
            !
            dip1kp1    = inv_dx*( field(i+1,j,  k+1)-field(i+1,j,  k  ))
            dip1       = inv_dx*( field(i+1,j,  k  )-field(i+1,j,  k-1))
            !
            djp1kp1    = inv_dx*( field(i,  j+1,k+1)-field(i,  j+1,k  ))
            djp1       = inv_dx*( field(i,  j+1,k  )-field(i,  j+1,k-1))
            !
            dip1jp1kp1 = inv_dx*( field(i+1,j+1,k+1)-field(i+1,j+1,k  ))
            dip1jp1    = inv_dx*( field(i+1,j+1,k  )-field(i+1,j+1,k-1))
            !
            ! Stencil shifts in x-direction because of derivs:
            fz = fz + 0.5d0
         else if (fz.le.0.5 .and. k.eq.1) then
            ! At the boundary and so we could use back/forwards differncing,
            ! but instead just juse derivative a half-step from boundary:
            dkp1       = inv_dx*( field(i,  j,  k+1)-field(i,  j,  k  ))
            dijk       = dkp1
            !
            dip1kp1    = inv_dx*( field(i+1,j,  k+1)-field(i+1,j,  k  ))
            dip1       = dip1kp1
            !
            djp1kp1    = inv_dx*( field(i,  j+1,k+1)-field(i,  j+1,k  ))
            djp1       = djp1kp1
            !
            dip1jp1kp1 = inv_dx*( field(i+1,j+1,k+1)-field(i+1,j+1,k  ))
            dip1jp1    = dip1jp1kp1
            ! Value of fz does not matter now because
         else if (fz.gt.0.5 .and. k.lt.nz-1) then
            dkp1       = inv_dx*( field(i,  j,  k+2)-field(i,  j,  k+1))
            dijk       = inv_dx*( field(i,  j,  k+1)-field(i,  j,  k  ))
            !
            dip1kp1    = inv_dx*( field(i+1,j,  k+2)-field(i+1,j,  k+1))
            dip1       = inv_dx*( field(i+1,j,  k+1)-field(i+1,j,  k  ))
            !
            djp1kp1    = inv_dx*( field(i,  j+1,k+2)-field(i,  j+1,k+1))
            djp1       = inv_dx*( field(i,  j+1,k+1)-field(i,  j+1,k  ))
            !
            dip1jp1kp1 = inv_dx*( field(i+1,j+1,k+2)-field(i+1,j+1,k+1))
            dip1jp1    = inv_dx*( field(i+1,j+1,k+1)-field(i+1,j+1,k  ))
            !
            ! Stencil shifts in z-direction because of derivs:
            fz = fz - 0.5d0
         else if (fz.gt.0.5 .and. k.eq.nz-1) then
            ! At the boundary and so we could use back/forwards differncing,
            ! but instead just juse derivative a half-step from boundary:
            dijk       = inv_dx*( field(i,  j,  k+1)-field(i,  j,  k  ))
            dkp1       = dijk
            !
            dip1       = inv_dx*( field(i+1,j,  k+1)-field(i+1,j,  k  ))
            dip1kp1    = dip1
            !
            djp1       = inv_dx*( field(i,  j+1,k+1)-field(i,  j+1,k  ))
            djp1kp1    = djp1
            !
            dip1jp1    = inv_dx*( field(i+1,j+1,k+1)-field(i+1,j+1,k  ))
            dip1jp1kp1 = dip1jp1
            ! Value of fz does not matter now because
         else 
            write(*,*) 'interpderiv: something very wrong w/ z-deriv'
            write(*,*) 'interpderiv:  i/j/k:',i, j, k
            write(*,*) 'interpderiv: fx/y/z:',fx,fy,fz
            write(*,*) 'interpderiv: nx/y/z:',nx,ny,nz
            somethingwrong = .true.
         end if
      end if
      !...................................................
      if (somethingwrong) then
      interpderiv = 0.d0
      else
      interpderiv = 
     *      + (1.d0-fx) * (1.d0-fy) * (1.d0-fz) * dijk
      !
     *      + (1.d0-fx) *       fy  *       fz  * djp1kp1
     *      +       fx  * (1.d0-fy) *       fz  * dip1kp1
     *      +       fx  *       fy  * (1.d0-fz) * dip1jp1
      !
     *      +       fx  * (1.d0-fy) * (1.d0-fz) * dip1
     *      + (1.d0-fx) *       fy  * (1.d0-fz) * djp1
     *      + (1.d0-fx) * (1.d0-fy) *       fz  * dkp1
      !
     *      +       fx  *       fy  *       fz  * dip1jp1kp1
      end if

      return
      end    ! END: interpderiv


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interp 4point Lagrange                                                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function interp4( field, fx, fy, fz,
     *                                  i,j,k, nx,ny,nz)
      implicit     none
      integer      i,j,k, nx, ny, nz, li,lm,lj,lk
      real(kind=8) field(nx,ny,nz),
     *             fx, fy, fz, F2int(4,4), F1int(4),half,sixth

      logical     ltrace
      parameter ( ltrace = .false. )

      if (i .ge. nx-1 .or. j .ge. ny-1 .or. k .ge. nz-1 .or.
     *    i .le. 1    .or. j .le. 1    .or. k .le. 1         ) return

      !
      half = 0.5d0
      sixth = .1666666666666666d0
      
      F2int = 0.0d0
      
c intepolate in x      
      do li=1,4
        do lm =1,4
	
	lj = j-2+li
	lk = k-2+lm
	
      F2int(li,lm)= -field(i-1,lj,lk) * fx*(fx-1.0)*(fx-2.0)      *sixth 
     *              +field(i,lj,lk)   * (fx+1.0)*(fx-1.0)*(fx-2.0)*half
     *              -field(i+1,lj,lk) * (fx+1.0)*(fx)*(fx-2.0)    *half
     *              +field(i+2,lj,lk) * (fx+1.0)*(fx)*(fx-1.0)    *sixth
        end do
       end do	
      
 
c interpolate in y      
      F1int = 0.0
      
      do lm = 1,4 
       F1int(lm) = - F2int(1,lm) * fy*(fy-1.0)*(fy-2.0)       *sixth
     *             + F2int(2,lm) * (fy+1.0)*(fy-1.0)*(fy-2.0) *half
     *             - F2int(3,lm) * (fy+1.0)*(fy)*(fy-2.0)     *half
     *             + F2int(4,lm) * (fy+1.0)*(fy)*(fy-1.0)     *sixth 
      end do

c interpolate in z      
      interp4 = 0.0
      
      do lm = 1,4 
       interp4   = - F1int(1) * fz*(fz-1.0)*(fz-2.0)       *sixth
     *             + F1int(2) * (fz+1.0)*(fz-1.0)*(fz-2.0) *half
     *             - F1int(3) * (fz+1.0)*(fz)*(fz-2.0)     *half
     *             + F1int(4) * (fz+1.0)*(fz)*(fz-1.0)     *sixth 
      end do


      return
      end    ! END: interp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interp2d:                                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function interp2d( field, fx, fy,
     *                                  i,j,nx,ny)
      implicit     none
      integer      i,j,nx, ny
      real(kind=8) field(nx,ny),
     *             fx, fy

      logical     ltrace
      parameter ( ltrace = .false. )

      if (i .eq. nx .or. j .eq. ny) then
         interp2d = 
     *                                field(i,  j  )
      else if (i .gt. nx .or. j .gt. ny .or.
     *    i .lt. 1    .or. j .lt. 1         ) then
         return
      else
         interp2d = 
     *      + (1.d0-fx) * (1.d0-fy) * field(i,  j  )
     *      +       fx  *       fy  * field(i+1,j+1)
     *      +       fx  * (1.d0-fy) * field(i+1,j  )
     *      + (1.d0-fx) *       fy  * field(i,  j+1)
      end if

      return
      end    ! END: interp2d
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mask_return_closestdist:  compute the radius of a given point from a      cc
cc                            black hole                                      cc
cc                            NB: if no black hole currently being excised,   cc
cc                                then returns coordinate distance from origincc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mask_return_closestdist(x,y,z,dist)

      implicit    none
      integer     l
      real(kind=8) x,y,z,dist
      real(kind=8) centerx,centery,centerz,tmpdist
      include     'mask.inc'
      include     'largesmall.inc'
      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'mask_return_closestdist x,y,z = ',x,y,z
      end if

      dist = LARGENUMBER
      do l=1,max_num_masks
        if (bh_true(l)) then
          centerx = 0.5*(mask_bbox(1,l)+mask_bbox(2,l))
          centery = 0.5*(mask_bbox(3,l)+mask_bbox(4,l))
          centerz = 0.5*(mask_bbox(5,l)+mask_bbox(6,l))
          tmpdist = sqrt( (centerx-x)**2 
     *                  + (centery-y)**2 
     *                  + (centerz-z)**2)
          if ( tmpdist .lt. dist ) dist = tmpdist
        end if  
      end do

      if (dist.eq.LARGENUMBER) dist = sqrt(x**2+y**2+z**2)

      if (ltrace) write(*,*) 'mask_return_closestdist: Done.'

      return
      end    ! END: mask_return_closestdist

      !-----------------------------------------------------------------
      !
      !  interp_from_parentC
      !
      !       This routine interpolates data from a parent array, up,
      !       to a child array, uc.  The arrays overlap over some 
      !       region, and the points within the region R
      !
      !            R := [mini, maxi] x [minj, maxj] x [mink, maxk]
      !
      !       of the array uc need to be filled from the up array.
      !       Thus the outcome of this routine is to modify the R of uc.
      !
      !       This routine can be used for interpolation of arbitrary
      !       order.  It is currently set up for cubic and quintic
      !       interpolation.  
      !
      !       The only allowed interpolation stencil is a centered
      !       stencil.  For the stencil below, the points "o" exist
      !       on the parent grid, and + is a point on the child grid.
      !       As the stencil is centered, the point + is in [x_i, x_{i+1}].
      !       
      !           o         o         o    +    o         o         o
      !         i + sl                i        i+1                 i + su
      !                          
      !       The end points of the stencil are labeld in reference
      !       to the point x_i.  The first point has the index (i+sl)
      !       and the last point is (i+su).  Clearly, sl < 0 and su > 0.
      !       The number of points in the stencil is ns.  Finally, to 
      !       simplify some loop variables, we define 
      !            s0 = 1 - sl.
      !
      !       The order of interpolation is determined (partially) by
      !       the stencil size.  For cubic interpolation, ns = 4, 
      !       sl = -1, and su = 2.  For quintic interpolation, ns = 6
      !       sl = -2, and su = 3.  Higher order interpolation stencils
      !       can be contructed by setting ns, sl, su, and s0.
      !
      !
      !       IMPORTANT NOTES
      !       -------------------------------------------------------
      !
      !       (1) The interpolation order is set by the CPP variable
      !           INTERP_ORDER.  Currently INTERP_ORDER can be 3 or 5.
      ! 
      !       (2) This routine is written only for r == 2.  This is
      !           the regridding ratio.  I did not attempt to generalize
      !           the routine for different ratios.
      !
      !-----------------------------------------------------------------
      subroutine interp_from_parentC(up,      uc,      chr,     mask,
     &                               nxp,     nyp,     nzp,
     &                               mini,    minj,    mink,
     &                               maxi,    maxj,    maxk,
     &                               bi,      bj,      bk,
     &                               nxc,     nyc,     nzc,     r,
     &                               gw)

        implicit none
        CCTK_INT                           :: nxp, nyp, nzp
        CCTK_INT                           :: nxc, nyc, nzc
        CCTK_INT                           :: bi,  bj,  bk
        CCTK_INT                           :: r,   gw
        CCTK_INT                           :: mini,minj,mink
        CCTK_INT                           :: maxi,maxj,maxk
        CCTK_REAL, dimension(nxp,nyp,nzp)  :: up
        CCTK_REAL, dimension(nxc,nyc,nzc)  :: uc
        CCTK_REAL, dimension(nxc,nyc,nzc)  :: chr
        CCTK_REAL, dimension(nxc,nyc,nzc)  :: mask

        ! Set interpolation order using CPP.

!>>>>>>>>>>>>>>>>>>  SET INTERPOLATION ORDER HERE
! INTERP_ORDER can be either 3 or 5
#define INTERP_ORDER 5
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#ifdef INTERP_ORDER
#if INTERP_ORDER == 3
#define INTERP interp3
        ! third order interpolation
        integer, parameter                 :: sl = -1
        integer, parameter                 :: su =  2
        integer, parameter                 :: ns =  4
        integer, parameter                 :: s0 =  2
        CCTK_REAL                          :: interp3
        external                           :: interp3
#elif INTERP_ORDER == 5
#define INTERP interp5
        ! fifth order interpolation
        integer, parameter                 :: sl = -2
        integer, parameter                 :: su =  3
        integer, parameter                 :: ns =  6
        integer, parameter                 :: s0 =  3
        CCTK_REAL                          :: interp5
        external                           :: interp5
#endif
#endif

        ! local vars
        CCTK_INT                           :: ic, jc, kc, ip, jp, kp
        CCTK_INT                           :: xic, xjc, xkc
        CCTK_INT                           :: xip, xjp, xkp
        CCTK_INT                           :: ii,  jj,  kk
        CCTK_REAL, dimension(ns)  :: utemp, utemp2, utemp3
        include 'largesmall.inc'
        include 'chr.inc'

        logical, parameter                 :: ltrace = .false.

        if (ltrace) write(*,*)'Begin interp_from_parentC'

        if (maxi .lt. mini .or. maxj .lt. minj .or. maxk .lt. mink) then
          write(*,*)'interp_from_parentC:  max < min, returning'
          write(*,*)'interp_from_parentC:  mini, maxi ',mini, maxi
          write(*,*)'interp_from_parentC:  minj, maxj ',minj, maxj
          write(*,*)'interp_from_parentC:  mink, maxk ',mink, maxk
          return
        end if

        ! Is there sufficient data to interpolate only symmetrically?
        kp = (maxk - mink)/r + 1 + bk
        if (kp .gt. nzp) then
          write(*,*) 'interp_from_parentC: Not sufficient for large z'
          write(*,*) '  kp   = ',kp
          write(*,*) '  maxk = ',maxk
        end if
        jp = (maxj - minj)/r + 1 + bj
        if (jp .gt. nyp) then
          write(*,*) 'interp_from_parentC: Not sufficient for large y'
          write(*,*) '  jp   = ',jp
          write(*,*) '  maxj = ',maxj
        end if
        ip = (maxi - mini)/r + 1 + bi
        if (ip .gt. nxp) then
          write(*,*) 'interp_from_parentC: Not sufficient for large x'
          write(*,*) '  ip   = ',ip
          write(*,*) '  maxi = ',maxi
        end if

        ! Interpolate in x, aligned in y and z
        do kc = mink, maxk, r
          kp = (kc - mink)/r + 1 + bk
          do jc = minj, maxj, r
            jp = (jc - minj)/r + 1 + bj
            ! Copy all points on refined grid that coincide with points
            ! on the course grid.
            do ic = mini, maxi
              ip = (ic - mini)/r + 1 + bi
              xic = r*(ip - 1 - bi) + mini
              if (xic .eq. ic) then
                ! The points on the child grid and parent grid align exactly
                uc(ic,jc,kc) = up(ip,jp,kp)
!               write(*,*) '%% point in common ',uc(ic,jc,kc)
              else
                ! The points on the child grid and parent grid 
                ! do not align in x.
                if (ip+sl .ge. 1 .and. ip+su .le. nxp) then
                  do ii = sl, su
                    utemp(ii+s0) = up(ip+ii,jp,kp)
                  end do
                  uc(ic,jc,kc) = INTERP(utemp)
!                 write(*,*) '%% interpolate from parent ',uc(ic,jc,kc)
!                 write(*,*) '%% ',utemp
                else
                  write(*,*) 'interp_from_parentC: ',
     &                      'I probably need a point here.'
                end if
              end if
            end do
          end do
        end do

        ! Interpolate in y, aligned in z
        do kc = mink, maxk, r
          kp = (kc - mink)/r + 1 + bk
          do jc = minj+1, maxj-1, r  ! go only along rows not aligned in y
            jp = (jc - minj)/r + 1 + bj
            do ic = mini, maxi
              ip = (ic - mini)/r + 1 + bi
              xic = r*(ip - 1 - bi) + mini
              if (ic .eq. xic) then
                ! aligned in x
                if (jp+sl .ge. 1 .and. jp+su .le. nyp) then
                  do jj = sl, su
                    utemp(jj+s0) = up(ip,jp+jj,kp)
                  end do
                  uc(ic,jc,kc) = INTERP(utemp)
                else
                  write(*,*)'trying to interpolate in y, and missing'
                  write(*,*)'some points.'
                end if
              else
                ! not aligned in x...need to build the stencil
                do jj = sl, su
                  xjp = jp + jj
                  xjc = r*(xjp - 1 - bj) + minj
                  if (xjc .ge. minj .and. xjc .le. maxj) then
                    ! This point was already computed, so we can use it
                    utemp(s0+jj) = uc(ic,xjc,kc)
                  else
                    if (ip+sl .ge. 1 .and. ip+su .le. nxp) then
                      ! We can build the point needed.
                      do ii = sl, su
                        utemp2(s0+ii) = up(ip+ii,xjp,kp)
                      end do
                      utemp(s0+jj) = INTERP(utemp2)
                    else
                      ! Not enough info to build this point.  Jump out.
                      write(*,*)'interp_from_parentC:',
     &                          'Can not interpolate in x to get points'
                      utemp(s0+jj) = LARGENUMBER
                    end if
                  end if
                end do
                uc(ic,jc,kc) = INTERP(utemp)
              end if
            end do
          end do
        end do

        ! Interpolate in z, 
        do kc = mink+1, maxk-1, r
          kp = (kc - mink)/r + 1 + bk
          do jc = minj, maxj
            jp = (jc - minj)/r + 1 + bj
            xjc = r*(jp - 1 - bj) + minj
            do ic = mini, maxi
              ip = (ic - mini)/r + 1 + bi
              xic = r*(ip - 1 - bi) + mini

              do kk = sl, su
                xkp = kp + kk
                xkc = r*(xkp - 1 - bk) + mink
                if (xkc .ge. mink .and. xkc .le. maxk) then
                  ! This is a point that has already been calculated
                  utemp(s0+kk) = uc(ic,jc,xkc)
                else 
                  ! We need to build the interpolated value.
                  if (xic .eq. ic .and. xjc .eq. jc) then
                    ! Here we are lucky, and x and y are aligned.
                    utemp(s0+kk) = up(ip,jp,xkp)
                  else if (xic .eq. ic .and. xjc .ne. jc) then
                    ! Aligned in x, interpolate in y
                    do jj = sl, su
                      xjp = jp + jj
                      utemp2(s0+jj) = up(ip,xjp,xkp)
                    end do 
                    utemp(s0+kk) = INTERP(utemp2)
                  else if (xic .ne. ic .and. xjc .eq. jc) then
                    ! Aligned in y, interpolate in x
                    do ii = sl, su
                      xip = ip + ii
                      utemp2(s0+ii) = up(xip,jp,xkp)
                    end do 
                    utemp(s0+kk) = INTERP(utemp2)
                  else if (xic .ne. ic .and. xjc .ne. jc) then
                    ! No alignment, interpolate in both x and y
                    do ii = sl, su
                      xip = ip + ii
                      do jj = sl, su
                        xjp = jp + jj
                        utemp3(s0+jj) = up(xip,xjp,xkp)
                      end do 
                      utemp2(s0+ii) = INTERP(utemp3)
                    end do 
                    utemp(s0+kk) = INTERP(utemp2)
                  else
                    write(*,*)'interp_from_parentC: ',
     &                        'Should not be here'
                  end if 
            
                end if
              end do
              uc(ic,jc,kc) = INTERP(utemp)
            end do
          end do
        end do

        if (ltrace) write(*,*)'End interp_from_parentC'

        return
      end

      !-----------------------------------------------------------------
      !
      !   interp3
      !
      !     cubic interpolation for centered interpolation stencil.  The
      !     "o" points are used to construct the function at the point "+"
      !      
      !            o         o    +    o         o 
      !                      i        i+1         
      !
      !-----------------------------------------------------------------
      CCTK_REAL function interp3(f)
        implicit none
        CCTK_REAL, dimension(4) :: f

        interp3 = 0.0625d0 * (9.0d0*(f(2) + f(3)) - (f(1) + f(4)))
        return 
      end

      !-----------------------------------------------------------------
      !
      !   interp5
      !
      !      quintic interpolation for a centered stencil.  The "o" points
      !      are used to construct the function at the point "+"
      !
      !           o         o         o    +    o         o         o
      !         i + sl                i        i+1                 i + su
      !
      !-----------------------------------------------------------------
      CCTK_REAL function interp5(f)
        implicit none
        CCTK_REAL, dimension(6) :: f

        interp5 = 0.00390625d0 * ( 150.0d0*(f(3) + f(4)) 
     &                 - 25.0d0*(f(2) + f(5)) + 3.0d0*(f(1)+f(6)) )

        return 
      end

      !-----------------------------------------------------------------
      !
      !
      !
      !-----------------------------------------------------------------
      subroutine interp_compare(fold, fnew, mini, minj, mink,
     &                          maxi, maxj, maxk,
     &                          nx, ny, nz)
        implicit none
        integer  :: nx, ny, nz
        integer  :: mini, minj, mink, maxi, maxj, maxk 
        real(kind=8), dimension(nx,ny,nz) :: fold, fnew

        ! local vars
        integer      :: i, j, k
        real(kind=8) :: df
        logical, parameter :: ltrace = .false.

        if (ltrace) write(*,*)'Begin interp_compare'

        if (mini .lt. 1 .or. mini .gt. nx) then
          write(*,*) 'interp_compare:  mini out of bounds: ',mini
        end if
        if (maxi .lt. 1 .or. maxi .gt. nx) then
          write(*,*) 'interp_compare:  maxi out of bounds: ',maxi
        end if
        if (minj .lt. 1 .or. minj .gt. ny) then
          write(*,*) 'interp_compare:  minj out of bounds: ',minj
        end if
        if (maxj .lt. 1 .or. maxj .gt. ny) then
          write(*,*) 'interp_compare:  maxj out of bounds: ',maxj
        end if
        if (mink .lt. 1 .or. mink .gt. nz) then
          write(*,*) 'interp_compare:  mink out of bounds: ',mink
        end if
        if (maxk .lt. 1 .or. maxk .gt. nz) then
          write(*,*) 'interp_compare:  maxk out of bounds: ',maxk
        end if

        if (ltrace) write(*,*)'Begin interp_compare LOOP'

        do k = mink, maxk
        do j = minj, maxj
        do i = mini, maxi
          df = abs(fold(i,j,k) - fnew(i,j,k))
          if (df .gt. 1.0d-12) then
            write(*,*)'<<<<<<'
            write(*,*)'   Problem at ',i,j,k
            write(*,*)'   fold       ',fold(i,j,k)
            write(*,*)'   fnew       ',fnew(i,j,k)
            stop
          end if
        end do
        end do
        end do

        if (ltrace) write(*,*)'End interp_compare'
      
        return
      end subroutine

      !-----------------------------------------------------------------
      !
      !
      !
      !-----------------------------------------------------------------
