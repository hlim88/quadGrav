#include "cctk.h"

!-------------------------------------------------------------------------
!   fdata1d       - data from sdf file.
!   coords        - coord array from sdf file
!   dnr           - size of fdata1d
!   f             - function used in the code
!   nx, ny, nz    - size of f
!   xf            - 1d coord x-array 
!   yf            - 1d coord y-array 
!   zf            - 1d coord z-array 
!
!-------------------------------------------------------------------------
      subroutine read_sdf_file1d                                    &
                           (filename, f, xf, yf, zf, nx, ny, nz,    &
                            level, fdata1d, coords1d, dnr,          &
                            interp, xc, yc, zc)
      implicit none

      CCTK_INT               nx, ny, nz, level
      CCTK_INT               dnr
      CCTK_INT               interp
      CCTK_REAL              f(nx,ny,nz)
      CCTK_REAL              xf(nx), yf(ny), zf(nz)
      CCTK_REAL              fdata1d(dnr)
      CCTK_REAL              coords1d(dnr)
      character*(*)          filename
      CCTK_REAL           ::   xc,  yc,  zc


      CCTK_INT               rank, shp(1), rc
      CCTK_INT               i, j, k, istat
      logical, parameter  :: FLT = .false.
      logical, parameter  :: ltrace = .false.

      CCTK_REAL              time
      character*64           cnames

      CCTK_INT            :: order
      CCTK_REAL           :: ddr, drmin, drmax
      CCTK_INT            :: m,   mb,  mc
      CCTK_REAL           :: xp,  yp,  zp,  rp
      CCTK_REAL           :: cr(10)

      ! functions
      CCTK_INT   gft_read_rank, gft_read_shape, gft_read_full
      CCTK_REAL, external ::  myl2norm1d
      logical    double_equal

      interface
         subroutine int4coef(c, xi, x)
           implicit none
           CCTK_REAL c(5), xi, x(*)
         end subroutine int4coef
      end interface

#ifdef AIX
      rc = gft_read_rank(trim(filename)//CHAR(0),%VAL(level), rank)
#else
      rc = gft_read_rank(filename, level, rank)
#endif
      if (FLT .or. ltrace) then
         write(*,*) '...Attempt to read ',filename
      end if
      if (rc .ne. 1) then
         write(*,*) 'Unable to read as SDF file ',filename
         stop
      end if
      if ( rank .ne. 1 ) then
         write(*,*) 'Trying to read an sdf which is NOT 1D' 
         stop
      endif
      if (FLT) then
         write(*,*) '...rank = ', rank
      end if
#ifdef AIX
      rc = gft_read_shape(trim(filename)//CHAR(0), %VAL(level), shp)
#else
      rc = gft_read_shape(filename, level, shp)
#endif
      if (rc .ne. 1) then
         write(*,*) 'Unable to read shape for data in ',filename
         stop
      end if
      if (FLT) then
         write(*,*) '...shp = ', shp(1)
      end if
      if ( shp(1) .gt. dnr ) then
         write(*,*) 'fdata1d is not large enough to read in this file'
         write(*,*) 'file = ', filename
         write(*,*) 'size of fdata1d = ',dnr
         write(*,*) 'size of data in file = ',shp(1)
         stop
      end if

      if (FLT) write(*,*) '...calling gft_read_full'
#ifdef AIX
      rc = gft_read_full(trim(filename)//CHAR(0), %VAL(level), shp,cnames, rank, &
#else
      rc = gft_read_full(filename, level, shp, cnames, rank, &
#endif
                         time, coords1d, fdata1d)
      if (FLT) write(*,*) '...end gft_read_full'

      if (rc .ne. 1) then
         write(*,*) 'Unable to read data in ',filename
         stop
      end if

      if (double_equal(myl2norm1d(fdata1d,shp(1)),0.d0))  then
         if (ltrace) write(*,*) '...initializing to zero: ',filename
         call load_scal3d( f, 0.d0, nx, ny, nz)
         goto 100
      end if

      ! Extract local data from the global data set

      if (interp .eq. 0) then
         drmin = coords1d(1)
         ddr   = coords1d(2) - coords1d(1)
         !
         ! load data in array
         if (FLT) write(*,*) '...straight loading data'
         do k = 1, nz
            zp  =  zf(k) - zc 
            do j = 1, ny
               yp  =  yf(j) - yc 
               do i = 1, nx
                  xp  =  xf(i) - xc 
                  rp  = sqrt(xp**2 + yp**2 + zp**2)
                  m   = int((rp-drmin)/ddr) + 1
                  f(i,j,k) = fdata1d(m)
               end do
            end do
         end do
      else if (interp .eq. 1) then
         ! interpolation
         if (FLT) write(*,*) '...interpolating'
         dnr   = shp(1)
         ddr   = coords1d(2) - coords1d(1)
         drmin = coords1d(1)
         drmax = coords1d(dnr)
         order = 4
         if (FLT) then
            write(*,*) 'dnr   = ',dnr
            write(*,*) 'drmin = ',drmin
            write(*,*) 'drmax = ',drmax
            write(*,*) 'ddr   = ',ddr
         end if

         do k = 1, nz
            zp = zf(k) - zc
            do j = 1, ny
               yp = yf(j) - yc
               do i = 1, nx
                  xp = xf(i) - xc
                  rp  = sqrt(xp**2 + yp**2 + zp**2)
!                  if (rp .lt. drmin) then
!                      write(*,*) 'Problem: ',xp,yp,zp
!                      stop
!                  end if
                  ! estimate the indices for the point in fdata closest to 
                  ! (xp,yp,zp)
                  m  = int(rp/ddr) + 1 

                  ! find the first point in the interpolation stencil.
                  ! we are careful about boundaries.

                  mb = min(max(m-2,1),dnr-4)
                  call int4coef(cr, rp, coords1d(mb))

                  f(i,j,k) = 0.0d0
                  do mc = 1, 5
                     f(i,j,k) = f(i,j,k)  + cr(mc)*fdata1d(mb+mc-1)
                  end do

               end do
            end do
         end do

      else
         write(*,*)'...interp has illegal value in read_sdf_file1d'
         write(*,*)'...interp = ',interp
         stop
      end if
      
100  continue

      return
   end subroutine read_sdf_file1d

!-------------------------------------------------------------------------------
      real*8 function myl2norm1d(u, nr)
!-------------------------------------------------------------------------------
      implicit none
      integer          nr
      real(kind=8)   u(nr)
      integer        i


      myl2norm1d = 0.d0

      if ((nr.le.0)) then
         write(*,*) 'myl2norm1d: ERROR'
      end if

      do i = 1, nr
         myl2norm1d = myl2norm1d + u(i)**2
      end do

      myl2norm1d = sqrt(myl2norm1d/nr)
      return
      end


