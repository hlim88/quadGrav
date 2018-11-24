!-------------------------------------------------------------------------
!
!  $Id$
!
!  This routine reads an SDF file.  In a small attempt to keep this
!  generic, this routine reads a single file for which the name, storage,
!  etc., etc., must be passed in.  This file should be called from a 
!  routine in the equation thorn.
!
!  WARNING! NOTE! ACHTUNG! Pay attention here:  Because the SDF libraries
!  are not installed everywhere, this routine is wrapped by a CPP macro,
!  and is only activated *IF* you specifically set the macro variable
!  HAVE_SDF_LIB before compiling.
!
!-------------------------------------------------------------------------

#include "cctk.h"

!-------------------------------------------------------------------------
!
!  This routine simply returns the shape of the sdf file filename.
!  This is done so that memory can be allocated to read the file.
!  This should be called before read_sdf_file.
!
!-------------------------------------------------------------------------
      subroutine sdf_file_mem(filename, level, shp)
      implicit none

      character*(*)   filename
      CCTK_INT         ::  shp(3), level

      CCTK_INT         ::  rc, rank
      logical, parameter                        :: FLT = .false.
      logical, parameter                        :: ltrace = .false.

      ! functions
      CCTK_INT   gft_read_rank, gft_read_shape

      level = 1
#ifdef AIX
      rc = gft_read_rank(trim(filename)//CHAR(0), %VAL(level), rank) 
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
         write(*,*) '...shp = ', shp(1),shp(2), shp(3)
      end if

      return
   end subroutine sdf_file_mem

!-------------------------------------------------------------------------
!   fdata         - data from sdf file.
!   coords        - coord array from sdf file
!   dnx, dny, dnz - size of fdata
!   f             - function used in the code
!   nx, ny, nz    - size of f
!   xf            - 1d coord x-array 
!   yf            - 1d coord y-array 
!   zf            - 1d coord z-array 
!
!-------------------------------------------------------------------------
      subroutine read_sdf_file(filename, f, xf, yf, zf, nx, ny, nz, &
                            level, fdata, coords, dnx, dny, dnz,    &
                            interp, cctk_llb, cctk_gshp)
      implicit none
      CCTK_INT        nx, ny, nz, level
      CCTK_INT        dnx, dny, dnz
      CCTK_INT        interp, cctk_llb(3), cctk_gshp(3)
      CCTK_REAL       f(nx,ny,nz)
      CCTK_REAL       xf(nx), yf(ny), zf(nz)
      CCTK_REAL       fdata(dnx,dny,dnz)
      CCTK_REAL       coords(dnx+dny+dnz)
      character*(*)   filename


      CCTK_INT        rank, shp(3), rc, ilb, jlb, klb
      CCTK_INT        i, j, k, istat
      logical, parameter                        :: FLT = .false.
      logical, parameter                        :: ltrace = .false.
      CCTK_REAL       time
      character*64    cnames
      CCTK_INT   ::   order
      CCTK_REAL  :: ddx, ddy, ddz, dxmin, dymin, dzmin
      CCTK_REAL  :: x1, y1, z1
      CCTK_INT   :: ii,jj, kk, ib, jb, kb, ic, jc, kc
      CCTK_REAL  :: xp, yp, zp
      CCTK_REAL  :: cx(10), cy(10), cz(10)

      ! functions
      CCTK_INT   gft_read_rank, gft_read_shape, gft_read_full
      CCTK_REAL  myl2norm3d
      logical    double_equal

      interface
         subroutine int4coef(c, xi, x)
           implicit none
           CCTK_REAL c(5), xi, x(*)
         end subroutine int4coef
      end interface

#ifdef AIX
      rc = gft_read_rank(trim(filename)//CHAR(0), %VAL(level), rank)
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
         write(*,*) '...shp = ', shp(1),shp(2), shp(3)
      end if
      if (shp(1)*shp(2)*shp(3) .gt. dnx*dny*dnz) then
         write(*,*) 'fdata is not large enough to read in this file'
         write(*,*) 'file = ', filename
         write(*,*) 'size of fdata = ',dnx, dny, dnz
         write(*,*) 'size of data in file = ',shp(1),shp(2),shp(3)
         stop
      end if

      if (FLT) write(*,*) '...calling gft_read_full'
#ifdef AIX
      rc = gft_read_full(trim(filename)//CHAR(0),         &
                         %VAL(level), shp, cnames, rank,   &
                         time, coords, fdata)
#else
      rc = gft_read_full(filename, level, shp, cnames, rank, &
                         time, coords, fdata)
#endif
      if (FLT) write(*,*) '...end gft_read_full'

      if (rc .ne. 1) then
         write(*,*) 'Unable to read data in ',filename
         stop
      end if

      if (double_equal(myl2norm3d(fdata,shp(1),shp(2),shp(3)),0.d0)) &
         then
         if (ltrace) write(*,*) '...initializing to zero: ',filename
         call load_scal3d( f, 0.d0, nx, ny, nz)
         goto 100
      end if

      ! Extract local data from the global data set
      ilb = cctk_llb(1)
      jlb = cctk_llb(2)
      klb = cctk_llb(3)

      if (interp .eq. 0) then
         do i = 1, rank
           if (cctk_gshp(i) .ne. shp(i)) then
           write(*,*)'...Inconsistent data and parameter file.'
           write(*,*)'...file has ',shp(i), ', but expecting ',cctk_gshp(i)
           stop
           end if
         end do
         !
         ! load data in array
         if (FLT) write(*,*) '...straight loading data'
         do k = 1, nz
         do j = 1, ny
         do i = 1, nx
             f(i,j,k) = fdata(i+ilb, j+jlb, k+klb)
         end do
         end do
         end do
      else if (interp .eq. 1) then
         ! interpolation
         if (FLT) write(*,*) '...interpolating'
         dnx = shp(1)
         dny = shp(2)
         dnz = shp(3)
         ddx = coords(2) - coords(1)
         ddy = coords(dnx + 2) - coords(dnx + 1)
         ddz = coords(dnx + dny + 2) - coords(dnx + dny + 1)
         dxmin = coords(1)
         dymin = coords(dnx+1)
         dzmin = coords(dnx+dny+1)
         order = 4
         if (FLT) then
            write(*,*) 'dnx/y/z   = ',dnx,dny,dnz
            write(*,*) 'dx/y/zmin = ',dxmin,dymin,dzmin
            write(*,*) 'ddx/y/z   = ',ddx,ddy,ddz
         end if

         !write(*,*) 'delta x = ',xf(2)-xf(1)
         do k = 1, nz
            zp = zf(k)
            do j = 1, ny
               yp = yf(j)
               do i = 1, nx
                  xp = xf(i)
                  if (xp .lt. dxmin .or. yp .lt. dymin .or. zp .lt.dzmin) then
                      write(*,*) 'Problem: ',xp,yp,zp
                      stop
                  end if
                  ! estimate the indices for the point in fdata closest to 
                  ! (xp,yp,zp)
                  ii = int((xp - dxmin)/ddx) + 1
                  jj = int((yp - dymin)/ddy) + 1
                  kk = int((zp - dzmin)/ddz) + 1

                  ! find the first point in the interpolation stencil.
                  ! we are careful about boundaries.
                  ib = min(max(ii-2,1),dnx-4)
                  jb = min(max(jj-2,1),dny-4)
                  kb = min(max(Kk-2,1),dnz-4)

                  call int4coef(cx, xp, coords(ib))
                  call int4coef(cy, yp, coords(dnx + jb ))
                  call int4coef(cz, zp, coords(dnx + dny + kb))

                  f(i,j,k) = 0.0d0
                  do kc = 1, 5
                  do jc = 1, 5
                  do ic = 1, 5
                     f(i,j,k) = f(i,j,k)  &
                     + cx(ic)*cy(jc)*cz(kc)*fdata(ib+ic-1,jb+jc-1,kb+kc-1)
                  end do
                  end do
                  end do

               end do
            end do
         end do
      else
         write(*,*)'...interp has illegal value in read_file_sdf'
         write(*,*)'...interp = ',interp
         stop
      end if
      
100  continue

      return
   end subroutine read_sdf_file


!-------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------
   subroutine int4coef(c, xi, x)
     implicit none
     CCTK_REAL c(5), xi, x(*)


     ! local vars
     CCTK_INT :: i

     ! initialize all coefs to zero
     do i = 1, 5
        c(i) = 0.0d0
     end do
     do i = 1, 5
        if (xi .eq. x(i)) then
          c(i) = 1.0d0
          return
        end if
     end do

  c(1) =    (xi-x(2))*(xi-x(3))*(xi-x(4))*(xi-x(5)) &
         / ((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))*(x(1)-x(5)))

  c(2) =    (xi-x(1))*(xi-x(3))*(xi-x(4))*(xi-x(5)) &
         / ((x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))*(x(2)-x(5)))

  c(3) =    (xi-x(1))*(xi-x(2))*(xi-x(4))*(xi-x(5)) &
         / ((x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))*(x(3)-x(5)))

  c(4) =    (xi-x(1))*(xi-x(2))*(xi-x(3))*(xi-x(5)) &
         / ((x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))*(x(4)-x(5)))

  c(5) =    (xi-x(1))*(xi-x(2))*(xi-x(3))*(xi-x(4)) &
         / ((x(5)-x(1))*(x(5)-x(2))*(x(5)-x(3))*(x(5)-x(4)))

     return
   end subroutine int4coef
