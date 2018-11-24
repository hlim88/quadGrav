c-----------------------------------------------------------------------
c
c   $Id: futil.f,v 1.1 2010-07-17 22:09:42 carlos Exp $
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c
c 
c
c-----------------------------------------------------------------------
      subroutine check_derivs_1(rc, f, dxf, dyf, dzf, fname,
     &                          time, myproc, x1d, y1d, z1d, nx, ny, nz)
        implicit none
        integer rc, nx, ny, nz, myproc
        real*8  f(nx,ny,nz), dxf(nx,ny,nz), dyf(nx,ny,nz), dzf(nx,ny,nz)
        real*8  x1d(nx), y1d(ny), z1d(nz)
        real*8  time
        character*(*) fname

c       local vars
        integer nd, ok, rc1, rc2, rc3, rc4, i, j
        integer rank
        integer shp(3)
        character*128 sdfname
        character*16  cnames
        integer       maxc
        parameter(    maxc = 512   )
        real*8        coords(maxc)

        ok = 1
        nd = nx*ny*nz

        call check_isfinite_1d(f,   rc1, nd)
        ok = ok*rc1
        call check_isfinite_1d(dxf, rc2, nd)
        ok = ok*rc2
        call check_isfinite_1d(dyf, rc3, nd)
        ok = ok*rc3
        call check_isfinite_1d(dzf, rc4, nd)
        ok = ok*rc4

        if (ok .eq. 0) then
          write(*,*)'-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
          write(*,*)'CHECK_DERIVS_1: problem with ',fname
          write(*,*)'  rc[1-4] ',rc1,rc2,rc3,rc4
          shp(1) = nx
          shp(2) = ny
          shp(3) = nz
          rank   = 3
          cnames = 'x|y|z'
  
          i = nx + ny + nz
          if (i .gt. maxc) then
            write(0,*)'check_derivs_1: Not enough memory to write ',
     &                 'coordinates'
            write(0,*) 'to sdf output.  Will not write SDFs.'
            return
          end if

          j = 0
          do i = 1, nx
            j = j + 1
            coords(j) = x1d(i)
          end do
          do i = 1, ny
            j = j + 1
            coords(j) = y1d(i)
          end do
          do i = 1, nz
            j = j + 1
            coords(j) = z1d(i)
          end do


          if (myproc .lt.10) then
            write(sdfname,95) myproc, fname
 95         format(i1,'P_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,100) myproc, fname
 100        format(i2,'P_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,105) myproc, fname
 105        format(i3,'P_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,f)

          if (myproc .lt.10) then
            write(sdfname,110) myproc, fname
 110        format(i1,'P_DX_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,115) myproc, fname
 115        format(i2,'P_DX_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,120) myproc, fname
 120        format(i3,'P_DX_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,dxf)

          if (myproc .lt.10) then
            write(sdfname,125) myproc, fname
 125        format(i1,'P_DY_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,130) myproc, fname
 130        format(i2,'P_DY_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,135) myproc, fname
 135        format(i3,'P_DY_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,dyf)

          if (myproc .lt.10) then
            write(sdfname,140) myproc, fname
 140        format(i1,'P_DZ_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,145) myproc, fname
 145        format(i2,'P_DZ_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,150) myproc, fname
 150        format(i3,'P_DZ_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,dzf)
        end if

        rc = ok

        return
      end subroutine
c-----------------------------------------------------------------------
c
c 
c
c-----------------------------------------------------------------------
      subroutine check_derivs_2(rc, f, dxxf, dxyf, dxzf, 
     &                          dyyf, dyzf, dzzf,  fname,
     &                          time, myproc, x1d, y1d, z1d, nx, ny, nz)
        implicit none
        integer rc, nx, ny, nz, myproc
        real*8  f(nx,ny,nz)
        real*8  dxxf(nx,ny,nz), dxyf(nx,ny,nz), dxzf(nx,ny,nz)
        real*8  dyyf(nx,ny,nz), dyzf(nx,ny,nz), dzzf(nx,ny,nz)
        real*8  x1d(nx), y1d(ny), z1d(nz)
        real*8  time
        character*(*) fname

c       local vars
        integer nd, ok, rc1, rc2, rc3, rc4, rc5, rc6, rc7
        integer rank, i, j
        integer shp(3)
        character*128 sdfname
        character*16  cnames
        integer       maxc
        parameter(    maxc = 512   )
        real*8        coords(maxc)

        ok = 1
        nd = nx*ny*nz
        

        call check_isfinite_1d(f,   rc1, nd)
        ok = ok*rc1
        call check_isfinite_1d(dxxf, rc2, nd)
        ok = ok*rc2
        call check_isfinite_1d(dxyf, rc3, nd)
        ok = ok*rc3
        call check_isfinite_1d(dxzf, rc4, nd)
        ok = ok*rc4
        call check_isfinite_1d(dyyf, rc5, nd)
        ok = ok*rc5
        call check_isfinite_1d(dyzf, rc6, nd)
        ok = ok*rc6
        call check_isfinite_1d(dzzf, rc7, nd)
        ok = ok*rc7


        if (ok .eq. 0) then
          write(*,*)'-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
          write(*,*)'CHECK_DERIVS_2: problem with ',fname
          write(*,*)'  rc[1-4] ',rc1,rc2,rc3,rc4

          i = nx + ny + nz
          if (i .gt. maxc) then
            write(0,*)'check_derivs_2: Not enough memory to write ',
     &                 'coordinates'
            write(0,*) 'to sdf output.  Will not write SDFs.'
            return
          end if

          shp(1) = nx
          shp(2) = ny
          shp(3) = nz
          rank   = 3
          cnames = 'x|y|z'
          j = 0
          do i = 1, nx
            j = j + 1
            coords(j) = x1d(i)
          end do
          do i = 1, ny
            j = j + 1
            coords(j) = y1d(i)
          end do
          do i = 1, nz
            j = j + 1
            coords(j) = z1d(i)
          end do

          if (myproc .lt.10) then
            write(sdfname,195) myproc, fname
 195        format(i1,'Q_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,200) myproc, fname
 200        format(i2,'Q_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,205) myproc, fname
 205        format(i3,'Q_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,f)

          if (myproc .lt.10) then
            write(sdfname,210) myproc, fname
 210        format(i1,'Q_DXX_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,215) myproc, fname
 215        format(i2,'Q_DXX_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,220) myproc, fname
 220        format(i3,'Q_DXX_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,dxxf)

          if (myproc .lt.10) then
            write(sdfname,225) myproc, fname
 225        format(i1,'Q_DXY_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,230) myproc, fname
 230        format(i2,'Q_DXY_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,235) myproc, fname
 235        format(i3,'Q_DXY_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,dxyf)

          if (myproc .lt.10) then
            write(sdfname,240) myproc, fname
 240        format(i1,'Q_DXZ_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,245) myproc, fname
 245        format(i2,'Q_DXZ_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,250) myproc, fname
 250        format(i3,'Q_DXZ_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,dxzf)

          if (myproc .lt.10) then
            write(sdfname,255) myproc, fname
 255        format(i1,'Q_DYY_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,260) myproc, fname
 260        format(i2,'Q_DYY_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,265) myproc, fname
 265        format(i3,'Q_DYY_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,dyyf)

          if (myproc .lt.10) then
            write(sdfname,270) myproc, fname
 270        format(i1,'Q_DYZ_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,275) myproc, fname
 275        format(i2,'Q_DYZ_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,280) myproc, fname
 280        format(i3,'Q_DYZ_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,dyzf)

          if (myproc .lt.10) then
            write(sdfname,285) myproc, fname
 285        format(i1,'Q_DZZ_',a)
          else if (myproc .gt. 9 .and. myproc .lt. 100) then
            write(sdfname,290) myproc, fname
 290        format(i2,'Q_DZZ_',a)
          else if (myproc .gt. 99 .and. myproc .lt. 1000) then
            write(sdfname,295) myproc, fname
 295        format(i3,'Q_DZZ_',a)
          else
          end if
          call gft_out_full(sdfname,time,shp,cnames,rank,coords,dzzf)


        end if

        rc = ok

        return
      end subroutine
c-----------------------------------------------------------------------
c
c 
c
c-----------------------------------------------------------------------
      subroutine check_magnetic(rc, Bx, By, Bz, nx, ny, nz, string)
        implicit none

        integer       nx, ny, nz, rc
        real*8        Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
        character*(*) string

c       local vars
        integer       err, i, j, k
        real*8        myl2norm3d
        external      myl2norm3d

        err = 0

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          if (abs(Bx(i,j,k)) .gt. 0.0d0) then 
            err = 1
          end if
        end do
        end do
        end do
        if (err .eq. 1) then
          write(*,*)'check_magnetic QQ'
          write(*,*) string
          write(*,*)' QQ Bx is not zero.  |Bx|=',
     &                     myl2norm3d(Bx, nx, ny, nz)
        end if

        err = 0
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          if (abs(By(i,j,k)) .gt. 0.0d0) then 
            err = 1
          end if
        end do
        end do
        end do
        if (err .eq. 1) then
          write(*,*)'check_magnetic QQ'
          write(*,*) string
          write(*,*)' QQ By is not zero.  |By|=',
     &                     myl2norm3d(By, nx, ny, nz)
        end if

        err = 0
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
          if (abs(Bz(i,j,k)) .gt. 0.0d0) then 
            err = 1
          end if
        end do
        end do
        end do
        if (err .eq. 1) then
          write(*,*)'check_magnetic QQ'
          write(*,*) string
          write(*,*)' QQ Bz is not zero.  |Bz|=',
     &                     myl2norm3d(Bz, nx, ny, nz)

        end if

        return
      end subroutine

c-----------------------------------------------------------------------
c
c     subroutine check_isfinite
c
c-----------------------------------------------------------------------
      subroutine check_isfinite(x, string, rc)
        implicit none
        real*8        x
        character*(*) string
        integer       rc

        call for_isfinite(x, rc)

        if (rc .eq. 0) then
          write(*,*)'check_isfinite: value not finite',x
          write(*,*)'check_isfinite: message: ',string
        end if

        return
      end subroutine


c-----------------------------------------------------------------------
c
c 
c
c-----------------------------------------------------------------------
