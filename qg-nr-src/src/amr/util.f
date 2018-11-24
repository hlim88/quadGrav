c. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .
c.   This file, util.f, contains utility     .
c.   routines used by this code but does     .
c.   not include any files, and takes all    .
c.   input as arguments. Hence, these        .
c.   routines don't know what AMR or grids   .
c.   are. Hence, the name utility routines.  .
c. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  paste_boundariesFV:                                                       cc
cc                                                                            cc
cc      NB: the loops can't go up to maxi/j/k, but only to maxi/j/k-1         cc
cc      NB: customized version of paste_boundaries() for finite volume        cc
cc          fields. No change to small x/y/z boundaries. For large            cc
cc          x/y/z boundaries, the change is that we really only have          cc
cc          nx/y/z-1 cells instead of nx/y/z points.                          cc
cc      copy data from work space into boundary data of a field               cc
cc      assuming that the necessary interpolation has happened everywhere.    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine paste_boundariesFV(  values,field,     chr,
     *                              mini,  minj,      mink,
     *                              maxi,  maxj,      maxk,
     *                              nx,    ny,        nz,
     *                              r, gw)
      implicit    none
      integer     nx,ny,nz, r, gw
      integer     mini, minj, mink,  maxi, maxj, maxk
      integer     minip,minjp,minkp, maxip,maxjp,maxkp
      real(kind=8)      values(nx,ny,nz),
     *            field(nx,ny,nz), chr(nx,ny,nz)
      include    'chr.inc'
      integer     i,j,k, width
      !
      ! To keep things as close to the regular version of this routine,
      ! I'll just subtract one from each of the bounds when setting 
      ! the boundary of the "high" side...because for cell centered,
      ! one has nx-1 instead of nx cells.
      !
      !
      logical     ltrace
      parameter ( ltrace  = .false. )

      ! DM This routine was not covering all the points, there were zeroes
      ! turning up on FV AMR boundaries that shouldn't be there. I'm not
      ! sure if my correction is correct, either, as it may take points
      ! from parents that are in the last row of cells. But it appears to 
      ! be working for now.

      if (ltrace) then
         write(*,*) 'paste_boundariesFV: nx/y/z:   ',nx, ny, nz
         write(*,*) 'paste_boundariesFV: mini/j/k: ',mini,minj,mink
         write(*,*) 'paste_boundariesFV: maxi/j/k: ',maxi,maxj,maxk
      end if

      !
      ! Low Z Boundary:
      !
      if (mink .ne. 1) goto 10
      width = gw
      if (width .gt. maxk) width = maxk
      if (ltrace) write(*,*) 'paste_boundariesFV: Low  Z k: ',1,width
      do k = 1, width
         do j = minj, maxj
            do i = mini, maxi
               if ( NINT(chr(i,j,1)) .ne. CHR_deco_bdy ) then
               !if (.not.double_equal(chr(i,j,1),CHR_deco_bdy)) then
                  field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

 10   continue

      !
      ! High Z Boundary:
      !
      if (maxk .ne. nz) goto 20
      width = nz-gw
      if (width .lt. mink) width = mink
      if (ltrace) write(*,*) 'paste_boundariesFV: High Z k: ',nz, width
      do k= nz, width, -1
         do j = minj, maxj
            do i = mini, maxi
               if ( NINT(chr(i,j,nz)) .ne. CHR_deco_bdy ) then
                  field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

 20   continue

      !
      ! Low Y Boundary:
      !
      if (minj .ne. 1) goto 30
      width = gw
      if (width .gt. maxj) width = maxj
      if (ltrace) write(*,*) 'paste_boundariesFV: Low  Y j: ',1,width
      do j = 1, width
         do k = mink, maxk
            do i = mini, maxi
               if ( NINT(chr(i,1,k)) .ne. CHR_deco_bdy ) then
                  field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

 30   continue

      !
      ! High Y Boundary:
      !
      if (maxj .ne. ny) goto 40
      width = ny-gw
      if (width .lt. minj) width = minj
      if (ltrace) write(*,*) 'paste_boundariesFV: High Y j: ',ny,width
      do j= ny, width, -1
         do k = mink, maxk
            do i = mini, maxi
               if ( NINT(chr(i,ny,k)) .ne. CHR_deco_bdy ) then
                  field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

 40   continue

      !
      ! Low X Boundary:
      !
      if (mini .ne. 1) goto 50
      width = gw
      if (width .gt. maxi) width = maxi
      if (ltrace) write(*,*) 'paste_boundariesFV: Low  X i: ',1,width
      do i = 1, width
         do k = mink, maxk
            do j = minj, maxj
               if ( NINT(chr(1,j,k)) .ne. CHR_deco_bdy ) then
                  field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

 50   continue

      !
      ! High X Boundary:
      !
      if (maxi .ne. nx) goto 60
      width = nx-gw
      if (width .lt. mini) width = mini
      if (ltrace) write(*,*) 'paste_boundariesFV: High X i: ',nx,width
      do i= nx, width, -1
         do k = mink, maxk
            do j = minj, maxj
               if ( NINT(chr(nx,j,k)) .ne. CHR_deco_bdy ) then
                  field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

 60   continue

      return
      end    ! END: paste_boundariesFV

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  init_fieldFV:                                                             cc
cc                                                                            cc
cc           Initialize a finite volume field over at least part of its       cc
cc           domain from data on the parent.                                  cc
cc               Input:   u_c --- coarse grid data in the region of overlap   cc
cc               Output:  u_f --- fine   grid data properly initialized via   cc
cc                                striaght copies of parent data in the cell  cc
cc           NB: for FV fields only                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_fieldFV( u_c, u_f, nxc,nyc,nzc,
     *                         mini,minj,mink,
     *                         maxi,maxj,maxk,
     *                         nxf, nyf, nzf, r, gw)
      implicit    none
      integer     nxc,nyc,nzc, nxf,nyf,nzf, r, gw
      integer     mini,minj,mink
      integer     maxi,maxj,maxk
      real(kind=8)      u_c(nxc,nyc,nzc), u_f(nxf,nyf,nzf)
      integer     ic,jc,kc, if,jf,kf, i
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      !
      ! (1) This is restricited to a refine
      !     factor of two, at least  for right now
      ! (2) It is not clear yet whether this covers
      !     all the points yet
      !
      !

      ! DM This routine was not covering all the points, there were zeroes
      ! turning up on FV AMR boundaries that shouldn't be there. I'm not
      ! sure if my correction is correct, either, as it may take points
      ! from parents that are in the last row of cells. But it appears to 
      ! be working for now.

c     .....copy points in common to both grids
      do kf = mink, maxk
         kc = (kf - mink)/r + 1
         do jf = minj, maxj
            jc = (jf - minj)/r + 1
            do if = mini, maxi
               ic = (if - mini)/r + 1
               u_f(if,  jf,  kf  ) = u_c(ic,jc,kc)
            end do
         end do
      end do


      return
      end    ! END: init_fieldFV


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc int2str:                                                                   cc
cc           Converts a non-negative integer to a string.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine int2str( myint, mystring )
      implicit none
      character*(*) mystring
      integer       myint
      integer       tmp, numdigits, the_int, poweroften, digit
      character     dig2char
      external      dig2char

      mystring = ''
      if (myint .lt. 0) then
         write(*,*) 'int2str: Cannot convert negative integer.'
         return
      else if (myint .lt. 10) then
         mystring = dig2char(myint)
         return
      end if

      tmp       = 10
      numdigits = 1
 10   if (myint .ge. tmp) then
         tmp       = 10 * tmp
         numdigits = numdigits + 1
         goto 10
      end if

      tmp               = 1
      the_int           = myint
      poweroften        = 10**(numdigits-1)

 20   digit             = the_int / poweroften
      mystring(tmp:tmp) = dig2char(digit)
      tmp               = tmp + 1
      the_int           = the_int - digit * poweroften
      poweroften        =  poweroften / 10
      if (tmp.le.numdigits) goto 20

      mystring(tmp:tmp) = ''

      return
      end        ! END: subroutine int2str

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc dig2char:                                                                  cc
cc           Converts a single digit [0..9] to a character.                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character function dig2char( digit )
      implicit none
      integer       digit

      if      (digit .eq. 0) then
              dig2char = '0'
      else if (digit .eq. 1) then
              dig2char = '1'
      else if (digit .eq. 2) then
              dig2char = '2'
      else if (digit .eq. 3) then
              dig2char = '3'
      else if (digit .eq. 4) then
              dig2char = '4'
      else if (digit .eq. 5) then
              dig2char = '5'
      else if (digit .eq. 6) then
              dig2char = '6'
      else if (digit .eq. 7) then
              dig2char = '7'
      else if (digit .eq. 8) then
              dig2char = '8'
      else if (digit .eq. 9) then
              dig2char = '9'
      end if

      return
      end        ! END: function dig2char

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc mystringlength:                                                            cc
cc                 Returns length of beginning non-space characters.          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function mystringlength( the_string )
      implicit none
      character*(*) the_string
      integer       the_string_length

      the_string_length = len(the_string)

      mystringlength    = index(the_string,' ') - 1

      if (mystringlength .lt. 0) mystringlength = the_string_length

      return
      end        ! END: function mystringlength

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interp_cubicshort                                                         cc
cc                                                                            cc
cc    From interp_cubic(): assumes two things:                                cc
cc        (1) y-values are equally spaced                                     cc
cc        (2) interpolated value in center between the other four             cc
cc    Finds cubic interpolated point at x based on points 1,2,3,4             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function interp_cubicshort(y1,y2,y3,y4)
      implicit    none
      real(kind=8)      y1, y2, y3, y4

      interp_cubicshort =  -        y1 / 16.d0
     *
     *                     + 9.d0 * y2 / 16.d0
     *
     *                     + 9.d0 * y3 / 16.d0
     *
     *                     -        y4 / 16.d0

c     interp_cubic =    xx2    *xx3    *xx4     * y1
c    *                / ( (x1-x2)*(x1-x3)*(x1-x4) )
c    *
c    *              +   xx1    *xx3    *xx4     * y2
c    *                / ( (x2-x1)*(x2-x3)*(x2-x4) )
c    *
c    *              +   xx1    *xx2    *xx4     * y3
c    *                / ( (x3-x1)*(x3-x2)*(x3-x4) )
c    *
c    *              +   xx1    *xx2    *xx3     * y4
c    *                / ( (x4-x1)*(x4-x2)*(x4-x3) )

      return
      end    ! END function interp_cubicshort

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interp_cubicb:                                                            cc
cc                                                                            cc
cc    Finds cubic interpolated point at x based on points 1,2,3,4             cc
cc    NB: this version checks for any points that are 'masked' out            cc
cc        by having been set to a LARENUMBER. In that event, reduced          cc
cc        order interpolation is done.                                        cc
cc    NB: since this is used sometimes for masked regions, do only            cc
cc        one-sided interpolations, never skipping a point.                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function interp_cubicb(y1,y2,y3,y4,x,x1,x2,x3,x4)
      implicit          none
      real(kind=8)      y1, y2, y3, y4, x, x1, x2, x3, x4
      real(kind=8)      xx1, xx2, xx3, xx4
      include          'largesmall.inc'

      logical     ltrace
      parameter ( ltrace = .false. )

      xx1  = x  - x1
      xx2  = x  - x2
      xx3  = x  - x3
      xx4  = x  - x4

      if (y1 .ge. LARGE_M_SMALL) then
         if (y2 .ge. LARGE_M_SMALL) then
            if (y3 .ge. LARGE_M_SMALL) then
               interp_cubicb=                        y4
            else
               if (y4 .ge. LARGE_M_SMALL) then
                  interp_cubicb=                     y3
               else
                  interp_cubicb=    
     *                                     xx4     * y3
     *                   / (                 (x3-x4) )
     *                 +                   xx3     * y4
     *                   / (                 (x4-x3) )
               end if
            end if
         else
            if (y3 .ge. LARGE_M_SMALL) then
               if (y4 .ge. LARGE_M_SMALL) then
                  interp_cubicb=                     y2
               else
                  interp_cubicb=                     y2
c                 interp_cubicb=    
c    *                                     xx4     * y2
c    *                   / (                 (x2-x4) )
c    *                 +                   xx2     * y4
c    *                   / (                 (x4-x2) )
               end if
            else
               if (y4 .ge. LARGE_M_SMALL) then
                  interp_cubicb=
     *                 +           xx3             * y2
     *                   / (         (x2-x3)         )
     *                 +           xx2             * y3
     *                   / (         (x3-x2)         )
               else
                  interp_cubicb=
     *                 +           xx3    *xx4     * y2
     *                   / (         (x2-x3)*(x2-x4) )
     *                 +           xx2    *xx4     * y3
     *                   / (         (x3-x2)*(x3-x4) )
     *                 +           xx2    *xx3     * y4
     *                   / (         (x4-x2)*(x4-x3) )
               end if
            end if
         end if
      else
         if (y2 .ge. LARGE_M_SMALL) then
            if (y3 .ge. LARGE_M_SMALL) then
               if (y4 .ge. LARGE_M_SMALL) then
                  interp_cubicb=                     y1
               else
                  interp_cubicb=                     y1
c                 interp_cubicb=
c    *                 +           xx4             * y1
c    *                   / (         (x1-x4)         )
c    *                 +           xx1             * y4
c    *                   / (         (x4-x1)         )
               end if
            else
               if (y4 .ge. LARGE_M_SMALL) then
                  interp_cubicb=                     y1
c                 interp_cubicb=
c    *                 +           xx3             * y1
c    *                   / (         (x1-x3)         )
c    *                 +           xx1             * y3
c    *                   / (         (x3-x1)         )
               else
                  ! could use y1 here:
                  interp_cubicb=
     *                 +           xx4             * y3
     *                   / (         (x3-x4)         )
     *                 +           xx3             * y4
     *                   / (         (x4-x3)         )
               end if
            end if
         else
            if (y3 .ge. LARGE_M_SMALL) then
               if (y4 .ge. LARGE_M_SMALL) then
                  interp_cubicb=
     *                 +           xx2             * y1
     *                   / (         (x1-x2)         )
     *                 +           xx1             * y2
     *                   / (         (x2-x1)         )
               else
                  ! could use y4 here
                  interp_cubicb=
     *                 +           xx2             * y1
     *                   / (         (x1-x2)         )
     *                 +           xx1             * y2
     *                   / (         (x2-x1)         )
               end if
            else
               if (y4 .ge. LARGE_M_SMALL) then
                  interp_cubicb=
     *                 +           xx2    *xx3     * y1
     *                   / (         (x1-x2)*(x1-x3) )
     *                 +           xx1    *xx3     * y2
     *                   / (         (x2-x1)*(x2-x3) )
     *                 +           xx1    *xx2     * y3
     *                   / (         (x3-x1)*(x3-x2) )
               else
                  interp_cubicb = 
     *                     xx2    *xx3    *xx4     * y1
     *                   / ( (x1-x2)*(x1-x3)*(x1-x4) )
     *
     *                 +   xx1    *xx3    *xx4     * y2
     *                   / ( (x2-x1)*(x2-x3)*(x2-x4) )
     *
     *                 +   xx1    *xx2    *xx4     * y3
     *                   / ( (x3-x1)*(x3-x2)*(x3-x4) )
     *
     *                 +   xx1    *xx2    *xx3     * y4
     *                   / ( (x4-x1)*(x4-x2)*(x4-x3) )
               end if
            end if
         end if
      end if

         if (ltrace .and. interp_cubicb .ge. 1.d80) then
            write(*,*) 'interp_cubicb: Invalid interpolation:'
            write(*,*) 'interp_cubicb: y1 = ',y1
            write(*,*) 'interp_cubicb: y2 = ',y2
            write(*,*) 'interp_cubicb: y3 = ',y3
            write(*,*) 'interp_cubicb: y4 = ',y4
         end if

      if (ltrace) then
         write(*,*) 'interp_cubicb: y1 = ', y1
         write(*,*) 'interp_cubicb: y2 = ', y2
         write(*,*) 'interp_cubicb: y3 = ', y3
         write(*,*) 'interp_cubicb: y4 = ', y4
         write(*,*) 'interp_cubicb: x1 = ', x1
         write(*,*) 'interp_cubicb: x2 = ', x2
         write(*,*) 'interp_cubicb: x3 = ', x3
         write(*,*) 'interp_cubicb: x4 = ', x4
         write(*,*) 'interp_cubicb: x  = ', x
         write(*,*) 'interp_cubicb: y  = ', interp_cubicb
      end if


      return
      end    ! END function interp_cubicb


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  interp_cubic:                                                             cc
cc                                                                            cc
cc    Finds cubic interpolated point at x based on points 1,2,3,4             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function interp_cubic(y1,y2,y3,y4, x,x1,x2,x3,x4)
      implicit    none
      real(kind=8)      y1, y2, y3, y4, x, x1, x2, x3, x4
      real(kind=8)      xx1, xx2, xx3, xx4


      logical     ltrace
      parameter ( ltrace = .false. )

      xx1  = x  - x1
      xx2  = x  - x2
      xx3  = x  - x3
      xx4  = x  - x4

      interp_cubic =    xx2    *xx3    *xx4     * y1
     *                / ( (x1-x2)*(x1-x3)*(x1-x4) )
     *
     *              +   xx1    *xx3    *xx4     * y2
     *                / ( (x2-x1)*(x2-x3)*(x2-x4) )
     *
     *              +   xx1    *xx2    *xx4     * y3
     *                / ( (x3-x1)*(x3-x2)*(x3-x4) )
     *
     *              +   xx1    *xx2    *xx3     * y4
     *                / ( (x4-x1)*(x4-x2)*(x4-x3) )
      if (ltrace) then
         write(*,*) 'interp_cubic: y1 = ', y1
         write(*,*) 'interp_cubic: y2 = ', y2
         write(*,*) 'interp_cubic: y3 = ', y3
         write(*,*) 'interp_cubic: y4 = ', y4
         write(*,*) 'interp_cubic: x1 = ', x1
         write(*,*) 'interp_cubic: x2 = ', x2
         write(*,*) 'interp_cubic: x3 = ', x3
         write(*,*) 'interp_cubic: x4 = ', x4
c        write(*,*) 'interp_cubic: xx1 = ', xx1
c        write(*,*) 'interp_cubic: xx2 = ', xx2
c        write(*,*) 'interp_cubic: xx3 = ', xx3
c        write(*,*) 'interp_cubic: xx4 = ', xx4
         write(*,*) 'interp_cubic: x  = ', x
         write(*,*) 'interp_cubic: y  = ', interp_cubic
      end if


      return
      end    ! END function interp_cubic

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  copy_interp_section:                                                      cc
cc                                                                            cc
cc     Copies data into a field and interpolates around the boundary          cc
cc     of this region. The reason for the interpolation is to smooth          cc
cc     out any irregularities. The points interpolated were almost            cc
cc     surely gotten originally by interpolating from data from a parent,     cc
cc     and hence this is just reinterpolating based on the new information.   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copy_interp_section( field_to, field_from, ox,oy,oz,
     *                                nx_f,ny_f,nz_f, nx_t,ny_t,nz_t,r)
      implicit    none
      integer     ox,oy,oz, nx_f,ny_f,nz_f, nx_t,ny_t,nz_t, r
      real(kind=8)      field_to(nx_t,ny_t,nz_t)
      real(kind=8)      field_from(nx_f,ny_f,nz_f)
      integer     i,j,k, if,jf,kf, fi,fj,fk

      !
      ! does not play well with domain decomposition:
      !
      logical     interp_boundary
      parameter ( interp_boundary = .false. )
 
      logical     ltrace
      parameter ( ltrace = .false. )

      logical     ltrace2
      parameter ( ltrace2 = .false. )

      if (ltrace) write(*,*) 'copy_interp_section: Enter'

      !
      ! Simple bounds checking:
      !
      if (     (nz_t.lt.nz_f+oz)
     *     .or.(ny_t.lt.ny_f+oy)
     *     .or.(nx_t.lt.nx_f+ox) ) then
          write(*,*) '                  nx_f  = ',nx_f
          write(*,*) '                  ny_f  = ',ny_f
          write(*,*) '                  nz_f  = ',nz_f
          write(*,*) '                  nx_t  = ',nx_t
          write(*,*) '                  ny_t  = ',ny_t
          write(*,*) '                  nz_t  = ',nz_t
          write(*,*) '                    ox  = ',ox
          write(*,*) '                    oy  = ',oy
          write(*,*) '                    oz  = ',oz
          write(*,*) '       refine factor: r = ',r
         write(*,*) 'copy_interp_section: PROBLEM!!!!'
         stop
      end if

      !
      ! With bounding box, copy values:
      !
      do kf = 1, nz_f
         k  = kf + oz
         do jf = 1, ny_f
            j  = jf + oy
            do if = 1, nx_f
               i  = if + ox
               field_to(i,j,k) = field_from(if,jf,kf)
            end do
         end do
      end do

      ! Upper bounds of region copied:
      fi = ox + nx_f
      fj = oy + ny_f
      fk = oz + nz_f
      !
      ! Check for complete overlap:
      !

      if (      (ox.eq.0) .and. (fi.eq.nx_t)
     *     .and.(oy.eq.0) .and. (fj.eq.ny_t)
     *     .and.(oz.eq.0) .and. (fk.eq.nz_t) ) then
         if (ltrace2) write(*,*) 'copy_interp: complete overlap'
         goto 100
      end if

      if (.not. interp_boundary) goto 100

      !
      ! Reinterpolate around "x" boundary:
      !
      if (ox .gt. 1) then
         if (ltrace2) write(*,*) 'copy_interp: reinterp at xmin'
          i = ox
          do k = 1, fk
             do j = 1, fj
               field_to(i,j,k) = 0.5d0 * 
     *                         ( field_to(i+1,j,k) + field_to(i-1,j,k) )
             end do
          end do
      end if
      if (fi .lt. nx_t) then
         if (ltrace2) write(*,*) 'copy_interp: reinterp at xmax'
          i = fi + 1
          do k = 1, fk
             do j = 1, fj
               field_to(i,j,k) = 0.5d0 * 
     *                         ( field_to(i+1,j,k) + field_to(i-1,j,k) )
             end do
          end do
      end if

      !
      ! Reinterpolate around "y" boundary:
      !
      if (oy .gt. 1) then
         if (ltrace2) write(*,*) 'copy_interp: reinterp at ymin'
          j = oy
          do k = 1, fk
             do i = 1, fi
               field_to(i,j,k) = 0.5d0 * 
     *                         ( field_to(i,j+1,k) + field_to(i,j-1,k) )
             end do
          end do
      end if
      if (fj .lt. ny_t) then
         if (ltrace2) write(*,*) 'copy_interp: reinterp at ymax'
          j = fj + 1
          do k = 1, fk
             do i = 1, fi
               field_to(i,j,k) = 0.5d0 * 
     *                         ( field_to(i,j+1,k) + field_to(i,j-1,k) )
             end do
          end do
      end if

      !
      ! Reinterpolate around "z" boundary:
      !
      if (oz .gt. 1) then
         if (ltrace2) write(*,*) 'copy_interp: reinterp at zmin'
          k = oz
          do j = 1, fj
             do i = 1, fi
               field_to(i,j,k) = 0.5d0 * 
     *                         ( field_to(i,j,k+1) + field_to(i,j,k-1) )
             end do
          end do
      end if
      if (fk .lt. nz_t) then
         if (ltrace2) write(*,*) 'copy_interp: reinterp at zmax'
          k = fk + 1
          do j = 1, fj
             do i = 1, fi
               field_to(i,j,k) = 0.5d0 * 
     *                         ( field_to(i,j,k+1) + field_to(i,j,k-1) )
             end do
          end do
      end if

 100  if (ltrace) then
          write(*,*) '                  nx_f  = ',nx_f
          write(*,*) '                  ny_f  = ',ny_f
          write(*,*) '                  nz_f  = ',nz_f
          write(*,*) '                  nx_t  = ',nx_t
          write(*,*) '                  ny_t  = ',ny_t
          write(*,*) '                  nz_t  = ',nz_t
          write(*,*) '                    ox  = ',ox
          write(*,*) '                    oy  = ',oy
          write(*,*) '                    oz  = ',oz
          write(*,*) '       refine factor: r = ',r
          write(*,*) 'copy_interp_section: Done.'
      end if

      return
      end    ! END: copy_interp_section

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  apply_dissGW:                                                             cc
cc                                                                            cc
cc     NB: ghostwidth aware. Only applies dissipation within ghostregion.     cc
cc                                                                            cc
cc     Applies Kreiss-Oliger dissipation to a field.                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine apply_dissGW(field, work, epsdis, nx,   ny,  nz, gw)
      implicit    none
      integer     nx,ny,nz, gw
      real(kind=8)      field(nx,ny,nz), work(nx,ny,nz), epsdis
      integer     i,j,k, itmp

      !
      ! Make copy of "field":
      !
      call mat_copy3d(field, work, nx,ny,nz)

      if ( gw .ge. 2) then
         itmp = gw + 1
      else
         itmp = 3
      end if

      do k = itmp, nz-itmp+1
         do j = itmp, ny-itmp+1
            do i = itmp, nx-itmp+1
               field(i,j,k) = work(i,j,k) - (epsdis/16.d0)*(
     *             6.d0 *  work(i,j,k)   + work(i+2,j,k) + work(i-2,j,k)
     *           - 4.d0 * (work(i+1,j,k) + work(i-1,j,k) )
     *           + 6.d0 *  work(i,j,k)   + work(i,j+2,k) + work(i,j-2,k)
     *           - 4.d0 * (work(i,j+1,k) + work(i,j-1,k) )
     *           + 6.d0 *  work(i,j,k)   + work(i,j,k+2) + work(i,j,k-2)
     *           - 4.d0 * (work(i,j,k+1) + work(i,j,k-1) )
     *           )
            end do
         end do
      end do

      return
      end    ! END: apply_diss

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  apply_diss:                                                               cc
cc                                                                            cc
cc     Applies Kreiss-Oliger dissipation to a field.                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine apply_diss(   field, work, epsdis, nx,   ny,  nz)
      implicit    none
      integer     nx,ny,nz
      real(kind=8)      field(nx,ny,nz), work(nx,ny,nz), epsdis
      integer     i,j,k

      !
      ! Make copy of "field":
      !
      call mat_copy3d(field, work, nx,ny,nz)

      do k = 3, nz-2
         do j = 3, ny-2
            do i = 3, nx-2
               field(i,j,k) = work(i,j,k) - (epsdis/16.d0)*(
     *             6.d0 *  work(i,j,k)   + work(i+2,j,k) + work(i-2,j,k)
     *           - 4.d0 * (work(i+1,j,k) + work(i-1,j,k) )
     *           + 6.d0 *  work(i,j,k)   + work(i,j+2,k) + work(i,j-2,k)
     *           - 4.d0 * (work(i,j+1,k) + work(i,j-1,k) )
     *           + 6.d0 *  work(i,j,k)   + work(i,j,k+2) + work(i,j,k-2)
     *           - 4.d0 * (work(i,j,k+1) + work(i,j,k-1) )
     *           )
            end do
         end do
      end do

      return
      end    ! END: apply_diss

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  apply_dissnh:                                                             cc
cc                                                                            cc
cc     Applies Kreiss-Oliger dissipation to a field using a wide stencil      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine apply_dissnh(   field, work, epsdis, n, nx,   ny,  nz)
      implicit    none
      integer     nx,ny,nz, n
      real(kind=8)      field(nx,ny,nz), work(nx,ny,nz), epsdis
      integer     i,j,k

      !
      ! Make copy of "field":
      !
      call mat_copy3d(field, work, nx,ny,nz)

      do k = 1+2*n, nz-2*n
         do j = 1+2*n, ny-2*n
            do i = 1+2*n, nx-2*n
               field(i,j,k) = work(i,j,k) - (epsdis/16.d0)*(
     *             6.d0* work(i,j,k)  +work(i+2*n,j,k) + work(i-2*n,j,k)
     *           - 4.d0*(work(i+n,j,k)+work(i-n,j,k) )
     *           + 6.d0* work(i,j,k)  +work(i,j+2*n,k) + work(i,j-2*n,k)
     *           - 4.d0*(work(i,j+n,k)+work(i,j-n,k) )
     *           + 6.d0* work(i,j,k)  +work(i,j,k+2*n) + work(i,j,k-2*n)
     *           - 4.d0*(work(i,j,k+n)+work(i,j,k-n) )
     *           )
            end do
         end do
      end do

      return
      end    ! END: apply_dissnh

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  apply_diss_torhs:                                                         cc
cc                                                                            cc
cc     Computes Kreiss-Oliger dissipation of a field and                      cc
cc     adds it to the RHS.  For use in had's rk3.f90 file.                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine apply_diss_torhs( RHS, u, epsdis, nx,   ny,  nz)
      implicit    none
      integer     nx,ny,nz
      real(kind=8)      RHS(nx,ny,nz), u(nx,ny,nz), epsdis
      integer     i,j,k

      do k = 3, nz-2
         do j = 3, ny-2
            do i = 3, nx-2
               RHS(i,j,k) = RHS(i,j,k) - (epsdis/16.d0)*(
     *             6.d0 *  u(i,j,k)   + u(i+2,j,k) + u(i-2,j,k)
     *           - 4.d0 * (u(i+1,j,k) + u(i-1,j,k) )
     *           + 6.d0 *  u(i,j,k)   + u(i,j+2,k) + u(i,j-2,k)
     *           - 4.d0 * (u(i,j+1,k) + u(i,j-1,k) )
     *           + 6.d0 *  u(i,j,k)   + u(i,j,k+2) + u(i,j,k-2)
     *           - 4.d0 * (u(i,j,k+1) + u(i,j,k-1) )
     *           )
            end do
         end do
      end do

      return
      end    ! END: apply_diss_torhs

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  apply_diss2D:                                                             cc
cc                                                                            cc
cc     Applies Kreiss-Oliger dissipation to a field on a 2D surface           cc
cc     in spherical coordinates and hence has no boundaries.                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine apply_diss2D(field, work, epsdis, nx,   ny)
      implicit    none
      integer     nx,ny
      real(kind=8)      field(nx,ny), work(nx,ny), epsdis
      integer     j,k
      !real(kind=8) work_jp2, work_jp1, work_jm1, work_jm2
      !real(kind=8) work_kp2, work_kp1, work_km1, work_km2
      integer     km1,km2,kp1,kp2
      integer     jm1,jm2,jp1,jp2

      !
      ! Make copy of "field":
      !
      call mat_copy1d(field, work, nx*ny)

      do k = 1, ny
         km1 = k-1
         km2 = k-2
         kp1 = k+1
         kp2 = k+2
         if (k.eq.1) then
            km1 = ny
            km2 = ny-1
         else if (k.eq.2) then
            km2 = ny
         else if (k.eq.ny) then
            kp1 = 1
            kp2 = 2
         else if (k.eq.ny-1) then
            kp2 = 1
         end if
         !
         do j = 1, nx
            jm1 = j-1
            jm2 = j-2
            jp1 = j+1
            jp2 = j+2
            if (j.eq.1) then
               jm1 = nx
               jm2 = nx-1
            else if (j.eq.2) then
               jm2 = nx
            else if (j.eq.nx) then
               jp1 = 1
               jp2 = 2
            else if (j.eq.nx-1) then
               jp2 = 1
            end if
            !
            !
            field(j,k) =   work(j,k) - (epsdis/16.d0)*(
     *           + 6.d0 *  work(j,k) + work(jp2,k) + work(jm2,k)
     *           - 4.d0 * (work(jp1,k) + work(jm1,k) )
     *           + 6.d0 *  work(j,k) + work(j,kp2) + work(j,km2)
     *           - 4.d0 * (work(j,kp1) + work(j,km1) )
     *           )
            !
         end do
      end do

      return
      end    ! END: apply_diss2D

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  smooth_nearGW:                                                            cc
cc                                                                            cc
cc     NB: ghostwidth aware.                                                  cc
cc                                                                            cc
cc     Smooth points next to boundary                                         cc
cc     Now checks chr() function to ensure boundary is an amr_bdy,            cc
cc     (see chr.inc).                                                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smooth_nearGW(   field, chr, nx,   ny,  nz, gw)
      implicit    none
      integer     nx,ny,nz, gw
      real(kind=8)      field(nx,ny,nz), chr(nx,ny,nz)
      integer     i,j,k
      real(kind=8)      temp0, temp2
      real(kind=8)      interp_cubicshort
      external          interp_cubicshort
      include    'chr.inc'
      !
      ! Use linear interpolation?
      !
      logical     linear

      logical     ltrace
      parameter ( ltrace = .false. )


      if (gw.lt.2.or.nx.lt.5.or.ny.lt.5.or.nz.lt.5) then
         !
         ! Linear smoothing:
         !
         linear = .true.
         temp0  = 0.5d0
         temp2  = 0.5d0
      else
         linear = .false.
      end if

      if (ltrace) then
         write(*,*) 'smooth_nearGW: using linear interpolation:',linear
         write(*,*) 'smooth_nearGW: gw                        =',gw
      end if

      if ( (nx.lt.3).or.(ny.lt.3).or.(nz.lt.3) ) then
         write(*,*) 'smooth_nearGW: Grid too small to be smoothed'
         write(*,*) 'smooth_nearGW: nx,ny,nz = ',nx,ny,nz
         write(*,*) 'smooth_nearGW: Quitting...'
         stop
      end if

      
      k = nz-gw
         do j = 1+gw, ny-gw
            do i = 1+gw, nx-gw
                if (       chr(i,j,nz  ) .eq. CHR_amr_bdy
     *               .and. chr(i,j,k)    .ne. CHR_deco_bdy ) then
                if (linear) then
                field(i,j,k) = temp0*field(i,j,k+1)+temp2*field(i,j,k-1)
                else
                field(i,j,k) = interp_cubicshort(  field(i,j,k+2),
     *                                             field(i,j,k+1),
     *                                             field(i,j,k-1),
     *                                             field(i,j,k-2)  )
                end if
                end if
            end do
         end do
      k = 1+gw
         do j = 1+gw, ny-gw
            do i = 1+gw, nx-gw
                if (       chr(i,j,1   ) .eq. CHR_amr_bdy
     *               .and. chr(i,j,k)    .ne. CHR_deco_bdy ) then
                if (linear) then
                field(i,j,k) = temp0*field(i,j,k+1)+temp2*field(i,j,k-1)
                else
                field(i,j,k) = interp_cubicshort(  field(i,j,k+2),
     *                                             field(i,j,k+1),
     *                                             field(i,j,k-1),
     *                                             field(i,j,k-2)  )
                end if
                end if
            end do
         end do
      j = ny-gw
         do k = 1+gw, nz-gw
            do i = 1+gw, nx-gw
                if (       chr(i,ny  ,k) .eq. CHR_amr_bdy
     *               .and. chr(i,j,k)    .ne. CHR_deco_bdy ) then
                if (linear) then
                field(i,j,k) = temp0*field(i,j+1,k)+temp2*field(i,j-1,k)
                else
                field(i,j,k) = interp_cubicshort(  field(i,j+2,k),
     *                                             field(i,j+1,k),
     *                                             field(i,j-1,k),
     *                                             field(i,j-2,k)  )
                end if
                end if
            end do
         end do
      j = 1+gw
         do k = 1+gw, nz-gw
            do i = 1+gw, nx-gw
                if (       chr(i,1  ,k) .eq. CHR_amr_bdy
     *               .and. chr(i,j,k)    .ne. CHR_deco_bdy ) then
                if (linear) then
                field(i,j,k) = temp0*field(i,j+1,k)+temp2*field(i,j-1,k)
                else
                field(i,j,k) = interp_cubicshort(  field(i,j+2,k),
     *                                             field(i,j+1,k),
     *                                             field(i,j-1,k),
     *                                             field(i,j-2,k)  )
                end if
                end if
            end do
         end do
      i = nx-gw
         do k = 1+gw, nz-gw
            do j = 1+gw, ny-gw
                if (       chr(nx  ,j,k) .eq. CHR_amr_bdy
     *               .and. chr(i,j,k)    .ne. CHR_deco_bdy ) then
                if (linear) then
                field(i,j,k) = temp0*field(i+1,j,k)+temp2*field(i-1,j,k)
                else
                field(i,j,k) = interp_cubicshort(  field(i+2,j,k),
     *                                             field(i+1,j,k),
     *                                             field(i-1,j,k),
     *                                             field(i-2,j,k)  )
                end if
                end if
            end do
         end do
      i = 1+gw
         do k = 1+gw, nz-gw
            do j = 1+gw, ny-gw
                if (       chr(1   ,j,k) .eq. CHR_amr_bdy
     *               .and. chr(i,j,k)    .ne. CHR_deco_bdy ) then
                if (linear) then
                field(i,j,k) = temp0*field(i+1,j,k)+temp2*field(i-1,j,k)
                else
                field(i,j,k) = interp_cubicshort(  field(i+2,j,k),
     *                                             field(i+1,j,k),
     *                                             field(i-1,j,k),
     *                                             field(i-2,j,k)  )
                end if
                end if
            end do
         end do

      return
      end    ! END: smooth_nearGW

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  smooth_near:                                                              cc
cc                                                                            cc
cc     Smooth points next to boundary                                         cc
cc     Now checks chr() function to ensure boundary is an amr_bdy,            cc
cc     (see chr.inc).                                                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smooth_near(   field, chr, nx,   ny,  nz)
      implicit    none
      integer     nx,ny,nz
      real(kind=8)      field(nx,ny,nz), chr(nx,ny,nz)
      real(kind=8)      temp0, temp2
      integer     i,j,k
      include    'chr.inc'
      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'smooth_near:'
      end if

      if ( (nx.lt.3).or.(ny.lt.3).or.(nz.lt.3) ) then
         write(*,*) 'smooth_near: Grid too small to be smoothed'
         write(*,*) 'smooth_near: nx,ny,nz = ',nx,ny,nz
         write(*,*) 'smooth_near: Quitting...'
         stop
      end if
      !
      ! Quadratic fit:
      !
      temp0 = 0.75d0
      temp2 = 0.25d0
      !
      ! Linear smoothing:
      !
      temp0 = 0.5d0
      temp2 = 0.5d0

      k = nz-1
         do j = 2, ny-1
            do i = 2, nx-1
         if (chr(i,j,k+1).eq.CHR_amr_bdy.and.chr(i,j,k).ne.CHR_deco_bdy)
     *          field(i,j,k) = temp0*field(i,j,k+1)+temp2*field(i,j,k-1)
     *          ! Not implemented yet
c    *          field(i,j,k) = interp_quartic( field(i,j,k+1)
c    *                                         field(i,j,k-1)
c    *                                         field(i,j,k-2)
c    *                                         2.d0,
c    *                                         0.d0,
c    *                                         1.d0,
c    *                                         3.d0,           )
            end do
         end do
      k = 2
         do j = 2, ny-1
            do i = 2, nx-1
         if (chr(i,j,k-1).eq.CHR_amr_bdy.and.chr(i,j,k).ne.CHR_deco_bdy)
     #         field(i,j,k) = temp0*field(i,j,k+1)+temp2*field(i,j,k-1)
            end do
         end do
      j = ny-1
         do k = 2, nz-1
            do i = 2, nx-1
         if (chr(i,j+1,k).eq.CHR_amr_bdy.and.chr(i,j,k).ne.CHR_deco_bdy)
     #         field(i,j,k) = temp0*field(i,j+1,k)+temp2*field(i,j-1,k)
            end do
         end do
      j = 2
         do k = 2, nz-1
            do i = 2, nx-1
         if (chr(i,j-1,k).eq.CHR_amr_bdy.and.chr(i,j,k).ne.CHR_deco_bdy)
     #         field(i,j,k) = temp0*field(i,j+1,k)+temp2*field(i,j-1,k)
            end do
         end do
      i = nx-1
         do k = 2, nz-1
            do j = 2, ny-1
         if (chr(i+1,j,k).eq.CHR_amr_bdy.and.chr(i,j,k).ne.CHR_deco_bdy)
     #         field(i,j,k) = temp0*field(i+1,j,k)+temp2*field(i-1,j,k)
            end do
         end do
      i = 2
         do k = 2, nz-1
            do j = 2, ny-1
         if (chr(i-1,j,k).eq.CHR_amr_bdy.and.chr(i,j,k).ne.CHR_deco_bdy)
     #         field(i,j,k) = temp0*field(i+1,j,k)+temp2*field(i-1,j,k)
            end do
         end do

      return
      end    ! END: smooth_near

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  inc_mask_overlapB:                                                        cc
cc                  Increase by one the mask in region of overlap.            cc
cc                  Applies to the parent in region where grid overlaps.      cc
cc             NB:  Can't be done in place in replace_part_field_avgGWa()     cc
cc                  because that is called for each field on a grid.          cc
cc                  Hence, loop structure mirrors that routine.               cc
cc             NB:  decrements mask to use negative numbers so that           cc
cc                  positive numbers still reflect *type* of point as         cc
cc                  described in chr.inc.                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inc_mask_overlapB( mask, fieldc,
     *                             ax,    ay,        az,
     *                             nxp,   nyp,       nzp,
     *                             nxs,   nys,       nzs,
     *                             nxc,   nyc,       nzc,
     *                             min_ic,max_ic,
     *                             min_jc,max_jc,
     *                             min_kc,max_kc,
     *                             r, gw   )
      implicit    none
      integer     nxp,nyp,nzp, ax,ay,az, gw,r,
     *            nxc,nyc,nzc, nxs,nys,nzs
      integer     min_ip,min_jp,min_kp
      integer     min_ic,max_ic, min_jc,max_jc,min_kc,max_kc
      real(kind=8) mask(nxp,nyp,nzp), fieldc(nxs,nys,nzs)
      integer     i,j,k,   ip,jp,kp, ic,jc,kc, gwc
      integer     gwc_i, gwc_j, gwc_k
      integer     ib,ie, jb,je, kb,ke, n
      include       'largesmall.inc'

      logical     ltrace
      parameter ( ltrace = .false. )

      !gwc  = (gw + 1)/r

      if (ltrace) write(*,*) 'inc_mask_overlapB: Enter'
      do kc = min_kc, max_kc, r
         kp = (kc - min_kc) / r + 1
         k  = kp + az - 1
         do jc = min_jc, max_jc, r
            jp = (jc - min_jc) / r + 1
            j  = jp + ay - 1
            do ic = min_ic, max_ic, r
               ip = (ic - min_ic) / r + 1
               i  = ip + ax - 1
               !
               n  = - NINT(mask(i,j,k))
               if (n .ge. 0 .and. fieldc(ip,jp,kp).lt.LARGE_M_SMALL)then
                  !
                  mask(i,j,k) = mask(i,j,k) - 1.d0
                  !
               end if
            end do
         end do
      end do
      if (ltrace) write(*,*) 'inc_mask_overlapB: Done'

      return
      end        ! END: inc_mask_overlapB

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  inc_mask_overlapshadow:                                                   cc
cc                  Increase by one the mask in region of overlap.            cc
cc                  Uses negative numbers in regions where the TRE has        cc
cc                  already been computed. Parallels logic for injection.     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inc_mask_overlapshadow( mask, field, error,
     *                             mini, minj, mink,
     *                             maxi, maxj, maxk,
     *                             nx,   ny,   nz     )
      implicit     none
      integer      mini, minj, mink
      integer      maxi, maxj, maxk
      integer      nx,   ny,   nz
      real(kind=8) mask(nx,ny,nz), field(nx,ny,nz), error(nx,ny,nz)
      integer      i,j,k
      integer      ichr
      include     'chr.inc'

      logical     ltrace
      parameter ( ltrace = .false. )


      if (ltrace) write(*,*) 'inc_mask_overlapshadow: Enter'
      do k = mink, maxk
         do j = minj, maxj
            do i = mini, maxi
               !
               ichr = NINT( mask(i,j,k) ) 
               ! All values in chr.inc are strictly positive, so as long
               ! as we do not have one of those points, computer TRE:
               if (ichr .le. 0) then
               !if (ichr .ne. CHR_deco_bdy .and. ichr.ne.CHR_amr_bdy)then
                  mask(i,j,k) = mask(i,j,k) - 1.d0
                  ichr = -NINT( mask(i,j,k) ) 
                  if (ichr .gt. 1) then
                     !write(*,*) i,j,k,' ichr: ',ichr
                     error(i,j,k) = 0.5d0 * error(i,j,k) 
                  end if
               end if
               !
            end do
         end do
      end do
      if (ltrace) write(*,*) 'inc_mask_overlapshadow: Done'


      return
      end        ! END: inc_mask_overlapshadow

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  replace_part_field_avgGWc:                                                cc
cc                  NB: Replace points even where DD boundary.                cc
cc                  NB: Now ghostwidth aware.                                 cc
cc                  Replace part of a field with entirity of another.         cc
cc                  Ignores boundaries!                                       cc
cc                  Averages according to the mask function so that           cc
cc                  overlapping grids can be accomadated.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine replace_part_field_avgGWc(fieldp, maskp,fieldc,
     *                             ax,    ay,        az,
     *                             nxp,   nyp,       nzp,
     *                             nxs,   nys,       nzs,
     *                             nxc,   nyc,       nzc,
     *                             min_ic,max_ic,
     *                             min_jc,max_jc,
     *                             min_kc,max_kc,
     *                             r, gw, bwidth   )
      implicit    none
      integer     nxp,nyp,nzp, nxc,nyc,nzc, ax,ay,az, gw,r,
     *            nxs,nys,nzs, bwidth
      integer     min_ip,min_jp,min_kp
      integer     min_ic,max_ic, min_jc,max_jc,min_kc,max_kc
      real(kind=8)      fieldp(nxp,nyp,nzp), maskp(nxp,nyp,nzp),
     *            fieldc(nxs,nys,nzs)
      include    'chr.inc'
      integer     i,j,k,   ip,jp,kp, ic,jc,kc, gwc
      integer     gwc_i, gwc_j, gwc_k
      integer     ib,ie, jb,je, kb,ke, n
      real(kind=8)      child_pt
      include       'largesmall.inc'
      logical           ltrace
      parameter       ( ltrace      = .false. )
      logical           ltrace2
      parameter       ( ltrace2     = .false. )


      do kc = min_kc, max_kc, r
         kp = (kc - min_kc) / r + 1
         k  = kp + az - 1
         do jc = min_jc, max_jc, r
            jp = (jc - min_jc) / r + 1
            j  = jp + ay - 1
            do ic = min_ic, max_ic, r
               ip = (ic - min_ic) / r + 1
               i  = ip + ax - 1
               if (ltrace2) then
                  if (i.ge.nxp) then
                     write(*,*) 'ic,ip,i: ',ic, ip, i
                  end if
               end if
               !
               n        = - NINT(maskp(i,j,k))
               child_pt = fieldc(ip,jp,kp)
               if ( child_pt .lt. LARGE_M_SMALL ) then
                  if (n .ge. 0) then
                     !
                     ! Regular interior point:
                     !
                     fieldp(i,j,k) =
     *                         (     n * fieldp(i,j,k)
     *                           +       child_pt
     *                          ) / ( 1.d0 + n)
                  else if (abs(n).ne.CHR_amr_bdy) then
                     !
                     ! These points represent points either
                     ! on the coarse grid boundary or
                     ! on a   domain deco boundary
                     ! ...no averaging here:
                     !
                     fieldp(i,j,k) = child_pt
                     !
                  end if
               end if
            end do
         end do
      end do

      if (ltrace) then
         write(*,*) 'replace_part_field_avgGWc: ax/y/z =',ax,ay,az
c        write(*,*) 'replace_part_field_avgGWc: i/j/kp =',ip,jp,kp
         write(*,*) 'replace_part_field_avgGWc: n[xyz]p=',nxp,nyp,nzp
         write(*,*) 'replace_part_field_avgGWc: n[xyz]c=',nxc,nyc,nzc
         write(*,*) 'replace_part_field_avgGWc: n[xyz]s=',nxs,nys,nzs
         write(*,*) 'replace_part_field_avgGWc:mini/j/kc',
     *                                           min_ic,min_jc,min_kc
         write(*,*) 'replace_part_field_avgGWc:maxi/j/kc',
     *                                           max_ic,max_jc,max_kc
c        write(*,*) 'replace_part_field_avgGWc:    gw  =',gw
         write(*,*) 'replace_part_field_avgGWc:    r   =',r
c        write(*,*) 'replace_part_field_avgGWc:    gwc =',gwc
c        write(*,*)
c    * 'replace_part_field_avgGWc:gwci/j/k=',gwc_i,gwc_j,gwc_k
      end if


      return
      end        ! END: subroutine replace_part_field_avgGWc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  restrict_fieldGWb:                                                        cc
cc                  Restrict a Child grid to a Parent grid.                   cc
cc               NB: Do away with ghostwidth. Just get all data possibly      cc
cc                   necessary and let replace_part_field_avgGWa() figure     cc
cc                   out what to do with it and what to ignore.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine restrict_fieldGWb( fieldp,  fieldc, maskc,
     *                             nxp,   nyp,       nzp,
     *                             nxc,   nyc,       nzc,
     *                             min_i, max_i,
     *                             min_j, max_j,
     *                             min_k, max_k,
     *                             r, gw)
      implicit    none
      integer     nxp,nyp,nzp, nxc,nyc,nzc, r, gw
      real(kind=8)      fieldp(nxp,nyp,nzp), fieldc(nxc,nyc,nzc),
     *            maskc(nxc,nyc,nzc)
      include    'chr.inc'
      integer     ic,jc,kc, ip,jp,kp, gwc, min_i,max_i,
     *                    min_j,max_j,min_k,max_k
      real(kind=8)      temp1, temp2
      include       'largesmall.inc'
      logical           HALFWEIGHTED
      parameter       ( HALFWEIGHTED = .false. )
      integer     AMR_BDY

      logical     ltrace
      parameter ( ltrace = .false. )


      AMR_BDY = NINT(CHR_amr_bdy)

      if (ltrace) then
         write(*,*) 'restrict_fieldGWb: Start....'
         call field_dump_info(fieldp,nxp,nyp,nzp)
      end if

      do kc = min_k, max_k, r
         kp = (kc - min_k) / r + 1
         do jc = min_j, max_j, r
            jp = (jc - min_j) / r + 1
            do ic = min_i, max_i, r
               ip = (ic - min_i) / r + 1
               !
               ! Determine if this data should get injected to parent:
               !
               if (
     *        (ic.le.gw      .and.NINT(maskc(1,jc,kc)  ).eq. AMR_BDY)
     *    .or.(jc.le.gw      .and.NINT(maskc(ic,1,kc)  ).eq. AMR_BDY)
     *    .or.(kc.le.gw      .and.NINT(maskc(ic,jc,1)  ).eq. AMR_BDY)
     *    .or.(ic.ge.nxc-gw+1.and.NINT(maskc(nxc,jc,kc)).eq. AMR_BDY)
     *    .or.(jc.ge.nyc-gw+1.and.NINT(maskc(ic,nyc,kc)).eq. AMR_BDY)
     *    .or.(kc.ge.nzc-gw+1.and.NINT(maskc(ic,jc,nzc)).eq. AMR_BDY)
     *      )  then
                  fieldp(ip,jp,kp) = LARGENUMBER
               else
                  if (HALFWEIGHTED.and.
     *                ic.gt.1  .and.jc.gt.1  .and.kc.gt.1  .and.
     *                ic.lt.nxc.and.jc.lt.nyc.and.kc.lt.nzc ) then
                     fieldp(ip,jp,kp) = 0.5d0 * fieldc(ic,  jc,  kc  )
     *                                 +(1.d0/12.d0)*(
     *                                        + fieldc(ic+1,jc,  kc  )
     *                                        + fieldc(ic-1,jc,  kc  )
     *                                        + fieldc(ic,  jc+1,kc  )
     *                                        + fieldc(ic,  jc-1,kc  )
     *                                        + fieldc(ic,  jc,  kc+1)
     *                                        + fieldc(ic,  jc,  kc-1)
     *                                                )
                  else
                     fieldp(ip,jp,kp) = fieldc(ic,  jc,  kc  )
                  end if
               end if
            end do
         end do
      end do

      if (ltrace) then
         call field_dump_info(fieldp,nxp,nyp,nzp)
         write(*,*) 'restrict_fieldGWb: i/j/kp =',ip,jp,kp
         write(*,*) 'restrict_fieldGWb: '
         write(*,*) 'restrict_fieldGWb: r   = ', r
         write(*,*) 'restrict_fieldGWb: gw  = ', gw
         write(*,*) 'restrict_fieldGWb: ',1+gwc*r,nzc-gwc*r
         write(*,*) 'restrict_fieldGWb: ip,jp,kp = ',ip,jp,kp
         write(*,*) 'restrict_fieldGWb: nx/y/zp = ',nxp,nyp,nzp
         write(*,*) 'restrict_fieldGWb: nx/y/zc = ',nxc,nyc,nzc
         write(*,*) 'restrict_fieldGWb: i/j/kc  = ',ic,jc,kc
         write(*,*) 'restrict_fieldGWb: min_i/j/k ',min_i,min_j,min_k
         write(*,*) 'restrict_fieldGWb: max_i/j/k ',max_i,max_j,max_k
      end if

      return
      end    ! END: restrict_fieldGWb


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  extrap_n_np1:                                                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine extrap_n_np1( field_np1, field_n, field_nm1,
     *                         nx, ny, nz)
      implicit    none
      integer     nx,ny,nz
      real(kind=8)      field_np1(nx,ny,nz), field_n(nx,ny,nz),
     *            field_nm1(nx,ny,nz)
      integer     i,j,k
      logical     ltrace
      parameter ( ltrace = .false. )

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               field_np1(i,j,k) = 2.d0*field_n(i,j,k) - field_nm1(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: extrap_n_np1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  int_bounds_time:                                                          cc
cc                                                                            cc
cc      Interpolates in time between the _np1 and _n levels some region       cc
cc      of a grid into a storage array.                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine int_bounds_time( field, field_np1, field_n,
     *                              t,     t_np1,     t_n,
     *                              mini,  minj,      mink,
     *                              maxi,  maxj,      maxk,
     *                              nxs,   nys,       nzs,
     *                              nx,    ny,        nz     )
      implicit          none
      integer           nx,ny,nz, nxs,nys,nzs
      integer           mini, minj, mink, maxi, maxj, maxk
      real(kind=8)      field(nxs,nys,nzs),
     *                  field_np1(nx,ny,nz),
     *                  field_n(nx,ny,nz),
     *                  t, t_np1, t_n
      integer           i,j,k, is, js, ks
      real(kind=8)      alpha_n, alpha_np1
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )


      if (t_np1 .eq. t_n) then
         write(*,*) 'int_bounds_time: Problem!'
         stop
      end if

      alpha_n   = ( t_np1 - t )/( t_np1 - t_n )
      alpha_np1 = 1.d0 - alpha_n

      if (ltrace) then
         write(*,*) 'int_bounds_time: t/t_np1/t_n: ',t,t_np1,t_n
         write(*,*) 'int_bounds_time: alpha_n/np1: ',alpha_n,alpha_np1
      end if

      do k = mink, maxk
         ks = (k-mink) + 1
         do j = minj, maxj
            js = (j-minj) + 1
            do i = mini, maxi
               is = (i-mini) + 1
               field(is,js,ks) =   alpha_np1 * field_np1(i,j,k)
     *                           + alpha_n   *   field_n(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: int_bounds_time


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  int_field_bounds_timeGWa:                                                 cc
cc                                                                            cc
cc      NB: Can handle multiple parents.                                      cc
cc      NB: modified from int_field_bounds() to interpolate within            cc
cc      a ghost region of width gw.                                           cc
cc                                                                            cc
cc      interpolate (in time) boundary values from parent grid                cc
cc      to a temporary array, stored linearly.                                cc
cc      Takes general refinement factor, r, on input.                         cc
cc      Does *not* do any spatial interpolation. Just interpolate             cc
cc      the boundary in time. Then this information can passed to a finer     cc
cc      grid and interpolated in space.                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine int_field_bounds_timeGWa(  field, field_np1, field_n,
     *                              t,     t_np1,     t_n,
     *                              mini,  minj,      mink,
     *                              maxi,  maxj,      maxk,
     *                              nx,    ny,        nz,
     *                              length, r,        gw)
      implicit    none
      integer     nx,ny,nz, r, length,gw
      integer     mini, minj, mink, maxi, maxj, maxk
      real(kind=8)      field(length),
     *            field_np1(nx,ny,nz),
     *            field_n(nx,ny,nz),
     *            t, t_np1, t_n
      integer     ic,jc,kc, index, gwc
      real(kind=8)      alpha_n, alpha_np1
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )


      if (t_np1.eq.t_n) then
         write(*,*) 'int_field_bounds_timeGWa: Problem!!!!'
         write(*,*) 'int_field_bounds_timeGWa: t_np1 = t_n'
         write(*,*) 'int_field_bounds_timeGWa:     t_np1 =', t_np1
         write(*,*) 'int_field_bounds_timeGWa:     t_n   =', t_n
         write(*,*) 'int_field_bounds_timeGWa:     t     =', t  
         write(*,*) 'int_field_bounds_timeGWa: Quitting...'
         stop
      end if

      alpha_n   = ( t_np1 - t )/( t_np1 - t_n )
      alpha_np1 = 1.d0 - alpha_n

      ! Find the number of points to interpolate..remember
      ! gw = 1   r = 2  --> gwc = 1
      ! gw = 2   r = 2  --> gwc = 2
      ! gw = 3   r = 2  --> gwc = 2
      ! gw = 4   r = 2  --> gwc = 3
      ! gw = 5   r = 2  --> gwc = 3
      ! gw = 6   r = 2  --> gwc = 4
      ! For general value of r, then we need:
      !
      gwc  = gw / r + 1

      !
      ! interpolate  points in space on the boundary points of parent:
      !
      index = 0

c     if (mink .gt. 1 .and. maxk-mink+1.ge.gwc) then
         if (ltrace) write(*,*)'int_field_bounds_timeGWa: Lower z bound'
         do kc = mink, mink+gwc-1
            do jc = minj, maxj
               do ic = mini, maxi
                  index = index + 1
                  field(index   ) =   alpha_np1 * field_np1(ic,jc,kc)
     *                              + alpha_n   * field_n(ic,jc,kc)
                  if (ltrace2) write(*,*) ic,jc,kc,index
               end do
            end do
         end do
c     end if
      if (ltrace)write(*,*)'int_field_bounds_timeGWa: index     =',index
      !
c     if (maxk .lt. nz .and. maxk-mink+1.ge.gwc) then
         if (ltrace) write(*,*)'int_field_bounds_timeGWa: Upper z bound'
         do kc = maxk, maxk-gwc+1, -1
            do jc = minj, maxj
               do ic = mini, maxi
                  index = index + 1
                  field(index   ) =   alpha_np1 * field_np1(ic,jc,kc)
     *                              + alpha_n   * field_n(ic,jc,kc)
                  if (ltrace2) write(*,*) ic,jc,kc,index
               end do
            end do
         end do
c     end if
      if (ltrace)write(*,*)'int_field_bounds_timeGWa: index     =',index
      !
c     if (minj .gt. 1 .and. maxj-minj+1.ge.gwc) then
         if (ltrace) write(*,*)'int_field_bounds_timeGWa: Lower y bound'
         do jc = minj, minj+gwc-1
            do kc = mink, maxk
               do ic = mini, maxi
                  index = index + 1
                  field(index   ) =   alpha_np1 * field_np1(ic,jc,kc)
     *                              + alpha_n   * field_n(ic,jc,kc)
                  if (ltrace2) write(*,*) ic,jc,kc,index
               end do
            end do
         end do
c     end if
      if (ltrace)write(*,*)'int_field_bounds_timeGWa: index     =',index
      !
c     if (maxj .lt. ny .and. maxj-minj+1.ge.gwc) then
         if (ltrace) write(*,*)'int_field_bounds_timeGWa: Upper y bound'
         do jc = maxj, maxj-gwc+1, -1
            do kc = mink, maxk
               do ic = mini, maxi
                  index = index + 1
                  field(index   ) =   alpha_np1 * field_np1(ic,jc,kc)
     *                              + alpha_n   * field_n(ic,jc,kc)
                  if (ltrace2) write(*,*) ic,jc,kc,index
               end do
            end do
         end do
c     end if
      if (ltrace)write(*,*)'int_field_bounds_timeGWa: index     =',index
      !
c     if (mini .gt. 1 .and. maxi-mini+1.ge.gwc) then
         if (ltrace) write(*,*)'int_field_bounds_timeGWa: Lower x bound'
         do ic = mini, mini+gwc-1
            do kc = mink, maxk
               do jc = minj, maxj
                  index = index + 1
                  field(index   ) =   alpha_np1 * field_np1(ic,jc,kc)
     *                              + alpha_n   * field_n(ic,jc,kc)
                  if (ltrace2) write(*,*) ic,jc,kc,index
               end do
            end do
         end do
c     end if
      if (ltrace)write(*,*)'int_field_bounds_timeGWa: index     =',index
      !
c     if (maxi .lt. nx .and. maxi-mini+1.ge.gwc) then
         if (ltrace)write(*,*) 'int_field_bounds_timeGWa: Upper x bound'
         do ic = maxi, maxi-gwc+1, -1
            do kc = mink, maxk
               do jc = minj, maxj
                  index = index + 1
                  field(index   ) =   alpha_np1 * field_np1(ic,jc,kc)
     *                              + alpha_n   * field_n(ic,jc,kc)
                  if (ltrace2) write(*,*) ic,jc,kc,index
               end do
            end do
         end do
c     end if
      if (ltrace)write(*,*)'int_field_bounds_timeGWa: index     =',index


      if (ltrace) then
         write(*,*) 'int_field_bounds_timeGWa:     t_np1 =', t_np1
         write(*,*) 'int_field_bounds_timeGWa:     t_n   =', t_n
         write(*,*) 'int_field_bounds_timeGWa:     t     =', t
         write(*,*) 'int_field_bounds_timeGWa:     nx    =', nx
         write(*,*) 'int_field_bounds_timeGWa:     ny    =', ny
         write(*,*) 'int_field_bounds_timeGWa:     nz    =', nz
         write(*,*) 'int_field_bounds_timeGWa:     r     =', r  
         write(*,*) 'int_field_bounds_timeGWa: alpha_np1 =', alpha_np1
         write(*,*) 'int_field_bounds_timeGWa: alpha_n   =', alpha_n
         write(*,*) 'int_field_bounds_timeGWa:     gw,gwc =',gw,gwc
         write(*,*) 'int_field_bounds_timeGWa:mini/j/k =',mini,minj,mink
         write(*,*) 'int_field_bounds_timeGWa:maxi/j/k =',maxi,maxj,maxk
         write(*,*) 'int_field_bounds_timeGWa: length     =',length
         write(*,*) 'int_field_bounds_timeGWa: index      =',index 
         write(*,*) ''
      end if

      if (index .gt. length) then
         write(*,*) 'int_field_bounds_timeGWa: length     =',length
         write(*,*) 'int_field_bounds_timeGWa: index      =',index 
         write(*,*) 'int_field_bounds_timeGWa: Problem with storage'
         stop
      end if

      return
      end    ! END: int_field_bounds_timeGWa

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  int_field_bounds_spaceGWa:                                                cc
cc                                                                            cc
cc      NB: modified from int_field_bounds() to interpolate within            cc
cc      a ghost region of width gw.                                           cc
cc                                                                            cc
cc      interpolate (in space) along the boundaries of a field                cc
cc      by inserting values passed into the routine and then interpolating.   cc
cc      Takes general refinement factor, r, on input.                         cc
cc      Unpacking of "values" must match packing as found in                  cc
cc      "int_field_bounds_time".                                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine int_field_bounds_spaceGWa(  values, field, chr,
     *                              mini,  minj,      mink,
     *                              maxi,  maxj,      maxk,
     *                              nx,    ny,        nz,
     *                              length, r, gw)
      implicit    none
      integer     nx,ny,nz, r, length, gw
      integer     mini,minj,mink, maxi,maxj,maxk
      real(kind=8)      values(length),
     *            field(nx,ny,nz), chr(nx,ny,nz)
      include    'chr.inc'
      integer     ic,jc,kc, i, if,jf,kf, index, gwc
      integer     itmp, tmpgw, indexL, indexR
      integer     nxcp, nycp, nzcp
      !
      ! need to keep track of the "index":
      !
      integer     LOWER, UPPER, FRONT, BACK, LEFT
      real(kind=8)      temp, tmp
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      !
      ! Initialize to zero:
      !
      LOWER  = 0
      UPPER  = 0
      FRONT  = 0
      BACK   = 0
      LEFT   = 0

      !
      ! If gw is even, then we have to loop up to gw+1:
      !   For example, if r=2 and g=6, then we need to loop
      !         over kf = 1, 3, 5, 7
      !   *not* over kf = 1, 3, 5
      ! For general value of r, then we need:
      !
      gwc   = gw / r + 1
      itmp  = mod(gw-1,r)
      tmpgw = gw + itmp

      nxcp = (maxi-mini)/r + 1
      nycp = (maxj-minj)/r + 1
      nzcp = (maxk-mink)/r + 1

c     copy  points on the boundary from passed in values
c       and then interpolate over grid between child points:
c
c       NB: the values passed in may have more data than we
c           want to copy in directly. Instead, they're used
c           just for interpolating (this is the case when gw is odd)
      !-----------
      ! Lower Face
      !-----------
      LOWER = gwc * nxcp * nycp
      UPPER = gwc * nxcp * nycp
      FRONT = gwc * nxcp * nzcp
      LEFT  = gwc * nycp * nzcp
      BACK  = gwc * nxcp * nzcp
      if (nzcp .lt. gwc) goto 20
      if (mink .ne. 1  ) goto 10
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: Lower face'
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: index=',index
      do kf = 1, gw, r
         kc = (kf - 1)/r + 1
         !
         ! Copy values directly:
         !
         do jf = minj, maxj, r
            jc = (jf - minj)/r + 1
            do if = mini, maxi, r
               ic = (if - mini)/r + 1
               index = ic + (jc-1)*nxcp + (kc-1)*nxcp*nycp
               if ( chr(if,jf,1 ).ne.CHR_deco_bdy)
     *         field(if,jf,kf) =   values(index)
               if (ltrace2) write(*,*) ic,jc,kc,index
            end do
         end do
         !
         ! Interpolate in z direction:
         !
         if (kf.lt.gw) then
         do jf = minj, maxj, r
            jc = (jf - minj)/r + 1
            do if = mini, maxi, r
               ic = (if - mini)/r + 1
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  indexL = ic + (jc-1)*nxcp + (kc-1)*nxcp*nycp
                  indexR = ic + (jc-1)*nxcp + (kc  )*nxcp*nycp
                  if ( chr(if,jf,1 ).ne.CHR_deco_bdy 
     *                .and. kf+i .le. gw)
     *            field(if,jf,kf+i) =   values(indexL)*(1.d0-temp)
     *                                + values(indexR)*temp
               end do
            end do
         end do
         end if
      end do
      do kf = 1, gw, 1
         !
         ! Interpolate in x direction:
         !
         do jf = minj, maxj, r
            do if = mini, maxi-r, r
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(if,jf,1 ).ne.CHR_deco_bdy )
     *            field(if+i,jf,kf) =   field(if,  jf,kf)*(1.d0-temp)
     *                                + field(if+r,jf,kf)*temp
               end do
            end do
         end do
         !
         ! Interpolate in y direction:
         !
         do jf = minj, maxj-r, r
            do if = mini, maxi
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(if,jf,1 ).ne.CHR_deco_bdy )
     *            field(if,jf+i,kf) =   field(if,jf  ,kf)* (1.d0-temp)
     *                                + field(if,jf+r,kf)* temp
               end do
            end do
         end do
      end do

 10   continue

      !-----------
      ! Upper face:
      !-----------
      if (maxk .ne. nz) goto 20
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: Upper face'
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: index=',LOWER
      tmp = (nz-1)/r
      do kf = nz, nz-gw, -r
         kc = (kf - 1)/r + 1
         !
         ! Copy values directly:
         !
         if (kf.ge.nz-gw+1) then
         do jf = minj, maxj, r
            jc = (jf - minj)/r + 1
            do if = mini, maxi, r
               ic = (if - mini)/r + 1
               index = ic + (jc-1)*nxcp - (kc-1-tmp)*nxcp*nycp
     *                 + LOWER
               if ( chr(if,jf,nz).ne.CHR_deco_bdy )
     *         field(if,jf,kf) =   values(index)
               if (ltrace2) write(*,*) ic,jc,kc,index
            end do
         end do
         end if
         !
         ! Interpolate in z direction:
         !
         if (kf.lt.nz) then
         do jf = minj, maxj, r
            jc = (jf - minj)/r + 1
            do if = mini, maxi, r
               ic = (if - mini)/r + 1
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  indexL = ic + (jc-1)*nxcp - (kc-1-tmp)*nxcp*nycp
     *                 + LOWER
                  indexR = ic + (jc-1)*nxcp - (kc  -tmp)*nxcp*nycp
     *                 + LOWER
                  if ( chr(if,jf,nz).ne.CHR_deco_bdy 
     *                .and. kf+i .ge. nz-gw+1)
     *            field(if,jf,kf+i) =   values(indexL)*(1.d0-temp)
     *                                + values(indexR)*temp
               end do
            end do
         end do
         end if
      end do
      do kf = nz, nz-gw+1, -1
         !
         ! Interpolate in x direction:
         !
         do jf = minj, maxj, r
            do if = mini, maxi-r, r
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(if,jf,nz).ne.CHR_deco_bdy )
     *            field(if+i,jf,kf) =   field(if,  jf,kf)*(1.d0-temp)
     *                                + field(if+r,jf,kf)*temp
               end do
            end do
         end do
         !
         ! Interpolate in y direction:
         !
         do jf = minj, maxj-r, r
            do if = mini, maxi
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(if,jf,nz).ne.CHR_deco_bdy )
     *            field(if,jf+i,kf) =   field(if,jf  ,kf)* (1.d0-temp)
     *                                + field(if,jf+r,kf)* temp
               end do
            end do
         end do
      end do

 20   continue

      !-----------
      ! Front face:
      !-----------
      if (nycp .lt. gwc) goto 40
      if (minj .ne. 1  ) goto 30
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: Front face'
      if(ltrace)write(*,*)'int_field_bounds_spaceGWa:index=',LOWER+UPPER
      do jf = 1, gw, r
         jc = (jf - 1)/r + 1
         !
         ! Copy values directly:
         !
         do kf = mink, maxk, r
            kc = (kf - mink)/r + 1
            do if = mini, maxi, r
               ic = (if - mini)/r + 1
               index = ic + (kc-1)*nxcp + (jc-1)*nxcp*nzcp
     *                 + LOWER + UPPER
               if ( chr(if,1,kf).ne.CHR_deco_bdy )
     *         field(if,jf,kf) =   values(index)
               if (ltrace2) write(*,*) ic,jc,kc,index
            end do
         end do
         !
         ! Interpolate in y direction:
         !
         if (jf.lt.gw) then
         do kf = mink, maxk, r
            kc = (kf - mink)/r + 1
            do if = mini, maxi, r
               ic = (if - mini)/r + 1
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  indexL = ic + (kc-1)*nxcp + (jc-1)*nxcp*nzcp
     *                 + LOWER + UPPER
                  indexR = ic + (kc-1)*nxcp + (jc  )*nxcp*nzcp
     *                 + LOWER + UPPER
                  if ( chr(if,1,kf).ne.CHR_deco_bdy 
     *                .and. jf+i .le. gw)
     *            field(if,jf+i,kf) =   values(indexL)*(1.d0-temp)
     *                                + values(indexR)*temp
               end do
            end do
         end do
         end if
      end do
      do jf = 1, gw, 1
         jc = (jf - 1)/r + 1
         !
         ! Interpolate in x direction:
         !
         do kf = mink, maxk, r
            do if = mini, maxi-r, r
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(if,1,kf).ne.CHR_deco_bdy )
     *            field(if+i,jf,kf) =   field(if,  jf,kf)*(1.d0-temp)
     *                                + field(if+r,jf,kf)*temp
               end do
            end do
         end do
         !
         ! Interpolate in z direction:
         !
         do kf = mink, maxk-r, r
            do if = mini, maxi
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(if,1,kf).ne.CHR_deco_bdy )
     *            field(if,jf,kf+i) =   field(if,jf,kf  )* (1.d0-temp)
     *                                + field(if,jf,kf+r)* temp
               end do
            end do
         end do
      end do

 30   continue

      !-----------
      ! Back face:
      !-----------
      if (maxj .ne. ny) goto 40
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: Back  face'
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: index=',LOWER+
     *     UPPER + FRONT
      tmp = (ny-1)/r
      do jf = ny, ny-gw, -r
         jc = (jf - 1)/r + 1
         !
         ! Copy values directly:
         !
         if (jf .ge. ny-gw+1) then
         do kf = mink, maxk, r
            kc = (kf - mink)/r + 1
            do if = mini, maxi, r
               ic = (if - mini)/r + 1
               index = ic + (kc-1)*nxcp - (jc-1-tmp)*nxcp*nzcp
     *                 + LOWER + UPPER + FRONT
               if ( chr(if,ny,kf).ne.CHR_deco_bdy )
     *         field(if,jf,kf) =   values(index)
               if (ltrace2) write(*,*) ic,jc,kc,index
            end do
         end do
         end if
         !
         ! Interpolate in y direction:
         !
         if (jf.lt.ny) then
         do kf = mink, maxk, r
            kc = (kf - mink)/r + 1
            do if = mini, maxi, r
               ic = (if - mini)/r + 1
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  indexL = ic + (kc-1)*nxcp - (jc-1-tmp)*nxcp*nzcp
     *                 + LOWER + UPPER + FRONT
                  indexR = ic + (kc-1)*nxcp - (jc  -tmp)*nxcp*nzcp
     *                 + LOWER + UPPER + FRONT
                  if ( chr(if,ny,kf).ne.CHR_deco_bdy 
     *                .and. jf+i .ge. ny-gw+1)
     *            field(if,jf+i,kf) =   values(indexL)*(1.d0-temp)
     *                                + values(indexR)*temp
               end do
            end do
         end do
         end if
      end do
      do jf = ny, ny-gw+1, -1
         !
         ! Interpolate in x direction:
         !
         do kf = mink, maxk, r
            do if = mini, maxi-r, r
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(if,ny,kf).ne.CHR_deco_bdy )
     *            field(if+i,jf,kf) =   field(if,  jf,kf)*(1.d0-temp)
     *                                + field(if+r,jf,kf)*temp
               end do
            end do
         end do
         !
         ! Interpolate in z direction:
         !
         do kf = mink, maxk-r, r
            do if = mini, maxi
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(if,ny,kf).ne.CHR_deco_bdy )
     *            field(if,jf,kf+i) =   field(if,jf,kf  )* (1.d0-temp)
     *                                + field(if,jf,kf+r)* temp
               end do
            end do
         end do
      end do

 40   continue

      !-----------
      ! Left face:
      !-----------
      if (nxcp .lt. gwc) goto 60
      if (mini .ne. 1  ) goto 50
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: Left  face'
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: index=',LOWER+
     *     UPPER + FRONT + BACK
      do if = 1, gw, r
         ic = (if - 1)/r + 1
         !
         ! Copy values directly:
         !
         do kf = mink, maxk, r
            kc = (kf - mink)/r + 1
            do jf = minj, maxj, r
               jc = (jf - minj)/r + 1
               index = jc + (kc-1)*nycp + (ic-1)*nycp*nzcp
     *                 + LOWER + UPPER + FRONT + BACK
               if ( chr(1, jf,kf).ne.CHR_deco_bdy )
     *         field(if,jf,kf) =   values(index)
               if (ltrace2) write(*,*) ic,jc,kc,index
            end do
         end do
         !
         ! Interpolate in x direction:
         !
         if (if.lt.gw) then
         do kf = mink, maxk, r
            kc = (kf - mink)/r + 1
            do jf = minj, maxj, r
               jc = (jf - minj)/r + 1
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  indexL = jc + (kc-1)*nycp + (ic-1)*nycp*nzcp
     *                 + LOWER + UPPER + FRONT + BACK
                  indexR = jc + (kc-1)*nycp + (ic  )*nycp*nzcp
     *                 + LOWER + UPPER + FRONT + BACK
                  if ( chr(1, jf,kf).ne.CHR_deco_bdy 
     *                .and. if+i .le. gw)
     *            field(if+i,jf,kf) =   values(indexL)*(1.d0-temp)
     *                                + values(indexR)*temp
               end do
            end do
         end do
         end if
      end do
      do if = 1, gw, 1
         ic = (if - 1)/r + 1
         !
         ! Interpolate in y direction:
         !
         do kf = mink, maxk, r
            do jf = minj, maxj-r, r
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(1, jf,kf).ne.CHR_deco_bdy )
     *            field(if,jf+i,kf) =   field(if,jf  ,kf)*(1.d0-temp)
     *                                + field(if,jf+r,kf)*temp
               end do
            end do
         end do
         !
         ! Interpolate in z direction:
         !
         do kf = mink, maxk-r, r
            do jf = minj, maxj
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(1, jf,kf).ne.CHR_deco_bdy )
     *            field(if,jf,kf+i) =   field(if,jf,kf  )* (1.d0-temp)
     *                                + field(if,jf,kf+r)* temp
               end do
            end do
         end do
      end do

 50   continue

      !-----------
      ! Right face:
      !-----------
      if (maxi .ne. nx) goto 60
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: Right face'
      if (ltrace) write(*,*) 'int_field_bounds_spaceGWa: index=',LOWER+
     *     UPPER + FRONT + BACK + LEFT
      tmp = (nx-1)/r
      do if = nx, nx-gw, -r
         ic = (if - 1)/r + 1
         !
         ! Copy values directly:
         !
         if (if .ge. nx-gw+1) then
         do kf = mink, maxk, r
            kc = (kf - mink)/r + 1
            do jf = minj, maxj, r
               jc = (jf - minj)/r + 1
               index = jc + (kc-1)*nycp - (ic-1-tmp)*nycp*nzcp
     *                 + LOWER + UPPER + FRONT + BACK + LEFT
               if ( chr(nx,jf,kf).ne.CHR_deco_bdy )
     *         field(if,jf,kf) =   values(index)
               if (ltrace2) write(*,*) ic,jc,kc,index
            end do
         end do
         end if
         !
         ! Interpolate in x direction:
         !
         if (if.lt.nx) then
         do kf = mink, maxk, r
            kc = (kf - mink)/r + 1
            do jf = minj, maxj, r
               jc = (jf - minj)/r + 1
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  indexL = jc + (kc-1)*nycp - (ic-1-tmp)*nycp*nzcp
     *                 + LOWER + UPPER + FRONT + BACK + LEFT
                  indexR = jc + (kc-1)*nycp - (ic  -tmp)*nycp*nzcp
     *                 + LOWER + UPPER + FRONT + BACK + LEFT
                  if ( chr(nx,jf,kf).ne.CHR_deco_bdy 
     *                .and. if+i .ge. nx-gw+1)
     *            field(if+i,jf,kf) =   values(indexL)*(1.d0-temp)
     *                                + values(indexR)*temp
               end do
            end do
         end do
         end if
      end do
      do if = nx, nx-gw+1, -1
         !
         ! Interpolate in y direction:
         !
         do kf = mink, maxk, r
            do jf = minj, maxj-r, r
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(nx,jf,kf).ne.CHR_deco_bdy )
     *            field(if,jf+i,kf) =   field(if,jf  ,kf)*(1.d0-temp)
     *                                + field(if,jf+r,kf)*temp
               end do
            end do
         end do
         !
         ! Interpolate in z direction:
         !
         do kf = mink, maxk-r, r
            do jf = minj, maxj
               do i = 1, r-1
                  temp = (1.d0*i) / r
                  if ( chr(nx,jf,kf).ne.CHR_deco_bdy )
     *            field(if,jf,kf+i) =   field(if,jf,kf  )* (1.d0-temp)
     *                                + field(if,jf,kf+r)* temp
               end do
            end do
         end do
      end do

 60   continue

      if (ltrace) then
         write(*,*) 'int_field_bounds_spaceGWa:     nx    =', nx
         write(*,*) 'int_field_bounds_spaceGWa:     ny    =', ny
         write(*,*) 'int_field_bounds_spaceGWa:     nz    =', nz
         write(*,*) 'int_field_bounds_spaceGWa:     r     =', r  
         write(*,*) 'int_field_bounds_spaceGWa:     gw    =', gw
         write(*,*) 'int_field_bounds_spaceGWa:     gwc   =', gwc
         write(*,*) 'int_field_bounds_spaceGWa:     tmpgw =', tmpgw
         write(*,*) 'int_field_bounds_spaceGWa: mini/j/k  =', 
     *                                                mini,minj,mink
         write(*,*) 'int_field_bounds_spaceGWa: maxi/j/k  =', 
     *                                                maxi,maxj,maxk
         write(*,*) 'int_field_bounds_spaceGWa:     length=', length
         write(*,*) 'int_field_bounds_spaceGWa:     index =', index
         write(*,*) 'int_field_bounds_spaceGWa:     nxcp  =', nxcp
         write(*,*) 'int_field_bounds_spaceGWa:     nycp  =', nycp
         write(*,*) 'int_field_bounds_spaceGWa:     nzcp  =', nzcp
      end if

      return
      end    ! END: int_field_bounds_spaceGWa

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  linear_ramp:                                                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine linear_ramp( field, min, h, num )
      implicit    none
      integer     num
      real(kind=8)      field(num), min, h
      integer     i

      do i = 1, num
         field(i) = min + (i-1.d0)*h
      end do

      return
      end    ! END: linear_ramp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  integrate3d:                                                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function integrate3d( field, h, nx, ny, nz )
      implicit    none
      integer     nx, ny, nz
      real(kind=8)      field(nx,ny,nz), h
      integer     i,  j,  k

      integrate3d = 0.0d0
      do k = 2, nz
         do j = 2, ny
            do i = 2, nx
               integrate3d = integrate3d
     .      + 0.125d0*(
     .           field(i  ,j  ,k  )
     .         + field(i-1,j  ,k  )
     .         + field(i  ,j-1,k  )
     .         + field(i-1,j-1,k  )
     .         + field(i  ,j  ,k-1)
     .         + field(i-1,j  ,k-1)
     .         + field(i  ,j-1,k-1)
     .         + field(i-1,j-1,k-1)
     .                               )*h**3
            end do
         end do
      end do

      return
      end    ! END: integrate3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  pintegrate3d:                                                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function pintegrate3d( field,
     *                                        istart, ifinish,
     *                                        jstart, jfinish,
     *                                        kstart, kfinish,
     *                                        h, nx, ny, nz )
      implicit    none
      integer     istart, ifinish,
     *            jstart, jfinish,
     *            kstart, kfinish,
     *            nx,     ny,      nz
      real(kind=8)      field(nx,ny,nz), h
      integer     i,  j,  k
      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) then
         write(*,*) 'pintegrated3d: istart, ifinish =',istart,ifinish
         write(*,*) 'pintegrated3d: jstart, jfinish =',jstart,jfinish
         write(*,*) 'pintegrated3d: kstart, kfinish =',kstart,kfinish
         write(*,*) 'pintegrated3d: h               =',h
         write(*,*) 'pintegrated3d: nx,ny,nz        =',nx,ny,nz
      end if

      pintegrate3d = 0.0d0
      do k = kstart+1, kfinish
         do j = jstart+1, jfinish
            do i = jstart+1, ifinish
               pintegrate3d = pintegrate3d
     .      + 0.125d0*(
     .           field(i  ,j  ,k  )
     .         + field(i-1,j  ,k  )
     .         + field(i  ,j-1,k  )
     .         + field(i-1,j-1,k  )
     .         + field(i  ,j  ,k-1)
     .         + field(i-1,j  ,k-1)
     .         + field(i  ,j-1,k-1)
     .         + field(i-1,j-1,k-1)
     .                               )*h**3
            end do
         end do
      end do

      return
      end    ! END: pintegrate3d


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mat_mat_sub3d:                                                            cc
cc                 Subtract two matricies and puts them in third.             cc
cc                   C = A + B                                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mat_mat_sub3d( field_c, field_a, field_b, nx,ny,nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field_c(nx,ny,nz), field_a(nx,ny,nz)
      real(kind=8)     field_b(nx,ny,nz)
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               field_c(i,j,k) = field_a(i,j,k) - field_b(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: mat_mat_sub3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mat_mat_add3d:                                                            cc
cc                 Adds two matricies and puts them in third.                 cc
cc                   C = A + B                                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mat_mat_add3d( field_c, field_a, field_b, nx,ny,nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field_c(nx,ny,nz), field_a(nx,ny,nz)
      real(kind=8)     field_b(nx,ny,nz)
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               field_c(i,j,k) = field_a(i,j,k) + field_b(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: mat_mat_add3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mat_copy1d:                                                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mat_copy1d( field_n, field_np1, nx)
      implicit   none
      integer    nx
      real(kind=8)     field_n(nx), field_np1(nx)
      integer    i

            do i = 1, nx
               field_np1(i) = field_n(i)
            end do

      return
      end    ! END: mat_copy1d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mat_copy3d:                                                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mat_copy3d( field_n, field_np1, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field_n(nx,ny,nz), field_np1(nx,ny,nz)
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               field_np1(i,j,k) = field_n(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: mat_copy3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  subtract_field3dL2:                                                       cc
cc                    Compute L2 norm of diff.                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function subtract_field3dL2(pfield, mfield, nx,
     .                                                           ny,nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     pfield(nx,ny,nz)
      real(kind=8)     mfield(nx,ny,nz)
      real(kind=8)     tmp
      integer    i,j,k

      tmp = 0.d0
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               tmp = tmp + (pfield(i,j,k) - mfield(i,j,k))**2
            end do
         end do
      end do

      subtract_field3dL2 = nx*ny*nz
      if (subtract_field3dL2.gt.0) then
         subtract_field3dL2 = sqrt(tmp/subtract_field3dL2)
      else
         subtract_field3dL2 = 0.d0
      end if

      return
      end    ! END: subtract_field3dL2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  subtract_field3d:                                                         cc
cc                    subtract two 3d vectors.                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine subtract_field3d( fieldout, pfield, mfield, nx,ny,nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     fieldout(nx,ny,nz), pfield(nx,ny,nz)
      real(kind=8)     mfield(nx,ny,nz)
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               fieldout(i,j,k) = pfield(i,j,k) - mfield(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: subtract_field3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  vec_scalmult1D:                                                           cc
cc                  Multiple existing vector by a scalar.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vec_scalmult1D(field, scalar, nx)
      implicit   none
      integer    nx
      real(kind=8)     field(nx), scalar
      integer    i

            do i = 1, nx
               field(i) = scalar * field(i)
            end do

      return
      end    ! END: vec_scalmult1D

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scal_mult3d:                                                         cc
cc                    Load scalar*fieldin to fieldout.                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scal_mult3d( fieldout, fieldin, scalar, nx,ny,nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     fieldout(nx,ny,nz), fieldin(nx,ny,nz), scalar
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               fieldout(i,j,k) = scalar * fieldin(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: load_scal_mult3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scal_add1d:                                                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scal_add1d( field, scalar, nx)
      implicit   none
      integer    nx
      real(kind=8)     field(nx), scalar
      integer    i

            do i = 1, nx
               field(i) = field(i) + scalar
            end do

      return
      end    ! END: load_scal_add1d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scal1d:                                                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scal1d( field, scalar, nx)
      implicit   none
      integer    nx
      real(kind=8)     field(nx), scalar
      integer    i

            do i = 1, nx
               field(i) = scalar
            end do

      return
      end    ! END: load_scal1d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scal_bound3d:                                                        cc
cc                    Loads scalar into boundary                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scal_bound3d( field, scalar, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), scalar
      integer    i,j,k

      !
      !   X boundaries
      !
      do i = 1, nx, nx-1
         do k = 1, nz
            do j = 1, ny
               field(i,j,k) = scalar
            end do
         end do
      end do
      !
      !   Y boundaries
      !
      do j = 1, ny, ny-1
         do k = 1, nz
            do i = 1, nx
               field(i,j,k) = scalar
            end do
         end do
      end do
      !
      !   Z boundaries
      !
      do k = 1, nz, nz-1
         do j = 1, ny
            do i = 1, nx
               field(i,j,k) = scalar
            end do
         end do
      end do

      return
      end    ! END: load_scal_bound3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scal_bound3dW:                                                       cc
cc                    Loads scalar into boundary of some width                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scal_bound3dW( field, scalar, width, sym,nx,ny,nz)
      implicit   none
      integer    nx, ny, nz, width, sym
      real(kind=8)     field(nx,ny,nz), scalar
      integer    i,j,k

      !
      !   X boundaries
      !
         if (sym .ne.1.and.sym.ne.6) then
         do k = 1, nz
            do j = 1, ny
               do i = 1, width
                  field(i,j,k) = scalar
               end do
            end do
         end do
         end if
         do k = 1, nz
            do j = 1, ny
               do i = nx-width+1, nx
                  field(i,j,k) = scalar
               end do
            end do
         end do
      !
      !   Y boundaries
      !
         if (sym .ne.2.and.sym.ne.6) then
         do k = 1, nz
            do j = 1, width
               do i = 1, nx
                  field(i,j,k) = scalar
               end do
            end do
         end do
         end if
         do k = 1, nz
            do j = ny-width+1, ny
               do i = 1, nx
                  field(i,j,k) = scalar
               end do
            end do
         end do
      !
      !   Z boundaries
      !
         if (sym .ne.3.and.sym.ne.6) then
            do k = 1, width
               do j = 1, ny
                  do i = 1, nx
                     field(i,j,k) = scalar
                  end do
               end do
            end do
         end if
         do k = nz-width+1, nz
            do j = 1, ny
               do i = 1, nx
                  field(i,j,k) = scalar
               end do
            end do
         end do

      return
      end    ! END: load_scal_bound3dW

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scal_add3d:                                                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scal_add3d( field, scalar, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), scalar
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               field(i,j,k) = field(i,j,k) + scalar
            end do
         end do
      end do

      return
      end    ! END: load_scal_add3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scal3d:                                                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scal3d( field, scalar, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), scalar
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               field(i,j,k) = scalar
            end do
         end do
      end do

      return
      end    ! END: load_scal3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scalallz:                                                            cc
cc                 Specialized routine to flag points for all of z, but       cc
cc                 only a fraction of x and y.                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scalallz( field, level, h,mx,my,mz,
     *                            scalar, nx,ny,nz,nxc,nyc,nzc)
      implicit   none
      integer    nx, ny, nz, level
      ! Coarsest grid values:
      integer    nxc,nyc,nzc
      ! grid spacing and minimums for when symmetry is assumed:
      real(kind=8)     h, mx,my,mz
      real(kind=8)     field(nx,ny,nz), scalar
      integer    i,j,k
      integer    centeri, centerj, centerk
      integer    widthi,  widthj,  widthk
      integer    bi, ei, bj, ej, bk, ek
      integer    izero, jzero, kzero
      real(kind=8)     cr
      logical          includeregionsatZextremes
      parameter      ( includeregionsatZextremes = .false. )
      logical          ltrace
      parameter      ( ltrace = .false. )

      !cr   = linfract
      !
      ! Hardcoded to refine only half the x- and y- directions:
      !
      cr   = 2.d0

      centeri = NINT( 0.5d0 * (nx+1) )
      centerj = NINT( 0.5d0 * (ny+1) )
      centerk = NINT( 0.5d0 * (nz+1) )

      if (.not.includeregionsatZextremes) then
         !This doesn't work when includeregionsatZextremes==true:
         widthi  = centeri / cr
         widthj  = centerj / cr
         widthk  = centerk / cr
      else
         !
         ! Arbitrary fixed choice:
         !widthi  = 10
         !widthj  = 10
         !widthk  = 10
         !
         ! Arbitrary fraction of grid (in terms of coarsest)
         widthi  = (nxc+1)/4
         widthj  = (nyc+1)/4
         widthk  = (nzc+1)/4
      end if

      bi      = centeri - widthi
      bj      = centerj - widthj
      bk      = centerk - widthk

      ei      = centeri + widthi
      ej      = centerj + widthj
      ek      = centerk + widthk

      if (bi.lt.1) bi=1
      if (bj.lt.1) bj=1
      if (bk.lt.1) bk=1
      !
      if (ei.gt.nx) ei=nx
      if (ej.gt.ny) ej=ny
      if (ek.gt.nz) ek=nz

      !
      ! Refine entire z-axis:
      !
      bk = 1
      ek = nz

      if (ltrace) then
        write(*,*)'load_scalallz:         cr:', cr
        write(*,*)'load_scalallz:     scalar:', scalar
        write(*,*)'load_scalallz:centeri/j/k:',centeri,centerj,centerk
        write(*,*)'load_scalallz: widthi/j/k: ',widthi,widthj,widthk
        write(*,*)'load_scalallz:      kzero:', kzero
        write(*,*)'load_scalallz:      level:', level
        write(*,*)'load_scalallz:          h:', h 
        write(*,*)'load_scalallz:    m/x/y/z:', mx, my, mz
        write(*,*)'load_scalallz:          i:', bi, ei, nx
        write(*,*)'load_scalallz:          j:', bj, ej, ny
        write(*,*)'load_scalallz:          k:', bk, ek, nz
      end if


      do k = bk, ek
      do j = bj, ej
      do i = bi, ei
         field(i,j,k) = scalar
      end do
      end do
      end do

      if (includeregionsatZextremes) then
         !widthk = NINT(0.2d0 * nz)  --> does not leave buffer regions
         !                               between multiple refined regions
         ! The following value is arbitrary...choose something 
         ! sufficient to enforce whatever condition you want:
         !widthk = 8
         widthk = (nzc+1)/5
         if (ltrace) then
           write(*,*)'load_scalallz: Including regions at extremes of z'
           write(*,*)'load_scalallz: of width (in points): ',widthk
         end if
         kzero = 0
         do k = 1, widthk
         do j = 1, ny
         do i = 1, nx
            kzero = kzero + 1
            field(i,j,k) = scalar
            !write(*,*)'load_scalallz: (i,j,k):',i,j,k,field(i,j,k),kzero
         end do
         end do
         end do
         do k = nz-widthk+1, nz
         do j = 1, ny
         do i = 1, nx
            kzero = kzero + 1
            field(i,j,k) = scalar
         end do
         end do
         end do
         if(ltrace) write(*,*) 'load_scalallz: kzero = ',kzero
      end if

      if (ltrace) call field_dump_stats(field, nx, ny, nz, 'flagk')
   
      return
      end    ! END: load_scalallz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scalallz2:                                                           cc
cc                 Specialized routine to flag points for                     cc
cc                 only half the extents of x and y and an equivalent amount  cc
cc                 in z ....to maintain cubical refined regions in a coarse   cc
cc                 level that is not cubical.                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scalallz2( field, linfract, h,mx,my,mz,
     *                            scalar, nx,ny,nz)
      implicit   none
      integer    nx, ny, nz
      ! grid spacing and minimums for when symmetry is assumed:
      real(kind=8)     h, mx,my,mz
      real(kind=8)     field(nx,ny,nz), scalar, linfract
      integer    i,j,k
      integer    centeri, centerj, centerk
      integer    widthi,  widthj,  widthk
      integer    bi, ei, bj, ej, bk, ek
      integer    izero, jzero, kzero
      real(kind=8)     cr
      logical          ltrace
      parameter      ( ltrace = .false. )

      !cr   = linfract
      !
      ! Hardcoded to refine only half the x- and y- directions:
      !
      cr   = 2.d0

      centeri = NINT( 0.5d0 * (nx+1) )
      centerj = NINT( 0.5d0 * (ny+1) )
      centerk = NINT( 0.5d0 * (nz+1) )

      widthi  = centeri / cr
      widthj  = centerj / cr
      widthk  = centerk / cr
      widthk  = max(widthi,widthj)

      bi      = centeri - widthi
      bj      = centerj - widthj
      bk      = centerk - widthk

      ei      = centeri + widthi
      ej      = centerj + widthj
      ek      = centerk + widthk

      if (bi.lt.1) bi=1
      if (bj.lt.1) bj=1
      if (bk.lt.1) bk=1
      !
      if (ei.gt.nx) ei=nx
      if (ej.gt.ny) ej=ny
      if (ek.gt.nz) ek=nz

      !!
      !! Refine entire z-axis:
      !!
      !bk = 1
      !ek = nz

      if (ltrace) then
        write(*,*)'load_scalallz2:         cr:', cr
        write(*,*)'load_scalallz2:     scalar:', scalar
        write(*,*)'load_scalallz2:centeri/j/k:',centeri,centerj,centerk
        write(*,*)'load_scalallz2: widthi/j/k: ',widthi,widthj,widthk
        write(*,*)'load_scalallz2:      kzero:', kzero
        write(*,*)'load_scalallz2:          h:', h 
        write(*,*)'load_scalallz2:    m/x/y/z:', mx, my, mz
        write(*,*)'load_scalallz2:          i:', bi, ei, nx
        write(*,*)'load_scalallz2:          j:', bj, ej, ny
        write(*,*)'load_scalallz2:          k:', bk, ek, nz
      end if


      do k = bk, ek
      do j = bj, ej
      do i = bi, ei
         field(i,j,k) = scalar
      end do
      end do
      end do

      return
      end    ! END: load_scalallz2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_scalfrac3d:                                                          cc
cc                    Load scalar in fraction volume of fractional linear     cc
cc                    dimension of linfract.                                  cc
cc                NB: sym--> symmetry conditions (if 0, no such conditions)   cc
cc                                               (if 1, reflection about x=0) cc
cc                                               (if 2, reflection about y=0) cc
cc                                               (if 3, reflection about z=0) cc
cc                                               (if 6, reflection about all) cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_scalfrac3d( field, linfract, h,mx,my,mz,
     *                            scalar, sym,extension, nx,ny,nz)
      implicit   none
      integer    nx, ny, nz, sym, extension
      ! grid spacing and minimums for when symmetry is assumed:
      real(kind=8)     h, mx,my,mz
      real(kind=8)     field(nx,ny,nz), scalar, linfract
      integer    i,j,k
      integer    centeri, centerj, centerk
      integer    widthi,  widthj,  widthk
      integer    bi, ei, bj, ej, bk, ek
      integer    izero, jzero, kzero
      real(kind=8)     cr
      logical          ltrace
      parameter      ( ltrace = .false. )

      cr   = linfract

      centeri = NINT( 0.5d0 * (nx+1) )
      centerj = NINT( 0.5d0 * (ny+1) )
      centerk = NINT( 0.5d0 * (nz+1) )

      widthi  = centeri / cr
      widthj  = centerj / cr
      widthk  = centerk / cr

      bi      = centeri - widthi
      bj      = centerj - widthj
      bk      = centerk - widthk

      ei      = centeri + widthi
      ej      = centerj + widthj
      ek      = centerk + widthk

      if (bi.lt.1) bi=1
      if (bj.lt.1) bj=1
      if (bk.lt.1) bk=1
      !
      if (ei.gt.nx) ei=nx
      if (ej.gt.ny) ej=ny
      if (ek.gt.nz) ek=nz

      if (mz.le.0.and.(sym.eq.3.or.sym.eq.6)) then
         kzero   = -NINT(mz/h) + 1
         centerk = kzero
         widthk  = (nz-kzero) / cr
         bk      = kzero - extension
         if (bk.lt.1) bk=1
         ek      = centerk + widthk
      end if
      if (my.le.0.and.(sym.eq.2.or.sym.eq.6)) then
         jzero   = -NINT(my/h) + 1
         centerj = jzero
         widthj  = (ny-jzero) / cr
         bj      = jzero - extension
         if (bj.lt.1) bj=1
         ej      = centerj + widthj
      end if
      if (mz.le.0.and.(sym.eq.1.or.sym.eq.6)) then
         izero   = -NINT(mx/h) + 1
         centeri = izero
         widthi  = (nz-izero) / cr
         bi      = izero - extension
         if (bi.lt.1) bi=1
         ei      = centeri + widthi
      end if

      if (ltrace) then
        write(*,*)'load_scalfrac3d:         cr:', cr
        write(*,*)'load_scalfrac3d:     scalar:', scalar
        write(*,*)'load_scalfrac3d:centeri/j/k:',centeri,centerj,centerk
        write(*,*)'load_scalfrac3d: widthi/j/k: ',widthi,widthj,widthk
        write(*,*)'load_scalfrac3d:  extension:', extension
        write(*,*)'load_scalfrac3d:      kzero:', kzero
        write(*,*)'load_scalfrac3d:          h:', h 
        write(*,*)'load_scalfrac3d:    m/x/y/z:', mx, my, mz
        write(*,*)'load_scalfrac3d:          i:', bi, ei, nx
        write(*,*)'load_scalfrac3d:          j:', bj, ej, ny
        write(*,*)'load_scalfrac3d:          k:', bk, ek, nz
      end if


      do k = bk, ek
      do j = bj, ej
      do i = bi, ei
         field(i,j,k) = scalar
      end do
      end do
      end do

      return
      end    ! END: load_scalfrac3d


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mult_field_r:                                                             cc
cc                    Multiply a field by radial distance "r".                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mult_field_r( field, x, y, z, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), x(nx), y(ny), z(nz)
      real(kind=8)     r
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               r = sqrt( x(i)**2 + y(j)**2 + z(k)**2 )
               field(i,j,k) = field(i,j,k) * r
            end do
         end do
      end do

      return
      end    ! END: mult_field_r


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mult_field_x:                                                             cc
cc                    Multiply a field by coordinate "x".                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mult_field_x( field, x, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), x(nx)
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               field(i,j,k) = field(i,j,k) * x(i)
            end do
         end do
      end do

      return
      end    ! END: mult_field_x


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mult_field_y:                                                             cc
cc                    Multiply a field by coordinate "y".                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mult_field_y( field, y, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), y(ny)
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               field(i,j,k) = field(i,j,k) * y(j)
            end do
         end do
      end do

      return
      end    ! END: mult_field_y

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mult_field_z:                                                             cc
cc                    Multiply a field by coordinate "z".                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mult_field_z( field, z, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     field(nx,ny,nz), z(nz)
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               field(i,j,k) = field(i,j,k) * z(k)
            end do
         end do
      end do

      return
      end    ! END: mult_field_z

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  vec_vec_mult:                                                             cc
cc                Multiply two fields together:                               cc
cc                A = B * C                                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vec_vec_mult( A, B, C, nx,ny,nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     A(nx,ny,nz), B(nx,ny,nz), C(nx,ny,nz)
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               A(i,j,k) = B(i,j,k) * C(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: vec_vec_mult

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  field_exp3d:                                                              cc
cc                    Compute and store fieldout = exp(fieldin)               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine field_exp3d( fieldout, fieldin, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     fieldout(nx,ny,nz), fieldin(nx,ny,nz)
      integer    i,j,k

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               fieldout(i,j,k) = exp(fieldin(i,j,k))
            end do
         end do
      end do

      return
      end    ! END: field_exp3d

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  deriv2_y:                                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deriv2_y( fieldy, field, h, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     fieldy(nx,ny,nz), field(nx,ny,nz), h
      integer    i,j,k

      do i = 1, nx
         do j = 2, ny-1
            do k = 1, nz
               fieldy(i,j,k) = 0.5d0*(field(i,j+1,k) - field(i,j-1,k))/h
c              if(fieldy(i,j,k).gt.0.00001)write(*,*)i,j,k,fieldy(i,j,k)
            end do
         end do
      end do

      j = 1
         do i = 1, nx
            do k = 1, nz
               fieldy(i,j,k) = ( -3.d0*field(i,j,k)
     *                           +4.d0*field(i,j+1,k)
     *                           -     field(i,j+2,k)
     *                           ) * 0.5d0 / h
            end do
         end do

      j = ny
         do i = 1, nx
            do k = 1, nz
               fieldy(i,j,k) = (  3.d0*field(i,j,k)
     *                           -4.d0*field(i,j-1,k)
     *                           +     field(i,j-2,k)
     *                           ) * 0.5d0 / h
            end do
         end do

      return
      end    ! END: deriv2_y

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  deriv2_z:                                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deriv2_z( fieldz, field, h, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     fieldz(nx,ny,nz), field(nx,ny,nz), h
      integer    i,j,k

      do i = 1, nx
         do j = 1, ny
            do k = 2, nz-1
               fieldz(i,j,k) = 0.5d0*(field(i,j,k+1) - field(i,j,k-1))/h
            end do
         end do
      end do

      k = 1
         do i = 1, nx
            do j = 1, ny
               fieldz(i,j,k) = ( -3.d0*field(i,j,k)
     *                           +4.d0*field(i,j,k+1)
     *                           -     field(i,j,k+2)
     *                           ) * 0.5d0 / h
            end do
         end do

      k = nz
         do i = 1, nx
            do j = 1, ny
               fieldz(i,j,k) = (  3.d0*field(i,j,k)
     *                           -4.d0*field(i,j,k-1)
     *                           +     field(i,j,k-2)
     *                           ) * 0.5d0 / h
            end do
         end do

      return
      end    ! END: deriv2_z

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  deriv2_x:                                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deriv2_x( fieldx, field, h, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     fieldx(nx,ny,nz), field(nx,ny,nz), h
      integer    i,j,k

      do i = 2, nx-1
         do j = 1, ny
            do k = 1, nz
               fieldx(i,j,k) = 0.5d0*(field(i+1,j,k) - field(i-1,j,k))/h
            end do
         end do
      end do

      i = 1
         do j = 1, ny
            do k = 1, nz
               fieldx(i,j,k) = ( -3.d0*field(i,j,k)
     *                           +4.d0*field(i+1,j,k)
     *                           -     field(i+2,j,k)
     *                           ) * 0.5d0 / h
            end do
         end do

      i = nx
         do j = 1, ny
            do k = 1, nz
               fieldx(i,j,k) = (  3.d0*field(i,j,k)
     *                           -4.d0*field(i-1,j,k)
     *                           +     field(i-2,j,k)
     *                           ) * 0.5d0 / h
            end do
         end do

      return
      end    ! END: deriv2_x

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  gaussian3d_add:                                                           cc
cc                   Adds Gaussian pulse to field.                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gaussian3d_add( field, h, minx, miny, minz,nx,ny,nz,
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
         write(*,*) 'gaussian3d_add:   h = ',h
         write(*,*) 'gaussian3d_add: minx = ',minx
         write(*,*) 'gaussian3d_add: miny = ',miny
         write(*,*) 'gaussian3d_add: minz = ',minz
         write(*,*) 'gaussian3d_add:    a = ',a
         write(*,*) 'gaussian3d_add:   ex = ',ex
         write(*,*) 'gaussian3d_add:   ey = ',ey
         write(*,*) 'gaussian3d_add:    r = ',r
         write(*,*) 'gaussian3d_add:    d = ',d
         write(*,*) 'gaussian3d_add:   xc = ',xc
         write(*,*) 'gaussian3d_add:   yc = ',yc
         write(*,*) 'gaussian3d_add:   zc = ',zc
      end if

      do k = 1, nz
         z = minz + h*(k-1.d0)
         do j = 1, ny
            y = miny + h*(j-1.d0)
            do i = 1, nx
               x = minx + h*(i-1.d0)
               field(i,j,k) = field(i,j,k)
     *                      + a*exp( -(sqrt(
     *             ex*(x-xc)**2 + ey*(y-yc)**2 + (z-zc)**2)-r)**2/d**2 )
c    *             ex*(x-xc)**2 + ey*(y-yc)**2            )-r)**2/d**2 )
            end do
         end do
      end do

      if (ltrace) then
         write(*,*) 'gaussian3d_add: maxx = ',x
         write(*,*) 'gaussian3d_add: maxy = ',y
         write(*,*) 'gaussian3d_add: maxz = ',z
      end if

      return
      end    ! END: gaussian3d_add

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  paste_boundaries:                                                         cc
cc                                                                            cc
cc      copy data from work space into boundary data of a field               cc
cc      assuming that the necessary interpolation has happened everywhere.    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine paste_boundaries(  values,field,     chr,
     *                              mini,  minj,      mink,
     *                              maxi,  maxj,      maxk,
     *                              nx,    ny,        nz,
     *                              r, gw)
      implicit    none
      integer     nx,ny,nz, r, gw
      integer     mini, minj, mink,  maxi, maxj, maxk
      integer     minip,minjp,minkp, maxip,maxjp,maxkp
      real(kind=8)      values(nx,ny,nz),
     *            field(nx,ny,nz), chr(nx,ny,nz)
      include    'chr.inc'
      integer     i,j,k, width
      integer     ii,jj,kk, tmpi
      integer     imid,jmid,kmid,  AMRBOUNDARY,DECOBOUNDARY,REFLBOUNDARY
      logical     nearamr
      integer     myid, proc_return_myid
      !
      logical     ltrace
      parameter ( ltrace   = .false. )
      logical     ltrace2
      parameter ( ltrace2  = .false. )

      !myid         = proc_return_myid()
      !ltrace       = .false.
      !ltrace2      = .false.
      !if (myid.eq.3) then
      !ltrace       = .true.
      !ltrace2      = .true.
      !end if

      AMRBOUNDARY  = NINT(CHR_amr_bdy)
      DECOBOUNDARY = NINT(CHR_deco_bdy)
      REFLBOUNDARY = NINT(CHR_refl_bdy)

      if (ltrace) then
         write(*,*) 'paste_boundaries: nx/y/z:   ',nx, ny, nz
         write(*,*) 'paste_boundaries: mini/j/k: ',mini,minj,mink
         write(*,*) 'paste_boundaries: maxi/j/k: ',maxi,maxj,maxk
         write(*,*) 'paste_boundaries: gw: ',gw
      end if

      imid = (mini+maxi)/2
      jmid = (minj+maxj)/2
      kmid = (mink+maxk)/2

      !
      ! Low Z Boundary:
      !
      if (mink .ne. 1) goto 10
      width = gw
      if (width .gt. maxk) width = maxk
      if (NINT(chr(imid,jmid,1)) .eq. REFLBOUNDARY ) goto 10
      if (ltrace) write(*,*) 'paste_boundaries: Low  Z k: ',1,width
      do j = minj, maxj
         do i = mini, maxi
            if (NINT(chr(i,j,1)) .eq. AMRBOUNDARY) then
               do k = 1, width
                  field(i,j,k) = values(i,j,k)
               end do
            else
               !
               ! for non-triviallly intersecting grids, at the
               ! inside corners, points need to be set by 
               ! interpolation which are near the corner,
               ! so check points within ghostwidth of this one
               !   So if this bdy point is not an AMR bdy pt,
               ! but is within gw of an AMR bdy pt, then we
               ! should still set the points.
               if(ltrace2)write(*,*)'pasteLZ:',tmpi,i,j
               tmpi = 1
               nearamr = .false.
 23            continue
               if (i.gt. 1+tmpi .and. j.gt.1+tmpi .and.
     .             i.lt. nx-tmpi.and. j.lt.ny-tmpi      ) then
                  if (     NINT(chr(i+tmpi,j,     1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. j.gt.1+tmpi .and.
     .                  i.lt. nx-tmpi                        ) then
                  if (     NINT(chr(i+tmpi,j,     1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. j.gt.1+tmpi .and.
     .                                     j.lt.ny-tmpi      ) then
                  if (     NINT(chr(i,     j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and.
     .                  i.lt. nx-tmpi.and. j.lt.ny-tmpi      ) then
                  if (     NINT(chr(i+tmpi,j,     1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (                   j.gt.1+tmpi .and.
     .                  i.lt. nx-tmpi.and. j.lt.ny-tmpi      ) then
                  if (     NINT(chr(i+tmpi,j,     1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. j.gt.1+tmpi ) then
                  if (     NINT(chr(i-tmpi,j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. j.lt.ny-tmpi ) then
                  if (     NINT(chr(i,     j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if ( i.lt. nx-tmpi.and. j.lt.ny-tmpi) then
                  if (     NINT(chr(i+tmpi,j,     1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j+tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j+tmpi,1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.lt. nx-tmpi.and. j.gt.1+tmpi) then
                  if (     NINT(chr(i+tmpi,j,     1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j-tmpi,1)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1 +tmpi) then
                  if (     NINT(chr(i-tmpi,j,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.lt. nx-tmpi) then
                  if (     NINT(chr(i+tmpi,j,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1 +tmpi) then
                  if (     NINT(chr(i,j-tmpi,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.lt. ny-tmpi) then
                  if (     NINT(chr(i,j+tmpi,     1)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               end if
               !
               if (nearamr) then
                  ! tmpi is the distance to an AMR point:
                  ! fill in values according to distance:
                  do k = 1, min(gw-tmpi,width)
                  !do k = 1, gw-tmpi
                     field(i,j,k) = values(i,j,k)
                  end do
               else
                  if (tmpi .lt. gw-1) then
                     tmpi = tmpi + 1
                     goto 23
                  end if
               end if
               !
            end if
         end do
      end do

 10   continue

      !
      ! High Z Boundary:
      !
      if (maxk .ne. nz) goto 20
      width = nz-gw+1
      if (width .lt. mink) width = mink
      if (ltrace) write(*,*) 'paste_boundaries: High Z k: ',nz, width
      do j = minj, maxj
         do i = mini, maxi
            if (NINT(chr(i,j,nz)) .eq. AMRBOUNDARY) then
               do k = nz, width, -1
                  field(i,j,k) = values(i,j,k)
               end do
            else
               !
               ! for non-triviallly intersecting grids, at the
               ! inside corners, points need to be set by 
               ! interpolation which are near the corner,
               ! so check points within ghostwidth of this one
               !   So if this bdy point is not an AMR bdy pt,
               ! but is within gw of an AMR bdy pt, then we
               ! should still set the points.
               if(ltrace2)write(*,*)'pasteHZ:',tmpi,i,j
               tmpi = 1
               nearamr = .false.
 33            continue
               if (i.gt. 1+tmpi .and. j.gt.1+tmpi .and.
     .             i.lt. nx-tmpi.and. j.lt.ny-tmpi      ) then
                  if (     NINT(chr(i+tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. j.gt.1+tmpi .and.
     .                  i.lt. nx-tmpi                        ) then
                  if (     NINT(chr(i+tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. j.gt.1+tmpi .and.
     .                                     j.lt.ny-tmpi      ) then
                  if (     NINT(chr(i,     j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and.
     .                  i.lt. nx-tmpi.and. j.lt.ny-tmpi      ) then
                  if (     NINT(chr(i+tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (                   j.gt.1+tmpi .and.
     .                  i.lt. nx-tmpi.and. j.lt.ny-tmpi      ) then
                  if (     NINT(chr(i+tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. j.gt.1+tmpi ) then
                  if (     NINT(chr(i-tmpi,j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. j.lt.ny-tmpi ) then
                  if (     NINT(chr(i,     j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if ( i.lt. nx-tmpi.and. j.lt.ny-tmpi) then
                  if (     NINT(chr(i+tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j+tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j+tmpi,nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.lt. nx-tmpi.and. j.gt.1+tmpi) then
                  if (     NINT(chr(i+tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,j-tmpi,nz)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     j-tmpi,nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1 +tmpi) then
                  if (     NINT(chr(i-tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.lt. nx-tmpi) then
                  if (     NINT(chr(i+tmpi,j,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1 +tmpi) then
                  if (     NINT(chr(i,j-tmpi,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.lt. ny-tmpi) then
                  if (     NINT(chr(i,j+tmpi,     nz)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               end if
               !
               if (nearamr) then
                  ! tmpi is the distance to an AMR point:
                  ! fill in values according to distance:
                  do k = nz, max(nz-gw+tmpi,width), -1
                  !do k = nz, nz-gw+tmpi, -1
                     field(i,j,k) = values(i,j,k)
                  end do
               else
                  if (tmpi .lt. gw-1) then
                     tmpi = tmpi + 1
                     goto 33
                  end if
               end if
               !
            end if
         end do
      end do

 20   continue

      !
      ! Low Y Boundary:
      !
      if (minj .ne. 1) goto 30
      width = gw
      if (width .gt. maxj) width = maxj
      if (NINT(chr(imid,1,kmid)) .eq. REFLBOUNDARY) goto 30
      if (ltrace) write(*,*) 'paste_boundaries: Low  Y j: ',1,width
      do k = mink, maxk
         do i = mini, maxi
            if (NINT(chr(i,1,k)) .eq. AMRBOUNDARY) then
               do j = 1, width
                  field(i,j,k) = values(i,j,k)
               end do
            else
               !
               ! for non-triviallly intersecting grids, at the
               ! inside corners, points need to be set by 
               ! interpolation which are near the corner,
               ! so check points within ghostwidth of this one
               !   So if this bdy point is not an AMR bdy pt,
               ! but is within gw of an AMR bdy pt, then we
               ! should still set the points.
               if(ltrace2)write(*,*)'pasteLY:',tmpi,i,k
               tmpi = 1
               nearamr = .false.
 43            continue
               if (i.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .             i.lt. nx-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(i+tmpi,1,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .                  i.lt. nx-tmpi                        ) then
                  if (     NINT(chr(i+tmpi,1,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .                                     k.lt.nz-tmpi      ) then
                  if (     NINT(chr(i,     1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and.
     .                  i.lt. nx-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(i+tmpi,1,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (                   k.gt.1+tmpi .and.
     .                  i.lt. nx-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(i+tmpi,1,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. k.gt.1+tmpi ) then
                  if (     NINT(chr(i-tmpi,1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. k.lt.nz-tmpi ) then
                  if (     NINT(chr(i,     1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,1,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if ( i.lt. nx-tmpi.and. k.lt.nz-tmpi) then
                  if (     NINT(chr(i+tmpi,1,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,1,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k+tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.lt. nx-tmpi.and. k.gt.1+tmpi) then
                  if (     NINT(chr(i+tmpi,1,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,1,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     1,k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1 +tmpi) then
                  if (     NINT(chr(i-tmpi,1,k)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.lt. nx-tmpi) then
                  if (     NINT(chr(i+tmpi,1,k)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (k.gt. 1 +tmpi) then
                  if (     NINT(chr(i,1,k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (k.lt. ny-tmpi) then
                  if (     NINT(chr(i,1,k+tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               end if
               !
               !
               if (nearamr) then
                  ! tmpi is the distance to an AMR point:
                  ! fill in values according to distance:
                  do j = 1, min(gw-tmpi,width)
                  !do j = 1, gw-tmpi
                     field(i,j,k) = values(i,j,k)
                  end do
               else
                  if (tmpi .lt. gw-1) then
                     tmpi = tmpi + 1
                     goto 43
                  end if
               end if
               !
            end if
         end do
      end do

 30   continue

      !
      ! High Y Boundary:
      !
      if (maxj .ne. ny) goto 40
      width = ny-gw+1
      if (width .lt. minj) width = minj
      if (ltrace) write(*,*) 'paste_boundaries: High Y j: ',ny,width
      do k = mink, maxk
         do i = mini, maxi
            if (NINT(chr(i,ny,k)) .eq. AMRBOUNDARY) then
               do j = ny, width, -1
                  field(i,j,k) = values(i,j,k)
               end do
            else
               !
               ! for non-triviallly intersecting grids, at the
               ! inside corners, points need to be set by 
               ! interpolation which are near the corner,
               ! so check points within ghostwidth of this one
               !   So if this bdy point is not an AMR bdy pt,
               ! but is within gw of an AMR bdy pt, then we
               ! should still set the points.
               if(ltrace2)write(*,*)'pasteHY:',tmpi,i,k
               tmpi = 1
               nearamr = .false.
 53            continue
               if (i.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .             i.lt. nx-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(i+tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .                  i.lt. nx-tmpi                        ) then
                  if (     NINT(chr(i+tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .                                     k.lt.nz-tmpi      ) then
                  if (     NINT(chr(i,     ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and.
     .                  i.lt. nx-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(i+tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (                   k.gt.1+tmpi .and.
     .                  i.lt. nx-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(i+tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. k.gt.1+tmpi ) then
                  if (     NINT(chr(i-tmpi,ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1+tmpi .and. k.lt.nz-tmpi ) then
                  if (     NINT(chr(i,     ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i-tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if ( i.lt. nx-tmpi.and. k.lt.nz-tmpi) then
                  if (     NINT(chr(i+tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,ny,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k+tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.lt. nx-tmpi.and. k.gt.1+tmpi) then
                  if (     NINT(chr(i+tmpi,ny,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i+tmpi,ny,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(i,     ny,k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.gt. 1 +tmpi) then
                  if (     NINT(chr(i-tmpi,ny,k)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (i.lt. nx-tmpi) then
                  if (     NINT(chr(i+tmpi,ny,k)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (k.gt. 1 +tmpi) then
                  if (     NINT(chr(i,ny,k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (k.lt. ny-tmpi) then
                  if (     NINT(chr(i,ny,k+tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               end if
               !
               if (nearamr) then
                  ! tmpi is the distance to an AMR point:
                  ! fill in values according to distance:
                  do j = ny, max(ny-gw+tmpi,width), -1
                  !do j = ny, ny-gw+tmpi, -1
                     field(i,j,k) = values(i,j,k)
                  end do
               else
                  if (tmpi .lt. gw-1) then
                  !if (tmpi .lt. gw) then
                     tmpi = tmpi + 1
                     goto 53
                  end if
               end if
               !
            end if
         end do
      end do

 40   continue

      !
      ! Low X Boundary:
      !
      if (mini .ne. 1) goto 50
      width = gw
      if (width .gt. maxi) width = maxi
      if (NINT(chr(1,jmid,kmid)) .eq. REFLBOUNDARY) goto 50
      if (ltrace) write(*,*) 'paste_boundaries: Low  X i: ',1,width
      do k = mink, maxk
         do j = minj, maxj
            if (NINT(chr(1,j,k)) .eq. AMRBOUNDARY) then
               do i = 1, width
                  field(i,j,k) = values(i,j,k)
               end do
            else
               !
               ! for non-triviallly intersecting grids, at the
               ! inside corners, points need to be set by 
               ! interpolation which are near the corner,
               ! so check points within ghostwidth of this one
               !   So if this bdy point is not an AMR bdy pt,
               ! but is within gw of an AMR bdy pt, then we
               ! should still set the points.
               if(ltrace2)write(*,*)'pasteLX:',tmpi,j,k
               tmpi = 1
               nearamr = .false.
 63            continue
               if (j.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .             j.lt. ny-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(1,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j+tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j+tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .                  j.lt. ny-tmpi                        ) then
                  if (     NINT(chr(1,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j+tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .                                     k.lt.nz-tmpi      ) then
                  if (     NINT(chr(1,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and.
     .                  j.lt. ny-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(1,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j+tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (                   k.gt.1+tmpi .and.
     .                  j.lt. ny-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(1,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j+tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j+tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and. k.gt.1+tmpi ) then
                  if (     NINT(chr(1,j-tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and. k.lt.nz-tmpi ) then
                  if (     NINT(chr(1,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if ( j.lt. ny-tmpi.and. k.lt.nz-tmpi) then
                  if (     NINT(chr(1,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j+tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.lt. ny-tmpi.and. k.gt.1+tmpi) then
                  if (     NINT(chr(1,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j+tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(1,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1 +tmpi) then
                  if (     NINT(chr(1,j-tmpi,k)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.lt. ny-tmpi) then
                  if (     NINT(chr(1,j+tmpi,k)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (k.gt. 1 +tmpi) then
                  if (     NINT(chr(1,j,k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (k.lt. nz-tmpi) then
                  if (     NINT(chr(1,j,k+tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               end if
               !
               !
               if (nearamr) then
                  ! tmpi is the distance to an AMR point:
                  ! fill in values according to distance:
                  do i = 1, min(gw-tmpi,width)
                  !do i = 1, gw-tmpi
                     field(i,j,k) = values(i,j,k)
                  end do
               else
                  if (tmpi .lt. gw-1) then
                  !if (tmpi .lt. gw) then
                     tmpi = tmpi + 1
                     goto 63
                  end if
               end if
               !
            end if
         end do
      end do

 50   continue

      !
      ! High X Boundary:
      !
      if (maxi .ne. nx) goto 60
      width = nx-gw+1
      if (width .lt. mini) width = mini
      if (ltrace) write(*,*) 'paste_boundaries: High X i: ',nx,width
      do k = mink, maxk
         do j = minj, maxj
            if (NINT(chr(nx,j,k)) .eq. AMRBOUNDARY) then
               do i = nx, width, -1
                  field(i,j,k) = values(i,j,k)
               end do
            else
               !
               ! for non-triviallly intersecting grids, at the
               ! inside corners, points need to be set by 
               ! interpolation which are near the corner,
               ! so check points within ghostwidth of this one
               !   So if this bdy point is not an AMR bdy pt,
               ! but is within gw of an AMR bdy pt, then we
               ! should still set the points.
               if(ltrace2)write(*,*)'pasteHX:',tmpi,j,k
               tmpi = 1
               nearamr = .false.
 73            continue
               if (j.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .             j.lt. ny-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(nx,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j+tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j+tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .                  j.lt. ny-tmpi                        ) then
                  if (     NINT(chr(nx,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j+tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and. k.gt.1+tmpi .and.
     .                                     k.lt.nz-tmpi      ) then
                  if (     NINT(chr(nx,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and.
     .                  j.lt. ny-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(nx,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j+tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (                   k.gt.1+tmpi .and.
     .                  j.lt. ny-tmpi.and. k.lt.nz-tmpi      ) then
                  if (     NINT(chr(nx,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j+tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j+tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and. k.gt.1+tmpi ) then
                  if (     NINT(chr(nx,j-tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1+tmpi .and. k.lt.nz-tmpi ) then
                  if (     NINT(chr(nx,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j-tmpi,k     )) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if ( j.lt. ny-tmpi.and. k.lt.nz-tmpi) then
                  if (     NINT(chr(nx,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j+tmpi,k+tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k+tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.lt. ny-tmpi.and. k.gt.1+tmpi) then
                  if (     NINT(chr(nx,j+tmpi,k     )) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j+tmpi,k-tmpi)) .eq. AMRBOUNDARY
     .                .or. NINT(chr(nx,j,     k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.gt. 1 +tmpi) then
                  if (     NINT(chr(nx,j-tmpi,k)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (j.lt. ny-tmpi) then
                  if (     NINT(chr(nx,j+tmpi,k)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (k.gt. 1 +tmpi) then
                  if (     NINT(chr(nx,j,k-tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               else if (k.lt. nz-tmpi) then
                  if (     NINT(chr(nx,j,k+tmpi)) .eq. AMRBOUNDARY
     .                ) then
                     nearamr = .true.
                  end if
               end if
               !
               if (nearamr) then
                  ! tmpi is the distance to an AMR point:
                  ! fill in values according to distance:
                  do i = nx, max(nx-gw+tmpi,width), -1
                  !do i = nx, nx-gw+tmpi, -1
                     field(i,j,k) = values(i,j,k)
                  end do
               else
                  if (tmpi .lt. gw-1) then
                  !if (tmpi .lt. gw) then
                     tmpi = tmpi + 1
                     goto 73
                  end if
               end if
               !
            end if
         end do
      end do

 60   continue

      return
      end    ! END: paste_boundaries

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  paste_boundariesNON                                                       cc
cc                                                                            cc
cc      copy data from work space into boundary data of a field               cc
cc      assuming that the necessary interpolation has happened everywhere.    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine paste_boundariesNON(  values,field,     chr,
     *                              mini,  minj,      mink,
     *                              maxi,  maxj,      maxk,
     *                              nx,    ny,        nz,
     *                              r, gw)
      implicit    none
      integer     nx,ny,nz, r, gw
      integer     mini, minj, mink,  maxi, maxj, maxk
      integer     minip,minjp,minkp, maxip,maxjp,maxkp
      real(kind=8)      values(nx,ny,nz),
     *            field(nx,ny,nz), chr(nx,ny,nz)
      include    'chr.inc'
      integer     i,j,k, width
      integer     imid,jmid,kmid,  AMRBOUNDARY
      integer     num
      !
      logical     ltrace
      parameter ( ltrace   = .true. )
      logical     ltrace2
      parameter ( ltrace2  = .false. )

      AMRBOUNDARY = NINT(CHR_amr_bdy)

      if (ltrace) then
         write(*,*) 'paste_boundaries: nx/y/z:   ',nx, ny, nz
         write(*,*) 'paste_boundaries: mini/j/k: ',mini,minj,mink
         write(*,*) 'paste_boundaries: maxi/j/k: ',maxi,maxj,maxk
      end if

      imid = (mini+maxi)/2
      jmid = (minj+maxj)/2
      kmid = (mink+maxk)/2

      !
      ! Low Z Boundary:
      !
      if (mink .ne. 1) goto 10
      width = gw
      if (width .gt. maxk) width = maxk
      if (NINT(chr(imid,jmid,1)) .eq. NINT(CHR_Refl_bdy)) goto 10
      if (ltrace) write(*,*) 'paste_boundaries: Low  Z k: ',1,width
      num = 0
      do k = 1, width
         do j = minj, maxj
            do i = mini, maxi
               if (NINT(chr(i,j,1)) .eq. AMRBOUNDARY) then
                  if(ltrace2)write(*,*)'Low  Z i,j,k:',i,j,k
                  num = num + 1
                  !field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do
      if (ltrace) write(*,*) 'paste_boundaries: num = ',num
 10   continue

      !
      ! High Z Boundary:
      !
      if (maxk .ne. nz) goto 20
      width = nz-gw+1
      if (width .lt. mink) width = mink
      if (ltrace) write(*,*) 'paste_boundaries: High Z k: ',nz, width
      num = 0
      do k = nz, width, -1
         do j = minj, maxj
            do i = mini, maxi
               !if (NINT(chr(i,j,nz)) .ne. CHR_deco_bdy) then
               if (NINT(chr(i,j,nz)) .eq. AMRBOUNDARY) then
                  if(ltrace2)write(*,*)'High Z i,j,k:',i,j,k
                  num = num + 1
                  !field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

      if (ltrace) write(*,*) 'paste_boundaries: num = ',num
 20   continue

      !
      ! Low Y Boundary:
      !
      if (minj .ne. 1) goto 30
      width = gw
      if (width .gt. maxj) width = maxj
      if (NINT(chr(imid,1,kmid)) .eq. NINT(CHR_Refl_bdy)) goto 30
      if (ltrace) write(*,*) 'paste_boundaries: Low  Y j: ',1,width
      num = 0
      do j = 1, width
         do k = mink, maxk
            do i = mini, maxi
               !if (NINT(chr(i,1,k)) .ne. CHR_deco_bdy) then
               if (NINT(chr(i,1,k)) .eq. AMRBOUNDARY) then
                  if(ltrace2)write(*,*)'Low  Y i,j,k:',i,j,k
                  num = num + 1
                  !field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

      if (ltrace) write(*,*) 'paste_boundaries: num = ',num
 30   continue

      !
      ! High Y Boundary:
      !
      if (maxj .ne. ny) goto 40
      width = ny-gw+1
      if (width .lt. minj) width = minj
      if (ltrace) write(*,*) 'paste_boundaries: High Y j: ',ny,width
      num = 0
      do j = ny, width, -1
         do k = mink, maxk
            do i = mini, maxi
               if (NINT(chr(i,ny,k)) .eq. AMRBOUNDARY) then
               !if (NINT(chr(i,ny,k)) .ne. CHR_deco_bdy) then
                  if(ltrace2)write(*,*)'High Y i,j,k:',i,j,k
                  num = num + 1
                  !field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

      if (ltrace) write(*,*) 'paste_boundaries: num = ',num
 40   continue

      !
      ! Low X Boundary:
      !
      if (mini .ne. 1) goto 50
      width = gw
      if (width .gt. maxi) width = maxi
      if (NINT(chr(1,jmid,kmid)) .eq. NINT(CHR_Refl_bdy)) goto 50
      if (ltrace) write(*,*) 'paste_boundaries: Low  X i: ',1,width
      num = 0
      do i = 1, width
         do k = mink, maxk
            do j = minj, maxj
               if (NINT(chr(1,j,k)) .eq. AMRBOUNDARY) then
               !if (NINT(chr(1,j,k)) .ne. CHR_deco_bdy) then
                  if(ltrace2)write(*,*)'Low  X i,j,k:',i,j,k
                  num = num + 1
                  !field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

      if (ltrace) write(*,*) 'paste_boundaries: num = ',num
 50   continue

      !
      ! High X Boundary:
      !
      if (maxi .ne. nx) goto 60
      width = nx-gw+1
      if (width .lt. mini) width = mini
      if (ltrace) write(*,*) 'paste_boundaries: High X i: ',nx,width
      num = 0
      do i = nx, width, -1
         do k = mink, maxk
            do j = minj, maxj
               if (NINT(chr(nx,j,k)) .eq. AMRBOUNDARY) then
               !if (NINT(chr(nx,j,k)) .ne. CHR_deco_bdy) then
                  if(ltrace2)write(*,*)'High X i,j,k:',i,j,k
                  num = num + 1
                  !field(i,j,k) = values(i,j,k)
               end if
            end do
         end do
      end do

      if (ltrace) write(*,*) 'paste_boundaries: num = ',num
 60   continue

      return
      end    ! END: paste_boundariesNON

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  copyFV:                                                                   cc
cc          Take a cell-centered, finite-volume field defined using           cc
cc          the normal vertex-centered memory storage nx ny nz and copy       cc
cc          into an array of dimensions (nx-1)(ny-1)(nz-1).                   cc
cc          This routine simply ignores the last row,column,face, etc         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copyFV( fieldFV, fieldFD, nxFD, nyFD, nzFD)
      implicit   none
      integer    nxFD, nyFD, nzFD
      real(kind=8)     fieldFV(nxFD-1,nyFD-1,nzFD-1)
      real(kind=8)     fieldFD(nxFD,  nyFD,  nzFD  )
      integer    i,j,k

      do k = 1, nzFD-1
         do j = 1, nyFD-1
            do i = 1, nxFD-1
               fieldFV(i,j,k) = fieldFD(i,j,k)
            end do
         end do
      end do

      return
      end    ! END: copyFV


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  restrict_fieldFV:                                                         cc
cc                  Restrict a Child grid to a Parent grid.                   cc
cc               NB: For FINITE VOLUME fields only.                           cc
cc               NB: Do away with ghostwidth. Just get all data possibly      cc
cc                   necessary and let replace_part_field_avgGWa() figure     cc
cc                   out what to do with it and what to ignore.               cc
cc               NB: aware of ghostwidth. Removed unused code from old versioncc
cc               NB: allows for multiple parents.                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine restrict_fieldFV( fieldp,  fieldc, maskc,
     *                             nxp,   nyp,       nzp,
     *                             nxc,   nyc,       nzc,
     *                             min_i, max_i,
     *                             min_j, max_j,
     *                             min_k, max_k,
     *                             r, gw)
      implicit    none
      integer     nxp,nyp,nzp, nxc,nyc,nzc, r, gw
      real(kind=8)      fieldp(nxp,nyp,nzp), fieldc(nxc,nyc,nzc),
     *            maskc(nxc,nyc,nzc)
      include    'chr.inc'
      integer     ic,jc,kc, ip,jp,kp, gwc, min_i,max_i,
     *                    min_j,max_j,min_k,max_k
      include       'largesmall.inc'
      !real(kind=8)      LARGENUMBER
      !parameter       ( LARGENUMBER = 9.d98 )

      logical     ltrace
      parameter ( ltrace = .false. )


      if (ltrace) then
         write(*,*) 'restrict_fieldFV: Start....'
         call field_dump_info(fieldp,nxp,nyp,nzp)
      end if

      do kc = min_k, max_k-r, r
      !do kc = min_k, max_k, r
         kp = (kc - min_k) / r + 1
         do jc = min_j, max_j-r, r
            jp = (jc - min_j) / r + 1
            do ic = min_i, max_i-r, r
               ip = (ic - min_i) / r + 1
               !
               ! Determine if this data should get injected to parent:
               !
               if (
     *      (ic .le. gw      .and.NINT(maskc(1,jc,kc)  ).eq.CHR_amr_bdy)
     *  .or.(jc .le. gw      .and.NINT(maskc(ic,1,kc)  ).eq.CHR_amr_bdy)
     *  .or.(kc .le. gw      .and.NINT(maskc(ic,jc,1)  ).eq.CHR_amr_bdy)
     *  .or.(ic .ge. nxc-gw+1.and.NINT(maskc(nxc,jc,kc)).eq.CHR_amr_bdy)
     *  .or.(jc .ge. nyc-gw+1.and.NINT(maskc(ic,nyc,kc)).eq.CHR_amr_bdy)
     *  .or.(kc .ge. nzc-gw+1.and.NINT(maskc(ic,jc,nzc)).eq.CHR_amr_bdy)
     *      )  then
                  fieldp(ip,jp,kp) = LARGENUMBER
               else
                  !fieldp(ip,jp,kp) = fieldc(ic,  jc,  kc  )
                  ! Average child cells to get value to be injected:
                  fieldp(ip,jp,kp) = 0.125d0*(
     *                                        + fieldc(ic,  jc,  kc  )
     *                                        + fieldc(ic+1,jc,  kc  )
     *                                        + fieldc(ic,  jc+1,kc  )
     *                                        + fieldc(ic,  jc,  kc+1)
     *                                        + fieldc(ic+1,jc+1,kc  )
     *                                        + fieldc(ic,  jc+1,kc+1)
     *                                        + fieldc(ic+1,jc,  kc+1)
     *                                        + fieldc(ic+1,jc+1,kc+1) )
               end if
            end do
         end do
      end do

      if (ltrace) then
         call field_dump_info(fieldp,nxp,nyp,nzp)
         write(*,*) 'restrict_fieldFV: i/j/kp =',ip,jp,kp
         write(*,*) 'restrict_fieldFV: '
         write(*,*) 'restrict_fieldFV: r   = ', r
         write(*,*) 'restrict_fieldFV: gw  = ', gw
         write(*,*) 'restrict_fieldFV: ',1+gwc*r,nzc-gwc*r
         write(*,*) 'restrict_fieldFV: ip,jp,kp = ',ip,jp,kp
         write(*,*) 'restrict_fieldFV: nx/y/zp = ',nxp,nyp,nzp
         write(*,*) 'restrict_fieldFV: nx/y/zc = ',nxc,nyc,nzc
         write(*,*) 'restrict_fieldFV: i/j/kc  = ',ic,jc,kc
         write(*,*) 'restrict_fieldFV: min_i/j/k ',min_i,min_j,min_k
         write(*,*) 'restrict_fieldFV: max_i/j/k ',max_i,max_j,max_k
      end if

      return
      end    ! END: restrict_fieldFV

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  replace_part_field_avgFV:                                                 cc
cc                  NB: For FINITE VOLUME fields only.                        cc
cc                  Replace part of a field with entirity of another.         cc
cc                  Ignores boundaries!                                       cc
cc                  Averages according to the mask function so that           cc
cc                  overlapping grids can be accomadated.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine replace_part_field_avgFV(fieldp, maskp,fieldc,
     *                             ax,    ay,        az,
     *                             nxp,   nyp,       nzp,
     *                             nxs,   nys,       nzs,
     *                             nxc,   nyc,       nzc,
     *                             min_ic,max_ic,
     *                             min_jc,max_jc,
     *                             min_kc,max_kc,
     *                             r, gw, bwidth   )
      implicit    none
      integer     nxp,nyp,nzp, nxc,nyc,nzc, ax,ay,az, gw,r,
     *            nxs,nys,nzs, bwidth
      integer     min_ip,min_jp,min_kp
      integer     min_ic,max_ic, min_jc,max_jc,min_kc,max_kc
      real(kind=8)      fieldp(nxp,nyp,nzp), maskp(nxp,nyp,nzp),
     *            fieldc(nxs,nys,nzs)
      include    'chr.inc'
      integer     i,j,k,   ip,jp,kp, ic,jc,kc, gwc
      integer     gwc_i, gwc_j, gwc_k
      integer     ib,ie, jb,je, kb,ke, n
c     logical     dd_bdy
c     logical     double_equal
c     external    double_equal
      include       'largesmall.inc'
      !real(kind=8)      LARGENUMBER
      !parameter       ( LARGENUMBER = 9.d98   )
      logical           ltrace
      parameter       ( ltrace      = .false. )


      do kc = min_kc, max_kc-r, r
         kp = (kc - min_kc) / r + 1
         k  = kp + az - 1
         do jc = min_jc, max_jc-r, r
            jp = (jc - min_jc) / r + 1
            j  = jp + ay - 1
            do ic = min_ic, max_ic-r, r
               ip = (ic - min_ic) / r + 1
               i  = ip + ax - 1
               if ( fieldc(ip,jp,kp).lt.LARGE_M_SMALL ) then
               !if ( fieldc(ip,jp,kp).lt.LARGENUMBER ) then
                  fieldp(i,j,k) = fieldc(ip,jp,kp)
               end if
c              n  = - NINT(maskp(i,j,k))
               !
               ! Since FV fields chop off the last point,
               ! you need to look at the "advanced" mask point
               ! to see if you're at DD boundary there:
               !
c              dd_bdy =    NINT(maskp(i,j,k)  ).eq. CHR_deco_bdy
c              if (.not. dd_bdy .and. i.ne.nxp)
c    *         dd_bdy =    NINT(maskp(i+1,j,k)).eq. CHR_deco_bdy
c              if (.not. dd_bdy .and. j.ne.nyp)
c    *         dd_bdy =    NINT(maskp(i,j+1,k)).eq. CHR_deco_bdy
c              if (.not. dd_bdy .and. k.ne.nzp)
c    *         dd_bdy =    NINT(maskp(i,j,k+1)).eq. CHR_deco_bdy

c              dd_bdy =    NINT(maskp(i,j,k)  ).eq. CHR_deco_bdy .or.
c    *                     NINT(maskp(i+1,j,k)).eq. CHR_deco_bdy .or.
c    *                     NINT(maskp(i,j+1,k)).eq. CHR_deco_bdy .or.
c    *                     NINT(maskp(i,j,k+1)).eq. CHR_deco_bdy
               !
               ! Needs to be in agreement with how FD fields
               ! are handled with inc_mask_overlapb():
               !
c              if ( n .ge. 0 .and.
c    *              fieldc(ip,jp,kp).lt.LARGENUMBER .and.
c    *              .not. dd_bdy ) then
c              if ( n .ge. 0 .and.
c    *              fieldc(ip,jp,kp).lt.LARGENUMBER) then
c                 fieldp(i,j,k) =
c    *                      (     n * fieldp(i,j,k)
c    *                        +       fieldc(ip,jp,kp)
c    *                       ) / ( 1.d0 + n)
c              else if ( fieldc(ip,jp,kp).lt.LARGENUMBER .and.
c    *                   dd_bdy ) then
c                 !
c                 ! These represent points on the domain deco.
c                 ! boundary. No averaging done here, but
c                 ! we still want to restrict to keep the coarse
c                 ! grid properly syncing to child:
c                 !
c                 fieldp(i,j,k) = fieldc(ip,jp,kp)
c              end if
c     if (i.eq.6.and.j.eq.32.and.k.eq.32)
c    *write(*,*) 'After: ',6,fieldp(6,32,32),nint(maskp(6,32,32))
            end do
         end do
      end do

      if (ltrace) then
         write(*,*) 'replace_part_field_avgFV: ax/y/z =',ax,ay,az
         write(*,*) 'replace_part_field_avgFV: i/j/kp =',ip,jp,kp
         write(*,*) 'replace_part_field_avgFV: n[xyz]p=',nxp,nyp,nzp
         write(*,*) 'replace_part_field_avgFV: n[xyz]c=',nxc,nyc,nzc
         write(*,*) 'replace_part_field_avgFV: n[xyz]s=',nxs,nys,nzs
         write(*,*) 'replace_part_field_avgFV:mini/j/kc',
     *                                           min_ic,min_jc,min_kc
         write(*,*) 'replace_part_field_avgFV:maxi/j/kc',
     *                                           max_ic,max_jc,max_kc
         write(*,*) 'replace_part_field_avgFV:    gw  =',gw
         write(*,*) 'replace_part_field_avgFV:    r   =',r
         write(*,*) 'replace_part_field_avgFV:    gwc =',gwc
         write(*,*)
     * 'replace_part_field_avgFV:gwci/j/k=',gwc_i,gwc_j,gwc_k
      end if


      return
      end        ! END: subroutine replace_part_field_avgFV


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  tre_subtract:                                                             cc
cc                                                                            cc
cc      With interpolated (in space) data from parent, subtract from          cc
cc      resident field and add to the error field for the purposes of         cc
cc      computing the truncation error estimation.                            cc
cc          numfields--- the number of fields for which this is computed      cc
cc                       and the error is divided/rescaled by this number     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tre_subtract(      parent, field,    error, chr,weight,
     *                              mini,  minj,      mink,
     *                              maxi,  maxj,      maxk,
     *                              nx,    ny,        nz,
     *                              numfields, r, gw)
      implicit    none
      integer     nx,ny,nz, r, gw
      integer     mini, minj, mink,  maxi, maxj, maxk
      integer     minip,minjp,minkp, maxip,maxjp,maxkp, numfields
      real(kind=8)      weight
      real(kind=8)      parent(nx,ny,nz),
     *            field(nx,ny,nz), error(nx,ny,nz), chr(nx,ny,nz)
      include    'chr.inc'
      real(kind=8)      tmp
      integer     i,j,k, ichr
      logical  isanan
      external isanan
      !
      logical     ltrace
      parameter ( ltrace  = .false. )

      if (ltrace) then
         write(*,*) 'tre_subtract: nx/y/z:   ',nx, ny, nz
         write(*,*) 'tre_subtract: mini/j/k: ',mini,minj,mink
         write(*,*) 'tre_subtract: maxi/j/k: ',maxi,maxj,maxk
      end if

      do k = mink, maxk
         do j = minj, maxj
            do i = mini, maxi
               !
               ichr = NINT(chr(i,j,k))
               !
               if (ichr .ne. CHR_deco_bdy .and. ichr.ne.CHR_amr_bdy)then
                  tmp          = abs(parent(i,j,k)-field(i,j,k))
     *                             / (1.d0*numfields)
                  if (isanan(tmp)) then
                  if(ltrace)write(*,*)'tre_subtract:ignoring Nan in TRE'
                  else
                     error(i,j,k) = error(i,j,k) + weight * tmp
                  end if
                  !error(i,j,k) = error(i,j,k) + tmp
               end if
               !
            end do
         end do
      end do

      return
      end    ! END: tre_subtract

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine myapply_outflow(f, chr, nx, ny, nz, ng)
      implicit   none
      include 'chr.inc'
      integer                        :: nx, ny, nz, ng
      real*8 f(nx,ny,nz),chr(nx,ny,nz)
      ! local vars
      integer                        :: i, j, k, count, mychr
      logical     ltrace
      parameter ( ltrace  = .true.)
      logical     ltrace2
      parameter ( ltrace2 = .false.)
      logical     previous
      parameter ( previous = .false.)

      count = 0
      if (ltrace) write(*,*) 'myapply_outflow: nx/y/z/g: ',nx,ny,nz,ng
      ! xmin
      do k = 1, nz
      do j = 1, ny
       mychr = NINT(chr(1,j,k))
       if (previous) mychr = 0
       if (  mychr .ne. NINT(CHR_amr_bdy)  .and.
     *       mychr .ne. NINT(CHR_deco_bdy) .and.
     *       mychr .ne. NINT(CHR_Refl_bdy)  ) then
           do i = 1, ng
             f(i,j,k) = f(ng+1,j,k)
           end do
       else
           count = count + 1
           if(ltrace2)write(*,*) 'myapply: i=1:',j,k,mychr
       end if
      end do
      end do

      ! xmax
      do k = 1, nz
      do j = 1, ny
       mychr = NINT(chr(nx,j,k))
       if (previous) mychr = 0
       if (  mychr .ne. NINT(CHR_amr_bdy)  .and.
     *       mychr .ne. NINT(CHR_deco_bdy) .and.
     *       mychr .ne. NINT(CHR_Refl_bdy)  ) then
           do i = nx-ng+1, nx
             f(i,j,k) = f(nx-ng,j,k)
           end do
       else
           if(ltrace2)write(*,*) 'myapply: i=nx:',j,k,mychr
            count = count + 1
        end if
      end do
      end do

      ! ymin
      do k = 1, nz
      do i = 1, nx
       mychr = NINT(chr(i,1,k))
       if (previous) mychr = 0
       if (  mychr .ne. NINT(CHR_amr_bdy)  .and.
     *       mychr .ne. NINT(CHR_deco_bdy) .and.
     *       mychr .ne. NINT(CHR_Refl_bdy)  ) then
          do j = 1, ng
             f(i,j,k) = f(i,ng+1,k)
          end do
       else
          if(ltrace2)write(*,*) 'myapply: j=1:',i,k,mychr
          count = count + 1
       end if
      end do
      end do

      ! ymax
      do k = 1, nz
      do i = 1, nx
        mychr = NINT(chr(i,ny,k))
       if (previous) mychr = 0
       if (  mychr .ne. NINT(CHR_amr_bdy)  .and.
     *       mychr .ne. NINT(CHR_deco_bdy) .and.
     *       mychr .ne. NINT(CHR_Refl_bdy)  ) then
          do j = ny-ng+1, ny
             f(i,j,k) = f(i,ny-ng,k)
          end do
        else
            if(ltrace2)write(*,*) 'myapply: j=ny:',i,k,mychr
            count = count + 1
        end if
      end do
      end do

      ! zmin
      do j = 1, ny
      do i = 1, nx
        mychr = NINT(chr(i,j,1))
       if (previous) mychr = 0
       if (  mychr .ne. NINT(CHR_amr_bdy)  .and.
     *       mychr .ne. NINT(CHR_deco_bdy) .and.
     *       mychr .ne. NINT(CHR_Refl_bdy)  ) then
          do k = 1, ng
             f(i,j,k) = f(i,j,ng+1)
          end do
        else
            if(ltrace2)write(*,*) 'myapply: k=1:',i,j,mychr
            count = count + 1
        end if
      end do
      end do

      ! zmax
      do j = 1, ny
      do i = 1, nx
        mychr = NINT(chr(i,j,nz))
       if (previous) mychr = 0
       if (  mychr .ne. NINT(CHR_amr_bdy)  .and.
     *       mychr .ne. NINT(CHR_deco_bdy) .and.
     *       mychr .ne. NINT(CHR_Refl_bdy)  ) then
          do k = nz-ng+1, nz
             f(i,j,k) = f(i,j,nz-ng)
          end do
        else
            if(ltrace2)write(*,*) 'myapply: k=nz:',i,j,mychr
            count = count + 1
        end if
      end do
      end do

      if (previous)write(*,*)'myapply_outflow: previous=.true.'
      if(ltrace)write(*,*)'myapply_outflow: %bound untouched: ',
     *                      100*count/(2*(nx*ny+ny*nz+nx*nz))

      return
      end 



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc    Enforce an assumed symmetry about z=0 by populating z<=0 points         cc
cc    according to whether a field is even or odd about the z=0 plane.        cc
cc    odd_or_even =  1 -- EVEN                                                cc
cc                = -1 -- ODD       (as specified in had/include/fields.inc   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine apply_symmetryz(f, z, odd_or_even, nx, ny, nz)
      implicit   none
      integer    nx,ny,nz, odd_or_even
      real*8     f(nx,ny,nz), z(nz)
      ! local vars
      integer    i,j,k, kzero
      integer    find_index
      external   find_index

      kzero =  find_index(z, 0.d0, nz )
      
      if (kzero .lt. 0) then
         write(*,*) 'apply_symmetry: PROBLEM, z=0 not present'
         write(*,*) 'apply_symmetry: min/max z ',z(1),z(nz)
         return
      end if

      do k =1, kzero-1
         do j = 1, ny
            do i = 1, nx
               f(i,j,k) = (1.d0*odd_or_even) * f(i,j,2*kzero-k)
            end do
         end do
      end do
      
      if (odd_or_even .lt. 0) then
         k = kzero
         do j = 1, ny
            do i = 1, nx
               f(i,j,k) = 0.d0
            end do
         end do
      end if

      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc    Enforce an assumed symmetry about z=0 by populating z<=0 points         cc
cc    according to whether a field is even or odd about the z=0 plane.        cc
cc    odd_or_even =  1 -- EVEN                                                cc
cc                = -1 -- ODD       (as specified in had/include/fields.inc   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine apply_symmetryy(f, y, odd_or_even, nx, ny, nz)
      implicit   none
      integer    nx,ny,nz, odd_or_even
      real*8     f(nx,ny,nz), y(ny)
      ! local vars
      integer    i,j,k, jzero
      integer    find_index
      external   find_index

      jzero =  find_index(y, 0.d0, ny )
      
      if (jzero .lt. 0) then
         write(*,*) 'apply_symmetryy: PROBLEM, y=0 not present'
         write(*,*) 'apply_symmetryy: min/max y ',y(1),y(nz)
         return
      end if

      do j =1, jzero-1
         do k = 1, nz
            do i = 1, nx
               f(i,j,k) = (1.d0*odd_or_even) * f(i,2*jzero-j,k)
            end do
         end do
      end do
      
      if (odd_or_even .lt. 0) then
         j = jzero
         do k = 1, nz
            do i = 1, nx
               f(i,j,k) = 0.d0
            end do
         end do
      end if

      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc    Enforce an assumed symmetry about z=0 by populating z<=0 points         cc
cc    according to whether a field is even or odd about the z=0 plane.        cc
cc    odd_or_even =  1 -- EVEN                                                cc
cc                = -1 -- ODD       (as specified in had/include/fields.inc   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine apply_symmetryx(f, x, odd_or_even, nx, ny, nz)
      implicit   none
      integer    nx,ny,nz, odd_or_even
      real*8     f(nx,ny,nz), x(nx)
      ! local vars
      integer    i,j,k, izero
      integer    find_index
      external   find_index

      izero =  find_index(x, 0.d0, nx )
      
      if (izero .lt. 0) then
         write(*,*) 'apply_symmetryx: PROBLEM, x=0 not present'
         write(*,*) 'apply_symmetryx: min/max x ',x(1),x(nx)
         return
      end if

      do i =1, izero-1
         do k = 1, nz
            do j = 1, ny
               f(i,j,k) = (1.d0*odd_or_even) * f(2*izero-i,j,k)
            end do
         end do
      end do
      
      if (odd_or_even .lt. 0) then
         i = izero
         do k = 1, nz
            do j = 1, ny
               f(i,j,k) = 0.d0
            end do
         end do
      end if

      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  deriv34_x:                                                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deriv34_x( du, u, h, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     du(nx,ny,nz), u(nx,ny,nz), h
      integer    i,j,k
      real*8    d00,d01,d02,d03, 
     &    d10,d11,d12,d13,d14,d15, 
     &    d20,d21,d22,d23,d24,d25, 
     &    d30,d31,d32,d33,d34,d35,d36, 
     &    d40,d41,d42,d43,d44,d45,d46, 
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


       do k = 1, nz
       do j = 1, ny
         du( 1,j,k) =  (  d00*u( 1,  j,k) + d01*u(   2,j,k)
     .                   +d02*u( 3,  j,k) + d03*u(   4,j,k)  )/h
         du(nx,j,k) = -(  d00*u(nx,  j,k) + d01*u(nx-1,j,k)
     .                   +d02*u(nx-2,j,k) + d03*u(nx-3,j,k)  )/h
         !
         du(   2,j,k) =  (  d10*u(   1,j,k) + d11*u(   2,j,k)
     .                     +d12*u(   3,j,k) + d13*u(   4,j,k)
     .                     +d14*u(   5,j,k) + d15*u(   6,j,k)  )/h
         du(nx-1,j,k) = -(  d10*u(  nx,j,k) + d11*u(nx-1,j,k)
     .                     +d12*u(nx-2,j,k) + d13*u(nx-3,j,k)
     .                     +d14*u(nx-4,j,k) + d15*u(nx-5,j,k)  )/h
         !
         du(   3,j,k) =  (  d20*u(   1,j,k) + d21*u(   2,j,k)
     .                     +d22*u(   3,j,k) + d23*u(   4,j,k)
     .                     +d24*u(   5,j,k) + d25*u(   6,j,k)  )/h
         du(nx-2,j,k) = -(  d20*u(  nx,j,k) + d21*u(nx-1,j,k)
     .                     +d22*u(nx-2,j,k) + d23*u(nx-3,j,k)
     .                     +d24*u(nx-4,j,k) + d25*u(nx-5,j,k)  )/h
         !
         du(   4,j,k) =  (  d30*u(   1,j,k) + d31*u(   2,j,k)
     .                     +d32*u(   3,j,k) + d33*u(   4,j,k)
     .                     +d34*u(   5,j,k) + d35*u(   6,j,k)
     .                     +d36*u(   7,j,k)                    )/h
         du(nx-3,j,k) = -(  d30*u(  nx,j,k) + d31*u(nx-1,j,k)
     .                     +d32*u(nx-2,j,k) + d33*u(nx-3,j,k)
     .                     +d34*u(nx-4,j,k) + d35*u(nx-5,j,k)
     .                     +d36*u(nx-6,j,k)                    )/h
         !
         du(   5,j,k) =  (  d40*u(   1,j,k) + d41*u(   2,j,k)
     .                     +d42*u(   3,j,k) + d43*u(   4,j,k)
     .                     +d44*u(   5,j,k) + d45*u(   6,j,k)
     .                     +d46*u(   7,j,k)                    )/h
         du(nx-4,j,k) = -(  d40*u(  nx,j,k) + d41*u(nx-1,j,k)
     .                     +d42*u(nx-2,j,k) + d43*u(nx-3,j,k)
     .                     +d44*u(nx-4,j,k) + d45*u(nx-5,j,k)
     .                     +d46*u(nx-6,j,k)                    )/h
         !
       end do
       end do

       do k = 1, nz
       do j = 1, ny
       do i = 6, nx-5
          du(i,j,k) = ( -d2*u(i-2,j,k) - d1*u(i-1,j,k)
     .                  +d1*u(i+1,j,k) + d2*u(i+2,j,k)  )/h
       end do
       end do
       end do

      return
      end    ! END: deriv34_x

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  deriv34_y:                                                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deriv34_y( du, u, h, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     du(nx,ny,nz), u(nx,ny,nz), h
      integer    i,j,k
      real*8    d00,d01,d02,d03, 
     &    d10,d11,d12,d13,d14,d15, 
     &    d20,d21,d22,d23,d24,d25, 
     &    d30,d31,d32,d33,d34,d35,d36, 
     &    d40,d41,d42,d43,d44,d45,d46, 
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


       do k = 1, nz
       do j = 1, nx
         du(j, 1,k) =  (  d00*u(j, 1,  k) + d01*u(j,   2,k)
     .                   +d02*u(j, 3,  k) + d03*u(j,   4,k)  )/h
         du(j,ny,k) = -(  d00*u(j,ny,  k) + d01*u(j,ny-1,k)
     .                   +d02*u(j,ny-2,k) + d03*u(j,ny-3,k)  )/h
         !
         du(j,   2,k) =  (  d10*u(j,   1,k) + d11*u(j,   2,k)
     .                     +d12*u(j,   3,k) + d13*u(j,   4,k)
     .                     +d14*u(j,   5,k) + d15*u(j,   6,k)  )/h
         du(j,ny-1,k) = -(  d10*u(j,  ny,k) + d11*u(j,ny-1,k)
     .                     +d12*u(j,ny-2,k) + d13*u(j,ny-3,k)
     .                     +d14*u(j,ny-4,k) + d15*u(j,ny-5,k)  )/h
         !
         du(j,   3,k) =  (  d20*u(j,   1,k) + d21*u(j,   2,k)
     .                     +d22*u(j,   3,k) + d23*u(j,   4,k)
     .                     +d24*u(j,   5,k) + d25*u(j,   6,k)  )/h
         du(j,ny-2,k) = -(  d20*u(j,  ny,k) + d21*u(j,ny-1,k)
     .                     +d22*u(j,ny-2,k) + d23*u(j,ny-3,k)
     .                     +d24*u(j,ny-4,k) + d25*u(j,ny-5,k)  )/h
         !
         du(j,   4,k) =  (  d30*u(j,   1,k) + d31*u(j,   2,k)
     .                     +d32*u(j,   3,k) + d33*u(j,   4,k)
     .                     +d34*u(j,   5,k) + d35*u(j,   6,k)
     .                     +d36*u(j,   7,k)                    )/h
         du(j,ny-3,k) = -(  d30*u(j,  ny,k) + d31*u(j,ny-1,k)
     .                     +d32*u(j,ny-2,k) + d33*u(j,ny-3,k)
     .                     +d34*u(j,ny-4,k) + d35*u(j,ny-5,k)
     .                     +d36*u(j,ny-6,k)                    )/h
         !
         du(j,   5,k) =  (  d40*u(j,   1,k) + d41*u(j,   2,k)
     .                     +d42*u(j,   3,k) + d43*u(j,   4,k)
     .                     +d44*u(j,   5,k) + d45*u(j,   6,k)
     .                     +d46*u(j,   7,k)                    )/h
         du(j,ny-4,k) = -(  d40*u(j,  ny,k) + d41*u(j,ny-1,k)
     .                     +d42*u(j,ny-2,k) + d43*u(j,ny-3,k)
     .                     +d44*u(j,ny-4,k) + d45*u(j,ny-5,k)
     .                     +d46*u(j,ny-6,k)                    )/h
         !
       end do
       end do

       do k = 1, nz
       do i = 6, ny-5
       do j = 1, nx
          du(j,i,k) = ( -d2*u(j,i-2,k) - d1*u(j,i-1,k)
     .                  +d1*u(j,i+1,k) + d2*u(j,i+2,k)  )/h
       end do
       end do
       end do

      return
      end    ! END: deriv34_y

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  deriv34_z:                                                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deriv34_z( du, u, h, nx, ny, nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)     du(nx,ny,nz), u(nx,ny,nz), h
      integer    i,j,k
      real*8    d00,d01,d02,d03, 
     &    d10,d11,d12,d13,d14,d15, 
     &    d20,d21,d22,d23,d24,d25, 
     &    d30,d31,d32,d33,d34,d35,d36, 
     &    d40,d41,d42,d43,d44,d45,d46, 
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


       do k = 1, ny
       do j = 1, nx
         du(j,k, 1) =  (  d00*u(j,k, 1  ) + d01*u(j,k,   2)
     .                   +d02*u(j,k, 3  ) + d03*u(j,k,   4)  )/h
         du(j,k,nz) = -(  d00*u(j,k,nz  ) + d01*u(j,k,nz-1)
     .                   +d02*u(j,k,nz-2) + d03*u(j,k,nz-3)  )/h
         !
         du(j,k,   2) =  (  d10*u(j,k,   1) + d11*u(j,k,   2)
     .                     +d12*u(j,k,   3) + d13*u(j,k,   4)
     .                     +d14*u(j,k,   5) + d15*u(j,k,   6)  )/h
         du(j,k,nz-1) = -(  d10*u(j,k,  nz) + d11*u(j,k,nz-1)
     .                     +d12*u(j,k,nz-2) + d13*u(j,k,nz-3)
     .                     +d14*u(j,k,nz-4) + d15*u(j,k,nz-5)  )/h
         !
         du(j,k,   3) =  (  d20*u(j,k,   1) + d21*u(j,k,   2)
     .                     +d22*u(j,k,   3) + d23*u(j,k,   4)
     .                     +d24*u(j,k,   5) + d25*u(j,k,   6)  )/h
         du(j,k,nz-2) = -(  d20*u(j,k,  nz) + d21*u(j,k,nz-1)
     .                     +d22*u(j,k,nz-2) + d23*u(j,k,nz-3)
     .                     +d24*u(j,k,nz-4) + d25*u(j,k,nz-5)  )/h
         !
         du(j,k,   4) =  (  d30*u(j,k,   1) + d31*u(j,k,   2)
     .                     +d32*u(j,k,   3) + d33*u(j,k,   4)
     .                     +d34*u(j,k,   5) + d35*u(j,k,   6)
     .                     +d36*u(j,k,   7)                    )/h
         du(j,k,nz-3) = -(  d30*u(j,k,  nz) + d31*u(j,k,nz-1)
     .                     +d32*u(j,k,nz-2) + d33*u(j,k,nz-3)
     .                     +d34*u(j,k,nz-4) + d35*u(j,k,nz-5)
     .                     +d36*u(j,k,nz-6)                    )/h
         !
         du(j,k,   5) =  (  d40*u(j,k,   1) + d41*u(j,k,   2)
     .                     +d42*u(j,k,   3) + d43*u(j,k,   4)
     .                     +d44*u(j,k,   5) + d45*u(j,k,   6)
     .                     +d46*u(j,k,   7)                    )/h
         du(j,k,nz-4) = -(  d40*u(j,k,  nz) + d41*u(j,k,nz-1)
     .                     +d42*u(j,k,nz-2) + d43*u(j,k,nz-3)
     .                     +d44*u(j,k,nz-4) + d45*u(j,k,nz-5)
     .                     +d46*u(j,k,nz-6)                    )/h
         !
       end do
       end do

       do i = 6, nz-5
       do k = 1, ny
       do j = 1, nx
          du(j,k,i) = ( -d2*u(j,k,i-2) - d1*u(j,k,i-1)
     .                  +d1*u(j,k,i+1) + d2*u(j,k,i+2)  )/h
       end do
       end do
       end do

      return
      end    ! END: deriv34_z

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  add_diss: compute and dissipation to the field                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine add_diss( u, uRK, work, eps, h, nx,ny,nz)
      implicit   none
      integer    nx, ny, nz
      real(kind=8)  u(nx,ny,nz), uRK(nx,ny,nz), eps, h
      real(kind=8)  work(nx,ny,nz)
      integer    i,j,k

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if (   i.gt.4 .and. i.lt.nx-4 .and.
     .          j.gt.4 .and. j.lt.ny-4 .and.
     .          k.gt.4 .and. k.lt.nz-4  ) then
            work(i,j,k) = (
     .     ! x-direction:
     .                   u(i+3,j,k) + u(i-3,j,k)
     .           -6.0d0*(u(i+2,j,k) + u(i-2,j,k))
     .         + 15.0d0*(u(i+1,j,k) + u(i-1,j,k))
     .          -20.0d0* u(i,j,k)
     .     ! y-direction:
     .                 + u(i,j+3,k) + u(i,j-3,k)
     .           -6.0d0*(u(i,j+2,k) + u(i,j-2,k))
     .         + 15.0d0*(u(i,j+1,k) + u(i,j-1,k))
     .          -20.0d0* u(i,j,k)
     .     ! z-direction:
     .                 + u(i,j,k+3) + u(i,j,k-3)
     .           -6.0d0*(u(i,j,k+2) + u(i,j,k-2))
     .         + 15.0d0*(u(i,j,k+1) + u(i,j,k-1))
     .          -20.0d0* u(i,j,k)
     .                     ) / h
         else
            ! not implementing these right now:
            work(i,j,k) = 0.d0
         end if
         uRK(i,j,k) = uRK(i,j,k) + eps*work(i,j,k)
      end do
      end do
      end do

      return
      end    ! END: add_diss

