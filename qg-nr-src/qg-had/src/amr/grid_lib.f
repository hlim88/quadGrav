cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_kill_children:                                                       cc
cc                                                                            cc
cc      During regridding, grids need to be turned off or killed              cc
cc      but still kept in memory so that their information can                cc
cc      be re-used or even possibly for them to come back to life             cc
cc      if the clusterer calls for a grid with identical dimensions           cc
cc      to be produced. This routine turns off the children of a given        cc
cc      grid. Another routine will then come along to get rid of              cc
cc      them from memory when the time comes.                                 cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function  grid_kill_children(grid)
      implicit none
      include 'glob.inc'
      integer  grid
      integer  gi,     previous,   next, sibling, lev_i

      logical  grid_return_existence
      integer  grid_return_child,
     *         grid_return_parent,
     *         grid_return_sibling
      external grid_return_existence,
     *         grid_return_child,
     *         grid_return_parent,
     *         grid_return_sibling

      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) write(*,*) 'kill_children: Children of grid: ',grid

      grid_kill_children = 0

      if (.not.grid_return_existence(grid)) then
        write(*,*) 'kill_children: WARNING!!!'
        write(*,*) 'kill_children: grid does not exist, ',grid
        return
      end if

      gi       = grid_return_child(grid)

      if (.not.grid_return_existence(grid)) then
         write(*,*) 'kill_children: Grid has no children',grid,gi
         return
      end if

      !
      ! All children of the inputted grid should be consecutive:
      !
  9   if (grid_return_parent(gi) .eq. grid) then
 10      call grid_mk_dead(gi)
         grid_kill_children = grid_kill_children + 1
         if (ltrace) write(*,*) '   Marking as dead grid: ',gi
         gi                 = grid_return_sibling(gi)
         if (grid_return_existence(gi)) then
            if (grid_return_parent(gi) .eq. grid) goto 10
         end if
      end if

      if (ltrace) write(*,*) 'kill_children: Done.'

      return
      end       ! END: function grid_kill_children()

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_is_local:                                                            cc
cc                  Is the owner of this grid the same processor asking?      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_is_local( gi )
      implicit    none
      integer     gi
      include    'mpif.h'
      include    'mpi_stuff.inc'
      integer     grid_return_owner
      external    grid_return_owner
      
      grid_is_local = myid .eq. grid_return_owner(gi)

      return
      end       ! END: grid_is_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_is_childof:                                                          cc
cc                  Does grid gi exist and is the child of "gp"?              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_is_childof( gi, gp )
      implicit    none
      integer     gi, gp
c     include    'grid.inc'
      integer     grid_return_parent
      logical     grid_return_existence
      logical     ltrace
      parameter ( ltrace = .false. )
      
      grid_is_childof = grid_return_existence(gi)

      if (grid_is_childof) then
         grid_is_childof = grid_return_parent(gi).eq.gp
      end if

      return
      end       ! END: grid_is_childof

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_ana:                                                          cc
cc                  Returns the value of gfunc_ana() for the given field      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_ana(field)
      implicit    none
      integer    field
      include    'grid.inc'

      grid_return_ana = gfunc_ana(field)

      return
      end       ! END: grid_return_ana

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_numgfuncs:                                                    cc
cc                  Returns the number of grid functions.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_numgfuncs()
      implicit    none
      include    'grid.inc'

      grid_return_numgfuncs = num_gfuncs

      return
      end       ! END: grid_return_numgfuncs

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_finest:                                                       cc
cc                  Returns pointer to (a) grid at finest level.              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_finest()
      implicit    none
c     include    'grid.inc'
      include    'glob.inc'
      integer     level_return_finest,
     *            level_return_start

c     grid_return_finest = Levelp(level_return_finest(Levelp,maxlev))
      grid_return_finest = level_return_start(
     *                            level_return_finest(Levelp,maxlev))

      if(grid_return_finest.le.0) then
         write(*,*)'grid_return_finest: Problem Warning:******'
         write(*,*)'grid_return_finest:grid:',grid_return_finest
         write(*,*)'grid_return_finest:level:',
     .                        level_return_finest(Levelp,maxlev)
      end if

      return
      end       ! END: grid_return_finest

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_at_t0:                                                               cc
cc               Returns true if grid is at time t=0.                         cc
cc               To determine this look at lev_count() for all                cc
cc               higher levels. If all are zero, then return true.            cc
cc               If any are nonzero, return false.                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_at_t0(grid)
      implicit    none
      integer     grid
c     include    'grid.inc'
      include    'glob.inc'
      integer     grid_return_level
      integer     level

      do level = grid_return_level(grid), 0, -1
         if ( lev_count(level) .ne. 0 ) then
            grid_at_t0 = .false.
            return
         end if
      end do

      grid_at_t0 = .true.

      return
      end       ! END: grid_at_t0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_init_chr:                                                            cc
cc                 Initialize characteristic chr() array for a given grid.    cc
cc                 Valid values taken from "chr.inc".                         cc
cc                 chr() values used in update routines applying              cc
cc                 different (or none at all) numerical updates depending     cc
cc                 on type of grid point described by characteristic value.   cc
cc             NB: now handles assumptions of symmetry.                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_init_chr( chr_array, minx, miny, minz,
     *                            hi, nx,ny,nz )
      implicit    none
      include    'chr.inc'
      include    'glob.inc'
      integer    nx, ny, nz
      real(kind=8)  chr_array(nx,ny,nz), minx, miny, minz, hi
      real(kind=8)  maxx, maxy, maxz
      integer    i,  j,  k
      logical     double_equal
      external    double_equal
      logical     atminx, atmaxx,
     *            atminy, atmaxy,
     *            atminz, atmaxz
      logical     ltrace
      parameter ( ltrace  = .false. )

      if (ltrace) write(*,*) 'grid_init_chr: Enter'

      maxx = minx + (nx-1.d0)*hi
      maxy = miny + (ny-1.d0)*hi
      maxz = minz + (nz-1.d0)*hi

      atminx = double_equal(minx,minx0)
      atminy = double_equal(miny,miny0)
      atminz = double_equal(minz,minz0)

      atmaxx = double_equal(maxx,maxx0)
      atmaxy = double_equal(maxy,maxy0)
      atmaxz = double_equal(maxz,maxz0)

      if (ltrace) then
         write(*,*) 'grid_init_chr: min/axx0    ',minx0,maxx0
         write(*,*) 'grid_init_chr: min/axy0    ',miny0,maxy0
         write(*,*) 'grid_init_chr: min/axz0    ',minz0,maxz0
         write(*,*) 'grid_init_chr: min/axx:    ',minx,maxx
         write(*,*) 'grid_init_chr: min/axy:    ',miny,maxy
         write(*,*) 'grid_init_chr: min/axz:    ',minz,maxz
         write(*,*) 'grid_init_chr: atmin/maxx: ',atminx,atmaxx
         write(*,*) 'grid_init_chr: atmin/maxy: ',atminy,atmaxy
         write(*,*) 'grid_init_chr: atmin/maxz: ',atminz,atmaxz
         write(*,13) 'grid_init_chr: diffs min x: ',abs(minx0-minx0)
         write(*,13) 'grid_init_chr: diffs min y: ',abs(miny0-miny0)
         write(*,13) 'grid_init_chr: diffs min z: ',abs(minz0-minz0)
         write(*,13) 'grid_init_chr: diffs max x: ',abs(maxx0-maxx0)
         write(*,13) 'grid_init_chr: diffs max y: ',abs(maxy0-maxy0)
         write(*,13) 'grid_init_chr: diffs max z: ',abs(maxz0-maxz0)
      end if

 13   format(A,G20.13)

      ! set all of interior
      do k = 2, nz-1
         do j = 2, ny-1
            do i = 2, nx-1
               chr_array(i,j,k) = CHR_interior
            end do
         end do
      end do

      ! set all AMR boundaries first:
      i = 1
      if (.not.atminx) then
         do k = 1, nz
            do j = 1, ny
               chr_array(i,j,k) = CHR_amr_bdy
            end do
         end do
      end if
      j = 1
      if (.not.atminy) then
         do k = 1, nz
            do i = 1, nx
               chr_array(i,j,k) = CHR_amr_bdy
            end do
         end do
      end if
      k = 1
      if (.not.atminz) then
         do j = 1, ny
            do i = 1, nx
               chr_array(i,j,k) = CHR_amr_bdy
            end do
         end do
      end if

      i = nx
      if (.not.atmaxx) then
         do k = 1, nz
            do j = 1, ny
               chr_array(i,j,k) = CHR_amr_bdy
            end do
         end do
      end if
      j = ny
      if (.not.atmaxy) then
         do k = 1, nz
            do i = 1, nx
               chr_array(i,j,k) = CHR_amr_bdy
            end do
         end do
      end if
      k = nz
      if (.not.atmaxz) then
         do j = 1, ny
            do i = 1, nx
               chr_array(i,j,k) = CHR_amr_bdy
            end do
         end do
      end if

      ! Go back, and re-set all outer boundaries
      i = 1
      if (atminx) then
         do k = 1, nz
            do j = 1, ny
               chr_array(i,j,k) = CHR_xmin_bdy
            end do
         end do
      end if
      j = 1
      if (atminy) then
         do k = 1, nz
            do i = 1, nx
               chr_array(i,j,k) = CHR_ymin_bdy
            end do
         end do
      end if
      k = 1
      if (atminz) then
         do j = 1, ny
            do i = 1, nx
               chr_array(i,j,k) = CHR_zmin_bdy
            end do
         end do
      end if
      i = nx
      if (atmaxx) then
         do k = 1, nz
            do j = 1, ny
               chr_array(i,j,k) = CHR_xmax_bdy
            end do
         end do
      end if
      j = ny
      if (atmaxy) then
         do k = 1, nz
            do i = 1, nx
               chr_array(i,j,k) = CHR_ymax_bdy
            end do
         end do
      end if
      k = nz
      if (atmaxz) then
         do j = 1, ny
            do i = 1, nx
               chr_array(i,j,k) = CHR_zmax_bdy
            end do
         end do
      end if


      if (atminx .and. atminy) then
         i = 1
         j = 1
         do k = 1, nz
            chr_array(i,j,k) = 1.d0*CHR_xmin_ymin_bdy
         end do
      end if
      if (atminx .and. atmaxy) then
         i = 1
         j = ny
         do k = 1, nz
            chr_array(i,j,k) = 1.d0*CHR_xmin_ymax_bdy
         end do
      end if
      if (atmaxx .and. atmaxy) then
         i = nx
         j = ny
         do k = 1, nz
            chr_array(i,j,k) = 1.d0*CHR_xmax_ymax_bdy
         end do
      end if
      if (atmaxx .and. atminy) then
         i = nx
         j = 1
         do k = 1, nz
            chr_array(i,j,k) = 1.d0*CHR_xmax_ymin_bdy
         end do
      end if
      !
      !
      !
      if (atminx .and. atminz) then
         i = 1
         k = 1
         do j = 1, ny
            chr_array(i,j,k) = 1.d0*CHR_xmin_zmin_bdy
         end do
      end if
      if (atminx .and. atmaxz) then
         i = 1
         k = nz
         do j = 1, ny
            chr_array(i,j,k) = 1.d0*CHR_xmin_zmax_bdy
         end do
      end if
      if (atmaxx .and. atmaxz) then
         i = nx
         k = nz
         do j = 1, ny
            chr_array(i,j,k) = 1.d0*CHR_xmax_zmax_bdy
         end do
      end if
      if (atmaxx .and. atminz) then
         i = nx
         k = 1
         do j = 1, ny
            chr_array(i,j,k) = 1.d0*CHR_xmax_zmin_bdy
         end do
      end if
      !
      !
      !
      if (atminy .and. atminz) then
         j = 1
         k = 1
         do i = 1, nx
            chr_array(i,j,k) = 1.d0*CHR_ymin_zmin_bdy
         end do
      end if
      if (atminy .and. atmaxz) then
         j = 1
         k = nz
         do i = 1, nx
            chr_array(i,j,k) = 1.d0*CHR_ymin_zmax_bdy
         end do
      end if
      if (atmaxy .and. atmaxz) then
         j = ny
         k = nz
         do i = 1, nx
            chr_array(i,j,k) = 1.d0*CHR_ymax_zmax_bdy
         end do
      end if
      if (atmaxy .and. atminz) then
         j = ny
         k = 1
         do i = 1, nx
            chr_array(i,j,k) = 1.d0*CHR_ymax_zmin_bdy
         end do
      end if

      !
      ! Corners:
      !
      if (       double_equal(minx,minx0)
     *     .and. double_equal(miny,miny0)
     *     .and. double_equal(minz,minz0) ) then
         chr_array( 1,  1,  1) = 1.d0*CHR_xmin_ymin_zmin_bdy
      end if
      if (       double_equal(minx,minx0)
     *     .and. double_equal(miny,miny0)
     *     .and. double_equal(maxz,maxz0) ) then
         chr_array( 1,  1, nz) = 1.d0*CHR_xmin_ymin_zmax_bdy
      end if
      if (       double_equal(minx,minx0)
     *     .and. double_equal(maxy,maxy0)
     *     .and. double_equal(minz,minz0) ) then
         chr_array( 1, ny,  1) = 1.d0*CHR_xmin_ymax_zmin_bdy
      end if
      if (       double_equal(minx,minx0)
     *     .and. double_equal(maxy,maxy0)
     *     .and. double_equal(maxz,maxz0) ) then
         chr_array( 1, ny, nz) = 1.d0*CHR_xmin_ymax_zmax_bdy
      end if
      if (       double_equal(maxx,maxx0)
     *     .and. double_equal(miny,miny0)
     *     .and. double_equal(minz,minz0) ) then
         chr_array(nx,  1,  1) = 1.d0*CHR_xmax_ymin_zmin_bdy
      end if
      if (       double_equal(maxx,maxx0)
     *     .and. double_equal(miny,miny0)
     *     .and. double_equal(maxz,maxz0) ) then
         chr_array(nx,  1, nz) = 1.d0*CHR_xmax_ymin_zmax_bdy
      end if
      if (       double_equal(maxx,maxx0)
     *     .and. double_equal(maxy,maxy0)
     *     .and. double_equal(minz,minz0) ) then
         chr_array(nx, ny,  1) = 1.d0*CHR_xmax_ymax_zmin_bdy
      end if
      if (       double_equal(maxx,maxx0)
     *     .and. double_equal(maxy,maxy0)
     *     .and. double_equal(maxz,maxz0) ) then
         chr_array(nx, ny, nz) = 1.d0*CHR_xmax_ymax_zmax_bdy
      end if

      if (assume_symmetry .gt. 0) then
         if(ltrace)write(*,*)'grid_init_chr:Symmetry: ',assume_symmetry
         if (assume_symmetry .eq. 3.or.assume_symmetry.eq.6) then
            if ( minz .le. 0 ) then
            !if ( double_equal(minz,0.d0) ) then
               if(ltrace)write(*,*)'grid_init_chr:Reflect on Low Z'
               k = 1
               do j = 1, ny
                  do i = 1, nx
                     chr_array(i,j,k) = 1.d0*CHR_REFL_bdy
                  end do
               end do
            end if
         end if
         if (assume_symmetry .eq. 2.or.assume_symmetry.eq.6) then
            if ( miny .le. 0 ) then
               if(ltrace)write(*,*)'grid_init_chr:Reflect on Low Y'
               j = 1
               do k = 1, nz
                  do i = 1, nx
                     chr_array(i,j,k) = 1.d0*CHR_REFL_bdy
                  end do
               end do
            end if
         end if
         if (assume_symmetry .eq. 1.or.assume_symmetry.eq.6) then
            if ( minx .le. 0 ) then
               if(ltrace)write(*,*)'grid_init_chr:Reflect on Low X'
               i = 1
               do k = 1, nz
                  do j = 1, ny
                     chr_array(i,j,k) = 1.d0*CHR_REFL_bdy
                  end do
               end do
            end if
         end if
         if (.not. (assume_symmetry.eq.1 .or.
     *              assume_symmetry.eq.2 .or.
     *              assume_symmetry.eq.3 .or.
     *              assume_symmetry.eq.6     ) ) then
            write(*,*)'grid_init_chr:Unknown value of assume_symmetry:',
     *                                                assume_symmetry
         end if
      end if

      if (ltrace) write(*,*) 'grid_init_chr: Done'

      return
      end    ! END: grid_init_chr

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_contains_pt:                                                         cc
cc                  Does this grid contain the given point?                   cc
cc                  NB: For now, strictly in the *interior* of grid.          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_contains_pt( gi, x,y,z )
      implicit     none
      integer      gi
      real(kind=8) x,y,z
      real(kind=8) minx,maxx, miny,maxy, minz,maxz, margin
      real(kind=8) grid_return_h
      external     grid_return_h
      
      !margin = 3.d0 * grid_return_h(gi)
      margin = 0.d0 * grid_return_h(gi)

      call grid_find_bounds(gi,minx,maxx,miny,maxy,minz,maxz)

      grid_contains_pt =     ( x .gt. minx+margin )
     *                  .and.( x .lt. maxx+margin )
     *                  .and.( y .gt. miny+margin )
     *                  .and.( y .lt. maxy+margin )
     *                  .and.( z .gt. minz+margin )
     *                  .and.( z .lt. maxz+margin )

!     grid_contains_pt =     ( x .gt. minx ) .and. ( x .lt. maxx)
!    *                  .and.( y .gt. miny ) .and. ( y .lt. maxy)
!    *                  .and.( z .gt. minz ) .and. ( z .lt. maxz)

      return
      end       ! END: grid_contains_pt

